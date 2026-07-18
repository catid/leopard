#!/usr/bin/env python3
"""Create and audit deterministic nested-thread-safe Leopard2 fuzz campaigns.

The campaign uses ``leopard2_lab.py`` for content-addressed executables,
stable seeds, affinity, memory limits, timeouts, resumable per-job results, and
deterministic merge.  Each sanitizer replay is deliberately a one-CPU,
one-thread job; the lab runner overrides inherited native-runtime defaults and
rejects observed oversubscription before the result can be accepted.  Leak
checking is deliberately a separate sanitizer phase because LeakSanitizer's
Linux stop-the-world implementation may briefly clone a helper process during
normal process teardown.
"""

from __future__ import print_function

import argparse
import hashlib
import json
import math
import os
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
from pathlib import Path

_TOOLS_DIRECTORY = str(Path(__file__).resolve().parent)
if _TOOLS_DIRECTORY not in sys.path:
    sys.path.insert(0, _TOOLS_DIRECTORY)

import leopard2_lab as lab  # noqa: E402


CAMPAIGN_SCHEMA = "leopard2-fuzz-campaign/v4"
AUDIT_SCHEMA = "leopard2-fuzz-campaign-audit/v4"
TARGETS = ("api", "pruned")
CAMPAIGN_NAME = "api-and-pruned-asan-ubsan-no-lsan"
SANITIZER_ENVIRONMENT = {
    "ASAN_OPTIONS":
        "abort_on_error=1:detect_leaks=0:halt_on_error=1:suppressions=",
    "LSAN_OPTIONS": "detect_leaks=0:suppressions=",
    "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1:suppressions=",
}
SANITIZER_SCOPE = {
    "address_sanitizer": True,
    "undefined_behavior_sanitizer": True,
    "leak_sanitizer": False,
    "leak_check_phase": "separate",
}
INSTRUMENTATION_SCHEMA = "leopard2-sanitizer-attestation/v1"
INSTRUMENTATION_ARGUMENT = "--leopard2-sanitizer-attestation-v1"
INSTRUMENTATION_SYMBOLS = (
    "__asan_address_is_poisoned",
    "__asan_init",
    "__ubsan_handle_type_mismatch_v1",
)
LEAK_CAMPAIGN_SCHEMA = "leopard2-fuzz-leak-campaign/v3"
LEAK_RESULT_SCHEMA = "leopard2-fuzz-leak-result/v1"
LEAK_MERGE_SCHEMA = "leopard2-fuzz-leak-merge/v1"
LEAK_AUDIT_SCHEMA = "leopard2-fuzz-leak-campaign-audit/v3"
LEAK_CANARY_SCHEMA = "leopard2-lsan-canary-attestation/v2"
LEAK_CANARY_ARGUMENT = "--leopard2-lsan-canary-v1"
LEAK_CANARY_BYTES = 12345
LEAK_CANARY_EXIT_CODE = 86
LEAK_CANARY_STDOUT = (
    "leopard2_lsan_canary armed bytes={}\n".format(LEAK_CANARY_BYTES))
LEAK_CANARY_DIAGNOSTICS = (
    "ERROR: LeakSanitizer: detected memory leaks",
    "Direct leak of {} byte(s)".format(LEAK_CANARY_BYTES),
    "leo2_test_core_leak_sanitizer_canary",
    "SUMMARY: AddressSanitizer: {} byte(s) leaked in 1 allocation(s).".format(
        LEAK_CANARY_BYTES),
)
LSAN_CONTROL_SYMBOLS = (
    "__lsan_is_turned_off",
    "__lsan_default_suppressions",
)
LEAK_SANITIZER_ENVIRONMENT = {
    "ASAN_OPTIONS":
        "abort_on_error=1:detect_leaks=1:halt_on_error=1:suppressions=",
    "LSAN_OPTIONS": "detect_leaks=1:suppressions=",
    "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1:suppressions=",
}
LEAK_CANARY_SANITIZER_ENVIRONMENT = {
    "ASAN_OPTIONS":
        "abort_on_error=0:detect_leaks=1:halt_on_error=1:suppressions=",
    "LSAN_OPTIONS": "detect_leaks=1:exitcode={}:suppressions=".format(
        LEAK_CANARY_EXIT_CODE),
    "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1:suppressions=",
}
LEAK_SANITIZER_SCOPE = {
    "address_sanitizer": True,
    "undefined_behavior_sanitizer": True,
    "leak_sanitizer": True,
    "process_count_evidence": False,
}


class CampaignError(Exception):
    pass


def _expected_lsan_canary_observed_result():
    stderr_evidence = {
        "allocation_bytes": LEAK_CANARY_BYTES,
        "diagnostic_markers": list(LEAK_CANARY_DIAGNOSTICS),
    }
    return {
        "schema": "leopard2-lsan-canary-result/v1",
        "exit_code": LEAK_CANARY_EXIT_CODE,
        "stdout": {
            "sha256": hashlib.sha256(
                LEAK_CANARY_STDOUT.encode("utf-8")).hexdigest(),
            "size_bytes": len(LEAK_CANARY_STDOUT.encode("utf-8")),
        },
        "stderr_evidence": stderr_evidence,
        "stderr_evidence_sha256": lab._digest(stderr_evidence),
    }


def _nm_symbol_type(line, symbol):
    fields = line.split()
    if len(fields) < 2:
        return None
    name = fields[-1].split("@", 1)[0]
    if name != symbol:
        return None
    type_code = fields[-2]
    if len(type_code) != 1:
        raise CampaignError(
            "cannot classify nm symbol line for {}: {!r}".format(
                symbol, line[-256:]))
    return type_code


def _validate_lsan_control_symbols(symbol_output):
    """Reject linked hooks that can disable or selectively suppress LSan."""
    states = {}
    for symbol in LSAN_CONTROL_SYMBOLS:
        matching_types = [
            type_code for type_code in (
                _nm_symbol_type(line, symbol)
                for line in symbol_output.splitlines())
            if type_code is not None
        ]
        # GNU and LLVM nm use U for an undefined strong reference and
        # lowercase w/v for undefined weak function/object references.  Every
        # other type is a definition (including uppercase W/V weak
        # definitions) and can disable the process or hide selected frames.
        defined_types = sorted(set(matching_types) - {"U", "w", "v"})
        if defined_types:
            raise CampaignError(
                "target defines {} (nm type {})".format(
                    symbol, ",".join(defined_types)))
        states[symbol] = "undefined" if matching_types else "absent"
    return states


def _validate_lsan_canary_process(completed, role):
    stdout = completed.stdout.decode("utf-8", errors="replace")
    stderr = completed.stderr.decode("utf-8", errors="replace")
    if len(stdout) > 4096 or len(stderr) > 2 * 1024 * 1024:
        raise CampaignError(
            "{} LSan canary emitted excessive output".format(role))
    if completed.returncode != LEAK_CANARY_EXIT_CODE:
        raise CampaignError(
            "{} LSan canary exited {}, expected {}".format(
                role, completed.returncode, LEAK_CANARY_EXIT_CODE))
    if stdout != LEAK_CANARY_STDOUT:
        raise CampaignError(
            "{} LSan canary marker is missing or mismatched".format(role))
    mismatched = [
        marker for marker in LEAK_CANARY_DIAGNOSTICS
        if stderr.count(marker) != 1
    ]
    if mismatched:
        raise CampaignError(
            "{} LSan canary requires exactly one diagnostic {!r}".format(
                role, mismatched[0]))
    if (stderr.count("Direct leak of ") != 1 or
            "Indirect leak of " in stderr or
            stderr.count("SUMMARY: AddressSanitizer:") != 1):
        raise CampaignError(
            "{} LSan canary did not report exactly one direct leak".format(
                role))
    if ("ERROR: AddressSanitizer:" in stderr or
            "runtime error:" in stderr or
            "UndefinedBehaviorSanitizer" in stderr):
        raise CampaignError(
            "{} LSan canary emitted an unrelated sanitizer error".format(role))
    return _expected_lsan_canary_observed_result()


def _run_lsan_canary_process(executable, environment, role):
    try:
        process = subprocess.Popen(
            [executable, LEAK_CANARY_ARGUMENT], env=environment,
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, start_new_session=True)
    except (OSError, subprocess.SubprocessError) as error:
        raise CampaignError(
            "{} LSan canary could not launch: {}".format(role, error))
    try:
        stdout, stderr = process.communicate(timeout=30.0)
    except subprocess.TimeoutExpired:
        _signal_leak_process_group(process, signal.SIGTERM)
        try:
            process.communicate(timeout=0.5)
        except subprocess.TimeoutExpired:
            _signal_leak_process_group(process, signal.SIGKILL)
            process.communicate()
        raise CampaignError("{} LSan canary timed out".format(role))
    return subprocess.CompletedProcess(
        args=[executable, LEAK_CANARY_ARGUMENT],
        returncode=process.returncode, stdout=stdout, stderr=stderr)


def _validate_lsan_canary_record(
        value, role, executable_identity, instrumentation_record):
    expected_fields = {
        "schema", "role", "executable_sha256", "argument", "environment",
        "expected_allocation_bytes", "expected_stdout", "diagnostic_markers",
        "exit_policy", "observed_result", "control_hook_probe",
    }
    if (not isinstance(value, dict) or set(value) != expected_fields or
            value.get("schema") != LEAK_CANARY_SCHEMA or
            value.get("role") != role or
            value.get("executable_sha256") != executable_identity["sha256"] or
            value.get("argument") != LEAK_CANARY_ARGUMENT or
            value.get("environment") != _lsan_canary_environment() or
            value.get("expected_allocation_bytes") != LEAK_CANARY_BYTES or
            value.get("expected_stdout") != LEAK_CANARY_STDOUT or
            value.get("diagnostic_markers") !=
                list(LEAK_CANARY_DIAGNOSTICS) or
            value.get("exit_policy") != "exact-lsan-exit-{}".format(
                LEAK_CANARY_EXIT_CODE) or
            value.get("observed_result") !=
                _expected_lsan_canary_observed_result()):
        raise CampaignError(
            "{} target has invalid LSan canary evidence".format(role))
    probe = value.get("control_hook_probe")
    expected_tool = instrumentation_record["binary_probe"]["tool"]
    if (not isinstance(probe, dict) or set(probe) != {
            "schema", "tool", "symbols"} or
            probe.get("schema") != "elf-nm-lsan-control-hooks/v2" or
            probe.get("tool") != expected_tool or
            not isinstance(probe.get("symbols"), dict) or
            set(probe["symbols"]) != set(LSAN_CONTROL_SYMBOLS) or
            any(state not in {"absent", "undefined"}
                for state in probe["symbols"].values())):
        raise CampaignError(
            "{} target has invalid LSan disable-hook evidence".format(role))
    _validate_executable_identity(
        probe["tool"], "{} LSan instrumentation tool".format(role))
    return value


def _collect_lsan_canary(
        executable, role, executable_identity, instrumentation_record):
    try:
        current_identity = lab._file_identity(executable)
    except lab.LabError as error:
        raise CampaignError(
            "{} LSan target identity cannot be checked: {}".format(
                role, error))
    if current_identity != executable_identity:
        raise CampaignError(
            "{} target changed before LSan canary".format(role))
    environment = os.environ.copy()
    environment.update(_lsan_canary_environment())
    nm_path = shutil.which("nm")
    if nm_path is None:
        raise CampaignError("nm is required for LSan canary evidence")
    try:
        symbols_completed = subprocess.run(
            [nm_path, "-g", executable], env=environment,
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=30.0, check=False)
    except (OSError, subprocess.SubprocessError) as error:
        raise CampaignError(
            "{} LSan disable-hook probe could not run: {}".format(
                role, error))
    symbols_stderr = symbols_completed.stderr.decode(
        "utf-8", errors="replace").strip()
    if symbols_completed.returncode != 0 or symbols_stderr:
        raise CampaignError(
            "{} LSan disable-hook probe failed (exit {}, stderr={!r})".format(
                role, symbols_completed.returncode, symbols_stderr[-512:]))
    symbols_stdout = symbols_completed.stdout.decode(
        "utf-8", errors="replace")
    control_hook_states = _validate_lsan_control_symbols(symbols_stdout)
    nm_identity = lab._file_identity(nm_path)
    if nm_identity != instrumentation_record["binary_probe"]["tool"]:
        raise CampaignError(
            "{} LSan probe tool differs from the signed instrumentation tool".format(
                role))
    completed = _run_lsan_canary_process(executable, environment, role)
    observed_result = _validate_lsan_canary_process(completed, role)
    if lab._file_identity(executable) != executable_identity:
        raise CampaignError(
            "{} target changed during LSan canary".format(role))
    record = {
        "schema": LEAK_CANARY_SCHEMA,
        "role": role,
        "executable_sha256": executable_identity["sha256"],
        "argument": LEAK_CANARY_ARGUMENT,
        "environment": _lsan_canary_environment(),
        "expected_allocation_bytes": LEAK_CANARY_BYTES,
        "expected_stdout": LEAK_CANARY_STDOUT,
        "diagnostic_markers": list(LEAK_CANARY_DIAGNOSTICS),
        "exit_policy": "exact-lsan-exit-{}".format(LEAK_CANARY_EXIT_CODE),
        "observed_result": observed_result,
        "control_hook_probe": {
            "schema": "elf-nm-lsan-control-hooks/v2",
            "tool": nm_identity,
            "symbols": control_hook_states,
        },
    }
    return _validate_lsan_canary_record(
        record, role, executable_identity, instrumentation_record)


def _validate_instrumentation_record(value, role, executable_identity):
    expected_attestation_fields = {
        "schema", "role", "driver", "address_compile", "address_runtime",
        "undefined_compile", "undefined_runtime", "core_address_compile",
        "core_undefined_compile",
    }
    if (not isinstance(value, dict) or
            set(value) != {
                "executable_sha256", "attestation", "binary_probe"} or
            value.get("executable_sha256") != executable_identity["sha256"]):
        raise CampaignError(
            "{} target has invalid instrumentation identity".format(role))
    attestation = value.get("attestation")
    if (not isinstance(attestation, dict) or
            set(attestation) != expected_attestation_fields or
            attestation.get("schema") != INSTRUMENTATION_SCHEMA or
            attestation.get("role") != role or
            attestation.get("driver") != "deterministic-replay-v1" or
            any(attestation.get(field) is not True for field in (
                "address_compile", "address_runtime",
                "undefined_compile", "undefined_runtime",
                "core_address_compile", "core_undefined_compile"))):
        raise CampaignError(
            "{} target is not proven ASan+UBSan instrumented".format(role))
    binary_probe = value.get("binary_probe")
    if (not isinstance(binary_probe, dict) or set(binary_probe) != {
            "schema", "tool", "required_symbols", "symbol_table_digest"} or
            binary_probe.get("schema") != "elf-nm-global/v1" or
            binary_probe.get("required_symbols") !=
                list(INSTRUMENTATION_SYMBOLS) or
            not isinstance(binary_probe.get("symbol_table_digest"), str) or
            len(binary_probe["symbol_table_digest"]) != 64 or
            any(character not in "0123456789abcdef"
                for character in binary_probe["symbol_table_digest"] or "")):
        raise CampaignError(
            "{} target has invalid independent binary evidence".format(role))
    _validate_executable_identity(
        binary_probe.get("tool"), "{} instrumentation tool".format(role))
    return value


def _collect_instrumentation(executable, role, executable_identity):
    environment = os.environ.copy()
    environment.update(SANITIZER_ENVIRONMENT)
    try:
        completed = subprocess.run(
            [executable, INSTRUMENTATION_ARGUMENT],
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=10.0,
            check=False)
    except (OSError, subprocess.SubprocessError) as error:
        raise CampaignError(
            "{} sanitizer attestation could not run: {}".format(role, error))
    stdout = completed.stdout.decode("utf-8", errors="replace").strip()
    stderr = completed.stderr.decode("utf-8", errors="replace").strip()
    if (completed.returncode != 0 or stderr or len(stdout) > 4096 or
            len(stdout.splitlines()) != 1):
        raise CampaignError(
            "{} sanitizer attestation failed (exit {}, stderr={!r})".format(
                role, completed.returncode, stderr[-512:]))
    try:
        attestation = json.loads(stdout)
    except ValueError as error:
        raise CampaignError(
            "{} sanitizer attestation is not JSON: {}".format(role, error))
    nm_path = shutil.which("nm")
    if nm_path is None:
        raise CampaignError("nm is required for sanitizer binary evidence")
    try:
        symbols_completed = subprocess.run(
            [nm_path, "-g", executable],
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30.0,
            check=False)
    except (OSError, subprocess.SubprocessError) as error:
        raise CampaignError(
            "{} sanitizer symbol probe could not run: {}".format(role, error))
    symbols_stderr = symbols_completed.stderr.decode(
        "utf-8", errors="replace").strip()
    symbols_stdout = symbols_completed.stdout.decode(
        "utf-8", errors="replace")
    if symbols_completed.returncode != 0 or symbols_stderr:
        raise CampaignError(
            "{} sanitizer symbol probe failed (exit {}, stderr={!r})".format(
                role, symbols_completed.returncode, symbols_stderr[-512:]))
    symbols = sorted(set(
        line.split()[-1] for line in symbols_stdout.splitlines()
        if line.split()))
    missing_symbols = sorted(set(INSTRUMENTATION_SYMBOLS) - set(symbols))
    if missing_symbols:
        raise CampaignError(
            "{} target lacks sanitizer symbols: {}".format(
                role, ", ".join(missing_symbols)))
    record = {
        "executable_sha256": executable_identity["sha256"],
        "attestation": attestation,
        "binary_probe": {
            "schema": "elf-nm-global/v1",
            "tool": lab._file_identity(nm_path),
            "required_symbols": list(INSTRUMENTATION_SYMBOLS),
            "symbol_table_digest": lab._digest(symbols),
        },
    }
    return _validate_instrumentation_record(
        record, role, executable_identity)


def _positive_int(value, label, allow_zero=False):
    if (isinstance(value, bool) or not isinstance(value, int) or
            value < 0 or (value == 0 and not allow_zero)):
        raise CampaignError("{} must be {} integer".format(
            label, "a non-negative" if allow_zero else "a positive"))
    return value


def _campaign_spec(arguments, instrumentation_override=None):
    topology = lab.detect_topology()
    worker_count = min(len(topology["allowed_cpus"]), 128)
    seeds_per_target = (
        arguments.seeds_per_target
        if arguments.seeds_per_target is not None else worker_count)
    _positive_int(worker_count, "visible worker count")
    _positive_int(seeds_per_target, "seeds per target")
    _positive_int(arguments.iterations, "iterations")
    _positive_int(arguments.minimum_memory_mb, "minimum memory", True)
    _positive_int(arguments.rss_limit_mb, "RSS limit")
    if arguments.minimum_memory_mb > arguments.rss_limit_mb:
        raise CampaignError("minimum memory may not exceed the RSS limit")

    executables = {
        "api": str(Path(arguments.api_executable).resolve()),
        "pruned": str(Path(arguments.pruned_executable).resolve()),
    }
    for target, executable in executables.items():
        path = Path(executable)
        if not path.is_file() or not os.access(str(path), os.X_OK):
            raise CampaignError(
                "{} fuzz executable is missing or not executable: {}".format(
                    target, executable))
    executable_identities = {
        target: lab._file_identity(executable)
        for target, executable in executables.items()
    }
    if instrumentation_override is None:
        instrumentation = {
            target: _collect_instrumentation(
                executables[target], target, executable_identities[target])
            for target in TARGETS
        }
    else:
        if (not isinstance(instrumentation_override, dict) or
                set(instrumentation_override) != set(TARGETS)):
            raise CampaignError("test instrumentation roles are invalid")
        instrumentation = {
            target: _validate_instrumentation_record(
                instrumentation_override[target], target,
                executable_identities[target])
            for target in TARGETS
        }

    jobs = []
    for target in TARGETS:
        for index in range(seeds_per_target):
            jobs.append({
                "id": "{}.{:03d}".format(target, index),
                "command": [
                    executables[target], "{seed}", str(arguments.iterations)],
                "env": dict(SANITIZER_ENVIRONMENT),
            })
    return {
        "schema": CAMPAIGN_SCHEMA,
        "root": str(Path(arguments.root).resolve()),
        "workers": worker_count,
        "base_seed": arguments.base_seed,
        "metadata": {
            "campaign": CAMPAIGN_NAME,
            "iterations_per_seed": arguments.iterations,
            "seeds_per_target": seeds_per_target,
            "targets": list(TARGETS),
            "target_executables": executable_identities,
            "target_instrumentation": instrumentation,
            "sanitizer_environment": dict(SANITIZER_ENVIRONMENT),
            "sanitizer_scope": dict(SANITIZER_SCOPE),
            "timeout_seconds": arguments.timeout_seconds,
            "memory_policy": {
                "address_space_limit_mb": 0,
                "minimum_memory_mb": arguments.minimum_memory_mb,
                "rss_limit_mb": arguments.rss_limit_mb,
            },
        },
        "defaults": {
            "timeout_seconds": arguments.timeout_seconds,
            # ASan reserves a very large virtual shadow mapping, so RLIMIT_AS
            # is disabled and the lab's sampled process-session RSS cap is the
            # enforceable per-job bound.
            "memory_mb": 0,
            "minimum_memory_mb": arguments.minimum_memory_mb,
            "rss_limit_mb": arguments.rss_limit_mb,
            "cpu_count": 1,
            "thread_runtime": {
                "max_threads": 1,
                "allow_internal_team": False,
            },
        },
        "jobs": jobs,
    }


def _resolved_path(path, label):
    try:
        return path.resolve()
    except (OSError, RuntimeError) as error:
        raise CampaignError("cannot resolve {}: {}".format(label, error))


def _refuse_executable_destination(destination, protected_paths, label):
    resolved_destination = _resolved_path(destination, label)
    if resolved_destination in protected_paths:
        raise CampaignError("{} may not overwrite campaign tooling".format(label))
    try:
        if (resolved_destination.is_file() and
                os.access(str(resolved_destination), os.X_OK)):
            raise CampaignError(
                "{} may not overwrite an executable".format(label))
    except OSError as error:
        raise CampaignError(
            "cannot inspect {} {}: {}".format(label, destination, error))


def _prepare_create_destination(arguments):
    destination = Path(arguments.output).absolute()
    executable_paths = {
        _resolved_path(Path(arguments.api_executable).absolute(), "API target"),
        _resolved_path(
            Path(arguments.pruned_executable).absolute(), "pruned target"),
    }
    _refuse_executable_destination(
        destination, executable_paths, "manifest output")
    try:
        destination.unlink()
    except FileNotFoundError:
        pass
    except OSError as error:
        raise CampaignError(
            "cannot clear prior manifest {}: {}".format(destination, error))
    return destination


def create_manifest(arguments):
    destination = _prepare_create_destination(arguments)
    specification = _campaign_spec(arguments)
    manifest = lab.build_manifest(specification)
    _validate_campaign_manifest(manifest)
    lab._atomic_write_json(destination, manifest)
    print("wrote {} nested-safe fuzz jobs to {} (digest {})".format(
        len(manifest["jobs"]), destination,
        manifest["manifest_digest"]))


def _load_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise CampaignError("cannot read {}: {}".format(path, error))


def _positive_number(value, label):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or
            not math.isfinite(value) or value <= 0):
        raise CampaignError("{} must be a positive number".format(label))
    return value


def _validate_executable_identity(value, label):
    if (not isinstance(value, dict) or
            set(value) != {"path", "sha256", "size_bytes"} or
            not isinstance(value.get("path"), str) or
            not Path(value["path"]).is_absolute() or
            not isinstance(value.get("sha256"), str) or
            len(value["sha256"]) != 64 or
            any(character not in "0123456789abcdef"
                for character in value["sha256"]) or
            isinstance(value.get("size_bytes"), bool) or
            not isinstance(value.get("size_bytes"), int) or
            value["size_bytes"] < 0):
        raise CampaignError("{} has an invalid executable identity".format(label))
    return value


def _validate_campaign_manifest(value):
    """Bind an audit to the exact canonical API/pruned sanitizer campaign."""
    manifest = lab.validate_manifest(value)
    source_spec = manifest.get("source_spec", {})
    if source_spec.get("schema") != CAMPAIGN_SCHEMA:
        raise CampaignError("manifest is not a versioned fuzz campaign")
    metadata = source_spec.get("metadata")
    expected_metadata_fields = {
        "campaign", "iterations_per_seed", "seeds_per_target", "targets",
        "target_executables", "target_instrumentation",
        "sanitizer_environment", "timeout_seconds", "sanitizer_scope",
        "memory_policy",
    }
    if (not isinstance(metadata, dict) or
            set(metadata) != expected_metadata_fields or
            metadata.get("campaign") != CAMPAIGN_NAME or
            metadata.get("targets") != list(TARGETS) or
            metadata.get("sanitizer_environment") != SANITIZER_ENVIRONMENT or
            metadata.get("sanitizer_scope") != SANITIZER_SCOPE):
        raise CampaignError("manifest campaign metadata is not canonical")

    seeds_per_target = metadata.get("seeds_per_target")
    iterations = metadata.get("iterations_per_seed")
    _positive_int(seeds_per_target, "manifest seeds per target")
    _positive_int(iterations, "manifest iterations per seed")
    timeout_seconds = _positive_number(
        metadata.get("timeout_seconds"), "manifest timeout")
    memory_policy = metadata.get("memory_policy")
    if (not isinstance(memory_policy, dict) or set(memory_policy) != {
            "address_space_limit_mb", "minimum_memory_mb", "rss_limit_mb"}):
        raise CampaignError("manifest campaign memory policy is not canonical")
    address_space_limit_mb = memory_policy.get("address_space_limit_mb")
    minimum_memory_mb = memory_policy.get("minimum_memory_mb")
    rss_limit_mb = memory_policy.get("rss_limit_mb")
    _positive_int(address_space_limit_mb, "manifest address-space limit", True)
    _positive_int(minimum_memory_mb, "manifest minimum memory", True)
    _positive_int(rss_limit_mb, "manifest RSS limit")
    if (address_space_limit_mb != 0 or
            minimum_memory_mb > rss_limit_mb):
        raise CampaignError("manifest campaign memory policy is invalid")

    executable_identities = metadata.get("target_executables")
    if (not isinstance(executable_identities, dict) or
            set(executable_identities) != set(TARGETS)):
        raise CampaignError("manifest campaign executable roles are invalid")
    for target in TARGETS:
        _validate_executable_identity(
            executable_identities[target], "{} target".format(target))
    instrumentation = metadata.get("target_instrumentation")
    if (not isinstance(instrumentation, dict) or
            set(instrumentation) != set(TARGETS)):
        raise CampaignError("manifest instrumentation roles are invalid")
    for target in TARGETS:
        _validate_instrumentation_record(
            instrumentation[target], target, executable_identities[target])

    base_seed = manifest.get("base_seed")
    _positive_int(base_seed, "manifest base seed", True)
    topology = manifest.get("topology")
    allowed_cpus = topology.get("allowed_cpus") if isinstance(topology, dict) else None
    if (not isinstance(allowed_cpus, list) or not allowed_cpus or
            not all(isinstance(cpu, int) and not isinstance(cpu, bool) and
                    cpu >= 0 for cpu in allowed_cpus) or
            allowed_cpus != sorted(set(allowed_cpus))):
        raise CampaignError("manifest campaign CPU topology is invalid")
    expected_workers = min(len(allowed_cpus), 128)
    if (manifest.get("workers") != expected_workers or
            manifest.get("cpu_policy") != "physical-first"):
        raise CampaignError("manifest campaign worker policy is not canonical")
    try:
        cpu_order = lab._cpu_order(topology, "physical-first")
    except (KeyError, TypeError, lab.LabError) as error:
        raise CampaignError(
            "manifest campaign CPU topology is invalid: {}".format(error))
    if (not all(isinstance(cpu, int) and not isinstance(cpu, bool) and
                cpu >= 0 for cpu in cpu_order) or
            len(cpu_order) != len(allowed_cpus) or
            set(cpu_order) != set(allowed_cpus)):
        raise CampaignError("manifest campaign CPU order is invalid")

    expected_roles = {}
    for target in TARGETS:
        for index in range(seeds_per_target):
            expected_roles["{}.{:03d}".format(target, index)] = target
    jobs = manifest["jobs"]
    if [job["id"] for job in jobs] != sorted(expected_roles):
        raise CampaignError("manifest campaign job IDs are incomplete or noncanonical")
    expected_job_fields = {
        "id", "command", "cwd", "env", "timeout_seconds", "memory_mb",
        "rss_limit_mb", "minimum_memory_mb", "cpu_set", "thread_runtime",
        "seed", "executable", "job_digest",
    }
    expected_thread_runtime = {
        "schema": lab.THREAD_POLICY_SCHEMA,
        "max_threads": 1,
        "allow_internal_team": False,
        "effective_env": lab._thread_environment(1),
    }
    for position, job in enumerate(jobs):
        job_id = job["id"]
        target = expected_roles[job_id]
        identity = executable_identities[target]
        expected_cpu_set = [cpu_order[position % len(cpu_order)]]
        if (set(job) != expected_job_fields or
                job.get("command") != [
                    identity["path"], "{seed}", str(iterations)] or
                job.get("cwd") != "." or
                job.get("env") != SANITIZER_ENVIRONMENT or
                job.get("timeout_seconds") != timeout_seconds or
                job.get("memory_mb") != address_space_limit_mb or
                job.get("minimum_memory_mb") != minimum_memory_mb or
                job.get("rss_limit_mb") != rss_limit_mb or
                job.get("cpu_set") != expected_cpu_set or
                job.get("thread_runtime") != expected_thread_runtime or
                job.get("seed") != lab._stable_seed(base_seed, job_id) or
                job.get("executable") != identity):
            raise CampaignError(
                "job {} is not a canonical {} sanitizer replay".format(
                    job_id, target))
    return manifest


def _prepare_audit_destination(arguments):
    """Clear stale output before validation without deleting input evidence."""
    destination = Path(arguments.output).absolute()
    resolved_destination = _resolved_path(destination, "audit output")
    resolved_manifest = _resolved_path(
        Path(arguments.manifest).absolute(), "audit manifest")
    resolved_results = _resolved_path(
        Path(arguments.output_dir).absolute(), "audit results")
    if resolved_destination == resolved_manifest:
        raise CampaignError("audit output may not overwrite its manifest")
    if (resolved_destination == resolved_results or
            resolved_results in resolved_destination.parents):
        raise CampaignError("audit output may not overwrite per-job evidence")
    protected_paths = set()
    try:
        untrusted = _load_json(arguments.manifest)
    except CampaignError:
        untrusted = {}
    metadata = (untrusted.get("source_spec", {}).get("metadata", {})
                if isinstance(untrusted, dict) else {})
    executable_roles = metadata.get("target_executables", {})
    if isinstance(executable_roles, dict):
        for identity in executable_roles.values():
            if isinstance(identity, dict) and isinstance(identity.get("path"), str):
                protected_paths.add(_resolved_path(
                    Path(identity["path"]).absolute(),
                    "manifest target executable"))
    instrumentation_roles = metadata.get("target_instrumentation", {})
    if isinstance(instrumentation_roles, dict):
        for record in instrumentation_roles.values():
            probe = record.get("binary_probe", {}) if isinstance(record, dict) else {}
            tool = probe.get("tool", {}) if isinstance(probe, dict) else {}
            if isinstance(tool, dict) and isinstance(tool.get("path"), str):
                protected_paths.add(_resolved_path(
                    Path(tool["path"]).absolute(),
                    "manifest instrumentation tool"))
    _refuse_executable_destination(
        destination, protected_paths, "audit output")
    try:
        destination.unlink()
    except FileNotFoundError:
        pass
    except OSError as error:
        raise CampaignError(
            "cannot clear prior audit output {}: {}".format(
                destination, error))
    return destination


def _verify_live_instrumentation(manifest):
    metadata = manifest["source_spec"]["metadata"]
    for target in TARGETS:
        identity = metadata["target_executables"][target]
        try:
            current_identity = lab._file_identity(identity["path"])
        except lab.LabError as error:
            raise CampaignError(
                "{} target identity cannot be rechecked: {}".format(
                    target, error))
        if current_identity != identity:
            raise CampaignError(
                "{} target changed after campaign creation".format(target))
        current = _collect_instrumentation(
            identity["path"], target, identity)
        if current != metadata["target_instrumentation"][target]:
            raise CampaignError(
                "{} target instrumentation changed after campaign creation".format(
                    target))


def _collect_live_lsan_canaries(source_manifest):
    metadata = source_manifest["source_spec"]["metadata"]
    return {
        target: _collect_lsan_canary(
            metadata["target_executables"][target]["path"], target,
            metadata["target_executables"][target],
            metadata["target_instrumentation"][target])
        for target in TARGETS
    }


def _verify_live_lsan_canaries(leak_manifest):
    current = _collect_live_lsan_canaries(leak_manifest["source_manifest"])
    if current != leak_manifest["target_lsan_canaries"]:
        raise CampaignError(
            "LSan canary evidence changed after companion creation")


def audit_campaign(arguments):
    destination = _prepare_audit_destination(arguments)
    manifest = _validate_campaign_manifest(_load_json(arguments.manifest))
    _verify_live_instrumentation(manifest)
    metadata = manifest.get("source_spec", {}).get("metadata", {})
    expected_per_target = metadata.get("seeds_per_target")
    counts = {target: 0 for target in TARGETS}
    seeds = set()
    for job in manifest["jobs"]:
        target = job["id"].split(".", 1)[0]
        if target not in counts:
            raise CampaignError("unexpected campaign job {}".format(job["id"]))
        counts[target] += 1
        if job["seed"] in seeds:
            raise CampaignError("campaign seeds are not distinct")
        seeds.add(job["seed"])
        policy = job["thread_runtime"]
        if (len(job["cpu_set"]) != 1 or policy["max_threads"] != 1 or
                policy["allow_internal_team"] or
                policy["effective_env"].get("OMP_NUM_THREADS") != "1" or
                policy["effective_env"].get("OMP_DYNAMIC") != "FALSE"):
            raise CampaignError(
                "job {} is not a one-CPU, one-thread workload".format(
                    job["id"]))
        if job.get("rss_limit_mb", 0) <= 0:
            raise CampaignError("job {} has no RSS bound".format(job["id"]))
    if any(count != expected_per_target for count in counts.values()):
        raise CampaignError("campaign target seed counts are incomplete")

    # Keep the lab's ordinary deterministic merge separate.  The requested
    # campaign audit artifact is written only after the stricter live-sample
    # gates below pass, so a failed audit cannot leave a newly forged-looking
    # "audited" output behind.
    merged = lab.merge_results(
        manifest, arguments.output_dir, allow_missing=False)
    if merged["summary"] != {"missing": 0, "success": len(manifest["jobs"])}:
        raise CampaignError(
            "campaign has non-success results: {}".format(merged["summary"]))
    for result in merged["results"]:
        observation = result["thread_runtime"]["observation"]
        # The generic lab runner seeds one launch observation after its
        # affinity-setting pre-exec hook succeeds so sub-millisecond commands
        # can still be resumed.  Fuzz evidence is held to a stronger gate: at
        # least one subsequent /proc session sample must observe affinity and
        # resident memory while the workload is live.
        if (observation["sample_count"] < 2 or
                observation["affinity_sample_count"] < 2 or
                observation["peak_thread_count"] > 1 or
                observation["peak_process_count"] > 1 or
                observation["oversubscribed"] or
                observation["outside_cpu_set"] or
                observation["rss_exceeded"]):
            raise CampaignError(
                "job {} violated its runtime allocation".format(
                    result["job_id"]))
    lab._atomic_write_json(destination, {
        "schema": AUDIT_SCHEMA,
        "manifest_digest": manifest["manifest_digest"],
        "job_count": len(manifest["jobs"]),
        "distinct_seed_count": len(seeds),
        "sanitizer_scope": dict(SANITIZER_SCOPE),
        "target_instrumentation": metadata["target_instrumentation"],
        "summary": merged["summary"],
        "merged_results": merged,
    })
    print("audited {} jobs, {} distinct seeds, one thread per CPU".format(
        len(manifest["jobs"]), len(seeds)))


def _leak_environment():
    environment = dict(LEAK_SANITIZER_ENVIRONMENT)
    environment.update(lab._thread_environment(1))
    return dict(sorted(environment.items()))


def _lsan_canary_environment():
    # An intentional canary must never use the production abort policy: SIGABRT
    # can create cores or invoke host crash reporters.  Compiler-rt's fixed LSan
    # exit code gives the probe one precise, side-effect-free success contract.
    environment = dict(LEAK_CANARY_SANITIZER_ENVIRONMENT)
    environment.update(lab._thread_environment(1))
    return dict(sorted(environment.items()))


def _derive_leak_jobs(source_manifest):
    metadata = source_manifest["source_spec"]["metadata"]
    iterations = metadata["iterations_per_seed"]
    jobs = []
    for source_job in source_manifest["jobs"]:
        target = source_job["id"].split(".", 1)[0]
        identity = metadata["target_executables"][target]
        job = {
            "id": source_job["id"],
            "source_job_digest": source_job["job_digest"],
            "seed": source_job["seed"],
            "iterations": iterations,
            "command": [
                identity["path"], str(source_job["seed"]), str(iterations)],
            "cwd": source_job["cwd"],
            "cpu_set": list(source_job["cpu_set"]),
            "timeout_seconds": source_job["timeout_seconds"],
            "environment": _leak_environment(),
            "executable": lab._json_copy(identity),
        }
        job["job_digest"] = lab._digest(job)
        jobs.append(job)
    return jobs


def _build_leak_manifest(source_manifest, target_lsan_canaries):
    source_manifest = _validate_campaign_manifest(source_manifest)
    metadata = source_manifest["source_spec"]["metadata"]
    if (not isinstance(target_lsan_canaries, dict) or
            set(target_lsan_canaries) != set(TARGETS)):
        raise CampaignError("LSan canary target roles are invalid")
    target_lsan_canaries = {
        target: _validate_lsan_canary_record(
            target_lsan_canaries[target], target,
            metadata["target_executables"][target],
            metadata["target_instrumentation"][target])
        for target in TARGETS
    }
    manifest = {
        "schema": LEAK_CAMPAIGN_SCHEMA,
        "source_manifest_digest": source_manifest["manifest_digest"],
        "source_manifest": lab._json_copy(source_manifest),
        "target_lsan_canaries": lab._json_copy(target_lsan_canaries),
        "execution_policy": {
            "mode": "serial",
            "native_max_threads": 1,
            "effective_environment": _leak_environment(),
            "process_count_evidence": False,
            "lsan_helper_processes": "permitted-not-classified",
        },
        "sanitizer_scope": dict(LEAK_SANITIZER_SCOPE),
        "jobs": _derive_leak_jobs(source_manifest),
    }
    manifest["manifest_digest"] = lab._digest(manifest)
    return _validate_leak_manifest(manifest)


def _validate_leak_manifest(value):
    expected_fields = {
        "schema", "source_manifest_digest", "source_manifest",
        "target_lsan_canaries", "execution_policy", "sanitizer_scope",
        "jobs", "manifest_digest",
    }
    if (not isinstance(value, dict) or set(value) != expected_fields or
            value.get("schema") != LEAK_CAMPAIGN_SCHEMA):
        raise CampaignError("unsupported or malformed leak campaign")
    unsigned = dict(value)
    expected_digest = unsigned.pop("manifest_digest", None)
    if (not isinstance(expected_digest, str) or
            expected_digest != lab._digest(unsigned)):
        raise CampaignError("leak campaign digest does not match its contents")
    source_manifest = _validate_campaign_manifest(value.get("source_manifest"))
    if value.get("source_manifest_digest") != source_manifest["manifest_digest"]:
        raise CampaignError("leak campaign source manifest identity is invalid")
    expected_policy = {
        "mode": "serial",
        "native_max_threads": 1,
        "effective_environment": _leak_environment(),
        "process_count_evidence": False,
        "lsan_helper_processes": "permitted-not-classified",
    }
    if (value.get("execution_policy") != expected_policy or
            value.get("sanitizer_scope") != LEAK_SANITIZER_SCOPE):
        raise CampaignError("leak campaign execution policy is not canonical")
    metadata = source_manifest["source_spec"]["metadata"]
    canaries = value.get("target_lsan_canaries")
    if not isinstance(canaries, dict) or set(canaries) != set(TARGETS):
        raise CampaignError("leak campaign LSan canary roles are invalid")
    for target in TARGETS:
        _validate_lsan_canary_record(
            canaries[target], target,
            metadata["target_executables"][target],
            metadata["target_instrumentation"][target])
    expected_jobs = _derive_leak_jobs(source_manifest)
    if value.get("jobs") != expected_jobs:
        raise CampaignError(
            "leak campaign does not reproduce the exact source seed matrix")
    if len({job["seed"] for job in expected_jobs}) != len(expected_jobs):
        raise CampaignError("leak campaign source seeds are not distinct")
    return value


def _untrusted_campaign_protected_paths(value):
    protected = set()
    if not isinstance(value, dict):
        return protected
    source = value.get("source_manifest", value)
    metadata = (source.get("source_spec", {}).get("metadata", {})
                if isinstance(source, dict) else {})
    executables = metadata.get("target_executables", {})
    if isinstance(executables, dict):
        for identity in executables.values():
            path = identity.get("path") if isinstance(identity, dict) else None
            if isinstance(path, str):
                try:
                    protected.add(Path(path).absolute().resolve())
                except (OSError, RuntimeError, ValueError):
                    pass
    instrumentation = metadata.get("target_instrumentation", {})
    if isinstance(instrumentation, dict):
        for record in instrumentation.values():
            probe = record.get("binary_probe", {}) if isinstance(record, dict) else {}
            tool = probe.get("tool", {}) if isinstance(probe, dict) else {}
            path = tool.get("path") if isinstance(tool, dict) else None
            if isinstance(path, str):
                try:
                    protected.add(Path(path).absolute().resolve())
                except (OSError, RuntimeError, ValueError):
                    pass
    return protected


def _prepare_leak_create_destination(arguments):
    destination = Path(arguments.output).absolute()
    source_path = _resolved_path(
        Path(arguments.manifest).absolute(), "source fuzz manifest")
    if _resolved_path(destination, "leak manifest output") == source_path:
        raise CampaignError("leak manifest may not overwrite its source manifest")
    try:
        untrusted = _load_json(arguments.manifest)
    except CampaignError:
        untrusted = {}
    _refuse_executable_destination(
        destination, _untrusted_campaign_protected_paths(untrusted),
        "leak manifest output")
    try:
        destination.unlink()
    except FileNotFoundError:
        pass
    except OSError as error:
        raise CampaignError(
            "cannot clear prior leak manifest {}: {}".format(
                destination, error))
    return destination


def create_leak_manifest(arguments):
    destination = _prepare_leak_create_destination(arguments)
    source_manifest = _validate_campaign_manifest(_load_json(arguments.manifest))
    _verify_live_instrumentation(source_manifest)
    canaries = _collect_live_lsan_canaries(source_manifest)
    manifest = _build_leak_manifest(source_manifest, canaries)
    lab._atomic_write_json(destination, manifest)
    print("wrote {} exact LSan replay jobs to {} (digest {})".format(
        len(manifest["jobs"]), destination, manifest["manifest_digest"]))


def _leak_job_directory(output_dir, job_id):
    return Path(output_dir) / "jobs" / job_id


def _expected_leak_stdout(job):
    return "leopard2_fuzz_replay seed={} iterations={} passed\n".format(
        job["seed"], job["iterations"]).encode("ascii")


def _validate_leak_result(result_path, result, job):
    required = {
        "schema", "state", "job_id", "job_digest", "outcome", "exit_code",
        "duration_seconds", "seed", "cpu_set", "command", "environment",
        "executable_before", "executable_after", "stdout", "stderr",
        "outputs", "result_digest",
    }
    if (not isinstance(result, dict) or
            not required.issubset(result) or
            not set(result).issubset(required | {"detail"}) or
            result.get("schema") != LEAK_RESULT_SCHEMA or
            result.get("state") != "complete" or
            result.get("job_id") != job["id"] or
            result.get("job_digest") != job["job_digest"] or
            result.get("seed") != job["seed"] or
            result.get("cpu_set") != job["cpu_set"] or
            result.get("command") != job["command"] or
            result.get("environment") != job["environment"] or
            result.get("stdout") != "stdout.txt" or
            result.get("stderr") != "stderr.txt" or
            result.get("outcome") not in {
                "success", "failed", "timeout", "launch_error",
                "diagnostic", "evidence_invalid"} or
            isinstance(result.get("duration_seconds"), bool) or
            not isinstance(result.get("duration_seconds"), (int, float)) or
            not math.isfinite(float(result["duration_seconds"])) or
            result["duration_seconds"] < 0 or
            (result.get("exit_code") is not None and
             (isinstance(result.get("exit_code"), bool) or
              not isinstance(result.get("exit_code"), int))) or
            ("detail" in result and
             (not isinstance(result["detail"], str) or not result["detail"]))):
        raise CampaignError(
            "invalid leak replay result for job {}".format(job["id"]))
    unsigned = dict(result)
    digest = unsigned.pop("result_digest", None)
    if digest != lab._digest(unsigned):
        raise CampaignError(
            "leak replay result digest is invalid for job {}".format(job["id"]))
    outputs = result.get("outputs")
    job_dir = Path(result_path).parent
    for name in ("stdout.txt", "stderr.txt", "result.json"):
        _assert_private_regular_file(
            job_dir / name,
            "LSan job {} {}".format(job["id"], name))
    if (not isinstance(outputs, dict) or set(outputs) != {"stdout", "stderr"}):
        raise CampaignError(
            "leak replay outputs are invalid for job {}".format(job["id"]))
    for name in ("stdout", "stderr"):
        if (not isinstance(outputs.get(name), dict) or
                lab._content_identity(job_dir / (name + ".txt")) !=
                    outputs[name]):
            raise CampaignError(
                "leak replay {} changed for job {}".format(name, job["id"]))
    if result["outcome"] == "success":
        try:
            stdout = (job_dir / "stdout.txt").read_bytes()
            stderr = (job_dir / "stderr.txt").read_bytes()
        except OSError as error:
            raise CampaignError(
                "cannot read leak replay output for job {}: {}".format(
                    job["id"], error))
        if (result.get("exit_code") != 0 or "detail" in result or
                result.get("executable_before") != job["executable"] or
                result.get("executable_after") != job["executable"] or
                stdout != _expected_leak_stdout(job) or stderr):
            raise CampaignError(
                "successful leak replay evidence is invalid for job {}".format(
                    job["id"]))
    return result


def _write_leak_result(job_dir, result):
    for name in ("stdout.txt", "stderr.txt"):
        _assert_private_regular_file(
            Path(job_dir) / name,
            "LSan result {}".format(name))
    _assert_private_regular_file(
        Path(job_dir) / "result.json", "LSan result.json",
        allow_missing=True)
    result["outputs"] = {
        "stdout": lab._content_identity(Path(job_dir) / "stdout.txt"),
        "stderr": lab._content_identity(Path(job_dir) / "stderr.txt"),
    }
    unsigned = dict(result)
    unsigned.pop("result_digest", None)
    result["result_digest"] = lab._digest(unsigned)
    lab._atomic_write_json(Path(job_dir) / "result.json", result)
    return result


def _read_resumable_leak_result(output_dir, job):
    job_dir = _prepare_leak_job_directory(output_dir, job)
    result_path = job_dir / "result.json"
    if not _assert_private_regular_file(
            result_path, "LSan job {} result.json".format(job["id"]),
            allow_missing=True):
        return None
    try:
        result = _load_json(result_path)
        _validate_leak_result(result_path, result, job)
    except (CampaignError, lab.LabError, OSError):
        return None
    return result if result.get("outcome") == "success" else None


def _signal_leak_process_group(process, sig):
    try:
        os.killpg(process.pid, sig)
    except (ProcessLookupError, OSError):
        try:
            process.send_signal(sig)
        except OSError:
            pass


def _assert_private_regular_file(path, label, allow_missing=False):
    path = Path(path)
    try:
        status = path.lstat()
    except FileNotFoundError:
        if allow_missing:
            return False
        raise CampaignError("{} is missing".format(label))
    except OSError as error:
        raise CampaignError("cannot inspect {}: {}".format(label, error))
    if not stat.S_ISREG(status.st_mode):
        raise CampaignError("{} must be a regular non-symlink file".format(label))
    if status.st_nlink != 1:
        raise CampaignError("{} may not be a hard-linked alias".format(label))
    return True


def _prepare_leak_job_directory(output_dir, job):
    jobs_root = Path(output_dir) / "jobs"
    if jobs_root.exists() or jobs_root.is_symlink():
        if jobs_root.is_symlink() or not jobs_root.is_dir():
            raise CampaignError("LSan jobs output must be a real directory")
    else:
        jobs_root.mkdir(parents=True)
    job_dir = jobs_root / job["id"]
    if job_dir.exists() or job_dir.is_symlink():
        if job_dir.is_symlink() or not job_dir.is_dir():
            raise CampaignError(
                "LSan job {} output must be a real directory".format(job["id"]))
    else:
        job_dir.mkdir()
    for name in ("stdout.txt", "stderr.txt", "result.json"):
        _assert_private_regular_file(
            job_dir / name, "LSan job {} {}".format(job["id"], name),
            allow_missing=True)
    return job_dir


def _reset_leak_job_evidence(job_dir, job):
    for name in ("stdout.txt", "stderr.txt", "result.json"):
        path = Path(job_dir) / name
        if _assert_private_regular_file(
                path, "LSan job {} {}".format(job["id"], name),
                allow_missing=True):
            try:
                path.unlink()
            except OSError as error:
                raise CampaignError(
                    "cannot clear LSan job {} {}: {}".format(
                        job["id"], name, error))


def _open_new_leak_output(path, label):
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = None
    try:
        descriptor = os.open(str(path), flags, 0o600)
        status = os.fstat(descriptor)
        if not stat.S_ISREG(status.st_mode) or status.st_nlink != 1:
            raise CampaignError("{} is not a private regular file".format(label))
        return os.fdopen(descriptor, "wb")
    except BaseException:
        try:
            if descriptor is not None:
                os.close(descriptor)
        except OSError:
            pass
        if descriptor is not None:
            try:
                Path(path).unlink()
            except OSError:
                pass
        raise


def _execute_leak_job(manifest, job, output_dir):
    job_dir = _prepare_leak_job_directory(output_dir, job)
    _reset_leak_job_evidence(job_dir, job)
    stdout_path = job_dir / "stdout.txt"
    stderr_path = job_dir / "stderr.txt"
    started = time.monotonic()
    before = None
    after = None
    exit_code = None
    outcome = "launch_error"
    detail = None
    process = None
    stdout_handle = None
    stderr_handle = None
    try:
        stdout_handle = _open_new_leak_output(
            stdout_path, "LSan job {} stdout".format(job["id"]))
        stderr_handle = _open_new_leak_output(
            stderr_path, "LSan job {} stderr".format(job["id"]))
        before = lab._file_identity(job["executable"]["path"])
        if before != job["executable"]:
            raise CampaignError("executable changed before leak replay")
        cwd = Path(job["cwd"])
        if not cwd.is_absolute():
            cwd = Path(manifest["source_manifest"]["root"]) / cwd
        environment = os.environ.copy()
        environment.update(job["environment"])
        process = subprocess.Popen(
            job["command"], cwd=str(cwd), env=environment,
            stdin=subprocess.DEVNULL, stdout=stdout_handle,
            stderr=stderr_handle, start_new_session=True,
            preexec_fn=lambda: lab._child_setup(job["cpu_set"], 0))
        try:
            exit_code = process.wait(timeout=job["timeout_seconds"])
        except subprocess.TimeoutExpired:
            outcome = "timeout"
            detail = "leak replay exceeded its signed timeout"
            _signal_leak_process_group(process, signal.SIGTERM)
            try:
                exit_code = process.wait(timeout=0.5)
            except subprocess.TimeoutExpired:
                _signal_leak_process_group(process, signal.SIGKILL)
                exit_code = process.wait()
        if outcome != "timeout":
            outcome = "success" if exit_code == 0 else "failed"
            if outcome == "failed":
                detail = "leak replay exited {}".format(exit_code)
    except (OSError, subprocess.SubprocessError, lab.LabError,
            CampaignError) as error:
        detail = str(error)
        outcome = "launch_error" if process is None else "evidence_invalid"
    finally:
        if process is not None and process.poll() is None:
            _signal_leak_process_group(process, signal.SIGTERM)
            try:
                process.wait(timeout=0.5)
            except subprocess.TimeoutExpired:
                _signal_leak_process_group(process, signal.SIGKILL)
                process.wait()
        if stdout_handle is not None:
            stdout_handle.close()
        if stderr_handle is not None:
            stderr_handle.close()
    try:
        after = lab._file_identity(job["executable"]["path"])
    except lab.LabError as error:
        after = None
        outcome = "evidence_invalid"
        detail = str(error)
    if before != job["executable"] or after != job["executable"]:
        outcome = "evidence_invalid"
        detail = "executable identity changed during leak replay"
    if outcome == "success":
        try:
            stdout = stdout_path.read_bytes()
            stderr = stderr_path.read_bytes()
        except OSError as error:
            outcome = "evidence_invalid"
            detail = "cannot read leak replay output: {}".format(error)
        else:
            if stderr:
                outcome = "diagnostic"
                detail = "leak replay emitted sanitizer or unexpected stderr"
            elif stdout != _expected_leak_stdout(job):
                outcome = "evidence_invalid"
                detail = "leak replay success marker is missing or mismatched"
    result = {
        "schema": LEAK_RESULT_SCHEMA,
        "state": "complete",
        "job_id": job["id"],
        "job_digest": job["job_digest"],
        "outcome": outcome,
        "exit_code": exit_code,
        "duration_seconds": round(max(0.0, time.monotonic() - started), 6),
        "seed": job["seed"],
        "cpu_set": job["cpu_set"],
        "command": job["command"],
        "environment": job["environment"],
        "executable_before": before,
        "executable_after": after,
        "stdout": "stdout.txt",
        "stderr": "stderr.txt",
    }
    if detail:
        result["detail"] = detail
    return _write_leak_result(job_dir, result)


def _verify_leak_runtime_cpus(manifest):
    if not hasattr(os, "sched_getaffinity"):
        raise CampaignError("LSan companion execution requires Linux affinity")
    allowed = set(os.sched_getaffinity(0))
    unavailable = sorted({
        cpu for job in manifest["jobs"] for cpu in job["cpu_set"]
        if cpu not in allowed
    })
    if unavailable:
        raise CampaignError(
            "LSan companion CPUs are unavailable: {}".format(
                lab.format_cpu_list(unavailable)))


def _prepare_leak_output_directory(manifest_path, manifest, output_dir):
    output_dir = _resolved_path(
        Path(output_dir).absolute(), "LSan result output directory")
    if output_dir.exists() and not output_dir.is_dir():
        raise CampaignError("LSan result output must be a directory")
    output_dir.mkdir(parents=True, exist_ok=True)
    retained_manifest = output_dir / "manifest.json"
    resolved_manifest_path = _resolved_path(
        Path(manifest_path).absolute(), "LSan manifest")
    if _resolved_path(retained_manifest, "retained LSan manifest") == \
            resolved_manifest_path:
        _assert_private_regular_file(
            retained_manifest, "retained LSan manifest")
        retained = _load_json(retained_manifest)
        if retained != manifest:
            raise CampaignError("LSan result manifest differs from its input")
        return output_dir
    if _assert_private_regular_file(
            retained_manifest, "retained LSan manifest", allow_missing=True):
        retained = _load_json(retained_manifest)
        if retained != manifest:
            raise CampaignError("LSan result directory belongs to another manifest")
    else:
        lab._atomic_write_json(retained_manifest, manifest)
        _assert_private_regular_file(
            retained_manifest, "retained LSan manifest")
    return output_dir


def _run_leak_manifest(manifest, manifest_path, output_dir, verify_live=True):
    manifest = _validate_leak_manifest(manifest)
    if verify_live:
        _verify_live_instrumentation(manifest["source_manifest"])
        _verify_live_lsan_canaries(manifest)
    _verify_leak_runtime_cpus(manifest)
    output_dir = _prepare_leak_output_directory(
        manifest_path, manifest, output_dir)
    resumed = 0
    executed = 0
    results = []
    for job in manifest["jobs"]:
        result = _read_resumable_leak_result(output_dir, job)
        if result is not None:
            resumed += 1
        else:
            result = _execute_leak_job(manifest, job, output_dir)
            executed += 1
        results.append(result)
    failures = [
        result["job_id"] for result in results
        if result.get("outcome") != "success"
    ]
    summary = {
        "total": len(results),
        "executed": executed,
        "resumed": resumed,
        "success": len(results) - len(failures),
        "failed": len(failures),
    }
    print(json.dumps(summary, sort_keys=True))
    if failures:
        raise CampaignError(
            "LSan companion has failed jobs: {}".format(", ".join(failures)))
    return summary


def run_leak_manifest(arguments):
    manifest = _validate_leak_manifest(_load_json(arguments.manifest))
    return _run_leak_manifest(
        manifest, arguments.manifest, arguments.output_dir, verify_live=True)


def _prepare_leak_audit_destination(arguments):
    destination = Path(arguments.output).absolute()
    resolved_destination = _resolved_path(destination, "LSan audit output")
    resolved_manifest = _resolved_path(
        Path(arguments.manifest).absolute(), "LSan manifest")
    resolved_results = _resolved_path(
        Path(arguments.output_dir).absolute(), "LSan result output directory")
    if resolved_destination == resolved_manifest:
        raise CampaignError("LSan audit may not overwrite its manifest")
    if (resolved_destination == resolved_results or
            resolved_results in resolved_destination.parents):
        raise CampaignError("LSan audit may not overwrite replay evidence")
    try:
        untrusted = _load_json(arguments.manifest)
    except CampaignError:
        untrusted = {}
    _refuse_executable_destination(
        destination, _untrusted_campaign_protected_paths(untrusted),
        "LSan audit output")
    try:
        destination.unlink()
    except FileNotFoundError:
        pass
    except OSError as error:
        raise CampaignError(
            "cannot clear prior LSan audit {}: {}".format(destination, error))
    return destination


def _audit_leak_campaign(arguments, verify_live=True):
    destination = _prepare_leak_audit_destination(arguments)
    manifest = _validate_leak_manifest(_load_json(arguments.manifest))
    if verify_live:
        _verify_live_instrumentation(manifest["source_manifest"])
        _verify_live_lsan_canaries(manifest)
    output_dir = _resolved_path(
        Path(arguments.output_dir).absolute(), "LSan result output directory")
    if not output_dir.is_dir():
        raise CampaignError("LSan result output directory is missing")
    _assert_private_regular_file(
        output_dir / "manifest.json", "retained LSan manifest")
    retained_manifest = _load_json(output_dir / "manifest.json")
    if retained_manifest != manifest:
        raise CampaignError("LSan result directory manifest is stale or mismatched")
    expected_root_entries = {"manifest.json", "jobs"}
    try:
        root_entries = {entry.name for entry in output_dir.iterdir()}
    except OSError as error:
        raise CampaignError("cannot inspect LSan result directory: {}".format(error))
    if root_entries != expected_root_entries:
        raise CampaignError("LSan result directory has unexpected or missing files")
    jobs_root = output_dir / "jobs"
    if jobs_root.is_symlink() or not jobs_root.is_dir():
        raise CampaignError("LSan jobs evidence must be a real directory")
    expected_ids = [job["id"] for job in manifest["jobs"]]
    try:
        job_entries = list(jobs_root.iterdir())
    except OSError as error:
        raise CampaignError("cannot inspect LSan job evidence: {}".format(error))
    if any(entry.is_symlink() or not entry.is_dir() for entry in job_entries):
        raise CampaignError("LSan job evidence contains a non-directory alias")
    actual_ids = sorted(entry.name for entry in job_entries)
    if actual_ids != sorted(expected_ids):
        raise CampaignError("LSan result seed matrix is incomplete or has extras")
    results = []
    for job in manifest["jobs"]:
        job_dir = _leak_job_directory(output_dir, job["id"])
        if job_dir.is_symlink() or not job_dir.is_dir():
            raise CampaignError(
                "LSan job {} evidence must be a real directory".format(
                    job["id"]))
        try:
            entries = {entry.name for entry in job_dir.iterdir()}
        except OSError as error:
            raise CampaignError(
                "cannot inspect LSan job {}: {}".format(job["id"], error))
        if entries != {"stdout.txt", "stderr.txt", "result.json"}:
            raise CampaignError(
                "LSan job {} has unexpected or missing evidence".format(job["id"]))
        for name in ("stdout.txt", "stderr.txt", "result.json"):
            _assert_private_regular_file(
                job_dir / name,
                "LSan job {} {}".format(job["id"], name))
        result_path = job_dir / "result.json"
        result = _validate_leak_result(
            result_path, _load_json(result_path), job)
        if result["outcome"] != "success":
            raise CampaignError(
                "LSan job {} did not succeed".format(job["id"]))
        results.append(result)
    seeds = {job["seed"] for job in manifest["jobs"]}
    if len(seeds) != len(manifest["jobs"]):
        raise CampaignError("LSan source campaign seeds are not distinct")
    merged = {
        "schema": LEAK_MERGE_SCHEMA,
        "leak_manifest_digest": manifest["manifest_digest"],
        "summary": {"missing": 0, "success": len(results)},
        "results": results,
    }
    audit = {
        "schema": LEAK_AUDIT_SCHEMA,
        "leak_manifest_digest": manifest["manifest_digest"],
        "source_manifest_digest": manifest["source_manifest_digest"],
        "job_count": len(results),
        "distinct_seed_count": len(seeds),
        "total_iterations": sum(job["iterations"] for job in manifest["jobs"]),
        "execution_policy": lab._json_copy(manifest["execution_policy"]),
        "sanitizer_scope": dict(LEAK_SANITIZER_SCOPE),
        "target_lsan_canaries": lab._json_copy(
            manifest["target_lsan_canaries"]),
        "target_instrumentation": lab._json_copy(
            manifest["source_manifest"]["source_spec"]["metadata"]
            ["target_instrumentation"]),
        "summary": merged["summary"],
        "merged_results": merged,
    }
    lab._atomic_write_json(destination, audit)
    print("audited {} exact LSan jobs and {} iterations".format(
        len(results), audit["total_iterations"]))
    return audit


def audit_leak_campaign(arguments):
    return _audit_leak_campaign(arguments, verify_live=True)


def _leak_companion_self_test(directory, nm_identity):
    weak_symbols = "".join(
        "                 w {}\n".format(symbol)
        for symbol in LSAN_CONTROL_SYMBOLS)
    if (_validate_lsan_control_symbols(weak_symbols) != {
            symbol: "undefined" for symbol in LSAN_CONTROL_SYMBOLS} or
            _validate_lsan_control_symbols("") != {
                symbol: "absent" for symbol in LSAN_CONTROL_SYMBOLS}):
        raise CampaignError("LSan control-hook classification is invalid")
    for symbol in LSAN_CONTROL_SYMBOLS:
        for type_code in ("T", "D", "W", "V"):
            try:
                _validate_lsan_control_symbols(
                    "0000000000001000 {} {}\n".format(
                        type_code, symbol))
            except CampaignError:
                pass
            else:
                raise CampaignError(
                    "LSan control hook {} {} definition was accepted".format(
                        symbol, type_code))

    valid_canary_process = subprocess.CompletedProcess(
        args=["modeled", LEAK_CANARY_ARGUMENT],
        returncode=LEAK_CANARY_EXIT_CODE,
        stdout=LEAK_CANARY_STDOUT.encode("utf-8"),
        stderr=("\n".join(LEAK_CANARY_DIAGNOSTICS) + "\n").encode(
            "utf-8"))
    _validate_lsan_canary_process(valid_canary_process, "modeled")
    for label, returncode, stdout, stderr in (
            ("silent success", 0, valid_canary_process.stdout,
             valid_canary_process.stderr),
            ("arbitrary abort", -signal.SIGABRT, valid_canary_process.stdout,
             valid_canary_process.stderr),
            ("wrong linked-core size", LEAK_CANARY_EXIT_CODE,
             b"leopard2_lsan_canary armed bytes=4096\n",
             valid_canary_process.stderr),
            ("duplicate leak report", LEAK_CANARY_EXIT_CODE,
             valid_canary_process.stdout,
             valid_canary_process.stderr * 2),
            ("extra indirect leak", LEAK_CANARY_EXIT_CODE,
             valid_canary_process.stdout,
             valid_canary_process.stderr +
             b"Indirect leak of 1 byte(s) in 1 object(s) allocated from:\n"),
            ("missing compiler-rt diagnostic", LEAK_CANARY_EXIT_CODE,
             valid_canary_process.stdout, b"")):
        try:
            _validate_lsan_canary_process(
                subprocess.CompletedProcess(
                    args=["modeled", LEAK_CANARY_ARGUMENT],
                    returncode=returncode, stdout=stdout, stderr=stderr),
                "modeled")
        except CampaignError:
            pass
        else:
            raise CampaignError("LSan canary accepted {}".format(label))

    helper = Path(directory) / "fake-exact-leak-replay.py"
    helper.write_text(
        "#!{}\n"
        "import os, sys, time\n"
        "if len(sys.argv) != 3:\n"
        "    raise SystemExit(2)\n"
        "seed = int(sys.argv[1], 0)\n"
        "iterations = int(sys.argv[2], 0)\n"
        "asan = dict(part.split('=', 1) for part in "
        "os.environ.get('ASAN_OPTIONS', '').split(':') if '=' in part)\n"
        "lsan = dict(part.split('=', 1) for part in "
        "os.environ.get('LSAN_OPTIONS', '').split(':') if '=' in part)\n"
        "if lsan.get('detect_leaks', asan.get('detect_leaks', '1')) != '1':\n"
        "    raise SystemExit(3)\n"
        "pid = os.fork()\n"
        "if pid == 0:\n"
        "    time.sleep(0.08)\n"
        "    os._exit(0)\n"
        "os.waitpid(pid, 0)\n"
        "print('leopard2_fuzz_replay seed={{}} iterations={{}} passed'.format("
        "seed, iterations))\n".format(sys.executable),
        encoding="utf-8")
    helper.chmod(0o700)
    identity = lab._file_identity(helper)
    instrumentation = {
        target: {
            "executable_sha256": identity["sha256"],
            "attestation": {
                "schema": INSTRUMENTATION_SCHEMA,
                "role": target,
                "driver": "deterministic-replay-v1",
                "address_compile": True,
                "address_runtime": True,
                "undefined_compile": True,
                "undefined_runtime": True,
                "core_address_compile": True,
                "core_undefined_compile": True,
            },
            "binary_probe": {
                "schema": "elf-nm-global/v1",
                "tool": nm_identity,
                "required_symbols": list(INSTRUMENTATION_SYMBOLS),
                "symbol_table_digest": "2" * 64,
            },
        }
        for target in TARGETS
    }
    source_arguments = argparse.Namespace(
        api_executable=str(helper),
        pruned_executable=str(helper),
        output=str(Path(directory) / "unused-source.json"),
        root=directory,
        seeds_per_target=2,
        iterations=19,
        base_seed=29,
        timeout_seconds=5.0,
        minimum_memory_mb=0,
        rss_limit_mb=128,
    )
    source = lab.build_manifest(
        _campaign_spec(source_arguments, instrumentation))
    canaries = {
        target: {
            "schema": LEAK_CANARY_SCHEMA,
            "role": target,
            "executable_sha256": identity["sha256"],
            "argument": LEAK_CANARY_ARGUMENT,
            "environment": _lsan_canary_environment(),
            "expected_allocation_bytes": LEAK_CANARY_BYTES,
            "expected_stdout": LEAK_CANARY_STDOUT,
            "diagnostic_markers": list(LEAK_CANARY_DIAGNOSTICS),
            "exit_policy": "exact-lsan-exit-{}".format(
                LEAK_CANARY_EXIT_CODE),
            "observed_result": _expected_lsan_canary_observed_result(),
            "control_hook_probe": {
                "schema": "elf-nm-lsan-control-hooks/v2",
                "tool": nm_identity,
                "symbols": {
                    symbol: "absent" for symbol in LSAN_CONTROL_SYMBOLS
                },
            },
        }
        for target in TARGETS
    }
    leak = _build_leak_manifest(source, canaries)
    if (leak["schema"] != LEAK_CAMPAIGN_SCHEMA or
            [job["id"] for job in leak["jobs"]] !=
                [job["id"] for job in source["jobs"]] or
            [job["seed"] for job in leak["jobs"]] !=
                [job["seed"] for job in source["jobs"]] or
            any(job["iterations"] != 19 for job in leak["jobs"]) or
            any(job["environment"] != _leak_environment()
                for job in leak["jobs"])):
        raise CampaignError("LSan companion changed the source seed matrix")

    coverage_arguments = argparse.Namespace(**vars(source_arguments))
    coverage_arguments.seeds_per_target = 30
    coverage_arguments.iterations = 8192
    coverage_source = lab.build_manifest(
        _campaign_spec(coverage_arguments, instrumentation))
    coverage_leak = _build_leak_manifest(coverage_source, canaries)
    if (len(coverage_leak["jobs"]) != 60 or
            len({job["seed"] for job in coverage_leak["jobs"]}) != 60 or
            sum(job["iterations"] for job in coverage_leak["jobs"]) !=
                60 * 8192 or
            [job["id"] for job in coverage_leak["jobs"]] !=
                [job["id"] for job in coverage_source["jobs"]]):
        raise CampaignError("LSan companion did not preserve the 60-job gate")

    def clone(value):
        return json.loads(json.dumps(value))

    def resign_manifest(value, job=None):
        if job is not None:
            unsigned_job = dict(job)
            unsigned_job.pop("job_digest", None)
            job["job_digest"] = lab._digest(unsigned_job)
        unsigned = dict(value)
        unsigned.pop("manifest_digest", None)
        value["manifest_digest"] = lab._digest(unsigned)
        return value

    def expect_invalid_manifest(value, label):
        try:
            _validate_leak_manifest(value)
        except (CampaignError, lab.LabError):
            return
        raise CampaignError("{} leak mutation was accepted".format(label))

    missing_job = clone(leak)
    missing_job["jobs"].pop()
    expect_invalid_manifest(resign_manifest(missing_job), "missing job")
    wrong_seed = clone(leak)
    wrong_seed_job = wrong_seed["jobs"][0]
    wrong_seed_job["seed"] += 1
    wrong_seed_job["command"][1] = str(wrong_seed_job["seed"])
    expect_invalid_manifest(
        resign_manifest(wrong_seed, wrong_seed_job), "seed matrix")
    wrong_lsan = clone(leak)
    wrong_lsan_job = wrong_lsan["jobs"][0]
    wrong_lsan_job["environment"]["LSAN_OPTIONS"] = "detect_leaks=0"
    expect_invalid_manifest(
        resign_manifest(wrong_lsan, wrong_lsan_job), "LSan environment")

    missing_canary = clone(leak)
    missing_canary["target_lsan_canaries"].pop("api")
    expect_invalid_manifest(
        resign_manifest(missing_canary), "missing LSan canary")

    wrong_canary_size = clone(leak)
    wrong_canary_size["target_lsan_canaries"]["api"][
        "expected_allocation_bytes"] = 4096
    expect_invalid_manifest(
        resign_manifest(wrong_canary_size), "wrong linked-core canary size")

    disabled_canary = clone(leak)
    disabled_canary["target_lsan_canaries"]["pruned"][
        "control_hook_probe"]["symbols"]["__lsan_is_turned_off"] = \
        "strong-defined:T"
    expect_invalid_manifest(
        resign_manifest(disabled_canary), "defined LSan turn-off hook")

    leak_path = Path(directory) / "leak-manifest.json"
    output_dir = Path(directory) / "leak-results"
    audit_path = Path(directory) / "leak-audit.json"
    lab._atomic_write_json(leak_path, leak)
    summary = _run_leak_manifest(
        leak, str(leak_path), str(output_dir), verify_live=False)
    if summary != {
            "total": 4, "executed": 4, "resumed": 0,
            "success": 4, "failed": 0}:
        raise CampaignError("LSan companion execution summary is invalid")
    audit_arguments = argparse.Namespace(
        manifest=str(leak_path), output_dir=str(output_dir),
        output=str(audit_path))

    def resign_result(path, update):
        value = _load_json(path)
        update(value)
        unsigned_value = dict(value)
        unsigned_value.pop("result_digest", None)
        value["result_digest"] = lab._digest(unsigned_value)
        lab._atomic_write_json(path, value)

    def expect_audit_failure(label):
        lab._atomic_write_json(audit_path, {"stale": label})
        try:
            _audit_leak_campaign(audit_arguments, verify_live=False)
        except (CampaignError, lab.LabError):
            pass
        else:
            raise CampaignError("LSan audit accepted {}".format(label))
        if audit_path.exists():
            raise CampaignError(
                "failed LSan audit retained stale output for {}".format(label))

    def expect_one_repair(label):
        repaired_summary = _run_leak_manifest(
            leak, str(leak_path), str(output_dir), verify_live=False)
        if (repaired_summary["executed"] != 1 or
                repaired_summary["resumed"] != 3):
            raise CampaignError(
                "LSan companion did not repair {}".format(label))
        _audit_leak_campaign(audit_arguments, verify_live=False)

    audit = _audit_leak_campaign(audit_arguments, verify_live=False)
    if (audit.get("schema") != LEAK_AUDIT_SCHEMA or
            audit.get("job_count") != 4 or
            audit.get("distinct_seed_count") != 4 or
            audit.get("total_iterations") != 76 or
            audit.get("sanitizer_scope") != LEAK_SANITIZER_SCOPE or
            audit.get("target_lsan_canaries") != canaries or
            audit.get("execution_policy", {}).get(
                "process_count_evidence") is not False or
            audit.get("execution_policy", {}).get(
                "native_max_threads") != 1):
        raise CampaignError("LSan companion audit summary is invalid")
    resumed = _run_leak_manifest(
        leak, str(leak_path), str(output_dir), verify_live=False)
    if resumed["resumed"] != 4 or resumed["executed"] != 0:
        raise CampaignError("LSan companion did not resume valid jobs")

    first_job = leak["jobs"][0]
    first_dir = _leak_job_directory(output_dir, first_job["id"])
    (first_dir / "stderr.txt").write_text(
        "ERROR: LeakSanitizer: modeled diagnostic\n", encoding="utf-8")
    resign_result(
        first_dir / "result.json",
        lambda value: value["outputs"].update({
            "stderr": lab._content_identity(first_dir / "stderr.txt")}))
    expect_audit_failure("a sanitizer diagnostic")
    expect_one_repair("diagnostic evidence")

    def make_failed(value):
        value["outcome"] = "failed"
        value["exit_code"] = 1
        value["detail"] = "modeled nonzero exit"

    resign_result(first_dir / "result.json", make_failed)
    expect_audit_failure("a nonzero replay result")
    expect_one_repair("nonzero evidence")

    def change_identity(value):
        value["executable_after"]["sha256"] = "0" * 64

    resign_result(first_dir / "result.json", change_identity)
    expect_audit_failure("a changed executable identity")
    expect_one_repair("identity-mutated evidence")

    (first_dir / "result.json").unlink()
    expect_audit_failure("an incomplete result matrix")
    expect_one_repair("incomplete evidence")

    helper_before = lab._file_identity(helper)
    (first_dir / "stdout.txt").unlink()
    (first_dir / "stdout.txt").symlink_to(helper)
    expect_audit_failure("an aliased per-job output")
    try:
        _run_leak_manifest(
            leak, str(leak_path), str(output_dir), verify_live=False)
    except (CampaignError, lab.LabError):
        pass
    else:
        raise CampaignError("LSan runner accepted an aliased per-job output")
    if lab._file_identity(helper) != helper_before:
        raise CampaignError("LSan runner modified an aliased executable")
    (first_dir / "stdout.txt").unlink()
    expect_one_repair("aliased output evidence")

    (first_dir / "stdout.txt").unlink()
    os.link(str(helper), str(first_dir / "stdout.txt"))
    expect_audit_failure("a hard-linked per-job output")
    try:
        _run_leak_manifest(
            leak, str(leak_path), str(output_dir), verify_live=False)
    except (CampaignError, lab.LabError):
        pass
    else:
        raise CampaignError("LSan runner accepted a hard-linked output")
    if lab._file_identity(helper) != helper_before:
        raise CampaignError("LSan runner modified a hard-linked executable")
    (first_dir / "stdout.txt").unlink()
    expect_one_repair("hard-linked output evidence")

    second_dir = _leak_job_directory(output_dir, leak["jobs"][1]["id"])
    second_real_dir = Path(str(second_dir) + ".saved")
    second_dir.rename(second_real_dir)
    second_dir.symlink_to(second_real_dir, target_is_directory=True)
    expect_audit_failure("an aliased job directory")
    try:
        _run_leak_manifest(
            leak, str(leak_path), str(output_dir), verify_live=False)
    except (CampaignError, lab.LabError):
        pass
    else:
        raise CampaignError("LSan runner accepted an aliased job directory")
    second_dir.unlink()
    second_real_dir.rename(second_dir)
    restored = _run_leak_manifest(
        leak, str(leak_path), str(output_dir), verify_live=False)
    if restored["resumed"] != 4 or restored["executed"] != 0:
        raise CampaignError("LSan companion lost restored directory evidence")
    _audit_leak_campaign(audit_arguments, verify_live=False)

    alias_arguments = argparse.Namespace(
        manifest=str(leak_path), output_dir=str(output_dir), output=str(helper))
    try:
        _prepare_leak_audit_destination(alias_arguments)
    except CampaignError:
        pass
    else:
        raise CampaignError("LSan audit accepted an executable output alias")
    if lab._file_identity(helper) != helper_before:
        raise CampaignError("LSan audit modified its aliased executable")


def self_test():
    with tempfile.TemporaryDirectory(
            prefix="leopard2-fuzz-campaign-self-test-") as directory:
        arguments = argparse.Namespace(
            api_executable=sys.executable,
            pruned_executable=sys.executable,
            output=str(Path(directory) / "manifest.json"),
            root=directory,
            seeds_per_target=2,
            iterations=17,
            base_seed=19,
            timeout_seconds=5.0,
            minimum_memory_mb=1,
            rss_limit_mb=2,
        )
        python_identity = lab._file_identity(sys.executable)
        Path(arguments.output).write_text("stale\n", encoding="utf-8")
        if (_prepare_create_destination(arguments) !=
                Path(arguments.output).absolute() or
                Path(arguments.output).exists()):
            raise CampaignError("create retained a stale manifest")
        alias_arguments = argparse.Namespace(**vars(arguments))
        alias_arguments.output = sys.executable
        try:
            _prepare_create_destination(alias_arguments)
        except CampaignError:
            pass
        else:
            raise CampaignError("create accepted output/executable alias")
        if lab._file_identity(sys.executable) != python_identity:
            raise CampaignError("alias rejection modified the executable")
        nm_path = shutil.which("nm")
        if nm_path is None:
            raise CampaignError("self-test requires nm")
        nm_identity = lab._file_identity(nm_path)
        _leak_companion_self_test(directory, nm_identity)
        test_instrumentation = {
            target: {
                "executable_sha256": python_identity["sha256"],
                "attestation": {
                    "schema": INSTRUMENTATION_SCHEMA,
                    "role": target,
                    "driver": "deterministic-replay-v1",
                    "address_compile": True,
                    "address_runtime": True,
                    "undefined_compile": True,
                    "undefined_runtime": True,
                    "core_address_compile": True,
                    "core_undefined_compile": True,
                },
                "binary_probe": {
                    "schema": "elf-nm-global/v1",
                    "tool": nm_identity,
                    "required_symbols": list(INSTRUMENTATION_SYMBOLS),
                    "symbol_table_digest": "1" * 64,
                },
            }
            for target in TARGETS
        }
        first = lab.build_manifest(
            _campaign_spec(arguments, test_instrumentation))
        second = lab.build_manifest(
            _campaign_spec(arguments, test_instrumentation))
        if lab._canonical_json_bytes(first) != lab._canonical_json_bytes(second):
            raise CampaignError("campaign manifest is not deterministic")
        alias_manifest = Path(directory) / "alias-manifest.json"
        lab._atomic_write_json(alias_manifest, first)
        for label, protected in (
                ("target", Path(sys.executable)),
                ("instrumentation tool", Path(nm_path))):
            protected_before = lab._file_identity(protected)
            alias_audit = argparse.Namespace(
                manifest=str(alias_manifest),
                output_dir=str(Path(directory) / "alias-results"),
                output=str(protected),
            )
            try:
                _prepare_audit_destination(alias_audit)
            except CampaignError:
                pass
            else:
                raise CampaignError(
                    "audit accepted {} output alias".format(label))
            if lab._file_identity(protected) != protected_before:
                raise CampaignError(
                    "audit modified aliased {}".format(label))
        symlink_alias = Path(directory) / "executable-output-link"
        symlink_alias.symlink_to(sys.executable)
        symlink_audit = argparse.Namespace(
            manifest=str(alias_manifest),
            output_dir=str(Path(directory) / "symlink-results"),
            output=str(symlink_alias),
        )
        try:
            _prepare_audit_destination(symlink_audit)
        except CampaignError:
            pass
        else:
            raise CampaignError("audit accepted executable symlink alias")
        if not symlink_alias.is_symlink():
            raise CampaignError("audit removed executable symlink alias")
        _validate_campaign_manifest(first)
        if len(first["jobs"]) != 4 or len({
                job["seed"] for job in first["jobs"]}) != 4:
            raise CampaignError("campaign did not assign distinct stable seeds")
        for job in first["jobs"]:
            if (len(job["cpu_set"]) != 1 or
                    job["thread_runtime"]["max_threads"] != 1 or
                    job["thread_runtime"]["effective_env"][
                        "OMP_NUM_THREADS"] != "1" or
                    job["env"] != SANITIZER_ENVIRONMENT or
                    "detect_leaks=0" not in job["env"]["ASAN_OPTIONS"] or
                    job["env"].get("LSAN_OPTIONS") !=
                        "detect_leaks=0:suppressions=" or
                    job["rss_limit_mb"] != 2):
                raise CampaignError("campaign job is not nested-thread safe")

        def clone(value):
            return json.loads(json.dumps(value))

        def resign(value, job=None):
            if job is not None:
                unsigned_job = dict(job)
                unsigned_job.pop("job_digest", None)
                job["job_digest"] = lab._digest(unsigned_job)
            unsigned_manifest = dict(value)
            unsigned_manifest.pop("manifest_digest", None)
            value["manifest_digest"] = lab._digest(unsigned_manifest)
            return value

        def expect_invalid(value, label):
            try:
                _validate_campaign_manifest(value)
            except (CampaignError, lab.LabError):
                return
            raise CampaignError("{} mutation was accepted".format(label))

        wrong_schema = clone(first)
        wrong_schema["source_spec"]["schema"] = CAMPAIGN_SCHEMA + ".forged"
        expect_invalid(resign(wrong_schema), "campaign schema")

        wrong_command = clone(first)
        wrong_command_job = wrong_command["jobs"][0]
        wrong_command_job["command"] = [
            sys.executable, "-c", "import time; time.sleep(1)"]
        expect_invalid(
            resign(wrong_command, wrong_command_job), "fuzzer command")

        wrong_environment = clone(first)
        wrong_environment_job = wrong_environment["jobs"][0]
        wrong_environment_job["env"]["ASAN_OPTIONS"] = (
            "abort_on_error=1:detect_leaks=1:halt_on_error=1")
        expect_invalid(
            resign(wrong_environment, wrong_environment_job),
            "sanitizer environment")

        wrong_lsan_environment = clone(first)
        wrong_lsan_job = wrong_lsan_environment["jobs"][0]
        wrong_lsan_job["env"]["LSAN_OPTIONS"] = "detect_leaks=1"
        expect_invalid(
            resign(wrong_lsan_environment, wrong_lsan_job),
            "LSan environment")

        wrong_scope = clone(first)
        wrong_scope["source_spec"]["metadata"]["sanitizer_scope"][
            "leak_sanitizer"] = True
        expect_invalid(resign(wrong_scope), "sanitizer scope")

        wrong_memory = clone(first)
        wrong_memory_job = wrong_memory["jobs"][0]
        wrong_memory_job["rss_limit_mb"] += 1
        expect_invalid(
            resign(wrong_memory, wrong_memory_job), "memory policy")

        wrong_identity = clone(first)
        wrong_identity["source_spec"]["metadata"]["target_executables"][
            "api"]["sha256"] = "0" * 64
        expect_invalid(resign(wrong_identity), "executable role identity")

        wrong_instrumentation_hash = clone(first)
        wrong_instrumentation_hash["source_spec"]["metadata"][
            "target_instrumentation"]["api"]["executable_sha256"] = "0" * 64
        expect_invalid(
            resign(wrong_instrumentation_hash), "instrumentation identity")

        missing_ubsan = clone(first)
        missing_ubsan["source_spec"]["metadata"]["target_instrumentation"][
            "pruned"]["attestation"]["undefined_runtime"] = False
        expect_invalid(resign(missing_ubsan), "missing UBSan runtime")

        unsanitized_core = clone(first)
        unsanitized_core["source_spec"]["metadata"][
            "target_instrumentation"]["api"]["attestation"][
                "core_address_compile"] = False
        expect_invalid(resign(unsanitized_core), "unsanitized codec core")

        swapped_role = clone(first)
        swapped_role["source_spec"]["metadata"]["target_instrumentation"][
            "api"]["attestation"]["role"] = "pruned"
        expect_invalid(resign(swapped_role), "instrumentation role")

        try:
            _collect_instrumentation(
                sys.executable, "api", python_identity)
        except CampaignError:
            pass
        else:
            raise CampaignError(
                "unsanitized non-replay executable passed attestation")

        forged_metadata = clone(first["source_spec"]["metadata"])
        forged_metadata["seeds_per_target"] = 1
        forged_manifest = lab.build_manifest({
            "schema": CAMPAIGN_SCHEMA,
            "root": directory,
            "base_seed": arguments.base_seed,
            "metadata": forged_metadata,
            "defaults": {
                "timeout_seconds": arguments.timeout_seconds,
                "memory_mb": 0,
                "minimum_memory_mb": arguments.minimum_memory_mb,
                "rss_limit_mb": arguments.rss_limit_mb,
                "cpu_count": 1,
                "thread_runtime": {
                    "max_threads": 1,
                    "allow_internal_team": False,
                },
            },
            "jobs": [{
                "id": "{}.{:03d}".format(target, 0),
                "command": [
                    sys.executable, "-c", "import time; time.sleep(1)"],
                "env": dict(SANITIZER_ENVIRONMENT),
            } for target in TARGETS],
        })
        expect_invalid(forged_manifest, "non-fuzzer campaign")

        # Model the Linux LeakSanitizer teardown shape without depending on a
        # particular compiler-rt build.  The canonical campaign must suppress
        # that helper so its one-process observation remains meaningful.  A
        # generic leak-enabled job must still trip the unchanged lab runtime
        # gate; otherwise a real application fork could be misclassified as a
        # sanitizer implementation detail.
        if (sys.platform.startswith("linux") and hasattr(os, "fork") and
                Path("/proc/self/stat").is_file()):
            helper = Path(directory) / "fake-sanitizer-replay.py"
            helper.write_text(
                "#!{}\n"
                "import os, time\n"
                "asan_options = dict(part.split('=', 1) for part in "
                "os.environ.get('ASAN_OPTIONS', '').split(':') "
                "if '=' in part)\n"
                "lsan_options = dict(part.split('=', 1) for part in "
                "os.environ.get('LSAN_OPTIONS', '').split(':') "
                "if '=' in part)\n"
                "detect_leaks = lsan_options.get('detect_leaks', "
                "asan_options.get('detect_leaks', '1'))\n"
                "time.sleep(0.15)\n"
                "if detect_leaks != '0':\n"
                "    pid = os.fork()\n"
                "    if pid == 0:\n"
                "        time.sleep(0.8)\n"
                "        os._exit(0)\n"
                "    os.waitpid(pid, 0)\n".format(sys.executable),
                encoding="utf-8")
            helper.chmod(0o700)
            runtime_manifest = lab.build_manifest({
                "root": directory,
                "workers": 2,
                "base_seed": 23,
                "defaults": {
                    "timeout_seconds": 5.0,
                    "memory_mb": 0,
                    "rss_limit_mb": 128,
                    "cpu_count": 1,
                    "thread_runtime": {
                        "max_threads": 1,
                        "allow_internal_team": False,
                    },
                },
                "jobs": [{
                    "id": "modeled-lsan-disabled-{}".format(target),
                    "command": [str(helper), "1", "17"],
                    "env": dict(SANITIZER_ENVIRONMENT),
                } for target in TARGETS],
            })
            runtime_output = Path(directory) / "runtime-results"
            inherited_lsan = os.environ.get("LSAN_OPTIONS")
            os.environ["LSAN_OPTIONS"] = "detect_leaks=1:exitcode=91"
            try:
                runtime_summary = lab.run_manifest(
                    runtime_manifest, runtime_output, quiet=True)
            finally:
                if inherited_lsan is None:
                    os.environ.pop("LSAN_OPTIONS", None)
                else:
                    os.environ["LSAN_OPTIONS"] = inherited_lsan
            if runtime_summary["outcomes"] != {"missing": 0, "success": 2}:
                raise CampaignError(
                    "LSan-disabled campaign replay was rejected: {}".format(
                        runtime_summary["outcomes"]))

            leak_manifest = lab.build_manifest({
                "root": directory,
                "workers": 1,
                "defaults": {
                    "memory_mb": 0,
                    "rss_limit_mb": 128,
                    "cpu_count": 1,
                },
                "jobs": [{
                    "id": "modeled-lsan-helper",
                    "command": [str(helper), "1", "17"],
                    "env": {
                        "ASAN_OPTIONS": (
                            "abort_on_error=1:detect_leaks=1:halt_on_error=1"),
                        "LSAN_OPTIONS": "detect_leaks=1",
                    },
                }],
            })
            leak_output = Path(directory) / "leak-helper-results"
            inherited_lsan = os.environ.get("LSAN_OPTIONS")
            os.environ["LSAN_OPTIONS"] = "detect_leaks=0"
            try:
                leak_summary = lab.run_manifest(
                    leak_manifest, leak_output, quiet=True)
            finally:
                if inherited_lsan is None:
                    os.environ.pop("LSAN_OPTIONS", None)
                else:
                    os.environ["LSAN_OPTIONS"] = inherited_lsan
            leak_result = _load_json(
                leak_output / "jobs" / "modeled-lsan-helper" /
                "result.json")
            leak_observation = leak_result["thread_runtime"]["observation"]
            if (leak_summary["outcomes"] != {
                    "evidence_invalid": 1, "missing": 0} or
                    leak_observation["peak_process_count"] < 2 or
                    not leak_observation["oversubscribed"]):
                raise CampaignError(
                    "generic runtime gate accepted a modeled LSan helper")

        invalid_manifest_path = Path(directory) / "invalid-manifest.json"
        stale_audit_path = Path(directory) / "stale-audit.json"
        lab._atomic_write_json(invalid_manifest_path, {"schema": "invalid"})
        lab._atomic_write_json(stale_audit_path, {"stale": True})
        invalid_audit_arguments = argparse.Namespace(
            manifest=str(invalid_manifest_path),
            output_dir=str(Path(directory) / "unused-results"),
            output=str(stale_audit_path),
        )
        try:
            audit_campaign(invalid_audit_arguments)
        except (CampaignError, lab.LabError):
            pass
        else:
            raise CampaignError("invalid campaign audit unexpectedly passed")
        if stale_audit_path.exists():
            raise CampaignError("failed audit retained stale output")
    print("leopard2_fuzz_campaign self-test: PASS")


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create", help="create a signed lab manifest")
    create.add_argument("--api-executable", required=True)
    create.add_argument("--pruned-executable", required=True)
    create.add_argument("--output", required=True)
    create.add_argument("--root", default=os.getcwd())
    create.add_argument("--seeds-per-target", type=int)
    create.add_argument("--iterations", type=int, default=8192)
    create.add_argument("--base-seed", type=int, default=0x79A1501)
    create.add_argument("--timeout-seconds", type=float, default=3600.0)
    create.add_argument("--minimum-memory-mb", type=int, default=512)
    create.add_argument("--rss-limit-mb", type=int, default=2048)

    audit = subparsers.add_parser("audit", help="validate and merge a campaign")
    audit.add_argument("--manifest", required=True)
    audit.add_argument("--output-dir", required=True)
    audit.add_argument("--output", required=True)
    leak_create = subparsers.add_parser(
        "leak-create", help="derive an exact LSan replay from a v4 manifest")
    leak_create.add_argument("--manifest", required=True)
    leak_create.add_argument("--output", required=True)
    leak_run = subparsers.add_parser(
        "leak-run", help="run or resume the serial exact LSan replay")
    leak_run.add_argument("--manifest", required=True)
    leak_run.add_argument("--output-dir", required=True)
    leak_audit = subparsers.add_parser(
        "leak-audit", help="audit the complete exact LSan replay")
    leak_audit.add_argument("--manifest", required=True)
    leak_audit.add_argument("--output-dir", required=True)
    leak_audit.add_argument("--output", required=True)
    subparsers.add_parser("self-test", help="verify deterministic safe specs")
    return parser


def main(argv=None):
    arguments = _parser().parse_args(argv)
    try:
        if arguments.command == "create":
            create_manifest(arguments)
        elif arguments.command == "audit":
            audit_campaign(arguments)
        elif arguments.command == "leak-create":
            create_leak_manifest(arguments)
        elif arguments.command == "leak-run":
            run_leak_manifest(arguments)
        elif arguments.command == "leak-audit":
            audit_leak_campaign(arguments)
        else:
            self_test()
        return 0
    except (CampaignError, lab.LabError) as error:
        print("leopard2_fuzz_campaign: error: {}".format(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
