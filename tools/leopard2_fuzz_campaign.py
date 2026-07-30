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

try:
    import resource
except ImportError:  # pragma: no cover - LSan companion requires Linux.
    resource = None

_TOOLS_DIRECTORY = str(Path(__file__).resolve().parent)
if _TOOLS_DIRECTORY not in sys.path:
    sys.path.insert(0, _TOOLS_DIRECTORY)

import leopard2_lab as lab  # noqa: E402


HISTORICAL_CAMPAIGN_SCHEMA = "leopard2-fuzz-campaign/v4"
CAMPAIGN_SCHEMA = "leopard2-fuzz-campaign/v5"
AUDIT_SCHEMA = "leopard2-fuzz-campaign-audit/v5"
PROBE_EXECUTION_POLICY_SCHEMA = "leopard2-probe-execution-policy/v1"
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
LEAK_CAMPAIGN_SCHEMA = "leopard2-fuzz-leak-campaign/v4"
LEAK_RESULT_SCHEMA = "leopard2-fuzz-leak-result/v2"
LEAK_MERGE_SCHEMA = "leopard2-fuzz-leak-merge/v2"
LEAK_AUDIT_SCHEMA = "leopard2-fuzz-leak-campaign-audit/v4"
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
LEAK_CAPTURE_LIMITS_BYTES = {
    # RLIMIT_FSIZE is a per-file limit, so use the same hard ceiling for both
    # inherited regular-file descriptors.  Successful stdout remains bound to
    # its exact deterministic marker by the result validator.
    "stdout": 2 * 1024 * 1024,
    "stderr": 2 * 1024 * 1024,
}
LEAK_CAPTURE_ENFORCEMENT = {
    "child_file_size_limit_bytes": max(LEAK_CAPTURE_LIMITS_BYTES.values()),
    "policy": "rlimit-fsize-and-retained-fd-v1",
    "sigxfsz_disposition": "default",
}
LEAK_RESULT_JSON_LIMIT_BYTES = 1024 * 1024
PROBE_ATTESTATION_STDOUT_LIMIT = 4096
PROBE_STDERR_LIMIT = 2 * 1024 * 1024
PROBE_NM_STDOUT_LIMIT = 64 * 1024 * 1024


def _probe_capture_profile(stdout_limit, stderr_limit):
    """Describe the exact limits enforced by one retained-file probe."""
    return {
        "stdout_limit_bytes": stdout_limit,
        "stderr_limit_bytes": stderr_limit,
        "rlimit_fsize_soft_bytes": max(stdout_limit, stderr_limit) + 1,
        "rlimit_fsize_hard_bytes": max(stdout_limit, stderr_limit) + 1,
    }


def _probe_execution_policy():
    """Return the closed-world policy used for live instrumentation probes."""
    return {
        "schema": PROBE_EXECUTION_POLICY_SCHEMA,
        "platform_contract": "linux-procfs-elf",
        "retained_input_descriptors": {
            "roles": ["target-executable", "nm-tool"],
            "open": "single-open-read-only-no-follow",
            "content_identity":
                "sha256-and-size-before-and-after-full-lifecycle",
            "child_path": "proc-self-fd",
            "inheritance": "explicit-pass-fds-only",
        },
        "retained_capture_descriptors": {
            "streams": ["stdout", "stderr"],
            "storage": "private-exclusive-regular-files",
            "directory": "retained-no-follow-descriptor",
            "published_name_binding":
                "name-to-retained-device-and-inode",
            "directory_aba_detection":
                "device-inode-mode-links-size-mtime-ctime",
            "capture_change_detection":
                "device-inode-mode-links-size-mtime-ctime",
            "read": "bounded-pread-from-retained-descriptor-after-cleanup",
        },
        "capture_profiles": {
            "sanitizer_attestation": _probe_capture_profile(
                PROBE_ATTESTATION_STDOUT_LIMIT, PROBE_STDERR_LIMIT),
            "sanitizer_symbol_nm": _probe_capture_profile(
                PROBE_NM_STDOUT_LIMIT, PROBE_STDERR_LIMIT),
            "lsan_control_hook_nm": _probe_capture_profile(
                PROBE_NM_STDOUT_LIMIT, PROBE_STDERR_LIMIT),
            "lsan_canary": _probe_capture_profile(
                4096, 2 * 1024 * 1024),
        },
        "capture_enforcement": {
            "rlimit": "RLIMIT_FSIZE",
            "soft_and_hard_limits": True,
            "limit_formula": "max-stream-limit-plus-one",
            "sigxfsz_disposition": "default",
            "oversize_read": "signed-limit-plus-one-byte",
        },
        "descendant_containment": {
            "child_subreaper": True,
            "containment_token":
                "sha256-of-random-32-bytes-and-probe-label",
            "containment_token_environment":
                "LEO2_LAB_CONTAINMENT_TOKEN",
            "direct_child_baseline": True,
            "retained_process_identity": "pid-and-starttime-ticks",
            "leader_signaling": "popen-sigterm-then-sigkill",
            "descendant_signaling":
                "pidfd-or-verified-pid-and-starttime",
            "raw_process_group_signaling": False,
            "bounded_cleanup_before_capture_read": True,
        },
    }


class CampaignError(Exception):
    pass


class _LeakResumeIdentityError(CampaignError):
    pass


def _leak_resume_checkpoint(phase):
    """Internal deterministic fault-injection seam for resume self-tests."""
    del phase


def _probe_terminal_exception(error):
    return isinstance(error, (KeyboardInterrupt, SystemExit))


def _merge_probe_exception(primary, later, context):
    """Preserve the earliest terminal exception across bounded teardown."""
    if primary is None:
        return later
    if _probe_terminal_exception(primary):
        if hasattr(primary, "add_note"):
            primary.add_note(
                "{}: {}: {}".format(
                    context, type(later).__name__, later))
        return primary
    if _probe_terminal_exception(later):
        if hasattr(later, "add_note"):
            later.add_note(
                "earlier probe failure: {}: {}".format(
                    type(primary).__name__, primary))
        return later
    if hasattr(primary, "add_note"):
        primary.add_note(
            "{}: {}: {}".format(context, type(later).__name__, later))
    return primary


def _raise_probe_exception(error, label):
    if _probe_terminal_exception(error) or isinstance(error, CampaignError):
        raise error
    raise CampaignError(
        "{} failed: {}: {}".format(
            label, type(error).__name__, error)) from error


def _retained_probe_descriptor_identity(descriptor, path, label):
    """Hash one retained regular-file descriptor without changing its offset."""
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise CampaignError(
                "{} is not a retained regular file".format(label))
        hasher = hashlib.sha256()
        offset = 0
        while offset < before.st_size:
            block = os.pread(
                descriptor, min(1024 * 1024, before.st_size - offset), offset)
            if not block:
                break
            hasher.update(block)
            offset += len(block)
        after = os.fstat(descriptor)
    except CampaignError:
        raise
    except OSError as error:
        raise CampaignError(
            "cannot hash retained {}: {}".format(label, error))
    if (offset != before.st_size or before.st_dev != after.st_dev or
            before.st_ino != after.st_ino or before.st_size != after.st_size or
            before.st_mtime_ns != after.st_mtime_ns):
        raise CampaignError(
            "{} changed while its retained descriptor was hashed".format(
                label))
    return {
        "path": str(Path(path).resolve()),
        "sha256": hasher.hexdigest(),
        "size_bytes": after.st_size,
    }


def _open_retained_probe_file(path, expected_identity, label):
    _validate_executable_identity(expected_identity, label)
    try:
        current_identity = lab._file_identity(path)
    except lab.LabError as error:
        raise CampaignError(
            "{} identity cannot be checked: {}".format(label, error))
    if current_identity != expected_identity:
        raise CampaignError("{} changed before probe launch".format(label))
    resolved = Path(expected_identity["path"]).resolve()
    flags = os.O_RDONLY
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(str(resolved), flags)
    except OSError as error:
        raise CampaignError(
            "cannot retain {}: {}".format(label, error))
    record = {
        "descriptor": descriptor,
        "path": str(resolved),
        "expected_identity": expected_identity,
        "label": label,
        "proc_path": "/proc/self/fd/{}".format(descriptor),
    }
    try:
        if (_retained_probe_descriptor_identity(
                descriptor, resolved, label) != expected_identity or
                lab._file_identity(resolved) != expected_identity):
            raise CampaignError(
                "{} changed while its descriptor was retained".format(label))
    except BaseException as error:
        primary_error = error
        try:
            os.close(descriptor)
        except BaseException as close_error:
            primary_error = _merge_probe_exception(
                primary_error, close_error,
                "{} retained descriptor close failed".format(label))
        _raise_probe_exception(
            primary_error, "{} retained descriptor setup".format(label))
    return record


def _verify_retained_probe_file(record):
    expected = record["expected_identity"]
    label = record["label"]
    if (_retained_probe_descriptor_identity(
            record["descriptor"], record["path"], label) != expected or
            lab._file_identity(record["path"]) != expected):
        raise CampaignError(
            "{} changed during the complete probe lifecycle".format(label))


def _retained_probe_argument(name):
    return ("retained-probe-file", name)


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


def _probe_capture_checkpoint(phase, directory_descriptor=None):
    """Internal deterministic fault-injection seam for capture self-tests."""
    del phase, directory_descriptor


def _probe_capture_descriptor_state(descriptor, label, regular):
    try:
        value = os.fstat(descriptor)
    except OSError as error:
        raise CampaignError(
            "cannot inspect retained {}: {}".format(label, error))
    expected_type = stat.S_ISREG if regular else stat.S_ISDIR
    if not expected_type(value.st_mode):
        raise CampaignError(
            "retained {} changed file type".format(label))
    if regular and value.st_nlink != 1:
        raise CampaignError(
            "retained {} must not be hard-linked".format(label))
    return {
        "descriptor": descriptor,
        "device": value.st_dev,
        "inode": value.st_ino,
        "mode": value.st_mode,
        "links": value.st_nlink,
        "size": value.st_size,
        "mtime_ns": value.st_mtime_ns,
        "ctime_ns": value.st_ctime_ns,
    }


def _verify_probe_capture_evidence(
        directory_descriptor, directory_state, captures,
        capture_states=None):
    current_directory = _probe_capture_descriptor_state(
        directory_descriptor, "probe capture directory", regular=False)
    if current_directory != directory_state:
        raise CampaignError(
            "probe capture directory changed during execution")
    for name, descriptor in captures.items():
        label = "probe capture {}".format(name)
        try:
            lab._verify_named_fd(
                directory_descriptor, name, descriptor, label, regular=True)
        except lab.LabError as error:
            raise CampaignError(str(error)) from error
        if capture_states is not None:
            current = _probe_capture_descriptor_state(
                descriptor, label, regular=True)
            if current != capture_states[name]:
                raise CampaignError(
                    "{} changed during evidence read".format(label))


def _read_bounded_probe_capture(descriptor, maximum_bytes, label):
    try:
        return _bounded_leak_capture_from_descriptor(
            descriptor, maximum_bytes, label)
    except CampaignError as error:
        if "exceeds its signed" in str(error):
            raise CampaignError(
                "{} emitted excessive output".format(label)) from error
        raise


def _lsan_cleanup_checkpoint(phase):
    """Internal deterministic fault-injection seam for cleanup self-tests."""
    del phase


def _new_lsan_containment(process, containment_token):
    return {
        "process": process,
        "containment_token": containment_token,
        "contained_identities": {},
        "direct_child_baseline": None,
    }


def _current_lsan_direct_children():
    """Return exact children currently owned by this isolated subreaper."""
    if not sys.platform.startswith("linux"):
        return []
    root = Path("/proc")
    try:
        entries = list(root.iterdir())
    except OSError as error:
        raise CampaignError(
            "cannot enumerate direct children for LSan containment: {}".format(
                error))
    records = []
    for entry in entries:
        if not entry.name.isdigit():
            continue
        parsed = lab._parse_linux_proc_stat(
            lab._read_text(entry / "stat") or "")
        if parsed is not None and parsed["ppid"] == os.getpid():
            records.append(parsed)
    return records


def _begin_lsan_direct_child_capture(active):
    """Snapshot peers before launch so later direct children are attributable."""
    if active.get("direct_child_baseline") is not None:
        raise CampaignError("LSan direct-child capture began twice")
    active["direct_child_baseline"] = {
        (record["pid"], record["starttime_ticks"]): record
        for record in _current_lsan_direct_children()
    }


def _remember_lsan_direct_children(active):
    """Retain exact new direct/adopted children, including unpublished leaders."""
    baseline = active.get("direct_child_baseline")
    if baseline is None:
        return
    known = active["contained_identities"]
    for record in _current_lsan_direct_children():
        key = (record["pid"], record["starttime_ticks"])
        if key not in baseline:
            known[key] = record


def _remember_lsan_descendants(active):
    """Retain exact descendants while the unreaped leader reserves its PID.

    Session IDs are safe attribution evidence only while the original Popen
    child is still unreaped.  Once it is reaped, that PID can be reused as an
    unrelated session ID.  Later cleanup therefore uses only these retained
    PID/start-time identities and the unguessable inherited token.
    """
    _remember_lsan_direct_children(active)
    process = active.get("process")
    if process is None or process.returncode is not None:
        return
    known = active["contained_identities"]
    for record in lab._contained_processes(active):
        if record["pid"] == process.pid:
            continue
        known[(record["pid"], record["starttime_ticks"])] = record


def _wait_lsan_leader(active, timeout_seconds):
    """Poll one leader boundedly, sampling descendants before every reap."""
    process = active["process"]
    deadline = time.monotonic() + max(0.0, float(timeout_seconds))
    while True:
        _remember_lsan_descendants(active)
        returncode = process.poll()
        if returncode is not None:
            return returncode
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return None
        time.sleep(min(0.02, remaining))


def _current_exact_lsan_descendants(active):
    """Return retained/token descendants without a stale session-ID scan."""
    _remember_lsan_direct_children(active)
    known = active["contained_identities"]
    process = active.get("process")
    for record in lab._token_process_records(
            active["containment_token"]):
        # Popen owns the exact leader until it is reaped.  It is terminated
        # separately with Popen.send_signal, never through a token PID scan.
        if (process is not None and process.returncode is None and
                record["pid"] == process.pid):
            continue
        known[(record["pid"], record["starttime_ticks"])] = record
    current = []
    for key, record in list(known.items()):
        observed = lab._same_process_identity(record)
        if observed is None:
            del known[key]
            continue
        known[key] = observed
        current.append(observed)
    return current


def _cleanup_exact_lsan_descendants(active):
    """Boundedly terminate exact retained/token descendants."""
    initial = []
    for attempt in range(2):
        initial = _current_exact_lsan_descendants(active)
        if initial:
            break
        if attempt == 0:
            time.sleep(0.02)
    if not initial:
        return {"detected": [], "remaining": []}

    detected = {
        (record["pid"], record["starttime_ticks"]): record
        for record in initial
    }
    current = initial
    deadline = time.monotonic() + 0.5
    while current and time.monotonic() < deadline:
        for record in current:
            if record["state"] == "Z":
                lab._reap_process_identity(record)
            else:
                lab._signal_process_identity(record, signal.SIGTERM)
        time.sleep(0.02)
        current = _current_exact_lsan_descendants(active)
        for record in current:
            detected[(record["pid"], record["starttime_ticks"])] = record

    deadline = time.monotonic() + 0.75
    while current and time.monotonic() < deadline:
        for record in current:
            if record["state"] == "Z":
                lab._reap_process_identity(record)
            else:
                lab._signal_process_identity(record, signal.SIGKILL)
        time.sleep(0.02)
        current = _current_exact_lsan_descendants(active)
        for record in current:
            detected[(record["pid"], record["starttime_ticks"])] = record

    return {
        "detected": sorted({record["pid"] for record in detected.values()}),
        "remaining": sorted({record["pid"] for record in current}),
    }


def _cleanup_lsan_process_tree(active):
    """Boundedly stop one exact leader and all attributable descendants."""
    _lsan_cleanup_checkpoint("begin")
    process = active.get("process")
    errors = []
    leader_remaining = False
    if process is not None and process.returncode is None:
        _remember_lsan_descendants(active)
        try:
            process.send_signal(signal.SIGTERM)
        except (OSError, ProcessLookupError):
            pass
        if _wait_lsan_leader(active, 0.5) is None:
            try:
                process.send_signal(signal.SIGKILL)
            except (OSError, ProcessLookupError):
                pass
            if _wait_lsan_leader(active, 1.0) is None:
                leader_remaining = True
                errors.append(
                    "leader {} remained live after bounded SIGKILL".format(
                        process.pid))

    _lsan_cleanup_checkpoint("before-descendants")
    containment = _cleanup_exact_lsan_descendants(active)
    if containment["remaining"]:
        errors.append(
            "contained processes remained after cleanup: {}".format(
                ",".join(str(pid) for pid in containment["remaining"])))
    return {
        "detected": containment["detected"],
        "remaining": containment["remaining"],
        "leader_remaining": leader_remaining,
        "errors": errors,
    }


def _cleanup_lsan_process_tree_with_retry(active):
    """Retry teardown once after any cleanup exception.

    KeyboardInterrupt and SystemExit remain authoritative, but teardown gets a
    second bounded attempt before either is re-raised.
    """
    try:
        return _cleanup_lsan_process_tree(active)
    except BaseException as first_error:
        try:
            cleanup = _cleanup_lsan_process_tree(active)
        except BaseException as retry_error:
            if (not isinstance(first_error, (KeyboardInterrupt, SystemExit)) and
                    isinstance(retry_error, (KeyboardInterrupt, SystemExit))):
                if hasattr(retry_error, "add_note"):
                    retry_error.add_note(
                        "initial LSan cleanup failed: {}: {}".format(
                            type(first_error).__name__, first_error))
                raise retry_error
            if hasattr(first_error, "add_note"):
                first_error.add_note(
                    "LSan cleanup retry failed: {}: {}".format(
                        type(retry_error).__name__, retry_error))
            raise first_error
        if isinstance(first_error, (KeyboardInterrupt, SystemExit)):
            cleanup_detail = "; ".join(cleanup["errors"])
            if cleanup_detail and hasattr(first_error, "add_note"):
                first_error.add_note(cleanup_detail)
            raise first_error
        cleanup["errors"].append(
            "initial cleanup attempt failed: {}: {}".format(
                type(first_error).__name__, first_error))
        return cleanup


def _probe_child_setup(capture_file_limit):
    """Apply one inherited regular-file ceiling before executing a probe."""
    if resource is None or not hasattr(resource, "RLIMIT_FSIZE"):
        raise RuntimeError("RLIMIT_FSIZE is required for bounded probes")
    signal.signal(signal.SIGXFSZ, signal.SIG_DFL)
    unused_soft, hard = resource.getrlimit(resource.RLIMIT_FSIZE)
    del unused_soft
    if hard != resource.RLIM_INFINITY and hard < capture_file_limit:
        raise RuntimeError(
            "inherited RLIMIT_FSIZE is below the probe capture limit")
    # The probed executable is not trusted to preserve a soft limit.  Lower
    # the inherited hard limit as well so it cannot raise RLIMIT_FSIZE again
    # after exec and grow either retained capture without bound.
    resource.setrlimit(
        resource.RLIMIT_FSIZE,
        (capture_file_limit, capture_file_limit))


def _run_bounded_probe_process(
        command, environment, label, timeout_seconds,
        stdout_limit, stderr_limit, pass_fds=()):
    """Run one identity-bound probe with bounded files and exact containment."""
    if (not isinstance(stdout_limit, int) or stdout_limit <= 0 or
            not isinstance(stderr_limit, int) or stderr_limit <= 0):
        raise CampaignError("{} has invalid capture limits".format(label))
    capture_file_limit = max(stdout_limit, stderr_limit) + 1
    containment_token = hashlib.sha256(
        os.urandom(32) + label.encode("utf-8")).hexdigest()
    probe_environment = dict(environment)
    probe_environment["LEO2_LAB_CONTAINMENT_TOKEN"] = containment_token
    process = None
    timed_out = False
    primary_error = None
    cleanup = None
    active = _new_lsan_containment(None, containment_token)
    stdout = None
    stderr = None

    try:
        # Adopting orphaned setsid descendants makes their final zombie state
        # waitable by this runner.  Baseline capture also catches a tokenless
        # child forked immediately before Popen itself raises.
        with lab._ChildSubreaperLease():
            _begin_lsan_direct_child_capture(active)
            with tempfile.TemporaryDirectory(
                    prefix="leopard2-bounded-probe-") as temporary:
                stdout_path = Path(temporary) / "stdout.txt"
                stderr_path = Path(temporary) / "stderr.txt"
                directory_descriptor = None
                try:
                    directory_descriptor = os.open(
                        temporary,
                        os.O_RDONLY | os.O_CLOEXEC |
                        os.O_DIRECTORY | os.O_NOFOLLOW)
                    with stdout_path.open("x+b", buffering=0) as stdout_capture, \
                            stderr_path.open("x+b", buffering=0) as stderr_capture:
                        captures = {
                            "stdout.txt": stdout_capture.fileno(),
                            "stderr.txt": stderr_capture.fileno(),
                        }
                        directory_state = _probe_capture_descriptor_state(
                            directory_descriptor,
                            "probe capture directory", regular=False)
                        try:
                            process = subprocess.Popen(
                                command, env=probe_environment,
                                stdin=subprocess.DEVNULL,
                                stdout=stdout_capture, stderr=stderr_capture,
                                start_new_session=True,
                                pass_fds=tuple(pass_fds),
                                preexec_fn=lambda: _probe_child_setup(
                                    capture_file_limit))
                            active["process"] = process
                            if _wait_lsan_leader(
                                    active, timeout_seconds) is None:
                                timed_out = True
                        except BaseException as error:
                            primary_error = error
                        try:
                            cleanup = (
                                _cleanup_lsan_process_tree_with_retry(active))
                        except BaseException as error:
                            primary_error = _merge_probe_exception(
                                primary_error, error,
                                "{} teardown failed".format(label))
                        cleanup_detail = ""
                        if cleanup is not None:
                            cleanup_detail = "; ".join(cleanup["errors"])
                            if cleanup["detected"]:
                                residual_detail = (
                                    "probe left residual descendants {}; "
                                    "cleanup remaining {}".format(
                                        ",".join(
                                            str(pid)
                                            for pid in cleanup["detected"]),
                                        ",".join(
                                            str(pid)
                                            for pid in cleanup["remaining"]) or
                                        "none"))
                                cleanup_detail = (
                                    residual_detail if not cleanup_detail
                                    else cleanup_detail + "; " +
                                    residual_detail)
                        if cleanup_detail:
                            primary_error = _merge_probe_exception(
                                primary_error,
                                CampaignError(cleanup_detail),
                                "{} containment failed".format(label))
                        if timed_out:
                            primary_error = _merge_probe_exception(
                                primary_error,
                                CampaignError("{} timed out".format(label)),
                                "{} timeout".format(label))

                        try:
                            _probe_capture_checkpoint(
                                "after-cleanup", directory_descriptor)
                            _verify_probe_capture_evidence(
                                directory_descriptor, directory_state,
                                captures)
                            capture_states = {
                                name: _probe_capture_descriptor_state(
                                    descriptor,
                                    "probe capture {}".format(name),
                                    regular=True)
                                for name, descriptor in captures.items()
                            }
                            _probe_capture_checkpoint(
                                "before-read", directory_descriptor)
                            _verify_probe_capture_evidence(
                                directory_descriptor, directory_state,
                                captures, capture_states)
                            if primary_error is None:
                                stdout = _read_bounded_probe_capture(
                                    captures["stdout.txt"], stdout_limit,
                                    "{} stdout".format(label))
                                stderr = _read_bounded_probe_capture(
                                    captures["stderr.txt"], stderr_limit,
                                    "{} stderr".format(label))
                            _probe_capture_checkpoint(
                                "after-read", directory_descriptor)
                            _verify_probe_capture_evidence(
                                directory_descriptor, directory_state,
                                captures, capture_states)
                        except BaseException as error:
                            primary_error = _merge_probe_exception(
                                primary_error, error,
                                "{} retained capture verification failed"
                                .format(label))
                except BaseException as error:
                    primary_error = _merge_probe_exception(
                        primary_error, error,
                        "{} capture lifecycle failed".format(label))
                    if cleanup is None:
                        try:
                            cleanup = (
                                _cleanup_lsan_process_tree_with_retry(active))
                        except BaseException as cleanup_error:
                            primary_error = _merge_probe_exception(
                                primary_error, cleanup_error,
                                "{} fallback teardown failed".format(label))
                finally:
                    if directory_descriptor is not None:
                        os.close(directory_descriptor)
    except BaseException as error:
        primary_error = _merge_probe_exception(
            primary_error, error,
            "{} outer lifecycle failed".format(label))

    if primary_error is not None:
        _raise_probe_exception(primary_error, label)

    return subprocess.CompletedProcess(
        args=command, returncode=process.returncode,
        stdout=stdout, stderr=stderr)


def _run_retained_probe(
        command_parts, file_specs, environment, label, timeout_seconds,
        stdout_limit, stderr_limit):
    """Run a probe against retained files and revalidate paths after teardown."""
    retained = {}
    primary_error = None
    result = None
    try:
        for name, spec in file_specs.items():
            path, expected_identity, file_label = spec
            retained[name] = _open_retained_probe_file(
                path, expected_identity, file_label)
        command = []
        for part in command_parts:
            if (isinstance(part, tuple) and len(part) == 2 and
                    part[0] == "retained-probe-file"):
                if part[1] not in retained:
                    raise CampaignError(
                        "{} references an unknown retained file {}".format(
                            label, part[1]))
                command.append(retained[part[1]]["proc_path"])
            elif isinstance(part, str):
                command.append(part)
            else:
                raise CampaignError(
                    "{} has an invalid command argument".format(label))
        result = _run_bounded_probe_process(
            command, environment, label, timeout_seconds,
            stdout_limit, stderr_limit,
            pass_fds=tuple(
                record["descriptor"] for record in retained.values()))
    except BaseException as error:
        primary_error = error
    finally:
        for record in retained.values():
            try:
                _verify_retained_probe_file(record)
            except BaseException as error:
                primary_error = _merge_probe_exception(
                    primary_error, error,
                    "{} retained identity verification failed".format(label))
        for record in retained.values():
            try:
                os.close(record["descriptor"])
            except BaseException as error:
                primary_error = _merge_probe_exception(
                    primary_error, error,
                    "{} retained descriptor close failed".format(label))
    if primary_error is not None:
        _raise_probe_exception(primary_error, label)
    return result


def _run_lsan_canary_process(
        executable, environment, role, timeout_seconds=30.0, pass_fds=()):
    """Run one canary without pipes or an unbounded post-timeout wait."""
    return _run_bounded_probe_process(
        [executable, LEAK_CANARY_ARGUMENT], environment,
        "{} LSan canary".format(role), timeout_seconds,
        4096, 2 * 1024 * 1024, pass_fds=pass_fds)


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
    environment = os.environ.copy()
    environment.update(_lsan_canary_environment())
    nm_path = shutil.which("nm")
    if nm_path is None:
        raise CampaignError("nm is required for LSan canary evidence")
    nm_identity = lab._file_identity(nm_path)
    if nm_identity != instrumentation_record["binary_probe"]["tool"]:
        raise CampaignError(
            "{} LSan probe tool differs from the signed instrumentation tool".format(
                role))

    def collect(retained):
        symbols_completed = _run_bounded_probe_process(
            [retained["tool"]["proc_path"], "-g",
             retained["target"]["proc_path"]],
            environment, "{} LSan disable-hook probe".format(role), 30.0,
            PROBE_NM_STDOUT_LIMIT, PROBE_STDERR_LIMIT,
            pass_fds=(
                retained["tool"]["descriptor"],
                retained["target"]["descriptor"]))
        symbols_stderr = symbols_completed.stderr.decode(
            "utf-8", errors="replace").strip()
        if symbols_completed.returncode != 0 or symbols_stderr:
            raise CampaignError(
                "{} LSan disable-hook probe failed "
                "(exit {}, stderr={!r})".format(
                    role, symbols_completed.returncode,
                    symbols_stderr[-512:]))
        control_hook_states = _validate_lsan_control_symbols(
            symbols_completed.stdout.decode(
                "utf-8", errors="replace"))
        completed = _run_lsan_canary_process(
            retained["target"]["proc_path"], environment, role,
            pass_fds=(retained["target"]["descriptor"],))
        observed_result = _validate_lsan_canary_process(completed, role)
        return {
            "schema": LEAK_CANARY_SCHEMA,
            "role": role,
            "executable_sha256": executable_identity["sha256"],
            "argument": LEAK_CANARY_ARGUMENT,
            "environment": _lsan_canary_environment(),
            "expected_allocation_bytes": LEAK_CANARY_BYTES,
            "expected_stdout": LEAK_CANARY_STDOUT,
            "diagnostic_markers": list(LEAK_CANARY_DIAGNOSTICS),
            "exit_policy": "exact-lsan-exit-{}".format(
                LEAK_CANARY_EXIT_CODE),
            "observed_result": observed_result,
            "control_hook_probe": {
                "schema": "elf-nm-lsan-control-hooks/v2",
                "tool": nm_identity,
                "symbols": control_hook_states,
            },
        }

    retained = {}
    primary_error = None
    record = None
    try:
        retained["target"] = _open_retained_probe_file(
            executable, executable_identity, "{} LSan target".format(role))
        retained["tool"] = _open_retained_probe_file(
            nm_path, nm_identity, "{} LSan instrumentation tool".format(role))
        record = collect(retained)
    except BaseException as error:
        primary_error = error
    finally:
        for retained_record in retained.values():
            try:
                _verify_retained_probe_file(retained_record)
            except BaseException as error:
                primary_error = _merge_probe_exception(
                    primary_error, error,
                    "{} LSan retained identity verification failed".format(
                        role))
        for retained_record in retained.values():
            try:
                os.close(retained_record["descriptor"])
            except BaseException as error:
                primary_error = _merge_probe_exception(
                    primary_error, error,
                    "{} LSan retained descriptor close failed".format(role))
    if primary_error is not None:
        _raise_probe_exception(
            primary_error, "{} LSan evidence lifecycle".format(role))
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
    nm_path = shutil.which("nm")
    if nm_path is None:
        raise CampaignError("nm is required for sanitizer binary evidence")
    nm_identity = lab._file_identity(nm_path)
    retained = {}
    primary_error = None
    record = None
    try:
        retained["target"] = _open_retained_probe_file(
            executable, executable_identity,
            "{} sanitizer target".format(role))
        retained["tool"] = _open_retained_probe_file(
            nm_path, nm_identity,
            "{} sanitizer instrumentation tool".format(role))
        completed = _run_bounded_probe_process(
            [retained["target"]["proc_path"],
             INSTRUMENTATION_ARGUMENT],
            environment, "{} sanitizer attestation".format(role), 10.0,
            PROBE_ATTESTATION_STDOUT_LIMIT, PROBE_STDERR_LIMIT,
            pass_fds=(retained["target"]["descriptor"],))
        stdout = completed.stdout.decode(
            "utf-8", errors="replace").strip()
        stderr = completed.stderr.decode(
            "utf-8", errors="replace").strip()
        if (completed.returncode != 0 or stderr or
                len(stdout.splitlines()) != 1):
            raise CampaignError(
                "{} sanitizer attestation failed "
                "(exit {}, stderr={!r})".format(
                    role, completed.returncode, stderr[-512:]))
        try:
            attestation = json.loads(stdout)
        except ValueError as error:
            raise CampaignError(
                "{} sanitizer attestation is not JSON: {}".format(
                    role, error))
        symbols_completed = _run_bounded_probe_process(
            [retained["tool"]["proc_path"], "-g",
             retained["target"]["proc_path"]],
            environment, "{} sanitizer symbol probe".format(role), 30.0,
            PROBE_NM_STDOUT_LIMIT, PROBE_STDERR_LIMIT,
            pass_fds=(
                retained["tool"]["descriptor"],
                retained["target"]["descriptor"]))
        symbols_stderr = symbols_completed.stderr.decode(
            "utf-8", errors="replace").strip()
        symbols_stdout = symbols_completed.stdout.decode(
            "utf-8", errors="replace")
        if symbols_completed.returncode != 0 or symbols_stderr:
            raise CampaignError(
                "{} sanitizer symbol probe failed "
                "(exit {}, stderr={!r})".format(
                    role, symbols_completed.returncode,
                    symbols_stderr[-512:]))
        symbols = sorted(set(
            line.split()[-1] for line in symbols_stdout.splitlines()
            if line.split()))
        missing_symbols = sorted(
            set(INSTRUMENTATION_SYMBOLS) - set(symbols))
        if missing_symbols:
            raise CampaignError(
                "{} target lacks sanitizer symbols: {}".format(
                    role, ", ".join(missing_symbols)))
        record = {
            "executable_sha256": executable_identity["sha256"],
            "attestation": attestation,
            "binary_probe": {
                "schema": "elf-nm-global/v1",
                "tool": nm_identity,
                "required_symbols": list(INSTRUMENTATION_SYMBOLS),
                "symbol_table_digest": lab._digest(symbols),
            },
        }
    except BaseException as error:
        primary_error = error
    finally:
        for retained_record in retained.values():
            try:
                _verify_retained_probe_file(retained_record)
            except BaseException as error:
                primary_error = _merge_probe_exception(
                    primary_error, error,
                    "{} sanitizer retained identity verification failed".format(
                        role))
        for retained_record in retained.values():
            try:
                os.close(retained_record["descriptor"])
            except BaseException as error:
                primary_error = _merge_probe_exception(
                    primary_error, error,
                    "{} sanitizer retained descriptor close failed".format(
                        role))
    if primary_error is not None:
        _raise_probe_exception(
            primary_error, "{} sanitizer evidence lifecycle".format(role))
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
            "probe_execution_policy": _probe_execution_policy(),
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


def _validate_campaign_manifest(value, allow_historical_live_replay=False):
    """Bind an audit to one exact API/pruned sanitizer campaign contract.

    Historical v4 manifests did not sign the probe execution policy.  They are
    accepted only by callers that immediately repeat every live
    instrumentation probe using the current hardened runner.  Offline
    consumers and all newly produced companion evidence require v5.
    """
    manifest = lab.validate_manifest(value)
    source_spec = manifest.get("source_spec", {})
    source_schema = source_spec.get("schema")
    historical = source_schema == HISTORICAL_CAMPAIGN_SCHEMA
    if (source_schema != CAMPAIGN_SCHEMA and
            not (allow_historical_live_replay and historical)):
        raise CampaignError("manifest is not a versioned fuzz campaign")
    metadata = source_spec.get("metadata")
    expected_metadata_fields = {
        "campaign", "iterations_per_seed", "seeds_per_target", "targets",
        "target_executables", "target_instrumentation",
        "sanitizer_environment", "timeout_seconds", "sanitizer_scope",
        "memory_policy",
    }
    if not historical:
        expected_metadata_fields.add("probe_execution_policy")
    if (not isinstance(metadata, dict) or
            set(metadata) != expected_metadata_fields or
            metadata.get("campaign") != CAMPAIGN_NAME or
            metadata.get("targets") != list(TARGETS) or
            metadata.get("sanitizer_environment") != SANITIZER_ENVIRONMENT or
            metadata.get("sanitizer_scope") != SANITIZER_SCOPE):
        raise CampaignError("manifest campaign metadata is not canonical")
    if (not historical and
            metadata.get("probe_execution_policy") !=
                _probe_execution_policy()):
        raise CampaignError(
            "manifest probe execution policy is not canonical")

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


def _campaign_probe_policy_binding(source_schema):
    if source_schema == CAMPAIGN_SCHEMA:
        return "source-manifest-and-live-revalidation"
    if source_schema == HISTORICAL_CAMPAIGN_SCHEMA:
        return "live-revalidation-only-historical-v4"
    raise CampaignError("campaign audit source schema is invalid")


def _validate_campaign_audit(value, source_manifest=None):
    """Validate the closed-world terminal v5 campaign evidence."""
    expected_fields = {
        "schema", "source_campaign_schema", "probe_policy_binding",
        "probe_execution_policy", "manifest_digest", "job_count",
        "distinct_seed_count", "sanitizer_scope",
        "target_instrumentation", "summary", "merged_results",
        "audit_digest",
    }
    if (not isinstance(value, dict) or set(value) != expected_fields or
            value.get("schema") != AUDIT_SCHEMA):
        raise CampaignError("campaign audit is malformed or unsupported")
    unsigned = dict(value)
    audit_digest = unsigned.pop("audit_digest", None)
    if (not isinstance(audit_digest, str) or
            audit_digest != lab._digest(unsigned)):
        raise CampaignError("campaign audit digest does not match its contents")
    source_schema = value.get("source_campaign_schema")
    if (source_schema not in
            (CAMPAIGN_SCHEMA, HISTORICAL_CAMPAIGN_SCHEMA) or
            value.get("probe_policy_binding") !=
                _campaign_probe_policy_binding(source_schema) or
            value.get("probe_execution_policy") !=
                _probe_execution_policy() or
            value.get("sanitizer_scope") != SANITIZER_SCOPE):
        raise CampaignError("campaign audit policy is not canonical")
    if (isinstance(value.get("job_count"), bool) or
            not isinstance(value.get("job_count"), int) or
            value["job_count"] <= 0 or
            isinstance(value.get("distinct_seed_count"), bool) or
            not isinstance(value.get("distinct_seed_count"), int) or
            value["distinct_seed_count"] != value["job_count"]):
        raise CampaignError("campaign audit counts are invalid")
    merged = value.get("merged_results")
    if (not isinstance(merged, dict) or
            value.get("summary") != merged.get("summary") or
            not isinstance(merged.get("results"), list) or
            len(merged["results"]) != value["job_count"] or
            value["summary"] != {
                "missing": 0, "success": value["job_count"]}):
        raise CampaignError("campaign audit merge summary is invalid")
    if source_manifest is not None:
        source_manifest = _validate_campaign_manifest(
            source_manifest,
            allow_historical_live_replay=True)
        metadata = source_manifest["source_spec"]["metadata"]
        if (source_schema != source_manifest["source_spec"]["schema"] or
                value.get("manifest_digest") !=
                    source_manifest["manifest_digest"] or
                value.get("job_count") != len(source_manifest["jobs"]) or
                value.get("target_instrumentation") !=
                    metadata["target_instrumentation"]):
            raise CampaignError(
                "campaign audit differs from its source manifest")
    return value


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
    manifest = _validate_campaign_manifest(
        _load_json(arguments.manifest),
        allow_historical_live_replay=True)
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
    source_schema = manifest["source_spec"]["schema"]
    audit = {
        "schema": AUDIT_SCHEMA,
        "source_campaign_schema": source_schema,
        "probe_policy_binding":
            _campaign_probe_policy_binding(source_schema),
        "probe_execution_policy": _probe_execution_policy(),
        "manifest_digest": manifest["manifest_digest"],
        "job_count": len(manifest["jobs"]),
        "distinct_seed_count": len(seeds),
        "sanitizer_scope": dict(SANITIZER_SCOPE),
        "target_instrumentation": metadata["target_instrumentation"],
        "summary": merged["summary"],
        "merged_results": merged,
    }
    audit["audit_digest"] = lab._digest(audit)
    _validate_campaign_audit(audit, manifest)
    lab._atomic_write_json(destination, audit)
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


def _leak_capture_limits():
    return dict(sorted(LEAK_CAPTURE_LIMITS_BYTES.items()))


def _leak_execution_policy():
    return {
        "mode": "serial",
        "native_max_threads": 1,
        "effective_environment": _leak_environment(),
        "process_count_evidence": False,
        "lsan_helper_processes": "permitted-not-classified",
        "capture_limits_bytes": _leak_capture_limits(),
        "capture_enforcement": dict(LEAK_CAPTURE_ENFORCEMENT),
    }


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
            "capture_limits_bytes": _leak_capture_limits(),
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
        "execution_policy": _leak_execution_policy(),
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
    expected_policy = _leak_execution_policy()
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


def _validate_leak_result(
        result_path, result, job, retained_artifacts=None):
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
    job_dir = None
    if retained_artifacts is None:
        job_dir = Path(result_path).parent
        for name in ("stdout.txt", "stderr.txt", "result.json"):
            _assert_private_regular_file(
                job_dir / name,
                "LSan job {} {}".format(job["id"], name))
    else:
        expected_streams = {"stdout.txt", "stderr.txt", "result.json"}
        if set(retained_artifacts.streams) != expected_streams:
            raise _LeakResumeIdentityError(
                "LSan job {} retained artifact set is incomplete".format(
                    job["id"]))
        retained_artifacts.verify()
        for name in sorted(expected_streams):
            try:
                lab._verify_named_fd(
                    retained_artifacts.fd, name,
                    retained_artifacts.streams[name],
                    "LSan job {} {}".format(job["id"], name),
                    regular=True)
            except lab.LabError as error:
                raise _LeakResumeIdentityError(str(error)) from error
    if (not isinstance(outputs, dict) or set(outputs) != {"stdout", "stderr"}):
        raise CampaignError(
            "leak replay outputs are invalid for job {}".format(job["id"]))
    for name in ("stdout", "stderr"):
        limit = _leak_capture_limit(job, name)
        if retained_artifacts is None:
            capture_path = job_dir / (name + ".txt")
            try:
                capture_size = capture_path.stat().st_size
            except OSError as error:
                raise CampaignError(
                    "cannot inspect leak replay {} for job {}: {}".format(
                        name, job["id"], error))
            capture_identity = lab._content_identity(capture_path)
        else:
            descriptor = retained_artifacts.streams[name + ".txt"]
            try:
                capture_size = os.fstat(descriptor).st_size
                capture_identity = retained_artifacts.content_identity(
                    name + ".txt")
            except (OSError, lab.LabError) as error:
                raise _LeakResumeIdentityError(
                    "cannot inspect retained leak replay {} for job {}: {}"
                    .format(name, job["id"], error)) from error
        if (capture_size > limit or
                not isinstance(outputs.get(name), dict) or
                outputs[name].get("size_bytes") != capture_size or
                capture_identity != outputs[name]):
            raise CampaignError(
                "leak replay {} changed for job {}".format(name, job["id"]))
    if result["outcome"] == "success":
        if retained_artifacts is None:
            stdout = _bounded_leak_capture_from_path(
                job_dir / "stdout.txt",
                _leak_capture_limit(job, "stdout"),
                "LSan job {} stdout".format(job["id"]))
            stderr = _bounded_leak_capture_from_path(
                job_dir / "stderr.txt",
                _leak_capture_limit(job, "stderr"),
                "LSan job {} stderr".format(job["id"]))
        else:
            try:
                stdout = _read_retained_leak_capture(
                    retained_artifacts, "stdout.txt",
                    _leak_capture_limit(job, "stdout"))
                stderr = _read_retained_leak_capture(
                    retained_artifacts, "stderr.txt",
                    _leak_capture_limit(job, "stderr"))
            except (CampaignError, lab.LabError, OSError) as error:
                raise _LeakResumeIdentityError(
                    "cannot read retained leak replay output for job {}: {}"
                    .format(job["id"], error)) from error
        if (result.get("exit_code") != 0 or "detail" in result or
                result.get("executable_before") != job["executable"] or
                result.get("executable_after") != job["executable"] or
                stdout != _expected_leak_stdout(job) or stderr):
            raise CampaignError(
                "successful leak replay evidence is invalid for job {}".format(
                    job["id"]))
    if retained_artifacts is not None:
        retained_artifacts.verify()
    return result


def _write_leak_result(job_dir, result, job_artifacts=None):
    result_path = Path(job_dir) / "result.json"
    if job_artifacts is None:
        for name in ("stdout.txt", "stderr.txt"):
            _assert_private_regular_file(
                Path(job_dir) / name,
                "LSan result {}".format(name))
        _assert_private_regular_file(
            Path(job_dir) / "result.json", "LSan result.json",
            allow_missing=True)
        outputs = {
            "stdout": lab._content_identity(Path(job_dir) / "stdout.txt"),
            "stderr": lab._content_identity(Path(job_dir) / "stderr.txt"),
        }
    else:
        outputs = {
            "stdout": job_artifacts.content_identity("stdout.txt"),
            "stderr": job_artifacts.content_identity("stderr.txt"),
        }
    result["outputs"] = outputs
    unsigned = dict(result)
    unsigned.pop("result_digest", None)
    result["result_digest"] = lab._digest(unsigned)
    try:
        if job_artifacts is None:
            lab._atomic_write_json(result_path, result)
        else:
            job_artifacts.atomic_write_json("result.json", result)
        for name in ("stdout", "stderr"):
            current = (
                lab._content_identity(Path(job_dir) / (name + ".txt"))
                if job_artifacts is None else
                job_artifacts.content_identity(name + ".txt"))
            if current != outputs[name]:
                raise CampaignError(
                    "LSan {} changed during result publication".format(name))
    except BaseException as publication_error:
        invalidation_error = None
        try:
            if job_artifacts is None:
                try:
                    result_path.unlink()
                except FileNotFoundError:
                    pass
            else:
                job_artifacts.invalidate_result()
        except BaseException as error:
            invalidation_error = error
        if invalidation_error is not None:
            note = "LSan result invalidation failed: {}: {}".format(
                type(invalidation_error).__name__, invalidation_error)
            if isinstance(publication_error, (KeyboardInterrupt, SystemExit)):
                if hasattr(publication_error, "add_note"):
                    publication_error.add_note(note)
                raise publication_error
            if isinstance(invalidation_error, (KeyboardInterrupt, SystemExit)):
                if hasattr(invalidation_error, "add_note"):
                    invalidation_error.add_note(
                        "LSan result publication failure: {}: {}".format(
                            type(publication_error).__name__,
                            publication_error))
                raise invalidation_error
            if hasattr(publication_error, "add_note"):
                publication_error.add_note(note)
        raise publication_error
    return result


def _leak_resume_descriptor_state(descriptor, label):
    try:
        value = os.fstat(descriptor)
    except OSError as error:
        raise _LeakResumeIdentityError(
            "cannot inspect retained {}: {}".format(label, error))
    return {
        "descriptor": descriptor,
        "label": label,
        "device": value.st_dev,
        "inode": value.st_ino,
        "mode": value.st_mode,
        "links": value.st_nlink,
        "size": value.st_size,
        "mtime_ns": value.st_mtime_ns,
        "ctime_ns": value.st_ctime_ns,
    }


def _snapshot_leak_resume_owner(owner):
    tree = owner["tree"]
    artifacts = owner["artifacts"]
    records = []
    if tree._edges:
        records.append(_leak_resume_descriptor_state(
            tree._edges[-1][0], "LSan result-root parent"))
    records.extend((
        _leak_resume_descriptor_state(
            tree.root_fd, "LSan result root"),
        _leak_resume_descriptor_state(
            tree.jobs_fd, "LSan jobs directory"),
        _leak_resume_descriptor_state(
            artifacts.fd,
            "LSan job {} directory".format(artifacts.job_id)),
    ))
    for name in sorted(artifacts.streams):
        records.append(_leak_resume_descriptor_state(
            artifacts.streams[name],
            "LSan job {} {}".format(artifacts.job_id, name)))
    owner["snapshots"] = records


def _verify_leak_resume_owner(owner):
    tree = owner["tree"]
    artifacts = owner["artifacts"]
    try:
        artifacts.verify()
        for name in sorted(artifacts.streams):
            lab._verify_named_fd(
                artifacts.fd, name, artifacts.streams[name],
                "LSan job {} {}".format(artifacts.job_id, name),
                regular=True)
    except lab.LabError as error:
        raise _LeakResumeIdentityError(str(error)) from error
    for snapshot in owner.get("snapshots", ()):
        current = _leak_resume_descriptor_state(
            snapshot["descriptor"], snapshot["label"])
        if current != snapshot:
            raise _LeakResumeIdentityError(
                "{} changed during resume validation".format(
                    snapshot["label"]))
    tree.verify()


def _close_leak_resume_owner(owner):
    """Close every retained descriptor while preserving terminal precedence."""
    primary_error = None
    artifacts = owner.get("artifacts")
    tree = owner.get("tree")
    descriptors = []
    if artifacts is not None:
        descriptors.extend(
            (descriptor, "LSan job {} {}".format(artifacts.job_id, name))
            for name, descriptor in artifacts.streams.items())
        artifacts.streams.clear()
        if artifacts.fd is not None:
            descriptors.append((
                artifacts.fd,
                "LSan job {} directory".format(artifacts.job_id)))
            artifacts.fd = None
    if tree is not None:
        descriptors.extend(
            (descriptor, "LSan result-tree descriptor")
            for descriptor in reversed(tree._fds))
        tree._fds = []
        tree.jobs_fd = None
    closed = set()
    for descriptor, label in descriptors:
        if descriptor in closed:
            continue
        closed.add(descriptor)
        close_error = None
        try:
            os.close(descriptor)
        except BaseException as error:
            # close(2) may already have released the numeric descriptor before
            # Python dispatches an asynchronous terminal exception.  The same
            # number can then name an unrelated OFD, even one for the same
            # inode, so retrying by integer can close someone else's resource.
            close_error = error
        if close_error is not None:
            primary_error = _merge_probe_exception(
                primary_error, close_error,
                "{} close failed".format(label))
    owner["artifacts"] = None
    owner["tree"] = None
    owner["snapshots"] = []
    return primary_error


def _open_leak_resume_owner(output_dir, job):
    owner = {"tree": None, "artifacts": None, "snapshots": []}
    primary_error = None
    try:
        owner["tree"] = lab._ResultTree(output_dir, create=True)
        owner["artifacts"] = owner["tree"].open_job(
            job["id"], create=True)
        if owner["artifacts"] is None:
            raise _LeakResumeIdentityError(
                "cannot retain LSan job {} directory".format(job["id"]))
        complete = True
        for name in ("stdout.txt", "stderr.txt", "result.json"):
            descriptor = owner["artifacts"].open_existing(name)
            if descriptor is None:
                complete = False
            else:
                owner["artifacts"].streams[name] = descriptor
        if not complete:
            return owner, False
        _snapshot_leak_resume_owner(owner)
        _verify_leak_resume_owner(owner)
        return owner, True
    except _LeakResumeIdentityError as error:
        primary_error = error
    except (lab.LabError, OSError) as error:
        primary_error = _LeakResumeIdentityError(
            "cannot retain LSan job {} evidence: {}".format(
                job["id"], error))
    except BaseException as error:
        primary_error = error
    cleanup_error = _close_leak_resume_owner(owner)
    if cleanup_error is not None:
        primary_error = _merge_probe_exception(
            primary_error, cleanup_error,
            "LSan resume-owner constructor cleanup failed")
    _raise_probe_exception(primary_error, "LSan resume-owner construction")


def _bounded_leak_result_json(descriptor, job):
    label = "LSan job {} result.json".format(job["id"])
    data = _bounded_leak_capture_from_descriptor(
        descriptor, LEAK_RESULT_JSON_LIMIT_BYTES, label)
    try:
        return json.loads(data.decode("utf-8"))
    except (UnicodeError, ValueError) as error:
        raise CampaignError(
            "invalid JSON in {}: {}".format(label, error))


def _read_resumable_leak_result(output_dir, job):
    owner = None
    eligible = False
    result = None
    primary_error = None
    invalid_result = False
    identity_failure = False
    cleanup_error = None
    try:
        owner, eligible = _open_leak_resume_owner(output_dir, job)
        if eligible:
            _leak_resume_checkpoint("after-retain")
            artifacts = owner["artifacts"]
            result_before = artifacts.content_identity("result.json")
            result = _bounded_leak_result_json(
                artifacts.streams["result.json"], job)
            _leak_resume_checkpoint("after-result-read")
            _validate_leak_result(
                None, result, job, retained_artifacts=artifacts)
            if artifacts.content_identity("result.json") != result_before:
                raise _LeakResumeIdentityError(
                    "LSan job {} result.json changed during validation".format(
                        job["id"]))
            _leak_resume_checkpoint("before-accept")
            _verify_leak_resume_owner(owner)
            if result.get("outcome") != "success":
                invalid_result = True
    except _LeakResumeIdentityError as error:
        primary_error = error
        identity_failure = True
    except lab.LabError as error:
        primary_error = _LeakResumeIdentityError(str(error))
        identity_failure = True
    except (CampaignError, OSError) as error:
        primary_error = error
        invalid_result = True
    except BaseException as error:
        primary_error = error
    finally:
        try:
            _leak_resume_checkpoint("before-close")
        except BaseException as error:
            primary_error = _merge_probe_exception(
                primary_error, error,
                "LSan resume pre-close checkpoint failed")
        if owner is not None:
            if eligible:
                try:
                    _verify_leak_resume_owner(owner)
                except BaseException as error:
                    if not _probe_terminal_exception(error):
                        identity_failure = True
                    if (not _probe_terminal_exception(error) and
                            not isinstance(
                                error, _LeakResumeIdentityError)):
                        error = _LeakResumeIdentityError(str(error))
                    primary_error = _merge_probe_exception(
                        primary_error, error,
                        "LSan resume final identity verification failed")
            cleanup_error = _close_leak_resume_owner(owner)
            if cleanup_error is not None:
                primary_error = _merge_probe_exception(
                    primary_error, cleanup_error,
                    "LSan resume descriptor cleanup failed")
    if _probe_terminal_exception(primary_error):
        raise primary_error
    if cleanup_error is not None:
        _raise_probe_exception(
            primary_error, "LSan resume descriptor cleanup")
    if identity_failure:
        _raise_probe_exception(primary_error, "LSan resume identity")
    if primary_error is not None or invalid_result or not eligible:
        return None
    return result


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


def _leak_capture_limit(job, name):
    limits = job.get("capture_limits_bytes")
    if limits != _leak_capture_limits() or name not in limits:
        raise CampaignError(
            "LSan job {} capture limits are not canonical".format(job["id"]))
    limit = limits[name]
    if isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0:
        raise CampaignError(
            "LSan job {} {} capture limit is invalid".format(
                job["id"], name))
    return limit


def _leak_child_setup(cpu_set, capture_limits):
    """Apply inherited execution limits before the replay executable starts."""
    if resource is None or not hasattr(resource, "RLIMIT_FSIZE"):
        raise RuntimeError("RLIMIT_FSIZE is unavailable")
    if capture_limits != _leak_capture_limits():
        raise RuntimeError("noncanonical LSan capture limits")
    signed_limit = LEAK_CAPTURE_ENFORCEMENT[
        "child_file_size_limit_bytes"]
    unused_soft, inherited_hard = resource.getrlimit(resource.RLIMIT_FSIZE)
    del unused_soft
    if (inherited_hard != resource.RLIM_INFINITY and
            inherited_hard < signed_limit):
        raise RuntimeError(
            "inherited RLIMIT_FSIZE is below the signed capture limit")
    resource.setrlimit(
        resource.RLIMIT_FSIZE, (signed_limit, signed_limit))
    if hasattr(signal, "SIGXFSZ"):
        signal.signal(signal.SIGXFSZ, signal.SIG_DFL)
    lab._child_setup(cpu_set, 0)


def _bounded_leak_capture_from_descriptor(
        descriptor, maximum_bytes, label):
    """Read at most maximum_bytes + 1 while detecting concurrent growth."""
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            raise CampaignError(
                "{} must be a private regular file".format(label))
        if before.st_size > maximum_bytes:
            raise CampaignError(
                "{} exceeds its signed {}-byte limit".format(
                    label, maximum_bytes))
        data = bytearray()
        offset = 0
        while True:
            request = min(1024 * 1024, maximum_bytes + 1 - offset)
            if request <= 0:
                raise CampaignError(
                    "{} exceeds its signed {}-byte limit".format(
                        label, maximum_bytes))
            block = os.pread(descriptor, request, offset)
            if not block:
                break
            data.extend(block)
            offset += len(block)
        after = os.fstat(descriptor)
    except CampaignError:
        raise
    except OSError as error:
        raise CampaignError("cannot read {}: {}".format(label, error))
    if (before.st_dev != after.st_dev or before.st_ino != after.st_ino or
            before.st_size != after.st_size or
            before.st_mtime_ns != after.st_mtime_ns or
            before.st_ctime_ns != after.st_ctime_ns):
        raise CampaignError("{} changed while it was read".format(label))
    return bytes(data)


def _bounded_leak_capture_from_path(path, maximum_bytes, label):
    flags = os.O_RDONLY
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(str(path), flags)
    except OSError as error:
        raise CampaignError("cannot open {}: {}".format(label, error))
    try:
        return _bounded_leak_capture_from_descriptor(
            descriptor, maximum_bytes, label)
    finally:
        os.close(descriptor)


def _retained_leak_capture_size(job_artifacts, name, maximum_bytes):
    descriptor = job_artifacts.streams.get(name)
    if descriptor is None:
        raise CampaignError(
            "LSan job {} {} is not retained".format(
                job_artifacts.job_id, name))
    try:
        size = os.fstat(descriptor).st_size
    except OSError as error:
        raise CampaignError(
            "cannot inspect LSan job {} {}: {}".format(
                job_artifacts.job_id, name, error))
    if size > maximum_bytes:
        raise CampaignError(
            "LSan job {} {} exceeds its signed {}-byte limit".format(
                job_artifacts.job_id, name, maximum_bytes))
    return size


def _read_retained_leak_capture(
        job_artifacts, name, maximum_bytes):
    """Read a retained capture inode and reject concurrent replacement/change."""
    descriptor = job_artifacts.streams.get(name)
    if descriptor is None:
        raise CampaignError(
            "LSan job {} {} is not retained".format(
                job_artifacts.job_id, name))
    _retained_leak_capture_size(
        job_artifacts, name, maximum_bytes)
    before = job_artifacts.content_identity(name)
    data = _bounded_leak_capture_from_descriptor(
        descriptor, maximum_bytes,
        "LSan job {} {}".format(job_artifacts.job_id, name))
    after = job_artifacts.content_identity(name)
    if before != after:
        raise CampaignError(
            "LSan job {} {} changed while it was read".format(
                job_artifacts.job_id, name))
    return data


def _append_leak_detail(detail, addition):
    return addition if not detail else detail + "; " + addition


def _execute_leak_job(manifest, job, output_dir):
    capture_limits = {
        name: _leak_capture_limit(job, name)
        for name in ("stdout", "stderr")
    }
    job_dir = _leak_job_directory(output_dir, job["id"])
    started = time.monotonic()
    before = None
    after = None
    exit_code = None
    outcome = "launch_error"
    detail = None
    process = None
    stdout_handle = None
    stderr_handle = None
    result_tree = None
    job_artifacts = None
    cleanup = None
    primary_error = None
    containment_token = hashlib.sha256(
        os.urandom(32) + job["id"].encode("utf-8")).hexdigest()
    active = _new_lsan_containment(None, containment_token)
    try:
        result_tree = lab._ResultTree(output_dir, create=True)
        job_artifacts = result_tree.open_job(job["id"], create=True)
        if job_artifacts is None:
            raise CampaignError(
                "cannot create LSan job {} output".format(job["id"]))
        # Preserve the established fail-closed contract for pre-existing
        # symlink, FIFO, or hard-link artifacts before publishing replacements.
        for name in ("stdout.txt", "stderr.txt", "result.json"):
            existing = job_artifacts.open_existing(name)
            if existing is not None:
                os.close(existing)
        job_artifacts.invalidate_result()
        stdout_descriptor = job_artifacts.open_capture("stdout.txt")
        stderr_descriptor = job_artifacts.open_capture("stderr.txt")
        stdout_handle = os.fdopen(
            os.dup(stdout_descriptor), "wb", buffering=0)
        stderr_handle = os.fdopen(
            os.dup(stderr_descriptor), "wb", buffering=0)
        before = lab._file_identity(job["executable"]["path"])
        if before != job["executable"]:
            raise CampaignError("executable changed before leak replay")
        cwd = Path(job["cwd"])
        if not cwd.is_absolute():
            cwd = Path(manifest["source_manifest"]["root"]) / cwd
        environment = os.environ.copy()
        environment.update(job["environment"])
        environment["LEO2_LAB_CONTAINMENT_TOKEN"] = containment_token
        with lab._ChildSubreaperLease():
            _begin_lsan_direct_child_capture(active)
            try:
                process = subprocess.Popen(
                    job["command"], cwd=str(cwd), env=environment,
                    stdin=subprocess.DEVNULL, stdout=stdout_handle,
                    stderr=stderr_handle, start_new_session=True,
                    preexec_fn=lambda: _leak_child_setup(
                        job["cpu_set"], capture_limits))
                active["process"] = process
                exit_code = _wait_lsan_leader(
                    active, job["timeout_seconds"])
                if exit_code is None:
                    outcome = "timeout"
                    detail = "leak replay exceeded its signed timeout"
                elif (hasattr(signal, "SIGXFSZ") and
                        exit_code == -int(signal.SIGXFSZ)):
                    outcome = "evidence_invalid"
                    detail = (
                        "leak replay exceeded its signed capture byte limit "
                        "(SIGXFSZ)")
                else:
                    outcome = "success" if exit_code == 0 else "failed"
                    if outcome == "failed":
                        detail = "leak replay exited {}".format(exit_code)
            except BaseException as error:
                primary_error = error
            finally:
                try:
                    cleanup = _cleanup_lsan_process_tree_with_retry(active)
                except BaseException as cleanup_error:
                    primary_error = _merge_probe_exception(
                        primary_error, cleanup_error,
                        "LSan replay cleanup failed")

        if primary_error is not None:
            if isinstance(primary_error, (
                    OSError, subprocess.SubprocessError,
                    lab.LabError, CampaignError)):
                detail = str(primary_error)
                outcome = (
                    "launch_error" if process is None
                    else "evidence_invalid")
            else:
                raise primary_error
        if process is not None and process.returncode is not None:
            exit_code = process.returncode
        if (exit_code is not None and hasattr(signal, "SIGXFSZ") and
                exit_code == -int(signal.SIGXFSZ) and
                "(SIGXFSZ)" not in (detail or "")):
            outcome = "evidence_invalid"
            detail = _append_leak_detail(
                detail,
                "leak replay exceeded its signed capture byte limit "
                "(SIGXFSZ)")
        if cleanup is None:
            raise CampaignError("LSan replay cleanup did not complete")
        if cleanup["leader_remaining"]:
            outcome = "evidence_invalid"
            detail = _append_leak_detail(
                detail,
                "job leader remained live after bounded SIGKILL cleanup")
        if cleanup["errors"]:
            outcome = "evidence_invalid"
            detail = _append_leak_detail(
                detail, "; ".join(cleanup["errors"]))
        if cleanup["detected"]:
            residual_detail = (
                "job leader exited with residual descendants {}; cleanup "
                "remaining {}".format(
                    ",".join(str(pid) for pid in cleanup["detected"]),
                    ",".join(str(pid) for pid in cleanup["remaining"]) or
                    "none"))
            if outcome in ("success", "failed"):
                outcome = "evidence_invalid"
            detail = _append_leak_detail(detail, residual_detail)

        if stdout_handle is not None:
            stdout_handle.close()
            stdout_handle = None
        if stderr_handle is not None:
            stderr_handle.close()
            stderr_handle = None
        capture_limit_errors = []
        capture_sizes = {}
        for name in ("stdout", "stderr"):
            try:
                capture_sizes[name] = _retained_leak_capture_size(
                    job_artifacts, name + ".txt", capture_limits[name])
            except CampaignError as error:
                capture_limit_errors.append(str(error))
        if capture_limit_errors:
            outcome = "evidence_invalid"
            detail = _append_leak_detail(
                detail, "; ".join(capture_limit_errors))
        saturated = [
            name for name, size in capture_sizes.items()
            if size == LEAK_CAPTURE_ENFORCEMENT[
                "child_file_size_limit_bytes"]
        ]
        if saturated and "capture byte limit" not in (detail or ""):
            # A process may ignore SIGXFSZ and observe EFBIG instead.  An exact
            # hard-limit capture cannot prove that its final diagnostic was
            # retained, so it is invalid even when the leader exits normally.
            outcome = "evidence_invalid"
            detail = _append_leak_detail(
                detail,
                "{} capture reached its signed capture byte limit "
                "(SIGXFSZ/EFBIG)".format(",".join(saturated)))
        try:
            after = lab._file_identity(job["executable"]["path"])
        except lab.LabError as error:
            after = None
            outcome = "evidence_invalid"
            detail = _append_leak_detail(detail, str(error))
        if before != job["executable"] or after != job["executable"]:
            outcome = "evidence_invalid"
            detail = _append_leak_detail(
                detail, "executable identity changed during leak replay")
        if outcome == "success":
            try:
                stdout = _read_retained_leak_capture(
                    job_artifacts, "stdout.txt",
                    capture_limits["stdout"])
                stderr = _read_retained_leak_capture(
                    job_artifacts, "stderr.txt",
                    capture_limits["stderr"])
            except (CampaignError, lab.LabError, OSError) as error:
                outcome = "evidence_invalid"
                detail = _append_leak_detail(
                    detail,
                    "cannot read leak replay output: {}".format(error))
            else:
                if stderr:
                    outcome = "diagnostic"
                    detail = _append_leak_detail(
                        detail,
                        "leak replay emitted sanitizer or unexpected stderr")
                elif stdout != _expected_leak_stdout(job):
                    outcome = "evidence_invalid"
                    detail = _append_leak_detail(
                        detail,
                        "leak replay success marker is missing or mismatched")

        result = {
            "schema": LEAK_RESULT_SCHEMA,
            "state": "complete",
            "job_id": job["id"],
            "job_digest": job["job_digest"],
            "outcome": outcome,
            "exit_code": exit_code,
            "duration_seconds": round(
                max(0.0, time.monotonic() - started), 6),
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
        return _write_leak_result(
            job_dir, result, job_artifacts=job_artifacts)
    finally:
        if stdout_handle is not None:
            stdout_handle.close()
        if stderr_handle is not None:
            stderr_handle.close()
        if job_artifacts is not None:
            job_artifacts.close()
        if result_tree is not None:
            result_tree.close()


def _verify_leak_runtime_cpus(manifest):
    if not hasattr(os, "sched_getaffinity"):
        raise CampaignError("LSan companion execution requires Linux affinity")
    if resource is None or not hasattr(resource, "RLIMIT_FSIZE"):
        raise CampaignError(
            "LSan companion execution requires RLIMIT_FSIZE")
    unused_soft, hard_limit = resource.getrlimit(resource.RLIMIT_FSIZE)
    del unused_soft
    signed_limit = LEAK_CAPTURE_ENFORCEMENT[
        "child_file_size_limit_bytes"]
    if (hard_limit != resource.RLIM_INFINITY and
            hard_limit < signed_limit):
        raise CampaignError(
            "inherited RLIMIT_FSIZE is below the signed capture limit")
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


def _probe_runner_self_test(directory, nm_path, nm_identity):
    """Exercise bounded probe containment and terminal-error precedence."""
    if (not sys.platform.startswith("linux") or not hasattr(os, "fork") or
            not hasattr(os, "setsid") or
            not Path("/proc/self/stat").is_file()):
        return

    helper = Path(directory) / "fake-bounded-probe.py"
    helper.write_text(
        "#!{}\n"
        "import os, resource, sys, time\n"
        "mode = sys.argv[1]\n"
        "pid_path = sys.argv[2] if len(sys.argv) > 2 else ''\n"
        "if mode == 'normal':\n"
        "    if pid_path:\n"
        "        open(pid_path, 'w').write(str(os.getpid()))\n"
        "    sys.stdout.write('probe stdout\\n')\n"
        "    sys.stderr.write('probe stderr\\n')\n"
        "    raise SystemExit(7)\n"
        "if mode == 'detached-writer':\n"
        "    pid = os.fork()\n"
        "    if pid == 0:\n"
        "        os.setsid()\n"
        "        open(pid_path, 'w').write(str(os.getpid()))\n"
        "        while True:\n"
        "            os.write(1, b'detached writer\\n')\n"
        "            time.sleep(0.02)\n"
        "    while not os.path.exists(pid_path):\n"
        "        time.sleep(0.005)\n"
        "    os._exit(0)\n"
        "if mode == 'excessive':\n"
        "    open(pid_path, 'w').write(str(os.getpid()))\n"
        "    while True:\n"
        "        os.write(1, b'x' * 4096)\n"
        "if mode == 'raise-file-limit':\n"
        "    open(pid_path, 'w').write(str(os.getpid()))\n"
        "    expected = int(sys.argv[3])\n"
        "    soft, hard = resource.getrlimit(resource.RLIMIT_FSIZE)\n"
        "    if soft != expected or hard != expected:\n"
        "        sys.stderr.write('unexpected file limit {{}} {{}}\\n'.format("
        "soft, hard))\n"
        "        raise SystemExit(91)\n"
        "    try:\n"
        "        resource.setrlimit(resource.RLIMIT_FSIZE, "
        "(expected + 1, expected + 1))\n"
        "    except (OSError, ValueError):\n"
        "        pass\n"
        "    else:\n"
        "        raise SystemExit(92)\n"
        "    while True:\n"
        "        os.write(1, b'r' * 4096)\n"
        "if mode == 'mutate':\n"
        "    with open(sys.argv[0], 'ab') as output:\n"
        "        output.write(b'\\n# changed during probe\\n')\n"
        "    raise SystemExit(0)\n"
        "open(pid_path, 'w').write(str(os.getpid()))\n"
        "time.sleep(30)\n".format(sys.executable),
        encoding="utf-8")
    helper.chmod(0o700)
    helper_identity = lab._file_identity(helper)

    def run_helper(mode, pid_path, timeout_seconds=1.0,
                   stdout_limit=4096, stderr_limit=4096):
        return _run_retained_probe(
            [_retained_probe_argument("target"), mode, str(pid_path),
             str(max(stdout_limit, stderr_limit) + 1)],
            {"target": (
                str(helper), helper_identity, "modeled probe target")},
            {}, "modeled {} probe".format(mode), timeout_seconds,
            stdout_limit, stderr_limit)

    def require_process_gone(pid_path, label):
        try:
            pid = int(pid_path.read_text(encoding="ascii"))
        except (OSError, ValueError) as error:
            raise CampaignError(
                "{} did not publish a process identity: {}".format(
                    label, error))
        deadline = time.monotonic() + 1.0
        while (Path("/proc") / str(pid)).exists() and \
                time.monotonic() < deadline:
            time.sleep(0.01)
        if (Path("/proc") / str(pid)).exists():
            raise CampaignError(
                "{} left process {} live".format(label, pid))

    ordinary_pid = Path(directory) / "ordinary-probe.pid"
    ordinary = run_helper("normal", ordinary_pid)
    if (ordinary.returncode != 7 or
            ordinary.stdout != b"probe stdout\n" or
            ordinary.stderr != b"probe stderr\n"):
        raise CampaignError(
            "bounded probe changed ordinary captures: returncode={}, "
            "stdout={!r}, stderr={!r}".format(
                ordinary.returncode, ordinary.stdout, ordinary.stderr))
    require_process_gone(ordinary_pid, "ordinary bounded probe")

    def probe_resource_snapshot():
        return (
            len(os.listdir("/proc/self/fd")),
            {
                (record["pid"], record["starttime_ticks"])
                for record in _current_lsan_direct_children()
            })

    def expect_capture_mutation(label, phase, restore_name):
        original_checkpoint = globals()["_probe_capture_checkpoint"]
        mutation_count = []
        pid_path = Path(directory) / (
            "{}-capture-probe.pid".format(label))
        before_resources = probe_resource_snapshot()

        def mutate_capture(observed_phase, directory_descriptor=None):
            if observed_phase != phase or mutation_count:
                return
            mutation_count.append(observed_phase)
            os.rename(
                "stdout.txt", ".stdout.saved",
                src_dir_fd=directory_descriptor,
                dst_dir_fd=directory_descriptor)
            replacement = os.open(
                "stdout.txt",
                os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                os.O_CLOEXEC | os.O_NOFOLLOW,
                0o600, dir_fd=directory_descriptor)
            try:
                os.write(replacement, b"substituted probe capture\n")
            finally:
                os.close(replacement)
            if restore_name:
                os.unlink("stdout.txt", dir_fd=directory_descriptor)
                os.rename(
                    ".stdout.saved", "stdout.txt",
                    src_dir_fd=directory_descriptor,
                    dst_dir_fd=directory_descriptor)

        globals()["_probe_capture_checkpoint"] = mutate_capture
        try:
            try:
                run_helper("normal", pid_path)
            except CampaignError:
                pass
            else:
                raise CampaignError(
                    "bounded probe accepted {}".format(label))
        finally:
            globals()["_probe_capture_checkpoint"] = original_checkpoint
        if mutation_count != [phase]:
            raise CampaignError(
                "bounded probe did not exercise {}".format(label))
        require_process_gone(pid_path, label)
        if probe_resource_snapshot() != before_resources:
            raise CampaignError(
                "bounded probe leaked resources after {}".format(label))

    expect_capture_mutation(
        "capture-rename-replacement", "after-cleanup", False)
    expect_capture_mutation(
        "capture-name-ABA", "after-read", True)

    original_capture_checkpoint = globals()["_probe_capture_checkpoint"]
    for label, terminal, expected_type, expected_code in (
            ("capture-KeyboardInterrupt", KeyboardInterrupt(),
             KeyboardInterrupt, None),
            ("capture-SystemExit", SystemExit(47), SystemExit, 47)):
        pid_path = Path(directory) / (label + ".pid")
        before_resources = probe_resource_snapshot()
        interruptions = []

        def interrupt_capture(
                observed_phase, unused_directory_descriptor=None,
                exception=terminal):
            del unused_directory_descriptor
            if observed_phase == "after-cleanup" and not interruptions:
                interruptions.append(observed_phase)
                raise exception

        globals()["_probe_capture_checkpoint"] = interrupt_capture
        try:
            try:
                run_helper("normal", pid_path)
            except BaseException as observed:
                if not isinstance(observed, expected_type):
                    raise CampaignError(
                        "bounded capture propagated {} instead of {}"
                        .format(type(observed).__name__, label))
                if (expected_code is not None and
                        getattr(observed, "code", None) != expected_code):
                    raise CampaignError(
                        "bounded capture changed {} code".format(label))
            else:
                raise CampaignError(
                    "bounded capture swallowed {}".format(label))
        finally:
            globals()["_probe_capture_checkpoint"] = (
                original_capture_checkpoint)
        if interruptions != ["after-cleanup"]:
            raise CampaignError(
                "bounded capture did not inject {} exactly once".format(
                    label))
        require_process_gone(pid_path, label)
        if probe_resource_snapshot() != before_resources:
            raise CampaignError(
                "bounded capture leaked resources after {}".format(label))

    tool_version = _run_retained_probe(
        [_retained_probe_argument("tool"), "--version"],
        {"tool": (
            nm_path, nm_identity, "modeled instrumentation tool")},
        {}, "ordinary instrumentation tool probe", 5.0,
        1024 * 1024, PROBE_STDERR_LIMIT)
    if tool_version.returncode != 0 or not tool_version.stdout:
        raise CampaignError(
            "bounded instrumentation tool probe did not complete normally")
    python_identity = lab._file_identity(sys.executable)
    symbol_probe = _run_retained_probe(
        [_retained_probe_argument("tool"), "-g",
         _retained_probe_argument("target")],
        {
            "tool": (
                nm_path, nm_identity, "modeled instrumentation tool"),
            "target": (
                sys.executable, python_identity,
                "modeled instrumentation target"),
        },
        {}, "ordinary retained symbol probe", 5.0,
        PROBE_NM_STDOUT_LIMIT, PROBE_STDERR_LIMIT)
    if symbol_probe.returncode != 0:
        raise CampaignError(
            "bounded retained symbol probe did not complete normally: "
            "returncode={}, stderr={!r}".format(
                symbol_probe.returncode, symbol_probe.stderr[-512:]))

    detached_pid = Path(directory) / "detached-probe.pid"
    try:
        run_helper("detached-writer", detached_pid)
    except CampaignError as error:
        if "residual descendants" not in str(error):
            raise
    else:
        raise CampaignError(
            "bounded probe accepted a detached output writer")
    require_process_gone(detached_pid, "detached bounded probe")

    timeout_pid = Path(directory) / "timeout-probe.pid"
    timeout_started = time.monotonic()
    try:
        run_helper("sleep", timeout_pid, timeout_seconds=0.05)
    except CampaignError as error:
        if "timed out" not in str(error):
            raise
    else:
        raise CampaignError("bounded probe timeout was accepted")
    if time.monotonic() - timeout_started > 3.0:
        raise CampaignError("bounded probe timeout teardown was unbounded")
    require_process_gone(timeout_pid, "timed-out bounded probe")

    excessive_pid = Path(directory) / "excessive-probe.pid"
    try:
        run_helper(
            "excessive", excessive_pid, timeout_seconds=1.0,
            stdout_limit=1024, stderr_limit=1024)
    except CampaignError as error:
        if "excessive output" not in str(error):
            raise
    else:
        raise CampaignError("bounded probe accepted excessive output")
    require_process_gone(excessive_pid, "excessive-output bounded probe")

    raised_limit_pid = Path(directory) / "raised-limit-probe.pid"
    try:
        run_helper(
            "raise-file-limit", raised_limit_pid, timeout_seconds=1.0,
            stdout_limit=1024, stderr_limit=1024)
    except CampaignError as error:
        if "excessive output" not in str(error):
            raise
    else:
        raise CampaignError(
            "bounded probe could raise its inherited file-size hard limit")
    require_process_gone(
        raised_limit_pid, "raised-file-limit bounded probe")

    if hasattr(signal, "setitimer") and hasattr(signal, "ITIMER_REAL"):
        for label, exception_factory in (
                ("KeyboardInterrupt", KeyboardInterrupt),
                ("SystemExit", lambda: SystemExit(23))):
            pid_path = Path(directory) / (
                "{}-probe.pid".format(label.lower()))
            prior_handler = signal.getsignal(signal.SIGALRM)

            def interrupt_probe(_signum, _frame, factory=exception_factory):
                raise factory()

            signal.signal(signal.SIGALRM, interrupt_probe)
            signal.setitimer(signal.ITIMER_REAL, 0.15)
            try:
                try:
                    run_helper("sleep", pid_path, timeout_seconds=5.0)
                except BaseException as error:
                    if label == "KeyboardInterrupt":
                        if not isinstance(error, KeyboardInterrupt):
                            raise
                    elif not (isinstance(error, SystemExit) and
                              error.code == 23):
                        raise
                else:
                    raise CampaignError(
                        "bounded probe swallowed {}".format(label))
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0.0)
                signal.signal(signal.SIGALRM, prior_handler)
            require_process_gone(
                pid_path, "{} bounded probe".format(label))

    original_popen = subprocess.Popen
    original_cleanup = globals()["_cleanup_lsan_process_tree_with_retry"]

    def cleanup_then(error):
        def scripted(active):
            original_cleanup(active)
            raise error
        return scripted

    precedence_cases = (
        (CampaignError("ordinary launch"), KeyboardInterrupt(),
         KeyboardInterrupt, None),
        (KeyboardInterrupt(), SystemExit(31), KeyboardInterrupt, None),
        (SystemExit(29), KeyboardInterrupt(), SystemExit, 29),
    )
    for index, (launch_error, cleanup_error, expected_type,
                expected_code) in enumerate(precedence_cases):
        def fail_launch(*unused_args, error=launch_error, **unused_kwargs):
            del unused_args, unused_kwargs
            raise error

        subprocess.Popen = fail_launch
        globals()["_cleanup_lsan_process_tree_with_retry"] = cleanup_then(
            cleanup_error)
        try:
            try:
                _run_bounded_probe_process(
                    [str(helper), "normal", ""], {},
                    "modeled precedence probe {}".format(index), 1.0,
                    4096, 4096)
            except BaseException as observed:
                if not isinstance(observed, expected_type):
                    raise CampaignError(
                        "probe precedence case {} propagated {} "
                        "instead of {}".format(
                            index, type(observed).__name__,
                            expected_type.__name__))
                if (expected_code is not None and
                        getattr(observed, "code", None) != expected_code):
                    raise CampaignError(
                        "probe precedence case {} changed SystemExit code"
                        .format(index))
            else:
                raise CampaignError(
                    "probe precedence case {} unexpectedly passed".format(
                        index))
        finally:
            subprocess.Popen = original_popen
            globals()["_cleanup_lsan_process_tree_with_retry"] = (
                original_cleanup)

    mutating = Path(directory) / "mutating-bounded-probe.py"
    shutil.copy2(helper, mutating)
    mutating.chmod(0o700)
    mutating_identity = lab._file_identity(mutating)
    try:
        _run_retained_probe(
            [_retained_probe_argument("target"), "mutate", ""],
            {"target": (
                str(mutating), mutating_identity,
                "mutating modeled probe target")},
            {}, "mutating modeled probe", 1.0, 4096, 4096)
    except CampaignError as error:
        if "changed during the complete probe lifecycle" not in str(error):
            raise
    else:
        raise CampaignError(
            "retained probe accepted an executable identity mutation")


def _leak_companion_self_test(directory, nm_identity):
    def verify_cleanup_exception_precedence(
            first_error, retry_error, expected_error, label):
        original_cleanup = globals()["_cleanup_lsan_process_tree"]
        failures = [first_error, retry_error]

        def scripted_cleanup(unused_active):
            del unused_active
            raise failures.pop(0)

        globals()["_cleanup_lsan_process_tree"] = scripted_cleanup
        try:
            try:
                _cleanup_lsan_process_tree_with_retry({})
            except BaseException as observed:
                if observed is not expected_error:
                    raise CampaignError(
                        "{} propagated {} instead of {}".format(
                            label, type(observed).__name__,
                            type(expected_error).__name__))
            else:
                raise CampaignError(
                    "{} unexpectedly completed cleanup".format(label))
        finally:
            globals()["_cleanup_lsan_process_tree"] = original_cleanup
        if failures:
            raise CampaignError(
                "{} did not execute the cleanup retry".format(label))

    ordinary_then_interrupt = (
        CampaignError("initial ordinary cleanup failure"),
        KeyboardInterrupt())
    verify_cleanup_exception_precedence(
        ordinary_then_interrupt[0], ordinary_then_interrupt[1],
        ordinary_then_interrupt[1],
        "ordinary-to-KeyboardInterrupt cleanup")
    ordinary_then_exit = (
        CampaignError("initial ordinary cleanup failure"), SystemExit(17))
    verify_cleanup_exception_precedence(
        ordinary_then_exit[0], ordinary_then_exit[1],
        ordinary_then_exit[1], "ordinary-to-SystemExit cleanup")
    interrupt_then_ordinary = (
        KeyboardInterrupt(), CampaignError("retry ordinary cleanup failure"))
    verify_cleanup_exception_precedence(
        interrupt_then_ordinary[0], interrupt_then_ordinary[1],
        interrupt_then_ordinary[0],
        "KeyboardInterrupt-to-ordinary cleanup")
    exit_then_interrupt = (SystemExit(19), KeyboardInterrupt())
    verify_cleanup_exception_precedence(
        exit_then_interrupt[0], exit_then_interrupt[1],
        exit_then_interrupt[0], "SystemExit-to-KeyboardInterrupt cleanup")

    class RecycledCloseArtifacts(object):
        def __init__(self, descriptor):
            self.job_id = "modeled-recycled-close"
            self.streams = {"stdout.txt": descriptor}
            self.fd = None

    recycled_path = Path(directory) / "same-inode-close-recycle.bin"
    recycled_path.write_bytes(b"same inode, different open file description")
    descriptor_count_before = len(os.listdir("/proc/self/fd"))
    owned_descriptor = os.open(
        str(recycled_path), os.O_RDONLY | os.O_CLOEXEC)
    recycled_artifacts = RecycledCloseArtifacts(owned_descriptor)
    recycled_owner = {
        "tree": None,
        "artifacts": recycled_artifacts,
        "snapshots": [],
    }
    original_os_close = os.close
    recycled = {"descriptor": None, "calls": 0}

    def close_then_recycle(descriptor):
        if descriptor == owned_descriptor and recycled["calls"] == 0:
            recycled["calls"] += 1
            original_os_close(descriptor)
            replacement = os.open(
                str(recycled_path), os.O_RDONLY | os.O_CLOEXEC)
            recycled["descriptor"] = replacement
            if replacement != descriptor:
                raise CampaignError(
                    "same-inode close test did not recycle the numeric fd")
            raise KeyboardInterrupt()
        recycled["calls"] += 1
        original_os_close(descriptor)

    os.close = close_then_recycle
    try:
        close_error = _close_leak_resume_owner(recycled_owner)
    finally:
        os.close = original_os_close
    try:
        if (not isinstance(close_error, KeyboardInterrupt) or
                recycled["calls"] != 1 or
                recycled["descriptor"] != owned_descriptor):
            raise CampaignError(
                "resume owner retried an ambiguously closed numeric fd")
        retained_state = os.fstat(recycled["descriptor"])
        path_state = recycled_path.stat()
        if (retained_state.st_dev != path_state.st_dev or
                retained_state.st_ino != path_state.st_ino):
            raise CampaignError(
                "same-inode recycled descriptor changed identity")
    finally:
        if recycled["descriptor"] is not None:
            original_os_close(recycled["descriptor"])
    if len(os.listdir("/proc/self/fd")) != descriptor_count_before:
        raise CampaignError(
            "same-inode recycled descriptor test leaked an fd")

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

    class PublicationRaceArtifacts(object):
        def __init__(self, interrupt=False):
            self.changed = False
            self.interrupt = interrupt
            self.invalidated = False

        def content_identity(self, name):
            del name
            return {
                "sha256": ("b" if self.changed else "a") * 64,
                "size_bytes": 1,
            }

        def atomic_write_json(self, name, value):
            del name, value
            self.changed = True
            if self.interrupt:
                raise KeyboardInterrupt()

        def invalidate_result(self):
            self.invalidated = True

    publication_race = PublicationRaceArtifacts()
    try:
        _write_leak_result(
            Path(directory) / "unused-publication-race", {
                "schema": LEAK_RESULT_SCHEMA,
                "state": "complete",
            }, job_artifacts=publication_race)
    except CampaignError:
        pass
    else:
        raise CampaignError(
            "LSan result publication accepted changed capture bytes")
    if not publication_race.invalidated:
        raise CampaignError(
            "LSan result publication did not invalidate changed evidence")

    interrupted_publication = PublicationRaceArtifacts(interrupt=True)
    try:
        _write_leak_result(
            Path(directory) / "unused-interrupted-publication", {
                "schema": LEAK_RESULT_SCHEMA,
                "state": "complete",
            }, job_artifacts=interrupted_publication)
    except KeyboardInterrupt:
        pass
    else:
        raise CampaignError(
            "LSan result publication swallowed KeyboardInterrupt")
    if not interrupted_publication.invalidated:
        raise CampaignError(
            "interrupted LSan result publication retained stale authority")

    if (sys.platform.startswith("linux") and hasattr(os, "fork") and
            hasattr(os, "setsid") and Path("/proc/self/stat").is_file()):
        canary_helper = Path(directory) / "fake-lsan-canary-process.py"
        canary_helper.write_text(
            "#!{}\n"
            "import os, sys, time\n"
            "mode = os.environ.get('LEO2_TEST_CANARY_MODE', 'normal')\n"
            "pid_path = os.environ.get('LEO2_TEST_CANARY_PID')\n"
            "if mode == 'normal':\n"
            "    sys.stdout.write('canary stdout\\n')\n"
            "    sys.stderr.write('canary stderr\\n')\n"
            "    raise SystemExit(7)\n"
            "if mode == 'detached':\n"
            "    pid = os.fork()\n"
            "    if pid == 0:\n"
            "        os.setsid()\n"
            "        with open(pid_path, 'w') as output:\n"
            "            output.write(str(os.getpid()))\n"
            "        time.sleep(30)\n"
            "        os._exit(0)\n"
            "    while not os.path.exists(pid_path):\n"
            "        time.sleep(0.005)\n"
            "    os._exit(0)\n"
            "with open(pid_path, 'w') as output:\n"
            "    output.write(str(os.getpid()))\n"
            "time.sleep(30)\n".format(sys.executable),
            encoding="utf-8")
        canary_helper.chmod(0o700)

        normal = _run_lsan_canary_process(
            str(canary_helper), {
                "LEO2_TEST_CANARY_MODE": "normal",
            }, "modeled", timeout_seconds=1.0)
        if (normal.returncode != 7 or
                normal.stdout != b"canary stdout\n" or
                normal.stderr != b"canary stderr\n"):
            raise CampaignError(
                "bounded LSan canary runner changed ordinary captures")

        def require_process_gone(pid_path, label):
            try:
                pid = int(pid_path.read_text(encoding="ascii"))
            except (OSError, ValueError) as error:
                raise CampaignError(
                    "{} did not publish a process identity: {}".format(
                        label, error))
            deadline = time.monotonic() + 1.0
            while (Path("/proc") / str(pid)).exists() and \
                    time.monotonic() < deadline:
                time.sleep(0.01)
            if (Path("/proc") / str(pid)).exists():
                raise CampaignError(
                    "{} left process {} live".format(label, pid))

        tokenless_pid = Path(directory) / "tokenless-canary.pid"
        original_popen = subprocess.Popen
        fault_child = {"pid": None}

        def interrupt_after_fork(*unused_args, **unused_kwargs):
            del unused_args, unused_kwargs
            pid = os.fork()
            if pid == 0:
                try:
                    os.setsid()
                    time.sleep(30)
                finally:
                    os._exit(0)
            fault_child["pid"] = pid
            tokenless_pid.write_text(str(pid), encoding="ascii")
            raise KeyboardInterrupt()

        subprocess.Popen = interrupt_after_fork
        try:
            try:
                _run_lsan_canary_process(
                    str(canary_helper), {}, "modeled-tokenless",
                    timeout_seconds=1.0)
            except KeyboardInterrupt:
                pass
            else:
                raise CampaignError(
                    "LSan canary swallowed a post-fork KeyboardInterrupt")
        finally:
            subprocess.Popen = original_popen
        try:
            require_process_gone(
                tokenless_pid, "tokenless unpublished LSan child cleanup")
        except BaseException:
            pid = fault_child["pid"]
            if pid is not None:
                try:
                    waited, unused_status = os.waitpid(pid, os.WNOHANG)
                    del unused_status
                except ChildProcessError:
                    waited = pid
                if waited == 0:
                    try:
                        os.kill(pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    try:
                        os.waitpid(pid, 0)
                    except ChildProcessError:
                        pass
            raise

        detached_pid = Path(directory) / "detached-canary.pid"
        try:
            _run_lsan_canary_process(
                str(canary_helper), {
                    "LEO2_TEST_CANARY_MODE": "detached",
                    "LEO2_TEST_CANARY_PID": str(detached_pid),
                }, "modeled-detached", timeout_seconds=1.0)
        except CampaignError as error:
            if "residual descendants" not in str(error):
                raise
        else:
            raise CampaignError(
                "LSan canary accepted a detached residual descendant")
        require_process_gone(
            detached_pid, "detached LSan canary cleanup")

        timeout_pid = Path(directory) / "timeout-canary.pid"
        timeout_started = time.monotonic()
        try:
            _run_lsan_canary_process(
                str(canary_helper), {
                    "LEO2_TEST_CANARY_MODE": "sleep",
                    "LEO2_TEST_CANARY_PID": str(timeout_pid),
                }, "modeled-timeout", timeout_seconds=0.05)
        except CampaignError as error:
            if "timed out" not in str(error):
                raise
        else:
            raise CampaignError("LSan canary timeout was accepted")
        if time.monotonic() - timeout_started > 3.0:
            raise CampaignError(
                "LSan canary timeout recovery was not bounded")
        require_process_gone(timeout_pid, "timed-out LSan canary cleanup")

        if hasattr(signal, "setitimer") and hasattr(signal, "ITIMER_REAL"):
            interrupted_pid = Path(directory) / "interrupted-canary.pid"
            prior_handler = signal.getsignal(signal.SIGALRM)

            def interrupt_canary(_signum, _frame):
                raise KeyboardInterrupt()

            signal.signal(signal.SIGALRM, interrupt_canary)
            signal.setitimer(signal.ITIMER_REAL, 0.1)
            try:
                try:
                    _run_lsan_canary_process(
                        str(canary_helper), {
                            "LEO2_TEST_CANARY_MODE": "sleep",
                            "LEO2_TEST_CANARY_PID": str(interrupted_pid),
                        }, "modeled-interrupt", timeout_seconds=5.0)
                except KeyboardInterrupt:
                    pass
                else:
                    raise CampaignError(
                        "LSan canary swallowed KeyboardInterrupt")
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0.0)
                signal.signal(signal.SIGALRM, prior_handler)
            require_process_gone(
                interrupted_pid, "interrupted LSan canary cleanup")

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
                for job in leak["jobs"]) or
            any(job["capture_limits_bytes"] != _leak_capture_limits()
                for job in leak["jobs"]) or
            leak["execution_policy"] != _leak_execution_policy()):
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

    if (sys.platform.startswith("linux") and hasattr(os, "fork") and
            hasattr(os, "setsid") and Path("/proc/self/stat").is_file()):
        replay_helper = Path(directory) / "fake-adversarial-leak-replay.py"
        replay_helper.write_text(
            "#!{}\n"
            "import os, signal, sys, time\n"
            "seed = int(sys.argv[1], 0)\n"
            "iterations = int(sys.argv[2], 0)\n"
            "mode = os.environ['LEO2_TEST_LEAK_MODE']\n"
            "pid_path = os.environ['LEO2_TEST_LEAK_PID']\n"
            "capture_limit = int(os.environ['LEO2_TEST_CAPTURE_LIMIT'])\n"
            "def publish(pid):\n"
            "    with open(pid_path, 'w') as output:\n"
            "        output.write(str(pid))\n"
            "        output.flush()\n"
            "        os.fsync(output.fileno())\n"
            "def marker():\n"
            "    print('leopard2_fuzz_replay seed={{}} iterations={{}} "
            "passed'.format(seed, iterations), flush=True)\n"
            "def fill_capture(byte, count):\n"
            "    chunk = byte * 65536\n"
            "    written = 0\n"
            "    while written < count:\n"
            "        written += os.write(1, chunk[:min(len(chunk), "
            "count - written)])\n"
            "if mode == 'excessive-output':\n"
            "    publish(os.getpid())\n"
            "    signal.signal(signal.SIGXFSZ, signal.SIG_DFL)\n"
            "    fill_capture(b'X', capture_limit + 1)\n"
            "    raise SystemExit(99)\n"
            "if mode == 'timeout-output':\n"
            "    publish(os.getpid())\n"
            "    fill_capture(b'T', capture_limit - 1)\n"
            "    time.sleep(30)\n"
            "if mode == 'detached-output':\n"
            "    pid = os.fork()\n"
            "    if pid == 0:\n"
            "        os.setsid()\n"
            "        publish(os.getpid())\n"
            "        fill_capture(b'D', capture_limit - 1)\n"
            "        time.sleep(30)\n"
            "        os._exit(0)\n"
            "    while not os.path.exists(pid_path):\n"
            "        time.sleep(0.005)\n"
            "    marker()\n"
            "    raise SystemExit(0)\n"
            "if mode == 'detached':\n"
            "    pid = os.fork()\n"
            "    if pid == 0:\n"
            "        os.setsid()\n"
            "        publish(os.getpid())\n"
            "        time.sleep(30)\n"
            "        os._exit(0)\n"
            "    while not os.path.exists(pid_path):\n"
            "        time.sleep(0.005)\n"
            "    marker()\n"
            "    raise SystemExit(0)\n"
            "if mode == 'timeout-tree':\n"
            "    pid = os.fork()\n"
            "    if pid == 0:\n"
            "        os.setsid()\n"
            "        publish(os.getpid())\n"
            "        time.sleep(30)\n"
            "        os._exit(0)\n"
            "    while not os.path.exists(pid_path):\n"
            "        time.sleep(0.005)\n"
            "    time.sleep(30)\n"
            "publish(os.getpid())\n"
            "time.sleep(30)\n".format(sys.executable),
            encoding="utf-8")
        replay_helper.chmod(0o700)
        replay_identity = lab._file_identity(replay_helper)

        def adversarial_job(name, mode, pid_path, timeout_seconds):
            value = clone(leak["jobs"][0])
            value["id"] = "adversarial-{}".format(name)
            value["command"] = [
                str(replay_helper), str(value["seed"]),
                str(value["iterations"])]
            value["cwd"] = directory
            value["timeout_seconds"] = timeout_seconds
            value["environment"] = dict(value["environment"])
            value["environment"]["LEO2_TEST_LEAK_MODE"] = mode
            value["environment"]["LEO2_TEST_LEAK_PID"] = str(pid_path)
            value["environment"]["LEO2_TEST_CAPTURE_LIMIT"] = str(
                _leak_capture_limit(value, "stdout"))
            value["executable"] = replay_identity
            unsigned_job = dict(value)
            unsigned_job.pop("job_digest", None)
            value["job_digest"] = lab._digest(unsigned_job)
            return value

        original_execute_popen = subprocess.Popen
        original_execute_cleanup = globals()[
            "_cleanup_lsan_process_tree_with_retry"]

        def execution_resource_snapshot():
            return (
                len(os.listdir("/proc/self/fd")),
                {
                    (record["pid"], record["starttime_ticks"])
                    for record in _current_lsan_direct_children()
                })

        execute_precedence_cases = (
            (KeyboardInterrupt(), SystemExit(53), KeyboardInterrupt, None),
            (SystemExit(51), KeyboardInterrupt(), SystemExit, 51),
            (CampaignError("ordinary launch failure"), SystemExit(55),
             SystemExit, 55),
        )
        for index, (launch_error, cleanup_terminal, expected_type,
                    expected_code) in enumerate(execute_precedence_cases):
            precedence_job = adversarial_job(
                "cleanup-precedence-{}".format(index), "sleep",
                Path(directory) / "unused-cleanup-precedence.pid", 1.0)
            before_resources = execution_resource_snapshot()

            def fail_execute_launch(
                    *unused_args, error=launch_error, **unused_kwargs):
                del unused_args, unused_kwargs
                raise error

            def cleanup_then_terminal(
                    active, terminal=cleanup_terminal):
                original_execute_cleanup(active)
                raise terminal

            subprocess.Popen = fail_execute_launch
            globals()["_cleanup_lsan_process_tree_with_retry"] = (
                cleanup_then_terminal)
            try:
                try:
                    _execute_leak_job(
                        leak, precedence_job,
                        str(Path(directory) /
                            "cleanup-precedence-results-{}".format(index)))
                except BaseException as observed:
                    if not isinstance(observed, expected_type):
                        raise CampaignError(
                            "LSan execution cleanup precedence {} "
                            "propagated {} instead of {}".format(
                                index, type(observed).__name__,
                                expected_type.__name__))
                    if (expected_code is not None and
                            getattr(observed, "code", None) != expected_code):
                        raise CampaignError(
                            "LSan execution cleanup precedence {} changed "
                            "SystemExit code".format(index))
                else:
                    raise CampaignError(
                        "LSan execution cleanup precedence {} unexpectedly "
                        "completed".format(index))
            finally:
                subprocess.Popen = original_execute_popen
                globals()["_cleanup_lsan_process_tree_with_retry"] = (
                    original_execute_cleanup)
            if execution_resource_snapshot() != before_resources:
                raise CampaignError(
                    "LSan execution cleanup precedence {} leaked resources"
                    .format(index))

        excessive_replay_pid = Path(directory) / "excessive-replay.pid"
        excessive_job = adversarial_job(
            "excessive-output", "excessive-output",
            excessive_replay_pid, 1.0)
        excessive_result = _execute_leak_job(
            leak, excessive_job,
            str(Path(directory) / "excessive-replay-results"))
        if (excessive_result["outcome"] != "evidence_invalid" or
                excessive_result.get("exit_code") !=
                    -int(signal.SIGXFSZ) or
                "capture byte limit (SIGXFSZ)" not in
                    excessive_result.get("detail", "")):
            raise CampaignError(
                "LSan replay did not classify excessive output exactly: "
                "{}".format(excessive_result))
        excessive_stdout = (
            _leak_job_directory(
                Path(directory) / "excessive-replay-results",
                excessive_job["id"]) / "stdout.txt")
        if excessive_stdout.stat().st_size > _leak_capture_limit(
                excessive_job, "stdout"):
            raise CampaignError(
                "LSan replay retained stdout beyond its signed limit")
        require_process_gone(
            excessive_replay_pid, "excessive-output LSan replay cleanup")

        timeout_output_pid = Path(directory) / "timeout-output.pid"
        timeout_output_job = adversarial_job(
            "timeout-output", "timeout-output", timeout_output_pid, 0.05)
        timeout_output_result = _execute_leak_job(
            leak, timeout_output_job,
            str(Path(directory) / "timeout-output-results"))
        if (timeout_output_result["outcome"] != "timeout" or
                "signed timeout" not in
                    timeout_output_result.get("detail", "")):
            raise CampaignError(
                "LSan replay did not bound a timed-out capture writer")
        timeout_output_stdout = (
            _leak_job_directory(
                Path(directory) / "timeout-output-results",
                timeout_output_job["id"]) / "stdout.txt")
        if timeout_output_stdout.stat().st_size > _leak_capture_limit(
                timeout_output_job, "stdout"):
            raise CampaignError(
                "timed-out LSan replay exceeded its signed capture limit")
        require_process_gone(
            timeout_output_pid, "timed-out output-writer cleanup")

        detached_output_pid = Path(directory) / "detached-output.pid"
        detached_output_job = adversarial_job(
            "detached-output", "detached-output",
            detached_output_pid, 1.0)
        detached_output_result = _execute_leak_job(
            leak, detached_output_job,
            str(Path(directory) / "detached-output-results"))
        if (detached_output_result["outcome"] != "evidence_invalid" or
                "residual descendants" not in
                    detached_output_result.get("detail", "") or
                "capture byte limit" not in
                    detached_output_result.get("detail", "")):
            raise CampaignError(
                "LSan replay accepted a detached capture writer")
        detached_output_stdout = (
            _leak_job_directory(
                Path(directory) / "detached-output-results",
                detached_output_job["id"]) / "stdout.txt")
        if detached_output_stdout.stat().st_size > _leak_capture_limit(
                detached_output_job, "stdout"):
            raise CampaignError(
                "detached LSan writer exceeded its signed capture limit")
        require_process_gone(
            detached_output_pid, "detached output-writer cleanup")

        detached_replay_pid = Path(directory) / "detached-replay.pid"
        detached_job = adversarial_job(
            "detached", "detached", detached_replay_pid, 1.0)
        detached_started = time.monotonic()
        detached_result = _execute_leak_job(
            leak, detached_job,
            str(Path(directory) / "detached-replay-results"))
        if (detached_result["outcome"] != "evidence_invalid" or
                "residual descendants" not in
                detached_result.get("detail", "") or
                time.monotonic() - detached_started > 3.0):
            raise CampaignError(
                "LSan replay accepted or stalled on a detached child "
                "holding its capture descriptors")
        detached_stdout = (
            _leak_job_directory(
                Path(directory) / "detached-replay-results",
                detached_job["id"]) / "stdout.txt").read_bytes()
        if detached_stdout != _expected_leak_stdout(detached_job):
            raise CampaignError(
                "detached LSan replay changed its retained stdout capture")
        require_process_gone(
            detached_replay_pid, "detached LSan replay cleanup")

        timeout_replay_pid = Path(directory) / "timeout-replay.pid"
        timeout_job = adversarial_job(
            "timeout", "timeout-tree", timeout_replay_pid, 0.05)
        timeout_started = time.monotonic()
        timeout_result = _execute_leak_job(
            leak, timeout_job,
            str(Path(directory) / "timeout-replay-results"))
        if (timeout_result["outcome"] != "timeout" or
                "signed timeout" not in timeout_result.get("detail", "") or
                "residual descendants" not in
                timeout_result.get("detail", "") or
                time.monotonic() - timeout_started > 3.0):
            raise CampaignError(
                "LSan replay timeout teardown was not bounded and exact")
        require_process_gone(
            timeout_replay_pid, "timed-out LSan replay cleanup")

        if hasattr(signal, "setitimer") and hasattr(signal, "ITIMER_REAL"):
            interrupted_replay_pid = (
                Path(directory) / "interrupted-replay.pid")
            interrupted_job = adversarial_job(
                "interrupt", "sleep", interrupted_replay_pid, 5.0)
            prior_handler = signal.getsignal(signal.SIGALRM)

            def interrupt_replay(_signum, _frame):
                raise KeyboardInterrupt()

            signal.signal(signal.SIGALRM, interrupt_replay)
            signal.setitimer(signal.ITIMER_REAL, 0.1)
            try:
                try:
                    _execute_leak_job(
                        leak, interrupted_job,
                        str(Path(directory) /
                            "interrupted-replay-results"))
                except KeyboardInterrupt:
                    pass
                else:
                    raise CampaignError(
                        "LSan replay swallowed KeyboardInterrupt")
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0.0)
                signal.signal(signal.SIGALRM, prior_handler)
            require_process_gone(
                interrupted_replay_pid,
                "interrupted LSan replay cleanup")

        cleanup_interrupt_pid = Path(directory) / "cleanup-interrupt.pid"
        cleanup_interrupt_job = adversarial_job(
            "cleanup-interrupt", "sleep", cleanup_interrupt_pid, 0.05)
        original_cleanup_checkpoint = globals()["_lsan_cleanup_checkpoint"]
        cleanup_checkpoints = []

        def interrupt_first_cleanup(phase):
            cleanup_checkpoints.append(phase)
            if len(cleanup_checkpoints) == 1:
                raise KeyboardInterrupt()

        globals()["_lsan_cleanup_checkpoint"] = interrupt_first_cleanup
        try:
            try:
                _execute_leak_job(
                    leak, cleanup_interrupt_job,
                    str(Path(directory) /
                        "cleanup-interrupt-results"))
            except KeyboardInterrupt:
                pass
            else:
                raise CampaignError(
                    "LSan replay swallowed cleanup-time KeyboardInterrupt")
        finally:
            globals()["_lsan_cleanup_checkpoint"] = (
                original_cleanup_checkpoint)
        if len(cleanup_checkpoints) < 3:
            raise CampaignError(
                "LSan replay did not retry interrupted cleanup")
        require_process_gone(
            cleanup_interrupt_pid,
            "cleanup-interrupted LSan replay cleanup")

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

    wrong_capture_limit = clone(leak)
    wrong_capture_job = wrong_capture_limit["jobs"][0]
    wrong_capture_job["capture_limits_bytes"]["stdout"] += 1
    expect_invalid_manifest(
        resign_manifest(wrong_capture_limit, wrong_capture_job),
        "job capture byte limit")

    wrong_capture_policy = clone(leak)
    wrong_capture_policy["execution_policy"]["capture_limits_bytes"][
        "stderr"] += 1
    expect_invalid_manifest(
        resign_manifest(wrong_capture_policy),
        "manifest capture byte limit")

    historical_leak_schema = clone(leak)
    historical_leak_schema["schema"] = "leopard2-fuzz-leak-campaign/v3"
    expect_invalid_manifest(
        resign_manifest(historical_leak_schema),
        "historical leak schema under current semantics")

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

    def resume_resource_snapshot():
        descriptor_count = len(os.listdir("/proc/self/fd"))
        children = {
            (record["pid"], record["starttime_ticks"])
            for record in _current_lsan_direct_children()
        }
        return descriptor_count, children

    def expect_resume_identity_race(
            label, phase, mutate, restore):
        original_checkpoint = globals()["_leak_resume_checkpoint"]
        mutations = []
        before_resources = resume_resource_snapshot()

        def race_checkpoint(observed_phase):
            if observed_phase == phase and not mutations:
                mutate()
                mutations.append(observed_phase)

        globals()["_leak_resume_checkpoint"] = race_checkpoint
        observed_error = None
        try:
            try:
                _read_resumable_leak_result(output_dir, first_job)
            except (CampaignError, lab.LabError) as error:
                observed_error = error
            else:
                raise CampaignError(
                    "LSan resume accepted {}".format(label))
        finally:
            globals()["_leak_resume_checkpoint"] = original_checkpoint
            restore()
        if not mutations:
            raise CampaignError(
                "LSan resume did not exercise {}".format(label))
        if not isinstance(observed_error, _LeakResumeIdentityError):
            raise CampaignError(
                "LSan resume classified {} as {} instead of an identity "
                "failure".format(label, type(observed_error).__name__))
        if resume_resource_snapshot() != before_resources:
            raise CampaignError(
                "LSan resume leaked resources after {}".format(label))

    capture_path = first_dir / "stdout.txt"
    capture_saved = first_dir / ".stdout.saved"

    def replace_capture():
        capture_path.rename(capture_saved)
        capture_path.write_bytes(b"replacement capture\n")

    def restore_capture():
        try:
            capture_path.unlink()
        except FileNotFoundError:
            pass
        if capture_saved.exists():
            capture_saved.rename(capture_path)

    expect_resume_identity_race(
        "a retained capture rename/replacement",
        "after-result-read", replace_capture, restore_capture)

    result_path = first_dir / "result.json"
    result_saved = first_dir / ".result.saved"

    def replace_result():
        result_path.rename(result_saved)
        result_path.write_text("{}\n", encoding="ascii")

    def restore_result():
        try:
            result_path.unlink()
        except FileNotFoundError:
            pass
        if result_saved.exists():
            result_saved.rename(result_path)

    expect_resume_identity_race(
        "a retained result rename/replacement",
        "after-retain", replace_result, restore_result)

    expect_resume_identity_race(
        "a late retained capture rename/replacement",
        "before-close", replace_capture, restore_capture)

    def aba_capture():
        capture_path.rename(capture_saved)
        capture_path.write_bytes(b"temporary capture replacement\n")
        capture_path.unlink()
        capture_saved.rename(capture_path)

    expect_resume_identity_race(
        "a retained capture ABA replacement",
        "after-result-read", aba_capture, lambda: None)

    def aba_result():
        result_path.rename(result_saved)
        result_path.write_text("{}\n", encoding="ascii")
        result_path.unlink()
        result_saved.rename(result_path)

    expect_resume_identity_race(
        "a retained result ABA replacement",
        "after-result-read", aba_result, lambda: None)

    first_saved_dir = Path(str(first_dir) + ".resume-saved")

    def replace_job_directory():
        first_dir.rename(first_saved_dir)
        first_dir.mkdir()

    def restore_job_directory():
        if first_dir.exists():
            first_dir.rmdir()
        if first_saved_dir.exists():
            first_saved_dir.rename(first_dir)

    expect_resume_identity_race(
        "a retained job-directory descriptor/path disagreement",
        "after-retain", replace_job_directory, restore_job_directory)

    original_resume_checkpoint = globals()["_leak_resume_checkpoint"]
    for label, phase, terminal, expected_type, expected_code in (
            ("KeyboardInterrupt", "after-result-read", KeyboardInterrupt(),
             KeyboardInterrupt, None),
            ("SystemExit", "before-close", SystemExit(43),
             SystemExit, 43)):
        before_resources = resume_resource_snapshot()
        interruptions = []

        def interrupt_resume(
                observed_phase, expected_phase=phase,
                exception=terminal):
            if observed_phase == expected_phase and not interruptions:
                interruptions.append(observed_phase)
                raise exception

        globals()["_leak_resume_checkpoint"] = interrupt_resume
        try:
            try:
                _read_resumable_leak_result(output_dir, first_job)
            except BaseException as observed:
                if not isinstance(observed, expected_type):
                    raise CampaignError(
                        "LSan resume propagated {} instead of {}".format(
                            type(observed).__name__, label))
                if (expected_code is not None and
                        getattr(observed, "code", None) != expected_code):
                    raise CampaignError(
                        "LSan resume changed {} code".format(label))
            else:
                raise CampaignError(
                    "LSan resume swallowed {}".format(label))
        finally:
            globals()["_leak_resume_checkpoint"] = (
                original_resume_checkpoint)
        if interruptions != [phase]:
            raise CampaignError(
                "LSan resume did not inject {} exactly once".format(label))
        if resume_resource_snapshot() != before_resources:
            raise CampaignError(
                "LSan resume leaked resources after {}".format(label))

    if _read_resumable_leak_result(output_dir, first_job) is None:
        raise CampaignError(
            "descriptor-bound LSan resume lost valid evidence")

    resign_result(
        first_dir / "result.json",
        lambda value: value.update({
            "schema": "leopard2-fuzz-leak-result/v1"}))
    if _read_resumable_leak_result(output_dir, first_job) is not None:
        raise CampaignError(
            "current leak runner resumed a historical result schema")
    expect_one_repair("historical result-schema evidence")

    (first_dir / "stderr.txt").write_text(
        "ERROR: LeakSanitizer: modeled diagnostic\n", encoding="utf-8")
    resign_result(
        first_dir / "result.json",
        lambda value: value["outputs"].update({
            "stderr": lab._content_identity(first_dir / "stderr.txt")}))
    expect_audit_failure("a sanitizer diagnostic")
    expect_one_repair("diagnostic evidence")

    with (first_dir / "stdout.txt").open("r+b") as oversized:
        oversized.truncate(
            _leak_capture_limit(first_job, "stdout") + 1)
    resign_result(
        first_dir / "result.json",
        lambda value: value["outputs"].update({
            "stdout": lab._content_identity(first_dir / "stdout.txt")}))
    expect_audit_failure("an oversized retained capture")
    if _read_resumable_leak_result(output_dir, first_job) is not None:
        raise CampaignError(
            "LSan runner resumed an oversized retained capture")
    expect_one_repair("oversized capture evidence")

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
        _read_resumable_leak_result(output_dir, first_job)
    except _LeakResumeIdentityError:
        pass
    else:
        raise CampaignError(
            "descriptor-bound LSan resume did not fail closed on a symlink")
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
        _read_resumable_leak_result(output_dir, first_job)
    except _LeakResumeIdentityError:
        pass
    else:
        raise CampaignError(
            "descriptor-bound LSan resume did not fail closed on a hard link")
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
        _read_resumable_leak_result(output_dir, leak["jobs"][1])
    except _LeakResumeIdentityError:
        pass
    else:
        raise CampaignError(
            "descriptor-bound LSan resume accepted an aliased job directory")
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
        _probe_runner_self_test(directory, nm_path, nm_identity)
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

        wrong_probe_policy = clone(first)
        wrong_probe_policy["source_spec"]["metadata"][
            "probe_execution_policy"]["capture_profiles"][
                "sanitizer_attestation"]["stdout_limit_bytes"] += 1
        expect_invalid(
            resign(wrong_probe_policy), "probe execution policy")

        missing_probe_policy = clone(first)
        del missing_probe_policy["source_spec"]["metadata"][
            "probe_execution_policy"]
        expect_invalid(
            resign(missing_probe_policy), "missing probe execution policy")

        historical = clone(first)
        historical["source_spec"]["schema"] = HISTORICAL_CAMPAIGN_SCHEMA
        del historical["source_spec"]["metadata"]["probe_execution_policy"]
        resign(historical)
        expect_invalid(historical, "offline historical campaign")
        _validate_campaign_manifest(
            historical, allow_historical_live_replay=True)

        historical_with_extension = clone(historical)
        historical_with_extension["source_spec"]["metadata"][
            "probe_execution_policy"] = _probe_execution_policy()
        expect_invalid(
            resign(historical_with_extension),
            "extended historical campaign")

        audit = {
            "schema": AUDIT_SCHEMA,
            "source_campaign_schema": CAMPAIGN_SCHEMA,
            "probe_policy_binding":
                _campaign_probe_policy_binding(CAMPAIGN_SCHEMA),
            "probe_execution_policy": _probe_execution_policy(),
            "manifest_digest": first["manifest_digest"],
            "job_count": len(first["jobs"]),
            "distinct_seed_count": len(first["jobs"]),
            "sanitizer_scope": dict(SANITIZER_SCOPE),
            "target_instrumentation": clone(first["source_spec"]["metadata"][
                "target_instrumentation"]),
            "summary": {
                "missing": 0, "success": len(first["jobs"])},
            "merged_results": {
                "summary": {
                    "missing": 0, "success": len(first["jobs"])},
                "results": [
                    {"job_id": job["id"]} for job in first["jobs"]],
            },
        }
        audit["audit_digest"] = lab._digest(audit)
        _validate_campaign_audit(audit, first)

        forged_audit_policy = clone(audit)
        forged_audit_policy["probe_execution_policy"][
            "capture_enforcement"]["soft_and_hard_limits"] = False
        forged_audit_policy["audit_digest"] = lab._digest({
            key: value for key, value in forged_audit_policy.items()
            if key != "audit_digest"
        })
        try:
            _validate_campaign_audit(forged_audit_policy, first)
        except CampaignError:
            pass
        else:
            raise CampaignError(
                "campaign audit accepted a forged probe policy")

        historical_audit = clone(audit)
        historical_audit["source_campaign_schema"] = \
            HISTORICAL_CAMPAIGN_SCHEMA
        historical_audit["probe_policy_binding"] = \
            _campaign_probe_policy_binding(HISTORICAL_CAMPAIGN_SCHEMA)
        historical_audit["manifest_digest"] = historical["manifest_digest"]
        historical_audit["audit_digest"] = lab._digest({
            key: value for key, value in historical_audit.items()
            if key != "audit_digest"
        })
        _validate_campaign_audit(historical_audit, historical)

        old_audit = clone(audit)
        old_audit["schema"] = "leopard2-fuzz-campaign-audit/v4"
        old_audit["audit_digest"] = lab._digest({
            key: value for key, value in old_audit.items()
            if key != "audit_digest"
        })
        try:
            _validate_campaign_audit(old_audit)
        except CampaignError:
            pass
        else:
            raise CampaignError("current validator accepted a v4 audit")

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
                # A successful job must remain live long enough for one
                # complete system-wide /proc affinity/RSS sample even on a
                # loaded CI host.  This does not relax the generic fork gate:
                # the leak-enabled branch below still creates a second
                # process and must remain evidence-invalid.
                "time.sleep(1.5)\n"
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
