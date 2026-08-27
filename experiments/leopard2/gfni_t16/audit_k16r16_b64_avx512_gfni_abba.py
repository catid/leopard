#!/usr/bin/env python3

"""Independent, no-timing replay of the frozen K16/R16/B64 GFNI campaign.

This file deliberately does not import the acquisition runner.  It accepts a
completed report and its append-only journal, reopens every retained raw file,
executes the separately frozen benchmark JSON validator directly from source,
reconstructs the complete event stream, and recomputes the paired-log
inference and promotion gates.  It never executes the benchmark binary.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import stat
import statistics
import subprocess
import sys
import tarfile
import tempfile
import types
from pathlib import Path
from typing import Any


TARGET_ROUNDS = 25
INACTIVE_ROUNDS = 2
ITERATIONS = 31
WARMUP = 64
REUSE = 32768
MIN_RETAINED_TIMER_WINDOW_US = 250.0
MAX_ROUND_ATTEMPTS = 5
EXPECTED_CPU = 52
EXPECTED_SIBLING = 116
T_CRITICAL_24 = 2.0638985616280205
MAX_REPORT_BYTES = 64 * 1024 * 1024
MAX_JOURNAL_BYTES = 128 * 1024 * 1024
MAX_RESULT_BYTES = 8 * 1024 * 1024
MAX_COMMAND_BYTES = 64 * 1024
MAX_STREAM_BYTES = 1024 * 1024
MAX_SOURCE_ARCHIVE_BYTES = 512 * 1024 * 1024
MAX_VALIDATOR_BYTES = 8 * 1024 * 1024
REPORT_SCHEMA = "leopard2-k16r16-b64-avx512-gfni-frozen-abba/v1"
AUDIT_SCHEMA = "leopard2-k16r16-b64-avx512-gfni-abba-audit/v1"
BENCHMARK_SCHEMA = "leopard2-benchmark-v34"
DIAGNOSTIC_OPTION = "--k16r16-b64-avx512-gfni-mode"
CONTRACT = (
    "LEGACY_HIGH_V1,GF8,AUTO,K=16,R=16,T=16,B=64,native_layout,"
    "balanced_packed_terminal,runtime_AVX512F_BW_VL_GFNI,startup_KAT,"
    "calibrated_AMD_1A_08,one_shot_and_one_item_batch"
)
CHILD_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}
GIT_ENVIRONMENT = dict(CHILD_ENVIRONMENT, **{
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_CONFIG_SYSTEM": "/dev/null",
    "GIT_NO_REPLACE_OBJECTS": "1",
    "GIT_OPTIONAL_LOCKS": "0",
})

# id, K, R, bytes, rounds, role, profile, field, backend, rationale
CELLS = (
    ("target-k16-r16-b64-high-gf8-auto", 16, 16, 64, TARGET_ROUNDS,
     "target", "high", "gf8", "auto", "exact production selector cell"),
    ("inactive-b63", 16, 16, 63, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "lower byte boundary"),
    ("inactive-b65", 16, 16, 65, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "upper byte boundary"),
    ("inactive-k15", 15, 16, 64, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "lower K boundary"),
    ("inactive-k17", 17, 16, 64, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "upper K boundary"),
    ("inactive-r15", 16, 15, 64, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "lower R boundary"),
    ("inactive-r17", 16, 17, 64, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "upper R boundary"),
    ("inactive-low-gf8", 16, 16, 64, INACTIVE_ROUNDS,
     "inactive", "low", "gf8", "auto", "profile boundary"),
    ("inactive-high-gf16", 16, 16, 64, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "auto", "field boundary"),
    ("inactive-explicit-avx2", 16, 16, 64, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "avx2", "backend boundary"),
    ("inactive-b128", 16, 16, 128, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "extended byte control"),
    ("inactive-k8", 8, 16, 64, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "extended lower K control"),
    ("inactive-k32", 32, 16, 64, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "extended upper K control"),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def exact_keys(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    require(type(value) is dict and set(value) == keys,
            f"{label} keys changed")
    return value


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def reject_constant(value: str) -> Any:
    raise RuntimeError(f"non-finite JSON constant: {value}")


def finite_float_token(value: str) -> float:
    result = float(value)
    require(math.isfinite(result), f"non-finite JSON float: {value}")
    return result


def strict_json_bytes(data: bytes, label: str) -> Any:
    require(len(data) > 0, f"{label} is empty")
    text = data.decode("utf-8", errors="strict")
    result = json.loads(
        text, object_pairs_hook=strict_object,
        parse_constant=reject_constant, parse_float=finite_float_token)
    require_finite_tree(result, label)
    return result


def require_finite_tree(value: Any, label: str) -> None:
    if type(value) is float:
        require(math.isfinite(value), f"{label} contains a non-finite number")
    elif type(value) is list:
        for item in value:
            require_finite_tree(item, label)
    elif type(value) is dict:
        for item in value.values():
            require_finite_tree(item, label)
    else:
        require(type(value) in (type(None), bool, int, str),
                f"{label} contains an unsupported JSON value")


def metadata_tuple(value: os.stat_result) -> tuple[int, ...]:
    return (value.st_mode, value.st_nlink, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns, value.st_dev, value.st_ino)


def canonical_path(path: Path, *, must_exist: bool) -> Path:
    require(path.is_absolute(), f"path is not absolute: {path}")
    parent = path.parent
    resolved_parent = parent.resolve(strict=True)
    require(resolved_parent == parent and
            stat.S_ISDIR(parent.lstat().st_mode) and not parent.is_symlink(),
            f"path parent is aliased or symbolic: {path}")
    if must_exist:
        require(path.resolve(strict=True) == path,
                f"path is aliased or symbolic: {path}")
    return path


def stable_file(path: Path, *, maximum_bytes: int,
                return_bytes: bool = False,
                single_link: bool = True,
                read_only: bool = False) -> Any:
    canonical_path(path, must_exist=True)
    before = path.lstat()
    require(stat.S_ISREG(before.st_mode) and not path.is_symlink(),
            f"not one regular file: {path}")
    require(before.st_size <= maximum_bytes,
            f"file exceeds its byte bound: {path}")
    if single_link:
        require(before.st_nlink == 1, f"file is multiply linked: {path}")
    if read_only:
        require(before.st_mode & 0o222 == 0,
                f"frozen artifact remains writable: {path}")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    descriptor = os.open(path, flags)
    try:
        opened = os.fstat(descriptor)
        require(metadata_tuple(opened) == metadata_tuple(before),
                f"file changed before read: {path}")
        digest = hashlib.sha256()
        retained = bytearray() if return_bytes else None
        total = 0
        while True:
            block = os.read(descriptor, min(
                1024 * 1024, maximum_bytes + 1 - total))
            if not block:
                break
            total += len(block)
            require(total <= maximum_bytes,
                    f"file exceeds its byte bound: {path}")
            digest.update(block)
            if retained is not None:
                retained.extend(block)
        closed = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after = path.lstat()
    require(metadata_tuple(before) == metadata_tuple(closed) ==
            metadata_tuple(after) and total == before.st_size,
            f"file changed while reading: {path}")
    identity = {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size": before.st_size,
        "mode": stat.S_IMODE(before.st_mode),
        "links": before.st_nlink,
        "device": before.st_dev,
        "inode": before.st_ino,
        "mtime_ns": before.st_mtime_ns,
        "ctime_ns": before.st_ctime_ns,
    }
    return (bytes(retained), identity) if return_bytes else identity


IDENTITY_KEYS = {
    "path", "sha256", "size", "mode", "links", "device", "inode",
    "mtime_ns", "ctime_ns",
}


def validate_identity(identity: Any, label: str) -> dict[str, Any]:
    identity = exact_keys(identity, IDENTITY_KEYS, label)
    require(type(identity["path"]) is str and
            type(identity["sha256"]) is str and
            len(identity["sha256"]) == 64 and
            all(value in "0123456789abcdef"
                for value in identity["sha256"]),
            f"{label} has malformed path or SHA-256")
    for key in IDENTITY_KEYS - {"path", "sha256"}:
        require(type(identity[key]) is int and identity[key] >= 0,
                f"{label} has malformed {key}")
    return identity


def identity_matches_acquisition_or_seal(
        expected: dict[str, Any], observed: dict[str, Any]) -> bool:
    invariant = IDENTITY_KEYS - {"mode", "ctime_ns"}
    if any(observed[key] != expected[key] for key in invariant):
        return False
    if observed["mode"] == expected["mode"]:
        return observed["ctime_ns"] == expected["ctime_ns"]
    sealed_mode = expected["mode"] & ~0o222
    return (expected["mode"] & 0o222 != 0 and
            observed["mode"] == sealed_mode and
            observed["ctime_ns"] >= expected["ctime_ns"])


def reopen_identity(identity: Any, label: str, maximum_bytes: int,
                    *, read_only: bool = False,
                    return_bytes: bool = False,
                    allow_monotonic_seal: bool = False) -> Any:
    expected = validate_identity(identity, label)
    observed = stable_file(
        Path(expected["path"]), maximum_bytes=maximum_bytes,
        return_bytes=return_bytes, single_link=True, read_only=read_only)
    if return_bytes:
        data, observed_identity = observed
        require(observed_identity == expected or
                (allow_monotonic_seal and
                 identity_matches_acquisition_or_seal(
                     expected, observed_identity)),
                f"{label} identity changed outside the sealing policy")
        return data, observed_identity
    require(observed == expected or
            (allow_monotonic_seal and
             identity_matches_acquisition_or_seal(expected, observed)),
            f"{label} identity changed outside the sealing policy")
    return observed


def strict_json_file(path: Path, maximum_bytes: int, label: str) \
        -> tuple[Any, dict[str, Any], bytes]:
    data, identity = stable_file(
        path, maximum_bytes=maximum_bytes, return_bytes=True,
        single_link=True)
    return strict_json_bytes(data, label), identity, data


def parse_journal(path: Path) -> tuple[list[dict[str, Any]],
                                      dict[str, Any], bytes]:
    data, identity = stable_file(
        path, maximum_bytes=MAX_JOURNAL_BYTES, return_bytes=True,
        single_link=True)
    require(data.endswith(b"\n"), "journal lacks its terminal newline")
    lines = data.splitlines(keepends=True)
    events: list[dict[str, Any]] = []
    for index, line in enumerate(lines):
        require(line.endswith(b"\n") and line != b"\n",
                f"journal line {index} is empty or unterminated")
        value = strict_json_bytes(line, f"journal line {index}")
        require(type(value) is dict, f"journal line {index} is not an object")
        canonical = (json.dumps(
            value, sort_keys=True, separators=(",", ":")) + "\n").encode()
        require(line == canonical, f"journal line {index} is not canonical")
        events.append(value)
    require(events, "journal has no events")
    return events, identity, data


def finite_positive(value: Any, label: str) -> float:
    require(type(value) in (int, float), f"{label} is not numeric")
    result = float(value)
    require(math.isfinite(result) and result > 0,
            f"{label} is not finite and positive")
    return result


def timing(document: dict[str, Any], name: str) \
        -> tuple[float, list[float], float]:
    summary = document["metrics"][name]
    samples_raw = summary["samples_us_per_batch_call"]
    require(type(samples_raw) is list and len(samples_raw) == ITERATIONS,
            f"{name} sample count changed")
    samples = [finite_positive(item, f"{name} sample")
               for item in samples_raw]
    observed = statistics.median(samples)
    reported = finite_positive(
        summary["median_us_per_batch_call"], f"{name} median")
    require(abs(observed - reported) <= 0.000003,
            f"{name} median differs from retained samples")
    minimum = min(item * REUSE for item in samples)
    require(minimum >= MIN_RETAINED_TIMER_WINDOW_US,
            f"{name} retained timer window is below the frozen floor")
    return reported, samples, minimum


def load_validator(path: Path, expected: dict[str, Any]) -> Any:
    source, _ = reopen_identity(
        expected, "frozen validator", MAX_VALIDATOR_BYTES,
        read_only=True, return_bytes=True)
    require(path == Path(expected["path"]), "validator path changed")
    module_name = "leopard2_independent_frozen_benchmark_validator"
    module = types.ModuleType(module_name)
    module.__file__ = str(path)
    sys.modules[module_name] = module
    exec(compile(source, str(path), "exec", dont_inherit=True, optimize=0),
         module.__dict__)
    require(callable(getattr(module, "validate_common", None)) and
            callable(getattr(module, "validate_workload_digests", None)),
            "frozen validator lacks required functions")
    return module


def summarize(logs: list[float], inferential: bool) -> dict[str, Any]:
    require(logs and all(math.isfinite(item) for item in logs),
            "summary logs are empty or non-finite")
    mean = statistics.mean(logs)
    result: dict[str, Any] = {
        "rounds": len(logs),
        "log_mean": mean,
        "geometric_mean_speedup": math.exp(mean),
        "inferential": inferential,
    }
    if inferential:
        require(len(logs) == TARGET_ROUNDS,
                "target inferential round count changed")
        sample_sd = statistics.stdev(logs)
        radius = T_CRITICAL_24 * sample_sd / math.sqrt(len(logs))
        result.update({
            "log_sample_sd": sample_sd,
            "t_critical": T_CRITICAL_24,
            "ci95": [math.exp(mean - radius), math.exp(mean + radius)],
        })
    else:
        result.update({
            "log_sample_sd": statistics.stdev(logs)
            if len(logs) > 1 else None,
            "t_critical": None,
            "ci95": None,
        })
    return result


def isolation_passed(isolation: Any, label: str) -> bool:
    isolation = exact_keys(
        isolation, {"benchmark_cpu", "reserved_sibling"}, label)
    for name, expected_cpu in (("benchmark_cpu", EXPECTED_CPU),
                               ("reserved_sibling", EXPECTED_SIBLING)):
        value = exact_keys(isolation[name],
                           {"cpu", "values", "idle", "nonidle", "total"},
                           f"{label} {name}")
        require(type(value["cpu"]) is int and value["cpu"] == expected_cpu and
                type(value["values"]) is list and len(value["values"]) >= 8 and
                all(type(item) is int and item >= 0
                    for item in value["values"]),
                f"{label} {name} CPU vector changed")
        values = value["values"]
        require(value["idle"] == values[3] + values[4] and
                value["nonidle"] == sum(values[:3]) + sum(values[5:8]) and
                value["total"] == sum(values),
                f"{label} {name} derived counters changed")
    return (isolation["benchmark_cpu"]["nonidle"] > 0 and
            isolation["reserved_sibling"]["total"] > 0 and
            isolation["reserved_sibling"]["idle"] > 0 and
            isolation["reserved_sibling"]["nonidle"] == 0)


def parse_controller_command(report: dict[str, Any], report_path: Path,
                             journal_path: Path) -> dict[str, str]:
    command = report["controller_command"]
    require(type(command) is list and len(command) >= 3 and
            all(type(item) is str for item in command),
            "controller command is malformed")
    allowed = {
        "--binary", "--build-binary", "--binary-sha256", "--runner-sha256",
        "--validator", "--validator-sha256", "--source-archive",
        "--source-archive-sha256", "--source-commit", "--source-tree",
        "--repository", "--cpu", "--sibling", "--output", "--journal",
        "--invocations", "--lock-fd",
    }
    required = allowed - {"--lock-fd"}
    require((len(command) - 1) % 2 == 0,
            "controller command arguments are not flag/value pairs")
    parsed: dict[str, str] = {}
    for index in range(1, len(command), 2):
        flag, value = command[index:index + 2]
        require(flag in allowed and flag not in parsed,
                f"unknown or duplicate controller option: {flag}")
        parsed[flag] = value
    require(required <= set(parsed) and set(parsed) <= allowed,
            "controller command option census changed")
    require(Path(command[0]) == Path(report["provenance"]["runner_pre"]["path"]),
            "controller executable is not the frozen runner")
    require(Path(parsed["--output"]) == report_path and
            Path(parsed["--journal"]) == journal_path,
            "controller report/journal paths changed")
    return parsed


def validate_canonical_lock(value: Any,
                            controller: dict[str, str]) -> dict[str, Any]:
    lock = exact_keys(
        value,
        {"path", "mode", "scope", "descriptor", "device", "inode"},
        "canonical lock")
    require(lock["path"] == "/tmp/leopard-gf8-authoritative.lock" and
            lock["mode"] in {
                "runner-acquired", "inherited-across-build-copy-campaign"} and
            lock["scope"] in {
                "campaign-only", "wrapper-build-copy-campaign"} and
            all(type(lock[key]) is int and lock[key] >= 0
                for key in ("descriptor", "device", "inode")),
            "canonical lock contract changed")
    if "--lock-fd" in controller:
        require(controller["--lock-fd"] == str(lock["descriptor"]) and
                lock["mode"] == "inherited-across-build-copy-campaign" and
                lock["scope"] == "wrapper-build-copy-campaign",
                "inherited lock descriptor or scope changed")
    else:
        require(lock["mode"] == "runner-acquired" and
                lock["scope"] == "campaign-only",
                "runner-acquired lock metadata changed")
    return lock


def expected_benchmark_command(binary: Path, output: Path,
                               cell: tuple[Any, ...], mode: int,
                               seed: int) -> list[str]:
    _, k, r, shard_bytes, _, _, profile, field, backend, _ = cell
    return [
        "/usr/bin/taskset", "-c", str(EXPECTED_CPU), str(binary),
        "--k", str(k), "--r", str(r), "--profile", profile,
        "--field", field, "--backend", backend,
        "--bytes", str(shard_bytes), "--loss", "1", "--batch", "1",
        "--reuse", str(REUSE), "--iterations", str(ITERATIONS),
        "--warmup", str(WARMUP), "--threads", "1", "--seed", str(seed),
        "--skip-legacy", "--retain-samples", "--measure-one-shot-encode",
        "--attest-source", DIAGNOSTIC_OPTION, str(mode),
        "--json", str(output),
    ]


def validate_launch(launch: Any, *, cell: tuple[Any, ...], mode: int,
                    label: str, seed: int, artifact_root: Path,
                    invocations: Path, binary: Path, source_commit: str,
                    source_tree: str, validator: Any,
                    raw_manifest: list[dict[str, Any]]) -> dict[str, Any]:
    launch_keys = {
        "label", "mode", "elapsed_ns", "reuse", "encode_us",
        "one_shot_us", "encode_samples_us", "one_shot_samples_us",
        "retained_timer_window_us", "workload_digests",
        "kernel_available", "kernel_qualified", "selector_selected",
        "observed_call_count", "source", "raw",
    }
    launch = exact_keys(launch, launch_keys, f"launch {label}")
    require(launch["label"] == label and launch["mode"] == mode and
            type(launch["elapsed_ns"]) is int and launch["elapsed_ns"] > 0 and
            launch["reuse"] == REUSE,
            f"launch identity changed: {label}")
    raw = exact_keys(launch["raw"],
                     {"command", "stdout", "stderr", "result"},
                     f"launch raw evidence {label}")
    invocation = invocations / label
    require(invocation.is_absolute() and invocation.parent == invocations and
            invocation.resolve(strict=True) == invocation and
            stat.S_ISDIR(invocation.lstat().st_mode) and
            not invocation.is_symlink(),
            f"invocation directory escaped or is aliased: {label}")
    require(stat.S_IMODE(invocation.lstat().st_mode) in {0o700, 0o500},
            f"invocation directory permissions changed: {label}")
    expected_paths = {
        "command": invocation / "command.json",
        "stdout": invocation / "stdout",
        "stderr": invocation / "stderr",
        "result": invocation / "result.json",
    }
    opened: dict[str, tuple[bytes, dict[str, Any]]] = {}
    for kind, maximum in (("command", MAX_COMMAND_BYTES),
                          ("stdout", MAX_STREAM_BYTES),
                          ("stderr", MAX_STREAM_BYTES),
                          ("result", MAX_RESULT_BYTES)):
        identity = validate_identity(raw[kind], f"{label} {kind}")
        require(identity["mode"] & 0o022 == 0 and
                identity["mode"] & 0o111 == 0,
                f"{label} {kind} permissions are unsafe")
        path = Path(identity["path"])
        require(path == expected_paths[kind] and
                path.is_relative_to(artifact_root),
                f"{label} {kind} path escaped or changed")
        opened[kind] = reopen_identity(
            identity, f"{label} {kind}", maximum, return_bytes=True,
            allow_monotonic_seal=True)
        raw_manifest.append({
            "label": label, "kind": kind, "path": str(path),
            "sha256": identity["sha256"], "size": identity["size"],
        })
    command_bytes = opened["command"][0]
    command = strict_json_bytes(command_bytes, f"{label} command")
    require(type(command) is list and all(type(item) is str for item in command),
            f"{label} command is not a string list")
    expected = expected_benchmark_command(
        binary, expected_paths["result"], cell, mode, seed)
    require(command == expected and command_bytes ==
            (json.dumps(expected, separators=(",", ":")) + "\n").encode(),
            f"{label} retained command changed")
    require(opened["stdout"][0] == b"" and opened["stderr"][0] == b"",
            f"{label} retained terminal streams are nonempty")
    document = strict_json_bytes(opened["result"][0], f"{label} result")
    require(type(document) is dict, f"{label} result is not an object")
    validator.validate_common(document, True)
    validator.validate_workload_digests(document)
    cell_id, k, r, shard_bytes, _, role, profile, field, backend, _ = cell
    parameters = document["parameters"]
    expected_profile = "legacy_high_v1" if profile == "high" else "low_v1"
    require(document["schema"] == BENCHMARK_SCHEMA and
            (parameters["K"], parameters["R"],
             parameters["shard_bytes"]) == (k, r, shard_bytes) and
            parameters["k16r16_b64_avx512_gfni_mode"] == mode and
            parameters["requested_profile"] == expected_profile and
            parameters["requested_field"] == field and
            parameters["requested_backend"] == backend and
            parameters["skip_legacy"] is True and
            parameters["retain_samples"] is True and
            parameters["measure_one_shot_encode"] is True and
            parameters["attest_source"] is True and
            parameters["batch"] == 1 and parameters["thread_count"] == 1 and
            parameters["reuse"] == REUSE and
            parameters["iterations"] == ITERATIONS and
            parameters["warmup"] == WARMUP and parameters["loss_count"] == 1 and
            parameters["seed"] == seed,
            f"{label} benchmark contract changed")
    require(document["resolved"]["profile"] == expected_profile and
            document["resolved"]["field"] == field and
            document["resolved"]["backend"] == "avx2",
            f"{label} resolved identity changed")
    build = document["build"]
    available = backend == "auto"
    selected = role == "target" and mode == 1
    require(build["source_commit"] == source_commit and
            build["source_tree"] == source_tree and
            build["source_tracked_dirty"] is False and
            build["k16r16_b64_avx512_gfni_diagnostic_mode"] == mode and
            build["k16r16_b64_avx512_gfni_diagnostic_disabled"] is (mode == 0) and
            build["k16r16_b64_avx512_gfni_mode_latched"] == mode and
            build["k16r16_b64_avx512_gfni_kernel_available"] is available and
            build["k16r16_b64_avx512_gfni_kernel_qualified"] is available and
            build["k16r16_b64_avx512_gfni_selector_expected_selected"] is selected and
            build["k16r16_b64_avx512_gfni_selector_selected"] is selected and
            build["k16r16_b64_avx512_gfni_observed_call_count"] ==
                (2 if selected else 0) and
            build["k16r16_b64_avx512_gfni_selector_contract"] == CONTRACT and
            build["k16r16_b64_avx512_gfni_timed_ordinary_encode_api"] ==
                "leo2_encode_batch:item_count=1:no_preflight_scratch" and
            build["k16r16_b64_avx512_gfni_timed_one_shot_encode_api"] ==
                "leo2_encode" and
            document["correctness"]["leopard2_round_trip"] is True,
            f"{label} selector/correctness contract changed")
    encode_us, encode_samples, encode_minimum = timing(
        document, "encode_execution")
    one_shot_us, one_shot_samples, one_shot_minimum = timing(
        document, "one_shot_encode")
    digests = document["workload_digests"]
    require(digests["algorithm"] == "fnv1a64",
            f"{label} digest algorithm changed")
    derived = dict(launch)
    derived.update({
        "encode_us": encode_us,
        "one_shot_us": one_shot_us,
        "encode_samples_us": encode_samples,
        "one_shot_samples_us": one_shot_samples,
        "retained_timer_window_us": {
            "required_floor": MIN_RETAINED_TIMER_WINDOW_US,
            "encode_execution_min": encode_minimum,
            "one_shot_encode_min": one_shot_minimum,
        },
        "workload_digests": digests,
        "kernel_available": available,
        "kernel_qualified": available,
        "selector_selected": selected,
        "observed_call_count": 2 if selected else 0,
        "source": {
            "source_commit": source_commit,
            "source_tree": source_tree,
            "source_tracked_dirty": False,
        },
        "raw": raw,
    })
    require(launch == derived, f"{label} launch record is not raw-derived")
    return launch


def validate_source_archive_member(
        member: tarfile.TarInfo, seen: set[str]) -> None:
    is_root = member.name in {"source", "source/"}
    require(member.name not in seen and
            (is_root or member.name.startswith("source/")) and
            not member.name.startswith("/") and
            ".." not in Path(member.name).parts,
            "source archive contains a duplicate or unsafe path")
    seen.add(member.name)
    require(member.isfile() or member.isdir() or member.issym(),
            "source archive contains an unsupported entry type")
    require(not is_root or member.isdir(),
            "source archive root is not a directory")
    if member.issym():
        require(not member.linkname.startswith("/") and
                ".." not in Path(member.linkname).parts,
                "source archive contains an unsafe symlink")


def validate_source_closure(report: dict[str, Any], controller: dict[str, str],
                            artifact_root: Path,
                            archive_only_source_closure: bool) \
        -> dict[str, Any]:
    pre = report["provenance"]
    post = report["post_identities"]
    exact_keys(pre, {
        "frozen_binary", "mutable_build_binary_pre", "runner_pre",
        "validator_pre", "source_archive_pre", "taskset_pre",
        "source_binding_pre",
    }, "pre-run provenance")
    exact_keys(post, {
        "frozen_binary", "mutable_build_binary", "runner", "validator",
        "source_archive", "taskset", "source_binding", "canonical_lock",
        "retained_raw_evidence",
    }, "post-run provenance")
    pairs = (("frozen_binary", "frozen_binary"),
             ("runner_pre", "runner"),
             ("validator_pre", "validator"),
             ("source_archive_pre", "source_archive"),
             ("taskset_pre", "taskset"))
    for before, after in pairs:
        require(pre[before] == post[after],
                f"pre/post identity changed: {before}")
    require(pre["mutable_build_binary_pre"] == post["mutable_build_binary"],
            "mutable build binary identity changed during campaign")
    require(pre["source_binding_pre"] == post["source_binding"],
            "source binding changed during campaign")
    frozen_specs = (
        (pre["frozen_binary"], "frozen benchmark", MAX_RESULT_BYTES, True),
        (pre["runner_pre"], "frozen runner", MAX_VALIDATOR_BYTES, True),
        (pre["validator_pre"], "frozen validator", MAX_VALIDATOR_BYTES, True),
        (pre["source_archive_pre"], "source archive",
         MAX_SOURCE_ARCHIVE_BYTES, True),
    )
    for identity, label, maximum, read_only in frozen_specs:
        observed = reopen_identity(
            identity, label, maximum, read_only=read_only)
        require(Path(observed["path"]).parent == artifact_root,
                f"{label} escaped the artifact root")
    binary_mode = pre["frozen_binary"]["mode"]
    runner_mode = pre["runner_pre"]["mode"]
    require(binary_mode & 0o111 != 0 and runner_mode & 0o111 != 0,
            "frozen binary or runner is not executable")
    taskset = reopen_identity(
        pre["taskset_pre"], "taskset executable", 64 * 1024 * 1024)
    require(taskset["path"] == "/usr/bin/taskset" and
            taskset["mode"] & 0o111 != 0,
            "taskset executable provenance changed")
    lock_post = exact_keys(post["canonical_lock"],
                           {"device", "inode", "links", "affinity"},
                           "post-run canonical lock")
    require(lock_post["device"] == report["canonical_lock"]["device"] and
            lock_post["inode"] == report["canonical_lock"]["inode"] and
            lock_post["links"] == 1 and
            lock_post["affinity"] == [EXPECTED_CPU],
            "post-run canonical lock continuity changed")
    require(report["binary"] == pre["frozen_binary"]["path"] and
            report["binary_sha256_pre"] == pre["frozen_binary"]["sha256"] ==
                report["binary_sha256_post"] and
            controller["--binary-sha256"] == pre["frozen_binary"]["sha256"] and
            controller["--runner-sha256"] == pre["runner_pre"]["sha256"] and
            controller["--validator-sha256"] == pre["validator_pre"]["sha256"] and
            controller["--source-archive-sha256"] ==
                pre["source_archive_pre"]["sha256"],
            "controller/report provenance hashes changed")
    require(controller["--build-binary"] ==
            pre["mutable_build_binary_pre"]["path"],
            "controller mutable build-binary path changed")
    for flag, identity in (("--binary", pre["frozen_binary"]),
                           ("--validator", pre["validator_pre"]),
                           ("--source-archive", pre["source_archive_pre"])):
        require(controller[flag] == identity["path"],
                f"controller {flag} path changed")

    binding = exact_keys(pre["source_binding_pre"], {
        "repository", "git", "head", "tree", "status_porcelain",
        "committed_frozen_sources", "canonical_archive",
    }, "source binding")
    require(binding["head"] == report["source_commit"] ==
            controller["--source-commit"] and
            binding["tree"] == report["source_tree"] ==
            controller["--source-tree"] and
            report["source_tracked_dirty"] is False and
            binding["status_porcelain"] == "",
            "source commit/tree binding changed")
    archive_bytes, _ = reopen_identity(
        pre["source_archive_pre"], "source archive",
        MAX_SOURCE_ARCHIVE_BYTES, read_only=True, return_bytes=True)
    archive = binding["canonical_archive"]
    exact_keys(archive, {"sha256", "size", "format"},
               "canonical archive metadata")
    require(archive["sha256"] == pre["source_archive_pre"]["sha256"] and
            archive["size"] == pre["source_archive_pre"]["size"] and
            archive["format"] == "git archive --format=tar --prefix=source/",
            "canonical archive metadata changed")
    members: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(archive_bytes), mode="r:") as stream:
        comment = stream.pax_headers.get("comment")
        require(comment == report["source_commit"],
                "archive commit comment changed")
        seen: set[str] = set()
        for member in stream:
            validate_source_archive_member(member, seen)
            if member.isfile() and member.name in {
                    "source/experiments/leopard2/gfni_t16/"
                    "run_k16r16_b64_avx512_gfni_abba.py",
                    "source/tools/leopard2_benchmark_json_test.py"}:
                source = stream.extractfile(member)
                require(source is not None, "cannot read frozen source member")
                members[member.name] = source.read(MAX_VALIDATOR_BYTES + 1)
                require(len(members[member.name]) <= MAX_VALIDATOR_BYTES,
                        "frozen source member exceeds its bound")
    committed = binding["committed_frozen_sources"]
    expected_relative = {
        "experiments/leopard2/gfni_t16/"
        "run_k16r16_b64_avx512_gfni_abba.py": pre["runner_pre"],
        "tools/leopard2_benchmark_json_test.py": pre["validator_pre"],
    }
    require(type(committed) is dict and set(committed) == set(expected_relative),
            "committed frozen source census changed")
    for relative, identity in expected_relative.items():
        metadata = exact_keys(committed[relative], {"sha256", "size"},
                              f"committed source {relative}")
        data = members.get("source/" + relative)
        require(data is not None and hashlib.sha256(data).hexdigest() ==
                metadata["sha256"] == identity["sha256"] and
                len(data) == metadata["size"] == identity["size"],
                f"source archive does not bind frozen {relative}")

    # Acquisition already binds and rechecks the live repository before and
    # after timing.  Archive-only replay deliberately avoids making durable
    # verification depend on the later state of that mutable checkout; the
    # canonical source archive and committed controller/validator bytes remain
    # fully checked above.
    repository = Path(binding["repository"])
    require(controller["--repository"] == str(repository),
            "controller repository path changed")
    live_checked = False
    if repository.exists() and not archive_only_source_closure:
        canonical_path(repository / ".git", must_exist=True)
        git_identity = reopen_identity(
            binding["git"], "Git executable", 64 * 1024 * 1024)
        git = git_identity["path"]

        def git_output(arguments: list[str], label: str,
                       maximum: int = MAX_SOURCE_ARCHIVE_BYTES) -> bytes:
            completed = subprocess.run(
                [git, "-C", str(repository)] + arguments,
                stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, check=False, timeout=120,
                env=GIT_ENVIRONMENT)
            require(completed.returncode == 0 and not completed.stderr and
                    len(completed.stdout) <= maximum,
                    f"{label} failed or exceeded its bound")
            return completed.stdout

        require(git_output(["rev-parse", "HEAD"], "Git HEAD").decode().strip() ==
                report["source_commit"] and
                git_output(["rev-parse", "HEAD^{tree}"],
                           "Git tree").decode().strip() == report["source_tree"] and
                git_output(["status", "--porcelain=v1", "--untracked-files=normal"],
                           "Git status") == b"",
                "live source no longer matches the frozen source")
        regenerated = git_output(
            ["archive", "--format=tar", "--prefix=source/",
             report["source_commit"]], "Git archive")
        require(regenerated == archive_bytes,
                "live canonical Git archive differs byte-for-byte")
        live_checked = True
    return {
        "source_commit": report["source_commit"],
        "source_tree": report["source_tree"],
        "source_archive_sha256": pre["source_archive_pre"]["sha256"],
        "frozen_binary_sha256": pre["frozen_binary"]["sha256"],
        "runner_sha256": pre["runner_pre"]["sha256"],
        "validator_sha256": pre["validator_pre"]["sha256"],
        "live_repository_checked": live_checked,
        "source_closure_mode": (
            "archive-only" if archive_only_source_closure else
            "archive-plus-live-when-present"),
    }


def replay(report_path: Path, journal_path: Path,
           *, archive_only_source_closure: bool = False) -> dict[str, Any]:
    auditor_identity = stable_file(
        Path(__file__).resolve(), maximum_bytes=MAX_VALIDATOR_BYTES,
        single_link=True)
    report, report_identity, _ = strict_json_file(
        report_path, MAX_REPORT_BYTES, "campaign report")
    require(type(report) is dict, "campaign report is not an object")
    require(report_identity["mode"] in {0o600, 0o400},
            "campaign report permissions changed")
    report_keys = {
        "schema", "claim_scope", "binary", "binary_sha256_pre",
        "controller_command", "provenance", "canonical_lock",
        "source_commit", "source_tree", "source_tracked_dirty", "cpu",
        "sibling", "affinity", "iterations_per_launch", "warmup_per_launch",
        "child_environment", "git_environment", "reuse_per_sample",
        "retained_timer_duration_floor_us", "launches_per_round",
        "max_round_attempts", "target_rounds", "inactive_control_rounds",
        "benchmark_schema_required", "diagnostic_option", "co_primary_gate",
        "started_unix_ns", "platform", "quiet_presample",
        "active_incomplete_cell", "cells", "status", "claim_passed",
        "binary_sha256_post", "post_identities", "journal_closure",
        "observed_sibling_list", "completed_unix_ns", "completed_cell_count",
        "total_cell_count", "gate_policy", "gate_results",
    }
    exact_keys(report, report_keys, "complete report")
    require(report["schema"] == REPORT_SCHEMA and
            report["claim_scope"] ==
                "same-binary mode-0 control versus mode-1 candidate at the "
                "exact K16/R16/B64/high/GF8/AUTO cell; no exact-main claim" and
            report["status"] == "complete" and
            type(report["claim_passed"]) is bool and
            report["active_incomplete_cell"] is None and
            report["cpu"] == EXPECTED_CPU and
            report["sibling"] == EXPECTED_SIBLING and
            report["affinity"] == [EXPECTED_CPU] and
            report["observed_sibling_list"] ==
                f"{EXPECTED_CPU},{EXPECTED_SIBLING}" and
            report["iterations_per_launch"] == ITERATIONS and
            report["warmup_per_launch"] == WARMUP and
            report["reuse_per_sample"] == REUSE and
            report["retained_timer_duration_floor_us"] ==
                MIN_RETAINED_TIMER_WINDOW_US and
            report["launches_per_round"] == 4 and
            report["max_round_attempts"] == MAX_ROUND_ATTEMPTS and
            report["target_rounds"] == TARGET_ROUNDS and
            report["inactive_control_rounds"] == INACTIVE_ROUNDS and
            report["benchmark_schema_required"] == BENCHMARK_SCHEMA and
            report["diagnostic_option"] == DIAGNOSTIC_OPTION and
            report["child_environment"] == CHILD_ENVIRONMENT and
            report["git_environment"] == GIT_ENVIRONMENT and
            report["co_primary_gate"] == {
                "ordinary_encode_execution_lcb95_min": 1.05,
                "one_shot_encode_lcb95_min": 1.05,
            } and
            report["gate_policy"] == {
                "target": (
                    "ordinary encode_execution and one_shot_encode are "
                    "co-primary; both 25-round paired-log t-interval 95% "
                    "lower bounds must be >= 1.05"),
                "route": (
                    "the target selects exactly in mode 1 with two untimed "
                    "route-probe calls; mode 0 and every inactive cell select "
                    "false with zero calls"),
                "digests": (
                    "all retained workload digests must be identical across "
                    "modes, order slots, attempts, and rounds within a cell"),
                "timer_duration": (
                    "every retained encode_execution and one_shot_encode "
                    "sample multiplied by reuse must span at least 250 "
                    "microseconds"),
                "inactive_controls": (
                    "timing is retained as descriptive evidence only and is "
                    "never a promotion gate"),
            } and
            type(report["platform"]) is str and bool(report["platform"]) and
            type(report["started_unix_ns"]) is int and
            type(report["completed_unix_ns"]) is int and
            report["completed_unix_ns"] >= report["started_unix_ns"],
            "frozen report contract changed")
    require(type(report["source_commit"]) is str and
            len(report["source_commit"]) == 40 and
            type(report["source_tree"]) is str and
            len(report["source_tree"]) == 40 and
            all(item in "0123456789abcdef"
                for item in report["source_commit"] + report["source_tree"]),
            "source identities are malformed")
    artifact_root = report_path.parent
    canonical_path(report_path, must_exist=True)
    canonical_path(journal_path, must_exist=True)
    require(journal_path.parent == artifact_root,
            "journal escaped the report artifact root")
    controller = parse_controller_command(report, report_path, journal_path)
    lock = validate_canonical_lock(report["canonical_lock"], controller)
    require(controller["--cpu"] == str(EXPECTED_CPU) and
            controller["--sibling"] == str(EXPECTED_SIBLING),
            "controller CPU topology changed")
    invocations = Path(controller["--invocations"])
    require(invocations.is_absolute() and invocations.parent == artifact_root and
            invocations.resolve(strict=True) == invocations and
            stat.S_ISDIR(invocations.lstat().st_mode) and
            not invocations.is_symlink(),
            "invocation root escaped or is aliased")
    require(stat.S_IMODE(invocations.lstat().st_mode) in {0o700, 0o500},
            "invocation root permissions changed")
    binary = Path(report["binary"])
    source_closure = validate_source_closure(
        report, controller, artifact_root, archive_only_source_closure)
    validator_identity = report["provenance"]["validator_pre"]
    validator = load_validator(Path(validator_identity["path"]),
                               validator_identity)
    events, journal_identity, journal_bytes = parse_journal(journal_path)
    require(journal_identity["mode"] in {0o600, 0o400},
            "campaign journal permissions changed")
    journal_expected = validate_identity(
        report["journal_closure"], "journal closure")
    require(identity_matches_acquisition_or_seal(
                journal_expected, journal_identity),
            "journal closure identity changed outside the sealing policy")
    quiet_passed = isolation_passed(report["quiet_presample"],
                                    "quiet presample")
    require(report["quiet_presample"]["reserved_sibling"]["total"] > 0 and
            report["quiet_presample"]["reserved_sibling"]["idle"] > 0 and
            report["quiet_presample"]["reserved_sibling"]["nonidle"] == 0,
            "quiet presample sibling was not idle")
    del quiet_passed  # benchmark CPU activity is not required in this sample.

    expected_events: list[dict[str, Any]] = [{
        "event": "campaign_started",
        "binary_sha256": report["provenance"]["frozen_binary"]["sha256"],
        "runner_sha256": report["provenance"]["runner_pre"]["sha256"],
        "validator_sha256": report["provenance"]["validator_pre"]["sha256"],
        "source_archive_sha256":
            report["provenance"]["source_archive_pre"]["sha256"],
        "source_commit": report["source_commit"],
        "source_tree": report["source_tree"],
        "controller_command": report["controller_command"],
        "canonical_lock": report["canonical_lock"],
        "started_unix_ns": report["started_unix_ns"],
    }, {
        "event": "quiet_presample",
        "isolation": report["quiet_presample"],
    }]
    cells = report["cells"]
    require(type(cells) is list and len(cells) == len(CELLS) and
            report["completed_cell_count"] == report["total_cell_count"] ==
                len(CELLS),
            "cell completion census changed")
    raw_manifest: list[dict[str, Any]] = []
    replayed_cells: list[dict[str, Any]] = []
    accepted_launches = rejected_launches = rejected_attempt_count = 0
    all_paths: set[str] = set()
    route_gates: dict[str, bool] = {}
    digest_gates: dict[str, bool] = {}
    timer_gates: dict[str, bool] = {}
    for cell_index, (cell, expected_cell) in enumerate(zip(cells, CELLS)):
        cell_id, k, r, shard_bytes, rounds, role, profile, field, backend, \
            rationale = expected_cell
        cell_keys = {
            "id", "K", "R", "shard_bytes", "role", "profile", "field",
            "backend", "rationale", "rounds", "reuse", "workload_digests",
            "encode_execution", "one_shot_encode", "records",
            "rejected_contaminated_attempts",
        }
        exact_keys(cell, cell_keys, f"cell {cell_id}")
        require((cell["id"], cell["K"], cell["R"], cell["shard_bytes"],
                 cell["rounds"], cell["role"], cell["profile"],
                 cell["field"], cell["backend"], cell["rationale"],
                 cell["reuse"]) ==
                (cell_id, k, r, shard_bytes, rounds, role, profile, field,
                 backend, rationale, REUSE),
                f"cell contract changed: {cell_id}")
        records = cell["records"]
        rejected = cell["rejected_contaminated_attempts"]
        require(type(records) is list and len(records) == rounds and
                type(rejected) is list,
                f"cell attempts changed: {cell_id}")
        accepted_by_round: dict[int, dict[str, Any]] = {}
        rejected_by_round: dict[int, list[dict[str, Any]]] = {
            index: [] for index in range(rounds)}
        for record in records:
            require(type(record) is dict and type(record.get("round")) is int and
                    record["round"] not in accepted_by_round and
                    0 <= record["round"] < rounds,
                    f"accepted round identity changed: {cell_id}")
            accepted_by_round[record["round"]] = record
        for attempt in rejected:
            require(type(attempt) is dict and type(attempt.get("round")) is int and
                    0 <= attempt["round"] < rounds,
                    f"rejected round identity changed: {cell_id}")
            rejected_by_round[attempt["round"]].append(attempt)
        require(set(accepted_by_round) == set(range(rounds)),
                f"accepted round coverage changed: {cell_id}")
        encode_logs: list[float] = []
        one_shot_logs: list[float] = []
        cell_launches: list[dict[str, Any]] = []
        seed = 0x16164000 + cell_index
        for round_index in range(rounds):
            accepted = accepted_by_round[round_index]
            accepted_attempt = accepted.get("attempt")
            require(type(accepted_attempt) is int and
                    0 <= accepted_attempt < MAX_ROUND_ATTEMPTS,
                    f"accepted attempt index changed: {cell_id}/{round_index}")
            rejected_round = rejected_by_round[round_index]
            require([item.get("attempt") for item in rejected_round] ==
                    list(range(accepted_attempt)),
                    f"rejected attempt prefix changed: {cell_id}/{round_index}")
            order = [1, 0, 0, 1] if round_index % 2 == 0 else [0, 1, 1, 0]
            for attempt_index, attempt in enumerate(rejected_round + [accepted]):
                expected_attempt_keys = {
                    "round", "attempt", "order", "launches", "isolation"}
                if attempt is accepted:
                    expected_attempt_keys |= {
                        "log_control_over_candidate_encode",
                        "log_control_over_candidate_one_shot"}
                exact_keys(attempt, expected_attempt_keys,
                           f"attempt {cell_id}/{round_index}/{attempt_index}")
                require(attempt["round"] == round_index and
                        attempt["attempt"] == attempt_index and
                        attempt["order"] == order and
                        type(attempt["launches"]) is list and
                        len(attempt["launches"]) == 4,
                        f"attempt order changed: {cell_id}/{round_index}")
                accepted_isolation = isolation_passed(
                    attempt["isolation"],
                    f"isolation {cell_id}/{round_index}/{attempt_index}")
                require(accepted_isolation is (attempt is accepted),
                        f"attempt disposition differs from isolation: {cell_id}")
                launches: list[dict[str, Any]] = []
                for slot, (launch, mode) in enumerate(zip(
                        attempt["launches"], order)):
                    label = (f"{cell_id}-round{round_index}-attempt"
                             f"{attempt_index}-slot{slot}-mode{mode}")
                    validated = validate_launch(
                        launch, cell=expected_cell, mode=mode, label=label,
                        seed=seed, artifact_root=artifact_root,
                        invocations=invocations, binary=binary,
                        source_commit=report["source_commit"],
                        source_tree=report["source_tree"], validator=validator,
                        raw_manifest=raw_manifest)
                    for identity in validated["raw"].values():
                        require(identity["path"] not in all_paths,
                                "raw evidence path is reused")
                        all_paths.add(identity["path"])
                    expected_events.extend(({
                        "event": "launch_started", "label": label,
                        "cell": cell_id, "mode": mode,
                        "command": expected_benchmark_command(
                            binary, Path(validated["raw"]["result"]["path"]),
                            expected_cell, mode, seed),
                        "invocation": str(invocations / label),
                    }, {
                        "event": "launch_validated", "label": label,
                        "record": validated,
                    }))
                    launches.append(validated)
                    cell_launches.append(validated)
                expected_events.append({
                    "event": "round_attempt_complete", "cell": cell_id,
                    "round": round_index, "attempt": attempt_index,
                    "isolation": attempt["isolation"],
                    "accepted": attempt is accepted,
                })
                if attempt is not accepted:
                    rejected_attempt_count += 1
                    rejected_launches += 4
                    continue
                accepted_launches += 4
                controls = [item for item in launches if item["mode"] == 0]
                candidates = [item for item in launches if item["mode"] == 1]
                require(len(controls) == len(candidates) == 2,
                        "ABBA round is unbalanced")
                encode_log = statistics.mean(
                    math.log(item["encode_us"]) for item in controls) - \
                    statistics.mean(
                        math.log(item["encode_us"]) for item in candidates)
                one_shot_log = statistics.mean(
                    math.log(item["one_shot_us"]) for item in controls) - \
                    statistics.mean(
                        math.log(item["one_shot_us"]) for item in candidates)
                require(attempt["log_control_over_candidate_encode"] == encode_log and
                        attempt["log_control_over_candidate_one_shot"] ==
                            one_shot_log,
                        f"reported ABBA contrast changed: {cell_id}")
                encode_logs.append(encode_log)
                one_shot_logs.append(one_shot_log)
                expected_events.append({
                    "event": "round_complete", "cell": cell_id,
                    "round": round_index, "encode_log": encode_log,
                    "one_shot_log": one_shot_log,
                })
        all_digests = [item["workload_digests"] for item in cell_launches]
        require(all_digests and all(item == cell["workload_digests"]
                                    for item in all_digests),
                f"cross-mode digest identity failed: {cell_id}")
        encode_summary = summarize(encode_logs, role == "target")
        one_shot_summary = summarize(one_shot_logs, role == "target")
        require(cell["encode_execution"] == encode_summary and
                cell["one_shot_encode"] == one_shot_summary,
                f"cell summary is not retained-round-derived: {cell_id}")
        route_gate = all(
            item["selector_selected"] is
                (role == "target" and item["mode"] == 1) and
            item["observed_call_count"] ==
                (2 if role == "target" and item["mode"] == 1 else 0)
            for item in cell_launches)
        timer_gate = all(
            item["retained_timer_window_us"]["encode_execution_min"] >=
                MIN_RETAINED_TIMER_WINDOW_US and
            item["retained_timer_window_us"]["one_shot_encode_min"] >=
                MIN_RETAINED_TIMER_WINDOW_US for item in cell_launches)
        route_gates[cell_id] = route_gate
        digest_gates[cell_id] = True
        timer_gates[cell_id] = timer_gate
        expected_events.append({
            "event": "cell_complete", "cell": cell_id,
            "encode_execution": encode_summary,
            "one_shot_encode": one_shot_summary,
        })
        replayed_cells.append({
            "id": cell_id,
            "role": role,
            "accepted_rounds": rounds,
            "rejected_contaminated_attempts":
                sum(len(items) for items in rejected_by_round.values()),
            "encode_execution": encode_summary,
            "one_shot_encode": one_shot_summary,
            "route_gate": route_gate,
            "digest_gate": True,
            "timer_window_gate": timer_gate,
        })

    manifest_records = []
    for cell in cells:
        for attempts, disposition in ((cell["records"], "accepted"),
                                      (cell["rejected_contaminated_attempts"],
                                       "rejected")):
            for attempt in attempts:
                for launch in attempt["launches"]:
                    for kind in ("command", "stdout", "stderr", "result"):
                        identity = launch["raw"][kind]
                        manifest_records.append({
                            "launch": launch["label"],
                            "disposition": disposition,
                            "kind": kind,
                            "path": identity["path"],
                            "sha256": identity["sha256"],
                            "size": identity["size"],
                        })
    manifest_hash = hashlib.sha256(json.dumps(
        manifest_records, sort_keys=True,
        separators=(",", ":")).encode()).hexdigest()
    raw_closure = {
        "accepted_launches": accepted_launches,
        "rejected_contaminated_launches": rejected_launches,
        "file_count": len(manifest_records),
        "manifest_sha256": manifest_hash,
        "manifest_semantics": (
            "SHA-256 of canonical JSON records ordered by cell, disposition, "
            "attempt, launch, and command/stdout/stderr/result"),
    }
    require(raw_closure == report["post_identities"]["retained_raw_evidence"],
            "retained raw evidence closure changed")
    target = replayed_cells[0]
    ordinary_speed = target["encode_execution"]["ci95"][0] >= 1.05
    one_shot_speed = target["one_shot_encode"]["ci95"][0] >= 1.05
    inactive_routes = {key: value for key, value in route_gates.items()
                       if key != CELLS[0][0]}
    route_safety = route_gates[CELLS[0][0]] and all(inactive_routes.values())
    digest_safety = all(digest_gates.values())
    timer_safety = all(timer_gates.values())
    gates = {
        "target_ordinary_encode_lcb95_at_least_1_05": ordinary_speed,
        "target_one_shot_encode_lcb95_at_least_1_05": one_shot_speed,
        "target_route_and_call_count": route_gates[CELLS[0][0]],
        "inactive_route_and_call_count": inactive_routes,
        "per_cell_cross_mode_digest_identity": digest_gates,
        "per_cell_retained_timer_duration_floor": timer_gates,
        "global_route_safety": route_safety,
        "global_digest_safety": digest_safety,
        "global_retained_timer_duration_safety": timer_safety,
        "inactive_timing_used_as_gate": False,
    }
    claim_passed = (ordinary_speed and one_shot_speed and route_safety and
                    digest_safety and timer_safety)
    require(gates == report["gate_results"] and
            claim_passed is report["claim_passed"],
            "recomputed promotion gates differ from report")
    expected_events.append({
        "event": "campaign_complete",
        "claim_passed": claim_passed,
        "gate_results": gates,
        "post_identities": report["post_identities"],
        "retained_raw_evidence": raw_closure,
    })
    require(events == expected_events,
            "journal event order or event bytes differ from complete replay")
    require(len(raw_manifest) == len(manifest_records) and
            len(all_paths) == len(manifest_records),
            "raw replay file census changed")
    return {
        "schema": AUDIT_SCHEMA,
        "status": "complete",
        "audit_passed": True,
        "timing_performed": False,
        "benchmark_executed": False,
        "auditor": {
            "path": auditor_identity["path"],
            "sha256": auditor_identity["sha256"],
            "size": auditor_identity["size"],
        },
        "report": {
            "path": str(report_path),
            "sha256": report_identity["sha256"],
            "size": report_identity["size"],
        },
        "journal": {
            "path": str(journal_path),
            "sha256": journal_identity["sha256"],
            "size": journal_identity["size"],
            "event_count": len(events),
            "canonical_jsonl_sha256": hashlib.sha256(journal_bytes).hexdigest(),
        },
        "contract": {
            "benchmark_schema": BENCHMARK_SCHEMA,
            "selector_contract": CONTRACT,
            "target_rounds": TARGET_ROUNDS,
            "inactive_rounds": INACTIVE_ROUNDS,
            "iterations": ITERATIONS,
            "warmup": WARMUP,
            "reuse": REUSE,
            "timer_window_floor_us": MIN_RETAINED_TIMER_WINDOW_US,
            "cpu": EXPECTED_CPU,
            "sibling": EXPECTED_SIBLING,
            "cell_count": len(CELLS),
        },
        "replay": {
            "validator_invocations": accepted_launches + rejected_launches,
            "accepted_launches": accepted_launches,
            "rejected_contaminated_launches": rejected_launches,
            "rejected_contaminated_attempts": rejected_attempt_count,
            "raw_file_count": len(manifest_records),
            "raw_manifest_sha256": manifest_hash,
            "cells": replayed_cells,
        },
        "source_and_binary_closure": source_closure,
        "gate_results": gates,
        "claim_passed": claim_passed,
    }


def write_canonical(path: Path, value: dict[str, Any]) -> None:
    canonical_path(path, must_exist=False)
    require(not path.exists(), "refusing to overwrite audit output")
    data = (json.dumps(value, sort_keys=True,
                       separators=(",", ":")) + "\n").encode()
    descriptor = os.open(
        path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise
    directory = os.open(path.parent, os.O_RDONLY | os.O_CLOEXEC)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def self_test() -> None:
    for payload in (b'{"x":1,"x":2}', b'{"x":NaN}', b'\xff'):
        try:
            strict_json_bytes(payload, "adversarial JSON")
        except (RuntimeError, UnicodeDecodeError, json.JSONDecodeError):
            pass
        else:
            raise RuntimeError("strict JSON parser accepted adversarial input")
    logs = [math.log(2.0 + index / 100.0) for index in range(TARGET_ROUNDS)]
    result = summarize(logs, True)
    require(result["rounds"] == TARGET_ROUNDS and
            result["ci95"][0] < result["geometric_mean_speedup"] <
            result["ci95"][1], "Student-t self-test failed")
    archive_bytes = io.BytesIO()
    with tarfile.open(fileobj=archive_bytes, mode="w") as archive:
        root_member = tarfile.TarInfo("source/")
        root_member.type = tarfile.DIRTYPE
        archive.addfile(root_member)
        file_member = tarfile.TarInfo("source/file")
        file_member.size = 0
        archive.addfile(file_member, io.BytesIO())
    archive_seen: set[str] = set()
    with tarfile.open(fileobj=io.BytesIO(archive_bytes.getvalue()), mode="r:") \
            as archive:
        for member in archive:
            validate_source_archive_member(member, archive_seen)
    require(archive_seen == {"source", "source/file"},
            "canonical source archive root self-test failed")
    adversarial_members = []
    root_file = tarfile.TarInfo("source")
    adversarial_members.append(root_file)
    root_slash_file = tarfile.TarInfo("source/")
    adversarial_members.append(root_slash_file)
    outside_prefix = tarfile.TarInfo("source-other/file")
    adversarial_members.append(outside_prefix)
    traversal = tarfile.TarInfo("source/../outside")
    adversarial_members.append(traversal)
    unsafe_symlink = tarfile.TarInfo("source/link")
    unsafe_symlink.type = tarfile.SYMTYPE
    unsafe_symlink.linkname = "../outside"
    adversarial_members.append(unsafe_symlink)
    for member in adversarial_members:
        try:
            validate_source_archive_member(member, set())
        except RuntimeError:
            pass
        else:
            raise RuntimeError(
                "source archive self-test accepted an unsafe member")
    duplicate_seen = {"source/file"}
    try:
        validate_source_archive_member(file_member, duplicate_seen)
    except RuntimeError:
        pass
    else:
        raise RuntimeError("source archive self-test accepted a duplicate")
    inherited_lock = {
        "path": "/tmp/leopard-gf8-authoritative.lock",
        "mode": "inherited-across-build-copy-campaign",
        "scope": "wrapper-build-copy-campaign",
        "descriptor": 9,
        "device": 1,
        "inode": 2,
    }
    validate_canonical_lock(inherited_lock, {"--lock-fd": "9"})
    validate_canonical_lock({
        **inherited_lock,
        "mode": "runner-acquired",
        "scope": "campaign-only",
    }, {})
    for adversarial_lock in (
            {key: value for key, value in inherited_lock.items()
             if key != "scope"},
            {**inherited_lock, "scope": "campaign-only"},
            {**inherited_lock, "mode": "runner-acquired"}):
        try:
            validate_canonical_lock(adversarial_lock, {"--lock-fd": "9"})
        except RuntimeError:
            pass
        else:
            raise RuntimeError("canonical lock self-test accepted a mismatch")
    with tempfile.TemporaryDirectory(prefix="leopard-t16-audit-self-test.") \
            as directory:
        root = Path(directory).resolve()
        path = root / "canonical.json"
        write_canonical(path, {"b": 2, "a": 1})
        require(path.read_bytes() == b'{"a":1,"b":2}\n',
                "canonical writer self-test failed")
        source = root / "source"
        source.write_bytes(b"xx")
        hardlink = root / "hardlink"
        os.link(source, hardlink)
        symlink = root / "symlink"
        symlink.symlink_to(source)
        adversarial_files = (
            (hardlink, 8, "hard link"),
            (symlink, 8, "symbolic link"),
            (source, 1, "oversized file"),
            (root / ".." / root.name / "source", 8, "aliased path"),
        )
        for candidate, maximum, label in adversarial_files:
            try:
                stable_file(candidate, maximum_bytes=maximum)
            except RuntimeError:
                pass
            else:
                raise RuntimeError(f"self-test accepted {label}")
        journal = root / "journal.jsonl"
        journal.write_bytes(b'{ "event": "not-canonical" }\n')
        try:
            parse_journal(journal)
        except RuntimeError:
            pass
        else:
            raise RuntimeError("self-test accepted non-canonical journal")
        seal = root / "seal"
        seal.write_bytes(b"sealed")
        os.chmod(seal, 0o600)
        acquired = stable_file(seal, maximum_bytes=16)
        os.chmod(seal, 0o400)
        sealed = stable_file(seal, maximum_bytes=16)
        require(identity_matches_acquisition_or_seal(acquired, sealed),
                "monotonic seal self-test failed")
        tampered = dict(sealed, mtime_ns=sealed["mtime_ns"] + 1)
        require(not identity_matches_acquisition_or_seal(acquired, tampered),
                "monotonic seal accepted an mtime change")
    print("K16/R16/B64 GFNI no-timing auditor self-test passed")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path)
    parser.add_argument("--journal", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--archive-only-source-closure", action="store_true")
    options = parser.parse_args()
    if options.self_test:
        require(options.report is None and options.journal is None and
                options.output is None and
                not options.archive_only_source_closure,
                "--self-test cannot be combined with replay paths")
        self_test()
        return 0
    require(options.report is not None and options.journal is not None and
            options.output is not None,
            "--report, --journal, and --output are required")
    for path in (options.report, options.journal, options.output):
        require(path.is_absolute(), "replay paths must be absolute")
    require(options.report.parent == options.journal.parent,
            "report and journal must share one artifact directory")
    require(options.output != options.report and options.output != options.journal,
            "audit output aliases an input artifact")
    result = replay(
        options.report, options.journal,
        archive_only_source_closure=options.archive_only_source_closure)
    write_canonical(options.output, result)
    print(json.dumps({
        "audit_passed": True,
        "claim_passed": result["claim_passed"],
        "output": str(options.output),
        "report_sha256": result["report"]["sha256"],
        "journal_sha256": result["journal"]["sha256"],
        "timing_performed": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
