#!/usr/bin/env python3
"""Run and replay one authoritative, full-matrix C7 AVX2 benchmark.

This tool never builds C7.  It binds an already-built standalone executable and
its core archive to an exact core commit beneath a clean tooling revision,
validates the final v4 A/B build provenance live, holds both a coordinator
reservation and the shared per-user physical-core lease, pins the child to one
SMT thread, and rejects the run when the sibling records any non-idle work.
``verify`` and ``verify-failure`` are portable: they replay canonical retained
bytes and every derived value without requiring the original source, build,
executable, or reservation paths.  A verified failure remains a failure.
``self-test`` is synthetic and never executes benchmark timing.
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
import platform
import re
import resource
import signal
import socket
import stat
import statistics
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path
from typing import Any, Mapping, Sequence


RAW_SCHEMA = "leopard2-c7-authoritative-raw/v3"
MANIFEST_SCHEMA = "leopard2-c7-authoritative-manifest/v3"
FAILURE_SCHEMA = "leopard2-c7-authoritative-failure/v4"
FAILURE_STATE_SCHEMA = "leopard2-c7-authoritative-failure-state/v2"
PUBLICATION_STATE_SCHEMA = "leopard2-c7-authoritative-publication-state/v1"
LIFECYCLE_SCHEMA = "leopard2-c7-authoritative-lifecycle/v1"
INPUT_SCHEMA = "leopard2-c7-authoritative-input/v2"
BUILD_PROVENANCE_SCHEMA = "leopard2-c7-authoritative-build-provenance/v2"
RESERVATION_SCHEMA = "leopard2-cpu-reservation/v1"
PAIR_LEASE_SCHEMA = "leopard2-cpu-pair-lease/v1"
KERNEL_LEASE_SCHEMA = "leopard2-linux-abstract-lease/v1"
ISOLATION_SCHEMA = "leopard2-c7-authoritative-isolation/v1"
BUILD_MANIFEST_SCHEMA = "leopard2-c7-build-run-manifest/v4"
BUILD_ATTESTATION_SCHEMA = "leopard2-c7-authoritative-build-attestation/v1"
CHILD_SCHEMA = "leopard2-c7-exact-low/v1"
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
GIT_SHA_RE = re.compile(r"[0-9a-f]{40}\Z")
SAMPLE_COUNT = 7
MAX_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_LOG_BYTES = 1024 * 1024
MAX_TOP_LEVEL_JSON_BYTES = 16 * 1024 * 1024
MAX_TIMEOUT_SECONDS = 24 * 60 * 60
MAX_CPU_ID = (1 << 20) - 1
MAX_CPU_SET_SIZE = 1 << 16
MAX_ARTIFACT_COUNT = 256
MAX_ARTIFACT_DEPTH = 16
MAX_ARTIFACT_PATH_BYTES = 4096
MAX_REPORTED_SCRATCH_BYTES = 1 << 40
MAX_DIAGNOSTIC_BYTES = 16 * 1024
CHILD_REAP_TIMEOUT_SECONDS = 5.0
SOURCE_RELATIVE = Path("experiments/leopard2/non_power_of_two/c7/c7_exact_low.cpp")
RUNNER_RELATIVE = Path("experiments/leopard2/non_power_of_two/c7/run_authoritative.py")
BUILD_RUNNER_RELATIVE = Path(
    "experiments/leopard2/non_power_of_two/c7/run_matrix.py")
BUILD_VALIDATOR_RELATIVE = Path(
    "experiments/leopard2/non_power_of_two/c7/validate_evidence.py")
BUILD_NAMES = ("scalar", "ssse3", "avx2", "auto", "asan-ubsan")
BUILD_SCOPE = (
    "correctness plus one affinity-selected non-authoritative harness "
    "smoke; no promotion timing")
CORE_SOURCE_CLOSURE = (
    "CMakeLists.txt", "Leopard2Backend.cpp", "Leopard2Backend.h",
    "Leopard2BackendAVX2.cpp", "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp",
    "Leopard2Direct.h", "Leopard2Dispatch.h", "Leopard2Plan.cpp",
    "Leopard2Plan.h", "LeopardCommon.cpp", "LeopardCommon.h",
    "LeopardFF16.cpp", "LeopardFF16.h", "LeopardFF8.cpp", "LeopardFF8.h",
    "cmake/leopardConfig.cmake.in", "leopard.cpp", "leopard.h",
    "leopard2.cpp", "leopard2.h",
)
CHILD_ENVIRONMENT = {
    "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1", "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE", "PATH": "/usr/bin:/bin", "TZ": "UTC",
}
VALIDATOR_ENVIRONMENT = {
    **CHILD_ENVIRONMENT, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0",
}
CPU_STAT_FIELDS = (
    "user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal",
)
CPU_STAT_IDLE_FIELDS = ("idle", "iowait")
FAILURE_STAGES = (
    "arguments", "request", "host", "reservation", "lease", "inputs",
    "attested", "child-captured", "isolated", "inputs-after", "host-after",
    "validated",
)
FAILURE_CONTEXT_FIELDS = (
    "arguments", "request", "host_before", "reservation", "pair_lease",
    "inputs_before", "build_provenance", "child", "isolation", "inputs_after",
    "host_after", "validated_output",
)
FAILURE_CODE_STAGE = {
    "arguments-invalid": -1,
    "request-invalid": 0,
    "host-capture-failed": 1,
    "reservation-failed": 2,
    "lease-failed": 3,
    "inputs-capture-failed": 4,
    "attestation-failed": 5,
    "child-capture-failed": 6,
    "isolation-capture-failed": 7,
    "child-timeout": 8,
    "child-exit": 8,
    "child-result-missing": 8,
    "child-stdout-invalid": 8,
    "child-stderr-invalid": 8,
    "isolation-rejected": 8,
    "reservation-post-child-invalid": 8,
    "lease-post-child-invalid": 8,
    "child-result-json-invalid": 8,
    "inputs-after-capture-failed": 8,
    "inputs-drift": 9,
    "host-after-capture-failed": 9,
    "host-drift": 10,
    "child-result-invalid": 10,
    "manifest-validation-failed": 11,
    "manifest-write-failed": 11,
    "guard-exit-failed": 11,
    "affinity-restore-failed": 11,
    "publication-state-write-failed": 11,
    "raw-write-failed": 11,
    "final-snapshot-invalid": 11,
}
FAILURE_DIAGNOSTICS = {
    "arguments-invalid": "authoritative arguments are invalid",
    "request-invalid": "authoritative request is invalid",
    "host-capture-failed": "authoritative host capture failed",
    "reservation-failed": "CPU reservation acquisition failed",
    "lease-failed": "CPU pair lease acquisition failed",
    "inputs-capture-failed": "authoritative input capture failed",
    "attestation-failed": "authoritative build attestation failed",
    "child-capture-failed": "C7 child capture failed",
    "isolation-capture-failed": "C7 isolation capture failed",
    "child-timeout": "C7 child timed out",
    "child-result-missing": "C7 child did not write result JSON",
    "child-stdout-invalid": "C7 child unexpectedly wrote stdout",
    "child-stderr-invalid": "C7 child stderr progress changed",
    "isolation-rejected": "C7 isolation policy rejected the child",
    "reservation-post-child-invalid":
        "post-child reservation guard validation failed",
    "lease-post-child-invalid": "post-child pair guard validation failed",
    "child-result-json-invalid": "C7 child result JSON is invalid",
    "inputs-after-capture-failed": "post-child input capture failed",
    "inputs-drift": "source/archive/executable changed during C7 run",
    "host-after-capture-failed": "post-child host capture failed",
    "host-drift": "host topology/frequency policy changed during C7 run",
    "child-result-invalid": "C7 child result validation failed",
    "guard-exit-failed": "authoritative guard teardown failed",
    "affinity-restore-failed": "original process affinity was not restored exactly",
    "publication-state-write-failed": "publication state staging failed",
    "raw-write-failed": "raw evidence staging failed",
    "manifest-validation-failed": "manifest validation failed",
    "manifest-write-failed": "terminal manifest publication failed",
    "final-snapshot-invalid": "final staged snapshot validation failed",
}

EXPECTED_PROFILE = {
    "family": 3, "version": 1, "coordinate_map": 1,
    "systematic": "0..K-1", "parity": "K..K+R-1",
    "production_enabled": False,
}
EXPECTED_CORRECTNESS = {
    "gf8_cases": 9, "gf16_cases": 5, "coefficients": 118717,
    "gf16_vandermonde_coefficients": 1500, "encode_executions": 117,
    "encode_symbol_comparisons": 1030423,
    "subset_encode_executions": 117, "decode_plans": 403,
    "decode_executions": 403, "decode_symbol_comparisons": 272487,
    "maximum_loss_plans": 117, "unavailable_parity_plans": 175,
    "no_loss_null_calls": 14, "parity_rebuilds": 403,
    "odd_gf16_rejections": 10, "overlap_rejections": 59,
    "parity_output_overlap_rejections": 13,
    "restored_output_overlap_rejections": 12,
    "restored_input_overlap_rejections": 20,
    "selected_parity_null_rejections": 14, "survivor_null_rejections": 6,
    "atomic_rejection_bytes_checked": 61570,
    "read_only_input_alias_calls": 13,
    "read_only_input_alias_symbol_comparisons": 2139,
    "decode_read_only_input_alias_calls": 117,
    "decode_read_only_input_alias_symbol_comparisons": 6025,
    "hot_path_allocations": 0, "digest_fnv64": "0xec4179e9f2776a58",
}

# K, R, bytes, batch, losses, exact field, padded field.
EXPECTED_CELLS = (
    (3, 253, 64, 8, 3, 1, 2),
    (3, 253, 1024, 8, 3, 1, 2),
    (3, 253, 65536, 1, 3, 1, 2),
    (16, 240, 1024, 8, 8, 1, 1),
    (64, 192, 65536, 1, 8, 1, 1),
    (100, 156, 1024, 8, 8, 1, 2),
    (127, 129, 65536, 1, 8, 1, 2),
    (128, 128, 1024, 8, 8, 1, 1),
    (192, 64, 1024, 8, 8, 1, 2),
    (248, 8, 65536, 1, 8, 1, 2),
    (129, 100, 1024, 8, 8, 2, 2),
    (1000, 200, 1024, 1, 8, 2, 2),
)
CELL_KEYS = {
    "K", "R", "batch", "bytes", "exact_coefficients", "exact_decode",
    "exact_decode_samples_us", "exact_decode_setup",
    "exact_decode_setup_samples_us", "exact_decode_terms", "exact_encode",
    "exact_encode_samples_us", "exact_field", "exact_setup",
    "exact_setup_samples_us", "losses", "padded_decode",
    "padded_decode_samples_us", "padded_decode_scratch",
    "padded_decode_setup", "padded_decode_setup_samples_us", "padded_encode",
    "padded_encode_samples_us", "padded_encode_scratch", "padded_field",
    "padded_setup", "padded_setup_samples_us",
}
SUMMARY_SAMPLE_PAIRS = (
    ("exact_setup", "exact_setup_samples_us"),
    ("padded_setup", "padded_setup_samples_us"),
    ("exact_decode_setup", "exact_decode_setup_samples_us"),
    ("padded_decode_setup", "padded_decode_setup_samples_us"),
    ("exact_encode", "exact_encode_samples_us"),
    ("padded_encode", "padded_encode_samples_us"),
    ("exact_decode", "exact_decode_samples_us"),
    ("padded_decode", "padded_decode_samples_us"),
)


class EvidenceError(RuntimeError):
    """The run or evidence is not authoritative."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def bounded_diagnostic(value: object, label: str) -> str:
    require(isinstance(value, str) and value and
            len(value.encode("utf-8", errors="replace")) <= MAX_DIAGNOSTIC_BYTES,
            f"{label} is not a bounded diagnostic")
    return value


def empty_lifecycle() -> dict[str, Any]:
    def guard() -> dict[str, Any]:
        return {"entered": False, "exit": "not-entered",
                "error_type": None, "error": None}

    def check() -> dict[str, Any]:
        return {"status": "not-reached", "error_type": None, "error": None}

    return {
        "schema": LIFECYCLE_SCHEMA,
        "reservation": guard(), "pair_lease": guard(),
        "guard_checks": {
            "reservation_post_child": check(), "pair_post_child": check(),
            "reservation_final": check(), "pair_final": check(),
        },
        "affinity": {
            "captured": None, "restore_attempted": False,
            "restore_succeeded": False, "observed": None, "exact": False,
            "error_type": None, "error": None,
        },
    }


def _bounded_exception(error: BaseException) -> tuple[str, str]:
    name = type(error).__name__[:128] or "Exception"
    message = str(error)
    encoded = message.encode("utf-8", errors="replace")
    if len(encoded) > MAX_DIAGNOSTIC_BYTES:
        encoded = encoded[:MAX_DIAGNOSTIC_BYTES]
        message = encoded.decode("utf-8", errors="replace")
    return name, message or name


def validate_lifecycle(value: object, *, require_success: bool) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "reservation", "pair_lease", "guard_checks", "affinity"} and
        value["schema"] == LIFECYCLE_SCHEMA,
        "authoritative lifecycle schema changed")
    for name in ("reservation", "pair_lease"):
        guard = value[name]
        require(isinstance(guard, dict) and set(guard) == {
            "entered", "exit", "error_type", "error"} and
            type(guard["entered"]) is bool and
            guard["exit"] in ("not-entered", "pass", "failed"),
            f"{name} lifecycle changed")
        if guard["entered"]:
            require(guard["exit"] in ("pass", "failed"),
                    f"{name} entered without teardown")
        else:
            require(guard["exit"] == "not-entered",
                    f"{name} exited without entering")
        if guard["exit"] == "failed":
            bounded_diagnostic(guard["error_type"], f"{name} exit type")
            bounded_diagnostic(guard["error"], f"{name} exit error")
        else:
            require(guard["error_type"] is None and guard["error"] is None,
                    f"{name} lifecycle has a spurious error")
    checks = value["guard_checks"]
    require(isinstance(checks, dict) and set(checks) == {
        "reservation_post_child", "pair_post_child",
        "reservation_final", "pair_final"},
        "guard-check lifecycle changed")
    for name, check in checks.items():
        require(isinstance(check, dict) and set(check) == {
            "status", "error_type", "error"} and
            check["status"] in ("not-reached", "pass", "failed"),
            f"{name} guard-check lifecycle changed")
        if check["status"] == "failed":
            bounded_diagnostic(check["error_type"], f"{name} failure type")
            bounded_diagnostic(check["error"], f"{name} failure error")
        else:
            require(check["error_type"] is None and check["error"] is None,
                    f"{name} has a spurious guard-check error")
    require(checks["pair_post_child"]["status"] == "not-reached" or
            checks["reservation_post_child"]["status"] == "pass",
            "pair post-child check bypassed reservation check")
    require(checks["reservation_final"]["status"] == "not-reached" or
            (checks["reservation_post_child"]["status"] == "pass" and
             checks["pair_post_child"]["status"] == "pass"),
            "final guard checks bypassed post-child checks")
    require(checks["pair_final"]["status"] == "not-reached" or
            checks["reservation_final"]["status"] == "pass",
            "final pair check bypassed final reservation check")
    affinity = value["affinity"]
    require(isinstance(affinity, dict) and set(affinity) == {
        "captured", "restore_attempted", "restore_succeeded", "observed",
        "exact", "error_type", "error"} and
        type(affinity["restore_attempted"]) is bool and
        type(affinity["restore_succeeded"]) is bool and
        type(affinity["exact"]) is bool,
        "affinity lifecycle changed")
    for key in ("captured", "observed"):
        item = affinity[key]
        require(item is None or (isinstance(item, list) and item and
                len(item) <= MAX_CPU_SET_SIZE and
                all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in item) and
                item == sorted(set(item))),
                f"affinity {key} CPU set is invalid")
    require((affinity["captured"] is not None) ==
            affinity["restore_attempted"],
            "affinity restore attempt does not match capture")
    if affinity["restore_succeeded"]:
        require(affinity["restore_attempted"] and
                affinity["observed"] is not None and
                affinity["error_type"] is None and affinity["error"] is None,
                "successful affinity restore lifecycle is inconsistent")
        require(affinity["exact"] is
                typed_equal(affinity["observed"], affinity["captured"]),
                "affinity exact-readback predicate changed")
    else:
        require(affinity["exact"] is False,
                "failed affinity restore cannot be exact")
        if affinity["restore_attempted"]:
            bounded_diagnostic(affinity["error_type"], "affinity restore type")
            bounded_diagnostic(affinity["error"], "affinity restore error")
        else:
            require(affinity["error_type"] is None and affinity["error"] is None and
                    affinity["observed"] is None,
                    "uncaptured affinity has spurious lifecycle data")
    if require_success:
        require(all(value[name]["entered"] and value[name]["exit"] == "pass"
                    for name in ("reservation", "pair_lease")) and
                all(check["status"] == "pass" for check in checks.values()) and
                affinity["restore_succeeded"] and affinity["exact"],
                "successful evidence has incomplete teardown or affinity restoration")
    return value


def run_guard_check(lifecycle: dict[str, Any], name: str,
                    operation: Any) -> None:
    check = lifecycle["guard_checks"][name]
    require(check["status"] == "not-reached",
            f"{name} guard check ran more than once")
    try:
        operation()
    except BaseException as error:
        error_type, message = _bounded_exception(error)
        check.update(status="failed", error_type=error_type, error=message)
        raise
    check["status"] = "pass"


class ExitRecordingGuard:
    """Ensure teardown runs, while preserving the body failure as primary."""

    def __init__(self, guard: Any, lifecycle: dict[str, Any], name: str):
        self.guard = guard
        self.lifecycle = lifecycle
        self.name = name

    def __enter__(self) -> Any:
        retained = self.guard.__enter__()
        self.lifecycle[self.name]["entered"] = True
        return retained

    def __exit__(self, exc_type: object, exc: object, tb: object) -> bool:
        try:
            self.guard.__exit__(exc_type, exc, tb)
            self.lifecycle[self.name]["exit"] = "pass"
        except BaseException as error:  # teardown must be recorded, never hidden
            error_type, message = _bounded_exception(error)
            self.lifecycle[self.name].update(
                exit="failed", error_type=error_type, error=message)
        return False


def typed_equal(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return set(left) == set(right) and all(
            typed_equal(left[key], right[key]) for key in left)
    if isinstance(left, (list, tuple)):
        return len(left) == len(right) and all(
            typed_equal(a, b) for a, b in zip(left, right))
    return bool(left == right)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def validate_utc(value: object, label: str) -> str:
    require(isinstance(value, str) and value.endswith("Z"),
            f"{label} is not canonical UTC")
    try:
        parsed = dt.datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise EvidenceError(f"{label} is not valid UTC") from error
    require(parsed.tzinfo is not None and parsed.utcoffset() == dt.timedelta(0),
            f"{label} is not UTC")
    return value


def canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"),
                          allow_nan=False).encode("utf-8")
    except (TypeError, ValueError, RecursionError) as error:
        raise EvidenceError(f"value is not canonical JSON: {error}") from error


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _lexical_absolute(path: Path) -> Path:
    return Path(os.path.abspath(os.fspath(path)))


def _open_parent_nofollow(path: Path) -> tuple[int, str, Path]:
    absolute = _lexical_absolute(path)
    parts = absolute.parts
    require(absolute.is_absolute() and len(parts) >= 2 and
            all(part not in ("", ".", "..") for part in parts[1:]) and
            len(os.fsencode(str(absolute))) <= MAX_ARTIFACT_PATH_BYTES,
            "file path is not a bounded canonical absolute path")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
        getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(parts[0], flags)
    try:
        for component in parts[1:-1]:
            child = os.open(component, flags, dir_fd=descriptor)
            metadata = os.fstat(child)
            require(stat.S_ISDIR(metadata.st_mode),
                    "file path component is not a directory")
            os.close(descriptor)
            descriptor = child
        return descriptor, parts[-1], absolute
    except Exception:
        os.close(descriptor)
        raise


def _snapshot_regular_fd(descriptor: int, maximum_bytes: int,
                         label: str) -> tuple[bytes, os.stat_result]:
    before = os.fstat(descriptor)
    require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1,
            f"{label} is not a single-link regular file")
    require(0 <= before.st_size <= maximum_bytes,
            f"{label} exceeds the evidence bound")
    os.lseek(descriptor, 0, os.SEEK_SET)
    chunks: list[bytes] = []
    retained = 0
    while retained <= maximum_bytes:
        chunk = os.read(descriptor, min(1 << 20, maximum_bytes + 1 - retained))
        if not chunk:
            break
        chunks.append(chunk)
        retained += len(chunk)
    require(retained <= maximum_bytes, f"{label} exceeds the evidence bound")
    data = b"".join(chunks)
    after = os.fstat(descriptor)
    stable_fields = ("st_dev", "st_ino", "st_mode", "st_nlink", "st_size",
                     "st_mtime_ns", "st_ctime_ns")
    require(all(getattr(before, key) == getattr(after, key)
                for key in stable_fields) and len(data) == before.st_size,
            f"{label} changed while it was read")
    return data, after


def regular_path_snapshot(path: Path, maximum_bytes: int,
                          label: str) -> tuple[bytes, os.stat_result, Path]:
    try:
        parent, name, absolute = _open_parent_nofollow(path)
    except OSError as error:
        raise EvidenceError(f"cannot securely open {label}: {error}") from error
    descriptor: int | None = None
    try:
        descriptor = os.open(
            name, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0),
            dir_fd=parent)
        data, metadata = _snapshot_regular_fd(descriptor, maximum_bytes, label)
        current = os.stat(name, dir_fd=parent, follow_symlinks=False)
        require((current.st_dev, current.st_ino, current.st_mode,
                 current.st_nlink, current.st_size, current.st_mtime_ns,
                 current.st_ctime_ns) ==
                (metadata.st_dev, metadata.st_ino, metadata.st_mode,
                 metadata.st_nlink, metadata.st_size, metadata.st_mtime_ns,
                 metadata.st_ctime_ns), f"{label} path was replaced while read")
        live_parent, live_name, _live_absolute = _open_parent_nofollow(absolute)
        try:
            live_directory = os.fstat(live_parent)
            opened_directory = os.fstat(parent)
            live = os.stat(
                live_name, dir_fd=live_parent, follow_symlinks=False)
            require((live_directory.st_dev, live_directory.st_ino) ==
                    (opened_directory.st_dev, opened_directory.st_ino) and
                    (live.st_dev, live.st_ino) ==
                    (metadata.st_dev, metadata.st_ino),
                    f"{label} path chain was replaced while read")
        finally:
            os.close(live_parent)
        return data, metadata, absolute
    except OSError as error:
        raise EvidenceError(f"cannot securely read {label}: {error}") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
        os.close(parent)


def sha256_file(path: Path) -> str:
    data, _metadata, _absolute = regular_path_snapshot(
        path, MAX_ARTIFACT_BYTES, path.name)
    return sha256_bytes(data)


def signed(value: Mapping[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(dict(value))
    require("digest" not in result, "digest field already exists")
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def verify_signature(value: object, label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    payload = copy.deepcopy(value)
    digest = payload.pop("digest", None)
    require(isinstance(digest, str) and SHA256_RE.fullmatch(digest) is not None,
            f"{label} has no canonical digest")
    require(sha256_bytes(canonical_bytes(payload)) == digest,
            f"{label} digest mismatch")
    return value


def strict_json(data: bytes, label: str) -> Any:
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in items:
            if key in result:
                raise EvidenceError(f"{label} contains duplicate key {key!r}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise EvidenceError(f"{label} contains non-finite number {value}")

    def finite_float(value: str) -> float:
        if len(value) > 128:
            raise EvidenceError(f"{label} contains an oversized number")
        parsed = float(value)
        if not math.isfinite(parsed):
            raise EvidenceError(f"{label} contains non-finite number {value}")
        return parsed

    def bounded_int(value: str) -> int:
        if len(value) > 128:
            raise EvidenceError(f"{label} contains an oversized integer")
        return int(value)

    try:
        return json.loads(data.decode("utf-8"), object_pairs_hook=pairs,
                          parse_constant=reject_constant, parse_float=finite_float,
                          parse_int=bounded_int)
    except (UnicodeDecodeError, ValueError, RecursionError) as error:
        raise EvidenceError(f"invalid {label}: {error}") from error


def canonical_pretty_bytes(value: object) -> bytes:
    try:
        return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) +
                "\n").encode("utf-8")
    except (TypeError, ValueError, RecursionError) as error:
        raise EvidenceError(f"value is not canonical pretty JSON: {error}") from error


def read_bounded(path: Path, maximum_bytes: int = MAX_ARTIFACT_BYTES) -> bytes:
    data, _metadata, _absolute = regular_path_snapshot(
        path, maximum_bytes, path.name)
    return data


def write_exclusive(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    parent: int | None = None
    descriptor: int | None = None
    created_identity: tuple[int, int] | None = None
    completed = False
    name = path.name
    try:
        parent, name, _absolute = _open_parent_nofollow(path)
        descriptor = os.open(
            name, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
            getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
            0o600, dir_fd=parent)
        opened = os.fstat(descriptor)
        created_identity = (opened.st_dev, opened.st_ino)
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            require(written > 0, f"short write to {path}")
            offset += written
        os.fsync(descriptor)
        metadata = os.fstat(descriptor)
        current = os.stat(name, dir_fd=parent, follow_symlinks=False)
        require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
                (metadata.st_dev, metadata.st_ino, metadata.st_size) ==
                (current.st_dev, current.st_ino, current.st_size) and
                metadata.st_size == len(data), f"written path changed: {path}")
        live_parent, live_name, _live_absolute = _open_parent_nofollow(path)
        try:
            live_directory = os.fstat(live_parent)
            opened_directory = os.fstat(parent)
            live = os.stat(
                live_name, dir_fd=live_parent, follow_symlinks=False)
            require((live_directory.st_dev, live_directory.st_ino) ==
                    (opened_directory.st_dev, opened_directory.st_ino) and
                    (live.st_dev, live.st_ino) ==
                    (metadata.st_dev, metadata.st_ino),
                    f"written path chain changed: {path}")
        finally:
            os.close(live_parent)
        completed = True
    except FileExistsError as error:
        raise EvidenceError(f"refusing to replace {path}") from error
    except OSError as error:
        raise EvidenceError(f"cannot securely write {path}: {error}") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if not completed and parent is not None and created_identity is not None:
            try:
                current = os.stat(name, dir_fd=parent, follow_symlinks=False)
                if (current.st_dev, current.st_ino) == created_identity:
                    os.unlink(name, dir_fd=parent)
            except OSError:
                pass
        if parent is not None:
            os.close(parent)


def write_json_exclusive(path: Path, value: object) -> None:
    write_exclusive(path, canonical_bytes(value) + b"\n")


def run_checked(arguments: Sequence[str], cwd: Path | None = None) -> bytes:
    completed, stdout, stderr, timed_out = run_capture_bounded(
        arguments, cwd=cwd, timeout_seconds=30,
        environment=CHILD_ENVIRONMENT)
    require(not timed_out, f"command timed out: {' '.join(arguments)}")
    if completed.returncode:
        raise EvidenceError(
            f"command failed ({completed.returncode}): {' '.join(arguments)}: "
            + stderr.decode("utf-8", errors="replace").strip())
    return stdout


def file_identity(path: Path, kind: str) -> dict[str, Any]:
    data, info, absolute = regular_path_snapshot(path, MAX_ARTIFACT_BYTES, kind)
    require(info.st_size > 0, f"{kind} size is outside the evidence bound")
    mode = stat.S_IMODE(info.st_mode)
    if kind in ("runner", "executable", "taskset", "python"):
        require(mode & 0o111 != 0, f"{kind} is not executable")
    if kind == "archive":
        require(data.startswith(b"!<arch>\n"), "C7 archive is not ar format")
    return {"kind": kind, "name": absolute.name, "size": info.st_size,
            "mode": mode, "device": info.st_dev, "inode": info.st_ino,
            "sha256": sha256_bytes(data)}


def committed_file_identity(root: Path, commit: str, relative: Path,
                            kind: str) -> dict[str, Any]:
    root = _lexical_absolute(root)
    path = _lexical_absolute(root / relative)
    try:
        path.relative_to(root)
    except ValueError as error:
        raise EvidenceError(f"{kind} escapes the source checkout") from error
    record = file_identity(path, kind)
    committed = run_checked(("git", "-C", str(root), "show",
                             f"{commit}:{relative.as_posix()}"))
    require(record["sha256"] == sha256_bytes(committed),
            f"{kind} differs from committed source")
    tree_entry = run_checked(("git", "-C", str(root), "ls-tree", "-z", commit,
                              "--", relative.as_posix()))
    require(tree_entry.endswith(b"\0") and tree_entry.count(b"\0") == 1 and
            b"\t" in tree_entry, f"{kind} committed mode is absent")
    header, retained_path = tree_entry[:-1].split(b"\t", 1)
    fields = header.split()
    require(len(fields) == 3 and fields[1] == b"blob" and
            fields[0] in (b"100644", b"100755") and
            retained_path == relative.as_posix().encode("utf-8"),
            f"{kind} committed mode is invalid")
    require(bool(record["mode"] & 0o111) == (fields[0] == b"100755"),
            f"{kind} executable mode differs from committed source")
    record["path"] = relative.as_posix()
    return record


def git_identity(root: Path, tooling_commit: str,
                 core_commit: str) -> dict[str, Any]:
    root = root.resolve(strict=True)
    require(GIT_SHA_RE.fullmatch(tooling_commit) is not None and
            GIT_SHA_RE.fullmatch(core_commit) is not None,
            "expected tooling/core commit is not full SHA-1")
    head = run_checked(("git", "-C", str(root), "rev-parse", "HEAD")).decode().strip()
    require(head == tooling_commit,
            f"source HEAD {head} differs from tooling commit {tooling_commit}")
    status_text = run_checked(("git", "-C", str(root), "status", "--porcelain=v1",
                               "--untracked-files=all")).decode()
    require(status_text == "", f"source checkout is not clean: {status_text!r}")
    ancestor = subprocess.run(
        ["git", "-C", str(root), "merge-base", "--is-ancestor",
         core_commit, tooling_commit], stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL, check=False, timeout=30,
        env=CHILD_ENVIRONMENT)
    require(ancestor.returncode == 0,
            "core commit is not an ancestor of tooling commit")
    tooling_tree = run_checked(("git", "-C", str(root), "rev-parse",
                                f"{tooling_commit}^{{tree}}")).decode().strip()
    core_tree = run_checked(("git", "-C", str(root), "rev-parse",
                             f"{core_commit}^{{tree}}")).decode().strip()
    listing = run_checked(("git", "-C", str(root), "ls-tree", "-r", "-z", "HEAD"))
    return {"tooling_commit": head, "tooling_tree": tooling_tree,
            "core_commit": core_commit, "core_tree": core_tree,
            "core_is_ancestor": True, "clean": True,
            "tracked_tree_sha256": sha256_bytes(listing)}


def manifest_artifact(value: object, expected_path: Path,
                      label: str) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {"bytes", "path", "sha256"},
            f"{label} build artifact schema changed")
    require(type(value["bytes"]) is int and 0 < value["bytes"] <= MAX_ARTIFACT_BYTES and
            value["path"] == expected_path.as_posix() and
            isinstance(value["sha256"], str) and
            SHA256_RE.fullmatch(value["sha256"]) is not None,
            f"{label} build artifact identity changed")
    return value


def artifact_matches_identity(record: Mapping[str, Any],
                              identity: Mapping[str, Any], label: str) -> None:
    size = identity.get("size", identity.get("bytes"))
    expected_path = identity.get("path")
    expected_name = identity.get("name")
    path_matches = (record["path"] == expected_path if expected_path is not None
                    else Path(record["path"]).name == expected_name)
    require(record["bytes"] == size and
            record["sha256"] == identity.get("sha256") and path_matches,
            f"{label} differs from the supplied exact file")


def retained_artifact(identity: Mapping[str, Any], path: Path) -> dict[str, Any]:
    return {"bytes": identity["size"], "path": path.as_posix(),
            "sha256": identity["sha256"]}


def retained_program(identity: Mapping[str, Any], path: Path) -> dict[str, Any]:
    return {"path": str(path.resolve(strict=True)),
            "sha256": identity["sha256"]}


def resolve_build_artifact(root: Path, record: Mapping[str, Any], label: str) -> Path:
    relative = record["path"]
    require(isinstance(relative, str) and relative and "\\" not in relative and
            ":" not in relative and not os.path.isabs(relative),
            f"{label} build path is not canonical")
    pure = Path(relative)
    require(pure.as_posix() == relative and
            all(part not in ("", ".", "..") for part in pure.parts),
            f"{label} build path is not canonical")
    candidate = root.joinpath(*pure.parts).resolve(strict=True)
    try:
        candidate.relative_to(root.resolve())
    except ValueError as error:
        raise EvidenceError(f"{label} build path escapes the checkout") from error
    return candidate


def manifest_program(value: object, label: str) -> dict[str, str]:
    require(isinstance(value, dict) and set(value) == {
        "path", "sha256", "version"} and
        isinstance(value["path"], str) and value["path"] and
        os.path.isabs(value["path"]) and
        isinstance(value["sha256"], str) and
        SHA256_RE.fullmatch(value["sha256"]) is not None and
        isinstance(value["version"], str) and value["version"],
        f"{label} program record changed")
    return {"path": value["path"], "sha256": value["sha256"]}


def derive_build_attestation(
    value: object, manifest_identity: Mapping[str, Any],
    expected_tooling_commit: str, expected_core_commit: str,
    source_identity: Mapping[str, Any], build_runner_identity: Mapping[str, Any],
    build_validator_identity: Mapping[str, Any], archive_identity: Mapping[str, Any],
    executable_identity: Mapping[str, Any],
    core_source_closure: Sequence[Mapping[str, Any]],
    taskset_identity: Mapping[str, Any], python_identity: Mapping[str, Any], *,
    source_root: Path | None = None,
    archive_path: Path | None = None, executable_path: Path | None = None,
) -> dict[str, Any]:
    require(isinstance(value, dict), "C7 build manifest is not an object")
    require(value.get("schema") == BUILD_MANIFEST_SCHEMA and
            value.get("status") == "pass" and value.get("scope") == BUILD_SCOPE and
            value.get("tooling_git_sha") == expected_tooling_commit and
            value.get("core_git_sha") == expected_core_commit,
            "C7 build manifest version, status, scope, or Git identity changed")
    source = manifest_artifact(value.get("source"), SOURCE_RELATIVE, "C7 source")
    build_runner = manifest_artifact(
        value.get("runner"), BUILD_RUNNER_RELATIVE, "C7 build runner")
    build_validator = manifest_artifact(
        value.get("validator"), BUILD_VALIDATOR_RELATIVE, "C7 build validator")
    artifact_matches_identity(source, source_identity, "C7 source")
    artifact_matches_identity(build_runner, build_runner_identity, "C7 build runner")
    artifact_matches_identity(
        build_validator, build_validator_identity, "C7 build validator")

    builds = value.get("builds")
    require(isinstance(builds, list) and
            [item.get("name") if isinstance(item, dict) else None for item in builds] ==
            list(BUILD_NAMES), "C7 build matrix changed")
    require(isinstance(core_source_closure, list) and
            len(core_source_closure) == len(CORE_SOURCE_CLOSURE),
            "retained core source closure changed")
    build_record_sha256: dict[str, str] = {}
    for build, name in zip(builds, BUILD_NAMES):
        require(typed_equal(build.get("source_closure"), core_source_closure),
                f"C7 {name} build source closure differs")
        launcher_python = manifest_program(
            build.get("launcher_python"), f"C7 {name} launcher Python")
        require(typed_equal(launcher_python, python_identity),
                f"C7 {name} launcher Python differs from the runner")
        build_record_sha256[name] = sha256_bytes(canonical_bytes(build))
    taskset = manifest_program(value.get("taskset"), "C7 taskset")
    require(typed_equal(taskset, taskset_identity),
            "C7 taskset differs from the runner")
    avx2 = builds[2]
    require(avx2.get("name") == "avx2" and avx2.get("backend") == "avx2" and
            avx2.get("sanitizer") is False,
            "C7 authoritative build is not the non-sanitized AVX2 build")
    build_dir = avx2.get("build_dir")
    require(isinstance(build_dir, str) and
            build_dir.endswith("/core-avx2") and "\\" not in build_dir and
            ":" not in build_dir and not os.path.isabs(build_dir) and
            Path(build_dir).as_posix() == build_dir and
            all(part not in ("", ".", "..") for part in Path(build_dir).parts),
            "C7 AVX2 build directory changed")
    library = manifest_artifact(
        avx2.get("library"), Path(build_dir) / "liblibleopard.a", "C7 AVX2 archive")
    executable_record = manifest_artifact(
        avx2.get("executable"), Path(build_dir).parent / "c7-avx2",
        "C7 AVX2 executable")
    artifact_matches_identity(library, archive_identity, "C7 AVX2 archive")
    artifact_matches_identity(
        executable_record, executable_identity, "C7 AVX2 executable")
    if source_root is not None:
        require(archive_path is not None and executable_path is not None,
                "local C7 build paths are absent")
        require(resolve_build_artifact(source_root, library, "C7 AVX2 archive") ==
                archive_path.resolve(strict=True) and
                resolve_build_artifact(source_root, executable_record,
                                       "C7 AVX2 executable") ==
                executable_path.resolve(strict=True),
                "supplied C7 AVX2 paths differ from the verified build manifest")

    compile_argv = avx2.get("compile_argv")
    require(isinstance(compile_argv, list) and
            all(isinstance(item, str) for item in compile_argv) and
            compile_argv.count("-O2") == 1 and compile_argv.count("-fopenmp") == 1 and
            not any(re.fullmatch(r"-O(?:0|1|3|g|fast|s|z)", item)
                    for item in compile_argv) and
            not any(item.startswith("-fsanitize") or
                    item == "-DLEO2_C7_REQUIRE_ASAN_UBSAN=1"
                    for item in compile_argv) and
            source["path"] in compile_argv and library["path"] in compile_argv and
            "-o" in compile_argv and
            compile_argv.index("-o") + 1 < len(compile_argv) and
            compile_argv[compile_argv.index("-o") + 1] == executable_record["path"],
            "C7 AVX2 exact compile command is not optimized and non-sanitized")
    instrumentation = avx2.get("instrumentation")
    zero_counts = {"asan_lines": 0, "ubsan_lines": 0}
    require(isinstance(instrumentation, dict) and
            instrumentation.get("required_compile_macro") is False and
            typed_equal(instrumentation.get("executable_counts"), zero_counts) and
            typed_equal(instrumentation.get("core_archive_counts"), zero_counts),
            "C7 AVX2 sanitizer attribution changed")
    reproducibility = value.get("reproducibility")
    require(isinstance(reproducibility, dict) and
            isinstance(reproducibility.get("comparison"), dict) and
            reproducibility["comparison"].get("status") == "pass" and
            isinstance(reproducibility.get("fingerprints"), dict) and
            typed_equal(reproducibility["fingerprints"].get("avx2"), {
                "library_sha256": library["sha256"],
                "executable_sha256": executable_record["sha256"],
            }), "C7 v4 A/B reproducibility attestation is absent or changed")
    manifest_size = manifest_identity.get("size", manifest_identity.get("bytes"))
    require(type(manifest_size) is int and manifest_size > 0 and
            isinstance(manifest_identity.get("sha256"), str) and
            SHA256_RE.fullmatch(manifest_identity["sha256"]) is not None,
            "C7 build manifest file identity changed")
    return {
        "schema": BUILD_ATTESTATION_SCHEMA,
        "manifest": {"bytes": manifest_size,
                     "schema": BUILD_MANIFEST_SCHEMA,
                     "sha256": manifest_identity["sha256"]},
        "status": "pass", "scope": BUILD_SCOPE,
        "tooling_commit": expected_tooling_commit,
        "core_commit": expected_core_commit,
        "comparison_status": "pass",
        "source": copy.deepcopy(source),
        "build_runner": copy.deepcopy(build_runner),
        "build_validator": copy.deepcopy(build_validator),
        "taskset": copy.deepcopy(taskset),
        "python": copy.deepcopy(python_identity),
        "core_source_closure_sha256": sha256_bytes(
            canonical_bytes(core_source_closure)),
        "build_record_sha256": build_record_sha256,
        "avx2": {
            "name": "avx2", "backend": "avx2", "sanitizer": False,
            "library": copy.deepcopy(library),
            "executable": copy.deepcopy(executable_record),
            "compile_argv_sha256": sha256_bytes(canonical_bytes(compile_argv)),
            "optimization": "-O2", "openmp": "-fopenmp",
            "instrumentation": {
                "required_compile_macro": False,
                "executable_counts": copy.deepcopy(zero_counts),
                "core_archive_counts": copy.deepcopy(zero_counts),
            },
        },
    }


def run_build_validator(source_root: Path, build_manifest: Path) -> None:
    validator = _lexical_absolute(source_root / BUILD_VALIDATOR_RELATIVE)
    completed, stdout, stderr, timed_out = run_capture_bounded(
        [sys.executable, str(validator), str(build_manifest),
         "--source-root", str(source_root), "--evidence-root", str(source_root),
         "--live", "--require-checkout-head"],
        cwd=None, timeout_seconds=300, environment=VALIDATOR_ENVIRONMENT)
    require(not timed_out and completed.returncode == 0,
            "C7 v4 live build validation failed: " +
            stderr.decode("utf-8", errors="replace").strip())
    require(stdout == b"C7 evidence validation passed (live)\n" and
            stderr == b"", "C7 v4 live validator output changed")


def verified_build_attestation(
    source_root: Path, build_manifest: Path, expected_tooling_commit: str,
    expected_core_commit: str, source_identity: Mapping[str, Any],
    build_runner_identity: Mapping[str, Any],
    build_validator_identity: Mapping[str, Any], archive_identity: Mapping[str, Any],
    executable_identity: Mapping[str, Any],
    core_source_closure: Sequence[Mapping[str, Any]],
    taskset_identity: Mapping[str, Any], python_identity: Mapping[str, Any],
    archive_path: Path, executable_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    build_manifest = build_manifest.resolve(strict=True)
    try:
        build_manifest.relative_to(source_root.resolve())
    except ValueError as error:
        raise EvidenceError("C7 v4 build manifest is outside the source checkout") from error
    before = file_identity(build_manifest, "build-manifest")
    run_build_validator(source_root, build_manifest)
    after = file_identity(build_manifest, "build-manifest")
    require(typed_equal(before, after), "C7 build manifest changed during live validation")
    data = read_bounded(build_manifest)
    require(len(data) == before["size"] and sha256_bytes(data) == before["sha256"],
            "C7 build manifest changed after live validation")
    parsed = strict_json(data, "C7 v4 build manifest")
    require(data == canonical_pretty_bytes(parsed),
            "C7 v4 build manifest is not canonical pretty JSON")
    attestation = derive_build_attestation(
        parsed, before, expected_tooling_commit, expected_core_commit,
        source_identity, build_runner_identity, build_validator_identity,
        archive_identity, executable_identity, core_source_closure,
        taskset_identity, python_identity, source_root=source_root,
        archive_path=archive_path, executable_path=executable_path)
    return before, attestation


def input_snapshot(source_root: Path, expected_tooling_commit: str,
                   expected_core_commit: str, archive: Path,
                   executable: Path, taskset: Path,
                   build_manifest: Path) -> dict[str, Any]:
    source_root = source_root.resolve(strict=True)
    runner = Path(__file__).resolve(strict=True)
    require(runner == (source_root / RUNNER_RELATIVE).resolve(strict=True),
            "runner is not the committed runner in --source-root")
    core_source_identities = [
        committed_file_identity(
            source_root, expected_core_commit, Path(relative), "core-source")
        for relative in CORE_SOURCE_CLOSURE]
    core_source_closure = [
        retained_artifact(identity, Path(relative))
        for identity, relative in zip(core_source_identities, CORE_SOURCE_CLOSURE)]
    live_git = git_identity(
        source_root, expected_tooling_commit, expected_core_commit)
    git = {"tooling_commit": live_git["tooling_commit"],
           "core_commit": live_git["core_commit"]}
    runner_identity = committed_file_identity(
        source_root, expected_tooling_commit, RUNNER_RELATIVE, "runner")
    source_identity = committed_file_identity(
        source_root, expected_tooling_commit, SOURCE_RELATIVE, "source")
    build_runner_identity = committed_file_identity(
        source_root, expected_tooling_commit, BUILD_RUNNER_RELATIVE, "build-runner")
    build_validator_identity = committed_file_identity(
        source_root, expected_tooling_commit, BUILD_VALIDATOR_RELATIVE,
        "build-validator")
    archive_identity = file_identity(archive, "archive")
    executable_identity = file_identity(executable, "executable")
    taskset_identity = retained_program(
        file_identity(taskset, "taskset"), taskset)
    python_path = Path(sys.executable).resolve(strict=True)
    python_identity = retained_program(
        file_identity(python_path, "python"), python_path)
    build_manifest_identity, build_attestation = verified_build_attestation(
        source_root, build_manifest, expected_tooling_commit, expected_core_commit,
        source_identity, build_runner_identity, build_validator_identity,
        archive_identity, executable_identity, core_source_closure,
        taskset_identity, python_identity, archive, executable)
    result = {
        "schema": INPUT_SCHEMA, "git": git,
        "runner": retained_artifact(runner_identity, RUNNER_RELATIVE),
        "source": retained_artifact(source_identity, SOURCE_RELATIVE),
        "build_runner": retained_artifact(
            build_runner_identity, BUILD_RUNNER_RELATIVE),
        "build_validator": retained_artifact(
            build_validator_identity, BUILD_VALIDATOR_RELATIVE),
        "build_manifest": {
            "bytes": build_manifest_identity["size"],
            "schema": BUILD_MANIFEST_SCHEMA,
            "sha256": build_manifest_identity["sha256"],
        },
        "build_attestation": build_attestation,
        "archive": copy.deepcopy(build_attestation["avx2"]["library"]),
        "executable": copy.deepcopy(build_attestation["avx2"]["executable"]),
        "taskset": taskset_identity, "python": python_identity,
        "core_source_closure": core_source_closure,
    }
    result["binding_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def validate_input_snapshot(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "git", "runner", "source", "build_runner", "build_validator",
        "build_manifest", "build_attestation", "archive", "executable",
        "taskset", "python", "core_source_closure", "binding_sha256"},
        "input snapshot schema changed")
    require(value["schema"] == INPUT_SCHEMA, "input snapshot version changed")
    payload = {key: value[key] for key in value if key != "binding_sha256"}
    require(isinstance(value["binding_sha256"], str) and
            SHA256_RE.fullmatch(value["binding_sha256"]) is not None and
            value["binding_sha256"] == sha256_bytes(canonical_bytes(payload)),
            "input snapshot binding changed")
    git = value["git"]
    require(isinstance(git, dict) and set(git) == {
        "tooling_commit", "core_commit"} and
        all(isinstance(git[key], str) and GIT_SHA_RE.fullmatch(git[key]) is not None
            for key in ("tooling_commit", "core_commit")),
        "retained Git identity is invalid")
    artifact_paths = {
        "runner": RUNNER_RELATIVE, "source": SOURCE_RELATIVE,
        "build_runner": BUILD_RUNNER_RELATIVE,
        "build_validator": BUILD_VALIDATOR_RELATIVE,
    }
    for key, expected_path in artifact_paths.items():
        manifest_artifact(value[key], expected_path, f"retained {key}")
    for key in ("archive", "executable"):
        record = value[key]
        require(isinstance(record, dict) and set(record) == {
            "bytes", "path", "sha256"} and
            type(record["bytes"]) is int and 0 < record["bytes"] <=
                MAX_ARTIFACT_BYTES and
            isinstance(record["path"], str) and record["path"] and
            not os.path.isabs(record["path"]) and "\\" not in record["path"] and
            ":" not in record["path"] and
            Path(record["path"]).as_posix() == record["path"] and
            all(part not in ("", ".", "..") for part in Path(record["path"]).parts) and
            isinstance(record["sha256"], str) and
            SHA256_RE.fullmatch(record["sha256"]) is not None,
            f"retained {key} identity is invalid")
    build_manifest = value["build_manifest"]
    require(isinstance(build_manifest, dict) and set(build_manifest) == {
        "bytes", "schema", "sha256"} and
        type(build_manifest["bytes"]) is int and
        0 < build_manifest["bytes"] <= MAX_ARTIFACT_BYTES and
        build_manifest["schema"] == BUILD_MANIFEST_SCHEMA and
        isinstance(build_manifest["sha256"], str) and
        SHA256_RE.fullmatch(build_manifest["sha256"]) is not None,
        "retained build manifest identity is invalid")
    for key in ("taskset", "python"):
        program = value[key]
        require(isinstance(program, dict) and set(program) == {"path", "sha256"} and
                isinstance(program["path"], str) and os.path.isabs(program["path"]) and
                isinstance(program["sha256"], str) and
                SHA256_RE.fullmatch(program["sha256"]) is not None,
                f"retained {key} identity is invalid")
    attestation = value["build_attestation"]
    require(isinstance(attestation, dict) and set(attestation) == {
        "schema", "manifest", "status", "scope",
        "tooling_commit", "core_commit", "comparison_status", "source",
        "build_runner", "build_validator", "taskset", "python",
        "core_source_closure_sha256", "build_record_sha256", "avx2"} and
        attestation["schema"] == BUILD_ATTESTATION_SCHEMA and
        attestation["status"] == "pass" and attestation["scope"] == BUILD_SCOPE and
        attestation["tooling_commit"] == git["tooling_commit"] and
        attestation["core_commit"] == git["core_commit"] and
        attestation["comparison_status"] == "pass",
        "retained C7 v4 build attestation changed")
    manifest = attestation["manifest"]
    require(isinstance(manifest, dict) and set(manifest) == {
        "bytes", "schema", "sha256"} and
        manifest["schema"] == BUILD_MANIFEST_SCHEMA and
        manifest["bytes"] == value["build_manifest"]["bytes"] and
        manifest["sha256"] == value["build_manifest"]["sha256"],
        "retained C7 build manifest attestation changed")
    require(typed_equal(attestation["source"], {
        "bytes": value["source"]["bytes"], "path": SOURCE_RELATIVE.as_posix(),
        "sha256": value["source"]["sha256"]}) and
        typed_equal(attestation["build_runner"], {
            "bytes": value["build_runner"]["bytes"],
            "path": BUILD_RUNNER_RELATIVE.as_posix(),
            "sha256": value["build_runner"]["sha256"]}) and
        typed_equal(attestation["build_validator"], {
            "bytes": value["build_validator"]["bytes"],
            "path": BUILD_VALIDATOR_RELATIVE.as_posix(),
            "sha256": value["build_validator"]["sha256"]}) and
        typed_equal(attestation["taskset"], value["taskset"]) and
        typed_equal(attestation["python"], value["python"]),
        "retained build tooling attestation changed")
    require(isinstance(attestation["core_source_closure_sha256"], str) and
            attestation["core_source_closure_sha256"] ==
                sha256_bytes(canonical_bytes(value["core_source_closure"])) and
            isinstance(attestation["build_record_sha256"], dict) and
            set(attestation["build_record_sha256"]) == set(BUILD_NAMES) and
            all(isinstance(item, str) and SHA256_RE.fullmatch(item) is not None
                for item in attestation["build_record_sha256"].values()),
            "retained build closure attestation changed")
    avx2 = attestation["avx2"]
    require(isinstance(avx2, dict) and set(avx2) == {
        "name", "backend", "sanitizer", "library", "executable",
        "compile_argv_sha256", "optimization", "openmp", "instrumentation"} and
        avx2["name"] == "avx2" and avx2["backend"] == "avx2" and
        avx2["sanitizer"] is False and avx2["optimization"] == "-O2" and
        avx2["openmp"] == "-fopenmp" and
        isinstance(avx2["compile_argv_sha256"], str) and
        SHA256_RE.fullmatch(avx2["compile_argv_sha256"]) is not None and
        isinstance(avx2["library"], dict) and
        set(avx2["library"]) == {"bytes", "path", "sha256"} and
        typed_equal(avx2["library"], value["archive"]) and
        isinstance(avx2["executable"], dict) and
        set(avx2["executable"]) == {"bytes", "path", "sha256"} and
        typed_equal(avx2["executable"], value["executable"]) and
        typed_equal(avx2["instrumentation"], {
            "required_compile_macro": False,
            "executable_counts": {"asan_lines": 0, "ubsan_lines": 0},
            "core_archive_counts": {"asan_lines": 0, "ubsan_lines": 0},
        }), "retained non-sanitized AVX2 build attestation changed")
    closure = value["core_source_closure"]
    require(isinstance(closure, list) and len(closure) == len(CORE_SOURCE_CLOSURE),
            "core source closure length changed")
    for record, expected in zip(closure, CORE_SOURCE_CLOSURE):
        manifest_artifact(record, Path(expected), f"core source {expected}")
    return value


def parse_cpu_list(text: str) -> set[int]:
    require(isinstance(text, str) and len(text) <= MAX_CPU_SET_SIZE * 16 and
            re.fullmatch(r"[0-9]+(?:-[0-9]+)?(?:,[0-9]+(?:-[0-9]+)?)*",
                         text.strip()) is not None,
            "invalid CPU list syntax")
    result: set[int] = set()
    for item in text.strip().split(","):
        if "-" in item:
            bounds = item.split("-", 1)
            require(len(bounds) == 2 and all(
                re.fullmatch(r"[0-9]+", bound) is not None and
                len(bound) <= len(str(MAX_CPU_ID)) for bound in bounds),
                "invalid CPU list")
            first, last = map(int, bounds)
            require(0 <= first <= last <= MAX_CPU_ID and
                    len(result) + last - first + 1 <= MAX_CPU_SET_SIZE,
                    "CPU list exceeds the supported bound")
            result.update(range(first, last + 1))
        elif item:
            require(re.fullmatch(r"[0-9]+", item) is not None and
                    len(item) <= len(str(MAX_CPU_ID)),
                    "invalid CPU list")
            cpu = int(item)
            require(cpu <= MAX_CPU_ID and len(result) < MAX_CPU_SET_SIZE,
                    "CPU list exceeds the supported bound")
            result.add(cpu)
    return result


def read_optional(path: Path) -> str | None:
    try:
        return path.read_text(encoding="ascii").strip()
    except FileNotFoundError:
        return None


def cpu_record(cpu: int) -> dict[str, Any]:
    root = Path(f"/sys/devices/system/cpu/cpu{cpu}")
    topology = {name: read_optional(root / "topology" / name) for name in
                ("core_id", "physical_package_id", "die_id",
                 "thread_siblings_list", "core_siblings_list")}
    frequency = {name: read_optional(root / "cpufreq" / name) for name in
                 ("scaling_driver", "scaling_governor",
                  "energy_performance_preference", "scaling_min_freq",
                  "scaling_max_freq", "cpuinfo_min_freq", "cpuinfo_max_freq")}
    return {"cpu": cpu, "online": read_optional(root / "online"),
            "topology": topology, "frequency_policy": frequency}


def host_identity(cpu: int, sibling: int, allowed: set[int]) -> dict[str, Any]:
    uname = platform.uname()
    return {
        "system": {"system": uname.system, "release": uname.release,
                   "version": uname.version, "machine": uname.machine,
                   "python": platform.python_version()},
        "allowed_at_launch": sorted(allowed),
        "timing_cpu": cpu_record(cpu), "sibling_cpu": cpu_record(sibling),
        "turbo": {str(path): read_optional(path) for path in (
            Path("/sys/devices/system/cpu/intel_pstate/no_turbo"),
            Path("/sys/devices/system/cpu/amd_pstate/status"),
            Path("/sys/devices/system/cpu/cpufreq/boost"))},
    }


def validate_host(value: object, cpu: int, sibling: int) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "system", "allowed_at_launch", "timing_cpu", "sibling_cpu", "turbo"},
        "host identity schema changed")
    system = value["system"]
    require(isinstance(system, dict) and set(system) == {
        "system", "release", "version", "machine", "python"} and
        all(isinstance(item, str) and item for item in system.values()),
        "host system identity is invalid")
    allowed = value["allowed_at_launch"]
    require(isinstance(allowed, list) and
            len(allowed) <= MAX_CPU_SET_SIZE and
            all(type(item) is int and 0 <= item <= MAX_CPU_ID for item in allowed) and
            allowed == sorted(set(allowed)) and cpu in allowed and sibling in allowed,
            "host launch affinity is invalid")
    topology_keys = {
        "core_id", "physical_package_id", "die_id",
        "thread_siblings_list", "core_siblings_list",
    }
    frequency_keys = {
        "scaling_driver", "scaling_governor", "energy_performance_preference",
        "scaling_min_freq", "scaling_max_freq", "cpuinfo_min_freq",
        "cpuinfo_max_freq",
    }
    for key, expected_cpu in (("timing_cpu", cpu), ("sibling_cpu", sibling)):
        record = value[key]
        require(isinstance(record, dict) and set(record) == {
            "cpu", "online", "topology", "frequency_policy"} and
            type(record["cpu"]) is int and record["cpu"] == expected_cpu and
            (record["online"] is None or isinstance(record["online"], str)) and
            isinstance(record["topology"], dict) and
            set(record["topology"]) == topology_keys and
            isinstance(record["frequency_policy"], dict) and
            set(record["frequency_policy"]) == frequency_keys and
            all(item is None or isinstance(item, str)
                for item in (*record["topology"].values(),
                             *record["frequency_policy"].values())),
            f"host {key} identity is invalid")
        siblings = record["topology"]["thread_siblings_list"]
        require(isinstance(siblings, str) and
                parse_cpu_list(siblings) == {cpu, sibling},
                f"host {key} SMT topology changed")
    turbo_keys = {
        "/sys/devices/system/cpu/intel_pstate/no_turbo",
        "/sys/devices/system/cpu/amd_pstate/status",
        "/sys/devices/system/cpu/cpufreq/boost",
    }
    require(isinstance(value["turbo"], dict) and
            set(value["turbo"]) == turbo_keys and
            all(item is None or isinstance(item, str)
                for item in value["turbo"].values()),
            "host turbo identity is invalid")
    return value


def validate_topology(cpu: int, sibling: int) -> tuple[set[int], set[int]]:
    require(type(cpu) is int and type(sibling) is int and
            0 <= cpu <= MAX_CPU_ID and 0 <= sibling <= MAX_CPU_ID and
            cpu != sibling, "invalid CPU pair")
    allowed = set(os.sched_getaffinity(0))
    require(cpu in allowed and sibling in allowed,
            "timing CPU and sibling must both be in launch affinity")
    text = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list") \
        .read_text(encoding="ascii")
    require(parse_cpu_list(text) == {cpu, sibling}, "requested CPUs are not one SMT pair")
    housekeeping = allowed - {cpu, sibling}
    require(housekeeping, "no housekeeping CPU remains")
    return allowed, housekeeping


def cpu_stat_snapshot(cpu: int) -> dict[str, Any]:
    require(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID,
            "CPU stat identity is invalid")
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            tokens = line.split()
            require(tokens[0] == f"cpu{cpu}" and
                    len(tokens) >= 1 + len(CPU_STAT_FIELDS),
                    f"CPU {cpu} has an incomplete /proc/stat record")
            try:
                values = [int(item) for item in
                          tokens[1:1 + len(CPU_STAT_FIELDS)]]
            except ValueError as error:
                raise EvidenceError(
                    f"CPU {cpu} has a non-integer /proc/stat record") from error
            require(all(value >= 0 for value in values),
                    f"CPU {cpu} has a negative /proc/stat counter")
            fields = dict(zip(CPU_STAT_FIELDS, values))
            return {"cpu": cpu, "fields": fields,
                    "total_jiffies": sum(values)}
    raise EvidenceError(f"CPU {cpu} is absent from /proc/stat")


def validate_cpu_stat_snapshot(value: object, expected_cpu: int,
                               label: str) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "cpu", "fields", "total_jiffies"} and
        type(value["cpu"]) is int and value["cpu"] == expected_cpu and
        isinstance(value["fields"], dict) and
        set(value["fields"]) == set(CPU_STAT_FIELDS) and
        all(type(value["fields"][name]) is int and value["fields"][name] >= 0
            for name in CPU_STAT_FIELDS) and
        type(value["total_jiffies"]) is int and
        value["total_jiffies"] == sum(value["fields"].values()),
        f"{label} CPU stat snapshot is invalid")
    return value


def cpu_stat_delta(before: Mapping[str, Any],
                   after: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(before, dict) and isinstance(after, dict) and
            type(before.get("cpu")) is int and type(after.get("cpu")) is int and
            before.get("cpu") == after.get("cpu"),
            "CPU stat snapshots refer to different CPUs")
    expected_cpu = before["cpu"]
    before = validate_cpu_stat_snapshot(before, expected_cpu, "before")
    after = validate_cpu_stat_snapshot(after, expected_cpu, "after")
    before_fields, after_fields = before["fields"], after["fields"]
    deltas: dict[str, int] = {}
    for name in CPU_STAT_FIELDS:
        first, last = before_fields[name], after_fields[name]
        require(type(first) is int and first >= 0 and
                type(last) is int and last >= first,
                f"CPU stat counter {name} moved backwards")
        deltas[name] = last - first
    idle = sum(deltas[name] for name in CPU_STAT_IDLE_FIELDS)
    nonidle = sum(value for name, value in deltas.items()
                  if name not in CPU_STAT_IDLE_FIELDS)
    return {"cpu": before["cpu"], "fields": deltas,
            "idle_jiffies": idle, "nonidle_jiffies": nonidle,
            "total_jiffies": idle + nonidle}


def kernel_lease_namespace(kind: str, payload: Mapping[str, Any]) -> str:
    require(kind in ("pair", "reservation", "terminal"),
            "kernel lease kind is invalid")
    digest = sha256_bytes(canonical_bytes(payload))
    return f"@leopard2-{kind}-{digest}"


def validate_kernel_lease_identity(value: object, kind: str,
                                   payload: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "authority", "namespace", "device", "inode", "payload_sha256"} and
        value["schema"] == KERNEL_LEASE_SCHEMA and
        value["authority"] == "linux-abstract-unix-bind" and
        value["namespace"] == kernel_lease_namespace(kind, payload) and
        type(value["device"]) is int and value["device"] >= 0 and
        type(value["inode"]) is int and value["inode"] > 0 and
        value["payload_sha256"] == sha256_bytes(canonical_bytes(payload)),
        "kernel lease identity changed")
    return value


class KernelNamespaceLease:
    """An unlinked, kernel-owned exclusive name that path replacement cannot split."""

    def __init__(self, kind: str, payload: Mapping[str, Any]):
        self.kind = kind
        self.payload = copy.deepcopy(dict(payload))
        self.namespace = kernel_lease_namespace(kind, self.payload)
        self.socket: socket.socket | None = None
        self.identity: dict[str, Any] | None = None
        self._identity_bytes: bytes | None = None

    def acquire(self) -> dict[str, Any]:
        require(self.socket is None, "kernel lease is already acquired")
        retained = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        try:
            retained.bind("\0" + self.namespace[1:])
        except OSError as error:
            retained.close()
            raise EvidenceError(
                f"{self.kind} kernel lease is already held or unavailable") from error
        metadata = os.fstat(retained.fileno())
        require(stat.S_ISSOCK(metadata.st_mode) and metadata.st_nlink == 1 and
                retained.getsockname() ==
                    b"\0" + self.namespace[1:].encode("ascii"),
                "kernel lease socket identity is invalid")
        self.socket = retained
        self.identity = {
            "schema": KERNEL_LEASE_SCHEMA,
            "authority": "linux-abstract-unix-bind",
            "namespace": self.namespace, "device": metadata.st_dev,
            "inode": metadata.st_ino,
            "payload_sha256": sha256_bytes(canonical_bytes(self.payload)),
        }
        self._identity_bytes = canonical_bytes(self.identity)
        return copy.deepcopy(self.identity)

    def validate_current(self) -> None:
        require(self.socket is not None and self.identity is not None,
                "kernel lease is not held")
        require(self._identity_bytes == canonical_bytes(self.identity),
                "kernel lease immutable identity changed")
        metadata = os.fstat(self.socket.fileno())
        require(stat.S_ISSOCK(metadata.st_mode) and metadata.st_nlink == 1 and
                (metadata.st_dev, metadata.st_ino) ==
                    (self.identity["device"], self.identity["inode"]) and
                self.socket.getsockname() ==
                    b"\0" + self.namespace[1:].encode("ascii"),
                "kernel lease descriptor changed")
        validate_kernel_lease_identity(self.identity, self.kind, self.payload)

    def release(self) -> None:
        if self.socket is not None:
            self.socket.close()
            self.socket = None


def pair_lease_payload(cpu: int, sibling: int,
                       uid: int | None = None) -> dict[str, Any]:
    require(type(cpu) is int and type(sibling) is int and
            0 <= cpu <= MAX_CPU_ID and 0 <= sibling <= MAX_CPU_ID and
            cpu != sibling,
            "CPU pair lease CPUs are invalid")
    retained_uid = os.getuid() if uid is None else uid
    require(type(retained_uid) is int and retained_uid >= 0,
            "CPU pair lease UID is invalid")
    return {"cpus": sorted((cpu, sibling)), "schema": PAIR_LEASE_SCHEMA,
            "uid": retained_uid}


def pair_lease_name(cpu: int, sibling: int, uid: int | None = None) -> str:
    require(type(cpu) is int and type(sibling) is int and
            0 <= cpu <= MAX_CPU_ID and 0 <= sibling <= MAX_CPU_ID and
            cpu != sibling,
            "CPU pair lease CPUs are invalid")
    retained_uid = os.getuid() if uid is None else uid
    first, second = sorted((cpu, sibling))
    return f"leopard2-cpu-pair-{retained_uid}-{first}-{second}.lock"


def pair_lease_runtime_root(uid: int | None = None) -> Path:
    retained_uid = os.getuid() if uid is None else uid
    return Path("/run/user") / str(retained_uid)


def pair_lease_directory(uid: int | None = None) -> Path:
    return pair_lease_runtime_root(uid) / "leopard2-cpu-leases"


class StableLeaseAnchor:
    """Shared non-replaceable ancestor lock for all Leopard2 evidence runners."""

    def __init__(self, path: Path | None = None):
        self.path = pair_lease_runtime_root() if path is None else _lexical_absolute(path)
        self.descriptor: int | None = None
        self.identity: tuple[int, int, int, int] | None = None

    @staticmethod
    def _metadata(metadata: os.stat_result) -> tuple[int, int, int, int]:
        mode = stat.S_IMODE(metadata.st_mode)
        require(stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and mode & 0o022 == 0,
                "stable runner lease anchor has unsafe ownership or mode")
        return metadata.st_dev, metadata.st_ino, metadata.st_uid, mode

    def validate_current(self) -> None:
        require(self.descriptor is not None and self.identity is not None,
                "stable runner lease anchor is not held")
        descriptor = os.fstat(self.descriptor)
        path = os.lstat(self.path)
        require(self._metadata(descriptor) == self.identity and
                (descriptor.st_dev, descriptor.st_ino) ==
                    (path.st_dev, path.st_ino),
                "stable runner lease anchor path was replaced")

    def __enter__(self) -> StableLeaseAnchor:
        require(self.descriptor is None,
                "stable runner lease anchor is already held")
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
            getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
        try:
            self.descriptor = os.open(self.path, flags)
            fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            descriptor = os.fstat(self.descriptor)
            path = os.lstat(self.path)
            self.identity = self._metadata(descriptor)
            require((descriptor.st_dev, descriptor.st_ino) ==
                    (path.st_dev, path.st_ino),
                    "stable runner lease anchor changed during acquisition")
            self.validate_current()
            return self
        except BaseException as error:
            self.__exit__(None, None, None)
            if isinstance(error, OSError):
                raise EvidenceError(
                    f"cannot acquire stable runner lease anchor: {error}") from error
            raise

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, exc, tb
        try:
            if self.descriptor is not None:
                descriptor = self.descriptor
                self.descriptor = None
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(descriptor)
        finally:
            self.identity = None


class PairLease:
    """Serialize Leopard2 evidence runners by one physical CPU pair."""

    def __init__(self, cpu: int, sibling: int,
                 runtime_root: Path | None = None,
                 anchor: StableLeaseAnchor | None = None):
        require(type(cpu) is int and type(sibling) is int and
                0 <= cpu <= MAX_CPU_ID and 0 <= sibling <= MAX_CPU_ID and
                cpu != sibling,
                "pair lease requires two distinct non-negative CPUs")
        self.cpu = cpu
        self.sibling = sibling
        self.runtime_root = (pair_lease_runtime_root() if runtime_root is None
                             else runtime_root)
        self.root = self.runtime_root / "leopard2-cpu-leases"
        self.path = self.root / pair_lease_name(cpu, sibling)
        self.descriptor: int | None = None
        self.directory_fd: int | None = None
        self.identity: dict[str, Any] | None = None
        self._identity_bytes: bytes | None = None
        self.anchor = anchor
        self.owns_anchor = anchor is None
        if self.owns_anchor:
            self.anchor = StableLeaseAnchor(self.runtime_root)
        self.kernel = KernelNamespaceLease(
            "pair", pair_lease_payload(cpu, sibling))

    def _validate_runtime_root(self) -> os.stat_result:
        metadata = os.lstat(self.runtime_root)
        require(stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o700,
                "CPU lease runtime root is not an owned mode-0700 directory")
        return metadata

    def _validate_directory(self) -> os.stat_result:
        self._validate_runtime_root()
        metadata = os.lstat(self.root)
        require(stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o700,
                "CPU pair lease directory is not an owned mode-0700 directory")
        return metadata

    def validate_current(self, expected_identity: Mapping[str, Any] | None = None) -> None:
        require(self.descriptor is not None and self.directory_fd is not None and
                self.identity is not None and self.anchor is not None,
                "CPU pair lease is not held")
        self.anchor.validate_current()
        require(self._identity_bytes == canonical_bytes(self.identity),
                "CPU pair lease immutable identity changed")
        if expected_identity is not None:
            require(typed_equal(expected_identity, self.identity),
                    "CPU pair lease retained identity differs from its guard")
        self.kernel.validate_current()
        directory = self._validate_directory()
        require(self.directory_fd is not None,
                "CPU pair lease directory descriptor was lost")
        opened_directory = os.fstat(self.directory_fd)
        descriptor = os.fstat(self.descriptor)
        path = os.stat(
            self.path.name, dir_fd=self.directory_fd, follow_symlinks=False)
        require(stat.S_ISREG(descriptor.st_mode) and
                descriptor.st_uid == os.getuid() and descriptor.st_nlink == 1 and
                stat.S_IMODE(descriptor.st_mode) == 0o600 and
                (descriptor.st_dev, descriptor.st_ino) == (path.st_dev, path.st_ino) and
                (descriptor.st_dev, descriptor.st_ino) ==
                    (self.identity["device"], self.identity["inode"]),
                "CPU pair lease path was replaced or its metadata changed")
        require((directory.st_dev, directory.st_ino) ==
                (opened_directory.st_dev, opened_directory.st_ino) ==
                (self.identity["directory_device"],
                 self.identity["directory_inode"]),
                "CPU pair lease directory was replaced")
        require(typed_equal(self.identity["authority"], self.kernel.identity),
                "CPU pair lease authority identity changed")
        expected = canonical_bytes(pair_lease_payload(self.cpu, self.sibling))
        retained, _metadata = _snapshot_regular_fd(
            self.descriptor, 4096, "CPU pair lease")
        require(retained == expected, "CPU pair lease contents changed while held")

    def __enter__(self) -> dict[str, Any]:
        require(self.anchor is not None, "CPU pair stable anchor is unavailable")
        if self.owns_anchor:
            self.anchor.__enter__()
        else:
            self.anchor.validate_current()
        try:
            authority = self.kernel.acquire()
            self._validate_runtime_root()
            try:
                self.root.mkdir(mode=0o700)
                os.chmod(self.root, 0o700)
            except FileExistsError:
                pass
            self._validate_directory()
            self.directory_fd, _unused = _open_directory_nofollow(
                self.root, "CPU pair lease directory")
            flags = os.O_RDWR | getattr(os, "O_CLOEXEC", 0) | \
                getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
            try:
                self.descriptor = os.open(
                    self.path.name, flags | os.O_CREAT | os.O_EXCL, 0o600,
                    dir_fd=self.directory_fd)
                os.fchmod(self.descriptor, 0o600)
            except FileExistsError:
                self.descriptor = os.open(
                    self.path.name, flags, dir_fd=self.directory_fd)
            metadata = os.fstat(self.descriptor)
            require(stat.S_ISREG(metadata.st_mode) and
                    metadata.st_uid == os.getuid() and metadata.st_nlink == 1 and
                    stat.S_IMODE(metadata.st_mode) == 0o600,
                    "CPU pair lease file has unsafe ownership, type, links, or permissions")
            try:
                fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise EvidenceError(
                    f"physical CPU pair is already leased: {self.path}") from error
            expected = canonical_bytes(pair_lease_payload(self.cpu, self.sibling))
            retained, _opened = _snapshot_regular_fd(
                self.descriptor, 4096, "CPU pair lease")
            if not retained:
                require(os.write(self.descriptor, expected) == len(expected),
                        "short write to CPU pair lease")
                os.fsync(self.descriptor)
                retained, _opened = _snapshot_regular_fd(
                    self.descriptor, 4096, "CPU pair lease")
            require(retained == expected,
                    "CPU pair lease has unexpected or noncanonical contents")
            directory = os.fstat(self.directory_fd)
            self.identity = {
                "schema": PAIR_LEASE_SCHEMA,
                "authority": authority,
                "device": metadata.st_dev,
                "directory_device": directory.st_dev,
                "directory_inode": directory.st_ino,
                "inode": metadata.st_ino,
                "lock": "exclusive_nonblocking_pair_wide",
                "path": str(_lexical_absolute(self.path)),
                "payload": pair_lease_payload(self.cpu, self.sibling),
                "sha256": sha256_bytes(retained),
            }
            self._identity_bytes = canonical_bytes(self.identity)
            self.validate_current()
            return copy.deepcopy(self.identity)
        except Exception:
            try:
                try:
                    if self.descriptor is not None:
                        descriptor = self.descriptor
                        self.descriptor = None
                        try:
                            fcntl.flock(descriptor, fcntl.LOCK_UN)
                        finally:
                            os.close(descriptor)
                finally:
                    if self.directory_fd is not None:
                        directory = self.directory_fd
                        self.directory_fd = None
                        os.close(directory)
            finally:
                try:
                    self.kernel.release()
                finally:
                    if self.owns_anchor and self.anchor is not None:
                        self.anchor.__exit__(None, None, None)
            raise

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        try:
            try:
                if self.descriptor is not None:
                    descriptor = self.descriptor
                    self.descriptor = None
                    try:
                        fcntl.flock(descriptor, fcntl.LOCK_UN)
                    finally:
                        os.close(descriptor)
            finally:
                if self.directory_fd is not None:
                    directory = self.directory_fd
                    self.directory_fd = None
                    os.close(directory)
        finally:
            try:
                self.kernel.release()
            finally:
                if self.owns_anchor and self.anchor is not None:
                    self.anchor.__exit__(exc_type, exc, tb)


def validate_pair_lease_identity(value: object, cpu: int,
                                 sibling: int) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "authority", "device", "directory_device", "directory_inode",
        "inode", "lock",
        "path", "payload", "sha256"}, "CPU pair lease identity is incomplete")
    payload = value["payload"]
    require(isinstance(payload, dict) and set(payload) == {"cpus", "schema", "uid"} and
            payload["schema"] == PAIR_LEASE_SCHEMA and
            typed_equal(payload["cpus"], sorted((cpu, sibling))) and
            type(payload["uid"]) is int and payload["uid"] >= 0,
            "CPU pair lease payload does not match the campaign")
    expected_path = pair_lease_directory(payload["uid"]) / \
        pair_lease_name(cpu, sibling, payload["uid"])
    require(value["path"] == str(expected_path),
            "CPU pair lease path does not identify the reserved pair")
    require(value["schema"] == PAIR_LEASE_SCHEMA,
            "CPU pair lease schema changed")
    validate_kernel_lease_identity(value["authority"], "pair", payload)
    require(all(type(value[name]) is int and value[name] >= 0 for name in
                ("device", "directory_device", "directory_inode", "inode")),
            "CPU pair lease filesystem identity is invalid")
    require(value["lock"] == "exclusive_nonblocking_pair_wide" and
            value["sha256"] == sha256_bytes(canonical_bytes(payload)),
            "CPU pair lease lock or digest is invalid")
    return value


def isolation_record(cpu: int, sibling: int, pair_lease: Mapping[str, Any],
                     before_cpu: Mapping[str, Any], after_cpu: Mapping[str, Any],
                     before_sibling: Mapping[str, Any],
                     after_sibling: Mapping[str, Any], started_ns: int,
                     ended_ns: int) -> dict[str, Any]:
    before_cpu = validate_cpu_stat_snapshot(
        before_cpu, cpu, "before timing CPU")
    after_cpu = validate_cpu_stat_snapshot(after_cpu, cpu, "after timing CPU")
    before_sibling = validate_cpu_stat_snapshot(
        before_sibling, sibling, "before sibling CPU")
    after_sibling = validate_cpu_stat_snapshot(
        after_sibling, sibling, "after sibling CPU")
    require(type(started_ns) is int and type(ended_ns) is int,
            "isolation monotonic timestamps are invalid")
    benchmark_delta = cpu_stat_delta(before_cpu, after_cpu)
    sibling_delta = cpu_stat_delta(before_sibling, after_sibling)
    accepted = (ended_ns > started_ns and benchmark_delta["nonidle_jiffies"] > 0 and
                sibling_delta["total_jiffies"] > 0 and
                sibling_delta["nonidle_jiffies"] == 0)
    return {
        "schema": ISOLATION_SCHEMA, "accepted": accepted,
        "timing_cpu": cpu, "sibling_cpu": sibling,
        "before": {"timing_cpu": dict(before_cpu),
                   "sibling_cpu": dict(before_sibling),
                   "monotonic_ns": started_ns},
        "after": {"timing_cpu": dict(after_cpu),
                  "sibling_cpu": dict(after_sibling),
                  "monotonic_ns": ended_ns},
        "delta": {"timing_cpu": benchmark_delta,
                  "sibling_cpu": sibling_delta},
        "duration_ns": ended_ns - started_ns,
        "pair_lease": dict(pair_lease),
        "policy": {"counter_source": "/proc/stat",
                   "idle_fields": list(CPU_STAT_IDLE_FIELDS),
                   "nonidle_fields": [name for name in CPU_STAT_FIELDS
                                      if name not in CPU_STAT_IDLE_FIELDS],
                   "sibling_max_nonidle_jiffies": 0,
                   "timing_min_nonidle_jiffies": 1,
                   "sibling_min_total_jiffies": 1},
    }


def validate_isolation(value: object, cpu: int, sibling: int, *,
                       require_accepted: bool = True) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "accepted", "timing_cpu", "sibling_cpu", "before", "after",
        "delta", "duration_ns", "pair_lease", "policy"},
        "isolation schema changed")
    require(value["schema"] == ISOLATION_SCHEMA and
            type(value["timing_cpu"]) is int and
            type(value["sibling_cpu"]) is int and
            value["timing_cpu"] == cpu and value["sibling_cpu"] == sibling,
            "isolation CPU pair changed")
    validate_pair_lease_identity(value["pair_lease"], cpu, sibling)
    before, after = value["before"], value["after"]
    require(isinstance(before, dict) and isinstance(after, dict) and
            set(before) == {"timing_cpu", "sibling_cpu", "monotonic_ns"} and
            set(after) == {"timing_cpu", "sibling_cpu", "monotonic_ns"} and
            type(before["monotonic_ns"]) is int and
            type(after["monotonic_ns"]) is int,
            "isolation snapshots changed")
    validate_cpu_stat_snapshot(before["timing_cpu"], cpu,
                               "before timing CPU")
    validate_cpu_stat_snapshot(after["timing_cpu"], cpu,
                               "after timing CPU")
    validate_cpu_stat_snapshot(before["sibling_cpu"], sibling,
                               "before sibling CPU")
    validate_cpu_stat_snapshot(after["sibling_cpu"], sibling,
                               "after sibling CPU")
    expected = isolation_record(
        cpu, sibling, value["pair_lease"], before["timing_cpu"],
        after["timing_cpu"], before["sibling_cpu"], after["sibling_cpu"],
        before["monotonic_ns"], after["monotonic_ns"])
    require(typed_equal(value, expected), "isolation derivation changed")
    if require_accepted:
        require(value["accepted"] is True,
                "timing CPU was inactive or SMT sibling was active during timing")
    return value


def parse_reservation(raw: bytes, cpu: int, sibling: int) -> dict[str, Any]:
    payload = strict_json(raw, "CPU reservation")
    require(isinstance(payload, dict) and set(payload) == {
        "benchmark_cpu", "nonce", "owner", "reserved_sibling", "schema", "status"},
        "CPU reservation fields changed")
    require(raw == canonical_bytes(payload),
            "CPU reservation is not canonical JSON without trailing newline")
    require(payload["schema"] == RESERVATION_SCHEMA and payload["status"] == "held",
            "CPU reservation is stale or not held")
    require(type(payload["benchmark_cpu"]) is int and
            type(payload["reserved_sibling"]) is int and
            payload["benchmark_cpu"] == cpu and payload["reserved_sibling"] == sibling,
            "CPU reservation pair differs")
    require(isinstance(payload["owner"], str) and payload["owner"].strip(),
            "CPU reservation owner is empty")
    require(isinstance(payload["nonce"], str) and 8 <= len(payload["nonce"]) <= 256,
            "CPU reservation nonce is invalid")
    return payload


def reservation_authority_payload(path: Path, cpu: int, sibling: int,
                                  reservation_sha256: str | None = None,
                                  uid: int | None = None) -> dict[str, Any]:
    retained_uid = os.getuid() if uid is None else uid
    result: dict[str, Any] = {
        "cpus": sorted((cpu, sibling)), "path": str(_lexical_absolute(path)),
        "uid": retained_uid,
    }
    if reservation_sha256 is not None:
        result["reservation_sha256"] = reservation_sha256
    return result


class Reservation:
    def __init__(self, path: Path, cpu: int, sibling: int,
                 anchor: StableLeaseAnchor | None = None):
        self.path = _lexical_absolute(path)
        self.cpu = cpu
        self.sibling = sibling
        self.fd: int | None = None
        self.parent_fd: int | None = None
        self.name = self.path.name
        self.raw = b""
        self.identity: dict[str, Any] | None = None
        self._identity_bytes: bytes | None = None
        self.anchor = anchor
        self.owns_anchor = anchor is None
        if self.owns_anchor:
            self.anchor = StableLeaseAnchor()
        self.kernel = KernelNamespaceLease(
            "reservation", reservation_authority_payload(
                self.path, cpu, sibling))

    def __enter__(self) -> dict[str, Any]:
        require(self.anchor is not None, "reservation stable anchor is unavailable")
        if self.owns_anchor:
            self.anchor.__enter__()
        else:
            self.anchor.validate_current()
        try:
            self.kernel.acquire()
            self.parent_fd, self.name, _absolute = _open_parent_nofollow(self.path)
            directory = os.fstat(self.parent_fd)
            require(stat.S_ISDIR(directory.st_mode),
                    "CPU reservation parent is not a directory")
            self.fd = os.open(
                self.name, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0),
                dir_fd=self.parent_fd)
            descriptor = os.fstat(self.fd)
            current = os.stat(
                self.name, dir_fd=self.parent_fd, follow_symlinks=False)
            require(stat.S_ISREG(descriptor.st_mode) and
                    descriptor.st_uid == os.getuid() and descriptor.st_nlink == 1 and
                    stat.S_IMODE(descriptor.st_mode) & 0o022 == 0 and
                    (descriptor.st_dev, descriptor.st_ino) ==
                    (current.st_dev, current.st_ino),
                    "CPU reservation has unsafe ownership, type, links, or permissions")
            try:
                fcntl.flock(self.fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise EvidenceError("CPU reservation is already locked") from error
            self.raw, descriptor = _snapshot_regular_fd(
                self.fd, MAX_LOG_BYTES, "CPU reservation")
            payload = parse_reservation(self.raw, self.cpu, self.sibling)
            reservation_sha = sha256_bytes(self.raw)
            authority_payload = reservation_authority_payload(
                self.path, self.cpu, self.sibling, reservation_sha)
            # Upgrade the pre-open namespace binding with the immutable byte digest.
            # The path/pair namespace remains held; this second name binds its content.
            content_kernel = KernelNamespaceLease("reservation", authority_payload)
            content_authority = content_kernel.acquire()
            self.content_kernel = content_kernel
            self.identity = {
                "schema": RESERVATION_SCHEMA, "bytes": len(self.raw),
                "sha256": reservation_sha, "payload": payload,
                "lock": "linux_abstract_bind_plus_flock",
                "path": str(self.path), "uid": os.getuid(),
                "device": descriptor.st_dev, "inode": descriptor.st_ino,
                "directory_device": directory.st_dev,
                "directory_inode": directory.st_ino,
                "authority": {"path": copy.deepcopy(self.kernel.identity),
                              "content": content_authority},
            }
            self._identity_bytes = canonical_bytes(self.identity)
            self.validate_current()
            return copy.deepcopy(self.identity)
        except Exception:
            self._release()
            raise

    def validate_current(self, expected_identity: Mapping[str, Any] | None = None) -> None:
        require(self.fd is not None and self.parent_fd is not None and
                self.identity is not None and hasattr(self, "content_kernel") and
                self.anchor is not None,
                "CPU reservation lock was lost")
        self.anchor.validate_current()
        require(self._identity_bytes == canonical_bytes(self.identity),
                "CPU reservation immutable identity changed")
        if expected_identity is not None:
            require(typed_equal(expected_identity, self.identity),
                    "CPU reservation retained identity differs from its guard")
        self.kernel.validate_current()
        self.content_kernel.validate_current()
        descriptor = os.fstat(self.fd)
        directory = os.fstat(self.parent_fd)
        current = os.stat(self.name, dir_fd=self.parent_fd, follow_symlinks=False)
        live_parent, live_name, _absolute = _open_parent_nofollow(self.path)
        try:
            live_directory = os.fstat(live_parent)
            live_path = os.stat(
                live_name, dir_fd=live_parent, follow_symlinks=False)
        finally:
            os.close(live_parent)
        require((descriptor.st_dev, descriptor.st_ino) ==
                (current.st_dev, current.st_ino) ==
                (live_path.st_dev, live_path.st_ino) ==
                (self.identity["device"], self.identity["inode"]) and
                (directory.st_dev, directory.st_ino) ==
                (live_directory.st_dev, live_directory.st_ino) ==
                (self.identity["directory_device"],
                 self.identity["directory_inode"]) and
                descriptor.st_uid == os.getuid() and descriptor.st_nlink == 1 and
                stat.S_ISREG(descriptor.st_mode) and
                stat.S_IMODE(descriptor.st_mode) & 0o022 == 0,
                "CPU reservation path was replaced or its metadata changed")
        raw, _metadata = _snapshot_regular_fd(
            self.fd, MAX_LOG_BYTES, "CPU reservation")
        require(raw == self.raw, "CPU reservation changed while locked")
        parse_reservation(raw, self.cpu, self.sibling)
        expected_authority_payload = reservation_authority_payload(
            self.path, self.cpu, self.sibling, self.identity["sha256"])
        require(typed_equal(self.identity["authority"]["path"],
                            self.kernel.identity) and
                typed_equal(self.identity["authority"]["content"],
                            self.content_kernel.identity),
                "CPU reservation authority identity changed")
        validate_kernel_lease_identity(
            self.identity["authority"]["path"], "reservation",
            reservation_authority_payload(self.path, self.cpu, self.sibling))
        validate_kernel_lease_identity(
            self.identity["authority"]["content"], "reservation",
            expected_authority_payload)

    def _release(self) -> None:
        try:
            try:
                if self.fd is not None:
                    descriptor = self.fd
                    self.fd = None
                    try:
                        fcntl.flock(descriptor, fcntl.LOCK_UN)
                    finally:
                        os.close(descriptor)
            finally:
                if self.parent_fd is not None:
                    parent = self.parent_fd
                    self.parent_fd = None
                    os.close(parent)
        finally:
            try:
                try:
                    if hasattr(self, "content_kernel"):
                        self.content_kernel.release()
                finally:
                    self.kernel.release()
            finally:
                if self.owns_anchor and self.anchor is not None:
                    self.anchor.__exit__(None, None, None)

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        self._release()


def validate_reservation_record(value: object, cpu: int, sibling: int) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "bytes", "sha256", "payload", "lock", "path", "uid",
        "device", "inode", "directory_device", "directory_inode", "authority"} and
        value["schema"] == RESERVATION_SCHEMA and
        value["lock"] == "linux_abstract_bind_plus_flock",
        "reservation record changed")
    require(type(value["bytes"]) is int and 0 < value["bytes"] <= MAX_LOG_BYTES and
            isinstance(value["sha256"], str) and
            SHA256_RE.fullmatch(value["sha256"]) is not None,
            "reservation byte identity is invalid")
    raw = canonical_bytes(value["payload"])
    require(value["bytes"] == len(raw) and value["sha256"] == sha256_bytes(raw),
            "reservation bytes or checksum changed")
    require(parse_reservation(raw, cpu, sibling) == value["payload"],
            "reservation semantics changed")
    require(isinstance(value["path"], str) and os.path.isabs(value["path"]) and
            type(value["uid"]) is int and value["uid"] >= 0 and
            all(type(value[name]) is int and value[name] >= 0 for name in
                ("device", "inode", "directory_device", "directory_inode")),
            "reservation opened-object identity is invalid")
    authority_payload = reservation_authority_payload(
        Path(value["path"]), cpu, sibling, value["sha256"], value["uid"])
    require(isinstance(value["authority"], dict) and
            set(value["authority"]) == {"path", "content"},
            "reservation authority schema changed")
    validate_kernel_lease_identity(
        value["authority"]["path"], "reservation",
        reservation_authority_payload(
            Path(value["path"]), cpu, sibling, uid=value["uid"]))
    validate_kernel_lease_identity(
        value["authority"]["content"], "reservation", authority_payload)
    return value


def finite(value: object, label: str) -> float:
    return bounded_positive_number(value, label, 1e300)


def validate_summary(cell: Mapping[str, Any], summary_key: str,
                     samples_key: str) -> dict[str, float]:
    samples_value = cell.get(samples_key)
    require(isinstance(samples_value, list) and len(samples_value) == SAMPLE_COUNT,
            f"{samples_key} does not contain seven samples")
    samples = [finite(value, samples_key) for value in samples_value]
    median = statistics.median(samples)
    mad = statistics.median(abs(value - median) for value in samples)
    summary = cell.get(summary_key)
    require(isinstance(summary, dict) and set(summary) == {"median_us", "mad_us"},
            f"{summary_key} schema changed")
    require(type(summary["median_us"]) in (int, float) and
            type(summary["mad_us"]) in (int, float) and
            float(summary["median_us"]) == median and float(summary["mad_us"]) == mad,
            f"{summary_key} differs from raw samples")
    return {"median_us": median, "mad_us": mad}


def validate_child_result(value: object, cpu: int, inputs: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "affinity", "allocation_tracking", "benchmarks", "core_git_sha",
        "correctness", "library_sha256", "omp_dynamic", "omp_num_threads",
        "production_constructor_rejected", "profile", "requested_backend",
        "runtime_backend", "sanitizer", "sanitizer_features", "schema",
        "source_sha256", "status", "timing_scope"}, "C7 result schema changed")
    require(value["schema"] == CHILD_SCHEMA and value["status"] == "pass" and
            typed_equal(value["profile"], EXPECTED_PROFILE) and
            value["production_constructor_rejected"] is True and
            value["timing_scope"] == "candidate-authoritative-requires-runner-manifest" and
            value["requested_backend"] == "avx2" and
            value["runtime_backend"] == "avx2" and
            typed_equal(value["affinity"], [cpu]) and
            value["omp_num_threads"] == "1" and value["omp_dynamic"] == "FALSE" and
            value["sanitizer"] == "none" and
            typed_equal(value["sanitizer_features"], {"address": False, "undefined": False}) and
            value["allocation_tracking"] == "global-new" and
            typed_equal(value["correctness"], EXPECTED_CORRECTNESS),
            "C7 result identity, backend, affinity, or correctness changed")
    require(value["source_sha256"] == inputs["source"]["sha256"] and
            value["library_sha256"] == inputs["archive"]["sha256"] and
            value["core_git_sha"] == inputs["git"]["core_commit"],
            "C7 embedded source/archive/Git identity changed")
    cells = value["benchmarks"]
    require(isinstance(cells, list) and len(cells) == len(EXPECTED_CELLS),
            "C7 result does not contain the full 12-cell matrix")
    normalized = []
    for index, (cell, expected) in enumerate(zip(cells, EXPECTED_CELLS)):
        require(isinstance(cell, dict) and set(cell) == CELL_KEYS,
                f"C7 cell {index} schema changed")
        k, r, byte_count, batch, losses, exact_field, padded_field = expected
        geometry = [cell[key] for key in
                    ("K", "R", "bytes", "batch", "losses", "exact_field", "padded_field")]
        require(typed_equal(geometry, list(expected)), f"C7 cell {index} order/geometry changed")
        require(type(cell["exact_coefficients"]) is int and
                cell["exact_coefficients"] == k * r and
                type(cell["exact_decode_terms"]) is int and
                cell["exact_decode_terms"] == k * losses,
                f"C7 cell {index} exact accounting changed")
        require(type(cell["padded_encode_scratch"]) is int and
                0 <= cell["padded_encode_scratch"] <=
                    MAX_REPORTED_SCRATCH_BYTES and
                type(cell["padded_decode_scratch"]) is int and
                0 <= cell["padded_decode_scratch"] <=
                    MAX_REPORTED_SCRATCH_BYTES,
                f"C7 cell {index} scratch accounting is invalid")
        summaries = {name: validate_summary(cell, name, samples)
                     for name, samples in SUMMARY_SAMPLE_PAIRS}
        normalized.append({"geometry": list(expected), "summaries": summaries})
    return {"schema": "leopard2-c7-authoritative-validation/v1",
            "cell_count": len(cells), "sample_count_per_metric": SAMPLE_COUNT,
            "correctness_checksum": value["correctness"]["digest_fnv64"],
            "cells": normalized}


def artifact_relative_path(relative: object) -> Path:
    require(isinstance(relative, str) and relative and
            "\\" not in relative and ":" not in relative and
            not os.path.isabs(relative),
            "artifact path is not relative")
    relative_path = Path(relative)
    require(relative_path.as_posix() == relative,
            "artifact path is not canonical")
    parts = relative_path.parts
    require(len(parts) <= MAX_ARTIFACT_DEPTH and
            len(os.fsencode(relative)) <= MAX_ARTIFACT_PATH_BYTES and
            all(part not in ("", ".", "..") for part in parts),
            "artifact path is not canonical")
    return relative_path


def _open_directory_nofollow(path: Path, label: str) -> tuple[int, Path]:
    absolute = _lexical_absolute(path)
    try:
        parent, name, _unused = _open_parent_nofollow(absolute)
        try:
            descriptor = os.open(
                name, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=parent)
        finally:
            os.close(parent)
    except OSError as error:
        raise EvidenceError(f"cannot securely open {label}: {error}") from error
    metadata = os.fstat(descriptor)
    require(stat.S_ISDIR(metadata.st_mode), f"{label} is not a directory")
    return descriptor, absolute


def _artifact_snapshot(root: Path, relative: object, maximum_bytes: int,
                       label: str) -> tuple[bytes, os.stat_result]:
    relative_path = artifact_relative_path(relative)
    root_descriptor, absolute_root = _open_directory_nofollow(
        root, "artifact root")
    directories = [root_descriptor]
    components: list[tuple[str, os.stat_result]] = []
    descriptor: int | None = None
    try:
        for component in relative_path.parts[:-1]:
            child = os.open(
                component, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=directories[-1])
            before = os.fstat(child)
            current = os.stat(
                component, dir_fd=directories[-1], follow_symlinks=False)
            require(stat.S_ISDIR(before.st_mode) and
                    (before.st_dev, before.st_ino) ==
                    (current.st_dev, current.st_ino),
                    f"{label} directory component was replaced")
            components.append((component, before))
            directories.append(child)
        descriptor = os.open(
            relative_path.parts[-1], os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0),
            dir_fd=directories[-1])
        data, metadata = _snapshot_regular_fd(
            descriptor, maximum_bytes, label)
        current = os.stat(
            relative_path.parts[-1], dir_fd=directories[-1],
            follow_symlinks=False)
        require((current.st_dev, current.st_ino, current.st_mode,
                 current.st_nlink, current.st_size, current.st_mtime_ns,
                 current.st_ctime_ns) ==
                (metadata.st_dev, metadata.st_ino, metadata.st_mode,
                 metadata.st_nlink, metadata.st_size, metadata.st_mtime_ns,
                 metadata.st_ctime_ns), f"{label} path was replaced while read")
        for index, (component, opened) in enumerate(components):
            live = os.stat(
                component, dir_fd=directories[index], follow_symlinks=False)
            require((live.st_dev, live.st_ino, live.st_mode,
                     live.st_mtime_ns, live.st_ctime_ns) ==
                    (opened.st_dev, opened.st_ino, opened.st_mode,
                     opened.st_mtime_ns, opened.st_ctime_ns),
                    f"{label} directory path was replaced while read")
        live_root = os.stat(absolute_root, follow_symlinks=False)
        opened_root = os.fstat(root_descriptor)
        require((live_root.st_dev, live_root.st_ino) ==
                (opened_root.st_dev, opened_root.st_ino),
                f"{label} artifact root was replaced while read")
        return data, metadata
    except OSError as error:
        raise EvidenceError(f"cannot securely read {label}: {error}") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
        for directory in reversed(directories):
            os.close(directory)


def artifact_record(root: Path, path: Path) -> dict[str, Any]:
    root = _lexical_absolute(root)
    path = _lexical_absolute(path)
    try:
        relative = path.relative_to(root).as_posix()
    except ValueError as error:
        raise EvidenceError("artifact path escapes evidence root") from error
    data, metadata = _artifact_snapshot(
        root, relative, MAX_ARTIFACT_BYTES, "retained artifact")
    return {"path": relative, "size": metadata.st_size,
            "sha256": sha256_bytes(data)}


def validate_artifact_record(record: object, label: str,
                             maximum_bytes: int = MAX_ARTIFACT_BYTES) -> dict[str, Any]:
    require(isinstance(record, dict) and set(record) == {"path", "size", "sha256"},
            f"{label} artifact record changed")
    require(type(record["size"]) is int and 0 <= record["size"] <= maximum_bytes and
            isinstance(record["sha256"], str) and
            SHA256_RE.fullmatch(record["sha256"]) is not None,
            f"{label} artifact identity is invalid")
    artifact_relative_path(record["path"])
    return record


def read_artifact(root: Path, record: object, label: str,
                  maximum_bytes: int = MAX_ARTIFACT_BYTES) -> bytes:
    record = validate_artifact_record(record, label, maximum_bytes)
    data, _metadata = _artifact_snapshot(
        root, record["path"], maximum_bytes, label)
    require(record["size"] == len(data) and record["sha256"] == sha256_bytes(data),
            f"retained {label} checksum changed")
    return data


def artifact_inventory(root: Path) -> list[dict[str, Any]]:
    """Return the exact staged inventory, excluding either terminal record."""
    descriptor, absolute_root = _open_directory_nofollow(root, "artifact root")
    records: list[dict[str, Any]] = []
    seen_directories: set[tuple[int, int]] = set()
    node_count = 0

    def scan(directory: int, prefix: tuple[str, ...]) -> None:
        nonlocal node_count
        directory_before = os.fstat(directory)
        identity = (directory_before.st_dev, directory_before.st_ino)
        require(identity not in seen_directories,
                "artifact directory graph contains a loop")
        seen_directories.add(identity)
        names: list[str] = []
        with os.scandir(directory) as iterator:
            for entry in iterator:
                node_count += 1
                require(node_count <= MAX_ARTIFACT_COUNT * 2,
                        "artifact tree contains too many entries")
                names.append(entry.name)
        for name in sorted(names):
            require(name not in ("", ".", "..") and "/" not in name and
                    "\\" not in name,
                    "artifact tree contains a noncanonical name")
            relative_parts = (*prefix, name)
            relative = "/".join(relative_parts)
            artifact_relative_path(relative)
            before = os.stat(name, dir_fd=directory, follow_symlinks=False)
            require(not stat.S_ISLNK(before.st_mode),
                    f"retained artifact is a symbolic link: {relative}")
            if stat.S_ISDIR(before.st_mode):
                child = os.open(
                    name, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_DIRECTORY", 0) |
                    getattr(os, "O_NOFOLLOW", 0), dir_fd=directory)
                try:
                    opened = os.fstat(child)
                    require((opened.st_dev, opened.st_ino) ==
                            (before.st_dev, before.st_ino),
                            "artifact directory was replaced while scanning")
                    scan(child, relative_parts)
                    after = os.stat(
                        name, dir_fd=directory, follow_symlinks=False)
                    require((after.st_dev, after.st_ino, after.st_mtime_ns,
                             after.st_ctime_ns) ==
                            (before.st_dev, before.st_ino, before.st_mtime_ns,
                             before.st_ctime_ns),
                            "artifact directory changed while scanning")
                finally:
                    os.close(child)
                continue
            require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1,
                    f"retained artifact is not a single-link regular file: {relative}")
            if relative in ("failure.json", "manifest.json"):
                continue
            require(len(records) < MAX_ARTIFACT_COUNT,
                    "artifact inventory contains too many files")
            child = os.open(
                name, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0),
                dir_fd=directory)
            try:
                data, opened = _snapshot_regular_fd(
                    child, MAX_ARTIFACT_BYTES, relative)
                after = os.stat(name, dir_fd=directory, follow_symlinks=False)
                require((after.st_dev, after.st_ino, after.st_mode,
                         after.st_nlink, after.st_size, after.st_mtime_ns,
                         after.st_ctime_ns) ==
                        (opened.st_dev, opened.st_ino, opened.st_mode,
                         opened.st_nlink, opened.st_size, opened.st_mtime_ns,
                         opened.st_ctime_ns),
                        "artifact changed while scanning")
                records.append({"path": relative, "size": len(data),
                                "sha256": sha256_bytes(data)})
            finally:
                os.close(child)
        directory_after = os.fstat(directory)
        require((directory_after.st_dev, directory_after.st_ino,
                 directory_after.st_mtime_ns, directory_after.st_ctime_ns) ==
                (directory_before.st_dev, directory_before.st_ino,
                 directory_before.st_mtime_ns, directory_before.st_ctime_ns),
                "artifact directory changed while scanning")

    try:
        scan(descriptor, ())
        current = os.stat(absolute_root, follow_symlinks=False)
        opened = os.fstat(descriptor)
        require((current.st_dev, current.st_ino) ==
                (opened.st_dev, opened.st_ino),
                "artifact root was replaced while scanning")
        return sorted(records, key=lambda item: item["path"])
    except OSError as error:
        raise EvidenceError(f"cannot securely inventory artifacts: {error}") from error
    finally:
        os.close(descriptor)


def validate_artifact_inventory(value: object, root: Path | None,
                                check_files: bool) -> list[dict[str, Any]]:
    require(isinstance(value, list), "failure artifact inventory is not a list")
    records: list[dict[str, Any]] = []
    previous = ""
    for index, item in enumerate(value):
        record = validate_artifact_record(item, f"failure inventory item {index}")
        path = record["path"]
        require(path not in ("failure.json", "manifest.json") and
                (index == 0 or previous < path),
                "failure artifact inventory is not strictly sorted and unique")
        previous = path
        records.append(record)
    if check_files:
        require(root is not None, "failure artifact root is absent")
        require(typed_equal(records, artifact_inventory(root)),
                "failure artifact inventory differs from retained files")
    return records


def inventory_contains(inventory: Sequence[Mapping[str, Any]],
                       record: Mapping[str, Any], label: str) -> None:
    require(any(typed_equal(item, record) for item in inventory),
            f"{label} is absent from the exact artifact inventory")


def retain_build_provenance(root: Path, build_manifest: Path,
                            inputs: Mapping[str, Any]) -> dict[str, Any]:
    data = read_bounded(_lexical_absolute(build_manifest))
    identity = inputs["build_manifest"]
    require(len(data) == identity["bytes"] and
            sha256_bytes(data) == identity["sha256"],
            "C7 build manifest changed before retention")
    manifest_path = root / "snapshot/provenance/build-run-manifest-v4.json"
    write_exclusive(manifest_path, data)
    runner_data = read_bounded(_lexical_absolute(Path(__file__)))
    require(len(runner_data) == inputs["runner"]["bytes"] and
            sha256_bytes(runner_data) == inputs["runner"]["sha256"],
            "authoritative runner changed before retention")
    runner_path = root / "snapshot/provenance/run_authoritative.py"
    write_exclusive(runner_path, runner_data)
    return {"schema": BUILD_PROVENANCE_SCHEMA,
            "manifest": artifact_record(root, manifest_path),
            "runner": artifact_record(root, runner_path),
            "attestation": copy.deepcopy(inputs["build_attestation"])}


def validate_build_provenance(value: object, root: Path,
                              inputs: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "manifest", "runner", "attestation"} and
        value["schema"] == BUILD_PROVENANCE_SCHEMA,
            "retained C7 build provenance schema changed")
    require(isinstance(value["manifest"], dict) and
            isinstance(value["runner"], dict) and
            value["manifest"].get("path") ==
                "snapshot/provenance/build-run-manifest-v4.json" and
            value["runner"].get("path") ==
                "snapshot/provenance/run_authoritative.py",
            "retained C7 build provenance paths changed")
    data = read_artifact(root, value["manifest"], "C7 v4 build manifest")
    require(data == canonical_pretty_bytes(
        strict_json(data, "retained C7 v4 build manifest")),
        "retained C7 v4 build manifest is not canonical pretty JSON")
    require(value["manifest"]["size"] == inputs["build_manifest"]["bytes"] and
            value["manifest"]["sha256"] == inputs["build_manifest"]["sha256"],
            "retained C7 build manifest differs from the live-validated input")
    runner_data = read_artifact(root, value["runner"], "authoritative runner")
    require(len(runner_data) == inputs["runner"]["bytes"] and
            value["runner"]["sha256"] == inputs["runner"]["sha256"],
            "retained authoritative runner differs from its input identity")
    parsed = strict_json(data, "retained C7 v4 build manifest")
    expected = derive_build_attestation(
        parsed, inputs["build_manifest"], inputs["git"]["tooling_commit"],
        inputs["git"]["core_commit"], inputs["source"], inputs["build_runner"],
        inputs["build_validator"], inputs["archive"], inputs["executable"],
        inputs["core_source_closure"], inputs["taskset"], inputs["python"])
    require(typed_equal(value["attestation"], expected) and
            typed_equal(inputs["build_attestation"], expected),
            "retained C7 v4 build attestation changed")
    return value


def expected_stderr() -> bytes:
    return b"".join(f"C7 benchmark {index}/12\n".encode("ascii")
                    for index in range(1, 13))


def _child_resource_limits() -> None:
    resource.setrlimit(
        resource.RLIMIT_FSIZE, (MAX_LOG_BYTES, MAX_LOG_BYTES))


def _kill_process_group_and_reap(process: subprocess.Popen[bytes]) -> None:
    """Kill the isolated child session and bound the leader reap operation."""
    deadline = time.monotonic() + CHILD_REAP_TIMEOUT_SECONDS
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    except OSError as error:
        raise EvidenceError(f"cannot kill child process group: {error}") from error
    if process.poll() is None:
        try:
            process.wait(timeout=max(0.001, deadline - time.monotonic()))
        except subprocess.TimeoutExpired as error:
            raise EvidenceError(
                "child process group leader did not terminate after SIGKILL") from error
    while True:
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            break
        except OSError as error:
            raise EvidenceError(
                f"cannot verify child process-group teardown: {error}") from error
        remaining = deadline - time.monotonic()
        require(remaining > 0,
                "child process group retained descendants after SIGKILL")
        time.sleep(min(0.01, remaining))


def _wait_isolated_child(
        process: subprocess.Popen[bytes], timeout_seconds: float) -> tuple[int, bool]:
    timed_out = False
    try:
        try:
            returncode = process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            returncode = 124
    finally:
        _kill_process_group_and_reap(process)
    return returncode, timed_out


def run_child_bounded(command: Sequence[str], stdout_path: Path,
                      stderr_path: Path, timeout_seconds: float,
                      environment: Mapping[str, str]) -> tuple[
                          subprocess.CompletedProcess[bytes], bool, int, int]:
    """Run without PIPE accumulation; inherited RLIMIT_FSIZE bounds every output."""
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | \
        getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    stdout_fd = os.open(stdout_path, flags, 0o600)
    stderr_fd: int | None = None
    process: subprocess.Popen[bytes] | None = None
    cleanup_attempted = False
    started = time.monotonic_ns()
    timed_out = False
    try:
        stderr_fd = os.open(stderr_path, flags, 0o600)
        process = subprocess.Popen(
            list(command), stdin=subprocess.DEVNULL, stdout=stdout_fd,
            stderr=stderr_fd, env=dict(environment), close_fds=True,
            preexec_fn=_child_resource_limits, start_new_session=True)
        cleanup_attempted = True
        returncode, timed_out = _wait_isolated_child(process, timeout_seconds)
        os.fsync(stdout_fd)
        os.fsync(stderr_fd)
        ended = time.monotonic_ns()
        completed = subprocess.CompletedProcess(
            list(command), returncode, stdout=b"", stderr=b"")
        return completed, timed_out, started, ended
    finally:
        if process is not None and not cleanup_attempted:
            cleanup_attempted = True
            _kill_process_group_and_reap(process)
        if stderr_fd is not None:
            os.close(stderr_fd)
        os.close(stdout_fd)


def run_capture_bounded(
        command: Sequence[str], *, cwd: Path | None, timeout_seconds: float,
        environment: Mapping[str, str]) -> tuple[
            subprocess.CompletedProcess[bytes], bytes, bytes, bool]:
    """Bound diagnostic command output without ever accumulating a PIPE."""
    with tempfile.TemporaryDirectory(prefix="leopard2-c7-capture-") as directory:
        root = Path(directory)
        stdout_path = root / "stdout.bin"
        stderr_path = root / "stderr.bin"
        retained_command = list(command)
        if cwd is not None:
            # Popen's cwd is needed for Git commands.  Use a tiny wrapper-free
            # local implementation so the same file-backed capture policy holds.
            flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | \
                getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
            stdout_fd = os.open(stdout_path, flags, 0o600)
            stderr_fd = os.open(stderr_path, flags, 0o600)
            process: subprocess.Popen[bytes] | None = None
            cleanup_attempted = False
            timed_out = False
            try:
                process = subprocess.Popen(
                    retained_command, cwd=str(cwd), stdin=subprocess.DEVNULL,
                    stdout=stdout_fd, stderr=stderr_fd,
                    env=dict(environment), close_fds=True,
                    preexec_fn=_child_resource_limits, start_new_session=True)
                cleanup_attempted = True
                returncode, timed_out = _wait_isolated_child(
                    process, timeout_seconds)
                os.fsync(stdout_fd)
                os.fsync(stderr_fd)
            finally:
                if process is not None and not cleanup_attempted:
                    cleanup_attempted = True
                    _kill_process_group_and_reap(process)
                os.close(stderr_fd)
                os.close(stdout_fd)
            completed = subprocess.CompletedProcess(
                retained_command, returncode, stdout=b"", stderr=b"")
        else:
            completed, timed_out, _started, _ended = run_child_bounded(
                retained_command, stdout_path, stderr_path, timeout_seconds,
                environment)
        stdout = read_bounded(stdout_path, MAX_LOG_BYTES)
        stderr = read_bounded(stderr_path, MAX_LOG_BYTES)
        return completed, stdout, stderr, timed_out


def bounded_positive_number(value: object, label: str,
                            maximum: float) -> float:
    require(type(value) in (int, float), f"{label} is not numeric")
    try:
        result = float(value)
    except (OverflowError, TypeError, ValueError) as error:
        raise EvidenceError(f"{label} cannot be represented") from error
    require(math.isfinite(result) and 0 < result <= maximum,
            f"{label} is outside the supported positive range")
    return result


def validate_request(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "backend", "cell_geometry", "child_environment", "command", "cpu",
        "sibling", "timeout_seconds"} and value["backend"] == "avx2" and
        typed_equal(value["cell_geometry"],
                    [list(cell) for cell in EXPECTED_CELLS]) and
        typed_equal(value["child_environment"], CHILD_ENVIRONMENT),
        "authoritative request changed")
    cpu, sibling = value["cpu"], value["sibling"]
    require(type(cpu) is int and type(sibling) is int and
            0 <= cpu <= MAX_CPU_ID and 0 <= sibling <= MAX_CPU_ID and
            cpu != sibling,
            "authoritative request CPU is invalid")
    bounded_positive_number(
        value["timeout_seconds"], "authoritative timeout", MAX_TIMEOUT_SECONDS)
    require(typed_equal(value["command"], [
        "${TASKSET}", "-c", str(cpu), "${C7_EXECUTABLE}",
        "--backend", "avx2", "${RESULT_JSON}"]),
        "authoritative child command changed")
    return value


PUBLICATION_CONTEXT_FIELDS = (
    "arguments", "request", "inputs_before", "inputs_after", "host_before",
    "host_after", "reservation", "pair_lease", "isolation",
    "build_provenance", "child", "validated_output", "lifecycle",
)


def publication_state_payload(created_utc: str,
                              context: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema": PUBLICATION_STATE_SCHEMA, "created_utc": created_utc,
        **{key: copy.deepcopy(context[key]) for key in PUBLICATION_CONTEXT_FIELDS},
    }


def validate_publication_state(value: object) -> dict[str, Any]:
    state = verify_signature(value, "C7 publication state")
    require(set(state) == {"schema", "created_utc", "digest",
                           *PUBLICATION_CONTEXT_FIELDS} and
            state["schema"] == PUBLICATION_STATE_SCHEMA,
            "publication state schema changed")
    validate_utc(state["created_utc"], "publication state creation time")
    validate_lifecycle(state["lifecycle"], require_success=True)
    return state


def validate_raw(value: object, root: Path, check_files: bool = True,
                 result_override: object | None = None) -> dict[str, Any]:
    raw = verify_signature(value, "C7 raw evidence")
    require(set(raw) == {"schema", "created_utc", "arguments", "request",
                         "publication_state", "pair_lease", "lifecycle", "inputs_before",
                         "inputs_after", "host_before", "host_after", "reservation",
                         "isolation", "build_provenance", "child",
                         "validated_output", "digest"} and
            raw["schema"] == RAW_SCHEMA, "raw evidence schema changed")
    validate_utc(raw["created_utc"], "raw creation time")
    arguments = validate_arguments(raw["arguments"])
    request = validate_request(raw["request"])
    cpu, sibling = request["cpu"], request["sibling"]
    require(cpu == arguments["cpu"] and sibling == arguments["sibling"] and
            typed_equal(request["timeout_seconds"],
                        arguments["timeout_seconds"]),
            "raw arguments differ from the request")
    validate_lifecycle(raw["lifecycle"], require_success=True)
    inputs = validate_input_snapshot(raw["inputs_before"])
    require(typed_equal(raw["inputs_after"], inputs), "input identity changed during run")
    validate_build_provenance(raw["build_provenance"], root, inputs)
    validate_host(raw["host_before"], cpu, sibling)
    validate_host(raw["host_after"], cpu, sibling)
    require(typed_equal(raw["host_before"], raw["host_after"]),
            "host topology/frequency policy changed during run")
    validate_reservation_record(raw["reservation"], cpu, sibling)
    pair_lease = validate_pair_lease_identity(raw["pair_lease"], cpu, sibling)
    validate_isolation(raw["isolation"], cpu, sibling)
    require(typed_equal(raw["isolation"]["pair_lease"], pair_lease),
            "raw isolation differs from its pair lease")
    child = raw["child"]
    require(isinstance(child, dict) and set(child) == {
        "command", "environment", "returncode", "timed_out", "duration_ns",
        "stdout", "stderr", "result"} and child["command"] == request["command"] and
        typed_equal(child["environment"], CHILD_ENVIRONMENT) and
        type(child["returncode"]) is int and child["returncode"] == 0 and
        child["timed_out"] is False and
        type(child["duration_ns"]) is int and child["duration_ns"] > 0 and
        child["duration_ns"] == raw["isolation"]["duration_ns"],
        "C7 child failed, timed out, or changed invocation")
    if check_files:
        stdout = read_artifact(root, child["stdout"], "stdout", MAX_LOG_BYTES)
        stderr = read_artifact(root, child["stderr"], "stderr", MAX_LOG_BYTES)
        result_bytes = read_artifact(root, child["result"], "result")
        require(stdout == b"", "C7 child unexpectedly wrote stdout")
        require(stderr == expected_stderr(), "C7 child stderr progress changed")
        parsed = strict_json(result_bytes, "C7 result")
    else:
        require(isinstance(result_override, dict), "parsed C7 result is absent")
        parsed = result_override
    normalized = validate_child_result(parsed, cpu, inputs)
    require(typed_equal(raw["validated_output"], normalized),
            "retained C7 validation summary changed")
    if check_files:
        state_record = validate_artifact_record(
            raw["publication_state"], "publication state",
            MAX_TOP_LEVEL_JSON_BYTES)
        require(state_record["path"] == "snapshot/publication-state.json",
                "publication state path changed")
        state_bytes = read_artifact(
            root, state_record, "publication state", MAX_TOP_LEVEL_JSON_BYTES)
        state = strict_json(state_bytes, "C7 publication state")
        require(state_bytes == canonical_bytes(state) + b"\n",
                "publication state is not canonical JSON")
        validate_publication_state(state)
        expected = signed(publication_state_payload(
            raw["created_utc"], {key: raw[key]
                                 for key in PUBLICATION_CONTEXT_FIELDS}))
        require(typed_equal(state, expected),
                "raw evidence differs from immutable publication state")
    return normalized


def validate_manifest(value: object, root: Path) -> dict[str, Any]:
    _require_terminal_absent(root, "failure.json")
    manifest = verify_signature(value, "C7 manifest")
    require(set(manifest) == {"schema", "created_utc", "valid", "raw",
                              "publication_state", "artifact_inventory",
                              "arguments", "request", "inputs", "reservation",
                              "pair_lease", "isolation", "build_provenance",
                              "lifecycle", "validated_output", "digest"} and
            manifest["schema"] == MANIFEST_SCHEMA and manifest["valid"] is True,
            "C7 manifest schema changed")
    validate_utc(manifest["created_utc"], "manifest creation time")
    raw_bytes = read_artifact(root, manifest["raw"], "raw")
    raw = strict_json(raw_bytes, "C7 raw evidence")
    require(raw_bytes == canonical_bytes(raw) + b"\n",
            "raw evidence is not canonical JSON")
    normalized = validate_raw(raw, root, check_files=True)
    inventory = validate_artifact_inventory(
        manifest["artifact_inventory"], root, check_files=True)
    inventory_contains(inventory, manifest["raw"], "raw evidence")
    inventory_contains(
        inventory, manifest["publication_state"], "publication state")
    require(manifest["created_utc"] == raw["created_utc"],
            "manifest creation time differs from raw evidence")
    for key, expected in (("arguments", raw["arguments"]),
                          ("request", raw["request"]),
                          ("inputs", raw["inputs_before"]),
                          ("reservation", raw["reservation"]),
                          ("pair_lease", raw["pair_lease"]),
                          ("isolation", raw["isolation"]),
                          ("build_provenance", raw["build_provenance"]),
                          ("lifecycle", raw["lifecycle"]),
                          ("validated_output", normalized)):
        require(typed_equal(manifest[key], expected), f"manifest {key} differs from raw")
    require(typed_equal(manifest["publication_state"], raw["publication_state"]),
            "manifest publication state differs from raw")
    return normalized


def validate_failure_child(value: object, request: Mapping[str, Any],
                           root: Path | None, check_files: bool) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "command", "environment", "returncode", "timed_out", "duration_ns",
        "stdout", "stderr", "result"},
        "failure child record changed")
    require(typed_equal(value["command"], request["command"]) and
            typed_equal(value["environment"], CHILD_ENVIRONMENT) and
            type(value["returncode"]) is int and
            type(value["timed_out"]) is bool and
            type(value["duration_ns"]) is int and value["duration_ns"] > 0,
            "failure child invocation or outcome is invalid")
    validate_artifact_record(value["stdout"], "failure stdout", MAX_LOG_BYTES)
    validate_artifact_record(value["stderr"], "failure stderr", MAX_LOG_BYTES)
    require(value["stdout"]["path"] == "snapshot/child/stdout.bin" and
            value["stderr"]["path"] == "snapshot/child/stderr.bin",
            "failure child stream paths changed")
    require(value["result"] is None or isinstance(value["result"], dict),
            "failure result artifact record changed")
    if value["result"] is not None:
        validate_artifact_record(value["result"], "failure result")
        require(value["result"]["path"] == "snapshot/child/result.json",
                "failure child result path changed")
    if check_files:
        require(root is not None, "failure artifact root is absent")
        read_artifact(root, value["stdout"], "failure stdout", MAX_LOG_BYTES)
        read_artifact(root, value["stderr"], "failure stderr", MAX_LOG_BYTES)
        if value["result"] is not None:
            read_artifact(root, value["result"], "failure result")
    return value


def completed_failure_stage(index: int) -> dict[str, Any]:
    require(type(index) is int and -1 <= index < len(FAILURE_STAGES),
            "failure stage index is invalid")
    return {"index": index,
            "name": "none" if index == -1 else FAILURE_STAGES[index]}


def validate_completed_failure_stage(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {"index", "name"} and
            type(value["index"]) is int and
            -1 <= value["index"] < len(FAILURE_STAGES) and
            value["name"] == ("none" if value["index"] == -1 else
                              FAILURE_STAGES[value["index"]]),
            "failure completed stage changed")
    return value


def validate_arguments(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "cpu", "sibling", "timeout_seconds", "expected_tooling_commit",
        "expected_core_commit"}, "failure argument record changed")
    require(type(value["cpu"]) is int and type(value["sibling"]) is int and
            0 <= value["cpu"] <= MAX_CPU_ID and
            0 <= value["sibling"] <= MAX_CPU_ID and
            value["cpu"] != value["sibling"] and
            isinstance(value["expected_tooling_commit"], str) and
            GIT_SHA_RE.fullmatch(value["expected_tooling_commit"]) is not None and
            isinstance(value["expected_core_commit"], str) and
            GIT_SHA_RE.fullmatch(value["expected_core_commit"]) is not None,
            "failure arguments are invalid")
    bounded_positive_number(
        value["timeout_seconds"], "failure timeout", MAX_TIMEOUT_SECONDS)
    return value


def failure_state_payload(completed_stage: Mapping[str, Any],
                          failure_code: str,
                          context: Mapping[str, Any],
                          lifecycle: Mapping[str, Any] | None = None,
                          publication: Mapping[str, Any] | None = None,
                          snapshot_inventory: Mapping[str, Any] | None = None) -> dict[str, Any]:
    return {
        "schema": FAILURE_STATE_SCHEMA,
        "completed_stage": copy.deepcopy(completed_stage),
        "failure_code": failure_code,
        "lifecycle": copy.deepcopy(
            empty_lifecycle() if lifecycle is None else lifecycle),
        "publication": copy.deepcopy(
            {"state": None, "raw": None} if publication is None else publication),
        "snapshot_inventory": copy.deepcopy(snapshot_inventory),
        **{key: copy.deepcopy(context[key]) for key in FAILURE_CONTEXT_FIELDS},
    }


def validate_failure_state(value: object) -> dict[str, Any]:
    state = verify_signature(value, "C7 failure state")
    require(set(state) == {
        "schema", "completed_stage", "failure_code", "lifecycle",
        "publication", "snapshot_inventory", "digest",
        *FAILURE_CONTEXT_FIELDS} and state["schema"] == FAILURE_STATE_SCHEMA and
        isinstance(state["failure_code"], str) and
        state["failure_code"] in FAILURE_CODE_STAGE,
        "failure state schema changed")
    stage = validate_completed_failure_stage(state["completed_stage"])
    require(FAILURE_CODE_STAGE[state["failure_code"]] == stage["index"],
            "failure code does not identify the completed stage")
    validate_lifecycle(state["lifecycle"], require_success=False)
    validate_publication_records(state["publication"])
    validate_snapshot_inventory(state["snapshot_inventory"])
    for index, key in enumerate(FAILURE_CONTEXT_FIELDS):
        require((state[key] is not None) == (stage["index"] >= index),
                f"failure state presence mask changed at {key}")
    return state


def validate_publication_records(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {"state", "raw"},
            "failure publication record changed")
    if value["state"] is not None:
        validate_artifact_record(
            value["state"], "failed publication state", MAX_TOP_LEVEL_JSON_BYTES)
        require(value["state"]["path"] == "snapshot/publication-state.json",
                "failed publication state path changed")
    if value["raw"] is not None:
        require(value["state"] is not None,
                "failed raw publication has no publication state")
        validate_artifact_record(
            value["raw"], "failed raw evidence", MAX_TOP_LEVEL_JSON_BYTES)
        require(value["raw"]["path"] == "snapshot/raw.json",
                "failed raw evidence path changed")
    return value


def validate_snapshot_inventory(value: object) -> dict[str, Any] | None:
    if value is None:
        return None
    require(isinstance(value, dict) and set(value) == {
        "expected", "observed", "predicate"} and
        value["predicate"] in (
            "inventory-mismatch", "reservation-guard-validation",
            "pair-guard-validation", "snapshot-revalidation"),
            "snapshot mismatch inventory changed")
    for name in ("expected", "observed"):
        records = value[name]
        require(isinstance(records, list) and len(records) <= MAX_ARTIFACT_COUNT,
                f"snapshot {name} inventory is invalid")
        previous = ""
        for index, record in enumerate(records):
            validate_artifact_record(record, f"snapshot {name} item {index}")
            require((index == 0 or previous < record["path"]) and
                    record["path"] not in (
                        "manifest.json", "failure.json", "failure/state-v2.json"),
                    f"snapshot {name} inventory order changed")
            previous = record["path"]
    if value["predicate"] == "inventory-mismatch":
        require(not typed_equal(value["expected"], value["observed"]),
                "snapshot inventory-mismatch predicate is false")
    else:
        require(typed_equal(value["expected"], value["observed"]),
                "guard/semantic snapshot predicate contains an inventory mismatch")
    return value


def _failure_result_bytes(root: Path, child: Mapping[str, Any]) -> bytes:
    require(child["result"] is not None, "failure result artifact is absent")
    return read_artifact(root, child["result"], "failure result")


def _validate_failure(value: object, root: Path | None,
                      check_files: bool) -> dict[str, Any]:
    require(check_files and root is not None,
            "v4 failure validation requires its retained artifact directory")
    _require_terminal_absent(root, "manifest.json")
    failure = verify_signature(value, "C7 failure evidence")
    require(set(failure) == {
        "schema", "status", "created_utc", "failure_code", "completed_stage",
        "error_type", "error", "traceback", "state", "artifact_inventory",
        "lifecycle", "publication", "snapshot_inventory", "digest",
        *FAILURE_CONTEXT_FIELDS} and
        failure["schema"] == FAILURE_SCHEMA and failure["status"] == "failed" and
        isinstance(failure["failure_code"], str) and
        failure["failure_code"] in FAILURE_CODE_STAGE and
        isinstance(failure["error_type"], str) and failure["error_type"] and
        isinstance(failure["error"], str) and failure["error"] and
        isinstance(failure["traceback"], str) and failure["traceback"] and
        len(failure["traceback"].encode("utf-8", errors="replace")) <=
            MAX_DIAGNOSTIC_BYTES,
        "failure evidence schema changed")
    validate_utc(failure["created_utc"], "failure creation time")
    stage = validate_completed_failure_stage(failure["completed_stage"])
    code = failure["failure_code"]
    require(FAILURE_CODE_STAGE[code] == stage["index"],
            "failure code does not identify the completed stage")
    require(failure["error_type"] == "EvidenceError",
            "failure error type is not canonical")
    if code == "child-exit":
        pass  # bound to the retained return code after child validation below
    else:
        require(code in FAILURE_DIAGNOSTICS and
                failure["error"] == FAILURE_DIAGNOSTICS[code],
                "failure diagnostic does not match its code")
    lifecycle = validate_lifecycle(
        failure["lifecycle"], require_success=False)
    publication = validate_publication_records(failure["publication"])
    snapshot_inventory = validate_snapshot_inventory(
        failure["snapshot_inventory"])
    require((snapshot_inventory is not None) ==
            (code == "final-snapshot-invalid"),
            "snapshot mismatch presence does not match its failure code")
    for index, key in enumerate(FAILURE_CONTEXT_FIELDS):
        require((failure[key] is not None) == (stage["index"] >= index),
                f"failure context presence mask changed at {key}")

    inventory = validate_artifact_inventory(
        failure["artifact_inventory"], root, check_files=True)
    state_record = validate_artifact_record(
        failure["state"], "failure state", MAX_TOP_LEVEL_JSON_BYTES)
    require(state_record["path"] == "failure/state-v2.json",
            "failure state path changed")
    inventory_contains(inventory, state_record, "failure state")
    state_bytes = read_artifact(
        root, state_record, "failure state", MAX_TOP_LEVEL_JSON_BYTES)
    state = strict_json(state_bytes, "C7 failure state")
    require(state_bytes == canonical_bytes(state) + b"\n",
            "failure state is not canonical JSON")
    validate_failure_state(state)
    expected_state = signed(failure_state_payload(
        stage, code, {key: failure[key] for key in FAILURE_CONTEXT_FIELDS},
        lifecycle, publication, snapshot_inventory))
    require(typed_equal(state, expected_state),
            "failure top-level context differs from retained state")

    expected_inventory = [state_record]
    if failure["build_provenance"] is not None:
        expected_inventory.extend((failure["build_provenance"]["manifest"],
                                   failure["build_provenance"]["runner"]))
    if failure["child"] is not None:
        expected_inventory.extend((failure["child"]["stdout"],
                                   failure["child"]["stderr"]))
        if failure["child"]["result"] is not None:
            expected_inventory.append(failure["child"]["result"])
    for retained in (publication["state"], publication["raw"]):
        if retained is not None:
            expected_inventory.append(retained)
    if snapshot_inventory is None:
        require(typed_equal(
            sorted(expected_inventory, key=lambda item: item["path"]), inventory),
            "failure inventory is not the exact stage-derived artifact set")
    else:
        require(typed_equal(
            sorted([*snapshot_inventory["observed"], state_record],
                   key=lambda item: item["path"]), inventory),
            "final-snapshot failure inventory differs from its observation")

    arguments = None
    request = None
    inputs_before = None
    reservation = None
    pair_lease = None
    provenance = None
    child = None
    isolation = None
    parsed_result: object | None = None

    if stage["index"] >= 0:
        arguments = validate_arguments(failure["arguments"])
    if stage["index"] >= 1:
        request = validate_request(failure["request"])
        require(request["cpu"] == arguments["cpu"] and
                request["sibling"] == arguments["sibling"] and
                typed_equal(request["timeout_seconds"],
                            arguments["timeout_seconds"]),
                "failure request differs from validated arguments")
    if stage["index"] >= 2:
        validate_host(failure["host_before"], request["cpu"], request["sibling"])
    if stage["index"] >= 3:
        reservation = validate_reservation_record(
            failure["reservation"], request["cpu"], request["sibling"])
    if stage["index"] >= 4:
        pair_lease = validate_pair_lease_identity(
            failure["pair_lease"], request["cpu"], request["sibling"])
    if stage["index"] >= 5:
        inputs_before = validate_input_snapshot(failure["inputs_before"])
        require(inputs_before["git"]["tooling_commit"] ==
                    arguments["expected_tooling_commit"] and
                inputs_before["git"]["core_commit"] ==
                    arguments["expected_core_commit"],
                "failure inputs differ from expected Git commits")
    if stage["index"] >= 6:
        provenance = validate_build_provenance(
            failure["build_provenance"], root, inputs_before)
        inventory_contains(inventory, provenance["manifest"],
                           "failure build manifest")
        inventory_contains(inventory, provenance["runner"],
                           "failure authoritative runner")
    if stage["index"] >= 7:
        child = validate_failure_child(
            failure["child"], request, root,
            check_files=code != "final-snapshot-invalid")
        if code != "final-snapshot-invalid":
            inventory_contains(inventory, child["stdout"], "failure child stdout")
            inventory_contains(inventory, child["stderr"], "failure child stderr")
            if child["result"] is not None:
                inventory_contains(inventory, child["result"], "failure child result")
    if code == "final-snapshot-invalid":
        require(stage["index"] == 11 and
                lifecycle["reservation"]["exit"] == "pass" and
                lifecycle["pair_lease"]["exit"] == "pass" and
                lifecycle["affinity"]["restore_succeeded"] and
                lifecycle["affinity"]["exact"],
                "final-snapshot lifecycle predicate is false")
        require(typed_equal(
            snapshot_inventory["expected"],
            sorted(expected_inventory[1:], key=lambda item: item["path"])),
            "final-snapshot expected inventory is not context-derived")
        checks = lifecycle["guard_checks"]
        predicate = snapshot_inventory["predicate"]
        if predicate == "reservation-guard-validation":
            require(checks["reservation_final"]["status"] == "failed" and
                    checks["pair_final"]["status"] == "not-reached",
                    "final reservation-guard predicate is false")
        elif predicate == "pair-guard-validation":
            require(checks["reservation_final"]["status"] == "pass" and
                    checks["pair_final"]["status"] == "failed",
                    "final pair-guard predicate is false")
        else:
            require(checks["reservation_final"]["status"] == "pass" and
                    checks["pair_final"]["status"] == "pass",
                    "final snapshot revalidation predicate is false")
        return failure
    if stage["index"] >= 8:
        isolation = validate_isolation(
            failure["isolation"], request["cpu"], request["sibling"],
            require_accepted=False)
        require(typed_equal(isolation["pair_lease"], pair_lease) and
                child["duration_ns"] == isolation["duration_ns"],
                "failure child/isolation/lease relationship changed")
        checks = lifecycle["guard_checks"]

        if code == "child-timeout":
            require(checks["reservation_post_child"]["status"] == "not-reached" and
                    checks["pair_post_child"]["status"] == "not-reached",
                    "timeout failure ran post-child guard checks")
            require(child["timed_out"] is True and child["returncode"] == 124 and
                    failure["error"] == FAILURE_DIAGNOSTICS[code],
                    "timeout failure predicate changed")
            return failure
        require(child["timed_out"] is False,
                "non-timeout failure retains a timed-out child")
        if code == "child-exit":
            require(checks["reservation_post_child"]["status"] == "not-reached" and
                    checks["pair_post_child"]["status"] == "not-reached",
                    "child-exit failure ran post-child guard checks")
            require(child["returncode"] != 0 and
                    failure["error_type"] == "EvidenceError" and
                    failure["error"] == f"C7 child exited {child['returncode']}",
                    "child-exit failure predicate changed")
            return failure
        require(child["returncode"] == 0,
                "post-success failure retains a nonzero child exit")
        if code == "child-result-missing":
            require(checks["reservation_post_child"]["status"] == "not-reached" and
                    checks["pair_post_child"]["status"] == "not-reached",
                    "missing-result failure ran post-child guard checks")
            require(child["result"] is None and
                    failure["error"] == FAILURE_DIAGNOSTICS[code],
                    "missing-result failure predicate changed")
            return failure
        require(child["result"] is not None,
                "post-success failure has no retained result")

        stdout = read_artifact(root, child["stdout"], "failure stdout", MAX_LOG_BYTES)
        if code == "child-stdout-invalid":
            require(checks["reservation_post_child"]["status"] == "not-reached" and
                    checks["pair_post_child"]["status"] == "not-reached",
                    "stdout failure ran post-child guard checks")
            require(stdout != b"", "stdout failure predicate is false")
            return failure
        require(stdout == b"", "post-success failure has unexpected stdout")
        stderr = read_artifact(root, child["stderr"], "failure stderr", MAX_LOG_BYTES)
        if code == "child-stderr-invalid":
            require(checks["reservation_post_child"]["status"] == "not-reached" and
                    checks["pair_post_child"]["status"] == "not-reached",
                    "stderr failure ran post-child guard checks")
            require(stderr != expected_stderr(), "stderr failure predicate is false")
            return failure
        require(stderr == expected_stderr(),
                "post-success failure has invalid stderr progress")
        if code == "isolation-rejected":
            require(checks["reservation_post_child"]["status"] == "not-reached" and
                    checks["pair_post_child"]["status"] == "not-reached",
                    "isolation failure ran post-child guard checks")
            require(isolation["accepted"] is False,
                    "isolation failure predicate is false")
            return failure
        require(isolation["accepted"] is True,
                "post-success failure retains rejected isolation")

        if code == "reservation-post-child-invalid":
            require(checks["reservation_post_child"]["status"] == "failed" and
                    checks["pair_post_child"]["status"] == "not-reached",
                    "reservation post-child guard predicate is false")
            return failure
        if code == "lease-post-child-invalid":
            require(checks["reservation_post_child"]["status"] == "pass" and
                    checks["pair_post_child"]["status"] == "failed",
                    "pair post-child guard predicate is false")
            return failure

        require(checks["reservation_post_child"]["status"] == "pass" and
                checks["pair_post_child"]["status"] == "pass",
                "post-guard failure lacks completed post-child checks")

        result_bytes = _failure_result_bytes(root, child)
        if code == "child-result-json-invalid":
            try:
                strict_json(result_bytes, "C7 result")
            except EvidenceError:
                return failure
            raise EvidenceError("result-JSON failure predicate is false")
        parsed_result = strict_json(result_bytes, "C7 result")

    if stage["index"] >= 9:
        inputs_after = validate_input_snapshot(failure["inputs_after"])
        if code == "inputs-drift":
            require(not typed_equal(inputs_after, inputs_before),
                    "input-drift failure predicate is false")
            return failure
        require(typed_equal(inputs_after, inputs_before),
                "post-input failure retains input drift")
    if stage["index"] >= 10:
        validate_host(failure["host_after"], request["cpu"], request["sibling"])
        if code == "host-drift":
            require(not typed_equal(failure["host_after"], failure["host_before"]),
                    "host-drift failure predicate is false")
            return failure
        require(typed_equal(failure["host_after"], failure["host_before"]),
                "post-host failure retains host drift")
        if code == "child-result-invalid":
            try:
                validate_child_result(parsed_result, request["cpu"], inputs_before)
            except EvidenceError:
                return failure
            raise EvidenceError("child-result failure predicate is false")
    if stage["index"] >= 11:
        normalized = validate_child_result(
            parsed_result, request["cpu"], inputs_before)
        require(typed_equal(failure["validated_output"], normalized),
                "failure validated output differs from the child result")
        guard_failed = any(lifecycle[name]["exit"] == "failed"
                           for name in ("reservation", "pair_lease"))
        affinity = lifecycle["affinity"]
        require(all(check["status"] == "pass" for check in
                    lifecycle["guard_checks"].values()),
                "post-validation failure lacks all guard checks")
        if code == "guard-exit-failed":
            require(guard_failed and affinity["restore_succeeded"] and
                    affinity["exact"], "guard-exit failure predicate is false")
            return failure
        require(not guard_failed,
                "non-guard failure retains a failed guard teardown")
        if code == "affinity-restore-failed":
            require(affinity["restore_attempted"] and
                    (not affinity["restore_succeeded"] or not affinity["exact"]),
                    "affinity-restore failure predicate is false")
            return failure
        require(affinity["restore_succeeded"] and affinity["exact"],
                "publication failure lacks exact affinity restoration")
        if code == "publication-state-write-failed":
            require(publication["state"] is None and publication["raw"] is None,
                    "publication-state failure predicate is false")
            return failure
        require(publication["state"] is not None,
                "post-publication-state failure lacks its staged state")
        if code == "raw-write-failed":
            require(publication["raw"] is None,
                    "raw-write failure predicate is false")
            return failure
        require(publication["raw"] is not None,
                "post-raw failure lacks staged raw evidence")
        require(code in ("manifest-validation-failed",
                         "final-snapshot-invalid", "manifest-write-failed"),
                "stage-11 failure has no exact lifecycle predicate")
    return failure


def validate_failure(value: object, root: Path | None = None,
                     check_files: bool = False) -> dict[str, Any]:
    try:
        return _validate_failure(value, root, check_files)
    except EvidenceError:
        raise
    except (AttributeError, IndexError, KeyError, OverflowError, RecursionError,
            TypeError) as error:
        raise EvidenceError(
            f"malformed C7 failure evidence: {type(error).__name__}: {error}") from error


def _failure_diagnostic(code: str, child: Mapping[str, Any] | None) -> str:
    if code == "child-exit":
        require(child is not None and type(child.get("returncode")) is int,
                "child-exit has no retained return code")
        return f"C7 child exited {child['returncode']}"
    require(code in FAILURE_DIAGNOSTICS, "failure code has no canonical diagnostic")
    return FAILURE_DIAGNOSTICS[code]


def _bounded_traceback(value: str) -> str:
    encoded = (value or "failure captured without a Python traceback\n").encode(
        "utf-8", errors="replace")
    if len(encoded) > MAX_DIAGNOSTIC_BYTES:
        encoded = encoded[:MAX_DIAGNOSTIC_BYTES]
    return encoded.decode("utf-8", errors="replace")


def _restore_affinity(lifecycle: dict[str, Any], original: set[int] | None) -> None:
    if original is None:
        return
    affinity = lifecycle["affinity"]
    affinity["restore_attempted"] = True
    try:
        os.sched_setaffinity(0, original)
        observed = set(os.sched_getaffinity(0))
        affinity["observed"] = sorted(observed)
        affinity["restore_succeeded"] = True
        affinity["exact"] = observed == original
        if observed != original:
            affinity["restore_succeeded"] = False
            affinity["error_type"] = "EvidenceError"
            affinity["error"] = "affinity readback differs from captured CPU set"
    except BaseException as error:
        error_type, message = _bounded_exception(error)
        affinity.update(restore_succeeded=False, observed=None, exact=False,
                        error_type=error_type, error=message)


def _root_identity(descriptor: int, root: Path, *,
                   require_private_owner: bool = False) -> tuple[int, int]:
    opened = os.fstat(descriptor)
    current = os.stat(_lexical_absolute(root), follow_symlinks=False)
    mode = stat.S_IMODE(opened.st_mode)
    require(stat.S_ISDIR(opened.st_mode) and
            (opened.st_dev, opened.st_ino) == (current.st_dev, current.st_ino),
            "evidence root path or inode changed")
    if require_private_owner:
        require(opened.st_uid == os.getuid() and mode & 0o022 == 0,
                "evidence root ownership or mode is unsafe for publication")
    return opened.st_dev, opened.st_ino


def _terminal_absent(descriptor: int, name: str) -> None:
    try:
        os.stat(name, dir_fd=descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return
    raise EvidenceError(f"terminal artifact already exists: {name}")


def _require_terminal_absent(root: Path, name: str) -> None:
    descriptor, absolute = _open_directory_nofollow(root, "evidence root")
    try:
        _root_identity(descriptor, absolute)
        _terminal_absent(descriptor, name)
        _root_identity(descriptor, absolute)
    finally:
        os.close(descriptor)


def _write_json_exclusive_at(
        descriptor: int, name: str, value: Mapping[str, Any]) -> tuple[int, int]:
    require(name in ("manifest.json", "failure.json"),
            "terminal artifact name is invalid")
    data = canonical_bytes(value) + b"\n"
    child: int | None = None
    identity: tuple[int, int] | None = None
    completed = False
    try:
        child = os.open(
            name, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
            getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
            0o600, dir_fd=descriptor)
        opened = os.fstat(child)
        identity = opened.st_dev, opened.st_ino
        offset = 0
        while offset < len(data):
            written = os.write(child, data[offset:])
            require(written > 0, f"short terminal write: {name}")
            offset += written
        os.fsync(child)
        retained = os.fstat(child)
        current = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
        require(stat.S_ISREG(retained.st_mode) and retained.st_nlink == 1 and
                (retained.st_dev, retained.st_ino, retained.st_size) ==
                    (current.st_dev, current.st_ino, current.st_size) and
                retained.st_size == len(data),
                f"terminal artifact changed while written: {name}")
        completed = True
        return identity
    except FileExistsError as error:
        raise EvidenceError(f"refusing to replace terminal artifact: {name}") from error
    except OSError as error:
        raise EvidenceError(f"cannot securely write terminal artifact {name}: {error}") \
            from error
    finally:
        if child is not None:
            os.close(child)
        if not completed and identity is not None:
            try:
                current = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
                if (current.st_dev, current.st_ino) == identity:
                    os.unlink(name, dir_fd=descriptor)
            except OSError:
                pass


def _remove_created_terminal(
        descriptor: int, name: str, identity: tuple[int, int]) -> None:
    try:
        current = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
        if (current.st_dev, current.st_ino) == identity:
            os.unlink(name, dir_fd=descriptor)
            os.fsync(descriptor)
    except OSError:
        pass


def _publish_terminal(root: Path, name: str, value: Mapping[str, Any],
                      expected_root: tuple[int, int] | None = None) -> None:
    require(name in ("manifest.json", "failure.json"),
            "terminal artifact name is invalid")
    other = "failure.json" if name == "manifest.json" else "manifest.json"
    descriptor, absolute = _open_directory_nofollow(root, "evidence root")
    authority = KernelNamespaceLease(
        "terminal", {"path": str(_lexical_absolute(root)), "uid": os.getuid()})
    created: tuple[int, int] | None = None
    acquired = False
    try:
        opened_root = _root_identity(
            descriptor, absolute, require_private_owner=True)
        if expected_root is not None:
            require(opened_root == expected_root,
                    "evidence root differs from its creation inode")
        authority.acquire()
        acquired = True
        require(_root_identity(
                    descriptor, absolute, require_private_owner=True) == opened_root,
                "evidence root changed before terminal publication")
        _terminal_absent(descriptor, name)
        _terminal_absent(descriptor, other)
        created = _write_json_exclusive_at(descriptor, name, value)
        _terminal_absent(descriptor, other)
        require(_root_identity(
                    descriptor, absolute, require_private_owner=True) == opened_root,
                "evidence root changed during terminal publication")
        os.fsync(descriptor)
        require(_root_identity(
                    descriptor, absolute, require_private_owner=True) == opened_root,
                "evidence root changed after terminal publication")
    except BaseException:
        if created is not None:
            _remove_created_terminal(descriptor, name, created)
        raise
    finally:
        if acquired:
            try:
                authority.release()
            except OSError:
                pass
        os.close(descriptor)


def run_campaign(options: argparse.Namespace) -> int:
    output = _lexical_absolute(options.output)
    require(not output.exists(), f"output already exists: {output}")
    output.mkdir(parents=True, mode=0o700)
    os.chmod(output, 0o700)
    output_descriptor, _output_absolute = _open_directory_nofollow(
        output, "evidence root")
    try:
        output_identity = _root_identity(
            output_descriptor, output, require_private_owner=True)
    finally:
        os.close(output_descriptor)
    lifecycle = empty_lifecycle()
    publication: dict[str, Any] = {"state": None, "raw": None}
    arguments: dict[str, Any] | None = None
    child_record: dict[str, Any] | None = None
    isolation: dict[str, Any] | None = None
    build_provenance: dict[str, Any] | None = None
    request: dict[str, Any] | None = None
    inputs_before: dict[str, Any] | None = None
    inputs_after: dict[str, Any] | None = None
    host_before: dict[str, Any] | None = None
    host_after: dict[str, Any] | None = None
    reservation_record: dict[str, Any] | None = None
    pair_lease_record: dict[str, Any] | None = None
    validated_output: dict[str, Any] | None = None
    original_affinity: set[int] | None = None
    completed_stage_index = -1
    failure_code = "arguments-invalid"
    body_error: BaseException | None = None
    body_traceback = ""
    reservation_guard: Reservation | None = None
    pair_guard: PairLease | None = None
    final_snapshot_predicate = "snapshot-revalidation"
    try:
        arguments = validate_arguments({
            "cpu": options.cpu, "sibling": options.sibling,
            "timeout_seconds": options.timeout,
            "expected_tooling_commit": options.expected_tooling_commit,
            "expected_core_commit": options.expected_core_commit,
        })
        completed_stage_index = 0
        failure_code = "request-invalid"
        request = {
            "backend": "avx2", "cell_geometry": [list(cell) for cell in EXPECTED_CELLS],
            "child_environment": dict(CHILD_ENVIRONMENT), "cpu": options.cpu,
            "sibling": options.sibling, "timeout_seconds": options.timeout,
            "command": ["${TASKSET}", "-c", str(options.cpu),
                        "${C7_EXECUTABLE}", "--backend", "avx2", "${RESULT_JSON}"],
        }
        validate_request(request)
        completed_stage_index = 1
        failure_code = "host-capture-failed"
        original_affinity = set(os.sched_getaffinity(0))
        require(original_affinity, "launch affinity is empty")
        lifecycle["affinity"]["captured"] = sorted(original_affinity)
        allowed, housekeeping = validate_topology(options.cpu, options.sibling)
        host_before = host_identity(options.cpu, options.sibling, allowed)
        validate_host(host_before, options.cpu, options.sibling)
        completed_stage_index = 2
        failure_code = "reservation-failed"
        reservation_guard = Reservation(
            options.reservation_file, options.cpu, options.sibling)
        pair_guard = PairLease(
            options.cpu, options.sibling, anchor=reservation_guard.anchor)
        with ExitRecordingGuard(reservation_guard, lifecycle, "reservation") as retained_reservation:
            require(typed_equal(retained_reservation, reservation_guard.identity),
                    "reservation guard returned a mutable or foreign identity")
            validate_reservation_record(retained_reservation, options.cpu, options.sibling)
            reservation_record = copy.deepcopy(retained_reservation)
            completed_stage_index = 3
            failure_code = "lease-failed"
            with ExitRecordingGuard(pair_guard, lifecycle, "pair_lease") as retained_pair:
                require(typed_equal(retained_pair, pair_guard.identity),
                        "pair guard returned a mutable or foreign identity")
                validate_pair_lease_identity(retained_pair, options.cpu, options.sibling)
                pair_lease_record = copy.deepcopy(retained_pair)
                completed_stage_index = 4
                failure_code = "inputs-capture-failed"
                os.sched_setaffinity(0, housekeeping)
                taskset = _lexical_absolute(Path("/usr/bin/taskset"))
                executable = _lexical_absolute(options.executable)
                archive = _lexical_absolute(options.archive)
                build_manifest = _lexical_absolute(options.build_manifest)
                source_root = _lexical_absolute(options.source_root)
                inputs_before = input_snapshot(
                    source_root, options.expected_tooling_commit,
                    options.expected_core_commit, archive, executable, taskset,
                    build_manifest)
                validate_input_snapshot(inputs_before)
                completed_stage_index = 5
                failure_code = "attestation-failed"
                build_provenance = retain_build_provenance(
                    output, build_manifest, inputs_before)
                validate_build_provenance(build_provenance, output, inputs_before)
                completed_stage_index = 6
                failure_code = "child-capture-failed"
                result_path = output / "snapshot/child/result.json"
                stdout_path = output / "snapshot/child/stdout.bin"
                stderr_path = output / "snapshot/child/stderr.bin"
                result_path.parent.mkdir(parents=True, exist_ok=True)
                before_cpu = cpu_stat_snapshot(options.cpu)
                before_sibling = cpu_stat_snapshot(options.sibling)
                command = [str(taskset), "-c", str(options.cpu), str(executable),
                           "--backend", "avx2", str(result_path)]
                completed, timed_out, started, ended = run_child_bounded(
                    command, stdout_path, stderr_path, options.timeout,
                    CHILD_ENVIRONMENT)
                after_cpu = cpu_stat_snapshot(options.cpu)
                after_sibling = cpu_stat_snapshot(options.sibling)
                try:
                    result_record = artifact_record(output, result_path)
                except EvidenceError:
                    result_record = None
                child_record = {
                    "command": request["command"], "environment": dict(CHILD_ENVIRONMENT),
                    "returncode": completed.returncode, "timed_out": timed_out,
                    "duration_ns": ended - started,
                    "stdout": artifact_record(output, stdout_path),
                    "stderr": artifact_record(output, stderr_path),
                    "result": result_record,
                }
                validate_failure_child(child_record, request, output, check_files=True)
                completed_stage_index = 7
                failure_code = "isolation-capture-failed"
                isolation = isolation_record(
                    options.cpu, options.sibling, pair_lease_record,
                    before_cpu, after_cpu, before_sibling, after_sibling,
                    started, ended)
                validate_isolation(isolation, options.cpu, options.sibling,
                                   require_accepted=False)
                completed_stage_index = 8
                failure_code = "child-timeout"
                require(not timed_out, FAILURE_DIAGNOSTICS[failure_code])
                failure_code = "child-exit"
                require(completed.returncode == 0,
                        f"C7 child exited {completed.returncode}")
                failure_code = "child-result-missing"
                require(result_record is not None, FAILURE_DIAGNOSTICS[failure_code])
                stdout = read_artifact(output, child_record["stdout"], "stdout", MAX_LOG_BYTES)
                stderr = read_artifact(output, child_record["stderr"], "stderr", MAX_LOG_BYTES)
                failure_code = "child-stdout-invalid"
                require(stdout == b"", FAILURE_DIAGNOSTICS[failure_code])
                failure_code = "child-stderr-invalid"
                require(stderr == expected_stderr(), FAILURE_DIAGNOSTICS[failure_code])
                failure_code = "isolation-rejected"
                validate_isolation(isolation, options.cpu, options.sibling)
                failure_code = "reservation-post-child-invalid"
                run_guard_check(
                    lifecycle, "reservation_post_child",
                    lambda: reservation_guard.validate_current(reservation_record))
                failure_code = "lease-post-child-invalid"
                run_guard_check(
                    lifecycle, "pair_post_child",
                    lambda: pair_guard.validate_current(pair_lease_record))
                failure_code = "child-result-json-invalid"
                result_bytes = read_artifact(output, result_record, "result")
                parsed = strict_json(result_bytes, "C7 result")
                failure_code = "inputs-after-capture-failed"
                inputs_after = input_snapshot(
                    source_root, options.expected_tooling_commit,
                    options.expected_core_commit, archive, executable, taskset,
                    build_manifest)
                validate_input_snapshot(inputs_after)
                completed_stage_index = 9
                failure_code = "inputs-drift"
                require(typed_equal(inputs_after, inputs_before),
                        FAILURE_DIAGNOSTICS[failure_code])
                failure_code = "host-after-capture-failed"
                host_after = host_identity(options.cpu, options.sibling, allowed)
                validate_host(host_after, options.cpu, options.sibling)
                completed_stage_index = 10
                failure_code = "host-drift"
                require(typed_equal(host_after, host_before),
                        FAILURE_DIAGNOSTICS[failure_code])
                failure_code = "child-result-invalid"
                validated_output = validate_child_result(parsed, options.cpu, inputs_before)
                completed_stage_index = 11

                # This is the final lease call.  Everything it protects is read
                # again after it returns, closing a validate-then-replace window.
                failure_code = "final-snapshot-invalid"
                final_snapshot_predicate = "reservation-guard-validation"
                run_guard_check(
                    lifecycle, "reservation_final",
                    lambda: reservation_guard.validate_current(reservation_record))
                final_snapshot_predicate = "pair-guard-validation"
                run_guard_check(
                    lifecycle, "pair_final",
                    lambda: pair_guard.validate_current(pair_lease_record))
                final_snapshot_predicate = "snapshot-revalidation"
                require(read_artifact(output, result_record, "result") == result_bytes,
                        FAILURE_DIAGNOSTICS[failure_code])
                require(read_artifact(output, child_record["stdout"], "stdout",
                                      MAX_LOG_BYTES) == b"",
                        FAILURE_DIAGNOSTICS[failure_code])
                require(read_artifact(output, child_record["stderr"], "stderr",
                                      MAX_LOG_BYTES) == expected_stderr(),
                        FAILURE_DIAGNOSTICS[failure_code])
                validate_child_result(strict_json(result_bytes, "C7 result"),
                                      options.cpu, inputs_before)
                expected_paths = sorted(
                    [build_provenance["manifest"], build_provenance["runner"],
                     child_record["stdout"], child_record["stderr"], result_record],
                    key=lambda item: item["path"])
                require(typed_equal(artifact_inventory(output), expected_paths),
                        FAILURE_DIAGNOSTICS[failure_code])
    except BaseException as error:
        body_error = error
        body_traceback = traceback.format_exc()

    # Context managers have exited before this point.  Restore and read back the
    # exact launch mask before selecting either terminal outcome.
    _restore_affinity(lifecycle, original_affinity)
    if body_error is None and any(
            lifecycle[name]["exit"] == "failed"
            for name in ("reservation", "pair_lease")):
        body_error = EvidenceError(FAILURE_DIAGNOSTICS["guard-exit-failed"])
        body_traceback = "guard teardown failed after successful body\n"
        failure_code = "guard-exit-failed"
    if body_error is None and original_affinity is not None and not (
            lifecycle["affinity"]["restore_succeeded"] and
            lifecycle["affinity"]["exact"]):
        body_error = EvidenceError(FAILURE_DIAGNOSTICS["affinity-restore-failed"])
        body_traceback = "affinity restore/readback failed after successful body\n"
        failure_code = "affinity-restore-failed"

    context = {
        "arguments": arguments, "request": request, "host_before": host_before,
        "reservation": reservation_record, "pair_lease": pair_lease_record,
        "inputs_before": inputs_before, "build_provenance": build_provenance,
        "child": child_record, "isolation": isolation,
        "inputs_after": inputs_after, "host_after": host_after,
        "validated_output": validated_output,
    }

    if body_error is None:
        try:
            validate_lifecycle(lifecycle, require_success=True)
            created = utc_now()
            success_context = {
                "arguments": arguments, "request": request,
                "inputs_before": inputs_before, "inputs_after": inputs_after,
                "host_before": host_before, "host_after": host_after,
                "reservation": reservation_record, "pair_lease": pair_lease_record,
                "isolation": isolation, "build_provenance": build_provenance,
                "child": child_record, "validated_output": validated_output,
                "lifecycle": lifecycle,
            }
            failure_code = "publication-state-write-failed"
            publication_state = signed(publication_state_payload(
                created, success_context))
            validate_publication_state(publication_state)
            state_path = output / "snapshot/publication-state.json"
            write_json_exclusive(state_path, publication_state)
            publication["state"] = artifact_record(output, state_path)
            failure_code = "raw-write-failed"
            raw = signed({
                "schema": RAW_SCHEMA, "created_utc": created,
                "publication_state": publication["state"], **success_context,
            })
            validate_raw(raw, output, check_files=True)
            raw_path = output / "snapshot/raw.json"
            write_json_exclusive(raw_path, raw)
            publication["raw"] = artifact_record(output, raw_path)
            failure_code = "manifest-validation-failed"
            manifest = signed({
                "schema": MANIFEST_SCHEMA, "created_utc": created, "valid": True,
                "publication_state": publication["state"], "raw": publication["raw"],
                "artifact_inventory": artifact_inventory(output),
                "arguments": arguments, "request": request, "inputs": inputs_before,
                "reservation": reservation_record, "pair_lease": pair_lease_record,
                "isolation": isolation, "build_provenance": build_provenance,
                "lifecycle": lifecycle, "validated_output": validated_output,
            })
            validate_manifest(manifest, output)
            failure_code = "final-snapshot-invalid"
            observed = set(os.sched_getaffinity(0))
            if observed != original_affinity:
                lifecycle["affinity"].update(
                    observed=sorted(observed), restore_succeeded=True,
                    exact=False, error_type=None, error=None)
                raise EvidenceError(FAILURE_DIAGNOSTICS["affinity-restore-failed"])
            require(typed_equal(manifest["artifact_inventory"],
                                artifact_inventory(output)),
                    FAILURE_DIAGNOSTICS[failure_code])
            validate_manifest(manifest, output)
            failure_code = "manifest-write-failed"
            _publish_terminal(
                output, "manifest.json", manifest, output_identity)
        except BaseException as error:
            body_error = error
            body_traceback = traceback.format_exc()
            if (failure_code == "final-snapshot-invalid" and
                    not lifecycle["affinity"]["exact"]):
                failure_code = "affinity-restore-failed"

    if body_error is not None:
        require(FAILURE_CODE_STAGE.get(failure_code) == completed_stage_index,
                "internal failure-code/stage mismatch")
        stage = completed_failure_stage(completed_stage_index)
        snapshot_inventory: dict[str, Any] | None = None
        if failure_code == "final-snapshot-invalid":
            expected_snapshot: list[dict[str, Any]] = []
            if build_provenance is not None:
                expected_snapshot.extend((build_provenance["manifest"],
                                          build_provenance["runner"]))
            if child_record is not None:
                expected_snapshot.extend((child_record["stdout"],
                                          child_record["stderr"]))
                if child_record["result"] is not None:
                    expected_snapshot.append(child_record["result"])
            for retained in (publication["state"], publication["raw"]):
                if retained is not None:
                    expected_snapshot.append(retained)
            snapshot_inventory = {
                "expected": sorted(expected_snapshot,
                                   key=lambda item: item["path"]),
                "observed": artifact_inventory(output),
            }
            snapshot_inventory["predicate"] = (
                "inventory-mismatch" if not typed_equal(
                    snapshot_inventory["expected"],
                    snapshot_inventory["observed"])
                else final_snapshot_predicate)
        state = signed(failure_state_payload(
            stage, failure_code, context, lifecycle, publication,
            snapshot_inventory))
        state_path = output / "failure/state-v2.json"
        write_json_exclusive(state_path, state)
        error_text = _failure_diagnostic(failure_code, child_record)
        failure = signed({
            "schema": FAILURE_SCHEMA, "status": "failed",
            "created_utc": utc_now(), "failure_code": failure_code,
            "completed_stage": stage, "error_type": "EvidenceError",
            "error": error_text, "traceback": _bounded_traceback(body_traceback),
            "state": artifact_record(output, state_path),
            "artifact_inventory": artifact_inventory(output),
            "lifecycle": lifecycle, "publication": publication,
            "snapshot_inventory": snapshot_inventory, **context,
        })
        validate_failure(failure, output, check_files=True)
        _publish_terminal(
            output, "failure.json", failure, output_identity)
        raise EvidenceError(error_text) from body_error

    try:
        print(output / "manifest.json")
    except OSError:
        pass
    return 0


def verify_campaign(options: argparse.Namespace) -> int:
    manifest_path = _lexical_absolute(options.manifest)
    root = manifest_path.parent
    manifest_bytes = read_bounded(manifest_path, MAX_TOP_LEVEL_JSON_BYTES)
    manifest = strict_json(manifest_bytes, "C7 manifest")
    require(manifest_bytes == canonical_bytes(manifest) + b"\n",
            "manifest is not canonical JSON")
    normalized = validate_manifest(manifest, root)
    print(json.dumps({"status": "PASS", "cells": normalized["cell_count"],
                      "manifest_sha256": sha256_bytes(manifest_bytes)},
                     sort_keys=True))
    return 0


def verify_failure_campaign(options: argparse.Namespace) -> int:
    failure_path = _lexical_absolute(options.failure)
    require(failure_path.name == "failure.json",
            "failure replay requires the canonical failure.json filename")
    root = failure_path.parent
    failure_bytes = read_bounded(failure_path, MAX_TOP_LEVEL_JSON_BYTES)
    failure = strict_json(failure_bytes, "C7 failure evidence")
    require(failure_bytes == canonical_bytes(failure) + b"\n",
            "failure evidence is not canonical JSON")
    validated = validate_failure(failure, root, check_files=True)
    print(json.dumps({
        "status": "VERIFIED_FAILURE",
        "failure_code": validated["failure_code"],
        "failure_sha256": sha256_bytes(failure_bytes),
    }, sort_keys=True))
    return 2


def synthetic_result(cpu: int, inputs: Mapping[str, Any]) -> dict[str, Any]:
    cells = []
    for expected in EXPECTED_CELLS:
        k, r, byte_count, batch, losses, exact_field, padded_field = expected
        cell: dict[str, Any] = {
            "K": k, "R": r, "bytes": byte_count, "batch": batch,
            "losses": losses, "exact_field": exact_field,
            "padded_field": padded_field, "exact_coefficients": k * r,
            "exact_decode_terms": k * losses, "padded_encode_scratch": 64,
            "padded_decode_scratch": 128,
        }
        for index, (summary, samples) in enumerate(SUMMARY_SAMPLE_PAIRS, 1):
            values = [float(index)] * SAMPLE_COUNT
            cell[samples] = values
            cell[summary] = {"median_us": float(index), "mad_us": 0.0}
        cells.append(cell)
    return {"schema": CHILD_SCHEMA, "status": "pass",
            "profile": copy.deepcopy(EXPECTED_PROFILE),
            "production_constructor_rejected": True,
            "timing_scope": "candidate-authoritative-requires-runner-manifest",
            "requested_backend": "avx2", "runtime_backend": "avx2",
            "affinity": [cpu], "omp_num_threads": "1", "omp_dynamic": "FALSE",
            "source_sha256": inputs["source"]["sha256"],
            "core_git_sha": inputs["git"]["core_commit"],
            "library_sha256": inputs["archive"]["sha256"], "sanitizer": "none",
            "sanitizer_features": {"address": False, "undefined": False},
            "allocation_tracking": "global-new",
            "correctness": copy.deepcopy(EXPECTED_CORRECTNESS), "benchmarks": cells}


def synthetic_build_manifest(inputs: Mapping[str, Any]) -> dict[str, Any]:
    def artifact(key: str, path: Path) -> dict[str, Any]:
        return {"bytes": inputs[key]["bytes"], "path": path.as_posix(),
                "sha256": inputs[key]["sha256"]}

    build_dir = Path(".research/fixture/build/core-avx2")
    library = {"bytes": inputs["archive"]["bytes"],
               "path": (build_dir / "liblibleopard.a").as_posix(),
               "sha256": inputs["archive"]["sha256"]}
    executable = {"bytes": inputs["executable"]["bytes"],
                  "path": (build_dir.parent / "c7-avx2").as_posix(),
                  "sha256": inputs["executable"]["sha256"]}
    builds: list[dict[str, Any]] = []
    for name in BUILD_NAMES:
        common = {
            "name": name,
            "source_closure": copy.deepcopy(inputs["core_source_closure"]),
            "launcher_python": {
                **copy.deepcopy(inputs["python"]), "version": "fixture-python"},
        }
        if name == "avx2":
            builds.append({
                **common, "backend": "avx2", "sanitizer": False,
                "build_dir": build_dir.as_posix(), "library": library,
                "executable": executable,
                "compile_argv": [
                    "/usr/bin/g++", "-O2", SOURCE_RELATIVE.as_posix(),
                    library["path"], "-pthread", "-fopenmp", "-o",
                    executable["path"],
                ],
                "instrumentation": {
                    "required_compile_macro": False,
                    "executable_counts": {"asan_lines": 0, "ubsan_lines": 0},
                    "core_archive_counts": {"asan_lines": 0, "ubsan_lines": 0},
                },
            })
        else:
            builds.append(common)
    return {
        "schema": BUILD_MANIFEST_SCHEMA, "status": "pass", "scope": BUILD_SCOPE,
        "tooling_git_sha": inputs["git"]["tooling_commit"],
        "core_git_sha": inputs["git"]["core_commit"],
        "source": artifact("source", SOURCE_RELATIVE),
        "runner": artifact("build_runner", BUILD_RUNNER_RELATIVE),
        "validator": artifact("build_validator", BUILD_VALIDATOR_RELATIVE),
        "taskset": {**copy.deepcopy(inputs["taskset"]),
                    "version": "fixture-taskset"},
        "builds": builds,
        "reproducibility": {
            "comparison": {"status": "pass"},
            "fingerprints": {"avx2": {
                "library_sha256": library["sha256"],
                "executable_sha256": executable["sha256"],
            }},
        },
    }


def synthetic_inputs() -> dict[str, Any]:
    synthetic_runner = read_bounded(Path(__file__).resolve(strict=True))

    def artifact(path: str, digest: str = "a" * 64) -> dict[str, Any]:
        return {"bytes": 1, "path": path, "sha256": digest}

    result = {
        "schema": INPUT_SCHEMA,
        "git": {"tooling_commit": "b" * 40, "core_commit": "e" * 40},
        "runner": {
            "bytes": len(synthetic_runner), "path": RUNNER_RELATIVE.as_posix(),
            "sha256": sha256_bytes(synthetic_runner)},
        "source": artifact(SOURCE_RELATIVE.as_posix()),
        "build_runner": artifact(BUILD_RUNNER_RELATIVE.as_posix()),
        "build_validator": artifact(BUILD_VALIDATOR_RELATIVE.as_posix()),
        "archive": artifact(
            ".research/fixture/build/core-avx2/liblibleopard.a"),
        "executable": artifact(".research/fixture/build/c7-avx2"),
        "taskset": {"path": "/usr/bin/taskset", "sha256": "a" * 64},
        "python": {"path": "/usr/bin/python3", "sha256": "a" * 64},
        "core_source_closure": [
            artifact(path)
            for path in CORE_SOURCE_CLOSURE
        ],
    }
    build_data = canonical_pretty_bytes(synthetic_build_manifest(result))
    result["build_manifest"] = {
        "bytes": len(build_data), "schema": BUILD_MANIFEST_SCHEMA,
        "sha256": sha256_bytes(build_data),
    }
    result["build_attestation"] = derive_build_attestation(
        synthetic_build_manifest(result), result["build_manifest"],
        result["git"]["tooling_commit"], result["git"]["core_commit"],
        result["source"], result["build_runner"], result["build_validator"],
        result["archive"], result["executable"],
        result["core_source_closure"], result["taskset"], result["python"])
    result["binding_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def synthetic_cpu_stat(cpu: int, *, user: int, idle: int) -> dict[str, Any]:
    fields = {"user": user, "nice": 0, "system": 0, "idle": idle,
              "iowait": 0, "irq": 0, "softirq": 0, "steal": 0}
    return {"cpu": cpu, "fields": fields, "total_jiffies": sum(fields.values())}


def synthetic_pair_lease(cpu: int, sibling: int) -> dict[str, Any]:
    payload = pair_lease_payload(cpu, sibling, uid=1000)
    authority = {
        "schema": KERNEL_LEASE_SCHEMA,
        "authority": "linux-abstract-unix-bind",
        "namespace": kernel_lease_namespace("pair", payload),
        "device": 1, "inode": 4,
        "payload_sha256": sha256_bytes(canonical_bytes(payload)),
    }
    return {"schema": PAIR_LEASE_SCHEMA, "authority": authority,
            "device": 1, "directory_device": 1, "directory_inode": 2,
            "inode": 3, "lock": "exclusive_nonblocking_pair_wide",
            "path": str(pair_lease_directory(1000) /
                        pair_lease_name(cpu, sibling, 1000)),
            "payload": payload, "sha256": sha256_bytes(canonical_bytes(payload))}


def synthetic_host(cpu: int, sibling: int) -> dict[str, Any]:
    topology = {
        "core_id": "0", "physical_package_id": "0", "die_id": "0",
        "thread_siblings_list": f"{cpu},{sibling}",
        "core_siblings_list": f"{cpu},{sibling},2",
    }
    frequency = {
        "scaling_driver": "fixture", "scaling_governor": "performance",
        "energy_performance_preference": "performance", "scaling_min_freq": "1",
        "scaling_max_freq": "2", "cpuinfo_min_freq": "1",
        "cpuinfo_max_freq": "2",
    }
    def record(index: int) -> dict[str, Any]:
        return {"cpu": index, "online": "1", "topology": copy.deepcopy(topology),
                "frequency_policy": copy.deepcopy(frequency)}
    return {
        "system": {"system": "Linux", "release": "fixture",
                   "version": "fixture", "machine": "x86_64", "python": "3.11"},
        "allowed_at_launch": sorted((cpu, sibling, 2)),
        "timing_cpu": record(cpu), "sibling_cpu": record(sibling),
        "turbo": {
            "/sys/devices/system/cpu/intel_pstate/no_turbo": "0",
            "/sys/devices/system/cpu/amd_pstate/status": None,
            "/sys/devices/system/cpu/cpufreq/boost": None,
        },
    }


def synthetic_bundle(root: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    cpu, sibling = 0, 1
    inputs = synthetic_inputs()
    result = synthetic_result(cpu, inputs)
    result_path = root / "snapshot/child/result.json"
    stdout_path = root / "snapshot/child/stdout.bin"
    stderr_path = root / "snapshot/child/stderr.bin"
    write_json_exclusive(result_path, result)
    write_exclusive(stdout_path, b"")
    write_exclusive(stderr_path, expected_stderr())
    retained_build_manifest = root / "snapshot/provenance/build-run-manifest-v4.json"
    write_exclusive(retained_build_manifest,
                    canonical_pretty_bytes(synthetic_build_manifest(inputs)))
    retained_runner = root / "snapshot/provenance/run_authoritative.py"
    write_exclusive(retained_runner, read_bounded(Path(__file__).resolve(strict=True)))
    build_provenance = {
        "schema": BUILD_PROVENANCE_SCHEMA,
        "manifest": artifact_record(root, retained_build_manifest),
        "runner": artifact_record(root, retained_runner),
        "attestation": copy.deepcopy(inputs["build_attestation"]),
    }
    payload = {"benchmark_cpu": cpu, "nonce": "fixture-nonce", "owner": "self-test",
               "reserved_sibling": sibling, "schema": RESERVATION_SCHEMA,
               "status": "held"}
    reservation_bytes = canonical_bytes(payload)
    reservation_sha = sha256_bytes(reservation_bytes)
    reservation_path = Path("/fixture/cpu-reservation.json")
    path_authority_payload = reservation_authority_payload(
        reservation_path, cpu, sibling, uid=1000)
    content_authority_payload = reservation_authority_payload(
        reservation_path, cpu, sibling, reservation_sha, uid=1000)
    def authority(identity: int, authority_payload: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "schema": KERNEL_LEASE_SCHEMA,
            "authority": "linux-abstract-unix-bind",
            "namespace": kernel_lease_namespace("reservation", authority_payload),
            "device": 1, "inode": identity,
            "payload_sha256": sha256_bytes(canonical_bytes(authority_payload)),
        }
    reservation = {
        "schema": RESERVATION_SCHEMA, "bytes": len(reservation_bytes),
        "sha256": reservation_sha, "payload": payload,
        "lock": "linux_abstract_bind_plus_flock", "path": str(reservation_path),
        "uid": 1000, "device": 1, "inode": 5,
        "directory_device": 1, "directory_inode": 6,
        "authority": {"path": authority(7, path_authority_payload),
                      "content": authority(8, content_authority_payload)},
    }
    isolation = isolation_record(
        cpu, sibling, synthetic_pair_lease(cpu, sibling),
        synthetic_cpu_stat(cpu, user=0, idle=0),
        synthetic_cpu_stat(cpu, user=1, idle=0),
        synthetic_cpu_stat(sibling, user=0, idle=0),
        synthetic_cpu_stat(sibling, user=0, idle=1), 1, 2)
    request = {
        "backend": "avx2", "cell_geometry": [list(cell) for cell in EXPECTED_CELLS],
        "child_environment": copy.deepcopy(CHILD_ENVIRONMENT), "cpu": cpu,
        "sibling": sibling, "timeout_seconds": 10.0,
        "command": ["${TASKSET}", "-c", "0", "${C7_EXECUTABLE}",
                    "--backend", "avx2", "${RESULT_JSON}"],
    }
    normalized = validate_child_result(result, cpu, inputs)
    created = "2026-07-17T00:00:00Z"
    arguments = {
        "cpu": cpu, "sibling": sibling, "timeout_seconds": 10.0,
        "expected_tooling_commit": inputs["git"]["tooling_commit"],
        "expected_core_commit": inputs["git"]["core_commit"],
    }
    lifecycle = empty_lifecycle()
    for name in ("reservation", "pair_lease"):
        lifecycle[name].update(entered=True, exit="pass")
    for check in lifecycle["guard_checks"].values():
        check["status"] = "pass"
    lifecycle["affinity"].update(
        captured=[0, 1, 2], restore_attempted=True,
        restore_succeeded=True, observed=[0, 1, 2], exact=True)
    child = {"command": request["command"],
             "environment": copy.deepcopy(CHILD_ENVIRONMENT),
             "returncode": 0, "timed_out": False, "duration_ns": 1,
             "stdout": artifact_record(root, stdout_path),
             "stderr": artifact_record(root, stderr_path),
             "result": artifact_record(root, result_path)}
    success_context = {
        "arguments": arguments, "request": request,
        "inputs_before": inputs, "inputs_after": copy.deepcopy(inputs),
        "host_before": synthetic_host(cpu, sibling),
        "host_after": synthetic_host(cpu, sibling), "reservation": reservation,
        "pair_lease": isolation["pair_lease"], "isolation": isolation,
        "build_provenance": build_provenance, "child": child,
        "validated_output": normalized, "lifecycle": lifecycle,
    }
    state = signed(publication_state_payload(created, success_context))
    state_path = root / "snapshot/publication-state.json"
    write_json_exclusive(state_path, state)
    raw = signed({
        "schema": RAW_SCHEMA, "created_utc": created,
        "publication_state": artifact_record(root, state_path),
        **success_context,
    })
    raw_path = root / "snapshot/raw.json"
    write_json_exclusive(raw_path, raw)
    manifest = signed({
        "schema": MANIFEST_SCHEMA, "created_utc": created, "valid": True,
        "raw": artifact_record(root, raw_path),
        "publication_state": artifact_record(root, state_path),
        "artifact_inventory": artifact_inventory(root),
        "arguments": arguments, "request": request,
        "inputs": inputs, "reservation": reservation,
        "pair_lease": isolation["pair_lease"], "isolation": isolation,
        "build_provenance": build_provenance, "lifecycle": lifecycle,
        "validated_output": normalized,
    })
    return manifest, raw, result


def self_test() -> int:
    inputs = synthetic_inputs()
    validate_input_snapshot(inputs)
    result = synthetic_result(0, inputs)
    validate_child_result(result, 0, inputs)
    rejected = 0

    def expect_rejected(function: Any, label: str) -> None:
        nonlocal rejected
        try:
            function()
        except EvidenceError:
            rejected += 1
        else:
            raise EvidenceError(f"{label} mutation was accepted")

    for label, mutate in (
        ("summary", lambda value: value["benchmarks"][0]["exact_encode"].update(
            median_us=99.0)),
        ("sample count", lambda value: value["benchmarks"][0][
            "exact_encode_samples_us"].pop()),
        ("cell order", lambda value: value["benchmarks"].reverse()),
        ("backend", lambda value: value.update(runtime_backend="scalar")),
        ("embedded source", lambda value: value.update(source_sha256="f" * 64)),
        ("Boolean scratch", lambda value: value["benchmarks"][0].update(
            padded_encode_scratch=False)),
    ):
        changed = copy.deepcopy(result)
        mutate(changed)
        expect_rejected(lambda changed=changed: validate_child_result(
            changed, 0, inputs), label)

    stale = {"benchmark_cpu": 0, "nonce": "fixture-nonce", "owner": "self-test",
             "reserved_sibling": 1, "schema": RESERVATION_SCHEMA,
             "status": "released"}
    expect_rejected(lambda: parse_reservation(canonical_bytes(stale), 0, 1),
                    "stale reservation")

    busy = isolation_record(
        0, 1, synthetic_pair_lease(0, 1),
        synthetic_cpu_stat(0, user=0, idle=0),
        synthetic_cpu_stat(0, user=1, idle=0),
        synthetic_cpu_stat(1, user=0, idle=0),
        synthetic_cpu_stat(1, user=1, idle=1), 1, 2)
    expect_rejected(lambda: validate_isolation(busy, 0, 1), "busy sibling")

    with tempfile.TemporaryDirectory(prefix="leopard2-c7-authoritative-") as directory:
        root = Path(directory)
        manifest, raw, _fixture_result = synthetic_bundle(root)
        validate_manifest(manifest, root)
        changed_raw = copy.deepcopy(raw)
        changed_raw["request"]["cell_geometry"].reverse()
        changed_raw = signed({key: value for key, value in changed_raw.items()
                              if key != "digest"})
        expect_rejected(lambda: validate_raw(changed_raw, root), "raw schedule")
        changed_raw = copy.deepcopy(raw)
        changed_raw["host_after"]["timing_cpu"]["frequency_policy"][
            "scaling_governor"] = "powersave"
        changed_raw = signed({key: value for key, value in changed_raw.items()
                              if key != "digest"})
        expect_rejected(lambda: validate_raw(changed_raw, root), "host drift")
        result_path = root / "snapshot/child/result.json"
        original = result_path.read_bytes()
        result_path.write_bytes(b"{}\n")
        expect_rejected(lambda: validate_manifest(manifest, root), "artifact bytes")
        result_path.write_bytes(original)
        changed_manifest = copy.deepcopy(manifest)
        changed_manifest["validated_output"]["cell_count"] = 11
        changed_manifest = signed({key: value for key, value in changed_manifest.items()
                                   if key != "digest"})
        expect_rejected(lambda: validate_manifest(
            changed_manifest, root), "manifest summary")

    print(json.dumps({"status": "PASS", "cells": len(EXPECTED_CELLS),
                      "samples_per_metric": SAMPLE_COUNT,
                      "mutations_rejected": rejected}, sort_keys=True))
    return 0


class EvidenceArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise EvidenceError(f"invalid command line: {message}")


def parser() -> argparse.ArgumentParser:
    result = EvidenceArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    run = commands.add_parser(
        "run", help="execute one fixed 12-cell authoritative AVX2 campaign")
    run.add_argument("--source-root", required=True, type=Path,
                     help="clean checkout containing the committed C7 tools")
    run.add_argument("--expected-tooling-commit", required=True,
                     help="full SHA-1 that must equal clean checkout HEAD")
    run.add_argument("--expected-core-commit", required=True,
                     help="full ancestor SHA-1 embedded in the C7 executable")
    run.add_argument("--archive", required=True, type=Path,
                     help="exact AVX2 liblibleopard.a linked into the harness")
    run.add_argument("--executable", required=True, type=Path,
                     help="exact non-sanitized AVX2 c7_exact_low executable")
    run.add_argument("--build-manifest", required=True, type=Path,
                     help=("final v4 A/B build manifest validated live and bound "
                           "to the supplied AVX2 archive/executable"))
    run.add_argument("--reservation-file", required=True, type=Path,
                     help="canonical held v1 CPU-pair reservation to lock")
    run.add_argument("--output", required=True, type=Path,
                     help="new directory for raw, manifest, logs, or failure")
    run.add_argument("--cpu", required=True, type=int,
                     help="allowed logical CPU used for all timing")
    run.add_argument("--sibling", required=True, type=int,
                     help="the timing CPU's allowed SMT sibling to reserve")
    run.add_argument("--timeout", type=float, default=1800.0,
                     help="positive child timeout in seconds (default: 1800)")
    run.set_defaults(function=run_campaign)
    verify = commands.add_parser(
        "verify", help="portably replay a retained canonical evidence bundle")
    verify.add_argument("--manifest", required=True, type=Path,
                        help="retained manifest.json")
    verify.set_defaults(function=verify_campaign)
    verify_failure = commands.add_parser(
        "verify-failure",
        help="portably authenticate a failed campaign without treating it as success")
    verify_failure.add_argument("--failure", required=True, type=Path,
                                help="retained failure.json")
    verify_failure.set_defaults(function=verify_failure_campaign)
    test = commands.add_parser(
        "self-test", help="run synthetic validation and mutation tests only")
    test.set_defaults(function=lambda unused: self_test())
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    try:
        options = parser().parse_args(arguments)
        return int(options.function(options))
    except (EvidenceError, OSError, ValueError, OverflowError, RecursionError,
            AttributeError, IndexError, KeyError, TypeError,
            subprocess.SubprocessError) as error:
        print(f"C7 authoritative evidence error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
