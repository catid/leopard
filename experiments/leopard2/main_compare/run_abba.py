#!/usr/bin/env python3
"""Fail-closed Leopard2 versus exact Leopard-main ABBA comparison runner.

The runner deliberately does not build either implementation.  It executes two
already-built, independently linked benchmark processes on one pinned CPU,
retains their byte-for-byte output, and refuses to analyze evidence unless the
workload and build identities match the signed request.
"""

from __future__ import annotations

import argparse
import base64
import binascii
import copy
import ctypes
import datetime as dt
import errno
import fcntl
import hashlib
import json
import math
import os
import platform
import re
import selectors
import shlex
import signal
import socket
import stat
import statistics
import subprocess
import sys
import time
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

_DECODER_DISPATCH_DIR = Path(__file__).resolve().parents[1] / "decoder_dispatch"
if str(_DECODER_DISPATCH_DIR) not in sys.path:
    sys.path.insert(0, str(_DECODER_DISPATCH_DIR))
import balanced_evidence_common as link_common


MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
RAW_SCHEMA_V1 = "leopard2-main-compare-raw/v1"
RAW_SCHEMA_V2 = "leopard2-main-compare-raw/v2"
RAW_SCHEMA_V3 = "leopard2-main-compare-raw/v3"
RAW_SCHEMA_V4 = "leopard2-main-compare-raw/v4"
RAW_SCHEMA_V5 = "leopard2-main-compare-raw/v5"
RAW_SCHEMA = "leopard2-main-compare-raw/v6"
HARDENED_HISTORICAL_BUILD_SCHEMA = \
    "leopard2-main-compare-build/hardened-historical-v1"
MANIFEST_SCHEMA_V1 = "leopard2-main-compare-manifest/v1"
MANIFEST_SCHEMA_V2 = "leopard2-main-compare-manifest/v2"
MANIFEST_SCHEMA_V3 = "leopard2-main-compare-manifest/v3"
MANIFEST_SCHEMA_V4 = "leopard2-main-compare-manifest/v4"
MANIFEST_SCHEMA_V5 = "leopard2-main-compare-manifest/v5"
MANIFEST_SCHEMA = "leopard2-main-compare-manifest/v6"
FAILURE_SCHEMA_V2 = "leopard2-main-compare-failure/v2"
FAILURE_SCHEMA_V3 = "leopard2-main-compare-failure/v3"
FAILURE_SCHEMA_V4 = "leopard2-main-compare-failure/v4"
FAILURE_SCHEMA_V5 = "leopard2-main-compare-failure/v5"
FAILURE_SCHEMA = "leopard2-main-compare-failure/v6"
RESERVATION_SCHEMA = "leopard2-cpu-reservation/v1"
PAIR_LEASE_SCHEMA = "leopard2-cpu-pair-lease/v1"
ISOLATION_SCHEMA = "leopard2-main-compare-isolation/v1"
SUPERVISION_SCHEMA = "leopard2-main-supervision/v1"
SUPERVISION_NONCE_ENV = "LEO2_AFFINITY_EXECUTION_NONCE"
CANONICAL_LDD_SCHEMA = "leopard2-main-compare-canonical-ldd/v1"
CANONICAL_LDD_NORMALIZATION = "terminal-aslr-load-address/v1"
CANONICAL_LDD_ADDRESS = "<ASLR_LOAD_ADDRESS>"

# CMake target and archive identity is evidence, not an interchangeable build
# detail.  Historical v1/v2 records predate the canonical target rename and
# must continue to replay against their exact old names.  Version 3 and later bind
# the canonical target/archive. Verification selects one exact identity from
# the signed schema; runtime backend selection cannot alter it.
HISTORICAL_CMAKE_IDENTITY = {
    "target": "libleopard",
    "archive": "liblibleopard.a",
    "target_directory": "libleopard.dir",
}
CANONICAL_CMAKE_IDENTITY = {
    "target": "leopard",
    "archive": "libleopard.a",
    "target_directory": "leopard.dir",
}
RAW_TO_CMAKE_IDENTITY = {
    RAW_SCHEMA_V1: HISTORICAL_CMAKE_IDENTITY,
    RAW_SCHEMA_V2: HISTORICAL_CMAKE_IDENTITY,
    RAW_SCHEMA_V3: CANONICAL_CMAKE_IDENTITY,
    RAW_SCHEMA_V4: CANONICAL_CMAKE_IDENTITY,
    RAW_SCHEMA_V5: CANONICAL_CMAKE_IDENTITY,
    RAW_SCHEMA: CANONICAL_CMAKE_IDENTITY,
}
# This internal build-only schema lets another evidence family authenticate an
# exact pre-rename CMake graph with current recipe-content hardening.  It is
# deliberately absent from RAW_TO_CMAKE_IDENTITY and every manifest/failure
# version map, so main-comparison evidence can never claim this private schema.
BUILD_SCHEMA_TO_CMAKE_IDENTITY = {
    **RAW_TO_CMAKE_IDENTITY,
    HARDENED_HISTORICAL_BUILD_SCHEMA: HISTORICAL_CMAKE_IDENTITY,
}
HARDENED_BUILD_SCHEMAS = frozenset((
    RAW_SCHEMA_V3,
    RAW_SCHEMA_V4,
    RAW_SCHEMA_V5,
    RAW_SCHEMA,
    HARDENED_HISTORICAL_BUILD_SCHEMA,
))
MANIFEST_TO_RAW_SCHEMA = {
    MANIFEST_SCHEMA_V1: RAW_SCHEMA_V1,
    MANIFEST_SCHEMA_V2: RAW_SCHEMA_V2,
    MANIFEST_SCHEMA_V3: RAW_SCHEMA_V3,
    MANIFEST_SCHEMA_V4: RAW_SCHEMA_V4,
    MANIFEST_SCHEMA_V5: RAW_SCHEMA_V5,
    MANIFEST_SCHEMA: RAW_SCHEMA,
}
FAILURE_TO_RAW_SCHEMA = {
    FAILURE_SCHEMA_V2: RAW_SCHEMA_V2,
    FAILURE_SCHEMA_V3: RAW_SCHEMA_V3,
    FAILURE_SCHEMA_V4: RAW_SCHEMA_V4,
    FAILURE_SCHEMA_V5: RAW_SCHEMA_V5,
    FAILURE_SCHEMA: RAW_SCHEMA,
}
CANDIDATE_MODE_SCHEMAS = frozenset((RAW_SCHEMA_V4, RAW_SCHEMA_V5, RAW_SCHEMA))
WORKSPACE_SELECTOR_SCHEMAS = frozenset((
    RAW_SCHEMA_V3, RAW_SCHEMA_V4, RAW_SCHEMA_V5, RAW_SCHEMA,
))
ISOLATION_SCHEMAS = frozenset((
    RAW_SCHEMA_V2, RAW_SCHEMA_V3, RAW_SCHEMA_V4, RAW_SCHEMA_V5, RAW_SCHEMA,
))
COMPLETE_EVIDENCE_SCHEMAS = frozenset((RAW_SCHEMA_V5, RAW_SCHEMA))
SUPERVISION_SCHEMAS = COMPLETE_EVIDENCE_SCHEMAS
CPU_STAT_FIELDS = (
    "user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal",
)
CPU_STAT_IDLE_FIELDS = ("idle", "iowait")
MAX_SIBLING_NONIDLE_JIFFIES = 0
MAX_CPU_ID = 1_048_575
MAX_CPU_LIST_ENTRIES = 4096
MAX_CPU_LIST_TEXT_BYTES = 65_536
ROUNDS = 3
ORDER = ("baseline", "candidate", "candidate", "baseline")
CANDIDATE_MODES = ("auto", "generic", "materialized", "tiled")
FNV_OFFSET = 14695981039346656037
FNV_PRIME = 1099511628211
MASK64 = (1 << 64) - 1
HEX64 = re.compile(r"^[0-9a-f]{16}$")
HEX256 = re.compile(r"^[0-9a-f]{64}$")
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
MAX_COMMAND_STDOUT_BYTES = 128 * 1024 * 1024
MAX_COMMAND_STDERR_BYTES = 8 * 1024 * 1024
MAX_COMMAND_TIMEOUT_SECONDS = 3600.0
MAX_IDENTITY_FILE_BYTES = 256 * 1024 * 1024
MAX_LINK_RECIPE_BYTES = 1024 * 1024
MAX_GIT_COMMIT_BYTES = 1024 * 1024
MAX_CAMPAIGN_COUNT = 1_000_000
MAX_CAMPAIGN_CELLS = 4096
MAX_SHARD_BYTES = 1 << 40
BASELINE_LIBRARY_SOURCES = (
    "leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp",
)
CANDIDATE_LIBRARY_SOURCES = (
    "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp", "Leopard2Plan.cpp",
    "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp",
    "Leopard2BackendSSSE3.cpp", "Leopard2BackendAVX2.cpp",
    "Leopard2BackendAVX2Xor.cpp",
    "Leopard2BackendAVX512.cpp", "Leopard2BackendGFNI.cpp",
)
CANDIDATE_NON_LIBRARY_COMPILE_TARGETS = {
    "tests/benchmark.cpp": "bench_leopard.dir",
    "bench/leopard2/benchmark.cpp": "bench_leopard2.dir",
    "bench/leopard2/batch_preflight_benchmark.cpp":
        "bench_leopard2_batch_preflight.dir",
    "bench/leopard2/locator_benchmark.cpp": "bench_leopard2_locator.dir",
    "bench/leopard2/sparse_encode_benchmark.cpp":
        "leopard2_sparse_encode_benchmark_object.dir",
    "tests/leopard2/direct_oracle.cpp":
        "leopard2_sparse_encode_oracle_object.dir",
    "experiments/leopard2/gf16_high_encode/attribution.cpp":
        "bench_leopard2_gf16_high_attribution.dir",
    "tests/experiments.cpp": "experiment_leopard.dir",
}
CANDIDATE_CONFIGURED_SOURCES = (
    *CANDIDATE_LIBRARY_SOURCES, *CANDIDATE_NON_LIBRARY_COMPILE_TARGETS,
)
BASELINE_EXPECTED_COMPILE_COMMAND_COUNT = 5
CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT = len(CANDIDATE_CONFIGURED_SOURCES)
COMPILE_COMMANDS_SCHEMA = "leopard2-main-compare-compile-commands/v2"
BASELINE_COMPILE_PROFILE = \
    "gnu-compatible-cxx11-native-x86_64-release/v1"
CANDIDATE_COMPILE_PROFILE = \
    "gnu-compatible-cxx11-runtime-dispatch-x86_64-release/v1"
CHILD_REAP_TIMEOUT_SECONDS = 5.0
PR_SET_CHILD_SUBREAPER = 36
PR_GET_CHILD_SUBREAPER = 37


class EvidenceError(RuntimeError):
    """The requested run or retained evidence is not authoritative."""


@dataclass(frozen=True)
class Cell:
    identifier: str
    k: int
    r: int
    shard_bytes: int
    losses: int
    seed: int


REPRESENTATIVE_CELLS = (
    Cell("xor-gf8", 129, 1, 65536, 1, 101),
    Cell("gf8-high-one", 240, 16, 65536, 1, 103),
    Cell("gf8-high-full", 240, 16, 65536, 16, 107),
    Cell("gf8-balanced-full", 128, 128, 65536, 128, 109),
    Cell("gf16-inflation-eight", 200, 50, 65536, 8, 113),
    Cell("gf16-high-one", 1000, 200, 65536, 1, 127),
    Cell("gf16-high-full", 1000, 200, 65536, 200, 131),
    Cell("gf16-large-eight", 4096, 512, 4096, 8, 137),
)
SMOKE_CELLS = (Cell("smoke", 8, 4, 64, 1, 1),)


def statistics_policy() -> dict[str, Any]:
    return {
        "method": "one mean log contrast per independent ABBA round",
        "confidence": 0.95,
        "critical_distribution": "Student-t",
        "independent_round_count_per_metric": ROUNDS,
        "constituent_pair_count_per_metric": 2 * ROUNDS,
        "degrees_of_freedom": ROUNDS - 1,
        "child_estimator": "median of retained per-invocation samples",
        "decode_first_use_semantics": (
            "derived median plan-create plus median one-execution with codec already "
            "created; separate timing loops; excludes codec setup"),
        "decode_reuse_amortized_semantics": (
            "derived median execution plus median plan-create divided by reuse; "
            "separate timing loops; excludes codec setup"),
    }


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise EvidenceError(f"value is not canonical JSON: {error}") from error


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def bounded_file_snapshot(
    path: Path, limit: int = MAX_IDENTITY_FILE_BYTES
) -> tuple[os.stat_result, str, bytes]:
    require(type(limit) is int and 0 <= limit <= MAX_IDENTITY_FILE_BYTES,
            "file identity limit is invalid")
    resolved = path.resolve(strict=True)
    before = os.lstat(resolved)
    require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
            0 <= before.st_size <= limit,
            f"file identity is not a bounded single-link regular file: {resolved}")
    descriptor = os.open(
        resolved, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))
    try:
        initial = os.fstat(descriptor)
        path_metadata = os.lstat(resolved)
        require(stat.S_ISREG(initial.st_mode) and initial.st_nlink == 1 and
                (initial.st_dev, initial.st_ino) ==
                (before.st_dev, before.st_ino) ==
                (path_metadata.st_dev, path_metadata.st_ino) and
                0 <= initial.st_size <= limit,
                f"file identity is not a bounded inode-bound single-link regular file: "
                f"{resolved}")
        digest = hashlib.sha256()
        prefix = b""
        retained = 0
        offset = 0
        while True:
            block = os.pread(
                descriptor, min(1024 * 1024, limit + 1 - retained), offset)
            if not block:
                break
            if not prefix:
                prefix = block[:8]
            digest.update(block)
            retained += len(block)
            offset += len(block)
            require(retained <= limit, f"file identity exceeds {limit} bytes")
        final = os.fstat(descriptor)
        final_path = os.lstat(resolved)
        require(retained == initial.st_size and final.st_nlink == 1 and
                stat.S_ISREG(final_path.st_mode) and final_path.st_nlink == 1 and
                (final.st_dev, final.st_ino) ==
                (initial.st_dev, initial.st_ino) ==
                (final_path.st_dev, final_path.st_ino) and
                (final.st_size, final.st_mtime_ns, final.st_ctime_ns,
                 final.st_dev, final.st_ino) ==
                (initial.st_size, initial.st_mtime_ns, initial.st_ctime_ns,
                 initial.st_dev, initial.st_ino),
                f"file identity changed while hashing: {resolved}")
        return initial, digest.hexdigest(), prefix
    finally:
        os.close(descriptor)


def sha256_file(path: Path, limit: int = MAX_IDENTITY_FILE_BYTES) -> str:
    return bounded_file_snapshot(path, limit)[1]


def bounded_file_contents_snapshot(
    path: Path, limit: int = MAX_IDENTITY_FILE_BYTES
) -> tuple[os.stat_result, bytes]:
    """Read one inode-bound, bounded snapshot and return its exact bytes."""
    require(type(limit) is int and 0 <= limit <= MAX_IDENTITY_FILE_BYTES,
            "file-content snapshot limit is invalid")
    resolved = path.resolve(strict=True)
    before = os.lstat(resolved)
    require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
            0 <= before.st_size <= limit,
            f"file-content snapshot is not a bounded single-link regular file: "
            f"{resolved}")
    descriptor = os.open(
        resolved, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))
    try:
        initial = os.fstat(descriptor)
        path_metadata = os.lstat(resolved)
        require(stat.S_ISREG(initial.st_mode) and initial.st_nlink == 1 and
                (initial.st_dev, initial.st_ino) ==
                (before.st_dev, before.st_ino) ==
                (path_metadata.st_dev, path_metadata.st_ino) and
                0 <= initial.st_size <= limit,
                "file-content snapshot is not an inode-bound regular file: "
                f"{resolved}")
        blocks: list[bytes] = []
        retained = 0
        while retained < initial.st_size:
            block = os.pread(
                descriptor, min(1024 * 1024, initial.st_size - retained),
                retained)
            require(block, f"file-content snapshot was truncated: {resolved}")
            blocks.append(block)
            retained += len(block)
        require(not os.pread(descriptor, 1, retained),
                f"file-content snapshot grew while reading: {resolved}")
        final = os.fstat(descriptor)
        final_path = os.lstat(resolved)
        require(retained == initial.st_size and final.st_nlink == 1 and
                stat.S_ISREG(final_path.st_mode) and final_path.st_nlink == 1 and
                (final.st_dev, final.st_ino) ==
                (initial.st_dev, initial.st_ino) ==
                (final_path.st_dev, final_path.st_ino) and
                (final.st_size, final.st_mtime_ns, final.st_ctime_ns,
                 final.st_dev, final.st_ino) ==
                (initial.st_size, initial.st_mtime_ns, initial.st_ctime_ns,
                 initial.st_dev, initial.st_ino),
                f"file-content snapshot changed while reading: {resolved}")
        return initial, b"".join(blocks)
    finally:
        os.close(descriptor)


def json_file_snapshot(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    """Parse JSON from exactly the bytes accepted by the bounded snapshot."""
    _, data = bounded_file_contents_snapshot(path)
    try:
        value = json.loads(data.decode("utf-8", errors="strict"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"{label} JSON is invalid: {error}") from error
    require(isinstance(value, dict), f"{label} is not an object")
    return value, data


def signed(value: Mapping[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(dict(value))
    require("digest" not in result, "digest field already exists")
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def verify_signature(value: object, what: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{what} is not an object")
    payload = copy.deepcopy(value)
    digest = payload.pop("digest", None)
    require(isinstance(digest, str) and HEX256.fullmatch(digest) is not None,
            f"{what} has no canonical SHA-256 digest")
    require(sha256_bytes(canonical_bytes(payload)) == digest,
            f"{what} digest mismatch")
    return value


def write_json_exclusive(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = canonical_bytes(value) + b"\n"
    try:
        with path.open("xb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise EvidenceError(f"refusing to replace evidence file {path}") from error


def write_bytes_exclusive(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("xb") as stream:
            stream.write(value)
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise EvidenceError(f"refusing to replace evidence file {path}") from error


def _linux_prctl(option: int, argument: object) -> None:
    """Invoke the Linux prctl needed for temporary child-subreaper ownership."""
    require(sys.platform.startswith("linux"),
            "child descendant containment requires Linux")
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        prctl = libc.prctl
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            f"Linux child-subreaper prctl is unavailable: {error}") from error
    ctypes.set_errno(0)
    result = prctl(ctypes.c_int(option), argument,
                   ctypes.c_ulong(0), ctypes.c_ulong(0), ctypes.c_ulong(0))
    if result != 0:
        error_number = ctypes.get_errno()
        raise EvidenceError(
            "Linux child-subreaper prctl failed: " +
            os.strerror(error_number or 1))


def _get_child_subreaper() -> int:
    value = ctypes.c_int(-1)
    _linux_prctl(PR_GET_CHILD_SUBREAPER, ctypes.byref(value))
    require(value.value in (0, 1),
            "Linux returned an invalid child-subreaper state")
    return value.value


def _set_child_subreaper(value: int) -> None:
    require(value in (0, 1), "invalid child-subreaper state")
    _linux_prctl(PR_SET_CHILD_SUBREAPER, ctypes.c_ulong(value))
    require(_get_child_subreaper() == value,
            "Linux did not apply the requested child-subreaper state")


def _linux_pidfd_open(pid: int) -> int | None:
    """Open a race-free Linux process handle or report that the PID is gone."""
    try:
        pidfd_open = ctypes.CDLL(None, use_errno=True).pidfd_open
    except (AttributeError, OSError) as error:
        raise EvidenceError(f"Linux pidfd_open is unavailable: {error}") from error
    ctypes.set_errno(0)
    descriptor = pidfd_open(ctypes.c_int(pid), ctypes.c_uint(0))
    if descriptor >= 0:
        return descriptor
    error_number = ctypes.get_errno()
    if error_number == errno.ESRCH:
        return None
    raise EvidenceError(
        f"cannot open Linux pidfd for process {pid}: " +
        os.strerror(error_number or errno.EPERM))


def _linux_pidfd_signal(descriptor: int, signal_number: int) -> None:
    """Signal the exact process referenced by a pidfd, never a reused PID."""
    try:
        send_signal = ctypes.CDLL(None, use_errno=True).pidfd_send_signal
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            f"Linux pidfd_send_signal is unavailable: {error}") from error
    ctypes.set_errno(0)
    result = send_signal(
        ctypes.c_int(descriptor), ctypes.c_int(signal_number), None,
        ctypes.c_uint(0))
    if result != 0:
        error_number = ctypes.get_errno()
        if error_number == errno.ESRCH:
            return
        raise EvidenceError(
            "cannot signal contained Linux process through pidfd: " +
            os.strerror(error_number or errno.EPERM))


def _validate_linux_pidfd_support() -> None:
    descriptor = _linux_pidfd_open(os.getpid())
    require(descriptor is not None,
            "Linux pidfd support cannot identify the runner process")
    try:
        _linux_pidfd_signal(descriptor, 0)
    finally:
        os.close(descriptor)


def _proc_process_record(pid: int) -> tuple[int, int, int, int, str] | None:
    """Return (ppid, pgrp, session, starttime, state) from Linux procfs."""
    try:
        data = (Path("/proc") / str(pid) / "stat").read_bytes()
    except (FileNotFoundError, ProcessLookupError):
        return None
    except OSError as error:
        if error.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise EvidenceError(f"cannot inspect Linux process {pid}: {error}") from error
    closing = data.rfind(b")")
    require(closing > 0 and closing + 2 < len(data),
            f"Linux process {pid} has malformed procfs stat data")
    fields = data[closing + 2:].split()
    require(len(fields) >= 20,
            f"Linux process {pid} has truncated procfs stat data")
    try:
        state = fields[0].decode("ascii")
        ppid = int(fields[1])
        pgrp = int(fields[2])
        session = int(fields[3])
        starttime = int(fields[19])
    except (UnicodeDecodeError, ValueError) as error:
        raise EvidenceError(
            f"Linux process {pid} has invalid procfs stat fields") from error
    require(len(state) == 1 and ppid >= 0 and pgrp >= 0 and
            session >= 0 and starttime >= 0,
            f"Linux process {pid} has invalid procfs process identity")
    return ppid, pgrp, session, starttime, state


def _proc_process_snapshot() -> dict[int, tuple[int, int, int, int, str]]:
    """Snapshot all procfs-visible processes, including every same-UID child."""
    proc = Path("/proc")
    require(proc.is_dir() and (proc / "self/stat").is_file(),
            "child descendant containment requires mounted Linux procfs")
    try:
        names = os.listdir(proc)
    except OSError as error:
        raise EvidenceError(f"cannot enumerate Linux procfs: {error}") from error
    result: dict[int, tuple[int, int, int, int, str]] = {}
    for name in names:
        if not name.isascii() or not name.isdigit():
            continue
        pid = int(name)
        try:
            record = _proc_process_record(pid)
        except EvidenceError:
            # hidepid may make unrelated processes unreadable.  Every same-UID
            # child of this runner must remain inspectable or the run fails.
            try:
                owner = (proc / name).stat().st_uid
            except OSError:
                continue
            if owner == os.getuid():
                raise
            continue
        if record is not None:
            result[pid] = record
    self_record = _proc_process_record(os.getpid())
    require(self_record is not None and os.getpid() in result,
            "Linux procfs does not expose the runner process")
    return result


def _process_identity(
    pid: int, snapshot: Mapping[int, tuple[int, int, int, int, str]]
) -> tuple[int, int] | None:
    record = snapshot.get(pid)
    return None if record is None else (pid, record[3])


def _emergency_linux_prctl(option: int, argument: object) -> None:
    """Independent prctl path reserved for post-spawn emergency cleanup."""
    if not sys.platform.startswith("linux"):
        raise EvidenceError("emergency child cleanup requires Linux")
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        prctl = libc.prctl
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            f"emergency child-subreaper prctl is unavailable: {error}") from error
    ctypes.set_errno(0)
    result = prctl(ctypes.c_int(option), argument,
                   ctypes.c_ulong(0), ctypes.c_ulong(0), ctypes.c_ulong(0))
    if result != 0:
        number = ctypes.get_errno()
        raise EvidenceError(
            "emergency child-subreaper prctl failed: " +
            os.strerror(number or errno.EPERM))


def _emergency_get_child_subreaper() -> int:
    value = ctypes.c_int(-1)
    _emergency_linux_prctl(PR_GET_CHILD_SUBREAPER, ctypes.byref(value))
    if value.value not in (0, 1):
        raise EvidenceError("emergency child-subreaper state is invalid")
    return value.value


def _emergency_restore_child_subreaper(value: int) -> None:
    if value not in (0, 1):
        raise EvidenceError("emergency restore state is invalid")
    _emergency_linux_prctl(PR_SET_CHILD_SUBREAPER, ctypes.c_ulong(value))
    if _emergency_get_child_subreaper() != value:
        raise EvidenceError("emergency child-subreaper restore did not persist")


def _emergency_proc_process_record(
    pid: int,
) -> tuple[int, int, int, int, str] | None:
    """Read a Linux process identity without using the normal procfs helpers."""
    path = f"/proc/{pid}/stat"
    descriptor: int | None = None
    try:
        descriptor = os.open(
            path, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))
        chunks = bytearray()
        while len(chunks) <= 65536:
            block = os.read(descriptor, min(4096, 65537 - len(chunks)))
            if not block:
                break
            chunks.extend(block)
        if len(chunks) > 65536:
            raise EvidenceError(f"emergency procfs record {pid} is oversized")
        data = bytes(chunks)
    except OSError as error:
        if error.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise EvidenceError(
            f"cannot inspect emergency Linux process {pid}: {error}") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
    closing = data.rfind(b")")
    if closing <= 0 or closing + 2 >= len(data):
        raise EvidenceError(f"emergency Linux process {pid} has malformed stat data")
    fields = data[closing + 2:].split()
    if len(fields) < 20:
        raise EvidenceError(f"emergency Linux process {pid} has truncated stat data")
    try:
        state = fields[0].decode("ascii")
        ppid, pgrp, session = int(fields[1]), int(fields[2]), int(fields[3])
        starttime = int(fields[19])
    except (UnicodeDecodeError, ValueError) as error:
        raise EvidenceError(
            f"emergency Linux process {pid} has invalid stat fields") from error
    if (len(state) != 1 or min(ppid, pgrp, session, starttime) < 0):
        raise EvidenceError(f"emergency Linux process {pid} identity is invalid")
    return ppid, pgrp, session, starttime, state


def _emergency_proc_process_snapshot(
) -> dict[int, tuple[int, int, int, int, str]]:
    """Independent same-UID procfs snapshot for exception cleanup."""
    try:
        names = os.listdir("/proc")
    except OSError as error:
        raise EvidenceError(f"cannot enumerate emergency Linux procfs: {error}") from error
    result: dict[int, tuple[int, int, int, int, str]] = {}
    for name in names:
        if not name.isascii() or not name.isdigit():
            continue
        pid = int(name)
        try:
            record = _emergency_proc_process_record(pid)
        except EvidenceError:
            try:
                owner = os.stat(f"/proc/{name}", follow_symlinks=False).st_uid
            except OSError:
                continue
            if owner == os.getuid():
                raise
            continue
        if record is not None:
            result[pid] = record
    if os.getpid() not in result:
        raise EvidenceError("emergency procfs does not expose the runner")
    return result


def _emergency_pidfd_open(pid: int) -> int | None:
    """Independent pidfd open used only after the normal path faults."""
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_open
    except (AttributeError, OSError) as error:
        raise EvidenceError(f"emergency pidfd_open is unavailable: {error}") from error
    ctypes.set_errno(0)
    descriptor = function(ctypes.c_int(pid), ctypes.c_uint(0))
    if descriptor >= 0:
        return descriptor
    number = ctypes.get_errno()
    if number == errno.ESRCH:
        return None
    raise EvidenceError(
        f"emergency pidfd_open failed for {pid}: " +
        os.strerror(number or errno.EPERM))


def _emergency_pidfd_signal(descriptor: int, signal_number: int) -> None:
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_send_signal
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            f"emergency pidfd_send_signal is unavailable: {error}") from error
    ctypes.set_errno(0)
    result = function(ctypes.c_int(descriptor), ctypes.c_int(signal_number),
                      None, ctypes.c_uint(0))
    if result != 0:
        number = ctypes.get_errno()
        if number == errno.ESRCH:
            return
        raise EvidenceError(
            "emergency pidfd_send_signal failed: " +
            os.strerror(number or errno.EPERM))


def _emergency_signal_identity(identity: tuple[int, int]) -> None:
    """SIGKILL exactly one still-matching process identity."""
    pid, starttime = identity
    record = _emergency_proc_process_record(pid)
    if record is None or record[3] != starttime:
        return
    descriptor = _emergency_pidfd_open(pid)
    if descriptor is None:
        return
    try:
        record = _emergency_proc_process_record(pid)
        if record is None or record[3] != starttime:
            return
        _emergency_pidfd_signal(descriptor, signal.SIGKILL)
    finally:
        os.close(descriptor)


class LinuxDescendantContainment:
    """Own, kill, reap, and prove teardown of one complete Linux child tree.

    A process can escape a start_new_session process group with setsid().  The
    runner therefore becomes a temporary child subreaper.  Escaped and
    double-forked descendants are adopted by this process and remain visible
    through procfs until killed and explicitly reaped.
    """

    def __init__(self) -> None:
        self.runner_pid = os.getpid()
        self.previous_subreaper: int | None = None
        self.baseline_children: set[tuple[int, int]] = set()
        self.leader: tuple[int, int] | None = None
        self.known: set[tuple[int, int]] = set()
        self.spawned_process: subprocess.Popen[bytes] | None = None
        self.active = False
        self.proven_empty = False

    @staticmethod
    def _direct_children(
        snapshot: Mapping[int, tuple[int, int, int, int, str]], parent: int
    ) -> set[tuple[int, int]]:
        return {(pid, record[3]) for pid, record in snapshot.items()
                if record[0] == parent}

    def __enter__(self) -> "LinuxDescendantContainment":
        require(not self.active, "child descendant containment is already active")
        require(sys.platform.startswith("linux"),
                "child descendant containment requires Linux")
        task_root = Path("/proc/self/task")
        require(task_root.is_dir(),
                "child descendant containment requires Linux procfs task data")
        try:
            task_count = sum(1 for name in os.listdir(task_root)
                             if name.isascii() and name.isdigit())
        except OSError as error:
            raise EvidenceError(
                f"cannot inspect runner threads for child containment: {error}") from error
        require(task_count == 1,
                "child descendant containment requires a single-threaded runner")
        _validate_linux_pidfd_support()
        self.previous_subreaper = _emergency_get_child_subreaper()
        require(_get_child_subreaper() == self.previous_subreaper,
                "normal and emergency child-subreaper reads disagree")
        try:
            _set_child_subreaper(1)
            snapshot = _proc_process_snapshot()
            self.baseline_children = self._direct_children(
                snapshot, self.runner_pid)
            require(not self.baseline_children,
                    "child descendant containment found pre-existing children")
            self.active = True
            return self
        except BaseException:
            if self.previous_subreaper is not None:
                _emergency_restore_child_subreaper(self.previous_subreaper)
            self.previous_subreaper = None
            raise

    def observe_spawn(self, process: subprocess.Popen[bytes]) -> None:
        """Record the Popen handle before any injectable normal attachment."""
        require(self.active and self.spawned_process is process and process.pid > 0,
                "invalid emergency child observation")
        record = _emergency_proc_process_record(process.pid)
        if record is not None and record[0] == self.runner_pid:
            identity = (process.pid, record[3])
            self.known.add(identity)

    def attach(self, pid: int) -> None:
        require(self.active and self.leader is None and type(pid) is int and pid > 0,
                "invalid child descendant containment attachment")
        record = _proc_process_record(pid)
        identity = None if record is None else (pid, record[3])
        require(identity is not None and record is not None and
                record[0] == self.runner_pid and
                identity not in self.baseline_children,
                "spawned process is not an owned direct child")
        self.leader = identity
        self.known.add(identity)

    def _discover(
        self, snapshot: Mapping[int, tuple[int, int, int, int, str]]
    ) -> set[tuple[int, int]]:
        targets = {identity for identity in self.known
                   if _process_identity(identity[0], snapshot) == identity}
        targets.update(
            identity for identity in self._direct_children(
                snapshot, self.runner_pid)
            if identity not in self.baseline_children)
        changed = True
        while changed:
            changed = False
            parent_pids = {pid for pid, _starttime in targets}
            for pid, record in snapshot.items():
                identity = (pid, record[3])
                if record[0] in parent_pids and identity not in targets:
                    targets.add(identity)
                    changed = True
        self.known.update(targets)
        return targets

    @staticmethod
    def _signal_identity(identity: tuple[int, int]) -> None:
        pid, starttime = identity
        record = _proc_process_record(pid)
        if record is None or record[3] != starttime:
            return
        descriptor = _linux_pidfd_open(pid)
        if descriptor is None:
            return
        try:
            record = _proc_process_record(pid)
            if record is None or record[3] != starttime:
                return
            _linux_pidfd_signal(descriptor, signal.SIGKILL)
        finally:
            os.close(descriptor)

    def terminate_and_reap(self, process: subprocess.Popen[bytes]) -> None:
        require(self.active and self.leader is not None and
                self.leader[0] == process.pid,
                "child descendant containment is not attached to this process")
        deadline = time.monotonic() + CHILD_REAP_TIMEOUT_SECONDS
        empty_scans = 0
        while True:
            snapshot = _proc_process_snapshot()
            targets = self._discover(snapshot)
            for identity in sorted(targets, reverse=True):
                self._signal_identity(identity)

            if process.poll() is None:
                try:
                    process.wait(timeout=max(
                        0.001, min(0.05, deadline - time.monotonic())))
                except subprocess.TimeoutExpired:
                    pass

            # Subreaper adoption makes every orphan a waitable direct child.
            for identity in sorted(self.known):
                if identity == self.leader:
                    continue
                record = _proc_process_record(identity[0])
                if record is None or record[3] != identity[1]:
                    continue
                try:
                    os.waitpid(identity[0], os.WNOHANG)
                except (ChildProcessError, ProcessLookupError):
                    pass
                except OSError as error:
                    raise EvidenceError(
                        f"cannot reap contained child {identity[0]}: {error}") from error

            snapshot = _proc_process_snapshot()
            live_targets = self._discover(snapshot)
            if process.poll() is not None and not live_targets:
                empty_scans += 1
                if empty_scans >= 2:
                    self.proven_empty = True
                    return
            else:
                empty_scans = 0
            remaining = deadline - time.monotonic()
            require(remaining > 0,
                    "contained child descendants remained after SIGKILL")
            time.sleep(min(0.01, remaining))

    def _emergency_discover(
        self, snapshot: Mapping[int, tuple[int, int, int, int, str]]
    ) -> set[tuple[int, int]]:
        targets = {identity for identity in self.known
                   if _process_identity(identity[0], snapshot) == identity}
        targets.update(
            identity for identity in self._direct_children(
                snapshot, self.runner_pid)
            if identity not in self.baseline_children)
        changed = True
        while changed:
            changed = False
            parents = {pid for pid, _starttime in targets}
            for pid, record in snapshot.items():
                identity = (pid, record[3])
                if record[0] in parents and identity not in targets:
                    targets.add(identity)
                    changed = True
        self.known.update(targets)
        return targets

    def emergency_terminate_and_reap(self) -> None:
        """Exception-independent, bounded cleanup for every spawned child.

        This deliberately does not call the normal procfs or pidfd helpers, so
        a fault injected into attachment or primary teardown cannot disable the
        final cleanup path.
        """
        process = self.spawned_process
        if process is None:
            return
        deadline = time.monotonic() + CHILD_REAP_TIMEOUT_SECONDS
        empty_scans = 0
        while True:
            snapshot = _emergency_proc_process_snapshot()
            record = snapshot.get(process.pid)
            if (record is not None and record[0] == self.runner_pid and
                    process.poll() is None):
                self.known.add((process.pid, record[3]))
            targets = self._emergency_discover(snapshot)
            for identity in sorted(targets, reverse=True):
                _emergency_signal_identity(identity)

            if process.poll() is None:
                try:
                    process.wait(timeout=max(
                        0.001, min(0.05, deadline - time.monotonic())))
                except subprocess.TimeoutExpired:
                    pass

            # No pre-existing direct child is permitted at entry, so every
            # adopted child identity discovered above is safe to reap here.
            for identity in sorted(self.known):
                if identity[0] == process.pid:
                    continue
                record = _emergency_proc_process_record(identity[0])
                if record is None or record[3] != identity[1]:
                    continue
                try:
                    os.waitpid(identity[0], os.WNOHANG)
                except (ChildProcessError, ProcessLookupError):
                    pass
                except OSError as error:
                    raise EvidenceError(
                        f"emergency child reap failed for {identity[0]}: {error}") \
                        from error

            snapshot = _emergency_proc_process_snapshot()
            live = self._emergency_discover(snapshot)
            if process.poll() is not None and not live:
                empty_scans += 1
                if empty_scans >= 2:
                    self.proven_empty = True
                    return
            else:
                empty_scans = 0
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise EvidenceError(
                    "emergency child descendants remained after bounded SIGKILL")
            time.sleep(min(0.01, remaining))

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del tb
        if not self.active:
            return
        cleanup_error: BaseException | None = None
        restore_error: BaseException | None = None
        previous = self.previous_subreaper
        try:
            if (self.spawned_process is not None and
                    (exc_type is not None or not self.proven_empty)):
                try:
                    self.emergency_terminate_and_reap()
                except BaseException as error:
                    cleanup_error = error
            elif self.spawned_process is None:
                # Popen was never reached.  A new direct child would still be
                # an invariant violation, but it cannot be silently ignored.
                try:
                    snapshot = _emergency_proc_process_snapshot()
                    current = self._direct_children(snapshot, self.runner_pid)
                    if current != self.baseline_children:
                        cleanup_error = EvidenceError(
                            "unattached child appeared during containment")
                except BaseException as error:
                    cleanup_error = error
        finally:
            self.active = False
            self.previous_subreaper = None
            if previous is None:
                restore_error = EvidenceError(
                    "previous child-subreaper state was lost")
            else:
                try:
                    _emergency_restore_child_subreaper(previous)
                except BaseException as error:
                    restore_error = error
        if cleanup_error is not None or restore_error is not None:
            parts = []
            if cleanup_error is not None:
                parts.append(
                    "emergency cleanup failed: " +
                    f"{type(cleanup_error).__name__}: {cleanup_error}")
            if restore_error is not None:
                parts.append(
                    "subreaper restore failed: " +
                    f"{type(restore_error).__name__}: {restore_error}")
            if exc is not None:
                parts.append(f"primary failure: {type(exc).__name__}: {exc}")
            raise EvidenceError("; ".join(parts)) from (
                cleanup_error or restore_error)


def terminate_process_group_bounded(
    process: subprocess.Popen[bytes], timeout: float = 5.0
) -> tuple[bool, int]:
    """Bounded legacy process-group primitive; runners use full containment."""
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    try:
        return True, process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, (
            process.returncode if isinstance(process.returncode, int)
            else -int(signal.SIGKILL)
        )


def run_process_bounded(
    arguments: Sequence[str],
    cwd: Path | None = None,
    environment: Mapping[str, str] | None = None,
    timeout: float = 30.0,
    max_stdout: int = MAX_COMMAND_STDOUT_BYTES,
    max_stderr: int = MAX_COMMAND_STDERR_BYTES,
) -> subprocess.CompletedProcess[bytes]:
    require(isinstance(arguments, Sequence) and 1 <= len(arguments) <= 512 and
            all(isinstance(item, str) and item and
                len(item.encode("utf-8")) <= 65536 for item in arguments),
            "subprocess argument vector is invalid or oversized")
    require(isinstance(timeout, (int, float)) and not isinstance(timeout, bool) and
            math.isfinite(float(timeout)) and 0 < timeout <= MAX_COMMAND_TIMEOUT_SECONDS,
            "subprocess timeout is invalid")
    require(type(max_stdout) is int and 0 <= max_stdout <= MAX_COMMAND_STDOUT_BYTES and
            type(max_stderr) is int and 0 <= max_stderr <= MAX_COMMAND_STDERR_BYTES,
            "subprocess output limits are invalid")
    process: subprocess.Popen[bytes] | None = None
    selector = selectors.DefaultSelector()
    stdout_fd = -1
    stderr_fd = -1
    outputs: dict[int, bytearray] = {}
    returncode = -int(signal.SIGKILL)
    failure: str | None = None
    try:
        with LinuxDescendantContainment() as containment:
            process = subprocess.Popen(
                list(arguments), cwd=None if cwd is None else str(cwd),
                stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, start_new_session=True,
                env=None if environment is None else dict(environment))
            # Store the handle with no intervening injectable helper call.
            containment.spawned_process = process
            containment.observe_spawn(process)
            containment.attach(process.pid)
            require(process.stdout is not None and process.stderr is not None,
                    "cannot capture subprocess output")
            stdout_fd = process.stdout.fileno()
            stderr_fd = process.stderr.fileno()
            outputs = {stdout_fd: bytearray(), stderr_fd: bytearray()}
            limits = {stdout_fd: max_stdout, stderr_fd: max_stderr}
            for stream in (process.stdout, process.stderr):
                os.set_blocking(stream.fileno(), False)
                selector.register(stream, selectors.EVENT_READ)
            deadline = time.monotonic() + float(timeout)
            try:
                while selector.get_map():
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        failure = f"command exceeded {float(timeout):.3f} seconds"
                        break
                    for key, _ in selector.select(min(remaining, 0.1)):
                        descriptor = key.fileobj.fileno()
                        try:
                            block = os.read(descriptor, 65536)
                        except BlockingIOError:
                            continue
                        if not block:
                            selector.unregister(key.fileobj)
                            continue
                        outputs[descriptor].extend(block)
                        if len(outputs[descriptor]) > limits[descriptor]:
                            failure = "command output exceeded its byte limit"
                            break
                    if failure is not None:
                        break
                if failure is None:
                    try:
                        returncode = process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        failure = "command closed its output but did not terminate"
            finally:
                # This is required after successful commands too: a command may
                # exit zero after daemonizing a descendant into a new session.
                containment.terminate_and_reap(process)
                if isinstance(process.returncode, int):
                    returncode = process.returncode
    finally:
        selector.close()
        if process is not None:
            if process.stdout is not None:
                process.stdout.close()
            if process.stderr is not None:
                process.stderr.close()
    stdout = bytes(outputs[stdout_fd])
    stderr = bytes(outputs[stderr_fd])
    if failure is not None:
        raise EvidenceError(failure)
    return subprocess.CompletedProcess(list(arguments), returncode, stdout, stderr)


def run_checked(
    arguments: Sequence[str],
    cwd: Path | None = None,
    environment: Mapping[str, str] | None = None,
    timeout: float = 30.0,
    max_stdout: int = MAX_COMMAND_STDOUT_BYTES,
    max_stderr: int = MAX_COMMAND_STDERR_BYTES,
) -> bytes:
    completed = run_process_bounded(
        arguments, cwd, environment, timeout, max_stdout, max_stderr)
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace").strip()
        raise EvidenceError(
            f"command failed ({completed.returncode}): {' '.join(arguments)}: {detail}")
    return completed.stdout


def git_commit_object_identity(raw: bytes, expected_commit: str) -> dict[str, Any]:
    require(isinstance(raw, bytes) and 0 < len(raw) <= MAX_GIT_COMMIT_BYTES,
            "Git commit object is outside the retained byte bound")
    object_id = hashlib.sha1(
        f"commit {len(raw)}\0".encode("ascii") + raw).hexdigest()
    require(object_id == expected_commit,
            "Git commit object bytes do not match the retained commit ID")
    encoded = base64.b64encode(raw).decode("ascii")
    return {
        "encoding": "base64", "size": len(raw),
        "sha256": sha256_bytes(raw), "object_id": object_id,
        "base64": encoded,
    }


def validate_git_commit_object_identity(
    value: object, expected_commit: str, expected_tree: str, label: str,
) -> dict[str, Any]:
    require(isinstance(expected_commit, str) and
            re.fullmatch(r"[0-9a-f]{40}", expected_commit) is not None and
            isinstance(expected_tree, str) and
            re.fullmatch(r"[0-9a-f]{40}", expected_tree) is not None,
            f"{label} expected Git identity is invalid")
    require(isinstance(value, dict) and set(value) == {
                "encoding", "size", "sha256", "object_id", "base64"} and
            value.get("encoding") == "base64" and
            isinstance(value.get("base64"), str),
            f"{label} Git commit-object identity shape differs")
    try:
        raw = base64.b64decode(value["base64"], validate=True)
    except (ValueError, binascii.Error) as error:
        raise EvidenceError(f"{label} Git commit object is not canonical base64") from error
    require(base64.b64encode(raw).decode("ascii") == value["base64"] and
            type(value.get("size")) is int and value["size"] == len(raw) and
            0 < len(raw) <= MAX_GIT_COMMIT_BYTES and
            value.get("sha256") == sha256_bytes(raw) and
            value.get("object_id") == expected_commit and
            git_commit_object_identity(raw, expected_commit) == value,
            f"{label} Git commit-object byte identity differs")
    require(b"\n\n" in raw,
            f"{label} Git commit object has no header/message boundary")
    header_lines = raw.split(b"\n\n", 1)[0].splitlines()
    expected_tree_line = b"tree " + expected_tree.encode("ascii")
    trees = [line for line in header_lines if line.startswith(b"tree ")]
    require(header_lines and header_lines[0] == expected_tree_line and
            trees == [expected_tree_line],
            f"{label} Git commit object names a different tree")
    return value


def git_identity(
    root: Path, expected_commit: str, detached: bool,
    include_commit_object: bool = False,
) -> dict[str, Any]:
    root = root.resolve(strict=True)
    git = Path("/usr/bin/git").resolve(strict=True)
    head = run_checked((str(git), "-C", str(root), "rev-parse", "HEAD")) \
        .decode().strip()
    require(head == expected_commit,
            f"source {root} is {head}, expected {expected_commit}")
    tree = run_checked((str(git), "-C", str(root), "rev-parse", "HEAD^{tree}")) \
        .decode().strip()
    status = run_checked((
        str(git), "-C", str(root), "status", "--porcelain=v1",
        "--untracked-files=normal")).decode()
    require(status == "", f"source {root} is not clean: {status!r}")
    symbolic = run_process_bounded(
        (str(git), "-C", str(root), "symbolic-ref", "-q", "HEAD"),
        max_stdout=65536, max_stderr=65536)
    is_detached = symbolic.returncode != 0
    if detached:
        require(is_detached, f"exact-main source {root} is not detached")
    tracked = run_checked((
        str(git), "-C", str(root), "ls-tree", "-r", "-z", "HEAD"))
    result = {
        "path": str(root),
        "head": head,
        "tree": tree,
        "detached": is_detached,
        "tracked_tree_listing_sha256": sha256_bytes(tracked),
        "tracked_status": "clean",
    }
    if include_commit_object:
        commit_object = run_checked((
            str(git), "-C", str(root), "cat-file", "commit", head),
            max_stdout=MAX_GIT_COMMIT_BYTES, max_stderr=65536)
        result["commit_object"] = git_commit_object_identity(
            commit_object, head)
        validate_git_commit_object_identity(
            result["commit_object"], head, tree, f"source {root}")
    return result


def artifact_identity(path: Path, kind: str) -> dict[str, Any]:
    path = path.resolve(strict=True)
    metadata, digest, prefix = bounded_file_snapshot(path)
    if kind == "executable":
        require(os.access(path, os.X_OK), f"benchmark is not executable: {path}")
    if kind == "archive":
        require(prefix == b"!<arch>\n", f"not an ar archive: {path}")
    if kind == "shared_library":
        require(prefix.startswith(b"\x7fELF"),
                f"not an ELF shared library: {path}")
    return {
        "path": str(path),
        "kind": kind,
        "size": metadata.st_size,
        "mode": metadata.st_mode & 0o7777,
        "mtime_ns": metadata.st_mtime_ns,
        "sha256": digest,
    }


def parse_cmake_cache(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        name_and_type, value = line.split("=", 1)
        if ":" not in name_and_type:
            continue
        name, _ = name_and_type.split(":", 1)
        values[name] = value
    return values


def compile_command_tokens(entry: Mapping[str, Any]) -> list[str]:
    has_arguments = "arguments" in entry
    has_command = "command" in entry
    require(has_arguments != has_command,
            "compile command must contain exactly one argv representation")
    arguments = entry.get("arguments")
    if has_arguments:
        require(isinstance(arguments, list) and arguments and
                all(isinstance(item, str) and item for item in arguments),
                "compile command arguments are invalid")
        return list(arguments)
    command = entry.get("command")
    require(isinstance(command, str) and command,
            "compile command text is invalid")
    try:
        tokens = shlex.split(command, posix=True)
    except ValueError as error:
        raise EvidenceError(f"cannot parse compile command: {error}") from error
    require(tokens and all(tokens), "compile command token stream is empty")
    return tokens


def validate_effective_flags(
    tokens: Sequence[str], what: str, policy: str,
) -> None:
    try:
        link_common.validate_effective_flags(tokens, what, policy)
    except link_common.EvidenceError as error:
        raise EvidenceError(str(error)) from error


EXTERNAL_LINK_INPUT_ROLES = link_common.EXTERNAL_LINK_INPUT_ROLES
EXTERNAL_LINK_INPUT_ORDER = link_common.EXTERNAL_LINK_INPUT_ORDER


def _validate_external_link_operand_path(
    operand: object, role: str, label: str,
) -> str:
    try:
        return link_common.validate_external_link_operand_path(
            operand, role, label)
    except link_common.EvidenceError as error:
        raise EvidenceError(str(error)) from error


def capture_external_link_inputs(
    tokens: Sequence[str], label: str,
) -> list[dict[str, Any]]:
    """Snapshot every accepted non-project file operand before execution.

    CMake's GCC OpenMP target emits an explicit libgomp shared-library operand
    and an explicit (often empty) libpthread archive.  Their lexical operands
    are part of the recipe while their resolved regular-file bytes are the
    durable identity.  No other external file input is inferred or accepted.
    """
    records: list[dict[str, Any]] = []
    for token in tokens:
        if not isinstance(token, str):
            continue
        descriptor = EXTERNAL_LINK_INPUT_ROLES.get(Path(token).name)
        if descriptor is None:
            continue
        role, kind = descriptor
        operand = _validate_external_link_operand_path(token, role, label)
        try:
            metadata, digest, snapshot, resolved, symlink_chain = \
                link_common.current_external_file_snapshot(
                    Path(operand), f"{label} {role}")
        except (OSError, RuntimeError, link_common.EvidenceError) as error:
            raise EvidenceError(
                f"{label} cannot snapshot {role}: {error}") from error
        require((kind == "archive" and snapshot.startswith(b"!<arch>\n")) or
                (kind == "shared_library" and
                 snapshot.startswith(b"\x7fELF")),
                f"{label} {role} has the wrong file format")
        records.append({
            "operand": operand,
            "role": role,
            "lexical_symlink_chain": symlink_chain,
            "artifact": {
                "path": str(resolved),
                "kind": kind,
                "size": metadata.st_size,
                "mode": stat.S_IMODE(metadata.st_mode),
                "mtime_ns": metadata.st_mtime_ns,
                "sha256": digest,
            },
        })
    require([record["role"] for record in records] ==
                list(EXTERNAL_LINK_INPUT_ORDER),
            f"{label} external OpenMP/pthread input closure differs")
    return records


def validate_external_link_inputs(
    value: object, label: str, *, normalized: bool = False,
) -> list[dict[str, Any]]:
    """Validate the retained exact identities of accepted external inputs."""
    require(isinstance(value, list) and len(value) ==
                len(EXTERNAL_LINK_INPUT_ORDER),
            f"{label} external-link identity closure differs")
    records: list[dict[str, Any]] = []
    for index, expected_role in enumerate(EXTERNAL_LINK_INPUT_ORDER):
        record = value[index]
        require(isinstance(record, dict) and set(record) == {
                    "operand", "role", "artifact",
                    "lexical_symlink_chain"} and
                record.get("role") == expected_role,
                f"{label} external-link identity {index} shape differs")
        operand = _validate_external_link_operand_path(
            record.get("operand"), expected_role, label)
        _, expected_kind = EXTERNAL_LINK_INPUT_ROLES[Path(operand).name]
        artifact = record.get("artifact")
        if normalized:
            require(isinstance(artifact, dict) and set(artifact) == {
                        "path", "kind", "size", "mode", "sha256"} and
                    artifact.get("kind") == expected_kind and
                    isinstance(artifact.get("path"), str) and
                    Path(artifact["path"]).is_absolute() and
                    type(artifact.get("size")) is int and artifact["size"] >= 0 and
                    type(artifact.get("mode")) is int and
                    0 <= artifact["mode"] <= 0o7777 and
                    isinstance(artifact.get("sha256"), str) and
                    HEX256.fullmatch(artifact["sha256"]) is not None,
                    f"{label} normalized external artifact {index} is invalid")
        else:
            artifact = validate_complete_artifact_identity(
                artifact, f"{label} external artifact {index}", expected_kind)
        records.append(record)
    try:
        link_common.validate_external_link_input_shape(records, label)
    except link_common.EvidenceError as error:
        raise EvidenceError(str(error)) from error
    return records


def validate_executable_link_semantics(
    tokens: Sequence[str], *, compiler_invocation: str, archive_name: str,
    executable_name: str, benchmark_object: str,
    external_link_inputs: object, label: str,
    normalized_external_inputs: bool = False,
) -> None:
    """Validate one executable-link command with a fail-closed grammar.

    Every token is consumed by one exact production.  This avoids the former
    suffix blacklist, under which fused -Wl controls and many compiler-driver
    inputs were invisible to the archive/object closure checks.
    """
    external = validate_external_link_inputs(
        external_link_inputs, label, normalized=normalized_external_inputs)
    try:
        link_common.validate_executable_link_semantics(
            tokens, compiler_invocation=compiler_invocation,
            archive_name=archive_name, executable_name=executable_name,
            benchmark_object=benchmark_object,
            external_link_inputs=external, label=label)
    except link_common.EvidenceError as error:
        raise EvidenceError(str(error)) from error


def _validate_historical_declared_archive_operands(
    tokens: Sequence[str], archive_name: str, label: str,
) -> None:
    """Preserve replay of pre-v5 producers without broadening complete schemas."""
    same_name_operands = [
        token for token in tokens
        if isinstance(token, str) and Path(token).name == archive_name
    ]
    require(same_name_operands == [archive_name] and
            not any(isinstance(token, str) and
                    ("@" in token or token.startswith("-l") or ",-l" in token)
                    for token in tokens),
            f"{label} historical archive/library closure differs")
    external_archives = [
        token for token in tokens
        if isinstance(token, str) and token.endswith(".a") and
        token != archive_name
    ]
    allowed_system_roots = ("/lib/", "/lib64/", "/usr/lib/", "/usr/lib64/")
    require(all(Path(token).is_absolute() and
                os.path.normpath(token) == token and
                Path(token).name == "libpthread.a" and
                token.startswith(allowed_system_roots)
                for token in external_archives),
            f"{label} historical static archive closure differs")


def compile_profile_for_implementation(implementation: str) -> str:
    require(implementation in {"baseline", "candidate"},
            "compile profile implementation is invalid")
    return (BASELINE_COMPILE_PROFILE if implementation == "baseline" else
            CANDIDATE_COMPILE_PROFILE)


def expected_compile_output(
    implementation: str, source: Path, specification: Mapping[str, Any],
) -> str:
    """Return the one CMake object operand permitted for a configured source."""
    baseline_root = Path(specification["baseline_source_root"])
    candidate_root = Path(specification["candidate_source_root"])
    if implementation == "baseline":
        adapter = candidate_root / \
            "experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        if source == adapter:
            return ("CMakeFiles/leopard_main_benchmark.dir/"
                    "legacy_main_benchmark.cpp.o")
        require(source.is_relative_to(baseline_root),
                "baseline compile source escapes exact-main root")
        return ("CMakeFiles/leopard_main_exact.dir/" +
                source.as_posix().removeprefix("/") + ".o")

    require(implementation == "candidate" and
            source.is_relative_to(candidate_root),
            "candidate compile source escapes source root")
    relative = source.relative_to(candidate_root).as_posix()
    backend_targets = {
        "Leopard2BackendSSSE3.cpp": "leopard2_backend_ssse3.dir",
        "Leopard2BackendAVX2.cpp": "leopard2_backend_avx2.dir",
        "Leopard2BackendAVX2Xor.cpp": "leopard2_backend_avx2.dir",
        "Leopard2BackendAVX512.cpp": "leopard2_backend_avx512.dir",
        "Leopard2BackendGFNI.cpp": "leopard2_backend_gfni.dir",
    }
    if relative in CANDIDATE_LIBRARY_SOURCES:
        target = backend_targets.get(relative, "leopard.dir")
    else:
        target = CANDIDATE_NON_LIBRARY_COMPILE_TARGETS.get(relative)
        require(target is not None,
                f"candidate compile source is not configured: {source}")
    return f"CMakeFiles/{target}/{relative}.o"


def expected_compile_argv(
    implementation: str, source: Path, specification: Mapping[str, Any],
    compiler_invocation: str,
) -> list[str]:
    """Construct the complete ordered argv for one configured translation unit."""
    output = expected_compile_output(implementation, source, specification)
    baseline_root = Path(specification["baseline_source_root"])
    candidate_root = Path(specification["candidate_source_root"])
    if implementation == "baseline":
        adapter = candidate_root / \
            "experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        definitions = ([] if source != adapter else [
            f'-DLEOPARD_MAIN_SOURCE_COMMIT="{MAIN_COMMIT}"',
        ])
        return [
            compiler_invocation, *definitions, f"-I{baseline_root}",
            "-g", "-O0", "-O3", "-std=gnu++11", "-march=native",
            "-Wall", "-Wextra", "-fopenmp",
            "-o", output, "-c", str(source),
        ]

    relative = source.relative_to(candidate_root).as_posix()
    require(relative in CANDIDATE_CONFIGURED_SOURCES,
            f"candidate compile source is not configured: {source}")
    isolated_flags = {
        "Leopard2BackendSSSE3.cpp": ["-mssse3", "-mno-avx"],
        "Leopard2BackendAVX2.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2BackendAVX2Xor.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2BackendAVX512.cpp": [
            "-mavx2", "-mavx512f", "-mavx512bw", "-mavx512vl",
            "-mprefer-vector-width=256", "-falign-functions=64"],
        "Leopard2BackendGFNI.cpp": [
            "-mavx2", "-mgfni", "-mno-avx512f", "-falign-functions=64"],
    }
    if relative == "bench/leopard2/locator_benchmark.cpp":
        candidate_commit = specification.get("candidate_commit")
        require(isinstance(candidate_commit, str) and
                re.fullmatch(r"[0-9a-f]{40}", candidate_commit) is not None,
                "candidate compile profile lacks its exact source commit")
        definitions = [
            "-DLEO2_LOCATOR_SOURCE_DIRTY=0",
            f'-DLEO2_LOCATOR_SOURCE_GIT_SHA="{candidate_commit}"',
        ]
        includes = [f"-I{candidate_root}"]
        propagated_openmp = ["-fopenmp"]
    elif relative == "bench/leopard2/sparse_encode_benchmark.cpp":
        candidate_commit = specification.get("candidate_commit")
        require(isinstance(candidate_commit, str) and
                re.fullmatch(r"[0-9a-f]{40}", candidate_commit) is not None,
                "candidate compile profile lacks its exact source commit")
        definitions = [
            "-DLEO2_SPARSE_ENCODE_LIBRARY_TEST_HOOKS=0",
            "-DLEO2_SPARSE_ENCODE_SOURCE_DIRTY=0",
            f'-DLEO2_SPARSE_ENCODE_SOURCE_GIT_SHA="{candidate_commit}"',
        ]
        includes = [
            f"-I{candidate_root}",
            f"-I{candidate_root / 'tests/leopard2'}",
        ]
        propagated_openmp = []
    elif relative == "tests/leopard2/direct_oracle.cpp":
        definitions = []
        includes = [f"-I{candidate_root / 'tests/leopard2'}"]
        propagated_openmp = []
    elif relative == "tests/experiments.cpp":
        definitions = []
        includes = []
        propagated_openmp = []
    elif relative == "Leopard2BackendAVX512.cpp":
        definitions = ["-DLEO2_HAVE_AVX2_BACKEND=1"]
        includes = [f"-I{candidate_root}"]
        propagated_openmp = []
    elif relative == "Leopard2BackendGFNI.cpp":
        definitions = [
            "-DLEO2_HAVE_AVX2_BACKEND=1", "-DLEO2_HAVE_GFNI_BACKEND=1"]
        includes = [f"-I{candidate_root}"]
        propagated_openmp = []
    elif relative in isolated_flags:
        definitions = []
        includes = [f"-I{candidate_root}"]
        propagated_openmp = []
    elif relative in CANDIDATE_LIBRARY_SOURCES:
        definitions = [
            "-DLEO2_DISABLE_AVX2_CODEGEN=1",
            "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_HAVE_AVX512_BACKEND=1",
            "-DLEO2_HAVE_GFNI_BACKEND=1",
            "-DLEO2_HAVE_SSSE3_BACKEND=1",
        ]
        includes = [f"-I{candidate_root}"]
        propagated_openmp = ["-fopenmp"]
    else:
        definitions = []
        includes = [f"-I{candidate_root}"]
        propagated_openmp = ["-fopenmp"]
    isa = isolated_flags.get(relative, [])
    return [
        compiler_invocation, *definitions, *includes,
        "-Wall", "-Wextra", "-fopenmp", "-O3", "-DNDEBUG", "-O3",
        "-std=gnu++11", *isa, *propagated_openmp,
        "-o", output, "-c", str(source),
    ]


def validate_compile_entry_io(
    entry: Mapping[str, Any], tokens: Sequence[str], source: Path,
    build: Path, compiler: Path, compiler_invocation: str,
) -> Path:
    """Consume the compiler, output, and sole positional source operands."""
    expected_keys = {"directory", "file", "output"} | (
        {"arguments"} if "arguments" in entry else {"command"})
    require(set(entry) == expected_keys and
            entry.get("directory") == str(build) and
            entry.get("file") == str(source) and
            isinstance(entry.get("output"), str) and entry["output"] and
            not Path(entry["output"]).is_absolute() and
            Path(entry["output"]).as_posix() == entry["output"] and
            ".." not in Path(entry["output"]).parts and
            "\\" not in entry["output"] and "@" not in entry["output"],
            f"compile command entry metadata differs for {source}")
    require(tokens and tokens[0] == compiler_invocation and
            Path(tokens[0]).resolve(strict=True) == compiler and
            all(isinstance(token, str) and token and "@" not in token
                for token in tokens),
            f"compile command for {source} has a compiler or response-file operand")
    output_positions = [
        index for index, token in enumerate(tokens) if token == "-o"]
    require(len(output_positions) == 1 and
            output_positions[0] + 1 < len(tokens) and
            tokens[output_positions[0] + 1] == entry["output"] and
            tokens.count("-c") == 1,
            f"compile command for {source} has invalid compile/output semantics")
    positional: list[str] = []
    index = 1
    while index < len(tokens):
        token = tokens[index]
        if token == "-o":
            index += 2
            continue
        if token == "-" or not token.startswith("-"):
            positional.append(token)
        index += 1
    require(positional == [str(source)],
            f"compile command for {source} does not have one exact source operand")
    # A targeted production build need not materialize objects for unrelated
    # configured targets that still appear in compile_commands.json.  Retain
    # fail-closed path containment for every entry, while leaving strict file
    # existence and byte identity to the required provenance-TU loop below.
    output = (build / entry["output"]).resolve(strict=False)
    require(output.is_relative_to(build),
            f"compile command output escapes build directory: {output}")
    return output


def validate_compile_commands(
    path: Path,
    implementation: str,
    specification: Mapping[str, Any],
    compiler: Path,
    raw_schema: str = RAW_SCHEMA,
    compiler_invocation: str | None = None,
) -> dict[str, Any]:
    try:
        entries = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise EvidenceError(f"invalid compile_commands.json: {error}") from error
    require(isinstance(entries, list) and entries,
            f"{implementation} compile_commands.json is empty")
    build = Path(specification[f"{implementation}_build_dir"]).resolve(
        strict=True)
    if compiler_invocation is None:
        compiler_invocation = str(compiler)
    require(isinstance(compiler_invocation, str) and compiler_invocation,
            f"{implementation} compiler invocation is invalid")
    by_source: dict[Path, tuple[list[str], Mapping[str, Any], Path]] = {}
    for entry in entries:
        require(isinstance(entry, dict) and isinstance(entry.get("file"), str),
                "compile command entry is malformed")
        source = Path(entry["file"]).resolve(strict=True)
        require(source not in by_source, f"duplicate compile command for {source}")
        tokens = compile_command_tokens(entry)
        output = validate_compile_entry_io(
            entry, tokens, source, build, compiler, compiler_invocation)
        by_source[source] = (tokens, entry, output)

    baseline_root = Path(specification["baseline_source_root"]).resolve(strict=True)
    candidate_root = Path(specification["candidate_source_root"]).resolve(strict=True)
    if implementation == "baseline":
        required = {
            *(baseline_root / name for name in BASELINE_LIBRARY_SOURCES),
            candidate_root / "experiments/leopard2/main_compare/legacy_main_benchmark.cpp",
        }
    else:
        required = {
            *(candidate_root / name for name in CANDIDATE_LIBRARY_SOURCES),
            candidate_root / "bench/leopard2/benchmark.cpp",
        }
    required = {path.resolve(strict=True) for path in required}
    expected_entry_count = (
        BASELINE_EXPECTED_COMPILE_COMMAND_COUNT
        if implementation == "baseline" else
        CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT)
    if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
        configured = (required if implementation == "baseline" else {
            candidate_root / name for name in CANDIDATE_CONFIGURED_SOURCES
        })
        configured = {source.resolve(strict=True) for source in configured}
        require(len(entries) == expected_entry_count and
                set(by_source) == configured,
                f"{implementation} compile-command entry closure differs")
        for source in configured:
            tokens, entry, _ = by_source[source]
            expected_tokens = expected_compile_argv(
                implementation, source, specification, compiler_invocation)
            require(tokens == expected_tokens and
                    entry["output"] == expected_compile_output(
                        implementation, source, specification),
                    f"compile command for {source} differs from the exact "
                    f"{compile_profile_for_implementation(implementation)} profile")
    missing = sorted(str(path) for path in required - set(by_source))
    require(not missing, f"{implementation} compile commands miss sources: {missing}")
    object_records: list[dict[str, Any]] = []
    required_entries: list[dict[str, Any]] = []
    for source in required:
        tokens, entry, output = by_source[source]
        if raw_schema not in COMPLETE_EVIDENCE_SCHEMAS:
            # Historical/private build schemas predate the exact ordered argv
            # record.  Preserve replay without reintroducing broad -D/-I
            # exceptions into the current production validator.
            historical_flags = [
                token for token in tokens
                if not token.startswith(("-D", "-I"))]
            validate_effective_flags(
                historical_flags, f"historical compile command for {source}",
                "compile")
            require("-fopenmp" in tokens,
                    f"historical {source} was not compiled with OpenMP")
        source_identity = artifact_identity(source, "source_file")
        object_identity = artifact_identity(output, "object_file")
        require(object_identity["mtime_ns"] >= source_identity["mtime_ns"],
                f"object {output} predates source {source}")
        object_records.append({
            "source": source_identity,
            "object": object_identity,
        })
        if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
            required_entries.append({
                "directory": entry["directory"],
                "file": entry["file"],
                "output": entry["output"],
                "arguments": list(tokens),
            })
    result = {
        "entry_count": len(entries),
        "required_sources": sorted(str(path) for path in required),
        "validated_optimization": "-O3",
        "validated_openmp": True,
        "required_source_object_pairs": sorted(
            object_records, key=lambda record: record["source"]["path"]),
        "isa_policy": (
            "whole-build -march=native" if implementation == "baseline" else
            "portable core with ISA flags only on SSSE3, AVX2, and "
            "AVX-512VL translation units"),
    }
    if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
        result.update({
            "schema": COMPILE_COMMANDS_SCHEMA,
            "implementation": implementation,
            "profile": compile_profile_for_implementation(implementation),
            "required_entries": sorted(
                required_entries, key=lambda entry: entry["file"]),
        })
    return result


def cmake_identity_for_raw_schema(raw_schema: str) -> Mapping[str, str]:
    require(isinstance(raw_schema, str),
            "main-comparison build schema is not a string")
    identity = RAW_TO_CMAKE_IDENTITY.get(raw_schema)
    require(identity is not None, "unsupported main-comparison build schema")
    return identity


def cmake_identity_for_build_schema(build_schema: str) -> Mapping[str, str]:
    require(isinstance(build_schema, str),
            "main-comparison build schema is not a string")
    identity = BUILD_SCHEMA_TO_CMAKE_IDENTITY.get(build_schema)
    require(identity is not None, "unsupported main-comparison build schema")
    return identity


def exact_text_content(text: str, label: str) -> dict[str, Any]:
    require(isinstance(text, str), f"{label} is not text")
    encoded = text.encode("utf-8")
    require(0 < len(encoded) <= MAX_LINK_RECIPE_BYTES and "\x00" not in text,
            f"{label} is outside the retained byte bound")
    return {"encoding": "utf-8", "size": len(encoded),
            "sha256": sha256_bytes(encoded), "text": text}


def validate_archive_link_recipe_content(
    content: object,
    recipe_identity: object,
    expected_archive: str,
    expected_target_directory: str,
    label: str,
    *,
    expected_objects: Sequence[str],
    expected_archiver: str,
    expected_ranlib: str,
) -> str:
    """Validate retained CMake archive-link bytes and their path semantics.

    This deliberately parses the retained bytes rather than trusting path
    labels or a normalized derivative that can be coherently re-authored.
    Backend object-library directories are distinct CMake targets and remain
    valid, while every ordinary library object must use the declared target.
    """
    require(isinstance(content, dict) and set(content) == {
                "encoding", "sha256", "size", "text"},
            f"{label} retained content is incomplete")
    require(content.get("encoding") == "utf-8",
            f"{label} retained content has the wrong encoding")
    text = content.get("text")
    require(isinstance(text, str), f"{label} retained content is not text")
    expected_content = exact_text_content(text, label)
    require(content == expected_content,
            f"{label} retained content identity is invalid")
    require(isinstance(recipe_identity, dict) and
            recipe_identity.get("size") == expected_content["size"] and
            recipe_identity.get("sha256") == expected_content["sha256"],
            f"{label} retained bytes differ from the recipe file identity")

    commands: list[list[str]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        try:
            tokens = shlex.split(line, posix=True)
        except ValueError as error:
            raise EvidenceError(f"cannot parse {label}: {error}") from error
        require(tokens, f"{label} contains an empty command")
        commands.append(tokens)
    require(len(commands) == 2,
            f"{label} must contain exactly archive and ranlib commands")

    require(isinstance(expected_archiver, str) and expected_archiver and
            isinstance(expected_ranlib, str) and expected_ranlib and
            isinstance(expected_objects, Sequence) and
            not isinstance(expected_objects, (str, bytes)) and
            expected_objects and all(isinstance(item, str) and item
                                     for item in expected_objects),
            f"{label} expected command closure is invalid")
    archive_tokens, ranlib_tokens = commands
    require(len(archive_tokens) >= 4 and
            archive_tokens[0] == expected_archiver and
            archive_tokens[1] in {"qc", "rc", "rcs"},
            f"{label} archive tool or mode differs from its build identity")
    archive_output = archive_tokens[2]
    require(archive_output == expected_archive,
            f"{label} archive output is not the canonical relative path")
    objects = [token.replace("\\", "/") for token in archive_tokens[3:]]
    require(objects == list(expected_objects) and
            all(token.endswith(".o") and not token.startswith("/") and
                "\\" not in token and "@" not in token for token in objects),
            f"{label} object order or closure differs from compile provenance")
    ordinary_prefix = f"CMakeFiles/{expected_target_directory}/"
    ordinary = [token for token in objects if
                "/leopard2_backend_" not in f"/{token}" and
                "/leopard_backend_" not in f"/{token}"]
    require(ordinary and all(ordinary_prefix in token for token in ordinary),
            f"{label} ordinary objects use the wrong CMake target directory")
    require(all(
                ordinary_prefix in token or
                re.search(r"(?:^|/)CMakeFiles/leopard2?_backend_[^/]+\.dir/", token)
                is not None
                for token in objects),
            f"{label} contains an object from an undeclared CMake target")
    require(len(ranlib_tokens) == 2 and
            ranlib_tokens[0] == expected_ranlib and
            ranlib_tokens[1] == archive_output,
            f"{label} ranlib command does not identify the produced archive")
    return text


def build_provenance(
    implementation: str, specification: Mapping[str, Any],
    raw_schema: str = RAW_SCHEMA,
) -> dict[str, Any]:
    build = Path(specification[f"{implementation}_build_dir"]).resolve(strict=True)
    require(build.is_dir(), f"{implementation} build path is not a directory: {build}")
    cmake_identity = cmake_identity_for_build_schema(raw_schema)
    names = ({
        "executable": "leopard_main_benchmark",
        "archive": "libleopard_main_exact.a",
        "executable_link": "CMakeFiles/leopard_main_benchmark.dir/link.txt",
        "archive_link": "CMakeFiles/leopard_main_exact.dir/link.txt",
    } if implementation == "baseline" else {
        "executable": "bench_leopard2",
        "archive": cmake_identity["archive"],
        "executable_link": "CMakeFiles/bench_leopard2.dir/link.txt",
        "archive_link": "CMakeFiles/{}/link.txt".format(
            cmake_identity["target_directory"]),
    })
    expected_executable = (build / names["executable"]).resolve(strict=True)
    expected_archive = (build / names["archive"]).resolve(strict=True)
    require(expected_executable ==
            Path(specification[f"{implementation}_executable"]).resolve(strict=True),
            f"{implementation} executable is not the declared build artifact")
    require(expected_archive ==
            Path(specification[f"{implementation}_archive"]).resolve(strict=True),
            f"{implementation} archive is not the declared build artifact")

    cache_path = build / "CMakeCache.txt"
    commands_path = build / "compile_commands.json"
    executable_link_path = build / names["executable_link"]
    archive_link_path = build / names["archive_link"]
    cache = parse_cmake_cache(cache_path)
    require(cache.get("CMAKE_BUILD_TYPE") == "Release",
            f"{implementation} build is not CMake Release")
    validate_effective_flags(
        shlex.split(cache.get("CMAKE_CXX_FLAGS_RELEASE", "")),
        f"{implementation} CMake Release flags", "release")
    compiler = Path(cache.get("CMAKE_CXX_COMPILER", "")).resolve(strict=True)
    candidate_root = Path(specification["candidate_source_root"]).resolve(strict=True)
    baseline_root = Path(specification["baseline_source_root"]).resolve(strict=True)
    if implementation == "baseline":
        expected_home = candidate_root / "experiments/leopard2/main_compare"
        require(Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve() ==
                expected_home.resolve(), "baseline CMake source is not the adapter directory")
        require(Path(cache.get("LEOPARD_MAIN_SOURCE_DIR", "")).resolve() == baseline_root,
                "baseline CMake cache points at the wrong exact-main source")
        required_cache = {"LEO_MAIN_HAS_MARCH_NATIVE": "1"}
        expected_archive_name = "libleopard_main_exact.a"
    else:
        require(Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve() == candidate_root,
                "candidate CMake source differs from candidate source root")
        required_cache = {
            "ENABLE_OPENMP": "ON",
            "LEO2_BACKEND_VARIANT": "auto",
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_BUILD_TESTS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
        }
        expected_archive_name = cmake_identity["archive"]
    for name, expected in required_cache.items():
        require(cache.get(name) == expected,
                f"{implementation} CMake cache {name} is {cache.get(name)!r}, "
                f"expected {expected!r}")
    executable_link = executable_link_path.read_text(encoding="utf-8")
    archive_link_bytes = archive_link_path.read_bytes()
    require(0 < len(archive_link_bytes) <= MAX_LINK_RECIPE_BYTES,
            f"{implementation} archive recipe is outside the retained byte bound")
    try:
        archive_link = archive_link_bytes.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise EvidenceError(
            f"{implementation} archive recipe is not strict UTF-8") from error
    require(names["archive"] in archive_link,
            f"{implementation} archive recipe does not produce its declared archive")
    executable_link_tokens = shlex.split(executable_link)
    require(executable_link_tokens and
            Path(executable_link_tokens[0]).resolve(strict=True) == compiler,
            f"{implementation} link recipe uses a different compiler")
    external_link_inputs: list[dict[str, Any]] = []
    if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
        external_link_inputs = capture_external_link_inputs(
            executable_link_tokens,
            f"{implementation} benchmark link recipe")
    else:
        _validate_historical_declared_archive_operands(
            executable_link_tokens, expected_archive_name,
            f"{implementation} benchmark link recipe")
    validate_effective_flags(
        executable_link_tokens, f"{implementation} executable link recipe",
        "link")
    semantics = validate_compile_commands(
        commands_path, implementation, specification, compiler, raw_schema,
        cache["CMAKE_CXX_COMPILER"])
    records = semantics["required_source_object_pairs"]
    benchmark_suffix = (
        "/experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        if implementation == "baseline" else "/bench/leopard2/benchmark.cpp")
    benchmark_records = [
        record for record in records
        if record["source"]["path"].endswith(benchmark_suffix)
    ]
    require(len(benchmark_records) == 1,
            f"{implementation} has no unique benchmark object")
    archive_records = [record for record in records if record not in benchmark_records]

    def linked_objects(recipe: str) -> list[Path]:
        result: list[Path] = []
        for line in recipe.splitlines():
            for token in shlex.split(line):
                if token.endswith(".o"):
                    item = Path(token)
                    if not item.is_absolute():
                        item = build / item
                    result.append(item.resolve(strict=True))
        return result

    archive_recipe_objects = linked_objects(archive_link)
    executable_recipe_objects = linked_objects(executable_link)
    expected_archive_objects = [
        Path(record["object"]["path"]) for record in archive_records
    ]
    expected_executable_objects = [
        Path(record["object"]["path"]) for record in benchmark_records
    ]
    require(len(archive_recipe_objects) == len(set(archive_recipe_objects)) and
            set(archive_recipe_objects) == set(expected_archive_objects),
            f"{implementation} archive object closure differs from compile commands")
    require(executable_recipe_objects == expected_executable_objects,
            f"{implementation} executable object closure differs from compile commands")
    if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
        validate_executable_link_semantics(
            executable_link_tokens,
            compiler_invocation=cache["CMAKE_CXX_COMPILER"],
            archive_name=expected_archive_name,
            executable_name=names["executable"],
            benchmark_object=Path(expected_executable_objects[0]).relative_to(
                build).as_posix(),
            external_link_inputs=external_link_inputs,
            label=f"{implementation} benchmark link recipe")
    archive_identity = artifact_identity(expected_archive, "archive")
    executable_identity = artifact_identity(expected_executable, "executable")
    require(all(archive_identity["mtime_ns"] >= record["object"]["mtime_ns"]
                for record in archive_records),
            f"{implementation} archive predates one of its object files")
    require(executable_identity["mtime_ns"] >= archive_identity["mtime_ns"] and
            all(executable_identity["mtime_ns"] >= record["object"]["mtime_ns"]
                for record in benchmark_records),
            f"{implementation} executable predates its link inputs")
    if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
        require(all(
                    executable_identity["mtime_ns"] >=
                        record["artifact"]["mtime_ns"]
                    for record in external_link_inputs),
                f"{implementation} executable predates an external link input")
    ar = Path(cache.get("CMAKE_AR", "")).resolve(strict=True)
    ranlib = Path(cache.get("CMAKE_RANLIB", "")).resolve(strict=True)
    members = run_checked((str(ar), "t", str(expected_archive)),
                          environment=CHILD_ENVIRONMENT).decode().splitlines()
    require(members == [path.name for path in archive_recipe_objects],
            f"{implementation} archive members differ from its link recipe")
    compiler_version = run_checked(
        (str(compiler), "--version"), environment=CHILD_ENVIRONMENT)
    archive_link_identity = artifact_identity(
        archive_link_path, "build_metadata")
    result = {
        "build_dir": str(build),
        "cmake_cache": artifact_identity(cache_path, "build_metadata"),
        "compile_commands": artifact_identity(commands_path, "build_metadata"),
        "executable_link_recipe": artifact_identity(
            executable_link_path, "build_metadata"),
        "archive_link_recipe": archive_link_identity,
        "compiler": artifact_identity(compiler, "compiler"),
        "compiler_version_stdout": {
            "sha256": sha256_bytes(compiler_version),
            "text": compiler_version.decode("utf-8", errors="strict"),
        },
        "archiver": artifact_identity(ar, "archiver"),
        "validated_archive_members": members,
        "validated_executable": executable_identity,
        "validated_archive": archive_identity,
        "validated_cache": {
            "CMAKE_BUILD_TYPE": cache["CMAKE_BUILD_TYPE"],
            "CMAKE_CXX_COMPILER": cache.get("CMAKE_CXX_COMPILER"),
            "CMAKE_CXX_FLAGS_RELEASE": cache["CMAKE_CXX_FLAGS_RELEASE"],
            **required_cache,
        },
        "validated_compile_commands": semantics,
    }
    if raw_schema in HARDENED_BUILD_SCHEMAS:
        content = exact_text_content(
            archive_link, f"{implementation} archive link recipe")
        require(content["size"] == archive_link_identity["size"] and
                content["sha256"] == archive_link_identity["sha256"],
                f"{implementation} archive recipe changed between reads")
        result["archive_link_recipe_content"] = content
        result["ranlib"] = artifact_identity(ranlib, "ranlib")
        result["archive_link_tool_invocations"] = {
            "archiver": {"invocation": cache["CMAKE_AR"],
                         "resolved_path": str(ar)},
            "ranlib": {"invocation": cache["CMAKE_RANLIB"],
                       "resolved_path": str(ranlib)},
        }
    if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
        executable_content = exact_text_content(
            executable_link, f"{implementation} executable link recipe")
        require(executable_content["size"] ==
                    result["executable_link_recipe"]["size"] and
                executable_content["sha256"] ==
                    result["executable_link_recipe"]["sha256"],
                f"{implementation} executable recipe changed between reads")
        result["executable_link_recipe_content"] = executable_content
        result["compiler_invocation"] = {
            "invocation": cache["CMAKE_CXX_COMPILER"],
            "resolved_path": str(compiler),
        }
        result["validated_external_link_inputs"] = external_link_inputs
    return result


def _parse_ldd_output_lines(
    output: str, label: str, *, canonical: bool,
) -> list[dict[str, Any]]:
    require(isinstance(output, str) and output and
            len(output.encode("utf-8")) <= MAX_LINK_RECIPE_BYTES,
            f"{label} ldd output is outside the retained byte bound")
    address = (re.escape(CANONICAL_LDD_ADDRESS) if canonical else
               r"0x[0-9A-Fa-f]+")
    shared = re.compile(
        rf"(?P<soname>[^\s=()]+)\s+=>\s+(?P<path>/\S+)\s+"
        rf"\((?P<address>{address})\)")
    direct = re.compile(
        rf"(?P<path>/\S+)\s+\((?P<address>{address})\)")
    virtual = re.compile(
        rf"(?P<soname>linux-vdso\.so\.1)\s+"
        rf"\((?P<address>{address})\)")
    entries: list[dict[str, Any]] = []
    for raw_line in output.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        match = shared.fullmatch(line)
        if match is not None:
            target = match.group("path")
            require(os.path.normpath(target) == target,
                    f"non-canonical runtime dependency: {line}")
            entries.append({
                "soname": match.group("soname"),
                "loader_path": target,
                "file_kind": "shared_library",
            })
            continue
        if "=>" in line and "not found" in line:
            raise EvidenceError(f"runtime dependency is missing: {line}")
        match = direct.fullmatch(line)
        if match is not None:
            target = match.group("path")
            require(os.path.normpath(target) == target,
                    f"non-canonical dynamic-loader path: {line}")
            entries.append({
                "soname": Path(target).name,
                "loader_path": target,
                "file_kind": "dynamic_loader",
            })
            continue
        match = virtual.fullmatch(line)
        if match is not None:
            entries.append({"soname": match.group("soname"), "virtual": True})
            continue
        raise EvidenceError(f"unrecognized ldd output: {line}")
    require(entries, f"{label} contains no runtime closure")
    require(len({entry["soname"] for entry in entries}) == len(entries),
            f"duplicate runtime dependency in {label}")
    require(sum(entry.get("file_kind") == "dynamic_loader" for entry in entries) == 1,
            f"{label} does not contain exactly one dynamic loader")
    require(any(entry.get("file_kind") == "shared_library" for entry in entries),
            f"{label} contains no shared-library dependency")
    require(sum(entry.get("virtual") is True for entry in entries) <= 1,
            f"{label} contains duplicate virtual loader records")
    return entries


def _render_canonical_ldd_output(entries: Sequence[Mapping[str, Any]]) -> str:
    lines: list[str] = []
    for entry in entries:
        if entry.get("virtual") is True:
            lines.append(f"{entry['soname']} ({CANONICAL_LDD_ADDRESS})")
        elif entry.get("file_kind") == "dynamic_loader":
            lines.append(f"{entry['loader_path']} ({CANONICAL_LDD_ADDRESS})")
        else:
            lines.append(
                f"{entry['soname']} => {entry['loader_path']} "
                f"({CANONICAL_LDD_ADDRESS})")
    return "\n".join(lines) + "\n"


def parse_ldd_output_text(output: str, label: str) -> list[dict[str, Any]]:
    """Strictly parse one real ldd transcript, including every load address."""
    entries = _parse_ldd_output_lines(output, label, canonical=False)
    return sorted(entries, key=lambda item: item["soname"])


def canonical_ldd_output(output: str, label: str) -> dict[str, Any]:
    """Retain all loader semantics while removing only ASLR-variant addresses."""
    entries = _parse_ldd_output_lines(output, label, canonical=False)
    text = _render_canonical_ldd_output(entries)
    identity = exact_text_content(text, f"canonical {label}")
    return {
        "schema": CANONICAL_LDD_SCHEMA,
        "normalization": CANONICAL_LDD_NORMALIZATION,
        **identity,
    }


def parse_canonical_ldd_output_text(
    output: str, label: str,
) -> list[dict[str, Any]]:
    entries = _parse_ldd_output_lines(output, label, canonical=True)
    require(output == _render_canonical_ldd_output(entries),
            f"{label} is not the exact canonical ldd rendering")
    return sorted(entries, key=lambda item: item["soname"])


def runtime_closure(
    ldd: Path, executable: Path, retained_output_schema: str | None = None,
) -> dict[str, Any]:
    ldd = ldd.resolve(strict=True)
    executable = executable.resolve(strict=True)
    output = run_checked(
        (str(ldd), str(executable)), environment=CHILD_ENVIRONMENT,
        max_stdout=MAX_LINK_RECIPE_BYTES,
    ).decode("utf-8", errors="strict")
    parsed = parse_ldd_output_text(output, f"ldd output for {executable}")
    entries = []
    for entry in parsed:
        if entry.get("virtual") is True:
            entries.append(dict(entry))
            continue
        retained = dict(entry)
        kind = retained.pop("file_kind")
        loader_path = retained["loader_path"]
        require(os.path.normpath(loader_path) == loader_path,
                f"ldd returned a non-canonical loader path: {loader_path}")
        file_identity = artifact_identity(Path(loader_path), kind)
        # The loader may report a symlinked namespace such as /lib while
        # Path.resolve() names /usr/lib.  Hash through the reported path, but
        # retain that exact path so offline validation can bind this file record
        # to its originating ldd line rather than accepting a swapped record.
        file_identity["path"] = loader_path
        retained["file"] = file_identity
        entries.append(retained)
    result = {
        "executable": str(executable),
        "dependencies": entries,
    }
    if retained_output_schema == RAW_SCHEMA_V5:
        result["raw_ldd_output"] = exact_text_content(
            output, f"ldd output for {executable}")
    elif retained_output_schema == RAW_SCHEMA:
        result["canonical_ldd_output"] = canonical_ldd_output(
            output, f"ldd output for {executable}")
    else:
        require(retained_output_schema is None,
                "unsupported retained ldd-output schema")
    return result


def input_snapshot(
    specification: Mapping[str, Any], raw_schema: str = RAW_SCHEMA,
) -> dict[str, Any]:
    ldd = Path(specification["ldd"])
    baseline_build = build_provenance("baseline", specification, raw_schema)
    candidate_build = build_provenance("candidate", specification, raw_schema)
    require(baseline_build["compiler"] == candidate_build["compiler"] and
            baseline_build["compiler_version_stdout"] ==
            candidate_build["compiler_version_stdout"],
            "baseline and candidate use different compiler binaries or versions")
    return {
        "runner": artifact_identity(Path(specification["runner"]), "file"),
        "taskset": artifact_identity(Path(specification["taskset"]), "executable"),
        "ldd": artifact_identity(ldd, "executable"),
        "baseline_executable": artifact_identity(
            Path(specification["baseline_executable"]), "executable"),
        "candidate_executable": artifact_identity(
            Path(specification["candidate_executable"]), "executable"),
        "baseline_archive": artifact_identity(
            Path(specification["baseline_archive"]), "archive"),
        "candidate_archive": artifact_identity(
            Path(specification["candidate_archive"]), "archive"),
        "baseline_build": baseline_build,
        "candidate_build": candidate_build,
        "baseline_runtime_closure": runtime_closure(
            ldd, Path(specification["baseline_executable"]),
            retained_output_schema=(raw_schema if raw_schema in
                                    COMPLETE_EVIDENCE_SCHEMAS else None)),
        "candidate_runtime_closure": runtime_closure(
            ldd, Path(specification["candidate_executable"]),
            retained_output_schema=(raw_schema if raw_schema in
                                    COMPLETE_EVIDENCE_SCHEMAS else None)),
        "baseline_source": git_identity(
            Path(specification["baseline_source_root"]), MAIN_COMMIT, True,
            include_commit_object=raw_schema in COMPLETE_EVIDENCE_SCHEMAS),
        "candidate_source": git_identity(
            Path(specification["candidate_source_root"]),
            str(specification["candidate_commit"]), False,
            include_commit_object=raw_schema in COMPLETE_EVIDENCE_SCHEMAS),
    }


def input_snapshots_equal(
    current: object, retained: object, raw_schema: str,
) -> bool:
    """Compare inputs while honoring the declared historical transcript schema."""
    if raw_schema != RAW_SCHEMA_V5:
        return current == retained
    require(isinstance(current, dict) and isinstance(retained, dict),
            "schema-v5 input snapshots are not objects")

    def normalized(snapshot: dict[str, Any]) -> dict[str, Any]:
        result = copy.deepcopy(snapshot)
        for role in ("baseline", "candidate"):
            closure = result.get(f"{role}_runtime_closure")
            require(isinstance(closure, dict),
                    f"schema-v5 {role} runtime closure is missing")
            raw = validate_complete_text_identity(
                closure.pop("raw_ldd_output", None),
                f"schema-v5 {role} raw ldd output")
            closure["canonical_ldd_output"] = canonical_ldd_output(
                raw["text"], f"schema-v5 {role} raw ldd output")
        return result

    return normalized(current) == normalized(retained)


def validate_complete_artifact_identity(
    value: object, label: str, expected_kind: str | None = None,
) -> dict[str, Any]:
    """Validate the complete, portable file identity retained since schema v5."""
    require(isinstance(value, dict) and set(value) == {
                "path", "kind", "size", "mode", "mtime_ns", "sha256"},
            f"{label} file identity shape differs")
    path = value["path"]
    kind = value["kind"]
    require(isinstance(path, str) and path and Path(path).is_absolute() and
            isinstance(kind, str) and kind and
            (expected_kind is None or kind == expected_kind) and
            type(value["size"]) is int and value["size"] >= 0 and
            type(value["mode"]) is int and 0 <= value["mode"] <= 0o7777 and
            type(value["mtime_ns"]) is int and value["mtime_ns"] >= 0 and
            isinstance(value["sha256"], str) and
            HEX256.fullmatch(value["sha256"]) is not None,
            f"{label} file identity is invalid")
    return value


def validate_complete_text_identity(value: object, label: str) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "encoding", "size", "sha256", "text"},
            f"{label} retained text shape differs")
    require(value.get("encoding") == "utf-8" and
            isinstance(value.get("text"), str),
            f"{label} retained text is invalid")
    require(value == exact_text_content(value["text"], label),
            f"{label} retained text identity differs")
    return value


def validate_canonical_ldd_output(
    value: object, label: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    require(isinstance(value, dict) and set(value) == {
                "schema", "normalization", "encoding", "size", "sha256", "text"},
            f"{label} canonical ldd-output shape differs")
    require(value.get("schema") == CANONICAL_LDD_SCHEMA and
            value.get("normalization") == CANONICAL_LDD_NORMALIZATION,
            f"{label} canonical ldd-output contract differs")
    retained = validate_complete_text_identity(
        {key: value[key] for key in ("encoding", "size", "sha256", "text")},
        f"{label} canonical ldd output")
    parsed = parse_canonical_ldd_output_text(
        retained["text"], f"{label} retained canonical ldd output")
    return value, parsed


def validate_complete_git_identity(
    value: object, label: str, expected_path: str, expected_head: str,
    *, require_detached: bool,
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "path", "head", "tree", "detached",
                "tracked_tree_listing_sha256", "tracked_status",
                "commit_object"},
            f"{label} source identity shape differs")
    require(Path(expected_path).is_absolute() and
            value.get("path") == expected_path and
            value.get("head") == expected_head and
            isinstance(value.get("tree"), str) and
            re.fullmatch(r"[0-9a-f]{40}", value["tree"]) is not None and
            type(value.get("detached")) is bool and
            (not require_detached or value["detached"] is True) and
            isinstance(value.get("tracked_tree_listing_sha256"), str) and
            HEX256.fullmatch(value["tracked_tree_listing_sha256"]) is not None and
            value.get("tracked_status") == "clean",
            f"{label} source identity is invalid")
    validate_git_commit_object_identity(
        value.get("commit_object"), expected_head, value["tree"], label)
    return value


def validate_complete_runtime_closure(
    value: object, label: str, expected_executable: str, raw_schema: str,
) -> dict[str, Any]:
    output_key = (
        "raw_ldd_output" if raw_schema == RAW_SCHEMA_V5 else
        "canonical_ldd_output" if raw_schema == RAW_SCHEMA else None)
    require(output_key is not None,
            f"{label} runtime closure uses an unsupported evidence schema")
    require(isinstance(value, dict) and set(value) == {
                "executable", "dependencies", output_key} and
            value.get("executable") == expected_executable and
            isinstance(value.get("dependencies"), list) and
            value["dependencies"],
            f"{label} runtime closure shape differs")
    if raw_schema == RAW_SCHEMA_V5:
        raw_output = validate_complete_text_identity(
            value.get(output_key), f"{label} raw ldd output")
        parsed = parse_ldd_output_text(
            raw_output["text"], f"{label} retained raw ldd output")
    else:
        _, parsed = validate_canonical_ldd_output(
            value.get(output_key), label)
    sonames: list[str] = []
    loader_paths: list[str] = []
    file_paths: list[str] = []
    for index, dependency in enumerate(value["dependencies"]):
        require(isinstance(dependency, dict),
                f"{label} runtime dependency {index} is not an object")
        soname = dependency.get("soname")
        require(isinstance(soname, str) and soname,
                f"{label} runtime dependency {index} has no soname")
        sonames.append(soname)
        if set(dependency) == {"soname", "virtual"}:
            require(soname == "linux-vdso.so.1" and
                    dependency.get("virtual") is True,
                    f"{label} virtual runtime dependency is invalid")
            continue
        require(set(dependency) == {"soname", "loader_path", "file"},
                f"{label} runtime dependency {soname} variant differs")
        loader_path = dependency.get("loader_path")
        require(isinstance(loader_path, str) and
                Path(loader_path).is_absolute(),
                f"{label} runtime dependency {soname} loader path is invalid")
        file_record = validate_complete_artifact_identity(
            dependency.get("file"), f"{label} runtime dependency {soname}")
        require(file_record["kind"] in {"shared_library", "dynamic_loader"} and
                file_record["path"] == loader_path,
                f"{label} runtime dependency {soname} has the wrong kind or loader path")
        loader_paths.append(loader_path)
        file_paths.append(file_record["path"])
    require(sonames == sorted(sonames) and len(sonames) == len(set(sonames)) and
            len(loader_paths) == len(set(loader_paths)) and
            len(file_paths) == len(set(file_paths)),
            f"{label} runtime dependency closure is not sorted and unique")
    normalized_dependencies = []
    for dependency in value["dependencies"]:
        if dependency.get("virtual") is True:
            normalized_dependencies.append(dict(dependency))
        else:
            normalized_dependencies.append({
                "soname": dependency["soname"],
                "loader_path": dependency["loader_path"],
                "file_kind": dependency["file"]["kind"],
            })
    require(normalized_dependencies == parsed,
            f"{label} dependency summary differs from retained ldd output")
    return value


def validate_complete_executable_recipe(
    content: object, identity: object, compiler_invocation: str,
    archive_name: str, executable_name: str, benchmark_object: str, label: str,
    external_link_inputs: object,
) -> None:
    retained = validate_complete_text_identity(content, label)
    recipe = validate_complete_artifact_identity(
        identity, f"{label} file", "build_metadata")
    require(recipe["size"] == retained["size"] and
            recipe["sha256"] == retained["sha256"],
            f"{label} bytes differ from its file identity")
    try:
        tokens = shlex.split(retained["text"], posix=True)
    except ValueError as error:
        raise EvidenceError(f"cannot parse {label}: {error}") from error
    validate_executable_link_semantics(
        tokens, compiler_invocation=compiler_invocation,
        archive_name=archive_name, executable_name=executable_name,
        benchmark_object=benchmark_object,
        external_link_inputs=external_link_inputs, label=label)


def validate_complete_build_identity(
    build: object, implementation: str, specification: Mapping[str, Any],
) -> dict[str, Any]:
    require(implementation in {"baseline", "candidate"} and
            isinstance(build, dict),
            f"{implementation} build identity is not an object")
    require(set(build) == {
                "build_dir", "cmake_cache", "compile_commands",
                "executable_link_recipe", "archive_link_recipe", "compiler",
                "compiler_version_stdout", "archiver", "ranlib",
                "validated_archive_members", "validated_executable",
                "validated_archive", "validated_cache",
                "validated_compile_commands", "archive_link_recipe_content",
                "executable_link_recipe_content",
                "archive_link_tool_invocations", "compiler_invocation",
                "validated_external_link_inputs"},
            f"{implementation} build identity shape differs")
    build_dir = build.get("build_dir")
    require(isinstance(build_dir, str) and Path(build_dir).is_absolute() and
            build_dir == specification[f"{implementation}_build_dir"],
            f"{implementation} build directory identity differs")
    cmake = CANONICAL_CMAKE_IDENTITY
    if implementation == "baseline":
        archive_name = "libleopard_main_exact.a"
        executable_name = "leopard_main_benchmark"
        archive_target_directory = "leopard_main_exact.dir"
        executable_recipe_relative = \
            "CMakeFiles/leopard_main_benchmark.dir/link.txt"
        required_cache = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_CXX_COMPILER": None,
            "CMAKE_CXX_FLAGS_RELEASE": None,
            "LEO_MAIN_HAS_MARCH_NATIVE": "1",
        }
        benchmark_suffix = \
            "/experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        isa_policy = "whole-build -march=native"
        compile_profile = BASELINE_COMPILE_PROFILE
        library_sources = BASELINE_LIBRARY_SOURCES
        expected_entry_count = BASELINE_EXPECTED_COMPILE_COMMAND_COUNT
    else:
        archive_name = cmake["archive"]
        executable_name = "bench_leopard2"
        archive_target_directory = cmake["target_directory"]
        executable_recipe_relative = "CMakeFiles/bench_leopard2.dir/link.txt"
        required_cache = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_CXX_COMPILER": None,
            "CMAKE_CXX_FLAGS_RELEASE": None,
            "ENABLE_OPENMP": "ON",
            "LEO2_BACKEND_VARIANT": "auto",
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_BUILD_TESTS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
        }
        benchmark_suffix = "/bench/leopard2/benchmark.cpp"
        isa_policy = (
            "portable core with ISA flags only on SSSE3, AVX2, and "
            "AVX-512VL translation units")
        compile_profile = CANDIDATE_COMPILE_PROFILE
        library_sources = CANDIDATE_LIBRARY_SOURCES
        expected_entry_count = CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT

    expected_metadata = {
        "cmake_cache": str(Path(build_dir) / "CMakeCache.txt"),
        "compile_commands": str(Path(build_dir) / "compile_commands.json"),
        "executable_link_recipe": str(
            Path(build_dir) / executable_recipe_relative),
        "archive_link_recipe": str(
            Path(build_dir) / "CMakeFiles" / archive_target_directory / "link.txt"),
    }
    for name, expected_path in expected_metadata.items():
        record = validate_complete_artifact_identity(
            build.get(name), f"{implementation} {name}", "build_metadata")
        require(record["path"] == expected_path,
                f"{implementation} {name} path differs")
    compiler = validate_complete_artifact_identity(
        build.get("compiler"), f"{implementation} compiler", "compiler")
    archiver = validate_complete_artifact_identity(
        build.get("archiver"), f"{implementation} archiver", "archiver")
    ranlib = validate_complete_artifact_identity(
        build.get("ranlib"), f"{implementation} ranlib", "ranlib")
    executable = validate_complete_artifact_identity(
        build.get("validated_executable"),
        f"{implementation} validated executable", "executable")
    archive = validate_complete_artifact_identity(
        build.get("validated_archive"),
        f"{implementation} validated archive", "archive")
    require(executable["path"] == str(Path(build_dir) / executable_name) and
            archive["path"] == str(Path(build_dir) / archive_name),
            f"{implementation} validated output paths differ")
    external_link_inputs = validate_external_link_inputs(
        build.get("validated_external_link_inputs"),
        f"{implementation} executable link recipe")
    require(all(executable["mtime_ns"] >= record["artifact"]["mtime_ns"]
                for record in external_link_inputs),
            f"{implementation} executable predates an external link input")

    version = build.get("compiler_version_stdout")
    require(isinstance(version, dict) and set(version) == {"sha256", "text"} and
            isinstance(version.get("text"), str) and version["text"] and
            isinstance(version.get("sha256"), str) and
            version["sha256"] == sha256_bytes(version["text"].encode("utf-8")),
            f"{implementation} compiler version identity differs")
    cache = build.get("validated_cache")
    require(isinstance(cache, dict) and set(cache) == set(required_cache),
            f"{implementation} validated CMake cache shape differs")
    require(cache.get("CMAKE_BUILD_TYPE") == "Release" and
            isinstance(cache.get("CMAKE_CXX_COMPILER"), str) and
            cache["CMAKE_CXX_COMPILER"] and
            isinstance(cache.get("CMAKE_CXX_FLAGS_RELEASE"), str) and
            all(cache.get(name) == expected for name, expected in
                required_cache.items() if expected is not None),
            f"{implementation} validated CMake cache semantics differ")
    try:
        cache_flags = shlex.split(cache["CMAKE_CXX_FLAGS_RELEASE"], posix=True)
    except ValueError as error:
        raise EvidenceError(
            f"cannot parse {implementation} CMake Release flags: {error}") from error
    validate_effective_flags(
        cache_flags, f"{implementation} retained CMake Release flags",
        "release")
    compiler_invocation = build.get("compiler_invocation")
    require(isinstance(compiler_invocation, dict) and
            set(compiler_invocation) == {"invocation", "resolved_path"} and
            compiler_invocation.get("invocation") ==
                cache["CMAKE_CXX_COMPILER"] and
            compiler_invocation.get("resolved_path") == compiler["path"],
            f"{implementation} compiler invocation identity differs")

    compile_record = build.get("validated_compile_commands")
    require(isinstance(compile_record, dict) and set(compile_record) == {
                "entry_count", "required_sources", "validated_optimization",
                "validated_openmp", "required_source_object_pairs", "isa_policy",
                "schema", "implementation", "profile", "required_entries"},
            f"{implementation} compile-command identity shape differs")
    pairs = compile_record.get("required_source_object_pairs")
    sources = compile_record.get("required_sources")
    retained_entries = compile_record.get("required_entries")
    require(type(compile_record.get("entry_count")) is int and
            compile_record["entry_count"] == expected_entry_count and
            compile_record.get("schema") == COMPILE_COMMANDS_SCHEMA and
            compile_record.get("implementation") == implementation and
            compile_record.get("profile") == compile_profile and
            compile_record.get("validated_optimization") == "-O3" and
            compile_record.get("validated_openmp") is True and
            compile_record.get("isa_policy") == isa_policy and
            isinstance(sources, list) and sources and
            sources == sorted(set(sources)) and
            isinstance(pairs, list) and pairs and
            isinstance(retained_entries, list) and retained_entries and
            len(retained_entries) == len(pairs),
            f"{implementation} compile-command identity is invalid")
    entries_by_source: dict[str, dict[str, Any]] = {}
    for index, entry in enumerate(retained_entries):
        require(isinstance(entry, dict) and set(entry) == {
                    "directory", "file", "output", "arguments"} and
                isinstance(entry.get("directory"), str) and
                isinstance(entry.get("file"), str) and
                isinstance(entry.get("output"), str) and entry["output"] and
                isinstance(entry.get("arguments"), list) and
                all(isinstance(token, str) and token and "@" not in token
                    for token in entry["arguments"]) and
                entry["file"] not in entries_by_source,
                f"{implementation} retained compile entry {index} differs")
        entries_by_source[entry["file"]] = entry
    pair_sources: list[str] = []
    pair_objects: list[str] = []
    for index, pair in enumerate(pairs):
        require(isinstance(pair, dict) and set(pair) == {"source", "object"},
                f"{implementation} compile pair {index} shape differs")
        source = validate_complete_artifact_identity(
            pair["source"], f"{implementation} source {index}", "source_file")
        obj = validate_complete_artifact_identity(
            pair["object"], f"{implementation} object {index}", "object_file")
        require(Path(obj["path"]).is_relative_to(Path(build_dir)),
                f"{implementation} object {index} escapes its build directory")
        retained_entry = entries_by_source.get(source["path"])
        require(retained_entry is not None,
                f"{implementation} source {index} lacks its full compiler argv")
        expected_output = expected_compile_output(
            implementation, Path(source["path"]), specification)
        expected_arguments = expected_compile_argv(
            implementation, Path(source["path"]), specification,
            cache["CMAKE_CXX_COMPILER"])
        require(retained_entry == {
                    "directory": build_dir,
                    "file": source["path"],
                    "output": expected_output,
                    "arguments": expected_arguments,
                } and
                str(Path(build_dir) / expected_output) == obj["path"],
                f"{implementation} source {index} compiler argv/output differs")
        pair_sources.append(source["path"])
        pair_objects.append(obj["path"])
    expected_sources = sorted([
        *(str(Path(specification[f"{implementation}_source_root"]) / name)
          for name in library_sources),
        str(Path(specification["candidate_source_root"]) /
            benchmark_suffix.removeprefix("/")),
    ])
    require(pair_sources == expected_sources and
            pair_sources == sources and len(pair_sources) == len(set(pair_sources)) and
            len(pair_objects) == len(set(pair_objects)) and
            set(entries_by_source) == set(pair_sources) and
            compile_record["entry_count"] >= len(pair_sources),
            f"{implementation} compile-command source/object closure differs")
    benchmark_pairs = [pair for pair in pairs
                       if pair["source"]["path"].endswith(benchmark_suffix)]
    require(len(benchmark_pairs) == 1,
            f"{implementation} benchmark object is not unique")
    archive_pairs = [pair for pair in pairs if pair not in benchmark_pairs]
    require(archive_pairs,
            f"{implementation} archive compile closure is empty")
    candidate_source_root = Path(specification["candidate_source_root"])
    implementation_source_root = Path(
        specification[f"{implementation}_source_root"])
    require(all(Path(pair["source"]["path"]).is_relative_to(
                    candidate_source_root)
                for pair in benchmark_pairs) and
            all(Path(pair["source"]["path"]).is_relative_to(
                    implementation_source_root)
                for pair in archive_pairs),
            f"{implementation} compile sources escape their declared source roots")
    members = build.get("validated_archive_members")
    expected_members = [Path(name).name + ".o" for name in library_sources]
    require(isinstance(members, list) and members == expected_members and
            len(members) == len(set(members)),
            f"{implementation} archive member identity is invalid")
    objects_by_member: dict[str, str] = {}
    for pair in archive_pairs:
        object_path = Path(pair["object"]["path"])
        require(object_path.name not in objects_by_member,
                f"{implementation} archive member basenames are ambiguous")
        objects_by_member[object_path.name] = object_path.relative_to(
            Path(build_dir)).as_posix()
    require(set(members) == set(objects_by_member),
            f"{implementation} archive members differ from compile closure")
    require(all(archive["mtime_ns"] >= pair["object"]["mtime_ns"]
                for pair in archive_pairs),
            f"{implementation} archive predates one of its retained objects")
    require(executable["mtime_ns"] >= archive["mtime_ns"] and
            all(executable["mtime_ns"] >= pair["object"]["mtime_ns"]
                for pair in benchmark_pairs),
            f"{implementation} executable predates a retained project input")

    tools = build.get("archive_link_tool_invocations")
    require(isinstance(tools, dict) and set(tools) == {"archiver", "ranlib"},
            f"{implementation} archive tool invocation shape differs")
    for name, artifact in (("archiver", archiver), ("ranlib", ranlib)):
        invocation = tools.get(name)
        require(isinstance(invocation, dict) and set(invocation) == {
                    "invocation", "resolved_path"} and
                isinstance(invocation.get("invocation"), str) and
                invocation["invocation"] and
                invocation.get("resolved_path") == artifact["path"],
                f"{implementation} {name} invocation identity differs")
    validate_archive_link_recipe_content(
        build.get("archive_link_recipe_content"),
        build["archive_link_recipe"], archive_name, archive_target_directory,
        f"{implementation} archive link recipe",
        expected_objects=[objects_by_member[member] for member in members],
        expected_archiver=tools["archiver"]["invocation"],
        expected_ranlib=tools["ranlib"]["invocation"])
    validate_complete_executable_recipe(
        build.get("executable_link_recipe_content"),
        build["executable_link_recipe"], cache["CMAKE_CXX_COMPILER"],
        archive_name, executable_name,
        Path(benchmark_pairs[0]["object"]["path"]).relative_to(
            Path(build_dir)).as_posix(),
        f"{implementation} executable link recipe",
        external_link_inputs)
    return build


def validate_complete_input_snapshot(
    specification: Mapping[str, Any], snapshot: object, raw_schema: str = RAW_SCHEMA,
) -> dict[str, Any]:
    require(raw_schema in COMPLETE_EVIDENCE_SCHEMAS,
            "complete input snapshot uses an unsupported schema")
    require(isinstance(snapshot, dict) and set(snapshot) == {
                "runner", "taskset", "ldd", "baseline_executable",
                "candidate_executable", "baseline_archive", "candidate_archive",
                "baseline_build", "candidate_build", "baseline_runtime_closure",
                "candidate_runtime_closure", "baseline_source", "candidate_source"},
            "complete input identity shape differs")
    for name, kind in (("runner", "file"), ("taskset", "executable"),
                       ("ldd", "executable")):
        record = validate_complete_artifact_identity(
            snapshot.get(name), name, kind)
        require(record["path"] == specification[name],
                f"{name} identity differs from the input specification")
    for role in ("baseline", "candidate"):
        executable = validate_complete_artifact_identity(
            snapshot.get(f"{role}_executable"), f"{role} executable", "executable")
        archive = validate_complete_artifact_identity(
            snapshot.get(f"{role}_archive"), f"{role} archive", "archive")
        require(executable["path"] == specification[f"{role}_executable"] and
                archive["path"] == specification[f"{role}_archive"],
                f"{role} output identity differs from the input specification")
        build = validate_complete_build_identity(
            snapshot.get(f"{role}_build"), role, specification)
        require(executable == build["validated_executable"] and
                archive == build["validated_archive"],
                f"{role} top-level/build output identity differs")
        validate_complete_runtime_closure(
            snapshot.get(f"{role}_runtime_closure"), role,
            executable["path"], raw_schema)
    validate_complete_git_identity(
        snapshot.get("baseline_source"), "baseline",
        specification["baseline_source_root"], MAIN_COMMIT,
        require_detached=True)
    validate_complete_git_identity(
        snapshot.get("candidate_source"), "candidate",
        specification["candidate_source_root"], specification["candidate_commit"],
        require_detached=False)
    require(snapshot["baseline_build"]["compiler"] ==
                snapshot["candidate_build"]["compiler"] and
            snapshot["baseline_build"]["compiler_version_stdout"] ==
                snapshot["candidate_build"]["compiler_version_stdout"],
            "baseline and candidate compiler identities differ")
    require(snapshot["baseline_build"]["validated_external_link_inputs"] ==
                snapshot["candidate_build"]["validated_external_link_inputs"],
            "baseline and candidate external link-input identities differ")
    return snapshot


def validate_candidate_cmake_identity(
    specification: Mapping[str, Any], snapshot: object, raw_schema: str,
) -> None:
    """Bind portable evidence to the schema's exact CMake identity.

    Current evidence also retains and parses exact bounded recipe bytes.
    Historical replay intentionally keeps its original record shape and does
    not require old build paths to remain present on the current machine.
    """
    identity = cmake_identity_for_raw_schema(raw_schema)
    require(isinstance(snapshot, dict), "input identity is not an object")
    archive = snapshot.get("candidate_archive")
    build = snapshot.get("candidate_build")
    require(isinstance(archive, dict) and isinstance(build, dict),
            "candidate archive/build identity is missing")
    validated_archive = build.get("validated_archive")
    archive_recipe = build.get("archive_link_recipe")
    require(isinstance(validated_archive, dict) and
            isinstance(archive_recipe, dict),
            "candidate CMake archive provenance is incomplete")
    expected_archive = identity["archive"]
    expected_recipe_suffix = "/CMakeFiles/{}/link.txt".format(
        identity["target_directory"])
    paths = (
        specification.get("candidate_archive"),
        archive.get("path"),
        validated_archive.get("path"),
    )
    require(all(isinstance(path, str) and Path(path).name == expected_archive
                for path in paths),
            "candidate archive name differs from its evidence schema")
    recipe_path = archive_recipe.get("path")
    require(isinstance(recipe_path, str) and
            recipe_path.replace("\\", "/").endswith(expected_recipe_suffix),
            "candidate CMake target directory differs from its evidence schema")
    require(archive == validated_archive,
            "candidate archive identity differs from validated build artifact")
    if raw_schema in HARDENED_BUILD_SCHEMAS:
        compiler_records = build.get("validated_compile_commands")
        members = build.get("validated_archive_members")
        build_dir = build.get("build_dir")
        archiver = build.get("archiver")
        ranlib = build.get("ranlib")
        tool_invocations = build.get("archive_link_tool_invocations")
        require(isinstance(compiler_records, dict) and
                isinstance(compiler_records.get("required_source_object_pairs"), list) and
                isinstance(members, list) and members and
                len(members) == len(set(members)) and
                all(isinstance(member, str) and member for member in members) and
                isinstance(build_dir, str) and build_dir and
                isinstance(archiver, dict) and isinstance(archiver.get("path"), str) and
                isinstance(ranlib, dict) and isinstance(ranlib.get("path"), str) and
                isinstance(tool_invocations, dict) and
                set(tool_invocations) == {"archiver", "ranlib"} and
                all(isinstance(tool_invocations.get(role), dict) and
                    set(tool_invocations[role]) == {"invocation", "resolved_path"} and
                    isinstance(tool_invocations[role]["invocation"], str) and
                    tool_invocations[role]["invocation"] and
                    tool_invocations[role]["resolved_path"] ==
                        (archiver if role == "archiver" else ranlib)["path"]
                    for role in ("archiver", "ranlib")),
                "candidate archive command closure is incomplete")
        objects_by_member: dict[str, str] = {}
        for record in compiler_records["required_source_object_pairs"]:
            require(isinstance(record, dict) and
                    isinstance(record.get("source"), dict) and
                    isinstance(record.get("object"), dict) and
                    isinstance(record["source"].get("path"), str) and
                    isinstance(record["object"].get("path"), str),
                    "candidate compile-command closure is malformed")
            if record["source"]["path"].endswith("/bench/leopard2/benchmark.cpp"):
                continue
            object_path = Path(record["object"]["path"])
            try:
                relative_object = object_path.relative_to(Path(build_dir)).as_posix()
            except ValueError as error:
                raise EvidenceError(
                    "candidate archive object escapes its build directory") from error
            require(object_path.name not in objects_by_member,
                    "candidate archive object basenames are ambiguous")
            objects_by_member[object_path.name] = relative_object
        require(set(objects_by_member) == set(members),
                "candidate archive members differ from compile-command closure")
        validate_archive_link_recipe_content(
            build.get("archive_link_recipe_content"), archive_recipe,
            expected_archive, identity["target_directory"],
            "candidate archive link recipe",
            expected_objects=[objects_by_member[member] for member in members],
            expected_archiver=tool_invocations["archiver"]["invocation"],
            expected_ranlib=tool_invocations["ranlib"]["invocation"])
    else:
        require("archive_link_recipe_content" not in build and
                "ranlib" not in build and
                "archive_link_tool_invocations" not in build,
                "historical build identity has a current-only recipe-content field")


def parse_cpu_list(text: str) -> set[int]:
    require(isinstance(text, str), "CPU list is not text")
    try:
        encoded = text.encode("ascii", errors="strict")
    except UnicodeEncodeError as error:
        raise EvidenceError("CPU list is not ASCII text") from error
    require(len(encoded) <= MAX_CPU_LIST_TEXT_BYTES,
            "CPU list is not bounded ASCII text")
    max_digits = len(str(MAX_CPU_ID))
    result: set[int] = set()
    for component in text.strip().split(","):
        if not component:
            continue
        if "-" in component:
            require(re.fullmatch(r"[0-9]+-[0-9]+", component) is not None,
                    f"invalid CPU range {component!r}")
            bounds = component.split("-", 1)
            require(all(len(item) <= max_digits for item in bounds),
                    f"CPU range is out of bounds: {component!r}")
            first, last = (int(item) for item in bounds)
            require(0 <= first <= last <= MAX_CPU_ID and
                    last - first + 1 <= MAX_CPU_LIST_ENTRIES,
                    f"CPU range is out of bounds: {component!r}")
            result.update(range(first, last + 1))
        else:
            require(re.fullmatch(r"[0-9]+", component) is not None,
                    f"invalid CPU identity {component!r}")
            require(len(component) <= max_digits,
                    f"CPU identity is out of bounds: {component!r}")
            value = int(component)
            require(0 <= value <= MAX_CPU_ID,
                    f"CPU identity is out of bounds: {component!r}")
            result.add(value)
        require(len(result) <= MAX_CPU_LIST_ENTRIES,
                "CPU list contains too many entries")
    return result


def validate_topology(cpu: int, sibling: int) -> tuple[set[int], set[int]]:
    require(type(cpu) is int and type(sibling) is int and
            0 <= cpu <= MAX_CPU_ID and 0 <= sibling <= MAX_CPU_ID and
            cpu != sibling,
            "benchmark CPU and reserved sibling must be distinct non-negative CPUs")
    require(hasattr(os, "sched_getaffinity") and hasattr(os, "sched_setaffinity"),
            "Linux scheduling affinity is required")
    allowed = set(os.sched_getaffinity(0))
    require(cpu in allowed and sibling in allowed,
            f"CPU pair {cpu}/{sibling} is not in initial affinity {sorted(allowed)}; "
            "launch the runner with an affinity containing both CPUs")
    topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
    require(topology.is_file(), f"missing topology file {topology}")
    siblings = parse_cpu_list(topology.read_text(encoding="ascii"))
    require(siblings == {cpu, sibling},
            f"physical core must have exactly the reserved pair {cpu}/{sibling}; "
            f"topology reports {sorted(siblings)}")
    housekeeping = allowed - {cpu, sibling}
    require(housekeeping, "no housekeeping CPU remains after reserving the physical core")
    return allowed, housekeeping


def read_text_optional(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        return None


def cpuinfo_identity(cpu: int) -> dict[str, str]:
    blocks = Path("/proc/cpuinfo").read_text(encoding="utf-8").split("\n\n")
    selected: dict[str, str] | None = None
    for block in blocks:
        values: dict[str, str] = {}
        for line in block.splitlines():
            if ":" in line:
                name, value = line.split(":", 1)
                values[name.strip()] = value.strip()
        if values.get("processor") == str(cpu):
            selected = values
            break
    require(selected is not None, f"CPU {cpu} is absent from /proc/cpuinfo")
    retained = {
        name: selected[name] for name in (
            "processor", "vendor_id", "cpu family", "model", "model name",
            "stepping", "microcode", "flags", "Features", "CPU implementer",
            "CPU architecture", "CPU variant", "CPU part", "CPU revision")
        if name in selected
    }
    require(any(name in retained for name in ("model name", "CPU part")),
            f"CPU {cpu} has no retained model identity")
    return retained


def sysfs_numeric_directory_inventory(
    root: Path, prefix: str, label: str,
) -> dict[str, Any]:
    """Retain a canonical raw directory listing independent of summaries."""
    require(root.is_dir(), f"{label} sysfs directory is absent")
    names = [name for name in os.listdir(root) if name.startswith(prefix)]
    require(names and all(name.removeprefix(prefix).isdigit() for name in names),
            f"{label} sysfs directory has an invalid indexed entry")
    names.sort(key=lambda name: int(name.removeprefix(prefix)))
    text = "".join(name + "\n" for name in names)
    return exact_text_content(text, f"{label} sysfs directory inventory")


def validate_sysfs_numeric_directory_inventory(
    value: object, prefix: str, label: str,
) -> list[str]:
    retained = validate_complete_text_identity(value, label)
    text = retained["text"]
    names = text.splitlines()
    require(names and text == "".join(name + "\n" for name in names) and
            all(name.startswith(prefix) and
                name.removeprefix(prefix).isdigit() for name in names) and
            names == sorted(set(names),
                            key=lambda name: int(name.removeprefix(prefix))),
            f"{label} is not a canonical numeric directory inventory")
    return names


def cpu_policy_identity(cpu: int) -> dict[str, Any]:
    root = Path(f"/sys/devices/system/cpu/cpu{cpu}")
    topology_root = root / "topology"
    topology = {
        name: read_text_optional(topology_root / name) for name in (
            "core_id", "physical_package_id", "die_id", "cluster_id",
            "thread_siblings_list", "core_siblings_list")
    }
    require(topology["thread_siblings_list"] is not None,
            f"CPU {cpu} has no thread sibling topology")
    frequency_root = root / "cpufreq"
    frequency = {
        name: read_text_optional(frequency_root / name) for name in (
            "scaling_driver", "scaling_governor", "energy_performance_preference",
            "scaling_min_freq", "scaling_max_freq", "cpuinfo_min_freq",
            "cpuinfo_max_freq")
    }
    cache_fields = (
        "level", "type", "size", "coherency_line_size", "number_of_sets",
        "ways_of_associativity", "physical_line_partition",
        "shared_cpu_list", "shared_cpu_map", "allocation_policy",
        "write_policy",
    )
    cache_directory_inventory_text = sysfs_numeric_directory_inventory(
        root / "cache", "index", f"CPU {cpu} cache")
    cache_roots = list((root / "cache").glob("index*"))
    require(all(path.name.removeprefix("index").isdigit()
                for path in cache_roots),
            f"CPU {cpu} has an invalid cache-index inventory")
    cache_roots.sort(
        key=lambda path: int(path.name.removeprefix("index")))
    cache_index_inventory = [path.name for path in cache_roots]
    require(cache_index_inventory ==
                cache_directory_inventory_text["text"].splitlines(),
            f"CPU {cpu} cache directory changed during identity capture")
    caches = []
    for index_root in cache_roots:
        suffix = index_root.name.removeprefix("index")
        require(suffix.isdigit(), f"CPU {cpu} has an invalid cache index")
        record = {
            "index": int(suffix),
            **{name: read_text_optional(index_root / name)
               for name in cache_fields},
        }
        require(all(record[name] is not None for name in (
                    "level", "type", "size", "coherency_line_size",
                    "shared_cpu_list", "shared_cpu_map")),
                f"CPU {cpu} cache index {suffix} lacks domain identity")
        require(cpu in parse_cpu_list(str(record["shared_cpu_list"])),
                f"CPU {cpu} cache index {suffix} excludes its own CPU")
        caches.append(record)
    require(caches, f"CPU {cpu} has no retained cache hierarchy")
    numa_node_paths = [
        path for path in root.glob("node*")
        if path.name.removeprefix("node").isdigit()]
    numa_node_paths.sort(
        key=lambda path: int(path.name.removeprefix("node")))
    numa_node_inventory = [path.name for path in numa_node_paths]
    numa_nodes = [int(name.removeprefix("node"))
                  for name in numa_node_inventory]
    core_class = {
        "core_type": read_text_optional(topology_root / "core_type"),
        "cpu_capacity": read_text_optional(root / "cpu_capacity"),
    }
    return {
        "cpu": cpu,
        "online": read_text_optional(root / "online"),
        "cpuinfo": cpuinfo_identity(cpu),
        "topology": topology,
        "frequency_policy": frequency,
        "cache_hierarchy": caches,
        "cache_index_inventory": cache_index_inventory,
        "cache_directory_inventory_text": cache_directory_inventory_text,
        "numa_nodes": numa_nodes,
        "numa_node_inventory": numa_node_inventory,
        "core_class": core_class,
    }


def host_identity(cpu: int, sibling: int, allowed_at_launch: set[int]) -> dict[str, Any]:
    require(cpu in allowed_at_launch and sibling in allowed_at_launch,
            "host identity launch CPU set does not contain the reserved pair")
    turbo_paths = (
        Path("/sys/devices/system/cpu/intel_pstate/no_turbo"),
        Path("/sys/devices/system/cpu/amd_pstate/status"),
        Path("/sys/devices/system/cpu/cpufreq/boost"),
    )
    uname = platform.uname()
    online_path = Path("/sys/devices/system/cpu/online")
    node_online_path = Path("/sys/devices/system/node/online")
    online = online_path.read_text(encoding="ascii")
    node_online = node_online_path.read_text(encoding="ascii")
    return {
        "system": {
            "system": uname.system,
            "node": uname.node,
            "release": uname.release,
            "version": uname.version,
            "machine": uname.machine,
            "python": platform.python_version(),
            "page_size": os.sysconf("SC_PAGE_SIZE"),
        },
        "allowed_cpu_set_at_launch": sorted(allowed_at_launch),
        "online_cpu_set": sorted(parse_cpu_list(online)),
        "online_cpu_list_text": exact_text_content(
            online, "sysfs online CPU list"),
        "online_node_list_text": exact_text_content(
            node_online, "sysfs online NUMA-node list"),
        "benchmark_cpu": cpu_policy_identity(cpu),
        "reserved_sibling": cpu_policy_identity(sibling),
        "turbo_and_pstate": {
            str(path): read_text_optional(path) for path in turbo_paths
        },
    }


def cpu_stat_snapshot(cpu: int) -> dict[str, Any]:
    """Read non-double-counted Linux scheduler counters for one logical CPU."""
    require(cpu >= 0, "CPU stat identity must be non-negative")
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if not line.startswith(prefix):
            continue
        tokens = line.split()
        require(tokens[0] == f"cpu{cpu}" and
                len(tokens) >= 1 + len(CPU_STAT_FIELDS),
                f"CPU {cpu} has an incomplete /proc/stat record")
        try:
            values = [int(value) for value in tokens[1:1 + len(CPU_STAT_FIELDS)]]
        except ValueError as error:
            raise EvidenceError(f"CPU {cpu} has a non-integer /proc/stat record") from error
        require(all(value >= 0 for value in values),
                f"CPU {cpu} has a negative /proc/stat counter")
        fields = dict(zip(CPU_STAT_FIELDS, values))
        return {"cpu": cpu, "fields": fields, "total_jiffies": sum(values)}
    raise EvidenceError(f"CPU {cpu} is absent from /proc/stat")


def cpu_stat_delta(before: Mapping[str, Any], after: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(before, dict) and isinstance(after, dict) and
            before.get("cpu") == after.get("cpu") and
            isinstance(before.get("cpu"), int) and
            not isinstance(before.get("cpu"), bool),
            "CPU stat snapshots refer to different CPUs")
    before_fields = before.get("fields")
    after_fields = after.get("fields")
    require(isinstance(before_fields, dict) and isinstance(after_fields, dict) and
            set(before_fields) == set(CPU_STAT_FIELDS) and
            set(after_fields) == set(CPU_STAT_FIELDS),
            "CPU stat snapshot fields are incomplete")
    deltas: dict[str, int] = {}
    for name in CPU_STAT_FIELDS:
        first = before_fields[name]
        last = after_fields[name]
        require(isinstance(first, int) and not isinstance(first, bool) and first >= 0 and
                isinstance(last, int) and not isinstance(last, bool) and last >= first,
                f"CPU stat counter {name} moved backwards")
        deltas[name] = last - first
    idle = sum(deltas[name] for name in CPU_STAT_IDLE_FIELDS)
    nonidle = sum(value for name, value in deltas.items()
                  if name not in CPU_STAT_IDLE_FIELDS)
    total = idle + nonidle
    require(after.get("total_jiffies") == sum(after_fields.values()) and
            before.get("total_jiffies") == sum(before_fields.values()),
            "CPU stat total does not match its fields")
    return {
        "cpu": before["cpu"],
        "fields": deltas,
        "idle_jiffies": idle,
        "nonidle_jiffies": nonidle,
        "total_jiffies": total,
    }


def pair_lease_payload(
    cpu: int, sibling: int, uid: int | None = None
) -> dict[str, Any]:
    retained_uid = os.getuid() if uid is None else uid
    return {
        "cpus": sorted((cpu, sibling)),
        "schema": PAIR_LEASE_SCHEMA,
        "uid": retained_uid,
    }


def pair_lease_name(cpu: int, sibling: int, uid: int | None = None) -> str:
    first, second = sorted((cpu, sibling))
    retained_uid = os.getuid() if uid is None else uid
    return f"leopard2-cpu-pair-{retained_uid}-{first}-{second}.lock"


def pair_lease_directory(uid: int | None = None) -> Path:
    retained_uid = os.getuid() if uid is None else uid
    return Path("/run/user") / str(retained_uid) / "leopard2-cpu-leases"


class PairLease:
    """Serialize normal Leopard2 evidence runners by physical CPU pair."""

    def __init__(self, cpu: int, sibling: int, root: Path | None = None):
        require(cpu >= 0 and sibling >= 0 and cpu != sibling,
                "pair lease requires two distinct non-negative CPUs")
        self.cpu = cpu
        self.sibling = sibling
        self.production_root = root is None
        self.root = pair_lease_directory() if root is None else root
        self.path = self.root / \
            pair_lease_name(cpu, sibling)
        self.descriptor: int | None = None
        self.identity: dict[str, Any] | None = None
        self.kernel_socket: socket.socket | None = None
        material = canonical_bytes({
            "cpus": sorted((cpu, sibling)),
            "root": os.path.abspath(self.root),
            "schema": PAIR_LEASE_SCHEMA,
            "uid": os.getuid(),
        })
        self.kernel_name = b"\0leopard2-pair-v1-" + \
            hashlib.sha256(material).hexdigest()[:40].encode("ascii")

    def _acquire_kernel_lease(self) -> None:
        require(sys.platform.startswith("linux") and hasattr(socket, "AF_UNIX"),
                "Linux abstract Unix sockets are required for stable CPU leases")
        lease = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        lease.set_inheritable(False)
        try:
            lease.bind(self.kernel_name)
        except OSError as error:
            lease.close()
            if error.errno == errno.EADDRINUSE:
                raise EvidenceError(
                    "physical CPU pair already has a kernel lease") from error
            raise EvidenceError(f"cannot bind stable CPU pair lease: {error}") from error
        self.kernel_socket = lease

    def _release_kernel_lease(self) -> None:
        if self.kernel_socket is not None:
            self.kernel_socket.close()
            self.kernel_socket = None

    def _validate_directory(self) -> os.stat_result:
        if self.production_root:
            runtime = os.lstat(self.root.parent)
            require(stat.S_ISDIR(runtime.st_mode) and
                    runtime.st_uid == os.getuid() and
                    stat.S_IMODE(runtime.st_mode) == 0o700,
                    "CPU pair runtime directory is not an owned mode-0700 directory")
        metadata = os.lstat(self.root)
        require(stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o700,
                "CPU pair lease directory is not an owned mode-0700 directory")
        return metadata

    def validate_current(self) -> None:
        require(self.descriptor is not None and self.identity is not None and
                self.kernel_socket is not None and
                self.kernel_socket.fileno() >= 0 and
                self.kernel_socket.getsockname() == self.kernel_name,
                "CPU pair lease is not held")
        try:
            directory = self._validate_directory()
            descriptor = os.fstat(self.descriptor)
            path = os.lstat(self.path)
        except OSError as error:
            raise EvidenceError(
                f"CPU pair lease path/descriptor revalidation failed: {error}") from error
        require(stat.S_ISREG(descriptor.st_mode) and
                descriptor.st_uid == os.getuid() and descriptor.st_nlink == 1 and
                stat.S_IMODE(descriptor.st_mode) == 0o600 and
                (descriptor.st_dev, descriptor.st_ino) == (path.st_dev, path.st_ino) and
                (descriptor.st_dev, descriptor.st_ino) ==
                    (self.identity["device"], self.identity["inode"]),
                "CPU pair lease path was replaced or its metadata changed")
        require((directory.st_dev, directory.st_ino) ==
                (self.identity["directory_device"],
                 self.identity["directory_inode"]),
                "CPU pair lease directory was replaced")
        expected = canonical_bytes(pair_lease_payload(self.cpu, self.sibling))
        os.lseek(self.descriptor, 0, os.SEEK_SET)
        retained = os.read(self.descriptor, 4096)
        require(retained == expected,
                "CPU pair lease contents changed while held")

    def __enter__(self) -> dict[str, Any]:
        self._acquire_kernel_lease()
        created_directory = False
        try:
            self.root.mkdir(mode=0o700)
            created_directory = True
        except FileExistsError:
            pass
        except OSError as error:
            self._release_kernel_lease()
            raise EvidenceError(f"cannot create CPU pair lease directory: {error}") from error
        try:
            if created_directory:
                os.chmod(self.root, 0o700)
            self._validate_directory()
        except Exception:
            self._release_kernel_lease()
            raise
        flags = os.O_RDWR
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        if hasattr(os, "O_NONBLOCK"):
            flags |= os.O_NONBLOCK
        created_file = False
        before: os.stat_result | None = None
        try:
            try:
                self.descriptor = os.open(
                    self.path, flags | os.O_CREAT | os.O_EXCL, 0o600)
                created_file = True
            except FileExistsError:
                before = os.lstat(self.path)
                require(stat.S_ISREG(before.st_mode),
                        "CPU pair lease path is not a regular file")
                self.descriptor = os.open(self.path, flags)
            if created_file:
                os.fchmod(self.descriptor, 0o600)
        except OSError as error:
            if self.descriptor is not None:
                os.close(self.descriptor)
                self.descriptor = None
            self._release_kernel_lease()
            raise EvidenceError(f"cannot open CPU pair lease: {error}") from error
        try:
            metadata = os.fstat(self.descriptor)
            path_metadata = os.lstat(self.path)
            require(stat.S_ISREG(metadata.st_mode) and metadata.st_uid == os.getuid() and
                    metadata.st_nlink == 1 and stat.S_IMODE(metadata.st_mode) == 0o600 and
                    (metadata.st_dev, metadata.st_ino) ==
                    (path_metadata.st_dev, path_metadata.st_ino) and
                    (before is None or (metadata.st_dev, metadata.st_ino) ==
                     (before.st_dev, before.st_ino)),
                    "CPU pair lease file has unsafe ownership, type, links, or permissions")
            try:
                fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise EvidenceError(
                    f"physical CPU pair is already leased: {self.path}") from error
            expected = canonical_bytes(pair_lease_payload(self.cpu, self.sibling))
            os.lseek(self.descriptor, 0, os.SEEK_SET)
            retained = os.read(self.descriptor, 4096)
            if not retained:
                written = os.write(self.descriptor, expected)
                require(written == len(expected), "short write to CPU pair lease")
                os.fsync(self.descriptor)
                retained = expected
            require(retained == expected,
                    "CPU pair lease has unexpected or noncanonical contents")
            directory = self._validate_directory()
            self.identity = {
                "device": metadata.st_dev,
                "directory_device": directory.st_dev,
                "directory_inode": directory.st_ino,
                "inode": metadata.st_ino,
                "lock": "exclusive_nonblocking_pair_wide",
                "path": str(self.path.resolve()),
                "payload": pair_lease_payload(self.cpu, self.sibling),
                "sha256": sha256_bytes(retained),
            }
            self.validate_current()
            return self.identity
        except Exception:
            if self.descriptor is not None:
                try:
                    fcntl.flock(self.descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(self.descriptor)
                    self.descriptor = None
            self._release_kernel_lease()
            raise

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        descriptor = self.descriptor
        self.descriptor = None
        try:
            if descriptor is not None:
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(descriptor)
        finally:
            self._release_kernel_lease()


def validate_pair_lease_identity(
    value: object, cpu: int, sibling: int
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "device", "directory_device", "directory_inode", "inode", "lock",
        "path", "payload", "sha256"},
        "CPU pair lease identity is incomplete")
    payload = value.get("payload")
    require(isinstance(payload, dict) and set(payload) == {"cpus", "schema", "uid"} and
            payload.get("schema") == PAIR_LEASE_SCHEMA and
            isinstance(payload.get("cpus"), list) and
            all(isinstance(item, int) and not isinstance(item, bool) and
                0 <= item <= MAX_CPU_ID for item in payload["cpus"]) and
            payload.get("cpus") == sorted((cpu, sibling)) and
            isinstance(payload.get("uid"), int) and
            not isinstance(payload.get("uid"), bool) and payload["uid"] >= 0,
            "CPU pair lease payload does not match the campaign")
    expected_name = pair_lease_name(cpu, sibling, payload["uid"])
    path = value.get("path")
    expected_path = pair_lease_directory(payload["uid"]) / expected_name
    require(isinstance(path, str) and Path(path) == expected_path,
            "CPU pair lease path does not identify the reserved pair")
    require(all(isinstance(value.get(name), int) and
                not isinstance(value.get(name), bool) and value[name] >= 0
                for name in ("device", "directory_device", "directory_inode", "inode")),
            "CPU pair lease filesystem identity is invalid")
    require(value.get("lock") == "exclusive_nonblocking_pair_wide" and
            value.get("sha256") == sha256_bytes(canonical_bytes(payload)),
            "CPU pair lease lock or digest is invalid")
    return value


def isolation_record(
    cpu: int,
    sibling: int,
    pair_lease: Mapping[str, Any],
    before_monotonic_ns: int,
    after_monotonic_ns: int,
    before_cpu: Mapping[str, Any],
    after_cpu: Mapping[str, Any],
    before_sibling: Mapping[str, Any],
    after_sibling: Mapping[str, Any],
) -> dict[str, Any]:
    require(all(isinstance(snapshot, dict) for snapshot in (
                before_cpu, after_cpu, before_sibling, after_sibling)),
            "isolation snapshots are not objects")
    require(isinstance(before_monotonic_ns, int) and
            not isinstance(before_monotonic_ns, bool) and before_monotonic_ns >= 0 and
            isinstance(after_monotonic_ns, int) and
            not isinstance(after_monotonic_ns, bool) and after_monotonic_ns >= 0,
            "isolation monotonic timestamps are invalid")
    require(
        before_cpu.get("cpu") == cpu and after_cpu.get("cpu") == cpu and
        before_sibling.get("cpu") == sibling and
        after_sibling.get("cpu") == sibling,
        "isolation snapshots do not match the reserved CPU pair")
    benchmark_delta = cpu_stat_delta(before_cpu, after_cpu)
    sibling_delta = cpu_stat_delta(before_sibling, after_sibling)
    require(benchmark_delta["cpu"] == cpu and sibling_delta["cpu"] == sibling,
            "isolation deltas do not match the reserved CPU pair")
    accepted = (
        after_monotonic_ns > before_monotonic_ns and
        benchmark_delta["nonidle_jiffies"] > 0 and
        sibling_delta["total_jiffies"] > 0 and
        sibling_delta["nonidle_jiffies"] <= MAX_SIBLING_NONIDLE_JIFFIES
    )
    return {
        "accepted": accepted,
        "after": {
            "benchmark_cpu": dict(after_cpu),
            "monotonic_ns": after_monotonic_ns,
            "reserved_sibling": dict(after_sibling),
        },
        "before": {
            "benchmark_cpu": dict(before_cpu),
            "monotonic_ns": before_monotonic_ns,
            "reserved_sibling": dict(before_sibling),
        },
        "benchmark_cpu": cpu,
        "delta": {
            "benchmark_cpu": benchmark_delta,
            "reserved_sibling": sibling_delta,
        },
        "pair_lease": dict(pair_lease),
        "policy": {
            "counter_source": "/proc/stat",
            "idle_fields": list(CPU_STAT_IDLE_FIELDS),
            "nonidle_fields": [
                name for name in CPU_STAT_FIELDS if name not in CPU_STAT_IDLE_FIELDS
            ],
            "reserved_sibling_max_nonidle_jiffies":
                MAX_SIBLING_NONIDLE_JIFFIES,
        },
        "reserved_sibling": sibling,
        "schema": ISOLATION_SCHEMA,
    }


def validate_isolation(
    value: object, cpu: int, sibling: int, require_accepted: bool = True
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "accepted", "after", "before", "benchmark_cpu", "delta", "pair_lease",
        "policy", "reserved_sibling", "schema"},
        "isolation evidence is incomplete")
    require(value.get("schema") == ISOLATION_SCHEMA and
            isinstance(value.get("benchmark_cpu"), int) and
            not isinstance(value.get("benchmark_cpu"), bool) and
            isinstance(value.get("reserved_sibling"), int) and
            not isinstance(value.get("reserved_sibling"), bool) and
            value.get("benchmark_cpu") == cpu and
            value.get("reserved_sibling") == sibling,
            "isolation evidence CPU identity is invalid")
    require(value.get("policy") == {
        "counter_source": "/proc/stat",
        "idle_fields": list(CPU_STAT_IDLE_FIELDS),
        "nonidle_fields": [
            name for name in CPU_STAT_FIELDS if name not in CPU_STAT_IDLE_FIELDS
        ],
        "reserved_sibling_max_nonidle_jiffies": MAX_SIBLING_NONIDLE_JIFFIES,
    }, "isolation evidence policy was edited")
    validate_pair_lease_identity(value.get("pair_lease"), cpu, sibling)
    before = value.get("before")
    after = value.get("after")
    require(isinstance(before, dict) and isinstance(after, dict) and
            set(before) == {"benchmark_cpu", "monotonic_ns", "reserved_sibling"} and
            set(after) == {"benchmark_cpu", "monotonic_ns", "reserved_sibling"} and
            isinstance(before.get("monotonic_ns"), int) and
            not isinstance(before.get("monotonic_ns"), bool) and
            isinstance(after.get("monotonic_ns"), int) and
            not isinstance(after.get("monotonic_ns"), bool) and
            all(isinstance(record, dict) for record in (
                before.get("benchmark_cpu"), before.get("reserved_sibling"),
                after.get("benchmark_cpu"), after.get("reserved_sibling"))),
            "isolation snapshots are incomplete")
    expected = isolation_record(
        cpu, sibling, value["pair_lease"], before["monotonic_ns"],
        after["monotonic_ns"], before["benchmark_cpu"], after["benchmark_cpu"],
        before["reserved_sibling"], after["reserved_sibling"])
    require(value == expected, "isolation deltas or acceptance were edited")
    if require_accepted:
        require(value["accepted"] is True,
                "reserved SMT sibling performed non-idle work during the campaign")
    return value


def supervision_record(
    nonce: str,
    runner_started_ns: int,
    runner_finished_ns: int,
    campaign: Mapping[str, Any],
    reservation: Mapping[str, Any],
    isolation: Mapping[str, Any],
) -> dict[str, Any]:
    require(HEX256.fullmatch(nonce) is not None,
            "affinity execution nonce is not 256-bit lowercase hex")
    require(type(runner_started_ns) is int and type(runner_finished_ns) is int and
            0 <= runner_started_ns <= runner_finished_ns,
            "runner supervision interval is invalid")
    payload = reservation.get("payload")
    require(isinstance(payload, dict), "supervision reservation payload is missing")
    before = isolation.get("before")
    after = isolation.get("after")
    require(isinstance(before, dict) and isinstance(after, dict),
            "supervision isolation interval is missing")
    return {
        "schema": SUPERVISION_SCHEMA,
        "execution_nonce": nonce,
        "runner_pid": os.getpid(),
        "runner_started_monotonic_ns": runner_started_ns,
        "runner_finished_monotonic_ns": runner_finished_ns,
        "launch_cpus": list(campaign["allowed_cpu_set_at_launch"]),
        "reserved_cpus": sorted((campaign["benchmark_cpu"],
                                 campaign["reserved_sibling"])),
        "campaign_sha256": sha256_bytes(canonical_bytes(campaign)),
        "reservation_sha256": reservation["sha256"],
        "reservation_nonce": payload["nonce"],
        "isolation_before_monotonic_ns": before["monotonic_ns"],
        "isolation_after_monotonic_ns": after["monotonic_ns"],
    }


def validate_supervision(
    value: object,
    campaign: Mapping[str, Any],
    reservation: Mapping[str, Any],
    isolation: Mapping[str, Any],
) -> dict[str, Any]:
    expected_keys = {
        "campaign_sha256", "execution_nonce",
        "isolation_after_monotonic_ns", "isolation_before_monotonic_ns",
        "launch_cpus", "reservation_nonce", "reservation_sha256",
        "reserved_cpus", "runner_finished_monotonic_ns", "runner_pid",
        "runner_started_monotonic_ns", "schema",
    }
    require(isinstance(value, dict) and set(value) == expected_keys,
            "supervision handshake has unexpected or missing fields")
    require(value.get("schema") == SUPERVISION_SCHEMA and
            isinstance(value.get("execution_nonce"), str) and
            HEX256.fullmatch(value["execution_nonce"]) is not None and
            type(value.get("runner_pid")) is int and value["runner_pid"] > 0,
            "supervision handshake identity is invalid")
    expected = supervision_record(
        value["execution_nonce"], value["runner_started_monotonic_ns"],
        value["runner_finished_monotonic_ns"], campaign, reservation, isolation)
    expected["runner_pid"] = value["runner_pid"]
    require(value == expected, "supervision handshake semantics were edited")
    require(value["runner_started_monotonic_ns"] <=
            value["isolation_before_monotonic_ns"] <=
            value["isolation_after_monotonic_ns"] <=
            value["runner_finished_monotonic_ns"],
            "supervision interval does not enclose scheduler evidence")
    return value


def pair_lease_runtime_root(uid: int | None = None) -> Path:
    """Return the user-owned, root-anchored runtime directory shared by runners."""
    retained_uid = os.getuid() if uid is None else uid
    return Path("/run/user") / str(retained_uid)


class StableLeaseAnchor:
    """Serialize current Leopard2 evidence runners across replaceable files."""

    def __init__(self, path: Path | None = None):
        self.path = (pair_lease_runtime_root() if path is None else
                     Path(os.path.abspath(os.fspath(path))))
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


class Reservation:
    """Hold the coordinator-created canonical reservation for the whole run."""

    def __init__(self, path: Path, cpu: int, sibling: int,
                 runtime_root: Path | None = None):
        self.requested_path = Path(os.path.abspath(os.fspath(path)))
        self.path = self.requested_path
        self.cpu = cpu
        self.sibling = sibling
        self.descriptor: int | None = None
        self.identity: dict[str, Any] | None = None
        self.anchor = StableLeaseAnchor(runtime_root)

    def __enter__(self) -> dict[str, Any]:
        try:
            self.anchor.__enter__()
            self.path = self.requested_path.resolve(strict=True)
            flags = os.O_RDONLY
            if hasattr(os, "O_NOFOLLOW"):
                flags |= os.O_NOFOLLOW
            self.descriptor = os.open(self.path, flags)
            try:
                fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise EvidenceError(
                    f"CPU reservation is already locked: {self.path}") from error
            raw = os.read(self.descriptor, 1024 * 1024)
            payload = parse_reservation(raw, self.cpu, self.sibling)
            self.identity = {
                "path": str(self.path),
                "sha256": sha256_bytes(raw),
                "payload": payload,
                "lock": "exclusive_nonblocking",
            }
            self.anchor.validate_current()
            return self.identity
        except BaseException:
            self._release()
            raise

    def _release(self) -> None:
        try:
            if self.descriptor is not None:
                descriptor = self.descriptor
                self.descriptor = None
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(descriptor)
        finally:
            self.anchor.__exit__(None, None, None)

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, exc, tb
        self._release()


def parse_reservation(raw: bytes, cpu: int, sibling: int) -> dict[str, Any]:
    require(type(cpu) is int and type(sibling) is int and
            0 <= cpu <= MAX_CPU_ID and 0 <= sibling <= MAX_CPU_ID and
            cpu != sibling,
            "CPU reservation request pair is invalid")
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"invalid CPU reservation JSON: {error}") from error
    expected_keys = {
        "benchmark_cpu", "nonce", "owner", "reserved_sibling", "schema", "status",
    }
    require(isinstance(payload, dict) and set(payload) == expected_keys,
            "CPU reservation has unexpected or missing fields")
    require(raw == canonical_bytes(payload),
            "CPU reservation is not canonical JSON without a trailing newline")
    require(payload["schema"] == RESERVATION_SCHEMA and payload["status"] == "held",
            "CPU reservation is not a held v1 reservation")
    require(type(payload["benchmark_cpu"]) is int and
            type(payload["reserved_sibling"]) is int and
            0 <= payload["benchmark_cpu"] <= MAX_CPU_ID and
            0 <= payload["reserved_sibling"] <= MAX_CPU_ID and
            payload["benchmark_cpu"] != payload["reserved_sibling"] and
            payload["benchmark_cpu"] == cpu and
            payload["reserved_sibling"] == sibling,
            "CPU reservation pair does not match the run request")
    require(isinstance(payload["owner"], str) and payload["owner"].strip(),
            "CPU reservation owner is empty")
    require(isinstance(payload["nonce"], str) and 8 <= len(payload["nonce"]) <= 256,
            "CPU reservation nonce is invalid")
    return payload


def validate_reservation_current(identity: Mapping[str, Any]) -> None:
    path = Path(str(identity.get("path", "")))
    require(path.is_file(), "CPU reservation disappeared during the campaign")
    raw = path.read_bytes()
    payload = identity.get("payload")
    require(isinstance(payload, dict), "retained CPU reservation payload is invalid")
    parsed = parse_reservation(
        raw, int(payload["benchmark_cpu"]), int(payload["reserved_sibling"]))
    require(parsed == payload and sha256_bytes(raw) == identity.get("sha256"),
            "CPU reservation changed during the campaign")


class XorShift64:
    def __init__(self, seed: int):
        self.state = seed if seed else 0x9E3779B97F4A7C15

    def next(self) -> int:
        value = self.state
        value ^= (value << 13) & MASK64
        value ^= value >> 7
        value ^= (value << 17) & MASK64
        self.state = value & MASK64
        return self.state


def expected_losses(cell: Cell) -> list[int]:
    order = list(range(cell.k))
    random = XorShift64(cell.seed ^ 0xD1B54A32D192ED03)
    for remaining in range(len(order), 1, -1):
        selected = random.next() % remaining
        order[remaining - 1], order[selected] = order[selected], order[remaining - 1]
    return sorted(order[:cell.losses])


def ceil_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def validate_cell(cell: Cell) -> None:
    require(type(cell.identifier) is str and
            all(type(value) is int for value in (
                cell.k, cell.r, cell.shard_bytes, cell.losses, cell.seed)),
            "cell fields have invalid types")
    require(re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", cell.identifier) is not None,
            f"invalid cell identifier {cell.identifier!r}")
    require(cell.k > 0 and cell.r > 0 and cell.r <= cell.k,
            f"cell {cell.identifier} is outside exact-main R <= K")
    require(0 < cell.shard_bytes <= MAX_SHARD_BYTES and
            cell.shard_bytes % 64 == 0,
            f"cell {cell.identifier} shard bytes are not a positive multiple of 64")
    require(0 <= cell.losses <= cell.r,
            f"cell {cell.identifier} has invalid loss count")
    parent = ceil_power_of_two(cell.k + ceil_power_of_two(cell.r))
    require(parent <= 65536, f"cell {cell.identifier} exceeds the legacy field")
    require(0 <= cell.seed <= MASK64, f"cell {cell.identifier} seed exceeds uint64")


def candidate_mode_for_campaign(campaign: Mapping[str, Any]) -> str:
    """Return the explicit v4+ mode, or AUTO for replay-only older bundles."""
    mode = campaign.get("candidate_mode", "auto")
    require(isinstance(mode, str) and mode in CANDIDATE_MODES,
            "campaign candidate mode is invalid")
    return mode


def validate_candidate_mode_schema(
    campaign: Mapping[str, Any], raw_schema: str
) -> str:
    if raw_schema in CANDIDATE_MODE_SCHEMAS:
        require("candidate_mode" in campaign,
                "current campaign does not bind its candidate mode")
    else:
        require("candidate_mode" not in campaign,
                "historical campaign contains an unversioned candidate mode")
    return candidate_mode_for_campaign(campaign)


def candidate_mode_flags(mode: str) -> dict[str, bool]:
    require(mode in CANDIDATE_MODES, "candidate mode is invalid")
    return {
        "force_generic_decode": mode == "generic",
        "force_specialized_decode": mode in {"materialized", "tiled"},
        "force_tiled_decode": mode == "tiled",
        "force_materialized_decode": mode == "materialized",
    }


def candidate_mode_arguments(mode: str) -> tuple[str, ...]:
    require(mode in CANDIDATE_MODES, "candidate mode is invalid")
    if mode == "auto":
        return ()
    if mode == "generic":
        return ("--force-generic",)
    if mode == "materialized":
        return ("--force-specialized", "--force-materialized")
    return ("--force-specialized", "--force-tiled")


def finite_number(value: object, what: str, positive: bool = True) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{what} is not numeric")
    result = float(value)
    require(math.isfinite(result), f"{what} is not finite")
    if positive:
        require(result > 0, f"{what} is not positive")
    else:
        require(result >= 0, f"{what} is negative")
    return result


def close_enough(actual: float, expected: float) -> bool:
    return abs(actual - expected) <= max(0.000002, abs(expected) * 0.000002)


def median(values: Sequence[float]) -> float:
    return statistics.median(values)


def validate_summary(
    summary: object, iterations: int, setup: bool = False
) -> list[float]:
    require(isinstance(summary, dict), "timing summary is not an object")
    prefix = "" if setup else "_per_batch_call"
    names = {
        "median": f"median_us{prefix}",
        "mad": f"mad_us{prefix}",
        "minimum": f"minimum_us{prefix}",
        "maximum": f"maximum_us{prefix}",
        "samples": "samples_us" if setup else "samples_us_per_batch_call",
    }
    samples_value = summary.get(names["samples"])
    require(isinstance(samples_value, list) and len(samples_value) == iterations,
            f"{names['samples']} does not contain exactly {iterations} samples")
    samples = [finite_number(value, names["samples"]) for value in samples_value]
    derived_median = median(samples)
    deviations = [abs(value - derived_median) for value in samples]
    expected = {
        names["median"]: derived_median,
        names["mad"]: median(deviations),
        names["minimum"]: min(samples),
        names["maximum"]: max(samples),
    }
    for name, derived in expected.items():
        actual = finite_number(summary.get(name), name, positive=name != names["mad"])
        require(close_enough(actual, derived), f"{name} is not derived from raw samples")
    return samples


def validate_digest_object(value: object) -> dict[str, str]:
    require(isinstance(value, dict), "workload_digests is not an object")
    require(value.get("algorithm") == "fnv1a64", "wrong workload digest algorithm")
    result: dict[str, str] = {}
    for name in ("original_data", "transmitted_parity", "recovered_originals"):
        digest = value.get(name)
        require(isinstance(digest, str) and HEX64.fullmatch(digest) is not None,
                f"invalid FNV-1a digest {name}")
        result[name] = digest
    return result


def expected_parameters(cell: Cell, campaign: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "K": cell.k,
        "R": cell.r,
        "shard_bytes": cell.shard_bytes,
        "loss_count": cell.losses,
        "missing_original_indices": expected_losses(cell),
        "batch": 1,
        "reuse": campaign["reuse"],
        "iterations": campaign["iterations"],
        "warmup": campaign["warmup"],
        "thread_count": 1,
        "seed": cell.seed,
    }


def validate_result(
    implementation: str,
    value: object,
    cell: Cell,
    campaign: Mapping[str, Any],
    raw_schema: str = RAW_SCHEMA,
) -> dict[str, Any]:
    require(raw_schema in RAW_TO_CMAKE_IDENTITY,
            "unsupported result-container schema")
    require(isinstance(value, dict), "benchmark output is not a JSON object")
    expected_schema = (
        "leopard-main-benchmark-v1" if implementation == "baseline"
        else "leopard2-benchmark-v2")
    require(value.get("schema") == expected_schema,
            f"{implementation} returned wrong schema")
    parameters = value.get("parameters")
    require(isinstance(parameters, dict), "parameters is not an object")
    for name, expected in expected_parameters(cell, campaign).items():
        require(parameters.get(name) == expected,
                f"{implementation} parameter {name} differs: "
                f"{parameters.get(name)!r} != {expected!r}")
    resolved = value.get("resolved")
    require(isinstance(resolved, dict), "resolved identity is not an object")
    padded = ceil_power_of_two(cell.r)
    parent = ceil_power_of_two(cell.k + padded)
    field = "gf8" if parent <= 256 else "gf16"
    for name, expected in {
        "profile": "legacy_high_v1", "field": field,
        "parent_count": parent, "padded_side": padded,
    }.items():
        require(resolved.get(name) == expected,
                f"{implementation} resolved {name} differs")
    if implementation == "baseline":
        require(resolved.get("thread_count") == 1,
                "baseline resolved more than one thread")
        require(value.get("build", {}).get("main_source_commit") == MAIN_COMMIT,
                "baseline did not attest the exact main commit")
        require(value.get("correctness", {}).get("round_trip") is True,
                "baseline round trip failed")
    else:
        candidate_mode = validate_candidate_mode_schema(campaign, raw_schema)
        mode_flags = candidate_mode_flags(candidate_mode)
        required_candidate = {
            "requested_profile": "legacy_high_v1",
            "requested_field": "auto",
            "requested_backend": "auto",
            "skip_legacy": True,
            "retain_samples": True,
            "force_generic_decode": mode_flags["force_generic_decode"],
            "force_specialized_decode": mode_flags["force_specialized_decode"],
        }
        for name, expected in required_candidate.items():
            require(parameters.get(name) == expected,
                    f"candidate option {name} is not comparison-safe")
        selector_names = (
            "force_tiled_decode", "force_materialized_decode")
        if raw_schema in WORKSPACE_SELECTOR_SCHEMAS:
            for name in selector_names:
                require(type(parameters.get(name)) is bool and
                        parameters[name] is mode_flags[name],
                        f"candidate option {name} is not current-schema safe")
        else:
            require(all(name not in parameters for name in selector_names),
                    "historical candidate result contains unversioned workspace selectors")
        require(resolved.get("thread_count") == 1,
                "candidate resolved more than one thread")
        require(resolved.get("backend") in {
                    "scalar", "ssse3", "avx2", "avx512", "neon"},
                "candidate did not resolve a production backend")
        require(value.get("correctness", {}).get("leopard2_round_trip") is True,
                "candidate round trip failed")
        legacy = value.get("legacy")
        require(isinstance(legacy, dict) and legacy.get("available") is False and
                legacy.get("unavailable_reason") == "disabled by --skip-legacy",
                "candidate silently ran the in-tree legacy comparison")
    digests = validate_digest_object(value.get("workload_digests"))
    metrics = value.get("metrics")
    require(isinstance(metrics, dict), "metrics is not an object")
    iterations = int(campaign["iterations"])
    encode = validate_summary(metrics.get("encode_execution"), iterations)
    if implementation == "baseline":
        require(metrics.get("codec_setup") is None and
                metrics.get("decode_timing_includes_setup") is True,
                "baseline decode setup semantics are ambiguous")
        decode = validate_summary(metrics.get("decode_including_setup"), iterations)
        return {
            "digests": digests,
            "backend": "exact_main_native",
            "encode": encode,
            "decode": decode,
        }
    codec_setup = validate_summary(metrics.get("codec_setup"), iterations, setup=True)
    plan_setup = validate_summary(metrics.get("decode_plan_setup"), iterations, setup=True)
    decode = validate_summary(metrics.get("decode_execution"), iterations)
    amortized = metrics.get("decode_amortized_at_reuse")
    require(isinstance(amortized, dict) and
            amortized.get("reuse_count") == campaign["reuse"],
            "candidate amortized decode has wrong reuse")
    expected_amortized = median(decode) + median(plan_setup) / campaign["reuse"]
    actual_amortized = finite_number(
        amortized.get("derived_median_us_per_batch_call"), "amortized decode median")
    require(close_enough(actual_amortized, expected_amortized),
            "candidate amortized decode is not derived from plan and execution")
    return {
        "digests": digests,
        "backend": resolved["backend"],
        "encode": encode,
        "codec_setup": codec_setup,
        "plan_setup": plan_setup,
        "decode": decode,
    }


def child_time(normalized: Mapping[str, Any], metric: str, reuse: int) -> float:
    if metric == "encode":
        return median(normalized["encode"])
    if "plan_setup" not in normalized:
        return median(normalized["decode"])
    divisor = 1 if metric == "decode_first_use" else reuse
    # Plan creation and execution are measured in separate loops.  Combine their
    # retained medians rather than pretending sample i from each loop was one
    # jointly timed first-use call.  Codec creation is intentionally excluded.
    return median(normalized["decode"]) + median(normalized["plan_setup"]) / divisor


def confidence_interval(log_ratios: Sequence[float]) -> dict[str, Any]:
    # The two contrasts within one ABBA round share drift and are not independent.
    # Each independent observation is therefore the mean log contrast for a whole
    # ABBA round; three rounds use Student-t with df=2.
    require(len(log_ratios) == ROUNDS,
            f"paired analysis requires {ROUNDS} independent round contrasts")
    mean = statistics.fmean(log_ratios)
    deviation = statistics.stdev(log_ratios)
    # Student-t 97.5th percentile with df=2.
    critical = 4.302652729911275
    half_width = critical * deviation / math.sqrt(len(log_ratios))
    lower = math.exp(mean - half_width)
    upper = math.exp(mean + half_width)
    speedup = math.exp(mean)
    return {
        "ratio_orientation": "baseline_time_over_candidate_time",
        "independent_round_count": len(log_ratios),
        "degrees_of_freedom": len(log_ratios) - 1,
        "independent_round_log_contrasts": list(log_ratios),
        "geometric_speedup": speedup,
        "ci95_lower": lower,
        "ci95_upper": upper,
        "faster_ci_excludes_one": lower > 1.0,
        "promotion_lower_bound_at_least_1_05": lower >= 1.05,
        "performance_result_does_not_affect_evidence_validity": True,
    }


def analyze_cell(records: Sequence[Mapping[str, Any]], reuse: int) -> dict[str, Any]:
    require(len(records) == ROUNDS * len(ORDER), "cell has wrong invocation count")
    metrics = ("encode", "decode_first_use", "decode_reuse_amortized")
    observations: dict[str, list[float]] = {name: [] for name in metrics}
    raw_pairs: dict[str, list[list[float]]] = {name: [] for name in metrics}
    for round_index in range(ROUNDS):
        group = records[round_index * 4:(round_index + 1) * 4]
        require(tuple(item["implementation"] for item in group) == ORDER,
                f"round {round_index} is not ABBA")
        round_pairs: dict[str, list[float]] = {name: [] for name in metrics}
        for baseline_index, candidate_index in ((0, 1), (3, 2)):
            baseline = group[baseline_index]["normalized"]
            candidate = group[candidate_index]["normalized"]
            for metric in metrics:
                baseline_time = child_time(baseline, metric, reuse)
                candidate_time = child_time(candidate, metric, reuse)
                round_pairs[metric].append(math.log(baseline_time / candidate_time))
        for metric in metrics:
            require(len(round_pairs[metric]) == 2,
                    f"round {round_index} does not contain two ABBA contrasts")
            raw_pairs[metric].append(round_pairs[metric])
            observations[metric].append(statistics.fmean(round_pairs[metric]))
    result: dict[str, Any] = {}
    for name, values in observations.items():
        result[name] = confidence_interval(values)
        result[name]["within_round_log_contrasts"] = raw_pairs[name]
        result[name]["constituent_pair_count"] = 2 * ROUNDS
    return result


def analyze(invocations: Sequence[Mapping[str, Any]], campaign: Mapping[str, Any]) -> dict[str, Any]:
    by_cell: dict[str, list[Mapping[str, Any]]] = {
        cell["identifier"]: [] for cell in campaign["cells"]
    }
    for invocation in invocations:
        identifier = invocation.get("cell_id")
        require(identifier in by_cell, f"unknown invocation cell {identifier!r}")
        by_cell[identifier].append(invocation)
    return {
        identifier: analyze_cell(records, int(campaign["reuse"]))
        for identifier, records in by_cell.items()
    }


def safe_evidence_path(root: Path, relative: object) -> Path:
    require(isinstance(relative, str) and relative and not os.path.isabs(relative),
            "evidence path is not relative")
    path = (root / relative).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as error:
        raise EvidenceError(f"evidence path escapes output directory: {relative}") from error
    return path


def validate_complete_cpu_policy_record(
    record: object, label: str, expected_cpu: int, expected_sibling: int,
) -> dict[str, Any]:
    require(isinstance(record, dict) and set(record) == {
                "cpu", "online", "cpuinfo", "topology", "frequency_policy",
                "cache_hierarchy", "cache_index_inventory",
                "cache_directory_inventory_text", "numa_nodes",
                "numa_node_inventory", "core_class"},
            f"host {label} policy identity shape differs")
    require(type(record.get("cpu")) is int and record["cpu"] == expected_cpu and
            (record.get("online") is None or isinstance(record["online"], str)),
            f"host {label} CPU identity is invalid")
    cpuinfo = record.get("cpuinfo")
    allowed_cpuinfo = {
        "processor", "vendor_id", "cpu family", "model", "model name",
        "stepping", "microcode", "flags", "Features", "CPU implementer",
        "CPU architecture", "CPU variant", "CPU part", "CPU revision",
    }
    require(isinstance(cpuinfo, dict) and cpuinfo and
            set(cpuinfo).issubset(allowed_cpuinfo) and
            cpuinfo.get("processor") == str(expected_cpu) and
            any(key in cpuinfo for key in ("model name", "CPU part")) and
            all(isinstance(item, str) for item in cpuinfo.values()),
            f"host {label} processor identity is invalid")
    topology = record.get("topology")
    topology_names = {
        "core_id", "physical_package_id", "die_id", "cluster_id",
        "thread_siblings_list", "core_siblings_list",
    }
    require(isinstance(topology, dict) and set(topology) == topology_names and
            all(item is None or isinstance(item, str)
                for item in topology.values()) and
            isinstance(topology["thread_siblings_list"], str) and
            parse_cpu_list(topology["thread_siblings_list"]) == {
                expected_cpu, expected_sibling} and
            isinstance(topology["core_siblings_list"], str) and
            {expected_cpu, expected_sibling}.issubset(
                parse_cpu_list(topology["core_siblings_list"])),
            f"host {label} topology is not the exact reserved SMT pair")
    policy = record.get("frequency_policy")
    policy_names = {
        "scaling_driver", "scaling_governor", "energy_performance_preference",
        "scaling_min_freq", "scaling_max_freq", "cpuinfo_min_freq",
        "cpuinfo_max_freq",
    }
    require(isinstance(policy, dict) and set(policy) == policy_names and
            all(item is None or isinstance(item, str) for item in policy.values()),
            f"host {label} frequency policy identity differs")
    caches = record.get("cache_hierarchy")
    cache_names = {
        "index", "level", "type", "size", "coherency_line_size",
        "number_of_sets", "ways_of_associativity", "physical_line_partition",
        "shared_cpu_list", "shared_cpu_map", "allocation_policy", "write_policy",
    }
    require(isinstance(caches, list) and caches,
            f"host {label} cache hierarchy is absent")
    indices: list[int] = []
    for cache in caches:
        require(isinstance(cache, dict) and set(cache) == cache_names and
                type(cache.get("index")) is int and cache["index"] >= 0 and
                all(cache.get(name) is None or isinstance(cache[name], str)
                    for name in cache_names - {"index"}) and
                all(isinstance(cache.get(name), str) and cache[name]
                    for name in ("level", "type", "size",
                                 "coherency_line_size", "shared_cpu_list",
                                 "shared_cpu_map")) and
                expected_cpu in parse_cpu_list(cache["shared_cpu_list"]),
                f"host {label} cache identity is invalid")
        indices.append(cache["index"])
    require(indices == sorted(indices) and len(indices) == len(set(indices)),
            f"host {label} cache indices are not sorted and unique")
    cache_inventory = record.get("cache_index_inventory")
    raw_cache_inventory = validate_sysfs_numeric_directory_inventory(
        record.get("cache_directory_inventory_text"), "index",
        f"host {label} raw cache-directory inventory")
    require(isinstance(cache_inventory, list) and
            cache_inventory == [f"index{index}" for index in indices] and
            cache_inventory == raw_cache_inventory,
            f"host {label} cache index inventory differs from retained caches")
    nodes = record.get("numa_nodes")
    node_inventory = record.get("numa_node_inventory")
    require(isinstance(nodes, list) and nodes and nodes == sorted(set(nodes)) and
            all(type(node) is int and node >= 0 for node in nodes),
            f"host {label} NUMA identity is invalid")
    require(isinstance(node_inventory, list) and
            node_inventory == [f"node{node}" for node in nodes],
            f"host {label} NUMA-node inventory differs")
    core_class = record.get("core_class")
    require(isinstance(core_class, dict) and set(core_class) == {
                "core_type", "cpu_capacity"} and
            all(item is None or isinstance(item, str)
                for item in core_class.values()),
            f"host {label} core-class identity differs")
    return record


def validate_host_record(
    value: object, cpu: int, sibling: int, allowed: Sequence[int], raw_schema: str,
) -> None:
    require(isinstance(value, dict), "host identity is not an object")
    if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
        require(set(value) == {
                    "system", "allowed_cpu_set_at_launch", "online_cpu_set",
                    "online_cpu_list_text", "online_node_list_text",
                    "benchmark_cpu", "reserved_sibling", "turbo_and_pstate"},
                "complete host identity shape differs")
        system = value.get("system")
        require(isinstance(system, dict) and set(system) == {
                    "system", "node", "release", "version", "machine",
                    "python", "page_size"} and
                all(isinstance(system.get(name), str) and system[name]
                    for name in ("system", "node", "release", "version",
                                 "machine", "python")) and
                type(system.get("page_size")) is int and system["page_size"] > 0,
                "complete host system identity differs")
        launch = value.get("allowed_cpu_set_at_launch")
        online = value.get("online_cpu_set")
        require(isinstance(launch, list) and launch == list(allowed) and
                launch == sorted(set(launch)) and
                3 <= len(launch) <= MAX_CPU_LIST_ENTRIES and
                all(type(item) is int and 0 <= item <= MAX_CPU_ID
                    for item in launch) and
                isinstance(online, list) and online == sorted(set(online)) and
                1 <= len(online) <= MAX_CPU_LIST_ENTRIES and
                all(type(item) is int and 0 <= item <= MAX_CPU_ID
                    for item in online),
                "complete host CPU-set identity differs")
        online_cpu_text = validate_complete_text_identity(
            value.get("online_cpu_list_text"), "host online CPU list")
        online_node_text = validate_complete_text_identity(
            value.get("online_node_list_text"), "host online NUMA-node list")
        online_nodes = parse_cpu_list(online_node_text["text"])
        require(parse_cpu_list(online_cpu_text["text"]) == set(online) and
                online_nodes,
                "complete host online CPU/node summaries differ from sysfs text")
        benchmark_record = validate_complete_cpu_policy_record(
            value.get("benchmark_cpu"), "benchmark_cpu", cpu, sibling)
        sibling_record = validate_complete_cpu_policy_record(
            value.get("reserved_sibling"), "reserved_sibling", sibling, cpu)
        require(benchmark_record["cache_hierarchy"] ==
                    sibling_record["cache_hierarchy"] and
                benchmark_record["numa_nodes"] == sibling_record["numa_nodes"] and
                set(benchmark_record["numa_nodes"]).issubset(online_nodes) and
                benchmark_record["core_class"] == sibling_record["core_class"],
                "reserved SMT pair differs in cache, NUMA, or core-class identity")
        turbo = value.get("turbo_and_pstate")
        turbo_paths = {
            "/sys/devices/system/cpu/intel_pstate/no_turbo",
            "/sys/devices/system/cpu/amd_pstate/status",
            "/sys/devices/system/cpu/cpufreq/boost",
        }
        require(isinstance(turbo, dict) and set(turbo) == turbo_paths and
                all(item is None or isinstance(item, str)
                    for item in turbo.values()),
                "complete host turbo/pstate identity differs")
    require(value.get("allowed_cpu_set_at_launch") == list(allowed),
            "host identity has the wrong launch affinity")
    online = value.get("online_cpu_set")
    require(type(online) is list and online == sorted(set(online)) and
            all(type(item) is int and 0 <= item <= MAX_CPU_ID for item in online) and
            cpu in online and sibling in online,
            "reserved CPUs were not retained as online")
    for name, expected_cpu, expected_sibling in (
        ("benchmark_cpu", cpu, sibling), ("reserved_sibling", sibling, cpu)):
        record = value.get(name)
        require(isinstance(record, dict) and
                type(record.get("cpu")) is int and
                record.get("cpu") == expected_cpu,
                f"host {name} identity is invalid")
        cpuinfo = record.get("cpuinfo")
        topology = record.get("topology")
        policy = record.get("frequency_policy")
        require(isinstance(cpuinfo, dict) and
                any(key in cpuinfo for key in ("model name", "CPU part")),
                f"host {name} model identity is missing")
        require(isinstance(topology, dict) and
                isinstance(topology.get("thread_siblings_list"), str) and
                parse_cpu_list(topology["thread_siblings_list"]) ==
                {expected_cpu, expected_sibling},
                f"host {name} topology is not exactly the reserved SMT pair")
        require(isinstance(policy, dict) and
                all(key in policy for key in (
                    "scaling_driver", "scaling_governor",
                    "energy_performance_preference")),
                f"host {name} frequency policy is missing")
    require(isinstance(value.get("turbo_and_pstate"), dict),
            "host turbo/pstate identity is missing")
    require(isinstance(value.get("system"), dict) and
            isinstance(value["system"].get("release"), str),
            "host system identity is missing")


def validate_raw(
    raw: object,
    output: Path | None,
    check_files: bool,
    check_current_inputs: bool,
) -> dict[str, Any]:
    raw = verify_signature(raw, "raw bundle")
    raw_schema = raw.get("schema")
    require(isinstance(raw_schema, str) and
            raw_schema in RAW_TO_CMAKE_IDENTITY,
            "wrong raw bundle schema")
    campaign = raw.get("campaign")
    require(isinstance(campaign, dict), "campaign is not an object")
    validate_candidate_mode_schema(campaign, raw_schema)
    require(type(campaign.get("rounds")) is int and
            campaign["rounds"] == ROUNDS and
            tuple(campaign.get("order", ())) == ORDER,
            "campaign does not contain exactly three ABBA rounds")
    require(type(campaign.get("batch")) is int and campaign["batch"] == 1 and
            type(campaign.get("threads")) is int and campaign["threads"] == 1,
            "campaign is not a one-stripe, one-thread comparison")
    require(campaign.get("child_environment") == CHILD_ENVIRONMENT,
            "campaign child environment is not the strict comparison environment")
    cpu = campaign.get("benchmark_cpu")
    sibling = campaign.get("reserved_sibling")
    require(type(cpu) is int and type(sibling) is int and
            0 <= cpu <= MAX_CPU_ID and 0 <= sibling <= MAX_CPU_ID and
            cpu != sibling,
            "campaign has no valid reserved CPU pair")
    require(type(campaign.get("iterations")) is int and
            3 <= campaign["iterations"] <= MAX_CAMPAIGN_COUNT,
            "campaign has too few timing samples")
    require(type(campaign.get("reuse")) is int and
            1 <= campaign["reuse"] <= MAX_CAMPAIGN_COUNT,
            "campaign reuse is invalid")
    require(type(campaign.get("warmup")) is int and
            1 <= campaign["warmup"] <= MAX_CAMPAIGN_COUNT,
            "campaign warmup is invalid")
    timeout = campaign.get("timeout_seconds")
    require(isinstance(timeout, (int, float)) and not isinstance(timeout, bool) and
            math.isfinite(timeout) and 0 < timeout <= MAX_COMMAND_TIMEOUT_SECONDS,
            "campaign timeout is invalid")
    require(campaign.get("statistics") == statistics_policy(),
            "campaign statistics policy is not the authoritative clustered ABBA policy")
    allowed = campaign.get("allowed_cpu_set_at_launch")
    require(type(allowed) is list and allowed == sorted(set(allowed)) and
            3 <= len(allowed) <= MAX_CPU_LIST_ENTRIES and
            all(type(item) is int and 0 <= item <= MAX_CPU_ID for item in allowed) and
            cpu in allowed and sibling in allowed and len(set(allowed) - {cpu, sibling}) > 0,
            "campaign launch affinity is invalid")
    host_initial = raw.get("host_initial")
    host_final = raw.get("host_final")
    require(host_initial == host_final, "host policy/topology changed during campaign")
    validate_host_record(host_initial, cpu, sibling, allowed, raw_schema)
    if raw_schema in ISOLATION_SCHEMAS:
        validate_isolation(raw.get("isolation"), cpu, sibling)
    else:
        require("isolation" not in raw,
                "legacy raw schema contains unversioned isolation evidence")
    cells_value = campaign.get("cells")
    require(type(cells_value) is list and
            1 <= len(cells_value) <= MAX_CAMPAIGN_CELLS and
            all(isinstance(value, dict) and set(value) == {
                "identifier", "k", "r", "shard_bytes", "losses", "seed"}
                for value in cells_value),
            "campaign cell list is incomplete or invalid")
    cells = [Cell(**value) for value in cells_value]
    require(len({cell.identifier for cell in cells}) == len(cells),
            "campaign cell identifiers are not unique")
    for cell in cells:
        validate_cell(cell)
    input_spec = raw.get("input_specification")
    initial = raw.get("identities_initial")
    final = raw.get("identities_final")
    require(isinstance(input_spec, dict) and isinstance(initial, dict) and
            isinstance(final, dict), "input identity is missing")
    expected_specification_keys = {
        "runner", "taskset", "ldd", "baseline_executable", "candidate_executable",
        "baseline_archive", "candidate_archive", "baseline_build_dir",
        "candidate_build_dir", "baseline_source_root", "candidate_source_root",
        "candidate_commit",
    }
    require(set(input_spec) == expected_specification_keys and
            all(isinstance(value, str) and value for value in input_spec.values()),
            "input specification is incomplete or has unexpected fields")
    require(re.fullmatch(r"[0-9a-f]{40}", input_spec["candidate_commit"]) is not None,
            "candidate commit is not a full lowercase SHA-1")
    require(initial == final, "input identities changed during the campaign")
    if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
        validate_complete_input_snapshot(input_spec, initial, raw_schema)
    validate_candidate_cmake_identity(input_spec, initial, raw_schema)
    reservation = raw.get("reservation")
    require(isinstance(reservation, dict) and
            reservation.get("lock") == "exclusive_nonblocking" and
            isinstance(reservation.get("path"), str) and
            isinstance(reservation.get("sha256"), str) and
            HEX256.fullmatch(reservation["sha256"]) is not None,
            "retained CPU reservation identity is invalid")
    reservation_payload = reservation.get("payload")
    require(isinstance(reservation_payload, dict),
            "retained CPU reservation payload is missing")
    require(parse_reservation(canonical_bytes(reservation_payload), cpu, sibling) ==
            reservation_payload,
            "retained CPU reservation semantics differ from the campaign")
    require(reservation["sha256"] ==
            sha256_bytes(canonical_bytes(reservation_payload)),
            "retained CPU reservation hash does not match its canonical payload")
    if raw_schema in SUPERVISION_SCHEMAS:
        require("supervision" in raw,
                "complete raw bundle omits its supervision handshake field")
        supervision = raw.get("supervision")
        if supervision is not None:
            validate_supervision(supervision, campaign, reservation, raw["isolation"])
    else:
        require("supervision" not in raw,
                "historical raw bundle contains unversioned supervision data")
    if check_current_inputs:
        require(input_snapshots_equal(
                    input_snapshot(input_spec, raw_schema), initial, raw_schema),
                "current executable/archive/source identity differs from retained evidence")
        validate_reservation_current(reservation)
    invocations = raw.get("invocations")
    expected_count = len(cells) * ROUNDS * len(ORDER)
    require(isinstance(invocations, list) and len(invocations) == expected_count,
            f"campaign has {0 if not isinstance(invocations, list) else len(invocations)} "
            f"invocations, expected {expected_count}")
    digest_by_cell: dict[str, dict[str, str]] = {}
    candidate_backend_by_cell: dict[str, str] = {}
    cell_by_id = {cell.identifier: cell for cell in cells}
    expected_sequence = [
        (cell.identifier, round_index, slot, implementation)
        for cell in cells for round_index in range(ROUNDS)
        for slot, implementation in enumerate(ORDER)
    ]
    total_child_duration_ns = 0
    for invocation, expected in zip(invocations, expected_sequence):
        require(isinstance(invocation, dict), "invocation is not an object")
        observed = (
            invocation.get("cell_id"), invocation.get("round"),
            invocation.get("slot"), invocation.get("implementation"))
        require(observed == expected,
                f"invocation order/relabel mismatch: {observed!r} != {expected!r}")
        require(invocation.get("identity_before") == initial and
                invocation.get("identity_after") == initial,
                "an invocation observed a changed input identity")
        require(invocation.get("reservation_before") == reservation and
                invocation.get("reservation_after") == reservation,
                "an invocation observed a changed CPU reservation")
        require(invocation.get("returncode") == 0,
                "benchmark child did not exit successfully")
        duration_ns = invocation.get("duration_ns")
        require(isinstance(duration_ns, int) and not isinstance(duration_ns, bool) and
                duration_ns > 0,
                "benchmark child duration is not a positive integer")
        total_child_duration_ns += duration_ns
        cell = cell_by_id[expected[0]]
        expected_command = [
            input_spec["taskset"], "-c", str(cpu),
            *benchmark_arguments(expected[3], Path(input_spec[
                f"{expected[3]}_executable"]), cell, campaign),
        ]
        require(invocation.get("command") == expected_command and
                invocation.get("pinned_cpu") == cpu,
                "benchmark command or CPU pinning was edited")
        require(invocation.get("environment") == CHILD_ENVIRONMENT,
                "benchmark invocation inherited or retained an unsafe environment")
        result = invocation.get("result")
        if check_files:
            require(output is not None, "output root is required for file verification")
            for stream_name in ("stdout", "stderr"):
                stream = invocation.get(stream_name)
                require(isinstance(stream, dict), f"missing {stream_name} evidence")
                path = safe_evidence_path(output, stream.get("path"))
                require(path.is_file(), f"missing retained {stream_name}: {path}")
                data = path.read_bytes()
                require(stream.get("size") == len(data) and
                        stream.get("sha256") == sha256_bytes(data),
                        f"retained {stream_name} identity mismatch")
                if stream_name == "stdout":
                    try:
                        parsed = json.loads(data.decode("utf-8"))
                    except (UnicodeDecodeError, json.JSONDecodeError) as error:
                        raise EvidenceError(f"retained stdout is not JSON: {error}") from error
                    require(parsed == result, "parsed retained stdout differs from raw result")
        normalized = validate_result(
            expected[3], result, cell, campaign, raw_schema)
        require(invocation.get("normalized") == normalized,
                "retained normalized benchmark data was edited")
        digests = normalized["digests"]
        if expected[0] in digest_by_cell:
            require(digests == digest_by_cell[expected[0]],
                    f"workload FNV digests differ in cell {expected[0]}")
        else:
            digest_by_cell[expected[0]] = digests
        if expected[3] == "candidate":
            backend = normalized["backend"]
            if expected[0] in candidate_backend_by_cell:
                require(backend == candidate_backend_by_cell[expected[0]],
                        f"candidate backend changed within cell {expected[0]}")
            else:
                candidate_backend_by_cell[expected[0]] = backend
    if raw_schema in ISOLATION_SCHEMAS:
        isolation = raw["isolation"]
        elapsed_ns = isolation["after"]["monotonic_ns"] - \
            isolation["before"]["monotonic_ns"]
        require(elapsed_ns >= total_child_duration_ns,
                "isolation interval does not cover all benchmark child durations")
    calculated = analyze(invocations, campaign)
    require(raw.get("analysis") == calculated, "paired-log analysis was edited")
    return calculated


def benchmark_arguments(
    implementation: str, executable: Path, cell: Cell, campaign: Mapping[str, Any]
) -> list[str]:
    arguments = [
        str(executable), "--k", str(cell.k), "--r", str(cell.r),
        "--bytes", str(cell.shard_bytes), "--loss", str(cell.losses),
        "--batch", "1", "--reuse", str(campaign["reuse"]),
        "--iterations", str(campaign["iterations"]),
        "--warmup", str(campaign["warmup"]), "--threads", "1",
        "--seed", str(cell.seed),
    ]
    if implementation == "candidate":
        arguments.extend((
            "--profile", "high", "--field", "auto", "--backend", "auto",
            "--skip-legacy", "--retain-samples",
        ))
        arguments.extend(candidate_mode_arguments(
            candidate_mode_for_campaign(campaign)))
    arguments.extend(("--json", "-"))
    return arguments


def run_child(
    implementation: str,
    cell: Cell,
    round_index: int,
    slot: int,
    campaign: Mapping[str, Any],
    specification: Mapping[str, Any],
    initial_identity: Mapping[str, Any],
    reservation: Mapping[str, Any],
    output: Path,
    cpu: int,
    timeout: float,
) -> dict[str, Any]:
    executable = Path(specification[f"{implementation}_executable"])
    child_arguments = benchmark_arguments(implementation, executable, cell, campaign)
    command = [specification["taskset"], "-c", str(cpu), *child_arguments]
    before = input_snapshot(specification)
    require(before == initial_identity, "input identity changed before benchmark launch")
    validate_reservation_current(reservation)
    environment = dict(CHILD_ENVIRONMENT)
    start_utc = utc_now()
    start = time.monotonic_ns()
    completed = run_process_bounded(
        command, environment=environment, timeout=timeout,
        max_stdout=8 * 1024 * 1024, max_stderr=1024 * 1024)
    duration_ns = time.monotonic_ns() - start
    stem = f"invocations/{cell.identifier}/round-{round_index}/slot-{slot}-{implementation}"
    stdout_path = output / f"{stem}.stdout"
    stderr_path = output / f"{stem}.stderr"
    write_bytes_exclusive(stdout_path, completed.stdout)
    write_bytes_exclusive(stderr_path, completed.stderr)
    after = input_snapshot(specification)
    require(after == initial_identity, "input identity changed after benchmark launch")
    validate_reservation_current(reservation)
    require(completed.returncode == 0,
            f"{implementation} exited {completed.returncode}; see {stderr_path}")
    try:
        result = json.loads(completed.stdout.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"{implementation} stdout is not one JSON value: {error}") from error
    normalized = validate_result(
        implementation, result, cell, campaign, RAW_SCHEMA)
    return {
        "cell_id": cell.identifier,
        "round": round_index,
        "slot": slot,
        "implementation": implementation,
        "command": command,
        "environment": environment,
        "pinned_cpu": cpu,
        "started_utc": start_utc,
        "duration_ns": duration_ns,
        "returncode": completed.returncode,
        "stdout": {
            "path": str(stdout_path.relative_to(output)),
            "size": len(completed.stdout),
            "sha256": sha256_bytes(completed.stdout),
        },
        "stderr": {
            "path": str(stderr_path.relative_to(output)),
            "size": len(completed.stderr),
            "sha256": sha256_bytes(completed.stderr),
        },
        "result": result,
        "normalized": normalized,
        "identity_before": before,
        "identity_after": after,
        "reservation_before": reservation,
        "reservation_after": reservation,
    }


def parse_cell(text: str) -> Cell:
    components = text.split(":")
    require(len(components) == 6,
            "--cell must be ID:K:R:BYTES:LOSSES:SEED")
    try:
        cell = Cell(
            components[0], *(int(component, 10) for component in components[1:]))
    except ValueError as error:
        raise EvidenceError(f"invalid --cell {text!r}: {error}") from error
    validate_cell(cell)
    return cell


def cells_from_options(options: argparse.Namespace) -> tuple[Cell, ...]:
    if options.cell:
        cells = tuple(parse_cell(value) for value in options.cell)
    else:
        cells = REPRESENTATIVE_CELLS if options.preset == "representative" else SMOKE_CELLS
    require(len({cell.identifier for cell in cells}) == len(cells),
            "cell identifiers must be unique")
    return cells


def retained_file_records(output: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    failure_path = output / "failure.json"
    for path in sorted(candidate for candidate in output.rglob("*")
                       if candidate.is_file() and candidate != failure_path):
        records.append({
            "path": str(path.relative_to(output)),
            "size": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    return records


def validate_failure(
    value: object, output: Path, check_files: bool = True
) -> dict[str, Any]:
    failure = verify_signature(value, "failed campaign")
    failure_schema = failure.get("schema")
    expected_fields = {
        "campaign", "created_utc", "digest", "error", "error_type",
        "host_initial", "identities_initial", "input_specification",
        "invocations", "isolation", "pair_lease", "reservation",
        "retained_files", "schema", "status", "traceback", "valid"}
    if failure_schema in (FAILURE_SCHEMA_V5, FAILURE_SCHEMA):
        expected_fields.add("supervision")
    require(set(failure) == expected_fields,
        "failed campaign has unexpected or missing fields")
    require(isinstance(failure_schema, str) and
            failure_schema in FAILURE_TO_RAW_SCHEMA and
            failure.get("status") == "failed" and failure.get("valid") is False,
            "failed campaign status is invalid")
    require(all(isinstance(failure.get(name), str) and failure[name]
                for name in ("created_utc", "error", "error_type", "traceback")),
            "failed campaign diagnostic fields are invalid")
    campaign = failure.get("campaign")
    require(isinstance(campaign, dict), "failed campaign metadata is missing")
    raw_schema = FAILURE_TO_RAW_SCHEMA[failure_schema]
    validate_candidate_mode_schema(campaign, raw_schema)
    require(type(campaign.get("rounds")) is int and
            campaign["rounds"] == ROUNDS and
            tuple(campaign.get("order", ())) == ORDER and
            type(campaign.get("batch")) is int and campaign["batch"] == 1 and
            type(campaign.get("threads")) is int and campaign["threads"] == 1 and
            type(campaign.get("iterations")) is int and
            3 <= campaign["iterations"] <= MAX_CAMPAIGN_COUNT and
            type(campaign.get("reuse")) is int and
            1 <= campaign["reuse"] <= MAX_CAMPAIGN_COUNT and
            type(campaign.get("warmup")) is int and
            1 <= campaign["warmup"] <= MAX_CAMPAIGN_COUNT,
            "failed campaign scalar counts are invalid")
    timeout = campaign.get("timeout_seconds")
    require(isinstance(timeout, (int, float)) and not isinstance(timeout, bool) and
            math.isfinite(timeout) and 0 < timeout <= MAX_COMMAND_TIMEOUT_SECONDS,
            "failed campaign timeout is invalid")
    cpu = campaign.get("benchmark_cpu")
    sibling = campaign.get("reserved_sibling")
    require(type(cpu) is int and type(sibling) is int and
            0 <= cpu <= MAX_CPU_ID and 0 <= sibling <= MAX_CPU_ID and
            cpu != sibling,
            "failed campaign CPU pair is invalid")
    allowed_value = campaign.get("allowed_cpu_set_at_launch")
    if allowed_value is not None:
        require(type(allowed_value) is list and
                allowed_value == sorted(set(allowed_value)) and
                3 <= len(allowed_value) <= MAX_CPU_LIST_ENTRIES and
                all(type(item) is int and 0 <= item <= MAX_CPU_ID
                    for item in allowed_value) and
                cpu in allowed_value and sibling in allowed_value and
                len(set(allowed_value) - {cpu, sibling}) > 0,
                "failed campaign launch affinity is invalid")
    host_initial = failure.get("host_initial")
    if host_initial is not None:
        require(allowed_value is not None,
                "failed campaign host identity lacks launch affinity")
        validate_host_record(
            host_initial, cpu, sibling, allowed_value, raw_schema)
    pair_lease = failure.get("pair_lease")
    isolation = failure.get("isolation")
    if pair_lease is not None:
        validate_pair_lease_identity(pair_lease, cpu, sibling)
    if isolation is not None:
        validate_isolation(isolation, cpu, sibling, require_accepted=False)
        require(pair_lease == isolation["pair_lease"],
                "failed campaign isolation uses another pair lease")
    reservation = failure.get("reservation")
    if reservation is not None:
        require(isinstance(reservation, dict) and
                reservation.get("lock") == "exclusive_nonblocking" and
                isinstance(reservation.get("payload"), dict),
                "failed campaign reservation is invalid")
        payload = reservation["payload"]
        require(parse_reservation(canonical_bytes(payload), cpu, sibling) == payload and
                reservation.get("sha256") == sha256_bytes(canonical_bytes(payload)),
            "failed campaign reservation identity is invalid")
    supervision = failure.get("supervision")
    if supervision is not None:
        require(failure_schema in (FAILURE_SCHEMA_V5, FAILURE_SCHEMA) and
                reservation is not None and
                isolation is not None and isinstance(campaign, dict),
                "failed supervision handshake has no complete-schema context")
        validate_supervision(supervision, campaign, reservation, isolation)
    invocations = failure.get("invocations")
    require(isinstance(invocations, list), "failed invocation prefix is not a list")
    cells_value = campaign.get("cells")
    require(type(cells_value) is list and
            1 <= len(cells_value) <= MAX_CAMPAIGN_CELLS and
            all(isinstance(item, dict) and set(item) == {
                "identifier", "k", "r", "shard_bytes", "losses", "seed"}
                for item in cells_value),
            "failed campaign cells are missing")
    cells = [Cell(**item) for item in cells_value]
    for cell in cells:
        validate_cell(cell)
    expected_sequence = [
        (cell.identifier, round_index, slot, implementation)
        for cell in cells for round_index in range(ROUNDS)
        for slot, implementation in enumerate(ORDER)
    ]
    require(len(invocations) <= len(expected_sequence),
            "failed invocation prefix is too long")
    specification = failure.get("input_specification")
    initial = failure.get("identities_initial")
    if isinstance(initial, dict):
        require(isinstance(specification, dict),
                "failed campaign identity lacks input specification")
        if raw_schema in COMPLETE_EVIDENCE_SCHEMAS:
            validate_complete_input_snapshot(specification, initial, raw_schema)
        validate_candidate_cmake_identity(
            specification, initial, raw_schema)
    if invocations:
        require(isinstance(specification, dict) and isinstance(initial, dict) and
                reservation is not None,
                "failed invocation prefix lacks build or reservation identity")
    cell_by_id = {cell.identifier: cell for cell in cells}
    for invocation, expected in zip(invocations, expected_sequence):
        require(isinstance(invocation, dict) and
                (invocation.get("cell_id"), invocation.get("round"),
                 invocation.get("slot"), invocation.get("implementation")) == expected and
                invocation.get("returncode") == 0,
                "failed campaign invocation prefix was edited")
        duration_ns = invocation.get("duration_ns")
        require(isinstance(duration_ns, int) and not isinstance(duration_ns, bool) and
                duration_ns > 0,
                "failed campaign invocation duration is invalid")
        implementation = expected[3]
        cell = cell_by_id[expected[0]]
        expected_command = [
            specification["taskset"], "-c", str(cpu),
            *benchmark_arguments(implementation, Path(specification[
                f"{implementation}_executable"]), cell, campaign),
        ]
        require(invocation.get("command") == expected_command and
                invocation.get("environment") == CHILD_ENVIRONMENT and
                invocation.get("pinned_cpu") == cpu and
                invocation.get("identity_before") == initial and
                invocation.get("identity_after") == initial and
                invocation.get("reservation_before") == reservation and
                invocation.get("reservation_after") == reservation,
                "failed campaign invocation execution identity was edited")
        normalized = validate_result(
            implementation, invocation.get("result"), cell, campaign,
            raw_schema)
        require(invocation.get("normalized") == normalized,
                "failed campaign invocation result was edited")
    retained = failure.get("retained_files")
    require(isinstance(retained, list), "failed retained-file list is invalid")
    retained_paths: set[str] = set()
    retained_by_path: dict[str, dict[str, Any]] = {}
    for record in retained:
        require(isinstance(record, dict) and set(record) == {"path", "sha256", "size"} and
                isinstance(record.get("path"), str) and
                isinstance(record.get("size"), int) and
                not isinstance(record.get("size"), bool) and record["size"] >= 0 and
                isinstance(record.get("sha256"), str) and
                HEX256.fullmatch(record["sha256"]) is not None and
                record["path"] not in retained_paths,
                "failed retained-file identity is invalid")
        retained_paths.add(record["path"])
        retained_by_path[record["path"]] = record
        if check_files:
            path = safe_evidence_path(output, record["path"])
            require(path.is_file() and path.stat().st_size == record["size"] and
                    sha256_file(path) == record["sha256"],
                    "failed retained file is missing or changed")
    for invocation in invocations:
        for name in ("stdout", "stderr"):
            stream = invocation.get(name)
            require(isinstance(stream, dict) and
                    retained_by_path.get(stream.get("path")) == stream,
                    "failed invocation stream is not retained")
            if check_files and name == "stdout":
                path = safe_evidence_path(output, stream["path"])
                try:
                    parsed = json.loads(path.read_text(encoding="utf-8"))
                except (UnicodeDecodeError, json.JSONDecodeError) as error:
                    raise EvidenceError(
                        f"failed invocation stdout is not JSON: {error}") from error
                require(parsed == invocation.get("result"),
                        "failed invocation stdout differs from embedded result")
    if check_files:
        require(retained == retained_file_records(output),
                "failed output directory contains unbound or missing files")
    return failure


def run_campaign(options: argparse.Namespace) -> int:
    runner_started_monotonic_ns = time.monotonic_ns()
    execution_nonce = os.environ.get(SUPERVISION_NONCE_ENV)
    if execution_nonce is not None:
        require(HEX256.fullmatch(execution_nonce) is not None,
                "supervisor execution nonce is malformed")
    output = options.output.resolve()
    require(not output.exists(), f"output path already exists: {output}")
    output.mkdir(parents=True)
    cells = cells_from_options(options)
    for cell in cells:
        validate_cell(cell)
    taskset = Path(options.taskset).resolve(strict=True)
    ldd = Path(options.ldd).resolve(strict=True)
    require(taskset == Path("/usr/bin/taskset").resolve(strict=True),
            "authoritative comparison requires /usr/bin/taskset")
    require(ldd == Path("/usr/bin/ldd").resolve(strict=True),
            "authoritative comparison requires /usr/bin/ldd")
    specification = {
        "runner": str(Path(__file__).resolve()),
        "taskset": str(taskset),
        "ldd": str(ldd),
        "baseline_executable": str(options.baseline.resolve(strict=True)),
        "candidate_executable": str(options.candidate.resolve(strict=True)),
        "baseline_archive": str(options.baseline_archive.resolve(strict=True)),
        "candidate_archive": str(options.candidate_archive.resolve(strict=True)),
        "baseline_build_dir": str(options.baseline_build_dir.resolve(strict=True)),
        "candidate_build_dir": str(options.candidate_build_dir.resolve(strict=True)),
        "baseline_source_root": str(options.baseline_source_root.resolve(strict=True)),
        "candidate_source_root": str(options.candidate_source_root.resolve(strict=True)),
        "candidate_commit": options.candidate_commit,
    }
    campaign = {
        "rounds": ROUNDS,
        "order": list(ORDER),
        "cells": [asdict(cell) for cell in cells],
        "candidate_mode": options.candidate_mode,
        "batch": 1,
        "reuse": options.reuse,
        "iterations": options.iterations,
        "warmup": options.warmup,
        "threads": 1,
        "child_environment": dict(CHILD_ENVIRONMENT),
        "benchmark_cpu": options.cpu,
        "reserved_sibling": options.reserved_sibling,
        "timeout_seconds": options.timeout,
        "statistics": statistics_policy(),
    }
    require(options.iterations >= 3, "--iterations must be at least 3")
    require(options.reuse >= 1 and options.warmup >= 1,
            "--reuse and --warmup must be positive")
    require(math.isfinite(options.timeout) and options.timeout > 0,
            "--timeout must be positive and finite")
    isolation: dict[str, Any] | None = None
    host_initial: dict[str, Any] | None = None
    reservation: dict[str, Any] | None = None
    pair_lease: dict[str, Any] | None = None
    initial: dict[str, Any] | None = None
    invocations: list[dict[str, Any]] = []
    before_monotonic_ns: int | None = None
    before_cpu: dict[str, Any] | None = None
    before_sibling: dict[str, Any] | None = None
    supervision: dict[str, Any] | None = None
    try:
        allowed_at_launch, housekeeping = validate_topology(
            options.cpu, options.reserved_sibling)
        campaign["allowed_cpu_set_at_launch"] = sorted(allowed_at_launch)
        host_initial = host_identity(
            options.cpu, options.reserved_sibling, allowed_at_launch)
        pair_guard = PairLease(options.cpu, options.reserved_sibling)
        with Reservation(
            options.reservation_file, options.cpu, options.reserved_sibling
        ) as reservation, pair_guard as pair_lease:
            os.sched_setaffinity(0, housekeeping)
            initial = input_snapshot(specification)
            before_monotonic_ns = time.monotonic_ns()
            before_cpu = cpu_stat_snapshot(options.cpu)
            before_sibling = cpu_stat_snapshot(options.reserved_sibling)
            try:
                for cell in cells:
                    for round_index in range(ROUNDS):
                        for slot, implementation in enumerate(ORDER):
                            invocation = run_child(
                                implementation, cell, round_index, slot, campaign,
                                specification, initial, reservation, output,
                                options.cpu, options.timeout)
                            pair_guard.validate_current()
                            invocations.append(invocation)
            finally:
                if before_monotonic_ns is not None and before_cpu is not None and \
                        before_sibling is not None:
                    after_cpu = cpu_stat_snapshot(options.cpu)
                    after_sibling = cpu_stat_snapshot(options.reserved_sibling)
                    after_monotonic_ns = time.monotonic_ns()
                    isolation = isolation_record(
                        options.cpu, options.reserved_sibling, pair_lease,
                        before_monotonic_ns, after_monotonic_ns,
                        before_cpu, after_cpu, before_sibling, after_sibling)
                    pair_guard.validate_current()
            require(isolation is not None,
                    "campaign produced no scheduler isolation evidence")
            require(isolation["accepted"] is True,
                    "reserved SMT sibling performed non-idle work during the campaign")
            final = input_snapshot(specification)
            require(final == initial, "input identity changed during campaign")
            host_final = host_identity(
                options.cpu, options.reserved_sibling, allowed_at_launch)
            require(host_final == host_initial,
                    "host topology/frequency policy changed during campaign")
            analysis = analyze(invocations, campaign)
            if execution_nonce is not None:
                supervision = supervision_record(
                    execution_nonce, runner_started_monotonic_ns,
                    time.monotonic_ns(), campaign, reservation, isolation)
            raw = signed({
                "schema": RAW_SCHEMA,
                "created_utc": utc_now(),
                "validity_is_independent_of_speed": True,
                "campaign": campaign,
                "host_initial": host_initial,
                "isolation": isolation,
                "reservation": reservation,
                "supervision": supervision,
                "input_specification": specification,
                "identities_initial": initial,
                "invocations": invocations,
                "identities_final": final,
                "host_final": host_final,
                "analysis": analysis,
            })
            validate_raw(raw, output, check_files=True, check_current_inputs=True)
            raw_path = output / "raw.json"
            write_json_exclusive(raw_path, raw)
            manifest = signed({
                "schema": MANIFEST_SCHEMA,
                "created_utc": utc_now(),
                "valid": True,
                "validity_is_independent_of_speed": True,
                "raw": {
                    "path": "raw.json",
                    "size": raw_path.stat().st_size,
                    "sha256": sha256_file(raw_path),
                    "payload_digest": raw["digest"],
                },
                "campaign": campaign,
                "host": host_initial,
                "isolation": isolation,
                "reservation": reservation,
                "supervision": supervision,
                "identities": initial,
                "analysis": analysis,
            })
            write_json_exclusive(output / "manifest.json", manifest)
    except Exception as error:
        failure = signed({
            "schema": FAILURE_SCHEMA,
            "created_utc": utc_now(),
            "status": "failed",
            "valid": False,
            "error_type": type(error).__name__,
            "error": str(error),
            "campaign": campaign,
            "host_initial": host_initial,
            "reservation": reservation,
            "supervision": supervision,
            "pair_lease": pair_lease,
            "isolation": isolation,
            "input_specification": specification,
            "identities_initial": initial,
            "invocations": invocations,
            "retained_files": retained_file_records(output),
            "traceback": traceback.format_exc(),
        })
        failure_path = output / "failure.json"
        if not failure_path.exists():
            write_json_exclusive(failure_path, failure)
        validate_failure(failure, output, check_files=True)
        raise
    print(output / "manifest.json")
    return 0


def verified_campaign_bundle(
    manifest_path: Path, no_current_input_check: bool = False,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Verify and return exact evidence plus the accepted manifest identity."""
    manifest_path = manifest_path.resolve(strict=True)
    manifest, manifest_bytes = json_file_snapshot(manifest_path, "manifest")
    manifest_snapshot = {
        "size": len(manifest_bytes),
        "sha256": sha256_bytes(manifest_bytes),
    }
    verify_signature(manifest, "manifest")
    manifest_schema = manifest.get("schema")
    require(isinstance(manifest_schema, str) and
            manifest_schema in MANIFEST_TO_RAW_SCHEMA and
            manifest.get("valid") is True,
            "manifest is not valid main-comparison evidence")
    output = manifest_path.parent
    raw_info = manifest.get("raw")
    require(isinstance(raw_info, dict), "manifest has no raw bundle identity")
    raw_path = safe_evidence_path(output, raw_info.get("path"))
    require(raw_path.is_file(), "retained raw bundle is missing")
    raw, raw_bytes = json_file_snapshot(raw_path, "retained raw bundle")
    require(raw_info.get("size") == len(raw_bytes) and
            raw_info.get("sha256") == sha256_bytes(raw_bytes),
            "raw bundle file identity mismatch")
    expected_raw_schema = MANIFEST_TO_RAW_SCHEMA[manifest_schema]
    require(raw.get("schema") == expected_raw_schema,
            "manifest/raw schema versions do not match")
    analysis = validate_raw(
        raw, output, check_files=True,
        check_current_inputs=not no_current_input_check)
    require(raw_info.get("payload_digest") == raw.get("digest"),
            "manifest/raw payload identity mismatch")
    names = ["campaign", "host", "reservation", "identities", "analysis"]
    if manifest_schema in (
        MANIFEST_SCHEMA_V2, MANIFEST_SCHEMA_V3, MANIFEST_SCHEMA_V4,
        MANIFEST_SCHEMA_V5, MANIFEST_SCHEMA
    ):
        names.append("isolation")
    else:
        require("isolation" not in manifest,
                "legacy manifest contains unversioned isolation evidence")
    if manifest_schema in (MANIFEST_SCHEMA_V5, MANIFEST_SCHEMA):
        names.append("supervision")
    else:
        require("supervision" not in manifest,
                "historical manifest contains unversioned supervision data")
    for name in names:
        if name == "identities":
            expected = raw["identities_initial"]
        elif name == "host":
            expected = raw["host_initial"]
        else:
            expected = raw[name]
        require(manifest.get(name) == expected,
                f"manifest {name} differs from retained raw bundle")
    require(manifest.get("analysis") == analysis, "manifest analysis was edited")
    return manifest, raw, analysis, manifest_snapshot


def verify_campaign(options: argparse.Namespace) -> int:
    manifest_path = options.manifest.resolve(strict=True)
    manifest, _, _, manifest_snapshot = verified_campaign_bundle(
        manifest_path, options.no_current_input_check)
    manifest_schema = manifest["schema"]
    affinity_binding = getattr(options, "affinity_binding", None)
    if affinity_binding is not None:
        affinity_binding = Path(os.path.abspath(os.fspath(affinity_binding)))
        supervisor = Path(__file__).resolve().parents[3] / \
            "tools/leopard2_affinity_supervisor.py"
        require(supervisor.is_file(), "affinity binding verifier is missing")
        completed = run_process_bounded(
            [sys.executable, "-I", "-S", str(supervisor),
             "verify-binding", "--binding",
             str(affinity_binding), "--manifest", str(manifest_path),
             "--manifest-sha256", manifest_snapshot["sha256"]],
            timeout=30.0, max_stdout=1024 * 1024, max_stderr=1024 * 1024)
        require(completed.returncode == 0,
                "affinity binding verification failed: {}".format(
                    completed.stderr.decode("utf-8", errors="replace")))
    if manifest_schema == MANIFEST_SCHEMA_V1:
        print("legacy exact-main v1 bundle verified; it has no v2 CPU-isolation "
              "qualification")
    elif options.no_current_input_check:
        print("exact-main ABBA bundle structurally verified only; current build/source "
              "closure was not revalidated")
    else:
        print("exact-main ABBA evidence and current build/source closure verified")
    return 0


def verify_failed_campaign(options: argparse.Namespace) -> int:
    failure_path = options.failure.resolve(strict=True)
    try:
        failure = json.loads(failure_path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"failed campaign JSON is invalid: {error}") from error
    validate_failure(failure, failure_path.parent, check_files=True)
    print("failed exact-main campaign diagnostics and retained files verified; "
          "this is not valid performance evidence")
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="run a new, isolated three-round ABBA campaign")
    run.add_argument("--baseline", required=True, type=Path)
    run.add_argument("--candidate", required=True, type=Path)
    run.add_argument("--baseline-archive", required=True, type=Path)
    run.add_argument("--candidate-archive", required=True, type=Path)
    run.add_argument("--baseline-build-dir", required=True, type=Path)
    run.add_argument("--candidate-build-dir", required=True, type=Path)
    run.add_argument("--baseline-source-root", required=True, type=Path)
    run.add_argument("--candidate-source-root", required=True, type=Path)
    run.add_argument("--candidate-commit", required=True)
    run.add_argument("--candidate-mode", choices=CANDIDATE_MODES,
                     default="auto",
                     help="candidate decoder path: auto, generic, materialized, or tiled")
    run.add_argument("--reservation-file", required=True, type=Path)
    run.add_argument("--output", required=True, type=Path)
    run.add_argument("--cpu", required=True, type=int)
    run.add_argument("--reserved-sibling", required=True, type=int)
    run.add_argument("--taskset", default="/usr/bin/taskset", type=Path)
    run.add_argument("--ldd", default="/usr/bin/ldd", type=Path)
    run.add_argument("--preset", choices=("representative", "smoke"),
                     default="representative")
    run.add_argument("--cell", action="append",
                     help="override preset with ID:K:R:BYTES:LOSSES:SEED")
    run.add_argument("--reuse", type=int, default=8)
    run.add_argument("--iterations", type=int, default=9)
    run.add_argument("--warmup", type=int, default=2)
    run.add_argument("--timeout", type=float, default=120.0)
    run.set_defaults(function=run_campaign)
    verify = commands.add_parser("verify", help="verify retained evidence and raw outputs")
    verify.add_argument("--manifest", required=True, type=Path)
    verify.add_argument("--no-current-input-check", action="store_true",
                        help="structural-only replay without revalidating original build paths")
    verify.add_argument(
        "--affinity-binding", type=Path,
        help="also require a verified accepted affinity-supervisor binding")
    verify.set_defaults(function=verify_campaign)
    verify_failure = commands.add_parser(
        "verify-failure", help="verify a retained failed campaign bundle")
    verify_failure.add_argument("--failure", required=True, type=Path)
    verify_failure.set_defaults(function=verify_failed_campaign)
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    options = parser().parse_args(arguments)
    try:
        return int(options.function(options))
    except EvidenceError as error:
        print(f"main comparison evidence error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
