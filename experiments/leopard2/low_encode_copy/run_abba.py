#!/usr/bin/env python3
"""Provenance-bound isolated A/B evidence for the low-encode copy removal.

Exit status 0 means structurally valid evidence passed the promotion policy,
2 means structurally valid evidence failed that policy, and 1 means the request,
run, or retained evidence was invalid.  A policy failure is not an integrity
failure and is retained as reviewable negative evidence.
"""

from __future__ import annotations

import argparse
import copy
import ctypes
import datetime as dt
import errno
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import re
import selectors
import shlex
import shutil
import signal
import socket
import stat
import statistics
import subprocess
import sys
import tempfile
import time
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


CONTROL_COMMIT = "4070e4e527935026fb87593567587558f0a08d51"
CANDIDATE_COMMIT = "6d3afee213b94d486cf5f8145ac18078883ebc20"
RAW_SCHEMA = "leopard2-low-encode-copy-raw/v3"
MANIFEST_SCHEMA = "leopard2-low-encode-copy-manifest/v3"
FAILURE_SCHEMA = "leopard2-low-encode-copy-failure/v3"
EXECUTION_SCHEMA = "leopard2-low-encode-copy-execution/v1"
RESERVATION_SCHEMA = "leopard2-cpu-reservation/v1"
ROUNDS = 5
SEQUENCE = (
    ("A1", "control"),
    ("B1", "candidate"),
    ("B2", "candidate"),
    ("A2", "control"),
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
MAX_JSON_BYTES = 64 * 1024 * 1024
MAX_STDOUT_BYTES = 8 * 1024 * 1024
MAX_STDERR_BYTES = 1024 * 1024
MAX_JSON_DEPTH = 48
MAX_JSON_NODES = 300_000
MAX_JSON_STRING_BYTES = 8 * 1024 * 1024
MAX_BUILD_METADATA_BYTES = 16 * 1024 * 1024
MAX_BUILD_ARTIFACT_BYTES = 128 * 1024 * 1024
MAX_INVENTORY_ENTRIES = 2048
MAX_JSON_INTEGER = (1 << 64) - 1
MAX_JSON_FLOAT = 1.0e15
MAX_ITERATIONS = 1000
MAX_REUSE = 1_000_000_000
MAX_WARMUP = 1000
MAX_TIMEOUT_SECONDS = 3600.0
MAX_PATH_BYTES = 4096
MAX_PATH_COMPONENT_BYTES = 255
HEX64 = re.compile(r"^[0-9a-f]{16}$")
HEX256 = re.compile(r"^[0-9a-f]{64}$")
TARGET_SPEEDUP = 1.05
NEIGHBOR_REGRESSION_FLOOR = 1.0 / 1.02
_RELINK_CACHE: dict[str, str] = {}
_RECOMPILE_CACHE: dict[str, dict[str, Any]] = {}

# These vectors were generated from the independently compiled C++
# bench/leopard2/benchmark.cpp XorShift64 implementation.  They deliberately
# do not call the Python mock, so copying the same wrong PRNG into the runner
# and mock cannot make the self-test pass again.
LOSS_VECTOR_FIXTURES = (
    (0x26100, 8, 1, (6,)),
    (0x26101, 16, 1, (1,)),
    (0x26102, 32, 1, (29,)),
    (0x26103, 64, 1, (33,)),
    (0x26104, 100, 1, (65,)),
    (0x26105, 127, 1, (83,)),
    (0x26106, 128, 1, (52,)),
    (0x26107, 129, 1, (125,)),
    (0x123456, 17, 4, (1, 8, 10, 11)),
    (0x26117, 129, 8, (16, 24, 48, 57, 93, 110, 113, 128)),
)

FAILURE_PHASES = (
    "initialized",
    "topology_validated",
    "locks_acquired",
    "affinity_isolated",
    "inputs_attested",
    "executables_staged",
    "measurement_started",
    "measurement_completed",
    "evidence_prepared",
)


def _load_support() -> tuple[Any, Path]:
    path = Path(__file__).resolve().parents[1] / "main_compare" / "run_abba.py"
    require_path = path.resolve(strict=True)
    specification = importlib.util.spec_from_file_location(
        "leopard2_exact_main_evidence_support", require_path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load exact-main evidence support {require_path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module, require_path


SUPPORT, SUPPORT_PATH = _load_support()
EvidenceError = SUPPORT.EvidenceError


@dataclass(frozen=True)
class Cell:
    identifier: str
    backend: str
    field: str
    k: int
    r: int
    shard_bytes: int
    role: str
    seed: int
    losses: int = 1


SHAPES = (
    ("gf8-target-p8", "gf8", 8, 248, 64, "target"),
    ("gf8-neighbor-p16", "gf8", 16, 240, 256, "neighbor"),
    ("gf8-target-p32", "gf8", 32, 224, 1024, "target"),
    ("gf8-target-p64", "gf8", 64, 192, 4096, "target"),
    ("gf16-target-p128", "gf16", 100, 156, 16384, "target"),
    ("gf16-neighbor-rgt-k", "gf16", 127, 129, 65536, "neighbor"),
    ("gf16-neighbor-balanced", "gf16", 128, 128, 262144, "neighbor"),
    ("gf16-neighbor-p-jump", "gf16", 129, 127, 1048576, "neighbor"),
)


def fixed_cells() -> tuple[Cell, ...]:
    result: list[Cell] = []
    seed = 0x26_100
    for backend in ("scalar", "ssse3", "avx2"):
        for name, field, k, r, shard_bytes, role in SHAPES:
            result.append(Cell(
                identifier=f"{backend}-{name}-{shard_bytes}",
                backend=backend,
                field=field,
                k=k,
                r=r,
                shard_bytes=shard_bytes,
                role=role,
                seed=seed,
            ))
            seed += 1
    return tuple(result)


FIXED_CELLS = fixed_cells()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def exact_json_equal(actual: object, expected: object) -> bool:
    if isinstance(expected, bool):
        return isinstance(actual, bool) and actual is expected
    if isinstance(expected, int):
        return isinstance(actual, int) and not isinstance(actual, bool) and actual == expected
    if isinstance(expected, float):
        return isinstance(actual, float) and actual == expected
    if isinstance(expected, str) or expected is None:
        return type(actual) is type(expected) and actual == expected
    if isinstance(expected, list):
        return isinstance(actual, list) and len(actual) == len(expected) and all(
            exact_json_equal(left, right) for left, right in zip(actual, expected))
    if isinstance(expected, dict):
        return isinstance(actual, dict) and set(actual) == set(expected) and all(
            exact_json_equal(actual[name], expected[name]) for name in expected)
    return type(actual) is type(expected) and actual == expected


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def validate_utc(value: object, what: str) -> str:
    require(isinstance(value, str) and 20 <= len(value) <= 40 and value.endswith("Z"),
            f"{what} is not a bounded UTC timestamp")
    try:
        parsed = dt.datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise EvidenceError(f"{what} is not an ISO-8601 timestamp: {error}") from error
    require(parsed.tzinfo is not None and parsed.utcoffset() == dt.timedelta(0),
            f"{what} is not UTC")
    return value


def canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise EvidenceError(f"value is not canonical JSON: {error}") from error


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def open_regular_readonly(
    path: Path, what: str, limit: int, *, require_single_link: bool = True
) -> tuple[int, os.stat_result]:
    """Open a bounded regular file without ever blocking on a special file."""
    require(type(limit) is int and 0 <= limit <= MAX_BUILD_ARTIFACT_BYTES,
            f"{what} byte limit is invalid")
    try:
        before = os.lstat(path)
    except OSError as error:
        raise EvidenceError(f"cannot inspect {what}: {error}") from error
    require(stat.S_ISREG(before.st_mode) and
            (not require_single_link or before.st_nlink == 1) and
            0 <= before.st_size <= limit,
            f"{what} is not a bounded " +
            ("single-link " if require_single_link else "") + "regular file")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise EvidenceError(f"cannot safely open {what}: {error}") from error
    try:
        opened = os.fstat(descriptor)
        after = os.lstat(path)
        require(stat.S_ISREG(opened.st_mode) and
                (not require_single_link or opened.st_nlink == 1) and
                (opened.st_dev, opened.st_ino) ==
                (before.st_dev, before.st_ino) == (after.st_dev, after.st_ino) and
                opened.st_size == before.st_size == after.st_size and
                opened.st_size <= limit,
                f"{what} changed type, identity, links, or size while opening")
        return descriptor, opened
    except Exception:
        os.close(descriptor)
        raise


def sha256_file(path: Path, limit: int = MAX_BUILD_ARTIFACT_BYTES) -> str:
    require(type(limit) is int and 0 <= limit <= MAX_BUILD_ARTIFACT_BYTES,
            "file digest limit is invalid")
    descriptor, initial = open_regular_readonly(
        path, "file digest input", limit, require_single_link=False)
    try:
        path_metadata = os.lstat(path)
        require(stat.S_ISREG(initial.st_mode) and
                (initial.st_dev, initial.st_ino) ==
                (path_metadata.st_dev, path_metadata.st_ino) and
                initial.st_size <= limit,
                "file digest input is not a bounded inode-bound regular file")
        digest, size = sha256_descriptor(descriptor, limit)
        final = os.fstat(descriptor)
        require(size == initial.st_size and
                (final.st_size, final.st_mtime_ns, final.st_ctime_ns) ==
                (initial.st_size, initial.st_mtime_ns, initial.st_ctime_ns),
                "file changed while hashing")
        return digest
    finally:
        os.close(descriptor)


def sha256_descriptor(descriptor: int, limit: int) -> tuple[str, int]:
    digest = hashlib.sha256()
    retained = 0
    offset = 0
    while True:
        block = os.pread(descriptor, min(1024 * 1024, limit + 1 - retained), offset)
        if not block:
            break
        digest.update(block)
        retained += len(block)
        offset += len(block)
        require(retained <= limit, f"file descriptor exceeds {limit} bytes")
    return digest.hexdigest(), retained


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


def validate_json_shape(value: object, what: str) -> None:
    nodes = 0
    stack: list[tuple[object, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        nodes += 1
        require(nodes <= MAX_JSON_NODES, f"{what} exceeds the JSON node limit")
        require(depth <= MAX_JSON_DEPTH, f"{what} exceeds the JSON depth limit")
        if isinstance(current, str):
            require(len(current.encode("utf-8")) <= MAX_JSON_STRING_BYTES,
                    f"{what} contains an oversized string")
        elif isinstance(current, list):
            stack.extend((item, depth + 1) for item in current)
        elif isinstance(current, dict):
            for key, item in current.items():
                require(isinstance(key, str), f"{what} contains a non-string key")
                require(len(key.encode("utf-8")) <= 4096,
                        f"{what} contains an oversized key")
                stack.append((item, depth + 1))
        else:
            require(current is None or isinstance(current, (bool, int, float)),
                    f"{what} contains a non-JSON value")
            if isinstance(current, int) and not isinstance(current, bool):
                require(abs(current) <= MAX_JSON_INTEGER,
                        f"{what} contains an oversized integer")
            if isinstance(current, float):
                require(math.isfinite(current) and abs(current) <= MAX_JSON_FLOAT,
                        f"{what} contains an invalid or oversized float")


def parse_json_bytes(value: bytes, what: str, limit: int = MAX_JSON_BYTES) -> object:
    require(len(value) <= limit, f"{what} exceeds {limit} bytes")
    try:
        decoded = value.decode("utf-8")
        result = json.loads(decoded)
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
        raise EvidenceError(f"{what} is not valid UTF-8 JSON: {error}") from error
    validate_json_shape(result, what)
    return result


def read_json_limited(path: Path, what: str, limit: int = MAX_JSON_BYTES) -> object:
    descriptor, metadata = open_regular_readonly(path, what, limit)
    try:
        path_metadata = os.lstat(path)
        require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
                (metadata.st_dev, metadata.st_ino) ==
                (path_metadata.st_dev, path_metadata.st_ino),
                f"{what} is not an inode-bound single-link regular file")
        require(metadata.st_size <= limit, f"{what} exceeds {limit} bytes")
        chunks: list[bytes] = []
        retained = 0
        while True:
            block = os.read(descriptor, min(1024 * 1024, limit + 1 - retained))
            if not block:
                break
            chunks.append(block)
            retained += len(block)
            require(retained <= limit, f"{what} exceeds {limit} bytes")
        final = os.fstat(descriptor)
        final_path = os.lstat(path)
        require(stat.S_ISREG(final.st_mode) and final.st_nlink == 1 and
                stat.S_ISREG(final_path.st_mode) and final_path.st_nlink == 1 and
                (final.st_dev, final.st_ino) ==
                (metadata.st_dev, metadata.st_ino) ==
                (final_path.st_dev, final_path.st_ino) and
                (final.st_size, final.st_mtime_ns, final.st_ctime_ns) ==
                (metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns),
                f"{what} changed while it was read")
        return parse_json_bytes(b"".join(chunks), what, limit)
    finally:
        os.close(descriptor)


def write_bytes_exclusive(path: Path, value: bytes) -> None:
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    try:
        view = memoryview(value)
        while view:
            count = os.write(descriptor, view)
            require(count > 0, f"short write to {path}")
            view = view[count:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_json_exclusive(path: Path, value: object) -> None:
    write_bytes_exclusive(path, canonical_bytes(value) + b"\n")


def validate_path_text(value: object, what: str, absolute: bool | None = None) -> str:
    require(isinstance(value, str) and value and "\0" not in value and
            len(value.encode("utf-8")) <= MAX_PATH_BYTES,
            f"{what} is not a bounded path")
    path = Path(value)
    require(absolute is None or path.is_absolute() is absolute,
            f"{what} has the wrong absolute/relative form")
    require(len(path.parts) <= 128 and all(
        len(part.encode("utf-8")) <= MAX_PATH_COMPONENT_BYTES for part in path.parts),
        f"{what} has too many or oversized components")
    return value


def top_level_evidence_path(path: Path, what: str) -> Path:
    """Reject a non-regular top-level input without following that entry."""
    # abspath removes lexical '.'/'..' components but, unlike resolve(), never
    # follows the top-level directory entry.  The subsequent descriptor open
    # uses O_NOFOLLOW and rebinds the inode, covering a replacement race.
    absolute = Path(os.path.abspath(os.fspath(path)))
    validate_path_text(str(absolute), what, absolute=True)
    try:
        metadata = os.lstat(absolute)
    except OSError as error:
        raise EvidenceError(f"cannot inspect {what}: {error}") from error
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
            f"{what} is not a single-link regular file")
    return absolute


def safe_evidence_path(root: Path, relative: object) -> Path:
    relative = validate_path_text(relative, "evidence path", absolute=False)
    require("\\" not in relative, "evidence path contains a backslash")
    canonical_root = root.resolve(strict=True)
    candidate = canonical_root.joinpath(relative)
    require(all(part not in ("", ".", "..") for part in Path(relative).parts),
            "evidence path contains an unsafe component")
    resolved = candidate.resolve(strict=True)
    try:
        resolved.relative_to(canonical_root)
    except ValueError as error:
        raise EvidenceError(f"evidence path escapes output directory: {relative}") from error
    require(resolved == candidate,
            f"evidence path contains a symlink or noncanonical component: {relative}")
    return resolved


def artifact_record(
    root: Path, path: Path, limit: int = MAX_JSON_BYTES
) -> dict[str, Any]:
    path = path.resolve(strict=True)
    path.relative_to(root.resolve())
    metadata = os.lstat(path)
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
            metadata.st_size <= limit,
            f"retained artifact is not a single-link regular file: {path}")
    return {
        "path": str(path.relative_to(root.resolve())),
        "sha256": sha256_file(path, limit),
        "size": metadata.st_size,
    }


def verify_artifact(
    root: Path, value: object, what: str, limit: int = MAX_JSON_BYTES
) -> Path:
    require(isinstance(value, dict) and set(value) == {"path", "sha256", "size"},
            f"{what} identity is incomplete")
    require(isinstance(value.get("sha256"), str) and
            HEX256.fullmatch(value["sha256"]) is not None and
            isinstance(value.get("size"), int) and not isinstance(value["size"], bool) and
            0 <= value["size"] <= limit,
            f"{what} identity is invalid")
    path = safe_evidence_path(root, value.get("path"))
    descriptor, metadata = open_regular_readonly(path, what, limit)
    try:
        path_metadata = os.lstat(path)
        digest, size = sha256_descriptor(descriptor, limit)
        final = os.fstat(descriptor)
        final_path = os.lstat(path)
        require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
                (metadata.st_dev, metadata.st_ino) ==
                (path_metadata.st_dev, path_metadata.st_ino) and
                stat.S_ISREG(final.st_mode) and final.st_nlink == 1 and
                stat.S_ISREG(final_path.st_mode) and final_path.st_nlink == 1 and
                (metadata.st_dev, metadata.st_ino) ==
                (final.st_dev, final.st_ino) ==
                (final_path.st_dev, final_path.st_ino) and
                (metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns) ==
                (final.st_size, final.st_mtime_ns, final.st_ctime_ns) and
                size == metadata.st_size == value["size"] and
                digest == value["sha256"],
                f"{what} is missing or changed")
    finally:
        os.close(descriptor)
    return path


def ceil_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def validate_cell(cell: Cell) -> None:
    require(re.fullmatch(r"[a-z0-9][a-z0-9-]{0,95}", cell.identifier) is not None,
            f"invalid cell identifier {cell.identifier!r}")
    require(cell.backend in {"scalar", "ssse3", "avx2"},
            f"invalid backend in {cell.identifier}")
    require(cell.field in {"gf8", "gf16"} and cell.role in {"target", "neighbor"},
            f"invalid field/role in {cell.identifier}")
    for name in ("k", "r", "shard_bytes", "seed", "losses"):
        item = getattr(cell, name)
        require(isinstance(item, int) and not isinstance(item, bool),
                f"{name} is not an exact integer in {cell.identifier}")
    require(cell.k > 0 and cell.r > 0 and 0 <= cell.losses <= min(cell.k, cell.r),
            f"invalid counts in {cell.identifier}")
    require(cell.shard_bytes in {64, 256, 1024, 4096, 16384, 65536, 262144, 1048576},
            f"invalid aligned size in {cell.identifier}")
    padded = ceil_power_of_two(cell.k)
    parent = ceil_power_of_two(padded + cell.r)
    require(parent <= 65536 and
            ((cell.field == "gf8" and parent <= 256) or
             cell.field == "gf16"),
            f"field/parent mismatch in {cell.identifier}")
    require(0 <= cell.seed < (1 << 64), f"seed overflow in {cell.identifier}")


def statistics_policy() -> dict[str, Any]:
    return {
        "method": "mean log contrast per independent A1/B1/B2/A2 round",
        "confidence": 0.95,
        "critical_distribution": "Student-t",
        "independent_round_count_per_metric": ROUNDS,
        "degrees_of_freedom": ROUNDS - 1,
        "constituent_pair_count_per_metric": 2 * ROUNDS,
        "child_estimator": "median of retained benchmark-v2 samples",
        "setup_execution_separation": (
            "codec_setup and encode_execution are measured and analyzed separately"
        ),
        "promotion": {
            "target_encode_ci95_lower": TARGET_SPEEDUP,
            "neighbor_encode_ci95_lower": NEIGHBOR_REGRESSION_FLOOR,
            "setup_is_reported_but_not_a_promotion_gate": True,
        },
    }


def finite_number(value: object, what: str, positive: bool = True) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{what} is not numeric")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise EvidenceError(f"{what} is not a bounded number") from error
    require(math.isfinite(result) and abs(result) <= MAX_JSON_FLOAT,
            f"{what} is not finite or bounded")
    require(result >= 1.0e-12 if positive else result >= 0,
            f"{what} is {'not positive' if positive else 'negative'}")
    return result


def close_enough(actual: float, expected: float) -> bool:
    return abs(actual - expected) <= max(0.000002, abs(expected) * 0.000002)


def validate_summary(value: object, iterations: int, setup: bool) -> list[float]:
    require(isinstance(value, dict), "timing summary is not an object")
    suffix = "" if setup else "_per_batch_call"
    names = {
        "median": f"median_us{suffix}",
        "mad": f"mad_us{suffix}",
        "minimum": f"minimum_us{suffix}",
        "maximum": f"maximum_us{suffix}",
        "samples": "samples_us" if setup else "samples_us_per_batch_call",
    }
    statistic_keys = set(names.values())
    if setup:
        require(set(value) == statistic_keys,
                "setup timing summary has unexpected or missing fields")
    retained = value.get(names["samples"])
    require(isinstance(retained, list) and len(retained) == iterations,
            f"{names['samples']} does not contain exactly {iterations} samples")
    samples = [finite_number(item, names["samples"]) for item in retained]
    derived_median = statistics.median(samples)
    expected = {
        names["median"]: derived_median,
        names["mad"]: statistics.median(
            [abs(item - derived_median) for item in samples]),
        names["minimum"]: min(samples),
        names["maximum"]: max(samples),
    }
    for name, derived in expected.items():
        actual = finite_number(value.get(name), name, positive=name != names["mad"])
        require(close_enough(actual, derived), f"{name} is not derived from raw samples")
    return samples


def validate_execution_summary(
    value: object,
    iterations: int,
    input_name: str,
    output_name: str,
    input_bytes: int,
    output_bytes: int,
) -> list[float]:
    require(isinstance(value, dict) and set(value) == {
        "mad_us_per_batch_call", "maximum_us_per_batch_call",
        "median_us_per_batch_call", "minimum_us_per_batch_call",
        "samples_us_per_batch_call", input_name, output_name},
        "execution timing summary has unexpected or missing fields")
    samples = validate_summary(value, iterations, False)
    median = statistics.median(samples)
    for name, byte_count in ((input_name, input_bytes), (output_name, output_bytes)):
        actual = finite_number(value.get(name), name)
        expected = byte_count / (median * 1000.0)
        require(close_enough(actual, expected),
                f"{name} is not derived from bytes and retained median")
    return samples


def validate_support_contract() -> dict[str, Any]:
    required_callables = (
        "artifact_identity", "build_provenance", "canonical_bytes",
        "cpu_stat_snapshot", "git_identity", "host_identity", "isolation_record",
        "PairLease", "parse_cpu_list", "runtime_closure", "validate_host_record",
        "validate_isolation", "validate_topology", "run_process_bounded",
        "terminate_process_group_bounded",
    )
    for name in required_callables:
        require(callable(getattr(SUPPORT, name, None)),
                f"exact-main evidence support lacks callable {name}")
    require(SUPPORT.RESERVATION_SCHEMA == RESERVATION_SCHEMA,
            "support reservation schema changed")
    require(SUPPORT.PAIR_LEASE_SCHEMA == "leopard2-cpu-pair-lease/v1",
            "support pair-lease schema changed")
    require(tuple(SUPPORT.CPU_STAT_FIELDS) == (
        "user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal"),
        "support /proc/stat field contract changed")
    require(tuple(SUPPORT.CPU_STAT_IDLE_FIELDS) == ("idle", "iowait") and
            SUPPORT.MAX_SIBLING_NONIDLE_JIFFIES == 0,
            "support sibling-idle contract changed")
    require(SUPPORT.MAX_CPU_ID == 1_048_575 and
            SUPPORT.MAX_CPU_LIST_ENTRIES == 4096 and
            SUPPORT.MAX_CPU_LIST_TEXT_BYTES == 65_536,
            "support bounded CPU-list contract changed")
    require(SUPPORT.MAX_COMMAND_STDOUT_BYTES == 128 * 1024 * 1024 and
            SUPPORT.MAX_COMMAND_STDERR_BYTES == 8 * 1024 * 1024 and
            SUPPORT.MAX_COMMAND_TIMEOUT_SECONDS == 3600.0,
            "support bounded subprocess contract changed")
    probe = {"b": 2, "a": 1}
    require(SUPPORT.canonical_bytes(probe) == canonical_bytes(probe),
            "support canonical JSON contract changed")
    payload = {
        "benchmark_cpu": 3,
        "nonce": "0123456789abcdef0123456789abcdef",
        "owner": "self-test",
        "reserved_sibling": 7,
        "schema": RESERVATION_SCHEMA,
        "status": "held",
    }
    require(SUPPORT.parse_reservation(canonical_bytes(payload), 3, 7) == payload,
            "support reservation semantic contract changed")
    return {
        "path": str(SUPPORT_PATH),
        "sha256": sha256_file(SUPPORT_PATH),
        "contracts": {
            "canonical_json": "sort_keys compact UTF-8 without newline",
            "cpu_stat_fields": list(SUPPORT.CPU_STAT_FIELDS),
            "max_cpu_id": SUPPORT.MAX_CPU_ID,
            "max_cpu_list_entries": SUPPORT.MAX_CPU_LIST_ENTRIES,
            "pair_lease_schema": SUPPORT.PAIR_LEASE_SCHEMA,
            "reservation_schema": SUPPORT.RESERVATION_SCHEMA,
            "sibling_max_nonidle_jiffies": SUPPORT.MAX_SIBLING_NONIDLE_JIFFIES,
            "subprocess_stderr_limit": SUPPORT.MAX_COMMAND_STDERR_BYTES,
            "subprocess_stdout_limit": SUPPORT.MAX_COMMAND_STDOUT_BYTES,
            "subprocess_timeout_limit": SUPPORT.MAX_COMMAND_TIMEOUT_SECONDS,
        },
    }


def side_build_specification(side: str, specification: Mapping[str, Any]) -> dict[str, Any]:
    root = specification[f"{side}_source_root"]
    return {
        "baseline_source_root": root,
        "candidate_source_root": root,
        "candidate_executable": specification[f"{side}_executable"],
        "candidate_archive": specification[f"{side}_archive"],
        "candidate_build_dir": specification[f"{side}_build_dir"],
    }


def require_bounded_regular(path: Path, what: str, limit: int) -> None:
    metadata = os.lstat(path)
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
            0 <= metadata.st_size <= limit,
            f"{what} is not a bounded single-link regular file")


def strict_compile_recipes(
    side: str, specification: Mapping[str, Any], provenance: Mapping[str, Any]
) -> dict[str, list[str]]:
    root = Path(specification[f"{side}_source_root"]).resolve(strict=True)
    build = Path(specification[f"{side}_build_dir"]).resolve(strict=True)
    compiler = Path(provenance["compiler"]["path"])
    commands_path = build / "compile_commands.json"
    entries = read_json_limited(
        commands_path, f"{side} compile commands", MAX_BUILD_METADATA_BYTES)
    require(isinstance(entries, list) and
            len(entries) == provenance["validated_compile_commands"]["entry_count"] and
            len(entries) <= 256,
            f"{side} compile-command entry count is invalid")
    by_source: dict[Path, Mapping[str, Any]] = {}
    for entry in entries:
        require(isinstance(entry, dict) and set(entry) == {
            "command", "directory", "file", "output"},
            f"{side} compile command has unexpected fields")
        source = Path(entry["file"]).resolve(strict=True)
        require(source not in by_source, f"{side} duplicate compile command for {source}")
        by_source[source] = entry

    allowed_literals = {
        "-Wall", "-Wextra", "-fopenmp", "-g", "-O0", "-O3",
        "-std=gnu++11", "-mavx2", "-mno-avx512f",
        "-falign-functions=64", "-mssse3", "-mno-avx",
        "-DLEO2_DISABLE_AVX2_CODEGEN=1", "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
        "-DLEO2_HAVE_AVX2_BACKEND=1", "-DLEO2_HAVE_SSSE3_BACKEND=1",
        "-c", "-o",
    }
    result: dict[str, list[str]] = {}
    pairs = provenance["validated_compile_commands"][
        "required_source_object_pairs"]
    for pair in pairs:
        source = Path(pair["source"]["path"])
        obj = Path(pair["object"]["path"])
        entry = by_source.get(source)
        require(entry is not None and
                Path(entry["directory"]).resolve(strict=True) == build,
                f"{side} compile command has the wrong build directory")
        tokens = SUPPORT.compile_command_tokens(entry)
        require(tokens and Path(tokens[0]).resolve(strict=True) == compiler,
                f"{side} compile command uses another compiler")
        normalized = ["<compiler>"]
        index = 1
        saw_source = False
        saw_output = False
        while index < len(tokens):
            token = tokens[index]
            if token == "-o":
                require(not saw_output and index + 1 < len(tokens),
                        f"{side} compile command has an invalid output")
                candidate = Path(tokens[index + 1])
                if not candidate.is_absolute():
                    candidate = build / candidate
                require(candidate.resolve(strict=True) == obj,
                        f"{side} compile command output differs from its object")
                normalized.extend(("-o", "<object>"))
                saw_output = True
                index += 2
                continue
            if token == str(source):
                require(not saw_source, f"{side} compile source is duplicated")
                normalized.append("<source>")
                saw_source = True
            elif token == f"-I{root}":
                normalized.append("-I<SOURCE_ROOT>")
            else:
                require(token in allowed_literals,
                        f"{side} compile command has a noncanonical token: {token!r}")
                normalized.append(token)
            index += 1
        require(saw_source and saw_output and normalized.count("-c") == 1 and
                normalized.count("-O3") == 1,
                f"{side} compile command is incomplete")
        relative = str(source.relative_to(root))
        require(entry["output"] and
                (build / entry["output"]).resolve(strict=True) == obj,
                f"{side} compile-command output metadata differs")
        result[relative] = normalized
    require(len(result) == len(pairs), f"{side} strict compile closure is incomplete")
    return dict(sorted(result.items()))


def clean_recompile_proof(
    side: str, specification: Mapping[str, Any], provenance: Mapping[str, Any],
    recipes: Mapping[str, Sequence[str]],
) -> list[dict[str, Any]]:
    root = Path(specification[f"{side}_source_root"]).resolve(strict=True)
    build = Path(specification[f"{side}_build_dir"]).resolve(strict=True)
    entries = read_json_limited(
        build / "compile_commands.json", f"{side} compile commands",
        MAX_BUILD_METADATA_BYTES)
    require(isinstance(entries, list), f"{side} compile commands are not a list")
    by_source: dict[Path, Mapping[str, Any]] = {}
    for entry in entries:
        require(isinstance(entry, dict) and isinstance(entry.get("file"), str),
                f"{side} clean-recompile entry is invalid")
        source = Path(entry["file"]).resolve(strict=True)
        require(source not in by_source,
                f"{side} clean-recompile source is duplicated")
        by_source[source] = entry
    result: list[dict[str, Any]] = []
    pairs = provenance["validated_compile_commands"]["required_source_object_pairs"]
    with tempfile.TemporaryDirectory(prefix=f"leo2-{side}-recompile-") as temporary:
        for index, pair in enumerate(pairs):
            source = Path(pair["source"]["path"])
            relative = str(source.relative_to(root))
            entry = by_source.get(source)
            require(entry is not None, f"{side} clean-recompile source is missing")
            tokens = SUPPORT.compile_command_tokens(entry)
            output_positions = [
                position for position, token in enumerate(tokens) if token == "-o"
            ]
            require(len(output_positions) == 1 and output_positions[0] + 1 < len(tokens),
                    f"{side} clean-recompile output is invalid")
            cache_key = sha256_bytes(canonical_bytes({
                "compiler": provenance["compiler"],
                "normalized_recipe": recipes[relative],
                "object": pair["object"],
                "source": pair["source"],
            }))
            cached = _RECOMPILE_CACHE.get(cache_key)
            if cached is None:
                output = Path(temporary) / f"object-{index}.o"
                command = list(tokens)
                command[output_positions[0] + 1] = str(output)
                completed = SUPPORT.run_process_bounded(
                    command, cwd=build, environment=CHILD_ENVIRONMENT,
                    timeout=120.0, max_stdout=MAX_STDOUT_BYTES,
                    max_stderr=MAX_STDERR_BYTES)
                require(completed.returncode == 0,
                        f"{side} clean recompilation failed for {relative}: " +
                        completed.stderr[:4096].decode("utf-8", errors="replace"))
                require_bounded_regular(
                    output, f"{side} clean-recompiled object",
                    MAX_BUILD_ARTIFACT_BYTES)
                digest = sha256_file(output, MAX_BUILD_ARTIFACT_BYTES)
                size = os.lstat(output).st_size
                require(digest == pair["object"]["sha256"] and
                        size == pair["object"]["size"],
                        f"{side} object is not reproducible from its retained source/recipe: "
                        f"{relative}")
                cached = {
                    "cache_key": cache_key,
                    "object_sha256": digest,
                    "object_size": size,
                    "source": relative,
                    "status": "byte_identical_clean_recompile",
                }
                _RECOMPILE_CACHE[cache_key] = copy.deepcopy(cached)
            result.append(copy.deepcopy(cached))
    return result


def strict_archive_and_link_provenance(
    side: str, specification: Mapping[str, Any], provenance: Mapping[str, Any]
) -> dict[str, Any]:
    build = Path(specification[f"{side}_build_dir"]).resolve(strict=True)
    archive = Path(specification[f"{side}_archive"]).resolve(strict=True)
    executable = Path(specification[f"{side}_executable"]).resolve(strict=True)
    compiler = Path(provenance["compiler"]["path"])
    archive_recipe_path = Path(provenance["archive_link_recipe"]["path"])
    executable_recipe_path = Path(provenance["executable_link_recipe"]["path"])
    pairs = provenance["validated_compile_commands"][
        "required_source_object_pairs"]
    benchmark_pairs = [
        pair for pair in pairs
        if pair["source"]["path"].endswith("/bench/leopard2/benchmark.cpp")
    ]
    require(len(benchmark_pairs) == 1, f"{side} has no unique benchmark object")
    benchmark_object = Path(benchmark_pairs[0]["object"]["path"])
    archive_pairs = [pair for pair in pairs if pair not in benchmark_pairs]
    archive_objects = [Path(pair["object"]["path"]) for pair in archive_pairs]

    ar = Path(provenance["archiver"]["path"])
    ranlib = Path("/usr/bin/ranlib").resolve(strict=True)
    archive_lines = [
        shlex.split(line) for line in archive_recipe_path.read_text(
            encoding="utf-8").splitlines() if line.strip()
    ]
    require(len(archive_lines) == 2 and len(archive_lines[0]) >= 3 and
            Path(archive_lines[0][0]).resolve(strict=True) == ar and
            archive_lines[0][1] == "qc" and archive_lines[1] and
            Path(archive_lines[1][0]).resolve(strict=True) == ranlib,
            f"{side} archive recipe is not canonical ar/ranlib")

    def resolve_build_token(token: str) -> Path:
        path = Path(token)
        if not path.is_absolute():
            path = build / path
        return path.resolve(strict=True)

    recipe_archive_objects = [
        resolve_build_token(token) for token in archive_lines[0][3:]
    ]
    require(resolve_build_token(archive_lines[0][2]) == archive and
            len(recipe_archive_objects) == len(set(recipe_archive_objects)) and
            set(recipe_archive_objects) == set(archive_objects) and
            len(archive_lines[1]) == 2 and
            resolve_build_token(archive_lines[1][1]) == archive,
            f"{side} archive recipe differs from its exact object closure")
    members = provenance["validated_archive_members"]
    require(members == [path.name for path in recipe_archive_objects],
            f"{side} archive member order differs from its objects")
    pair_by_object = {
        Path(pair["object"]["path"]): pair for pair in archive_pairs
    }
    ordered_archive_pairs = [pair_by_object[path] for path in recipe_archive_objects]
    member_content: list[dict[str, Any]] = []
    for member, pair in zip(members, ordered_archive_pairs):
        content = SUPPORT.run_checked(
            (str(ar), "p", str(archive), member), environment=CHILD_ENVIRONMENT)
        require(len(content) <= MAX_BUILD_ARTIFACT_BYTES and
                sha256_bytes(content) == pair["object"]["sha256"] and
                len(content) == pair["object"]["size"],
                f"{side} archive member {member} differs from its object file")
        member_content.append({
            "member": member,
            "object_path": pair["object"]["path"],
            "sha256": pair["object"]["sha256"],
            "size": pair["object"]["size"],
        })

    link_lines = [
        line for line in executable_recipe_path.read_text(
            encoding="utf-8").splitlines() if line.strip()
    ]
    require(len(link_lines) == 1, f"{side} executable link recipe is not one command")
    link_tokens = shlex.split(link_lines[0])
    allowed_flags = {"-Wall", "-Wextra", "-fopenmp", "-g", "-O0", "-O3"}
    normalized: list[str] = ["<compiler>"]
    external_inputs: list[dict[str, Any]] = []
    output_index: int | None = None
    saw_object = False
    saw_archive = False
    require(link_tokens and Path(link_tokens[0]).resolve(strict=True) == compiler,
            f"{side} executable link recipe uses another compiler")
    index = 1
    while index < len(link_tokens):
        token = link_tokens[index]
        if token == "-o":
            require(output_index is None and index + 1 < len(link_tokens),
                    f"{side} executable link output is invalid")
            require(resolve_build_token(link_tokens[index + 1]) == executable,
                    f"{side} executable link output path differs")
            normalized.extend(("-o", "<executable>"))
            output_index = index + 1
            index += 2
            continue
        if token in allowed_flags:
            normalized.append(token)
            index += 1
            continue
        resolved = resolve_build_token(token)
        if resolved == benchmark_object:
            require(not saw_object, f"{side} benchmark object is duplicated")
            normalized.append("<benchmark-object>")
            saw_object = True
        elif resolved == archive:
            require(not saw_archive, f"{side} archive is duplicated")
            normalized.append("<archive>")
            saw_archive = True
        else:
            require_bounded_regular(resolved, f"{side} external link input",
                                    MAX_BUILD_ARTIFACT_BYTES)
            identity = SUPPORT.artifact_identity(resolved, "link_input")
            external_inputs.append(identity)
            normalized.append(f"<external:{identity['path']}>")
        index += 1
    require(output_index is not None and saw_object and saw_archive and
            normalized.count("-O3") == 1,
            f"{side} executable link recipe is incomplete")

    cache_key = sha256_bytes(canonical_bytes({
        "archive": provenance["validated_archive"],
        "benchmark_object": benchmark_pairs[0]["object"],
        "compiler": provenance["compiler"],
        "external_inputs": external_inputs,
        "normalized_link": normalized,
    }))
    expected_sha = _RELINK_CACHE.get(cache_key)
    if expected_sha is None:
        with tempfile.TemporaryDirectory(prefix="leo2-low-copy-relink-") as temporary:
            relinked = Path(temporary) / "bench_leopard2"
            command = list(link_tokens)
            command[output_index] = str(relinked)
            completed = SUPPORT.run_process_bounded(
                command, cwd=build, environment=CHILD_ENVIRONMENT,
                timeout=30.0, max_stdout=MAX_STDOUT_BYTES,
                max_stderr=MAX_STDERR_BYTES)
            require(completed.returncode == 0 and
                    len(completed.stdout) <= MAX_STDOUT_BYTES and
                    len(completed.stderr) <= MAX_STDERR_BYTES,
                    f"{side} reproducible relink failed: " +
                    completed.stderr[:4096].decode("utf-8", errors="replace"))
            require_bounded_regular(relinked, f"{side} relinked executable",
                                    MAX_BUILD_ARTIFACT_BYTES)
            expected_sha = sha256_file(relinked)
        _RELINK_CACHE[cache_key] = expected_sha
    validate_relinked_executable_sha(
        provenance["validated_executable"]["sha256"], expected_sha, side)
    return {
        "archive_member_content": member_content,
        "external_link_inputs": external_inputs,
        "normalized_executable_link": normalized,
        "ranlib": SUPPORT.artifact_identity(ranlib, "ranlib"),
        "relink_cache_key": cache_key,
        "reproducible_executable_sha256": expected_sha,
    }


def validate_relinked_executable_sha(
    declared_sha: object, relinked_sha: object, side: str
) -> None:
    require(isinstance(declared_sha, str) and isinstance(relinked_sha, str) and
            HEX256.fullmatch(declared_sha) is not None and
            HEX256.fullmatch(relinked_sha) is not None and
            declared_sha == relinked_sha,
            f"{side} executable bytes are not the result of the retained link inputs")


def validate_build_provenance(
    side: str, specification: Mapping[str, Any]
) -> dict[str, Any]:
    build = Path(specification[f"{side}_build_dir"]).resolve(strict=True)
    for relative in (
        "CMakeCache.txt", "compile_commands.json",
        "CMakeFiles/bench_leopard2.dir/link.txt",
        "CMakeFiles/libleopard.dir/link.txt",
    ):
        require_bounded_regular(build / relative, f"{side} {relative}",
                                MAX_BUILD_METADATA_BYTES)
    require_bounded_regular(Path(specification[f"{side}_executable"]),
                            f"{side} executable", MAX_BUILD_ARTIFACT_BYTES)
    require_bounded_regular(Path(specification[f"{side}_archive"]),
                            f"{side} archive", MAX_BUILD_ARTIFACT_BYTES)
    result = SUPPORT.build_provenance(
        "candidate", side_build_specification(side, specification))
    cache_path = Path(specification[f"{side}_build_dir"]) / "CMakeCache.txt"
    cache = SUPPORT.parse_cmake_cache(cache_path)
    require(cache.get("CMAKE_GENERATOR") == "Unix Makefiles",
            f"{side} build does not use the declared Unix Makefiles recipe")
    require(cache.get("CMAKE_EXPORT_COMPILE_COMMANDS") == "ON",
            f"{side} build did not retain compile commands")
    result["validated_cache"]["CMAKE_GENERATOR"] = cache["CMAKE_GENERATOR"]
    result["validated_cache"]["CMAKE_EXPORT_COMPILE_COMMANDS"] = cache[
        "CMAKE_EXPORT_COMPILE_COMMANDS"]
    for name in (
        "CMAKE_CXX_FLAGS", "CMAKE_EXE_LINKER_FLAGS",
        "CMAKE_EXE_LINKER_FLAGS_RELEASE", "CMAKE_STATIC_LINKER_FLAGS",
        "CMAKE_STATIC_LINKER_FLAGS_RELEASE",
    ):
        require(cache.get(name) == "", f"{side} build has noncanonical {name}")
        result["validated_cache"][name] = ""
    for pair in result["validated_compile_commands"]["required_source_object_pairs"]:
        require_bounded_regular(Path(pair["source"]["path"]), f"{side} source",
                                MAX_BUILD_METADATA_BYTES)
        require_bounded_regular(Path(pair["object"]["path"]), f"{side} object",
                                MAX_BUILD_ARTIFACT_BYTES)
    result["strict_compile_recipes"] = strict_compile_recipes(
        side, specification, result)
    result["clean_recompile_proof"] = clean_recompile_proof(
        side, specification, result, result["strict_compile_recipes"])
    result["strict_link_provenance"] = strict_archive_and_link_provenance(
        side, specification, result)
    return result


def validate_equal_build_policy(
    control_build: Mapping[str, Any], candidate_build: Mapping[str, Any]
) -> None:
    require(control_build["compiler"] == candidate_build["compiler"] and
            control_build["compiler_version_stdout"] ==
            candidate_build["compiler_version_stdout"],
            "control and candidate use different compiler identities")
    control_cache = control_build["validated_cache"]
    candidate_cache = candidate_build["validated_cache"]
    for name in (
        "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS_RELEASE",
        "CMAKE_GENERATOR", "CMAKE_EXPORT_COMPILE_COMMANDS", "ENABLE_OPENMP",
        "LEO2_BACKEND_VARIANT", "LEO2_BUILD_BENCHMARKS", "LEO2_BUILD_FUZZERS",
        "LEO2_BUILD_TESTS", "LEO2_ENABLE_CUDA", "CMAKE_CXX_FLAGS",
        "CMAKE_EXE_LINKER_FLAGS", "CMAKE_EXE_LINKER_FLAGS_RELEASE",
        "CMAKE_STATIC_LINKER_FLAGS", "CMAKE_STATIC_LINKER_FLAGS_RELEASE",
    ):
        require(control_cache.get(name) == candidate_cache.get(name),
                f"control/candidate CMake setting differs: {name}")
    require(control_build["strict_compile_recipes"] ==
            candidate_build["strict_compile_recipes"],
            "control/candidate effective compile recipes differ")
    for name in ("normalized_executable_link", "external_link_inputs", "ranlib"):
        require(control_build["strict_link_provenance"][name] ==
                candidate_build["strict_link_provenance"][name],
                f"control/candidate link provenance differs: {name}")


def input_snapshot(specification: Mapping[str, Any]) -> dict[str, Any]:
    require(Path(sys.executable).resolve(strict=True) ==
            Path(specification["python"]).resolve(strict=True),
            "current Python interpreter differs from the declared coordinator")
    support_contract = validate_support_contract()
    control_build = validate_build_provenance("control", specification)
    candidate_build = validate_build_provenance("candidate", specification)
    validate_equal_build_policy(control_build, candidate_build)
    ldd = Path(specification["ldd"])
    control_root = Path(specification["control_source_root"])
    candidate_root = Path(specification["candidate_source_root"])
    return {
        "runner": SUPPORT.artifact_identity(Path(specification["runner"]), "file"),
        "support": {
            "artifact": SUPPORT.artifact_identity(SUPPORT_PATH, "file"),
            "validated_contract": support_contract,
        },
        "taskset": SUPPORT.artifact_identity(
            Path(specification["taskset"]), "executable"),
        "ldd": SUPPORT.artifact_identity(ldd, "executable"),
        "git": SUPPORT.artifact_identity(
            Path(specification["git"]), "executable"),
        "python": SUPPORT.artifact_identity(
            Path(specification["python"]), "executable"),
        "control_executable": SUPPORT.artifact_identity(
            Path(specification["control_executable"]), "executable"),
        "candidate_executable": SUPPORT.artifact_identity(
            Path(specification["candidate_executable"]), "executable"),
        "control_archive": SUPPORT.artifact_identity(
            Path(specification["control_archive"]), "archive"),
        "candidate_archive": SUPPORT.artifact_identity(
            Path(specification["candidate_archive"]), "archive"),
        "control_build": control_build,
        "candidate_build": candidate_build,
        "control_runtime_closure": SUPPORT.runtime_closure(
            ldd, Path(specification["control_executable"])),
        "candidate_runtime_closure": SUPPORT.runtime_closure(
            ldd, Path(specification["candidate_executable"])),
        "control_source": SUPPORT.git_identity(
            control_root, CONTROL_COMMIT, True),
        "candidate_source": SUPPORT.git_identity(
            candidate_root, CANDIDATE_COMMIT, False),
    }


def validate_file_identity(
    value: object, what: str, expected_path: str | None = None,
    expected_kind: str | None = None,
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) ==
            {"kind", "mode", "mtime_ns", "path", "sha256", "size"},
            f"retained {what} identity is incomplete")
    path = value.get("path")
    require(isinstance(path, str), f"retained {what} path is invalid")
    validate_path_text(path, f"retained {what} path", absolute=True)
    require(Path(path).is_absolute() and
            (expected_path is None or path == str(Path(expected_path).resolve())),
            f"retained {what} path is invalid")
    require(isinstance(value.get("kind"), str) and value["kind"] and
            (expected_kind is None or value["kind"] == expected_kind),
            f"retained {what} kind is invalid")
    require(isinstance(value.get("size"), int) and not isinstance(value["size"], bool) and
            0 <= value["size"] <= MAX_BUILD_ARTIFACT_BYTES and
            isinstance(value.get("mode"), int) and not isinstance(value["mode"], bool) and
            0 <= value["mode"] <= 0o7777 and
            isinstance(value.get("mtime_ns"), int) and
            not isinstance(value["mtime_ns"], bool) and value["mtime_ns"] >= 0 and
            isinstance(value.get("sha256"), str) and
            HEX256.fullmatch(value["sha256"]) is not None,
            f"retained {what} metadata or digest is invalid")
    return value


def validate_source_identity(
    value: object, what: str, expected_path: str, expected_head: str,
    detached: bool,
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "detached", "head", "path", "tracked_status",
        "tracked_tree_listing_sha256", "tree"},
        f"retained {what} source identity is incomplete")
    require(value.get("path") == str(Path(expected_path).resolve()) and
            value.get("head") == expected_head and value.get("detached") is detached and
            value.get("tracked_status") == "clean" and
            isinstance(value.get("tree"), str) and
            re.fullmatch(r"[0-9a-f]{40}", value["tree"]) is not None and
            isinstance(value.get("tracked_tree_listing_sha256"), str) and
            HEX256.fullmatch(value["tracked_tree_listing_sha256"]) is not None,
            f"retained {what} source commit/tree/cleanliness is invalid")
    return value


def validate_clean_recompile_record(
    proof: object, source_name: str, pair: Mapping[str, Any],
    compiler: Mapping[str, Any], normalized_recipe: Sequence[str],
) -> dict[str, Any]:
    require(isinstance(proof, dict) and set(proof) == {
        "cache_key", "object_sha256", "object_size", "source", "status"},
        "retained clean-recompile record is invalid")
    expected_key = sha256_bytes(canonical_bytes({
        "compiler": compiler,
        "normalized_recipe": list(normalized_recipe),
        "object": pair["object"],
        "source": pair["source"],
    }))
    require(proof.get("source") == source_name and
            proof.get("status") == "byte_identical_clean_recompile" and
            proof.get("cache_key") == expected_key and
            proof.get("object_sha256") == pair["object"]["sha256"] and
            proof.get("object_size") == pair["object"]["size"],
            "retained clean-recompile digest is invalid")
    return proof


def validate_build_identity(
    value: object, side: str, specification: Mapping[str, Any],
    outer_executable: Mapping[str, Any], outer_archive: Mapping[str, Any],
) -> dict[str, Any]:
    expected_keys = {
        "archive_link_recipe", "archiver", "build_dir", "cmake_cache",
        "compile_commands", "compiler", "compiler_version_stdout",
        "executable_link_recipe", "validated_archive", "validated_archive_members",
        "validated_cache", "validated_compile_commands", "validated_executable",
        "strict_compile_recipes", "clean_recompile_proof", "strict_link_provenance",
    }
    require(isinstance(value, dict) and set(value) == expected_keys,
            f"retained {side} build identity is incomplete")
    build_dir = str(Path(specification[f"{side}_build_dir"]).resolve())
    require(value.get("build_dir") == build_dir,
            f"retained {side} build directory is invalid")
    for name, kind in (
        ("cmake_cache", "build_metadata"),
        ("compile_commands", "build_metadata"),
        ("executable_link_recipe", "build_metadata"),
        ("archive_link_recipe", "build_metadata"),
        ("compiler", "compiler"),
        ("archiver", "archiver"),
        ("validated_executable", "executable"),
        ("validated_archive", "archive"),
    ):
        validate_file_identity(value.get(name), f"{side} {name}", expected_kind=kind)
    require(value["validated_executable"] == outer_executable and
            value["validated_archive"] == outer_archive,
            f"retained {side} build artifacts differ from the top-level closure")
    for name in ("cmake_cache", "compile_commands", "executable_link_recipe",
                 "archive_link_recipe"):
        require(Path(value[name]["path"]).is_relative_to(Path(build_dir)),
                f"retained {side} {name} is outside its build directory")
    version = value.get("compiler_version_stdout")
    require(isinstance(version, dict) and set(version) == {"sha256", "text"} and
            isinstance(version.get("text"), str) and
            len(version["text"].encode("utf-8")) <= 65536 and
            version.get("sha256") == sha256_bytes(version["text"].encode("utf-8")),
            f"retained {side} compiler version is invalid")
    cache = value.get("validated_cache")
    expected_cache = {
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_CXX_COMPILER": "/usr/bin/c++",
        "CMAKE_CXX_FLAGS": "",
        "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
        "CMAKE_EXE_LINKER_FLAGS": "",
        "CMAKE_EXE_LINKER_FLAGS_RELEASE": "",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_GENERATOR": "Unix Makefiles",
        "ENABLE_OPENMP": "ON",
        "LEO2_BACKEND_VARIANT": "auto",
        "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_BUILD_TESTS": "OFF",
        "LEO2_ENABLE_CUDA": "OFF",
        "CMAKE_STATIC_LINKER_FLAGS": "",
        "CMAKE_STATIC_LINKER_FLAGS_RELEASE": "",
    }
    require(cache == expected_cache,
            f"retained {side} CMake production configuration is invalid")
    validate_file_identity(value.get("compiler"), f"{side} compiler",
                           expected_path=cache["CMAKE_CXX_COMPILER"],
                           expected_kind="compiler")
    validate_file_identity(value.get("archiver"), f"{side} archiver",
                           expected_path="/usr/bin/ar", expected_kind="archiver")
    members = value.get("validated_archive_members")
    require(isinstance(members, list) and members and len(members) == len(set(members)) and
            all(isinstance(item, str) and item and "/" not in item and "\\" not in item
                for item in members),
            f"retained {side} archive member list is invalid")
    commands = value.get("validated_compile_commands")
    require(isinstance(commands, dict) and set(commands) == {
        "entry_count", "isa_policy", "required_source_object_pairs",
        "required_sources", "validated_openmp", "validated_optimization"},
        f"retained {side} compile-command proof is incomplete")
    sources = commands.get("required_sources")
    pairs = commands.get("required_source_object_pairs")
    require(isinstance(sources, list) and sources == sorted(set(sources)) and
            all(isinstance(item, str) and Path(item).is_absolute() for item in sources) and
            isinstance(pairs, list) and len(pairs) == len(sources) and
            isinstance(commands.get("entry_count"), int) and
            not isinstance(commands["entry_count"], bool) and
            commands["entry_count"] >= len(sources) and
            commands.get("validated_openmp") is True and
            commands.get("validated_optimization") == "-O3" and
            commands.get("isa_policy") ==
            "portable core with ISA flags only on SSSE3/AVX2 translation units",
            f"retained {side} compile-command semantics are invalid")
    observed_sources: list[str] = []
    for pair in pairs:
        require(isinstance(pair, dict) and set(pair) == {"object", "source"},
                f"retained {side} source/object pair is invalid")
        source = validate_file_identity(pair.get("source"), f"{side} source",
                                        expected_kind="source_file")
        validate_file_identity(pair.get("object"), f"{side} object",
                               expected_kind="object_file")
        require(Path(source["path"]).is_relative_to(
            Path(specification[f"{side}_source_root"]).resolve()),
            f"retained {side} compiled source is outside its source tree")
        require(Path(pair["object"]["path"]).is_relative_to(Path(build_dir)),
                f"retained {side} object is outside its build tree")
        observed_sources.append(source["path"])
    require(sorted(observed_sources) == sources,
            f"retained {side} source/object closure differs from required sources")
    strict_recipes = value.get("strict_compile_recipes")
    require(isinstance(strict_recipes, dict) and
            set(strict_recipes) == {
                str(Path(path).relative_to(
                    Path(specification[f"{side}_source_root"]).resolve()))
                for path in sources
            } and
            all(isinstance(tokens, list) and tokens and
                all(isinstance(token, str) for token in tokens)
                for tokens in strict_recipes.values()),
            f"retained {side} strict compile recipes are invalid")
    proofs = value.get("clean_recompile_proof")
    require(isinstance(proofs, list) and len(proofs) == len(pairs),
            f"retained {side} clean-recompile proof is incomplete")
    pair_by_relative = {
        str(Path(pair["source"]["path"]).relative_to(
            Path(specification[f"{side}_source_root"]).resolve())): pair
        for pair in pairs
    }
    seen_proofs: set[str] = set()
    for proof in proofs:
        require(isinstance(proof, dict),
                f"retained {side} clean-recompile record is invalid")
        source_name = proof.get("source")
        require(isinstance(source_name, str) and source_name in pair_by_relative and
                source_name not in seen_proofs and
                proof.get("status") == "byte_identical_clean_recompile",
                f"retained {side} clean-recompile source/status is invalid")
        seen_proofs.add(source_name)
        pair = pair_by_relative[source_name]
        validate_clean_recompile_record(
            proof, source_name, pair, value["compiler"],
            strict_recipes[source_name])
    require(seen_proofs == set(pair_by_relative),
            f"retained {side} clean-recompile closure differs")
    link = value.get("strict_link_provenance")
    require(isinstance(link, dict) and set(link) == {
        "archive_member_content", "external_link_inputs",
        "normalized_executable_link", "ranlib", "relink_cache_key",
        "reproducible_executable_sha256"},
        f"retained {side} strict link provenance is incomplete")
    require(isinstance(link.get("normalized_executable_link"), list) and
            all(isinstance(token, str) for token in link["normalized_executable_link"]) and
            isinstance(link.get("external_link_inputs"), list) and
            isinstance(link.get("archive_member_content"), list) and
            isinstance(link.get("relink_cache_key"), str) and
            HEX256.fullmatch(link["relink_cache_key"]) is not None and
            link.get("reproducible_executable_sha256") ==
            value["validated_executable"]["sha256"],
            f"retained {side} strict link values are invalid")
    validate_file_identity(link.get("ranlib"), f"{side} ranlib",
                           expected_path="/usr/bin/ranlib", expected_kind="ranlib")
    object_by_path = {
        pair["object"]["path"]: pair["object"] for pair in pairs
    }
    require(len(link["archive_member_content"]) == len(members),
            f"retained {side} archive content proof is incomplete")
    for item, member in zip(link["archive_member_content"], members):
        require(isinstance(item, dict) and set(item) == {
            "member", "object_path", "sha256", "size"} and
            item.get("member") == member and item.get("object_path") in object_by_path and
            item.get("sha256") == object_by_path[item["object_path"]]["sha256"] and
            item.get("size") == object_by_path[item["object_path"]]["size"],
            f"retained {side} archive member content is invalid")
    for item in link["external_link_inputs"]:
        validate_file_identity(item, f"{side} external link input",
                               expected_kind="link_input")
    return value


def validate_runtime_closure(
    value: object, what: str, executable: str
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {"dependencies", "executable"} and
            value.get("executable") == str(Path(executable).resolve()) and
            isinstance(value.get("dependencies"), list) and value["dependencies"],
            f"retained {what} runtime closure is invalid")
    sonames: set[str] = set()
    for dependency in value["dependencies"]:
        require(isinstance(dependency, dict) and isinstance(dependency.get("soname"), str) and
                dependency["soname"] not in sonames,
                f"retained {what} runtime dependency is invalid or duplicated")
        sonames.add(dependency["soname"])
        if dependency.get("virtual") is True:
            require(set(dependency) == {"soname", "virtual"},
                    f"retained {what} virtual dependency is invalid")
        else:
            require(set(dependency) == {"file", "loader_path", "soname"} and
                    isinstance(dependency.get("loader_path"), str) and
                    Path(dependency["loader_path"]).is_absolute(),
                    f"retained {what} dependency path is invalid")
            validate_file_identity(
                dependency.get("file"), f"{what} dependency",
                expected_path=dependency["loader_path"])
    return value


def validate_retained_snapshot(
    value: object, specification: Mapping[str, Any]
) -> dict[str, Any]:
    expected_keys = {
        "candidate_archive", "candidate_build", "candidate_executable",
        "candidate_runtime_closure", "candidate_source", "control_archive",
        "control_build", "control_executable", "control_runtime_closure",
        "control_source", "git", "ldd", "python", "runner", "support", "taskset",
    }
    require(isinstance(value, dict) and set(value) == expected_keys,
            "retained source/build snapshot is incomplete")
    for name, spec_name, kind in (
        ("runner", "runner", "file"), ("taskset", "taskset", "executable"),
        ("ldd", "ldd", "executable"), ("git", "git", "executable"),
        ("python", "python", "executable"),
        ("control_executable", "control_executable", "executable"),
        ("candidate_executable", "candidate_executable", "executable"),
        ("control_archive", "control_archive", "archive"),
        ("candidate_archive", "candidate_archive", "archive"),
    ):
        validate_file_identity(
            value.get(name), name, specification[spec_name], kind)
    support = value.get("support")
    require(isinstance(support, dict) and set(support) ==
            {"artifact", "validated_contract"},
            "retained support identity is incomplete")
    expected_support = Path(specification["runner"]).parent.parent / \
        "main_compare" / "run_abba.py"
    support_artifact = validate_file_identity(
        support.get("artifact"), "support", str(expected_support), "file")
    contract = support.get("validated_contract")
    require(isinstance(contract, dict) and set(contract) ==
            {"contracts", "path", "sha256"} and
            contract.get("path") == support_artifact["path"] and
            contract.get("sha256") == support_artifact["sha256"] and
            contract.get("contracts") == validate_support_contract()["contracts"],
            "retained support contract is invalid")
    validate_source_identity(
        value.get("control_source"), "control",
        specification["control_source_root"], CONTROL_COMMIT, True)
    validate_source_identity(
        value.get("candidate_source"), "candidate",
        specification["candidate_source_root"], CANDIDATE_COMMIT, False)
    validate_build_identity(
        value.get("control_build"), "control", specification,
        value["control_executable"], value["control_archive"])
    validate_build_identity(
        value.get("candidate_build"), "candidate", specification,
        value["candidate_executable"], value["candidate_archive"])
    require(value["control_build"]["compiler"] ==
            value["candidate_build"]["compiler"] and
            value["control_build"]["compiler_version_stdout"] ==
            value["candidate_build"]["compiler_version_stdout"],
            "retained builds use different compilers")
    validate_runtime_closure(
        value.get("control_runtime_closure"), "control",
        specification["control_executable"])
    validate_runtime_closure(
        value.get("candidate_runtime_closure"), "candidate",
        specification["candidate_executable"])
    return value


class HardenedReservation:
    """Hold and continuously bind a coordinator-created CPU reservation."""

    def __init__(self, path: Path, cpu: int, sibling: int):
        absolute = path.absolute()
        resolved = path.resolve(strict=True)
        require(absolute == resolved,
                "CPU reservation path must be canonical and contain no symlink")
        self.path = resolved
        self.cpu = cpu
        self.sibling = sibling
        self.descriptor: int | None = None
        self.identity: dict[str, Any] | None = None
        self.kernel_socket: socket.socket | None = None
        material = canonical_bytes({
            "cpu": cpu,
            "path": str(self.path),
            "schema": RESERVATION_SCHEMA,
            "sibling": sibling,
            "uid": os.getuid(),
        })
        self.kernel_name = b"\0leopard2-reservation-v1-" + \
            hashlib.sha256(material).hexdigest()[:36].encode("ascii")

    def _acquire_kernel_lease(self) -> None:
        require(sys.platform.startswith("linux") and hasattr(socket, "AF_UNIX"),
                "Linux abstract Unix sockets are required for stable reservations")
        lease = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        lease.set_inheritable(False)
        try:
            lease.bind(self.kernel_name)
        except OSError as error:
            lease.close()
            if error.errno == errno.EADDRINUSE:
                raise EvidenceError("CPU reservation already has a kernel lease") from error
            raise EvidenceError(f"cannot bind stable CPU reservation: {error}") from error
        self.kernel_socket = lease

    def _release_kernel_lease(self) -> None:
        if self.kernel_socket is not None:
            self.kernel_socket.close()
            self.kernel_socket = None

    @staticmethod
    def parse(raw: bytes, cpu: int, sibling: int) -> dict[str, Any]:
        require(len(raw) <= 4096, "CPU reservation exceeds 4096 bytes")
        payload = parse_json_bytes(raw, "CPU reservation", 4096)
        expected = {
            "benchmark_cpu", "nonce", "owner", "reserved_sibling", "schema", "status",
        }
        require(isinstance(payload, dict) and set(payload) == expected,
                "CPU reservation has unexpected or missing fields")
        require(raw == canonical_bytes(payload),
                "CPU reservation must be canonical JSON without a trailing newline")
        require(payload.get("schema") == RESERVATION_SCHEMA and
                payload.get("status") == "held" and
                isinstance(payload.get("benchmark_cpu"), int) and
                not isinstance(payload.get("benchmark_cpu"), bool) and
                isinstance(payload.get("reserved_sibling"), int) and
                not isinstance(payload.get("reserved_sibling"), bool) and
                payload.get("benchmark_cpu") == cpu and
                payload.get("reserved_sibling") == sibling,
                "CPU reservation does not bind the requested pair")
        require(isinstance(payload.get("owner"), str) and
                1 <= len(payload["owner"].encode("utf-8")) <= 256,
                "CPU reservation owner is invalid")
        require(isinstance(payload.get("nonce"), str) and
                re.fullmatch(r"[0-9a-f]{32,128}", payload["nonce"]) is not None,
                "CPU reservation nonce is invalid")
        return payload

    def __enter__(self) -> dict[str, Any]:
        self._acquire_kernel_lease()
        try:
            parent = os.lstat(self.path.parent)
            path_before = os.lstat(self.path)
            require(stat.S_ISDIR(parent.st_mode) and parent.st_uid == os.getuid() and
                    stat.S_IMODE(parent.st_mode) == 0o700,
                    "CPU reservation parent must be an owned mode-0700 directory")
            require(stat.S_ISREG(path_before.st_mode) and
                    path_before.st_uid == os.getuid() and
                    path_before.st_nlink == 1 and
                    stat.S_IMODE(path_before.st_mode) == 0o600,
                    "CPU reservation must be an owned single-link mode-0600 file")
        except Exception:
            self._release_kernel_lease()
            raise
        flags = os.O_RDONLY | getattr(os, "O_NONBLOCK", 0)
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        try:
            self.descriptor = os.open(self.path, flags)
            metadata = os.fstat(self.descriptor)
            path_metadata = os.lstat(self.path)
            require(stat.S_ISREG(metadata.st_mode) and metadata.st_uid == os.getuid() and
                    metadata.st_nlink == 1 and stat.S_IMODE(metadata.st_mode) == 0o600 and
                    (metadata.st_dev, metadata.st_ino) ==
                    (path_before.st_dev, path_before.st_ino) ==
                    (path_metadata.st_dev, path_metadata.st_ino),
                    "CPU reservation must be an owned single-link mode-0600 file")
            try:
                fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise EvidenceError(
                    f"CPU reservation is already locked: {self.path}") from error
            raw = os.pread(self.descriptor, 4097, 0)
            payload = self.parse(raw, self.cpu, self.sibling)
            self.identity = {
                "device": metadata.st_dev,
                "directory_device": parent.st_dev,
                "directory_inode": parent.st_ino,
                "inode": metadata.st_ino,
                "lock": "exclusive_nonblocking_inode_bound",
                "mode": stat.S_IMODE(metadata.st_mode),
                "path": str(self.path),
                "payload": payload,
                "sha256": sha256_bytes(raw),
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

    def validate_current(self) -> None:
        require(self.descriptor is not None and self.identity is not None and
                self.kernel_socket is not None and
                self.kernel_socket.fileno() >= 0 and
                self.kernel_socket.getsockname() == self.kernel_name,
                "CPU reservation is not held")
        descriptor = os.fstat(self.descriptor)
        path = os.lstat(self.path)
        parent = os.lstat(self.path.parent)
        require(stat.S_ISREG(descriptor.st_mode) and descriptor.st_uid == os.getuid() and
                descriptor.st_nlink == 1 and stat.S_IMODE(descriptor.st_mode) == 0o600 and
                (descriptor.st_dev, descriptor.st_ino) == (path.st_dev, path.st_ino) ==
                (self.identity["device"], self.identity["inode"]),
                "CPU reservation path, inode, ownership, links, or mode changed")
        require(stat.S_ISDIR(parent.st_mode) and parent.st_uid == os.getuid() and
                stat.S_IMODE(parent.st_mode) == 0o700 and
                (parent.st_dev, parent.st_ino) ==
                (self.identity["directory_device"], self.identity["directory_inode"]),
                "CPU reservation parent changed")
        raw = os.pread(self.descriptor, 4097, 0)
        require(sha256_bytes(raw) == self.identity["sha256"] and
                self.parse(raw, self.cpu, self.sibling) == self.identity["payload"],
                "CPU reservation contents changed")

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


def validate_reservation_identity(
    value: object, cpu: int, sibling: int
) -> dict[str, Any]:
    expected = {
        "device", "directory_device", "directory_inode", "inode", "lock", "mode",
        "path", "payload", "sha256",
    }
    require(isinstance(value, dict) and set(value) == expected,
            "retained CPU reservation identity is incomplete")
    require(all(isinstance(value.get(name), int) and not isinstance(value[name], bool) and
                value[name] >= 0 for name in
                ("device", "directory_device", "directory_inode", "inode")),
            "retained CPU reservation filesystem identity is invalid")
    require(value.get("mode") == 0o600 and
            value.get("lock") == "exclusive_nonblocking_inode_bound" and
            isinstance(value.get("path"), str) and Path(value["path"]).is_absolute(),
            "retained CPU reservation mode, lock, or path is invalid")
    payload = value.get("payload")
    require(isinstance(payload, dict) and
            HardenedReservation.parse(canonical_bytes(payload), cpu, sibling) == payload,
            "retained CPU reservation payload is invalid")
    require(value.get("sha256") == sha256_bytes(canonical_bytes(payload)),
            "retained CPU reservation digest is invalid")
    return value


def benchmark_arguments(
    executable: Path, cell: Cell, campaign: Mapping[str, Any]
) -> list[str]:
    return [
        str(executable),
        "--k", str(cell.k),
        "--r", str(cell.r),
        "--profile", "low",
        "--field", cell.field,
        "--backend", cell.backend,
        "--bytes", str(cell.shard_bytes),
        "--loss", str(cell.losses),
        "--batch", "1",
        "--reuse", str(campaign["reuse"]),
        "--iterations", str(campaign["iterations"]),
        "--warmup", str(campaign["warmup"]),
        "--threads", "1",
        "--seed", str(cell.seed),
        "--skip-legacy",
        "--retain-samples",
        "--json", "-",
    ]


class BoundedChildError(EvidenceError):
    def __init__(
        self, message: str, returncode: int, stdout: bytes, stderr: bytes,
        duration_ns: int,
    ):
        super().__init__(message)
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.duration_ns = duration_ns


class ImmutableExecutables:
    """Copy attested executables to private unlinked read-only descriptors."""

    def __init__(
        self, stage: Path, specification: Mapping[str, Any],
        initial: Mapping[str, Any],
    ):
        self.stage = stage
        self.specification = specification
        self.initial = initial
        self.descriptors: dict[str, int] = {}
        self.identities: dict[str, dict[str, Any]] = {}

    def __enter__(self) -> tuple[dict[str, dict[str, Any]], dict[str, int]]:
        temporary = Path(tempfile.mkdtemp(prefix=".execution-", dir=self.stage))
        os.chmod(temporary, 0o700)
        try:
            for side in ("control", "candidate"):
                source_path = Path(self.specification[f"{side}_executable"])
                expected = self.initial[f"{side}_executable"]
                source_flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
                    getattr(os, "O_NOFOLLOW", 0)
                source = os.open(source_path, source_flags)
                try:
                    metadata = os.fstat(source)
                    path_metadata = os.lstat(source_path)
                    digest, size = sha256_descriptor(source, MAX_BUILD_ARTIFACT_BYTES)
                    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
                            (metadata.st_dev, metadata.st_ino) ==
                            (path_metadata.st_dev, path_metadata.st_ino) and
                            digest == expected["sha256"] and size == expected["size"] and
                            stat.S_IMODE(metadata.st_mode) == expected["mode"] and
                            os.access(source_path, os.X_OK),
                            f"{side} executable changed before immutable staging")
                    staged = temporary / side
                    descriptor = os.open(
                        staged, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                        getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
                        0o600)
                    try:
                        offset = 0
                        while offset < size:
                            block = os.pread(source, min(1024 * 1024, size - offset), offset)
                            require(block, f"{side} executable shortened during staging")
                            view = memoryview(block)
                            while view:
                                written = os.write(descriptor, view)
                                require(written > 0, f"short immutable {side} executable write")
                                view = view[written:]
                            offset += len(block)
                        os.fchmod(descriptor, 0o500)
                        os.fsync(descriptor)
                    finally:
                        os.close(descriptor)
                finally:
                    os.close(source)
                read_descriptor = os.open(
                    staged, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_NOFOLLOW", 0))
                staged_metadata = os.fstat(read_descriptor)
                staged_digest, staged_size = sha256_descriptor(
                    read_descriptor, MAX_BUILD_ARTIFACT_BYTES)
                require(staged_digest == expected["sha256"] and
                        staged_size == expected["size"] and
                        stat.S_IMODE(staged_metadata.st_mode) == 0o500,
                        f"private {side} executable copy differs from attestation")
                staged.unlink()
                self.descriptors[side] = read_descriptor
                self.identities[side] = {
                    "command_path": f"/proc/self/fd/{read_descriptor}",
                    "descriptor": read_descriptor,
                    "logical_path": expected["path"],
                    "mode": 0o500,
                    "schema": EXECUTION_SCHEMA,
                    "sha256": expected["sha256"],
                    "size": expected["size"],
                    "strategy": "private_unlinked_readonly_fd",
                }
            temporary.rmdir()
            return copy.deepcopy(self.identities), dict(self.descriptors)
        except Exception as primary:
            shutil.rmtree(temporary, ignore_errors=True)
            try:
                self.__exit__(None, None, None)
            except Exception as cleanup:
                raise EvidenceError(
                    "immutable executable staging failed and cleanup also failed: "
                    f"{type(primary).__name__}: {primary}; "
                    f"{type(cleanup).__name__}: {cleanup}") from primary
            raise

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        descriptors = list(self.descriptors.values())
        self.descriptors.clear()
        failures: list[str] = []
        for descriptor in descriptors:
            try:
                os.close(descriptor)
            except OSError as error:
                failures.append(f"fd {descriptor}: {error}")
        require(not failures,
                "cannot close immutable executable descriptors: " +
                "; ".join(failures))


def validate_execution_identity(
    value: object, initial: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    require(isinstance(value, dict) and set(value) == {"control", "candidate"},
            "immutable execution identity is incomplete")
    descriptors: set[int] = set()
    for side in ("control", "candidate"):
        item = value.get(side)
        expected = initial[f"{side}_executable"]
        require(isinstance(item, dict) and set(item) == {
            "command_path", "descriptor", "logical_path", "mode", "schema",
            "sha256", "size", "strategy"},
            f"immutable {side} execution identity is incomplete")
        descriptor = item.get("descriptor")
        require(item.get("schema") == EXECUTION_SCHEMA and
                item.get("strategy") == "private_unlinked_readonly_fd" and
                isinstance(descriptor, int) and not isinstance(descriptor, bool) and
                0 <= descriptor <= 1_048_575 and descriptor not in descriptors and
                item.get("command_path") == f"/proc/self/fd/{descriptor}" and
                item.get("logical_path") == expected["path"] and
                item.get("sha256") == expected["sha256"] and
                item.get("size") == expected["size"] and item.get("mode") == 0o500,
                f"immutable {side} execution identity differs from attestation")
        descriptors.add(descriptor)
    return value


def terminate_process_group_bounded(
    process: subprocess.Popen[bytes], timeout: float = 5.0
) -> tuple[bool, int]:
    """SIGKILL a child process group and make one bounded reap attempt."""
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


def run_bounded(
    command: Sequence[str], environment: Mapping[str, str], timeout: float,
    pass_fds: Sequence[int] = (),
) -> tuple[int, bytes, bytes, int]:
    require(isinstance(timeout, (int, float)) and not isinstance(timeout, bool) and
            math.isfinite(float(timeout)) and 0 < timeout <= MAX_TIMEOUT_SECONDS,
            "benchmark timeout is invalid")
    require(isinstance(command, Sequence) and 1 <= len(command) <= 64 and
            all(isinstance(item, str) and item and
                len(item.encode("utf-8")) <= MAX_PATH_BYTES for item in command),
            "benchmark command is invalid or oversized")
    require(len(pass_fds) <= 8 and all(
        isinstance(item, int) and not isinstance(item, bool) and 0 <= item <= 1_048_575
        for item in pass_fds), "benchmark inherited descriptor list is invalid")
    started = time.monotonic_ns()
    process = subprocess.Popen(
        list(command), cwd="/", env=dict(environment), stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True,
        pass_fds=tuple(pass_fds))
    require(process.stdout is not None and process.stderr is not None,
            "cannot capture benchmark pipes")
    selector = selectors.DefaultSelector()
    stdout_fd = process.stdout.fileno()
    stderr_fd = process.stderr.fileno()
    outputs: dict[int, bytearray] = {stdout_fd: bytearray(), stderr_fd: bytearray()}
    limits = {
        process.stdout.fileno(): MAX_STDOUT_BYTES,
        process.stderr.fileno(): MAX_STDERR_BYTES,
    }
    for stream in (process.stdout, process.stderr):
        os.set_blocking(stream.fileno(), False)
        selector.register(stream, selectors.EVENT_READ)
    deadline = time.monotonic() + timeout
    failure: str | None = None
    try:
        while selector.get_map():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                failure = f"benchmark exceeded {timeout:.3f} seconds"
                break
            events = selector.select(min(remaining, 0.1))
            for key, _ in events:
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
                    failure = "benchmark output exceeded its retained byte limit"
                    break
            if failure is not None:
                break
        if failure is not None:
            reaped, returncode = terminate_process_group_bounded(process)
            if not reaped:
                failure += "; process group was not reapable within five seconds"
        else:
            try:
                returncode = process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                reaped, returncode = terminate_process_group_bounded(process)
                failure = (
                    "benchmark closed its output but did not terminate within five seconds"
                )
                if not reaped:
                    failure += "; process group was not reapable within five seconds"
    except subprocess.TimeoutExpired:
        reaped, returncode = terminate_process_group_bounded(process)
        failure = "benchmark process did not terminate within five seconds after SIGKILL"
        if not reaped:
            failure += "; bounded reap failed"
    except BaseException:
        terminate_process_group_bounded(process)
        raise
    finally:
        selector.close()
        process.stdout.close()
        process.stderr.close()
    duration = time.monotonic_ns() - started
    stdout = bytes(outputs[stdout_fd])
    stderr = bytes(outputs[stderr_fd])
    if failure is not None:
        raise BoundedChildError(
            failure, returncode,
            stdout[:MAX_STDOUT_BYTES], stderr[:MAX_STDERR_BYTES], duration)
    return returncode, stdout, stderr, duration


class XorShift64:
    def __init__(self, seed: int):
        self.state = seed & ((1 << 64) - 1)
        if self.state == 0:
            self.state = 0x9E3779B97F4A7C15

    def next(self) -> int:
        value = self.state
        value ^= (value << 13) & ((1 << 64) - 1)
        value ^= value >> 7
        value ^= (value << 17) & ((1 << 64) - 1)
        self.state = value & ((1 << 64) - 1)
        return self.state


def expected_losses(cell: Cell) -> list[int]:
    order = list(range(cell.k))
    random = XorShift64(cell.seed ^ 0xD1B54A32D192ED03)
    for remaining in range(cell.k, 1, -1):
        selected = random.next() % remaining
        order[remaining - 1], order[selected] = order[selected], order[remaining - 1]
    return sorted(order[:cell.losses])


def validate_digests(value: object) -> dict[str, str]:
    require(isinstance(value, dict) and set(value) == {
        "algorithm", "original_data", "recovered_originals", "transmitted_parity"} and
            value.get("algorithm") == "fnv1a64",
            "workload digests are missing or use another algorithm")
    result: dict[str, str] = {}
    for name in ("original_data", "transmitted_parity", "recovered_originals"):
        digest = value.get(name)
        require(isinstance(digest, str) and HEX64.fullmatch(digest) is not None,
                f"invalid workload digest {name}")
        result[name] = digest
    # benchmark-v2 deliberately uses different domains here: original_data
    # covers every K source shard, while recovered_originals covers only the L
    # missing shards.  Round-trip truth is asserted by the benchmark's checked
    # decoder, and validate_invocations independently requires all three
    # digests (including parity and repaired output) to match across A/B runs.
    return result


def validate_result(
    value: object, cell: Cell, campaign: Mapping[str, Any]
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "build", "correctness", "legacy", "memory", "metrics", "parameters",
        "resolved", "schema", "workload_digests"} and
            value.get("schema") == "leopard2-benchmark-v2",
            "benchmark returned the wrong JSON schema")
    build = value.get("build")
    require(isinstance(build, dict) and set(build) ==
            {"compiler", "compiler_version", "cplusplus"} and
            isinstance(build.get("compiler"), str) and build["compiler"] and
            isinstance(build.get("compiler_version"), str) and
            len(build["compiler_version"].encode("utf-8")) <= 4096 and
            isinstance(build.get("cplusplus"), int) and
            not isinstance(build["cplusplus"], bool) and build["cplusplus"] >= 201103,
            "benchmark build identity is invalid")
    parameters = value.get("parameters")
    require(isinstance(parameters, dict), "benchmark parameters are missing")
    expected = {
        "K": cell.k,
        "R": cell.r,
        "requested_profile": "low_v1",
        "requested_field": cell.field,
        "requested_backend": cell.backend,
        "force_generic_decode": False,
        "force_specialized_decode": False,
        "skip_legacy": True,
        "retain_samples": True,
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
    require(set(parameters) == set(expected),
            "benchmark parameters have unexpected or missing fields")
    for name, expected_value in expected.items():
        require(exact_json_equal(parameters.get(name), expected_value),
                f"benchmark parameter {name} differs: "
                f"{parameters.get(name)!r} != {expected_value!r}")
    resolved = value.get("resolved")
    padded = ceil_power_of_two(cell.k)
    parent = ceil_power_of_two(padded + cell.r)
    require(isinstance(resolved, dict) and set(resolved) == {
            "backend", "field", "padded_side", "parent_count", "profile", "thread_count"} and
            resolved.get("profile") == "low_v1" and
            resolved.get("field") == cell.field and
            resolved.get("backend") == cell.backend and
            type(resolved.get("thread_count")) is int and
            resolved.get("thread_count") == 1 and
            type(resolved.get("parent_count")) is int and
            resolved.get("parent_count") == parent and
            type(resolved.get("padded_side")) is int and
            resolved.get("padded_side") == padded,
            "benchmark resolved a different profile, field, backend, or parent")
    correctness = value.get("correctness")
    require(isinstance(correctness, dict) and set(correctness) ==
            {"legacy_comparison", "leopard2_round_trip"} and
            correctness.get("leopard2_round_trip") is True and
            correctness.get("legacy_comparison") is None,
            "benchmark round trip failed")
    memory = value.get("memory")
    require(isinstance(memory, dict) and set(memory) == {
        "decode_scratch_bytes_batch", "decode_scratch_bytes_per_stripe",
        "encode_scratch_bytes_batch", "encode_scratch_bytes_per_stripe",
        "scratch_alignment"} and
        all(isinstance(memory.get(name), int) and not isinstance(memory[name], bool) and
            memory[name] >= 0 for name in memory) and
        memory["scratch_alignment"] > 0 and
        memory["scratch_alignment"] & (memory["scratch_alignment"] - 1) == 0 and
        memory["encode_scratch_bytes_batch"] ==
        memory["encode_scratch_bytes_per_stripe"] and
        memory["decode_scratch_bytes_batch"] ==
        memory["decode_scratch_bytes_per_stripe"],
        "benchmark scratch identity is invalid")
    legacy = value.get("legacy")
    require(isinstance(legacy, dict) and set(legacy) == {
        "available", "codec_setup", "decode_including_setup",
        "decode_timing_includes_setup", "encode_execution", "unavailable_reason"} and
            legacy.get("available") is False and
            legacy.get("unavailable_reason") == "disabled by --skip-legacy",
            "benchmark silently ran legacy comparison")
    require(legacy.get("codec_setup") is None and
            legacy.get("decode_including_setup") is None and
            legacy.get("encode_execution") is None and
            legacy.get("decode_timing_includes_setup") is True,
            "disabled legacy metrics are not null with explicit setup semantics")
    metrics = value.get("metrics")
    require(isinstance(metrics, dict) and set(metrics) == {
        "codec_setup", "decode_amortized_at_reuse", "decode_execution",
        "decode_plan_setup", "encode_execution", "rate_semantics"},
        "benchmark metrics are missing or contain unexpected fields")
    require(metrics.get("rate_semantics") ==
            "offered_received counts all non-null shard pointers supplied; "
            "a plan may read a deterministic subset",
            "benchmark rate semantics changed")
    iterations = campaign["iterations"]
    setup = validate_summary(metrics.get("codec_setup"), iterations, True)
    encode = validate_execution_summary(
        metrics.get("encode_execution"), iterations,
        "input_GB_per_s", "parity_output_GB_per_s",
        cell.k * cell.shard_bytes, cell.r * cell.shard_bytes)
    plan = validate_summary(metrics.get("decode_plan_setup"), iterations, True)
    decode = validate_execution_summary(
        metrics.get("decode_execution"), iterations,
        "offered_received_GB_per_s", "repaired_output_GB_per_s",
        (cell.k - cell.losses + cell.r) * cell.shard_bytes,
        cell.losses * cell.shard_bytes)
    amortized = metrics.get("decode_amortized_at_reuse")
    require(isinstance(amortized, dict) and set(amortized) == {
        "derived_median_us_per_batch_call", "offered_received_GB_per_s",
        "repaired_output_GB_per_s", "reuse_count"} and
            amortized.get("reuse_count") == campaign["reuse"],
            "benchmark amortized decode reuse differs")
    expected_amortized = statistics.median(decode) + \
        statistics.median(plan) / campaign["reuse"]
    require(close_enough(finite_number(
        amortized.get("derived_median_us_per_batch_call"), "amortized decode"),
        expected_amortized),
        "amortized decode is not derived from retained setup/execution samples")
    for name, byte_count in (
        ("offered_received_GB_per_s", (cell.k - cell.losses + cell.r) * cell.shard_bytes),
        ("repaired_output_GB_per_s", cell.losses * cell.shard_bytes),
    ):
        require(close_enough(
            finite_number(amortized.get(name), f"amortized {name}"),
            byte_count / (expected_amortized * 1000.0)),
            f"amortized {name} is not derived from bytes and timing")
    return {
        "backend": cell.backend,
        "codec_setup_samples_us": setup,
        "decode_execution_samples_us_per_batch_call": decode,
        "decode_plan_setup_samples_us": plan,
        "digests": validate_digests(value.get("workload_digests")),
        "encode_execution_samples_us_per_batch_call": encode,
    }


def confidence_interval(round_log_contrasts: Sequence[float]) -> dict[str, Any]:
    require(len(round_log_contrasts) == ROUNDS,
            f"analysis requires {ROUNDS} independent round contrasts")
    mean = statistics.fmean(round_log_contrasts)
    deviation = statistics.stdev(round_log_contrasts)
    critical = 2.7764451051977987  # two-sided 95%, Student-t df=4
    half_width = critical * deviation / math.sqrt(ROUNDS)
    return {
        "ratio_orientation": "control_time_over_candidate_time",
        "independent_round_count": ROUNDS,
        "degrees_of_freedom": ROUNDS - 1,
        "independent_round_log_contrasts": list(round_log_contrasts),
        "geometric_speedup": math.exp(mean),
        "ci95_lower": math.exp(mean - half_width),
        "ci95_upper": math.exp(mean + half_width),
        "constituent_pair_count": 2 * ROUNDS,
    }


def normalized_median(value: Mapping[str, Any], metric: str) -> float:
    name = (
        "codec_setup_samples_us" if metric == "codec_setup"
        else "encode_execution_samples_us_per_batch_call"
    )
    samples = value[name]
    require(isinstance(samples, list) and samples,
            f"normalized {metric} samples are missing")
    return statistics.median(samples)


def analyze_cell(records: Sequence[Mapping[str, Any]], cell: Cell) -> dict[str, Any]:
    require(len(records) == ROUNDS * len(SEQUENCE),
            f"cell {cell.identifier} has the wrong invocation count")
    observations = {"codec_setup": [], "encode_execution": []}
    within = {"codec_setup": [], "encode_execution": []}
    for round_index in range(ROUNDS):
        group = records[round_index * 4:(round_index + 1) * 4]
        require(tuple((item.get("label"), item.get("implementation")) for item in group) ==
                SEQUENCE,
                f"cell {cell.identifier} round {round_index} is not A1/B1/B2/A2")
        for metric in observations:
            contrasts = []
            for control_index, candidate_index in ((0, 1), (3, 2)):
                control = normalized_median(group[control_index]["normalized"], metric)
                candidate = normalized_median(
                    group[candidate_index]["normalized"], metric)
                contrasts.append(math.log(control) - math.log(candidate))
            within[metric].append(contrasts)
            observations[metric].append(statistics.fmean(contrasts))
    result = {
        metric: confidence_interval(values)
        for metric, values in observations.items()
    }
    for metric in result:
        result[metric]["within_round_log_contrasts"] = within[metric]
    threshold = TARGET_SPEEDUP if cell.role == "target" else NEIGHBOR_REGRESSION_FLOOR
    result["encode_gate"] = {
        "cell_role": cell.role,
        "required_ci95_lower": threshold,
        "passed": result["encode_execution"]["ci95_lower"] >= threshold,
    }
    return result


def analyze(
    invocations: Sequence[Mapping[str, Any]], cells: Sequence[Cell]
) -> dict[str, Any]:
    by_cell: dict[str, list[Mapping[str, Any]]] = {
        cell.identifier: [] for cell in cells
    }
    for invocation in invocations:
        identifier = invocation.get("cell_id")
        require(identifier in by_cell, f"invocation uses unknown cell {identifier!r}")
        by_cell[identifier].append(invocation)
    cells_analysis = {
        cell.identifier: analyze_cell(by_cell[cell.identifier], cell)
        for cell in cells
    }
    failures = sorted(
        identifier for identifier, result in cells_analysis.items()
        if not result["encode_gate"]["passed"]
    )
    return {
        "cells": cells_analysis,
        "policy": {
            "failed_cells": failures,
            "passed": not failures,
            "policy": statistics_policy()["promotion"],
        },
    }


def run_child(
    implementation: str,
    label: str,
    cell: Cell,
    round_index: int,
    slot: int,
    campaign: Mapping[str, Any],
    specification: Mapping[str, Any],
    initial_identity: Mapping[str, Any],
    reservation: Mapping[str, Any],
    execution: Mapping[str, Mapping[str, Any]],
    execution_descriptors: Mapping[str, int],
    output: Path,
    cpu: int,
    timeout: float,
    snapshot: Callable[[], Mapping[str, Any]],
    validate_guards: Callable[[], None],
) -> dict[str, Any]:
    validate_guards()
    identity_before = snapshot()
    require(identity_before == initial_identity,
            "input identity changed before benchmark invocation")
    executable = Path(execution[implementation]["command_path"])
    command = [
        specification["taskset"], "-c", str(cpu),
        *benchmark_arguments(executable, cell, campaign),
    ]
    stem = f"{cell.identifier}/r{round_index:02d}-{slot:02d}-{label}-{implementation}"
    stdout_path = output / "children" / f"{stem}.stdout.json"
    stderr_path = output / "children" / f"{stem}.stderr.txt"
    children_root = output / "children"
    children_root.mkdir(mode=0o700, exist_ok=True)
    os.chmod(children_root, 0o700)
    stdout_path.parent.mkdir(mode=0o700, exist_ok=True)
    os.chmod(stdout_path.parent, 0o700)
    try:
        returncode, stdout, stderr, duration_ns = run_bounded(
            command, CHILD_ENVIRONMENT, timeout,
            (execution_descriptors[implementation],))
    except BoundedChildError as error:
        write_bytes_exclusive(stdout_path, error.stdout)
        write_bytes_exclusive(stderr_path, error.stderr)
        raise EvidenceError(
            f"benchmark {implementation}/{cell.identifier}/{label}: {error}") from error
    write_bytes_exclusive(stdout_path, stdout)
    write_bytes_exclusive(stderr_path, stderr)
    stdout_record = artifact_record(output, stdout_path)
    stderr_record = artifact_record(output, stderr_path)
    require(returncode == 0,
            f"benchmark {implementation}/{cell.identifier}/{label} exited {returncode}")
    result = parse_json_bytes(stdout, "benchmark stdout", MAX_STDOUT_BYTES)
    normalized = validate_result(result, cell, campaign)
    validate_guards()
    identity_after = snapshot()
    require(identity_after == initial_identity,
            "input identity changed after benchmark invocation")
    return {
        "cell_id": cell.identifier,
        "command": command,
        "duration_ns": duration_ns,
        "environment": dict(CHILD_ENVIRONMENT),
        "identity_after_digest": sha256_bytes(canonical_bytes(identity_after)),
        "identity_before_digest": sha256_bytes(canonical_bytes(identity_before)),
        "implementation": implementation,
        "label": label,
        "normalized": normalized,
        "pinned_cpu": cpu,
        "reservation_after": reservation,
        "reservation_before": reservation,
        "result": result,
        "returncode": returncode,
        "round": round_index,
        "slot": slot,
        "stderr": stderr_record,
        "stdout": stdout_record,
    }


def campaign_cells(
    campaign: Mapping[str, Any], allow_self_test: bool
) -> tuple[Cell, ...]:
    values = campaign.get("cells")
    require(isinstance(values, list) and values, "campaign has no cells")
    cells: list[Cell] = []
    for value in values:
        require(isinstance(value, dict) and set(value) == set(asdict(FIXED_CELLS[0])),
                "campaign cell has unexpected or missing fields")
        try:
            cell = Cell(**value)
        except (TypeError, ValueError) as error:
            raise EvidenceError(f"campaign cell is invalid: {error}") from error
        validate_cell(cell)
        cells.append(cell)
    require(len({cell.identifier for cell in cells}) == len(cells),
            "campaign cell identifiers are not unique")
    if allow_self_test:
        require((campaign.get("mode") == "self_test" and len(cells) == 1) or
                (campaign.get("mode") == "self_test_full" and
                 tuple(cells) == FIXED_CELLS),
                "self-test campaign identity is invalid")
    else:
        require(campaign.get("mode") == "authoritative" and
                tuple(cells) == FIXED_CELLS,
                "authoritative campaign differs from the fixed matrix")
    return tuple(cells)


def validate_campaign(
    value: object, allow_self_test: bool = False
) -> tuple[dict[str, Any], tuple[Cell, ...]]:
    require(isinstance(value, dict), "campaign is not an object")
    expected_keys = {
        "allowed_cpu_set_at_launch", "batch", "benchmark_cpu", "cells",
        "child_environment", "iterations", "mode", "reserved_sibling", "reuse",
        "rounds", "sequence", "statistics", "threads", "timeout_seconds", "warmup",
    }
    require(set(value) == expected_keys, "campaign has unexpected or missing fields")
    require(isinstance(value.get("rounds"), int) and
            not isinstance(value.get("rounds"), bool) and
            value.get("rounds") == ROUNDS and
            value.get("sequence") == [list(item) for item in SEQUENCE],
            "campaign does not use five A1/B1/B2/A2 rounds")
    require(type(value.get("batch")) is int and value.get("batch") == 1 and
            type(value.get("threads")) is int and value.get("threads") == 1 and
            value.get("child_environment") == CHILD_ENVIRONMENT,
            "campaign is not a strict one-stripe, one-thread run")
    for name, minimum, maximum in (
        ("iterations", 5, MAX_ITERATIONS), ("reuse", 1, MAX_REUSE),
        ("warmup", 1, MAX_WARMUP),
    ):
        item = value.get(name)
        require(isinstance(item, int) and not isinstance(item, bool) and
                minimum <= item <= maximum,
                f"campaign {name} is invalid")
    timeout = value.get("timeout_seconds")
    require(isinstance(timeout, (int, float)) and not isinstance(timeout, bool) and
            0 < timeout <= MAX_TIMEOUT_SECONDS and
            math.isfinite(float(timeout)),
            "campaign timeout is invalid")
    require(value.get("statistics") == statistics_policy(),
            "campaign statistics or promotion policy was edited")
    cpu = value.get("benchmark_cpu")
    sibling = value.get("reserved_sibling")
    allowed = value.get("allowed_cpu_set_at_launch")
    require(isinstance(cpu, int) and not isinstance(cpu, bool) and cpu >= 0 and
            isinstance(sibling, int) and not isinstance(sibling, bool) and sibling >= 0 and
            cpu != sibling and isinstance(allowed, list) and
            allowed == sorted(set(allowed)) and
            all(isinstance(item, int) and not isinstance(item, bool) and item >= 0
                and item <= SUPPORT.MAX_CPU_ID for item in allowed) and
            len(allowed) <= SUPPORT.MAX_CPU_LIST_ENTRIES and
            cpu in allowed and sibling in allowed and
            bool(set(allowed) - {cpu, sibling}),
            "campaign CPU pair, allowed set, or housekeeping set is invalid")
    return value, campaign_cells(value, allow_self_test)


def validate_input_specification(value: object) -> dict[str, str]:
    expected = {
        "candidate_archive", "candidate_build_dir", "candidate_executable",
        "candidate_source_root", "control_archive", "control_build_dir",
        "control_executable", "control_source_root", "git", "ldd", "python",
        "runner", "taskset",
    }
    require(isinstance(value, dict) and set(value) == expected and
            all(isinstance(item, str) and item for item in value.values()),
            "input specification is incomplete or has unexpected fields")
    for name, item in value.items():
        validate_path_text(item, f"input specification {name}", absolute=True)
    return value


def expected_invocation_sequence(cells: Sequence[Cell]) -> list[tuple[Any, ...]]:
    return [
        (cell.identifier, round_index, slot, label, implementation)
        for cell in cells
        for round_index in range(ROUNDS)
        for slot, (label, implementation) in enumerate(SEQUENCE)
    ]


def validate_invocations(
    invocations: object,
    cells: Sequence[Cell],
    campaign: Mapping[str, Any],
    specification: Mapping[str, Any],
    initial: Mapping[str, Any],
    reservation: Mapping[str, Any],
    execution: Mapping[str, Mapping[str, Any]],
    output: Path,
    complete: bool,
) -> list[dict[str, Any]]:
    require(isinstance(invocations, list), "invocations are not a list")
    expected = expected_invocation_sequence(cells)
    require(len(invocations) == len(expected) if complete else len(invocations) <= len(expected),
            "invocation count is invalid")
    cell_by_id = {cell.identifier: cell for cell in cells}
    digest_by_cell: dict[str, dict[str, str]] = {}
    validated: list[dict[str, Any]] = []
    referenced_artifacts: list[dict[str, Any]] = []
    cpu = campaign["benchmark_cpu"]
    for invocation, sequence in zip(invocations, expected):
        require(isinstance(invocation, dict) and set(invocation) == {
            "cell_id", "command", "duration_ns", "environment",
            "identity_after_digest", "identity_before_digest", "implementation",
            "label", "normalized", "pinned_cpu", "reservation_after",
            "reservation_before", "result", "returncode", "round", "slot",
            "stderr", "stdout"},
            "invocation is not an exact evidence object")
        observed = (
            invocation.get("cell_id"), invocation.get("round"), invocation.get("slot"),
            invocation.get("label"), invocation.get("implementation"),
        )
        require(type(invocation.get("round")) is int and
                type(invocation.get("slot")) is int and
                observed == sequence,
                f"invocation order/relabel mismatch: {observed!r} != {sequence!r}")
        cell = cell_by_id[sequence[0]]
        implementation = sequence[4]
        expected_command = [
            specification["taskset"], "-c", str(cpu),
            *benchmark_arguments(
                Path(execution[implementation]["command_path"]), cell, campaign),
        ]
        require(invocation.get("command") == expected_command and
                invocation.get("environment") == CHILD_ENVIRONMENT and
                type(invocation.get("pinned_cpu")) is int and
                invocation.get("pinned_cpu") == cpu and
                type(invocation.get("returncode")) is int and
                invocation.get("returncode") == 0,
                "invocation execution contract was edited")
        duration = invocation.get("duration_ns")
        require(isinstance(duration, int) and not isinstance(duration, bool) and
                0 < duration <= int((campaign["timeout_seconds"] + 6.0) * 1e9),
                "invocation duration is invalid")
        initial_digest = sha256_bytes(canonical_bytes(initial))
        require(invocation.get("identity_before_digest") == initial_digest and
                invocation.get("identity_after_digest") == initial_digest and
                invocation.get("reservation_before") == reservation and
                invocation.get("reservation_after") == reservation,
                "invocation source/build/reservation identity was edited")
        stdout_path = verify_artifact(
            output, invocation.get("stdout"), "stdout", MAX_STDOUT_BYTES)
        verify_artifact(
            output, invocation.get("stderr"), "stderr", MAX_STDERR_BYTES)
        referenced_artifacts.extend((invocation["stdout"], invocation["stderr"]))
        parsed = read_json_limited(stdout_path, "retained benchmark stdout", MAX_STDOUT_BYTES)
        require(parsed == invocation.get("result"),
                "retained stdout differs from embedded benchmark result")
        normalized = validate_result(parsed, cell, campaign)
        require(invocation.get("normalized") == normalized,
                "normalized benchmark result was edited")
        digests = normalized["digests"]
        if cell.identifier in digest_by_cell:
            require(digests == digest_by_cell[cell.identifier],
                    f"control/candidate wire or recovery digests differ in {cell.identifier}")
        else:
            digest_by_cell[cell.identifier] = digests
        validated.append(invocation)
    if complete:
        require(
            sorted(referenced_artifacts, key=lambda item: item["path"]) ==
            retained_child_records(output),
            "successful output contains unbound or missing child artifacts")
    return validated


def validate_raw(
    raw: object,
    output: Path,
    check_current_inputs: bool,
    allow_self_test: bool = False,
) -> dict[str, Any]:
    raw = verify_signature(raw, "raw bundle")
    expected_keys = {
        "analysis", "campaign", "created_utc", "digest", "host_final", "host_initial",
        "identities_final", "identities_initial", "input_specification", "invocations",
        "execution",
        "isolation", "reservation", "schema", "validity_is_independent_of_speed",
    }
    require(set(raw) == expected_keys and raw.get("schema") == RAW_SCHEMA and
            raw.get("validity_is_independent_of_speed") is True,
            "raw bundle schema or keys are invalid")
    validate_utc(raw.get("created_utc"), "raw created_utc")
    campaign, cells = validate_campaign(raw.get("campaign"), allow_self_test)
    cpu = campaign["benchmark_cpu"]
    sibling = campaign["reserved_sibling"]
    allowed = campaign["allowed_cpu_set_at_launch"]
    require(raw.get("host_initial") == raw.get("host_final"),
            "host policy/topology changed during campaign")
    SUPPORT.validate_host_record(raw.get("host_initial"), cpu, sibling, allowed)
    isolation = SUPPORT.validate_isolation(raw.get("isolation"), cpu, sibling)
    reservation = validate_reservation_identity(raw.get("reservation"), cpu, sibling)
    specification = validate_input_specification(raw.get("input_specification"))
    initial = raw.get("identities_initial")
    require(isinstance(initial, dict) and raw.get("identities_final") == initial,
            "input source/build identities changed during campaign")
    if not allow_self_test:
        validate_retained_snapshot(initial, specification)
    if check_current_inputs:
        require(input_snapshot(specification) == initial,
                "current source/build/artifact closure differs from retained evidence")
    execution = validate_execution_identity(raw.get("execution"), initial)
    invocations = validate_invocations(
        raw.get("invocations"), cells, campaign, specification, initial,
        reservation, execution, output, complete=True)
    total_duration = sum(item["duration_ns"] for item in invocations)
    elapsed = isolation["after"]["monotonic_ns"] - isolation["before"]["monotonic_ns"]
    require(elapsed >= total_duration,
            "CPU isolation interval does not cover all child durations")
    calculated = analyze(invocations, cells)
    require(raw.get("analysis") == calculated, "raw paired analysis was edited")
    return calculated


def retained_child_records(output: Path) -> list[dict[str, Any]]:
    root = output / "children"
    if not root.exists():
        return []
    root_metadata = os.lstat(root)
    require(stat.S_ISDIR(root_metadata.st_mode) and
            root_metadata.st_uid == os.getuid() and
            stat.S_IMODE(root_metadata.st_mode) == 0o700,
            "children artifact root is not an owned mode-0700 directory")
    records: list[dict[str, Any]] = []
    directories: set[str] = set()
    stack = [root]
    entries = 0
    while stack:
        directory = stack.pop()
        with os.scandir(directory) as iterator:
            children = sorted(iterator, key=lambda entry: entry.name)
        for entry in children:
            entries += 1
            require(entries <= MAX_INVENTORY_ENTRIES,
                    "child artifact inventory contains too many entries")
            path = Path(entry.path)
            metadata = entry.stat(follow_symlinks=False)
            relative = str(path.relative_to(root))
            if stat.S_ISDIR(metadata.st_mode):
                require(metadata.st_uid == os.getuid() and
                        stat.S_IMODE(metadata.st_mode) == 0o700,
                        f"child artifact directory is not owned mode-0700: {path}")
                directories.add(relative)
                stack.append(path)
                continue
            require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
                    f"unexpected retained artifact type: {path}")
            if path.name.endswith(".stdout.json"):
                limit = MAX_STDOUT_BYTES
            elif path.name.endswith(".stderr.txt"):
                limit = MAX_STDERR_BYTES
            else:
                raise EvidenceError(f"unexpected retained artifact name: {path}")
            records.append(artifact_record(output, path, limit))
    expected_directories: set[str] = set()
    for record in records:
        relative = Path(record["path"]).relative_to("children")
        parent = relative.parent
        while parent != Path("."):
            expected_directories.add(str(parent))
            parent = parent.parent
    require(directories == expected_directories,
            "child artifact tree contains an empty or unbound directory")
    return sorted(records, key=lambda item: item["path"])


def validate_failure_lifecycle(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "completed_phases", "failed_phase", "teardown_error", "teardown_status"},
        "failure lifecycle is incomplete")
    completed = value.get("completed_phases")
    require(isinstance(completed, list) and completed and
            all(isinstance(item, str) for item in completed) and
            completed == list(FAILURE_PHASES[:len(completed)]),
            "failure lifecycle phases are not an exact prefix")
    failed_phase = value.get("failed_phase")
    expected_next = (
        FAILURE_PHASES[len(completed)] if len(completed) < len(FAILURE_PHASES)
        else "teardown"
    )
    teardown_status = value.get("teardown_status")
    teardown_error = value.get("teardown_error")
    require(teardown_status in {"completed", "failed"} and
            ((teardown_status == "completed" and teardown_error is None) or
             (teardown_status == "failed" and isinstance(teardown_error, str) and
              1 <= len(teardown_error.encode("utf-8")) <= 8192)),
            "failure lifecycle teardown result is invalid")
    require(
        (teardown_status == "failed" and
         (failed_phase == expected_next or
          (failed_phase == "teardown" and completed == list(FAILURE_PHASES)))) or
        (teardown_status == "completed" and failed_phase == expected_next) or
        (teardown_status == "completed" and failed_phase == "publication" and
         completed == list(FAILURE_PHASES)),
        "failure lifecycle failed phase is inconsistent")
    return value


def validate_failure(
    failure: object, output: Path, allow_self_test: bool = False
) -> dict[str, Any]:
    failure = verify_signature(failure, "failure bundle")
    expected_keys = {
        "campaign", "created_utc", "digest", "error", "error_type", "host_initial",
        "identities_initial", "input_specification", "invocations", "isolation",
        "execution",
        "lifecycle", "pair_lease", "reservation", "retained_files", "schema",
        "status", "traceback",
        "valid",
    }
    require(set(failure) == expected_keys and failure.get("schema") == FAILURE_SCHEMA and
            failure.get("status") == "failed" and failure.get("valid") is False,
            "failure bundle schema, status, or keys are invalid")
    validate_utc(failure.get("created_utc"), "failure created_utc")
    for name, maximum in (("error_type", 256), ("error", 8192), ("traceback", 65536)):
        item = failure.get(name)
        require(isinstance(item, str) and 1 <= len(item.encode("utf-8")) <= maximum,
                f"failure {name} is invalid or oversized")
    lifecycle = validate_failure_lifecycle(failure.get("lifecycle"))
    completed = set(lifecycle["completed_phases"])
    campaign, cells = validate_campaign(failure.get("campaign"), allow_self_test)
    cpu = campaign["benchmark_cpu"]
    sibling = campaign["reserved_sibling"]
    specification = validate_input_specification(failure.get("input_specification"))
    host = failure.get("host_initial")
    if "topology_validated" in completed:
        require(host is not None, "failure topology phase lacks host identity")
    if host is not None:
        SUPPORT.validate_host_record(
            host, cpu, sibling, campaign["allowed_cpu_set_at_launch"])
    pair_lease = failure.get("pair_lease")
    if "locks_acquired" in completed:
        require(pair_lease is not None, "failure lock phase lacks pair lease")
    if pair_lease is not None:
        SUPPORT.validate_pair_lease_identity(pair_lease, cpu, sibling)
    reservation = failure.get("reservation")
    if "locks_acquired" in completed:
        require(reservation is not None, "failure lock phase lacks reservation")
    if reservation is not None:
        validate_reservation_identity(reservation, cpu, sibling)
    isolation = failure.get("isolation")
    if isolation is not None:
        SUPPORT.validate_isolation(isolation, cpu, sibling, require_accepted=False)
        require(pair_lease == isolation["pair_lease"],
                "failure isolation uses another pair lease")
    invocations = failure.get("invocations")
    initial = failure.get("identities_initial")
    execution = failure.get("execution")
    if "inputs_attested" in completed:
        require(isinstance(initial, dict),
                "failure input-attestation phase lacks identities")
    if isinstance(initial, dict) and not allow_self_test:
        validate_retained_snapshot(initial, specification)
    if "executables_staged" in completed:
        require(execution is not None,
                "failure executable-staging phase lacks execution identity")
    if "measurement_started" not in completed:
        require(invocations == [],
                "failure has invocations before measurement started")
    if "measurement_completed" in completed:
        require(isinstance(isolation, dict) and isolation.get("accepted") is True,
                "completed measurement lacks accepted isolation evidence")
    if invocations:
        require(isinstance(initial, dict) and reservation is not None,
                "failed invocation prefix lacks source/build/reservation identity")
        validated_execution = validate_execution_identity(execution, initial)
        validate_invocations(
            invocations, cells, campaign, specification, initial,
            reservation, validated_execution, output, complete=False)
    else:
        require(invocations == [], "failure invocation prefix is invalid")
        require(initial is None or isinstance(initial, dict),
                "failure initial identity is invalid")
        require(execution is None or
                (isinstance(initial, dict) and
                 validate_execution_identity(execution, initial)),
                "failure execution identity is invalid")
    retained = failure.get("retained_files")
    require(isinstance(retained, list), "failure retained file list is invalid")
    seen: set[str] = set()
    for record in retained:
        require(isinstance(record, dict) and isinstance(record.get("path"), str) and
                record["path"] not in seen,
                "failure retained file identities are invalid or duplicated")
        seen.add(record["path"])
        verify_artifact(output, record, "failure retained file")
    require(retained == retained_child_records(output),
            "failure output has unbound, missing, or reordered child artifacts")
    top_level = sorted(path.name for path in output.iterdir())
    expected_top = ["failure.json"] + (["children"] if retained else [])
    require(top_level == sorted(expected_top),
            "failure output contains unbound top-level artifacts")
    return failure


def raw_file_identity(path: Path, payload: Mapping[str, Any]) -> dict[str, Any]:
    metadata = os.lstat(path)
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
            "raw bundle is not a single-link regular file")
    return {
        "path": "raw.json",
        "payload_digest": payload["digest"],
        "sha256": sha256_file(path, MAX_JSON_BYTES),
        "size": metadata.st_size,
    }


def manifest_exit_status(
    expected_status: str, check_current_inputs: bool, allow_self_test: bool
) -> int:
    require(expected_status in {"passed", "policy_failed"},
            "manifest policy status is invalid")
    if not check_current_inputs and not allow_self_test:
        return 1
    return 0 if expected_status == "passed" else 2


def validate_manifest(
    manifest: object,
    output: Path,
    check_current_inputs: bool,
    allow_self_test: bool = False,
) -> tuple[dict[str, Any], int]:
    manifest = verify_signature(manifest, "manifest")
    expected_keys = {
        "analysis", "campaign", "created_utc", "digest", "host", "identities",
        "isolation", "raw", "reservation", "schema", "status", "valid", "execution",
        "validity_is_independent_of_speed",
    }
    require(set(manifest) == expected_keys and manifest.get("schema") == MANIFEST_SCHEMA and
            manifest.get("valid") is True and
            manifest.get("validity_is_independent_of_speed") is True and
            manifest.get("status") in {"passed", "policy_failed"},
            "manifest schema, keys, or status are invalid")
    validate_utc(manifest.get("created_utc"), "manifest created_utc")
    raw_info = manifest.get("raw")
    require(isinstance(raw_info, dict) and set(raw_info) ==
            {"path", "payload_digest", "sha256", "size"},
            "manifest raw identity is invalid")
    raw_path = verify_artifact(output, {
        "path": raw_info.get("path"),
        "sha256": raw_info.get("sha256"),
        "size": raw_info.get("size"),
    }, "raw bundle")
    raw = read_json_limited(raw_path, "raw bundle")
    analysis = validate_raw(raw, output, check_current_inputs, allow_self_test)
    require(raw_info.get("payload_digest") == raw.get("digest"),
            "manifest raw payload digest differs")
    for manifest_name, raw_name in (
        ("campaign", "campaign"), ("host", "host_initial"),
        ("identities", "identities_initial"), ("isolation", "isolation"),
        ("reservation", "reservation"), ("analysis", "analysis"),
        ("execution", "execution"),
    ):
        require(manifest.get(manifest_name) == raw.get(raw_name),
                f"manifest {manifest_name} differs from the raw bundle")
    require(manifest.get("analysis") == analysis,
            "manifest analysis differs from recomputation")
    expected_status = "passed" if analysis["policy"]["passed"] else "policy_failed"
    require(manifest.get("status") == expected_status,
            "manifest status differs from the authenticated promotion policy")
    require(sorted(path.name for path in output.iterdir()) ==
            ["children", "manifest.json", "raw.json"],
            "successful output contains unbound top-level artifacts")
    return manifest, manifest_exit_status(
        expected_status, check_current_inputs, allow_self_test)


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_tree(root: Path) -> None:
    directories = [root]
    directories.extend(path for path in root.rglob("*") if path.is_dir())
    for path in sorted(directories, key=lambda item: len(item.parts), reverse=True):
        fsync_directory(path)


def publish_no_replace(stage: Path, output: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    require(renameat2 is not None,
            "Linux renameat2 is required for no-replace evidence publication")
    renameat2.argtypes = [
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint,
    ]
    renameat2.restype = ctypes.c_int
    at_fdcwd = -100
    rename_noreplace = 1
    result = renameat2(
        at_fdcwd, os.fsencode(stage), at_fdcwd, os.fsencode(output), rename_noreplace)
    if result != 0:
        code = ctypes.get_errno()
        if code == errno.EEXIST:
            raise EvidenceError(f"refusing to replace evidence output {output}")
        raise EvidenceError(f"cannot atomically publish evidence: {os.strerror(code)}")


def create_stage(output: Path) -> Path:
    validate_path_text(str(output.absolute()), "evidence output", absolute=True)
    output = output.resolve()
    require(not output.exists(), f"output path already exists: {output}")
    output.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    require(output.parent.is_dir(), "evidence output parent is not a directory")
    for _ in range(16):
        nonce = hashlib.sha256(os.urandom(32)).hexdigest()[:16]
        stage = output.parent / f".{output.name}.staging-{os.getpid()}-{nonce}"
        try:
            os.mkdir(stage, 0o700)
            os.chmod(stage, 0o700)
            return stage
        except FileExistsError:
            continue
    raise EvidenceError("cannot allocate a unique evidence staging directory")


def publish_manifest(
    stage: Path, output: Path, manifest: Mapping[str, Any],
    check_current_inputs: bool = True, allow_self_test: bool = False,
) -> int:
    _, status = validate_manifest(
        manifest, stage, check_current_inputs=check_current_inputs,
        allow_self_test=allow_self_test)
    fsync_tree(stage)
    fsync_directory(output.parent)
    publish_no_replace(stage, output)
    return status


def publish_failure(
    stage: Path, output: Path, failure: Mapping[str, Any]
) -> None:
    validate_failure(failure, stage)
    fsync_tree(stage)
    fsync_directory(output.parent)
    publish_no_replace(stage, output)


def authoritative_specification(options: argparse.Namespace) -> dict[str, str]:
    taskset = options.taskset.resolve(strict=True)
    ldd = options.ldd.resolve(strict=True)
    git = Path("/usr/bin/git").resolve(strict=True)
    python = Path(sys.executable).resolve(strict=True)
    require(taskset == Path("/usr/bin/taskset").resolve(strict=True),
            "authoritative evidence requires /usr/bin/taskset")
    require(ldd == Path("/usr/bin/ldd").resolve(strict=True),
            "authoritative evidence requires /usr/bin/ldd")
    specification = {
        "candidate_archive": str(options.candidate_archive.resolve(strict=True)),
        "candidate_build_dir": str(options.candidate_build_dir.resolve(strict=True)),
        "candidate_executable": str(options.candidate.resolve(strict=True)),
        "candidate_source_root": str(options.candidate_source_root.resolve(strict=True)),
        "control_archive": str(options.control_archive.resolve(strict=True)),
        "control_build_dir": str(options.control_build_dir.resolve(strict=True)),
        "control_executable": str(options.control.resolve(strict=True)),
        "control_source_root": str(options.control_source_root.resolve(strict=True)),
        "git": str(git),
        "ldd": str(ldd),
        "python": str(python),
        "runner": str(Path(__file__).resolve(strict=True)),
        "taskset": str(taskset),
    }
    validate_input_specification(specification)
    return specification


def validate_exact_affinity(expected: set[int]) -> None:
    require(expected and set(os.sched_getaffinity(0)) == expected,
            "coordinator affinity left the exact housekeeping set")


class TrackedContextExit:
    """Preserve a primary failure while retaining an independent exit failure."""

    def __init__(self, context: Any, label: str):
        self.context = context
        self.label = label
        self.exit_error: str | None = None

    def __enter__(self) -> Any:
        return self.context.__enter__()

    def __exit__(self, exc_type: object, exc: object, tb: object) -> bool:
        try:
            self.context.__exit__(exc_type, exc, tb)
        except Exception as error:
            self.exit_error = f"{self.label}: {type(error).__name__}: {error}"
        return False


def run_campaign(options: argparse.Namespace) -> int:
    require(type(options.iterations) is int and
            5 <= options.iterations <= MAX_ITERATIONS,
            "--iterations is out of bounds")
    require(type(options.warmup) is int and 1 <= options.warmup <= MAX_WARMUP and
            type(options.reuse) is int and 1 <= options.reuse <= MAX_REUSE,
            "--warmup or --reuse is out of bounds")
    require(isinstance(options.timeout, float) and math.isfinite(options.timeout) and
            0 < options.timeout <= MAX_TIMEOUT_SECONDS,
            "--timeout must be positive, finite, and bounded")
    require(type(options.cpu) is int and type(options.reserved_sibling) is int and
            0 <= options.cpu <= SUPPORT.MAX_CPU_ID and
            0 <= options.reserved_sibling <= SUPPORT.MAX_CPU_ID and
            options.cpu != options.reserved_sibling,
            "benchmark CPU pair is invalid")
    allowed = set(os.sched_getaffinity(0))
    require(options.cpu in allowed and options.reserved_sibling in allowed and
            bool(allowed - {options.cpu, options.reserved_sibling}),
            "launch affinity must contain the pair and at least one housekeeping CPU")
    for cell in FIXED_CELLS:
        validate_cell(cell)
    specification = authoritative_specification(options)
    campaign: dict[str, Any] = {
        "allowed_cpu_set_at_launch": sorted(allowed),
        "batch": 1,
        "benchmark_cpu": options.cpu,
        "cells": [asdict(cell) for cell in FIXED_CELLS],
        "child_environment": dict(CHILD_ENVIRONMENT),
        "iterations": options.iterations,
        "mode": "authoritative",
        "reserved_sibling": options.reserved_sibling,
        "reuse": options.reuse,
        "rounds": ROUNDS,
        "sequence": [list(item) for item in SEQUENCE],
        "statistics": statistics_policy(),
        "threads": 1,
        "timeout_seconds": options.timeout,
        "warmup": options.warmup,
    }
    validate_campaign(campaign)
    output = options.output.resolve()
    stage = create_stage(output)
    original_affinity = set(os.sched_getaffinity(0))
    host_initial: dict[str, Any] | None = None
    reservation: dict[str, Any] | None = None
    pair_lease: dict[str, Any] | None = None
    initial: dict[str, Any] | None = None
    isolation: dict[str, Any] | None = None
    execution_guard: ImmutableExecutables | None = None
    execution: dict[str, dict[str, Any]] | None = None
    execution_descriptors: dict[str, int] = {}
    invocations: list[dict[str, Any]] = []
    manifest: dict[str, Any] | None = None
    completed_phases = ["initialized"]
    caught: Exception | None = None
    caught_traceback = ""
    teardown_errors: list[str] = []
    pair_context: TrackedContextExit | None = None
    reservation_context: TrackedContextExit | None = None
    try:
        allowed_checked, housekeeping = SUPPORT.validate_topology(
            options.cpu, options.reserved_sibling)
        require(allowed_checked == allowed,
                "launch affinity changed during topology validation")
        host_initial = SUPPORT.host_identity(
            options.cpu, options.reserved_sibling, allowed)
        completed_phases.append("topology_validated")
        pair_guard = SUPPORT.PairLease(options.cpu, options.reserved_sibling)
        reservation_guard = HardenedReservation(
            options.reservation_file, options.cpu, options.reserved_sibling)
        pair_context = TrackedContextExit(pair_guard, "CPU pair lease teardown")
        reservation_context = TrackedContextExit(
            reservation_guard, "CPU reservation teardown")
        with pair_context as pair_lease, reservation_context as reservation:
            completed_phases.append("locks_acquired")
            os.sched_setaffinity(0, housekeeping)
            completed_phases.append("affinity_isolated")

            def validate_guards() -> None:
                pair_guard.validate_current()
                reservation_guard.validate_current()
                validate_exact_affinity(housekeeping)

            validate_guards()
            initial = input_snapshot(specification)
            completed_phases.append("inputs_attested")
            execution_guard = ImmutableExecutables(stage, specification, initial)
            execution, execution_descriptors = execution_guard.__enter__()
            validate_execution_identity(execution, initial)
            completed_phases.append("executables_staged")

            def snapshot() -> Mapping[str, Any]:
                return input_snapshot(specification)

            completed_phases.append("measurement_started")
            before_monotonic_ns = time.monotonic_ns()
            before_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = SUPPORT.cpu_stat_snapshot(options.reserved_sibling)
            try:
                for cell in FIXED_CELLS:
                    for round_index in range(ROUNDS):
                        for slot, (label, implementation) in enumerate(SEQUENCE):
                            invocation = run_child(
                                implementation, label, cell, round_index, slot,
                                campaign, specification, initial, reservation,
                                execution, execution_descriptors, stage,
                                options.cpu, options.timeout, snapshot, validate_guards)
                            validate_guards()
                            invocations.append(invocation)
            finally:
                after_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
                after_sibling = SUPPORT.cpu_stat_snapshot(options.reserved_sibling)
                after_monotonic_ns = time.monotonic_ns()
                isolation = SUPPORT.isolation_record(
                    options.cpu, options.reserved_sibling, pair_lease,
                    before_monotonic_ns, after_monotonic_ns,
                    before_cpu, after_cpu, before_sibling, after_sibling)
                validate_guards()
            require(isolation["accepted"] is True,
                    "reserved SMT sibling performed non-idle work during the campaign")
            completed_phases.append("measurement_completed")
            final = input_snapshot(specification)
            require(final == initial, "input identity changed during campaign")
            host_final = SUPPORT.host_identity(
                options.cpu, options.reserved_sibling, allowed)
            require(host_final == host_initial,
                    "host topology or frequency policy changed during campaign")
            analysis = analyze(invocations, FIXED_CELLS)
            raw = signed({
                "analysis": analysis,
                "campaign": campaign,
                "created_utc": utc_now(),
                "host_final": host_final,
                "host_initial": host_initial,
                "identities_final": final,
                "identities_initial": initial,
                "execution": execution,
                "input_specification": specification,
                "invocations": invocations,
                "isolation": isolation,
                "reservation": reservation,
                "schema": RAW_SCHEMA,
                "validity_is_independent_of_speed": True,
            })
            validate_guards()
            validate_raw(raw, stage, check_current_inputs=True)
            raw_path = stage / "raw.json"
            write_json_exclusive(raw_path, raw)
            status = "passed" if analysis["policy"]["passed"] else "policy_failed"
            manifest = signed({
                "analysis": analysis,
                "campaign": campaign,
                "created_utc": utc_now(),
                "host": host_initial,
                "identities": initial,
                "execution": execution,
                "isolation": isolation,
                "raw": raw_file_identity(raw_path, raw),
                "reservation": reservation,
                "schema": MANIFEST_SCHEMA,
                "status": status,
                "valid": True,
                "validity_is_independent_of_speed": True,
            })
            write_json_exclusive(stage / "manifest.json", manifest)
            validate_guards()
            validate_manifest(
                manifest, stage, check_current_inputs=True, allow_self_test=False)
            completed_phases.append("evidence_prepared")
    except Exception as error:
        caught = error
        caught_traceback = traceback.format_exc()[-65536:]
    finally:
        if execution_guard is not None:
            try:
                execution_guard.__exit__(None, None, None)
            except Exception as error:
                teardown_errors.append(
                    f"immutable executable teardown: {type(error).__name__}: {error}")
            execution_guard = None
        for tracked in (reservation_context, pair_context):
            if tracked is not None and tracked.exit_error is not None:
                teardown_errors.append(tracked.exit_error)
        try:
            os.sched_setaffinity(0, original_affinity)
            validate_exact_affinity(original_affinity)
        except Exception as error:
            teardown_errors.append(
                f"affinity restoration: {type(error).__name__}: {error}")

    if caught is None and (
        manifest is None or completed_phases != list(FAILURE_PHASES)
    ):
        caught = EvidenceError(
            "campaign teardown reached publication without complete prepared evidence")
        caught_traceback = str(caught)

    primary_failed_phase = (
        FAILURE_PHASES[len(completed_phases)]
        if caught is not None and len(completed_phases) < len(FAILURE_PHASES)
        else None
    )
    if teardown_errors:
        teardown_error = "; ".join(teardown_errors)[:8192]
        if caught is None:
            caught = EvidenceError(teardown_error)
            caught_traceback = teardown_error
        else:
            caught_traceback = (caught_traceback + "\n" + teardown_error)[-65536:]
    else:
        teardown_error = None

    def remove_staged_summary_files() -> None:
        for name in ("manifest.json", "raw.json", "failure.json"):
            path = stage / name
            if not path.exists():
                continue
            metadata = os.lstat(path)
            require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
                    f"cannot clean unsafe staged artifact {path}")
            path.unlink()

    if caught is not None:
        try:
            remove_staged_summary_files()
            failed_phase = primary_failed_phase or (
                "teardown" if teardown_errors else
                (FAILURE_PHASES[len(completed_phases)]
                 if len(completed_phases) < len(FAILURE_PHASES) else "teardown"))
            failure = signed({
                "campaign": campaign,
                "created_utc": utc_now(),
                "error": (str(caught) or repr(caught))[:8192],
                "error_type": type(caught).__name__[:256] or "Exception",
                "host_initial": host_initial,
                "identities_initial": initial,
                "execution": execution,
                "input_specification": specification,
                "invocations": invocations,
                "isolation": isolation,
                "lifecycle": {
                    "completed_phases": completed_phases,
                    "failed_phase": failed_phase,
                    "teardown_error": teardown_error,
                    "teardown_status": "failed" if teardown_errors else "completed",
                },
                "pair_lease": pair_lease,
                "reservation": reservation,
                "retained_files": retained_child_records(stage),
                "schema": FAILURE_SCHEMA,
                "status": "failed",
                "traceback": (caught_traceback or repr(caught))[-65536:],
                "valid": False,
            })
            write_json_exclusive(stage / "failure.json", failure)
            publish_failure(stage, output, failure)
            stage = Path()
        except Exception:
            if stage != Path() and stage.exists():
                shutil.rmtree(stage)
            raise
        if isinstance(caught, EvidenceError):
            raise caught
        raise EvidenceError(
            f"campaign failed: {type(caught).__name__}: {caught}") from caught

    require(manifest is not None,
            "campaign publication lost its prepared manifest")
    try:
        result = publish_manifest(stage, output, manifest)
    except Exception as error:
        if output.exists():
            if stage.exists():
                shutil.rmtree(stage)
            raise EvidenceError(
                "publication failed because the output path is already visible") from error
        remove_staged_summary_files()
        failure = signed({
            "campaign": campaign,
            "created_utc": utc_now(),
            "error": (str(error) or repr(error))[:8192],
            "error_type": type(error).__name__[:256] or "Exception",
            "host_initial": host_initial,
            "identities_initial": initial,
            "execution": execution,
            "input_specification": specification,
            "invocations": invocations,
            "isolation": isolation,
            "lifecycle": {
                "completed_phases": completed_phases,
                "failed_phase": "publication",
                "teardown_error": None,
                "teardown_status": "completed",
            },
            "pair_lease": pair_lease,
            "reservation": reservation,
            "retained_files": retained_child_records(stage),
            "schema": FAILURE_SCHEMA,
            "status": "failed",
            "traceback": traceback.format_exc()[-65536:],
            "valid": False,
        })
        write_json_exclusive(stage / "failure.json", failure)
        publish_failure(stage, output, failure)
        stage = Path()
        raise EvidenceError(f"campaign publication failed: {error}") from error
    stage = Path()
    return result


def verify_campaign(options: argparse.Namespace) -> int:
    manifest_path = top_level_evidence_path(options.manifest, "manifest path")
    manifest = read_json_limited(manifest_path, "manifest")
    _, status = validate_manifest(
        manifest, manifest_path.parent,
        check_current_inputs=not options.no_current_input_check)
    if options.no_current_input_check:
        print(
            "low-encode-copy bundle is structurally consistent only; live input "
            "provenance was not rebound, so this is not authoritative evidence")
        return 1
    if status == 0:
        print("low-encode-copy evidence verified and promotion policy passed")
    else:
        print("low-encode-copy evidence verified; promotion policy failed")
    return status


def verify_failed_campaign(options: argparse.Namespace) -> int:
    failure_path = top_level_evidence_path(options.failure, "failure path")
    failure = read_json_limited(failure_path, "failure bundle")
    validate_failure(failure, failure_path.parent)
    print("failed low-encode-copy campaign diagnostics verified; not performance evidence")
    return 1


def find_self_test_pair() -> tuple[int, int, set[int]]:
    allowed = set(os.sched_getaffinity(0))
    for cpu in sorted(allowed):
        path = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
        if not path.is_file():
            continue
        siblings = SUPPORT.parse_cpu_list(path.read_text(encoding="ascii"))
        if len(siblings) == 2 and siblings <= allowed and allowed - siblings:
            first, second = sorted(siblings)
            return first, second, allowed
    raise EvidenceError("self-test needs an allowed two-thread SMT pair and housekeeping CPU")


def synthetic_pair_lease(cpu: int, sibling: int) -> dict[str, Any]:
    payload = SUPPORT.pair_lease_payload(cpu, sibling)
    return {
        "device": 1,
        "directory_device": 1,
        "directory_inode": 2,
        "inode": 3,
        "lock": "exclusive_nonblocking_pair_wide",
        "path": str(SUPPORT.pair_lease_directory() /
                    SUPPORT.pair_lease_name(cpu, sibling)),
        "payload": payload,
        "sha256": sha256_bytes(canonical_bytes(payload)),
    }


def synthetic_reservation(root: Path, cpu: int, sibling: int) -> dict[str, Any]:
    payload = {
        "benchmark_cpu": cpu,
        "nonce": "0123456789abcdef0123456789abcdef",
        "owner": "low-copy-self-test",
        "reserved_sibling": sibling,
        "schema": RESERVATION_SCHEMA,
        "status": "held",
    }
    return {
        "device": 1,
        "directory_device": 1,
        "directory_inode": 2,
        "inode": 3,
        "lock": "exclusive_nonblocking_inode_bound",
        "mode": 0o600,
        "path": str((root / "reservation.json").resolve()),
        "payload": payload,
        "sha256": sha256_bytes(canonical_bytes(payload)),
    }


def synthetic_isolation(
    cpu: int, sibling: int, pair_lease: Mapping[str, Any], elapsed_ns: int
) -> dict[str, Any]:
    def snapshot(selected: int, user: int, idle: int) -> dict[str, Any]:
        fields = {name: 100 for name in SUPPORT.CPU_STAT_FIELDS}
        fields["user"] += user
        fields["idle"] += idle
        return {"cpu": selected, "fields": fields, "total_jiffies": sum(fields.values())}

    return SUPPORT.isolation_record(
        cpu, sibling, pair_lease, 1000, 1000 + max(elapsed_ns, 1),
        snapshot(cpu, 0, 0), snapshot(cpu, 10, 10),
        snapshot(sibling, 0, 0), snapshot(sibling, 0, 20))


def self_test_campaign(
    cells: Sequence[Cell], cpu: int, sibling: int, allowed: set[int]
) -> dict[str, Any]:
    return {
        "allowed_cpu_set_at_launch": sorted(allowed),
        "batch": 1,
        "benchmark_cpu": cpu,
        "cells": [asdict(cell) for cell in cells],
        "child_environment": dict(CHILD_ENVIRONMENT),
        "iterations": 5,
        "mode": "self_test" if len(cells) == 1 else "self_test_full",
        "reserved_sibling": sibling,
        "reuse": 2,
        "rounds": ROUNDS,
        "sequence": [list(item) for item in SEQUENCE],
        "statistics": statistics_policy(),
        "threads": 1,
        "timeout_seconds": 10.0,
        "warmup": 1,
    }


def build_self_test_bundle(
    root: Path,
    candidate: Path,
    cpu: int,
    sibling: int,
    allowed: set[int],
    cells: Sequence[Cell] = (FIXED_CELLS[0],),
) -> tuple[dict[str, Any], dict[str, Any]]:
    output = root
    output.mkdir(mode=0o700)
    retained_cells = tuple(copy.deepcopy(tuple(cells)))
    campaign = self_test_campaign(retained_cells, cpu, sibling, allowed)
    validate_campaign(campaign, allow_self_test=True)
    control = root.parent / "control-mock"
    specification = {
        "candidate_archive": str(root.parent / "candidate.a"),
        "candidate_build_dir": str(root.parent / "candidate-build"),
        "candidate_executable": str(candidate),
        "candidate_source_root": str(root.parent / "candidate-source"),
        "control_archive": str(root.parent / "control.a"),
        "control_build_dir": str(root.parent / "control-build"),
        "control_executable": str(control),
        "control_source_root": str(root.parent / "control-source"),
        "git": str(Path("/usr/bin/git").resolve(strict=True)),
        "ldd": "/usr/bin/ldd",
        "python": str(Path(sys.executable).resolve(strict=True)),
        "runner": str(Path(__file__).resolve()),
        "taskset": "/usr/bin/taskset",
    }
    initial = {
        "mock_identity": "stable",
        "control_executable": SUPPORT.artifact_identity(control, "executable"),
        "candidate_executable": SUPPORT.artifact_identity(candidate, "executable"),
    }
    reservation = synthetic_reservation(root.parent, cpu, sibling)
    invocations: list[dict[str, Any]] = []
    snapshot = lambda: initial
    validate_guards = lambda: None
    execution_guard = ImmutableExecutables(output, specification, initial)
    with execution_guard as (execution, execution_descriptors):
        for cell in retained_cells:
            for round_index in range(ROUNDS):
                for slot, (label, implementation) in enumerate(SEQUENCE):
                    invocations.append(run_child(
                        implementation, label, cell, round_index, slot, campaign,
                        specification, initial, reservation, execution,
                        execution_descriptors, output, cpu, 10.0,
                        snapshot, validate_guards))
    pair_lease = synthetic_pair_lease(cpu, sibling)
    isolation = synthetic_isolation(
        cpu, sibling, pair_lease,
        sum(item["duration_ns"] for item in invocations) + 1)
    host = SUPPORT.host_identity(cpu, sibling, allowed)
    analysis = analyze(invocations, retained_cells)
    raw = signed({
        "analysis": analysis,
        "campaign": campaign,
        "created_utc": "2000-01-01T00:00:00Z",
        "host_final": host,
        "host_initial": host,
        "identities_final": initial,
        "identities_initial": initial,
        "execution": execution,
        "input_specification": specification,
        "invocations": invocations,
        "isolation": isolation,
        "reservation": reservation,
        "schema": RAW_SCHEMA,
        "validity_is_independent_of_speed": True,
    })
    validate_raw(raw, output, check_current_inputs=False, allow_self_test=True)
    raw_path = output / "raw.json"
    write_json_exclusive(raw_path, raw)
    status = "passed" if analysis["policy"]["passed"] else "policy_failed"
    manifest = signed({
        "analysis": analysis,
        "campaign": campaign,
        "created_utc": "2000-01-01T00:00:00Z",
        "host": host,
        "identities": initial,
        "execution": execution,
        "isolation": isolation,
        "raw": raw_file_identity(raw_path, raw),
        "reservation": reservation,
        "schema": MANIFEST_SCHEMA,
        "status": status,
        "valid": True,
        "validity_is_independent_of_speed": True,
    })
    write_json_exclusive(output / "manifest.json", manifest)
    _, manifest_status = validate_manifest(
        manifest, output, check_current_inputs=False, allow_self_test=True)
    require(manifest_status == (0 if status == "passed" else 2),
            "self-test manifest returned the wrong status")
    return raw, manifest


def expect_evidence_error(action: Callable[[], object], name: str) -> None:
    try:
        action()
    except (EvidenceError, FileNotFoundError, NotADirectoryError, OSError, ValueError,
            TypeError, KeyError, AttributeError, RuntimeError):
        return
    raise EvidenceError(f"adversarial mutation was accepted: {name}")


def resign(value: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(value))
    payload.pop("digest", None)
    return signed(payload)


def run_self_test(_options: argparse.Namespace) -> int:
    validate_support_contract()
    for seed, k, count, expected in LOSS_VECTOR_FIXTURES:
        fixture = Cell(
            identifier="loss-vector", backend="scalar", field="gf16",
            k=k, r=max(1, count), shard_bytes=64, role="target",
            seed=seed, losses=count)
        require(tuple(expected_losses(fixture)) == expected,
                f"benchmark-compatible loss vector changed for seed {seed:#x}")
    for cell in FIXED_CELLS:
        validate_cell(cell)
    cpu, sibling, allowed = find_self_test_pair()
    with tempfile.TemporaryDirectory(prefix="leopard2-low-copy-selftest-") as temporary:
        root = Path(temporary)
        os.chmod(root, 0o700)
        mock = Path(__file__).with_name("mock_benchmark.py")
        control = root / "control-mock"
        fast = root / "candidate-fast-mock"
        slow = root / "candidate-slow-mock"
        for destination, role in (
            (control, "control"), (fast, "candidate-fast"),
            (slow, "candidate-slow"),
        ):
            destination.write_bytes(
                mock.read_bytes() + f"\n# MOCK_ROLE={role}\n".encode("ascii"))
            os.chmod(destination, 0o700)
        fast_raw, fast_manifest = build_self_test_bundle(
            root / "fast", fast, cpu, sibling, allowed)
        slow_raw, slow_manifest = build_self_test_bundle(
            root / "slow", slow, cpu, sibling, allowed)
        full_raw, full_manifest = build_self_test_bundle(
            root / "full-matrix", fast, cpu, sibling, allowed,
            cells=FIXED_CELLS)
        require(fast_raw["analysis"]["policy"]["passed"] is True and
                slow_raw["analysis"]["policy"]["passed"] is False,
                "mock pass/fail policy separation failed")
        require(len(full_raw["invocations"]) ==
                len(FIXED_CELLS) * ROUNDS * len(SEQUENCE) and
                full_raw["analysis"]["policy"]["passed"] is True,
                "full 24-cell mock campaign did not execute or pass")
        _, fast_status = validate_manifest(
            fast_manifest, root / "fast", False, allow_self_test=True)
        _, slow_status = validate_manifest(
            slow_manifest, root / "slow", False, allow_self_test=True)
        _, full_status = validate_manifest(
            full_manifest, root / "full-matrix", False, allow_self_test=True)
        require((fast_status, slow_status) == (0, 2),
                "pass/policy-failure exit semantics are not distinct")
        require(full_status == 0, "full mock matrix has the wrong status")
        require(manifest_exit_status("passed", False, False) == 1 and
                manifest_exit_status("policy_failed", False, False) == 1 and
                manifest_exit_status("passed", True, False) == 0 and
                manifest_exit_status("policy_failed", True, False) == 2,
                "portable replay was not downgraded to non-authoritative status")

        failed_root = root / "failed"
        failed_root.mkdir(mode=0o700)
        failure = signed({
            "campaign": fast_raw["campaign"],
            "created_utc": "2000-01-01T00:00:00Z",
            "error": "intentional self-test failure",
            "error_type": "EvidenceError",
            "execution": None,
            "host_initial": None,
            "identities_initial": None,
            "input_specification": fast_raw["input_specification"],
            "invocations": [],
            "isolation": None,
            "lifecycle": {
                "completed_phases": ["initialized"],
                "failed_phase": "topology_validated",
                "teardown_error": None,
                "teardown_status": "completed",
            },
            "pair_lease": None,
            "reservation": None,
            "retained_files": [],
            "schema": FAILURE_SCHEMA,
            "status": "failed",
            "traceback": "intentional self-test traceback",
            "valid": False,
        })
        write_json_exclusive(failed_root / "failure.json", failure)
        validate_failure(failure, failed_root, allow_self_test=True)
        bad_lifecycle = copy.deepcopy(failure)
        bad_lifecycle["lifecycle"]["completed_phases"] = [
            "initialized", "locks_acquired"]
        bad_lifecycle = resign(bad_lifecycle)
        expect_evidence_error(
            lambda: validate_failure(
                bad_lifecycle, failed_root, allow_self_test=True),
            "failure lifecycle prefix")
        bad_teardown = copy.deepcopy(failure)
        bad_teardown["lifecycle"]["teardown_status"] = "failed"
        bad_teardown["lifecycle"]["teardown_error"] = None
        bad_teardown = resign(bad_teardown)
        expect_evidence_error(
            lambda: validate_failure(
                bad_teardown, failed_root, allow_self_test=True),
            "failure teardown state")
        combined_failure = copy.deepcopy(failure)
        combined_failure["lifecycle"]["teardown_status"] = "failed"
        combined_failure["lifecycle"]["teardown_error"] = \
            "intentional teardown failure"
        combined_failure = resign(combined_failure)
        validate_failure(combined_failure, failed_root, allow_self_test=True)

        class BrokenExit:
            def __enter__(self) -> None:
                return None

            def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
                raise EvidenceError("intentional context-exit failure")

        tracked_exit = TrackedContextExit(BrokenExit(), "self-test teardown")
        try:
            with tracked_exit:
                raise EvidenceError("intentional primary failure")
        except EvidenceError as error:
            require(str(error) == "intentional primary failure" and
                    tracked_exit.exit_error is not None and
                    "context-exit failure" in tracked_exit.exit_error,
                    "tracked teardown obscured the primary failure")
        else:
            raise EvidenceError("tracked context suppressed a primary failure")

        mutations: list[tuple[str, Callable[[dict[str, Any]], None]]] = [
            ("command", lambda value: value["invocations"][0]["command"].__setitem__(
                2, str(sibling))),
            ("boolean round", lambda value: value["invocations"][0].__setitem__(
                "round", False)),
            ("duration", lambda value: value["invocations"][0].__setitem__(
                "duration_ns", 0)),
            ("wire digest", lambda value: value["invocations"][0]["result"]
             ["workload_digests"].__setitem__("transmitted_parity", "f" * 16)),
            ("recovery digest", lambda value: value["invocations"][0]["result"]
             ["workload_digests"].__setitem__("recovered_originals", "e" * 16)),
            ("analysis", lambda value: value["analysis"]["cells"]
             [FIXED_CELLS[0].identifier]["encode_execution"].__setitem__(
                 "ci95_lower", 99.0)),
            ("reservation", lambda value: value["reservation"]["payload"].__setitem__(
                "benchmark_cpu", sibling)),
            ("isolation", lambda value: value["isolation"].__setitem__(
                "accepted", False)),
            ("path traversal", lambda value: value["invocations"][0]["stdout"]
             .__setitem__("path", "../escape")),
            ("nested null", lambda value: value.__setitem__("campaign", None)),
        ]
        for name, mutation in mutations:
            changed = copy.deepcopy(fast_raw)
            mutation(changed)
            changed = resign(changed)
            expect_evidence_error(
                lambda changed=changed: validate_raw(
                    changed, root / "fast", False, allow_self_test=True), name)
        unsigned = copy.deepcopy(fast_raw)
        unsigned["analysis"] = None
        expect_evidence_error(
            lambda: validate_raw(unsigned, root / "fast", False, True),
            "raw signature")
        changed_manifest = copy.deepcopy(fast_manifest)
        changed_manifest["status"] = "policy_failed"
        changed_manifest = resign(changed_manifest)
        expect_evidence_error(
            lambda: validate_manifest(
                changed_manifest, root / "fast", False, allow_self_test=True),
            "manifest status")
        expect_evidence_error(
            lambda: parse_json_bytes(b'"' + b'x' * (MAX_JSON_STRING_BYTES + 1) + b'"',
                                     "oversized string"),
            "bounded JSON string")
        expect_evidence_error(
            lambda: parse_json_bytes(b"[" * 64 + b"0" + b"]" * 64,
                                     "deep JSON"),
            "bounded JSON depth")
        huge_integer = ("{\"value\":" + "9" * 100 + "}").encode("ascii")
        expect_evidence_error(
            lambda: parse_json_bytes(huge_integer, "oversized integer"),
            "bounded JSON integer")
        expect_evidence_error(
            lambda: SUPPORT.parse_cpu_list("0-1000000000"),
            "bounded CPU range")
        bad_cell = Cell(**{
            **asdict(FIXED_CELLS[0]), "losses": True,
        })
        expect_evidence_error(lambda: validate_cell(bad_cell), "boolean cell integer")
        bad_reservation = {
            "benchmark_cpu": False,
            "nonce": "f" * 32,
            "owner": "self-test",
            "reserved_sibling": 16,
            "schema": RESERVATION_SCHEMA,
            "status": "held",
        }
        expect_evidence_error(
            lambda: HardenedReservation.parse(
                canonical_bytes(bad_reservation), 0, 16),
            "boolean reservation CPU")
        bad_lease = synthetic_pair_lease(0, 16)
        bad_lease["payload"]["cpus"][0] = False
        bad_lease["sha256"] = sha256_bytes(canonical_bytes(bad_lease["payload"]))
        expect_evidence_error(
            lambda: SUPPORT.validate_pair_lease_identity(
                bad_lease, 0, 16), "boolean pair-lease CPU")
        distinct_domains = validate_digests({
            "algorithm": "fnv1a64", "original_data": "0" * 16,
            "transmitted_parity": "1" * 16,
            "recovered_originals": "2" * 16,
        })
        require(len(set(distinct_domains.values())) == 3,
                "benchmark-v2 digest domains were conflated")
        expect_evidence_error(
            lambda: validate_digests({
                "algorithm": "fnv1a64", "original_data": "0" * 16,
                "transmitted_parity": "1" * 16,
                "recovered_originals": "not-a-digest",
            }), "malformed recovery digest")
        expect_evidence_error(
            lambda: validate_relinked_executable_sha("0" * 64, "1" * 64, "mock"),
            "stale executable relink")
        synthetic_pair = {
            "source": {"path": "/mock/source.cpp", "sha256": "1" * 64},
            "object": {"path": "/mock/object.o", "sha256": "2" * 64, "size": 7},
        }
        synthetic_compiler = {"path": "/mock/c++", "sha256": "3" * 64}
        synthetic_recipe = ["<compiler>", "-O3", "<source>", "-o", "<object>"]
        synthetic_proof = {
            "cache_key": sha256_bytes(canonical_bytes({
                "compiler": synthetic_compiler,
                "normalized_recipe": synthetic_recipe,
                "object": synthetic_pair["object"],
                "source": synthetic_pair["source"],
            })),
            "object_sha256": "2" * 64,
            "object_size": 7,
            "source": "source.cpp",
            "status": "byte_identical_clean_recompile",
        }
        validate_clean_recompile_record(
            synthetic_proof, "source.cpp", synthetic_pair,
            synthetic_compiler, synthetic_recipe)
        bad_proof = copy.deepcopy(synthetic_proof)
        bad_proof["object_sha256"] = "4" * 64
        expect_evidence_error(
            lambda: validate_clean_recompile_record(
                bad_proof, "source.cpp", synthetic_pair,
                synthetic_compiler, synthetic_recipe),
            "clean-recompile digest mismatch")
        control_policy = {
            "compiler": {"sha256": "0" * 64},
            "compiler_version_stdout": {"text": "mock"},
            "validated_cache": {},
            "strict_compile_recipes": {"mock.cpp": ["-O3"]},
            "strict_link_provenance": {
                "normalized_executable_link": ["-O3"],
                "external_link_inputs": [], "ranlib": {"sha256": "0" * 64},
            },
        }
        candidate_policy = copy.deepcopy(control_policy)
        candidate_policy["strict_compile_recipes"]["mock.cpp"].append("-fno-inline")
        expect_evidence_error(
            lambda: validate_equal_build_policy(control_policy, candidate_policy),
            "mismatched effective compile flags")
        symlink = root / "artifact-link"
        symlink.symlink_to(root / "fast" / "raw.json")
        expect_evidence_error(
            lambda: safe_evidence_path(root, "artifact-link"),
            "retained artifact symlink")
        top_level = root / "top-level.json"
        write_bytes_exclusive(top_level, b"{}")
        top_symlink = root / "top-level-symlink.json"
        top_symlink.symlink_to(top_level)
        expect_evidence_error(
            lambda: top_level_evidence_path(top_symlink, "top-level symlink"),
            "top-level evidence symlink")
        top_hardlink = root / "top-level-hardlink.json"
        os.link(top_level, top_hardlink)
        expect_evidence_error(
            lambda: top_level_evidence_path(top_hardlink, "top-level hardlink"),
            "top-level evidence hardlink")
        top_fifo = root / "top-level.fifo"
        os.mkfifo(top_fifo, 0o600)
        expect_evidence_error(
            lambda: read_json_limited(top_fifo, "top-level FIFO"),
            "nonblocking FIFO input")
        expect_evidence_error(
            lambda: top_level_evidence_path(root, "top-level directory"),
            "top-level evidence directory")
        expect_evidence_error(
            lambda: top_level_evidence_path(Path("/dev/null"), "device input"),
            "top-level evidence device")
        empty = root / "fast" / "children" / "unbound-empty-directory"
        empty.mkdir(mode=0o700)
        expect_evidence_error(
            lambda: validate_manifest(
                fast_manifest, root / "fast", False, allow_self_test=True),
            "unbound empty artifact directory")
        empty.rmdir()
        oversized_stderr = root / "oversized.stderr.txt"
        write_bytes_exclusive(oversized_stderr, b"x" * (MAX_STDERR_BYTES + 1))
        oversized_record = artifact_record(root, oversized_stderr)
        expect_evidence_error(
            lambda: verify_artifact(
                root, oversized_record, "oversized stderr", MAX_STDERR_BYTES),
            "retained stderr cap")
        try:
            run_bounded(
                ("/usr/bin/python3", "-c",
                 f"import sys;sys.stdout.write('x'*{MAX_STDOUT_BYTES + 65536})"),
                CHILD_ENVIRONMENT, 5.0)
        except BoundedChildError as error:
            require(len(error.stdout) == MAX_STDOUT_BYTES,
                    "bounded child did not retain the capped stdout prefix")
        else:
            raise EvidenceError("oversized child stdout was accepted")
        try:
            run_bounded(
                ("/usr/bin/python3", "-c", "while True: pass"),
                CHILD_ENVIRONMENT, 0.05)
        except BoundedChildError as error:
            require("exceeded" in str(error), "child timeout used the wrong error")
        else:
            raise EvidenceError("child timeout was accepted")
        descendant_marker = root / "escaped-descendant"
        descendant_program = (
            "import os,pathlib,time\n"
            "child=os.fork()\n"
            "if child == 0:\n"
            " time.sleep(.25);pathlib.Path(" + repr(str(descendant_marker)) +
            ").write_text('escaped')\n"
            "else:\n"
            " while True: pass\n"
        )
        try:
            run_bounded(
                (str(Path(sys.executable).resolve()), "-c", descendant_program),
                CHILD_ENVIRONMENT, 0.05)
        except BoundedChildError:
            pass
        else:
            raise EvidenceError("descendant timeout was accepted")
        time.sleep(0.4)
        require(not descendant_marker.exists(),
                "timed-out benchmark descendant escaped process-group cleanup")

        class NeverReaps:
            pid = 999_999_999
            returncode: int | None = None

            def __init__(self) -> None:
                self.wait_timeouts: list[float | None] = []

            def wait(self, timeout: float | None = None) -> int:
                self.wait_timeouts.append(timeout)
                raise subprocess.TimeoutExpired(("never-reaps",), timeout)

        never_reaps = NeverReaps()
        reaped, _ = terminate_process_group_bounded(never_reaps)  # type: ignore[arg-type]
        require(not reaped and never_reaps.wait_timeouts == [5.0] and
                all(value is not None for value in never_reaps.wait_timeouts),
                "benchmark cleanup performed an unbounded wait")
        support_never_reaps = NeverReaps()
        support_reaped, _ = SUPPORT.terminate_process_group_bounded(
            support_never_reaps)  # type: ignore[arg-type]
        require(not support_reaped and
                support_never_reaps.wait_timeouts == [5.0] and
                all(value is not None for value in support_never_reaps.wait_timeouts),
                "support cleanup performed an unbounded wait")
        expect_evidence_error(
            lambda: SUPPORT.run_process_bounded(
                ("/usr/bin/python3", "-c",
                 "import sys;sys.stdout.buffer.write(b'x'*2048)"),
                environment=CHILD_ENVIRONMENT, timeout=5.0,
                max_stdout=1024, max_stderr=1024),
            "support subprocess output cap")

        reservation_root = root / "reservation-root"
        reservation_root.mkdir(mode=0o700)
        reservation_path = reservation_root / "reservation.json"
        payload = {
            "benchmark_cpu": cpu,
            "nonce": "fedcba9876543210fedcba9876543210",
            "owner": "self-test",
            "reserved_sibling": sibling,
            "schema": RESERVATION_SCHEMA,
            "status": "held",
        }
        write_bytes_exclusive(reservation_path, canonical_bytes(payload))
        reservation_guard = HardenedReservation(reservation_path, cpu, sibling)
        with reservation_guard:
            replacement = reservation_root / "replacement.json"
            write_bytes_exclusive(replacement, canonical_bytes(payload))
            os.replace(replacement, reservation_path)
            expect_evidence_error(
                lambda: HardenedReservation.validate_current(
                    reservation_guard), "reservation replacement")
            replacement_guard = HardenedReservation(
                reservation_path, cpu, sibling)
            expect_evidence_error(
                lambda: replacement_guard.__enter__(),
                "reservation kernel-lease replacement overlap")

        lease_root = root / "lease-root"
        lease_root.mkdir(mode=0o700)
        first = SUPPORT.PairLease(cpu, sibling, root=lease_root)
        second = SUPPORT.PairLease(sibling, cpu, root=lease_root)
        with first:
            expect_evidence_error(lambda: second.__enter__(), "pair lease overlap")

        replacement_lease_root = root / "replacement-lease-root"
        replacement_lease_root.mkdir(mode=0o700)
        replacement_first = SUPPORT.PairLease(
            cpu, sibling, root=replacement_lease_root)
        with replacement_first:
            replacement_first.path.rename(replacement_lease_root / "old.lock")
            expect_evidence_error(
                lambda: SUPPORT.PairLease(
                    cpu, sibling, root=replacement_lease_root).__enter__(),
                "pair lease file replacement overlap")

        directory_lease_root = root / "directory-lease-root"
        directory_lease_root.mkdir(mode=0o700)
        directory_first = SUPPORT.PairLease(cpu, sibling, root=directory_lease_root)
        with directory_first:
            moved_root = root / "directory-lease-root-moved"
            directory_lease_root.rename(moved_root)
            directory_lease_root.mkdir(mode=0o700)
            expect_evidence_error(
                lambda: SUPPORT.PairLease(
                    cpu, sibling, root=directory_lease_root).__enter__(),
                "pair lease directory replacement overlap")

        butterfly_path = Path(__file__).resolve().parents[1] / \
            "backend_butterfly" / "run_abba.py"
        butterfly_spec = importlib.util.spec_from_file_location(
            "leopard2_low_copy_butterfly_lease_test", butterfly_path)
        require(butterfly_spec is not None and butterfly_spec.loader is not None,
                "cannot load butterfly lease interoperability oracle")
        butterfly = importlib.util.module_from_spec(butterfly_spec)
        sys.modules[butterfly_spec.name] = butterfly
        butterfly_spec.loader.exec_module(butterfly)

        def expect_butterfly_rejection(action: Callable[[], object], name: str) -> None:
            try:
                action()
            except butterfly.EvidenceError:
                return
            raise EvidenceError(f"butterfly interoperability accepted: {name}")

        cross_root = root / "butterfly-cross-lease-root"
        cross_root.mkdir(mode=0o700)
        with SUPPORT.PairLease(cpu, sibling, root=cross_root) as identity:
            Path(identity["path"]).rename(cross_root / "old.lock")
            expect_butterfly_rejection(
                lambda: butterfly.PairLease(
                    sibling, cpu, root=cross_root).__enter__(),
                "butterfly file-replacement lease overlap")
        reverse_root = root / "butterfly-reverse-lease-root"
        reverse_root.mkdir(mode=0o700)
        with butterfly.PairLease(cpu, sibling, root=reverse_root):
            moved = root / "butterfly-reverse-lease-root-moved"
            reverse_root.rename(moved)
            reverse_root.mkdir(mode=0o700)
            expect_evidence_error(
                lambda: SUPPORT.PairLease(
                    sibling, cpu, root=reverse_root).__enter__(),
                "butterfly directory-replacement lease overlap")

        immutable_root = root / "immutable"
        immutable_root.mkdir(mode=0o700)
        immutable_control = immutable_root / "control"
        immutable_candidate = immutable_root / "candidate"
        for destination in (immutable_control, immutable_candidate):
            shutil.copy2("/bin/true", destination)
            os.chmod(destination, 0o700)
        immutable_initial = {
            "control_executable": SUPPORT.artifact_identity(
                immutable_control, "executable"),
            "candidate_executable": SUPPORT.artifact_identity(
                immutable_candidate, "executable"),
        }
        immutable_specification = {
            "control_executable": str(immutable_control),
            "candidate_executable": str(immutable_candidate),
        }
        immutable_guard = ImmutableExecutables(
            immutable_root, immutable_specification, immutable_initial)
        with immutable_guard as (immutable_execution, immutable_descriptors):
            shutil.copy2(mock, immutable_candidate)
            os.chmod(immutable_candidate, 0o700)
            returncode, stdout, stderr, _ = run_bounded(
                ("/usr/bin/taskset", "-c", str(cpu),
                 immutable_execution["candidate"]["command_path"]),
                CHILD_ENVIRONMENT, 5.0,
                (immutable_descriptors["candidate"],))
            require(returncode == 0 and stdout == b"" and stderr == b"",
                    "immutable executable descriptor followed a transient path replacement")

        expect_evidence_error(
            lambda: validate_exact_affinity(set(allowed) | {SUPPORT.MAX_CPU_ID}),
            "housekeeping affinity mutation")

        publish_root = root / "publish"
        publish_root.mkdir(mode=0o700)
        stage = publish_root / "stage"
        stage.mkdir(mode=0o700)
        final = publish_root / "final"
        publish_no_replace(stage, final)
        another = publish_root / "another"
        another.mkdir(mode=0o700)
        expect_evidence_error(
            lambda: publish_no_replace(another, final), "no-replace publication")

        loop_root = root / "loop"
        loop_root.mkdir(mode=0o700)
        (loop_root / "a").symlink_to("b")
        (loop_root / "b").symlink_to("a")
        loop_result = SUPPORT.run_process_bounded(
            (str(Path(sys.executable).resolve()), str(Path(__file__).resolve()),
             "verify", "--manifest", str(loop_root / "a")),
            environment=CHILD_ENVIRONMENT, timeout=10.0,
            max_stdout=65536, max_stderr=65536)
        require(loop_result.returncode == 1 and b"Traceback" not in loop_result.stderr,
                "symlink-loop input was not rejected cleanly")

        invalid_publish = copy.deepcopy(fast_manifest)
        invalid_publish["status"] = "policy_failed"
        invalid_publish = resign(invalid_publish)
        rejected_output = root / "rejected-publication"
        expect_evidence_error(
            lambda: publish_manifest(
                root / "fast", rejected_output, invalid_publish,
                check_current_inputs=False, allow_self_test=True),
            "invalid manifest publication preflight")
        require((root / "fast").is_dir() and not rejected_output.exists(),
                "invalid publication changed the staged/output paths")
        published_output = root / "published-fast"
        published_status = publish_manifest(
            root / "fast", published_output, fast_manifest,
            check_current_inputs=False, allow_self_test=True)
        require(published_status == 0 and not (root / "fast").exists() and
                published_output.is_dir(),
                "valid publication did not perform one final no-replace rename")
        _, replay_status = validate_manifest(
            fast_manifest, published_output,
            check_current_inputs=False, allow_self_test=True)
        require(replay_status == 0,
                "published self-test evidence did not replay")
    print("low-encode-copy runner self-test passed: full 24-cell mock ABBA and hardened adversarial gates")
    return 0


def run_production_smoke(options: argparse.Namespace) -> int:
    """Exercise real benchmark-v2 semantics without producing timing evidence."""
    backends = tuple(options.backend or ("scalar",))
    require(len(backends) == len(set(backends)),
            "production-smoke backends must be unique")
    require(isinstance(options.timeout, float) and math.isfinite(options.timeout) and
            0 < options.timeout <= 60.0,
            "production-smoke timeout must be in (0,60]")
    executables = {
        "control": options.control.absolute(),
        "candidate": options.candidate.absolute(),
    }
    for side, executable in executables.items():
        require_bounded_regular(
            executable, f"production-smoke {side} executable",
            MAX_BUILD_ARTIFACT_BYTES)
        require(os.access(executable, os.X_OK),
                f"production-smoke {side} executable is not executable")
    campaign = {"iterations": 1, "reuse": 1, "warmup": 0}
    for backend in backends:
        base = FIXED_CELLS[0]
        cell = Cell(
            identifier=f"production-smoke-{backend}", backend=backend,
            field=base.field, k=base.k, r=base.r,
            shard_bytes=base.shard_bytes, role="target", seed=base.seed,
            losses=base.losses)
        validated: dict[str, dict[str, Any]] = {}
        for side, executable in executables.items():
            returncode, stdout, stderr, _ = run_bounded(
                benchmark_arguments(executable, cell, campaign),
                CHILD_ENVIRONMENT, options.timeout)
            require(returncode == 0 and not stderr,
                    f"production-smoke {side}/{backend} failed or wrote stderr")
            parsed = parse_json_bytes(
                stdout, f"production-smoke {side}/{backend}", MAX_STDOUT_BYTES)
            validated[side] = validate_result(parsed, cell, campaign)
        require(validated["control"]["digests"] ==
                validated["candidate"]["digests"],
                f"production-smoke control/candidate original, parity, or repaired "
                f"digests differ for {backend}")
    print(
        "low-encode-copy production smoke passed real benchmark-v2 round-trip and "
        "A/B wire/recovery semantics for: " + ",".join(backends))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    run = commands.add_parser(
        "run", help="run the fixed isolated low-encode-copy A/B campaign")
    run.add_argument(
        "--control",
        default="/home/catid/leopard-wt-low-copy-baseline/build/"
                "low-copy-authoritative/bench_leopard2",
        type=Path)
    run.add_argument(
        "--candidate",
        default="/home/catid/leopard-wt-low-copy/build/"
                "low-copy-authoritative/bench_leopard2",
        type=Path)
    run.add_argument(
        "--control-archive",
        default="/home/catid/leopard-wt-low-copy-baseline/build/"
                "low-copy-authoritative/liblibleopard.a",
        type=Path)
    run.add_argument(
        "--candidate-archive",
        default="/home/catid/leopard-wt-low-copy/build/"
                "low-copy-authoritative/liblibleopard.a",
        type=Path)
    run.add_argument(
        "--control-build-dir",
        default="/home/catid/leopard-wt-low-copy-baseline/build/low-copy-authoritative",
        type=Path)
    run.add_argument(
        "--candidate-build-dir",
        default="/home/catid/leopard-wt-low-copy/build/low-copy-authoritative",
        type=Path)
    run.add_argument(
        "--control-source-root",
        default="/home/catid/leopard-wt-low-copy-baseline",
        type=Path)
    run.add_argument(
        "--candidate-source-root",
        default="/home/catid/leopard-wt-low-copy",
        type=Path)
    run.add_argument("--reservation-file", required=True, type=Path)
    run.add_argument("--output", required=True, type=Path)
    run.add_argument("--cpu", required=True, type=int)
    run.add_argument("--reserved-sibling", required=True, type=int)
    run.add_argument("--iterations", default=5, type=int)
    run.add_argument("--warmup", default=1, type=int)
    run.add_argument("--reuse", default=1, type=int)
    run.add_argument("--timeout", default=180.0, type=float)
    run.add_argument("--taskset", default="/usr/bin/taskset", type=Path)
    run.add_argument("--ldd", default="/usr/bin/ldd", type=Path)
    run.set_defaults(function=run_campaign)

    verify = commands.add_parser(
        "verify", help="verify a retained pass or policy-failure manifest")
    verify.add_argument("--manifest", required=True, type=Path)
    verify.add_argument(
        "--no-current-input-check", action="store_true",
        help="replay structure/files without reopening original source/build paths")
    verify.set_defaults(function=verify_campaign)

    verify_failure = commands.add_parser(
        "verify-failure", help="verify retained invalid-run diagnostics")
    verify_failure.add_argument("--failure", required=True, type=Path)
    verify_failure.set_defaults(function=verify_failed_campaign)

    self_test = commands.add_parser(
        "self-test", help="run fast mock/adversarial tests without real benchmarking")
    self_test.set_defaults(function=run_self_test)
    smoke = commands.add_parser(
        "production-smoke",
        help="run a non-authoritative correctness/schema smoke against real binaries")
    smoke.add_argument(
        "--control",
        default="/home/catid/leopard-wt-low-copy-baseline/build/"
                "low-copy-authoritative/bench_leopard2",
        type=Path)
    smoke.add_argument(
        "--candidate",
        default="/home/catid/leopard-wt-low-copy/build/"
                "low-copy-authoritative/bench_leopard2",
        type=Path)
    smoke.add_argument(
        "--backend", action="append", choices=("scalar", "ssse3", "avx2"),
        help="backend to smoke; repeat for more than one (default: scalar)")
    smoke.add_argument("--timeout", default=30.0, type=float)
    smoke.set_defaults(function=run_production_smoke)
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    try:
        options = parser().parse_args(arguments)
        return int(options.function(options))
    except SystemExit as error:
        return 0 if error.code == 0 else 1
    except EvidenceError as error:
        print(f"low-encode-copy evidence error: {error}", file=sys.stderr)
        return 1
    except (OSError, ValueError, TypeError, KeyError, AttributeError,
            OverflowError, RecursionError, RuntimeError, MemoryError,
            subprocess.SubprocessError,
            json.JSONDecodeError) as error:
        print(
            f"low-encode-copy evidence input error: {type(error).__name__}: {error}",
            file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
