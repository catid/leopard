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
import shutil
import signal
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
RAW_SCHEMA = "leopard2-low-encode-copy-raw/v1"
MANIFEST_SCHEMA = "leopard2-low-encode-copy-manifest/v1"
FAILURE_SCHEMA = "leopard2-low-encode-copy-failure/v1"
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
HEX64 = re.compile(r"^[0-9a-f]{16}$")
HEX256 = re.compile(r"^[0-9a-f]{64}$")
TARGET_SPEEDUP = 1.05
NEIGHBOR_REGRESSION_FLOOR = 1.0 / 1.02


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    metadata = os.lstat(path)
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
            f"{what} is not a single-link regular file")
    require(metadata.st_size <= limit, f"{what} exceeds {limit} bytes")
    return parse_json_bytes(path.read_bytes(), what, limit)


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


def safe_evidence_path(root: Path, relative: object) -> Path:
    require(isinstance(relative, str) and relative and not os.path.isabs(relative),
            "evidence path is not relative")
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


def artifact_record(root: Path, path: Path) -> dict[str, Any]:
    path = path.resolve(strict=True)
    path.relative_to(root.resolve())
    metadata = os.lstat(path)
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
            f"retained artifact is not a single-link regular file: {path}")
    return {
        "path": str(path.relative_to(root.resolve())),
        "sha256": sha256_file(path),
        "size": metadata.st_size,
    }


def verify_artifact(root: Path, value: object, what: str) -> Path:
    require(isinstance(value, dict) and set(value) == {"path", "sha256", "size"},
            f"{what} identity is incomplete")
    require(isinstance(value.get("sha256"), str) and
            HEX256.fullmatch(value["sha256"]) is not None and
            isinstance(value.get("size"), int) and not isinstance(value["size"], bool) and
            0 <= value["size"] <= MAX_JSON_BYTES,
            f"{what} identity is invalid")
    path = safe_evidence_path(root, value.get("path"))
    metadata = os.lstat(path)
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
            metadata.st_size == value["size"] and
            sha256_file(path) == value["sha256"],
            f"{what} is missing or changed")
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
    require(cell.k > 0 and cell.r > 0 and 0 <= cell.losses <= min(cell.k, cell.r),
            f"invalid counts in {cell.identifier}")
    require(cell.shard_bytes in {64, 256, 1024, 4096, 16384, 65536, 262144, 1048576},
            f"invalid aligned size in {cell.identifier}")
    padded = ceil_power_of_two(cell.k)
    parent = ceil_power_of_two(padded + cell.r)
    require(parent <= 65536 and
            ((cell.field == "gf8" and parent <= 256) or
             (cell.field == "gf16" and parent > 256)),
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
    result = float(value)
    require(math.isfinite(result), f"{what} is not finite")
    require(result > 0 if positive else result >= 0,
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
        "validate_isolation", "validate_topology",
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
            "pair_lease_schema": SUPPORT.PAIR_LEASE_SCHEMA,
            "reservation_schema": SUPPORT.RESERVATION_SCHEMA,
            "sibling_max_nonidle_jiffies": SUPPORT.MAX_SIBLING_NONIDLE_JIFFIES,
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


def validate_build_provenance(
    side: str, specification: Mapping[str, Any]
) -> dict[str, Any]:
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
    return result


def input_snapshot(specification: Mapping[str, Any]) -> dict[str, Any]:
    support_contract = validate_support_contract()
    control_build = validate_build_provenance("control", specification)
    candidate_build = validate_build_provenance("candidate", specification)
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
        "LEO2_BUILD_TESTS", "LEO2_ENABLE_CUDA",
    ):
        require(control_cache.get(name) == candidate_cache.get(name),
                f"control/candidate CMake setting differs: {name}")
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
    require(isinstance(path, str) and Path(path).is_absolute() and
            (expected_path is None or path == str(Path(expected_path).resolve())),
            f"retained {what} path is invalid")
    require(isinstance(value.get("kind"), str) and value["kind"] and
            (expected_kind is None or value["kind"] == expected_kind),
            f"retained {what} kind is invalid")
    require(isinstance(value.get("size"), int) and not isinstance(value["size"], bool) and
            value["size"] >= 0 and
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


def validate_build_identity(
    value: object, side: str, specification: Mapping[str, Any],
    outer_executable: Mapping[str, Any], outer_archive: Mapping[str, Any],
) -> dict[str, Any]:
    expected_keys = {
        "archive_link_recipe", "archiver", "build_dir", "cmake_cache",
        "compile_commands", "compiler", "compiler_version_stdout",
        "executable_link_recipe", "validated_archive", "validated_archive_members",
        "validated_cache", "validated_compile_commands", "validated_executable",
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
        "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_GENERATOR": "Unix Makefiles",
        "ENABLE_OPENMP": "ON",
        "LEO2_BACKEND_VARIANT": "auto",
        "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_BUILD_TESTS": "OFF",
        "LEO2_ENABLE_CUDA": "OFF",
    }
    require(cache == expected_cache,
            f"retained {side} CMake production configuration is invalid")
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
            validate_file_identity(dependency.get("file"), f"{what} dependency")
    return value


def validate_retained_snapshot(
    value: object, specification: Mapping[str, Any]
) -> dict[str, Any]:
    expected_keys = {
        "candidate_archive", "candidate_build", "candidate_executable",
        "candidate_runtime_closure", "candidate_source", "control_archive",
        "control_build", "control_executable", "control_runtime_closure",
        "control_source", "ldd", "runner", "support", "taskset",
    }
    require(isinstance(value, dict) and set(value) == expected_keys,
            "retained source/build snapshot is incomplete")
    for name, spec_name, kind in (
        ("runner", "runner", "file"), ("taskset", "taskset", "executable"),
        ("ldd", "ldd", "executable"),
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
        parent = os.lstat(self.path.parent)
        require(stat.S_ISDIR(parent.st_mode) and parent.st_uid == os.getuid() and
                stat.S_IMODE(parent.st_mode) == 0o700,
                "CPU reservation parent must be an owned mode-0700 directory")
        flags = os.O_RDONLY
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        self.descriptor = os.open(self.path, flags)
        try:
            metadata = os.fstat(self.descriptor)
            path_metadata = os.lstat(self.path)
            require(stat.S_ISREG(metadata.st_mode) and metadata.st_uid == os.getuid() and
                    metadata.st_nlink == 1 and stat.S_IMODE(metadata.st_mode) == 0o600 and
                    (metadata.st_dev, metadata.st_ino) ==
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
            raise

    def validate_current(self) -> None:
        require(self.descriptor is not None and self.identity is not None,
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
        if self.descriptor is not None:
            fcntl.flock(self.descriptor, fcntl.LOCK_UN)
            os.close(self.descriptor)
            self.descriptor = None


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


def run_bounded(
    command: Sequence[str], environment: Mapping[str, str], timeout: float
) -> tuple[int, bytes, bytes, int]:
    started = time.monotonic_ns()
    process = subprocess.Popen(
        list(command), cwd="/", env=dict(environment), stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
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
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        returncode = process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
        returncode = process.returncode
        failure = "benchmark process did not terminate within five seconds after SIGKILL"
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


class XorShift64Star:
    def __init__(self, seed: int):
        self.state = seed & ((1 << 64) - 1)
        if self.state == 0:
            self.state = 0x9E3779B97F4A7C15

    def next(self) -> int:
        value = self.state
        value ^= value >> 12
        value ^= (value << 25) & ((1 << 64) - 1)
        value ^= value >> 27
        self.state = value & ((1 << 64) - 1)
        return (self.state * 2685821657736338717) & ((1 << 64) - 1)


def expected_losses(cell: Cell) -> list[int]:
    order = list(range(cell.k))
    random = XorShift64Star(cell.seed ^ 0xD1B54A32D192ED03)
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
        require(parameters.get(name) == expected_value,
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
            resolved.get("thread_count") == 1 and
            resolved.get("parent_count") == parent and
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
                contrasts.append(math.log(control / candidate))
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
    executable = Path(specification[f"{implementation}_executable"])
    command = [
        specification["taskset"], "-c", str(cpu),
        *benchmark_arguments(executable, cell, campaign),
    ]
    stem = f"{cell.identifier}/r{round_index:02d}-{slot:02d}-{label}-{implementation}"
    stdout_path = output / "children" / f"{stem}.stdout.json"
    stderr_path = output / "children" / f"{stem}.stderr.txt"
    try:
        returncode, stdout, stderr, duration_ns = run_bounded(
            command, CHILD_ENVIRONMENT, timeout)
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
        require(campaign.get("mode") == "self_test" and len(cells) == 1,
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
    require(value.get("rounds") == ROUNDS and
            value.get("sequence") == [list(item) for item in SEQUENCE],
            "campaign does not use five A1/B1/B2/A2 rounds")
    require(value.get("batch") == 1 and value.get("threads") == 1 and
            value.get("child_environment") == CHILD_ENVIRONMENT,
            "campaign is not a strict one-stripe, one-thread run")
    for name, minimum in (("iterations", 5), ("reuse", 1), ("warmup", 1)):
        item = value.get(name)
        require(isinstance(item, int) and not isinstance(item, bool) and item >= minimum,
                f"campaign {name} is invalid")
    timeout = value.get("timeout_seconds")
    require(isinstance(timeout, (int, float)) and not isinstance(timeout, bool) and
            math.isfinite(timeout) and timeout > 0,
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
                for item in allowed) and cpu in allowed and sibling in allowed and
            bool(set(allowed) - {cpu, sibling}),
            "campaign CPU pair, allowed set, or housekeeping set is invalid")
    return value, campaign_cells(value, allow_self_test)


def validate_input_specification(value: object) -> dict[str, str]:
    expected = {
        "candidate_archive", "candidate_build_dir", "candidate_executable",
        "candidate_source_root", "control_archive", "control_build_dir",
        "control_executable", "control_source_root", "ldd", "runner", "taskset",
    }
    require(isinstance(value, dict) and set(value) == expected and
            all(isinstance(item, str) and item for item in value.values()),
            "input specification is incomplete or has unexpected fields")
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
        require(observed == sequence,
                f"invocation order/relabel mismatch: {observed!r} != {sequence!r}")
        cell = cell_by_id[sequence[0]]
        implementation = sequence[4]
        expected_command = [
            specification["taskset"], "-c", str(cpu),
            *benchmark_arguments(
                Path(specification[f"{implementation}_executable"]), cell, campaign),
        ]
        require(invocation.get("command") == expected_command and
                invocation.get("environment") == CHILD_ENVIRONMENT and
                invocation.get("pinned_cpu") == cpu and
                invocation.get("returncode") == 0,
                "invocation execution contract was edited")
        duration = invocation.get("duration_ns")
        require(isinstance(duration, int) and not isinstance(duration, bool) and duration > 0,
                "invocation duration is invalid")
        initial_digest = sha256_bytes(canonical_bytes(initial))
        require(invocation.get("identity_before_digest") == initial_digest and
                invocation.get("identity_after_digest") == initial_digest and
                invocation.get("reservation_before") == reservation and
                invocation.get("reservation_after") == reservation,
                "invocation source/build/reservation identity was edited")
        stdout_path = verify_artifact(output, invocation.get("stdout"), "stdout")
        verify_artifact(output, invocation.get("stderr"), "stderr")
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
    invocations = validate_invocations(
        raw.get("invocations"), cells, campaign, specification, initial,
        reservation, output, complete=True)
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
    records: list[dict[str, Any]] = []
    for path in sorted(root.rglob("*")):
        metadata = os.lstat(path)
        if stat.S_ISDIR(metadata.st_mode):
            continue
        require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
                f"unexpected retained artifact type: {path}")
        records.append(artifact_record(output, path))
    return records


def validate_failure(
    failure: object, output: Path, allow_self_test: bool = False
) -> dict[str, Any]:
    failure = verify_signature(failure, "failure bundle")
    expected_keys = {
        "campaign", "created_utc", "digest", "error", "error_type", "host_initial",
        "identities_initial", "input_specification", "invocations", "isolation",
        "pair_lease", "reservation", "retained_files", "schema", "status", "traceback",
        "valid",
    }
    require(set(failure) == expected_keys and failure.get("schema") == FAILURE_SCHEMA and
            failure.get("status") == "failed" and failure.get("valid") is False,
            "failure bundle schema, status, or keys are invalid")
    validate_utc(failure.get("created_utc"), "failure created_utc")
    for name, maximum in (("error_type", 256), ("error", 8192), ("traceback", 65536)):
        item = failure.get(name)
        require(isinstance(item, str) and len(item.encode("utf-8")) <= maximum,
                f"failure {name} is invalid or oversized")
    campaign, cells = validate_campaign(failure.get("campaign"), allow_self_test)
    cpu = campaign["benchmark_cpu"]
    sibling = campaign["reserved_sibling"]
    specification = validate_input_specification(failure.get("input_specification"))
    host = failure.get("host_initial")
    if host is not None:
        SUPPORT.validate_host_record(
            host, cpu, sibling, campaign["allowed_cpu_set_at_launch"])
    pair_lease = failure.get("pair_lease")
    if pair_lease is not None:
        SUPPORT.validate_pair_lease_identity(pair_lease, cpu, sibling)
    reservation = failure.get("reservation")
    if reservation is not None:
        validate_reservation_identity(reservation, cpu, sibling)
    isolation = failure.get("isolation")
    if isolation is not None:
        SUPPORT.validate_isolation(isolation, cpu, sibling, require_accepted=False)
        require(pair_lease == isolation["pair_lease"],
                "failure isolation uses another pair lease")
    invocations = failure.get("invocations")
    initial = failure.get("identities_initial")
    if invocations:
        require(isinstance(initial, dict) and reservation is not None,
                "failed invocation prefix lacks source/build/reservation identity")
        validate_invocations(
            invocations, cells, campaign, specification, initial,
            reservation, output, complete=False)
    else:
        require(invocations == [], "failure invocation prefix is invalid")
        require(initial is None or isinstance(initial, dict),
                "failure initial identity is invalid")
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
        "sha256": sha256_file(path),
        "size": metadata.st_size,
    }


def validate_manifest(
    manifest: object,
    output: Path,
    check_current_inputs: bool,
    allow_self_test: bool = False,
) -> tuple[dict[str, Any], int]:
    manifest = verify_signature(manifest, "manifest")
    expected_keys = {
        "analysis", "campaign", "created_utc", "digest", "host", "identities",
        "isolation", "raw", "reservation", "schema", "status", "valid",
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
    return manifest, 0 if expected_status == "passed" else 2


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
    fsync_directory(output.parent)


def create_stage(output: Path) -> Path:
    output = output.resolve()
    require(not output.exists(), f"output path already exists: {output}")
    output.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    require(output.parent.is_dir(), "evidence output parent is not a directory")
    for _ in range(16):
        nonce = hashlib.sha256(os.urandom(32)).hexdigest()[:16]
        stage = output.parent / f".{output.name}.staging-{os.getpid()}-{nonce}"
        try:
            os.mkdir(stage, 0o700)
            return stage
        except FileExistsError:
            continue
    raise EvidenceError("cannot allocate a unique evidence staging directory")


def publish_manifest(stage: Path, output: Path, manifest: Mapping[str, Any]) -> int:
    _, status = validate_manifest(
        manifest, stage, check_current_inputs=True, allow_self_test=False)
    fsync_tree(stage)
    publish_no_replace(stage, output)
    retained = read_json_limited(output / "manifest.json", "published manifest")
    _, replay_status = validate_manifest(
        retained, output, check_current_inputs=True, allow_self_test=False)
    require(replay_status == status, "published manifest status changed")
    return status


def publish_failure(
    stage: Path, output: Path, failure: Mapping[str, Any]
) -> None:
    validate_failure(failure, stage)
    fsync_tree(stage)
    publish_no_replace(stage, output)
    retained = read_json_limited(output / "failure.json", "published failure")
    validate_failure(retained, output)


def authoritative_specification(options: argparse.Namespace) -> dict[str, str]:
    taskset = options.taskset.resolve(strict=True)
    ldd = options.ldd.resolve(strict=True)
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
        "ldd": str(ldd),
        "runner": str(Path(__file__).resolve(strict=True)),
        "taskset": str(taskset),
    }
    validate_input_specification(specification)
    return specification


def run_campaign(options: argparse.Namespace) -> int:
    require(options.iterations >= 5, "--iterations must be at least 5")
    require(options.warmup >= 1 and options.reuse >= 1,
            "--warmup and --reuse must be positive")
    require(math.isfinite(options.timeout) and options.timeout > 0,
            "--timeout must be positive and finite")
    require(options.cpu >= 0 and options.reserved_sibling >= 0 and
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
    invocations: list[dict[str, Any]] = []
    try:
        allowed_checked, housekeeping = SUPPORT.validate_topology(
            options.cpu, options.reserved_sibling)
        require(allowed_checked == allowed,
                "launch affinity changed during topology validation")
        host_initial = SUPPORT.host_identity(
            options.cpu, options.reserved_sibling, allowed)
        pair_guard = SUPPORT.PairLease(options.cpu, options.reserved_sibling)
        reservation_guard = HardenedReservation(
            options.reservation_file, options.cpu, options.reserved_sibling)
        with pair_guard as pair_lease, reservation_guard as reservation:
            os.sched_setaffinity(0, housekeeping)

            def validate_guards() -> None:
                pair_guard.validate_current()
                reservation_guard.validate_current()

            validate_guards()
            initial = input_snapshot(specification)

            def snapshot() -> Mapping[str, Any]:
                return input_snapshot(specification)

            before_monotonic_ns = time.monotonic_ns()
            before_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = SUPPORT.cpu_stat_snapshot(options.reserved_sibling)
            try:
                for cell in FIXED_CELLS:
                    for round_index in range(ROUNDS):
                        for slot, (label, implementation) in enumerate(SEQUENCE):
                            invocation = run_child(
                                implementation, label, cell, round_index, slot,
                                campaign, specification, initial, reservation, stage,
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
            result = publish_manifest(stage, output, manifest)
            stage = Path()
            print(output / "manifest.json")
            return result
    except Exception as error:
        if stage != Path() and stage.exists():
            for name in ("manifest.json", "raw.json", "failure.json"):
                path = stage / name
                if path.exists():
                    metadata = os.lstat(path)
                    require(stat.S_ISREG(metadata.st_mode),
                            f"cannot clean unsafe staged artifact {path}")
                    path.unlink()
            failure = signed({
                "campaign": campaign,
                "created_utc": utc_now(),
                "error": str(error)[:8192],
                "error_type": type(error).__name__[:256],
                "host_initial": host_initial,
                "identities_initial": initial,
                "input_specification": specification,
                "invocations": invocations,
                "isolation": isolation,
                "pair_lease": pair_lease,
                "reservation": reservation,
                "retained_files": retained_child_records(stage),
                "schema": FAILURE_SCHEMA,
                "status": "failed",
                "traceback": traceback.format_exc()[-65536:],
                "valid": False,
            })
            try:
                write_json_exclusive(stage / "failure.json", failure)
                publish_failure(stage, output, failure)
                stage = Path()
            except Exception:
                if stage.exists():
                    shutil.rmtree(stage)
                raise
        if isinstance(error, EvidenceError):
            raise
        raise EvidenceError(f"campaign failed: {type(error).__name__}: {error}") from error
    finally:
        try:
            os.sched_setaffinity(0, original_affinity)
        except OSError:
            pass
        if stage != Path() and stage.exists():
            shutil.rmtree(stage)


def verify_campaign(options: argparse.Namespace) -> int:
    manifest_path = options.manifest.resolve(strict=True)
    manifest = read_json_limited(manifest_path, "manifest")
    _, status = validate_manifest(
        manifest, manifest_path.parent,
        check_current_inputs=not options.no_current_input_check)
    if status == 0:
        print("low-encode-copy evidence verified and promotion policy passed")
    else:
        print("low-encode-copy evidence verified; promotion policy failed")
    return status


def verify_failed_campaign(options: argparse.Namespace) -> int:
    failure_path = options.failure.resolve(strict=True)
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
    cell: Cell, cpu: int, sibling: int, allowed: set[int]
) -> dict[str, Any]:
    return {
        "allowed_cpu_set_at_launch": sorted(allowed),
        "batch": 1,
        "benchmark_cpu": cpu,
        "cells": [asdict(cell)],
        "child_environment": dict(CHILD_ENVIRONMENT),
        "iterations": 5,
        "mode": "self_test",
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
) -> tuple[dict[str, Any], dict[str, Any]]:
    output = root
    output.mkdir(mode=0o700)
    cell = copy.deepcopy(FIXED_CELLS[0])
    campaign = self_test_campaign(cell, cpu, sibling, allowed)
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
        "ldd": "/usr/bin/ldd",
        "runner": str(Path(__file__).resolve()),
        "taskset": "/usr/bin/taskset",
    }
    initial = {"mock_identity": "stable"}
    reservation = synthetic_reservation(root.parent, cpu, sibling)
    invocations: list[dict[str, Any]] = []
    snapshot = lambda: initial
    validate_guards = lambda: None
    for round_index in range(ROUNDS):
        for slot, (label, implementation) in enumerate(SEQUENCE):
            invocations.append(run_child(
                implementation, label, cell, round_index, slot, campaign,
                specification, initial, reservation, output, cpu, 10.0,
                snapshot, validate_guards))
    pair_lease = synthetic_pair_lease(cpu, sibling)
    isolation = synthetic_isolation(
        cpu, sibling, pair_lease,
        sum(item["duration_ns"] for item in invocations) + 1)
    host = SUPPORT.host_identity(cpu, sibling, allowed)
    analysis = analyze(invocations, (cell,))
    raw = signed({
        "analysis": analysis,
        "campaign": campaign,
        "created_utc": "2000-01-01T00:00:00Z",
        "host_final": host,
        "host_initial": host,
        "identities_final": initial,
        "identities_initial": initial,
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
            TypeError, KeyError, AttributeError):
        return
    raise EvidenceError(f"adversarial mutation was accepted: {name}")


def resign(value: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(value))
    payload.pop("digest", None)
    return signed(payload)


def run_self_test(_options: argparse.Namespace) -> int:
    validate_support_contract()
    cpu, sibling, allowed = find_self_test_pair()
    with tempfile.TemporaryDirectory(prefix="leopard2-low-copy-selftest-") as temporary:
        root = Path(temporary)
        os.chmod(root, 0o700)
        mock = Path(__file__).with_name("mock_benchmark.py")
        control = root / "control-mock"
        fast = root / "candidate-fast-mock"
        slow = root / "candidate-slow-mock"
        for destination in (control, fast, slow):
            shutil.copy2(mock, destination)
            os.chmod(destination, 0o700)
        fast_raw, fast_manifest = build_self_test_bundle(
            root / "fast", fast, cpu, sibling, allowed)
        slow_raw, slow_manifest = build_self_test_bundle(
            root / "slow", slow, cpu, sibling, allowed)
        require(fast_raw["analysis"]["policy"]["passed"] is True and
                slow_raw["analysis"]["policy"]["passed"] is False,
                "mock pass/fail policy separation failed")
        _, fast_status = validate_manifest(
            fast_manifest, root / "fast", False, allow_self_test=True)
        _, slow_status = validate_manifest(
            slow_manifest, root / "slow", False, allow_self_test=True)
        require((fast_status, slow_status) == (0, 2),
                "pass/policy-failure exit semantics are not distinct")

        failed_root = root / "failed"
        failed_root.mkdir(mode=0o700)
        failure = signed({
            "campaign": fast_raw["campaign"],
            "created_utc": "2000-01-01T00:00:00Z",
            "error": "intentional self-test failure",
            "error_type": "EvidenceError",
            "host_initial": None,
            "identities_initial": None,
            "input_specification": fast_raw["input_specification"],
            "invocations": [],
            "isolation": None,
            "pair_lease": None,
            "reservation": None,
            "retained_files": [],
            "schema": FAILURE_SCHEMA,
            "status": "failed",
            "traceback": "",
            "valid": False,
        })
        write_json_exclusive(failed_root / "failure.json", failure)
        validate_failure(failure, failed_root, allow_self_test=True)

        mutations: list[tuple[str, Callable[[dict[str, Any]], None]]] = [
            ("command", lambda value: value["invocations"][0]["command"].__setitem__(
                2, str(sibling))),
            ("duration", lambda value: value["invocations"][0].__setitem__(
                "duration_ns", 0)),
            ("wire digest", lambda value: value["invocations"][0]["result"]
             ["workload_digests"].__setitem__("transmitted_parity", "f" * 16)),
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
        symlink = root / "artifact-link"
        symlink.symlink_to(root / "fast" / "raw.json")
        expect_evidence_error(
            lambda: safe_evidence_path(root, "artifact-link"),
            "retained artifact symlink")
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

        lease_root = root / "lease-root"
        lease_root.mkdir(mode=0o700)
        first = SUPPORT.PairLease(cpu, sibling, root=lease_root)
        second = SUPPORT.PairLease(sibling, cpu, root=lease_root)
        with first:
            expect_evidence_error(lambda: second.__enter__(), "pair lease overlap")

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
    print("low-encode-copy runner self-test passed: mock ABBA, policy exits, 19 adversarial gates")
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
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    try:
        options = parser().parse_args(arguments)
        return int(options.function(options))
    except EvidenceError as error:
        print(f"low-encode-copy evidence error: {error}", file=sys.stderr)
        return 1
    except (OSError, ValueError, TypeError, KeyError, AttributeError,
            json.JSONDecodeError) as error:
        print(
            f"low-encode-copy evidence input error: {type(error).__name__}: {error}",
            file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
