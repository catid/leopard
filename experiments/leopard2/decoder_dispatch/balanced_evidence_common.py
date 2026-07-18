#!/usr/bin/env python3
"""Shared, fail-closed helpers for balanced forced-decoder evidence."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import stat
import subprocess
import sys
from typing import Any


MATRIX_SCHEMA = "leopard2-balanced-forced-matrix/v1"
RUN_SCHEMA = "leopard2-balanced-forced-abba/v2"
CHECKPOINT_SCHEMA = "leopard2-balanced-forced-checkpoint/v1"
RECORD_SCHEMA = "leopard2-balanced-forced-record/v1"
SUMMARY_SCHEMA = "leopard2-balanced-forced-summary/v2"
BENCHMARK_SCHEMA = "leopard2-benchmark-v2"
ROUND_ORDERS = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
MODES = ("generic", "materialized", "tiled")
BACKENDS = ("scalar", "ssse3", "avx2")
MODE_SELECTORS = {
    "generic": ("--force-generic",),
    "materialized": ("--force-specialized", "--force-materialized"),
    "tiled": ("--force-specialized", "--force-tiled"),
}
MODE_PARAMETERS = {
    "generic": {
        "force_generic_decode": True,
        "force_specialized_decode": False,
        "force_tiled_decode": False,
        "force_materialized_decode": False,
    },
    "materialized": {
        "force_generic_decode": False,
        "force_specialized_decode": True,
        "force_tiled_decode": False,
        "force_materialized_decode": True,
    },
    "tiled": {
        "force_generic_decode": False,
        "force_specialized_decode": True,
        "force_tiled_decode": True,
        "force_materialized_decode": False,
    },
}
CPU_FIELDS = ("user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal")
EMPTY_SHA256 = hashlib.sha256(b"").hexdigest()
RUNNER_RELATIVE = "experiments/leopard2/decoder_dispatch/run_balanced_abba.py"
ANALYZER_RELATIVE = "experiments/leopard2/decoder_dispatch/analyze_balanced.py"
COMMON_RELATIVE = "experiments/leopard2/decoder_dispatch/balanced_evidence_common.py"
LEASE_RELATIVE = "tools/leopard2_jerasure_compare.py"


class EvidenceError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise EvidenceError(f"value is not canonical JSON: {error}") from error


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    try:
        raw = path.read_bytes()
        value = json.loads(raw.decode("utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise EvidenceError(f"cannot read {label} {path}: {error}") from error
    require(isinstance(value, dict), f"{label} is not a JSON object")
    return value, raw


def require_keys(value: object, keys: set[str], label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    require(set(value) == keys,
            f"{label} keys changed: expected {sorted(keys)}, got {sorted(value)}")
    return value


def strict_int(value: object, label: str, minimum: int = 0) -> int:
    require(type(value) is int and value >= minimum,
            f"{label} is not an integer >= {minimum}")
    return value


def hex_digest(value: object, label: str, length: int = 64) -> str:
    require(isinstance(value, str) and len(value) == length and value == value.lower(),
            f"{label} is not lowercase {length}-digit hexadecimal")
    try:
        int(value, 16)
    except ValueError as error:
        raise EvidenceError(f"{label} is not hexadecimal") from error
    return value


def finite(value: object, label: str, minimum: float = 0.0) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{label} is not numeric")
    result = float(value)
    require(math.isfinite(result) and result >= minimum,
            f"{label} is not finite and >= {minimum}")
    return result


def checked_output(arguments: list[str]) -> str:
    try:
        result = subprocess.run(
            arguments, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise EvidenceError(f"command failed: {arguments!r}: {error}") from error
    require(not result.stderr, f"command emitted stderr: {arguments!r}")
    return result.stdout.strip()


def file_identity(path: Path, root: Path | None, label: str) -> dict[str, Any]:
    path = path.absolute()
    require(path.is_file(), f"missing {label}: {path}")
    descriptor = path.stat()
    relative: str | None = None
    if root is not None:
        try:
            relative = str(path.resolve().relative_to(root.resolve()))
        except ValueError:
            relative = None
    return {
        "path": str(path),
        "realpath": str(path.resolve()),
        "relative_path": relative,
        "sha256": sha256_file(path),
        "size": descriptor.st_size,
        "mode": stat.S_IMODE(descriptor.st_mode),
        "device": descriptor.st_dev,
        "inode": descriptor.st_ino,
        "mtime_ns": descriptor.st_mtime_ns,
    }


def git_identity(root: Path, expected_commit: str) -> dict[str, Any]:
    root = root.resolve()
    hex_digest(expected_commit, "expected commit", 40)
    head = checked_output(["git", "-C", str(root), "rev-parse", "HEAD"])
    tree = checked_output(["git", "-C", str(root), "rev-parse", "HEAD^{tree}"])
    expected_tree = checked_output(
        ["git", "-C", str(root), "rev-parse", expected_commit + "^{tree}"])
    status = checked_output(
        ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=all"])
    require(head == expected_commit, f"source HEAD {head} != expected {expected_commit}")
    require(tree == expected_tree, "source HEAD tree differs from expected commit tree")
    require(not status, f"source tree is dirty: {status!r}")
    return {
        "root": str(root), "head": head, "tree": tree,
        "expected_commit": expected_commit, "status": "clean",
        "status_sha256": EMPTY_SHA256,
    }


def find_build_root(binary: Path) -> Path:
    for candidate in (binary.parent, *binary.parents):
        if (candidate / "CMakeCache.txt").is_file():
            return candidate
    raise EvidenceError(f"cannot find CMakeCache.txt above {binary}")


def refresh_build(binary: Path) -> list[str]:
    root = find_build_root(binary)
    command = ["cmake", "--build", str(root), "--target", "bench_leopard2", "-j", "1"]
    try:
        completed = subprocess.run(
            command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise EvidenceError(f"benchmark build refresh failed: {error}") from error
    return command


def _unique_file(candidates: list[Path], label: str) -> Path:
    unique = sorted({item.resolve() for item in candidates if item.is_file()})
    require(len(unique) == 1,
            f"expected exactly one {label}, found {[str(item) for item in unique]}")
    return unique[0]


def build_identity(source_root: Path, binary_relative: str) -> dict[str, Any]:
    source_root = source_root.resolve()
    relative = Path(binary_relative)
    require(not relative.is_absolute() and ".." not in relative.parts,
            "benchmark binary path must be source-root-relative")
    binary = source_root / relative
    require(binary.is_file() and os.access(binary, os.X_OK),
            f"benchmark executable is absent or not executable: {binary}")
    build_root = find_build_root(binary)
    cache = build_root / "CMakeCache.txt"
    homes = [line.split("=", 1)[1] for line in cache.read_text(
        encoding="utf-8", errors="strict").splitlines()
             if line.startswith("CMAKE_HOME_DIRECTORY:INTERNAL=")]
    require(len(homes) == 1 and Path(homes[0]).resolve() == source_root,
            "CMake cache does not name the clean source root")
    benchmark_source = source_root / "bench/leopard2/benchmark.cpp"
    decoder_source = source_root / "leopard2.cpp"
    dispatch_source = source_root / "Leopard2Dispatch.h"
    benchmark_object = _unique_file([
        item for item in build_root.rglob("benchmark.cpp.o")
        if "bench_leopard2.dir" in str(item.parent)
    ] + [
        item for item in build_root.rglob("benchmark.cpp.obj")
        if "bench_leopard2.dir" in str(item.parent)
    ], "bench_leopard2 benchmark object")
    decoder_object = _unique_file([
        item for item in build_root.rglob("leopard2.cpp.o")
        if "leopard.dir" in str(item.parent) and "test_hooks" not in str(item)
    ] + [
        item for item in build_root.rglob("leopard2.cpp.obj")
        if "leopard.dir" in str(item.parent) and "test_hooks" not in str(item)
    ], "production leopard2 object")
    archive_names = ("libleopard.a", "leopard.lib")
    archive = _unique_file(
        [item for name in archive_names for item in build_root.rglob(name)],
        "production Leopard archive")
    graph_candidates = [build_root / "build.ninja", build_root / "Makefile"]
    graph = _unique_file(graph_candidates, "CMake build graph")
    identities = {
        "benchmark": file_identity(benchmark_source, source_root, "benchmark source"),
        "decoder": file_identity(decoder_source, source_root, "decoder source"),
        "dispatch": file_identity(dispatch_source, source_root, "dispatch source"),
    }
    objects = {
        "benchmark": file_identity(benchmark_object, source_root, "benchmark object"),
        "decoder": file_identity(decoder_object, source_root, "decoder object"),
    }
    archive_identity = file_identity(archive, source_root, "production archive")
    binary_identity = file_identity(binary, source_root, "benchmark executable")
    require(objects["benchmark"]["mtime_ns"] <= binary_identity["mtime_ns"],
            "benchmark executable predates its benchmark object")
    require(objects["decoder"]["mtime_ns"] <= archive_identity["mtime_ns"] <=
            binary_identity["mtime_ns"], "executable/archive/object timestamps are stale")
    return {
        "root": str(build_root),
        "cmake_home": str(source_root),
        "cache": file_identity(cache, source_root, "CMake cache"),
        "graph": file_identity(graph, source_root, "CMake build graph"),
        "sources": identities,
        "objects": objects,
        "archive": archive_identity,
        "binary": binary_identity,
    }


def _external_file_identity(path: Path, label: str) -> dict[str, Any]:
    return file_identity(path, None, label)


def stable_cpuinfo(cpus: list[int]) -> list[dict[str, str]]:
    """Retain CPU identity fields while excluding frequency/load telemetry."""
    wanted = {
        "processor", "vendor_id", "cpu family", "model", "model name",
        "stepping", "microcode", "flags", "bugs", "address sizes",
        "CPU implementer", "CPU architecture", "CPU variant", "CPU part",
        "CPU revision", "Features",
    }
    blocks = Path("/proc/cpuinfo").read_text(
        encoding="utf-8", errors="strict").strip().split("\n\n")
    by_processor: dict[int, dict[str, str]] = {}
    for block in blocks:
        record: dict[str, str] = {}
        for line in block.splitlines():
            if ":" not in line:
                continue
            key, value = (item.strip() for item in line.split(":", 1))
            if key in wanted:
                record[key] = value
        try:
            processor = int(record.get("processor", "-1"))
        except ValueError:
            continue
        by_processor[processor] = record
    result = []
    for cpu in cpus:
        require(cpu in by_processor, f"CPU {cpu} is absent from stable /proc/cpuinfo")
        result.append(by_processor[cpu])
    return result


def runtime_identity(binary: Path, cpu: int, sibling: int,
                     affinity: set[int]) -> dict[str, Any]:
    require(cpu in affinity and sibling in affinity and cpu != sibling,
            "runtime CPU pair is outside affinity or not distinct")
    records = []
    expected = sorted((cpu, sibling))
    for logical in expected:
        topology = Path(f"/sys/devices/system/cpu/cpu{logical}/topology")
        try:
            sibling_text = (topology / "thread_siblings_list").read_text(
                encoding="ascii").strip()
            core_id = (topology / "core_id").read_text(encoding="ascii").strip()
            package_id = (topology / "physical_package_id").read_text(
                encoding="ascii").strip()
        except OSError as error:
            raise EvidenceError(f"cannot read CPU {logical} topology: {error}") from error
        parsed: list[int] = []
        for part in sibling_text.split(","):
            bounds = part.split("-", 1)
            start = int(bounds[0])
            stop = int(bounds[-1])
            parsed.extend(range(start, stop + 1))
        require(sorted(set(parsed)) == expected,
                f"CPU {logical} is not in the exact SMT pair {expected}")
        records.append({
            "cpu": logical, "thread_siblings_list": sibling_text,
            "thread_siblings": sorted(set(parsed)), "core_id": core_id,
            "physical_package_id": package_id,
        })
    require(records[0]["core_id"] == records[1]["core_id"] and
            records[0]["physical_package_id"] == records[1]["physical_package_id"],
            "CPU pair does not share one physical core")
    try:
        ldd = subprocess.run(
            ["ldd", str(binary)], check=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise EvidenceError(f"cannot inspect benchmark runtime libraries: {error}") from error
    require(not ldd.stderr, "ldd emitted stderr")
    library_paths = []
    for line in ldd.stdout.splitlines():
        words = line.strip().split()
        candidate = None
        if len(words) >= 3 and words[1] == "=>" and words[2].startswith("/"):
            candidate = words[2]
        elif words and words[0].startswith("/"):
            candidate = words[0]
        if candidate:
            library_paths.append(Path(candidate))
    libraries = [_external_file_identity(path, "runtime library")
                 for path in sorted(set(library_paths))]
    normalized_ldd = re.sub(r"\s+\(0x[0-9a-fA-F]+\)", "", ldd.stdout)
    uname = platform.uname()
    environment = {key: os.environ.get(key) for key in (
        "LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT", "OMP_NUM_THREADS",
        "OMP_DYNAMIC", "OMP_PROC_BIND")}
    require(not environment["LD_PRELOAD"] and not environment["LD_AUDIT"],
            "LD_PRELOAD and LD_AUDIT must be unset for authoritative evidence")
    return {
        "platform": {
            "system": uname.system, "node": uname.node, "release": uname.release,
            "version": uname.version, "machine": uname.machine,
        },
        "python": {
            "version": sys.version, "implementation": platform.python_implementation(),
            "executable": _external_file_identity(
                Path(sys.executable), "Python executable"),
            "byteorder": sys.byteorder,
        },
        "affinity": sorted(affinity),
        "cpu": cpu,
        "sibling": sibling,
        "topology": records,
        "clock_ticks_per_second": os.sysconf("SC_CLK_TCK"),
        "cpuinfo": stable_cpuinfo(expected),
        "runtime_libraries": libraries,
        "ldd_normalized_sha256": sha256_bytes(normalized_ldd.encode("utf-8")),
        "environment": environment,
    }


def normalize_matrix(value: object) -> tuple[str, list[dict[str, Any]]]:
    matrix = require_keys(value, {"schema", "name", "cases"}, "matrix")
    require(matrix["schema"] == MATRIX_SCHEMA, "unexpected matrix schema")
    name = matrix["name"]
    require(isinstance(name, str) and name, "matrix name is empty")
    raw_cases = matrix["cases"]
    require(isinstance(raw_cases, list) and raw_cases, "matrix cases are empty")
    cases = []
    names = set()
    keys = {
        "name", "K", "R", "profile", "field", "backend", "shard_bytes",
        "loss_count", "batch", "reuse", "iterations", "warmup", "seed",
        "control_mode", "candidate_mode",
    }
    for index, raw in enumerate(raw_cases):
        case = require_keys(raw, keys, f"matrix case {index}")
        case_name = case["name"]
        require(isinstance(case_name, str) and case_name and case_name not in names,
                f"matrix case {index} name is empty or duplicated")
        names.add(case_name)
        for key in ("K", "R", "shard_bytes", "loss_count", "batch", "reuse",
                    "iterations", "seed"):
            strict_int(case[key], f"case {case_name} {key}", 1)
        strict_int(case["warmup"], f"case {case_name} warmup")
        require(case["profile"] == "legacy_high_v1" and case["field"] == "gf8",
                f"case {case_name} is outside balanced legacy-high GF8")
        k = case["K"]
        r = case["R"]
        require(k == r and 5 <= k <= 128 and case["loss_count"] == k,
                f"case {case_name} is not balanced full-original recovery")
        side = 1 << (r - 1).bit_length()
        parent = 1 << (k + side - 1).bit_length()
        require(side >= 8 and parent == 2 * side and parent <= 256,
                f"case {case_name} has an invalid balanced GF8 parent")
        require(case["backend"] in BACKENDS,
                f"case {case_name} uses an unsupported backend")
        control = case["control_mode"]
        candidate = case["candidate_mode"]
        require(control in MODES and candidate in MODES and control != candidate,
                f"case {case_name} has invalid or identical forced modes")
        normalized = dict(case)
        normalized["padded_side"] = side
        normalized["parent_count"] = parent
        cases.append(normalized)
    return name, cases


def role_mode(case: dict[str, Any], role: str) -> str:
    require(role in ("control", "candidate"), f"unknown role: {role}")
    return case[role + "_mode"]


def benchmark_command(binary: Path, case: dict[str, Any], output: Path,
                      cpu: int, role: str) -> list[str]:
    mode = role_mode(case, role)
    return [
        "/usr/bin/taskset", "-c", str(cpu),
        "/usr/bin/env", "OMP_NUM_THREADS=1", "OMP_DYNAMIC=false",
        "OMP_PROC_BIND=close", str(binary.absolute()),
        "--k", str(case["K"]), "--r", str(case["R"]), "--profile", "high",
        "--field", "gf8", "--backend", case["backend"],
        "--bytes", str(case["shard_bytes"]), "--loss", str(case["loss_count"]),
        "--batch", str(case["batch"]), "--reuse", str(case["reuse"]),
        "--iterations", str(case["iterations"]), "--warmup", str(case["warmup"]),
        "--threads", "1", "--seed", str(case["seed"]),
        *MODE_SELECTORS[mode], "--skip-legacy", "--retain-samples",
        "--json", str(output.absolute()),
    ]


def cpu_snapshot(cpu: int) -> dict[str, int]:
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            values = [int(item) for item in line.split()[1:1 + len(CPU_FIELDS)]]
            require(len(values) == len(CPU_FIELDS), f"CPU {cpu} stat row is incomplete")
            return dict(zip(CPU_FIELDS, values))
    raise EvidenceError(f"CPU {cpu} is absent from /proc/stat")


def cpu_delta(before: dict[str, int], after: dict[str, int]) -> dict[str, int]:
    require(set(before) == set(CPU_FIELDS) and set(after) == set(CPU_FIELDS),
            "CPU snapshot fields changed")
    result = {key: after[key] - before[key] for key in CPU_FIELDS}
    require(all(value >= 0 for value in result.values()), "CPU counters moved backwards")
    result["idle_total"] = result["idle"] + result["iowait"]
    result["nonidle_total"] = sum(
        result[key] for key in CPU_FIELDS if key not in {"idle", "iowait"})
    result["total"] = result["idle_total"] + result["nonidle_total"]
    return result


def isolation(before_cpu: dict[str, int], after_cpu: dict[str, int],
              before_sibling: dict[str, int], after_sibling: dict[str, int]) -> dict[str, Any]:
    timed = cpu_delta(before_cpu, after_cpu)
    sibling = cpu_delta(before_sibling, after_sibling)
    accepted = (timed["nonidle_total"] >= 1 and sibling["total"] >= 1 and
                sibling["nonidle_total"] == 0)
    return {
        "policy": {
            "source": "/proc/stat", "timed_min_nonidle_jiffies": 1,
            "sibling_min_total_jiffies": 1, "sibling_max_nonidle_jiffies": 0,
        },
        "accepted": accepted,
        "timed_before": before_cpu, "timed_after": after_cpu, "timed_delta": timed,
        "sibling_before": before_sibling, "sibling_after": after_sibling,
        "sibling_delta": sibling,
    }
