#!/usr/bin/env python3
"""Authenticate and summarize materialized-versus-tiled high decode evidence.

The current v2 schema compares explicit modes in one binary from one clean Git
tree.  Validation is fail-closed: exact argv, source/build closure, retained
artifact bytes, selector fields, raw samples, derived statistics/rates, workload
digests, CPU accounting, and the CPU-pair lease are all recomputed.

Historical v1 two-worktree manifests remain replayable under their original,
explicitly weaker contract.  A v1 artifact can never be relabelled as v2.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
from typing import Any


CURRENT_SCHEMA = "leopard2-tiling-followup/v2"
HISTORICAL_SCHEMA = "leopard2-tiling-followup/v1"
BENCHMARK_SCHEMA = "leopard2-benchmark-v2"
CURRENT_SUMMARY_SCHEMA = "leopard2-tiled-high-dispatch-evidence/v2"
HISTORICAL_SUMMARY_SCHEMA = "leopard2-tiled-high-dispatch-evidence/v1"
BENCHMARK_SEED = 424242
BACKENDS = {"scalar", "ssse3", "avx2"}
CPU_FIELDS = ("user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal")
ROUND_ORDERS = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
ROLE_MODE = {"control": "materialized", "candidate": "tiled"}
MODE_SELECTOR = {
    "materialized": "--force-materialized",
    "tiled": "--force-tiled",
}
EMPTY_SHA256 = hashlib.sha256(b"").hexdigest()
RUNNER_RELATIVE = "experiments/leopard2/decoder_dispatch/run_tiled_high_abba.py"
ANALYZER_RELATIVE = "experiments/leopard2/decoder_dispatch/analyze_tiled_high.py"
LEASE_RELATIVE = "tools/leopard2_jerasure_compare.py"
PAIR_LEASE_SCHEMA = "leopard2-cpu-pair-lease/v1"
KERNEL_LEASE_SCHEMA = "leopard2-kernel-cpu-pair-lease/v1"
PAIR_LEASE_LOCK = (
    "dual_linux_abstract_af_unix_bind_and_fcntl_flock_"
    "exclusive_nonblocking_pair_wide"
)
KERNEL_LEASE_MECHANISM = "linux-abstract-af-unix-bind-exclusive"
# Two-sided Student-t 95% critical values.
T95 = {
    1: 12.706204736, 2: 4.302652730, 3: 3.182446305,
    4: 2.776445105, 5: 2.570581836, 6: 2.446911851,
    7: 2.364624252, 8: 2.306004135, 9: 2.262157163,
    10: 2.228138852, 11: 2.200985160, 12: 2.178812830,
    13: 2.160368656, 14: 2.144786688, 15: 2.131449546,
}


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


def load_json_bytes(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    try:
        raw = path.read_bytes()
        value = json.loads(raw.decode("utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise EvidenceError(f"cannot read {label} {path}: {error}") from error
    require(isinstance(value, dict), f"{label} is not a JSON object")
    return value, raw


def require_keys(value: object, expected: set[str], label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    require(set(value) == expected,
            f"{label} keys changed: expected {sorted(expected)}, got {sorted(value)}")
    return value


def hex_digest(value: object, label: str, length: int = 64) -> str:
    require(isinstance(value, str) and len(value) == length,
            f"{label} is not a {length * 4}-bit hexadecimal digest")
    try:
        int(value, 16)
    except ValueError as error:
        raise EvidenceError(f"{label} is not hexadecimal") from error
    require(value == value.lower(), f"{label} is not lowercase canonical hexadecimal")
    return value


def git_oid(value: object, label: str) -> str:
    return hex_digest(value, label, 40)


def strict_int(value: object, label: str, minimum: int = 0) -> int:
    require(isinstance(value, int) and not isinstance(value, bool) and value >= minimum,
            f"{label} is not an integer >= {minimum}")
    return value


def finite(value: object, label: str, minimum: float = 0.0) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{label} is not numeric")
    result = float(value)
    require(math.isfinite(result) and result >= minimum,
            f"{label} is not finite and >= {minimum}")
    return result


def close_six(actual: object, expected: float, label: str) -> float:
    value = finite(actual, label)
    require(abs(value - expected) <= 2.1e-6,
            f"{label}={value!r}, recomputed {expected!r}")
    return value


def require_exact_scalar(actual: object, expected: object, label: str) -> None:
    if isinstance(expected, bool):
        require(type(actual) is bool and actual is expected,
                f"{label}={actual!r}, expected exact boolean {expected!r}")
        return
    if isinstance(expected, int):
        require(type(actual) is int and actual == expected,
                f"{label}={actual!r}, expected exact integer {expected!r}")
        return
    if expected is None:
        require(actual is None, f"{label}={actual!r}, expected null")
        return
    require(type(actual) is type(expected) and actual == expected,
            f"{label}={actual!r}, expected exact {expected!r}")


def normalize_case(raw: object) -> dict[str, Any]:
    require(isinstance(raw, list), "case is not an array")
    if len(raw) == 6:
        name, k, r, byte_count, warmup, reuse = raw
        loss, backend, batch = 8, "avx2", 1
    elif len(raw) == 9:
        name, k, r, byte_count, warmup, reuse, loss, backend, batch = raw
    else:
        raise EvidenceError("case must contain six historical or nine calibrated fields")
    require(isinstance(name, str) and name, "case name is empty")
    for label, value in (("K", k), ("R", r), ("bytes", byte_count),
                         ("warmup", warmup), ("reuse", reuse), ("loss", loss),
                         ("batch", batch)):
        strict_int(value, f"case {name} {label}", 1)
    require(isinstance(backend, str) and backend in BACKENDS,
            f"case {name} has an unsupported backend")
    require(k + r <= 256 and loss <= min(k, r), f"case {name} has invalid GF8 counts")
    padded_side = 1 << (r - 1).bit_length()
    parent_count = 1 << (k + padded_side - 1).bit_length()
    return {
        "name": name, "K": k, "R": r, "shard_bytes": byte_count,
        "warmup": warmup, "reuse": reuse, "loss": loss,
        "backend": backend, "batch": batch, "padded_side": padded_side,
        "parent_count": parent_count,
    }


def versioned_checkpoint_cases() -> list[list[Any]]:
    cases: list[list[Any]] = []
    for kib in (16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96):
        for loss in (1, 4, 8, 16, 32):
            cases.append([f"avx2-t32-b{kib}k-l{loss}", 224, 32, kib * 1024,
                          8, 8, loss, "avx2", 1])
    for kib in (24, 32, 40, 48, 56, 64, 80, 96, 112):
        for loss in (1, 8, 16):
            cases.append([f"ssse3-t32-b{kib}k-l{loss}", 224, 32, kib * 1024,
                          8, 8, loss, "ssse3", 1])
    for kib in (32, 48, 64, 80, 96):
        for loss in (1, 8, 16):
            cases.append([f"scalar-t32-b{kib}k-l{loss}", 224, 32, kib * 1024,
                          8, 8, loss, "scalar", 1])
    for kib in (32, 48, 64, 80):
        for loss in (1, 8, 16):
            cases.append([f"avx2-t32-b{kib}k-l{loss}-batch2", 224, 32,
                          kib * 1024, 4, 4, loss, "avx2", 2])
    cases.append(["avx2-t32-b64k-l8-reuse1", 224, 32, 64 * 1024,
                  8, 1, 8, "avx2", 1])
    cases.append(["avx2-t32-b64k-l8-reuse64", 224, 32, 64 * 1024,
                  8, 64, 8, "avx2", 1])
    for recovery in (16, 64):
        original = 256 - recovery
        for kib in (32, 64, 96):
            cases.append([f"avx2-t{recovery}-b{kib}k-l8", original, recovery,
                          kib * 1024, 8, 8, 8, "avx2", 1])
    return cases


def versioned_smoke_cases() -> list[list[Any]]:
    return [
        ["avx2-t32-b32k-l8", 224, 32, 32 * 1024, 2, 2, 8, "avx2", 1],
        ["avx2-t32-b96k-l8", 224, 32, 96 * 1024, 2, 2, 8, "avx2", 1],
    ]


def clustered_summary(round_logs: list[float]) -> dict[str, Any]:
    require(len(round_logs) >= 2, "at least two independent rounds are required")
    degrees = len(round_logs) - 1
    require(degrees in T95, "unsupported independent-round count")
    mean = statistics.fmean(round_logs)
    deviation = statistics.stdev(round_logs)
    half_width = T95[degrees] * deviation / math.sqrt(len(round_logs))
    lower = math.exp(mean - half_width)
    upper = math.exp(mean + half_width)
    speedup = math.exp(mean)
    return {
        "independent_round_count": len(round_logs),
        "degrees_of_freedom": degrees,
        "independent_round_log_contrasts": round_logs,
        "geometric_control_over_tiled": speedup,
        "improvement_percent": 100.0 * (speedup - 1.0),
        "ci95_lower_percent": 100.0 * (lower - 1.0),
        "ci95_upper_percent": 100.0 * (upper - 1.0),
        "credible_regression_over_2_percent": upper < 0.98,
        "credible_gain_at_least_5_percent": lower >= 1.05,
    }


def cpu_snapshot(value: object, label: str) -> dict[str, int]:
    result = require_keys(value, set(CPU_FIELDS), label)
    return {key: strict_int(result[key], f"{label}.{key}") for key in CPU_FIELDS}


def recompute_cpu_delta(before: dict[str, int], after: dict[str, int]) -> dict[str, int]:
    result = {key: after[key] - before[key] for key in CPU_FIELDS}
    require(all(value >= 0 for value in result.values()), "CPU counters moved backwards")
    result["idle_total"] = result["idle"] + result["iowait"]
    result["nonidle_total"] = sum(
        result[key] for key in CPU_FIELDS if key not in {"idle", "iowait"})
    result["total"] = result["idle_total"] + result["nonidle_total"]
    return result


def validate_cpu_triplet(before_value: object, after_value: object,
                         delta_value: object, label: str) -> dict[str, int]:
    before = cpu_snapshot(before_value, label + " before")
    after = cpu_snapshot(after_value, label + " after")
    expected = recompute_cpu_delta(before, after)
    delta = require_keys(delta_value, set(expected), label + " delta")
    for key, expected_value in expected.items():
        require(strict_int(delta[key], f"{label} delta.{key}") == expected_value,
                f"{label} delta.{key} is not recomputed")
    return expected


def validate_source(value: object, label: str) -> dict[str, Any]:
    source = require_keys(value, {
        "root", "head", "tree", "expected_commit", "status", "status_sha256"
    }, label)
    require(isinstance(source["root"], str) and Path(source["root"]).is_absolute(),
            f"{label} root is not absolute")
    require(git_oid(source["head"], label + " head") ==
            git_oid(source["expected_commit"], label + " expected commit"),
            f"{label} head/expected commit differ")
    git_oid(source["tree"], label + " tree")
    require(source["status"] == "clean" and source["status_sha256"] == EMPTY_SHA256,
            f"{label} is not a clean source snapshot")
    return source


def validate_file_identity(value: object, label: str,
                           expected_relative: str | None = None) -> dict[str, Any]:
    identity = require_keys(value, {
        "path", "realpath", "relative_path", "sha256", "size", "mode", "device", "inode"
    }, label)
    require(isinstance(identity["path"], str) and Path(identity["path"]).is_absolute(),
            f"{label} path is not absolute")
    require(isinstance(identity["realpath"], str) and Path(identity["realpath"]).is_absolute(),
            f"{label} realpath is not absolute")
    require(isinstance(identity["relative_path"], str) and identity["relative_path"] and
            not Path(identity["relative_path"]).is_absolute() and
            ".." not in Path(identity["relative_path"]).parts,
            f"{label} relative path is unsafe")
    if expected_relative is not None:
        require(identity["relative_path"] == expected_relative,
                f"{label} relative path changed")
    hex_digest(identity["sha256"], label + " SHA-256")
    strict_int(identity["size"], label + " size", 1)
    strict_int(identity["mode"], label + " mode", 1)
    strict_int(identity["device"], label + " device")
    strict_int(identity["inode"], label + " inode", 1)
    return identity


def validate_build(value: object, root: str, label: str) -> dict[str, Any]:
    build = require_keys(value, {"binary", "cache", "cmake_home"}, label)
    binary = validate_file_identity(build["binary"], label + " binary")
    cache = validate_file_identity(build["cache"], label + " cache")
    require(build["cmake_home"] == root, f"{label} CMake source root changed")
    require(cache["relative_path"] == str(Path(binary["relative_path"]).parent /
                                          "CMakeCache.txt"),
            f"{label} binary/cache directories differ")
    require(binary["path"] == str(Path(root) / binary["relative_path"]) and
            cache["path"] == str(Path(root) / cache["relative_path"]),
            f"{label} paths are not rooted at the source snapshot")
    return build


def checked_git(root: Path, arguments: list[str]) -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(root), *arguments], check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise EvidenceError(f"cannot validate live Git source {root}: {error}") from error
    require(not result.stderr, f"live Git validation emitted stderr: {result.stderr}")
    return result.stdout.strip()


def validate_live_identity(source: dict[str, Any], build: dict[str, Any],
                           runner: dict[str, Any], lease_source: dict[str, Any],
                           source_override: Path | None,
                           require_canonical_analyzer: bool = True) -> None:
    stored_root = Path(source["root"])
    root = stored_root if source_override is None else source_override.resolve()
    require(root.is_dir(), f"live source root is absent: {root}")
    require(checked_git(root, ["rev-parse", "HEAD"]) == source["head"],
            "live source HEAD differs from evidence")
    require(checked_git(root, ["rev-parse", "HEAD^{tree}"]) == source["tree"],
            "live source tree differs from evidence")
    require(not checked_git(root, ["status", "--porcelain", "--untracked-files=all"]),
            "live source is dirty")
    if require_canonical_analyzer:
        require(Path(__file__).resolve() == (root / ANALYZER_RELATIVE).resolve(),
                "current evidence must be analyzed by the canonical clean-tree analyzer")

    def live_file(identity: dict[str, Any], executable: bool = False) -> Path:
        path = root / identity["relative_path"]
        require(path.is_file(), f"live identity file is absent: {path}")
        require(path.stat().st_size == identity["size"] and
                sha256_file(path) == identity["sha256"],
                f"live identity bytes differ: {path}")
        require((path.stat().st_mode & 0o7777) == identity["mode"],
                f"live identity mode differs: {path}")
        if executable:
            require(os.access(path, os.X_OK), f"live benchmark is not executable: {path}")
        if source_override is None:
            require(str(path) == identity["path"] and str(path.resolve()) == identity["realpath"] and
                    path.stat().st_dev == identity["device"] and
                    path.stat().st_ino == identity["inode"],
                    f"same-root live identity metadata differs: {path}")
        return path

    live_file(build["binary"], executable=True)
    cache = live_file(build["cache"])
    prefix = "CMAKE_HOME_DIRECTORY:INTERNAL="
    homes = [line[len(prefix):] for line in cache.read_text(encoding="utf-8").splitlines()
             if line.startswith(prefix)]
    require(len(homes) == 1 and Path(homes[0]).resolve() == stored_root,
            "retained CMake cache does not name its recorded source root")
    live_file(runner)
    live_file(lease_source)


def validate_lease(value: object, cpu: int, sibling: int) -> dict[str, Any]:
    lease = require_keys(value, {
        "lock", "path", "payload", "sha256", "device", "inode",
        "directory_device", "directory_inode", "kernel_lease"
    }, "pair lease")
    require(lease["lock"] == PAIR_LEASE_LOCK, "pair lease mechanism changed")
    payload = require_keys(lease["payload"], {"cpus", "schema", "uid"},
                           "pair lease payload")
    require(payload["schema"] == PAIR_LEASE_SCHEMA and
            payload["cpus"] == sorted((cpu, sibling)),
            "pair lease payload does not bind the CPU pair")
    uid = strict_int(payload["uid"], "pair lease UID")
    require(lease["sha256"] == canonical_sha256(payload),
            "pair lease payload SHA-256 mismatch")
    expected_path = (Path("/run/user") / str(uid) / "leopard2-cpu-leases" /
                     f"leopard2-cpu-pair-{uid}-{min(cpu, sibling)}-{max(cpu, sibling)}.lock")
    require(lease["path"] == str(expected_path), "pair lease path/payload mismatch")
    for key in ("device", "inode", "directory_device", "directory_inode"):
        strict_int(lease[key], "pair lease " + key, 1)
    kernel = require_keys(lease["kernel_lease"], {"schema", "mechanism", "name_sha256"},
                          "kernel pair lease")
    require(kernel["schema"] == KERNEL_LEASE_SCHEMA and
            kernel["mechanism"] == KERNEL_LEASE_MECHANISM,
            "kernel pair lease mechanism changed")
    hex_digest(kernel["name_sha256"], "kernel pair lease name SHA-256")
    return lease


def validate_stat_summary(value: object, label: str, iterations: int,
                          execution: bool, rates: dict[str, float] | None = None) -> dict[str, Any]:
    suffix = "_us_per_batch_call" if execution else "_us"
    sample_key = "samples" + suffix
    expected_keys = {
        "median" + suffix, "mad" + suffix, "minimum" + suffix,
        "maximum" + suffix, sample_key,
    }
    if rates:
        expected_keys.update(rates)
    summary = require_keys(value, expected_keys, label)
    samples = summary[sample_key]
    require(isinstance(samples, list) and len(samples) == iterations,
            f"{label} does not retain exactly {iterations} samples")
    numeric = [finite(item, f"{label} sample", 1e-300) for item in samples]
    median = statistics.median(numeric)
    mad = statistics.median(abs(item - median) for item in numeric)
    close_six(summary["median" + suffix], median, label + " median")
    close_six(summary["mad" + suffix], mad, label + " MAD")
    close_six(summary["minimum" + suffix], min(numeric), label + " minimum")
    close_six(summary["maximum" + suffix], max(numeric), label + " maximum")
    if rates:
        for name, byte_count in rates.items():
            close_six(summary[name], byte_count / (median * 1000.0), label + " " + name)
    return {"median": median, "mad": mad, "minimum": min(numeric),
            "maximum": max(numeric), "samples": numeric}


def validate_digests(value: object, label: str) -> dict[str, str]:
    digests = require_keys(value, {
        "algorithm", "original_data", "transmitted_parity", "recovered_originals"
    }, label)
    require(digests["algorithm"] == "fnv1a64", f"{label} algorithm changed")
    for key in ("original_data", "transmitted_parity", "recovered_originals"):
        hex_digest(digests[key], f"{label}.{key}", 16)
    return digests


def validate_inactive_selector_pair(parameters: dict[str, Any], label: str) -> None:
    tiled_present = "force_tiled_decode" in parameters
    materialized_present = "force_materialized_decode" in parameters
    require(tiled_present == materialized_present,
            f"{label}: selector fields are only partially present")
    if not tiled_present:
        return
    tiled = parameters["force_tiled_decode"]
    materialized = parameters["force_materialized_decode"]
    require(isinstance(tiled, bool) and isinstance(materialized, bool),
            f"{label}: selector fields are not boolean")
    require(not tiled and not materialized,
            f"{label}: historical/mixed selector fields are active")


def expected_command(binary: str, case: dict[str, Any], json_path: str,
                     cpu: int, iterations: int, variant: str) -> list[str]:
    selector = MODE_SELECTOR[ROLE_MODE[variant]]
    return [
        "/usr/bin/taskset", "-c", str(cpu),
        "/usr/bin/env", "OMP_NUM_THREADS=1", "OMP_DYNAMIC=false",
        "OMP_PROC_BIND=close", binary,
        "--k", str(case["K"]), "--r", str(case["R"]), "--profile", "high",
        "--field", "gf8", "--backend", case["backend"],
        "--bytes", str(case["shard_bytes"]), "--loss", str(case["loss"]),
        "--batch", str(case["batch"]), "--reuse", str(case["reuse"]),
        "--iterations", str(iterations), "--warmup", str(case["warmup"]),
        "--threads", "1", "--seed", str(BENCHMARK_SEED),
        "--force-specialized", selector, "--skip-legacy", "--retain-samples",
        "--json", json_path,
    ]


def validate_current_result(path: Path, case: dict[str, Any], variant: str,
                            iterations: int, expected_sha: str,
                            expected_digest_sha: str) -> dict[str, Any]:
    document, raw = load_json_bytes(path, "benchmark result")
    require(sha256_bytes(raw) == expected_sha, f"benchmark result hash changed: {path}")
    require_keys(document, {
        "schema", "build", "parameters", "resolved", "correctness",
        "workload_digests", "memory", "metrics", "legacy"
    }, f"benchmark result {path}")
    require(document["schema"] == BENCHMARK_SCHEMA, f"unexpected benchmark schema in {path}")
    require_keys(document["build"], {"compiler", "compiler_version", "cplusplus"},
                 f"{path} build")
    parameters = require_keys(document["parameters"], {
        "K", "R", "requested_profile", "requested_field", "requested_backend",
        "force_generic_decode", "force_specialized_decode", "force_tiled_decode",
        "force_materialized_decode", "skip_legacy", "retain_samples", "shard_bytes",
        "loss_count", "missing_original_indices", "batch", "reuse", "iterations",
        "warmup", "thread_count", "seed"
    }, f"{path} parameters")
    expected = {
        "K": case["K"], "R": case["R"], "requested_profile": "legacy_high_v1",
        "requested_field": "gf8", "requested_backend": case["backend"],
        "force_generic_decode": False, "force_specialized_decode": True,
        "force_tiled_decode": variant == "candidate",
        "force_materialized_decode": variant == "control",
        "skip_legacy": True, "retain_samples": True,
        "shard_bytes": case["shard_bytes"], "loss_count": case["loss"],
        "batch": case["batch"], "reuse": case["reuse"], "iterations": iterations,
        "warmup": case["warmup"], "thread_count": 1, "seed": BENCHMARK_SEED,
    }
    for key, expected_value in expected.items():
        require_exact_scalar(parameters[key], expected_value, f"{path}: {key}")
    missing = parameters["missing_original_indices"]
    require(isinstance(missing, list) and len(missing) == case["loss"] and
            all(isinstance(item, int) and not isinstance(item, bool) and
                0 <= item < case["K"] for item in missing) and
            missing == sorted(set(missing)), f"{path}: invalid missing-original set")
    resolved = require_keys(document["resolved"], {
        "profile", "field", "backend", "thread_count", "parent_count", "padded_side"
    }, f"{path} resolved")
    expected_resolved = {
        "profile": "legacy_high_v1", "field": "gf8", "backend": case["backend"],
        "thread_count": 1, "parent_count": case["parent_count"],
        "padded_side": case["padded_side"],
    }
    for key, expected_value in expected_resolved.items():
        require_exact_scalar(resolved[key], expected_value, f"{path}: resolved.{key}")
    correctness = require_keys(document["correctness"], {
        "leopard2_round_trip", "legacy_comparison"
    }, f"{path} correctness")
    require_exact_scalar(correctness["leopard2_round_trip"], True,
                         f"{path}: leopard2_round_trip")
    require_exact_scalar(correctness["legacy_comparison"], None,
                         f"{path}: legacy_comparison")
    digests = validate_digests(document["workload_digests"], f"{path} digests")
    require(canonical_sha256(digests) == expected_digest_sha,
            f"{path}: workload digest snapshot changed")
    memory = require_keys(document["memory"], {
        "scratch_alignment", "encode_scratch_bytes_per_stripe",
        "decode_scratch_bytes_per_stripe", "encode_scratch_bytes_batch",
        "decode_scratch_bytes_batch"
    }, f"{path} memory")
    for key, value in memory.items():
        strict_int(value, f"{path} memory.{key}", 1)
    require(memory["encode_scratch_bytes_batch"] ==
            memory["encode_scratch_bytes_per_stripe"] * case["batch"] and
            memory["decode_scratch_bytes_batch"] ==
            memory["decode_scratch_bytes_per_stripe"] * case["batch"],
            f"{path}: batch scratch is not derived from per-stripe scratch")

    metrics = require_keys(document["metrics"], {
        "codec_setup", "encode_execution", "decode_plan_setup", "decode_execution",
        "decode_amortized_at_reuse", "rate_semantics"
    }, f"{path} metrics")
    encode_input = case["K"] * case["shard_bytes"] * case["batch"]
    encode_output = case["R"] * case["shard_bytes"] * case["batch"]
    decode_input = (case["K"] - case["loss"] + case["R"]) * \
        case["shard_bytes"] * case["batch"]
    decode_output = case["loss"] * case["shard_bytes"] * case["batch"]
    validate_stat_summary(metrics["codec_setup"], f"{path} codec setup",
                          iterations, False)
    validate_stat_summary(metrics["encode_execution"], f"{path} encode execution",
                          iterations, True,
                          {"input_GB_per_s": encode_input,
                           "parity_output_GB_per_s": encode_output})
    setup = validate_stat_summary(metrics["decode_plan_setup"], f"{path} plan setup",
                                  iterations, False)
    execution = validate_stat_summary(
        metrics["decode_execution"], f"{path} decode execution", iterations, True,
        {"offered_received_GB_per_s": decode_input,
         "repaired_output_GB_per_s": decode_output})
    amortized = require_keys(metrics["decode_amortized_at_reuse"], {
        "reuse_count", "derived_median_us_per_batch_call",
        "offered_received_GB_per_s", "repaired_output_GB_per_s"
    }, f"{path} amortized decode")
    require(amortized["reuse_count"] == case["reuse"], f"{path}: reuse changed")
    amortized_us = execution["median"] + setup["median"] / case["reuse"]
    close_six(amortized["derived_median_us_per_batch_call"], amortized_us,
              f"{path} amortized median")
    close_six(amortized["offered_received_GB_per_s"],
              decode_input / (amortized_us * 1000.0), f"{path} amortized input rate")
    close_six(amortized["repaired_output_GB_per_s"],
              decode_output / (amortized_us * 1000.0), f"{path} amortized output rate")
    require(metrics["rate_semantics"] ==
            "offered_received counts all non-null shard pointers supplied; "
            "a plan may read a deterministic subset", f"{path}: rate semantics changed")
    legacy = require_keys(document["legacy"], {
        "available", "unavailable_reason", "codec_setup", "decode_timing_includes_setup",
        "encode_execution", "decode_including_setup"
    }, f"{path} legacy")
    expected_legacy = {
        "available": False, "unavailable_reason": "disabled by --skip-legacy",
        "codec_setup": None, "decode_timing_includes_setup": True,
        "encode_execution": None, "decode_including_setup": None,
    }
    for key, expected_value in expected_legacy.items():
        require_exact_scalar(legacy[key], expected_value, f"{path}: legacy.{key}")
    return {
        "variant": variant, "decode_us": execution["median"],
        "decode_mad_us": execution["mad"], "setup_us": setup["median"],
        "scratch_bytes": memory["decode_scratch_bytes_per_stripe"],
        "digests": digests, "raw_sha256": expected_sha,
    }


def validate_artifact(value: object, expected_name: str, path: Path,
                      label: str, empty: bool = False) -> str:
    artifact = require_keys(value, {"name", "size", "sha256"}, label)
    require(artifact["name"] == expected_name, f"{label} name changed")
    size = strict_int(artifact["size"], label + " size")
    digest = hex_digest(artifact["sha256"], label + " SHA-256")
    require(path.is_file() and path.stat().st_size == size and sha256_file(path) == digest,
            f"{label} retained bytes changed")
    if empty:
        require(size == 0 and digest == EMPTY_SHA256, f"{label} is not empty")
    else:
        require(size > 0, f"{label} is empty")
    return digest


def analyze_current(manifest_path: Path, require_rounds: int,
                    source_override: Path | None = None,
                    validate_live: bool = True,
                    allow_self_test: bool = False) -> dict[str, Any]:
    manifest, manifest_raw = load_json_bytes(manifest_path, "manifest")
    require_keys(manifest, {
        "schema", "status", "valid", "source", "build", "runner", "lease_source",
        "campaign", "records", "final", "content_sha256"
    }, "current manifest")
    require(manifest["schema"] == CURRENT_SCHEMA and manifest["status"] == "complete" and
            manifest["valid"] is True, "current manifest schema/status/validity is wrong")
    content_digest = hex_digest(manifest["content_sha256"], "manifest content SHA-256")
    unsigned = dict(manifest)
    del unsigned["content_sha256"]
    require(canonical_sha256(unsigned) == content_digest,
            "manifest canonical content SHA-256 mismatch")
    source = validate_source(manifest["source"], "source")
    build = validate_build(manifest["build"], source["root"], "build")
    runner = validate_file_identity(manifest["runner"], "runner", RUNNER_RELATIVE)
    lease_source = validate_file_identity(
        manifest["lease_source"], "lease source", LEASE_RELATIVE)
    require(runner["path"] == str(Path(source["root"]) / RUNNER_RELATIVE) and
            lease_source["path"] == str(Path(source["root"]) / LEASE_RELATIVE),
            "runner/lease source paths are not rooted at the source snapshot")
    final = require_keys(manifest["final"], {"source", "build", "runner", "lease_source"},
                         "final identity")
    require(final == {"source": source, "build": build, "runner": runner,
                      "lease_source": lease_source},
            "source/build/runner/lease identity changed during collection")
    if validate_live:
        validate_live_identity(
            source, build, runner, lease_source, source_override,
            require_canonical_analyzer=not allow_self_test)

    campaign = require_keys(manifest["campaign"], {
        "cpu", "sibling", "original_affinity", "housekeeping_affinity", "matrix",
        "iterations", "seed", "cases", "round_orders", "lease",
        "lease_identity_sha256", "start_time_ns", "end_time_ns", "cpu_before",
        "cpu_after", "cpu_delta", "sibling_before", "sibling_after", "sibling_delta"
    }, "campaign")
    cpu = strict_int(campaign["cpu"], "campaign CPU")
    sibling = strict_int(campaign["sibling"], "campaign sibling")
    require(cpu != sibling, "campaign CPU and sibling are identical")
    for key in ("original_affinity", "housekeeping_affinity"):
        affinity = campaign[key]
        require(isinstance(affinity, list) and affinity and
                all(isinstance(item, int) and not isinstance(item, bool) and item >= 0
                    for item in affinity) and affinity == sorted(set(affinity)),
                f"campaign {key} is invalid")
    require(cpu in campaign["original_affinity"] and sibling in campaign["original_affinity"] and
            cpu not in campaign["housekeeping_affinity"] and
            sibling not in campaign["housekeeping_affinity"] and
            campaign["housekeeping_affinity"] == sorted(
                set(campaign["original_affinity"]) - {cpu, sibling}),
            "campaign housekeeping affinity does not exclude exactly the timed pair")
    matrix = campaign["matrix"]
    require(matrix in {"checkpoint", "smoke"} or
            (allow_self_test and matrix == "self-test"), "unknown campaign matrix")
    iterations = strict_int(campaign["iterations"], "campaign iterations", 3)
    require(campaign["seed"] == BENCHMARK_SEED, "campaign seed changed")
    start_ns = strict_int(campaign["start_time_ns"], "campaign start time", 1)
    end_ns = strict_int(campaign["end_time_ns"], "campaign end time", 1)
    require(end_ns > start_ns, "campaign timestamps are not increasing")
    campaign_cpu_delta = validate_cpu_triplet(
        campaign["cpu_before"], campaign["cpu_after"], campaign["cpu_delta"],
        "campaign CPU")
    campaign_sibling_delta = validate_cpu_triplet(
        campaign["sibling_before"], campaign["sibling_after"],
        campaign["sibling_delta"], "campaign sibling")
    require(campaign_cpu_delta["nonidle_total"] > 0,
            "timed CPU recorded no non-idle work")
    lease = validate_lease(campaign["lease"], cpu, sibling)
    lease_digest = hex_digest(campaign["lease_identity_sha256"], "lease identity SHA-256")
    require(lease_digest == canonical_sha256(lease), "lease identity digest mismatch")

    if matrix == "checkpoint":
        require(campaign["cases"] == versioned_checkpoint_cases(),
                "checkpoint case matrix changed, reordered, or was truncated")
    elif matrix == "smoke":
        require(campaign["cases"] == versioned_smoke_cases(),
                "smoke case matrix changed, reordered, or was truncated")
    cases = [normalize_case(value) for value in campaign["cases"]]
    require(cases and len({case["name"] for case in cases}) == len(cases),
            "case names are empty or duplicated")
    case_by_name = {case["name"]: case for case in cases}
    orders = campaign["round_orders"]
    require(require_rounds == len(ROUND_ORDERS),
            "current schema requires all three versioned independent rounds")
    require(isinstance(orders, list) and len(orders) == require_rounds,
            "manifest must contain exactly three independent rounds")
    require(orders == [list(order) for order in ROUND_ORDERS],
            "current round schedule is not the versioned ABBA/BAAB sequence")

    records = manifest["records"]
    require(isinstance(records, list), "manifest records are missing")
    expected_count = len(cases) * len(orders) * 4
    require(len(records) == expected_count,
            f"manifest has {len(records)} records, expected {expected_count}")
    raw_root = manifest_path.parent / "raw"
    grouped: dict[str, dict[int, dict[int, dict[str, Any]]]] = {}
    raw_hashes = []
    seen = set()
    sequence = 0
    for case in cases:
        for round_index, order in enumerate(orders):
            for slot, expected_variant in enumerate(order):
                record = records[sequence]
                require_keys(record, {
                    "sequence", "case", "round", "slot", "variant", "mode", "selector",
                    "command", "command_sha256", "benchmark_argv0", "json_path_at_run",
                    "iterations", "seed", "lease_identity_sha256", "start_time_ns",
                    "end_time_ns", "cpu_before", "cpu_after", "cpu_delta",
                    "sibling_before", "sibling_after", "sibling_delta", "artifacts",
                    "workload_digests_sha256"
                }, f"record {sequence}")
                record_sequence = strict_int(record["sequence"],
                                             f"record {sequence} sequence")
                record_round = strict_int(record["round"], f"record {sequence} round")
                record_slot = strict_int(record["slot"], f"record {sequence} slot")
                require(record_sequence == sequence and record["case"] == case["name"] and
                        record_round == round_index and record_slot == slot,
                        f"record {sequence} coordinate/relabel mismatch")
                variant = record["variant"]
                require(variant == expected_variant, f"record {sequence} variant/order mismatch")
                mode = ROLE_MODE[variant]
                require(record["mode"] == mode and record["selector"] == MODE_SELECTOR[mode],
                        f"record {sequence} role/mode/selector mapping changed")
                require(record["benchmark_argv0"] == build["binary"]["path"],
                        f"record {sequence} benchmark argv[0] differs from bound binary")
                require(strict_int(record["iterations"], f"record {sequence} iterations", 3) ==
                        iterations and
                        strict_int(record["seed"], f"record {sequence} seed") == BENCHMARK_SEED,
                        f"record {sequence} iteration/seed snapshot changed")
                require(record["lease_identity_sha256"] == lease_digest,
                        f"record {sequence} lease identity changed")
                json_at_run = record["json_path_at_run"]
                require(isinstance(json_at_run, str) and Path(json_at_run).is_absolute() and
                        Path(json_at_run).parent.name == "raw",
                        f"record {sequence} JSON path-at-run is invalid")
                expected = expected_command(
                    build["binary"]["path"], case, json_at_run, cpu, iterations, variant)
                require(record["command"] == expected,
                        f"record {sequence} exact argv mismatch")
                require(record["command_sha256"] == canonical_sha256(expected),
                        f"record {sequence} argv SHA-256 mismatch")
                record_start = strict_int(record["start_time_ns"],
                                          f"record {sequence} start", 1)
                record_end = strict_int(record["end_time_ns"], f"record {sequence} end", 1)
                require(start_ns <= record_start < record_end <= end_ns,
                        f"record {sequence} timestamps escape the campaign")
                record_cpu = validate_cpu_triplet(
                    record["cpu_before"], record["cpu_after"], record["cpu_delta"],
                    f"record {sequence} CPU")
                validate_cpu_triplet(
                    record["sibling_before"], record["sibling_after"],
                    record["sibling_delta"], f"record {sequence} sibling")
                require(record_cpu["nonidle_total"] > 0,
                        f"record {sequence} timed CPU recorded no work")
                artifacts = require_keys(record["artifacts"], {"json", "stdout", "stderr"},
                                         f"record {sequence} artifacts")
                prefix = f"{case['name']}-round{round_index}-slot{slot}-{variant}"
                json_name = prefix + ".json"
                json_path = raw_root / json_name
                json_sha = validate_artifact(
                    artifacts["json"], json_name, json_path,
                    f"record {sequence} JSON")
                validate_artifact(artifacts["stdout"], prefix + ".stdout",
                                  raw_root / (prefix + ".stdout"),
                                  f"record {sequence} stdout", empty=True)
                validate_artifact(artifacts["stderr"], prefix + ".stderr",
                                  raw_root / (prefix + ".stderr"),
                                  f"record {sequence} stderr", empty=True)
                digest_sha = hex_digest(record["workload_digests_sha256"],
                                        f"record {sequence} workload digest SHA-256")
                result = validate_current_result(
                    json_path, case, variant, iterations, json_sha, digest_sha)
                coordinate = (case["name"], round_index, slot)
                require(coordinate not in seen, "duplicate manifest record")
                seen.add(coordinate)
                grouped.setdefault(case["name"], {}).setdefault(round_index, {})[slot] = result
                raw_hashes.append(json_sha)
                sequence += 1

    cells = []
    for case in cases:
        round_logs = []
        setup_logs = []
        scratch = {"control": set(), "candidate": set()}
        digest_reference = None
        cell_raw_hashes = []
        for round_index, order in enumerate(orders):
            slots = grouped[case["name"]][round_index]
            require(set(slots) == {0, 1, 2, 3}, "incomplete ABBA round")
            for result in slots.values():
                scratch[result["variant"]].add(result["scratch_bytes"])
                if digest_reference is None:
                    digest_reference = result["digests"]
                require(result["digests"] == digest_reference,
                        f"{case['name']}: control/candidate workload digests differ")
                cell_raw_hashes.append(result["raw_sha256"])
            pair_logs = []
            pair_setup_logs = []
            for left, right in ((0, 1), (2, 3)):
                pair = (slots[left], slots[right])
                control = pair[0] if pair[0]["variant"] == "control" else pair[1]
                candidate = pair[1] if pair[1]["variant"] == "candidate" else pair[0]
                pair_logs.append(math.log(control["decode_us"] / candidate["decode_us"]))
                pair_setup_logs.append(math.log(control["setup_us"] / candidate["setup_us"]))
            round_logs.append(statistics.fmean(pair_logs))
            setup_logs.append(statistics.fmean(pair_setup_logs))
        require(all(len(values) == 1 for values in scratch.values()),
                f"{case['name']}: scratch size changed within a role")
        materialized = next(iter(scratch["control"]))
        tiled = next(iter(scratch["candidate"]))
        require(materialized > 0 and 0 < tiled <= materialized,
                f"{case['name']}: explicit-mode scratch relation is invalid")
        cells.append({
            **case,
            "decode_execution": clustered_summary(round_logs),
            "decode_plan_setup": clustered_summary(setup_logs),
            "scratch": {
                "materialized_bytes": materialized,
                "tiled_bytes": tiled,
                "reduction_percent": 100.0 * (materialized - tiled) / materialized,
            },
            "workload_digests": digest_reference,
            "raw_result_set_sha256": canonical_sha256(cell_raw_hashes),
        })

    sibling_fraction = (campaign_sibling_delta["nonidle_total"] /
                        campaign_sibling_delta["total"]
                        if campaign_sibling_delta["total"] else 0.0)
    summary: dict[str, Any] = {
        "schema": CURRENT_SUMMARY_SCHEMA,
        "source_manifest_sha256": sha256_bytes(manifest_raw),
        "source_manifest_content_sha256": content_digest,
        "source": source,
        "build": build,
        "runner_sha256": runner["sha256"],
        "analyzer_sha256": sha256_file(Path(__file__)),
        "lease_source_sha256": lease_source["sha256"],
        "raw_result_set_sha256": canonical_sha256(raw_hashes),
        "method": {
            "binary_relation": "same executable and source identity",
            "control_mode": "--force-materialized",
            "candidate_mode": "--force-tiled",
            "cpu": cpu, "sibling": sibling,
            "lease_identity_sha256": lease_digest,
            "round_orders": orders,
            "confidence_intervals": "clustered paired-log Student-t 95%",
            "ratio_orientation": "materialized_time_over_tiled_time",
            "raw_statistics_recomputed": True,
            "sibling_nonidle_fraction": sibling_fraction,
            "validity_is_independent_of_speed": True,
        },
        "cells": cells,
        "status": "pass",
    }
    summary["content_sha256"] = canonical_sha256(summary)
    return summary


# Historical v1 replay -----------------------------------------------------

def historical_command_options(command: object) -> dict[str, str]:
    require(isinstance(command, list) and all(isinstance(x, str) for x in command),
            "historical record command is not a string array")
    options: dict[str, str] = {}
    boolean = {"--force-specialized", "--skip-legacy", "--retain-samples"}
    index = 0
    while index < len(command):
        item = command[index]
        if item.startswith("--"):
            require(item not in options, f"duplicate historical command option {item}")
            if item in boolean:
                options[item] = "true"
            else:
                require(index + 1 < len(command), f"missing value after {item}")
                options[item] = command[index + 1]
                index += 1
        index += 1
    require("--force-tiled" not in options and "--force-materialized" not in options,
            "historical command contains a current selector")
    return options


def validate_historical_result(path: Path, case: dict[str, Any], variant: str,
                               expected_sha: str) -> dict[str, Any]:
    require(path.is_file() and sha256_file(path) == expected_sha,
            f"historical benchmark result hash changed: {path}")
    document, _ = load_json_bytes(path, "historical benchmark result")
    require(document.get("schema") == BENCHMARK_SCHEMA,
            f"unexpected historical benchmark schema in {path}")
    parameters = document.get("parameters")
    resolved = document.get("resolved")
    require(isinstance(parameters, dict) and isinstance(resolved, dict),
            f"missing historical benchmark identity in {path}")
    validate_inactive_selector_pair(parameters, str(path))
    expected = {
        "K": case["K"], "R": case["R"], "requested_profile": "legacy_high_v1",
        "requested_field": "gf8", "requested_backend": case["backend"],
        "force_generic_decode": False, "force_specialized_decode": True,
        "skip_legacy": True, "retain_samples": True,
        "shard_bytes": case["shard_bytes"], "loss_count": case["loss"],
        "batch": case["batch"], "reuse": case["reuse"],
        "warmup": case["warmup"], "thread_count": 1, "seed": BENCHMARK_SEED,
    }
    for key, expected_value in expected.items():
        require(parameters.get(key) == expected_value,
                f"{path}: historical {key} mismatch")
    require(resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and resolved.get("backend") == case["backend"] and
            resolved.get("parent_count") == case["parent_count"] and
            resolved.get("padded_side") == case["padded_side"],
            f"{path}: historical resolved codec mismatch")
    require(document.get("correctness", {}).get("leopard2_round_trip") is True,
            f"{path}: historical round trip failed")
    execution = document.get("metrics", {}).get("decode_execution", {})
    setup = document.get("metrics", {}).get("decode_plan_setup", {})
    memory = document.get("memory", {})
    digests = validate_digests(document.get("workload_digests"), f"{path} digests")
    return {
        "variant": variant,
        "decode_us": finite(execution.get("median_us_per_batch_call"),
                            f"{path} decode median", 1e-300),
        "setup_us": finite(setup.get("median_us"), f"{path} setup median", 1e-300),
        "scratch_bytes": strict_int(memory.get("decode_scratch_bytes_per_stripe"),
                                    f"{path} scratch", 1),
        "digests": digests,
    }


def analyze_historical(manifest_path: Path, require_rounds: int) -> dict[str, Any]:
    manifest, raw_manifest = load_json_bytes(manifest_path, "historical manifest")
    require(manifest.get("schema") == HISTORICAL_SCHEMA and manifest.get("valid") is True,
            "historical manifest schema/validity is wrong")
    git_oid(manifest.get("control_commit"), "historical control commit")
    git_oid(manifest.get("candidate_commit"), "historical candidate commit")
    identities = manifest.get("identities")
    require(isinstance(identities, dict) and identities, "missing historical identities")
    for key, value in identities.items():
        hex_digest(value, f"historical identity {key}")
    cases = [normalize_case(value) for value in manifest.get("cases", [])]
    require(cases and len({case["name"] for case in cases}) == len(cases),
            "historical case names are empty or duplicated")
    case_by_name = {case["name"]: case for case in cases}
    orders = manifest.get("round_orders")
    require(isinstance(orders, list) and len(orders) == require_rounds,
            f"historical manifest must contain exactly {require_rounds} rounds")
    for order in orders:
        require(order in (["control", "candidate", "candidate", "control"],
                          ["candidate", "control", "control", "candidate"]),
                "historical round is not ABBA or BAAB")
    records = manifest.get("records")
    require(isinstance(records, list) and len(records) == len(cases) * len(orders) * 4,
            "historical record count mismatch")
    raw_root = manifest_path.parent / "raw"
    grouped: dict[str, dict[int, dict[int, dict[str, Any]]]] = {}
    seen = set()
    for record in records:
        require(isinstance(record, dict), "historical record is not an object")
        name, round_index, slot, variant = (record.get("case"), record.get("round"),
                                            record.get("slot"), record.get("variant"))
        require(name in case_by_name and isinstance(round_index, int) and
                isinstance(slot, int) and 0 <= round_index < len(orders) and 0 <= slot < 4,
                "historical record coordinate is invalid")
        require(variant == orders[round_index][slot],
                "historical record variant/order mismatch")
        coordinate = (name, round_index, slot)
        require(coordinate not in seen, "duplicate historical record")
        seen.add(coordinate)
        options = historical_command_options(record.get("command"))
        case = case_by_name[name]
        expected = {
            "--k": str(case["K"]), "--r": str(case["R"]), "--profile": "high",
            "--field": "gf8", "--backend": case["backend"],
            "--bytes": str(case["shard_bytes"]), "--loss": str(case["loss"]),
            "--batch": str(case["batch"]), "--reuse": str(case["reuse"]),
            "--iterations": str(manifest.get("iterations")),
            "--warmup": str(case["warmup"]), "--threads": "1",
            "--seed": str(BENCHMARK_SEED), "--force-specialized": "true",
            "--skip-legacy": "true", "--retain-samples": "true",
        }
        for key, expected_value in expected.items():
            require(options.get(key) == expected_value,
                    f"historical record command {key} mismatch")
        relative = record.get("json")
        require(isinstance(relative, str) and Path(relative).name == relative,
                "historical result path is not a basename")
        json_path = raw_root / relative
        for suffix in (".stdout", ".stderr"):
            stream = json_path.with_suffix(suffix)
            require(stream.is_file() and stream.stat().st_size == 0,
                    f"historical benchmark stream is missing/nonempty: {stream}")
        result = validate_historical_result(
            json_path, case, variant,
            hex_digest(record.get("json_sha256"), "historical result digest"))
        grouped.setdefault(name, {}).setdefault(round_index, {})[slot] = result
    cells = []
    for case in cases:
        round_logs, setup_logs = [], []
        scratch = {"control": set(), "candidate": set()}
        digest_reference = None
        for round_index, order in enumerate(orders):
            slots = grouped[case["name"]][round_index]
            require(set(slots) == {0, 1, 2, 3}, "incomplete historical ABBA round")
            for result in slots.values():
                scratch[result["variant"]].add(result["scratch_bytes"])
                if digest_reference is None:
                    digest_reference = result["digests"]
                require(result["digests"] == digest_reference,
                        "historical control/candidate workload digests differ")
            contrasts, setup_contrasts = [], []
            for left, right in ((0, 1), (2, 3)):
                pair = (slots[left], slots[right])
                control = pair[0] if pair[0]["variant"] == "control" else pair[1]
                candidate = pair[1] if pair[1]["variant"] == "candidate" else pair[0]
                contrasts.append(math.log(control["decode_us"] / candidate["decode_us"]))
                setup_contrasts.append(math.log(control["setup_us"] / candidate["setup_us"]))
            round_logs.append(statistics.fmean(contrasts))
            setup_logs.append(statistics.fmean(setup_contrasts))
        require(all(len(item) == 1 for item in scratch.values()),
                "historical scratch changed within a variant")
        control_scratch = next(iter(scratch["control"]))
        candidate_scratch = next(iter(scratch["candidate"]))
        require(control_scratch > 0 and 0 < candidate_scratch <= control_scratch,
                "historical candidate scratch is invalid")
        cells.append({
            **case,
            "decode_execution": clustered_summary(round_logs),
            "decode_plan_setup": clustered_summary(setup_logs),
            "scratch": {
                "materialized_bytes": control_scratch,
                "tiled_bytes": candidate_scratch,
                "reduction_percent": 100.0 *
                    (control_scratch - candidate_scratch) / control_scratch,
            },
            "workload_digests": digest_reference,
        })
    summary: dict[str, Any] = {
        "schema": HISTORICAL_SUMMARY_SCHEMA,
        "source_manifest_sha256": sha256_bytes(raw_manifest),
        "control_commit": manifest["control_commit"],
        "candidate_commit": manifest["candidate_commit"],
        "identities": identities,
        "method": {
            "cpu": manifest.get("cpu"), "sibling": manifest.get("sibling"),
            "lease": manifest.get("lease"), "round_orders": orders,
            "confidence_intervals": "clustered paired-log Student-t 95%",
            "ratio_orientation": "materialized_time_over_tiled_time",
            "validity_is_independent_of_speed": True,
        },
        "cells": cells,
        "status": "pass",
    }
    summary["content_sha256"] = canonical_sha256(summary)
    return summary


def analyze(manifest_path: Path, require_rounds: int,
            source_override: Path | None = None,
            validate_live: bool = True) -> dict[str, Any]:
    manifest, _ = load_json_bytes(manifest_path, "manifest")
    schema = manifest.get("schema")
    if schema == CURRENT_SCHEMA:
        return analyze_current(manifest_path, require_rounds, source_override, validate_live)
    if schema == HISTORICAL_SCHEMA:
        require(source_override is None, "historical replay does not accept source override")
        return analyze_historical(manifest_path, require_rounds)
    raise EvidenceError(f"unsupported tiled-high evidence schema: {schema!r}")


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True, allow_nan=False)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


# Self-test fixtures -------------------------------------------------------

def fixture_stats(samples: list[float], execution: bool,
                  rates: dict[str, int] | None = None) -> dict[str, Any]:
    median = statistics.median(samples)
    mad = statistics.median(abs(value - median) for value in samples)
    suffix = "_us_per_batch_call" if execution else "_us"
    result: dict[str, Any] = {
        "median" + suffix: median,
        "mad" + suffix: mad,
        "minimum" + suffix: min(samples),
        "maximum" + suffix: max(samples),
        "samples" + suffix: samples,
    }
    if rates:
        for name, byte_count in rates.items():
            result[name] = byte_count / (median * 1000.0)
    return result


def fixture_file(path: Path, root: Path, relative: str) -> dict[str, Any]:
    return {
        "path": str(path), "realpath": str(path), "relative_path": relative,
        "sha256": sha256_file(path), "size": path.stat().st_size,
        "mode": path.stat().st_mode & 0o7777, "device": path.stat().st_dev,
        "inode": path.stat().st_ino,
    }


def fixture_delta(before: dict[str, int], after: dict[str, int]) -> dict[str, int]:
    return recompute_cpu_delta(before, after)


def make_fixture(root: Path) -> Path:
    source_root = root / "source"
    raw_root = root / "evidence" / "raw"
    (source_root / Path(RUNNER_RELATIVE).parent).mkdir(parents=True)
    (source_root / Path(LEASE_RELATIVE).parent).mkdir(parents=True)
    binary = source_root / "build-audit/bench_leopard2"
    cache = source_root / "build-audit/CMakeCache.txt"
    binary.parent.mkdir(parents=True)
    binary.write_bytes(b"fixture benchmark\n")
    binary.chmod(0o755)
    cache.write_text(f"CMAKE_HOME_DIRECTORY:INTERNAL={source_root}\n", encoding="utf-8")
    runner = source_root / RUNNER_RELATIVE
    runner.write_text("fixture runner\n", encoding="utf-8")
    lease_path = source_root / LEASE_RELATIVE
    lease_path.write_text("fixture lease\n", encoding="utf-8")
    for command in (
            ["git", "init", "-q"],
            ["git", "add", "."],
            ["git", "-c", "user.name=Leopard2 Self Test",
             "-c", "user.email=leopard2-self-test@example.invalid",
             "commit", "-q", "-m", "fixture"]):
        try:
            subprocess.run(command, cwd=source_root, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except (OSError, subprocess.CalledProcessError) as error:
            raise EvidenceError(f"cannot create self-test Git fixture: {error}") from error
    head = checked_git(source_root, ["rev-parse", "HEAD"])
    tree = checked_git(source_root, ["rev-parse", "HEAD^{tree}"])
    raw_root.mkdir(parents=True)
    source = {
        "root": str(source_root), "head": head, "tree": tree,
        "expected_commit": head, "status": "clean", "status_sha256": EMPTY_SHA256,
    }
    build = {
        "binary": fixture_file(binary, source_root, "build-audit/bench_leopard2"),
        "cache": fixture_file(cache, source_root, "build-audit/CMakeCache.txt"),
        "cmake_home": str(source_root),
    }
    runner_id = fixture_file(runner, source_root, RUNNER_RELATIVE)
    lease_source = fixture_file(lease_path, source_root, LEASE_RELATIVE)
    cpu, sibling = 2, 3
    payload = {"cpus": [2, 3], "schema": PAIR_LEASE_SCHEMA, "uid": os.getuid()}
    lease = {
        "lock": PAIR_LEASE_LOCK,
        "path": str(Path("/run/user") / str(os.getuid()) / "leopard2-cpu-leases" /
                    f"leopard2-cpu-pair-{os.getuid()}-2-3.lock"),
        "payload": payload, "sha256": canonical_sha256(payload),
        "device": 1, "inode": 2, "directory_device": 1, "directory_inode": 3,
        "kernel_lease": {"schema": KERNEL_LEASE_SCHEMA,
                         "mechanism": KERNEL_LEASE_MECHANISM,
                         "name_sha256": "3" * 64},
    }
    lease_digest = canonical_sha256(lease)
    cases = [["fixture", 224, 32, 1024, 1, 2, 8, "avx2", 1]]
    case = normalize_case(cases[0])
    workload = {"algorithm": "fnv1a64", "original_data": "1" * 16,
                "transmitted_parity": "2" * 16, "recovered_originals": "3" * 16}
    records = []
    campaign_start, campaign_end = 1000, 100000
    sequence = 0
    for round_index, order in enumerate(ROUND_ORDERS):
        for slot, variant in enumerate(order):
            prefix = f"fixture-round{round_index}-slot{slot}-{variant}"
            json_path = raw_root / (prefix + ".json")
            stdout_path = raw_root / (prefix + ".stdout")
            stderr_path = raw_root / (prefix + ".stderr")
            decode_samples = ([20.0, 21.0, 22.0] if variant == "control"
                              else [17.0, 18.0, 19.0])
            encode_input = case["K"] * case["shard_bytes"]
            encode_output = case["R"] * case["shard_bytes"]
            decode_input = (case["K"] - case["loss"] + case["R"]) * case["shard_bytes"]
            decode_output = case["loss"] * case["shard_bytes"]
            setup_samples = [1.8, 2.0, 2.2]
            setup_median = statistics.median(setup_samples)
            decode_median = statistics.median(decode_samples)
            amortized = decode_median + setup_median / case["reuse"]
            raw_document = {
                "schema": BENCHMARK_SCHEMA,
                "build": {"compiler": "fixture", "compiler_version": "1", "cplusplus": 201103},
                "parameters": {
                    "K": case["K"], "R": case["R"],
                    "requested_profile": "legacy_high_v1", "requested_field": "gf8",
                    "requested_backend": "avx2", "force_generic_decode": False,
                    "force_specialized_decode": True,
                    "force_tiled_decode": variant == "candidate",
                    "force_materialized_decode": variant == "control",
                    "skip_legacy": True, "retain_samples": True,
                    "shard_bytes": case["shard_bytes"], "loss_count": case["loss"],
                    "missing_original_indices": list(range(case["loss"])),
                    "batch": 1, "reuse": 2, "iterations": 3, "warmup": 1,
                    "thread_count": 1, "seed": BENCHMARK_SEED,
                },
                "resolved": {"profile": "legacy_high_v1", "field": "gf8",
                             "backend": "avx2", "thread_count": 1,
                             "parent_count": 256, "padded_side": 32},
                "correctness": {"leopard2_round_trip": True, "legacy_comparison": None},
                "workload_digests": workload,
                "memory": {"scratch_alignment": 64,
                           "encode_scratch_bytes_per_stripe": 1024,
                           "decode_scratch_bytes_per_stripe":
                               4096 if variant == "control" else 2048,
                           "encode_scratch_bytes_batch": 1024,
                           "decode_scratch_bytes_batch":
                               4096 if variant == "control" else 2048},
                "metrics": {
                    "codec_setup": fixture_stats([0.9, 1.0, 1.1], False),
                    "encode_execution": fixture_stats(
                        [9.0, 10.0, 11.0], True,
                        {"input_GB_per_s": encode_input,
                         "parity_output_GB_per_s": encode_output}),
                    "decode_plan_setup": fixture_stats(setup_samples, False),
                    "decode_execution": fixture_stats(
                        decode_samples, True,
                        {"offered_received_GB_per_s": decode_input,
                         "repaired_output_GB_per_s": decode_output}),
                    "decode_amortized_at_reuse": {
                        "reuse_count": 2,
                        "derived_median_us_per_batch_call": amortized,
                        "offered_received_GB_per_s": decode_input / (amortized * 1000.0),
                        "repaired_output_GB_per_s": decode_output / (amortized * 1000.0),
                    },
                    "rate_semantics":
                        "offered_received counts all non-null shard pointers supplied; "
                        "a plan may read a deterministic subset",
                },
                "legacy": {"available": False,
                           "unavailable_reason": "disabled by --skip-legacy",
                           "codec_setup": None, "decode_timing_includes_setup": True,
                           "encode_execution": None, "decode_including_setup": None},
            }
            json_path.write_text(json.dumps(raw_document, indent=2) + "\n", encoding="utf-8")
            stdout_path.write_bytes(b"")
            stderr_path.write_bytes(b"")
            json_at_run = str(json_path)
            command = expected_command(str(binary), case, json_at_run, cpu, 3, variant)
            before = {key: 1000 + sequence * 10 for key in CPU_FIELDS}
            after = dict(before)
            after["user"] += 2
            after["idle"] += 3
            sibling_before = {key: 2000 + sequence * 10 for key in CPU_FIELDS}
            sibling_after = dict(sibling_before)
            sibling_after["idle"] += 5
            artifacts = {
                "json": {"name": json_path.name, "size": json_path.stat().st_size,
                         "sha256": sha256_file(json_path)},
                "stdout": {"name": stdout_path.name, "size": 0, "sha256": EMPTY_SHA256},
                "stderr": {"name": stderr_path.name, "size": 0, "sha256": EMPTY_SHA256},
            }
            records.append({
                "sequence": sequence, "case": "fixture", "round": round_index,
                "slot": slot, "variant": variant, "mode": ROLE_MODE[variant],
                "selector": MODE_SELECTOR[ROLE_MODE[variant]], "command": command,
                "command_sha256": canonical_sha256(command),
                "benchmark_argv0": str(binary), "json_path_at_run": json_at_run,
                "iterations": 3, "seed": BENCHMARK_SEED,
                "lease_identity_sha256": lease_digest,
                "start_time_ns": 2000 + sequence * 100,
                "end_time_ns": 2050 + sequence * 100,
                "cpu_before": before, "cpu_after": after,
                "cpu_delta": fixture_delta(before, after),
                "sibling_before": sibling_before, "sibling_after": sibling_after,
                "sibling_delta": fixture_delta(sibling_before, sibling_after),
                "artifacts": artifacts,
                "workload_digests_sha256": canonical_sha256(workload),
            })
            sequence += 1
    campaign_before = {key: 100 for key in CPU_FIELDS}
    campaign_after = dict(campaign_before)
    campaign_after["user"] += 100
    campaign_after["idle"] += 1000
    sibling_before = {key: 200 for key in CPU_FIELDS}
    sibling_after = dict(sibling_before)
    sibling_after["idle"] += 1100
    manifest: dict[str, Any] = {
        "schema": CURRENT_SCHEMA, "status": "complete", "valid": True,
        "source": source, "build": build, "runner": runner_id,
        "lease_source": lease_source,
        "campaign": {
            "cpu": cpu, "sibling": sibling, "original_affinity": [0, 1, 2, 3],
            "housekeeping_affinity": [0, 1], "matrix": "self-test", "iterations": 3,
            "seed": BENCHMARK_SEED, "cases": cases,
            "round_orders": [list(order) for order in ROUND_ORDERS],
            "lease": lease, "lease_identity_sha256": lease_digest,
            "start_time_ns": campaign_start, "end_time_ns": campaign_end,
            "cpu_before": campaign_before, "cpu_after": campaign_after,
            "cpu_delta": fixture_delta(campaign_before, campaign_after),
            "sibling_before": sibling_before, "sibling_after": sibling_after,
            "sibling_delta": fixture_delta(sibling_before, sibling_after),
        },
        "records": records,
        "final": {"source": copy.deepcopy(source), "build": copy.deepcopy(build),
                  "runner": copy.deepcopy(runner_id),
                  "lease_source": copy.deepcopy(lease_source)},
    }
    manifest["content_sha256"] = canonical_sha256(manifest)
    path = root / "evidence/manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def rewrite_manifest(path: Path, manifest: dict[str, Any]) -> None:
    manifest.pop("content_sha256", None)
    manifest["content_sha256"] = canonical_sha256(manifest)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def refresh_record_json(path: Path, manifest: dict[str, Any], index: int,
                        document: dict[str, Any]) -> None:
    record = manifest["records"][index]
    artifact = record["artifacts"]["json"]
    raw = path.parent / "raw" / artifact["name"]
    raw.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    artifact["size"] = raw.stat().st_size
    artifact["sha256"] = sha256_file(raw)
    record["workload_digests_sha256"] = canonical_sha256(document["workload_digests"])


def self_test() -> None:
    require(len(versioned_checkpoint_cases()) == 117 and
            len({case[0] for case in versioned_checkpoint_cases()}) == 117,
            "versioned checkpoint matrix changed")
    require(len(versioned_smoke_cases()) == 2,
            "versioned smoke matrix changed")
    validate_inactive_selector_pair({}, "absent historical selector pair")
    validate_inactive_selector_pair(
        {"force_tiled_decode": False, "force_materialized_decode": False},
        "inactive mixed selector pair")
    for label, selector_pair in (
            ("partial-tiled", {"force_tiled_decode": False}),
            ("partial-materialized", {"force_materialized_decode": False}),
            ("active-tiled", {"force_tiled_decode": True,
                              "force_materialized_decode": False}),
            ("active-materialized", {"force_tiled_decode": False,
                                     "force_materialized_decode": True}),
            ("active-both", {"force_tiled_decode": True,
                             "force_materialized_decode": True}),
            ("non-boolean", {"force_tiled_decode": 0,
                             "force_materialized_decode": False})):
        try:
            validate_inactive_selector_pair(selector_pair, label)
        except EvidenceError:
            pass
        else:
            raise EvidenceError("historical/mixed selector mutation was accepted: " + label)
    summary = clustered_summary([math.log(1.10), math.log(1.11), math.log(1.09)])
    require(summary["credible_gain_at_least_5_percent"], "gain classification failed")
    regression = clustered_summary([math.log(0.95), math.log(0.96), math.log(0.95)])
    require(regression["credible_regression_over_2_percent"],
            "regression classification failed")
    with tempfile.TemporaryDirectory(prefix="tiled-high-evidence-self-test-") as temporary:
        root = Path(temporary)
        base = make_fixture(root / "base")
        result = analyze_current(base, 3, validate_live=True, allow_self_test=True)
        require(len(result["cells"]) == 1 and result["status"] == "pass",
                "valid current fixture was rejected")
        try:
            analyze_current(base, 2, validate_live=False, allow_self_test=True)
        except EvidenceError:
            pass
        else:
            raise EvidenceError("current evidence accepted a truncated round requirement")

        def rejected(label: str, mutation, live: bool = False) -> None:
            attack_root = root / label
            path = make_fixture(attack_root)
            manifest, _ = load_json_bytes(path, label + " manifest")
            mutation(path, manifest)
            rewrite_manifest(path, manifest)
            try:
                analyze_current(path, 3, validate_live=live, allow_self_test=True)
            except EvidenceError:
                return
            raise EvidenceError("adversarial fixture was accepted: " + label)

        def coherent_role_relabel(path, m):
            swap = {"control": "candidate", "candidate": "control"}
            m["campaign"]["round_orders"] = [
                [swap[item] for item in order]
                for order in m["campaign"]["round_orders"]]
            case = normalize_case(m["campaign"]["cases"][0])
            for record in m["records"]:
                old_variant = record["variant"]
                new_variant = swap[old_variant]
                old_prefix = (f"fixture-round{record['round']}-slot{record['slot']}-"
                              f"{old_variant}")
                new_prefix = (f"fixture-round{record['round']}-slot{record['slot']}-"
                              f"{new_variant}")
                for kind, suffix in (("json", ".json"), ("stdout", ".stdout"),
                                     ("stderr", ".stderr")):
                    old_path = path.parent / "raw" / (old_prefix + suffix)
                    new_path = path.parent / "raw" / (new_prefix + suffix)
                    old_path.rename(new_path)
                    record["artifacts"][kind]["name"] = new_path.name
                raw = path.parent / "raw" / (new_prefix + ".json")
                document, _ = load_json_bytes(raw, "coherent relabel raw")
                document["parameters"]["force_tiled_decode"] = new_variant == "candidate"
                document["parameters"]["force_materialized_decode"] = new_variant == "control"
                raw.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
                record["artifacts"]["json"]["size"] = raw.stat().st_size
                record["artifacts"]["json"]["sha256"] = sha256_file(raw)
                record["variant"] = new_variant
                record["mode"] = ROLE_MODE[new_variant]
                record["selector"] = MODE_SELECTOR[ROLE_MODE[new_variant]]
                record["json_path_at_run"] = str(raw)
                record["command"] = expected_command(
                    record["benchmark_argv0"], case, str(raw), 2, 3, new_variant)
                record["command_sha256"] = canonical_sha256(record["command"])
        rejected("coherent-role-relabel", coherent_role_relabel)
        rejected("truncated-smoke-matrix", lambda _p, m:
                 m["campaign"].__setitem__("matrix", "smoke"))

        def swapped_executable(_p, m):
            record = m["records"][0]
            record["benchmark_argv0"] = "/tmp/swapped-bench"
            record["command"][7] = "/tmp/swapped-bench"
            record["command_sha256"] = canonical_sha256(record["command"])
        rejected("swapped-executable", swapped_executable)

        def command_rehash(_p, m):
            record = m["records"][0]
            index = record["command"].index("--seed") + 1
            record["command"][index] = "9"
            record["seed"] = 9
            record["command_sha256"] = canonical_sha256(record["command"])
        rejected("command-rehash", command_rehash)

        def raw_identity_rehash(path, m):
            raw = path.parent / "raw" / m["records"][0]["artifacts"]["json"]["name"]
            document, _ = load_json_bytes(raw, "attack raw")
            document["parameters"]["K"] = 223
            refresh_record_json(path, m, 0, document)
        rejected("raw-identity-rehash", raw_identity_rehash)

        def sample_attack(path, m):
            raw = path.parent / "raw" / m["records"][0]["artifacts"]["json"]["name"]
            document, _ = load_json_bytes(raw, "attack raw")
            document["metrics"]["decode_execution"]["samples_us_per_batch_call"][0] = 50.0
            refresh_record_json(path, m, 0, document)
        rejected("sample-rehash", sample_attack)

        def summary_attack(path, m):
            raw = path.parent / "raw" / m["records"][0]["artifacts"]["json"]["name"]
            document, _ = load_json_bytes(raw, "attack raw")
            document["metrics"]["decode_execution"]["median_us_per_batch_call"] += 1.0
            refresh_record_json(path, m, 0, document)
        rejected("summary-rehash", summary_attack)

        def rate_attack(path, m):
            raw = path.parent / "raw" / m["records"][0]["artifacts"]["json"]["name"]
            document, _ = load_json_bytes(raw, "attack raw")
            document["metrics"]["decode_execution"]["offered_received_GB_per_s"] += 1.0
            refresh_record_json(path, m, 0, document)
        rejected("derived-rate-rehash", rate_attack)

        def coherent_lease_rehash(_p, m):
            lease = m["campaign"]["lease"]
            lease["payload"]["cpus"] = [2, 4]
            lease["path"] = str(Path(lease["path"]).with_name(
                f"leopard2-cpu-pair-{lease['payload']['uid']}-2-4.lock"))
            lease["sha256"] = canonical_sha256(lease["payload"])
            digest = canonical_sha256(lease)
            m["campaign"]["lease_identity_sha256"] = digest
            for record in m["records"]:
                record["lease_identity_sha256"] = digest
        rejected("coherent-lease-rehash", coherent_lease_rehash)

        def coherent_source_rehash(_p, m):
            for source_key in (m["source"], m["final"]["source"]):
                source_key["head"] = "4" * 40
                source_key["expected_commit"] = "4" * 40
                source_key["tree"] = "5" * 40
        rejected("coherent-source-rehash", coherent_source_rehash, live=True)
        rejected("cpu-delta-rehash", lambda _p, m:
                 m["records"][0]["cpu_delta"].__setitem__("nonidle_total", 99))

        def digest_attack(path, m):
            raw = path.parent / "raw" / m["records"][0]["artifacts"]["json"]["name"]
            document, _ = load_json_bytes(raw, "attack raw")
            document["workload_digests"]["recovered_originals"] = "f" * 16
            refresh_record_json(path, m, 0, document)
        rejected("workload-digest-rehash", digest_attack)

        def both_selectors(path, m):
            raw = path.parent / "raw" / m["records"][0]["artifacts"]["json"]["name"]
            document, _ = load_json_bytes(raw, "attack raw")
            document["parameters"]["force_tiled_decode"] = True
            document["parameters"]["force_materialized_decode"] = True
            refresh_record_json(path, m, 0, document)
        rejected("both-selectors", both_selectors)

        def numeric_selector(path, m, key, value):
            raw = path.parent / "raw" / m["records"][0]["artifacts"]["json"]["name"]
            document, _ = load_json_bytes(raw, "attack raw")
            document["parameters"][key] = value
            refresh_record_json(path, m, 0, document)
        rejected("numeric-false-tiled-selector", lambda p, m:
                 numeric_selector(p, m, "force_tiled_decode", 0))
        rejected("numeric-true-materialized-selector", lambda p, m:
                 numeric_selector(p, m, "force_materialized_decode", 1))

        def partial_selector(path, m):
            raw = path.parent / "raw" / m["records"][0]["artifacts"]["json"]["name"]
            document, _ = load_json_bytes(raw, "attack raw")
            del document["parameters"]["force_tiled_decode"]
            refresh_record_json(path, m, 0, document)
        rejected("partial-selector", partial_selector)

        def stdout_attack(path, m):
            record = m["records"][0]
            artifact = record["artifacts"]["stdout"]
            stream = path.parent / "raw" / artifact["name"]
            stream.write_bytes(b"noise")
            artifact["size"] = 5
            artifact["sha256"] = sha256_file(stream)
        rejected("stdout-rehash", stdout_attack)

        path = make_fixture(root / "unsigned-content")
        manifest, _ = load_json_bytes(path, "unsigned attack")
        manifest["records"][0]["seed"] = 9
        path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        try:
            analyze_current(path, 3, validate_live=False, allow_self_test=True)
        except EvidenceError:
            pass
        else:
            raise EvidenceError("unsigned manifest mutation was accepted")
    print("tiled-high evidence analyzer self-test passed: current v2 + 18 attacks")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest")
    parser.add_argument("--output")
    parser.add_argument("--require-rounds", type=int, default=3)
    parser.add_argument("--source-root",
                        help="relocate current source/build closure for replay on another machine")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.manifest or not args.output:
        parser.error("--manifest and --output are required unless --self-test is used")
    source_override = Path(args.source_root) if args.source_root else None
    summary = analyze(Path(args.manifest), args.require_rounds, source_override)
    write_json(Path(args.output), summary)
    print(json.dumps({"cells": len(summary["cells"]),
                      "content_sha256": summary["content_sha256"],
                      "schema": summary["schema"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except EvidenceError as error:
        print(f"analyze_tiled_high.py: {error}", file=os.sys.stderr)
        raise SystemExit(1)
