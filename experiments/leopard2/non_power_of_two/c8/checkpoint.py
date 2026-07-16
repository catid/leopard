#!/usr/bin/env python3
"""Fail-closed evidence merger for the default-off Leopard2 C8 experiment."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any


EXPECTED_CORRECTNESS = {
    "cases": 19,
    "coefficients": 78362,
    "byte_comparisons": 99420,
    "legacy_identity_cases": 66,
    "legacy_difference_cases": 8,
    "hot_path_allocations": 0,
    "digest": "0x4bd56b387de0757e",
}
EXPECTED_BACKENDS = ("auto", "scalar", "ssse3", "avx2")
EXPECTED_DIFFERENTIAL = (
    ("gf8",1,1,4,1,True,True,False),
    ("gf8",3,1,4,3,True,True,False),
    ("gf8",5,3,4,15,False,True,True),
    ("gf8",12,4,4,48,True,True,False),
    ("gf8",24,8,4,192,True,True,False),
    ("gf8",31,3,4,93,False,True,True),
    ("gf8",48,16,4,768,True,True,False),
    ("gf8",60,4,4,240,True,True,False),
    ("gf8",63,1,4,63,True,True,False),
    ("gf8",96,32,4,3072,True,True,False),
    ("gf8",120,7,4,840,False,True,True),
    ("gf8",120,8,4,960,True,True,False),
    ("gf8",120,9,4,1080,False,True,True),
    ("gf8",192,64,4,12288,True,False,False),
    ("gf8",224,31,4,6944,False,True,True),
    ("gf16",255,1,4,255,True,True,False),
    ("gf16",500,3,4,1500,False,True,True),
    ("gf16",1000,17,4,17000,False,True,True),
    ("gf16",1000,33,4,33000,False,True,True),
)
EXPECTED_BENCHMARKS = (
    ("gf8",31,2,64,8,8), ("gf8",31,2,1024,8,4),
    ("gf8",31,3,64,8,8), ("gf8",31,3,1024,8,4),
    ("gf8",31,3,65536,1,1),
    ("gf8",31,4,64,8,8), ("gf8",31,4,1024,8,4),
    ("gf8",60,4,1024,8,4),
    ("gf8",60,4,65536,1,1), ("gf8",120,7,64,8,8),
    ("gf8",120,7,1024,8,4), ("gf8",120,7,65536,1,1),
    ("gf8",120,8,1024,8,4), ("gf8",120,9,1024,8,4),
    ("gf8",224,31,64,8,4), ("gf8",224,31,1024,1,1),
    ("gf16",500,2,64,8,4), ("gf16",500,2,1024,1,1),
    ("gf16",500,3,64,8,4), ("gf16",500,3,1024,1,1),
    ("gf16",500,4,64,8,4), ("gf16",500,4,1024,1,1),
    ("gf16",500,5,64,8,4), ("gf16",500,5,1024,1,1),
    ("gf16",1000,17,64,1,1), ("gf16",1000,17,1024,1,1),
    ("gf16",1000,33,64,1,1),
)


def fail(message: str) -> None:
    raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        fail(f"{path}: top-level JSON is not an object")
    return value


def parse_mapping(values: list[str], option: str) -> dict[str, Path]:
    output: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            fail(f"{option}: expected LABEL=PATH")
        label, path = value.split("=", 1)
        if not label or label in output:
            fail(f"{option}: duplicate/empty label {label!r}")
        output[label] = Path(path)
    return output


def require_commit(repository: Path, revision: str) -> None:
    result = subprocess.run(
        ["git", "cat-file", "-e", f"{revision}^{{commit}}"],
        cwd=repository, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, check=False,
    )
    if result.returncode:
        fail(f"core_git_sha is not a commit: {revision}")


def validate_executable(
    path: Path, source_hash: str, library: Path, expected_label: str,
    repository: Path, sanitizer: bool = False,
) -> dict[str, Any]:
    value = read_json(path)
    if value.get("schema") != "leopard2-c8-executable-v1":
        fail(f"{path}: wrong schema")
    expected_mode = "correctness"
    if value.get("mode") != expected_mode:
        fail(f"{path}: wrong mode")
    if value.get("backend_label") != expected_label:
        fail(f"{path}: wrong backend label")
    if value.get("profile") != "exact_high_prefix_v1_candidate" or \
            value.get("default_off") is not True:
        fail(f"{path}: profile/default-off claim changed")
    if value.get("source_sha256") != source_hash:
        fail(f"{path}: stale source binding")
    if value.get("library_sha256") != sha256(library):
        fail(f"{path}: library hash mismatch")
    revision = value.get("core_git_sha")
    if not isinstance(revision, str) or len(revision) != 40:
        fail(f"{path}: invalid core revision")
    require_commit(repository, revision)
    if value.get("correctness") != EXPECTED_CORRECTNESS:
        fail(f"{path}: correctness projection changed")
    if value.get("benchmarks") != []:
        fail(f"{path}: correctness artifact contains timing")
    if value.get("openmp_max_threads") != 1:
        fail(f"{path}: correctness did not force one OpenMP thread")
    if sanitizer:
        if value.get("sanitizer") != "asan-ubsan" or \
                value.get("allocation_tracking") != "disabled-for-sanitizer":
            fail(f"{path}: wrong sanitizer provenance")
    else:
        if value.get("sanitizer") != "none" or \
                value.get("allocation_tracking") != "global-new":
            fail(f"{path}: normal artifact has wrong instrumentation")
        runtime = value.get("runtime_backend")
        if expected_label != "auto" and runtime != expected_label:
            fail(f"{path}: forced backend did not execute")
    return value


def close_enough(observed: float, expected: float, scale: float = 1.) -> bool:
    return math.isfinite(observed) and abs(observed - expected) <= \
        max(5e-6, abs(scale) * 5e-6)


def summarize_samples(values: object, expected_count: int,
                      label: str) -> tuple[float, float]:
    if not isinstance(values, list) or len(values) != expected_count:
        fail(f"benchmark: {label} must contain {expected_count} samples")
    samples: list[float] = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or \
                not math.isfinite(value) or value <= 0:
            fail(f"benchmark: invalid {label} sample")
        samples.append(float(value))
    ordered = sorted(samples)
    median = ordered[len(ordered) // 2]
    deviations = sorted(abs(value - median) for value in samples)
    return median, deviations[len(deviations) // 2]


def validate_benchmark(
    path: Path, source_hash: str, library: Path, repository: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    value = read_json(path)
    if value.get("schema") != "leopard2-c8-executable-v1" or \
            value.get("mode") != "benchmark":
        fail("benchmark: wrong schema/mode")
    if value.get("backend_label") != "avx2" or \
            value.get("runtime_backend") != "avx2":
        fail("benchmark: timing did not execute forced AVX2")
    if value.get("sanitizer") != "none" or \
            value.get("allocation_tracking") != "global-new":
        fail("benchmark: timing was instrumented")
    if value.get("source_sha256") != source_hash or \
            value.get("library_sha256") != sha256(library):
        fail("benchmark: stale source/library binding")
    revision = value.get("core_git_sha")
    if not isinstance(revision, str):
        fail("benchmark: missing core revision")
    require_commit(repository, revision)
    if value.get("openmp_max_threads") != 1:
        fail("benchmark: OpenMP max threads was not one")
    affinity = value.get("affinity")
    if not isinstance(affinity, list) or len(affinity) != 1 or \
            not isinstance(affinity[0], int):
        fail("benchmark: process was not pinned to one CPU")
    rows = value.get("benchmarks")
    if not isinstance(rows, list) or len(rows) != len(EXPECTED_BENCHMARKS):
        fail("benchmark: incomplete matrix")
    observed_keys = []
    validated = []
    for row in rows:
        if not isinstance(row, dict):
            fail("benchmark: row is not an object")
        key = (row.get("field"), row.get("K"), row.get("R"), row.get("bytes"),
               row.get("batch"), row.get("reuse"))
        observed_keys.append(key)
        for name in ("exact_median_us", "exact_mad_us", "padded_median_us",
                     "padded_mad_us", "exact_setup_median_us",
                     "exact_setup_mad_us", "padded_setup_median_us",
                     "padded_setup_mad_us"):
            number = row.get(name)
            if not isinstance(number, (int, float)) or not math.isfinite(number) or \
                    number < 0 or ("median" in name and number == 0):
                fail(f"benchmark: invalid {name}")
        sample_fields = (
            ("exact_samples_us", 11, "exact_median_us", "exact_mad_us"),
            ("padded_samples_us", 11, "padded_median_us", "padded_mad_us"),
            ("exact_setup_samples_us", 7, "exact_setup_median_us",
             "exact_setup_mad_us"),
            ("padded_setup_samples_us", 7, "padded_setup_median_us",
             "padded_setup_mad_us"),
        )
        for samples_name, count, median_name, mad_name in sample_fields:
            median, mad = summarize_samples(row.get(samples_name), count,
                                            samples_name)
            if not close_enough(float(row[median_name]), median, median) or \
                    not close_enough(float(row[mad_name]), mad,
                                     max(median, mad)):
                fail(f"benchmark: {samples_name} summary mismatch")
        expected_ratio = row["padded_median_us"] / row["exact_median_us"]
        expected_gain = ((row["padded_median_us"] - row["padded_mad_us"]) /
                         (row["exact_median_us"] + row["exact_mad_us"]) - 1.)
        executions = row["batch"] * row["reuse"]
        amortized_exact = (row["exact_median_us"] + row["exact_mad_us"] +
                           (row["exact_setup_median_us"] +
                            row["exact_setup_mad_us"]) / executions)
        amortized_padded = (row["padded_median_us"] - row["padded_mad_us"] +
                            max(0., row["padded_setup_median_us"] -
                                row["padded_setup_mad_us"]) / executions)
        amortized_gain = amortized_padded / amortized_exact - 1.
        if not close_enough(float(row.get("ratio", math.nan)), expected_ratio,
                            expected_ratio):
            fail("benchmark: forged ratio")
        if not close_enough(float(row.get("credible_gain", math.nan)), expected_gain,
                            expected_gain):
            fail("benchmark: forged credible gain")
        width = 1 if row["field"] == "gf8" else 2
        expected_table = row["K"] * row["R"] * width
        expected_traffic = row["R"] * (1 + 3 * (row["K"] - 1)) * row["bytes"]
        if row.get("exact_table_bytes") != expected_table or \
                row.get("exact_scratch_bytes") != 0 or \
                row.get("exact_logical_bytes") != expected_traffic:
            fail("benchmark: exact memory accounting changed")
        if not isinstance(row.get("padded_scratch_bytes"), int) or \
                row["padded_scratch_bytes"] <= 0:
            fail("benchmark: invalid padded scratch")
        validated.append({**row, "recomputed_ratio": expected_ratio,
                          "recomputed_credible_gain": expected_gain,
                          "credible_amortized_gain": amortized_gain})
    if tuple(observed_keys) != EXPECTED_BENCHMARKS:
        fail("benchmark: geometry/order changed")
    if value.get("correctness") is not None:
        fail("benchmark: timing artifact unexpectedly contains correctness")
    return value, validated


def validate_isolation(path: Path, timing_cpu: int) -> dict[str, Any]:
    value = read_json(path)
    if value.get("schema") != "leopard2-c8-isolation-v1" or \
            value.get("exit_code") != 0:
        fail("isolation: wrong schema or failed child")
    if value.get("timing_cpu") != timing_cpu or \
            value.get("sibling_cpu") == timing_cpu:
        fail("isolation: timing CPU/sibling mismatch")
    if value.get("sibling_idle") is not True or \
            value.get("sibling_nonidle_jiffies") != 0:
        fail("isolation: sibling was not idle")
    if value.get("environment") != {
            "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1"}:
        fail("isolation: threading environment changed")
    elapsed = value.get("elapsed_seconds")
    if not isinstance(elapsed, (int, float)) or not math.isfinite(elapsed) or elapsed <= 0:
        fail("isolation: invalid duration")
    command = value.get("command")
    if not isinstance(command, list) or "--mode" not in command or \
            "benchmark" not in command or "--backend-label" not in command or \
            "avx2" not in command:
        fail("isolation: wrong benchmark child")
    delta = value.get("delta", {}).get(str(timing_cpu), {})
    if not isinstance(delta, dict) or sum(
            int(delta.get(name, 0)) for name in
            ("user", "nice", "system", "irq", "softirq", "steal")) <= 0:
        fail("isolation: timing CPU recorded no execution")
    return value


def validate_algebra(directory: Path, source: Path) -> dict[str, Any]:
    manifest = read_json(directory / "algebra-manifest.json")
    if manifest.get("schema") != "leopard2-c8-manifest-v1":
        fail("algebra: wrong manifest schema")
    for name, expected in manifest.get("files", {}).items():
        if sha256(directory / name) != expected:
            fail(f"algebra: hash mismatch for {name}")
    value = read_json(directory / "algebra.json")
    if value.get("schema") != "leopard2-c8-algebra-v1" or \
            value.get("source_sha256") != sha256(source):
        fail("algebra: stale source/schema")
    profile = value.get("profile", {})
    if profile != {
            "name": "exact_high_prefix_v1_candidate",
            "status": "unserialized_default_off_candidate",
            "coordinate_map_version": 1,
            "polynomial_degree_bound": "degree < K",
            "wire_order": ["parity evaluations 0..R-1",
                           "systematic evaluations R..R+K-1"],
    }:
        fail("algebra: profile identity changed")
    gf4 = value.get("gf4", {})
    required_gf4 = {
        "public_geometries": 120,
        "k_coordinate_subsets": 131038,
        "solver_vectors": 57670,
        "parity_symbols": 403993,
        "full_parent_wire_identity_geometries": 10,
        "nonidentity_wire_witness_geometries": 75,
    }
    if gf4 != required_gf4:
        fail("algebra: GF4 projection changed")
    model = value.get("gf8_model", {})
    if model.get("public_pairs") != 32640 or \
            model.get("method_rows") != 228480 or \
            model.get("winner_counts") != {
                "padded_t_lch": 29359,
                "precomputed_dense_r_solve": 3281,
            }:
        fail("algebra: incomplete GF8 model")
    differential = value.get("differential", [])
    if not isinstance(differential, list) or \
            not all(isinstance(row, dict) for row in differential):
        fail("algebra: invalid GF8/GF16 differential set")
    projection = tuple(
        (row.get("field"), row.get("K"), row.get("R"), row.get("vectors"),
         row.get("coefficients"), row.get("legacy_wire_identity"),
         row.get("dense_constraint_checked"), row.get("schur_checked"))
        for row in differential)
    if projection != EXPECTED_DIFFERENTIAL:
        fail("algebra: GF8/GF16 differential projection changed")
    return value


def validate_default_off(repository: Path) -> None:
    root_cmake = (repository / "CMakeLists.txt").read_text(encoding="utf-8")
    forbidden = ("c8_exact_high.cpp", "non_power_of_two/c8")
    if any(token in root_cmake for token in forbidden):
        fail("default root CMake graph references C8")
    header = (repository / "leopard2.h").read_text(encoding="utf-8")
    if "EXACT_EXPERIMENTAL_V1 reserves a distinct code identity" not in header or \
            "currently returns" not in header:
        fail("stable exact-profile unsupported contract changed")


def canonical_write(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n",
                    encoding="utf-8")


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    backends = parse_mapping(arguments.backend, "--backend")
    libraries = parse_mapping(arguments.library, "--library")
    if tuple(sorted(backends)) != tuple(sorted(EXPECTED_BACKENDS)) or \
            set(backends) != set(libraries):
        fail("backend/library labels must be auto, scalar, ssse3, avx2")
    source_hash = sha256(arguments.source)
    algebra = validate_algebra(arguments.algebra_dir, arguments.algebra_source)
    validated_backends = {}
    core_revisions = set()
    for label in EXPECTED_BACKENDS:
        value = validate_executable(backends[label], source_hash, libraries[label],
                                    label, arguments.repository)
        core_revisions.add(value["core_git_sha"])
        validated_backends[label] = {
            "path": str(backends[label]), "sha256": sha256(backends[label]),
            "library_sha256": sha256(libraries[label]),
            "runtime_backend": value["runtime_backend"],
        }
    sanitizer = validate_executable(
        arguments.sanitizer, source_hash, arguments.sanitizer_library,
        "asan-ubsan", arguments.repository, sanitizer=True)
    core_revisions.add(sanitizer["core_git_sha"])
    benchmark, rows = validate_benchmark(
        arguments.benchmark, source_hash, arguments.benchmark_library,
        arguments.repository)
    core_revisions.add(benchmark["core_git_sha"])
    isolation = validate_isolation(arguments.isolation, benchmark["affinity"][0])
    if len(core_revisions) != 1:
        fail("artifacts use different core revisions")
    validate_default_off(arguments.repository)

    positive = [row for row in rows if row["recomputed_credible_gain"] >= .10]
    amortized_positive = [row for row in rows
                          if row["credible_amortized_gain"] >= .10]
    regressions = [row for row in rows if row["recomputed_credible_gain"] < -.02]
    # The GF16 R=2..4 reused tiny-shard region is real, but an encoder-only,
    # unserialized profile cannot be promoted before C9/C10 establish decode
    # behavior and a deterministic region dispatcher.  The broad candidate is
    # rejected; the narrow region remains an explicit inconclusive follow-up.
    disposition = "inconclusive_no_promotion"
    reason = ("dense exact-R execution has a contiguous GF16 tiny-R win, but "
              "setup amortization, GF8 losses, the R>=5 crossover, missing "
              "exact-profile decoding/serialization, and absent C10 dispatch "
              "evidence prevent production promotion")
    output = {
        "schema": "leopard2-c8-checkpoint-v1",
        "profile": "exact_high_prefix_v1_candidate",
        "profile_status": "unserialized_default_off_candidate",
        "core_git_sha": next(iter(core_revisions)),
        "source_sha256": source_hash,
        "checkpoint_source_sha256": sha256(Path(__file__)),
        "algebra": {
            "manifest_sha256": sha256(arguments.algebra_dir / "algebra-manifest.json"),
            "gf4": algebra["gf4"], "gf8_model": algebra["gf8_model"],
            "differential_cases": len(algebra["differential"]),
        },
        "backends": validated_backends,
        "sanitizer": {
            "path": str(arguments.sanitizer),
            "sha256": sha256(arguments.sanitizer),
            "library_sha256": sha256(arguments.sanitizer_library),
            "mode": sanitizer["sanitizer"],
        },
        "benchmark": {
            "path": str(arguments.benchmark),
            "sha256": sha256(arguments.benchmark),
            "library_sha256": sha256(arguments.benchmark_library),
            "isolation_sha256": sha256(arguments.isolation),
            "sibling_cpu": isolation["sibling_cpu"],
            "sibling_nonidle_jiffies": isolation["sibling_nonidle_jiffies"],
            "timing_cpu": benchmark["affinity"][0],
            "cells": len(rows),
            "cells_at_least_10_percent": len(positive),
            "amortized_cells_at_least_10_percent": len(amortized_positive),
            "cells_regressing_over_2_percent": len(regressions),
            "credible_gain_minimum": min(row["recomputed_credible_gain"] for row in rows),
            "credible_gain_median": sorted(
                row["recomputed_credible_gain"] for row in rows)[len(rows) // 2],
            "credible_gain_maximum": max(row["recomputed_credible_gain"] for row in rows),
            "rows": rows,
        },
        "disposition": disposition,
        "reason": reason,
        "promoted": [],
        "rejected": [
            "dense exact-R as a universal encoder",
            "factorized R-by-R execution",
            "Newton execution",
            "dyadic Schur execution",
        ],
        "retained_inconclusive": [
            "precomputed dense GF16 execution for R=2..4 with codec reuse",
            "partial/transposed exact-R LCH schedule",
            "Tang-Han epsilon adaptation to the suffix coordinate set",
        ],
        "production_changes": False,
        "public_profile_changes": False,
        "default_build_changes": False,
    }
    canonical_write(arguments.output, output)
    return output


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--algebra-dir", type=Path, required=True)
    result.add_argument("--algebra-source", type=Path, required=True)
    result.add_argument("--backend", action="append", default=[], required=True)
    result.add_argument("--library", action="append", default=[], required=True)
    result.add_argument("--sanitizer", type=Path, required=True)
    result.add_argument("--sanitizer-library", type=Path, required=True)
    result.add_argument("--benchmark", type=Path, required=True)
    result.add_argument("--benchmark-library", type=Path, required=True)
    result.add_argument("--isolation", type=Path, required=True)
    result.add_argument("--source", type=Path, required=True)
    result.add_argument("--repository", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    return result


def main() -> int:
    run(parser().parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
