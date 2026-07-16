#!/usr/bin/env python3
"""Validate and merge Leopard2 C5 dyadic-block C++ evidence.

The merger is intentionally fail-closed: raw results are accepted only when
their source/library bindings match the supplied files, every backend has the
same deterministic correctness projection, the sanitizer ran the same matrix,
and the timing artifact contains the declared ABBA matrix.  The output is a
deterministic checkpoint over those retained inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import statistics
import subprocess
from typing import Any


SCHEMA = "leopard2-c5-cpp-v1"
EXPECTED_BACKENDS = ("auto", "scalar", "ssse3", "avx2")
EXPECTED_CORRECTNESS_CASES = 144
EXPECTED_DIRECT_CHECKS = 388
EXPECTED_FACTOR_CHECKS = 1_260_048
BENCHMARK_Q = (33, 65, 129, 257, 513)


class CheckpointError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CheckpointError(message)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: pathlib.Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CheckpointError(f"cannot read {path}: {error}") from error
    require(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def parse_binding(values: list[str], option: str) -> dict[str, pathlib.Path]:
    result: dict[str, pathlib.Path] = {}
    for value in values:
        require("=" in value, f"{option} requires LABEL=PATH")
        label, raw_path = value.split("=", 1)
        require(label and raw_path, f"{option} requires LABEL=PATH")
        require(label not in result, f"duplicate {option} label: {label}")
        result[label] = pathlib.Path(raw_path)
    return result


def validate_common(
    path: pathlib.Path,
    payload: dict[str, Any],
    source_hash: str,
    expected_mode: str,
) -> None:
    require(payload.get("schema_version") == SCHEMA,
            f"{path}: wrong schema")
    require(payload.get("status") == "pass", f"{path}: status is not pass")
    require(payload.get("wire_identity") == "existing padded dyadic parent",
            f"{path}: wrong wire identity")
    require(payload.get("exact_profile_implemented") is False,
            f"{path}: exact profile claim changed")
    require(payload.get("default_build_integration") is False,
            f"{path}: default-build claim changed")
    require(payload.get("mode") == expected_mode,
            f"{path}: expected mode {expected_mode}")
    require(payload.get("kernel_quantum_bytes") == 64,
            f"{path}: kernel quantum changed")
    binding = payload.get("build_binding")
    require(isinstance(binding, dict), f"{path}: missing build binding")
    require(binding.get("source_sha256") == source_hash,
            f"{path}: source hash mismatch")
    require(isinstance(binding.get("core_git_sha"), str) and
            len(binding["core_git_sha"]) == 40,
            f"{path}: invalid core git SHA")
    require(isinstance(binding.get("linked_library_sha256"), str) and
            len(binding["linked_library_sha256"]) == 64,
            f"{path}: invalid library hash")
    require(len(payload.get("correctness_cases", [])) ==
            EXPECTED_CORRECTNESS_CASES,
            f"{path}: correctness matrix size changed")
    require(payload.get("direct_symbol_checks") == EXPECTED_DIRECT_CHECKS,
            f"{path}: direct check count changed")
    require(payload.get("factor_checks") == EXPECTED_FACTOR_CHECKS,
            f"{path}: factor check count changed")
    require(isinstance(payload.get("correctness_digest_fnv1a64"), str),
            f"{path}: missing correctness digest")

    saw_tail = False
    saw_gf8 = False
    saw_gf16 = False
    saw_large_gf16 = False
    for case in payload["correctness_cases"]:
        require(case["field_bits"] in (8, 16), f"{path}: bad field")
        saw_gf8 |= case["field_bits"] == 8
        saw_gf16 |= case["field_bits"] == 16
        saw_large_gf16 |= case["field_bits"] == 16 and case["q"] >= 8191
        processed = (case["valid_bytes"] + 63) // 64 * 64
        require(case["processed_bytes"] == processed,
                f"{path}: tail rounding changed")
        saw_tail |= case["valid_bytes"] != processed
        require(case["compared_bytes"] == case["parent"] * processed,
                f"{path}: compared-byte count is inconsistent")
        require(case["q"] <= case["parent"], f"{path}: q exceeds parent")
        require(case["parent"] & (case["parent"] - 1) == 0,
                f"{path}: non-dyadic parent")
        require(case["candidate_scratch_bytes"] >= 0 and
                case["padded_scratch_bytes"] >= 0,
                f"{path}: invalid scratch accounting")
        require(case["candidate_traffic_bytes"] > 0 and
                case["padded_traffic_bytes"] > 0,
                f"{path}: invalid traffic accounting")
    require(saw_tail and saw_gf8 and saw_gf16 and saw_large_gf16,
            f"{path}: field/tail coverage incomplete")


def correctness_projection(payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "digest": payload["correctness_digest_fnv1a64"],
        "direct": payload["direct_symbol_checks"],
        "factors": payload["factor_checks"],
        "cases": payload["correctness_cases"],
    }


def validate_library_binding(
    label: str,
    result_path: pathlib.Path,
    payload: dict[str, Any],
    library_path: pathlib.Path,
) -> None:
    require(result_path.exists(), f"missing result: {result_path}")
    require(library_path.exists(), f"missing library: {library_path}")
    require(payload.get("requested_backend") == label,
            f"{result_path}: backend label mismatch")
    require(payload["build_binding"]["linked_library_sha256"] ==
            sha256(library_path),
            f"{result_path}: linked library hash mismatch")


def validate_benchmarks(path: pathlib.Path, payload: dict[str, Any]) -> None:
    rows = payload.get("benchmarks")
    require(isinstance(rows, list), f"{path}: benchmarks are absent")
    require(len(rows) == 54, f"{path}: expected 54 benchmark cells")
    seen_q = {row["q"] for row in rows}
    require(seen_q == set(BENCHMARK_Q), f"{path}: benchmark q matrix changed")
    require({row["valid_bytes"] for row in rows} >=
            {64, 1024, 65536, 1024 * 1024},
            f"{path}: byte matrix incomplete")
    require({row["batch"] for row in rows} >= {1, 8},
            f"{path}: batch matrix incomplete")
    require({row["reuse"] for row in rows} >= {1, 8},
            f"{path}: reuse matrix incomplete")
    for row in rows:
        require(row["samples"] >= 7 and row["samples"] % 2 == 1,
                f"{path}: invalid ABBA sample count")
        for key in (
            "setup_median_us", "candidate_median_us", "padded_median_us",
            "candidate_scratch_bytes_per_stripe",
            "padded_scratch_bytes_per_stripe", "resident_batch_bytes",
            "candidate_traffic_bytes_per_execution",
            "padded_traffic_bytes_per_execution",
        ):
            require(math.isfinite(row[key]) and row[key] >= 0,
                    f"{path}: invalid {key}")
        require(row["candidate_median_us"] > 0 and
                row["padded_median_us"] > 0,
                f"{path}: zero execution time")


def summarize_benchmarks(rows: list[dict[str, Any]]) -> dict[str, Any]:
    per_q: list[dict[str, Any]] = []
    qualified_q: list[int] = []
    for q in BENCHMARK_Q:
        selected = [row for row in rows if row["q"] == q]
        gains = [row["credible_gain_percent"] for row in selected]
        credible_wins = [row for row in selected
                         if row["credible_gain_percent"] >= 10.0]
        byte_sizes = {row["valid_bytes"] for row in credible_wins}
        meaningful = len(credible_wins) >= 3 and len(byte_sizes) >= 2 and all(
            row["credible_gain_percent"] >= -2.0 for row in selected)
        if meaningful:
            qualified_q.append(q)
        per_q.append({
            "q": q,
            "cells": len(selected),
            "credible_gain_percent_min": min(gains),
            "credible_gain_percent_median": statistics.median(gains),
            "credible_gain_percent_max": max(gains),
            "credible_wins_at_least_10_percent": len(credible_wins),
            "meaningful_region_gate": meaningful,
            "modeled_traffic_reduction_percent_min": min(
                100.0 * (1.0 - row["candidate_traffic_bytes_per_execution"] /
                         row["padded_traffic_bytes_per_execution"])
                for row in selected),
            "modeled_traffic_reduction_percent_max": max(
                100.0 * (1.0 - row["candidate_traffic_bytes_per_execution"] /
                         row["padded_traffic_bytes_per_execution"])
                for row in selected),
        })
    return {
        "promotion_threshold_percent": 10.0,
        "neighbor_regression_limit_percent": 2.0,
        "qualified_q": qualified_q,
        "production_promotion": False,
        "experimental_followup": bool(qualified_q),
        "disposition": (
            "inconclusive-followup: a bounded q region passed the measured gate; "
            "no default integration was performed"
            if qualified_q else
            "reject: no bounded q region passed the measured 10% gate without "
            "neighboring regressions"
        ),
        "per_q": per_q,
    }


def git_status(repository: pathlib.Path) -> str:
    completed = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return completed.stdout.strip()


def run(args: argparse.Namespace) -> dict[str, Any]:
    source = pathlib.Path(args.source)
    require(source.exists(), f"missing source: {source}")
    source_hash = sha256(source)
    backend_results = parse_binding(args.backend, "--backend")
    backend_libraries = parse_binding(args.library, "--library")
    require(set(backend_results) == set(EXPECTED_BACKENDS),
            "backend results must be auto, scalar, ssse3, and avx2")
    require(set(backend_libraries) == set(EXPECTED_BACKENDS),
            "backend libraries must be auto, scalar, ssse3, and avx2")

    loaded_backends: dict[str, dict[str, Any]] = {}
    projections: list[dict[str, Any]] = []
    raw_hashes: dict[str, str] = {}
    core_shas: set[str] = set()
    for label in EXPECTED_BACKENDS:
        path = backend_results[label]
        payload = load_json(path)
        validate_common(path, payload, source_hash, "backend")
        validate_library_binding(label, path, payload, backend_libraries[label])
        require(len(payload["benchmarks"]) == len(BENCHMARK_Q),
                f"{path}: backend timing matrix changed")
        loaded_backends[label] = payload
        projections.append(correctness_projection(payload))
        raw_hashes[f"backend_{label}"] = sha256(path)
        core_shas.add(payload["build_binding"]["core_git_sha"])
    require(all(item == projections[0] for item in projections[1:]),
            "backend correctness projections differ")

    sanitizer_path = pathlib.Path(args.sanitizer)
    sanitizer_library = pathlib.Path(args.sanitizer_library)
    sanitizer = load_json(sanitizer_path)
    validate_common(sanitizer_path, sanitizer, source_hash, "correctness")
    require(sanitizer["requested_backend"] == "asan-ubsan",
            "sanitizer label mismatch")
    require(sanitizer["build_binding"]["sanitizer_mode"] == "asan-ubsan",
            "sanitizer mode mismatch")
    require(sanitizer_library.exists(), "sanitizer library is missing")
    require(sanitizer["build_binding"]["linked_library_sha256"] ==
            sha256(sanitizer_library), "sanitizer library hash mismatch")
    require(correctness_projection(sanitizer) == projections[0],
            "sanitizer correctness projection differs")
    raw_hashes["sanitizer"] = sha256(sanitizer_path)
    core_shas.add(sanitizer["build_binding"]["core_git_sha"])

    benchmark_path = pathlib.Path(args.benchmark)
    benchmark = load_json(benchmark_path)
    validate_common(benchmark_path, benchmark, source_hash, "all")
    require(benchmark["requested_backend"] == "avx2",
            "benchmark must use the fastest available AVX2 library")
    require(benchmark["build_binding"]["linked_library_sha256"] ==
            sha256(backend_libraries["avx2"]),
            "benchmark AVX2 library hash mismatch")
    require(correctness_projection(benchmark) == projections[0],
            "benchmark correctness projection differs")
    validate_benchmarks(benchmark_path, benchmark)
    raw_hashes["benchmark"] = sha256(benchmark_path)
    core_shas.add(benchmark["build_binding"]["core_git_sha"])
    require(len(core_shas) == 1, "raw artifacts use different core revisions")

    current_head = git_status(pathlib.Path(args.repository))
    return {
        "schema_version": "leopard2-c5-checkpoint-v1",
        "status": "pass",
        "wire_identity": "existing padded dyadic parent",
        "exact_profile_implemented": False,
        "default_build_integration": False,
        "build_binding": {
            "experiment_source_sha256": source_hash,
            "core_baseline_git_sha": next(iter(core_shas)),
            "checkpoint_generation_head": current_head,
        },
        "correctness": {
            "backend_variants": list(EXPECTED_BACKENDS),
            "sanitizers": ["address", "undefined"],
            "cases_per_variant": EXPECTED_CORRECTNESS_CASES,
            "backend_case_executions":
                EXPECTED_CORRECTNESS_CASES * len(EXPECTED_BACKENDS),
            "sanitizer_case_executions": EXPECTED_CORRECTNESS_CASES,
            "direct_symbol_checks_per_variant": EXPECTED_DIRECT_CHECKS,
            "factor_checks_per_variant": EXPECTED_FACTOR_CHECKS,
            "digest": projections[0]["digest"],
            "tail_policy": benchmark["tail_policy"],
        },
        "benchmark": {
            "cells": len(benchmark["benchmarks"]),
            "setup_reported_separately": True,
            "execution_order": "alternating ABBA",
            "backend": benchmark["runtime_backend"],
            "thread_count": 1,
            "affinity_cpus":
                benchmark["runtime_environment"]["process_affinity_cpus"],
            "shard_bytes": sorted({row["valid_bytes"]
                                    for row in benchmark["benchmarks"]}),
            "batch": sorted({row["batch"] for row in benchmark["benchmarks"]}),
            "reuse": sorted({row["reuse"] for row in benchmark["benchmarks"]}),
            "setup_median_us_range": [
                min(row["setup_median_us"] for row in benchmark["benchmarks"]),
                max(row["setup_median_us"] for row in benchmark["benchmarks"]),
            ],
        },
        "promotion": summarize_benchmarks(benchmark["benchmarks"]),
        "limitations": [
            "The candidate operates on existing parent LCH coefficients; it is not a systematic encoder API.",
            "The only fused route is q=2^a+1; general multi-block factors use scalar multiplication.",
            "Tail bytes are zero-padded to the legacy 64-byte kernel quantum.",
            "No exact-size coordinate profile or serialized code identity is introduced.",
            "Timing is authoritative only for the recorded pinned host/backend.",
        ],
        "artifacts_sha256": raw_hashes,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", action="append", default=[],
                        help="LABEL=raw-result.json")
    parser.add_argument("--library", action="append", default=[],
                        help="LABEL=linked-libleopard.a")
    parser.add_argument("--sanitizer", required=True)
    parser.add_argument("--sanitizer-library", required=True)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--repository", default=".")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> int:
    try:
        args = parse_args()
        checkpoint = run(args)
        output = pathlib.Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(checkpoint, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        return 0
    except (CheckpointError, subprocess.CalledProcessError) as error:
        print(f"FAIL c5 checkpoint: {error}", file=__import__("sys").stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
