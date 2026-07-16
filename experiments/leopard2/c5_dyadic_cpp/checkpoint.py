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
import itertools
import json
import math
import pathlib
import re
import statistics
import subprocess
from typing import Any


SCHEMA = "leopard2-c5-cpp-v1"
EXPECTED_BACKENDS = ("auto", "scalar", "ssse3", "avx2")
EXPECTED_CORRECTNESS_CASES = 144
EXPECTED_DIRECT_CHECKS = 388
EXPECTED_FACTOR_CHECKS = 1_260_048
EXPECTED_CORRECTNESS_DIGEST = "0x56af7bbc2fb6a888"
CORE_BASELINE_GIT_SHA = "37e852774bce6f9effb1acf1fcde99a758ecfe6e"
TAIL_POLICY = "zero-pad to kernel quantum; compare valid and padded bytes"
NORMAL_SANITIZER_MODE = "none"
SANITIZED_MODE = "asan-ubsan"
RETAINED_POINTER_BYTES = 8
GUARD_BYTES_PER_SHARD = 64
CORRECTNESS_GF8_Q = (
    1, 3, 5, 7, 9, 15, 17, 31, 33, 63, 65, 127, 129, 191, 255,
)
CORRECTNESS_GF16_Q = (257, 259, 511, 513, 1000, 1025, 4097, 8191)
CORRECTNESS_TAIL_BYTES = (1, 17, 63, 64, 65, 129)
CORRECTNESS_BOUNDARY_BYTES = (4095, 4096, 4097)
BENCHMARK_GEOMETRIES = (
    (8, 33), (8, 65), (8, 129), (16, 257), (16, 513),
)
BENCHMARK_Q = tuple(q for _field_bits, q in BENCHMARK_GEOMETRIES)

# The C++ result writer uses fixed precision with six digits after the decimal
# point.  A displayed primitive or derived value can therefore differ from its
# in-memory value by at most half a unit in the last place.  Derived-field
# validation below propagates this interval instead of relying on a broad,
# empirically chosen relative tolerance.
JSON_DECIMAL_HALF_ULP = 0.5e-6
FLOAT_INTERVAL_SLOP = 1e-12


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


def require_int(
    record: dict[str, Any],
    key: str,
    path: pathlib.Path,
    minimum: int = 0,
) -> int:
    value = record.get(key)
    require(isinstance(value, int) and not isinstance(value, bool),
            f"{path}: {key} is not an integer")
    require(value >= minimum, f"{path}: invalid {key}")
    return value


def require_finite_nonnegative(
    record: dict[str, Any],
    key: str,
    path: pathlib.Path,
) -> float:
    value = record.get(key)
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{path}: {key} is not numeric")
    result = float(value)
    require(math.isfinite(result) and result >= 0.0,
            f"{path}: invalid {key}")
    return result


def ceil_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def canonical_binary_blocks(q: int) -> tuple[tuple[int, int], ...]:
    require(q >= 1, "binary prefix must be nonempty")
    result: list[tuple[int, int]] = []
    offset = 0
    remaining = q
    while remaining:
        block_size = 1 << (remaining.bit_length() - 1)
        require(offset & (block_size - 1) == 0,
                "internal binary block is not aligned")
        result.append((offset, block_size))
        offset += block_size
        remaining -= block_size
    require(offset == q, "internal binary blocks do not cover prefix")
    return tuple(result)


def active_jobs_for_block(offset: int, size: int, parent: int) -> int:
    """Count nonzero LCH block factors in the legacy coordinate order.

    A non-base block factor contains the subspace polynomial for the highest
    set bit of ``offset``.  Its roots are exactly the coordinate prefix below
    that power of two, so all aligned shifts at or above it are active.
    """
    require(size >= 1 and parent % size == 0,
            "internal job geometry is invalid")
    if offset == 0:
        return parent // size
    first_active_shift = 1 << (offset.bit_length() - 1)
    require(first_active_shift % size == 0 and first_active_shift < parent,
            "internal active-job boundary is invalid")
    return (parent - first_active_shift) // size


def expected_execution_accounting(q: int, processed_bytes: int) -> dict[str, int]:
    """Derive every non-timing execution counter written by c5_dyadic.cpp."""
    require(processed_bytes >= 1 and processed_bytes % 64 == 0,
            "processed bytes must use the retained kernel quantum")
    parent = ceil_power_of_two(q)
    blocks = canonical_binary_blocks(q)
    largest_block = blocks[0][1]
    fused_power_plus_one = q == largest_block + 1
    largest_tail_block = blocks[1][1] if len(blocks) > 1 else 0

    jobs = sum(active_jobs_for_block(offset, size, parent)
               for offset, size in blocks)
    candidate_scratch = largest_block * RETAINED_POINTER_BYTES
    if not fused_power_plus_one and largest_tail_block > 1:
        candidate_scratch += largest_tail_block * (
            processed_bytes + GUARD_BYTES_PER_SHARD)
    padded_scratch = parent * RETAINED_POINTER_BYTES

    candidate_traffic = 0
    for block_index, (offset, size) in enumerate(blocks):
        active_jobs = active_jobs_for_block(offset, size, parent)
        base = block_index == 0
        if fused_power_plus_one and not base:
            continue
        shard_bytes = size * processed_bytes
        butterfly_bytes = 2 * processed_bytes * size * (size.bit_length() - 1)
        if base:
            candidate_traffic += active_jobs * (
                2 * shard_bytes + butterfly_bytes)
            if fused_power_plus_one:
                # The upper base coset injects the one-shard tail: one input
                # read plus one output read and write.
                require(active_jobs == 2,
                        "internal fused base must have two output cosets")
                candidate_traffic += 3 * processed_bytes
        elif size > 1:
            candidate_traffic += active_jobs * (
                5 * shard_bytes + butterfly_bytes)
        else:
            candidate_traffic += active_jobs * 3 * shard_bytes

    padded_traffic = processed_bytes * (
        q + parent + 2 * parent * (parent.bit_length() - 1))
    return {
        "jobs": jobs,
        "candidate_scratch_bytes": candidate_scratch,
        "padded_scratch_bytes": padded_scratch,
        "candidate_traffic_bytes": candidate_traffic,
        "padded_traffic_bytes": padded_traffic,
    }


def expected_resident_batch_bytes(
    q: int,
    processed_bytes: int,
    batch: int,
) -> int:
    parent = ceil_power_of_two(q)
    accounting = expected_execution_accounting(q, processed_bytes)
    shard_storage = (q + parent) * (
        processed_bytes + GUARD_BYTES_PER_SHARD)
    return batch * (
        shard_storage + accounting["candidate_scratch_bytes"] +
        accounting["padded_scratch_bytes"])


def expected_correctness_cells() -> frozenset[tuple[int, int, int]]:
    result = {
        (8, q, valid_bytes)
        for q in CORRECTNESS_GF8_Q
        for valid_bytes in CORRECTNESS_TAIL_BYTES
    }
    result.update({
        (16, q, valid_bytes)
        for q in CORRECTNESS_GF16_Q
        for valid_bytes in CORRECTNESS_TAIL_BYTES
    })
    for valid_bytes in CORRECTNESS_BOUNDARY_BYTES:
        result.add((8, 129, valid_bytes))
        result.add((16, 257, valid_bytes))
    require(len(result) == EXPECTED_CORRECTNESS_CASES,
            "internal correctness-cell definition is inconsistent")
    return frozenset(result)


def expected_benchmark_cells(
    backend_only: bool = False,
) -> frozenset[tuple[int, int, int, int, int]]:
    if backend_only:
        return frozenset(
            (field_bits, q, 1024, 1, 8)
            for field_bits, q in BENCHMARK_GEOMETRIES
        )

    result: set[tuple[int, int, int, int, int]] = set()
    for field_bits, q in BENCHMARK_GEOMETRIES:
        for valid_bytes in (64, 1024):
            for batch in (1, 8):
                for reuse in (1, 8):
                    result.add((field_bits, q, valid_bytes, batch, reuse))
        for reuse in (1, 8):
            result.add((field_bits, q, 65536, 1, reuse))
    for field_bits, q in BENCHMARK_GEOMETRIES[:-1]:
        result.add((field_bits, q, 1024 * 1024, 1, 1))
    require(len(result) == 54,
            "internal benchmark-cell definition is inconsistent")
    return frozenset(result)


def parse_binding(values: list[str], option: str) -> dict[str, pathlib.Path]:
    result: dict[str, pathlib.Path] = {}
    for value in values:
        require("=" in value, f"{option} requires LABEL=PATH")
        label, raw_path = value.split("=", 1)
        require(label and raw_path, f"{option} requires LABEL=PATH")
        require(label not in result, f"duplicate {option} label: {label}")
        result[label] = pathlib.Path(raw_path)
    return result


def expected_direct_checks(field_bits: int, q: int, valid_bytes: int) -> int:
    if valid_bytes != 64 and not (
        valid_bytes == 1 and q in (129, 257)
    ):
        return 0
    points = 1
    parent = ceil_power_of_two(q)
    if parent > 1:
        points += 1
    if parent > 2:
        points += 1
    if parent > 3:
        points += 1
    # Both retained fields sample four lanes.
    require(field_bits in (8, 16), "internal direct-check field is invalid")
    return points * 4


def expected_factor_checks(q: int) -> int:
    parent = ceil_power_of_two(q)
    result = 0
    for offset, block_size in canonical_binary_blocks(q):
        # One normalizer comparison per canonical block, then one subspace
        # comparison per set bit of the block offset on every output coset.
        result += 1 + (parent // block_size) * offset.bit_count()
    return result


def validate_correctness_cases(
    path: pathlib.Path,
    payload: dict[str, Any],
) -> None:
    cases = payload.get("correctness_cases")
    require(isinstance(cases, list), f"{path}: correctness cases are absent")
    expected = expected_correctness_cells()
    require(len(cases) == len(expected),
            f"{path}: correctness matrix size changed")
    seen: set[tuple[int, int, int]] = set()
    direct_total = 0
    factor_total = 0
    for index, case in enumerate(cases):
        require(isinstance(case, dict),
                f"{path}: correctness case {index} is not an object")
        field_bits = require_int(case, "field_bits", path, 1)
        q = require_int(case, "q", path, 1)
        valid_bytes = require_int(case, "valid_bytes", path, 1)
        key = (field_bits, q, valid_bytes)
        require(key not in seen,
                f"{path}: duplicate correctness cell {key}")
        require(key in expected,
                f"{path}: unexpected correctness cell {key}")
        seen.add(key)

        parent = require_int(case, "parent", path, 1)
        processed_bytes = require_int(case, "processed_bytes", path, 1)
        require(parent == ceil_power_of_two(q),
                f"{path}: wrong parent for correctness cell {key}")
        require(processed_bytes == (valid_bytes + 63) // 64 * 64,
                f"{path}: wrong tail rounding for correctness cell {key}")
        require(require_int(case, "compared_bytes", path, 1) ==
                parent * processed_bytes,
                f"{path}: wrong compared-byte count for correctness cell {key}")
        require(require_int(case, "blocks", path, 1) == q.bit_count(),
                f"{path}: wrong binary-block count for correctness cell {key}")
        accounting = expected_execution_accounting(q, processed_bytes)
        require(require_int(case, "jobs", path, 1) == accounting["jobs"],
                f"{path}: wrong active-job count for correctness cell {key}")
        factor_checks = require_int(case, "factor_checks", path, 1)
        require(factor_checks == expected_factor_checks(q),
                f"{path}: wrong factor-check count for correctness cell {key}")
        factor_total += factor_checks
        direct = require_int(case, "direct_symbol_checks", path)
        require(direct == expected_direct_checks(field_bits, q, valid_bytes),
                f"{path}: wrong direct-check count for correctness cell {key}")
        direct_total += direct
        for accounting_key in (
            "candidate_scratch_bytes", "padded_scratch_bytes",
            "candidate_traffic_bytes", "padded_traffic_bytes",
        ):
            minimum = 1 if "traffic" in accounting_key else 0
            require(require_int(case, accounting_key, path, minimum) ==
                    accounting[accounting_key],
                    f"{path}: wrong {accounting_key} for correctness cell "
                    f"{key}")

    require(seen == expected, f"{path}: correctness matrix is incomplete")
    require(direct_total == EXPECTED_DIRECT_CHECKS,
            f"{path}: per-cell direct checks do not match the declared total")
    require(factor_total == EXPECTED_FACTOR_CHECKS,
            f"{path}: per-cell factor checks do not match the declared total")


def validate_common(
    path: pathlib.Path,
    payload: dict[str, Any],
    source_hash: str,
    expected_mode: str,
    expected_sanitizer_mode: str,
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
    require(payload.get("tail_policy") == TAIL_POLICY,
            f"{path}: tail policy changed")
    binding = payload.get("build_binding")
    require(isinstance(binding, dict), f"{path}: missing build binding")
    require(binding.get("source_sha256") == source_hash,
            f"{path}: source hash mismatch")
    require(binding.get("core_git_sha") == CORE_BASELINE_GIT_SHA,
            f"{path}: core git SHA is not the retained baseline")
    require(isinstance(binding.get("linked_library_sha256"), str) and
            len(binding["linked_library_sha256"]) == 64,
            f"{path}: invalid library hash")
    require(binding.get("sanitizer_mode") == expected_sanitizer_mode,
            f"{path}: wrong sanitizer mode")
    require(payload.get("direct_symbol_checks") == EXPECTED_DIRECT_CHECKS,
            f"{path}: direct check count changed")
    require(payload.get("factor_checks") == EXPECTED_FACTOR_CHECKS,
            f"{path}: factor check count changed")
    digest = payload.get("correctness_digest_fnv1a64")
    require(isinstance(digest, str) and
            re.fullmatch(r"0x[0-9a-f]{16}", digest) is not None and
            digest == EXPECTED_CORRECTNESS_DIGEST,
            f"{path}: correctness digest is not the retained oracle value")
    validate_correctness_cases(path, payload)


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


def rounded_interval(value: float) -> tuple[float, float]:
    return (
        max(0.0, value - JSON_DECIMAL_HALF_ULP),
        value + JSON_DECIMAL_HALF_ULP,
    )


def recomputed_benchmark_metrics(
    row: dict[str, Any],
) -> tuple[float, float]:
    candidate = float(row["candidate_median_us"])
    candidate_mad = float(row["candidate_mad_us"])
    padded = float(row["padded_median_us"])
    padded_mad = float(row["padded_mad_us"])
    return (
        padded / candidate,
        100.0 * ((padded - padded_mad) /
                 (candidate + candidate_mad) - 1.0),
    )


def derived_metric_intervals(
    row: dict[str, Any],
) -> tuple[tuple[float, float], tuple[float, float]]:
    candidate = rounded_interval(float(row["candidate_median_us"]))
    candidate_mad = rounded_interval(float(row["candidate_mad_us"]))
    padded = rounded_interval(float(row["padded_median_us"]))
    padded_mad = rounded_interval(float(row["padded_mad_us"]))
    require(candidate[0] > 0.0,
            "rounded candidate timing interval includes zero")

    ratios = [
        padded_value / candidate_value
        for padded_value, candidate_value in itertools.product(
            padded, candidate)
    ]
    gains = [
        100.0 * ((padded_value - padded_mad_value) /
                 (candidate_value + candidate_mad_value) - 1.0)
        for padded_value, padded_mad_value, candidate_value,
            candidate_mad_value in itertools.product(
                padded, padded_mad, candidate, candidate_mad)
    ]
    return ((min(ratios), max(ratios)), (min(gains), max(gains)))


def require_rounded_derived_value(
    path: pathlib.Path,
    key: str,
    claimed: float,
    interval: tuple[float, float],
) -> None:
    # The claimed derived value was independently rounded by the same writer.
    # Its half-ULP interval must overlap the complete interval propagated from
    # the four displayed primitives.
    require(
        claimed + JSON_DECIMAL_HALF_ULP + FLOAT_INTERVAL_SLOP >= interval[0]
        and claimed - JSON_DECIMAL_HALF_ULP - FLOAT_INTERVAL_SLOP <= interval[1],
        f"{path}: {key} does not match the displayed timing primitives",
    )


def validate_benchmark_rows(
    path: pathlib.Path,
    rows: Any,
    expected: frozenset[tuple[int, int, int, int, int]],
) -> list[dict[str, Any]]:
    require(isinstance(rows, list), f"{path}: benchmarks are absent")
    require(len(rows) == len(expected),
            f"{path}: expected {len(expected)} benchmark cells")
    seen: set[tuple[int, int, int, int, int]] = set()
    normalized: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        require(isinstance(row, dict),
                f"{path}: benchmark cell {index} is not an object")
        field_bits = require_int(row, "field_bits", path, 1)
        q = require_int(row, "q", path, 1)
        valid_bytes = require_int(row, "valid_bytes", path, 1)
        batch = require_int(row, "batch", path, 1)
        reuse = require_int(row, "reuse", path, 1)
        key = (field_bits, q, valid_bytes, batch, reuse)
        require(key not in seen, f"{path}: duplicate benchmark cell {key}")
        require(key in expected, f"{path}: unexpected benchmark cell {key}")
        seen.add(key)

        parent = require_int(row, "parent", path, 1)
        processed_bytes = require_int(row, "processed_bytes", path, 1)
        require(parent == ceil_power_of_two(q),
                f"{path}: wrong parent for benchmark cell {key}")
        require(processed_bytes == (valid_bytes + 63) // 64 * 64,
                f"{path}: wrong tail rounding for benchmark cell {key}")
        expected_samples = 7 if valid_bytes >= 65536 else 11
        require(require_int(row, "samples", path, 1) == expected_samples,
                f"{path}: wrong ABBA sample count for benchmark cell {key}")

        for numeric_key in (
            "setup_median_us", "setup_mad_us",
            "candidate_median_us", "candidate_mad_us",
            "padded_median_us", "padded_mad_us",
        ):
            require_finite_nonnegative(row, numeric_key, path)
        require(float(row["candidate_median_us"]) > 0.0 and
                float(row["padded_median_us"]) > 0.0,
                f"{path}: zero execution time in benchmark cell {key}")
        for storage_key in (
            "candidate_scratch_bytes_per_stripe",
            "padded_scratch_bytes_per_stripe", "resident_batch_bytes",
            "candidate_traffic_bytes_per_execution",
            "padded_traffic_bytes_per_execution",
        ):
            require_int(row, storage_key, path)
        accounting = expected_execution_accounting(q, processed_bytes)
        benchmark_accounting = {
            "candidate_scratch_bytes_per_stripe":
                accounting["candidate_scratch_bytes"],
            "padded_scratch_bytes_per_stripe":
                accounting["padded_scratch_bytes"],
            "resident_batch_bytes": expected_resident_batch_bytes(
                q, processed_bytes, batch),
            "candidate_traffic_bytes_per_execution":
                accounting["candidate_traffic_bytes"],
            "padded_traffic_bytes_per_execution":
                accounting["padded_traffic_bytes"],
        }
        for storage_key, expected_value in benchmark_accounting.items():
            require(row[storage_key] == expected_value,
                    f"{path}: wrong {storage_key} for benchmark cell {key}")

        claimed_ratio = require_finite_nonnegative(
            row, "padded_over_candidate", path)
        claimed_gain_value = row.get("credible_gain_percent")
        require(isinstance(claimed_gain_value, (int, float)) and
                not isinstance(claimed_gain_value, bool) and
                math.isfinite(float(claimed_gain_value)),
                f"{path}: invalid credible_gain_percent")
        ratio_interval, gain_interval = derived_metric_intervals(row)
        require_rounded_derived_value(
            path, "padded_over_candidate", claimed_ratio, ratio_interval)
        require_rounded_derived_value(
            path, "credible_gain_percent", float(claimed_gain_value),
            gain_interval)

        ratio, gain = recomputed_benchmark_metrics(row)
        canonical = dict(row)
        canonical["padded_over_candidate"] = ratio
        canonical["credible_gain_percent"] = gain
        normalized.append(canonical)

    require(seen == expected, f"{path}: benchmark matrix is incomplete")
    return normalized


def validate_benchmarks(
    path: pathlib.Path,
    payload: dict[str, Any],
) -> list[dict[str, Any]]:
    rows = payload.get("benchmarks")
    return validate_benchmark_rows(path, rows, expected_benchmark_cells())


def validate_backend_benchmarks(
    path: pathlib.Path,
    payload: dict[str, Any],
) -> list[dict[str, Any]]:
    return validate_benchmark_rows(
        path, payload.get("benchmarks"), expected_benchmark_cells(True))


def validate_forced_backend_runtime(
    path: pathlib.Path,
    label: str,
    payload: dict[str, Any],
) -> None:
    expected_runtime = "scalar" if label == "auto" else label
    require(payload.get("runtime_backend") == expected_runtime,
            f"{path}: requested {label} but runtime backend is not "
            f"{expected_runtime}")


def validate_authoritative_runtime(
    path: pathlib.Path,
    payload: dict[str, Any],
) -> None:
    require(payload.get("requested_backend") == "avx2",
            f"{path}: authoritative benchmark did not request AVX2")
    require(payload.get("runtime_backend") == "avx2",
            f"{path}: authoritative benchmark did not run AVX2")
    environment = payload.get("runtime_environment")
    require(isinstance(environment, dict),
            f"{path}: runtime environment is absent")
    affinity = environment.get("process_affinity_cpus")
    require(isinstance(affinity, list) and len(affinity) == 1 and
            isinstance(affinity[0], int) and not isinstance(affinity[0], bool)
            and affinity[0] >= 0,
            f"{path}: authoritative timing must have one-CPU affinity")
    require(environment.get("omp_num_threads_env") == "1",
            f"{path}: OMP_NUM_THREADS must be 1")
    require(environment.get("openmp_max_threads") == 1,
            f"{path}: OpenMP maximum thread count must be 1")


def summarize_benchmarks(rows: list[dict[str, Any]]) -> dict[str, Any]:
    per_q: list[dict[str, Any]] = []
    qualified_q: list[int] = []
    for q in BENCHMARK_Q:
        selected = [row for row in rows if row["q"] == q]
        require(selected, f"cannot summarize absent q={q} benchmark cells")
        gains = [recomputed_benchmark_metrics(row)[1] for row in selected]
        credible_wins = [
            row for row in selected
            if recomputed_benchmark_metrics(row)[1] >= 10.0
        ]
        byte_sizes = {row["valid_bytes"] for row in credible_wins}
        meaningful = len(credible_wins) >= 3 and len(byte_sizes) >= 2 and all(
            gain >= -2.0 for gain in gains)
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


def validate_core_revision(repository: pathlib.Path) -> None:
    try:
        completed = subprocess.run(
            ["git", "-C", str(repository), "cat-file", "-t",
             CORE_BASELINE_GIT_SHA],
            check=True, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    except (OSError, subprocess.CalledProcessError) as error:
        raise CheckpointError(
            "retained core baseline is not present in the repository") from error
    require(completed.stdout.strip() == "commit",
            "retained core baseline does not name a commit")


def run(args: argparse.Namespace) -> dict[str, Any]:
    source = pathlib.Path(args.source)
    require(source.exists(), f"missing source: {source}")
    source_hash = sha256(source)
    repository = pathlib.Path(args.repository)
    validate_core_revision(repository)
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
        validate_common(path, payload, source_hash, "backend",
                        NORMAL_SANITIZER_MODE)
        validate_library_binding(label, path, payload, backend_libraries[label])
        validate_forced_backend_runtime(path, label, payload)
        validate_backend_benchmarks(path, payload)
        loaded_backends[label] = payload
        projections.append(correctness_projection(payload))
        raw_hashes[f"backend_{label}"] = sha256(path)
        core_shas.add(payload["build_binding"]["core_git_sha"])
    require(all(item == projections[0] for item in projections[1:]),
            "backend correctness projections differ")

    sanitizer_path = pathlib.Path(args.sanitizer)
    sanitizer_library = pathlib.Path(args.sanitizer_library)
    sanitizer = load_json(sanitizer_path)
    validate_common(sanitizer_path, sanitizer, source_hash, "correctness",
                    SANITIZED_MODE)
    require(sanitizer["requested_backend"] == "asan-ubsan",
            "sanitizer label mismatch")
    require(sanitizer["build_binding"]["sanitizer_mode"] == SANITIZED_MODE,
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
    validate_common(benchmark_path, benchmark, source_hash, "all",
                    NORMAL_SANITIZER_MODE)
    validate_authoritative_runtime(benchmark_path, benchmark)
    require(benchmark["build_binding"]["linked_library_sha256"] ==
            sha256(backend_libraries["avx2"]),
            "benchmark AVX2 library hash mismatch")
    require(correctness_projection(benchmark) == projections[0],
            "benchmark correctness projection differs")
    benchmark_rows = validate_benchmarks(benchmark_path, benchmark)
    raw_hashes["benchmark"] = sha256(benchmark_path)
    core_shas.add(benchmark["build_binding"]["core_git_sha"])
    require(len(core_shas) == 1, "raw artifacts use different core revisions")

    return {
        "schema_version": "leopard2-c5-checkpoint-v1",
        "status": "pass",
        "wire_identity": "existing padded dyadic parent",
        "exact_profile_implemented": False,
        "default_build_integration": False,
        "build_binding": {
            "experiment_source_sha256": source_hash,
            "checkpoint_merger_sha256": sha256(
                pathlib.Path(__file__).resolve()),
            "core_baseline_git_sha": next(iter(core_shas)),
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
            "cells": len(benchmark_rows),
            "setup_reported_separately": True,
            "execution_order": "alternating ABBA",
            "backend": benchmark["runtime_backend"],
            "thread_count": 1,
            "affinity_cpus":
                benchmark["runtime_environment"]["process_affinity_cpus"],
            "shard_bytes": sorted({row["valid_bytes"]
                                    for row in benchmark_rows}),
            "batch": sorted({row["batch"] for row in benchmark_rows}),
            "reuse": sorted({row["reuse"] for row in benchmark_rows}),
            "setup_median_us_range": [
                min(row["setup_median_us"] for row in benchmark_rows),
                max(row["setup_median_us"] for row in benchmark_rows),
            ],
        },
        "promotion": summarize_benchmarks(benchmark_rows),
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
