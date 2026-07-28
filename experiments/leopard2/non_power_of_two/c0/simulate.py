#!/usr/bin/env python3
"""Leopard2 C0 symbolic cost and dependency simulator.

This program ranks *experiments*.  It is deliberately not a cycle model and
must not be used as a dispatch table without measured validation.  The parent
preserving radix-2 counts are exact for the structural dependency DAG defined
below.  Counts for exact-length algebra are estimates and carry an explicit
confidence label in every output row.

All arithmetic and memory counts are per shard payload byte.  A ``load`` or
``store`` means one byte moved between an application/scratch shard and the
kernel.  A fixed multiplication is multiplication of one field symbol by a
known coefficient.  ``temporary_shard_slots`` is independent of shard bytes.

The compact command evaluates all 1 <= K,R <= 256 using independent worker
processes and writes deterministic summaries and heat maps.  The matrix
command additionally writes the complete long-form matrix as deterministic
gzip-compressed CSV.  No timestamps or host-specific paths enter the output.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, fields, replace
from functools import lru_cache
from pathlib import Path
from typing import Mapping, Sequence


SCHEMA_VERSION = "leopard2-c0-v2"
PROFILES = ("legacy_high_v1", "low_v1")
METHODS = (
    "padded_dyadic_parent",
    "prefix_pruned_parent",
    "dependency_pruned_parent",
    "recursive_truncated_parent",
    "binary_dyadic_exact",
    "full_block_tail_exact",
    "direct_matrix_exact",
    "subproduct_tree_exact",
)
PARENT_METHODS = frozenset(METHODS[:4])
FOCUS_COUNTS = (3, 5, 7, 9, 15, 17, 31, 33, 63, 65, 127, 129, 255)
GF16_PAIRS = (
    (257, 1), (1, 257), (257, 255), (255, 257),
    (513, 511), (511, 513), (1000, 200), (200, 1000),
    (1025, 1023), (1023, 1025), (4096, 512), (512, 4096),
    (4097, 4095), (4095, 4097), (16384, 2048), (2048, 16384),
    (30000, 22500), (22500, 30000), (32768, 4096), (4096, 32768),
    (32768, 32768),
)
GAIN_THRESHOLDS = (1.0, 1.1, 1.25, 1.5, 2.0, 4.0)
TOP_COUNT = 24


@dataclass(frozen=True)
class Geometry:
    k: int
    r: int
    profile: str
    transform_side: int
    parent_size: int
    transmitted_size: int
    padded_coordinates: int
    shortened_coordinates: int
    punctured_coordinates: int
    parent_field: str
    exact_field: str
    padding_forces_gf16: bool


@dataclass(frozen=True)
class Metrics:
    butterfly_equivalents: int = 0
    xor_byte_ops: int = 0
    fixed_multiplications: int = 0
    load_byte_ops: int = 0
    store_byte_ops: int = 0
    temporary_shard_slots: int = 0
    irregular_boundary_ops: int = 0
    transform_passes: int = 0
    copied_shard_bytes: int = 0
    zero_filled_shard_bytes: int = 0

    def plus(self, other: "Metrics") -> "Metrics":
        values = {}
        for item in fields(Metrics):
            name = item.name
            if name == "temporary_shard_slots":
                values[name] = max(getattr(self, name), getattr(other, name))
            else:
                values[name] = getattr(self, name) + getattr(other, name)
        return Metrics(**values)

    def repeated(self, count: int) -> "Metrics":
        if count < 0:
            raise ValueError("repeat count cannot be negative")
        values = {}
        for item in fields(Metrics):
            name = item.name
            value = getattr(self, name)
            values[name] = value if name == "temporary_shard_slots" else value * count
        return Metrics(**values)

    @property
    def score(self) -> int:
        """Architecture-neutral prioritization score, not predicted cycles.

        Butterflies are reported independently and are not added because their
        component XOR, multiplication, load and store work is already counted.
        Irregular operations get a modest four-unit scheduling penalty.
        """
        return (
            self.xor_byte_ops
            + 2 * self.fixed_multiplications
            + self.load_byte_ops
            + self.store_byte_ops
            + 4 * self.irregular_boundary_ops
        )


@dataclass(frozen=True)
class CostRow:
    schema_version: str
    operation: str
    metric_unit: str
    k: int
    r: int
    profile: str
    method: str
    wire_scope: str
    confidence: str
    transform_side: int
    candidate_transform_side: int
    parent_size: int
    candidate_size: int
    transmitted_size: int
    padded_coordinates: int
    shortened_coordinates: int
    punctured_coordinates: int
    parent_inflation: float
    parent_field: str
    candidate_field: str
    padding_forces_gf16: bool
    gf8_boundary_rescue_candidate: bool
    butterfly_equivalents: int
    xor_byte_ops: int
    fixed_multiplications: int
    load_byte_ops: int
    store_byte_ops: int
    temporary_shard_slots: int
    irregular_boundary_ops: int
    transform_passes: int
    copied_shard_bytes: int
    zero_filled_shard_bytes: int
    symbolic_score: int
    gain_vs_padded: float


def ceil_pow2(value: int) -> int:
    if value <= 0:
        raise ValueError("ceil_pow2 requires a positive integer")
    return 1 << (value - 1).bit_length()


def ceil_log2(value: int) -> int:
    if value <= 0:
        raise ValueError("ceil_log2 requires a positive integer")
    return (value - 1).bit_length()


def field_for_size(value: int) -> str:
    if value <= 256:
        return "gf8"
    if value <= 65536:
        return "gf16"
    return "unsupported"


def make_geometry(k: int, r: int, profile: str) -> Geometry:
    if k <= 0 or r <= 0:
        raise ValueError("K and R must be positive")
    if profile == "legacy_high_v1":
        side = ceil_pow2(r)
        parent = ceil_pow2(k + side)
        dimension = parent - side
        shortened = dimension - k
        punctured = side - r
    elif profile == "low_v1":
        side = ceil_pow2(k)
        parent = ceil_pow2(side + r)
        shortened = side - k
        punctured = parent - side - r
    else:
        raise ValueError(f"unknown profile: {profile}")
    transmitted = k + r
    exact_field = field_for_size(transmitted)
    parent_field = field_for_size(parent)
    return Geometry(
        k=k,
        r=r,
        profile=profile,
        transform_side=side,
        parent_size=parent,
        transmitted_size=transmitted,
        padded_coordinates=parent - transmitted,
        shortened_coordinates=shortened,
        punctured_coordinates=punctured,
        parent_field=parent_field,
        exact_field=exact_field,
        padding_forces_gf16=(exact_field == "gf8" and parent_field == "gf16"),
    )


def chunks(total: int, width: int) -> tuple[int, ...]:
    if total < 0 or width <= 0:
        raise ValueError("invalid chunk dimensions")
    full, tail = divmod(total, width)
    values = [width] * full
    if tail:
        values.append(tail)
    return tuple(values)


def full_transform(length: int) -> Metrics:
    if length & (length - 1):
        raise ValueError("full transform length must be a power of two")
    layers = ceil_log2(length)
    butterflies = length * layers // 2
    return Metrics(
        butterfly_equivalents=butterflies,
        xor_byte_ops=2 * butterflies,
        fixed_multiplications=butterflies,
        load_byte_ops=2 * butterflies,
        store_byte_ops=2 * butterflies,
        temporary_shard_slots=length,
        transform_passes=layers,
    )


@lru_cache(maxsize=None)
def dependency_counts(
    length: int, active_input_prefix: int, requested_output_prefix: int
) -> tuple[int, int, int]:
    """Return input-pruned, bidirectionally-pruned, active-stage counts.

    The simulated network is the regular in-place radix-2 butterfly DAG.  A
    generic butterfly makes both outputs depend on both inputs.  Known factors
    zero and one are intentionally not guessed here; C1 must derive those from
    the actual skew schedule.
    """
    if length & (length - 1):
        raise ValueError("dependency network length must be a power of two")
    if not 0 <= active_input_prefix <= length:
        raise ValueError("active input prefix outside network")
    if not 0 <= requested_output_prefix <= length:
        raise ValueError("requested output prefix outside network")
    if length == 1 or active_input_prefix == 0 or requested_output_prefix == 0:
        return (0, 0, 0)

    live = [index < active_input_prefix for index in range(length)]
    live_before: list[list[bool]] = []
    widths: list[int] = []
    width = 1
    while width < length:
        before = live[:]
        after = live[:]
        for base in range(0, length, 2 * width):
            for offset in range(width):
                left = base + offset
                right = left + width
                value = before[left] or before[right]
                after[left] = value
                after[right] = value
        live_before.append(before)
        widths.append(width)
        live = after
        width *= 2

    needed = [index < requested_output_prefix for index in range(length)]
    needed_after: list[list[bool]] = [[False] * length for _ in widths]
    for stage in range(len(widths) - 1, -1, -1):
        needed_after[stage] = needed[:]
        before = needed[:]
        width = widths[stage]
        for base in range(0, length, 2 * width):
            for offset in range(width):
                left = base + offset
                right = left + width
                value = needed[left] or needed[right]
                before[left] = value
                before[right] = value
        needed = before

    input_pruned = 0
    bidirectional = 0
    active_stages = 0
    for stage, width in enumerate(widths):
        stage_count = 0
        for base in range(0, length, 2 * width):
            for offset in range(width):
                left = base + offset
                right = left + width
                input_live = live_before[stage][left] or live_before[stage][right]
                output_needed = (
                    needed_after[stage][left] or needed_after[stage][right]
                )
                input_pruned += int(input_live)
                stage_count += int(input_live and output_needed)
        bidirectional += stage_count
        active_stages += int(stage_count != 0)
    return input_pruned, bidirectional, active_stages


def metrics_from_butterflies(
    butterflies: int,
    passes: int,
    temp_slots: int,
    irregular_ops: int = 0,
) -> Metrics:
    return Metrics(
        butterfly_equivalents=butterflies,
        xor_byte_ops=2 * butterflies,
        fixed_multiplications=butterflies,
        load_byte_ops=2 * butterflies,
        store_byte_ops=2 * butterflies,
        temporary_shard_slots=temp_slots,
        irregular_boundary_ops=irregular_ops,
        transform_passes=passes,
    )


def transform_shapes(g: Geometry, exact: bool = False) -> tuple[tuple[int, int], ...]:
    side = (g.r if g.profile == "legacy_high_v1" else g.k) if exact else g.transform_side
    if g.profile == "legacy_high_v1":
        return tuple((count, side) for count in chunks(g.k, side)) + ((side, g.r),)
    return ((g.k, side),) + tuple((side, count) for count in chunks(g.r, side))


def common_transform_io(g: Geometry, shapes: Sequence[tuple[int, int]], padded: bool) -> Metrics:
    side = g.transform_side
    copied = sum(active for active, _ in shapes)
    zeroed = sum(side - active for active, _ in shapes) if padded else 0
    # Requested parity is scattered once into application storage.  Transform
    # stores are already accounted for above.
    result = Metrics(
        load_byte_ops=copied + g.r,
        store_byte_ops=copied + zeroed + g.r,
        copied_shard_bytes=copied + g.r,
        zero_filled_shard_bytes=zeroed,
    )
    if g.profile == "legacy_high_v1":
        block_count = len(chunks(g.k, side))
        accumulated = side * max(0, block_count - 1)
        result = result.plus(Metrics(
            xor_byte_ops=accumulated,
            load_byte_ops=2 * accumulated,
            store_byte_ops=accumulated,
        ))
    return result


def parent_metrics(g: Geometry, method: str) -> Metrics:
    side = g.transform_side
    shapes = transform_shapes(g)
    result = Metrics()
    padded = method == "padded_dyadic_parent"
    for active, requested in shapes:
        if padded:
            item = full_transform(side)
        else:
            prefix, bidirectional, passes = dependency_counts(side, active, requested)
            if method == "prefix_pruned_parent":
                count = prefix
                # Prefix-only execution still walks every live stage.
                passes = ceil_log2(side) if count else 0
                irregular = int(active != side)
                temp = 2 * side
            elif method == "dependency_pruned_parent":
                count = bidirectional
                irregular = int(active != side) + int(requested != side)
                temp = 2 * side
            elif method == "recursive_truncated_parent":
                count = bidirectional
                ragged_paths = (
                    int(active not in (0, side)) + int(requested not in (0, side))
                )
                irregular = ragged_paths * ceil_log2(side)
                # A lower bound for a packed recursive frontier.  C2 must
                # measure the real schedule and scratch requirement.
                temp = max(1, active, requested)
            else:
                raise ValueError(f"not a parent method: {method}")
            item = metrics_from_butterflies(count, passes, temp, irregular)
        result = result.plus(item)
    # The regular implementation reuses two side-wide workspaces.
    if method != "recursive_truncated_parent":
        result = replace(result, temporary_shard_slots=2 * side)
    return result.plus(common_transform_io(g, shapes, padded))


def binary_parts(value: int) -> tuple[int, ...]:
    if value <= 0:
        raise ValueError("binary decomposition requires a positive value")
    return tuple(
        1 << bit
        for bit in range(value.bit_length() - 1, -1, -1)
        if value & (1 << bit)
    )


def binary_dyadic_one(length: int) -> Metrics:
    parts = binary_parts(length)
    result = Metrics()
    for part in parts:
        result = result.plus(full_transform(part))
    joins = max(0, len(parts) - 1)
    # Optimistic O(q) work per aligned-coset join.  If C5 derives a dense
    # cross-block map, this lower bound must be replaced by that measured map.
    cross = length * joins
    return result.plus(Metrics(
        xor_byte_ops=cross,
        fixed_multiplications=cross,
        load_byte_ops=2 * cross,
        store_byte_ops=cross,
        temporary_shard_slots=length,
        irregular_boundary_ops=joins,
        transform_passes=joins,
    ))


def dense_tail(length: int) -> Metrics:
    head = 1 << (length.bit_length() - 1)
    tail = length - head
    if tail == 0:
        return full_transform(length)
    head_metrics = full_transform(head)
    # Tail interpolation plus cross-block factors.  This intentionally charges
    # the potentially dense head/tail join instead of assuming it is free.
    tail_products = tail * tail
    cross_products = head * tail
    products = tail_products + cross_products
    xors = max(0, tail - 1) * tail + cross_products
    tail_metrics = Metrics(
        xor_byte_ops=xors,
        fixed_multiplications=products,
        load_byte_ops=products + xors,
        store_byte_ops=tail + cross_products,
        temporary_shard_slots=tail,
        irregular_boundary_ops=1,
        transform_passes=1,
    )
    return head_metrics.plus(tail_metrics).plus(Metrics(
        temporary_shard_slots=length,
        irregular_boundary_ops=1,
        transform_passes=1,
    ))


def exact_transform_metrics(g: Geometry, method: str) -> Metrics:
    side = g.r if g.profile == "legacy_high_v1" else g.k
    shapes = transform_shapes(g, exact=True)
    if method == "binary_dyadic_exact":
        one = binary_dyadic_one(side)
    elif method == "full_block_tail_exact":
        one = dense_tail(side)
    else:
        raise ValueError(f"unknown exact transform method: {method}")
    result = one.repeated(len(shapes))
    copied = sum(active for active, _ in shapes) + g.r
    result = result.plus(Metrics(
        load_byte_ops=copied,
        store_byte_ops=copied,
        copied_shard_bytes=copied,
    ))
    if g.profile == "legacy_high_v1":
        block_count = len(chunks(g.k, side))
        accumulated = side * max(0, block_count - 1)
        result = result.plus(Metrics(
            xor_byte_ops=accumulated,
            load_byte_ops=2 * accumulated,
            store_byte_ops=accumulated,
        ))
    return replace(result, temporary_shard_slots=side)


def direct_matrix_metrics(g: Geometry) -> Metrics:
    products = g.k * g.r
    xors = max(0, g.k - 1) * g.r
    return Metrics(
        xor_byte_ops=xors,
        fixed_multiplications=products,
        load_byte_ops=products,
        store_byte_ops=g.r,
        temporary_shard_slots=1,
        transform_passes=1,
    )


def subproduct_tree_metrics(g: Geometry) -> Metrics:
    # Three generic polynomial phases (tree/interpolate/evaluate), using the
    # conservative M(q) log(q) proxy q*ceil(log2(q))^2 per phase.
    total = g.transmitted_size
    layers = ceil_log2(total)
    units = 3 * total * layers * layers
    return Metrics(
        butterfly_equivalents=units,
        xor_byte_ops=units,
        fixed_multiplications=units,
        load_byte_ops=2 * units,
        store_byte_ops=units,
        temporary_shard_slots=4 * total,
        irregular_boundary_ops=layers,
        transform_passes=3 * layers,
    )


def method_metrics(g: Geometry) -> dict[str, Metrics]:
    return {
        "padded_dyadic_parent": parent_metrics(g, "padded_dyadic_parent"),
        "prefix_pruned_parent": parent_metrics(g, "prefix_pruned_parent"),
        "dependency_pruned_parent": parent_metrics(g, "dependency_pruned_parent"),
        "recursive_truncated_parent": parent_metrics(g, "recursive_truncated_parent"),
        "binary_dyadic_exact": exact_transform_metrics(g, "binary_dyadic_exact"),
        "full_block_tail_exact": exact_transform_metrics(g, "full_block_tail_exact"),
        "direct_matrix_exact": direct_matrix_metrics(g),
        "subproduct_tree_exact": subproduct_tree_metrics(g),
    }


def make_row(g: Geometry, method: str, metrics: Metrics, padded_score: int) -> CostRow:
    parent = method in PARENT_METHODS
    candidate_size = g.parent_size if parent else g.transmitted_size
    candidate_field = g.parent_field if parent else g.exact_field
    if method in ("padded_dyadic_parent", "prefix_pruned_parent", "dependency_pruned_parent"):
        confidence = "exact_structural_radix2_dag"
    elif method == "recursive_truncated_parent":
        confidence = "exact_arithmetic_lower_bound_scratch"
    elif method == "direct_matrix_exact":
        confidence = "exact_dense_generator_arithmetic"
    elif method == "binary_dyadic_exact":
        confidence = "heuristic_optimistic_cross_block"
    elif method == "full_block_tail_exact":
        confidence = "heuristic_conservative_dense_join"
    else:
        confidence = "asymptotic_generic_polynomial_proxy"
    score = metrics.score
    gain = math.inf if score == 0 and padded_score else (
        1.0 if score == 0 else padded_score / score
    )
    return CostRow(
        schema_version=SCHEMA_VERSION,
        operation="systematic_encode",
        metric_unit="per_payload_byte",
        k=g.k,
        r=g.r,
        profile=g.profile,
        method=method,
        wire_scope="parent_preserving" if parent else "new_exact_profile",
        confidence=confidence,
        transform_side=g.transform_side,
        candidate_transform_side=(
            g.transform_side
            if parent
            else (g.r if g.profile == "legacy_high_v1" else g.k)
        ),
        parent_size=g.parent_size,
        candidate_size=candidate_size,
        transmitted_size=g.transmitted_size,
        padded_coordinates=g.padded_coordinates,
        shortened_coordinates=g.shortened_coordinates,
        punctured_coordinates=g.punctured_coordinates,
        parent_inflation=round(g.parent_size / g.transmitted_size, 9),
        parent_field=g.parent_field,
        candidate_field=candidate_field,
        padding_forces_gf16=g.padding_forces_gf16,
        gf8_boundary_rescue_candidate=(g.padding_forces_gf16 and candidate_field == "gf8"),
        butterfly_equivalents=metrics.butterfly_equivalents,
        xor_byte_ops=metrics.xor_byte_ops,
        fixed_multiplications=metrics.fixed_multiplications,
        load_byte_ops=metrics.load_byte_ops,
        store_byte_ops=metrics.store_byte_ops,
        temporary_shard_slots=metrics.temporary_shard_slots,
        irregular_boundary_ops=metrics.irregular_boundary_ops,
        transform_passes=metrics.transform_passes,
        copied_shard_bytes=metrics.copied_shard_bytes,
        zero_filled_shard_bytes=metrics.zero_filled_shard_bytes,
        symbolic_score=score,
        gain_vs_padded=round(gain, 9),
    )


def rows_for(k: int, r: int) -> tuple[CostRow, ...]:
    rows: list[CostRow] = []
    for profile in PROFILES:
        g = make_geometry(k, r, profile)
        costs = method_metrics(g)
        padded_score = costs["padded_dyadic_parent"].score
        rows.extend(make_row(g, method, costs[method], padded_score) for method in METHODS)
    return tuple(rows)


def row_dict(row: CostRow) -> dict[str, object]:
    return asdict(row)


def flat_row(row: CostRow) -> dict[str, object]:
    result = asdict(row)
    result["padding_forces_gf16"] = str(row.padding_forces_gf16).lower()
    result["gf8_boundary_rescue_candidate"] = str(
        row.gf8_boundary_rescue_candidate
    ).lower()
    result["parent_inflation"] = f"{row.parent_inflation:.9f}"
    result["gain_vs_padded"] = f"{row.gain_vs_padded:.9f}"
    return result


def _top_insert(items: list[tuple[float, int, int]], value: tuple[float, int, int]) -> None:
    items.append(value)
    items.sort(key=lambda item: (-item[0], item[1], item[2]))
    del items[TOP_COUNT:]


def analyze_k(args: tuple[int, int, bool]) -> dict[str, object]:
    k, r_max, include_heatmaps = args
    winners: Counter[str] = Counter()
    thresholds: Counter[str] = Counter()
    score_sums: Counter[str] = Counter()
    gain_sums: defaultdict[str, float] = defaultdict(float)
    method_cells: Counter[str] = Counter()
    top: defaultdict[str, list[tuple[float, int, int]]] = defaultdict(list)
    focus: list[dict[str, object]] = []
    heatmaps: defaultdict[str, list[float]] = defaultdict(list)
    rescue_profile_cells = 0
    inflation_sum = 0.0
    inflation_max = 0.0
    profile_cells = 0
    for r in range(1, r_max + 1):
        cell_rows = rows_for(k, r)
        for profile in PROFILES:
            selected = [row for row in cell_rows if row.profile == profile]
            winner = min(selected, key=lambda row: (row.symbolic_score, METHODS.index(row.method)))
            winners[f"{profile}|{winner.method}"] += 1
            profile_cells += 1
            inflation_sum += selected[0].parent_inflation
            inflation_max = max(inflation_max, selected[0].parent_inflation)
            if selected[0].padding_forces_gf16:
                rescue_profile_cells += 1
            for row in selected:
                key = f"{profile}|{row.method}"
                method_cells[key] += 1
                score_sums[key] += row.symbolic_score
                gain_sums[key] += row.gain_vs_padded
                for threshold in GAIN_THRESHOLDS:
                    if row.gain_vs_padded >= threshold:
                        thresholds[f"{key}|{threshold:g}"] += 1
                _top_insert(top[key], (row.gain_vs_padded, k, r))
                if include_heatmaps:
                    heatmaps[key].append(row.gain_vs_padded)
                if k in FOCUS_COUNTS and r in FOCUS_COUNTS:
                    focus.append(row_dict(row))
    return {
        "k": k,
        "winners": dict(winners),
        "thresholds": dict(thresholds),
        "score_sums": dict(score_sums),
        "gain_sums": dict(gain_sums),
        "method_cells": dict(method_cells),
        "top": {key: values for key, values in top.items()},
        "focus": focus,
        "heatmaps": dict(heatmaps),
        "rescue_profile_cells": rescue_profile_cells,
        "inflation_sum": inflation_sum,
        "inflation_max": inflation_max,
        "profile_cells": profile_cells,
    }


def merge_analysis(batches: Sequence[Mapping[str, object]]) -> dict[str, object]:
    winners: Counter[str] = Counter()
    thresholds: Counter[str] = Counter()
    score_sums: Counter[str] = Counter()
    method_cells: Counter[str] = Counter()
    gain_sums: defaultdict[str, float] = defaultdict(float)
    top: defaultdict[str, list[tuple[float, int, int]]] = defaultdict(list)
    focus: list[dict[str, object]] = []
    rescue = 0
    inflation_sum = 0.0
    inflation_max = 0.0
    profile_cells = 0
    for batch in batches:
        winners.update(batch["winners"])  # type: ignore[arg-type]
        thresholds.update(batch["thresholds"])  # type: ignore[arg-type]
        score_sums.update(batch["score_sums"])  # type: ignore[arg-type]
        method_cells.update(batch["method_cells"])  # type: ignore[arg-type]
        for key, value in batch["gain_sums"].items():  # type: ignore[union-attr]
            gain_sums[key] += float(value)
        for key, values in batch["top"].items():  # type: ignore[union-attr]
            for value in values:
                _top_insert(top[key], tuple(value))  # type: ignore[arg-type]
        focus.extend(batch["focus"])  # type: ignore[arg-type]
        rescue += int(batch["rescue_profile_cells"])
        inflation_sum += float(batch["inflation_sum"])
        inflation_max = max(inflation_max, float(batch["inflation_max"]))
        profile_cells += int(batch["profile_cells"])

    method_summary = []
    for profile in PROFILES:
        for method in METHODS:
            key = f"{profile}|{method}"
            cells = method_cells[key]
            method_summary.append({
                "profile": profile,
                "method": method,
                "cells": cells,
                "average_symbolic_score": round(score_sums[key] / cells, 6),
                "average_gain_vs_padded": round(gain_sums[key] / cells, 9),
                "gain_threshold_cells": {
                    f"{threshold:g}": thresholds[f"{key}|{threshold:g}"]
                    for threshold in GAIN_THRESHOLDS
                },
                "top_gain_cells": [
                    {"gain": gain, "k": k, "r": r}
                    for gain, k, r in top[key]
                ],
            })
    winner_summary = []
    for profile in PROFILES:
        for method in METHODS:
            count = winners[f"{profile}|{method}"]
            if count:
                winner_summary.append({
                    "profile": profile,
                    "method": method,
                    "cells": count,
                })
    focus.sort(key=lambda row: (
        int(row["k"]), int(row["r"]), PROFILES.index(str(row["profile"])),
        METHODS.index(str(row["method"])),
    ))
    return {
        "schema_version": SCHEMA_VERSION,
        "operation": "systematic_encode",
        "metric_unit": "per_payload_byte",
        "input_pair_count": profile_cells // len(PROFILES),
        "profile_cell_count": profile_cells,
        "method_row_count": profile_cells * len(METHODS),
        "padding_forces_gf16_profile_cell_count": rescue,
        "average_parent_inflation": round(inflation_sum / profile_cells, 9),
        "maximum_parent_inflation": round(inflation_max, 9),
        "winner_counts": winner_summary,
        "method_summary": method_summary,
        "focus_counts": list(FOCUS_COUNTS),
        "focus_rows": focus,
    }


def analyze_pairs(pairs: Sequence[tuple[int, int]]) -> dict[str, object]:
    rows = [row for k, r in sorted(set(pairs)) for row in rows_for(k, r)]
    winners: Counter[str] = Counter()
    for k, r in sorted(set(pairs)):
        selected_cell = [row for row in rows if row.k == k and row.r == r]
        for profile in PROFILES:
            selected = [row for row in selected_cell if row.profile == profile]
            winner = min(selected, key=lambda row: (row.symbolic_score, METHODS.index(row.method)))
            winners[f"{profile}|{winner.method}"] += 1
    return {
        "pairs": [{"k": k, "r": r} for k, r in sorted(set(pairs))],
        "profile_cell_count": len(set(pairs)) * len(PROFILES),
        "method_row_count": len(rows),
        "winner_counts": [
            {"profile": profile, "method": method, "cells": winners[f"{profile}|{method}"]}
            for profile in PROFILES for method in METHODS
            if winners[f"{profile}|{method}"]
        ],
        "rows": [row_dict(row) for row in rows],
    }


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            stream.write(text)
        os.replace(temp_name, path)
    except BaseException:
        Path(temp_name).unlink(missing_ok=True)
        raise


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(
    path: Path,
    rows: Sequence[Mapping[str, object]],
    field_names: Sequence[str] | None = None,
) -> None:
    if field_names is None:
        if not rows:
            raise ValueError("empty CSV requires explicit field names")
        field_names = list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=field_names, lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temp_name, path)
    except BaseException:
        Path(temp_name).unlink(missing_ok=True)
        raise


def write_focus_csv(path: Path, focus_rows: Sequence[Mapping[str, object]]) -> None:
    rows = []
    for row in focus_rows:
        converted = dict(row)
        converted["padding_forces_gf16"] = str(converted["padding_forces_gf16"]).lower()
        converted["gf8_boundary_rescue_candidate"] = str(
            converted["gf8_boundary_rescue_candidate"]
        ).lower()
        rows.append(converted)
    write_csv(path, rows, [item.name for item in fields(CostRow)])


def write_heatmaps(
    output_dir: Path, batches: Sequence[Mapping[str, object]], r_max: int
) -> list[Path]:
    heatmap_dir = output_dir / "heatmaps"
    heatmap_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for profile in PROFILES:
        for method in METHODS:
            key = f"{profile}|{method}"
            path = heatmap_dir / f"{profile}__{method}.csv"
            descriptor, temp_name = tempfile.mkstemp(
                prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, text=True
            )
            try:
                with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
                    writer = csv.writer(stream, lineterminator="\n")
                    writer.writerow(["k\\r"] + list(range(1, r_max + 1)))
                    for batch in batches:
                        writer.writerow([batch["k"]] + batch["heatmaps"][key])  # type: ignore[index]
                os.replace(temp_name, path)
            except BaseException:
                Path(temp_name).unlink(missing_ok=True)
                raise
            paths.append(path)
    return paths


def full_matrix_part(args: tuple[int, int, str]) -> str:
    k, r_max, directory = args
    path = Path(directory) / f"part-{k:04d}.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[item.name for item in fields(CostRow)],
            lineterminator="\n",
        )
        for r in range(1, r_max + 1):
            for row in rows_for(k, r):
                writer.writerow(flat_row(row))
    return str(path)


def write_full_matrix(path: Path, k_max: int, r_max: int, workers: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="leopard2-c0-matrix-", dir=path.parent) as directory:
        inputs = [(k, r_max, directory) for k in range(1, k_max + 1)]
        with ProcessPoolExecutor(max_workers=workers) as executor:
            parts = list(executor.map(full_matrix_part, inputs, chunksize=1))
        descriptor, temp_name = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
        )
        os.close(descriptor)
        try:
            with Path(temp_name).open("wb") as raw:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                    header = ",".join(item.name for item in fields(CostRow)) + "\n"
                    compressed.write(header.encode("utf-8"))
                    for part_name in parts:
                        with Path(part_name).open("rb") as part:
                            shutil.copyfileobj(part, compressed)
            os.replace(temp_name, path)
        except BaseException:
            Path(temp_name).unlink(missing_ok=True)
            raise


def output_manifest(output_dir: Path, paths: Sequence[Path], workers: int) -> dict[str, object]:
    files = []
    for path in sorted(paths):
        files.append({
            "path": str(path.relative_to(output_dir)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    return {
        "schema_version": SCHEMA_VERSION,
        "workers": workers,
        "deterministic_output": True,
        "files": files,
    }


def command_compact(args: argparse.Namespace) -> int:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    inputs = [(k, args.r_max, args.heatmaps) for k in range(1, args.k_max + 1)]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        batches = list(executor.map(analyze_k, inputs, chunksize=1))
    summary = merge_analysis(batches)
    focus_rows = summary.pop("focus_rows")
    gf16 = analyze_pairs(GF16_PAIRS)
    summary_path = output_dir / "summary.json"
    focus_path = output_dir / "focus.csv"
    gf16_path = output_dir / "gf16_samples.json"
    atomic_write_text(summary_path, json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_focus_csv(focus_path, focus_rows)  # type: ignore[arg-type]
    atomic_write_text(gf16_path, json.dumps(gf16, indent=2, sort_keys=True) + "\n")
    paths = [summary_path, focus_path, gf16_path]
    if args.heatmaps:
        paths.extend(write_heatmaps(output_dir, batches, args.r_max))
    if args.full_matrix:
        matrix_path = output_dir / "cost_matrix.csv.gz"
        write_full_matrix(matrix_path, args.k_max, args.r_max, args.workers)
        paths.append(matrix_path)
    manifest = output_manifest(output_dir, paths, args.workers)
    manifest_path = output_dir / "manifest.json"
    atomic_write_text(manifest_path, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "schema_version": SCHEMA_VERSION,
        "input_pair_count": summary["input_pair_count"],
        "profile_cell_count": summary["profile_cell_count"],
        "method_row_count": summary["method_row_count"],
        "padding_forces_gf16_profile_cell_count": summary[
            "padding_forces_gf16_profile_cell_count"
        ],
        "workers": args.workers,
        "manifest": str(manifest_path),
    }, indent=2, sort_keys=True))
    return 0


def brute_dependency_counts(length: int, active: int, requested: int) -> tuple[int, int]:
    """Independent graph trace used only by self-test for small networks."""
    nodes: list[list[tuple[int, int]]] = []
    live = [index < active for index in range(length)]
    live_states = []
    width = 1
    while width < length:
        pairs = []
        next_live = live[:]
        for base in range(0, length, 2 * width):
            for offset in range(width):
                pair = (base + offset, base + offset + width)
                pairs.append(pair)
                value = live[pair[0]] or live[pair[1]]
                next_live[pair[0]] = value
                next_live[pair[1]] = value
        nodes.append(pairs)
        live_states.append(live)
        live = next_live
        width *= 2
    prefix = sum(
        int(state[left] or state[right])
        for state, pairs in zip(live_states, nodes)
        for left, right in pairs
    )
    required = {index for index in range(requested)}
    selected: list[set[int]] = [set() for _ in nodes]
    for stage in range(len(nodes) - 1, -1, -1):
        previous = set(required)
        for index, (left, right) in enumerate(nodes[stage]):
            if left in required or right in required:
                selected[stage].add(index)
                previous.add(left)
                previous.add(right)
        required = previous
    bidirectional = sum(
        int(index in selected[stage] and (live_states[stage][left] or live_states[stage][right]))
        for stage, pairs in enumerate(nodes)
        for index, (left, right) in enumerate(pairs)
    )
    return prefix, bidirectional


def self_test() -> None:
    check_count = 0

    def check(condition: bool, label: str) -> None:
        nonlocal check_count
        check_count += 1
        if not condition:
            raise RuntimeError("C0 simulator self-test failed: " + label)

    def assert_line_numbers(source: str) -> tuple[int, ...]:
        return tuple(
            node.lineno for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Assert)
        )

    check(ceil_pow2(1) == 1, "ceil_pow2 one")
    check(ceil_pow2(255) == 256, "ceil_pow2 rounds upward")
    check(ceil_pow2(257) == 512, "ceil_pow2 crosses GF8 boundary")
    check(chunks(10, 4) == (4, 4, 2), "chunk decomposition")
    check(binary_parts(13) == (8, 4, 1), "binary decomposition")

    high = make_geometry(240, 16, "legacy_high_v1")
    check(
        (high.transform_side, high.parent_size) == (16, 256),
        "legacy-high geometry",
    )
    check(
        (high.shortened_coordinates, high.punctured_coordinates) == (0, 0),
        "legacy-high exact public geometry",
    )
    low = make_geometry(16, 240, "low_v1")
    check(
        (low.transform_side, low.parent_size) == (16, 256),
        "low geometry",
    )
    check(
        (low.shortened_coordinates, low.punctured_coordinates) == (0, 0),
        "low exact public geometry",
    )
    rescue = make_geometry(129, 100, "legacy_high_v1")
    check(
        rescue.parent_size == 512 and rescue.padding_forces_gf16,
        "GF8 field-boundary rescue geometry",
    )

    # Exhaust every structural prefix pair through length 32 against a second
    # graph representation, and check monotonic/bounded invariants to 256.
    checked_networks = 0
    for length in (1, 2, 4, 8, 16, 32):
        full = full_transform(length).butterfly_equivalents
        for active in range(1, length + 1):
            previous = -1
            for requested in range(1, length + 1):
                prefix, dependency, passes = dependency_counts(length, active, requested)
                brute_prefix, brute_dependency = brute_dependency_counts(
                    length, active, requested
                )
                label = "length={} active={} requested={}".format(
                    length, active, requested)
                check(
                    (prefix, dependency) ==
                    (brute_prefix, brute_dependency),
                    "independent dependency trace " + label,
                )
                check(
                    0 <= dependency <= prefix <= full,
                    "dependency bounds " + label,
                )
                check(
                    dependency >= previous,
                    "dependency monotonicity " + label,
                )
                check(
                    0 <= passes <= ceil_log2(length),
                    "dependency pass bounds " + label,
                )
                previous = dependency
                checked_networks += 1
    for length in (64, 128, 256):
        full = full_transform(length).butterfly_equivalents
        for active in (1, 3, length // 2, length - 1, length):
            previous = -1
            for requested in (1, 3, length // 2, length - 1, length):
                prefix, dependency, _ = dependency_counts(length, active, requested)
                label = "length={} active={} requested={}".format(
                    length, active, requested)
                check(
                    0 <= dependency <= prefix <= full,
                    "large dependency bounds " + label,
                )
                check(
                    dependency >= previous,
                    "large dependency monotonicity " + label,
                )
                previous = dependency
                checked_networks += 1

    checked_rows = 0
    for k in (1,) + FOCUS_COUNTS + (256, 257, 1000, 4096):
        for r in (1,) + FOCUS_COUNTS + (256, 257, 1000):
            metrics_by_profile = {
                profile: method_metrics(make_geometry(k, r, profile))
                for profile in PROFILES
            }
            for profile in PROFILES:
                profile_metrics = metrics_by_profile[profile]
                padded = profile_metrics["padded_dyadic_parent"]
                prefix = profile_metrics["prefix_pruned_parent"]
                dependency = profile_metrics["dependency_pruned_parent"]
                recursive = profile_metrics["recursive_truncated_parent"]
                label = f"profile={profile} K={k} R={r}"
                check(
                    dependency.butterfly_equivalents <=
                    prefix.butterfly_equivalents,
                    "dependency versus prefix butterflies " + label,
                )
                check(
                    prefix.butterfly_equivalents <=
                    padded.butterfly_equivalents,
                    "prefix versus padded butterflies " + label,
                )
                check(
                    recursive.butterfly_equivalents ==
                    dependency.butterfly_equivalents,
                    "recursive dependency butterflies " + label,
                )
                check(
                    recursive.temporary_shard_slots <=
                    2 * make_geometry(k, r, profile).transform_side,
                    "recursive scratch bound " + label,
                )
            for row in rows_for(k, r):
                label = (
                    f"profile={row.profile} method={row.method} K={k} R={r}")
                check(
                    row.k == k and row.r == r,
                    "row public geometry " + label,
                )
                check(
                    row.method in METHODS and row.profile in PROFILES,
                    "row method/profile identity " + label,
                )
                check(row.symbolic_score >= 0, "row score " + label)
                check(
                    row.parent_size >= row.transmitted_size,
                    "row parent bound " + label,
                )
                check(
                    row.padded_coordinates ==
                    row.parent_size - row.transmitted_size,
                    "row padding identity " + label,
                )
                for name in fields(Metrics):
                    value = getattr(
                        metrics_by_profile[row.profile][row.method], name.name)
                    check(
                        value >= 0,
                        "row metric {} {}".format(name.name, label),
                    )
                if row.method == "padded_dyadic_parent" and row.symbolic_score:
                    check(
                        row.gain_vs_padded == 1.0,
                        "padded baseline gain " + label,
                    )
                if row.wire_scope == "parent_preserving":
                    check(
                        row.candidate_size == row.parent_size,
                        "parent-preserving size " + label,
                    )
                checked_rows += 1

    exact_side_row = next(
        row
        for row in rows_for(129, 100)
        if row.profile == "legacy_high_v1" and row.method == "direct_matrix_exact"
    )
    check(exact_side_row.transform_side == 128, "exact rescue padded side")
    check(
        exact_side_row.candidate_transform_side == 100,
        "exact rescue candidate side",
    )
    check(exact_side_row.parent_field == "gf16", "exact rescue parent field")
    check(
        exact_side_row.candidate_field == "gf8",
        "exact rescue candidate field",
    )
    check(
        exact_side_row.gf8_boundary_rescue_candidate,
        "exact rescue classification",
    )

    # Dense systematic generator work is monotonic in K and R.
    for k in range(1, 32):
        for r in range(1, 32):
            here = direct_matrix_metrics(make_geometry(k, r, "low_v1"))
            check(
                direct_matrix_metrics(
                    make_geometry(k + 1, r, "low_v1")).
                fixed_multiplications >= here.fixed_multiplications,
                f"direct K monotonicity K={k} R={r}",
            )
            check(
                direct_matrix_metrics(
                    make_geometry(k, r + 1, "low_v1")).
                fixed_multiplications >= here.fixed_multiplications,
                f"direct R monotonicity K={k} R={r}",
            )

    # Stable serialization and deterministic parallel worker payload.
    first = [row_dict(row) for row in rows_for(7, 9)]
    second = [row_dict(row) for row in rows_for(7, 9)]
    check(first == second, "stable row serialization")
    batch_first = analyze_k((3, 9, True))
    batch_second = analyze_k((3, 9, True))
    check(batch_first == batch_second, "stable worker payload")
    merged = merge_analysis([analyze_k((k, 3, False)) for k in range(1, 4)])
    check(merged["input_pair_count"] == 9, "merged input-pair count")
    check(merged["profile_cell_count"] == 18, "merged profile-cell count")
    check(merged["method_row_count"] == 144, "merged method-row count")
    check(
        sum(item["cells"] for item in merged["winner_counts"]
            if item["profile"] == "low_v1") == 9,
        "merged low-profile winner count",
    )
    check(checked_networks == 1440, "dependency network coverage count")
    check(checked_rows == 4896, "cost-row coverage count")
    check(
        assert_line_numbers("value = 1\n") == (),
        "assert scanner accepts assert-free source",
    )
    check(
        assert_line_numbers("value = 1\nassert value\n") == (2,),
        "assert scanner detects optimized-away checks",
    )
    check(
        assert_line_numbers(Path(__file__).read_text(encoding="utf-8")) == (),
        "C0 simulator contains no Python assert statements",
    )
    check(check_count > 1, "self-test executed a nontrivial check set")
    print(json.dumps({
        "status": "PASS",
        "checked_dependency_networks": checked_networks,
        "checked_cost_rows": checked_rows,
        "checks": check_count,
        "schema_version": SCHEMA_VERSION,
    }, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("self-test", help="run deterministic algebra/model invariants")
    compact = subparsers.add_parser("compact", help="run sweep and write deterministic evidence")
    compact.add_argument("--output-dir", required=True)
    compact.add_argument("--k-max", type=int, default=256)
    compact.add_argument("--r-max", type=int, default=256)
    compact.add_argument("--workers", type=int, default=min(128, os.cpu_count() or 1))
    compact.add_argument("--heatmaps", action="store_true")
    compact.add_argument("--full-matrix", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "self-test":
        self_test()
        return 0
    if args.k_max <= 0 or args.r_max <= 0:
        parser.error("--k-max and --r-max must be positive")
    if args.workers <= 0:
        parser.error("--workers must be positive")
    return command_compact(args)


if __name__ == "__main__":
    raise SystemExit(main())
