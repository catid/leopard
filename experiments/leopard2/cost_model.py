#!/usr/bin/env python3
"""Symbolic cost model for Leopard2 non-power-of-two experiments.

This is an implementation-prioritization model, not a performance predictor.
Every count is expressed per payload byte so that it is independent of shard
size.  The model deliberately keeps two classes of candidates separate:

* parent-preserving methods compute the existing shortened/punctured dyadic
  code and therefore may be wire-compatible;
* exact-profile methods use exactly K+R evaluation coordinates and require a
  separately versioned profile unless generator equivalence is later proved.

The radix-2 network counts for the padded and pruning variants are symbolic.
The irregular exact-transform counts are explicitly marked as estimates.  The
integer ranking score is only a stable way to prioritize experiments:

    4 * butterflies + 2 * fixed multiplications + XOR byte units
      + load byte units + store byte units + 4 * irregular boundaries

It intentionally does not encode an assumed CPU's instruction throughput.

The generated long-form CSV/JSONL schema is versioned.  Output order and JSON
key order are deterministic, timestamps are omitted, and files are replaced
atomically.  A full 1<=K,R<=256 sweep streams its rows and uses bounded memory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import os
import tempfile
from dataclasses import asdict, dataclass, fields, replace
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Iterator, Sequence, TextIO


SCHEMA_VERSION = "leopard2-c0-cost-v1"
METHOD_ORDER = (
    "padded_dyadic",
    "prefix_pruned_parent",
    "dependency_pruned_parent",
    "recursive_truncated_parent",
    "direct_matrix_exact",
    "recursive_truncated_exact",
    "binary_dyadic_exact",
    "full_block_tail_exact",
    "generic_subproduct_exact",
)
PROFILE_ORDER = ("legacy_high_v1", "low_v1")
EXACT_GF8_LIMIT = 256
GF16_LIMIT = 65536


@dataclass(frozen=True)
class Geometry:
    k: int
    r: int
    profile: str
    parent_size: int
    transform_side: int
    stream_side: int
    transmitted_length: int
    padded_coordinates: int
    baseline_field: str
    field_boundary_rescue_possible: bool


@dataclass(frozen=True)
class Metrics:
    butterfly_equivalents: int = 0
    xor_byte_units: int = 0
    fixed_multiplications: int = 0
    load_byte_units: int = 0
    store_byte_units: int = 0
    temporary_shard_slots: int = 0
    irregular_boundaries: int = 0
    transform_passes: int = 0
    zero_fill_slots: int = 0
    copy_slots: int = 0

    def plus(self, other: "Metrics") -> "Metrics":
        return Metrics(**{
            field.name: getattr(self, field.name) + getattr(other, field.name)
            for field in fields(Metrics)
        })

    def scaled(self, factor: int) -> "Metrics":
        return Metrics(**{
            field.name: getattr(self, field.name) * factor
            for field in fields(Metrics)
        })

    @property
    def weighted_cost(self) -> int:
        return (
            4 * self.butterfly_equivalents
            + 2 * self.fixed_multiplications
            + self.xor_byte_units
            + self.load_byte_units
            + self.store_byte_units
            + 4 * self.irregular_boundaries
        )


@dataclass(frozen=True)
class CostRow:
    schema_version: str
    operation_scope: str
    metric_unit: str
    k: int
    r: int
    profile: str
    method: str
    estimate_kind: str
    wire_compatible_candidate: bool
    baseline_field: str
    candidate_field: str
    profile_parent_size: int
    candidate_parent_size: int
    transform_side: int
    transmitted_length: int
    padded_coordinates: int
    parent_inflation_numerator: int
    parent_inflation_denominator: int
    field_boundary_rescue_possible: bool
    field_boundary_rescue_candidate: bool
    butterfly_equivalents: int
    xor_byte_units: int
    fixed_multiplications: int
    load_byte_units: int
    store_byte_units: int
    temporary_shard_slots: int
    irregular_boundaries: int
    transform_passes: int
    zero_fill_slots: int
    copy_slots: int
    weighted_cost: int
    gain_vs_padded: float


def ceil_pow2(value: int) -> int:
    if value <= 0:
        raise ValueError("value must be positive")
    return 1 << (value - 1).bit_length()


def ceil_log2(value: int) -> int:
    if value <= 0:
        raise ValueError("value must be positive")
    return (value - 1).bit_length()


def field_for_length(length: int) -> str:
    if length <= EXACT_GF8_LIMIT:
        return "gf8"
    if length <= GF16_LIMIT:
        return "gf16"
    return "unsupported"


def geometry(k: int, r: int, profile: str) -> Geometry:
    if k <= 0 or r <= 0:
        raise ValueError("K and R must be positive")
    if profile == "legacy_high_v1":
        side = ceil_pow2(r)
        parent = ceil_pow2(k + side)
        stream = k
    elif profile == "low_v1":
        side = ceil_pow2(k)
        parent = ceil_pow2(side + r)
        stream = r
    else:
        raise ValueError(f"unknown profile: {profile}")
    transmitted = k + r
    baseline_field = field_for_length(parent)
    return Geometry(
        k=k,
        r=r,
        profile=profile,
        parent_size=parent,
        transform_side=side,
        stream_side=stream,
        transmitted_length=transmitted,
        padded_coordinates=parent - transmitted,
        baseline_field=baseline_field,
        field_boundary_rescue_possible=(
            transmitted <= EXACT_GF8_LIMIT and parent > EXACT_GF8_LIMIT
        ),
    )


def chunks(total: int, width: int) -> tuple[int, ...]:
    full, tail = divmod(total, width)
    result = [width] * full
    if tail:
        result.append(tail)
    return tuple(result)


def full_transform_metrics(length: int) -> Metrics:
    """One complete radix-2 transform, per payload byte."""
    layers = ceil_log2(length)
    butterflies = length * layers // 2
    return Metrics(
        butterfly_equivalents=butterflies,
        xor_byte_units=2 * butterflies,
        fixed_multiplications=butterflies,
        load_byte_units=2 * butterflies,
        store_byte_units=2 * butterflies,
        temporary_shard_slots=length,
        transform_passes=layers,
    )


@lru_cache(maxsize=None)
def network_counts(
    length: int, active_input_count: int, requested_output_count: int
) -> tuple[int, int, int]:
    """Return (prefix-only, bidirectional, active stages) butterfly counts.

    A radix-2 butterfly is structurally live when either input is live.  For
    bidirectional dependency pruning it must additionally feed a requested
    output.  Fixed factors equal to zero or one are deliberately not assumed;
    those savings belong to later coefficient-aware schedule simulation.
    """
    if length & (length - 1):
        raise ValueError("network length must be a power of two")
    if not 0 <= active_input_count <= length:
        raise ValueError("active input count outside transform")
    if not 0 <= requested_output_count <= length:
        raise ValueError("requested output count outside transform")
    if length == 1 or active_input_count == 0 or requested_output_count == 0:
        return (0, 0, 0)

    active = [index < active_input_count for index in range(length)]
    active_before: list[list[bool]] = []
    active_after: list[list[bool]] = []
    widths: list[int] = []
    width = 1
    while width < length:
        before = active[:]
        after = active[:]
        span = width * 2
        for base in range(0, length, span):
            for offset in range(width):
                left = base + offset
                right = left + width
                live = before[left] or before[right]
                after[left] = live
                after[right] = live
        active_before.append(before)
        active_after.append(after)
        widths.append(width)
        active = after
        width = span

    needed = [index < requested_output_count for index in range(length)]
    needed_after: list[list[bool]] = [[False] * length for _ in widths]
    for stage in range(len(widths) - 1, -1, -1):
        needed_after[stage] = needed[:]
        before = needed[:]
        width = widths[stage]
        span = width * 2
        for base in range(0, length, span):
            for offset in range(width):
                left = base + offset
                right = left + width
                live = needed[left] or needed[right]
                before[left] = live
                before[right] = live
        needed = before

    prefix_count = 0
    dependency_count = 0
    active_stages = 0
    for stage, width in enumerate(widths):
        span = width * 2
        stage_dependency = 0
        for base in range(0, length, span):
            for offset in range(width):
                left = base + offset
                right = left + width
                input_live = (
                    active_before[stage][left] or active_before[stage][right]
                )
                output_needed = (
                    needed_after[stage][left] or needed_after[stage][right]
                )
                prefix_count += int(input_live)
                stage_dependency += int(input_live and output_needed)
        dependency_count += stage_dependency
        active_stages += int(stage_dependency != 0)
    return prefix_count, dependency_count, active_stages


def transform_metrics_from_count(
    butterflies: int,
    passes: int,
    temporary_slots: int,
    irregular_boundaries: int = 0,
) -> Metrics:
    return Metrics(
        butterfly_equivalents=butterflies,
        xor_byte_units=2 * butterflies,
        fixed_multiplications=butterflies,
        load_byte_units=2 * butterflies,
        store_byte_units=2 * butterflies,
        temporary_shard_slots=temporary_slots,
        irregular_boundaries=irregular_boundaries,
        transform_passes=passes,
    )


def combine_transforms(parts: Iterable[Metrics]) -> Metrics:
    result = Metrics()
    max_temp = 0
    for part in parts:
        max_temp = max(max_temp, part.temporary_shard_slots)
        result = result.plus(replace(part, temporary_shard_slots=0))
    return replace(result, temporary_shard_slots=max_temp)


def transform_shapes(g: Geometry) -> tuple[tuple[int, int], ...]:
    """(active inputs, requested outputs) for each padded-side transform."""
    side = g.transform_side
    if g.profile == "legacy_high_v1":
        # Message-block IFFTs construct all side coefficients.  The final FFT
        # has all accumulated coefficients but only R transmitted outputs.
        return tuple((count, side) for count in chunks(g.k, side)) + ((side, g.r),)
    # One systematic IFFT constructs all coefficients.  Shifted parity FFTs
    # request one transmitted block at a time.
    return ((g.k, side),) + tuple((side, count) for count in chunks(g.r, side))


def preparation_metrics(g: Geometry, method: str) -> Metrics:
    shapes = transform_shapes(g)
    copied = sum(active for active, _ in shapes)
    if method == "padded_dyadic":
        zeroed = sum(g.transform_side - active for active, _ in shapes)
    else:
        zeroed = 0
    # Copies and zero fills are memory operations, not transform butterflies.
    return Metrics(
        load_byte_units=copied,
        store_byte_units=copied + zeroed,
        zero_fill_slots=zeroed,
        copy_slots=copied,
    )


def padded_parent_metrics(g: Geometry) -> Metrics:
    shapes = transform_shapes(g)
    transforms = full_transform_metrics(g.transform_side).scaled(len(shapes))
    transforms = replace(
        transforms,
        # Scratch is reused across transforms rather than multiplied.
        temporary_shard_slots=2 * g.transform_side,
    )
    return transforms.plus(preparation_metrics(g, "padded_dyadic"))


def pruned_parent_metrics(g: Geometry, dependency: bool) -> Metrics:
    parts = []
    boundaries = 0
    for active, requested in transform_shapes(g):
        prefix, bidirectional, passes = network_counts(
            g.transform_side, active, requested
        )
        butterflies = bidirectional if dependency else prefix
        boundary = int(active != g.transform_side or requested != g.transform_side)
        boundaries += boundary
        parts.append(transform_metrics_from_count(
            butterflies, passes, g.transform_side, boundary
        ))
    combined = combine_transforms(parts).plus(preparation_metrics(g, "pruned"))
    return replace(
        combined,
        temporary_shard_slots=2 * g.transform_side,
        irregular_boundaries=boundaries,
    )


def recursive_parent_metrics(g: Geometry) -> Metrics:
    """Wire-compatible recursive execution of the bidirectional live DAG.

    Arithmetic counts equal the dependency DAG.  The different experiment is
    represented by eliminating suffix materialization and by charging one
    irregular recursive boundary per ragged input/output path.
    """
    parts = []
    boundaries = 0
    maximum_materialized = 1
    for active, requested in transform_shapes(g):
        _, butterflies, passes = network_counts(g.transform_side, active, requested)
        input_boundary = ceil_log2(g.transform_side) if active != g.transform_side else 0
        output_boundary = (
            ceil_log2(g.transform_side) if requested != g.transform_side else 0
        )
        boundary = input_boundary + output_boundary
        boundaries += boundary
        maximum_materialized = max(maximum_materialized, active, requested)
        parts.append(transform_metrics_from_count(
            butterflies, passes, max(active, requested), boundary
        ))
    combined = combine_transforms(parts).plus(preparation_metrics(g, "recursive"))
    return replace(
        combined,
        temporary_shard_slots=maximum_materialized,
        irregular_boundaries=boundaries,
    )


def exact_transform_repetitions(g: Geometry) -> int:
    exact_side = g.r if g.profile == "legacy_high_v1" else g.k
    stream = g.k if g.profile == "legacy_high_v1" else g.r
    return 1 + (stream + exact_side - 1) // exact_side


def exact_side(g: Geometry) -> int:
    return g.r if g.profile == "legacy_high_v1" else g.k


def recursive_exact_one(length: int) -> Metrics:
    layers = ceil_log2(length)
    butterflies = (length * layers + 1) // 2
    boundary = layers if length != ceil_pow2(length) else 0
    return transform_metrics_from_count(
        butterflies, layers, length, boundary
    )


def recursive_exact_metrics(g: Geometry) -> Metrics:
    side = exact_side(g)
    repeats = exact_transform_repetitions(g)
    result = recursive_exact_one(side).scaled(repeats)
    return replace(
        result,
        temporary_shard_slots=side,
        copy_slots=g.k + g.r,
        load_byte_units=result.load_byte_units + g.k + g.r,
        store_byte_units=result.store_byte_units + g.k + g.r,
    )


def binary_parts(value: int) -> tuple[int, ...]:
    return tuple(
        1 << bit
        for bit in range(value.bit_length() - 1, -1, -1)
        if value & (1 << bit)
    )


def binary_exact_one(length: int) -> Metrics:
    parts = binary_parts(length)
    transforms = combine_transforms(full_transform_metrics(part) for part in parts)
    # Joining s aligned additive cosets requires a candidate derivation.  Until
    # then charge one linear cross-block pass per join, rather than pretending
    # the operation is free or assuming a dense matrix.
    joins = max(0, len(parts) - 1)
    join_units = length * joins
    return Metrics(
        butterfly_equivalents=transforms.butterfly_equivalents,
        xor_byte_units=transforms.xor_byte_units + join_units,
        fixed_multiplications=transforms.fixed_multiplications + join_units,
        load_byte_units=transforms.load_byte_units + 2 * join_units,
        store_byte_units=transforms.store_byte_units + join_units,
        temporary_shard_slots=length,
        irregular_boundaries=joins,
        transform_passes=transforms.transform_passes + joins,
    )


def binary_exact_metrics(g: Geometry) -> Metrics:
    side = exact_side(g)
    repeats = exact_transform_repetitions(g)
    result = binary_exact_one(side).scaled(repeats)
    return replace(
        result,
        temporary_shard_slots=side,
        copy_slots=g.k + g.r,
        load_byte_units=result.load_byte_units + g.k + g.r,
        store_byte_units=result.store_byte_units + g.k + g.r,
    )


def tail_exact_one(length: int) -> Metrics:
    if length == ceil_pow2(length):
        return full_transform_metrics(length)
    largest = 1 << (length.bit_length() - 1)
    tail = length - largest
    head = full_transform_metrics(largest)
    tail_parent = ceil_pow2(tail)
    padded_tail = full_transform_metrics(tail_parent)
    dense_tail_cost = tail * length
    padded_tail_score = padded_tail.weighted_cost
    if dense_tail_cost * 6 <= padded_tail_score:
        tail_part = Metrics(
            xor_byte_units=max(0, tail - 1) * tail,
            fixed_multiplications=dense_tail_cost,
            load_byte_units=dense_tail_cost,
            store_byte_units=tail,
            temporary_shard_slots=tail,
            irregular_boundaries=1,
            transform_passes=1,
        )
    else:
        tail_part = padded_tail
    # Cross-block factors touch each tail element once per head layer.
    cross = tail * ceil_log2(largest)
    result = combine_transforms((head, tail_part))
    return replace(
        result,
        xor_byte_units=result.xor_byte_units + cross,
        fixed_multiplications=result.fixed_multiplications + cross,
        load_byte_units=result.load_byte_units + 2 * cross,
        store_byte_units=result.store_byte_units + cross,
        temporary_shard_slots=length,
        irregular_boundaries=result.irregular_boundaries + 1,
        transform_passes=result.transform_passes + 1,
    )


def tail_exact_metrics(g: Geometry) -> Metrics:
    side = exact_side(g)
    repeats = exact_transform_repetitions(g)
    result = tail_exact_one(side).scaled(repeats)
    return replace(
        result,
        temporary_shard_slots=side,
        copy_slots=g.k + g.r,
        load_byte_units=result.load_byte_units + g.k + g.r,
        store_byte_units=result.store_byte_units + g.k + g.r,
    )


def direct_matrix_metrics(g: Geometry) -> Metrics:
    products = g.k * g.r
    xors = max(0, g.k - 1) * g.r
    return Metrics(
        xor_byte_units=xors,
        fixed_multiplications=products,
        load_byte_units=products,
        store_byte_units=g.r,
        temporary_shard_slots=1,
        transform_passes=1,
        copy_slots=g.k,
    )


def subproduct_metrics(g: Geometry) -> Metrics:
    """Conservative M(n) log(n) estimate for generic exact interpolation.

    This is not an LCH transform count.  Three polynomial phases (tree,
    interpolation, evaluation) are each charged n*ceil(log2(n))^2 units.
    """
    total = g.transmitted_length
    layers = ceil_log2(total)
    units = 3 * total * layers * layers
    return Metrics(
        butterfly_equivalents=units,
        xor_byte_units=units,
        fixed_multiplications=units,
        load_byte_units=2 * units,
        store_byte_units=units,
        temporary_shard_slots=4 * total,
        irregular_boundaries=layers,
        transform_passes=3 * layers,
        copy_slots=g.k + g.r,
    )


def all_method_metrics(g: Geometry) -> dict[str, Metrics]:
    return {
        "padded_dyadic": padded_parent_metrics(g),
        "prefix_pruned_parent": pruned_parent_metrics(g, dependency=False),
        "dependency_pruned_parent": pruned_parent_metrics(g, dependency=True),
        "recursive_truncated_parent": recursive_parent_metrics(g),
        "direct_matrix_exact": direct_matrix_metrics(g),
        "recursive_truncated_exact": recursive_exact_metrics(g),
        "binary_dyadic_exact": binary_exact_metrics(g),
        "full_block_tail_exact": tail_exact_metrics(g),
        "generic_subproduct_exact": subproduct_metrics(g),
    }


def row_for(g: Geometry, method: str, metrics: Metrics, padded_cost: int) -> CostRow:
    parent_candidate = method.endswith("_parent") or method == "padded_dyadic"
    candidate_length = g.parent_size if parent_candidate else g.transmitted_length
    candidate_field = g.baseline_field if parent_candidate else field_for_length(candidate_length)
    if method in ("padded_dyadic", "prefix_pruned_parent", "dependency_pruned_parent"):
        estimate_kind = "symbolic_radix2"
    elif method == "direct_matrix_exact":
        estimate_kind = "symbolic_dense"
    else:
        estimate_kind = "heuristic_estimate"
    gain = 0.0 if metrics.weighted_cost == 0 else padded_cost / metrics.weighted_cost
    return CostRow(
        schema_version=SCHEMA_VERSION,
        operation_scope="systematic_encoding",
        metric_unit="per_payload_byte",
        k=g.k,
        r=g.r,
        profile=g.profile,
        method=method,
        estimate_kind=estimate_kind,
        wire_compatible_candidate=parent_candidate,
        baseline_field=g.baseline_field,
        candidate_field=candidate_field,
        profile_parent_size=g.parent_size,
        candidate_parent_size=candidate_length,
        transform_side=g.transform_side,
        transmitted_length=g.transmitted_length,
        padded_coordinates=g.padded_coordinates,
        parent_inflation_numerator=g.parent_size,
        parent_inflation_denominator=g.transmitted_length,
        field_boundary_rescue_possible=g.field_boundary_rescue_possible,
        field_boundary_rescue_candidate=(
            g.field_boundary_rescue_possible and candidate_field == "gf8"
        ),
        butterfly_equivalents=metrics.butterfly_equivalents,
        xor_byte_units=metrics.xor_byte_units,
        fixed_multiplications=metrics.fixed_multiplications,
        load_byte_units=metrics.load_byte_units,
        store_byte_units=metrics.store_byte_units,
        temporary_shard_slots=metrics.temporary_shard_slots,
        irregular_boundaries=metrics.irregular_boundaries,
        transform_passes=metrics.transform_passes,
        zero_fill_slots=metrics.zero_fill_slots,
        copy_slots=metrics.copy_slots,
        weighted_cost=metrics.weighted_cost,
        gain_vs_padded=round(gain, 6),
    )


def rows_for_pair(k: int, r: int, profiles: Sequence[str]) -> Iterator[CostRow]:
    for profile in profiles:
        g = geometry(k, r, profile)
        metrics = all_method_metrics(g)
        padded_cost = metrics["padded_dyadic"].weighted_cost
        for method in METHOD_ORDER:
            yield row_for(g, method, metrics[method], padded_cost)


def parse_range(spec: str) -> tuple[int, ...]:
    """Parse inclusive START:STOP[:STEP] or a comma-separated value list."""
    if ":" not in spec:
        values = tuple(int(part) for part in spec.split(",") if part)
    else:
        parts = tuple(int(part) for part in spec.split(":"))
        if len(parts) not in (2, 3):
            raise ValueError("range must be START:STOP[:STEP]")
        start, stop = parts[:2]
        step = parts[2] if len(parts) == 3 else 1
        if step <= 0:
            raise ValueError("range step must be positive")
        values = tuple(range(start, stop + 1, step))
    if not values or min(values) <= 0:
        raise ValueError("ranges must contain positive integers")
    return tuple(sorted(set(values)))


def iter_pairs(
    k_values: Sequence[int],
    r_values: Sequence[int],
    gf16_k_values: Sequence[int],
    gf16_r_values: Sequence[int],
) -> Iterator[tuple[int, int]]:
    pairs = {(k, r) for k in k_values for r in r_values}
    if gf16_k_values or gf16_r_values:
        if not gf16_k_values or not gf16_r_values:
            raise ValueError("both GF16 K and R ranges must be provided")
        pairs.update((k, r) for k in gf16_k_values for r in gf16_r_values)
    yield from sorted(pairs)


class AtomicTextFile:
    def __init__(self, destination: Path):
        self.destination = destination
        self.temp_path: Path | None = None
        self.stream: TextIO | None = None

    def __enter__(self) -> TextIO:
        self.destination.parent.mkdir(parents=True, exist_ok=True)
        descriptor, name = tempfile.mkstemp(
            dir=self.destination.parent,
            prefix=f".{self.destination.name}.",
            suffix=".tmp",
            text=True,
        )
        self.temp_path = Path(name)
        self.stream = os.fdopen(descriptor, "w", encoding="utf-8", newline="")
        return self.stream

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        assert self.stream is not None
        assert self.temp_path is not None
        self.stream.close()
        if exc_type is None:
            os.replace(self.temp_path, self.destination)
        else:
            self.temp_path.unlink(missing_ok=True)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class Summary:
    def __init__(self, top_count: int = 32, threshold: float = 1.10):
        self.top_count = top_count
        self.threshold = threshold
        self.heaps: dict[tuple[str, str], list[tuple[float, int, int]]] = {}
        self.winner_counts: dict[tuple[str, str], int] = {}
        self._current_pair: tuple[int, int, str] | None = None
        self._current_rows: list[CostRow] = []
        self.intervals: list[dict[str, object]] = []
        self._open_interval: dict[tuple[str, int], dict[str, object]] = {}
        self.rows = 0
        self.pairs = 0
        self.rescue_pairs: set[tuple[int, int, str]] = set()

    def add(self, row: CostRow) -> None:
        self.rows += 1
        pair = (row.k, row.r, row.profile)
        if self._current_pair is not None and pair != self._current_pair:
            self._finish_pair()
        if self._current_pair is None or pair != self._current_pair:
            self._current_pair = pair
            self._current_rows = []
        self._current_rows.append(row)
        gain = row.gain_vs_padded
        key = (row.profile, row.method)
        heap = self.heaps.setdefault(key, [])
        item = (gain, -row.k, -row.r)
        if len(heap) < self.top_count:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)
        if row.field_boundary_rescue_candidate:
            self.rescue_pairs.add(pair)

    def _finish_pair(self) -> None:
        if not self._current_rows:
            return
        self.pairs += 1
        winner = min(
            self._current_rows,
            key=lambda row: (row.weighted_cost, METHOD_ORDER.index(row.method)),
        )
        self.winner_counts[(winner.profile, winner.method)] = (
            self.winner_counts.get((winner.profile, winner.method), 0) + 1
        )
        gain = winner.gain_vs_padded
        interval_key = (winner.profile, winner.k)
        previous = self._open_interval.get(interval_key)
        if (
            gain >= self.threshold
            and previous is not None
            and previous["method"] == winner.method
            and previous["r_end"] == winner.r - 1
        ):
            previous["r_end"] = winner.r
            previous["gain_sum"] = float(previous["gain_sum"]) + gain
            previous["max_gain"] = max(float(previous["max_gain"]), gain)
            previous["cells"] = int(previous["cells"]) + 1
        else:
            if previous is not None:
                self.intervals.append(previous)
            if gain >= self.threshold:
                self._open_interval[interval_key] = {
                    "profile": winner.profile,
                    "k": winner.k,
                    "r_start": winner.r,
                    "r_end": winner.r,
                    "method": winner.method,
                    "gain_sum": gain,
                    "max_gain": gain,
                    "cells": 1,
                }
            else:
                self._open_interval.pop(interval_key, None)
        self._current_pair = None
        self._current_rows = []

    def finish(self) -> dict[str, object]:
        self._finish_pair()
        self.intervals.extend(self._open_interval.values())
        top_cells = []
        for (profile, method), heap in sorted(self.heaps.items()):
            cells = [
                {"gain": round(gain, 6), "k": -neg_k, "r": -neg_r}
                for gain, neg_k, neg_r in sorted(heap, reverse=True)
            ]
            top_cells.append({
                "profile": profile,
                "method": method,
                "cells": cells,
            })
        regions = []
        for item in self.intervals:
            cells = int(item["cells"])
            regions.append({
                "profile": item["profile"],
                "k": item["k"],
                "r_start": item["r_start"],
                "r_end": item["r_end"],
                "method": item["method"],
                "cells": cells,
                "average_gain": round(float(item["gain_sum"]) / cells, 6),
                "max_gain": round(float(item["max_gain"]), 6),
            })
        regions.sort(key=lambda item: (
            -int(item["cells"]),
            -float(item["average_gain"]),
            str(item["profile"]),
            int(item["k"]),
            int(item["r_start"]),
        ))
        winners = [
            {"profile": profile, "method": method, "cells": count}
            for (profile, method), count in sorted(self.winner_counts.items())
        ]
        return {
            "row_count": self.rows,
            "profile_pair_count": self.pairs,
            "field_boundary_rescue_profile_pair_count": len(self.rescue_pairs),
            "winner_counts": winners,
            "top_cells": top_cells,
            "top_regions": regions[: self.top_count],
        }


def bool_for_csv(value: bool) -> str:
    return "true" if value else "false"


def csv_dict(row: CostRow) -> dict[str, object]:
    result = asdict(row)
    for key, value in tuple(result.items()):
        if isinstance(value, bool):
            result[key] = bool_for_csv(value)
    result["gain_vs_padded"] = f"{row.gain_vs_padded:.6f}"
    return result


def write_matrix(
    output_dir: Path,
    pairs: Iterable[tuple[int, int]],
    profiles: Sequence[str],
    formats: Sequence[str],
    top_count: int,
) -> dict[str, object]:
    csv_path = output_dir / "cost_matrix.csv"
    jsonl_path = output_dir / "cost_matrix.jsonl"
    summary = Summary(top_count=top_count)
    field_names = [field.name for field in fields(CostRow)]
    pair_list = tuple(pairs)

    # A result directory represents one invocation.  Do not leave an old
    # matrix in a format that the current invocation did not request.
    for format_name, path in (("csv", csv_path), ("jsonl", jsonl_path)):
        if format_name not in formats:
            path.unlink(missing_ok=True)

    csv_context = AtomicTextFile(csv_path) if "csv" in formats else None
    json_context = AtomicTextFile(jsonl_path) if "jsonl" in formats else None
    csv_stream = csv_context.__enter__() if csv_context else None
    json_stream = json_context.__enter__() if json_context else None
    try:
        writer = None
        if csv_stream:
            writer = csv.DictWriter(csv_stream, fieldnames=field_names, lineterminator="\n")
            writer.writeheader()
        for k, r in pair_list:
            for row in rows_for_pair(k, r, profiles):
                summary.add(row)
                if writer:
                    writer.writerow(csv_dict(row))
                if json_stream:
                    json.dump(asdict(row), json_stream, separators=(",", ":"))
                    json_stream.write("\n")
    except BaseException:
        if json_context:
            json_context.__exit__(*__import__("sys").exc_info())
        if csv_context:
            csv_context.__exit__(*__import__("sys").exc_info())
        raise
    else:
        if json_context:
            json_context.__exit__(None, None, None)
        if csv_context:
            csv_context.__exit__(None, None, None)

    result = summary.finish()
    files = []
    for path in (csv_path, jsonl_path):
        if path.exists() and path.suffix.lstrip(".") in formats:
            files.append({
                "name": path.name,
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            })
    result.update({
        "schema_version": SCHEMA_VERSION,
        "columns": field_names,
        "methods": list(METHOD_ORDER),
        "profiles": list(profiles),
        "input_pair_count": len(pair_list),
        "files": files,
    })
    with AtomicTextFile(output_dir / "summary.json") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


def self_test() -> None:
    assert ceil_pow2(1) == 1
    assert ceil_pow2(3) == 4
    assert ceil_pow2(256) == 256
    assert chunks(10, 4) == (4, 4, 2)
    assert binary_parts(13) == (8, 4, 1)

    high = geometry(240, 16, "legacy_high_v1")
    assert (high.transform_side, high.parent_size) == (16, 256)
    low = geometry(16, 240, "low_v1")
    assert (low.transform_side, low.parent_size) == (16, 256)
    rescue = geometry(129, 100, "legacy_high_v1")
    assert rescue.parent_size == 512
    assert rescue.transmitted_length == 229
    assert rescue.field_boundary_rescue_possible

    for length in (1, 2, 4, 8, 16, 256):
        prefix, dependency, passes = network_counts(length, length, length)
        expected = full_transform_metrics(length).butterfly_equivalents
        assert prefix == expected
        assert dependency == expected
        assert passes == ceil_log2(length)

    for length in (2, 4, 8, 16, 32):
        padded = full_transform_metrics(length).butterfly_equivalents
        for active in range(1, length + 1):
            previous = -1
            for requested in range(1, length + 1):
                prefix, dependency, _ = network_counts(length, active, requested)
                assert 0 <= dependency <= prefix <= padded
                assert dependency >= previous
                previous = dependency

    for profile in PROFILE_ORDER:
        for k, r in ((1, 1), (3, 5), (16, 16), (129, 100), (255, 1)):
            g = geometry(k, r, profile)
            metrics = all_method_metrics(g)
            assert tuple(metrics) == METHOD_ORDER
            padded = metrics["padded_dyadic"].weighted_cost
            assert padded >= 0
            for method, cost in metrics.items():
                assert cost.weighted_cost >= 0, method
                for field in fields(Metrics):
                    assert getattr(cost, field.name) >= 0
            padded_row = row_for(g, "padded_dyadic", metrics["padded_dyadic"], padded)
            if padded:
                assert padded_row.gain_vs_padded == 1.0

    # Dense generator work is monotone in either public dimension.
    for k in range(1, 16):
        for r in range(1, 16):
            here = direct_matrix_metrics(geometry(k, r, "low_v1")).fixed_multiplications
            next_k = direct_matrix_metrics(geometry(k + 1, r, "low_v1")).fixed_multiplications
            next_r = direct_matrix_metrics(geometry(k, r + 1, "low_v1")).fixed_multiplications
            assert next_k >= here
            assert next_r >= here

    parent_row = next(
        row for row in rows_for_pair(129, 100, ("legacy_high_v1",))
        if row.method == "dependency_pruned_parent"
    )
    exact_row = next(
        row for row in rows_for_pair(129, 100, ("legacy_high_v1",))
        if row.method == "direct_matrix_exact"
    )
    assert parent_row.candidate_field == "gf16"
    assert not parent_row.field_boundary_rescue_candidate
    assert exact_row.candidate_field == "gf8"
    assert exact_row.field_boundary_rescue_candidate

    # Stable construction: identical inputs yield identical rows and key order.
    first = [asdict(row) for row in rows_for_pair(7, 9, PROFILE_ORDER)]
    second = [asdict(row) for row in rows_for_pair(7, 9, PROFILE_ORDER)]
    assert first == second
    assert list(first[0]) == [field.name for field in fields(CostRow)]

    # Exercise the streaming/atomic output path, not only the pure model.
    with tempfile.TemporaryDirectory(prefix="leopard2-cost-self-test-") as directory:
        output_dir = Path(directory)
        pairs = ((3, 5), (7, 9))
        expected_rows = len(pairs) * len(PROFILE_ORDER) * len(METHOD_ORDER)
        result = write_matrix(
            output_dir, pairs, PROFILE_ORDER, ("csv", "jsonl"), top_count=4
        )
        assert result["row_count"] == expected_rows
        assert result["profile_pair_count"] == len(pairs) * len(PROFILE_ORDER)
        assert len((output_dir / "cost_matrix.csv").read_text(
            encoding="utf-8"
        ).splitlines()) == expected_rows + 1
        json_lines = (output_dir / "cost_matrix.jsonl").read_text(
            encoding="utf-8"
        ).splitlines()
        assert len(json_lines) == expected_rows
        assert all(json.loads(line)["schema_version"] == SCHEMA_VERSION
                   for line in json_lines)
        first_hashes = {
            path.name: file_sha256(path)
            for path in (
                output_dir / "cost_matrix.csv",
                output_dir / "cost_matrix.jsonl",
                output_dir / "summary.json",
            )
        }
        write_matrix(
            output_dir, pairs, PROFILE_ORDER, ("csv", "jsonl"), top_count=4
        )
        second_hashes = {
            path.name: file_sha256(path)
            for path in (
                output_dir / "cost_matrix.csv",
                output_dir / "cost_matrix.jsonl",
                output_dir / "summary.json",
            )
        }
        assert first_hashes == second_hashes
        write_matrix(output_dir, pairs, PROFILE_ORDER, ("csv",), top_count=4)
        assert not (output_dir / "cost_matrix.jsonl").exists()


def command_generate(args: argparse.Namespace) -> int:
    k_values = parse_range(args.k_range)
    r_values = parse_range(args.r_range)
    gf16_k_values = parse_range(args.gf16_k_range) if args.gf16_k_range else ()
    gf16_r_values = parse_range(args.gf16_r_range) if args.gf16_r_range else ()
    pairs = tuple(iter_pairs(k_values, r_values, gf16_k_values, gf16_r_values))
    profiles = PROFILE_ORDER if args.profile == "all" else (args.profile,)
    formats = ("csv", "jsonl") if args.format == "both" else (args.format,)
    summary = write_matrix(
        Path(args.output_dir), pairs, profiles, formats, args.top_count
    )
    print(json.dumps({
        "schema_version": summary["schema_version"],
        "input_pair_count": summary["input_pair_count"],
        "profile_pair_count": summary["profile_pair_count"],
        "row_count": summary["row_count"],
        "files": summary["files"],
    }, indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("self-test", help="run deterministic model invariants")
    generate = subparsers.add_parser("generate", help="stream a cost matrix")
    generate.add_argument(
        "--k-range", default="1:256",
        help="inclusive START:STOP[:STEP] or comma-separated values",
    )
    generate.add_argument(
        "--r-range", default="1:256",
        help="inclusive START:STOP[:STEP] or comma-separated values",
    )
    generate.add_argument(
        "--gf16-k-range", default="",
        help="optional additional GF16 K values/range",
    )
    generate.add_argument(
        "--gf16-r-range", default="",
        help="optional additional GF16 R values/range",
    )
    generate.add_argument(
        "--profile", choices=("all",) + PROFILE_ORDER, default="all",
    )
    generate.add_argument(
        "--format", choices=("csv", "jsonl", "both"), default="csv",
    )
    generate.add_argument("--output-dir", required=True)
    generate.add_argument("--top-count", type=int, default=32)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "self-test":
        self_test()
        print("cost_model self-test: PASS")
        return 0
    if args.top_count <= 0:
        parser.error("--top-count must be positive")
    return command_generate(args)


if __name__ == "__main__":
    raise SystemExit(main())
