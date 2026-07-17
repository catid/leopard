#!/usr/bin/env python3
"""Experiment F: deterministic received-coordinate subset selection.

This is a scalar correctness and symbolic-cost checkpoint.  It does not alter
the production decoder, profile identity, or dispatcher, and it deliberately
does not report wall-clock performance.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import hashlib
import itertools
import json
import math
import os
import sys
from pathlib import Path
from typing import Iterable, Sequence


SCHEMA = "leopard2-received-subset-checkpoint-v2"
FIELD_ORDER = 16
FIELD_POLYNOMIAL = 0x13  # x^4 + x + 1
POLICIES = (
    "lowest_parent_coordinate",
    "prefer_systematic",
    "block_aligned_greedy",
    "exact_block_dp",
)
FULL_COMPARISON_KEYS = (
    "strictly_better_than_prefer_systematic",
    "equal_to_prefer_systematic",
    "strictly_worse_than_prefer_systematic",
    "strictly_worse_than_exact",
    "equal_to_exact",
)
# Unprefixed comparison counters cover SelectionCost.total_key: all three
# schedule metrics followed by the selected-coordinate tuple.  Metric-prefixed
# counters deliberately omit only that final deterministic tie-break.
METRIC_COMPARISON_KEYS = tuple(
    f"metric_{key}" for key in FULL_COMPARISON_KEYS
)
TIEBREAK_COMPARISON_KEY = "metric_equal_but_tiebreak_worse_than_exact"
COMPARISON_KEYS = (
    FULL_COMPARISON_KEYS + METRIC_COMPARISON_KEYS +
    (TIEBREAK_COMPARISON_KEY,)
)


def gf_multiply_slow(left: int, right: int) -> int:
    if not 0 <= left < FIELD_ORDER or not 0 <= right < FIELD_ORDER:
        raise ValueError("GF(2^4) symbol outside the field")
    result = 0
    value = left
    multiplier = right
    while multiplier:
        if multiplier & 1:
            result ^= value
        multiplier >>= 1
        value <<= 1
        if value & FIELD_ORDER:
            value ^= FIELD_POLYNOMIAL
    return result


GF_MULTIPLY = tuple(
    tuple(gf_multiply_slow(left, right) for right in range(FIELD_ORDER))
    for left in range(FIELD_ORDER)
)


def gf_power(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = GF_MULTIPLY[result][value]
        value = GF_MULTIPLY[value][value]
        exponent >>= 1
    return result


GF_INVERSE = (0,) + tuple(gf_power(value, FIELD_ORDER - 2)
                          for value in range(1, FIELD_ORDER))


def ceil_power_of_two(value: int) -> int:
    if value <= 0:
        raise ValueError("power-of-two input must be positive")
    return 1 << (value - 1).bit_length()


@dataclasses.dataclass(frozen=True)
class Geometry:
    profile: str
    k: int
    r: int
    side: int
    parent: int
    dimension: int
    original_coordinates: tuple[int, ...]
    recovery_coordinates: tuple[int, ...]
    source_coordinates: tuple[int, ...]

    @property
    def public_coordinates(self) -> tuple[int, ...]:
        return self.original_coordinates + self.recovery_coordinates


def make_geometry(profile: str, k: int, r: int) -> Geometry:
    if k <= 0 or r <= 0:
        raise ValueError("K and R must be positive")
    if profile == "legacy_high_v1":
        side = ceil_power_of_two(r)
        parent = ceil_power_of_two(k + side)
        dimension = parent - side
        originals = tuple(range(side, side + k))
        recoveries = tuple(range(r))
        sources = tuple(range(side, parent))
    elif profile == "low_v1":
        side = ceil_power_of_two(k)
        parent = ceil_power_of_two(side + r)
        dimension = side
        originals = tuple(range(k))
        recoveries = tuple(range(side, side + r))
        sources = tuple(range(side))
    else:
        raise ValueError("unknown profile")
    if parent > FIELD_ORDER:
        raise ValueError("parent does not fit GF(2^4)")
    if len(sources) != dimension:
        raise AssertionError("source-coordinate geometry is inconsistent")
    if tuple(sources[:k]) != originals:
        raise AssertionError("actual originals are not the source prefix")
    return Geometry(profile, k, r, side, parent, dimension,
                    originals, recoveries, sources)


def matrix_inverse(matrix: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    work = [list(row) + [int(row_index == column) for column in range(size)]
            for row_index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise ValueError("singular matrix")
        work[column], work[pivot] = work[pivot], work[column]
        scale = GF_INVERSE[work[column][column]]
        work[column] = [GF_MULTIPLY[value][scale] for value in work[column]]
        for row in range(size):
            factor = work[row][column]
            if row == column or factor == 0:
                continue
            work[row] = [left ^ GF_MULTIPLY[factor][right]
                         for left, right in zip(work[row], work[column])]
    return tuple(tuple(row[size:]) for row in work)


def matrix_multiply(left: Sequence[Sequence[int]],
                    right: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("matrix dimensions differ")
    columns = len(right[0])
    output = []
    for left_row in left:
        row = []
        for column in range(columns):
            value = 0
            for index, coefficient in enumerate(left_row):
                value ^= GF_MULTIPLY[coefficient][right[index][column]]
            row.append(value)
        output.append(tuple(row))
    return tuple(output)


def matrix_vector(matrix: Sequence[Sequence[int]],
                  vector: Sequence[int]) -> tuple[int, ...]:
    if any(len(row) != len(vector) for row in matrix):
        raise ValueError("matrix/vector dimensions differ")
    return tuple(_xor_products(row, vector) for row in matrix)


def _xor_products(row: Sequence[int], vector: Sequence[int]) -> int:
    value = 0
    for coefficient, source in zip(row, vector):
        value ^= GF_MULTIPLY[coefficient][source]
    return value


def identity(size: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(row == column) for column in range(size))
                 for row in range(size))


def lagrange_row(source_coordinates: Sequence[int], point: int,
                 active_count: int) -> tuple[int, ...]:
    row = []
    for source_index, source in enumerate(source_coordinates[:active_count]):
        numerator = 1
        denominator = 1
        for other_index, other in enumerate(source_coordinates):
            if other_index == source_index:
                continue
            numerator = GF_MULTIPLY[numerator][point ^ other]
            denominator = GF_MULTIPLY[denominator][source ^ other]
        row.append(GF_MULTIPLY[numerator][GF_INVERSE[denominator]])
    return tuple(row)


def generator_rows(geometry: Geometry) -> dict[int, tuple[int, ...]]:
    rows = {
        coordinate: tuple(int(row == column) for column in range(geometry.k))
        for row, coordinate in enumerate(geometry.original_coordinates)
    }
    for coordinate in geometry.recovery_coordinates:
        rows[coordinate] = lagrange_row(
            geometry.source_coordinates, coordinate, geometry.k)
    if len(rows) != geometry.k + geometry.r:
        raise AssertionError("public generator rows overlap")
    return rows


def vandermonde_generator_rows(
    geometry: Geometry,
) -> dict[int, tuple[int, ...]]:
    """Independent monomial-basis oracle for the shortened parent code."""
    source_matrix = tuple(
        tuple(gf_power(point, power) for power in range(geometry.dimension))
        for point in geometry.source_coordinates
    )
    source_inverse = matrix_inverse(source_matrix)
    rows = {}
    for coordinate in geometry.public_coordinates:
        evaluation = tuple(gf_power(coordinate, power)
                           for power in range(geometry.dimension))
        full_row = matrix_multiply((evaluation,), source_inverse)[0]
        rows[coordinate] = tuple(full_row[:geometry.k])
    return rows


def ifft_prefix_butterflies(size: int, active_prefix: int) -> int:
    """Current fused DIT-IFFT schedule in radix-2 butterfly equivalents."""
    if size <= 0 or size & (size - 1):
        raise ValueError("transform size must be a power of two")
    if not 0 <= active_prefix <= size:
        raise ValueError("active prefix outside transform")
    if active_prefix == 0 or size == 1:
        return 0
    count = 0
    dist = 1
    dist4 = 4
    while dist4 <= size:
        count += ((active_prefix + dist4 - 1) // dist4) * dist4
        dist = dist4
        dist4 <<= 2
    if dist < size:
        count += size // 2
    return count


@dataclasses.dataclass(frozen=True)
class SelectionCost:
    active_blocks: int
    ifft_butterflies: int
    prefix_slots: int
    selected: tuple[int, ...]

    @property
    def metric_key(self) -> tuple[int, int, int]:
        return self.active_blocks, self.ifft_butterflies, self.prefix_slots

    @property
    def total_key(self) -> tuple[int, int, int, tuple[int, ...]]:
        return self.metric_key + (self.selected,)


def selection_cost(selected: Iterable[int], side: int) -> SelectionCost:
    coordinates = tuple(sorted(selected))
    blocks: dict[int, list[int]] = {}
    for coordinate in coordinates:
        blocks.setdefault(coordinate // side, []).append(coordinate % side)
    prefixes = [max(local) + 1 for local in blocks.values()]
    return SelectionCost(
        active_blocks=len(blocks),
        ifft_butterflies=sum(ifft_prefix_butterflies(side, prefix)
                             for prefix in prefixes),
        prefix_slots=sum(prefixes),
        selected=coordinates,
    )


def select_lowest_parent(geometry: Geometry,
                         available: Sequence[int]) -> tuple[int, ...]:
    return tuple(sorted(available)[:geometry.k])


def select_prefer_systematic(geometry: Geometry,
                             available: Sequence[int]) -> tuple[int, ...]:
    present = set(available)
    selected = [coordinate for coordinate in geometry.original_coordinates
                if coordinate in present]
    selected.extend(coordinate for coordinate in geometry.recovery_coordinates
                    if coordinate in present and len(selected) < geometry.k)
    if len(selected) != geometry.k:
        raise ValueError("not enough received coordinates")
    return tuple(sorted(selected))


def select_block_greedy(geometry: Geometry,
                        available: Sequence[int]) -> tuple[int, ...]:
    grouped: dict[int, tuple[int, ...]] = {}
    for coordinate in sorted(available):
        block = coordinate // geometry.side
        grouped.setdefault(block, tuple())
        grouped[block] += (coordinate,)
    selected: list[int] = []
    remaining = geometry.k
    unused = dict(grouped)
    while remaining:
        candidates = []
        for block, coordinates in unused.items():
            take = min(remaining, len(coordinates))
            local_selection = coordinates[:take]
            cost = selection_cost(local_selection, geometry.side)
            candidates.append((-take, cost.ifft_butterflies,
                               cost.prefix_slots, block, coordinates))
        if not candidates:
            raise ValueError("not enough received coordinates")
        _, _, _, block, coordinates = min(candidates)
        take = min(remaining, len(coordinates))
        selected.extend(coordinates[:take])
        remaining -= take
        del unused[block]
    return tuple(sorted(selected))


def select_exact_block_dp(geometry: Geometry,
                          available: Sequence[int]) -> tuple[int, ...]:
    """Exact DP for the documented lexicographic symbolic objective."""
    grouped: dict[int, tuple[int, ...]] = {}
    for coordinate in sorted(available):
        block = coordinate // geometry.side
        grouped.setdefault(block, tuple())
        grouped[block] += (coordinate,)
    # The value is the best complete SelectionCost for exactly `count` shards.
    states = {0: selection_cost((), geometry.side)}
    for block in sorted(grouped):
        coordinates = grouped[block]
        next_states: dict[int, SelectionCost] = {}
        for count, prior in states.items():
            limit = min(len(coordinates), geometry.k - count)
            for take in range(limit + 1):
                selected = prior.selected + coordinates[:take]
                candidate = selection_cost(selected, geometry.side)
                new_count = count + take
                incumbent = next_states.get(new_count)
                if incumbent is None or candidate.total_key < incumbent.total_key:
                    next_states[new_count] = candidate
        states = next_states
    if geometry.k not in states:
        raise ValueError("not enough received coordinates")
    return states[geometry.k].selected


def select(policy: str, geometry: Geometry,
           available: Sequence[int]) -> tuple[int, ...]:
    available_set = set(available)
    if len(available_set) != len(available):
        raise ValueError("received coordinates must be distinct")
    if not available_set.issubset(geometry.public_coordinates):
        raise ValueError("received coordinate is outside the public code")
    if len(available_set) < geometry.k:
        raise ValueError("not enough received coordinates")
    functions = {
        "lowest_parent_coordinate": select_lowest_parent,
        "prefer_systematic": select_prefer_systematic,
        "block_aligned_greedy": select_block_greedy,
        "exact_block_dp": select_exact_block_dp,
    }
    try:
        selected = functions[policy](geometry, available)
    except KeyError as error:
        raise ValueError("unknown selection policy") from error
    if len(selected) != geometry.k or len(set(selected)) != geometry.k:
        raise AssertionError("selection did not return K distinct coordinates")
    if not set(selected).issubset(available_set):
        raise AssertionError("selection used an unavailable coordinate")
    return selected


def subset_mask(geometry: Geometry, selected: Iterable[int]) -> int:
    indices = {coordinate: index
               for index, coordinate in enumerate(geometry.public_coordinates)}
    mask = 0
    for coordinate in selected:
        mask |= 1 << indices[coordinate]
    return mask


def all_available_patterns(geometry: Geometry):
    public = geometry.public_coordinates
    for erased_count in range(geometry.r + 1):
        for erased_indices in itertools.combinations(range(len(public)),
                                                     erased_count):
            erased = set(erased_indices)
            yield tuple(coordinate for index, coordinate in enumerate(public)
                        if index not in erased)


def proof_cell(arguments: tuple[str, int, int]) -> dict:
    profile, k, r = arguments
    geometry = make_geometry(profile, k, r)
    rows = generator_rows(geometry)
    vandermonde_rows = vandermonde_generator_rows(geometry)
    if rows != vandermonde_rows:
        raise AssertionError("Lagrange and monomial generator oracles differ")
    public = geometry.public_coordinates
    expected_identity = identity(k)
    mds_masks: set[int] = set()
    mds_subsets = 0
    basis_decode_vectors = 0
    exhaustive_message_decodes = 0
    algebra_hash = hashlib.sha256()

    for indices in itertools.combinations(range(k + r), k):
        selected_rows = tuple(rows[public[index]] for index in indices)
        inverse = matrix_inverse(selected_rows)
        if matrix_multiply(inverse, selected_rows) != expected_identity:
            raise AssertionError("direct inverse failed a basis decode")
        mask = sum(1 << index for index in indices)
        mds_masks.add(mask)
        mds_subsets += 1
        basis_decode_vectors += k
        algebra_hash.update(mask.to_bytes(2, "little"))
        algebra_hash.update(bytes(value for row in inverse for value in row))
        if k <= 2:
            for message in itertools.product(range(FIELD_ORDER), repeat=k):
                encoded = matrix_vector(selected_rows, message)
                decoded = matrix_vector(inverse, encoded)
                if decoded != message:
                    raise AssertionError("exhaustive message decode mismatch")
                exhaustive_message_decodes += 1

    policy_sums = {
        policy: {"active_blocks": 0, "ifft_butterflies": 0,
                 "prefix_slots": 0, "selected_subset_changes": 0}
        for policy in POLICIES
    }
    comparison = {
        policy: {key: 0 for key in COMPARISON_KEYS}
        for policy in POLICIES
    }
    availability_patterns = 0
    policy_decode_cases = 0
    maximum_active_block_reduction = 0
    maximum_butterfly_reduction = 0
    best_example = None
    best_example_key = (-1, -1, -1)

    for available in all_available_patterns(geometry):
        availability_patterns += 1
        costs = {}
        selections = {}
        for policy in POLICIES:
            selected = select(policy, geometry, available)
            if subset_mask(geometry, selected) not in mds_masks:
                raise AssertionError("policy selected a subset not covered by MDS proof")
            policy_decode_cases += 1
            selections[policy] = selected
            costs[policy] = selection_cost(selected, geometry.side)
            sums = policy_sums[policy]
            sums["active_blocks"] += costs[policy].active_blocks
            sums["ifft_butterflies"] += costs[policy].ifft_butterflies
            sums["prefix_slots"] += costs[policy].prefix_slots

        baseline_selection = selections["prefer_systematic"]
        for policy in POLICIES:
            if selections[policy] != baseline_selection:
                policy_sums[policy]["selected_subset_changes"] += 1

        baseline = costs["prefer_systematic"]
        exact = costs["exact_block_dp"]
        for policy, cost in costs.items():
            if cost.total_key < baseline.total_key:
                comparison[policy]["strictly_better_than_prefer_systematic"] += 1
            elif cost.total_key == baseline.total_key:
                comparison[policy]["equal_to_prefer_systematic"] += 1
            else:
                comparison[policy]["strictly_worse_than_prefer_systematic"] += 1
            if cost.total_key == exact.total_key:
                comparison[policy]["equal_to_exact"] += 1
            elif cost.total_key > exact.total_key:
                comparison[policy]["strictly_worse_than_exact"] += 1
            else:
                raise AssertionError("exact DP lost its declared objective")

            if cost.metric_key < baseline.metric_key:
                comparison[policy][
                    "metric_strictly_better_than_prefer_systematic"] += 1
            elif cost.metric_key == baseline.metric_key:
                comparison[policy]["metric_equal_to_prefer_systematic"] += 1
            else:
                comparison[policy][
                    "metric_strictly_worse_than_prefer_systematic"] += 1
            if cost.metric_key == exact.metric_key:
                comparison[policy]["metric_equal_to_exact"] += 1
                if cost.total_key > exact.total_key:
                    comparison[policy][TIEBREAK_COMPARISON_KEY] += 1
            elif cost.metric_key > exact.metric_key:
                comparison[policy]["metric_strictly_worse_than_exact"] += 1
            else:
                raise AssertionError("exact DP lost its metric objective")

        block_reduction = baseline.active_blocks - exact.active_blocks
        butterfly_reduction = baseline.ifft_butterflies - exact.ifft_butterflies
        example_key = (block_reduction, butterfly_reduction,
                       baseline.prefix_slots - exact.prefix_slots)
        maximum_active_block_reduction = max(
            maximum_active_block_reduction, block_reduction)
        maximum_butterfly_reduction = max(
            maximum_butterfly_reduction, butterfly_reduction)
        if example_key > best_example_key:
            best_example_key = example_key
            best_example = {
                "available": list(available),
                "prefer_systematic": list(selections["prefer_systematic"]),
                "exact_block_dp": list(selections["exact_block_dp"]),
                "prefer_systematic_cost": dataclasses.asdict(baseline),
                "exact_block_dp_cost": dataclasses.asdict(exact),
            }

    return {
        "profile": profile,
        "K": k,
        "R": r,
        "side": geometry.side,
        "parent": geometry.parent,
        "dimension": geometry.dimension,
        "generator_entries_compared": (geometry.k + geometry.r) * geometry.k,
        "mds_subsets": mds_subsets,
        "basis_decode_vectors": basis_decode_vectors,
        "exhaustive_message_decodes": exhaustive_message_decodes,
        "availability_patterns": availability_patterns,
        "policy_decode_cases": policy_decode_cases,
        "algebra_sha256": algebra_hash.hexdigest(),
        "policy_sums": policy_sums,
        "comparison": comparison,
        "maximum_active_block_reduction": maximum_active_block_reduction,
        "maximum_butterfly_reduction": maximum_butterfly_reduction,
        "best_example": best_example,
    }


def valid_cells() -> tuple[tuple[str, int, int], ...]:
    cells = []
    for profile in ("legacy_high_v1", "low_v1"):
        for k in range(1, FIELD_ORDER):
            for r in range(1, FIELD_ORDER + 1 - k):
                try:
                    make_geometry(profile, k, r)
                except ValueError:
                    continue
                cells.append((profile, k, r))
    return tuple(cells)


def merge(cells: Sequence[dict], source_sha256: str) -> dict:
    ordered = sorted(cells, key=lambda row: (row["profile"], row["K"], row["R"]))
    totals = {
        "profile_cells": len(ordered),
        "generator_entries_compared": sum(
            row["generator_entries_compared"] for row in ordered),
        "mds_subsets": sum(row["mds_subsets"] for row in ordered),
        "basis_decode_vectors": sum(row["basis_decode_vectors"] for row in ordered),
        "exhaustive_message_decodes": sum(
            row["exhaustive_message_decodes"] for row in ordered),
        "availability_patterns": sum(row["availability_patterns"] for row in ordered),
        "policy_decode_cases": sum(row["policy_decode_cases"] for row in ordered),
    }


    policy_totals = {}
    for policy in POLICIES:
        sums = {key: sum(row["policy_sums"][policy][key] for row in ordered)
                for key in ("active_blocks", "ifft_butterflies", "prefix_slots",
                            "selected_subset_changes")}
        comparisons = {
            key: sum(row["comparison"][policy][key] for row in ordered)
            for key in COMPARISON_KEYS
        }
        policy_totals[policy] = {"cost_sums": sums, "comparisons": comparisons}
    best = sorted(
        (row for row in ordered if row["best_example"] is not None),
        key=lambda row: (-row["maximum_active_block_reduction"],
                         -row["maximum_butterfly_reduction"],
                         row["profile"], row["K"], row["R"]),
    )[:16]
    return {
        "schema": SCHEMA,
        "source_sha256": source_sha256,
        "field": {"name": "GF(2^4)", "order": FIELD_ORDER,
                  "polynomial": "x^4+x+1", "polynomial_bits": FIELD_POLYNOMIAL},
        "scope": {
            "production_enabled": False,
            "timing_claim": False,
            "objective_order": ["active_blocks", "ifft_butterflies",
                                "prefix_slots", "lexicographic_coordinates"],
            "profiles": ["legacy_high_v1", "low_v1"],
            "maximum_parent": FIELD_ORDER,
        },
        "policies": list(POLICIES),
        "totals": totals,
        "policy_totals": policy_totals,
        "best_examples": [
            {key: row[key] for key in (
                "profile", "K", "R", "side", "parent",
                "maximum_active_block_reduction", "maximum_butterfly_reduction",
                "best_example")}
            for row in best
        ],
        "cells": ordered,
    }


def validate_checkpoint(path: Path) -> dict:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != SCHEMA:
        raise ValueError("checkpoint schema differs")
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if document.get("source_sha256") != source_hash:
        raise ValueError("checkpoint source hash differs")
    cells = document.get("cells")
    if not isinstance(cells, list):
        raise ValueError("checkpoint cells are absent")
    identities = [(row.get("profile"), row.get("K"), row.get("R"))
                  for row in cells]
    if len(identities) != len(set(identities)):
        raise ValueError("checkpoint has duplicate cells")
    if set(identities) != set(valid_cells()):
        raise ValueError("checkpoint cell coverage differs")
    for row in cells:
        geometry = make_geometry(row["profile"], row["K"], row["R"])
        expected_mds = math.comb(geometry.k + geometry.r, geometry.k)
        expected_patterns = sum(
            math.comb(geometry.k + geometry.r, erased)
            for erased in range(geometry.r + 1))
        expected = {
            "side": geometry.side,
            "parent": geometry.parent,
            "dimension": geometry.dimension,
            "generator_entries_compared": (
                geometry.k + geometry.r) * geometry.k,
            "mds_subsets": expected_mds,
            "basis_decode_vectors": expected_mds * geometry.k,
            "exhaustive_message_decodes": (
                expected_mds * FIELD_ORDER ** geometry.k
                if geometry.k <= 2 else 0),
            "availability_patterns": expected_patterns,
            "policy_decode_cases": expected_patterns * len(POLICIES),
        }
        for key, value in expected.items():
            if row.get(key) != value:
                raise ValueError(f"checkpoint cell {key} differs")
        digest = row.get("algebra_sha256")
        if (not isinstance(digest, str) or len(digest) != 64 or
                any(character not in "0123456789abcdef" for character in digest)):
            raise ValueError("checkpoint algebra digest is malformed")
        for policy in POLICIES:
            sums = row["policy_sums"][policy]
            if any(not isinstance(value, int) or value < 0
                   for value in sums.values()):
                raise ValueError("checkpoint policy sum is malformed")
            comparison = row["comparison"][policy]
            if set(comparison) != set(COMPARISON_KEYS):
                raise ValueError("checkpoint comparison keys differ")
            if any(not isinstance(value, int) or value < 0
                   for value in comparison.values()):
                raise ValueError("checkpoint comparison count is malformed")
            if sum(comparison[key] for key in (
                "strictly_better_than_prefer_systematic",
                "equal_to_prefer_systematic",
                "strictly_worse_than_prefer_systematic",
            )) != expected_patterns:
                raise ValueError("baseline comparison partition differs")
            if sum(comparison[key] for key in (
                "strictly_worse_than_exact", "equal_to_exact",
            )) != expected_patterns:
                raise ValueError("exact comparison partition differs")
            if sum(comparison[key] for key in (
                "metric_strictly_better_than_prefer_systematic",
                "metric_equal_to_prefer_systematic",
                "metric_strictly_worse_than_prefer_systematic",
            )) != expected_patterns:
                raise ValueError("metric baseline comparison partition differs")
            if sum(comparison[key] for key in (
                "metric_strictly_worse_than_exact", "metric_equal_to_exact",
            )) != expected_patterns:
                raise ValueError("metric exact comparison partition differs")
            if comparison[TIEBREAK_COMPARISON_KEY] != (
                comparison["metric_equal_to_exact"] -
                comparison["equal_to_exact"]
            ):
                raise ValueError("metric-tie/tiebreak partition differs")
            if comparison["equal_to_exact"] > comparison["metric_equal_to_exact"]:
                raise ValueError("full exact matches exceed metric matches")
            if sums["selected_subset_changes"] != (
                expected_patterns -
                comparison["equal_to_prefer_systematic"]
            ):
                raise ValueError("selection-change/full-objective partition differs")
            if policy == "prefer_systematic" and any((
                comparison["strictly_better_than_prefer_systematic"] != 0,
                comparison["equal_to_prefer_systematic"] != expected_patterns,
                comparison["strictly_worse_than_prefer_systematic"] != 0,
                comparison[
                    "metric_strictly_better_than_prefer_systematic"] != 0,
                comparison["metric_equal_to_prefer_systematic"] !=
                expected_patterns,
                comparison[
                    "metric_strictly_worse_than_prefer_systematic"] != 0,
            )):
                raise ValueError("baseline self-comparison differs")
            if policy == "exact_block_dp" and any((
                comparison["strictly_worse_than_exact"] != 0,
                comparison["equal_to_exact"] != expected_patterns,
                comparison["metric_strictly_worse_than_exact"] != 0,
                comparison["metric_equal_to_exact"] != expected_patterns,
                comparison[TIEBREAK_COMPARISON_KEY] != 0,
            )):
                raise ValueError("exact self-comparison differs")
    if document != merge(cells, source_hash):
        raise ValueError("checkpoint aggregates or ordering differ")
    return {
        **document["totals"],
        "checkpoint_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def validate_field() -> int:
    checks = 0
    for left in range(FIELD_ORDER):
        for right in range(FIELD_ORDER):
            if GF_MULTIPLY[left][right] != GF_MULTIPLY[right][left]:
                raise AssertionError("field multiplication is not commutative")
            checks += 1
            for third in range(FIELD_ORDER):
                if GF_MULTIPLY[left][right ^ third] != (
                    GF_MULTIPLY[left][right] ^ GF_MULTIPLY[left][third]
                ):
                    raise AssertionError("field multiplication is not distributive")
                checks += 1
    for value in range(1, FIELD_ORDER):
        if GF_MULTIPLY[value][GF_INVERSE[value]] != 1:
            raise AssertionError("field inverse is wrong")
        checks += 1
    return checks


def brute_force_select(geometry: Geometry,
                       available: Sequence[int]) -> tuple[int, ...]:
    return min(
        (selection_cost(candidate, geometry.side)
         for candidate in itertools.combinations(sorted(available), geometry.k)),
        key=lambda cost: cost.total_key,
    ).selected


def self_test() -> dict:
    field_checks = validate_field()
    dp_checks = 0
    for profile, k, r in (
        ("legacy_high_v1", 1, 3),
        ("legacy_high_v1", 3, 5),
        ("low_v1", 3, 5),
        ("low_v1", 5, 3),
    ):
        geometry = make_geometry(profile, k, r)
        public = geometry.public_coordinates
        for mask in range(1 << len(public)):
            available = tuple(public[index] for index in range(len(public))
                              if mask & (1 << index))
            if len(available) < k:
                continue
            if select_exact_block_dp(geometry, available) != brute_force_select(
                geometry, available
            ):
                raise AssertionError("exact DP differs from exhaustive search")
            for policy in POLICIES:
                select(policy, geometry, available)
            dp_checks += 1
    representative = [proof_cell(cell) for cell in (
        ("legacy_high_v1", 4, 4), ("low_v1", 4, 4))]
    if any(row["mds_subsets"] != 70 for row in representative):
        raise AssertionError("representative MDS subset count is wrong")
    return {
        "schema": "leopard2-received-subset-self-test-v1",
        "field_checks": field_checks,
        "dp_vs_bruteforce_patterns": dp_checks,
        "representative_mds_subsets": sum(row["mds_subsets"]
                                           for row in representative),
        "representative_policy_decode_cases": sum(
            row["policy_decode_cases"] for row in representative),
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def run(output: Path, workers: int) -> dict:
    if workers <= 0:
        raise ValueError("workers must be positive")
    cells = valid_cells()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(proof_cell, cells, chunksize=1))
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    checkpoint = merge(results, source_hash)
    write_json(output, checkpoint)
    return {
        "output": str(output),
        "workers": workers,
        "allowed_cpus": len(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity") else os.cpu_count(),
        **checkpoint["totals"],
        "checkpoint_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("self-test")
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument(
        "--input", type=Path,
        default=Path(__file__).resolve().parent / "results/checkpoint.json")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--workers", type=int,
                            default=min(128, len(os.sched_getaffinity(0)))
                            if hasattr(os, "sched_getaffinity") else
                            min(128, os.cpu_count() or 1))
    run_parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).resolve().parent / "results/checkpoint.json")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_args(sys.argv[1:] if argv is None else argv)
    if arguments.command == "self-test":
        print(json.dumps(self_test(), indent=2, sort_keys=True))
    elif arguments.command == "validate":
        print(json.dumps(validate_checkpoint(arguments.input),
                         indent=2, sort_keys=True))
    else:
        print(json.dumps(run(arguments.output, arguments.workers),
                         indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
