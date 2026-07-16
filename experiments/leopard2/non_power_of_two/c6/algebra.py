#!/usr/bin/env python3
"""Independent algebra checkpoint for C6 exact-GF256 boundary rescue.

This program does not call Leopard.  It uses carryless polynomial arithmetic
behind the declared coordinate bases to check the proposed exact prefix code:

    systematic evaluations: 0 .. K-1
    parity evaluations:      K .. K+R-1

For non-dyadic K this is a new code profile.  It is never described as a
wire-compatible implementation of a padded GF16 parent.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import importlib.util
import itertools
import json
import random
import struct
import sys
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[4]
C3B = ROOT / "experiments/leopard2/non_power_of_two/c3b/fast_inverse.py"
SCHEMA = "leopard2-c6-algebra/v2"


def load_c3b():
    name = "_leopard2_c6_c3b_reference"
    spec = importlib.util.spec_from_file_location(name, C3B)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load C3b independent field reference")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


c3b = load_c3b()


def ceil_pow2(value: int) -> int:
    if value <= 0:
        raise ValueError("count must be positive")
    return 1 << (value - 1).bit_length()


def parent_size(profile: str, k: int, r: int) -> int:
    if profile == "legacy_high_v1":
        return ceil_pow2(k + ceil_pow2(r))
    if profile == "low_v1":
        return ceil_pow2(ceil_pow2(k) + r)
    raise ValueError("unknown profile")


def affected(profile: str, k: int, r: int) -> bool:
    return k > 0 and r > 0 and k + r <= 256 and parent_size(profile, k, r) > 256


def generator_parity(field, k: int, r: int) -> tuple[tuple[int, ...], ...]:
    """Direct Lagrange rows for exact prefix evaluations."""
    weights = []
    for i in range(k):
        denominator = 1
        for other in range(k):
            if other != i:
                denominator = field.multiply(denominator, i ^ other)
        weights.append(field.inverse(denominator))
    rows = []
    for parity in range(r):
        point = k + parity
        vanishing = 1
        for systematic in range(k):
            vanishing = field.multiply(vanishing, point ^ systematic)
        rows.append(tuple(
            field.multiply(
                field.multiply(vanishing, field.inverse(point ^ systematic)),
                weights[systematic],
            )
            for systematic in range(k)
        ))
    return tuple(rows)


def rank(field, matrix: Sequence[Sequence[int]]) -> int:
    if not matrix:
        return 0
    rows = [list(row) for row in matrix]
    columns = len(rows[0])
    if any(len(row) != columns for row in rows):
        raise ValueError("ragged matrix")
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, len(rows))
                      if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        inverse = field.inverse(rows[pivot_row][column])
        rows[pivot_row] = [field.multiply(value, inverse)
                           for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [left ^ field.multiply(factor, right)
                         for left, right in zip(rows[row], rows[pivot_row])]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def invert_matrix(field, matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    augmented = [list(row) + [int(row_index == column)
                              for column in range(size)]
                 for row_index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if augmented[row][column]), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = field.inverse(augmented[column][column])
        augmented[column] = [field.multiply(value, inverse)
                             for value in augmented[column]]
        for row in range(size):
            if row == column or augmented[row][column] == 0:
                continue
            factor = augmented[row][column]
            augmented[row] = [left ^ field.multiply(factor, right)
                              for left, right in zip(augmented[row],
                                                     augmented[column])]
    return [row[size:] for row in augmented]


def deterministic_missing(k: int, losses: int) -> tuple[int, ...]:
    if not 0 < losses <= k:
        raise ValueError("invalid loss count")
    result = tuple((index * k) // losses for index in range(losses))
    if len(set(result)) != losses:
        raise AssertionError("deterministic loss set contains a duplicate")
    return result


def direct_decode_terms(field, rows: Sequence[Sequence[int]], k: int,
                        missing: Sequence[int]) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    losses = len(missing)
    selected = tuple(range(losses))
    inverse = invert_matrix(field, [[rows[parity][original]
                                     for original in missing]
                                    for parity in selected])
    missing_set = set(missing)
    outputs = []
    for output in range(losses):
        terms = []
        for equation, parity in enumerate(selected):
            coefficient = inverse[output][equation]
            if coefficient:
                terms.append((1, parity, coefficient))
        for original in range(k):
            if original in missing_set:
                continue
            coefficient = 0
            for equation, parity in enumerate(selected):
                coefficient ^= field.multiply(inverse[output][equation],
                                              rows[parity][original])
            if coefficient:
                terms.append((0, original, coefficient))
        if not terms:
            raise AssertionError("decode output contains no terms")
        outputs.append(tuple(terms))
    return tuple(outputs)


def systematic_generator(field, k: int, r: int) -> tuple[tuple[int, ...], ...]:
    identity = tuple(tuple(int(row == column) for column in range(k))
                     for row in range(k))
    return identity + generator_parity(field, k, r)


def apply(field, matrix: Sequence[Sequence[int]], values: Sequence[int]) -> list[int]:
    output = []
    for row in matrix:
        value = 0
        for coefficient, source in zip(row, values):
            if coefficient and source:
                value ^= field.multiply(coefficient, source)
        output.append(value)
    return output


def lagrange_evaluate(field, values: Sequence[int], point: int) -> int:
    value = 0
    for i, source in enumerate(values):
        numerator = 1
        denominator = 1
        for other in range(len(values)):
            if other == i:
                continue
            numerator = field.multiply(numerator, point ^ other)
            denominator = field.multiply(denominator, i ^ other)
        value ^= field.multiply(source,
                                field.multiply(numerator, field.inverse(denominator)))
    return value


def exact_identity_digest(k: int, r: int) -> str:
    canonical = bytearray(b"LEO2-EXACT-PREFIX-GF8-UNFROZEN\0")
    canonical += struct.pack(">II", k, r)
    for role, points in ((1, range(k)), (2, range(k, k + r))):
        canonical += struct.pack(">BI", role, len(points))
        canonical += bytes(points)
    return hashlib.sha256(canonical).hexdigest()


def legacy_gf16_coefficient(field, profile: str, k: int, r: int,
                            systematic: int, parity: int) -> int:
    if profile == "legacy_high_v1":
        side = ceil_pow2(r)
        parent = ceil_pow2(k + side)
        systematic_points = range(side, parent)
        source_point = side + systematic
        parity_point = parity
    elif profile == "low_v1":
        side = ceil_pow2(k)
        parent = ceil_pow2(side + r)
        systematic_points = range(side)
        source_point = systematic
        parity_point = side + parity
    else:
        raise ValueError("unknown profile")
    numerator = 1
    denominator = 1
    for point in systematic_points:
        if point != source_point:
            numerator = field.multiply(numerator, parity_point ^ point)
            denominator = field.multiply(denominator, source_point ^ point)
    return field.multiply(numerator, field.inverse(denominator))


def exhaustive_gf4() -> dict[str, int]:
    field = c3b.field_named("gf4")
    subsets = 0
    generators = 0
    for k in range(1, 16):
        generator = systematic_generator(field, k, 16 - k)
        generators += 1
        for selected in itertools.combinations(range(16), k):
            subsets += 1
            if rank(field, [generator[row] for row in selected]) != k:
                raise AssertionError(f"GF4 non-MDS subset: K={k} rows={selected}")
    if subsets != 65534:
        raise AssertionError(f"unexpected GF4 subset count {subsets}")
    return {"full_length_generators": generators,
            "k_coordinate_subsets": subsets,
            "implied_public_geometries": 120}


def boundary_cases() -> tuple[tuple[str, int, int], ...]:
    high = ((1, 129), (3, 129), (17, 129), (31, 129), (63, 129),
            (65, 129), (127, 129), (1, 193), (3, 193), (17, 193),
            (31, 193), (63, 193))
    low = tuple((k, r) for r, k in high)
    return tuple(("legacy_high_v1", k, r) for k, r in high) + tuple(
        ("low_v1", k, r) for k, r in low)


def deterministic_indices(count: int, amount: int, salt: int) -> tuple[int, ...]:
    if amount > count:
        raise ValueError("sample exceeds population")
    rng = random.Random((count << 20) ^ (amount << 8) ^ salt)
    return tuple(sorted(rng.sample(range(count), amount)))


def validate_gf8_boundaries() -> tuple[list[dict], dict[str, int]]:
    gf8 = c3b.field_named("gf8")
    gf16 = c3b.field_named("gf16")
    records = []
    rank_checks = 0
    coefficient_checks = 0
    evaluation_checks = 0
    witnesses = 0
    for case_index, (profile, k, r) in enumerate(boundary_cases()):
        if not affected(profile, k, r):
            raise AssertionError(f"case is not a field-boundary case: {profile} {k} {r}")
        parity_rows = generator_parity(gf8, k, r)
        if any(coefficient == 0 for row in parity_rows for coefficient in row):
            raise AssertionError("exact prefix RS parity row contains zero")
        coefficient_checks += k * r

        local_rank_checks = 0
        for loss in sorted(set((1, min(2, k, r), min(4, k, r), min(8, k, r)))):
            missing = deterministic_indices(k, loss, case_index * 17 + loss)
            parity = deterministic_indices(r, loss, case_index * 29 + loss)
            square = [[parity_rows[row][column] for column in missing]
                      for row in parity]
            if rank(gf8, square) != loss:
                raise AssertionError("GF8 repair minor is singular")
            local_rank_checks += 1
        rank_checks += local_rank_checks

        values = [((index * 73 + case_index * 19 + 1) & 255) for index in range(k)]
        matrix_values = apply(gf8, parity_rows, values)
        for parity in sorted(set((0, r // 2, r - 1))):
            expected = lagrange_evaluate(gf8, values, k + parity)
            if matrix_values[parity] != expected:
                raise AssertionError("GF8 matrix/direct evaluation mismatch")
            evaluation_checks += 1

        witness = None
        for parity in sorted(set((0, r // 2, r - 1))):
            for systematic in sorted(set((0, k // 2, k - 1))):
                exact = parity_rows[parity][systematic]
                padded = legacy_gf16_coefficient(
                    gf16, profile, k, r, systematic, parity)
                if exact != padded:
                    witness = {
                        "systematic": systematic,
                        "parity": parity,
                        "exact_gf8_coefficient": exact,
                        "padded_gf16_coefficient": padded,
                    }
                    break
            if witness:
                break
        if witness is None:
            raise AssertionError("no explicit field/profile difference witness")
        witnesses += 1
        records.append({
            "profile": profile,
            "K": k,
            "R": r,
            "transmitted": k + r,
            "padded_parent": parent_size(profile, k, r),
            "exact_coordinate_digest_sha256": exact_identity_digest(k, r),
            "parity_coefficients": k * r,
            "repair_minor_rank_checks": local_rank_checks,
            "wire_difference_witness": witness,
        })
    return records, {
        "boundary_cases": len(records),
        "parity_coefficients_nonzero": coefficient_checks,
        "repair_minor_rank_checks": rank_checks,
        "direct_evaluation_checks": evaluation_checks,
        "explicit_wire_difference_witnesses": witnesses,
    }


def decoder_patterns() -> tuple[tuple[str, int, int, int], ...]:
    return (
        ("legacy_high_v1", 2, 129, 2),
        ("legacy_high_v1", 3, 129, 1),
        ("legacy_high_v1", 3, 129, 3),
        ("legacy_high_v1", 3, 130, 3),
        ("legacy_high_v1", 4, 129, 4),
        ("legacy_high_v1", 17, 129, 1),
        ("legacy_high_v1", 17, 129, 4),
        ("legacy_high_v1", 17, 129, 17),
        ("low_v1", 129, 2, 2),
        ("low_v1", 129, 3, 1),
        ("low_v1", 129, 3, 3),
        ("low_v1", 130, 3, 3),
        ("low_v1", 129, 4, 4),
        ("low_v1", 129, 17, 1),
        ("low_v1", 129, 17, 4),
        ("low_v1", 129, 17, 17),
    )


def validate_decoder_plans() -> tuple[list[dict], dict[str, int]]:
    field = c3b.field_named("gf8")
    records = []
    total_terms = 0
    recovered_values = 0
    for pattern_index, (profile, k, r, losses) in enumerate(decoder_patterns()):
        rows = generator_parity(field, k, r)
        missing = deterministic_missing(k, losses)
        outputs = direct_decode_terms(field, rows, k, missing)
        message = [((index * 67 + pattern_index * 31 + 9) & 255)
                   for index in range(k)]
        parity = apply(field, rows, message)
        for output_index, terms in enumerate(outputs):
            recovered = 0
            for is_parity, index, coefficient in terms:
                source = parity[index] if is_parity else message[index]
                recovered ^= field.multiply(coefficient, source)
            if recovered != message[missing[output_index]]:
                raise AssertionError("independent direct decoder recovered incorrectly")
            recovered_values += 1
        encoded = bytearray()
        for terms in outputs:
            for is_parity, index, coefficient in terms:
                encoded += struct.pack(">BHB", is_parity, index, coefficient)
        term_count = sum(len(terms) for terms in outputs)
        total_terms += term_count
        records.append({
            "profile": profile,
            "K": k,
            "R": r,
            "losses": losses,
            "missing_originals": list(missing),
            "selected_parities": list(range(losses)),
            "term_count": term_count,
            "term_digest_sha256": hashlib.sha256(encoded).hexdigest(),
        })
    return records, {
        "plans": len(records),
        "folded_nonzero_terms": total_terms,
        "recovered_values": recovered_values,
        "no_loss_execution_terms": 0,
    }


def prior_artifact(path: str) -> dict[str, object]:
    full = ROOT / path
    data = full.read_bytes()
    return {"path": path, "sha256": hashlib.sha256(data).hexdigest()}


def run() -> dict:
    affected_counts = {}
    for profile in ("legacy_high_v1", "low_v1"):
        affected_counts[profile] = sum(
            affected(profile, k, r)
            for k in range(1, 256) for r in range(1, 257 - k)
        )
    if affected_counts != {"legacy_high_v1": 10795, "low_v1": 10795}:
        raise AssertionError(f"affected geometry changed: {affected_counts}")
    records, gf8_counts = validate_gf8_boundaries()
    decoder_records, decoder_counts = validate_decoder_plans()
    return {
        "schema": SCHEMA,
        "status": "pass",
        "code_definition": {
            "field": "legacy_gf8_representation_v1",
            "systematic_coordinates": "0..K-1",
            "parity_coordinates": "K..K+R-1",
            "degree_bound": "degree < K",
            "identity_status": "unfrozen-new-profile; no serialized production identity",
            "legacy_wire_compatible": False,
        },
        "affected_profile_cells": affected_counts,
        "gf4_exhaustive": exhaustive_gf4(),
        "gf8": gf8_counts,
        "boundary_records": records,
        "decoder": decoder_counts,
        "decoder_plan_records": decoder_records,
        "method_scope": {
            "legacy_padded_gf16": "measured by the C6 C++ public-API baseline",
            "exact_direct_gf8": "implemented and measured by the C6 C++ candidate",
            "wire_compatible_truncated_gf8": (
                "impossible in the target cells: the declared parent has 512 "
                "GF16 coordinates, exceeding GF8 order; C2 remains a same-field parent executor"
            ),
            "tang_han_epsilon_gf8": (
                "C3b exact inverse algebra only; scalar best gain 0.862 and no "
                "end-to-end exact encoder/profile identity"
            ),
            "binary_block_gf8": (
                "C5 found the materialized join dense/rejected; a fully fused "
                "exact join is the direct evaluator measured here"
            ),
            "generic_gf8": (
                "same exact prefix polynomial and wire definition; generic "
                "per-stripe interpolation is retained as an oracle, not relabeled as a distinct plan"
            ),
        },
        "prior_evidence": [
            prior_artifact("experiments/leopard2/c2_truncated_cpp/results/checkpoint.json"),
            prior_artifact("experiments/leopard2/non_power_of_two/c3b/results/checkpoint.json"),
            prior_artifact("experiments/leopard2/c5_dyadic_cpp/results/checkpoint.json"),
            prior_artifact("experiments/leopard2/non_power_of_two/c0/results/summary.json"),
        ],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
