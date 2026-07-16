#!/usr/bin/env python3
"""Independent algebra for the C7 exact-prefix low-rate profile.

This file intentionally does not import Leopard or C6.  It constructs binary
extension fields from their polynomial and coordinate-basis definitions, then
derives the systematic generator in two independent ways:

* barycentric Lagrange products; and
* inversion of a monomial Vandermonde matrix.

The exhaustive GF(2^4) checks cover MDS, density, affine invariance, and a
bounded coordinate-set search.  GF8/GF16 samples bind the resulting profile to
Leopard's declared legacy field representations without using production
tables as an oracle.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


SCHEMA = "leopard2-c7-algebra/v1"


@dataclass(frozen=True)
class Field:
    name: str
    bits: int
    polynomial: int
    coordinate_basis: tuple[int, ...]

    def __post_init__(self) -> None:
        if len(self.coordinate_basis) != self.bits:
            raise ValueError("coordinate basis has the wrong width")
        order = 1 << self.bits
        seen = set()
        for coordinate in range(order):
            value = self.to_polynomial(coordinate)
            if value >= order or value in seen:
                raise ValueError("coordinate basis is not independent")
            seen.add(value)

    @property
    def order(self) -> int:
        return 1 << self.bits

    @functools.lru_cache(maxsize=None)
    def to_polynomial(self, coordinate: int) -> int:
        value = 0
        for bit, basis in enumerate(self.coordinate_basis):
            if coordinate & (1 << bit):
                value ^= basis
        return value

    @functools.lru_cache(maxsize=None)
    def from_polynomial(self, polynomial: int) -> int:
        try:
            return coordinate_reverse(self)[polynomial]
        except KeyError as error:
            raise AssertionError(
                "polynomial element is outside the coordinate map") from error

    @functools.lru_cache(maxsize=None)
    def multiply(self, left: int, right: int) -> int:
        a = self.to_polynomial(left)
        b = self.to_polynomial(right)
        product = 0
        while b:
            if b & 1:
                product ^= a
            b >>= 1
            a <<= 1
        for bit in range(2 * self.bits - 2, self.bits - 1, -1):
            if product & (1 << bit):
                product ^= self.polynomial << (bit - self.bits)
        return self.from_polynomial(product)

    @functools.lru_cache(maxsize=None)
    def power(self, value: int, exponent: int) -> int:
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            exponent >>= 1
            if exponent:
                value = self.multiply(value, value)
        return result

    @functools.lru_cache(maxsize=None)
    def inverse(self, value: int) -> int:
        if value == 0:
            raise ZeroDivisionError("inverse of zero")
        return self.power(value, self.order - 2)


@functools.lru_cache(maxsize=None)
def coordinate_reverse(field: Field) -> dict[int, int]:
    # Built from the declared basis, independently of Leopard's log tables.
    return {field.to_polynomial(coordinate): coordinate
            for coordinate in range(field.order)}


GF4 = Field("gf4-test", 4, 0x13, (1, 2, 4, 8))
GF8 = Field("legacy-gf8-v1", 8, 0x11D,
            (1, 214, 152, 146, 86, 200, 88, 230))
GF16 = Field("legacy-gf16-v1", 16, 0x1002D,
             (0x0001, 0xACCA, 0x3C0E, 0x163E,
              0xC582, 0xED2E, 0x914C, 0x4012,
              0x6C98, 0x10D8, 0x6A72, 0xB900,
              0xFDB8, 0xFB34, 0xFF38, 0x991E))


def matrix_inverse(field: Field, matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    rows = [list(row) + [int(row_index == column) for column in range(size)]
            for row_index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if rows[row][column]), None)
        if pivot is None:
            raise ValueError("singular matrix")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = field.inverse(rows[column][column])
        rows[column] = [field.multiply(value, scale) for value in rows[column]]
        for row in range(size):
            if row == column or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [left ^ field.multiply(factor, right)
                         for left, right in zip(rows[row], rows[column])]
    return [row[size:] for row in rows]


def matrix_rank(field: Field, matrix: Sequence[Sequence[int]]) -> int:
    if not matrix:
        return 0
    rows = [list(row) for row in matrix]
    columns = len(rows[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, len(rows))
                      if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = field.inverse(rows[pivot_row][column])
        rows[pivot_row] = [field.multiply(value, scale)
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


def barycentric_rows(field: Field, systematic: Sequence[int],
                     parity: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    if len(set(systematic) | set(parity)) != len(systematic) + len(parity):
        raise ValueError("evaluation coordinates must be distinct")
    weights = []
    for index, point in enumerate(systematic):
        denominator = 1
        for other_index, other in enumerate(systematic):
            if other_index != index:
                denominator = field.multiply(denominator, point ^ other)
        weights.append(field.inverse(denominator))
    rows = []
    for point in parity:
        vanishing = 1
        for source in systematic:
            vanishing = field.multiply(vanishing, point ^ source)
        rows.append(tuple(
            field.multiply(
                field.multiply(vanishing, field.inverse(point ^ source)),
                weights[index],
            )
            for index, source in enumerate(systematic)
        ))
    return tuple(rows)


def vandermonde_rows(field: Field, systematic: Sequence[int],
                     parity: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Independent monomial interpolation/evaluation construction."""
    size = len(systematic)
    vandermonde = []
    for point in systematic:
        row = [1]
        for _ in range(1, size):
            row.append(field.multiply(row[-1], point))
        vandermonde.append(row)
    inverse = matrix_inverse(field, vandermonde)
    rows = []
    for point in parity:
        powers = [1]
        for _ in range(1, size):
            powers.append(field.multiply(powers[-1], point))
        rows.append(tuple(
            reduce_xor(field.multiply(powers[degree], inverse[degree][source])
                       for degree in range(size))
            for source in range(size)
        ))
    return tuple(rows)


def reduce_xor(values: Iterable[int]) -> int:
    result = 0
    for value in values:
        result ^= value
    return result


def systematic_generator(rows: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    k = len(rows[0]) if rows else 0
    identity = tuple(tuple(int(row == column) for column in range(k))
                     for row in range(k))
    return identity + tuple(tuple(row) for row in rows)


def affine_points(field: Field, points: Sequence[int], scale: int,
                  translation: int) -> tuple[int, ...]:
    if scale == 0:
        raise ValueError("affine scale must be nonzero")
    return tuple(field.multiply(scale, point) ^ translation for point in points)


def is_global_affine_image(field: Field, target: Sequence[int]) -> bool:
    """Whether ordered target points equal phi(0),phi(1),... for one phi."""
    if len(target) < 2:
        return True
    translation = target[0]
    scale = target[1] ^ translation
    if scale == 0:
        return False
    return affine_points(field, tuple(range(len(target))), scale,
                         translation) == tuple(target)


def aligned_union_points(k: int, r: int, order: int) -> tuple[tuple[int, ...],
                                                               tuple[int, ...]] | None:
    """Exact low map obtained by removing the padded-low shortened interval.

    Systematic points remain the prefix 0..K-1.  Parity begins at
    P=ceil_pow2(K), an aligned additive coset boundary.  It is a real
    non-prefix union when K is not a power of two and is unavailable if P+R
    crosses the field boundary.
    """
    p = 1 << (k - 1).bit_length()
    if p + r > order:
        return None
    return tuple(range(k)), tuple(range(p, p + r))


def dyadic_fragment_count(points: Sequence[int], order: int) -> int:
    """Minimal complete aligned-dyadic-block cover of an exact point set."""
    selected = set(points)

    def visit(start: int, size: int) -> int:
        present = sum(point in selected for point in range(start, start + size))
        if present == 0:
            return 0
        if present == size:
            return 1
        if size == 1:
            return 1
        half = size // 2
        return visit(start, half) + visit(start + half, half)

    return visit(0, order)


def coordinate_digest(field: Field, systematic: Sequence[int],
                      parity: Sequence[int]) -> str:
    encoded = bytearray(b"LEO2-EXACT-LOW-PREFIX-V1\0")
    encoded += struct.pack(">BII", field.bits, len(systematic), len(parity))
    width = (field.bits + 7) // 8
    for role, points in ((1, systematic), (2, parity)):
        encoded += struct.pack(">BI", role, len(points))
        for point in points:
            encoded += point.to_bytes(width, "big")
    return hashlib.sha256(encoded).hexdigest()


def exhaustive_small_field() -> dict[str, object]:
    subsets = 0
    geometries = 0
    dense_coefficients = 0
    oracle_coefficients = 0
    affine_maps = 0
    affine_coefficients = 0
    aligned_available = 0
    aligned_affine = 0
    aligned_changed_coefficients = 0
    aligned_dense_coefficients = 0
    aligned_prefix_one_coefficients = 0
    aligned_union_one_coefficients = 0
    aligned_prefix_fragments = 0
    aligned_union_fragments = 0
    for transmitted in range(2, GF4.order + 1):
        points = tuple(range(transmitted))
        for k in range(1, transmitted):
            geometries += 1
            systematic = points[:k]
            parity = points[k:]
            rows = barycentric_rows(GF4, systematic, parity)
            oracle = vandermonde_rows(GF4, systematic, parity)
            if rows != oracle:
                raise AssertionError("barycentric/Vandermonde mismatch")
            if any(coefficient == 0 for row in rows for coefficient in row):
                raise AssertionError("exact parity generator is not dense")
            dense_coefficients += k * (transmitted - k)
            oracle_coefficients += k * (transmitted - k)
            generator = systematic_generator(rows)
            for selected in itertools.combinations(range(transmitted), k):
                subsets += 1
                if matrix_rank(GF4, [generator[index] for index in selected]) != k:
                    raise AssertionError("exact-low small-field MDS failure")

            # Exhaust every invertible affine map of GF(16), including maps
            # that move this prefix onto non-prefix and aligned-coset points.
            for scale in range(1, GF4.order):
                for translation in range(GF4.order):
                    transformed_systematic = affine_points(
                        GF4, systematic, scale, translation)
                    transformed_parity = affine_points(
                        GF4, parity, scale, translation)
                    transformed = barycentric_rows(
                        GF4, transformed_systematic, transformed_parity)
                    if transformed != rows:
                        raise AssertionError("affine generator invariance failed")
                    affine_maps += 1
                    affine_coefficients += k * (transmitted - k)

            aligned = aligned_union_points(k, transmitted - k, GF4.order)
            if aligned is not None:
                aligned_available += 1
                aligned_systematic, aligned_parity = aligned
                aligned_rows = barycentric_rows(
                    GF4, aligned_systematic, aligned_parity)
                if any(value == 0 for row in aligned_rows for value in row):
                    raise AssertionError("aligned-union generator is not dense")
                aligned_dense_coefficients += k * (transmitted - k)
                target = aligned_systematic + aligned_parity
                globally_affine = is_global_affine_image(GF4, target)
                aligned_affine += int(globally_affine)
                if globally_affine and aligned_rows != rows:
                    raise AssertionError("affine aligned-union changed generator")
                if not globally_affine and aligned_rows != rows:
                    aligned_changed_coefficients += sum(
                        left != right
                        for left_row, right_row in zip(rows, aligned_rows)
                        for left, right in zip(left_row, right_row)
                    )
                aligned_prefix_one_coefficients += sum(
                    value == 1 for row in rows for value in row)
                aligned_union_one_coefficients += sum(
                    value == 1 for row in aligned_rows for value in row)
                aligned_prefix_fragments += (
                    dyadic_fragment_count(systematic, GF4.order) +
                    dyadic_fragment_count(parity, GF4.order)
                )
                aligned_union_fragments += (
                    dyadic_fragment_count(aligned_systematic, GF4.order) +
                    dyadic_fragment_count(aligned_parity, GF4.order)
                )

    # Exhaust every full-field systematic/parity partition.  A searched set can
    # change how many coefficients equal one, but distinct evaluation points
    # cannot create zero coefficients: all candidates remain dense K*R maps.
    search = []
    all_points = tuple(range(GF4.order))
    searched_partitions = 0
    for k in range(1, GF4.order):
        canonical_rows = barycentric_rows(GF4, all_points[:k], all_points[k:])
        canonical_ones = sum(value == 1 for row in canonical_rows for value in row)
        best_ones = -1
        best_systematic: tuple[int, ...] | None = None
        for systematic in itertools.combinations(all_points, k):
            parity = tuple(point for point in all_points if point not in systematic)
            rows = barycentric_rows(GF4, systematic, parity)
            searched_partitions += 1
            if any(value == 0 for row in rows for value in row):
                raise AssertionError("searched generator unexpectedly sparse")
            ones = sum(value == 1 for row in rows for value in row)
            if ones > best_ones:
                best_ones = ones
                best_systematic = systematic
        assert best_systematic is not None
        search.append({
            "K": k,
            "R": GF4.order - k,
            "canonical_one_coefficients": canonical_ones,
            "best_one_coefficients": best_ones,
            "best_systematic_coordinates": list(best_systematic),
            "dense_coefficients": k * (GF4.order - k),
        })

    # Offline transform-aware search keeps the frozen systematic prefix but
    # considers every nonempty parity subset of the remaining GF(16) points.
    # It first minimizes complete dyadic fragments, then maximizes coefficient
    # ones, then uses lexicographic coordinates for deterministic tie-breaking.
    transform_search_candidates = 0
    transform_search_wins = 0
    transform_search_focus = []
    for k in range(1, GF4.order):
        systematic = tuple(range(k))
        remaining = tuple(range(k, GF4.order))
        for r in range(1, len(remaining) + 1):
            prefix_parity = remaining[:r]
            prefix_rows = barycentric_rows(GF4, systematic, prefix_parity)
            prefix_score = (
                dyadic_fragment_count(prefix_parity, GF4.order),
                -sum(value == 1 for row in prefix_rows for value in row),
                prefix_parity,
            )
            best_score = prefix_score
            best_parity = prefix_parity
            for parity in itertools.combinations(remaining, r):
                rows = barycentric_rows(GF4, systematic, parity)
                transform_search_candidates += 1
                if any(value == 0 for row in rows for value in row):
                    raise AssertionError("searched parity generator is not dense")
                score = (
                    dyadic_fragment_count(parity, GF4.order),
                    -sum(value == 1 for row in rows for value in row),
                    parity,
                )
                if score < best_score:
                    best_score = score
                    best_parity = parity
            if best_parity != prefix_parity:
                transform_search_wins += 1
            if (k, r) in ((3, 3), (3, 5), (5, 3), (7, 5), (9, 3)):
                aligned = aligned_union_points(k, r, GF4.order)
                aligned_record = None
                if aligned is not None:
                    aligned_rows = barycentric_rows(GF4, *aligned)
                    aligned_record = {
                        "parity": list(aligned[1]),
                        "global_affine": is_global_affine_image(
                            GF4, aligned[0] + aligned[1]),
                        "fragments": dyadic_fragment_count(
                            aligned[1], GF4.order),
                        "one_coefficients": sum(
                            value == 1 for row in aligned_rows for value in row),
                    }
                transform_search_focus.append({
                    "K": k,
                    "R": r,
                    "prefix": {
                        "parity": list(prefix_parity),
                        "fragments": prefix_score[0],
                        "one_coefficients": -prefix_score[1],
                    },
                    "aligned_union": aligned_record,
                    "searched": {
                        "parity": list(best_parity),
                        "fragments": best_score[0],
                        "one_coefficients": -best_score[1],
                    },
                })

    return {
        "geometries": geometries,
        "mds_subsets": subsets,
        "dense_coefficients": dense_coefficients,
        "vandermonde_oracle_coefficients": oracle_coefficients,
        "affine_maps": affine_maps,
        "affine_coefficients": affine_coefficients,
        "aligned_union": {
            "available_geometries": aligned_available,
            "globally_affine_geometries": aligned_affine,
            "nonaffine_changed_coefficients": aligned_changed_coefficients,
            "dense_coefficients": aligned_dense_coefficients,
            "prefix_one_coefficients": aligned_prefix_one_coefficients,
            "aligned_one_coefficients": aligned_union_one_coefficients,
            "prefix_symbolic_fragments": aligned_prefix_fragments,
            "aligned_symbolic_fragments": aligned_union_fragments,
        },
        "searched_full_field_partitions": searched_partitions,
        "search": search,
        "transform_search_candidates": transform_search_candidates,
        "transform_search_wins": transform_search_wins,
        "transform_search_focus": transform_search_focus,
    }


def sampled_large_fields() -> dict[str, object]:
    cases = (
        (GF8, 3, 253), (GF8, 7, 249), (GF8, 16, 240),
        (GF8, 64, 192), (GF8, 127, 129), (GF8, 248, 8),
        (GF16, 3, 500), (GF16, 129, 100), (GF16, 256, 64),
        (GF16, 1000, 17),
    )
    records = []
    coefficient_checks = 0
    affine_checks = 0
    aligned_records = []
    for case_index, (field, k, r) in enumerate(cases):
        systematic = tuple(range(k))
        parity = tuple(range(k, k + r))
        # Full barycentric rows are intentionally exercised for the requested
        # boundary geometries.  The independent Vandermonde construction is
        # sampled at bounded K because its cubic setup is an oracle, not a
        # candidate implementation.
        rows = barycentric_rows(field, systematic, parity)
        if any(value == 0 for row in rows for value in row):
            raise AssertionError("large-field exact generator is not dense")
        coefficient_checks += k * r
        oracle_rows = 0
        if k <= 64:
            oracle_count = r if field is GF16 and k == 3 and r == 500 else min(r, 8)
            sample_parity = parity[:oracle_count]
            oracle = vandermonde_rows(field, systematic, sample_parity)
            if tuple(rows[:len(sample_parity)]) != oracle:
                raise AssertionError("large-field Vandermonde mismatch")
            oracle_rows = len(sample_parity)

        scales = (1, 2, 3, (case_index * 37 + 5) % field.order)
        translations = (0, 1, (case_index * 53 + 17) % field.order)
        for scale in scales:
            if scale == 0:
                continue
            for translation in translations:
                # Sample complete transformed rows for large K without hiding
                # the formula behind a random coefficient comparison.
                transformed = barycentric_rows(
                    field,
                    affine_points(field, systematic, scale, translation),
                    affine_points(field, parity[:min(r, 8)], scale, translation),
                )
                if transformed != tuple(rows[:min(r, 8)]):
                    raise AssertionError("large-field affine invariance failed")
                affine_checks += min(r, 8) * k
        records.append({
            "field": field.name,
            "K": k,
            "R": r,
            "dense_coefficients": k * r,
            "vandermonde_rows": oracle_rows,
            "coordinate_digest_sha256": coordinate_digest(
                field, systematic, parity),
        })
        aligned = aligned_union_points(k, r, field.order)
        if aligned is None:
            aligned_records.append({
                "field": field.name, "K": k, "R": r,
                "available": False,
                "reason": "ceil_pow2(K)+R exceeds the field order",
            })
        else:
            aligned_rows = barycentric_rows(field, *aligned)
            changed = sum(
                left != right
                for left_row, right_row in zip(rows, aligned_rows)
                for left, right in zip(left_row, right_row)
            )
            aligned_records.append({
                "field": field.name, "K": k, "R": r,
                "available": True,
                "globally_affine": is_global_affine_image(
                    field, aligned[0] + aligned[1]),
                "changed_coefficients": changed,
                "dense_coefficients": k * r,
                "prefix_one_coefficients": sum(
                    value == 1 for row in rows for value in row),
                "aligned_one_coefficients": sum(
                    value == 1 for row in aligned_rows for value in row),
                "prefix_symbolic_fragments": (
                    dyadic_fragment_count(systematic, field.order) +
                    dyadic_fragment_count(parity, field.order)
                ),
                "aligned_symbolic_fragments": (
                    dyadic_fragment_count(aligned[0], field.order) +
                    dyadic_fragment_count(aligned[1], field.order)
                ),
            })
    return {
        "cases": len(records),
        "dense_coefficients": coefficient_checks,
        "affine_coefficients": affine_checks,
        "records": records,
        "aligned_union_records": aligned_records,
    }


def run() -> dict[str, object]:
    small = exhaustive_small_field()
    large = sampled_large_fields()
    return {
        "schema": SCHEMA,
        "status": "pass",
        "profile": {
            "identity_family": 3,
            "profile_version": 1,
            "coordinate_map_version": 1,
            "systematic_coordinates": "omega_0 .. omega_(K-1)",
            "parity_coordinates": "omega_K .. omega_(K+R-1)",
            "parent_count": "K+R",
            "padded_side": "K",
            "shortening": "none",
            "puncturing": "none",
            "future_exact_high_identity_family": 4,
        },
        "derivation": {
            "row_formula": "G[j,i] = Z(q_j) / ((q_j+x_i) * Z'(x_i))",
            "density": (
                "all evaluation points are distinct, so every numerator and "
                "denominator factor is nonzero and all K*R coefficients are nonzero"
            ),
            "affine_invariance": (
                "under phi(x)=a*x+b, a!=0, the K numerator factors contribute "
                "a^K while (q+x_i) and Z'(x_i)'s K-1 factors contribute a^K; "
                "the ratio and therefore every systematic generator byte is unchanged"
            ),
        },
        "gf4_exhaustive": small,
        "large_field": large,
        "coordinate_comparison": {
            "prefix": "selected as canonical map 1 and inherits C6 execution evidence",
            "aligned_coset": (
                "constructed as systematic 0..K-1 plus parity P..P+R-1. It "
                "is generally a non-affine dyadic union with different coefficients, "
                "and is unavailable when P+R exceeds the field; see measured symbolic "
                "fragment and coefficient-one comparisons"
            ),
            "affine_searched": (
                "exactly generator-identical; it cannot change coefficient or execution cost"
            ),
            "general_searched": (
                "can change the count of coefficient-one specializations but remains dense; "
                "it would require another coordinate-map identity and has no C6 evidence"
            ),
            "decision": "freeze the simple prefix map; do not add search tables to V1",
        },
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = run()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
