#!/usr/bin/env python3
"""Independent scalar algebra and cost evidence for Leopard2 experiment C8.

This is a default-off research program.  It defines the candidate coordinate
map explicitly and never changes the public Leopard2 profile or dispatcher.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import dataclasses
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
import os
import random
import struct
import sys
from pathlib import Path
from typing import Iterable, Sequence


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
C3B_PATH = ROOT / "experiments/leopard2/non_power_of_two/c3b/fast_inverse.py"
PROFILE_NAME = "exact_high_prefix_v1_candidate"
COORDINATE_MAP_VERSION = 1


def load_c3b():
    specification = importlib.util.spec_from_file_location("leo2_c3b", C3B_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load the C3b field oracle")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


c3b = load_c3b()


@dataclasses.dataclass
class Counter:
    multiplications: int = 0
    inversions: int = 0
    xors: int = 0
    loads: int = 0
    stores: int = 0

    def as_dict(self) -> dict[str, int]:
        return dataclasses.asdict(self)


def checked_counts(field, k: int, r: int) -> None:
    if k <= 0 or r <= 0 or k + r > field.order:
        raise ValueError("invalid exact-high counts")


def multiply(field, left: int, right: int, counter: Counter | None = None) -> int:
    if counter is not None:
        counter.multiplications += 1
    return field.multiply(left, right)


def inverse(field, value: int, counter: Counter | None = None) -> int:
    if counter is not None:
        counter.inversions += 1
    return field.inverse(value)


def matrix_apply(
    field, matrix: Sequence[Sequence[int]], values: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    if any(len(row) != len(values) for row in matrix):
        raise ValueError("matrix/vector dimensions differ")
    output = []
    for row in matrix:
        value = 0
        first = True
        for coefficient, source in zip(row, values):
            if counter is not None:
                counter.loads += 2
            if coefficient and source:
                term = multiply(field, coefficient, source, counter)
                if not first:
                    if counter is not None:
                        counter.xors += 1
                    value ^= term
                else:
                    value = term
                    first = False
        if counter is not None:
            counter.stores += 1
        output.append(value)
    return output


def matrix_multiply(field, left, right, counter: Counter | None = None):
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("invalid matrix dimensions")
    columns = len(right[0])
    if any(len(row) != columns for row in right):
        raise ValueError("ragged right matrix")
    output = []
    for left_row in left:
        row = []
        for column in range(columns):
            value = 0
            first = True
            for index, coefficient in enumerate(left_row):
                other = right[index][column]
                if coefficient and other:
                    product = multiply(field, coefficient, other, counter)
                    if not first:
                        if counter is not None:
                            counter.xors += 1
                        value ^= product
                    else:
                        value = product
                        first = False
            row.append(value)
        output.append(row)
    return output


def matrix_inverse(field, matrix, counter: Counter | None = None):
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
        scale = inverse(field, work[column][column], counter)
        for index in range(2 * size):
            work[column][index] = multiply(field, work[column][index], scale,
                                           counter)
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            factor = work[row][column]
            for index in range(2 * size):
                work[row][index] ^= multiply(
                    field, factor, work[column][index], counter)
                if counter is not None:
                    counter.xors += 1
    return [row[size:] for row in work]


def matrix_rank(field, matrix) -> int:
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


def lagrange_rows(field, k: int, r: int,
                  counter: Counter | None = None) -> tuple[tuple[int, ...], ...]:
    """Generator rows for parity points 0:R and message points R:R+K."""
    checked_counts(field, k, r)
    message_points = tuple(range(r, r + k))
    weights = []
    for point in message_points:
        denominator = 1
        for other in message_points:
            if other != point:
                denominator = multiply(field, denominator, point ^ other,
                                       counter)
        weights.append(inverse(field, denominator, counter))
    rows = []
    for parity_point in range(r):
        vanishing = 1
        for point in message_points:
            vanishing = multiply(field, vanishing, parity_point ^ point,
                                 counter)
        row = []
        for index, point in enumerate(message_points):
            coefficient = multiply(
                field,
                multiply(field, vanishing,
                         inverse(field, parity_point ^ point, counter), counter),
                weights[index], counter)
            row.append(coefficient)
        rows.append(tuple(row))
    return tuple(rows)


def direct_lagrange_evaluate(field, values: Sequence[int], nodes: Sequence[int],
                             point: int) -> int:
    if len(values) != len(nodes):
        raise ValueError("value/node size mismatch")
    output = 0
    for source_index, source in enumerate(values):
        numerator = 1
        denominator = 1
        for other_index, other in enumerate(nodes):
            if other_index == source_index:
                continue
            numerator = field.multiply(numerator, point ^ other)
            denominator = field.multiply(
                denominator, nodes[source_index] ^ other)
        output ^= field.multiply(
            source, field.multiply(numerator, field.inverse(denominator)))
    return output


def barycentric_weights(field, nodes: Sequence[int]) -> tuple[int, ...]:
    weights = []
    for node in nodes:
        denominator = 1
        for other in nodes:
            if other != node:
                denominator = field.multiply(denominator, node ^ other)
        weights.append(field.inverse(denominator))
    return tuple(weights)


def barycentric_evaluate(field, values: Sequence[int], nodes: Sequence[int],
                         weights: Sequence[int], point: int) -> int:
    if len(values) != len(nodes) or len(weights) != len(nodes):
        raise ValueError("barycentric dimensions differ")
    vanishing = 1
    for node in nodes:
        vanishing = field.multiply(vanishing, point ^ node)
    output = 0
    for source, node, weight in zip(values, nodes, weights):
        coefficient = field.multiply(
            field.multiply(vanishing, field.inverse(point ^ node)), weight)
        output ^= field.multiply(source, coefficient)
    return output


def dual_constraint_matrices(field, k: int, r: int,
                             counter: Counter | None = None):
    """Build the independent GRS dual constraints A*p + B*m = 0.

    For prefix points x_i=omega_i, dual weight u_i is the inverse derivative of
    the full transmitted-set vanishing polynomial.  Rows u_i*x_i^j,
    j=0..R-1, annihilate every evaluation vector of degree less than K.
    """
    checked_counts(field, k, r)
    n = k + r
    weights = []
    for point in range(n):
        derivative = 1
        for other in range(n):
            if other != point:
                derivative = multiply(field, derivative, point ^ other, counter)
        weights.append(inverse(field, derivative, counter))
    rows = []
    for exponent in range(r):
        row = []
        for point in range(n):
            power = 1
            for _ in range(exponent):
                power = multiply(field, power, point, counter)
            row.append(multiply(field, weights[point], power, counter))
        rows.append(row)
    a = [row[:r] for row in rows]
    b = [row[r:] for row in rows]
    return a, b


def dense_solve_rows(field, k: int, r: int,
                     counter: Counter | None = None):
    a, b = dual_constraint_matrices(field, k, r, counter)
    return tuple(tuple(value for value in row) for row in
                 matrix_multiply(field, matrix_inverse(field, a, counter), b,
                                 counter))


def newton_coefficients(field, nodes: Sequence[int], values: Sequence[int],
                        counter: Counter | None = None) -> list[int]:
    if len(nodes) != len(values):
        raise ValueError("node/value size mismatch")
    coefficients = list(values)
    for order in range(1, len(nodes)):
        for index in range(len(nodes) - 1, order - 1, -1):
            coefficients[index] ^= coefficients[index - 1]
            if counter is not None:
                counter.xors += 1
            coefficients[index] = multiply(
                field, coefficients[index],
                inverse(field, nodes[index] ^ nodes[index - order], counter),
                counter)
    return coefficients


def newton_evaluate(field, nodes: Sequence[int], coefficients: Sequence[int],
                    point: int, counter: Counter | None = None) -> int:
    value = coefficients[-1]
    for index in range(len(coefficients) - 2, -1, -1):
        value = multiply(field, value, point ^ nodes[index], counter)
        value ^= coefficients[index]
        if counter is not None:
            counter.xors += 1
    return value


def newton_encode(field, k: int, r: int, values: Sequence[int],
                  counter: Counter | None = None) -> list[int]:
    nodes = tuple(range(r, r + k))
    coefficients = newton_coefficients(field, nodes, values, counter)
    return [newton_evaluate(field, nodes, coefficients, point, counter)
            for point in range(r)]


@dataclasses.dataclass(frozen=True)
class SchurPlan:
    split: int
    a11_inverse: tuple[tuple[int, ...], ...]
    a12: tuple[tuple[int, ...], ...]
    a21: tuple[tuple[int, ...], ...]
    schur_inverse: tuple[tuple[int, ...], ...]
    b: tuple[tuple[int, ...], ...]


def schur_plan(field, k: int, r: int,
               counter: Counter | None = None) -> SchurPlan:
    if r < 3 or r & (r - 1) == 0:
        raise ValueError("Schur plan requires a non-power-of-two R >= 3")
    a, b = dual_constraint_matrices(field, k, r, counter)
    split = 1 << (r.bit_length() - 1)
    a11 = [row[:split] for row in a[:split]]
    a12 = [row[split:] for row in a[:split]]
    a21 = [row[:split] for row in a[split:]]
    a22 = [row[split:] for row in a[split:]]
    a11_inverse = matrix_inverse(field, a11, counter)
    lower_left = matrix_multiply(field, a21, a11_inverse, counter)
    correction = matrix_multiply(field, lower_left, a12, counter)
    schur = [[left ^ right for left, right in zip(row_left, row_right)]
             for row_left, row_right in zip(a22, correction)]
    if counter is not None:
        counter.xors += sum(len(row) for row in schur)
    return SchurPlan(
        split,
        tuple(tuple(row) for row in a11_inverse),
        tuple(tuple(row) for row in a12),
        tuple(tuple(row) for row in a21),
        tuple(tuple(row) for row in matrix_inverse(field, schur, counter)),
        tuple(tuple(row) for row in b),
    )


def schur_encode(field, plan: SchurPlan, values: Sequence[int],
                 counter: Counter | None = None) -> list[int]:
    right = matrix_apply(field, plan.b, values, counter)
    split = plan.split
    y1 = matrix_apply(field, plan.a11_inverse, right[:split], counter)
    a21_y1 = matrix_apply(field, plan.a21, y1, counter)
    right2 = [left ^ other for left, other in zip(right[split:], a21_y1)]
    if counter is not None:
        counter.xors += len(right2)
    x2 = matrix_apply(field, plan.schur_inverse, right2, counter)
    a12_x2 = matrix_apply(field, plan.a12, x2, counter)
    corrected1 = [left ^ other for left, other in zip(right[:split], a12_x2)]
    if counter is not None:
        counter.xors += len(corrected1)
    x1 = matrix_apply(field, plan.a11_inverse, corrected1, counter)
    return x1 + x2


def ceil_pow2(value: int) -> int:
    return 1 << (value - 1).bit_length()


def legacy_high_rows(field, k: int, r: int):
    """Direct oracle rows for LEGACY_HIGH_V1, including shortened zeros."""
    t = ceil_pow2(r)
    n = ceil_pow2(k + t)
    if n > field.order:
        raise ValueError("legacy parent exceeds field")
    systematic_points = tuple(range(t, n))
    weights = []
    for point in systematic_points:
        denominator = 1
        for other in systematic_points:
            if other != point:
                denominator = field.multiply(denominator, point ^ other)
        weights.append(field.inverse(denominator))
    rows = []
    for parity in range(r):
        vanishing = 1
        for point in systematic_points:
            vanishing = field.multiply(vanishing, parity ^ point)
        rows.append(tuple(
            field.multiply(
                field.multiply(vanishing, field.inverse(parity ^ (t + index))),
                weights[index])
            for index in range(k)
        ))
    return tuple(rows)


def legacy_high_coefficient(field, k: int, r: int,
                            parity: int, systematic: int) -> int:
    t = ceil_pow2(r)
    n = ceil_pow2(k + t)
    if n > field.order or systematic >= k or parity >= r:
        raise ValueError("invalid legacy coefficient")
    source = t + systematic
    numerator = 1
    denominator = 1
    for point in range(t, n):
        if point != source:
            numerator = field.multiply(numerator, parity ^ point)
            denominator = field.multiply(denominator, source ^ point)
    return field.multiply(numerator, field.inverse(denominator))


def exact_identity_payload(field_name: str, k: int, r: int) -> bytes:
    encoded = bytearray(b"LEO2-EXACT-HIGH-PREFIX-V1-CANDIDATE\0")
    encoded += struct.pack(">BIII", 8 if field_name == "gf8" else 16,
                           k, r, COORDINATE_MAP_VERSION)
    encoded += struct.pack(">I", r)
    for point in range(r):
        encoded += struct.pack(">BI", 1, point)
    encoded += struct.pack(">I", k)
    for point in range(r, r + k):
        encoded += struct.pack(">BI", 2, point)
    return bytes(encoded)


def exact_identity_digest(field_name: str, k: int, r: int) -> str:
    return hashlib.sha256(exact_identity_payload(field_name, k, r)).hexdigest()


def systematic_generator(rows, k: int):
    return tuple(rows) + tuple(tuple(int(row == column) for column in range(k))
                               for row in range(k))


def stable_values(field, k: int, r: int, salt: int) -> list[int]:
    rng = random.Random((field.bits << 56) ^ (k << 32) ^ (r << 16) ^ salt)
    return [rng.randrange(field.order) for _ in range(k)]


def exhaustive_gf4() -> dict[str, int]:
    field = c3b.field_named("gf4")
    geometries = subsets = solver_vectors = parity_symbols = 0
    full_parent_identity = nonidentity_witnesses = 0
    for n in range(2, 17):
        for k in range(1, n):
            r = n - k
            rows = lagrange_rows(field, k, r)
            dense = dense_solve_rows(field, k, r)
            if rows != dense:
                raise AssertionError(f"dual solve mismatch GF4 K={k} R={r}")
            generator = systematic_generator(rows, k)
            for selected in itertools.combinations(range(n), k):
                subsets += 1
                if matrix_rank(field, [generator[index] for index in selected]) != k:
                    raise AssertionError(
                        f"non-MDS exact-high GF4 K={k} R={r} rows={selected}")
            geometries += 1
            vectors: Iterable[Sequence[int]]
            if k <= 3:
                vectors = itertools.product(range(16), repeat=k)
            else:
                vectors = ([int(index == basis) for index in range(k)]
                           for basis in range(k))
            plan = schur_plan(field, k, r) if r >= 3 and r & (r - 1) else None
            for values_tuple in vectors:
                values = list(values_tuple)
                expected = matrix_apply(field, rows, values)
                observed = [direct_lagrange_evaluate(
                    field, values, tuple(range(r, r + k)), point)
                    for point in range(r)]
                if expected != observed or newton_encode(field, k, r, values) != expected:
                    raise AssertionError("independent GF4 encoder mismatch")
                if plan is not None and schur_encode(field, plan, values) != expected:
                    raise AssertionError("GF4 Schur encoder mismatch")
                solver_vectors += 1
                parity_symbols += r
            legacy_possible = ceil_pow2(k + ceil_pow2(r)) <= field.order
            if legacy_possible:
                legacy = legacy_high_rows(field, k, r)
                should_match = (r & (r - 1) == 0 and k + r == ceil_pow2(k + r))
                if should_match:
                    if legacy != rows:
                        raise AssertionError("full-parent legacy identity failed")
                    full_parent_identity += 1
                elif legacy != rows:
                    nonidentity_witnesses += 1
                else:
                    raise AssertionError("unexpected legacy identity outside full parent")
    return {
        "public_geometries": geometries,
        "k_coordinate_subsets": subsets,
        "solver_vectors": solver_vectors,
        "parity_symbols": parity_symbols,
        "full_parent_wire_identity_geometries": full_parent_identity,
        "nonidentity_wire_witness_geometries": nonidentity_witnesses,
    }


def differential_specs() -> tuple[tuple[str, int, int], ...]:
    return (
        ("gf8", 1, 1), ("gf8", 3, 1), ("gf8", 5, 3),
        ("gf8", 12, 4), ("gf8", 24, 8), ("gf8", 31, 3),
        ("gf8", 48, 16), ("gf8", 60, 4), ("gf8", 63, 1),
        ("gf8", 96, 32), ("gf8", 120, 7), ("gf8", 120, 8),
        ("gf8", 120, 9), ("gf8", 192, 64), ("gf8", 224, 31),
        ("gf16", 255, 1), ("gf16", 500, 3), ("gf16", 1000, 17),
        ("gf16", 1000, 33),
    )


def differential_job(specification: tuple[str, int, int]) -> dict[str, object]:
    field_name, k, r = specification
    field = c3b.field_named(field_name)
    rows = lagrange_rows(field, k, r)
    coefficient_checks = k * r
    dense_checked = r <= 33 and k + r <= 1100
    if dense_checked and rows != dense_solve_rows(field, k, r):
        raise AssertionError("dual constraint rows differ")
    schur_checked = dense_checked and r >= 3 and r & (r - 1) != 0
    plan = schur_plan(field, k, r) if schur_checked else None
    digest = hashlib.sha256()
    vector_count = 0
    nodes = tuple(range(r, r + k))
    independent_weights = barycentric_weights(field, nodes)
    for salt in range(1 if k > 1000 else 4):
        values = stable_values(field, k, r, salt)
        expected = matrix_apply(field, rows, values)
        direct = [barycentric_evaluate(
            field, values, nodes, independent_weights, point)
                  for point in range(r)]
        if expected != direct or newton_encode(field, k, r, values) != expected:
            raise AssertionError("GF8/GF16 direct/Newton mismatch")
        if plan is not None and schur_encode(field, plan, values) != expected:
            raise AssertionError("GF8/GF16 Schur mismatch")
        width = 1 if field.bits == 8 else 2
        for value in expected:
            digest.update(value.to_bytes(width, "little"))
        vector_count += 1
    should_match = r & (r - 1) == 0 and k + r == ceil_pow2(k + r)
    if should_match:
        if legacy_high_rows(field, k, r) != rows:
            raise AssertionError("full-parent legacy wire identity failed")
    else:
        witnesses = ((0, 0), (r // 2, k // 2), (r - 1, k - 1))
        if not any(rows[parity][systematic] != legacy_high_coefficient(
                field, k, r, parity, systematic)
                   for parity, systematic in witnesses):
            raise AssertionError("missing nonidentity legacy wire witness")
    return {
        "field": field_name,
        "K": k,
        "R": r,
        "identity_sha256": exact_identity_digest(field_name, k, r),
        "coefficients": coefficient_checks,
        "vectors": vector_count,
        "dense_constraint_checked": dense_checked,
        "schur_checked": schur_checked,
        "legacy_wire_identity": should_match,
        "parity_digest_sha256": digest.hexdigest(),
    }


def execution_counts(k: int, r: int, method: str) -> dict[str, int | str | bool]:
    t = ceil_pow2(r)
    log_t = t.bit_length() - 1
    if method == "precomputed_dense_r_solve":
        return {"implemented": True, "multiplications": k * r,
                "xors": r * max(0, k - 1), "temporary_symbols": 0,
                "logical_loads": 2 * k * r,
                "logical_stores": k * r, "irregular": 0}
    if method == "factorized_r_solve":
        return {"implemented": True,
                "multiplications": k * r + r * r,
                "xors": r * max(0, k - 1) + r * max(0, r - 1),
                "temporary_symbols": 2 * r, "logical_loads": 2 * (k * r + r * r),
                "logical_stores": k * r + r * r, "irregular": r * r}
    if method == "newton":
        work = k * (k - 1) // 2 + k * r
        return {"implemented": True, "multiplications": work,
                "xors": work, "temporary_symbols": k,
                "logical_loads": 2 * work, "logical_stores": work,
                "irregular": k + r}
    if method == "dyadic_schur":
        if r < 3 or r & (r - 1) == 0:
            return {"implemented": False, "reason": "R is a power of two or tiny"}
        b = 1 << (r.bit_length() - 1)
        d = r - b
        work = k * r + 2 * b * b + 2 * b * d + d * d
        return {"implemented": True, "multiplications": work,
                "xors": work - r, "temporary_symbols": 2 * r,
                "logical_loads": 2 * work, "logical_stores": work,
                "irregular": b * d + d * d}
    if method == "padded_t_lch":
        blocks = (k + t - 1) // t
        butterflies = (blocks + 1) * (t // 2) * log_t
        return {"implemented": True, "multiplications": butterflies,
                "xors": 2 * butterflies, "temporary_symbols": 2 * t,
                "logical_loads": 2 * butterflies,
                "logical_stores": 2 * butterflies, "irregular": t - r}
    if method == "partial_transposed_lch_lower_bound":
        # This deliberately optimistic bound omits an executable irregular
        # transposed-LCH schedule.  It may prioritize later work, never promote.
        work = (k + r) * max(1, math.ceil(math.log2(max(2, r))))
        return {"implemented": False, "lower_bound": True,
                "multiplications": work, "xors": 2 * work,
                "temporary_symbols": r, "logical_loads": 2 * work,
                "logical_stores": 2 * work, "irregular": t - r,
                "reason": "no proven executable transposed exact-R LCH schedule"}
    if method == "tang_han_epsilon_control":
        p = ceil_pow2(k)
        work = (p // 2) * max(1, p.bit_length() - 1) + k * r
        return {"implemented": False, "control_profile_only": True,
                "multiplications": work, "xors": 2 * work,
                "temporary_symbols": p, "logical_loads": 2 * work,
                "logical_stores": 2 * work, "irregular": p - k,
                "reason": "published prefix IFFT does not interpolate the C8 suffix coordinate set"}
    raise ValueError("unknown method")


METHODS = (
    "precomputed_dense_r_solve", "factorized_r_solve", "newton",
    "dyadic_schur", "padded_t_lch", "partial_transposed_lch_lower_bound",
    "tang_han_epsilon_control",
)


def model_row(pair: tuple[int, int]) -> list[dict[str, object]]:
    k, r = pair
    rows = []
    for method in METHODS:
        counts = execution_counts(k, r, method)
        numeric = {name: value for name, value in counts.items()
                   if isinstance(value, int) and not isinstance(value, bool)}
        score = (3 * numeric.get("multiplications", 0) +
                 numeric.get("xors", 0) + numeric.get("logical_loads", 0) +
                 2 * numeric.get("logical_stores", 0) +
                 4 * numeric.get("irregular", 0))
        if method == "padded_t_lch":
            wire_scope = "legacy_high_v1"
        elif method == "tang_han_epsilon_control":
            wire_scope = "different_prefix_coordinate_control"
        else:
            wire_scope = PROFILE_NAME
        rows.append({"K": k, "R": r, "method": method,
                     "wire_scope": wire_scope,
                     **counts, "score": score})
    return rows


def canonical_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n",
                    encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def source_sha256() -> str:
    return sha256(Path(__file__))


def analyze(workers: int, output_dir: Path) -> dict[str, object]:
    if workers <= 0:
        raise ValueError("workers must be positive")
    gf4 = exhaustive_gf4()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        differential = list(executor.map(differential_job, differential_specs()))
    pairs = [(k, r) for k in range(1, 256) for r in range(1, 257 - k)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        modeled = list(executor.map(model_row, pairs, chunksize=64))
    flat = [row for group in modeled for row in group]
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = output_dir / "gf8_model.csv.gz"
    columns = sorted({key for row in flat for key in row})
    with matrix_path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
            import io
            wrapper = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
            writer = csv.DictWriter(wrapper, fieldnames=columns, extrasaction="raise")
            writer.writeheader()
            writer.writerows(flat)
            wrapper.flush()
    winner_counts: dict[str, int] = {}
    gains_over_padded: dict[str, list[float]] = {method: [] for method in METHODS}
    for group in modeled:
        usable = [row for row in group if row.get("implemented")]
        winner = min(usable, key=lambda row: (int(row["score"]), str(row["method"])))
        winner_counts[str(winner["method"])] = winner_counts.get(str(winner["method"]), 0) + 1
        padded = next(row for row in group if row["method"] == "padded_t_lch")
        for row in group:
            if int(row["score"]):
                gains_over_padded[str(row["method"])].append(
                    int(padded["score"]) / int(row["score"]))
    summary = {
        "schema": "leopard2-c8-algebra-v1",
        "source_sha256": source_sha256(),
        "profile": {
            "name": PROFILE_NAME,
            "status": "unserialized_default_off_candidate",
            "coordinate_map_version": COORDINATE_MAP_VERSION,
            "wire_order": ["parity evaluations 0..R-1",
                           "systematic evaluations R..R+K-1"],
            "polynomial_degree_bound": "degree < K",
        },
        "workers_requested": workers,
        "workers_used": min(workers, len(differential_specs())),
        "gf4": gf4,
        "differential": differential,
        "gf8_model": {
            "public_pairs": len(pairs), "method_rows": len(flat),
            "winner_counts": winner_counts,
            "gain_over_padded": {
                method: {
                    "maximum": max(values),
                    "mean": sum(values) / len(values),
                    "cells_at_least_1_10": sum(value >= 1.10 for value in values),
                } for method, values in gains_over_padded.items()
            },
        },
        "method_status": {
            "precomputed_dense_r_solve": "executable scalar algebra and C++ candidate",
            "factorized_r_solve": "executable scalar algebra",
            "newton": "executable independent scalar control",
            "dyadic_schur": "executable scalar algebra for non-power R",
            "padded_t_lch": "production legacy baseline, different wire except identity gate",
            "partial_transposed_lch_lower_bound": "non-executable optimistic model",
            "tang_han_epsilon_control": "not applicable to suffix coordinate set without a new derivation",
        },
        "disposition": "pending executable benchmark",
    }
    canonical_json(output_dir / "algebra.json", summary)
    manifest = {
        "schema": "leopard2-c8-manifest-v1",
        "files": {
            "algebra.json": sha256(output_dir / "algebra.json"),
            "gf8_model.csv.gz": sha256(matrix_path),
        },
    }
    canonical_json(output_dir / "algebra-manifest.json", manifest)
    return summary


def verify(output_dir: Path) -> None:
    manifest_path = output_dir / "algebra-manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "leopard2-c8-manifest-v1":
        raise ValueError("wrong C8 manifest schema")
    for name, expected in manifest["files"].items():
        observed = sha256(output_dir / name)
        if observed != expected:
            raise ValueError(f"hash mismatch for {name}")
    result = json.loads((output_dir / "algebra.json").read_text(encoding="utf-8"))
    if result.get("schema") != "leopard2-c8-algebra-v1":
        raise ValueError("wrong C8 algebra schema")
    if result.get("source_sha256") != source_sha256():
        raise ValueError("stale C8 algebra source hash")
    gf4 = result.get("gf4", {})
    if gf4.get("public_geometries") != 120:
        raise ValueError("incomplete GF4 geometry coverage")
    if gf4.get("k_coordinate_subsets") != 131038:
        raise ValueError("incomplete GF4 MDS coverage")
    model = result.get("gf8_model", {})
    if model.get("public_pairs") != 32640 or model.get("method_rows") != 228480:
        raise ValueError("incomplete GF8 model matrix")
    differential = result.get("differential", [])
    if len(differential) != len(differential_specs()):
        raise ValueError("incomplete differential matrix")
    for row in differential:
        expected = (int(row["R"]) & (int(row["R"]) - 1) == 0 and
                    int(row["K"]) + int(row["R"]) ==
                    ceil_pow2(int(row["K"]) + int(row["R"])))
        if row.get("legacy_wire_identity") is not expected:
            raise ValueError("forged wire identity gate")


def self_test() -> None:
    field = c3b.field_named("gf4")
    for k, r in ((1, 1), (3, 1), (3, 2), (5, 3), (8, 7)):
        rows = lagrange_rows(field, k, r)
        if rows != dense_solve_rows(field, k, r):
            raise AssertionError("self-test dense solve mismatch")
        values = stable_values(field, k, r, 9)
        expected = matrix_apply(field, rows, values)
        if newton_encode(field, k, r, values) != expected:
            raise AssertionError("self-test Newton mismatch")
        if r >= 3 and r & (r - 1):
            if schur_encode(field, schur_plan(field, k, r), values) != expected:
                raise AssertionError("self-test Schur mismatch")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser("self-test")
    analyze_parser = commands.add_parser("analyze")
    analyze_parser.add_argument("--workers", type=int, default=min(128, os.cpu_count() or 1))
    analyze_parser.add_argument("--output-dir", type=Path, default=HERE / "results")
    verify_parser = commands.add_parser("verify")
    verify_parser.add_argument("--output-dir", type=Path, default=HERE / "results")
    return result


def main() -> int:
    arguments = parser().parse_args()
    if arguments.command == "self-test":
        self_test()
    elif arguments.command == "analyze":
        analyze(arguments.workers, arguments.output_dir)
    elif arguments.command == "verify":
        verify(arguments.output_dir)
    else:
        raise AssertionError("unhandled command")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
