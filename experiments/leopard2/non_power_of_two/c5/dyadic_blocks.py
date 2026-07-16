#!/usr/bin/env python3
"""Scalar checkpoint for binary dyadic-block decomposition (C5).

This is deliberately an algebra and operation-count experiment, not a timing
program and not a production codec.  It studies two different constructions:

* a parent-wire-compatible decomposition of an LCH coefficient prefix into
  aligned dyadic blocks; and
* a separately identified exact prefix-evaluation Reed-Solomon code whose K
  source points are interpolated without a dyadic parent.

The two constructions must not be conflated.  The first is checked against the
full padded LCH evaluation map.  The second is checked with direct Lagrange
algebra and exhaustive GF(2^4) MDS tests.
"""

from __future__ import annotations

import argparse
import functools
import itertools
import json
from pathlib import Path
from typing import Iterable, Sequence, cast


SCHEMA_VERSION = "leopard2-c5-v1"


class BinaryField:
    """Table-backed GF(2^m) with an explicit additive coordinate basis."""

    def __init__(
        self, bits: int, polynomial: int, coordinate_basis: Sequence[int]
    ) -> None:
        if bits <= 0 or bits > 8 or len(coordinate_basis) != bits:
            raise ValueError("this scalar checkpoint supports 1 <= m <= 8")
        if not polynomial & (1 << bits):
            raise ValueError("polynomial degree does not match the field")
        self.bits = bits
        self.order = 1 << bits
        self.polynomial = polynomial
        self.coordinate_basis = tuple(coordinate_basis)

        to_polynomial = [0] * self.order
        to_coordinate = [-1] * self.order
        for coordinate in range(self.order):
            value = 0
            for bit, basis_value in enumerate(self.coordinate_basis):
                if coordinate & (1 << bit):
                    value ^= basis_value
            if value >= self.order or to_coordinate[value] >= 0:
                raise ValueError("coordinate basis is not independent")
            to_polynomial[coordinate] = value
            to_coordinate[value] = coordinate
        self.coordinate_to_polynomial = tuple(to_polynomial)
        self.polynomial_to_coordinate = tuple(to_coordinate)
        self.mul_rows = tuple(
            bytes(self._multiply_slow(left, right)
                  for right in range(self.order))
            for left in range(self.order)
        )

    def _multiply_slow(self, left: int, right: int) -> int:
        a = self.coordinate_to_polynomial[left]
        b = self.coordinate_to_polynomial[right]
        product = 0
        while b:
            if b & 1:
                product ^= a
            b >>= 1
            a <<= 1
        for bit in range(2 * self.bits - 2, self.bits - 1, -1):
            if product & (1 << bit):
                product ^= self.polynomial << (bit - self.bits)
        return self.polynomial_to_coordinate[product]

    def multiply(self, left: int, right: int) -> int:
        return self.mul_rows[left][right]

    def power(self, value: int, exponent: int) -> int:
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            exponent >>= 1
            if exponent:
                value = self.multiply(value, value)
        return result

    def inverse(self, value: int) -> int:
        if value == 0:
            raise ZeroDivisionError("zero has no inverse")
        return self.power(value, self.order - 2)


def make_gf2_4() -> BinaryField:
    return BinaryField(4, 0x13, (1, 2, 4, 8))


def make_legacy_gf8() -> BinaryField:
    # Leopard's polynomial and public Cantor-coordinate basis.
    return BinaryField(8, 0x11D, (1, 214, 152, 146, 86, 200, 88, 230))


Polynomial = list[int]
Matrix = list[list[int]]


def polynomial_multiply(
    field: BinaryField, left: Sequence[int], right: Sequence[int]
) -> Polynomial:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a == 0:
            continue
        for j, b in enumerate(right):
            if b:
                result[i + j] ^= field.multiply(a, b)
    return result


def polynomial_evaluate(
    field: BinaryField, polynomial: Sequence[int], point: int
) -> int:
    result = 0
    for coefficient in reversed(polynomial):
        result = field.multiply(result, point) ^ coefficient
    return result


def subspace_polynomial(field: BinaryField, dimension: int) -> Polynomial:
    result = [1]
    for point in range(1 << dimension):
        result = polynomial_multiply(field, result, [point, 1])
    return result


@functools.lru_cache(maxsize=None)
def normalized_lch_polynomials(
    field: BinaryField,
) -> tuple[tuple[int, ...], ...]:
    """Return all normalized LCH basis polynomials for this field.

    X_i(x) is the product of normalized subspace polynomials selected by the
    one bits of i.  Integer coordinates are additive/Cantor coordinates, so
    field addition is integer XOR.
    """
    subspaces = tuple(subspace_polynomial(field, bit)
                      for bit in range(field.bits))
    result: list[tuple[int, ...]] = []
    for index in range(field.order):
        polynomial = [1]
        normalizer = 1
        for bit, subspace in enumerate(subspaces):
            if not index & (1 << bit):
                continue
            polynomial = polynomial_multiply(field, polynomial, subspace)
            normalizer = field.multiply(
                normalizer,
                polynomial_evaluate(field, subspace, 1 << bit),
            )
        inverse = field.inverse(normalizer)
        result.append(tuple(field.multiply(value, inverse)
                            for value in polynomial))
    return tuple(result)


@functools.lru_cache(maxsize=None)
def lch_evaluation_table(
    field: BinaryField,
) -> tuple[tuple[int, ...], ...]:
    basis = normalized_lch_polynomials(field)
    return tuple(
        tuple(polynomial_evaluate(field, polynomial, point)
              for polynomial in basis)
        for point in range(field.order)
    )


def binary_prefix_blocks(length: int) -> list[tuple[int, int]]:
    """Canonical q=sum(2^a) partition into aligned additive cosets."""
    if length <= 0:
        raise ValueError("length must be positive")
    blocks: list[tuple[int, int]] = []
    offset = 0
    for exponent in range(length.bit_length() - 1, -1, -1):
        size = 1 << exponent
        if length & size:
            if offset & (size - 1):
                raise AssertionError("canonical block is not aligned")
            blocks.append((offset, size))
            offset += size
    if offset != length:
        raise AssertionError("binary partition did not cover the prefix")
    return blocks


def dyadic_interval_blocks(start: int, length: int) -> list[tuple[int, int]]:
    """Partition [start,start+length) into maximal aligned dyadic cosets."""
    if start < 0 or length <= 0:
        raise ValueError("interval must be nonempty and nonnegative")
    end = start + length
    blocks: list[tuple[int, int]] = []
    while start < end:
        size = 1 << ((end - start).bit_length() - 1)
        while start & (size - 1):
            size >>= 1
        blocks.append((start, size))
        start += size
    return blocks


def identity_matrix(size: int) -> Matrix:
    return [[int(row == column) for column in range(size)]
            for row in range(size)]


def invert_matrix(field: BinaryField, matrix: Sequence[Sequence[int]]) -> Matrix:
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("matrix must be square")
    work = [list(row) + [int(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_inverse = field.inverse(work[column][column])
        work[column] = [field.multiply(value, pivot_inverse)
                        for value in work[column]]
        pivot_row = work[column]
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            factor = work[row][column]
            target = work[row]
            work[row] = [left ^ field.multiply(factor, right)
                         for left, right in zip(target, pivot_row)]
    return [row[size:] for row in work]


def matrix_multiply(
    field: BinaryField,
    left: Sequence[Sequence[int]],
    right: Sequence[Sequence[int]],
) -> Matrix:
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("incompatible matrix dimensions")
    columns = len(right[0])
    result: Matrix = []
    for left_row in left:
        row: list[int] = []
        for column in range(columns):
            value = 0
            for inner, coefficient in enumerate(left_row):
                other = right[inner][column]
                if coefficient and other:
                    value ^= field.multiply(coefficient, other)
            row.append(value)
        result.append(row)
    return result


def matrix_rank(field: BinaryField, matrix: Sequence[Sequence[int]]) -> int:
    if not matrix:
        return 0
    rows = [list(row) for row in matrix]
    column_count = len(rows[0])
    rank = 0
    for column in range(column_count):
        pivot = next((row for row in range(rank, len(rows))
                      if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = field.inverse(rows[rank][column])
        rows[rank] = [field.multiply(value, inverse)
                      for value in rows[rank]]
        for row in range(rank + 1, len(rows)):
            factor = rows[row][column]
            if factor:
                rows[row] = [left ^ field.multiply(factor, right)
                             for left, right in zip(rows[row], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def evaluation_matrix(
    field: BinaryField, points: Iterable[int], coefficient_count: int
) -> Matrix:
    result: Matrix = []
    for point in points:
        row = [1]
        for _ in range(1, coefficient_count):
            row.append(field.multiply(row[-1], point))
        result.append(row)
    return result


def newton_evaluation_matrix(
    field: BinaryField, points: Sequence[int], roots: Sequence[int]
) -> Matrix:
    result: Matrix = []
    for point in points:
        row = [1]
        for root in roots[:-1]:
            row.append(field.multiply(row[-1], point ^ root))
        result.append(row)
    return result


def direct_lagrange_matrix(
    field: BinaryField,
    source_points: Sequence[int],
    output_points: Sequence[int],
) -> Matrix:
    """Independent direct interpolation/evaluation matrix."""
    denominators: list[int] = []
    for index, point in enumerate(source_points):
        denominator = 1
        for other_index, other in enumerate(source_points):
            if index != other_index:
                denominator = field.multiply(denominator, point ^ other)
        denominators.append(field.inverse(denominator))
    rows: Matrix = []
    for output in output_points:
        row: list[int] = []
        for index in range(len(source_points)):
            numerator = 1
            for other_index, other in enumerate(source_points):
                if index != other_index:
                    numerator = field.multiply(numerator, output ^ other)
            row.append(field.multiply(numerator, denominators[index]))
        rows.append(row)
    return rows


def block_diagonal(matrices: Sequence[Sequence[Sequence[int]]]) -> Matrix:
    total = sum(len(matrix) for matrix in matrices)
    result = [[0] * total for _ in range(total)]
    offset = 0
    for matrix in matrices:
        size = len(matrix)
        for row in range(size):
            result[offset + row][offset:offset + size] = matrix[row]
        offset += size
    return result


def matrix_cost(matrix: Sequence[Sequence[int]]) -> dict[str, int]:
    nonzero = 0
    nontrivial = 0
    xors = 0
    for row in matrix:
        row_nonzero = sum(value != 0 for value in row)
        nonzero += row_nonzero
        nontrivial += sum(value not in (0, 1) for value in row)
        xors += max(0, row_nonzero - 1)
    return {
        "entries": sum(len(row) for row in matrix),
        "nonzero": nonzero,
        "nontrivial_fixed_multiplications": nontrivial,
        "row_accumulation_xors": xors,
    }


def cross_block_cost(
    matrix: Sequence[Sequence[int]], blocks: Sequence[tuple[int, int]]
) -> dict[str, int]:
    groups = [0] * len(matrix)
    for group, (offset, size) in enumerate(blocks):
        for index in range(offset, offset + size):
            groups[index] = group
    entries = 0
    nonzero = 0
    nontrivial = 0
    for row, values in enumerate(matrix):
        for column, value in enumerate(values):
            if groups[row] == groups[column]:
                continue
            entries += 1
            nonzero += value != 0
            nontrivial += value not in (0, 1)
    return {
        "entries": entries,
        "nonzero": nonzero,
        "nontrivial_fixed_multiplications": nontrivial,
    }


def input_output_block_coupling(
    matrix: Sequence[Sequence[int]],
    input_blocks: Sequence[tuple[int, int]],
    output_blocks: Sequence[tuple[int, int]],
    output_origin: int,
) -> dict[str, int]:
    coupled = 0
    fully_dense = 0
    coefficient_entries = 0
    nonzero_coefficients = 0
    for output_offset, output_size in output_blocks:
        row_begin = output_offset - output_origin
        for input_offset, input_size in input_blocks:
            values = [
                matrix[row][column]
                for row in range(row_begin, row_begin + output_size)
                for column in range(input_offset, input_offset + input_size)
            ]
            nonzero = sum(value != 0 for value in values)
            coefficient_entries += len(values)
            nonzero_coefficients += nonzero
            coupled += nonzero != 0
            fully_dense += nonzero == len(values)
    return {
        "block_pairs": len(input_blocks) * len(output_blocks),
        "coupled_block_pairs": coupled,
        "fully_dense_block_pairs": fully_dense,
        "coefficient_entries": coefficient_entries,
        "nonzero_coefficients": nonzero_coefficients,
    }


def basis_matrices(
    field: BinaryField, q: int, basis: str
) -> tuple[Matrix, Matrix, Matrix]:
    points = list(range(q))
    blocks = binary_prefix_blocks(q)
    if basis == "lagrange":
        global_evaluation = identity_matrix(q)
        local_evaluation = identity_matrix(q)
    elif basis == "newton":
        global_evaluation = newton_evaluation_matrix(field, points, points)
        local_evaluation = block_diagonal([
            newton_evaluation_matrix(
                field, list(range(size)), list(range(size))
            )
            for _, size in blocks
        ])
    elif basis == "lch":
        table = lch_evaluation_table(field)
        global_evaluation = [list(row[:q]) for row in table[:q]]
        local_evaluation = block_diagonal([
            # This is the existing shifted-coset LCH map: coefficients remain
            # in the global X_0..X_(size-1) basis while points are in the
            # aligned coset beginning at offset.
            [list(row[:size]) for row in table[offset:offset + size]]
            for offset, size in blocks
        ])
    else:
        raise ValueError(f"unsupported basis: {basis}")
    join = matrix_multiply(field, invert_matrix(field, global_evaluation),
                           local_evaluation)
    return global_evaluation, local_evaluation, join


def output_basis_matrix(
    field: BinaryField, q: int, outputs: Sequence[int], basis: str
) -> Matrix:
    if basis == "lagrange":
        return direct_lagrange_matrix(field, list(range(q)), outputs)
    if basis == "newton":
        return newton_evaluation_matrix(field, list(outputs), list(range(q)))
    if basis == "lch":
        table = lch_evaluation_table(field)
        return [list(table[point][:q]) for point in outputs]
    raise ValueError(f"unsupported basis: {basis}")


def join_case(
    field: BinaryField, q: int, output_count: int, basis: str
) -> tuple[dict[str, object], int]:
    if q + output_count > field.order:
        raise ValueError("exact profile exceeds the field")
    blocks = binary_prefix_blocks(q)
    _, local_evaluation, join = basis_matrices(field, q, basis)
    outputs = list(range(q, q + output_count))
    output_matrix = output_basis_matrix(field, q, outputs, basis)
    effective = matrix_multiply(field, output_matrix, join)
    output_blocks = dyadic_interval_blocks(q, output_count)

    direct = direct_lagrange_matrix(field, list(range(q)), outputs)
    independent_effective = matrix_multiply(field, direct, local_evaluation)
    if effective != independent_effective:
        raise AssertionError(f"{basis} join differs from direct algebra at q={q}")
    if any(value == 0 for row in output_matrix for value in row):
        raise AssertionError(f"{basis} parity evaluation was not dense at q={q}")

    result: dict[str, object] = {
        "field_bits": field.bits,
        "q": q,
        "output_count": output_count,
        "blocks": [{"offset": offset, "size": size}
                   for offset, size in blocks],
        "basis": basis,
        "join": matrix_cost(join),
        "join_cross_block": cross_block_cost(join, blocks),
        "parity_evaluation": matrix_cost(output_matrix),
        "effective_local_to_parity": matrix_cost(effective),
        "input_output_block_coupling": input_output_block_coupling(
            effective, blocks, output_blocks, q
        ),
        "output_blocks": [
            {"offset": offset, "size": size}
            for offset, size in output_blocks
        ],
        # This is the simple immutable-input/out-of-place implementation, not
        # a lower bound: triangular or identity joins can overwrite/reuse more.
        "straightforward_out_of_place_join_slots": q,
        "largest_local_block_slots": blocks[0][1],
    }
    compared_entries = len(effective) * len(effective[0])
    return result, compared_entries


def parent_wire_case(
    field: BinaryField, q: int
) -> tuple[dict[str, object], int]:
    """Count and verify block accumulation of an LCH coefficient prefix."""
    parent = 1 << (q - 1).bit_length()
    blocks = binary_prefix_blocks(q)
    table = lch_evaluation_table(field)

    active_cosets = 0
    active_point_couplings = 0
    scale_multiplications = 0
    butterflies = 0
    comparisons = 0
    point_contributions = [0] * parent

    for offset, size in blocks:
        log_size = size.bit_length() - 1
        for coset in range(0, parent, size):
            factors = table[coset][offset:offset + 1]
            factor = factors[0]
            # X_offset is constant on a coset of V_log2(size), because every
            # set bit in offset is at least log2(size).
            for point in range(coset, coset + size):
                if table[point][offset] != factor:
                    raise AssertionError("block factor is not coset-constant")
                for local in range(size):
                    expected = table[point][offset + local]
                    # The existing shifted LCH kernel evaluates the global
                    # basis X_local at this coset; it does not translate the
                    # polynomial variable to point XOR coset.
                    actual = field.multiply(factor, table[point][local])
                    comparisons += 1
                    if actual != expected:
                        raise AssertionError(
                            "normalized LCH block product identity failed"
                        )
            if factor == 0:
                continue
            active_cosets += 1
            active_point_couplings += size
            butterflies += (size // 2) * log_size
            if factor != 1:
                scale_multiplications += size
            for point in range(coset, coset + size):
                point_contributions[point] += 1

    accumulation_xors = sum(max(0, count - 1)
                            for count in point_contributions)
    padded_butterflies = (parent // 2) * (parent.bit_length() - 1)
    result: dict[str, object] = {
        "q": q,
        "parent": parent,
        "blocks": [{"offset": offset, "size": size}
                   for offset, size in blocks],
        "block_count": len(blocks),
        "active_block_cosets": active_cosets,
        "active_point_couplings": active_point_couplings,
        "kernel_butterfly_equivalents": butterflies,
        "padded_kernel_butterfly_equivalents": padded_butterflies,
        "block_factor_fixed_multiplications": scale_multiplications,
        "cross_block_accumulation_xors": accumulation_xors,
        "kernel_logical_loads": 2 * butterflies,
        "kernel_logical_stores": 2 * butterflies,
        "largest_local_scratch_slots": blocks[0][1],
        "padded_scratch_slots": parent,
    }
    return result, comparisons


def exhaustive_gf2_4_mds(field: BinaryField) -> dict[str, int]:
    """Check every K-subset of the full 16-coordinate exact prefix code.

    Any shorter transmitted prefix is a row subset of this generator, so this
    also exhausts every valid public (K,R) geometry in GF(2^4).
    """
    subsets = 0
    generator_entries = 0
    for q in range(1, field.order):
        source = evaluation_matrix(field, range(q), q)
        source_inverse = invert_matrix(field, source)
        generator = matrix_multiply(
            field,
            evaluation_matrix(field, range(field.order), q),
            source_inverse,
        )
        if generator[:q] != identity_matrix(q):
            raise AssertionError("exact prefix generator is not systematic")
        generator_entries += field.order * q
        for selected in itertools.combinations(range(field.order), q):
            subsets += 1
            if matrix_rank(field, [generator[row] for row in selected]) != q:
                raise AssertionError(
                    f"exact prefix code is not MDS: q={q}, rows={selected}"
                )
    return {
        "systematic_generators": field.order - 1,
        "generator_entries_checked": generator_entries,
        "full_length_coordinate_subsets": subsets,
        "implied_public_k_r_geometries": (
            field.order * (field.order - 1) // 2
        ),
    }


def self_test() -> dict[str, object]:
    gf4 = make_gf2_4()
    gf8 = make_legacy_gf8()
    bases = ("lch", "newton", "lagrange")

    gf4_join_entries = 0
    gf4_join_cases = 0
    gf4_block_pairs = 0
    gf4_coupled_block_pairs = 0
    gf4_representative: list[dict[str, object]] = []
    gf4_wire_rows: list[dict[str, object]] = []
    gf4_wire_identity_comparisons = 0
    for q in range(1, gf4.order):
        output_count = gf4.order - q
        for basis in bases:
            case, entries = join_case(gf4, q, output_count, basis)
            gf4_join_cases += 1
            gf4_join_entries += entries
            coupling = cast(
                dict[str, int], case["input_output_block_coupling"]
            )
            gf4_block_pairs += coupling["block_pairs"]
            gf4_coupled_block_pairs += coupling["coupled_block_pairs"]
            if q in (3, 5, 7, 9, 13, 15):
                gf4_representative.append(case)
        wire_row, wire_comparisons = parent_wire_case(gf4, q)
        gf4_wire_rows.append(wire_row)
        gf4_wire_identity_comparisons += wire_comparisons

    # Matrix inversion is intentionally bounded in GF8.  The all-q sweep below
    # independently covers block geometry, product identities, and the complete
    # parity-output density theorem over every prefix point pair.
    gf8_join_q = (3, 5, 7, 9, 15, 17, 31, 33, 63, 65)
    gf8_join_entries = 0
    gf8_block_pairs = 0
    gf8_coupled_block_pairs = 0
    gf8_join_cases: list[dict[str, object]] = []
    for q in gf8_join_q:
        output_count = min(17, gf8.order - q)
        for basis in bases:
            case, entries = join_case(gf8, q, output_count, basis)
            gf8_join_entries += entries
            coupling = cast(
                dict[str, int], case["input_output_block_coupling"]
            )
            gf8_block_pairs += coupling["block_pairs"]
            gf8_coupled_block_pairs += coupling["coupled_block_pairs"]
            gf8_join_cases.append(case)

    if gf4_coupled_block_pairs != gf4_block_pairs:
        raise AssertionError("a GF(2^4) input/output dyadic block pair decoupled")
    if gf8_coupled_block_pairs != gf8_block_pairs:
        raise AssertionError("a selected GF8 input/output dyadic block pair decoupled")

    wire_rows: list[dict[str, object]] = []
    wire_identity_comparisons = 0
    exact_dense_coefficients = 0
    exact_dense_coefficients_checked = 0
    table = lch_evaluation_table(gf8)
    for q in range(1, gf8.order):
        row, comparisons = parent_wire_case(gf8, q)
        wire_rows.append(row)
        wire_identity_comparisons += comparisons

        # All three global bases are dense at parity points beta >= q:
        # Lagrange factors beta+alpha, Newton factors beta+alpha, and every LCH
        # subspace factor have roots strictly below q.  Check LCH explicitly;
        # direct join cases above check all three bases algebraically.
        coefficient_count = q * (gf8.order - q)
        exact_dense_coefficients += coefficient_count * len(bases)
        for point in range(q, gf8.order):
            for coefficient in range(q):
                exact_dense_coefficients_checked += 1
                if table[point][coefficient] == 0:
                    raise AssertionError("GF8 LCH parity map was not dense")

    mds = exhaustive_gf2_4_mds(gf4)
    selected_wire_q = {
        3, 5, 7, 9, 15, 17, 31, 33, 63, 65, 127, 129, 191, 223,
        239, 247, 249, 255,
    }
    representative_wire = [row for row in wire_rows if row["q"] in selected_wire_q]

    return {
        "schema": SCHEMA_VERSION,
        "scope": {
            "timing_performed": False,
            "production_code": False,
            "wire_compatible_candidate_changes_code_identity": False,
            "exact_prefix_profile_requires_new_code_identity": True,
        },
        "coordinate_systems": {
            "gf2_4": {
                "polynomial": "0x13",
                "coordinate_basis": [1, 2, 4, 8],
            },
            "legacy_gf8": {
                "polynomial": "0x11d",
                "cantor_coordinate_basis": [1, 214, 152, 146, 86, 200, 88, 230],
            },
        },
        "gf2_4": {
            "prefix_geometries": gf4.order - 1,
            "basis_join_cases": gf4_join_cases,
            "effective_matrix_entries_compared": gf4_join_entries,
            "input_output_block_pairs": gf4_block_pairs,
            "coupled_input_output_block_pairs": gf4_coupled_block_pairs,
            "all_prefix_parent_rows": gf4_wire_rows,
            "parent_block_identity_entries_compared": (
                gf4_wire_identity_comparisons
            ),
            "mds": mds,
            "representative_join_costs": gf4_representative,
        },
        "legacy_gf8": {
            "all_prefix_geometries": gf8.order - 1,
            "all_prefix_parent_rows": wire_rows,
            "parent_block_identity_entries_compared": wire_identity_comparisons,
            "join_q_values": list(gf8_join_q),
            "basis_join_cases": len(gf8_join_cases),
            "effective_matrix_entries_compared": gf8_join_entries,
            "input_output_block_pairs": gf8_block_pairs,
            "coupled_input_output_block_pairs": gf8_coupled_block_pairs,
            "join_costs": gf8_join_cases,
            "global_parity_matrix_coefficients_proved_nonzero_all_bases": (
                exact_dense_coefficients
            ),
            "lch_parity_matrix_coefficients_checked_nonzero": (
                exact_dense_coefficients_checked
            ),
            "representative_parent_costs": representative_wire,
        },
        "accounting_contract": {
            "matrix_fixed_multiplications": (
                "nonzero coefficients other than one; setup excluded"
            ),
            "matrix_xors": "row nonzero count minus one; setup excluded",
            "butterfly_equivalent": (
                "one complete radix-2 LCH butterfly; zero/one skew "
                "specialization is not estimated"
            ),
            "kernel_load_store": (
                "two logical operands/results per butterfly; block packing, "
                "external scatter, caches, and Python objects excluded"
            ),
            "scratch_slots": (
                "field-symbol/shard slots, excluding caller-owned input and output"
            ),
        },
        "disposition": {
            "checkpoint_scope": "bounded_scalar_nonpromotion_result",
            "explicit_same_basis_join_plus_parent": (
                "do_not_promote_without_a_new_factorization_and_measurement"
            ),
            "broader_binary_block_exact_factorization": (
                "inconclusive_parent_c5_remains_open"
            ),
            "reason": (
                "the explicit local-to-global basis join adds cross-block work before "
                "parent evaluation, and every sampled input/output dyadic-block pair "
                "remains coupled; density and coupling do not rule out an unstudied "
                "factored or pruned-output schedule, so this checkpoint rejects only "
                "promotion of the naive explicit join-plus-parent route"
            ),
            "parent_wire_block_accumulation": (
                "correct_but_do_not_promote_separately_from_c1_c2"
            ),
            "parent_c5_bead": (
                "leave_open_for_representative_gf16_and_measured_byte_batch_reuse_gates"
            ),
            "next_gate": (
                "measure any surviving optimized factorization in representative "
                "GF16 and byte/batch/reuse cells; only promote at >=10% end-to-end"
            ),
        },
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="deterministic JSON evidence path",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = self_test()
    write_json(args.output, result)
    print(json.dumps({
        "schema": result["schema"],
        "gf2_4_mds_subsets": result["gf2_4"]["mds"]["full_length_coordinate_subsets"],
        "gf8_prefix_geometries": result["legacy_gf8"]["all_prefix_geometries"],
        "output": str(args.output),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
