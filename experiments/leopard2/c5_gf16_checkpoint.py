#!/usr/bin/env python3
"""Deterministic GF16 and memory checkpoint for Leopard2 experiment C5.

This program evaluates two deliberately separate constructions:

* ``parent_wire_block_accumulation`` computes the existing padded LCH map by
  splitting an active coefficient prefix into aligned dyadic blocks.  It does
  not define a new code or change parity bytes.
* ``exact_prefix_join`` interpolates exactly K prefix points and therefore is a
  new, experimental code identity unless generator equivalence is proved.

The program is a scalar algebra and traffic simulator.  It does not report
Python wall-clock time as codec performance evidence.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
from pathlib import Path
from typing import Iterable, Sequence


SCHEMA = "leopard2-c5-gf16-checkpoint/v1"
GF16_POLYNOMIAL = 0x1002D
GF16_CANTOR_BASIS = (
    0x0001, 0xACCA, 0x3C0E, 0x163E,
    0xC582, 0xED2E, 0x914C, 0x4012,
    0x6C98, 0x10D8, 0x6A72, 0xB900,
    0xFDB8, 0xFB34, 0xFF38, 0x991E,
)
PARENT_Q = (
    257, 511, 513, 1000, 1023, 1025, 2047, 2049,
    4095, 4096, 4097, 8191, 8193, 16383, 16385,
    32767, 32769, 65535,
)
JOIN_Q = (17, 33, 63, 65, 129, 191, 255, 257)
SHARD_BYTES = (64, 1024, 64 * 1024, 1024 * 1024, 16 * 1024 * 1024)
BATCHES = (1, 8, 64, 1024)
REUSES = (1, 8, 64, 1024)


class BinaryField:
    """GF(2^m) with additive coordinates and O(1) log multiplication."""

    def __init__(
        self, bits: int, polynomial: int, coordinate_basis: Sequence[int]
    ) -> None:
        self.bits = bits
        self.order = 1 << bits
        self.modulus = self.order - 1
        self.polynomial = polynomial
        if len(coordinate_basis) != bits or not polynomial & self.order:
            raise ValueError("field definition has the wrong dimension")

        coordinate_to_polynomial = [0] * self.order
        polynomial_to_coordinate = [-1] * self.order
        for coordinate in range(self.order):
            value = 0
            for bit, basis_value in enumerate(coordinate_basis):
                if coordinate & (1 << bit):
                    value ^= basis_value
            if value >= self.order or polynomial_to_coordinate[value] >= 0:
                raise ValueError("coordinate basis is not independent")
            coordinate_to_polynomial[coordinate] = value
            polynomial_to_coordinate[value] = coordinate

        exponent = [0] * self.modulus
        logarithm = [-1] * self.order
        state = 1
        for power in range(self.modulus):
            if logarithm[state] >= 0:
                raise ValueError("field polynomial does not provide a primitive x")
            exponent[power] = state
            logarithm[state] = power
            state <<= 1
            if state & self.order:
                state ^= polynomial
        if state != 1 or any(value < 0 for value in logarithm[1:]):
            raise ValueError("incomplete multiplicative cycle")

        self.coordinate_to_polynomial = coordinate_to_polynomial
        self.polynomial_to_coordinate = polynomial_to_coordinate
        self.exponent = exponent
        self.logarithm = logarithm

    def multiply(self, left: int, right: int) -> int:
        if left == 0 or right == 0:
            return 0
        left_polynomial = self.coordinate_to_polynomial[left]
        right_polynomial = self.coordinate_to_polynomial[right]
        product = self.exponent[
            (self.logarithm[left_polynomial] + self.logarithm[right_polynomial])
            % self.modulus
        ]
        return self.polynomial_to_coordinate[product]

    def multiply_shift_reduce(self, left: int, right: int) -> int:
        """Independent carryless polynomial multiply used only by validation."""
        left_polynomial = self.coordinate_to_polynomial[left]
        right_polynomial = self.coordinate_to_polynomial[right]
        product = 0
        while right_polynomial:
            if right_polynomial & 1:
                product ^= left_polynomial
            left_polynomial <<= 1
            right_polynomial >>= 1
        for bit in range(2 * self.bits - 2, self.bits - 1, -1):
            if product & (1 << bit):
                product ^= self.polynomial << (bit - self.bits)
        return self.polynomial_to_coordinate[product]

    def inverse(self, value: int) -> int:
        if value == 0:
            raise ZeroDivisionError("zero has no inverse")
        polynomial = self.coordinate_to_polynomial[value]
        inverse = self.exponent[-self.logarithm[polynomial] % self.modulus]
        return self.polynomial_to_coordinate[inverse]

    def divide(self, numerator: int, denominator: int) -> int:
        if numerator == 0:
            return 0
        return self.multiply(numerator, self.inverse(denominator))


def make_gf16() -> BinaryField:
    return BinaryField(16, GF16_POLYNOMIAL, GF16_CANTOR_BASIS)


def ceil_power_of_two(value: int) -> int:
    if value <= 0:
        raise ValueError("value must be positive")
    return 1 << (value - 1).bit_length()


def binary_prefix_blocks(length: int) -> list[tuple[int, int]]:
    blocks: list[tuple[int, int]] = []
    offset = 0
    for bit in range(length.bit_length() - 1, -1, -1):
        size = 1 << bit
        if length & size:
            if offset & (size - 1):
                raise AssertionError("prefix block is not aligned")
            blocks.append((offset, size))
            offset += size
    if offset != length:
        raise AssertionError("prefix decomposition is incomplete")
    return blocks


def dyadic_interval_blocks(start: int, length: int) -> list[tuple[int, int]]:
    end = start + length
    blocks: list[tuple[int, int]] = []
    while start < end:
        size = 1 << ((end - start).bit_length() - 1)
        while start & (size - 1):
            size >>= 1
        blocks.append((start, size))
        start += size
    return blocks


def direct_subspace_value(field: BinaryField, size: int, point: int) -> int:
    result = 1
    for root in range(size):
        result = field.multiply(result, point ^ root)
    return result


def normalizer_factors(field: BinaryField) -> list[int]:
    """Return b_j=s_j(v_j) using only the additive recurrence."""
    factors: list[int] = []
    for bit in range(field.bits):
        subspace = 1 << bit
        for lower in range(bit):
            subspace = field.multiply(subspace, subspace ^ factors[lower])
        factors.append(subspace)
    return factors


def normalized_subspaces(
    field: BinaryField, factors: Sequence[int], point: int
) -> list[int]:
    result = [0] * field.bits
    subspace = point
    for bit in range(field.bits):
        result[bit] = field.divide(subspace, factors[bit])
        subspace = field.multiply(subspace, subspace ^ factors[bit])
    return result


def lch_values(
    field: BinaryField, factors: Sequence[int], point: int, count: int
) -> list[int]:
    normalized = normalized_subspaces(field, factors, point)
    values = [1] * count
    for index in range(1, count):
        low_bit = index & -index
        values[index] = field.multiply(
            values[index ^ low_bit], normalized[low_bit.bit_length() - 1]
        )
    return values


def lch_value(
    field: BinaryField, factors: Sequence[int], point: int, index: int
) -> int:
    normalized = normalized_subspaces(field, factors, point)
    result = 1
    while index:
        low_bit = index & -index
        result = field.multiply(
            result, normalized[low_bit.bit_length() - 1]
        )
        index ^= low_bit
    return result


Matrix = list[list[int]]


def identity_matrix(size: int) -> Matrix:
    return [[int(row == column) for column in range(size)] for row in range(size)]


def invert_matrix(field: BinaryField, matrix: Sequence[Sequence[int]]) -> Matrix:
    size = len(matrix)
    work = [list(row) + identity for row, identity in zip(matrix, identity_matrix(size))]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            raise AssertionError("matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        inverse = field.inverse(work[column][column])
        if inverse != 1:
            work[column] = [field.multiply(value, inverse) for value in work[column]]
        pivot_row = work[column]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            if factor == 0:
                continue
            target = work[row]
            work[row] = [
                left ^ field.multiply(factor, right)
                for left, right in zip(target, pivot_row)
            ]
    return [row[size:] for row in work]


def matrix_multiply(
    field: BinaryField,
    left: Sequence[Sequence[int]],
    right: Sequence[Sequence[int]],
) -> Matrix:
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("incompatible matrices")
    result = [[0] * len(right[0]) for _ in left]
    for row_index, left_row in enumerate(left):
        output = result[row_index]
        for inner, coefficient in enumerate(left_row):
            if coefficient == 0:
                continue
            right_row = right[inner]
            if coefficient == 1:
                for column, value in enumerate(right_row):
                    output[column] ^= value
            else:
                for column, value in enumerate(right_row):
                    if value:
                        output[column] ^= field.multiply(coefficient, value)
    return result


def matrix_cost(matrix: Sequence[Sequence[int]]) -> dict[str, int]:
    nonzero = sum(value != 0 for row in matrix for value in row)
    nontrivial = sum(value not in (0, 1) for row in matrix for value in row)
    xors = sum(max(0, sum(value != 0 for value in row) - 1) for row in matrix)
    return {
        "entries": sum(len(row) for row in matrix),
        "nonzero": nonzero,
        "nontrivial_fixed_multiplications": nontrivial,
        "row_accumulation_xors": xors,
    }


def direct_lagrange_matrix(
    field: BinaryField, source: Sequence[int], outputs: Sequence[int]
) -> Matrix:
    inverse_denominators: list[int] = []
    for index, point in enumerate(source):
        denominator = 1
        for other_index, other in enumerate(source):
            if index != other_index:
                denominator = field.multiply(denominator, point ^ other)
        inverse_denominators.append(field.inverse(denominator))
    result: Matrix = []
    for output in outputs:
        row: list[int] = []
        for index in range(len(source)):
            numerator = 1
            for other_index, other in enumerate(source):
                if index != other_index:
                    numerator = field.multiply(numerator, output ^ other)
            row.append(field.multiply(numerator, inverse_denominators[index]))
        result.append(row)
    return result


def cross_block_cost(
    matrix: Sequence[Sequence[int]], blocks: Sequence[tuple[int, int]]
) -> dict[str, int]:
    groups = [0] * len(matrix)
    for group, (offset, size) in enumerate(blocks):
        for index in range(offset, offset + size):
            groups[index] = group
    entries = nonzero = nontrivial = 0
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


def coupling_cost(
    matrix: Sequence[Sequence[int]],
    input_blocks: Sequence[tuple[int, int]],
    output_blocks: Sequence[tuple[int, int]],
    output_origin: int,
) -> dict[str, int]:
    pairs = coupled = fully_dense = entries = nonzero = 0
    for output_offset, output_size in output_blocks:
        row_begin = output_offset - output_origin
        for input_offset, input_size in input_blocks:
            pairs += 1
            block_entries = output_size * input_size
            block_nonzero = sum(
                matrix[row][column] != 0
                for row in range(row_begin, row_begin + output_size)
                for column in range(input_offset, input_offset + input_size)
            )
            entries += block_entries
            nonzero += block_nonzero
            coupled += block_nonzero != 0
            fully_dense += block_nonzero == block_entries
    return {
        "block_pairs": pairs,
        "coupled_block_pairs": coupled,
        "fully_dense_block_pairs": fully_dense,
        "coefficient_entries": entries,
        "nonzero_coefficients": nonzero,
    }


def exact_join_case(q: int) -> dict[str, object]:
    field = make_gf16()
    factors = normalizer_factors(field)
    output_count = 17
    source = list(range(q))
    outputs = list(range(q, q + output_count))
    blocks = binary_prefix_blocks(q)

    global_evaluation = [lch_values(field, factors, point, q) for point in source]
    local_evaluation = [[0] * q for _ in range(q)]
    for offset, size in blocks:
        for point in range(offset, offset + size):
            local_evaluation[point][offset:offset + size] = lch_values(
                field, factors, point, size
            )
    join = matrix_multiply(field, invert_matrix(field, global_evaluation), local_evaluation)
    parity = [lch_values(field, factors, point, q) for point in outputs]
    if any(value == 0 for row in parity for value in row):
        raise AssertionError("global exact-profile parity map is not dense")
    effective = matrix_multiply(field, parity, join)
    direct = direct_lagrange_matrix(field, source, outputs)
    independent = matrix_multiply(field, direct, local_evaluation)
    if effective != independent:
        raise AssertionError("GF16 LCH join differs from direct Lagrange algebra")

    coupling = coupling_cost(
        effective, blocks, dyadic_interval_blocks(q, output_count), q
    )
    if coupling["coupled_block_pairs"] != coupling["block_pairs"]:
        raise AssertionError("a GF16 input/output dyadic block pair decoupled")
    join_cost = matrix_cost(join)
    effective_cost = matrix_cost(effective)
    return {
        "construction": "new_exact_prefix_profile_not_parent_wire_compatible",
        "q": q,
        "output_count": output_count,
        "blocks": [{"offset": offset, "size": size} for offset, size in blocks],
        "join": join_cost,
        "join_cross_block": cross_block_cost(join, blocks),
        "global_parity_evaluation": matrix_cost(parity),
        "effective_local_to_parity": effective_cost,
        "input_output_block_coupling": coupling,
        "straightforward_join_coefficient_bytes_u16": 2 * q * q,
        "sparse_join_record_bytes_at_8_bytes_each": 8 * join_cost["nonzero"],
        "effective_parity_coefficient_bytes_u16": 2 * output_count * q,
        "largest_local_block_slots": blocks[0][1],
        "straightforward_out_of_place_join_slots": q,
    }


def parent_wire_case(q: int) -> dict[str, object]:
    field = make_gf16()
    factors = normalizer_factors(field)
    parent = ceil_power_of_two(q)
    blocks = binary_prefix_blocks(q)
    point_contributions = [0] * parent
    active_cosets = factor_zero = factor_one = factor_nontrivial = 0
    active_points = butterflies = scale_multiplications = 0
    factor_evaluations = setup_multiplications = identity_checks = 0

    for offset, size in blocks:
        log_size = size.bit_length() - 1
        offset_products = max(0, offset.bit_count() - 1)
        for coset in range(0, parent, size):
            factor = lch_value(field, factors, coset, offset)
            factor_evaluations += 1
            setup_multiplications += 2 * field.bits + offset_products
            if factor == 0:
                factor_zero += 1
                continue
            active_cosets += 1
            active_points += size
            butterflies += (size // 2) * log_size
            if factor == 1:
                factor_one += 1
            else:
                factor_nontrivial += 1
                scale_multiplications += size
            for point in range(coset, coset + size):
                point_contributions[point] += 1

        # Directly validate representative instances of X_(o+t)=X_o*X_t.
        sample_cosets = sorted({0, (parent // (2 * size)) * size, parent - size})
        sample_locals = sorted({0, size // 2, size - 1})
        for coset in sample_cosets:
            factor = lch_value(field, factors, coset, offset)
            for local in sample_locals:
                point = coset ^ (local % size)
                expected = lch_value(field, factors, point, offset + local)
                actual = field.multiply(
                    factor, lch_value(field, factors, point, local)
                )
                if actual != expected:
                    raise AssertionError("parent LCH dyadic block identity failed")
                identity_checks += 1

    accumulation_xors = sum(max(0, value - 1) for value in point_contributions)
    padded_butterflies = (parent // 2) * (parent.bit_length() - 1)
    padded_transfers = 4 * padded_butterflies
    materialized_transfers = (
        4 * butterflies + 2 * scale_multiplications + 3 * accumulation_xors
    )
    ideal_fused_transfers = 4 * butterflies
    schedule_bytes = 32 + 16 * len(blocks) + 16 * active_cosets
    materialized_improvement = (
        0.0 if padded_transfers == 0
        else 100.0 * (padded_transfers - materialized_transfers) / padded_transfers
    )
    ideal_fused_improvement = (
        0.0 if padded_transfers == 0
        else 100.0 * (padded_transfers - ideal_fused_transfers) / padded_transfers
    )
    return {
        "construction": "existing_parent_wire_compatible_map",
        "q": q,
        "parent": parent,
        "blocks": [{"offset": offset, "size": size} for offset, size in blocks],
        "block_count": len(blocks),
        "factor_cosets": {
            "evaluated": factor_evaluations,
            "zero": factor_zero,
            "one": factor_one,
            "nontrivial": factor_nontrivial,
        },
        "active_block_cosets": active_cosets,
        "active_point_couplings": active_points,
        "kernel_butterfly_equivalents": butterflies,
        "padded_kernel_butterfly_equivalents": padded_butterflies,
        "block_factor_fixed_multiplications": scale_multiplications,
        "cross_block_accumulation_xors": accumulation_xors,
        "materialized_candidate_logical_shard_transfers": materialized_transfers,
        "ideal_fused_lower_bound_logical_shard_transfers": ideal_fused_transfers,
        "padded_logical_shard_transfers": padded_transfers,
        "materialized_memory_improvement_percent": round(
            materialized_improvement, 9
        ),
        "ideal_fused_upper_bound_improvement_percent": round(
            ideal_fused_improvement, 9
        ),
        "largest_local_scratch_slots": blocks[0][1],
        "padded_scratch_slots": parent,
        "modeled_schedule_bytes": schedule_bytes,
        "schedule_factor_evaluation_field_multiplications": setup_multiplications,
        "block_identity_samples_checked": identity_checks,
    }


def validation_summary() -> dict[str, object]:
    field = make_gf16()
    factors = normalizer_factors(field)
    direct_checks = 0
    for bit, factor in enumerate(factors):
        expected = direct_subspace_value(field, 1 << bit, 1 << bit)
        if factor != expected:
            raise AssertionError("GF16 subspace recurrence differs from direct product")
        direct_checks += 1
    recurrence_checks = 0
    for point in (0xA5A5, 0xFFFF):
        normalized = normalized_subspaces(field, factors, point)
        for bit, factor in enumerate(factors):
            recurrence_value = field.multiply(normalized[bit], factor)
            direct_value = direct_subspace_value(field, 1 << bit, point)
            if recurrence_value != direct_value:
                raise AssertionError(
                    "GF16 subspace recurrence differs away from its normalizer"
                )
            recurrence_checks += 1

    samples = (1, 2, 3, 0x1234, 0x8000, 0xFFFF)
    field_checks = 0
    shift_reduce_checks = 0
    for value in samples:
        if field.multiply(value, field.inverse(value)) != 1:
            raise AssertionError("GF16 inverse identity failed")
        for other in samples:
            if field.multiply(value, other) != field.multiply(other, value):
                raise AssertionError("GF16 multiplication is not commutative")
            if field.multiply(value, other) != field.multiply_shift_reduce(
                value, other
            ):
                raise AssertionError("GF16 log multiply differs from shift/reduce")
            field_checks += 1
            shift_reduce_checks += 1
    return {
        "normalizer_direct_product_checks": direct_checks,
        "subspace_recurrence_direct_product_checks": recurrence_checks,
        "field_multiply_inverse_checks": field_checks,
        "shift_reduce_multiply_checks": shift_reduce_checks,
        "coordinate_samples": list(samples),
    }


def traffic_grid(parent_cases: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for case in parent_cases:
        q = int(case["q"])
        schedule_bytes = int(case["modeled_schedule_bytes"])
        baseline_transfers = int(case["padded_logical_shard_transfers"])
        materialized_transfers = int(
            case["materialized_candidate_logical_shard_transfers"]
        )
        ideal_fused_transfers = int(
            case["ideal_fused_lower_bound_logical_shard_transfers"]
        )
        candidate_scratch_slots = int(case["largest_local_scratch_slots"])
        padded_scratch_slots = int(case["padded_scratch_slots"])
        for shard_bytes in SHARD_BYTES:
            baseline_stripe = baseline_transfers * shard_bytes
            materialized_stripe = materialized_transfers * shard_bytes
            ideal_fused_stripe = ideal_fused_transfers * shard_bytes
            candidate_scratch = candidate_scratch_slots * shard_bytes
            padded_scratch = padded_scratch_slots * shard_bytes
            for batch in BATCHES:
                for reuse in REUSES:
                    stripes = batch * reuse
                    baseline_total = baseline_stripe * stripes
                    materialized_total = (
                        materialized_stripe * stripes + schedule_bytes
                    )
                    ideal_fused_total = ideal_fused_stripe * stripes + schedule_bytes
                    materialized_improvement = 100.0 * (
                        baseline_total - materialized_total
                    ) / baseline_total
                    ideal_fused_improvement = 100.0 * (
                        baseline_total - ideal_fused_total
                    ) / baseline_total
                    result.append({
                        "q": q,
                        "shard_bytes": shard_bytes,
                        "batch": batch,
                        "plan_reuse_count": reuse,
                        "baseline_execution_bytes": baseline_total,
                        "materialized_execution_plus_incremental_setup_bytes": (
                            materialized_total
                        ),
                        "ideal_fused_execution_plus_incremental_setup_bytes": (
                            ideal_fused_total
                        ),
                        "candidate_scratch_bytes_per_stripe": candidate_scratch,
                        "padded_scratch_bytes_per_stripe": padded_scratch,
                        "candidate_scratch_bytes_if_batch_resident": (
                            candidate_scratch * batch
                        ),
                        "padded_scratch_bytes_if_batch_resident": padded_scratch * batch,
                        "incremental_setup_bytes_per_stripe": round(
                            schedule_bytes / stripes, 9
                        ),
                        "materialized_memory_improvement_percent": round(
                            materialized_improvement, 9
                        ),
                        "ideal_fused_upper_bound_improvement_percent": round(
                            ideal_fused_improvement, 9
                        ),
                    })
    return result


def run_task(task: tuple[str, int]) -> tuple[str, int, dict[str, object]]:
    kind, q = task
    if kind == "parent":
        return kind, q, parent_wire_case(q)
    if kind == "join":
        return kind, q, exact_join_case(q)
    raise ValueError("unknown task")


def allowed_cpu_count() -> int:
    if hasattr(os, "sched_getaffinity"):
        try:
            return len(os.sched_getaffinity(0))
        except OSError:
            pass
    return os.cpu_count() or 1


def experiment(workers: int) -> dict[str, object]:
    validation = validation_summary()
    tasks = [("parent", q) for q in PARENT_Q] + [("join", q) for q in JOIN_Q]
    rows: list[tuple[str, int, dict[str, object]]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for row in executor.map(run_task, tasks):
            rows.append(row)
    parent_cases = [row for kind, _, row in rows if kind == "parent"]
    exact_cases = [row for kind, _, row in rows if kind == "join"]
    parent_cases.sort(key=lambda row: int(row["q"]))
    exact_cases.sort(key=lambda row: int(row["q"]))
    grid = traffic_grid(parent_cases)

    best_materialized_geometry = max(
        parent_cases,
        key=lambda row: float(row["materialized_memory_improvement_percent"]),
    )
    best_ideal_geometry = max(
        parent_cases,
        key=lambda row: float(row["ideal_fused_upper_bound_improvement_percent"]),
    )
    best_materialized_grid = max(
        grid, key=lambda row: float(row["materialized_memory_improvement_percent"])
    )
    worst_materialized_grid = min(
        grid, key=lambda row: float(row["materialized_memory_improvement_percent"])
    )
    best_ideal_grid = max(
        grid, key=lambda row: float(row["ideal_fused_upper_bound_improvement_percent"])
    )
    all_pairs = sum(
        int(case["input_output_block_coupling"]["block_pairs"])
        for case in exact_cases
    )
    coupled_pairs = sum(
        int(case["input_output_block_coupling"]["coupled_block_pairs"])
        for case in exact_cases
    )
    if coupled_pairs != all_pairs:
        raise AssertionError("representative GF16 block coupling is incomplete")

    return {
        "schema": SCHEMA,
        "scope": {
            "field": "legacy GF(2^16)",
            "timing_performed": False,
            "production_code_changed": False,
            "parent_candidate_changes_wire_profile": False,
            "exact_prefix_candidate_requires_new_profile_identifier": True,
            "worker_count_excluded_from_deterministic_identity": True,
        },
        "coordinate_system": {
            "polynomial": hex(GF16_POLYNOMIAL),
            "cantor_coordinate_basis": list(GF16_CANTOR_BASIS),
        },
        "validation": validation,
        "accounting_contract": {
            "butterfly_logical_transfers": "two shard reads plus two shard writes",
            "fixed_multiply_logical_transfers": "one shard read plus one shard write",
            "accumulation_xor_logical_transfers": "source and destination reads plus destination write",
            "materialized_candidate": "optimistic transform-stage traffic plus separate fixed-multiply and accumulation passes; repeated input staging, final output placement, packing, scatter, allocator, and cache-line effects excluded",
            "ideal_fused_lower_bound": "butterfly traffic only; all factor, tail-input, and accumulation traffic is unrealistically free and therefore only an upper bound on possible improvement",
            "padded_baseline": "optimistic transform-stage traffic under the same exclusions",
            "schedule_record_model": "32-byte header plus 16 bytes per block and active coset",
            "setup_amortization": "one incremental candidate schedule divided over batch * plan_reuse_count stripes",
            "scratch_slots": "full shard slots excluding caller input/output",
        },
        "axes": {
            "parent_q": list(PARENT_Q),
            "exact_join_q": list(JOIN_Q),
            "shard_bytes": list(SHARD_BYTES),
            "batch": list(BATCHES),
            "plan_reuse_count": list(REUSES),
        },
        "parent_wire_block_accumulation": parent_cases,
        "exact_prefix_join": {
            "cases": exact_cases,
            "input_output_block_pairs": all_pairs,
            "coupled_input_output_block_pairs": coupled_pairs,
        },
        "byte_batch_reuse_memory_grid": grid,
        "summary": {
            "best_materialized_parent_geometry": best_materialized_geometry,
            "best_ideal_fusion_parent_geometry": best_ideal_geometry,
            "best_materialized_amortized_memory_cell": best_materialized_grid,
            "worst_materialized_amortized_memory_cell": worst_materialized_grid,
            "best_ideal_fusion_amortized_memory_cell": best_ideal_grid,
            "traffic_cells": len(grid),
            "exact_gf16_join_cases": len(exact_cases),
        },
        "disposition": {
            "materialized_parent_wire_block_accumulation": "kill",
            "parent_wire_reason": (
                "no representative GF16 geometry reaches the 10 percent promotion "
                "threshold with separate materialized accumulation; the only symbolic "
                "ceiling above 10 percent assumes free fusion, factor, tail-input, and "
                "accumulation traffic and belongs in the general C1/C2 parent scheduler"
            ),
            "fused_parent_wire_schedule": "not_a_distinct_c5_route_send_to_c1_c2",
            "parent_wire_low_memory_note": (
                "largest-block scratch can be smaller than the parent, but C1/C2 "
                "pruned or tiled parent schedules should own that opportunity"
            ),
            "naive_exact_join_plus_parent": "kill",
            "exact_reason": (
                "all representative GF16 input/output block pairs remain coupled, "
                "the global parity map is dense, and explicit join storage/execution "
                "is additional work; a different factored exact algorithm is not disproved"
            ),
            "broader_factored_exact_profile": "inconclusive_route_to_c6_c7_c8",
            "promotion": "none",
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
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--workers", type=int,
        default=min(128, allowed_cpu_count(), len(PARENT_Q) + len(JOIN_Q)),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.workers <= 0:
        raise SystemExit("--workers must be positive")
    result = experiment(args.workers)
    write_json(args.output, result)
    print(json.dumps({
        "output": str(args.output),
        "schema": result["schema"],
        "traffic_cells": result["summary"]["traffic_cells"],
        "exact_gf16_join_cases": result["summary"]["exact_gf16_join_cases"],
        "promotion": result["disposition"]["promotion"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
