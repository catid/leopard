#!/usr/bin/env python3
"""Compact, parent-preserving truncated LCH transform experiment.

This is a scalar research model, not a production kernel.  It evaluates the
same dyadic parent transform as Leopard while accepting only declared active
input coordinates and returning only requested output coordinates.  Boundary
nodes use compact maps; a dense vector is created only for a complete dyadic
subtransform whose input and output masks are both full.

The forward and inverse recurrences are independently checked against a direct
LCH evaluation matrix (and its Gaussian-elimination inverse), as well as an
explicit padded butterfly transform.  No timing command is provided: C2 is a
correctness and operation/memory-count checkpoint.
"""

from __future__ import annotations

import argparse
import functools
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence


SCHEMA_VERSION = "leopard2-c2-v1"


class BinaryField:
    """Table-backed GF(2^m) with an explicit public coordinate basis."""

    def __init__(
        self, bits: int, polynomial: int, coordinate_basis: Sequence[int]
    ) -> None:
        if bits <= 0 or bits > 8 or len(coordinate_basis) != bits:
            raise ValueError("the scalar checkpoint supports GF(2^m), m <= 8")
        if polynomial & (1 << bits) == 0:
            raise ValueError("the polynomial has the wrong degree")
        self.bits = bits
        self.order = 1 << bits
        self.polynomial = polynomial
        self.coordinate_basis = tuple(coordinate_basis)

        coordinate_to_polynomial = [0] * self.order
        polynomial_to_coordinate = [-1] * self.order
        for coordinate in range(self.order):
            value = 0
            for bit, basis_value in enumerate(self.coordinate_basis):
                if coordinate & (1 << bit):
                    value ^= basis_value
            if value >= self.order or polynomial_to_coordinate[value] >= 0:
                raise ValueError("coordinate basis is not independent")
            coordinate_to_polynomial[coordinate] = value
            polynomial_to_coordinate[value] = coordinate
        self.coordinate_to_polynomial = tuple(coordinate_to_polynomial)
        self.polynomial_to_coordinate = tuple(polynomial_to_coordinate)

        rows: list[bytes] = []
        for left in range(self.order):
            rows.append(bytes(self._multiply_slow(left, right)
                              for right in range(self.order)))
        self.mul_rows = tuple(rows)

    def _multiply_slow(self, left: int, right: int) -> int:
        a = self.coordinate_to_polynomial[left]
        b = self.coordinate_to_polynomial[right]
        product = 0
        while b:
            if b & 1:
                product ^= a
            b >>= 1
            a <<= 1
        for bit in range(self.bits * 2 - 2, self.bits - 1, -1):
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
    # Leopard's polynomial representation and public Cantor-coordinate basis.
    return BinaryField(8, 0x11D, (1, 214, 152, 146, 86, 200, 88, 230))


Polynomial = list[int]


def polynomial_multiply(
    field: BinaryField, left: Sequence[int], right: Sequence[int]
) -> Polynomial:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
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
    field: BinaryField, length: int
) -> tuple[tuple[int, ...], ...]:
    dimension = _validate_length(field, length)
    subspaces = tuple(subspace_polynomial(field, bit)
                      for bit in range(dimension))
    result: list[tuple[int, ...]] = []
    for index in range(length):
        polynomial = [1]
        normalizer = 1
        for bit, subspace in enumerate(subspaces):
            if not index & (1 << bit):
                continue
            polynomial = polynomial_multiply(field, polynomial, subspace)
            normalizer = field.multiply(
                normalizer, polynomial_evaluate(field, subspace, 1 << bit)
            )
        inverse = field.inverse(normalizer)
        result.append(tuple(field.multiply(value, inverse)
                            for value in polynomial))
    return tuple(result)


@functools.lru_cache(maxsize=None)
def direct_lch_matrix(
    field: BinaryField, length: int, shift: int
) -> tuple[tuple[int, ...], ...]:
    """Independent LCH-to-Lagrange matrix built from monomial polynomials."""
    _validate_geometry(field, length, shift)
    basis = normalized_lch_polynomials(field, length)
    rows: list[tuple[int, ...]] = []
    for output in range(length):
        point = shift ^ output
        rows.append(tuple(polynomial_evaluate(field, polynomial, point)
                          for polynomial in basis))
    return tuple(rows)


def invert_matrix(
    field: BinaryField, matrix: Sequence[Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    size = len(matrix)
    work = [list(row) + [int(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    if any(len(row) != 2 * size for row in work):
        raise ValueError("matrix must be square")
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        inverse = field.inverse(work[column][column])
        work[column] = [field.multiply(value, inverse)
                        for value in work[column]]
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            factor = work[row][column]
            scaled = [field.multiply(factor, value)
                      for value in work[column]]
            work[row] = [left ^ right
                         for left, right in zip(work[row], scaled)]
    return tuple(tuple(row[size:]) for row in work)


@functools.lru_cache(maxsize=None)
def direct_inverse_lch_matrix(
    field: BinaryField, length: int, shift: int
) -> tuple[tuple[int, ...], ...]:
    return invert_matrix(field, direct_lch_matrix(field, length, shift))


def matrix_vector(
    field: BinaryField,
    matrix: Sequence[Sequence[int]],
    vector: Sequence[int],
) -> list[int]:
    result: list[int] = []
    for row in matrix:
        value = 0
        for coefficient, item in zip(row, vector):
            if coefficient and item:
                value ^= item if coefficient == 1 else field.multiply(
                    coefficient, item
                )
        result.append(value)
    return result


def make_fft_skew(field: BinaryField) -> tuple[int, ...]:
    """Reconstruct Leopard's FFTSkew values before logarithm conversion."""
    skew = [0] * (field.order - 1)
    temp = [1 << bit for bit in range(1, field.bits)]
    for layer in range(field.bits - 1):
        step = 1 << (layer + 1)
        first = (1 << layer) - 1
        skew[first] = 0
        for index in range(layer, field.bits - 1):
            offset = 1 << (index + 1)
            for position in range(first, offset, step):
                skew[position + offset] = skew[position] ^ temp[index]
        scale = field.inverse(field.multiply(temp[layer], temp[layer] ^ 1))
        for index in range(layer + 1, field.bits - 1):
            temp[index] = field.multiply(
                field.multiply(temp[index], temp[index] ^ 1), scale
            )
    return tuple(skew)


def _validate_length(field: BinaryField, length: int) -> int:
    if length <= 0 or length & (length - 1):
        raise ValueError("parent length must be a positive power of two")
    dimension = length.bit_length() - 1
    if dimension > field.bits:
        raise ValueError("parent length exceeds the field")
    return dimension


def _validate_geometry(field: BinaryField, length: int, shift: int) -> None:
    _validate_length(field, length)
    if shift < 0 or shift & (length - 1) or shift + length > field.order:
        raise ValueError("shift must name an aligned in-field additive coset")


def _mask(value: int | Sequence[bool], length: int, name: str) -> tuple[bool, ...]:
    if isinstance(value, int):
        if value < 0 or value > length:
            raise ValueError(f"{name} prefix is outside the parent")
        return tuple(index < value for index in range(length))
    result = tuple(bool(item) for item in value)
    if len(result) != length:
        raise ValueError(f"{name} mask has the wrong length")
    return result


@dataclass
class OperationCounter:
    fixed_multiplications: int = 0
    xors: int = 0
    loads: int = 0
    stores: int = 0
    complete_blocks: int = 0
    complete_butterflies: int = 0
    boundary_nodes: int = 0
    boundary_pairs: int = 0
    compact_values_written: int = 0
    maximum_dense_block: int = 0
    partial_dense_materializations: int = 0


def _constant_product(
    field: BinaryField, coefficient: int, value: int,
    counter: OperationCounter | None,
) -> int:
    if coefficient == 0:
        return 0
    if coefficient == 1:
        return value
    if counter is not None:
        counter.fixed_multiplications += 1
    return field.multiply(coefficient, value)


def _xor(left: int, right: int, counter: OperationCounter | None) -> int:
    if counter is not None:
        counter.xors += 1
    return left ^ right


def _linear_pair(
    field: BinaryField,
    left: int,
    right: int,
    left_coefficient: int,
    right_coefficient: int,
    counter: OperationCounter | None,
) -> int:
    terms: list[int] = []
    if left_coefficient:
        terms.append(_constant_product(field, left_coefficient, left, counter))
    if right_coefficient:
        terms.append(_constant_product(field, right_coefficient, right, counter))
    if not terms:
        return 0
    return terms[0] if len(terms) == 1 else _xor(terms[0], terms[1], counter)


def _padded_forward(
    field: BinaryField,
    values: Sequence[int],
    shift: int,
    skew: Sequence[int],
    counter: OperationCounter | None = None,
) -> list[int]:
    length = len(values)
    if length == 1:
        return list(values)
    half = length // 2
    multiplier = skew[shift + half - 1]
    left: list[int] = []
    right: list[int] = []
    for offset in range(half):
        a, b = values[offset], values[half + offset]
        if counter is not None:
            counter.loads += 2
            counter.stores += 2
            counter.complete_butterflies += 1
        if multiplier == 0:
            u = a
            v = _xor(a, b, counter)
        elif multiplier == 1:
            u = _xor(a, b, counter)
            v = a
        else:
            product = _constant_product(field, multiplier, b, counter)
            u = _xor(a, product, counter)
            v = _xor(b, u, counter)
        left.append(u)
        right.append(v)
    return (_padded_forward(field, left, shift, skew, counter) +
            _padded_forward(field, right, shift + half, skew, counter))


def _padded_inverse(
    field: BinaryField,
    values: Sequence[int],
    shift: int,
    skew: Sequence[int],
    counter: OperationCounter | None = None,
) -> list[int]:
    length = len(values)
    if length == 1:
        return list(values)
    half = length // 2
    multiplier = skew[shift + half - 1]
    left = _padded_inverse(field, values[:half], shift, skew, counter)
    right = _padded_inverse(field, values[half:], shift + half, skew, counter)
    result = [0] * length
    for offset, (u, v) in enumerate(zip(left, right)):
        if counter is not None:
            counter.loads += 2
            counter.stores += 2
            counter.complete_butterflies += 1
        if multiplier == 0:
            x = u
            y = _xor(u, v, counter)
        elif multiplier == 1:
            x = v
            y = _xor(u, v, counter)
        else:
            y = _xor(u, v, counter)
            x = _xor(u, _constant_product(field, multiplier, y, counter), counter)
        result[offset] = x
        result[half + offset] = y
    return result


@dataclass
class TruncatedNode:
    length: int
    shift: int
    direction: str
    input_mask: tuple[bool, ...]
    output_mask: tuple[bool, ...]
    required_inputs: tuple[bool, ...]
    kind: str
    multiplier: int = 0
    left: "TruncatedNode | None" = None
    right: "TruncatedNode | None" = None
    left_inputs: tuple[bool, ...] = ()
    right_inputs: tuple[bool, ...] = ()
    live_outputs: tuple[bool, ...] = ()


@dataclass
class TruncatedPlan:
    field: BinaryField
    parent_size: int
    shift: int
    direction: str
    input_mask: tuple[bool, ...]
    output_mask: tuple[bool, ...]
    root: TruncatedNode
    skew: tuple[int, ...]


def _compile_forward_node(
    field: BinaryField,
    skew: Sequence[int],
    length: int,
    shift: int,
    input_mask: tuple[bool, ...],
    output_mask: tuple[bool, ...],
) -> TruncatedNode:
    if not any(input_mask) or not any(output_mask):
        return TruncatedNode(length, shift, "forward", input_mask, output_mask,
                             (False,) * length, "empty",
                             live_outputs=(False,) * length)
    if length == 1:
        required = (bool(input_mask[0] and output_mask[0]),)
        return TruncatedNode(length, shift, "forward", input_mask, output_mask,
                             required, "leaf", live_outputs=required)
    if all(input_mask) and all(output_mask):
        return TruncatedNode(length, shift, "forward", input_mask, output_mask,
                             input_mask, "complete", live_outputs=output_mask)

    half = length // 2
    multiplier = skew[shift + half - 1]
    a_live, b_live = input_mask[:half], input_mask[half:]
    u_live = tuple(a or (b and multiplier != 0)
                   for a, b in zip(a_live, b_live))
    v_live = tuple(a or (b and (multiplier ^ 1) != 0)
                   for a, b in zip(a_live, b_live))
    left = _compile_forward_node(
        field, skew, half, shift, u_live, output_mask[:half]
    )
    right = _compile_forward_node(
        field, skew, half, shift + half, v_live, output_mask[half:]
    )
    required_a: list[bool] = []
    required_b: list[bool] = []
    for index in range(half):
        need_u = left.required_inputs[index]
        need_v = right.required_inputs[index]
        required_a.append(a_live[index] and (need_u or need_v))
        required_b.append(b_live[index] and (
            (need_u and multiplier != 0) or
            (need_v and (multiplier ^ 1) != 0)
        ))
    required = tuple(required_a + required_b)
    live_outputs = left.live_outputs + right.live_outputs
    if not any(required):
        return TruncatedNode(
            length, shift, "forward", input_mask, output_mask, required,
            "empty", live_outputs=(False,) * length
        )
    return TruncatedNode(
        length, shift, "forward", input_mask, output_mask,
        required, "boundary", multiplier, left, right,
        left.required_inputs, right.required_inputs, live_outputs
    )


def _compile_inverse_node(
    field: BinaryField,
    skew: Sequence[int],
    length: int,
    shift: int,
    input_mask: tuple[bool, ...],
    output_mask: tuple[bool, ...],
) -> TruncatedNode:
    if not any(input_mask) or not any(output_mask):
        return TruncatedNode(length, shift, "inverse", input_mask, output_mask,
                             (False,) * length, "empty",
                             live_outputs=(False,) * length)
    if length == 1:
        required = (bool(input_mask[0] and output_mask[0]),)
        return TruncatedNode(length, shift, "inverse", input_mask, output_mask,
                             required, "leaf", live_outputs=required)
    if all(input_mask) and all(output_mask):
        return TruncatedNode(length, shift, "inverse", input_mask, output_mask,
                             input_mask, "complete", live_outputs=output_mask)

    half = length // 2
    multiplier = skew[shift + half - 1]
    need_x, need_y = output_mask[:half], output_mask[half:]
    need_u = tuple(y or (x and (multiplier ^ 1) != 0)
                   for x, y in zip(need_x, need_y))
    need_v = tuple(y or (x and multiplier != 0)
                   for x, y in zip(need_x, need_y))
    left = _compile_inverse_node(
        field, skew, half, shift, input_mask[:half], need_u
    )
    right = _compile_inverse_node(
        field, skew, half, shift + half, input_mask[half:], need_v
    )
    required = left.required_inputs + right.required_inputs
    live_x = tuple(
        requested and (
            (left.live_outputs[index] and (multiplier ^ 1) != 0) or
            (right.live_outputs[index] and multiplier != 0)
        )
        for index, requested in enumerate(need_x)
    )
    live_y = tuple(
        requested and (
            left.live_outputs[index] or right.live_outputs[index]
        )
        for index, requested in enumerate(need_y)
    )
    live_outputs = live_x + live_y
    if not any(required):
        return TruncatedNode(
            length, shift, "inverse", input_mask, output_mask, required,
            "empty", live_outputs=(False,) * length
        )
    return TruncatedNode(
        length, shift, "inverse", input_mask, output_mask,
        required, "boundary", multiplier, left, right,
        left.live_outputs, right.live_outputs, live_outputs
    )


def compile_truncated_plan(
    field: BinaryField,
    parent_size: int,
    active_inputs: int | Sequence[bool],
    requested_outputs: int | Sequence[bool],
    *,
    shift: int = 0,
    direction: str = "forward",
) -> TruncatedPlan:
    _validate_geometry(field, parent_size, shift)
    input_mask = _mask(active_inputs, parent_size, "active input")
    output_mask = _mask(requested_outputs, parent_size, "requested output")
    skew = make_fft_skew(field)
    if direction == "forward":
        root = _compile_forward_node(
            field, skew, parent_size, shift, input_mask, output_mask
        )
    elif direction == "inverse":
        root = _compile_inverse_node(
            field, skew, parent_size, shift, input_mask, output_mask
        )
    else:
        raise ValueError("direction must be forward or inverse")
    return TruncatedPlan(
        field, parent_size, shift, direction, input_mask, output_mask, root,
        tuple(skew)
    )


def _forward_boundary_pair(
    field: BinaryField,
    multiplier: int,
    a: int,
    b: int,
    have_a: bool,
    have_b: bool,
    need_u: bool,
    need_v: bool,
    counter: OperationCounter | None,
) -> tuple[int | None, int | None]:
    if counter is not None:
        counter.loads += int(have_a) + int(have_b)
        counter.stores += int(need_u) + int(need_v)
        counter.compact_values_written += int(need_u) + int(need_v)
        counter.boundary_pairs += 1
    if need_u and need_v and have_a and have_b:
        if multiplier == 0:
            return a, _xor(a, b, counter)
        if multiplier == 1:
            return _xor(a, b, counter), a
        u = _xor(a, _constant_product(field, multiplier, b, counter), counter)
        return u, _xor(b, u, counter)
    u = (_linear_pair(field, a, b, int(have_a),
                      multiplier if have_b else 0, counter)
         if need_u else None)
    v = (_linear_pair(field, a, b, int(have_a),
                      (multiplier ^ 1) if have_b else 0, counter)
         if need_v else None)
    return u, v


def _inverse_boundary_pair(
    field: BinaryField,
    multiplier: int,
    u: int,
    v: int,
    have_u: bool,
    have_v: bool,
    need_x: bool,
    need_y: bool,
    counter: OperationCounter | None,
) -> tuple[int | None, int | None]:
    if counter is not None:
        counter.loads += int(have_u) + int(have_v)
        counter.stores += int(need_x) + int(need_y)
        counter.compact_values_written += int(need_x) + int(need_y)
        counter.boundary_pairs += 1
    if need_x and need_y and have_u and have_v:
        if multiplier == 0:
            return u, _xor(u, v, counter)
        if multiplier == 1:
            return v, _xor(u, v, counter)
        y = _xor(u, v, counter)
        return _xor(u, _constant_product(field, multiplier, y, counter), counter), y
    if need_x and not need_y and have_u and have_v:
        if multiplier == 0:
            return u, None
        if multiplier == 1:
            return v, None
        delta = _xor(u, v, counter)
        return _xor(
            u, _constant_product(field, multiplier, delta, counter), counter
        ), None
    x = (_linear_pair(field, u, v,
                      (multiplier ^ 1) if have_u else 0,
                      multiplier if have_v else 0, counter)
         if need_x else None)
    y = (_linear_pair(field, u, v, int(have_u), int(have_v), counter)
         if need_y else None)
    return x, y


def _execute_node(
    plan: TruncatedPlan,
    node: TruncatedNode,
    values: dict[int, int],
    counter: OperationCounter | None,
) -> dict[int, int]:
    if node.kind == "empty":
        return {}
    if node.kind == "leaf":
        return {0: values.get(0, 0)} if node.output_mask[0] else {}
    if node.kind == "complete":
        dense = [values[index] for index in range(node.length)]
        if counter is not None:
            counter.complete_blocks += 1
            counter.maximum_dense_block = max(
                counter.maximum_dense_block, node.length
            )
            if not all(node.input_mask) or not all(node.output_mask):
                counter.partial_dense_materializations += 1
        transformed = (
            _padded_forward(plan.field, dense, node.shift, plan.skew, counter)
            if node.direction == "forward" else
            _padded_inverse(plan.field, dense, node.shift, plan.skew, counter)
        )
        return dict(enumerate(transformed))

    if counter is not None:
        counter.boundary_nodes += 1
    half = node.length // 2
    assert node.left is not None and node.right is not None
    if node.direction == "forward":
        left_values: dict[int, int] = {}
        right_values: dict[int, int] = {}
        for index, (need_u, need_v) in enumerate(zip(
            node.left_inputs, node.right_inputs
        )):
            if not need_u and not need_v:
                continue
            have_a = node.required_inputs[index]
            have_b = node.required_inputs[half + index]
            u, v = _forward_boundary_pair(
                plan.field, node.multiplier,
                values.get(index, 0), values.get(half + index, 0),
                have_a, have_b, need_u, need_v, counter
            )
            if need_u:
                left_values[index] = 0 if u is None else u
            if need_v:
                right_values[index] = 0 if v is None else v
        left_output = _execute_node(plan, node.left, left_values, counter)
        right_output = _execute_node(plan, node.right, right_values, counter)
        result = dict(left_output)
        result.update((half + index, value)
                      for index, value in right_output.items())
        return result

    left_input = {index: values[index]
                  for index in range(half) if index in values}
    right_input = {index: values[half + index]
                   for index in range(half) if half + index in values}
    left_output = _execute_node(plan, node.left, left_input, counter)
    right_output = _execute_node(plan, node.right, right_input, counter)
    result: dict[int, int] = {}
    for index, (need_x, need_y) in enumerate(zip(
        node.live_outputs[:half], node.live_outputs[half:]
    )):
        if not need_x and not need_y:
            continue
        have_u = node.left_inputs[index]
        have_v = node.right_inputs[index]
        x, y = _inverse_boundary_pair(
            plan.field, node.multiplier,
            left_output.get(index, 0), right_output.get(index, 0),
            have_u, have_v, need_x, need_y, counter
        )
        if need_x:
            result[index] = 0 if x is None else x
        if need_y:
            result[half + index] = 0 if y is None else y
    return result


def execute_truncated(
    plan: TruncatedPlan,
    active_values: Sequence[int],
    *,
    counter: OperationCounter | None = None,
) -> tuple[int, ...]:
    expected = sum(plan.input_mask)
    if len(active_values) != expected:
        raise ValueError("compact input count does not match the active mask")
    if any(value < 0 or value >= plan.field.order for value in active_values):
        raise ValueError("input symbol is outside the field")
    # Values irrelevant to the requested outputs never enter the recursive
    # executor.  This is plan-time projection, not a byte-loop mask branch.
    required: dict[int, int] = {}
    compact_index = 0
    for index, active in enumerate(plan.input_mask):
        if not active:
            continue
        if plan.root.required_inputs[index]:
            required[index] = active_values[compact_index]
        compact_index += 1
    output = _execute_node(plan, plan.root, required, counter)
    return tuple(output.get(index, 0)
                 for index, requested in enumerate(plan.output_mask)
                 if requested)


def lch_truncated_forward(
    field: BinaryField,
    active_values: Sequence[int],
    active_inputs: int | Sequence[bool],
    requested_outputs: int | Sequence[bool],
    parent_size: int,
    shift: int = 0,
) -> tuple[int, ...]:
    plan = compile_truncated_plan(
        field, parent_size, active_inputs, requested_outputs,
        shift=shift, direction="forward"
    )
    return execute_truncated(plan, active_values)


def lch_truncated_inverse(
    field: BinaryField,
    active_values: Sequence[int],
    active_inputs: int | Sequence[bool],
    requested_outputs: int | Sequence[bool],
    parent_size: int,
    shift: int = 0,
) -> tuple[int, ...]:
    plan = compile_truncated_plan(
        field, parent_size, active_inputs, requested_outputs,
        shift=shift, direction="inverse"
    )
    return execute_truncated(plan, active_values)


@dataclass
class TestCounts:
    gf2_4_prefix_plans: int = 0
    gf2_4_sparse_plans: int = 0
    gf8_prefix_plans: int = 0
    gf8_sparse_plans: int = 0
    basis_vectors: int = 0
    zero_vectors: int = 0
    padded_vectors: int = 0
    direct_symbol_comparisons: int = 0
    padded_symbol_comparisons: int = 0
    compact_contract_rejections: int = 0
    partial_dense_materializations: int = 0


def _indices(mask: Sequence[bool]) -> tuple[int, ...]:
    return tuple(index for index, enabled in enumerate(mask) if enabled)


def _verify_plan(
    field: BinaryField,
    length: int,
    shift: int,
    input_mask: tuple[bool, ...],
    output_mask: tuple[bool, ...],
    direction: str,
    counts: TestCounts,
    *,
    check_padded: bool = True,
) -> None:
    plan = compile_truncated_plan(
        field, length, input_mask, output_mask, shift=shift,
        direction=direction
    )
    inputs = _indices(input_mask)
    outputs = _indices(output_mask)
    matrix = (direct_lch_matrix(field, length, shift)
              if direction == "forward" else
              direct_inverse_lch_matrix(field, length, shift))

    zero = execute_truncated(plan, [0] * len(inputs))
    if zero != (0,) * len(outputs):
        raise AssertionError("zero vector did not remain zero")
    counts.zero_vectors += 1
    counts.direct_symbol_comparisons += len(outputs)

    # Equality on every active unit vector proves equality of these GF-linear
    # maps for every possible field-valued input, without enumerating 16^N.
    for compact_index, input_index in enumerate(inputs):
        values = [0] * len(inputs)
        values[compact_index] = 1
        actual = execute_truncated(plan, values)
        expected = tuple(matrix[output][input_index] for output in outputs)
        if actual != expected:
            raise AssertionError(
                f"direct oracle mismatch: GF(2^{field.bits}), N={length}, "
                f"shift={shift}, {direction}, input={input_index}"
            )
        counts.basis_vectors += 1
        counts.direct_symbol_comparisons += len(outputs)

    if check_padded:
        compact = [((index + 1) * 7 + shift * 3 + len(inputs)) % field.order
                   for index in range(len(inputs))]
        full = [0] * length
        for index, value in zip(inputs, compact):
            full[index] = value
        padded = (_padded_forward(field, full, shift, plan.skew)
                  if direction == "forward" else
                  _padded_inverse(field, full, shift, plan.skew))
        actual = execute_truncated(plan, compact)
        expected = tuple(padded[index] for index in outputs)
        if actual != expected:
            raise AssertionError("padded transform mismatch")
        counts.padded_vectors += 1
        counts.padded_symbol_comparisons += len(outputs)

    # Dense allocations are permitted only for full/full dyadic subtrees.
    counter = OperationCounter()
    execute_truncated(plan, [1] * len(inputs), counter=counter)
    counts.partial_dense_materializations += counter.partial_dense_materializations
    if counter.partial_dense_materializations:
        raise AssertionError("an incomplete subtree was materialized densely")


def _representative_metrics(
    field: BinaryField,
    length: int,
    active: int,
    requested: int,
    direction: str,
) -> dict[str, object]:
    plan = compile_truncated_plan(
        field, length, active, requested, direction=direction
    )
    candidate = OperationCounter()
    execute_truncated(plan, [1] * active, counter=candidate)
    padded = OperationCounter()
    padded.complete_blocks = 1
    padded.maximum_dense_block = length
    full = [1 if index < active else 0 for index in range(length)]
    if direction == "forward":
        _padded_forward(field, full, 0, plan.skew, padded)
    else:
        _padded_inverse(field, full, 0, plan.skew, padded)
    candidate_dict = asdict(candidate)
    padded_dict = asdict(padded)
    return {
        "direction": direction,
        "parent_size": length,
        "active_inputs": active,
        "requested_outputs": requested,
        "candidate": candidate_dict,
        "padded": padded_dict,
        "ratios": {
            "fixed_multiplications": (
                candidate.fixed_multiplications / padded.fixed_multiplications
                if padded.fixed_multiplications else 0.0
            ),
            "xors": candidate.xors / padded.xors if padded.xors else 0.0,
            "loads": candidate.loads / padded.loads if padded.loads else 0.0,
            "stores": candidate.stores / padded.stores if padded.stores else 0.0,
        },
        "workspace_contract": {
            "padded_parent_slots": length,
            "largest_dense_candidate_block": candidate.maximum_dense_block,
            "incomplete_dense_blocks": candidate.partial_dense_materializations,
            "compact_boundary_values_written": candidate.compact_values_written,
        },
    }


def run_self_test() -> dict[str, object]:
    counts = TestCounts()
    gf2_4 = make_gf2_4()
    gf8 = make_legacy_gf8()

    # Every prefix geometry, aligned coset, and direction in GF(2^4).  Unit
    # vectors establish equality of each complete linear map.
    for length in (1, 2, 4, 8, 16):
        for shift in range(0, gf2_4.order, length):
            for active in range(length + 1):
                input_mask = tuple(index < active for index in range(length))
                for requested in range(length + 1):
                    output_mask = tuple(index < requested
                                        for index in range(length))
                    for direction in ("forward", "inverse"):
                        _verify_plan(
                            gf2_4, length, shift, input_mask, output_mask,
                            direction, counts
                        )
                        counts.gf2_4_prefix_plans += 1

    # Every sparse input/output mask through N=8 in every aligned GF(2^4)
    # coset.  At N=8 this is 256 x 256 masks x two cosets x two directions.
    # Prefix and sparse tests overlap intentionally: their contracts and
    # result counters are independent.  The independent matrix is exhaustive;
    # the secondary padded-vector check is retained through N=4 because broad
    # padded coverage is also supplied by the prefix and GF8 sweeps.
    for length in (1, 2, 4, 8):
        for shift in range(0, gf2_4.order, length):
            for input_bits in range(1 << length):
                input_mask = tuple(bool(input_bits & (1 << index))
                                   for index in range(length))
                for output_bits in range(1 << length):
                    output_mask = tuple(bool(output_bits & (1 << index))
                                        for index in range(length))
                    for direction in ("forward", "inverse"):
                        _verify_plan(
                            gf2_4, length, shift, input_mask, output_mask,
                            direction, counts, check_padded=(length <= 4)
                        )
                        counts.gf2_4_sparse_plans += 1

    # Broad GF8 shifted-coset coverage around every dyadic boundary.  Direct
    # matrices are cached per (N, shift), so both directions remain independent
    # of the butterfly candidate.
    for length in (2, 4, 8, 16, 32, 64, 128, 256):
        shifts = sorted({0, gf8.order - length,
                         (gf8.order // (2 * length)) * length})
        candidates = sorted({
            0, 1, max(0, length // 2 - 1), length // 2,
            min(length, length // 2 + 1), max(0, length - 1), length,
        })
        for shift in shifts:
            for active in candidates:
                input_mask = tuple(index < active for index in range(length))
                for requested in candidates:
                    output_mask = tuple(index < requested
                                        for index in range(length))
                    for direction in ("forward", "inverse"):
                        _verify_plan(
                            gf8, length, shift, input_mask, output_mask,
                            direction, counts,
                            check_padded=(active in (0, 1, length // 2, length)
                                          and requested in
                                          (0, 1, length // 2, length)),
                        )
                        counts.gf8_prefix_plans += 1

        sparse_inputs = (
            tuple(index % 2 == 0 for index in range(length)),
            tuple(index in (0, length // 2, length - 1)
                  for index in range(length)),
            tuple((index * 5 + 3) % 11 < 4 for index in range(length)),
        )
        sparse_outputs = (
            tuple(index % 3 == 0 for index in range(length)),
            tuple(index in (0, length // 2, length - 1)
                  for index in range(length)),
            tuple((index * 7 + 1) % 13 < 5 for index in range(length)),
        )
        for shift in shifts:
            for input_mask in sparse_inputs:
                for output_mask in sparse_outputs:
                    for direction in ("forward", "inverse"):
                        _verify_plan(
                            gf8, length, shift, input_mask, output_mask,
                            direction, counts, check_padded=True
                        )
                        counts.gf8_sparse_plans += 1

    # Fail-closed compact interface checks.
    plan = compile_truncated_plan(gf8, 8, 3, 5)
    for invalid in ((1, 2), (1, 2, 3, 4), (1, 2, 256)):
        try:
            execute_truncated(plan, invalid)
        except ValueError:
            counts.compact_contract_rejections += 1
        else:
            raise AssertionError("compact interface accepted invalid input")

    representative = [
        _representative_metrics(gf8, length, active, requested, direction)
        for direction in ("forward", "inverse")
        for length, active, requested in (
            (16, 9, 7), (64, 33, 17), (256, 129, 65)
        )
    ]
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "pass",
        "oracle": {
            "forward": "direct normalized-LCH polynomial evaluation matrix",
            "inverse": "GF Gaussian inverse of the direct evaluation matrix",
            "secondary": "explicit full padded Leopard butterfly recurrence",
        },
        "coverage": asdict(counts),
        "operation_count_scope": (
            "butterfly and compact boundary combinations; plan setup, "
            "Python container conversion, and final scatter excluded"
        ),
        "representative_operation_and_memory_counts": representative,
        "timing": "not measured",
    }


def write_json(path: Path | None, payload: object) -> None:
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if path is None:
        sys.stdout.write(encoded)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(encoded, encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args(argv)
    payload = run_self_test()
    write_json(arguments.output, payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
