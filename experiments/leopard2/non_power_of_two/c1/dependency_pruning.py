#!/usr/bin/env python3
"""Wire-preserving dependency pruning for Leopard's dyadic LCH transform.

This is an isolated scalar research model, not a production kernel.  It keeps
the enclosing power-of-two parent, normalized LCH basis, coordinate order, and
butterfly factors unchanged.  Only operations proven irrelevant from known-zero
inputs and requested outputs are removed.

The four execution forms are:

* recursive: a pruned recursive transform tree with runtime node branches;
* flat: a precompiled list containing only required scalar operations;
* hybrid: complete dyadic subtransforms plus scalar boundary operations;
* generated: an experiment-only Python function with the flat constants baked
  into source.  It is never imported by the library or default build.

Run ``self-test`` before ``benchmark``.  JSON output is deterministic except
for the explicitly host-dependent timing and platform sections of benchmark.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import statistics
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence


SCHEMA_VERSION = "leopard2-c1-v1"


class BinaryField:
    """Small table-backed binary extension field in a public coordinate basis."""

    def __init__(
        self, bits: int, polynomial: int, coordinate_basis: Sequence[int]
    ) -> None:
        if bits <= 0 or bits > 8 or len(coordinate_basis) != bits:
            raise ValueError("C1 supports test fields through GF(2^8)")
        if polynomial & (1 << bits) == 0:
            raise ValueError("irreducible polynomial has the wrong degree")
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


def make_gf4() -> BinaryField:
    return BinaryField(4, 0x13, (1, 2, 4, 8))


def make_legacy_gf8() -> BinaryField:
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


def normalized_lch_polynomials(
    field: BinaryField, length: int
) -> tuple[tuple[int, ...], ...]:
    if length <= 0 or length & (length - 1):
        raise ValueError("LCH length must be a power of two")
    dimension = length.bit_length() - 1
    if dimension > field.bits:
        raise ValueError("LCH length exceeds field")
    subspaces = [subspace_polynomial(field, bit)
                 for bit in range(dimension)]
    result: list[tuple[int, ...]] = []
    for index in range(length):
        polynomial = [1]
        normalizer = 1
        for bit in range(dimension):
            if not index & (1 << bit):
                continue
            polynomial = polynomial_multiply(field, polynomial, subspaces[bit])
            factor = polynomial_evaluate(field, subspaces[bit], 1 << bit)
            normalizer = field.multiply(normalizer, factor)
        inverse = field.inverse(normalizer)
        result.append(tuple(field.multiply(value, inverse)
                            for value in polynomial))
    return tuple(result)


def direct_lch_evaluate(
    field: BinaryField,
    coefficients: Sequence[int],
    shift: int,
) -> list[int]:
    """Independent monomial-polynomial evaluation oracle."""
    length = len(coefficients)
    basis = normalized_lch_polynomials(field, length)
    polynomial = [0] * length
    for coefficient, basis_polynomial in zip(coefficients, basis):
        for degree, value in enumerate(basis_polynomial):
            polynomial[degree] ^= field.multiply(coefficient, value)
    return [polynomial_evaluate(field, polynomial, shift ^ output)
            for output in range(length)]


def make_fft_skew(field: BinaryField) -> tuple[int, ...]:
    """Reconstruct Leopard FFTSkew before its logarithm conversion.

    Leopard stores logarithms for its multiply kernels.  The scalar model keeps
    the equivalent field multiplier, which also exposes exact 0/1 factors.
    """
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

        product = field.multiply(temp[layer], temp[layer] ^ 1)
        scale = field.inverse(product)
        for index in range(layer + 1, field.bits - 1):
            temp[index] = field.multiply(
                field.multiply(temp[index], temp[index] ^ 1), scale
            )
    return tuple(skew)


@dataclass(frozen=True)
class RawOperation:
    identifier: int
    x: int
    y: int
    a: int
    b: int
    c: int
    d: int
    inverse: bool
    node_start: int
    node_length: int


@dataclass
class TransformNode:
    start: int
    length: int
    shift: int
    local_ids: tuple[int, ...] = ()
    left: "TransformNode | None" = None
    right: "TransformNode | None" = None
    active_count: int = 0


def make_operation_tree(
    field: BinaryField, length: int, shift: int, direction: str
) -> tuple[TransformNode, tuple[RawOperation, ...]]:
    if length <= 0 or length & (length - 1):
        raise ValueError("transform length must be a power of two")
    if shift < 0 or shift & (length - 1) or shift + length > field.order:
        raise ValueError("shift must name an aligned in-field additive coset")
    if direction not in ("forward", "inverse"):
        raise ValueError("direction must be forward or inverse")
    skew = make_fft_skew(field)
    operations: list[RawOperation] = []

    def build(start: int, size: int, coset: int) -> TransformNode:
        node = TransformNode(start, size, coset)
        if size == 1:
            return node
        half = size // 2
        multiplier = skew[coset + half - 1]
        local: list[RawOperation] = []
        for offset in range(half):
            identifier = len(operations) + len(local)
            if direction == "forward":
                matrix = (1, multiplier, 1, multiplier ^ 1)
            else:
                matrix = (multiplier ^ 1, multiplier, 1, 1)
            local.append(RawOperation(
                identifier, start + offset, start + half + offset,
                *matrix, direction == "inverse", start, size
            ))

        if direction == "forward":
            operations.extend(local)
            node.local_ids = tuple(operation.identifier for operation in local)
            node.left = build(start, half, coset)
            node.right = build(start + half, half, coset + half)
        else:
            node.left = build(start, half, coset)
            node.right = build(start + half, half, coset + half)
            # IDs are assigned only when appended so that identifier == flat index.
            local = []
            for offset in range(half):
                identifier = len(operations) + len(local)
                multiplier = skew[coset + half - 1]
                local.append(RawOperation(
                    identifier, start + offset, start + half + offset,
                    multiplier ^ 1, multiplier, 1, 1, True, start, size
                ))
            operations.extend(local)
            node.local_ids = tuple(operation.identifier for operation in local)
        return node

    root = build(0, length, shift)
    return root, tuple(operations)


@dataclass(frozen=True)
class PlannedOperation:
    raw: RawOperation
    live_x: bool
    live_y: bool
    need_x: bool
    need_y: bool

    @property
    def write_x(self) -> bool:
        if not self.need_x:
            return False
        raw = self.raw
        other_term = self.live_y and raw.b != 0
        return not ((self.live_x and raw.a == 1 and not other_term) or
                    (not self.live_x and not other_term))

    @property
    def write_y(self) -> bool:
        if not self.need_y:
            return False
        raw = self.raw
        other_term = self.live_x and raw.c != 0
        return not ((self.live_y and raw.d == 1 and not other_term) or
                    (not self.live_y and not other_term))

    @property
    def executable(self) -> bool:
        return self.write_x or self.write_y


@dataclass(frozen=True)
class PlanMetrics:
    length: int
    full_butterflies: int
    active_butterflies: int
    skipped_butterflies: int
    one_output_butterflies: int
    input_zero_specializations: int
    multiplier_zero_butterflies: int
    multiplier_one_butterflies: int
    fixed_multiplications_per_byte: int
    xors_per_byte: int
    loads_per_byte: int
    stores_per_byte: int
    padded_fixed_multiplications_per_byte: int
    padded_xors_per_byte: int
    padded_loads_per_byte: int
    padded_stores_per_byte: int
    flat_schedule_bytes: int
    recursive_nodes_visited: int
    recursive_nodes_skipped: int
    hybrid_complete_subtransforms: int
    hybrid_complete_butterflies: int
    hybrid_boundary_butterflies: int
    hybrid_schedule_items: int


@dataclass
class TransformPlan:
    field: BinaryField
    length: int
    shift: int
    direction: str
    input_mask: tuple[bool, ...]
    output_mask: tuple[bool, ...]
    root: TransformNode
    raw_operations: tuple[RawOperation, ...]
    by_identifier: tuple[PlannedOperation, ...]
    flat_operations: tuple[PlannedOperation, ...]
    hybrid_items: tuple[tuple[str, object], ...]
    metrics: PlanMetrics


def _linear_dependencies(
    left_live: bool, right_live: bool, left_coefficient: int,
    right_coefficient: int
) -> bool:
    return ((left_live and left_coefficient != 0) or
            (right_live and right_coefficient != 0))


def _operation_cost(operation: PlannedOperation) -> tuple[int, int, int, int]:
    """Return fixed multiplications, XORs, loads, stores per payload byte."""
    raw = operation.raw
    used_inputs: set[int] = set()
    multiplications = 0
    xors = 0
    stores = 0
    write_x, write_y = operation.write_x, operation.write_y
    # When both destinations change, preserve Leopard's two-step butterfly so
    # the second output reuses the first and only one fixed multiplication is
    # needed.  A single requested output uses its direct matrix row.
    if write_x and write_y:
        multiplier = raw.b
        if operation.live_x and operation.live_y:
            return ((1 if multiplier else 0), (2 if multiplier else 1),
                    2, (2 if multiplier else 1))
        if operation.live_x:
            if raw.inverse:
                factor = multiplier ^ 1
                return ((1 if factor not in (0, 1) else 0), 0, 1, 2)
            return (0, 0, 1, 1)  # x stays in place; y is a copy.
        if operation.live_y:
            if raw.inverse:
                return ((1 if multiplier not in (0, 1) else 0), 0, 1,
                        1 if multiplier else 0)
            return ((1 if multiplier else 0), (1 if multiplier else 0), 1,
                    2 if multiplier else 0)
        return (0, 0, 0, 0)

    for needed, coefficients in (
        (write_x, (raw.a, raw.b)),
        (write_y, (raw.c, raw.d)),
    ):
        if not needed:
            continue
        terms = 0
        for live, coefficient, slot in zip(
            (operation.live_x, operation.live_y), coefficients, (raw.x, raw.y)
        ):
            if live and coefficient:
                terms += 1
                used_inputs.add(slot)
                if coefficient != 1:
                    multiplications += 1
        xors += max(0, terms - 1)
        stores += 1
    return multiplications, xors, len(used_inputs), stores


def compile_plan(
    field: BinaryField,
    length: int,
    active_inputs: int | Sequence[bool],
    requested_outputs: int | Sequence[bool],
    *,
    shift: int = 0,
    direction: str = "forward",
) -> TransformPlan:
    def mask(value: int | Sequence[bool], name: str) -> tuple[bool, ...]:
        if isinstance(value, int):
            if value < 0 or value > length:
                raise ValueError(f"{name} prefix is outside the transform")
            return tuple(index < value for index in range(length))
        result = tuple(bool(item) for item in value)
        if len(result) != length:
            raise ValueError(f"{name} mask has the wrong length")
        return result

    input_mask = mask(active_inputs, "input")
    output_mask = mask(requested_outputs, "output")
    root, operations = make_operation_tree(field, length, shift, direction)

    live = list(input_mask)
    live_before: list[tuple[bool, bool]] = []
    for operation in operations:
        left, right = live[operation.x], live[operation.y]
        live_before.append((left, right))
        live[operation.x] = _linear_dependencies(left, right,
                                                  operation.a, operation.b)
        live[operation.y] = _linear_dependencies(left, right,
                                                  operation.c, operation.d)

    # A requested output that is structurally zero requires no input and no
    # store: the caller-provided/locally-created padded slot already contains
    # zero.  Masking here also prevents a false dependency on an unavailable
    # known-zero input in the all-zero plan.
    needed = [requested and alive for requested, alive in zip(output_mask, live)]
    planned_reverse: list[PlannedOperation] = []
    for operation, before in zip(reversed(operations), reversed(live_before)):
        need_x, need_y = needed[operation.x], needed[operation.y]
        left, right = before
        planned_reverse.append(PlannedOperation(
            operation, left, right, need_x, need_y
        ))
        needed[operation.x] = left and (
            (need_x and operation.a != 0) or (need_y and operation.c != 0)
        )
        needed[operation.y] = right and (
            (need_x and operation.b != 0) or (need_y and operation.d != 0)
        )
    by_identifier = tuple(reversed(planned_reverse))
    if any(required and not present
           for required, present in zip(needed, input_mask)):
        raise AssertionError("pruning analysis requested an unavailable input")
    flat = tuple(operation for operation in by_identifier if operation.executable)

    def mark_counts(node: TransformNode) -> int:
        count = sum(by_identifier[item].executable for item in node.local_ids)
        if node.left:
            count += mark_counts(node.left)
            count += mark_counts(node.right)  # type: ignore[arg-type]
        node.active_count = count
        return count

    mark_counts(root)

    def subtree_ids(node: TransformNode) -> Iterable[int]:
        if direction == "forward":
            yield from node.local_ids
        if node.left:
            yield from subtree_ids(node.left)
            yield from subtree_ids(node.right)  # type: ignore[arg-type]
        if direction == "inverse":
            yield from node.local_ids

    hybrid: list[tuple[str, object]] = []

    def complete(node: TransformNode) -> bool:
        identifiers = tuple(subtree_ids(node))
        return bool(identifiers) and all(
            by_identifier[item].executable and
            by_identifier[item].need_x and by_identifier[item].need_y and
            by_identifier[item].live_x and by_identifier[item].live_y
            for item in identifiers
        )

    def emit_hybrid(node: TransformNode) -> None:
        if node.active_count == 0:
            return
        if node.length >= 4 and complete(node):
            hybrid.append(("complete", tuple(
                operations[item] for item in subtree_ids(node)
            )))
            return
        if direction == "forward":
            hybrid.extend(("boundary", by_identifier[item])
                          for item in node.local_ids
                          if by_identifier[item].executable)
        if node.left:
            emit_hybrid(node.left)
            emit_hybrid(node.right)  # type: ignore[arg-type]
        if direction == "inverse":
            hybrid.extend(("boundary", by_identifier[item])
                          for item in node.local_ids
                          if by_identifier[item].executable)

    emit_hybrid(root)

    multiplications = xors = loads = stores = 0
    for operation in flat:
        values = _operation_cost(operation)
        multiplications += values[0]
        xors += values[1]
        loads += values[2]
        stores += values[3]

    padded_multiplications = padded_xors = padded_loads = padded_stores = 0
    for raw in operations:
        values = _operation_cost(PlannedOperation(raw, True, True, True, True))
        padded_multiplications += values[0]
        padded_xors += values[1]
        padded_loads += values[2]
        padded_stores += values[3]

    def recursive_counts(node: TransformNode) -> tuple[int, int]:
        if node.active_count == 0:
            return (0, 1)
        visited, skipped = 1, 0
        if node.left:
            left_counts = recursive_counts(node.left)
            right_counts = recursive_counts(node.right)  # type: ignore[arg-type]
            visited += left_counts[0] + right_counts[0]
            skipped += left_counts[1] + right_counts[1]
        return visited, skipped

    recursive_visited, recursive_skipped = recursive_counts(root)
    metrics = PlanMetrics(
        length=length,
        full_butterflies=len(operations),
        active_butterflies=len(flat),
        skipped_butterflies=len(operations) - len(flat),
        one_output_butterflies=sum(operation.need_x != operation.need_y
                                   for operation in flat),
        input_zero_specializations=sum(operation.live_x != operation.live_y
                                       for operation in flat),
        multiplier_zero_butterflies=sum(operation.raw.b == 0
                                        for operation in flat),
        multiplier_one_butterflies=sum(operation.raw.b == 1
                                       for operation in flat),
        fixed_multiplications_per_byte=multiplications,
        xors_per_byte=xors,
        loads_per_byte=loads,
        stores_per_byte=stores,
        padded_fixed_multiplications_per_byte=padded_multiplications,
        padded_xors_per_byte=padded_xors,
        padded_loads_per_byte=padded_loads,
        padded_stores_per_byte=padded_stores,
        # Stable C-like estimate: 2x uint16 indices, four uint8 factors/flags,
        # and one uint32 node tag rounded to 16 bytes.
        flat_schedule_bytes=16 * len(flat),
        recursive_nodes_visited=recursive_visited,
        recursive_nodes_skipped=recursive_skipped,
        hybrid_complete_subtransforms=sum(kind == "complete"
                                          for kind, _ in hybrid),
        hybrid_complete_butterflies=sum(len(item) for kind, item in hybrid
                                        if kind == "complete"),
        hybrid_boundary_butterflies=sum(kind == "boundary"
                                        for kind, _ in hybrid),
        hybrid_schedule_items=len(hybrid),
    )
    return TransformPlan(
        field, length, shift, direction, input_mask, output_mask, root,
        operations, by_identifier, flat, tuple(hybrid), metrics
    )


def _linear_value(
    field: BinaryField, left: int, right: int, a: int, b: int,
    live_left: bool = True, live_right: bool = True
) -> int:
    result = 0
    if live_left and a:
        result ^= left if a == 1 else field.mul_rows[a][left]
    if live_right and b:
        result ^= right if b == 1 else field.mul_rows[b][right]
    return result


def _apply_symbol_operation(
    field: BinaryField, work: list[int], operation: PlannedOperation
) -> None:
    raw = operation.raw
    left, right = work[raw.x], work[raw.y]
    write_x, write_y = operation.write_x, operation.write_y
    multiplier = raw.b
    if write_x and write_y:
        if raw.inverse and operation.live_x and not operation.live_y:
            factor = multiplier ^ 1
            new_right = left
            new_left = (0 if factor == 0 else left if factor == 1 else
                        field.mul_rows[factor][left])
        elif raw.inverse:
            new_right = right ^ left
            new_left = left ^ (0 if multiplier == 0 else
                               field.mul_rows[multiplier][new_right])
        else:
            new_left = left ^ (0 if multiplier == 0 else
                               field.mul_rows[multiplier][right])
            new_right = right ^ new_left
        work[raw.x], work[raw.y] = new_left, new_right
        return
    if write_x:
        work[raw.x] = _linear_value(
            field, left, right, raw.a, raw.b,
            operation.live_x, operation.live_y
        )
    if write_y:
        work[raw.y] = _linear_value(
            field, left, right, raw.c, raw.d,
            operation.live_x, operation.live_y
        )


def _check_input(plan: TransformPlan, values: Sequence[int]) -> None:
    if len(values) != plan.length:
        raise ValueError("input has the wrong transform length")
    if any(value and not live for value, live in zip(values, plan.input_mask)):
        raise ValueError("input violates the plan's known-zero mask")


def execute_flat(plan: TransformPlan, values: Sequence[int]) -> list[int]:
    _check_input(plan, values)
    work = list(values)
    for operation in plan.flat_operations:
        _apply_symbol_operation(plan.field, work, operation)
    return work


def execute_recursive(plan: TransformPlan, values: Sequence[int]) -> list[int]:
    _check_input(plan, values)
    work = list(values)

    def visit(node: TransformNode) -> None:
        if node.active_count == 0:
            return
        if plan.direction == "forward":
            for identifier in node.local_ids:
                operation = plan.by_identifier[identifier]
                if operation.executable:
                    _apply_symbol_operation(plan.field, work, operation)
        if node.left:
            visit(node.left)
            visit(node.right)  # type: ignore[arg-type]
        if plan.direction == "inverse":
            for identifier in node.local_ids:
                operation = plan.by_identifier[identifier]
                if operation.executable:
                    _apply_symbol_operation(plan.field, work, operation)
    visit(plan.root)
    return work


def execute_hybrid(plan: TransformPlan, values: Sequence[int]) -> list[int]:
    _check_input(plan, values)
    work = list(values)
    for kind, item in plan.hybrid_items:
        if kind == "boundary":
            _apply_symbol_operation(plan.field, work, item)  # type: ignore[arg-type]
        else:
            for raw in item:  # type: ignore[union-attr]
                operation = PlannedOperation(raw, True, True, True, True)
                _apply_symbol_operation(plan.field, work, operation)
    return work


def execute_full(plan: TransformPlan, values: Sequence[int]) -> list[int]:
    if len(values) != plan.length:
        raise ValueError("input has the wrong transform length")
    work = list(values)
    for raw in plan.raw_operations:
        left, right = work[raw.x], work[raw.y]
        work[raw.x] = _linear_value(plan.field, left, right, raw.a, raw.b)
        work[raw.y] = _linear_value(plan.field, left, right, raw.c, raw.d)
    return work


def compile_generated(
    plan: TransformPlan, *, validate_input: bool = True
) -> object:
    """Compile an experiment-only scalar kernel with embedded constants.

    Validation is on by default, matching the other executors.  The benchmark
    may disable it only after constructing and checking its own valid input so
    it can isolate schedule-interpreter overhead from common validation cost.
    """
    lines = [
        "def generated(values, rows=_rows):",
        "    if len(values) != %d:" % plan.length,
        "        raise ValueError('input has the wrong transform length')",
        "    w = list(values)",
    ]

    def term(slot: int, coefficient: int, live: bool) -> str:
        if not live or coefficient == 0:
            return "0"
        if coefficient == 1:
            return f"v{slot}"
        return f"rows[{coefficient}][v{slot}]"

    for index, operation in enumerate(plan.flat_operations):
        raw = operation.raw
        lines.append(f"    v{raw.x}, v{raw.y} = w[{raw.x}], w[{raw.y}]")
        write_x, write_y = operation.write_x, operation.write_y
        if write_x and write_y:
            multiplier_term = term(raw.y, raw.b, operation.live_y)
            if raw.inverse and operation.live_x and not operation.live_y:
                factor = raw.b ^ 1
                factor_term = ("0" if factor == 0 else
                               f"v{raw.x}" if factor == 1 else
                               f"rows[{factor}][v{raw.x}]")
                lines.append(f"    n{index}y = v{raw.x}")
                lines.append(f"    n{index}x = {factor_term}")
            elif raw.inverse:
                lines.append(f"    n{index}y = v{raw.y} ^ v{raw.x}")
                multiplier_term = ("0" if raw.b == 0 else
                                   (f"n{index}y" if raw.b == 1 else
                                    f"rows[{raw.b}][n{index}y]"))
                lines.append(f"    n{index}x = v{raw.x} ^ {multiplier_term}")
            else:
                lines.append(f"    n{index}x = v{raw.x} ^ {multiplier_term}")
                lines.append(f"    n{index}y = v{raw.y} ^ n{index}x")
            lines.append(f"    w[{raw.x}], w[{raw.y}] = n{index}x, n{index}y")
            continue
        if write_x:
            lines.append(
                f"    w[{raw.x}] = "
                f"{term(raw.x, raw.a, operation.live_x)} ^ "
                f"{term(raw.y, raw.b, operation.live_y)}"
            )
        if write_y:
            lines.append(
                f"    w[{raw.y}] = "
                f"{term(raw.x, raw.c, operation.live_x)} ^ "
                f"{term(raw.y, raw.d, operation.live_y)}"
            )
    lines.append("    return w")
    namespace: dict[str, object] = {"_rows": plan.field.mul_rows}
    source = "\n".join(lines) + "\n"
    exec(compile(source, "<leopard2-c1-generated>", "exec"), namespace)
    function = namespace["generated"]
    source_sha256 = hashlib.sha256(source.encode()).hexdigest()
    source_bytes = len(source.encode())
    setattr(function, "source_sha256", source_sha256)
    setattr(function, "source_bytes", source_bytes)
    if not validate_input:
        return function

    def validated(values: Sequence[int]) -> list[int]:
        _check_input(plan, values)
        return function(values)  # type: ignore[operator]

    setattr(validated, "source_sha256", source_sha256)
    setattr(validated, "source_bytes", source_bytes)
    return validated


def _apply_byte_operation(
    field: BinaryField, work: list[bytearray], operation: PlannedOperation
) -> None:
    raw = operation.raw
    left, right = work[raw.x], work[raw.y]
    rows = field.mul_rows
    write_x, write_y = operation.write_x, operation.write_y
    multiplier = raw.b
    byte_count = len(left)
    if write_x and write_y:
        if raw.inverse and operation.live_x and not operation.live_y:
            factor = multiplier ^ 1
            if factor == 0:
                for index in range(byte_count):
                    right[index], left[index] = left[index], 0
            elif factor == 1:
                for index in range(byte_count):
                    right[index] = left[index]
            else:
                row = rows[factor]
                for index in range(byte_count):
                    value = left[index]
                    right[index], left[index] = value, row[value]
        elif raw.inverse:
            row = rows[multiplier]
            for index in range(byte_count):
                x_value = left[index]
                new_right = right[index] ^ x_value
                right[index] = new_right
                left[index] = x_value ^ row[new_right]
        else:
            row = rows[multiplier]
            for index in range(byte_count):
                y_value = right[index]
                new_left = left[index] ^ row[y_value]
                left[index] = new_left
                right[index] = y_value ^ new_left
        return

    def write_linear(
        destination: bytearray, left_coefficient: int, right_coefficient: int
    ) -> None:
        left_coefficient = left_coefficient if operation.live_x else 0
        right_coefficient = right_coefficient if operation.live_y else 0
        if left_coefficient == 0:
            row = rows[right_coefficient]
            if right_coefficient == 1:
                destination[:] = right
            else:
                for index in range(byte_count):
                    destination[index] = row[right[index]]
        elif right_coefficient == 0:
            row = rows[left_coefficient]
            if left_coefficient == 1:
                destination[:] = left
            else:
                for index in range(byte_count):
                    destination[index] = row[left[index]]
        else:
            left_row = rows[left_coefficient]
            right_row = rows[right_coefficient]
            if left_coefficient == 1 and right_coefficient == 1:
                for index in range(byte_count):
                    destination[index] = left[index] ^ right[index]
            elif left_coefficient == 1:
                for index in range(byte_count):
                    destination[index] = left[index] ^ right_row[right[index]]
            elif right_coefficient == 1:
                for index in range(byte_count):
                    destination[index] = left_row[left[index]] ^ right[index]
            else:
                for index in range(byte_count):
                    destination[index] = (left_row[left[index]] ^
                                          right_row[right[index]])

    if write_x:
        write_linear(left, raw.a, raw.b)
    elif write_y:
        write_linear(right, raw.c, raw.d)


def execute_bytes(plan: TransformPlan, shards: Sequence[bytes], form: str) -> list[bytes]:
    if len(shards) != plan.length or not shards:
        raise ValueError("byte shards have the wrong transform length")
    byte_count = len(shards[0])
    if any(len(shard) != byte_count for shard in shards):
        raise ValueError("byte shards have inconsistent sizes")
    if any(any(shard) and not live for shard, live in zip(shards, plan.input_mask)):
        raise ValueError("byte input violates the plan's known-zero mask")
    work = [bytearray(shard) for shard in shards]

    if form == "padded":
        for raw in plan.raw_operations:
            _apply_byte_operation(plan.field, work,
                                  PlannedOperation(raw, True, True, True, True))
    elif form == "flat":
        for operation in plan.flat_operations:
            _apply_byte_operation(plan.field, work, operation)
    elif form == "recursive":
        def visit(node: TransformNode) -> None:
            if node.active_count == 0:
                return
            if plan.direction == "forward":
                for identifier in node.local_ids:
                    operation = plan.by_identifier[identifier]
                    if operation.executable:
                        _apply_byte_operation(plan.field, work, operation)
            if node.left:
                visit(node.left)
                visit(node.right)  # type: ignore[arg-type]
            if plan.direction == "inverse":
                for identifier in node.local_ids:
                    operation = plan.by_identifier[identifier]
                    if operation.executable:
                        _apply_byte_operation(plan.field, work, operation)
        visit(plan.root)
    elif form == "hybrid":
        for kind, item in plan.hybrid_items:
            if kind == "boundary":
                _apply_byte_operation(plan.field, work, item)  # type: ignore[arg-type]
            else:
                for raw in item:  # type: ignore[union-attr]
                    _apply_byte_operation(plan.field, work,
                                          PlannedOperation(raw, True, True, True, True))
    else:
        raise ValueError("unknown byte execution form")
    return [bytes(shard) for shard in work]


def _assert_requested_equal(
    plan: TransformPlan, actual: Sequence[int], expected: Sequence[int]
) -> None:
    for index, requested in enumerate(plan.output_mask):
        if requested and actual[index] != expected[index]:
            raise AssertionError(
                f"{plan.direction} output {index} differs: "
                f"{actual[index]} != {expected[index]}"
            )


@dataclass
class TestCounts:
    direct_oracle_vectors: int = 0
    inverse_vectors: int = 0
    gf4_prefix_plans: int = 0
    gf4_sparse_plans: int = 0
    gf4_sparse_input_plans: int = 0
    gf8_swept_plans: int = 0
    gf8_sparse_input_plans: int = 0
    execution_form_comparisons: int = 0
    requested_symbol_comparisons: int = 0
    generated_kernel_comparisons: int = 0
    generated_mask_rejections: int = 0


def run_self_test() -> dict[str, object]:
    counts = TestCounts()
    gf4 = make_gf4()
    gf8 = make_legacy_gf8()
    rng = random.Random(0xC1D3E3D)

    # Direct polynomial evaluation is independent of the skew/butterfly model.
    for field in (gf4, gf8):
        maximum = 16 if field.bits == 4 else 256
        lengths = (1, 2, 4, 8, 16) if field.bits == 4 else (2, 4, 8, 16, 32, 64, 128, 256)
        for length in lengths:
            shifts = range(0, maximum, length) if field.bits == 4 else (
                (0, maximum - length) if length < maximum else (0,)
            )
            for shift in shifts:
                full = compile_plan(field, length, length, length, shift=shift)
                inverse = compile_plan(field, length, length, length,
                                       shift=shift, direction="inverse")
                vectors: list[list[int]] = []
                for index in range(length):
                    vector = [0] * length
                    vector[index] = 1
                    vectors.append(vector)
                if length <= 2:
                    if length == 1:
                        vectors.extend([[value] for value in range(field.order)])
                    else:
                        vectors.extend([[left, right]
                                        for left in range(field.order)
                                        for right in range(field.order)])
                vectors.extend([[rng.randrange(field.order) for _ in range(length)]
                                for _ in range(4)])
                for vector in vectors:
                    direct = direct_lch_evaluate(field, vector, shift)
                    transformed = execute_full(full, vector)
                    if transformed != direct:
                        raise AssertionError(
                            f"direct LCH oracle differs for GF{field.order}, "
                            f"N={length}, shift={shift}"
                        )
                    counts.direct_oracle_vectors += 1
                    recovered = execute_full(inverse, transformed)
                    if recovered != vector:
                        raise AssertionError("inverse transform is not an inverse")
                    counts.inverse_vectors += 1

    def verify_plan(
        plan: TransformPlan, values: list[int], expected: list[int],
        *, generated: bool = False
    ) -> None:
        forms: tuple[list[int], ...] = (
            execute_recursive(plan, values),
            execute_flat(plan, values),
            execute_hybrid(plan, values),
        )
        if generated:
            forms += (compile_generated(plan)(values),)
            counts.generated_kernel_comparisons += 1
        for actual in forms:
            _assert_requested_equal(plan, actual, expected)
            counts.execution_form_comparisons += 1
            counts.requested_symbol_comparisons += sum(plan.output_mask)

    # Every GF4 prefix geometry in every supported active parent, both ways.
    for length in (1, 2, 4, 8, 16):
        for active in range(length + 1):
            values = [rng.randrange(1, gf4.order) if index < active else 0
                      for index in range(length)]
            for requested in range(length + 1):
                for direction in ("forward", "inverse"):
                    plan = compile_plan(gf4, length, active, requested,
                                        direction=direction)
                    expected = execute_full(plan, values)
                    verify_plan(plan, values, expected,
                                generated=(length <= 4 and active == requested))
                    counts.gf4_prefix_plans += 1

    # Exhaust every sparse GF4 output mask through N=8 for every input prefix.
    for length in (1, 2, 4, 8):
        for active in range(length + 1):
            values = [rng.randrange(1, gf4.order) if index < active else 0
                      for index in range(length)]
            for bits in range(1 << length):
                requested = tuple(bool(bits & (1 << index))
                                  for index in range(length))
                plan = compile_plan(gf4, length, active, requested)
                verify_plan(plan, values, execute_full(plan, values))
                counts.gf4_sparse_plans += 1

    # Exhaust arbitrary input and output masks through N=4, then sweep all N=8
    # input masks against representative irregular outputs.  This catches DAG
    # mistakes hidden by the encoder's usual known-nonzero prefix convention.
    for length in (1, 2, 4):
        for input_bits in range(1 << length):
            active_mask = tuple(bool(input_bits & (1 << index))
                                for index in range(length))
            values = [rng.randrange(1, gf4.order) if active else 0
                      for active in active_mask]
            for output_bits in range(1 << length):
                requested_mask = tuple(bool(output_bits & (1 << index))
                                       for index in range(length))
                for direction in ("forward", "inverse"):
                    plan = compile_plan(gf4, length, active_mask,
                                        requested_mask, direction=direction)
                    verify_plan(plan, values, execute_full(plan, values))
                    counts.gf4_sparse_input_plans += 1
    length = 8
    output_masks = [
        tuple(index < prefix for index in range(length))
        for prefix in (0, 1, 3, 4, 5, 7, 8)
    ] + [
        tuple(index % 3 == residue for index in range(length))
        for residue in range(3)
    ]
    for input_bits in range(1 << length):
        active_mask = tuple(bool(input_bits & (1 << index))
                            for index in range(length))
        values = [rng.randrange(1, gf4.order) if active else 0
                  for active in active_mask]
        for requested_mask in output_masks:
            for direction in ("forward", "inverse"):
                plan = compile_plan(gf4, length, active_mask, requested_mask,
                                    direction=direction)
                verify_plan(plan, values, execute_full(plan, values))
                counts.gf4_sparse_input_plans += 1

    # GF8 boundaries just below/at/above dyadic sizes, prefixes and sparse masks.
    for length in (2, 4, 8, 16, 32, 64, 128, 256):
        candidates = sorted({0, 1, max(0, length // 2 - 1), length // 2,
                             min(length, length // 2 + 1), max(0, length - 1), length})
        masks = [tuple(index < requested for index in range(length))
                 for requested in candidates]
        masks.extend(
            tuple(index == 0 or index + 1 == length or index % 7 == residue
                  for index in range(length))
            for residue in (0, 3)
        )
        for active in candidates:
            values = [rng.randrange(1, gf8.order) if index < active else 0
                      for index in range(length)]
            for requested in masks:
                for direction in ("forward", "inverse"):
                    plan = compile_plan(gf8, length, active, requested,
                                        direction=direction)
                    generated = (
                        direction == "forward" and
                        length in (16, 64, 256) and
                        active == length // 2 + 1 and
                        requested == masks[3]
                    )
                    verify_plan(plan, values, execute_full(plan, values),
                                generated=generated)
                    counts.gf8_swept_plans += 1

        sparse_inputs = (
            tuple(index % 2 == 0 for index in range(length)),
            tuple(index in (0, length // 2, length - 1) for index in range(length)),
            tuple((index * 5 + 3) % 11 < 4 for index in range(length)),
        )
        sparse_outputs = (masks[1], masks[-2], masks[-1])
        for active_mask in sparse_inputs:
            values = [rng.randrange(1, gf8.order) if active else 0
                      for active in active_mask]
            for requested_mask in sparse_outputs:
                for direction in ("forward", "inverse"):
                    plan = compile_plan(gf8, length, active_mask,
                                        requested_mask, direction=direction)
                    verify_plan(plan, values, execute_full(plan, values))
                    counts.gf8_sparse_input_plans += 1

    # Byte execution exercises tails and validates the three runtime forms.
    byte_cases = 0
    byte_comparisons = 0
    for length, active, requested, byte_count in (
        (8, 5, 3, 1), (16, 9, 7, 7), (32, 17, 19, 31),
        (64, 33, 17, 65), (128, 65, 63, 257), (256, 129, 65, 1025),
    ):
        plan = compile_plan(gf8, length, active, requested)
        shards = [bytes(rng.randrange(256) for _ in range(byte_count))
                  if index < active else bytes(byte_count)
                  for index in range(length)]
        # Full byte oracle is evaluated independently one stripe at a time.
        expected = [[0] * byte_count for _ in range(length)]
        for byte in range(byte_count):
            vector = [shard[byte] for shard in shards]
            output = execute_full(plan, vector)
            for index, value in enumerate(output):
                expected[index][byte] = value
        for form in ("recursive", "flat", "hybrid"):
            actual = execute_bytes(plan, shards, form)
            for index, requested_output in enumerate(plan.output_mask):
                if requested_output and actual[index] != bytes(expected[index]):
                    raise AssertionError(f"{form} byte execution differs")
                if requested_output:
                    byte_comparisons += byte_count
        byte_cases += 1

    generated_plan = compile_plan(gf8, 8, 3, 8)
    generated = compile_generated(generated_plan)
    invalid = [1, 2, 3, 0, 0, 0, 0, 1]
    try:
        generated(invalid)  # type: ignore[operator]
    except ValueError:
        counts.generated_mask_rejections += 1
    else:
        raise AssertionError("generated kernel accepted a known-zero violation")

    payload = {
        "schema_version": SCHEMA_VERSION,
        "seed": "0x0c1d3e3d",
        "counts": asdict(counts),
        "byte_cases": byte_cases,
        "byte_comparisons": byte_comparisons,
        "status": "pass",
    }
    return payload


def _benchmark_case(
    field: BinaryField,
    length: int,
    active: int,
    requested: int,
    byte_count: int,
    seed: int,
) -> dict[str, object]:
    plan = compile_plan(field, length, active, requested)
    rng = random.Random(seed)
    shards = [bytes(rng.randrange(256) for _ in range(byte_count))
              if index < active else bytes(byte_count)
              for index in range(length)]
    reference = execute_bytes(plan, shards, "padded")
    repetitions = 9 if byte_count <= 1024 else 7
    forms = ("padded", "recursive", "flat", "hybrid")
    samples_by_form: dict[str, list[int]] = {form: [] for form in forms}
    # Warm interpreter specialization, page mappings, and input copies for each
    # form.  Then rotate execution order so frequency/thermal drift cannot
    # consistently favor one schedule representation.
    for form in forms:
        for _ in range(2):
            execute_bytes(plan, shards, form)
    for repetition in range(repetitions):
        ordered = forms[repetition % len(forms):] + forms[:repetition % len(forms)]
        for form in ordered:
            started = time.perf_counter_ns()
            actual = execute_bytes(plan, shards, form)
            samples_by_form[form].append(time.perf_counter_ns() - started)
            if any(actual[index] != reference[index]
                   for index, needed in enumerate(plan.output_mask) if needed):
                raise AssertionError(f"benchmark {form} differs from flat")
    timings: dict[str, dict[str, object]] = {}
    for form in forms:
        samples = samples_by_form[form]
        timings[form] = {
            "median_ns": int(statistics.median(samples)),
            "mad_ns": int(statistics.median(
                abs(sample - statistics.median(samples)) for sample in samples
            )),
            "repetitions": repetitions,
        }
    # Inputs were constructed from this plan's active prefix above.  Time the
    # prevalidated generated body so this micro-row measures schedule dispatch,
    # while the default public experiment callable remains fail-closed.
    generated = compile_generated(plan, validate_input=False)
    scalar = [shard[0] if shard else 0 for shard in shards]
    generated_samples: list[int] = []
    flat_scalar_samples: list[int] = []
    generated_repetitions = max(101, 100000 // max(1, len(plan.flat_operations)))
    expected_scalar = execute_flat(plan, scalar)
    for _ in range(100):
        generated(scalar)
        execute_flat(plan, scalar)
    for repetition in range(generated_repetitions):
        if repetition & 1:
            started = time.perf_counter_ns()
            actual_flat_scalar = execute_flat(plan, scalar)
            flat_scalar_samples.append(time.perf_counter_ns() - started)
        started = time.perf_counter_ns()
        actual_scalar = generated(scalar)
        generated_samples.append(time.perf_counter_ns() - started)
        if not repetition & 1:
            started = time.perf_counter_ns()
            actual_flat_scalar = execute_flat(plan, scalar)
            flat_scalar_samples.append(time.perf_counter_ns() - started)
        _assert_requested_equal(plan, actual_scalar, expected_scalar)
        _assert_requested_equal(plan, actual_flat_scalar, expected_scalar)
    timings["generated_scalar"] = {
        "median_ns": int(statistics.median(generated_samples)),
        "mad_ns": int(statistics.median(
            abs(sample - statistics.median(generated_samples))
            for sample in generated_samples
        )),
        "repetitions": generated_repetitions,
        "source_bytes": getattr(generated, "source_bytes"),
        "source_sha256": getattr(generated, "source_sha256"),
    }
    timings["flat_scalar"] = {
        "median_ns": int(statistics.median(flat_scalar_samples)),
        "mad_ns": int(statistics.median(
            abs(sample - statistics.median(flat_scalar_samples))
            for sample in flat_scalar_samples
        )),
        "repetitions": generated_repetitions,
    }
    return {
        "direction": "forward",
        "length": length,
        "active_inputs": active,
        "requested_outputs": requested,
        "shard_bytes": byte_count,
        "metrics": asdict(plan.metrics),
        "timings": timings,
    }


def run_benchmark() -> dict[str, object]:
    field = make_legacy_gf8()
    cases: list[tuple[int, int, int, int]] = []
    # The Python scalar model covers tiny through 64 KiB.  Larger N uses a
    # smaller payload cap so the validation remains suitable for CI/laptops.
    for length, active, requested, sizes in (
        (16, 9, 7, (64, 1024, 16384, 65536)),
        (64, 33, 17, (64, 1024, 16384)),
        (256, 129, 65, (64, 1024)),
    ):
        cases.extend((length, active, requested, byte_count)
                     for byte_count in sizes)
    results = [
        _benchmark_case(field, *case, 0xC1B30000 + index)
        for index, case in enumerate(cases)
    ]
    cpu_model = "unknown"
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("model name"):
                cpu_model = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "pass",
        "clock": "time.perf_counter_ns",
        "timing_scope": "Python scalar experiment; setup excluded",
        "seed": "0x0c1b30000 + case_index",
        "host": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "cpu_model": cpu_model,
            "online_cpus": os.cpu_count(),
            "allowed_cpus": len(os.sched_getaffinity(0))
                if hasattr(os, "sched_getaffinity") else os.cpu_count(),
        },
        "cases": results,
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
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("self-test", "benchmark"):
        child = subparsers.add_parser(command)
        child.add_argument("--output", type=Path)
    arguments = parser.parse_args(argv)
    if arguments.command == "self-test":
        write_json(arguments.output, run_self_test())
    else:
        write_json(arguments.output, run_benchmark())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
