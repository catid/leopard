#!/usr/bin/env python3
"""C3 truncated-inverse and basis-conversion scalar experiment.

This file is deliberately isolated from the production library.  It studies
two different mathematical questions and never conflates them:

* ``parent_*`` methods invert the same dyadic parent evaluation vector as
  Leopard.  The unprovided suffix is fixed to zero, so these methods must match
  the full padded inverse byte-for-byte.
* ``exact_*`` methods interpolate the unique degree-<q polynomial through
  exactly q prefix points.  For non-power-of-two q this is a different code
  definition and would require a separately versioned wire profile.

The dense parent-column construction is a correctness model for a
Lagrange-to-LCH basis conversion, not an implementation of Coxon's fast
recursive conversion.  The Newton path is likewise a transparent scalar
conversion model, not a transcription of Tang-Han's epsilon inverse.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import functools
import hashlib
import json
import os
import platform
import random
import statistics
import sys
import time
from pathlib import Path
from typing import Iterable, Sequence


SCHEMA = "leopard2-c3-truncated-inverse/v1"
TIMING_SCHEMA = "leopard2-c3-truncated-inverse-timing/v1"
METHODS = (
    "parent_padded_lch",
    "parent_dense_lagrange_to_lch",
    "exact_dense_lch",
    "exact_newton_to_lch",
    "exact_generic_lagrange",
)


@dataclasses.dataclass
class Counter:
    xors: int = 0
    multiplications: int = 0
    inversions: int = 0
    loads: int = 0
    stores: int = 0

    @property
    def arithmetic(self) -> int:
        return self.xors + self.multiplications

    def add(self, other: "Counter") -> None:
        self.xors += other.xors
        self.multiplications += other.multiplications
        self.inversions += other.inversions
        self.loads += other.loads
        self.stores += other.stores


class BinaryField:
    """GF(2^m) in an explicit public coordinate basis.

    GF4/GF8 use an independently generated multiplication table.  GF16 keeps
    only the two basis-conversion tables and performs carryless multiplication
    directly, avoiding a prohibitively large product table.
    """

    def __init__(
        self, name: str, bits: int, polynomial: int,
        coordinate_basis: Sequence[int],
    ) -> None:
        if bits <= 0 or bits > 16 or len(coordinate_basis) != bits:
            raise ValueError("unsupported binary field")
        if polynomial & (1 << bits) == 0:
            raise ValueError("irreducible polynomial has the wrong degree")
        self.name = name
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

        if self.order <= 256:
            self.mul_rows: tuple[bytes, ...] | None = tuple(
                bytes(self._multiply_slow(left, right)
                      for right in range(self.order))
                for left in range(self.order)
            )
        else:
            self.mul_rows = None

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
        if not (0 <= left < self.order and 0 <= right < self.order):
            raise ValueError("field element out of range")
        if self.mul_rows is not None:
            return self.mul_rows[left][right]
        return self._multiply_slow(left, right)

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


@functools.lru_cache(maxsize=None)
def field_named(name: str) -> BinaryField:
    if name == "gf4":
        return BinaryField("gf4", 4, 0x13, (1, 2, 4, 8))
    if name == "gf8":
        return BinaryField(
            "gf8", 8, 0x11D, (1, 214, 152, 146, 86, 200, 88, 230)
        )
    if name == "gf16":
        return BinaryField(
            "gf16", 16, 0x1002D,
            (
                0x0001, 0xACCA, 0x3C0E, 0x163E,
                0xC582, 0xED2E, 0x914C, 0x4012,
                0x6C98, 0x10D8, 0x6A72, 0xB900,
                0xFDB8, 0xFB34, 0xFF38, 0x991E,
            ),
        )
    raise ValueError("unknown field")


def ceil_power_of_two(value: int) -> int:
    if value <= 0:
        raise ValueError("length must be positive")
    return 1 << (value - 1).bit_length()


def validate_geometry(field: BinaryField, parent: int, shift: int) -> None:
    if parent <= 0 or parent & (parent - 1):
        raise ValueError("parent must be a power of two")
    if shift < 0 or shift & (parent - 1) or shift + parent > field.order:
        raise ValueError("shift must name an aligned in-field additive coset")


Polynomial = list[int]


def trim(polynomial: Polynomial) -> Polynomial:
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def polynomial_add(left: Sequence[int], right: Sequence[int]) -> Polynomial:
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] ^= value
    for index, value in enumerate(right):
        result[index] ^= value
    return trim(result)


def polynomial_multiply(
    field: BinaryField, left: Sequence[int], right: Sequence[int],
    counter: Counter | None = None,
) -> Polynomial:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if counter is not None:
                counter.loads += 2
                counter.stores += 1
                if a != 0 and b not in (0, 1):
                    counter.multiplications += 1
                if result[i + j] != 0 and a != 0 and b != 0:
                    counter.xors += 1
            if a and b:
                result[i + j] ^= a if b == 1 else field.multiply(a, b)
    return trim(result)


def polynomial_scale_add(
    field: BinaryField, destination: Polynomial, source: Sequence[int],
    scalar: int, counter: Counter | None = None,
) -> None:
    if len(destination) < len(source):
        destination.extend([0] * (len(source) - len(destination)))
    for index, coefficient in enumerate(source):
        if counter is not None:
            counter.loads += 2
            counter.stores += 1
            if scalar not in (0, 1) and coefficient != 0:
                counter.multiplications += 1
            if destination[index] != 0 and scalar != 0 and coefficient != 0:
                counter.xors += 1
        if scalar and coefficient:
            destination[index] ^= (
                coefficient if scalar == 1
                else field.multiply(scalar, coefficient)
            )


def polynomial_evaluate(
    field: BinaryField, polynomial: Sequence[int], point: int,
    counter: Counter | None = None,
) -> int:
    value = 0
    for coefficient in reversed(polynomial):
        if counter is not None:
            counter.loads += 1
            counter.stores += 1
            if point not in (0, 1) and value != 0:
                counter.multiplications += 1
            if value != 0 and coefficient != 0:
                counter.xors += 1
        value = (value if point == 1 else field.multiply(value, point)) ^ coefficient
    return value


@functools.lru_cache(maxsize=None)
def subspace_polynomial(field_name: str, dimension: int) -> tuple[int, ...]:
    field = field_named(field_name)
    result: Polynomial = [1]
    for point in range(1 << dimension):
        result = polynomial_multiply(field, result, (point, 1))
    return tuple(result)


@functools.lru_cache(maxsize=None)
def normalized_lch_polynomials(
    field_name: str, parent: int,
) -> tuple[tuple[int, ...], ...]:
    field = field_named(field_name)
    if parent <= 0 or parent & (parent - 1) or parent > field.order:
        raise ValueError("invalid LCH parent")
    dimension = parent.bit_length() - 1
    subspaces = tuple(
        subspace_polynomial(field_name, bit) for bit in range(dimension)
    )
    basis: list[tuple[int, ...]] = []
    for index in range(parent):
        polynomial: Polynomial = [1]
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
        basis.append(tuple(field.multiply(value, inverse)
                           for value in polynomial))
    return tuple(basis)


def lch_evaluation_matrix(
    field: BinaryField, count: int, shift: int, parent: int,
    counter: Counter | None = None,
) -> tuple[tuple[int, ...], ...]:
    validate_geometry(field, parent, shift)
    if not (0 < count <= parent):
        raise ValueError("invalid exact count")
    basis = normalized_lch_polynomials(field.name, parent)
    return tuple(
        tuple(
            polynomial_evaluate(
                field, basis[column], shift ^ row, counter
            )
            for column in range(count)
        )
        for row in range(count)
    )


def invert_matrix(
    field: BinaryField, matrix: Sequence[Sequence[int]],
    counter: Counter | None = None,
) -> tuple[tuple[int, ...], ...]:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    work = [list(row) + [int(row_index == column)
                         for column in range(size)]
            for row_index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        if counter is not None:
            counter.inversions += 1
        scale = field.inverse(work[column][column])
        for index in range(2 * size):
            if counter is not None:
                counter.loads += 1
                counter.stores += 1
                if scale not in (0, 1) and work[column][index] != 0:
                    counter.multiplications += 1
            work[column][index] = field.multiply(work[column][index], scale)
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            factor = work[row][column]
            for index in range(2 * size):
                source = work[column][index]
                if counter is not None:
                    counter.loads += 2
                    counter.stores += 1
                    if factor not in (0, 1) and source != 0:
                        counter.multiplications += 1
                    if work[row][index] != 0 and source != 0:
                        counter.xors += 1
                if source:
                    work[row][index] ^= (
                        source if factor == 1
                        else field.multiply(factor, source)
                    )
    return tuple(tuple(row[size:]) for row in work)


def apply_matrix(
    field: BinaryField, matrix: Sequence[Sequence[int]],
    values: Sequence[int], counter: Counter | None = None,
) -> list[int]:
    if any(len(row) != len(values) for row in matrix):
        raise ValueError("matrix/vector dimensions differ")
    output: list[int] = []
    for row in matrix:
        result = 0
        have_term = False
        for coefficient, value in zip(row, values):
            if coefficient == 0:
                continue
            if counter is not None:
                counter.loads += 2
                if coefficient != 1:
                    counter.multiplications += 1
                if have_term:
                    counter.xors += 1
            product = value if coefficient == 1 else field.multiply(coefficient, value)
            result ^= product
            have_term = True
        if counter is not None:
            counter.stores += 1
        output.append(result)
    return output


@functools.lru_cache(maxsize=None)
def fft_skew(field_name: str) -> tuple[int, ...]:
    field = field_named(field_name)
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


def padded_lch_inverse(
    field: BinaryField, values: Sequence[int], shift: int,
    counter: Counter | None = None,
) -> list[int]:
    length = len(values)
    validate_geometry(field, length, shift)
    if length == 1:
        return list(values)
    half = length // 2
    left = padded_lch_inverse(field, values[:half], shift, counter)
    right = padded_lch_inverse(field, values[half:], shift + half, counter)
    multiplier = fft_skew(field.name)[shift + half - 1]
    output = [0] * length
    for index, (u, v) in enumerate(zip(left, right)):
        if counter is not None:
            counter.loads += 2
            counter.stores += 2
        if multiplier == 0:
            x = u
            y = u ^ v
            if counter is not None:
                counter.xors += 1
        elif multiplier == 1:
            x = v
            y = u ^ v
            if counter is not None:
                counter.xors += 1
        else:
            y = u ^ v
            x = u ^ field.multiply(multiplier, y)
            if counter is not None:
                counter.xors += 2
                counter.multiplications += 1
        output[index] = x
        output[half + index] = y
    return output


@dataclasses.dataclass(frozen=True)
class ParentDensePlan:
    parent: int
    q: int
    shift: int
    columns: tuple[tuple[int, ...], ...]


def build_parent_dense_plan(
    field: BinaryField, q: int, shift: int,
    counter: Counter | None = None,
) -> ParentDensePlan:
    parent = ceil_power_of_two(q)
    full = lch_evaluation_matrix(field, parent, shift, parent, counter)
    inverse = invert_matrix(field, full, counter)
    columns = tuple(tuple(row[column] for column in range(q))
                    for row in inverse)
    return ParentDensePlan(parent, q, shift, columns)


@dataclasses.dataclass(frozen=True)
class ExactDensePlan:
    parent: int
    q: int
    shift: int
    inverse: tuple[tuple[int, ...], ...]


def build_exact_dense_plan(
    field: BinaryField, q: int, shift: int,
    counter: Counter | None = None,
) -> ExactDensePlan:
    parent = ceil_power_of_two(q)
    matrix = lch_evaluation_matrix(field, q, shift, parent, counter)
    return ExactDensePlan(parent, q, shift, invert_matrix(field, matrix, counter))


@dataclasses.dataclass(frozen=True)
class NewtonPlan:
    parent: int
    q: int
    shift: int
    points: tuple[int, ...]
    inverse_denominators: tuple[tuple[int, ...], ...]
    basis_polynomials: tuple[tuple[int, ...], ...]
    lch_leading_inverses: tuple[int, ...]


def build_newton_plan(
    field: BinaryField, q: int, shift: int,
    counter: Counter | None = None,
) -> NewtonPlan:
    parent = ceil_power_of_two(q)
    validate_geometry(field, parent, shift)
    points = tuple(shift ^ index for index in range(q))
    denominators: list[tuple[int, ...]] = [tuple()]
    for order in range(1, q):
        row: list[int] = []
        for index in range(order, q):
            difference = points[index] ^ points[index - order]
            if counter is not None:
                counter.inversions += 1
            row.append(field.inverse(difference))
        denominators.append(tuple(row))
    basis_polynomials: list[tuple[int, ...]] = []
    basis: Polynomial = [1]
    for index in range(q):
        basis_polynomials.append(tuple(basis))
        if index + 1 < q:
            basis = polynomial_multiply(
                field, basis, (points[index], 1), counter
            )
    lch_basis = normalized_lch_polynomials(field.name, parent)
    leading_inverses: list[int] = []
    for index in range(q):
        if counter is not None:
            counter.inversions += 1
        leading_inverses.append(field.inverse(lch_basis[index][index]))
    return NewtonPlan(
        parent, q, shift, points, tuple(denominators),
        tuple(basis_polynomials), tuple(leading_inverses)
    )


def monomial_to_lch(
    field: BinaryField, polynomial: Sequence[int], parent: int,
    count: int, leading_inverses: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    if len(leading_inverses) != count:
        raise ValueError("LCH leading-inverse count mismatch")
    residual = list(polynomial) + [0] * max(0, count - len(polynomial))
    residual = residual[:count]
    basis = normalized_lch_polynomials(field.name, parent)
    output = [0] * count
    for index in range(count - 1, -1, -1):
        inverse = leading_inverses[index]
        value = field.multiply(residual[index], inverse)
        if counter is not None and inverse != 1:
            counter.multiplications += 1
        output[index] = value
        for degree, coefficient in enumerate(basis[index]):
            if counter is not None:
                counter.loads += 2
                counter.stores += 1
                if value not in (0, 1) and coefficient != 0:
                    counter.multiplications += 1
                if residual[degree] != 0 and value != 0 and coefficient != 0:
                    counter.xors += 1
            if value and coefficient:
                residual[degree] ^= (
                    coefficient if value == 1
                    else field.multiply(value, coefficient)
                )
    if any(residual):
        raise AssertionError("LCH triangular conversion left a remainder")
    return output


def exact_newton_to_lch(
    field: BinaryField, plan: NewtonPlan, values: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    if len(values) != plan.q:
        raise ValueError("Newton value count mismatch")
    divided = list(values)
    for order in range(1, plan.q):
        inverses = plan.inverse_denominators[order]
        for index in range(plan.q - 1, order - 1, -1):
            difference = divided[index] ^ divided[index - 1]
            inverse = inverses[index - order]
            divided[index] = field.multiply(difference, inverse)
            if counter is not None:
                counter.loads += 3
                counter.stores += 1
                counter.xors += 1
                if inverse != 1:
                    counter.multiplications += 1

    polynomial: Polynomial = [0]
    for coefficient, basis in zip(divided, plan.basis_polynomials):
        polynomial_scale_add(field, polynomial, basis, coefficient, counter)
    return monomial_to_lch(
        field, polynomial, plan.parent, plan.q,
        plan.lch_leading_inverses, counter
    )


@dataclasses.dataclass(frozen=True)
class LagrangePlan:
    parent: int
    q: int
    shift: int
    points: tuple[int, ...]
    cardinal_polynomials: tuple[tuple[int, ...], ...]
    lch_leading_inverses: tuple[int, ...]


def build_lagrange_plan(
    field: BinaryField, q: int, shift: int,
    counter: Counter | None = None,
) -> LagrangePlan:
    parent = ceil_power_of_two(q)
    validate_geometry(field, parent, shift)
    points = tuple(shift ^ index for index in range(q))
    cardinals: list[tuple[int, ...]] = []
    for index, point in enumerate(points):
        polynomial: Polynomial = [1]
        denominator = 1
        for other_index, other in enumerate(points):
            if other_index == index:
                continue
            polynomial = polynomial_multiply(
                field, polynomial, (other, 1), counter
            )
            difference = point ^ other
            denominator = field.multiply(denominator, difference)
            if counter is not None:
                counter.multiplications += int(difference != 1)
                counter.loads += 2
                counter.stores += 1
        if counter is not None:
            counter.inversions += 1
        scale = field.inverse(denominator)
        cardinals.append(tuple(field.multiply(value, scale)
                               for value in polynomial))
        if counter is not None:
            counter.multiplications += sum(value != 0 and scale != 1
                                           for value in polynomial)
    lch_basis = normalized_lch_polynomials(field.name, parent)
    leading_inverses: list[int] = []
    for index in range(q):
        if counter is not None:
            counter.inversions += 1
        leading_inverses.append(field.inverse(lch_basis[index][index]))
    return LagrangePlan(
        parent, q, shift, points, tuple(cardinals), tuple(leading_inverses)
    )


def exact_generic_lagrange(
    field: BinaryField, plan: LagrangePlan, values: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    if len(values) != plan.q:
        raise ValueError("Lagrange value count mismatch")
    polynomial: Polynomial = [0] * plan.q
    for value, cardinal in zip(values, plan.cardinal_polynomials):
        polynomial_scale_add(field, polynomial, cardinal, value, counter)
    return monomial_to_lch(
        field, polynomial, plan.parent, plan.q,
        plan.lch_leading_inverses, counter
    )


def lch_coefficients_to_polynomial(
    field: BinaryField, coefficients: Sequence[int], parent: int,
) -> Polynomial:
    polynomial: Polynomial = [0]
    for coefficient, basis in zip(
        coefficients, normalized_lch_polynomials(field.name, parent)
    ):
        polynomial_scale_add(field, polynomial, basis, coefficient)
    return trim(polynomial)


def evaluate_lch(
    field: BinaryField, coefficients: Sequence[int], parent: int,
    points: Sequence[int],
) -> list[int]:
    polynomial = lch_coefficients_to_polynomial(field, coefficients, parent)
    return [polynomial_evaluate(field, polynomial, point) for point in points]


def stable_seed(field_name: str, q: int, shift: int) -> int:
    payload = f"c3:{field_name}:{q}:{shift}".encode("ascii")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def value_vectors(field: BinaryField, q: int, shift: int) -> Iterable[tuple[int, ...]]:
    if field.name == "gf4" and q <= 3:
        total = field.order ** q
        for encoded in range(total):
            value = encoded
            vector = []
            for _ in range(q):
                vector.append(value % field.order)
                value //= field.order
            yield tuple(vector)
        return
    yield (0,) * q
    yield (1,) * q
    for index in range(q):
        yield tuple(1 if column == index else 0 for column in range(q))
    yield tuple((column * 29 + q * 7 + shift) % field.order
                for column in range(q))
    random = __import__("random").Random(stable_seed(field.name, q, shift))
    for _ in range(6):
        yield tuple(random.randrange(field.order) for _ in range(q))


def correctness_job(specification: tuple[str, int, int]) -> dict[str, object]:
    field_name, q, shift = specification
    field = field_named(field_name)
    parent = ceil_power_of_two(q)
    validate_geometry(field, parent, shift)
    parent_plan = build_parent_dense_plan(field, q, shift)
    exact_plan = build_exact_dense_plan(field, q, shift)
    newton_plan = build_newton_plan(field, q, shift)
    lagrange_plan = build_lagrange_plan(field, q, shift)
    points = tuple(shift ^ index for index in range(q))
    parent_points = tuple(shift ^ index for index in range(parent))
    digest = hashlib.sha256()
    vectors = 0
    parent_exact_differences = 0
    exact_suffix_nonzero = 0
    comparisons = 0

    for values_tuple in value_vectors(field, q, shift):
        values = list(values_tuple)
        padded_input = values + [0] * (parent - q)
        parent_padded = padded_lch_inverse(field, padded_input, shift)
        parent_dense = apply_matrix(field, parent_plan.columns, values)
        if parent_padded != parent_dense:
            raise AssertionError("parent dense conversion differs from padded inverse")
        if evaluate_lch(field, parent_padded, parent, parent_points) != padded_input:
            raise AssertionError("parent inverse fails direct polynomial evaluation")

        exact_dense = apply_matrix(field, exact_plan.inverse, values)
        exact_newton = exact_newton_to_lch(field, newton_plan, values)
        exact_lagrange = exact_generic_lagrange(field, lagrange_plan, values)
        if not (exact_dense == exact_newton == exact_lagrange):
            raise AssertionError("exact inverse methods disagree")
        if evaluate_lch(field, exact_dense, parent, points) != values:
            raise AssertionError("exact coefficients fail direct evaluation")

        if q == parent:
            if parent_padded != exact_dense:
                raise AssertionError("full parent and exact inverse disagree")
        elif parent_padded[:q] != exact_dense or any(parent_padded[q:]):
            parent_exact_differences += 1
        if q < parent:
            suffix = evaluate_lch(
                field, exact_dense, parent, parent_points[q:]
            )
            exact_suffix_nonzero += int(any(suffix))

        for vector in (parent_padded, parent_dense, exact_dense,
                       exact_newton, exact_lagrange):
            for item in vector:
                digest.update(int(item).to_bytes(2, "little"))
        vectors += 1
        comparisons += 5 * q + 2 * parent

    return {
        "comparisons": comparisons,
        "digest": digest.hexdigest(),
        "exact_suffix_nonzero_vectors": exact_suffix_nonzero,
        "field": field_name,
        "parent": parent,
        "parent_exact_difference_vectors": parent_exact_differences,
        "q": q,
        "shift": shift,
        "vectors": vectors,
    }


def shifts_for(field: BinaryField, parent: int) -> tuple[int, ...]:
    candidates = (
        0,
        (field.order // 2) & ~(parent - 1),
        field.order - parent,
    )
    return tuple(sorted(set(value for value in candidates
                            if 0 <= value <= field.order - parent)))


def correctness_specs() -> list[tuple[str, int, int]]:
    selections = {
        "gf4": tuple(range(1, 17)),
        "gf8": (1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31),
        "gf16": (3, 5, 9, 17),
    }
    specs: list[tuple[str, int, int]] = []
    for field_name, counts in selections.items():
        field = field_named(field_name)
        for q in counts:
            parent = ceil_power_of_two(q)
            for shift in shifts_for(field, parent):
                specs.append((field_name, q, shift))
    return specs


def source_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def operation_row(field_name: str, q: int, shift: int) -> dict[str, object]:
    field = field_named(field_name)
    parent = ceil_power_of_two(q)
    rows: list[dict[str, object]] = []

    padded_counter = Counter()
    padded_lch_inverse(field, [0] * parent, shift, padded_counter)
    rows.append({
        "method": METHODS[0], "setup": dataclasses.asdict(Counter()),
        "execution": dataclasses.asdict(padded_counter),
        "table_elements": len(fft_skew(field_name)),
        "output_coefficients": parent,
    })

    if parent <= 32:
        parent_setup = Counter()
        parent_plan = build_parent_dense_plan(field, q, shift, parent_setup)
        parent_execution = Counter()
        apply_matrix(field, parent_plan.columns, [0] * q, parent_execution)
        exact_setup = Counter()
        exact_plan = build_exact_dense_plan(field, q, shift, exact_setup)
        exact_execution = Counter()
        apply_matrix(field, exact_plan.inverse, [0] * q, exact_execution)
    else:
        parent_matrix = parent * parent * (parent + 1) // 2
        parent_setup = Counter(
            inversions=parent,
            multiplications=parent ** 3 + parent_matrix,
            xors=parent ** 3 + parent_matrix,
            loads=2 * parent ** 3 + parent_matrix,
            stores=parent ** 3 + parent_matrix,
        )
        parent_execution = Counter(
            multiplications=parent * q, xors=parent * max(0, q - 1),
            loads=2 * parent * q, stores=parent,
        )
        exact_matrix = q * q * (q + 1) // 2
        exact_setup = Counter(
            inversions=q,
            multiplications=q ** 3 + exact_matrix,
            xors=q ** 3 + exact_matrix,
            loads=2 * q ** 3 + exact_matrix,
            stores=q ** 3 + exact_matrix,
        )
        exact_execution = Counter(
            multiplications=q * q, xors=q * max(0, q - 1),
            loads=2 * q * q, stores=q,
        )
    rows.append({
        "method": METHODS[1], "setup": dataclasses.asdict(parent_setup),
        "execution": dataclasses.asdict(parent_execution),
        "table_elements": parent * q, "output_coefficients": parent,
    })
    rows.append({
        "method": METHODS[2], "setup": dataclasses.asdict(exact_setup),
        "execution": dataclasses.asdict(exact_execution),
        "table_elements": q * q, "output_coefficients": q,
    })

    if parent <= 32:
        probe = [((index + 1) * 29) % field.order or 1
                 for index in range(q)]
        newton_setup = Counter()
        newton_plan = build_newton_plan(field, q, shift, newton_setup)
        newton_execution = Counter()
        exact_newton_to_lch(field, newton_plan, probe, newton_execution)
    else:
        divided = q * max(0, q - 1) // 2
        # Setup precomputes every successive Newton-basis polynomial and the
        # triangular LCH leading-coefficient inverses.  Execution only applies
        # those immutable tables.  This remains a dense scalar conversion, not
        # Coxon's fast recursive algorithm.
        basis_work = q * max(0, q - 1)
        triangular = q * (q + 1) // 2
        newton_setup = Counter(
            xors=basis_work, multiplications=basis_work,
            inversions=divided + q, loads=2 * basis_work,
            stores=basis_work,
        )
        newton_execution = Counter(
            xors=divided + 2 * triangular,
            multiplications=divided + 2 * triangular + q,
            loads=3 * divided + 4 * triangular,
            stores=divided + 2 * triangular,
        )
    rows.append({
        "method": METHODS[3], "setup": dataclasses.asdict(newton_setup),
        "execution": dataclasses.asdict(newton_execution),
        "table_elements": (
            q * max(0, q - 1) // 2 + q * (q + 1) // 2 + q
        ),
        "output_coefficients": q,
    })

    if parent <= 32:
        lagrange_setup = Counter()
        lagrange_plan = build_lagrange_plan(field, q, shift, lagrange_setup)
        lagrange_execution = Counter()
        exact_generic_lagrange(field, lagrange_plan, probe, lagrange_execution)
    else:
        cardinal_work = q * q * max(0, q - 1)
        dense = q * q
        triangular = q * (q + 1) // 2
        lagrange_setup = Counter(
            xors=cardinal_work, multiplications=cardinal_work,
            inversions=2 * q, loads=2 * cardinal_work,
            stores=cardinal_work,
        )
        lagrange_execution = Counter(
            xors=dense + triangular,
            multiplications=dense + triangular,
            loads=2 * (dense + triangular), stores=dense + triangular,
        )
    rows.append({
        "method": METHODS[4], "setup": dataclasses.asdict(lagrange_setup),
        "execution": dataclasses.asdict(lagrange_execution),
        "table_elements": q * q + q, "output_coefficients": q,
    })
    return {"field": field_name, "parent": parent, "q": q,
            "shift": shift, "methods": rows}


def cost(counter: dict[str, int]) -> int:
    return (
        counter["multiplications"] * 3 + counter["xors"] +
        counter["loads"] + counter["stores"] * 2 +
        counter["inversions"] * 12
    )


def modeled_cells(operation_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    cells: list[dict[str, object]] = []
    for operation in operation_rows:
        field_name = str(operation["field"])
        bytes_per_symbol = 1 if field_name == "gf8" else 2
        by_method = {str(row["method"]): row for row in operation["methods"]}
        for shard_bytes in (64, 1024, 65536):
            symbols = shard_bytes // bytes_per_symbol
            for batch in (1, 8, 64):
                for reuse in (1, 8, 64, 1024):
                    baseline = by_method[METHODS[0]]
                    baseline_score = (
                        cost(baseline["execution"]) * symbols * batch +
                        cost(baseline["setup"]) / reuse
                    )
                    for method in METHODS:
                        row = by_method[method]
                        score = (
                            cost(row["execution"]) * symbols * batch +
                            cost(row["setup"]) / reuse
                        )
                        cells.append({
                            "batch": batch,
                            "field": field_name,
                            "gain_over_padded": round(
                                baseline_score / score, 9
                            ) if score else 1.0,
                            "method": method,
                            "parent": operation["parent"],
                            "q": operation["q"],
                            "reuse": reuse,
                            "score": round(score, 3),
                            "shard_bytes": shard_bytes,
                        })
    return cells


def deterministic_analysis(workers: int) -> dict[str, object]:
    specs = correctness_specs()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(correctness_job, specs))
    results.sort(key=lambda row: (str(row["field"]), int(row["q"]),
                                  int(row["shift"])))

    operation_rows: list[dict[str, object]] = []
    for field_name, counts in (
        ("gf8", (3, 5, 7, 9, 15, 17, 31, 33, 63, 65, 127, 129, 255)),
        ("gf16", (3, 5, 9, 17, 33, 65, 129, 257)),
    ):
        field = field_named(field_name)
        for q in counts:
            parent = ceil_power_of_two(q)
            operation_rows.append(operation_row(field_name, q,
                                                  shifts_for(field, parent)[0]))
    cells = modeled_cells(operation_rows)
    compatible_candidates = [cell for cell in cells
                             if cell["method"] == METHODS[1]]
    exact_candidates = [cell for cell in cells
                        if cell["method"] in METHODS[2:]]
    best_compatible = max(compatible_candidates,
                          key=lambda cell: float(cell["gain_over_padded"]))
    best_exact = max(exact_candidates,
                     key=lambda cell: float(cell["gain_over_padded"]))

    summary = {
        "comparisons": sum(int(row["comparisons"]) for row in results),
        "correctness_jobs": len(results),
        "exact_suffix_nonzero_vectors": sum(
            int(row["exact_suffix_nonzero_vectors"]) for row in results
        ),
        "modeled_cells": len(cells),
        "parent_exact_difference_vectors": sum(
            int(row["parent_exact_difference_vectors"]) for row in results
        ),
        "vectors": sum(int(row["vectors"]) for row in results),
    }
    return {
        "best_modeled_compatible": best_compatible,
        "best_modeled_exact": best_exact,
        "correctness": results,
        "disposition": {
            "exact_profile_methods": "inconclusive-new-profile",
            "parent_dense_lagrange_to_lch": "inconclusive-tiny-q-only",
            "production_promotion": "none",
            "true_fast_coxon_or_epsilon_inverse": "not-implemented",
        },
        "method_class": {
            METHODS[0]: "wire-compatible-parent",
            METHODS[1]: "wire-compatible-parent",
            METHODS[2]: "new-exact-profile",
            METHODS[3]: "new-exact-profile",
            METHODS[4]: "new-exact-profile",
        },
        "modeled_cells": cells,
        "operation_rows": operation_rows,
        "schema": SCHEMA,
        "source_sha256": source_sha256(),
        "summary": summary,
    }


def canonical_write(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, sort_keys=True, indent=2, separators=(",", ": ")) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def run_analysis(arguments: argparse.Namespace) -> int:
    if arguments.workers <= 0:
        raise ValueError("workers must be positive")
    value = deterministic_analysis(arguments.workers)
    canonical_write(Path(arguments.output), value)
    print(json.dumps(value["summary"], sort_keys=True))
    return 0


def verify_result(arguments: argparse.Namespace) -> int:
    value = json.loads(Path(arguments.result).read_text(encoding="utf-8"))
    if value.get("schema") != SCHEMA:
        raise ValueError("result schema mismatch")
    if value.get("source_sha256") != source_sha256():
        raise ValueError("result is stale for this source")
    if value.get("disposition", {}).get("production_promotion") != "none":
        raise ValueError("unexpected production promotion")
    if not value.get("correctness") or not value.get("modeled_cells"):
        raise ValueError("result is incomplete")
    if any(row.get("vectors", 0) <= 0 for row in value["correctness"]):
        raise ValueError("empty correctness job")
    expected = deterministic_analysis(1)
    if value != expected:
        raise ValueError("result content differs from deterministic analysis")
    print("leopard2 C3 deterministic result verified")
    return 0


def median_mad(values: Sequence[int]) -> tuple[int, int]:
    median = int(statistics.median(values))
    mad = int(statistics.median(abs(value - median) for value in values))
    return median, mad


def time_call(function, iterations: int) -> int:
    start = time.perf_counter_ns()
    for _ in range(iterations):
        function()
    return time.perf_counter_ns() - start


def benchmark_case(field: BinaryField, q: int, samples: int,
                   iterations: int) -> dict[str, object]:
    parent = ceil_power_of_two(q)
    shift = 0
    random_source = random.Random(stable_seed(field.name, q, shift))
    values = [random_source.randrange(field.order) for _ in range(q)]
    padded_values = values + [0] * (parent - q)

    builders = {
        METHODS[1]: lambda: build_parent_dense_plan(field, q, shift),
        METHODS[2]: lambda: build_exact_dense_plan(field, q, shift),
        METHODS[3]: lambda: build_newton_plan(field, q, shift),
        METHODS[4]: lambda: build_lagrange_plan(field, q, shift),
    }
    setup: dict[str, dict[str, int]] = {}
    plans: dict[str, object] = {}
    for method, builder in builders.items():
        times = []
        for _ in range(samples):
            start = time.perf_counter_ns()
            plan = builder()
            times.append(time.perf_counter_ns() - start)
        plans[method] = plan
        median, mad = median_mad(times)
        setup[method] = {"mad_ns": mad, "median_ns": median}
    setup[METHODS[0]] = {"mad_ns": 0, "median_ns": 0}

    functions = {
        METHODS[0]: lambda: padded_lch_inverse(field, padded_values, shift),
        METHODS[1]: lambda: apply_matrix(
            field, plans[METHODS[1]].columns, values
        ),
        METHODS[2]: lambda: apply_matrix(
            field, plans[METHODS[2]].inverse, values
        ),
        METHODS[3]: lambda: exact_newton_to_lch(
            field, plans[METHODS[3]], values
        ),
        METHODS[4]: lambda: exact_generic_lagrange(
            field, plans[METHODS[4]], values
        ),
    }
    execution: dict[str, dict[str, int]] = {}
    for method, function in functions.items():
        function()
        times = [time_call(function, iterations) // iterations
                 for _ in range(samples)]
        median, mad = median_mad(times)
        execution[method] = {"mad_ns": mad, "median_ns": median}
    return {
        "execution": execution,
        "field": field.name,
        "iterations_per_sample": iterations,
        "parent": parent,
        "q": q,
        "samples": samples,
        "setup": setup,
        "shift": shift,
    }


def affinity() -> list[int]:
    if hasattr(os, "sched_getaffinity"):
        return sorted(os.sched_getaffinity(0))
    return []


def run_benchmark(arguments: argparse.Namespace) -> int:
    if arguments.samples < 3 or arguments.iterations <= 0:
        raise ValueError("benchmark requires samples>=3 and iterations>0")
    before = affinity()
    if len(before) > 1:
        raise ValueError("benchmark must be pinned to one allowed CPU")
    field = field_named(arguments.field)
    cases = [benchmark_case(field, q, arguments.samples, arguments.iterations)
             for q in arguments.counts]
    after = affinity()
    if after != before:
        raise RuntimeError("process affinity changed during benchmark")
    value = {
        "affinity": before,
        "cases": cases,
        "machine": {
            "architecture": platform.machine(),
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python": platform.python_version(),
        },
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "schema": TIMING_SCHEMA,
        "source_sha256": source_sha256(),
    }
    canonical_write(Path(arguments.output), value)
    print(f"wrote {arguments.output}")
    return 0


def self_test() -> int:
    field = field_named("gf4")
    if field.multiply(7, 11) != field.multiply(11, 7):
        raise AssertionError("GF4 multiplication is not commutative")
    for q in range(1, 8):
        parent = ceil_power_of_two(q)
        values = [(index * 3 + q) % field.order for index in range(q)]
        padded = padded_lch_inverse(field, values + [0] * (parent - q), 0)
        parent_plan = build_parent_dense_plan(field, q, 0)
        exact_plan = build_exact_dense_plan(field, q, 0)
        newton = build_newton_plan(field, q, 0)
        lagrange = build_lagrange_plan(field, q, 0)
        if padded != apply_matrix(field, parent_plan.columns, values):
            raise AssertionError("parent dense self-test differs")
        exact = apply_matrix(field, exact_plan.inverse, values)
        if exact != exact_newton_to_lch(field, newton, values):
            raise AssertionError("Newton self-test differs")
        if exact != exact_generic_lagrange(field, lagrange, values):
            raise AssertionError("Lagrange self-test differs")
        if q == parent and exact != padded:
            raise AssertionError("full parent self-test differs")
    print("leopard2 C3 self-test passed")
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)
    subparsers.add_parser("self-test")
    analyze = subparsers.add_parser("analyze")
    analyze.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 128))
    analyze.add_argument("--output", required=True)
    verify = subparsers.add_parser("verify")
    verify.add_argument("result")
    benchmark = subparsers.add_parser("benchmark")
    benchmark.add_argument("--field", choices=("gf8", "gf16"), default="gf8")
    benchmark.add_argument("--counts", type=int, nargs="+", default=(3, 5, 9, 17))
    benchmark.add_argument("--samples", type=int, default=9)
    benchmark.add_argument("--iterations", type=int, default=32)
    benchmark.add_argument("--output", required=True)
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        if arguments.command == "self-test":
            return self_test()
        if arguments.command == "analyze":
            return run_analysis(arguments)
        if arguments.command == "verify":
            return verify_result(arguments)
        if arguments.command == "benchmark":
            return run_benchmark(arguments)
        raise AssertionError("unreachable command")
    except (AssertionError, OSError, ValueError, ZeroDivisionError) as error:
        print(f"C3 error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
