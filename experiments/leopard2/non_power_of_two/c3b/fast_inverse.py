#!/usr/bin/env python3
"""C3b executable Coxon and Tang--Han inverse-transform research.

This module is deliberately isolated from the production library.  It contains
literal scalar implementations of these published algorithms:

* Coxon Algorithm 1 (normalised Newton to LCH), using a reduction tree;
* Coxon Algorithm 3 (truncated Lagrange to LCH), including its mixed recursive
  states; and
* Tang--Han Appendix-A Algorithm 8 (epsilon-point inverse LCH transform).

All three solve exact-prefix interpolation.  For non-power-of-two ``count``
this is a different mathematical object from Leopard's wire-compatible padded
parent inverse, which fixes the missing *evaluations* to zero.  The distinction
is represented in every result row and is tested rather than inferred.

The field implementation uses carryless polynomial multiplication behind an
explicit coordinate-basis map.  It does not call Leopard or its transform code.
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
from typing import Callable, Iterable, Sequence


SCHEMA = "leopard2-c3b-fast-inverse/v1"
TIMING_SCHEMA = "leopard2-c3b-fast-inverse-timing/v1"
METHOD_PADDED = "parent_padded_lch"
METHOD_DENSE = "exact_dense_inverse"
METHOD_GENERIC = "exact_generic_interpolation"
METHOD_COXON_L = "exact_coxon_algorithm3_lagrange_to_lch"
METHOD_COXON_N = "exact_coxon_algorithm1_newton_to_lch"
METHOD_TANG = "exact_tang_han_algorithm8_epsilon_ifft"
METHODS = (
    METHOD_PADDED,
    METHOD_DENSE,
    METHOD_GENERIC,
    METHOD_COXON_L,
    METHOD_COXON_N,
    METHOD_TANG,
)


@dataclasses.dataclass
class Counter:
    xors: int = 0
    multiplications: int = 0
    inversions: int = 0
    loads: int = 0
    stores: int = 0
    recursive_calls: int = 0
    full_transform_symbols: int = 0

    def add(self, other: "Counter") -> None:
        for field in dataclasses.fields(self):
            setattr(self, field.name,
                    getattr(self, field.name) + getattr(other, field.name))


class BinaryField:
    """GF(2^m) in the public coordinate representation under test."""

    def __init__(
        self,
        name: str,
        bits: int,
        polynomial: int,
        coordinate_basis: Sequence[int],
    ) -> None:
        if not (1 <= bits <= 16) or len(coordinate_basis) != bits:
            raise ValueError("invalid binary field")
        if polynomial & (1 << bits) == 0:
            raise ValueError("field polynomial has the wrong degree")
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

        self.mul_rows: tuple[bytes, ...] | None = None
        if self.order <= 256:
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
        if not (0 <= left < self.order and 0 <= right < self.order):
            raise ValueError("field element outside representation")
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

    def divide(self, numerator: int, denominator: int) -> int:
        return self.multiply(numerator, self.inverse(denominator))

    def frobenius(self, value: int, dimension: int) -> int:
        for _ in range(dimension):
            value = self.multiply(value, value)
        return value


@functools.lru_cache(maxsize=None)
def field_named(name: str) -> BinaryField:
    if name == "gf4":
        # Test-only field.  Its polynomial-coordinate basis intentionally is
        # not assumed to be Cantor; the Coxon tree builder handles it generally.
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
        raise ValueError("count must be positive")
    return 1 << (value - 1).bit_length()


Polynomial = list[int]


def trim(polynomial: Polynomial) -> Polynomial:
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def polynomial_multiply(
    field: BinaryField,
    left: Sequence[int],
    right: Sequence[int],
    counter: Counter | None = None,
) -> Polynomial:
    output = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if counter is not None:
                counter.loads += 2
                counter.stores += 1
                counter.multiplications += int(a != 0 and b not in (0, 1))
                counter.xors += int(output[i + j] != 0 and a != 0 and b != 0)
            if a and b:
                output[i + j] ^= a if b == 1 else field.multiply(a, b)
    return trim(output)


def polynomial_evaluate(
    field: BinaryField,
    polynomial: Sequence[int],
    point: int,
    counter: Counter | None = None,
) -> int:
    value = 0
    for coefficient in reversed(polynomial):
        if counter is not None:
            counter.loads += 1
            counter.stores += 1
            counter.multiplications += int(point not in (0, 1) and value != 0)
            counter.xors += int(value != 0 and coefficient != 0)
        value = (value if point == 1 else field.multiply(value, point)) ^ coefficient
    return value


@functools.lru_cache(maxsize=None)
def subspace_polynomial(field_name: str, dimension: int) -> tuple[int, ...]:
    field = field_named(field_name)
    polynomial: Polynomial = [1]
    for point in range(1 << dimension):
        polynomial = polynomial_multiply(field, polynomial, (point, 1))
    return tuple(polynomial)


@functools.lru_cache(maxsize=None)
def normalized_lch_polynomials(
    field_name: str, parent: int
) -> tuple[tuple[int, ...], ...]:
    field = field_named(field_name)
    if parent <= 0 or parent & (parent - 1) or parent > field.order:
        raise ValueError("invalid LCH parent")
    dimension = parent.bit_length() - 1
    subspaces = tuple(subspace_polynomial(field_name, bit)
                      for bit in range(dimension))
    output: list[tuple[int, ...]] = []
    for index in range(parent):
        polynomial: Polynomial = [1]
        normalizer = 1
        for bit, subspace in enumerate(subspaces):
            if index & (1 << bit):
                polynomial = polynomial_multiply(field, polynomial, subspace)
                normalizer = field.multiply(
                    normalizer,
                    polynomial_evaluate(field, subspace, 1 << bit),
                )
        scale = field.inverse(normalizer)
        output.append(tuple(field.multiply(value, scale)
                            for value in polynomial))
    return tuple(output)


@functools.lru_cache(maxsize=None)
def _subspace_multiplier_value(
    field_name: str, dimension: int, shift: int
) -> int:
    field = field_named(field_name)
    if dimension < 0 or dimension >= field.bits:
        raise ValueError("invalid subspace dimension")
    polynomial = subspace_polynomial(field_name, dimension)
    numerator = polynomial_evaluate(field, polynomial, shift)
    denominator = polynomial_evaluate(
        field, polynomial, 1 << dimension
    )
    return field.divide(numerator, denominator)


def subspace_multiplier(
    field: BinaryField, dimension: int, shift: int,
    counter: Counter | None = None,
) -> int:
    """Return precomputed s_dimension(shift) / s_dimension(v_dimension).

    The cached value models Leopard's immutable skew/normalizer table.  When a
    setup counter is supplied, account for constructing that constant once;
    transform execution itself performs only the fixed multiplication.
    """
    value = _subspace_multiplier_value(field.name, dimension, shift)
    if counter is not None:
        polynomial = subspace_polynomial(field.name, dimension)
        polynomial_evaluate(field, polynomial, shift, counter)
        polynomial_evaluate(field, polynomial, 1 << dimension, counter)
        counter.inversions += 1
        counter.multiplications += int(value != 0)
    return value


def _fixed_product(
    field: BinaryField, coefficient: int, value: int,
    counter: Counter | None,
) -> int:
    if coefficient == 0 or value == 0:
        return 0
    if coefficient == 1:
        return value
    if counter is not None:
        counter.multiplications += 1
    return field.multiply(coefficient, value)


def lch_forward(
    field: BinaryField,
    coefficients: Sequence[int],
    shift: int,
    counter: Counter | None = None,
) -> list[int]:
    length = len(coefficients)
    if length <= 0 or length & (length - 1) or length > field.order:
        raise ValueError("forward length is not a supported power of two")
    if counter is not None:
        counter.recursive_calls += 1
        counter.full_transform_symbols += length
    if length == 1:
        return list(coefficients)
    half = length // 2
    dimension = half.bit_length() - 1
    multiplier = subspace_multiplier(field, dimension, shift)
    left_coefficients = [0] * half
    right_coefficients = [0] * half
    for index in range(half):
        low = coefficients[index]
        high = coefficients[half + index]
        first = low ^ _fixed_product(field, multiplier, high, counter)
        second = first ^ high
        if counter is not None:
            counter.loads += 2
            counter.stores += 2
            counter.xors += 2
        left_coefficients[index] = first
        right_coefficients[index] = second
    return (
        lch_forward(field, left_coefficients, shift, counter) +
        lch_forward(field, right_coefficients, shift ^ half, counter)
    )


def lch_inverse(
    field: BinaryField,
    values: Sequence[int],
    shift: int,
    counter: Counter | None = None,
) -> list[int]:
    length = len(values)
    if length <= 0 or length & (length - 1) or length > field.order:
        raise ValueError("inverse length is not a supported power of two")
    if counter is not None:
        counter.recursive_calls += 1
        counter.full_transform_symbols += length
    if length == 1:
        return list(values)
    half = length // 2
    dimension = half.bit_length() - 1
    low = lch_inverse(field, values[:half], shift, counter)
    high = lch_inverse(field, values[half:], shift ^ half, counter)
    multiplier = subspace_multiplier(field, dimension, shift)
    output = [0] * length
    for index, (first, second) in enumerate(zip(low, high)):
        upper = first ^ second
        lower = first ^ _fixed_product(field, multiplier, upper, counter)
        if counter is not None:
            counter.loads += 2
            counter.stores += 2
            counter.xors += 2
        output[index] = lower
        output[half + index] = upper
    return output


def lch_to_monomial(
    field: BinaryField, coefficients: Sequence[int], parent: int
) -> Polynomial:
    polynomial: Polynomial = [0] * parent
    for value, basis in zip(
        coefficients, normalized_lch_polynomials(field.name, parent)
    ):
        if value == 0:
            continue
        for degree, coefficient in enumerate(basis):
            if coefficient:
                polynomial[degree] ^= (
                    coefficient if value == 1
                    else field.multiply(value, coefficient)
                )
    return trim(polynomial)


def evaluate_lch(
    field: BinaryField,
    coefficients: Sequence[int],
    parent: int,
    points: Sequence[int],
) -> list[int]:
    polynomial = lch_to_monomial(field, coefficients, parent)
    return [polynomial_evaluate(field, polynomial, point) for point in points]


def evaluation_matrix(
    field: BinaryField,
    count: int,
    shift: int,
    parent: int,
    counter: Counter | None = None,
) -> tuple[tuple[int, ...], ...]:
    basis = normalized_lch_polynomials(field.name, parent)
    return tuple(
        tuple(polynomial_evaluate(field, basis[column], shift ^ row, counter)
              for column in range(count))
        for row in range(count)
    )


def invert_matrix(
    field: BinaryField,
    matrix: Sequence[Sequence[int]],
    counter: Counter | None = None,
) -> tuple[tuple[int, ...], ...]:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix is not square")
    work = [list(row) + [int(row_index == column)
                         for column in range(size)]
            for row_index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            raise ValueError("singular matrix")
        work[column], work[pivot] = work[pivot], work[column]
        scale = field.inverse(work[column][column])
        if counter is not None:
            counter.inversions += 1
        for index in range(2 * size):
            work[column][index] = _fixed_product(
                field, scale, work[column][index], counter
            )
            if counter is not None:
                counter.loads += 1
                counter.stores += 1
        for row in range(size):
            factor = work[row][column]
            if row == column or factor == 0:
                continue
            for index in range(2 * size):
                product = _fixed_product(
                    field, factor, work[column][index], counter
                )
                work[row][index] ^= product
                if counter is not None:
                    counter.loads += 2
                    counter.stores += 1
                    counter.xors += 1
    return tuple(tuple(row[size:]) for row in work)


def apply_matrix(
    field: BinaryField,
    matrix: Sequence[Sequence[int]],
    values: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    if any(len(row) != len(values) for row in matrix):
        raise ValueError("matrix/vector dimensions differ")
    output: list[int] = []
    for row in matrix:
        accumulator = 0
        first = True
        for coefficient, value in zip(row, values):
            if coefficient == 0:
                continue
            product = _fixed_product(field, coefficient, value, counter)
            if not first and counter is not None:
                counter.xors += 1
            accumulator ^= product
            first = False
            if counter is not None:
                counter.loads += 2
        if counter is not None:
            counter.stores += 1
        output.append(accumulator)
    return output


@dataclasses.dataclass(frozen=True)
class DensePlan:
    parent: int
    count: int
    shift: int
    inverse: tuple[tuple[int, ...], ...]


def build_dense_plan(
    field: BinaryField, count: int, shift: int,
    counter: Counter | None = None,
) -> DensePlan:
    parent = ceil_power_of_two(count)
    if parent > field.order:
        raise ValueError("parent exceeds field")
    matrix = evaluation_matrix(field, count, shift, parent, counter)
    return DensePlan(parent, count, shift, invert_matrix(field, matrix, counter))


def dense_inverse(
    field: BinaryField, plan: DensePlan, values: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    return apply_matrix(field, plan.inverse, values, counter) + [0] * (
        plan.parent - plan.count
    )


def monomial_to_lch(
    field: BinaryField,
    polynomial: Sequence[int],
    parent: int,
    count: int,
    leading_inverses: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    residual = list(polynomial) + [0] * max(0, count - len(polynomial))
    residual = residual[:count]
    basis = normalized_lch_polynomials(field.name, parent)
    output = [0] * count
    for index in range(count - 1, -1, -1):
        value = _fixed_product(
            field, leading_inverses[index], residual[index], counter
        )
        output[index] = value
        for degree, coefficient in enumerate(basis[index]):
            product = _fixed_product(field, value, coefficient, counter)
            residual[degree] ^= product
            if counter is not None:
                counter.loads += 2
                counter.stores += 1
                counter.xors += 1
    if any(residual):
        raise AssertionError("triangular LCH conversion left a remainder")
    return output + [0] * (parent - count)


@dataclasses.dataclass(frozen=True)
class GenericPlan:
    parent: int
    count: int
    shift: int
    cardinals: tuple[tuple[int, ...], ...]
    leading_inverses: tuple[int, ...]


def build_generic_plan(
    field: BinaryField, count: int, shift: int,
    counter: Counter | None = None,
) -> GenericPlan:
    parent = ceil_power_of_two(count)
    points = tuple(shift ^ index for index in range(count))
    cardinals: list[tuple[int, ...]] = []
    for index, point in enumerate(points):
        polynomial: Polynomial = [1]
        denominator = 1
        for other_index, other in enumerate(points):
            if index == other_index:
                continue
            polynomial = polynomial_multiply(
                field, polynomial, (other, 1), counter
            )
            denominator = field.multiply(denominator, point ^ other)
            if counter is not None:
                counter.multiplications += 1
        scale = field.inverse(denominator)
        if counter is not None:
            counter.inversions += 1
        cardinals.append(tuple(_fixed_product(field, scale, value, counter)
                               for value in polynomial))
    basis = normalized_lch_polynomials(field.name, parent)
    leading = []
    for index in range(count):
        leading.append(field.inverse(basis[index][index]))
        if counter is not None:
            counter.inversions += 1
    return GenericPlan(parent, count, shift, tuple(cardinals), tuple(leading))


def generic_inverse(
    field: BinaryField, plan: GenericPlan, values: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    polynomial = [0] * plan.count
    for value, cardinal in zip(values, plan.cardinals):
        for degree, coefficient in enumerate(cardinal):
            polynomial[degree] ^= _fixed_product(
                field, value, coefficient, counter
            )
            if counter is not None:
                counter.loads += 2
                counter.stores += 1
                counter.xors += 1
    return monomial_to_lch(
        field, polynomial, plan.parent, plan.count,
        plan.leading_inverses, counter
    )


@dataclasses.dataclass(frozen=True)
class ReductionNode:
    beta: tuple[int, ...]
    d: int
    alpha: "ReductionNode | None"
    delta: "ReductionNode | None"
    sigma_phi_alpha: tuple[tuple[int, ...], ...]

    @property
    def dimension(self) -> int:
        return len(self.beta)

    @property
    def size(self) -> int:
        return 1 << self.dimension

    @property
    def leaf(self) -> bool:
        return self.dimension == 1


def _valid_reduction(
    field: BinaryField, beta: Sequence[int], d: int
) -> bool:
    beta0_inverse = field.inverse(beta[0])
    for value in beta[:d]:
        ratio = field.multiply(value, beta0_inverse)
        if field.frobenius(ratio, d) != ratio:
            return False
    return True


def _phi_for(
    field: BinaryField, node: ReductionNode, shift: int,
    counter: Counter | None = None,
) -> tuple[int, ...]:
    if node.leaf:
        if counter is not None:
            counter.inversions += 1
            counter.multiplications += int(shift != 0)
        return (field.divide(shift, node.beta[0]),)
    assert node.alpha is not None and node.delta is not None
    beta0_inverse = field.inverse(node.beta[0])
    ratio = field.multiply(shift, beta0_inverse)
    delta_shift = field.frobenius(ratio, node.d) ^ ratio
    if counter is not None:
        counter.inversions += 1
        counter.multiplications += (1 + node.d) * int(shift != 0)
        counter.xors += 1
    return (
        _phi_for(field, node.alpha, shift, counter) +
        _phi_for(field, node.delta, delta_shift, counter)
    )


def build_reduction_tree(
    field: BinaryField,
    beta: Sequence[int],
    counter: Counter | None = None,
) -> ReductionNode:
    beta_tuple = tuple(beta)
    if not beta_tuple:
        raise ValueError("empty Coxon basis")
    if len(beta_tuple) == 1:
        return ReductionNode(beta_tuple, 0, None, None, tuple())

    candidates = [value for value in range(1, len(beta_tuple))
                  if value & (value - 1) == 0]
    valid = [value for value in candidates
             if _valid_reduction(field, beta_tuple, value)]
    if not valid:
        raise ValueError("no Coxon reduction exists")
    d = max(valid)
    beta0_inverse = field.inverse(beta_tuple[0])
    if counter is not None:
        counter.inversions += 1
    delta_beta = []
    for value in beta_tuple[d:]:
        ratio = field.multiply(value, beta0_inverse)
        delta_beta.append(field.frobenius(ratio, d) ^ ratio)
        if counter is not None:
            counter.multiplications += 1 + d
            counter.xors += 1
    alpha = build_reduction_tree(field, beta_tuple[:d], counter)
    delta = build_reduction_tree(field, delta_beta, counter)

    # Coxon sigma_(v,i) is the prefix XOR of gamma = beta[d:].  Algorithm
    # 3 only consumes phi values on the alpha leaves.  These values are plan
    # constants and are therefore excluded from byte-heavy execution.
    sigma_rows: list[tuple[int, ...]] = []
    sigma = 0
    for value in beta_tuple[d:]:
        sigma ^= value
        provisional = ReductionNode(
            beta_tuple, d, alpha, delta, tuple()
        )
        sigma_rows.append(_phi_for(field, provisional, sigma, counter)[:d])
        if counter is not None:
            counter.xors += 1
    return ReductionNode(beta_tuple, d, alpha, delta, tuple(sigma_rows))


def trailing_ones(value: int) -> int:
    result = 0
    while value & 1:
        result += 1
        value >>= 1
    return result


def _xor_vector(left: Sequence[int], right: Sequence[int],
                counter: Counter | None) -> tuple[int, ...]:
    if len(left) != len(right):
        raise ValueError("Coxon phi dimensions differ")
    if counter is not None:
        counter.loads += 2 * len(left)
        counter.stores += len(left)
        counter.xors += len(left)
    return tuple(a ^ b for a, b in zip(left, right))


def _leaf_multiply_add(
    field: BinaryField, values: list[int], destination: int,
    source: int, coefficient: int, counter: Counter | None,
) -> None:
    values[destination] ^= _fixed_product(
        field, coefficient, values[source], counter
    )
    if counter is not None:
        counter.loads += 2
        counter.stores += 1
        counter.xors += 1


def coxon_l2x(
    field: BinaryField,
    node: ReductionNode,
    phi: Sequence[int],
    c: int,
    ell: int,
    b: int,
    values: list[int],
    indices: Sequence[int],
    counter: Counter | None = None,
) -> None:
    """Coxon 2019 Algorithm 3, Lagrange to LCH, in-place."""
    if counter is not None:
        counter.recursive_calls += 1
    if len(phi) != node.dimension or len(indices) != node.size:
        raise ValueError("invalid Coxon Algorithm 3 subproblem")
    if not (0 <= c <= ell <= node.size and b in (0, 1)
            and 1 <= b + c <= node.size):
        raise ValueError("invalid Coxon Algorithm 3 state")
    if node.leaf:
        i0, i1 = indices
        coefficient = phi[0]
        if c == 2:
            values[i1] ^= values[i0]
            if counter is not None:
                counter.loads += 2
                counter.stores += 1
                counter.xors += 1
            _leaf_multiply_add(
                field, values, i0, i1, coefficient, counter
            )
        elif c == 1 and ell == 2 and b == 1:
            product = _fixed_product(field, coefficient, values[i1], counter)
            old0 = values[i0]
            values[i1] ^= old0
            values[i0] = old0 ^ product
            if counter is not None:
                counter.loads += 3
                counter.stores += 2
                counter.xors += 2
        elif c == 1 and ell == 2 and b == 0:
            _leaf_multiply_add(
                field, values, i0, i1, coefficient, counter
            )
        elif c == 0 and ell == 2:
            _leaf_multiply_add(
                field, values, i0, i1, coefficient, counter
            )
        elif c == 1 and ell == 1 and b == 1:
            values[i1] = values[i0]
            if counter is not None:
                counter.loads += 1
                counter.stores += 1
        return

    assert node.alpha is not None and node.delta is not None
    block = 1 << node.d
    c1, c2 = divmod(c, block)
    ell1, ell2 = divmod(ell, block)
    ell2_prime = min(block, ell)
    b_prime = min(b + c2, 1)
    s = min(c2, ell2)
    t = max(c2, ell2)
    mu = tuple(phi[:node.d])
    nu = tuple(phi[node.d:])

    for i in range(max(0, c1 + b_prime - 1)):
        row = indices[block * i:block * (i + 1)]
        coxon_l2x(field, node.alpha, mu, block, block, 0,
                  values, row, counter)
        delta_index = trailing_ones(i)
        mu = _xor_vector(mu, node.sigma_phi_alpha[delta_index], counter)
    if b_prime == 0:
        row = indices[block * (c1 - 1):block * c1]
        coxon_l2x(field, node.alpha, mu, block, block, 0,
                  values, row, counter)

    for j in range(c2, t):
        column = indices[j::block]
        coxon_l2x(field, node.delta, nu, c1, ell1 + 1, b_prime,
                  values, column, counter)
    for j in range(t, ell2_prime):
        column = indices[j::block]
        coxon_l2x(field, node.delta, nu, c1, ell1, b_prime,
                  values, column, counter)
    if b_prime == 1:
        row = indices[block * c1:block * (c1 + 1)]
        coxon_l2x(field, node.alpha, mu, c2, ell2_prime, b,
                  values, row, counter)
    for j in range(0, s):
        column = indices[j::block]
        coxon_l2x(field, node.delta, nu, c1 + 1, ell1 + 1, 0,
                  values, column, counter)
    for j in range(s, c2):
        column = indices[j::block]
        coxon_l2x(field, node.delta, nu, c1 + 1, ell1, 0,
                  values, column, counter)


def coxon_n2x(
    field: BinaryField,
    node: ReductionNode,
    phi: Sequence[int],
    ell: int,
    values: list[int],
    indices: Sequence[int],
    counter: Counter | None = None,
) -> None:
    """Coxon 2019 Algorithm 1, normalised Newton to LCH, in-place."""
    if counter is not None:
        counter.recursive_calls += 1
    if not (1 <= ell <= node.size):
        raise ValueError("invalid Coxon Algorithm 1 length")
    if node.leaf:
        if ell == 2:
            _leaf_multiply_add(
                field, values, indices[0], indices[1], phi[0], counter
            )
        return
    assert node.alpha is not None and node.delta is not None
    block = 1 << node.d
    ell1 = (ell + block - 1) // block - 1
    ell2 = ell - block * ell1
    ell2_prime = min(block, ell)
    mu = tuple(phi[:node.d])
    nu = tuple(phi[node.d:])
    for i in range(ell1):
        row = indices[block * i:block * (i + 1)]
        coxon_n2x(field, node.alpha, mu, block, values, row, counter)
        mu = _xor_vector(
            mu, node.sigma_phi_alpha[trailing_ones(i)], counter
        )
    # The subproblem's storage vector still has block entries.  Only ell2 are
    # semantically live, but Algorithm 1's recursive addressing requires the
    # complete row view.
    full_row = indices[block * ell1:block * (ell1 + 1)]
    coxon_n2x(field, node.alpha, mu, ell2, values, full_row, counter)
    for j in range(ell2):
        coxon_n2x(field, node.delta, nu, ell1 + 1,
                  values, indices[j::block], counter)
    for j in range(ell2, ell2_prime):
        coxon_n2x(field, node.delta, nu, ell1,
                  values, indices[j::block], counter)


@dataclasses.dataclass(frozen=True)
class CoxonPlan:
    parent: int
    count: int
    shift: int
    root: ReductionNode
    root_phi: tuple[int, ...]
    inverse_differences: tuple[tuple[int, ...], ...]
    newton_normalizers: tuple[int, ...]


def build_coxon_plan(
    field: BinaryField, count: int, shift: int,
    counter: Counter | None = None,
    include_newton: bool = True,
) -> CoxonPlan:
    parent = ceil_power_of_two(count)
    dimension = parent.bit_length() - 1
    if dimension == 0:
        # Coxon's tree is defined for a non-empty beta.  A synthetic beta=(1)
        # gives a two-entry root, while count=1 uses only its base case.
        root = build_reduction_tree(field, (1,), counter)
    else:
        root = build_reduction_tree(
            field, tuple(1 << bit for bit in range(dimension)), counter
        )
    root_phi = _phi_for(field, root, shift, counter)
    inverse_differences: list[tuple[int, ...]] = [tuple()]
    normalizers: list[int] = [1]
    if include_newton:
        points = tuple(shift ^ index for index in range(count))
        for order in range(1, count):
            row = []
            for index in range(order, count):
                row.append(field.inverse(points[index] ^ points[index - order]))
                if counter is not None:
                    counter.inversions += 1
            inverse_differences.append(tuple(row))
        for index in range(1, count):
            value = 1
            for earlier in range(index):
                value = field.multiply(value, points[index] ^ points[earlier])
                if counter is not None:
                    counter.multiplications += 1
            normalizers.append(value)
    return CoxonPlan(
        parent, count, shift, root, root_phi,
        tuple(inverse_differences), tuple(normalizers)
    )


def coxon_lagrange_inverse(
    field: BinaryField, plan: CoxonPlan, values: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    if len(values) != plan.count:
        raise ValueError("Coxon Lagrange input length mismatch")
    if plan.count == 1:
        return [values[0]]
    work = list(values) + [0] * (plan.parent - plan.count)
    coxon_l2x(
        field, plan.root, plan.root_phi, plan.count, plan.count, 0,
        work, tuple(range(plan.parent)), counter
    )
    return work[:plan.count] + [0] * (plan.parent - plan.count)


def normalized_newton_coefficients(
    field: BinaryField, plan: CoxonPlan, values: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    if (len(plan.inverse_differences) != plan.count or
            len(plan.newton_normalizers) != plan.count):
        raise ValueError("Coxon plan does not include Newton constants")
    divided = list(values)
    for order in range(1, plan.count):
        inverses = plan.inverse_differences[order]
        for index in range(plan.count - 1, order - 1, -1):
            difference = divided[index] ^ divided[index - 1]
            divided[index] = _fixed_product(
                field, inverses[index - order], difference, counter
            )
            if counter is not None:
                counter.loads += 3
                counter.stores += 1
                counter.xors += 1
    return [
        _fixed_product(field, normalizer, value, counter)
        for normalizer, value in zip(plan.newton_normalizers, divided)
    ]


def coxon_newton_inverse(
    field: BinaryField, plan: CoxonPlan, values: Sequence[int],
    counter: Counter | None = None,
) -> list[int]:
    if plan.count == 1:
        return [values[0]]
    coefficients = normalized_newton_coefficients(
        field, plan, values, counter
    )
    work = coefficients + [0] * (plan.parent - plan.count)
    coxon_n2x(
        field, plan.root, plan.root_phi, plan.count, work,
        tuple(range(plan.parent)), counter
    )
    return work[:plan.count] + [0] * (plan.parent - plan.count)


@dataclasses.dataclass(frozen=True)
class TangHanPlan:
    parent: int
    count: int
    shift: int
    split_multipliers: tuple[int, ...]


def build_tang_han_plan(
    field: BinaryField, count: int, shift: int,
    counter: Counter | None = None,
) -> TangHanPlan:
    parent = ceil_power_of_two(count)
    multipliers = []
    local_shift = shift
    local_count = count
    local_parent = parent
    while local_parent > 1:
        half = local_parent // 2
        multipliers.append(subspace_multiplier(
            field, half.bit_length() - 1, local_shift, counter
        ))
        if local_count > half:
            local_count -= half
            local_shift ^= half
        local_parent = half
    return TangHanPlan(parent, count, shift, tuple(multipliers))


def _tang_han_recursive(
    field: BinaryField,
    values: Sequence[int],
    count: int,
    parent: int,
    shift: int,
    counter: Counter | None,
) -> tuple[list[int], list[int]]:
    """Tang--Han 2022 Appendix-A Algorithm 8."""
    if counter is not None:
        counter.recursive_calls += 1
    if parent == 1:
        return [values[0]], [values[0]]
    half = parent // 2
    dimension = half.bit_length() - 1
    if count <= half:
        coefficients0, completed0 = _tang_han_recursive(
            field, values, count, half, shift, counter
        )
        completed1 = lch_forward(
            field, coefficients0, shift ^ half, counter
        )
        return coefficients0 + [0] * half, completed0 + completed1

    # Paper Algorithm 8 lines 9--12.
    w = lch_inverse(field, values[:half], shift, counter)
    w_prime = lch_forward(field, w, shift ^ half, counter)
    tail_count = count - half
    adjusted = [values[half + index] ^ w_prime[index]
                for index in range(tail_count)]
    if counter is not None:
        counter.loads += 2 * tail_count
        counter.stores += tail_count
        counter.xors += tail_count
    coefficients1, adjusted_completed1 = _tang_han_recursive(
        field, adjusted, tail_count, half, shift ^ half, counter
    )
    multiplier = subspace_multiplier(field, dimension, shift)
    coefficients0 = [
        low ^ _fixed_product(field, multiplier, high, counter)
        for low, high in zip(w, coefficients1)
    ]
    completed1 = [left ^ right for left, right
                  in zip(adjusted_completed1, w_prime)]
    if counter is not None:
        counter.loads += 4 * half
        counter.stores += 2 * half
        counter.xors += 2 * half
    return coefficients0 + coefficients1, list(values[:half]) + completed1


def tang_han_inverse(
    field: BinaryField, plan: TangHanPlan, values: Sequence[int],
    counter: Counter | None = None,
) -> tuple[list[int], list[int]]:
    if len(values) != plan.count:
        raise ValueError("Tang--Han input length mismatch")
    return _tang_han_recursive(
        field, values, plan.count, plan.parent, plan.shift, counter
    )


def stable_seed(field_name: str, count: int, shift: int) -> int:
    digest = hashlib.sha256(
        f"leopard2-c3b:{field_name}:{count}:{shift}".encode("ascii")
    ).digest()
    return int.from_bytes(digest[:8], "little")


def value_vectors(
    field: BinaryField, count: int, shift: int
) -> Iterable[tuple[int, ...]]:
    if field.name == "gf4" and count <= 3:
        for encoded in range(field.order ** count):
            work = encoded
            output = []
            for _ in range(count):
                output.append(work % field.order)
                work //= field.order
            yield tuple(output)
        return
    yield (0,) * count
    yield (1,) * count
    for index in range(count):
        yield tuple(int(column == index) for column in range(count))
    yield tuple((column * 29 + count * 7 + shift) % field.order
                for column in range(count))
    generator = random.Random(stable_seed(field.name, count, shift))
    for _ in range(6):
        yield tuple(generator.randrange(field.order) for _ in range(count))


def shifts_for(field: BinaryField, parent: int) -> tuple[int, ...]:
    candidates = (0, 1, 3, field.order // 2 + 1, field.order - 1)
    return tuple(sorted(set(value for value in candidates
                            if 0 <= value < field.order)))


def correctness_specs() -> list[tuple[str, int, int]]:
    counts = {
        "gf4": tuple(range(1, 17)),
        "gf8": (1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31),
        "gf16": (1, 3, 5, 9, 17),
    }
    output = []
    for field_name, selections in counts.items():
        field = field_named(field_name)
        for count in selections:
            parent = ceil_power_of_two(count)
            for shift in shifts_for(field, parent):
                output.append((field_name, count, shift))
    return output


def correctness_job(specification: tuple[str, int, int]) -> dict[str, object]:
    field_name, count, shift = specification
    field = field_named(field_name)
    parent = ceil_power_of_two(count)
    dense_plan = build_dense_plan(field, count, shift)
    generic_plan = build_generic_plan(field, count, shift)
    coxon_plan = build_coxon_plan(field, count, shift)
    tang_plan = build_tang_han_plan(field, count, shift)
    points = tuple(shift ^ index for index in range(count))
    parent_points = tuple(shift ^ index for index in range(parent))
    digest = hashlib.sha256()
    vectors = 0
    comparisons = 0
    padded_exact_differences = 0
    exact_suffix_nonzero = 0
    for vector in value_vectors(field, count, shift):
        values = list(vector)
        padded = lch_inverse(
            field, values + [0] * (parent - count), shift
        )
        dense = dense_inverse(field, dense_plan, values)
        generic = generic_inverse(field, generic_plan, values)
        coxon_lagrange = coxon_lagrange_inverse(field, coxon_plan, values)
        coxon_newton = coxon_newton_inverse(field, coxon_plan, values)
        tang, completed = tang_han_inverse(field, tang_plan, values)
        if not (dense == generic == coxon_lagrange == coxon_newton == tang):
            raise AssertionError(
                f"exact methods differ for {field_name} K={count} shift={shift}"
            )
        if count < parent:
            scratch = list(values) + [
                (field.order - 1 - index) % field.order
                for index in range(parent - count)
            ]
            coxon_l2x(
                field, coxon_plan.root, coxon_plan.root_phi,
                count, count, 0, scratch, tuple(range(parent))
            )
            if scratch[:count] + [0] * (parent - count) != dense:
                raise AssertionError(
                    "Coxon Algorithm 3 depends on unspecified scratch"
                )
        if evaluate_lch(field, dense, parent, points) != values:
            raise AssertionError("direct exact evaluation failed")
        if lch_forward(field, dense, shift) != completed:
            raise AssertionError("Tang--Han completed evaluations are wrong")
        if completed[:count] != values:
            raise AssertionError("Tang--Han changed supplied evaluations")
        if evaluate_lch(field, padded, parent, parent_points) != (
            values + [0] * (parent - count)
        ):
            raise AssertionError("wire-compatible padded inverse failed")
        if count == parent:
            if dense != padded:
                raise AssertionError("full-parent exact and padded inverses differ")
        else:
            padded_exact_differences += int(dense != padded)
            exact_suffix_nonzero += int(any(completed[count:]))
        for result in (padded, dense, generic, coxon_lagrange,
                       coxon_newton, tang, completed):
            for item in result:
                digest.update(int(item).to_bytes(2, "little"))
        vectors += 1
        comparisons += 7 * parent + 3 * count
    return {
        "comparisons": comparisons,
        "count": count,
        "digest": digest.hexdigest(),
        "exact_suffix_nonzero_vectors": exact_suffix_nonzero,
        "field": field_name,
        "padded_exact_difference_vectors": padded_exact_differences,
        "parent": parent,
        "shift": shift,
        "vectors": vectors,
    }


def method_class(method: str) -> str:
    return (
        "wire-compatible-parent" if method == METHOD_PADDED
        else "new-exact-profile"
    )


def operation_case(field_name: str, count: int, shift: int) -> dict[str, object]:
    field = field_named(field_name)
    parent = ceil_power_of_two(count)
    values = [((index + 1) * 29 + shift) % field.order
              for index in range(count)]
    setup: dict[str, Counter] = {METHOD_PADDED: Counter()}
    execution: dict[str, Counter] = {}

    dense_setup = Counter()
    dense_plan = build_dense_plan(field, count, shift, dense_setup)
    setup[METHOD_DENSE] = dense_setup
    generic_setup = Counter()
    generic_plan = build_generic_plan(field, count, shift, generic_setup)
    setup[METHOD_GENERIC] = generic_setup
    coxon_l_setup = Counter()
    coxon_l_plan = build_coxon_plan(
        field, count, shift, coxon_l_setup, include_newton=False
    )
    setup[METHOD_COXON_L] = coxon_l_setup
    coxon_n_setup = Counter()
    coxon_n_plan = build_coxon_plan(
        field, count, shift, coxon_n_setup, include_newton=True
    )
    setup[METHOD_COXON_N] = coxon_n_setup
    tang_setup = Counter()
    tang_plan = build_tang_han_plan(field, count, shift, tang_setup)
    setup[METHOD_TANG] = tang_setup

    padded_counter = Counter()
    lch_inverse(field, values + [0] * (parent - count), shift, padded_counter)
    execution[METHOD_PADDED] = padded_counter
    dense_counter = Counter()
    dense_inverse(field, dense_plan, values, dense_counter)
    execution[METHOD_DENSE] = dense_counter
    generic_counter = Counter()
    generic_inverse(field, generic_plan, values, generic_counter)
    execution[METHOD_GENERIC] = generic_counter
    coxon_l_counter = Counter()
    coxon_lagrange_inverse(field, coxon_l_plan, values, coxon_l_counter)
    execution[METHOD_COXON_L] = coxon_l_counter
    coxon_n_counter = Counter()
    coxon_newton_inverse(field, coxon_n_plan, values, coxon_n_counter)
    execution[METHOD_COXON_N] = coxon_n_counter
    tang_counter = Counter()
    tang_han_inverse(field, tang_plan, values, tang_counter)
    execution[METHOD_TANG] = tang_counter

    return {
        "count": count,
        "field": field_name,
        "methods": [
            {
                "class": method_class(method),
                "execution": dataclasses.asdict(execution[method]),
                "method": method,
                "setup": dataclasses.asdict(setup[method]),
            }
            for method in METHODS
        ],
        "parent": parent,
        "shift": shift,
    }


def counter_score(counter: dict[str, int]) -> int:
    return (
        3 * counter["multiplications"] + counter["xors"] +
        counter["loads"] + 2 * counter["stores"] +
        12 * counter["inversions"]
    )


def modeled_cells(cases: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for case in cases:
        methods = {row["method"]: row for row in case["methods"]}
        bytes_per_symbol = 1 if case["field"] == "gf8" else 2
        for shard_bytes in (64, 1024, 65536):
            symbols = max(1, shard_bytes // bytes_per_symbol)
            for batch in (1, 8, 64):
                for reuse in (1, 8, 64, 1024):
                    for method in METHODS:
                        row = methods[method]
                        execution = counter_score(row["execution"])
                        setup = counter_score(row["setup"])
                        score = execution * symbols * batch + setup / reuse
                        baseline = methods[METHOD_PADDED]
                        baseline_score = (
                            counter_score(baseline["execution"]) *
                            symbols * batch +
                            counter_score(baseline["setup"]) / reuse
                        )
                        output.append({
                            "batch": batch,
                            "class": row["class"],
                            "count": case["count"],
                            "field": case["field"],
                            "gain_over_padded": round(
                                baseline_score / score, 9
                            ) if score else 1.0,
                            "method": method,
                            "parent": case["parent"],
                            "reuse": reuse,
                            "score": round(score, 3),
                            "shard_bytes": shard_bytes,
                            "shift": case["shift"],
                        })
    return output


def source_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def deterministic_analysis(workers: int = 1) -> dict[str, object]:
    if not (1 <= workers <= 8):
        raise ValueError("C3b workers must be in 1..8")
    specifications = correctness_specs()
    if workers == 1:
        correctness = [correctness_job(specification)
                       for specification in specifications]
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=workers
        ) as executor:
            correctness = list(executor.map(
                correctness_job, specifications
            ))
    correctness.sort(key=lambda row: (
        str(row["field"]), int(row["count"]), int(row["shift"])
    ))
    cases = []
    for field_name, counts in (
        ("gf8", (3, 5, 7, 9, 15, 17, 31)),
        ("gf16", (3, 5, 9, 17)),
    ):
        for count in counts:
            cases.append(operation_case(field_name, count, 3))
    cells = modeled_cells(cases)
    exact_cells = [cell for cell in cells if cell["class"] == "new-exact-profile"]
    best_exact = max(exact_cells, key=lambda row: float(row["gain_over_padded"]))
    by_method = {}
    for method in METHODS[1:]:
        selections = [cell for cell in exact_cells if cell["method"] == method]
        by_method[method] = max(
            selections, key=lambda row: float(row["gain_over_padded"])
        )
    summary = {
        "comparisons": sum(int(row["comparisons"]) for row in correctness),
        "correctness_jobs": len(correctness),
        "exact_suffix_nonzero_vectors": sum(
            int(row["exact_suffix_nonzero_vectors"]) for row in correctness
        ),
        "modeled_cells": len(cells),
        "padded_exact_difference_vectors": sum(
            int(row["padded_exact_difference_vectors"]) for row in correctness
        ),
        "vectors": sum(int(row["vectors"]) for row in correctness),
    }
    return {
        "best_modeled_exact": best_exact,
        "best_modeled_exact_by_method": by_method,
        "correctness": correctness,
        "disposition": {
            "exact_profile_identity": "required-before-any-production-use",
            "production_promotion": "none",
            "wire_compatible_candidate": "none",
        },
        "method_classes": {method: method_class(method) for method in METHODS},
        "modeled_cells": cells,
        "operation_cases": cases,
        "schema": SCHEMA,
        "source_sha256": source_sha256(),
        "summary": summary,
    }


def canonical_write(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, separators=(",", ": ")) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def run_analysis(arguments: argparse.Namespace) -> int:
    result = deterministic_analysis(arguments.workers)
    canonical_write(Path(arguments.output), result)
    print(json.dumps(result["summary"], sort_keys=True))
    return 0


def verify_result(arguments: argparse.Namespace) -> int:
    actual = json.loads(Path(arguments.result).read_text(encoding="utf-8"))
    if actual.get("schema") != SCHEMA:
        raise ValueError("C3b schema mismatch")
    if actual.get("source_sha256") != source_sha256():
        raise ValueError("C3b result is stale for its source")
    expected = deterministic_analysis(1)
    if actual != expected:
        raise ValueError("C3b result differs from deterministic recomputation")
    print("leopard2 C3b deterministic result verified")
    return 0


def median_mad(values: Sequence[int]) -> tuple[int, int]:
    median = int(statistics.median(values))
    mad = int(statistics.median(abs(item - median) for item in values))
    return median, mad


def benchmark_case(
    field: BinaryField, count: int, shift: int, samples: int, iterations: int
) -> dict[str, object]:
    parent = ceil_power_of_two(count)
    generator = random.Random(stable_seed(field.name, count, shift))
    values = [generator.randrange(field.order) for _ in range(count)]
    padded_values = values + [0] * (parent - count)
    builders: dict[str, Callable[[], object]] = {
        METHOD_DENSE: lambda: build_dense_plan(field, count, shift),
        METHOD_GENERIC: lambda: build_generic_plan(field, count, shift),
        METHOD_COXON_L: lambda: build_coxon_plan(
            field, count, shift, include_newton=False
        ),
        METHOD_COXON_N: lambda: build_coxon_plan(
            field, count, shift, include_newton=True
        ),
        METHOD_TANG: lambda: build_tang_han_plan(field, count, shift),
    }
    plans: dict[str, object] = {}
    setup = {METHOD_PADDED: {"mad_ns": 0, "median_ns": 0}}
    for method, builder in builders.items():
        times = []
        for _ in range(samples):
            start = time.perf_counter_ns()
            plan = builder()
            times.append(time.perf_counter_ns() - start)
        plans[method] = plan
        median, mad = median_mad(times)
        setup[method] = {"mad_ns": mad, "median_ns": median}

    functions: dict[str, Callable[[], object]] = {
        METHOD_PADDED: lambda: lch_inverse(field, padded_values, shift),
        METHOD_DENSE: lambda: dense_inverse(
            field, plans[METHOD_DENSE], values
        ),
        METHOD_GENERIC: lambda: generic_inverse(
            field, plans[METHOD_GENERIC], values
        ),
        METHOD_COXON_L: lambda: coxon_lagrange_inverse(
            field, plans[METHOD_COXON_L], values
        ),
        METHOD_COXON_N: lambda: coxon_newton_inverse(
            field, plans[METHOD_COXON_N], values
        ),
        METHOD_TANG: lambda: tang_han_inverse(
            field, plans[METHOD_TANG], values
        )[0],
    }
    execution = {}
    reference = functions[METHOD_DENSE]()
    for method, function in functions.items():
        result = function()
        if method == METHOD_PADDED:
            if count == parent and result != reference:
                raise AssertionError("timing full-parent result mismatch")
        elif result != reference:
            raise AssertionError(f"timing result mismatch for {method}")
        times = []
        for _ in range(samples):
            start = time.perf_counter_ns()
            for _ in range(iterations):
                function()
            times.append((time.perf_counter_ns() - start) // iterations)
        median, mad = median_mad(times)
        execution[method] = {"mad_ns": mad, "median_ns": median}
    return {
        "count": count,
        "execution_per_symbol_position": execution,
        "field": field.name,
        "iterations_per_sample": iterations,
        "parent": parent,
        "samples": samples,
        "setup": setup,
        "shift": shift,
    }


def current_affinity() -> list[int]:
    if hasattr(os, "sched_getaffinity"):
        return sorted(os.sched_getaffinity(0))
    return []


def run_benchmark(arguments: argparse.Namespace) -> int:
    if arguments.samples < 3 or arguments.iterations <= 0:
        raise ValueError("benchmark requires samples>=3 and iterations>0")
    before = current_affinity()
    if len(before) > 1:
        raise ValueError("benchmark must be pinned to exactly one allowed CPU")
    cases = []
    for field_name in arguments.fields:
        field = field_named(field_name)
        for count in arguments.counts:
            if ceil_power_of_two(count) <= field.order:
                cases.append(benchmark_case(
                    field, count, arguments.shift,
                    arguments.samples, arguments.iterations
                ))
    if current_affinity() != before:
        raise RuntimeError("CPU affinity changed during benchmark")
    result = {
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
    canonical_write(Path(arguments.output), result)
    print(f"wrote {arguments.output}")
    return 0


def summarize_timing(arguments: argparse.Namespace) -> int:
    timing = json.loads(Path(arguments.timing).read_text(encoding="utf-8"))
    if timing.get("schema") != TIMING_SCHEMA:
        raise ValueError("timing schema mismatch")
    if timing.get("source_sha256") != source_sha256():
        raise ValueError("timing is stale for the source")
    cells = []
    for case in timing["cases"]:
        bytes_per_symbol = 1 if case["field"] == "gf8" else 2
        for shard_bytes in (64, 1024, 65536):
            symbols = max(1, shard_bytes // bytes_per_symbol)
            for batch in (1, 8, 64):
                for reuse in (1, 8, 64, 1024):
                    padded_ns = (
                        case["execution_per_symbol_position"][METHOD_PADDED]["median_ns"] *
                        symbols * batch
                    )
                    for method in METHODS:
                        total_ns = (
                            case["execution_per_symbol_position"][method]["median_ns"] *
                            symbols * batch +
                            case["setup"][method]["median_ns"] / reuse
                        )
                        cells.append({
                            "batch": batch,
                            "class": method_class(method),
                            "count": case["count"],
                            "field": case["field"],
                            "gain_over_padded": round(padded_ns / total_ns, 9),
                            "method": method,
                            "parent": case["parent"],
                            "reuse": reuse,
                            "shard_bytes": shard_bytes,
                            "total_ns_extrapolated": round(total_ns, 3),
                        })
    exact = [cell for cell in cells if cell["class"] == "new-exact-profile"]
    best_by_method = {}
    for method in METHODS[1:]:
        selection = [cell for cell in exact if cell["method"] == method]
        best_by_method[method] = max(
            selection, key=lambda row: float(row["gain_over_padded"])
        )
    result = {
        "best_exact_by_method": best_by_method,
        "cells": cells,
        "disposition": {
            "production_promotion": "none",
            "reason": (
                "Python scalar exact-profile timing cannot promote a wire-incompatible "
                "profile; retain candidates for C7/C8 C++ whole-codec evaluation"
            ),
            "threshold": "credible >=10% whole-codec gain in a meaningful region",
        },
        "schema": "leopard2-c3b-timing-summary/v1",
        "source_sha256": source_sha256(),
        "timing_sha256": hashlib.sha256(
            Path(arguments.timing).read_bytes()
        ).hexdigest(),
    }
    canonical_write(Path(arguments.output), result)
    print(f"wrote {arguments.output}")
    return 0


def self_test() -> int:
    for field_name in ("gf4", "gf8", "gf16"):
        field = field_named(field_name)
        probes = (0, 1, 2, 3, field.order - 1)
        for left in probes:
            for right in probes:
                if field.multiply(left, right) != field.multiply(right, left):
                    raise AssertionError("field multiplication is not commutative")
        for value in probes[1:]:
            if field.multiply(value, field.inverse(value)) != 1:
                raise AssertionError("field inverse is wrong")
        for count in (1, 2, 3, 5):
            if ceil_power_of_two(count) > field.order:
                continue
            for shift in (0, 1, field.order - 1):
                result = correctness_job((field_name, count, shift))
                if int(result["vectors"]) <= 0:
                    raise AssertionError("empty self-test")
    print("leopard2 C3b fast inverse self-test passed")
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser("self-test")
    analyze = commands.add_parser("analyze")
    analyze.add_argument("--output", required=True)
    analyze.add_argument("--workers", type=int, default=8)
    verify = commands.add_parser("verify")
    verify.add_argument("result")
    benchmark = commands.add_parser("benchmark")
    benchmark.add_argument("--fields", nargs="+", choices=("gf8", "gf16"),
                           default=("gf8", "gf16"))
    benchmark.add_argument("--counts", nargs="+", type=int,
                           default=(3, 5, 9, 17))
    benchmark.add_argument("--shift", type=int, default=3)
    benchmark.add_argument("--samples", type=int, default=9)
    benchmark.add_argument("--iterations", type=int, default=16)
    benchmark.add_argument("--output", required=True)
    summarize = commands.add_parser("summarize-timing")
    summarize.add_argument("timing")
    summarize.add_argument("--output", required=True)
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
        if arguments.command == "summarize-timing":
            return summarize_timing(arguments)
        raise AssertionError("unknown command")
    except (AssertionError, OSError, ValueError, ZeroDivisionError) as error:
        print(f"C3b error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
