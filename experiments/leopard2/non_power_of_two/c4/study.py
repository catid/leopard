#!/usr/bin/env python3
"""C4 full-dyadic-block plus tail scalar experiment.

This program deliberately studies two different maps:

* ``parent_*`` methods compute the existing power-of-two parent inverse with
  the untransmitted suffix fixed to zero.  They must match the padded Leopard
  parent byte-for-byte.
* ``exact_*`` methods interpolate a degree-<q polynomial from exactly q
  points.  For non-power-of-two q this is a new code definition, even when a
  full 2^a transform is reused for the head block.

The direct, smaller-padded, and recursive prefix tail executors are scalar
research models.  The Newton tail is an exact-profile oracle, not a claim to
implement Coxon's fast recursive conversion or Tang-Han epsilon transforms.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import functools
import hashlib
import importlib.util
import itertools
import json
import os
import platform
import random
import statistics
import sys
import time
from pathlib import Path
from typing import Iterable, Sequence


SCHEMA = "leopard2-c4-full-block-tail/v1"
TIMING_SCHEMA = "leopard2-c4-full-block-tail-timing/v1"
ROOT = Path(__file__).resolve().parents[4]
C3_PATH = ROOT / "experiments/leopard2/non_power_of_two/c3/study.py"


def _load_c3():
    name = "_leopard2_c3_reference"
    existing = sys.modules.get(name)
    if existing is not None:
        return existing
    specification = importlib.util.spec_from_file_location(name, C3_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load the C3 scalar reference")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


c3 = _load_c3()


@dataclasses.dataclass
class Counter:
    xors: int = 0
    multiplications: int = 0
    inversions: int = 0
    loads: int = 0
    stores: int = 0

    def add(self, other: "Counter", scale: int = 1) -> None:
        self.xors += other.xors * scale
        self.multiplications += other.multiplications * scale
        self.inversions += other.inversions * scale
        self.loads += other.loads * scale
        self.stores += other.stores * scale


def counter_from_c3(value) -> Counter:
    return Counter(
        value.xors, value.multiplications, value.inversions,
        value.loads, value.stores,
    )


def ceil_power_of_two(value: int) -> int:
    if value <= 0:
        raise ValueError("length must be positive")
    return 1 << (value - 1).bit_length()


def split_geometry(q: int) -> tuple[int, int, int]:
    if q < 3 or q & (q - 1) == 0:
        raise ValueError("q must be a non-power-of-two count of at least three")
    block = 1 << (q.bit_length() - 1)
    tail = q - block
    return block, tail, 2 * block


def matrix_multiply(field, left: Sequence[Sequence[int]],
                    right: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    if not left or not right or any(len(row) != len(right) for row in left):
        raise ValueError("matrix dimensions differ")
    columns = len(right[0])
    if any(len(row) != columns for row in right):
        raise ValueError("ragged matrix")
    result: list[tuple[int, ...]] = []
    for left_row in left:
        row: list[int] = []
        for column in range(columns):
            value = 0
            for index, coefficient in enumerate(left_row):
                other = right[index][column]
                if coefficient and other:
                    value ^= (other if coefficient == 1 else
                              field.multiply(coefficient, other))
            row.append(value)
        result.append(tuple(row))
    return tuple(result)


def matrix_xor(left: Sequence[Sequence[int]],
               right: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    if len(left) != len(right) or any(len(a) != len(b)
                                      for a, b in zip(left, right)):
        raise ValueError("matrix dimensions differ")
    return tuple(tuple(a ^ b for a, b in zip(left_row, right_row))
                 for left_row, right_row in zip(left, right))


def apply(field, matrix: Sequence[Sequence[int]], values: Sequence[int],
          counter=None) -> list[int]:
    return c3.apply_matrix(field, matrix, values, counter)


def combine_inverse(field, left: Sequence[int], right: Sequence[int],
                    shift: int, counter=None) -> list[int]:
    if len(left) != len(right) or not left:
        raise ValueError("inverse children must have equal nonzero size")
    half = len(left)
    multiplier = c3.fft_skew(field.name)[shift + half - 1]
    output = [0] * (2 * half)
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
class DirectTailPlan:
    block: int
    tail: int
    shift: int
    columns: tuple[tuple[int, ...], ...]
    construction: str


def build_direct_tail_plan(field, block: int, tail: int, shift: int,
                           independent_limit: int = 32) -> DirectTailPlan:
    if not 0 < tail < block:
        raise ValueError("tail must be strictly smaller than its block")
    if block <= independent_limit:
        matrix = c3.lch_evaluation_matrix(field, block, shift, block)
        inverse = c3.invert_matrix(field, matrix)
        columns = tuple(tuple(row[column] for column in range(tail))
                        for row in inverse)
        construction = "direct-polynomial-matrix-inverse"
    else:
        columns_by_input: list[list[int]] = []
        for column in range(tail):
            unit = [0] * block
            unit[column] = 1
            columns_by_input.append(c3.padded_lch_inverse(field, unit, shift))
        columns = tuple(tuple(columns_by_input[column][row]
                              for column in range(tail))
                        for row in range(block))
        construction = "unit-columns-from-padded-reference"
    return DirectTailPlan(block, tail, shift, columns, construction)


def direct_tail(field, plan: DirectTailPlan, values: Sequence[int],
                counter=None) -> list[int]:
    if len(values) != plan.tail:
        raise ValueError("direct tail value count differs from plan")
    return apply(field, plan.columns, values, counter)


def smaller_padded_tail(field, values: Sequence[int], block: int, shift: int,
                        counter=None) -> list[int]:
    tail = len(values)
    if not 0 < tail < block:
        raise ValueError("invalid smaller-padded tail")
    small = ceil_power_of_two(tail)
    output = c3.padded_lch_inverse(
        field, list(values) + [0] * (small - tail), shift, counter
    )
    size = small
    while size < block:
        output = combine_inverse(field, output, [0] * size, shift, counter)
        size *= 2
    return output


def recursive_prefix_tail(field, values: Sequence[int], block: int, shift: int,
                          counter=None) -> list[int]:
    """C2-style prefix recursion specialized to all requested outputs."""
    active = len(values)
    if not 0 < active <= block or block & (block - 1):
        raise ValueError("invalid recursive tail geometry")
    if active == block:
        return c3.padded_lch_inverse(field, values, shift, counter)
    half = block // 2
    if active <= half:
        left = recursive_prefix_tail(field, values, half, shift, counter)
        right = [0] * half
    else:
        left = c3.padded_lch_inverse(field, values[:half], shift, counter)
        right = recursive_prefix_tail(
            field, values[half:], half, shift + half, counter
        )
    return combine_inverse(field, left, right, shift, counter)


def parent_full_block(field, values: Sequence[int], shift: int, method: str,
                      direct_plan: DirectTailPlan | None = None,
                      counter=None) -> list[int]:
    block, tail, parent = split_geometry(len(values))
    c3.validate_geometry(field, parent, shift)
    head = c3.padded_lch_inverse(field, values[:block], shift, counter)
    tail_shift = shift + block
    if method == "direct":
        if direct_plan is None:
            direct_plan = build_direct_tail_plan(field, block, tail, tail_shift)
        tail_coefficients = direct_tail(field, direct_plan, values[block:], counter)
    elif method == "smaller_padded":
        tail_coefficients = smaller_padded_tail(
            field, values[block:], block, tail_shift, counter
        )
    elif method == "recursive":
        tail_coefficients = recursive_prefix_tail(
            field, values[block:], block, tail_shift, counter
        )
    else:
        raise ValueError("unknown parent tail method")
    return combine_inverse(field, head, tail_coefficients, shift, counter)


@dataclasses.dataclass(frozen=True)
class ExactSplitPlan:
    q: int
    block: int
    tail: int
    parent: int
    shift: int
    head_to_tail: tuple[tuple[int, ...], ...]
    head_adjustment: tuple[tuple[int, ...], ...]
    tail_inverse: tuple[tuple[int, ...], ...]


def build_exact_split_plan(field, q: int, shift: int) -> ExactSplitPlan:
    block, tail, parent = split_geometry(q)
    c3.validate_geometry(field, parent, shift)
    matrix = c3.lch_evaluation_matrix(field, q, shift, parent)
    a00 = tuple(tuple(row[:block]) for row in matrix[:block])
    a01 = tuple(tuple(row[block:]) for row in matrix[:block])
    a10 = tuple(tuple(row[:block]) for row in matrix[block:])
    a11 = tuple(tuple(row[block:]) for row in matrix[block:])
    a00_inverse = c3.invert_matrix(field, a00)
    head_adjustment = matrix_multiply(field, a00_inverse, a01)
    schur = matrix_xor(
        a11, matrix_multiply(field, a10, head_adjustment)
    )
    return ExactSplitPlan(
        q, block, tail, parent, shift, a10, head_adjustment,
        c3.invert_matrix(field, schur),
    )


def exact_split_direct(field, plan: ExactSplitPlan, values: Sequence[int],
                       counter=None) -> list[int]:
    if len(values) != plan.q:
        raise ValueError("exact split value count differs from plan")
    head = c3.padded_lch_inverse(
        field, values[:plan.block], plan.shift, counter
    )
    predicted_tail = apply(field, plan.head_to_tail, head, counter)
    residual = [value ^ predicted for value, predicted
                in zip(values[plan.block:], predicted_tail)]
    if counter is not None:
        counter.loads += 2 * plan.tail
        counter.stores += plan.tail
        counter.xors += plan.tail
    tail = apply(field, plan.tail_inverse, residual, counter)
    if any(any(row) for row in plan.head_adjustment):
        adjustment = apply(field, plan.head_adjustment, tail, counter)
        head = [value ^ delta for value, delta in zip(head, adjustment)]
        if counter is not None:
            counter.loads += 2 * plan.block
            counter.stores += plan.block
            counter.xors += plan.block
    return head + tail


@dataclasses.dataclass(frozen=True)
class ExactNewtonTailPlan:
    split: ExactSplitPlan
    newton: object


def build_exact_newton_tail_plan(field, q: int) -> ExactNewtonTailPlan:
    split = build_exact_split_plan(field, q, 0)
    if any(any(row) for row in split.head_adjustment):
        raise AssertionError("zero-coset exact split unexpectedly needs head adjustment")
    newton = c3.build_newton_plan(field, split.tail, split.block)
    expected = c3.lch_evaluation_matrix(
        field, split.tail, split.block, ceil_power_of_two(split.tail)
    )
    if split.tail_inverse != c3.invert_matrix(field, expected):
        raise AssertionError("Newton tail basis does not match the exact Schur basis")
    return ExactNewtonTailPlan(split, newton)


def exact_split_newton(field, plan: ExactNewtonTailPlan,
                       values: Sequence[int], counter=None) -> list[int]:
    split = plan.split
    if len(values) != split.q:
        raise ValueError("exact Newton value count differs from plan")
    head = c3.padded_lch_inverse(field, values[:split.block], 0, counter)
    predicted_tail = apply(field, split.head_to_tail, head, counter)
    residual = [value ^ predicted for value, predicted
                in zip(values[split.block:], predicted_tail)]
    if counter is not None:
        counter.loads += 2 * split.tail
        counter.stores += split.tail
        counter.xors += split.tail
    tail = c3.exact_newton_to_lch(field, plan.newton, residual, counter)
    return head + tail


def stable_seed(label: str, field_name: str, q: int, shift: int) -> int:
    digest = hashlib.sha256(
        f"c4:{label}:{field_name}:{q}:{shift}".encode("ascii")
    ).digest()
    return int.from_bytes(digest[:8], "little")


def value_vectors(field, q: int, shift: int) -> Iterable[tuple[int, ...]]:
    yield (0,) * q
    if field.name == "gf4" and q <= 3:
        yield from itertools.product(range(field.order), repeat=q)
        return
    for index in range(min(q, 10)):
        value = [0] * q
        value[index] = 1
        yield tuple(value)
    yield tuple(1 for _ in range(q))
    yield tuple((index * 29 + q * 7 + shift) % field.order
                for index in range(q))
    random_source = random.Random(stable_seed("vectors", field.name, q, shift))
    for _ in range(12):
        yield tuple(random_source.randrange(field.order) for _ in range(q))


def correctness_job(specification: tuple[str, int, int]) -> dict[str, object]:
    field_name, q, shift = specification
    field = c3.field_named(field_name)
    block, tail, parent = split_geometry(q)
    direct_plan = build_direct_tail_plan(field, block, tail, shift + block)
    exact_checked = block <= 32 and tail <= 17
    exact_plan = build_exact_split_plan(field, q, shift) if exact_checked else None
    exact_oracle = c3.build_exact_dense_plan(field, q, shift) if exact_checked else None
    newton_plan = (build_exact_newton_tail_plan(field, q)
                   if exact_checked and shift == 0 else None)
    vectors = comparisons = parent_exact_differences = 0
    for values in value_vectors(field, q, shift):
        vectors += 1
        padded = c3.padded_lch_inverse(
            field, list(values) + [0] * (parent - q), shift
        )
        for method in ("direct", "smaller_padded", "recursive"):
            candidate = parent_full_block(
                field, values, shift, method, direct_plan
            )
            if candidate != padded:
                raise AssertionError(
                    f"{field_name} q={q} shift={shift} {method} changed parent"
                )
            comparisons += parent
        if exact_checked:
            assert exact_plan is not None and exact_oracle is not None
            exact = exact_split_direct(field, exact_plan, values)
            oracle = apply(field, exact_oracle.inverse, values)
            if exact != oracle:
                raise AssertionError("full-block exact direct differs from dense oracle")
            points = tuple(shift ^ index for index in range(q))
            if c3.evaluate_lch(field, exact, parent, points) != list(values):
                raise AssertionError("exact split does not interpolate its inputs")
            comparisons += 2 * q
            if exact != padded[:q] or any(padded[q:]):
                parent_exact_differences += 1
            if newton_plan is not None:
                newton = exact_split_newton(field, newton_plan, values)
                if newton != exact:
                    raise AssertionError("Newton tail differs from direct exact tail")
                comparisons += q
    return {
        "comparisons": comparisons,
        "direct_plan_construction": direct_plan.construction,
        "exact_checked": exact_checked,
        "field": field_name,
        "parent": parent,
        "parent_exact_difference_vectors": parent_exact_differences,
        "q": q,
        "shift": shift,
        "tail": tail,
        "vectors": vectors,
    }


def rank(field, matrix: Sequence[Sequence[int]]) -> int:
    if not matrix:
        return 0
    work = [list(row) for row in matrix]
    columns = len(work[0])
    row = 0
    for column in range(columns):
        pivot = next((index for index in range(row, len(work))
                      if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = field.inverse(work[row][column])
        work[row] = [field.multiply(value, inverse) for value in work[row]]
        for other in range(len(work)):
            if other == row or work[other][column] == 0:
                continue
            factor = work[other][column]
            work[other] = [
                value ^ field.multiply(factor, source)
                for value, source in zip(work[other], work[row])
            ]
        row += 1
        if row == len(work):
            break
    return row


def mds_job(q: int) -> dict[str, int]:
    field = c3.field_named("gf4")
    full_parent = 16
    evaluation = c3.lch_evaluation_matrix(field, full_parent, 0, full_parent)
    systematic = tuple(tuple(row[:q]) for row in evaluation[:q])
    inverse = c3.invert_matrix(field, systematic)
    generator = matrix_multiply(
        field, tuple(tuple(row[:q]) for row in evaluation), inverse
    )
    identity = tuple(tuple(int(row == column) for column in range(q))
                     for row in range(q))
    if generator[:q] != identity:
        raise AssertionError("exact profile generator is not systematic")
    subsets = 0
    for coordinates in itertools.combinations(range(full_parent), q):
        subsets += 1
        if rank(field, [generator[index] for index in coordinates]) != q:
            raise AssertionError(f"GF(2^4) exact q={q} is not MDS")
    return {"field_order": 16, "q": q, "subsets": subsets}


def correctness_specs() -> list[tuple[str, int, int]]:
    selections = {
        "gf4": (3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15),
        "gf8": (3, 5, 6, 7, 9, 15, 17, 31, 33, 63, 65, 127, 129, 255),
        "gf16": (3, 5, 7, 9, 15, 17, 31, 33, 63, 65, 127, 129),
    }
    result: list[tuple[str, int, int]] = []
    for field_name, counts in selections.items():
        field = c3.field_named(field_name)
        for q in counts:
            _, _, parent = split_geometry(q)
            shifts = (0,) if parent == field.order else (0, field.order - parent)
            result.extend((field_name, q, shift) for shift in shifts)
    return result


def root_counter(field_name: str, half: int, shift: int,
                 right_known_zero: bool = False) -> Counter:
    field = c3.field_named(field_name)
    multiplier = c3.fft_skew(field_name)[shift + half - 1]
    if right_known_zero:
        return Counter(
            xors=half if multiplier not in (0, 1) else 0,
            multiplications=half if multiplier not in (0, 1) else 0,
            loads=half, stores=2 * half,
        )
    return Counter(
        xors=half if multiplier in (0, 1) else 2 * half,
        multiplications=0 if multiplier in (0, 1) else half,
        loads=2 * half, stores=2 * half,
    )


@functools.lru_cache(maxsize=None)
def transform_counter(field_name: str, size: int, shift: int) -> Counter:
    if size == 1:
        return Counter()
    half = size // 2
    result = Counter()
    result.add(transform_counter(field_name, half, shift))
    result.add(transform_counter(field_name, half, shift + half))
    result.add(root_counter(field_name, half, shift))
    return result


def smaller_tail_counter(field_name: str, block: int, tail: int,
                         shift: int) -> Counter:
    small = ceil_power_of_two(tail)
    result = Counter()
    result.add(transform_counter(field_name, small, shift))
    size = small
    while size < block:
        result.add(root_counter(field_name, size, shift, True))
        size *= 2
    return result


def recursive_tail_counter(field_name: str, block: int, tail: int,
                           shift: int) -> Counter:
    if tail == block:
        copy = Counter()
        copy.add(transform_counter(field_name, block, shift))
        return copy
    half = block // 2
    result = Counter(stores=1)  # one compact schedule node at setup/execution boundary
    if tail <= half:
        result.add(recursive_tail_counter(field_name, half, tail, shift))
        result.add(root_counter(field_name, half, shift, True))
    else:
        result.add(transform_counter(field_name, half, shift))
        result.add(recursive_tail_counter(
            field_name, half, tail - half, shift + half
        ))
        result.add(root_counter(field_name, half, shift))
    return result


def dense_counter(rows: int, columns: int) -> Counter:
    return Counter(
        xors=rows * max(0, columns - 1),
        multiplications=rows * columns,
        loads=2 * rows * columns,
        stores=rows,
    )


def gaussian_setup_counter(size: int) -> Counter:
    cube = size ** 3
    return Counter(
        xors=cube, multiplications=cube, inversions=size,
        loads=2 * cube, stores=cube,
    )


def add_counters(*values: Counter) -> Counter:
    result = Counter()
    for value in values:
        result.add(value)
    return result


def method_row(name: str, profile: str, setup: Counter, execution: Counter,
               table_elements: int) -> dict[str, object]:
    return {
        "execution": dataclasses.asdict(execution),
        "method": name,
        "profile_class": profile,
        "setup": dataclasses.asdict(setup),
        "table_elements": table_elements,
    }


def operation_job(specification: tuple[str, int]) -> dict[str, object]:
    field_name, q = specification
    block, tail, parent = split_geometry(q)
    shift = 0
    head = transform_counter(field_name, block, shift)
    root = root_counter(field_name, block, shift)
    baseline = transform_counter(field_name, parent, shift)

    direct_tail_execution = dense_counter(block, tail)
    direct_setup = transform_counter(field_name, block, block)
    scaled_setup = Counter(stores=block * tail)
    scaled_setup.add(direct_setup, tail)
    direct_execution = add_counters(head, direct_tail_execution, root)

    smaller_execution = add_counters(
        head, smaller_tail_counter(field_name, block, tail, block), root
    )
    recursive_tail = recursive_tail_counter(field_name, block, tail, block)
    recursive_execution = add_counters(head, recursive_tail, root)
    recursive_setup = Counter(stores=max(1, tail.bit_count() + block.bit_length()))

    cross = dense_counter(tail, block)
    exact_tail_dense = dense_counter(tail, tail)
    exact_direct_execution = add_counters(
        head, cross, Counter(xors=tail, loads=2 * tail, stores=tail),
        exact_tail_dense,
    )
    exact_direct_setup = add_counters(
        gaussian_setup_counter(tail), Counter(stores=tail * (block + tail))
    )

    divided = tail * max(0, tail - 1) // 2
    triangular = tail * (tail + 1) // 2
    newton_execution = add_counters(
        head, cross, Counter(xors=tail, loads=2 * tail, stores=tail),
        Counter(
            xors=divided + 2 * triangular,
            multiplications=divided + 2 * triangular + tail,
            loads=3 * divided + 4 * triangular,
            stores=divided + 2 * triangular,
        ),
    )
    newton_setup = Counter(
        xors=tail * max(0, tail - 1),
        multiplications=tail * max(0, tail - 1),
        inversions=divided + tail,
        loads=2 * tail * max(0, tail - 1),
        stores=tail * max(0, tail - 1),
    )

    methods = (
        method_row("parent_full_padded", "wire-compatible-parent",
                   Counter(), baseline, len(c3.fft_skew(field_name))),
        method_row("parent_direct_dense_tail", "wire-compatible-parent",
                   scaled_setup, direct_execution, block * tail),
        method_row("parent_smaller_padded_tail", "wire-compatible-parent",
                   Counter(), smaller_execution, 0),
        method_row("parent_recursive_tail", "wire-compatible-parent",
                   recursive_setup, recursive_execution, 0),
        method_row("exact_direct_dense_tail", "new-exact-profile",
                   exact_direct_setup, exact_direct_execution,
                   tail * (block + tail)),
        method_row("exact_newton_tail", "new-exact-profile",
                   newton_setup, newton_execution,
                   divided + triangular + tail + block * tail),
    )
    return {
        "block": block, "field": field_name, "methods": methods,
        "parent": parent, "q": q, "shift": shift, "tail": tail,
    }


def operation_specs() -> list[tuple[str, int]]:
    result: list[tuple[str, int]] = []
    for block in (2, 4, 8, 16, 32, 64, 128):
        result.extend(("gf8", block + tail) for tail in range(1, block))
    for block in (2, 4, 8, 16, 32, 64):
        result.extend(("gf16", block + tail) for tail in range(1, block))
    for block in (128, 256, 512, 1024, 2048):
        tails = sorted({1, 2, 3, block // 8, block // 4,
                        block // 2, block - 1})
        result.extend(("gf16", block + tail) for tail in tails)
    return result


def counter_cost(counter: dict[str, int]) -> int:
    return (
        3 * counter["multiplications"] + counter["xors"] +
        counter["loads"] + 2 * counter["stores"] +
        12 * counter["inversions"]
    )


METHOD_IDS = (
    "parent_full_padded",
    "parent_direct_dense_tail",
    "parent_smaller_padded_tail",
    "parent_recursive_tail",
    "exact_direct_dense_tail",
    "exact_newton_tail",
)


def modeled_choices(rows: Sequence[dict[str, object]]) -> tuple[
    list[dict[str, object]], dict[str, object], dict[str, object]
]:
    """Return a compact deterministic offline table and both best cells.

    Each cell is ``[shard_bytes, batch, reuse, parent_method_id,
    parent_gain, exact_method_id, exact_gain]``.  Method IDs index
    ``METHOD_IDS``.  This avoids repeating two large dictionaries for every
    otherwise identical geometry/cell.
    """
    tables: list[dict[str, object]] = []
    best_parent: dict[str, object] | None = None
    best_exact: dict[str, object] | None = None
    for operation in rows:
        by_name = {str(method["method"]): method
                   for method in operation["methods"]}
        baseline = by_name["parent_full_padded"]
        cells: list[list[object]] = []
        for shard_bytes in (64, 1024, 65536, 1048576):
            symbol_bytes = 1 if operation["field"] == "gf8" else 2
            symbols = max(1, shard_bytes // symbol_bytes)
            for batch in (1, 8, 64):
                for reuse in (1, 8, 64, 1024):
                    def score(method):
                        return (
                            counter_cost(method["execution"]) * symbols * batch +
                            counter_cost(method["setup"]) / reuse
                        )

                    baseline_score = score(baseline)
                    parent_candidates = [
                        method for method in operation["methods"]
                        if method["profile_class"] == "wire-compatible-parent"
                    ]
                    exact_candidates = [
                        method for method in operation["methods"]
                        if method["profile_class"] == "new-exact-profile"
                    ]
                    parent_winner = min(parent_candidates, key=score)
                    exact_winner = min(exact_candidates, key=score)
                    parent_score = score(parent_winner)
                    exact_score = score(exact_winner)
                    parent_gain = round(baseline_score / parent_score, 9)
                    exact_gain = round(baseline_score / exact_score, 9)
                    cells.append([
                        shard_bytes, batch, reuse,
                        METHOD_IDS.index(str(parent_winner["method"])),
                        parent_gain,
                        METHOD_IDS.index(str(exact_winner["method"])),
                        exact_gain,
                    ])
                    parent_cell = {
                        "batch": batch, "field": operation["field"],
                        "gain_over_full_padded": parent_gain,
                        "parent": operation["parent"],
                        "profile_class": "wire-compatible-parent",
                        "q": operation["q"], "reuse": reuse,
                        "score": round(parent_score, 3),
                        "shard_bytes": shard_bytes,
                        "tail": operation["tail"],
                        "winner": parent_winner["method"],
                    }
                    exact_cell = {
                        "batch": batch, "field": operation["field"],
                        "gain_over_full_padded": exact_gain,
                        "parent": operation["parent"],
                        "profile_class": "new-exact-profile",
                        "q": operation["q"], "reuse": reuse,
                        "score": round(exact_score, 3),
                        "shard_bytes": shard_bytes,
                        "tail": operation["tail"],
                        "winner": exact_winner["method"],
                    }
                    if (best_parent is None or parent_gain >
                            float(best_parent["gain_over_full_padded"])):
                        best_parent = parent_cell
                    if (best_exact is None or exact_gain >
                            float(best_exact["gain_over_full_padded"])):
                        best_exact = exact_cell
        tables.append({
            "block": operation["block"], "cells": cells,
            "field": operation["field"], "parent": operation["parent"],
            "q": operation["q"], "tail": operation["tail"],
        })
    assert best_parent is not None and best_exact is not None
    return tables, best_parent, best_exact


def dependency_hashes() -> dict[str, str]:
    return {
        "c3_study.py": hashlib.sha256(C3_PATH.read_bytes()).hexdigest(),
    }


def source_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def deterministic_analysis(workers: int) -> dict[str, object]:
    if workers <= 0:
        raise ValueError("workers must be positive")
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        correctness = list(executor.map(correctness_job, correctness_specs(),
                                        chunksize=1))
        operations = list(executor.map(operation_job, operation_specs(),
                                       chunksize=1))
        mds = list(executor.map(
            mds_job, (3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15), chunksize=1
        ))
    correctness.sort(key=lambda row: (
        str(row["field"]), int(row["q"]), int(row["shift"])
    ))
    operations.sort(key=lambda row: (str(row["field"]), int(row["q"])))
    mds.sort(key=lambda row: int(row["q"]))
    choices, best_compatible, best_exact = modeled_choices(operations)
    return {
        "best_modeled_exact": best_exact,
        "best_modeled_wire_compatible": best_compatible,
        "choice_method_ids": METHOD_IDS,
        "correctness": correctness,
        "dependencies": dependency_hashes(),
        "disposition": {
            "exact_profile": "inconclusive-new-profile",
            "parent_direct_dense_tail": "inconclusive-tiny-tail",
            "parent_recursive_tail": "inconclusive-needs-fused-simd",
            "parent_smaller_padded_tail": "model-candidate-needs-cpp",
            "production_promotion": "none",
        },
        "mds": mds,
        "modeled_choices": choices,
        "operation_rows": operations,
        "schema": SCHEMA,
        "source_sha256": source_sha256(),
        "summary": {
            "correctness_comparisons": sum(int(row["comparisons"])
                                            for row in correctness),
            "correctness_jobs": len(correctness),
            "exact_difference_vectors": sum(
                int(row["parent_exact_difference_vectors"])
                for row in correctness
            ),
            "mds_coordinate_subsets": sum(int(row["subsets"]) for row in mds),
            "modeled_cells": sum(len(row["cells"]) for row in choices) * 2,
            "operation_geometries": len(operations),
            "vectors": sum(int(row["vectors"]) for row in correctness),
        },
    }


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True) + "\n").encode("ascii")


def seal(value: dict[str, object]) -> dict[str, object]:
    result = dict(value)
    result.pop("integrity_sha256", None)
    result["integrity_sha256"] = hashlib.sha256(canonical_bytes(result)).hexdigest()
    return result


def check_seal(value: dict[str, object]) -> bool:
    expected = value.get("integrity_sha256")
    candidate = dict(value)
    candidate.pop("integrity_sha256", None)
    return isinstance(expected, str) and expected == hashlib.sha256(
        canonical_bytes(candidate)
    ).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(canonical_bytes(value))
    temporary.replace(path)


def run_analysis(arguments: argparse.Namespace) -> int:
    allowed = sorted(os.sched_getaffinity(0))
    if arguments.workers > len(allowed):
        raise ValueError("workers exceed the allowed CPU set")
    analysis = deterministic_analysis(arguments.workers)
    artifact = seal({
        "analysis": analysis,
        "generation": {
            "allowed_cpus": allowed,
            "hostname": platform.node(),
            "platform": platform.platform(),
            "python": platform.python_version(),
            "workers": arguments.workers,
        },
        "schema": SCHEMA,
    })
    write_json(Path(arguments.output), artifact)
    print(json.dumps(analysis["summary"], sort_keys=True))
    return 0


def verify_result(arguments: argparse.Namespace) -> int:
    value = json.loads(Path(arguments.artifact).read_text(encoding="ascii"))
    if value.get("schema") != SCHEMA or not check_seal(value):
        raise ValueError("artifact schema or integrity digest is invalid")
    analysis = value.get("analysis")
    if not isinstance(analysis, dict):
        raise ValueError("artifact analysis is missing")
    if analysis.get("source_sha256") != source_sha256():
        raise ValueError("artifact source hash is stale")
    if analysis.get("dependencies") != dependency_hashes():
        raise ValueError("artifact dependency hash is stale")
    recomputed = deterministic_analysis(arguments.workers)
    if canonical_bytes(recomputed) != canonical_bytes(analysis):
        raise ValueError("artifact differs from deterministic recomputation")
    print(json.dumps(analysis["summary"], sort_keys=True))
    return 0


def affinity() -> list[int]:
    return sorted(os.sched_getaffinity(0))


def median_mad(values: Sequence[int]) -> tuple[int, int]:
    median = int(statistics.median(values))
    return median, int(statistics.median(abs(value - median) for value in values))


def time_call(function, iterations: int) -> int:
    started = time.perf_counter_ns()
    for _ in range(iterations):
        function()
    return (time.perf_counter_ns() - started) // iterations


def benchmark_case(field_name: str, q: int, samples: int,
                   iterations: int) -> dict[str, object]:
    field = c3.field_named(field_name)
    block, tail, parent = split_geometry(q)
    random_source = random.Random(stable_seed("timing", field_name, q, 0))
    values = tuple(random_source.randrange(field.order) for _ in range(q))
    direct_plan = build_direct_tail_plan(field, block, tail, block)
    exact_plan = build_exact_split_plan(field, q, 0)
    newton_plan = build_exact_newton_tail_plan(field, q)
    calls = {
        "parent_full_padded": lambda: c3.padded_lch_inverse(
            field, list(values) + [0] * (parent - q), 0
        ),
        "parent_direct_dense_tail": lambda: parent_full_block(
            field, values, 0, "direct", direct_plan
        ),
        "parent_smaller_padded_tail": lambda: parent_full_block(
            field, values, 0, "smaller_padded", direct_plan
        ),
        "parent_recursive_tail": lambda: parent_full_block(
            field, values, 0, "recursive", direct_plan
        ),
        "exact_direct_dense_tail": lambda: exact_split_direct(
            field, exact_plan, values
        ),
        "exact_newton_tail": lambda: exact_split_newton(
            field, newton_plan, values
        ),
    }
    parent_oracle = calls["parent_full_padded"]()
    if any(calls[name]() != parent_oracle for name in (
        "parent_direct_dense_tail", "parent_smaller_padded_tail",
        "parent_recursive_tail",
    )):
        raise AssertionError("benchmark parent candidates disagree")
    if calls["exact_direct_dense_tail"]() != calls["exact_newton_tail"]():
        raise AssertionError("benchmark exact candidates disagree")
    for call in calls.values():
        call()
    setup_calls = {
        "parent_direct_dense_tail": lambda: build_direct_tail_plan(
            field, block, tail, block
        ),
        "exact_direct_dense_tail": lambda: build_exact_split_plan(field, q, 0),
        "exact_newton_tail": lambda: build_exact_newton_tail_plan(field, q),
    }
    for call in setup_calls.values():
        call()
    timing: dict[str, dict[str, int]] = {}
    for name, call in calls.items():
        measurements = [time_call(call, iterations) for _ in range(samples)]
        median, mad = median_mad(measurements)
        timing[name] = {"mad_ns": mad, "median_ns": median}
    setup_timing: dict[str, dict[str, int]] = {}
    setup_iterations = max(1, iterations // 8)
    for name, call in setup_calls.items():
        measurements = [time_call(call, setup_iterations)
                        for _ in range(samples)]
        median, mad = median_mad(measurements)
        setup_timing[name] = {"mad_ns": mad, "median_ns": median}
    return {
        "block": block, "field": field_name, "iterations": iterations,
        "methods": timing, "parent": parent, "q": q,
        "samples": samples, "setup_iterations": setup_iterations,
        "setup_methods": setup_timing, "tail": tail,
    }


def run_benchmark(arguments: argparse.Namespace) -> int:
    if len(affinity()) != 1:
        raise ValueError("benchmark must be pinned to one allowed CPU")
    if arguments.samples < 3 or arguments.iterations <= 0:
        raise ValueError("benchmark requires at least three samples")
    cases = [benchmark_case(arguments.field, q, arguments.samples,
                            arguments.iterations) for q in arguments.counts]
    artifact = seal({
        "affinity": affinity(),
        "cases": cases,
        "dependencies": dependency_hashes(),
        "environment": {
            "hostname": platform.node(), "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "schema": TIMING_SCHEMA,
        "source_sha256": source_sha256(),
    })
    write_json(Path(arguments.output), artifact)
    print(json.dumps({"affinity": affinity(), "cases": len(cases)},
                     sort_keys=True))
    return 0


def verify_timing(arguments: argparse.Namespace) -> int:
    value = json.loads(Path(arguments.artifact).read_text(encoding="ascii"))
    if value.get("schema") != TIMING_SCHEMA or not check_seal(value):
        raise ValueError("timing schema or integrity digest is invalid")
    if value.get("source_sha256") != source_sha256():
        raise ValueError("timing source hash is stale")
    if value.get("dependencies") != dependency_hashes():
        raise ValueError("timing dependency hash is stale")
    recorded_affinity = value.get("affinity")
    if not isinstance(recorded_affinity, list) or len(recorded_affinity) != 1:
        raise ValueError("timing was not recorded under one-CPU affinity")
    cases = value.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("timing contains no cases")
    required = {
        "parent_full_padded", "parent_direct_dense_tail",
        "parent_smaller_padded_tail", "parent_recursive_tail",
        "exact_direct_dense_tail", "exact_newton_tail",
    }
    required_setup = {
        "parent_direct_dense_tail", "exact_direct_dense_tail",
        "exact_newton_tail",
    }
    for case in cases:
        if set(case.get("methods", {})) != required:
            raise ValueError("timing method set is incomplete")
        if set(case.get("setup_methods", {})) != required_setup:
            raise ValueError("timing setup-method set is incomplete")
        for measurement in case["methods"].values():
            if measurement.get("median_ns", 0) <= 0 or measurement.get("mad_ns", -1) < 0:
                raise ValueError("timing measurement is invalid")
        for measurement in case["setup_methods"].values():
            if measurement.get("median_ns", 0) <= 0 or measurement.get("mad_ns", -1) < 0:
                raise ValueError("timing setup measurement is invalid")
    print(json.dumps({"affinity": recorded_affinity, "cases": len(cases)},
                     sort_keys=True))
    return 0


def self_test() -> int:
    field = c3.field_named("gf4")
    for q in (3, 5, 7):
        row = correctness_job(("gf4", q, 0))
        if row["vectors"] <= 0:
            raise AssertionError("self-test ran no vectors")
    sample = seal({"schema": "test", "value": [1, 2, 3]})
    if not check_seal(sample):
        raise AssertionError("valid integrity seal failed")
    sample["value"][1] = 9
    if check_seal(sample):
        raise AssertionError("tampered integrity seal was accepted")
    try:
        split_geometry(4)
    except ValueError:
        pass
    else:
        raise AssertionError("power-of-two q was accepted")
    values = (0, 0, 1)
    exact = exact_split_direct(field, build_exact_split_plan(field, 3, 0), values)
    padded = c3.padded_lch_inverse(field, list(values) + [0], 0)
    if exact == padded[:3] and not padded[3]:
        raise AssertionError("self-test lacks a parent/exact distinction witness")
    print(json.dumps({"q": [3, 5, 7], "status": "ok"}, sort_keys=True))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser("self-test")
    analyze = commands.add_parser("analyze")
    analyze.add_argument("--workers", type=int,
                         default=min(len(os.sched_getaffinity(0)), 128))
    analyze.add_argument("--output", required=True)
    verify = commands.add_parser("verify")
    verify.add_argument("artifact")
    verify.add_argument("--workers", type=int, default=1)
    benchmark = commands.add_parser("benchmark")
    benchmark.add_argument("--field", choices=("gf8", "gf16"), default="gf8")
    benchmark.add_argument("--counts", nargs="+", type=int,
                           default=(3, 5, 9, 17))
    benchmark.add_argument("--samples", type=int, default=7)
    benchmark.add_argument("--iterations", type=int, default=32)
    benchmark.add_argument("--output", required=True)
    timing = commands.add_parser("verify-timing")
    timing.add_argument("artifact")
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
        if arguments.command == "verify-timing":
            return verify_timing(arguments)
        raise AssertionError("unreachable command")
    except (AssertionError, OSError, RuntimeError, ValueError) as error:
        print(f"c4: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
