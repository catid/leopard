#!/usr/bin/env python3
"""Generate the static FF8 XOR circuits used by the experimental backend.

The generated C++ deliberately contains named wire variables and literal XOR
statements.  The payload kernels include it after defining XorValue() for the
selected scalar or SIMD value type.
"""

from __future__ import print_function

import argparse
import collections
import difflib
import hashlib
import random
import sys
import time
from pathlib import Path


FIELD_BITS = 8
FIELD_ORDER = 1 << FIELD_BITS
FIELD_MODULUS = FIELD_ORDER - 1
FIELD_POLYNOMIAL = 0x11D
CANTOR_BASIS = (1, 214, 152, 146, 86, 200, 88, 230)
WIRE_COUNT_MULTIPLY = 8
WIRE_COUNT_BUTTERFLY = 16
GENERATOR_VERSION = b"LeopardFF8XorCircuits-v3-greedy\0"


def add_mod(a, b):
    """Match LeopardFF8.cpp AddMod(), including its representation of zero."""
    total = a + b
    return (total + (total >> FIELD_BITS)) & 0xFF


def make_field_tables():
    """Reproduce InitializeLogarithmTables() from LeopardFF8.cpp exactly."""
    log_lut = [0] * FIELD_ORDER
    exp_lut = [0] * FIELD_ORDER

    # At this point ExpLUT is temporarily indexed by polynomial-basis values
    # and contains their logarithms.
    state = 1
    for exponent in range(FIELD_MODULUS):
        exp_lut[state] = exponent
        state <<= 1
        if state >= FIELD_ORDER:
            state ^= FIELD_POLYNOMIAL
    exp_lut[0] = FIELD_MODULUS

    # Build the Cantor-coordinate-byte -> polynomial-basis-byte map, then turn
    # it into the final coordinate-byte -> logarithm table.
    log_lut[0] = 0
    for bit, basis in enumerate(CANTOR_BASIS):
        width = 1 << bit
        for value in range(width):
            log_lut[value + width] = log_lut[value] ^ basis

    for value in range(FIELD_ORDER):
        log_lut[value] = exp_lut[log_lut[value]]

    # Regenerate ExpLUT as logarithm -> Cantor-coordinate byte.  Leopard uses
    # index 255 as a redundant representation of exponent zero, not as zero.
    for value in range(FIELD_ORDER):
        exp_lut[log_lut[value]] = value
    exp_lut[FIELD_MODULUS] = exp_lut[0]

    assert log_lut[0] == FIELD_MODULUS
    assert exp_lut[0] == 1
    assert exp_lut[FIELD_MODULUS] == 1
    assert sorted(log_lut[1:]) == list(range(FIELD_MODULUS))
    return tuple(log_lut), tuple(exp_lut)


LOG_LUT, EXP_LUT = make_field_tables()


def scalar_multiply_log(value, log_multiplier):
    """Match LeopardFF8.cpp MultiplyLog(value, log_multiplier)."""
    if value == 0:
        return 0
    return EXP_LUT[add_mod(LOG_LUT[value], log_multiplier)]


def cantor_to_polynomial(value):
    """Convert a coordinate byte without using either generated field table."""
    polynomial_value = 0
    for bit, basis in enumerate(CANTOR_BASIS):
        if (value >> bit) & 1:
            polynomial_value ^= basis
    return polynomial_value


CANTOR_TO_POLYNOMIAL = tuple(
    cantor_to_polynomial(value) for value in range(FIELD_ORDER))
POLYNOMIAL_TO_CANTOR_LIST = [0] * FIELD_ORDER
for _coordinate, _polynomial_value in enumerate(CANTOR_TO_POLYNOMIAL):
    POLYNOMIAL_TO_CANTOR_LIST[_polynomial_value] = _coordinate
POLYNOMIAL_TO_CANTOR = tuple(POLYNOMIAL_TO_CANTOR_LIST)
del POLYNOMIAL_TO_CANTOR_LIST, _coordinate, _polynomial_value


def polynomial_multiply(left, right):
    """Independent shift/XOR GF(256) reference in polynomial coordinates."""
    product = 0
    for _ in range(FIELD_BITS):
        if right & 1:
            product ^= left
        right >>= 1
        left <<= 1
        if left & FIELD_ORDER:
            left ^= FIELD_POLYNOMIAL
    return product


POLYNOMIAL_EXP = [1] * FIELD_ORDER
for _exponent in range(1, FIELD_ORDER):
    POLYNOMIAL_EXP[_exponent] = polynomial_multiply(
        POLYNOMIAL_EXP[_exponent - 1], 2)
POLYNOMIAL_EXP = tuple(POLYNOMIAL_EXP)
del _exponent


def independent_scalar_multiply_log(value, log_multiplier):
    polynomial_product = polynomial_multiply(
        CANTOR_TO_POLYNOMIAL[value], POLYNOMIAL_EXP[log_multiplier])
    return POLYNOMIAL_TO_CANTOR[polynomial_product]


def linear_map(width, function):
    """Return row masks for a linear function using the documented convention."""
    rows = [0] * width
    for column in range(width):
        output = function(1 << column)
        if output < 0 or output >= (1 << width):
            raise ValueError("linear function returned an out-of-range value")
        for row in range(width):
            rows[row] |= ((output >> row) & 1) << column
    return tuple(rows)


def apply_matrix(rows, value):
    output = 0
    for row, mask in enumerate(rows):
        # Keep regeneration compatible with the older Python 3 versions that
        # may accompany this project's CMake 3.7 minimum.
        output |= (bin(mask & value).count("1") & 1) << row
    return output


def multiplication_matrix(log_multiplier):
    return linear_map(
        WIRE_COUNT_MULTIPLY,
        lambda value: scalar_multiply_log(value, log_multiplier))


def scalar_forward_butterfly(state, skew, multiplication_matrices):
    x = state & 0xFF
    y = (state >> 8) & 0xFF
    if skew != FIELD_MODULUS:
        x ^= apply_matrix(multiplication_matrices[skew], y)
    y ^= x
    return x | (y << 8)


def scalar_inverse_butterfly(state, skew, multiplication_matrices):
    x = state & 0xFF
    y = (state >> 8) & 0xFF
    y ^= x
    if skew != FIELD_MODULUS:
        x ^= apply_matrix(multiplication_matrices[skew], y)
    return x | (y << 8)


def cancel_adjacent_duplicate_gates(gates):
    """Remove adjacent equal CNOTs, including duplicates exposed by removal."""
    output = []
    for gate in gates:
        if output and output[-1] == gate:
            output.pop()
        else:
            output.append(gate)
    return tuple(output)


def synthesize_reversible_map(rows):
    """Synthesize an invertible binary row map using only dst ^= src gates."""
    width = len(rows)
    reduced = list(rows)
    reduction_gates = []

    for column in range(width):
        if ((reduced[column] >> column) & 1) == 0:
            pivot = None
            for candidate in range(column + 1, width):
                if (reduced[candidate] >> column) & 1:
                    pivot = candidate
                    break
            if pivot is None:
                raise ValueError("matrix is singular")

            # Do not introduce a physical swap: adding the lower row supplies
            # the missing pivot and is itself a reversible CNOT.
            reduced[column] ^= reduced[pivot]
            reduction_gates.append((column, pivot))

        for row in range(width):
            if row != column and ((reduced[row] >> column) & 1):
                reduced[row] ^= reduced[column]
                reduction_gates.append((row, column))

    if reduced != [1 << row for row in range(width)]:
        raise AssertionError("row reduction did not produce identity")

    # If E_n ... E_1 A = I, the implementation order for A is E_n ... E_1
    # as a list of sequential in-place wire operations, hence reverse the row
    # operations recorded during reduction.
    return cancel_adjacent_duplicate_gates(reversed(reduction_gates))


def circuit_matrix(width, gates):
    rows = [1 << row for row in range(width)]
    for destination, source in gates:
        if destination == source:
            raise AssertionError("a CNOT gate may not use the same wire twice")
        rows[destination] ^= rows[source]
    return tuple(rows)


def apply_circuit(state, gates):
    for destination, source in gates:
        state ^= ((state >> source) & 1) << destination
    return state


def circuit_depth(gates, width):
    """Count CNOT layers, conservatively forbidding shared wires per layer."""
    last_layer = [0] * width
    maximum = 0
    for destination, source in gates:
        layer = max(last_layer[destination], last_layer[source]) + 1
        last_layer[destination] = layer
        last_layer[source] = layer
        maximum = max(maximum, layer)
    return maximum


_POPCOUNT_16 = bytes(
    bin(value).count("1") for value in range(1 << WIRE_COUNT_BUTTERFLY))


def bit_count(value):
    # All emitted 8- and 16-wire maps use this compact lookup.  Keep the
    # fallback for the representative 32-wire synthesis proof and for older
    # Python versions that do not provide int.bit_count().
    if value < len(_POPCOUNT_16):
        return _POPCOUNT_16[value]
    return bin(value).count("1")


def transpose_matrix(rows):
    width = len(rows)
    return tuple(
        sum(((rows[row] >> column) & 1) << row for row in range(width))
        for column in range(width))


def inverse_matrix(rows):
    """Invert a binary matrix independently of the circuit synthesizer."""
    width = len(rows)
    left = list(rows)
    right = [1 << row for row in range(width)]
    for column in range(width):
        pivot = next((row for row in range(column, width)
                      if (left[row] >> column) & 1), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        if pivot != column:
            left[column], left[pivot] = left[pivot], left[column]
            right[column], right[pivot] = right[pivot], right[column]
        for row in range(width):
            if row != column and ((left[row] >> column) & 1):
                left[row] ^= left[column]
                right[row] ^= right[column]
    if left != [1 << row for row in range(width)]:
        raise AssertionError("matrix inversion did not produce identity")
    return tuple(right)


def transpose_circuit(gates):
    """Turn a circuit for A into one for transpose(A)."""
    return tuple((source, destination)
                 for destination, source in reversed(gates))


def cnot_gates_commute(left, right):
    left_destination, left_source = left
    right_destination, right_source = right
    return (left_destination != right_source and
            right_destination != left_source)


def cancel_commuting_duplicate_gates(gates):
    """Cancel equal CNOTs separated only by gates that commute with them."""
    output = list(gates)
    while True:
        changed = False
        for left in range(len(output)):
            for right in range(left + 1, len(output)):
                if output[right] == output[left]:
                    del output[right]
                    del output[left]
                    changed = True
                    break
                if not cnot_gates_commute(output[left], output[right]):
                    break
            if changed:
                break
        if not changed:
            return tuple(output)


def reverse_bits(value, width):
    result = 0
    for unused in range(width):
        result = (result << 1) | (value & 1)
        value >>= 1
    return result


def synthesis_column_orders(rows):
    """Return deterministic column orders for the synthesis portfolio."""
    width = len(rows)
    natural = list(range(width))
    result = []

    def add(name, order):
        order = tuple(order)
        if all(existing != order for unused_name, existing in result):
            result.append((name, order))

    add("forward", natural)
    add("reverse", reversed(natural))
    add("sparse-column", sorted(
        natural,
        key=lambda column: (
            sum((rows[row] >> column) & 1 for row in natural), column)))
    add("dense-column", sorted(
        natural,
        key=lambda column: (
            -sum((rows[row] >> column) & 1 for row in natural), column)))
    add("sparse-row", sorted(
        natural, key=lambda row: (bit_count(rows[row]), row)))
    add("dense-row", sorted(
        natural, key=lambda row: (-bit_count(rows[row]), row)))
    bit_width = (width - 1).bit_length()
    if 1 << bit_width == width:
        add("bit-reverse", sorted(
            natural, key=lambda value: reverse_bits(value, bit_width)))

    # All rotations are inexpensive and materially useful for the 8x8 maps.
    # The compact order set above is enough for the already-cheap 16x16 direct
    # butterflies and keeps checked-in regeneration comfortably finite.
    if width <= WIRE_COUNT_MULTIPLY:
        for shift in range(1, width):
            add("rotate-%d" % shift,
                natural[shift:] + natural[:shift])
    return tuple(result)


PIVOT_MODES = ("first", "last", "sparse", "dense")


def synthesize_reversible_map_ordered(
        rows, column_order, pivot_mode, reverse_elimination):
    """Gauss-Jordan CNOT synthesis with explicit deterministic choices."""
    width = len(rows)
    reduced = list(rows)
    reduction_gates = []
    processed = set()

    for column in column_order:
        if ((reduced[column] >> column) & 1) == 0:
            candidates = [
                row for row in range(width)
                if (row != column and row not in processed and
                    ((reduced[row] >> column) & 1))]
            if not candidates:
                raise ValueError("matrix is singular")
            if pivot_mode == "first":
                pivot = min(candidates)
            elif pivot_mode == "last":
                pivot = max(candidates)
            elif pivot_mode == "sparse":
                pivot = min(candidates, key=lambda row: (
                    bit_count(reduced[column] ^ reduced[row]),
                    bit_count(reduced[row]), row))
            elif pivot_mode == "dense":
                pivot = min(candidates, key=lambda row: (
                    -bit_count(reduced[column] ^ reduced[row]),
                    -bit_count(reduced[row]), row))
            else:
                raise ValueError("unknown pivot mode: %s" % pivot_mode)
            reduced[column] ^= reduced[pivot]
            reduction_gates.append((column, pivot))

        elimination_rows = [
            row for row in range(width)
            if row != column and ((reduced[row] >> column) & 1)]
        if reverse_elimination:
            elimination_rows.reverse()
        for row in elimination_rows:
            reduced[row] ^= reduced[column]
            reduction_gates.append((row, column))
        processed.add(column)

    if reduced != [1 << row for row in range(width)]:
        raise AssertionError("ordered row reduction did not produce identity")
    return cancel_adjacent_duplicate_gates(reversed(reduction_gates))


_EXACT_CNOT_TABLES = {}


def exact_cnot_table(width):
    """Return lexical shortest CNOT circuits for every map up to width four."""
    if width > 4:
        raise ValueError("bounded exact CNOT search supports at most four wires")
    cached = _EXACT_CNOT_TABLES.get(width)
    if cached is not None:
        return cached
    identity = tuple(1 << wire for wire in range(width))
    paths = {identity: ()}
    pending = collections.deque((identity,))
    gates = tuple(
        (destination, source)
        for destination in range(width)
        for source in range(width)
        if destination != source)
    while pending:
        rows = pending.popleft()
        path = paths[rows]
        for gate in gates:
            candidate = list(rows)
            candidate[gate[0]] ^= candidate[gate[1]]
            candidate = tuple(candidate)
            if candidate not in paths:
                paths[candidate] = path + (gate,)
                pending.append(candidate)
    _EXACT_CNOT_TABLES[width] = paths
    return paths


def optimize_exact_windows(gates, width, maximum_window=8):
    """Replace small contiguous CNOT windows with exact shortest circuits."""
    gates = tuple(gates)
    while True:
        best = None
        for begin in range(len(gates)):
            touched = set()
            for end in range(begin, min(len(gates), begin + maximum_window)):
                touched.update(gates[end])
                if len(touched) > 4:
                    break
                wires = tuple(sorted(touched))
                positions = {wire: index for index, wire in enumerate(wires)}
                local_gates = tuple(
                    (positions[destination], positions[source])
                    for destination, source in gates[begin:end + 1])
                local_rows = circuit_matrix(len(wires), local_gates)
                replacement = exact_cnot_table(len(wires))[local_rows]
                replacement = tuple(
                    (wires[destination], wires[source])
                    for destination, source in replacement)
                if len(replacement) >= end + 1 - begin:
                    continue
                candidate = (
                    gates[:begin] + replacement + gates[end + 1:])
                key = circuit_key(candidate, width)
                if best is None or key < best[0]:
                    best = (key, candidate)
        if best is None:
            return gates
        gates = best[1]


GREEDY_OUTPUT_SIDE = 0
GREEDY_INPUT_SIDE = 1
GREEDY_PAIR_BEAM_FACTOR = 2
_GREEDY_IDENTITY_DISTANCE_TABLES = {}
_GREEDY_CANONICAL_CIRCUITS = {}


def greedy_identity_distance_table(width):
    cached = _GREEDY_IDENTITY_DISTANCE_TABLES.get(width)
    if cached is not None:
        return cached
    if width > WIRE_COUNT_BUTTERFLY:
        return None
    table = tuple(bytes(
        _POPCOUNT_16[value ^ (1 << row)]
        for value in range(1 << width)) for row in range(width))
    _GREEDY_IDENTITY_DISTANCE_TABLES[width] = table
    return table


def greedy_synthesis_cost(rows, inverse_rows):
    """Hamming distance of a map and its inverse from identity."""
    distance_table = greedy_identity_distance_table(len(rows))
    if distance_table is not None:
        return sum(
            distance_table[index][row]
            for index, row in enumerate(rows)) + sum(
                distance_table[index][row]
                for index, row in enumerate(inverse_rows))
    return sum(
        bit_count(row ^ (1 << index))
        for index, row in enumerate(rows)) + sum(
            bit_count(row ^ (1 << index))
            for index, row in enumerate(inverse_rows))


def add_cnot_on_right(rows, destination, source):
    """Return rows * CNOT(destination, source).

    Left multiplication by this CNOT adds source row to destination row.
    Right multiplication instead adds destination column to source column.
    """
    return tuple(
        row ^ (((row >> destination) & 1) << source)
        for row in rows)


def apply_bidirectional_cnot(rows, inverse_rows, move):
    """Apply one output- or input-side CNOT to a map and its inverse."""
    side, destination, source = move
    if destination == source:
        raise ValueError("a CNOT gate may not use the same wire twice")

    if side == GREEDY_OUTPUT_SIDE:
        changed_rows = list(rows)
        changed_rows[destination] ^= changed_rows[source]
        return (tuple(changed_rows), add_cnot_on_right(
            inverse_rows, destination, source))
    if side == GREEDY_INPUT_SIDE:
        changed_inverse = list(inverse_rows)
        changed_inverse[destination] ^= changed_inverse[source]
        return (add_cnot_on_right(rows, destination, source),
                tuple(changed_inverse))
    raise ValueError("unknown greedy CNOT side: %r" % (side,))


def synthesize_reversible_map_greedy(rows):
    """Synthesize a map with deterministic bidirectional greedy reduction.

    Boyar-Peralta and Paar-style linear synthesis obtains short programs by
    greedily reusing profitable XOR relationships.  Their ordinary straight-
    line programs may allocate temporary signals, which this experiment's
    in-place kernels cannot do.  This candidate uses the reversible analogue:
    legal CNOTs are tried on both the input and output sides, and the selected
    move minimizes the Hamming distance of both the remainder matrix and its
    inverse from identity.  A bounded two-CNOT lookahead can cross a one-gate
    local minimum.

    For the eight-wire multipliers every step considers a pair; this produces
    materially smaller circuits at modest generation cost.  On wider maps a
    strictly improving single gate is preferred, and pair lookahead is used at
    local minima.  The pair search keeps the best 2*width first moves, ordered
    deterministically by cost and literal gate tuple.  If the greedy search
    still reaches a local minimum, ordinary row synthesis completes the map,
    so this function is total for every invertible input.
    """
    requested_rows = tuple(rows)
    width = len(requested_rows)
    identity = tuple(1 << wire for wire in range(width))
    requested_inverse = inverse_matrix(requested_rows)

    # A circuit and its reverse implement inverse maps.  Canonicalizing the
    # pair preserves determinism while halving work for the inverse multiplier
    # pairs and for matching FFT/IFFT butterflies.
    reverse_result = requested_inverse < requested_rows
    canonical_rows = requested_inverse if reverse_result else requested_rows
    remainder = canonical_rows
    inverse_remainder = requested_rows if reverse_result else requested_inverse
    cached = _GREEDY_CANONICAL_CIRCUITS.get(canonical_rows)
    if cached is not None:
        return tuple(reversed(cached)) if reverse_result else cached
    input_gates = []
    output_gates = []
    moves = tuple(
        (side, destination, source)
        for side in (GREEDY_OUTPUT_SIDE, GREEDY_INPUT_SIDE)
        for destination in range(width)
        for source in range(width)
        if destination != source)

    def append_move(move):
        side, destination, source = move
        if side == GREEDY_OUTPUT_SIDE:
            output_gates.append((destination, source))
        else:
            input_gates.append((destination, source))

    while remainder != identity:
        current_cost = greedy_synthesis_cost(
            remainder, inverse_remainder)
        first_candidates = []
        for move in moves:
            changed_rows, changed_inverse = apply_bidirectional_cnot(
                remainder, inverse_remainder, move)
            first_candidates.append((
                greedy_synthesis_cost(changed_rows, changed_inverse),
                move, changed_rows, changed_inverse))
        first_candidates.sort(key=lambda candidate: (
            candidate[0], candidate[1]))
        best_single = first_candidates[0]

        # A final single CNOT must not be padded with an unnecessary pair.
        # Wider matrices also take profitable single steps to bound runtime.
        if (best_single[2] == identity or
                (width > WIRE_COUNT_MULTIPLY and
                 best_single[0] < current_cost)):
            unused_cost, move, remainder, inverse_remainder = best_single
            append_move(move)
            continue

        best_pair = None
        beam_width = GREEDY_PAIR_BEAM_FACTOR * width
        for unused_cost, first_move, first_rows, first_inverse in \
                first_candidates[:beam_width]:
            for second_move in moves:
                if second_move == first_move:
                    continue
                changed_rows, changed_inverse = apply_bidirectional_cnot(
                    first_rows, first_inverse, second_move)
                changed_cost = greedy_synthesis_cost(
                    changed_rows, changed_inverse)
                key = (changed_cost, first_move, second_move)
                if (changed_cost < current_cost and
                        (best_pair is None or key < best_pair[0])):
                    best_pair = (
                        key, first_move, second_move,
                        changed_rows, changed_inverse)

        if best_pair is not None:
            (unused_key, first_move, second_move,
             remainder, inverse_remainder) = best_pair
            append_move(first_move)
            append_move(second_move)
            continue

        if best_single[0] < current_cost:
            unused_cost, move, remainder, inverse_remainder = best_single
            append_move(move)
            continue
        break

    # The maintained invariant is remainder = L * rows * R.  input_gates
    # implement R^-1 in their recorded order; reversed output_gates implement
    # L^-1.  A normal row circuit supplies any locally-minimal remainder.
    completion = synthesize_reversible_map(remainder)
    gates = cancel_commuting_duplicate_gates(
        tuple(input_gates) + completion + tuple(reversed(output_gates)))
    _GREEDY_CANONICAL_CIRCUITS[canonical_rows] = gates
    if reverse_result:
        gates = tuple(reversed(gates))
    if circuit_matrix(width, gates) != requested_rows:
        raise AssertionError("greedy synthesis candidate map mismatch")
    return gates


def synthesize_reversible_map_portfolio(rows):
    """Choose a verified circuit from deterministic synthesis variants."""
    width = len(rows)
    inverse_rows = inverse_matrix(rows)
    orientations = (
        ("direct", rows, lambda gates: gates),
        ("transpose", transpose_matrix(rows), transpose_circuit),
        ("inverse", inverse_rows, lambda gates: tuple(reversed(gates))),
        ("inverse-transpose", transpose_matrix(inverse_rows),
         lambda gates: tuple(reversed(transpose_circuit(gates)))),
    )
    candidates = {}
    for orientation_name, oriented_rows, transform in orientations:
        for order_name, column_order in synthesis_column_orders(oriented_rows):
            for pivot_mode in PIVOT_MODES:
                for reverse_elimination in (False, True):
                    gates = synthesize_reversible_map_ordered(
                        oriented_rows, column_order, pivot_mode,
                        reverse_elimination)
                    gates = cancel_commuting_duplicate_gates(transform(gates))
                    if circuit_matrix(width, gates) != rows:
                        raise AssertionError("portfolio candidate map mismatch")
                    name = "%s/%s/%s/%s" % (
                        orientation_name, order_name, pivot_mode,
                        "reverse" if reverse_elimination else "forward")
                    previous = candidates.get(gates)
                    if previous is None or name < previous:
                        candidates[gates] = name

    # The wider greedy butterfly schedules reduce source gate count but create
    # enough live-range pressure for GCC to spill the 16 named vector wires in
    # the portable/SIMD128 and AVX2 loops.  Keep the verified greedy candidate
    # for eight-wire multipliers, where the strict assembly census remains
    # spill-free.  Wider synthesis is still exercised directly by tests and is
    # available for a future register-pressure-aware/AVX-512-only portfolio.
    if width <= WIRE_COUNT_MULTIPLY:
        greedy_gates = synthesize_reversible_map_greedy(rows)
        if circuit_matrix(width, greedy_gates) != rows:
            raise AssertionError("greedy portfolio candidate map mismatch")
        previous = candidates.get(greedy_gates)
        if previous is None or "bidirectional-greedy" < previous:
            candidates[greedy_gates] = "bidirectional-greedy"

    # Exact peephole search is cheap for a bounded shortlist and can discover
    # CNOT templates that row reduction and commuting cancellation cannot.
    ranked = sorted(
        candidates.items(),
        key=lambda item: circuit_key(item[0], width) + (item[1],))
    for gates, name in ranked[:16]:
        optimized = optimize_exact_windows(gates, width)
        if circuit_matrix(width, optimized) != rows:
            raise AssertionError("exact-window rewrite changed the map")
        optimized_name = name + "/exact4x8"
        previous = candidates.get(optimized)
        if previous is None or optimized_name < previous:
            candidates[optimized] = optimized_name

    gates, name = min(
        candidates.items(),
        key=lambda item: circuit_key(item[0], width) + (item[1],))
    return gates, name


def direct_butterfly_circuit(skew, inverse, multiplication_matrices):
    gates = []
    matrix = multiplication_matrices[skew]

    if inverse:
        for bit in range(FIELD_BITS):
            gates.append((8 + bit, bit))

    if skew != FIELD_MODULUS:
        for output_bit in range(FIELD_BITS):
            for input_bit in range(FIELD_BITS):
                if (matrix[output_bit] >> input_bit) & 1:
                    gates.append((output_bit, 8 + input_bit))

    if not inverse:
        for bit in range(FIELD_BITS):
            gates.append((8 + bit, bit))

    return cancel_adjacent_duplicate_gates(gates)


def circuit_key(gates, width):
    return (len(gates), circuit_depth(gates, width), gates)


def choose_butterfly_circuit(rows, direct_gates):
    synthesized_gates, synthesized_name = \
        synthesize_reversible_map_portfolio(rows)
    direct_gates = cancel_commuting_duplicate_gates(direct_gates)
    direct_optimized = optimize_exact_windows(direct_gates, len(rows))
    if circuit_matrix(len(rows), synthesized_gates) != rows:
        raise AssertionError("synthesized butterfly circuit is incorrect")
    if circuit_matrix(len(rows), direct_gates) != rows:
        raise AssertionError("direct butterfly circuit is incorrect")
    if circuit_matrix(len(rows), direct_optimized) != rows:
        raise AssertionError("optimized direct butterfly circuit is incorrect")

    # The tuple comparison exactly implements gate count, then dependency
    # depth, then lexical gate ordering.  Identical gate lists are equivalent.
    candidates = (
        (synthesized_gates, "portfolio:" + synthesized_name),
        (direct_gates, "direct"),
        (direct_optimized, "direct/exact4x8"),
    )
    return min(candidates, key=lambda item:
               circuit_key(item[0], len(rows)) + (item[1],))


def validate_representative_wide_synthesis():
    """Prove the portfolio machinery on deterministic 32x32 maps."""
    width = 32
    deterministic_random = random.Random(0x32C107FF8)
    representative_maps = []
    for gate_count in (64, 127, 191):
        gates = []
        for unused in range(gate_count):
            destination = deterministic_random.randrange(width)
            source = deterministic_random.randrange(width - 1)
            if source >= destination:
                source += 1
            gates.append((destination, source))
        representative_maps.append(circuit_matrix(width, gates))

    for index, rows in enumerate(representative_maps):
        baseline = synthesize_reversible_map(rows)
        if circuit_matrix(width, baseline) != rows:
            raise AssertionError("32x32 baseline synthesis mismatch")

        # Exercise alternate order/pivot choices on every wide fixture.
        orders = synthesis_column_orders(rows)
        for order_name, order in (orders[0], orders[-1]):
            for pivot_mode in (PIVOT_MODES[0], PIVOT_MODES[-1]):
                gates = synthesize_reversible_map_ordered(
                    rows, order, pivot_mode, index % 2 != 0)
                if circuit_matrix(width, gates) != rows:
                    raise AssertionError(
                        "32x32 ordered synthesis mismatch: %s/%s" %
                        (order_name, pivot_mode))

        inverse_rows = inverse_matrix(rows)
        if transpose_matrix(transpose_matrix(rows)) != rows:
            raise AssertionError("32x32 transpose involution failed")
        validation_values = (
            0, (1 << width) - 1,
            1 << (index * 7),
            deterministic_random.randrange(1 << width),
            deterministic_random.randrange(1 << width),
        )
        for value in validation_values:
            if apply_matrix(inverse_rows, apply_matrix(rows, value)) != value:
                raise AssertionError("32x32 inverse validation failed")

    # Run the full selector for one nontrivial wide map.  This proves all
    # orientation transforms and post-passes at a width larger than emitted
    # FF8 circuits without making every --check invocation unnecessarily slow.
    selected, unused_name = synthesize_reversible_map_portfolio(
        representative_maps[0])
    if circuit_matrix(width, selected) != representative_maps[0]:
        raise AssertionError("32x32 portfolio selection mismatch")


def build_circuits():
    multiplication_matrices = tuple(
        multiplication_matrix(log_multiplier)
        for log_multiplier in range(FIELD_ORDER))
    multiplication_circuits = []
    forward_matrices = []
    inverse_matrices = []
    forward_circuits = []
    inverse_circuits = []
    multiplication_variants = []
    forward_variants = []
    inverse_variants = []

    # Exhaustive multiplier validation is fast and guards the column/row and
    # Cantor-coordinate conventions independently of the circuit synthesis.
    for log_multiplier, matrix in enumerate(multiplication_matrices):
        gates, variant = synthesize_reversible_map_portfolio(matrix)
        if circuit_matrix(WIRE_COUNT_MULTIPLY, gates) != matrix:
            raise AssertionError("incorrect multiplier circuit")
        for value in range(FIELD_ORDER):
            expected = scalar_multiply_log(value, log_multiplier)
            independent_expected = independent_scalar_multiply_log(
                value, log_multiplier)
            if expected != independent_expected:
                raise AssertionError("field table and polynomial references differ")
            if apply_matrix(matrix, value) != expected:
                raise AssertionError("incorrect multiplication matrix")
            if apply_circuit(value, gates) != expected:
                raise AssertionError("incorrect multiplication circuit")
        multiplication_circuits.append(gates)
        multiplication_variants.append(variant)

    deterministic_random = random.Random(0xFF8C1AC017)
    random_states = tuple(
        deterministic_random.randrange(1 << WIRE_COUNT_BUTTERFLY)
        for _ in range(64))
    validation_states = (
        (0, 0xFFFF)
        + tuple(1 << bit for bit in range(WIRE_COUNT_BUTTERFLY))
        + random_states)

    for skew in range(FIELD_ORDER):
        forward_function = lambda state, skew=skew: scalar_forward_butterfly(
            state, skew, multiplication_matrices)
        inverse_function = lambda state, skew=skew: scalar_inverse_butterfly(
            state, skew, multiplication_matrices)
        forward_matrix = linear_map(WIRE_COUNT_BUTTERFLY, forward_function)
        inverse_matrix = linear_map(WIRE_COUNT_BUTTERFLY, inverse_function)

        direct_forward = direct_butterfly_circuit(
            skew, False, multiplication_matrices)
        direct_inverse = direct_butterfly_circuit(
            skew, True, multiplication_matrices)
        chosen_forward, forward_variant = choose_butterfly_circuit(
            forward_matrix, direct_forward)
        chosen_inverse, inverse_variant = choose_butterfly_circuit(
            inverse_matrix, direct_inverse)

        for state in validation_states:
            expected_forward = forward_function(state)
            expected_inverse = inverse_function(state)
            if apply_circuit(state, chosen_forward) != expected_forward:
                raise AssertionError("incorrect forward butterfly circuit")
            if apply_circuit(state, chosen_inverse) != expected_inverse:
                raise AssertionError("incorrect inverse butterfly circuit")
            if apply_circuit(expected_forward, chosen_inverse) != state:
                raise AssertionError("inverse(forward(state)) failed")
            if apply_circuit(expected_inverse, chosen_forward) != state:
                raise AssertionError("forward(inverse(state)) failed")

        # Equality of the induced matrices proves the result for all 2^16
        # states, while the checks above independently exercise scalar order.
        if circuit_matrix(WIRE_COUNT_BUTTERFLY, chosen_forward) != forward_matrix:
            raise AssertionError("forward circuit map mismatch")
        if circuit_matrix(WIRE_COUNT_BUTTERFLY, chosen_inverse) != inverse_matrix:
            raise AssertionError("inverse circuit map mismatch")

        forward_matrices.append(forward_matrix)
        inverse_matrices.append(inverse_matrix)
        forward_circuits.append(chosen_forward)
        inverse_circuits.append(chosen_inverse)
        forward_variants.append(forward_variant)
        inverse_variants.append(inverse_variant)

    validate_representative_wide_synthesis()

    return {
        "multiply_matrices": tuple(multiplication_matrices),
        "forward_matrices": tuple(forward_matrices),
        "inverse_matrices": tuple(inverse_matrices),
        "multiply_circuits": tuple(multiplication_circuits),
        "forward_circuits": tuple(forward_circuits),
        "inverse_circuits": tuple(inverse_circuits),
        "multiply_variants": tuple(multiplication_variants),
        "forward_variants": tuple(forward_variants),
        "inverse_variants": tuple(inverse_variants),
    }


def circuit_checksum(circuit_data):
    digest = hashlib.sha256()
    digest.update(GENERATOR_VERSION)
    digest.update(bytes(LOG_LUT))
    digest.update(bytes(EXP_LUT))

    families = (
        (WIRE_COUNT_MULTIPLY, "multiply_matrices", "multiply_circuits"),
        (WIRE_COUNT_BUTTERFLY, "forward_matrices", "forward_circuits"),
        (WIRE_COUNT_BUTTERFLY, "inverse_matrices", "inverse_circuits"),
    )
    for width, matrix_name, circuit_name in families:
        digest.update(bytes((width,)))
        for rows, gates in zip(
                circuit_data[matrix_name], circuit_data[circuit_name]):
            for row in rows:
                digest.update(row.to_bytes(2, byteorder="little"))
            digest.update(len(gates).to_bytes(2, byteorder="little"))
            for destination, source in gates:
                digest.update(bytes((destination, source)))
    for variant_name in (
            "multiply_variants", "forward_variants", "inverse_variants"):
        for variant in circuit_data[variant_name]:
            encoded = variant.encode("ascii")
            digest.update(len(encoded).to_bytes(2, byteorder="little"))
            digest.update(encoded)
    return digest.hexdigest()


def format_integer_array(name, values):
    lines = ["static const uint8_t %s[256] = {" % name]
    for offset in range(0, len(values), 16):
        row = values[offset:offset + 16]
        lines.append("    " + ", ".join("%2d" % value for value in row) + ",")
    lines.append("};")
    return lines


def format_uint16_array(name, values):
    lines = ["static const uint16_t %s[256] = {" % name]
    for offset in range(0, len(values), 12):
        row = values[offset:offset + 12]
        lines.append("    " + ", ".join("%3d" % value for value in row) + ",")
    lines.append("};")
    return lines


def append_stats(lines, prefix, circuits, width):
    counts = [len(gates) for gates in circuits]
    depths = [circuit_depth(gates, width) for gates in circuits]
    lines.extend(format_integer_array("k%sGateCounts" % prefix, counts))
    lines.append("")
    lines.extend(format_integer_array("k%sDepths" % prefix, depths))
    lines.append("")
    lines.append("static const unsigned k%sMinGateCount = %d;" % (prefix, min(counts)))
    lines.append("static const unsigned k%sMaxGateCount = %d;" % (prefix, max(counts)))
    lines.append("static const double k%sAverageGateCount = %.9f;" % (
        prefix, float(sum(counts)) / len(counts)))
    lines.append("static const unsigned k%sMinDepth = %d;" % (prefix, min(depths)))
    lines.append("static const unsigned k%sMaxDepth = %d;" % (prefix, max(depths)))
    lines.append("static const double k%sAverageDepth = %.9f;" % (
        prefix, float(sum(depths)) / len(depths)))


def wire_names(width):
    if width == WIRE_COUNT_MULTIPLY:
        return tuple("x%d" % bit for bit in range(FIELD_BITS))
    return (
        tuple("x%d" % bit for bit in range(FIELD_BITS))
        + tuple("y%d" % bit for bit in range(FIELD_BITS)))


def append_circuit_specializations(lines, struct_name, circuits, width):
    names = wire_names(width)
    lines.append("template <unsigned Coefficient> struct %s;" % struct_name)
    lines.append("")
    for coefficient, gates in enumerate(circuits):
        lines.append("template <> struct %s<%d>" % (struct_name, coefficient))
        lines.append("{")
        lines.append("    template <typename Value>")
        lines.append("    static LEO_FORCE_INLINE void Apply(")
        for index, name in enumerate(names):
            comma = "," if index + 1 != len(names) else ")"
            lines.append("        Value& %s%s" % (name, comma))
        lines.append("    {")
        # This is compile-time-only and suppresses unused parameter warnings
        # for identity circuits without introducing payload instructions.
        for name in names:
            lines.append("        (void)%s;" % name)
        if gates:
            lines.append("")
        for destination, source in gates:
            destination_name = names[destination]
            source_name = names[source]
            lines.append("        %s = XorValue(%s, %s);" % (
                destination_name, destination_name, source_name))
        lines.append("    }")
        lines.append("};")
        lines.append("")


def generate_cpp(circuit_data):
    checksum = circuit_checksum(circuit_data)
    multiply_circuits = circuit_data["multiply_circuits"]
    forward_circuits = circuit_data["forward_circuits"]
    inverse_circuits = circuit_data["inverse_circuits"]
    variant_catalog = tuple(sorted(set(
        circuit_data["multiply_variants"] +
        circuit_data["forward_variants"] +
        circuit_data["inverse_variants"])))
    variant_ids = {name: index for index, name in enumerate(variant_catalog)}

    lines = [
        "// Generated by tools/generate_ff8_xor_circuits.py.  DO NOT EDIT.",
        "// Generator model: polynomial 0x11D; Cantor basis",
        "// { 1, 214, 152, 146, 86, 200, 88, 230 }.",
        "// Circuit checksum (SHA-256): %s" % checksum,
        "// CNOT depth forbids two gates in one layer from sharing a wire.",
        "// MultiplyCircuit<255> is identity; butterfly skew 255 omits multiply.",
        "",
        "#ifndef LEOPARD_FF8_XOR_CIRCUITS_GENERATED_INL",
        "#define LEOPARD_FF8_XOR_CIRCUITS_GENERATED_INL",
        "",
        "#include <stdint.h>",
        "",
        "#define LEOPARD_FF8XOR_FOR_EACH_LOG(M) \\",
    ]
    for coefficient in range(FIELD_ORDER):
        suffix = " \\" if coefficient + 1 != FIELD_ORDER else ""
        lines.append("    M(%d)%s" % (coefficient, suffix))

    lines.extend([
        "",
        "namespace leopard { namespace ff8xor { namespace generated {",
        "",
        "static const char kCircuitChecksum[] = \"%s\";" % checksum,
        "",
    ])
    lines.append("static const unsigned kSynthesisVariantCount = %d;" %
                 len(variant_catalog))
    lines.append(
        "static const char* const kSynthesisVariantNames[%d] = {" %
        len(variant_catalog))
    for name in variant_catalog:
        lines.append("    \"%s\"," % name)
    lines.append("};")
    lines.append("")
    lines.extend(format_uint16_array(
        "kMultiplyVariantIds",
        [variant_ids[name] for name in circuit_data["multiply_variants"]]))
    lines.append("")
    lines.extend(format_uint16_array(
        "kFFTVariantIds",
        [variant_ids[name] for name in circuit_data["forward_variants"]]))
    lines.append("")
    lines.extend(format_uint16_array(
        "kIFFTVariantIds",
        [variant_ids[name] for name in circuit_data["inverse_variants"]]))
    lines.append("")
    append_stats(
        lines, "Multiply", multiply_circuits, WIRE_COUNT_MULTIPLY)
    lines.append("")
    append_stats(
        lines, "FFT", forward_circuits, WIRE_COUNT_BUTTERFLY)
    lines.append("")
    append_stats(
        lines, "IFFT", inverse_circuits, WIRE_COUNT_BUTTERFLY)
    lines.append("")

    append_circuit_specializations(
        lines, "MultiplyCircuit", multiply_circuits, WIRE_COUNT_MULTIPLY)
    append_circuit_specializations(
        lines, "FFTCircuit", forward_circuits, WIRE_COUNT_BUTTERFLY)
    append_circuit_specializations(
        lines, "IFFTCircuit", inverse_circuits, WIRE_COUNT_BUTTERFLY)

    lines.extend([
        "}}} // namespace leopard::ff8xor::generated",
        "",
        "#endif // LEOPARD_FF8_XOR_CIRCUITS_GENERATED_INL",
        "",
    ])
    return "\n".join(lines)


def check_output(output_path, generated_text):
    try:
        current_text = output_path.read_text(encoding="utf-8")
    except FileNotFoundError:
        print("generated circuit file is missing: %s" % output_path, file=sys.stderr)
        return False

    if current_text == generated_text:
        print("FF8 XOR circuits are up to date: %s" % output_path)
        return True

    print("generated FF8 XOR circuits are stale: %s" % output_path, file=sys.stderr)
    difference = difflib.unified_diff(
        current_text.splitlines(),
        generated_text.splitlines(),
        fromfile=str(output_path),
        tofile="regenerated",
        lineterm="")
    for line in difference:
        print(line, file=sys.stderr)
    return False


def print_synthesis_summary(circuit_data, elapsed_seconds):
    families = (
        ("Multiply", "multiply_circuits", "multiply_variants",
         WIRE_COUNT_MULTIPLY),
        ("FFT", "forward_circuits", "forward_variants",
         WIRE_COUNT_BUTTERFLY),
        ("IFFT", "inverse_circuits", "inverse_variants",
         WIRE_COUNT_BUTTERFLY),
    )
    for label, circuit_name, variant_name, width in families:
        circuits = circuit_data[circuit_name]
        gates = [len(circuit) for circuit in circuits]
        depths = [circuit_depth(circuit, width) for circuit in circuits]
        print("%s synthesis: gates=%d (min=%d max=%d avg=%.6f) "
              "depth=%d (min=%d max=%d avg=%.6f) variants=%d" % (
                  label, sum(gates), min(gates), max(gates),
                  float(sum(gates)) / len(gates),
                  sum(depths), min(depths), max(depths),
                  float(sum(depths)) / len(depths),
                  len(set(circuit_data[variant_name]))))
    print("Synthesis generation time: %.3f seconds" % elapsed_seconds)


def parse_arguments():
    repository_root = Path(__file__).resolve().parents[1]
    default_output = repository_root / "generated" / "LeopardFF8XorCircuits.inl"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify that the checked-in generated source is current")
    parser.add_argument(
        "--output",
        type=Path,
        default=default_output,
        help="generated output path (default: %(default)s)")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    begin = time.monotonic()
    circuit_data = build_circuits()
    generated_text = generate_cpp(circuit_data)
    elapsed_seconds = time.monotonic() - begin

    if arguments.check:
        result = 0 if check_output(arguments.output, generated_text) else 1
        print_synthesis_summary(circuit_data, elapsed_seconds)
        return result

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(generated_text, encoding="utf-8")
    print("Generated %s" % arguments.output)
    print("Circuit checksum: %s" % circuit_checksum(circuit_data))
    print_synthesis_summary(circuit_data, elapsed_seconds)
    return 0


if __name__ == "__main__":
    sys.exit(main())
