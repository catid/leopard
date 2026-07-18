#!/usr/bin/env python3
"""Deterministic XOR2/XOR3 scheduling helpers for generated FF8 circuits.

The checked-in FF8 circuit generator represents a CNOT as ``(dst, src)``.
AVX-512 can combine two adjacent CNOTs with the same destination into one
VPTERNLOG operation::

    dst ^= src_a
    dst ^= src_b

becomes::

    dst = dst XOR src_a XOR src_b

This module deliberately has no dependency on the generator.  It provides a
small, independently testable scheduling layer that the generator can import
once its synthesis portfolio is stable.  Schedules remain tuples of literal
named-wire operations; this is not a runtime gate interpreter.
"""

from __future__ import print_function


OP_XOR2 = 2
OP_XOR3 = 3

VARIANT_XOR2 = 0
VARIANT_ADJACENT = 1
VARIANT_COMMUTE_BACKWARD = 2
VARIANT_COMMUTE_FORWARD = 3

VARIANT_NAMES = (
    "xor2",
    "adjacent",
    "commute-backward",
    "commute-forward",
)

# These are intentionally labelled estimates.  Exact encodings depend on the
# chosen registers and compiler, while these stable values are useful as a
# deterministic code-size tie-break and metadata field.
DEFAULT_XOR2_CODE_BYTES = 6
DEFAULT_XOR3_CODE_BYTES = 7


def _validate_width(width):
    if not isinstance(width, int) or width <= 0:
        raise ValueError("wire width must be a positive integer")


def validate_gates(gates, width):
    """Return *gates* as a validated immutable CNOT sequence."""
    _validate_width(width)
    result = tuple(tuple(gate) for gate in gates)
    for gate in result:
        if len(gate) != 2:
            raise ValueError("a CNOT gate must contain destination and source")
        destination, source = gate
        if not isinstance(destination, int) or not isinstance(source, int):
            raise ValueError("wire indices must be integers")
        if not (0 <= destination < width and 0 <= source < width):
            raise ValueError("CNOT wire index is outside the circuit")
        if destination == source:
            raise ValueError("a CNOT source and destination must differ")
    return result


def validate_operations(operations, width):
    """Return a validated immutable mixed XOR2/XOR3 operation sequence."""
    _validate_width(width)
    result = tuple(tuple(operation) for operation in operations)
    for operation in result:
        if not operation:
            raise ValueError("empty circuit operation")
        opcode = operation[0]
        expected_length = 3 if opcode == OP_XOR2 else 4 if opcode == OP_XOR3 else 0
        if not expected_length or len(operation) != expected_length:
            raise ValueError("invalid XOR2/XOR3 operation")
        destination = operation[1]
        sources = operation[2:]
        if not isinstance(destination, int) or not all(
                isinstance(source, int) for source in sources):
            raise ValueError("wire indices must be integers")
        if not (0 <= destination < width) or not all(
                0 <= source < width for source in sources):
            raise ValueError("operation wire index is outside the circuit")
        if destination in sources:
            raise ValueError("an XOR source may not equal its destination")
        if opcode == OP_XOR3 and sources[0] == sources[1]:
            raise ValueError("an XOR3 operation must have distinct sources")
    return result


def cnot_gates_commute(first, second):
    """Return whether two CNOT row additions can exchange positions."""
    destination_a, source_a = first
    destination_b, source_b = second
    return destination_a != source_b and destination_b != source_a


def xor2_schedule(gates, width):
    """Convert CNOTs to the canonical XOR2-only schedule."""
    return tuple((OP_XOR2, destination, source)
                 for destination, source in validate_gates(gates, width))


def _combined_operation(first, second):
    """Combine adjacent same-destination CNOTs, or cancel duplicates."""
    destination_a, source_a = first
    destination_b, source_b = second
    if destination_a != destination_b:
        raise ValueError("only equal-destination CNOTs can be combined")
    if source_a == source_b:
        return None
    # XOR is commutative.  Canonical source order ensures that schedules found
    # through different legal commuting directions compare and hash the same.
    source_a, source_b = sorted((source_a, source_b))
    return (OP_XOR3, destination_a, source_a, source_b)


def pair_adjacent(gates, width):
    """Greedily combine already adjacent equal-destination CNOTs."""
    pending = validate_gates(gates, width)
    result = []
    index = 0
    while index < len(pending):
        first = pending[index]
        if index + 1 < len(pending) and first[0] == pending[index + 1][0]:
            combined = _combined_operation(first, pending[index + 1])
            if combined is not None:
                result.append(combined)
            index += 2
        else:
            result.append((OP_XOR2, first[0], first[1]))
            index += 1
    return tuple(result)


def pair_commuting_backward(gates, width):
    """Move a later CNOT backward across commuting gates, then combine it.

    The earliest legal same-destination partner is selected.  Moving the
    partner, rather than the first gate, gives a deterministic schedule and
    preserves semantics because every crossed gate is checked explicitly.
    """
    pending = list(validate_gates(gates, width))
    result = []
    while pending:
        first = pending.pop(0)
        partner_index = None
        for index, candidate in enumerate(pending):
            if candidate[0] != first[0]:
                continue
            if all(cnot_gates_commute(candidate, crossed)
                   for crossed in pending[:index]):
                partner_index = index
                break

        if partner_index is None:
            result.append((OP_XOR2, first[0], first[1]))
            continue

        partner = pending.pop(partner_index)
        combined = _combined_operation(first, partner)
        if combined is not None:
            result.append(combined)
    return tuple(result)


def pair_commuting_forward(gates, width):
    """Move an earlier CNOT forward across commuting gates, then combine it."""
    validated = validate_gates(gates, width)

    def schedule(sequence):
        if not sequence:
            return ()

        first = sequence[0]
        partner_index = None
        for index in range(1, len(sequence)):
            candidate = sequence[index]
            if candidate[0] != first[0]:
                continue
            if all(cnot_gates_commute(first, crossed)
                   for crossed in sequence[1:index]):
                partner_index = index
                break

        if partner_index is None:
            return ((OP_XOR2, first[0], first[1]),) + schedule(sequence[1:])

        middle = sequence[1:partner_index]
        partner = sequence[partner_index]
        combined = _combined_operation(first, partner)
        pair = () if combined is None else (combined,)
        return schedule(middle) + pair + schedule(sequence[partner_index + 1:])

    return schedule(validated)


def expand_operations(operations, width):
    """Expand a mixed schedule back to its equivalent CNOT sequence."""
    expanded = []
    for operation in validate_operations(operations, width):
        if operation[0] == OP_XOR2:
            expanded.append((operation[1], operation[2]))
        else:
            expanded.append((operation[1], operation[2]))
            expanded.append((operation[1], operation[3]))
    return tuple(expanded)


def apply_gates(state, gates, width):
    """Apply a CNOT sequence to an integer bit-vector."""
    validated = validate_gates(gates, width)
    if not isinstance(state, int) or not (0 <= state < (1 << width)):
        raise ValueError("state is outside the circuit width")
    for destination, source in validated:
        state ^= ((state >> source) & 1) << destination
    return state


def apply_operations(state, operations, width):
    """Apply a mixed XOR2/XOR3 schedule to an integer bit-vector."""
    validated = validate_operations(operations, width)
    if not isinstance(state, int) or not (0 <= state < (1 << width)):
        raise ValueError("state is outside the circuit width")
    for operation in validated:
        destination = operation[1]
        value = (state >> operation[2]) & 1
        if operation[0] == OP_XOR3:
            value ^= (state >> operation[3]) & 1
        state ^= value << destination
    return state


def circuit_matrix(gates, width):
    """Return row masks for the binary linear map implemented by CNOTs."""
    rows = [1 << row for row in range(width)]
    for destination, source in validate_gates(gates, width):
        rows[destination] ^= rows[source]
    return tuple(rows)


def operation_matrix(operations, width):
    """Return row masks for the binary linear map implemented by a schedule."""
    rows = [1 << row for row in range(width)]
    for operation in validate_operations(operations, width):
        destination = operation[1]
        rows[destination] ^= rows[operation[2]]
        if operation[0] == OP_XOR3:
            rows[destination] ^= rows[operation[3]]
    return tuple(rows)


def verify_schedule(gates, operations, width):
    """Raise if a mixed schedule does not implement the original CNOT map."""
    if circuit_matrix(gates, width) != operation_matrix(operations, width):
        raise AssertionError("XOR3 schedule does not match its CNOT circuit")
    return True


def schedule_depth(operations, width):
    """Count conservative dependency layers, forbidding shared live wires."""
    last_layer = [0] * width
    maximum = 0
    for operation in validate_operations(operations, width):
        wires = operation[1:]
        layer = max(last_layer[wire] for wire in wires) + 1
        for wire in wires:
            last_layer[wire] = layer
        maximum = max(maximum, layer)
    return maximum


def schedule_metadata(operations, width,
                      xor2_code_bytes=DEFAULT_XOR2_CODE_BYTES,
                      xor3_code_bytes=DEFAULT_XOR3_CODE_BYTES):
    """Return deterministic cost and diagnostic metadata for a schedule.

    Peak liveness is the circuit width: these are reversible in-place maps and
    every named input wire remains a named output wire that must be stored.
    """
    validated = validate_operations(operations, width)
    xor2_count = sum(operation[0] == OP_XOR2 for operation in validated)
    xor3_count = sum(operation[0] == OP_XOR3 for operation in validated)
    return {
        "xor2_count": xor2_count,
        "xor3_count": xor3_count,
        "instruction_count": xor2_count + xor3_count,
        "depth": schedule_depth(validated, width),
        "peak_live_wires": width,
        "estimated_code_bytes": (
            xor2_count * xor2_code_bytes + xor3_count * xor3_code_bytes),
    }


def schedule_key(operations, width,
                 xor2_code_bytes=DEFAULT_XOR2_CODE_BYTES,
                 xor3_code_bytes=DEFAULT_XOR3_CODE_BYTES):
    """Return the stable AVX-512 schedule preference key."""
    validated = validate_operations(operations, width)
    metadata = schedule_metadata(
        validated, width, xor2_code_bytes, xor3_code_bytes)
    return (
        metadata["instruction_count"],
        metadata["depth"],
        metadata["peak_live_wires"],
        metadata["estimated_code_bytes"],
        validated,
    )


def schedule_variants(gates, width):
    """Build and independently verify the deterministic scheduling portfolio."""
    validated = validate_gates(gates, width)
    variants = (
        (VARIANT_XOR2, xor2_schedule(validated, width)),
        (VARIANT_ADJACENT, pair_adjacent(validated, width)),
        (VARIANT_COMMUTE_BACKWARD,
         pair_commuting_backward(validated, width)),
        (VARIANT_COMMUTE_FORWARD,
         pair_commuting_forward(validated, width)),
    )
    for unused_variant, operations in variants:
        verify_schedule(validated, operations, width)
    return variants


def choose_schedule(gates, width,
                    xor2_code_bytes=DEFAULT_XOR2_CODE_BYTES,
                    xor3_code_bytes=DEFAULT_XOR3_CODE_BYTES):
    """Return ``(variant_id, operations)`` for the best verified schedule."""
    variants = schedule_variants(gates, width)
    return min(
        variants,
        key=lambda item: schedule_key(
            item[1], width, xor2_code_bytes, xor3_code_bytes) + (item[0],))
