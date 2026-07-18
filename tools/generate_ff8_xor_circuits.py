#!/usr/bin/env python3
"""Generate the static FF8 XOR circuits used by the experimental backend.

The generated C++ deliberately contains named wire variables and literal XOR
statements.  The payload kernels include it after defining XorValue() for the
selected scalar or SIMD value type.
"""

from __future__ import print_function

import argparse
import difflib
import hashlib
import random
import sys
from pathlib import Path


FIELD_BITS = 8
FIELD_ORDER = 1 << FIELD_BITS
FIELD_MODULUS = FIELD_ORDER - 1
FIELD_POLYNOMIAL = 0x11D
CANTOR_BASIS = (1, 214, 152, 146, 86, 200, 88, 230)
WIRE_COUNT_MULTIPLY = 8
WIRE_COUNT_BUTTERFLY = 16
GENERATOR_VERSION = b"LeopardFF8XorCircuits-v1\0"


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
    synthesized_gates = synthesize_reversible_map(rows)
    if circuit_matrix(len(rows), synthesized_gates) != rows:
        raise AssertionError("synthesized butterfly circuit is incorrect")
    if circuit_matrix(len(rows), direct_gates) != rows:
        raise AssertionError("direct butterfly circuit is incorrect")

    # The tuple comparison exactly implements gate count, then dependency
    # depth, then lexical gate ordering.  Identical gate lists are equivalent.
    return min(
        (synthesized_gates, direct_gates),
        key=lambda gates: circuit_key(gates, len(rows)))


def build_circuits():
    multiplication_matrices = tuple(
        multiplication_matrix(log_multiplier)
        for log_multiplier in range(FIELD_ORDER))
    multiplication_circuits = []
    forward_matrices = []
    inverse_matrices = []
    forward_circuits = []
    inverse_circuits = []

    # Exhaustive multiplier validation is fast and guards the column/row and
    # Cantor-coordinate conventions independently of the circuit synthesis.
    for log_multiplier, matrix in enumerate(multiplication_matrices):
        gates = synthesize_reversible_map(matrix)
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
        chosen_forward = choose_butterfly_circuit(
            forward_matrix, direct_forward)
        chosen_inverse = choose_butterfly_circuit(
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

    return {
        "multiply_matrices": tuple(multiplication_matrices),
        "forward_matrices": tuple(forward_matrices),
        "inverse_matrices": tuple(inverse_matrices),
        "multiply_circuits": tuple(multiplication_circuits),
        "forward_circuits": tuple(forward_circuits),
        "inverse_circuits": tuple(inverse_circuits),
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
    return digest.hexdigest()


def format_integer_array(name, values):
    lines = ["static const uint8_t %s[256] = {" % name]
    for offset in range(0, len(values), 16):
        row = values[offset:offset + 16]
        lines.append("    " + ", ".join("%2d" % value for value in row) + ",")
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
    circuit_data = build_circuits()
    generated_text = generate_cpp(circuit_data)

    if arguments.check:
        return 0 if check_output(arguments.output, generated_text) else 1

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(generated_text, encoding="utf-8")
    print("Generated %s" % arguments.output)
    print("Circuit checksum: %s" % circuit_checksum(circuit_data))
    return 0


if __name__ == "__main__":
    sys.exit(main())
