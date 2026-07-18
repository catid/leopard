#!/usr/bin/env python3
"""Generate reachable AVX-512 four-buffer FF8 XOR circuits.

This is deliberately separate from the two-buffer circuit corpus.  It proves
and emits the two-layer radix-4 maps which an AVX-512 kernel can keep in its 32
architectural vector registers.  The generated payload bodies contain only
literal operations on named wires; tuple selection and dispatch belong outside
the vector-chunk loop.
"""

from __future__ import print_function

import argparse
import collections
import difflib
import hashlib
import json
import random
import sys
import time
from pathlib import Path


TOOLS_DIRECTORY = Path(__file__).resolve().parent
if str(TOOLS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIRECTORY))

# Reuse the independently tested field model and two-buffer synthesis
# portfolio, but construct and verify every 32-wire map in this file.
import generate_ff8_xor_circuits as two_buffer  # noqa: E402
import ff8_xor3_schedule as xor3  # noqa: E402


FIELD_BITS = 8
FIELD_ORDER = 256
FIELD_MODULUS = 255
WIRE_COUNT = 32
GENERATOR_ID = "ff8-four-buffer-v1"
RANDOM_SEED = 0x4FF8C17
RANDOM_STATE_COUNT = 128
WIRE_NAMES = tuple(
    "%s%d" % (buffer_name, bit)
    for buffer_name in ("a", "b", "c", "d")
    for bit in range(FIELD_BITS))


def make_fft_skew_padded():
    """Reproduce LeopardFF8Xor.cpp's padded FFT skew table exactly."""
    temporary = [1 << bit for bit in range(1, FIELD_BITS)]
    skew = [FIELD_MODULUS] + [0] * FIELD_MODULUS

    def get(logical_index):
        return skew[logical_index + 1]

    def put(logical_index, value):
        skew[logical_index + 1] = value

    for level in range(FIELD_BITS - 1):
        step = 1 << (level + 1)
        put((1 << level) - 1, 0)
        for bit in range(level, FIELD_BITS - 1):
            span = 1 << (bit + 1)
            for index in range((1 << level) - 1, span, step):
                put(index + span, get(index) ^ temporary[bit])

        temporary[level] = FIELD_MODULUS - two_buffer.LOG_LUT[
            two_buffer.scalar_multiply_log(
                temporary[level],
                two_buffer.LOG_LUT[temporary[level] ^ 1])]
        for bit in range(level + 1, FIELD_BITS - 1):
            logarithm = two_buffer.add_mod(
                two_buffer.LOG_LUT[temporary[bit] ^ 1], temporary[level])
            temporary[bit] = two_buffer.scalar_multiply_log(
                temporary[bit], logarithm)

    for index in range(FIELD_MODULUS):
        put(index, two_buffer.LOG_LUT[get(index)])
    return tuple(skew)


FFT_SKEW_PADDED = make_fft_skew_padded()


def next_power_of_two(value):
    result = 1
    while result < value:
        result <<= 1
    return result


def record_ifft_tuples(counter, count, count_truncated, skew_base):
    """Record radix-4 IFFT calls, including their inner-distance frequency."""
    distance = 1
    while distance * 4 <= count:
        span = distance * 4
        for range_start in range(0, count_truncated, span):
            # The following distance*2 layer is scheduled for this range as
            # soon as range_start is live, but a complete radix-4 unit also
            # requires the second distance-wide pair from this layer.  When
            # it is truncated away, production must retain the natural
            # two-way schedule instead of applying a four-buffer map that
            # includes an extra (c,d) butterfly.
            if range_start + 2 * distance >= count_truncated:
                continue
            index = skew_base + range_start
            tuple_value = (
                FFT_SKEW_PADDED[index + distance],
                FFT_SKEW_PADDED[index + 3 * distance],
                FFT_SKEW_PADDED[index + 2 * distance],
            )
            counter[tuple_value] += distance
        distance *= 4


def record_fft_tuples(counter, count, count_truncated, skew_base):
    """Record radix-4 FFT calls, including their inner-distance frequency."""
    distance4 = count
    distance = count >> 2
    while distance:
        for range_start in range(0, count_truncated, distance4):
            # Both distance-wide pairs in the lower layer must execute before
            # the preceding distance*2 layer can be represented by one
            # complete radix-4 circuit.
            if range_start + 2 * distance >= count_truncated:
                continue
            index = skew_base + range_start
            tuple_value = (
                FFT_SKEW_PADDED[index + distance],
                FFT_SKEW_PADDED[index + 3 * distance],
                FFT_SKEW_PADDED[index + 2 * distance],
            )
            counter[tuple_value] += distance
        distance4 = distance
        distance >>= 2


def enumerate_reachable_tuples():
    """Enumerate every tuple reachable by a valid public FF8 API shape.

    The loop deliberately ranges over all positive k/r pairs whose padded
    decoder transform fits in 256 elements.  Encode chunk skew bases, encode
    output truncation, and the decoder's full n-point transforms are all
    represented.  Frequencies count fused whole-buffer calls in this complete
    parameter census; they are provenance, not a real-world workload model.
    """
    fft_counts = collections.Counter()
    ifft_counts = collections.Counter()
    parameter_count = 0

    for original_count in range(1, FIELD_ORDER):
        for recovery_count in range(1, FIELD_ORDER + 1 - original_count):
            m = next_power_of_two(recovery_count)
            n = next_power_of_two(m + original_count)
            if n > FIELD_ORDER:
                continue
            parameter_count += 1

            for chunk_start in range(0, original_count, m):
                chunk_count = min(m, original_count - chunk_start)
                record_ifft_tuples(
                    ifft_counts, m, chunk_count, m + chunk_start)
            record_fft_tuples(fft_counts, m, recovery_count, 0)

            record_ifft_tuples(ifft_counts, n, m + original_count, 0)
            record_fft_tuples(fft_counts, n, m + original_count, 0)

    tuples = tuple(sorted(set(fft_counts) | set(ifft_counts)))
    if len(tuples) != 64:
        raise AssertionError(
            "expected 64 reachable skew tuples, observed %d" % len(tuples))
    if set(fft_counts) != set(ifft_counts):
        raise AssertionError("reachable FFT and IFFT tuple sets differ")
    if [value for value in tuples if FIELD_MODULUS in value] != [
            (FIELD_MODULUS, 85, FIELD_MODULUS)]:
        raise AssertionError("unexpected sentinel-bearing tuple set")
    return {
        "parameter_count": parameter_count,
        "tuples": tuples,
        "fft_counts": fft_counts,
        "ifft_counts": ifft_counts,
    }


def scalar_butterfly(left, right, skew, inverse, matrices):
    state = left | (right << FIELD_BITS)
    if inverse:
        state = two_buffer.scalar_inverse_butterfly(state, skew, matrices)
    else:
        state = two_buffer.scalar_forward_butterfly(state, skew, matrices)
    return state & 0xFF, (state >> FIELD_BITS) & 0xFF


def scalar_four_buffer(state, tuple_value, inverse, matrices):
    values = [(state >> (FIELD_BITS * index)) & 0xFF
              for index in range(4)]

    def apply(left_index, right_index, skew):
        values[left_index], values[right_index] = scalar_butterfly(
            values[left_index], values[right_index], skew, inverse, matrices)

    skew01, skew23, skew02 = tuple_value
    if inverse:
        apply(0, 1, skew01)
        apply(2, 3, skew23)
        apply(0, 2, skew02)
        apply(1, 3, skew02)
    else:
        apply(0, 2, skew02)
        apply(1, 3, skew02)
        apply(0, 1, skew01)
        apply(2, 3, skew23)
    return sum(value << (FIELD_BITS * index)
               for index, value in enumerate(values))


def map_two_buffer_gates(gates, left_index, right_index):
    offsets = (FIELD_BITS * left_index, FIELD_BITS * right_index)
    result = []
    for destination, source in gates:
        destination_half = 0 if destination < FIELD_BITS else 1
        source_half = 0 if source < FIELD_BITS else 1
        result.append((
            offsets[destination_half] + destination % FIELD_BITS,
            offsets[source_half] + source % FIELD_BITS,
        ))
    return tuple(result)


def four_buffer_sequence(tuple_value, inverse):
    skew01, skew23, skew02 = tuple_value
    if inverse:
        return (
            (0, 1, skew01), (2, 3, skew23),
            (0, 2, skew02), (1, 3, skew02),
        )
    return (
        (0, 2, skew02), (1, 3, skew02),
        (0, 1, skew01), (2, 3, skew23),
    )


def direct_four_buffer_circuit(tuple_value, inverse, matrices):
    result = []
    for left, right, skew in four_buffer_sequence(tuple_value, inverse):
        gates = two_buffer.direct_butterfly_circuit(
            skew, inverse, matrices)
        result.extend(map_two_buffer_gates(gates, left, right))
    return tuple(result)


def build_two_buffer_portfolio(skews, inverse, matrices):
    """Select the established verified two-buffer circuits for composition."""
    result = {}
    for skew in skews:
        function = (lambda state, skew=skew:
                    two_buffer.scalar_inverse_butterfly(
                        state, skew, matrices)) if inverse else (
                    lambda state, skew=skew:
                    two_buffer.scalar_forward_butterfly(
                        state, skew, matrices))
        rows = two_buffer.linear_map(16, function)
        direct = two_buffer.direct_butterfly_circuit(
            skew, inverse, matrices)
        gates, variant = two_buffer.choose_butterfly_circuit(
            rows, direct)
        if two_buffer.circuit_matrix(16, gates) != rows:
            raise AssertionError("two-buffer portfolio map mismatch")
        result[skew] = (gates, variant)
    return result


def composed_four_buffer_circuit(
        tuple_value, inverse, two_buffer_portfolio):
    result = []
    variants = []
    for left, right, skew in four_buffer_sequence(tuple_value, inverse):
        gates, variant = two_buffer_portfolio[skew]
        result.extend(map_two_buffer_gates(gates, left, right))
        variants.append(variant)
    return tuple(result), tuple(variants)


def whole_map_candidates(rows):
    """Return deterministic full-32x32 Gauss-Jordan synthesis candidates."""
    inverse_rows = two_buffer.inverse_matrix(rows)
    transpose_gates = two_buffer.synthesize_reversible_map(
        two_buffer.transpose_matrix(rows))
    inverse_transpose_gates = two_buffer.synthesize_reversible_map(
        two_buffer.transpose_matrix(inverse_rows))
    candidates = (
        ("whole-map/direct",
         two_buffer.synthesize_reversible_map(rows)),
        ("whole-map/inverse",
         tuple(reversed(two_buffer.synthesize_reversible_map(inverse_rows)))),
        ("whole-map/transpose",
         two_buffer.transpose_circuit(transpose_gates)),
        ("whole-map/inverse-transpose",
         tuple(reversed(two_buffer.transpose_circuit(
             inverse_transpose_gates)))),
    )
    result = []
    for name, gates in candidates:
        gates = two_buffer.cancel_commuting_duplicate_gates(gates)
        if two_buffer.circuit_matrix(WIRE_COUNT, gates) != rows:
            raise AssertionError("%s map mismatch" % name)
        result.append((name, gates))
    return tuple(result)


def operation_literal(operation):
    destination = WIRE_NAMES[operation[1]]
    if operation[0] == xor3.OP_XOR2:
        source = WIRE_NAMES[operation[2]]
        return "%s = XorValue(%s, %s);" % (
            destination, destination, source)
    source_a = WIRE_NAMES[operation[2]]
    source_b = WIRE_NAMES[operation[3]]
    return "%s = Xor3Value(%s, %s, %s);" % (
        destination, destination, source_a, source_b)


def operation_metadata(operations):
    metadata = xor3.schedule_metadata(operations, WIRE_COUNT)
    metadata["literal_statement_bytes"] = sum(
        len(operation_literal(operation).encode("ascii")) + 9
        for operation in operations)
    return metadata


def candidate_record(name, gates):
    if two_buffer.circuit_matrix(WIRE_COUNT, gates) != \
            xor3.circuit_matrix(gates, WIRE_COUNT):
        raise AssertionError("independent CNOT evaluators disagree")
    xor2_operations = xor3.xor2_schedule(gates, WIRE_COUNT)
    xor3_variant, xor3_operations = xor3.choose_schedule(gates, WIRE_COUNT)
    xor3.verify_schedule(gates, xor3_operations, WIRE_COUNT)
    return {
        "name": name,
        "gates": tuple(gates),
        "gate_count": len(gates),
        "cnot_depth": two_buffer.circuit_depth(gates, WIRE_COUNT),
        "xor2_operations": xor2_operations,
        "xor2_metadata": operation_metadata(xor2_operations),
        "xor3_variant": xor3.VARIANT_NAMES[xor3_variant],
        "xor3_operations": xor3_operations,
        "xor3_metadata": operation_metadata(xor3_operations),
    }


def xor2_candidate_key(candidate):
    return (
        candidate["xor2_metadata"]["instruction_count"],
        candidate["xor2_metadata"]["depth"],
        candidate["xor2_metadata"]["estimated_code_bytes"],
        candidate["gates"],
        candidate["name"],
    )


def xor3_candidate_key(candidate):
    operations = candidate["xor3_operations"]
    return xor3.schedule_key(operations, WIRE_COUNT) + (
        candidate["gate_count"], candidate["name"])


def validate_selected_circuits(
        rows, inverse_rows, xor2_candidate, xor3_candidate,
        tuple_value, inverse, matrices, random_states):
    if xor3.operation_matrix(
            xor2_candidate["xor2_operations"], WIRE_COUNT) != rows:
        raise AssertionError("selected XOR2 map mismatch")
    if xor3.operation_matrix(
            xor3_candidate["xor3_operations"], WIRE_COUNT) != rows:
        raise AssertionError("selected XOR3 map mismatch")

    states = (
        (0, (1 << WIRE_COUNT) - 1)
        + tuple(1 << bit for bit in range(WIRE_COUNT))
        + random_states)
    for state in states:
        expected = scalar_four_buffer(
            state, tuple_value, inverse, matrices)
        if two_buffer.apply_matrix(rows, state) != expected:
            raise AssertionError("32-wire matrix does not match scalar map")
        if xor3.apply_operations(
                state, xor2_candidate["xor2_operations"], WIRE_COUNT) != expected:
            raise AssertionError("XOR2 circuit does not match scalar map")
        if xor3.apply_operations(
                state, xor3_candidate["xor3_operations"], WIRE_COUNT) != expected:
            raise AssertionError("XOR3 circuit does not match scalar map")
        inverse_value = scalar_four_buffer(
            expected, tuple_value, not inverse, matrices)
        if inverse_value != state:
            raise AssertionError("four-buffer FFT/IFFT scalar inverse failed")
        if two_buffer.apply_matrix(inverse_rows, expected) != state:
            raise AssertionError("four-buffer matrix inverse failed")


def build_artifact():
    reachability = enumerate_reachable_tuples()
    tuples = reachability["tuples"]
    matrices = tuple(two_buffer.multiplication_matrix(skew)
                     for skew in range(FIELD_ORDER))
    skews = tuple(sorted(set(value for item in tuples for value in item)))
    portfolios = {
        False: build_two_buffer_portfolio(skews, False, matrices),
        True: build_two_buffer_portfolio(skews, True, matrices),
    }
    deterministic_random = random.Random(RANDOM_SEED)
    random_states = tuple(
        deterministic_random.randrange(1 << WIRE_COUNT)
        for unused in range(RANDOM_STATE_COUNT))

    families = {"fft": [], "ifft": []}
    for tuple_index, tuple_value in enumerate(tuples):
        direction_rows = {}
        for inverse, family_name in ((False, "fft"), (True, "ifft")):
            rows = two_buffer.linear_map(
                WIRE_COUNT,
                lambda state, tuple_value=tuple_value, inverse=inverse:
                scalar_four_buffer(state, tuple_value, inverse, matrices))
            direction_rows[inverse] = rows
            direct = direct_four_buffer_circuit(
                tuple_value, inverse, matrices)
            composed, component_variants = composed_four_buffer_circuit(
                tuple_value, inverse, portfolios[inverse])
            candidates = [
                candidate_record("direct", direct),
                candidate_record("composed-two-buffer-portfolio", composed),
            ]
            for name, gates in whole_map_candidates(rows):
                candidates.append(candidate_record(name, gates))
            for candidate in candidates:
                if xor3.circuit_matrix(
                        candidate["gates"], WIRE_COUNT) != rows:
                    raise AssertionError("candidate does not implement tuple map")

            xor2_selected = min(candidates, key=xor2_candidate_key)
            xor3_selected = min(candidates, key=xor3_candidate_key)
            record = {
                "tuple_index": tuple_index,
                "tuple": tuple_value,
                "rows": rows,
                "component_variants": component_variants,
                "candidates": candidates,
                "xor2_selected": xor2_selected,
                "xor3_selected": xor3_selected,
            }
            families[family_name].append(record)

        fft_rows = direction_rows[False]
        if two_buffer.inverse_matrix(fft_rows) != direction_rows[True]:
            raise AssertionError("FFT/IFFT full maps are not inverses")
        validate_selected_circuits(
            direction_rows[False], direction_rows[True],
            families["fft"][-1]["xor2_selected"],
            families["fft"][-1]["xor3_selected"],
            tuple_value, False, matrices, random_states)
        validate_selected_circuits(
            direction_rows[True], direction_rows[False],
            families["ifft"][-1]["xor2_selected"],
            families["ifft"][-1]["xor3_selected"],
            tuple_value, True, matrices, random_states)

    return {
        "generator": GENERATOR_ID,
        "reachability": reachability,
        "families": families,
        "skew_count": len(skews),
        "random_seed": RANDOM_SEED,
        "random_state_count": RANDOM_STATE_COUNT,
    }


def checksum_payload(artifact):
    payload = {
        "generator": artifact["generator"],
        "tuples": artifact["reachability"]["tuples"],
        "families": {
            family_name: [
                {
                    "tuple": record["tuple"],
                    "rows": record["rows"],
                    "xor2_name": record["xor2_selected"]["name"],
                    "xor2": record["xor2_selected"]["xor2_operations"],
                    "xor3_name": record["xor3_selected"]["name"],
                    "xor3": record["xor3_selected"]["xor3_operations"],
                }
                for record in records]
            for family_name, records in artifact["families"].items()
        },
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def append_uint_array(lines, declaration, values, columns=16):
    lines.append("static const uint16_t %s[%d] = {" % (declaration, len(values)))
    for start in range(0, len(values), columns):
        lines.append("    " + ", ".join(
            str(value) for value in values[start:start + columns]) + ",")
    lines.append("};")


def append_circuit(lines, struct_name, tuple_index, operations):
    lines.append("template <> struct %s<%d>" % (struct_name, tuple_index))
    lines.append("{")
    lines.append("    template <typename Value>")
    lines.append("    static LEO_FORCE_INLINE void Apply(")
    for index, name in enumerate(WIRE_NAMES):
        suffix = "," if index + 1 != len(WIRE_NAMES) else ")"
        lines.append("        Value& %s%s" % (name, suffix))
    lines.append("    {")
    for name in WIRE_NAMES:
        lines.append("        (void)%s;" % name)
    lines.append("")
    for operation in operations:
        lines.append("        " + operation_literal(operation))
    lines.append("    }")
    lines.append("};")
    lines.append("")


def generate_cpp(artifact):
    checksum = checksum_payload(artifact)
    tuples = artifact["reachability"]["tuples"]
    lines = [
        "// Generated by tools/generate_ff8_xor_four_buffer_circuits.py.  DO NOT EDIT.",
        "// Reachable FF8 radix-4 tuple count: %d." % len(tuples),
        "// Circuit checksum (SHA-256): %s" % checksum,
        "// Xor3Value must implement parity truth table immediate 0x96.",
        "// Tuple dispatch must occur outside the vector-chunk loop.",
        "// Use destructive +v XOR wrappers for 32-live-ZMM integration;",
        "// tools/inspect_ff8xor_four_buffer_avx512.py proves spill status.",
        "",
        "#ifndef LEOPARD_FF8_XOR_FOUR_BUFFER_CIRCUITS_GENERATED_INL",
        "#define LEOPARD_FF8_XOR_FOUR_BUFFER_CIRCUITS_GENERATED_INL",
        "",
        "#include <stdint.h>",
        "",
        "#define LEOPARD_FF8XOR_FOR_EACH_FOUR_TUPLE(M) \\",
    ]
    for index, tuple_value in enumerate(tuples):
        suffix = " \\" if index + 1 != len(tuples) else ""
        lines.append("    M(%d, %d, %d, %d)%s" % (
            index, tuple_value[0], tuple_value[1], tuple_value[2], suffix))
    lines.extend([
        "",
        "namespace leopard { namespace ff8xor { namespace generated4 {",
        "",
        "struct FourBufferSkewTuple",
        "{",
        "    uint8_t Skew01;",
        "    uint8_t Skew23;",
        "    uint8_t Skew02;",
        "};",
        "",
        "static const char kFourBufferCircuitChecksum[] = \"%s\";" % checksum,
        "static const FourBufferSkewTuple kFourBufferSkewTuples[%d] = {" %
        len(tuples),
    ])
    for tuple_value in tuples:
        lines.append("    { %d, %d, %d }," % tuple_value)
    lines.extend(["};", ""])

    for family_name, prefix in (("fft", "FFT4"), ("ifft", "IFFT4")):
        records = artifact["families"][family_name]
        append_uint_array(lines, "k%sXor2Counts" % prefix, [
            record["xor2_selected"]["xor2_metadata"]["instruction_count"]
            for record in records])
        append_uint_array(lines, "k%sXor3Xor2Counts" % prefix, [
            record["xor3_selected"]["xor3_metadata"]["xor2_count"]
            for record in records])
        append_uint_array(lines, "k%sXor3Counts" % prefix, [
            record["xor3_selected"]["xor3_metadata"]["xor3_count"]
            for record in records])
        append_uint_array(lines, "k%sXor3InstructionCounts" % prefix, [
            record["xor3_selected"]["xor3_metadata"]["instruction_count"]
            for record in records])
        append_uint_array(lines, "k%sXor2Depths" % prefix, [
            record["xor2_selected"]["xor2_metadata"]["depth"]
            for record in records])
        append_uint_array(lines, "k%sXor3Depths" % prefix, [
            record["xor3_selected"]["xor3_metadata"]["depth"]
            for record in records])
        lines.append("")

    lines.extend([
        "template <unsigned TupleIndex> struct FFT4CircuitXor2;",
        "template <unsigned TupleIndex> struct FFT4CircuitXor3;",
        "template <unsigned TupleIndex> struct IFFT4CircuitXor2;",
        "template <unsigned TupleIndex> struct IFFT4CircuitXor3;",
        "",
    ])
    for family_name, prefix in (("fft", "FFT4"), ("ifft", "IFFT4")):
        for record in artifact["families"][family_name]:
            append_circuit(
                lines, prefix + "CircuitXor2", record["tuple_index"],
                record["xor2_selected"]["xor2_operations"])
            append_circuit(
                lines, prefix + "CircuitXor3", record["tuple_index"],
                record["xor3_selected"]["xor3_operations"])

    lines.extend([
        "}}} // namespace leopard::ff8xor::generated4",
        "",
        "#endif // LEOPARD_FF8_XOR_FOUR_BUFFER_CIRCUITS_GENERATED_INL",
        "",
    ])
    return "\n".join(lines)


def summarize_numbers(values):
    return {
        "total": sum(values),
        "min": min(values),
        "max": max(values),
        "average": float(sum(values)) / len(values),
    }


def public_candidate(candidate):
    return {
        "name": candidate["name"],
        "gate_count": candidate["gate_count"],
        "cnot_depth": candidate["cnot_depth"],
        "xor2": candidate["xor2_metadata"],
        "xor3_variant": candidate["xor3_variant"],
        "xor3": candidate["xor3_metadata"],
    }


def generate_metadata(artifact, generated_source_bytes):
    checksum = checksum_payload(artifact)
    result = {
        "schema": "leopard-ff8xor-four-buffer-circuits-v1",
        "generator": artifact["generator"],
        "checksum_sha256": checksum,
        "field_polynomial": "0x11D",
        "cantor_basis": list(two_buffer.CANTOR_BASIS),
        "wire_count": WIRE_COUNT,
        "peak_live_wires": WIRE_COUNT,
        "reachable_tuple_count": len(artifact["reachability"]["tuples"]),
        "reachable_skew_count": artifact["skew_count"],
        "valid_parameter_pairs_enumerated":
            artifact["reachability"]["parameter_count"],
        "random_seed": artifact["random_seed"],
        "random_states_per_tuple_and_direction": artifact["random_state_count"],
        "generated_source_bytes": generated_source_bytes,
        "production_integration_note": (
            "The 32-live-ZMM kernel requires destructive +v inline-assembly "
            "XorValue/Xor3Value wrappers; ordinary intrinsic wrappers were "
            "observed to spill after compiler DAG reassociation.  Re-run the "
            "strict GCC/Clang inspector for the production translation unit."),
        "tuple_frequency_note": (
            "Counts are fused calls in the exhaustive valid-k/r API-shape "
            "census, not a production workload model."),
        "families": {},
    }

    for family_name, records in artifact["families"].items():
        public_records = []
        for record in records:
            public_records.append({
                "tuple_index": record["tuple_index"],
                "tuple": list(record["tuple"]),
                "parameter_census_calls": artifact["reachability"][
                    family_name + "_counts"][record["tuple"]],
                "component_variants": list(record["component_variants"]),
                "candidates": [public_candidate(candidate)
                               for candidate in record["candidates"]],
                "xor2_selected": public_candidate(record["xor2_selected"]),
                "xor3_selected": public_candidate(record["xor3_selected"]),
            })
        xor2_counts = [record["xor2_selected"]["xor2_metadata"][
            "instruction_count"] for record in records]
        xor3_counts = [record["xor3_selected"]["xor3_metadata"][
            "instruction_count"] for record in records]
        xor3_ternaries = [record["xor3_selected"]["xor3_metadata"][
            "xor3_count"] for record in records]
        result["families"][family_name] = {
            "xor2_instruction_summary": summarize_numbers(xor2_counts),
            "xor3_instruction_summary": summarize_numbers(xor3_counts),
            "xor3_ternary_summary": summarize_numbers(xor3_ternaries),
            "records": public_records,
        }
    return result


def check_file(path, expected):
    try:
        current = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        print("generated file is missing: %s" % path, file=sys.stderr)
        return False
    if current == expected:
        print("FF8 four-buffer generated file is current: %s" % path)
        return True
    print("generated file is stale: %s" % path, file=sys.stderr)
    for line in difflib.unified_diff(
            current.splitlines(), expected.splitlines(),
            fromfile=str(path), tofile="regenerated", lineterm=""):
        print(line, file=sys.stderr)
    return False


def print_summary(artifact, elapsed_seconds):
    print("Reachable tuples: %d across %d valid k/r pairs and %d skews" % (
        len(artifact["reachability"]["tuples"]),
        artifact["reachability"]["parameter_count"], artifact["skew_count"]))
    print("Circuit checksum: %s" % checksum_payload(artifact))
    for family_name in ("fft", "ifft"):
        records = artifact["families"][family_name]
        xor2_values = [record["xor2_selected"]["xor2_metadata"][
            "instruction_count"] for record in records]
        xor3_values = [record["xor3_selected"]["xor3_metadata"][
            "instruction_count"] for record in records]
        ternaries = [record["xor3_selected"]["xor3_metadata"][
            "xor3_count"] for record in records]
        print("%s: XOR2=%d (min=%d max=%d avg=%.6f); "
              "mixed=%d (min=%d max=%d avg=%.6f; XOR3=%d)" % (
                  family_name.upper(), sum(xor2_values), min(xor2_values),
                  max(xor2_values), float(sum(xor2_values)) / len(xor2_values),
                  sum(xor3_values), min(xor3_values), max(xor3_values),
                  float(sum(xor3_values)) / len(xor3_values), sum(ternaries)))
    print("Generation and proof time: %.3f seconds" % elapsed_seconds)


def parse_arguments():
    repository_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--output", type=Path,
        default=repository_root / "generated" /
        "LeopardFF8XorFourBufferCircuits.inl")
    parser.add_argument(
        "--metadata", type=Path,
        default=repository_root / "generated" /
        "FF8XorFourBufferCircuits.json")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    begin = time.monotonic()
    artifact = build_artifact()
    generated_cpp = generate_cpp(artifact)
    metadata = generate_metadata(
        artifact, len(generated_cpp.encode("utf-8")))
    generated_metadata = json.dumps(
        metadata, sort_keys=True, indent=2) + "\n"
    elapsed = time.monotonic() - begin

    if arguments.check:
        success = check_file(arguments.output, generated_cpp)
        success = check_file(arguments.metadata, generated_metadata) and success
        print_summary(artifact, elapsed)
        return 0 if success else 1

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.metadata.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(generated_cpp, encoding="utf-8")
    arguments.metadata.write_text(generated_metadata, encoding="utf-8")
    print("Generated %s" % arguments.output)
    print("Generated %s" % arguments.metadata)
    print_summary(artifact, elapsed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
