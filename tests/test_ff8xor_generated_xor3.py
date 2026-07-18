#!/usr/bin/env python3
"""Independently verify every emitted FF8 XOR3 named-register schedule."""

from __future__ import print_function

import hashlib
import json
import math
import random
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "generated" / "LeopardFF8XorCircuits.inl"
XOR3_PATH = ROOT / "generated" / "LeopardFF8XorCircuitsXor3.inl"
EVIDENCE_PATH = (ROOT / "tools" / "profiles" /
                 "FF8XorGeneratedXor3EvaluationZen5Gcc13.json")

FIELD_BITS = 8
FIELD_ORDER = 1 << FIELD_BITS
FIELD_MODULUS = FIELD_ORDER - 1
FIELD_POLYNOMIAL = 0x11D
CANTOR_BASIS = (1, 214, 152, 146, 86, 200, 88, 230)

XOR2 = 2
XOR3 = 3


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_checksum(value):
    checked = dict(value)
    checked.pop("checksum_sha256", None)
    encoded = json.dumps(
        checked, sort_keys=True, separators=(",", ":"),
        ensure_ascii=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def validate_evidence(base_text, xor3_text):
    evidence = json.loads(EVIDENCE_PATH.read_text(encoding="utf-8"))
    require(evidence.get("schema") ==
            "leopard.ff8xor.generated-xor3-evaluation.v1",
            "generated XOR3 evidence has the wrong schema")
    require(evidence.get("checksum_sha256") == canonical_checksum(evidence),
            "generated XOR3 evidence checksum is stale")

    artifacts = evidence["artifacts"]
    require(artifacts["xor2"]["path"] ==
            "generated/LeopardFF8XorCircuits.inl" and
            artifacts["xor3"]["path"] ==
            "generated/LeopardFF8XorCircuitsXor3.inl",
            "evidence names the wrong generated artifacts")
    require(artifacts["xor2"]["file_sha256"] == file_sha256(BASE_PATH) and
            artifacts["xor3"]["file_sha256"] == file_sha256(XOR3_PATH),
            "generated artifact changed without refreshing XOR3 evidence")
    base_checksum = re.search(
        r"Circuit checksum \(SHA-256\): ([0-9a-f]{64})", base_text)
    xor3_checksum = re.search(
        r"Explicit XOR3 schedule checksum \(SHA-256\): ([0-9a-f]{64})",
        xor3_text)
    require(base_checksum and xor3_checksum and
            artifacts["xor2"]["circuit_checksum_sha256"] ==
            base_checksum.group(1) and
            artifacts["xor3"]["base_circuit_checksum_sha256"] ==
            base_checksum.group(1) and
            artifacts["xor3"]["schedule_checksum_sha256"] ==
            xor3_checksum.group(1),
            "evidence checksum provenance does not match generated source")

    environment = evidence["environment"]
    require(environment["cpu_model"] ==
            "AMD Ryzen Threadripper PRO 9985WX 64-Cores" and
            environment["pinned_logical_cpu"] == 24 and
            environment["compiler"] == "GCC 13.3.0" and
            environment["build_type"] == "Release (NDEBUG)" and
            environment["source_base_commit"] ==
            "fdd1ebd2a7fe992890723c9673adc200e838efed" and
            environment["configure_command"] ==
            "cmake -S . -B build -DCMAKE_BUILD_TYPE=Release "
            "-DENABLE_OPENMP=OFF" and
            environment["build_command"] ==
            "cmake --build build -j 30 --target bench_leopard_ff8xor" and
            environment["circuit_object_flags"] == [
                "-march=x86-64", "-mavx2", "-mavx512f", "-mavx512vl",
                "-fno-tree-reassoc"],
            "generated XOR3 build provenance changed unexpectedly")

    methodology = evidence["methodology"]
    require(methodology["benchmark_command"] ==
            "taskset -c 24 $BINARY --quick --json --abba "
            "--ff8xor-mode $MODE --cpu 24" and
            methodology["modes"] == ["avx512vl", "avx512zmm"] and
            methodology["process_repetitions_per_variant_per_mode"] == 9 and
            methodology["inner_measurement_order"] == "ABBA" and
            methodology["warmup_iterations"] == 1 and
            methodology["measured_iterations"] == 3 and
            methodology["minimum_sample_usec"] == 250 and
            methodology["end_to_end_row_identities_per_process"] == 14 and
            methodology["end_to_end_rows_per_variant_per_mode"] == 126 and
            methodology["normalized_circuit_row_identities_per_process"] ==
            9 and
            methodology["normalized_circuit_rows_per_variant_per_mode"] ==
            81 and
            methodology["transpose_included"] is False and
            methodology["allocation_included"] is False,
            "generated XOR3 benchmark methodology is incomplete")
    require(methodology["raw_samples_retained"] is False and
            methodology["raw_sample_manifest_definition"] ==
            "SHA-256 of the concatenated binary SHA-256 digests of nine "
            "JSONL files sorted by filename within one variant and mode." and
            "unavailable" in methodology["raw_samples_note"] and
            "do not make the raw data recoverable" in
            methodology["raw_samples_note"],
            "evidence must not imply unavailable raw samples are reparsable")

    require(evidence["generated_schedule_totals"] == {
        "multiply": {
            "source_xor2_gates": 4903,
            "xor2_operations": 3113,
            "xor3_operations": 895,
            "scheduled_instructions": 4008,
        },
        "fft": {
            "source_xor2_gates": 10240,
            "xor2_operations": 3072,
            "xor3_operations": 3584,
            "scheduled_instructions": 6656,
        },
        "ifft": {
            "source_xor2_gates": 10240,
            "xor2_operations": 3072,
            "xor3_operations": 3584,
            "scheduled_instructions": 6656,
        },
    }, "generated schedule totals changed unexpectedly")

    exact_speedups = {
        ("full_explicit", "avx512vl"): 0.9812084393895245,
        ("full_explicit", "avx512zmm"): 0.9768714400996034,
        ("fft_only", "avx512vl"): 1.0036487693501197,
        ("fft_only", "avx512zmm"): 0.9928193776180321,
    }
    exact_family_speedups = {
        ("full_explicit", "avx512vl"): {
            "fft": 1.004570038756457,
            "ifft": 0.9956878052426169,
            "multiply": 0.9987391169402606,
        },
        ("full_explicit", "avx512zmm"): {
            "fft": 0.9888663803348539,
            "ifft": 0.9961432213563884,
            "multiply": 0.9910613430083693,
        },
        ("fft_only", "avx512vl"): {
            "fft": 0.9940766941655861,
            "ifft": 1.0101652961560683,
            "multiply": 0.998067289413305,
        },
        ("fft_only", "avx512zmm"): {
            "fft": 0.9921804810779209,
            "ifft": 1.0028381416581607,
            "multiply": 0.9969145382482083,
        },
    }
    expected_groups = {
        ("decode", 1024): 5,
        ("decode", 65536): 5,
        ("encode", 1024): 2,
        ("encode", 65536): 2,
    }
    for variant in ("full_explicit", "fft_only"):
        for mode in ("avx512vl", "avx512zmm"):
            item = evidence["results"][variant][mode]
            require(re.fullmatch(
                r"[0-9a-f]{64}", item["control_manifest_sha256"]) and
                re.fullmatch(
                    r"[0-9a-f]{64}", item["candidate_manifest_sha256"]),
                    "raw-sample manifest digest is malformed")
            end_to_end = item["normalized_end_to_end"]
            require(end_to_end["geomean_speedup"] ==
                    exact_speedups[(variant, mode)] and
                    end_to_end["minimum_speedup"] <=
                    end_to_end["geomean_speedup"] <=
                    end_to_end["maximum_speedup"],
                    "normalized end-to-end result changed unexpectedly")
            groups = {
                (entry["operation"], entry["buffer_bytes"]):
                entry["row_identities"]
                for entry in end_to_end["by_operation_and_buffer"]
            }
            require(groups == expected_groups and
                    sum(groups.values()) == 14,
                    "end-to-end benchmark result set is incomplete")
            families = {
                entry["family"]: entry["geomean_speedup"]
                for entry in item["normalized_circuit_families"]
            }
            require(families == exact_family_speedups[(variant, mode)] and
                    all(entry["coefficients"] == 3 and
                        math.isfinite(entry["minimum_speedup"]) and
                        entry["minimum_speedup"] <=
                        entry["geomean_speedup"] <=
                        entry["maximum_speedup"]
                        for entry in item["normalized_circuit_families"]),
                    "circuit-family benchmark result changed unexpectedly")

    expected_control = {
        "multiply_avx512vl": {
            "instructions": 18092, "code_bytes": 72644,
            "xor2": 2265, "xor3": 924},
        "multiply_avx512zmm": {
            "instructions": 18182, "code_bytes": 83268,
            "xor2": 2265, "xor3": 924},
        "fft_avx512vl": {
            "instructions": 27254, "code_bytes": 123758,
            "xor2": 3059, "xor3": 3438},
        "fft_avx512zmm": {
            "instructions": 27247, "code_bytes": 143558,
            "xor2": 3059, "xor3": 3438},
        "ifft_avx512vl": {
            "instructions": 27021, "code_bytes": 121750,
            "xor2": 3076, "xor3": 3308},
        "ifft_avx512zmm": {
            "instructions": 25342, "code_bytes": 132214,
            "xor2": 3076, "xor3": 3308},
    }
    expected_full = {
        "multiply_avx512vl": {
            "instructions": 19258, "code_bytes": 77340,
            "xor2": 2297, "xor3": 1147},
        "multiply_avx512zmm": {
            "instructions": 19397, "code_bytes": 87964,
            "xor2": 2297, "xor3": 1147},
        "fft_avx512vl": {
            "instructions": 26809, "code_bytes": 119574,
            "xor2": 2913, "xor3": 3584},
        "fft_avx512zmm": {
            "instructions": 26810, "code_bytes": 135974,
            "xor2": 2913, "xor3": 3584},
        "ifft_avx512vl": {
            "instructions": 26979, "code_bytes": 121950,
            "xor2": 3056, "xor3": 3584},
        "ifft_avx512zmm": {
            "instructions": 26973, "code_bytes": 138998,
            "xor2": 3056, "xor3": 3584},
    }
    assembly = evidence["assembly"]
    require(assembly["command"] ==
            "python3 tools/inspect_ff8xor_assembly.py "
            "build/liblibleopard.a --strict --fail-on-spills "
            "--require-avx2 --require-avx512 --json" and
            assembly["raw_reports_retained"] is False and all(
                re.fullmatch(r"[0-9a-f]{64}", assembly[field])
                for field in (
                    "control_report_sha256",
                    "full_explicit_report_sha256",
                    "fft_only_report_sha256")),
            "assembly command or unavailable-report hashes are incomplete")
    require(assembly["control"] == expected_control and
            assembly["full_explicit"] == expected_full,
            "exact production assembly totals changed unexpectedly")
    expected_fft_only = dict(expected_control)
    expected_fft_only["fft_avx512vl"] = expected_full["fft_avx512vl"]
    expected_fft_only["fft_avx512zmm"] = expected_full["fft_avx512zmm"]
    require(assembly["fft_only"] == expected_fft_only and
            assembly["all_variants"] == {
                "strict_pass": True,
                "spill_check_pass": True,
                "hard_violation_count": 0,
                "spill_warning_count": 0,
                "inner_loop_calls": 0,
                "forbidden_instructions": 0,
                "non_xor3_ternary_immediates": 0,
            },
            "FFT-only assembly or structural checks changed unexpectedly")

    decision = evidence["decision"]
    require(decision["enable_full_explicit_in_production"] is False and
            decision["enable_fft_only_in_production"] is False and
            decision["production_avx512_uses_historical_xor2_source"] is True,
            "negative XOR3 result was accidentally promoted")


def coordinate_to_polynomial(value):
    result = 0
    for bit, basis in enumerate(CANTOR_BASIS):
        if (value >> bit) & 1:
            result ^= basis
    return result


COORDINATE_TO_POLYNOMIAL = tuple(
    coordinate_to_polynomial(value) for value in range(FIELD_ORDER))
POLYNOMIAL_TO_COORDINATE_LIST = [0] * FIELD_ORDER
for _coordinate, _polynomial in enumerate(COORDINATE_TO_POLYNOMIAL):
    POLYNOMIAL_TO_COORDINATE_LIST[_polynomial] = _coordinate
POLYNOMIAL_TO_COORDINATE = tuple(POLYNOMIAL_TO_COORDINATE_LIST)
del POLYNOMIAL_TO_COORDINATE_LIST, _coordinate, _polynomial


def polynomial_multiply(left, right):
    result = 0
    for unused_bit in range(FIELD_BITS):
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & FIELD_ORDER:
            left ^= FIELD_POLYNOMIAL
    return result


POLYNOMIAL_EXP_LIST = [1] * FIELD_ORDER
for _logarithm in range(1, FIELD_ORDER):
    POLYNOMIAL_EXP_LIST[_logarithm] = polynomial_multiply(
        POLYNOMIAL_EXP_LIST[_logarithm - 1], 2)
POLYNOMIAL_EXP = tuple(POLYNOMIAL_EXP_LIST)
del POLYNOMIAL_EXP_LIST, _logarithm


def multiply_log(value, logarithm):
    product = polynomial_multiply(
        COORDINATE_TO_POLYNOMIAL[value], POLYNOMIAL_EXP[logarithm])
    return POLYNOMIAL_TO_COORDINATE[product]


def forward_butterfly(state, skew):
    x = state & 0xff
    y = state >> 8
    if skew != FIELD_MODULUS:
        x ^= multiply_log(y, skew)
    y ^= x
    return x | (y << 8)


def inverse_butterfly(state, skew):
    x = state & 0xff
    y = state >> 8
    y ^= x
    if skew != FIELD_MODULUS:
        x ^= multiply_log(y, skew)
    return x | (y << 8)


def wire_names(width):
    names = tuple("x%d" % index for index in range(FIELD_BITS))
    if width == 16:
        names += tuple("y%d" % index for index in range(FIELD_BITS))
    return names


def parse_schedules(text, struct_name, width):
    pattern = re.compile(
        r"template <> struct %s<(\d+)>\n"
        r"\{.*?static LEO_FORCE_INLINE void Apply\(.*?\n"
        r"    \{(.*?)\n    \}\n\};" % re.escape(struct_name),
        re.DOTALL)
    names = wire_names(width)
    name_to_index = {name: index for index, name in enumerate(names)}
    schedules = {}
    xor2_pattern = re.compile(
        r"([xy]\d) = XorValue\(([xy]\d), ([xy]\d)\);")
    xor3_pattern = re.compile(
        r"([xy]\d) = Xor3Value\(([xy]\d), ([xy]\d), ([xy]\d)\);")
    void_pattern = re.compile(r"\(void\)([xy]\d);")

    for coefficient_text, body in pattern.findall(text):
        coefficient = int(coefficient_text)
        require(coefficient not in schedules,
                "%s coefficient %d is duplicated" %
                (struct_name, coefficient))
        operations = []
        for line in body.splitlines():
            statement = line.strip()
            if not statement:
                continue
            void_match = void_pattern.fullmatch(statement)
            if void_match:
                require(void_match.group(1) in name_to_index,
                        "invalid named wire in void cast")
                continue
            xor2_match = xor2_pattern.fullmatch(statement)
            xor3_match = xor3_pattern.fullmatch(statement)
            require(xor2_match is not None or xor3_match is not None,
                    "non-literal operation in %s<%d>: %s" %
                    (struct_name, coefficient, statement))
            match = xor2_match or xor3_match
            require(match.group(1) == match.group(2),
                    "destination is not the first XOR input")
            require(all(name in name_to_index for name in match.groups()),
                    "operation uses an invalid named wire")
            destination = name_to_index[match.group(1)]
            first_source = name_to_index[match.group(3)]
            require(destination != first_source,
                    "XOR source equals its destination")
            if xor2_match:
                operations.append((XOR2, destination, first_source))
            else:
                second_source = name_to_index[match.group(4)]
                require(destination != second_source and
                        first_source != second_source,
                        "XOR3 must use three distinct named wires")
                operations.append((
                    XOR3, destination, first_source, second_source))
        schedules[coefficient] = tuple(operations)

    require(set(schedules) == set(range(FIELD_ORDER)),
            "%s does not contain exactly 256 coefficients" % struct_name)
    return tuple(schedules[index] for index in range(FIELD_ORDER))


def apply_schedule(state, operations):
    for operation in operations:
        value = (state >> operation[2]) & 1
        if operation[0] == XOR3:
            value ^= (state >> operation[3]) & 1
        state ^= value << operation[1]
    return state


def schedule_rows(operations, width):
    rows = [0] * width
    for input_wire in range(width):
        output = apply_schedule(1 << input_wire, operations)
        for output_wire in range(width):
            rows[output_wire] |= (
                ((output >> output_wire) & 1) << input_wire)
    return tuple(rows)


def schedule_metadata(operations, width):
    xor2_count = sum(operation[0] == XOR2 for operation in operations)
    xor3_count = sum(operation[0] == XOR3 for operation in operations)
    last_layer = [0] * width
    depth = 0
    for operation in operations:
        wires = operation[1:]
        layer = max(last_layer[wire] for wire in wires) + 1
        for wire in wires:
            last_layer[wire] = layer
        depth = max(depth, layer)
    return {
        "xor2": xor2_count,
        "xor3": xor3_count,
        "instructions": len(operations),
        "depth": depth,
        "live": width,
        "code": xor2_count * 6 + xor3_count * 7,
    }


def parse_integer_array(text, name):
    match = re.search(
        r"static const uint(?:8|16)_t %s\[256\] = \{(.*?)\n\};" %
        re.escape(name), text, re.DOTALL)
    require(match is not None, "missing generated metadata array %s" % name)
    values = tuple(int(value) for value in re.findall(r"\d+", match.group(1)))
    require(len(values) == FIELD_ORDER,
            "%s does not have exactly 256 values" % name)
    return values


def verify_metadata(text, prefix, schedules, width):
    expected = tuple(schedule_metadata(operations, width)
                     for operations in schedules)
    fields = (
        ("Xor2Counts", "xor2"),
        ("Xor3Counts", "xor3"),
        ("InstructionCounts", "instructions"),
        ("Depths", "depth"),
        ("PeakLiveWireCounts", "live"),
        ("EstimatedCodeBytes", "code"),
    )
    for suffix, key in fields:
        observed = parse_integer_array(text, "k%sXor3%s" % (prefix, suffix))
        require(observed == tuple(entry[key] for entry in expected),
                "%s %s metadata does not match emitted bodies" %
                (prefix, key))
    variants = parse_integer_array(
        text, "k%sXor3ScheduleVariantIds" % prefix)
    require(all(0 <= variant < 4 for variant in variants),
            "%s has an invalid schedule variant id" % prefix)


def verify_multiplier_schedules(schedules):
    for logarithm, operations in enumerate(schedules):
        rows = schedule_rows(operations, 8)
        for value in range(FIELD_ORDER):
            expected = multiply_log(value, logarithm)
            require(apply_schedule(value, operations) == expected,
                    "multiplier %d failed input %d" % (logarithm, value))
        require(all(0 <= row < FIELD_ORDER for row in rows),
                "multiplier row is outside its width")


def verify_butterfly_schedules(forward, inverse):
    deterministic_random = random.Random(0xFF8A5123)
    random_states = tuple(
        deterministic_random.randrange(1 << 16) for unused in range(128))
    states = ((0, 0xffff) + tuple(1 << bit for bit in range(16)) +
              random_states)
    for skew in range(FIELD_ORDER):
        forward_rows = schedule_rows(forward[skew], 16)
        inverse_rows = schedule_rows(inverse[skew], 16)
        for input_wire in range(16):
            state = 1 << input_wire
            require(apply_schedule(state, forward[skew]) ==
                    forward_butterfly(state, skew),
                    "forward map mismatch for skew %d" % skew)
            require(apply_schedule(state, inverse[skew]) ==
                    inverse_butterfly(state, skew),
                    "inverse map mismatch for skew %d" % skew)
        require(len(forward_rows) == len(inverse_rows) == 16,
                "butterfly row width mismatch")
        for state in states:
            expected_forward = forward_butterfly(state, skew)
            expected_inverse = inverse_butterfly(state, skew)
            observed_forward = apply_schedule(state, forward[skew])
            observed_inverse = apply_schedule(state, inverse[skew])
            require(observed_forward == expected_forward,
                    "forward schedule mismatch for skew %d" % skew)
            require(observed_inverse == expected_inverse,
                    "inverse schedule mismatch for skew %d" % skew)
            require(apply_schedule(observed_forward, inverse[skew]) == state,
                    "inverse(forward) mismatch for skew %d" % skew)
            require(apply_schedule(observed_inverse, forward[skew]) == state,
                    "forward(inverse) mismatch for skew %d" % skew)


def main():
    base_text = BASE_PATH.read_text(encoding="utf-8")
    text = XOR3_PATH.read_text(encoding="utf-8")
    validate_evidence(base_text, text)
    require("registers[" not in text and "for each gate" not in text,
            "generated XOR3 hot code contains an indexed interpreter")

    base_checksum = re.search(
        r"Circuit checksum \(SHA-256\): ([0-9a-f]{64})", base_text)
    xor3_base_checksum = re.search(
        r"Base XOR2 circuit checksum \(SHA-256\): ([0-9a-f]{64})", text)
    xor3_checksum = re.search(
        r"Explicit XOR3 schedule checksum \(SHA-256\): ([0-9a-f]{64})",
        text)
    require(base_checksum and xor3_base_checksum and xor3_checksum,
            "generated checksum metadata is missing")
    require(base_checksum.group(1) == xor3_base_checksum.group(1),
            "XOR3 schedules were generated from a different base circuit")
    require(re.search(
        r'kXor3CircuitChecksum\[\] = "%s";' % xor3_checksum.group(1), text),
        "XOR3 checksum constant does not match the file header")

    multiply = parse_schedules(text, "MultiplyCircuitXor3", 8)
    forward = parse_schedules(text, "FFTCircuitXor3", 16)
    inverse = parse_schedules(text, "IFFTCircuitXor3", 16)
    verify_metadata(text, "Multiply", multiply, 8)
    verify_metadata(text, "FFT", forward, 16)
    verify_metadata(text, "IFFT", inverse, 16)
    verify_multiplier_schedules(multiply)
    verify_butterfly_schedules(forward, inverse)

    require(not multiply[FIELD_MODULUS],
            "multiply log 255 must remain the identity")
    for label, sentinel in (("FFT", forward[FIELD_MODULUS]),
                            ("IFFT", inverse[FIELD_MODULUS])):
        metadata = schedule_metadata(sentinel, 16)
        require(metadata["xor2"] == 8 and metadata["xor3"] == 0,
                "%s skew 255 must contain only eight sentinel XORs" % label)

    print("Generated FF8 XOR3 circuits independently verified: "
          "checksum=%s multiply=%d/%d FFT=%d/%d IFFT=%d/%d" % (
              xor3_checksum.group(1),
              sum(schedule_metadata(item, 8)["xor2"] for item in multiply),
              sum(schedule_metadata(item, 8)["xor3"] for item in multiply),
              sum(schedule_metadata(item, 16)["xor2"] for item in forward),
              sum(schedule_metadata(item, 16)["xor3"] for item in forward),
              sum(schedule_metadata(item, 16)["xor2"] for item in inverse),
              sum(schedule_metadata(item, 16)["xor3"] for item in inverse)))


if __name__ == "__main__":
    main()
