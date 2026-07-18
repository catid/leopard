#!/usr/bin/env python3
"""Independent finite tests for generated FF8 four-buffer circuits."""

from __future__ import print_function

import hashlib
import importlib.util
import json
import random
import re
import sys
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

SPEC = importlib.util.spec_from_file_location(
    "ff8_four_buffer_generator",
    str(TOOLS / "generate_ff8_xor_four_buffer_circuits.py"))
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)
INSPECTOR_SPEC = importlib.util.spec_from_file_location(
    "ff8_four_buffer_inspector",
    str(TOOLS / "inspect_ff8xor_four_buffer_avx512.py"))
INSPECTOR = importlib.util.module_from_spec(INSPECTOR_SPEC)
INSPECTOR_SPEC.loader.exec_module(INSPECTOR)


FIELD_POLYNOMIAL = 0x11D
CANTOR_BASIS = (1, 214, 152, 146, 86, 200, 88, 230)
TUPLE_CHECKSUM = "b8733d765f8ff4663d0615c898b01791f04b4c4536bd14c3a99fb08e37e6f67d"


def polynomial_multiply(left, right):
    product = 0
    for unused in range(8):
        if right & 1:
            product ^= left
        right >>= 1
        left <<= 1
        if left & 0x100:
            left ^= FIELD_POLYNOMIAL
    return product


def cantor_to_polynomial(value):
    result = 0
    for bit, basis in enumerate(CANTOR_BASIS):
        if (value >> bit) & 1:
            result ^= basis
    return result


POLYNOMIAL_TO_CANTOR = [0] * 256
for _coordinate in range(256):
    POLYNOMIAL_TO_CANTOR[cantor_to_polynomial(_coordinate)] = _coordinate
POLYNOMIAL_EXP = [1]
for _exponent in range(1, 256):
    POLYNOMIAL_EXP.append(polynomial_multiply(POLYNOMIAL_EXP[-1], 2))


def independent_multiply_log(value, logarithm):
    product = polynomial_multiply(
        cantor_to_polynomial(value), POLYNOMIAL_EXP[logarithm])
    return POLYNOMIAL_TO_CANTOR[product]


def independent_butterfly(left, right, skew, inverse):
    if inverse:
        right ^= left
        if skew != 255:
            left ^= independent_multiply_log(right, skew)
    else:
        if skew != 255:
            left ^= independent_multiply_log(right, skew)
        right ^= left
    return left, right


def independent_four_buffer(state, tuple_value, inverse):
    values = [(state >> (8 * index)) & 0xFF for index in range(4)]

    def apply(left, right, skew):
        values[left], values[right] = independent_butterfly(
            values[left], values[right], skew, inverse)

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
    return sum(value << (8 * index) for index, value in enumerate(values))


def local_apply_operations(state, operations):
    for operation in operations:
        opcode, destination = operation[:2]
        value = (state >> operation[2]) & 1
        if opcode == 3:
            value ^= (state >> operation[3]) & 1
        elif opcode != 2:
            raise AssertionError("unknown circuit opcode")
        state ^= value << destination
    return state


def local_matrix(operations):
    rows = [0] * 32
    for column in range(32):
        output = local_apply_operations(1 << column, operations)
        for row in range(32):
            rows[row] |= ((output >> row) & 1) << column
    return tuple(rows)


class FourBufferCircuitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact = GENERATOR.build_artifact()
        cls.generated = GENERATOR.generate_cpp(cls.artifact)
        cls.metadata = GENERATOR.generate_metadata(
            cls.artifact, len(cls.generated.encode("utf-8")))

    def test_reachability_census_is_frozen_and_complete(self):
        reachability = self.artifact["reachability"]
        tuples = reachability["tuples"]
        self.assertEqual(reachability["parameter_count"], 21845)
        self.assertEqual(len(tuples), 64)
        self.assertEqual(set(reachability["fft_counts"]), set(tuples))
        self.assertEqual(set(reachability["ifft_counts"]), set(tuples))
        self.assertEqual(sum(reachability["fft_counts"].values()), 4747804)
        self.assertEqual(sum(reachability["ifft_counts"].values()), 4945408)
        self.assertEqual(len({tuple_value[0] for tuple_value in tuples}), 64)
        self.assertEqual(
            [value for value in tuples if 255 in value], [(255, 85, 255)])
        encoded = json.dumps(tuples, separators=(",", ":")).encode("ascii")
        self.assertEqual(hashlib.sha256(encoded).hexdigest(), TUPLE_CHECKSUM)

    def test_truncated_partial_radix4_units_are_not_counted(self):
        for recorder in (
                GENERATOR.record_ifft_tuples,
                GENERATOR.record_fft_tuples):
            for count_truncated in (1, 2):
                counter = Counter()
                recorder(counter, 4, count_truncated, 0)
                self.assertEqual(counter, Counter())
            for count_truncated in (3, 4):
                counter = Counter()
                recorder(counter, 4, count_truncated, 0)
                self.assertEqual(sum(counter.values()), 1)

            counter = Counter()
            recorder(counter, 16, 3, 0)
            self.assertEqual(sum(counter.values()), 1)
            counter = Counter()
            recorder(counter, 16, 9, 0)
            self.assertEqual(sum(counter.values()), 6)

    def test_reachability_matches_independent_schedule_corpus(self):
        corpus = json.loads((
            ROOT / "generated" / "FF8XorScheduleCorpus.json").read_text(
                encoding="utf-8"))

        def parse_histogram(record, field):
            return {
                tuple(int(value) for value in key.split(",")): count
                for key, count in record.get(field, {}).items()
            }

        fft_weighted = {}
        ifft_weighted = {}
        combined_weighted = {}
        for record in corpus["records"]:
            fft = parse_histogram(record, "fft_four_buffer_tuple_frequency")
            ifft = parse_histogram(record, "ifft_four_buffer_tuple_frequency")
            combined = parse_histogram(record, "four_buffer_tuple_frequency")
            observed_combined = dict(fft)
            for tuple_value, count in ifft.items():
                observed_combined[tuple_value] = (
                    observed_combined.get(tuple_value, 0) + count)
            self.assertEqual(observed_combined, combined, record["id"])
            for source, destination in (
                    (fft, fft_weighted), (ifft, ifft_weighted),
                    (combined, combined_weighted)):
                for tuple_value, count in source.items():
                    destination[tuple_value] = (
                        destination.get(tuple_value, 0) + count)

        expected = set(self.artifact["reachability"]["tuples"])
        self.assertEqual(set(fft_weighted), expected)
        self.assertEqual(set(ifft_weighted), expected)
        self.assertEqual(set(combined_weighted), expected)
        self.assertEqual(sum(fft_weighted.values()), 2728)
        self.assertEqual(sum(ifft_weighted.values()), 9188)
        self.assertEqual(sum(combined_weighted.values()), 11916)
        for tuple_value in expected:
            self.assertEqual(
                combined_weighted[tuple_value],
                fft_weighted[tuple_value] + ifft_weighted[tuple_value])

    def test_all_maps_match_independent_scalar_and_are_mutual_inverses(self):
        deterministic_random = random.Random(0xC17F04B)
        random_states = tuple(
            deterministic_random.randrange(1 << 32) for unused in range(128))
        states = (
            (0, (1 << 32) - 1)
            + tuple(1 << bit for bit in range(32))
            + random_states)

        for tuple_index, tuple_value in enumerate(
                self.artifact["reachability"]["tuples"]):
            fft = self.artifact["families"]["fft"][tuple_index]
            ifft = self.artifact["families"]["ifft"][tuple_index]
            self.assertEqual(
                local_matrix(fft["xor2_selected"]["xor2_operations"]),
                fft["rows"])
            self.assertEqual(
                local_matrix(fft["xor3_selected"]["xor3_operations"]),
                fft["rows"])
            self.assertEqual(
                local_matrix(ifft["xor2_selected"]["xor2_operations"]),
                ifft["rows"])
            self.assertEqual(
                local_matrix(ifft["xor3_selected"]["xor3_operations"]),
                ifft["rows"])

            for state in states:
                expected_fft = independent_four_buffer(
                    state, tuple_value, False)
                expected_ifft = independent_four_buffer(
                    state, tuple_value, True)
                self.assertEqual(
                    local_apply_operations(
                        state, fft["xor2_selected"]["xor2_operations"]),
                    expected_fft)
                self.assertEqual(
                    local_apply_operations(
                        state, fft["xor3_selected"]["xor3_operations"]),
                    expected_fft)
                self.assertEqual(
                    local_apply_operations(
                        state, ifft["xor2_selected"]["xor2_operations"]),
                    expected_ifft)
                self.assertEqual(
                    local_apply_operations(
                        state, ifft["xor3_selected"]["xor3_operations"]),
                    expected_ifft)
                self.assertEqual(
                    independent_four_buffer(expected_fft, tuple_value, True),
                    state)
                self.assertEqual(
                    independent_four_buffer(expected_ifft, tuple_value, False),
                    state)

    def test_every_candidate_was_map_verified_and_metadata_is_exact(self):
        for family in ("fft", "ifft"):
            for record in self.artifact["families"][family]:
                for candidate in record["candidates"]:
                    self.assertEqual(
                        local_matrix(candidate["xor2_operations"]),
                        record["rows"])
                    expanded = GENERATOR.xor3.expand_operations(
                        candidate["xor3_operations"], 32)
                    observed_xor2 = sum(
                        operation[0] == 2
                        for operation in candidate["xor3_operations"])
                    observed_xor3 = sum(
                        operation[0] == 3
                        for operation in candidate["xor3_operations"])
                    self.assertEqual(
                        candidate["xor3_metadata"]["xor2_count"],
                        observed_xor2)
                    self.assertEqual(
                        candidate["xor3_metadata"]["xor3_count"],
                        observed_xor3)
                    self.assertEqual(
                        candidate["xor3_metadata"]["instruction_count"],
                        observed_xor2 + observed_xor3)
                    self.assertEqual(
                        candidate["xor3_metadata"]["peak_live_wires"], 32)
                    self.assertEqual(
                        local_matrix(candidate["xor3_operations"]),
                        record["rows"])
                    self.assertEqual(
                        GENERATOR.xor3.circuit_matrix(expanded, 32),
                        record["rows"])

    def test_candidate_comparison_retains_negative_whole_map_result(self):
        for family in ("fft", "ifft"):
            records = self.artifact["families"][family]
            self.assertEqual(len(records), 64)
            for record in records:
                names = {candidate["name"]
                         for candidate in record["candidates"]}
                self.assertIn("direct", names)
                self.assertIn("composed-two-buffer-portfolio", names)
                self.assertTrue(any(name.startswith("whole-map/")
                                    for name in names))
                self.assertFalse(
                    record["xor2_selected"]["name"].startswith("whole-map/"))
                self.assertFalse(
                    record["xor3_selected"]["name"].startswith("whole-map/"))

    def test_generated_hot_bodies_use_only_named_literal_operations(self):
        self.assertEqual(self.generated.count("template <> struct FFT4"), 128)
        self.assertEqual(self.generated.count("template <> struct IFFT4"), 128)
        self.assertNotIn("registers[", self.generated)
        self.assertNotIn("wire[", self.generated)
        self.assertNotIn("for (", self.generated)
        self.assertNotIn("while (", self.generated)
        bodies = re.findall(r"static LEO_FORCE_INLINE void Apply\(.*?\n    \}\n\};",
                            self.generated, re.DOTALL)
        self.assertEqual(len(bodies), 256)
        xor2_pattern = re.compile(
            r"^([abcd][0-7]) = XorValue\(([abcd][0-7]), "
            r"([abcd][0-7])\);$")
        xor3_pattern = re.compile(
            r"^([abcd][0-7]) = Xor3Value\(([abcd][0-7]), "
            r"([abcd][0-7]), ([abcd][0-7])\);$")

        def parse_literal(statement):
            xor2_match = xor2_pattern.fullmatch(statement)
            xor3_match = xor3_pattern.fullmatch(statement)
            if (xor2_match is None) == (xor3_match is None):
                raise ValueError("not exactly one literal XOR operation")
            match = xor2_match if xor2_match is not None else xor3_match
            names = match.groups()
            if names[0] != names[1]:
                raise ValueError("the assignment destination must be updated")
            operation_names = (names[0],) + names[2:]
            if len(set(operation_names)) != len(operation_names):
                raise ValueError("literal operation wires must be distinct")
            opcode = 2 if xor2_match is not None else 3
            return (opcode,) + tuple(
                GENERATOR.WIRE_NAMES.index(name) for name in operation_names)

        expected_operations = []
        for family in ("fft", "ifft"):
            for record in self.artifact["families"][family]:
                expected_operations.append(
                    record["xor2_selected"]["xor2_operations"])
                expected_operations.append(
                    record["xor3_selected"]["xor3_operations"])

        for body, expected in zip(bodies, expected_operations):
            statements = [line.strip() for line in body.splitlines()
                          if " = Xor" in line]
            self.assertGreater(len(statements), 0)
            self.assertEqual(
                tuple(parse_literal(statement) for statement in statements),
                expected)
        for malformed in (
                "registers[0] = XorValue(registers[0], registers[1]);",
                "a8 = XorValue(a8, b0);",
                "a0 = XorValue(a0, wires[index]);",
                "a0 ^= b0;",
                "a0 = Xor3Value(a0, b0);",
                "a0 = XorValue(a0, b0, c0);",
                "a0 = XorValue(b0, c0);",
                "a0 = XorValue(a0, a0);",
                "a0 = Xor3Value(a0, b0, b0);"):
            with self.assertRaises(ValueError, msg=malformed):
                parse_literal(malformed)

    def test_checked_in_outputs_are_deterministic_and_current(self):
        output_path = ROOT / "generated" / "LeopardFF8XorFourBufferCircuits.inl"
        metadata_path = ROOT / "generated" / "FF8XorFourBufferCircuits.json"
        self.assertEqual(
            output_path.read_text(encoding="utf-8"), self.generated)
        expected_metadata = json.dumps(
            self.metadata, sort_keys=True, indent=2) + "\n"
        self.assertEqual(
            metadata_path.read_text(encoding="utf-8"), expected_metadata)
        self.assertEqual(
            self.metadata["checksum_sha256"],
            GENERATOR.checksum_payload(self.artifact))
        self.assertEqual(
            self.metadata["generated_source_bytes"],
            len(self.generated.encode("utf-8")))

    def test_strict_inspector_rejects_stress_and_representative_spills(self):
        for selection in (
                "frequency-representative", "max-instruction-stress"):
            kernel = {
                "symbol": "fixture",
                "selection": selection,
                "lowering": "xor2",
                "call_count": 0,
                "vector_stack_reference_count": 1,
                "forbidden_field_instruction_count": 0,
                "vpternlog_count": 0,
                "vector_xor_count": 1,
            }
            failures = INSPECTOR.validate_strict(({
                "compiler": "/fixture/g++", "kernels": (kernel,),
            },))
            self.assertEqual(
                failures, ["g++:fixture contains vector stack spills"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
