#!/usr/bin/env python3
"""Independent semantic tests for the FF8 AVX-512 XOR3 scheduler."""

from __future__ import print_function

import os
import random
import sys
import unittest


REPOSITORY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPOSITORY_ROOT, "tools"))

import ff8_xor3_schedule as xor3  # noqa: E402


class Xor3ScheduleTests(unittest.TestCase):
    def assert_schedule_equivalent(self, gates, operations, width):
        self.assertTrue(xor3.verify_schedule(gates, operations, width))
        if width <= 6:
            for state in range(1 << width):
                self.assertEqual(
                    xor3.apply_gates(state, gates, width),
                    xor3.apply_operations(state, operations, width))

    def test_validation_rejects_illegal_wires_and_operations(self):
        with self.assertRaises(ValueError):
            xor3.validate_gates(((0, 0),), 2)
        with self.assertRaises(ValueError):
            xor3.validate_gates(((2, 0),), 2)
        with self.assertRaises(ValueError):
            xor3.validate_operations(((xor3.OP_XOR3, 0, 1, 1),), 3)
        with self.assertRaises(ValueError):
            xor3.validate_operations(((4, 0, 1),), 2)

    def test_commutation_rule_matches_exhaustive_linear_maps(self):
        width = 4
        gates = tuple((destination, source)
                      for destination in range(width)
                      for source in range(width)
                      if destination != source)
        for first in gates:
            for second in gates:
                actually_commutes = (
                    xor3.circuit_matrix((first, second), width) ==
                    xor3.circuit_matrix((second, first), width))
                self.assertEqual(
                    xor3.cnot_gates_commute(first, second),
                    actually_commutes,
                    (first, second))

    def test_adjacent_pair_and_duplicate_cancellation(self):
        gates = ((0, 2), (0, 1), (3, 1), (3, 1))
        operations = xor3.pair_adjacent(gates, 4)
        self.assertEqual(operations, ((xor3.OP_XOR3, 0, 1, 2),))
        self.assert_schedule_equivalent(gates, operations, 4)

    def test_backward_commuting_pair(self):
        gates = ((0, 1), (2, 3), (0, 3))
        operations = xor3.pair_commuting_backward(gates, 4)
        self.assertEqual(operations, (
            (xor3.OP_XOR3, 0, 1, 3),
            (xor3.OP_XOR2, 2, 3),
        ))
        self.assert_schedule_equivalent(gates, operations, 4)

    def test_backward_scheduler_does_not_cross_dependency(self):
        gates = ((0, 1), (2, 3), (0, 2))
        operations = xor3.pair_commuting_backward(gates, 4)
        self.assertEqual(operations, (
            (xor3.OP_XOR2, 0, 1),
            (xor3.OP_XOR2, 2, 3),
            (xor3.OP_XOR2, 0, 2),
        ))
        self.assert_schedule_equivalent(gates, operations, 4)

    def test_forward_variant_finds_pair_blocked_in_backward_direction(self):
        gates = ((0, 1), (2, 3), (0, 2))
        backward = xor3.pair_commuting_backward(gates, 4)
        forward = xor3.pair_commuting_forward(gates, 4)
        self.assertEqual(len(backward), 3)
        self.assertEqual(forward, (
            (xor3.OP_XOR2, 2, 3),
            (xor3.OP_XOR3, 0, 1, 2),
        ))
        self.assert_schedule_equivalent(gates, forward, 4)

    def test_random_portfolio_variants_preserve_every_map(self):
        generator = random.Random(0xF8C10C17)
        # Eight-wire multipliers and sixteen-wire butterflies are both covered.
        for width in tuple(range(2, 9)) + (16,):
            for unused_case in range(100):
                gates = []
                for unused_gate in range(generator.randrange(0, 41)):
                    destination = generator.randrange(width)
                    source = generator.randrange(width - 1)
                    if source >= destination:
                        source += 1
                    gates.append((destination, source))
                gates = tuple(gates)
                for unused_variant, operations in xor3.schedule_variants(
                        gates, width):
                    self.assert_schedule_equivalent(gates, operations, width)

    def test_metadata_is_exact_and_explicit_about_liveness(self):
        operations = (
            (xor3.OP_XOR3, 0, 1, 2),
            (xor3.OP_XOR2, 3, 4),
        )
        metadata = xor3.schedule_metadata(operations, 5)
        self.assertEqual(metadata, {
            "xor2_count": 1,
            "xor3_count": 1,
            "instruction_count": 2,
            "depth": 1,
            "peak_live_wires": 5,
            "estimated_code_bytes": 13,
        })

    def test_selection_is_deterministic_and_uses_profitable_xor3(self):
        gates = ((0, 1), (2, 3), (0, 2), (1, 3), (1, 2))
        first = xor3.choose_schedule(gates, 4)
        for unused_iteration in range(100):
            self.assertEqual(xor3.choose_schedule(gates, 4), first)
        metadata = xor3.schedule_metadata(first[1], 4)
        self.assertGreater(metadata["xor3_count"], 0)
        self.assertLess(metadata["instruction_count"], len(gates))
        self.assert_schedule_equivalent(gates, first[1], 4)

    def test_expand_operations_preserves_ordered_xor_semantics(self):
        operations = (
            (xor3.OP_XOR3, 0, 2, 1),
            (xor3.OP_XOR2, 3, 0),
        )
        expanded = xor3.expand_operations(operations, 4)
        self.assertEqual(expanded, ((0, 2), (0, 1), (3, 0)))
        self.assertEqual(
            xor3.circuit_matrix(expanded, 4),
            xor3.operation_matrix(operations, 4))


if __name__ == "__main__":
    unittest.main(verbosity=2)
