#!/usr/bin/env python3
"""Pure protocol tests; never launches a codec or reads a benchmark clock."""
import copy
import unittest
from run_gfni_boundary_screen import (ORDERS, PROTOCOL, analyze, strict_equal,
                                      validate, validate_plan)


def rows():
    result = []
    for cell in range(2):
        for round_id in range(3):
            for comparison, order in ORDERS.items():
                for slot, variant in enumerate(order):
                    result.append({"cell": cell, "round": round_id,
                        "comparison": comparison, "slot": slot, "variant": variant,
                        "sibling_delta": 0, "record": {"samples_ns":
                            [120 if variant in ("auto", "main") else 100] * 21}})
    return result


class ProtocolTests(unittest.TestCase):
    def test_positive_analysis(self):
        result = analyze(rows())
        self.assertTrue(result["controls_pass"])
        for cell in result["cells"]:
            self.assertAlmostEqual(cell["ratios"]["auto_vs_gfni"], 1.2)
            self.assertAlmostEqual(cell["ratios"]["main_vs_gfni"], 1.2)
            self.assertEqual(cell["decision"], "qualify_bounded_auto_candidate")
        self.assertFalse(result["production_promotion"])
        self.assertFalse(result["neighbor_qualification"])

    def test_control_failure_suppresses_every_cell(self):
        changed = rows()
        for row in changed:
            if row["cell"] == 1 and row["variant"] == "gfni_a":
                row["record"]["samples_ns"] = [110] * 21
        result = analyze(changed)
        self.assertFalse(result["controls_pass"])
        self.assertTrue(all(cell["decision"] == "inconclusive_controls"
                            for cell in result["cells"]))

    def test_five_percent_gate(self):
        changed = rows()
        for row in changed:
            if row["variant"] == "auto":
                row["record"]["samples_ns"] = [104] * 21
        self.assertTrue(all(cell["decision"] == "reject_for_this_screen"
                            for cell in analyze(changed)["cells"]))

    def test_rejects_invalid_rows(self):
        for key, value in (("cell", 1), ("round", 1), ("slot", 1),
                           ("sibling_delta", 1), ("sibling_delta", False),
                           ("comparison", "main_vs_gfni"), ("variant", "gfni")):
            changed = rows()
            changed[0][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                analyze(changed)
        with self.assertRaises(ValueError):
            analyze(rows()[:-1])
        for samples in ([True] * 21, [0] * 21, [1] * 20, [1.0] * 21):
            changed = rows()
            changed[0]["record"]["samples_ns"] = samples
            with self.assertRaises(ValueError):
                analyze(changed)

    def test_record_validation(self):
        expected = {"codec_commit": "pinned", "cell": 0, "requested": "gfni"}
        validate(dict(expected, samples_ns=[]), expected, False)
        validate(dict(expected, samples_ns=[100] * 21), expected, True)
        for key, value in (("codec_commit", "other"), ("cell", False),
                           ("requested", "auto"), ("samples_ns", [True] * 21)):
            changed = dict(expected, samples_ns=[100] * 21)
            changed[key] = value
            with self.assertRaises(ValueError):
                validate(changed, expected, True)

    def test_frozen_plan(self):
        plan = copy.deepcopy(PROTOCOL)
        plan["artifact_sha256"] = {key: "a" * 64 for key in
                                  ("main", "current", "main.a", "current.a", "expected.json")}
        validate_plan(plan)
        for key, value in (("cpu", 22), ("sibling", 86), ("controller_cpu", 1),
                           ("attempt_budget", 2), ("rounds", 4), ("rounds", 3.0),
                           ("samples_per_process", 20), ("passive_seconds", 0),
                           ("condition", "old-condition"), ("minimum_gain", 1.04),
                           ("control_bound", 1.03), ("production_promotion", 0),
                           ("neighbor_qualification", True), ("host", {}),
                           ("cells", list(reversed(plan["cells"]))), ("orders", {}),
                           ("artifact_sha256", {"main": "invalid"})):
            changed = copy.deepcopy(plan)
            changed[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_plan(changed)

    def test_strict_identity(self):
        with self.assertRaises(ValueError):
            strict_equal({"cell": 0}, {"cell": False}, "type")


if __name__ == "__main__":
    unittest.main()
