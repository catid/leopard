#!/usr/bin/env python3
"""Pure protocol tests; no native programs or benchmark clocks."""
import copy
import unittest

from run_r199_boundary_screen import (ARTIFACTS, FROZEN_FILES, ORDERS, PROFILES, PROTOCOL, analyze,
                                 identity, same, schedule, validate, validate_pins, validate_plan)


def fixtures():
    expected = {p: [dict(profile=p, cell=c, output_hash="fixed") for c in range(1)]
                for p in PROFILES}
    durations = {"native": 100, "auto": 110, "gfni": 80}
    rows = [dict(item, sibling_delta=0,
                 record=dict(expected[item["profile"]][item["cell"]],
                             encode_calls=26, samples_ns=[durations[item["profile"]]] * 21))
            for item in schedule()]
    return expected, rows


class ProtocolTests(unittest.TestCase):
    def test_ratios_and_no_promotion(self):
        expected, rows = fixtures()
        report = analyze(rows, expected)
        self.assertTrue(report["controls_pass"])
        for cell in report["cells"]:
            self.assertAlmostEqual(cell["ratios"]["auto_over_native"], 100 / 110)
            self.assertAlmostEqual(cell["ratios"]["gfni_over_native"], 100 / 80)
            self.assertAlmostEqual(cell["ratios"]["gfni_over_auto"], 110 / 80)
            for profile in PROFILES:
                self.assertAlmostEqual(cell["ratios"]["same_" + profile], 1)
        self.assertFalse(report["production_promotion"])
        self.assertFalse(report["authoritative_v19"])

    def test_each_control_can_suppress_conclusions(self):
        for cell in range(1):
            for profile in PROFILES:
                expected, rows = fixtures()
                for row in rows:
                    if row["cell"] == cell and row["comparison"] == "same_" + profile and row["slot"] in (0, 3):
                        row["record"]["samples_ns"] = [200] * 21
                result = analyze(rows, expected)
                self.assertFalse(result["controls_pass"])
                self.assertEqual(result["decision"], "inconclusive_controls")

    def test_incomplete_or_reordered(self):
        expected, rows = fixtures()
        self.assertEqual(len(rows), 72)
        for changed in (rows[:-1], rows + [rows[0]], rows[1:] + rows[:1]):
            with self.assertRaises(ValueError):
                analyze(changed, expected)

    def test_row_mutations(self):
        expected, rows = fixtures()
        for field, value in (("cell", False), ("cell", 1), ("round", 1),
                             ("slot", 1), ("profile", "auto"),
                             ("comparison", "same_native"), ("sibling_delta", 1),
                             ("sibling_delta", False)):
            changed = copy.deepcopy(rows)
            changed[0][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                analyze(changed, expected)
        for field, value in (("profile", "other"), ("output_hash", "changed"),
                             ("encode_calls", 1), ("encode_calls", 26.0),
                             ("samples_ns", [True] * 21), ("samples_ns", [0] * 21),
                             ("samples_ns", [1.0] * 21), ("samples_ns", [1] * 20)):
            changed = copy.deepcopy(rows)
            changed[0]["record"][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                analyze(changed, expected)

    def test_untimed_modes(self):
        expected = dict(profile="auto", cell=0)
        validate(dict(expected, samples_ns=[], encode_calls=1), expected, False)
        validate(dict(expected, samples_ns=[], encode_calls=26), expected, False, True)
        for calls in (1, 25, 27, True):
            with self.assertRaises(ValueError):
                validate(dict(expected, samples_ns=[], encode_calls=calls), expected, False, True)

    def test_plan_mutations(self):
        plan = copy.deepcopy(PROTOCOL)
        plan["artifact_sha256"] = {k: "a" * 64 for k in ARTIFACTS}
        validate_plan(plan)
        for key, value in (("cpu", 4), ("sibling", 68), ("controller_cpu", 1),
                           ("attempt_budget", 2), ("rounds", 4), ("rounds", 3.0),
                           ("samples_per_process", 20), ("passive_seconds", 0),
                           ("condition", "old"), ("control_bound", 1.03),
                           ("future_candidate_minimum_gain", 1.04), ("orders", {}),
                           ("cells", []), ("host", {}),
                           ("production_promotion", 0), ("confidence_intervals", True),
                           ("artifact_sha256", {"native": "bad"})):
            changed = copy.deepcopy(plan)
            changed[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_plan(changed)

    def test_fixed_gain_and_every_round_gates(self):
        expected, rows = fixtures()
        for row in rows:
            if row["profile"] == "auto":
                row["record"]["samples_ns"] = [104] * 21
            elif row["profile"] == "gfni":
                row["record"]["samples_ns"] = [100] * 21
        self.assertEqual(analyze(rows, expected)["decision"], "reject_for_this_screen")
        expected, rows = fixtures()
        for row in rows:
            if row["round"] == 1 and row["comparison"] == "gfni_over_auto" and row["profile"] == "auto":
                row["record"]["samples_ns"] = [70] * 21
        self.assertEqual(analyze(rows, expected)["decision"], "reject_for_this_screen")

    def test_strict_types(self):
        with self.assertRaises(ValueError):
            same(0, False)
        self.assertEqual(identity(dict(cell=0, encode_calls=1, samples_ns=[])), dict(cell=0))

    def test_frozen_inventory(self):
        good = dict(schema='leopard-r199-boundary-pins/v1', files={k:'a'*64 for k in FROZEN_FILES})
        validate_pins(good)
        for name in FROZEN_FILES:
            bad = copy.deepcopy(good)
            del bad['files'][name]
            with self.assertRaises(ValueError): validate_pins(bad)
        for name,value in (('../escape','a'*64),('native',None),('gfni','not-a-hash')):
            bad = copy.deepcopy(good)
            bad['files'][name] = value
            with self.assertRaises(ValueError): validate_pins(bad)


if __name__ == "__main__":
    unittest.main()
