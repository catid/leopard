"""Small pure tests; no benchmark subprocess or timing is executed."""
import copy
import unittest

from run_current_route_screen import ORDERS, analyze, validate


class CurrentRouteScreenTests(unittest.TestCase):
    def setUp(self):
        self.rows = []
        for cell in range(6):
            for round_id in range(3):
                for comparison, order in ORDERS.items():
                    for slot, variant in enumerate(order):
                        self.rows.append({"cell": cell, "round": round_id,
                            "comparison": comparison, "slot": slot, "variant": variant,
                            "sibling_delta": 0, "record": {"samples_ns": [100] * 21}})

    def test_equal_is_uncertain_and_never_authoritative(self):
        result = analyze(self.rows)
        self.assertEqual(result["decision"], "diagnostic_complete")
        self.assertTrue(all(x["interpretation"] == "near_parity_or_uncertain"
                            for x in result["cells"]))
        for flag in ("confidence_intervals", "authoritative_v19", "production_promotion",
                     "historical_exact_main_gap_closed"):
            self.assertFalse(result[flag])

    def test_clear_directions_need_clean_controls(self):
        for row in self.rows:
            if row["comparison"] == "main_vs_current" and row["variant"] == "main":
                row["record"]["samples_ns"] = [90 if row["cell"] == 0 else 110] * 21
        result = analyze(self.rows)
        self.assertEqual(result["cells"][0]["interpretation"], "investigate_current_deficit")
        self.assertEqual(result["cells"][1]["interpretation"], "current_advantage")
        for row in self.rows:
            if row["comparison"] == "same_current" and row["variant"] == "current_a":
                row["record"]["samples_ns"] = [103] * 21
        result = analyze(self.rows)
        self.assertEqual(result["decision"], "inconclusive_controls")
        self.assertTrue(all(x["interpretation"] == "no_inference_controls"
                            for x in result["cells"]))

    def test_mixed_rounds_or_small_effect_do_not_classify(self):
        for row in self.rows:
            if row["comparison"] == "main_vs_current" and row["variant"] == "main":
                row["record"]["samples_ns"] = [99] * 21
        self.assertEqual(analyze(self.rows)["cells"][0]["interpretation"],
                         "near_parity_or_uncertain")
        for row in self.rows:
            if row["comparison"] == "main_vs_current" and row["variant"] == "main":
                row["record"]["samples_ns"] = [101 if row["round"] == 0 else 80] * 21
        self.assertEqual(analyze(self.rows)["cells"][0]["interpretation"],
                         "near_parity_or_uncertain")

    def test_partial_contaminated_and_reordered_fail(self):
        with self.assertRaises(ValueError):
            analyze(self.rows[:-1])
        for key, value in (("comparison", "same_current"), ("sibling_delta", 1),
                           ("variant", "current"), ("slot", 2), ("cell", 1)):
            changed = copy.deepcopy(self.rows)
            changed[0][key] = value
            with self.assertRaises(ValueError):
                analyze(changed)
        for value in ([True] * 21, [0] * 21, [1] * 20):
            changed = copy.deepcopy(self.rows)
            changed[0]["record"]["samples_ns"] = value
            with self.assertRaises(ValueError):
                analyze(changed)

    def test_source_route_scratch_and_samples_are_bound(self):
        expected = {"schema": "leopard-gf16-current-route-screen/v1", "codec_commit": "c",
                    "execution_route": "gfni", "scratch_bytes": 42, "output_hash": "a"}
        validate(dict(expected, samples_ns=[]), expected, False)
        validate(dict(expected, samples_ns=[1] * 21), expected, True)
        for key in expected:
            with self.assertRaises(ValueError):
                validate(dict(expected, samples_ns=[], **{key: "changed"}), expected, False)
        with self.assertRaises(ValueError):
            validate(dict(expected, samples_ns=[True] * 21), expected, True)


if __name__ == "__main__":
    unittest.main()
