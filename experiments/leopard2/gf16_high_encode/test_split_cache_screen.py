"""Pure, small screen checks: no benchmark child or clock is invoked."""
import copy
import json
from pathlib import Path
import unittest

from run_split_cache_screen import analyze, validate


class ScreenTests(unittest.TestCase):
    def setUp(self):
        self.plan = json.loads(Path(__file__).with_name(
            "split_cache_screen_plan.json").read_text())
        self.rows = []
        for cell in self.plan["cells"]:
            for round_id in range(3):
                for slot, variant in enumerate(self.plan["order"]):
                    self.rows.append({"cell": cell["id"], "round": round_id,
                        "slot": slot, "variant": variant, "sibling_delta": 0,
                        "record": {"samples_ns": [100] * 21}})

    def test_neutral_is_not_a_gain(self):
        result = analyze(self.plan, self.rows)
        self.assertEqual(result["decision"], "reject_for_this_screen")
        self.assertFalse(result["production_promotion"])
        self.assertFalse(result["exact_leopard1_claim"])

    def test_gain_requires_controls_and_neighbors(self):
        for row in self.rows:
            if row["cell"] == 0 and row["variant"] == "on":
                row["record"]["samples_ns"] = [90] * 21
        self.assertEqual(analyze(self.plan, self.rows)["decision"],
                         "continue_to_future_qualification")
        for row in self.rows:
            if row["cell"] == 5 and row["variant"] == "on":
                row["record"]["samples_ns"] = [110] * 21
        self.assertEqual(analyze(self.plan, self.rows)["decision"],
                         "inconclusive_controls")
        for row in self.rows:
            if row["cell"] == 5:
                row["record"]["samples_ns"] = [100] * 21
            if row["cell"] == 2 and row["variant"] == "on":
                row["record"]["samples_ns"] = [110] * 21
        self.assertEqual(analyze(self.plan, self.rows)["decision"],
                         "reject_for_this_screen")

    def test_incomplete_contaminated_or_wrong_order_rejected(self):
        with self.assertRaises(ValueError):
            analyze(self.plan, self.rows[:-1])
        for field, value in (("sibling_delta", 1), ("variant", "on"),
                             ("slot", 2)):
            rows = copy.deepcopy(self.rows)
            rows[0][field] = value
            with self.assertRaises(ValueError):
                analyze(self.plan, rows)

    def test_record_shape_route_identity_and_samples(self):
        cell = self.plan["cells"][0]
        expected = {"schema": "leopard2-gf16-split-screen/v1", "cell": 0,
                    "k": 1000, "r": 200, "bytes": 32768,
                    "execution_route": "avx512", "input_hash": "a",
                    "output_hash": "b", "scratch_bytes": 42}
        record = dict(expected, samples_ns=[100] * 21)
        validate(record, cell, expected, True)
        validate(dict(expected, samples_ns=[]), cell, expected, False)
        for field, value in (("execution_route", "gfni"), ("output_hash", "x"),
                             ("scratch_bytes", 43), ("samples_ns", [True] * 21),
                             ("samples_ns", [0] * 21), ("samples_ns", [])):
            with self.assertRaises(ValueError):
                validate(dict(record, **{field: value}), cell, expected, True)


if __name__ == "__main__":
    unittest.main()
