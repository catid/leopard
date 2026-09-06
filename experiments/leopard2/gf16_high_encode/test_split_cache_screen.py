"""Pure, small screen checks: no benchmark child or clock is invoked."""
import copy
import json
from pathlib import Path
import unittest
from unittest import mock

from run_split_cache_screen import analyze, check_passive, validate, validate_plan


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

    def test_only_fixed_profiles_and_unchanged_method(self):
        for name in ("split_cache_screen_plan.json",
                     "split_cache_screen_foureyes_plan.json"):
            plan = json.loads(Path(__file__).with_name(name).read_text())
            validate_plan(plan, name)
            for key, value in (("cpu", 23), ("sibling", 87),
                               ("passive_seconds", 11), ("attempt_budget", 2),
                               ("rounds", 2), ("samples_per_process", 20),
                               ("order", ["on", "off", "off", "on"]),
                               ("cell_order", list(reversed(range(6))))):
                with self.assertRaises(ValueError):
                    validate_plan(dict(plan, **{key: value}), name)
            changed = copy.deepcopy(plan)
            changed["cells"][0]["bytes"] = 65536
            with self.assertRaises(ValueError):
                validate_plan(changed, name)
            if "host" in plan:
                changed = copy.deepcopy(plan)
                changed["host"]["hostname"] = "ripper"
                with self.assertRaises(ValueError):
                    validate_plan(changed, name)
        with self.assertRaises(ValueError):
            validate_plan(self.plan, "unregistered.json")

    def test_passive_records_success_contamination_and_short_window(self):
        plan = {"passive_seconds": 10, "sibling": 86}
        for after, elapsed, accepted in ((100, 10000000001, True),
                                         (101, 10000000001, False),
                                         (99, 10000000001, False),
                                         (100, 9999999999, False)):
            state = {}
            with mock.patch("run_split_cache_screen.sibling_ticks",
                            side_effect=[100, after]) as ticks, \
                 mock.patch("run_split_cache_screen.time.sleep") as sleep, \
                 mock.patch("run_split_cache_screen.time.monotonic_ns",
                            side_effect=[1, 1 + elapsed]):
                if accepted:
                    check_passive(state, plan)
                else:
                    with self.assertRaises(ValueError):
                        check_passive(state, plan)
                self.assertEqual(state["passive"], {
                    "before": 100, "after": after, "elapsed_ns": elapsed})
                self.assertEqual(ticks.call_args_list, [mock.call(86)] * 2)
                sleep.assert_called_once_with(10)
        state = {}
        with mock.patch("run_split_cache_screen.sibling_ticks") as ticks:
            check_passive(state, self.plan)
            ticks.assert_not_called()
            self.assertEqual(state, {})

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
