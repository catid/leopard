"""Pure protocol tests. No benchmark or clock is executed."""
import copy
import unittest

from run_gfni_source_stage_screen import ARTIFACTS, CELLS, ORDERS, analyze, validate_plan, validate_trace


def trace(cell, enabled, measured, exercise=False):
    encodes = 26 if measured or exercise else 1
    matches = encodes * 2 if cell == 0 else 0
    return {"schema": "gfni-source-stage-timing/v1", "cell": cell, "enabled": enabled,
        "encodes": encodes, "calls": encodes * (2, 2, 1, 1, 2, 1)[cell],
        "matches": matches, "changed": matches if enabled else 0,
        "timed": measured, "exercise": exercise}


class StageScreenTests(unittest.TestCase):
    def setUp(self):
        self.rows = []
        for cell in range(6):
            for round_id in range(3):
                for comparison, order in ORDERS.items():
                    for slot, variant in enumerate(order):
                        self.rows.append({"cell": cell, "round": round_id,
                            "comparison": comparison, "slot": slot, "variant": variant,
                            "sibling_delta": 0, "record": {"samples_ns": [100] * 21},
                            "trace": trace(cell, variant == "on", True)})

    def test_equal_rejects_and_never_promotes(self):
        result = analyze(self.rows)
        self.assertEqual(result["decision"], "reject_for_this_screen")
        for flag in ("confidence_intervals", "production_promotion", "exact_leopard1_claim",
                     "authoritative_v19"):
            self.assertFalse(result[flag])

    def test_target_and_control_gates(self):
        for row in self.rows:
            if row["cell"] == 0 and row["variant"] == "off":
                row["record"]["samples_ns"] = [106] * 21
        self.assertEqual(analyze(self.rows)["decision"], "continue_to_future_qualification")
        for comparison, variant in (("same_current", "off_a"), ("current_vs_staged", "off")):
            for sample in (97, 103):
                changed = copy.deepcopy(self.rows)
                for row in changed:
                    if row["cell"] == 1 and row["comparison"] == comparison and row["variant"] == variant:
                        row["record"]["samples_ns"] = [sample] * 21
                self.assertEqual(analyze(changed)["decision"], "inconclusive_controls")

    def test_small_or_mixed_target_rejected(self):
        for row in self.rows:
            if row["cell"] == 0 and row["variant"] == "off":
                row["record"]["samples_ns"] = [104] * 21
        self.assertEqual(analyze(self.rows)["decision"], "reject_for_this_screen")
        for row in self.rows:
            if row["cell"] == 0 and row["variant"] == "off":
                row["record"]["samples_ns"] = [99 if row["round"] == 0 else 120] * 21
        self.assertEqual(analyze(self.rows)["decision"], "reject_for_this_screen")

    def test_partial_order_isolation_and_samples_rejected(self):
        with self.assertRaises(ValueError):
            analyze(self.rows[:-1])
        for key, value in (("cell", 1), ("round", 1), ("slot", 1), ("variant", "on"),
                           ("comparison", "same_current"), ("sibling_delta", 1)):
            changed = copy.deepcopy(self.rows)
            changed[0][key] = value
            with self.assertRaises(ValueError):
                analyze(changed)
        for value in ([True] * 21, [0] * 21, [1] * 20, [1.0] * 21):
            changed = copy.deepcopy(self.rows)
            changed[0]["record"]["samples_ns"] = value
            with self.assertRaises(ValueError):
                analyze(changed)

    def test_trace_every_cell_mode_phase_and_field(self):
        for cell in range(6):
            for enabled in (False, True):
                for measured, exercise in ((False, False), (True, False), (False, True)):
                    expected = trace(cell, enabled, measured, exercise)
                    validate_trace(expected, cell, enabled, measured, exercise)
                    for key in expected:
                        changed = dict(expected, **{key: "changed"})
                        with self.assertRaises(ValueError):
                            validate_trace(changed, cell, enabled, measured, exercise)
                    for key, value in (("calls", True), ("enabled", int(enabled)),
                                       ("timed", int(measured)), ("extra", 0)):
                        with self.assertRaises(ValueError):
                            validate_trace(dict(expected, **{key: value}), cell, enabled, measured, exercise)

    def test_trace_corruption_blocks_analysis(self):
        self.rows[1]["trace"]["changed"] = 0
        with self.assertRaises(ValueError):
            analyze(self.rows)

    def test_fixed_plan_and_complete_artifact_inventory(self):
        plan = {"schema": "gfni-source-stage-screen-plan/v1", "bead": "leopard-79h.38.5.4.11",
            "artifact_sha256": {name: "0" * 64 for name in ARTIFACTS}, "cpu": 22, "sibling": 86,
            "controller_cpu": 0, "passive_seconds": 10, "attempt_budget": 1, "rounds": 3,
            "samples_per_process": 21, "orders": ORDERS,
            "cells": [dict(id=i, k=k, r=r, bytes=size, current_route=route)
                      for i, (k, r, size, route) in enumerate(CELLS)],
            "host": {"hostname": "foureyes", "kernel": "6.8.0-138-generic",
                "vendor_id": "AuthenticAMD", "cpu family": "26", "model": "8",
                "model name": "AMD Ryzen Threadripper PRO 9985WX 64-Cores"}}
        validate_plan(plan)
        for key, value in (("cpu", 23), ("sibling", 87), ("controller_cpu", 1),
                           ("passive_seconds", 0), ("attempt_budget", True), ("rounds", 4),
                           ("samples_per_process", 20), ("schema", "other"), ("bead", "other"),
                           ("artifact_sha256", {}), ("orders", {}), ("cells", []), ("host", {})):
            with self.assertRaises(ValueError):
                validate_plan(dict(plan, **{key: value}))
        for name in ARTIFACTS:
            changed = copy.deepcopy(plan)
            del changed["artifact_sha256"][name]
            with self.assertRaises(ValueError):
                validate_plan(changed)


if __name__ == "__main__":
    unittest.main()
