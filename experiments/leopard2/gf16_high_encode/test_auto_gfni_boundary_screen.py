#!/usr/bin/env python3
"""Pure tests for the new AUTO qualification protocol, without codec execution."""
import copy
import unittest
from run_auto_gfni_boundary_screen import PROTOCOL, analyze, orders, validate, validate_plan


def rows():
    result = []
    for cell in range(8):
        for round_id in range(3):
            for comparison, order in orders(cell).items():
                for slot, variant in enumerate(order):
                    result.append({"cell":cell,"round":round_id,"comparison":comparison,
                        "slot":slot,"variant":variant,"sibling_delta":0,"record":{"samples_ns":
                            [140 if (cell < 4 and variant == "off") or variant == "main" else 100]*21}})
    return result


class Tests(unittest.TestCase):
    def test_positive(self):
        data = rows()
        self.assertEqual(len(data),216)
        result = analyze(data)
        self.assertEqual(result["decision"],"continue_to_production_integration")
        self.assertTrue(result["targets_pass"] and result["neighbors_pass"] and result["controls_pass"])
        self.assertFalse(result["production_promotion"])
        for cell in result["cells"][:4]: self.assertAlmostEqual(cell["ratios"]["off_on"],1.4)
        self.assertNotIn("main_on",result["cells"][2]["ratios"])

    def test_control_failure(self):
        data = rows()
        for row in data:
            if row["cell"] == 7 and row["variant"] == "on_a": row["record"]["samples_ns"] = [105]*21
        self.assertEqual(analyze(data)["decision"],"inconclusive_controls")

    def test_neighbor_equivalence(self):
        for value in (95,105):
            data = rows()
            for row in data:
                if row["cell"] == 5 and row["variant"] == "off": row["record"]["samples_ns"] = [value]*21
            self.assertEqual(analyze(data)["decision"],"reject_neighbor_gate")

    def test_all_four_targets_required(self):
        for cell in range(4):
            data = rows()
            for row in data:
                if row["cell"] == cell and row["variant"] == "off": row["record"]["samples_ns"] = [104]*21
            self.assertEqual(analyze(data)["decision"],"reject_target_gate")

    def test_rows_fail_closed(self):
        with self.assertRaises(ValueError): analyze(rows()[:-1])
        for key,value in (("cell",1),("round",1),("slot",1),("sibling_delta",False),
                          ("sibling_delta",1),("comparison","main_on"),("variant","on")):
            data = rows(); data[0][key] = value
            with self.assertRaises(ValueError): analyze(data)
        for values in ([True]*21,[0]*21,[100.0]*21,[100]*20):
            data = rows(); data[0]["record"]["samples_ns"] = values
            with self.assertRaises(ValueError): analyze(data)

    def test_workload_identity(self):
        expected = {"cell":0,"codec_commit":"pinned","boundary_mode":1,"untimed_route_calls":1,"api":"leo2_encode"}
        validate(dict(expected,samples_ns=[]),expected,False)
        validate(dict(expected,samples_ns=[100]*21),expected,True)
        for key,value in (("cell",False),("codec_commit","other"),("boundary_mode",0),
                          ("untimed_route_calls",0),("api","leo_encode")):
            data = dict(expected,samples_ns=[100]*21); data[key] = value
            with self.assertRaises(ValueError): validate(data,expected,True)

    def test_plan(self):
        plan = copy.deepcopy(PROTOCOL)
        plan["artifact_sha256"] = {name:"a"*64 for name in ("main","current","main.a","current.a","expected.json")}
        validate_plan(plan)
        for key,value in (("cpu",22),("sibling",86),("controller_cpu",1),("attempt_budget",2),
                          ("rounds",4),("rounds",3.0),("samples_per_process",20),("passive_seconds",0),
                          ("minimum_gain",1.04),("equivalence_bound",1.03),("condition","other"),
                          ("codec_commit","other"),("core_sha256","other"),("header_sha256","other"),
                          ("main_cells",[0,1,2]),("main_order",[]),("orders",{}),("host",{}),
                          ("cells",list(reversed(plan["cells"]))),("production_promotion",0),
                          ("artifact_sha256",{})):
            changed = copy.deepcopy(plan); changed[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): validate_plan(changed)


if __name__ == "__main__": unittest.main()
