"""Pure four-mode protocol tests: no benchmark, clock or subprocess execution."""
import copy
import unittest
from run_gfni_combined_screen import ARTIFACTS, CELLS, analyze, schedule, validate_plan, validate_trace

ORDERS = ((0,1,2,3,3,2,1,0), (1,2,3,0,0,3,2,1), (2,3,0,1,1,0,3,2))


def trace(cell, mode, measured, exercise=False):
    encodes = 26 if measured or exercise else 1
    matches = encodes * 2 if cell == 0 else 0
    return dict(schema="gfni-combined-timing/v1", cell=cell, mode=mode, encodes=encodes,
        calls=encodes * (2,2,1,1,2,1)[cell], matches=matches,
        first=matches if mode in (1,3) else 0, terminal=matches if mode in (2,3) else 0,
        timed=measured, exercise=exercise)


def plan_fixture():
    return dict(schema="gfni-combined-screen-plan/v1", bead="leopard-79h.38.5.4.16",
        artifact_sha256={name:"0"*64 for name in ARTIFACTS}, cpu=22, sibling=86,
        controller_cpu=0, passive_seconds=10, attempt_budget=1, rounds=3,
        samples_per_process=21, timed_invocations=288, seed=20260906,
        initial_untimed_encodes=1, additional_warmups=4, round_orders=ORDERS,
        comparisons=["factorial","same_off"],
        thresholds=dict(target_ratio=1.05, control_factor=1.02, positive_target_rounds=3),
        cells=[dict(id=i,k=k,r=r,bytes=size,current_route=route)
               for i,(k,r,size,route) in enumerate(CELLS)],
        host={"hostname":"foureyes", "kernel":"6.8.0-138-generic", "vendor_id":"AuthenticAMD",
              "cpu family":"26", "model":"8", "model name":"AMD Ryzen Threadripper PRO 9985WX 64-Cores"})


class CombinedScreenTests(unittest.TestCase):
    def setUp(self):
        self.rows=[]
        for cell in range(6):
            for round_id,order in enumerate(ORDERS):
                for comparison in ("factorial","same_off"):
                    for slot,variant in enumerate(order):
                        mode=variant if comparison=="factorial" else 0
                        self.rows.append(dict(cell=cell,round=round_id,comparison=comparison,
                            slot=slot,variant=variant,mode=mode,sibling_delta=0,
                            record={"samples_ns":[100]*21},trace=trace(cell,mode,True)))

    def test_fixed_schedule_and_mirror_balance(self):
        actual=list(schedule())
        self.assertEqual(len(actual),288)
        self.assertEqual(actual,[{key:row[key] for key in actual[0]} for row in self.rows])
        for order in ORDERS:
            for mode in range(4):
                slots=[i for i,value in enumerate(order) if value==mode]
                self.assertEqual(len(slots),2)
                self.assertEqual(sum(slots),7)

    def test_equal_rejects_without_promotion(self):
        result=analyze(self.rows)
        self.assertEqual(result["decision"],"reject_for_this_screen")
        self.assertEqual(result["aggregate_controls"],33)
        self.assertEqual(result["cells"][0]["interaction_factor"],1)
        for flag in ("confidence_intervals","production_promotion","exact_leopard1_claim","authoritative_v19"):
            self.assertFalse(result[flag])

    def test_target_and_every_control_mode(self):
        for row in self.rows:
            if row["cell"]==0 and row["comparison"]=="factorial" and row["mode"]==3:
                row["record"]["samples_ns"]=[94]*21
        self.assertEqual(analyze(self.rows)["decision"],"continue_to_future_qualification")
        for comparison in ("factorial","same_off"):
            for mode in (1,2,3):
                for value in (97,103):
                    changed=copy.deepcopy(self.rows)
                    for row in changed:
                        if row["cell"]==1 and row["comparison"]==comparison and row["variant"]==mode:
                            row["record"]["samples_ns"]=[value]*21
                    self.assertEqual(analyze(changed)["decision"],"inconclusive_controls")

    def test_small_or_mixed_target_rejected(self):
        for row in self.rows:
            if row["cell"]==0 and row["comparison"]=="factorial" and row["mode"]==3:
                row["record"]["samples_ns"]=[97]*21
        self.assertEqual(analyze(self.rows)["decision"],"reject_for_this_screen")
        for row in self.rows:
            if row["cell"]==0 and row["comparison"]=="factorial" and row["mode"]==3:
                row["record"]["samples_ns"]=[101 if row["round"]==0 else 80]*21
        self.assertEqual(analyze(self.rows)["decision"],"reject_for_this_screen")

    def test_interaction_from_same_fresh_rounds(self):
        for both in (900,1000,1100):
            for row in self.rows:
                if row["cell"]==0 and row["comparison"]=="factorial":
                    row["record"]["samples_ns"]=[(1320,1200,1100,both)[row["mode"]]]*21
            result=analyze(self.rows)["cells"][0]
            self.assertAlmostEqual(result["ratios"]["factorial"]["1"],1.1)
            self.assertAlmostEqual(result["ratios"]["factorial"]["2"],1.2)
            self.assertAlmostEqual(result["ratios"]["factorial"]["3"],1320/both)
            self.assertAlmostEqual(result["interaction_factor"],1000/both)
            for value in result["interaction_rounds"]:
                self.assertAlmostEqual(value,1000/both)

    def test_partial_order_isolation_and_samples(self):
        with self.assertRaises(ValueError): analyze(self.rows[:-1])
        with self.assertRaises(ValueError): analyze(self.rows+self.rows[:1])
        for key,value in (("cell",1),("round",1),("comparison","same_off"),("slot",1),
                          ("variant",True),("mode",True),("sibling_delta",True),("sibling_delta",1)):
            changed=copy.deepcopy(self.rows); changed[0][key]=value
            with self.assertRaises(ValueError): analyze(changed)
        for values in ([True]*21,[0]*21,[1]*20,[1.0]*21):
            changed=copy.deepcopy(self.rows); changed[0]["record"]["samples_ns"]=values
            with self.assertRaises(ValueError): analyze(changed)

    def test_all_trace_fields_and_types(self):
        for cell in range(6):
            for mode in range(4):
                for measured,exercise in ((False,False),(True,False),(False,True)):
                    original=trace(cell,mode,measured,exercise)
                    validate_trace(original,cell,mode,measured,exercise)
                    for key in original:
                        changed=dict(original,**{key:None})
                        with self.assertRaises(ValueError): validate_trace(changed,cell,mode,measured,exercise)
                    for key,value in (("mode",bool(mode)),("calls",True),("timed",int(measured)),("extra",0)):
                        changed=dict(original,**{key:value})
                        with self.assertRaises(ValueError): validate_trace(changed,cell,mode,measured,exercise)
        for mode in (True,-1,4):
            with self.assertRaises(ValueError): validate_trace(trace(0,0,False),0,mode,False)
        with self.assertRaises(ValueError): validate_trace(trace(0,0,True),0,0,True,True)

    def test_bad_trace_blocks_analysis(self):
        self.rows[1]["trace"]["first"]=0
        with self.assertRaises(ValueError): analyze(self.rows)

    def test_fixed_plan_and_full_inventory(self):
        plan=plan_fixture(); validate_plan(plan)
        for key,value in (("cpu",23),("sibling",87),("controller_cpu",1),("passive_seconds",0),
            ("attempt_budget",True),("rounds",4),("samples_per_process",20),("timed_invocations",144),
            ("seed",0),("initial_untimed_encodes",0),("additional_warmups",0),("schema","other"),
            ("bead","other"),("artifact_sha256",{}),("round_orders",[]),("comparisons",[]),
            ("cells",[]),("host",{}),("thresholds",{})):
            with self.assertRaises(ValueError): validate_plan(dict(plan,**{key:value}))
        for key in plan["thresholds"]:
            changed=copy.deepcopy(plan); changed["thresholds"][key]=0
            with self.assertRaises(ValueError): validate_plan(changed)
        for name in ARTIFACTS:
            changed=copy.deepcopy(plan); del changed["artifact_sha256"][name]
            with self.assertRaises(ValueError): validate_plan(changed)

    def test_independent_log_replay_matches_collector(self):
        from replay_gfni_combined_screen import derive,close
        for target in ((100,100,100,100),(1320,1200,1100,900),(1320,1200,1100,1100)):
            for contaminated in (False,True):
                rows=copy.deepcopy(self.rows)
                medians={}
                for row in rows:
                    value=target[row["mode"]] if row["cell"]==0 and row["comparison"]=="factorial" else 100
                    if contaminated and row["cell"]==1 and row["comparison"]=="same_off" and row["variant"]==1:
                        value=103
                    row["record"]["samples_ns"]=[value]*21
                    medians[row["cell"],row["round"],row["comparison"],row["slot"]]=value
                close(analyze(rows),derive(medians))

    def test_independent_comparison_rejects_bad_fields(self):
        from replay_gfni_combined_screen import close
        original=analyze(self.rows)
        for key,value in (("decision","continue_to_future_qualification"),("aggregate_controls",True),
                           ("production_promotion",True),("cells",[]),("extra",0)):
            with self.assertRaises(ValueError): close(dict(original,**{key:value}),original)
        for value in (True,1.1,float("nan"),float("inf")):
            changed=copy.deepcopy(original); changed["cells"][0]["interaction_factor"]=value
            with self.assertRaises(ValueError): close(changed,original)
if __name__=="__main__": unittest.main()
