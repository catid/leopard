#!/usr/bin/env python3
"""Single local same-binary AUTO extension qualification; no old samples."""
import fcntl
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

from run_gfni_boundary_screen import HOST, strict_equal, validate
from run_split_cache_screen import check_passive, digest, host_identity, require, sibling_ticks

PLAN = "auto_gfni_boundary_screen_plan.json"
ORDERS = {"off_on": ["off", "on", "on", "off"],
          "same_on": ["on_a", "on_b", "on_b", "on_a"]}
MAIN_ORDER = ["main", "on", "on", "main"]
CELLS = [
    {"k":1000,"r":200,"bytes":32768,"api":"encode","role":"target"},
    {"k":1000,"r":199,"bytes":65536,"api":"encode","role":"target"},
    {"k":1000,"r":200,"bytes":32768,"api":"one_item_batch","role":"target"},
    {"k":1000,"r":199,"bytes":65536,"api":"one_item_batch","role":"target"},
    {"k":1000,"r":200,"bytes":65536,"api":"encode","role":"unchanged_neighbor"},
    {"k":1000,"r":199,"bytes":32768,"api":"encode","role":"unchanged_neighbor"},
    {"k":1000,"r":200,"bytes":32768,"api":"explicit_avx2_encode","role":"unchanged_neighbor"},
    {"k":4096,"r":512,"bytes":4096,"api":"encode","role":"unchanged_neighbor"},
]
PROTOCOL = {"schema":"leopard-auto-gfni-boundary-plan/v1", "bead":"leopard-79h.38.5.4.17.1",
            "host":HOST, "cpu":26, "sibling":90, "controller_cpu":0,
            "condition":"slipgate-disabled-20260909", "attempt_budget":1,
            "rounds":3, "samples_per_process":21, "passive_seconds":10,
            "orders":ORDERS, "main_order":MAIN_ORDER, "main_cells":[0,1], "cells":CELLS,
            "minimum_gain":1.05, "equivalence_bound":1.02,
            "codec_commit":"35f53fb5c12f3f336f5a0e6235c5c243993a604a",
            "main_commit":"6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198",
            "core_sha256":"0e4cdc4485e96e6d1ca711ad0b5192925df01e7a4951bca3d40b35e580c061b7",
            "header_sha256":"0c6ff3efdfbb8754cd5abebaf520f05241f16e4f5487f27b093aae4dbd2ca065",
            "production_promotion":False, "authoritative_v19":False}


def orders(cell):
    return dict(ORDERS, main_on=MAIN_ORDER) if cell < 2 else ORDERS


def validate_plan(plan):
    require(set(plan) == set(PROTOCOL) | {"artifact_sha256"}, "plan fields")
    strict_equal({key:plan[key] for key in PROTOCOL}, PROTOCOL, "changed protocol")
    require(set(plan["artifact_sha256"]) == {"main","current","main.a","current.a","expected.json"},
            "artifact set")
    require(all(type(value) is str and len(value) == 64 and
                all(c in "0123456789abcdef" for c in value)
                for value in plan["artifact_sha256"].values()), "artifact hashes")


def analyze(rows):
    require(type(rows) is list and len(rows) == 216, "partial attempts have no analysis")
    cells, cursor = [], 0
    for cell in range(8):
        rounds = {name:[] for name in orders(cell)}
        for round_id in range(3):
            for comparison, order in orders(cell).items():
                medians = []
                for slot, variant in enumerate(order):
                    row = rows[cursor]
                    cursor += 1
                    strict_equal({key:row[key] for key in
                        ("cell","round","comparison","slot","variant","sibling_delta")},
                        {"cell":cell,"round":round_id,"comparison":comparison,"slot":slot,
                         "variant":variant,"sibling_delta":0}, "order/isolation")
                    values = row["record"]["samples_ns"]
                    require(type(values) is list and len(values) == 21 and
                            all(type(x) is int and x > 0 for x in values), "samples")
                    medians.append(statistics.median(values))
                rounds[comparison].append(math.sqrt((medians[0]/medians[1])*(medians[3]/medians[2])))
        cells.append({"cell":cell,"role":CELLS[cell]["role"],"round_ratios":rounds,
                      "ratios":{name:math.exp(statistics.mean(math.log(x) for x in values))
                                for name,values in rounds.items()}})
    controls = all(1/1.02 <= cell["ratios"]["same_on"] <= 1.02 for cell in cells)
    neighbors = all(1/1.02 <= cell["ratios"]["off_on"] <= 1.02 for cell in cells[4:])
    targets = all(cell["ratios"]["off_on"] >= 1.05 and
                  min(cell["round_ratios"]["off_on"]) > 1 for cell in cells[:4])
    decision = ("inconclusive_controls" if not controls else
                "reject_neighbor_gate" if not neighbors else
                "reject_target_gate" if not targets else "continue_to_production_integration")
    return {"cells":cells,"controls_pass":controls,"neighbors_pass":neighbors,
            "targets_pass":targets,"decision":decision,"confidence_intervals":False,
            "production_promotion":False,"authoritative_v19":False}


def run(bundle, output):
    output.mkdir(mode=0o700)  # Never overwrite, resume, retry or pool.
    plan = json.loads((bundle/PLAN).read_text())
    pins = json.loads((bundle/"pins.json").read_text())
    expected = json.loads((bundle/"expected.json").read_text())
    validate_plan(plan)
    strict_equal(host_identity(), HOST, "host changed")
    for cpu in (26,90):
        require(Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
                .read_text().strip() == "26,90", "topology changed")
    os.sched_setaffinity(0,{0})
    state = {"schema":"leopard-auto-gfni-boundary-attempt/v1","pins":pins,
             "plan_sha256":digest(bundle/PLAN),"host":host_identity(),
             "preflight":[],"invocations":[],"complete":False}
    locks = []

    def verify():
        for name,sha in pins["files"].items():
            require(Path(name).name == name, "unsafe input name")
            path = bundle/name
            require(path.is_file() and not path.is_symlink() and
                    not path.stat().st_mode & 0o222 and digest(path) == sha, "changed input: "+name)
        for name,sha in plan["artifact_sha256"].items():
            require(pins["files"][name] == sha, "plan artifact")
        require(pins["files"]["leopard2.cpp"] == plan["core_sha256"] and
                pins["files"]["Leopard2Direct.h"] == plan["header_sha256"], "codec source pins")

    try:
        leases = Path(f"/run/user/{os.getuid()}/leopard2-cpu-leases")
        require(leases.is_dir() and not leases.is_symlink() and leases.stat().st_uid == os.getuid()
                and leases.stat().st_mode & 0o777 == 0o700, "safe lease directory")
        for path in (Path("/tmp/leopard-gf8-authoritative.lock"),
                     leases/f"leopard2-cpu-pair-{os.getuid()}-26-90.lock"):
            fd = os.open(path,os.O_RDONLY|os.O_CREAT|os.O_NOFOLLOW,0o600)
            locks.append(fd)
            fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        verify()
        env = {"PATH":"/usr/bin:/bin","LANG":"C","LC_ALL":"C",
               "OMP_NUM_THREADS":"1","OMP_DYNAMIC":"FALSE","OMP_THREAD_LIMIT":"1"}

        def condition(name):
            with (output/(name+".stdout")).open("xb") as out, (output/(name+".stderr")).open("xb") as err:
                result = subprocess.run(["/bin/bash",str(bundle/"check-condition.sh")],
                                        stdout=out,stderr=err,timeout=30,check=False)
            require(result.returncode == 0 and (output/(name+".stderr")).stat().st_size == 0,
                    "shutdown condition: "+name)

        def invoke(variant,cell,name,measured):
            executable,mode = ("main","0") if variant == "main" else ("current","0" if variant == "off" else "1")
            command = ["/usr/bin/taskset","-c","26","/usr/bin/prlimit","--cpu=30:30",
                       "--fsize=1048576:1048576","--",str(bundle/executable),
                       "--measure" if measured else "--check",str(cell),mode]
            before = sibling_ticks(90)
            stdout,stderr = output/(name+".stdout"),output/(name+".stderr")
            with stdout.open("xb") as out, stderr.open("xb") as err:
                result = subprocess.run(command,stdout=out,stderr=err,env=env,timeout=60,check=False)
            delta = sibling_ticks(90)-before
            require(result.returncode == 0 and stderr.stat().st_size == 0, "child: "+name)
            require(stdout.stat().st_size <= 1048576, "record size")
            record = json.loads(stdout.read_text())
            verify()
            return record,delta

        condition("condition-before")
        for cell in range(8):
            for variant in ("main","off","on"):
                record,_ = invoke(variant,cell,f"check-{cell}-{variant}",False)
                validate(record,expected[variant][cell],False)
                state["preflight"].append(record)
        check_passive(state,plan)
        for cell in range(8):
            for round_id in range(3):
                for comparison,order in orders(cell).items():
                    for slot,variant in enumerate(order):
                        actual = "on" if variant.startswith("on") else variant
                        record,delta = invoke(actual,cell,f"cell-{cell}-round-{round_id}-{comparison}-slot-{slot}",True)
                        state["invocations"].append({"cell":cell,"round":round_id,"comparison":comparison,
                            "slot":slot,"variant":variant,"sibling_delta":delta,"record":record})
                        validate(record,expected[actual][cell],True)
                        require(delta == 0, "sibling activity; attempt stopped")
                print(f"cell {cell} round {round_id}: comparisons and control retained",flush=True)
        condition("condition-after")
        require((output/"condition-before.stdout").read_bytes() ==
                (output/"condition-after.stdout").read_bytes(), "shutdown state changed")
        verify()
        state["analysis"] = analyze(state["invocations"])
        state["complete"] = True
    except Exception as error:
        state["failure"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        (output/"attempt.json").write_text(json.dumps(state,indent=2)+"\n")
        for fd in reversed(locks): os.close(fd)


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: run_auto_gfni_boundary_screen.py frozen output")
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
