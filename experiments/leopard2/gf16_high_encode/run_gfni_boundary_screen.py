#!/usr/bin/env python3
"""One local explicit-GFNI boundary screen, not AUTO selector promotion."""
import fcntl
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

from run_split_cache_screen import (check_passive, digest, host_identity,
                                    identity, require, sibling_ticks)

PLAN = "gfni_boundary_screen_plan.json"
ORDERS = {"auto_vs_gfni": ["auto", "gfni", "gfni", "auto"],
          "main_vs_gfni": ["main", "gfni", "gfni", "main"],
          "same_gfni": ["gfni_a", "gfni_b", "gfni_b", "gfni_a"]}
HOST = {"hostname": "work", "kernel": "6.8.0-137-generic",
        "vendor_id": "AuthenticAMD", "cpu family": "26", "model": "8",
        "model name": "AMD Ryzen Threadripper 9980X 64-Cores"}
PROTOCOL = {"schema": "leopard-gfni-boundary-plan/v1", "bead": "leopard-79h.38.5.4.17",
            "host": HOST, "cpu": 26, "sibling": 90, "controller_cpu": 0,
            "condition": "slipgate-disabled-20260909", "attempt_budget": 1,
            "rounds": 3, "samples_per_process": 21, "passive_seconds": 10,
            "orders": ORDERS,
            "cells": [{"k": 1000, "r": 200, "bytes": 32768},
                      {"k": 1000, "r": 199, "bytes": 65536}],
            "minimum_gain": 1.05, "control_bound": 1.02,
            "production_promotion": False, "neighbor_qualification": False}


def strict_equal(actual, expected, message):
    require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), message)


def validate_plan(plan):
    require(set(plan) == set(PROTOCOL) | {"artifact_sha256"}, "plan fields")
    strict_equal({key: plan[key] for key in PROTOCOL}, PROTOCOL, "changed protocol")
    require(set(plan["artifact_sha256"]) ==
            {"main", "current", "main.a", "current.a", "expected.json"}, "artifact set")
    require(all(type(value) is str and len(value) == 64 and
                all(c in "0123456789abcdef" for c in value)
                for value in plan["artifact_sha256"].values()), "artifact digests")


def validate(record, expected, measured):
    strict_equal(identity(record), expected, "changed workload/route/source/scratch")
    samples = record.get("samples_ns")
    require(type(samples) is list and len(samples) == (21 if measured else 0) and
            all(type(x) is int and x > 0 for x in samples), "sample values/count")


def analyze(rows):
    require(type(rows) is list and len(rows) == 72, "partial attempts have no analysis")
    cells = []
    cursor = 0
    for cell in range(2):
        rounds = {name: [] for name in ORDERS}
        for round_id in range(3):
            for comparison, order in ORDERS.items():
                medians = []
                for slot, variant in enumerate(order):
                    row = rows[cursor]
                    cursor += 1
                    strict_equal({key: row[key] for key in
                        ("cell", "round", "comparison", "slot", "variant", "sibling_delta")},
                        {"cell": cell, "round": round_id, "comparison": comparison,
                         "slot": slot, "variant": variant, "sibling_delta": 0}, "order/isolation")
                    samples = row["record"]["samples_ns"]
                    require(type(samples) is list and len(samples) == 21 and
                            all(type(x) is int and x > 0 for x in samples), "analysis samples")
                    medians.append(statistics.median(samples))
                rounds[comparison].append(math.sqrt(
                    (medians[0] / medians[1]) * (medians[3] / medians[2])))
        cells.append({"cell": cell, "round_ratios": rounds,
                      "ratios": {name: math.exp(statistics.mean(math.log(x) for x in values))
                                 for name, values in rounds.items()}})
    controls_ok = all(1 / 1.02 <= cell["ratios"]["same_gfni"] <= 1.02 for cell in cells)
    for cell in cells:
        cell["decision"] = ("inconclusive_controls" if not controls_ok else
            "qualify_bounded_auto_candidate" if cell["ratios"]["auto_vs_gfni"] >= 1.05
            and min(cell["round_ratios"]["auto_vs_gfni"]) > 1 else "reject_for_this_screen")
    return {"controls_pass": controls_ok, "cells": cells, "confidence_intervals": False,
            "production_promotion": False, "neighbor_qualification": False,
            "authoritative_v19": False}


def run(bundle, output):
    output.mkdir(mode=0o700)  # Exactly one attempt; never resume or overwrite.
    plan = json.loads((bundle / PLAN).read_text())
    pins = json.loads((bundle / "pins.json").read_text())
    expected = json.loads((bundle / "expected.json").read_text())
    validate_plan(plan)
    strict_equal(host_identity(), HOST, "host changed")
    for cpu in (26, 90):
        require(Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
                .read_text().strip() == "26,90", "topology changed")
    os.sched_setaffinity(0, {0})
    state = {"schema": "leopard-gfni-boundary-attempt/v1", "pins": pins,
             "plan_sha256": digest(bundle / PLAN), "host": host_identity(),
             "preflight": [], "invocations": [], "complete": False}
    locks = []

    def verify():
        for name, sha in pins["files"].items():
            require(Path(name).name == name, "unsafe frozen name")
            path = bundle / name
            require(path.is_file() and not path.is_symlink() and
                    not path.stat().st_mode & 0o222 and digest(path) == sha,
                    "input changed: " + name)
        for name, sha in plan["artifact_sha256"].items():
            require(pins["files"][name] == sha, "plan artifact: " + name)

    try:
        lease_root = Path(f"/run/user/{os.getuid()}/leopard2-cpu-leases")
        require(lease_root.is_dir() and not lease_root.is_symlink() and
                lease_root.stat().st_uid == os.getuid() and
                lease_root.stat().st_mode & 0o777 == 0o700, "unsafe lease directory")
        for path in (Path("/tmp/leopard-gf8-authoritative.lock"),
                     lease_root / f"leopard2-cpu-pair-{os.getuid()}-26-90.lock"):
            fd = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_NOFOLLOW, 0o600)
            locks.append(fd)
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        verify()
        env = {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
               "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "1"}

        def condition(label):
            with (output / (label + ".stdout")).open("xb") as out, \
                    (output / (label + ".stderr")).open("xb") as err:
                result = subprocess.run(["/bin/bash", str(bundle / "check-condition.sh")],
                                        stdout=out, stderr=err, timeout=30, check=False)
            require(result.returncode == 0 and (output / (label + ".stderr")).stat().st_size == 0,
                    "shutdown condition: " + label)

        def invoke(variant, cell, label, measured):
            mode = "--measure" if measured else "--check"
            argv = [str(bundle / "main"), mode, str(cell + 3)] if variant == "main" else \
                [str(bundle / "current"), mode, str(cell), variant]
            command = ["/usr/bin/taskset", "-c", "26", "/usr/bin/prlimit", "--cpu=30:30",
                       "--fsize=1048576:1048576", "--"] + argv
            before = sibling_ticks(90)
            stdout, stderr = output / (label + ".stdout"), output / (label + ".stderr")
            with stdout.open("xb") as out, stderr.open("xb") as err:
                result = subprocess.run(command, stdout=out, stderr=err, env=env,
                                        timeout=60, check=False)
            delta = sibling_ticks(90) - before
            require(result.returncode == 0 and stderr.stat().st_size == 0, "child: " + label)
            require(stdout.stat().st_size <= 1048576, "oversized record")
            record = json.loads(stdout.read_text())
            verify()
            return record, delta

        condition("condition-before")
        for cell in range(2):
            for variant in ("main", "auto", "avx2", "gfni"):
                record, _ = invoke(variant, cell, f"check-{cell}-{variant}", False)
                validate(record, expected[variant][cell], False)
                state["preflight"].append(record)
        check_passive(state, plan)
        for cell in range(2):
            for round_id in range(3):
                for comparison, order in ORDERS.items():
                    for slot, variant in enumerate(order):
                        route = "gfni" if variant.startswith("gfni") else variant
                        record, delta = invoke(route, cell,
                            f"cell-{cell}-round-{round_id}-{comparison}-slot-{slot}", True)
                        state["invocations"].append({"cell": cell, "round": round_id,
                            "comparison": comparison, "slot": slot, "variant": variant,
                            "sibling_delta": delta, "record": record})
                        validate(record, expected[route][cell], True)
                        require(delta == 0, "sibling activity; attempt stopped")
                print(f"cell {cell} round {round_id}: AUTO, Leopard1 and same-GFNI control retained",
                      flush=True)
        condition("condition-after")
        require((output / "condition-before.stdout").read_bytes() ==
                (output / "condition-after.stdout").read_bytes(), "condition changed")
        verify()
        state["analysis"] = analyze(state["invocations"])
        state["complete"] = True
    except Exception as error:
        state["failure"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        (output / "attempt.json").write_text(json.dumps(state, indent=2) + "\n")
        for fd in reversed(locks):
            os.close(fd)


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: run_gfni_boundary_screen.py frozen output")
    run(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
