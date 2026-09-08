#!/usr/bin/env python3
"""One fixed current-route diagnostic; never imports the v19 controller."""
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

ORDERS = {"main_vs_current": ["main", "current", "current", "main"],
          "same_current": ["current_a", "current_b", "current_b", "current_a"]}
CELLS = [(1000, 200, 65536, "gfni"), (1000, 200, 65536, "avx2"),
         (1000, 200, 65536, "avx512"), (1000, 200, 32768, "avx2"),
         (1000, 199, 65536, "avx2"), (4096, 512, 4096, "avx2")]


def validate(record, expected, measured):
    require(record.get("schema") == "leopard-gf16-current-route-screen/v1", "schema")
    require(identity(record) == expected, "changed workload, route, source or scratch")
    samples = record.get("samples_ns")
    require(type(samples) is list and len(samples) == (21 if measured else 0), "samples")
    require(all(type(x) is int and x > 0 for x in samples), "sample values")


def analyze(rows):
    require(len(rows) == 144, "partial attempts have no analysis")
    cells = []
    cursor = 0
    for cell in range(6):
        ratios = {name: [] for name in ORDERS}
        for round_id in range(3):
            for comparison, order in ORDERS.items():
                medians = []
                for slot, variant in enumerate(order):
                    row = rows[cursor]
                    cursor += 1
                    require((row["cell"], row["round"], row["comparison"],
                             row["slot"], row["variant"], row["sibling_delta"]) ==
                            (cell, round_id, comparison, slot, variant, 0), "order/isolation")
                    values = row["record"]["samples_ns"]
                    require(len(values) == 21 and all(type(x) is int and x > 0
                                                     for x in values), "analysis samples")
                    medians.append(statistics.median(values))
                ratios[comparison].append(math.sqrt(
                    (medians[0] / medians[1]) * (medians[3] / medians[2])))
        cells.append({"cell": cell, "round_ratios": ratios,
            "ratios": {name: math.exp(statistics.mean(math.log(x) for x in values))
                       for name, values in ratios.items()}})
    controls_ok = all(1 / 1.02 <= cell["ratios"]["same_current"] <= 1.02 for cell in cells)
    for cell in cells:
        ratio = cell["ratios"]["main_vs_current"]
        rounds = cell["round_ratios"]["main_vs_current"]
        cell["interpretation"] = (
            "no_inference_controls" if not controls_ok else
            "investigate_current_deficit" if ratio < 1 / 1.02 and max(rounds) < 1 else
            "current_advantage" if ratio > 1.02 and min(rounds) > 1 else
            "near_parity_or_uncertain")
    return {"decision": "diagnostic_complete" if controls_ok else "inconclusive_controls",
            "cells": cells, "confidence_intervals": False,
            "authoritative_v19": False, "production_promotion": False,
            "historical_exact_main_gap_closed": False}


def validate_plan(plan, name):
    profiles = {
        "current_route_screen_plan.json": (22, 86, "foureyes", "6.8.0-138-generic",
            "AMD Ryzen Threadripper PRO 9985WX 64-Cores"),
        "current_route_screen_work_plan.json": (26, 90, "work", "6.8.0-137-generic",
            "AMD Ryzen Threadripper 9980X 64-Cores"),
    }
    require(name in profiles, "unsupported plan name")
    cpu, sibling, hostname, kernel, model = profiles[name]
    for key in ("cpu", "sibling", "controller_cpu", "passive_seconds",
                "attempt_budget", "rounds", "samples_per_process"):
        require(type(plan[key]) is int, "protocol integer: " + key)
    require((plan["cpu"], plan["sibling"], plan["controller_cpu"], plan["passive_seconds"],
             plan["attempt_budget"], plan["rounds"], plan["samples_per_process"]) ==
            (cpu, sibling, 0, 10, 1, 3, 21), "changed protocol")
    require(plan["orders"] == ORDERS and
            [(x["k"], x["r"], x["bytes"], x["current_route"]) for x in plan["cells"]]
                == CELLS, "changed cases/order")
    require(plan["host"] == {"hostname": hostname, "kernel": kernel,
        "vendor_id": "AuthenticAMD", "cpu family": "26", "model": "8",
        "model name": model}, "changed host profile")


def run(bundle, output, plan_name="current_route_screen_plan.json"):
    validate_name = plan_name in ("current_route_screen_plan.json",
                                 "current_route_screen_work_plan.json")
    require(validate_name, "unsupported plan name")
    output.mkdir(mode=0o700)  # No resume, overwrite, retry or partial pooling.
    plan = json.loads((bundle / plan_name).read_text())
    pins = json.loads((bundle / "pins.json").read_text())
    expected = json.loads((bundle / "expected.json").read_text())
    validate_plan(plan, plan_name)
    cpu, sibling = plan["cpu"], plan["sibling"]
    require(plan["host"] == host_identity(), "host changed")
    for member in (cpu, sibling):
        require(Path(f"/sys/devices/system/cpu/cpu{member}/topology/thread_siblings_list")
                .read_text().strip() == f"{cpu},{sibling}", "physical topology changed")
    os.sched_setaffinity(0, {0})
    state = {"schema": "leopard-gf16-current-route-attempt/v1", "pins": pins,
             "plan_sha256": digest(bundle / plan_name),
             "host": host_identity(), "preflight": [], "invocations": [], "complete": False}
    lock_fds = []

    def verify():
        for name, sha in pins["files"].items():
            path = bundle / name
            require(path.is_file() and not path.is_symlink() and
                    not path.stat().st_mode & 0o222, "unfrozen input")
            require(digest(path) == sha, "input changed: " + name)
        for name in ("main", "current", "main.a", "current.a", "expected.json"):
            require(pins["files"][name] == plan["artifact_sha256"][name], "plan artifact")

    try:
        root = Path(f"/run/user/{os.getuid()}/leopard2-cpu-leases")
        require(root.is_dir() and not root.is_symlink() and root.stat().st_uid == os.getuid()
                and root.stat().st_mode & 0o777 == 0o700, "unsafe lease directory")
        for path in (Path("/tmp/leopard-gf8-authoritative.lock"),
                     root / f"leopard2-cpu-pair-{os.getuid()}-{cpu}-{sibling}.lock"):
            fd = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_NOFOLLOW, 0o600)
            lock_fds.append(fd)
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        verify()
        env = {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
               "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "1"}

        def invoke(executable, cell, label, measured):
            command = ["/usr/bin/taskset", "-c", str(cpu), "/usr/bin/prlimit", "--cpu=30:30",
                       "--fsize=1048576:1048576", "--", str(bundle / executable),
                       "--measure" if measured else "--check", str(cell)]
            before = sibling_ticks(sibling)
            stdout, stderr = output / (label + ".stdout"), output / (label + ".stderr")
            with stdout.open("xb") as out, stderr.open("xb") as err:
                child = subprocess.run(command, stdout=out, stderr=err, env=env,
                                       timeout=60, check=False)
            after = sibling_ticks(sibling)
            require(child.returncode == 0 and stderr.stat().st_size == 0, "child: " + label)
            require(stdout.stat().st_size <= 1048576, "output size")
            record = json.loads(stdout.read_text())
            verify()
            return record, after - before

        for cell in range(6):
            for variant in ("main", "current"):
                record, _ = invoke(variant, cell, f"check-{cell}-{variant}", False)
                validate(record, expected[variant][cell], False)
                state["preflight"].append(record)
        check_passive(state, plan)
        for cell in range(6):
            for round_id in range(3):
                for comparison, order in ORDERS.items():
                    for slot, variant in enumerate(order):
                        executable = "main" if variant == "main" else "current"
                        record, delta = invoke(executable, cell,
                            f"cell-{cell}-round-{round_id}-{comparison}-slot-{slot}", True)
                        state["invocations"].append({"cell": cell, "round": round_id,
                            "comparison": comparison, "slot": slot, "variant": variant,
                            "sibling_delta": delta, "record": record})
                        validate(record, expected[executable][cell], True)
                        require(delta == 0, "sibling activity; attempt stopped")
                print(f"cell {cell} round {round_id} and same-binary control retained", flush=True)
        verify()
        state["analysis"] = analyze(state["invocations"])
        state["complete"] = True
    except Exception as error:
        state["failure"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        (output / "attempt.json").write_text(json.dumps(state, indent=2) + "\n")
        for fd in reversed(lock_fds):
            os.close(fd)


if __name__ == "__main__":
    require(len(sys.argv) in (3, 4),
            "usage: run_current_route_screen.py frozen output [plan-name]")
    run(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(),
        sys.argv[3] if len(sys.argv) == 4 else "current_route_screen_plan.json")
