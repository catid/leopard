#!/usr/bin/env python3
"""One .38.5.4.14 fused-ON/overlay-OFF filter; not a retry of the .10 L1 plan."""
import fcntl
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

from run_current_route_screen import CELLS, validate
from run_split_cache_screen import (check_passive, digest, host_identity,
                                    require, sibling_ticks)

ORDERS = {"off_vs_fused": ["off", "on", "on", "off"],
          "same_off": ["off_a", "off_b", "off_b", "off_a"]}
PLAN_NAME = "gfni_terminal_screen_plan.json"
ARTIFACTS = {"timing", "candidate.a", "expected.json", "gfni_terminal_timing.cpp",
             "gfni_terminal.h", "gfni_terminal.cpp", "gfni_terminal.patch",
             "LeopardFF16.cpp", "gfni_terminal_backend.cpp", "Leopard2BackendGFNI.cpp",
             "Leopard2BackendAVX2.cpp", "current_route_screen.cpp", "run_gfni_terminal_screen.py",
             "run_current_route_screen.py", "run_split_cache_screen.py"}


def validate_trace(trace, cell, enabled, measured, exercise=False):
    encodes = 26 if measured or exercise else 1
    passes = (2, 2, 1, 1, 2, 1)[cell]
    matches = encodes * 2 if cell == 0 else 0
    expected = {"schema": "gfni-terminal-timing/v1", "cell": cell,
                "enabled": enabled, "encodes": encodes, "calls": encodes * passes,
                "matches": matches, "changed": matches if enabled else 0,
                "timed": measured, "exercise": exercise}
    require(json.dumps(trace, sort_keys=True) == json.dumps(expected, sort_keys=True),
            "terminal pass trace differs")


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
                    require((row["cell"], row["round"], row["comparison"], row["slot"],
                             row["variant"], row["sibling_delta"]) ==
                            (cell, round_id, comparison, slot, variant, 0), "order/isolation")
                    validate_trace(row["trace"], cell, variant == "on", True)
                    values = row["record"]["samples_ns"]
                    require(type(values) is list and len(values) == 21 and
                            all(type(x) is int and x > 0 for x in values), "analysis samples")
                    medians.append(statistics.median(values))
                ratios[comparison].append(math.sqrt(
                    (medians[0] / medians[1]) * (medians[3] / medians[2])))
        cells.append({"cell": cell, "round_ratios": ratios,
            "ratios": {name: math.exp(statistics.mean(math.log(x) for x in values))
                       for name, values in ratios.items()}})
    controls_ok = all(1 / 1.02 <= cell["ratios"]["same_off"] <= 1.02 for cell in cells)
    # Every non-target is mechanically unchanged; demand equivalence, not
    # merely no slowdown, so implausible off/on gains also invalidate the filter.
    neighbors_ok = all(1 / 1.02 <= cell["ratios"]["off_vs_fused"] <= 1.02
                       for cell in cells[1:])
    target_ok = (cells[0]["ratios"]["off_vs_fused"] >= 1.05 and
                 min(cells[0]["round_ratios"]["off_vs_fused"]) > 1)
    decision = ("inconclusive_controls" if not controls_ok or not neighbors_ok else
                "continue_to_future_qualification" if target_ok else "reject_for_this_screen")
    return {"decision": decision, "cells": cells, "confidence_intervals": False,
            "production_promotion": False, "exact_leopard1_claim": False,
            "authoritative_v19": False}


def validate_plan(plan):
    require(plan["schema"] == "gfni-terminal-screen-plan/v1" and
            plan["bead"] == "leopard-79h.38.5.4.14", "plan identity")
    require(set(plan["artifact_sha256"]) == ARTIFACTS and
            all(type(sha) is str and len(sha) == 64 and
                all(char in "0123456789abcdef" for char in sha)
                for sha in plan["artifact_sha256"].values()), "artifact inventory")
    require(all(type(plan[key]) is int for key in ("cpu", "sibling", "controller_cpu",
                "passive_seconds", "attempt_budget", "rounds", "samples_per_process")),
            "integer protocol fields")
    require((plan["cpu"], plan["sibling"], plan["controller_cpu"], plan["passive_seconds"],
             plan["attempt_budget"], plan["rounds"], plan["samples_per_process"]) ==
            (22, 86, 0, 10, 1, 3, 21), "changed protocol")
    require(plan["orders"] == ORDERS and
            [(x["k"], x["r"], x["bytes"], x["current_route"]) for x in plan["cells"]]
                == CELLS and [x["id"] for x in plan["cells"]] == list(range(6)),
            "changed cases/order")
    require(plan["host"] == {"hostname": "foureyes", "kernel": "6.8.0-138-generic",
        "vendor_id": "AuthenticAMD", "cpu family": "26", "model": "8",
        "model name": "AMD Ryzen Threadripper PRO 9985WX 64-Cores"}, "changed host")


def run(bundle, output):
    output.mkdir(mode=0o700)  # No resume, overwrite, retries or partial pooling.
    plan = json.loads((bundle / PLAN_NAME).read_text())
    pins = json.loads((bundle / "pins.json").read_text())
    expected = json.loads((bundle / "expected.json").read_text())
    state = {"schema": "gfni-terminal-screen-attempt/v1", "pins": pins,
             "plan_sha256": digest(bundle / PLAN_NAME), "host": host_identity(),
             "preflight": [], "invocations": [], "complete": False}
    lock_fds = []
    executable_identity = None

    def verify():
        nonlocal executable_identity
        require(pins["files"] == dict(plan["artifact_sha256"],
                **{PLAN_NAME: state["plan_sha256"]}), "plan/pin inventory differs")
        for name, sha in pins["files"].items():
            require(Path(name).name == name, "flat frozen input name")
            path = bundle / name
            require(path.is_file() and not path.is_symlink() and
                    not path.stat().st_mode & 0o222, "unfrozen input")
            require(digest(path) == sha, "input changed: " + name)
        info = (bundle / "timing").stat()
        identity = [info.st_dev, info.st_ino, info.st_size, info.st_mode]
        if executable_identity is None:
            executable_identity = identity
            state["executable_identity"] = identity
        require(identity == executable_identity, "executable identity changed")

    try:
        validate_plan(plan)
        require(plan["host"] == state["host"], "host changed")
        for cpu in (22, 86):
            require(Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
                    .read_text().strip() == "22,86", "physical topology changed")
        os.sched_setaffinity(0, {0})
        root = Path(f"/run/user/{os.getuid()}/leopard2-cpu-leases")
        require(root.is_dir() and not root.is_symlink() and root.stat().st_uid == os.getuid()
                and root.stat().st_mode & 0o777 == 0o700, "unsafe lease directory")
        for path in (Path("/tmp/leopard-gf8-authoritative.lock"),
                     root / f"leopard2-cpu-pair-{os.getuid()}-22-86.lock"):
            fd = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_NOFOLLOW, 0o600)
            lock_fds.append(fd)
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        verify()
        env = {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
               "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "1"}

        def invoke(cell, enabled, label, measured):
            verify()
            command = ["/usr/bin/taskset", "-c", "22", "/usr/bin/prlimit", "--cpu=30:30",
                       "--fsize=1048576:1048576", "--", str(bundle / "timing"),
                       "--measure" if measured else "--check", str(cell),
                       "--terminal=1" if enabled else "--terminal=0"]
            before = sibling_ticks(86)
            stdout, stderr = output / (label + ".stdout"), output / (label + ".stderr")
            with stdout.open("xb") as out, stderr.open("xb") as err:
                child = subprocess.run(command, stdout=out, stderr=err, env=env,
                                       timeout=60, check=False)
            after = sibling_ticks(86)
            require(child.returncode == 0, "child failed: " + label)
            require(stdout.stat().st_size <= 1048576 and stderr.stat().st_size <= 1048576,
                    "output size")
            record, trace = json.loads(stdout.read_text()), json.loads(stderr.read_text())
            verify()
            validate(record, expected[cell], measured)
            validate_trace(trace, cell, enabled, measured)
            return record, trace, after - before

        for cell in range(6):
            for enabled in (False, True):
                record, trace, _ = invoke(cell, enabled, f"check-{cell}-{int(enabled)}", False)
                state["preflight"].append({"record": record, "trace": trace})
        check_passive(state, plan)
        for cell in range(6):
            for round_id in range(3):
                for comparison, order in ORDERS.items():
                    for slot, variant in enumerate(order):
                        record, trace, delta = invoke(cell, variant == "on",
                            f"cell-{cell}-round-{round_id}-{comparison}-slot-{slot}", True)
                        state["invocations"].append({"cell": cell, "round": round_id,
                            "comparison": comparison, "slot": slot, "variant": variant,
                            "sibling_delta": delta, "record": record, "trace": trace})
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
    require(len(sys.argv) == 3, "usage: run_gfni_terminal_screen.py frozen output")
    run(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
