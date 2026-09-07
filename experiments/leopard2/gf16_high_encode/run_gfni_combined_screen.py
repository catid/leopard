#!/usr/bin/env python3
"""One bounded four-mode .38.5.4.16 screen; no retry of earlier experiments."""
import fcntl
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

from run_current_route_screen import CELLS, validate
from run_split_cache_screen import check_passive, digest, host_identity, require, sibling_ticks

PLAN_NAME = "gfni_combined_screen_plan.json"
ARTIFACTS = {"timing", "candidate.a", "expected.json", "gfni_combined_timing.cpp",
    "gfni_combined.h", "gfni_combined.cpp", "gfni_combined.patch", "LeopardFF16.cpp",
    "gfni_terminal_backend.cpp", "Leopard2BackendGFNI.cpp", "Leopard2BackendAVX2.cpp",
    "current_route_screen.cpp", "run_gfni_combined_screen.py",
    "run_current_route_screen.py", "run_split_cache_screen.py"}
ROUND_ORDERS = ((0, 1, 2, 3, 3, 2, 1, 0),
                (1, 2, 3, 0, 0, 3, 2, 1),
                (2, 3, 0, 1, 1, 0, 3, 2))
COMPARISONS = ("factorial", "same_off")


def equal(actual, expected):
    require(json.dumps(actual, sort_keys=True) == json.dumps(expected, sort_keys=True), "record mismatch")


def schedule():
    for cell in range(6):
        for round_id, order in enumerate(ROUND_ORDERS):
            for comparison in COMPARISONS:
                for slot, variant in enumerate(order):
                    yield dict(cell=cell, round=round_id, comparison=comparison, slot=slot,
                        variant=variant, mode=variant if comparison == "factorial" else 0)


def validate_trace(trace, cell, mode, measured, exercise=False):
    require(type(cell) is int and 0 <= cell < 6 and type(mode) is int and 0 <= mode < 4 and
            type(measured) is bool and type(exercise) is bool and not (measured and exercise), "trace scope")
    encodes = 26 if measured or exercise else 1
    matches = encodes * 2 if cell == 0 else 0
    equal(trace, dict(schema="gfni-combined-timing/v1", cell=cell, mode=mode,
        encodes=encodes, calls=encodes * (2, 2, 1, 1, 2, 1)[cell], matches=matches,
        first=matches if mode & 1 else 0, terminal=matches if mode & 2 else 0,
        timed=measured, exercise=exercise))


def analyze(rows):
    require(type(rows) is list and len(rows) == 288, "partial attempts have no analysis")
    times = {}
    for row, position in zip(rows, schedule()):
        equal({key:row[key] for key in position}, position)
        require(type(row["sibling_delta"]) is int and row["sibling_delta"] == 0, "sibling isolation")
        validate_trace(row["trace"], position["cell"], position["mode"], True)
        values = row["record"]["samples_ns"]
        require(type(values) is list and len(values) == 21 and
                all(type(value) is int and value > 0 for value in values), "analysis samples")
        key = (row["cell"], row["round"], row["comparison"], row["variant"])
        times.setdefault(key, []).append(statistics.median(values))
    cells = []
    for cell in range(6):
        rounds = {comparison: {str(mode):[] for mode in (1, 2, 3)} for comparison in COMPARISONS}
        interactions = []
        for round_id in range(3):
            for comparison in COMPARISONS:
                means = []
                for mode in range(4):
                    values = times[cell, round_id, comparison, mode]
                    require(len(values) == 2, "two mirrored observations per mode")
                    means.append(math.sqrt(values[0] * values[1]))
                for mode in (1, 2, 3):
                    rounds[comparison][str(mode)].append(means[0] / means[mode])
                if comparison == "factorial":
                    # >1 means the combined speedup exceeds the product of
                    # the two single-mode speedups within this same fresh round.
                    interactions.append(means[1] * means[2] / (means[0] * means[3]))
        geometric = lambda values: math.exp(statistics.mean(math.log(value) for value in values))
        cells.append(dict(cell=cell, round_ratios=rounds,
            ratios={comparison:{mode:geometric(values) for mode, values in modes.items()}
                    for comparison, modes in rounds.items()},
            interaction_rounds=interactions, interaction_factor=geometric(interactions)))
    controls = [value for cell in cells for value in cell["ratios"]["same_off"].values()]
    controls += [value for cell in cells[1:] for value in cell["ratios"]["factorial"].values()]
    controls_ok = all(1 / 1.02 <= value <= 1.02 for value in controls)
    target_ok = (cells[0]["ratios"]["factorial"]["3"] >= 1.05 and
                 min(cells[0]["round_ratios"]["factorial"]["3"]) > 1)
    decision = ("inconclusive_controls" if not controls_ok else
                "continue_to_future_qualification" if target_ok else "reject_for_this_screen")
    return dict(decision=decision, cells=cells, aggregate_controls=len(controls),
        confidence_intervals=False, production_promotion=False,
        exact_leopard1_claim=False, authoritative_v19=False)


def validate_plan(plan):
    require(plan["schema"] == "gfni-combined-screen-plan/v1" and
            plan["bead"] == "leopard-79h.38.5.4.16", "plan identity")
    require(set(plan["artifact_sha256"]) == ARTIFACTS and
            all(type(sha) is str and len(sha) == 64 and all(c in "0123456789abcdef" for c in sha)
                for sha in plan["artifact_sha256"].values()), "artifact inventory")
    fields = ("cpu", "sibling", "controller_cpu", "passive_seconds", "attempt_budget",
              "rounds", "samples_per_process", "timed_invocations", "seed",
              "initial_untimed_encodes", "additional_warmups")
    equal([plan[key] for key in fields], [22, 86, 0, 10, 1, 3, 21, 288, 20260906, 1, 4])
    equal(plan["round_orders"], ROUND_ORDERS)
    equal(plan["comparisons"], COMPARISONS)
    equal(plan["thresholds"], dict(target_ratio=1.05, control_factor=1.02, positive_target_rounds=3))
    equal([[x["k"], x["r"], x["bytes"], x["current_route"]] for x in plan["cells"]], CELLS)
    equal([x["id"] for x in plan["cells"]], list(range(6)))
    equal(plan["host"], {"hostname":"foureyes", "kernel":"6.8.0-138-generic",
        "vendor_id":"AuthenticAMD", "cpu family":"26", "model":"8",
        "model name":"AMD Ryzen Threadripper PRO 9985WX 64-Cores"})


def run(bundle, output):
    output.mkdir(mode=0o700)  # One attempt; no resume, overwrite or pooling.
    plan = json.loads((bundle / PLAN_NAME).read_text())
    pins = json.loads((bundle / "pins.json").read_text())
    expected = json.loads((bundle / "expected.json").read_text())
    state = dict(schema="gfni-combined-screen-attempt/v1", pins=pins,
        plan_sha256=digest(bundle / PLAN_NAME), host=host_identity(),
        preflight=[], invocations=[], complete=False)
    lock_fds = []
    executable_identity = None

    def verify():
        nonlocal executable_identity
        equal(pins["files"], dict(plan["artifact_sha256"], **{PLAN_NAME:state["plan_sha256"]}))
        for name, sha in pins["files"].items():
            require(Path(name).name == name, "flat frozen input name")
            path = bundle / name
            require(path.is_file() and not path.is_symlink() and not path.stat().st_mode & 0o222,
                    "unfrozen input")
            require(digest(path) == sha, "input changed: " + name)
        info = (bundle / "timing").stat()
        identity = [info.st_dev, info.st_ino, info.st_size, info.st_mode]
        if executable_identity is None:
            executable_identity = identity
            state["executable_identity"] = identity
        equal(identity, executable_identity)

    try:
        validate_plan(plan)
        equal(plan["host"], state["host"])
        for cpu in (22, 86):
            require(Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
                .read_text().strip() == "22,86", "physical topology changed")
        os.sched_setaffinity(0, {0})
        root = Path(f"/run/user/{os.getuid()}/leopard2-cpu-leases")
        require(root.is_dir() and not root.is_symlink() and root.stat().st_uid == os.getuid() and
                root.stat().st_mode & 0o777 == 0o700, "unsafe lease directory")
        for path in (Path("/tmp/leopard-gf8-authoritative.lock"),
                     root / f"leopard2-cpu-pair-{os.getuid()}-22-86.lock"):
            fd = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_NOFOLLOW, 0o600)
            lock_fds.append(fd)
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        verify()
        env = {"PATH":"/usr/bin:/bin", "LANG":"C", "LC_ALL":"C", "OMP_NUM_THREADS":"1",
               "OMP_DYNAMIC":"FALSE", "OMP_THREAD_LIMIT":"1"}

        def invoke(cell, mode, label, measured):
            verify()
            command = ["/usr/bin/taskset", "-c", "22", "/usr/bin/prlimit", "--cpu=30:30",
                "--fsize=1048576:1048576", "--", str(bundle / "timing"),
                "--measure" if measured else "--check", str(cell), "--mode=%d" % mode]
            before = sibling_ticks(86)
            stdout, stderr = output / (label + ".stdout"), output / (label + ".stderr")
            with stdout.open("xb") as out, stderr.open("xb") as err:
                child = subprocess.run(command, stdout=out, stderr=err, env=env, timeout=60, check=False)
            after = sibling_ticks(86)
            require(child.returncode == 0, "child failed: " + label)
            require(stdout.stat().st_size <= 1048576 and stderr.stat().st_size <= 1048576, "output size")
            record, trace = json.loads(stdout.read_text()), json.loads(stderr.read_text())
            verify()
            validate(record, expected[cell], measured)
            validate_trace(trace, cell, mode, measured)
            return record, trace, after - before

        for cell in range(6):
            for mode in range(4):
                record, trace, _ = invoke(cell, mode, f"check-{cell}-{mode}", False)
                state["preflight"].append(dict(record=record, trace=trace))
        check_passive(state, plan)
        for position in schedule():
            cell, round_id, comparison, slot = (position[key] for key in ("cell", "round", "comparison", "slot"))
            label = f"cell-{cell}-round-{round_id}-{comparison}-slot-{slot}"
            record, trace, delta = invoke(cell, position["mode"], label, True)
            state["invocations"].append(dict(position, sibling_delta=delta, record=record, trace=trace))
            require(delta == 0, "sibling activity; attempt stopped")
            if comparison == "same_off" and slot == 7:
                print(f"cell {cell} round {round_id} four modes and matched control retained", flush=True)
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
    require(len(sys.argv) == 3, "usage: run_gfni_combined_screen.py frozen output")
    run(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
