#!/usr/bin/env python3
"""Single local native/pure-L1/L2 AVX2 attribution attempt, not promotion."""
import fcntl
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

from run_split_cache_screen import check_passive, digest, host_identity, require, sibling_ticks

PLAN = "avx2_isa_screen_plan.json"
PROFILES = ("native", "pure", "current")
ORDERS = {
    "pure_over_native": ["native", "pure", "pure", "native"],
    "current_over_native": ["native", "current", "current", "native"],
    "current_over_pure": ["pure", "current", "current", "pure"],
    "same_native": ["native"] * 4,
    "same_pure": ["pure"] * 4,
    "same_current": ["current"] * 4,
}
HOST = {"hostname": "work", "kernel": "6.8.0-137-generic", "vendor_id": "AuthenticAMD",
        "cpu family": "26", "model": "8", "model name": "AMD Ryzen Threadripper 9980X 64-Cores"}
PROTOCOL = {
    "schema": "leopard-avx2-isa-plan/v1", "bead": "leopard-79h.38.5.4.18.1",
    "host": HOST, "cpu": 26, "sibling": 90, "controller_cpu": 0,
    "condition": "slipgate-disabled-20260909", "attempt_budget": 1,
    "rounds": 3, "samples_per_process": 21, "passive_seconds": 10,
    "orders": ORDERS, "cells": [
        {"k": 1000, "r": 200, "bytes": 65536},
        {"k": 1000, "r": 200, "bytes": 32768},
        {"k": 1000, "r": 199, "bytes": 65536},
        {"k": 4096, "r": 512, "bytes": 4096}],
    "control_bound": 1.02, "future_candidate_minimum_gain": 1.05,
    "production_promotion": False, "confidence_intervals": False,
    "attribution": "native target/tuning versus generic AVX2 ceiling; not single-instruction ablation",
}
ARTIFACTS = {*PROFILES, *(x + ".a" for x in PROFILES), "expected.json"}
FROZEN_FILES = ARTIFACTS | {PLAN, "avx2_isa_screen.cpp", "run_avx2_isa_screen.py",
    "run_split_cache_screen.py", "test_avx2_isa_screen.py", "replay_avx2_isa_screen.py",
    "test_replay_avx2_isa_screen.py", "check-condition.sh", "check.sh", "build.json", "preflight.json"}


def same(actual, expected):
    require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), "value/type mismatch")


def identity(record):
    return {k: v for k, v in record.items() if k not in ("samples_ns", "encode_calls")}


def validate(record, expected, measured, exercise=False):
    same(identity(record), expected)
    same(record["encode_calls"], 26 if measured or exercise else 1)
    samples = record["samples_ns"]
    require(type(samples) is list and len(samples) == (21 if measured else 0) and
            all(type(v) is int and v > 0 for v in samples), "invalid samples")


def validate_plan(plan):
    same(set_as_list(plan), sorted([*PROTOCOL, "artifact_sha256"]))
    same({k: plan[k] for k in PROTOCOL}, PROTOCOL)
    same(sorted(plan["artifact_sha256"]), sorted(ARTIFACTS))
    require(all(type(v) is str and len(v) == 64 and all(c in "0123456789abcdef" for c in v)
                for v in plan["artifact_sha256"].values()), "invalid artifact hashes")


def set_as_list(value):
    require(type(value) is dict, "expected object")
    return sorted(value)


def validate_pins(pins):
    same(set_as_list(pins), ["files", "schema"])
    same(pins["schema"], "leopard-avx2-isa-pins/v1")
    same(set_as_list(pins["files"]), sorted(FROZEN_FILES))
    require(all(type(v) is str and len(v) == 64 and all(c in "0123456789abcdef" for c in v)
                for v in pins["files"].values()), "invalid frozen hashes")


def schedule():
    for cell in range(4):
        for round_id in range(3):
            for comparison, order in ORDERS.items():
                for slot, profile in enumerate(order):
                    yield dict(cell=cell, round=round_id, comparison=comparison,
                               slot=slot, profile=profile)


def analyze(rows, expected):
    require(type(rows) is list and len(rows) == 288, "partial attempt has no analysis")
    ratios = [{key: [] for key in ORDERS} for _ in range(4)]
    medians = []
    for row, wanted in zip(rows, schedule()):
        same({k: row[k] for k in wanted}, wanted)
        same(row["sibling_delta"], 0)
        validate(row["record"], expected[wanted["profile"]][wanted["cell"]], True)
        medians.append(statistics.median(row["record"]["samples_ns"]))
        if wanted["slot"] == 3:
            ratios[wanted["cell"]][wanted["comparison"]].append(
                math.sqrt((medians[0] / medians[1]) * (medians[3] / medians[2])))
            medians = []
    cells = [dict(cell=i, round_ratios=values, ratios={k: math.exp(statistics.mean(
                  math.log(v) for v in data)) for k, data in values.items()})
             for i, values in enumerate(ratios)]
    controls = all(1 / 1.02 <= cell["ratios"]["same_" + p] <= 1.02
                   for cell in cells for p in PROFILES)
    # Same-current at cell3 is the unchanged neighbor control. Native/pure
    # differences in that cell are measured attribution effects, not regressions
    # in an unmodified L2 candidate. No candidate promotion is decided here.
    return dict(controls_pass=controls, cells=cells,
                decision="attribution_available" if controls else "inconclusive_controls",
                production_promotion=False, confidence_intervals=False, authoritative_v19=False)


def run(bundle, output, preregistration):
    require(len(preregistration) == 40 and all(c in "0123456789abcdef" for c in preregistration),
            "preregistration commit")
    require(output == bundle.parent / "attempt1", "fixed single-attempt location")
    output.mkdir(mode=0o700)
    plan = json.loads((bundle / PLAN).read_text())
    pins = json.loads((bundle / "pins.json").read_text())
    expected = json.loads((bundle / "expected.json").read_text())
    validate_plan(plan)
    validate_pins(pins)
    repo_path = "experiments/leopard2/gf16_high_encode/"
    require(subprocess.check_output(["git", "show", preregistration + ":" + repo_path + PLAN]) ==
            (bundle / PLAN).read_bytes(), "plan differs from preregistration")
    subprocess.run(["git", "merge-base", "--is-ancestor", preregistration,
                    "origin/codex/claude-fable-5-1-audit"], check=True)
    for name in ("run_avx2_isa_screen.py", "run_split_cache_screen.py",
                 "avx2_isa_screen.cpp", "test_avx2_isa_screen.py"):
        require(subprocess.check_output(["git", "show", preregistration + ":" + repo_path + name]) ==
                (bundle / name).read_bytes(), "source differs from preregistration")
    same(host_identity(), HOST)
    for cpu in (26, 90):
        require(Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
                .read_text().strip() == "26,90", "topology")
    os.sched_setaffinity(0, {0})
    state = dict(schema="leopard-avx2-isa-attempt/v1", preregistration=preregistration,
                 pins=pins, plan_sha256=digest(bundle / PLAN), host=host_identity(),
                 preflight=[], invocations=[], complete=False)
    locks = []

    def verify():
        for name, value in pins["files"].items():
            require(Path(name).name == name, "unsafe pin name")
            path = bundle / name
            require(path.is_file() and not path.is_symlink() and
                    not path.stat().st_mode & 0o222 and digest(path) == value, "input drift: " + name)
        for name, value in plan["artifact_sha256"].items():
            same(pins["files"][name], value)

    try:
        lease = Path(f"/run/user/{os.getuid()}/leopard2-cpu-leases")
        require(lease.is_dir() and not lease.is_symlink() and lease.stat().st_uid == os.getuid()
                and lease.stat().st_mode & 0o777 == 0o700, "unsafe lease directory")
        for path in (Path("/tmp/leopard-gf8-authoritative.lock"),
                     lease / f"leopard2-cpu-pair-{os.getuid()}-26-90.lock"):
            fd = os.open(path, os.O_RDONLY | os.O_CREAT | os.O_NOFOLLOW, 0o600)
            locks.append(fd)
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        verify()
        env = dict(PATH="/usr/bin:/bin", LANG="C", LC_ALL="C",
                   OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", OMP_THREAD_LIMIT="1")

        def condition(label):
            with (output / (label + ".stdout")).open("xb") as out, \
                    (output / (label + ".stderr")).open("xb") as err:
                result = subprocess.run(["/bin/bash", str(bundle / "check-condition.sh")],
                                        stdout=out, stderr=err, timeout=30)
            require(result.returncode == 0 and (output / (label + ".stderr")).stat().st_size == 0,
                    "condition: " + label)

        def invoke(profile, cell, label, measured):
            command = ["/usr/bin/taskset", "-c", "26", "/usr/bin/prlimit", "--cpu=30:30",
                       "--fsize=1048576:1048576", "--", str(bundle / profile),
                       "--measure" if measured else "--check", str(cell)]
            before = sibling_ticks(90)
            stdout, stderr = output / (label + ".stdout"), output / (label + ".stderr")
            with stdout.open("xb") as out, stderr.open("xb") as err:
                result = subprocess.run(command, stdout=out, stderr=err, env=env, timeout=60)
            delta = sibling_ticks(90) - before
            require(result.returncode == 0 and stderr.stat().st_size == 0, "child: " + label)
            require(stdout.stat().st_size <= 1048576, "oversized output")
            record = json.loads(stdout.read_text())
            verify()
            return record, delta

        condition("condition-before")
        for cell in range(4):
            for profile in PROFILES:
                record, _ = invoke(profile, cell, f"check-{cell}-{profile}", False)
                validate(record, expected[profile][cell], False)
                state["preflight"].append(record)
        check_passive(state, plan)
        for item in schedule():
            cell, profile = item["cell"], item["profile"]
            label = f"cell-{cell}-round-{item['round']}-{item['comparison']}-slot-{item['slot']}"
            record, delta = invoke(profile, cell, label, True)
            state["invocations"].append(dict(item, record=record, sibling_delta=delta))
            validate(record, expected[profile][cell], True)
            same(delta, 0)
            if item["slot"] == 3 and item["comparison"] == "same_current":
                print(f"cell {cell} round {item['round']}: three comparisons and three controls retained",
                      flush=True)
        condition("condition-after")
        same((output / "condition-before.stdout").read_text(),
             (output / "condition-after.stdout").read_text())
        verify()
        state["analysis"] = analyze(state["invocations"], expected)
        state["complete"] = True
    except Exception as error:
        state["failure"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        (output / "attempt.json").write_text(json.dumps(state, indent=2) + "\n")
        for fd in reversed(locks):
            os.close(fd)


if __name__ == "__main__":
    require(len(sys.argv) == 4, "usage: run_avx2_isa_screen.py frozen attempt1 pushed_commit")
    run(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), sys.argv[3])
