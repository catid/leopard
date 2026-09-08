#!/usr/bin/env python3
"""Offline validation of the 2026-09-08 local failure; executes no codec.

This standalone replay intentionally produces no ratios from the partial run.
It imports neither the collector nor production code. Historical correctness
is separately replayed from the unchanged referenced archive.
"""
import copy
import hashlib
import json
from pathlib import Path
import sys


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(65536), b""):
            result.update(block)
    return result.hexdigest()


def read_json(path):
    require(path.stat().st_size < 1048576, "JSON bound")
    return json.loads(path.read_text())


def validate_attempt(state, expected, pins, host):
    require(state["schema"] == "leopard-gf16-current-route-attempt/v1", "schema")
    require(state["pins"] == pins and state["host"] == host, "pins/host")
    require(state["plan_sha256"] == pins["files"]["current_route_screen_work_plan.json"],
            "plan binding")
    require(state["complete"] is False and "analysis" not in state, "partial inference")
    require(state["failure"] == "ValueError: sibling activity; attempt stopped", "failure")
    require(state["passive"] == {"before": 569301, "after": 569301,
                                "elapsed_ns": 10000076453}, "passive gate")
    require(len(state["preflight"]) == 12, "check count")
    for cell in range(6):
        for index, variant in enumerate(("main", "current")):
            require(state["preflight"][cell * 2 + index] ==
                    dict(expected[variant][cell], samples_ns=[]), "check identity")
    rows = state["invocations"]
    require(len(rows) == 24, "partial invocation count")
    orders = (("main_vs_current", ("main", "current", "current", "main")),
              ("same_current", ("current_a", "current_b", "current_b", "current_a")))
    for index, row in enumerate(rows):
        comparison, order = orders[(index % 8) // 4]
        variant = order[index % 4]
        require((row["cell"], row["round"], row["comparison"], row["slot"], row["variant"])
                == (0, index // 8, comparison, index % 4, variant), "invocation order")
        require(type(row["sibling_delta"]) is int and
                row["sibling_delta"] == (4 if index == 23 else 0), "isolation failure")
        record = row["record"]
        samples = record["samples_ns"]
        require(type(samples) is list and len(samples) == 21 and
                all(type(x) is int and x > 0 for x in samples), "samples")
        identity = {key: value for key, value in record.items() if key != "samples_ns"}
        require(identity == expected["main" if variant == "main" else "current"][0],
                "measured identity")


def replay(root):
    fixed = {
        "attempt1/attempt.json": "c9885e8e70f5ce71aca405d253961fe0c9bf1a1b202b74b946ff237365126e2f",
        "attempt1.log": "67644df0e8eccc1a43908bdb72fc05b793ee57ac03a058996e8dab2fee7ad523",
        "frozen/pins.json": "eb1ada8ff7665924884a350b792ab41e518f1baf437d0ff5489bfbea65bbf316",
    }
    for name, digest in fixed.items():
        require(sha(root / name) == digest, "fixed evidence: " + name)
    pins = read_json(root / "frozen/pins.json")
    require(pins["source_commit"] == "b8b9b1a9d8074b14f75fcf08a88646eefa8f8c2a" and
            len(pins["files"]) == 14, "preregistration/inventory")
    for name, digest in pins["files"].items():
        path = root / "frozen" / name
        require(Path(name).name == name and not path.is_symlink() and
                not path.stat().st_mode & 0o222 and sha(path) == digest, "frozen: " + name)
    plan = read_json(root / "frozen/current_route_screen_work_plan.json")
    require((plan["cpu"], plan["sibling"], plan["controller_cpu"], plan["attempt_budget"])
            == (26, 90, 0, 1), "local protocol")
    require(plan["host"]["hostname"] == "work", "local host")
    for name, digest in plan["artifact_sha256"].items():
        require(pins["files"][name] == digest, "plan artifact")
    state = read_json(root / "attempt1/attempt.json")
    expected = read_json(root / "frozen/expected.json")
    validate_attempt(state, expected, pins, plan["host"])
    for cell in range(6):
        for index, variant in enumerate(("main", "current")):
            prefix = root / "attempt1" / f"check-{cell}-{variant}"
            require(read_json(prefix.with_suffix(".stdout")) ==
                    state["preflight"][cell * 2 + index], "raw check")
            require(prefix.with_suffix(".stderr").stat().st_size == 0, "check stderr")
    for row in state["invocations"]:
        prefix = root / "attempt1" / (
            f"cell-{row['cell']}-round-{row['round']}-{row['comparison']}-slot-{row['slot']}")
        require(read_json(prefix.with_suffix(".stdout")) == row["record"], "raw timing")
        require(prefix.with_suffix(".stderr").stat().st_size == 0, "timing stderr")
    log = (root / "attempt1.log").read_text()
    require("memory.peak\n132247552\nmemory.max\n268435456\n" in log and
            "memory.events\nlow 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\n"
            "memory.swap.current\n0\nmemory.swap.max\n0\n" in log, "resource evidence")

    # Adversarial semantic checks bypass fixed hashes deliberately; no files change.
    rejected = 0
    for mutation in range(11):
        changed = copy.deepcopy(state)
        if mutation == 0:
            changed["complete"] = True
        elif mutation == 1:
            changed["analysis"] = {}
        elif mutation == 2:
            changed["passive"]["after"] += 1
        elif mutation == 3:
            changed["invocations"].pop()
        elif mutation == 4:
            changed["invocations"][-1]["sibling_delta"] = 0
        elif mutation == 5:
            changed["invocations"][0]["sibling_delta"] = 1
        elif mutation == 6:
            changed["invocations"][0]["record"]["samples_ns"][0] = True
        elif mutation == 7:
            changed["invocations"][0]["slot"] = 1
        elif mutation == 8:
            changed["preflight"].pop()
        elif mutation == 9:
            changed["host"]["hostname"] = "foureyes"
        else:
            changed["plan_sha256"] = "0" * 64
        try:
            validate_attempt(changed, expected, pins, plan["host"])
        except ValueError:
            rejected += 1
    require(rejected == 11, "mutation escaped")
    return {"checks": 12, "timed_invocations": 24, "failure_invocation": 24,
            "sibling_nonidle_jiffies": 4, "complete": False,
            "performance_inference": False, "rejected_mutations": rejected}


if __name__ == "__main__":
    require(len(sys.argv) == 2, "usage: replay_current_route_work_failure.py ROOT")
    print(json.dumps(replay(Path(sys.argv[1])), sort_keys=True))
