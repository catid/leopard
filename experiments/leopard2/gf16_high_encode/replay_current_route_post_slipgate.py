#!/usr/bin/env python3
"""Standalone offline replay of the complete post-Slipgate diagnostic.

Imports no collector or codec; executes no benchmark. Recomputes every ratio
from raw, ordered samples and checks the fixed evidence and host condition.
"""
import copy
import hashlib
import json
import math
from pathlib import Path
import sys

PLAN = "current_route_screen_post_slipgate_plan.json"
ORDERS = (("main_vs_current", ("main", "current", "current", "main")),
          ("same_current", ("current_a", "current_b", "current_b", "current_a")))


def require(condition, message):
    if not condition:
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


def derive(rows, expected):
    require(type(rows) is list and len(rows) == 144, "incomplete attempt")
    medians = []
    for index, row in enumerate(rows):
        cell, round_id = index // 24, index % 24 // 8
        comparison, order = ORDERS[index % 8 // 4]
        slot = index % 4
        variant = order[slot]
        require(all(type(row[key]) is int for key in
                    ("cell", "round", "slot", "sibling_delta")), "integer metadata")
        require((row["cell"], row["round"], row["slot"], row["comparison"],
                 row["variant"], row["sibling_delta"]) ==
                (cell, round_id, slot, comparison, variant, 0), "order/isolation")
        record = row["record"]
        samples = record["samples_ns"]
        require(type(samples) is list and len(samples) == 21 and
                all(type(x) is int and x > 0 for x in samples), "sample values")
        identity = {k: v for k, v in record.items() if k != "samples_ns"}
        require(identity == expected["main" if variant == "main" else "current"][cell],
                "workload/source/route/scratch changed")
        medians.append(sorted(samples)[10])
    cells = []
    for cell in range(6):
        round_ratios = {}
        for comparison_index, (name, _) in enumerate(ORDERS):
            values = []
            for round_id in range(3):
                start = cell * 24 + round_id * 8 + comparison_index * 4
                a, b, c, d = medians[start:start + 4]
                values.append(math.sqrt((a * d) / (b * c)))
            round_ratios[name] = values
        ratios = {name: math.prod(values) ** (1 / 3)
                  for name, values in round_ratios.items()}
        cells.append({"cell": cell, "round_ratios": round_ratios, "ratios": ratios})
    controls = all(1 / 1.02 <= cell["ratios"]["same_current"] <= 1.02 for cell in cells)
    for cell in cells:
        ratio = cell["ratios"]["main_vs_current"]
        rounds = cell["round_ratios"]["main_vs_current"]
        if not controls:
            classification = "no_inference_controls"
        elif ratio < 1 / 1.02 and all(x < 1 for x in rounds):
            classification = "investigate_current_deficit"
        elif ratio > 1.02 and all(x > 1 for x in rounds):
            classification = "current_advantage"
        else:
            classification = "near_parity_or_uncertain"
        cell["interpretation"] = classification
    return {"decision": "diagnostic_complete" if controls else "inconclusive_controls",
            "cells": cells, "confidence_intervals": False, "authoritative_v19": False,
            "production_promotion": False, "historical_exact_main_gap_closed": False}


def same_analysis(actual, derived):
    require(set(actual) == set(derived), "analysis fields")
    for key in derived.keys() - {"cells"}:
        require(type(actual[key]) is type(derived[key]) and actual[key] == derived[key],
                "analysis decision/claim")
    require(len(actual["cells"]) == 6, "cell count")
    for observed, rebuilt in zip(actual["cells"], derived["cells"]):
        require(set(observed) == set(rebuilt) and observed["cell"] == rebuilt["cell"] and
                observed["interpretation"] == rebuilt["interpretation"], "cell classification")
        for key in ("round_ratios", "ratios"):
            require(set(observed[key]) == set(rebuilt[key]), "comparison keys")
            for name in rebuilt[key]:
                a = observed[key][name] if key == "round_ratios" else [observed[key][name]]
                b = rebuilt[key][name] if key == "round_ratios" else [rebuilt[key][name]]
                require(len(a) == len(b) and all(math.isclose(x, y, rel_tol=2e-14,
                        abs_tol=0) for x, y in zip(a, b)), "recomputed ratios")


def replay(root):
    for name, digest in {
        "attempt1/attempt.json": "a1ebcc14caf096152317175f4cebd5eaaf83d29804b9b0aadd94e314b441ffc4",
        "attempt1.log": "491f30599bf60e1fdbaba1c3060aee778a734b69b5265588d3ad7d53d538310c",
        "frozen/pins.json": "28eaea50b852f35c9ba8ba92cd693be70559d911070b4dcf67416cd93a2df17d",
    }.items():
        require(sha(root / name) == digest, "fixed evidence: " + name)
    pins = read_json(root / "frozen/pins.json")
    require(pins["source_commit"] == "23bcdd7f45bc6dc00eeaf51b30090d4e901d6f3c" and
            len(pins["files"]) == 16, "preregistration/inventory")
    for name, digest in pins["files"].items():
        path = root / "frozen" / name
        require(Path(name).name == name and not path.is_symlink() and
                not path.stat().st_mode & 0o222 and sha(path) == digest, "frozen: " + name)
    plan = read_json(root / "frozen" / PLAN)
    require(plan["condition"] == "slipgate-disabled-20260909" and
            (plan["cpu"], plan["sibling"], plan["controller_cpu"], plan["attempt_budget"])
            == (26, 90, 0, 1), "local post-shutdown plan")
    for name, digest in plan["artifact_sha256"].items():
        require(pins["files"][name] == digest, "plan artifact")
    state = read_json(root / "attempt1/attempt.json")
    expected = read_json(root / "frozen/expected.json")
    require(state["pins"] == pins and state["plan_sha256"] == pins["files"][PLAN] and
            state["host"] == plan["host"], "journal bindings")
    require(state["complete"] is True and "failure" not in state, "completion")
    require(state["passive"] == {"before": 570054, "after": 570054,
                                "elapsed_ns": 10000241142}, "passive evidence")
    require(len(state["preflight"]) == 12, "preflight count")
    for cell in range(6):
        for index, variant in enumerate(("main", "current")):
            record = dict(expected[variant][cell], samples_ns=[])
            prefix = root / "attempt1" / f"check-{cell}-{variant}"
            require(state["preflight"][cell * 2 + index] == record and
                    read_json(prefix.with_suffix(".stdout")) == record, "raw check")
            require(prefix.with_suffix(".stderr").stat().st_size == 0, "check stderr")
    for row in state["invocations"]:
        prefix = root / "attempt1" / (
            f"cell-{row['cell']}-round-{row['round']}-{row['comparison']}-slot-{row['slot']}")
        require(read_json(prefix.with_suffix(".stdout")) == row["record"], "raw timing")
        require(prefix.with_suffix(".stderr").stat().st_size == 0, "timing stderr")
    condition = []
    for unit in ("slipgate-catstream.service", "slipgate-obs.service", "slipgate-headless-x.service"):
        condition.extend(("MainPID=0", "Id=" + unit, "ActiveState=inactive", "UnitFileState=disabled"))
    for container in plan["shutdown_condition"]["containers"]:
        condition.append(container + " exited|false|no")
    for name in ("condition-before.log", "condition-after.log"):
        require((root / name).read_text().splitlines() == condition, "shutdown condition")
    log = (root / "attempt1.log").read_text()
    require("\tExit status: 0\n" in log and
            "memory.peak\n133763072\nmemory.max\n268435456\n" in log and
            "memory.events\nlow 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\n"
            "memory.swap.current\n0\nmemory.swap.max\n0\n" in log, "resource evidence")
    derived = derive(state["invocations"], expected)
    same_analysis(state["analysis"], derived)
    rejected = 0
    for key, value in (("cell", 1), ("round", 1), ("slot", 1), ("sibling_delta", 1),
                       ("sibling_delta", False), ("variant", "current"),
                       ("comparison", "same_current"), ("partial", None),
                       ("samples", [0] * 21), ("samples", [True] * 21),
                       ("samples", [1] * 20), ("identity", "changed")):
        changed = copy.deepcopy(state["invocations"])
        if key == "partial":
            changed.pop()
        elif key == "samples":
            changed[0]["record"]["samples_ns"] = value
        elif key == "identity":
            changed[0]["record"]["codec_commit"] = value
        else:
            changed[0][key] = value
        try:
            derive(changed, expected)
        except ValueError:
            rejected += 1
    require(rejected == 12, "semantic mutation escaped")
    for key, value in (("decision", "wrong"), ("production_promotion", True)):
        changed = copy.deepcopy(state["analysis"])
        changed[key] = value
        try:
            same_analysis(changed, derived)
        except ValueError:
            rejected += 1
    changed = copy.deepcopy(state["analysis"])
    changed["cells"][0]["ratios"]["main_vs_current"] *= 1.001
    try:
        same_analysis(changed, derived)
    except ValueError:
        rejected += 1
    require(rejected == 15, "analysis mutation escaped")
    contaminated_control = copy.deepcopy(state["invocations"])
    for row in contaminated_control:
        if row["comparison"] == "same_current" and row["variant"] == "current_a":
            row["record"]["samples_ns"] = [x * 2 for x in row["record"]["samples_ns"]]
    failed_control = derive(contaminated_control, expected)
    require(failed_control["decision"] == "inconclusive_controls" and
            all(cell["interpretation"] == "no_inference_controls"
                for cell in failed_control["cells"]), "control gate self-check")
    return {"checks": 12, "timed_invocations": 144, "sibling_nonidle_jiffies": 0,
            "rejected_mutations": rejected, "control_failure_selfcheck": True,
            "analysis": derived}


if __name__ == "__main__":
    require(len(sys.argv) == 2, "usage: replay_current_route_post_slipgate.py ROOT")
    print(json.dumps(replay(Path(sys.argv[1])), sort_keys=True))
