#!/usr/bin/env python3
"""Independent stdlib replay of the fixed .11 attempt; executes no benchmark."""
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    require(path.stat().st_size <= 1048576, "bounded record")
    return json.loads(path.read_text())


def equal(actual, expected):
    require(json.dumps(actual, sort_keys=True) == json.dumps(expected, sort_keys=True), "record mismatch")


def checked_pair(root, label, cell, enabled, measured, expected):
    record = read(root / (label + ".stdout"))
    trace = read(root / (label + ".stderr"))
    equal({key: value for key, value in record.items() if key != "samples_ns"}, expected[cell])
    samples = record.get("samples_ns")
    require(type(samples) is list and len(samples) == (21 if measured else 0) and
            all(type(x) is int and x > 0 for x in samples), "raw samples")
    encodes = 26 if measured else 1
    matches = encodes * 2 if cell == 0 else 0
    equal(trace, {"schema": "gfni-source-stage-timing/v1", "cell": cell, "enabled": enabled,
        "encodes": encodes, "calls": encodes * (2, 2, 1, 1, 2, 1)[cell],
        "matches": matches, "changed": matches if enabled else 0, "timed": measured,
        "exercise": False})
    return record, trace


def resource_peak(path, complete):
    lines = path.read_text().splitlines()
    require(lines.count("memory.peak") == 1, "resource block")
    index = lines.index("memory.peak")
    peak = int(lines[index + 1])
    require(0 < peak <= 268435456 and lines[index + 2:] == ["memory.max", "268435456",
        "memory.events", "low 0", "high 0", "max 0", "oom 0", "oom_kill 0", "oom_group_kill 0",
        "memory.swap.current", "0", "memory.swap.max", "0"], "resource envelope")
    require("\tExit status: %d" % (0 if complete else 1) in lines, "scope exit")
    return peak


def replay(root):
    frozen, attempt = root / "frozen", root / "attempt"
    require(sha(frozen / "gfni_source_stage_screen_plan.json") ==
            "fee19236152f027e2319f022a243e8d28e3fb7298526a260f4e550e97ecd0873", "preregistered plan")
    require(sha(frozen / "pins.json") ==
            "c60f03b1f030aba96f03472acd532d46af3921f6e67de9fd10f28038330e24e8", "frozen inventory")
    pins = read(frozen / "pins.json")
    require(len(pins["files"]) == 11, "pin count")
    for name, expected_sha in pins["files"].items():
        path = frozen / name
        require(Path(name).name == name and path.is_file() and not path.is_symlink() and
                not path.stat().st_mode & 0o222 and sha(path) == expected_sha, "frozen artifact")
    plan = read(frozen / "gfni_source_stage_screen_plan.json")
    expected = read(frozen / "expected.json")
    state = read(attempt / "attempt.json")
    equal(state["pins"], pins)
    equal(state["host"], plan["host"])
    require(state["schema"] == "gfni-source-stage-screen-attempt/v1" and
            state["plan_sha256"] == sha(frozen / "gfni_source_stage_screen_plan.json"), "attempt identity")
    require(len(state["preflight"]) == 12, "preflight incomplete")
    for cell in range(6):
        for enabled in (False, True):
            record, trace = checked_pair(attempt, "check-%d-%d" % (cell, enabled),
                                        cell, enabled, False, expected)
            equal(state["preflight"][2 * cell + enabled], {"record": record, "trace": trace})
    passive = state["passive"]
    require(passive["elapsed_ns"] >= 10000000000, "short passive observation")
    rows = state["invocations"]
    require(len(rows) <= 144 and type(state["complete"]) is bool, "attempt bounds")
    orders = (("current_vs_staged", ("off", "on", "on", "off")),
              ("same_current", ("off_a", "off_b", "off_b", "off_a")))
    cursor = 0
    medians = {}
    for cell in range(6):
        for round_id in range(3):
            for comparison, order in orders:
                for slot, variant in enumerate(order):
                    if cursor >= len(rows):
                        continue
                    row = rows[cursor]
                    equal([row[k] for k in ("cell", "round", "comparison", "slot", "variant")],
                          [cell, round_id, comparison, slot, variant])
                    require(type(row["sibling_delta"]) is int and row["sibling_delta"] >= 0 and
                            (row["sibling_delta"] == 0 or cursor == len(rows) - 1), "sibling history")
                    label = "cell-%d-round-%d-%s-slot-%d" % (cell, round_id, comparison, slot)
                    record, trace = checked_pair(attempt, label, cell, variant == "on", True, expected)
                    equal(row["record"], record)
                    equal(row["trace"], trace)
                    medians[cell, round_id, comparison, slot] = statistics.median(record["samples_ns"])
                    cursor += 1
    peak = resource_peak(root / "attempt.scope.log", state["complete"])
    if not state["complete"]:
        require("analysis" not in state, "partial attempt analyzed")
        passive_failure = passive["after"] > passive["before"] and not rows
        timed_failure = (passive["after"] == passive["before"] and rows and rows[-1]["sibling_delta"] > 0)
        require(passive_failure or timed_failure, "unsupported failure requires separate investigation")
        require(state["failure"] == "ValueError: %ssibling activity; attempt stopped" %
                ("passive " if passive_failure else ""), "terminal failure")
        return {"complete": False, "timed_invocations": len(rows), "performance_inference": False,
                "failure": state["failure"], "peak_bytes": peak}
    require(len(rows) == 144 and passive["after"] == passive["before"] and
            all(row["sibling_delta"] == 0 for row in rows) and "failure" not in state,
            "complete isolation")
    derived = []
    for cell in range(6):
        result = {"cell": cell, "round_ratios": {}, "ratios": {}}
        for comparison, _ in orders:
            logs = []
            for round_id in range(3):
                log_values = [math.log(medians[cell, round_id, comparison, slot]) for slot in range(4)]
                logs.append((log_values[0] + log_values[3] - log_values[1] - log_values[2]) / 2)
            result["round_ratios"][comparison] = [math.exp(value) for value in logs]
            result["ratios"][comparison] = math.exp(statistics.mean(logs))
        derived.append(result)
    controls = all(1 / 1.02 <= cell["ratios"]["same_current"] <= 1.02 for cell in derived)
    neighbors = all(1 / 1.02 <= cell["ratios"]["current_vs_staged"] <= 1.02 for cell in derived[1:])
    target = (derived[0]["ratios"]["current_vs_staged"] >= 1.05 and
              min(derived[0]["round_ratios"]["current_vs_staged"]) > 1)
    decision = ("inconclusive_controls" if not controls or not neighbors else
                "continue_to_future_qualification" if target else "reject_for_this_screen")
    analysis = state["analysis"]
    require(set(analysis) == {"decision", "cells", "confidence_intervals", "production_promotion",
                              "exact_leopard1_claim", "authoritative_v19"}, "analysis schema")
    require(analysis["decision"] == decision and all(analysis[k] is False for k in
            ("confidence_intervals", "production_promotion", "exact_leopard1_claim", "authoritative_v19")),
            "decision or unsupported claim")
    require(len(analysis["cells"]) == 6, "analysis cells")
    for actual, calculated in zip(analysis["cells"], derived):
        require(set(actual) == set(calculated) and actual["cell"] == calculated["cell"], "cell schema")
        for field in ("ratios", "round_ratios"):
            require(set(actual[field]) == set(calculated[field]), "comparison schema")
        for comparison, _ in orders:
            require(math.isclose(actual["ratios"][comparison], calculated["ratios"][comparison],
                                 rel_tol=1e-12, abs_tol=1e-12), "aggregate ratio")
            require(len(actual["round_ratios"][comparison]) == 3 and all(
                math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12) for a, b in
                zip(actual["round_ratios"][comparison], calculated["round_ratios"][comparison])), "round ratios")
    return {"complete": True, "timed_invocations": 144, "preflights": 12,
            "decision": decision, "cells": derived, "peak_bytes": peak,
            "production_promotion": False, "exact_leopard1_claim": False,
            "performance_inference": decision != "inconclusive_controls"}


if __name__ == "__main__":
    require(len(sys.argv) == 2, "usage: replay_gfni_source_stage_screen.py ROOT")
    print(json.dumps(replay(Path(sys.argv[1])), sort_keys=True))
