"""Fixed release encode timing payload and complete-ABBA analysis contracts."""
import math
import statistics

from check_encode import CELLS, BASELINE, validate, require

CANDIDATE = "e35b1f0c3a5ef9f241ad3b174b0fbfdb61bbccab"
GROUPS = (4194304, 1048576, 256, 256, 256, 16, 16, 128)
SAMPLES = 9
ROUNDS = 3
MIN_GROUP_NS = 20000000
COMPARISONS = (
    ("native_vs_current", ("main", "current", "current", "main")),
    ("same_native", ("main", "main", "main", "main")),
    ("same_current", ("current", "current", "current", "current")),
)


def validate_timing(record, index, implementation, expected, mode="measure"):
    require(mode in {"measure", "synthetic", "exercise"}, "unsupported timing mode")
    require(type(record) is dict and record.get("schema") ==
            "leopard-native-release-encode-timing/v1", "timing schema")
    base = dict(record)
    timing = base.pop("grouped_timing", None)
    total = 6 + SAMPLES * GROUPS[index]
    require(type(base.get("public_encode_calls")) is int and
            base["public_encode_calls"] == total, "total public call count")
    base["schema"] = "leopard-native-release-encode-check/v1"
    base["public_encode_calls"] = 2
    validate(base, index, implementation, BASELINE if implementation == "main" else CANDIDATE)
    require(base == expected, "workload, source, output or scratch changed")
    require(type(timing) is dict and set(timing) == {
        "schema", "mode", "clock_kind", "group_calls", "warmup_calls",
        "total_public_calls", "clock_calls", "elapsed_ns", "public_calls_at_clock"},
        "timing fields")
    exact = {"schema": "native-encode-groups/v1", "mode": mode,
             "clock_kind": "synthetic" if mode == "synthetic" else
                           "abort" if mode == "exercise" else "steady",
             "group_calls": GROUPS[index], "warmup_calls": 4,
             "total_public_calls": total, "clock_calls": 0 if mode == "exercise" else 18}
    for key, value in exact.items():
        require(type(timing[key]) is type(value) and timing[key] == value,
                "timing mismatch: " + key)
    samples = timing["elapsed_ns"]
    require(type(samples) is list and len(samples) == SAMPLES and
            all(type(value) is int and 0 <= value <= 9007199254740991 for value in samples),
            "invalid sample array")
    if mode == "measure":
        require(min(samples) >= MIN_GROUP_NS, "retained timer window below 20ms")
    else:
        require(samples == ([31000000 + i * 100 for i in range(SAMPLES)]
                            if mode == "synthetic" else [0] * SAMPLES), "synthetic durations")
    endpoints = timing["public_calls_at_clock"]
    expected_endpoints = [6 + ((endpoint + 1) // 2) * GROUPS[index]
                          for endpoint in range(18)] if mode == "synthetic" else []
    require(type(endpoints) is list and all(type(value) is int for value in endpoints) and
            endpoints == expected_endpoints, "clock/public-call boundary mismatch")
    return [value / GROUPS[index] for value in samples]


def schedule():
    return [(cell, round_id, comparison, slot, implementation)
            for cell in range(len(CELLS)) for round_id in range(ROUNDS)
            for comparison, order in COMPARISONS for slot, implementation in enumerate(order)]


def analyze(rows, expected):
    order = schedule()
    require(len(rows) == len(order), "incomplete attempts have no performance analysis")
    medians = []
    for row, identity in zip(rows, order):
        require(all(type(row[key]) is int for key in ("cell", "round", "slot", "sibling_delta")),
                "invalid invocation metadata type")
        require((row["cell"], row["round"], row["comparison"], row["slot"],
                 row["implementation"]) == identity and row["sibling_delta"] == 0,
                "invocation order or sibling contamination")
        cell, _, _, _, implementation = identity
        samples = validate_timing(row["record"], cell, implementation,
                                  expected[cell][implementation])
        medians.append(statistics.median(samples))
    cells = []
    cursor = 0
    for cell in range(len(CELLS)):
        ratios = {name: [] for name, _ in COMPARISONS}
        comparison_medians = {"main": [], "current": []}
        for _ in range(ROUNDS):
            for comparison, order in COMPARISONS:
                values = medians[cursor:cursor + 4]
                cursor += 4
                ratios[comparison].append(math.sqrt(values[0] * values[3] /
                                                    (values[1] * values[2])))
                if comparison == "native_vs_current":
                    comparison_medians["main"].extend((values[0], values[3]))
                    comparison_medians["current"].extend((values[1], values[2]))
        aggregates = {key: math.exp(statistics.mean(math.log(v) for v in values))
                      for key, values in ratios.items()}
        cells.append({"cell": cell, "id": CELLS[cell][0], "round_ratios": ratios,
                      "ratios": aggregates, "comparison_process_medians_ns": comparison_medians})
    controls_ok = all(1 / 1.02 <= cell["ratios"][key] <= 1.02
                      for cell in cells for key in ("same_native", "same_current"))
    for cell in cells:
        ratio = cell["ratios"]["native_vs_current"]
        rounds = cell["round_ratios"]["native_vs_current"]
        cell["classification"] = (
            "inconclusive_controls" if not controls_ok else
            "current_advantage" if ratio > 1.02 and min(rounds) > 1 else
            "current_deficit" if ratio < 1 / 1.02 and max(rounds) < 1 else
            "near_parity_or_uncertain")
    return {"decision": "complete" if controls_ok else "inconclusive_controls",
            "cells": cells, "confidence_intervals": False, "promotion": False,
            "boundary": "ordinary public full encode, setup and allocation excluded"}
