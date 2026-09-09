#!/usr/bin/env python3
"""Independent retained-only replay; imports no collector and runs no codec."""
import copy
import hashlib
import json
import math
from pathlib import Path
import sys

ORDERS = (("auto_vs_gfni", ("auto", "gfni", "gfni", "auto")),
          ("main_vs_gfni", ("main", "gfni", "gfni", "main")),
          ("same_gfni", ("gfni_a", "gfni_b", "gfni_b", "gfni_a")))


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            result.update(block)
    return result.hexdigest()


def read(path):
    require(path.stat().st_size < 1048576, "bounded JSON")
    return json.loads(path.read_text())


def same(actual, expected):
    require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), "exact typed identity")


def derive(rows, expected):
    require(type(rows) is list and len(rows) == 72, "complete ordered attempt")
    medians = []
    for index, row in enumerate(rows):
        cell, round_id = index // 36, index % 36 // 12
        comparison, order = ORDERS[index % 12 // 4]
        slot = index % 4
        variant = order[slot]
        same({key: row[key] for key in
              ("cell", "round", "comparison", "slot", "variant", "sibling_delta")},
             {"cell": cell, "round": round_id, "comparison": comparison,
              "slot": slot, "variant": variant, "sibling_delta": 0})
        record = row["record"]
        route = "gfni" if variant.startswith("gfni") else variant
        same({key: value for key, value in record.items() if key != "samples_ns"},
             expected[route][cell])
        samples = record["samples_ns"]
        require(type(samples) is list and len(samples) == 21 and
                all(type(x) is int and x > 0 for x in samples), "positive integer samples")
        medians.append(sorted(samples)[10])
    cells = []
    for cell in range(2):
        rounds = {}
        for comparison_index, (name, _) in enumerate(ORDERS):
            values = []
            for round_id in range(3):
                start = cell * 36 + round_id * 12 + comparison_index * 4
                a, b, c, d = medians[start:start + 4]
                values.append(math.sqrt((a * d) / (b * c)))
            rounds[name] = values
        cells.append({"cell": cell, "round_ratios": rounds,
                      "ratios": {name: math.prod(values) ** (1 / 3)
                                 for name, values in rounds.items()}})
    controls_ok = all(1 / 1.02 <= cell["ratios"]["same_gfni"] <= 1.02 for cell in cells)
    for cell in cells:
        cell["decision"] = ("inconclusive_controls" if not controls_ok else
            "qualify_bounded_auto_candidate" if cell["ratios"]["auto_vs_gfni"] >= 1.05 and
            all(x > 1 for x in cell["round_ratios"]["auto_vs_gfni"]) else "reject_for_this_screen")
    return {"controls_pass": controls_ok, "cells": cells, "confidence_intervals": False,
            "production_promotion": False, "neighbor_qualification": False,
            "authoritative_v19": False}


def analysis_equal(actual, derived):
    require(set(actual) == set(derived), "analysis fields")
    for key in derived.keys() - {"cells"}:
        same(actual[key], derived[key])
    require(len(actual["cells"]) == 2, "two cells")
    for a, b in zip(actual["cells"], derived["cells"]):
        require(set(a) == set(b), "cell fields")
        same(a["cell"], b["cell"])
        same(a["decision"], b["decision"])
        for key in ("round_ratios", "ratios"):
            require(set(a[key]) == set(b[key]), "comparison fields")
            for name in b[key]:
                left = a[key][name] if key == "round_ratios" else [a[key][name]]
                right = b[key][name] if key == "round_ratios" else [b[key][name]]
                require(len(left) == len(right) and all(math.isclose(x, y,
                        rel_tol=2e-14, abs_tol=0) for x, y in zip(left, right)), "raw ratio replay")


def resource(path, peak, maximum):
    lines = path.read_text().splitlines()
    require(lines.count("memory.peak") == 1 and "\tExit status: 0" in lines, "scope exit")
    index = lines.index("memory.peak")
    same(lines[index:], ["memory.peak", str(peak), "memory.max", str(maximum),
        "memory.events", "low 0", "high 0", "max 0", "oom 0", "oom_kill 0",
        "oom_group_kill 0", "memory.swap.current", "0", "memory.swap.max", "0"])


def replay(root):
    for name, digest in {
        "attempt1/attempt.json": "337d847c992973d5b0d54a05914a4cc09f890b7faf073dc3b8f9e401da72b6be",
        "attempt1.log": "ca4621f0fe80f989e7bfe6f89e59864055679b84698afc816fa571f5a5b2cbba",
        "frozen/pins.json": "9e6e63281427faacb51fb4c581ba239cfd24f81936b11de8e600105c87233862",
    }.items():
        require(sha(root / name) == digest, "fixed evidence: " + name)
    pins = read(root / "frozen/pins.json")
    require(pins["source_commit"] == "d1d3cbfbe6e2fa94460c32bc374f1a0950aaf2e3" and
            len(pins["files"]) == 12, "preregistration and inventory")
    for name, digest in pins["files"].items():
        path = root / "frozen" / name
        require(Path(name).name == name and not path.is_symlink() and
                not path.stat().st_mode & 0o222 and sha(path) == digest, "frozen input: " + name)
    plan = read(root / "frozen/gfni_boundary_screen_plan.json")
    state = read(root / "attempt1/attempt.json")
    expected = read(root / "frozen/expected.json")
    same(state["pins"], pins)
    same(state["host"], plan["host"])
    same(state["plan_sha256"], pins["files"]["gfni_boundary_screen_plan.json"])
    require(state["complete"] is True and "failure" not in state, "completion")
    same(state["passive"], {"before": 570109, "after": 570109, "elapsed_ns": 10000091881})
    for name, digest in plan["artifact_sha256"].items():
        same(pins["files"][name], digest)
    build_artifacts = 0
    for line in (root / "build.log").read_text().splitlines():
        words = line.split()
        if len(words) == 2 and len(words[0]) == 64 and \
                words[1].startswith("/tmp/leopard-gfni-boundary.iIwWWx/stage/"):
            require(all(c in "0123456789abcdef" for c in words[0]), "build hash format")
            require(sha(root / "stage" / Path(words[1]).name) == words[0], "build artifact changed")
            build_artifacts += 1
    require(build_artifacts == 10, "complete build artifact inventory")
    require(len(state["preflight"]) == 8, "eight preflights")
    parity_bytes = 0
    for cell in range(2):
        for index, route in enumerate(("main", "auto", "avx2", "gfni")):
            record = dict(expected[route][cell], samples_ns=[])
            same(state["preflight"][cell * 4 + index], record)
            same(read(root / "attempt1" / f"check-{cell}-{route}.stdout"), record)
            require((root / "attempt1" / f"check-{cell}-{route}.stderr").stat().st_size == 0,
                    "preflight stderr")
            same(read(root / "checks" / f"{cell}-{route}.json"), record)
            if route != "main":
                same(read(root / "checks" / f"{cell}-{route}-sanitize.json"), record)
                actual = root / "checks" / f"{cell}-{route}.parity"
                original = root / "checks" / f"{cell}-main.parity"
                size = plan["cells"][cell]["r"] * plan["cells"][cell]["bytes"]
                require(actual.stat().st_size == original.stat().st_size == size, "parity size")
                with actual.open("rb") as a, original.open("rb") as b:
                    while True:
                        block = a.read(65536)
                        require(block == b.read(65536), "full Leopard1 parity")
                        if not block:
                            break
                parity_bytes += size
    for cell in range(8):
        r = (199 if cell % 2 else 200) if cell < 6 else 7
        size = ((65536 if cell % 2 else 32768) + (2 if cell >= 4 else 0)) if cell < 6 else \
            (65 if cell == 6 else 66)
        record = {"schema": "leopard-gfni-boundary-guards/v1", "cell": cell,
                  "k": 1000 if cell < 6 else 17, "r": r, "bytes": size,
                  "misalignment": 0 if cell < 2 else 1 if cell < 4 or cell >= 6 else 2,
                  "field": 1 if cell == 6 else 2, "subset_masks": 6,
                  "scratch_bytes": read(root / "checks" / f"guard-{cell}-release.json")["scratch_bytes"],
                  "timed": False}
        require(type(record["scratch_bytes"]) is int and record["scratch_bytes"] > 0, "guard scratch")
        for mode in ("release", "sanitize"):
            same(read(root / "checks" / f"guard-{cell}-{mode}.json"), record)
    for row in state["invocations"]:
        prefix = root / "attempt1" / (
            f"cell-{row['cell']}-round-{row['round']}-{row['comparison']}-slot-{row['slot']}")
        same(read(prefix.with_suffix(".stdout")), row["record"])
        require(prefix.with_suffix(".stderr").stat().st_size == 0, "timed stderr")
    condition = []
    for unit in ("slipgate-catstream.service", "slipgate-obs.service", "slipgate-headless-x.service"):
        condition.extend(("MainPID=0", "Id=" + unit, "ActiveState=inactive", "UnitFileState=disabled"))
    for container in ("3fad3dfb247b3d4f3991e0a6801fba0805aae0901bdc89018841d4707a2d888c",
                      "002ba0a05bf860b12f2d799b49da1976aed663cb1494b2f8236751e18dc50182"):
        condition.append(container + " exited|false|no")
    for name in ("condition-before", "condition-after"):
        same((root / "attempt1" / (name + ".stdout")).read_text().splitlines(), condition)
        require((root / "attempt1" / (name + ".stderr")).stat().st_size == 0, "condition stderr")
    resource(root / "build.log", 130965504, 536870912)
    resource(root / "check.log", 256364544, 268435456)
    resource(root / "attempt1.log", 127823872, 268435456)
    derived = derive(state["invocations"], expected)
    analysis_equal(state["analysis"], derived)
    rejected = 0
    for key, value in (("cell", 1), ("round", 1), ("slot", 1), ("sibling_delta", 1),
                       ("sibling_delta", False), ("variant", "gfni"),
                       ("comparison", "same_gfni"), ("partial", None),
                       ("samples", [True] * 21), ("samples", [0] * 21),
                       ("samples", [1] * 20), ("identity", "changed")):
        changed = copy.deepcopy(state["invocations"])
        if key == "partial":
            changed.pop()
        elif key == "samples":
            changed[0]["record"]["samples_ns"] = value
        elif key == "identity":
            changed[0]["record"]["requested"] = value
        else:
            changed[0][key] = value
        try:
            derive(changed, expected)
        except ValueError:
            rejected += 1
    require(rejected == 12, "row mutation escaped")
    for key in ("production_promotion", "neighbor_qualification", "authoritative_v19"):
        changed = copy.deepcopy(derived)
        changed[key] = True
        try:
            analysis_equal(changed, derived)
        except ValueError:
            rejected += 1
    require(rejected == 15, "claim mutation escaped")
    bad_control = copy.deepcopy(state["invocations"])
    for row in bad_control:
        if row["variant"] == "gfni_a":
            row["record"]["samples_ns"] = [2 * x for x in row["record"]["samples_ns"]]
    rejected_control = derive(bad_control, expected)
    require(rejected_control["controls_pass"] is False and all(
        cell["decision"] == "inconclusive_controls" for cell in rejected_control["cells"]),
        "bad control must suppress both decisions")
    return {"preflights": 8, "timed_invocations": 72, "sibling_nonidle_jiffies": 0,
            "compared_parity_bytes": parity_bytes, "guarded_shapes_per_build": 8,
            "rejected_mutations": rejected, "control_failure_selfcheck": True,
            "analysis": derived}


if __name__ == "__main__":
    require(len(sys.argv) == 2, "usage: replay_gfni_boundary_screen.py ROOT")
    print(json.dumps(replay(Path(sys.argv[1])), sort_keys=True))
