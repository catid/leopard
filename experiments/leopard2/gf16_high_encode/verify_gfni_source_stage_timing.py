#!/usr/bin/env python3
"""Retained-only timing-front-end validation; never executes codecs or clocks."""
import json
from pathlib import Path
import sys

from verify_gfni_source_stage import equal, require, resource_peak


def replay(root, baseline, original_stage):
    peaks = []
    compared = 0
    for flavor in ("release", "sanitize"):
        for capacity in (16, 64):
            path = root / "checks" / ("%s-capacity-%d.log" % (flavor, capacity))
            peaks.append(resource_peak(path))
            line = "source-stage capacity %d: both modes, overflow delegation, reset, neighbor passed" % capacity
            require(path.read_text().splitlines().count(line) == 1, "capacity test")
        for mode in (0, 1):
            for cell in range(6):
                for operation in (("check", "exercise") if cell == 0 else ("check",)):
                    prefix = root / "checks" / ("%s-%d-%d-%s" % (flavor, mode, cell, operation))
                    peaks.append(resource_peak(prefix.with_suffix(".log")))
                    expected = (baseline / ("%d-current.json" % cell)).read_bytes()
                    require(prefix.with_suffix(".jsonl").read_bytes() ==
                            expected * (26 if operation == "exercise" else 1), "workload records")
                    encodes = 26 if operation == "exercise" else 1
                    matches = encodes * 2 if cell == 0 else 0
                    equal(json.loads(prefix.with_suffix(".trace.json").read_text()), {
                        "schema": "gfni-source-stage-timing/v1", "cell": cell,
                        "enabled": bool(mode), "encodes": encodes,
                        "calls": encodes * (2, 2, 1, 1, 2, 1)[cell], "matches": matches,
                        "changed": matches if mode else 0, "timed": False,
                        "exercise": operation == "exercise"})
                    if flavor == "release" and operation == "check":
                        parity = prefix.with_suffix(".parity")
                        reference = baseline / ("%d-main.parity" % cell)
                        require(parity.stat().st_size == reference.stat().st_size, "parity size")
                        with parity.open("rb") as actual, reference.open("rb") as original:
                            while True:
                                block = actual.read(65536)
                                require(block == original.read(65536), "exact Leopard1 parity")
                                if not block:
                                    break
                        compared += parity.stat().st_size
    for cell in range(6):
        peaks.append(resource_peak(root / "checks" / ("default-%d.log" % cell)))
        for suffix in ("json", "probe.jsonl", "parity"):
            actual_path = root / "checks" / ("default-%d.%s" % (cell, suffix))
            original_path = original_stage / "checks" / ("release-1-%d.%s" % (cell, suffix))
            require(actual_path.stat().st_size == original_path.stat().st_size, "default probe size")
            with actual_path.open("rb") as actual, original_path.open("rb") as original:
                while True:
                    block = actual.read(65536)
                    require(block == original.read(65536), "default probe bytes")
                    if not block:
                        break
    peaks.append(resource_peak(root / "guards.log"))
    for flavor in ("release", "sanitize"):
        for index in range(8):
            prefix = root / "guards" / ("%s-%d" % (flavor, index))
            require(prefix.with_suffix(".stdout").read_bytes() == b"", "guard stdout")
            lines = prefix.with_suffix(".stderr").read_text().splitlines()
            require(len(lines) == 1 and lines[0].startswith("source-stage timing driver: "),
                    "guard diagnostic")
    require(compared == 122028032, "parity byte total")
    return {"fixed_records": 24, "clock_free_exercise_records": 104,
            "capacity_tests": 4, "unchanged_default_probe_cases": 6,
            "malformed_request_guards": 16, "compared_parity_bytes": compared,
            "default_probe_parity_bytes": 61014016,
            "validated_scopes": len(peaks), "peak_bytes": max(peaks),
            "timed_invocations": 0, "performance_inference": False}


if __name__ == "__main__":
    require(len(sys.argv) == 4, "usage: verify_gfni_source_stage_timing.py ROOT BASELINE ORIGINAL_STAGE")
    print(json.dumps(replay(*(Path(value) for value in sys.argv[1:])), sort_keys=True))
