#!/usr/bin/env python3
"""Offline deterministic replay of .38.5.4.10.1; no codec execution/timings."""
import json
from pathlib import Path
import sys


def require(value, message):
    if not value:
        raise ValueError(message)


def expected(cell, variant):
    k, r, size = ((1000, 200, 65536), (1000, 200, 65536),
                  (1000, 200, 65536), (1000, 200, 32768),
                  (1000, 199, 65536), (4096, 512, 4096))[cell]
    passes = 2 if variant != "main" and cell in (0, 1, 4) else 1
    copied = variant == "main" or cell != 5
    return {"cell": cell, "encode_calls": 1,
            "input_copy_calls": k * passes if copied else 0,
            "input_copy_bytes": k * size if copied else 0,
            "output_copy_calls": r * passes if copied else 0,
            "output_copy_bytes": r * size if copied else 0,
            "other_copy_calls": 0, "other_copy_bytes": 0,
            "zero_calls": 24 * passes if cell != 5 else 0,
            "zero_bytes": 24 * size if cell != 5 else 0,
            "other_set_calls": 0, "other_set_bytes": 0,
            "checked_copy_calls": 0, "checked_set_calls": 0,
            "external_calls_only": True, "timed": False}


def same_record(actual, want):
    # bool/int equality must not silently admit an altered schema.
    require(set(actual) == set(want) and all(
        type(actual[key]) is type(value) and actual[key] == value
        for key, value in want.items()), "copy count/schema mismatch")


def resource_peak(path):
    lines = path.read_text().splitlines()
    require(lines.count("memory.peak") == 1, "resource block count")
    index = lines.index("memory.peak")
    peak = int(lines[index + 1])
    require(0 < peak <= 268435456 and lines[index + 2:] == [
        "memory.max", "268435456", "memory.events", "low 0", "high 0",
        "max 0", "oom 0", "oom_kill 0", "oom_group_kill 0",
        "memory.swap.current", "0", "memory.swap.max", "0"], "resource limits")
    require("\tExit status: 0" in lines, "scope exit")
    return peak


def replay(root, baseline):
    negative_checks = 0
    sample = expected(0, "current")
    for key, value in sample.items():
        bad = dict(sample)
        bad[key] = not value if type(value) is bool else value + 1
        try:
            same_record(bad, sample)
        except ValueError:
            negative_checks += 1
        else:
            raise ValueError("mutation accepted")
    peaks = []
    parity_bytes = 0
    for variant in ("main", "current", "sanitize"):
        for cell in range(6):
            prefix = "%d-%s" % (cell, variant)
            same_record(json.loads((root / (prefix + ".copy.json")).read_text()),
                        expected(cell, variant))
            ref = "%d-%s.json" % (cell, "main" if variant == "main" else "current")
            require((root / (prefix + ".json")).read_bytes() ==
                    (baseline / ref).read_bytes(), "workload byte identity")
            peaks.append(resource_peak(root / (prefix + ".scope.log")))
            if variant == "sanitize":
                continue
            parity = root / (prefix + ".parity")
            reference = baseline / ("%d-main.parity" % cell)
            record = json.loads((baseline / ref).read_text())
            size = record["r"] * record["bytes"]
            require(parity.stat().st_size == reference.stat().st_size == size,
                    "full parity length")
            with parity.open("rb") as actual, reference.open("rb") as original:
                while True:
                    block = actual.read(65536)
                    require(block == original.read(65536), "full parity bytes")
                    if not block:
                        break
            parity_bytes += size
    return {"copy_records": 18, "workload_records": 18, "parity_files": 12,
            "compared_parity_bytes": parity_bytes, "peak_bytes": max(peaks),
            "rejected_count_mutations": negative_checks, "performance_inference": False}


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: verify_current_route_copies.py CHECKS BASELINE")
    print(json.dumps(replay(Path(sys.argv[1]), Path(sys.argv[2])), sort_keys=True))
