#!/usr/bin/env python3
"""Retained-only .38.5.4.11 replay; does not launch codecs or infer speed."""
import json
from pathlib import Path
import re
import sys


def require(value, message):
    if not value:
        raise ValueError(message)


def equal(actual, expected):
    require(json.dumps(actual, sort_keys=True) == json.dumps(expected, sort_keys=True),
            "record values/schema differ")


def resource_peak(path):
    lines = path.read_text().splitlines()
    require(lines.count("memory.peak") == 1, "resource block count")
    index = lines.index("memory.peak")
    peak = int(lines[index + 1])
    require(0 < peak <= 268435456 and lines[index + 2:] == [
        "memory.max", "268435456", "memory.events", "low 0", "high 0", "max 0",
        "oom 0", "oom_kill 0", "oom_group_kill 0", "memory.swap.current", "0",
        "memory.swap.max", "0"] and "\tExit status: 0" in lines, "resource envelope")
    return peak


def replay(root, baseline, copies):
    matrix = ((1000, 200, 65536, 6, 32768, 2), (1000, 200, 65536, 3, 32768, 2),
              (1000, 200, 65536, 5, 65536, 1), (1000, 200, 32768, 3, 32768, 1),
              (1000, 199, 65536, 3, 32768, 2), (4096, 512, 4096, 3, 4096, 1))
    peaks = []
    compared = 0
    for flavor in ("release", "sanitize"):
        for mode in (0, 1):
            for cell, (k, r, size, kind, tile, passes) in enumerate(matrix):
                prefix = root / "checks" / ("%s-%d-%d" % (flavor, mode, cell))
                require(prefix.with_suffix(".json").read_bytes() ==
                        (baseline / ("%d-current.json" % cell)).read_bytes(), "public workload")
                raw = prefix.with_suffix(".probe.jsonl").read_text().splitlines()
                require(len(raw) == 2, "two exact probe records")
                actual_copy, actual_stage = map(json.loads, raw)
                expected_copy = json.loads((copies / ("%d-current.copy.json" % cell)).read_text())
                changed = cell == 0 and mode == 1
                if changed:
                    for key in ("input_copy_calls", "input_copy_bytes", "output_copy_calls", "output_copy_bytes"):
                        expected_copy[key] = 0
                equal(actual_copy, expected_copy)
                record = {"kind": kind, "k": k, "r": r, "requested": r,
                          "side": 512 if cell == 5 else 256, "sparse_blocks": 0,
                          "bytes": tile, "source_policy": size,
                          "effective_policy": 16384 if changed else size}
                equal(actual_stage, {"schema": "gfni-source-stage/v1", "enabled": bool(mode),
                    "matches": 2 if cell == 0 else 0, "changed": 2 if changed else 0,
                    "calls": [record] * passes, "timed": False})
                peaks.append(resource_peak(prefix.with_suffix(".scope.log")))
                if flavor == "release":
                    parity = prefix.with_suffix(".parity")
                    reference = baseline / ("%d-main.parity" % cell)
                    require(parity.stat().st_size == reference.stat().st_size == r * size,
                            "parity length")
                    with parity.open("rb") as actual, reference.open("rb") as original:
                        while True:
                            block = actual.read(65536)
                            require(block == original.read(65536), "exact Leopard1 parity")
                            if not block:
                                break
                    compared += r * size
    directed_records = {}
    for flavor in ("release", "sanitize"):
        for item in ("kernel", *(str(i) for i in range(14))):
            path = root / "directed-final" / ("%s-%s.log" % (flavor, item))
            peaks.append(resource_peak(path))
            lines = path.read_text().splitlines()
            output = [line for line in lines if line.startswith(("public case ", "predicate:", "GFNI first-stage"))]
            if item == "kernel":
                require(output == ["predicate: exact match and 16 negative neighbors passed",
                    "GFNI first-stage kernel: 224 exact-end/unaligned/scalar cases passed"], "kernel matrix")
            else:
                index = int(item)
                require(len(output) == 1 and output[0].startswith("public case %d:" % index), "directed identity")
                if index in (1, 2):
                    require(output[0] == "public case %d: native odd B%d rejected in both modes before transform" %
                            (index, 65535 if index == 1 else 65537), "odd rejection")
                else:
                    match = re.search(r" calls(\d+) changed(\d+) parity([0-9a-f]{16})$", output[0])
                    require(match is not None, "directed result")
                    calls = 3 if index in (4, 13) else 1 if index in (10, 11) else 2
                    require(int(match[1]) == calls and int(match[2]) == (2 if index in (0, 3, 4) else 0),
                            "directed pass guard")
            if flavor == "release":
                directed_records[item] = output
            else:
                require(output == directed_records[item], "sanitizer/release directed identity")
    require(compared == 122028032, "full parity byte total")
    return {"fixed_public_records": 24, "full_parity_files": 12,
            "compared_parity_bytes": compared, "kernel_cases_per_build": 224,
            "public_cases_per_build": 14, "predicate_negatives_per_build": 16,
            "validated_scopes": len(peaks), "peak_bytes": max(peaks),
            "target_explicit_input_copy_bytes_removed": 65536000, "performance_inference": False}


if __name__ == "__main__":
    require(len(sys.argv) == 4, "usage: verify_gfni_source_stage.py RUN BASELINE COPIES")
    print(json.dumps(replay(*(Path(value) for value in sys.argv[1:])), sort_keys=True))
