#!/usr/bin/env python3
"""Audit retained clock-free qualification; never execute a codec or read a timer."""
import argparse
import json
from pathlib import Path

from check_encode import require, sha256, equal_files
from frozen_inputs import FrozenInputs, TIMING_ARTIFACTS
from timing_contract import validate_timing


def audit(bundle, previous, output):
    reference = previous / "records.json"
    require(sha256(reference) == "dc97ce966bc489d8ac90a1f90b61df346fc817333a518dc26979a9eeb6609687",
            "reference records changed")
    expected = [{key: row[key] for key in ("main", "current")}
                for row in json.loads(reference.read_text())]
    require(json.loads((output / "expected.json").read_text()) == expected, "expected records differ")
    frozen = FrozenInputs(bundle, TIMING_ARTIFACTS)
    records = json.loads((output / "invocations.json").read_text())
    cursor = 0
    counts = dict.fromkeys(("unit", "check", "synthetic", "exercise", "abort", "fault", "cli"), 0)

    def check(name, arguments, category, code=0, fault=None):
        nonlocal cursor
        require(cursor < len(records), "missing invocation")
        row = records[cursor]
        label = f"{cursor:03d}-{category}-{name}"
        require(row == {"name": name, "arguments": arguments, "category": category,
                        "exit_code": code, "stdout": label + ".stdout",
                        "stderr": label + ".stderr", "fault": fault}, "invocation differs: " + label)
        require(name in frozen.pins, "unverified executable")
        out, err = output / row["stdout"], output / row["stderr"]
        require(out.is_file() and err.is_file(), "missing raw output")
        if code == 0:
            require(err.stat().st_size == 0, "successful child stderr")
            result = json.loads(out.read_text())
        else:
            require(out.stat().st_size == 0 and err.stat().st_size > 0, "failure evidence")
            if code == 86:
                require(err.read_text() == "unexpected timing clock in clock-free probe\n", "abort witness")
            result = None
        cursor += 1
        counts[category] += 1
        return result

    for name in ("timing-unit", "timing-unit-sanitized"):
        require(check(name, [], "unit") == {
            "schema": "native-encode-timing-unit/v1", "cases": 26, "real_clocks": 0}, "unit result")
    for implementation in ("main", "current"):
        for kind in ("steady", "synthetic", "abort"):
            for cell in range(8):
                require(check(f"{implementation}-{kind}", ["--check", str(cell)], "check") ==
                        expected[cell][implementation], "clock-free record differs")
        for arguments in ([], ["--synthetic", "0"], ["--measure", "8"],
                          ["--measure", "0", "forbidden"], ["--check", "00"]):
            check(implementation + "-steady", arguments, "cli", 1)
        check(implementation + "-synthetic", ["--measure", "0"], "cli", 1)
        check(implementation + "-abort", ["--measure", "0"], "abort", 86)
        for fault in ("negative", "equal", "reverse", "huge", "unknown"):
            check(implementation + "-synthetic", ["--synthetic", "0"], "fault", 1, fault)
        for cell in (0, 5):
            value = check(implementation + "-abort", ["--exercise", str(cell)], "exercise")
            validate_timing(value, cell, implementation, expected[cell][implementation], "exercise")
        for cell in range(8):
            name = f"{implementation}-{cell}.parity"
            arguments = records[cursor]["arguments"]
            require(len(arguments) == 3 and arguments[:2] == ["--synthetic", str(cell)] and
                    Path(arguments[2]).name == name, "parity output argument")
            value = check(implementation + "-synthetic", arguments, "synthetic")
            validate_timing(value, cell, implementation, expected[cell][implementation], "synthetic")
            equal_files(output / name, previous / f"cell-{cell}-main.parity",
                        expected[cell][implementation]["parity_bytes"])
    require(cursor == len(records) == 94, "qualification invocation inventory")
    require(json.loads((output / "result.json").read_text()) == {
        "complete": True, "counts": counts, "artifact_sha256": frozen.pins,
        "real_benchmark_clocks": 0, "performance_conclusion": None}, "qualification summary")
    frozen.verify()
    return {"audited": True, "invocations": cursor, "counts": counts,
            "manifest_sha256": frozen.manifest_sha256, "new_codec_executions": 0}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--previous-checks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.bundle, args.previous_checks, args.output), indent=2))
