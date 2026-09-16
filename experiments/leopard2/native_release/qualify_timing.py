#!/usr/bin/env python3
"""Qualify timing frontends with synthetic/aborting clocks, never real timings."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from check_encode import sha256, equal_files, require
from timing_contract import validate_timing
from frozen_inputs import FrozenInputs, TIMING_ARTIFACTS


def run(bundle, previous, output):
    expected_path = previous / "records.json"
    require(sha256(expected_path) == "dc97ce966bc489d8ac90a1f90b61df346fc817333a518dc26979a9eeb6609687",
            "previous native qualification records changed")
    expected = [{name: row[name] for name in ("main", "current")}
                for row in json.loads(expected_path.read_text())]
    frozen = FrozenInputs(bundle, TIMING_ARTIFACTS)
    pins = frozen.pins
    verify = frozen.verify

    output.mkdir(mode=0o700)
    verify()
    environment = {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
                   "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "1",
                   "ASAN_OPTIONS": "detect_leaks=1:halt_on_error=1:quarantine_size_mb=8:thread_local_quarantine_size_kb=64",
                   "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1"}
    invocations = []
    counts = {"unit": 0, "check": 0, "synthetic": 0, "exercise": 0, "abort": 0, "fault": 0, "cli": 0}

    def invoke(name, arguments, category, code=0, fault=None):
        label = f"{len(invocations):03d}-{category}-{name}"
        stdout, stderr = (output / (label + suffix) for suffix in (".stdout", ".stderr"))
        env = dict(environment)
        if fault: env["LEO_NATIVE_TEST_CLOCK_FAULT"] = fault
        verify()
        with stdout.open("xb") as out, stderr.open("xb") as err:
            child = subprocess.run([frozen.executable(name), *arguments], env=env,
                                   stdout=out, stderr=err, timeout=180)
        invocations.append({"name": name, "arguments": arguments, "category": category,
                            "exit_code": child.returncode, "stdout": stdout.name,
                            "stderr": stderr.name, "fault": fault})
        (output / "invocations.json").write_text(json.dumps(invocations, indent=2) + "\n")
        verify()
        require(child.returncode == code, "unexpected exit: " + label)
        if code == 0:
            require(stderr.stat().st_size == 0, "unexpected stderr: " + label)
            record = json.loads(stdout.read_text())
        else:
            require(stdout.stat().st_size == 0 and stderr.stat().st_size > 0,
                    "invalid failure evidence: " + label)
            if code == 86:
                require(stderr.read_text() == "unexpected timing clock in clock-free probe\n",
                        "abort guard did not fire")
            record = None
        counts[category] += 1
        return record

    for name in ("timing-unit", "timing-unit-sanitized"):
        require(invoke(name, [], "unit") == {
            "schema": "native-encode-timing-unit/v1", "cases": 26, "real_clocks": 0}, "unit result")
    for implementation in ("main", "current"):
        for kind in ("steady", "synthetic", "abort"):
            for index in range(8):
                actual = invoke(f"{implementation}-{kind}", ["--check", str(index)], "check")
                require(actual == expected[index][implementation], "clock-free check changed")
        for arguments in ([], ["--synthetic", "0"], ["--measure", "8"],
                          ["--measure", "0", "forbidden"], ["--check", "00"]):
            invoke(f"{implementation}-steady", arguments, "cli", 1)
        invoke(f"{implementation}-synthetic", ["--measure", "0"], "cli", 1)
        invoke(f"{implementation}-abort", ["--measure", "0"], "abort", 86)
        for fault in ("negative", "equal", "reverse", "huge", "unknown"):
            invoke(f"{implementation}-synthetic", ["--synthetic", "0"], "fault", 1, fault)
        for index in (0, 5):
            record = invoke(f"{implementation}-abort", ["--exercise", str(index)], "exercise")
            validate_timing(record, index, implementation, expected[index][implementation], "exercise")
        for index in range(8):
            parity = output / f"{implementation}-{index}.parity"
            record = invoke(f"{implementation}-synthetic", ["--synthetic", str(index), str(parity)], "synthetic")
            validate_timing(record, index, implementation, expected[index][implementation], "synthetic")
            equal_files(parity, previous / f"cell-{index}-main.parity", expected[index][implementation]["parity_bytes"])
    verify()
    require(sha256(expected_path) == "dc97ce966bc489d8ac90a1f90b61df346fc817333a518dc26979a9eeb6609687",
            "reference records changed")
    (output / "expected.json").write_text(json.dumps(expected, indent=2) + "\n")
    (output / "result.json").write_text(json.dumps({
        "complete": True, "counts": counts, "artifact_sha256": pins,
        "real_benchmark_clocks": 0, "performance_conclusion": None,
    }, indent=2) + "\n")
    print(json.dumps(counts, sort_keys=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--previous-checks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.bundle.resolve(), args.previous_checks.resolve(), args.output.resolve())
