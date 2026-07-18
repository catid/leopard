#!/usr/bin/env python3
"""Regression checks for sanitizer-aware assembly-test registration."""

from __future__ import print_function

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CMAKE = shutil.which("cmake")
CTEST = shutil.which("ctest")

ARTIFACT_TESTS = (
    "ff8xor_baseline_isa_census",
    "ff8xor_baseline_isa_custom_target",
    "ff8xor_xor3_assembly",
    "ff8xor_xor_batch_assembly",
    "ff8xor_assembly_census",
)
ARTIFACT_TARGETS = (
    "check_ff8xor_baseline_isa",
    "check_ff8xor_xor3_assembly",
    "inspect_ff8xor_assembly",
    "check_ff8xor_assembly",
    "check_ff8xor_xor_batch_assembly",
)
PARSER_TESTS = (
    "ff8xor_assembly_inspector",
    "ff8xor_xor_batch_inspector",
    "ff8xor_baseline_isa_inspector",
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def run(command):
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True)
    if result.returncode != 0:
        raise AssertionError(
            "command failed (%d): %s\n%s" %
            (result.returncode, " ".join(command), result.stdout))
    return result.stdout


def check_configuration(label, generator, flags, configuration=None):
    with tempfile.TemporaryDirectory(
            prefix="leopard-ff8xor-sanitizer-gate-") as temporary:
        build = Path(temporary) / "build"
        output = run([
            CMAKE, "-S", str(ROOT), "-B", str(build),
            "-G", generator, "-DENABLE_OPENMP=OFF",
        ] + list(flags))
        require("FF8 XOR production assembly checks: skipped for "
                "sanitizer-instrumented objects" in output,
                "%s configure did not report sanitizer gating" % label)

        ctest_command = [CTEST, "--test-dir", str(build), "-N"]
        if configuration:
            ctest_command.extend(("-C", configuration))
        tests = run(ctest_command)
        for name in ARTIFACT_TESTS:
            require(name not in tests,
                    "%s retained sanitizer artifact test %s" %
                    (label, name))
        for name in PARSER_TESTS:
            require(name in tests,
                    "%s lost sanitizer-safe parser test %s" %
                    (label, name))

        targets = run([CMAKE, "--build", str(build), "--target", "help"])
        for name in ARTIFACT_TARGETS:
            require(name not in targets,
                    "%s retained sanitizer artifact target %s" %
                    (label, name))


def main():
    require(CMAKE is not None and CTEST is not None,
            "cmake and ctest are required")
    help_text = run([CMAKE, "--help"])
    ninja = shutil.which("ninja")
    single_generator = "Ninja" if ninja else "Unix Makefiles"
    check_configuration(
        "single-config Release", single_generator,
        ("-DCMAKE_BUILD_TYPE=Release",
         "-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined"))

    if ninja and "Ninja Multi-Config" in help_text:
        check_configuration(
            "multi-config Debug", "Ninja Multi-Config",
            ("-DCMAKE_CXX_FLAGS_DEBUG=-fsanitize=address,undefined",),
            "Debug")

    print("FF8 XOR sanitizer assembly-gating tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
