#!/usr/bin/env python3
"""Unit checks for the FF8 XOR assembly-census regression rules."""

from __future__ import print_function

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import inspect_ff8xor_assembly as inspector  # noqa: E402


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def test_required_families():
    violations = inspector.hard_violations(
        [], {}, inspector.BASE_FAMILIES)
    observed = {item["family"] for item in violations
                if item["rule"] == "all_256_specializations_present"}
    require(observed == set(inspector.BASE_FAMILIES),
            "empty input did not fail every mandatory base family")

    summary = {
        "multiply": {"coefficient_count": 256},
        "fft": {"coefficient_count": 255},
        "ifft": {"coefficient_count": 256},
    }
    violations = inspector.hard_violations(
        [], summary, inspector.BASE_FAMILIES)
    require(len(violations) == 1 and violations[0]["family"] == "fft",
            "partial specialization family was not rejected exactly")


def test_payload_lookup_and_forbidden_shapes():
    # The backwards branch makes addresses 0x10..0x24 a cyclic payload block.
    # The table base is formed once outside the loop, copied, and then used by
    # a narrow SIB load.  The loop also contains direct RIP data and a vector
    # permute that a static named-register XOR circuit must never need.
    instructions = [
        (0x00, "lea", "0x0(%rip),%rax"),
        (0x07, "mov", "%rax,%r9"),
        (0x0a, "jmp", "10"),
        (0x10, "movzbl", "(%r9,%rcx,1),%edx"),
        (0x14, "movq", "0x0(%rip),%r8"),
        (0x1b, "vpermd", "%ymm0,%ymm1,%ymm2"),
        (0x20, "add", "$0x1,%rcx"),
        (0x24, "jne", "10"),
        (0x26, "ret", ""),
    ]
    census = inspector.instruction_census(
        "void MultiplyWholeBuffer<0u>", 64, instructions)
    counts = census["counts"]
    require(counts["loop_narrow_indexed_loads"] == 1,
            "narrow indexed lookup was not detected")
    require(counts["static_table_indexed_refs"] == 1,
            "RIP-derived static table base was not detected")
    require(counts["loop_rip_memory_refs"] == 1,
            "direct RIP-relative loop access was not detected")
    require(counts["forbidden_instructions"] == 1,
            "vector permute was not detected")


def test_allowed_xor3():
    instructions = [
        (0x00, "vpternlogd", "$0x96,%zmm1,%zmm2,%zmm0"),
        (0x07, "ret", ""),
    ]
    census = inspector.instruction_census(
        "void MultiplyWholeBuffer<1u>", 8, instructions)
    require(census["counts"]["vector_xor3"] == 1,
            "XOR3 ternary instruction was not counted")
    require(census["counts"]["non_xor_ternary"] == 0,
            "valid XOR3 immediate was rejected")


def main():
    test_required_families()
    test_payload_lookup_and_forbidden_shapes()
    test_allowed_xor3()
    print("FF8 XOR assembly inspector tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
