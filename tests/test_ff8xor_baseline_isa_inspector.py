#!/usr/bin/env python3
"""Finite regression tests for the FF8 XOR baseline-ISA assembly checker."""

from __future__ import print_function

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "tools" / "inspect_ff8xor_baseline_isa.py"
SPEC = importlib.util.spec_from_file_location(
    "inspect_ff8xor_baseline_isa", str(MODULE_PATH))
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def disassembly(symbol, mnemonic, operands=""):
    suffix = (" " + operands) if operands else ""
    return "0000000000000000 <%s>:\n   0: %s%s\n" % (
        symbol, mnemonic, suffix)


def expect_pass(symbol, mnemonic, operands=""):
    report = MODULE.inspect_disassembly(
        disassembly(symbol, mnemonic, operands), "synthetic.o")
    if report["violations"]:
        raise AssertionError("unexpected rejection: %r" % report["violations"])


def expect_failure(symbol, mnemonic, reason, operands=""):
    report = MODULE.inspect_disassembly(
        disassembly(symbol, mnemonic, operands), "synthetic.o")
    observed = [item["reason"] for item in report["violations"]]
    if reason not in observed:
        raise AssertionError(
            "expected %s for %s, observed %r" % (reason, mnemonic, observed))


def main():
    # x86-64-v1 scalar and SSE2 instructions used by the named-register tails.
    expect_pass("portable", "xor", "%rax,%rdx")
    expect_pass("simd128", "pxor", "%xmm1,%xmm0")
    expect_pass("simd128", "pshufd", "$0x4e,%xmm0,%xmm1")
    expect_pass("simd128", "pmuludq", "%xmm1,%xmm0")
    expect_pass("probe", "cpuid")
    expect_pass("leopard::ff8xor::GetX86VectorCapabilities()", "xgetbv")
    expect_pass("leopard::ff8xor::DetectX86VectorCapabilities()", "xgetbv")
    expect_pass("leopard::InitializeCPUArch()", "xgetbv")
    expect_pass("leopard::_xgetbv0()", "xgetbv")

    # Representative SSSE3, SSE4, AVX/AVX-512, GFNI and BMI regressions.
    expect_failure("bad", "pshufb", "post_x86_64_v1_instruction",
                   "%xmm1,%xmm0")
    expect_failure("bad", "addsubps", "post_x86_64_v1_instruction",
                   "%xmm1,%xmm0")
    expect_failure("bad", "haddpd", "post_x86_64_v1_instruction",
                   "%xmm1,%xmm0")
    expect_failure("bad", "lddqu", "post_x86_64_v1_instruction",
                   "(%rax),%xmm0")
    expect_failure("bad", "movddup", "post_x86_64_v1_instruction",
                   "%xmm1,%xmm0")
    expect_failure("bad", "pblendw", "post_x86_64_v1_instruction",
                   "$0x3,%xmm1,%xmm0")
    expect_failure("bad", "pmuldq", "post_x86_64_v1_instruction",
                   "%xmm1,%xmm0")
    expect_failure("bad", "movntdqa", "post_x86_64_v1_instruction",
                   "(%rax),%xmm0")
    expect_failure("bad", "crc32l", "post_x86_64_v1_instruction",
                   "(%rax),%edx")
    expect_failure("bad", "fisttpl", "post_x86_64_v1_instruction",
                   "(%rax)")
    expect_failure("bad", "popcnt", "post_x86_64_v1_instruction",
                   "%rax,%rdx")
    expect_failure("bad", "vpxor", "vex_or_evex_instruction",
                   "%ymm2,%ymm1,%ymm0")
    expect_failure("bad", "vpternlogd", "vex_or_evex_instruction",
                   "$0x96,%zmm2,%zmm1,%zmm0")
    expect_failure("bad", "gf2p8mulb", "post_x86_64_v1_instruction",
                   "%xmm1,%xmm0")
    expect_failure("bad", "pext", "post_x86_64_v1_instruction",
                   "%rax,%rdx,%rcx")
    expect_failure("bad", "mov", "ymm_zmm_or_opmask_register",
                   "%ymm0,(%rax)")
    expect_failure("unrelated", "xgetbv",
                   "xgetbv_outside_guarded_capability_probe")
    expect_failure("leopard::ff8xor::NotGetX86VectorCapabilitiesWrapper()",
                   "xgetbv", "xgetbv_outside_guarded_capability_probe")

    # Parsing must retain function provenance so the XGETBV exception cannot
    # silently expand to every symbol in an object.
    combined = (
        disassembly("leopard::ff8xor::GetX86VectorCapabilities()", "xgetbv") +
        disassembly("payload", "xgetbv"))
    report = MODULE.inspect_disassembly(combined, "synthetic.o")
    if report["violation_count"] != 1:
        raise AssertionError("function provenance regression: %r" % report)

    filtered = MODULE.inspect_disassembly(
        disassembly("baseline", "pxor", "%xmm1,%xmm0") +
        disassembly("optimized", "vpxor", "%ymm2,%ymm1,%ymm0"),
        "synthetic.o", MODULE.re.compile(r"baseline"))
    if filtered["violation_count"] != 0 or \
            filtered["matched_symbols"] != ["baseline"]:
        raise AssertionError("member symbol filter regression: %r" % filtered)

    print("FF8 XOR baseline ISA inspector tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
