#!/usr/bin/env python3
"""Unit checks for the plain multi-buffer XOR assembly census."""

from __future__ import print_function

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import inspect_ff8xor_xor_batch as inspector  # noqa: E402


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def rules(result):
    return {item["rule"] for item in result["violations"]}


def test_symbol_classification_and_parsing():
    avx2_symbol = (
        "leopard::ff8xor::avx2::Xor2(void*, void const*, void*, "
        "void const*, unsigned long, unsigned long)")
    avx512_symbol = (
        "leopard::ff8xor::avx512::Xor4_512(void*, void const*, void*, "
        "void const*, void*, void const*, void*, void const*, unsigned long)")
    require(inspector.descriptor(avx2_symbol) ==
            ("avx2_xor2", 2, "ymm"),
            "AVX2 XOR2 symbol was not classified")
    require(inspector.descriptor(avx512_symbol) ==
            ("avx512zmm_xor4", 4, "zmm"),
            "AVX-512 ZMM XOR4 symbol was not classified")
    require(inspector.descriptor("unrelated::Xor2()") is None,
            "unrelated symbol was classified as a batch kernel")

    sizes, symbols = inspector.parse_sizes(
        "00000010 00000020 T %s\n"
        "00000010 00000030 T %s\n" % (avx2_symbol, avx2_symbol))
    require(sizes == {"avx2_xor2": 0x30} and
            symbols["avx2_xor2"] == avx2_symbol,
            "nm parsing did not retain the largest matching body")

    disassembly = (
        "00000010 <%s>:\n"
        "  10: vmovdqu (%%rsi,%%rax,1),%%ymm0\n"
        "  16: vpxor (%%rdx,%%rax,1),%%ymm0,%%ymm0\n"
        "  1c: ret\n" % avx2_symbol)
    functions = inspector.parse_functions(disassembly)
    require(len(functions.get("avx2_xor2", ())) == 3,
            "objdump parsing lost batch-kernel instructions")


def legal_loop(streams, register):
    instructions = []
    address = 0x10
    for stream in range(streams):
        vector = "%%%s%d" % (register, stream)
        instructions.extend((
            (address, "vmovdqu", "(%rsi,%rax,1),{}".format(vector)),
            (address + 4, "vpxor",
             "(%rdx,%rax,1),{},{}".format(vector, vector)),
            (address + 8, "vmovdqu",
             "{},(%rdi,%rax,1)".format(vector)),
        ))
        address += 12
    instructions.extend((
        (address, "add", "$0x20,%rax"),
        (address + 4, "cmp", "%rax,%r11"),
        (address + 8, "jne", "10"),
        (address + 10, "ret", ""),
    ))
    return instructions


def test_legal_named_stream_loops_pass():
    for key, streams, register in (
            ("avx2_xor2", 2, "ymm"),
            ("avx2_xor4", 4, "ymm"),
            ("avx512zmm_xor4", 4, "zmm")):
        result = inspector.analyze(
            key, key, 128, legal_loop(streams, register),
            streams, register)
        require(result["pass"],
                "%s legal loop failed: %r" % (key, result["violations"]))
        require(result["loop_vector_xor_count"] == streams,
                "%s independent XOR count changed" % key)


def test_bad_loop_shapes_are_rejected():
    instructions = legal_loop(2, "ymm")
    # Keep the loop cyclic while adding every structural shape that this small
    # inspector exists to reject.
    instructions[4] = (instructions[4][0], "vpshufb",
                       "%ymm2,%ymm1,%ymm1")
    instructions = instructions[:-4] + [
        (0x25, "vmovdqu", "%ymm0,-0x20(%rsp)"),
        (0x26, "vmovdqu", "0x0(%rip),%ymm3"),
        (0x27, "call", "external_helper"),
    ] + instructions[-4:]
    result = inspector.analyze(
        "avx2_xor2", "bad", 128, instructions, 2, "ymm")
    observed = rules(result)
    require("independent_xors_per_loop" in observed,
            "missing independent XOR was not rejected")
    require("no_calls_in_batch_kernel" in observed and
            "no_calls_in_payload_loop" in observed,
            "payload-loop call was not rejected")
    require("no_vector_stack_spills" in observed,
            "vector stack spill was not rejected")
    require("no_static_table_reads" in observed,
            "RIP-relative loop read was not rejected")
    require("xor_load_store_only_vector_shape" in observed,
            "forbidden shuffle was not rejected")

    wrong_width = inspector.analyze(
        "avx512zmm_xor2", "wrong-width", 128,
        legal_loop(2, "ymm"), 2, "zmm")
    require("required_vector_width_present" in rules(wrong_width) and
            "unexpected_vector_width_absent" in rules(wrong_width),
            "wrong AVX-512 vector width was not rejected")


def main():
    test_symbol_classification_and_parsing()
    test_legal_named_stream_loops_pass()
    test_bad_loop_shapes_are_rejected()
    print("FF8 XOR batch assembly inspector tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
