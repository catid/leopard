#!/usr/bin/env python3
"""Inspect the static two/four-stream FF8 XOR loops.

The generated field circuits have a separate assembly census.  This smaller
gate covers the plain-buffer accumulation kernels used by encoder chunking and
the formal derivative: each required function must contain a vector loop with
the requested number of independent XORs, no calls, no vector spills, no
dynamic stack indexing, and no lookup/shuffle/multiply instructions.
"""

from __future__ import print_function

import argparse
import json
import os
import re
import sys
from pathlib import Path

import inspect_ff8xor_assembly as core


SCHEMA = "leopard.ff8xor.xor-batch-assembly.v1"
AVX2_NAME = re.compile(r"leopard::ff8xor::avx2::Xor([24])\(")
AVX512_NAME = re.compile(
    r"leopard::ff8xor::avx512::Xor([24])_(256|512)\(")


def descriptor(symbol):
    match = AVX2_NAME.search(symbol)
    if match:
        return "avx2_xor%s" % match.group(1), int(match.group(1)), "ymm"
    match = AVX512_NAME.search(symbol)
    if match:
        width = "vl" if match.group(2) == "256" else "zmm"
        register = "ymm" if width == "vl" else "zmm"
        return ("avx512%s_xor%s" % (width, match.group(1)),
                int(match.group(1)), register)
    return None


def parse_sizes(text):
    sizes = {}
    symbols = {}
    for line in text.splitlines():
        match = core.NM_SYMBOL.match(line)
        if not match:
            continue
        item = descriptor(match.group(3))
        if item is None:
            continue
        key = item[0]
        size = int(match.group(2), 16)
        if key not in sizes or size > sizes[key]:
            sizes[key] = size
            symbols[key] = match.group(3)
    return sizes, symbols


def parse_functions(text):
    functions = {}
    current = None
    for line in text.splitlines():
        header = core.FUNCTION_HEADER.match(line)
        if header:
            item = descriptor(header.group(1))
            current = item[0] if item else None
            if current is not None:
                functions.setdefault(current, [])
            continue
        if current is None:
            continue
        match = core.INSTRUCTION.match(line)
        if match:
            functions[current].append((
                int(match.group(1), 16),
                match.group(2).lower(),
                (match.group(3) or "").lower()))
    return functions


def is_vector_xor(mnemonic, operands):
    return mnemonic in core.XOR2_MNEMONICS and bool(
        core.VECTOR_REGISTER.search(operands))


def analyze(key, symbol, size, instructions, streams, required_register):
    cycles = core.cyclic_instruction_components(instructions)
    loop_addresses = set().union(*cycles) if cycles else set()
    loop = [item for item in instructions if item[0] in loop_addresses]
    vector_xors = sum(is_vector_xor(mnemonic, operands)
                      for unused, mnemonic, operands in instructions)
    loop_vector_xors = sum(is_vector_xor(mnemonic, operands)
                           for unused, mnemonic, operands in loop)
    calls = sum(mnemonic.startswith("call")
                for unused, mnemonic, operands in instructions)
    inner_calls = sum(mnemonic.startswith("call")
                      for unused, mnemonic, operands in loop)
    vector_stack_refs = 0
    scaled_stack_refs = 0
    loop_rip_refs = 0
    forbidden = []
    for address, mnemonic, operands in instructions:
        operands_list = core.split_operands(operands)
        if core.VECTOR_REGISTER.search(operands):
            if any(core.is_stack_memory(operand, True)
                   for operand in operands_list):
                vector_stack_refs += 1
        for operand in operands_list:
            parts = core.memory_address_parts(operand)
            if parts is not None:
                base, index, scale = parts
                if base in ("%rsp", "%rbp") and index:
                    scaled_stack_refs += 1
        if address in loop_addresses:
            if "(%rip)" in operands:
                loop_rip_refs += 1
            if (mnemonic in core.FORBIDDEN_EXACT or
                    any(mnemonic.startswith(prefix)
                        for prefix in core.FORBIDDEN_PREFIXES)):
                forbidden.append(mnemonic)

    required_register_pattern = re.compile(
        r"%%%s\d+" % required_register)
    required_width_in_loop = any(required_register_pattern.search(operands)
                                 for unused, mnemonic, operands in loop)
    wrong_width_pattern = (re.compile(r"%zmm\d+")
                           if required_register == "ymm"
                           else re.compile(r"%[xy]mm\d+"))
    wrong_width_in_loop = any(wrong_width_pattern.search(operands)
                              for unused, mnemonic, operands in loop)

    violations = []
    checks = (
        (bool(cycles), "vector_payload_loop_present", len(cycles)),
        (loop_vector_xors >= streams,
         "independent_xors_per_loop", loop_vector_xors),
        (calls == 0, "no_calls_in_batch_kernel", calls),
        (inner_calls == 0, "no_calls_in_payload_loop", inner_calls),
        (vector_stack_refs == 0, "no_vector_stack_spills", vector_stack_refs),
        (scaled_stack_refs == 0,
         "no_dynamic_stack_array_indexing", scaled_stack_refs),
        (loop_rip_refs == 0, "no_static_table_reads", loop_rip_refs),
        (not forbidden, "xor_load_store_only_vector_shape", forbidden),
        (required_width_in_loop, "required_vector_width_present", False),
        (not wrong_width_in_loop, "unexpected_vector_width_absent", True),
    )
    for passed, rule, observed in checks:
        if not passed:
            violations.append({"rule": rule, "observed": observed})

    return {
        "key": key,
        "symbol": symbol,
        "code_bytes": size,
        "streams": streams,
        "required_register": required_register,
        "instruction_count": len(instructions),
        "loop_component_count": len(cycles),
        "vector_xor_count": vector_xors,
        "loop_vector_xor_count": loop_vector_xors,
        "call_count": calls,
        "inner_call_count": inner_calls,
        "vector_stack_ref_count": vector_stack_refs,
        "scaled_stack_ref_count": scaled_stack_refs,
        "loop_rip_ref_count": loop_rip_refs,
        "forbidden_mnemonics": sorted(set(forbidden)),
        "violations": violations,
        "pass": not violations,
    }


def build_report(arguments):
    artifact = arguments.artifact.resolve()
    if not artifact.is_file():
        raise RuntimeError("artifact does not exist: %s" % artifact)
    nm_tool = core.find_tool(arguments.nm, ("nm", "llvm-nm"))
    objdump_tool = core.find_tool(
        arguments.objdump, ("objdump", "llvm-objdump"))
    sizes, symbols = parse_sizes(core.run_tool([
        nm_tool, "-S", "-C", "--defined-only", str(artifact)]))
    functions = parse_functions(core.run_tool([
        objdump_tool, "-drwC", "--no-show-raw-insn", str(artifact)]))

    required = []
    if arguments.require_avx2:
        required.extend(("avx2_xor2", "avx2_xor4"))
    if arguments.require_avx512:
        required.extend((
            "avx512vl_xor2", "avx512vl_xor4",
            "avx512zmm_xor2", "avx512zmm_xor4"))
    if not required:
        required = sorted(set(sizes) & set(functions))

    results = []
    violations = []
    for key in required:
        item = descriptor(symbols[key]) if key in symbols else None
        if key not in sizes or key not in functions or item is None:
            violations.append({
                "key": key,
                "rule": "required_batch_kernel_present",
                "observed": False,
            })
            continue
        result = analyze(
            key, symbols[key], sizes[key], functions[key], item[1], item[2])
        results.append(result)
        for violation in result["violations"]:
            record = dict(violation)
            record["key"] = key
            violations.append(record)

    return {
        "schema": SCHEMA,
        "artifact": artifact.name,
        "tools": {
            "nm": os.path.basename(nm_tool),
            "objdump": os.path.basename(objdump_tool),
        },
        "required": required,
        "kernels": results,
        "violation_count": len(violations),
        "violations": violations,
        "strict_pass": not violations,
    }


def human_report(report):
    lines = [
        "FF8 XOR batch assembly census: %s" % report["artifact"],
        "strict batch-loop checks: %s" %
        ("PASS" if report["strict_pass"] else "FAIL"),
    ]
    for item in report["kernels"]:
        lines.append(
            "%s: %d bytes, %d instructions, loop-xor=%d, calls=%d "
            "(inner=%d), vector-stack=%d, forbidden=%s" % (
                item["key"], item["code_bytes"], item["instruction_count"],
                item["loop_vector_xor_count"], item["call_count"],
                item["inner_call_count"], item["vector_stack_ref_count"],
                ",".join(item["forbidden_mnemonics"]) or "none"))
    for violation in report["violations"]:
        lines.append("  %s" % json.dumps(violation, sort_keys=True))
    return "\n".join(lines) + "\n"


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--nm")
    parser.add_argument("--objdump")
    parser.add_argument("--require-avx2", action="store_true")
    parser.add_argument("--require-avx512", action="store_true")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    try:
        report = build_report(arguments)
    except (KeyError, RuntimeError) as error:
        print("batch assembly census unavailable: %s" % error,
              file=sys.stderr)
        return 2
    output = (json.dumps(report, sort_keys=True, indent=2) + "\n"
              if arguments.json else human_report(report))
    sys.stdout.write(output)
    return 1 if arguments.strict and not report["strict_pass"] else 0


if __name__ == "__main__":
    sys.exit(main())
