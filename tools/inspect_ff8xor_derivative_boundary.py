#!/usr/bin/env python3
"""Inspect the production AVX-512 formal-derivative boundary row loops.

The source-arity switch is metadata dispatch performed once before one of
eight statically named payload loops.  A compiler may lower it either as a
comparison tree or as one pre-loop jump table.  This inspection permits both
shapes, then rejects calls, any payload-time indirect jump, vector stack
operands, non-XOR ternary immediates, and field-multiplication opcodes.
"""

from __future__ import print_function

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import inspect_ff8xor_four_buffer_production as common  # noqa: E402


OBJECT_BASENAME = "LeopardFF8XorAVX512.cpp.o"
SYMBOL_RE = re.compile(
    r"derivative_detail::ApplyRow<leopard::ff8xor::avx512::"
    r"(Ymm|Zmm)DerivativeOps>")
INDIRECT_JUMP_RE = re.compile(r"(?:^|\s)jmp\s+\*", re.IGNORECASE)


def select_archive_member(member_names):
    exact = [name for name in member_names
             if Path(name).name == OBJECT_BASENAME]
    if len(exact) == 1:
        return exact[0]
    fallback = [name for name in member_names
                if "LeopardFF8XorAVX512.cpp" in Path(name).name and
                name.endswith((".o", ".obj"))]
    if len(fallback) == 1:
        return fallback[0]
    raise RuntimeError("archive does not contain one unambiguous %s" %
                       OBJECT_BASENAME)


def extract_member(ar, archive, destination):
    members = [line.strip() for line in common.run_text(
        [ar, "t", str(archive)]).splitlines() if line.strip()]
    member = select_archive_member(members)
    result = subprocess.run(
        [ar, "p", str(archive), member], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    destination.write_bytes(result.stdout)
    return member


def inspect_body(body):
    records = common.instruction_records(body)
    rbp_frame = common.establishes_rbp_frame(records)
    first_vector = None
    indirect_jumps = []
    calls = []
    vector_stack = []
    bad_ternary = []
    forbidden = []
    ternary_count = 0
    vector_count = 0

    for index, record in enumerate(records):
        line = record["line"]
        mnemonic = record["mnemonic"]
        operands = record["operands"]
        has_vector = common.VECTOR_RE.search(operands) is not None
        if has_vector:
            vector_count += 1
            if first_vector is None:
                first_vector = index
        if INDIRECT_JUMP_RE.search("%s %s" % (mnemonic, operands)):
            indirect_jumps.append(index)
        if mnemonic.startswith("call"):
            calls.append(line)
        if has_vector and common.RSP_MEMORY_RE.search(operands):
            vector_stack.append(line)
        if (has_vector and rbp_frame and
                common.RBP_MEMORY_RE.search(operands)):
            vector_stack.append(line)
        if mnemonic.startswith("vpternlog"):
            ternary_count += 1
            if not common.ternary_immediate_is_xor3(operands):
                bad_ternary.append(line)
        if ("pshufb" in mnemonic or "gf2p8" in mnemonic or
                "pclmul" in mnemonic or mnemonic.startswith("vpmull") or
                mnemonic.startswith(("vpand", "pand")) or
                mnemonic in ("andps", "andpd", "andnps", "andnpd",
                             "vandps", "vandpd")):
            forbidden.append(line)

    failures = []
    if len(indirect_jumps) > 1:
        failures.append("expected at most one source-arity indirect jump, "
                        "found %d" % len(indirect_jumps))
    elif (len(indirect_jumps) == 1 and first_vector is not None and
          indirect_jumps[0] >= first_vector):
        failures.append("source-arity indirect jump is inside payload code")
    if calls:
        failures.append("payload implementation contains calls")
    if vector_stack:
        failures.append("payload vectors reference the stack")
    if bad_ternary:
        failures.append("VPTERNLOG immediate is not XOR3 0x96")
    if forbidden:
        failures.append("forbidden field-multiplication opcode found")
    if ternary_count == 0:
        failures.append("no VPTERNLOG XOR3 lowering found")
    if vector_count == 0:
        failures.append("no vector payload instructions found")
    return {
        "instruction_count": len(records),
        "vector_instruction_count": vector_count,
        "vpternlog_xor3_count": ternary_count - len(bad_ternary),
        "indirect_jump_count": len(indirect_jumps),
        "indirect_jump_before_payload": bool(
            len(indirect_jumps) == 1 and first_vector is not None and
            indirect_jumps[0] < first_vector),
        "call_count": len(calls),
        "vector_stack_reference_count": len(vector_stack),
        "bad_ternary_count": len(bad_ternary),
        "forbidden_instruction_count": len(forbidden),
        "failures": failures,
    }


def inspect_object(path, objdump):
    output = common.run_text([
        objdump, "-d", "-C", "--no-show-raw-insn", str(path)])
    bodies = common.parse_disassembly(output)
    observations = {}
    failures = []
    for symbol, body in bodies.items():
        match = SYMBOL_RE.search(symbol)
        if not match:
            continue
        width = match.group(1).lower()
        if width in observations:
            failures.append("duplicate %s ApplyRow specialization" % width)
            continue
        observation = inspect_body(body)
        observation["symbol"] = symbol
        observations[width] = observation
        failures.extend("%s: %s" % (width, failure)
                        for failure in observation["failures"])
    for width in ("ymm", "zmm"):
        if width not in observations:
            failures.append("missing %s ApplyRow specialization" % width)
    return observations, failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True)
    parser.add_argument("--ar")
    parser.add_argument("--objdump")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    archive = Path(args.archive).resolve()
    ar = args.ar or shutil.which("llvm-ar") or shutil.which("ar")
    objdump = (args.objdump or shutil.which("llvm-objdump") or
               shutil.which("objdump"))
    if not archive.is_file() or not ar or not objdump:
        print("archive, ar, and objdump are required", file=sys.stderr)
        return 2

    try:
        with tempfile.TemporaryDirectory(
                prefix="ff8xor-derivative-boundary-") as directory:
            path = Path(directory) / OBJECT_BASENAME
            member = extract_member(ar, archive, path)
            observations, failures = inspect_object(path, objdump)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print("derivative-boundary inspection failed: %s" % error,
              file=sys.stderr)
        return 2

    report = {
        "archive_member": member,
        "status": "clean" if not failures else "failed",
        "observations": observations,
        "failures": failures,
    }
    if args.json:
        print(json.dumps(report, sort_keys=True, indent=2))
    else:
        print("derivative boundary member: %s" % member)
        for width in sorted(observations):
            item = observations[width]
            print("  %s: instructions=%d vectors=%d xor3=%d "
                  "dispatch=%d calls=%d spills=%d" % (
                      width,
                      item["instruction_count"],
                      item["vector_instruction_count"],
                      item["vpternlog_xor3_count"],
                      item["indirect_jump_count"],
                      item["call_count"],
                      item["vector_stack_reference_count"]))
        for failure in failures:
            print("FAIL: %s" % failure, file=sys.stderr)
    return 1 if args.strict and failures else 0


if __name__ == "__main__":
    sys.exit(main())
