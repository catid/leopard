#!/usr/bin/env python3
"""Verify the exact instruction shape of the FF8 XOR3 benchmark probes."""

from __future__ import print_function

import argparse
import collections
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


SCHEMA = "leopard.ff8xor.xor3-assembly.v1"
PROBES = {
    "ff8xor_xor3_probe_forced_xmm": ("forced-xor2", "xmm"),
    "ff8xor_xor3_probe_forced_ymm": ("forced-xor2", "ymm"),
    "ff8xor_xor3_probe_forced_zmm": ("forced-xor2", "zmm"),
    "ff8xor_xor3_probe_explicit_xmm": ("explicit-xor3", "xmm"),
    "ff8xor_xor3_probe_explicit_ymm": ("explicit-xor3", "ymm"),
    "ff8xor_xor3_probe_explicit_zmm": ("explicit-xor3", "zmm"),
}
NM_SYMBOL = re.compile(
    r"^\s*([0-9a-fA-F]+)\s+([0-9a-fA-F]+)\s+[A-Za-z]\s+(\S+)\s*$")
FUNCTION_HEADER = re.compile(r"^\s*[0-9a-fA-F]+ <([^>]+)>:$")
INSTRUCTION = re.compile(
    r"^\s*([0-9a-fA-F]+):\s+([A-Za-z0-9_.]+)(?:\s+(.*))?$")
STACK_REFERENCE = re.compile(r"\([^)]*%(?:r|e)?(?:sp|bp)[^)]*\)")
XOR2_MNEMONICS = {"pxor", "vpxor", "vpxord", "vpxorq"}
XOR3_MNEMONICS = {"vpternlogd", "vpternlogq"}


def find_tool(explicit, candidates):
    if explicit:
        path = shutil.which(explicit) or explicit
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
        raise RuntimeError("tool is not executable: %s" % explicit)
    for candidate in candidates:
        path = shutil.which(candidate)
        if path:
            return path
    raise RuntimeError("none of these tools is available: %s" %
                       ", ".join(candidates))


def run_tool(command):
    try:
        return subprocess.check_output(
            command, stderr=subprocess.STDOUT, universal_newlines=True)
    except subprocess.CalledProcessError as error:
        raise RuntimeError("command failed (%s):\n%s" %
                           (" ".join(command), error.output))


def parse_sizes(text):
    sizes = {}
    for line in text.splitlines():
        match = NM_SYMBOL.match(line)
        if match and match.group(3) in PROBES:
            sizes[match.group(3)] = (
                int(match.group(1), 16), int(match.group(2), 16))
    return sizes


def parse_functions(text):
    functions = {}
    current = None
    for line in text.splitlines():
        header = FUNCTION_HEADER.match(line)
        if header:
            current = header.group(1)
            if current not in PROBES:
                current = None
            else:
                functions[current] = []
            continue
        if current is None:
            continue
        instruction = INSTRUCTION.match(line)
        if instruction:
            functions[current].append((
                int(instruction.group(1), 16),
                instruction.group(2).lower(), instruction.group(3) or ""))
    return functions


def inspect_probe(name, location, instructions):
    operation, width = PROBES[name]
    begin, size = location
    instructions = [item for item in instructions
                    if begin <= item[0] < begin + size]
    mnemonics = collections.Counter()
    xor2 = 0
    xor3 = 0
    bad_xor3_immediates = 0
    wrong_width = 0
    calls = 0
    stack_refs = 0
    width_pattern = re.compile(r"%%%s\d+" % width)

    for unused_address, mnemonic, operands in instructions:
        mnemonics[mnemonic] += 1
        if mnemonic in XOR2_MNEMONICS:
            xor2 += 1
            if not width_pattern.search(operands):
                wrong_width += 1
        if mnemonic in XOR3_MNEMONICS:
            if re.search(r"(?:\$)?0x96\b", operands.lower()):
                xor3 += 1
            else:
                bad_xor3_immediates += 1
            if not width_pattern.search(operands):
                wrong_width += 1
        if mnemonic.startswith("call"):
            calls += 1
        if STACK_REFERENCE.search(operands):
            stack_refs += 1

    violations = []
    if size <= 0:
        violations.append("missing_code_size")
    if not instructions:
        violations.append("missing_disassembly")
    if operation == "forced-xor2":
        if xor2 != 2:
            violations.append("forced_control_requires_exactly_two_xors")
        if xor3 != 0 or bad_xor3_immediates != 0:
            violations.append("forced_control_contains_ternary_logic")
    else:
        if xor3 != 1:
            violations.append("explicit_form_requires_exactly_one_xor3")
        if xor2 != 0:
            violations.append("explicit_form_contains_xor2")
        if bad_xor3_immediates != 0:
            violations.append("vpternlog_immediate_is_not_0x96")
    if wrong_width:
        violations.append("vector_width_mismatch")
    if calls:
        violations.append("probe_contains_call")
    if stack_refs:
        violations.append("probe_contains_stack_reference")

    return {
        "name": name,
        "operation": operation,
        "width": width,
        "code_bytes": size,
        "instruction_count": len(instructions),
        "xor2_count": xor2,
        "xor3_count": xor3,
        "bad_xor3_immediate_count": bad_xor3_immediates,
        "call_count": calls,
        "stack_reference_count": stack_refs,
        "mnemonics": dict(sorted(mnemonics.items())),
        "violations": violations,
    }


def build_report(arguments):
    artifact = arguments.artifact.resolve()
    if not artifact.is_file():
        raise RuntimeError("artifact does not exist: %s" % artifact)
    nm_tool = find_tool(arguments.nm, ("nm", "llvm-nm"))
    objdump_tool = find_tool(arguments.objdump, ("objdump", "llvm-objdump"))
    sizes = parse_sizes(run_tool([
        nm_tool, "-S", "--defined-only", str(artifact)]))
    disassembly = run_tool([
        objdump_tool, "-drw", "--no-show-raw-insn", str(artifact)])
    functions = parse_functions(disassembly)

    probes = []
    for name in sorted(PROBES):
        probes.append(inspect_probe(
            name, sizes.get(name, (0, 0)), functions.get(name, [])))
    violations = []
    for probe in probes:
        for violation in probe["violations"]:
            violations.append({"probe": probe["name"], "rule": violation})
    return {
        "schema": SCHEMA,
        "artifact": artifact.name,
        "tools": {
            "nm": os.path.basename(nm_tool),
            "objdump": os.path.basename(objdump_tool),
        },
        "probes": probes,
        "hard_violation_count": len(violations),
        "hard_violations": violations,
        "strict_pass": not violations,
    }


def human_report(report):
    lines = ["FF8 XOR3 assembly probes: %s" % report["artifact"]]
    for probe in report["probes"]:
        lines.append(
            "  %(width)s %(operation)s: bytes=%(code_bytes)d "
            "instructions=%(instruction_count)d xor2=%(xor2_count)d "
            "xor3=%(xor3_count)d calls=%(call_count)d stack=%(stack_reference_count)d" %
            probe)
    lines.append("  strict=%s violations=%d" % (
        "pass" if report["strict_pass"] else "FAIL",
        report["hard_violation_count"]))
    for violation in report["hard_violations"]:
        lines.append("    %s: %s" %
                     (violation["probe"], violation["rule"]))
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--nm")
    parser.add_argument("--objdump")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--strict", action="store_true")
    arguments = parser.parse_args()
    try:
        report = build_report(arguments)
    except RuntimeError as error:
        print("XOR3 assembly inspection failed: %s" % error, file=sys.stderr)
        return 2
    if arguments.json:
        print(json.dumps(report, sort_keys=True))
    else:
        print(human_report(report))
    return 1 if arguments.strict and not report["strict_pass"] else 0


if __name__ == "__main__":
    sys.exit(main())
