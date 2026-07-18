#!/usr/bin/env python3
"""Inspect retained AVX-512 packed/plane transpose archive members."""

from __future__ import print_function

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


SCHEMA = "leopard.ff8xor.transpose-avx512-assembly.v1"
PROBES = {
    "bitalg": {
        "member": "LeopardFF8XorTransposeBITALG.cpp.o",
        "symbol": "PackedToPlaneBitalg",
        "required": {"vpshufbitqmb": 8, "kmov": 8},
        "allowed_vector": (
            "vmov", "vpaddb", "vpshufbitqmb", "vzeroupper"),
        "max_code_bytes": 1024,
    },
    "vbmi": {
        "member": "LeopardFF8XorTransposeVBMI.cpp.o",
        "symbol": "PlaneToPackedVbmi",
        "required": {"vpermt2b": 24},
        "allowed_vector": (
            "vmov", "vpbroadcast", "vpand", "vpermt2b", "vpsllw",
            "vpsrlw", "vpternlog", "vpxor", "vzeroupper"),
        "max_code_bytes": 4096,
    },
}
FUNCTION_HEADER = re.compile(r"^\s*[0-9a-fA-F]+ <(.+)>:$")
INSTRUCTION = re.compile(
    r"^\s*([0-9a-fA-F]+):\s+([A-Za-z0-9_.]+)(?:\s+(.*))?$")
NM_SYMBOL = re.compile(
    r"^\s*([0-9a-fA-F]+)\s+([0-9a-fA-F]+)\s+[A-Za-z]\s+(\S+)\s*$")
STACK_REFERENCE = re.compile(
    r"(?:\([^)]*%(?:r|e)?(?:sp|bp)[^)]*\)|\[(?:r|e)?(?:sp|bp)[^]]*\])",
    re.IGNORECASE)


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


def run_text(command):
    try:
        return subprocess.check_output(
            command, stderr=subprocess.STDOUT, universal_newlines=True)
    except subprocess.CalledProcessError as error:
        raise RuntimeError("command failed (%s):\n%s" %
                           (" ".join(command), error.output))


def archive_payload(ar_tool, artifact, member):
    try:
        return subprocess.check_output(
            [ar_tool, "p", str(artifact), member], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        raise RuntimeError("could not extract %s: %s" %
                           (member, error.output.decode("utf-8", "replace")))


def inspect_probe(probe, artifact, tools):
    payload = archive_payload(tools["ar"], artifact, probe["member"])
    with tempfile.NamedTemporaryFile(suffix=".o") as temporary:
        temporary.write(payload)
        temporary.flush()
        nm_text = run_text([
            tools["nm"], "-S", "--defined-only", temporary.name])
        disassembly = run_text([
            tools["objdump"], "-drwC", "--no-show-raw-insn",
            temporary.name])

    code_bytes = 0
    for line in nm_text.splitlines():
        match = NM_SYMBOL.match(line)
        if match and probe["symbol"] in match.group(3):
            code_bytes = int(match.group(2), 16)
            break

    instructions = []
    active = False
    for line in disassembly.splitlines():
        header = FUNCTION_HEADER.match(line)
        if header:
            active = probe["symbol"] in header.group(1)
            continue
        if active:
            match = INSTRUCTION.match(line)
            if match:
                instructions.append((
                    match.group(2).lower(), (match.group(3) or "").lower()))

    counts = {}
    calls = 0
    stack_references = 0
    zmm_references = 0
    wrong_width_references = 0
    unexpected_vectors = set()
    for mnemonic, operands in instructions:
        counts[mnemonic] = counts.get(mnemonic, 0) + 1
        if mnemonic.startswith("call"):
            calls += 1
        if STACK_REFERENCE.search(operands):
            stack_references += 1
        if re.search(r"\bzmm\d+\b", operands):
            zmm_references += 1
        if re.search(r"\b[xy]mm\d+\b", operands):
            wrong_width_references += 1
        if (mnemonic.startswith("v") and
                not any(mnemonic.startswith(prefix)
                        for prefix in probe["allowed_vector"])):
            unexpected_vectors.add(mnemonic)

    violations = []
    if not instructions:
        violations.append("missing_disassembly")
    if code_bytes <= 0:
        violations.append("missing_code_size")
    elif code_bytes > probe["max_code_bytes"]:
        violations.append("code_size_exceeds_limit")
    if calls:
        violations.append("call_inside_kernel")
    if stack_references:
        violations.append("stack_spill_or_reference")
    if not zmm_references:
        violations.append("missing_zmm_operations")
    if wrong_width_references:
        violations.append("unexpected_xmm_or_ymm_reference")
    for mnemonic in sorted(unexpected_vectors):
        violations.append("unexpected_vector_instruction_%s" % mnemonic)
    for prefix, minimum in probe["required"].items():
        actual = sum(value for mnemonic, value in counts.items()
                     if mnemonic.startswith(prefix))
        if actual < minimum:
            violations.append("requires_%s_at_least_%d" % (prefix, minimum))

    return {
        "member": probe["member"],
        "symbol": probe["symbol"],
        "code_bytes": code_bytes,
        "instruction_count": len(instructions),
        "call_count": calls,
        "stack_reference_count": stack_references,
        "zmm_reference_count": zmm_references,
        "wrong_width_reference_count": wrong_width_references,
        "unexpected_vector_instructions": sorted(unexpected_vectors),
        "mnemonics": dict(sorted(counts.items())),
        "violations": violations,
    }


def build_report(arguments):
    artifact = arguments.artifact.resolve()
    if not artifact.is_file():
        raise RuntimeError("artifact does not exist: %s" % artifact)
    tools = {
        "ar": find_tool(arguments.ar, ("ar", "llvm-ar")),
        "nm": find_tool(arguments.nm, ("nm", "llvm-nm")),
        "objdump": find_tool(arguments.objdump, ("objdump", "llvm-objdump")),
    }
    probes = {}
    violations = []
    for name in sorted(PROBES):
        probes[name] = inspect_probe(PROBES[name], artifact, tools)
        for violation in probes[name]["violations"]:
            violations.append({"probe": name, "rule": violation})
    return {
        "schema": SCHEMA,
        "artifact": artifact.name,
        "tools": dict((name, os.path.basename(path))
                      for name, path in tools.items()),
        "probes": probes,
        "hard_violation_count": len(violations),
        "hard_violations": violations,
        "strict_pass": not violations,
    }


def human_report(report):
    lines = ["FF8 XOR AVX-512 transpose assembly: %s" % report["artifact"]]
    for name in sorted(report["probes"]):
        probe = report["probes"][name]
        lines.append(
            "  %s: bytes=%d instructions=%d calls=%d stack=%d zmm_refs=%d" %
            (name, probe["code_bytes"], probe["instruction_count"],
             probe["call_count"], probe["stack_reference_count"],
             probe["zmm_reference_count"]))
    lines.append("  strict=%s violations=%d" %
                 ("pass" if report["strict_pass"] else "FAIL",
                  report["hard_violation_count"]))
    for violation in report["hard_violations"]:
        lines.append("    %s: %s" %
                     (violation["probe"], violation["rule"]))
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--ar")
    parser.add_argument("--nm")
    parser.add_argument("--objdump")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--strict", action="store_true")
    arguments = parser.parse_args()
    try:
        report = build_report(arguments)
    except RuntimeError as error:
        print("AVX-512 transpose inspection failed: %s" % error,
              file=sys.stderr)
        return 2
    print(json.dumps(report, sort_keys=True) if arguments.json
          else human_report(report))
    return 1 if arguments.strict and not report["strict_pass"] else 0


if __name__ == "__main__":
    sys.exit(main())
