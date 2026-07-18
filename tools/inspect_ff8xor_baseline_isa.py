#!/usr/bin/env python3
"""Reject post-baseline instructions in FF8 XOR dispatcher/portable objects.

The optimized AVX2 and AVX-512 circuit kernels intentionally live in separate
archive members.  This checker extracts only the public API gateways and the
portable/128-bit dispatcher member, then verifies that their disassembly stays
within the x86-64-v1/SSE2 baseline.  XGETBV is the sole intentional exception:
it is permitted only in the CPUID/OSXSAVE-gated capability probe.
"""

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


SCHEMA = "leopard.ff8xor.baseline-isa-census.v1"
FUNCTION_HEADER = re.compile(r"^\s*[0-9a-fA-F]+ <(.+)>:$")
INSTRUCTION = re.compile(
    r"^\s*([0-9a-fA-F]+):\s+([A-Za-z0-9_.]+)(?:\s+(.*))?$")
WIDE_REGISTER = re.compile(r"%(?:[yz]mm\d+|k[0-7])\b", re.IGNORECASE)

# Instructions introduced after the x86-64-v1/SSE2 baseline.  VEX/EVEX forms
# are caught separately by their leading "v" mnemonic.  Keep this explicit so
# a failure identifies the exact extension class instead of relying on byte
# prefix heuristics that can mistake inline data for instructions.
FORBIDDEN_EXACT = {
    # SSE3 (x86-64-v2, not the v1/SSE2 baseline)
    "addsubpd", "addsubps", "fisttp", "haddpd", "haddps", "hsubpd",
    "hsubps", "lddqu", "movddup", "movshdup", "movsldup", "monitor",
    "mwait",
    # SSSE3
    "pabsb", "pabsw", "pabsd", "palignr", "phaddw", "phaddd",
    "phaddsw", "phsubw", "phsubd", "phsubsw", "pmaddubsw",
    "pmulhrsw", "pshufb", "psignb", "psignw", "psignd",
    # SSE4.x and POPCNT
    "blendpd", "blendps", "blendvpd", "blendvps", "crc32", "dppd",
    "dpps", "extractps", "insertps", "mpsadbw", "packusdw", "pblendvb",
    "movntdqa", "pblendw", "pcmpeqq", "pcmpestri", "pcmpestrm",
    "pcmpgtq", "pcmpistri", "pcmpistrm", "pextrb", "pextrd", "pextrq",
    "phminposuw", "pinsrb",
    "pinsrd", "pinsrq", "pmaxsb", "pmaxsd", "pmaxud", "pmaxuw",
    "pminsb", "pminsd", "pminud", "pminuw", "pmovsxbd", "pmovsxbq",
    "pmovsxbw", "pmovsxdq", "pmovsxwd", "pmovsxwq", "pmovzxbd",
    "pmovzxbq", "pmovzxbw", "pmovzxdq", "pmovzxwd", "pmovzxwq",
    "pmuldq", "pmulld", "pmullq", "ptest", "roundpd", "roundps", "roundsd",
    "roundss", "popcnt",
    # AMD SSE4a
    "extrq", "insertq", "movntsd", "movntss",
    # AES, carry-less multiply, SHA, GFNI
    "aesdec", "aesdeclast", "aesenc", "aesenclast", "aesimc",
    "aeskeygenassist", "pclmulqdq", "sha1msg1", "sha1msg2", "sha1nexte",
    "sha1rnds4", "sha256msg1", "sha256msg2", "sha256rnds2",
    "gf2p8affineinvqb", "gf2p8affineqb", "gf2p8mulb",
    # BMI, BMI2, ADX, LZCNT and later scalar extensions
    "adcx", "adox", "andn", "bextr", "blsi", "blsmsk", "blsr", "bzhi",
    "cmpxchg16b", "lahf", "lzcnt", "movbe", "mulx", "pdep", "pext",
    "prefetchw", "rdfsbase", "rdgsbase", "rdpid", "rdrand", "rdseed",
    "rdtscp", "rorx", "sahf", "sarx", "shlx", "shrx", "tzcnt",
    "wrfsbase", "wrgsbase", "xrstor", "xrstor64", "xrstors", "xrstors64",
    "xsave", "xsave64", "xsavec", "xsavec64", "xsaveopt", "xsaveopt64",
    "xsaves", "xsaves64", "xsetbv",
}
FORBIDDEN_PREFIXES = (
    # GNU objdump appends operand-size suffixes to these mnemonics.
    "crc32", "fisttp",
    "clflushopt", "clwb", "enqcmd", "movdir", "prefetchwt1", "serialize",
    "tile", "tpause", "umonitor", "umwait",
)
V_PREFIX_BASELINE_EXCEPTIONS = {"verr", "verw"}
GUARDED_XGETBV_PROBE = re.compile(
    r"(?:^|::)(?:DetectX86VectorCapabilities|GetX86VectorCapabilities|"
    r"DetectAVX512TransposeCapabilities|InitializeCPUArch|_xgetbv0)"
    r"\([^)]*\)(?: \[clone[^]]*\])?$"
)


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


def parse_disassembly(text):
    current = "<unknown>"
    instructions = []
    for line in text.splitlines():
        header = FUNCTION_HEADER.match(line)
        if header:
            current = header.group(1)
            continue
        match = INSTRUCTION.match(line)
        if not match:
            continue
        instructions.append({
            "address": int(match.group(1), 16),
            "symbol": current,
            "mnemonic": match.group(2).lower(),
            "operands": (match.group(3) or "").lower(),
        })
    return instructions


def classify_instruction(instruction):
    mnemonic = instruction["mnemonic"]
    operands = instruction["operands"]
    symbol = instruction["symbol"]
    if mnemonic.startswith("v") and mnemonic not in V_PREFIX_BASELINE_EXCEPTIONS:
        return "vex_or_evex_instruction"
    if WIDE_REGISTER.search(operands):
        return "ymm_zmm_or_opmask_register"
    if mnemonic in FORBIDDEN_EXACT:
        return "post_x86_64_v1_instruction"
    if any(mnemonic.startswith(prefix) for prefix in FORBIDDEN_PREFIXES):
        return "post_x86_64_v1_instruction"
    if mnemonic == "xgetbv" and not GUARDED_XGETBV_PROBE.search(symbol):
        return "xgetbv_outside_guarded_capability_probe"
    return None


def inspect_disassembly(text, member, symbol_filter=None):
    instructions = parse_disassembly(text)
    if symbol_filter is not None:
        instructions = [instruction for instruction in instructions
                        if symbol_filter.search(instruction["symbol"])]
    violations = []
    xgetbv_count = 0
    cpuid_count = 0
    for instruction in instructions:
        if instruction["mnemonic"] == "xgetbv":
            xgetbv_count += 1
        elif instruction["mnemonic"] == "cpuid":
            cpuid_count += 1
        reason = classify_instruction(instruction)
        if reason:
            item = dict(instruction)
            item["member"] = member
            item["reason"] = reason
            item["address"] = "0x%x" % item["address"]
            violations.append(item)
    return {
        "member": member,
        "symbol_filter": symbol_filter.pattern if symbol_filter else None,
        "matched_symbols": sorted(set(
            instruction["symbol"] for instruction in instructions)),
        "instruction_count": len(instructions),
        "cpuid_count": cpuid_count,
        "guarded_xgetbv_count": xgetbv_count,
        "violation_count": len(violations),
        "violations": violations,
    }


def archive_members(ar_tool, artifact):
    return [line.strip() for line in
            run_text([ar_tool, "t", str(artifact)]).splitlines()
            if line.strip()]


def disassemble_member(ar_tool, objdump_tool, artifact, member):
    try:
        payload = subprocess.check_output(
            [ar_tool, "p", str(artifact), member], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        raise RuntimeError("could not extract %s from %s: %s" %
                           (member, artifact, error.output.decode(
                               "utf-8", "replace")))
    suffix = os.path.splitext(member)[1] or ".o"
    with tempfile.NamedTemporaryFile(suffix=suffix) as temporary:
        temporary.write(payload)
        temporary.flush()
        return run_text([
            objdump_tool, "-dC", "--no-show-raw-insn", temporary.name])


def build_report(arguments):
    artifact = arguments.artifact.resolve()
    if not artifact.is_file():
        raise RuntimeError("artifact does not exist: %s" % artifact)
    objdump_tool = find_tool(arguments.objdump, ("objdump", "llvm-objdump"))
    reports = []
    filters = {}
    for specification in arguments.member_symbol:
        if "=" not in specification:
            raise RuntimeError(
                "--member-symbol requires MEMBER=REGEX: %s" % specification)
        member_name, pattern = specification.split("=", 1)
        if not member_name or not pattern:
            raise RuntimeError(
                "--member-symbol requires MEMBER=REGEX: %s" % specification)
        try:
            filters[member_name] = re.compile(pattern)
        except re.error as error:
            raise RuntimeError("invalid symbol regex %s: %s" %
                               (pattern, error))
    if arguments.member:
        ar_tool = find_tool(arguments.ar, ("ar", "llvm-ar"))
        available = archive_members(ar_tool, artifact)
        for requested in arguments.member:
            matches = [name for name in available
                       if name == requested or os.path.basename(name) == requested]
            if len(matches) != 1:
                raise RuntimeError(
                    "archive member %s matched %d entries (available: %s)" %
                    (requested, len(matches), ", ".join(available)))
            text = disassemble_member(
                ar_tool, objdump_tool, artifact, matches[0])
            symbol_filter = filters.get(requested, filters.get(matches[0]))
            report = inspect_disassembly(text, matches[0], symbol_filter)
            if symbol_filter is not None and not report["matched_symbols"]:
                raise RuntimeError(
                    "symbol filter %s matched nothing in %s" %
                    (symbol_filter.pattern, matches[0]))
            reports.append(report)
    else:
        text = run_text([
            objdump_tool, "-dC", "--no-show-raw-insn", str(artifact)])
        reports.append(inspect_disassembly(text, artifact.name))
    violations = [item for report in reports for item in report["violations"]]
    return {
        "schema": SCHEMA,
        "artifact": artifact.name,
        "members": reports,
        "violation_count": len(violations),
        "violations": violations,
        "strict_pass": not violations,
    }


def human_report(report):
    lines = [
        "FF8 XOR baseline ISA census: %s" % report["artifact"],
        "x86-64-v1/SSE2 check: %s" %
        ("PASS" if report["strict_pass"] else "FAIL"),
    ]
    for member in report["members"]:
        lines.append(
            "%s%s: instructions=%d cpuid=%d guarded-xgetbv=%d violations=%d" %
            (member["member"],
             (" [%s]" % member["symbol_filter"]
              if member["symbol_filter"] else ""),
             member["instruction_count"],
             member["cpuid_count"], member["guarded_xgetbv_count"],
             member["violation_count"]))
    for violation in report["violations"][:20]:
        lines.append("  %s" % json.dumps(violation, sort_keys=True))
    return "\n".join(lines) + "\n"


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path,
                        help="object file or static Leopard archive")
    parser.add_argument("--member", action="append", default=[],
                        help="archive member to inspect (repeatable)")
    parser.add_argument(
        "--member-symbol", action="append", default=[],
        help="limit one archive member to demangled symbols: MEMBER=REGEX")
    parser.add_argument("--ar", help="ar/llvm-ar executable")
    parser.add_argument("--objdump", help="objdump/llvm-objdump executable")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--strict", action="store_true",
                        help="exit nonzero when a forbidden instruction exists")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    try:
        report = build_report(arguments)
    except RuntimeError as error:
        print("baseline ISA census unavailable: %s" % error, file=sys.stderr)
        return 2
    output = (json.dumps(report, sort_keys=True, indent=2) + "\n"
              if arguments.json else human_report(report))
    sys.stdout.write(output)
    return 1 if arguments.strict and not report["strict_pass"] else 0


if __name__ == "__main__":
    sys.exit(main())
