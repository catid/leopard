#!/usr/bin/env python3
"""Create a machine-readable census of generated FF8 XOR kernels.

The inspector accepts either LeopardFF8Xor.cpp.o or the static Leopard archive.
It recognizes every coefficient-specialized whole-buffer kernel and reports
instruction shape, code size, calls, stack/vector spills, scaled stack indexing,
and instructions that would violate the table-free XOR experiment.  `--strict`
turns hard hot-loop regressions into a nonzero exit status.
"""

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


SCHEMA = "leopard.ff8xor.assembly-census.v1"
FUNCTION_HEADER = re.compile(r"^\s*[0-9a-fA-F]+ <(.+)>:$")
INSTRUCTION = re.compile(
    r"^\s*([0-9a-fA-F]+):\s+([A-Za-z0-9_.]+)(?:\s+(.*))?$")
NM_SYMBOL = re.compile(
    r"^\s*([0-9a-fA-F]+)\s+([0-9a-fA-F]+)\s+[A-Za-z]\s+(.+)$")
MULTIPLY_NAME = re.compile(r"MultiplyWholeBuffer<(\d+)u?>")
BUTTERFLY_NAME = re.compile(
    r"ButterflyWholeBuffer<(\d+)u?,\s*(false|true)>")
AVX_MULTIPLY_NAME = re.compile(r"avx512::Multiply(256|512)<(\d+)u?>")
AVX_BUTTERFLY_NAME = re.compile(
    r"avx512::Butterfly(256|512)<(\d+)u?,\s*(false|true)>")
STACK_MEMORY = re.compile(r"\([^)]*%(?:r|e)?(?:sp|bp)[^)]*\)")
SCALED_MEMORY = re.compile(r"\([^,]+,[^,]+,[1248]\)")
SCALED_STACK_MEMORY = re.compile(
    r"\([^)]*%(?:r|e)?(?:sp|bp)[^,]*,[^,]+,[1248]\)")
VECTOR_REGISTER = re.compile(r"%(?:[xyz]mm\d+|mm\d+)")

XOR2_MNEMONICS = {
    "xor", "xorb", "xorw", "xorl", "xorq", "pxor", "vpxor",
    "vpxord", "vpxorq",
}
TERNARY_MNEMONICS = {"vpternlogd", "vpternlogq"}
FORBIDDEN_PREFIXES = (
    "pshufb", "vpshufb", "vpermb", "vpermw", "vpgather", "vgather",
    "gf2p8", "vgf2p8", "pclmul", "vpclmul", "pmull", "vpmull",
)
FORBIDDEN_EXACT = {
    "imul", "imulb", "imulw", "imull", "imulq", "mul", "mulb",
    "mulw", "mull", "mulq", "pand", "vpand", "vpandd", "vpandq",
}


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


def parse_sizes(nm_output):
    sizes = {}
    for line in nm_output.splitlines():
        match = NM_SYMBOL.match(line)
        if not match:
            continue
        name = match.group(3)
        if (MULTIPLY_NAME.search(name) or BUTTERFLY_NAME.search(name) or
                AVX_MULTIPLY_NAME.search(name) or
                AVX_BUTTERFLY_NAME.search(name)):
            sizes[name] = int(match.group(2), 16)
    return sizes


def split_operands(text):
    # Relocation annotations and explanatory comments are not operands.
    text = text.split("\t", 1)[0].split("#", 1)[0].strip()
    operands = []
    begin = 0
    depth = 0
    for index, character in enumerate(text):
        if character == '(':
            depth += 1
        elif character == ')':
            depth = max(0, depth - 1)
        elif character == ',' and depth == 0:
            operand = text[begin:index].strip()
            if operand:
                operands.append(operand)
            begin = index + 1
    operand = text[begin:].strip()
    if operand:
        operands.append(operand)
    return operands


def family_and_coefficient(name):
    match = AVX_MULTIPLY_NAME.search(name)
    if match:
        width = "vl" if match.group(1) == "256" else "zmm"
        return "multiply_avx512" + width, int(match.group(2))
    match = AVX_BUTTERFLY_NAME.search(name)
    if match:
        width = "vl" if match.group(1) == "256" else "zmm"
        operation = "ifft" if match.group(3) == "true" else "fft"
        return operation + "_avx512" + width, int(match.group(2))
    match = MULTIPLY_NAME.search(name)
    if match:
        return "multiply", int(match.group(1))
    match = BUTTERFLY_NAME.search(name)
    if match:
        return ("ifft" if match.group(2) == "true" else "fft",
                int(match.group(1)))
    return None, None


def direct_target(operand_text):
    operands = split_operands(operand_text)
    if not operands:
        return None
    match = re.match(r"^(?:0x)?([0-9a-fA-F]+)", operands[0])
    return int(match.group(1), 16) if match else None


def loop_instruction_addresses(instructions):
    """Find instructions in cyclic CFG components, not address intervals.

    Compilers commonly place cold dispatch blocks inside the address range of
    a backwards branch.  Treating that whole range as a loop falsely labels
    one-per-buffer AVX-512 dispatch calls as inner-loop calls, so build a small
    direct-branch CFG and identify its cyclic strongly connected components.
    """
    if not instructions:
        return set()
    addresses = [item[0] for item in instructions]
    address_set = set(addresses)
    starts = {addresses[0]}
    for index, (unused_address, mnemonic, operands) in enumerate(instructions):
        mnemonic = mnemonic.lower()
        target = direct_target(operands) if mnemonic.startswith("j") else None
        if target in address_set:
            starts.add(target)
        if (mnemonic.startswith("j") or mnemonic.startswith("ret")) and \
                index + 1 < len(instructions):
            starts.add(instructions[index + 1][0])

    ordered_starts = sorted(starts)
    block_for_address = {}
    blocks = []
    start_index = {address: index for index, address in enumerate(addresses)}
    for block_index, start in enumerate(ordered_starts):
        begin = start_index[start]
        end = (start_index[ordered_starts[block_index + 1]]
               if block_index + 1 < len(ordered_starts) else len(instructions))
        members = instructions[begin:end]
        blocks.append(members)
        for address, unused_mnemonic, unused_operands in members:
            block_for_address[address] = block_index

    edges = [set() for unused in blocks]
    for block_index, members in enumerate(blocks):
        unused_address, mnemonic, operands = members[-1]
        mnemonic = mnemonic.lower()
        target = direct_target(operands) if mnemonic.startswith("j") else None
        if target in block_for_address:
            edges[block_index].add(block_for_address[target])
        unconditional = mnemonic in ("jmp", "jmpq")
        terminal = unconditional or mnemonic.startswith("ret")
        if not terminal and block_index + 1 < len(blocks):
            edges[block_index].add(block_index + 1)

    index_counter = [0]
    indices = [-1] * len(blocks)
    lowlink = [0] * len(blocks)
    stack = []
    on_stack = [False] * len(blocks)
    cyclic_blocks = set()

    def visit(vertex):
        indices[vertex] = index_counter[0]
        lowlink[vertex] = index_counter[0]
        index_counter[0] += 1
        stack.append(vertex)
        on_stack[vertex] = True
        for successor in edges[vertex]:
            if indices[successor] < 0:
                visit(successor)
                lowlink[vertex] = min(lowlink[vertex], lowlink[successor])
            elif on_stack[successor]:
                lowlink[vertex] = min(lowlink[vertex], indices[successor])
        if lowlink[vertex] == indices[vertex]:
            component = []
            while True:
                member = stack.pop()
                on_stack[member] = False
                component.append(member)
                if member == vertex:
                    break
            if len(component) > 1 or vertex in edges[vertex]:
                cyclic_blocks.update(component)

    for vertex in range(len(blocks)):
        if indices[vertex] < 0:
            visit(vertex)

    return set(address for block_index in cyclic_blocks
               for address, unused_mnemonic, unused_operands in blocks[block_index])


def instruction_census(name, size, instructions):
    family, coefficient = family_and_coefficient(name)
    counts = collections.Counter()
    mnemonic_counts = collections.Counter()

    loop_addresses = loop_instruction_addresses(instructions)

    for address, mnemonic, operand_text in instructions:
        mnemonic = mnemonic.lower()
        mnemonic_counts[mnemonic] += 1
        counts["instructions"] += 1
        operands = split_operands(operand_text)
        has_memory = "(" in operand_text
        has_vector = bool(VECTOR_REGISTER.search(operand_text))

        if mnemonic in XOR2_MNEMONICS:
            if has_vector or mnemonic.startswith(("p", "v")):
                counts["vector_xor2"] += 1
            else:
                counts["scalar_xor"] += 1
        if mnemonic in TERNARY_MNEMONICS:
            immediate = operands[0].lower() if operands else ""
            if immediate in ("$0x96", "0x96"):
                counts["vector_xor3"] += 1
            else:
                counts["non_xor_ternary"] += 1
        if mnemonic.startswith("call"):
            counts["calls"] += 1
            if address in loop_addresses:
                counts["inner_loop_calls"] += 1
        if STACK_MEMORY.search(operand_text):
            counts["stack_memory_refs"] += 1
            if has_vector:
                counts["vector_stack_refs"] += 1
        if SCALED_MEMORY.search(operand_text):
            counts["scaled_memory_refs"] += 1
        if SCALED_STACK_MEMORY.search(operand_text):
            counts["scaled_stack_refs"] += 1

        if (mnemonic in FORBIDDEN_EXACT or
                any(mnemonic.startswith(prefix) for prefix in FORBIDDEN_PREFIXES)):
            counts["forbidden_instructions"] += 1

        if has_memory and "%rip" in operand_text and has_vector:
            counts["vector_rip_memory_refs"] += 1

        if has_memory and mnemonic.startswith(("mov", "vmov")):
            destination = operands[-1] if operands else ""
            if "(" in destination:
                counts["memory_stores"] += 1
            else:
                counts["memory_loads"] += 1

    return {
        "symbol": name,
        "family": family,
        "coefficient": coefficient,
        "code_bytes": size,
        "counts": {key: counts.get(key, 0) for key in (
            "instructions", "vector_xor2", "vector_xor3", "scalar_xor",
            "memory_loads", "memory_stores", "calls", "stack_memory_refs",
            "inner_loop_calls",
            "vector_stack_refs", "scaled_memory_refs", "scaled_stack_refs",
            "vector_rip_memory_refs", "non_xor_ternary",
            "forbidden_instructions")},
        "mnemonics": dict(sorted(mnemonic_counts.items())),
    }


def parse_disassembly(disassembly, sizes):
    functions = []
    current_name = None
    current_instructions = []

    def finish():
        if current_name is None:
            return
        family, unused_coefficient = family_and_coefficient(current_name)
        if family:
            functions.append(instruction_census(
                current_name, sizes.get(current_name, 0), current_instructions))

    for line in disassembly.splitlines():
        header = FUNCTION_HEADER.match(line)
        if header:
            finish()
            current_name = header.group(1)
            current_instructions = []
            continue
        if current_name is None:
            continue
        instruction = INSTRUCTION.match(line)
        if instruction:
            current_instructions.append((
                int(instruction.group(1), 16), instruction.group(2),
                instruction.group(3) or ""))
    finish()
    return functions


def aggregate(functions):
    result = {}
    families = sorted(set(function["family"] for function in functions))
    for family in families:
        members = [function for function in functions
                   if function["family"] == family]
        totals = collections.Counter()
        for member in members:
            totals["code_bytes"] += member["code_bytes"]
            totals.update(member["counts"])
        result[family] = {
            "function_count": len(members),
            "coefficient_count": len(set(
                member["coefficient"] for member in members)),
            "totals": dict(sorted(totals.items())),
        }
    return result


def hard_violations(functions, summary):
    violations = []
    for family in sorted(summary):
        if summary[family]["coefficient_count"] != 256:
            violations.append({
                "family": family,
                "rule": "all_256_specializations_present",
                "observed": summary[family]["coefficient_count"],
            })
    rules = (
        ("inner_loop_calls", "call_inside_payload_loop"),
        ("vector_stack_refs", "possible_vector_spill_or_reload"),
        ("scaled_stack_refs", "dynamic_stack_array_indexing"),
        ("vector_rip_memory_refs", "possible_payload_lookup_or_vector_constant"),
        ("non_xor_ternary", "vpternlog_immediate_is_not_xor3_0x96"),
        ("forbidden_instructions", "shuffle_lookup_multiply_or_vector_and"),
    )
    for function in functions:
        for field, rule in rules:
            if function["counts"][field]:
                violations.append({
                    "family": function["family"],
                    "coefficient": function["coefficient"],
                    "symbol": function["symbol"],
                    "rule": rule,
                    "observed": function["counts"][field],
                })
    return violations


def representative(functions):
    selected = []
    for family in sorted(set(function["family"] for function in functions)):
        candidates = [function for function in functions
                      if function["family"] == family and
                      function["coefficient"] == 42]
        if candidates:
            # Const-propagated clones may coexist; report the largest body.
            selected.append(max(candidates, key=lambda item: item["code_bytes"]))
    return selected


def build_report(arguments):
    nm_tool = find_tool(arguments.nm, ("nm", "llvm-nm"))
    objdump_tool = find_tool(arguments.objdump, ("objdump", "llvm-objdump"))
    artifact = arguments.artifact.resolve()
    if not artifact.is_file():
        raise RuntimeError("artifact does not exist: %s" % artifact)

    nm_output = run_tool([nm_tool, "-S", "-C", "--defined-only", str(artifact)])
    sizes = parse_sizes(nm_output)
    disassembly = run_tool([
        objdump_tool, "-drwC", "--no-show-raw-insn", str(artifact)])
    functions = parse_disassembly(disassembly, sizes)
    summary = aggregate(functions)
    violations = hard_violations(functions, summary)
    return {
        "schema": SCHEMA,
        "artifact": artifact.name,
        "tools": {
            "nm": os.path.basename(nm_tool),
            "objdump": os.path.basename(objdump_tool),
        },
        "summary": summary,
        "representative_coefficient_42": representative(functions),
        "hard_violation_count": len(violations),
        "hard_violations": violations,
        "strict_pass": not violations,
    }


def human_report(report):
    lines = [
        "FF8 XOR assembly census: %s" % report["artifact"],
        "strict hot-loop checks: %s" %
        ("PASS" if report["strict_pass"] else "FAIL"),
    ]
    for family in sorted(report["summary"]):
        item = report["summary"][family]
        totals = item["totals"]
        lines.append(
            "%s: %d coefficients, %d bytes, xor2=%d xor3=%d, "
            "vector-stack=%d calls=%d (inner=%d) forbidden=%d" % (
                family, item["coefficient_count"], totals.get("code_bytes", 0),
                totals.get("vector_xor2", 0), totals.get("vector_xor3", 0),
                totals.get("vector_stack_refs", 0), totals.get("calls", 0),
                totals.get("inner_loop_calls", 0),
                totals.get("forbidden_instructions", 0)))
    for item in report["representative_coefficient_42"]:
        counts = item["counts"]
        lines.append(
            "representative %s<42>: %d bytes, %d instructions, xor2=%d "
            "xor3=%d, loads=%d stores=%d, vector-stack=%d" % (
                item["family"], item["code_bytes"], counts["instructions"],
                counts["vector_xor2"], counts["vector_xor3"],
                counts["memory_loads"], counts["memory_stores"],
                counts["vector_stack_refs"]))
    if report["hard_violations"]:
        lines.append("hard violations (first 20 of %d):" %
                     report["hard_violation_count"])
        for violation in report["hard_violations"][:20]:
            lines.append("  %s" % json.dumps(violation, sort_keys=True))
    return "\n".join(lines) + "\n"


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path,
                        help="LeopardFF8Xor object file or static archive")
    parser.add_argument("--nm", help="nm/llvm-nm executable")
    parser.add_argument("--objdump", help="objdump/llvm-objdump executable")
    parser.add_argument("--json", action="store_true",
                        help="emit JSON instead of the human summary")
    parser.add_argument("--output", type=Path,
                        help="write the selected output to this path")
    parser.add_argument("--strict", action="store_true",
                        help="fail on spills, calls, lookups, or missing kernels")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    try:
        report = build_report(arguments)
    except RuntimeError as error:
        print("assembly census unavailable: %s" % error, file=sys.stderr)
        return 2
    output = (json.dumps(report, sort_keys=True, indent=2) + "\n"
              if arguments.json else human_report(report))
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(output, encoding="utf-8")
    else:
        sys.stdout.write(output)
    return 1 if arguments.strict and not report["strict_pass"] else 0


if __name__ == "__main__":
    sys.exit(main())
