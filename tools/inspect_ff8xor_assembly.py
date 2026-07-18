#!/usr/bin/env python3
"""Create a machine-readable census of generated FF8 XOR kernels.

The inspector accepts either LeopardFF8Xor.cpp.o or the static Leopard archive.
It recognizes every coefficient-specialized whole-buffer kernel and reports
instruction shape, code size, calls, stack/vector spills, scaled stack indexing,
and instructions that would violate the table-free XOR experiment.  `--strict`
turns structural hot-loop regressions into a nonzero exit status;
`--fail-on-spills` optionally enforces compiler code-generation quality too.
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


SCHEMA = "leopard.ff8xor.assembly-census.v2"
BASE_FAMILIES = ("multiply", "fft", "ifft")
AVX2_FAMILIES = ("multiply_avx2", "fft_avx2", "ifft_avx2")
AVX512_FAMILIES = (
    "multiply_avx512vl", "fft_avx512vl", "ifft_avx512vl",
    "multiply_avx512zmm", "fft_avx512zmm", "ifft_avx512zmm",
)
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
AVX2_MULTIPLY_NAME = re.compile(r"avx2::Multiply<(\d+)u?>")
AVX2_BUTTERFLY_NAME = re.compile(
    r"avx2::Butterfly<(\d+)u?,\s*(false|true)>")
VECTOR_REGISTER = re.compile(r"%(?:[xyz]mm\d+|mm\d+)")

XOR2_MNEMONICS = {
    "xor", "xorb", "xorw", "xorl", "xorq", "pxor", "vpxor",
    "vpxord", "vpxorq", "xorps", "xorpd", "vxorps", "vxorpd",
}
TERNARY_MNEMONICS = {"vpternlogd", "vpternlogq"}
FORBIDDEN_PREFIXES = (
    "pshuf", "vpshuf", "shuf", "vshuf", "perm", "vperm",
    "palignr", "vpalignr", "pblend", "vpblend", "blend", "vblend",
    "pand", "vpand", "andps", "andpd", "andnps", "andnpd", "vand",
    "kand", "por", "vpor", "orps", "orpd", "vorps", "vorpd",
    "psll", "vpsll", "psrl", "vpsrl", "psra", "vpsra", "vprol",
    "vpror", "vpgather", "vgather",
    "gf2p8", "vgf2p8", "pclmul", "vpclmul", "pmull", "vpmull",
    "mulx",
)
FORBIDDEN_EXACT = {
    "imul", "imulb", "imulw", "imull", "imulq", "mul", "mulb",
    "mulw", "mull", "mulq",
}
NARROW_LOAD_PREFIXES = ("movzb", "movzw", "movsb", "movsw")


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
                AVX2_MULTIPLY_NAME.search(name) or
                AVX2_BUTTERFLY_NAME.search(name) or
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


def memory_address_parts(operand):
    """Return AT&T base/index/scale parts for a real memory operand.

    Demangled direct-branch targets contain parentheses too, for example a
    function signature ending in ``(void*, unsigned long)``.  Requiring a
    register base or index prevents those annotations from being mistaken for
    payload memory references while still accepting index-only SIB operands.
    """
    for match in re.finditer(r"\(([^()]*)\)", operand):
        parts = [part.strip() for part in match.group(1).split(",")]
        if not parts or len(parts) > 3:
            continue
        base = parts[0]
        index = parts[1] if len(parts) >= 2 else ""
        if base and not re.match(r"^%[A-Za-z0-9]+$", base):
            continue
        if index and not re.match(r"^%[A-Za-z0-9]+$", index):
            continue
        if not base and not index:
            continue
        scale = 1
        if len(parts) >= 3:
            if not index:
                continue
            try:
                scale = int(parts[2], 0)
            except ValueError:
                continue
            if scale not in (1, 2, 4, 8):
                continue
        return base, index, scale
    return None


def canonical_gpr(register):
    """Return the 64-bit family name for an x86 GPR spelling."""
    register = register.strip().lower()
    legacy = {
        "%rax": "%rax", "%eax": "%rax", "%ax": "%rax",
        "%al": "%rax", "%ah": "%rax",
        "%rbx": "%rbx", "%ebx": "%rbx", "%bx": "%rbx",
        "%bl": "%rbx", "%bh": "%rbx",
        "%rcx": "%rcx", "%ecx": "%rcx", "%cx": "%rcx",
        "%cl": "%rcx", "%ch": "%rcx",
        "%rdx": "%rdx", "%edx": "%rdx", "%dx": "%rdx",
        "%dl": "%rdx", "%dh": "%rdx",
        "%rsi": "%rsi", "%esi": "%rsi", "%si": "%rsi",
        "%sil": "%rsi",
        "%rdi": "%rdi", "%edi": "%rdi", "%di": "%rdi",
        "%dil": "%rdi",
        "%rbp": "%rbp", "%ebp": "%rbp", "%bp": "%rbp",
        "%bpl": "%rbp",
        "%rsp": "%rsp", "%esp": "%rsp", "%sp": "%rsp",
        "%spl": "%rsp",
    }
    result = legacy.get(register)
    if result:
        return result
    match = re.match(r"^%r(8|9|1[0-5])(?:d|w|b)?$", register)
    return ("%s%s" % ("%r", match.group(1))) if match else None


def memory_address_registers(operand):
    parts = memory_address_parts(operand)
    if parts is None:
        return None, None
    return canonical_gpr(parts[0]), canonical_gpr(parts[1])


def is_stack_memory(operand, frame_pointer_active=False):
    base, unused_index = memory_address_registers(operand)
    return base == "%rsp" or (frame_pointer_active and base == "%rbp")


def family_and_coefficient(name):
    match = AVX2_MULTIPLY_NAME.search(name)
    if match:
        return "multiply_avx2", int(match.group(1))
    match = AVX2_BUTTERFLY_NAME.search(name)
    if match:
        operation = "ifft" if match.group(2) == "true" else "fft"
        return operation + "_avx2", int(match.group(1))
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


def cyclic_instruction_components(instructions):
    """Return instruction-address sets for cyclic CFG components.

    Compilers commonly place cold dispatch blocks inside the address range of
    a backwards branch.  Treating that whole range as a loop falsely labels
    one-per-buffer AVX-512 dispatch calls as inner-loop calls, so build a small
    direct-branch CFG and identify its cyclic strongly connected components.
    """
    if not instructions:
        return []
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
    cyclic_components = []

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
                cyclic_components.append(set(
                    address for block_index in component
                    for address, unused_mnemonic, unused_operands
                    in blocks[block_index]))

    for vertex in range(len(blocks)):
        if indices[vertex] < 0:
            visit(vertex)

    return cyclic_components


def loop_instruction_addresses(instructions):
    return set(address for component in
               cyclic_instruction_components(instructions)
               for address in component)


def frame_pointer_states(instructions):
    """Return whether %rbp is an established frame base at each address."""
    states = {}
    active = False
    for address, mnemonic, operand_text in instructions:
        states[address] = active
        mnemonic = mnemonic.lower()
        operands = split_operands(operand_text)
        destination = operands[-1] if operands else ""
        destination_gpr = canonical_gpr(destination)
        source_gpr = canonical_gpr(operands[0]) if operands else None
        if destination_gpr == "%rbp":
            frame_source = (memory_address_registers(operands[0])
                            if operands else (None, None))
            if ((mnemonic.startswith("mov") and source_gpr == "%rsp") or
                    (mnemonic.startswith("lea") and
                     frame_source[0] == "%rsp")):
                active = True
            elif mnemonic.startswith(("pop", "leave")):
                active = False
            elif not mnemonic.startswith("push"):
                active = False
        elif mnemonic.startswith("leave"):
            active = False
    return states


def instruction_successors(instructions):
    """Build the direct instruction-level CFG used by loop taint analysis."""
    successors = {address: set() for address, unused_mnemonic,
                  unused_operands in instructions}
    addresses = [item[0] for item in instructions]
    address_set = set(addresses)
    for index, (address, mnemonic, operand_text) in enumerate(instructions):
        mnemonic = mnemonic.lower()
        target = direct_target(operand_text) if mnemonic.startswith("j") else None
        if target in address_set:
            successors[address].add(target)
        terminal = mnemonic in ("jmp", "jmpq") or mnemonic.startswith("ret")
        if not terminal and index + 1 < len(instructions):
            successors[address].add(addresses[index + 1])
    return successors


def static_address_transfer(mnemonic, operands, tainted):
    """Transfer RIP-derived static-address provenance through one instruction.

    This is a may-analysis: read-only or unfamiliar instructions retain the
    incoming set.  Provenance is removed only for an explicit overwrite whose
    result is known not to depend on a static address.  That conservative
    default prevents comparisons, tests, pushes, and other flag/control
    instructions from hiding a later static-table read.
    """
    result = set(tainted)
    mnemonic = mnemonic.lower()
    destination = operands[-1] if operands else ""
    destination_gpr = canonical_gpr(destination)
    source_operands = operands[:-1]
    source_gprs = set(filter(None, (
        canonical_gpr(operand) for operand in source_operands)))

    def source_depends_on_static():
        if source_gprs.intersection(tainted):
            return True
        for operand in source_operands:
            base, index = memory_address_registers(operand)
            if base in tainted or index in tainted:
                return True
        return False

    if not destination_gpr:
        return result

    if mnemonic.startswith("lea") and operands:
        if "%rip" in operands[0]:
            result.add(destination_gpr)
        else:
            base, index = memory_address_registers(operands[0])
            if base in tainted or index in tainted:
                result.add(destination_gpr)
            else:
                result.discard(destination_gpr)
    elif mnemonic.startswith("mov"):
        # A full-width RIP load may fetch a stored static pointer.  Do not give
        # the same provenance to 8/16/32-bit scalar globals such as feature
        # flags or enum selectors merely because their storage is static.
        rip_pointer_load = (
            destination.lower() == destination_gpr and
            any("%rip" in operand and
                memory_address_parts(operand) is not None
                for operand in source_operands))
        if source_depends_on_static() or rip_pointer_load:
            result.add(destination_gpr)
        else:
            result.discard(destination_gpr)
    elif mnemonic.startswith("cmov"):
        # A conditional move may retain the old destination or take its source.
        if destination_gpr in tainted or source_depends_on_static():
            result.add(destination_gpr)
        else:
            result.discard(destination_gpr)
    elif mnemonic.startswith("xchg") and len(operands) == 2:
        left = canonical_gpr(operands[0])
        right = canonical_gpr(operands[1])
        if left and right:
            left_static = left in tainted
            right_static = right in tainted
            result.discard(left)
            result.discard(right)
            if right_static:
                result.add(left)
            if left_static:
                result.add(right)
    elif (mnemonic in ("pop", "popq", "popl", "popw") or
          mnemonic.startswith("set")):
        # These overwrite the explicit destination with a stack value or a
        # boolean.  Stack-slot provenance is outside this register analysis.
        result.discard(destination_gpr)
    elif mnemonic.startswith((
            "add", "adc", "sub", "sbb", "and", "or", "xor",
            "shl", "shr", "sal", "sar", "rol", "ror", "inc",
            "dec", "neg", "not", "imul")):
        if (mnemonic.startswith(("xor", "sub")) and
                len(source_gprs) == 1 and
                destination_gpr in source_gprs):
            result.discard(destination_gpr)
        elif (destination_gpr in tainted or
              source_depends_on_static()):
            result.add(destination_gpr)
        else:
            result.discard(destination_gpr)
    return result


def static_address_entry_states(instructions):
    """Return RIP-derived static-address provenance at each instruction.

    A may-analysis over the complete direct CFG is required because a compiler
    may merge a static-table path and a payload-pointer path before one loop.
    A lexical scan can incorrectly forget the static path when the payload
    overwrite happens to appear later in disassembly order.
    """
    if not instructions:
        return {}
    by_address = {item[0]: item for item in instructions}
    successors = instruction_successors(instructions)
    predecessors = {address: set() for address in by_address}
    for address, targets in successors.items():
        for target in targets:
            predecessors[target].add(address)

    entry = {address: set() for address in by_address}
    outgoing = {address: set() for address in by_address}
    worklist = list(item[0] for item in instructions)
    queued = set(worklist)
    while worklist:
        address = worklist.pop(0)
        queued.discard(address)
        incoming = set()
        for predecessor in predecessors[address]:
            incoming.update(outgoing[predecessor])
        unused_address, mnemonic, operand_text = by_address[address]
        transferred = static_address_transfer(
            mnemonic, split_operands(operand_text), incoming)
        if incoming == entry[address] and transferred == outgoing[address]:
            continue
        entry[address] = incoming
        outgoing[address] = transferred
        for successor in successors.get(address, ()):
            if successor not in queued:
                worklist.append(successor)
                queued.add(successor)
    return entry


def loaded_data_transfer(mnemonic, operands, tainted, frame_pointer_active):
    """Transfer function for scalar values loaded from runtime memory."""
    result = set(tainted)
    mnemonic = mnemonic.lower()
    destination = operands[-1] if operands else ""
    destination_gpr = canonical_gpr(destination)
    source_gpr = canonical_gpr(operands[0]) if operands else None
    source_memory_operands = [
        operand for operand in operands[:-1]
        if memory_address_parts(operand) is not None]
    if destination_gpr and mnemonic.startswith("mov"):
        if source_memory_operands:
            if any(not is_stack_memory(operand, frame_pointer_active)
                   for operand in source_memory_operands):
                result.add(destination_gpr)
            else:
                # Stack reloads normally hold loop bounds or payload pointers,
                # not gate descriptors.  Tracking tainted spill slots would be
                # a separate, more expensive alias analysis.
                result.discard(destination_gpr)
        elif source_gpr in tainted:
            result.add(destination_gpr)
        else:
            result.discard(destination_gpr)
    elif destination_gpr and mnemonic.startswith("lea"):
        source_base, source_index = (
            memory_address_registers(operands[0])
            if operands else (None, None))
        if source_base in tainted or source_index in tainted:
            result.add(destination_gpr)
        else:
            result.discard(destination_gpr)
    elif destination_gpr and mnemonic.startswith((
            "add", "adc", "sub", "sbb", "and", "or", "xor",
            "shl", "shr", "sal", "sar", "rol", "ror", "inc",
            "dec", "neg", "not", "bswap", "cmov")):
        if (mnemonic.startswith(("xor", "sub")) and
                source_gpr == destination_gpr):
            result.discard(destination_gpr)
        elif (destination_gpr in tainted or source_gpr in tainted or
              any(not is_stack_memory(operand, frame_pointer_active)
                  for operand in source_memory_operands)):
            result.add(destination_gpr)
        else:
            result.discard(destination_gpr)
    return result


def loaded_value_address_refs(instructions, components, frame_states):
    """Count loop memory operands addressed by loop-loaded scalar data.

    This is a forward may-analysis solved to a fixed point inside each cyclic
    CFG component.  Fixed-point propagation is required because compilers can
    load and decode the next gate descriptor at a loop latch, then consume its
    src/dst indices at the next iteration's header.
    """
    by_address = {item[0]: item for item in instructions}
    successors = instruction_successors(instructions)
    total = 0
    for component in components:
        predecessors = {address: set() for address in component}
        for address in component:
            for successor in successors.get(address, ()):
                if successor in component:
                    predecessors[successor].add(address)

        entry = {address: set() for address in component}
        outgoing = {address: set() for address in component}
        worklist = list(sorted(component))
        queued = set(worklist)
        while worklist:
            address = worklist.pop(0)
            queued.discard(address)
            incoming = set()
            for predecessor in predecessors[address]:
                incoming.update(outgoing[predecessor])
            unused_address, mnemonic, operand_text = by_address[address]
            operands = split_operands(operand_text)
            transferred = loaded_data_transfer(
                mnemonic, operands, incoming,
                frame_states.get(address, False))
            if incoming == entry[address] and transferred == outgoing[address]:
                continue
            entry[address] = incoming
            outgoing[address] = transferred
            for successor in successors.get(address, ()):
                if successor in component and successor not in queued:
                    worklist.append(successor)
                    queued.add(successor)

        for address in component:
            unused_address, unused_mnemonic, operand_text = by_address[address]
            for operand in split_operands(operand_text):
                if memory_address_parts(operand) is None:
                    continue
                base, index = memory_address_registers(operand)
                if base in entry[address] or index in entry[address]:
                    total += 1
    return total


def instruction_census(name, size, instructions):
    family, coefficient = family_and_coefficient(name)
    counts = collections.Counter()
    mnemonic_counts = collections.Counter()

    cyclic_components = cyclic_instruction_components(instructions)
    loop_addresses = set(address for component in cyclic_components
                         for address in component)
    frame_states = frame_pointer_states(instructions)
    static_address_states = static_address_entry_states(instructions)

    for address, mnemonic, operand_text in instructions:
        mnemonic = mnemonic.lower()
        mnemonic_counts[mnemonic] += 1
        counts["instructions"] += 1
        operands = split_operands(operand_text)
        memory_operands = [operand for operand in operands
                           if memory_address_parts(operand) is not None]
        source_memory_operands = [operand for operand in operands[:-1]
                                  if memory_address_parts(operand) is not None]
        frame_pointer_active = frame_states.get(address, False)
        stack_memory_operands = [
            operand for operand in memory_operands
            if is_stack_memory(operand, frame_pointer_active)]
        has_memory = bool(memory_operands)
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
        if stack_memory_operands:
            counts["stack_memory_refs"] += 1
            if has_vector:
                counts["vector_stack_refs"] += 1
                if address in loop_addresses:
                    counts["loop_vector_stack_refs"] += 1
        if any(memory_address_parts(operand)[1]
               for operand in memory_operands):
            counts["scaled_memory_refs"] += 1
        if any(memory_address_parts(operand)[1]
               for operand in stack_memory_operands):
            counts["scaled_stack_refs"] += 1

        if (mnemonic in FORBIDDEN_EXACT or
                any(mnemonic.startswith(prefix) for prefix in FORBIDDEN_PREFIXES)):
            counts["forbidden_instructions"] += 1

        if has_memory and "%rip" in operand_text and has_vector:
            counts["vector_rip_memory_refs"] += 1
        if (address in loop_addresses and has_memory and
                "%rip" in operand_text):
            counts["loop_rip_memory_refs"] += 1
        if (address in loop_addresses and
                mnemonic.startswith(NARROW_LOAD_PREFIXES) and
                any(memory_address_parts(operand)[1]
                    for operand in source_memory_operands)):
            counts["loop_narrow_indexed_loads"] += 1

        # x86 may fold a vector load into an arithmetic instruction.  Keep
        # these reads separate from explicit mov/vmov loads so representative
        # instruction-shape reports do not understate payload memory traffic.
        if (has_vector and source_memory_operands and
                not mnemonic.startswith(("mov", "vmov"))):
            counts["folded_vector_memory_reads"] += len(
                source_memory_operands)
        elif (source_memory_operands and
              not mnemonic.startswith(("lea", "mov", "vmov"))):
            counts["folded_scalar_memory_reads"] += len(
                source_memory_operands)

        # A generated payload loop should never form a RIP-relative table base
        # and then index it by a runtime value; literal named-register circuits
        # need only payload pointers.  The CFG may-analysis preserves static
        # provenance across branch merges instead of trusting lexical order.
        if address in loop_addresses:
            for operand in memory_operands:
                base, index = memory_address_registers(operand)
                static_addresses = static_address_states.get(address, set())
                if (base in static_addresses or index in static_addresses):
                    counts["static_table_indexed_refs"] += 1

        if has_memory and mnemonic.startswith(("mov", "vmov")):
            destination = operands[-1] if operands else ""
            if memory_address_parts(destination) is not None:
                counts["memory_stores"] += 1
            else:
                counts["memory_loads"] += 1

    counts["loop_loaded_value_address_refs"] = loaded_value_address_refs(
        instructions, cyclic_components, frame_states)

    return {
        "symbol": name,
        "family": family,
        "coefficient": coefficient,
        "code_bytes": size,
        "counts": {key: counts.get(key, 0) for key in (
            "instructions", "vector_xor2", "vector_xor3", "scalar_xor",
            "memory_loads", "folded_vector_memory_reads",
            "folded_scalar_memory_reads", "memory_stores",
            "calls", "stack_memory_refs",
            "inner_loop_calls",
            "vector_stack_refs", "loop_vector_stack_refs",
            "scaled_memory_refs", "scaled_stack_refs",
            "vector_rip_memory_refs", "non_xor_ternary",
            "loop_rip_memory_refs", "loop_narrow_indexed_loads",
            "loop_loaded_value_address_refs",
            "static_table_indexed_refs", "forbidden_instructions")},
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


def hard_violations(functions, summary, required_families):
    violations = []
    all_families = set(summary) | set(required_families)
    for family in sorted(all_families):
        observed = (summary[family]["coefficient_count"]
                    if family in summary else 0)
        if observed != 256:
            violations.append({
                "family": family,
                "rule": "all_256_specializations_present",
                "observed": observed,
            })
    rules = (
        ("inner_loop_calls", "call_inside_payload_loop"),
        ("scaled_stack_refs", "dynamic_stack_array_indexing"),
        ("loop_rip_memory_refs", "rip_relative_read_inside_payload_loop"),
        ("loop_narrow_indexed_loads", "narrow_indexed_payload_lookup"),
        ("loop_loaded_value_address_refs",
         "payload_or_gate_data_used_as_memory_address"),
        ("static_table_indexed_refs", "static_table_index_inside_payload_loop"),
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


def spill_warnings(functions):
    warnings = []
    for function in functions:
        observed = function["counts"]["vector_stack_refs"]
        if observed:
            warnings.append({
                "family": function["family"],
                "coefficient": function["coefficient"],
                "symbol": function["symbol"],
                "rule": "possible_vector_spill_or_reload",
                "observed": observed,
                "inside_loop":
                    function["counts"]["loop_vector_stack_refs"],
            })
    return warnings


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
    required_families = set(BASE_FAMILIES)
    required_families.update(arguments.require_family)
    if arguments.require_avx2 or any(
            family in summary for family in AVX2_FAMILIES):
        required_families.update(AVX2_FAMILIES)
    if arguments.require_avx512 or any(
            family in summary for family in AVX512_FAMILIES):
        required_families.update(AVX512_FAMILIES)
    violations = hard_violations(functions, summary, required_families)
    spills = spill_warnings(functions)
    return {
        "schema": SCHEMA,
        "artifact": artifact.name,
        "tools": {
            "nm": os.path.basename(nm_tool),
            "objdump": os.path.basename(objdump_tool),
        },
        "summary": summary,
        "required_families": sorted(required_families),
        "representative_coefficient_42": representative(functions),
        "hard_violation_count": len(violations),
        "hard_violations": violations,
        "strict_pass": not violations,
        "spill_warning_count": len(spills),
        "spill_warnings": spills,
        "spill_check_pass": not spills,
    }


def human_report(report):
    lines = [
        "FF8 XOR assembly census: %s" % report["artifact"],
        "strict hot-loop checks: %s" %
        ("PASS" if report["strict_pass"] else "FAIL"),
        "vector spill check: %s (%d specialized functions affected)" % (
            "PASS" if report["spill_check_pass"] else "WARNING",
            report["spill_warning_count"]),
    ]
    for family in sorted(report["summary"]):
        item = report["summary"][family]
        totals = item["totals"]
        lines.append(
            "%s: %d coefficients, %d bytes, xor2=%d xor3=%d, "
            "vector-stack=%d (loop=%d) scaled-stack=%d vector-rip=%d "
            "calls=%d (inner=%d) forbidden=%d" % (
                family, item["coefficient_count"], totals.get("code_bytes", 0),
                totals.get("vector_xor2", 0), totals.get("vector_xor3", 0),
                totals.get("vector_stack_refs", 0),
                totals.get("loop_vector_stack_refs", 0),
                totals.get("scaled_stack_refs", 0),
                totals.get("vector_rip_memory_refs", 0),
                totals.get("calls", 0),
                totals.get("inner_loop_calls", 0),
                totals.get("forbidden_instructions", 0)))
    for item in report["representative_coefficient_42"]:
        counts = item["counts"]
        lines.append(
            "representative %s<42>: %d bytes, %d instructions, xor2=%d "
            "xor3=%d, explicit-loads=%d folded-vector-reads=%d "
            "folded-scalar-reads=%d total-load-reads=%d stores=%d, "
            "vector-stack=%d (loop=%d)" % (
                item["family"], item["code_bytes"], counts["instructions"],
                counts["vector_xor2"], counts["vector_xor3"],
                counts["memory_loads"],
                counts["folded_vector_memory_reads"],
                counts["folded_scalar_memory_reads"],
                counts["memory_loads"] +
                counts["folded_vector_memory_reads"] +
                counts["folded_scalar_memory_reads"],
                counts["memory_stores"],
                counts["vector_stack_refs"],
                counts["loop_vector_stack_refs"]))
    if report["hard_violations"]:
        lines.append("hard violations (first 20 of %d):" %
                     report["hard_violation_count"])
        for violation in report["hard_violations"][:20]:
            lines.append("  %s" % json.dumps(violation, sort_keys=True))
    if report["spill_warnings"]:
        lines.append("spill warnings (first 20 of %d):" %
                     report["spill_warning_count"])
        for warning in report["spill_warnings"][:20]:
            lines.append("  %s" % json.dumps(warning, sort_keys=True))
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
                        help="fail on structural hot-loop violations")
    parser.add_argument(
        "--fail-on-spills", action="store_true",
        help="also fail when the selected compiler spills vector values")
    parser.add_argument(
        "--require-family", action="append", default=[],
        help="require a named family with all 256 specializations")
    parser.add_argument(
        "--require-avx2", action="store_true",
        help="require complete isolated AVX2 multiply/FFT/IFFT families")
    parser.add_argument(
        "--require-avx512", action="store_true",
        help="require complete AVX-512VL and ZMM multiply/FFT/IFFT families")
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
    failed = arguments.strict and not report["strict_pass"]
    failed = failed or (arguments.fail_on_spills and
                        not report["spill_check_pass"])
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
