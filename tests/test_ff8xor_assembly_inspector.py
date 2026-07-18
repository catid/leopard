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

    # Index-only SIB addressing has an empty base field and must still count
    # as a dynamically indexed narrow lookup.
    index_only = [
        (0x10, "movzbl", "0x0(,%rax,1),%edx"),
        (0x17, "add", "$0x1,%rax"),
        (0x1b, "jne", "10"),
        (0x1d, "ret", ""),
    ]
    index_census = inspector.instruction_census(
        "void MultiplyWholeBuffer<6u>", 30, index_only)
    require(index_census["counts"]["loop_narrow_indexed_loads"] == 1,
            "index-only narrow lookup was not detected")


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


def test_runtime_gate_table_interpreter_is_rejected():
    # The gate table bases arrive through runtime registers rather than a
    # RIP-relative static address.  Decode packed src/dst fields and use them
    # to index a runtime wire array: this is the dynamic interpreter shape the
    # named-register generator must never regress to.  Exercise both ordinary
    # 32-bit and 64-bit descriptor loads.
    fixtures = (
        (2, "movl", "(%r10,%rax,4),%r8d", "movl", "%r8d,%r9d",
         "andl", "$0xf,%r8d", "shrl", "$0x4,%r9d"),
        (3, "movq", "(%r10,%rax,8),%r8", "movq", "%r8,%r9",
         "andq", "$0xf,%r8", "shrq", "$0x20,%r9"),
    )
    for (coefficient, load, load_operands, copy, copy_operands,
         mask, mask_operands, shift, shift_operands) in fixtures:
        instructions = [
            (0x00, "mov", "%rdx,%r10"),
            (0x03, "mov", "%rcx,%r11"),
            (0x06, "xor", "%eax,%eax"),
            (0x10, load, load_operands),
            (0x15, copy, copy_operands),
            (0x18, mask, mask_operands),
            (0x1c, shift, shift_operands),
            (0x20, "movq", "(%r11,%r8,8),%r12"),
            (0x25, "xorq", "%r12,(%r11,%r9,8)"),
            (0x2a, "add", "$0x1,%rax"),
            (0x2e, "jne", "10"),
            (0x30, "ret", ""),
        ]
        census = inspector.instruction_census(
            "void MultiplyWholeBuffer<%du>" % coefficient,
            49, instructions)
        counts = census["counts"]
        require(counts["loop_narrow_indexed_loads"] == 0,
                "fixture unexpectedly relied on the old narrow-load rule")
        require(counts["static_table_indexed_refs"] == 0,
                "fixture unexpectedly acquired static table provenance")
        require(counts["loop_loaded_value_address_refs"] == 2,
                "%s-bit gate interpreter was not detected" %
                (32 if load == "movl" else 64))

        summary = {"multiply": {"coefficient_count": 256}}
        violations = inspector.hard_violations(
            [census], summary, ("multiply",))
        matching = [item for item in violations if item["rule"] ==
                    "payload_or_gate_data_used_as_memory_address"]
        require(len(matching) == 1 and matching[0]["observed"] == 2,
                "runtime gate interpreter was not a strict violation")


def test_loop_carried_gate_indices_are_rejected():
    # The next descriptor is decoded at the loop latch and used at the next
    # iteration's header.  A single lexical taint pass misses this back-edge;
    # fixed-point CFG propagation must mark both src and dst address uses.
    instructions = [
        (0x00, "xor", "%eax,%eax"),
        (0x10, "movq", "(%r11,%r8,8),%r12"),
        (0x15, "xorq", "%r12,(%r11,%r9,8)"),
        (0x1a, "add", "$0x1,%rax"),
        (0x1e, "movl", "(%r10,%rax,4),%r8d"),
        (0x23, "movl", "%r8d,%r9d"),
        (0x26, "andl", "$0xf,%r8d"),
        (0x2a, "shrl", "$0x4,%r9d"),
        (0x2e, "jne", "10"),
        (0x30, "ret", ""),
    ]
    census = inspector.instruction_census(
        "void MultiplyWholeBuffer<9u>", 49, instructions)
    require(census["counts"]["loop_loaded_value_address_refs"] == 2,
            "loop-carried gate src/dst indices escaped fixed-point taint")
    summary = {"multiply": {"coefficient_count": 256}}
    violations = inspector.hard_violations(
        [census], summary, ("multiply",))
    matching = [item for item in violations if item["rule"] ==
                "payload_or_gate_data_used_as_memory_address"]
    require(len(matching) == 1 and matching[0]["observed"] == 2,
            "loop-carried gate interpreter was not a strict violation")


def test_adjusted_static_table_base_is_rejected():
    # Static table provenance must survive ordinary pointer adjustment before
    # the loop, including both add-immediate and a non-RIP lea copy.
    instructions = [
        (0x00, "lea", "0x0(%rip),%r10"),
        (0x07, "add", "$0x20,%r10"),
        (0x0b, "lea", "0x20(%r10),%r11"),
        (0x10, "vmovdqu", "(%r11,%rax,1),%ymm0"),
        (0x16, "add", "$0x20,%rax"),
        (0x1a, "jne", "10"),
        (0x1c, "ret", ""),
    ]
    census = inspector.instruction_census(
        "void MultiplyWholeBuffer<10u>", 29, instructions)
    require(census["counts"]["static_table_indexed_refs"] == 1,
            "adjusted RIP-derived table base lost provenance")
    summary = {"multiply": {"coefficient_count": 256}}
    violations = inspector.hard_violations(
        [census], summary, ("multiply",))
    matching = [item for item in violations if item["rule"] ==
                "static_table_index_inside_payload_loop"]
    require(len(matching) == 1 and matching[0]["observed"] == 1,
            "adjusted static table read was not a strict violation")


def test_branch_merged_static_table_base_is_rejected():
    # One branch retains the RIP-derived table base while the lexical
    # fallthrough overwrites that register with a legal payload pointer.  A
    # linear provenance scan sees only the overwrite; CFG may-analysis must
    # preserve the static possibility at the shared loop header.
    instructions = [
        (0x00, "lea", "0x0(%rip),%r10"),
        (0x07, "test", "%edi,%edi"),
        (0x09, "je", "10"),
        (0x0b, "mov", "%rsi,%r10"),
        (0x10, "vmovdqu", "(%r10,%rax,1),%ymm0"),
        (0x16, "add", "$0x20,%rax"),
        (0x1a, "jne", "10"),
        (0x1c, "ret", ""),
    ]
    census = inspector.instruction_census(
        "void MultiplyWholeBuffer<11u>", 29, instructions)
    require(census["counts"]["static_table_indexed_refs"] == 1,
            "branch-merged static table base lost CFG provenance")
    summary = {"multiply": {"coefficient_count": 256}}
    violations = inspector.hard_violations(
        [census], summary, ("multiply",))
    matching = [item for item in violations if item["rule"] ==
                "static_table_index_inside_payload_loop"]
    require(len(matching) == 1 and matching[0]["observed"] == 1,
            "branch-merged static table read was not a strict violation")


def test_read_only_static_base_uses_are_rejected():
    # cmp/test read their explicit operands but do not overwrite them.  Treating
    # the final AT&T operand as a destination loses static provenance and lets
    # the subsequent indexed table loop escape the strict check.
    fixtures = (
        (12, "cmp", "$0x0,%r10"),
        (13, "test", "%r10,%r10"),
    )
    for coefficient, mnemonic, operands in fixtures:
        instructions = [
            (0x00, "lea", "0x0(%rip),%r10"),
            (0x07, mnemonic, operands),
            (0x10, "vmovdqu", "(%r10,%rax,1),%ymm0"),
            (0x16, "add", "$0x20,%rax"),
            (0x1a, "jne", "10"),
            (0x1c, "ret", ""),
        ]
        census = inspector.instruction_census(
            "void MultiplyWholeBuffer<%du>" % coefficient,
            29, instructions)
        require(census["counts"]["static_table_indexed_refs"] == 1,
                "%s incorrectly cleared static table provenance" % mnemonic)
        summary = {"multiply": {"coefficient_count": 256}}
        violations = inspector.hard_violations(
            [census], summary, ("multiply",))
        matching = [item for item in violations if item["rule"] ==
                    "static_table_index_inside_payload_loop"]
        require(len(matching) == 1 and matching[0]["observed"] == 1,
                "%s-hidden static table read was not strict" % mnemonic)


def test_rip_load_width_distinguishes_pointers_from_scalars():
    # Full-width RIP loads can fetch a stored table pointer.  Narrow globals
    # such as the runtime kernel enum and feature flags are scalar metadata and
    # must not taint an address register merely because their storage is RIP-
    # relative.
    pointer_instructions = [
        (0x00, "movq", "0x0(%rip),%r10"),
        (0x07, "cmp", "$0x0,%r10"),
        (0x10, "vmovdqu", "(%r10,%rax,1),%ymm0"),
        (0x16, "add", "$0x20,%rax"),
        (0x1a, "jne", "10"),
        (0x1c, "ret", ""),
    ]
    pointer_states = inspector.static_address_entry_states(
        pointer_instructions)
    require("%r10" in pointer_states[0x07],
            "full-width RIP pointer load lost provenance")
    pointer_census = inspector.instruction_census(
        "void MultiplyWholeBuffer<14u>", 29, pointer_instructions)
    require(pointer_census["counts"]["static_table_indexed_refs"] == 1,
            "full-width RIP pointer table read escaped strict analysis")

    scalar_instructions = [
        (0x00, "movl", "0x0(%rip),%r10d"),
        (0x07, "cmp", "$0x1,%r10d"),
        (0x0b, "ret", ""),
    ]
    scalar_states = inspector.static_address_entry_states(
        scalar_instructions)
    require("%r10" not in scalar_states[0x07],
            "narrow RIP scalar load was mistaken for a static pointer")


def test_legal_payload_addressing_and_folded_reads():
    # Compilers may represent the byte offset as a word or vector element index,
    # so scale four/eight is legal for payload.  Loaded payload values are XORed
    # and stored, but never become addresses.  The arithmetic memory operands
    # are real folded reads and supplement explicit mov/vmov load counts.
    instructions = [
        (0x00, "xor", "%eax,%eax"),
        (0x10, "movl", "(%rsi,%rcx,4),%eax"),
        (0x15, "xorl", "(%rdi,%rcx,4),%eax"),
        (0x1a, "movl", "%eax,(%rdi,%rcx,4)"),
        (0x1f, "movq", "(%rsi,%rcx,8),%rdx"),
        (0x24, "xorq", "(%rdi,%rcx,8),%rdx"),
        (0x29, "movq", "%rdx,(%rdi,%rcx,8)"),
        (0x2e, "vmovups", "(%rsi,%rcx,8),%ymm0"),
        # %rbp is a legal general payload pointer when the function did not
        # establish a frame pointer.
        (0x34, "vxorps", "(%rbp,%rcx,8),%ymm0,%ymm1"),
        (0x3a, "vpternlogq", "$0x96,(%rsi,%rcx,8),%ymm1,%ymm0"),
        (0x42, "vmovups", "%ymm0,(%rdi,%rcx,8)"),
        (0x48, "add", "$0x1,%rcx"),
        (0x4c, "jne", "10"),
        (0x4e, "ret", ""),
    ]
    census = inspector.instruction_census(
        "void MultiplyWholeBuffer<4u>", 79, instructions)
    counts = census["counts"]
    require(counts["loop_loaded_value_address_refs"] == 0,
            "legal scaled payload addressing was rejected")
    require(counts["memory_loads"] == 3,
            "explicit vector payload load count changed")
    require(counts["folded_vector_memory_reads"] == 2,
            "folded vector payload reads were not counted")
    require(counts["folded_scalar_memory_reads"] == 2,
            "folded scalar payload reads were not counted")
    require(counts["memory_stores"] == 3,
            "explicit vector payload store count changed")
    require(counts["vector_stack_refs"] == 0 and
            counts["scaled_stack_refs"] == 0,
            "general-purpose %rbp payload addressing was called a spill")
    require(counts["vector_xor2"] == 1,
            "floating-spelling vector XOR was not counted")
    require(counts["vector_xor3"] == 1,
            "memory-source XOR3 was not recognized")

    summary = {"multiply": {"coefficient_count": 256}}
    violations = inspector.hard_violations(
        [census], summary, ("multiply",))
    require(not violations,
            "legal payload addressing produced a strict violation: %r" %
            violations)

    require(inspector.memory_address_parts("0x4000(,%rax,4)") ==
            ("", "%rax", 4), "index-only SIB operand was not parsed")
    require(inspector.memory_address_parts("callee(void*, unsigned long)")
            is None, "demangled signature was mistaken for memory")


def test_real_frame_pointer_spill_is_counted():
    instructions = [
        (0x00, "push", "%rbp"),
        (0x01, "mov", "%rsp,%rbp"),
        (0x04, "xor", "%eax,%eax"),
        (0x10, "vmovdqu", "-0x40(%rbp,%rax,1),%ymm0"),
        (0x17, "vpxor", "%ymm1,%ymm0,%ymm0"),
        (0x1b, "add", "$0x20,%rax"),
        (0x1f, "jne", "10"),
        (0x21, "leave", ""),
        (0x22, "ret", ""),
    ]
    census = inspector.instruction_census(
        "void MultiplyWholeBuffer<5u>", 35, instructions)
    require(census["counts"]["vector_stack_refs"] == 1,
            "frame-pointer vector spill was not counted")
    require(census["counts"]["loop_vector_stack_refs"] == 1,
            "in-loop frame-pointer vector spill was not counted")
    require(census["counts"]["scaled_stack_refs"] == 1,
            "indexed frame-pointer stack slot was not counted")
    warnings = inspector.spill_warnings([census])
    require(len(warnings) == 1 and warnings[0]["observed"] == 1 and
            warnings[0]["inside_loop"] == 1,
            "compiler vector spill was not reported as a warning")

    summary = {"multiply": {"coefficient_count": 256}}
    violations = inspector.hard_violations(
        [census], summary, ("multiply",))
    matching = [item for item in violations
                if item["rule"] == "possible_vector_spill_or_reload"]
    require(not matching,
            "compiler spill warning was incorrectly made structural")


def main():
    test_required_families()
    test_payload_lookup_and_forbidden_shapes()
    test_allowed_xor3()
    test_runtime_gate_table_interpreter_is_rejected()
    test_loop_carried_gate_indices_are_rejected()
    test_adjusted_static_table_base_is_rejected()
    test_branch_merged_static_table_base_is_rejected()
    test_read_only_static_base_uses_are_rejected()
    test_rip_load_width_distinguishes_pointers_from_scalars()
    test_legal_payload_addressing_and_folded_reads()
    test_real_frame_pointer_spill_is_counted()
    print("FF8 XOR assembly inspector tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
