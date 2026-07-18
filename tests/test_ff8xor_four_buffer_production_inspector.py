#!/usr/bin/env python3
"""Finite parser tests for the production four-buffer assembly inspector."""

from __future__ import print_function

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOL_PATH = ROOT / "tools" / "inspect_ff8xor_four_buffer_production.py"
SPEC = importlib.util.spec_from_file_location(
    "ff8xor_four_buffer_production_inspector", str(TOOL_PATH))
INSPECTOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(INSPECTOR)


class ProductionInspectorParserTests(unittest.TestCase):
    def test_archive_member_selection_is_exact_and_unambiguous(self):
        self.assertEqual(
            INSPECTOR.select_archive_member([
                "LeopardFF8Xor.cpp.o",
                "LeopardFF8XorAVX512Four.cpp.o",
                "LeopardFF8XorAVX512.cpp.o",
            ]),
            "LeopardFF8XorAVX512Four.cpp.o")
        with self.assertRaises(RuntimeError):
            INSPECTOR.select_archive_member(["LeopardFF8Xor.cpp.o"])
        with self.assertRaises(RuntimeError):
            INSPECTOR.select_archive_member([
                "a/LeopardFF8XorAVX512Four.cpp.o",
                "b/LeopardFF8XorAVX512Four.cpp.o",
            ])

    def test_template_symbols_are_not_truncated_at_inner_angle_bracket(self):
        symbol = (
            "void leopard::ff8xor::avx512four::"
            "ApplyWhole<17u, true, false>(void*, void*, void*, void*, "
            "unsigned long)")
        text = (
            "0000 <%s>:\n"
            "   0:\tvpxord %%zmm1,%%zmm0,%%zmm0\n"
            "0006 <ordinary>:\n"
            "   6:\tret\n" % symbol)
        bodies = INSPECTOR.parse_disassembly(text)
        self.assertIn(symbol, bodies)
        self.assertEqual(
            INSPECTOR.specialization_key(symbol), (17, True, False))
        self.assertIn("vpxord", bodies[symbol])

    def test_rsp_is_spill_but_payload_rbp_is_not(self):
        payload_rbp = "\n".join([
            "   0: push %rbp",
            "   1: vmovups 0x40(%rbp,%rax,1),%zmm3",
            "   8: vpxord %zmm2,%zmm3,%zmm3",
            "   e: vmovups %zmm3,0x40(%rbp,%rax,1)",
            "  15: ret",
        ])
        observed = INSPECTOR.inspect_body(payload_rbp)
        self.assertFalse(observed["rbp_is_frame_pointer"])
        self.assertEqual(observed["vector_stack_reference_count"], 0)
        self.assertEqual(observed["payload_rbp_vector_reference_count"], 2)
        self.assertEqual(observed["payload_vector_load_count"], 1)
        self.assertEqual(observed["payload_vector_store_count"], 1)

        real_stack = "\n".join([
            "   0: push %rbp",
            "   1: mov %rsp,%rbp",
            "   4: vmovdqa64 %zmm3,-0x40(%rbp)",
            "   b: vmovdqa64 -0x40(%rbp),%zmm3",
            "  12: vmovdqa64 %zmm4,0x20(%rsp)",
            "  19: ret",
        ])
        observed = INSPECTOR.inspect_body(real_stack)
        self.assertTrue(observed["rbp_is_frame_pointer"])
        self.assertEqual(observed["vector_stack_reference_count"], 3)

    def test_calls_tables_and_forbidden_field_instructions_are_detected(self):
        body = "\n".join([
            "   0: call 20 <helper>",
            "   5: vmovdqu64 0x10(%rip),%zmm0",
            "   c: vpshufb %zmm1,%zmm2,%zmm3",
            "  12: vpclmulqdq $0,%xmm1,%xmm2,%xmm3",
            "  19: imul %rax,%rdx",
            "  1c: vpternlogd $0xe8,%zmm1,%zmm2,%zmm3",
        ])
        observed = INSPECTOR.inspect_body(body)
        self.assertEqual(observed["call_count"], 1)
        self.assertEqual(observed["rip_relative_memory_reference_count"], 1)
        # Integer IMUL is permitted for address arithmetic; the two vector
        # field/table multiplication shapes are not.
        self.assertEqual(observed["forbidden_field_instruction_count"], 2)
        self.assertEqual(observed["non_xor3_ternary_count"], 1)

    def test_memory_source_xor_and_displacement_cannot_spoof_checks(self):
        body = "\n".join([
            "   0: vpxord 0x40(%rax),%zmm1,%zmm1",
            "   8: vpternlogd $0xe8,0x96(%rax),%zmm2,%zmm3",
        ])
        observed = INSPECTOR.inspect_body(body)
        self.assertEqual(observed["non_move_vector_memory_count"], 2)
        self.assertEqual(observed["non_xor3_ternary_count"], 1)
        self.assertFalse(INSPECTOR.ternary_immediate_is_xor3(
            "$0xe8,0x96(%rax),%zmm2,%zmm3"))
        self.assertTrue(INSPECTOR.ternary_immediate_is_xor3(
            "$0x96,%zmm1,%zmm2,%zmm3"))
        self.assertTrue(INSPECTOR.ternary_immediate_is_xor3(
            "zmm3,zmm2,zmm1,150"))

    def test_nm_sizes_preserve_demangled_spaces(self):
        name = (
            "void leopard::ff8xor::avx512four::"
            "ApplyWhole<0u, false, true>(void*, void*, void*, void*, "
            "unsigned long)")
        parsed = INSPECTOR.parse_nm_sizes("0000 01af t %s\n" % name)
        self.assertEqual(parsed[name], 0x1af)

    def test_validation_rejects_missing_specializations_and_shape_errors(self):
        kernel = {
            "family": "fft",
            "lowering": "xor3",
            "tuple_index": 0,
            "code_bytes": 0,
            "call_count": 1,
            "vector_stack_reference_count": 1,
            "forbidden_field_instruction_count": 1,
            "non_xor3_ternary_count": 1,
            "non_move_vector_memory_count": 1,
            "rip_relative_memory_reference_count": 1,
            "payload_vector_load_count": 31,
            "payload_vector_store_count": 33,
            "vector_xor_count": 2,
            "expected_vector_xor_count": 3,
            "vpternlog_count": 4,
            "expected_vpternlog_count": 5,
        }
        failures = INSPECTOR.validate_kernels(
            [kernel], ["missing ApplyWhole specialization"])
        self.assertGreaterEqual(len(failures), 10)
        self.assertTrue(any("contains calls" in item for item in failures))
        self.assertTrue(any("non-XOR3" in item for item in failures))
        self.assertTrue(any("non-move vector memory" in item
                            for item in failures))
        self.assertTrue(any("payload loads" in item for item in failures))
        self.assertTrue(any("payload stores" in item for item in failures))
        self.assertTrue(any("VPTERNLOGs" in item for item in failures))

    def test_summary_reports_excess_reloads_and_permits_elided_stores(self):
        kernels = []
        for index, (loads, stores) in enumerate(((32, 32), (35, 24))):
            kernels.append({
                "family": "fft",
                "lowering": "xor2",
                "tuple_index": index,
                "code_bytes": 100 + index,
                "instruction_count": 200 + index,
                "vector_move_count": loads + stores,
                "payload_vector_load_count": loads,
                "payload_vector_store_count": stores,
            })
        summary = INSPECTOR.summarize(kernels)["fft_xor2"]
        self.assertEqual(summary["excess_payload_reload_count"], 3)
        self.assertEqual(summary["clean_memory_traffic_specialization_count"], 1)
        self.assertEqual(summary["payload_vector_store_count_min"], 24)


if __name__ == "__main__":
    unittest.main()
