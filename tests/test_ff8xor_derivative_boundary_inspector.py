#!/usr/bin/env python3

"""Parser regressions for the derivative-boundary assembly gate."""

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import inspect_ff8xor_derivative_boundary as inspector  # noqa: E402


def body(payload):
    return """
   0: mov    %%rsp,%%rbp
   3: lea    0x0(%%rip),%%rax
   a: notrack jmp *%%rax
%s
""" % payload


class DerivativeBoundaryInspectorTests(unittest.TestCase):
    def test_clean_metadata_dispatch(self):
        result = inspector.inspect_body(body("""
  10: vmovdqu (%%rsi),%%zmm0
  16: vpternlogd $0x96,(%%rdx),%%zmm1,%%zmm0
  1d: vmovdqu %%zmm0,(%%rdi)
  23: jne    10
  25: ret
"""))
        self.assertEqual([], result["failures"])
        self.assertTrue(result["indirect_jump_before_payload"])

    def test_clean_comparison_tree_dispatch(self):
        result = inspector.inspect_body("""
   0: cmp    $0x3,%edi
   3: jne    20
  10: vmovdqu (%rsi),%zmm0
  16: vpternlogd $0x96,(%rdx),%zmm1,%zmm0
  1d: vmovdqu %zmm0,(%rdi)
  23: jne    10
  25: ret
""")
        self.assertEqual([], result["failures"])
        self.assertFalse(result["indirect_jump_before_payload"])

    def test_dynamic_payload_dispatch_rejected(self):
        result = inspector.inspect_body(body("""
  10: vmovdqu (%%rsi),%%zmm0
  16: jmp    *%%rdx
  18: vpternlogd $0x96,%%zmm2,%%zmm1,%%zmm0
"""))
        self.assertTrue(any("one source-arity" in value
                            for value in result["failures"]))

    def test_call_rejected(self):
        result = inspector.inspect_body(body("""
  10: vmovdqu (%%rsi),%%ymm0
  16: vpternlogd $0x96,%%ymm2,%%ymm1,%%ymm0
  1d: call   helper
"""))
        self.assertTrue(any("calls" in value for value in result["failures"]))

    def test_vector_stack_operand_rejected(self):
        result = inspector.inspect_body(body("""
  10: vmovdqu 0x20(%%rsp),%%zmm0
  18: vpternlogd $0x96,%%zmm2,%%zmm1,%%zmm0
"""))
        self.assertTrue(any("stack" in value for value in result["failures"]))

    def test_bad_ternary_immediate_rejected(self):
        result = inspector.inspect_body(body("""
  10: vmovdqu (%%rsi),%%zmm0
  18: vpternlogd $0x69,%%zmm2,%%zmm1,%%zmm0
"""))
        self.assertTrue(any("0x96" in value for value in result["failures"]))

    def test_field_shuffle_rejected(self):
        result = inspector.inspect_body(body("""
  10: vmovdqu (%%rsi),%%ymm0
  18: vpternlogd $0x96,%%ymm2,%%ymm1,%%ymm0
  20: vpshufb %%ymm3,%%ymm0,%%ymm0
"""))
        self.assertTrue(any("field-multiplication" in value
                            for value in result["failures"]))

    def test_scalar_stack_bookkeeping_allowed(self):
        result = inspector.inspect_body(body("""
  10: mov    0x20(%%rbp),%%r12
  18: vmovdqu (%%rsi),%%ymm0
  20: vpternlogd $0x96,(%%r12),%%ymm1,%%ymm0
"""))
        self.assertEqual([], result["failures"])

    def test_folded_payload_source_allowed(self):
        result = inspector.inspect_body(body("""
  10: vmovdqu (%%rsi),%%zmm0
  18: vpternlogd $0x96,0x40(%%rdx),%%zmm1,%%zmm0
  20: vmovdqu %%zmm0,(%%rdi)
"""))
        self.assertEqual([], result["failures"])


if __name__ == "__main__":
    unittest.main()
