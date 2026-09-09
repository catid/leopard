#!/usr/bin/env python3
"""Static-audit parser checks; no codecs or timings (leopard-79h.38.5.4.18)."""
import unittest

from avx2_isa_attribution import select_recipe, summarize


class AuditTests(unittest.TestCase):
    def test_vex_is_not_evex(self):
        result = summarize("0000 <f>:\n 0:\tc5 85 ef c0 \tvpxor %ymm0,%ymm15,%ymm0\n")
        self.assertEqual(result["counts"]["instructions"], 1)
        self.assertEqual(result["counts"]["evex"], 0)
        self.assertEqual(result["counts"]["high_ymm"], 0)

    def test_evex_low_registers_still_detected(self):
        result = summarize("0000 <f>:\n 0:\t62 f3 dd 08 25 c1 96 \tvpternlogq $0x96,%xmm1,%xmm4,%xmm0\n")
        self.assertEqual(result["counts"]["evex"], 1)
        self.assertEqual(result["counts"]["ternary_logic"], 1)
        self.assertEqual(result["counts"]["zmm"], 0)

    def test_high_register_and_function_scope(self):
        result = summarize("0000 <f>:\n 0:\t62 e1 7f 28 6f 42 02 \tvmovdqu8 0x40(%rdx),%ymm16\n"
                           "0007 <g>:\n 7:\tc3 \tret\n")
        self.assertEqual(result["counts"]["high_ymm"], 1)
        self.assertEqual(result["functions"]["g"]["counts"]["evex"], 0)
        self.assertIn("%ymm16", result["functions"]["f"]["examples"]["high_ymm"])

    def test_invalid_disassembly(self):
        for raw in ("", "0:\tc3\tret\n", "0000 <f>:\n0001 <f>:\n"):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                summarize(raw)

    def test_production_recipe_not_hooks(self):
        production = {"file": "/src/a.cpp", "output": "CMakeFiles/prod.dir/a.cpp.o"}
        hooks = {"file": "/src/a.cpp", "output": "CMakeFiles/prod_test_hooks.dir/a.cpp.o"}
        self.assertEqual(select_recipe([hooks, production], "a.cpp", "prod"), production)
        for recipes in ([], [hooks], [production, production]):
            with self.subTest(recipes=recipes), self.assertRaises(ValueError):
                select_recipe(recipes, "a.cpp", "prod")


if __name__ == "__main__":
    unittest.main()
