#!/usr/bin/env python3
"""Pure adversarial replay checks; uses saved disassembly, never a codec."""
from pathlib import Path
import sys
import unittest

from replay_tower_algebra import codegen, field_proof


class ReplayTests(unittest.TestCase):
    def test_actual_object(self):
        result = codegen(ASSEMBLY)
        self.assertEqual(result['tower_product_blocks']['shuffles'], 6)
        self.assertEqual(result['tower_convert_blocks']['shuffles'], 4)

    def test_missing_shuffle(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('vpshufb', 'vpxor', 1))

    def test_evex(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('c4 62 7d 5a', '62 62 7d 5a', 1))

    def test_high_register(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('ymm10', 'ymm30', 1))

    def test_stack_reference(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('[rdi]', '[rsp]', 1))

    def test_missing_loop(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('jne', 'je', 1))

    def test_missing_function(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY[:ASSEMBLY.index('00000000000000e0 <tower_convert_blocks>:')])

    def test_noninvertible_basis(self):
        with self.assertRaises(ValueError):
            field_proof([1] * 16)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: test_replay_tower_algebra.py CODEGEN_TEXT')
    ASSEMBLY = Path(sys.argv.pop()).read_text()
    unittest.main()
