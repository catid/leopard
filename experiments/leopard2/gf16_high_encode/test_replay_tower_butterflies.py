#!/usr/bin/env python3
"""Adversarial tests of retained codegen/record checks; no codec execution."""
from pathlib import Path
import sys
import unittest
from replay_tower_butterflies import codegen, conversion_proof, native_record, stack


class ReplayTests(unittest.TestCase):
    def test_actual(self):
        result = codegen(ASSEMBLY)
        self.assertEqual(len(result), 4)
        self.assertEqual(result['tower_convert_involution']['shuffles'], 2)
        self.assertTrue(all(x['whole_function_stack_accesses'] > 0 for x in result.values()))

    def test_shuffle(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('vpshufb', 'vpxor', 1))

    def test_stack(self):
        for instruction in ('push', 'pop', 'pushfq', 'popfq', 'call', 'lcall', 'enter', 'leave'):
            with self.subTest(instruction=instruction), self.assertRaises(ValueError):
                codegen(ASSEMBLY.replace('vpsrlw', instruction, 1))
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('ymm12,YMMWORD PTR [rsi+rax*1]', 'ymm12,YMMWORD PTR [rax+rsp*1]', 1))

    def test_ret_is_not_spill(self):
        self.assertFalse(stack((0, 'ret', '')))
        self.assertTrue(stack((0, 'mov', 'rax,[rbp-8]')))

    def test_evex(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('c5 fe 6f', '62 fe 6f', 1))

    def test_higher_register(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('ymm12', 'ymm30', 1))

    def test_missing_loop(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('ja ', 'jb ', 1))

    def test_missing_dispatch(self):
        with self.assertRaises(ValueError):
            codegen(ASSEMBLY.replace('<tower_fft_out>:', '<wrong_name>:'))

    def test_exaggerated_native_claim(self):
        import json
        record = json.loads(NATIVE)
        native_record(record)
        for key, value in [('public_codec_qualified', True), ('timed', True),
                           ('exhaustive_butterfly_cases', 65536), ('boundary_butterfly_cases', 1561)]:
            with self.subTest(key=key), self.assertRaises(ValueError):
                native_record(dict(record, **{key: value}))

    def test_bad_basis(self):
        with self.assertRaises(ValueError):
            conversion_proof([1] * 16)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: test_replay_tower_butterflies.py ROOT')
    root = Path(sys.argv.pop())
    ASSEMBLY = (root / 'codegen.txt').read_text()
    NATIVE = (root / 'release.stdout').read_text()
    unittest.main()
