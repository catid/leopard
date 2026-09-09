import unittest
from audit_avx2_adjacent_schedule import FUNCTIONS, selected_loop, isa_ceiling
from audit_avx2_pair_schedule import FUNCTION
from test_audit_avx2_pair_schedule import fixture


class AdjacentAuditTests(unittest.TestCase):
    def test_selects_only_requested_function(self):
        text = ''.join(fixture(key == 'accumulating').replace(FUNCTION, name)
                       for key, name in FUNCTIONS.items())
        for key, name in FUNCTIONS.items():
            self.assertEqual(len(selected_loop(text, name)['stack_references']),
                             int(key == 'accumulating'))

    def test_missing_duplicate_and_wrong_body(self):
        name = FUNCTIONS['forward']
        text = fixture().replace(FUNCTION, name)
        for bad in ('', fixture(), text + text, text.replace('vpshufb', 'vpxor')):
            with self.assertRaises(ValueError):
                selected_loop(bad, name)

    def test_rejects_each_excluded_isa_family(self):
        isa_ceiling(fixture())
        for raw, assembly in [('62 01', 'nop'), ('90', 'vpxor %ymm16,%ymm0,%ymm0'),
                              ('90', 'vmovdqu %zmm0,(%rax)'), ('90', 'vpternlogd $0,%ymm0,%ymm0,%ymm0'),
                              ('90', 'vgf2p8affineqb $0,%ymm0,%ymm0,%ymm0')]:
            with self.assertRaises(ValueError):
                isa_ceiling(f' 0:\t{raw}\t{assembly}\n')

    def test_nested_range_selects_inner_but_disjoint_loops_reject(self):
        nested = fixture() + ' a:\t90\tnop\n b:\t77 f3\tjne 0 <outer>\n'
        self.assertEqual(selected_loop(nested, FUNCTION)['end'], '0x9')
        second = ''.join(f' {i + 16:x}:\tc4 e2 75 00 c2\tvpshufb %ymm2,%ymm1,%ymm0\n'
                         for i in range(8)) + ' 19:\t77 f5\tja 10 <second>\n'
        with self.assertRaises(ValueError):
            selected_loop(fixture() + second, FUNCTION)


if __name__ == '__main__':
    unittest.main()
