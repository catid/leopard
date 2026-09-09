import unittest
from audit_avx2_pair_schedule import FUNCTION, instructions, pair_loop


def fixture(stack=False):
    result = '0000000000000000 <'+FUNCTION+'>:\n'
    for i in range(8):
        result += f' {i:x}:\tc4 e2 75 00 c2\tvpshufb %ymm2,%ymm1,%ymm0\n'
    if stack: result += ' 8:\tc5 fd 6f 44 24 e0\tvmovdqa -0x20(%rsp),%ymm0\n'
    result += ' 9:\t77 f5\tja 0 <'+FUNCTION+'>\n'
    return result


class AuditTests(unittest.TestCase):
    def test_loop_and_stack(self):
        self.assertEqual(pair_loop(fixture())['stack_references'], [])
        self.assertEqual(len(pair_loop(fixture(True))['stack_references']), 1)
        self.assertEqual(pair_loop(fixture())['instructions'], 9)

    def test_missing_or_duplicate_function(self):
        for text in ('',fixture().replace(FUNCTION,'other'),fixture()+fixture()):
            with self.assertRaises(ValueError): pair_loop(text)

    def test_loop_structure_rejections(self):
        for old,new in [('vpshufb','vpxor'),('ja 0 ','jmp 0 '),('ja 0 ','ja a ')]:
            with self.assertRaises(ValueError): pair_loop(fixture().replace(old,new))

    def test_raw_instruction_validation(self):
        for text in ('', ' 0:\t\tnop\n', ' 0:\txy\tnop\n', ' 0:\t90\t\n'):
            with self.assertRaises(ValueError): instructions(text)


if __name__ == '__main__':
    unittest.main()
