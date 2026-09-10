"""Pure fixtures for the retained-only GF8 diagnostic; no native execution."""
import copy
import unittest
from unittest import mock

import audit_gf8_shift as audit


def archive(rows):
    result = b'!<arch>\n'
    for name, data in rows:
        header = f'{name:<16}{0:<12}{0:<6}{0:<6}{"100644":<8}{len(data):<10}`\n'.encode()
        result += header + data + (b'\n' if len(data) % 2 else b'')
    return result


class AuditTests(unittest.TestCase):
    def test_strict_raw_json(self):
        for raw in ('{"a":1,"a":2}', '{"outer":{"a":1,"a":1}}',
                    '{"a":NaN}', '{"a":Infinity}', '{"a":-Infinity}'):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                audit.strict_json(raw)
        self.assertFalse(audit.typed_equal({'a': True}, {'a': 1}))
        self.assertFalse(audit.typed_equal({'a': 1.0}, {'a': 1}))
        self.assertTrue(audit.typed_equal({'a': 1, 'b': 2}, {'b': 2, 'a': 1}))

    def test_instruction_tail_differences(self):
        prefix = [(0, 'nop')]
        extra = prefix + [(1, 'ret')]
        self.assertEqual(audit.differences(prefix, extra),
                         [{'index': 1, 'original': None, 'current': (1, 'ret')}])
        self.assertEqual(audit.differences(extra, prefix),
                         [{'index': 1, 'original': (1, 'ret'), 'current': None}])

    def test_archive_members_and_names(self):
        data = archive([('/', b'symbols'), ('//', b'long-name.cpp.o/\n'),
                        ('/0', b'object'), ('short.o/', b'x')])
        self.assertEqual(audit.members(data), {'long-name.cpp.o': audit.digest(b'object'),
                                              'short.o': audit.digest(b'x')})

    def test_archive_rejects_corruption(self):
        good = archive([('one.o/', b'1')])
        for bad in (b'!<thin>\n', good[:-1], good + b'x', good.replace(b'`\n', b'xx'),
                    archive([('one.o/', b'1'), ('one.o/', b'2')]),
                    archive([('../one.o/', b'1')]), archive([('/999', b'1')]),
                    archive([('//', b'one.o/\n'), ('/999', b'1')])):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                audit.members(bad)

    def test_only_annotated_addresses_normalize(self):
        self.assertEqual(audit.normalize('call   2c000 <leo2_encode>'), 'call <leo2_encode>')
        self.assertEqual(audit.normalize('je  2c034 <leo2_encode+0x34>'), 'je <leo2_encode+0x34>')
        self.assertEqual(audit.normalize('mov 0x123(%rip),%eax # 12345 <mode>'),
                         'mov <mode>(%rip),%eax')
        self.assertEqual(audit.normalize('mov -0x123(%rip),%eax # 54321 <mode>'),
                         'mov <mode>(%rip),%eax')
        for before in ('mov $0x123,%eax', 'call *%rax', 'mov 0x123(%rip),%eax',
                       'mov 0x123(%rax),%eax # 12345 <mode>'):
            self.assertEqual(audit.normalize(before), before)

    def test_real_operand_changes_survive(self):
        pairs = [('call 10 <a>', 'call 20 <b>'),
                 ('mov $0x1,%eax', 'mov $0x2,%eax'),
                 ('mov 0x1(%rip),%eax # 10 <mode>', 'mov 0x2(%rip),%eax # 20 <other>'),
                 ('je 10 <a+0x1>', 'je 20 <a+0x2>'),
                 ('mov 0x1(%rip),%eax # 10 <mode>', 'mov 0x2(%rip),%ebx # 20 <mode>')]
        for a, b in pairs:
            self.assertNotEqual(audit.normalize(a), audit.normalize(b))

    def test_disassembly_coverage(self):
        raw = '  1000:\t90 \tnop\n  1001:\tc3 \tret\n'
        row, instructions = audit.parse_disassembly(raw, 0x1000, 2)
        self.assertEqual(instructions, [(0, 'nop'), (1, 'ret')])
        self.assertEqual(row['raw_sha256'], audit.digest(b'\x90\xc3'))
        for bad, start, size in ((raw, 0x1001, 2), (raw, 0x1000, 3),
                                 (raw.replace('1001:', '1002:'), 0x1000, 2),
                                 ('nothing', 0x1000, 2)):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                audit.parse_disassembly(bad, start, size)

    def test_public_loop_layout(self):
        instructions = [(0, 'nop'), (16, 'mov %rax,%rdi'), (19, 'call <leo2_encode>'),
                        (24, 'test %eax,%eax'), (26, 'jne <main.cold+0x1>'),
                        (32, 'add $0x1,%ebx'), (35, 'jne <main+0x10>')]
        a = audit.grouped_loop(instructions, 0x48d0)
        b = audit.grouped_loop(instructions, 0x4970)
        self.assertEqual(a['normalized_instructions'], b['normalized_instructions'])
        self.assertNotEqual(a['start_mod64'], b['start_mod64'])
        for bad in (instructions[:-1], instructions + instructions):
            with self.assertRaises(ValueError):
                audit.grouped_loop(bad, 0x48d0)

    def test_raw_groups_keep_all_samples(self):
        record = {'group': 256, 'samples': [[256 * (100 + i), 100 + i] for i in range(84)]}
        result = audit.process_summary(record)
        self.assertEqual(result['retained_spans'], 84)
        self.assertEqual(result['median'], 141.5)
        self.assertEqual(result['minimum'], 100)
        self.assertEqual(result['maximum'], 183)
        self.assertEqual(result['slot_medians'], [140, 141, 142, 143])
        for field, value in (('group', 1), ('samples', record['samples'][:-1])):
            bad = copy.deepcopy(record); bad[field] = value
            with self.assertRaises(ValueError):
                audit.process_summary(bad)
        bad = copy.deepcopy(record); bad['samples'][0][1] += 0.1
        with self.assertRaises(ValueError):
            audit.process_summary(bad)

    def test_codec_execution_refused(self):
        with mock.patch.object(audit.subprocess, 'run') as run:
            with self.assertRaises(ValueError):
                audit.tool(['/tmp/current', '--check'])
            run.assert_not_called()


if __name__ == '__main__':
    unittest.main()
