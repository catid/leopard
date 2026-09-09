#!/usr/bin/env python3
"""Pure adversarial protocol checks; no codec execution and no clocks."""
import copy
import unittest

from verify_paired_r19932 import ARCHIVES, bad_arguments, cases, equal, expected, parse, valid, witness


class Protocol(unittest.TestCase):
    def test_matrix(self):
        rows = cases()
        self.assertEqual(len(rows), 198)
        self.assertEqual(sum(row[-1] == 0 for row in rows), 180)
        self.assertEqual(sum(row[-1] == 86 for row in rows), 18)
        self.assertEqual(len({r[0] for r in rows}), len(rows))
        self.assertEqual(sum(len(bad_arguments(p)) for p in ARCHIVES), 39)

    def test_exact_public_calls(self):
        for profile in ARCHIVES:
            schedule = 'NNNN' if profile=='native' else '0110'
            for group,total in ((1,104),(256,25604)):
                row = expected(profile,8,schedule,group,True)
                self.assertEqual(row['encode_calls'],total)
                self.assertEqual(row['per_slot_calls'],[total//4]*4)
                self.assertEqual(row['warmup_passes'],4)
                self.assertEqual(row['exercise_passes'],21)
                self.assertEqual(witness(profile,8,schedule,group,'exercise')['calls'],total)
                self.assertEqual(witness(profile,8,schedule,group,'clock-guard')['calls'],4+16*group)
                self.assertEqual(expected(profile,8,schedule,group,False)['encode_calls'],4)

    def test_route_boundaries(self):
        for cell in range(9):
            forward = expected('release',cell,'0110',1,True)
            reverse = expected('sanitize',cell,'1001',1,True)
            self.assertEqual(forward['probes'], [0,1,1,0] if cell<2 else [int(2<=cell<=4)]*4)
            self.assertEqual(reverse['probes'], [1,0,0,1] if cell<2 else [int(2<=cell<=4)]*4)
            self.assertFalse(forward['timed'])
            self.assertFalse(forward['default_enabled'])

    def test_one_item_batch_not_group_batch(self):
        row = witness('release',1,'0110',1,'exercise')
        self.assertEqual(row['apis'],[0,104,0])
        self.assertEqual(row['states'],[52,52,0])
        with self.assertRaises(ValueError):
            valid('release',1,'0110',256)

    def test_order_not_just_histogram(self):
        a = witness('release',0,'0110',1,'exercise')
        b = witness('release',0,'1001',1,'exercise')
        self.assertEqual(a['states'], b['states'])
        self.assertNotEqual(a['order_hash'], b['order_hash'])
        # Independent byte-string reference for sample-major order.
        stream = bytes([0,1,1,0]) + bytes([0]*256+[1]*512+[0]*256)*25
        h = 0xcbf29ce484222325
        for byte in stream:
            h ^= byte
            h = h * 0x100000001b3 % (2**64)
        self.assertEqual(witness('release',8,'0110',256,'exercise')['order_hash'],f'{h:016x}')

    def test_same_path_controls(self):
        for schedule,counts in (('0000',[104,0,0]),('1111',[0,104,0])):
            self.assertEqual(witness('release',0,schedule,1,'exercise')['states'],counts)
        self.assertEqual(witness('native',0,'NNNN',1,'exercise')['apis'],[0,0,104])

    def test_reject_bad_spec(self):
        for args in [('other',0,'0110',1),('release',True,'0110',1),('release',9,'0110',1),
                     ('release',0,'0101',1),('native',0,'0110',1),('release',0,'NNNN',1),
                     ('release',8,'0110',True),('release',8,'0110',256.0),
                     ('release',8,'0110',0),('release',8,'0110',257),('release',0,'0110',256)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                valid(*args)

    def test_record_mutations(self):
        good = expected('release',0,'0110',1,True)
        mutations = dict(encode_calls=100,per_slot_calls=[25]*4,group=256,
                         warmup_passes=0,exercise_passes=25,selections=4,probes=[0]*4,
                         timed=True,default_enabled=True,schedule='1001',scratch_bytes=0,
                         input_hash='0'*16,output_hash='0'*16,api='leo_encode')
        equal(copy.deepcopy(good),good)
        for key,value in mutations.items():
            bad = dict(good,**{key:value})
            with self.subTest(key=key), self.assertRaises(ValueError):
                equal(bad,good)

    def test_types_keys_and_duplicates(self):
        for a,b in ((False,0),(True,1),(1.0,1),({'a':1,'b':2},{'a':1})):
            with self.assertRaises(ValueError): equal(a,b)
        with self.assertRaises(ValueError): parse('{"calls":103,"calls":104}')
        with self.assertRaises(ValueError): parse('{"outer":{"timed":true,"timed":false}}')


if __name__=='__main__':
    unittest.main()
