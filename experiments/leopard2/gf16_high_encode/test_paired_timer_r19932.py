#!/usr/bin/env python3
"""Adversarial replay expectations; no native execution or real clocks."""
import math
from pathlib import Path
import unittest
from verify_paired_r19932 import equal
from verify_paired_timer_r19932 import clock_expected, driver_expected, fault_witness, inventory, bad_args


class TimerProtocol(unittest.TestCase):
    def test_inventory(self):
        rows = inventory()
        self.assertEqual(len(rows),254)
        self.assertEqual(len({r[0] for r in rows}),254)
        self.assertEqual(sum(r[6]==0 for r in rows),182)
        self.assertEqual(sum(r[5]=='synthetic' for r in rows),90)
        self.assertEqual(sum(r[5]=='abort' for r in rows),18)
        self.assertEqual(sum(r[5]=='fault' for r in rows),12)
        self.assertEqual(sum(r[5]=='bad' for r in rows),42)

    def test_boundaries(self):
        for group in (1,256):
            row = clock_expected(group)
            calls = row['public_calls_at_clock']
            self.assertEqual(row['clock_calls'],168)
            self.assertEqual(calls[0],4+4*4*group)
            self.assertEqual(calls[-1],4+25*4*group)
            for i in range(84):
                self.assertEqual(calls[2*i+1]-calls[2*i],group)
                if i<83: self.assertEqual(calls[2*i+1],calls[2*i+2])

    def test_fractional_normalization(self):
        row = driver_expected('release',8,'0110',256,True)
        self.assertEqual(len(row['samples']),84)
        self.assertEqual(row['samples'][0],[257,1.00390625])
        for i,(elapsed,average) in enumerate(row['samples']):
            self.assertEqual(elapsed,257+17*i)
            self.assertEqual(average*256,elapsed)
            self.assertTrue(math.isfinite(average))
        self.assertFalse(row['timed'])
        self.assertEqual(row['clock_source'],'synthetic')

    def test_plain_no_samples(self):
        for cell in range(9):
            row = driver_expected('release',cell,'0110',1,False)
            self.assertEqual(row['samples'],[])
            self.assertEqual(row['clock_source'],'steady')
            self.assertEqual(row['encode_calls'],104)
            self.assertFalse(row['timed'])

    def test_fault_stops_after_first_group(self):
        for p in ('release','sanitize','native'):
            row = fault_witness(p,256)
            self.assertEqual(row['calls'],4+17*256)
            self.assertEqual(sum(row['apis']),row['calls'])
            self.assertEqual(sum(row['states']),row['calls'])
            clock = clock_expected(256,1)
            self.assertEqual(clock['public_calls_at_clock'],[4100,4356])

    def test_no_measure_invocations(self):
        for p in ('release','sanitize','native'):
            self.assertTrue(all('--measure' not in args for args in bad_args(p)))

    def test_allocation_guards_identical_to_qualified_prototype(self):
        source = Path(__file__).resolve().parent
        original = (source/'paired_r19932.cpp').read_text()
        current = (source/'paired_timer_r19932.cpp').read_text()
        def buffer(text):
            return text.split('struct Buffer\n',1)[1].split('uint64_t Hash(',1)[0]
        self.assertEqual(buffer(original),buffer(current))

    def test_plain_guard_cli_refusal_is_in_matrix(self):
        for p in ('release','sanitize','native'):
            s = 'NNNN' if p=='native' else '0110'
            self.assertIn(['--clock-guard','8',s,'256'],bad_args(p))

    def test_reject_clock_boundary_mutations(self):
        good = clock_expected(256)
        for index in (0,1,2,83,167):
            for delta in (-256,-1,1,256):
                bad = dict(good,public_calls_at_clock=list(good['public_calls_at_clock']))
                bad['public_calls_at_clock'][index] += delta
                with self.assertRaises(ValueError): equal(bad,good)
        with self.assertRaises(ValueError): equal(dict(good,clock_calls=167),good)

    def test_reject_wrong_cost_and_provenance(self):
        good = driver_expected('release',8,'0110',256,True)
        for average in (257.,1.,0.,257/255):
            samples = list(good['samples']); samples[0] = [257,average]
            with self.assertRaises(ValueError): equal(dict(good,samples=samples),good)
        for key,value in (('clock_source','steady'),('timed',True),('group',1),('default_enabled',True)):
            with self.assertRaises(ValueError): equal(dict(good,**{key:value}),good)


if __name__=='__main__': unittest.main()
