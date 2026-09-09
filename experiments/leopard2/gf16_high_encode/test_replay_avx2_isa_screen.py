#!/usr/bin/env python3
"""Analytical and adversarial raw-replay tests; no native execution."""
import copy
import unittest
from replay_avx2_isa_screen import compare_analysis,derive
from test_avx2_isa_screen import fixtures
from run_avx2_isa_screen import analyze


class ReplayTests(unittest.TestCase):
    def test_known_ratios(self):
        expected,rows=fixtures()
        result=derive(rows,expected)
        compare_analysis(analyze(rows,expected),result)
        for cell in result['cells']:
            self.assertAlmostEqual(cell['ratios']['pure_over_native'],.8)
            self.assertAlmostEqual(cell['ratios']['current_over_native'],100/105)
            self.assertAlmostEqual(cell['ratios']['current_over_pure'],125/105)

    def test_reject_bad_rows(self):
        expected,rows=fixtures()
        for field,value in (('cell',1),('round',1),('slot',1),('profile','pure'),
                            ('comparison','same_native'),('sibling_delta',1),('sibling_delta',False)):
            bad=copy.deepcopy(rows)
            bad[0][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError): derive(bad,expected)
        for field,value in (('profile','other'),('output_hash','wrong'),('encode_calls',1),
                            ('encode_calls',True),('samples_ns',[True]*21),('samples_ns',[0]*21),
                            ('samples_ns',[1.0]*21),('samples_ns',[1]*20)):
            bad=copy.deepcopy(rows)
            bad[0]['record'][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError): derive(bad,expected)
        for bad in (rows[:-1],rows+[rows[0]],rows[1:]+rows[:1]):
            with self.assertRaises(ValueError): derive(bad,expected)

    def test_wrong_analysis(self):
        expected,rows=fixtures()
        correct=derive(rows,expected)
        bad=copy.deepcopy(correct)
        bad['cells'][0]['ratios']['current_over_native']+=.01
        with self.assertRaises(ValueError): compare_analysis(bad,correct)
        bad=copy.deepcopy(correct)
        bad['production_promotion']=True
        with self.assertRaises(ValueError): compare_analysis(bad,correct)

    def test_control_failure(self):
        expected,rows=fixtures()
        for row in rows:
            if row['cell']==3 and row['comparison']=='same_current' and row['slot'] in (0,3):
                row['record']['samples_ns']=[200]*21
        result=derive(rows,expected)
        self.assertFalse(result['controls_pass'])
        self.assertEqual(result['decision'],'inconclusive_controls')


if __name__=='__main__': unittest.main()
