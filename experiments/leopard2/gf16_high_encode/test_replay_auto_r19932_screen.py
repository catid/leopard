import copy
import unittest
from test_auto_r19932_screen import fixtures
from run_auto_r19932_screen import analyze
from replay_auto_r19932_screen import compare_analysis,derive


class ReplayTests(unittest.TestCase):
    def test_independent_ratios(self):
        expected,rows = fixtures()
        result = derive(rows,expected); compare_analysis(analyze(rows,expected),result)
        self.assertEqual(result['decision'],'qualify_default_on_artifact')
        self.assertAlmostEqual(result['cells'][0]['ratios']['off_on'],1.4)
        self.assertAlmostEqual(result['cells'][1]['ratios']['native_on'],1.3)

    def test_reject_bad_records(self):
        expected,rows = fixtures()
        for field,value in [('cell',1),('round',1),('profile','on'),('sibling_delta',False),
                            ('sibling_delta',1),('slot',1),('comparison','same_off')]:
            bad = copy.deepcopy(rows); bad[0][field]=value
            with self.assertRaises(ValueError): derive(bad,expected)
        for field,value in [('output_hash','wrong'),('encode_calls',True),('encode_calls',1),
                            ('samples_ns',[True]*21),('samples_ns',[0]*21),('samples_ns',[1]*20)]:
            bad = copy.deepcopy(rows); bad[0]['record'][field]=value
            with self.assertRaises(ValueError): derive(bad,expected)
        with self.assertRaises(ValueError): derive(rows[:-1],expected)

    def test_claim_and_ratio_mutations(self):
        expected,rows = fixtures(); good = derive(rows,expected)
        for key in ('production_promotion','authoritative_v19','confidence_intervals'):
            bad = copy.deepcopy(good); bad[key]=True
            with self.assertRaises(ValueError): compare_analysis(bad,good)
        bad = copy.deepcopy(good); bad['cells'][0]['ratios']['native_on']+=0.01
        with self.assertRaises(ValueError): compare_analysis(bad,good)

    def test_decision_gates(self):
        expected,rows = fixtures()
        for comp,cell,profile,duration,decision in [('same_off',8,'off',200,'inconclusive_controls'),
            ('off_on',8,'off',105,'reject_neighbor_gate'),('off_on',1,'off',104,'reject_target_gate'),
            ('native_on',1,'native',104,'reject_native_gate')]:
            bad = copy.deepcopy(rows)
            for r in bad:
                if r['cell']==cell and r['comparison']==comp and r['profile']==profile and r['slot'] in (0,3):
                    r['record']['samples_ns']=[duration]*21
            derived = derive(bad,expected)
            self.assertEqual(derived['decision'],decision)
            compare_analysis(analyze(bad,expected),derived)


if __name__=='__main__': unittest.main()
