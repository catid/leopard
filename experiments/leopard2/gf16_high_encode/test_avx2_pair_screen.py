"""Analytical/adversarial protocol and independent-replay checks, no codecs."""
import copy
import unittest
from run_avx2_pair_screen import (ARTIFACTS,FROZEN_FILES,PROTOCOL,analyze,
                                  schedule,validate,validate_plan,validate_pins)
from replay_avx2_pair_screen import derive,compare_analysis


def fixtures():
    expected={p:[dict(profile='native' if p=='native' else 'current',mode=p,
                     cell=c,traced=False,output_hash='fixed') for c in range(8)]
              for p in ('off','on','native')}
    rows=[]
    for item in schedule():
        p,c=item['profile'],item['cell']
        duration=95 if p=='native' else 90 if p=='on' and c<5 else 100
        rows.append(dict(item,sibling_delta=0,record=dict(expected[p][c],
            encode_calls=26,samples_ns=[duration]*21,route_counts=[0,0,0,0])))
    return expected,rows


class Tests(unittest.TestCase):
    def test_analytical_ratios_and_no_promotion(self):
        expected,rows=fixtures()
        self.assertEqual(len(rows),360)
        for function in (analyze,derive):
            result=function(rows,expected)
            self.assertEqual(result['decision'],'candidate_pass')
            self.assertFalse(result['production_promotion'])
            self.assertFalse(result['authoritative_v19'])
            for cell in result['cells']:
                self.assertAlmostEqual(cell['ratios']['on_over_off'],100/90 if cell['cell']<5 else 1)
                if cell['cell']<3: self.assertAlmostEqual(cell['ratios']['on_over_native'],95/90)
        compare_analysis(analyze(rows,expected),derive(rows,expected))

    def test_every_control(self):
        controls=[(c,p) for c in range(8) for p in (('off','on','native') if c<3 else ('off','on'))]
        self.assertEqual(len(controls),19)
        for c,p in controls:
            expected,rows=fixtures()
            for row in rows:
                if row['cell']==c and row['comparison']=='same_'+p and row['slot'] in (0,3):
                    row['record']['samples_ns']=[200]*21
            for function in (analyze,derive):
                self.assertEqual(function(rows,expected)['decision'],'inconclusive_controls')

    def test_neighbor_and_target_gates(self):
        for cell,duration,decision in [(c,110,'inconclusive_controls') for c in (5,6,7)]+\
                                      [(c,110,'reject_neighbor_regression') for c in (3,4)]+\
                                      [(c,98,'below_target_gate') for c in (0,1,2)]:
            expected,rows=fixtures()
            for row in rows:
                if row['cell']==cell and row['comparison']=='on_over_off' and row['profile']=='on':
                    row['record']['samples_ns']=[duration]*21
            for function in (analyze,derive): self.assertEqual(function(rows,expected)['decision'],decision)

    def test_partial_and_reordered(self):
        expected,rows=fixtures()
        for bad in (rows[:-1],rows+[rows[0]],rows[1:]+rows[:1]):
            for function in (analyze,derive):
                with self.assertRaises(ValueError): function(bad,expected)

    def test_mutations(self):
        expected,rows=fixtures()
        for where,key,value in [('row','cell',False),('row','round',1),('row','slot',2),
            ('row','profile','on'),('row','comparison','same_off'),('row','sibling_delta',False),
            ('row','sibling_delta',1),('record','output_hash','wrong'),('record','mode','on'),
            ('record','traced',True),('record','route_counts',[0,1,0,1]),
            ('record','route_counts',[False,0,0,0]),('record','encode_calls',26.0),
            ('record','samples_ns',[True]*21),('record','samples_ns',[0]*21),
            ('record','samples_ns',[1.0]*21),('record','samples_ns',[1]*20)]:
            bad=copy.deepcopy(rows)
            target=bad[0] if where=='row' else bad[0]['record']
            target[key]=value
            for function in (analyze,derive):
                with self.subTest(where=where,key=key),self.assertRaises(ValueError): function(bad,expected)

    def test_plan_and_inventory(self):
        good=dict(copy.deepcopy(PROTOCOL),artifact_sha256={k:'a'*64 for k in ARTIFACTS})
        validate_plan(good)
        for key,value in [('cpu',4),('sibling',68),('controller_cpu',1),('rounds',3.0),
            ('attempt_budget',2),('minimum_target_gain',1.04),('control_bound',1.03),
            ('affected_neighbors',[]),('target_cells',[0]),('unchanged_neighbors',[]),
            ('orders',{}),('samples_per_process',20),('production_promotion',True)]:
            bad=copy.deepcopy(good); bad[key]=value
            with self.assertRaises(ValueError): validate_plan(bad)
        pins=dict(schema='leopard-avx2-pair-pins/v1',files={k:'a'*64 for k in FROZEN_FILES})
        validate_pins(pins)
        for key in FROZEN_FILES:
            bad=copy.deepcopy(pins); del bad['files'][key]
            with self.assertRaises(ValueError): validate_pins(bad)

    def test_exercise_calls_and_analysis(self):
        expected,rows=fixtures()
        value=dict(expected['off'][0],samples_ns=[],encode_calls=26,route_counts=[0,0,0,0])
        validate(value,expected['off'][0],False,True)
        for calls in (1,25,27,True):
            with self.assertRaises(ValueError): validate(dict(value,encode_calls=calls),expected['off'][0],False,True)
        result=derive(rows,expected)
        bad=copy.deepcopy(result); bad['cells'][0]['ratios']['on_over_off']+=0.01
        with self.assertRaises(ValueError): compare_analysis(bad,result)


if __name__=='__main__': unittest.main()
