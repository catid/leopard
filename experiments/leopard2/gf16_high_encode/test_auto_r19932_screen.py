"""Pure protocol tests, no codec execution or sample collection."""
import copy
import unittest
from run_auto_r19932_screen import (ARTIFACTS,FROZEN_FILES,PROTOCOL,analyze,identity,schedule,
                                     validate,validate_pins,validate_plan)
from verify_auto_r19932_checks import public_expected


def fixtures():
    expected = {v:[identity(public_expected('native' if v=='native' else 'release',int(v=='on'),c))
                   for c in range(9)] for v in ('native','off','on')}
    rows = []
    for item in schedule():
        cell,profile = item['cell'],item['profile']
        duration = 140 if cell<2 and profile=='off' else 130 if profile=='native' else 100
        rows.append(dict(item,sibling_delta=0,record=dict(expected[profile][cell],
                         encode_calls=26,samples_ns=[duration]*21)))
    return expected,rows


class ProtocolTests(unittest.TestCase):
    def test_positive(self):
        expected,rows = fixtures(); self.assertEqual(len(rows),372)
        result = analyze(rows,expected)
        self.assertEqual(result['decision'],'qualify_default_on_artifact')
        self.assertFalse(result['production_promotion'])
        for c in result['cells'][:2]:
            self.assertAlmostEqual(c['ratios']['off_on'],1.4)
            self.assertAlmostEqual(c['ratios']['native_on'],1.3)
        self.assertEqual(sum(k.startswith('same_') for c in result['cells'] for k in c['ratios']),20)

    def test_all_control_gates(self):
        expected,rows = fixtures()
        for cell,comparison in {(r['cell'],r['comparison']) for r in rows if r['comparison'].startswith('same_')}:
            bad = copy.deepcopy(rows)
            for r in bad:
                if r['cell']==cell and r['comparison']==comparison and r['slot'] in (0,3):
                    r['record']['samples_ns']=[200]*21
            self.assertEqual(analyze(bad,expected)['decision'],'inconclusive_controls')

    def test_seven_neighbors_both_directions(self):
        expected,rows = fixtures()
        for cell in range(2,9):
            for value in (95,105):
                bad = copy.deepcopy(rows)
                for r in bad:
                    if r['cell']==cell and r['comparison']=='off_on' and r['profile']=='off':
                        r['record']['samples_ns']=[value]*21
                self.assertEqual(analyze(bad,expected)['decision'],'reject_neighbor_gate')

    def test_each_target_and_native_gate(self):
        expected,rows = fixtures()
        for cell in range(2):
            for comp,profile,decision in [('off_on','off','reject_target_gate'),
                                          ('native_on','native','reject_native_gate')]:
                bad = copy.deepcopy(rows)
                for r in bad:
                    if r['cell']==cell and r['comparison']==comp and r['profile']==profile:
                        r['record']['samples_ns']=[104]*21
                self.assertEqual(analyze(bad,expected)['decision'],decision)

    def test_each_target_round_must_improve(self):
        expected,rows = fixtures()
        for comp,profile,decision in [('off_on','off','reject_target_gate'),('native_on','native','reject_native_gate')]:
            bad = copy.deepcopy(rows)
            for r in bad:
                if r['cell']==0 and r['round']==0 and r['comparison']==comp and r['profile']==profile:
                    r['record']['samples_ns']=[99]*21
            self.assertEqual(analyze(bad,expected)['decision'],decision)

    def test_order_isolation_and_samples(self):
        expected,rows = fixtures()
        for key,value in [('cell',True),('round',1),('comparison','same_on'),('profile','on'),
                          ('slot',1),('sibling_delta',False),('sibling_delta',1)]:
            bad = copy.deepcopy(rows); bad[0][key]=value
            with self.assertRaises(ValueError): analyze(bad,expected)
        for values in ([True]*21,[0]*21,[1.0]*21,[100]*20):
            bad = copy.deepcopy(rows); bad[0]['record']['samples_ns']=values
            with self.assertRaises(ValueError): analyze(bad,expected)
        for bad in (rows[:-1],rows+[rows[0]],rows[1:]+rows[:1]):
            with self.assertRaises(ValueError): analyze(bad,expected)

    def test_workload_identity(self):
        expected,rows = fixtures(); good = rows[0]['record']; wanted = expected['off'][0]
        validate(good,wanted,True)
        validate(dict(wanted,encode_calls=1,samples_ns=[]),wanted,False)
        for key,value in [('execution_route','gfni'),('codec_commit','other'),('cell',False),
                          ('encode_calls',True),('encode_calls',1),('boundary_mode',1),
                          ('untimed_route_calls',1),('output_hash','other'),('api','leo_encode')]:
            bad = copy.deepcopy(good); bad[key]=value
            with self.assertRaises(ValueError): validate(bad,wanted,True)

    def test_fixed_plan(self):
        plan = copy.deepcopy(PROTOCOL); plan['artifact_sha256']={n:'a'*64 for n in ARTIFACTS}
        validate_plan(plan)
        for key,value in [('cpu',22),('sibling',86),('controller_cpu',1),('rounds',3.0),
                          ('attempt_budget',2),('samples_per_process',20),('passive_seconds',0),
                          ('minimum_gain',1.04),('native_minimum_gain',1),('equivalence_bound',1.03),
                          ('codec_commit','other'),('native_commit','other'),('orders',{}),
                          ('native_orders',{}),('cells',[]),('native_cells',[0]),
                          ('host',{}),('attempt_root','elsewhere'),('production_promotion',0),
                          ('core_sha256','other'),('header_sha256','other'),('artifact_sha256',{})]:
            bad = copy.deepcopy(plan); bad[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError): validate_plan(bad)

    def test_exact_pins(self):
        pins = dict(schema='leopard-auto-r19932-pins/v1',files={n:'a'*64 for n in FROZEN_FILES})
        validate_pins(pins)
        for field in ('run_auto_r19932_screen.py','current','expected.json','leopard2.cpp'):
            bad = copy.deepcopy(pins); del bad['files'][field]
            with self.assertRaises(ValueError): validate_pins(bad)
        bad = copy.deepcopy(pins); bad['files']['unexpected']='a'*64
        with self.assertRaises(ValueError): validate_pins(bad)


if __name__=='__main__': unittest.main()
