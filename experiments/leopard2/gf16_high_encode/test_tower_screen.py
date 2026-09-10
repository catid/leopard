"""Synthetic protocol/adversarial tests: no native execution or actual clocks."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import run_tower_screen as collect
import replay_tower_screen as replay


def samples(row, averages):
    group = row['record']['group']
    row['record']['samples'] = [[int(v*group),float(int(v*group)/group)] for v in averages]


def fixtures():
    rows = []
    for item in collect.schedule():
        cell, order = item['cell'],item['order']
        row = dict(item,sibling_delta=0,record=collect.expected_record(cell,order,True))
        values = [130000 if s=='N' else 120000 if s=='P' and cell in (0,1,2,4,8) else
                  140000 if s=='0' and cell in (0,1,2,4,8) else 100000 for s in order]*21
        samples(row,values); rows.append(row)
    return rows


def result(rows):
    a,b = collect.analyze(rows),replay.derive(rows)
    replay.compare_analysis(a,b)
    return a


class ProtocolTests(unittest.TestCase):
    def test_positive_and_inventory(self):
        rows = fixtures(); self.assertEqual(len(rows),714)
        a = result(rows); self.assertEqual(a['decision'],'qualify_production_candidate')
        self.assertFalse(a['production_promotion'])
        self.assertEqual(sum(k.startswith('same_') for c in a['cells'] for k in c['ratios']),32)
        self.assertEqual(sum(k.startswith('within_') for c in a['cells'] for k in c['ratios']),128)
        self.assertEqual(sum(len(c['ratios']) for c in a['cells']),201)
        for c in [a['cells'][i] for i in (0,1,2,4,8)]:
            for order in ('0110','1001'): self.assertAlmostEqual(c['ratios']['paired_'+order],1.4)
            self.assertAlmostEqual(c['ratios']['native_on'],1.3)
            self.assertAlmostEqual(c['ratios']['original_on'],1.2)
            self.assertAlmostEqual(c['ratios']['original_off'],1.2/1.4)
        self.assertEqual(len({collect.label(r) for r in rows}),714)

    def test_every_cross_control(self):
        original = fixtures()
        controls = {(r['cell'],r['comparison']) for r in original if r['comparison'].startswith('same_')}
        self.assertEqual(len(controls),32)
        for cell,comp in controls:
            rows = copy.deepcopy(original)
            for row in rows:
                if row['cell']==cell and row['comparison']==comp and row['slot'] in (0,3):
                    samples(row,[v[1]*1.1 for v in row['record']['samples']])
            a = result(rows)
            self.assertGreater(a['cells'][cell]['ratios'][comp],1.02)
            self.assertEqual(a['decision'],'inconclusive_controls')

    def test_every_within_control_independently(self):
        original = fixtures()
        controls = {(r['cell'],r['comparison'],r['slot']) for r in original if r['comparison'].startswith('same_')}
        self.assertEqual(len(controls),128)
        for cell,comp,slot in controls:
            rows = copy.deepcopy(original)
            for row in rows:
                if (row['cell'],row['comparison'],row['slot'])==(cell,comp,slot):
                    samples(row,[v[1]*(1.1 if i%4 in (0,3) else 1/1.1)
                                 for i,v in enumerate(row['record']['samples'])])
            a = result(rows)
            self.assertGreater(a['cells'][cell]['ratios'][f'within_{comp}_{slot}'],1.02)
            self.assertTrue(all(1/1.02<=v<=1.02 for c in a['cells'] for k,v in c['ratios'].items()
                                if k.startswith('same_')))
            self.assertEqual(a['decision'],'inconclusive_controls')

    def test_all_neighbor_comparisons(self):
        original = fixtures()
        for cell in (3,5,6,7):
            for comp in ('paired_0110','paired_1001','original_on'):
                for factor in (.95,1.05):
                    rows = copy.deepcopy(original)
                    for row in rows:
                        if row['cell']==cell and row['comparison']==comp:
                            changed = 'P' if comp=='original_on' else '0'
                            samples(row,[v[1]*(factor if row['order'][i%4]==changed else 1)
                                         for i,v in enumerate(row['record']['samples'])])
                    self.assertEqual(result(rows)['decision'],'reject_unchanged_neighbor')

    def test_both_targets_each_order_native_and_each_round(self):
        original = fixtures()
        for cell in (0,1,2,4,8):
            for comp in ('paired_0110','paired_1001','original_on','native_on'):
                for rnd in (None,0,1,2):
                    rows = copy.deepcopy(original)
                    for row in rows:
                        if row['cell']!=cell or row['comparison']!=comp or (rnd is not None and row['round']!=rnd):
                            continue
                        order = row['order']; changed = 'N' if comp=='native_on' else 'P' if comp=='original_on' else '0'
                        samples(row,[(104000 if rnd is None else 99000) if s==changed else 100000 for s in order]*21)
                    decision = ('reject_native_gate' if comp=='native_on' else
                                'reject_original_gate' if comp=='original_on' else 'reject_target_gate')
                    self.assertEqual(result(rows)['decision'],decision)

    def test_medians_use_all_spans_and_sample_major_pairing(self):
        rows = fixtures(); row = rows[0]
        samples(row,[400000,100000,400000,100000]*11+[100000,400000,100000,400000]*10)
        self.assertEqual(collect.within(row['record']),1)
        self.assertEqual(replay.internal(row['record']),1)
        # Lower 42 values and upper 42 values both participate in median84.
        for row in rows:
            if row['cell']==0 and row['comparison']=='native_on' and row['order']=='NNNN':
                samples(row,[120000]*42+[160000]*42)
        self.assertAlmostEqual(result(rows)['cells'][0]['ratios']['native_on'],1.4)

    def test_reversed_order_is_not_inverted_gain(self):
        row = fixtures()[1]
        self.assertAlmostEqual(collect.within(row['record'],True),1.4)
        self.assertAlmostEqual(replay.internal(row['record'],True),1.4)
        self.assertAlmostEqual(collect.within(row['record']),1/1.4)

    def test_exact_target_threshold_and_adjacent_integer(self):
        for comparison,changed,rejection in (
                ('paired_0110','0','reject_target_gate'),
                ('paired_1001','0','reject_target_gate'),
                ('original_on','P','reject_original_gate'),
                ('native_on','N','reject_native_gate')):
            for cost in (110000,109999,105000):
                rows = fixtures()
                for row in rows:
                    if row['cell'] in (0,1,2,4,8) and row['comparison']==comparison:
                        samples(row,[cost if s==changed else 100000 for s in row['order']]*21)
                self.assertEqual(result(rows)['decision'],
                    'qualify_production_candidate' if cost==110000 else rejection)

    def test_bad_record_types_counts_routes_and_normalization(self):
        row = fixtures()[0]['record']
        for key,value in [('timed',False),('timed',1),('clock_source','synthetic'),('default_enabled',True),
                          ('cell',True),('group',256),('api','leo_encode'),('traced',True),('backend',0),('field',1),('encode_calls',100),
                          ('selections',100),('per_slot_calls',[25]*4),('probes',[1]*4),
                          ('schedule','1001'),('scratch_bytes',1),('output_hash','bad'),('extra',1)]:
            bad = dict(row,**{key:value})
            for check in (collect.validate,replay.record):
                with self.subTest(key=key,check=check.__name__),self.assertRaises(ValueError):
                    check(bad,0,'0110',True)
        for values in ([],[[1,1.0]]*83,[[True,1.0]]*84,[[0,0.0]]*84,[[1.0,1.0]]*84,
                       [[2**53,float(2**53)]]*84,[[1,1]]*84,[[1,float('nan')]]*84,
                       [[1,float('inf')]]*84,[[1,1.1]]*84,[[1,1.0,1]]*84):
            for check in (collect.validate,replay.record):
                with self.assertRaises(ValueError): check(dict(row,samples=values),0,'0110',True)
        small = collect.expected_record(7,'0000',True)
        small['samples'] = [[257,257/256]]*84
        collect.validate(small,7,'0000',True); replay.record(small,7,'0000',True)

    def test_preflight_records(self):
        for cell in range(9):
            for order in ('NNNN','PPPP','0000','1111'):
                row = dict(collect.expected_record(cell,order,False),samples=[])
                collect.validate(row,cell,order,False); replay.record(row,cell,order,False)
                with self.assertRaises(ValueError): collect.validate(row,cell,order,True)

    def test_bad_schedule_isolation_inventory(self):
        original = fixtures()
        for key,value in [('cell',True),('round',1),('comparison','same_off'),('slot',1),
                          ('order','1001'),('sibling_delta',False),('sibling_delta',1),('extra',1)]:
            rows = copy.deepcopy(original); rows[0][key] = value
            for check in (collect.analyze,replay.derive):
                with self.assertRaises(ValueError): check(rows)
        for rows in (original[:-1],original+[original[0]],original[1:]+original[:1]):
            for check in (collect.analyze,replay.derive):
                with self.assertRaises(ValueError): check(rows)

    def test_fixed_plan_and_pins(self):
        collect.validate_plan(collect.PROTOCOL)
        for key in collect.PROTOCOL:
            plan = copy.deepcopy(collect.PROTOCOL); plan[key] = None
            with self.subTest(key=key),self.assertRaises(ValueError): collect.validate_plan(plan)
        pins = dict(schema='leopard-tower-pins/v1',files={n:collect.PINS.get(n,'a'*64) for n in collect.FROZEN_FILES})
        collect.validate_pins(pins)
        for name in pins['files']:
            bad = copy.deepcopy(pins); del bad['files'][name]
            with self.assertRaises(ValueError): collect.validate_pins(bad)
        for name in collect.PINS:
            bad = copy.deepcopy(pins); bad['files'][name] = 'b'*64
            with self.assertRaises(ValueError): collect.validate_pins(bad)
        self.assertEqual(collect.FROZEN_FILES,replay.FILES)
        self.assertEqual(collect.SOURCE_FILES,replay.SOURCES)

    def test_replay_rejects_altered_analysis(self):
        good = result(fixtures())
        for key,value in [('decision','other'),('controls_pass',False),('production_promotion',True)]:
            with self.assertRaises(ValueError): replay.compare_analysis(dict(good,**{key:value}),good)
        bad = copy.deepcopy(good); bad['cells'][0]['ratios']['paired_0110'] *= 1.01
        with self.assertRaises(ValueError): replay.compare_analysis(bad,good)

    def test_duplicate_json_and_resource_failures(self):
        with tempfile.TemporaryDirectory(prefix='paired-protocol-') as temp:
            path = Path(temp)/'data'
            path.write_text('{"x":1,"x":2}')
            for reader in (collect.read,replay.read):
                with self.assertRaises(ValueError): reader(path)
            good = '\tExit status: 0\nmemory.peak\n123456\nmemory.max\n268435456\nmemory.events\nlow 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\nmemory.swap.current\n0\nmemory.swap.max\n0\n'
            path.write_text(good); self.assertEqual(replay.resource(path),123456)
            for text in (good.replace('Exit status: 0','Exit status: 1'),good.replace('max 0','max 1'),
                         good.replace('123456','268435457'),good.replace('oom 0','oom 1'),
                         good.replace('memory.swap.current\n0','memory.swap.current\n1')):
                path.write_text(text)
                with self.assertRaises(ValueError): replay.resource(path)

    def test_no_native_execution_in_pure_analysis(self):
        with mock.patch.object(collect.subprocess,'run',side_effect=AssertionError('execution')), \
             mock.patch.object(collect.subprocess,'check_output',side_effect=AssertionError('execution')):
            result(fixtures())

    def test_full_raw_file_replay_and_mutations(self):
        # Synthetic files test raw decoding/terminal analysis. Input byte provenance
        # is checked separately against the real immutable qualified bundle.
        rows = fixtures(); commit = 'a'*40
        plan = collect.PROTOCOL; pins = dict(files={collect.PLAN:'b'*64})
        state = dict(schema='leopard-tower-attempt/v1',preregistration=commit,
            pins=pins,host=plan['host'],plan_sha256='b'*64,complete=True,
            passive=dict(before=123,after=123,elapsed_ns=10000000000),
            preflight=[dict(collect.expected_record(c,o,False),samples=[]) for c in range(9)
                       for o in ('NNNN','PPPP','0000','1111')],invocations=rows,analysis=collect.analyze(rows))
        with tempfile.TemporaryDirectory(prefix='paired-raw-fixture-') as temp:
            root = Path(temp); attempt = root/'attempt1'; attempt.mkdir()
            def put(name,value): (attempt/name).write_text(json.dumps(value,separators=(',',':'))+'\n')
            put('attempt.json',state)
            for c in range(9):
                for i,o in enumerate(('NNNN','PPPP','0000','1111')):
                    name = f'check-{c}-{o}'
                    put(name+'.stdout',state['preflight'][c*4+i]); (attempt/(name+'.stderr')).write_text('')
            for row in rows:
                name = collect.label(row)
                put(name+'.stdout',row['record']); (attempt/(name+'.stderr')).write_text('')
            condition = ''.join(f'MainPID=0\nId={u}\nActiveState=inactive\nUnitFileState=disabled\n' for u in
                ('slipgate-catstream.service','slipgate-obs.service','slipgate-headless-x.service'))
            condition += '3fad3dfb247b3d4f3991e0a6801fba0805aae0901bdc89018841d4707a2d888c exited|false|no\n'
            condition += '002ba0a05bf860b12f2d799b49da1976aed663cb1494b2f8236751e18dc50182 exited|false|no\n'
            for name in ('condition-before','condition-after'):
                (attempt/(name+'.stdout')).write_text(condition); (attempt/(name+'.stderr')).write_text('')
            (root/'attempt1.log').write_text('\tExit status: 0\nmemory.peak\n123456\nmemory.max\n268435456\n'
                'memory.events\nlow 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\n'
                'memory.swap.current\n0\nmemory.swap.max\n0\n')
            with mock.patch.object(replay,'inputs',return_value=(plan,pins)), \
                 mock.patch.object(replay.subprocess,'run',side_effect=AssertionError('execution')), \
                 mock.patch.object(replay.subprocess,'check_output',side_effect=AssertionError('execution')):
                actual = replay.replay(root,commit)
                self.assertEqual(actual['timed_invocations'],714)
                self.assertEqual(actual['analysis']['decision'],'qualify_production_candidate')
                for key,value in [('complete',False),('complete',1),('preregistration','c'*40),
                                  ('preflight',state['preflight'][:-1]),('invocations',rows[:-1]),
                                  ('passive',dict(before=123,after=124,elapsed_ns=10000000000)),
                                  ('failure','test')]:
                    put('attempt.json',dict(state,**{key:value}))
                    with self.assertRaises(ValueError): replay.replay(root,commit)
                put('attempt.json',state)
                for name,replacement in [('check-0-NNNN.stdout','{}\n'),
                                          (collect.label(rows[0])+'.stdout','{}\n'),
                                          (collect.label(rows[0])+'.stderr','failure\n'),
                                          ('condition-after.stdout',condition.replace('disabled','enabled'))]:
                    path = attempt/name; saved = path.read_text(); path.write_text(replacement)
                    with self.assertRaises(ValueError): replay.replay(root,commit)
                    path.write_text(saved)
                (attempt/'unexpected.stdout').write_text('{}')
                with self.assertRaises(ValueError): replay.replay(root,commit)

    def test_original_off_is_not_borrowed_gain(self):
        rows = fixtures()
        for row in rows:
            if row['comparison']=='original_off':
                samples(row,[200000 if s=='P' else 100000 for s in row['order']]*21)
        a = result(rows)
        self.assertEqual(a['decision'],'qualify_production_candidate')
        for c in a['cells']:
            self.assertAlmostEqual(c['ratios']['original_off'],2)

    def test_tower_counters_and_original_identity(self):
        for order in ('0110','PPPP','NNNN'):
            row = dict(collect.expected_record(0,order,True),samples=[[100000,100000.0]]*84)
            for target in ('tower_probes','tower_totals'):
                bad = copy.deepcopy(row)
                counter = bad[target][0] if target=='tower_probes' else bad[target]
                counter['values'][0] = 1
                for check in (collect.validate,replay.record):
                    with self.assertRaises(ValueError): check(bad,0,order,True)
            for key,value in (('traced',1),('traced',True),('codec','bad'),('backend',True),('field',False)):
                for check in (collect.validate,replay.record):
                    with self.assertRaises(ValueError): check(dict(row,**{key:value}),0,order,True)

    def test_lifetime_initialization_all_cells_and_orders(self):
        for cell in range(9):
            for order in ('0110','1001','0000','1111','PPPP','NNNN'):
                for measured in (False,True):
                    row = dict(collect.expected_record(cell,order,measured),
                        samples=[[100000,100000.0/collect.PROTOCOL['groups'][cell]]]*84 if measured else [])
                    eligible = cell in (0,1,2,4,8) and order not in ('PPPP','NNNN')
                    self.assertEqual([r['initializations'] for r in row['tower_probes']],
                        [int(eligible and '1' in order[:i+1]) for i in range(4)])
                    self.assertEqual(row['tower_totals']['initializations'],int(eligible and '1' in order))
                    collect.validate(row,cell,order,measured); replay.record(row,cell,order,measured)
                    for index in range(5):
                        for key in ('values','initializations'):
                            bad = copy.deepcopy(row)
                            snapshot = bad['tower_probes'][index] if index<4 else bad['tower_totals']
                            if key=='values': snapshot[key][index] = 1
                            else: snapshot[key] = 1-snapshot[key]
                            for check in (collect.validate,replay.record):
                                with self.assertRaises(ValueError): check(bad,cell,order,measured)
                    bad = copy.deepcopy(row); bad['tower_probes'].pop()
                    for check in (collect.validate,replay.record):
                        with self.assertRaises(ValueError): check(bad,cell,order,measured)

    def test_independent_fixed_spec_and_original_controls(self):
        collect.validate_plan(replay.SPEC)
        for key in replay.SPEC:
            bad = copy.deepcopy(replay.SPEC); bad[key] = None
            with self.assertRaises(ValueError): replay.same(bad,replay.SPEC)
        self.assertEqual(replay.SPEC['target_cells'],[0,1,2,4,8])
        self.assertEqual(replay.SPEC['affected_neighbors'],[])
        self.assertEqual(replay.SPEC['unchanged_neighbors'],[3,5,6,7])
        self.assertEqual(replay.SPEC['groups'],[1]*7+[256,1])
        self.assertEqual(sum(r['order']=='PPPP' for r in fixtures()),216)

    def test_complete_maximum_interval_record_fits_existing_bound(self):
        rows = fixtures()
        for row in rows:
            g = row['record']['group']
            row['record']['samples'] = [[2**53-1,(2**53-1)/g]]*84
            collect.validate(row['record'],row['cell'],row['order'],True)
            replay.record(row['record'],row['cell'],row['order'],True)
        state = dict(schema='leopard-tower-attempt/v1',preregistration='a'*40,
            pins=dict(schema='leopard-tower-pins/v1',files={n:'a'*64 for n in collect.FROZEN_FILES}),
            host=collect.HOST,plan_sha256='a'*64,complete=True,
            passive=dict(before=2**63,after=2**63,elapsed_ns=10**10),
            preflight=[dict(collect.expected_record(c,o,False),samples=[]) for c in range(9)
                       for o in ('NNNN','PPPP','0000','1111')],
            invocations=rows,analysis=collect.analyze(rows))
        payload = json.dumps(state,separators=(',',':'),allow_nan=False)+'\n'
        self.assertLess(len(payload.encode()),4*1024**2)
        with tempfile.TemporaryDirectory(prefix='adjacent-max-record-') as temp:
            path = Path(temp)/'record.json'; path.write_text(payload)
            collect.equal(collect.read(path),state); replay.same(replay.read(path),state)


if __name__=='__main__': unittest.main()
