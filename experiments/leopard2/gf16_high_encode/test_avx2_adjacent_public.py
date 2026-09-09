import copy
import json
import unittest

from verify_avx2_adjacent_public import (ARCHIVES, CELLS, expected, structural, pair_counts,
    schedules, valid, inventory, verify_record, witness, clocks)


class PublicFrontendTests(unittest.TestCase):
    def test_exact_cells_and_structural_oracle(self):
        from replay_avx2_pair_screen import CELLS as previous
        self.assertEqual(CELLS, previous+[previous[0]])
        for c,(f,a) in enumerate(((665,384),(1330,768),(1330,768),(1793,1792),(665,384),
                                  (0,0),(0,0),(0,0),(665,384))):
            model = structural(c)
            self.assertEqual([model['forward_pairs'],model['accumulating_pairs']],[f,a])
            blocks = 64 if c==3 else 512
            self.assertEqual([model['forward_blocks'],model['accumulating_blocks']],[f*blocks,a*blocks])

    def test_isolated_probes_and_initialization_exclusion(self):
        row = expected('trace',0,'0110',1,'--exercise','steady')
        self.assertEqual(row['pair_probes'],[pair_counts('trace',0,s) for s in '110'])
        self.assertEqual(row['pair_totals']['calls'],[[33250,33250],[19200,19200]])
        self.assertEqual(expected('trace',0,'0110',1,'--check','steady')['pair_totals'],
                         dict(calls=[[0,0],[0,0]],blocks=[[0,0],[0,0]]))
        for p in ('native','original','release'):
            row = expected(p,0,schedules(p)[0],1,'--exercise','steady')
            self.assertFalse(row['traced'])
            self.assertEqual(row['pair_totals'],dict(calls=[[0,0],[0,0]],blocks=[[0,0],[0,0]]))

    def test_inventory_unique_and_no_benchmark_clocks(self):
        rows = inventory()
        self.assertEqual(len({r['label'] for r in rows}),len(rows))
        self.assertEqual(sum(r['code']==0 for r in rows),307)
        self.assertEqual(sum(r['code']==86 for r in rows),42)
        self.assertEqual(sum(r['fault'] is not None for r in rows),20)
        self.assertTrue(all('--measure' not in r['args'] for r in rows))

    def test_argument_scope(self):
        for p in ARCHIVES:
            s = schedules(p)[0]
            for c,g in ((-1,1),(9,1),(True,1),(0,True),(0,256),(8,256),(7,2)):
                with self.assertRaises(ValueError): valid(p,c,s,g)
            for s in ('0101','111','', 'NPPP'):
                with self.assertRaises(ValueError): valid(p,0,s,1)
        with self.assertRaises(ValueError): schedules('invented')
        with self.assertRaises(ValueError): structural(True)

    def test_synthetic_endpoint_and_batch_accounting(self):
        for g in (1,256):
            c = 7 if g==256 else 8
            row = expected('release',c,'1001',g,'--clock-exercise','synthetic')
            self.assertEqual(row['encode_calls'],4+100*g)
            self.assertEqual(len(row['samples']),84)
            self.assertEqual(clocks(g)['public_calls_at_clock'][0],4+16*g)
            self.assertEqual(clocks(g)['public_calls_at_clock'][-1],4+100*g)
            w = witness('release',c,'1001',g)
            self.assertEqual(w['apis'][int(c==8)],4+100*g)
            self.assertEqual(w['states'],[2+50*g,2+50*g,0])

    def test_parser_mutations(self):
        row = next(r for r in inventory() if r['label']=='trace-0-0110-1-synthetic')
        correct = expected('trace',0,'0110',1,'--clock-exercise','synthetic')
        extra = [witness('trace',0,'0110',1),clocks(1)]
        def output(value): return '\n'.join(map(json.dumps,[value,*extra]))
        verify_record(row,output(correct),'')
        for key,value in (('traced',False),('timed',True),('default_enabled',True),
                          ('clock_source','steady'),('backend',0),('field',1),('encode_calls',105),
                          ('codec','original:'+ARCHIVES['original']),('samples',[])):
            bad = copy.deepcopy(correct); bad[key] = value
            with self.assertRaises(ValueError): verify_record(row,output(bad),'')
        for key in ('calls','blocks'):
            for target in ('pair_totals','pair_probes'):
                bad = copy.deepcopy(correct)
                item = bad[target][0] if target=='pair_probes' else bad[target]
                item[key][0][1] += 1
                with self.assertRaises(ValueError): verify_record(row,output(bad),'')
        bad = copy.deepcopy(correct); bad['encode_calls'] = True
        with self.assertRaises(ValueError): verify_record(row,output(bad),'')
        with self.assertRaises(ValueError): verify_record(row,output(correct)+'\n'+json.dumps(correct),'')
        with self.assertRaises(ValueError): verify_record(row,output(correct),'sanitizer warning')

    def test_fault_and_abort_records(self):
        for row in inventory():
            if row['code']==86:
                _,c,s,g = row['args']
                verify_record(row,json.dumps(witness(row['profile'],int(c),s,int(g),4)),
                              'unexpected driver benchmark clock\n')
            if row['fault']:
                p = row['profile']; s = schedules(p)[0]
                records = [witness(p,7,s,256,4,256),clocks(256,1)]
                error = ('group duration exceeds exact binary64 integer range' if row['fault']=='huge'
                         else 'nonpositive or reversed grouped clock interval')+'\n'
                verify_record(row,'\n'.join(map(json.dumps,records)),error)
                records[0]['calls'] -= 1
                with self.assertRaises(ValueError): verify_record(row,'\n'.join(map(json.dumps,records)),error)

    def test_external_witness_and_clock_mutations(self):
        row = next(r for r in inventory() if r['label']=='release-7-1001-256-synthetic')
        records = [expected('release',7,'1001',256,'--clock-exercise','synthetic'),
                   witness('release',7,'1001',256),clocks(256)]
        verify_record(row,'\n'.join(map(json.dumps,records)),'')
        for kind in ('clock_count','begin','end','state','api','order'):
            bad = copy.deepcopy(records)
            if kind=='clock_count': bad[2]['clock_calls'] -= 1
            elif kind=='begin': bad[2]['public_calls_at_clock'][0] += 1
            elif kind=='end': bad[2]['public_calls_at_clock'][-1] -= 1
            elif kind=='state': bad[1]['states'][0] += 1
            elif kind=='api': bad[1]['apis'][0] -= 1
            else: bad[1]['order_hash'] = '0'*16
            with self.assertRaises(ValueError): verify_record(row,'\n'.join(map(json.dumps,bad)),'')

    def test_check_mode_and_empty_cli_refusals(self):
        for p in ARCHIVES:
            row = next(r for r in inventory() if r['label']==p+'-0-check')
            output = expected(p,0,schedules(p)[0],1,'--check','steady')
            verify_record(row,json.dumps(output),'')
            output['encode_calls'] += 1
            with self.assertRaises(ValueError): verify_record(row,json.dumps(output),'')
            row = next(r for r in inventory() if r['label']==p+'-bad-0')
            verify_record(row,'','usage\n')
            with self.assertRaises(ValueError): verify_record(row,'','')


if __name__=='__main__': unittest.main()
