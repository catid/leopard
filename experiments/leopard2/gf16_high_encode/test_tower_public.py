import copy
import json
from pathlib import Path
import tempfile
import unittest

from tower_public_overlay import adapt
from verify_tower_public import (ARCHIVES,CELLS,FIELDS,SELECTED,counts,expected,inventory,
                                input_snapshot,schedules,scope_footer,split_scope,structural,valid,verify_record,witness)
from build_tower_public import sha
from verify_avx2_adjacent_public import clocks


class TowerPublicTests(unittest.TestCase):
    def test_pinned_overlay_preserves_timing_and_public_calls(self):
        folder = Path(__file__).resolve().parent
        base = (folder/'avx2_adjacent_public.cpp').read_text()
        target = adapt(base,'driver')
        for begin,end in (('        // Sample-major order:','        const unsigned expected'),
                          ('        const auto encode = [&]() {','        const auto check_buffers')):
            a = base[base.index(begin):base.index(end,base.index(begin))]
            b = target[target.index(begin):target.index(end,target.index(begin))]
            self.assertEqual(a.replace('LeoAdjacent','LeoTower').replace('select adjacent state','select tower state'),b)
        self.assertIn('LeoTowerCounts tower_probes[4]',target)
        self.assertIn('i < 4; ++i) { if (i)',target)
        self.assertNotIn('.overflow',target)
        self.assertIn('tower_probes[slot] = LeoTowerPublicCounts()',target)
        for text,kind in ((base+'\n','driver'),(base,'witness'),(base,'unknown')):
            with self.assertRaises(ValueError): adapt(text,kind)

    def test_structural_zero_skew_edges_and_policy_exclusions(self):
        self.assertEqual(SELECTED,(0,1,2,4,8))
        self.assertEqual(structural(0),[1,1000,32768000,200,6553600,3672,917,384])
        self.assertEqual(structural(2),[2*n for n in structural(0)])
        self.assertEqual(structural(1),[2,2000,65536000,398,13041664,7344,1834,768])
        self.assertEqual(structural(4),[1,1000,32768000,199,6520832,3672,917,384])
        self.assertEqual(structural(8),structural(0))
        for c in (3,5,6,7): self.assertEqual(structural(c),[0]*8)

    def test_cache_lifetime_and_off_counts(self):
        for p in ('release','trace','sanitize'):
            for s,initializations in (('0110',[0,1,1,1]),('1001',[1,1,1,1]),
                                      ('0000',[0,0,0,0]),('1111',[1,1,1,1])):
                row = expected(p,0,s,1,'--exercise','steady')
                self.assertEqual([r['initializations'] for r in row['tower_probes']],initializations)
                for state,c in zip(s,row['tower_probes']):
                    self.assertEqual(c['values'],structural(0) if state=='1' and p!='release' else [0]*8)
                self.assertEqual(row['tower_totals'],counts(p,0,s*25,any(initializations)))
                self.assertEqual(expected(p,3,s,1,'--exercise','steady')['tower_totals'],
                                 dict(values=[0]*8,initializations=0))
        for p in ('native','original','original-sanitize'):
            row = expected(p,0,schedules(p)[0],1,'--exercise','steady')
            self.assertFalse(row['traced'])
            self.assertEqual(row['tower_probes'],[dict(values=[0]*8,initializations=0)]*4)

    def test_inventory_and_all_successful_records(self):
        rows = inventory()
        self.assertEqual(len({r['label'] for r in rows}),len(rows))
        self.assertEqual(set(r['profile'] for r in rows),set(ARCHIVES))
        self.assertFalse(any('--measure' in r['args'] for r in rows))
        for row in rows:
            if row['code'] or not row['parity']: continue
            p = row['profile']; mode,c,s,g = row['args']; c,g = int(c),int(g)
            binary = row['binary']; clock = 'steady' if binary=='plain' else binary
            records = [expected(p,c,s,g,mode,clock)]
            if binary!='plain': records.append(witness(p,c,s,g))
            if clock=='synthetic': records.append(clocks(g))
            verify_record(row,'\n'.join(map(json.dumps,records)),'')

    def test_scope_exit_and_resource_mutations(self):
        stderr = 'Running as unit: run-u123.scope; invocation ID: '+'0'*32+'\n'
        resources = ('TOWER_CHILD_EXIT=0\nmemory.peak\n1000000\nmemory.max\n268435456\n'
            'memory.events\nlow 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\n'
            'memory.swap.current\n0\nmemory.swap.max\n0\n')
        self.assertEqual(split_scope('{}\n'+resources,stderr,0),('{}\n','',1000000))
        for old,new in (('TOWER_CHILD_EXIT=0','TOWER_CHILD_EXIT=1'),('max 0','max 1'),
                        ('1000000','268435456'),('268435456','536870912'),
                        ('memory.swap.current\n0','memory.swap.current\n1')):
            with self.assertRaises(ValueError): split_scope(resources.replace(old,new),stderr,0)
        with self.assertRaises(ValueError): split_scope(resources,stderr+'bad header',1)

    def test_all_tower_counter_and_metadata_mutations(self):
        row = next(r for r in inventory() if r['label']=='trace-0-0110-1-synthetic')
        records = [expected('trace',0,'0110',1,'--clock-exercise','synthetic'),
                   witness('trace',0,'0110',1),clocks(1)]
        for target in ('tower_probes','tower_totals'):
            for i in range(len(FIELDS)+1):
                bad = copy.deepcopy(records)
                item = bad[0][target][1] if target=='tower_probes' else bad[0][target]
                if i == 8: item['initializations'] += 1
                else: item['values'][i] += 1
                with self.assertRaises(ValueError): verify_record(row,'\n'.join(map(json.dumps,bad)),'')
        for key,value in (('timed',True),('traced',False),('codec','original:'+ARCHIVES['original']),
                          ('encode_calls',105),('scratch_bytes',0),('probes',[1]*4),('samples',[])):
            bad = copy.deepcopy(records); bad[0][key] = value
            with self.assertRaises(ValueError): verify_record(row,'\n'.join(map(json.dumps,bad)),'')

    def test_schedule_and_group_refusals(self):
        for p in ARCHIVES:
            for c,g in ((True,1),(9,1),(0,256),(8,256),(7,True)):
                with self.assertRaises(ValueError): valid(p,c,schedules(p)[0],g)
            with self.assertRaises(ValueError): valid(p,0,'0101',1)

    def test_public_and_clock_witness_mutations(self):
        row = next(r for r in inventory() if r['label']=='release-7-1001-256-synthetic')
        records = [expected('release',7,'1001',256,'--clock-exercise','synthetic'),
                   witness('release',7,'1001',256),clocks(256)]
        self.assertEqual(records[0]['encode_calls'],25604)
        self.assertEqual(records[2]['clock_calls'],168)
        for key in ('calls','states','apis','order_hash','clock_calls','public_calls_at_clock'):
            bad = copy.deepcopy(records); obj = bad[1] if key in bad[1] else bad[2]
            if isinstance(obj[key],list): obj[key][0] += 1
            elif isinstance(obj[key],str): obj[key] = '0'*16
            else: obj[key] += 1
            with self.assertRaises(ValueError): verify_record(row,'\n'.join(map(json.dumps,bad)),'')

    def test_build_inputs_resolve_to_retained_not_live_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root/'build').mkdir(); (root/'tools').mkdir()
            source = root/'build/tower_public_link.cpp'; source.write_text('retained-source')
            name = '/home/catid/leopard/experiments/leopard2/gf16_high_encode/tower_public_link.cpp'
            self.assertEqual(input_snapshot(root,name,sha(source)),source)
            with self.assertRaises(ValueError): input_snapshot(root,name,'0'*64)
            builder = root/'tools/build_tower_public.py'; builder.write_text('retained-builder')
            self.assertEqual(input_snapshot(root,'/gone/build_tower_public.py',sha(builder)),builder)
            with self.assertRaises(ValueError): input_snapshot(root,'/gone/unknown.input','0'*64)

    def test_build_resource_footer_is_not_native_cap(self):
        resource = ('TOWER_CHILD_EXIT=0\nmemory.peak\n326787072\nmemory.max\n536870912\n'
            'memory.events\nlow 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\n'
            'memory.swap.current\n0\nmemory.swap.max\n0\n')
        self.assertEqual(scope_footer('build complete\n'+resource,0,512*1024**2),
                         ('build complete\n',326787072))
        with self.assertRaises(ValueError): scope_footer(resource,0,256*1024**2)


if __name__ == '__main__': unittest.main()
