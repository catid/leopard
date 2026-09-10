"""Pure metadata/schema/overlay tests; no codec or benchmark execution."""
import copy
from pathlib import Path
import unittest
import tempfile
import contextlib
import io
import retain_paired_metadata as retention

from paired_metadata_overlay import adapt
from verify_paired_metadata import (BEAD, META, CELLS, expected, validate_metadata,
                                   end, inventory, parse, check_inventory, scope_command)
from verify_paired_r19932 import sha, equal


def fixture(profile='release', cell=8, schedule='0110', exercise=True):
    native=profile=='native'; k,r,size=CELLS[cell]
    old=expected(profile,cell,schedule,1,exercise)
    bias=0 if profile=='sanitize' else 0x500000000000
    base=0x400000 if profile=='sanitize' else 0
    image=dict(type=2 if profile=='sanitize' else 3,
        segments=[dict(address=base,bytes=0x20000,flags=5)],
        functions={name:dict(value=base+0x1000*(i+1),size=64) for i,name in enumerate(
            ('main','encode','metadata_anchor') if native else ('main','encode','batch','metadata_anchor'))})
    allocations={}
    for i,(name,amount) in enumerate(dict(source=k*size,reference=r*size,scratch=old['scratch_bytes'],output=r*size).items()):
        address=0x100000000+i*0x10000000
        allocations[name]=dict(raw=address,data=address+64,bytes=amount)
    if native: allocations['output']=dict(raw=0,data=0,bytes=0)
    parity=allocations['scratch' if native else 'output']['data']
    count=old['scratch_bytes']//size if native else r
    functions={n:bias+s['value'] for n,s in image['functions'].items()}
    if native: functions['batch']=0
    snapshot=dict(endpoint=0,allocations=allocations,
        spans=dict(parity=dict(address=parity,bytes=r*size),input_array=dict(address=0x200000000,bytes=k*8),
                   output_array=dict(address=0x300000000,bytes=count*8)),
        inputs=[allocations['source']['data']+i*size for i in range(k)],
        outputs=[parity+i*size for i in range(count)],functions=functions,
        dladdr_image_base=bias+base,load_bias=bias,page_bytes=4096,
        segments=[dict(address=bias+base,bytes=0x20000,flags=5)])
    other=copy.deepcopy(snapshot); other['endpoint']=1
    selections=[]
    for i in range(104 if exercise else 4):
        slot=i%4
        selections.append({'index':i,'phase':'preflight' if i<4 else 'exercise','pass':-1 if i<4 else (i-4)//4,
            'slot':slot,'requested_backend':-1 if native else 3 if cell==6 else 0,'context_backend':-1 if native else 3,
            'candidate_state':-1 if native else int(schedule[slot]),
            'operation_gfni':-1 if native else int(2<=cell<=4 or (cell<2 and schedule[slot]=='1'))})
    meta=dict(schema=META,bead=BEAD,timed=False,observation='new_frontend_endpoints',pointer_bytes=8,
        preflight_gfni_counts=old['probes'],selections=selections,
        native_route_label='original_native_compiler_policy' if native else 'not_native',snapshots=[snapshot,other])
    return meta,image


class MetadataTests(unittest.TestCase):
    def test_all_profiles_cells_schedules(self):
        for profile in ('native','release','sanitize'):
            for cell in range(9):
                for schedule in (('NNNN',) if profile=='native' else ('0110','1001','0000','1111')):
                    for exercise in (False,True):
                        meta,image=fixture(profile,cell,schedule,exercise)
                        validate_metadata(meta,profile,cell,schedule,exercise,image)

    def test_every_leaf_is_validated(self):
        meta,image=fixture()
        def leaves(value,path=()):
            if type(value) is dict:
                for k,v in value.items(): yield from leaves(v,path+(k,))
            elif type(value) is list:
                for k,v in enumerate(value): yield from leaves(v,path+(k,))
            else: yield path
        # Mutate every recorded field, including every one of the 104 selections.
        for path in leaves(meta):
            changed=copy.deepcopy(meta); at=changed
            for step in path[:-1]: at=at[step]
            at[path[-1]]=None
            with self.subTest(path=path),self.assertRaises((ValueError,TypeError)):
                validate_metadata(changed,'release',8,'0110',True,image)

    def test_consistent_endpoint_corruption(self):
        meta,image=fixture()
        for field in ('load_bias','dladdr_image_base','page_bytes'):
            bad=copy.deepcopy(meta)
            for snapshot in bad['snapshots']: snapshot[field]+=1
            with self.assertRaises(ValueError): validate_metadata(bad,'release',8,'0110',True,image)
        for name in ('source','scratch','reference','output'):
            bad=copy.deepcopy(meta)
            for snapshot in bad['snapshots']: snapshot['allocations'][name]['bytes']+=64
            with self.assertRaises(ValueError): validate_metadata(bad,'release',8,'0110',True,image)

    def test_pointer_arrays_and_native_alias(self):
        for profile,schedule in (('release','0110'),('native','NNNN')):
            meta,image=fixture(profile,7,schedule)
            self.assertEqual(len(meta['snapshots'][0]['outputs']),1024 if profile=='native' else 512)
            for field in ('inputs','outputs'):
                bad=copy.deepcopy(meta)
                for snapshot in bad['snapshots']: snapshot[field][-1]+=8
                with self.assertRaises(ValueError): validate_metadata(bad,profile,7,schedule,True,image)
            bad=copy.deepcopy(meta)
            for snapshot in bad['snapshots']: snapshot['spans']['parity']['address']+=64
            with self.assertRaises(ValueError): validate_metadata(bad,profile,7,schedule,True,image)

    def test_capacity_and_missing_extra_fields(self):
        meta,image=fixture()
        for field in ('selections','snapshots'):
            for remove in (True,False):
                bad=copy.deepcopy(meta)
                if remove: bad[field].pop()
                else: bad[field].append(copy.deepcopy(bad[field][-1]))
                with self.assertRaises(ValueError): validate_metadata(bad,'release',8,'0110',True,image)
        for field in meta:
            bad=copy.deepcopy(meta); del bad[field]
            with self.assertRaises(ValueError): validate_metadata(bad,'release',8,'0110',True,image)
        bad=copy.deepcopy(meta); bad['cause']='ASLR'
        with self.assertRaises(ValueError): validate_metadata(bad,'release',8,'0110',True,image)

    def test_typed_numbers_and_overflow(self):
        self.assertEqual(end(2**64-2,1),2**64-1)
        for a,b in ((2**64-1,1),(-1,0),(1,-1),(True,1),(1,1.0),(2**64,0)):
            with self.assertRaises(ValueError): end(a,b)
        for value in ('{"x":1,"x":2}','{"x":NaN}','{"x":Infinity}'):
            with self.assertRaises(ValueError): parse(value)

    def test_inventory_and_auto_request(self):
        rows=inventory()
        self.assertEqual(sum(r['code']==0 for r in rows),270)
        self.assertEqual(sum(r['code']==86 for r in rows),18)
        self.assertEqual(sum(bool(r['fault']) for r in rows),12)
        self.assertEqual(sum('-bad-' in r['label'] for r in rows),42)
        meta,image=fixture(); meta['selections'][0]['requested_backend']=3
        with self.assertRaises(ValueError): validate_metadata(meta,'release',8,'0110',True,image)

    def test_overlay_fails_closed(self):
        source=(Path(__file__).parent/'paired_timer_r19932.cpp').read_text()
        out=adapt(source)
        self.assertEqual(out.count('snapshot(0);'),1); self.assertEqual(out.count('snapshot(1);'),1)
        self.assertLess(out.index('metadata refuses real timing'),out.index('Buffer source('))
        self.assertIn('options.backend = index == 6 ? LEO2_BACKEND_AVX2 : LEO2_BACKEND_AUTO;',out)
        self.assertEqual(out.count('paired_metadata::Select('),2)
        for change in (source+'\n',source.replace('{17,7,64}','{17,7,128}'),source.replace('encode();','encode(); encode();')):
            with self.assertRaises(ValueError): adapt(change)

    def test_inventory_omissions(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); path=root/'abort'; path.write_bytes(b'pinned executable fixture')
            pins={'abort':sha(path)}
            check_inventory(root,pins)
            with self.assertRaises(ValueError): check_inventory(root,{})
            with self.assertRaises(ValueError): check_inventory(root,{'other':pins['abort']})
            with self.assertRaises(ValueError): check_inventory(root,dict(pins,extra=pins['abort']))

    def test_scope_command_mutations(self):
        command=scope_command('/tmp/qualification','release','abort',['--check','8','0110','256'])
        for original in ('/tmp/leopard-gf8-authoritative.lock','--cpu=60:60','--core=0:0',
                         'MemoryMax=256M','MemorySwapMax=0','120','/tmp/qualification/build/release/abort'):
            bad=list(command); bad[bad.index(original)]='wrong'
            with self.assertRaises(ValueError): equal(bad,scope_command('/tmp/qualification','release','abort',
                                                                       ['--check','8','0110','256']))


class MetadataRetentionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='paired-metadata-retention-test-')
        self.root=Path(self.temp.name); self.source=self.root/'source'; self.source.mkdir()
        self.output=self.root/'output'; self.output.mkdir(mode=0o700)
        for name in ('checks','units'):
            (self.source/name).mkdir()
            (self.source/name/'checks.json').write_text('{"completed":true,"timed":false}\n')
        for mode in ('normal','optimized'): (self.source/('replay-'+mode+'.json')).write_text('{}\n')
        (self.source/'artifact').write_bytes(b'qualification fixture')

    def tearDown(self):
        for path in self.root.rglob('*'):
            if path.is_dir() and not path.is_symlink(): path.chmod(0o700)
        self.temp.cleanup()

    def test_readonly_private_copy(self):
        output=io.StringIO()
        with contextlib.redirect_stdout(output): retention.retain(self.source,self.output)
        self.assertEqual(parse(output.getvalue())['files'],6)
        for source in self.source.rglob('*'):
            if source.is_file():
                target=self.output/source.relative_to(self.source)
                self.assertEqual(sha(source),sha(target))
                self.assertEqual(target.stat().st_mode&0o777,0o444)
                self.assertNotEqual((source.stat().st_dev,source.stat().st_ino),(target.stat().st_dev,target.stat().st_ino))
        self.assertEqual(self.output.stat().st_mode&0o777,0o555)

    def test_reserved_manifest_and_nonempty_destination(self):
        (self.source/'SHA256SUMS').write_text('preserve')
        with self.assertRaises(ValueError): retention.retain(self.source,self.output)
        self.assertFalse(any(self.output.iterdir()))
        (self.output/'owned').write_text('preserve')
        with self.assertRaises(ValueError): retention.retain(self.source,self.output)
        self.assertEqual((self.output/'owned').read_text(),'preserve')

    def test_incomplete_or_disagreeing_results(self):
        (self.source/'units/checks.json').write_text('{"completed":false,"timed":false}\n')
        with self.assertRaises(ValueError): retention.retain(self.source,self.output)
        (self.source/'units/checks.json').write_text('{"completed":true,"timed":false}\n')
        (self.source/'replay-optimized.json').write_text('{"different":true}\n')
        with self.assertRaises(ValueError): retention.retain(self.source,self.output)

    def test_links_and_overlap(self):
        alias=self.root/'alias'; alias.symlink_to(self.source,target_is_directory=True)
        with self.assertRaises(ValueError): retention.retain(alias,self.output)
        nested=self.source/'nested'; nested.mkdir(mode=0o700)
        with self.assertRaises(ValueError): retention.retain(self.source,nested)
        (self.source/'link').symlink_to(self.source/'artifact')
        with self.assertRaises(ValueError): retention.retain(self.source,self.output)


if __name__=='__main__': unittest.main()
