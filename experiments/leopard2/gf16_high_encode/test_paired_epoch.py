"""Pure epoch-overlay and accounting tests; never execute a codec or real clock."""
from pathlib import Path
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import paired_epoch_overlay as overlay
import verify_paired_epoch as verify

ROOT = Path(__file__).resolve().parent

class OverlayTests(unittest.TestCase):
    def test_driver(self):
        original = (ROOT/'paired_timer_r19932.cpp').read_text()
        result = overlay.driver(original)
        self.assertIn('samples.reserve(252);',result)
        self.assertIn('selections == (exercise ? 312U : 12U)',result)
        self.assertIn('epoch == 0 && slot == 0',result)
        self.assertIn('paired_metadata::RequireQuiescentProbe(probes[3]);',result)
        self.assertEqual(result.count('for (unsigned epoch = 0; epoch < 3; ++epoch)'),2)
        self.assertLess(result.index('metadata refuses real timing'), result.index('Buffer source('))
        self.assertIn('index == 6 ? LEO2_BACKEND_AVX2 : LEO2_BACKEND_AUTO',result)
        with self.assertRaises(ValueError): overlay.driver(original+'\n')

    def test_header(self):
        result = overlay.header((ROOT/'PairedRuntimeMetadata.h').read_text())
        self.assertIn('snapshots[6]',result)
        self.assertIn('kSelections = 312',result)
        self.assertIn('if (endpoint != 0) Require(Same(records.snapshots[0], s)',result)
        self.assertIn('local < 4',result)
        self.assertIn('selections == 12 || selections == 312',result)
        self.assertIn('paired-epoch-runtime-metadata/v1',result)
        self.assertIn('paired-epoch-progress/v1',result)
        self.assertLess(result.index('static Store records'), result.index('static Progress progress'))
        self.assertIn('!d::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics()',result)

    def test_clock(self):
        result = overlay.clock((ROOT/'paired_timer_clock.cpp').read_text())
        self.assertIn('calls[504]',result)
        self.assertIn('trace.count >= 504',result)
        self.assertIn('index / 168 == target_epoch',result)
        self.assertIn('LEO_PAIRED_TEST_CLOCK_FAULT_EPOCH',result)

    def test_witness(self):
        source = (ROOT/'paired_timer_witness.cpp').read_text()
        result = overlay.witness(source)
        self.assertTrue(result.startswith(source))
        self.assertIn('Mark marks[6]',result)
        self.assertNotIn('witness =',result)
        self.assertIn('endpoint >= 6',result)
        for name, fn in (('paired_timer_witness.cpp',overlay.witness),
                         ('paired_timer_clock.cpp',overlay.clock),('PairedRuntimeMetadata.h',overlay.header)):
            with self.assertRaises(ValueError): fn((ROOT/name).read_text()+'\n')

class OracleTests(unittest.TestCase):
    def test_inventory(self):
        rows = verify.inventory()
        self.assertEqual(sum(r['code']==0 for r in rows),270)
        self.assertEqual(sum(r['code']==86 for r in rows),18)
        self.assertEqual(sum(r['fault'] is not None for r in rows),36)
        self.assertEqual(sum('-bad-' in r['label'] for r in rows),42)
        for epoch in range(3): self.assertEqual(sum(r['fault_epoch']==epoch for r in rows),12)

    def test_public_counts_and_marks(self):
        for p in verify.ARCHIVES:
            for g in (1,256):
                s = 'NNNN' if p=='native' else '0110'
                for exercise in (False, True):
                    c = 4+100*g if exercise else 4
                    w,m = verify.public_oracle(p,8,s,g,exercise)
                    self.assertEqual(w['calls'],3*c)
                    self.assertEqual([r['calls'] for r in m['marks']],[0,c,c,2*c,2*c,3*c])
                    self.assertEqual(w['calls'],sum(w['states']))
                    self.assertEqual(w['calls'],sum(w['apis']))
                for e in range(3):
                    w,m = verify.public_oracle(p,8,s,g,True,e,'fault')
                    self.assertEqual(w['calls'],e*(4+100*g)+4+17*g)
                    self.assertEqual(len(m['marks']),2*e+1)
                w,m = verify.public_oracle(p,8,s,g,True,0,'abort')
                self.assertEqual(w['calls'],4+16*g)
                self.assertEqual(len(m['marks']),1)

    def test_clock_epoch_boundaries(self):
        for g,starts in ((1,[20,124,228]),(256,[4100,29704,55308])):
            row = verify.clocks(g)
            self.assertEqual(row['clock_calls'],504)
            self.assertEqual(row['public_calls_at_clock'][::168],starts)
            self.assertEqual(row['public_calls_at_clock'][-1],3*(4+100*g))
        for n in (-1,505,True):
            with self.assertRaises(ValueError): verify.clocks(1,n)

    def test_native_order_oracle_independently(self):
        # All native calls carry symbol10: state2 plus4*API2. Recompute the
        # twelve-call check digest without the oracle's epoch/group loops.
        digest = 14695981039346656037
        for _ in range(12): digest = ((digest^10)*1099511628211) % (2**64)
        row,_ = verify.public_oracle('native',0,'NNNN',1,False)
        self.assertEqual(row['order_hash'],f'{digest:016x}')

    def test_duplicate_and_nonfinite_json(self):
        for value in ('{"a":1,"a":2}','{"a":NaN}','{"a":Infinity}'):
            with self.assertRaises(ValueError): verify.parse(value)
        with self.assertRaises(ValueError): verify.records('{"schema":"x"}\n{"schema":"x"}')

    def test_missing_local_import_and_digest(self):
        with tempfile.TemporaryDirectory() as raw:
            folder=Path(raw)
            (folder/'entry.py').write_text('import helper\n')
            (folder/'helper.py').write_text('VALUE = 1\n')
            pins={p.name:verify.sha(p) for p in folder.iterdir()}
            required={'entry.py','helper.py'}
            verify.verify_tools(folder,pins,required,['entry.py'])
            (folder/'helper.py').unlink(); pins.pop('helper.py')
            # The former optional-discovery closure now reports just entry.py;
            # exact required names must reject the coordinated removal.
            self.assertEqual(verify.tool_closure(folder,['entry.py']),{'entry.py'})
            with self.assertRaises(ValueError): verify.verify_tools(folder,pins,required,['entry.py'])

    def test_fault_record_mutations(self):
        for row in [r for r in verify.inventory() if r['fault']]:
            _,c,s,g = row['args']; c,g = int(c),int(g); e = row['fault_epoch']
            w,m = verify.public_oracle(row['profile'],c,s,g,True,e,'fault')
            records = [w,m,verify.clocks(g,168*e+2),dict(schema=verify.PROGRESS,
                       selection_count=104*e+21,snapshot_count=2*e+1,timed=False)]
            error = ('group duration exceeds exact binary64 integer range' if row['fault']=='huge'
                     else 'nonpositive or reversed grouped clock interval')+'\n'
            def check(value): verify.verify_record(row,'\n'.join(json.dumps(r) for r in value),error,None)
            check(records)
            mutations=[]
            bad=copy.deepcopy(records); bad[-1]['selection_count']-=1; mutations.append(bad)
            bad=copy.deepcopy(records); bad[0]['calls']=0; mutations.append(bad)
            bad=copy.deepcopy(records); bad[1]['marks'].pop(); mutations.append(bad)
            bad=copy.deepcopy(records); bad[2]['public_calls_at_clock'][-1]-=1; mutations.append(bad)
            mutations.append(records[:-1])
            mutations.append(records+[dict(schema=verify.DRIVER,timed=False)])
            for bad in mutations:
                with self.assertRaises(ValueError): check(bad)

class FinalToolTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name); self.folder=self.root/'final_tools'; self.folder.mkdir()
        self.pins={}
        baseline=self.root/'build/baseline'; baseline.mkdir(parents=True)
        units=self.root/'units'; units.mkdir()
        for name in verify.FINAL_PYTHON | set(verify.BASE_NAMES) | set(verify.UNIT_SOURCES):
            if name in verify.FINAL_PYTHON: shutil.copyfile(ROOT/name,self.folder/name)
            else: (self.folder/name).write_text('fixture asset: '+name+'\n')
            self.pins[name]=verify.sha(self.folder/name)
            if name in verify.BASE_NAMES: shutil.copyfile(self.folder/name,baseline/name)
            if name in verify.UNIT_SOURCES: shutil.copyfile(self.folder/name,units/name)
        (self.root/'build/build.json').write_text(json.dumps(dict(
            baseline={n:self.pins[n] for n in verify.BASE_NAMES})))
        (units/'build.json').write_text(json.dumps(dict(
            source_sha256={n:self.pins[n] for n in verify.UNIT_SOURCES})))
        self.manifest()

    def manifest(self):
        (self.root/'final-tools.json').write_text(json.dumps(dict(bead=verify.BEAD,timed=False,files=self.pins)))

    def changed(self,name):
        path=self.folder/name; path.write_text(path.read_text()+'\n')
        self.pins[name]=verify.sha(path); self.manifest()

    def test_complete_final_inventory(self):
        verify.verify_final_tools(self.root)
        self.assertEqual(len(self.pins),24)

    def test_missing_dependency_and_digest(self):
        name='verify_auto_gfni_boundary_checks.py'
        (self.folder/name).unlink(); self.pins.pop(name); self.manifest()
        with self.assertRaises(ValueError): verify.verify_final_tools(self.root)

    def test_changed_asset_with_updated_digest(self):
        for name in ('paired_timer_r19932.cpp','test_paired_epoch_native.cpp'):
            original=(self.folder/name).read_bytes(); self.changed(name)
            with self.assertRaises(ValueError): verify.verify_final_tools(self.root)
            (self.folder/name).write_bytes(original); self.pins[name]=verify.sha(self.folder/name); self.manifest()

    def test_executing_verifier_identity(self):
        self.changed('verify_paired_epoch.py')
        with self.assertRaises(ValueError): verify.verify_final_tools(self.root)

    def executing_check(self, source, extra=''):
        program = ('import sys; from pathlib import Path; '
                   'sys.path.insert(0, sys.argv[1]); '
                   'import verify_paired_epoch as verify; import verify_paired_epoch_units; '
                   + extra + 'verify.verify_final_tools(Path(sys.argv[2]))')
        return subprocess.run([sys.executable, '-I', '-B', *([] if __debug__ else ['-O']),
                               '-c', program, str(source), str(self.root)],
                              capture_output=True, text=True, timeout=10)

    def test_drifted_executing_dependency(self):
        source=self.root/'executing'; shutil.copytree(self.folder,source)
        self.assertEqual(self.executing_check(source).returncode,0)
        for name in ('verify_paired_epoch_units.py','verify_paired_metadata.py','paired_epoch_overlay.py'):
            with self.subTest(name=name):
                original=(source/name).read_text()
                (source/name).write_text(original+'\nDRIFTED_DEPENDENCY = True\n')
                result=self.executing_check(source)
                self.assertNotEqual(result.returncode,0)
                self.assertIn('executing dependency hash: '+name,result.stderr)
                (source/name).write_text(original)

    def test_drifted_duplicate_main_module(self):
        source=self.root/'executing'; shutil.copytree(self.folder,source)
        alias=source/'alias'; alias.mkdir()
        (alias/'verify_paired_epoch.py').write_text(
            (source/'verify_paired_epoch.py').read_text()+'\nDRIFTED_DEPENDENCY = True\n')
        extra = ("import importlib.util; "
                 "spec=importlib.util.spec_from_file_location('epoch_alias', "
                 "Path(sys.argv[1])/'alias/verify_paired_epoch.py'); "
                 "module=importlib.util.module_from_spec(spec); "
                 "sys.modules['epoch_alias']=module; spec.loader.exec_module(module); ")
        result=self.executing_check(source,extra)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('executing dependency hash: verify_paired_epoch.py',result.stderr)

    def test_unlisted_directory(self):
        (self.folder/'unlisted').mkdir()
        with self.assertRaises(ValueError): verify.verify_final_tools(self.root)


class RetainedRecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw=os.environ.get('LEO_PAIRED_EPOCH_EVIDENCE')
        if not raw: raise unittest.SkipTest('real-record tests require LEO_PAIRED_EPOCH_EVIDENCE')
        cls.root=Path(raw).resolve(strict=True)
        state=verify.parse((cls.root/'checks/checks.json').read_text())
        verify.equal(state['completed'],True)

    def record(self,profile='release',cell=8,mode='--clock-exercise'):
        row=next(r for r in verify.inventory() if r['profile']==profile and r['code']==0 and
                 r['args']==[mode,str(cell),'NNNN' if profile=='native' else '0110','1'])
        label=row['label']; folder=self.root/'checks'
        output,error,_=verify.split_scope((folder/(label+'.stdout')).read_text(),(folder/(label+'.stderr')).read_text(),0)
        image=verify.elf(self.root/'build'/profile/row['binary'],profile=='native')
        verify.verify_record(row,output,error,image)
        return row,verify.records(output),image

    def reject(self,row,records,image):
        with self.assertRaises(ValueError):
            verify.verify_record(row,'\n'.join(json.dumps(r) for r in records.values()),'',image)

    def test_epoch_accounting_mutations(self):
        row,good,image=self.record()
        mutations=[]
        bad=copy.deepcopy(good); bad[verify.DRIVER]['encode_calls']//=3; mutations.append(bad)
        bad=copy.deepcopy(good); bad[verify.ACCOUNT]['epochs'].pop(); mutations.append(bad)
        bad=copy.deepcopy(good); bad[verify.ACCOUNT]['epochs'][1]['sample_begin']=0; mutations.append(bad)
        bad=copy.deepcopy(good); bad[verify.MARKS]['marks'][2]['calls']=0; mutations.append(bad)
        bad=copy.deepcopy(good); bad[verify.MARKS]['marks'][2:4]=reversed(bad[verify.MARKS]['marks'][2:4]); mutations.append(bad)
        bad=copy.deepcopy(good); bad[verify.WITNESS]['order_hash']='0000000000000000'; mutations.append(bad)
        bad=copy.deepcopy(good); bad.pop(verify.PROGRESS); mutations.append(bad)
        bad=copy.deepcopy(good); bad[verify.PROGRESS]['snapshot_count']=True; mutations.append(bad)
        for endpoint in (168,336):
            bad=copy.deepcopy(good); bad[verify.CLOCK]['public_calls_at_clock'][endpoint]-=4; mutations.append(bad)
        for bad in mutations: self.reject(row,bad,image)

    def test_check_and_exercise_metadata_mutations(self):
        for mode in ('--check','--exercise','--clock-exercise'):
            row,good,image=self.record(mode=mode); stride=4 if mode=='--check' else 104
            mutations=[]
            bad=copy.deepcopy(good); bad[verify.META]['selections'][stride]['phase']='exercise'; mutations.append(bad)
            bad=copy.deepcopy(good); bad[verify.META]['selections'][stride]['index']=0; mutations.append(bad)
            bad=copy.deepcopy(good); bad[verify.META]['preflight_gfni_counts'][1][0]=1; mutations.append(bad)
            bad=copy.deepcopy(good); bad[verify.META]['snapshots'][2]['epoch']=0; mutations.append(bad)
            bad=copy.deepcopy(good); bad[verify.META]['snapshots'].pop(); mutations.append(bad)
            for bad in mutations: self.reject(row,bad,image)

    def test_later_pair_equal_but_cross_epoch_changed(self):
        row,good,image=self.record(mode='--check')
        bad=copy.deepcopy(good)
        for snapshot in bad[verify.META]['snapshots'][2:4]:
            snapshot['spans']['input_array']['address']+=8
        meta=bad[verify.META]
        # Both changed snapshots still form a valid pair under the old
        # two-endpoint geometry checker. Only cross-epoch equality rejects it.
        normal=copy.deepcopy(meta)
        normal.update(schema=verify.prior.META,bead=verify.prior.BEAD,observation='new_frontend_endpoints',
                      preflight_gfni_counts=meta['preflight_gfni_counts'][1],selections=[],snapshots=[])
        for i,selection in enumerate(meta['selections'][4:8]):
            selection=dict(selection); selection.pop('epoch'); selection['index']=i
            normal['selections'].append(selection)
        for i,snapshot in enumerate(meta['snapshots'][2:4]):
            normal['snapshots'].append(dict({k:v for k,v in snapshot.items() if k not in ('epoch','phase','endpoint')},endpoint=i))
        verify.prior.validate_metadata(normal,'release',8,'0110',False,image)
        self.reject(row,bad,image)

    def test_native_full_output_geometry(self):
        row,good,image=self.record(profile='native',cell=7,mode='--check')
        self.assertEqual(len(good[verify.META]['snapshots'][0]['outputs']),1024)
        bad=copy.deepcopy(good); bad[verify.META]['snapshots'][5]['outputs'].pop()
        self.reject(row,bad,image)

if __name__ == '__main__': unittest.main()
