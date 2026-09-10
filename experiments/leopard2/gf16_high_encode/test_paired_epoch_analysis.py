"""Synthetic analysis/mutation tests; no clocks or native processes."""
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest

import paired_epoch_analysis as collect
import replay_paired_epoch_analysis as replay
import test_paired_r19932_screen as prior_tests


def fixture():
    return [prior_tests.fixtures() for _ in range(3)]


def checked(rows):
    actual,derived = collect.analyze_projected(rows),replay.derive_projected(rows)
    replay.compare(actual,derived)
    return actual


class AnalysisTests(unittest.TestCase):
    def test_complete_inventory_and_diagnostic_only(self):
        self.assertEqual(list(collect.inherited.schedule()),list(replay.schedule()))
        result=checked(fixture())
        self.assertTrue(result['all_epochs_pass'])
        self.assertFalse(result['production_promotion'])
        self.assertFalse(result['default_enabled'])
        self.assertEqual(len(result['homogeneous_processes']),264)
        self.assertEqual([e['epoch'] for e in result['epochs']],[0,1,2])
        self.assertEqual(sum(k.startswith(('same_','within_')) for e in result['epochs']
            for c in e['cells'] for k in c['ratios']),300)
        for e in result['epochs']:
            self.assertEqual(e['decision'],'diagnostic_all_gates_pass')
            self.assertTrue(e['diagnostic_only'])

    def test_every_epoch_retained_not_pooled(self):
        for epoch_id in range(3):
            rows=fixture()
            for row in rows[epoch_id]:
                if row['cell']==8 and row['comparison']=='same_off' and row['slot'] in (0,3):
                    prior_tests.samples(row,[v[1]*1.10 for v in row['record']['samples']])
            result=checked(rows)
            self.assertFalse(result['all_epochs_pass'])
            self.assertEqual([r['decision'] for r in result['epochs']],
                ['inconclusive_controls' if i==epoch_id else 'diagnostic_all_gates_pass' for i in range(3)])

    def test_all_300_control_locations(self):
        # Shared inherited estimators are already independently tested. Check
        # every per-epoch control remains represented and individually gated;
        # only the changed epoch needs its estimators reevaluated each case.
        pristine=fixture()
        for e in range(3):
            keys={(r['cell'],r['comparison'],r['slot']) for r in pristine[e] if r['comparison'].startswith('same_')}
            cross={(c,k) for c,k,_ in keys}
            self.assertEqual((len(keys),len(cross)),(80,20))
            for cell,comp,slot in sorted(keys) + [(c,k,None) for c,k in sorted(cross)]:
                rows=copy.deepcopy(pristine[e])
                for row in rows:
                    if row['cell']!=cell or row['comparison']!=comp: continue
                    if slot is None and row['slot'] in (0,3):
                        prior_tests.samples(row,[v[1]*1.1 for v in row['record']['samples']])
                    elif row['slot']==slot:
                        prior_tests.samples(row,[v[1]*(1.1 if i%4 in (0,3) else 1/1.1)
                                                for i,v in enumerate(row['record']['samples'])])
                actual,derived=collect.inherited.analyze(rows),replay.prior.derive(rows)
                replay.prior.compare_analysis(actual,derived)
                self.assertEqual(actual['decision'],'inconclusive_controls')

    def test_neighbors_targets_native_each_epoch(self):
        for e in range(3):
            for cell,comparison,outcome in ((8,'paired_1001','reject_neighbor_gate'),
                                           (0,'paired_0110','reject_target_gate'),(1,'native_on','reject_native_gate')):
                rows=fixture()
                for row in rows[e]:
                    if row['cell']!=cell or row['comparison']!=comparison: continue
                    if cell==8: values=[110000 if s=='0' else 100000 for s in row['order']]*21
                    else: values=[100000]*84
                    prior_tests.samples(row,values)
                result=checked(rows)
                self.assertEqual(result['epochs'][e]['decision'],outcome)
                self.assertFalse(result['all_epochs_pass'])

    def test_homogeneous_epoch_ratios_and_mixed_exclusion(self):
        rows=fixture()
        for e,factor in enumerate((1,2,3)):
            for row in rows[e]: prior_tests.samples(row,[v[1]*factor for v in row['record']['samples']])
        result=checked(rows)
        self.assertTrue(result['all_epochs_pass'])
        for row in result['homogeneous_processes']:
            self.assertIn(row['order'],('0000','1111','NNNN'))
            self.assertEqual(row['epoch1_over_epoch0'],2)
            self.assertEqual(row['epoch2_over_epoch0'],3)

    def test_missing_duplicate_and_reordered_processes(self):
        for change in ('epoch','process','order'):
            rows=fixture()
            if change=='epoch': rows.pop()
            elif change=='process': rows[2].pop()
            else: rows[1][0],rows[1][1]=rows[1][1],rows[1][0]
            for fn in (collect.analyze_projected,replay.derive_projected):
                with self.assertRaises(ValueError): fn(rows)


class FullRecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw=os.environ.get('LEO_PAIRED_EPOCH_EVIDENCE')
        if not raw: raise unittest.SkipTest('retained fixtures require LEO_PAIRED_EPOCH_EVIDENCE')
        cls.root=Path(raw).resolve(strict=True)

    def record(self):
        folder=self.root/'checks'; label='release-8-0110-256-clock-exercise'
        stdout,error,_=collect.epoch.split_scope((folder/(label+'.stdout')).read_text(),
                                                (folder/(label+'.stderr')).read_text(),0)
        records=collect.epoch.records(stdout)
        img=collect.epoch.elf(self.root/'build/release/synthetic',False)
        base=records.pop(collect.epoch.DRIVER)
        base.update(schema=collect.SCHEMA,clock_source='steady',timed=True)
        records[collect.SCHEMA]=base
        for key in (collect.epoch.WITNESS,collect.epoch.MARKS,collect.epoch.CLOCK): records.pop(key)
        return records,img

    def test_complete_before_projection(self):
        good,img=self.record()
        a=collect.projections(good,8,'0110',img); b=replay.projections(good,8,'0110',img)
        self.assertEqual(a,b)
        self.assertEqual([len(r['samples']) for r in a],[84,84,84])
        self.assertEqual(a[2]['samples'][0],good[collect.SCHEMA]['samples'][168])
        for e in range(3):
            for field in ('sample_begin','calls','selections'):
                bad=copy.deepcopy(good); bad[collect.epoch.ACCOUNT]['epochs'][e][field]+=1
                for fn in (collect.projections,replay.projections):
                    with self.assertRaises(ValueError): fn(bad,8,'0110',img)
        for key in good:
            bad=copy.deepcopy(good); bad.pop(key)
            for fn in (collect.projections,replay.projections):
                with self.assertRaises((ValueError,KeyError)): fn(bad,8,'0110',img)

    def test_fake_trace_metadata_and_samples_cannot_be_discarded(self):
        good,img=self.record(); mutations=[]
        bad=copy.deepcopy(good); bad[collect.epoch.CLOCK]={'schema':collect.epoch.CLOCK}; mutations.append(bad)
        bad=copy.deepcopy(good); bad[collect.epoch.META]['snapshots'][5]['epoch']=0; mutations.append(bad)
        bad=copy.deepcopy(good); bad[collect.epoch.META]['selections'][208]['index']=0; mutations.append(bad)
        bad=copy.deepcopy(good); bad[collect.SCHEMA]['samples'].pop(); mutations.append(bad)
        for value in (True,0,2**53):
            bad=copy.deepcopy(good); bad[collect.SCHEMA]['samples'][251][0]=value; mutations.append(bad)
        bad=copy.deepcopy(good); bad[collect.SCHEMA]['samples'][84][1]=float('nan'); mutations.append(bad)
        bad=copy.deepcopy(good); bad[collect.SCHEMA]['timed']=1; mutations.append(bad)
        for bad in mutations:
            for fn in (collect.projections,replay.projections):
                with self.assertRaises(ValueError): fn(bad,8,'0110',img)


class FileBackedTests(unittest.TestCase):
    """Shape/analysis fixtures, deliberately NOT authentic timing evidence."""
    @classmethod
    def setUpClass(cls):
        raw=os.environ.get('LEO_PAIRED_EPOCH_EVIDENCE')
        if not raw: raise unittest.SkipTest('retained fixtures require LEO_PAIRED_EPOCH_EVIDENCE')
        cls.root=Path(raw).resolve(strict=True)

    def make_files(self,folder):
        epochs=fixture(); entries=[]; images={}; templates={}
        for item in collect.inherited.schedule():
            profile='native' if item['order']=='NNNN' else 'release'
            if profile not in images:
                images[profile]=collect.epoch.elf(self.root/'build'/profile/'synthetic',profile=='native')
            key=(profile,item['cell'],item['order'])
            if key not in templates:
                group=256 if item['cell']==8 else 1
                name=f"{profile}-{item['cell']}-{item['order']}-{group}-clock-exercise"
                output,_,_=collect.epoch.split_scope((self.root/'checks'/(name+'.stdout')).read_text(),
                    (self.root/'checks'/(name+'.stderr')).read_text(),0)
                records=collect.epoch.records(output)
                base=records.pop(collect.epoch.DRIVER)
                base.update(schema=collect.SCHEMA,clock_source='steady',timed=True)
                records[collect.SCHEMA]=base
                for schema in (collect.epoch.CLOCK,collect.epoch.MARKS,collect.epoch.WITNESS): records.pop(schema)
                templates[key]=records
            records=copy.deepcopy(templates[key]); index=len(entries)
            records[collect.SCHEMA]['samples']=[s for e in epochs for s in e[index]['record']['samples']]
            name=collect.inherited.label(item)+'.stdout'; path=folder/name
            path.write_text('\n'.join(json.dumps(r,allow_nan=False) for r in records.values())+'\n')
            entries.append(dict(item,sibling_delta=0,stdout=name,stdout_sha256=collect.sha(path)))
        return entries,images

    def test_complete_file_backed_analysis_and_mutations(self):
        with tempfile.TemporaryDirectory(prefix='leopard-epoch-analysis-fixture-') as raw:
            folder=Path(raw); entries,images=self.make_files(folder)
            actual=collect.analyze(folder,entries,images); derived=replay.derive(folder,entries,images)
            replay.compare(actual,derived)
            self.assertTrue(actual['all_epochs_pass']); self.assertFalse(actual['production_promotion'])
            mutations=[]
            bad=copy.deepcopy(entries); bad.pop(); mutations.append(bad)
            bad=copy.deepcopy(entries); bad[1]=bad[0]; mutations.append(bad)
            bad=copy.deepcopy(entries); bad[0],bad[1]=bad[1],bad[0]; mutations.append(bad)
            for key,value in (('stdout','../elsewhere'),('stdout_sha256','0'*64),('sibling_delta',1)):
                bad=copy.deepcopy(entries); bad[0][key]=value; mutations.append(bad)
            for bad in mutations:
                for fn in (collect.analyze,replay.derive):
                    with self.assertRaises(ValueError): fn(folder,bad,images)
            first=folder/entries[0]['stdout']; first.rename(folder/'saved.stdout'); first.symlink_to('saved.stdout')
            for fn in (collect.analyze,replay.derive):
                with self.assertRaises(ValueError): fn(folder,entries,images)
            first.unlink(); (folder/'saved.stdout').rename(first)
            # Valid updated file digest cannot conceal corruption in the last
            # process's final epoch. Both complete streaming paths must refuse.
            last=folder/entries[-1]['stdout']; records=collect.epoch.records(last.read_text())
            records[collect.epoch.ACCOUNT]['epochs'][2]['calls']+=1
            last.write_text('\n'.join(json.dumps(r) for r in records.values())+'\n')
            entries[-1]['stdout_sha256']=collect.sha(last)
            for fn in (collect.analyze,replay.derive):
                with self.assertRaises(ValueError): fn(folder,entries,images)


if __name__ == '__main__': unittest.main()
