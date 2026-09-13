"""Synthetic full-campaign fixtures and refusal tests; no native execution.

The native metadata templates are retained check/synthetic records. Their ELF
placement fields are deliberately adapted for a TEST fixture, never evidence.
"""
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import paired_epoch_campaign as c
import freeze_paired_epoch_campaign as freeze
import run_paired_epoch_campaign as collector
import replay_paired_epoch_campaign as replay
import test_paired_epoch_analysis as old_tests

PREREG='b64efb7957e1d99bea70435dd509f3bf4e80a48f'  # Deliberately NOT a campaign preregistration.
TEST_QUALIFICATION={'synthetic_qualification_standin':True}
OLD=Path('/home/catid/leopard/.research/leopard-79h/paired-epoch-qualified.unzse9')


def write(path,value):
    path.write_text(json.dumps(value,sort_keys=True,allow_nan=False)+'\n')


def rewrite_pin(bundle,name,data):
    file=bundle/name; file.chmod(0o600); file.write_bytes(data)
    file.chmod(0o555 if name in ('native','current') else 0o444)
    pinfile=bundle/'pins.json'; pins=c.read(pinfile); pins['files'][name]=c.sha(file)
    pinfile.chmod(0o600); write(pinfile,pins); pinfile.chmod(0o444)


def metadata(records,image,native):
    for snapshot in records[c.qualified.epoch.META]['snapshots']:
        bias=snapshot['load_bias'] if image['type']==3 else 0
        snapshot['load_bias']=bias
        snapshot['segments']=[dict(s,address=bias+s['address']) for s in image['segments']]
        snapshot['dladdr_image_base']=bias+min(s['address']//4096*4096 for s in image['segments'])
        snapshot['functions']={name:bias+s['value'] for name,s in image['functions'].items()}
        if native: snapshot['functions']['batch']=0


def records_text(records):
    return '\n'.join(json.dumps(r,allow_nan=False) for r in records.values())+'\n'


def resource_text(bundle=None,output=None,qualification=False):
    bundle=bundle or c.ROOT/'frozen'; output=output or c.ROOT/'attempt1'
    return '\n'.join(['EPOCH_SCOPE=/user.slice/synthetic-test.scope',
        'EPOCH_WRAPPER='+str(bundle/'paired_epoch_scope.sh'),'EPOCH_COMMAND_BEGIN',
        *c.controller_command(bundle,output,qualification,PREREG),'EPOCH_COMMAND_END',
        'TOWER_CHILD_EXIT=0','memory.peak','120000000','memory.max','268435456',
        'memory.events','low 0','high 0','max 0','oom 0','oom_kill 0','oom_group_kill 0',
        'memory.swap.current','0','memory.swap.max','0','\tExit status: 0',''])


def scope():
    return dict(host=c.HOST,topology='26,90',controller_affinity=[0],memory_max=268435456,
                swap_max=0,initial_memory_events=[0]*6,cgroup='/user.slice/synthetic-test.scope')


class CampaignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory(prefix='leopard-epoch-campaign-TEST-')
        cls.root=Path(cls.tmp.name); cls.bundle=cls.root/'frozen'
        freeze.copy_bundle(cls.bundle,Path(__file__).resolve().parent)
        cls.root_patch=patch.object(c,'ROOT',cls.root); cls.root_patch.start()
        rewrite_pin(cls.bundle,c.PLAN,(json.dumps(c.plan(True))+'\n').encode())
        _,cls.pins,cls.images=c.inputs(cls.bundle)
        cls.fixtures=cls.root/'synthetic-inputs'; cls.fixtures.mkdir()
        helper=old_tests.FileBackedTests(); helper.root=OLD
        entries,_=helper.make_files(cls.fixtures)
        for entry in entries:
            path=cls.fixtures/entry['stdout']; records=c.qualified.epoch.records(path.read_text())
            native=entry['order']=='NNNN'; metadata(records,cls.images['native' if native else 'release'],native)
            path.write_text(records_text(records))
        for item in c.preflights():
            native=item['order']=='NNNN'; profile='native' if native else 'release'; cell=item['cell']
            group=256 if cell==8 else 1
            label=f"{profile}-{cell}-{item['order']}-{group}-check"
            output,error,_=c.qualified.epoch.split_scope((OLD/'checks'/(label+'.stdout')).read_text(),
                (OLD/'checks'/(label+'.stderr')).read_text(),0)
            c.equal(error,''); records=c.qualified.epoch.records(output)
            base=records.pop(c.qualified.epoch.DRIVER)
            base.update(schema=collector.analysis.SCHEMA,clock_source='steady')
            records[collector.analysis.SCHEMA]=base
            for name in (c.qualified.epoch.CLOCK,c.qualified.epoch.MARKS,c.qualified.epoch.WITNESS): records.pop(name,None)
            metadata(records,cls.images[profile],native)
            text=records_text(records); c.check_record(text,cell,item['order'],cls.images[profile])
            (cls.fixtures/(c.name(item)+'.stdout')).write_text(text)
        cls.resource=cls.root/'resource.log'; cls.resource.write_text(resource_text())
        cls.attempt=cls.root/'attempt1'
        def fake_run(command,stdout,stderr,env,timeout):
            # Actual collector launch writes the intent BEFORE this stand-in.
            intent=c.parse((cls.attempt/'launches.jsonl').read_text().splitlines()[-1])
            c.equal(intent['command'],command); c.equal(intent['environment'],env)
            label=intent['label']
            text='\n'.join(c.condition_lines())+'\n' if label.startswith('condition-') else (cls.fixtures/(label+'.stdout')).read_text()
            stdout.write(text.encode()); return SimpleNamespace(returncode=0)
        def passive(state,value): state['passive']=dict(before=123,after=123,elapsed_ns=10**10)
        def locks(): return [os.open('/dev/null',os.O_RDONLY),os.open('/dev/null',os.O_RDONLY)]
        with patch.object(collector,'preregistration'),patch.object(collector,'qualification_gate',return_value=TEST_QUALIFICATION), \
             patch.object(collector,'acquire_locks',side_effect=locks), \
             patch.object(collector,'scope_and_host',side_effect=scope),patch.object(collector,'check_passive',side_effect=passive), \
             patch.object(collector,'sibling_ticks',return_value=123),patch.object(collector.subprocess,'run',side_effect=fake_run):
            collector.run(cls.bundle,PREREG)
        cls.state=c.read(cls.attempt/'attempt.json')

    @classmethod
    def tearDownClass(cls):
        cls.root_patch.stop(); cls.tmp.cleanup()

    def replay(self):
        with patch.object(replay,'publication'),patch.object(replay,'qualification_gate',return_value=TEST_QUALIFICATION):
            return replay.verify(self.attempt,self.bundle,self.resource)

    def test_complete_actual_collector_and_independent_replayer(self):
        result=self.replay()
        self.assertEqual(result['timed_processes'],318)
        self.assertEqual(result['preflights'],27)
        self.assertTrue(result['analysis']['all_epochs_pass'])
        self.assertFalse(result['production_promotion'])
        self.assertEqual(len(result['analysis']['homogeneous_processes']),264)

    def test_clock_free_qualification_and_relocated_full_replay(self):
        bundle=Path(c.plan()['qualification_path']); shutil.copytree(self.bundle,bundle)
        rewrite_pin(bundle,c.PLAN,(json.dumps(c.plan(False))+'\n').encode())
        output=Path(c.plan()['qualification_output']); resource=Path(c.plan()['qualification_resource'])
        resource.write_text(resource_text(bundle,output,True))
        def fake_run(command,stdout,stderr,env,timeout):
            intent=c.parse((output/'launches.jsonl').read_text().splitlines()[-1])
            self.assertEqual(intent['command'],command); self.assertNotIn('--measure',command)
            label=intent['label']
            text='\n'.join(c.condition_lines())+'\n' if label.startswith('condition-') else (self.fixtures/(label+'.stdout')).read_text()
            stdout.write(text.encode()); return SimpleNamespace(returncode=0)
        def locks(): return [os.open('/dev/null',os.O_RDONLY),os.open('/dev/null',os.O_RDONLY)]
        with patch.object(collector,'acquire_locks',side_effect=locks),patch.object(collector,'scope_and_host',side_effect=scope), \
             patch.object(collector,'sibling_ticks',return_value=123),patch.object(collector.subprocess,'run',side_effect=fake_run), \
             patch.object(collector,'check_passive') as passive:
            freeze.check(bundle,output); passive.assert_not_called()
        original=replay.verify(output,bundle,resource,qualification=True)
        self.assertEqual(original['timed_processes'],0)
        self.assertEqual(original['preflights'],27)
        for path in output.iterdir(): path.chmod(0o444)
        output.chmod(0o555); resource.chmod(0o444)
        self.assertEqual(replay.qualification_gate(self.pins),original)
        # Matching readiness cannot bypass missing or writable evidence.
        resource.chmod(0o644)
        with self.assertRaisesRegex(ValueError,'readonly qualification'): replay.qualification_gate(self.pins)
        resource.chmod(0o444)
        relocated=self.root/'relocated-qualification'; shutil.copytree(output,relocated)
        self.assertEqual(replay.verify(relocated,bundle,resource,qualification=True),original)
        state=c.read(relocated/'attempt.json'); state['output_path']=str(self.root/'wrong-output')
        (relocated/'attempt.json').chmod(0o644); write(relocated/'attempt.json',state)
        with self.assertRaises(ValueError): replay.verify(relocated,bundle,resource,qualification=True)

    def test_journal_mutations(self):
        mutations=[]
        for key,value in [('complete',1),('plan_sha256','0'*64),('child_environment',{}),('failure','late error')]:
            bad=copy.deepcopy(self.state); bad[key]=value; mutations.append(bad)
        bad=copy.deepcopy(self.state); bad['scope']['initial_memory_events'][2]=1; mutations.append(bad)
        bad=copy.deepcopy(self.state); bad['passive']['elapsed_ns']=10**10-1; mutations.append(bad)
        bad=copy.deepcopy(self.state); bad['preflight'].pop(); mutations.append(bad)
        bad=copy.deepcopy(self.state); bad['invocations'][-1]['sibling_delta']=1; mutations.append(bad)
        bad=copy.deepcopy(self.state); bad['invocations'][0]['command'][2]='27'; mutations.append(bad)
        bad=copy.deepcopy(self.state); bad['conditions'][1]['returncode']=1; mutations.append(bad)
        bad=copy.deepcopy(self.state); bad['invocations'][-1]['timed_out']=True; mutations.append(bad)
        try:
            for index,bad in enumerate(mutations):
                with self.subTest(index=index):
                    write(self.attempt/'attempt.json',bad)
                    with self.assertRaises(ValueError): self.replay()
        finally: write(self.attempt/'attempt.json',self.state)

    def test_updated_digest_late_metadata_cannot_hide_corruption(self):
        row=self.state['invocations'][-1]; file=self.attempt/row['stdout']; saved=file.read_bytes()
        outcome=self.attempt/(row['label']+'.outcome.json'); saved_outcome=outcome.read_bytes()
        try:
            records=c.qualified.epoch.records(saved.decode()); records[c.qualified.epoch.ACCOUNT]['epochs'][2]['calls']+=1
            file.write_text(records_text(records)); bad=copy.deepcopy(self.state)
            bad['invocations'][-1]['stdout_sha256']=c.sha(file); write(self.attempt/'attempt.json',bad)
            result=c.parse(saved_outcome.decode()); result['stdout_sha256']=c.sha(file); write(outcome,result)
            with self.assertRaises(ValueError): self.replay()
        finally:
            file.write_bytes(saved); outcome.write_bytes(saved_outcome); write(self.attempt/'attempt.json',self.state)

    def test_files_intents_conditions_and_resource_mutations(self):
        for target in (self.attempt/'launches.jsonl',self.resource,
                       self.attempt/'condition-after.stdout',self.attempt/'condition-after.stderr'):
            original=target.read_bytes()
            try:
                target.write_bytes(original+b'corrupt\n')
                # A trailing resource line is not inherently invalid: mutate
                # the actual required event counter in that case.
                if target==self.resource: target.write_bytes(original.replace(b'oom 0',b'oom 1'))
                with self.assertRaises((ValueError,json.JSONDecodeError)): self.replay()
            finally: target.write_bytes(original)
        extra=self.attempt/'unexpected'; extra.mkdir()
        try:
            with self.assertRaises(ValueError): self.replay()
        finally: extra.rmdir()
        first=self.attempt/self.state['preflight'][0]['stdout']; backup=first.with_suffix('.saved')
        first.rename(backup); first.symlink_to(backup.name)
        try:
            with self.assertRaises(ValueError): self.replay()
        finally: first.unlink(); backup.rename(first)

    def test_false_readiness_precedes_all_acquisition(self):
        saved=(self.bundle/c.PLAN).read_bytes()
        try:
            rewrite_pin(self.bundle,c.PLAN,(json.dumps(c.plan(False))+'\n').encode())
            with patch.object(collector.subprocess,'run') as process,patch.object(collector,'preregistration') as publication, \
                 patch.object(collector,'acquire_locks') as locks,patch.object(collector,'check_passive') as passive, \
                 patch.object(collector,'qualification_gate') as qualification:
                with self.assertRaisesRegex(ValueError,'not qualified'): collector.run(self.bundle,PREREG)
                for fn in (process,publication,locks,passive,qualification): fn.assert_not_called()
        finally: rewrite_pin(self.bundle,c.PLAN,saved)

    def test_consumed_attempt_never_reused(self):
        original=(self.attempt/'attempt.json').read_bytes()
        with patch.object(collector,'preregistration'),patch.object(collector,'qualification_gate',return_value=TEST_QUALIFICATION), \
             patch.object(collector,'acquire_locks') as locks:
            with self.assertRaises(FileExistsError): collector.run(self.bundle,PREREG)
            locks.assert_not_called()
        self.assertEqual((self.attempt/'attempt.json').read_bytes(),original)

    def test_nonpublished_and_drifted_preregistration(self):
        for fn in (collector.preregistration,replay.publication):
            with patch.object(subprocess,'check_output',return_value=b'wrong source'), \
                 patch.object(subprocess,'run',side_effect=subprocess.CalledProcessError(1,['git'])):
                with self.assertRaises((ValueError,subprocess.CalledProcessError)): fn(self.bundle,PREREG)
            with self.assertRaises(ValueError): fn(self.bundle,'not-a-commit')
        with patch.object(collector,'preregistration',side_effect=ValueError('not pushed')), \
             patch.object(collector,'qualification_gate',return_value=TEST_QUALIFICATION), \
             patch.object(collector,'acquire_locks') as locks:
            with self.assertRaisesRegex(ValueError,'not pushed'): collector.run(self.bundle,PREREG)
            locks.assert_not_called()

    def test_immutable_source_and_real_binary_pins(self):
        binary=self.bundle/'current'; binary.chmod(0o444)
        try:
            with self.assertRaises(ValueError): c.inputs(self.bundle)
        finally: binary.chmod(0o555)
        for name in ('current','qualification.tools.json','verify_paired_metadata.py'):
            data=(self.bundle/name).read_bytes()
            try:
                replacement=b'fake binary' if name=='current' else data+b'\n'
                rewrite_pin(self.bundle,name,replacement)
                with self.assertRaises(ValueError): c.inputs(self.bundle)
            finally:
                rewrite_pin(self.bundle,name,data)
                if name=='current': (self.bundle/name).chmod(0o555)
        pins=copy.deepcopy(self.pins); pins['files']['run_paired_epoch_campaign.py']='0'*64
        with self.assertRaisesRegex(ValueError,'executing dependency'): c.executing_sources(pins)
        saved=collector.__file__
        try:
            collector.__file__=str(self.bundle/'paired_epoch_campaign.py')
            with self.assertRaises(ValueError): c.executing_sources(self.pins)
        finally: collector.__file__=saved

    def test_exact_frozen_namespace(self):
        extra=self.bundle/'unexpected'; self.bundle.chmod(0o755)
        try:
            extra.mkdir(); self.bundle.chmod(0o555)
            with self.assertRaises(ValueError): c.inputs(self.bundle)
        finally:
            self.bundle.chmod(0o755); extra.rmdir(); self.bundle.chmod(0o555)

    def test_independent_module_does_not_import_campaign_collector(self):
        text=(Path(__file__).parent/'replay_paired_epoch_campaign.py').read_text()
        self.assertNotIn('import run_paired_epoch_campaign',text)
        self.assertNotIn('import freeze_paired_epoch_campaign',text)


class LaunchTests(unittest.TestCase):
    def test_condition_tools_available_in_restricted_path(self):
        for program in ('grep','systemctl','docker'):
            self.assertIsNotNone(shutil.which(program,path=c.condition_env()['PATH']),program)

    def test_missing_qualification_refuses_before_publication_or_attempt(self):
        with tempfile.TemporaryDirectory(prefix='leopard-epoch-gate-TEST-') as raw, \
             patch.object(c,'ROOT',Path(raw)):
            root=Path(raw); value=c.plan(True); pins={'files':{c.PLAN:'0'*64}}
            # Supply only the timing input result; leave the actual qualification
            # gate to discover the missing readiness-false bundle itself.
            inputs=c.inputs
            def supplied(bundle,*args,**kwargs):
                if bundle==root/'frozen': return value,pins,{}
                return inputs(bundle,*args,**kwargs)
            with patch.object(c,'inputs',side_effect=supplied),patch.object(collector,'preregistration') as publication, \
                 patch.object(collector,'acquire_locks') as locks:
                with self.assertRaisesRegex(ValueError,'immutable bundle'): collector.run(root/'frozen',PREREG)
                publication.assert_not_called(); locks.assert_not_called()
            self.assertFalse((root/'attempt1').exists())
            plan_path=Path(freeze.__file__).resolve().parent/c.PLAN
            read=c.read
            def read_plan(path): return value if path==plan_path else read(path)
            with patch.object(c,'read',side_effect=read_plan),patch.object(freeze,'copy_bundle') as copying:
                with self.assertRaisesRegex(ValueError,'immutable bundle'): freeze.freeze()
                copying.assert_not_called()
            self.assertFalse((root/'frozen').exists())

    def test_late_validation_failure_records_incomplete_and_closes_leases(self):
        # Focused terminal-state test, not a replacement for the complete
        # 318-process synthetic fixture above. No child is executed here.
        with tempfile.TemporaryDirectory(prefix='leopard-epoch-terminal-TEST-') as raw, \
             patch.object(c,'ROOT',Path(raw)):
            root=Path(raw); bundle=root/'frozen'; value=c.plan(True)
            pins={'files':{c.PLAN:'0'*64}}; leases=[]
            def locks():
                leases.extend(os.open('/dev/null',os.O_RDONLY) for _ in range(2)); return leases
            def condition_run(command,stdout,stderr,env,timeout):
                stdout.write(('\n'.join(c.condition_lines())+'\n').encode())
                return SimpleNamespace(returncode=0)
            with patch.object(c,'inputs',return_value=(value,pins,{})), \
                 patch.object(collector,'qualification_gate',return_value=TEST_QUALIFICATION), \
                 patch.object(collector,'preregistration'),patch.object(collector,'acquire_locks',side_effect=locks), \
                 patch.object(collector,'scope_and_host',side_effect=scope),patch.object(c,'preflights',return_value=[]), \
                 patch.object(collector,'check_passive'),patch.object(collector.analysis.inherited,'schedule',return_value=[]), \
                 patch.object(collector.analysis,'analyze',return_value={'synthetic':True}), \
                 patch.object(collector.subprocess,'run',side_effect=condition_run), \
                 patch.object(c,'executing_sources',side_effect=ValueError('late executing drift')):
                with self.assertRaisesRegex(ValueError,'late executing drift'): collector.run(bundle,PREREG)
            state=c.read(root/'attempt1/attempt.json')
            self.assertIs(state['complete'],False)
            self.assertIn('late executing drift',state['failure'])
            self.assertEqual(state['analysis'],{'synthetic':True})
            self.assertEqual(len(leases),2)
            for fd in leases:
                with self.assertRaises(OSError): os.fstat(fd)

    def test_qualification_cannot_consume_reserved_attempt(self):
        for path in (c.ROOT/'attempt1',c.ROOT/'attempt1'/'nested',c.ROOT):
            with self.assertRaisesRegex(ValueError,'reserved timing attempt'): c.disjoint_attempt(path)
        bundle=Path(c.plan()['qualification_path'])
        with patch.object(c,'inputs',return_value=(c.plan(False),{},{})), \
             patch.object(Path,'mkdir') as mkdir,patch.object(collector,'acquire_locks') as locks:
            with self.assertRaisesRegex(ValueError,'reserved timing attempt'):
                freeze.check(bundle,c.ROOT/'attempt1')
            mkdir.assert_not_called(); locks.assert_not_called()

    def test_launch_errors_and_timeouts_have_durable_intent_and_outcome(self):
        for error in (FileNotFoundError('missing'),subprocess.TimeoutExpired(['unused'],60)):
            with tempfile.TemporaryDirectory(prefix='leopard-epoch-launch-TEST-') as raw:
                root=Path(raw)
                with (root/'launches.jsonl').open('x') as intents, \
                     patch.object(collector.subprocess,'run',side_effect=error), \
                     patch.object(collector,'sibling_ticks',return_value=123):
                    row=collector.launch(root,intents,'failure',['unused'],c.ENV)
                self.assertEqual(c.read(root/'failure.outcome.json'),row)
                self.assertEqual(c.parse((root/'launches.jsonl').read_text())['command'],['unused'])
                self.assertTrue(row['failure'] is not None or row['timed_out'])
                with self.assertRaises(ValueError): collector.successful(root,row)

    def test_bounded_reads_and_resource_counters(self):
        with tempfile.TemporaryDirectory(prefix='leopard-epoch-bounds-TEST-') as raw:
            path=Path(raw)/'raw'; path.write_bytes(b'x'*128)
            with self.assertRaisesRegex(ValueError,'bounded'): c.bounded_text(path,128)
            args=('/user.slice/synthetic-test.scope',c.controller_command(c.ROOT/'frozen',c.ROOT/'attempt1',False,PREREG),
                  c.ROOT/'frozen'/'paired_epoch_scope.sh')
            path.write_text(resource_text()); self.assertEqual(c.resource(path,*args),120000000)
            for before,after in [('max 0','max 1'),('268435456','536870912'),('memory.swap.max\n0','memory.swap.max\n1'),
                                 ('synthetic-test.scope','unrelated.scope'),('--cpu=600:600','--cpu=900:900'),
                                 ('TOWER_CHILD_EXIT=0','TOWER_CHILD_EXIT=0\nTOWER_CHILD_EXIT=1')]:
                path.write_text(resource_text().replace(before,after))
                with self.assertRaises(ValueError): c.resource(path,*args)


if __name__=='__main__': unittest.main()
