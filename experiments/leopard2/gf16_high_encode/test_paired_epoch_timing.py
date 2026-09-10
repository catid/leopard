"""Clock-free tests of the separate diagnostic adapter and binary policies."""
from pathlib import Path
import copy
import json
import os
import unittest
import paired_epoch_overlay as epoch
import paired_epoch_timing_overlay as overlay
import paired_epoch_timing_elf as image
import verify_paired_epoch_timing as verify

ROOT = Path(__file__).resolve().parent


class AdapterTests(unittest.TestCase):
    def test_driver_delta(self):
        original = epoch.driver((ROOT/'paired_timer_r19932.cpp').read_text())
        actual = overlay.driver(original)
        self.assertNotIn('metadata refuses real timing', actual)
        self.assertIn('unknown diagnostic clock binding', actual)
        for line in ('Require(!synthetic || !std::strcmp(clock_kind, "synthetic"), "synthetic clock required");',
                     'Require(!measured || !std::strcmp(clock_kind, "steady"), "steady clock required");',
                     'Require(!clock_guard || !std::strcmp(clock_kind, "abort"), "abort clock required");'):
            self.assertIn(line,actual)
        self.assertIn(overlay.SCHEMA,actual)
        self.assertIn('samples.reserve(252);',actual)
        self.assertIn('selections == (exercise ? 312U : 12U)',actual)
        self.assertEqual(actual.count('for (unsigned epoch = 0; epoch < 3; ++epoch)'),2)
        with self.assertRaises(ValueError): overlay.driver(original+'\n')
        with self.assertRaises(ValueError): overlay.driver(actual)

    def test_fake_clock_delta(self):
        original = epoch.clock((ROOT/'paired_timer_clock.cpp').read_text())
        actual = overlay.fake_steady_clock(original)
        self.assertEqual(actual,original.replace('return "synthetic";', 'return "steady";'))
        self.assertIn('calls[504]',actual)
        self.assertIn('index / 168 == target_epoch',actual)
        with self.assertRaises(ValueError): overlay.fake_steady_clock(original+'\n')

    def test_plain_mark_is_empty(self):
        self.assertEqual(overlay.PLAIN_MARK.splitlines()[-1], 'extern "C" void LeoPairedWitnessMark(unsigned) {}')


class BindingTests(unittest.TestCase):
    def symbols(self,native,variant):
        defined = {'LeoPairedWitnessMark','LeoPairedClockKind'}
        undefined = set()
        if variant == 'steady': undefined.add(image.CLOCK)
        else: defined.add('__wrap_'+image.CLOCK)
        if variant in ('synthetic','fake-steady'):
            defined.add('LeoPairedWitnessCalls')
            defined.update({'__wrap_leo_encode'} if native else {'__wrap_leo2_encode','__wrap_leo2_encode_batch'})
        return defined,undefined

    def test_all_bindings(self):
        for native in (False,True):
            for variant in image.VARIANTS:
                defined,undefined = self.symbols(native,variant)
                image.bindings(defined,undefined,native,variant)
                for name in defined:
                    with self.assertRaises(ValueError):
                        image.bindings(defined-{name},undefined,native,variant)

    def test_fake_cannot_be_real(self):
        for native in (False,True):
            for variant in ('abort','synthetic','fake-steady'):
                defined,undefined = self.symbols(native,variant)
                with self.assertRaises(ValueError): image.bindings(defined,undefined,native,'steady')
                with self.assertRaises(ValueError): image.bindings(defined,undefined|{image.CLOCK},native,variant)

    def test_real_cannot_be_fake_or_observed(self):
        for native in (False,True):
            defined,undefined = self.symbols(native,'steady')
            for variant in ('abort','synthetic','fake-steady'):
                with self.assertRaises(ValueError): image.bindings(defined,undefined,native,variant)
            for name in ('__wrap_leo_encode','LeoPairedWitnessCalls','__wrap_'+image.CLOCK,image.CLOCK):
                with self.assertRaises(ValueError): image.bindings(defined|{name},undefined,native,'steady')


class QualificationTests(unittest.TestCase):
    def test_inventory_and_clock_firewall(self):
        rows = verify.inventory()
        self.assertEqual(len(rows),189)
        self.assertEqual(sum(r['code']==0 for r in rows),126)
        self.assertEqual(sum(bool(r['fault']) for r in rows),36)
        self.assertEqual(sum(bool(r['reason']) for r in rows),24)
        for r in rows:
            command = verify.scope_command('/owned',r['profile'],r['variant'],r['args'])
            self.assertIn('/owned/build/'+verify.executable(r['profile'],r['variant']),command)
            if r['variant']=='steady': self.assertIn(r['args'][0],('--check','--exercise'))
        for mode in ('--measure','--clock-exercise','--clock-guard','--unknown'):
            with self.assertRaises(ValueError): verify.qualification_args('steady',[mode,'0','0110','1'])

    def test_failure_oracles(self):
        for row in verify.inventory():
            if row['code']==0: continue
            wanted = {verify.epoch.PROGRESS:dict(schema=verify.epoch.PROGRESS,selection_count=0,
                                                snapshot_count=0,timed=False)}
            p,v = row['profile'],row['variant']
            calls,marks = verify.epoch.public_oracle(p,0,'NNNN' if p=='native' else '0110',1,False,completed=0)
            if row['reason']:
                error = row['reason']+'\n'
                if v in ('synthetic','fake-steady'): wanted[verify.epoch.CLOCK] = verify.epoch.clocks(1,0)
            elif row['code']==86:
                error = 'unexpected driver benchmark clock\n'
                wanted[verify.epoch.PROGRESS].update(selection_count=21,snapshot_count=1)
            else:
                _,c,s,g = row['args']; c,g = int(c),int(g); e = row['fault_epoch']
                calls,marks = verify.epoch.public_oracle(p,c,s,g,True,e,'fault')
                wanted[verify.epoch.CLOCK] = verify.epoch.clocks(g,168*e+2)
                wanted[verify.epoch.PROGRESS].update(selection_count=104*e+21,snapshot_count=2*e+1)
                error = ('group duration exceeds exact binary64 integer range' if row['fault']=='huge' else
                         'nonpositive or reversed grouped clock interval')+'\n'
            if v in ('synthetic','fake-steady'): wanted.update({verify.epoch.WITNESS:calls,verify.epoch.MARKS:marks})
            def check(value): verify.verify_record(row,'\n'.join(json.dumps(x) for x in value.values()),error,None)
            check(wanted)
            bad=copy.deepcopy(wanted); bad[verify.epoch.PROGRESS]['selection_count']+=1
            with self.assertRaises(ValueError): check(bad)
            bad=copy.deepcopy(wanted); bad[overlay.SCHEMA]={'schema':overlay.SCHEMA}
            with self.assertRaises(ValueError): check(bad)


class RecordedFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw = os.environ.get('LEO_PAIRED_EPOCH_EVIDENCE')
        if not raw: raise unittest.SkipTest('retained fixtures require LEO_PAIRED_EPOCH_EVIDENCE')
        cls.root=Path(raw).resolve(strict=True)

    def fixture(self,mode):
        old=next(r for r in verify.epoch.inventory() if r['profile']=='release' and r['code']==0 and
                 r['args']==[mode,'8','0110','256'])
        folder=self.root/'checks'; label=old['label']
        output,error,_=verify.split_scope((folder/(label+'.stdout')).read_text(),(folder/(label+'.stderr')).read_text(),0)
        img=verify.epoch.elf(self.root/'build/release'/old['binary'],False)
        verify.epoch.verify_record(old,output,error,img)
        records=verify.epoch.records(output)
        driver=records.pop(verify.epoch.DRIVER); driver['schema']=overlay.SCHEMA
        records[overlay.SCHEMA]=driver
        return records,img

    def test_plain_record_shape(self):
        records,img=self.fixture('--exercise')
        for key in (verify.epoch.WITNESS,verify.epoch.MARKS): records.pop(key)
        records[overlay.SCHEMA]['clock_source']='steady'
        row=next(r for r in verify.inventory() if r['profile']=='release' and r['variant']=='steady' and
                 r['args']==['--exercise','8','0110','256'])
        def check(value): verify.verify_record(row,'\n'.join(json.dumps(x) for x in value.values()),'',img)
        check(records)
        bad=copy.deepcopy(records); bad[verify.epoch.WITNESS]={'schema':verify.epoch.WITNESS}
        with self.assertRaises(ValueError): check(bad)

    def test_fake_measure_shape_and_later_epochs(self):
        records,img=self.fixture('--clock-exercise')
        records[overlay.SCHEMA].update(clock_source='steady',timed=True)
        row=next(r for r in verify.inventory() if r['profile']=='release' and r['variant']=='fake-steady' and
                 r['args']==['--measure','8','0110','256'] and r['code']==0)
        def check(value): verify.verify_record(row,'\n'.join(json.dumps(x) for x in value.values()),'',img)
        check(records)
        bad=copy.deepcopy(records); bad[verify.epoch.ACCOUNT]['epochs'][2]['sample_begin']=0
        with self.assertRaises(ValueError): check(bad)
        bad=copy.deepcopy(records); bad[verify.epoch.META]['snapshots'][5]['epoch']=1
        with self.assertRaises(ValueError): check(bad)
        bad=copy.deepcopy(records); bad[overlay.SCHEMA]['samples'].pop()
        with self.assertRaises(ValueError): check(bad)


class NewRecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw=os.environ.get('LEO_PAIRED_EPOCH_TIMING_EVIDENCE')
        if not raw: raise unittest.SkipTest('new native fixtures require LEO_PAIRED_EPOCH_TIMING_EVIDENCE')
        cls.root=Path(raw).resolve(strict=True)

    def test_actual_variant_records_and_late_mutations(self):
        for profile in verify.ARCHIVES:
            for variant in image.VARIANTS:
                row=next(r for r in verify.inventory() if r['profile']==profile and r['variant']==variant and
                         r['code']==0 and r['args'][1]=='8' and r['args'][3]=='256' and r['args'][0]!='--check')
                label=row['label']; folder=self.root/'checks'
                output,error,_=verify.split_scope((folder/(label+'.stdout')).read_text(),
                                                 (folder/(label+'.stderr')).read_text(),0)
                img=image.elf(self.root/'build'/verify.executable(profile,variant),profile=='native',variant)
                def check(records): verify.verify_record(row,'\n'.join(json.dumps(r) for r in records.values()),error,img)
                records=verify.epoch.records(output); check(records)
                bad=copy.deepcopy(records); bad[verify.epoch.ACCOUNT]['epochs'][2]['calls']+=1
                with self.assertRaises(ValueError): check(bad)
                bad=copy.deepcopy(records); bad[verify.epoch.META]['snapshots'][5]['epoch']=0
                with self.assertRaises(ValueError): check(bad)
                bad=copy.deepcopy(records); bad[verify.epoch.PROGRESS]['selection_count']=104
                with self.assertRaises(ValueError): check(bad)
                if variant in ('synthetic','fake-steady'):
                    bad=copy.deepcopy(records); bad[overlay.SCHEMA]['samples'][251][0]+=1
                    with self.assertRaises(ValueError): check(bad)
                    bad=copy.deepcopy(records); bad[verify.epoch.CLOCK]['clock_calls']=168
                    with self.assertRaises(ValueError): check(bad)


if __name__ == '__main__': unittest.main()
