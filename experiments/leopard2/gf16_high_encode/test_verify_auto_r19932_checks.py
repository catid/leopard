import copy
import unittest
from verify_auto_r19932_checks import arguments, guard_expected, inventory, public_expected, same, scope


class QualificationTests(unittest.TestCase):
    def test_inventory(self):
        rows = inventory()
        self.assertEqual(len(rows),259)
        self.assertEqual(len(set(n for n,_ in rows)),259)
        self.assertEqual([sum(c==n for _,c in rows) for n in (0,1,86)],[223,31,5])
        for name,_ in rows: self.assertTrue(arguments(name))

    def test_target_and_unchanged_routes(self):
        for cell in range(9):
            old,new = [public_expected('release',m,cell) for m in (0,1)]
            if cell<2:
                self.assertEqual(old['execution_route'],'avx2')
                self.assertEqual(new['execution_route'],'gfni')
            else:
                old['boundary_mode']=1
                self.assertEqual(old,new)
            for mode in range(2):
                value = public_expected('sanitize',mode,cell,True)
                self.assertEqual(value['encode_calls'],26)
                self.assertEqual(value['samples_ns'],[])

    def test_reject_record_mutations(self):
        expected = public_expected('release',1,0)
        for key,value in [('cell',False),('k',999),('r',200),('bytes',65536),
                          ('boundary_mode',True),('field',1),('api','leo_encode'),
                          ('execution_route','avx2'),('scratch_bytes',16777216),
                          ('untimed_route_calls',0),('input_hash','changed'),
                          ('output_hash','changed'),('samples_ns',[1]),('encode_calls',26),
                          ('codec_commit',public_expected('sanitize',1,0)['codec_commit'])]:
            changed = copy.deepcopy(expected); changed[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError): same(changed,expected)

    def test_native_mode_and_cell_refusals(self):
        for profile,mode,cell in [('native',1,0),('release',2,0),('release',True,0),
                                  ('other',0,0),('release',0,9),('release',0,False)]:
            with self.assertRaises(ValueError): public_expected(profile,mode,cell)

    def test_exact_cli(self):
        self.assertEqual(arguments('release-old-fault-1-kat')[1:],['--fault','1','kat'])
        self.assertEqual(arguments('release-new-fault-oom')[1:],['--fault','oom'])
        self.assertEqual(arguments('release-new-guard-1-3')[1:],['--guards','3','1'])
        self.assertEqual(arguments('release-1-0-exercise')[1:],['--exercise','0','1'])
        self.assertEqual(arguments('native-0-clock-guard')[1:],['--measure','0','0'])

    def test_guards(self):
        for old in (True,False):
            for cell in range(8):
                value = guard_expected(old,cell)
                self.assertEqual(value['subset_masks'],6)
                self.assertFalse(value['timed'])
        self.assertEqual(guard_expected(False,3)['scratch_bytes'],16839744)
        self.assertEqual(guard_expected(True,3)['bytes'],65536)
        self.assertEqual(guard_expected(False,3)['bytes'],32766)

    def test_resource_pressure_cannot_be_hidden(self):
        text = '\tExit status: 0\nmemory.peak\n256\nmemory.max\n256\nmemory.events\n'
        text += 'low 0\nhigh 0\nmax 612\noom 0\noom_kill 0\noom_group_kill 0\n'
        text += 'memory.swap.current\n0\nmemory.swap.max\n0\n'
        result = scope(text,256,612)
        self.assertEqual(result['max_events'],612)
        with self.assertRaises(ValueError): scope(text,256)
        for before,after in [('oom 0','oom 1'),('oom_kill 0','oom_kill 1'),
                             ('status: 0','status: 1'),('peak\n256','peak\n257'),
                             ('swap.current\n0','swap.current\n1')]:
            with self.assertRaises(ValueError): scope(text.replace(before,after),256,612)


if __name__=='__main__': unittest.main()
