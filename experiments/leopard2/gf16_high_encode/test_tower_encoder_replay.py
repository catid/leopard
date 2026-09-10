import unittest
from pathlib import Path
from replay_tower_encoder import audit_isa, resource
from tower_encoder_overlay import overlay


class QualificationTests(unittest.TestCase):
    @staticmethod
    def field_source():
        retained = Path(__file__).resolve().parent.parent/'build/source/LeopardFF16.original.cpp'
        return (retained if retained.exists() else Path(__file__).resolve().parents[3]/'LeopardFF16.cpp').read_text()

    def test_resource_success(self):
        self.assertEqual(resource(self.good_resources()), 1024)

    @staticmethod
    def good_resources():
        return ['memory.peak','1024','memory.max','268435456','memory.events','low 0','high 0','max 0',
                'oom 0','oom_kill 0','oom_group_kill 0','memory.swap.current','0','memory.swap.max','0']

    def test_resource_failure(self):
        for index in range(5,11):
            bad = self.good_resources(); bad[index] = bad[index].replace(' 0',' 1')
            with self.assertRaises(ValueError): resource(bad)
        for peak in ('0','-1','268435456'):
            bad = self.good_resources(); bad[1] = peak
            with self.assertRaises(ValueError): resource(bad)
        for index in (3,12,14):
            bad = self.good_resources(); bad[index] = '1'
            with self.assertRaises(ValueError): resource(bad)

    def test_baseline_and_avx2(self):
        self.assertEqual(audit_isa(' 0: 66 0f ef c0 pxor xmm0,xmm0',True),1)
        self.assertEqual(audit_isa(' 0: c5 fd ef c0 vpxor ymm0,ymm0,ymm0',False),1)
        with self.assertRaises(ValueError): audit_isa(' 0: c5 fd ef c0 vpxor ymm0,ymm0,ymm0',True)

    def test_wide_and_special_isa(self):
        for text in (' 0: 62 01 02 03 vpxor ymm0,ymm0,ymm0',
                     ' 0: c5 01 02 03 vpxor ymm16,ymm0,ymm0',
                     ' 0: c5 01 02 03 vpxor zmm0,zmm0,zmm0',
                     ' 0: c5 01 02 03 vpxor k1,k2,k3',
                     ' 0: c5 01 02 03 vgf2p8affineqb ymm0,ymm1,ymm2,0',
                     ' 0: c5 01 02 03 vpternlogd ymm0,ymm1,ymm2,0',
                     ' 0: c5 01 02 03 vpclmulqdq ymm0,ymm1,ymm2,0'):
            with self.assertRaises(ValueError): audit_isa(text,False)

    def test_empty_disassembly(self):
        with self.assertRaises(ValueError): audit_isa('missing function',False)

    def test_overlay_scope(self):
        source = self.field_source()
        changed = overlay(source)
        self.assertEqual(changed.count('tower_encoder::CopySource(ops, work[i], data[i], bytes);'),2)
        self.assertEqual(changed.count('tower_encoder::Select('),1)
        self.assertEqual(changed.count('tower_encoder::Finish(work, recovery_count, buffer_bytes);'),1)
        self.assertEqual(source[source.index('void ReedSolomonEncodeLow('):],
                         changed[changed.index('void ReedSolomonEncodeLow('):])

    def test_overlay_refuses_source_drift(self):
        source = self.field_source()
        for bad in (source+'\n',source.replace('memcpy(work[i], data[i], bytes);','bad();',1),''):
            with self.assertRaises(ValueError): overlay(bad)


if __name__ == '__main__':
    unittest.main()
