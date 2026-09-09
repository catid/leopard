import copy
import unittest
from verify_r199_boundary_checks import codec_identity, inventory, public_record, same


def record(profile='release',variant='auto'):
    native = profile=='native'
    return dict(schema='leopard-r199-boundary/v1',profile='native' if native else 'current',
        codec_commit=codec_identity(profile),cell=0,k=1000,r=199,bytes=32768,requested=variant,
        execution_route='native_l1' if native else 'avx2' if variant=='auto' else 'gfni',
        scratch_bytes=16777216 if native else 16808512,
        output_semantics='first_r_work_buffers' if native else 'separate_output_buffers',
        input_hash='8b78decd04e67d27',output_hash='00648d9dbf0f2b20',outer_guards=True,
        encode_calls=1,samples_ns=[])


class QualificationTests(unittest.TestCase):
    def test_inventory(self):
        rows = inventory()
        self.assertEqual(len(rows),82)
        self.assertEqual(len(set(k for k,_ in rows)),82)
        self.assertEqual([sum(code==n for _,code in rows) for n in (0,1,86)],[49,28,5])

    def test_all_routes_and_calls(self):
        for profile in ('native','release','sanitize'):
            for variant in (('native',) if profile=='native' else ('auto','gfni')):
                value = record(profile,variant)
                public_record(value,profile,variant)
                value['encode_calls'] = 26
                public_record(value,profile,variant,True)
                with self.assertRaises(ValueError): public_record(value,profile,variant)

    def test_reject_identity_and_timing_mutations(self):
        for key,value in (('cell',False),('k',999),('r',200),('bytes',65536),
                          ('execution_route','gfni'),('requested','gfni'),('profile','native'),
                          ('codec_commit',codec_identity('sanitize')),('scratch_bytes',16777216),
                          ('output_semantics','first_r_work_buffers'),('input_hash','changed'),
                          ('output_hash','changed'),('outer_guards',1),('encode_calls',True),
                          ('samples_ns',[1]),('schema','old')):
            bad = copy.deepcopy(record()); bad[key] = value
            with self.subTest(key=key),self.assertRaises(ValueError):
                public_record(bad,'release','auto')

    def test_strict_type_comparison(self):
        for a,b in ((1,True),(0,False),(1,1.0),([],{})):
            with self.assertRaises(ValueError): same(a,b)


if __name__ == '__main__': unittest.main()
