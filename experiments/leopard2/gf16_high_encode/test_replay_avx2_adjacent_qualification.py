import copy
import unittest
from replay_avx2_adjacent_qualification import parity_labels, public_identity, same


def record():
    return dict(cell=0,k=1000,r=200,bytes=32768,field=2,requested_backend=3,
                input_hash='input',output_hash='output',outer_guards=True,codec_commit='candidate',
                mode='off',traced=False,encode_calls=1,samples_ns=[],route_counts=[0,0,0,0])


class ReplayTests(unittest.TestCase):
    def test_parity_inventory_has_each_profile_and_cell_once(self):
        labels = parity_labels()
        self.assertEqual(len(labels),32)
        self.assertEqual(len(set(labels)),32)
        self.assertEqual(labels[0],'release-no-clock-0')
        self.assertEqual(labels[-1],'sanitize-observer-7')

    def test_public_binding_and_every_observed_field(self):
        base = record()
        public_identity('release-no-clock-0',base,base,'candidate')
        for key,value in (('cell',1),('cell',False),('k',999),('r',199),('bytes',65536),
                          ('field',1),('requested_backend',6),('input_hash','changed'),
                          ('output_hash','changed'),('outer_guards',False),('codec_commit','baseline'),
                          ('mode','on'),('traced',True),('encode_calls',True),('samples_ns',[1]),
                          ('route_counts',[1,0,0,0])):
            bad = copy.deepcopy(base); bad[key] = value
            with self.subTest(key=key,value=value), self.assertRaises(ValueError):
                public_identity('release-no-clock-0',bad,base,'candidate')

    def test_swapped_matching_public_and_native_cells_still_reject(self):
        wrong = record(); wrong['cell'] = 1
        with self.assertRaises(ValueError):
            public_identity('release-no-clock-0',wrong,wrong,'candidate')

    def test_type_aware_equality(self):
        for a,b in ((0,False),(1,True),(1,1.0),([],{}),({'value':1},{'value':True})):
            with self.assertRaises(ValueError): same(a,b)


if __name__ == '__main__':
    unittest.main()
