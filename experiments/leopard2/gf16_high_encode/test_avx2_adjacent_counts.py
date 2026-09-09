import copy
import unittest
from avx2_adjacent_counts import SHAPES, validate
from verify_gf16_callback_probe import model_shape


def record(cell):
    if cell == 7:
        return dict(schema='gf16-callback-counts/v1',timed=False,calls=0,passes=[],buckets=[])
    k,r,public_bytes,kind,tile,passes = SHAPES[cell]
    model = model_shape(SHAPES[cell])
    return dict(schema='gf16-callback-counts/v1',timed=False,calls=sum(model.values()),
        passes=[dict(kind=kind,k=k,r=r,requested=r,side=1 << (r-1).bit_length(),sparse_blocks=0,
                     bytes=tile,source_policy=public_bytes) for _ in range(passes)],
        buckets=[dict(op=op,distance=distance,zero_mask=zero,prefer_fused=fused,bytes=size,calls=calls)
                 for (op,distance,zero,fused,size),calls in model.items()])


class AdjacentCountTests(unittest.TestCase):
    def test_expected_target_and_neighbor_pairs(self):
        for cell, pairs in enumerate(((665,384),(1330,768),(1330,768),(1793,1792),(665,384))):
            actual = validate(record(cell),cell)
            self.assertEqual((actual['forward_pairs'],actual['accumulating_pairs']),pairs)
            blocks = 64 if cell == 3 else 512
            self.assertEqual(actual['forward_blocks'],pairs[0]*blocks)
            self.assertEqual(actual['accumulating_blocks'],pairs[1]*blocks)

    def test_unchanged_backend_and_field_neighbors(self):
        for cell in (5,6,7):
            actual = validate(record(cell),cell)
            self.assertFalse(actual['affected'])
            self.assertEqual(actual['forward_blocks'] + actual['accumulating_blocks'],0)

    def test_native_bucket_corruption_and_types(self):
        for cell in range(7):
            for key,value in (('distance',0),('calls',True),('zero_mask',7),('prefer_fused',True),('op','invented')):
                bad = copy.deepcopy(record(cell)); bad['buckets'][0][key] = value
                with self.assertRaises(ValueError): validate(bad,cell)
            bad = record(cell); bad['passes'][0]['kind'] = 5
            with self.assertRaises(ValueError): validate(bad,cell)
        for key,value in (('timed',0),('calls',False),('calls',1),('buckets',[{}])):
            bad = record(7); bad[key] = value
            with self.assertRaises(ValueError): validate(bad,7)

    def test_bad_cell_and_dropped_bucket(self):
        for cell in (-1,8,True,'0'):
            with self.assertRaises(ValueError): validate(record(0),cell)
        bad = record(0); bad['buckets'].pop()
        with self.assertRaises(ValueError): validate(bad,0)


if __name__ == '__main__':
    unittest.main()
