"""Pure schedule/replay checks; no codec, clock or subprocess is executed."""
import copy
import unittest
from verify_gf16_callback_probe import CELLS, model, validate_counts


def record(cell):
    k,r,public_bytes,kind,tile,passes = CELLS[cell]
    return {"schema":"gf16-callback-counts/v1","timed":False,"calls":sum(model(cell).values()),
        "passes":[{"kind":kind,"k":k,"r":r,"requested":r,"side":1 << (r-1).bit_length(),
                   "sparse_blocks":0,"bytes":tile,"source_policy":public_bytes} for _ in range(passes)],
        "buckets":[dict(op=op,distance=distance,zero_mask=mask,prefer_fused=hint,bytes=size,calls=n)
                   for (op,distance,mask,hint,size),n in model(cell).items()]}


class CallbackModelTests(unittest.TestCase):
    def test_exact_target_stage_arithmetic(self):
        result = validate_counts(record(0),0)
        self.assertEqual(sum(model(0).values()),4132)
        self.assertEqual(result["ifft2"]["calls"],2000+768)
        self.assertEqual(result["ifft2_xor"]["calls"],768)
        self.assertEqual(result["ifft4_range"],{"calls":160,"lane_groups":1144})
        self.assertEqual(result["fft4_range"],{"calls":36,"lane_groups":360})
        self.assertEqual(result["fft2"]["calls"]+result["xor"]["calls"],400)
        self.assertEqual(4*(1144+360)+2768+768+400,9952)

    def test_all_fixed_shapes_and_odd_layer_side(self):
        for cell,total in enumerate((4132,4132,2066,2066,4132,3749)):
            self.assertEqual(record(cell)["calls"],total)
            validate_counts(record(cell),cell)
        result = validate_counts(record(5),5)
        self.assertEqual(result["ifft4_out"]["calls"],4096//4)
        self.assertEqual(result["ifft2_xor"]["calls"],7*256)
        self.assertEqual(result["ifft4_range"],{"calls":8*42,"lane_groups":8*3*128})
        self.assertEqual(result["fft4_range"],{"calls":85,"lane_groups":512})

    def test_every_pass_field_is_bound(self):
        for cell in range(6):
            for key in record(cell)["passes"][0]:
                changed = record(cell)
                changed["passes"][0][key] += 1
                with self.assertRaises(ValueError): validate_counts(changed,cell)

    def test_bucket_corruption_and_metadata_reject(self):
        for cell in range(6):
            base = record(cell)
            for key,value in (("schema","other"),("timed",True),("calls",True),("calls",0)):
                with self.assertRaises(ValueError): validate_counts(dict(base,**{key:value}),cell)
            for key,value in (("op","other"),("distance",0),("distance",True),
                              ("zero_mask",7),("prefer_fused",True),("bytes",3),("calls",False)):
                changed = copy.deepcopy(base)
                changed["buckets"][0][key] = value
                with self.assertRaises(ValueError): validate_counts(changed,cell)
            changed = copy.deepcopy(base)
            changed["buckets"].append(dict(changed["buckets"][0]))
            with self.assertRaises(ValueError): validate_counts(changed,cell)
            changed = copy.deepcopy(base)
            changed["buckets"].pop()
            with self.assertRaises(ValueError): validate_counts(changed,cell)


if __name__ == "__main__":
    unittest.main()
