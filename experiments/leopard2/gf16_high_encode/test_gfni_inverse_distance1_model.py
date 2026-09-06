import copy
import unittest
from test_gf16_callback_model import record
from verify_gfni_inverse_distance1 import validate


def candidate():
    result = record(0)
    for b in result["buckets"]:
        if b["op"] == "ifft2": b["calls"] -= 2000
    result["buckets"].append(dict(op="ifft4_range",distance=1,zero_mask=0,
        prefer_fused=False,bytes=32768,calls=500))
    result["calls"] -= 1500
    return result


class InverseDistanceModelTests(unittest.TestCase):
    def test_exact_substitution(self):
        validate(candidate(),0,True)
        self.assertEqual(candidate()["calls"],2632)
        for cell in range(6): validate(record(cell),cell,False)
        for cell in range(1,6): validate(record(cell),cell,True)

    def test_no_implicit_activation(self):
        with self.assertRaises(ValueError): validate(candidate(),0,False)
        with self.assertRaises(ValueError): validate(record(0),0,True)

    def test_exact_new_bucket(self):
        for key,value in (("op","ifft2"),("distance",4),("distance",True),
                          ("zero_mask",1),("prefer_fused",True),("bytes",32770),
                          ("calls",499),("calls",True)):
            changed = candidate(); changed["buckets"][-1][key] = value
            with self.assertRaises(ValueError): validate(changed,0,True)
        changed = candidate(); changed["buckets"].append(copy.deepcopy(changed["buckets"][-1]))
        with self.assertRaises(ValueError): validate(changed,0,True)

    def test_remaining_callbacks_unchanged(self):
        for index in range(len(candidate()["buckets"])-1):
            changed = candidate(); changed["buckets"][index]["calls"] += 1
            with self.assertRaises(ValueError): validate(changed,0,True)


if __name__ == "__main__": unittest.main()
