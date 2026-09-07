import copy
from pathlib import Path
import tempfile
import unittest
from test_gf16_callback_model import record
from verify_gfni_terminal import validate, terminal_shape, resource_result


def kernel(selected=False):
    return dict(schema="gfni-terminal-kernel-counts/v1", timed=False,
                calls=6 if selected else 0, lane_groups=384 if selected else 0)


def candidate():
    result = record(0)
    for b in result["buckets"]:
        if b["op"] == "ifft2": b["calls"] -= 768
    result["buckets"] = [b for b in result["buckets"] if b["op"] != "ifft2_xor"]
    result["calls"] -= 1536
    return result


class TerminalModelTests(unittest.TestCase):
    def test_exact_substitution(self):
        validate(candidate(),kernel(True),0,True)
        self.assertEqual(candidate()["calls"],2596)
        self.assertEqual(terminal_shape(0),(6,384,100663296))
        self.assertEqual(next(b["calls"] for b in candidate()["buckets"] if b["op"]=="ifft2"),2000)
        for cell in range(6): validate(record(cell),kernel(),cell,False)
        for cell in range(1,6): validate(record(cell),kernel(),cell,True)

    def test_no_implicit_activation(self):
        for c,k,enabled in ((candidate(),kernel(True),False), (record(0),kernel(),True),
                            (candidate(),kernel(),True), (record(0),kernel(True),False)):
            with self.assertRaises(ValueError): validate(c,k,0,enabled)

    def test_exact_external_counts(self):
        for key,value in (("schema","other"),("timed",True),("calls",5),("calls",True),
                          ("lane_groups",383),("lane_groups",True),("extra",0)):
            changed = kernel(True); changed[key] = value
            with self.assertRaises(ValueError): validate(candidate(),changed,0,True)

    def test_all_remaining_callbacks_unchanged(self):
        for index in range(len(candidate()["buckets"])):
            for key,value in (("calls",True),("calls",1),("distance",True),
                              ("zero_mask",7),("prefer_fused",True),("bytes",2),("extra",1)):
                changed = candidate(); changed["buckets"][index][key] = value
                with self.assertRaises(ValueError): validate(changed,kernel(True),0,True)
        changed = candidate(); changed["calls"] = True
        with self.assertRaises(ValueError): validate(changed,kernel(True),0,True)

    def test_duplicate_and_restored_accumulating_bucket(self):
        for bucket in (candidate()["buckets"][0],
                       next(b for b in record(0)["buckets"] if b["op"]=="ifft2_xor")):
            changed=candidate(); changed["buckets"].append(copy.deepcopy(bucket))
            with self.assertRaises(ValueError): validate(changed,kernel(True),0,True)

    def test_resource_failures_are_not_passes(self):
        body = "\tExit status: 0\nmemory.peak\n1000\nmemory.max\n268435456\nmemory.events\nlow 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\nmemory.swap.current\n0\nmemory.swap.max\n0\n"
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"resource.log"
            path.write_text(body)
            self.assertEqual(resource_result(path),1000)
            with self.assertRaises(ValueError): resource_result(path,1)
            path.write_text(body.replace("Exit status: 0","Exit status: 1"))
            self.assertEqual(resource_result(path,1),1000)
            with self.assertRaises(ValueError): resource_result(path)
            for old,new in (("oom 0","oom 1"),("high 0","high 1"),("1000","268435457"),
                            ("memory.swap.current\n0","memory.swap.current\n1"),
                            ("memory.max\n268435456","memory.max\n536870912")):
                path.write_text(body.replace(old,new))
                with self.assertRaises(ValueError): resource_result(path)


if __name__ == "__main__": unittest.main()
