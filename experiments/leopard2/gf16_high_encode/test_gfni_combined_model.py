"""Independent fixtures and malformed-evidence checks for the four-mode replay."""
import copy
import unittest
from test_gf16_callback_model import record
from verify_gfni_combined import expected_hook, validate_hook, validate


def candidate(cell, mode):
    result = record(cell)
    if cell != 0:
        return result
    removed_pairs = (2000 if mode & 1 else 0) + (768 if mode & 2 else 0)
    for bucket in result["buckets"]:
        if bucket["op"] == "ifft2":
            bucket["calls"] -= removed_pairs
    result["buckets"] = [b for b in result["buckets"] if b["calls"] and
                         not (mode & 2 and b["op"] == "ifft2_xor")]
    if mode & 1:
        result["buckets"].append(dict(op="ifft4_range", distance=1, zero_mask=0,
            prefer_fused=False, bytes=32768, calls=500))
    result["calls"] -= (1500 if mode & 1 else 0) + (1536 if mode & 2 else 0)
    return result


def kernel(cell, mode):
    selected = cell == 0 and bool(mode & 2)
    return dict(schema="gfni-terminal-kernel-counts/v1", timed=False,
        calls=6 if selected else 0, lane_groups=384 if selected else 0)


class CombinedModelTests(unittest.TestCase):
    def test_all_modes_and_neighbors(self):
        for cell in range(6):
            for mode in range(4):
                validate(candidate(cell, mode), kernel(cell, mode), cell, mode)
                validate_hook(expected_hook(cell, mode), cell, mode)
        self.assertEqual([candidate(0, mode)["calls"] for mode in range(4)], [4132, 2632, 2596, 1096])

    def test_combined_has_no_pair_buckets(self):
        self.assertFalse(any(b["op"] in ("ifft2", "ifft2_xor") for b in candidate(0, 3)["buckets"]))
        validate(candidate(0, 3), kernel(0, 3), 0, 3)

    def test_all_wrong_mode_combinations_refused(self):
        for actual in range(4):
            for requested in range(4):
                if actual != requested:
                    with self.assertRaises(ValueError):
                        validate(candidate(0, actual), kernel(0, actual), 0, requested)

    def test_scope_types_and_bounds(self):
        for mode in (-1, 4, True, 1.0, "1", None):
            with self.assertRaises(ValueError): validate(record(0), kernel(0, 0), 0, mode)
        for cell in (-1, 6, False, 0.0, "0", None):
            with self.assertRaises(ValueError): validate(record(0), kernel(0, 0), cell, 0)

    def test_bucket_corruption_and_boolean_types(self):
        for mode in range(4):
            original = candidate(0, mode)
            for index in range(len(original["buckets"])):
                for key, value in (("op", "copy"), ("distance", True), ("zero_mask", 7),
                    ("prefer_fused", 0), ("bytes", 2), ("calls", True), ("calls", 0), ("extra", 1)):
                    changed = copy.deepcopy(original)
                    changed["buckets"][index][key] = value
                    with self.assertRaises(ValueError): validate(changed, kernel(0, mode), 0, mode)

    def test_duplicate_missing_and_zero_bucket_refused(self):
        for mode in range(4):
            original = candidate(0, mode)
            duplicate = copy.deepcopy(original)
            duplicate["buckets"].append(copy.deepcopy(duplicate["buckets"][0]))
            missing = copy.deepcopy(original); missing["buckets"].pop()
            zero = copy.deepcopy(original)
            zero["buckets"].append(dict(op="ifft2", distance=1, zero_mask=0,
                                       prefer_fused=False, bytes=32768, calls=0))
            for changed in (duplicate, missing, zero):
                with self.assertRaises(ValueError): validate(changed, kernel(0, mode), 0, mode)

    def test_callback_metadata_exact(self):
        for mode in range(4):
            original = candidate(0, mode)
            for key, value in (("schema", "other"), ("timed", 0), ("calls", True), ("extra", 1),
                                ("buckets", {}), ("passes", [])):
                changed = dict(original, **{key:value})
                with self.assertRaises(ValueError): validate(changed, kernel(0, mode), 0, mode)
            for key in original["passes"][0]:
                changed = copy.deepcopy(original)
                changed["passes"][0][key] += 1
                with self.assertRaises(ValueError): validate(changed, kernel(0, mode), 0, mode)

    def test_kernel_counts_exact(self):
        for mode in range(4):
            for key, value in (("schema", "other"), ("timed", 0), ("calls", True),
                                ("calls", 5), ("lane_groups", 383), ("extra", 0)):
                changed = kernel(0, mode); changed[key] = value
                with self.assertRaises(ValueError): validate(candidate(0, mode), changed, 0, mode)

    def test_hook_metadata_predicates_and_modes(self):
        for cell in range(6):
            for mode in range(4):
                original = expected_hook(cell, mode)
                for key in original:
                    changed = copy.deepcopy(original); changed[key] = None
                    with self.assertRaises(ValueError): validate_hook(changed, cell, mode)
                for key in original["records"][0]:
                    changed = copy.deepcopy(original); changed["records"][0][key] = None
                    with self.assertRaises(ValueError): validate_hook(changed, cell, mode)
                for key in ("mode", "calls", "matches", "first", "terminal"):
                    changed = copy.deepcopy(original); changed[key] = bool(changed[key])
                    with self.assertRaises(ValueError): validate_hook(changed, cell, mode)

    def test_order_irrelevant_inputs_not_mutated(self):
        for mode in range(4):
            changed = candidate(0, mode); changed["buckets"].reverse()
            before = copy.deepcopy(changed)
            validate(changed, kernel(0, mode), 0, mode)
            self.assertEqual(changed, before)


if __name__ == "__main__":
    unittest.main()
