#!/usr/bin/env python3
"""Pure contract/mutation checks; no codec execution or timing."""
import copy
from pathlib import Path
import tempfile
import unittest

import check_encode as check


class ContractTests(unittest.TestCase):
    def fixture(self, index=0, main=True):
        name, k, r, size = check.CELLS[index]
        count = r if k == 1 else 1 if r == 1 else 2 * (1 << (r - 1).bit_length())
        return {
            "schema": "leopard-native-release-encode-check/v1",
            "codec_commit": check.BASELINE,
            "implementation": "leopard1-native" if main else "leopard2",
            "cell": index, "id": name, "k": k, "r": r, "bytes": size,
            "requested_backend": "native" if main else "avx2" if index == 6 else "auto",
            "context_backend": -1 if main else 3, "threads": 1,
            "input_bytes": k * size, "parity_bytes": r * size,
            "workspace_bytes": count * size if main else 0,
            "separate_output_bytes": 0 if main else r * size,
            "output_layout": "parity_is_first_r_work_rows" if main else
                             "separate_parity_and_scratch",
            "input_unchanged": True, "repeated_encode_equal": True,
            "public_encode_calls": 2, "input_hash": "a" * 16, "parity_hash": "b" * 16,
        }

    def test_all_fixed_records(self):
        for index in range(8):
            for main in (True, False):
                check.validate(self.fixture(index, main), index,
                               "main" if main else "current", check.BASELINE)

    def test_every_field_is_bound_or_validated(self):
        for main in (True, False):
            record = self.fixture(6, main)
            for key, value in record.items():
                bad = copy.deepcopy(record)
                bad[key] = "invalid" if isinstance(value, str) else None
                with self.subTest(main=main, field=key), self.assertRaises(ValueError):
                    check.validate(bad, 6, "main" if main else "current", check.BASELINE)

    def test_boolean_is_not_integer(self):
        record = self.fixture()
        record["k"] = True
        with self.assertRaises(ValueError):
            check.validate(record, 0, "main", check.BASELINE)

    def test_no_timing_fields(self):
        record = self.fixture()
        record["samples_ns"] = [100]
        with self.assertRaises(ValueError):
            check.validate(record, 0, "main", check.BASELINE)

    def test_native_geometry_and_explicit_backend(self):
        record = self.fixture(6)
        record["workspace_bytes"] += 64
        with self.assertRaises(ValueError):
            check.validate(record, 6, "main", check.BASELINE)
        record = self.fixture(6, False)
        record["context_backend"] = 6
        with self.assertRaises(ValueError):
            check.validate(record, 6, "current", check.BASELINE)

    def test_full_file_comparison(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = (Path(directory) / name for name in ("a", "b"))
            a.write_bytes(b"a" * 131073)
            b.write_bytes(a.read_bytes())
            check.equal_files(a, b, 131073)
            b.write_bytes(b"a" * 131072 + b"b")
            with self.assertRaises(ValueError):
                check.equal_files(a, b, 131073)
            b.write_bytes(b"a" * 131072)
            with self.assertRaises(ValueError):
                check.equal_files(a, b, 131073)


if __name__ == "__main__":
    unittest.main()
