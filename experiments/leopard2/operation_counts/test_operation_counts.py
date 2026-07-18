#!/usr/bin/env python3
"""Independent tests for tools/leopard2_operation_counts.py."""

from __future__ import annotations

import csv
import importlib.util
import io
import json
import pathlib
import subprocess
import sys
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[3]
TOOL = ROOT / "tools" / "leopard2_operation_counts.py"
SPEC = importlib.util.spec_from_file_location("leopard2_operation_counts", TOOL)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load operation-count tool")
COUNTS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = COUNTS
SPEC.loader.exec_module(COUNTS)


def independent_full_transform_count(size: int) -> int:
    operations = 0
    span = 2
    while span <= size:
        for _block in range(0, size, span):
            for _lane in range(span // 2):
                operations += 1
        span *= 2
    return operations


class OperationCountTests(unittest.TestCase):
    def test_full_schedules_against_explicit_layer_enumeration(self) -> None:
        for size in (1, 2, 4, 8, 16, 32, 64, 128, 256):
            expected = independent_full_transform_count(size)
            self.assertEqual(COUNTS.ifft_prefix_butterflies(size, size), expected)
            self.assertEqual(COUNTS.fft_prefix_butterflies(size, size), expected)
            self.assertEqual(
                COUNTS.fft_mask_butterflies(size, set(range(size))), expected
            )

    def test_mask_schedule_is_prefix_schedule_for_prefix_masks(self) -> None:
        for size in (2, 4, 8, 16, 32, 64):
            for prefix in range(size + 1):
                self.assertEqual(
                    COUNTS.fft_mask_butterflies(size, set(range(prefix))),
                    COUNTS.fft_prefix_butterflies(size, prefix),
                )

    def test_high_and_low_representative_cells(self) -> None:
        high = COUNTS.model_high_encode(240, 16, 16, 1024, set(range(16)))
        self.assertEqual(high.butterflies, 512)
        self.assertEqual(high.nontransform_xor_vectors, 224)
        low = COUNTS.model_low_encode(8, 248, 8, 1024, set(range(248)))
        self.assertEqual(low.butterflies, 384)
        self.assertEqual(low.details["active_parity_blocks"], 31)
        self.assertEqual(low.copies, 8 + 248)
        self.assertTrue(low.details["out_of_place_first_fft_layer"])

    def test_low_encode_source_guard_rejects_old_copy_loop(self) -> None:
        old_copy = """
            for (unsigned i = 0; i < p; ++i)
                memcpy(work[p + i], work[i], buffer_bytes);
        """
        for filename in ("LeopardFF8.cpp", "LeopardFF16.cpp"):
            source = (ROOT / filename).read_text(encoding="utf-8")
            COUNTS.verify_low_encode_no_copy_source(source, filename)
            with self.assertRaises(COUNTS.ModelError):
                COUNTS.verify_low_encode_no_copy_source(
                    source + old_copy, filename + " mutation"
                )

    def test_sparse_requested_parity_reduces_low_work(self) -> None:
        all_outputs = COUNTS.model_low_encode(8, 248, 8, 1024, set(range(248)))
        sparse = COUNTS.model_low_encode(8, 248, 8, 1024, {0, 247})
        self.assertLess(sparse.butterflies, all_outputs.butterflies)
        self.assertEqual(sparse.details["active_parity_blocks"], 2)
        self.assertEqual(sparse.copies, 8 + 2)

    def test_deterministic_received_subset(self) -> None:
        high_data, parity = COUNTS.deterministic_decode_coordinates(
            9, 7, "high", 8, {0, 3}
        )
        self.assertEqual(parity, {0, 1})
        self.assertEqual(len(high_data), 9)
        self.assertIn(8 + 1, high_data)
        self.assertNotIn(8 + 3, high_data)

    def test_direct_bound_and_dispatch_gates(self) -> None:
        model = COUNTS.model_direct_repair(16, 8, 16, 16, "low", set(range(4)))
        self.assertEqual(model.fixed_multiply_vectors, 64)
        self.assertEqual(model.nontransform_xor_vectors, 60)
        self.assertEqual(model.api_scratch_slots, 0)
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.model_direct_repair(17, 8, 32, 32, "low", {0})

    def test_json_and_csv_are_deterministic(self) -> None:
        command = [
            sys.executable, str(TOOL), "report",
            "--path", "generic_decode", "--profile", "high",
            "--k", "9", "--r", "7", "--field", "gf8",
            "--shard-bytes", "65", "--loss-mask", "0,3",
        ]
        first = subprocess.check_output(command, text=True)
        second = subprocess.check_output(command, text=True)
        self.assertEqual(first, second)
        document = json.loads(first)
        self.assertEqual(document["schema_version"], 1)
        self.assertEqual(
            document["metrics"]["transform_butterflies"]["classification"],
            "exact_schedule",
        )
        csv_text = subprocess.check_output(command + ["--format", "csv"], text=True)
        rows = list(csv.DictReader(io.StringIO(csv_text)))
        self.assertGreater(len(rows), 10)
        self.assertEqual([row["metric"] for row in rows],
                         sorted(row["metric"] for row in rows))

    def test_gf16_odd_length_is_rejected(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(TOOL), "report", "--path", "low_encode",
             "--k", "8", "--r", "8", "--field", "gf16",
             "--shard-bytes", "3"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("even shard byte count", completed.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
