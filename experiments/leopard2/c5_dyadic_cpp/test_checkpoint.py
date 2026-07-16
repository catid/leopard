#!/usr/bin/env python3
"""Fail-closed unit tests for the C5 checkpoint merger."""

from __future__ import annotations

import copy
import importlib.util
import pathlib
import tempfile
import unittest


HERE = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("c5_checkpoint", HERE / "checkpoint.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def benchmark_row(q: int, valid_bytes: int, batch: int, reuse: int) -> dict:
    return {
        "field_bits": 8 if q <= 129 else 16,
        "q": q,
        "parent": 1 << (q - 1).bit_length(),
        "valid_bytes": valid_bytes,
        "processed_bytes": (valid_bytes + 63) // 64 * 64,
        "batch": batch,
        "reuse": reuse,
        "samples": 11,
        "setup_median_us": 1.0,
        "setup_mad_us": 0.1,
        "candidate_median_us": 10.0,
        "candidate_mad_us": 0.1,
        "padded_median_us": 10.0,
        "padded_mad_us": 0.1,
        "padded_over_candidate": 1.0,
        "credible_gain_percent": -1.99,
        "candidate_scratch_bytes_per_stripe": 64,
        "padded_scratch_bytes_per_stripe": 64,
        "resident_batch_bytes": 1024,
        "candidate_traffic_bytes_per_execution": 900,
        "padded_traffic_bytes_per_execution": 1000,
    }


def full_benchmark_rows() -> list[dict]:
    rows: list[dict] = []
    for q in MODULE.BENCHMARK_Q:
        for valid_bytes in (64, 1024):
            for batch in (1, 8):
                for reuse in (1, 8):
                    rows.append(benchmark_row(q, valid_bytes, batch, reuse))
        for reuse in (1, 8):
            rows.append(benchmark_row(q, 65536, 1, reuse))
    for q in (33, 65, 129, 257):
        rows.append(benchmark_row(q, 1024 * 1024, 1, 1))
    assert len(rows) == 54
    return rows


class CheckpointTests(unittest.TestCase):
    def test_parse_binding_rejects_duplicates(self) -> None:
        with self.assertRaises(MODULE.CheckpointError):
            MODULE.parse_binding(["auto=a", "auto=b"], "--backend")

    def test_benchmark_matrix_accepts_declared_shape(self) -> None:
        MODULE.validate_benchmarks(pathlib.Path("synthetic"), {
            "benchmarks": full_benchmark_rows(),
        })

    def test_benchmark_matrix_rejects_even_samples(self) -> None:
        rows = full_benchmark_rows()
        rows[0]["samples"] = 10
        with self.assertRaises(MODULE.CheckpointError):
            MODULE.validate_benchmarks(pathlib.Path("synthetic"), {
                "benchmarks": rows,
            })

    def test_noisy_equal_candidate_is_rejected(self) -> None:
        summary = MODULE.summarize_benchmarks(full_benchmark_rows())
        self.assertFalse(summary["production_promotion"])
        self.assertFalse(summary["experimental_followup"])
        self.assertEqual(summary["qualified_q"], [])
        self.assertTrue(summary["disposition"].startswith("reject:"))

    def test_three_credible_cells_across_two_sizes_qualify_followup(self) -> None:
        rows = full_benchmark_rows()
        selected = [row for row in rows if row["q"] == 33]
        winners = [
            next(row for row in selected if row["valid_bytes"] == 64),
            next(row for row in selected if row["valid_bytes"] == 1024),
            next(row for row in selected
                 if row["valid_bytes"] == 1024 and row["batch"] == 8),
        ]
        for row in winners:
            row["credible_gain_percent"] = 12.0
        summary = MODULE.summarize_benchmarks(rows)
        self.assertEqual(summary["qualified_q"], [33])
        self.assertFalse(summary["production_promotion"])
        self.assertTrue(summary["experimental_followup"])

    def test_neighbor_regression_blocks_followup(self) -> None:
        rows = full_benchmark_rows()
        selected = [row for row in rows if row["q"] == 33]
        for row in selected[:3]:
            row["credible_gain_percent"] = 12.0
        selected[-1]["credible_gain_percent"] = -2.01
        summary = MODULE.summarize_benchmarks(rows)
        self.assertNotIn(33, summary["qualified_q"])

    def test_sha256_is_content_bound(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "input"
            path.write_bytes(b"c5\n")
            first = MODULE.sha256(path)
            path.write_bytes(b"C5\n")
            second = MODULE.sha256(path)
            self.assertNotEqual(first, second)

    def test_correctness_projection_detects_one_byte_equivalent_mutation(self) -> None:
        payload = {
            "correctness_digest_fnv1a64": "0x1",
            "direct_symbol_checks": 1,
            "factor_checks": 2,
            "correctness_cases": [{"q": 33}],
        }
        changed = copy.deepcopy(payload)
        changed["correctness_cases"][0]["q"] = 34
        self.assertNotEqual(
            MODULE.correctness_projection(payload),
            MODULE.correctness_projection(changed),
        )


if __name__ == "__main__":
    unittest.main()
