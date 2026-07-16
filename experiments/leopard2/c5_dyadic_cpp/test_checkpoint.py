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
    row = {
        "field_bits": 8 if q <= 129 else 16,
        "q": q,
        "parent": 1 << (q - 1).bit_length(),
        "valid_bytes": valid_bytes,
        "processed_bytes": (valid_bytes + 63) // 64 * 64,
        "batch": batch,
        "reuse": reuse,
        "samples": 7 if valid_bytes >= 65536 else 11,
        "setup_median_us": 1.0,
        "setup_mad_us": 0.1,
        "candidate_median_us": 10.0,
        "candidate_mad_us": 0.1,
        "padded_median_us": 10.0,
        "padded_mad_us": 0.1,
        "candidate_scratch_bytes_per_stripe": 64,
        "padded_scratch_bytes_per_stripe": 64,
        "resident_batch_bytes": 1024,
        "candidate_traffic_bytes_per_execution": 900,
        "padded_traffic_bytes_per_execution": 1000,
    }
    ratio, gain = MODULE.recomputed_benchmark_metrics(row)
    row["padded_over_candidate"] = round(ratio, 6)
    row["credible_gain_percent"] = round(gain, 6)
    return row


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


def correctness_case(field_bits: int, q: int, valid_bytes: int) -> dict:
    parent = MODULE.ceil_power_of_two(q)
    processed_bytes = (valid_bytes + 63) // 64 * 64
    return {
        "field_bits": field_bits,
        "q": q,
        "parent": parent,
        "valid_bytes": valid_bytes,
        "processed_bytes": processed_bytes,
        "blocks": q.bit_count(),
        "jobs": 1,
        "factor_checks": MODULE.expected_factor_checks(q),
        "compared_bytes": parent * processed_bytes,
        "direct_symbol_checks": MODULE.expected_direct_checks(
            field_bits, q, valid_bytes),
        "candidate_scratch_bytes": 0,
        "padded_scratch_bytes": 0,
        "candidate_traffic_bytes": 1,
        "padded_traffic_bytes": 1,
    }


def full_correctness_cases() -> list[dict]:
    return [
        correctness_case(field_bits, q, valid_bytes)
        for field_bits, q, valid_bytes in
        sorted(MODULE.expected_correctness_cells())
    ]


def authoritative_payload() -> dict:
    return {
        "requested_backend": "avx2",
        "runtime_backend": "avx2",
        "runtime_environment": {
            "process_affinity_cpus": [15],
            "omp_num_threads_env": "1",
            "openmp_max_threads": 1,
        },
    }


def set_timing(
    row: dict,
    candidate: float,
    candidate_mad: float,
    padded: float,
    padded_mad: float,
) -> None:
    row["candidate_median_us"] = candidate
    row["candidate_mad_us"] = candidate_mad
    row["padded_median_us"] = padded
    row["padded_mad_us"] = padded_mad
    ratio, gain = MODULE.recomputed_benchmark_metrics(row)
    row["padded_over_candidate"] = round(ratio, 6)
    row["credible_gain_percent"] = round(gain, 6)


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

    def test_benchmark_matrix_rejects_forged_ratio(self) -> None:
        rows = full_benchmark_rows()
        rows[0]["padded_over_candidate"] += 0.01
        with self.assertRaises(MODULE.CheckpointError):
            MODULE.validate_benchmarks(pathlib.Path("synthetic"), {
                "benchmarks": rows,
            })

    def test_benchmark_matrix_rejects_forged_credible_gain(self) -> None:
        rows = full_benchmark_rows()
        rows[0]["credible_gain_percent"] = 42.0
        with self.assertRaises(MODULE.CheckpointError):
            MODULE.validate_benchmarks(pathlib.Path("synthetic"), {
                "benchmarks": rows,
            })

    def test_benchmark_matrix_rejects_negative_mad(self) -> None:
        rows = full_benchmark_rows()
        rows[0]["candidate_mad_us"] = -0.001
        with self.assertRaises(MODULE.CheckpointError):
            MODULE.validate_benchmarks(pathlib.Path("synthetic"), {
                "benchmarks": rows,
            })

    def test_benchmark_matrix_rejects_duplicate_cell(self) -> None:
        rows = full_benchmark_rows()
        rows[-1] = copy.deepcopy(rows[0])
        with self.assertRaises(MODULE.CheckpointError):
            MODULE.validate_benchmarks(pathlib.Path("synthetic"), {
                "benchmarks": rows,
            })

    def test_benchmark_matrix_rejects_missing_cell(self) -> None:
        rows = full_benchmark_rows()
        rows.pop()
        with self.assertRaises(MODULE.CheckpointError):
            MODULE.validate_benchmarks(pathlib.Path("synthetic"), {
                "benchmarks": rows,
            })

    def test_correctness_matrix_is_exact_and_unique(self) -> None:
        payload = {"correctness_cases": full_correctness_cases()}
        MODULE.validate_correctness_cases(pathlib.Path("synthetic"), payload)
        duplicate = copy.deepcopy(payload)
        duplicate["correctness_cases"][-1] = copy.deepcopy(
            duplicate["correctness_cases"][0])
        with self.assertRaises(MODULE.CheckpointError):
            MODULE.validate_correctness_cases(
                pathlib.Path("synthetic"), duplicate)

    def test_correctness_matrix_rejects_redistributed_factor_checks(self) -> None:
        payload = {"correctness_cases": full_correctness_cases()}
        payload["correctness_cases"][0]["factor_checks"] += 1
        payload["correctness_cases"][1]["factor_checks"] -= 1
        with self.assertRaises(MODULE.CheckpointError):
            MODULE.validate_correctness_cases(
                pathlib.Path("synthetic"), payload)

    def test_authoritative_runtime_metadata_is_fail_closed(self) -> None:
        MODULE.validate_authoritative_runtime(
            pathlib.Path("synthetic"), authoritative_payload())
        mutations = (
            ("requested", lambda value: value.update(
                {"requested_backend": "scalar"})),
            ("runtime", lambda value: value.update(
                {"runtime_backend": "scalar"})),
            ("affinity", lambda value: value["runtime_environment"].update(
                {"process_affinity_cpus": [15, 31]})),
            ("omp-env", lambda value: value["runtime_environment"].update(
                {"omp_num_threads_env": "2"})),
            ("omp-runtime", lambda value: value["runtime_environment"].update(
                {"openmp_max_threads": 2})),
        )
        for label, mutate in mutations:
            with self.subTest(label=label):
                payload = authoritative_payload()
                mutate(payload)
                with self.assertRaises(MODULE.CheckpointError):
                    MODULE.validate_authoritative_runtime(
                        pathlib.Path("synthetic"), payload)

    def test_forced_backend_must_match_runtime(self) -> None:
        for label, runtime in (
            ("auto", "scalar"),
            ("scalar", "scalar"),
            ("ssse3", "ssse3"),
            ("avx2", "avx2"),
        ):
            with self.subTest(label=label):
                MODULE.validate_forced_backend_runtime(
                    pathlib.Path("synthetic"), label,
                    {"runtime_backend": runtime})
                with self.assertRaises(MODULE.CheckpointError):
                    MODULE.validate_forced_backend_runtime(
                        pathlib.Path("synthetic"), label,
                        {"runtime_backend": "neon"})

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
            set_timing(row, 10.0, 0.1, 12.0, 0.1)
        summary = MODULE.summarize_benchmarks(rows)
        self.assertEqual(summary["qualified_q"], [33])
        self.assertFalse(summary["production_promotion"])
        self.assertTrue(summary["experimental_followup"])

    def test_neighbor_regression_blocks_followup(self) -> None:
        rows = full_benchmark_rows()
        selected = [row for row in rows if row["q"] == 33]
        for row in selected[:3]:
            set_timing(row, 10.0, 0.1, 12.0, 0.1)
        set_timing(selected[-1], 10.0, 0.1, 9.9, 0.1)
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
            "correctness_cases": [{
                "q": 33,
                "jobs": 2,
                "candidate_scratch_bytes": 64,
                "candidate_traffic_bytes": 1024,
            }],
        }
        for key in (
            "q", "jobs", "candidate_scratch_bytes",
            "candidate_traffic_bytes",
        ):
            with self.subTest(key=key):
                changed = copy.deepcopy(payload)
                changed["correctness_cases"][0][key] += 1
                self.assertNotEqual(
                    MODULE.correctness_projection(payload),
                    MODULE.correctness_projection(changed),
                )


if __name__ == "__main__":
    unittest.main()
