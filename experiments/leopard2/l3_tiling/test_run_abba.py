#!/usr/bin/env python3
"""Unit tests for the L3-tiling evidence contract (no timings are run)."""

from __future__ import annotations

import importlib.util
import json
import math
import os
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("run_abba.py")
SPEC = importlib.util.spec_from_file_location("leopard2_l3_tiling_runner", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load the L3-tiling evidence runner")
RUNNER = importlib.util.module_from_spec(SPEC)
import sys
sys.modules[SPEC.name] = RUNNER
SPEC.loader.exec_module(RUNNER)


class ContractTests(unittest.TestCase):
    def test_exact_matrix(self) -> None:
        observed = [
            (cell.profile, cell.k, cell.r, cell.shard_bytes, cell.losses)
            for cell in RUNNER.CELLS
        ]
        self.assertEqual(observed, [
            ("high", 300, 100, 393216, 9),
            ("high", 300, 100, 1048576, 9),
            ("low", 100, 193, 524288, 9),
            ("low", 100, 193, 1048576, 9),
            ("high", 1000, 200, 65536, 9),
            ("high", 2000, 500, 65536, 9),
            ("high", 2000, 500, 131072, 9),
            ("low", 200, 800, 65536, 9),
            ("high", 2000, 512, 65536, 512),
        ])
        self.assertEqual(len({cell.identifier for cell in RUNNER.CELLS}), 9)
        self.assertEqual(len({cell.seed for cell in RUNNER.CELLS}), 9)

    def test_abba_orders(self) -> None:
        self.assertEqual(RUNNER.ORDER_SIDES, (
            ("left", "right", "right", "left"),
            ("right", "left", "left", "right"),
            ("left", "right", "right", "left"),
        ))
        self.assertEqual(RUNNER.RETAINED_SAMPLES, 9)
        self.assertEqual(RUNNER.ROUNDS, 3)
        self.assertLessEqual(
            RUNNER.MAX_BUILD_OUTPUT,
            RUNNER.SUPPORT.MAX_COMMAND_STDOUT_BYTES)
        self.assertLessEqual(
            RUNNER.MAX_BUILD_OUTPUT,
            RUNNER.SUPPORT.MAX_COMMAND_STDERR_BYTES)

    def test_ratio_summary(self) -> None:
        result = RUNNER.ratio_summary([1.1, 1.1, 1.1])
        self.assertAlmostEqual(result["speedup_of_right_over_left"], 1.1)
        self.assertEqual(result["confidence_interval_95"], [1.1, 1.1])
        self.assertTrue(
            result["credible_right_improvement_at_least_5_percent"])
        self.assertFalse(
            result["point_estimate_right_regression_over_2_percent"])

    def test_parse_requires_exact_named_pair(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            # parse_lane resolves roots but source validation is intentionally
            # deferred until the campaign owns the canonical lock.
            commit_a = "1" * 40
            commit_b = "2" * 40
            options = RUNNER.parse_args([
                "--lane", "fixed", str(root), commit_a, "build/a/bench_leopard2",
                "--lane", "new", str(root), commit_b, "build/b/bench_leopard2",
                "--comparison", "fixed", "new",
                "--output", str(root / "out"),
                "--cpu", "4", "--sibling", "20",
                "--cells", "high-side-512-maxloss", "high-live-96m",
            ])
            self.assertEqual(
                [request.label for request in options.requests],
                ["fixed", "new"])
            self.assertEqual(
                [cell.identifier for cell in options.selected_cells],
                ["high-live-96m", "high-side-512-maxloss"])
            with self.assertRaises(RUNNER.EvidenceError):
                RUNNER.parse_args([
                    "--lane", "fixed", str(root), commit_a,
                    "build/a/bench_leopard2",
                    "--lane", "new", str(root), commit_b,
                    "build/b/bench_leopard2",
                    "--comparison", "fixed", "other",
                    "--output", str(root / "out"),
                    "--cpu", "4", "--sibling", "20",
                ])

    def test_freeze_and_tamper_detection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.write_bytes(b"#!/bin/sh\nexit 0\n")
            source.chmod(0o700)
            expected = RUNNER.sha256(source)
            destination = root / "artifacts" / "lane" / "bench_leopard2"
            record = RUNNER.freeze_binary(source, destination, expected)
            RUNNER.verify_file_identity(record["identity"], executable=True)
            destination.chmod(0o700)
            with destination.open("ab") as stream:
                stream.write(b"# tamper\n")
            with self.assertRaises(RUNNER.EvidenceError):
                RUNNER.verify_file_identity(
                    record["identity"], executable=True)
            # Restore directory permissions for TemporaryDirectory cleanup.
            destination.parent.chmod(0o700)

    def test_timing_result_retains_samples_and_path(self) -> None:
        cell = RUNNER.CELLS[0]
        metric = {
            "median_us_per_batch_call": 2.0,
            "samples_us_per_batch_call": [2.0] * RUNNER.RETAINED_SAMPLES,
        }
        result = {
            "schema": "leopard2-benchmark-v3",
            "parameters": {
                "K": cell.k, "R": cell.r,
                "shard_bytes": cell.shard_bytes,
                "loss_count": cell.losses,
                "reuse": RUNNER.REUSE,
                "iterations": RUNNER.RETAINED_SAMPLES,
                "warmup": RUNNER.WARMUP,
                "seed": cell.seed,
                "batch": 1, "thread_count": 1,
                "requested_profile": "legacy_high_v1",
                "requested_field": "gf16",
                "requested_backend": "avx2",
                "skip_legacy": True,
                "retain_samples": True,
                "report_decode_path": True,
            },
            "resolved": {
                "profile": "legacy_high_v1",
                "field": "gf16",
                "backend": "avx2",
                "thread_count": 1,
                "selected_decode_path": "tiled",
                "selected_decode_rule": "workspace_tiled",
                "decode_required_work_slots": 264,
            },
            "correctness": {"leopard2_round_trip": True},
            "workload_digests": {
                "algorithm": "fnv1a64",
                "original_data": "1" * 16,
                "transmitted_parity": "2" * 16,
                "recovered_originals": "3" * 16,
            },
            "metrics": {
                "encode_execution": dict(metric),
                "decode_execution": dict(metric),
            },
        }
        normalized = RUNNER.validate_timing_result(result, cell)
        self.assertEqual(
            len(normalized["decode_samples_us"]),
            RUNNER.RETAINED_SAMPLES)
        bad = json.loads(json.dumps(result))
        bad["resolved"]["selected_decode_path"] = "direct"
        with self.assertRaises(RUNNER.EvidenceError):
            RUNNER.validate_timing_result(bad, cell)


if __name__ == "__main__":
    unittest.main()
