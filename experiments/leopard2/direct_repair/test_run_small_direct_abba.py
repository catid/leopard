#!/usr/bin/env python3
"""Fail-closed scheduler-isolation tests for run_small_direct_abba.py."""

from __future__ import annotations

import argparse
import copy
import errno
import fcntl
import importlib.util
import json
import os
import signal
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock


RUNNER_PATH = Path(__file__).with_name("run_small_direct_abba.py")
SPEC = importlib.util.spec_from_file_location(
    "leopard2_test_direct_odd_runner", RUNNER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot import direct-repair ABBA runner")
RUNNER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUNNER
SPEC.loader.exec_module(RUNNER)


def record() -> dict:
    return {
        "pid": 101, "tid": 102, "command": "worker",
        "affinity": [0, 1, 2], "start_time_ticks": 12345,
        "task_inode": 777,
    }


def runner_record(affinity: list[int]) -> dict:
    return {
        "pid": 1001, "tid": 1001, "affinity": affinity,
        "start_time_ticks": 55555, "task_inode": 999,
    }


class IsolationTests(unittest.TestCase):
    def test_small_direct_matrices_are_complete_and_stable(self) -> None:
        core = RUNNER.make_matrix()
        self.assertEqual(core["schema"], RUNNER.MATRIX_SCHEMA)
        self.assertEqual(core["cell_count"], 160)
        self.assertEqual(
            core["matrix_sha256"],
            "745c53312d2da5a1e01ecc5c507f882d4849cf1017a6c782e95d0483e2585b67")
        self.assertEqual(
            core["matrix_sha256"], RUNNER.make_matrix()["matrix_sha256"])
        self.assertEqual(
            {cell["K"] for cell in core["cells"]}, {5, 8, 9, 12, 16})
        self.assertEqual(
            {cell["bytes"] for cell in core["cells"]},
            {64, 2048, 4096, 65536})
        self.assertEqual(
            [cell["index"] for cell in core["cells"]],
            list(range(core["cell_count"])))
        self.assertEqual(
            len({cell["id"] for cell in core["cells"]}),
            core["cell_count"])
        self.assertTrue(all(
            cell["loss"] <= min(cell["K"], cell["R"])
            for cell in core["cells"]))

        full = RUNNER.make_large_matrix()
        self.assertEqual(full["schema"], RUNNER.LARGE_MATRIX_SCHEMA)
        self.assertEqual(full["cell_count"], 3792)
        self.assertEqual(
            full["matrix_sha256"],
            "80fe9c95bd4e9908eba08f95aadd32e5783394d02a8ddc02190814feb32bdabf")
        self.assertEqual(
            full["matrix_sha256"],
            RUNNER.make_large_matrix()["matrix_sha256"])
        self.assertEqual(
            [cell["index"] for cell in full["cells"]],
            list(range(full["cell_count"])))
        self.assertEqual(
            len({cell["id"] for cell in full["cells"]}),
            full["cell_count"])
        expected = {
            (k, r, loss, byte_count, reuse)
            for k in range(5, 17)
            for r in range(5, 9)
            for loss in range(4, min(k, r) + 1)
            for byte_count in
                (64, 65, 256, 1024, 2048, 2049, 4096, 65536)
            for reuse in (1, 8, 64)
        }
        self.assertEqual({
            (cell["K"], cell["R"], cell["loss"], cell["bytes"],
             cell["reuse"])
            for cell in full["cells"]
        }, expected)

        tiny = RUNNER.make_tiny_matrix()
        self.assertEqual(tiny["schema"], RUNNER.TINY_MATRIX_SCHEMA)
        self.assertEqual(tiny["cell_count"], 96)
        self.assertEqual(
            tiny["matrix_sha256"],
            "97251d268a5d2a59e078ca5aff1ec064728176872d42213bad2218d58f7c9b1f")
        self.assertEqual(
            tiny["matrix_sha256"],
            RUNNER.make_tiny_matrix()["matrix_sha256"])
        self.assertEqual(
            [cell["index"] for cell in tiny["cells"]],
            list(range(tiny["cell_count"])))
        self.assertEqual(
            len({cell["id"] for cell in tiny["cells"]}),
            tiny["cell_count"])
        expected_tiny = {
            (k, r, loss, byte_count)
            for k, r, loss in (
                (5, 5, 4), (5, 5, 5),
                (8, 8, 4), (8, 8, 5), (8, 8, 8),
                (16, 8, 4), (16, 8, 5), (16, 8, 8))
            for byte_count in
                (1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63)
        }
        self.assertEqual({
            (cell["K"], cell["R"], cell["loss"], cell["bytes"])
            for cell in tiny["cells"]
        }, expected_tiny)
        self.assertFalse(any(
            cell["exact_main_required"] for cell in tiny["cells"]))
        for matrix in (core, full, tiny):
            self.assertEqual(
                matrix["directional_scope"],
                RUNNER.frozen_directional_scope())
            policy = matrix["same_source_promotion"]
            self.assertEqual(
                policy["neighbor_ci95_low_at_least"], 0.98)
            self.assertNotIn("neighbor_ci95_high_at_least", policy)

    def test_neighbor_promotion_uses_lower_confidence_endpoint(self) -> None:
        policy = {
            "target": "whole decode_execution",
            "candidate_ci95_low_at_least": 1.05,
            "neighbor_ci95_low_at_least": 0.98,
            "orientation": "baseline_time_over_candidate_time",
        }
        neighbor = {"id": "neighbor", "role": "loss4_neighbor"}
        target = {"id": "target", "role": "loss5_to_8_target"}
        matrix = {
            "same_source_promotion": policy,
            "directional_scope": {
                "baseline_candidate_pairs": [
                    ["transform", "output"],
                    ["transform", "source"],
                ],
                "direct_head_to_head": "directional comparison only",
            },
            "cells": [neighbor, target],
        }

        def summary(cell: dict, low: float, high: float) -> dict:
            return {
                "cell": cell,
                "metric_ratios": {
                    "execution": {"ci95": [low, high]},
                },
            }

        rejected = RUNNER.analyze_same_source_promotion(matrix, [
            summary(neighbor, 0.90, 1.20),
            summary(target, 1.06, 1.10),
        ], {"baseline": "transform", "candidate": "output"})
        self.assertFalse(rejected["checks"][0]["passes"])
        self.assertFalse(rejected["promotion_eligible"])

        accepted = RUNNER.analyze_same_source_promotion(matrix, [
            summary(neighbor, 0.98, 1.20),
            summary(target, 1.06, 1.10),
        ], {"baseline": "transform", "candidate": "output"})
        self.assertTrue(accepted["checks"][0]["passes"])
        self.assertTrue(accepted["same_source_screen_eligible"])
        self.assertTrue(accepted["same_source_direction_authorized"])
        self.assertFalse(accepted["promotion_eligible"])
        self.assertFalse(accepted["promotion_authorized"])
        self.assertEqual(
            accepted["exact_leopard1_evidence"]["status"], "not_bound")
        self.assertEqual(
            accepted["exact_leopard1_evidence"]["required_commit"],
            RUNNER.EXACT_LEOPARD1_COMMIT)

        reversed_result = RUNNER.analyze_same_source_promotion(matrix, [
            summary(neighbor, 1.20, 1.30),
            summary(target, 1.20, 1.30),
        ], {"baseline": "source", "candidate": "transform"})
        self.assertFalse(
            reversed_result["same_source_direction_authorized"])
        self.assertFalse(reversed_result["promotion_authorized"])
        self.assertFalse(reversed_result["promotion_eligible"])
        self.assertEqual(
            reversed_result["requested_mode_pair"],
            ["source", "transform"])
        self.assertIsNotNone(
            reversed_result["promotion_ineligibility_reason"])

        head_to_head = RUNNER.analyze_same_source_promotion(matrix, [
            summary(neighbor, 1.20, 1.30),
            summary(target, 1.20, 1.30),
        ], {"baseline": "output", "candidate": "source"})
        self.assertFalse(head_to_head["promotion_authorized"])
        self.assertFalse(
            head_to_head["same_source_direction_authorized"])
        self.assertFalse(head_to_head["promotion_eligible"])
        self.assertEqual(
            head_to_head["requested_mode_pair"], ["output", "source"])

        forged = json.loads(json.dumps(matrix))
        forged["same_source_promotion"][
            "neighbor_ci95_high_at_least"] = 0.98
        del forged["same_source_promotion"][
            "neighbor_ci95_low_at_least"]
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "invalid.*promotion contract"):
            RUNNER.analyze_same_source_promotion(
                forged, [],
                {"baseline": "transform", "candidate": "output"})

    def test_campaign_rejects_nonfinite_timeouts_before_use(self) -> None:
        for lock_timeout in (
                float("nan"), float("inf"), float("-inf")):
            with self.subTest(lock_timeout=lock_timeout), \
                    mock.patch.object(RUNNER, "acquire_lock") as acquire:
                options = argparse.Namespace(
                    lock_timeout=lock_timeout, lock=RUNNER.DEFAULT_LOCK)
                with self.assertRaisesRegex(
                        RUNNER.EvidenceError, "finite and nonnegative"):
                    RUNNER.run_campaign(options)
                acquire.assert_not_called()

        for timeout in (float("nan"), float("inf"), float("-inf")):
            with self.subTest(timeout=timeout):
                options = argparse.Namespace(
                    rounds=1, max_retries=0, timeout=timeout,
                    lock_timeout=0.0)
                with self.assertRaisesRegex(
                        RUNNER.EvidenceError, "finite and positive"):
                    RUNNER.run_campaign_locked(
                        options, None, {}, set())

    def test_small_direct_modes_are_explicit_and_directional(self) -> None:
        self.assertEqual(
            RUNNER.comparison_modes("transform", "output"),
            {"baseline": "transform", "candidate": "output"})
        self.assertEqual(
            RUNNER.comparison_modes("transform", "source"),
            {"baseline": "transform", "candidate": "source"})
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "distinct and known"):
            RUNNER.comparison_modes("source", "source")
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "distinct and known"):
            RUNNER.comparison_modes("unknown", "source")

        loss4 = {"loss": 4}
        loss5 = {"loss": 5}
        self.assertEqual(
            RUNNER.expected_direct_executor("transform", loss4),
            "output_major")
        self.assertEqual(
            RUNNER.expected_direct_executor("output", loss5),
            "output_major")
        self.assertEqual(
            RUNNER.expected_direct_executor("source", loss5),
            "source_major")
        self.assertEqual(
            RUNNER.expected_direct_executor("transform", loss5), "none")

        arguments = ["c++", "-O3", RUNNER.MODE_COMPILE_DEFINITIONS["source"]]
        self.assertEqual(
            RUNNER.strip_mode_definition(arguments, "source", "test"),
            ["c++", "-O3"])
        self.assertEqual(
            RUNNER.strip_mode_definition(
                ["c++", "-O3"], "transform", "test"),
            ["c++", "-O3"])
        self.assertEqual(RUNNER.mode_compile_arguments("transform"), [])
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "contains a diagnostic mode"):
            RUNNER.strip_mode_definition(
                arguments, "transform", "test")
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "exact mode definition once"):
            RUNNER.strip_mode_definition(["c++", "-O3"], "source", "test")
        alternate_spellings = (
            ["c++", "-D", RUNNER.MODE_CACHE_VARIABLE + "=2"],
            ["c++", "-D" + RUNNER.MODE_CACHE_VARIABLE + "=02"],
            ["c++", "-U" + RUNNER.MODE_CACHE_VARIABLE,
             RUNNER.MODE_COMPILE_DEFINITIONS["source"]],
            ["c++", "-Wp,-D" + RUNNER.MODE_CACHE_VARIABLE + "=2"],
            ["c++", RUNNER.MODE_COMPILE_DEFINITIONS["source"],
             "-U", RUNNER.MODE_CACHE_VARIABLE],
        )
        for alternate in alternate_spellings:
            with self.subTest(arguments=alternate), self.assertRaisesRegex(
                    RUNNER.EvidenceError, "exact mode definition once"):
                RUNNER.strip_mode_definition(alternate, "source", "test")

        benchmark_definitions = RUNNER.benchmark_attestation_definitions(
            "/tmp/leo2-build", "0" * 64)
        archive = {
            "role": "archive",
            "source": {"path": "/source/leopard2.cpp"},
            "compile_entry": {
                "output": "CMakeFiles/leopard.dir/leopard2.cpp.o",
                "arguments": [
                    "c++", "-O3",
                    RUNNER.MODE_COMPILE_DEFINITIONS["source"],
                ],
            },
        }
        benchmark = {
            "role": "benchmark",
            "compile_entry": {
                "output":
                    "CMakeFiles/bench_leopard2.dir/benchmark.cpp.o",
                "arguments": [
                    "c++", "-O3", *benchmark_definitions,
                ],
            },
        }
        backend_object = {
            "role": "archive",
            "source": {"path": "/source/Leopard2BackendAVX2.cpp"},
            "compile_entry": {
                "output":
                    "CMakeFiles/leopard2_backend_avx2.dir/"
                    "Leopard2BackendAVX2.cpp.o",
                "arguments": ["c++", "-O3"],
            },
        }
        self.assertEqual(
            RUNNER.normalized_compile_arguments(archive, "source"),
            ["c++", "-O3"])
        self.assertEqual(
            RUNNER.normalized_compile_arguments(
                benchmark, "source", benchmark_definitions),
            ["c++", "-O3"])
        self.assertEqual(
            RUNNER.normalized_compile_arguments(backend_object, "source"),
            ["c++", "-O3"])
        bad_benchmark = {
            "role": "benchmark",
            "compile_entry": {
                "output":
                    "CMakeFiles/bench_leopard2.dir/benchmark.cpp.o",
                "arguments": [
                    "c++", RUNNER.MODE_COMPILE_DEFINITIONS["source"],
                    *benchmark_definitions,
                ],
            },
        }
        with self.assertRaisesRegex(
                RUNNER.EvidenceError,
                "backend-object or benchmark.*diagnostic mode"):
            RUNNER.normalized_compile_arguments(
                bad_benchmark, "source", benchmark_definitions)
        for variable in (
                RUNNER.BENCHMARK_CONFIGURATION_DEFINITION,
                RUNNER.BENCHMARK_SOURCE_HEADER_DEFINITION):
            alternate_attestation = {
                "role": "benchmark",
                "compile_entry": {
                    "output":
                        "CMakeFiles/bench_leopard2.dir/benchmark.cpp.o",
                    "arguments": [
                        "c++", "-O3", *benchmark_definitions,
                        "-U", variable,
                    ],
                },
            }
            with self.subTest(
                    alternate_attestation=variable), self.assertRaisesRegex(
                    RUNNER.EvidenceError, "does not bind exact"):
                RUNNER.normalized_compile_arguments(
                    alternate_attestation, "source",
                    benchmark_definitions)

    def test_benchmark_v6_binds_direct_executor(self) -> None:
        cell = next(
            value for value in RUNNER.make_tiny_matrix()["cells"]
            if value["K"] == 8 and value["R"] == 8 and
            value["loss"] == 5 and value["bytes"] == 33)

        def result_for(mode: str) -> dict:
            executor = RUNNER.expected_direct_executor(mode, cell)
            direct = executor != "none"
            iterations = cell["iterations"]
            self.assertEqual(iterations, 15)
            factors = (
                0.70, 0.80, 0.85, 0.90, 0.94, 0.97, 0.99, 1.00,
                1.02, 1.05, 1.09, 1.14, 1.20, 1.30, 1.50,
            )

            def setup_summary(value: float) -> dict:
                return {
                    "median_us": value,
                    "mad_us": value * 0.10,
                    "minimum_us": value * 0.70,
                    "maximum_us": value * 1.50,
                    "samples_us": [
                        value * factor for factor in factors],
                }

            def execution_summary(
                    value: float, input_name: str, output_name: str,
                    input_bytes: int, output_bytes: int) -> dict:
                return {
                    "median_us_per_batch_call": value,
                    "mad_us_per_batch_call": value * 0.10,
                    "minimum_us_per_batch_call": value * 0.70,
                    "maximum_us_per_batch_call": value * 1.50,
                    "samples_us_per_batch_call": [
                        value * factor for factor in factors],
                    input_name: input_bytes / (value * 1000.0),
                    output_name: output_bytes / (value * 1000.0),
                }

            decode_us = 1.25
            plan_us = 0.5
            amortized_us = decode_us + plan_us / cell["reuse"]
            decode_input = (
                cell["K"] - cell["loss"] + cell["R"]) * cell["bytes"]
            decode_output = cell["loss"] * cell["bytes"]
            return {
                "schema": "leopard2-benchmark-v6",
                "parameters": {
                    "K": cell["K"], "R": cell["R"],
                    "shard_bytes": cell["bytes"],
                    "loss_count": cell["loss"],
                    "batch": cell["batch"], "reuse": cell["reuse"],
                    "iterations": cell["iterations"],
                    "warmup": cell["warmup"],
                    "thread_count": 1, "seed": cell["seed"],
                    "requested_profile": "legacy_high_v1",
                    "requested_field": "gf8",
                    "requested_backend": "avx2",
                    "force_generic_decode": False,
                    "force_specialized_decode": False,
                    "force_tiled_decode": False,
                    "force_materialized_decode": False,
                    "skip_legacy": True, "retain_samples": True,
                    "report_decode_path": True,
                    "report_direct_executor": True,
                    "missing_original_indices":
                        list(range(cell["loss"])),
                },
                "resolved": {
                    "profile": "legacy_high_v1", "field": "gf8",
                    "backend": "avx2", "thread_count": 1,
                    "parent_count": 16, "padded_side": 8,
                    "selected_decode_path":
                        "direct" if direct else "materialized",
                    "selected_decode_rule":
                        "direct" if direct else "high_specialized",
                    "decode_required_work_slots": 0 if direct else 16,
                    "decode_aligned_prefix_bytes": 0,
                    "decode_tail_bytes": 33,
                    "decode_rounded_bytes": 64,
                    "decode_multi_item_batch": False,
                    "selected_direct_executor": executor,
                },
                "correctness": {
                    "leopard2_round_trip": True,
                    "legacy_comparison": None,
                },
                "workload_digests": {
                    "algorithm": "fnv1a64",
                    "original_data": "0123456789abcdef",
                    "transmitted_parity": "123456789abcdef0",
                    "recovered_originals": "23456789abcdef01",
                },
                "metrics": {
                    "codec_setup": setup_summary(0.25),
                    "encode_execution": execution_summary(
                        2.0, "input_GB_per_s", "parity_output_GB_per_s",
                        cell["K"] * cell["bytes"],
                        cell["R"] * cell["bytes"]),
                    "decode_execution": execution_summary(
                        decode_us, "offered_received_GB_per_s",
                        "repaired_output_GB_per_s",
                        decode_input, decode_output),
                    "decode_plan_setup": setup_summary(plan_us),
                    "decode_amortized_at_reuse": {
                        "reuse_count": cell["reuse"],
                        "derived_median_us_per_batch_call": amortized_us,
                        "offered_received_GB_per_s":
                            decode_input / (amortized_us * 1000.0),
                        "repaired_output_GB_per_s":
                            decode_output / (amortized_us * 1000.0),
                    },
                    "rate_semantics":
                        "offered_received counts all non-null shard pointers "
                        "supplied; a plan may read a deterministic subset",
                },
            }

        self.assertIn(
            "--report-direct-executor", RUNNER.benchmark_arguments(cell))
        transform = RUNNER.validate_result(
            result_for("transform"), cell, "transform")
        source = RUNNER.validate_result(
            result_for("source"), cell, "source")
        self.assertEqual(transform["build_bound_executor"], "none")
        self.assertEqual(source["build_bound_executor"], "source_major")
        self.assertEqual(
            RUNNER.statistics.median(
                source["decode_execution_samples_us_per_batch_call"]),
            1.25)
        self.assertNotEqual(
            source["decode_execution_samples_us_per_batch_call"],
            [1.25] * cell["iterations"])

        stale = result_for("source")
        stale["schema"] = "leopard2-benchmark-v5"
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "direct-executor v6"):
            RUNNER.validate_result(stale, cell, "source")

        for label, mutate in (
                ("execution median", lambda value: value["metrics"][
                    "decode_execution"].update({
                        "median_us_per_batch_call": 0.01})),
                ("execution MAD", lambda value: value["metrics"][
                    "decode_execution"].update({
                        "mad_us_per_batch_call": 0.5})),
                ("setup median", lambda value: value["metrics"][
                    "decode_plan_setup"].update({"median_us": 9.0})),
                ("encode rate", lambda value: value["metrics"][
                    "encode_execution"].update({"input_GB_per_s": 999.0})),
                ("amortized median", lambda value: value["metrics"][
                    "decode_amortized_at_reuse"].update({
                        "derived_median_us_per_batch_call": 0.1}))):
            with self.subTest(forgery=label):
                forged = result_for("source")
                mutate(forged)
                with self.assertRaisesRegex(
                        RUNNER.EvidenceError,
                        "derived|retained|timing"):
                    RUNNER.validate_result(forged, cell, "source")

        nonfinite = result_for("source")
        nonfinite["metrics"]["decode_execution"][
            "median_us_per_batch_call"] = float("nan")
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "finite positive"):
            RUNNER.validate_result(nonfinite, cell, "source")

        overflowed_rate = result_for("source")
        encode = overflowed_rate["metrics"]["encode_execution"]
        smallest_positive = float.fromhex("0x0.0000000000001p-1022")
        encode.update({
            "median_us_per_batch_call": smallest_positive,
            "mad_us_per_batch_call": 0.0,
            "minimum_us_per_batch_call": smallest_positive,
            "maximum_us_per_batch_call": smallest_positive,
            "samples_us_per_batch_call":
                [smallest_positive] * cell["iterations"],
            "input_GB_per_s": 1.0,
            "parity_output_GB_per_s": 1.0,
        })
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "not derived"):
            RUNNER.validate_result(overflowed_rate, cell, "source")

        exact_type_mutations = (
            ("K", "K", float(cell["K"])),
            ("thread count", "thread_count", True),
            ("Boolean selector", "force_generic_decode", 0),
            ("reuse", "reuse", float(cell["reuse"])),
        )
        for label, key, replacement in exact_type_mutations:
            with self.subTest(exact_parameter=label):
                forged = result_for("source")
                forged["parameters"][key] = replacement
                with self.assertRaisesRegex(
                        RUNNER.EvidenceError,
                        "parameter .* does not match"):
                    RUNNER.validate_result(forged, cell, "source")
        forged = result_for("source")
        forged["parameters"]["unexpected"] = True
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "parameters have unexpected"):
            RUNNER.validate_result(forged, cell, "source")
        for label, key, replacement in (
                ("thread count", "thread_count", True),
                ("parent count", "parent_count", 16.0),
                ("multi-item batch", "decode_multi_item_batch", 0)):
            with self.subTest(exact_resolved=label):
                forged = result_for("source")
                forged["resolved"][key] = replacement
                with self.assertRaises(RUNNER.EvidenceError):
                    RUNNER.validate_result(forged, cell, "source")
        forged = result_for("source")
        forged["resolved"]["unexpected"] = 1
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "resolved fields"):
            RUNNER.validate_result(forged, cell, "source")

    def test_summary_reports_setup_and_first_use_separately(self) -> None:
        cell = RUNNER.make_matrix()["cells"][0]
        with tempfile.TemporaryDirectory(
                prefix="leo2-small-direct-summary-") as directory:
            root = Path(directory)
            invocations = []
            for index, implementation in enumerate(RUNNER.ORDER):
                envelope = root / ("%d.json" % index)
                envelope.write_text("{}\n")
                baseline = implementation == "baseline"
                invocations.append({
                    "implementation": implementation,
                    "execution_identity": {
                        "schema": RUNNER.EXECUTION_SCHEMA,
                        "command_path": "/proc/self/fd/9",
                    },
                    "normalized": {
                        "median_us": 10.0 if baseline else 5.0,
                        "plan_setup_us": 4.0 if baseline else 1.0,
                        "digests": {"workload": "same"},
                        "missing_original_indices": [0, 1, 2, 3],
                        "decode_path": "direct",
                        "decode_rule": "direct",
                        "build_bound_executor": "output_major",
                    },
                    "envelope_path": str(envelope),
                    "reserved_sibling_nonidle_jiffies": 0,
                    "same_user_pair_affinity_before": [],
                    "same_user_pair_affinity_after": [],
                    "target_runtime": {},
                    "target_interrupts": {},
                    "isolation_epoch_digest": "0" * 64,
                })
            summary = RUNNER.summarize_cell(
                cell, invocations, 1, {"digest": "a" * 64})
        self.assertEqual(summary["binary_identity_digest"], "a" * 64)
        self.assertNotIn("binary_identity", summary)
        self.assertEqual(
            summary["plan_setup_us"]["baseline"]["samples"], [4.0, 4.0])
        self.assertEqual(
            summary["plan_setup_us"]["candidate"]["median"], 1.0)
        self.assertAlmostEqual(
            summary["metric_ratios"]["execution"]
                ["geometric_mean_ratio"],
            2.0)
        self.assertAlmostEqual(
            summary["metric_ratios"]["first_use_plan_plus_execution"]
                ["geometric_mean_ratio"],
            14.0 / 6.0)

    def test_default_mode_source_contract_is_exact(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-default-mode-") as directory:
            root = Path(directory)
            source = root / "leopard2.cpp"
            source.write_text(
                "#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                "#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0\n"
                "#endif\n")
            RUNNER.validate_default_mode_source(root)
            source.write_text(
                "#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                "#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\n"
                "#endif\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "does not default"):
                RUNNER.validate_default_mode_source(root)
            source.write_text(
                "/*\n"
                "#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                "#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0\n"
                "#endif\n"
                "*/\n"
                "#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "does not default"):
                RUNNER.validate_default_mode_source(root)
            default_block = (
                "#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                "#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0\n"
                "#endif\n")
            source.write_text(
                default_block +
                "#un\\\ndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                "#def\\\nine LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "does not default"):
                RUNNER.validate_default_mode_source(root)
            for whitespace in (" ", "\t", "\v", "\f"):
                source.write_text(
                    default_block +
                    "#un\\" + whitespace +
                    "\ndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                    "#def\\" + whitespace +
                    "\nine LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\n")
                with self.subTest(
                        splice_whitespace=ord(whitespace)), \
                        self.assertRaisesRegex(
                            RUNNER.EvidenceError,
                            "whitespace after a line splice"):
                    RUNNER.validate_default_mode_source(root)
            source.write_text(
                default_block +
                "%:undef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                "%:define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "does not default"):
                RUNNER.validate_default_mode_source(root)
            source.write_text(
                default_block +
                "??=undef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                "??=define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "does not default"):
                RUNNER.validate_default_mode_source(root)
            for whitespace in ("\v", "\f"):
                source.write_text(
                    default_block +
                    "#" + whitespace +
                    "undef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                    "#" + whitespace +
                    "define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\n")
                with self.subTest(
                        directive_whitespace=ord(whitespace)), \
                        self.assertRaisesRegex(
                            RUNNER.EvidenceError, "does not default"):
                    RUNNER.validate_default_mode_source(root)
            source.write_text(
                default_block +
                "#undef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\r"
                "#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\r")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "does not default"):
                RUNNER.validate_default_mode_source(root)
            source.write_bytes(
                default_block.encode("ascii") +
                b"#\0undef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\n"
                b"#\0define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "control characters"):
                RUNNER.validate_default_mode_source(root)
            source.write_text(
                "\ufeff#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 2\n" +
                default_block)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "does not default"):
                RUNNER.validate_default_mode_source(root)
            source.write_text(
                "#if 0\n" + default_block + "#endif\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "does not default"):
                RUNNER.validate_default_mode_source(root)

    def test_effective_configuration_normalizes_only_mode(self) -> None:
        def attestation(mode: str) -> dict:
            entries = {
                name: "" for name in RUNNER.BUILD_CONFIGURATION_VARIABLES
            }
            entries.update({
                "CMAKE_BUILD_TYPE": "Release",
                RUNNER.MODE_CACHE_VARIABLE: mode,
                **RUNNER.REQUIRED_DISABLED_EXPERIMENTS,
            })
            return {
                "entries": entries,
                "path": "/fixture/effective-configuration.txt",
                "schema": RUNNER.BUILD_CONFIGURATION_ATTESTATION_SCHEMA,
                "sha256": RUNNER.build_configuration_digest(
                    entries, RUNNER.BUILD_CONFIGURATION_VARIABLES),
            }

        transform = attestation("0")
        source = attestation("2")
        self.assertEqual(
            RUNNER.normalized_effective_configuration(
                transform, "transform"),
            RUNNER.normalized_effective_configuration(source, "source"))
        self.assertEqual(
            RUNNER.normalized_validated_cache(
                transform["entries"], "transform"),
            RUNNER.normalized_validated_cache(source["entries"], "source"))

        forged = copy.deepcopy(source)
        forged["entries"][RUNNER.MODE_CACHE_VARIABLE] = "1"
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "wrong diagnostic mode"):
            RUNNER.normalized_effective_configuration(forged, "source")
        missing = copy.deepcopy(source)
        missing["entries"].pop(RUNNER.MODE_CACHE_VARIABLE)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "wrong diagnostic mode"):
            RUNNER.normalized_effective_configuration(missing, "source")
        for name in RUNNER.REQUIRED_DISABLED_EXPERIMENTS:
            forged = copy.deepcopy(source)
            forged["entries"][name] = "ON"
            with self.subTest(name=name), self.assertRaisesRegex(
                    RUNNER.EvidenceError,
                    "does not disable unrelated"):
                RUNNER.normalized_effective_configuration(forged, "source")
            with self.subTest(
                    validated_cache_name=name), self.assertRaisesRegex(
                    RUNNER.EvidenceError,
                    "does not disable unrelated"):
                RUNNER.normalized_validated_cache(
                    forged["entries"], "source")

    def test_historical_v2_attestation_replays_but_current_v3_is_strict(
            self) -> None:
        def make_attestation(
                schema: str, variables: tuple[str, ...],
                disabled: dict[str, str]) -> dict:
            entries = {name: "" for name in variables}
            entries.update({
                RUNNER.MODE_CACHE_VARIABLE: "2",
                **disabled,
            })
            return {
                "entries": entries,
                "path": "/retained/effective-configuration.txt",
                "schema": schema,
                "sha256": RUNNER.build_configuration_digest(
                    entries, variables),
            }

        historical = make_attestation(
            RUNNER.BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2,
            RUNNER.BUILD_CONFIGURATION_VARIABLES_V2,
            RUNNER.REQUIRED_DISABLED_EXPERIMENTS_V2)
        normalized = RUNNER.normalized_effective_configuration(
            historical, "source", RUNNER.BUILD_SCHEMA_V4)
        self.assertNotIn(
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT", normalized)
        historical_cache = dict(historical["entries"])
        RUNNER.normalized_validated_cache(
            historical_cache, "source", RUNNER.BUILD_SCHEMA_V4)

        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "framing"):
            RUNNER.normalized_effective_configuration(
                historical, "source", RUNNER.BUILD_SCHEMA)

        current = make_attestation(
            RUNNER.BUILD_CONFIGURATION_ATTESTATION_SCHEMA,
            RUNNER.BUILD_CONFIGURATION_VARIABLES,
            RUNNER.REQUIRED_DISABLED_EXPERIMENTS)
        general = "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"
        missing = copy.deepcopy(current)
        missing["entries"].pop(general)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "selector keys"):
            RUNNER.normalized_effective_configuration(
                missing, "source", RUNNER.BUILD_SCHEMA)
        missing_cache = dict(current["entries"])
        missing_cache.pop(general)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "versioned selector keys"):
            RUNNER.normalized_validated_cache(
                missing_cache, "source", RUNNER.BUILD_SCHEMA)

        enabled = copy.deepcopy(current)
        enabled["entries"][general] = "ON"
        enabled["sha256"] = RUNNER.build_configuration_digest(
            enabled["entries"], RUNNER.BUILD_CONFIGURATION_VARIABLES)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "does not disable unrelated"):
            RUNNER.normalized_effective_configuration(
                enabled, "source", RUNNER.BUILD_SCHEMA)

        extended_historical = copy.deepcopy(historical)
        extended_historical["entries"][general] = "OFF"
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "selector keys"):
            RUNNER.normalized_effective_configuration(
                extended_historical, "source", RUNNER.BUILD_SCHEMA_V4)
        extended_historical_cache = dict(historical_cache)
        extended_historical_cache[general] = "OFF"
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "versioned selector keys"):
            RUNNER.normalized_validated_cache(
                extended_historical_cache, "source",
                RUNNER.BUILD_SCHEMA_V4)

    def test_source_closure_binds_attestation_generator(self) -> None:
        self.assertIn(
            "cmake/GenerateBenchmarkSourceAttestation.cmake",
            RUNNER.SOURCE_FILES)

    def test_provenance_binds_c_and_cxx_compile_drivers(self) -> None:
        provenance, unused = RUNNER.load_current_provenance_modules(
            RUNNER_PATH.parents[3])
        with tempfile.TemporaryDirectory(
                prefix="leo2-compile-language-") as directory:
            root = Path(directory)
            c_compiler = root / "cc"
            cxx_compiler = root / "c++"
            launcher = root / "launcher"
            for path in (c_compiler, cxx_compiler, launcher):
                path.write_text("test\n")
            c_source = root / "abi.c"
            cxx_source = root / "codec.cpp"
            upper_cxx_source = root / "codec.C"
            unknown_source = root / "kernel.cu"
            for path in (
                    c_source, cxx_source, upper_cxx_source, unknown_source):
                path.write_text("test\n")
            c_compiler = c_compiler.resolve(strict=True)
            cxx_compiler = cxx_compiler.resolve(strict=True)

            provenance._require_compile_driver(
                [str(c_compiler)], c_source, c_compiler, cxx_compiler)
            provenance._require_compile_driver(
                [str(cxx_compiler)], cxx_source, c_compiler, cxx_compiler)
            provenance._require_compile_driver(
                [str(cxx_compiler)], upper_cxx_source,
                c_compiler, cxx_compiler)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "another command driver"):
                provenance._require_compile_driver(
                    [str(cxx_compiler)], c_source,
                    c_compiler, cxx_compiler)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "another command driver"):
                provenance._require_compile_driver(
                    [str(launcher), str(c_compiler)], c_source,
                    c_compiler, cxx_compiler)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError, "response file"):
                provenance._require_compile_driver(
                    [str(c_compiler), "@arguments.rsp"], c_source,
                    c_compiler, cxx_compiler)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError, "omits"):
                provenance._require_compile_driver(
                    [], c_source, c_compiler, cxx_compiler)
            with self.assertRaisesRegex(
                    provenance.BuildProvenanceError,
                    "unsupported source language"):
                provenance._require_compile_driver(
                    [str(cxx_compiler)], unknown_source,
                    c_compiler, cxx_compiler)

    def test_provenance_cmake_cache_is_exact_lf_and_typed(self) -> None:
        provenance, unused = RUNNER.load_current_provenance_modules(
            RUNNER_PATH.parents[3])
        valid = (
            "# generated cache\n"
            "CMAKE_AR:FILEPATH=/usr/bin/ar\n"
            "CMAKE_BUILD_TYPE:STRING=Release\n"
            "CMAKE_C_COMPILER:FILEPATH=/usr/bin/cc\n"
            "CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++\n"
            "CMAKE_GENERATOR:INTERNAL=Unix Makefiles\n"
            "LEOPARD_ENABLE_GF8:BOOL=ON\n"
            "LEO2_BACKEND_VARIANT:STRING=avx2\n"
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256:INTERNAL="
                + "0" * 64 + "\n"
            "LEO2_FLAG_MAVX2:INTERNAL=1\n"
            "CUSTOM_VALUE:STRING=before-after\n"
        )
        parsed = provenance.parse_cmake_cache(valid.encode("utf-8"))
        self.assertEqual(parsed["CUSTOM_VALUE"], "before-after")

        type_confusions = {
            "archiver": (
                "CMAKE_AR:STRING=/usr/bin/ar",
                "CMAKE_AR:FILEPATH=/usr/bin/ar"),
            "build profile": (
                "CMAKE_BUILD_TYPE:BOOL=Release",
                "CMAKE_BUILD_TYPE:STRING=Release"),
            "C compiler": (
                "CMAKE_C_COMPILER:STRING=/usr/bin/cc",
                "CMAKE_C_COMPILER:FILEPATH=/usr/bin/cc"),
            "C++ compiler": (
                "CMAKE_CXX_COMPILER:STRING=/usr/bin/c++",
                "CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++"),
            "field option": (
                "LEOPARD_ENABLE_GF8:STRING=ON",
                "LEOPARD_ENABLE_GF8:BOOL=ON"),
            "backend profile": (
                "LEO2_BACKEND_VARIANT:BOOL=avx2",
                "LEO2_BACKEND_VARIANT:STRING=avx2"),
            "attestation": (
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256:STRING=" +
                    "0" * 64,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256:INTERNAL=" +
                    "0" * 64),
            "compiler probe": (
                "LEO2_FLAG_MAVX2:BOOL=1",
                "LEO2_FLAG_MAVX2:INTERNAL=1"),
        }
        for label, (forged, canonical) in type_confusions.items():
            with self.subTest(label=label):
                payload = valid.replace(canonical, forged)
                with self.assertRaisesRegex(
                        provenance.BuildProvenanceError, "has type"):
                    provenance.parse_cmake_cache(payload.encode("utf-8"))

        malformed = {
            "CR delimiter": valid.replace("\n", "\r\n", 1).encode(),
            "NUL delimiter": valid.replace(
                "Release", "Release\0", 1).encode(),
            "duplicate": (
                valid + "CMAKE_BUILD_TYPE:STRING=Release\n").encode(),
            "unsupported type": valid.replace(
                "CUSTOM_VALUE:STRING", "CUSTOM_VALUE:EVIL").encode(),
            "missing type": valid.replace(
                "CUSTOM_VALUE:STRING", "CUSTOM_VALUE").encode(),
            "extra type delimiter": valid.replace(
                "CUSTOM_VALUE:STRING", "CUSTOM_VALUE:STRING:MORE").encode(),
            "invalid key": valid.replace(
                "CUSTOM_VALUE:STRING", "CUSTOM VALUE:STRING").encode(),
            "Unicode line separator": valid.replace(
                "before-after", "before\u2028after").encode(),
        }
        for label, payload in malformed.items():
            with self.subTest(label=label):
                with self.assertRaises(provenance.BuildProvenanceError):
                    provenance.parse_cmake_cache(payload)

    def test_clean_replay_proof_uses_shared_complete_validator(self) -> None:
        closure = {"source_root": "/source"}
        proof = {"schema": "leopard2-reproducible-build-proof/v2"}
        calls = []

        class SharedValidator:
            @staticmethod
            def validate_reproducible_build_proof(
                    value, retained, *, label: str) -> None:
                calls.append((value, retained, label))
                if value != proof:
                    raise RuntimeError("immutable replay proof differs")

        RUNNER.validate_reproducible_build_proof(
            proof, closure, "fixture",
            provenance_module=SharedValidator)
        self.assertEqual(calls, [(proof, closure, "fixture")])
        with self.assertRaisesRegex(
                RUNNER.EvidenceError,
                "reproducible-build proof.*immutable replay"):
            RUNNER.validate_reproducible_build_proof(
                {"schema": "leopard2-reproducible-build-proof/v1"},
                closure, "fixture",
                provenance_module=SharedValidator)

    def test_capture_current_build_runs_clean_replay_with_inherited_lock(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-clean-replay-binding-") as directory:
            root = Path(directory).resolve()
            source = root / "source"
            build = root / "build"
            source.mkdir()
            build.mkdir()
            executable = build / "bench_leopard2"
            executable.write_bytes(b"executable")
            executable.chmod(0o755)
            for relative in (
                    "cmake/Leopard2BenchmarkAttestation.cmake",
                    "cmake/GenerateBenchmarkSourceAttestation.cmake",
                    "tools/leopard2_build_provenance.py",
                    "tools/leopard2_direct_encode_crossover.py"):
                path = source / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(relative + "\n", encoding="utf-8")
            cache_path = build / "CMakeCache.txt"
            cache_path.write_text("fixture\n", encoding="utf-8")
            cache_identity = RUNNER.file_identity(cache_path)
            entries = {
                name: "" for name in RUNNER.BUILD_CONFIGURATION_VARIABLES
            }
            entries.update({
                RUNNER.MODE_CACHE_VARIABLE: "2",
                "LEO2_BACKEND_VARIANT": "avx2",
                **RUNNER.REQUIRED_DISABLED_EXPERIMENTS,
            })
            configuration_digest = RUNNER.build_configuration_digest(
                entries, RUNNER.BUILD_CONFIGURATION_VARIABLES)
            executable_identity = RUNNER.file_identity(executable)
            closure = {
                "cmake_cache": cache_identity,
                "source_root": str(source),
                "executable_target": "bench_leopard2",
                "compiler": {"sha256": "1" * 64},
                "source_object_compile_closure": [],
                "archive_member_identities": [],
                "archive": {"sha256": "2" * 64},
                "executable": executable_identity,
                "validated_cache": entries,
            }
            proof = {
                "schema": "leopard2-reproducible-build-proof/v2",
                "method":
                    "runner-owned-empty-directory-configure-build-byte-compare",
                "source_root": str(source),
                "executable_target": "bench_leopard2",
                "compiler_sha256": "1" * 64,
                "objects": [], "archive_members": [],
                "archive_sha256": "2" * 64,
                "executable_sha256": executable_identity["sha256"],
            }
            calls = []

            class Provenance:
                @staticmethod
                def candidate_build_provenance(
                        unused_build, unused_source, unused_binary,
                        unused_target, *, inherited_descriptors=()):
                    calls.append(("capture", inherited_descriptors))
                    return closure

                @staticmethod
                def file_snapshot(unused_path, unused_label):
                    return cache_identity, b"cache"

                @staticmethod
                def parse_cmake_cache(unused_payload):
                    return entries

                @staticmethod
                def verify_reproducible_candidate_build(
                        unused_closure, *, inherited_descriptors=()):
                    calls.append(("replay", inherited_descriptors))
                    return proof

                @staticmethod
                def validate_reproducible_build_proof(
                        value, retained, *, label: str):
                    self.assertIs(value, proof)
                    self.assertIs(retained, closure)
                    self.assertEqual(label, "source")

            class Attestation:
                BUILD_CONFIGURATION_RELATIVE_PATH = Path("effective.json")

                @staticmethod
                def cmake_build_metadata(unused_binary):
                    return {
                        "build_root": str(build),
                        "executable": executable_identity,
                        "effective_configuration_attestation": {
                            "entries": entries,
                            "path": str(build / "effective.json"),
                            "schema":
                                RUNNER.BUILD_CONFIGURATION_ATTESTATION_SCHEMA,
                            "sha256": configuration_digest,
                        },
                    }

                @staticmethod
                def validate_build_configuration_attestation(
                        unused_record, unused_path):
                    return None

            with mock.patch.object(
                    RUNNER, "load_current_provenance_modules",
                    return_value=(Provenance, Attestation)), \
                 mock.patch.object(
                    RUNNER, "benchmark_attestation_definitions",
                    return_value=("unused-one", "unused-two")):
                result = RUNNER.capture_current_build(
                    executable, source, "source",
                    inherited_descriptors=(9,))
            self.assertEqual(result["schema"], RUNNER.BUILD_SCHEMA)
            self.assertEqual(result["reproducible_build"], proof)
            self.assertEqual(
                calls, [("capture", (9,)), ("replay", (9,))])

    def test_offline_frozen_binaries_are_lane_bound_and_rehashed(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-retained-binaries-") as directory:
            root = Path(directory).resolve()
            records = {}
            for label, content in (
                    ("baseline", b"baseline"), ("candidate", b"candidate")):
                path = root / "artifacts" / label / "bench_leopard2"
                path.parent.mkdir(parents=True)
                path.write_bytes(content)
                path.chmod(0o555)
                records[label] = {
                    "executable": {
                        "frozen": RUNNER.file_identity(path),
                    },
                }
            identity = {
                "baseline": records["baseline"],
                "candidate": records["candidate"],
            }
            RUNNER.validate_retained_frozen_executables(identity, root)

            candidate = root / "artifacts/candidate/bench_leopard2"
            original_metadata = candidate.stat()
            candidate.chmod(0o755)
            candidate.write_bytes(b"candidaXe")
            os.utime(candidate, ns=(
                original_metadata.st_atime_ns, original_metadata.st_mtime_ns))
            candidate.chmod(0o555)
            self.assertEqual(
                candidate.stat().st_size,
                identity["candidate"]["executable"]["frozen"]["size"])
            self.assertEqual(
                candidate.stat().st_mtime_ns, original_metadata.st_mtime_ns)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "retained.*changed"):
                RUNNER.validate_retained_frozen_executables(identity, root)
            candidate.chmod(0o755)
            candidate.write_bytes(b"candidate")
            candidate.chmod(0o555)

            baseline = root / "artifacts/baseline/bench_leopard2"
            baseline.unlink()
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "retained.*changed"):
                RUNNER.validate_retained_frozen_executables(identity, root)
            baseline.write_bytes(b"baseline")
            baseline.chmod(0o555)

            outside = root / "outside"
            outside.write_bytes(b"baseline")
            outside.chmod(0o555)
            escaped = json.loads(json.dumps(identity))
            escaped["baseline"]["executable"]["frozen"] = \
                RUNNER.file_identity(outside)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "escaped.*artifact lane"):
                RUNNER.validate_retained_frozen_executables(escaped, root)

            hardlink = root / "hardlinked-baseline"
            os.link(baseline, hardlink)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "retained.*changed"):
                RUNNER.validate_retained_frozen_executables(identity, root)
            hardlink.unlink()

    def test_immutable_descriptor_ignores_frozen_path_replacement(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-immutable-direct-") as directory:
            root = Path(directory).resolve()
            identity = {}
            for label, output in (
                    ("baseline", "baseline-original"),
                    ("candidate", "candidate-original")):
                path = root / "artifacts" / label / "bench_leopard2"
                path.parent.mkdir(parents=True)
                path.write_text(
                    "#!/bin/sh\nprintf '%s\\n'\n" % output,
                    encoding="utf-8")
                path.chmod(0o555)
                identity[label] = {
                    "executable": {
                        "frozen": RUNNER.file_identity(path),
                    },
                }
            guard = RUNNER.ImmutableFrozenExecutables(root, identity)
            with guard as (execution, descriptors):
                RUNNER.validate_execution_identity(execution, identity)
                candidate_path = Path(
                    identity["candidate"]["executable"]["frozen"]["path"])
                candidate_path.unlink()
                candidate_path.write_text(
                    "#!/bin/sh\nprintf 'replacement-evil\\n'\n",
                    encoding="utf-8")
                candidate_path.chmod(0o555)
                descriptor = descriptors["candidate"]
                self.assertEqual(
                    fcntl.fcntl(descriptor, RUNNER.LINUX_F_GET_SEALS) &
                        RUNNER.LINUX_REQUIRED_EXECUTABLE_SEALS,
                    RUNNER.LINUX_REQUIRED_EXECUTABLE_SEALS)
                os.fchmod(descriptor, 0o700)
                writer = os.open(
                    execution["candidate"]["command_path"], os.O_WRONLY)
                try:
                    with self.assertRaises(OSError) as mutation:
                        os.pwrite(writer, b"replacement-evil", 0)
                    self.assertEqual(mutation.exception.errno, errno.EPERM)
                finally:
                    os.close(writer)
                    os.fchmod(descriptor, 0o500)
                completed = subprocess.run(
                    [execution["candidate"]["command_path"]],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    pass_fds=(descriptor,), check=False)
                self.assertEqual(completed.returncode, 0)
                self.assertEqual(
                    completed.stdout, b"candidate-original\n")
                self.assertNotIn(b"replacement-evil", completed.stdout)
                RUNNER.validate_execution_descriptor(
                    descriptors["candidate"], execution["candidate"])

    def test_gated_invocation_executes_bound_descriptor_not_swapped_path(
            self) -> None:
        source_root = RUNNER_PATH.parents[3]
        cpu, sibling = self.free_physical_pair(source_root)
        PairLease, unused = RUNNER.load_pair_lease(source_root)
        with tempfile.TemporaryDirectory(
                prefix="leo2-gated-immutable-") as directory:
            root = Path(directory).resolve()
            os.chmod(root, 0o700)
            identity = {}
            for label in ("baseline", "candidate"):
                path = root / "artifacts" / label / "bench_leopard2"
                path.parent.mkdir(parents=True)
                path.write_text(
                    "#!/bin/sh\nprintf '%s-original\\n'\n" % label,
                    encoding="utf-8")
                path.chmod(0o555)
                identity[label] = {
                    "executable": {
                        "frozen": RUNNER.file_identity(path),
                    },
                }
            with RUNNER.ImmutableFrozenExecutables(
                    root, identity) as (execution, descriptors):
                candidate_path = Path(
                    identity["candidate"]["executable"]["frozen"]["path"])
                candidate_path.unlink()
                candidate_path.write_text(
                    "#!/bin/sh\nprintf 'swapped-path\\n'\n",
                    encoding="utf-8")
                candidate_path.chmod(0o555)
                with PairLease(cpu, sibling):
                    result = RUNNER.run_gated_benchmark(
                        [execution["candidate"]["command_path"]],
                        cpu, sibling, 10.0,
                        root / "stdout", root / "stderr",
                        executable_descriptor=descriptors["candidate"],
                        execution_identity=execution["candidate"])
                self.assertEqual(result["return_code"], 0)
                self.assertEqual(
                    (root / "stdout").read_bytes(),
                    b"candidate-original\n")
                self.assertNotIn(
                    b"swapped-path", (root / "stdout").read_bytes())

    def test_timed_child_retains_flock_after_coordinator_sigkill(self) -> None:
        source_root = RUNNER_PATH.parents[3]
        cpu, unused_sibling = self.free_physical_pair(source_root)
        with tempfile.TemporaryDirectory(
                prefix="leo2-gated-lock-sigkill-") as directory:
            root = Path(directory).resolve()
            os.chmod(root, 0o700)
            lock_path = root / "campaign.lock"
            marker = root / "timed-child-started"
            coordinator = os.fork()
            if coordinator == 0:
                try:
                    descriptor = os.open(
                        lock_path, os.O_CREAT | os.O_RDWR | os.O_CLOEXEC,
                        0o600)
                    fcntl.flock(descriptor, fcntl.LOCK_EX)
                    code = (
                        "from pathlib import Path; import os,time; "
                        f"os.fstat({descriptor}); "
                        f"Path({str(marker)!r}).write_text(str(os.getpid())); "
                        "time.sleep(1.5)"
                    )
                    RUNNER.run_gated_benchmark(
                        [sys.executable, "-c", code],
                        cpu, unused_sibling, 5.0,
                        root / "stdout", root / "stderr",
                        campaign_lock_descriptor=descriptor)
                    os.close(descriptor)
                    os._exit(0)
                except BaseException:
                    os._exit(91)

            contender = -1
            try:
                deadline = time.monotonic() + 5.0
                while not marker.is_file() and time.monotonic() < deadline:
                    time.sleep(0.01)
                self.assertTrue(marker.is_file(), "timed child did not start")
                contender = os.open(lock_path, os.O_RDWR | os.O_CLOEXEC)
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(
                        contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                os.kill(coordinator, signal.SIGKILL)
                unused_pid, status = os.waitpid(coordinator, 0)
                coordinator = -1
                self.assertTrue(os.WIFSIGNALED(status))
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(
                        contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = False
                deadline = time.monotonic() + 5.0
                while time.monotonic() < deadline:
                    try:
                        fcntl.flock(
                            contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        acquired = True
                        break
                    except BlockingIOError:
                        time.sleep(0.02)
                self.assertTrue(
                    acquired, "timed child leaked the campaign lock")
                fcntl.flock(contender, fcntl.LOCK_UN)
            finally:
                if coordinator > 0:
                    try:
                        os.kill(coordinator, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    os.waitpid(coordinator, 0)
                if contender >= 0:
                    os.close(contender)

    def test_normal_gated_completion_releases_inherited_flock(self) -> None:
        source_root = RUNNER_PATH.parents[3]
        cpu, sibling = self.free_physical_pair(source_root)
        PairLease, unused = RUNNER.load_pair_lease(source_root)
        with tempfile.TemporaryDirectory(
                prefix="leo2-gated-lock-cleanup-") as directory:
            root = Path(directory).resolve()
            os.chmod(root, 0o700)
            lock_path = root / "campaign.lock"
            owner = os.open(
                lock_path, os.O_CREAT | os.O_RDWR | os.O_CLOEXEC, 0o600)
            contender = -1
            try:
                fcntl.flock(owner, fcntl.LOCK_EX)
                with PairLease(cpu, sibling):
                    result = RUNNER.run_gated_benchmark(
                        ["/bin/true"], cpu, sibling, 5.0,
                        root / "stdout", root / "stderr",
                        campaign_lock_descriptor=owner)
                self.assertEqual(result["return_code"], 0)
                contender = os.open(
                    lock_path, os.O_RDWR | os.O_CLOEXEC)
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(
                        contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                os.close(owner)
                owner = -1
                fcntl.flock(contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(contender, fcntl.LOCK_UN)
            finally:
                if owner >= 0:
                    os.close(owner)
                if contender >= 0:
                    os.close(contender)

    def test_invalid_execution_binding_creates_no_retained_outputs(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-invalid-execution-binding-") as directory:
            root = Path(directory).resolve()
            os.chmod(root, 0o700)
            stdout = root / "stdout"
            stderr = root / "stderr"
            for arguments in (
                    {"executable_descriptor": 3},
                    {"execution_identity": {"forged": True}}):
                with self.subTest(arguments=arguments), \
                        self.assertRaisesRegex(
                            RUNNER.EvidenceError,
                            "descriptor and identity must be supplied "
                            "together"):
                    RUNNER.run_gated_benchmark(
                        ["/bin/true"], 0, 1, 1.0, stdout, stderr,
                        **arguments)
                self.assertFalse(stdout.exists())
                self.assertFalse(stderr.exists())

    def test_gated_resource_setup_failures_restore_mask_and_fds(self) -> None:
        def descriptors() -> set[str]:
            return set(os.listdir("/proc/self/fd"))

        before_mask = signal.pthread_sigmask(signal.SIG_BLOCK, set())
        with tempfile.TemporaryDirectory(
                prefix="leo2-gated-selector-failure-") as directory:
            root = Path(directory).resolve()
            os.chmod(root, 0o700)
            before_fds = descriptors()
            with mock.patch.object(
                    RUNNER.selectors, "DefaultSelector",
                    side_effect=OSError(errno.EMFILE, "selector unavailable")), \
                 mock.patch.object(RUNNER.os, "fork") as fork, \
                 self.assertRaisesRegex(OSError, "selector unavailable"):
                RUNNER.run_gated_benchmark(
                    ["/bin/true"], 0, 1, 1.0,
                    root / "stdout", root / "stderr")
            fork.assert_not_called()
            self.assertEqual(descriptors(), before_fds)
            self.assertEqual(
                signal.pthread_sigmask(signal.SIG_BLOCK, set()),
                before_mask)
            self.assertFalse((root / "stdout").exists())
            self.assertFalse((root / "stderr").exists())

        with tempfile.TemporaryDirectory(
                prefix="leo2-gated-pipe-failure-") as directory:
            root = Path(directory).resolve()
            os.chmod(root, 0o700)
            before_fds = descriptors()
            with mock.patch.object(
                    RUNNER.os, "pipe2",
                    side_effect=OSError(errno.EMFILE, "pipe unavailable")), \
                 mock.patch.object(RUNNER.os, "fork") as fork, \
                 self.assertRaisesRegex(OSError, "pipe unavailable"):
                RUNNER.run_gated_benchmark(
                    ["/bin/true"], 0, 1, 1.0,
                    root / "stdout", root / "stderr")
            fork.assert_not_called()
            self.assertEqual(descriptors(), before_fds)
            self.assertEqual(
                signal.pthread_sigmask(signal.SIG_BLOCK, set()),
                before_mask)

        real_pthread_sigmask = signal.pthread_sigmask
        calls = 0

        def fail_after_block(how, mask):
            nonlocal calls
            calls += 1
            result = real_pthread_sigmask(how, mask)
            if calls == 2:
                raise OSError(errno.EINVAL, "mask update failed")
            return result

        with tempfile.TemporaryDirectory(
                prefix="leo2-gated-mask-failure-") as directory:
            root = Path(directory).resolve()
            os.chmod(root, 0o700)
            before_fds = descriptors()
            with mock.patch.object(
                    RUNNER.signal, "pthread_sigmask",
                    side_effect=fail_after_block), \
                 mock.patch.object(RUNNER.os, "fork") as fork, \
                 self.assertRaisesRegex(OSError, "mask update failed"):
                RUNNER.run_gated_benchmark(
                    ["/bin/true"], 0, 1, 1.0,
                    root / "stdout", root / "stderr")
            fork.assert_not_called()
            self.assertGreaterEqual(calls, 3)
            self.assertEqual(descriptors(), before_fds)
            self.assertEqual(
                signal.pthread_sigmask(signal.SIG_BLOCK, set()),
                before_mask)

    def test_current_identity_rejects_rebuild_with_unchanged_binary_records(
            self) -> None:
        def file_record(path: str, byte: str) -> dict:
            return {"path": path, "size": 1, "sha256": byte * 64}

        runner = file_record("/runner", "1")
        process_support = file_record("/process-support", "5")
        pair_lease = file_record("/pair-lease", "2")
        baseline_origin = file_record("/baseline-origin", "3")
        baseline_frozen = file_record("/baseline-frozen", "3")
        candidate_origin = file_record("/candidate-origin", "4")
        candidate_frozen = file_record("/candidate-frozen", "4")
        baseline_build = {
            "digest": "baseline-build", "reproducible_build": {}}
        candidate_build = {
            "digest": "candidate-build", "reproducible_build": {}}
        recorded_source = {"root": "/source", "snapshot": "unchanged"}
        identity = {
            "source": recorded_source,
            "runner": runner,
            "process_containment_source": process_support,
            "pair_lease_source": pair_lease,
            "comparison_modes": {
                "baseline": "transform", "candidate": "source",
            },
            "baseline": {
                "executable": {
                    "origin": baseline_origin, "frozen": baseline_frozen,
                },
                "build": baseline_build,
            },
            "candidate": {
                "executable": {
                    "origin": candidate_origin, "frozen": candidate_frozen,
                },
                "build": candidate_build,
            },
            "build_comparison": {"digest": "comparison"},
        }
        records = {
            value["path"]: value for value in (
                runner, process_support, pair_lease, baseline_origin,
                baseline_frozen, candidate_origin, candidate_frozen)
        }

        def current_file(path: Path) -> dict:
            return records[str(path)]

        with mock.patch.object(
                RUNNER, "validate_binary_identity_structure"), \
             mock.patch.object(
                 RUNNER, "source_identity", return_value=recorded_source), \
             mock.patch.object(
                 RUNNER, "file_identity", side_effect=current_file), \
             mock.patch.object(
                 RUNNER, "capture_current_build",
                 side_effect=[
                     baseline_build, {"digest": "silently-rebuilt"},
                 ]), \
             mock.patch.object(RUNNER, "compare_current_builds"), \
             self.assertRaisesRegex(
                 RUNNER.EvidenceError,
                 "current candidate build provenance changed"):
            RUNNER.validate_current_binary_identity(identity)

    def test_strict_evidence_json_rejects_duplicates_and_nonfinite(
            self) -> None:
        malformed = {
            "request duplicate":
                b'{"request":{"cpu":1,"cpu":2}}',
            "manifest duplicate":
                b'{"manifest":1,"manifest":2}',
            "accepted envelope duplicate":
                b'{"return_code":0,"return_code":1}',
            "nested non-finite":
                b'{"result":{"median_us":NaN}}',
            "positive infinity":
                b'{"value":Infinity}',
            "negative infinity":
                b'{"value":-Infinity}',
            "positive exponent overflow":
                b'{"value":1e9999}',
            "negative exponent overflow":
                b'{"value":-1e9999}',
            "invalid UTF-8":
                b'{"value":"' + bytes((0xff,)) + b'"}',
        }
        for label, payload in malformed.items():
            with self.subTest(label=label):
                with self.assertRaises(RUNNER.EvidenceError):
                    RUNNER.strict_json_bytes(payload, label)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "not strict JSON"):
            RUNNER.strict_json_bytes(
                b"[" * 10000 + b"0" + b"]" * 10000,
                "excessively nested")
        self.assertEqual(
            RUNNER.strict_json_bytes(
                b'{"nested":{"left":1,"right":2}}', "valid"),
            {"nested": {"left": 1, "right": 2}})

    def test_reservation_payload_requires_exact_keys_and_integer_cpus(
            self) -> None:
        valid = {
            "schema": "leopard2-cpu-reservation/v1",
            "status": "held",
            "benchmark_cpu": 1,
            "reserved_sibling": 17,
            "owner": "test",
            "nonce": "fixture",
        }
        RUNNER.validate_reservation_payload(valid, 1, 17)
        for label, mutate in (
                ("Boolean CPU",
                 lambda value: value.update({"benchmark_cpu": True})),
                ("floating sibling",
                 lambda value: value.update({"reserved_sibling": 17.0})),
                ("extra key",
                 lambda value: value.update({"extra": 1})),
                ("missing key", lambda value: value.pop("owner"))):
            with self.subTest(label=label):
                forged = dict(valid)
                mutate(forged)
                with self.assertRaises(RUNNER.EvidenceError):
                    RUNNER.validate_reservation_payload(forged, 1, 17)

    def test_accepted_envelope_rejects_boolean_return_code(self) -> None:
        keys = {
            "schema", "implementation", "command", "execution_identity",
            "started_ns", "ended_ns",
            "return_code", "stdout", "stderr", "reserved_sibling_before",
            "reserved_sibling_after", "reserved_sibling_nonidle_jiffies",
            "reservation_before", "reservation_after",
            "same_user_pair_affinity_before",
            "same_user_pair_affinity_after", "gate", "target_runtime",
            "target_interrupts", "sibling_runtime", "sibling_interrupts",
            "wait4_crosscheck", "isolation_epoch_digest", "accepted",
            "result", "normalized", "envelope_path",
        }
        value = {key: None for key in keys}
        envelope = Path("/canonical/cell/accepted.envelope.json")
        value.update({
            "schema": RUNNER.ENVELOPE_SCHEMA,
            "implementation": "baseline",
            "return_code": False,
            "accepted": True,
            "envelope_path": str(envelope),
        })
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "header is invalid"):
            RUNNER.validate_envelope(
                value, envelope, {}, "baseline", 0, 0, {}, "0" * 64,
                Path("/canonical"))

    def test_source_identity_rejects_coherently_redigested_mutations(
            self) -> None:
        source = {
            "root": "/canonical/source",
            "head": "1" * 40,
            "head_tree": "2" * 40,
            "status_short": [],
            "files": {
                relative: {"size": 1, "sha256": "3" * 64}
                for relative in RUNNER.SOURCE_FILES
            },
        }

        def seal(value: dict) -> dict:
            sealed = json.loads(json.dumps(value))
            sealed["snapshot_sha256"] = RUNNER.sha256_bytes(
                RUNNER.canonical_bytes(sealed))
            return sealed

        RUNNER.validate_source_identity_structure(seal(source))
        mutations = {}
        for label, mutate in (
                ("relative root",
                 lambda value: value.update(root="relative/source")),
                ("noncanonical root",
                 lambda value: value.update(root="/canonical/../escape")),
                ("malformed head",
                 lambda value: value.update(head="not-a-commit")),
                ("malformed tree",
                 lambda value: value.update(head_tree="2" * 41)),
                ("dirty status",
                 lambda value: value.update(status_short=[" M codec.cpp"]))):
            value = json.loads(json.dumps(source))
            mutate(value)
            mutations[label] = seal(value)
        unsafe = json.loads(json.dumps(source))
        removed = next(iter(unsafe["files"]))
        record = unsafe["files"].pop(removed)
        unsafe["files"]["../" + removed] = record
        mutations["unsafe source key"] = seal(unsafe)
        for label, value in mutations.items():
            with self.subTest(label=label):
                with self.assertRaises(RUNNER.EvidenceError):
                    RUNNER.validate_source_identity_structure(value)

    def test_canonical_json_rejects_nonfinite_and_huge_integer_input(
            self) -> None:
        for value in (
                float("nan"), float("inf"), float("-inf"),
                {"nested": [float("nan")]}):
            with self.subTest(value=value), self.assertRaisesRegex(
                    RUNNER.EvidenceError, "canonical finite JSON"):
                RUNNER.canonical_bytes(value)
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "not strict JSON"):
            RUNNER.strict_json_bytes(b"1" * 5000, "huge integer")

    def test_output_root_is_fresh_canonical_empty_and_0700(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-output-root-") as directory:
            root = Path(directory).resolve()
            source = root / "source"
            parent = root / "campaigns"
            source.mkdir()
            parent.mkdir()
            output = parent / "fresh"
            candidate, parent_identity = \
                RUNNER.validate_output_root_candidate(output, source)
            self.assertFalse(output.exists())
            created = RUNNER.create_output_root(candidate, parent_identity)
            self.assertEqual(created, output)
            self.assertEqual(created.resolve(strict=True), created)
            self.assertEqual(stat.S_IMODE(created.stat().st_mode), 0o700)
            self.assertEqual(list(created.iterdir()), [])

    def test_output_root_rejects_existing_data_without_mutation(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-output-existing-") as directory:
            root = Path(directory).resolve()
            source = root / "source"
            parent = root / "campaigns"
            source.mkdir()
            parent.mkdir()
            cases = {
                "empty": None,
                "matrix": ("matrix.json", b"caller matrix\n"),
                "extra": ("untracked.txt", b"caller evidence\n"),
            }
            for label, retained in cases.items():
                with self.subTest(label=label):
                    output = parent / label
                    output.mkdir()
                    if retained is not None:
                        name, payload = retained
                        (output / name).write_bytes(payload)
                    before = {
                        item.name: item.read_bytes()
                        for item in output.iterdir() if item.is_file()
                    }
                    with self.assertRaisesRegex(
                            RUNNER.EvidenceError, "must not already exist"):
                        RUNNER.validate_output_root_candidate(output, source)
                    after = {
                        item.name: item.read_bytes()
                        for item in output.iterdir() if item.is_file()
                    }
                    self.assertEqual(after, before)

    def test_output_root_rejects_symlinks_missing_parent_and_creation_race(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-output-path-") as directory:
            root = Path(directory).resolve()
            source = root / "source"
            parent = root / "campaigns"
            source.mkdir()
            parent.mkdir()

            target = root / "target"
            target.mkdir()
            output_link = parent / "output-link"
            output_link.symlink_to(target, target_is_directory=True)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "must not already exist"):
                RUNNER.validate_output_root_candidate(output_link, source)

            parent_link = root / "parent-link"
            parent_link.symlink_to(parent, target_is_directory=True)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "parent is not canonical"):
                RUNNER.validate_output_root_candidate(
                    parent_link / "escaped", source)
            self.assertFalse((parent / "escaped").exists())

            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "output path is invalid"):
                RUNNER.validate_output_root_candidate(
                    root / "missing-parent" / "campaign", source)
            self.assertFalse((root / "missing-parent").exists())

            raced = parent / "raced"
            candidate, parent_identity = \
                RUNNER.validate_output_root_candidate(raced, source)
            raced.mkdir()
            retained = raced / "matrix.json"
            retained.write_bytes(b"race winner\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "appeared before creation"):
                RUNNER.create_output_root(candidate, parent_identity)
            self.assertEqual(retained.read_bytes(), b"race winner\n")

    def test_atomic_json_never_replaces_and_enforces_bound(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-atomic-new-") as directory:
            root = Path(directory).resolve()
            retained = root / "retained.json"
            retained.write_bytes(b"caller bytes\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "destination already exists"):
                RUNNER.atomic_json(retained, {"replacement": True})
            self.assertEqual(retained.read_bytes(), b"caller bytes\n")
            self.assertEqual(
                [item.name for item in root.iterdir()], ["retained.json"])

            oversized = root / "oversized.json"
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "retained byte bound"):
                RUNNER.atomic_json(
                    oversized, {"payload": "x" * 100}, maximum_bytes=16)
            self.assertFalse(oversized.exists())

    def test_stable_snapshot_rejects_symlink_fifo_and_oversize_sparse_file(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-snapshot-types-") as directory:
            root = Path(directory).resolve()
            target = root / "target.json"
            target.write_text("{}\n")
            link = root / "link.json"
            link.symlink_to(target)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "not canonical"):
                RUNNER.stable_file_snapshot(
                    link, "symlink", 1024, require_canonical=True)

            fifo = root / "evidence.fifo"
            os.mkfifo(fifo)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "not a regular file"):
                RUNNER.stable_file_snapshot(
                    fifo, "FIFO", 1024, require_canonical=True)

            sparse = root / "oversize.json"
            with sparse.open("wb") as stream:
                stream.truncate(1025)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "retained byte bound"):
                RUNNER.strict_json_file(
                    sparse, "oversize sparse JSON", maximum_bytes=1024)

    def test_source_snapshot_rejects_symlinked_source_file(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-source-symlink-") as directory:
            root = Path(directory).resolve()
            target = root / "real.cpp"
            target.write_text("int source;\n")
            (root / "alias.cpp").symlink_to(target)
            with mock.patch.object(RUNNER, "SOURCE_FILES", ("alias.cpp",)), \
                 self.assertRaisesRegex(
                     RUNNER.EvidenceError, "not canonical"):
                RUNNER.source_identity(root)

    def test_stable_snapshot_rejects_replacement_during_read(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-snapshot-replace-") as directory:
            root = Path(directory).resolve()
            evidence = root / "evidence.bin"
            replacement = root / "replacement.bin"
            evidence.write_bytes(b"a" * ((2 << 20) + 1))
            replacement.write_bytes(b"b" * ((2 << 20) + 1))
            real_read = RUNNER.os.read
            replaced = False

            def replace_after_first_read(
                    descriptor: int, byte_count: int) -> bytes:
                nonlocal replaced
                payload = real_read(descriptor, byte_count)
                if payload and not replaced:
                    replaced = True
                    os.replace(replacement, evidence)
                return payload

            with mock.patch.object(
                    RUNNER.os, "read",
                    side_effect=replace_after_first_read), \
                 self.assertRaisesRegex(
                     RUNNER.EvidenceError, "changed|path was replaced"):
                RUNNER.stable_file_snapshot(
                    evidence, "replace race", 4 << 20,
                    require_canonical=True)
            self.assertTrue(replaced)

    def test_stable_snapshot_rejects_truncation_during_read(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-snapshot-truncate-") as directory:
            evidence = Path(directory).resolve() / "evidence.bin"
            evidence.write_bytes(b"a" * ((2 << 20) + 1))
            real_read = RUNNER.os.read
            truncated = False

            def truncate_after_first_read(
                    descriptor: int, byte_count: int) -> bytes:
                nonlocal truncated
                payload = real_read(descriptor, byte_count)
                if payload and not truncated:
                    truncated = True
                    os.truncate(evidence, 1)
                return payload

            with mock.patch.object(
                    RUNNER.os, "read",
                    side_effect=truncate_after_first_read), \
                 self.assertRaisesRegex(
                     RUNNER.EvidenceError, "truncated|changed"):
                RUNNER.stable_file_snapshot(
                    evidence, "truncate race", 4 << 20,
                    require_canonical=True)
            self.assertTrue(truncated)

    def test_campaign_inventory_rejects_unaudited_siblings(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-campaign-inventory-") as directory:
            root = Path(directory).resolve()
            for name in ("request.json", "matrix.json", "manifest.json"):
                (root / name).write_text("{}\n")
            for name in ("artifacts", "cells", "raw"):
                (root / name).mkdir()
                os.chmod(root / name, 0o700)
            for label in ("baseline", "candidate"):
                lane = root / "artifacts" / label
                lane.mkdir()
                os.chmod(lane, 0o700)
                (lane / "bench_leopard2").write_bytes(label.encode())
            (root / "cells" / "cell.json").write_text("{}\n")
            (root / "raw" / "cell").mkdir()
            os.chmod(root / "raw" / "cell", 0o700)
            RUNNER.validate_campaign_inventory(root, ["cell"])

            extra = root / "untracked.txt"
            extra.write_text("extra\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "root contains"):
                RUNNER.validate_campaign_inventory(root, ["cell"])
            extra.unlink()

            extra_raw = root / "raw" / "extra-cell"
            extra_raw.mkdir()
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "raw directory"):
                RUNNER.validate_campaign_inventory(root, ["cell"])
            extra_raw.rmdir()

            extra_lane = root / "artifacts" / "unreviewed"
            extra_lane.mkdir()
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "artifact lanes"):
                RUNNER.validate_campaign_inventory(root, ["cell"])

    def test_retry_chain_is_exact_and_declares_isolation_failure(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-retained-retries-") as directory:
            root = Path(directory).resolve()
            cell = {
                "id": "cell", "K": 5, "R": 5, "bytes": 64, "loss": 5,
                "batch": 1, "reuse": 8, "iterations": 3, "warmup": 1,
                "seed": 7, "role": "loss5_to_8_target",
            }
            raw = root / "raw" / cell["id"]
            raw.mkdir(parents=True)
            prefix = "round00-slot0-baseline-attempt"
            for attempt in (0, 1):
                stem = prefix + ("%02d" % attempt)
                (raw / (stem + ".envelope.json")).write_text("{}\n")
                (raw / (stem + ".stdout.json")).write_text("{}\n")
                (raw / (stem + ".stderr.txt")).write_text("")
            accepted = raw / (prefix + "01.envelope.json")
            attempts = RUNNER.retained_attempt_paths(
                accepted, cell, "baseline", 0, 0, 3, root)
            self.assertEqual(len(attempts), 2)

            missing = attempts[0][1]
            missing.unlink()
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "extra or missing"):
                RUNNER.retained_attempt_paths(
                    accepted, cell, "baseline", 0, 0, 3, root)
            missing.write_text("{}\n")

            extra = raw / (prefix + "02.envelope.json")
            extra.write_text("{}\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "extra or missing"):
                RUNNER.retained_attempt_paths(
                    accepted, cell, "baseline", 0, 0, 3, root)
            extra.unlink()

            outside = root / "escaped.envelope.json"
            outside.write_text("{}\n")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "escaped.*cell directory"):
                RUNNER.retained_attempt_paths(
                    outside, cell, "baseline", 0, 0, 3, root)

            retry_path, stdout_path, stderr_path = attempts[0]
            zero_stat = [100] * 10
            interrupts = json.loads(RUNNER.canonical_bytes(
                RUNNER.target_interrupt_evidence(
                    tuple(zero_stat), tuple(zero_stat))))
            gate = {
                "pid": 123, "session": 123, "cpu": 0,
                "start_time_ticks": 456, "task_inode": 789,
                "stopped_process_state": {
                    "state": "T", "pgrp": 123, "session": 123,
                    "processor": 0,
                },
                "zombie_process_state": {
                    "state": "Z", "pgrp": 123, "session": 123,
                    "processor": 0,
                },
                "retained_session_members": [{
                    "pid": 123, "state": "Z", "pgrp": 123,
                    "session": 123, "start_time_ticks": 456,
                }],
                "pre_target_child_cpu_time_ns": 0,
                "child_schedstat_before_ns": 100,
                "child_schedstat_after_ns": 110,
                "affinity": [0],
            }
            target_runtime = RUNNER.target_runtime_evidence(
                1000, 21000, 10, RUNNER.TARGET_RUNTIME_TOLERANCE_NS)
            self.assertFalse(target_runtime["accepted"])
            retry = {
                "schema": RUNNER.ENVELOPE_SCHEMA,
                "implementation": "baseline",
                "execution_identity": {
                    "command_path": "/proc/self/fd/9",
                    "descriptor": 9,
                },
                "command": RUNNER.benchmark_command(
                    Path("/proc/self/fd/9"), cell, 0),
                "started_ns": 1, "ended_ns": 2, "return_code": 0,
                "stdout": RUNNER.file_identity(stdout_path),
                "stderr": RUNNER.file_identity(stderr_path),
                "reserved_sibling_before": zero_stat,
                "reserved_sibling_after": zero_stat,
                "reserved_sibling_nonidle_jiffies": 0,
                "reservation_before": {"held": True},
                "reservation_after": {"held": True},
                "same_user_pair_affinity_before": [],
                "same_user_pair_affinity_after": [],
                "gate": gate,
                "target_runtime": target_runtime,
                "target_interrupts": interrupts,
                "sibling_runtime": {
                    "scheduler_before_ns": 200,
                    "scheduler_after_ns": 200,
                    "scheduler_delta_ns": 0,
                    "accepted": True,
                },
                "sibling_interrupts": interrupts,
                "wait4_crosscheck": {
                    "child_cpu_time_ns": 10,
                    "child_schedstat_delta_ns": 10,
                    "absolute_difference_ns": 0,
                    "tolerance_ns":
                        RUNNER.RUSAGE_CROSSCHECK_TOLERANCE_NS,
                    "accepted": True,
                },
                "isolation_epoch_digest": "1" * 64,
                "accepted": False,
                "envelope_path": str(retry_path),
            }
            retry_path.write_text(json.dumps(retry) + "\n")
            request = {
                "execution_identity": {
                    "baseline": {
                        "command_path": "/proc/self/fd/9",
                        "descriptor": 9,
                    },
                },
                "cpu": 0,
                "reservation": {"held": True},
            }
            RUNNER.validate_retry_envelope(
                retry, retry_path, cell, "baseline", 0, 0, 0, request,
                "1" * 64, root)

            retry["target_runtime"] = RUNNER.target_runtime_evidence(
                1000, 1010, 10, RUNNER.TARGET_RUNTIME_TOLERANCE_NS)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "no declared isolation failure"):
                RUNNER.validate_retry_envelope(
                    retry, retry_path, cell, "baseline", 0, 0, 0, request,
                    "1" * 64, root)

    @staticmethod
    def physical_pairs() -> list[tuple[int, int]]:
        allowed = RUNNER.cgroup_effective_cpus() & set(os.sched_getaffinity(0))
        pairs = []
        seen = set()
        for cpu in sorted(allowed):
            path = Path(
                "/sys/devices/system/cpu/cpu%d/topology/"
                "thread_siblings_list" % cpu)
            siblings = RUNNER.parse_cpu_list(path.read_text()) & allowed
            if len(siblings) == 2:
                sibling = next(value for value in siblings if value != cpu)
                pair = tuple(sorted((cpu, sibling)))
                if pair not in seen:
                    seen.add(pair)
                    pairs.append(pair)
        return pairs

    @classmethod
    def physical_pair(cls) -> tuple[int, int]:
        pairs = cls.physical_pairs()
        if pairs:
            return pairs[0]
        raise unittest.SkipTest("no process-visible two-thread SMT pair")

    @classmethod
    def free_physical_pair(cls, source_root: Path) -> tuple[int, int]:
        PairLease, unused = RUNNER.load_pair_lease(source_root)
        for cpu, sibling in cls.physical_pairs():
            try:
                with PairLease(cpu, sibling):
                    pass
                return cpu, sibling
            except Exception:
                continue
        raise unittest.SkipTest("no unleased process-visible SMT pair")

    def test_exact_masks_are_excluded_and_restored(self) -> None:
        state = {102: {0, 1, 2}}

        def get_affinity(tid: int) -> set[int]:
            return set(state[tid])

        def set_affinity(tid: int, cpus: set[int]) -> None:
            state[tid] = set(cpus)

        with mock.patch.object(
                RUNNER, "cgroup_effective_cpus",
                return_value={0, 1, 2, 3}), \
             mock.patch.object(
                RUNNER, "same_user_pair_affinity",
                side_effect=[[record()], [], []]), \
             mock.patch.object(
                RUNNER, "runner_affinity_record",
                side_effect=[runner_record([0, 1, 2, 3]),
                             runner_record([0, 2])]), \
             mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", side_effect=get_affinity), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity", side_effect=set_affinity):
            exclusion = RUNNER.exclude_same_user_from_pair(1, 3)
        self.assertEqual(exclusion["changed"][0]["affinity"], [0, 1, 2])
        self.assertEqual(exclusion["changed"][0]["after"], [0, 2])
        self.assertEqual(state[102], {0, 2})

        with mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER, "runner_affinity_record",
                return_value=runner_record([0, 1, 2, 3])), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", side_effect=get_affinity), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity", side_effect=set_affinity):
            restoration = RUNNER.restore_same_user_affinity(exclusion)
        self.assertEqual(state[102], {0, 1, 2})
        self.assertEqual(restoration["restored"][0]["affinity"], [0, 1, 2])
        self.assertEqual(restoration["failures"], [])

    def test_vanished_task_identity_is_rejected(self) -> None:
        with mock.patch.object(
                RUNNER, "current_task_identity",
                side_effect=FileNotFoundError("gone")):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "vanished before mutation"):
                RUNNER.require_current_task(record())

    def test_pid_reused_task_identity_is_rejected(self) -> None:
        with mock.patch.object(
                RUNNER, "current_task_identity", return_value=(99999, 777)):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "reused or replaced"):
                RUNNER.require_current_task(record())

    def test_same_tick_reused_task_inode_is_rejected(self) -> None:
        with mock.patch.object(
                RUNNER, "current_task_identity", return_value=(12345, 888)):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "inode was reused or replaced"):
                RUNNER.require_current_task(record())

    def test_failed_affinity_mutation_rolls_back_and_rejects(self) -> None:
        with mock.patch.object(
                RUNNER, "cgroup_effective_cpus",
                return_value={0, 1, 2, 3}), \
             mock.patch.object(
                RUNNER, "same_user_pair_affinity",
                return_value=[record()]), \
             mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", return_value={0, 1, 2}), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity",
                side_effect=PermissionError("denied")):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "failed to exclude"):
                RUNNER.exclude_same_user_from_pair(1, 3)

    def test_new_pair_eligible_task_is_captured_and_restored(self) -> None:
        state = {102: {0, 1, 2}, 202: {1}}
        foreign = {
            "pid": 201, "tid": 202, "command": "foreign",
            "affinity": [1], "start_time_ticks": 54321,
            "task_inode": 778,
        }

        def get_affinity(tid: int) -> set[int]:
            return set(state[tid])

        def set_affinity(tid: int, cpus: set[int]) -> None:
            state[tid] = set(cpus)

        with mock.patch.object(
                RUNNER, "cgroup_effective_cpus",
                return_value={0, 1, 2, 3}), \
             mock.patch.object(
                RUNNER, "same_user_pair_affinity",
                side_effect=[[record()], [foreign], [], []]), \
             mock.patch.object(
                RUNNER, "runner_affinity_record",
                side_effect=[runner_record([0, 1, 2, 3]),
                             runner_record([0, 2]),
                             runner_record([0, 1, 2, 3])]), \
             mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", side_effect=get_affinity), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity", side_effect=set_affinity):
            exclusion = RUNNER.exclude_same_user_from_pair(1, 3)
            self.assertEqual(len(exclusion["changed"]), 2)
            self.assertEqual(state[202], {0, 2})
            RUNNER.restore_same_user_affinity(exclusion)
        self.assertEqual(state[102], {0, 1, 2})
        self.assertEqual(state[202], {1})

    def test_restore_failure_is_rejected(self) -> None:
        changed = {**record(), "after": [0, 2]}
        with mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", return_value={0, 2}), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity",
                side_effect=PermissionError("denied")):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "restoration failed"):
                RUNNER.restore_same_user_affinity({"changed": [changed]})

    def test_cgroup_eligibility_is_distinct_from_inherited_mask(self) -> None:
        identity = RUNNER.validate_smt_pair_identity(
            1, 3, {0, 1, 3}, {0}, {1, 3})
        self.assertEqual(identity["cgroup_effective_cpus"], [0, 1, 3])
        self.assertEqual(identity["launch_affinity"], [0])
        self.assertEqual(identity["housekeeping_affinity"], [0])

    def test_transient_target_runtime_is_rejected(self) -> None:
        evidence = RUNNER.target_runtime_evidence(
            1_000_000, 2_000_000,
            1_000_000 - RUNNER.TARGET_RUNTIME_TOLERANCE_NS - 1)
        self.assertFalse(evidence["accepted"])
        self.assertEqual(
            evidence["unexplained_runtime_ns"],
            RUNNER.TARGET_RUNTIME_TOLERANCE_NS + 1)

    def test_runtime_accounting_tolerance_is_fixed(self) -> None:
        evidence = RUNNER.target_runtime_evidence(
            10, 10 + RUNNER.TARGET_RUNTIME_TOLERANCE_NS, 0)
        self.assertTrue(evidence["accepted"])
        self.assertEqual(evidence["tolerance_ns"], 5_000)

    def test_child_runtime_absent_from_target_is_rejected(self) -> None:
        evidence = RUNNER.target_runtime_evidence(
            100, 200, 200 + RUNNER.TARGET_RUNTIME_TOLERANCE_NS)
        self.assertFalse(evidence["accepted"])
        self.assertLess(evidence["signed_difference_ns"], 0)

    def test_target_interrupt_fields_are_rejected(self) -> None:
        before = (0,) * 10
        for index, name in ((5, "irq"), (6, "softirq"), (7, "steal"),
                            (8, "guest"), (9, "guest_nice")):
            with self.subTest(field=name):
                after = list(before)
                after[index] = 1
                evidence = RUNNER.target_interrupt_evidence(
                    before, tuple(after))
                self.assertFalse(evidence["accepted"])
                self.assertEqual(evidence["rejected_fields"][name], 1)

    def test_cleanup_never_signals_after_reap_even_if_pid_reused(self) -> None:
        with mock.patch.object(
                RUNNER, "require_child_identity_gone_after_reap",
                side_effect=RUNNER.EvidenceError("reused")), \
             mock.patch.object(RUNNER.os, "kill") as kill, \
             mock.patch.object(RUNNER.os, "killpg") as killpg:
            failures = RUNNER.cleanup_gated_child(1234, True, True, 0.01)
        self.assertTrue(any("reused" in failure for failure in failures))
        kill.assert_not_called()
        killpg.assert_not_called()

    def test_cleanup_reaps_when_retained_session_scan_fails(self) -> None:
        with mock.patch.object(
                RUNNER, "signal_retained_child_session", return_value=[]), \
             mock.patch.object(
                RUNNER, "waitid_exit_without_reap", return_value=object()), \
             mock.patch.object(
                RUNNER, "process_group_or_session_members",
                side_effect=PermissionError("proc denied")), \
             mock.patch.object(
                RUNNER.os, "wait4", return_value=(1234, 0, object())) as wait4, \
             mock.patch.object(
                RUNNER, "require_child_identity_gone_after_reap"):
            failures = RUNNER.cleanup_gated_child(1234, False, True, 0.01)
        self.assertTrue(any("proc denied" in failure for failure in failures))
        wait4.assert_called_once_with(1234, 0)
        self.assertFalse(any("leader was not reaped" in failure
                             for failure in failures))

    def test_cleanup_kills_descendant_before_releasing_leader(self) -> None:
        leader = {
            "pid": 1234, "state": "Z", "pgrp": 1234, "session": 1234,
            "start_time_ticks": 10,
        }
        descendant = {
            "pid": 1235, "state": "S", "pgrp": 1234, "session": 1234,
            "start_time_ticks": 11,
        }
        with mock.patch.object(
                RUNNER, "signal_retained_child_session",
                return_value=[]) as signal_child, \
             mock.patch.object(
                RUNNER, "waitid_exit_without_reap", return_value=object()), \
             mock.patch.object(
                RUNNER, "process_group_or_session_members",
                side_effect=[[leader, descendant], [leader]]), \
             mock.patch.object(
                RUNNER.os, "wait4", return_value=(1234, 0, object())), \
             mock.patch.object(
                RUNNER, "require_child_identity_gone_after_reap"):
            failures = RUNNER.cleanup_gated_child(1234, False, True, 0.01)
        self.assertEqual(failures, [])
        self.assertEqual(signal_child.call_count, 2)

    def test_waitid_proof_error_reaps_without_resignaling(self) -> None:
        with mock.patch.object(
                RUNNER, "signal_retained_child_session",
                return_value=[]) as signal_child, \
             mock.patch.object(
                RUNNER, "waitid_exit_without_reap",
                side_effect=RUNNER.EvidenceError("bad waitid")), \
             mock.patch.object(
                RUNNER, "wait4_until", return_value=(1234, 0, object())), \
             mock.patch.object(
                RUNNER, "require_child_identity_gone_after_reap"):
            failures = RUNNER.cleanup_gated_child(1234, False, True, 0.01)
        self.assertTrue(any("bad waitid" in failure for failure in failures))
        signal_child.assert_called_once_with(1234, True)

    def test_gated_output_limits_are_enforced_while_child_runs(self) -> None:
        source_root = RUNNER_PATH.parents[3]
        cpu, sibling = self.free_physical_pair(source_root)
        PairLease, unused = RUNNER.load_pair_lease(source_root)
        cases = (
            ("stdout", 1, RUNNER.MAX_BENCHMARK_STDOUT_BYTES),
            ("stderr", 2, RUNNER.MAX_BENCHMARK_STDERR_BYTES),
        )
        with PairLease(cpu, sibling):
            for label, descriptor, limit in cases:
                with self.subTest(stream=label), \
                        tempfile.TemporaryDirectory(
                            prefix="leo2-bounded-output-") as directory:
                    root = Path(directory).resolve()
                    os.chmod(root, 0o700)
                    stdout = root / "stdout"
                    stderr = root / "stderr"
                    program = (
                        "import os;"
                        "p=bytes(range(251));"
                        "s={limit}+65536;"
                        "b=(p*((s+len(p)-1)//len(p)))[:s];"
                        "n=0;"
                        "\nwhile n<len(b):"
                        "\n n+=os.write({descriptor},b[n:])"
                    ).format(limit=limit, descriptor=descriptor)
                    with self.assertRaisesRegex(
                            RUNNER.EvidenceError,
                            "%s exceeded retained byte limit" % label):
                        RUNNER.run_gated_benchmark(
                            [sys.executable, "-c", program],
                            cpu, sibling, 10.0, stdout, stderr)
                    retained = stdout if label == "stdout" else stderr
                    self.assertEqual(retained.stat().st_size, limit)
                    pattern = bytes(range(251))
                    expected = (
                        pattern * ((limit + len(pattern) - 1) //
                                   len(pattern)))[:limit]
                    self.assertEqual(retained.read_bytes(), expected)

    def test_gated_cleanup_contains_setsid_double_fork_on_all_exits(
            self) -> None:
        source_root = RUNNER_PATH.parents[3]
        cpu, sibling = self.free_physical_pair(source_root)
        PairLease, unused = RUNNER.load_pair_lease(source_root)
        program = r"""
import os,sys,time
mode,pidfile=sys.argv[1],sys.argv[2]
first=os.fork()
if first == 0:
    os.setsid()
    second=os.fork()
    if second > 0:
        os._exit(0)
    for descriptor in (0,1,2):
        try:
            os.close(descriptor)
        except OSError:
            pass
    with open(pidfile, "w", encoding="ascii") as stream:
        stream.write(str(os.getpid()))
        stream.flush()
        os.fsync(stream.fileno())
    time.sleep(30)
    os._exit(0)
os.waitpid(first, 0)
deadline=time.monotonic()+5
while not os.path.exists(pidfile) and time.monotonic() < deadline:
    time.sleep(0.001)
if mode == "timeout":
    time.sleep(30)
sys.exit(0 if mode == "success" else 7)
"""
        before = RUNNER.PROCESS_SUPPORT._get_child_subreaper()
        with PairLease(cpu, sibling):
            for mode in ("success", "failure", "timeout"):
                with self.subTest(mode=mode), \
                        tempfile.TemporaryDirectory(
                            prefix="leo2-double-fork-") as directory:
                    root = Path(directory).resolve()
                    os.chmod(root, 0o700)
                    pidfile = root / "daemon.pid"
                    timeout = 0.25 if mode == "timeout" else 10.0
                    if mode == "timeout":
                        with self.assertRaisesRegex(
                                RUNNER.EvidenceError,
                                "benchmark exceeded"):
                            RUNNER.run_gated_benchmark(
                                [sys.executable, "-c", program, mode,
                                 str(pidfile)],
                                cpu, sibling, timeout,
                                root / "stdout", root / "stderr")
                    else:
                        result = RUNNER.run_gated_benchmark(
                            [sys.executable, "-c", program, mode,
                             str(pidfile)],
                            cpu, sibling, timeout,
                            root / "stdout", root / "stderr")
                        self.assertEqual(
                            result["return_code"],
                            0 if mode == "success" else 7)
                    self.assertTrue(pidfile.is_file())
                    daemon = int(pidfile.read_text(encoding="ascii"))
                    self.assertFalse(
                        (Path("/proc") / str(daemon)).exists())
                    with self.assertRaises(ProcessLookupError):
                        os.kill(daemon, 0)
                    self.assertEqual(
                        RUNNER.PROCESS_SUPPORT._get_child_subreaper(), before)

    def test_pending_signal_aborts_after_child_before_next_slot(self) -> None:
        cell = RUNNER.make_matrix()["cells"][0]
        reservation = {"frozen": True}
        with tempfile.TemporaryDirectory(
                prefix="leo2-pending-slot-") as directory, \
             mock.patch.object(
                RUNNER, "reservation_identity", return_value=reservation), \
             mock.patch.object(
                RUNNER, "same_user_pair_affinity", return_value=[]), \
             mock.patch.object(
                RUNNER, "run_gated_benchmark", return_value={}) as gated, \
             mock.patch.object(
                RUNNER.signal, "sigpending", return_value={signal.SIGTERM}):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "pending control signal"):
                RUNNER.run_invocation(
                    Path("/bin/true"), "baseline", cell, "transform",
                    1, 17, 1.0,
                    Path(directory), 0, 0, 0, Path("/unused"), reservation,
                    "0" * 64, set(), 0)
        gated.assert_called_once()

    def test_control_signals_wait_for_exact_child_cleanup(self) -> None:
        source_root = RUNNER_PATH.parents[3]
        cpu, sibling = self.free_physical_pair(source_root)
        PairLease, unused = RUNNER.load_pair_lease(source_root)
        helper = (
            "import importlib.util,signal,sys;"
            "[signal.signal(s,signal.SIG_DFL) for s in "
            "(signal.SIGHUP,signal.SIGINT,signal.SIGTERM)];"
            "p=sys.argv[1];"
            "s=importlib.util.spec_from_file_location('runner_signal_test',p);"
            "m=importlib.util.module_from_spec(s);"
            "sys.modules[s.name]=m;"
            "s.loader.exec_module(m);"
            "m.run_gated_benchmark(['/bin/sleep','10'],int(sys.argv[2]),"
            "int(sys.argv[3]),0.25,m.Path(sys.argv[4]),m.Path(sys.argv[5]))"
        )
        with PairLease(cpu, sibling):
            for control_signal in (
                    signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
                with self.subTest(signal=control_signal), \
                        tempfile.TemporaryDirectory(
                            prefix="leo2-gated-signal-") as directory:
                    root = Path(directory)
                    process = subprocess.Popen([
                        sys.executable, "-c", helper, str(RUNNER_PATH),
                        str(cpu), str(sibling), str(root / "stdout"),
                        str(root / "stderr"),
                    ])
                    child_pid = None
                    deadline = time.monotonic() + 5.0
                    children_path = Path("/proc") / str(process.pid) / \
                        "task" / str(process.pid) / "children"
                    while time.monotonic() < deadline:
                        if process.poll() is not None:
                            break
                        try:
                            children = children_path.read_text().split()
                            if children:
                                candidate = int(children[0])
                                if os.getsid(candidate) == candidate:
                                    child_pid = candidate
                                    break
                        except (FileNotFoundError, ProcessLookupError):
                            pass
                        time.sleep(0.001)
                    self.assertIsNotNone(
                        child_pid, "gated child did not establish its session")
                    os.kill(process.pid, control_signal)
                    return_code = process.wait(timeout=5.0)
                    self.assertEqual(return_code, -control_signal)
                    self.assertFalse(
                        (Path("/proc") / str(child_pid)).exists())
                    self.assertEqual(
                        RUNNER.process_group_or_session_members(child_pid), [])

    def test_campaign_signal_restores_affinity_and_releases_locks(self) -> None:
        source_root = RUNNER_PATH.parents[3]
        cpu, sibling = self.free_physical_pair(source_root)
        helper = """
import importlib.util,json,os,signal,sys,time
from pathlib import Path
[signal.signal(s, signal.SIG_DFL) for s in
 (signal.SIGHUP, signal.SIGINT, signal.SIGTERM)]
runner_path = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location('runner_campaign_signal', runner_path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
cpu, sibling = int(sys.argv[2]), int(sys.argv[3])
module.DEFAULT_LOCK = Path(sys.argv[4])
status_path = Path(sys.argv[5])
source_root = Path(sys.argv[6])
stdout_path, stderr_path = Path(sys.argv[7]), Path(sys.argv[8])
read_descriptor, write_descriptor = os.pipe()
first = os.fork()
if first == 0:
    os.close(read_descriptor)
    os.setsid()
    second = os.fork()
    if second > 0:
        os._exit(0)
    os.write(write_descriptor, str(os.getpid()).encode('ascii'))
    os.close(write_descriptor)
    time.sleep(10)
    os._exit(0)
os.close(write_descriptor)
os.waitpid(first, 0)
victim_pid = int(os.read(read_descriptor, 64))
os.close(read_descriptor)
original = sorted(os.sched_getaffinity(victim_pid))
changed = sorted(set(original) - {cpu, sibling})
if not changed or changed == original:
    raise RuntimeError('victim affinity does not exercise restoration')
PairLease, unused = module.load_pair_lease(source_root)
class Options:
    lock_timeout = 5.0
    lock = module.DEFAULT_LOCK
def fake_campaign(options, lock_stream, lock_identity, child_signal_mask):
    guard = PairLease(cpu, sibling)
    guard.__enter__()
    try:
        os.sched_setaffinity(victim_pid, set(changed))
        module.atomic_json(status_path, {
            'stage': 'mask_changed', 'victim_pid': victim_pid,
            'original': original, 'changed': changed})
        return module.run_gated_benchmark(
            ['/bin/sleep', '10'], cpu, sibling, 0.25,
            stdout_path, stderr_path, child_signal_mask,
            campaign_lock_descriptor=lock_stream.fileno())
    finally:
        os.sched_setaffinity(victim_pid, set(original))
        guard.__exit__(None, None, None)
module.run_campaign_locked = fake_campaign
module.run_campaign(Options())
"""
        for control_signal in (signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
            with self.subTest(signal=control_signal), \
                    tempfile.TemporaryDirectory(
                        prefix="leo2-campaign-signal-") as directory:
                root = Path(directory)
                test_lock = root / "campaign.lock"
                status_path = root / "status.json"
                process = subprocess.Popen([
                    sys.executable, "-c", helper, str(RUNNER_PATH),
                    str(cpu), str(sibling), str(test_lock), str(status_path),
                    str(source_root), str(root / "stdout"),
                    str(root / "stderr"),
                ])
                victim_pid = None
                gated_pid = None
                try:
                    deadline = time.monotonic() + 5.0
                    while time.monotonic() < deadline:
                        if process.poll() is not None:
                            break
                        try:
                            status = json.loads(status_path.read_text())
                            if status.get("stage") == "mask_changed":
                                victim_pid = status["victim_pid"]
                                break
                        except (FileNotFoundError, json.JSONDecodeError):
                            pass
                        time.sleep(0.001)
                    self.assertIsNotNone(
                        victim_pid, "campaign did not mutate victim mask")
                    children_path = Path("/proc") / str(process.pid) / \
                        "task" / str(process.pid) / "children"
                    deadline = time.monotonic() + 5.0
                    while time.monotonic() < deadline:
                        if process.poll() is not None:
                            break
                        try:
                            children = children_path.read_text().split()
                        except FileNotFoundError:
                            break
                        for raw_pid in children:
                            candidate = int(raw_pid)
                            if candidate != victim_pid and \
                                    os.getsid(candidate) == candidate:
                                gated_pid = candidate
                                break
                        if gated_pid is not None:
                            break
                        time.sleep(0.001)
                    self.assertIsNotNone(
                        gated_pid, "campaign gated child is unavailable")
                    os.kill(process.pid, control_signal)
                    self.assertEqual(
                        process.wait(timeout=5.0), -control_signal)
                    self.assertEqual(
                        sorted(os.sched_getaffinity(victim_pid)),
                        status["original"])
                    self.assertFalse(
                        (Path("/proc") / str(gated_pid)).exists())
                    self.assertEqual(
                        RUNNER.process_group_or_session_members(gated_pid), [])

                    descriptor = os.open(
                        test_lock,
                        os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW)
                    try:
                        fcntl.flock(
                            descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        fcntl.flock(descriptor, fcntl.LOCK_UN)
                    finally:
                        os.close(descriptor)
                    PairLease, unused = RUNNER.load_pair_lease(source_root)
                    with PairLease(cpu, sibling):
                        pass
                finally:
                    if process.poll() is None:
                        process.kill()
                        process.wait(timeout=5.0)
                    for pid in (gated_pid, victim_pid):
                        if pid is None:
                            continue
                        try:
                            os.kill(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass

if __name__ == "__main__":
    unittest.main()
