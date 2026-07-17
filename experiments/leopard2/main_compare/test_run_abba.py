#!/usr/bin/env python3
"""Self and mutation tests for the exact-main ABBA evidence runner."""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import math
import os
import stat
import sys
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("run_abba.py")
SPEC = importlib.util.spec_from_file_location("main_compare_run_abba", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


CELL = runner.Cell("fixture", 8, 4, 64, 2, 7)
SPECIFICATION = {
    "runner": "/fixture/run_abba.py",
    "taskset": "/usr/bin/taskset",
    "ldd": "/usr/bin/ldd",
    "baseline_executable": "/fixture/baseline",
    "candidate_executable": "/fixture/candidate",
    "baseline_archive": "/fixture/baseline.a",
    "candidate_archive": "/fixture/candidate.a",
    "baseline_build_dir": "/fixture/baseline-build",
    "candidate_build_dir": "/fixture/candidate-build",
    "baseline_source_root": "/fixture/main",
    "candidate_source_root": "/fixture/candidate-source",
    "candidate_commit": "a" * 40,
}
CAMPAIGN = {
    "rounds": runner.ROUNDS,
    "order": list(runner.ORDER),
    "cells": [asdict(CELL)],
    "batch": 1,
    "reuse": 8,
    "iterations": 3,
    "warmup": 2,
    "threads": 1,
    "child_environment": copy.deepcopy(runner.CHILD_ENVIRONMENT),
    "benchmark_cpu": 0,
    "reserved_sibling": 1,
    "allowed_cpu_set_at_launch": [0, 1, 2],
    "timeout_seconds": 10.0,
    "statistics": runner.statistics_policy(),
}
DIGESTS = {
    "algorithm": "fnv1a64",
    "original_data": "0123456789abcdef",
    "transmitted_parity": "1111111111111111",
    "recovered_originals": "fedcba9876543210",
}
RESERVATION_PAYLOAD = {
    "benchmark_cpu": 0,
    "nonce": "fixture-nonce",
    "owner": "unit test",
    "reserved_sibling": 1,
    "schema": runner.RESERVATION_SCHEMA,
    "status": "held",
}
RESERVATION = {
    "path": "/fixture/reservation.json",
    "sha256": runner.sha256_bytes(runner.canonical_bytes(RESERVATION_PAYLOAD)),
    "payload": RESERVATION_PAYLOAD,
    "lock": "exclusive_nonblocking",
}


def host_cpu(cpu: int) -> dict:
    return {
        "cpu": cpu,
        "online": "1",
        "cpuinfo": {"processor": str(cpu), "model name": "fixture cpu"},
        "topology": {
            "core_id": "0", "physical_package_id": "0", "die_id": "0",
            "cluster_id": None, "thread_siblings_list": "0-1",
            "core_siblings_list": "0-2",
        },
        "frequency_policy": {
            "scaling_driver": "fixture", "scaling_governor": "performance",
            "energy_performance_preference": "performance",
            "scaling_min_freq": "1", "scaling_max_freq": "2",
            "cpuinfo_min_freq": "1", "cpuinfo_max_freq": "2",
        },
    }


HOST = {
    "system": {
        "system": "Linux", "node": "fixture", "release": "fixture",
        "version": "fixture", "machine": "x86_64", "python": "3",
        "page_size": 4096,
    },
    "allowed_cpu_set_at_launch": [0, 1, 2],
    "online_cpu_set": [0, 1, 2],
    "benchmark_cpu": host_cpu(0),
    "reserved_sibling": host_cpu(1),
    "turbo_and_pstate": {"fixture": "0"},
}


def cpu_stat(cpu: int, *, user: int, idle: int) -> dict:
    fields = {
        "user": user, "nice": 0, "system": 0, "idle": idle,
        "iowait": 0, "irq": 0, "softirq": 0, "steal": 0,
    }
    return {"cpu": cpu, "fields": fields, "total_jiffies": sum(fields.values())}


PAIR_LEASE_PAYLOAD = runner.pair_lease_payload(0, 1)
PAIR_LEASE = {
    "device": 1,
    "directory_device": 1,
    "directory_inode": 2,
    "inode": 3,
    "lock": "exclusive_nonblocking_pair_wide",
    "path": str(runner.pair_lease_directory() / runner.pair_lease_name(0, 1)),
    "payload": PAIR_LEASE_PAYLOAD,
    "sha256": runner.sha256_bytes(runner.canonical_bytes(PAIR_LEASE_PAYLOAD)),
}
ISOLATION = runner.isolation_record(
    0, 1, PAIR_LEASE, 1_000, 2_000,
    cpu_stat(0, user=100, idle=100), cpu_stat(0, user=110, idle=110),
    cpu_stat(1, user=100, idle=100), cpu_stat(1, user=100, idle=120),
)


def summary(samples: list[float], setup: bool = False) -> dict:
    middle = sorted(samples)[len(samples) // 2]
    deviations = sorted(abs(value - middle) for value in samples)
    suffix = "" if setup else "_per_batch_call"
    result = {
        f"median_us{suffix}": middle,
        f"mad_us{suffix}": deviations[len(deviations) // 2],
        f"minimum_us{suffix}": min(samples),
        f"maximum_us{suffix}": max(samples),
        "samples_us" if setup else "samples_us_per_batch_call": samples,
    }
    return result


def common_parameters() -> dict:
    return {
        "K": CELL.k,
        "R": CELL.r,
        "shard_bytes": CELL.shard_bytes,
        "loss_count": CELL.losses,
        "missing_original_indices": runner.expected_losses(CELL),
        "batch": 1,
        "reuse": CAMPAIGN["reuse"],
        "iterations": CAMPAIGN["iterations"],
        "warmup": CAMPAIGN["warmup"],
        "thread_count": 1,
        "seed": CELL.seed,
    }


def baseline_result(scale: float = 1.0) -> dict:
    encode = summary([10.0 * scale, 11.0 * scale, 12.0 * scale])
    decode = summary([20.0 * scale, 21.0 * scale, 22.0 * scale])
    encode.update({"input_GB_per_s": 1.0, "parity_output_GB_per_s": 0.5})
    decode.update({"offered_received_GB_per_s": 1.0,
                   "repaired_output_GB_per_s": 0.5})
    return {
        "schema": "leopard-main-benchmark-v1",
        "build": {"main_source_commit": runner.MAIN_COMMIT, "cplusplus": 201103},
        "parameters": common_parameters(),
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "parent_count": 16, "padded_side": 4, "thread_count": 1,
        },
        "correctness": {"round_trip": True},
        "workload_digests": copy.deepcopy(DIGESTS),
        "metrics": {
            "codec_setup": None,
            "decode_timing_includes_setup": True,
            "encode_execution": encode,
            "decode_including_setup": decode,
            "rate_semantics": "fixture",
        },
    }


def candidate_result(scale: float = 0.8) -> dict:
    parameters = common_parameters()
    parameters.update({
        "requested_profile": "legacy_high_v1",
        "requested_field": "auto",
        "requested_backend": "auto",
        "force_generic_decode": False,
        "force_specialized_decode": False,
        "skip_legacy": True,
        "retain_samples": True,
    })
    codec = summary([3.0, 3.1, 3.2], setup=True)
    plan = summary([4.0, 4.1, 4.2], setup=True)
    encode = summary([10.0 * scale, 11.0 * scale, 12.0 * scale])
    decode = summary([12.0 * scale, 13.0 * scale, 14.0 * scale])
    encode.update({"input_GB_per_s": 1.0, "parity_output_GB_per_s": 0.5})
    decode.update({"offered_received_GB_per_s": 1.0,
                   "repaired_output_GB_per_s": 0.5})
    amortized = decode["median_us_per_batch_call"] + \
        plan["median_us"] / CAMPAIGN["reuse"]
    return {
        "schema": "leopard2-benchmark-v2",
        "build": {"compiler": "fixture"},
        "parameters": parameters,
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8", "backend": "avx2",
            "thread_count": 1, "parent_count": 16, "padded_side": 4,
        },
        "correctness": {"leopard2_round_trip": True, "legacy_comparison": None},
        "workload_digests": copy.deepcopy(DIGESTS),
        "metrics": {
            "codec_setup": codec,
            "encode_execution": encode,
            "decode_plan_setup": plan,
            "decode_execution": decode,
            "decode_amortized_at_reuse": {
                "reuse_count": CAMPAIGN["reuse"],
                "derived_median_us_per_batch_call": amortized,
                "offered_received_GB_per_s": 1.0,
                "repaired_output_GB_per_s": 0.5,
            },
            "rate_semantics": "fixture",
        },
        "legacy": {
            "available": False,
            "unavailable_reason": "disabled by --skip-legacy",
            "codec_setup": None,
            "decode_timing_includes_setup": True,
            "encode_execution": None,
            "decode_including_setup": None,
        },
    }


def synthetic_raw(candidate_scale: float = 0.8) -> dict:
    identity = {"fixture": {"sha256": "a" * 64}}
    invocations = []
    for round_index in range(runner.ROUNDS):
        for slot, implementation in enumerate(runner.ORDER):
            result = (baseline_result() if implementation == "baseline"
                      else candidate_result(candidate_scale))
            normalized = runner.validate_result(
                implementation, result, CELL, CAMPAIGN)
            invocations.append({
                "cell_id": CELL.identifier,
                "round": round_index,
                "slot": slot,
                "implementation": implementation,
                "command": [
                    SPECIFICATION["taskset"], "-c", "0",
                    *runner.benchmark_arguments(
                        implementation,
                        Path(SPECIFICATION[f"{implementation}_executable"]),
                        CELL, CAMPAIGN),
                ],
                "environment": copy.deepcopy(runner.CHILD_ENVIRONMENT),
                "pinned_cpu": 0,
                "started_utc": "2026-07-16T00:00:00Z",
                "duration_ns": 1,
                "returncode": 0,
                "stdout": {"path": "unused", "size": 0, "sha256": "0" * 64},
                "stderr": {"path": "unused", "size": 0, "sha256": "0" * 64},
                "result": result,
                "normalized": normalized,
                "identity_before": identity,
                "identity_after": identity,
                "reservation_before": RESERVATION,
                "reservation_after": RESERVATION,
            })
    analysis = runner.analyze(invocations, CAMPAIGN)
    return runner.signed({
        "schema": runner.RAW_SCHEMA,
        "created_utc": "2026-07-16T00:00:00Z",
        "validity_is_independent_of_speed": True,
        "campaign": copy.deepcopy(CAMPAIGN),
        "host_initial": copy.deepcopy(HOST),
        "isolation": copy.deepcopy(ISOLATION),
        "reservation": RESERVATION,
        "input_specification": copy.deepcopy(SPECIFICATION),
        "identities_initial": identity,
        "invocations": invocations,
        "identities_final": identity,
        "host_final": copy.deepcopy(HOST),
        "analysis": analysis,
    })


def resign(value: dict) -> dict:
    payload = copy.deepcopy(value)
    payload.pop("digest", None)
    return runner.signed(payload)


class MainCompareRunnerTests(unittest.TestCase):
    def assert_rejected(self, value: dict) -> None:
        with self.assertRaises(runner.EvidenceError):
            runner.validate_raw(
                resign(value), None, check_files=False, check_current_inputs=False)

    def test_valid_fixture_and_analysis(self) -> None:
        value = synthetic_raw()
        analysis = runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)
        for metric in ("encode", "decode_first_use", "decode_reuse_amortized"):
            result = analysis[CELL.identifier][metric]
            self.assertEqual(result["independent_round_count"], 3)
            self.assertEqual(result["constituent_pair_count"], 6)
            self.assertEqual(result["degrees_of_freedom"], 2)
            self.assertTrue(math.isfinite(result["ci95_lower"]))
            self.assertTrue(result["performance_result_does_not_affect_evidence_validity"])

    def test_legacy_v1_raw_fixture_remains_replayable(self) -> None:
        value = synthetic_raw()
        value["schema"] = runner.RAW_SCHEMA_V1
        value.pop("isolation")
        value = resign(value)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_slower_candidate_is_valid_evidence(self) -> None:
        value = synthetic_raw(candidate_scale=2.0)
        analysis = runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)
        result = analysis[CELL.identifier]["encode"]
        self.assertLess(result["geometric_speedup"], 1.0)
        self.assertFalse(result["faster_ci_excludes_one"])

    def test_parameter_and_missing_mutations_rejected(self) -> None:
        for name, replacement in (
            ("K", 9),
            ("missing_original_indices", [0, 1]),
            ("requested_profile", "low_v1"),
        ):
            with self.subTest(name=name):
                value = synthetic_raw()
                value["invocations"][1]["result"]["parameters"][name] = replacement
                self.assert_rejected(value)

    def test_every_fnv_digest_mutation_rejected(self) -> None:
        for name in ("original_data", "transmitted_parity", "recovered_originals"):
            with self.subTest(name=name):
                value = synthetic_raw()
                value["invocations"][1]["result"]["workload_digests"][name] = "2" * 16
                self.assert_rejected(value)

    def test_sample_mutations_rejected(self) -> None:
        mutations = (
            lambda samples: samples.pop(),
            lambda samples: samples.__setitem__(0, 0.0),
            lambda samples: samples.__setitem__(0, float("nan")),
        )
        for mutation in mutations:
            value = synthetic_raw()
            samples = value["invocations"][0]["result"]["metrics"] \
                ["encode_execution"]["samples_us_per_batch_call"]
            mutation(samples)
            self.assert_rejected(value)

    def test_order_identity_and_analysis_mutations_rejected(self) -> None:
        value = synthetic_raw()
        value["invocations"][0]["implementation"] = "candidate"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["invocations"][0]["identity_after"] = {"fixture": {"sha256": "b" * 64}}
        self.assert_rejected(value)
        value = synthetic_raw()
        value["analysis"][CELL.identifier]["encode"]["geometric_speedup"] += 1.0
        self.assert_rejected(value)

    def test_environment_host_and_build_spec_mutations_rejected(self) -> None:
        value = synthetic_raw()
        value["invocations"][0]["environment"]["LD_PRELOAD"] = "/tmp/injected.so"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["campaign"]["child_environment"]["OMP_NUM_THREADS"] = "2"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["host_final"]["benchmark_cpu"]["frequency_policy"] \
            ["scaling_governor"] = "powersave"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["input_specification"].pop("candidate_build_dir")
        self.assert_rejected(value)
        value = synthetic_raw()
        value["invocations"][5]["result"]["resolved"]["backend"] = "scalar"
        value["invocations"][5]["normalized"]["backend"] = "scalar"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["reservation"]["sha256"] = "f" * 64
        for invocation in value["invocations"]:
            invocation["reservation_before"]["sha256"] = "f" * 64
            invocation["reservation_after"]["sha256"] = "f" * 64
        self.assert_rejected(value)

    def test_isolation_mutations_and_active_sibling_rejected(self) -> None:
        value = synthetic_raw()
        value["isolation"]["delta"]["reserved_sibling"]["idle_jiffies"] += 1
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["pair_lease"]["path"] = "/tmp/wrong-pair.lock"
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["accepted"] = False
        self.assert_rejected(value)

        value = synthetic_raw()
        before = value["isolation"]["before"]
        after = value["isolation"]["after"]
        active_after = copy.deepcopy(after["reserved_sibling"])
        active_after["fields"]["user"] += 1
        active_after["total_jiffies"] += 1
        value["isolation"] = runner.isolation_record(
            0, 1, value["isolation"]["pair_lease"],
            before["monotonic_ns"], after["monotonic_ns"],
            before["benchmark_cpu"], after["benchmark_cpu"],
            before["reserved_sibling"], active_after)
        self.assertFalse(value["isolation"]["accepted"])
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["before"]["benchmark_cpu"]["cpu"] = 98
        value["isolation"]["after"]["benchmark_cpu"]["cpu"] = 98
        value["isolation"]["delta"]["benchmark_cpu"]["cpu"] = 98
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["before"]["benchmark_cpu"] = []
        self.assert_rejected(value)

        value = synthetic_raw()
        before = value["isolation"]["before"]
        after = value["isolation"]["after"]
        before["benchmark_cpu"], before["reserved_sibling"] = (
            before["reserved_sibling"], before["benchmark_cpu"])
        after["benchmark_cpu"], after["reserved_sibling"] = (
            after["reserved_sibling"], after["benchmark_cpu"])
        value["isolation"]["delta"]["benchmark_cpu"], \
            value["isolation"]["delta"]["reserved_sibling"] = (
                value["isolation"]["delta"]["reserved_sibling"],
                value["isolation"]["delta"]["benchmark_cpu"])
        self.assert_rejected(value)

        value = synthetic_raw()
        value["invocations"][0]["duration_ns"] = 0
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["after"]["monotonic_ns"] = 1_001
        self.assert_rejected(value)

    def test_isolation_replay_does_not_require_original_uid(self) -> None:
        value = synthetic_raw()
        payload = runner.pair_lease_payload(0, 1, uid=12345)
        lease = {
            "device": 10,
            "directory_device": 11,
            "directory_inode": 12,
            "inode": 13,
            "lock": "exclusive_nonblocking_pair_wide",
            "path": str(runner.pair_lease_directory(12345) /
                runner.pair_lease_name(0, 1, uid=12345)),
            "payload": payload,
            "sha256": runner.sha256_bytes(runner.canonical_bytes(payload)),
        }
        before = value["isolation"]["before"]
        after = value["isolation"]["after"]
        value["isolation"] = runner.isolation_record(
            0, 1, lease, before["monotonic_ns"], after["monotonic_ns"],
            before["benchmark_cpu"], after["benchmark_cpu"],
            before["reserved_sibling"], after["reserved_sibling"])
        value = resign(value)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_signature_mutation_rejected(self) -> None:
        value = synthetic_raw()
        value["created_utc"] = "changed"
        with self.assertRaises(runner.EvidenceError):
            runner.validate_raw(
                value, None, check_files=False, check_current_inputs=False)

    def test_retained_stdout_mutation_rejected(self) -> None:
        value = synthetic_raw()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for index, invocation in enumerate(value["invocations"]):
                stdout = json.dumps(invocation["result"]).encode("utf-8")
                stderr = b""
                stdout_path = root / f"{index}.stdout"
                stderr_path = root / f"{index}.stderr"
                stdout_path.write_bytes(stdout)
                stderr_path.write_bytes(stderr)
                invocation["stdout"] = {
                    "path": stdout_path.name,
                    "size": len(stdout),
                    "sha256": runner.sha256_bytes(stdout),
                }
                invocation["stderr"] = {
                    "path": stderr_path.name,
                    "size": 0,
                    "sha256": runner.sha256_bytes(stderr),
                }
            value = resign(value)
            runner.validate_raw(
                value, root, check_files=True, check_current_inputs=False)
            (root / "0.stdout").write_bytes(b"{}")
            with self.assertRaises(runner.EvidenceError):
                runner.validate_raw(
                    value, root, check_files=True, check_current_inputs=False)

    def test_v2_manifest_binds_isolation_and_replays_portably(self) -> None:
        value = synthetic_raw()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for index, invocation in enumerate(value["invocations"]):
                stdout = json.dumps(invocation["result"]).encode("utf-8")
                stdout_path = root / f"{index}.stdout"
                stderr_path = root / f"{index}.stderr"
                stdout_path.write_bytes(stdout)
                stderr_path.write_bytes(b"")
                invocation["stdout"] = {
                    "path": stdout_path.name,
                    "size": len(stdout),
                    "sha256": runner.sha256_bytes(stdout),
                }
                invocation["stderr"] = {
                    "path": stderr_path.name,
                    "size": 0,
                    "sha256": runner.sha256_bytes(b""),
                }
            value = resign(value)
            raw_path = root / "raw.json"
            runner.write_json_exclusive(raw_path, value)
            manifest = runner.signed({
                "schema": runner.MANIFEST_SCHEMA,
                "created_utc": "2026-07-16T00:00:00Z",
                "valid": True,
                "validity_is_independent_of_speed": True,
                "raw": {
                    "path": raw_path.name,
                    "size": raw_path.stat().st_size,
                    "sha256": runner.sha256_file(raw_path),
                    "payload_digest": value["digest"],
                },
                "campaign": value["campaign"],
                "host": value["host_initial"],
                "isolation": value["isolation"],
                "reservation": value["reservation"],
                "identities": value["identities_initial"],
                "analysis": value["analysis"],
            })
            manifest_path = root / "manifest.json"
            runner.write_json_exclusive(manifest_path, manifest)
            options = argparse.Namespace(
                manifest=manifest_path, no_current_input_check=True)
            self.assertEqual(runner.verify_campaign(options), 0)

            edited = copy.deepcopy(manifest)
            edited["isolation"]["accepted"] = False
            edited = resign(edited)
            manifest_path.unlink()
            runner.write_json_exclusive(manifest_path, edited)
            with self.assertRaises(runner.EvidenceError):
                runner.verify_campaign(options)

    def test_failed_bundle_binds_retained_invocation_files(self) -> None:
        value = synthetic_raw()
        invocation = copy.deepcopy(value["invocations"][0])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stdout = root / "first.stdout"
            stderr = root / "first.stderr"
            stdout.write_text(json.dumps(invocation["result"]), encoding="utf-8")
            stderr.write_bytes(b"diagnostic")
            invocation["stdout"] = {
                "path": stdout.name, "size": stdout.stat().st_size,
                "sha256": runner.sha256_file(stdout),
            }
            invocation["stderr"] = {
                "path": stderr.name, "size": stderr.stat().st_size,
                "sha256": runner.sha256_file(stderr),
            }
            failure = runner.signed({
                "schema": runner.FAILURE_SCHEMA,
                "created_utc": "2026-07-16T00:00:00Z",
                "status": "failed", "valid": False,
                "error_type": "EvidenceError", "error": "fixture failure",
                "campaign": copy.deepcopy(CAMPAIGN),
                "host_initial": copy.deepcopy(HOST),
                "reservation": copy.deepcopy(RESERVATION),
                "pair_lease": copy.deepcopy(PAIR_LEASE),
                "isolation": copy.deepcopy(ISOLATION),
                "input_specification": copy.deepcopy(SPECIFICATION),
                "identities_initial": copy.deepcopy(invocation["identity_before"]),
                "invocations": [invocation],
                "retained_files": runner.retained_file_records(root),
                "traceback": "fixture traceback",
            })
            failure_path = root / "failure.json"
            runner.write_json_exclusive(failure_path, failure)
            runner.validate_failure(failure, root, check_files=True)
            self.assertEqual(runner.verify_failed_campaign(
                argparse.Namespace(failure=failure_path)), 0)

            stdout.write_bytes(b"{}")
            semantic_mismatch = copy.deepcopy(failure)
            replacement = {
                "path": stdout.name, "size": stdout.stat().st_size,
                "sha256": runner.sha256_file(stdout),
            }
            semantic_mismatch["invocations"][0]["stdout"] = replacement
            for index, retained in enumerate(semantic_mismatch["retained_files"]):
                if retained["path"] == stdout.name:
                    semantic_mismatch["retained_files"][index] = replacement
            semantic_mismatch = resign(semantic_mismatch)
            with self.assertRaises(runner.EvidenceError):
                runner.validate_failure(semantic_mismatch, root, check_files=True)

            stdout.write_bytes(b"edited")
            with self.assertRaises(runner.EvidenceError):
                runner.validate_failure(failure, root, check_files=True)

    def test_reservation_is_locked_and_canonical(self) -> None:
        payload = {
            "benchmark_cpu": 0,
            "nonce": "fixture-nonce",
            "owner": "unit test",
            "reserved_sibling": 1,
            "schema": runner.RESERVATION_SCHEMA,
            "status": "held",
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reservation.json"
            path.write_bytes(runner.canonical_bytes(payload))
            with runner.Reservation(path, 0, 1) as identity:
                self.assertEqual(identity["payload"], payload)
                with self.assertRaises(runner.EvidenceError):
                    with runner.Reservation(path, 0, 1):
                        pass
            path.write_bytes(runner.canonical_bytes(payload) + b"\n")
            with self.assertRaises(runner.EvidenceError):
                with runner.Reservation(path, 0, 1):
                    pass

    def test_pair_lease_serializes_different_reservation_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            with runner.PairLease(0, 1, root=root) as identity:
                self.assertEqual(identity["payload"], runner.pair_lease_payload(0, 1))
                with self.assertRaises(runner.EvidenceError):
                    with runner.PairLease(1, 0, root=root):
                        pass
            path = root / runner.pair_lease_name(0, 1)
            path.chmod(0o644)
            with self.assertRaises(runner.EvidenceError):
                with runner.PairLease(0, 1, root=root):
                    pass

    def test_pair_lease_detects_unlink_and_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            first = runner.PairLease(0, 1, root=root)
            with first:
                first.path.unlink()
                with self.assertRaises(runner.EvidenceError):
                    with runner.PairLease(0, 1, root=root):
                        pass
                with self.assertRaises(runner.EvidenceError):
                    first.validate_current()

    def test_pair_lease_interoperates_with_jerasure_runner(self) -> None:
        jerasure_path = MODULE_PATH.resolve().parents[3] / \
            "tools/leopard2_jerasure_compare.py"
        specification = importlib.util.spec_from_file_location(
            "jerasure_compare_pair_lease_test", jerasure_path)
        self.assertIsNotNone(specification)
        self.assertIsNotNone(specification.loader)
        jerasure = importlib.util.module_from_spec(specification)
        sys.modules[specification.name] = jerasure
        specification.loader.exec_module(jerasure)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lease_directory = root / jerasure.pair_lease_directory_name()
            with runner.PairLease(0, 1, root=lease_directory):
                with self.assertRaises(jerasure.ComparisonError):
                    with jerasure.PairLease(1, 0, root=root):
                        pass
            with jerasure.PairLease(0, 1, root=root):
                with self.assertRaises(runner.EvidenceError):
                    with runner.PairLease(1, 0, root=lease_directory):
                        pass

    def test_pair_lease_creation_ignores_restrictive_umask(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            previous = os.umask(0o777)
            try:
                with runner.PairLease(0, 1, root=root) as identity:
                    self.assertEqual(
                        stat.S_IMODE(os.lstat(Path(identity["path"])).st_mode),
                        0o600)
            finally:
                os.umask(previous)
    def test_stable_anchor_blocks_recreated_reservation_inode(self) -> None:
        payload = {
            "benchmark_cpu": 0,
            "nonce": "stable-anchor-fixture",
            "owner": "unit test",
            "reserved_sibling": 1,
            "schema": runner.RESERVATION_SCHEMA,
            "status": "held",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            runtime_root = root / "runtime"
            runtime_root.mkdir(mode=0o700)
            runtime_root.chmod(0o700)
            path = root / "reservation.json"
            old = root / "reservation.old"
            path.write_bytes(runner.canonical_bytes(payload))
            with runner.Reservation(
                    path, 0, 1, runtime_root=runtime_root):
                path.rename(old)
                path.write_bytes(runner.canonical_bytes(payload))
                try:
                    with self.assertRaises(runner.EvidenceError):
                        with runner.Reservation(
                                path, 0, 1, runtime_root=runtime_root):
                            pass
                finally:
                    path.unlink()
                    old.rename(path)

    def test_custom_cell_parser_rejects_non_wire_cases(self) -> None:
        good = runner.parse_cell("cell:240:16:65536:8:1")
        self.assertEqual((good.k, good.r, good.losses), (240, 16, 8))
        for bad in (
            "bad:8:9:64:1:1",
            "bad:8:4:65:1:1",
            "bad:8:4:64:5:1",
            "missing-fields",
        ):
            with self.subTest(cell=bad), self.assertRaises(runner.EvidenceError):
                runner.parse_cell(bad)


if __name__ == "__main__":
    unittest.main(verbosity=2)
