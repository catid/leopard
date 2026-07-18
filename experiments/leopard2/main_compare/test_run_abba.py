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
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock
from dataclasses import asdict
from pathlib import Path
from typing import Mapping


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
    "candidate_archive": "/fixture/candidate-build/libleopard.a",
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
    "candidate_mode": "auto",
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


def candidate_result(
    scale: float = 0.8, raw_schema: str = runner.RAW_SCHEMA,
    campaign: Mapping[str, object] = CAMPAIGN,
) -> dict:
    parameters = common_parameters()
    parameters.update({
        "requested_profile": "legacy_high_v1",
        "requested_field": "auto",
        "requested_backend": "auto",
        "skip_legacy": True,
        "retain_samples": True,
    })
    flags = runner.candidate_mode_flags(
        runner.candidate_mode_for_campaign(campaign))
    parameters.update({
        "force_generic_decode": flags["force_generic_decode"],
        "force_specialized_decode": flags["force_specialized_decode"],
    })
    if raw_schema in (runner.RAW_SCHEMA_V3, runner.RAW_SCHEMA):
        parameters.update({name: flags[name] for name in (
            "force_tiled_decode", "force_materialized_decode")})
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


def archive_recipe_fixture_text(cmake: dict | Mapping[str, str]) -> str:
    return (
        f"/usr/bin/ar qc {cmake['archive']} "
        f"CMakeFiles/{cmake['target_directory']}/LeopardCommon.cpp.o "
        "CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2.cpp.o\n"
        f"/usr/bin/ranlib {cmake['archive']}\n"
    )


def cmake_fixture_identity(raw_schema: str) -> tuple[dict, dict]:
    cmake = runner.cmake_identity_for_raw_schema(raw_schema)
    archive_path = f"/fixture/candidate-build/{cmake['archive']}"
    archive = {"path": archive_path, "sha256": "a" * 64}
    recipe_text = archive_recipe_fixture_text(cmake)
    recipe_content = runner.exact_text_content(
        recipe_text, "fixture archive link recipe")
    recipe = {
        "path": "/fixture/candidate-build/CMakeFiles/"
                f"{cmake['target_directory']}/link.txt",
        "size": recipe_content["size"],
        "sha256": recipe_content["sha256"],
    }
    build = {
        "build_dir": "/fixture/candidate-build",
        "validated_archive": copy.deepcopy(archive),
        "archive_link_recipe": recipe,
        "archiver": {"path": "/usr/bin/ar"},
        "validated_archive_members": [
            "LeopardCommon.cpp.o", "Leopard2BackendAVX2.cpp.o"],
        "validated_compile_commands": {
            "required_source_object_pairs": [
                {
                    "source": {"path": "/fixture/source/LeopardCommon.cpp"},
                    "object": {"path": "/fixture/candidate-build/CMakeFiles/"
                               f"{cmake['target_directory']}/LeopardCommon.cpp.o"},
                },
                {
                    "source": {"path": "/fixture/source/Leopard2BackendAVX2.cpp"},
                    "object": {"path": "/fixture/candidate-build/CMakeFiles/"
                               "leopard2_backend_avx2.dir/"
                               "Leopard2BackendAVX2.cpp.o"},
                },
            ],
        },
    }
    if raw_schema in runner.HARDENED_BUILD_SCHEMAS:
        build["archive_link_recipe_content"] = recipe_content
        build["ranlib"] = {"path": "/usr/bin/ranlib"}
        build["archive_link_tool_invocations"] = {
            "archiver": {"invocation": "/usr/bin/ar",
                         "resolved_path": "/usr/bin/ar"},
            "ranlib": {"invocation": "/usr/bin/ranlib",
                       "resolved_path": "/usr/bin/ranlib"},
        }
    identity = {
        "fixture": {"sha256": "a" * 64},
        "candidate_archive": copy.deepcopy(archive),
        "candidate_build": build,
    }
    specification = copy.deepcopy(SPECIFICATION)
    specification["candidate_archive"] = archive_path
    return identity, specification


def synthetic_raw(
    candidate_scale: float = 0.8, raw_schema: str = runner.RAW_SCHEMA,
    candidate_mode: str = "auto",
) -> dict:
    identity, specification = cmake_fixture_identity(raw_schema)
    campaign = copy.deepcopy(CAMPAIGN)
    if raw_schema == runner.RAW_SCHEMA:
        campaign["candidate_mode"] = candidate_mode
    else:
        campaign.pop("candidate_mode", None)
    invocations = []
    for round_index in range(runner.ROUNDS):
        for slot, implementation in enumerate(runner.ORDER):
            result = (baseline_result() if implementation == "baseline"
                      else candidate_result(
                          candidate_scale, raw_schema, campaign))
            normalized = runner.validate_result(
                implementation, result, CELL, campaign, raw_schema)
            invocations.append({
                "cell_id": CELL.identifier,
                "round": round_index,
                "slot": slot,
                "implementation": implementation,
                "command": [
                    specification["taskset"], "-c", "0",
                    *runner.benchmark_arguments(
                        implementation,
                        Path(specification[f"{implementation}_executable"]),
                        CELL, campaign),
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
    analysis = runner.analyze(invocations, campaign)
    return runner.signed({
        "schema": raw_schema,
        "created_utc": "2026-07-16T00:00:00Z",
        "validity_is_independent_of_speed": True,
        "campaign": campaign,
        "host_initial": copy.deepcopy(HOST),
        "isolation": copy.deepcopy(ISOLATION),
        "reservation": RESERVATION,
        "input_specification": specification,
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


def recursively_replace_strings(value: object, replacements: tuple[tuple[str, str], ...]
                                ) -> object:
    if isinstance(value, str):
        result = value
        for before, after in replacements:
            result = result.replace(before, after)
        return result
    if isinstance(value, list):
        return [recursively_replace_strings(item, replacements) for item in value]
    if isinstance(value, dict):
        return {
            key: recursively_replace_strings(item, replacements)
            for key, item in value.items()
        }
    return value


def attach_recipe_content(value: object, content: dict) -> None:
    if isinstance(value, list):
        for item in value:
            attach_recipe_content(item, content)
        return
    if not isinstance(value, dict):
        return
    build = value.get("candidate_build")
    if isinstance(build, dict) and isinstance(build.get("archive_link_recipe"), dict):
        build["archive_link_recipe"]["size"] = content["size"]
        build["archive_link_recipe"]["sha256"] = content["sha256"]
        build["archive_link_recipe_content"] = copy.deepcopy(content)
        build["ranlib"] = {"path": "/usr/bin/ranlib"}
        build["archive_link_tool_invocations"] = {
            "archiver": {"invocation": "/usr/bin/ar",
                         "resolved_path": "/usr/bin/ar"},
            "ranlib": {"invocation": "/usr/bin/ranlib",
                       "resolved_path": "/usr/bin/ranlib"},
        }
    for item in value.values():
        attach_recipe_content(item, content)


def replace_current_recipe_text(value: dict, text: str) -> None:
    content = runner.exact_text_content(text, "mutated archive link recipe")
    build = value["identities_initial"]["candidate_build"]
    build["archive_link_recipe_content"] = content
    build["archive_link_recipe"]["size"] = content["size"]
    build["archive_link_recipe"]["sha256"] = content["sha256"]
    value["identities_final"] = copy.deepcopy(value["identities_initial"])
    for invocation in value["invocations"]:
        invocation["identity_before"] = copy.deepcopy(value["identities_initial"])
        invocation["identity_after"] = copy.deepcopy(value["identities_initial"])


def synthetic_failure(raw_schema: str) -> dict:
    raw = synthetic_raw(raw_schema=raw_schema)
    failure_schema = {
        runner.RAW_SCHEMA_V2: runner.FAILURE_SCHEMA_V2,
        runner.RAW_SCHEMA_V3: runner.FAILURE_SCHEMA_V3,
        runner.RAW_SCHEMA: runner.FAILURE_SCHEMA,
    }[raw_schema]
    return runner.signed({
        "schema": failure_schema,
        "created_utc": "2026-07-16T00:00:00Z",
        "status": "failed",
        "valid": False,
        "error_type": "EvidenceError",
        "error": "fixture failure",
        "campaign": copy.deepcopy(raw["campaign"]),
        "host_initial": copy.deepcopy(raw["host_initial"]),
        "reservation": copy.deepcopy(raw["reservation"]),
        "pair_lease": copy.deepcopy(PAIR_LEASE),
        "isolation": copy.deepcopy(raw["isolation"]),
        "input_specification": copy.deepcopy(raw["input_specification"]),
        "identities_initial": copy.deepcopy(raw["identities_initial"]),
        "invocations": [],
        "retained_files": [],
        "traceback": "fixture traceback",
    })


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
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V1)
        value.pop("isolation")
        value = resign(value)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_legacy_v2_raw_fixture_remains_replayable(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V2)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_legacy_v3_raw_fixture_remains_replayable(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V3)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_candidate_modes_bind_flags_and_exact_argv(self) -> None:
        expected_arguments = {
            "auto": (),
            "generic": ("--force-generic",),
            "materialized": ("--force-specialized", "--force-materialized"),
            "tiled": ("--force-specialized", "--force-tiled"),
        }
        for mode in runner.CANDIDATE_MODES:
            with self.subTest(mode=mode):
                value = synthetic_raw(candidate_mode=mode)
                runner.validate_raw(
                    value, None, check_files=False,
                    check_current_inputs=False)
                flags = runner.candidate_mode_flags(mode)
                invocation = next(item for item in value["invocations"]
                                  if item["implementation"] == "candidate")
                parameters = invocation["result"]["parameters"]
                for name, expected in flags.items():
                    self.assertIs(parameters[name], expected)
                command = invocation["command"]
                for argument in expected_arguments[mode]:
                    self.assertIn(argument, command)

        value = synthetic_raw(candidate_mode="generic")
        value["campaign"]["candidate_mode"] = "auto"
        self.assert_rejected(value)

        value = synthetic_raw(candidate_mode="generic")
        invocation = next(item for item in value["invocations"]
                          if item["implementation"] == "candidate")
        invocation["command"].remove("--force-generic")
        self.assert_rejected(value)

        value = synthetic_raw(candidate_mode="tiled")
        invocation = next(item for item in value["invocations"]
                          if item["implementation"] == "candidate")
        invocation["result"]["parameters"]["force_tiled_decode"] = False
        self.assert_rejected(value)

        value = synthetic_raw()
        value["campaign"].pop("candidate_mode")
        self.assert_rejected(value)

        value = synthetic_raw()
        value["campaign"]["candidate_mode"] = ["auto"]
        self.assert_rejected(value)

        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V3)
        value["campaign"]["candidate_mode"] = "auto"
        self.assert_rejected(value)

    def test_workspace_selector_schema_boundary_is_fail_closed(self) -> None:
        for name in ("force_tiled_decode", "force_materialized_decode"):
            value = synthetic_raw()
            for invocation in value["invocations"]:
                if invocation["implementation"] == "candidate":
                    invocation["result"]["parameters"].pop(name)
            self.assert_rejected(value)

            value = synthetic_raw()
            for invocation in value["invocations"]:
                if invocation["implementation"] == "candidate":
                    invocation["result"]["parameters"][name] = True
            self.assert_rejected(value)

        for raw_schema in (runner.RAW_SCHEMA_V1, runner.RAW_SCHEMA_V2):
            value = synthetic_raw(raw_schema=raw_schema)
            for invocation in value["invocations"]:
                if invocation["implementation"] == "candidate":
                    invocation["result"]["parameters"].update({
                        "force_tiled_decode": False,
                        "force_materialized_decode": False,
                    })
            self.assert_rejected(value)

    def test_cmake_identity_and_cross_schema_relabels_are_rejected(self) -> None:
        value = synthetic_raw()
        value["input_specification"]["candidate_archive"] = \
            "/fixture/candidate-build/liblibleopard.a"
        self.assert_rejected(value)

        value = synthetic_raw()
        value["identities_initial"]["candidate_build"][
            "archive_link_recipe"]["path"] = \
            "/fixture/candidate-build/CMakeFiles/libleopard.dir/link.txt"
        value["identities_final"] = copy.deepcopy(value["identities_initial"])
        for invocation in value["invocations"]:
            invocation["identity_before"] = copy.deepcopy(
                value["identities_initial"])
            invocation["identity_after"] = copy.deepcopy(
                value["identities_initial"])
        self.assert_rejected(value)

        value = synthetic_raw()
        value["schema"] = runner.RAW_SCHEMA_V2
        self.assert_rejected(value)

        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V2)
        value["schema"] = runner.RAW_SCHEMA
        self.assert_rejected(value)

    def test_private_hardened_historical_build_schema_is_not_evidence(self) -> None:
        private = runner.HARDENED_HISTORICAL_BUILD_SCHEMA
        self.assertIn(private, runner.BUILD_SCHEMA_TO_CMAKE_IDENTITY)
        self.assertIn(private, runner.HARDENED_BUILD_SCHEMAS)
        self.assertNotIn(private, runner.RAW_TO_CMAKE_IDENTITY)
        self.assertNotIn(private, runner.MANIFEST_TO_RAW_SCHEMA.values())
        self.assertNotIn(private, runner.FAILURE_TO_RAW_SCHEMA.values())
        self.assertEqual(
            runner.HISTORICAL_CMAKE_IDENTITY,
            runner.cmake_identity_for_build_schema(private))
        with self.assertRaises(runner.EvidenceError):
            runner.cmake_identity_for_raw_schema(private)

        value = synthetic_raw()
        value["schema"] = private
        self.assert_rejected(value)

    def test_coherent_historical_recipe_relabel_is_rejected(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V2)
        value = recursively_replace_strings(value, (
            ("liblibleopard.a", "libleopard.a"),
            ("CMakeFiles/libleopard.dir", "CMakeFiles/leopard.dir"),
        ))
        self.assertIsInstance(value, dict)
        value["schema"] = runner.RAW_SCHEMA
        historical = runner.cmake_identity_for_raw_schema(runner.RAW_SCHEMA_V2)
        old_content = runner.exact_text_content(
            archive_recipe_fixture_text(historical),
            "historical fixture archive link recipe")
        attach_recipe_content(value, old_content)
        # Every path and every surrounding digest can be coherently relabeled,
        # but the retained old recipe bytes still describe the old CMake target.
        self.assert_rejected(value)

    def test_recipe_content_and_identity_mutations_are_rejected(self) -> None:
        value = synthetic_raw()
        value["identities_initial"]["candidate_build"][
            "archive_link_recipe_content"]["text"] += "\n"
        value["identities_final"] = copy.deepcopy(value["identities_initial"])
        for invocation in value["invocations"]:
            invocation["identity_before"] = copy.deepcopy(
                value["identities_initial"])
            invocation["identity_after"] = copy.deepcopy(
                value["identities_initial"])
        self.assert_rejected(value)

    def test_recipe_command_semantic_mutations_are_rejected(self) -> None:
        canonical = archive_recipe_fixture_text(
            runner.cmake_identity_for_raw_schema(runner.RAW_SCHEMA))
        mutations = {
            "noncanonical output path": canonical.replace(
                "ar qc libleopard.a", "ar qc nested/libleopard.a", 1),
            "different archiver": canonical.replace(
                "/usr/bin/ar", "/tmp/ar", 1),
            "response file": canonical.replace(
                "LeopardCommon.cpp.o ", "LeopardCommon.cpp.o @objects.rsp ", 1),
            "object order": canonical.replace(
                "CMakeFiles/leopard.dir/LeopardCommon.cpp.o "
                "CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2.cpp.o",
                "CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2.cpp.o "
                "CMakeFiles/leopard.dir/LeopardCommon.cpp.o"),
            "different ranlib": canonical.replace(
                "/usr/bin/ranlib", "/tmp/ranlib", 1),
        }
        for label, recipe in mutations.items():
            with self.subTest(label=label):
                value = synthetic_raw()
                replace_current_recipe_text(value, recipe)
                self.assert_rejected(value)

        value = synthetic_raw()
        value["identities_initial"]["candidate_build"][
            "archive_link_recipe"]["sha256"] = "b" * 64
        value["identities_final"] = copy.deepcopy(value["identities_initial"])
        for invocation in value["invocations"]:
            invocation["identity_before"] = copy.deepcopy(
                value["identities_initial"])
            invocation["identity_after"] = copy.deepcopy(
                value["identities_initial"])
        self.assert_rejected(value)

    def test_coherent_failed_historical_recipe_relabel_is_rejected(self) -> None:
        historical = synthetic_failure(runner.RAW_SCHEMA_V2)
        runner.validate_failure(historical, Path("/unused"), check_files=False)
        current = synthetic_failure(runner.RAW_SCHEMA)
        runner.validate_failure(current, Path("/unused"), check_files=False)

        relabeled = recursively_replace_strings(historical, (
            ("liblibleopard.a", "libleopard.a"),
            ("CMakeFiles/libleopard.dir", "CMakeFiles/leopard.dir"),
        ))
        self.assertIsInstance(relabeled, dict)
        relabeled["schema"] = runner.FAILURE_SCHEMA
        old_identity = runner.cmake_identity_for_raw_schema(runner.RAW_SCHEMA_V2)
        attach_recipe_content(relabeled, runner.exact_text_content(
            archive_recipe_fixture_text(old_identity),
            "historical failed fixture archive link recipe"))
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(relabeled), Path("/unused"), check_files=False)

    def test_non_string_schema_values_fail_closed(self) -> None:
        value = synthetic_raw()
        value["schema"] = {"unexpected": "object"}
        self.assert_rejected(value)
        failure = synthetic_failure(runner.RAW_SCHEMA)
        failure["schema"] = [runner.FAILURE_SCHEMA]
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(failure), Path("/unused"), check_files=False)

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

    def test_v4_manifest_binds_mode_isolation_and_replays_portably(self) -> None:
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

            cross_schema = copy.deepcopy(manifest)
            cross_schema["schema"] = runner.MANIFEST_SCHEMA_V2
            cross_schema = resign(cross_schema)
            manifest_path.unlink()
            runner.write_json_exclusive(manifest_path, cross_schema)
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

            legacy = copy.deepcopy(failure)
            legacy["schema"] = runner.FAILURE_SCHEMA_V2
            legacy["campaign"].pop("candidate_mode", None)
            old_identity, old_specification = cmake_fixture_identity(
                runner.RAW_SCHEMA_V2)
            legacy["input_specification"] = old_specification
            legacy["identities_initial"] = old_identity
            legacy["invocations"] = []
            legacy["retained_files"] = runner.retained_file_records(root)
            legacy = resign(legacy)
            runner.validate_failure(legacy, root, check_files=True)

    def test_bounded_file_snapshot_rejects_fifo_without_open_block(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fifo = Path(directory) / "identity.fifo"
            os.mkfifo(fifo, 0o600)
            with self.assertRaises(runner.EvidenceError):
                runner.bounded_file_snapshot(fifo)

    def test_process_group_reap_never_uses_unbounded_wait(self) -> None:
        class NeverReaps:
            pid = 123456
            returncode = None

            def __init__(self) -> None:
                self.calls: list[float | None] = []

            def wait(self, timeout: float | None = None) -> int:
                self.calls.append(timeout)
                raise subprocess.TimeoutExpired(("never",), timeout)

        process = NeverReaps()
        with mock.patch.object(os, "killpg"):
            reaped, returncode = runner.terminate_process_group_bounded(process)
        self.assertFalse(reaped)
        self.assertEqual(returncode, -9)
        self.assertEqual(process.calls, [5.0])
        self.assertNotIn(None, process.calls)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux descendant containment test")
    def test_timeout_kills_setsid_double_fork_descendant(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            marker = root / "escaped-marker"
            ready = root / "escaped-ready"
            subreaper_before = runner._get_child_subreaper()
            child = (
                "import os,sys,time\n"
                "pid=os.fork()\n"
                "if pid == 0:\n"
                " os.setsid()\n"
                " daemon=os.fork()\n"
                " if daemon != 0: os._exit(0)\n"
                " open(sys.argv[2], 'w').write(str(os.getpid()))\n"
                " time.sleep(1.5)\n"
                " open(sys.argv[1], 'w').write('escaped')\n"
                " os._exit(0)\n"
                "deadline=time.monotonic()+5\n"
                "while not os.path.exists(sys.argv[2]) and time.monotonic()<deadline:\n"
                " time.sleep(.01)\n"
                "while True: time.sleep(1)\n"
            )
            with self.assertRaisesRegex(runner.EvidenceError, "exceeded"):
                runner.run_process_bounded(
                    [sys.executable, "-c", child, str(marker), str(ready)],
                    timeout=0.8, max_stdout=1024, max_stderr=1024)
            escaped_pid = int(ready.read_text(encoding="utf-8"))
            self.assertEqual(runner._get_child_subreaper(), subreaper_before)
            time.sleep(1.0)
            self.assertFalse(marker.exists())
            self.assertFalse(Path("/proc", str(escaped_pid)).exists())

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux descendant containment test")
    def test_successful_leader_cannot_leave_daemon_descendant(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            marker = root / "daemon-marker"
            ready = root / "daemon-ready"
            child = (
                "import os,sys,time\n"
                "pid=os.fork()\n"
                "if pid == 0:\n"
                " os.setsid()\n"
                " daemon=os.fork()\n"
                " if daemon != 0: os._exit(0)\n"
                " null=os.open('/dev/null', os.O_WRONLY)\n"
                " os.dup2(null, 1);os.dup2(null, 2);os.close(null)\n"
                " open(sys.argv[2], 'w').write(str(os.getpid()))\n"
                " time.sleep(1.0)\n"
                " open(sys.argv[1], 'w').write('escaped')\n"
                " os._exit(0)\n"
                "deadline=time.monotonic()+5\n"
                "while not os.path.exists(sys.argv[2]) and time.monotonic()<deadline:\n"
                " time.sleep(.01)\n"
            )
            completed = runner.run_process_bounded(
                [sys.executable, "-c", child, str(marker), str(ready)],
                timeout=2.0, max_stdout=1024, max_stderr=1024)
            self.assertEqual(completed.returncode, 0)
            escaped_pid = int(ready.read_text(encoding="utf-8"))
            time.sleep(1.1)
            self.assertFalse(marker.exists())
            self.assertFalse(Path("/proc", str(escaped_pid)).exists())

    def test_descendant_containment_fails_closed_before_spawn(self) -> None:
        with mock.patch.object(runner.sys, "platform", "not-linux"), \
             mock.patch.object(runner.subprocess, "Popen") as popen, \
             self.assertRaisesRegex(runner.EvidenceError, "requires Linux"):
            runner.run_process_bounded(
                [sys.executable, "-c", "pass"], timeout=1.0,
                max_stdout=1024, max_stderr=1024)
        popen.assert_not_called()

        for unavailable in (
                "_validate_linux_pidfd_support", "_get_child_subreaper",
                "_proc_process_snapshot"):
            with self.subTest(unavailable=unavailable), \
                 mock.patch.object(
                     runner, unavailable,
                     side_effect=runner.EvidenceError(unavailable + " unavailable")), \
                 mock.patch.object(runner.subprocess, "Popen") as popen, \
                 self.assertRaisesRegex(runner.EvidenceError, "unavailable"):
                runner.run_process_bounded(
                    [sys.executable, "-c", "pass"], timeout=1.0,
                    max_stdout=1024, max_stderr=1024)
            popen.assert_not_called()

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
