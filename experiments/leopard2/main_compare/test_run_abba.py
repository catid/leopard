#!/usr/bin/env python3
"""Self and mutation tests for the exact-main ABBA evidence runner."""

from __future__ import annotations

import argparse
import base64
import copy
import hashlib
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


def load_plan_runner():
    path = MODULE_PATH.parents[1] / "decoder_dispatch" / \
        "plan_balanced_promotion.py"
    module_name = "main_compare_test_balanced_promotion"
    existing = sys.modules.get(module_name)
    if existing is not None:
        return existing
    module_spec = importlib.util.spec_from_file_location(module_name, path)
    assert module_spec is not None and module_spec.loader is not None
    module = importlib.util.module_from_spec(module_spec)
    decoder_directory = str(path.parent)
    inserted = decoder_directory not in sys.path
    if inserted:
        sys.path.insert(0, decoder_directory)
    try:
        sys.modules[module_name] = module
        module_spec.loader.exec_module(module)
    finally:
        if inserted:
            sys.path.remove(decoder_directory)
    return module


BASELINE_TREE = "b7c8830d96a978f6ec14fe747095f066e351ae72"
BASELINE_COMMIT_BASE64 = (
    "dHJlZSBiN2M4ODMwZDk2YTk3OGY2ZWMxNGZlNzQ3MDk1ZjA2NmUzNTFhZTcyCnBhcmVudCAy"
    "MmRkYzc4MDQ5OThkMzFjOGYxYTI2MTdlZTcyMGUwNjNiMWZhNmNkCnBhcmVudCAzNjQyN2Rk"
    "MjViZjY2NWY0MTA1MjVhMDNhMWYwYTBlYTkyNzUxNTBiCmF1dGhvciBDaHJpcyBUYXlsb3Ig"
    "PG1yY2F0aWRAZ21haWwuY29tPiAxNzExMTkyODE2IC0wNzAwCmNvbW1pdHRlciBHaXRIdWIg"
    "PG5vcmVwbHlAZ2l0aHViLmNvbT4gMTcxMTE5MjgxNiAtMDcwMApncGdzaWcgLS0tLS1CRUdJ"
    "TiBQR1AgU0lHTkFUVVJFLS0tLS0KIAogd3NGY0JBQUJDQUFRQlFKbC9ycndDUkMxYVE3dXU1"
    "VWhsQUFBWUVVUUFKQkh1ZjBtUWRxQW9mclZTR3ZHT2ZRZwogSllQQVNZTmh5azhoUEs0N3E3"
    "REhPWVI0UUN3eHZySC9VaVpYcGMyZndSK3FlNUxSRkhqUjFaNE43K3AvK0J6cgogMnd6cXBZ"
    "NmR1SUt5NzVvdjFNR0RGa0FocmpPdlBXbDFCR2RlcTNCYnRQL0NzQ0JVSVEyOCtNZG5OVGRq"
    "dGprVwogMkJMWUFYN3VnU2pKWGQ1OEE0eTJ2MHJvUDZKeVRWbHR3b0NRcUtlYUlTVkNDa0Y0"
    "amcvcmJTZHRUck9WeE1oZAogenlPdmtTMzlzcGVOUk9UaUlsTTB4N3VESkR6ZFZja2draEV2"
    "N2dDN3Z6VWdzRFRHYUh6WXkwTU5ISDlmaHhBaAogcEdNQXEzV3p4ZGgxdjJLQTF5cExvZ1VV"
    "ekwwSzNFRFpvRjREL1RZMG96SDdESm1TUytzV1BVRWpUSHNXQzg5eAogdFR2bk8vSzhLNjZj"
    "SkZDQ1N1LzBvSm4yVFdhQnZUdldnNXJ5WFN2RlNibmVjNGRzQlAzczBZVS9odGJRbDN4Tgog"
    "TTkvMmw3M21ra2NJNHIxOVRHQU0zV0JzYTVxKzBrUHk5QnE2SHloeTM5ampOV2pkTnVzZFFG"
    "a3p0dXNSVFJsNwogYUx5SlNRTXB4cWVzeHRacWVEMjRmQmZYUmFRMlhhYkxFUytsTThBd0ly"
    "dytMUFFScWcrQ20rb3Rxb21BeG9QaQogemVmcUU3OGxXVkJRenVESkpLZnNHaERGY1BFaWJJ"
    "UXYyN0hOT3BsaTBtWkFJbGFnckJiOEF4cjRkbytVS0d6cAogMkhYYlprTjZFTWVGYVRKNzY4"
    "MUhLZHlQRFdiWUM3STVPRFBJNktkUW1kaWFJNmZldUJ4QTgybjU2R3laQ2w1OQogclFHZnFM"
    "dGU3VFVvRnVMTE53c1YKID00MTBSCiAtLS0tLUVORCBQR1AgU0lHTkFUVVJFLS0tLS0KIAoK"
    "TWVyZ2UgcHVsbCByZXF1ZXN0ICMyMyBmcm9tIGdibGV0cjQyL21hc3RlcgoKSGFuZGxlIHRo"
    "ZSBjYXNlIHdoZW4gbm8gb3JpZ2luYWwgZGF0YSB3YXMgbG9zdA==")
CANDIDATE_TREE = "8" * 40
CANDIDATE_COMMIT_RAW = (
    f"tree {CANDIDATE_TREE}\nauthor Fixture <fixture@example.com> 1 +0000\n"
    "committer Fixture <fixture@example.com> 1 +0000\n\nfixture\n"
).encode("utf-8")
CANDIDATE_COMMIT = hashlib.sha1(
    f"commit {len(CANDIDATE_COMMIT_RAW)}\0".encode("ascii") +
    CANDIDATE_COMMIT_RAW).hexdigest()


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
    "candidate_commit": CANDIDATE_COMMIT,
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
        "cache_hierarchy": [
            {
                "index": 0, "level": "1", "type": "Data", "size": "32K",
                "coherency_line_size": "64", "number_of_sets": "64",
                "ways_of_associativity": "8", "physical_line_partition": "1",
                "shared_cpu_list": "0-1", "shared_cpu_map": "00000003",
                "allocation_policy": None, "write_policy": "WriteBack",
            },
            {
                "index": 1, "level": "3", "type": "Unified", "size": "8M",
                "coherency_line_size": "64", "number_of_sets": "16384",
                "ways_of_associativity": "8", "physical_line_partition": "1",
                "shared_cpu_list": "0-2", "shared_cpu_map": "00000007",
                "allocation_policy": None, "write_policy": "WriteBack",
            },
        ],
        "cache_index_inventory": ["index0", "index1"],
        "cache_directory_inventory_text": runner.exact_text_content(
            "index0\nindex1\n", "fixture cache-directory inventory"),
        "numa_nodes": [0],
        "numa_node_inventory": ["node0"],
        "core_class": {"core_type": None, "cpu_capacity": None},
    }


HOST = {
    "system": {
        "system": "Linux", "node": "fixture", "release": "fixture",
        "version": "fixture", "machine": "x86_64", "python": "3",
        "page_size": 4096,
    },
    "allowed_cpu_set_at_launch": [0, 1, 2],
    "online_cpu_set": [0, 1, 2],
    "online_cpu_list_text": runner.exact_text_content(
        "0-2\n", "fixture online CPU list"),
    "online_node_list_text": runner.exact_text_content(
        "0\n", "fixture online node list"),
    "benchmark_cpu": host_cpu(0),
    "reserved_sibling": host_cpu(1),
    "turbo_and_pstate": {
        "/sys/devices/system/cpu/intel_pstate/no_turbo": "0",
        "/sys/devices/system/cpu/amd_pstate/status": None,
        "/sys/devices/system/cpu/cpufreq/boost": None,
    },
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
    if raw_schema in runner.WORKSPACE_SELECTOR_SCHEMAS:
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


def complete_artifact(path: str, kind: str, character: str,
                      *, content: str | None = None) -> dict:
    encoded = None if content is None else content.encode("utf-8")
    return {
        "path": path,
        "kind": kind,
        "size": 1 if encoded is None else len(encoded),
        "mode": 0o755 if kind in {"executable", "compiler", "archiver",
                                   "ranlib"} else 0o644,
        "mtime_ns": 1,
        "sha256": (character * 64 if encoded is None else
                   runner.sha256_bytes(encoded)),
    }


def complete_build_fixture(role: str) -> dict:
    baseline = role == "baseline"
    build_dir = SPECIFICATION[f"{role}_build_dir"]
    source_root = SPECIFICATION[f"{role}_source_root"]
    if baseline:
        archive_name = "libleopard_main_exact.a"
        executable_name = "leopard_main_benchmark"
        target_directory = "leopard_main_exact.dir"
        benchmark_source = (SPECIFICATION["candidate_source_root"] +
            "/experiments/leopard2/main_compare/legacy_main_benchmark.cpp")
        benchmark_object = (build_dir +
            "/CMakeFiles/leopard_main_benchmark.dir/legacy_main_benchmark.cpp.o")
        cache = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
            "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
            "LEO_MAIN_HAS_MARCH_NATIVE": "1",
        }
        isa_policy = "whole-build -march=native"
        library_names = runner.BASELINE_LIBRARY_SOURCES
        entry_count = 5
    else:
        archive_name = runner.CANONICAL_CMAKE_IDENTITY["archive"]
        executable_name = "bench_leopard2"
        target_directory = runner.CANONICAL_CMAKE_IDENTITY["target_directory"]
        benchmark_source = source_root + "/bench/leopard2/benchmark.cpp"
        benchmark_object = (build_dir +
            "/CMakeFiles/bench_leopard2.dir/benchmark.cpp.o")
        cache = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
            "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
            "ENABLE_OPENMP": "ON",
            "LEO2_BACKEND_VARIANT": "auto",
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_BUILD_TESTS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
        }
        isa_policy = (
            "portable core with ISA flags only on SSSE3, AVX2, and "
            "AVX-512VL translation units")
        library_names = runner.CANDIDATE_LIBRARY_SOURCES
        entry_count = 20
    library_pairs = []
    archive_objects = []
    for index, name in enumerate(library_names):
        if not baseline and name in {
                "Leopard2BackendSSSE3.cpp", "Leopard2BackendAVX2.cpp",
                "Leopard2BackendAVX512.cpp"}:
            backend = name.removeprefix("Leopard2Backend").removesuffix(
                ".cpp").lower()
            object_relative = (
                f"CMakeFiles/leopard2_backend_{backend}.dir/{name}.o")
        else:
            object_relative = f"CMakeFiles/{target_directory}/{name}.o"
        archive_objects.append(object_relative)
        library_pairs.append({
            "source": complete_artifact(
                source_root + "/" + name, "source_file",
                format((index % 9) + 1, "x")),
            "object": complete_artifact(
                build_dir + "/" + object_relative, "object_file",
                format(((index + 3) % 9) + 1, "x")),
        })
    pairs = sorted([
        *library_pairs,
        {
            "source": complete_artifact(
                benchmark_source, "source_file", "5" if baseline else "6"),
            "object": complete_artifact(
                benchmark_object, "object_file", "7" if baseline else "8"),
        },
    ], key=lambda record: record["source"]["path"])
    archive_recipe_text = (
        f"/usr/bin/ar qc {archive_name} "
        + " ".join(archive_objects) + "\n"
        + f"/usr/bin/ranlib {archive_name}\n")
    archive_recipe_content = runner.exact_text_content(
        archive_recipe_text, f"fixture {role} archive recipe")
    executable_recipe_text = (
        f"/usr/bin/compiler -O3 -fopenmp {archive_name} "
        f"{Path(benchmark_object).relative_to(Path(build_dir)).as_posix()} "
        f"-o {executable_name}\n")
    executable_recipe_content = runner.exact_text_content(
        executable_recipe_text, f"fixture {role} executable recipe")
    archive_path = build_dir + "/" + archive_name
    executable_path = build_dir + "/" + executable_name
    executable_recipe_path = (
        build_dir + f"/CMakeFiles/{executable_name}.dir/link.txt")
    archive_recipe_path = (
        build_dir + f"/CMakeFiles/{target_directory}/link.txt")
    compiler_text = "fixture compiler 1.0\n"
    return {
        "build_dir": build_dir,
        "cmake_cache": complete_artifact(
            build_dir + "/CMakeCache.txt", "build_metadata", "9"),
        "compile_commands": complete_artifact(
            build_dir + "/compile_commands.json", "build_metadata", "a"),
        "executable_link_recipe": complete_artifact(
            executable_recipe_path, "build_metadata", "b",
            content=executable_recipe_text),
        "archive_link_recipe": complete_artifact(
            archive_recipe_path, "build_metadata", "c",
            content=archive_recipe_text),
        "compiler": complete_artifact("/usr/bin/compiler", "compiler", "d"),
        "compiler_invocation": {
            "invocation": "/usr/bin/compiler",
            "resolved_path": "/usr/bin/compiler",
        },
        "compiler_version_stdout": {
            "sha256": runner.sha256_bytes(compiler_text.encode("utf-8")),
            "text": compiler_text,
        },
        "archiver": complete_artifact("/usr/bin/ar", "archiver", "e"),
        "ranlib": complete_artifact("/usr/bin/ranlib", "ranlib", "f"),
        "validated_archive_members": [Path(name).name + ".o"
                                       for name in library_names],
        "validated_executable": complete_artifact(
            executable_path, "executable", "0" if baseline else "1"),
        "validated_archive": complete_artifact(
            archive_path, "archive", "2" if baseline else "3"),
        "validated_cache": cache,
        "validated_compile_commands": {
            "entry_count": entry_count,
            "required_sources": sorted(pair["source"]["path"] for pair in pairs),
            "validated_optimization": "-O3",
            "validated_openmp": True,
            "required_source_object_pairs": pairs,
            "isa_policy": isa_policy,
        },
        "archive_link_recipe_content": archive_recipe_content,
        "executable_link_recipe_content": executable_recipe_content,
        "archive_link_tool_invocations": {
            "archiver": {"invocation": "/usr/bin/ar",
                         "resolved_path": "/usr/bin/ar"},
            "ranlib": {"invocation": "/usr/bin/ranlib",
                       "resolved_path": "/usr/bin/ranlib"},
        },
    }


def complete_runtime_fixture(executable: str, character: str) -> dict:
    raw = (
        "linux-vdso.so.1 (0x00000000)\n"
        "libc.so.6 => /lib/libc.so.6 (0x00000000)\n"
        "libm.so.6 => /lib/libm.so.6 (0x00000000)\n"
        "/lib64/ld-linux-x86-64.so.2 (0x00000000)\n")
    return {
        "executable": executable,
        "raw_ldd_output": runner.exact_text_content(
            raw, "fixture raw ldd output"),
        "dependencies": [
            {
                "soname": "ld-linux-x86-64.so.2",
                "loader_path": "/lib64/ld-linux-x86-64.so.2",
                "file": complete_artifact(
                    "/lib64/ld-linux-x86-64.so.2",
                    "dynamic_loader", character),
            },
            {
                "soname": "libc.so.6",
                "loader_path": "/lib/libc.so.6",
                "file": complete_artifact(
                    "/lib/libc.so.6",
                    "shared_library", character),
            },
            {
                "soname": "libm.so.6",
                "loader_path": "/lib/libm.so.6",
                "file": complete_artifact(
                    "/lib/libm.so.6", "shared_library",
                    hex((int(character, 16) + 1) & 15)[2:]),
            },
            {"soname": "linux-vdso.so.1", "virtual": True},
        ],
    }


def complete_source_fixture(role: str) -> dict:
    baseline = role == "baseline"
    raw = (base64.b64decode(BASELINE_COMMIT_BASE64)
           if baseline else CANDIDATE_COMMIT_RAW)
    head = runner.MAIN_COMMIT if baseline else CANDIDATE_COMMIT
    tree = BASELINE_TREE if baseline else CANDIDATE_TREE
    commit_object = runner.git_commit_object_identity(raw, head)
    runner.validate_git_commit_object_identity(
        commit_object, head, tree, f"fixture {role}")
    return {
        "path": SPECIFICATION[f"{role}_source_root"],
        "head": head, "tree": tree, "detached": baseline,
        "tracked_tree_listing_sha256": "7" * 64 if baseline else "9" * 64,
        "tracked_status": "clean", "commit_object": commit_object,
    }


def cmake_fixture_identity(raw_schema: str) -> tuple[dict, dict]:
    if raw_schema == runner.RAW_SCHEMA:
        baseline_build = complete_build_fixture("baseline")
        candidate_build = complete_build_fixture("candidate")
        baseline_executable = baseline_build["validated_executable"]
        candidate_executable = candidate_build["validated_executable"]
        baseline_archive = baseline_build["validated_archive"]
        candidate_archive = candidate_build["validated_archive"]
        compiler = baseline_build["compiler"]
        candidate_build["compiler"] = copy.deepcopy(compiler)
        candidate_build["compiler_version_stdout"] = copy.deepcopy(
            baseline_build["compiler_version_stdout"])
        identity = {
            "runner": complete_artifact(
                SPECIFICATION["runner"], "file", "1"),
            "taskset": complete_artifact(
                SPECIFICATION["taskset"], "executable", "2"),
            "ldd": complete_artifact(
                SPECIFICATION["ldd"], "executable", "3"),
            "baseline_executable": copy.deepcopy(baseline_executable),
            "candidate_executable": copy.deepcopy(candidate_executable),
            "baseline_archive": copy.deepcopy(baseline_archive),
            "candidate_archive": copy.deepcopy(candidate_archive),
            "baseline_build": baseline_build,
            "candidate_build": candidate_build,
            "baseline_runtime_closure": complete_runtime_fixture(
                baseline_executable["path"], "4"),
            "candidate_runtime_closure": complete_runtime_fixture(
                candidate_executable["path"], "5"),
            "baseline_source": complete_source_fixture("baseline"),
            "candidate_source": complete_source_fixture("candidate"),
        }
        specification = copy.deepcopy(SPECIFICATION)
        specification.update({
            "baseline_executable": baseline_executable["path"],
            "candidate_executable": candidate_executable["path"],
            "baseline_archive": baseline_archive["path"],
            "candidate_archive": candidate_archive["path"],
        })
        return identity, specification

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
    if raw_schema in runner.CANDIDATE_MODE_SCHEMAS:
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


def synchronize_identity(value: dict) -> None:
    value["identities_final"] = copy.deepcopy(value["identities_initial"])
    for invocation in value["invocations"]:
        invocation["identity_before"] = copy.deepcopy(value["identities_initial"])
        invocation["identity_after"] = copy.deepcopy(value["identities_initial"])


def write_complete_evidence_bundle(
    root: Path, value: dict, manifest_schema: str = runner.MANIFEST_SCHEMA,
) -> Path:
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
        "schema": manifest_schema,
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
    return manifest_path


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


def replace_current_executable_recipe_text(value: dict, text: str) -> None:
    content = runner.exact_text_content(
        text, "mutated executable link recipe")
    build = value["identities_initial"]["candidate_build"]
    build["executable_link_recipe_content"] = content
    build["executable_link_recipe"]["size"] = content["size"]
    build["executable_link_recipe"]["sha256"] = content["sha256"]
    synchronize_identity(value)


def synthetic_failure(raw_schema: str) -> dict:
    raw = synthetic_raw(raw_schema=raw_schema)
    failure_schema = {
        runner.RAW_SCHEMA_V2: runner.FAILURE_SCHEMA_V2,
        runner.RAW_SCHEMA_V3: runner.FAILURE_SCHEMA_V3,
        runner.RAW_SCHEMA_V4: runner.FAILURE_SCHEMA_V4,
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

    def test_count_and_cpu_fields_require_bounded_non_boolean_integers(self) -> None:
        for key, replacement in (
            ("rounds", True), ("batch", True), ("threads", True),
            ("iterations", True), ("iterations", runner.MAX_CAMPAIGN_COUNT + 1),
            ("reuse", True), ("reuse", runner.MAX_CAMPAIGN_COUNT + 1),
            ("warmup", True), ("warmup", runner.MAX_CAMPAIGN_COUNT + 1),
            ("benchmark_cpu", False),
            ("benchmark_cpu", runner.MAX_CPU_ID + 1),
            ("reserved_sibling", True),
            ("reserved_sibling", runner.MAX_CPU_ID + 1),
        ):
            with self.subTest(kind="raw-campaign", key=key, value=replacement):
                value = synthetic_raw()
                value["campaign"][key] = replacement
                self.assert_rejected(value)

        for key, replacement in (
            ("rounds", True), ("batch", True), ("threads", True),
            ("iterations", True), ("iterations", runner.MAX_CAMPAIGN_COUNT + 1),
            ("reuse", True), ("reuse", runner.MAX_CAMPAIGN_COUNT + 1),
            ("warmup", True), ("warmup", runner.MAX_CAMPAIGN_COUNT + 1),
            ("benchmark_cpu", False),
            ("benchmark_cpu", runner.MAX_CPU_ID + 1),
            ("reserved_sibling", True),
            ("reserved_sibling", runner.MAX_CPU_ID + 1),
        ):
            with self.subTest(kind="failed-campaign", key=key, value=replacement):
                failure = synthetic_failure(runner.RAW_SCHEMA)
                failure["campaign"][key] = replacement
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_failure(
                        resign(failure), Path("/unused"), check_files=False)

        for invalid in (False, runner.MAX_CPU_ID + 1):
            with self.subTest(kind="topology", value=invalid), \
                 self.assertRaises(runner.EvidenceError):
                runner.validate_topology(invalid, 1)
            with self.subTest(kind="reservation-request", value=invalid), \
                 self.assertRaises(runner.EvidenceError):
                runner.parse_reservation(
                    runner.canonical_bytes(RESERVATION_PAYLOAD), invalid, 1)

        for key in ("benchmark_cpu", "reserved_sibling"):
            for invalid in (False, runner.MAX_CPU_ID + 1):
                with self.subTest(kind="reservation-payload", key=key,
                                  value=invalid):
                    payload = copy.deepcopy(RESERVATION_PAYLOAD)
                    payload[key] = invalid
                    with self.assertRaises(runner.EvidenceError):
                        runner.parse_reservation(
                            runner.canonical_bytes(payload), 0, 1)

        for field in ("allowed_cpu_set_at_launch", "online_cpu_set"):
            for invalid in (False, runner.MAX_CPU_ID + 1):
                with self.subTest(kind="host", field=field, value=invalid):
                    host = copy.deepcopy(HOST)
                    host[field][0] = invalid
                    with self.assertRaises(runner.EvidenceError):
                        runner.validate_host_record(
                            host, 0, 1, CAMPAIGN["allowed_cpu_set_at_launch"],
                            runner.RAW_SCHEMA)

        for invalid in (False, runner.MAX_CPU_ID + 1):
            with self.subTest(kind="campaign-allowed", value=invalid):
                value = synthetic_raw()
                value["campaign"]["allowed_cpu_set_at_launch"][0] = invalid
                self.assert_rejected(value)
            with self.subTest(kind="failed-campaign-allowed", value=invalid):
                failure = synthetic_failure(runner.RAW_SCHEMA)
                failure["campaign"]["allowed_cpu_set_at_launch"][0] = invalid
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_failure(
                        resign(failure), Path("/unused"), check_files=False)
            with self.subTest(kind="failed-host-online", value=invalid):
                failure = synthetic_failure(runner.RAW_SCHEMA)
                failure["host_initial"]["online_cpu_set"][0] = invalid
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_failure(
                        resign(failure), Path("/unused"), check_files=False)

        for field, replacement in (
            ("k", True), ("r", True), ("shard_bytes", True),
            ("shard_bytes", runner.MAX_SHARD_BYTES + 64),
            ("losses", True), ("seed", True), ("seed", runner.MASK64 + 1),
        ):
            with self.subTest(kind="cell", field=field, value=replacement):
                values = asdict(CELL)
                values[field] = replacement
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_cell(runner.Cell(**values))

    def test_uniformly_incomplete_v5_identities_fail_offline_verify_and_promotion(
        self,
    ) -> None:
        plan = load_plan_runner()

        with tempfile.TemporaryDirectory() as directory:
            valid_manifest = write_complete_evidence_bundle(
                Path(directory), synthetic_raw())
            runner.verify_campaign(argparse.Namespace(
                manifest=valid_manifest, no_current_input_check=True))
            document, scope = plan.verify_exact_manifest(valid_manifest)
            self.assertEqual(document["schema"], runner.MANIFEST_SCHEMA)
            self.assertEqual(
                plan.validate_evidence_scope(scope)["schema"],
                "leopard2-balanced-evidence-scope/v3")

        def missing_sources(value: dict) -> None:
            for role in ("baseline", "candidate"):
                value["identities_initial"][f"{role}_source"].pop("tree")
            synchronize_identity(value)

        def incomplete_runtime_files(value: dict) -> None:
            for role in ("baseline", "candidate"):
                dependency = value["identities_initial"][
                    f"{role}_runtime_closure"]["dependencies"][0]
                dependency["file"] = {"path": dependency["file"]["path"]}
            synchronize_identity(value)

        def path_only_compile_pairs(value: dict) -> None:
            for role in ("baseline", "candidate"):
                pairs = value["identities_initial"][f"{role}_build"][
                    "validated_compile_commands"]["required_source_object_pairs"]
                for pair in pairs:
                    pair["source"] = {"path": pair["source"]["path"]}
                    pair["object"] = {"path": pair["object"]["path"]}
            synchronize_identity(value)

        def empty_toolchain_records(value: dict) -> None:
            for role in ("baseline", "candidate"):
                build = value["identities_initial"][f"{role}_build"]
                for name in ("compiler", "cmake_cache", "compile_commands",
                             "executable_link_recipe", "archive_link_recipe"):
                    build[name] = {}
            synchronize_identity(value)

        def reduced_outputs(value: dict) -> None:
            identities = value["identities_initial"]
            for role in ("baseline", "candidate"):
                for output in ("archive", "executable"):
                    key = f"{role}_{output}"
                    old = identities[key]
                    reduced = {
                        "path": old["path"], "kind": old["kind"],
                        "sha256": old["sha256"],
                    }
                    identities[key] = copy.deepcopy(reduced)
                    identities[f"{role}_build"][f"validated_{output}"] = \
                        copy.deepcopy(reduced)
            synchronize_identity(value)

        def topology_only_host(value: dict) -> None:
            for name in ("benchmark_cpu", "reserved_sibling"):
                record = value["host_initial"][name]
                value["host_initial"][name] = {
                    "cpu": record["cpu"], "topology": record["topology"]}
            value["host_final"] = copy.deepcopy(value["host_initial"])

        def truncated_tu_closure(value: dict) -> None:
            for role, names in (("baseline", runner.BASELINE_LIBRARY_SOURCES),
                                ("candidate", runner.CANDIDATE_LIBRARY_SOURCES)):
                build = value["identities_initial"][f"{role}_build"]
                semantics = build["validated_compile_commands"]
                suffix = "/" + names[-1]
                pair = next(item for item in semantics[
                    "required_source_object_pairs"]
                    if item["source"]["path"].endswith(suffix))
                semantics["required_source_object_pairs"].remove(pair)
                semantics["required_sources"].remove(pair["source"]["path"])
                semantics["entry_count"] -= 1
                member = Path(pair["object"]["path"]).name
                build["validated_archive_members"].remove(member)
                relative = Path(pair["object"]["path"]).relative_to(
                    Path(build["build_dir"])).as_posix()
                text = build["archive_link_recipe_content"]["text"].replace(
                    " " + relative, "", 1)
                content = runner.exact_text_content(
                    text, f"truncated {role} archive recipe")
                build["archive_link_recipe_content"] = content
                build["archive_link_recipe"]["size"] = content["size"]
                build["archive_link_recipe"]["sha256"] = content["sha256"]
            synchronize_identity(value)

        def missing_dynamic_loader(value: dict) -> None:
            for role in ("baseline", "candidate"):
                closure = value["identities_initial"][f"{role}_runtime_closure"]
                closure["dependencies"] = [item for item in closure["dependencies"]
                    if item["soname"] != "ld-linux-x86-64.so.2"]
                text = "\n".join(line for line in
                    closure["raw_ldd_output"]["text"].splitlines()
                    if not line.startswith("/lib64/ld-linux")) + "\n"
                closure["raw_ldd_output"] = runner.exact_text_content(
                    text, f"truncated {role} ldd output")
            synchronize_identity(value)

        def swapped_runtime_file_records(value: dict) -> None:
            for role in ("baseline", "candidate"):
                dependencies = value["identities_initial"][
                    f"{role}_runtime_closure"]["dependencies"]
                libc = next(item for item in dependencies
                            if item["soname"] == "libc.so.6")
                libm = next(item for item in dependencies
                            if item["soname"] == "libm.so.6")
                libc["file"], libm["file"] = libm["file"], libc["file"]
            synchronize_identity(value)

        def noncanonical_runtime_loader_paths(value: dict) -> None:
            for role in ("baseline", "candidate"):
                closure = value["identities_initial"][f"{role}_runtime_closure"]
                libc = next(item for item in closure["dependencies"]
                            if item["soname"] == "libc.so.6")
                libc["loader_path"] = "/lib/./libc.so.6"
                libc["file"]["path"] = libc["loader_path"]
                text = closure["raw_ldd_output"]["text"].replace(
                    "/lib/libc.so.6", libc["loader_path"])
                closure["raw_ldd_output"] = runner.exact_text_content(
                    text, f"noncanonical {role} ldd output")
            synchronize_identity(value)

        def truncated_cache_inventory(value: dict) -> None:
            for name in ("benchmark_cpu", "reserved_sibling"):
                value["host_initial"][name]["cache_hierarchy"].pop()
                value["host_initial"][name]["cache_index_inventory"].pop()
            value["host_final"] = copy.deepcopy(value["host_initial"])

        def empty_numa_summary(value: dict) -> None:
            for name in ("benchmark_cpu", "reserved_sibling"):
                value["host_initial"][name]["numa_nodes"] = []
                value["host_initial"][name]["numa_node_inventory"] = []
            value["host_final"] = copy.deepcopy(value["host_initial"])

        def source_tree_commit_mismatch(value: dict) -> None:
            for role in ("baseline", "candidate"):
                value["identities_initial"][f"{role}_source"]["tree"] = "f" * 40
            synchronize_identity(value)

        def boolean_campaign_counts(value: dict) -> None:
            value["campaign"]["iterations"] = True

        mutations = {
            "sources": missing_sources,
            "runtime-files": incomplete_runtime_files,
            "path-only-compile-pairs": path_only_compile_pairs,
            "empty-toolchain-records": empty_toolchain_records,
            "reduced-outputs": reduced_outputs,
            "topology-only-host": topology_only_host,
            "truncated-tu-closure": truncated_tu_closure,
            "missing-dynamic-loader": missing_dynamic_loader,
            "swapped-runtime-file-records": swapped_runtime_file_records,
            "noncanonical-runtime-loader-paths":
                noncanonical_runtime_loader_paths,
            "truncated-cache-inventory": truncated_cache_inventory,
            "empty-numa-summary": empty_numa_summary,
            "source-tree-commit-mismatch": source_tree_commit_mismatch,
            "boolean-campaign-counts": boolean_campaign_counts,
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                value = synthetic_raw()
                mutate(value)
                manifest_path = write_complete_evidence_bundle(root, value)
                with self.assertRaises(runner.EvidenceError):
                    runner.verify_campaign(argparse.Namespace(
                        manifest=manifest_path, no_current_input_check=True))
                with self.assertRaises(plan.PlanError):
                    plan.verify_exact_manifest(manifest_path)

    def test_coherent_ordinary_runtime_rewrite_is_internal_consistency_only(
        self,
    ) -> None:
        plan = load_plan_runner()
        value = synthetic_raw()
        for role in ("baseline", "candidate"):
            closure = value["identities_initial"][f"{role}_runtime_closure"]
            closure["dependencies"] = [
                dependency for dependency in closure["dependencies"]
                if dependency["soname"] != "libm.so.6"]
            text = "\n".join(
                line for line in closure["raw_ldd_output"]["text"].splitlines()
                if not line.startswith("libm.so.6 =>")) + "\n"
            closure["raw_ldd_output"] = runner.exact_text_content(
                text, f"coherently rewritten {role} ldd output")
        synchronize_identity(value)
        with tempfile.TemporaryDirectory() as directory:
            manifest_path = write_complete_evidence_bundle(
                Path(directory), value)
            self.assertEqual(runner.verify_campaign(argparse.Namespace(
                manifest=manifest_path, no_current_input_check=True)), 0)
            document, scope = plan.verify_exact_manifest(manifest_path)
            self.assertEqual(document["schema"], runner.MANIFEST_SCHEMA)
            plan.validate_evidence_scope(scope)

    def test_coherent_cache_inventory_rewrite_is_internal_consistency_only(
        self,
    ) -> None:
        plan = load_plan_runner()
        value = synthetic_raw()
        for name in ("benchmark_cpu", "reserved_sibling"):
            record = value["host_initial"][name]
            record["cache_hierarchy"].pop()
            record["cache_index_inventory"].pop()
            text = "".join(
                f"{entry}\n" for entry in record["cache_index_inventory"])
            record["cache_directory_inventory_text"] = \
                runner.exact_text_content(
                    text, f"coherently rewritten {name} cache inventory")
        value["host_final"] = copy.deepcopy(value["host_initial"])
        with tempfile.TemporaryDirectory() as directory:
            manifest_path = write_complete_evidence_bundle(
                Path(directory), value)
            self.assertEqual(runner.verify_campaign(argparse.Namespace(
                manifest=manifest_path, no_current_input_check=True)), 0)
            document, scope = plan.verify_exact_manifest(manifest_path)
            self.assertEqual(document["schema"], runner.MANIFEST_SCHEMA)
            plan.validate_evidence_scope(scope)

    def test_promotion_consumes_the_exact_verified_manifest_and_raw_snapshots(
        self,
    ) -> None:
        plan = load_plan_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = write_complete_evidence_bundle(
                root, synthetic_raw())
            raw_path = root / "raw.json"
            exact_runner = plan.load_exact_main_runner()

            class InterleavingRunner:
                accepted = None

                def verified_campaign_bundle(
                    self, path: Path, no_current_input_check: bool = False,
                ):
                    self.accepted = exact_runner.verified_campaign_bundle(
                        path, no_current_input_check)
                    path.write_text("{}\n", encoding="utf-8")
                    raw_path.write_text("{}\n", encoding="utf-8")
                    return self.accepted

            interleaving = InterleavingRunner()
            with mock.patch.object(
                plan, "load_exact_main_runner", return_value=interleaving,
            ):
                document, scope = plan.verify_exact_manifest(manifest_path)
            self.assertIsNotNone(interleaving.accepted)
            self.assertEqual(document, interleaving.accepted[0])
            plan.validate_evidence_scope(scope)
            self.assertEqual(manifest_path.read_text(encoding="utf-8"), "{}\n")
            self.assertEqual(raw_path.read_text(encoding="utf-8"), "{}\n")

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

    def test_legacy_v4_raw_fixture_remains_replayable(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V4)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_legacy_v4_manifest_fixture_remains_replayable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manifest = write_complete_evidence_bundle(
                Path(directory), synthetic_raw(raw_schema=runner.RAW_SCHEMA_V4),
                runner.MANIFEST_SCHEMA_V4)
            self.assertEqual(runner.verify_campaign(argparse.Namespace(
                manifest=manifest, no_current_input_check=True)), 0)

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

    def test_executable_recipe_requires_the_declared_archive_operand(self) -> None:
        canonical = complete_build_fixture("candidate")[
            "executable_link_recipe_content"]["text"]
        mutations = {
            "absolute same-basename archive": canonical.replace(
                " libleopard.a ", " /tmp/libleopard.a ", 1),
            "sibling same-basename archive": canonical.replace(
                " libleopard.a ", " ../foreign/libleopard.a ", 1),
            "duplicate same-basename archive": canonical.replace(
                " libleopard.a ",
                " libleopard.a /tmp/libleopard.a ", 1),
        }
        for label, recipe in mutations.items():
            with self.subTest(label=label):
                value = synthetic_raw()
                replace_current_executable_recipe_text(value, recipe)
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

    def test_v5_manifest_binds_complete_identity_and_replays_portably(self) -> None:
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
                "campaign": copy.deepcopy(value["campaign"]),
                "host_initial": copy.deepcopy(HOST),
                "reservation": copy.deepcopy(RESERVATION),
                "pair_lease": copy.deepcopy(PAIR_LEASE),
                "isolation": copy.deepcopy(ISOLATION),
                "input_specification": copy.deepcopy(value["input_specification"]),
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
            with self.assertRaises(runner.EvidenceError):
                runner.bounded_file_contents_snapshot(fifo)

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
