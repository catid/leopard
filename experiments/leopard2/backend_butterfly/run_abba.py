#!/usr/bin/env python3
"""Run and replay portable, fail-closed Leopard2 butterfly ABBA evidence."""

from __future__ import print_function

import argparse
import base64
import binascii
import copy
import ctypes
import errno
import hashlib
import importlib.util
import json
import math
import os
import platform
import re
import selectors
import shlex
import shutil
import signal
import socket
import stat
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

try:
    import fcntl
except ImportError:  # Verification remains usable on non-POSIX hosts.
    fcntl = None


LEGACY_SCHEMA = "leopard2-backend-butterfly-abba/v6"
PRE_XOR_SCHEMA = "leopard2-backend-butterfly-abba/v7"
SCHEMA = "leopard2-backend-butterfly-abba/v8"
RAW_SCHEMA = "leopard2-backend-butterfly-raw/v1"
RESERVATION_SCHEMA = "leopard2-cpu-reservation/v1"
PAIR_LEASE_SCHEMA = "leopard2-cpu-pair-lease/v1"
MATRIX_SCHEMA_V1 = "leopard2-backend-matrix/v1"
MATRIX_SCHEMA = "leopard2-backend-matrix/v2"
SUPPORTED_BACKENDS = ("ssse3", "avx2")
SEQUENCES = (("A1", "baseline"), ("B1", "candidate"),
             ("B2", "candidate"), ("A2", "baseline"))
ROUNDS = tuple(range(1, 17))
# This campaign qualifies a correctness-driven backend-routing refactor, not a
# speculative speedup.  "target" identifies the primary rate/size regimes;
# both those cells and their neighbors use the bead's one-sided non-regression
# floor.  Positive gains remain visible in the replayed summary but are not a
# promotion claim.
TARGET_THRESHOLD = -2.0
NEIGHBOR_FLOOR = -2.0
INVOCATION_SAMPLES = 7
# One-sided 95% Student-t critical value with sixteen ABBA rounds (df=15).
ONE_SIDED_T95 = 1.7530503556925547
FRESH_BUILD_ENVIRONMENT = {
    "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin"}
MAX_BENCHMARK_STDOUT_BYTES = 64 * 1024 * 1024
MAX_BENCHMARK_STDERR_BYTES = 8 * 1024 * 1024
MAX_LINK_RECIPE_BYTES = 1024 * 1024
CHILD_REAP_TIMEOUT_SECONDS = 5.0
PR_SET_CHILD_SUBREAPER = 36
PR_GET_CHILD_SUBREAPER = 37
# Keep authenticated policy rejection distinct from both evidence failure (1)
# and argparse's command-line usage error (2).
EXPECTED_POLICY_FAILURE_EXIT = 3


def load_current_matrix_contract():
    """Load the exact v2 producer contract used by current evidence.

    Historical v1 constants remain frozen below.  Current evidence refuses a
    producer schema change instead of silently accepting a drifted matrix.
    """
    path = Path(__file__).resolve().parents[3] / \
        "tools/leopard2_backend_matrix.py"
    specification = importlib.util.spec_from_file_location(
        "leopard2_backend_matrix_contract_v2", str(path))
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load backend-matrix v2 contract")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    contract = module.evidence_contract()
    expected_keys = {
        "schema", "source_files", "expected_compile_source_counts",
        "compare_tests", "run_only_tests", "run_tests", "test_specs",
        "build_targets", "build_cache_keys", "base_backend_failure_tests",
        "avx512_backend_failure_tests", "backend_failure_ctest_regex",
        "portable_ctest_regex", "cuda_ctest_regex",
    }
    if not isinstance(contract, dict) or set(contract) != expected_keys or \
            contract.get("schema") != MATRIX_SCHEMA:
        raise RuntimeError("backend-matrix producer contract/schema changed")
    return contract


CURRENT_MATRIX_CONTRACT = load_current_matrix_contract()


def cell(name, k, r, profile, field, byte_count, loss, kind):
    resolved_profile = "legacy_high_v1" if profile == "high" else "low_v1"
    active_side = r if profile == "high" else k
    padded_side = 1 << max(0, (active_side - 1).bit_length())
    parent_input = k + padded_side if profile == "high" else padded_side + r
    parent = 1 << max(0, (parent_input - 1).bit_length())
    rounded = (byte_count + 63) & ~63
    return {
        "name": name,
        "K": k,
        "R": r,
        "profile": profile,
        "resolved_profile": resolved_profile,
        "field": field,
        "bytes": byte_count,
        "loss": loss,
        "kind": kind,
        "minimum_speedup_percent": (
            TARGET_THRESHOLD if kind == "target" else NEIGHBOR_FLOOR),
        "parent_count": parent,
        "padded_side": padded_side,
        "rounded_transform_bytes": rounded,
        "tail_mapping": (
            "direct complete 64-byte transform tiles" if rounded == byte_count
            else "public compact/partial shard staged into complete 64-byte transform tiles"),
    }


CELLS = (
    cell("high-gf8-64b", 240, 16, "high", "gf8", 64, 4, "neighbor"),
    cell("high-gf8-tail", 240, 16, "high", "gf8", 65, 4, "neighbor"),
    cell("high-gf8-1k", 240, 16, "high", "gf8", 1024, 4, "neighbor"),
    cell("high-gf8-1025", 240, 16, "high", "gf8", 1025, 4, "neighbor"),
    cell("high-gf8", 240, 16, "high", "gf8", 65536, 4, "target"),
    cell("low-gf8-64b", 32, 224, "low", "gf8", 64, 16, "neighbor"),
    cell("low-gf8-tail", 32, 224, "low", "gf8", 65, 16, "neighbor"),
    cell("low-gf8-1k", 32, 224, "low", "gf8", 1024, 16, "neighbor"),
    cell("low-gf8-1025", 32, 224, "low", "gf8", 1025, 16, "neighbor"),
    cell("low-gf8", 32, 224, "low", "gf8", 65536, 16, "target"),
    cell("high-gf16-64b", 1000, 200, "high", "gf16", 64, 8, "neighbor"),
    cell("high-gf16-tail", 1000, 200, "high", "gf16", 66, 8, "neighbor"),
    cell("high-gf16-130", 1000, 200, "high", "gf16", 130, 8, "neighbor"),
    cell("high-gf16-1k", 1000, 200, "high", "gf16", 1024, 8, "neighbor"),
    cell("high-gf16", 1000, 200, "high", "gf16", 16384, 8, "target"),
    cell("low-gf16-64b", 128, 1024, "low", "gf16", 64, 16, "neighbor"),
    cell("low-gf16-tail", 128, 1024, "low", "gf16", 66, 16, "neighbor"),
    cell("low-gf16-130", 128, 1024, "low", "gf16", 130, 16, "neighbor"),
    cell("low-gf16-1k", 128, 1024, "low", "gf16", 1024, 16, "neighbor"),
    cell("low-gf16", 128, 1024, "low", "gf16", 16384, 16, "target"),
    # Balanced-rate neighbors exercise the same kernels without extrapolating
    # the high/low-rate target result across a qualitatively different shape.
    cell("balanced-gf8", 128, 128, "high", "gf8", 65536, 8, "neighbor"),
    cell("balanced-gf16", 512, 512, "high", "gf16", 16384, 8, "neighbor"),
    # A second decode-loss geometry prevents a single erasure count from
    # becoming the accidental promotion oracle for each production regime.
    cell("high-gf8-loss1", 240, 16, "high", "gf8", 65536, 1, "neighbor"),
    cell("low-gf8-loss1", 32, 224, "low", "gf8", 65536, 1, "neighbor"),
    cell("high-gf16-loss1", 1000, 200, "high", "gf16", 16384, 1, "neighbor"),
    cell("low-gf16-loss1", 128, 1024, "low", "gf16", 16384, 1, "neighbor"),
)

PRE_XOR_BUILD_TRANSLATION_UNITS = (
    "leopard.cpp",
    "leopard2.cpp",
    "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
    "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp",
    "LeopardCommon.cpp",
    "LeopardFF8.cpp",
    "LeopardFF16.cpp",
    "bench/leopard2/benchmark.cpp",
)

BUILD_TRANSLATION_UNITS = (
    *PRE_XOR_BUILD_TRANSLATION_UNITS[:6],
    "Leopard2BackendAVX2Xor.cpp",
    "Leopard2BackendAVX512.cpp",
    *PRE_XOR_BUILD_TRANSLATION_UNITS[6:],
)

# CMake appends GF8/GF16 and then the optional object-library sources to the
# ordinary LIB_SOURCE_FILES list.  This is a semantic part of pre-XOR v7's
# retained archive recipe, not a sorted-set convention.
PRE_XOR_ARCHIVE_TRANSLATION_UNITS = (
    "leopard.cpp",
    "leopard2.cpp",
    "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp",
    "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp",
    "LeopardCommon.cpp",
    "LeopardFF8.cpp",
    "LeopardFF16.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
)

CURRENT_ARCHIVE_TRANSLATION_UNITS = (
    *PRE_XOR_ARCHIVE_TRANSLATION_UNITS,
    "Leopard2BackendAVX2Xor.cpp",
    "Leopard2BackendAVX512.cpp",
)

HISTORICAL_ARCHIVE_TRANSLATION_UNITS = (
    "leopard.cpp",
    "leopard2.cpp",
    "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp",
    "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp",
    "LeopardCommon.cpp",
    "LeopardFF16.cpp",
    "LeopardFF8.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
)

PRE_XOR_MATRIX_SOURCE_FILES = (
    "CMakeLists.txt",
    "LeopardCommon.cpp",
    "LeopardCommon.h",
    "Leopard2Backend.cpp",
    "Leopard2Backend.h",
    "Leopard2BackendScalar.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
    "Leopard2CpuFeatures.cpp",
    "LeopardFF8.cpp",
    "LeopardFF8.h",
    "LeopardFF16.cpp",
    "LeopardFF16.h",
    "Leopard2Plan.cpp",
    "Leopard2Plan.h",
    "Leopard2Direct.h",
    "Leopard2Dispatch.h",
    "leopard.cpp",
    "leopard.h",
    "leopard2.cpp",
    "leopard2.h",
    "tests/leopard2/test_legacy_golden.cpp",
    "tests/leopard2/test_backend_failures.cpp",
    "tests/leopard2/test_backend_ops.cpp",
    "tests/leopard2/test_context_backends.cpp",
    "tests/leopard2/test_r1_xor.cpp",
    "tests/leopard2/test_field_options.cpp",
    "tests/leopard2/test_concurrent_initialization.cpp",
    "tests/leopard2/test_legacy_simd_init_failure.cpp",
    "tests/leopard2/test_initialization_threads.cpp",
    "tests/leopard2/legacy_golden_vectors.h",
    "tests/leopard2/test_api.cpp",
    "tests/leopard2/test_public_api_contract.cpp",
    "tests/leopard2/test_batch_aliasing.cpp",
    "tests/leopard2/test_random.cpp",
    "tests/leopard2/test_locator.cpp",
    "tests/leopard2/test_boundaries.cpp",
    "tests/leopard2/test_active_lch.cpp",
    "tests/leopard2/test_gf16_tails.cpp",
    "tests/leopard2/test_gf16_padded_odd.cpp",
    "tests/leopard2/test_encoder_gf16_legacy_matrix.cpp",
    "tests/leopard2/test_low_gf16_direct_rows.cpp",
    "tests/leopard2/test_decode_high_acceptance.cpp",
    "tests/leopard2/test_high_pruned_legacy.cpp",
    "tests/leopard2/test_decode_low_acceptance.cpp",
    "tests/leopard2/test_decode_plan_schedule.cpp",
    "tests/leopard2/test_direct_encode.cpp",
    "tests/leopard2/test_arbitrary_counts_acceptance.cpp",
    "tests/leopard2/test_max_counts.cpp",
    "tests/leopard2/test_encode_concurrency.cpp",
    "tests/leopard2/test_codec_options_abi.c",
    "tests/leopard2/test_transform_differential.cpp",
    "tests/leopard2/test_pruned_transform.cpp",
    "tests/leopard2/direct_oracle.cpp",
    "tests/leopard2/direct_oracle.h",
    "tests/leopard2/test_direct_oracle.cpp",
    "tests/leopard2/direct_repair.cpp",
    "tests/leopard2/direct_repair.h",
    "tests/leopard2/test_direct_repair.cpp",
    "tests/leopard2/fuzz_api.cpp",
    "tests/leopard2/fuzz_pruned_transform.cpp",
    "tests/leopard2/fuzz_replay.cpp",
    "tests/cmake/test_cuda_optional.cmake",
    "cmake/leopardConfig.cmake.in",
    "tools/check_leopard2_portable_isa.sh",
    "tools/leopard2_backend_matrix.py",
)

PRE_XOR_MATRIX_COMPARE_TESTS = (
    "field_options",
    "direct_oracle",
    "backend_ops",
    "context_backends",
    "r1_xor",
    "legacy_golden",
    "api",
    "public_api_contract",
    "initialization_threads_legacy",
    "initialization_threads_explicit",
    "initialization_threads_default",
    "initialization_threads_gf16_high_codec",
    "initialization_threads_gf16_high_plan",
    "initialization_threads_gf16_low_plan",
    "random",
    "locator",
    "active_lch",
    "gf16_tails",
    "gf16_padded_odd",
    "gf16_legacy_encoder_matrix",
    "low_gf16_direct_rows",
    "decode_high_acceptance",
    "high_pruned_legacy",
    "decode_low_acceptance",
    "decode_plan_schedule",
    "direct_encode",
    "arbitrary_counts_acceptance",
    "max_counts",
    "encode_concurrency",
    "codec_options_abi",
    "direct_repair",
    "boundaries",
    "transform_differential",
    "fuzz_smoke",
    "pruned_fuzz_smoke",
)

# This test intentionally reports which qualified backends were enumerated, so
# its stdout differs between forced variants even though every byte check must
# pass.  Preserve execution and provenance without requiring identical output.
PRE_XOR_MATRIX_RUN_ONLY_TESTS = ("pruned_transform",)
PRE_XOR_MATRIX_RUN_TESTS = \
    PRE_XOR_MATRIX_COMPARE_TESTS + PRE_XOR_MATRIX_RUN_ONLY_TESTS

PRE_XOR_MATRIX_BACKEND_FAILURE_TESTS = (
    "leopard2_backend_failure_scalar_ff8_allocation",
    "leopard2_backend_failure_scalar_ff16_allocation",
    "leopard2_backend_failure_scalar_kat",
    "leopard2_backend_failure_ssse3_ff8_allocation",
    "leopard2_backend_failure_ssse3_ff16_allocation",
    "leopard2_backend_failure_ssse3_kat",
    "leopard2_backend_failure_avx2_ff8_allocation",
    "leopard2_backend_failure_avx2_ff16_allocation",
    "leopard2_backend_failure_avx2_kat",
)

PRE_XOR_MATRIX_TEST_SPECS = {
    "field_options": ("leopard2_field_options_test", []),
    "direct_oracle": ("leopard2_direct_oracle_test", []),
    "backend_ops": ("leopard2_backend_ops_test", []),
    "context_backends": ("leopard2_context_backends_test", []),
    "r1_xor": ("leopard2_r1_xor_test", []),
    "legacy_golden": ("leopard2_legacy_golden_test", []),
    "api": ("leopard2_api_test", []),
    "public_api_contract": ("leopard2_public_api_contract_test", []),
    "initialization_threads_legacy": (
        "leopard2_initialization_threads_test", ["legacy"]),
    "initialization_threads_explicit": (
        "leopard2_initialization_threads_test", ["explicit"]),
    "initialization_threads_default": (
        "leopard2_initialization_threads_test", ["default"]),
    "initialization_threads_gf16_high_codec": (
        "leopard2_initialization_threads_test", ["gf16-high-codec"]),
    "initialization_threads_gf16_high_plan": (
        "leopard2_initialization_threads_test", ["gf16-high-plan"]),
    "initialization_threads_gf16_low_plan": (
        "leopard2_initialization_threads_test", ["gf16-low-plan"]),
    "random": ("leopard2_random_test", [
        "--seed", "0x4c656f7061726432", "--cases", "64", "--threads", "1"]),
    "locator": ("leopard2_locator_test", []),
    "active_lch": ("leopard2_active_lch_test", []),
    "gf16_tails": ("leopard2_gf16_tails_test", []),
    "gf16_padded_odd": ("leopard2_gf16_padded_odd_test", []),
    "gf16_legacy_encoder_matrix": (
        "leopard2_gf16_legacy_encoder_matrix_test", []),
    "low_gf16_direct_rows": ("leopard2_low_gf16_direct_rows_test", []),
    "decode_high_acceptance": ("leopard2_decode_high_acceptance_test", []),
    "high_pruned_legacy": ("leopard2_high_pruned_legacy_test", []),
    "decode_low_acceptance": ("leopard2_decode_low_acceptance_test", []),
    "decode_plan_schedule": ("leopard2_decode_plan_schedule_test", []),
    "direct_encode": ("leopard2_direct_encode_test", []),
    "arbitrary_counts_acceptance": (
        "leopard2_arbitrary_counts_acceptance_test", []),
    "max_counts": ("leopard2_max_counts_test", []),
    "encode_concurrency": ("leopard2_encode_concurrency_test", []),
    "codec_options_abi": ("leopard2_codec_options_abi_test", []),
    "direct_repair": ("leopard2_direct_repair_test", []),
    "boundaries": ("leopard2_boundaries_test", []),
    "transform_differential": ("leopard2_transform_differential_test", []),
    "pruned_transform": ("leopard2_pruned_transform_test", []),
    "fuzz_smoke": ("leopard2_fuzz_smoke", []),
    "pruned_fuzz_smoke": ("leopard2_pruned_fuzz_smoke", []),
}

PRE_XOR_MATRIX_BUILD_TARGETS = (
    "leopard2_field_options_test",
    "leopard2_direct_oracle_test",
    "leopard2_backend_ops_test",
    "leopard2_backend_failures_test",
    "leopard2_context_backends_test",
    "leopard2_r1_xor_test",
    "leopard2_legacy_golden_test",
    "leopard2_api_test",
    "leopard2_public_api_contract_test",
    "leopard2_initialization_threads_test",
    "leopard2_random_test",
    "leopard2_locator_test",
    "leopard2_active_lch_test",
    "leopard2_gf16_tails_test",
    "leopard2_gf16_padded_odd_test",
    "leopard2_gf16_legacy_encoder_matrix_test",
    "leopard2_low_gf16_direct_rows_test",
    "leopard2_decode_high_acceptance_test",
    "leopard2_high_pruned_legacy_test",
    "leopard2_decode_low_acceptance_test",
    "leopard2_decode_plan_schedule_test",
    "leopard2_direct_encode_test",
    "leopard2_arbitrary_counts_acceptance_test",
    "leopard2_max_counts_test",
    "leopard2_encode_concurrency_test",
    "leopard2_codec_options_abi_test",
    "leopard2_direct_repair_test",
    "leopard2_boundaries_test",
    "leopard2_transform_differential_test",
    "leopard2_pruned_transform_test",
    "leopard2_fuzz_smoke",
    "leopard2_pruned_fuzz_smoke",
)

PRE_XOR_MATRIX_BUILD_CACHE_KEYS = (
    "CMAKE_BUILD_TYPE", "CMAKE_GENERATOR", "CMAKE_C_FLAGS",
    "CMAKE_C_FLAGS_RELEASE", "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE",
    "CMAKE_EXE_LINKER_FLAGS", "CMAKE_EXE_LINKER_FLAGS_RELEASE",
    "CMAKE_STATIC_LINKER_FLAGS", "CMAKE_STATIC_LINKER_FLAGS_RELEASE",
    "ENABLE_OPENMP", "LEO2_BACKEND_VARIANT", "LEO2_BUILD_TESTS",
    "LEO2_BUILD_BENCHMARKS", "LEO2_BUILD_FUZZERS", "LEO2_ENABLE_CUDA",
    "CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER",
)

MATRIX_COMPARE_TESTS = tuple(CURRENT_MATRIX_CONTRACT["compare_tests"])
MATRIX_RUN_ONLY_TESTS = tuple(CURRENT_MATRIX_CONTRACT["run_only_tests"])
MATRIX_RUN_TESTS = tuple(CURRENT_MATRIX_CONTRACT["run_tests"])
MATRIX_TEST_SPECS = {
    name: (specification[0], list(specification[1]))
    for name, specification in CURRENT_MATRIX_CONTRACT["test_specs"].items()
}
MATRIX_BUILD_TARGETS = tuple(CURRENT_MATRIX_CONTRACT["build_targets"])
MATRIX_BUILD_CACHE_KEYS = tuple(
    CURRENT_MATRIX_CONTRACT["build_cache_keys"])
MATRIX_BASE_BACKEND_FAILURE_TESTS = tuple(
    CURRENT_MATRIX_CONTRACT["base_backend_failure_tests"])
MATRIX_AVX512_BACKEND_FAILURE_TESTS = tuple(
    CURRENT_MATRIX_CONTRACT["avx512_backend_failure_tests"])
MATRIX_BACKEND_FAILURE_CTEST_REGEX = \
    CURRENT_MATRIX_CONTRACT["backend_failure_ctest_regex"]
MATRIX_PORTABLE_CTEST_REGEX = \
    CURRENT_MATRIX_CONTRACT["portable_ctest_regex"]
MATRIX_CUDA_CTEST_REGEX = CURRENT_MATRIX_CONTRACT["cuda_ctest_regex"]

PRE_XOR_MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS = {
    "Leopard2Backend.cpp": 2,
    "Leopard2BackendAVX2.cpp": 1,
    "Leopard2BackendSSSE3.cpp": 1,
    "Leopard2BackendScalar.cpp": 2,
    "Leopard2CpuFeatures.cpp": 2,
    "Leopard2Plan.cpp": 2,
    "LeopardCommon.cpp": 2,
    "LeopardFF16.cpp": 2,
    "LeopardFF8.cpp": 2,
    "leopard.cpp": 2,
    "leopard2.cpp": 2,
    "tests/leopard2/direct_oracle.cpp": 14,
    "tests/leopard2/direct_repair.cpp": 1,
    "tests/leopard2/fuzz_api.cpp": 1,
    "tests/leopard2/fuzz_pruned_transform.cpp": 1,
    "tests/leopard2/fuzz_replay.cpp": 2,
    "tests/leopard2/test_active_lch.cpp": 1,
    "tests/leopard2/test_api.cpp": 1,
    "tests/leopard2/test_arbitrary_counts_acceptance.cpp": 1,
    "tests/leopard2/test_backend_failures.cpp": 1,
    "tests/leopard2/test_backend_ops.cpp": 1,
    "tests/leopard2/test_batch_aliasing.cpp": 1,
    "tests/leopard2/test_concurrent_initialization.cpp": 1,
    "tests/leopard2/test_context_backends.cpp": 1,
    "tests/leopard2/test_r1_xor.cpp": 1,
    "tests/leopard2/test_boundaries.cpp": 1,
    "tests/leopard2/test_codec_options_abi.c": 1,
    "tests/leopard2/test_decode_high_acceptance.cpp": 1,
    "tests/leopard2/test_decode_low_acceptance.cpp": 1,
    "tests/leopard2/test_decode_plan_schedule.cpp": 1,
    "tests/leopard2/test_direct_encode.cpp": 1,
    "tests/leopard2/test_direct_oracle.cpp": 1,
    "tests/leopard2/test_direct_repair.cpp": 1,
    "tests/leopard2/test_encode_concurrency.cpp": 1,
    "tests/leopard2/test_encoder_gf16_legacy_matrix.cpp": 1,
    "tests/leopard2/test_field_options.cpp": 1,
    "tests/leopard2/test_gf16_padded_odd.cpp": 1,
    "tests/leopard2/test_gf16_tails.cpp": 1,
    "tests/leopard2/test_high_pruned_legacy.cpp": 1,
    "tests/leopard2/test_initialization_threads.cpp": 1,
    "tests/leopard2/test_legacy_golden.cpp": 1,
    "tests/leopard2/test_legacy_simd_init_failure.cpp": 1,
    "tests/leopard2/test_locator.cpp": 1,
    "tests/leopard2/test_low_gf16_direct_rows.cpp": 1,
    "tests/leopard2/test_max_counts.cpp": 1,
    "tests/leopard2/test_pruned_transform.cpp": 1,
    "tests/leopard2/test_public_api_contract.cpp": 1,
    "tests/leopard2/test_random.cpp": 1,
    "tests/leopard2/test_transform_differential.cpp": 1,
}

MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS = dict(
    CURRENT_MATRIX_CONTRACT["expected_compile_source_counts"])

CONFIGURATION_KEYS = (
    "CMAKE_BUILD_TYPE",
    "CMAKE_C_FLAGS",
    "CMAKE_C_FLAGS_RELEASE",
    "CMAKE_CXX_FLAGS",
    "CMAKE_CXX_FLAGS_RELEASE",
    "CMAKE_EXE_LINKER_FLAGS",
    "CMAKE_EXE_LINKER_FLAGS_RELEASE",
    "CMAKE_STATIC_LINKER_FLAGS",
    "CMAKE_STATIC_LINKER_FLAGS_RELEASE",
    "CMAKE_GENERATOR",
    "ENABLE_OPENMP",
    "LEO2_BACKEND_VARIANT",
    "LEO2_BUILD_BENCHMARKS",
    "LEO2_BUILD_FUZZERS",
    "LEO2_BUILD_TESTS",
    "LEO2_ENABLE_CUDA",
)

MATRIX_SOURCE_FILES = tuple(CURRENT_MATRIX_CONTRACT["source_files"])

TOOL_CACHE_KEYS = (
    ("cmake", "CMAKE_COMMAND"),
    ("cc", "CMAKE_C_COMPILER"),
    ("cxx", "CMAKE_CXX_COMPILER"),
    ("ar", "CMAKE_AR"),
    ("ranlib", "CMAKE_RANLIB"),
    ("cxx_ar", "CMAKE_CXX_COMPILER_AR"),
    ("cxx_ranlib", "CMAKE_CXX_COMPILER_RANLIB"),
    ("linker", "CMAKE_LINKER"),
    ("make", "CMAKE_MAKE_PROGRAM"),
)

PRE_XOR_CONFIGURED_TRANSLATION_UNITS = PRE_XOR_BUILD_TRANSLATION_UNITS + (
    "bench/leopard2/locator_benchmark.cpp",
    "tests/benchmark.cpp",
    "tests/experiments.cpp",
)

# Current v8 scopes the CMake File API proof to the benchmark and its complete
# target dependency graph.  Unrelated configured benchmarks cannot affect the
# linked executable and no longer make its evidence schema drift accidentally.
CONFIGURED_TRANSLATION_UNITS = BUILD_TRANSLATION_UNITS

HISTORICAL_CMAKE_IDENTITY = {
    "target": "libleopard",
    "archive": "liblibleopard.a",
    "target_directory": "libleopard.dir",
}
CANONICAL_CMAKE_IDENTITY = {
    "target": "leopard",
    "archive": "libleopard.a",
    "target_directory": "leopard.dir",
}
SCHEMA_TO_CMAKE_IDENTITY = {
    LEGACY_SCHEMA: HISTORICAL_CMAKE_IDENTITY,
    PRE_XOR_SCHEMA: CANONICAL_CMAKE_IDENTITY,
    SCHEMA: CANONICAL_CMAKE_IDENTITY,
}

HARDENED_SCHEMAS = frozenset((PRE_XOR_SCHEMA, SCHEMA))


def build_translation_units_for_schema(schema):
    return (BUILD_TRANSLATION_UNITS if schema == SCHEMA else
            PRE_XOR_BUILD_TRANSLATION_UNITS)


def configured_translation_units_for_schema(schema):
    return (CONFIGURED_TRANSLATION_UNITS if schema == SCHEMA else
            PRE_XOR_CONFIGURED_TRANSLATION_UNITS)


def archive_translation_units_for_schema(schema):
    if schema == SCHEMA:
        return CURRENT_ARCHIVE_TRANSLATION_UNITS
    if schema == PRE_XOR_SCHEMA:
        return PRE_XOR_ARCHIVE_TRANSLATION_UNITS
    require(schema == LEGACY_SCHEMA, "unsupported butterfly archive schema")
    return HISTORICAL_ARCHIVE_TRANSLATION_UNITS

PRE_XOR_RELEVANT_NON_LIBRARY_TARGETS = (
    "leopard2_backend_ssse3",
    "leopard2_backend_avx2",
    "bench_leopard2",
)
RELEVANT_NON_LIBRARY_TARGETS = (
    "leopard2_backend_ssse3",
    "leopard2_backend_avx2",
    "leopard2_backend_avx512",
    "bench_leopard2",
)


def cmake_identity_for_schema(schema):
    require(isinstance(schema, str),
            "butterfly manifest schema must be a string")
    identity = SCHEMA_TO_CMAKE_IDENTITY.get(schema)
    require(identity is not None, "unsupported butterfly manifest schema")
    return identity


def relevant_targets_for_schema(schema):
    non_library = (RELEVANT_NON_LIBRARY_TARGETS if schema == SCHEMA else
                   PRE_XOR_RELEVANT_NON_LIBRARY_TARGETS)
    return (cmake_identity_for_schema(schema)["target"],) + non_library

POWER_STATE_FILES = (
    ("scaling_driver", "cpu{cpu}/cpufreq/scaling_driver"),
    ("scaling_governor", "cpu{cpu}/cpufreq/scaling_governor"),
    ("energy_performance_preference",
     "cpu{cpu}/cpufreq/energy_performance_preference"),
    ("cpuinfo_min_freq_khz", "cpu{cpu}/cpufreq/cpuinfo_min_freq"),
    ("cpuinfo_max_freq_khz", "cpu{cpu}/cpufreq/cpuinfo_max_freq"),
    ("scaling_min_freq_khz", "cpu{cpu}/cpufreq/scaling_min_freq"),
    ("scaling_max_freq_khz", "cpu{cpu}/cpufreq/scaling_max_freq"),
    ("amd_pstate_status", "amd_pstate/status"),
    ("cpufreq_boost", "cpufreq/boost"),
    ("amd_pstate_no_turbo", "amd_pstate/no_turbo"),
    ("intel_pstate_no_turbo", "intel_pstate/no_turbo"),
)


class EvidenceError(Exception):
    pass


def require(condition, message):
    if not condition:
        raise EvidenceError(message)


def _linux_prctl(option, argument):
    require(sys.platform.startswith("linux"),
            "child descendant containment requires Linux")
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        prctl = libc.prctl
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            "Linux child-subreaper prctl is unavailable: {}".format(error))
    ctypes.set_errno(0)
    result = prctl(ctypes.c_int(option), argument,
                   ctypes.c_ulong(0), ctypes.c_ulong(0), ctypes.c_ulong(0))
    if result != 0:
        error_number = ctypes.get_errno()
        raise EvidenceError(
            "Linux child-subreaper prctl failed: " +
            os.strerror(error_number or 1))


def _get_child_subreaper():
    value = ctypes.c_int(-1)
    _linux_prctl(PR_GET_CHILD_SUBREAPER, ctypes.byref(value))
    require(value.value in (0, 1),
            "Linux returned an invalid child-subreaper state")
    return value.value


def _set_child_subreaper(value):
    require(value in (0, 1), "invalid child-subreaper state")
    _linux_prctl(PR_SET_CHILD_SUBREAPER, ctypes.c_ulong(value))
    require(_get_child_subreaper() == value,
            "Linux did not apply the requested child-subreaper state")


def _linux_pidfd_open(pid):
    try:
        pidfd_open = ctypes.CDLL(None, use_errno=True).pidfd_open
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            "Linux pidfd_open is unavailable: {}".format(error))
    ctypes.set_errno(0)
    descriptor = pidfd_open(ctypes.c_int(pid), ctypes.c_uint(0))
    if descriptor >= 0:
        return descriptor
    error_number = ctypes.get_errno()
    if error_number == errno.ESRCH:
        return None
    raise EvidenceError(
        "cannot open Linux pidfd for process {}: {}".format(
            pid, os.strerror(error_number or errno.EPERM)))


def _linux_pidfd_signal(descriptor, signal_number):
    try:
        send_signal = ctypes.CDLL(None, use_errno=True).pidfd_send_signal
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            "Linux pidfd_send_signal is unavailable: {}".format(error))
    ctypes.set_errno(0)
    result = send_signal(
        ctypes.c_int(descriptor), ctypes.c_int(signal_number), None,
        ctypes.c_uint(0))
    if result != 0:
        error_number = ctypes.get_errno()
        if error_number == errno.ESRCH:
            return
        raise EvidenceError(
            "cannot signal contained Linux process through pidfd: " +
            os.strerror(error_number or errno.EPERM))


def _validate_linux_pidfd_support():
    descriptor = _linux_pidfd_open(os.getpid())
    require(descriptor is not None,
            "Linux pidfd support cannot identify the runner process")
    try:
        _linux_pidfd_signal(descriptor, 0)
    finally:
        os.close(descriptor)


def _proc_process_record(pid):
    """Return (ppid, pgrp, session, starttime, state) from Linux procfs."""
    try:
        data = (Path("/proc") / str(pid) / "stat").read_bytes()
    except (FileNotFoundError, ProcessLookupError):
        return None
    except OSError as error:
        if error.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise EvidenceError(
            "cannot inspect Linux process {}: {}".format(pid, error))
    closing = data.rfind(b")")
    require(closing > 0 and closing + 2 < len(data),
            "Linux process {} has malformed procfs stat data".format(pid))
    fields = data[closing + 2:].split()
    require(len(fields) >= 20,
            "Linux process {} has truncated procfs stat data".format(pid))
    try:
        state = fields[0].decode("ascii")
        ppid = int(fields[1])
        pgrp = int(fields[2])
        session = int(fields[3])
        starttime = int(fields[19])
    except (UnicodeDecodeError, ValueError) as error:
        raise EvidenceError(
            "Linux process {} has invalid procfs stat fields: {}".format(
                pid, error))
    require(len(state) == 1 and ppid >= 0 and pgrp >= 0 and
            session >= 0 and starttime >= 0,
            "Linux process {} has invalid procfs process identity".format(pid))
    return ppid, pgrp, session, starttime, state


def _proc_process_snapshot():
    proc = Path("/proc")
    require(proc.is_dir() and (proc / "self/stat").is_file(),
            "child descendant containment requires mounted Linux procfs")
    try:
        names = os.listdir(str(proc))
    except OSError as error:
        raise EvidenceError("cannot enumerate Linux procfs: {}".format(error))
    result = {}
    for name in names:
        if not name.isascii() or not name.isdigit():
            continue
        pid = int(name)
        try:
            record = _proc_process_record(pid)
        except EvidenceError:
            try:
                owner = (proc / name).stat().st_uid
            except OSError:
                continue
            if owner == os.getuid():
                raise
            continue
        if record is not None:
            result[pid] = record
    self_record = _proc_process_record(os.getpid())
    require(self_record is not None and os.getpid() in result,
            "Linux procfs does not expose the runner process")
    return result


def _process_identity(pid, snapshot):
    record = snapshot.get(pid)
    return None if record is None else (pid, record[3])


def _emergency_linux_prctl(option, argument):
    if not sys.platform.startswith("linux"):
        raise EvidenceError("emergency child cleanup requires Linux")
    try:
        function = ctypes.CDLL(None, use_errno=True).prctl
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            "emergency child-subreaper prctl is unavailable: {}".format(error))
    ctypes.set_errno(0)
    result = function(ctypes.c_int(option), argument,
                      ctypes.c_ulong(0), ctypes.c_ulong(0), ctypes.c_ulong(0))
    if result != 0:
        number = ctypes.get_errno()
        raise EvidenceError(
            "emergency child-subreaper prctl failed: " +
            os.strerror(number or errno.EPERM))


def _emergency_get_child_subreaper():
    value = ctypes.c_int(-1)
    _emergency_linux_prctl(PR_GET_CHILD_SUBREAPER, ctypes.byref(value))
    if value.value not in (0, 1):
        raise EvidenceError("emergency child-subreaper state is invalid")
    return value.value


def _emergency_restore_child_subreaper(value):
    if value not in (0, 1):
        raise EvidenceError("emergency restore state is invalid")
    _emergency_linux_prctl(PR_SET_CHILD_SUBREAPER, ctypes.c_ulong(value))
    if _emergency_get_child_subreaper() != value:
        raise EvidenceError("emergency child-subreaper restore did not persist")


def _emergency_proc_process_record(pid):
    descriptor = None
    try:
        descriptor = os.open(
            "/proc/{}/stat".format(pid),
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))
        chunks = bytearray()
        while len(chunks) <= 65536:
            block = os.read(descriptor, min(4096, 65537 - len(chunks)))
            if not block:
                break
            chunks.extend(block)
        if len(chunks) > 65536:
            raise EvidenceError(
                "emergency procfs record {} is oversized".format(pid))
        data = bytes(chunks)
    except OSError as error:
        if error.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise EvidenceError(
            "cannot inspect emergency Linux process {}: {}".format(pid, error))
    finally:
        if descriptor is not None:
            os.close(descriptor)
    closing = data.rfind(b")")
    if closing <= 0 or closing + 2 >= len(data):
        raise EvidenceError(
            "emergency Linux process {} has malformed stat data".format(pid))
    fields = data[closing + 2:].split()
    if len(fields) < 20:
        raise EvidenceError(
            "emergency Linux process {} has truncated stat data".format(pid))
    try:
        state = fields[0].decode("ascii")
        ppid, pgrp, session = int(fields[1]), int(fields[2]), int(fields[3])
        starttime = int(fields[19])
    except (UnicodeDecodeError, ValueError) as error:
        raise EvidenceError(
            "emergency Linux process {} has invalid stat fields: {}".format(
                pid, error))
    if len(state) != 1 or min(ppid, pgrp, session, starttime) < 0:
        raise EvidenceError(
            "emergency Linux process {} identity is invalid".format(pid))
    return ppid, pgrp, session, starttime, state


def _emergency_proc_process_snapshot():
    try:
        names = os.listdir("/proc")
    except OSError as error:
        raise EvidenceError(
            "cannot enumerate emergency Linux procfs: {}".format(error))
    result = {}
    for name in names:
        if not name.isascii() or not name.isdigit():
            continue
        pid = int(name)
        try:
            record = _emergency_proc_process_record(pid)
        except EvidenceError:
            try:
                owner = os.stat(
                    "/proc/{}".format(name), follow_symlinks=False).st_uid
            except OSError:
                continue
            if owner == os.getuid():
                raise
            continue
        if record is not None:
            result[pid] = record
    if os.getpid() not in result:
        raise EvidenceError("emergency procfs does not expose the runner")
    return result


def _emergency_pidfd_open(pid):
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_open
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            "emergency pidfd_open is unavailable: {}".format(error))
    ctypes.set_errno(0)
    descriptor = function(ctypes.c_int(pid), ctypes.c_uint(0))
    if descriptor >= 0:
        return descriptor
    number = ctypes.get_errno()
    if number == errno.ESRCH:
        return None
    raise EvidenceError(
        "emergency pidfd_open failed for {}: {}".format(
            pid, os.strerror(number or errno.EPERM)))


def _emergency_pidfd_signal(descriptor, signal_number):
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_send_signal
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            "emergency pidfd_send_signal is unavailable: {}".format(error))
    ctypes.set_errno(0)
    result = function(ctypes.c_int(descriptor), ctypes.c_int(signal_number),
                      None, ctypes.c_uint(0))
    if result != 0:
        number = ctypes.get_errno()
        if number == errno.ESRCH:
            return
        raise EvidenceError(
            "emergency pidfd_send_signal failed: " +
            os.strerror(number or errno.EPERM))


def _emergency_signal_identity(identity):
    pid, starttime = identity
    record = _emergency_proc_process_record(pid)
    if record is None or record[3] != starttime:
        return
    descriptor = _emergency_pidfd_open(pid)
    if descriptor is None:
        return
    try:
        record = _emergency_proc_process_record(pid)
        if record is None or record[3] != starttime:
            return
        _emergency_pidfd_signal(descriptor, signal.SIGKILL)
    finally:
        os.close(descriptor)


class LinuxDescendantContainment(object):
    """Kill and reap descendants even after setsid() and a double fork."""

    def __init__(self):
        self.runner_pid = os.getpid()
        self.previous_subreaper = None
        self.baseline_children = set()
        self.leader = None
        self.known = set()
        self.spawned_process = None
        self.active = False
        self.proven_empty = False

    @staticmethod
    def _direct_children(snapshot, parent):
        return {(pid, record[3]) for pid, record in snapshot.items()
                if record[0] == parent}

    def __enter__(self):
        require(not self.active, "child descendant containment is already active")
        require(sys.platform.startswith("linux"),
                "child descendant containment requires Linux")
        task_root = Path("/proc/self/task")
        require(task_root.is_dir(),
                "child descendant containment requires Linux procfs task data")
        try:
            task_count = sum(1 for name in os.listdir(str(task_root))
                             if name.isascii() and name.isdigit())
        except OSError as error:
            raise EvidenceError(
                "cannot inspect runner threads for child containment: {}".format(
                    error))
        require(task_count == 1,
                "child descendant containment requires a single-threaded runner")
        _validate_linux_pidfd_support()
        self.previous_subreaper = _emergency_get_child_subreaper()
        require(_get_child_subreaper() == self.previous_subreaper,
                "normal and emergency child-subreaper reads disagree")
        try:
            _set_child_subreaper(1)
            snapshot = _proc_process_snapshot()
            self.baseline_children = self._direct_children(
                snapshot, self.runner_pid)
            require(not self.baseline_children,
                    "child descendant containment found pre-existing children")
            self.active = True
            return self
        except BaseException:
            if self.previous_subreaper is not None:
                _emergency_restore_child_subreaper(self.previous_subreaper)
            self.previous_subreaper = None
            raise

    def observe_spawn(self, process):
        require(self.active and self.spawned_process is process and process.pid > 0,
                "invalid emergency child observation")
        record = _emergency_proc_process_record(process.pid)
        if record is not None and record[0] == self.runner_pid:
            self.known.add((process.pid, record[3]))

    def attach(self, pid):
        require(self.active and self.leader is None and
                isinstance(pid, int) and pid > 0,
                "invalid child descendant containment attachment")
        record = _proc_process_record(pid)
        identity = None if record is None else (pid, record[3])
        require(identity is not None and record is not None and
                record[0] == self.runner_pid and
                identity not in self.baseline_children,
                "spawned process is not an owned direct child")
        self.leader = identity
        self.known.add(identity)

    def _discover(self, snapshot):
        targets = {identity for identity in self.known
                   if _process_identity(identity[0], snapshot) == identity}
        targets.update(
            identity for identity in self._direct_children(
                snapshot, self.runner_pid)
            if identity not in self.baseline_children)
        changed = True
        while changed:
            changed = False
            parent_pids = {pid for pid, _starttime in targets}
            for pid, record in snapshot.items():
                identity = (pid, record[3])
                if record[0] in parent_pids and identity not in targets:
                    targets.add(identity)
                    changed = True
        self.known.update(targets)
        return targets

    @staticmethod
    def _signal_identity(identity):
        pid, starttime = identity
        record = _proc_process_record(pid)
        if record is None or record[3] != starttime:
            return
        descriptor = _linux_pidfd_open(pid)
        if descriptor is None:
            return
        try:
            record = _proc_process_record(pid)
            if record is None or record[3] != starttime:
                return
            _linux_pidfd_signal(descriptor, signal.SIGKILL)
        finally:
            os.close(descriptor)

    def terminate_and_reap(self, process):
        require(self.active and self.leader is not None and
                self.leader[0] == process.pid,
                "child descendant containment is not attached to this process")
        deadline = time.monotonic() + CHILD_REAP_TIMEOUT_SECONDS
        empty_scans = 0
        while True:
            snapshot = _proc_process_snapshot()
            targets = self._discover(snapshot)
            for identity in sorted(targets, reverse=True):
                self._signal_identity(identity)
            if process.poll() is None:
                try:
                    process.wait(timeout=max(
                        0.001, min(0.05, deadline - time.monotonic())))
                except subprocess.TimeoutExpired:
                    pass
            for identity in sorted(self.known):
                if identity == self.leader:
                    continue
                record = _proc_process_record(identity[0])
                if record is None or record[3] != identity[1]:
                    continue
                try:
                    os.waitpid(identity[0], os.WNOHANG)
                except (ChildProcessError, ProcessLookupError):
                    pass
                except OSError as error:
                    raise EvidenceError(
                        "cannot reap contained child {}: {}".format(
                            identity[0], error))
            snapshot = _proc_process_snapshot()
            live_targets = self._discover(snapshot)
            if process.poll() is not None and not live_targets:
                empty_scans += 1
                if empty_scans >= 2:
                    self.proven_empty = True
                    return
            else:
                empty_scans = 0
            remaining = deadline - time.monotonic()
            require(remaining > 0,
                    "contained child descendants remained after SIGKILL")
            time.sleep(min(0.01, remaining))

    def _emergency_discover(self, snapshot):
        targets = {identity for identity in self.known
                   if _process_identity(identity[0], snapshot) == identity}
        targets.update(
            identity for identity in self._direct_children(
                snapshot, self.runner_pid)
            if identity not in self.baseline_children)
        changed = True
        while changed:
            changed = False
            parents = {pid for pid, _starttime in targets}
            for pid, record in snapshot.items():
                identity = (pid, record[3])
                if record[0] in parents and identity not in targets:
                    targets.add(identity)
                    changed = True
        self.known.update(targets)
        return targets

    def emergency_terminate_and_reap(self):
        process = self.spawned_process
        if process is None:
            return
        deadline = time.monotonic() + CHILD_REAP_TIMEOUT_SECONDS
        empty_scans = 0
        while True:
            snapshot = _emergency_proc_process_snapshot()
            record = snapshot.get(process.pid)
            if (record is not None and record[0] == self.runner_pid and
                    process.poll() is None):
                self.known.add((process.pid, record[3]))
            targets = self._emergency_discover(snapshot)
            for identity in sorted(targets, reverse=True):
                _emergency_signal_identity(identity)
            if process.poll() is None:
                try:
                    process.wait(timeout=max(
                        0.001, min(0.05, deadline - time.monotonic())))
                except subprocess.TimeoutExpired:
                    pass
            for identity in sorted(self.known):
                if identity[0] == process.pid:
                    continue
                record = _emergency_proc_process_record(identity[0])
                if record is None or record[3] != identity[1]:
                    continue
                try:
                    os.waitpid(identity[0], os.WNOHANG)
                except (ChildProcessError, ProcessLookupError):
                    pass
                except OSError as error:
                    raise EvidenceError(
                        "emergency child reap failed for {}: {}".format(
                            identity[0], error))
            snapshot = _emergency_proc_process_snapshot()
            live = self._emergency_discover(snapshot)
            if process.poll() is not None and not live:
                empty_scans += 1
                if empty_scans >= 2:
                    self.proven_empty = True
                    return
            else:
                empty_scans = 0
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise EvidenceError(
                    "emergency child descendants remained after bounded SIGKILL")
            time.sleep(min(0.01, remaining))

    def __exit__(self, exc_type, exc, tb):
        del tb
        if not self.active:
            return
        previous = self.previous_subreaper
        cleanup_error = None
        restore_error = None
        try:
            if (self.spawned_process is not None and
                    (exc_type is not None or not self.proven_empty)):
                try:
                    self.emergency_terminate_and_reap()
                except BaseException as error:
                    cleanup_error = error
            elif self.spawned_process is None:
                try:
                    snapshot = _emergency_proc_process_snapshot()
                    current = self._direct_children(snapshot, self.runner_pid)
                    if current != self.baseline_children:
                        cleanup_error = EvidenceError(
                            "unattached child appeared during containment")
                except BaseException as error:
                    cleanup_error = error
        finally:
            self.active = False
            self.previous_subreaper = None
            if previous is None:
                restore_error = EvidenceError(
                    "previous child-subreaper state was lost")
            else:
                try:
                    _emergency_restore_child_subreaper(previous)
                except BaseException as error:
                    restore_error = error
        if cleanup_error is not None or restore_error is not None:
            parts = []
            if cleanup_error is not None:
                parts.append("emergency cleanup failed: {}: {}".format(
                    type(cleanup_error).__name__, cleanup_error))
            if restore_error is not None:
                parts.append("subreaper restore failed: {}: {}".format(
                    type(restore_error).__name__, restore_error))
            if exc is not None:
                parts.append("primary failure: {}: {}".format(
                    type(exc).__name__, exc))
            raise EvidenceError("; ".join(parts))


def run_benchmark_bounded(argv, environment, timeout):
    """Capture one benchmark with bounded output and full-tree teardown."""
    require(isinstance(timeout, (int, float)) and
            not isinstance(timeout, bool) and math.isfinite(float(timeout)) and
            timeout > 0,
            "benchmark timeout must be positive and finite")
    process = None
    selector = selectors.DefaultSelector()
    stdout_fd = -1
    stderr_fd = -1
    outputs = {}
    returncode = -int(signal.SIGKILL)
    timed_out = False
    overflow = None
    try:
        with LinuxDescendantContainment() as containment:
            process = subprocess.Popen(
                [str(value) for value in argv], stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env=dict(environment), start_new_session=True)
            # Store the handle with no intervening injectable helper call.
            containment.spawned_process = process
            containment.observe_spawn(process)
            containment.attach(process.pid)
            require(process.stdout is not None and process.stderr is not None,
                    "cannot capture benchmark pipes")
            stdout_fd = process.stdout.fileno()
            stderr_fd = process.stderr.fileno()
            outputs = {stdout_fd: bytearray(), stderr_fd: bytearray()}
            limits = {stdout_fd: MAX_BENCHMARK_STDOUT_BYTES,
                      stderr_fd: MAX_BENCHMARK_STDERR_BYTES}
            for stream in (process.stdout, process.stderr):
                os.set_blocking(stream.fileno(), False)
                selector.register(stream, selectors.EVENT_READ)
            deadline = time.monotonic() + timeout
            try:
                while selector.get_map():
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        timed_out = True
                        break
                    for key, _events in selector.select(min(remaining, 0.1)):
                        descriptor = key.fileobj.fileno()
                        try:
                            block = os.read(descriptor, 65536)
                        except BlockingIOError:
                            continue
                        if not block:
                            selector.unregister(key.fileobj)
                            continue
                        outputs[descriptor].extend(block)
                        if len(outputs[descriptor]) > limits[descriptor]:
                            overflow = "benchmark output exceeded retained byte limit"
                            break
                    if overflow is not None:
                        break
                if not timed_out and overflow is None:
                    try:
                        returncode = process.wait(timeout=CHILD_REAP_TIMEOUT_SECONDS)
                    except subprocess.TimeoutExpired:
                        timed_out = True
            finally:
                # Run this even after returncode 0: a benchmark could otherwise
                # daemonize a child into a new session and escape the campaign.
                containment.terminate_and_reap(process)
                if isinstance(process.returncode, int):
                    returncode = process.returncode
    finally:
        selector.close()
        if process is not None:
            if process.stdout is not None:
                process.stdout.close()
            if process.stderr is not None:
                process.stderr.close()
    if overflow is not None:
        raise EvidenceError(overflow)
    if timed_out:
        raise EvidenceError("benchmark exceeded {} seconds".format(timeout))
    return subprocess.CompletedProcess(
        [str(value) for value in argv], returncode,
        bytes(outputs[stdout_fd]), bytes(outputs[stderr_fd]))


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False).encode("ascii")


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256_file(path):
    return sha256_bytes(Path(path).read_bytes())


def exact_utf8_file_content(path, label):
    path = Path(path)
    with path.open("rb") as input_file:
        raw = input_file.read(MAX_LINK_RECIPE_BYTES + 1)
    require(0 < len(raw) <= MAX_LINK_RECIPE_BYTES,
            "{} must contain 1..{} bytes".format(
                label, MAX_LINK_RECIPE_BYTES))
    try:
        text = raw.decode("utf-8", "strict")
    except UnicodeDecodeError as error:
        raise EvidenceError("{} is not strict UTF-8: {}".format(label, error))
    return exact_utf8_text_content(text, label)


def exact_utf8_text_content(text, label):
    require(isinstance(text, str), label + " text must be a string")
    try:
        raw = text.encode("utf-8", "strict")
    except UnicodeEncodeError as error:
        raise EvidenceError("{} text is not strict UTF-8: {}".format(
            label, error))
    require(0 < len(raw) <= MAX_LINK_RECIPE_BYTES and "\x00" not in text,
            "{} must contain 1..{} non-NUL UTF-8 bytes".format(
                label, MAX_LINK_RECIPE_BYTES))
    return {
        "encoding": "utf-8",
        "size": len(raw),
        "sha256": sha256_bytes(raw),
        "text": text,
    }


def parse_exact_recipe_content(content, expected_size, expected_sha256, label):
    require(isinstance(expected_size, int) and
            not isinstance(expected_size, bool) and
            0 < expected_size <= MAX_LINK_RECIPE_BYTES,
            label + " outer byte count")
    require(isinstance(content, dict) and set(content) == {
        "encoding", "size", "sha256", "text"},
        label + " content key set")
    require(content.get("encoding") == "utf-8" and
            isinstance(content.get("size"), int) and
            not isinstance(content.get("size"), bool) and
            0 < content["size"] <= MAX_LINK_RECIPE_BYTES and
            isinstance(content.get("text"), str),
            label + " content identity")
    try:
        raw = content["text"].encode("utf-8", "strict")
    except UnicodeEncodeError as error:
        raise EvidenceError("{} text is not strict UTF-8: {}".format(
            label, error))
    require("\x00" not in content["text"] and
            len(raw) == content["size"] == expected_size and
            re.fullmatch(r"[0-9a-f]{64}", content.get("sha256") or "") is
            not None and
            sha256_bytes(raw) == content["sha256"] == expected_sha256,
            label + " content size/SHA binding")
    try:
        lines = [shlex.split(value, posix=True)
                 for value in content["text"].splitlines() if value.strip()]
    except ValueError as error:
        raise EvidenceError("{} shell syntax: {}".format(label, error))
    require(lines and all(line and all(isinstance(value, str) and value
                                       for value in line) for line in lines),
            label + " parsed argv")
    return lines


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as output:
        json.dump(value, output, indent=2, sort_keys=True, ensure_ascii=True,
                  allow_nan=False)
        output.write("\n")
        output.flush()
        os.fsync(output.fileno())
    os.replace(str(temporary), str(path))


def checked_json_bytes(value, label):
    try:
        return canonical_bytes(value)
    except (TypeError, ValueError) as error:
        raise EvidenceError("{} is not canonical JSON: {}".format(label, error))


def parse_json_bytes(raw, label):
    def object_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise EvidenceError("{} has duplicate JSON key: {}".format(
                    label, key))
            result[key] = value
        return result

    def invalid_constant(value):
        raise EvidenceError("{} has non-finite JSON token: {}".format(
            label, value))

    try:
        return json.loads(raw.decode("utf-8"), object_pairs_hook=object_pairs,
                          parse_constant=invalid_constant)
    except UnicodeError as error:
        raise EvidenceError("{} is not UTF-8: {}".format(label, error))


def read_json(path, label):
    return parse_json_bytes(Path(path).read_bytes(), label)


def command_output(argv, cwd, label):
    completed = subprocess.run(
        argv, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False)
    require(completed.returncode == 0,
            "{} failed rc={}: {}".format(
                label, completed.returncode,
                completed.stderr.decode("utf-8", "replace").strip()))
    return completed.stdout


def git_blob(repo, commit, relative):
    return command_output(
        ["git", "show", "{}:{}".format(commit, relative)], repo,
        "git show " + relative)


def git_record(repo, commit, require_clean_checkout=False):
    require(re.fullmatch(r"[0-9a-f]{40}", commit or "") is not None,
            "commit must be a full lowercase SHA-1")
    command_output(["git", "cat-file", "-e", commit + "^{commit}"], repo,
                   "git commit lookup")
    tree = command_output(
        ["git", "rev-parse", commit + "^{tree}"], repo,
        "git tree lookup").decode("ascii").strip()
    require(re.fullmatch(r"[0-9a-f]{40}", tree) is not None,
            "git tree is malformed")
    tree_listing = command_output(
        ["git", "ls-tree", "-r", "-z", "--full-tree", commit], repo,
        "git tree manifest")
    tree_entries = [value for value in tree_listing.split(b"\0") if value]
    require(tree_entries, "git tree manifest is empty")
    if require_clean_checkout:
        head = command_output(
            ["git", "rev-parse", "HEAD"], repo,
            "git HEAD lookup").decode("ascii").strip()
        require(head == commit, "source checkout HEAD differs from declared commit")
        status = command_output(
            ["git", "status", "--porcelain=v1", "-z",
             "--untracked-files=all", "--ignore-submodules=none"], repo,
            "git status")
        require(status == b"",
                "source checkout has tracked, untracked, or submodule changes")
    return {
        "commit": commit,
        "tree": tree,
        "tree_entry_count": len(tree_entries),
        "tree_manifest_sha256": sha256_bytes(tree_listing),
    }


def validate_git_record(repo, record):
    require(set(record) == {
        "commit", "tree", "tree_entry_count", "tree_manifest_sha256"},
        "git record keys")
    actual = git_record(repo, record["commit"], False)
    require(actual == record, "git commit/tree binding mismatch")


def parse_cpu_list(text):
    cpus = set()
    for part in text.strip().split(","):
        if not part:
            continue
        if "-" in part:
            first, last = (int(value) for value in part.split("-", 1))
            require(first <= last, "descending CPU range")
            cpus.update(range(first, last + 1))
        else:
            cpus.add(int(part))
    return sorted(cpus)


def sibling_topology(cpu):
    path = Path("/sys/devices/system/cpu/cpu{}/topology/thread_siblings_list".format(cpu))
    try:
        raw = path.read_text(encoding="ascii").strip()
        parsed = parse_cpu_list(raw)
    except (OSError, ValueError) as error:
        raise EvidenceError("cannot read CPU sibling topology: {}".format(error))
    require(parsed, "empty CPU sibling topology")
    return {"raw": raw, "sha256": sha256_bytes(raw.encode("ascii")),
            "cpus": parsed}


def optional_text(path):
    try:
        value = Path(path).read_text(encoding="ascii").strip()
    except (OSError, UnicodeError):
        return None
    return value if value else None


def cpu_identity(cpu):
    raw = Path("/proc/cpuinfo").read_text(
        encoding="utf-8", errors="replace")
    selected = None
    for block in raw.split("\n\n"):
        fields = {}
        for line in block.splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                fields[key.strip()] = value.strip()
        if fields.get("processor") == str(cpu):
            selected = fields
            break
    require(selected is not None,
            "benchmark CPU is absent from /proc/cpuinfo")
    aliases = (
        ("vendor_id", "vendor_id"),
        ("model_name", "model name"),
        ("cpu_family", "cpu family"),
        ("model", "model"),
        ("stepping", "stepping"),
        ("microcode", "microcode"),
    )
    result = {logical: selected.get(source) for logical, source in aliases}
    require(all(isinstance(value, str) and value for value in result.values()),
            "x86 CPU identity is incomplete")
    result["processor"] = cpu
    return result


def power_state(cpu):
    root = Path("/sys/devices/system/cpu")
    return {key: optional_text(root / relative.format(cpu=cpu))
            for key, relative in POWER_STATE_FILES}


def current_frequency(cpu):
    return optional_text(Path(
        "/sys/devices/system/cpu/cpu{}/cpufreq/scaling_cur_freq".format(cpu)))


def host_record(cpu, pre_frequency, post_frequency, static_power):
    uname = os.uname()
    payload = {
        "cpu": cpu_identity(cpu),
        "platform": {
            "platform": platform.platform(),
            "uname": {"sysname": uname.sysname, "nodename": uname.nodename,
                      "release": uname.release, "version": uname.version,
                      "machine": uname.machine},
        },
        "power": {
            "state": static_power,
            "frequency_snapshots_khz": {
                "pre_campaign": pre_frequency,
                "post_campaign": post_frequency,
            },
        },
    }
    return {"record": payload, "sha256": sha256_bytes(canonical_bytes(payload))}


def validate_host_record(host, cpu):
    require(set(host) == {"record", "sha256"}, "host record keys")
    record = host["record"]
    require(sha256_bytes(canonical_bytes(record)) == host["sha256"],
            "host record digest")
    require(set(record) == {"cpu", "platform", "power"},
            "host payload keys")
    identity = record["cpu"]
    require(set(identity) == {"processor", "vendor_id", "model_name",
                             "cpu_family", "model", "stepping", "microcode"} and
            identity["processor"] == cpu and
            all(isinstance(identity[key], str) and identity[key]
                for key in set(identity) - {"processor"}),
            "host CPU identity")
    platform_record = record["platform"]
    require(set(platform_record) == {"platform", "uname"} and
            isinstance(platform_record["platform"], str) and
            platform_record["platform"], "host platform identity")
    uname = platform_record["uname"]
    require(set(uname) == {"sysname", "nodename", "release", "version", "machine"} and
            all(isinstance(value, str) and value for value in uname.values()),
            "host uname identity")
    power = record["power"]
    require(set(power) == {"state", "frequency_snapshots_khz"},
            "host power keys")
    state = power["state"]
    require(set(state) == {key for key, _ in POWER_STATE_FILES} and
            all(value is None or (isinstance(value, str) and value)
                for value in state.values()), "host power-state identity")
    for key in ("cpuinfo_min_freq_khz", "cpuinfo_max_freq_khz",
                "scaling_min_freq_khz", "scaling_max_freq_khz"):
        require(state[key] is None or state[key].isdigit(),
                "host frequency bound: " + key)
    snapshots = power["frequency_snapshots_khz"]
    require(set(snapshots) == {"pre_campaign", "post_campaign"} and
            all(value is None or (isinstance(value, str) and value.isdigit())
                for value in snapshots.values()),
            "host current-frequency snapshots")


def benchmark_arguments(item, backend):
    require(backend in SUPPORTED_BACKENDS,
            "unsupported campaign backend: " + str(backend))
    return [
        "--k", str(item["K"]), "--r", str(item["R"]),
        "--profile", item["profile"], "--field", item["field"],
        "--bytes", str(item["bytes"]), "--loss", str(item["loss"]),
        "--batch", "1", "--reuse", "8", "--iterations", "7",
        "--warmup", "3", "--threads", "1", "--backend", backend,
        "--seed", "42", "--json", "-",
    ]


def expected_jobs():
    jobs = []
    for item in CELLS:
        for round_number in ROUNDS:
            for sequence, build in SEQUENCES:
                name = "{}-r{}-{}-{}".format(
                    item["name"], round_number, sequence, build)
                jobs.append((name, item, round_number, sequence, build))
    return jobs


def logical_executable(build, digest):
    return "@{}:sha256:{}".format(build, digest)


def finite_number(value, label):
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            label + " must be numeric")
    result = float(value)
    require(math.isfinite(result), label + " must be finite")
    return result


def validate_metric(metric, label, rate_keys):
    required = {
        "mad_us_per_batch_call", "maximum_us_per_batch_call",
        "median_us_per_batch_call", "minimum_us_per_batch_call",
    } | set(rate_keys)
    require(set(metric) == required, label + " metric keys")
    values = {key: finite_number(metric[key], label + " " + key)
              for key in required}
    require(values["minimum_us_per_batch_call"] > 0.0,
            label + " minimum must be positive")
    require(values["minimum_us_per_batch_call"] <=
            values["median_us_per_batch_call"] <=
            values["maximum_us_per_batch_call"], label + " ordering")
    require(values["mad_us_per_batch_call"] >= 0.0, label + " MAD")
    for key in rate_keys:
        require(finite_number(metric[key], label + " " + key) > 0.0,
                label + " nonpositive rate")
    return {key: values[key] for key in (
        "minimum_us_per_batch_call", "median_us_per_batch_call",
        "maximum_us_per_batch_call", "mad_us_per_batch_call")}


def validate_setup_metric(metric, label):
    keys = {"median_us", "mad_us", "minimum_us", "maximum_us"}
    require(isinstance(metric, dict) and set(metric) == keys,
            label + " setup metric keys")
    values = {key: finite_number(metric[key], label + " " + key)
              for key in keys}
    require(0.0 <= values["minimum_us"] <= values["median_us"] <=
            values["maximum_us"] and values["mad_us"] >= 0.0,
            label + " setup metric ordering")


def check_raw(
        raw, item, label, requested_backend, missing_indices=None,
        evidence_schema=SCHEMA):
    require(requested_backend in SUPPORTED_BACKENDS,
            label + " unsupported requested backend")
    require(evidence_schema in SCHEMA_TO_CMAKE_IDENTITY,
            label + " unsupported evidence schema")
    require(set(raw) == {"schema", "build", "parameters", "resolved",
                         "correctness", "memory", "metrics", "legacy"},
            label + " top-level key set")
    require(raw.get("schema") == "leopard2-benchmark-v1", label + " schema")
    build = raw.get("build")
    require(isinstance(build, dict) and
            set(build) == {"compiler", "compiler_version", "cplusplus"} and
            isinstance(build["compiler"], str) and build["compiler"] and
            isinstance(build["compiler_version"], str) and build["compiler_version"] and
            isinstance(build["cplusplus"], int) and not isinstance(build["cplusplus"], bool) and
            build["cplusplus"] >= 201103,
            label + " build identity")
    parameters = raw.get("parameters")
    require(isinstance(parameters, dict), label + " parameters")
    required_parameter_keys = {
        "K", "R", "requested_profile", "requested_field",
        "requested_backend", "force_generic_decode",
        "force_specialized_decode", "shard_bytes", "loss_count",
        "missing_original_indices", "batch", "reuse", "iterations",
        "warmup", "thread_count", "seed",
    }
    if evidence_schema in HARDENED_SCHEMAS:
        required_parameter_keys.update({
            "force_tiled_decode", "force_materialized_decode"})
    require(set(parameters) == required_parameter_keys,
            label + " parameter key set")
    expected = {
        "K": item["K"], "R": item["R"],
        "requested_profile": item["resolved_profile"],
        "requested_field": item["field"],
        "requested_backend": requested_backend,
        "force_generic_decode": False, "force_specialized_decode": False,
        "shard_bytes": item["bytes"], "loss_count": item["loss"],
        "batch": 1, "reuse": 8, "iterations": 7, "warmup": 3,
        "thread_count": 1, "seed": 42,
    }
    if evidence_schema in HARDENED_SCHEMAS:
        expected.update({
            "force_tiled_decode": False,
            "force_materialized_decode": False,
        })
    for key, expected_value in expected.items():
        actual_value = parameters[key]
        matches = (
            type(actual_value) is bool and actual_value is expected_value
            if type(expected_value) is bool else
            actual_value == expected_value)
        require(matches,
                "{} parameter {} expected {!r}, got {!r}".format(
                    label, key, expected_value, actual_value))
    missing = parameters["missing_original_indices"]
    require(isinstance(missing, list) and len(missing) == item["loss"],
            label + " missing-index count")
    require(missing == sorted(set(missing)) and
            all(isinstance(index, int) and 0 <= index < item["K"] for index in missing),
            label + " missing-index geometry")
    if missing_indices is not None:
        require(missing == missing_indices, label + " missing indices changed")

    resolved = raw.get("resolved")
    require(isinstance(resolved, dict), label + " resolved")
    expected_resolved = {
        "profile": item["resolved_profile"], "field": item["field"],
        "backend": requested_backend, "thread_count": 1,
        "parent_count": item["parent_count"],
        "padded_side": item["padded_side"],
    }
    require(resolved == expected_resolved, label + " resolved identity")
    correctness = raw.get("correctness")
    legacy_comparison = None
    if item["profile"] == "high" and item["bytes"] % 64 == 0:
        legacy_comparison = "matched"
    expected_correctness = {
        "leopard2_round_trip": True,
        "legacy_comparison": legacy_comparison,
    }
    require(correctness == expected_correctness, label + " correctness identity")
    memory = raw.get("memory")
    memory_keys = {"scratch_alignment", "encode_scratch_bytes_per_stripe",
                   "decode_scratch_bytes_per_stripe", "encode_scratch_bytes_batch",
                   "decode_scratch_bytes_batch"}
    require(isinstance(memory, dict) and set(memory) == memory_keys and
            all(isinstance(memory[key], int) and not isinstance(memory[key], bool) and
                memory[key] >= 0 for key in memory_keys),
            label + " memory identity")
    alignment = memory["scratch_alignment"]
    require(alignment > 0 and alignment & (alignment - 1) == 0 and
            memory["encode_scratch_bytes_batch"] ==
            memory["encode_scratch_bytes_per_stripe"] and
            memory["decode_scratch_bytes_batch"] ==
            memory["decode_scratch_bytes_per_stripe"],
            label + " scratch geometry")
    metrics = raw.get("metrics")
    require(isinstance(metrics, dict) and set(metrics) == {
        "codec_setup", "encode_execution", "decode_plan_setup",
        "decode_execution", "decode_amortized_at_reuse", "rate_semantics"},
        label + " metrics")
    validate_setup_metric(metrics["codec_setup"], label + " codec setup")
    validate_setup_metric(metrics["decode_plan_setup"], label + " plan setup")
    encode = validate_metric(
        metrics["encode_execution"], label + " encode",
        ("input_GB_per_s", "parity_output_GB_per_s"))
    decode = validate_metric(
        metrics["decode_execution"], label + " decode",
        ("offered_received_GB_per_s", "repaired_output_GB_per_s"))
    amortized = metrics["decode_amortized_at_reuse"]
    require(isinstance(amortized, dict) and set(amortized) == {
        "reuse_count", "derived_median_us_per_batch_call",
        "offered_received_GB_per_s", "repaired_output_GB_per_s"} and
        amortized["reuse_count"] == 8 and
        all(finite_number(amortized[key], label + " amortized " + key) > 0.0
            for key in ("derived_median_us_per_batch_call",
                        "offered_received_GB_per_s", "repaired_output_GB_per_s")),
        label + " amortized decode")
    require(metrics["rate_semantics"] ==
            "offered_received counts all non-null shard pointers supplied; "
            "a plan may read a deterministic subset",
            label + " rate semantics")
    legacy = raw.get("legacy")
    require(isinstance(legacy, dict) and set(legacy) == {
        "available", "unavailable_reason", "codec_setup",
        "decode_timing_includes_setup", "encode_execution",
        "decode_including_setup"}, label + " legacy keys")
    legacy_available = legacy_comparison == "matched"
    require(legacy["available"] is legacy_available and
            legacy["codec_setup"] is None and
            legacy["decode_timing_includes_setup"] is True,
            label + " legacy availability")
    if legacy_available:
        require(legacy["unavailable_reason"] is None,
                label + " legacy reason")
        validate_metric(legacy["encode_execution"], label + " legacy encode",
                        ("input_GB_per_s", "parity_output_GB_per_s"))
        validate_metric(legacy["decode_including_setup"], label + " legacy decode",
                        ("offered_received_GB_per_s", "repaired_output_GB_per_s"))
    else:
        expected_reason = (
            "old Leopard only defines the legacy high wire profile"
            if item["profile"] != "high" else
            "old Leopard requires shard bytes divisible by 64")
        require(legacy["unavailable_reason"] == expected_reason and
                legacy["encode_execution"] is None and
                legacy["decode_including_setup"] is None,
                label + " legacy unavailable semantics")
    return encode, decode, missing, build


def invocation_log_standard_error(metric):
    """Conservative robust log-scale SE from one seven-sample invocation.

    The benchmark does not retain the seven underlying samples in v6, so this
    deliberately uses both the reported MAD and full min/max span.  It is not a
    claim that the original samples can be reconstructed.
    """
    median = metric["median_us_per_batch_call"]
    mad = metric["mad_us_per_batch_call"]
    span = (metric["maximum_us_per_batch_call"] -
            metric["minimum_us_per_batch_call"])
    mad_se = 1.4826 * mad / (math.sqrt(INVOCATION_SAMPLES) * median)
    span_se = span / (2.0 * math.sqrt(INVOCATION_SAMPLES) * median)
    return max(mad_se, span_se)


def paired_abba_statistics(entries, raw_by_name, item, metric_index):
    invocation = {}
    for entry in entries:
        if entry["cell"] == item["name"]:
            invocation[(entry["round"], entry["sequence"])] = \
                raw_by_name[entry["name"]][metric_index]
    expected = {(round_number, sequence) for round_number in ROUNDS
                for sequence, _ in SEQUENCES}
    require(set(invocation) == expected,
            item["name"] + " paired ABBA invocation geometry")

    round_logs = []
    round_within_variances = []
    round_speedups = []
    for round_number in ROUNDS:
        a1 = invocation[(round_number, "A1")]
        b1 = invocation[(round_number, "B1")]
        b2 = invocation[(round_number, "B2")]
        a2 = invocation[(round_number, "A2")]
        log_ratio = 0.5 * (
            math.log(a1["median_us_per_batch_call"]) +
            math.log(a2["median_us_per_batch_call"]) -
            math.log(b1["median_us_per_batch_call"]) -
            math.log(b2["median_us_per_batch_call"]))
        round_logs.append(log_ratio)
        round_speedups.append((math.exp(log_ratio) - 1.0) * 100.0)
        # Each round is the difference of two two-invocation log means.
        round_within_variances.append(sum(
            invocation_log_standard_error(value) ** 2
            for value in (a1, a2, b1, b2)) / 4.0)

    mean_log = statistics.mean(round_logs)
    between_mean_variance = statistics.variance(round_logs) / len(round_logs)
    within_mean_variance = sum(round_within_variances) / len(round_logs) ** 2
    standard_error = math.sqrt(between_mean_variance + within_mean_variance)
    lower_log = mean_log - ONE_SIDED_T95 * standard_error
    return {
        "method": "paired ABBA log-ratio; one-sided 95% t lower bound, df=15",
        "round_speedup_percent": round_speedups,
        "speedup_percent": (math.exp(mean_log) - 1.0) * 100.0,
        "lower_confidence_speedup_percent": (math.exp(lower_log) - 1.0) * 100.0,
        "log_standard_error": standard_error,
        "confidence_level": 0.95,
    }


def summarize(entries, raw_by_name, requested_backend):
    require(requested_backend in SUPPORTED_BACKENDS,
            "unsupported summary backend: " + str(requested_backend))
    summary = []
    for item in CELLS:
        cell_result = {
            "name": item["name"],
            "kind": item["kind"],
            "minimum_speedup_percent": item["minimum_speedup_percent"],
            "parameters": {
                "K": item["K"], "R": item["R"],
                "profile": item["resolved_profile"], "field": item["field"],
                "requested_backend": requested_backend,
                "resolved_backend": requested_backend,
                "shard_bytes": item["bytes"], "losses": item["loss"],
                "rounded_transform_bytes": item["rounded_transform_bytes"],
                "tail_mapping": item["tail_mapping"],
            },
        }
        for metric_index, metric_name in ((0, "encode"), (1, "decode")):
            result = {}
            for build in ("baseline", "candidate"):
                metrics = [raw_by_name[entry["name"]][metric_index]
                          for entry in entries
                          if entry["cell"] == item["name"] and
                          entry["build"] == build]
                expected_invocations = 2 * len(ROUNDS)
                require(len(metrics) == expected_invocations,
                        "{} {} {} requires {} invocations".format(
                            item["name"], metric_name, build,
                            expected_invocations))
                values = [value["median_us_per_batch_call"] for value in metrics]
                median = statistics.median(values)
                mad = statistics.median(abs(value - median) for value in values)
                result[build + "_median_us"] = median
                result[build + "_mad_us"] = mad
                result[build + "_minimum_us"] = min(values)
                result[build + "_maximum_us"] = max(values)
            paired = paired_abba_statistics(
                entries, raw_by_name, item, metric_index)
            result.update(paired)
            result["minimum_speedup_percent"] = item["minimum_speedup_percent"]
            result["accepted"] = (
                paired["lower_confidence_speedup_percent"] >=
                item["minimum_speedup_percent"])
            cell_result[metric_name] = result
        summary.append(cell_result)
    return summary


def stable_raw_digest(entries):
    evidence = [{
        "name": entry["name"],
        "command_sha256": entry["command_sha256"],
        "binary_sha256": entry["binary_sha256"],
        "stdout_sha256": entry["stdout_sha256"],
        "stderr_sha256": entry["stderr_sha256"],
        "returncode": entry["returncode"],
    } for entry in sorted(entries, key=lambda value: value["name"])]
    return sha256_bytes(canonical_bytes(evidence))


def cache_values(cache_path):
    result = {}
    for line in Path(cache_path).read_text(
            encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"([^:#=]+):[^=]+=(.*)$", line)
        if match:
            result[match.group(1)] = match.group(2)
    return result


def executable_identity(path, logical_name, cache_key):
    executable = Path(path).resolve()
    require(executable.is_file(), "cached tool is missing: " + cache_key)
    completed = subprocess.run(
        [str(executable), "--version"], cwd=str(executable.parent),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    require(completed.returncode == 0,
            "{} --version failed rc={}".format(cache_key, completed.returncode))
    raw_version = completed.stdout + completed.stderr
    require(raw_version, cache_key + " produced an empty version identity")
    version = raw_version.decode("utf-8", "replace")
    return executable, {
        "logical_name": logical_name,
        "cache_key": cache_key,
        "basename": executable.name,
        "binary_sha256": sha256_file(executable),
        "version": version,
        "version_sha256": sha256_bytes(version.encode("utf-8")),
    }


def cache_tools(cache_path):
    values = cache_values(cache_path)
    paths = {}
    records = {}
    for logical_name, cache_key in TOOL_CACHE_KEYS:
        raw = values.get(cache_key)
        require(raw, "CMake cache omits tool: " + cache_key)
        paths[logical_name], records[logical_name] = executable_identity(
            raw, logical_name, cache_key)
    return paths, records


def relative_to(path, root):
    try:
        return path.resolve().relative_to(root.resolve())
    except ValueError:
        return None


def tagged_path(path, source_root, build_root):
    path = path.resolve()
    relative = relative_to(path, build_root)
    if relative is not None:
        return "@build/" + relative.as_posix()
    relative = relative_to(path, source_root)
    if relative is not None:
        return "@source/" + relative.as_posix()
    return "@external/{}/{}".format(
        sha256_bytes(str(path).encode("utf-8")), path.name)


def command_executable(value, directory):
    candidate = Path(value)
    if candidate.is_absolute() or "/" in value:
        if not candidate.is_absolute():
            candidate = directory / candidate
        return candidate.resolve()
    found = shutil.which(value)
    require(found is not None, "compile command executable not found: " + value)
    return Path(found).resolve()


def command_output_path(command, directory):
    outputs = []
    for index, value in enumerate(command):
        if value == "-o":
            require(index + 1 < len(command), "compile command has bare -o")
            outputs.append(command[index + 1])
        elif value.startswith("-o") and len(value) > 2:
            outputs.append(value[2:])
    require(len(outputs) == 1, "compile command must contain exactly one output")
    output = Path(outputs[0])
    if not output.is_absolute():
        output = directory / output
    return output.resolve()


def normalize_command(command, source_root, build_root, tool_paths):
    replacements = [
        (str(source_root.resolve()), "@source"),
        (str(build_root.resolve()), "@build"),
    ]
    for logical_name, path in tool_paths.items():
        # Preserve both the spelling supplied to CMake (often a /usr/bin
        # symlink) and the canonical executable identity used for comparison.
        for spelling in {str(path), str(path.resolve())}:
            replacements.append((spelling, "@tool/" + logical_name))
    replacements.sort(key=lambda value: len(value[0]), reverse=True)
    normalized = []
    for value in command:
        result = value
        for actual, logical in replacements:
            result = result.replace(actual, logical)
        normalized.append(result)
    return normalized


def normalized_compile_entry(entry, source_root, build_root, tool_paths):
    require(isinstance(entry, dict), "compile command entry must be an object")
    require("file" in entry, "compile command entry omits file")
    source = Path(entry["file"])
    if not source.is_absolute():
        source = Path(entry.get("directory", source_root)) / source
    source = source.resolve()
    try:
        relative = source.relative_to(source_root.resolve()).as_posix()
    except ValueError:
        raise EvidenceError("compile source escapes declared source root: " + str(source))
    require(source.is_file(), "compile source is missing: " + relative)
    directory = Path(entry.get("directory", build_root)).resolve()
    require(relative_to(directory, build_root) is not None,
            "compile directory escapes declared build root: " + relative)
    if "arguments" in entry:
        command = [str(value) for value in entry["arguments"]]
    else:
        command = shlex.split(entry.get("command", ""))
    require(command, "empty compile command for " + relative)
    require(command_executable(command[0], directory) == tool_paths["cxx"].resolve(),
            "compile command/compiler mismatch: " + relative)
    source_arguments = []
    for value in command[1:]:
        if value.startswith("-"):
            continue
        candidate = Path(value)
        if not candidate.is_absolute():
            candidate = directory / candidate
        if candidate.resolve() == source:
            source_arguments.append(value)
    require(len(source_arguments) == 1,
            "compile command is not bound to declared source: " + relative)
    output = command_output_path(command, directory)
    require(relative_to(output, build_root) is not None,
            "compile output escapes declared build root: " + relative)
    normalized_argv = normalize_command(
        command, source_root, build_root, tool_paths)
    # compile_commands.json may retain a symlink spelling such as /usr/bin/c++
    # even though CMake canonicalized the compiler cache entry.  Resolution was
    # checked above; store the stable logical tool identity rather than either
    # host-specific spelling.
    normalized_argv[0] = "@tool/cxx"
    normalized = {
        "file": relative,
        "argv": normalized_argv,
        "directory": tagged_path(directory, source_root, build_root),
        "output": tagged_path(output, source_root, build_root),
        # CMake exports commands for configured targets that this focused
        # evidence build intentionally does not materialize.  The production
        # target subset is checked below and must have a concrete digest.
        "output_sha256": sha256_file(output) if output.is_file() else None,
    }
    normalized["command_sha256"] = sha256_bytes(canonical_bytes(normalized))
    return relative, (normalized, output, directory)


def parse_dependency_file(path):
    raw = path.read_text(encoding="utf-8", errors="surrogateescape")
    unfolded = raw.replace("\\\n", " ")
    require(":" in unfolded, "dependency file omits target: " + str(path))
    dependencies = shlex.split(unfolded.split(":", 1)[1], posix=True)
    require(dependencies, "dependency file has no inputs: " + str(path))
    return dependencies


def dependency_closure(compile_entries, source_root, build_root):
    dependencies = {}
    source_files = set()
    by_translation_unit = {}
    for relative, (_, output, directory) in compile_entries.items():
        depfile = Path(str(output) + ".d")
        require(depfile.is_file(),
                "compiler dependency file missing for " + relative)
        unit_dependencies = set()
        for value in parse_dependency_file(depfile):
            path = Path(value)
            if not path.is_absolute():
                path = directory / path
            path = path.resolve()
            require(path.is_file(), "dependency input missing: " + str(path))
            tag = tagged_path(path, source_root, build_root)
            digest = sha256_file(path)
            require(tag not in dependencies or dependencies[tag] == digest,
                    "dependency path collision: " + tag)
            dependencies[tag] = digest
            unit_dependencies.add(tag)
            if relative_to(path, source_root) is not None:
                source_files.add(tag)
        expected_source = "@source/" + relative
        require(expected_source in dependencies,
                "dependency closure omits source: " + relative)
        by_translation_unit[relative] = sorted(unit_dependencies)
    require(dependencies, "empty compiler dependency closure")
    manifest = [{"path": path, "sha256": dependencies[path]}
                for path in sorted(dependencies)]
    return {
        "file_count": len(manifest),
        "source_file_count": len(source_files),
        "manifest_sha256": sha256_bytes(canonical_bytes(manifest)),
        "manifest": manifest,
        "by_translation_unit": by_translation_unit,
    }


def command_evidence(argv, cwd, label, logical_argv, environment=None):
    completed = subprocess.run(
        [str(value) for value in argv], cwd=str(cwd), stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False, env=environment)
    require(completed.returncode == 0,
            "{} failed rc={}: {}".format(
                label, completed.returncode,
                completed.stderr.decode("utf-8", "replace").strip()))
    return {
        "argv": logical_argv,
        "command_sha256": sha256_bytes(canonical_bytes(logical_argv)),
        "stdout_sha256": sha256_bytes(completed.stdout),
        "stderr_sha256": sha256_bytes(completed.stderr),
        "returncode": completed.returncode,
    }


def file_api_targets(source_root, build_root, schema=SCHEMA):
    cmake_identity = cmake_identity_for_schema(schema)
    relevant_targets = relevant_targets_for_schema(schema)
    configured_translation_units = configured_translation_units_for_schema(
        schema)
    build_translation_units = build_translation_units_for_schema(schema)
    reply = build_root / ".cmake/api/v1/reply"
    indexes = sorted(reply.glob("index-*.json"))
    require(len(indexes) == 1, "CMake File API must produce exactly one index")
    index = read_json(indexes[0], "CMake File API index")
    response = index.get("reply", {}).get("codemodel-v2")
    require(isinstance(response, dict) and isinstance(response.get("jsonFile"), str),
            "CMake File API omits codemodel-v2")
    codemodel_path = reply / response["jsonFile"]
    codemodel = read_json(codemodel_path, "CMake File API codemodel")
    configurations = codemodel.get("configurations")
    require(isinstance(configurations, list) and len(configurations) == 1,
            "evidence build requires one CMake configuration")
    target_refs = configurations[0].get("targets")
    require(isinstance(target_refs, list), "CMake codemodel target list")
    id_to_name = {value.get("id"): value.get("name") for value in target_refs}
    documents = {}
    for reference in target_refs:
        name = reference.get("name")
        filename = reference.get("jsonFile")
        require(isinstance(name, str) and name and isinstance(filename, str),
                "malformed CMake target reference")
        require(name not in documents, "duplicate CMake target: " + name)
        documents[name] = read_json(reply / filename, "CMake target " + name)
    require(set(relevant_targets) <= set(documents),
            "CMake target graph omits evidence target")

    configured_units = set()
    configured_documents = (
        [documents[name] for name in relevant_targets]
        if schema == SCHEMA else documents.values())
    for document in configured_documents:
        for source in document.get("sources", []):
            if "compileGroupIndex" not in source:
                continue
            path = Path(source["path"])
            if not path.is_absolute():
                path = source_root / path
            require(relative_to(path, source_root) is not None,
                    "configured translation unit escapes source root: " + str(path))
            configured_units.add(path.resolve().relative_to(source_root).as_posix())
    require(configured_units == set(configured_translation_units),
            "configured translation-unit set mismatch: missing={} extra={}".format(
                sorted(set(configured_translation_units) - configured_units),
                sorted(configured_units - set(configured_translation_units))))

    targets = {}
    artifact_paths = {}
    relevant_units = set()
    for name in relevant_targets:
        document = documents[name]
        artifacts = document.get("artifacts", [])
        expected_artifact_count = (
            2 if schema == SCHEMA and
            name == "leopard2_backend_avx2" else 1)
        require(isinstance(artifacts, list) and
                len(artifacts) == expected_artifact_count,
                "evidence target artifact count changed: " + name)
        resolved_artifacts = []
        for artifact_record in artifacts:
            artifact = Path(artifact_record.get("path", ""))
            if not artifact.is_absolute():
                artifact = build_root / artifact
            artifact = artifact.resolve()
            require(relative_to(artifact, build_root) is not None,
                    "target artifact escapes fresh build: " + name)
            resolved_artifacts.append(artifact)
        artifact_paths[name] = resolved_artifacts[0]
        sources = []
        for source in document.get("sources", []):
            path = Path(source["path"])
            if not path.is_absolute():
                path = (source_root if not str(path).startswith("CMakeFiles/")
                        else build_root) / path
            tag = tagged_path(path, source_root, build_root)
            compiled = "compileGroupIndex" in source
            sources.append({"path": tag, "compiled": compiled})
            if compiled:
                require(tag.startswith("@source/"),
                        "relevant compiled source is not Git-backed: " + tag)
                relevant_units.add(tag[len("@source/"):])
        dependencies = [id_to_name.get(value.get("id"))
                        for value in document.get("dependencies", [])]
        require(all(isinstance(value, str) and value for value in dependencies),
                "unresolved CMake target dependency: " + name)
        dependencies.sort()
        link = document.get("link")
        normalized_link = None
        if link is not None:
            fragments = link.get("commandFragments")
            require(isinstance(fragments, list), "CMake link fragments: " + name)
            normalized_fragments = []
            for value in fragments:
                role = value.get("role")
                fragment = str(value.get("fragment", ""))
                if role == "libraries":
                    tokens = shlex.split(fragment)
                    require(len(tokens) == 1,
                            "CMake library fragment must be one path: " + name)
                    fragment = tagged_path(resolve_recipe_path(
                        tokens[0], build_root), source_root, build_root)
                else:
                    fragment = normalize_command(
                        [fragment], source_root, build_root, {})[0]
                normalized_fragments.append({"role": role,
                                             "fragment": fragment})
            normalized_link = {
                "language": link.get("language"),
                "fragments": normalized_fragments,
            }
        targets[name] = {
            "type": document.get("type"),
            "artifact": tagged_path(
                resolved_artifacts[0], source_root, build_root),
            "dependencies": dependencies,
            "sources": sorted(sources, key=lambda value: (
                value["path"], value["compiled"])),
            "link": normalized_link,
        }
        if schema == SCHEMA:
            targets[name]["artifacts"] = [
                tagged_path(value, source_root, build_root)
                for value in resolved_artifacts]
    require(relevant_units == set(build_translation_units),
            "evidence target translation-unit closure mismatch")
    require(targets["bench_leopard2"]["type"] == "EXECUTABLE" and
            targets["bench_leopard2"]["dependencies"] ==
            [cmake_identity["target"]],
            "benchmark target dependency identity")
    expected_library_dependencies = (
        ["leopard2_backend_avx2", "leopard2_backend_avx512",
         "leopard2_backend_ssse3"] if schema == SCHEMA else
        ["leopard2_backend_avx2", "leopard2_backend_ssse3"])
    expected_object_targets = (
        ("leopard2_backend_ssse3", "leopard2_backend_avx2",
         "leopard2_backend_avx512") if schema == SCHEMA else
        ("leopard2_backend_ssse3", "leopard2_backend_avx2"))
    require(targets[cmake_identity["target"]]["type"] == "STATIC_LIBRARY" and
            targets[cmake_identity["target"]]["dependencies"] ==
            expected_library_dependencies,
            "library target dependency identity")
    require(all(targets[name]["type"] == "OBJECT_LIBRARY" and
                targets[name]["dependencies"] == []
                for name in expected_object_targets),
            "backend object-target identity")
    record = {
        "configured_translation_units": sorted(configured_units),
        "targets": targets,
        "index_sha256": sha256_file(indexes[0]),
        "codemodel_sha256": sha256_file(codemodel_path),
    }
    record["digest"] = sha256_bytes(canonical_bytes(record))
    return artifact_paths, record


def fresh_rebuild(source_root, template_cache, fresh_root, jobs):
    source_root = Path(source_root).resolve()
    template_cache = Path(template_cache).resolve()
    build_root = Path(fresh_root).resolve()
    values = cache_values(template_cache)
    require(Path(values.get("CMAKE_HOME_DIRECTORY", "")).resolve() == source_root,
            "CMake cache/source-root mismatch")
    require(values.get("CMAKE_BUILD_TYPE") == "Release",
            "evidence build must be Release")
    require(values.get("LEO2_BUILD_BENCHMARKS", "ON") in ("ON", "1", "TRUE"),
            "evidence build omits benchmarks")
    require(values.get("CMAKE_EXPORT_COMPILE_COMMANDS") in ("ON", "1", "TRUE"),
            "evidence build omits compile commands")
    require(values.get("LEO2_BUILD_TESTS") in ("OFF", "0", "FALSE", ""),
            "evidence build must disable test hooks")
    require(values.get("LEO2_BUILD_FUZZERS") in ("OFF", "0", "FALSE", ""),
            "evidence build must disable fuzzers")
    require(isinstance(jobs, int) and 1 <= jobs <= 128,
            "build jobs must be in [1,128]")
    cmake_value = values.get("CMAKE_COMMAND")
    require(cmake_value, "CMake cache omits CMAKE_COMMAND")
    cmake = Path(cmake_value).resolve()
    require(cmake.is_file(), "cached CMake executable is missing")
    require(not build_root.exists(), "fresh build root already exists")
    query = build_root / ".cmake/api/v1/query"
    query.mkdir(parents=True)
    (query / "codemodel-v2").write_bytes(b"")
    generator = values.get("CMAKE_GENERATOR")
    require(generator, "template cache omits CMake generator")
    configure = [str(cmake), "-S", str(source_root), "-B", str(build_root),
                 "-G", generator]
    for key in CONFIGURATION_KEYS:
        if key == "CMAKE_GENERATOR":
            continue
        require(key in values, "template cache omits configuration key: " + key)
        configure.append("-D{}={}".format(key, values[key]))
    configure.extend([
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_C_COMPILER={}".format(values.get("CMAKE_C_COMPILER", "")),
        "-DCMAKE_CXX_COMPILER={}".format(values.get("CMAKE_CXX_COMPILER", "")),
        "-DCMAKE_AR={}".format(values.get("CMAKE_AR", "")),
        "-DCMAKE_RANLIB={}".format(values.get("CMAKE_RANLIB", "")),
        "-DCMAKE_MAKE_PROGRAM={}".format(values.get("CMAKE_MAKE_PROGRAM", "")),
    ])
    configure_tools = {"cmake": cmake}
    for logical_name, cache_key in TOOL_CACHE_KEYS:
        value = values.get(cache_key)
        require(value, "template cache omits tool: " + cache_key)
        configure_tools[logical_name] = Path(value)
    logical_configure = normalize_command(
        configure, source_root, build_root, configure_tools)
    logical_configure[0] = "@tool/cmake"
    configure_record = command_evidence(
        configure, source_root, "fresh CMake configure", logical_configure,
        FRESH_BUILD_ENVIRONMENT)
    artifacts, target_graph = file_api_targets(source_root, build_root, SCHEMA)
    cmake_identity = cmake_identity_for_schema(SCHEMA)
    argv = [str(cmake), "--build", str(build_root), "--parallel", str(jobs),
            "--target", cmake_identity["target"], "bench_leopard2"]
    logical_argv = ["@tool/cmake", "--build", "@build", "--parallel",
                    str(jobs), "--target", cmake_identity["target"],
                    "bench_leopard2"]
    build_record = command_evidence(
        argv, build_root, "fresh evidence build", logical_argv,
        FRESH_BUILD_ENVIRONMENT)
    version = command_output([str(cmake), "--version"], cmake.parent,
                             "cmake --version").decode("utf-8", "replace")
    rebuild = {
        "isolation": "new-empty-shadow-build",
        "environment": FRESH_BUILD_ENVIRONMENT,
        "configure": configure_record,
        "build": build_record,
        "cmake_binary_sha256": sha256_file(cmake),
        "cmake_version": version,
        "cmake_version_sha256": sha256_bytes(version.encode("utf-8")),
    }
    return {
        "build_root": build_root,
        "compile_commands": build_root / "compile_commands.json",
        "cmake_cache": build_root / "CMakeCache.txt",
        "library": artifacts[cmake_identity["target"]],
        "binary": artifacts["bench_leopard2"],
        "target_graph": target_graph,
        "rebuild": rebuild,
    }


def self_test_rebuild_record(jobs, cache_path, schema=SCHEMA):
    values = cache_values(cache_path)
    cmake_path = values["CMAKE_COMMAND"]
    configure_argv = ["@tool/cmake", "-S", "@source", "-B", "@build",
                      "-G", values["CMAKE_GENERATOR"]]
    for key in CONFIGURATION_KEYS:
        if key != "CMAKE_GENERATOR":
            configure_argv.append("-D{}={}".format(key, values[key]))
    configure_argv.extend([
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_C_COMPILER=@tool/cc",
        "-DCMAKE_CXX_COMPILER=@tool/cxx",
        "-DCMAKE_AR=@tool/ar",
        "-DCMAKE_RANLIB=@tool/ranlib",
        "-DCMAKE_MAKE_PROGRAM=@tool/make",
    ])
    cmake_identity = cmake_identity_for_schema(schema)
    build_argv = ["@tool/cmake", "--build", "@build", "--parallel",
                  str(jobs), "--target", cmake_identity["target"],
                  "bench_leopard2"]
    _, cmake = executable_identity(cmake_path, "cmake", "CMAKE_COMMAND")
    version = cmake["version"]
    empty_digest = sha256_bytes(b"")
    def record(argv):
        return {"argv": argv,
                "command_sha256": sha256_bytes(canonical_bytes(argv)),
                "stdout_sha256": empty_digest,
                "stderr_sha256": empty_digest,
                "returncode": 0}
    return {
        "isolation": "new-empty-shadow-build",
        "environment": FRESH_BUILD_ENVIRONMENT,
        "configure": record(configure_argv),
        "build": record(build_argv),
        "cmake_binary_sha256": cmake["binary_sha256"],
        "cmake_version": version,
        "cmake_version_sha256": sha256_bytes(version.encode("utf-8")),
    }


def normalized_configuration(values, source_root, build_root, tool_paths):
    configuration = {}
    for key in CONFIGURATION_KEYS:
        require(key in values, "CMake cache omits configuration key: " + key)
        configuration[key] = values[key]
    configuration["CMAKE_HOME_DIRECTORY"] = tagged_path(
        Path(values.get("CMAKE_HOME_DIRECTORY", "")), source_root, build_root)
    configuration["CMAKE_CACHEFILE_DIR"] = tagged_path(
        Path(values.get("CMAKE_CACHEFILE_DIR", "")), source_root, build_root)
    for logical_name, cache_key in TOOL_CACHE_KEYS:
        require(Path(values.get(cache_key, "")).resolve() ==
                tool_paths[logical_name].resolve(),
                "CMake tool cache identity mismatch: " + cache_key)
        configuration[cache_key] = "@tool/" + logical_name
    require(configuration["CMAKE_HOME_DIRECTORY"] == "@source/.",
            "normalized CMake source root")
    require(configuration["CMAKE_CACHEFILE_DIR"] == "@build/.",
            "normalized CMake build root")
    require(configuration["CMAKE_BUILD_TYPE"] == "Release",
            "evidence build must be Release")
    require(configuration["LEO2_BACKEND_VARIANT"] == "auto",
            "evidence build must use runtime backend dispatch")
    require(configuration["LEO2_BUILD_BENCHMARKS"] in ("ON", "1", "TRUE"),
            "evidence build omits benchmarks")
    require(configuration["LEO2_BUILD_TESTS"] in ("OFF", "0", "FALSE", "") and
            configuration["LEO2_BUILD_FUZZERS"] in ("OFF", "0", "FALSE", ""),
            "evidence build must not contain tests/fuzzers or test hooks")
    require(configuration["LEO2_ENABLE_CUDA"] in ("OFF", "0", "FALSE", ""),
            "butterfly evidence must not include optional CUDA")
    return configuration


def resolve_recipe_path(value, directory):
    path = Path(value)
    if not path.is_absolute():
        path = directory / path
    return path.resolve()


def literal_link_recipes(source_root, build_root, tool_paths, compile_entries,
                         library, binary, schema=SCHEMA):
    cmake_identity = cmake_identity_for_schema(schema)
    library_recipe_path = build_root / "CMakeFiles" / \
        cmake_identity["target_directory"] / "link.txt"
    benchmark_recipe_path = build_root / "CMakeFiles/bench_leopard2.dir/link.txt"
    require(library_recipe_path.is_file() and benchmark_recipe_path.is_file(),
            "evidence collection requires CMake literal link.txt recipes")
    library_content = exact_utf8_file_content(
        library_recipe_path, "static-library link recipe") \
        if schema in HARDENED_SCHEMAS else None
    benchmark_content = exact_utf8_file_content(
        benchmark_recipe_path, "benchmark link recipe") \
        if schema in HARDENED_SCHEMAS else None
    library_text = (library_content["text"] if library_content is not None else
                    library_recipe_path.read_text(encoding="utf-8"))
    benchmark_text = (benchmark_content["text"]
                      if benchmark_content is not None else
                      benchmark_recipe_path.read_text(encoding="utf-8"))
    library_lines = [shlex.split(value) for value in library_text.splitlines()
                     if value.strip()]
    benchmark_lines = [shlex.split(value) for value in
                       benchmark_text.splitlines()
                       if value.strip()]
    require(len(library_lines) == 2 and len(benchmark_lines) == 1,
            "unexpected CMake archive/link recipe shape")
    archive, index = library_lines
    require(command_executable(archive[0], build_root) == tool_paths["ar"] and
            len(archive) >= 4 and archive[1] in ("qc", "rc", "rcs"),
            "static-library recipe does not use captured CMAKE_AR")
    require(resolve_recipe_path(archive[2], build_root) == library,
            "static-library recipe output is not target artifact")
    expected_objects = {
        value[1].resolve() for relative, value in compile_entries.items()
        if relative != "bench/leopard2/benchmark.cpp"
    }
    recipe_objects = [resolve_recipe_path(value, build_root)
                      for value in archive[3:]]
    require(len(set(recipe_objects)) == len(recipe_objects) and
            set(recipe_objects) == expected_objects,
            "static-library recipe object set mismatch")
    require(len(index) == 2 and
            command_executable(index[0], build_root) == tool_paths["ranlib"] and
            resolve_recipe_path(index[1], build_root) == library,
            "static-library index recipe does not use captured CMAKE_RANLIB")

    link = benchmark_lines[0]
    require(command_executable(link[0], build_root) == tool_paths["cxx"],
            "benchmark link recipe does not use captured C++ driver")
    require(all(not value.startswith("@") and ",@" not in value for value in link),
            "response-file link recipes are not accepted")
    require(command_output_path(link, build_root) == binary,
            "benchmark recipe output is not target artifact")
    benchmark_object = compile_entries["bench/leopard2/benchmark.cpp"][1].resolve()
    resolved_file_inputs = []
    external_inputs = []
    output_indices = set()
    for index_value, value in enumerate(link):
        if value == "-o" and index_value + 1 < len(link):
            output_indices.add(index_value + 1)
        elif value.startswith("-o") and len(value) > 2:
            output_indices.add(index_value)
    for index_value, value in enumerate(link[1:], 1):
        if index_value in output_indices or value.startswith("-"):
            continue
        path = resolve_recipe_path(value, build_root)
        if not path.is_file():
            continue
        resolved_file_inputs.append(path)
        if relative_to(path, build_root) is None:
            require(relative_to(path, source_root) is None,
                    "benchmark links a source-tree file directly")
            external = {"path": tagged_path(path, source_root, build_root),
                        "sha256": sha256_file(path)}
            if schema in HARDENED_SCHEMAS:
                external["raw_path"] = value
            external_inputs.append(external)
    require(resolved_file_inputs.count(benchmark_object) == 1 and
            resolved_file_inputs.count(library) == 1,
            "benchmark link inputs omit or duplicate target object/archive")
    build_inputs = {path for path in resolved_file_inputs
                    if relative_to(path, build_root) is not None}
    require(build_inputs == {benchmark_object, library},
            "benchmark link recipe contains an unexpected build input")

    def normalized_recipe(line, logical_tool):
        result = list(line)
        result[0] = "@tool/" + logical_tool
        for index_value, value in enumerate(result[1:], 1):
            if value.startswith("-"):
                continue
            candidate = resolve_recipe_path(value, build_root)
            if candidate.exists() or (index_value > 0 and
                                      result[index_value - 1] == "-o"):
                result[index_value] = tagged_path(
                    candidate, source_root, build_root)
        return result

    recipes = {
        "library": [normalized_recipe(library_lines[0], "ar"),
                    normalized_recipe(library_lines[1], "ranlib")],
        "benchmark": [normalized_recipe(link, "cxx")],
        "external_link_inputs": sorted(
            external_inputs, key=lambda value: value["path"]),
        "library_recipe_sha256": (
            library_content["sha256"] if library_content is not None else
            sha256_file(library_recipe_path)),
        "benchmark_recipe_sha256": (
            benchmark_content["sha256"] if benchmark_content is not None else
            sha256_file(benchmark_recipe_path)),
    }
    if schema in HARDENED_SCHEMAS:
        recipes["library_recipe_size"] = library_content["size"]
        recipes["benchmark_recipe_size"] = benchmark_content["size"]
        recipes["library_recipe_content"] = library_content
        recipes["benchmark_recipe_content"] = benchmark_content
    recipes["digest"] = sha256_bytes(canonical_bytes(recipes))
    return recipes, recipe_objects


def archive_manifest(library, ar_tool, expected_objects):
    listing = command_output([str(ar_tool), "t", str(library)], library.parent,
                             "archive member listing").decode("utf-8").splitlines()
    require(listing and len(listing) == len(set(listing)),
            "archive is empty or has duplicate member names")
    expected_by_name = {}
    for path in expected_objects:
        require(path.name not in expected_by_name,
                "expected archive objects have duplicate basenames")
        expected_by_name[path.name] = path
    require(listing == [path.name for path in expected_objects],
            "archive member order differs from target object recipe")
    members = []
    for name in listing:
        raw = command_output([str(ar_tool), "p", str(library), name],
                             library.parent, "archive member extraction")
        object_path = expected_by_name[name]
        require(sha256_bytes(raw) == sha256_file(object_path),
                "archive member differs from freshly built object: " + name)
        members.append({"name": name,
                        "object": "@build/" +
                        object_path.relative_to(library.parent).as_posix(),
                        "sha256": sha256_bytes(raw)})
    result = {"members": members, "member_count": len(members)}
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def validate_retained_link_recipe_semantics(recipes, tools, units, schema):
    library_lines = parse_exact_recipe_content(
        recipes["library_recipe_content"],
        recipes["library_recipe_size"],
        recipes["library_recipe_sha256"], "static-library link recipe")
    benchmark_lines = parse_exact_recipe_content(
        recipes["benchmark_recipe_content"],
        recipes["benchmark_recipe_size"],
        recipes["benchmark_recipe_sha256"], "benchmark link recipe")
    require(len(library_lines) == 2 and len(benchmark_lines) == 1,
            "retained CMake archive/link recipe shape")
    archive, index = library_lines
    benchmark = benchmark_lines[0]
    require(archive[0] == tools["ar"]["recipe_argv0"] and
            index[0] == tools["ranlib"]["recipe_argv0"] and
            benchmark[0] == tools["cxx"]["recipe_argv0"],
            "retained recipe tool spelling differs from recorded tool")

    normalized_archive = recipes["library"][0]
    require(len(archive) == len(normalized_archive) and len(archive) >= 4 and
            archive[1] == normalized_archive[1] and
            archive[1] in ("qc", "rc", "rcs") and
            archive[2] == "libleopard.a" and
            normalized_archive[2] == "@build/libleopard.a",
            "retained archive command/output semantics")
    require(all(value.startswith("@build/")
                for value in normalized_archive[3:]),
            "normalized archive object location")
    expected_raw_objects = [value[len("@build/"):]
                            for value in normalized_archive[3:]]
    require(archive[3:] == expected_raw_objects and
            len(expected_raw_objects) == len(set(expected_raw_objects)),
            "retained archive object closure/order")
    for relative, entry in units.items():
        if relative == "bench/leopard2/benchmark.cpp":
            expected_prefix = "@build/CMakeFiles/bench_leopard2.dir/"
        elif relative == "Leopard2BackendSSSE3.cpp":
            expected_prefix = "@build/CMakeFiles/leopard2_backend_ssse3.dir/"
        elif relative in (
                "Leopard2BackendAVX2.cpp", "Leopard2BackendAVX2Xor.cpp"):
            expected_prefix = "@build/CMakeFiles/leopard2_backend_avx2.dir/"
        elif relative == "Leopard2BackendAVX512.cpp":
            expected_prefix = \
                "@build/CMakeFiles/leopard2_backend_avx512.dir/"
        else:
            expected_prefix = "@build/CMakeFiles/leopard.dir/"
        require(entry["output"].startswith(expected_prefix),
                "canonical CMake object directory: " + relative)
    expected_library_outputs = [
        units[relative]["output"]
        for relative in archive_translation_units_for_schema(schema)]
    require(normalized_archive[3:] == expected_library_outputs,
            "retained archive translation-unit closure/order")
    require(index == [tools["ranlib"]["recipe_argv0"], "libleopard.a"] and
            recipes["library"][1] ==
            ["@tool/ranlib", "@build/libleopard.a"],
            "retained ranlib archive identity")

    normalized_benchmark = recipes["benchmark"][0]
    require(len(benchmark) == len(normalized_benchmark) and
            all(not value.startswith("@") and ",@" not in value
                for value in benchmark),
            "retained benchmark response-file/argv shape")
    external_by_path = {
        value["path"]: value for value in recipes["external_link_inputs"]}
    require(len(external_by_path) == len(recipes["external_link_inputs"]),
            "retained benchmark external-input uniqueness")
    for offset, normalized in enumerate(normalized_benchmark):
        if offset == 0:
            require(normalized == "@tool/cxx",
                    "normalized benchmark tool identity")
        elif normalized.startswith("@build/"):
            require(benchmark[offset] == normalized[len("@build/"):],
                    "retained benchmark build-path identity")
        elif normalized.startswith("@external/"):
            external = external_by_path.get(normalized)
            raw_external = (Path(external["raw_path"])
                            if external is not None else None)
            require(external is not None and
                    benchmark[offset] == external["raw_path"] and
                    raw_external.is_absolute() and raw_external.is_file() and
                    raw_external.resolve().name == Path(normalized).name and
                    sha256_file(raw_external) == external["sha256"],
                    "retained benchmark external-path identity")
        else:
            require(benchmark[offset] == normalized,
                    "retained benchmark flag/argument identity")
    output_values = []
    for offset, value in enumerate(benchmark):
        if value == "-o":
            require(offset + 1 < len(benchmark),
                    "retained benchmark bare -o")
            output_values.append(benchmark[offset + 1])
        elif value.startswith("-o") and len(value) > 2:
            output_values.append(value[2:])
    benchmark_object = units["bench/leopard2/benchmark.cpp"]["output"]
    require(output_values == ["bench_leopard2"] and
            benchmark.count("libleopard.a") == 1 and
            benchmark.count(benchmark_object[len("@build/"):]) == 1 and
            normalized_benchmark.count("@build/libleopard.a") == 1,
            "retained benchmark output/object/archive semantics")


def build_record(source_root, source_identity, compile_commands, cmake_cache,
                 library, binary, rebuild, target_graph, schema=SCHEMA):
    configured_translation_units = configured_translation_units_for_schema(
        schema)
    build_translation_units = build_translation_units_for_schema(schema)
    source_root = Path(source_root).resolve()
    compile_commands = Path(compile_commands).resolve()
    cmake_cache = Path(cmake_cache).resolve()
    library = Path(library).resolve()
    binary = Path(binary).resolve()
    build_root = cmake_cache.parent
    require(compile_commands.parent == build_root,
            "compile commands and CMake cache have different build roots")
    require(relative_to(library, build_root) is not None and
            relative_to(binary, build_root) is not None,
            "build artifacts escape declared build root")
    for path, label in ((compile_commands, "compile_commands.json"),
                        (cmake_cache, "CMakeCache.txt"),
                        (library, "static library"), (binary, "benchmark binary")):
        require(path.is_file(), label + " missing")
    document = read_json(compile_commands, "compile_commands.json")
    require(isinstance(document, list), "compile_commands.json is not an array")
    values = cache_values(cmake_cache)
    require(Path(values.get("CMAKE_HOME_DIRECTORY", "")).resolve() == source_root,
            "CMake cache/source-root mismatch")
    tool_paths, tools = cache_tools(cmake_cache)
    all_by_file = {}
    for entry in document:
        relative, normalized = normalized_compile_entry(
            entry, source_root, build_root, tool_paths)
        require(relative not in all_by_file,
                "duplicate compile command for " + relative)
        all_by_file[relative] = normalized
    if schema == SCHEMA:
        # compile_commands.json is global to the configured build tree, while
        # the current evidence contract deliberately scopes the executable
        # closure to bench_leopard2 and its target dependencies.  Require every
        # relevant command exactly once, retain the digest of the complete
        # compile database below, and record only the commands that can reach
        # the benchmark.  Unrelated configured tools must not become an
        # accidental part of the benchmark's source identity.
        missing_commands = set(configured_translation_units) - set(all_by_file)
        require(not missing_commands,
                "configured compile-command set omits relevant units: {}".format(
                    sorted(missing_commands)))
        configured_by_file = {
            relative: all_by_file[relative]
            for relative in configured_translation_units
        }
    else:
        # Historical schemas described the entire configured compile database;
        # preserve that exact replay contract.
        require(set(all_by_file) == set(configured_translation_units),
                "configured compile-command set mismatch: missing={} extra={}".format(
                    sorted(set(configured_translation_units) - set(all_by_file)),
                    sorted(set(all_by_file) - set(configured_translation_units))))
        configured_by_file = all_by_file
    by_file = {relative: all_by_file[relative]
               for relative in build_translation_units}
    require(set(by_file) == set(build_translation_units),
            "compile command translation-unit set mismatch: missing={} extra={}".format(
                sorted(set(build_translation_units) - set(by_file)),
                sorted(set(by_file) - set(build_translation_units))))
    require(all(value[0]["output_sha256"] is not None for value in by_file.values()),
            "evidence target compile output is missing")
    commands = {relative: by_file[relative][0]
                for relative in sorted(by_file)}
    recipes, expected_objects = literal_link_recipes(
        source_root, build_root, tool_paths, by_file, library, binary, schema)
    if schema in HARDENED_SCHEMAS:
        retained_library = parse_exact_recipe_content(
            recipes["library_recipe_content"],
            recipes["library_recipe_size"],
            recipes["library_recipe_sha256"], "static-library link recipe")
        retained_benchmark = parse_exact_recipe_content(
            recipes["benchmark_recipe_content"],
            recipes["benchmark_recipe_size"],
            recipes["benchmark_recipe_sha256"], "benchmark link recipe")
        require(len(retained_library) == 2 and len(retained_benchmark) == 1,
                "retained link recipe shape")
        tools["ar"]["recipe_argv0"] = retained_library[0][0]
        tools["ranlib"]["recipe_argv0"] = retained_library[1][0]
        tools["cxx"]["recipe_argv0"] = retained_benchmark[0][0]
    archive = archive_manifest(library, tool_paths["ar"], expected_objects)
    record = {
        "source": source_identity,
        "configuration": normalized_configuration(
            values, source_root, build_root, tool_paths),
        "build_input_files": {
            "compile_commands.json": sha256_file(compile_commands),
            "CMakeCache.txt": sha256_file(cmake_cache),
        },
        "tools": tools,
        "translation_units": commands,
        "configured_translation_units": sorted(configured_by_file),
        "dependency_closure": dependency_closure(
            by_file, source_root, build_root),
        "target_graph": target_graph,
        "link_recipes": recipes,
        "archive": archive,
        "rebuild": rebuild,
        "artifacts": {
            "library_sha256": sha256_file(library),
            "benchmark_sha256": sha256_file(binary),
        },
    }
    record["digest"] = sha256_bytes(canonical_bytes(record))
    return record


def validate_build_record(record, repo, schema):
    cmake_identity = cmake_identity_for_schema(schema)
    relevant_targets = relevant_targets_for_schema(schema)
    configured_translation_units = configured_translation_units_for_schema(
        schema)
    build_translation_units = build_translation_units_for_schema(schema)
    library_target = cmake_identity["target"]
    library_artifact = "@build/" + cmake_identity["archive"]
    expected_keys = {
        "source", "configuration", "build_input_files", "tools",
        "translation_units", "configured_translation_units",
        "dependency_closure", "target_graph", "link_recipes", "archive",
        "rebuild", "artifacts", "digest",
    }
    require(set(record) == expected_keys, "build record keys")
    payload = dict(record)
    digest = payload.pop("digest")
    require(re.fullmatch(r"[0-9a-f]{64}", digest or "") is not None,
            "build-record digest format")
    require(sha256_bytes(canonical_bytes(payload)) == digest,
            "build-record digest mismatch")
    validate_git_record(repo, record["source"])
    configuration = record["configuration"]
    require(set(configuration) == set(CONFIGURATION_KEYS) | {
        "CMAKE_HOME_DIRECTORY", "CMAKE_CACHEFILE_DIR"} |
        {cache_key for _, cache_key in TOOL_CACHE_KEYS},
        "configuration key set")
    require(configuration["CMAKE_HOME_DIRECTORY"] == "@source/." and
            configuration["CMAKE_CACHEFILE_DIR"] == "@build/." and
            configuration["CMAKE_BUILD_TYPE"] == "Release" and
            configuration["LEO2_BACKEND_VARIANT"] == "auto" and
            configuration["LEO2_BUILD_BENCHMARKS"] in ("ON", "1", "TRUE") and
            configuration["LEO2_BUILD_TESTS"] in ("OFF", "0", "FALSE", "") and
            configuration["LEO2_BUILD_FUZZERS"] in ("OFF", "0", "FALSE", "") and
            configuration["LEO2_ENABLE_CUDA"] in ("OFF", "0", "FALSE", ""),
            "configuration identity")
    require(set(record["build_input_files"]) == {
        "compile_commands.json", "CMakeCache.txt"},
        "build input key set")
    require(set(record["tools"]) ==
            {logical_name for logical_name, _ in TOOL_CACHE_KEYS},
            "tool identity set")
    for logical_name, cache_key in TOOL_CACHE_KEYS:
        tool = record["tools"][logical_name]
        expected_tool_keys = {"logical_name", "cache_key", "basename",
                              "binary_sha256", "version", "version_sha256"}
        if schema in HARDENED_SCHEMAS and logical_name in (
                "ar", "ranlib", "cxx"):
            expected_tool_keys.add("recipe_argv0")
        require(set(tool) == expected_tool_keys,
                "tool identity keys: " + logical_name)
        require(tool["logical_name"] == logical_name and
                tool["cache_key"] == cache_key and
                isinstance(tool["basename"], str) and tool["basename"],
                "tool logical identity: " + logical_name)
        if "recipe_argv0" in expected_tool_keys:
            require(isinstance(tool["recipe_argv0"], str) and
                    tool["recipe_argv0"] and "\x00" not in tool["recipe_argv0"],
                    "recipe tool spelling: " + logical_name)
        require(configuration[cache_key] == "@tool/" + logical_name,
                "configuration/tool binding: " + logical_name)
        require(sha256_bytes(tool["version"].encode("utf-8")) ==
                tool["version_sha256"], "tool version digest: " + logical_name)
        require(re.fullmatch(r"[0-9a-f]{64}", tool["binary_sha256"] or "") is not None,
                "tool binary digest: " + logical_name)
    units = record["translation_units"]
    require(set(units) == set(build_translation_units),
            "build translation-unit set")
    for relative, entry in units.items():
        require(set(entry) == {
            "file", "argv", "directory", "output", "output_sha256",
            "command_sha256"},
            "translation-unit keys: " + relative)
        require(entry["file"] == relative and entry["directory"].startswith("@build/") and
                entry["output"].startswith("@build/") and
                isinstance(entry["argv"], list) and entry["argv"] and
                entry["argv"][0] == "@tool/cxx",
                "translation-unit identity: " + relative)
        require(re.fullmatch(r"[0-9a-f]{64}", entry["output_sha256"] or "")
                is not None, "translation-unit output digest: " + relative)
        require(all(isinstance(value, str) for value in entry["argv"]) and
                entry["argv"].count("@source/" + relative) == 1,
                "translation-unit source argument: " + relative)
        payload = dict(entry)
        command_digest = payload.pop("command_sha256")
        require(sha256_bytes(canonical_bytes(payload)) == command_digest,
                "translation-unit command digest: " + relative)
    require(record["configured_translation_units"] ==
            sorted(configured_translation_units),
            "configured translation-unit closure")
    require(set(record["artifacts"]) == {
        "library_sha256", "benchmark_sha256"}, "build artifact key set")
    closure = record["dependency_closure"]
    require(set(closure) == {
        "file_count", "source_file_count", "manifest_sha256", "manifest",
        "by_translation_unit"},
        "dependency-closure key set")
    manifest = closure["manifest"]
    require(isinstance(manifest, list) and manifest == sorted(
        manifest, key=lambda value: value.get("path", "")),
        "dependency manifest ordering")
    require(len(manifest) == closure["file_count"] and
            sha256_bytes(canonical_bytes(manifest)) == closure["manifest_sha256"],
            "dependency manifest digest/count")
    manifest_paths = {}
    source_count = 0
    for dependency in manifest:
        require(set(dependency) == {"path", "sha256"},
                "dependency manifest entry keys")
        path = dependency["path"]
        digest_value = dependency["sha256"]
        require(path not in manifest_paths and
                re.fullmatch(r"[0-9a-f]{64}", digest_value or "") is not None,
                "dependency manifest identity")
        manifest_paths[path] = digest_value
        if path.startswith("@source/"):
            source_count += 1
            relative = path[len("@source/"):]
            require(relative and not relative.startswith("/") and ".." not in Path(relative).parts,
                    "unsafe source dependency path")
            require(sha256_bytes(git_blob(
                repo, record["source"]["commit"], relative)) == digest_value,
                "source dependency differs from declared Git tree: " + relative)
        elif path.startswith("@build/"):
            require(not path.startswith("@build/../"),
                    "unsafe build dependency path")
        else:
            require(re.fullmatch(r"@external/[0-9a-f]{64}/[^/]+", path) is not None,
                    "unclassified external dependency path")
    require(source_count == closure["source_file_count"],
            "dependency source-file count")
    by_unit = closure["by_translation_unit"]
    require(set(by_unit) == set(build_translation_units),
            "dependency translation-unit set")
    for relative, dependencies in by_unit.items():
        require(isinstance(dependencies, list) and dependencies == sorted(set(dependencies)) and
                "@source/" + relative in dependencies and
                all(path in manifest_paths for path in dependencies),
                "dependency list identity: " + relative)
    require(isinstance(closure["file_count"], int) and
            closure["file_count"] >= len(build_translation_units),
            "dependency-closure file count")
    require(isinstance(closure["source_file_count"], int) and
            len(build_translation_units) <= closure["source_file_count"] <=
            closure["file_count"], "dependency-closure source count")

    target_graph = record["target_graph"]
    require(set(target_graph) == {"configured_translation_units", "targets",
                                 "index_sha256", "codemodel_sha256", "digest"},
            "target graph key set")
    graph_payload = dict(target_graph)
    graph_digest = graph_payload.pop("digest")
    require(sha256_bytes(canonical_bytes(graph_payload)) == graph_digest and
            target_graph["configured_translation_units"] ==
            sorted(configured_translation_units) and
            set(target_graph["targets"]) == set(relevant_targets),
            "target graph identity")
    targets = target_graph["targets"]
    expected_library_dependencies = (
        ["leopard2_backend_avx2", "leopard2_backend_avx512",
         "leopard2_backend_ssse3"] if schema == SCHEMA else
        ["leopard2_backend_avx2", "leopard2_backend_ssse3"])
    require(targets["bench_leopard2"]["type"] == "EXECUTABLE" and
            targets["bench_leopard2"]["artifact"] == "@build/bench_leopard2" and
            targets["bench_leopard2"]["dependencies"] == [library_target] and
            targets[library_target]["type"] == "STATIC_LIBRARY" and
            targets[library_target]["artifact"] == library_artifact and
            targets[library_target]["dependencies"] ==
            expected_library_dependencies,
            "target dependency graph")
    compiled_target_sources = set()
    for name, target in targets.items():
        expected_target_keys = {
            "type", "artifact", "dependencies", "sources", "link"}
        if schema == SCHEMA:
            expected_target_keys.add("artifacts")
        require(set(target) == expected_target_keys and
                target["artifact"].startswith("@build/") and
                isinstance(target["sources"], list),
                "target record: " + name)
        if schema == SCHEMA:
            expected_count = 2 if name == "leopard2_backend_avx2" else 1
            require(isinstance(target["artifacts"], list) and
                    len(target["artifacts"]) == expected_count and
                    target["artifacts"][0] == target["artifact"] and
                    len(set(target["artifacts"])) == expected_count and
                    all(value.startswith("@build/")
                        for value in target["artifacts"]),
                    "target artifact closure: " + name)
            if name == "leopard2_backend_avx2":
                require([Path(value).name for value in target["artifacts"]] == [
                            "Leopard2BackendAVX2.cpp.o",
                            "Leopard2BackendAVX2Xor.cpp.o"],
                        "AVX2 object artifact order changed")
        for source in target["sources"]:
            require(set(source) == {"path", "compiled"} and
                    isinstance(source["compiled"], bool),
                    "target source record: " + name)
            if source["compiled"]:
                require(source["path"].startswith("@source/"),
                        "compiled target source is not Git-backed")
                compiled_target_sources.add(source["path"][len("@source/"):])
    require(compiled_target_sources == set(build_translation_units),
            "target compiled-source closure")

    recipes = record["link_recipes"]
    expected_recipe_keys = {
        "library", "benchmark", "external_link_inputs",
        "library_recipe_sha256", "benchmark_recipe_sha256", "digest"}
    if schema in HARDENED_SCHEMAS:
        expected_recipe_keys.update({
            "library_recipe_content", "benchmark_recipe_content",
            "library_recipe_size", "benchmark_recipe_size"})
    require(set(recipes) == expected_recipe_keys, "link recipe key set")
    recipe_payload = dict(recipes)
    recipe_digest = recipe_payload.pop("digest")
    require(sha256_bytes(canonical_bytes(recipe_payload)) == recipe_digest and
            len(recipes["library"]) == 2 and len(recipes["benchmark"]) == 1 and
            recipes["library"][0][0] == "@tool/ar" and
            recipes["library"][1][0] == "@tool/ranlib" and
            recipes["benchmark"][0][0] == "@tool/cxx",
            "literal archive/link recipe identity")
    library_command = recipes["library"][0]
    require(len(library_command) >= 4 and library_command[1] in ("qc", "rc", "rcs") and
            library_command[2] == library_artifact,
            "literal archive output recipe")
    expected_library_objects = {
        entry["output"] for relative, entry in units.items()
        if relative != "bench/leopard2/benchmark.cpp"}
    require(len(library_command[3:]) == len(set(library_command[3:])) and
            set(library_command[3:]) == expected_library_objects and
            recipes["library"][1] ==
            ["@tool/ranlib", library_artifact],
            "literal archive object/index closure")
    benchmark_command = recipes["benchmark"][0]
    expected_benchmark_object = units[
        "bench/leopard2/benchmark.cpp"]["output"]
    for external in recipes["external_link_inputs"]:
        expected_external_keys = {"path", "sha256"}
        if schema in HARDENED_SCHEMAS:
            expected_external_keys.add("raw_path")
        require(set(external) == expected_external_keys and
                external["path"].startswith("@external/") and
                re.fullmatch(r"[0-9a-f]{64}", external["sha256"]) is not None and
                (schema not in HARDENED_SCHEMAS or
                 (isinstance(external["raw_path"], str) and
                  external["raw_path"] and
                  not external["raw_path"].startswith("@"))),
                "external link input identity")
        if schema in HARDENED_SCHEMAS:
            raw_external = Path(external["raw_path"])
            require(raw_external.is_absolute() and raw_external.is_file() and
                    raw_external.resolve().name == Path(external["path"]).name and
                    sha256_file(raw_external) == external["sha256"],
                    "external raw/resolved link identity")
    file_api_libraries = [
        value["fragment"]
        for value in targets["bench_leopard2"]["link"]["fragments"]
        if value["role"] == "libraries"]
    require(file_api_libraries.count(library_artifact) == 1,
            "CMake File API benchmark library target binding")
    expected_external_links = sorted(
        value for value in file_api_libraries if value.startswith("@external/"))
    require(expected_external_links == sorted(
        value["path"] for value in recipes["external_link_inputs"]),
        "literal link external inputs differ from CMake target metadata")
    link = targets["bench_leopard2"]["link"]
    require(isinstance(link, dict) and set(link) == {"language", "fragments"} and
            link["language"] == "CXX" and isinstance(link["fragments"], list),
            "CMake File API benchmark link identity")
    expected_benchmark_command = ["@tool/cxx"]
    for fragment in link["fragments"]:
        require(set(fragment) == {"role", "fragment"} and
                fragment["role"] in ("flags", "libraries") and
                isinstance(fragment["fragment"], str),
                "CMake File API benchmark link fragment")
        if fragment["role"] == "flags":
            expected_benchmark_command.extend(shlex.split(fragment["fragment"]))
    expected_benchmark_command.extend([
        expected_benchmark_object, "-o", "@build/bench_leopard2"])
    expected_benchmark_command.extend(file_api_libraries)
    require(benchmark_command == expected_benchmark_command,
            "literal benchmark recipe differs from CMake target metadata")
    if schema in HARDENED_SCHEMAS:
        validate_retained_link_recipe_semantics(
            recipes, record["tools"], units, schema)

    archive = record["archive"]
    require(set(archive) == {"members", "member_count", "digest"},
            "archive manifest key set")
    archive_payload = dict(archive)
    archive_digest = archive_payload.pop("digest")
    require(sha256_bytes(canonical_bytes(archive_payload)) == archive_digest and
            archive["member_count"] == len(archive["members"]) ==
            len(build_translation_units) - 1,
            "archive manifest identity")
    expected_archive_members = {
        entry["output"]: entry["output_sha256"]
        for relative, entry in units.items()
        if relative != "bench/leopard2/benchmark.cpp"
    }
    member_names = []
    member_objects = []
    for member in archive["members"]:
        require(set(member) == {"name", "object", "sha256"} and
                member["object"] in expected_archive_members and
                member["name"] == Path(member["object"]).name and
                member["sha256"] == expected_archive_members[member["object"]],
                "archive member identity")
        member_names.append(member["name"])
        member_objects.append(member["object"])
    require(len(member_names) == len(set(member_names)) and
            set(member_objects) == set(expected_archive_members) and
            member_objects == library_command[3:],
            "archive member/object closure")
    rebuild = record["rebuild"]
    require(set(rebuild) == {"isolation", "environment", "configure", "build",
                            "cmake_binary_sha256", "cmake_version",
                            "cmake_version_sha256"} and
            rebuild["isolation"] == "new-empty-shadow-build",
            "rebuild key/isolation identity")
    require(rebuild["environment"] == FRESH_BUILD_ENVIRONMENT,
            "fresh build environment identity")
    for phase in ("configure", "build"):
        command = rebuild[phase]
        require(set(command) == {"argv", "command_sha256", "stdout_sha256",
                                "stderr_sha256", "returncode"} and
                command["returncode"] == 0 and
                sha256_bytes(canonical_bytes(command["argv"])) ==
                command["command_sha256"],
                "rebuild {} command identity".format(phase))
    expected_configure = [
        "@tool/cmake", "-S", "@source", "-B", "@build", "-G",
        configuration["CMAKE_GENERATOR"],
    ]
    for key in CONFIGURATION_KEYS:
        if key != "CMAKE_GENERATOR":
            expected_configure.append("-D{}={}".format(key, configuration[key]))
    expected_configure.extend([
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_C_COMPILER=@tool/cc",
        "-DCMAKE_CXX_COMPILER=@tool/cxx",
        "-DCMAKE_AR=@tool/ar",
        "-DCMAKE_RANLIB=@tool/ranlib",
        "-DCMAKE_MAKE_PROGRAM=@tool/make",
    ])
    require(rebuild["configure"]["argv"] == expected_configure,
            "fresh configure argv")
    build_argv = rebuild["build"]["argv"]
    require(build_argv[:4] == ["@tool/cmake", "--build", "@build", "--parallel"] and
            build_argv[5:] == ["--target", library_target, "bench_leopard2"] and
            build_argv[4].isdigit() and 1 <= int(build_argv[4]) <= 128,
            "fresh build argv")
    require(sha256_bytes(rebuild["cmake_version"].encode("utf-8")) ==
            rebuild["cmake_version_sha256"], "CMake version digest")
    require(rebuild["cmake_binary_sha256"] ==
            record["tools"]["cmake"]["binary_sha256"] and
            rebuild["cmake_version"] == record["tools"]["cmake"]["version"],
            "rebuild/CMake tool binding")
    for value in ([entry["command_sha256"] for entry in units.values()] +
                  [entry["output_sha256"] for entry in units.values()] + [
            record["build_input_files"].get("compile_commands.json"),
            record["build_input_files"].get("CMakeCache.txt"),
            record["artifacts"].get("library_sha256"),
            record["artifacts"].get("benchmark_sha256"),
            closure.get("manifest_sha256"), target_graph.get("digest"),
            recipes.get("digest"), archive.get("digest"),
            rebuild["configure"].get("command_sha256"),
            rebuild["build"].get("command_sha256"),
            rebuild.get("cmake_binary_sha256"),
            rebuild.get("cmake_version_sha256"),
            rebuild["configure"].get("stdout_sha256"),
            rebuild["configure"].get("stderr_sha256"),
            rebuild["build"].get("stdout_sha256"),
            rebuild["build"].get("stderr_sha256")]):
        require(re.fullmatch(r"[0-9a-f]{64}", value or "") is not None,
                "build record SHA-256 format")


def validate_matrix_compiler(compiler, label):
    require(set(compiler) == {"executable", "binary_sha256", "version",
                             "version_sha256"} and
            isinstance(compiler["executable"], str) and
            Path(compiler["executable"]).is_absolute() and
            re.fullmatch(r"[0-9a-f]{64}", compiler["binary_sha256"]) is not None and
            isinstance(compiler["version"], str) and compiler["version"] and
            sha256_bytes(compiler["version"].encode("utf-8")) ==
            compiler["version_sha256"], label + " compiler identity")


def compiler_version_signature(version):
    """Remove only an argv[0]-dependent GCC driver label.

    GCC prints the invoked driver name before the parenthesized package
    identity.  CMake records the resolved executable while a matrix runner may
    have invoked the same binary through cc/c++.  The executable basename and
    binary digest are checked separately, so retaining the package/version and
    all following lines is the stable semantic comparison.
    """
    lines = version.strip().splitlines()
    _, marker, suffix = lines[0].partition(" (")
    if marker:
        lines[0] = "(" + suffix
    return "\n".join(lines)


def validate_matrix_command(command, label, extra_keys=()):
    keys = {"argv", "cwd", "label", "returncode", "stderr_log",
            "stdout_log", "timed_out", "stderr_sha256", "stdout_sha256",
            "command_sha256"} | set(extra_keys)
    require(set(command) == keys and command["label"] == label and
            command["returncode"] == 0 and command["timed_out"] is False and
            command["stdout_log"] == label + ".stdout.log" and
            command["stderr_log"] == label + ".stderr.log" and
            isinstance(command["argv"], list) and command["argv"] and
            all(isinstance(value, str) for value in command["argv"]) and
            isinstance(command["cwd"], str) and Path(command["cwd"]).is_absolute(),
            "matrix command schema: " + label)
    payload = dict(command)
    digest = payload.pop("command_sha256")
    require(sha256_bytes(canonical_bytes(payload)) == digest and
            all(re.fullmatch(r"[0-9a-f]{64}", command[key]) is not None
                for key in ("stdout_sha256", "stderr_sha256")),
            "matrix command digest: " + label)


def validate_matrix_document(
        document, repo, candidate_commit, evidence_schema=SCHEMA):
    current = evidence_schema == SCHEMA
    expected_matrix_schema = (
        MATRIX_SCHEMA if current else MATRIX_SCHEMA_V1)
    expected_source_files = (
        MATRIX_SOURCE_FILES if current else
        PRE_XOR_MATRIX_SOURCE_FILES)
    expected_compile_source_counts = (
        MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS
        if current else
        PRE_XOR_MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS)
    matrix_compare_tests = (
        MATRIX_COMPARE_TESTS if current else PRE_XOR_MATRIX_COMPARE_TESTS)
    matrix_run_tests = (
        MATRIX_RUN_TESTS if current else PRE_XOR_MATRIX_RUN_TESTS)
    matrix_test_specs = (
        MATRIX_TEST_SPECS if current else PRE_XOR_MATRIX_TEST_SPECS)
    matrix_build_targets = (
        MATRIX_BUILD_TARGETS if current else PRE_XOR_MATRIX_BUILD_TARGETS)
    matrix_build_cache_keys = (
        MATRIX_BUILD_CACHE_KEYS if current else
        PRE_XOR_MATRIX_BUILD_CACHE_KEYS)
    backend_failure_ctest_regex = (
        MATRIX_BACKEND_FAILURE_CTEST_REGEX if current else
        "^leopard2_backend_failure_")
    portable_ctest_regex = (
        MATRIX_PORTABLE_CTEST_REGEX if current else
        "^leopard2_portable_isa$")
    cuda_ctest_regex = (
        MATRIX_CUDA_CTEST_REGEX if current else "^leopard2_cuda_optional$")
    require(set(document) == {
        "c_compiler", "compiler", "generator", "jobs", "jobs_per_variant",
        "machine", "mismatches",
        "schema", "source_changed_during_run", "source_fingerprint",
        "status", "variant_workers", "variants"},
        "matrix top-level key set")
    require(document.get("schema") == expected_matrix_schema,
            "matrix schema")
    require(document.get("status") == "passed", "matrix status")
    require(document.get("source_changed_during_run") is False,
            "matrix source changed during run")
    fingerprint = document.get("source_fingerprint", {})
    files = fingerprint.get("files")
    require(isinstance(files, dict), "matrix source files")
    require(set(files) == set(expected_source_files),
            "matrix source fingerprint file set")
    for relative, expected_digest in files.items():
        require(re.fullmatch(r"[0-9a-f]{64}", expected_digest or "") is not None,
                "matrix source digest format")
        require(sha256_bytes(git_blob(repo, candidate_commit, relative)) == expected_digest,
                "matrix source digest mismatch: " + relative)
    require(sha256_bytes(canonical_bytes(files)) == fingerprint.get("digest"),
            "matrix source aggregate digest")
    c_compiler = document["c_compiler"]
    compiler = document["compiler"]
    validate_matrix_compiler(c_compiler, "matrix C")
    validate_matrix_compiler(compiler, "matrix C++")
    require(document["generator"] in ("Ninja", "Unix Makefiles"),
            "matrix generator identity")
    require(all(isinstance(document[key], int) and document[key] > 0
                for key in ("jobs", "jobs_per_variant", "variant_workers")),
            "matrix worker geometry")
    machine = document["machine"]
    require(set(machine) == {"allowed_cpu_list", "architecture", "cpu_flags",
                            "logical_cpus_allowed", "platform"} and
            isinstance(machine["logical_cpus_allowed"], int) and
            machine["logical_cpus_allowed"] > 0,
            "matrix machine identity")
    variants = document.get("variants")
    require(isinstance(variants, list), "matrix variants")
    variant_names = [value.get("variant") for value in variants]
    require(sorted(variant_names) == ["auto", "avx2", "scalar", "ssse3"],
            "matrix must contain exactly auto/scalar/ssse3/avx2")
    by_variant = {value.get("variant"): value for value in variants}
    for value in variants:
        variant = value.get("variant")
        expected_variant_keys = {
            "configuration_id", "resumed", "schema", "source_fingerprint",
            "variant", "fresh_build", "commands", "expected_runtime_backend",
            "pin_cpu", "reason", "selected_cache_variant", "status", "tests",
            "build_environment", "build_identity",
        }
        require(set(value) == expected_variant_keys,
                "matrix variant key set: " + str(variant))
        require(value.get("schema") == expected_matrix_schema and
                value.get("status") == "passed" and
                value.get("selected_cache_variant") == variant and
                value.get("reason") == "" and
                value.get("resumed") is False and
                isinstance(value.get("pin_cpu"), int) and
                isinstance(value.get("commands"), list) and value["commands"],
                "matrix variant failed or misconfigured")
        require(value.get("expected_runtime_backend") ==
                (None if variant == "auto" else variant),
                "matrix expected runtime backend: " + str(variant))
        expected_environment = {
            "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
            "OMP_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"}
        if variant != "auto":
            expected_environment["LEO2_EXPECT_BACKEND"] = variant
        require(value.get("build_environment") == expected_environment,
                "matrix build/test environment: " + str(variant))
        configuration_input = {
            "c_compiler": c_compiler, "compiler": compiler,
            "generator": document["generator"],
            "jobs_per_variant": document["jobs_per_variant"],
            "machine": machine, "source": fingerprint, "variant": variant,
            "environment": expected_environment,
        }
        expected_configuration_id = sha256_bytes(
            canonical_bytes(configuration_input))
        require(value["configuration_id"] == expected_configuration_id,
                "matrix configuration identity: " + str(variant))
        identity = value.get("build_identity")
        require(isinstance(identity, dict) and set(identity) == {
                    "cache", "cache_sha256", "compile_commands",
                    "compile_commands_sha256", "test_executables", "tools",
                    "digest"} and
                sha256_bytes(canonical_bytes(identity["cache"])) ==
                identity["cache_sha256"] and
                sha256_bytes(canonical_bytes(identity["compile_commands"])) ==
                identity["compile_commands_sha256"],
                "matrix exact build identity: " + str(variant))
        identity_payload = dict(identity)
        identity_digest = identity_payload.pop("digest")
        require(sha256_bytes(canonical_bytes(identity_payload)) == identity_digest,
                "matrix build identity digest: " + str(variant))
        cache = identity["cache"]
        require(set(cache) == set(matrix_build_cache_keys) and
                cache.get("CMAKE_BUILD_TYPE") == "Release" and
                cache.get("CMAKE_GENERATOR") == document["generator"] and
                cache.get("CMAKE_C_FLAGS") == "" and
                cache.get("CMAKE_C_FLAGS_RELEASE") == "-O3 -DNDEBUG" and
                cache.get("CMAKE_CXX_FLAGS") == "" and
                cache.get("CMAKE_CXX_FLAGS_RELEASE") == "-O3 -DNDEBUG" and
                cache.get("CMAKE_EXE_LINKER_FLAGS") == "" and
                cache.get("CMAKE_EXE_LINKER_FLAGS_RELEASE") == "" and
                cache.get("CMAKE_STATIC_LINKER_FLAGS") == "" and
                cache.get("CMAKE_STATIC_LINKER_FLAGS_RELEASE") == "" and
                cache.get("ENABLE_OPENMP") == "ON" and
                cache.get("CMAKE_C_COMPILER") == "@tool/cc" and
                cache.get("CMAKE_CXX_COMPILER") == "@tool/cxx" and
                cache.get("LEO2_BACKEND_VARIANT") == variant and
                cache.get("LEO2_BUILD_TESTS") in ("ON", "1", "TRUE") and
                cache.get("LEO2_BUILD_BENCHMARKS") in ("OFF", "0", "FALSE", "") and
                cache.get("LEO2_BUILD_FUZZERS") in ("OFF", "0", "FALSE", "") and
                cache.get("LEO2_ENABLE_CUDA") in ("OFF", "0", "FALSE", "") and
                isinstance(identity["compile_commands"], list) and
                identity["compile_commands"],
                "matrix normalized CMake identity: " + str(variant))
        tools = identity["tools"]
        require(set(tools) == {"cmake", "ctest", "cc", "cxx"},
                "matrix build tool identity set")
        for tool in tools.values():
            require(set(tool) == {"basename", "binary_sha256", "version_sha256"} and
                    all(re.fullmatch(r"[0-9a-f]{64}", tool[key]) is not None
                        for key in ("binary_sha256", "version_sha256")),
                    "matrix build tool identity")
        require(tools["cc"]["basename"] == Path(c_compiler["executable"]).name and
                tools["cc"]["binary_sha256"] == c_compiler["binary_sha256"] and
                tools["cc"]["version_sha256"] == c_compiler["version_sha256"] and
                tools["cxx"]["basename"] == Path(compiler["executable"]).name and
                tools["cxx"]["binary_sha256"] == compiler["binary_sha256"] and
                tools["cxx"]["version_sha256"] == compiler["version_sha256"],
                "matrix compiler/tool crosslink")

        compile_commands = identity["compile_commands"]
        counts = {}
        for command in compile_commands:
            require(set(command) == {"file", "language", "argv"} and
                    isinstance(command["file"], str) and
                    command["file"] in expected_compile_source_counts and
                    command["language"] ==
                    ("C" if command["file"].endswith(".c") else "CXX") and
                    isinstance(command["argv"], list) and command["argv"] and
                    command["argv"][0] ==
                    ("@tool/cc" if command["language"] == "C" else "@tool/cxx") and
                    command["argv"].count("@source/" + command["file"]) == 1,
                    "matrix compile command identity")
            counts[command["file"]] = counts.get(command["file"], 0) + 1
        require(counts == expected_compile_source_counts,
                "matrix compile source multiset")
        matrix_backend_failure_tests = (
            MATRIX_BASE_BACKEND_FAILURE_TESTS +
            (MATRIX_AVX512_BACKEND_FAILURE_TESTS
             if current and "Leopard2BackendAVX512.cpp" in counts else ())
            if current else PRE_XOR_MATRIX_BACKEND_FAILURE_TESTS)

        commands = value["commands"]
        require(len(commands) == 2 + len(matrix_run_tests) + 2 +
                (1 if variant == "auto" else 0),
                "matrix command count: " + variant)
        configure = commands[0]
        validate_matrix_command(configure, "configure")
        configure_argv = configure["argv"]
        require(len(configure_argv) == 16 and
                Path(configure_argv[0]).name == tools["cmake"]["basename"],
                "matrix configure tool/shape: " + variant)
        source_root = configure_argv[2]
        build_root = configure_argv[4]
        expected_configure = [
            configure_argv[0], "-S", source_root, "-B", build_root,
            "-G", document["generator"], "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
            "-DCMAKE_C_COMPILER={}".format(c_compiler["executable"]),
            "-DCMAKE_CXX_COMPILER={}".format(compiler["executable"]),
            "-DLEO2_BACKEND_VARIANT={}".format(variant),
            "-DLEO2_BUILD_TESTS=ON", "-DLEO2_BUILD_BENCHMARKS=OFF",
            "-DLEO2_BUILD_FUZZERS=OFF", "-DLEO2_ENABLE_CUDA=OFF"]
        require(configure_argv == expected_configure and
                configure["cwd"] == source_root,
                "matrix configure argv: " + variant)
        fresh = value.get("fresh_build")
        require(isinstance(fresh, dict) and set(fresh) == {
                    "configured_from_empty", "identity_sha256"} and
                fresh["configured_from_empty"] is True and
                fresh["identity_sha256"] == sha256_bytes(canonical_bytes({
                    "configuration_id": expected_configuration_id,
                    "configure_argv": expected_configure,
                    "environment": expected_environment})),
                "matrix fresh-build identity: " + variant)

        build_command = commands[1]
        validate_matrix_command(build_command, "build")
        require(build_command["cwd"] == source_root and
                build_command["argv"] == [
                    configure_argv[0], "--build", build_root, "--config", "Release",
                    "-j", str(document["jobs_per_variant"]), "--target"] +
                list(matrix_build_targets),
                "matrix build argv: " + variant)
        require(value.get("source_fingerprint") == fingerprint["digest"],
                "matrix variant source mismatch")
        tests = value.get("tests")
        required_tests = set(matrix_run_tests) | {
            "backend_failures", "portable_isa"}
        if value.get("variant") == "auto":
            required_tests.add("cuda_optional")
        require(isinstance(tests, dict) and set(tests) == required_tests,
                "matrix test set: " + str(value.get("variant")))
        executables = identity["test_executables"]
        require(set(executables) == set(matrix_run_tests),
                "matrix executable set: " + variant)
        command_index = 2
        for test_name in matrix_run_tests:
            test = tests[test_name]
            validate_matrix_command(test, "test_" + test_name,
                                    {"executable_sha256"})
            require(canonical_bytes(test) == canonical_bytes(commands[command_index]),
                    "matrix test/command crosslink: " + test_name)
            command_index += 1
            executable = executables[test_name]
            target, arguments = matrix_test_specs[test_name]
            require(set(executable) == {"path", "sha256"} and
                    executable["path"].startswith("@build/") and
                    Path(executable["path"]).name == target and
                    executable["sha256"] == test["executable_sha256"],
                    "matrix test executable identity: " + test_name)
            executable_path = str(Path(build_root) /
                                  executable["path"][len("@build/"):])
            require(test["cwd"] == source_root and len(test["argv"]) >= 4 and
                    Path(test["argv"][0]).name == "taskset" and
                    test["argv"][1:3] == ["-c", str(value["pin_cpu"])] and
                    test["argv"][3:] == [executable_path] + arguments,
                    "matrix test argv: " + test_name)

        failures = tests["backend_failures"]
        validate_matrix_command(
            failures, "test_backend_failures",
            {"ctest_executed", "ctest_executed_tests"})
        require(
            failures["ctest_executed"] is True and
            failures["ctest_executed_tests"] ==
                sorted(matrix_backend_failure_tests) and
            canonical_bytes(failures) ==
                canonical_bytes(commands[command_index]) and
            failures["cwd"] == source_root and len(failures["argv"]) >= 4 and
            Path(failures["argv"][0]).name == "taskset" and
            failures["argv"][1:3] == ["-c", str(value["pin_cpu"])] and
            Path(failures["argv"][3]).name == tools["ctest"]["basename"] and
            failures["argv"][4:] == [
                "--test-dir", build_root, "-C", "Release", "-R",
                backend_failure_ctest_regex, "--output-on-failure"],
            "matrix backend-failure CTest command")
        command_index += 1

        portable = tests["portable_isa"]
        validate_matrix_command(portable, "test_portable_isa",
                                {"ctest_executed"})
        require(portable["ctest_executed"] is True and
                canonical_bytes(portable) == canonical_bytes(commands[command_index]) and
                portable["cwd"] == source_root and
                Path(portable["argv"][0]).name == tools["ctest"]["basename"] and
                portable["argv"][1:] == ["--test-dir", build_root, "-C", "Release",
                    "-R", portable_ctest_regex, "--output-on-failure"],
                "matrix portable-ISA CTest command")
        command_index += 1
        if variant == "auto":
            cuda = tests["cuda_optional"]
            validate_matrix_command(cuda, "test_cuda_optional",
                                    {"ctest_executed"})
            require(cuda["ctest_executed"] is True and
                    canonical_bytes(cuda) == canonical_bytes(commands[command_index]) and
                    cuda["cwd"] == source_root and
                    cuda["argv"] == [portable["argv"][0], "--test-dir", build_root,
                        "-C", "Release", "-R", cuda_ctest_regex,
                        "--output-on-failure"],
                    "matrix CUDA-optional CTest command")
            command_index += 1
        require(command_index == len(commands), "matrix command closure: " + variant)

    recomputed_mismatches = []
    auto_tests = by_variant["auto"]["tests"]
    for variant in sorted(set(by_variant) - {"auto"}):
        for test_name in matrix_compare_tests:
            for stream in ("stdout", "stderr"):
                key = stream + "_sha256"
                expected_digest = auto_tests[test_name][key]
                actual_digest = by_variant[variant]["tests"][test_name][key]
                if actual_digest != expected_digest:
                    recomputed_mismatches.append({
                        "actual": actual_digest, "expected": expected_digest,
                        "stream": stream, "test": test_name,
                        "variant": variant,
                    })
    require(document.get("mismatches") == recomputed_mismatches,
            "matrix mismatch set does not replay from test hashes")
    require(recomputed_mismatches == [], "matrix backend outputs differ")


def validate_declared_template_paths(args, build, schema=SCHEMA):
    cache = Path(getattr(args, build + "_cmake_cache")).resolve()
    root = cache.parent
    cmake_identity = cmake_identity_for_schema(schema)
    expected = {
        build + "_compile_commands": root / "compile_commands.json",
        build + "_library": root / cmake_identity["archive"],
        build: root / "bench_leopard2",
    }
    for attribute, expected_path in expected.items():
        supplied = Path(getattr(args, attribute)).resolve()
        require(supplied == expected_path.resolve(),
                "{} path substitution: expected canonical target path {}".format(
                    attribute, expected_path))
        require(supplied.is_file(), attribute + " template artifact is missing")


def matrix_record(path, repo, candidate_commit, evidence_schema=SCHEMA):
    raw = Path(path).read_bytes()
    require(raw, "matrix artifact is empty")
    document = parse_json_bytes(raw, "backend matrix")
    validate_matrix_document(
        document, repo, candidate_commit, evidence_schema)
    return raw, {
        "sha256": sha256_bytes(raw),
        "source_fingerprint": document["source_fingerprint"]["digest"],
        "variant_count": len(document["variants"]),
    }


def pair_lease_payload(cpu, sibling, uid=None):
    retained_uid = os.getuid() if uid is None else uid
    require(isinstance(cpu, int) and not isinstance(cpu, bool) and cpu >= 0 and
            isinstance(sibling, int) and not isinstance(sibling, bool) and
            sibling >= 0 and cpu != sibling and
            isinstance(retained_uid, int) and not isinstance(retained_uid, bool) and
            retained_uid >= 0,
            "pair lease requires distinct non-negative CPUs and UID")
    return {
        "cpus": sorted((cpu, sibling)),
        "schema": PAIR_LEASE_SCHEMA,
        "uid": retained_uid,
    }


def pair_lease_name(cpu, sibling, uid=None):
    payload = pair_lease_payload(cpu, sibling, uid)
    return "leopard2-cpu-pair-{}-{}-{}.lock".format(
        payload["uid"], payload["cpus"][0], payload["cpus"][1])


def pair_lease_directory(uid=None):
    retained_uid = os.getuid() if uid is None else uid
    return pair_lease_runtime_root(retained_uid) / "leopard2-cpu-leases"


class PairLease(object):
    """Serialize legacy and current evidence runners by physical CPU pair."""

    def __init__(self, cpu, sibling, root=None):
        pair_lease_payload(cpu, sibling)
        self.cpu = cpu
        self.sibling = sibling
        self.production_root = root is None
        self.root = pair_lease_directory() if root is None else Path(root)
        self.path = self.root / pair_lease_name(cpu, sibling)
        self.descriptor = None
        self.identity = None
        self.kernel_socket = None
        material = canonical_bytes({
            "cpus": sorted((cpu, sibling)),
            "root": os.path.abspath(os.fspath(self.root)),
            "schema": PAIR_LEASE_SCHEMA,
            "uid": os.getuid(),
        })
        self.kernel_name = b"\0leopard2-pair-v1-" + \
            hashlib.sha256(material).hexdigest()[:40].encode("ascii")

    def _acquire_kernel_lease(self):
        require(sys.platform.startswith("linux") and hasattr(socket, "AF_UNIX"),
                "Linux abstract Unix sockets are required for stable CPU leases")
        lease = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        lease.set_inheritable(False)
        try:
            lease.bind(self.kernel_name)
        except OSError as error:
            lease.close()
            if error.errno == errno.EADDRINUSE:
                raise EvidenceError(
                    "physical CPU pair already has a kernel lease")
            raise EvidenceError(
                "cannot bind stable CPU pair lease: {}".format(error))
        self.kernel_socket = lease

    def _release_kernel_lease(self):
        if self.kernel_socket is not None:
            self.kernel_socket.close()
            self.kernel_socket = None

    def _validate_directory(self):
        if self.production_root:
            runtime = os.lstat(self.root.parent)
            require(stat.S_ISDIR(runtime.st_mode) and
                    runtime.st_uid == os.getuid() and
                    stat.S_IMODE(runtime.st_mode) == 0o700,
                    "CPU pair runtime directory is not owned mode-0700")
        directory = os.lstat(self.root)
        require(stat.S_ISDIR(directory.st_mode) and
                directory.st_uid == os.getuid() and
                stat.S_IMODE(directory.st_mode) == 0o700,
                "CPU pair lease directory is not owned mode-0700")
        return directory

    def validate_current(self):
        require(self.descriptor is not None and self.identity is not None and
                self.kernel_socket is not None and
                self.kernel_socket.fileno() >= 0 and
                self.kernel_socket.getsockname() == self.kernel_name,
                "CPU pair lease is not held")
        try:
            directory = self._validate_directory()
            descriptor = os.fstat(self.descriptor)
            path = os.lstat(self.path)
        except OSError as error:
            raise EvidenceError(
                "CPU pair lease path/descriptor revalidation failed: {}".format(
                    error))
        require(stat.S_ISREG(descriptor.st_mode) and
                descriptor.st_uid == os.getuid() and descriptor.st_nlink == 1 and
                stat.S_IMODE(descriptor.st_mode) == 0o600 and
                (descriptor.st_dev, descriptor.st_ino) ==
                (path.st_dev, path.st_ino) ==
                (self.identity["device"], self.identity["inode"]),
                "CPU pair lease path was replaced or changed")
        require((directory.st_dev, directory.st_ino) ==
                (self.identity["directory_device"],
                 self.identity["directory_inode"]),
                "CPU pair lease directory was replaced")
        expected = canonical_bytes(pair_lease_payload(self.cpu, self.sibling))
        os.lseek(self.descriptor, 0, os.SEEK_SET)
        require(os.read(self.descriptor, 4096) == expected,
                "CPU pair lease contents changed while held")
        return self.identity

    def __enter__(self):
        require(fcntl is not None and sys.platform.startswith("linux"),
                "stable CPU pair leases require Linux fcntl")
        self._acquire_kernel_lease()
        try:
            try:
                self.root.mkdir(mode=0o700)
                os.chmod(self.root, 0o700)
            except FileExistsError:
                pass
            self._validate_directory()
            flags = os.O_RDWR | getattr(os, "O_CLOEXEC", 0) | \
                getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
            created = False
            before = None
            try:
                before = os.lstat(self.path)
                require(stat.S_ISREG(before.st_mode),
                        "CPU pair lease path is not regular")
                self.descriptor = os.open(self.path, flags)
            except FileNotFoundError:
                self.descriptor = os.open(
                    self.path, flags | os.O_CREAT | os.O_EXCL, 0o600)
                created = True
            if created:
                os.fchmod(self.descriptor, 0o600)
            metadata = os.fstat(self.descriptor)
            path = os.lstat(self.path)
            require(stat.S_ISREG(metadata.st_mode) and
                    metadata.st_uid == os.getuid() and metadata.st_nlink == 1 and
                    stat.S_IMODE(metadata.st_mode) == 0o600 and
                    (metadata.st_dev, metadata.st_ino) ==
                    (path.st_dev, path.st_ino) and
                    (before is None or (metadata.st_dev, metadata.st_ino) ==
                     (before.st_dev, before.st_ino)),
                    "CPU pair lease file has unsafe type, links, mode, or identity")
            try:
                fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise EvidenceError("physical CPU pair is already leased") from error
            expected = canonical_bytes(pair_lease_payload(self.cpu, self.sibling))
            os.lseek(self.descriptor, 0, os.SEEK_SET)
            retained = os.read(self.descriptor, 4096)
            if not retained:
                require(os.write(self.descriptor, expected) == len(expected),
                        "short CPU pair lease write")
                os.fsync(self.descriptor)
                retained = expected
            require(retained == expected,
                    "CPU pair lease has unexpected contents")
            directory = self._validate_directory()
            self.identity = {
                "device": metadata.st_dev,
                "directory_device": directory.st_dev,
                "directory_inode": directory.st_ino,
                "inode": metadata.st_ino,
                "lock": "exclusive_nonblocking_pair_wide",
                "path": str(self.path.absolute()),
                "payload": pair_lease_payload(self.cpu, self.sibling),
                "sha256": sha256_bytes(retained),
            }
            return self.validate_current()
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, exc_type, exc, tb):
        descriptor = self.descriptor
        self.descriptor = None
        self.identity = None
        try:
            if descriptor is not None:
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(descriptor)
        finally:
            self._release_kernel_lease()


def pair_lease_runtime_root(uid=None):
    """Return the user-owned, root-anchored runtime directory shared by runners."""
    retained_uid = os.getuid() if uid is None else uid
    return Path("/run/user") / str(retained_uid)


class StableLeaseAnchor(object):
    """Serialize current Leopard2 evidence runners across replaceable files."""

    def __init__(self, path=None):
        self.path = (pair_lease_runtime_root() if path is None else
                     Path(os.path.abspath(os.fspath(path))))
        self.descriptor = None
        self.identity = None

    @staticmethod
    def _metadata(metadata):
        mode = stat.S_IMODE(metadata.st_mode)
        require(stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and mode & 0o022 == 0,
                "stable runner lease anchor has unsafe ownership or mode")
        return (metadata.st_dev, metadata.st_ino, metadata.st_uid, mode)

    def validate_current(self):
        require(self.descriptor is not None and self.identity is not None,
                "stable runner lease anchor is not held")
        descriptor = os.fstat(self.descriptor)
        path = os.lstat(self.path)
        require(self._metadata(descriptor) == self.identity and
                (descriptor.st_dev, descriptor.st_ino) ==
                (path.st_dev, path.st_ino),
                "stable runner lease anchor path was replaced")

    def acquire(self):
        require(fcntl is not None,
                "evidence collection requires Linux/POSIX fcntl file locking")
        require(self.descriptor is None,
                "stable runner lease anchor is already held")
        flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                 getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0))
        try:
            self.descriptor = os.open(self.path, flags)
            fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            descriptor = os.fstat(self.descriptor)
            path = os.lstat(self.path)
            self.identity = self._metadata(descriptor)
            require((descriptor.st_dev, descriptor.st_ino) ==
                    (path.st_dev, path.st_ino),
                    "stable runner lease anchor changed during acquisition")
            self.validate_current()
            return self
        except BaseException as error:
            self.release()
            if isinstance(error, OSError):
                raise EvidenceError(
                    "cannot acquire stable runner lease anchor: {}".format(error))
            raise

    def release(self):
        try:
            if self.descriptor is not None:
                descriptor = self.descriptor
                self.descriptor = None
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_UN)
                finally:
                    os.close(descriptor)
        finally:
            self.identity = None


class ReservationHandle(object):
    """Release the reservation, stable anchor, and legacy pair lease together."""

    def __init__(self, handle, anchor, pair):
        self.handle = handle
        self.anchor = anchor
        self.pair = pair

    def validate_current(self):
        require(self.handle is not None and not self.handle.closed,
                "CPU reservation is not held")
        require(self.anchor is not None and self.pair is not None,
                "CPU reservation guards are not held")
        self.anchor.validate_current()
        self.pair.validate_current()

    def close(self):
        try:
            if self.handle is not None:
                handle = self.handle
                self.handle = None
                try:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
                finally:
                    handle.close()
        finally:
            try:
                if self.anchor is not None:
                    anchor = self.anchor
                    self.anchor = None
                    anchor.release()
            finally:
                if self.pair is not None:
                    pair = self.pair
                    self.pair = None
                    pair.__exit__(None, None, None)


def reservation_record(path, cpu, sibling, runtime_root=None):
    require(fcntl is not None,
            "evidence collection requires Linux/POSIX fcntl file locking")
    pair_root = None if runtime_root is None else \
        Path(runtime_root) / "leopard2-cpu-leases"
    pair = PairLease(cpu, sibling, root=pair_root)
    anchor = StableLeaseAnchor(runtime_root)
    handle = None
    try:
        pair.__enter__()
        anchor.acquire()
        path = Path(path).resolve()
        require(path.is_file(), "reservation file missing")
        handle = path.open("r+")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            raise EvidenceError(
                "cannot acquire exclusive CPU reservation: {}".format(error))
        raw = handle.read().encode("utf-8")
        document = parse_json_bytes(raw, "reservation")
        canonical = checked_json_bytes(document, "reservation")
        require(raw == canonical,
                "reservation file must be canonical JSON with no trailing newline")
        require(set(document) == {"schema", "status", "benchmark_cpu",
                                 "reserved_sibling", "owner", "nonce"},
                "reservation key set")
        require(document.get("schema") == RESERVATION_SCHEMA,
                "reservation schema")
        require(document.get("status") == "held", "reservation status")
        require(document.get("benchmark_cpu") == cpu and
                document.get("reserved_sibling") == sibling,
                "reservation CPU identity")
        require(isinstance(document.get("owner"), str) and document["owner"].strip(),
                "reservation owner")
        require(isinstance(document.get("nonce"), str) and
                len(document["nonce"]) >= 16,
                "reservation nonce")
        anchor.validate_current()
        pair.validate_current()
        return (ReservationHandle(handle, anchor, pair),
                {"sha256": sha256_bytes(canonical), "document": document})
    except BaseException:
        if handle is not None:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            finally:
                handle.close()
        try:
            anchor.release()
        finally:
            pair.__exit__(None, None, None)
        raise


def validate_execution(execution):
    expected_keys = {
        "initial_affinity", "enforced_affinity", "benchmark_cpu",
        "reserved_sibling", "topology", "reservation", "child_environment",
        "timeout_seconds", "host",
    }
    require(set(execution) == expected_keys, "execution key set")
    cpu = execution["benchmark_cpu"]
    sibling = execution["reserved_sibling"]
    initial = execution["initial_affinity"]
    require(isinstance(cpu, int) and not isinstance(cpu, bool) and cpu >= 0 and
            isinstance(sibling, int) and not isinstance(sibling, bool) and sibling >= 0 and
            isinstance(initial, list) and
            all(isinstance(value, int) and not isinstance(value, bool) and value >= 0
                for value in initial) and
            initial == sorted(set(initial)) and cpu in initial and sibling in initial,
            "initial affinity does not contain the isolated core pair")
    require(execution["enforced_affinity"] == [cpu],
            "enforced affinity is not the benchmark CPU singleton")
    topology = execution["topology"]
    require(set(topology) == {"raw", "sha256", "cpus"}, "topology keys")
    require(sha256_bytes(topology["raw"].encode("ascii")) == topology["sha256"],
            "topology raw digest")
    require(parse_cpu_list(topology["raw"]) == topology["cpus"],
            "topology raw parse")
    require(cpu in topology["cpus"] and sibling in topology["cpus"] and cpu != sibling,
            "benchmark/reserved CPUs are not distinct topology siblings")
    reservation = execution["reservation"]
    require(set(reservation) == {"sha256", "document"}, "reservation keys")
    require(sha256_bytes(canonical_bytes(reservation["document"])) ==
            reservation["sha256"], "reservation digest")
    document = reservation["document"]
    require(set(document) == {"schema", "status", "benchmark_cpu",
                             "reserved_sibling", "owner", "nonce"} and
            document.get("schema") == RESERVATION_SCHEMA and
            document.get("status") == "held" and
            document.get("benchmark_cpu") == cpu and
            document.get("reserved_sibling") == sibling and
            isinstance(document.get("owner"), str) and document["owner"].strip() and
            isinstance(document.get("nonce"), str) and len(document["nonce"]) >= 16,
            "reservation document")
    require(execution["child_environment"] == {
        "LC_ALL": "C", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
        "PATH": "/usr/bin:/bin"},
        "benchmark child environment")
    validate_host_record(execution["host"], cpu)
    require(isinstance(execution["timeout_seconds"], int) and
            execution["timeout_seconds"] > 0, "timeout")


def encode_raw_record(stdout, stderr):
    return {
        "stdout_base64": base64.b64encode(stdout).decode("ascii"),
        "stdout_sha256": sha256_bytes(stdout),
        "stderr_base64": base64.b64encode(stderr).decode("ascii"),
        "stderr_sha256": sha256_bytes(stderr),
    }


def decode_raw_record(record, label):
    require(set(record) == {
        "stdout_base64", "stdout_sha256", "stderr_base64", "stderr_sha256"},
        label + " raw record keys")
    try:
        stdout = base64.b64decode(record["stdout_base64"], validate=True)
        stderr = base64.b64decode(record["stderr_base64"], validate=True)
    except (ValueError, TypeError) as error:
        raise EvidenceError(label + " invalid base64: " + str(error))
    require(sha256_bytes(stdout) == record["stdout_sha256"],
            label + " stdout mutation")
    require(sha256_bytes(stderr) == record["stderr_sha256"],
            label + " stderr mutation")
    return stdout, stderr


def validate_manifest(manifest_path, repo, raw_bundle_path=None,
                      binaries=None, matrix_path=None, allow_self_test=False):
    manifest_path = Path(manifest_path)
    manifest = read_json(manifest_path, "ABBA manifest")
    require(set(manifest) == {
        "schema", "status", "provenance", "campaign", "entries",
        "raw_bundle_file", "raw_bundle_sha256", "raw_evidence_sha256",
        "summary"}, "manifest top-level key set")
    manifest_schema = manifest.get("schema")
    require(isinstance(manifest_schema, str) and
            manifest_schema in SCHEMA_TO_CMAKE_IDENTITY,
            "unsupported manifest schema")
    require(manifest.get("status") in ("passed", "failed"),
            "campaign status")
    campaign = manifest.get("campaign", {})
    require(set(campaign) == {
        "order_per_round", "rounds", "cell_count", "entry_count",
        "samples_per_invocation", "warmups_per_invocation", "reuse_per_sample",
        "batch", "threads", "seed", "target_threshold_percent",
        "neighbor_floor_percent", "requested_backend", "resolved_backend",
        "self_test", "started_unix",
        "finished_unix"}, "campaign key set")
    is_self_test = campaign.get("self_test") is True
    require(not is_self_test or allow_self_test,
            "self-test manifest is not production evidence")
    requested_backend = campaign.get("requested_backend")
    resolved_backend = campaign.get("resolved_backend")
    require(requested_backend in SUPPORTED_BACKENDS and
            resolved_backend == requested_backend,
            "campaign requested/resolved backend identity")
    provenance = manifest.get("provenance", {})
    require(set(provenance) == {
        "runner_sha256", "git", "builds", "matrix", "execution"},
        "provenance key set")
    require(set(provenance["git"]) == {"baseline", "candidate"},
            "Git provenance key set")
    require(set(provenance["builds"]) == {"baseline", "candidate"},
            "build provenance key set")
    require(set(provenance["matrix"]) == {
        "sha256", "source_fingerprint", "variant_count"},
        "matrix provenance key set")
    for build in ("baseline", "candidate"):
        validate_git_record(repo, provenance["git"][build])
        validate_build_record(
            provenance["builds"][build], repo, manifest_schema)
        require(provenance["builds"][build]["source"] ==
                provenance["git"][build],
                build + " build/source Git binding")
    require(provenance["builds"]["baseline"]["tools"] ==
            provenance["builds"]["candidate"]["tools"],
            "baseline/candidate toolchain identity mismatch")
    require(provenance["builds"]["baseline"]["configuration"] ==
            provenance["builds"]["candidate"]["configuration"],
            "baseline/candidate configuration mismatch")
    if is_self_test:
        require(sha256_file(Path(__file__).resolve()) == provenance["runner_sha256"],
                "self-test runner hash")
    else:
        runner_blob = git_blob(
            repo, provenance["git"]["candidate"]["commit"],
            "experiments/leopard2/backend_butterfly/run_abba.py")
        require(sha256_bytes(runner_blob) == provenance["runner_sha256"],
                "runner is not bound to candidate commit")
    validate_execution(provenance["execution"])
    require(re.fullmatch(r"[0-9a-f]{64}", provenance["runner_sha256"] or "")
            is not None, "runner digest format")
    require(campaign["order_per_round"] == [value[0] for value in SEQUENCES] and
            campaign["rounds"] == len(ROUNDS) and
            campaign["cell_count"] == len(CELLS) and
            campaign["entry_count"] == len(expected_jobs()) and
            campaign["samples_per_invocation"] == 7 and
            campaign["warmups_per_invocation"] == 3 and
            campaign["reuse_per_sample"] == 8 and campaign["batch"] == 1 and
            campaign["threads"] == 1 and campaign["seed"] == 42 and
            campaign["target_threshold_percent"] == TARGET_THRESHOLD and
            campaign["neighbor_floor_percent"] == NEIGHBOR_FLOOR,
            "campaign geometry")
    started = finite_number(campaign["started_unix"], "campaign start")
    finished = finite_number(campaign["finished_unix"], "campaign finish")
    require(0.0 < started <= finished, "campaign time interval")

    require(manifest["raw_bundle_file"] == "abba_raw.json",
            "raw bundle must use the canonical adjacent filename")
    for key in ("raw_bundle_sha256", "raw_evidence_sha256"):
        require(re.fullmatch(r"[0-9a-f]{64}", manifest[key] or "") is not None,
                key + " format")
    if raw_bundle_path is None:
        raw_bundle_path = manifest_path.parent / manifest["raw_bundle_file"]
    raw_bundle_path = Path(raw_bundle_path)
    require(raw_bundle_path.is_file(), "adjacent raw evidence bundle missing")
    raw_bundle_bytes = raw_bundle_path.read_bytes()
    require(sha256_bytes(raw_bundle_bytes) == manifest["raw_bundle_sha256"],
            "raw bundle digest mismatch")
    bundle = parse_json_bytes(raw_bundle_bytes, "ABBA raw bundle")
    require(bundle.get("schema") == RAW_SCHEMA, "raw bundle schema")
    require(set(bundle) == {"schema", "matrix_base64", "matrix_sha256", "raw"},
            "raw bundle key set")
    matrix_raw = base64.b64decode(bundle["matrix_base64"], validate=True)
    require(sha256_bytes(matrix_raw) == bundle["matrix_sha256"] ==
            provenance["matrix"]["sha256"], "embedded matrix digest")
    matrix_document = parse_json_bytes(matrix_raw, "embedded backend matrix")
    validate_matrix_document(
        matrix_document, repo, provenance["git"]["candidate"]["commit"],
        manifest_schema)
    matrix_backend = next(
        (value for value in matrix_document["variants"]
         if value["variant"] == requested_backend), None)
    require(matrix_backend is not None and
            matrix_backend["status"] == "passed" and
            matrix_backend["expected_runtime_backend"] == resolved_backend,
            "campaign backend is not bound to a passing matrix variant")
    candidate_cc = provenance["builds"]["candidate"]["tools"]["cc"]
    candidate_cxx = provenance["builds"]["candidate"]["tools"]["cxx"]
    matrix_cc = matrix_document["c_compiler"]
    matrix_cxx = matrix_document["compiler"]
    require(Path(matrix_cc["executable"]).name == candidate_cc["basename"] and
            matrix_cc["binary_sha256"] == candidate_cc["binary_sha256"] and
            compiler_version_signature(matrix_cc["version"]) ==
            compiler_version_signature(candidate_cc["version"]) and
            Path(matrix_cxx["executable"]).name == candidate_cxx["basename"] and
            matrix_cxx["binary_sha256"] == candidate_cxx["binary_sha256"] and
            compiler_version_signature(matrix_cxx["version"]) ==
            compiler_version_signature(candidate_cxx["version"]),
            "backend matrix/C and C++ compiler build identity mismatch")
    require(matrix_document["source_fingerprint"]["digest"] ==
            provenance["matrix"]["source_fingerprint"],
            "matrix source fingerprint record")
    require(provenance["matrix"]["variant_count"] == 4 and
            provenance["matrix"]["variant_count"] ==
            len(matrix_document["variants"]), "matrix variant count")
    if matrix_path is not None:
        require(sha256_file(matrix_path) == provenance["matrix"]["sha256"],
                "supplied matrix artifact digest")

    entries = manifest.get("entries", [])
    jobs = expected_jobs()
    require(len(entries) == len(jobs), "manifest entry count")
    require([entry.get("name") for entry in entries] ==
            [job[0] for job in jobs], "manifest ABBA execution order")
    expected = {name: (item, round_number, sequence, build)
                for name, item, round_number, sequence, build in jobs}
    require(set(bundle["raw"]) == set(expected), "raw bundle job geometry")
    seen = set()
    raw_by_name = {}
    missing_by_cell = {}
    reported_builds = {}
    binary_hashes = {
        build: provenance["builds"][build]["artifacts"]["benchmark_sha256"]
        for build in ("baseline", "candidate")}
    for entry in entries:
        required_entry_keys = {
            "name", "cell", "round", "sequence", "build", "argv",
            "command_sha256", "binary_sha256", "stdout_sha256",
            "stderr_sha256", "returncode",
        }
        require(set(entry) == required_entry_keys, "entry key set")
        name = entry["name"]
        require(name in expected and name not in seen,
                "unexpected, duplicate, or relabelled entry: " + name)
        seen.add(name)
        item, round_number, sequence, build = expected[name]
        require((entry["cell"], entry["round"], entry["sequence"], entry["build"]) ==
                (item["name"], round_number, sequence, build),
                "entry geometry mismatch: " + name)
        require(entry["returncode"] == 0, "nonzero return code: " + name)
        require(entry["binary_sha256"] == binary_hashes[build],
                "binary relabel: " + name)
        expected_argv = [logical_executable(build, binary_hashes[build])] + \
            benchmark_arguments(item, requested_backend)
        require(entry["argv"] == expected_argv, "logical argv mismatch: " + name)
        command_record = {
            "affinity": provenance["execution"]["enforced_affinity"],
            "argv": expected_argv,
            "environment": provenance["execution"]["child_environment"],
        }
        require(entry["command_sha256"] ==
                sha256_bytes(canonical_bytes(command_record)),
                "command hash mismatch: " + name)
        stdout, stderr = decode_raw_record(bundle["raw"][name], name)
        require(entry["stdout_sha256"] == sha256_bytes(stdout) and
                entry["stderr_sha256"] == sha256_bytes(stderr),
                "entry/raw digest mismatch: " + name)
        try:
            raw = parse_json_bytes(stdout, "benchmark output " + name)
        except (UnicodeError, ValueError) as error:
            raise EvidenceError("invalid benchmark JSON {}: {}".format(name, error))
        parsed = check_raw(
            raw, item, name, requested_backend,
            missing_by_cell.get(item["name"]), manifest_schema)
        missing_by_cell.setdefault(item["name"], parsed[2])
        require(build not in reported_builds or reported_builds[build] == parsed[3],
                "benchmark build identity changed: " + name)
        reported_builds.setdefault(build, parsed[3])
        raw_by_name[name] = parsed[:2]
    require(seen == set(expected), "campaign geometry incomplete")
    require(reported_builds.get("baseline") == reported_builds.get("candidate"),
            "baseline/candidate reported compiler identity mismatch")
    require(stable_raw_digest(entries) == manifest["raw_evidence_sha256"],
            "stable raw evidence digest mismatch")
    recomputed = summarize(entries, raw_by_name, requested_backend)
    require(canonical_bytes(recomputed) == canonical_bytes(manifest["summary"]),
            "summary does not replay from raw evidence")
    policy_failures = []
    for item in recomputed:
        for metric in ("encode", "decode"):
            if not item[metric]["accepted"]:
                policy_failures.append(item["name"] + ":" + metric)
    if manifest["status"] == "passed":
        if policy_failures:
            raise EvidenceError(
                policy_failures[0] + " violates its promotion/neighbor floor")
    else:
        require(policy_failures,
                "failed campaign has no replayed promotion/neighbor-floor failure")
    if binaries is not None:
        for build in ("baseline", "candidate"):
            path = Path(binaries[build]).resolve()
            require(path.is_file(), build + " binary missing")
            require(sha256_file(path) == binary_hashes[build],
                    build + " supplied binary hash mismatch")
    return manifest


def run_campaign(
    args, repo, allow_dirty=False, self_test=False, evidence_schema=SCHEMA,
):
    require(evidence_schema == SCHEMA or self_test,
            "new production campaigns must use the current butterfly schema")
    cmake_identity_for_schema(evidence_schema)
    require(hasattr(os, "sched_getaffinity") and hasattr(os, "sched_setaffinity"),
            "Linux scheduler affinity APIs are required")
    require(isinstance(args.build_jobs, int) and 1 <= args.build_jobs <= 128,
            "build jobs must be in [1,128]")
    require(args.backend in SUPPORTED_BACKENDS,
            "campaign backend must be one of: " +
            ",".join(SUPPORTED_BACKENDS))
    initial_affinity = sorted(os.sched_getaffinity(0))
    require(args.cpu in initial_affinity and
            args.reserved_sibling in initial_affinity,
            "isolated core pair is not initially allowed")
    if not self_test:
        # Fresh provenance rebuilds precede the measured campaign.  Keep every
        # configure/compiler/link descendant off both reserved SMT siblings so
        # setup cannot perturb the core immediately before timing.  Restore the
        # declared launch mask after the rebuild; the measured phase below then
        # narrows it to the benchmark CPU as before.
        setup_affinity = set(initial_affinity) - {
            args.cpu, args.reserved_sibling}
        require(setup_affinity,
                "no housekeeping CPU remains outside the reserved pair")
        os.sched_setaffinity(0, setup_affinity)
        require(set(os.sched_getaffinity(0)) == setup_affinity,
                "failed to isolate provenance rebuilds from reserved pair")
    source_roots = {
        "baseline": Path(args.baseline_source_root).resolve(),
        "candidate": Path(args.candidate_source_root).resolve(),
    }
    commits = {"baseline": args.baseline_commit,
               "candidate": args.candidate_commit}
    git_records = {
        build: git_record(source_roots[build], commits[build], not allow_dirty)
        for build in ("baseline", "candidate")}
    for build in ("baseline", "candidate"):
        validate_declared_template_paths(args, build, evidence_schema)
    shadow_context = tempfile.TemporaryDirectory(
        prefix="leo2-butterfly-fresh-builds-")
    shadow_root = Path(shadow_context.name)
    if self_test:
        isolated = {}
        for build in ("baseline", "candidate"):
            artifacts = getattr(args, build + "_self_test_artifacts")
            isolated[build] = {
                "compile_commands": Path(getattr(args, build + "_compile_commands")),
                "cmake_cache": Path(getattr(args, build + "_cmake_cache")),
                "library": Path(getattr(args, build + "_library")),
                "binary": Path(getattr(args, build)),
                "target_graph": artifacts,
                "rebuild": self_test_rebuild_record(
                    args.build_jobs, getattr(args, build + "_cmake_cache"),
                    evidence_schema),
            }
    else:
        isolated = {
            "baseline": fresh_rebuild(
                source_roots["baseline"], args.baseline_cmake_cache,
                shadow_root / "baseline", args.build_jobs),
            "candidate": fresh_rebuild(
                source_roots["candidate"], args.candidate_cmake_cache,
                shadow_root / "candidate", args.build_jobs),
        }
        # A clean build must not have changed either declared source tree.
        for build in ("baseline", "candidate"):
            require(git_record(source_roots[build], commits[build], True) ==
                    git_records[build], build + " source changed during rebuild")
    builds = {
        "baseline": build_record(
            source_roots["baseline"], git_records["baseline"],
            isolated["baseline"]["compile_commands"],
            isolated["baseline"]["cmake_cache"], isolated["baseline"]["library"],
            isolated["baseline"]["binary"], isolated["baseline"]["rebuild"],
            isolated["baseline"]["target_graph"], evidence_schema),
        "candidate": build_record(
            source_roots["candidate"], git_records["candidate"],
            isolated["candidate"]["compile_commands"],
            isolated["candidate"]["cmake_cache"], isolated["candidate"]["library"],
            isolated["candidate"]["binary"], isolated["candidate"]["rebuild"],
            isolated["candidate"]["target_graph"], evidence_schema),
    }
    binaries = {build: isolated[build]["binary"].resolve()
                for build in ("baseline", "candidate")}
    matrix_raw, matrix_info = matrix_record(
        args.matrix, repo, commits["candidate"], evidence_schema)
    if not self_test:
        os.sched_setaffinity(0, set(initial_affinity))
        require(sorted(os.sched_getaffinity(0)) == initial_affinity,
                "failed to restore declared launch affinity after rebuilds")
    topology = sibling_topology(args.cpu)
    require(args.reserved_sibling in topology["cpus"] and
            args.reserved_sibling != args.cpu,
            "reserved CPU is not a distinct topology sibling")
    reservation_handle = None
    try:
        reservation_handle, reservation = reservation_record(
            args.reservation_file, args.cpu, args.reserved_sibling)
        output_directory = args.output.resolve()
        require(not output_directory.exists() or not any(output_directory.iterdir()),
                "output directory must be absent or empty")
        output_directory.mkdir(parents=True, exist_ok=True)
        os.sched_setaffinity(0, {args.cpu})
        enforced_affinity = sorted(os.sched_getaffinity(0))
        require(enforced_affinity == [args.cpu], "failed singleton affinity")
        static_power = power_state(args.cpu)
        pre_frequency = current_frequency(args.cpu)
        child_environment = {
            "LC_ALL": "C", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
            "PATH": "/usr/bin:/bin"}
        execution = {
            "initial_affinity": initial_affinity,
            "enforced_affinity": enforced_affinity,
            "benchmark_cpu": args.cpu,
            "reserved_sibling": args.reserved_sibling,
            "topology": topology,
            "reservation": reservation,
            "child_environment": child_environment,
            "timeout_seconds": args.timeout,
        }
        entries = []
        raw_records = {}
        raw_by_name = {}
        missing_by_cell = {}
        reported_builds = {}
        jobs = expected_jobs()
        campaign_start = time.time()
        for index, (name, item, round_number, sequence, build) in enumerate(jobs, 1):
            require(sorted(os.sched_getaffinity(0)) == enforced_affinity,
                    "runner affinity changed before " + name)
            actual_argv = [str(binaries[build])] + \
                benchmark_arguments(item, args.backend)
            logical_argv = [logical_executable(
                build, builds[build]["artifacts"]["benchmark_sha256"])] + \
                benchmark_arguments(item, args.backend)
            command_record = {"affinity": enforced_affinity,
                              "argv": logical_argv,
                              "environment": child_environment}
            print("[{}/{}] {}".format(index, len(jobs), name),
                  file=sys.stderr, flush=True)
            try:
                completed = run_benchmark_bounded(
                    actual_argv, child_environment, args.timeout)
            except EvidenceError as error:
                raise EvidenceError(
                    "benchmark execution failed {}: {}".format(name, error))
            require(completed.returncode == 0,
                    "benchmark failed {} rc={}".format(name, completed.returncode))
            raw = parse_json_bytes(completed.stdout, "benchmark output " + name)
            parsed = check_raw(
                raw, item, name, args.backend,
                missing_by_cell.get(item["name"]), evidence_schema)
            missing_by_cell.setdefault(item["name"], parsed[2])
            require(build not in reported_builds or
                    reported_builds[build] == parsed[3],
                    "benchmark build identity changed: " + name)
            reported_builds.setdefault(build, parsed[3])
            raw_by_name[name] = parsed[:2]
            raw_records[name] = encode_raw_record(completed.stdout, completed.stderr)
            entries.append({
                "name": name, "cell": item["name"], "round": round_number,
                "sequence": sequence, "build": build, "argv": logical_argv,
                "command_sha256": sha256_bytes(canonical_bytes(command_record)),
                "binary_sha256": builds[build]["artifacts"]["benchmark_sha256"],
                "stdout_sha256": sha256_bytes(completed.stdout),
                "stderr_sha256": sha256_bytes(completed.stderr),
                "returncode": completed.returncode,
            })
            reservation_handle.validate_current()
        campaign_end = time.time()
        require(reported_builds.get("baseline") == reported_builds.get("candidate"),
                "baseline/candidate reported compiler identity mismatch")
        require(sorted(os.sched_getaffinity(0)) == enforced_affinity,
                "runner affinity changed during campaign")
        post_frequency = current_frequency(args.cpu)
        require(power_state(args.cpu) == static_power,
                "CPU power/governor state changed during campaign")
        reservation_handle.validate_current()
        execution["host"] = host_record(
            args.cpu, pre_frequency, post_frequency, static_power)
        validate_execution(execution)
        for build in ("baseline", "candidate"):
            require(sha256_file(binaries[build]) ==
                    builds[build]["artifacts"]["benchmark_sha256"],
                    build + " binary changed during campaign")
        bundle = {
            "schema": RAW_SCHEMA,
            "matrix_base64": base64.b64encode(matrix_raw).decode("ascii"),
            "matrix_sha256": matrix_info["sha256"],
            "raw": raw_records,
        }
        raw_bundle_path = output_directory / "abba_raw.json"
        atomic_json(raw_bundle_path, bundle)
        summary = summarize(entries, raw_by_name, args.backend)
        manifest = {
            "schema": evidence_schema,
            "status": "pending",
            "provenance": {
                "runner_sha256": sha256_file(Path(__file__).resolve()),
                "git": git_records,
                "builds": builds,
                "matrix": matrix_info,
                "execution": execution,
            },
            "campaign": {
                "order_per_round": [value[0] for value in SEQUENCES],
                "rounds": len(ROUNDS), "cell_count": len(CELLS),
                "entry_count": len(jobs), "samples_per_invocation": 7,
                "warmups_per_invocation": 3, "reuse_per_sample": 8,
                "batch": 1, "threads": 1, "seed": 42,
                "target_threshold_percent": TARGET_THRESHOLD,
                "neighbor_floor_percent": NEIGHBOR_FLOOR,
                "requested_backend": args.backend,
                "resolved_backend": args.backend,
                "self_test": self_test,
                "started_unix": campaign_start, "finished_unix": campaign_end,
            },
            "entries": entries,
            "raw_bundle_file": raw_bundle_path.name,
            "raw_bundle_sha256": sha256_file(raw_bundle_path),
            "raw_evidence_sha256": stable_raw_digest(entries),
            "summary": summary,
        }
        manifest_path = output_directory / "abba_manifest.json"
        if not all(value[metric]["accepted"] for value in summary
                   for metric in ("encode", "decode")):
            manifest["status"] = "failed"
            atomic_json(manifest_path, manifest)
            raise EvidenceError(
                "campaign performance failed its paired confidence bound; "
                "a failed manifest was retained")
        validation_path = output_directory / ".abba_manifest.validation.json"
        manifest["status"] = "passed"
        atomic_json(validation_path, manifest)
        try:
            validate_manifest(
                validation_path, repo, raw_bundle_path, binaries, args.matrix,
                allow_self_test=self_test)
        except BaseException:
            manifest["status"] = "failed"
            atomic_json(manifest_path, manifest)
            try:
                validation_path.unlink()
            except OSError:
                pass
            raise
        os.replace(str(validation_path), str(manifest_path))
        return manifest
    finally:
        try:
            if reservation_handle is not None:
                reservation_handle.close()
        finally:
            try:
                os.sched_setaffinity(0, set(initial_affinity))
            finally:
                shadow_context.cleanup()


def git_file_hashes(repo, commit, relatives):
    files = {relative: sha256_bytes(git_blob(repo, commit, relative))
             for relative in relatives}
    return {"digest": sha256_bytes(canonical_bytes(files)), "files": files}


def pre_xor_matrix_fixture(document):
    """Downgrade a synthetic current matrix to the frozen v1 contract."""
    value = copy.deepcopy(document)
    value["schema"] = MATRIX_SCHEMA_V1
    current_files = value["source_fingerprint"]["files"]
    require(set(PRE_XOR_MATRIX_SOURCE_FILES) <= set(current_files),
            "current matrix fixture omits a historical v1 source")
    files = {
        relative: current_files[relative]
        for relative in PRE_XOR_MATRIX_SOURCE_FILES
    }
    value["source_fingerprint"]["files"] = files
    value["source_fingerprint"]["digest"] = \
        sha256_bytes(canonical_bytes(files))
    fingerprint = value["source_fingerprint"]

    def rehash_command(command):
        payload = dict(command)
        payload.pop("command_sha256", None)
        command["command_sha256"] = sha256_bytes(canonical_bytes(payload))

    for variant in value["variants"]:
        variant["schema"] = MATRIX_SCHEMA_V1
        variant["source_fingerprint"] = fingerprint["digest"]
        identity = variant["build_identity"]
        retained_commands = []
        retained_counts = {}
        for command in identity["compile_commands"]:
            relative = command["file"]
            expected = PRE_XOR_MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS.get(
                relative, 0)
            used = retained_counts.get(relative, 0)
            if used < expected:
                retained_commands.append(command)
                retained_counts[relative] = used + 1
        require(retained_counts ==
                PRE_XOR_MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS,
                "current matrix fixture cannot reconstruct v1 compile closure")
        identity["compile_commands"] = retained_commands
        identity["test_executables"] = {
            name: identity["test_executables"][name]
            for name in PRE_XOR_MATRIX_RUN_TESTS
        }
        identity["compile_commands_sha256"] = sha256_bytes(
            canonical_bytes(identity["compile_commands"]))
        identity_payload = dict(identity)
        identity_payload.pop("digest", None)
        identity["digest"] = sha256_bytes(canonical_bytes(identity_payload))

        tests = variant["tests"]
        require("auto_encode_backend" in tests,
                "current matrix fixture omits current-only auto test")
        tests.pop("auto_encode_backend")
        failures = tests["backend_failures"]
        failures["ctest_executed_tests"] = sorted(
            PRE_XOR_MATRIX_BACKEND_FAILURE_TESTS)
        failures["argv"][-2] = "^leopard2_backend_failure_"
        rehash_command(failures)
        configure = variant["commands"][0]
        build = variant["commands"][1]
        target_index = build["argv"].index("--target") + 1
        build["argv"][target_index:] = list(PRE_XOR_MATRIX_BUILD_TARGETS)
        rehash_command(build)
        commands = [configure, build]
        commands.extend(tests[name] for name in PRE_XOR_MATRIX_RUN_TESTS)
        commands.extend((failures, tests["portable_isa"]))
        if variant["variant"] == "auto":
            commands.append(tests["cuda_optional"])
        for command in commands:
            rehash_command(command)
        variant["commands"] = commands

        configuration_input = {
            "c_compiler": value["c_compiler"],
            "compiler": value["compiler"],
            "generator": value["generator"],
            "jobs_per_variant": value["jobs_per_variant"],
            "machine": value["machine"],
            "source": fingerprint,
            "variant": variant["variant"],
            "environment": variant["build_environment"],
        }
        variant["configuration_id"] = sha256_bytes(
            canonical_bytes(configuration_input))
        configure_argv = variant["commands"][0]["argv"]
        variant["fresh_build"]["identity_sha256"] = sha256_bytes(
            canonical_bytes({
                "configuration_id": variant["configuration_id"],
                "configure_argv": configure_argv,
                "environment": variant["build_environment"],
            }))
    return value


def write_mock(path, factor, historical=False):
    selector_fields = ("" if historical else
                       "'force_tiled_decode':False,"
                       "'force_materialized_decode':False,")
    source = """#!/usr/bin/env python3
import json,sys
a=sys.argv[1:]
def v(name): return a[a.index(name)+1]
k=int(v('--k')); r=int(v('--r')); profile=v('--profile'); field=v('--field'); backend=v('--backend')
byte_count=int(v('--bytes')); loss=int(v('--loss')); base=float(k+r+loss+byte_count/64.0+1)
factor=%.8f
missing=list(range(loss))
def metric(value):
 return {'minimum_us_per_batch_call':value*.98,'median_us_per_batch_call':value,'maximum_us_per_batch_call':value*1.02,'mad_us_per_batch_call':value*.01}
def rated(value,input_key,output_key):
 result=metric(value); result[input_key]=1.0/value; result[output_key]=2.0/value; return result
def setup(value):
 return {'minimum_us':value*.98,'median_us':value,'maximum_us':value*1.02,'mad_us':value*.01}
resolved_profile='legacy_high_v1' if profile=='high' else 'low_v1'
legacy_comparison='matched' if profile=='high' and byte_count%%64==0 else None
active=r if profile=='high' else k
padded=1<<(active-1).bit_length()
parent_input=k+padded if profile=='high' else padded+r
parent=1<<(parent_input-1).bit_length()
encode=rated(base*factor,'input_GB_per_s','parity_output_GB_per_s')
decode=rated(base*2*factor,'offered_received_GB_per_s','repaired_output_GB_per_s')
legacy_reason=None if legacy_comparison else ('old Leopard only defines the legacy high wire profile' if profile!='high' else 'old Leopard requires shard bytes divisible by 64')
out={'schema':'leopard2-benchmark-v1','build':{'compiler':'mock','compiler_version':'1','cplusplus':201103},'parameters':{'K':k,'R':r,'requested_profile':resolved_profile,'requested_field':field,'requested_backend':backend,'force_generic_decode':False,'force_specialized_decode':False,%s'shard_bytes':byte_count,'loss_count':loss,'missing_original_indices':missing,'batch':int(v('--batch')),'reuse':int(v('--reuse')),'iterations':int(v('--iterations')),'warmup':int(v('--warmup')),'thread_count':int(v('--threads')),'seed':int(v('--seed'))},'resolved':{'profile':resolved_profile,'field':field,'backend':backend,'thread_count':1,'parent_count':parent,'padded_side':padded},'correctness':{'leopard2_round_trip':True,'legacy_comparison':legacy_comparison},'memory':{'scratch_alignment':64,'encode_scratch_bytes_per_stripe':64,'decode_scratch_bytes_per_stripe':128,'encode_scratch_bytes_batch':64,'decode_scratch_bytes_batch':128},'metrics':{'codec_setup':setup(1.0),'encode_execution':encode,'decode_plan_setup':setup(2.0),'decode_execution':decode,'decode_amortized_at_reuse':{'reuse_count':8,'derived_median_us_per_batch_call':base*2*factor+.25,'offered_received_GB_per_s':1.0/(base*2*factor+.25),'repaired_output_GB_per_s':2.0/(base*2*factor+.25)},'rate_semantics':'offered_received counts all non-null shard pointers supplied; a plan may read a deterministic subset'},'legacy':{'available':legacy_comparison is not None,'unavailable_reason':legacy_reason,'codec_setup':None,'decode_timing_includes_setup':True,'encode_execution':rated(base*1.1,'input_GB_per_s','parity_output_GB_per_s') if legacy_comparison else None,'decode_including_setup':rated(base*2.2,'offered_received_GB_per_s','repaired_output_GB_per_s') if legacy_comparison else None}}
print(json.dumps(out,sort_keys=True,allow_nan=False))
""" % (factor, selector_fields)
    path.write_text(source, encoding="utf-8")
    path.chmod(0o755)


def write_self_test_build_files(root, source_root, binary, schema=SCHEMA):
    cmake_identity = cmake_identity_for_schema(schema)
    configured_translation_units = configured_translation_units_for_schema(
        schema)
    build_translation_units = build_translation_units_for_schema(schema)
    compiler = Path(shutil.which("c++") or "/usr/bin/c++").resolve()
    cc = Path(shutil.which("cc") or "/usr/bin/cc").resolve()
    ar = Path(shutil.which("ar") or "/usr/bin/ar").resolve()
    ranlib = Path(shutil.which("ranlib") or "/usr/bin/ranlib").resolve()
    compile_commands = []
    outputs = {}
    for index, relative in enumerate(configured_translation_units):
        if relative == "bench/leopard2/benchmark.cpp":
            output = root / "CMakeFiles/bench_leopard2.dir" / \
                (relative + ".o")
        elif relative == "Leopard2BackendSSSE3.cpp":
            output = root / "CMakeFiles/leopard2_backend_ssse3.dir" / \
                (relative + ".o")
        elif relative in (
                "Leopard2BackendAVX2.cpp", "Leopard2BackendAVX2Xor.cpp"):
            output = root / "CMakeFiles/leopard2_backend_avx2.dir" / \
                (relative + ".o")
        elif relative == "Leopard2BackendAVX512.cpp":
            output = root / "CMakeFiles/leopard2_backend_avx512.dir" / \
                (relative + ".o")
        else:
            output = root / "CMakeFiles" / \
                cmake_identity["target_directory"] / (relative + ".o")
        output.parent.mkdir(parents=True, exist_ok=True)
        outputs[relative] = output
        output.write_bytes((relative + " object").encode("utf-8"))
        Path(str(output) + ".d").write_text(
            "{}: {}\n".format(output, source_root / relative),
            encoding="utf-8")
        compile_commands.append({
            "directory": str(root),
            "file": str(source_root / relative),
            "command": "{} -O2 -MMD -c {} -o {}".format(
                shlex.quote(str(compiler)),
                shlex.quote(str(source_root / relative)),
                shlex.quote(str(output))),
        })
    compile_path = root / "compile_commands.json"
    atomic_json(compile_path, compile_commands)
    cache = root / "CMakeCache.txt"
    tool_values = {
        "CMAKE_COMMAND:INTERNAL": Path(shutil.which("cmake") or "/usr/bin/cmake").resolve(),
        "CMAKE_C_COMPILER:FILEPATH": cc,
        "CMAKE_CXX_COMPILER:FILEPATH": compiler,
        "CMAKE_AR:FILEPATH": ar,
        "CMAKE_RANLIB:FILEPATH": ranlib,
        "CMAKE_CXX_COMPILER_AR:FILEPATH": Path(
            shutil.which("gcc-ar") or shutil.which("ar")).resolve(),
        "CMAKE_CXX_COMPILER_RANLIB:FILEPATH": Path(
            shutil.which("gcc-ranlib") or shutil.which("ranlib")).resolve(),
        "CMAKE_LINKER:FILEPATH": Path(shutil.which("ld")).resolve(),
        "CMAKE_MAKE_PROGRAM:FILEPATH": Path(
            shutil.which("gmake") or shutil.which("make")).resolve(),
    }
    cache_lines = [
        "{}={}".format(key, value) for key, value in sorted(tool_values.items())]
    cache_lines.extend([
        "CMAKE_BUILD_TYPE:STRING=Release",
        "CMAKE_CACHEFILE_DIR:INTERNAL={}".format(root),
        "CMAKE_C_FLAGS:STRING=",
        "CMAKE_C_FLAGS_RELEASE:STRING=-O3 -DNDEBUG",
        "CMAKE_CXX_FLAGS:STRING=",
        "CMAKE_CXX_FLAGS_RELEASE:STRING=-O3 -DNDEBUG",
        "CMAKE_EXE_LINKER_FLAGS:STRING=",
        "CMAKE_EXE_LINKER_FLAGS_RELEASE:STRING=",
        "CMAKE_STATIC_LINKER_FLAGS:STRING=",
        "CMAKE_STATIC_LINKER_FLAGS_RELEASE:STRING=",
        "CMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON",
        "CMAKE_GENERATOR:INTERNAL=Unix Makefiles",
        "CMAKE_HOME_DIRECTORY:INTERNAL={}".format(source_root),
        "ENABLE_OPENMP:BOOL=ON",
        "LEO2_BACKEND_VARIANT:STRING=auto",
        "LEO2_BUILD_BENCHMARKS:BOOL=ON",
        "LEO2_BUILD_FUZZERS:BOOL=OFF",
        "LEO2_BUILD_TESTS:BOOL=OFF",
        "LEO2_ENABLE_CUDA:BOOL=OFF",
    ])
    cache.write_text("\n".join(cache_lines) + "\n", encoding="utf-8")
    library = root / cmake_identity["archive"]
    archive_order = archive_translation_units_for_schema(schema)
    archive_objects = [outputs[value] for value in archive_order]
    command_output([str(ar), "rcs", str(library)] +
                   [str(value) for value in archive_objects], root,
                   "self-test archive")
    external_library = root.parent / "butterfly-self-test-external.a"
    if not external_library.exists():
        external_library.write_bytes(b"butterfly external link fixture\n")
    library_link = root / "CMakeFiles" / \
        cmake_identity["target_directory"] / "link.txt"
    benchmark_link = root / "CMakeFiles/bench_leopard2.dir/link.txt"
    library_link.parent.mkdir(parents=True, exist_ok=True)
    benchmark_link.parent.mkdir(parents=True, exist_ok=True)
    library_link.write_text(
        "{} qc {} {}\n{} {}\n".format(
            shlex.quote(str(ar)), shlex.quote(library.name),
            " ".join(shlex.quote(value.relative_to(root).as_posix())
                     for value in archive_objects),
            shlex.quote(str(ranlib)), shlex.quote(library.name)),
        encoding="utf-8")
    benchmark_link.write_text(
        "{} -O3 {} -o {} {} {}\n".format(
            shlex.quote(str(compiler)),
            shlex.quote(outputs[
                "bench/leopard2/benchmark.cpp"].relative_to(root).as_posix()),
            shlex.quote(binary.relative_to(root).as_posix()),
            shlex.quote(library.name), shlex.quote(str(external_library))),
        encoding="utf-8")

    core = set(build_translation_units) - {
        "Leopard2BackendSSSE3.cpp", "Leopard2BackendAVX2.cpp",
        "Leopard2BackendAVX2Xor.cpp", "Leopard2BackendAVX512.cpp",
        "bench/leopard2/benchmark.cpp"}
    def target(target_type, artifact, dependencies, compiled):
        return {"type": target_type, "artifact": "@build/" + artifact.name,
                "dependencies": dependencies,
                "sources": [{"path": "@source/" + value, "compiled": True}
                            for value in sorted(compiled)],
                "link": None}
    library_dependencies = (
        ["leopard2_backend_avx2", "leopard2_backend_avx512",
         "leopard2_backend_ssse3"] if schema == SCHEMA else
        ["leopard2_backend_avx2", "leopard2_backend_ssse3"])
    targets = {
        cmake_identity["target"]: target(
            "STATIC_LIBRARY", library, library_dependencies, core),
        "leopard2_backend_ssse3": target(
            "OBJECT_LIBRARY", outputs["Leopard2BackendSSSE3.cpp"], [],
            {"Leopard2BackendSSSE3.cpp"}),
        "leopard2_backend_avx2": target(
            "OBJECT_LIBRARY", outputs["Leopard2BackendAVX2.cpp"], [],
            ({"Leopard2BackendAVX2.cpp", "Leopard2BackendAVX2Xor.cpp"}
             if schema == SCHEMA else {"Leopard2BackendAVX2.cpp"})),
        "bench_leopard2": target(
            "EXECUTABLE", binary, [cmake_identity["target"]],
            {"bench/leopard2/benchmark.cpp"}),
    }
    if schema == SCHEMA:
        targets["leopard2_backend_avx512"] = target(
            "OBJECT_LIBRARY", outputs["Leopard2BackendAVX512.cpp"], [],
            {"Leopard2BackendAVX512.cpp"})
        for target_record in targets.values():
            target_record["artifacts"] = [target_record["artifact"]]
        targets["leopard2_backend_avx2"]["artifacts"] = [
            "@build/" + outputs["Leopard2BackendAVX2.cpp"].name,
            "@build/" + outputs["Leopard2BackendAVX2Xor.cpp"].name,
        ]
    targets["bench_leopard2"]["link"] = {
        "language": "CXX",
        "fragments": [
            {"role": "flags", "fragment": "-O3"},
            {"role": "libraries",
             "fragment": "@build/" + cmake_identity["archive"]},
            {"role": "libraries",
             "fragment": tagged_path(external_library, source_root, root)},
        ]}
    graph = {"configured_translation_units": sorted(
                 configured_translation_units),
             "targets": targets, "index_sha256": sha256_bytes(b"index"),
             "codemodel_sha256": sha256_bytes(b"codemodel")}
    graph["digest"] = sha256_bytes(canonical_bytes(graph))
    return compile_path, cache, library, graph


def coordinated_manifest_mutation(source_manifest, source_bundle, root, mutate):
    manifest = copy.deepcopy(source_manifest)
    bundle = copy.deepcopy(source_bundle)
    mutate(manifest, bundle)
    bundle_path = root / "abba_raw.json"
    atomic_json(bundle_path, bundle)
    manifest["raw_bundle_sha256"] = sha256_file(bundle_path)
    manifest_path = root / "abba_manifest.json"
    atomic_json(manifest_path, manifest)
    return manifest_path, bundle_path


def rehash_build_record(record):
    payload = dict(record)
    payload.pop("digest", None)
    record["digest"] = sha256_bytes(canonical_bytes(payload))


def rehash_compile_entry(entry):
    payload = dict(entry)
    payload.pop("command_sha256", None)
    entry["command_sha256"] = sha256_bytes(canonical_bytes(payload))


def rehash_matrix_command(command):
    payload = dict(command)
    payload.pop("command_sha256", None)
    command["command_sha256"] = sha256_bytes(canonical_bytes(payload))


def rehash_matrix_build_identity(identity):
    identity["cache_sha256"] = sha256_bytes(
        canonical_bytes(identity["cache"]))
    identity["compile_commands_sha256"] = sha256_bytes(
        canonical_bytes(identity["compile_commands"]))
    rehash_nested_record(identity)


def rehash_dependency_closure(closure):
    closure["file_count"] = len(closure["manifest"])
    closure["source_file_count"] = sum(
        item["path"].startswith("@source/") for item in closure["manifest"])
    closure["manifest_sha256"] = sha256_bytes(
        canonical_bytes(closure["manifest"]))


def rehash_nested_record(record):
    payload = dict(record)
    payload.pop("digest", None)
    record["digest"] = sha256_bytes(canonical_bytes(payload))


def set_self_test_recipe_content(build, recipe, text):
    recipes = build["link_recipes"]
    content_key = recipe + "_recipe_content"
    sha_key = recipe + "_recipe_sha256"
    recipes[content_key] = exact_utf8_text_content(
        text, "self-test {} link recipe".format(recipe))
    recipes[recipe + "_recipe_size"] = recipes[content_key]["size"]
    recipes[sha_key] = recipes[content_key]["sha256"]
    rehash_nested_record(recipes)
    rehash_build_record(build)


def mutate_self_test_recipe_argv(build, recipe, callback):
    text = build["link_recipes"][recipe + "_recipe_content"]["text"]
    lines = [shlex.split(value, posix=True) for value in text.splitlines()
             if value.strip()]
    callback(lines)
    set_self_test_recipe_content(
        build, recipe, "\n".join(shlex.join(line) for line in lines) + "\n")


def replace_identity_strings(value, replacements):
    if isinstance(value, str):
        for old, new in replacements:
            if old in ("leopard", "libleopard"):
                if value == old:
                    return new
            else:
                value = value.replace(old, new)
        return value
    if isinstance(value, list):
        return [replace_identity_strings(item, replacements) for item in value]
    if isinstance(value, dict):
        return {
            replace_identity_strings(key, replacements):
                replace_identity_strings(item, replacements)
            for key, item in value.items()
        }
    return value


def rehash_relabelled_build(build):
    for entry in build["translation_units"].values():
        rehash_compile_entry(entry)
    rehash_dependency_closure(build["dependency_closure"])
    rehash_nested_record(build["target_graph"])
    rehash_nested_record(build["link_recipes"])
    rehash_nested_record(build["archive"])
    for phase in ("configure", "build"):
        command = build["rebuild"][phase]
        command["command_sha256"] = sha256_bytes(
            canonical_bytes(command["argv"]))
    rehash_build_record(build)


def self_test_recipe_texts(build):
    recipes = build["link_recipes"]
    return {
        recipe: recipes[recipe + "_recipe_content"]["text"]
        for recipe in ("library", "benchmark")
    }


def relabel_historical_manifest_as_pre_xor(
        manifest, retained_texts_by_build):
    replacements = (
        ("liblibleopard.a", "libleopard.a"),
        ("libleopard.dir", "leopard.dir"),
        ("libleopard", "leopard"),
    )
    manifest["schema"] = PRE_XOR_SCHEMA
    for name, build in list(manifest["provenance"]["builds"].items()):
        build = replace_identity_strings(build, replacements)
        manifest["provenance"]["builds"][name] = build
        desired_objects = [
            build["translation_units"][relative]["output"]
            for relative in PRE_XOR_ARCHIVE_TRANSLATION_UNITS]
        build["link_recipes"]["library"][0][3:] = desired_objects
        members_by_object = {
            member["object"]: member for member in build["archive"]["members"]}
        require(set(members_by_object) == set(desired_objects),
                "relabeled archive fixture closure")
        build["archive"]["members"] = [
            members_by_object[object_path] for object_path in desired_objects]
        texts = retained_texts_by_build[name]
        for recipe in ("library", "benchmark"):
            content = exact_utf8_text_content(
                texts[recipe], "retained historical {} recipe".format(recipe))
            build["link_recipes"][recipe + "_recipe_content"] = content
            build["link_recipes"][recipe + "_recipe_size"] = content["size"]
            require(content["sha256"] ==
                    build["link_recipes"][recipe + "_recipe_sha256"],
                    "historical recipe fixture SHA")
        retained_library = parse_exact_recipe_content(
            build["link_recipes"]["library_recipe_content"],
            build["link_recipes"]["library_recipe_size"],
            build["link_recipes"]["library_recipe_sha256"],
            "historical library fixture")
        retained_benchmark = parse_exact_recipe_content(
            build["link_recipes"]["benchmark_recipe_content"],
            build["link_recipes"]["benchmark_recipe_size"],
            build["link_recipes"]["benchmark_recipe_sha256"],
            "historical benchmark fixture")
        normalized_benchmark = build["link_recipes"]["benchmark"][0]
        external_by_path = {
            value["path"]: value
            for value in build["link_recipes"]["external_link_inputs"]}
        for offset, normalized in enumerate(normalized_benchmark):
            if normalized.startswith("@external/"):
                external_by_path[normalized]["raw_path"] = \
                    retained_benchmark[0][offset]
        build["tools"]["ar"]["recipe_argv0"] = retained_library[0][0]
        build["tools"]["ranlib"]["recipe_argv0"] = retained_library[1][0]
        build["tools"]["cxx"]["recipe_argv0"] = retained_benchmark[0][0]
        rehash_relabelled_build(build)


def downgrade_pre_xor_manifest_for_self_test(manifest, bundle):
    replacements = (
        ("libleopard.a", "liblibleopard.a"),
        ("leopard.dir", "libleopard.dir"),
        ("leopard", "libleopard"),
    )
    manifest["schema"] = LEGACY_SCHEMA
    retained_texts = {}
    for name, build in list(manifest["provenance"]["builds"].items()):
        build = replace_identity_strings(build, replacements)
        manifest["provenance"]["builds"][name] = build
        retained_texts[name] = self_test_recipe_texts(build)
        recipes = build["link_recipes"]
        for recipe in ("library", "benchmark"):
            recipes[recipe + "_recipe_content"] = exact_utf8_text_content(
                retained_texts[name][recipe],
                "downgraded historical {} recipe".format(recipe))
            recipes[recipe + "_recipe_size"] = \
                recipes[recipe + "_recipe_content"]["size"]
            recipes[recipe + "_recipe_sha256"] = \
                recipes[recipe + "_recipe_content"]["sha256"]
            del recipes[recipe + "_recipe_content"]
            del recipes[recipe + "_recipe_size"]
        for logical_name in ("ar", "ranlib", "cxx"):
            del build["tools"][logical_name]["recipe_argv0"]
        for external in recipes["external_link_inputs"]:
            external.pop("raw_path", None)
        rehash_relabelled_build(build)
    entries = {entry["name"]: entry for entry in manifest["entries"]}
    require(set(entries) == set(bundle["raw"]),
            "downgraded historical raw geometry")
    for name, record in bundle["raw"].items():
        stdout = base64.b64decode(record["stdout_base64"], validate=True)
        raw = parse_json_bytes(stdout, "downgraded historical raw " + name)
        parameters = raw.get("parameters")
        require(isinstance(parameters, dict),
                "downgraded historical parameters " + name)
        for selector in ("force_tiled_decode", "force_materialized_decode"):
            require(type(parameters.get(selector)) is bool and
                    parameters[selector] is False,
                    "downgraded current selector " + name)
            del parameters[selector]
        encoded = canonical_bytes(raw) + b"\n"
        digest = sha256_bytes(encoded)
        record["stdout_base64"] = base64.b64encode(encoded).decode("ascii")
        record["stdout_sha256"] = digest
        entries[name]["stdout_sha256"] = digest
    manifest["raw_evidence_sha256"] = stable_raw_digest(manifest["entries"])
    recompute_self_test_summary(manifest, bundle)
    return retained_texts


def recompute_self_test_summary(manifest, bundle):
    raw_by_name = {}
    entries = {value["name"]: value for value in manifest["entries"]}
    missing_by_cell = {}
    backend = manifest["campaign"]["requested_backend"]
    for name, item, _, _, _ in expected_jobs():
        stdout, _ = decode_raw_record(bundle["raw"][name], name)
        parsed = check_raw(
            parse_json_bytes(stdout, "self-test summary " + name), item, name,
            backend,
            missing_by_cell.get(item["name"]), manifest["schema"])
        missing_by_cell.setdefault(item["name"], parsed[2])
        raw_by_name[name] = parsed[:2]
    manifest["summary"] = summarize(
        list(entries.values()), raw_by_name, backend)


def mutate_raw_stdout(manifest, bundle, name, callback):
    record = bundle["raw"][name]
    raw = parse_json_bytes(
        base64.b64decode(record["stdout_base64"], validate=True),
        "self-test raw mutation")
    callback(raw)
    encoded = canonical_bytes(raw) + b"\n"
    record["stdout_base64"] = base64.b64encode(encoded).decode("ascii")
    record["stdout_sha256"] = sha256_bytes(encoded)
    entry = next(value for value in manifest["entries"] if value["name"] == name)
    entry["stdout_sha256"] = sha256_bytes(encoded)
    manifest["raw_evidence_sha256"] = stable_raw_digest(manifest["entries"])


def mutate_embedded_matrix(manifest, bundle, callback):
    raw = base64.b64decode(bundle["matrix_base64"], validate=True)
    document = parse_json_bytes(raw, "self-test embedded matrix mutation")
    callback(document)
    encoded = canonical_bytes(document)
    digest = sha256_bytes(encoded)
    bundle["matrix_base64"] = base64.b64encode(encoded).decode("ascii")
    bundle["matrix_sha256"] = digest
    manifest["provenance"]["matrix"]["sha256"] = digest
    manifest["provenance"]["matrix"]["source_fingerprint"] = \
        document.get("source_fingerprint", {}).get("digest")
    manifest["provenance"]["matrix"]["variant_count"] = len(
        document.get("variants", []))


def expect_failure(callback, label):
    try:
        callback()
    except (EvidenceError, KeyError, OSError, TypeError, ValueError,
            json.JSONDecodeError, base64.binascii.Error):
        return
    raise EvidenceError("mutation was accepted: " + label)


def self_test(repo):
    require(hasattr(os, "sched_getaffinity"), "affinity unavailable")
    allowed = sorted(os.sched_getaffinity(0))
    cpu = None
    sibling = None
    for candidate in allowed:
        siblings = sibling_topology(candidate)["cpus"]
        partner = next((value for value in siblings
                        if value != candidate and value in allowed), None)
        if partner is not None:
            cpu, sibling = candidate, partner
            break
    require(cpu is not None, "self-test needs an allowed SMT sibling pair")
    commit = command_output(
        ["git", "rev-parse", "HEAD"], repo, "self-test HEAD").decode("ascii").strip()
    with tempfile.TemporaryDirectory(prefix="leo2-butterfly-runner-") as temporary:
        root = Path(temporary)
        # Exercise the current contract against CMake's real codemodel.  The
        # synthetic graph below remains useful for adversarial replay, but it
        # must not be able to hide a newly configured ISA object or target.
        file_api_build = root / "real-file-api"
        query = file_api_build / ".cmake/api/v1/query"
        query.mkdir(parents=True)
        (query / "codemodel-v2").write_bytes(b"")
        cmake = Path(shutil.which("cmake") or "/usr/bin/cmake").resolve()
        command_output([
            str(cmake), "-S", str(repo), "-B", str(file_api_build),
            "-G", "Unix Makefiles", "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON", "-DLEO2_BUILD_TESTS=OFF",
            "-DLEO2_BUILD_BENCHMARKS=ON", "-DLEO2_BUILD_FUZZERS=OFF",
            "-DLEO2_ENABLE_CUDA=OFF",
        ], repo, "current CMake File API integration configure")
        file_api_artifacts, file_api_graph = file_api_targets(
            repo, file_api_build, SCHEMA)
        require(set(file_api_artifacts) == set(relevant_targets_for_schema(
                    SCHEMA)) and
                file_api_graph["configured_translation_units"] ==
                    sorted(CONFIGURED_TRANSLATION_UNITS),
                "current real CMake File API integration contract")
        real_archive_recipe = file_api_build / \
            "CMakeFiles/leopard.dir/link.txt"
        recipe_lines = [
            shlex.split(line) for line in
            real_archive_recipe.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        require(len(recipe_lines) == 2 and
                [Path(value).name for value in recipe_lines[0][3:]] == [
                    Path(relative).name + ".o"
                    for relative in CURRENT_ARCHIVE_TRANSLATION_UNITS],
                "current real CMake archive member order changed")
        # The File API graph and compile_commands.json have intentionally
        # different scopes.  Exercise the complete real collection path so an
        # unrelated configured benchmark cannot either enter the linked-source
        # closure or make every current campaign fail after a successful build.
        real_build = fresh_rebuild(
            repo, file_api_build / "CMakeCache.txt",
            root / "real-build-record", 2)
        compile_document = read_json(
            real_build["compile_commands"],
            "current real compile_commands integration")
        require(isinstance(compile_document, list) and
                len(compile_document) > len(CONFIGURED_TRANSLATION_UNITS),
                "current integration must cover a global compile database")
        real_record = build_record(
            repo, git_record(repo, commit), real_build["compile_commands"],
            real_build["cmake_cache"], real_build["library"],
            real_build["binary"], real_build["rebuild"],
            real_build["target_graph"], SCHEMA)
        validate_build_record(real_record, repo, SCHEMA)
        require(real_record["configured_translation_units"] ==
                sorted(CONFIGURED_TRANSLATION_UNITS),
                "current real build record leaked unrelated configured units")
        marker = root / "escaped-descendant-marker"
        ready = root / "escaped-descendant-ready"
        escaped_program = (
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
        subreaper_before = _get_child_subreaper()
        try:
            run_benchmark_bounded(
                [sys.executable, "-c", escaped_program,
                 str(marker), str(ready)], os.environ, 1)
        except EvidenceError as error:
            require("exceeded" in str(error),
                    "escaped-descendant self-test used the wrong failure")
        else:
            raise EvidenceError("escaped-descendant timeout was accepted")
        escaped_pid = int(ready.read_text(encoding="utf-8"))
        require(_get_child_subreaper() == subreaper_before,
                "escaped-descendant cleanup did not restore subreaper state")
        time.sleep(1.0)
        require(not marker.exists(),
                "setsid/double-fork descendant survived benchmark timeout")
        require(not Path("/proc", str(escaped_pid)).exists(),
                "setsid/double-fork descendant process survived benchmark timeout")
        baseline_root = root / "baseline-build"
        candidate_root = root / "candidate-build"
        baseline_root.mkdir()
        candidate_root.mkdir()
        baseline = baseline_root / "bench_leopard2"
        candidate = candidate_root / "bench_leopard2"
        write_mock(baseline, 1.0)
        write_mock(candidate, 0.8)
        baseline_build = write_self_test_build_files(
            baseline_root, repo, baseline)
        candidate_build = write_self_test_build_files(
            candidate_root, repo, candidate)
        fingerprint = git_file_hashes(repo, commit, MATRIX_SOURCE_FILES)
        empty_digest = sha256_bytes(b"")
        cache = cache_values(candidate_build[1])

        def matrix_compiler_identity(cache_key, logical_name):
            path = cache[cache_key]
            _, record = executable_identity(path, logical_name, cache_key)
            version = record["version"].strip()
            return {
                "executable": str(Path(path).resolve()),
                "binary_sha256": record["binary_sha256"],
                "version": version,
                "version_sha256": sha256_bytes(version.encode("utf-8")),
            }

        def matrix_tool_identity(path):
            executable = Path(path).resolve()
            completed = subprocess.run(
                [str(executable), "--version"], stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, check=False)
            raw = completed.stdout + completed.stderr
            require(completed.returncode == 0 and raw,
                    "self-test matrix tool identity")
            return {"basename": executable.name,
                    "binary_sha256": sha256_file(executable),
                    "version_sha256": sha256_bytes(raw)}

        def matrix_command(label, argv, extra=None):
            record = {
                "argv": [str(value) for value in argv],
                "cwd": str(repo.resolve()),
                "label": label,
                "returncode": 0,
                "stderr_log": label + ".stderr.log",
                "stderr_sha256": empty_digest,
                "stdout_log": label + ".stdout.log",
                "stdout_sha256": empty_digest,
                "timed_out": False,
            }
            if extra:
                record.update(extra)
            payload = dict(record)
            record["command_sha256"] = sha256_bytes(canonical_bytes(payload))
            return record

        matrix_c_compiler = matrix_compiler_identity(
            "CMAKE_C_COMPILER", "cc")
        matrix_cxx_compiler = matrix_compiler_identity(
            "CMAKE_CXX_COMPILER", "cxx")
        matrix_cmake = str(Path(cache["CMAKE_COMMAND"]).resolve())
        matrix_ctest = str(Path(shutil.which("ctest") or "/usr/bin/ctest").resolve())
        matrix_taskset = str(Path(shutil.which("taskset") or
                                  "/usr/bin/taskset").resolve())
        matrix_generator = "Unix Makefiles"
        matrix_machine = {
            "allowed_cpu_list": "{},{}".format(cpu, sibling),
            "architecture": "self-test", "cpu_flags": ["avx2"],
            "logical_cpus_allowed": 2, "platform": "self-test",
        }
        variants = []
        for value in ("auto", "scalar", "ssse3", "avx2"):
            build_environment = {
                "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
                "OMP_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"}
            if value != "auto":
                build_environment["LEO2_EXPECT_BACKEND"] = value
            matrix_build_root = str((root / "matrix-build" / value).resolve())
            compile_commands = []
            for relative in sorted(MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS):
                language = "C" if relative.endswith(".c") else "CXX"
                tool = "@tool/cc" if language == "C" else "@tool/cxx"
                for index in range(MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS[relative]):
                    output_name = re.sub(r"[^A-Za-z0-9_.-]", "_", relative)
                    compile_commands.append({
                        "file": relative,
                        "language": language,
                        "argv": [tool, "-O3", "-DNDEBUG", "-c",
                                 "@source/" + relative, "-o",
                                 "@build/mock/{}-{}.o".format(
                                     output_name, index)],
                    })
            build_cache = {
                "CMAKE_BUILD_TYPE": "Release",
                "CMAKE_GENERATOR": matrix_generator,
                "CMAKE_C_FLAGS": "", "CMAKE_C_FLAGS_RELEASE": "-O3 -DNDEBUG",
                "CMAKE_CXX_FLAGS": "",
                "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
                "CMAKE_EXE_LINKER_FLAGS": "",
                "CMAKE_EXE_LINKER_FLAGS_RELEASE": "",
                "CMAKE_STATIC_LINKER_FLAGS": "",
                "CMAKE_STATIC_LINKER_FLAGS_RELEASE": "",
                "ENABLE_OPENMP": "ON", "LEO2_BACKEND_VARIANT": value,
                "LEO2_BUILD_TESTS": "ON", "LEO2_BUILD_BENCHMARKS": "OFF",
                "LEO2_BUILD_FUZZERS": "OFF", "LEO2_ENABLE_CUDA": "OFF",
                "CMAKE_C_COMPILER": "@tool/cc",
                "CMAKE_CXX_COMPILER": "@tool/cxx",
            }
            tools = {
                "cmake": matrix_tool_identity(matrix_cmake),
                "ctest": matrix_tool_identity(matrix_ctest),
                "cc": {"basename": Path(matrix_c_compiler["executable"]).name,
                       "binary_sha256": matrix_c_compiler["binary_sha256"],
                       "version_sha256": matrix_c_compiler["version_sha256"]},
                "cxx": {"basename": Path(matrix_cxx_compiler["executable"]).name,
                        "binary_sha256": matrix_cxx_compiler["binary_sha256"],
                        "version_sha256": matrix_cxx_compiler["version_sha256"]},
            }
            test_executables = {
                name: {"path": "@build/" + MATRIX_TEST_SPECS[name][0],
                       "sha256": empty_digest}
                for name in MATRIX_RUN_TESTS
            }
            build_identity = {
                "cache": build_cache,
                "compile_commands": compile_commands,
                "test_executables": test_executables,
                "tools": tools,
            }
            build_identity["cache_sha256"] = sha256_bytes(
                canonical_bytes(build_identity["cache"]))
            build_identity["compile_commands_sha256"] = sha256_bytes(
                canonical_bytes(build_identity["compile_commands"]))
            build_identity["digest"] = sha256_bytes(canonical_bytes(build_identity))
            configuration_input = {
                "c_compiler": matrix_c_compiler,
                "compiler": matrix_cxx_compiler,
                "generator": matrix_generator,
                "jobs_per_variant": 1,
                "machine": matrix_machine,
                "source": fingerprint,
                "variant": value,
                "environment": build_environment,
            }
            configuration_id = sha256_bytes(canonical_bytes(configuration_input))
            configure_argv = [
                matrix_cmake, "-S", str(repo.resolve()), "-B", matrix_build_root,
                "-G", matrix_generator, "-DCMAKE_BUILD_TYPE=Release",
                "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
                "-DCMAKE_C_COMPILER={}".format(matrix_c_compiler["executable"]),
                "-DCMAKE_CXX_COMPILER={}".format(
                    matrix_cxx_compiler["executable"]),
                "-DLEO2_BACKEND_VARIANT={}".format(value),
                "-DLEO2_BUILD_TESTS=ON", "-DLEO2_BUILD_BENCHMARKS=OFF",
                "-DLEO2_BUILD_FUZZERS=OFF", "-DLEO2_ENABLE_CUDA=OFF",
            ]
            commands = [matrix_command("configure", configure_argv)]
            commands.append(matrix_command("build", [
                matrix_cmake, "--build", matrix_build_root, "--config", "Release",
                "-j", "1", "--target"] + list(MATRIX_BUILD_TARGETS)))
            tests = {}
            for name in MATRIX_RUN_TESTS:
                target, arguments = MATRIX_TEST_SPECS[name]
                test = matrix_command(
                    "test_" + name,
                    [matrix_taskset, "-c", str(cpu),
                     str(Path(matrix_build_root) / target)] + arguments,
                    {"executable_sha256": empty_digest})
                tests[name] = test
                commands.append(test)
            failures = matrix_command(
                "test_backend_failures",
                [matrix_taskset, "-c", str(cpu), matrix_ctest,
                 "--test-dir", matrix_build_root, "-C", "Release", "-R",
                 MATRIX_BACKEND_FAILURE_CTEST_REGEX,
                 "--output-on-failure"],
                {"ctest_executed": True,
                 "ctest_executed_tests":
                     sorted(MATRIX_BASE_BACKEND_FAILURE_TESTS +
                            MATRIX_AVX512_BACKEND_FAILURE_TESTS)})
            tests["backend_failures"] = failures
            commands.append(failures)
            portable = matrix_command(
                "test_portable_isa",
                [matrix_ctest, "--test-dir", matrix_build_root, "-C", "Release",
                 "-R", MATRIX_PORTABLE_CTEST_REGEX,
                 "--output-on-failure"],
                {"ctest_executed": True})
            tests["portable_isa"] = portable
            commands.append(portable)
            if value == "auto":
                cuda = matrix_command(
                    "test_cuda_optional",
                    [matrix_ctest, "--test-dir", matrix_build_root, "-C", "Release",
                     "-R", MATRIX_CUDA_CTEST_REGEX,
                     "--output-on-failure"],
                    {"ctest_executed": True})
                tests["cuda_optional"] = cuda
                commands.append(cuda)
            variants.append({
                "configuration_id": configuration_id, "resumed": False,
                "variant": value, "status": "passed", "reason": "",
                "schema": MATRIX_SCHEMA,
                "selected_cache_variant": value,
                "expected_runtime_backend": None if value == "auto" else value,
                "pin_cpu": cpu, "commands": commands,
                "build_environment": build_environment,
                "build_identity": build_identity,
                "fresh_build": {"configured_from_empty": True,
                                "identity_sha256": sha256_bytes(canonical_bytes({
                                    "configuration_id": configuration_id,
                                    "configure_argv": configure_argv,
                                    "environment": build_environment,
                                }))},
                "source_fingerprint": fingerprint["digest"],
                "tests": tests,
            })
        matrix = root / "matrix.json"
        matrix_document = {
            "c_compiler": matrix_c_compiler,
            "compiler": matrix_cxx_compiler,
            "generator": matrix_generator,
            "jobs": 1, "jobs_per_variant": 1,
            "machine": matrix_machine,
            "schema": MATRIX_SCHEMA, "status": "passed",
            "source_changed_during_run": False, "mismatches": [],
            "source_fingerprint": fingerprint, "variant_workers": 1,
            "variants": variants,
        }
        atomic_json(matrix, matrix_document)
        pre_xor_matrix = root / "matrix-v1-pre-xor.json"
        atomic_json(pre_xor_matrix, pre_xor_matrix_fixture(matrix_document))
        reservation = root / "reservation.json"
        reservation_document = {
            "schema": RESERVATION_SCHEMA, "status": "held",
            "benchmark_cpu": cpu, "reserved_sibling": sibling,
            "owner": "runner self-test", "nonce": "self-test-0123456789abcdef",
        }
        reservation.write_bytes(canonical_bytes(reservation_document))
        output = root / "evidence"
        args = argparse.Namespace(
            baseline=baseline, candidate=candidate,
            baseline_commit=commit, candidate_commit=commit,
            baseline_source_root=repo, candidate_source_root=repo,
            baseline_compile_commands=baseline_build[0],
            candidate_compile_commands=candidate_build[0],
            baseline_cmake_cache=baseline_build[1],
            candidate_cmake_cache=candidate_build[1],
            baseline_library=baseline_build[2], candidate_library=candidate_build[2],
            matrix=matrix, output=output, cpu=cpu, reserved_sibling=sibling,
            reservation_file=reservation, build_jobs=1, timeout=30,
            backend="avx2",
        )
        args.baseline_self_test_artifacts = baseline_build[3]
        args.candidate_self_test_artifacts = candidate_build[3]
        run_campaign(args, repo, allow_dirty=True, self_test=True)

        ssse3_args = copy.copy(args)
        ssse3_args.backend = "ssse3"
        ssse3_args.output = root / "evidence-ssse3"
        run_campaign(ssse3_args, repo, allow_dirty=True, self_test=True)
        validate_manifest(
            ssse3_args.output / "abba_manifest.json", repo,
            ssse3_args.output / "abba_raw.json", None, matrix,
            allow_self_test=True)

        pre_xor_baseline_root = root / "pre-xor-baseline-build"
        pre_xor_candidate_root = root / "pre-xor-candidate-build"
        pre_xor_baseline_root.mkdir()
        pre_xor_candidate_root.mkdir()
        pre_xor_baseline = pre_xor_baseline_root / "bench_leopard2"
        pre_xor_candidate = pre_xor_candidate_root / "bench_leopard2"
        write_mock(pre_xor_baseline, 1.0)
        write_mock(pre_xor_candidate, 0.8)
        pre_xor_baseline_build = write_self_test_build_files(
            pre_xor_baseline_root, repo, pre_xor_baseline, PRE_XOR_SCHEMA)
        pre_xor_candidate_build = write_self_test_build_files(
            pre_xor_candidate_root, repo, pre_xor_candidate, PRE_XOR_SCHEMA)
        pre_xor_args = copy.copy(args)
        pre_xor_args.baseline = pre_xor_baseline
        pre_xor_args.candidate = pre_xor_candidate
        pre_xor_args.baseline_compile_commands = pre_xor_baseline_build[0]
        pre_xor_args.candidate_compile_commands = pre_xor_candidate_build[0]
        pre_xor_args.baseline_cmake_cache = pre_xor_baseline_build[1]
        pre_xor_args.candidate_cmake_cache = pre_xor_candidate_build[1]
        pre_xor_args.baseline_library = pre_xor_baseline_build[2]
        pre_xor_args.candidate_library = pre_xor_candidate_build[2]
        pre_xor_args.baseline_self_test_artifacts = pre_xor_baseline_build[3]
        pre_xor_args.candidate_self_test_artifacts = pre_xor_candidate_build[3]
        pre_xor_args.matrix = pre_xor_matrix
        pre_xor_args.output = root / "pre-xor-v7-evidence"
        run_campaign(
            pre_xor_args, repo, allow_dirty=True, self_test=True,
            evidence_schema=PRE_XOR_SCHEMA)
        pre_xor_manifest_path = \
            pre_xor_args.output / "abba_manifest.json"
        pre_xor_bundle_path = pre_xor_args.output / "abba_raw.json"
        validate_manifest(
            pre_xor_manifest_path, repo, pre_xor_bundle_path, None,
            pre_xor_matrix, allow_self_test=True)
        relabeled_pre_xor = read_json(
            pre_xor_manifest_path, "pre-XOR self-test manifest")
        relabeled_pre_xor["schema"] = SCHEMA
        relabeled_pre_xor_path = root / "pre-xor-relabeled-v8.json"
        atomic_json(relabeled_pre_xor_path, relabeled_pre_xor)
        expect_failure(
            lambda: validate_manifest(
                relabeled_pre_xor_path, repo, pre_xor_bundle_path, None,
                pre_xor_matrix, allow_self_test=True),
            "pre-XOR v7 evidence relabeled as current v8")

        legacy_baseline_root = root / "legacy-baseline-build"
        legacy_candidate_root = root / "legacy-candidate-build"
        legacy_baseline_root.mkdir()
        legacy_candidate_root.mkdir()
        legacy_baseline = legacy_baseline_root / "bench_leopard2"
        legacy_candidate = legacy_candidate_root / "bench_leopard2"
        write_mock(legacy_baseline, 1.0, historical=True)
        write_mock(legacy_candidate, 0.8, historical=True)
        legacy_baseline_build = write_self_test_build_files(
            legacy_baseline_root, repo, legacy_baseline, LEGACY_SCHEMA)
        legacy_candidate_build = write_self_test_build_files(
            legacy_candidate_root, repo, legacy_candidate, LEGACY_SCHEMA)
        legacy_args = copy.copy(args)
        legacy_args.baseline = legacy_baseline
        legacy_args.candidate = legacy_candidate
        legacy_args.baseline_compile_commands = legacy_baseline_build[0]
        legacy_args.candidate_compile_commands = legacy_candidate_build[0]
        legacy_args.baseline_cmake_cache = legacy_baseline_build[1]
        legacy_args.candidate_cmake_cache = legacy_candidate_build[1]
        legacy_args.baseline_library = legacy_baseline_build[2]
        legacy_args.candidate_library = legacy_candidate_build[2]
        legacy_args.baseline_self_test_artifacts = legacy_baseline_build[3]
        legacy_args.candidate_self_test_artifacts = legacy_candidate_build[3]
        legacy_args.matrix = pre_xor_matrix
        legacy_args.output = root / "legacy-v6-evidence"
        run_campaign(
            legacy_args, repo, allow_dirty=True, self_test=True,
            evidence_schema=LEGACY_SCHEMA)
        legacy_manifest_path = legacy_args.output / "abba_manifest.json"
        legacy_bundle_path = legacy_args.output / "abba_raw.json"
        validate_manifest(
            legacy_manifest_path, repo, legacy_bundle_path, None,
            pre_xor_matrix,
            allow_self_test=True)
        legacy_manifest = read_json(
            legacy_manifest_path, "legacy self-test manifest")
        legacy_bundle = read_json(
            legacy_bundle_path, "legacy self-test raw bundle")
        selector_name, selector_item, _, _, _ = expected_jobs()[0]
        current_bundle = read_json(
            args.output / "abba_raw.json", "current self-test raw bundle")
        current_stdout, _ = decode_raw_record(
            current_bundle["raw"][selector_name], selector_name)
        current_selector_result = parse_json_bytes(
            current_stdout, "current selector self-test")
        for selector in ("force_tiled_decode", "force_materialized_decode"):
            changed = copy.deepcopy(current_selector_result)
            changed["parameters"].pop(selector)
            expect_failure(
                lambda changed=changed: check_raw(
                    changed, selector_item, "current selector omission",
                    "avx2", evidence_schema=SCHEMA),
                "current selector omission " + selector)
            changed = copy.deepcopy(current_selector_result)
            changed["parameters"][selector] = True
            expect_failure(
                lambda changed=changed: check_raw(
                    changed, selector_item, "current selector activation",
                    "avx2", evidence_schema=SCHEMA),
                "current selector activation " + selector)
        legacy_stdout, _ = decode_raw_record(
            legacy_bundle["raw"][selector_name], selector_name)
        legacy_selector_result = parse_json_bytes(
            legacy_stdout, "historical selector self-test")
        legacy_selector_result["parameters"].update({
            "force_tiled_decode": False,
            "force_materialized_decode": False,
        })
        expect_failure(
            lambda: check_raw(
                legacy_selector_result, selector_item,
                "historical selector injection", "avx2",
                evidence_schema=LEGACY_SCHEMA),
            "historical selector injection")
        legacy_recipe_texts = {
            "baseline": {
                "library": (legacy_baseline_root / "CMakeFiles" /
                            HISTORICAL_CMAKE_IDENTITY["target_directory"] /
                            "link.txt").read_text(encoding="utf-8"),
                "benchmark": (legacy_baseline_root /
                              "CMakeFiles/bench_leopard2.dir/link.txt").read_text(
                                  encoding="utf-8"),
            },
            "candidate": {
                "library": (legacy_candidate_root / "CMakeFiles" /
                            HISTORICAL_CMAKE_IDENTITY["target_directory"] /
                            "link.txt").read_text(encoding="utf-8"),
                "benchmark": (legacy_candidate_root /
                              "CMakeFiles/bench_leopard2.dir/link.txt").read_text(
                                  encoding="utf-8"),
            },
        }
        relabeled_root = root / "legacy-relabeled-current"
        relabeled_root.mkdir()
        relabeled_manifest_path, relabeled_bundle_path = \
            coordinated_manifest_mutation(
                legacy_manifest, legacy_bundle, relabeled_root,
                lambda manifest, _bundle:
                    relabel_historical_manifest_as_pre_xor(
                        manifest, legacy_recipe_texts))
        expect_failure(
            lambda: validate_manifest(
                relabeled_manifest_path, repo, relabeled_bundle_path,
                None, pre_xor_matrix, allow_self_test=True),
            "coherently re-signed historical v6 target/archive relabel retains old recipes")

        slow_root = root / "slow-candidate-build"
        slow_root.mkdir()
        slow_candidate = slow_root / "bench_leopard2"
        write_mock(slow_candidate, 1.25)
        slow_build = write_self_test_build_files(slow_root, repo, slow_candidate)
        slow_args = copy.copy(args)
        slow_args.candidate = slow_candidate
        slow_args.candidate_compile_commands = slow_build[0]
        slow_args.candidate_cmake_cache = slow_build[1]
        slow_args.candidate_library = slow_build[2]
        slow_args.candidate_self_test_artifacts = slow_build[3]
        slow_args.output = root / "failed-performance-evidence"
        expect_failure(lambda: run_campaign(
            slow_args, repo, allow_dirty=True, self_test=True),
            "negative-performance campaign")
        failed_manifest = read_json(
            slow_args.output / "abba_manifest.json",
            "negative-performance manifest")
        require(failed_manifest.get("status") == "failed",
                "negative-performance campaign published a passed manifest")
        failed_bundle_path = slow_args.output / "abba_raw.json"
        failed_replay = validate_manifest(
            slow_args.output / "abba_manifest.json", repo,
            failed_bundle_path, None, matrix, allow_self_test=True)
        require(failed_replay.get("status") == "failed",
                "negative-performance replay lost expected policy status")

        failed_bundle = read_json(
            failed_bundle_path, "negative-performance raw bundle")

        pre_xor_slow_root = root / "pre-xor-slow-candidate-build"
        pre_xor_slow_root.mkdir()
        pre_xor_slow_candidate = \
            pre_xor_slow_root / "bench_leopard2"
        write_mock(pre_xor_slow_candidate, 1.25)
        pre_xor_slow_build = write_self_test_build_files(
            pre_xor_slow_root, repo, pre_xor_slow_candidate,
            PRE_XOR_SCHEMA)
        pre_xor_slow_args = copy.copy(pre_xor_args)
        pre_xor_slow_args.candidate = pre_xor_slow_candidate
        pre_xor_slow_args.candidate_compile_commands = \
            pre_xor_slow_build[0]
        pre_xor_slow_args.candidate_cmake_cache = pre_xor_slow_build[1]
        pre_xor_slow_args.candidate_library = pre_xor_slow_build[2]
        pre_xor_slow_args.candidate_self_test_artifacts = \
            pre_xor_slow_build[3]
        pre_xor_slow_args.output = root / "failed-performance-v7"
        expect_failure(lambda: run_campaign(
            pre_xor_slow_args, repo, allow_dirty=True, self_test=True,
            evidence_schema=PRE_XOR_SCHEMA),
            "pre-XOR negative-performance campaign")
        pre_xor_failed_manifest = read_json(
            pre_xor_slow_args.output / "abba_manifest.json",
            "pre-XOR negative-performance manifest")
        pre_xor_failed_bundle = read_json(
            pre_xor_slow_args.output / "abba_raw.json",
            "pre-XOR negative-performance raw bundle")

        historical_failed_root = root / "failed-policy-v6"
        historical_failed_root.mkdir()
        historical_failed_texts = {}

        def downgrade_failed_policy(manifest_value, _bundle_value):
            historical_failed_texts.update(
                downgrade_pre_xor_manifest_for_self_test(
                    manifest_value, _bundle_value))

        historical_failed_manifest_path, historical_failed_bundle_path = \
            coordinated_manifest_mutation(
                pre_xor_failed_manifest, pre_xor_failed_bundle,
                historical_failed_root,
                downgrade_failed_policy)
        historical_failed_replay = validate_manifest(
            historical_failed_manifest_path, repo,
            historical_failed_bundle_path, None, pre_xor_matrix,
            allow_self_test=True)
        require(historical_failed_replay.get("status") == "failed" and
                historical_failed_replay.get("schema") == LEGACY_SCHEMA,
                "historical v6 failed-policy replay")
        historical_failed_manifest = read_json(
            historical_failed_manifest_path,
            "historical failed-policy manifest")
        historical_failed_bundle = read_json(
            historical_failed_bundle_path,
            "historical failed-policy raw bundle")
        relabeled_failed_root = root / "failed-policy-v6-relabeled-v7"
        relabeled_failed_root.mkdir()
        relabeled_failed_manifest_path, relabeled_failed_bundle_path = \
            coordinated_manifest_mutation(
                historical_failed_manifest, historical_failed_bundle,
                relabeled_failed_root,
                lambda manifest_value, _bundle_value:
                    relabel_historical_manifest_as_pre_xor(
                        manifest_value, historical_failed_texts))
        expect_failure(
            lambda: validate_manifest(
                relabeled_failed_manifest_path, repo,
                relabeled_failed_bundle_path, None, pre_xor_matrix,
                allow_self_test=True),
            "coherently re-signed failed v6 relabel retains old recipes")
        failed_mutations = (
            ("failed raw", lambda m, b:
             b["raw"][m["entries"][0]["name"]].__setitem__(
                 "stdout_base64", "e30K")),
            ("failed summary", lambda m, b:
             m["summary"][0]["encode"].__setitem__(
                 "speedup_percent", 999.0)),
            ("failed provenance", lambda m, b:
             m["provenance"].__setitem__("runner_sha256", "0" * 64)),
            ("failed status", lambda m, b:
             m.__setitem__("status", "passed")),
        )
        for index, (label, mutation) in enumerate(failed_mutations):
            mutated = root / ("failed-mutation-{}".format(index))
            mutated.mkdir()
            mp, bp = coordinated_manifest_mutation(
                failed_manifest, failed_bundle, mutated, mutation)
            expect_failure(
                lambda mp=mp, bp=bp: validate_manifest(
                    mp, repo, bp, None, matrix, allow_self_test=True),
                label)
        manifest_path = output / "abba_manifest.json"
        bundle_path = output / "abba_raw.json"
        manifest = read_json(manifest_path, "self-test manifest")
        bundle = read_json(bundle_path, "self-test raw bundle")
        validate = lambda mp, bp: validate_manifest(
            mp, repo, bp, None, None, allow_self_test=True)
        validate(manifest_path, bundle_path)

        mutations = []
        mutations.append(("current v8 target/archive relabeled as v6",
                          lambda m, b: m.__setitem__("schema", LEGACY_SCHEMA)))
        mutations.append(("passed status", lambda m, b:
                          m.__setitem__("status", "failed")))
        mutations.append(("return code", lambda m, b:
                          m["entries"][0].__setitem__("returncode", 17)))
        mutations.append(("ABBA entry order", lambda m, b:
                          m["entries"].__setitem__(slice(0, 2),
                                                   list(reversed(m["entries"][:2])))))
        mutations.append(("missing raw job", lambda m, b:
                          b["raw"].pop(m["entries"][0]["name"])))
        mutations.append(("logical argv", lambda m, b:
                          m["entries"][0]["argv"].__setitem__(0, "@wrong")))
        mutations.append(("initial affinity", lambda m, b:
                          m["provenance"]["execution"].__setitem__(
                              "initial_affinity", [9999])))
        mutations.append(("topology", lambda m, b:
                          m["provenance"]["execution"]["topology"].__setitem__(
                              "raw", str(cpu))))
        mutations.append(("reservation", lambda m, b:
                          m["provenance"]["execution"]["reservation"]["document"].__setitem__(
                              "nonce", "x")))
        mutations.append(("host power metadata", lambda m, b:
                          m["provenance"]["execution"]["host"]["record"]["power"][
                              "state"].__setitem__("scaling_governor", "forged")))
        mutations.append(("build record", lambda m, b:
                          m["provenance"]["builds"]["candidate"]["artifacts"].__setitem__(
                              "benchmark_sha256", "0" * 64)))

        def mutate_compile_command(m, b):
            build = m["provenance"]["builds"]["candidate"]
            entry = build["translation_units"]["Leopard2BackendAVX2.cpp"]
            source_index = entry["argv"].index(
                "@source/Leopard2BackendAVX2.cpp")
            entry["argv"][source_index] = "@source/Leopard2BackendScalar.cpp"
            rehash_compile_entry(entry)
            rehash_build_record(build)
        mutations.append(("coordinated compile command", mutate_compile_command))

        def mutate_xor_compile_command(m, b):
            build = m["provenance"]["builds"]["candidate"]
            entry = build["translation_units"][
                "Leopard2BackendAVX2Xor.cpp"]
            source_index = entry["argv"].index(
                "@source/Leopard2BackendAVX2Xor.cpp")
            entry["argv"][source_index] = \
                "@source/Leopard2BackendAVX2.cpp"
            rehash_compile_entry(entry)
            rehash_build_record(build)
        mutations.append(("AVX2 XOR compile source substitution",
                          mutate_xor_compile_command))

        def mutate_source_closure(m, b):
            build = m["provenance"]["builds"]["candidate"]
            closure = build["dependency_closure"]
            dependency = next(value for value in closure["manifest"]
                              if value["path"] ==
                              "@source/Leopard2BackendAVX2.cpp")
            dependency["sha256"] = "0" * 64
            rehash_dependency_closure(closure)
            rehash_build_record(build)
        mutations.append(("coordinated source dependency", mutate_source_closure))

        def mutate_xor_source_closure(m, b):
            build = m["provenance"]["builds"]["candidate"]
            closure = build["dependency_closure"]
            dependency = next(value for value in closure["manifest"]
                              if value["path"] ==
                              "@source/Leopard2BackendAVX2Xor.cpp")
            dependency["sha256"] = "0" * 64
            rehash_dependency_closure(closure)
            rehash_build_record(build)
        mutations.append(("AVX2 XOR source dependency substitution",
                          mutate_xor_source_closure))

        def mutate_toolchain(m, b):
            build = m["provenance"]["builds"]["candidate"]
            tool = build["tools"]["cxx"]
            tool["version"] = "forged compiler\n"
            tool["version_sha256"] = sha256_bytes(
                tool["version"].encode("utf-8"))
            rehash_build_record(build)
        mutations.append(("coordinated toolchain", mutate_toolchain))

        def mutate_configuration(m, b):
            build = m["provenance"]["builds"]["candidate"]
            build["configuration"]["CMAKE_CXX_FLAGS_RELEASE"] = \
                "-O0 -DNDEBUG"
            rehash_build_record(build)
        mutations.append(("coordinated configuration", mutate_configuration))

        def mutate_test_hooks(m, b):
            build = m["provenance"]["builds"]["candidate"]
            build["configuration"]["LEO2_BUILD_TESTS"] = "ON"
            rehash_build_record(build)
        mutations.append(("test-hook build configuration", mutate_test_hooks))

        def mutate_target_artifact(m, b):
            build = m["provenance"]["builds"]["candidate"]
            graph = build["target_graph"]
            graph["targets"]["bench_leopard2"]["artifact"] = \
                "@build/copied-bench_leopard2"
            rehash_nested_record(graph)
            rehash_build_record(build)
        mutations.append(("target artifact path substitution", mutate_target_artifact))

        def mutate_xor_target_artifact_omission(m, b):
            build = m["provenance"]["builds"]["candidate"]
            graph = build["target_graph"]
            graph["targets"]["leopard2_backend_avx2"]["artifacts"].pop()
            rehash_nested_record(graph)
            rehash_build_record(build)
        mutations.append(("AVX2 XOR File API artifact omission",
                          mutate_xor_target_artifact_omission))

        def mutate_avx512_target_omission(m, b):
            build = m["provenance"]["builds"]["candidate"]
            graph = build["target_graph"]
            graph["targets"].pop("leopard2_backend_avx512")
            graph["targets"]["leopard"]["dependencies"].remove(
                "leopard2_backend_avx512")
            rehash_nested_record(graph)
            rehash_build_record(build)
        mutations.append(("AVX512 File API target omission",
                          mutate_avx512_target_omission))

        def mutate_extra_translation_unit(m, b):
            build = m["provenance"]["builds"]["candidate"]
            build["configured_translation_units"].append("injected.cpp")
            build["configured_translation_units"].sort()
            rehash_build_record(build)
        mutations.append(("extra configured translation unit", mutate_extra_translation_unit))

        def mutate_archive_member(m, b):
            build = m["provenance"]["builds"]["candidate"]
            archive = build["archive"]
            archive["members"].append({"name": "injected.o",
                                       "object": "@build/injected.o",
                                       "sha256": "1" * 64})
            archive["member_count"] += 1
            rehash_nested_record(archive)
            rehash_build_record(build)
        mutations.append(("extra archive member", mutate_archive_member))

        def mutate_same_count_archive_member(m, b):
            build = m["provenance"]["builds"]["candidate"]
            member = build["archive"]["members"][0]
            member.update({"name": "injected.o", "object": "@build/injected.o",
                           "sha256": "1" * 64})
            rehash_nested_record(build["archive"])
            rehash_build_record(build)
        mutations.append(("same-count archive member substitution",
                          mutate_same_count_archive_member))

        def mutate_link_recipe(m, b):
            build = m["provenance"]["builds"]["candidate"]
            recipes = build["link_recipes"]
            recipes["benchmark"][0].append("@build/injected.o")
            rehash_nested_record(recipes)
            rehash_build_record(build)
        mutations.append(("extra benchmark link input", mutate_link_recipe))

        def mutate_link_library(m, b):
            build = m["provenance"]["builds"]["candidate"]
            build["link_recipes"]["benchmark"][0].append("-lmalicious")
            rehash_nested_record(build["link_recipes"])
            rehash_build_record(build)
        mutations.append(("unbound benchmark link library", mutate_link_library))

        def mutate_archive_tool(m, b):
            build = m["provenance"]["builds"]["candidate"]
            recipes = build["link_recipes"]
            recipes["library"][0][0] = "@tool/cxx_ar"
            rehash_nested_record(recipes)
            rehash_build_record(build)
        mutations.append(("archive recipe tool substitution", mutate_archive_tool))

        def mutate_recipe_content_size(m, b):
            build = m["provenance"]["builds"]["candidate"]
            recipes = build["link_recipes"]
            recipes["library_recipe_size"] += 1
            rehash_nested_record(recipes)
            rehash_build_record(build)
        mutations.append(("retained recipe outer size binding",
                          mutate_recipe_content_size))

        def mutate_recipe_content_sha(m, b):
            build = m["provenance"]["builds"]["candidate"]
            recipes = build["link_recipes"]
            recipes["library_recipe_sha256"] = "0" * 64
            rehash_nested_record(recipes)
            rehash_build_record(build)
        mutations.append(("retained recipe outer SHA binding",
                          mutate_recipe_content_sha))

        def mutate_retained_archive_tool(m, b):
            build = m["provenance"]["builds"]["candidate"]
            mutate_self_test_recipe_argv(
                build, "library",
                lambda lines: lines[0].__setitem__(0, "/forged/ar"))
        mutations.append(("retained archive tool spelling",
                          mutate_retained_archive_tool))

        def mutate_retained_archive_output(m, b):
            build = m["provenance"]["builds"]["candidate"]
            mutate_self_test_recipe_argv(
                build, "library",
                lambda lines: lines[0].__setitem__(2, "forged.a"))
        mutations.append(("retained canonical archive output",
                          mutate_retained_archive_output))

        def mutate_retained_archive_object(m, b):
            build = m["provenance"]["builds"]["candidate"]
            mutate_self_test_recipe_argv(
                build, "library",
                lambda lines: lines[0].__setitem__(
                    3, "CMakeFiles/leopard.dir/forged.cpp.o"))
        mutations.append(("retained archive object closure/order",
                          mutate_retained_archive_object))

        def mutate_retained_ranlib_output(m, b):
            build = m["provenance"]["builds"]["candidate"]
            mutate_self_test_recipe_argv(
                build, "library",
                lambda lines: lines[1].__setitem__(1, "forged.a"))
        mutations.append(("retained ranlib archive identity",
                          mutate_retained_ranlib_output))

        def mutate_retained_benchmark_archive(m, b):
            build = m["provenance"]["builds"]["candidate"]

            def replace_archive(lines):
                index = lines[0].index("libleopard.a")
                lines[0][index] = "forged.a"

            mutate_self_test_recipe_argv(
                build, "benchmark", replace_archive)
        mutations.append(("retained benchmark archive identity",
                          mutate_retained_benchmark_archive))

        def mutate_retained_benchmark_flag(m, b):
            build = m["provenance"]["builds"]["candidate"]

            def replace_flag(lines):
                index = lines[0].index("-O3")
                lines[0][index] = "-O0"

            mutate_self_test_recipe_argv(build, "benchmark", replace_flag)
        mutations.append(("retained benchmark flag identity",
                          mutate_retained_benchmark_flag))

        def mutate_retained_benchmark_external(m, b):
            build = m["provenance"]["builds"]["candidate"]
            raw_path = build["link_recipes"][
                "external_link_inputs"][0]["raw_path"]

            def replace_external(lines):
                index = lines[0].index(raw_path)
                lines[0][index] = "/forged/" + Path(raw_path).name

            mutate_self_test_recipe_argv(
                build, "benchmark", replace_external)
        mutations.append(("retained benchmark external-path identity",
                          mutate_retained_benchmark_external))

        def mutate_rebuild_target(m, b):
            build = m["provenance"]["builds"]["candidate"]
            command = build["rebuild"]["build"]
            command["argv"][-1] = "copied_benchmark"
            command["command_sha256"] = sha256_bytes(
                canonical_bytes(command["argv"]))
            rehash_build_record(build)
        mutations.append(("fresh rebuild target substitution", mutate_rebuild_target))

        def mutate_rebuild_configure(m, b):
            build = m["provenance"]["builds"]["candidate"]
            command = build["rebuild"]["configure"]
            command["argv"] = ["@tool/cmake", "-E", "true"]
            command["command_sha256"] = sha256_bytes(
                canonical_bytes(command["argv"]))
            rehash_build_record(build)
        mutations.append(("fresh configure no-op substitution",
                          mutate_rebuild_configure))

        mutations.append(("git commit", lambda m, b:
                          m["provenance"]["git"]["baseline"].__setitem__(
                              "commit", "f" * 40)))

        def mutate_matrix_source_set(m, b):
            def callback(document):
                files = document["source_fingerprint"]["files"]
                files.pop("Leopard2BackendAVX2.cpp")
                digest = sha256_bytes(canonical_bytes(files))
                document["source_fingerprint"]["digest"] = digest
                for variant in document["variants"]:
                    variant["source_fingerprint"] = digest
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("coordinated matrix source set", mutate_matrix_source_set))

        def mutate_matrix_xor_source_set(m, b):
            def callback(document):
                files = document["source_fingerprint"]["files"]
                files.pop("Leopard2BackendAVX2Xor.cpp")
                digest = sha256_bytes(canonical_bytes(files))
                document["source_fingerprint"]["digest"] = digest
                for variant in document["variants"]:
                    variant["source_fingerprint"] = digest
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("coordinated matrix AVX2 XOR source omission",
                          mutate_matrix_xor_source_set))

        def mutate_matrix_avx512_source_set(m, b):
            def callback(document):
                files = document["source_fingerprint"]["files"]
                files.pop("Leopard2BackendAVX512.cpp")
                digest = sha256_bytes(canonical_bytes(files))
                document["source_fingerprint"]["digest"] = digest
                for variant in document["variants"]:
                    variant["source_fingerprint"] = digest
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("coordinated matrix AVX512 source omission",
                          mutate_matrix_avx512_source_set))

        def mutate_matrix_test(m, b):
            def callback(document):
                document["variants"][0]["tests"]["backend_ops"]["returncode"] = 1
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("coordinated matrix test", mutate_matrix_test))

        def mutate_matrix_run_only_omission(m, b):
            def callback(document):
                document["variants"][0]["tests"].pop("pruned_transform")
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix run-only test omission",
                          mutate_matrix_run_only_omission))

        def mutate_matrix_auto_encode_omission(m, b):
            def callback(document):
                for variant in document["variants"]:
                    variant["tests"].pop("auto_encode_backend")
                    variant["build_identity"]["test_executables"].pop(
                        "auto_encode_backend")
                    variant["commands"] = [
                        command for command in variant["commands"]
                        if command.get("label") !=
                        "test_auto_encode_backend"]
                    rehash_matrix_build_identity(
                        variant["build_identity"])
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix auto-encode test omission",
                          mutate_matrix_auto_encode_omission))

        def mutate_matrix_output_mismatch(m, b):
            def callback(document):
                scalar = next(value for value in document["variants"]
                              if value["variant"] == "scalar")
                scalar["tests"]["backend_ops"]["stdout_sha256"] = "2" * 64
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("unreported backend output mismatch",
                          mutate_matrix_output_mismatch))

        def mutate_matrix_runtime_backend(m, b):
            def callback(document):
                scalar = next(value for value in document["variants"]
                              if value["variant"] == "scalar")
                scalar["expected_runtime_backend"] = "avx2"
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("forced runtime backend identity",
                          mutate_matrix_runtime_backend))

        def mutate_cuda_not_executed(m, b):
            def callback(document):
                auto = next(value for value in document["variants"]
                            if value["variant"] == "auto")
                auto["tests"]["cuda_optional"]["ctest_executed"] = False
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("CUDA optional CTest execution", mutate_cuda_not_executed))

        def mutate_cuda_source_binding(m, b):
            def callback(document):
                files = document["source_fingerprint"]["files"]
                files.pop("tests/cmake/test_cuda_optional.cmake")
                digest = sha256_bytes(canonical_bytes(files))
                document["source_fingerprint"]["digest"] = digest
                for variant in document["variants"]:
                    variant["source_fingerprint"] = digest
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("CUDA optional source closure", mutate_cuda_source_binding))

        def mutate_matrix_environment(m, b):
            def callback(document):
                document["variants"][0]["build_environment"]["LD_PRELOAD"] = \
                    "/tmp/injected.so"
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix ambient environment injection",
                          mutate_matrix_environment))

        def mutate_matrix_build_identity(m, b):
            def callback(document):
                identity = document["variants"][0]["build_identity"]
                identity["cache"]["CMAKE_CXX_FLAGS_RELEASE"] = "-O0"
                rehash_matrix_build_identity(identity)
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix coordinated build identity",
                          mutate_matrix_build_identity))

        def mutate_matrix_c_flags(m, b):
            def callback(document):
                identity = document["variants"][0]["build_identity"]
                identity["cache"]["CMAKE_C_FLAGS_RELEASE"] = "-O0"
                rehash_matrix_build_identity(identity)
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix coordinated C flags", mutate_matrix_c_flags))

        def mutate_matrix_c_compiler(m, b):
            def callback(document):
                document["c_compiler"]["version"] = "forged C compiler"
                document["c_compiler"]["version_sha256"] = sha256_bytes(
                    document["c_compiler"]["version"].encode("utf-8"))
                for variant in document["variants"]:
                    identity = variant["build_identity"]
                    identity["tools"]["cc"]["version_sha256"] = \
                        document["c_compiler"]["version_sha256"]
                    rehash_matrix_build_identity(identity)
                    configuration_input = {
                        "c_compiler": document["c_compiler"],
                        "compiler": document["compiler"],
                        "generator": document["generator"],
                        "jobs_per_variant": document["jobs_per_variant"],
                        "machine": document["machine"],
                        "source": document["source_fingerprint"],
                        "variant": variant["variant"],
                        "environment": variant["build_environment"],
                    }
                    variant["configuration_id"] = sha256_bytes(
                        canonical_bytes(configuration_input))
                    variant["fresh_build"]["identity_sha256"] = sha256_bytes(
                        canonical_bytes({
                            "configuration_id": variant["configuration_id"],
                            "configure_argv": variant["commands"][0]["argv"],
                            "environment": variant["build_environment"],
                        }))
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix coordinated C compiler",
                          mutate_matrix_c_compiler))

        def mutate_matrix_compile_source(m, b):
            def callback(document):
                identity = document["variants"][0]["build_identity"]
                command = identity["compile_commands"][0]
                old_source = "@source/" + command["file"]
                command["file"] = "injected.cpp"
                command["language"] = "CXX"
                command["argv"] = [
                    "@source/injected.cpp" if value == old_source else value
                    for value in command["argv"]]
                rehash_matrix_build_identity(identity)
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix injected compile source",
                          mutate_matrix_compile_source))

        def compile_action_omission(relative):
            def mutation(m, b):
                def callback(document):
                    identity = document["variants"][0]["build_identity"]
                    index = next(
                        index for index, command in
                        enumerate(identity["compile_commands"])
                        if command["file"] == relative)
                    identity["compile_commands"].pop(index)
                    rehash_matrix_build_identity(identity)
                mutate_embedded_matrix(m, b, callback)
            return mutation

        for relative in (
                "Leopard2BackendAVX2Xor.cpp",
                "Leopard2BackendAVX512.cpp",
                "tests/leopard2/test_auto_encode_backend.cpp",
                "tests/leopard2/test_decode_scratch_probe.cpp",
                "tests/leopard2/test_high_pruned_legacy.cpp"):
            mutations.append((
                "matrix compile action omission " + relative,
                compile_action_omission(relative)))

        def mutate_matrix_commands(m, b):
            def callback(document):
                document["variants"][0]["commands"] = [{"forged": True}]
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix forged command list", mutate_matrix_commands))

        def mutate_matrix_configuration_id(m, b):
            def callback(document):
                variant = document["variants"][0]
                variant["configuration_id"] = "4" * 64
                variant["fresh_build"]["identity_sha256"] = sha256_bytes(
                    canonical_bytes({
                        "configuration_id": variant["configuration_id"],
                        "configure_argv": variant["commands"][0]["argv"],
                        "environment": variant["build_environment"],
                    }))
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix coordinated configuration identity",
                          mutate_matrix_configuration_id))

        def mutate_matrix_fresh_identity(m, b):
            def callback(document):
                document["variants"][0]["fresh_build"]["identity_sha256"] = \
                    "5" * 64
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix fresh-build identity",
                          mutate_matrix_fresh_identity))

        def mutate_matrix_resumed(m, b):
            def callback(document):
                document["variants"][0]["resumed"] = True
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix resumed evidence", mutate_matrix_resumed))

        def mutate_matrix_executable_crosslink(m, b):
            def callback(document):
                variant = document["variants"][0]
                command = next(value for value in variant["commands"]
                               if value.get("label") == "test_backend_ops")
                command["executable_sha256"] = "6" * 64
                rehash_matrix_command(command)
                variant["tests"]["backend_ops"] = copy.deepcopy(command)
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("matrix test executable crosslink",
                          mutate_matrix_executable_crosslink))

        first_name = manifest["entries"][0]["name"]
        def mutate_raw_requested_backend(m, b):
            mutate_raw_stdout(
                m, b, first_name,
                lambda raw: raw["parameters"].__setitem__(
                    "requested_backend", "ssse3"))
        mutations.append(("raw requested backend",
                          mutate_raw_requested_backend))

        def mutate_raw_resolved_backend(m, b):
            mutate_raw_stdout(
                m, b, first_name,
                lambda raw: raw["resolved"].__setitem__(
                    "backend", "ssse3"))
        mutations.append(("raw resolved backend",
                          mutate_raw_resolved_backend))

        mutations.append((
            "campaign requested backend",
            lambda m, b: m["campaign"].__setitem__(
                "requested_backend", "ssse3")))
        mutations.append((
            "campaign resolved backend",
            lambda m, b: m["campaign"].__setitem__(
                "resolved_backend", "ssse3")))

        def mutate_profile(m, b):
            mutate_raw_stdout(
                m, b, first_name,
                lambda raw: raw["parameters"].__setitem__(
                    "requested_profile", "low_v1"))
        mutations.append(("requested profile", mutate_profile))

        tail_name = "high-gf8-tail-r1-A1-baseline"
        def mutate_tail_legacy(m, b):
            mutate_raw_stdout(
                m, b, tail_name,
                lambda raw: raw["correctness"].__setitem__(
                    "legacy_comparison", "matched"))
        mutations.append(("compact-tail legacy semantics", mutate_tail_legacy))

        def mutate_nonfinite(m, b):
            record = b["raw"][first_name]
            text = base64.b64decode(record["stdout_base64"]).decode("utf-8")
            raw = json.loads(text)
            raw["metrics"]["encode_execution"]["median_us_per_batch_call"] = 1e309
            encoded = json.dumps(raw, sort_keys=True).encode("utf-8")
            record["stdout_base64"] = base64.b64encode(encoded).decode("ascii")
            record["stdout_sha256"] = sha256_bytes(encoded)
            m["entries"][0]["stdout_sha256"] = sha256_bytes(encoded)
            m["raw_evidence_sha256"] = stable_raw_digest(m["entries"])
        mutations.append(("nonfinite timing", mutate_nonfinite))
        mutations.append(("summary", lambda m, b:
                          m["summary"][0]["encode"].__setitem__(
                              "speedup_percent", 999.0)))
        mutations.append(("promotion threshold", lambda m, b:
                          m["campaign"].__setitem__(
                              "target_threshold_percent", -99.0)))

        def mutate_high_variance(m, b):
            prefix = "high-gf8-r"
            for entry in m["entries"]:
                if (entry["name"].startswith(prefix) and
                        entry["build"] == "candidate"):
                    # Preserve the candidate's geometric mean while making
                    # the paired round ratios alternate sharply.  This
                    # remains an effective confidence-gate fixture as the
                    # campaign's round count changes.
                    scale = 0.4 if entry["round"] % 2 else 2.5
                    def widen(raw, scale=scale):
                        for metric_name in ("encode_execution", "decode_execution"):
                            metric = raw["metrics"][metric_name]
                            for key in ("minimum_us_per_batch_call",
                                        "median_us_per_batch_call",
                                        "maximum_us_per_batch_call",
                                        "mad_us_per_batch_call"):
                                metric[key] *= scale
                    mutate_raw_stdout(m, b, entry["name"], widen)
            recompute_self_test_summary(m, b)
            require(not next(value for value in m["summary"]
                             if value["name"] == "high-gf8")["encode"]["accepted"],
                    "high-variance fixture unexpectedly passed confidence gate")
        mutations.append(("paired ABBA high-variance overlap", mutate_high_variance))

        expect_failure(
            lambda: exact_utf8_text_content("\ud800", "surrogate recipe"),
            "retained recipe strict UTF-8")
        expect_failure(
            lambda: exact_utf8_text_content(
                "x" * (MAX_LINK_RECIPE_BYTES + 1), "oversized recipe"),
            "retained recipe byte bound")
        expect_failure(
            lambda: cmake_identity_for_schema({"schema": SCHEMA}),
            "non-string butterfly schema")

        for index, (label, mutation) in enumerate(mutations):
            mutated = root / ("mutation-{}".format(index))
            mutated.mkdir()
            mp, bp = coordinated_manifest_mutation(
                manifest, bundle, mutated, mutation)
            expect_failure(lambda mp=mp, bp=bp: validate(mp, bp), label)

        substituted = copy.copy(args)
        copied_binary = candidate_root / "copied-bench_leopard2"
        shutil.copy2(str(candidate), str(copied_binary))
        substituted.candidate = copied_binary
        expect_failure(lambda: validate_declared_template_paths(
            substituted, "candidate"), "caller binary path substitution")

        injected_object = candidate_root / "injected.o"
        injected_object.write_bytes(b"injected archive object")
        injected_archive = candidate_root / "injected-libleopard.a"
        shutil.copy2(str(candidate_build[2]), str(injected_archive))
        ar_tool = Path(cache_values(candidate_build[1])["CMAKE_AR"])
        command_output([str(ar_tool), "q", str(injected_archive),
                        str(injected_object)], candidate_root,
                       "self-test archive injection")
        compile_document = read_json(candidate_build[0], "self-test compile commands")
        tool_paths, _ = cache_tools(candidate_build[1])
        compiled = {}
        for entry in compile_document:
            relative, normalized = normalized_compile_entry(
                entry, repo, candidate_root, tool_paths)
            if relative in BUILD_TRANSLATION_UNITS:
                compiled[relative] = normalized
        expected_objects = {value[1] for relative, value in compiled.items()
                            if relative != "bench/leopard2/benchmark.cpp"}
        expect_failure(lambda: archive_manifest(
            injected_archive, ar_tool, expected_objects),
            "physical extra archive member")

        outside = candidate_root / "injected.cpp"
        outside.write_text("int injected;\n", encoding="utf-8")
        injected_compile = candidate_root / "compile_commands.injected.json"
        compile_document.append({
            "directory": str(candidate_root), "file": str(outside),
            "command": "{} -c {} -o {}".format(
                shlex.quote(str(tool_paths["cxx"])), shlex.quote(str(outside)),
                shlex.quote(str(candidate_root / "injected-extra.o"))),
        })
        atomic_json(injected_compile, compile_document)
        expect_failure(lambda: build_record(
            repo, manifest["provenance"]["git"]["candidate"],
            injected_compile, candidate_build[1], candidate_build[2], candidate,
            self_test_rebuild_record(1, candidate_build[1]), candidate_build[3]),
            "physical extra compile translation unit")

        missing = root / "missing-bundle"
        missing.mkdir()
        missing_manifest = missing / "abba_manifest.json"
        atomic_json(missing_manifest, manifest)
        expect_failure(lambda: validate(missing_manifest, missing / "absent.json"),
                       "missing adjacent raw bundle")

        noncanonical = root / "noncanonical-reservation.json"
        noncanonical.write_text(json.dumps(reservation_document, indent=2) + "\n",
                                encoding="utf-8")
        expect_failure(lambda: reservation_record(noncanonical, cpu, sibling),
                       "noncanonical reservation")
        anchor_root = root / "stable-anchor-runtime"
        anchor_root.mkdir(mode=0o700)
        anchor_root.chmod(0o700)
        anchored = root / "stable-anchor-reservation.json"
        anchored_old = root / "stable-anchor-reservation.old"
        anchored.write_bytes(canonical_bytes(reservation_document))
        anchored_handle, _anchored_identity = reservation_record(
            anchored, cpu, sibling, runtime_root=anchor_root)
        try:
            anchored.rename(anchored_old)
            anchored.write_bytes(canonical_bytes(reservation_document))
            expect_failure(lambda: reservation_record(
                anchored, cpu, sibling, runtime_root=anchor_root),
                "recreated reservation inode split")
        finally:
            anchored.unlink()
            anchored_old.rename(anchored)
            anchored_handle.close()

        def load_interop_module(name, path):
            specification = importlib.util.spec_from_file_location(name, path)
            require(specification is not None and
                    specification.loader is not None,
                    "cannot load cross-runner lease oracle: " + name)
            module = importlib.util.module_from_spec(specification)
            sys.modules[name] = module
            specification.loader.exec_module(module)
            return module

        def use_context(context):
            with context:
                pass

        def expect_cross_failure(callback, expected, label):
            try:
                callback()
            except expected:
                return
            raise EvidenceError("cross-runner lease overlap was accepted: " + label)

        jerasure = load_interop_module(
            "leopard2_butterfly_jerasure_lease_test",
            repo / "tools/leopard2_jerasure_compare.py")
        main_runner = load_interop_module(
            "leopard2_butterfly_main_lease_test",
            repo / "experiments/leopard2/main_compare/run_abba.py")
        c7_runner = load_interop_module(
            "leopard2_butterfly_c7_lease_test",
            repo / "experiments/leopard2/non_power_of_two/c7/run_authoritative.py")

        cross_runtime = root / "cross-runner-runtime"
        cross_runtime.mkdir(mode=0o700)
        cross_runtime.chmod(0o700)
        cross_reservation = root / "cross-runner-reservation.json"
        cross_reservation.write_bytes(canonical_bytes(reservation_document))

        # The production acquisition owns both layers at once: legacy
        # pair-only collectors are excluded by the file/socket lease, while
        # current collectors are excluded by the non-replaceable anchor.
        cross_handle, _cross_record = reservation_record(
            cross_reservation, cpu, sibling, runtime_root=cross_runtime)
        try:
            expect_cross_failure(
                lambda: use_context(jerasure.PairLease(
                    sibling, cpu, root=cross_runtime)),
                jerasure.ComparisonError,
                "butterfly then Jerasure on the same pair")
            expect_cross_failure(
                lambda: use_context(main_runner.StableLeaseAnchor(cross_runtime)),
                main_runner.EvidenceError,
                "butterfly then exact-main stable anchor")
            expect_cross_failure(
                lambda: use_context(c7_runner.StableLeaseAnchor(cross_runtime)),
                c7_runner.EvidenceError,
                "butterfly then C7 stable anchor")

            # A pair-only legacy collector on another core does not conflict
            # with this pair.  The global anchor intentionally remains more
            # conservative for current runners.
            use_context(jerasure.PairLease(
                cpu + 2000000, sibling + 2000000, root=cross_runtime))
            cross_handle.validate_current()
        finally:
            cross_handle.close()

        def acquire_backend(path, first, second):
            handle, _record = reservation_record(
                path, first, second, runtime_root=cross_runtime)
            handle.close()

        # Reverse acquisition is equally fail closed for both protocols.
        with jerasure.PairLease(cpu, sibling, root=cross_runtime):
            expect_cross_failure(
                lambda: acquire_backend(cross_reservation, cpu, sibling),
                EvidenceError,
                "Jerasure then butterfly on the same pair")

            disjoint_payload = dict(reservation_document)
            disjoint_payload["benchmark_cpu"] = cpu + 2000000
            disjoint_payload["reserved_sibling"] = sibling + 2000000
            disjoint_reservation = root / "cross-runner-disjoint.json"
            disjoint_reservation.write_bytes(canonical_bytes(disjoint_payload))
            acquire_backend(
                disjoint_reservation, cpu + 2000000, sibling + 2000000)

        for label, anchor_type in (
                ("exact-main", main_runner.StableLeaseAnchor),
                ("C7", c7_runner.StableLeaseAnchor)):
            with anchor_type(cross_runtime):
                expect_cross_failure(
                    lambda: acquire_backend(cross_reservation, cpu, sibling),
                    EvidenceError,
                    label + " stable anchor then butterfly")
            # A failed second-stage anchor acquisition must unwind the pair
            # lease; otherwise this legacy probe would still be rejected.
            use_context(jerasure.PairLease(cpu, sibling, root=cross_runtime))

        mutation_count = len(mutations) + len(failed_mutations) + 23
    print("butterfly ABBA v8 self-test passed: canonical and historical replay + {} adversarial mutations".format(
        mutation_count))


def add_run_arguments(parser):
    parser.add_argument("--backend", choices=SUPPORTED_BACKENDS, required=True,
                        help="explicit production runtime backend to qualify")
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--baseline-commit", required=True)
    parser.add_argument("--candidate-commit", required=True)
    parser.add_argument("--baseline-source-root", type=Path, required=True)
    parser.add_argument("--candidate-source-root", type=Path, required=True)
    parser.add_argument("--baseline-compile-commands", type=Path, required=True)
    parser.add_argument("--candidate-compile-commands", type=Path, required=True)
    parser.add_argument("--baseline-cmake-cache", type=Path, required=True)
    parser.add_argument("--candidate-cmake-cache", type=Path, required=True)
    parser.add_argument("--baseline-library", type=Path, required=True)
    parser.add_argument("--candidate-library", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--reserved-sibling", type=int, required=True)
    parser.add_argument("--reservation-file", type=Path, required=True)
    parser.add_argument("--build-jobs", type=int, default=min(os.cpu_count() or 1, 128))
    parser.add_argument("--timeout", type=int, default=120)


def main():
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command")
    run_parser = subparsers.add_parser("run")
    add_run_arguments(run_parser)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--manifest", type=Path, required=True)
    verify_parser.add_argument("--raw-bundle", type=Path)
    verify_parser.add_argument("--baseline", type=Path)
    verify_parser.add_argument("--candidate", type=Path)
    verify_parser.add_argument("--matrix", type=Path)
    subparsers.add_parser("self-test")
    args = parser.parse_args()
    try:
        if args.command == "run":
            run_campaign(args, repo)
            print("butterfly ABBA v8 campaign passed: backend={} cells={} entries={}".format(
                args.backend,
                len(CELLS), len(expected_jobs())))
        elif args.command == "verify":
            supplied = None
            if args.baseline is not None or args.candidate is not None:
                require(args.baseline is not None and args.candidate is not None,
                        "supply both binaries or neither")
                supplied = {"baseline": args.baseline, "candidate": args.candidate}
            manifest = validate_manifest(
                args.manifest, repo, args.raw_bundle, supplied, args.matrix)
            if manifest["status"] == "failed":
                print(
                    "butterfly ABBA evidence replay authenticated: "
                    "campaign failed its statistical policy",
                    file=sys.stderr)
                return EXPECTED_POLICY_FAILURE_EXIT
            print("butterfly ABBA path-independent evidence replay passed: schema={}".format(
                manifest["schema"]))
        elif args.command == "self-test":
            self_test(repo)
        else:
            parser.error("a command is required")
    except (EvidenceError, KeyError, OSError, TypeError, ValueError,
            json.JSONDecodeError, base64.binascii.Error) as error:
        print("butterfly ABBA evidence failed: {}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
