#!/usr/bin/env python3
"""Build and compare deterministic Leopard2 x86 backend variants.

The runner uses isolated build trees, pins each test process to an allowed CPU,
and writes resumable, machine-readable results.  It is intentionally standard
library only so that backend verification adds no runtime dependency.
"""

from __future__ import print_function

import argparse
import ast
import collections
import concurrent.futures
import hashlib
import json
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


HISTORICAL_SCHEMA = "leopard2-backend-matrix/v1"
PRE_GFNI_SCHEMA = "leopard2-backend-matrix/v2"
PRE_T16_SCHEMA = "leopard2-backend-matrix/v3"
SCHEMA = "leopard2-backend-matrix/v4"
VARIANTS = ("auto", "scalar", "ssse3", "avx2", "avx512")
COMPARE_TESTS = (
    "field_options",
    "direct_oracle",
    "backend_ops",
    "auto_encode_backend",
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
# This test intentionally enumerates every backend qualified by the selected
# build.  Forced scalar/SSSE3 variants therefore report different aggregate
# counters even though every byte comparison passes, so run it in every
# variant without requiring identical stdout.
RUN_ONLY_TESTS = ("pruned_transform",)
RUN_TESTS = COMPARE_TESTS + RUN_ONLY_TESTS
TEST_SPECS = {
    "field_options": ("leopard2_field_options_test", []),
    "direct_oracle": ("leopard2_direct_oracle_test", []),
    "backend_ops": ("leopard2_backend_ops_test", []),
    "auto_encode_backend": ("leopard2_auto_encode_backend_test", []),
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
        "--seed", "0x4c656f7061726432", "--cases", "64", "--threads", "1"
    ]),
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
BUILD_TARGETS = (
    "leopard2_field_options_test",
    "leopard2_direct_oracle_test",
    "leopard2_backend_ops_test",
    "leopard2_backend_failures_test",
    "leopard2_auto_encode_backend_test",
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
BASE_BACKEND_FAILURE_TESTS = (
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
AVX512_BACKEND_FAILURE_TESTS = (
    "leopard2_backend_failure_avx512_ff8_allocation",
    "leopard2_backend_failure_avx512_ff16_allocation",
    "leopard2_backend_failure_avx512_kat",
    "leopard2_backend_auto_avx512_kat_fallback",
)
# GFNI has no AUTO fallback case: AUTO never selects it, so there is no
# selected-backend qualification to fall back from.
GFNI_BACKEND_FAILURE_TESTS = (
    "leopard2_backend_failure_gfni_ff8_allocation",
    "leopard2_backend_failure_gfni_ff16_allocation",
    "leopard2_backend_failure_gfni_kat",
)
BACKEND_FAILURE_TESTS = (
    BASE_BACKEND_FAILURE_TESTS + AVX512_BACKEND_FAILURE_TESTS +
    GFNI_BACKEND_FAILURE_TESTS
)
BACKEND_FAILURE_CTEST_REGEX = \
    "^leopard2_backend_(failure_|auto_avx512_kat_fallback$)"
PORTABLE_CTEST_REGEX = "^leopard2_portable_isa$"
CUDA_CTEST_REGEX = "^leopard2_cuda_optional$"
PRE_T16_BUILD_CACHE_KEYS = (
    "CMAKE_BUILD_TYPE", "CMAKE_GENERATOR",
    "CMAKE_C_FLAGS", "CMAKE_C_FLAGS_RELEASE",
    "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE",
    "CMAKE_EXE_LINKER_FLAGS", "CMAKE_EXE_LINKER_FLAGS_RELEASE",
    "CMAKE_STATIC_LINKER_FLAGS", "CMAKE_STATIC_LINKER_FLAGS_RELEASE",
    "ENABLE_OPENMP", "LEO2_BACKEND_VARIANT", "LEO2_BUILD_TESTS",
    "LEO2_BUILD_BENCHMARKS", "LEO2_BUILD_FUZZERS", "LEO2_ENABLE_CUDA",
    "LEOPARD_ENABLE_GF8", "LEOPARD_ENABLE_GF16",
)
CURRENT_PRODUCTION_OBJECT_OPTIONS = (
    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED",
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED",
    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK",
    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL",
)
BUILD_CACHE_KEYS = PRE_T16_BUILD_CACHE_KEYS + \
    CURRENT_PRODUCTION_OBJECT_OPTIONS

PRE_GFNI_EXPECTED_COMPILE_SOURCE_COUNTS = {
    # The default dual-field test configuration builds the production
    # archive, the test-hook archive, and the legacy whole-TU SIMD
    # initialization fixture.  ISA sources appear in the production and
    # test-hook object libraries when their compiler probes succeed.
    "Leopard2Backend.cpp": 3,
    "Leopard2BackendAVX2.cpp": 2,
    "Leopard2BackendAVX2Xor.cpp": 2,
    "Leopard2BackendAVX512.cpp": 2,
    "Leopard2BackendSSSE3.cpp": 2,
    "Leopard2BackendScalar.cpp": 3,
    "Leopard2CpuFeatures.cpp": 3,
    "Leopard2Plan.cpp": 3,
    "LeopardCommon.cpp": 3,
    "LeopardFF16.cpp": 3,
    "LeopardFF8.cpp": 3,
    "leopard.cpp": 3,
    "leopard2.cpp": 3,
    "tests/leopard2/direct_oracle.cpp": 14,
    "tests/leopard2/direct_repair.cpp": 1,
    "tests/leopard2/fuzz_api.cpp": 1,
    "tests/leopard2/fuzz_pruned_transform.cpp": 1,
    "tests/leopard2/fuzz_replay.cpp": 2,
    "tests/leopard2/test_active_lch.cpp": 1,
    "tests/leopard2/test_api.cpp": 1,
    "tests/leopard2/test_arbitrary_counts_acceptance.cpp": 1,
    "tests/leopard2/test_auto_encode_backend.cpp": 1,
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
    "tests/leopard2/test_decode_scratch_probe.cpp": 1,
    "tests/leopard2/test_direct_encode.cpp": 1,
    "tests/leopard2/test_direct_oracle.cpp": 1,
    "tests/leopard2/test_direct_repair.cpp": 1,
    "tests/leopard2/test_encode_concurrency.cpp": 1,
    "tests/leopard2/test_encoder_gf16_legacy_matrix.cpp": 1,
    "tests/leopard2/test_field_options.cpp": 1,
    "tests/leopard2/test_gf16_padded_odd.cpp": 1,
    "tests/leopard2/test_gf16_tails.cpp": 1,
    "tests/leopard2/test_high_pruned_legacy.cpp": 2,
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
PRE_T16_EXPECTED_COMPILE_SOURCE_COUNTS = dict(
    PRE_GFNI_EXPECTED_COMPILE_SOURCE_COUNTS)
PRE_T16_EXPECTED_COMPILE_SOURCE_COUNTS["Leopard2BackendGFNI.cpp"] = 2
EXPECTED_COMPILE_SOURCE_COUNTS = dict(
    PRE_T16_EXPECTED_COMPILE_SOURCE_COUNTS)
EXPECTED_COMPILE_SOURCE_COUNTS["Leopard2BackendAVX2T2K4.cpp"] = 1
EXPECTED_COMPILE_SOURCE_COUNTS["Leopard2BackendAVX2T16B64.cpp"] = 1
EXPECTED_COMPILE_SOURCE_COUNTS["Leopard2BackendAVX2T32B256.cpp"] = 1
EXPECTED_COMPILE_SOURCE_COUNTS["Leopard2LowP32B64AVX2.cpp"] = 1
EXPECTED_COMPILE_SOURCE_COUNTS["tests/leopard2/test_t2_k4_production.cpp"] = 1
EXPECTED_COMPILE_SOURCE_COUNTS.update({
    "tests/leopard2/direct_oracle.cpp": 33,
    "tests/leopard2/test_balanced_b64_terminal.cpp": 1,
    "tests/leopard2/test_balanced_b64_terminal_production.cpp": 2,
    "tests/leopard2/test_dense_plan_policy.cpp": 1,
    "tests/leopard2/test_high_direct_production.cpp": 1,
    "tests/leopard2/test_high_low_duality.cpp": 1,
    "tests/leopard2/test_high_t32_tiny_schedule_policy.cpp": 1,
    "tests/leopard2/test_k16r8_b256_terminal.cpp": 1,
    "tests/leopard2/test_k5r4_b64_terminal.cpp": 1,
    "tests/leopard2/test_k5r5_b64_terminal.cpp": 1,
    "tests/leopard2/test_k8r4_one_shot_direct.cpp": 1,
    "tests/leopard2/test_k8r8_b64_terminal.cpp": 1,
    "tests/leopard2/test_k8r8_b64_terminal_production.cpp": 1,
    "tests/leopard2/test_k9r5_b256_terminal.cpp": 1,
    "tests/leopard2/test_low_p32_b64_terminal.cpp": 1,
    "tests/leopard2/test_one_shot_no_loss.cpp": 2,
    "tests/leopard2/test_small_direct_exhaustive.cpp": 1,
    "tests/leopard2/test_sparse_encode_production.cpp": 1,
    "tests/leopard2/test_t2_packed_terminal.cpp": 1,
    "tests/leopard2/test_t32_final_ifft_range.cpp": 1,
    "tests/leopard2/test_t4_packed_terminal_family.cpp": 1,
    "tests/leopard2/test_t4_packed_terminal_family_production.cpp": 1,
    "tests/leopard2/test_translated_metadata.cpp": 2,
})
PRE_GFNI_SOURCE_FILES = (
    "CMakeLists.txt",
    "LeopardCommon.cpp",
    "LeopardCommon.h",
    "Leopard2Backend.cpp",
    "Leopard2Backend.h",
    "Leopard2BackendScalar.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
    "Leopard2BackendAVX2Xor.cpp",
    "Leopard2BackendAVX512.cpp",
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
    "tests/leopard2/allocation_audit_config.h",
    "tests/leopard2/test_legacy_golden.cpp",
    "tests/leopard2/test_auto_encode_backend.cpp",
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
PRE_T16_SOURCE_FILES = PRE_GFNI_SOURCE_FILES + (
    "Leopard2BackendGFNI.cpp",)
SOURCE_FILES = PRE_T16_SOURCE_FILES + (
    "Leopard2BackendAVX2T2K4.cpp",
    "Leopard2BackendAVX2T16B64.cpp",
    "Leopard2BackendAVX2T32B256.cpp",
    "Leopard2LowP32B64AVX2.cpp",
    "tests/leopard2/test_t2_k4_production.cpp",
    "tests/leopard2/test_balanced_b64_terminal.cpp",
    "tests/leopard2/test_balanced_b64_terminal_production.cpp",
    "tests/leopard2/test_decode_scratch_probe.cpp",
    "tests/leopard2/test_dense_plan_policy.cpp",
    "tests/leopard2/test_high_direct_production.cpp",
    "tests/leopard2/test_high_low_duality.cpp",
    "tests/leopard2/test_high_t32_tiny_schedule_policy.cpp",
    "tests/leopard2/test_k16r8_b256_terminal.cpp",
    "tests/leopard2/test_k5r4_b64_terminal.cpp",
    "tests/leopard2/test_k5r5_b64_terminal.cpp",
    "tests/leopard2/test_k8r4_one_shot_direct.cpp",
    "tests/leopard2/test_k8r8_b64_terminal.cpp",
    "tests/leopard2/test_k8r8_b64_terminal_production.cpp",
    "tests/leopard2/test_k9r5_b256_terminal.cpp",
    "tests/leopard2/test_low_p32_b64_terminal.cpp",
    "tests/leopard2/test_one_shot_no_loss.cpp",
    "tests/leopard2/test_small_direct_exhaustive.cpp",
    "tests/leopard2/test_sparse_encode_production.cpp",
    "tests/leopard2/test_t2_packed_terminal.cpp",
    "tests/leopard2/test_t32_final_ifft_range.cpp",
    "tests/leopard2/test_t4_packed_terminal_family.cpp",
    "tests/leopard2/test_t4_packed_terminal_family_production.cpp",
    "tests/leopard2/test_translated_metadata.cpp",
)


def evidence_contract(schema=SCHEMA):
    """Return a frozen, versioned producer contract consumed by ABBA tools.

    Keeping this as data, rather than another hand-copied list, makes a matrix
    producer change fail the consumer's schema/contract binding immediately.
    Schema v2 remains available for historical butterfly-v9 replay, and v3
    retains the GFNI-only closure used by butterfly-v10.  Current schema v4
    additionally binds every promoted conditional AVX2/GF8 object, the
    dedicated T2/K4 test source, and the selectors that control that source
    closure.
    """
    if schema not in (PRE_GFNI_SCHEMA, PRE_T16_SCHEMA, SCHEMA):
        raise ValueError("unsupported backend-matrix evidence schema")
    current = schema == SCHEMA
    pre_t16 = schema == PRE_T16_SCHEMA
    contract = {
        "schema": schema,
        "source_files": tuple(
            SOURCE_FILES if current else
            PRE_T16_SOURCE_FILES if pre_t16 else PRE_GFNI_SOURCE_FILES),
        "expected_compile_source_counts": dict(
            EXPECTED_COMPILE_SOURCE_COUNTS if current else
            PRE_T16_EXPECTED_COMPILE_SOURCE_COUNTS if pre_t16 else
            PRE_GFNI_EXPECTED_COMPILE_SOURCE_COUNTS),
        "compare_tests": tuple(COMPARE_TESTS),
        "run_only_tests": tuple(RUN_ONLY_TESTS),
        "run_tests": tuple(RUN_TESTS),
        "test_specs": {
            name: (specification[0], list(specification[1]))
            for name, specification in TEST_SPECS.items()
        },
        "build_targets": tuple(BUILD_TARGETS),
        "build_cache_keys": tuple(
            BUILD_CACHE_KEYS if current else PRE_T16_BUILD_CACHE_KEYS) + (
            "CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER"),
        "base_backend_failure_tests": tuple(BASE_BACKEND_FAILURE_TESTS),
        "avx512_backend_failure_tests": tuple(
            AVX512_BACKEND_FAILURE_TESTS),
        "backend_failure_ctest_regex": BACKEND_FAILURE_CTEST_REGEX,
        "portable_ctest_regex": PORTABLE_CTEST_REGEX,
        "cuda_ctest_regex": CUDA_CTEST_REGEX,
    }
    if current or pre_t16:
        contract["gfni_backend_failure_tests"] = tuple(
            GFNI_BACKEND_FAILURE_TESTS)
    return contract


class MatrixError(Exception):
    """An actionable matrix configuration or execution error."""


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def digest_value(value):
    return digest_bytes(canonical_bytes(value))


def normalized_output(value):
    return value.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def cmake_cache_true(cache, name):
    return cache.get(name, "").upper() in ("1", "ON", "TRUE", "YES", "Y")


def expected_compile_source_counts(cache):
    """Bind optional ISA objects to the compiler probes from this build."""
    expected = dict(EXPECTED_COMPILE_SOURCE_COUNTS)

    # GCC, Clang, and AppleClang use the paired isolation flags below.  MSVC
    # always has the x64 SSSE3 object and records its AVX2 probe in
    # LEO2_FLAG_ARCH_AVX2.  AVX-512VL is deliberately unavailable in the
    # current MSVC graph, so an unsupported compiler remains a clean AUTO
    # build instead of failing the compile-command identity check.
    msvc_x86 = "LEO2_FLAG_ARCH_AVX2" in cache
    have_ssse3 = msvc_x86 or (
        cmake_cache_true(cache, "LEO2_FLAG_MSSSE3") and
        cmake_cache_true(cache, "LEO2_FLAG_MNO_AVX")
    )
    have_avx2 = cmake_cache_true(cache, "LEO2_FLAG_ARCH_AVX2") or (
        cmake_cache_true(cache, "LEO2_FLAG_MAVX2") and
        cmake_cache_true(cache, "LEO2_FLAG_MNO_AVX512F")
    )
    have_avx512 = (
        have_avx2 and
        cmake_cache_true(cache, "LEO2_FLAG_MAVX512F") and
        cmake_cache_true(cache, "LEO2_FLAG_MAVX512BW") and
        cmake_cache_true(cache, "LEO2_FLAG_MAVX512VL") and
        cmake_cache_true(cache, "LEO2_FLAG_MPREFER_VECTOR_WIDTH_256")
    )
    # The GFNI member is compiled with -mavx2 -mgfni -mno-avx512f, so its
    # probe is the AVX2 pair plus the affine flag.
    have_gfni = have_avx2 and cmake_cache_true(cache, "LEO2_FLAG_MGFNI")
    have_gf8 = cmake_cache_true(cache, "LEOPARD_ENABLE_GF8") \
        if "LEOPARD_ENABLE_GF8" in cache else True
    have_t16 = (
        have_avx2 and have_gf8 and
        (cmake_cache_true(cache, "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED")
         if "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED" in cache else True)
    )
    have_t32 = (
        have_avx2 and have_gf8 and (
            (cmake_cache_true(
                cache, "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED")
             if "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED" in cache else
             True) or
            (cmake_cache_true(
                cache, "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK")
             if "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK" in cache else
             True)))
    have_p32 = (
        have_avx2 and have_gf8 and
        (cmake_cache_true(cache, "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL")
         if "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL" in cache else True)
    )
    if not (have_avx2 and have_gf8):
        # The ordinary-archive T2/K4 gate contributes one additional direct
        # oracle compilation only when its AVX2/GF8 executable is configured.
        expected["tests/leopard2/direct_oracle.cpp"] = 14
    for source, available in (
            ("Leopard2BackendSSSE3.cpp", have_ssse3),
            ("Leopard2BackendAVX2.cpp", have_avx2),
            ("Leopard2BackendAVX2Xor.cpp", have_avx2 and have_gf8),
            ("Leopard2BackendAVX2T2K4.cpp", have_avx2 and have_gf8),
            ("Leopard2BackendAVX2T16B64.cpp", have_t16),
            ("Leopard2BackendAVX2T32B256.cpp", have_t32),
            ("Leopard2LowP32B64AVX2.cpp", have_p32),
            ("tests/leopard2/test_t2_k4_production.cpp",
             have_avx2 and have_gf8),
            ("Leopard2BackendAVX512.cpp", have_avx512),
            ("Leopard2BackendGFNI.cpp", have_gfni)):
        if not available:
            expected.pop(source)
    return expected


def require_compile_source_counts(counts, expected):
    """Fail closed when CMake's configured translation-unit graph changes."""
    if dict(counts) != dict(expected):
        raise MatrixError("compile-command source multiset mismatch")


CTEST_RESULT_LINE = re.compile(
    br"^\s*\d+/\d+\s+Test\s+#\d+:\s+(\S+)\s+\.+(?:\s|$)"
)


def ctest_executed_test_names(stdout, stderr=b""):
    names = []
    for line in normalized_output(stdout + stderr).splitlines():
        match = CTEST_RESULT_LINE.match(line)
        if match:
            names.append(match.group(1).decode("ascii", errors="replace"))
    return tuple(names)


def named_ctest_executed(stdout, test_name, stderr=b""):
    return ctest_executed_test_names(stdout, stderr) == (test_name,)


def backend_failure_ctest_executed(
    stdout, stderr=b"", expected=BACKEND_FAILURE_TESTS
):
    executed = ctest_executed_test_names(stdout, stderr)
    return (
        len(executed) == len(expected)
        and collections.Counter(executed) == collections.Counter(expected)
    )


def backend_failure_tests_for_build(build_identity):
    sources = {
        command["file"] for command in build_identity["compile_commands"]
    }
    tests = BASE_BACKEND_FAILURE_TESTS
    # The two explicit-only members are independent compiler probes, so each
    # contributes its own failure cases only when its object was compiled.
    if "Leopard2BackendAVX512.cpp" in sources:
        tests += AVX512_BACKEND_FAILURE_TESTS
    if "Leopard2BackendGFNI.cpp" in sources:
        tests += GFNI_BACKEND_FAILURE_TESTS
    return tests


def portable_ctest_executed(stdout, stderr=b""):
    return named_ctest_executed(stdout, "leopard2_portable_isa", stderr)


def cuda_ctest_executed(stdout, stderr=b""):
    return named_ctest_executed(stdout, "leopard2_cuda_optional", stderr)


def pinned_command(argv, taskset, cpu):
    command = [str(value) for value in argv]
    if taskset:
        return [str(taskset), "-c", str(cpu)] + command
    return command


def atomic_write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True, ensure_ascii=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, str(path))
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def compact_cpu_list(cpus):
    values = sorted(set(int(cpu) for cpu in cpus))
    if not values:
        return ""
    result = []
    first = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        result.append(str(first) if first == previous else "{}-{}".format(first, previous))
        first = previous = value
    result.append(str(first) if first == previous else "{}-{}".format(first, previous))
    return ",".join(result)


def allowed_cpus():
    if hasattr(os, "sched_getaffinity"):
        try:
            cpus = sorted(os.sched_getaffinity(0))
            if cpus:
                return cpus
        except OSError:
            pass
    return list(range(os.cpu_count() or 1))


def cpu_flags_for(allowed):
    try:
        text = Path("/proc/cpuinfo").read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        return []
    wanted = set(allowed)
    fallback = []
    for block in text.split("\n\n"):
        fields = {}
        for line in block.splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                fields[key.strip().lower()] = value.strip()
        flags = fields.get("flags", fields.get("features", "")).split()
        if flags and not fallback:
            fallback = flags
        try:
            processor = int(fields.get("processor", "-1"))
        except ValueError:
            processor = -1
        if processor in wanted:
            return sorted(set(flags))
    return sorted(set(fallback))


def compiler_identity(compiler):
    resolved = shutil.which(compiler)
    if not resolved:
        raise MatrixError("C++ compiler not found: {}".format(compiler))
    completed = subprocess.run(
        [resolved, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    output = normalized_output(completed.stdout + completed.stderr).decode(
        "utf-8", errors="replace"
    ).strip()
    if completed.returncode != 0:
        raise MatrixError("cannot query compiler {}: {}".format(resolved, output))
    # Driver spelling is semantically significant even when two frontends are
    # symlinks to one binary.  In particular, Clang's clang and clang++ names
    # select the C and C++ link drivers.  Resolving clang++ to clang here made
    # every C++ executable link with the C driver and omitted libstdc++.
    executable = str(Path(resolved).absolute())
    return {
        "executable": executable,
        "binary_sha256": digest_bytes(Path(resolved).read_bytes()),
        "version": output,
        "version_sha256": digest_bytes(output.encode("utf-8")),
    }


def variant_flags(variant):
    if variant == "ssse3":
        return ["-mssse3", "-mno-avx"]
    if variant == "avx2":
        return ["-mavx2", "-mno-avx512f"]
    if variant == "avx512":
        return [
            "-mavx2", "-mavx512f", "-mavx512bw", "-mavx512vl",
            "-mprefer-vector-width=256",
        ]
    return []


def compiler_accepts(compiler, flags):
    if not flags:
        return True, ""
    completed = subprocess.run(
        [compiler, "-x", "c++", "-std=c++11", "-fsyntax-only", "-"] + flags,
        input=b"int main() { return 0; }\n",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    message = normalized_output(completed.stdout + completed.stderr).decode(
        "utf-8", errors="replace"
    ).strip()
    return completed.returncode == 0, message


def availability(variant, machine, compiler):
    if variant == "auto":
        return True, ""
    architecture = machine["architecture"].lower()
    if architecture not in ("x86_64", "amd64", "i386", "i486", "i586", "i686"):
        return False, "forced variants are x86-only on this implementation"
    if variant == "scalar":
        return True, ""
    required = {
        "ssse3": ("ssse3",),
        "avx2": ("avx2",),
        "avx512": ("avx2", "avx512f", "avx512bw", "avx512vl"),
    }[variant]
    missing = [flag for flag in required if flag not in machine["cpu_flags"]]
    if missing:
        return False, "host CPU does not advertise {}".format(
            ",".join(missing)
        )
    supported, message = compiler_accepts(compiler["executable"], variant_flags(variant))
    if not supported:
        detail = ": " + message if message else ""
        return False, "compiler rejects {}{}".format(" ".join(variant_flags(variant)), detail)
    return True, ""


def source_fingerprint(source):
    files = {}
    for relative in SOURCE_FILES:
        path = source / relative
        if not path.is_file():
            raise MatrixError("required matrix input is missing: {}".format(path))
        files[relative] = digest_bytes(path.read_bytes())
    return {"digest": digest_value(files), "files": files}


def executable_path(build, name):
    direct = build / name
    if direct.is_file():
        return direct
    release = build / "Release" / (name + (".exe" if os.name == "nt" else ""))
    if release.is_file():
        return release
    windows = build / (name + ".exe")
    if windows.is_file():
        return windows
    raise MatrixError("built executable not found: {}".format(name))


def read_cache_variant(build):
    cache = build / "CMakeCache.txt"
    try:
        lines = cache.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise MatrixError("cannot read {}: {}".format(cache, error))
    prefix = "LEO2_BACKEND_VARIANT:STRING="
    for line in lines:
        if line.startswith(prefix):
            return line[len(prefix):]
    raise MatrixError("LEO2_BACKEND_VARIANT missing from {}".format(cache))


def seal_command(record):
    payload = dict(record)
    payload.pop("command_sha256", None)
    record["command_sha256"] = digest_value(payload)
    return record


def run_command(label, argv, cwd, log_dir, timeout, environment=None, hash_output=False):
    timed_out = False
    try:
        completed = subprocess.run(
            [str(item) for item in argv],
            cwd=str(cwd),
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
        returncode = completed.returncode
        stdout = normalized_output(completed.stdout)
        stderr = normalized_output(completed.stderr)
    except subprocess.TimeoutExpired as error:
        timed_out = True
        returncode = 124
        stdout = normalized_output(error.stdout or b"")
        stderr = normalized_output(error.stderr or b"")
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_name = label + ".stdout.log"
    stderr_name = label + ".stderr.log"
    (log_dir / stdout_name).write_bytes(stdout)
    (log_dir / stderr_name).write_bytes(stderr)
    record = {
        "argv": [str(item) for item in argv],
        "cwd": str(cwd),
        "label": label,
        "returncode": returncode,
        "stderr_log": stderr_name,
        "stdout_log": stdout_name,
        "timed_out": timed_out,
    }
    record["stderr_sha256"] = digest_bytes(stderr)
    record["stdout_sha256"] = digest_bytes(stdout)
    return record


def prepare_fresh_directory(path):
    path = Path(path)
    if path.exists():
        if path.is_symlink() or not path.is_dir():
            raise MatrixError("refusing non-directory build path: {}".format(path))
        shutil.rmtree(str(path))
    path.mkdir(parents=True, exist_ok=False)
    return True


def parse_cache(path):
    values = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("//") or line.startswith("#") or "=" not in line:
            continue
        key_type, value = line.split("=", 1)
        key = key_type.split(":", 1)[0]
        values[key] = value
    return values


def normalized_build_identity(build, source, c_compiler, compiler, cmake, ctest):
    cache_path = build / "CMakeCache.txt"
    compile_path = build / "compile_commands.json"
    if not cache_path.is_file() or not compile_path.is_file():
        raise MatrixError("fresh build omitted cache or compile commands")
    cache = parse_cache(cache_path)
    normalized_cache = {}
    for key in BUILD_CACHE_KEYS:
        if key not in cache:
            raise MatrixError("fresh cache omits {}".format(key))
        normalized_cache[key] = cache[key]
    compiler_cache = {
        "CMAKE_C_COMPILER": (c_compiler, "@tool/cc"),
        "CMAKE_CXX_COMPILER": (compiler, "@tool/cxx"),
    }
    for key, (identity_record, tag) in compiler_cache.items():
        cached = Path(cache.get(key, "")).resolve()
        if cached != Path(identity_record["executable"]).resolve():
            raise MatrixError("fresh cache/compiler mismatch: {}".format(key))
        normalized_cache[key] = tag
    commands = json.loads(compile_path.read_text(encoding="utf-8"))
    if not isinstance(commands, list) or not commands:
        raise MatrixError("compile_commands.json is empty or malformed")
    normalized_commands = []
    seen = set()
    replacements = ((str(source.resolve()), "@source"),
                    (str(build.resolve()), "@build"),
                    (str(Path(c_compiler["executable"]).resolve()), "@tool/cc"),
                    (str(Path(compiler["executable"]).resolve()), "@tool/cxx"))
    for entry in commands:
        path = Path(entry.get("file", ""))
        if not path.is_absolute():
            path = Path(entry.get("directory", build)) / path
        path = path.resolve()
        try:
            relative = path.relative_to(source.resolve()).as_posix()
        except ValueError:
            raise MatrixError("compile source escapes source root: {}".format(path))
        raw = entry.get("arguments")
        if raw is None:
            raw = shlex.split(entry.get("command", ""))
        argv = [str(value) for value in raw]
        if not argv:
            raise MatrixError("empty compile command: {}".format(relative))
        expected_compiler = c_compiler if relative.endswith(".c") else compiler
        executable = Path(argv[0])
        if not executable.is_absolute():
            resolved = shutil.which(argv[0], path="/usr/bin:/bin")
            if not resolved:
                raise MatrixError("compile executable not found: {}".format(argv[0]))
            executable = Path(resolved)
        if executable.resolve() != Path(expected_compiler["executable"]).resolve():
            raise MatrixError("compile language/tool mismatch: {}".format(relative))
        normalized_argv = []
        for argument in argv:
            for actual, logical in replacements:
                argument = argument.replace(actual, logical)
            normalized_argv.append(argument)
        normalized_argv[0] = "@tool/cc" if relative.endswith(".c") else "@tool/cxx"
        source_tag = "@source/" + relative
        if normalized_argv.count(source_tag) != 1:
            raise MatrixError("compile command source mismatch: {}".format(relative))
        language = "C" if relative.endswith(".c") else "CXX"
        identity = (relative, tuple(normalized_argv))
        if identity in seen:
            raise MatrixError("duplicate compile command: {}".format(relative))
        seen.add(identity)
        normalized_commands.append({"file": relative, "language": language,
                                    "argv": normalized_argv})
    normalized_commands.sort(key=lambda value: (value["file"], value["argv"]))
    counts = collections.Counter(value["file"] for value in normalized_commands)
    require_compile_source_counts(counts, expected_compile_source_counts(cache))

    def identity(path, name):
        resolved = Path(path).resolve()
        completed = subprocess.run(
            [str(resolved), "--version"], stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        raw = normalized_output(completed.stdout + completed.stderr)
        if completed.returncode != 0 or not raw:
            raise MatrixError("cannot identify {}".format(name))
        return {"basename": resolved.name,
                "binary_sha256": digest_bytes(resolved.read_bytes()),
                "version_sha256": digest_bytes(raw)}

    result = {
        "cache": normalized_cache,
        "cache_sha256": digest_value(normalized_cache),
        "compile_commands": normalized_commands,
        "compile_commands_sha256": digest_value(normalized_commands),
        "test_executables": {},
        "tools": {
            "cmake": identity(cmake, "cmake"),
            "ctest": identity(ctest, "ctest"),
            "cxx": {"basename": Path(compiler["executable"]).name,
                    "binary_sha256": compiler["binary_sha256"],
                    "version_sha256": compiler["version_sha256"]},
            "cc": {"basename": Path(c_compiler["executable"]).name,
                   "binary_sha256": c_compiler["binary_sha256"],
                   "version_sha256": c_compiler["version_sha256"]},
        },
    }
    result["digest"] = digest_value(result)
    return result


def run_variant(context, variant, index):
    result_dir = context["result_dir"] / variant
    result_path = context["result_dir"] / (variant + ".json")
    build = context["build_root"] / variant
    available, reason = availability(variant, context["machine"], context["compiler"])
    environment = {
        "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
        "OMP_NUM_THREADS": "1", "PATH": "/usr/bin:/bin",
    }
    if variant != "auto":
        environment["LEO2_EXPECT_BACKEND"] = variant
    identity_input = {
        "c_compiler": context["c_compiler"],
        "compiler": context["compiler"],
        "generator": context["generator"],
        "jobs_per_variant": context["jobs_per_variant"],
        "machine": context["machine"],
        "source": context["source_fingerprint"],
        "variant": variant,
        "environment": environment,
    }
    configuration_id = digest_value(identity_input)

    if context["resume"] and result_path.is_file():
        try:
            previous = json.loads(result_path.read_text(encoding="utf-8"))
            if (previous.get("configuration_id") == configuration_id and
                    previous.get("status") in ("passed", "unavailable")):
                return previous
        except (OSError, UnicodeError, ValueError):
            pass

    base = {
        "configuration_id": configuration_id,
        "resumed": False,
        "schema": SCHEMA,
        "source_fingerprint": context["source_fingerprint"]["digest"],
        "variant": variant,
    }
    if not available:
        base.update({"commands": [], "reason": reason, "status": "unavailable", "tests": {}})
        atomic_write_json(result_path, base)
        return base

    # A source fingerprint is not enough if an older build tree can retain an
    # object whose mtime is newer than the checkout.  Every non-resumed run is
    # therefore configured from an actually empty per-variant directory.
    prepare_fresh_directory(build)
    if result_dir.exists():
        shutil.rmtree(str(result_dir))
    result_dir.mkdir(parents=True, exist_ok=False)
    commands = []
    base["build_environment"] = dict(environment)
    configure = [
        context["cmake"], "-S", context["source"], "-B", build,
        "-G", context["generator"],
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_C_COMPILER={}".format(context["c_compiler"]["executable"]),
        "-DCMAKE_CXX_COMPILER={}".format(context["compiler"]["executable"]),
        "-DLEO2_BACKEND_VARIANT={}".format(variant),
        "-DLEO2_BUILD_TESTS=ON",
        "-DLEO2_BUILD_BENCHMARKS=OFF",
        "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_ENABLE_CUDA=OFF",
        "-DLEOPARD_ENABLE_GF8=ON",
        "-DLEOPARD_ENABLE_GF16=ON",
    ]
    base["fresh_build"] = {
        "configured_from_empty": True,
        "identity_sha256": digest_value({
            "configuration_id": configuration_id,
            "configure_argv": [str(value) for value in configure],
            "environment": environment,
        }),
    }
    command = run_command(
        "configure", configure, context["source"], result_dir,
        context["timeout"], environment
    )
    seal_command(command)
    commands.append(command)
    if command["returncode"] != 0:
        base.update({"commands": commands, "reason": "configure failed", "status": "failed", "tests": {}})
        atomic_write_json(result_path, base)
        return base

    base["build_identity"] = normalized_build_identity(
        build, context["source"], context["c_compiler"], context["compiler"],
        context["cmake"], context["ctest"])

    selected = read_cache_variant(build)
    if selected != variant:
        base.update({
            "commands": commands,
            "reason": "CMake selected {!r}, expected {!r}".format(selected, variant),
            "status": "failed",
            "tests": {},
        })
        atomic_write_json(result_path, base)
        return base

    build_command = [
        context["cmake"], "--build", build, "--config", "Release",
        "-j", str(context["jobs_per_variant"]), "--target",
    ] + list(BUILD_TARGETS)
    command = run_command(
        "build", build_command, context["source"], result_dir,
        context["timeout"], environment
    )
    seal_command(command)
    commands.append(command)
    if command["returncode"] != 0:
        base.update({"commands": commands, "reason": "build failed", "status": "failed", "tests": {}})
        atomic_write_json(result_path, base)
        return base

    tests = {}
    pin_cpu = context["allowed_cpus"][index % len(context["allowed_cpus"])]
    for name in RUN_TESTS:
        target, arguments = TEST_SPECS[name]
        executable = executable_path(build, target)
        argv = pinned_command(
            [str(executable)] + arguments, context["taskset"], pin_cpu
        )
        test_environment = environment
        if name == "context_backends":
            test_environment = dict(environment)
            test_environment["OMP_DYNAMIC"] = "FALSE"
            test_environment["OMP_NUM_THREADS"] = "4"
        command = run_command(
            "test_" + name, argv, context["source"], result_dir,
            context["timeout"], test_environment, hash_output=True
        )
        command["executable_sha256"] = digest_bytes(executable.read_bytes())
        base["build_identity"]["test_executables"][name] = {
            "path": "@build/" + executable.relative_to(build).as_posix(),
            "sha256": command["executable_sha256"],
        }
        seal_command(command)
        tests[name] = command
        commands.append(command)
        if command["returncode"] != 0:
            base.update({
                "commands": commands,
                "pin_cpu": pin_cpu,
                "reason": "{} test failed".format(name),
                "selected_cache_variant": selected,
                "status": "failed",
                "tests": tests,
            })
            atomic_write_json(result_path, base)
            return base

    failure_command = [
        context["ctest"], "--test-dir", build, "-C", "Release",
        "-R", BACKEND_FAILURE_CTEST_REGEX,
        "--output-on-failure",
    ]
    failure_command = pinned_command(
        failure_command, context["taskset"], pin_cpu
    )
    command = run_command(
        "test_backend_failures", failure_command, context["source"],
        result_dir, context["timeout"], environment, hash_output=True
    )
    tests["backend_failures"] = command
    commands.append(command)
    failure_stdout = (result_dir / command["stdout_log"]).read_bytes()
    failure_stderr = (result_dir / command["stderr_log"]).read_bytes()
    executed_failure_tests = ctest_executed_test_names(
        failure_stdout, failure_stderr
    )
    expected_failure_tests = backend_failure_tests_for_build(
        base["build_identity"]
    )
    failure_was_run = backend_failure_ctest_executed(
        failure_stdout, failure_stderr, expected_failure_tests
    )
    command["ctest_executed"] = failure_was_run
    command["ctest_executed_tests"] = sorted(executed_failure_tests)
    command["ctest_expected_tests"] = sorted(expected_failure_tests)
    seal_command(command)
    if command["returncode"] != 0 or not failure_was_run:
        base.update({
            "commands": commands,
            "pin_cpu": pin_cpu,
            "reason": (
                "backend failure matrix failed" if command["returncode"] != 0
                else "backend failure CTest set was incomplete or unexpected"
            ),
            "selected_cache_variant": selected,
            "status": "failed",
            "tests": tests,
        })
        atomic_write_json(result_path, base)
        return base

    if variant in VARIANTS:
        portable_command = [
            context["ctest"], "--test-dir", build, "-C", "Release",
            "-R", PORTABLE_CTEST_REGEX, "--output-on-failure",
        ]
        command = run_command(
            "test_portable_isa", portable_command, context["source"], result_dir,
            context["timeout"], environment, hash_output=True
        )
        tests["portable_isa"] = command
        commands.append(command)
        portable_was_run = portable_ctest_executed(
            (result_dir / command["stdout_log"]).read_bytes(),
            (result_dir / command["stderr_log"]).read_bytes(),
        )
        command["ctest_executed"] = portable_was_run
        seal_command(command)
        if command["returncode"] != 0 or not portable_was_run:
            if command["returncode"] == 0:
                reason = (
                    "portable-ISA test was not registered or executed; install "
                    "objdump or llvm-objdump and a POSIX sh, then reconfigure"
                )
            else:
                reason = "portable-ISA test failed"
            base.update({
                "commands": commands,
                "pin_cpu": pin_cpu,
                "reason": reason,
                "selected_cache_variant": selected,
                "status": "failed",
                "tests": tests,
            })
            atomic_write_json(result_path, base)
            return base

    if variant == "auto":
        cuda_command = [
            context["ctest"], "--test-dir", build, "-C", "Release",
            "-R", CUDA_CTEST_REGEX, "--output-on-failure",
        ]
        command = run_command(
            "test_cuda_optional", cuda_command, context["source"], result_dir,
            context["timeout"], environment, hash_output=True
        )
        tests["cuda_optional"] = command
        commands.append(command)
        cuda_was_run = cuda_ctest_executed(
            (result_dir / command["stdout_log"]).read_bytes(),
            (result_dir / command["stderr_log"]).read_bytes(),
        )
        command["ctest_executed"] = cuda_was_run
        seal_command(command)
        if command["returncode"] != 0 or not cuda_was_run:
            base.update({
                "commands": commands,
                "pin_cpu": pin_cpu,
                "reason": (
                    "optional-CUDA test failed" if command["returncode"] != 0
                    else "optional-CUDA test was not registered or executed"
                ),
                "selected_cache_variant": selected,
                "status": "failed",
                "tests": tests,
            })
            atomic_write_json(result_path, base)
            return base

    identity_payload = dict(base["build_identity"])
    identity_payload.pop("digest", None)
    base["build_identity"]["digest"] = digest_value(identity_payload)
    base.update({
        "commands": commands,
        "expected_runtime_backend": None if variant == "auto" else variant,
        "pin_cpu": pin_cpu,
        "reason": "",
        "selected_cache_variant": selected,
        "status": "passed",
        "tests": tests,
    })
    atomic_write_json(result_path, base)
    return base


def compare_results(results):
    passed = {result["variant"]: result for result in results if result["status"] == "passed"}
    mismatches = []
    if "auto" not in passed:
        return [{"reason": "auto variant did not pass"}]
    for variant in sorted(passed):
        if variant == "auto":
            continue
        for test in COMPARE_TESTS:
            for stream in ("stdout", "stderr"):
                key = stream + "_sha256"
                expected = passed["auto"]["tests"][test][key]
                actual = passed[variant]["tests"][test][key]
                if actual != expected:
                    mismatches.append({
                        "actual": actual,
                        "expected": expected,
                        "stream": stream,
                        "test": test,
                        "variant": variant,
                    })
    return mismatches


def matrix_run(arguments):
    source = Path(arguments.source).resolve()
    build_root = Path(arguments.build_root)
    if not build_root.is_absolute():
        build_root = source / build_root
    result_dir = Path(arguments.result_dir)
    if not result_dir.is_absolute():
        result_dir = source / result_dir
    requested = []
    for value in arguments.variants.split(","):
        variant = value.strip().lower()
        if variant and variant not in requested:
            requested.append(variant)
    invalid = [variant for variant in requested if variant not in VARIANTS]
    if invalid or not requested:
        raise MatrixError("variants must be a comma-separated subset of {}".format(",".join(VARIANTS)))

    cpus = allowed_cpus()
    jobs = min(max(1, arguments.jobs), len(cpus), 128)
    workers = min(max(1, arguments.variant_workers), len(requested), jobs)
    c_compiler = compiler_identity(arguments.c_compiler)
    compiler = compiler_identity(arguments.compiler)
    machine = {
        "allowed_cpu_list": compact_cpu_list(cpus),
        "architecture": platform.machine(),
        "cpu_flags": cpu_flags_for(cpus),
        "logical_cpus_allowed": len(cpus),
        "platform": platform.platform(),
    }
    source_state = source_fingerprint(source)
    cmake = shutil.which(arguments.cmake)
    ctest = shutil.which(arguments.ctest)
    if not cmake or not ctest:
        raise MatrixError("cmake and ctest must both be available")
    context = {
        "allowed_cpus": cpus,
        "build_root": build_root,
        "cmake": str(Path(cmake).resolve()),
        "c_compiler": c_compiler,
        "compiler": compiler,
        "ctest": str(Path(ctest).resolve()),
        "generator": arguments.generator,
        "jobs_per_variant": max(1, jobs // workers),
        "machine": machine,
        "result_dir": result_dir,
        "resume": not arguments.no_resume,
        "source": source,
        "source_fingerprint": source_state,
        "taskset": shutil.which("taskset") if os.name == "posix" else None,
        "timeout": arguments.timeout,
    }
    print(
        "backend matrix: variants={} workers={} jobs/variant={} cpus={}".format(
            ",".join(requested), workers, context["jobs_per_variant"],
            machine["allowed_cpu_list"]
        ),
        flush=True,
    )
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(run_variant, context, variant, index): variant
            for index, variant in enumerate(requested)
        }
        for future in concurrent.futures.as_completed(futures):
            variant = futures[future]
            result = future.result()
            results.append(result)
            print("{}: {}".format(variant, result["status"]), flush=True)
    results.sort(key=lambda item: VARIANTS.index(item["variant"]))

    final_state = source_fingerprint(source)
    source_changed = final_state["digest"] != source_state["digest"]
    mismatches = compare_results(results)
    failed = [result["variant"] for result in results if result["status"] == "failed"]
    summary = {
        "c_compiler": c_compiler,
        "compiler": compiler,
        "generator": arguments.generator,
        "jobs": jobs,
        "jobs_per_variant": context["jobs_per_variant"],
        "machine": machine,
        "mismatches": mismatches,
        "schema": SCHEMA,
        "source_changed_during_run": source_changed,
        "source_fingerprint": source_state,
        "status": "failed" if failed or mismatches or source_changed else "passed",
        "variant_workers": workers,
        "variants": results,
    }
    atomic_write_json(result_dir / "matrix.json", summary)
    print("matrix: {} ({})".format(summary["status"], result_dir / "matrix.json"), flush=True)
    return 0 if summary["status"] == "passed" else 1


def self_test():
    check_count = 0

    def check(condition, label):
        nonlocal check_count
        check_count += 1
        if not condition:
            raise MatrixError("self-test failed: " + label)

    def assert_line_numbers(source):
        return tuple(
            node.lineno for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Assert)
        )

    def failure_output_for(names):
        return b"".join(
            "{}/{} Test #{}: {} ... Passed\n".format(
                index, len(names), index, name
            ).encode("ascii")
            for index, name in enumerate(names, 1)
        )

    check(
        compact_cpu_list([3, 2, 1, 7, 9, 8]) == "1-3,7-9",
        "compact CPU ranges",
    )
    check(compact_cpu_list([]) == "", "empty CPU range")
    current_contract = evidence_contract()
    historical_contract = evidence_contract(PRE_GFNI_SCHEMA)
    pre_t16_contract = evidence_contract(PRE_T16_SCHEMA)
    check(current_contract["schema"] == SCHEMA,
          "current contract schema")
    check(len(current_contract["source_files"]) ==
          len(set(current_contract["source_files"])),
          "current contract source fingerprint is unique")
    check(set(current_contract["expected_compile_source_counts"]) <=
          set(current_contract["source_files"]),
          "current compile closure is covered by the source fingerprint")
    check(historical_contract["schema"] == PRE_GFNI_SCHEMA,
          "pre-GFNI contract schema")
    check(pre_t16_contract["schema"] == PRE_T16_SCHEMA,
          "pre-T16 contract schema")
    check("gfni_backend_failure_tests" in current_contract,
          "current contract GFNI failure tests")
    check("gfni_backend_failure_tests" not in historical_contract,
          "pre-GFNI contract excludes GFNI failure tests")
    check("gfni_backend_failure_tests" in pre_t16_contract,
          "pre-T16 contract retains GFNI failure tests")
    check("Leopard2BackendGFNI.cpp" in current_contract["source_files"],
          "current contract GFNI source")
    check("Leopard2BackendGFNI.cpp" not in
          historical_contract["source_files"],
          "pre-GFNI contract excludes GFNI source")
    check("Leopard2BackendGFNI.cpp" in
          current_contract["expected_compile_source_counts"],
          "current contract GFNI compile closure")
    check("Leopard2BackendGFNI.cpp" not in
          historical_contract["expected_compile_source_counts"],
          "pre-GFNI contract excludes GFNI compile closure")
    check("Leopard2BackendGFNI.cpp" in pre_t16_contract["source_files"],
          "pre-T16 contract retains GFNI source")
    check("Leopard2BackendAVX2T2K4.cpp" in current_contract["source_files"],
          "current contract T2/K4 source")
    check("Leopard2BackendAVX2T2K4.cpp" not in
          pre_t16_contract["source_files"],
          "pre-T16 contract excludes T2/K4 source")
    check("tests/leopard2/test_t2_k4_production.cpp" in
          current_contract["source_files"],
          "current contract T2/K4 production test source")
    check("Leopard2BackendAVX2T2K4.cpp" in
          current_contract["expected_compile_source_counts"],
          "current contract T2/K4 compile closure")
    check("tests/leopard2/test_t2_k4_production.cpp" in
          current_contract["expected_compile_source_counts"],
          "current contract T2/K4 test compile closure")
    check(current_contract["expected_compile_source_counts"].get(
              "tests/leopard2/direct_oracle.cpp") == 33,
          "current contract T2/K4 direct-oracle compile closure")
    check(historical_contract["expected_compile_source_counts"].get(
              "tests/leopard2/direct_oracle.cpp") == 14,
          "pre-GFNI contract direct-oracle compile closure")
    check(pre_t16_contract["expected_compile_source_counts"].get(
              "tests/leopard2/direct_oracle.cpp") == 14,
          "pre-T16 contract direct-oracle compile closure")
    check("Leopard2BackendAVX2T16B64.cpp" in
          current_contract["source_files"],
          "current contract T16/B64 source")
    check("Leopard2BackendAVX2T16B64.cpp" not in
          pre_t16_contract["source_files"],
          "pre-T16 contract excludes T16/B64 source")
    check("Leopard2BackendAVX2T16B64.cpp" in
          current_contract["expected_compile_source_counts"],
          "current contract T16/B64 compile closure")
    for source_name, label in (
            ("Leopard2BackendAVX2T32B256.cpp", "T32/B256"),
            ("Leopard2LowP32B64AVX2.cpp", "low P32/B64")):
        check(source_name in current_contract["source_files"],
              "current contract {} source".format(label))
        check(source_name in
              current_contract["expected_compile_source_counts"],
              "current contract {} compile closure".format(label))
        check(source_name not in pre_t16_contract["source_files"],
              "pre-T16 contract excludes {} source".format(label))
    for option in CURRENT_PRODUCTION_OBJECT_OPTIONS:
        check(option in current_contract["build_cache_keys"],
              "current contract object option provenance: " + option)
        check(option not in pre_t16_contract["build_cache_keys"],
              "pre-T16 contract excludes object option provenance: " +
              option)
    try:
        evidence_contract("leopard2-backend-matrix/v999")
    except ValueError:
        unknown_schema_rejected = True
    else:
        unknown_schema_rejected = False
    check(unknown_schema_rejected, "unknown evidence contract schema")
    check(
        normalized_output(b"a\r\nb\rc\n") == b"a\nb\nc\n",
        "normalized command output",
    )
    check(
        portable_ctest_executed(
            b"1/1 Test #1: leopard2_portable_isa ... Passed\n"
        ),
        "portable CTest recognized",
    )
    check(
        not portable_ctest_executed(b"No tests were found!!!\n"),
        "missing portable CTest rejected",
    )
    check(
        cuda_ctest_executed(
            b"1/1 Test #2: leopard2_cuda_optional ... Passed\n"
        ),
        "optional CUDA CTest recognized",
    )
    check(
        not cuda_ctest_executed(b"No tests were found!!!\n"),
        "missing optional CUDA CTest rejected",
    )
    failure_output = b"".join(
        "{}/{} Test #{}: {} ... Passed\n".format(
            index, len(BACKEND_FAILURE_TESTS), index, name
        ).encode("ascii")
        for index, name in enumerate(BACKEND_FAILURE_TESTS, 1)
    )
    check(
        ctest_executed_test_names(failure_output) == BACKEND_FAILURE_TESTS,
        "backend failure CTest names",
    )
    check(
        backend_failure_ctest_executed(failure_output),
        "complete backend failure CTest set",
    )
    check(
        backend_failure_ctest_executed(
            failure_output_for(BASE_BACKEND_FAILURE_TESTS),
            expected=BASE_BACKEND_FAILURE_TESTS,
        ),
        "base backend failure CTest set",
    )
    check(
        not backend_failure_ctest_executed(
            failure_output_for(BASE_BACKEND_FAILURE_TESTS),
        ),
        "incomplete current backend failure CTest set rejected",
    )
    check(
        backend_failure_tests_for_build({"compile_commands": []}) ==
        BASE_BACKEND_FAILURE_TESTS,
        "base backend failure tests selected",
    )
    check(
        backend_failure_tests_for_build({"compile_commands": [
            {"file": "Leopard2BackendAVX512.cpp"},
        ]}) == BASE_BACKEND_FAILURE_TESTS + AVX512_BACKEND_FAILURE_TESTS,
        "AVX-512 backend failure tests selected",
    )
    check(
        backend_failure_tests_for_build({"compile_commands": [
            {"file": "Leopard2BackendGFNI.cpp"},
        ]}) == BASE_BACKEND_FAILURE_TESTS + GFNI_BACKEND_FAILURE_TESTS,
        "GFNI backend failure tests selected",
    )
    check(
        backend_failure_tests_for_build({"compile_commands": [
            {"file": "Leopard2BackendAVX512.cpp"},
            {"file": "Leopard2BackendGFNI.cpp"},
        ]}) == BACKEND_FAILURE_TESTS,
        "all backend failure tests selected",
    )

    check(
        backend_failure_ctest_executed(
            failure_output_for(tuple(reversed(BACKEND_FAILURE_TESTS)))
        ),
        "backend failure CTest order independence",
    )
    mutations = (
        BACKEND_FAILURE_TESTS[:-1],
        BACKEND_FAILURE_TESTS[:4] + BACKEND_FAILURE_TESTS[5:],
        BACKEND_FAILURE_TESTS[:-1] + (
            "leopard2_backend_failure_avx2_kat_renamed",
        ),
        BACKEND_FAILURE_TESTS + (
            "leopard2_backend_failure_scalar_kat_extra",
        ),
        BACKEND_FAILURE_TESTS[:-1] + (BACKEND_FAILURE_TESTS[0],),
        (BACKEND_FAILURE_TESTS[0],),
    )
    for index, mutation in enumerate(mutations):
        check(
            not backend_failure_ctest_executed(failure_output_for(mutation)),
            "backend failure CTest mutation {}".format(index),
        )
    check(
        pinned_command(["ctest", "--test-dir", "build"], None, 7) == [
            "ctest", "--test-dir", "build",
        ],
        "unpinned command",
    )
    check(
        pinned_command(
            ["ctest", "--test-dir", "build"], "/usr/bin/taskset", 23
        ) == [
            "/usr/bin/taskset", "-c", "23", "ctest", "--test-dir", "build",
        ],
        "pinned command",
    )
    check(
        digest_value({"b": 2, "a": 1}) == digest_value({"a": 1, "b": 2}),
        "canonical mapping digest",
    )
    check(variant_flags("scalar") == [], "scalar variant flags")
    check(
        variant_flags("ssse3") == ["-mssse3", "-mno-avx"],
        "SSSE3 variant flags",
    )
    check(
        variant_flags("avx2") == ["-mavx2", "-mno-avx512f"],
        "AVX2 variant flags",
    )
    check(
        variant_flags("avx512") == [
            "-mavx2", "-mavx512f", "-mavx512bw", "-mavx512vl",
            "-mprefer-vector-width=256",
        ],
        "AVX-512 variant flags",
    )
    check(
        availability(
            "scalar",
            {"architecture": "x86_64", "cpu_flags": []},
            {"executable": "not-needed-for-empty-flags"},
        ) == (True, ""),
        "scalar availability",
    )
    check(
        availability(
            "avx512",
            {"architecture": "x86_64", "cpu_flags": ["avx2"]},
            {"executable": "not-reached-when-host-flags-are-missing"},
        ) == (
            False,
            "host CPU does not advertise avx512f,avx512bw,avx512vl",
        ),
        "AVX-512 missing-host-feature rejection",
    )
    require_compile_source_counts(
        EXPECTED_COMPILE_SOURCE_COUNTS,
        EXPECTED_COMPILE_SOURCE_COUNTS,
    )
    check(True, "unchanged compile-source contract")
    no_isa = expected_compile_source_counts({})
    check("Leopard2BackendSSSE3.cpp" not in no_isa,
          "SSSE3 object excluded without compiler probe")
    check("Leopard2BackendAVX2.cpp" not in no_isa,
          "AVX2 object excluded without compiler probe")
    check("Leopard2BackendAVX2Xor.cpp" not in no_isa,
          "AVX2 XOR object excluded without compiler probe")
    check("Leopard2BackendAVX2T2K4.cpp" not in no_isa,
          "AVX2 T2/K4 object excluded without compiler probe")
    check("Leopard2BackendAVX2T16B64.cpp" not in no_isa,
          "AVX2 T16/B64 object excluded without compiler probe")
    check("Leopard2BackendAVX2T32B256.cpp" not in no_isa,
          "AVX2 T32/B256 object excluded without compiler probe")
    check("Leopard2LowP32B64AVX2.cpp" not in no_isa,
          "AVX2 low P32/B64 object excluded without compiler probe")
    check("tests/leopard2/test_t2_k4_production.cpp" not in no_isa,
          "AVX2 T2/K4 test excluded without compiler probe")
    check(no_isa.get("tests/leopard2/direct_oracle.cpp") == 14,
          "T2/K4 direct-oracle compile excluded without compiler probe")
    check("Leopard2BackendAVX512.cpp" not in no_isa,
          "AVX-512 object excluded without compiler probe")
    check("Leopard2BackendGFNI.cpp" not in no_isa,
          "GFNI object excluded without compiler probe")
    avx2_only = expected_compile_source_counts({
        "LEO2_FLAG_MSSSE3": "1",
        "LEO2_FLAG_MNO_AVX": "1",
        "LEO2_FLAG_MAVX2": "1",
        "LEO2_FLAG_MNO_AVX512F": "1",
    })
    check(avx2_only.get("Leopard2BackendSSSE3.cpp") == 2,
          "SSSE3 compile-source count")
    check(avx2_only.get("Leopard2BackendAVX2.cpp") == 2,
          "AVX2 compile-source count")
    check(avx2_only.get("Leopard2BackendAVX2Xor.cpp") == 2,
          "AVX2 XOR compile-source count")
    check(avx2_only.get("Leopard2BackendAVX2T2K4.cpp") == 1,
          "AVX2 T2/K4 compile-source count")
    check(avx2_only.get("Leopard2BackendAVX2T16B64.cpp") == 1,
          "AVX2 T16/B64 compile-source count")
    check(avx2_only.get("Leopard2BackendAVX2T32B256.cpp") == 1,
          "AVX2 T32/B256 compile-source count")
    check(avx2_only.get("Leopard2LowP32B64AVX2.cpp") == 1,
          "AVX2 low P32/B64 compile-source count")
    check(avx2_only.get("tests/leopard2/test_t2_k4_production.cpp") == 1,
          "AVX2 T2/K4 test compile-source count")
    check(avx2_only.get("tests/leopard2/direct_oracle.cpp") == 33,
          "AVX2 T2/K4 direct-oracle compile-source count")
    check("Leopard2BackendAVX512.cpp" not in avx2_only,
          "AVX-512 object excluded from AVX2-only probes")
    check("Leopard2BackendGFNI.cpp" not in avx2_only,
          "GFNI object excluded from AVX2-only probes")
    t16_off = expected_compile_source_counts({
        "LEO2_FLAG_MSSSE3": "1",
        "LEO2_FLAG_MNO_AVX": "1",
        "LEO2_FLAG_MAVX2": "1",
        "LEO2_FLAG_MNO_AVX512F": "1",
        "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "OFF",
    })
    check("Leopard2BackendAVX2T16B64.cpp" not in t16_off,
          "T16/B64 object excluded when experiment is disabled")
    check(t16_off.get("Leopard2BackendAVX2T2K4.cpp") == 1,
          "T2/K4 object retained when T16/B64 is disabled")
    t32_off = expected_compile_source_counts({
        "LEO2_FLAG_MSSSE3": "1",
        "LEO2_FLAG_MNO_AVX": "1",
        "LEO2_FLAG_MAVX2": "1",
        "LEO2_FLAG_MNO_AVX512F": "1",
        "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
        "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "OFF",
    })
    check("Leopard2BackendAVX2T32B256.cpp" not in t32_off,
          "T32/B256 object excluded when both selectors are disabled")
    p32_off = expected_compile_source_counts({
        "LEO2_FLAG_MSSSE3": "1",
        "LEO2_FLAG_MNO_AVX": "1",
        "LEO2_FLAG_MAVX2": "1",
        "LEO2_FLAG_MNO_AVX512F": "1",
        "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "OFF",
    })
    check("Leopard2LowP32B64AVX2.cpp" not in p32_off,
          "low P32/B64 object excluded when its selector is disabled")
    gfni_only = expected_compile_source_counts({
        "LEO2_FLAG_MSSSE3": "1",
        "LEO2_FLAG_MNO_AVX": "1",
        "LEO2_FLAG_MAVX2": "1",
        "LEO2_FLAG_MNO_AVX512F": "1",
        "LEO2_FLAG_MGFNI": "1",
    })
    check(gfni_only.get("Leopard2BackendGFNI.cpp") == 2,
          "GFNI compile-source count")
    check("Leopard2BackendAVX512.cpp" not in gfni_only,
          "AVX-512 object excluded from GFNI-only probes")
    gf16_only = expected_compile_source_counts({
        "LEO2_FLAG_MSSSE3": "1",
        "LEO2_FLAG_MNO_AVX": "1",
        "LEO2_FLAG_MAVX2": "1",
        "LEO2_FLAG_MNO_AVX512F": "1",
        "LEOPARD_ENABLE_GF8": "OFF",
        "LEOPARD_ENABLE_GF16": "ON",
    })
    check("Leopard2BackendAVX2Xor.cpp" not in gf16_only,
          "GF8-only AVX2 XOR object excluded from GF16 build")
    check("Leopard2BackendAVX2T2K4.cpp" not in gf16_only,
          "GF8-only AVX2 T2/K4 object excluded from GF16 build")
    check("Leopard2BackendAVX2T16B64.cpp" not in gf16_only,
          "GF8-only AVX2 T16/B64 object excluded from GF16 build")
    check("Leopard2BackendAVX2T32B256.cpp" not in gf16_only,
          "GF8-only AVX2 T32/B256 object excluded from GF16 build")
    check("Leopard2LowP32B64AVX2.cpp" not in gf16_only,
          "GF8-only AVX2 low P32/B64 object excluded from GF16 build")
    check("tests/leopard2/test_t2_k4_production.cpp" not in gf16_only,
          "GF8-only AVX2 T2/K4 test excluded from GF16 build")
    check(gf16_only.get("tests/leopard2/direct_oracle.cpp") == 14,
          "GF8-only T2/K4 direct-oracle compile excluded from GF16 build")
    check(gf16_only.get("Leopard2BackendAVX2.cpp") == 2,
          "AVX2 object retained in GF16 build")
    all_isa = expected_compile_source_counts({
        "LEO2_FLAG_MSSSE3": "1",
        "LEO2_FLAG_MNO_AVX": "1",
        "LEO2_FLAG_MAVX2": "1",
        "LEO2_FLAG_MNO_AVX512F": "1",
        "LEO2_FLAG_MAVX512F": "1",
        "LEO2_FLAG_MAVX512BW": "1",
        "LEO2_FLAG_MAVX512VL": "1",
        "LEO2_FLAG_MPREFER_VECTOR_WIDTH_256": "1",
        "LEO2_FLAG_MGFNI": "1",
    })
    check(all_isa == EXPECTED_COMPILE_SOURCE_COUNTS,
          "all-ISA compile-source closure")
    contract_mutations = []
    missing = dict(EXPECTED_COMPILE_SOURCE_COUNTS)
    missing.pop("tests/leopard2/test_batch_aliasing.cpp")
    contract_mutations.append(missing)
    duplicated = dict(EXPECTED_COMPILE_SOURCE_COUNTS)
    duplicated["Leopard2Backend.cpp"] += 1
    contract_mutations.append(duplicated)
    unexpected = dict(EXPECTED_COMPILE_SOURCE_COUNTS)
    unexpected["tests/leopard2/untracked_backend_test.cpp"] = 1
    contract_mutations.append(unexpected)
    for index, mutation in enumerate(contract_mutations):
        try:
            require_compile_source_counts(
                mutation,
                EXPECTED_COMPILE_SOURCE_COUNTS,
            )
        except MatrixError as error:
            rejected = str(error) == \
                "compile-command source multiset mismatch"
        else:
            rejected = False
        check(
            rejected,
            "compile-source contract mutation {}".format(index),
        )
    with tempfile.TemporaryDirectory(prefix="leo2-backend-self-test-") as directory:
        path = Path(directory) / "result.json"
        value = {"z": [3, 2, 1], "a": "stable"}
        atomic_write_json(path, value)
        check(
            json.loads(path.read_text(encoding="utf-8")) == value,
            "atomic JSON write",
        )
        stale = Path(directory) / "stale-build"
        stale.mkdir()
        (stale / "copied-object.o").write_bytes(b"stale")
        check(prepare_fresh_directory(stale), "fresh build directory prepared")
        check(list(stale.iterdir()) == [], "stale build artifact removed")

        source = Path(directory) / "source"
        for relative in SOURCE_FILES:
            source_path = source / relative
            source_path.parent.mkdir(parents=True, exist_ok=True)
            source_path.write_bytes((relative + "\n").encode("utf-8"))
        source_before = source_fingerprint(source)
        check(
            source_fingerprint(source) == source_before,
            "stable source fingerprint",
        )
        changed_source = source / SOURCE_FILES[0]
        changed_source.write_bytes(changed_source.read_bytes() + b"mutated\n")
        source_after = source_fingerprint(source)
        check(
            source_after["digest"] != source_before["digest"],
            "source-content mutation changes fingerprint",
        )
        check(
            source_after["files"][SOURCE_FILES[0]] !=
            source_before["files"][SOURCE_FILES[0]],
            "mutated source file hash changes",
        )
        changed_source.unlink()
        try:
            source_fingerprint(source)
        except MatrixError as error:
            missing_source_rejected = \
                str(error).startswith("required matrix input is missing:")
        else:
            missing_source_rejected = False
        check(missing_source_rejected, "missing source input rejected")

        artifact = Path(directory) / "fixture-compiler"
        artifact.write_text(
            "#!/bin/sh\nprintf 'fixture compiler 1.0\\n'\n",
            encoding="utf-8",
        )
        artifact.chmod(0o755)
        artifact_before = compiler_identity(str(artifact))
        artifact.write_text(
            "#!/bin/sh\nprintf 'fixture compiler 1.0\\n'\n# mutation\n",
            encoding="utf-8",
        )
        artifact_after = compiler_identity(str(artifact))
        check(
            artifact_before["binary_sha256"] !=
            artifact_after["binary_sha256"],
            "compiler artifact mutation changes identity",
        )
        check(
            artifact_before["version_sha256"] ==
            artifact_after["version_sha256"],
            "artifact hash is independent of version text",
        )

        driver_link = Path(directory) / "python-driver"
        driver_link.symlink_to(sys.executable)
        driver_identity = compiler_identity(str(driver_link))
        check(
            driver_identity["executable"] == str(driver_link.absolute()),
            "compiler driver spelling retained",
        )
        check(
            Path(driver_identity["executable"]).resolve() ==
            Path(sys.executable).resolve(),
            "compiler driver symlink target",
        )
        timeout_record = run_command(
            "timeout", [sys.executable, "-c", "import time; time.sleep(1)"],
            directory, Path(directory), 0.01
        )
        check(timeout_record["returncode"] == 124, "timeout return code")
        check(timeout_record["timed_out"] is True, "timeout flag")

    check(
        assert_line_numbers("value = 1\n") == (),
        "assert scanner accepts assert-free source",
    )
    check(
        assert_line_numbers("value = 1\nassert value\n") == (2,),
        "assert scanner detects optimized-away check",
    )
    check(
        assert_line_numbers(Path(__file__).read_text(encoding="utf-8")) == (),
        "backend matrix contains no Python assert statements",
    )
    check(check_count > 1, "self-test executed a nontrivial check set")
    print(
        "leopard2 backend matrix self-test passed ({} checks)".format(
            check_count
        )
    )
    return 0


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command")
    run = subparsers.add_parser("run", help="build and compare backend variants")
    default_source = str(Path(__file__).resolve().parents[1])
    run.add_argument("--source", default=default_source)
    run.add_argument("--build-root", default="build/leopard2-backend-matrix")
    run.add_argument("--result-dir", default="results/leopard2/backend-matrix")
    run.add_argument("--variants", default=",".join(VARIANTS))
    run.add_argument("--jobs", type=int, default=min(os.cpu_count() or 1, 128))
    run.add_argument("--variant-workers", type=int, default=len(VARIANTS))
    run.add_argument("--timeout", type=int, default=900)
    run.add_argument("--compiler", default=os.environ.get("CXX", "c++"))
    run.add_argument("--c-compiler", default=os.environ.get("CC", "cc"))
    run.add_argument("--cmake", default="cmake")
    run.add_argument("--ctest", default="ctest")
    run.add_argument("--generator", default="Ninja" if shutil.which("ninja") else "Unix Makefiles")
    run.add_argument("--no-resume", action="store_true")
    subparsers.add_parser("self-test", help="run deterministic utility tests")
    return result


def main():
    arguments = parser().parse_args()
    try:
        if arguments.command == "self-test":
            return self_test()
        if arguments.command == "run":
            if arguments.jobs <= 0 or arguments.variant_workers <= 0 or arguments.timeout <= 0:
                raise MatrixError("jobs, variant-workers, and timeout must be positive")
            return matrix_run(arguments)
        parser().print_help()
        return 2
    except (MatrixError, OSError, subprocess.SubprocessError) as error:
        print("backend matrix error: {}".format(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
