#!/usr/bin/env python3
"""Generate and advance the exact-main-first balanced decoder matrix."""

from __future__ import annotations

import argparse
import base64
import binascii
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import stat
import sys
import tempfile
from typing import Any, Iterable, Mapping, Sequence

import balanced_evidence_common as common


PLAN_SCHEMA = "leopard2-balanced-promotion-plan/v2"
MAIN_CELL_SCHEMA = "leopard2-main-compare-cell-list/v2"
SURVIVOR_SCHEMA = "leopard2-balanced-promotion-survivors/v4"
STAGE_SCHEMA = "leopard2-balanced-promotion-stage/v5"
ATTESTATION_SCHEMA = "leopard2-balanced-auto-path-attestation/v5"
ATTESTATION_RESULT_SCHEMA_V5 = "leopard2-balanced-auto-path-result/v5"
ATTESTATION_RESULT_SCHEMA_V6 = "leopard2-balanced-auto-path-result/v6"
ATTESTATION_RESULT_SCHEMA_V7 = "leopard2-balanced-auto-path-result/v7"
ATTESTATION_RESULT_SCHEMA_V8 = "leopard2-balanced-auto-path-result/v8"
ATTESTATION_RESULT_SCHEMA_V9 = "leopard2-balanced-auto-path-result/v9"
ATTESTATION_RESULT_SCHEMA_V10 = "leopard2-balanced-auto-path-result/v10"
ATTESTATION_RESULT_SCHEMA_V11 = "leopard2-balanced-auto-path-result/v11"
ATTESTATION_RESULT_SCHEMA_V12 = "leopard2-balanced-auto-path-result/v12"
ATTESTATION_RESULT_SCHEMA_V13 = "leopard2-balanced-auto-path-result/v13"
ATTESTATION_RESULT_SCHEMA = "leopard2-balanced-auto-path-result/v14"
PROMOTION_TIMING_SCHEMA = "leopard2-balanced-promotion-timing-evidence/v1"
EXACT_MANIFEST_SCHEMA_V5 = "leopard2-main-compare-manifest/v5"
EXACT_MANIFEST_SCHEMA_V6 = "leopard2-main-compare-manifest/v6"
EXACT_MANIFEST_SCHEMA_V7 = "leopard2-main-compare-manifest/v7"
EXACT_MANIFEST_SCHEMA_V8 = "leopard2-main-compare-manifest/v8"
EXACT_MANIFEST_SCHEMA_V9 = "leopard2-main-compare-manifest/v9"
EXACT_MANIFEST_SCHEMA_V10 = "leopard2-main-compare-manifest/v10"
EXACT_MANIFEST_SCHEMA_V11 = "leopard2-main-compare-manifest/v11"
EXACT_MANIFEST_SCHEMA_V12 = "leopard2-main-compare-manifest/v12"
EXACT_MANIFEST_SCHEMA_V13 = "leopard2-main-compare-manifest/v13"
EXACT_MANIFEST_SCHEMA_V14 = "leopard2-main-compare-manifest/v14"
EXACT_MANIFEST_SCHEMA_V15 = "leopard2-main-compare-manifest/v15"
# v17 is an encode-only GFNI campaign.  It cannot satisfy this decoder
# promotion planner's paired decode evidence contract, so v16 remains the
# newest admitted exact-main generation.
EXACT_MANIFEST_SCHEMA = "leopard2-main-compare-manifest/v16"
EXACT_MANIFEST_SCHEMAS = frozenset((
    EXACT_MANIFEST_SCHEMA_V5, EXACT_MANIFEST_SCHEMA_V6,
    EXACT_MANIFEST_SCHEMA_V7, EXACT_MANIFEST_SCHEMA_V8,
    EXACT_MANIFEST_SCHEMA_V9, EXACT_MANIFEST_SCHEMA_V10,
    EXACT_MANIFEST_SCHEMA_V11, EXACT_MANIFEST_SCHEMA_V12,
    EXACT_MANIFEST_SCHEMA_V13, EXACT_MANIFEST_SCHEMA_V14,
    EXACT_MANIFEST_SCHEMA_V15, EXACT_MANIFEST_SCHEMA,
))
EXACT_MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
EXACT_MAIN_TREE = "b7c8830d96a978f6ec14fe747095f066e351ae72"
# Exact commit payload used only by deterministic scope self-tests.  Retaining
# the signed commit bytes makes the fixture exercise the same object-ID/tree
# proof as real complete-schema evidence rather than substituting a plausible hash.
EXACT_MAIN_COMMIT_BASE64 = (
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
PROMOTION_GATE = 1.05
NEIGHBOR_MAXIMUM_REGRESSION = 0.02
FORCED_SURVIVOR_MAXIMUM_RATIO = 1.0 + NEIGHBOR_MAXIMUM_REGRESSION
BOUNDARY_K = (
    5, 7, 8, 9, 14, 15, 16, 17, 29, 30, 31, 32, 33, 62, 63, 64,
    65, 96, 112, 120, 124, 125, 126, 127, 128,
)
EXACT_GATE_BYTES = (256, 4096, 65536)
ALIGNED_CONFIRM_BYTES = (192, 4032, 4096, 65536, 1024 * 1024)
TRUE_TAIL_BYTES = (193, 4033, 4097, 1024 * 1024 + 1)
BACKENDS = ("scalar", "ssse3", "avx2")
MODE_PAIRS = (("generic", "materialized"), ("generic", "tiled"))
CONFIRM_MODES = ("generic", "auto", "materialized", "tiled")
EXACT_ORDER = ("baseline", "candidate", "candidate", "baseline")
FORCED_ORDERS = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
# This promotion campaign has same-binary forced confirmation for exactly these
# backends.  AVX-512 and NEON are intentionally rejected rather than borrowing
# generic speed evidence collected under a kernel that the confirmation runner
# cannot force.  A later schema can add them with their own forced matrices.
CAMPAIGN_BACKENDS = ("scalar", "ssse3", "avx2")
AUTO_BACKENDS = CAMPAIGN_BACKENDS
EXACT_RAW_SCHEMA_V5 = "leopard2-main-compare-raw/v5"
EXACT_RAW_SCHEMA_V6 = "leopard2-main-compare-raw/v6"
EXACT_RAW_SCHEMA_V7 = "leopard2-main-compare-raw/v7"
EXACT_RAW_SCHEMA_V8 = "leopard2-main-compare-raw/v8"
EXACT_RAW_SCHEMA_V9 = "leopard2-main-compare-raw/v9"
EXACT_RAW_SCHEMA_V10 = "leopard2-main-compare-raw/v10"
EXACT_RAW_SCHEMA_V11 = "leopard2-main-compare-raw/v11"
EXACT_RAW_SCHEMA_V12 = "leopard2-main-compare-raw/v12"
EXACT_RAW_SCHEMA_V13 = "leopard2-main-compare-raw/v13"
EXACT_RAW_SCHEMA_V14 = "leopard2-main-compare-raw/v14"
EXACT_RAW_SCHEMA_V15 = "leopard2-main-compare-raw/v15"
EXACT_RAW_SCHEMA = "leopard2-main-compare-raw/v16"
EXACT_RAW_SCHEMAS = frozenset((
    EXACT_RAW_SCHEMA_V5, EXACT_RAW_SCHEMA_V6, EXACT_RAW_SCHEMA_V7,
    EXACT_RAW_SCHEMA_V8, EXACT_RAW_SCHEMA_V9, EXACT_RAW_SCHEMA_V10,
    EXACT_RAW_SCHEMA_V11, EXACT_RAW_SCHEMA_V12, EXACT_RAW_SCHEMA_V13,
    EXACT_RAW_SCHEMA_V14, EXACT_RAW_SCHEMA_V15, EXACT_RAW_SCHEMA,
))
EXACT_SCHEMA_PAIRS = frozenset((
    (EXACT_MANIFEST_SCHEMA_V5, EXACT_RAW_SCHEMA_V5),
    (EXACT_MANIFEST_SCHEMA_V6, EXACT_RAW_SCHEMA_V6),
    (EXACT_MANIFEST_SCHEMA_V7, EXACT_RAW_SCHEMA_V7),
    (EXACT_MANIFEST_SCHEMA_V8, EXACT_RAW_SCHEMA_V8),
    (EXACT_MANIFEST_SCHEMA_V9, EXACT_RAW_SCHEMA_V9),
    (EXACT_MANIFEST_SCHEMA_V10, EXACT_RAW_SCHEMA_V10),
    (EXACT_MANIFEST_SCHEMA_V11, EXACT_RAW_SCHEMA_V11),
    (EXACT_MANIFEST_SCHEMA_V12, EXACT_RAW_SCHEMA_V12),
    (EXACT_MANIFEST_SCHEMA_V13, EXACT_RAW_SCHEMA_V13),
    (EXACT_MANIFEST_SCHEMA_V14, EXACT_RAW_SCHEMA_V14),
    (EXACT_MANIFEST_SCHEMA_V15, EXACT_RAW_SCHEMA_V15),
    (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA),
))
CANONICAL_LDD_SCHEMA = "leopard2-main-compare-canonical-ldd/v1"
CANONICAL_LDD_NORMALIZATION = "terminal-aslr-load-address/v1"
CANONICAL_LDD_ADDRESS = "<ASLR_LOAD_ADDRESS>"
MAX_GIT_COMMIT_BYTES = 1024 * 1024
MAX_RETAINED_TEXT_BYTES = 1024 * 1024
MAX_BUILD_CONFIGURATION_BYTES = 64 * 1024
MAX_BUILD_TOOL_VERSION_BYTES = 64 * 1024
MAX_NINJA_GRAPH_FILES = 64
MAX_NINJA_GRAPH_TOTAL_BYTES = 8 * 1024 * 1024
MAX_JSON_DOCUMENT_BYTES = 64 * 1024 * 1024
NINJA_GRAPH_CLOSURE_SCHEMA = "leopard2-ninja-graph-closure/v1"
MAX_CPU_ID = 1_048_575
MAX_CPU_LIST_ENTRIES = 4096
BASELINE_LIBRARY_SOURCES = (
    "leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp",
)
CANDIDATE_LIBRARY_SOURCES_V9 = (
    "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp", "Leopard2Plan.cpp",
    "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp",
    "Leopard2BackendSSSE3.cpp", "Leopard2BackendAVX2.cpp",
    "Leopard2BackendAVX2Xor.cpp",
    "Leopard2BackendAVX512.cpp", "Leopard2BackendGFNI.cpp",
)
CANDIDATE_LIBRARY_SOURCES_V11 = tuple(
    source for source in CANDIDATE_LIBRARY_SOURCES_V9
    if source != "Leopard2BackendAVX512.cpp")
CANDIDATE_LIBRARY_SOURCES_V12 = (
    *CANDIDATE_LIBRARY_SOURCES_V11[:-1],
    "Leopard2BackendAVX2T32B256.cpp",
    "Leopard2BackendAVX2T16B64.cpp",
    "Leopard2LowP32B64AVX2.cpp",
    "Leopard2BackendAVX2T2K4.cpp",
    CANDIDATE_LIBRARY_SOURCES_V11[-1],
)
CANDIDATE_LIBRARY_SOURCES = (
    *CANDIDATE_LIBRARY_SOURCES_V12[:-1],
    "Leopard2BackendAVX2T8K8B1024.cpp",
    CANDIDATE_LIBRARY_SOURCES_V12[-1],
)
BASELINE_EXPECTED_COMPILE_COMMAND_COUNT = 5
CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT_V11 = 22
CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT_V12 = 26
CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT = 27
COMPILE_COMMANDS_SCHEMA_V2 = "leopard2-main-compare-compile-commands/v2"
COMPILE_COMMANDS_SCHEMA_V3 = "leopard2-main-compare-compile-commands/v3"
COMPILE_COMMANDS_SCHEMA_V4 = "leopard2-main-compare-compile-commands/v4"
COMPILE_COMMANDS_SCHEMA_V5 = "leopard2-main-compare-compile-commands/v5"
COMPILE_COMMANDS_SCHEMA_V6 = "leopard2-main-compare-compile-commands/v6"
COMPILE_COMMANDS_SCHEMA_V7 = "leopard2-main-compare-compile-commands/v7"
COMPILE_COMMANDS_SCHEMA_V8 = "leopard2-main-compare-compile-commands/v8"
COMPILE_COMMANDS_SCHEMA_V9 = "leopard2-main-compare-compile-commands/v9"
COMPILE_COMMANDS_SCHEMA_V10 = "leopard2-main-compare-compile-commands/v10"
COMPILE_COMMANDS_SCHEMA_V11 = "leopard2-main-compare-compile-commands/v11"
COMPILE_COMMANDS_SCHEMA = "leopard2-main-compare-compile-commands/v12"
GNU_CXX_DRIVER_BASENAME = re.compile(
    r"(?:g\+\+|(?:[A-Za-z0-9_.+]+-)+"
    r"(?:gnu|gnueabi(?:hf)?|eabi|elf|musl|mingw32)-g\+\+)"
    r"(?:-[0-9]+(?:\.[0-9]+)*)?")
BASELINE_COMPILE_PROFILE = \
    "gnu-compatible-cxx11-native-x86_64-release/v1"
BASELINE_PURE_AVX2_COMPILE_PROFILE = \
    "gnu-compatible-cxx11-pure-avx2-x86_64-release/v1"
CANDIDATE_COMPILE_PROFILE_V1 = \
    "gnu-compatible-cxx11-runtime-dispatch-x86_64-release/v1"
CANDIDATE_COMPILE_PROFILE_V2 = \
    "gnu-compatible-cxx11-runtime-dispatch-x86_64-release/v2"
CANDIDATE_COMPILE_PROFILE_V3 = \
    "gnu-compatible-cxx11-runtime-dispatch-effective-avx2-x86_64-release/v3"
CANDIDATE_COMPILE_PROFILE_V4 = \
    "gnu-compatible-cxx11-runtime-dispatch-effective-avx2-x86_64-release/v4"
CANDIDATE_COMPILE_PROFILE_V5 = \
    "gnu-compatible-cxx11-runtime-dispatch-effective-avx2-x86_64-release/v5"
CANDIDATE_COMPILE_PROFILE_V6 = \
    "gnu-compatible-cxx11-runtime-dispatch-effective-avx2-x86_64-release/v6"
CANDIDATE_COMPILE_PROFILE_V7 = \
    "gnu-compatible-cxx11-runtime-dispatch-effective-avx2-x86_64-release/v7"
CANDIDATE_COMPILE_PROFILE = \
    "gnu-compatible-cxx11-runtime-dispatch-effective-avx2-x86_64-release/v8"
BASELINE_NATIVE_ISA_POLICY = "whole-build -march=native"
BASELINE_PURE_AVX2_ISA_POLICY = (
    "whole-build -march=x86-64 -mtune=generic -mavx2 -mno-avx512f")
CANDIDATE_ISA_POLICY_V9 = (
    "portable core with ISA flags only on SSSE3, AVX2, and "
    "AVX-512VL translation units")
CANDIDATE_ISA_POLICY_V11 = (
    "portable core with ISA flags only on SSSE3, AVX2, and GFNI "
    "translation units; AVX-512 probes disabled; AUTO resolved AVX2")
CANDIDATE_ISA_POLICY_V12 = (
    "portable core with ISA flags only on SSSE3, AVX2, isolated generated "
    "AVX2, and GFNI translation units; AVX-512 probes disabled; AUTO "
    "resolved AVX2")
CANDIDATE_ISA_POLICY = (
    "portable core with ISA flags only on SSSE3, AVX2, isolated generated "
    "and fixed-shape K8 AVX2, and GFNI translation units; AVX-512 probes "
    "disabled; AUTO resolved AVX2")
BUILD_CONFIGURATION_RECORD_SCHEMA_V2 = \
    "leopard2-main-compare-build-configuration/v2"
BUILD_CONFIGURATION_RECORD_SCHEMA_V3 = \
    "leopard2-main-compare-build-configuration/v3"
BUILD_CONFIGURATION_RECORD_SCHEMA_V4 = \
    "leopard2-main-compare-build-configuration/v4"
BUILD_CONFIGURATION_RECORD_SCHEMA_V5 = \
    "leopard2-main-compare-build-configuration/v5"
BUILD_CONFIGURATION_RECORD_SCHEMA_V6 = \
    "leopard2-main-compare-build-configuration/v6"
BUILD_CONFIGURATION_RECORD_SCHEMA_V7 = \
    "leopard2-main-compare-build-configuration/v7"
BUILD_CONFIGURATION_RECORD_SCHEMA_V8 = \
    "leopard2-main-compare-build-configuration/v8"
BUILD_CONFIGURATION_RECORD_SCHEMA = \
    "leopard2-main-compare-build-configuration/v9"
BUILD_CONFIGURATION_FILE_SCHEMA_V2 = \
    "leopard2-benchmark-build-configuration/v2"
BUILD_CONFIGURATION_FILE_SCHEMA_V3 = \
    "leopard2-benchmark-build-configuration/v3"
BUILD_CONFIGURATION_FILE_SCHEMA_V4 = \
    "leopard2-benchmark-build-configuration/v4"
BUILD_CONFIGURATION_FILE_SCHEMA_V5 = \
    "leopard2-benchmark-build-configuration/v5"
BUILD_CONFIGURATION_FILE_SCHEMA_V6 = \
    "leopard2-benchmark-build-configuration/v6"
BUILD_CONFIGURATION_FILE_SCHEMA_V7 = \
    "leopard2-benchmark-build-configuration/v7"
BUILD_CONFIGURATION_FILE_SCHEMA_V8 = \
    "leopard2-benchmark-build-configuration/v8"
BUILD_CONFIGURATION_FILE_SCHEMA = \
    "leopard2-benchmark-build-configuration/v9"
BUILD_CONFIGURATION_RELATIVE_PATH = (
    "generated/leopard2-benchmark-attestation/"
    "leopard2_benchmark_build_configuration.txt"
)
BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH = \
    "cmake/Leopard2BenchmarkAttestation.cmake"
EVIDENCE_HELPER_RELATIVE_PATH = \
    "experiments/leopard2/decoder_dispatch/balanced_evidence_common.py"
BALANCED_ANALYZER_RELATIVE_PATH = \
    "experiments/leopard2/decoder_dispatch/analyze_balanced.py"
EVIDENCE_SCOPE_SCHEMA_V3 = "leopard2-balanced-evidence-scope/v3"
EVIDENCE_SCOPE_SCHEMA_V4 = "leopard2-balanced-evidence-scope/v4"
EVIDENCE_SCOPE_SCHEMA_V5 = "leopard2-balanced-evidence-scope/v5"
EVIDENCE_SCOPE_SCHEMA_V6 = "leopard2-balanced-evidence-scope/v6"
EVIDENCE_SCOPE_SCHEMA_V7 = "leopard2-balanced-evidence-scope/v7"
EVIDENCE_SCOPE_SCHEMA_V8 = "leopard2-balanced-evidence-scope/v8"
EVIDENCE_SCOPE_SCHEMA_V9 = "leopard2-balanced-evidence-scope/v9"
EVIDENCE_SCOPE_SCHEMA_V10 = "leopard2-balanced-evidence-scope/v10"
EVIDENCE_SCOPE_SCHEMA_V11 = "leopard2-balanced-evidence-scope/v11"
EVIDENCE_SCOPE_SCHEMA = "leopard2-balanced-evidence-scope/v12"
CANONICAL_BUILD_VALIDATOR_V2 = \
    "exact-main/run_abba.py build_provenance schema v8"
CANONICAL_BUILD_VALIDATOR_V3 = \
    "exact-main/run_abba.py build_provenance schema v9"
CANONICAL_BUILD_VALIDATOR_V4 = \
    "exact-main/run_abba.py build_provenance schema v10"
CANONICAL_BUILD_VALIDATOR_V5 = \
    "exact-main/run_abba.py build_provenance schema v11"
CANONICAL_BUILD_VALIDATOR_V6 = \
    "exact-main/run_abba.py build_provenance schema v12"
CANONICAL_BUILD_VALIDATOR_V7 = \
    "exact-main/run_abba.py build_provenance schema v13"
CANONICAL_BUILD_VALIDATOR_V8 = \
    "exact-main/run_abba.py build_provenance schema v14"
CANONICAL_BUILD_VALIDATOR_V9 = \
    "exact-main/run_abba.py build_provenance schema v15"
CANONICAL_BUILD_VALIDATOR = \
    "exact-main/run_abba.py build_provenance schema v16"
CANONICAL_PRODUCTION_BUILD_SCHEMA_V2 = \
    "leopard2-canonical-production-build/v2"
CANONICAL_PRODUCTION_BUILD_SCHEMA_V3 = \
    "leopard2-canonical-production-build/v3"
CANONICAL_PRODUCTION_BUILD_SCHEMA_V4 = \
    "leopard2-canonical-production-build/v4"
CANONICAL_PRODUCTION_BUILD_SCHEMA_V5 = \
    "leopard2-canonical-production-build/v5"
CANONICAL_PRODUCTION_BUILD_SCHEMA_V6 = \
    "leopard2-canonical-production-build/v6"
CANONICAL_PRODUCTION_BUILD_SCHEMA_V7 = \
    "leopard2-canonical-production-build/v7"
CANONICAL_PRODUCTION_BUILD_SCHEMA_V8 = \
    "leopard2-canonical-production-build/v8"
CANONICAL_PRODUCTION_BUILD_SCHEMA_V9 = \
    "leopard2-canonical-production-build/v9"
CANONICAL_PRODUCTION_BUILD_SCHEMA = \
    "leopard2-canonical-production-build/v10"
BUILD_CONFIGURATION_VARIABLES_V2 = (
    "CMAKE_BUILD_TYPE", "CMAKE_GENERATOR", "CMAKE_CONFIGURATION_TYPES",
    "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_DEBUG",
    "CMAKE_CXX_FLAGS_RELEASE", "CMAKE_CXX_FLAGS_RELWITHDEBINFO",
    "CMAKE_CXX_FLAGS_MINSIZEREL", "ENABLE_OPENMP", "LEOPARD_ENABLE_GF8",
    "LEOPARD_ENABLE_GF16", "LEO2_BACKEND_VARIANT",
    "LEO2_BENCHMARK_GIT_EXECUTABLE", "LEO2_BUILD_BENCHMARKS",
    "LEO2_BUILD_TESTS", "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE",
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE",
)
BUILD_CONFIGURATION_VARIABLES_V3 = (
    *BUILD_CONFIGURATION_VARIABLES_V2[:-1],
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
    BUILD_CONFIGURATION_VARIABLES_V2[-1],
)
BUILD_CONFIGURATION_VARIABLES_V4 = (
    *BUILD_CONFIGURATION_VARIABLES_V2[:-1],
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR",
    "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING",
    "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING",
    "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING",
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT",
    "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE",
    BUILD_CONFIGURATION_VARIABLES_V2[-1],
)
BUILD_CONFIGURATION_VARIABLES_V5 = (
    *BUILD_CONFIGURATION_VARIABLES_V4[:-4],
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED",
    *BUILD_CONFIGURATION_VARIABLES_V4[-4:],
)
BUILD_CONFIGURATION_VARIABLES_V6 = (
    *BUILD_CONFIGURATION_VARIABLES_V5[:16],
    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED",
    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK",
    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL",
    *BUILD_CONFIGURATION_VARIABLES_V5[16:],
)
BUILD_CONFIGURATION_VARIABLES_V7 = (
    *BUILD_CONFIGURATION_VARIABLES_V6,
    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
)
BUILD_CONFIGURATION_VARIABLES_V8 = BUILD_CONFIGURATION_VARIABLES_V7
BUILD_CONFIGURATION_VARIABLES = (
    *BUILD_CONFIGURATION_VARIABLES_V8,
    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK",
)
REQUIRED_BUILD_CONFIGURATION_ENTRIES_V2 = {
    "CMAKE_CXX_FLAGS": " -Wall -Wextra -fopenmp",
    "CMAKE_CXX_FLAGS_DEBUG": "-g -O0",
    "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG -O3",
    "CMAKE_CXX_FLAGS_RELWITHDEBINFO": "-O2 -g -DNDEBUG",
    "CMAKE_CXX_FLAGS_MINSIZEREL": "-Os -DNDEBUG",
    "ENABLE_OPENMP": "ON",
    "LEOPARD_ENABLE_GF8": "ON",
    "LEOPARD_ENABLE_GF16": "ON",
    "LEO2_BACKEND_VARIANT": "auto",
    "LEO2_BENCHMARK_GIT_EXECUTABLE": "/usr/bin/git",
    "LEO2_BUILD_BENCHMARKS": "ON",
    "LEO2_BUILD_TESTS": "OFF",
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
}
REQUIRED_BUILD_CONFIGURATION_ENTRIES_V3 = {
    **REQUIRED_BUILD_CONFIGURATION_ENTRIES_V2,
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
}
REQUIRED_BUILD_CONFIGURATION_ENTRIES_V4 = {
    **REQUIRED_BUILD_CONFIGURATION_ENTRIES_V2,
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
    "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": "ON",
    "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
    "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": "ON",
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
    "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
}
REQUIRED_BUILD_CONFIGURATION_ENTRIES_V5 = {
    **REQUIRED_BUILD_CONFIGURATION_ENTRIES_V4,
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
}
REQUIRED_BUILD_CONFIGURATION_ENTRIES_V6 = {
    **REQUIRED_BUILD_CONFIGURATION_ENTRIES_V5,
    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
}
REQUIRED_BUILD_CONFIGURATION_ENTRIES_V7 = {
    **REQUIRED_BUILD_CONFIGURATION_ENTRIES_V6,
    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": "ON",
}
REQUIRED_BUILD_CONFIGURATION_ENTRIES_V8 = \
    REQUIRED_BUILD_CONFIGURATION_ENTRIES_V7
REQUIRED_BUILD_CONFIGURATION_ENTRIES_V9 = {
    **REQUIRED_BUILD_CONFIGURATION_ENTRIES_V8,
    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": "OFF",
}
REQUIRED_BUILD_CONFIGURATION_ENTRIES = {
    **REQUIRED_BUILD_CONFIGURATION_ENTRIES_V9,
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
}
REQUIRED_LEGACY_CANDIDATE_CACHE = {
    "ENABLE_OPENMP": "ON",
    "LEO2_BACKEND_VARIANT": "auto",
    "LEO2_BUILD_BENCHMARKS": "ON",
    "LEO2_BUILD_FUZZERS": "OFF",
    "LEO2_BUILD_TESTS": "OFF",
    "LEO2_ENABLE_CUDA": "OFF",
}
REQUIRED_CANDIDATE_CACHE_V2 = {
    **REQUIRED_LEGACY_CANDIDATE_CACHE,
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
}
REQUIRED_CANDIDATE_CACHE_V3 = {
    **REQUIRED_CANDIDATE_CACHE_V2,
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
}
REQUIRED_CANDIDATE_CACHE_V4 = {
    **REQUIRED_CANDIDATE_CACHE_V2,
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
    "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": "ON",
    "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
    "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": "ON",
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
    "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
    "LEO2_FLAG_MSSSE3": "1",
    "LEO2_FLAG_MNO_AVX": "1",
    "LEO2_FLAG_MAVX2": "1",
    "LEO2_FLAG_MNO_AVX512F": "1",
    "LEO2_FLAG_FALIGN_FUNCTIONS_64": "1",
    "LEO2_FLAG_MGFNI": "1",
    "LEO2_FLAG_MAVX512F": "FALSE",
    "LEO2_FLAG_MAVX512BW": "FALSE",
    "LEO2_FLAG_MAVX512VL": "FALSE",
    "LEO2_FLAG_MPREFER_VECTOR_WIDTH_256": "FALSE",
}
REQUIRED_CANDIDATE_CACHE_V5 = {
    **REQUIRED_CANDIDATE_CACHE_V4,
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
}
REQUIRED_CANDIDATE_CACHE_V6 = {
    **REQUIRED_CANDIDATE_CACHE_V5,
    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
}
REQUIRED_CANDIDATE_CACHE_V7 = {
    **REQUIRED_CANDIDATE_CACHE_V6,
    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": "ON",
}
REQUIRED_CANDIDATE_CACHE_V8 = REQUIRED_CANDIDATE_CACHE_V7
REQUIRED_CANDIDATE_CACHE_V9 = {
    **REQUIRED_CANDIDATE_CACHE_V8,
    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": "OFF",
}
REQUIRED_CANDIDATE_CACHE = {
    **REQUIRED_CANDIDATE_CACHE_V9,
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
}
CANONICAL_NINJA_PATH = "/usr/bin/ninja"
REQUIRED_BASELINE_CACHE = {
    "LEO_MAIN_HAS_MARCH_NATIVE": "1",
}
REQUIRED_PURE_AVX2_BASELINE_CACHE = {
    **REQUIRED_BASELINE_CACHE,
    "LEO_MAIN_PURE_AVX2": "ON",
    "LEO_MAIN_HAS_MARCH_X86_64": "1",
    "LEO_MAIN_HAS_MTUNE_GENERIC": "1",
    "LEO_MAIN_HAS_MAVX2": "1",
    "LEO_MAIN_HAS_MNO_AVX512F": "1",
}
EXCLUDED_CAMPAIGN_BACKENDS = {
    "avx512": "no forced AVX-512 runner coverage in this campaign schema",
    "neon": "no native NEON forced-runner coverage in this campaign schema",
}
# Exact string pairs emitted by DecodePathName/DecodePathRuleName in
# Leopard2Dispatch.h.  This list is deliberately complete: it prevents a
# negative predicate ("not generic") from accepting unknown or cross-paired
# strings while the exact AUTO prediction below narrows each attestation cell
# to the one pair its plan and selector geometry can produce.
PRODUCTION_DECODE_PATH_RULE_PAIRS = frozenset({
    ("no_op", "no_op"),
    ("direct", "direct"),
    ("generic", "forced_generic"),
    ("materialized", "forced_materialized"),
    ("tiled", "forced_tiled"),
    ("generic", "balanced_generic"),
    ("tiled", "measured_batch_tiled"),
    ("materialized", "measured_materialized"),
    ("tiled", "workspace_tiled"),
    ("materialized", "workspace_materialized"),
    ("materialized", "unsupported_profile"),
})
# Production direct-plan limits mirrored from leopard2.cpp.  The attestation
# source/build identity binds that file, so changing these values requires an
# explicit evidence-protocol update rather than silently broadening acceptance.
DIRECT_MAX_ORIGINALS = 16
DIRECT_MAX_LOSSES = 4


class PlanError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PlanError(message)


def finite_real(value: object, label: str) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{label} is not numeric")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise PlanError(f"{label} is outside the finite float range") from error
    require(math.isfinite(result), f"{label} is not finite")
    return result


def finite_number(value: object, label: str, *, positive: bool = False) -> float:
    result = finite_real(value, label)
    require((result > 0 if positive else result >= 0),
            f"{label} is not finite and "
            f"{'positive' if positive else 'nonnegative'}")
    return result


def canonical_bytes(value: object) -> bytes:
    try:
        rendered = json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        )
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False,
            ensure_ascii=False).encode("utf-8", errors="strict")
        return rendered.encode("utf-8")
    except (TypeError, ValueError, UnicodeEncodeError) as error:
        raise PlanError(f"value is not canonical JSON: {error}") from error


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(
        stable_regular_file_bytes(path, MAX_JSON_DOCUMENT_BYTES, str(path))
    ).hexdigest()


def _stable_stat_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_mode, value.st_nlink, value.st_size, value.st_mtime_ns,
        value.st_ctime_ns, value.st_dev, value.st_ino,
    )


def stable_regular_file_bytes(path: Path, maximum: int, label: str) -> bytes:
    """Read one bounded, single-link regular file without following symlinks."""
    require(type(maximum) is int and maximum > 0,
            f"{label} byte bound is invalid")
    try:
        path = path.absolute()
        parent = path.parent
        require(parent.resolve(strict=True) == parent and
                stat.S_ISDIR(os.lstat(parent).st_mode),
                f"{label} parent traverses a symbolic link")
        before = os.lstat(path)
        require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
                0 < before.st_size <= maximum,
                f"{label} is not a bounded single-link regular file")
        descriptor = os.open(
            path, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))
        try:
            initial = os.fstat(descriptor)
            path_initial = os.lstat(path)
            require(_stable_stat_identity(initial) ==
                    _stable_stat_identity(before) ==
                    _stable_stat_identity(path_initial),
                    f"{label} changed before its identity read")
            blocks: list[bytes] = []
            retained = 0
            while True:
                block = os.read(descriptor, min(1024 * 1024, maximum + 1 - retained))
                if not block:
                    break
                blocks.append(block)
                retained += len(block)
                require(retained <= maximum, f"{label} exceeds its byte bound")
            final = os.fstat(descriptor)
            path_final = os.lstat(path)
            require(_stable_stat_identity(initial) ==
                    _stable_stat_identity(final) ==
                    _stable_stat_identity(path_final) and
                    retained == final.st_size,
                    f"{label} changed during its identity read")
            return b"".join(blocks)
        finally:
            os.close(descriptor)
    except (OSError, ValueError) as error:
        if isinstance(error, PlanError):
            raise
        raise PlanError(f"cannot snapshot {label}: {error}") from error


def deterministic_seed(namespace: str, *values: int) -> int:
    seed = int.from_bytes(hashlib.sha256(canonical_bytes({
        "namespace": namespace, "values": values,
    })).digest()[:8], "big")
    return seed or 1


def signed(value: dict[str, Any]) -> dict[str, Any]:
    result = dict(value)
    require("content_sha256" not in result, "document is already signed")
    result["content_sha256"] = canonical_sha256(result)
    return result


def unsigned(value: object, schema: str, label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    result = dict(value)
    digest = result.pop("content_sha256", None)
    require(result.get("schema") == schema, f"{label} schema differs")
    require(isinstance(digest, str) and len(digest) == 64 and
            canonical_sha256(result) == digest, f"{label} digest differs")
    return result


def _open_directory_chain(path: Path, *, create: bool) -> tuple[int, Path]:
    """Open one absolute directory component-by-component without symlinks."""
    require(os.name == "posix",
            "secure balanced-evidence publication requires POSIX dirfds")
    absolute = Path(os.path.abspath(os.fspath(path)))
    require(absolute.is_absolute() and absolute.name not in {".", ".."},
            f"invalid JSON output directory {path}")
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_DIRECTORY", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    descriptor = -1
    try:
        descriptor = os.open(absolute.anchor, flags)
        for component in absolute.parts[1:]:
            require(component not in {"", ".", ".."},
                    f"invalid JSON output directory component {component!r}")
            created = False
            if create:
                try:
                    os.mkdir(component, 0o700, dir_fd=descriptor)
                    created = True
                except FileExistsError:
                    pass
            child = os.open(component, flags, dir_fd=descriptor)
            try:
                identity = os.fstat(child)
                require(stat.S_ISDIR(identity.st_mode),
                        f"JSON output component is not a directory: "
                        f"{absolute}")
                if created:
                    os.fchmod(child, 0o700)
            except BaseException:
                os.close(child)
                raise
            os.close(descriptor)
            descriptor = child
        return descriptor, absolute
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
        raise


def _directory_fd_matches_path(descriptor: int, path: Path) -> bool:
    """Confirm that a lexical no-symlink path still names a retained dirfd."""
    reopened = -1
    try:
        reopened, _ = _open_directory_chain(path, create=False)
        held = os.fstat(descriptor)
        current = os.fstat(reopened)
        return (
            stat.S_ISDIR(held.st_mode) and stat.S_ISDIR(current.st_mode) and
            (held.st_dev, held.st_ino) == (current.st_dev, current.st_ino)
        )
    except (OSError, PlanError):
        return False
    finally:
        if reopened >= 0:
            os.close(reopened)


def _write_fd_all(descriptor: int, payload: bytes) -> None:
    offset = 0
    while offset < len(payload):
        written = os.write(descriptor, payload[offset:])
        require(written > 0, "JSON output write made no progress")
        offset += written


def _read_fd_all(descriptor: int, maximum: int) -> bytes:
    os.lseek(descriptor, 0, os.SEEK_SET)
    blocks: list[bytes] = []
    retained = 0
    while True:
        block = os.read(descriptor, min(1024 * 1024, maximum + 1 - retained))
        if not block:
            break
        blocks.append(block)
        retained += len(block)
        require(retained <= maximum, "published JSON exceeds its byte bound")
    return b"".join(blocks)


def write_json(path: Path, value: object) -> None:
    try:
        payload = (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8", errors="strict")
    except (TypeError, ValueError, UnicodeEncodeError) as error:
        raise PlanError(f"cannot encode canonical JSON for {path}: {error}") from error
    absolute = Path(os.path.abspath(os.fspath(path)))
    require(absolute.name not in {"", ".", ".."},
            f"invalid JSON output name {path}")
    parent_descriptor = -1
    temporary_descriptor = -1
    final_descriptor = -1
    temporary_name: str | None = None
    final_created = False
    completed = False
    primary_error: BaseException | None = None
    cleanup_errors: list[BaseException] = []
    try:
        parent_descriptor, parent = _open_directory_chain(
            absolute.parent, create=True)
        require(_directory_fd_matches_path(parent_descriptor, parent),
                f"JSON output parent changed before publication: {parent}")
        for _ in range(128):
            candidate = (
                f".{absolute.name}.{os.getpid()}."
                f"{os.urandom(16).hex()}.tmp"
            )
            try:
                temporary_descriptor = os.open(
                    candidate,
                    os.O_RDWR | os.O_CREAT | os.O_EXCL |
                    getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_NOFOLLOW", 0),
                    0o600, dir_fd=parent_descriptor)
                temporary_name = candidate
                break
            except FileExistsError:
                continue
        require(temporary_descriptor >= 0 and temporary_name is not None,
                f"cannot allocate temporary JSON output for {absolute}")
        os.fchmod(temporary_descriptor, 0o600)
        _write_fd_all(temporary_descriptor, payload)
        os.fsync(temporary_descriptor)
        staged = os.fstat(temporary_descriptor)
        require(stat.S_ISREG(staged.st_mode) and
                stat.S_IMODE(staged.st_mode) == 0o600 and
                staged.st_nlink == 1 and staged.st_size == len(payload),
                f"temporary JSON output has unsafe identity: {absolute}")
        require(_directory_fd_matches_path(parent_descriptor, parent),
                f"JSON output parent changed before commit: {parent}")
        try:
            os.link(
                temporary_name, absolute.name,
                src_dir_fd=parent_descriptor, dst_dir_fd=parent_descriptor,
                follow_symlinks=False)
        except FileExistsError as error:
            raise PlanError(
                f"refusing to replace existing output {absolute}") from error
        final_created = True
        final_descriptor = os.open(
            absolute.name,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0) |
            getattr(os, "O_NONBLOCK", 0),
            dir_fd=parent_descriptor)
        published_before = os.fstat(final_descriptor)
        require(
            stat.S_ISREG(published_before.st_mode) and
            stat.S_IMODE(published_before.st_mode) == 0o600 and
            published_before.st_nlink == 2 and
            published_before.st_size == len(payload) and
            (published_before.st_dev, published_before.st_ino) ==
                (staged.st_dev, staged.st_ino),
            f"published JSON output has unsafe identity: {absolute}")
        require(_read_fd_all(final_descriptor, len(payload)) == payload,
                f"published JSON bytes differ: {absolute}")
        published_after = os.fstat(final_descriptor)
        require(_stable_stat_identity(published_before) ==
                _stable_stat_identity(published_after),
                f"published JSON changed during readback: {absolute}")
        os.unlink(temporary_name, dir_fd=parent_descriptor)
        temporary_name = None
        os.fsync(parent_descriptor)
        linked = os.fstat(final_descriptor)
        require(linked.st_nlink == 1 and
                _read_fd_all(final_descriptor, len(payload)) == payload,
                f"published JSON changed after commit: {absolute}")
        require(_directory_fd_matches_path(parent_descriptor, parent),
                f"JSON output parent changed during publication: {parent}")
        completed = True
    except BaseException as error:
        if isinstance(error, (OSError, ValueError)) and \
                not isinstance(error, PlanError):
            primary_error = PlanError(
                f"cannot publish JSON output {absolute}: {error}")
        else:
            primary_error = error

    if final_descriptor >= 0:
        try:
            os.close(final_descriptor)
        except BaseException as error:
            cleanup_errors.append(error)
    if temporary_descriptor >= 0:
        try:
            os.close(temporary_descriptor)
        except BaseException as error:
            cleanup_errors.append(error)
    if temporary_name is not None and parent_descriptor >= 0:
        try:
            os.unlink(temporary_name, dir_fd=parent_descriptor)
        except FileNotFoundError:
            pass
        except BaseException as error:
            cleanup_errors.append(error)
    if final_created and not completed and parent_descriptor >= 0:
        try:
            os.unlink(absolute.name, dir_fd=parent_descriptor)
            os.fsync(parent_descriptor)
        except FileNotFoundError:
            pass
        except BaseException as error:
            cleanup_errors.append(error)
    if parent_descriptor >= 0:
        try:
            os.close(parent_descriptor)
        except BaseException as error:
            cleanup_errors.append(error)

    if primary_error is not None:
        if cleanup_errors:
            details = "; ".join(str(error) for error in cleanup_errors)
            raise PlanError(
                f"{primary_error}; JSON publication cleanup failed: {details}"
            ) from primary_error
        raise primary_error
    if cleanup_errors:
        details = "; ".join(str(error) for error in cleanup_errors)
        raise PlanError(
            f"JSON publication cleanup failed for {absolute}: {details}")


def decode_json_bytes(payload: bytes, label: str) -> object:
    def reject_duplicate_keys(
            pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise PlanError(
                    f"{label} contains duplicate JSON key {key!r}")
            result[key] = value
        return result

    def reject_constant(value: str) -> object:
        raise PlanError(
            f"{label} contains non-standard JSON constant {value!r}")

    def finite_float(value: str) -> float:
        result = float(value)
        if not math.isfinite(result):
            raise PlanError(
                f"{label} contains non-finite JSON number {value!r}")
        return result

    try:
        text = payload.decode("utf-8", errors="strict")
        return json.loads(
            text, object_pairs_hook=reject_duplicate_keys,
            parse_constant=reject_constant, parse_float=finite_float)
    except (OSError, UnicodeDecodeError, ValueError, OverflowError,
            RecursionError) as error:
        raise PlanError(f"cannot decode strict JSON {label}: {error}") from error


def load_json(path: Path) -> object:
    payload = stable_regular_file_bytes(
        path, MAX_JSON_DOCUMENT_BYTES, f"JSON document {path}")
    return decode_json_bytes(payload, str(path))


def main_cell(identifier: str, k: int, r: int, shard_bytes: int,
              losses: int, namespace: str) -> dict[str, Any]:
    require(re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", identifier) is not None,
            f"invalid exact-main identifier {identifier!r}")
    require(0 < r <= k and 0 < losses <= r and
            shard_bytes > 0 and shard_bytes % 64 == 0,
            f"invalid exact-main cell {identifier}")
    seed = deterministic_seed(namespace, k, r, shard_bytes, losses)
    return {
        "identifier": identifier, "K": k, "R": r,
        "shard_bytes": shard_bytes, "loss_count": losses, "seed": seed,
        "runner_cell": f"{identifier}:{k}:{r}:{shard_bytes}:{losses}:{seed}",
    }


def validate_main_cell(cell: object, label: str) -> dict[str, Any]:
    require(isinstance(cell, dict) and set(cell) == {
        "identifier", "K", "R", "shard_bytes", "loss_count", "seed",
        "runner_cell",
    }, f"{label} shape differs")
    identifier = cell["identifier"]
    require(isinstance(identifier, str) and
            re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", identifier) is not None,
            f"{label} identifier is invalid")
    for key in ("K", "R", "shard_bytes", "loss_count", "seed"):
        require(type(cell[key]) is int and cell[key] > 0,
                f"{label} {key} is invalid")
    require(cell["R"] <= cell["K"] and cell["loss_count"] <= cell["R"] and
            cell["shard_bytes"] % 64 == 0 and cell["seed"] < 1 << 64,
            f"{label} is outside exact-main limits")
    require(cell["runner_cell"] == (
        f"{identifier}:{cell['K']}:{cell['R']}:{cell['shard_bytes']}:"
        f"{cell['loss_count']}:{cell['seed']}"),
        f"{label} runner encoding differs")
    return cell


def main_document(name: str, mode: str, cells: list[dict[str, Any]],
                  purpose: str, run_condition: str) -> dict[str, Any]:
    require(mode in ("auto", "generic", "materialized", "tiled"),
            "invalid exact-main mode")
    return signed({
        "schema": MAIN_CELL_SCHEMA, "name": name, "candidate_mode": mode,
        "exact_main_commit": EXACT_MAIN_COMMIT, "reuse": 8,
        "iterations": 9, "warmup": 2, "purpose": purpose,
        "run_condition": run_condition, "cells": cells,
    })


def validate_main_document(value: object, label: str) -> dict[str, Any]:
    result = unsigned(value, MAIN_CELL_SCHEMA, label)
    require(set(result) == {
        "schema", "name", "candidate_mode", "exact_main_commit", "reuse",
        "iterations", "warmup", "purpose", "run_condition", "cells",
    }, f"{label} fields differ")
    require(type(result["candidate_mode"]) is str and
            result["candidate_mode"] in CONFIRM_MODES and
            result["exact_main_commit"] == EXACT_MAIN_COMMIT and
            type(result["reuse"]) is int and result["reuse"] == 8 and
            type(result["iterations"]) is int and result["iterations"] == 9 and
            type(result["warmup"]) is int and result["warmup"] == 2 and
            isinstance(result["name"], str) and bool(result["name"]) and
            isinstance(result["purpose"], str) and bool(result["purpose"]) and
            isinstance(result["run_condition"], str) and
            bool(result["run_condition"]),
            f"{label} semantics differ")
    cells = result["cells"]
    require(isinstance(cells, list) and cells, f"{label} has no cells")
    for index, cell in enumerate(cells):
        validate_main_cell(cell, f"{label} cell {index}")
    require(len({cell["identifier"] for cell in cells}) == len(cells),
            f"{label} identifiers are duplicated")
    return result


def transform_groups(counts: Iterable[int]) -> dict[int, list[int]]:
    result: dict[int, list[int]] = {}
    for k in sorted(set(counts)):
        result.setdefault(1 << (k - 1).bit_length(), []).append(k)
    return result


def artifact(path: Path, root: Path, case_count: int, child_count: int,
             kind: str) -> dict[str, Any]:
    snapshot = stable_regular_file_bytes(
        path, MAX_JSON_DOCUMENT_BYTES, f"generated artifact {path}")
    return {
        "path": str(path.relative_to(root)), "kind": kind,
        "case_count": case_count, "timed_child_count": child_count,
        "size": len(snapshot), "sha256": hashlib.sha256(snapshot).hexdigest(),
    }


def generate_plan(root: Path) -> dict[str, Any]:
    require(not root.exists(), f"refusing to replace existing output {root}")
    root.mkdir(parents=True)
    artifacts = []
    for side, counts in sorted(transform_groups(BOUNDARY_K).items()):
        cells = [
            main_cell(f"gate-k{k}-b{shard_bytes}", k, k, shard_bytes, k,
                      "exact-main-generic-gate")
            for k in counts for shard_bytes in EXACT_GATE_BYTES
        ]
        path = root / "exact-main-gate" / f"t{side}.json"
        write_json(path, main_document(
            f"balanced-generic-exact-main-gate-t{side}", "generic", cells,
            "External gate for balanced full-loss generic decode.",
            "Run before any internal forced-path promotion comparison."))
        artifacts.append(artifact(path, root, len(cells), len(cells) * 12,
                                  "exact_main_gate"))
    plan = signed({
        "schema": PLAN_SCHEMA,
        "name": "balanced-current-source-exact-main-first",
        "exact_main_commit": EXACT_MAIN_COMMIT,
        "promotion_gate": {
            "metrics": ["decode_first_use", "decode_reuse_amortized"],
            "minimum_ci95_lower": PROMOTION_GATE,
            "neighbor_maximum_regression": NEIGHBOR_MAXIMUM_REGRESSION,
        },
        "runner_ordering": {
            "exact_main": {
                "round_count": 3, "order_in_every_round": list(EXACT_ORDER),
                "labels": ["ABBA", "ABBA", "ABBA"],
            },
            "forced_same_binary": {
                "round_count": 3,
                "orders": [list(order) for order in FORCED_ORDERS],
                "labels": ["ABBA", "BAAB", "ABBA"],
            },
        },
        "workflow": {
            "select_command": (
                "select verifies every complete-schema exact-main gate manifest and "
                "derives an evidence-bound survivor set"),
            "refinement": (
                "opposite outcomes separated by an unmeasured K emit exact "
                "runner_cell refinements; advance refuses until none remain"),
            "advance_command": (
                "advance emits only surviving K/byte forced comparisons, "
                "aligned exact-main confirmations, true tails, and rejection controls"),
            "promotion": (
                "AUTO selected-path attestation and exact-main confirmation are "
                "required before a dispatcher change"),
        },
        "artifacts": artifacts,
        "initial_case_count": sum(item["case_count"] for item in artifacts),
        "initial_timed_child_count": sum(
            item["timed_child_count"] for item in artifacts),
    })
    write_json(root / "plan.json", plan)
    return plan


def json_tree(root: Path) -> list[Path]:
    return sorted(path.relative_to(root) for path in root.rglob("*.json"))


def compare_trees(actual: Path, expected: Path, label: str) -> None:
    actual_files = json_tree(actual)
    expected_files = json_tree(expected)
    require(actual_files == expected_files, f"{label} file set differs")
    for relative in expected_files:
        require(stable_regular_file_bytes(
                    actual / relative, MAX_JSON_DOCUMENT_BYTES,
                    f"{label} actual {relative}") ==
                stable_regular_file_bytes(
                    expected / relative, MAX_JSON_DOCUMENT_BYTES,
                    f"{label} expected {relative}"),
                f"{label} canonical bytes differ: {relative}")


def structural_plan(root: Path) -> dict[str, Any]:
    plan = unsigned(load_json(root / "plan.json"), PLAN_SCHEMA, "plan")
    require(set(plan) == {
        "schema", "name", "exact_main_commit", "promotion_gate",
        "runner_ordering", "workflow", "artifacts", "initial_case_count",
        "initial_timed_child_count",
    }, "plan fields differ")
    require(plan["exact_main_commit"] == EXACT_MAIN_COMMIT and
            plan["promotion_gate"] == {
                "metrics": ["decode_first_use", "decode_reuse_amortized"],
                "minimum_ci95_lower": PROMOTION_GATE,
                "neighbor_maximum_regression": NEIGHBOR_MAXIMUM_REGRESSION,
            }, "plan gate semantics differ")
    total_cases = 0
    total_children = 0
    for item in plan["artifacts"]:
        require(isinstance(item, dict) and set(item) == {
            "path", "kind", "case_count", "timed_child_count", "size", "sha256",
        } and item["kind"] == "exact_main_gate", "plan artifact differs")
        path = root / item["path"]
        snapshot = stable_regular_file_bytes(
            path, MAX_JSON_DOCUMENT_BYTES, f"plan artifact {item['path']}")
        require(len(snapshot) == item["size"] and
                hashlib.sha256(snapshot).hexdigest() == item["sha256"],
                "plan artifact bytes differ")
        document = validate_main_document(
            decode_json_bytes(snapshot, str(path)), str(path))
        require(document["candidate_mode"] == "generic" and
                len(document["cells"]) == item["case_count"] and
                item["timed_child_count"] == item["case_count"] * 12,
                "gate artifact accounting differs")
        total_cases += item["case_count"]
        total_children += item["timed_child_count"]
    require(total_cases == plan["initial_case_count"] and
            total_children == plan["initial_timed_child_count"],
            "plan totals differ")
    return plan


def validate_plan(root: Path) -> dict[str, Any]:
    plan = structural_plan(root)
    with tempfile.TemporaryDirectory(prefix="leopard2-plan-expected-") as temporary:
        expected = Path(temporary) / "plan"
        generate_plan(expected)
        compare_trees(root, expected, "plan")
    return plan


def planned_gate_cells(root: Path, plan: dict[str, Any]) -> dict[tuple[int, int], dict[str, Any]]:
    result = {}
    for item in plan["artifacts"]:
        document = validate_main_document(load_json(root / item["path"]), item["path"])
        for cell in document["cells"]:
            key = (cell["K"], cell["shard_bytes"])
            require(key not in result, "planned gate cell is duplicated")
            result[key] = cell
    return result


def _normalize_bound_paths(value: object, replacements: tuple[tuple[str, str], ...]) -> object:
    """Remove volatile filesystem metadata while retaining exact byte identity."""
    require(replacements and all(isinstance(original, str) and original and
                                 isinstance(marker, str) and marker
                                 for original, marker in replacements) and
            all(Path(original).is_absolute() and original != "/" and
                os.path.normpath(original) == original
                for original, _ in replacements) and
            len({original for original, _ in replacements}) == len(replacements) and
            len({marker for _, marker in replacements}) == len(replacements),
            "scope path replacements are invalid or ambiguous")
    # A baseline build may live below the candidate source tree.  Replacing the
    # longest roots first prevents a broad source marker from swallowing the
    # more specific, role-distinct build identity.
    ordered = tuple(sorted(replacements, key=lambda item: len(item[0]), reverse=True))

    token_prefix_delimiters = frozenset(" \t\r\n\"'`=(:,;[{")
    token_suffix_delimiters = frozenset(" \t\r\n\"'`,:;)]}")
    cmake_external_object_suffix = re.compile(
        r"CMakeFiles/[^/\s\"'`]+\.dir(?:/[^/\s\"'`]+)?$")

    def replace_root_tokens(text: str, original: str, marker: str) -> str:
        """Replace an absolute root only at a standalone path-token boundary."""
        pieces: list[str] = []
        cursor = 0
        while True:
            index = text.find(original, cursor)
            if index < 0:
                pieces.append(text[cursor:])
                return "".join(pieces)
            end = index + len(original)
            prefix = text[:index]
            cmake_match = cmake_external_object_suffix.search(prefix)
            cmake_external_source = cmake_match is not None and (
                cmake_match.start() == 0 or
                prefix[cmake_match.start() - 1] == "/" or
                prefix[cmake_match.start() - 1] in token_prefix_delimiters)
            fused_include = index >= 2 and text[index - 2:index] == "-I" and (
                index == 2 or text[index - 3] in token_prefix_delimiters)
            before_ok = index == 0 or \
                text[index - 1] in token_prefix_delimiters or \
                cmake_external_source or fused_include
            after_ok = end == len(text) or text[end] == "/" or \
                text[end] in token_suffix_delimiters
            pieces.append(text[cursor:index])
            if before_ok and after_ok:
                # CMake encodes an absolute external source below
                # CMakeFiles/<target>.dir by stripping no leading slash:
                # target.dir/home/user/source.cpp.o.  Preserve a separator
                # before the abstract root marker in that one typed context.
                pieces.append(("/" if cmake_external_source else "") + marker)
            else:
                pieces.append(original)
            cursor = end

    def visit(item: object) -> object:
        if isinstance(item, dict):
            result = {
                key: visit(child) for key, child in sorted(item.items())
                if key not in {"device", "inode", "mtime_ns"}
            }
            # Retained recipe text legitimately contains source/build roots.
            # Replacing those roots changes the normalized bytes, so recompute
            # their normalized content identity instead of leaving a stale
            # size/hash pair or retaining machine-specific absolute paths.
            if set(result) == {"encoding", "size", "sha256", "text"} and \
                    result.get("encoding") == "utf-8" and \
                    isinstance(result.get("text"), str):
                encoded = result["text"].encode("utf-8")
                result["size"] = len(encoded)
                result["sha256"] = hashlib.sha256(encoded).hexdigest()
            ninja_graph = result.get("multi_config_ninja_graph")
            if ninja_graph is None:
                for prefix in ("archive", "executable"):
                    content = result.get(f"{prefix}_link_recipe_content")
                    identity = result.get(f"{prefix}_link_recipe")
                    if isinstance(content, dict) and \
                            isinstance(identity, dict) and \
                            type(content.get("size")) is int and \
                            isinstance(content.get("sha256"), str):
                        identity["size"] = content["size"]
                        identity["sha256"] = content["sha256"]
            elif isinstance(ninja_graph, dict):
                entrypoint = ninja_graph.get("entrypoint")
                match = (re.fullmatch(r"build-([A-Za-z0-9_.+-]+)\.ninja",
                                      entrypoint)
                         if isinstance(entrypoint, str) else None)
                files = ninja_graph.get("files")
                require(match is not None and isinstance(files, list),
                        "normalized multi-config Ninja graph is malformed")
                link_relative = \
                    f"CMakeFiles/impl-{match.group(1)}.ninja"
                link_records = [
                    record for record in files
                    if isinstance(record, dict) and
                    record.get("relative_path") == link_relative
                ]
                require(len(link_records) == 1 and
                        isinstance(link_records[0].get("artifact"), dict),
                        "normalized multi-config Ninja link graph is absent")
                graph_identity = link_records[0]["artifact"]
                for prefix in ("archive", "executable"):
                    identity = result.get(f"{prefix}_link_recipe")
                    require(isinstance(identity, dict),
                            f"normalized {prefix} link identity is absent")
                    for field in ("path", "size", "mode", "sha256"):
                        require(field in graph_identity,
                                "normalized multi-config Ninja link graph "
                                f"omits {field}")
                        identity[field] = graph_identity[field]
            if set(result) == {"relative_path", "artifact", "content"}:
                content = result.get("content")
                artifact = result.get("artifact")
                if isinstance(content, dict) and isinstance(artifact, dict) and \
                        type(content.get("size")) is int and \
                        isinstance(content.get("sha256"), str):
                    artifact["size"] = content["size"]
                    artifact["sha256"] = content["sha256"]
            return result
        if isinstance(item, list):
            return [visit(child) for child in item]
        if isinstance(item, str):
            result = item
            for original, marker in ordered:
                result = replace_root_tokens(result, original, marker)
            return result
        return item

    return visit(value)


def _parse_scope_cpu_list(value: object, label: str) -> set[int]:
    require(isinstance(value, str) and value and
            len(value.encode("utf-8")) <= 65_536,
            f"{label} is not a bounded CPU-list string")
    result: set[int] = set()
    for component in value.split(","):
        component = component.strip()
        require(re.fullmatch(r"[0-9]+(?:-[0-9]+)?", component) is not None,
                f"{label} is malformed")
        bounds = [int(item) for item in component.split("-", 1)]
        first, last = bounds[0], bounds[-1]
        require(0 <= first <= last <= MAX_CPU_ID,
                f"{label} range is reversed or outside the CPU-ID bound")
        result.update(range(first, last + 1))
        require(len(result) <= MAX_CPU_LIST_ENTRIES,
                f"{label} expands beyond the CPU-count bound")
    require(result, f"{label} is empty")
    return result


def _validate_scope_artifact(
    value: object, label: str, expected_kind: str | None = None,
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "path", "kind", "size", "mode", "sha256"},
            f"{label} normalized artifact shape differs")
    require(isinstance(value.get("path"), str) and value["path"] and
            isinstance(value.get("kind"), str) and value["kind"] and
            (expected_kind is None or value["kind"] == expected_kind) and
            type(value.get("size")) is int and value["size"] >= 0 and
            type(value.get("mode")) is int and 0 <= value["mode"] <= 0o7777 and
            isinstance(value.get("sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is not None,
            f"{label} normalized artifact is invalid")
    if expected_kind == "build_tool":
        require(value["mode"] & 0o111,
                f"{label} normalized artifact is not executable")
    return value


def _validate_scope_build_tool_version(
    value: object, label: str,
) -> dict[str, str]:
    require(isinstance(value, dict) and set(value) == {"sha256", "text"} and
            isinstance(value.get("text"), str) and value["text"] and
            "\0" not in value["text"] and "\r" not in value["text"] and
            isinstance(value.get("sha256"), str),
            f"{label} shape differs")
    try:
        encoded = value["text"].encode("utf-8", errors="strict")
    except UnicodeEncodeError as error:
        raise PlanError(f"{label} is not strict UTF-8") from error
    require(len(encoded) <= MAX_BUILD_TOOL_VERSION_BYTES and
            value["sha256"] == hashlib.sha256(encoded).hexdigest(),
            f"{label} identity differs")
    return value


def _canonical_ninja_graph_relative_path(value: object, label: str) -> str:
    require(isinstance(value, str) and value and
            re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._/+@=-]*", value)
            is not None,
            f"{label} is not a literal Ninja graph path")
    path = Path(value)
    require(not path.is_absolute() and
            all(component not in {"", ".", ".."} for component in path.parts) and
            path.as_posix() == value,
            f"{label} is not a canonical relative Ninja graph path")
    return value


def _ninja_graph_includes(text: object, label: str) -> tuple[str, ...]:
    require(isinstance(text, str) and text.endswith("\n") and
            "\0" not in text and "\r" not in text,
            f"{label} is not canonical LF-terminated Ninja text")
    result: list[str] = []
    for line_number, line in enumerate(text[:-1].split("\n"), 1):
        match = re.match(r"^(include|subninja)(?:[ \t]|$)", line)
        if match is None:
            continue
        directive = re.fullmatch(
            r"(?:include|subninja)[ \t]+([^ \t]+)[ \t]*", line)
        require(directive is not None,
                f"{label}:{line_number} has a non-literal Ninja include")
        result.append(_canonical_ninja_graph_relative_path(
            directive.group(1), f"{label}:{line_number} include"))
    return tuple(result)


def _validate_scope_ninja_graph_closure(
    value: object, build_root: str, selected_configuration: str, label: str,
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "schema", "entrypoint", "files"} and
            value.get("schema") == NINJA_GRAPH_CLOSURE_SCHEMA and
            value.get("entrypoint") ==
                f"build-{selected_configuration}.ninja" and
            isinstance(value.get("files"), list) and
            0 < len(value["files"]) <= MAX_NINJA_GRAPH_FILES,
            f"{label} shape differs")
    records: dict[str, dict[str, Any]] = {}
    total_bytes = 0
    for record in value["files"]:
        require(isinstance(record, dict) and set(record) == {
                    "relative_path", "artifact", "content"},
                f"{label} record shape differs")
        relative = _canonical_ninja_graph_relative_path(
            record.get("relative_path"), f"{label} record path")
        require(relative not in records,
                f"{label} contains a duplicate graph path")
        artifact = _validate_scope_artifact(
            record.get("artifact"), f"{label} {relative}",
            "ninja_graph_input")
        content = _validate_scope_text(
            record.get("content"), f"{label} {relative} content")
        total_bytes += content["size"]
        require(
            artifact["path"] ==
                f"{build_root}/{relative}" and
            artifact["size"] == content["size"] and
            artifact["sha256"] == content["sha256"] and
            total_bytes <= MAX_NINJA_GRAPH_TOTAL_BYTES,
            f"{label} {relative} artifact/content identity differs")
        records[relative] = record
    require(list(records) == sorted(records),
            f"{label} records are not canonically ordered")
    pending = [value["entrypoint"]]
    visited: set[str] = set()
    while pending:
        relative = pending.pop()
        require(relative in records,
                f"{label} references an unretained graph input: {relative}")
        if relative in visited:
            continue
        visited.add(relative)
        pending.extend(_ninja_graph_includes(
            records[relative]["content"]["text"],
            f"{label} {relative} content"))
    require(visited == set(records),
            f"{label} contains inputs outside its entrypoint closure")
    return value


def _validate_scope_text(value: object, label: str) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "encoding", "size", "sha256", "text"} and
            value.get("encoding") == "utf-8" and
            isinstance(value.get("text"), str),
            f"{label} normalized retained text shape differs")
    try:
        encoded = value["text"].encode("utf-8", errors="strict")
    except UnicodeEncodeError as error:
        raise PlanError(f"{label} is not strict UTF-8 text") from error
    require(type(value.get("size")) is int and
            0 < value["size"] <= MAX_RETAINED_TEXT_BYTES and
            value["size"] == len(encoded) and "\x00" not in value["text"] and
            isinstance(value.get("sha256"), str) and
            value["sha256"] == hashlib.sha256(encoded).hexdigest(),
            f"{label} normalized retained text identity differs")
    return value


def _build_configuration_contract(
    file_schema: str,
) -> tuple[str, tuple[str, ...], Mapping[str, str]]:
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA:
        return (
            BUILD_CONFIGURATION_RECORD_SCHEMA,
            BUILD_CONFIGURATION_VARIABLES,
            REQUIRED_CANDIDATE_CACHE,
        )
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V8:
        return (
            BUILD_CONFIGURATION_RECORD_SCHEMA_V8,
            BUILD_CONFIGURATION_VARIABLES_V8,
            REQUIRED_CANDIDATE_CACHE_V8,
        )
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V7:
        return (
            BUILD_CONFIGURATION_RECORD_SCHEMA_V7,
            BUILD_CONFIGURATION_VARIABLES_V6,
            REQUIRED_CANDIDATE_CACHE_V6,
        )
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V6:
        return (
            BUILD_CONFIGURATION_RECORD_SCHEMA_V6,
            BUILD_CONFIGURATION_VARIABLES_V6,
            REQUIRED_CANDIDATE_CACHE_V6,
        )
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V5:
        return (
            BUILD_CONFIGURATION_RECORD_SCHEMA_V5,
            BUILD_CONFIGURATION_VARIABLES_V5,
            REQUIRED_CANDIDATE_CACHE_V5,
        )
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V3:
        return (
            BUILD_CONFIGURATION_RECORD_SCHEMA_V3,
            BUILD_CONFIGURATION_VARIABLES_V3,
            REQUIRED_CANDIDATE_CACHE_V3,
        )
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V4:
        return (
            BUILD_CONFIGURATION_RECORD_SCHEMA_V4,
            BUILD_CONFIGURATION_VARIABLES_V4,
            REQUIRED_CANDIDATE_CACHE_V4,
        )
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V2:
        return (
            BUILD_CONFIGURATION_RECORD_SCHEMA_V2,
            BUILD_CONFIGURATION_VARIABLES_V2,
            REQUIRED_CANDIDATE_CACHE_V2,
        )
    raise PlanError(
        "candidate normalized effective-configuration schema differs")


def _build_configuration_contract_for_compile_schema(
    compile_schema: str, file_schema: str,
) -> tuple[str, tuple[str, ...], Mapping[str, str]]:
    """Bind sidecar versions to the compile-command contract that names them."""
    require(compile_schema in {
                COMPILE_COMMANDS_SCHEMA_V2, COMPILE_COMMANDS_SCHEMA_V3,
                COMPILE_COMMANDS_SCHEMA_V4, COMPILE_COMMANDS_SCHEMA_V5,
                COMPILE_COMMANDS_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V7,
                COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
                COMPILE_COMMANDS_SCHEMA_V10,
                COMPILE_COMMANDS_SCHEMA_V11,
                COMPILE_COMMANDS_SCHEMA},
            "normalized compile-command schema differs")
    if compile_schema == COMPILE_COMMANDS_SCHEMA:
        require(file_schema == BUILD_CONFIGURATION_FILE_SCHEMA,
                "current compile-command schema requires the current "
                "effective-configuration schema")
    elif compile_schema == COMPILE_COMMANDS_SCHEMA_V11:
        require(file_schema == BUILD_CONFIGURATION_FILE_SCHEMA,
                "v15 compile-command schema requires the v15 "
                "effective-configuration schema")
        return (
            BUILD_CONFIGURATION_RECORD_SCHEMA,
            BUILD_CONFIGURATION_VARIABLES,
            REQUIRED_CANDIDATE_CACHE_V9,
        )
    elif compile_schema == COMPILE_COMMANDS_SCHEMA_V10:
        require(file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V8,
                "v14 compile-command schema requires the v14 "
                "effective-configuration schema")
    elif compile_schema == COMPILE_COMMANDS_SCHEMA_V9:
        require(file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V7,
                "v13 compile-command schema requires the v13 "
                "effective-configuration schema")
    elif compile_schema == COMPILE_COMMANDS_SCHEMA_V8:
        require(file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V6,
                "v12 compile-command schema requires the v12 "
                "effective-configuration schema")
    elif compile_schema == COMPILE_COMMANDS_SCHEMA_V7:
        require(file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V5,
                "v11 compile-command schema requires the v11 "
                "effective-configuration schema")
    elif compile_schema == COMPILE_COMMANDS_SCHEMA_V6:
        require(file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V4,
                "v10 compile-command schema requires the v10 "
                "effective-configuration schema")
    elif compile_schema == COMPILE_COMMANDS_SCHEMA_V5:
        require(file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V3,
                "v9 compile-command schema requires the v9 "
                "effective-configuration schema")
    else:
        require(file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V2,
                "historical compile-command schema requires the historical "
                "effective-configuration schema")
    return _build_configuration_contract(file_schema)


def _required_build_configuration_entries(
    file_schema: str, compile_schema: str = COMPILE_COMMANDS_SCHEMA,
) -> Mapping[str, str]:
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA:
        return (REQUIRED_BUILD_CONFIGURATION_ENTRIES_V9
                if compile_schema == COMPILE_COMMANDS_SCHEMA_V11 else
                REQUIRED_BUILD_CONFIGURATION_ENTRIES)
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V8:
        return REQUIRED_BUILD_CONFIGURATION_ENTRIES_V8
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V7:
        return REQUIRED_BUILD_CONFIGURATION_ENTRIES_V6
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V6:
        return REQUIRED_BUILD_CONFIGURATION_ENTRIES_V6
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V5:
        return REQUIRED_BUILD_CONFIGURATION_ENTRIES_V5
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V3:
        return REQUIRED_BUILD_CONFIGURATION_ENTRIES_V3
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V4:
        return REQUIRED_BUILD_CONFIGURATION_ENTRIES_V4
    if file_schema == BUILD_CONFIGURATION_FILE_SCHEMA_V2:
        return REQUIRED_BUILD_CONFIGURATION_ENTRIES_V2
    raise PlanError(
        "candidate normalized effective-configuration schema differs")


def _build_configuration_material(
    entries: Mapping[str, str],
    variables: Sequence[str] = BUILD_CONFIGURATION_VARIABLES,
) -> bytes:
    require(tuple(variables) in (
                BUILD_CONFIGURATION_VARIABLES,
                BUILD_CONFIGURATION_VARIABLES_V8,
                BUILD_CONFIGURATION_VARIABLES_V7,
                BUILD_CONFIGURATION_VARIABLES_V6,
                BUILD_CONFIGURATION_VARIABLES_V5,
                BUILD_CONFIGURATION_VARIABLES_V4,
                BUILD_CONFIGURATION_VARIABLES_V3,
                BUILD_CONFIGURATION_VARIABLES_V2) and
            set(entries) == set(variables),
            "candidate normalized effective-configuration variables differ")
    lines: list[str] = []
    for name in variables:
        value = entries.get(name)
        require(isinstance(value, str) and
                not any(character in value for character in ("\0", "\r", "\n")),
                f"candidate normalized effective-configuration {name} is invalid")
        lines.append(f"{name}={value}\n")
    try:
        return "".join(lines).encode("utf-8", errors="strict")
    except UnicodeEncodeError as error:
        raise PlanError(
            "candidate normalized effective configuration is not strict UTF-8") \
            from error


def _parse_build_configuration_bytes(retained: bytes) -> dict[str, Any]:
    require(isinstance(retained, bytes) and
            0 < len(retained) <= MAX_BUILD_CONFIGURATION_BYTES and
            b"\0" not in retained and b"\r" not in retained and
            retained.endswith(b"\n"),
            "candidate normalized effective-configuration bytes are invalid")
    try:
        text = retained.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise PlanError(
            "candidate normalized effective configuration is not strict UTF-8") \
            from error
    lines = text[:-1].split("\n")
    require(lines and lines[0].startswith("schema="),
            "candidate normalized effective-configuration framing differs")
    file_schema = lines[0][len("schema="):]
    unused_record_schema, variables, unused_required_cache = \
        _build_configuration_contract(file_schema)
    del unused_record_schema, unused_required_cache
    require(len(lines) == len(variables) + 2 and
            lines[0] == f"schema={file_schema}" and
            lines[1].startswith("sha256="),
            "candidate normalized effective-configuration framing differs")
    digest = lines[1][len("sha256="):]
    require(re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
            "candidate normalized effective-configuration digest is invalid")
    entries: dict[str, str] = {}
    for expected, line in zip(variables, lines[2:]):
        name, separator, value = line.partition("=")
        require(separator == "=" and name == expected and name not in entries,
                "candidate normalized effective-configuration order differs")
        entries[name] = value
    material = _build_configuration_material(entries, variables)
    require(retained == (
                f"schema={file_schema}\n"
                f"sha256={digest}\n").encode("ascii") + material and
            hashlib.sha256(material).hexdigest() == digest,
            "candidate normalized effective-configuration digest differs")
    return {
        "configuration_schema": file_schema,
        "configuration_sha256": digest,
        "entries": entries,
    }


def _validate_embedded_build_type(
    entries: Mapping[str, str], embedded_build_type: object,
) -> str:
    require(isinstance(embedded_build_type, str) and embedded_build_type and
            isinstance(entries, Mapping),
            "candidate normalized embedded build type is invalid")
    generator = entries.get("CMAKE_GENERATOR")
    configured = entries.get("CMAKE_BUILD_TYPE")
    encoded_types = entries.get("CMAKE_CONFIGURATION_TYPES")
    require(all(isinstance(value, str)
                for value in (generator, configured, encoded_types)) and
            generator,
            "candidate normalized configuration omits generator semantics")
    multi = _cmake_generator_is_multi_config(generator)
    configuration_types: tuple[str, ...] = ()
    if encoded_types:
        configuration_types = tuple(encoded_types.split(";"))
        require(all(configuration_types) and
                len(configuration_types) == len(set(configuration_types)),
                "candidate normalized CMAKE_CONFIGURATION_TYPES is malformed")
    if multi:
        require(configuration_types and
                embedded_build_type in configuration_types,
                "candidate normalized multi-config build type is outside "
                "CMAKE_CONFIGURATION_TYPES")
    else:
        require(embedded_build_type == configured,
                "candidate normalized single-config build type differs from "
                "CMAKE_BUILD_TYPE")
    require(embedded_build_type == "Release",
            "candidate normalized authoritative benchmark is not Release")
    return embedded_build_type


def _cmake_generator_is_multi_config(generator: object) -> bool:
    require(isinstance(generator, str) and generator,
            "normalized CMake generator identity is invalid")
    return (generator == "Xcode" or
            generator.startswith("Visual Studio") or
            "Multi-Config" in generator)


def _cmake_build_layout(
    entries: Mapping[str, Any],
) -> tuple[bool, tuple[str, ...], str | None]:
    require(isinstance(entries, Mapping),
            "normalized CMake build-layout identity is invalid")
    generator = entries.get("CMAKE_GENERATOR")
    encoded_types = entries.get("CMAKE_CONFIGURATION_TYPES", "")
    build_type = entries.get("CMAKE_BUILD_TYPE", "")
    require(isinstance(encoded_types, str) and isinstance(build_type, str),
            "normalized CMake configuration values are invalid")
    multi = _cmake_generator_is_multi_config(generator)
    if multi:
        types = tuple(encoded_types.split(";")) if encoded_types else ()
        require(types == ("Debug", "Release") and
                generator == "Ninja Multi-Config",
                "normalized multi-config build lacks the canonical "
                "Debug/Release closure")
        return True, types, "Release"
    require(build_type == "Release",
            "normalized single-config build is not Release")
    return False, ("Release",), None


def _benchmark_attestation_text(commit: str, tree: str) -> str:
    require(re.fullmatch(r"[0-9a-f]{40}", commit) is not None and
            re.fullmatch(r"[0-9a-f]{40}", tree) is not None,
            "normalized benchmark attestation source identity is invalid")
    return (
        "#ifndef LEOPARD2_BENCHMARK_SOURCE_ATTESTATION_GENERATED_H\n"
        "#define LEOPARD2_BENCHMARK_SOURCE_ATTESTATION_GENERATED_H\n"
        "\n"
        "#undef LEO2_BENCHMARK_SOURCE_COMMIT\n"
        "#undef LEO2_BENCHMARK_SOURCE_TREE\n"
        "#undef LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY\n"
        f"#define LEO2_BENCHMARK_SOURCE_COMMIT \"{commit}\"\n"
        f"#define LEO2_BENCHMARK_SOURCE_TREE \"{tree}\"\n"
        "#define LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY 0\n"
        "\n"
        "#endif\n")


def _validate_normalized_benchmark_attestation(
    semantics: Mapping[str, Any], commit: str, tree: str,
    build_root: str, source_root: str,
) -> dict[str, Any]:
    attestation = semantics.get("generated_attestation_header")
    require(isinstance(attestation, dict) and set(attestation) == {
                "artifact", "content", "source_commit", "source_tree",
                "source_tracked_dirty"} and
            attestation.get("source_commit") == commit and
            attestation.get("source_tree") == tree and
            attestation.get("source_tracked_dirty") is False,
            "candidate normalized generated-attestation source identity differs")
    artifact = _validate_scope_artifact(
        attestation.get("artifact"),
        "candidate normalized generated attestation",
        "generated_compile_input")
    expected_path = (
        f"{build_root}/generated/leopard2-benchmark-attestation/"
        "leopard2_benchmark_source_attestation.h")
    require(artifact["path"] == expected_path,
            "candidate normalized generated-attestation path differs")
    content = _validate_scope_text(
        attestation.get("content"),
        "candidate normalized generated attestation")
    expected_text = _benchmark_attestation_text(commit, tree)
    require(content["text"] == expected_text and
            content["size"] == artifact["size"] and
            content["sha256"] == artifact["sha256"],
            "candidate normalized generated-attestation bytes differ")

    entries = semantics.get("required_entries")
    require(isinstance(entries, list),
            "candidate normalized generated-attestation compiler entry is absent")
    benchmark_entries = [
        entry for entry in entries
        if isinstance(entry, Mapping) and
           entry.get("file") == f"{source_root}/bench/leopard2/benchmark.cpp"]
    require(len(benchmark_entries) == 1 and
            isinstance(benchmark_entries[0].get("arguments"), list),
            "candidate normalized generated-attestation compiler entry differs")
    arguments = benchmark_entries[0]["arguments"]
    enable_definition_name = "LEO2_BENCHMARK_SOURCE_ATTESTATION"
    exact_enable_definition = f"-D{enable_definition_name}=1"
    enable_definitions = [
        token for token in arguments
        if isinstance(token, str) and
           (token == f"-D{enable_definition_name}" or
            token.startswith(f"-D{enable_definition_name}="))]
    header_definition_name = "LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER"
    header_definition_prefix = f"-D{header_definition_name}="
    exact_header_definition = (
        f'{header_definition_prefix}"{expected_path}"')
    header_definitions = [
        token for token in arguments
        if isinstance(token, str) and
           (token == f"-D{header_definition_name}" or
            token.startswith(header_definition_prefix))]
    forbidden_identity_prefixes = (
        "-DLEO2_BENCHMARK_SOURCE_COMMIT=",
        "-DLEO2_BENCHMARK_SOURCE_TREE=",
        "-DLEO2_BENCHMARK_SOURCE_TRACKED_DIRTY=",
    )
    require(enable_definitions == [exact_enable_definition] and
            not any(isinstance(token, str) and
                    (token == f"-U{enable_definition_name}" or
                     token.startswith(f"-U{enable_definition_name}="))
                    for token in arguments) and
            header_definitions == [exact_header_definition] and
            not any(isinstance(token, str) and
                    (token == f"-U{header_definition_name}" or
                     token.startswith(f"-U{header_definition_name}="))
                    for token in arguments) and
            not any(isinstance(token, str) and
                    token.startswith(forbidden_identity_prefixes)
                    for token in arguments) and
            f"-I{build_root}/generated/leopard2-benchmark-attestation"
                not in arguments,
            "candidate normalized generated-attestation compile input differs")
    return attestation


def _validate_normalized_build_configuration(
    semantics: Mapping[str, Any], cache: Mapping[str, Any],
    build_root: str, source_root: str,
) -> dict[str, Any]:
    value = semantics.get("effective_build_configuration")
    require(isinstance(value, dict) and set(value) == {
                "schema", "artifact", "content", "configuration_schema",
                "configuration_sha256", "entries", "embedded_build_type",
                "helper_source"},
            "candidate normalized effective-configuration shape differs")
    artifact = _validate_scope_artifact(
        value.get("artifact"), "candidate normalized effective configuration",
        "generated_build_configuration")
    expected_path = f"{build_root}/{BUILD_CONFIGURATION_RELATIVE_PATH}"
    require(artifact["path"] == expected_path,
            "candidate normalized effective-configuration path differs")
    content = _validate_scope_text(
        value.get("content"),
        "candidate normalized effective configuration")
    parsed = _parse_build_configuration_bytes(
        content["text"].encode("utf-8"))
    record_schema, unused_variables, required_entries = \
        _build_configuration_contract_for_compile_schema(
            str(semantics.get("schema")), parsed["configuration_schema"])
    del unused_variables
    _validate_embedded_build_type(
        parsed["entries"], value.get("embedded_build_type"))
    unused_multi, unused_types, selected_configuration = \
        _cmake_build_layout(parsed["entries"])
    del unused_multi, unused_types
    expected_cache_keys = {
        "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS_RELEASE",
        "CMAKE_CONFIGURATION_TYPES", "CMAKE_GENERATOR",
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA",
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256",
        *required_entries,
    }
    if selected_configuration:
        expected_cache_keys.add("CMAKE_MAKE_PROGRAM")
    require(content["size"] == artifact["size"] and
            content["sha256"] == artifact["sha256"] and
            value.get("schema") == record_schema and
            value.get("configuration_schema") ==
                parsed["configuration_schema"] and
            value.get("configuration_sha256") ==
                parsed["configuration_sha256"] and
            value.get("entries") == parsed["entries"] and
            set(cache) == expected_cache_keys and
            cache.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") ==
                parsed["configuration_schema"] and
            cache.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256") ==
                parsed["configuration_sha256"] and
            cache.get("CMAKE_GENERATOR") ==
                parsed["entries"]["CMAKE_GENERATOR"] and
            cache.get("CMAKE_CONFIGURATION_TYPES") ==
                parsed["entries"]["CMAKE_CONFIGURATION_TYPES"] and
            (selected_configuration is not None or
             cache.get("CMAKE_BUILD_TYPE") ==
                parsed["entries"]["CMAKE_BUILD_TYPE"]) and
            parsed["entries"]["CMAKE_CXX_COMPILER"] ==
                cache.get("CMAKE_CXX_COMPILER") and
            all(parsed["entries"].get(name) == expected
                for name, expected in
                _required_build_configuration_entries(
                    parsed["configuration_schema"],
                    str(semantics.get("schema"))).items()) and
            all(cache.get(name) == expected
                for name, expected in required_entries.items()),
            "candidate normalized effective configuration differs from "
            "cache/sidecar")
    helper = _validate_scope_artifact(
        value.get("helper_source"),
        "candidate normalized benchmark attestation helper", "source_file")
    require(helper["path"] ==
                f"{source_root}/{BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH}",
            "candidate normalized benchmark attestation helper path differs")
    return value


def _validate_scope_numeric_directory_inventory(
    value: object, prefix: str, label: str,
) -> list[str]:
    retained = _validate_scope_text(value, label)
    text = retained["text"]
    names = text.splitlines()
    require(names and text == "".join(name + "\n" for name in names) and
            all(name.startswith(prefix) and
                name.removeprefix(prefix).isdigit() for name in names) and
            names == sorted(set(names),
                            key=lambda name: int(name.removeprefix(prefix))),
            f"{label} is not a canonical numeric directory inventory")
    return names


def _validate_scope_commit_object(
    value: object, expected_commit: str, expected_tree: str, label: str,
) -> dict[str, Any]:
    require(isinstance(expected_commit, str) and
            re.fullmatch(r"[0-9a-f]{40}", expected_commit) is not None and
            isinstance(expected_tree, str) and
            re.fullmatch(r"[0-9a-f]{40}", expected_tree) is not None,
            f"{label} normalized expected Git identity is invalid")
    require(isinstance(value, dict) and set(value) == {
                "encoding", "size", "sha256", "object_id", "base64"} and
            value.get("encoding") == "base64" and
            isinstance(value.get("base64"), str),
            f"{label} normalized Git commit-object shape differs")
    try:
        raw = base64.b64decode(value["base64"], validate=True)
    except (ValueError, binascii.Error) as error:
        raise PlanError(
            f"{label} normalized Git commit object is not canonical base64") from error
    canonical = base64.b64encode(raw).decode("ascii")
    object_id = hashlib.sha1(
        f"commit {len(raw)}\0".encode("ascii") + raw).hexdigest()
    require(value["base64"] == canonical and
            type(value.get("size")) is int and value["size"] == len(raw) and
            0 < len(raw) <= MAX_GIT_COMMIT_BYTES and
            value.get("sha256") == hashlib.sha256(raw).hexdigest() and
            value.get("object_id") == expected_commit and
            object_id == expected_commit,
            f"{label} normalized Git commit-object byte identity differs")
    require(b"\n\n" in raw,
            f"{label} normalized Git commit has no header/message boundary")
    header_lines = raw.split(b"\n\n", 1)[0].splitlines()
    expected_tree_line = b"tree " + expected_tree.encode("ascii")
    trees = [line for line in header_lines if line.startswith(b"tree ")]
    require(header_lines and header_lines[0] == expected_tree_line and
            trees == [expected_tree_line],
            f"{label} normalized Git commit names another tree")
    return value


def _validate_scope_source(value: object, role: str) -> dict[str, Any]:
    require(role in {"baseline", "candidate"} and isinstance(value, dict) and
            set(value) == {
                "path", "head", "tree", "detached",
                "tracked_tree_listing_sha256", "tracked_status",
                "commit_object"},
            f"{role} normalized source identity shape differs")
    expected_path = "$BASELINE_SOURCE" if role == "baseline" else "$CANDIDATE_SOURCE"
    expected_head = EXACT_MAIN_COMMIT if role == "baseline" else None
    require(value.get("path") == expected_path and
            isinstance(value.get("head"), str) and
            re.fullmatch(r"[0-9a-f]{40}", value["head"]) is not None and
            (expected_head is None or value["head"] == expected_head) and
            isinstance(value.get("tree"), str) and
            re.fullmatch(r"[0-9a-f]{40}", value["tree"]) is not None and
            type(value.get("detached")) is bool and
            (role != "baseline" or value["detached"] is True) and
            isinstance(value.get("tracked_tree_listing_sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}",
                         value["tracked_tree_listing_sha256"]) is not None and
            value.get("tracked_status") == "clean",
            f"{role} normalized source identity differs")
    _validate_scope_commit_object(
        value.get("commit_object"), value["head"], value["tree"], role)
    return value


def _normalized_compile_output(
    role: str, source: str, selected_configuration: str | None = None,
) -> str:
    baseline = role == "baseline"
    if baseline:
        if source == ("$CANDIDATE_SOURCE/experiments/leopard2/main_compare/"
                      "legacy_main_benchmark.cpp"):
            return ("CMakeFiles/leopard_main_benchmark.dir/" +
                    (f"{selected_configuration}/"
                     if selected_configuration else "") +
                    "legacy_main_benchmark.cpp.o")
        require(source.startswith("$BASELINE_SOURCE/"),
                "baseline normalized compiler source escapes its root")
        relative = source[len("$BASELINE_SOURCE/"):]
        return ("CMakeFiles/leopard_main_exact.dir/" +
                (f"{selected_configuration}/"
                 if selected_configuration else "") +
                "$BASELINE_SOURCE/" +
                relative + ".o")
    require(role == "candidate" and
            source.startswith("$CANDIDATE_SOURCE/"),
            "candidate normalized compiler source escapes its root")
    relative = source[len("$CANDIDATE_SOURCE/"):]
    if relative == "bench/leopard2/benchmark.cpp":
        configuration = (
            f"{selected_configuration}/" if selected_configuration else "")
        return (
            f"CMakeFiles/bench_leopard2.dir/{configuration}"
            "bench/leopard2/benchmark.cpp.o")
    backend_targets = {
        "Leopard2BackendSSSE3.cpp": "leopard2_backend_ssse3.dir",
        "Leopard2BackendAVX2.cpp": "leopard2_backend_avx2.dir",
        "Leopard2BackendAVX2Xor.cpp": "leopard2_backend_avx2.dir",
        "Leopard2BackendAVX2T32B256.cpp":
            "leopard2_backend_avx2_t32_b256.dir",
        "Leopard2BackendAVX2T16B64.cpp":
            "leopard2_backend_avx2_t16_b64.dir",
        "Leopard2LowP32B64AVX2.cpp": "leopard2_low_p32_b64_avx2.dir",
        "Leopard2BackendAVX2T2K4.cpp": "leopard2_backend_avx2_t2_k4.dir",
        "Leopard2BackendAVX2T8K8B1024.cpp":
            "leopard2_backend_avx2_t8_k8_b1024.dir",
        "Leopard2BackendAVX512.cpp": "leopard2_backend_avx512.dir",
        "Leopard2BackendGFNI.cpp": "leopard2_backend_gfni.dir",
    }
    configuration = (
        f"{selected_configuration}/" if selected_configuration else "")
    return (
        f"CMakeFiles/{backend_targets.get(relative, 'leopard.dir')}/"
        f"{configuration}{relative}.o")


def _resolved_compiler_is_gnu(compiler_path: Path | str | None) -> bool:
    """Classify only an unambiguous resolved GNU C++ driver pathname."""
    if compiler_path is None:
        return False
    path = Path(compiler_path)
    return path.is_absolute() and \
        GNU_CXX_DRIVER_BASENAME.fullmatch(path.name) is not None


def _normalized_compile_argv(
    role: str, source: str, compiler_invocation: str,
    compile_schema: str = COMPILE_COMMANDS_SCHEMA,
    candidate_commit: str | None = None,
    candidate_tree: str | None = None,
    build_configuration: Mapping[str, Any] | None = None,
    selected_configuration: str | None = None,
    compiler_path: Path | str | None = None,
) -> list[str]:
    require(compile_schema in {
                COMPILE_COMMANDS_SCHEMA_V2, COMPILE_COMMANDS_SCHEMA_V3,
                COMPILE_COMMANDS_SCHEMA_V4, COMPILE_COMMANDS_SCHEMA_V5,
                COMPILE_COMMANDS_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V7,
                COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
                COMPILE_COMMANDS_SCHEMA_V10,
                COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA},
            "normalized compile argv schema differs")
    output = _normalized_compile_output(
        role, source, selected_configuration)
    if role == "baseline":
        adapter = ("$CANDIDATE_SOURCE/experiments/leopard2/main_compare/"
                   "legacy_main_benchmark.cpp")
        definitions = ([] if source != adapter else [
            f'-DLEOPARD_MAIN_SOURCE_COMMIT="{EXACT_MAIN_COMMIT}"',
        ])
        configuration_definition = (
            [f'-DCMAKE_INTDIR="{selected_configuration}"']
            if selected_configuration else [])
        isa_flags = (
            ["-march=x86-64", "-mtune=generic", "-mavx2",
             "-mno-avx512f"]
            if compile_schema in (COMPILE_COMMANDS_SCHEMA_V6,
                                  COMPILE_COMMANDS_SCHEMA_V7,
                                  COMPILE_COMMANDS_SCHEMA_V8,
                                  COMPILE_COMMANDS_SCHEMA_V9,
                                  COMPILE_COMMANDS_SCHEMA_V10,
                                  COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
            ["-march=native"])
        return [
            compiler_invocation, *definitions, *configuration_definition,
            "-I$BASELINE_SOURCE",
            "-g", "-O0", "-O3", "-std=gnu++11", *isa_flags,
            "-Wall", "-Wextra", "-fopenmp",
            "-o", output, "-c", source,
        ]

    relative = source[len("$CANDIDATE_SOURCE/"):]
    isolated_flags = {
        "Leopard2BackendSSSE3.cpp": ["-mssse3", "-mno-avx"],
        "Leopard2BackendAVX2.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2BackendAVX2Xor.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2BackendAVX2T32B256.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2BackendAVX2T16B64.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2LowP32B64AVX2.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2BackendAVX2T2K4.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2BackendAVX2T8K8B1024.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2BackendAVX512.cpp": [
            "-mavx2", "-mavx512f", "-mavx512bw", "-mavx512vl",
            "-mprefer-vector-width=256", "-falign-functions=64"],
        "Leopard2BackendGFNI.cpp": [
            "-mavx2", "-mgfni", "-mno-avx512f", "-falign-functions=64"],
    }
    if compile_schema in (
            COMPILE_COMMANDS_SCHEMA_V9, COMPILE_COMMANDS_SCHEMA_V10,
            COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) and \
            relative == "Leopard2BackendAVX2T8K8B1024.cpp" and \
            _resolved_compiler_is_gnu(compiler_path):
        isolated_flags[relative].insert(-1, "-flive-range-shrinkage")
    effective_avx2 = compile_schema in (
        COMPILE_COMMANDS_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V7,
        COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
        COMPILE_COMMANDS_SCHEMA_V10,
        COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA)
    source_definitions: list[str] = []
    global_definitions = ([
        "-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=1",
        "-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=1",
    ] if effective_avx2 else [])
    if relative == "Leopard2BackendAVX512.cpp":
        require(not effective_avx2,
                "current effective-AVX2 profile contains an AVX-512 source")
        definitions = ["-DLEO2_HAVE_AVX2_BACKEND=1"]
    elif relative == "Leopard2BackendGFNI.cpp":
        definitions = [*global_definitions,
            "-DLEO2_HAVE_AVX2_BACKEND=1", "-DLEO2_HAVE_GFNI_BACKEND=1"]
    elif relative == "Leopard2BackendAVX2T32B256.cpp":
        definitions = sorted([
            *global_definitions,
            *(["-DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=0",
               "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1"]
              if compile_schema == COMPILE_COMMANDS_SCHEMA else []),
            "-DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK=0",
            "-DLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
        ])
    elif relative == "Leopard2BackendAVX2T16B64.cpp":
        definitions = sorted([
            *global_definitions,
            "-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
        ])
    elif relative == "Leopard2LowP32B64AVX2.cpp":
        definitions = sorted([
            *global_definitions,
            "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
        ])
    elif relative == "Leopard2BackendAVX2T2K4.cpp":
        definitions = sorted([
            *global_definitions,
            "-DLEO2_HAVE_AVX2_BACKEND=1",
        ])
    elif relative == "Leopard2BackendAVX2T8K8B1024.cpp":
        require(compile_schema in (
                    COMPILE_COMMANDS_SCHEMA_V9,
                    COMPILE_COMMANDS_SCHEMA_V10,
                    COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA),
                "historical compile profile contains the K8/R8/B1024 source")
        definitions = sorted([
            *global_definitions,
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1",
        ])
    elif relative == "bench/leopard2/benchmark.cpp":
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V2:
            require(isinstance(candidate_commit, str) and
                    re.fullmatch(r"[0-9a-f]{40}", candidate_commit) is not None and
                    isinstance(candidate_tree, str) and
                    re.fullmatch(r"[0-9a-f]{40}", candidate_tree) is not None,
                    "historical normalized benchmark source identity differs")
            definitions = [
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1",
                f'-DLEO2_BENCHMARK_SOURCE_COMMIT="{candidate_commit}"',
                "-DLEO2_BENCHMARK_SOURCE_TRACKED_DIRTY=0",
                f'-DLEO2_BENCHMARK_SOURCE_TREE="{candidate_tree}"',
            ]
        elif compile_schema == COMPILE_COMMANDS_SCHEMA_V3:
            definitions = [
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1",
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER="
                '"$CANDIDATE_BUILD/generated/leopard2-benchmark-attestation/'
                'leopard2_benchmark_source_attestation.h"',
            ]
        else:
            configuration_file_schema = (
                build_configuration.get("configuration_schema")
                if isinstance(build_configuration, Mapping) else None)
            expected_record_schema, unused_variables, unused_cache = \
                _build_configuration_contract_for_compile_schema(
                    compile_schema, str(configuration_file_schema))
            del unused_variables, unused_cache
            require(isinstance(build_configuration, Mapping) and
                    build_configuration.get("schema") ==
                        expected_record_schema and
                    re.fullmatch(r"[0-9a-f]{64}", str(
                        build_configuration.get(
                            "configuration_sha256"))) is not None and
                    build_configuration.get("embedded_build_type") == "Release",
                    "normalized compile profile lacks effective configuration")
            definitions = [
                "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256="
                f'"{build_configuration["configuration_sha256"]}"',
                '-DLEO2_BENCHMARK_BUILD_TYPE="Release"',
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1",
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER="
                '"$CANDIDATE_BUILD/generated/leopard2-benchmark-attestation/'
                'leopard2_benchmark_source_attestation.h"',
                *global_definitions,
            ]
    elif relative in isolated_flags:
        definitions = list(global_definitions)
        if effective_avx2 and relative in {
                "Leopard2BackendAVX2.cpp",
                "Leopard2BackendAVX2Xor.cpp"}:
            definitions.insert(
                1, "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1")
        if effective_avx2 and relative == "Leopard2BackendAVX2.cpp":
            source_definitions.append(
                "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1")
    else:
        definitions = [
            "-DLEO2_DISABLE_AVX2_CODEGEN=1",
            "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
            *global_definitions[:1],
            *(["-DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=1",
               "-DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=1",
               "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1"]
              if effective_avx2 else []),
            *global_definitions[1:],
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            *([] if effective_avx2 else ["-DLEO2_HAVE_AVX512_BACKEND=1"]),
            "-DLEO2_HAVE_GFNI_BACKEND=1",
            "-DLEO2_HAVE_SSSE3_BACKEND=1",
        ]
        if compile_schema in (COMPILE_COMMANDS_SCHEMA_V8,
                              COMPILE_COMMANDS_SCHEMA_V9,
                              COMPILE_COMMANDS_SCHEMA_V10,
                              COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA):
            definitions.insert(
                definitions.index("-DLEO2_HAVE_AVX2_BACKEND=1"),
                "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1")
        if compile_schema in (
                COMPILE_COMMANDS_SCHEMA_V9,
                COMPILE_COMMANDS_SCHEMA_V10,
                COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA):
            definitions.insert(
                definitions.index("-DLEO2_HAVE_AVX2_BACKEND=1"),
                "-DLEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1")
        if effective_avx2 and relative in {
                "leopard2.cpp", "Leopard2Backend.cpp"}:
            source_definitions.append(
                "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1")
        if compile_schema in (COMPILE_COMMANDS_SCHEMA_V8,
                              COMPILE_COMMANDS_SCHEMA_V9,
                              COMPILE_COMMANDS_SCHEMA_V10,
                              COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) and \
                relative == "leopard2.cpp":
            source_definitions.extend([
                "-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=1",
                *(["-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1"]
                  if compile_schema == COMPILE_COMMANDS_SCHEMA else []),
                "-DLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=1",
            ])
        if compile_schema in (
                COMPILE_COMMANDS_SCHEMA_V10,
                COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) and \
                relative == "leopard2.cpp":
            source_definitions.append(
                "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1")
        if compile_schema in (
                COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) and \
                relative == "leopard2.cpp":
            source_definitions.extend([
                "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=1",
                "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=0",
            ])
    # CMake emits the target-wide definitions in lexical order for the
    # effective-AVX2 build profiles.  Keep this reconstruction byte-for-byte
    # aligned with the retained compile command; source-property definitions
    # remain a separate, following group.
    if effective_avx2:
        definitions.sort()
        source_definitions.sort()
    includes = ["-I$CANDIDATE_SOURCE"]
    configuration_definition = (
        [f'-DCMAKE_INTDIR="{selected_configuration}"']
        if selected_configuration else [])
    return [
        compiler_invocation, *definitions, *configuration_definition,
        *source_definitions, *includes,
        "-Wall", "-Wextra", "-fopenmp", "-O3", "-DNDEBUG", "-O3",
        "-std=gnu++11", *isolated_flags.get(relative, []),
        *([] if relative in isolated_flags else ["-fopenmp"]),
        "-o", output, "-c", source,
    ]


def _validate_scope_build(
    build: object, role: str, source_identity: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    require(role in {"baseline", "candidate"} and isinstance(build, dict),
            f"{role} normalized build identity is invalid")
    semantics = build.get("validated_compile_commands")
    require(isinstance(semantics, dict),
            f"{role} normalized compile-command identity is absent")
    compile_schema = semantics.get("schema")
    require(compile_schema in {
                COMPILE_COMMANDS_SCHEMA_V2, COMPILE_COMMANDS_SCHEMA_V3,
                COMPILE_COMMANDS_SCHEMA_V4, COMPILE_COMMANDS_SCHEMA_V5,
                COMPILE_COMMANDS_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V7,
                COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
                COMPILE_COMMANDS_SCHEMA_V10,
                COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA},
            f"{role} normalized compile-command schema differs")
    complete_build = compile_schema in (
        COMPILE_COMMANDS_SCHEMA_V4, COMPILE_COMMANDS_SCHEMA_V5,
        COMPILE_COMMANDS_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V7,
        COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
        COMPILE_COMMANDS_SCHEMA_V10,
        COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA)
    expected_build_keys = {
                "build_dir", "cmake_cache", "compile_commands",
                "executable_link_recipe", "archive_link_recipe", "compiler",
                "compiler_version_stdout", "archiver", "ranlib",
                "validated_archive_members", "validated_executable",
                "validated_archive", "validated_cache",
                "validated_compile_commands", "archive_link_recipe_content",
                "executable_link_recipe_content",
                "archive_link_tool_invocations", "compiler_invocation",
                "validated_external_link_inputs"}
    if complete_build:
        expected_build_keys.update({
            "multi_config_build_tool",
            "multi_config_build_tool_version_stdout",
            "multi_config_ninja_graph",
        })
    require(set(build) == expected_build_keys,
            f"{role} normalized build identity shape differs")
    expected_root = "$BASELINE_BUILD" if role == "baseline" else "$CANDIDATE_BUILD"
    require(build.get("build_dir") == expected_root,
            f"{role} normalized build directory differs")
    baseline = role == "baseline"
    cache = build.get("validated_cache")
    require(isinstance(cache, dict),
            f"{role} validated CMake cache is absent")
    if baseline:
        required_cache = (
            REQUIRED_PURE_AVX2_BASELINE_CACHE
            if compile_schema in (COMPILE_COMMANDS_SCHEMA_V6,
                                  COMPILE_COMMANDS_SCHEMA_V7,
                                  COMPILE_COMMANDS_SCHEMA_V8,
                                  COMPILE_COMMANDS_SCHEMA_V9,
                                  COMPILE_COMMANDS_SCHEMA_V10,
                                  COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
            REQUIRED_BASELINE_CACHE)
    elif complete_build:
        unused_record_schema, unused_variables, required_cache = \
            _build_configuration_contract_for_compile_schema(
                compile_schema, str(cache.get(
                    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA")))
        del unused_record_schema, unused_variables
    else:
        required_cache = REQUIRED_LEGACY_CANDIDATE_CACHE
    expected_cache_keys = {
        "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS_RELEASE",
        *required_cache,
    }
    selected_configuration: str | None = None
    configured_variant_count = 1
    if complete_build:
        expected_cache_keys.update({
            "CMAKE_CONFIGURATION_TYPES", "CMAKE_GENERATOR",
        })
        _, configuration_types, selected_configuration = \
            _cmake_build_layout(cache)
        configured_variant_count = (
            len(configuration_types) if selected_configuration else 1)
        if selected_configuration:
            expected_cache_keys.add("CMAKE_MAKE_PROGRAM")
        if not baseline:
            expected_cache_keys.update({
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA",
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256",
            })
    require(set(cache) == expected_cache_keys and
            isinstance(cache.get("CMAKE_CXX_COMPILER"), str) and
            cache["CMAKE_CXX_COMPILER"] and
            isinstance(cache.get("CMAKE_CXX_FLAGS_RELEASE"), str) and
            all(cache.get(key) == expected
                for key, expected in required_cache.items()) and
            (complete_build or cache.get("CMAKE_BUILD_TYPE") == "Release") and
            (not selected_configuration or
             isinstance(cache.get("CMAKE_MAKE_PROGRAM"), str) and
             bool(cache["CMAKE_MAKE_PROGRAM"])) and
            (baseline or not complete_build or
             cache.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") ==
                (BUILD_CONFIGURATION_FILE_SCHEMA
                 if compile_schema == COMPILE_COMMANDS_SCHEMA else
                 BUILD_CONFIGURATION_FILE_SCHEMA
                 if compile_schema == COMPILE_COMMANDS_SCHEMA_V11 else
                 BUILD_CONFIGURATION_FILE_SCHEMA_V8
                 if compile_schema == COMPILE_COMMANDS_SCHEMA_V10 else
                 BUILD_CONFIGURATION_FILE_SCHEMA_V7
                 if compile_schema == COMPILE_COMMANDS_SCHEMA_V9 else
                 BUILD_CONFIGURATION_FILE_SCHEMA_V6
                 if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
                 BUILD_CONFIGURATION_FILE_SCHEMA_V5
                 if compile_schema == COMPILE_COMMANDS_SCHEMA_V7 else
                 BUILD_CONFIGURATION_FILE_SCHEMA_V4
                 if compile_schema == COMPILE_COMMANDS_SCHEMA_V6 else
                 BUILD_CONFIGURATION_FILE_SCHEMA_V3
                 if compile_schema == COMPILE_COMMANDS_SCHEMA_V5 else
                 BUILD_CONFIGURATION_FILE_SCHEMA_V2)) and
            (baseline or not complete_build or re.fullmatch(
                r"[0-9a-f]{64}", str(cache.get(
                    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256")))
                is not None),
            f"{role} validated CMake cache differs")
    archive_name = "libleopard_main_exact.a" if baseline else "libleopard.a"
    executable_name = "leopard_main_benchmark" if baseline else "bench_leopard2"
    archive_target = "leopard_main_exact.dir" if baseline else "leopard.dir"
    library_sources = (
        BASELINE_LIBRARY_SOURCES if baseline else
        CANDIDATE_LIBRARY_SOURCES
        if compile_schema in (
            COMPILE_COMMANDS_SCHEMA_V9, COMPILE_COMMANDS_SCHEMA_V10,
            COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
        CANDIDATE_LIBRARY_SOURCES_V12
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
        CANDIDATE_LIBRARY_SOURCES_V11
        if compile_schema in (COMPILE_COMMANDS_SCHEMA_V6,
                              COMPILE_COMMANDS_SCHEMA_V7) else
        CANDIDATE_LIBRARY_SOURCES_V9)
    expected_entry_count = (
        BASELINE_EXPECTED_COMPILE_COMMAND_COUNT if baseline else
        (CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT
         if compile_schema in (
             COMPILE_COMMANDS_SCHEMA_V9, COMPILE_COMMANDS_SCHEMA_V10,
             COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
         CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT_V12
         if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
         CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT_V11)) * configured_variant_count
    output_prefix = (
        f"{selected_configuration}/" if selected_configuration else "")
    link_metadata = (
        f"{expected_root}/CMakeFiles/impl-{selected_configuration}.ninja"
        if selected_configuration else None)
    metadata_paths = {
        "cmake_cache": f"{expected_root}/CMakeCache.txt",
        "compile_commands": f"{expected_root}/compile_commands.json",
        "executable_link_recipe": (
            link_metadata or
            f"{expected_root}/CMakeFiles/{executable_name}.dir/link.txt"),
        "archive_link_recipe": (
            link_metadata or
            f"{expected_root}/CMakeFiles/{archive_target}/link.txt"),
    }
    for name in ("cmake_cache", "compile_commands", "executable_link_recipe",
                 "archive_link_recipe"):
        record = _validate_scope_artifact(
            build.get(name), f"{role} {name}", "build_metadata")
        require(record["path"] == metadata_paths[name],
                f"{role} normalized {name} path differs")
    compiler = _validate_scope_artifact(
        build.get("compiler"), f"{role} compiler", "compiler")
    archiver = _validate_scope_artifact(
        build.get("archiver"), f"{role} archiver", "archiver")
    ranlib = _validate_scope_artifact(
        build.get("ranlib"), f"{role} ranlib", "ranlib")
    archive = _validate_scope_artifact(
        build.get("validated_archive"), f"{role} archive", "archive")
    executable = _validate_scope_artifact(
        build.get("validated_executable"), f"{role} executable", "executable")
    require(archive["path"] ==
                f"{expected_root}/{output_prefix}{archive_name}" and
            executable["path"] ==
                f"{expected_root}/{output_prefix}{executable_name}",
            f"{role} normalized output paths differ")
    version = build.get("compiler_version_stdout")
    require(isinstance(version, dict) and set(version) == {"sha256", "text"} and
            isinstance(version.get("text"), str) and version["text"] and
            isinstance(version.get("sha256"), str) and
            version["sha256"] == hashlib.sha256(
                version["text"].encode("utf-8")).hexdigest(),
            f"{role} compiler version identity differs")
    if complete_build:
        build_tool = build.get("multi_config_build_tool")
        build_tool_version = build.get(
            "multi_config_build_tool_version_stdout")
        ninja_graph = build.get("multi_config_ninja_graph")
        if selected_configuration:
            tool = _validate_scope_artifact(
                build_tool, f"{role} normalized multi-config build tool",
                "build_tool")
            validated_tool_version = _validate_scope_build_tool_version(
                build_tool_version,
                f"{role} normalized multi-config build-tool version")
            require(
                cache.get("CMAKE_MAKE_PROGRAM") == CANONICAL_NINJA_PATH and
                tool["path"] == CANONICAL_NINJA_PATH and
                validated_tool_version is build_tool_version,
                f"{role} normalized multi-config Ninja identity/version "
                "differs")
            graph = _validate_scope_ninja_graph_closure(
                ninja_graph, expected_root, selected_configuration,
                f"{role} normalized multi-config Ninja graph")
            graph_by_path = {
                record["relative_path"]: record["artifact"]
                for record in graph["files"]
            }
            link_graph = graph_by_path.get(
                f"CMakeFiles/impl-{selected_configuration}.ninja")
            require(isinstance(link_graph, dict) and all(
                        link_graph[key] ==
                            build["archive_link_recipe"][key] ==
                            build["executable_link_recipe"][key]
                        for key in ("path", "size", "mode", "sha256")),
                    f"{role} normalized link metadata differs from its "
                    "retained Ninja graph")
        else:
            require(build_tool is None and build_tool_version is None and
                    ninja_graph is None,
                    f"{role} normalized single-config build retained a "
                    "multi-config build tool/graph")
    try:
        release_flags = shlex.split(
            cache["CMAKE_CXX_FLAGS_RELEASE"], posix=True)
    except ValueError as error:
        raise PlanError(
            f"cannot parse {role} normalized CMake Release flags: {error}") \
            from error
    try:
        common.validate_effective_flags(
            release_flags, f"{role} normalized CMake Release flags",
            "release")
    except common.EvidenceError as error:
        raise PlanError(str(error)) from error
    compiler_invocation = build.get("compiler_invocation")
    require(isinstance(compiler_invocation, dict) and
            set(compiler_invocation) == {"invocation", "resolved_path"} and
            compiler_invocation.get("invocation") ==
                cache["CMAKE_CXX_COMPILER"] and
            compiler_invocation.get("resolved_path") == compiler["path"],
            f"{role} normalized compiler invocation differs")
    expected_semantics_keys = {
                "entry_count", "required_sources", "validated_optimization",
                "validated_openmp", "required_source_object_pairs", "isa_policy",
                "schema", "implementation", "profile", "required_entries"}
    if compile_schema in (
            COMPILE_COMMANDS_SCHEMA_V3, COMPILE_COMMANDS_SCHEMA_V4,
            COMPILE_COMMANDS_SCHEMA_V5, COMPILE_COMMANDS_SCHEMA_V6,
            COMPILE_COMMANDS_SCHEMA_V7,
            COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
            COMPILE_COMMANDS_SCHEMA_V10,
            COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA):
        expected_semantics_keys.add("generated_attestation_header")
    if compile_schema in (
            COMPILE_COMMANDS_SCHEMA_V4, COMPILE_COMMANDS_SCHEMA_V5,
            COMPILE_COMMANDS_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V7,
            COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
            COMPILE_COMMANDS_SCHEMA_V10,
            COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA):
        expected_semantics_keys.add("effective_build_configuration")
    expected_profile = (
        (BASELINE_PURE_AVX2_COMPILE_PROFILE
         if compile_schema in (COMPILE_COMMANDS_SCHEMA_V6,
                               COMPILE_COMMANDS_SCHEMA_V7,
                               COMPILE_COMMANDS_SCHEMA_V8,
                               COMPILE_COMMANDS_SCHEMA_V9,
                               COMPILE_COMMANDS_SCHEMA_V10,
                               COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
         BASELINE_COMPILE_PROFILE) if baseline else
        CANDIDATE_COMPILE_PROFILE
        if compile_schema == COMPILE_COMMANDS_SCHEMA else
        CANDIDATE_COMPILE_PROFILE_V7
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V11 else
        CANDIDATE_COMPILE_PROFILE_V6
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V10 else
        CANDIDATE_COMPILE_PROFILE_V5
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V9 else
        CANDIDATE_COMPILE_PROFILE_V4
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
        CANDIDATE_COMPILE_PROFILE_V3
        if compile_schema in (COMPILE_COMMANDS_SCHEMA_V6,
                              COMPILE_COMMANDS_SCHEMA_V7) else
        CANDIDATE_COMPILE_PROFILE_V2
        if compile_schema in (
            COMPILE_COMMANDS_SCHEMA_V4, COMPILE_COMMANDS_SCHEMA_V5) else
        CANDIDATE_COMPILE_PROFILE_V1)
    expected_isa_policy = (
        (BASELINE_PURE_AVX2_ISA_POLICY
         if compile_schema in (COMPILE_COMMANDS_SCHEMA_V6,
                               COMPILE_COMMANDS_SCHEMA_V7,
                               COMPILE_COMMANDS_SCHEMA_V8,
                               COMPILE_COMMANDS_SCHEMA_V9,
                               COMPILE_COMMANDS_SCHEMA_V10,
                               COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
         BASELINE_NATIVE_ISA_POLICY) if baseline else
        CANDIDATE_ISA_POLICY
        if compile_schema in (
            COMPILE_COMMANDS_SCHEMA_V9, COMPILE_COMMANDS_SCHEMA_V10,
            COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
        CANDIDATE_ISA_POLICY_V12
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
        CANDIDATE_ISA_POLICY_V11
        if compile_schema in (COMPILE_COMMANDS_SCHEMA_V6,
                              COMPILE_COMMANDS_SCHEMA_V7) else
        CANDIDATE_ISA_POLICY_V9)
    require(set(semantics) == expected_semantics_keys and
            type(semantics.get("entry_count")) is int and
            semantics["entry_count"] == expected_entry_count and
            semantics.get("implementation") == role and
            semantics.get("profile") == expected_profile and
            semantics.get("validated_optimization") == "-O3" and
            semantics.get("validated_openmp") is True and
            semantics.get("isa_policy") == expected_isa_policy and
            isinstance(semantics.get("required_sources"), list) and
            semantics["required_sources"] and
            len(semantics["required_sources"]) == len(set(
                semantics["required_sources"])) and
            isinstance(semantics.get("required_source_object_pairs"), list) and
            semantics["required_source_object_pairs"] and
            isinstance(semantics.get("required_entries"), list) and
            len(semantics["required_entries"]) ==
                len(semantics["required_source_object_pairs"]),
            f"{role} normalized compile-command identity differs")
    entries_by_source: dict[str, dict[str, Any]] = {}
    for index, entry in enumerate(semantics["required_entries"]):
        require(isinstance(entry, dict) and set(entry) == {
                    "directory", "file", "output", "arguments"} and
                entry.get("directory") == expected_root and
                isinstance(entry.get("file"), str) and entry["file"] and
                isinstance(entry.get("output"), str) and entry["output"] and
                isinstance(entry.get("arguments"), list) and
                all(isinstance(token, str) and token and "@" not in token
                    for token in entry["arguments"]) and
                entry["file"] not in entries_by_source,
                f"{role} normalized compiler entry {index} differs")
        entries_by_source[entry["file"]] = entry
    sources: list[str] = []
    objects: list[str] = []
    for index, pair in enumerate(semantics["required_source_object_pairs"]):
        require(isinstance(pair, dict) and set(pair) == {"source", "object"},
                f"{role} normalized compile pair {index} shape differs")
        source = _validate_scope_artifact(
            pair["source"], f"{role} source {index}", "source_file")
        obj = _validate_scope_artifact(
            pair["object"], f"{role} object {index}", "object_file")
        entry = entries_by_source.get(source["path"])
        expected_output = _normalized_compile_output(
            role, source["path"], selected_configuration)
        require(entry == {
                    "directory": expected_root,
                    "file": source["path"],
                    "output": expected_output,
                    "arguments": _normalized_compile_argv(
                        role, source["path"],
                        compiler_invocation["invocation"], compile_schema,
                        (source_identity.get("head")
                         if role == "candidate" and
                         isinstance(source_identity, Mapping) else None),
                        (source_identity.get("tree")
                         if role == "candidate" and
                         isinstance(source_identity, Mapping) else None),
                        semantics.get("effective_build_configuration"),
                        selected_configuration,
                        compiler_path=compiler["path"]),
                } and
                obj["path"] == f"{expected_root}/{expected_output}",
                f"{role} normalized compiler argv/output {index} differs")
        sources.append(source["path"])
        objects.append(obj["path"])
    require(sources == semantics["required_sources"] and
            len(sources) == len(set(sources)) and
            len(objects) == len(set(objects)) and
            set(entries_by_source) == set(sources) and
            semantics["entry_count"] >= len(sources),
            f"{role} normalized compile source/object closure differs")
    benchmark_suffix = (
        "/experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        if baseline else "/bench/leopard2/benchmark.cpp")
    implementation_source = "$BASELINE_SOURCE" if baseline else "$CANDIDATE_SOURCE"
    expected_sources = sorted([
        *(f"{implementation_source}/{name}" for name in library_sources),
        "$CANDIDATE_SOURCE" + benchmark_suffix,
    ])
    require(len(sources) == len(expected_sources) and
            set(sources) == set(expected_sources),
            f"{role} normalized translation-unit set differs from the producer")
    benchmark_pairs = [
        pair for pair in semantics["required_source_object_pairs"]
        if pair["source"]["path"].endswith(benchmark_suffix)]
    require(len(benchmark_pairs) == 1 and
            all(path.startswith(expected_root + "/") for path in objects),
            f"{role} normalized benchmark/build object closure differs")
    if compile_schema in (
            COMPILE_COMMANDS_SCHEMA_V3, COMPILE_COMMANDS_SCHEMA_V4,
            COMPILE_COMMANDS_SCHEMA_V5, COMPILE_COMMANDS_SCHEMA_V6,
            COMPILE_COMMANDS_SCHEMA_V7,
            COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
            COMPILE_COMMANDS_SCHEMA_V10,
            COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA):
        if baseline:
            require(semantics.get("generated_attestation_header") is None,
                    "baseline normalized build unexpectedly has an attestation")
        else:
            require(isinstance(source_identity, Mapping),
                    "candidate normalized source identity is absent")
            _validate_normalized_benchmark_attestation(
                semantics, source_identity.get("head"),
                source_identity.get("tree"),
                "$CANDIDATE_BUILD", "$CANDIDATE_SOURCE")
    if compile_schema in (
            COMPILE_COMMANDS_SCHEMA_V4, COMPILE_COMMANDS_SCHEMA_V5,
            COMPILE_COMMANDS_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V7,
            COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
            COMPILE_COMMANDS_SCHEMA_V10,
            COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA):
        if baseline:
            require(semantics.get("effective_build_configuration") is None,
                    "baseline normalized build has candidate configuration")
        else:
            _validate_normalized_build_configuration(
                semantics, cache, "$CANDIDATE_BUILD", "$CANDIDATE_SOURCE")
    archive_pairs = [pair for pair in semantics["required_source_object_pairs"]
                     if pair not in benchmark_pairs]
    require(archive_pairs and
            all(pair["source"]["path"].startswith("$CANDIDATE_SOURCE/")
                for pair in benchmark_pairs) and
            all(pair["source"]["path"].startswith(implementation_source + "/")
                for pair in archive_pairs),
            f"{role} normalized compile sources escape their source roots")
    members = build.get("validated_archive_members")
    archive_source_order = library_sources
    if selected_configuration and not baseline:
        object_library_source_order = (
            "Leopard2BackendSSSE3.cpp",
            "Leopard2BackendAVX2.cpp",
            "Leopard2BackendAVX2Xor.cpp",
            "Leopard2BackendAVX512.cpp",
            "Leopard2BackendAVX2T32B256.cpp",
            "Leopard2BackendAVX2T16B64.cpp",
            "Leopard2LowP32B64AVX2.cpp",
            "Leopard2BackendAVX2T2K4.cpp",
            "Leopard2BackendAVX2T8K8B1024.cpp",
            "Leopard2BackendGFNI.cpp",
        )
        object_library_sources = tuple(
            name for name in object_library_source_order
            if name in library_sources)
        archive_source_order = (
            *object_library_sources,
            *(name for name in library_sources
              if name not in object_library_sources),
        )
    expected_members = [
        Path(name).name + ".o" for name in archive_source_order]
    require(isinstance(members, list) and members == expected_members and
            len(members) == len(set(members)) and
            all(isinstance(member, str) and member and "/" not in member
                for member in members),
            f"{role} normalized archive members differ")
    objects_by_member = {
        pair["object"]["path"].rsplit("/", 1)[-1]: pair["object"]["path"]
        for pair in archive_pairs
    }
    require(len(objects_by_member) == len(archive_pairs) and
            set(objects_by_member) == set(members),
            f"{role} normalized archive member/object closure differs")
    archive_text = _validate_scope_text(
        build.get("archive_link_recipe_content"),
        f"{role} archive link recipe")
    executable_text = _validate_scope_text(
        build.get("executable_link_recipe_content"),
        f"{role} executable link recipe")
    if not selected_configuration:
        require(
            archive_text["size"] == build["archive_link_recipe"]["size"] and
            archive_text["sha256"] ==
                build["archive_link_recipe"]["sha256"] and
            executable_text["size"] ==
                build["executable_link_recipe"]["size"] and
            executable_text["sha256"] ==
                build["executable_link_recipe"]["sha256"],
            f"{role} normalized retained link bytes differ from file identity")
    tools = build.get("archive_link_tool_invocations")
    require(isinstance(tools, dict) and set(tools) == {"archiver", "ranlib"},
            f"{role} archive tool invocation shape differs")
    for name, artifact in (("archiver", archiver), ("ranlib", ranlib)):
        invocation = tools.get(name)
        require(isinstance(invocation, dict) and set(invocation) == {
                    "invocation", "resolved_path"} and
                isinstance(invocation.get("invocation"), str) and
                invocation["invocation"] and
                invocation.get("resolved_path") == artifact["path"],
                f"{role} normalized {name} invocation differs")
    archive_commands = [
        shlex.split(line, posix=True)
        for line in archive_text["text"].splitlines() if line.strip()
    ]
    expected_archive_objects = [
        objects_by_member[member][len(expected_root) + 1:]
        for member in members
    ]
    require(len(archive_commands) == 2 and
            len(archive_commands[0]) >= 4 and
            archive_commands[0][0] == tools["archiver"]["invocation"] and
            archive_commands[0][1] in {"qc", "rc", "rcs"} and
            archive_commands[0][2] == output_prefix + archive_name and
            archive_commands[0][3:] == expected_archive_objects and
            archive_commands[1] == [
                tools["ranlib"]["invocation"], output_prefix + archive_name],
            f"{role} normalized archive recipe semantics differ")
    try:
        executable_tokens = common.parse_single_executable_recipe(
            executable_text["text"],
            f"{role} normalized executable link recipe")
    except common.EvidenceError as error:
        raise PlanError(str(error)) from error
    expected_benchmark_object = benchmark_pairs[0]["object"]["path"][
        len(expected_root) + 1:]
    external_link_inputs = build.get("validated_external_link_inputs")
    require(isinstance(external_link_inputs, list),
            f"{role} normalized external-link identity is absent")
    for index, record in enumerate(external_link_inputs):
        require(isinstance(record, dict) and
                isinstance(record.get("artifact"), dict),
                f"{role} normalized external-link identity {index} differs")
        _validate_scope_artifact(
            record["artifact"],
            f"{role} normalized external-link artifact {index}")
    semantic_archive_name = output_prefix + archive_name
    semantic_executable_name = output_prefix + executable_name
    semantic_tokens = executable_tokens
    if selected_configuration:
        semantic_tokens = [
            archive_name if token == semantic_archive_name else
            executable_name if token == semantic_executable_name else token
            for token in executable_tokens
        ]
        semantic_archive_name = archive_name
        semantic_executable_name = executable_name
    try:
        common.validate_executable_link_semantics(
            semantic_tokens,
            compiler_invocation=compiler_invocation["invocation"],
            archive_name=semantic_archive_name,
            executable_name=semantic_executable_name,
            benchmark_object=expected_benchmark_object,
            external_link_inputs=external_link_inputs,
            label=f"{role} normalized executable link recipe")
    except common.EvidenceError as error:
        raise PlanError(
            f"{role} normalized executable recipe semantics differ: {error}") \
            from error
    require(compiler["path"], f"{role} normalized compiler path is empty")
    return build


def _parse_scope_ldd_output(
    value: str, label: str, *, canonical: bool,
) -> list[dict[str, Any]]:
    address = (re.escape(CANONICAL_LDD_ADDRESS) if canonical else
               r"0x[0-9A-Fa-f]+")
    shared = re.compile(
        rf"(?P<soname>[^\s=()]+)\s+=>\s+(?P<path>/\S+)\s+"
        rf"\((?P<address>{address})\)")
    direct = re.compile(
        rf"(?P<path>/\S+)\s+\((?P<address>{address})\)")
    virtual = re.compile(
        rf"(?P<soname>linux-vdso\.so\.1)\s+"
        rf"\((?P<address>{address})\)")
    entries: list[dict[str, Any]] = []
    for raw_line in value.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        match = shared.fullmatch(line)
        if match is not None:
            target = match.group("path")
            require(os.path.normpath(target) == target,
                    f"{label} retains a non-canonical dependency")
            entries.append({
                "soname": match.group("soname"), "loader_path": target,
                "file_kind": "shared_library",
            })
            continue
        if "=>" in line and "not found" in line:
            raise PlanError(f"{label} retains a missing dependency")
        match = direct.fullmatch(line)
        if match is not None:
            target = match.group("path")
            require(os.path.normpath(target) == target,
                    f"{label} retains a non-canonical dynamic-loader path")
            entries.append({
                "soname": Path(target).name, "loader_path": target,
                "file_kind": "dynamic_loader",
            })
            continue
        match = virtual.fullmatch(line)
        if match is not None:
            entries.append({"soname": match.group("soname"), "virtual": True})
            continue
        raise PlanError(f"{label} contains unrecognized ldd output: {line}")
    require(entries and
            len({entry["soname"] for entry in entries}) == len(entries) and
            sum(entry.get("file_kind") == "dynamic_loader"
                for entry in entries) == 1 and
            any(entry.get("file_kind") == "shared_library"
                for entry in entries) and
            sum(entry.get("virtual") is True for entry in entries) <= 1,
            f"{label} is not an internally consistent, unique loader record")
    if canonical:
        rendered = []
        for entry in entries:
            if entry.get("virtual") is True:
                rendered.append(
                    f"{entry['soname']} ({CANONICAL_LDD_ADDRESS})")
            elif entry.get("file_kind") == "dynamic_loader":
                rendered.append(
                    f"{entry['loader_path']} ({CANONICAL_LDD_ADDRESS})")
            else:
                rendered.append(
                    f"{entry['soname']} => {entry['loader_path']} "
                    f"({CANONICAL_LDD_ADDRESS})")
        require(value == "\n".join(rendered) + "\n",
                f"{label} is not the exact canonical ldd rendering")
    return sorted(entries, key=lambda item: item["soname"])


def _validate_runtime_closure(value: object, label: str) -> dict[str, Any]:
    require(isinstance(value, dict),
            f"{label} normalized runtime closure is not an object")
    canonical = "canonical_ldd_output" in value
    output_key = "canonical_ldd_output" if canonical else "raw_ldd_output"
    require(set(value) == {"executable", "dependencies", output_key} and
            isinstance(value.get("executable"), str) and value["executable"] and
            isinstance(value.get("dependencies"), list) and value["dependencies"],
            f"{label} normalized runtime closure shape differs")
    output = value.get(output_key)
    if canonical:
        require(isinstance(output, dict) and set(output) == {
                    "schema", "normalization", "encoding", "size", "sha256", "text"} and
                output.get("schema") == CANONICAL_LDD_SCHEMA and
                output.get("normalization") == CANONICAL_LDD_NORMALIZATION,
                f"{label} normalized canonical ldd-output contract differs")
        retained = _validate_scope_text(
            {key: output[key] for key in ("encoding", "size", "sha256", "text")},
            f"{label} normalized canonical ldd output")
    else:
        retained = _validate_scope_text(
            output, f"{label} normalized raw ldd output")
    parsed = _parse_scope_ldd_output(
        retained["text"], f"{label} normalized ldd output",
        canonical=canonical)
    sonames: list[str] = []
    loader_paths: list[str] = []
    file_paths: list[str] = []
    for dependency in value["dependencies"]:
        require(isinstance(dependency, dict) and
                isinstance(dependency.get("soname"), str) and
                dependency["soname"],
                f"{label} normalized runtime dependency is invalid")
        sonames.append(dependency["soname"])
        if set(dependency) == {"soname", "virtual"}:
            require(dependency["soname"] == "linux-vdso.so.1" and
                    dependency.get("virtual") is True,
                    f"{label} normalized virtual dependency differs")
            continue
        require(set(dependency) == {"soname", "loader_path", "file"} and
                isinstance(dependency.get("loader_path"), str) and
                dependency["loader_path"],
                f"{label} normalized runtime dependency variant differs")
        file_record = _validate_scope_artifact(
            dependency["file"], f"{label} runtime {dependency['soname']}")
        require(file_record["kind"] in {"shared_library", "dynamic_loader"} and
                file_record["path"] == dependency["loader_path"],
                f"{label} normalized runtime dependency kind/loader path differs")
        loader_paths.append(dependency["loader_path"])
        file_paths.append(file_record["path"])
    require(sonames == sorted(sonames) and len(sonames) == len(set(sonames)) and
            len(loader_paths) == len(set(loader_paths)) and
            len(file_paths) == len(set(file_paths)),
            f"{label} normalized runtime closure is not sorted and unique")
    normalized_dependencies = []
    for dependency in value["dependencies"]:
        if dependency.get("virtual") is True:
            normalized_dependencies.append(dict(dependency))
        else:
            normalized_dependencies.append({
                "soname": dependency["soname"],
                "loader_path": dependency["loader_path"],
                "file_kind": dependency["file"]["kind"],
            })
    require(normalized_dependencies == parsed,
            f"{label} normalized dependency summary differs from raw ldd bytes")
    return value


def _validate_scope_cpu_policy(
    value: object, label: str, expected_sibling: int,
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "cpu", "online", "cpuinfo", "topology", "frequency_policy",
                "cache_hierarchy", "cache_index_inventory",
                "cache_directory_inventory_text", "numa_nodes",
                "numa_node_inventory", "core_class"} and
            type(value.get("cpu")) is int and
            0 <= value["cpu"] <= MAX_CPU_ID and
            (value.get("online") is None or isinstance(value["online"], str)),
            f"{label} normalized CPU policy shape differs")
    cpu = value["cpu"]
    cpuinfo = value.get("cpuinfo")
    allowed_cpuinfo = {
        "processor", "vendor_id", "cpu family", "model", "model name",
        "stepping", "microcode", "flags", "Features", "CPU implementer",
        "CPU architecture", "CPU variant", "CPU part", "CPU revision",
    }
    require(isinstance(cpuinfo, dict) and cpuinfo and
            set(cpuinfo).issubset(allowed_cpuinfo) and
            cpuinfo.get("processor") == str(cpu) and
            any(name in cpuinfo for name in ("model name", "CPU part")) and
            all(isinstance(item, str) for item in cpuinfo.values()),
            f"{label} normalized CPU model identity differs")
    topology = value.get("topology")
    require(isinstance(topology, dict) and set(topology) == {
                "core_id", "physical_package_id", "die_id", "cluster_id",
                "thread_siblings_list", "core_siblings_list"} and
            all(item is None or isinstance(item, str)
                for item in topology.values()) and
            _parse_scope_cpu_list(topology.get("thread_siblings_list"),
                                  f"{label} thread siblings") == {
                cpu, expected_sibling} and
            {cpu, expected_sibling}.issubset(_parse_scope_cpu_list(
                topology.get("core_siblings_list"),
                f"{label} core siblings")),
            f"{label} normalized CPU topology differs")
    frequency = value.get("frequency_policy")
    require(isinstance(frequency, dict) and set(frequency) == {
                "scaling_driver", "scaling_governor",
                "energy_performance_preference", "scaling_min_freq",
                "scaling_max_freq", "cpuinfo_min_freq", "cpuinfo_max_freq"} and
            all(item is None or isinstance(item, str)
                for item in frequency.values()),
            f"{label} normalized frequency identity differs")
    caches = value.get("cache_hierarchy")
    require(isinstance(caches, list) and caches,
            f"{label} normalized cache hierarchy is absent")
    indices = []
    for cache in caches:
        require(isinstance(cache, dict) and set(cache) == {
                    "index", "level", "type", "size", "coherency_line_size",
                    "number_of_sets", "ways_of_associativity",
                    "physical_line_partition", "shared_cpu_list", "shared_cpu_map",
                    "allocation_policy", "write_policy"} and
                type(cache.get("index")) is int and cache["index"] >= 0 and
                all(cache.get(name) is None or isinstance(cache[name], str)
                    for name in set(cache) - {"index"}) and
                all(isinstance(cache.get(name), str) and cache[name]
                    for name in ("level", "type", "size", "coherency_line_size",
                                 "shared_cpu_list", "shared_cpu_map")) and
                cpu in _parse_scope_cpu_list(
                    cache["shared_cpu_list"], f"{label} cache shared CPUs"),
                f"{label} normalized cache identity differs")
        indices.append(cache["index"])
    require(indices == sorted(indices) and len(indices) == len(set(indices)),
            f"{label} normalized cache indices differ")
    cache_inventory = value.get("cache_index_inventory")
    raw_cache_inventory = _validate_scope_numeric_directory_inventory(
        value.get("cache_directory_inventory_text"), "index",
        f"{label} raw cache-directory inventory")
    require(isinstance(cache_inventory, list) and
            cache_inventory == [f"index{index}" for index in indices] and
            cache_inventory == raw_cache_inventory,
            f"{label} normalized cache inventory differs")
    numa = value.get("numa_nodes")
    node_inventory = value.get("numa_node_inventory")
    core_class = value.get("core_class")
    require(isinstance(numa, list) and numa and numa == sorted(set(numa)) and
            all(type(node) is int and 0 <= node <= MAX_CPU_ID for node in numa) and
            isinstance(node_inventory, list) and
            node_inventory == [f"node{node}" for node in numa] and
            isinstance(core_class, dict) and set(core_class) == {
                "core_type", "cpu_capacity"} and
            all(item is None or isinstance(item, str)
                for item in core_class.values()),
            f"{label} normalized NUMA/core-class identity differs")
    return value


def _validate_scope_host(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "system", "allowed_cpu_set_at_launch", "online_cpu_set",
                "online_cpu_list_text", "online_node_list_text",
                "benchmark_cpu", "reserved_sibling", "turbo_and_pstate"},
            "gate normalized host identity shape differs")
    system = value.get("system")
    require(isinstance(system, dict) and set(system) == {
                "system", "node", "release", "version", "machine", "python",
                "page_size"} and
            all(isinstance(system.get(name), str) and system[name]
                for name in ("system", "node", "release", "version", "machine",
                             "python")) and
            type(system.get("page_size")) is int and system["page_size"] > 0,
            "gate normalized host system identity differs")
    allowed = value.get("allowed_cpu_set_at_launch")
    online = value.get("online_cpu_set")
    require(isinstance(allowed, list) and allowed == sorted(set(allowed)) and
            3 <= len(allowed) <= MAX_CPU_LIST_ENTRIES and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in allowed) and
            isinstance(online, list) and online == sorted(set(online)) and
            1 <= len(online) <= MAX_CPU_LIST_ENTRIES and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in online),
            "gate normalized host CPU sets differ")
    raw_online = _validate_scope_text(
        value.get("online_cpu_list_text"), "gate normalized online CPU list")
    raw_nodes = _validate_scope_text(
        value.get("online_node_list_text"), "gate normalized online NUMA list")
    online_nodes = _parse_scope_cpu_list(
        raw_nodes["text"], "gate normalized online NUMA list")
    require(_parse_scope_cpu_list(
                raw_online["text"], "gate normalized online CPU list") ==
                set(online),
            "gate normalized online CPU summary differs from sysfs bytes")
    benchmark = value.get("benchmark_cpu")
    sibling = value.get("reserved_sibling")
    require(isinstance(benchmark, dict) and isinstance(sibling, dict) and
            type(benchmark.get("cpu")) is int and
            type(sibling.get("cpu")) is int and
            benchmark["cpu"] != sibling["cpu"],
            "gate normalized host reserved pair differs")
    _validate_scope_cpu_policy(
        benchmark, "benchmark CPU", sibling["cpu"])
    _validate_scope_cpu_policy(
        sibling, "reserved sibling", benchmark["cpu"])
    pair = {benchmark["cpu"], sibling["cpu"]}
    require(pair.issubset(set(allowed)) and pair.issubset(set(online)) and
            benchmark["cache_hierarchy"] == sibling["cache_hierarchy"] and
            benchmark["cache_index_inventory"] ==
                sibling["cache_index_inventory"] and
            benchmark["numa_nodes"] == sibling["numa_nodes"] and
            benchmark["numa_node_inventory"] ==
                sibling["numa_node_inventory"] and
            set(benchmark["numa_nodes"]).issubset(online_nodes) and
            benchmark["core_class"] == sibling["core_class"],
            "gate normalized pair is outside CPU sets or differs in cache/class/domain")
    turbo = value.get("turbo_and_pstate")
    require(isinstance(turbo, dict) and set(turbo) == {
                "/sys/devices/system/cpu/intel_pstate/no_turbo",
                "/sys/devices/system/cpu/amd_pstate/status",
                "/sys/devices/system/cpu/cpufreq/boost"} and
            all(item is None or isinstance(item, str) for item in turbo.values()),
            "gate normalized host turbo/pstate identity differs")
    return value


def selection_scope_from_verified_bundle(
    manifest: dict[str, Any], raw: dict[str, Any]
) -> dict[str, Any]:
    """Derive the one environment in which every gate cell is meaningful."""
    require(raw.get("schema") in EXACT_RAW_SCHEMAS,
            "exact-main raw bundle is not a supported complete-identity schema")
    identities = manifest.get("identities")
    host = manifest.get("host")
    require(isinstance(identities, dict) and isinstance(host, dict),
            "exact-main scope identity is incomplete")
    baseline_build = identities.get("baseline_build")
    candidate_build = identities.get("candidate_build")
    baseline_source = identities.get("baseline_source")
    candidate_source = identities.get("candidate_source")
    require(all(isinstance(item, dict) for item in (
        baseline_build, candidate_build, baseline_source, candidate_source)),
        "exact-main baseline/candidate build or source scope is absent")
    baseline_source_root = baseline_source.get("path")
    candidate_source_root = candidate_source.get("path")
    baseline_build_root = baseline_build.get("build_dir")
    candidate_build_root = candidate_build.get("build_dir")
    require(all(isinstance(item, str) and item for item in (
        baseline_source_root, candidate_source_root,
        baseline_build_root, candidate_build_root)),
        "exact-main role-specific source/build roots are absent")

    backends = {
        invocation.get("normalized", {}).get("backend")
        for invocation in raw.get("invocations", [])
        if isinstance(invocation, dict) and
           invocation.get("implementation") == "candidate" and
           isinstance(invocation.get("normalized"), dict)
    }
    require(len(backends) == 1, "exact-main candidate resolved multiple/no backends")
    resolved_backend = next(iter(backends))
    require(resolved_backend in CAMPAIGN_BACKENDS,
            "resolved AUTO backend lacks forced confirmation coverage; "
            "AVX-512 and NEON are outside this campaign")
    if raw.get("schema") in (
            EXACT_RAW_SCHEMA_V10, EXACT_RAW_SCHEMA_V11,
            EXACT_RAW_SCHEMA_V12, EXACT_RAW_SCHEMA_V13,
            EXACT_RAW_SCHEMA_V14, EXACT_RAW_SCHEMA_V15, EXACT_RAW_SCHEMA):
        require(resolved_backend == "avx2",
                "current exact-main candidate did not resolve the effective-AVX2 "
                "AUTO backend")

    replacements = (
        (baseline_source_root, "$BASELINE_SOURCE"),
        (candidate_source_root, "$CANDIDATE_SOURCE"),
        (baseline_build_root, "$BASELINE_BUILD"),
        (candidate_build_root, "$CANDIDATE_BUILD"),
    )
    raw_schema = raw.get("schema")
    rich_scope = raw_schema in {
        EXACT_RAW_SCHEMA_V7, EXACT_RAW_SCHEMA_V8, EXACT_RAW_SCHEMA_V9,
        EXACT_RAW_SCHEMA_V10, EXACT_RAW_SCHEMA_V11, EXACT_RAW_SCHEMA_V12,
        EXACT_RAW_SCHEMA_V13, EXACT_RAW_SCHEMA_V14, EXACT_RAW_SCHEMA_V15,
        EXACT_RAW_SCHEMA}
    tool_names = ["runner", "taskset", "ldd"]
    if rich_scope:
        tool_names.append("evidence_helper")
    source_scope_keys = (
        "path", "head", "tree", "detached",
        "tracked_tree_listing_sha256", "tracked_status", "commit_object",
    )
    # v8 carries the complete recursive Git/tree proof.  The producer verifies
    # that proof before returning the bundle; the selection scope intentionally
    # projects only its mathematical source identity so path-local Git metadata
    # and guard implementation details do not make otherwise identical gate
    # campaigns incomparable.  The complete manifest snapshot remains retained.
    baseline_source_scope = {
        key: baseline_source[key] for key in source_scope_keys}
    candidate_source_scope = {
        key: candidate_source[key] for key in source_scope_keys}
    scope = {
        "schema": (
            EVIDENCE_SCOPE_SCHEMA if raw_schema == EXACT_RAW_SCHEMA else
            EVIDENCE_SCOPE_SCHEMA_V11
            if raw_schema == EXACT_RAW_SCHEMA_V15 else
            EVIDENCE_SCOPE_SCHEMA_V10
            if raw_schema == EXACT_RAW_SCHEMA_V14 else
            EVIDENCE_SCOPE_SCHEMA_V9
            if raw_schema == EXACT_RAW_SCHEMA_V13 else
            EVIDENCE_SCOPE_SCHEMA_V8
            if raw_schema == EXACT_RAW_SCHEMA_V12 else
            EVIDENCE_SCOPE_SCHEMA_V7
            if raw_schema == EXACT_RAW_SCHEMA_V11 else
            EVIDENCE_SCOPE_SCHEMA_V6
            if raw_schema == EXACT_RAW_SCHEMA_V10 else
            EVIDENCE_SCOPE_SCHEMA_V5
            if raw_schema == EXACT_RAW_SCHEMA_V9 else
            EVIDENCE_SCOPE_SCHEMA_V4 if rich_scope else
            EVIDENCE_SCOPE_SCHEMA_V3),
        # Retain the exact pinned CPU pair.  This is intentionally stricter than
        # topology cardinality: the current evidence format lacks cache/NUMA and
        # heterogeneous-core descriptors needed to equate different pairs.
        "host": _normalize_bound_paths(host, replacements),
        "sources": {
            "baseline": _normalize_bound_paths(
                baseline_source_scope, replacements),
            "candidate": _normalize_bound_paths(
                candidate_source_scope, replacements),
        },
        "builds": {
            "baseline": _normalize_bound_paths(baseline_build, replacements),
            "candidate": _normalize_bound_paths(candidate_build, replacements),
        },
        "artifacts": {
            key: _normalize_bound_paths(identities.get(key), replacements)
            for key in (
                "baseline_archive", "baseline_executable",
                "candidate_archive", "candidate_executable")
        },
        "runtime_closures": {
            "baseline": _normalize_bound_paths(
                identities.get("baseline_runtime_closure"), replacements),
            "candidate": _normalize_bound_paths(
                identities.get("candidate_runtime_closure"), replacements),
        },
        "tools": {
            key: _normalize_bound_paths(identities.get(key), replacements)
            for key in tool_names
        },
        "resolved_auto_backend": resolved_backend,
        "forced_confirmation_backends": list(
            BACKENDS[:BACKENDS.index(resolved_backend) + 1]),
        "excluded_backends": dict(EXCLUDED_CAMPAIGN_BACKENDS),
    }
    return scope


def validate_evidence_scope(scope: object) -> dict[str, Any]:
    scope_schema = scope.get("schema") if isinstance(scope, dict) else None
    require(isinstance(scope, dict) and set(scope) == {
        "schema", "host", "sources", "builds", "artifacts",
        "runtime_closures", "tools",
        "resolved_auto_backend", "forced_confirmation_backends",
        "excluded_backends",
    } and scope_schema in {
        EVIDENCE_SCOPE_SCHEMA_V3, EVIDENCE_SCOPE_SCHEMA_V4,
        EVIDENCE_SCOPE_SCHEMA_V5, EVIDENCE_SCOPE_SCHEMA_V6,
        EVIDENCE_SCOPE_SCHEMA_V7, EVIDENCE_SCOPE_SCHEMA_V8,
        EVIDENCE_SCOPE_SCHEMA_V9, EVIDENCE_SCOPE_SCHEMA_V10,
        EVIDENCE_SCOPE_SCHEMA_V11, EVIDENCE_SCOPE_SCHEMA},
            "gate evidence scope shape differs")
    backend = scope.get("resolved_auto_backend")
    require(backend in CAMPAIGN_BACKENDS and
            (scope_schema not in (EVIDENCE_SCOPE_SCHEMA_V6,
                                  EVIDENCE_SCOPE_SCHEMA_V7,
                                  EVIDENCE_SCOPE_SCHEMA_V8,
                                  EVIDENCE_SCOPE_SCHEMA_V9,
                                  EVIDENCE_SCOPE_SCHEMA_V10,
                                  EVIDENCE_SCOPE_SCHEMA_V11, EVIDENCE_SCOPE_SCHEMA) or
             backend == "avx2") and
            scope.get("forced_confirmation_backends") ==
                list(BACKENDS[:BACKENDS.index(backend) + 1]) and
            scope.get("excluded_backends") == EXCLUDED_CAMPAIGN_BACKENDS,
            "gate evidence backend coverage declaration differs")
    host = scope.get("host")
    sources = scope.get("sources")
    builds = scope.get("builds")
    artifacts = scope.get("artifacts")
    closures = scope.get("runtime_closures")
    tools = scope.get("tools")
    _validate_scope_host(host)
    require(isinstance(sources, dict) and set(sources) == {"baseline", "candidate"} and
            isinstance(builds, dict) and set(builds) == {"baseline", "candidate"},
            "gate evidence role-specific source/build scope differs")
    _validate_scope_source(sources["baseline"], "baseline")
    _validate_scope_source(sources["candidate"], "candidate")
    _validate_scope_build(
        builds["baseline"], "baseline", sources["baseline"])
    _validate_scope_build(
        builds["candidate"], "candidate", sources["candidate"])
    baseline_compile = builds["baseline"]["validated_compile_commands"]
    candidate_compile = builds["candidate"]["validated_compile_commands"]
    require(baseline_compile["schema"] == candidate_compile["schema"],
            "gate evidence compile-command schema differs between roles")
    expected_scope_compile_schemas = {
        EVIDENCE_SCOPE_SCHEMA_V3: frozenset((
            COMPILE_COMMANDS_SCHEMA_V2, COMPILE_COMMANDS_SCHEMA_V3)),
        EVIDENCE_SCOPE_SCHEMA_V4: frozenset((COMPILE_COMMANDS_SCHEMA_V4,)),
        # v5 replays historical raw v9 evidence, where decode_first_use was
        # derived as plan setup plus execution.  v6 is deliberately distinct:
        # raw v10 measures the public one-shot decode call itself.
        EVIDENCE_SCOPE_SCHEMA_V5: frozenset((COMPILE_COMMANDS_SCHEMA_V5,)),
        EVIDENCE_SCOPE_SCHEMA_V6: frozenset((COMPILE_COMMANDS_SCHEMA_V6,)),
        EVIDENCE_SCOPE_SCHEMA_V7: frozenset((COMPILE_COMMANDS_SCHEMA_V7,)),
        EVIDENCE_SCOPE_SCHEMA_V8: frozenset((COMPILE_COMMANDS_SCHEMA_V8,)),
        EVIDENCE_SCOPE_SCHEMA_V9: frozenset((COMPILE_COMMANDS_SCHEMA_V9,)),
        EVIDENCE_SCOPE_SCHEMA_V10: frozenset((COMPILE_COMMANDS_SCHEMA_V10,)),
        EVIDENCE_SCOPE_SCHEMA_V11: frozenset((COMPILE_COMMANDS_SCHEMA_V11,)),
        EVIDENCE_SCOPE_SCHEMA: frozenset((COMPILE_COMMANDS_SCHEMA,)),
    }
    require(
        candidate_compile["schema"] in
            expected_scope_compile_schemas[scope_schema],
        "gate evidence scope schema differs from its compile-command schema")
    if candidate_compile["schema"] in (
        COMPILE_COMMANDS_SCHEMA_V3, COMPILE_COMMANDS_SCHEMA_V4,
        COMPILE_COMMANDS_SCHEMA_V5, COMPILE_COMMANDS_SCHEMA_V6,
        COMPILE_COMMANDS_SCHEMA_V7,
        COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
        COMPILE_COMMANDS_SCHEMA_V10,
        COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA,
    ):
        attestation = candidate_compile["generated_attestation_header"]
        require(attestation["source_commit"] == sources["candidate"]["head"] and
                attestation["source_tree"] == sources["candidate"]["tree"] and
                attestation["source_tracked_dirty"] is False and
                sources["candidate"]["tracked_status"] == "clean",
                "gate evidence generated attestation differs from source identity")
    require(builds["baseline"]["compiler"] == builds["candidate"]["compiler"] and
            builds["baseline"]["compiler_version_stdout"] ==
                builds["candidate"]["compiler_version_stdout"],
            "gate evidence baseline/candidate compiler identity differs")
    if scope_schema in (
            EVIDENCE_SCOPE_SCHEMA_V4, EVIDENCE_SCOPE_SCHEMA_V5,
            EVIDENCE_SCOPE_SCHEMA_V6, EVIDENCE_SCOPE_SCHEMA_V7,
            EVIDENCE_SCOPE_SCHEMA_V8,
            EVIDENCE_SCOPE_SCHEMA_V9,
            EVIDENCE_SCOPE_SCHEMA_V10,
            EVIDENCE_SCOPE_SCHEMA_V11, EVIDENCE_SCOPE_SCHEMA):
        require(
            builds["baseline"]["multi_config_build_tool"] ==
                builds["candidate"]["multi_config_build_tool"] and
            builds["baseline"][
                "multi_config_build_tool_version_stdout"] ==
                builds["candidate"][
                    "multi_config_build_tool_version_stdout"],
            "gate evidence baseline/candidate multi-config build-tool "
            "identity differs")
    require(builds["baseline"]["validated_external_link_inputs"] ==
                builds["candidate"]["validated_external_link_inputs"],
            "gate evidence baseline/candidate external link inputs differ")
    require(isinstance(artifacts, dict) and set(artifacts) == {
        "baseline_archive", "baseline_executable",
        "candidate_archive", "candidate_executable",
    }, "gate evidence artifact closure shape differs")
    for key, artifact_value in artifacts.items():
        _validate_scope_artifact(
            artifact_value, key,
            "archive" if key.endswith("archive") else "executable")
    require(artifacts["baseline_archive"] ==
                builds["baseline"]["validated_archive"] and
            artifacts["baseline_executable"] ==
                builds["baseline"]["validated_executable"] and
            artifacts["candidate_archive"] ==
                builds["candidate"]["validated_archive"] and
            artifacts["candidate_executable"] ==
                builds["candidate"]["validated_executable"],
            "gate evidence top-level/build artifact closure differs")
    require(isinstance(closures, dict) and set(closures) == {
        "baseline", "candidate",
    }, "gate evidence runtime closure shape differs")
    _validate_runtime_closure(closures["baseline"], "baseline")
    _validate_runtime_closure(closures["candidate"], "candidate")
    require(closures["baseline"]["executable"] ==
                artifacts["baseline_executable"].get("path") and
            closures["candidate"]["executable"] ==
                artifacts["candidate_executable"].get("path"),
            "gate evidence runtime closure names a different executable")
    expected_tools = {"runner", "taskset", "ldd"}
    if scope_schema in (
            EVIDENCE_SCOPE_SCHEMA_V4, EVIDENCE_SCOPE_SCHEMA_V5,
            EVIDENCE_SCOPE_SCHEMA_V6, EVIDENCE_SCOPE_SCHEMA_V7,
            EVIDENCE_SCOPE_SCHEMA_V8,
            EVIDENCE_SCOPE_SCHEMA_V9,
            EVIDENCE_SCOPE_SCHEMA_V10,
            EVIDENCE_SCOPE_SCHEMA_V11, EVIDENCE_SCOPE_SCHEMA):
        expected_tools.add("evidence_helper")
    require(isinstance(tools, dict) and set(tools) == expected_tools,
            "gate evidence tool scope shape differs")
    for key, tool in tools.items():
        _validate_scope_artifact(
            tool, key,
            "file" if key in {"runner", "evidence_helper"} else "executable")
    if scope_schema in (
            EVIDENCE_SCOPE_SCHEMA_V4, EVIDENCE_SCOPE_SCHEMA_V5,
            EVIDENCE_SCOPE_SCHEMA_V6, EVIDENCE_SCOPE_SCHEMA_V7,
            EVIDENCE_SCOPE_SCHEMA_V8,
            EVIDENCE_SCOPE_SCHEMA_V9,
            EVIDENCE_SCOPE_SCHEMA_V10,
            EVIDENCE_SCOPE_SCHEMA_V11, EVIDENCE_SCOPE_SCHEMA):
        require(
            tools["evidence_helper"]["path"] ==
                f"$CANDIDATE_SOURCE/{EVIDENCE_HELPER_RELATIVE_PATH}",
            "gate evidence helper is not bound to the candidate source")
    return scope


def _validate_exact_schema_pair(
    manifest: Mapping[str, Any], raw: Mapping[str, Any],
) -> None:
    """Reject coherent-looking cross-version or downgraded evidence pairs."""
    pair = (manifest.get("schema"), raw.get("schema"))
    require(pair in EXACT_SCHEMA_PAIRS,
            "exact-main manifest/raw schema versions do not form one "
            "supported evidence contract")


def verify_exact_manifest(
    path: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    try:
        document, raw, _, snapshot = \
            load_exact_main_runner().verified_campaign_bundle(
                path, no_current_input_check=True)
    except Exception as error:
        raise PlanError(f"exact-main verifier rejected {path}: {error}") from error
    require(isinstance(document, dict) and
            document.get("schema") in EXACT_MANIFEST_SCHEMAS and
            document.get("valid") is True,
            f"exact-main manifest is not a supported complete schema: {path}")
    require(isinstance(raw, dict),
            f"exact-main raw bundle is not an object: {path}")
    _validate_exact_schema_pair(document, raw)
    require(isinstance(snapshot, dict) and set(snapshot) == {"size", "sha256"} and
            type(snapshot.get("size")) is int and snapshot["size"] > 0 and
            isinstance(snapshot.get("sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", snapshot["sha256"]) is not None,
            f"exact-main manifest snapshot identity is invalid: {path}")
    return (document, validate_evidence_scope(
        selection_scope_from_verified_bundle(document, raw)), snapshot)


def manifest_cell(cell: object) -> dict[str, Any]:
    require(isinstance(cell, dict) and set(cell) == {
        "identifier", "k", "r", "shard_bytes", "losses", "seed",
    }, "exact-main campaign cell shape differs")
    require(isinstance(cell["identifier"], str),
            "exact-main campaign identifier is invalid")
    for key in ("k", "r", "shard_bytes", "losses", "seed"):
        require(type(cell[key]) is int and cell[key] > 0,
                f"exact-main campaign {key} is invalid")
    expected = main_cell(
        f"gate-k{cell['k']}-b{cell['shard_bytes']}",
        cell["k"], cell["r"], cell["shard_bytes"], cell["losses"],
                         "exact-main-generic-gate")
    require(cell == {
        "identifier": expected["identifier"], "k": expected["K"],
        "r": expected["R"], "shard_bytes": expected["shard_bytes"],
        "losses": expected["loss_count"], "seed": expected["seed"],
    }, "exact-main campaign cell differs from its canonical gate cell")
    return expected


def refinement_cells(evaluated: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: dict[tuple[int, int], dict[str, Any]] = {}
    grouped: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for cell in evaluated:
        side = 1 << (cell["K"] - 1).bit_length()
        grouped.setdefault((side, cell["shard_bytes"]), []).append(cell)
    for (_, shard_bytes), cells in grouped.items():
        cells.sort(key=lambda item: item["K"])
        for first, second in zip(cells, cells[1:]):
            if first["passes_gate"] == second["passes_gate"] or \
                    second["K"] == first["K"] + 1:
                continue
            for k in range(first["K"] + 1, second["K"]):
                candidate = main_cell(
                    f"gate-k{k}-b{shard_bytes}", k, k, shard_bytes, k,
                    "exact-main-generic-gate")
                result[(k, shard_bytes)] = candidate
    return [result[key] for key in sorted(result)]


def derive_survivors(plan_root: Path, manifests: list[dict[str, Any]],
                     references: list[dict[str, Any]],
                     scopes: list[dict[str, Any]]) -> dict[str, Any]:
    plan = validate_plan(plan_root)
    require(len(manifests) == len(references) == len(scopes) and manifests,
            "gate manifests, references, and evidence scopes differ")
    scope = validate_evidence_scope(scopes[0])
    for candidate_scope in scopes[1:]:
        validate_evidence_scope(candidate_scope)
    require(all(candidate == scope for candidate in scopes),
            "gate manifests do not share one normalized evidence environment")
    planned = planned_gate_cells(plan_root, plan)
    evaluated = []
    seen = set()
    candidate_commit = None
    for manifest in manifests:
        campaign = manifest.get("campaign")
        identities = manifest.get("identities")
        analysis = manifest.get("analysis")
        require(isinstance(campaign, dict) and isinstance(identities, dict) and
                isinstance(analysis, dict), "exact-main manifest is incomplete")
        require(campaign.get("candidate_mode") == "generic" and
                campaign.get("batch") == 1 and campaign.get("reuse") == 8 and
                campaign.get("iterations") == 9 and campaign.get("warmup") == 2 and
                campaign.get("threads") == 1 and campaign.get("rounds") == 3 and
                campaign.get("order") == list(EXACT_ORDER),
                "exact-main gate campaign semantics differ")
        baseline = identities.get("baseline_source", {})
        candidate = identities.get("candidate_source", {})
        require(baseline.get("head") == EXACT_MAIN_COMMIT and
                isinstance(candidate.get("head"), str),
                "exact-main source identity differs")
        if candidate_commit is None:
            candidate_commit = candidate["head"]
        require(candidate["head"] == candidate_commit,
                "gate manifests use different candidate commits")
        require(scope["sources"]["candidate"]["head"] == candidate_commit and
                scope["sources"]["baseline"]["head"] == EXACT_MAIN_COMMIT,
                "normalized evidence scope names a different candidate commit")
        raw_cells = campaign.get("cells")
        require(isinstance(raw_cells, list), "gate campaign cells are absent")
        for raw_cell in raw_cells:
            cell = manifest_cell(raw_cell)
            key = (cell["K"], cell["shard_bytes"])
            require(key not in seen and 5 <= cell["K"] <= 128 and
                    cell["R"] == cell["K"] and
                    cell["loss_count"] == cell["K"] and
                    cell["shard_bytes"] in EXACT_GATE_BYTES,
                    "gate manifest cell is duplicated or outside the refinement domain")
            seen.add(key)
            metrics = analysis.get(cell["identifier"])
            require(isinstance(metrics, dict), "gate analysis cell is absent")
            lower = {}
            for metric in ("decode_first_use", "decode_reuse_amortized"):
                result = metrics.get(metric)
                require(isinstance(result, dict) and
                        isinstance(result.get("ci95_lower"), (int, float)) and
                        not isinstance(result.get("ci95_lower"), bool) and
                        result.get("ratio_orientation") ==
                            "baseline_time_over_candidate_time",
                        "gate analysis metric is invalid")
                lower[metric] = finite_number(
                    result["ci95_lower"],
                    f"gate analysis {cell['identifier']} {metric} CI lower")
            evaluated.append({
                **cell, "ci95_lower": lower,
                "passes_gate": all(value >= PROMOTION_GATE for value in lower.values()),
            })
    require(set(planned).issubset(seen), "verified manifests omit planned gate cells")
    planned_evaluated = [
        cell for cell in evaluated
        if (cell["K"], cell["shard_bytes"]) in planned
    ]
    requested_refinements = {
        (cell["K"], cell["shard_bytes"])
        for cell in refinement_cells(planned_evaluated)
    }
    supplemental = seen - set(planned)
    require(supplemental.issubset(requested_refinements),
            "verified manifests contain an unrequested supplemental gate cell")
    evaluated.sort(key=lambda item: (item["shard_bytes"], item["K"]))
    survivors = [item for item in evaluated if item["passes_gate"]]
    rejected = [item for item in evaluated if not item["passes_gate"]]
    return signed({
        "schema": SURVIVOR_SCHEMA,
        "plan_content_sha256": (load_json(plan_root / "plan.json"))["content_sha256"],
        "exact_main_commit": EXACT_MAIN_COMMIT,
        "candidate_commit": candidate_commit,
        "promotion_minimum_ci95_lower": PROMOTION_GATE,
        "evidence_scope": scope,
        "evidence_scope_sha256": canonical_sha256(scope),
        "gate_manifests": references,
        "evaluated_cells": evaluated,
        "survivor_cells": survivors,
        "rejected_cells": rejected,
        "required_refinement_cells": refinement_cells(evaluated),
    })


def select_survivors(plan_root: Path, paths: list[Path], output: Path) -> dict[str, Any]:
    require(not output.exists(), f"refusing to replace {output}")
    manifests = []
    scopes = []
    references = []
    for path in paths:
        resolved = path.resolve(strict=True)
        manifest, scope, snapshot = verify_exact_manifest(resolved)
        manifests.append(manifest)
        scopes.append(scope)
        references.append({
            "path": str(resolved), **snapshot,
            "payload_digest": manifest["digest"],
        })
    result = derive_survivors(plan_root, manifests, references, scopes)
    write_json(output, result)
    return result


def validate_survivors(plan_root: Path, path: Path,
                       manifests: list[dict[str, Any]] | None = None,
                       scopes: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    retained = load_json(path)
    value = unsigned(retained, SURVIVOR_SCHEMA, "survivor set")
    require(set(value) == {
        "schema", "plan_content_sha256", "exact_main_commit", "candidate_commit",
        "promotion_minimum_ci95_lower", "gate_manifests", "evaluated_cells",
        "survivor_cells", "rejected_cells", "required_refinement_cells",
        "evidence_scope", "evidence_scope_sha256",
    }, "survivor fields differ")
    references = value["gate_manifests"]
    if manifests is None:
        manifests = []
        scopes = []
        require(isinstance(references, list) and references,
                "survivor gate references are absent")
        for reference in references:
            require(isinstance(reference, dict) and set(reference) == {
                "path", "size", "sha256", "payload_digest",
            }, "survivor gate reference differs")
            manifest_path = Path(reference["path"])
            manifest, scope, snapshot = verify_exact_manifest(manifest_path)
            require(snapshot == {
                        "size": reference["size"],
                        "sha256": reference["sha256"],
                    }, "survivor gate manifest bytes changed")
            require(manifest["digest"] == reference["payload_digest"],
                    "survivor gate payload changed")
            manifests.append(manifest)
            scopes.append(scope)
    require(scopes is not None, "verified evidence scopes are absent")
    expected = derive_survivors(plan_root, manifests, references, scopes)
    require(canonical_bytes(retained) == canonical_bytes(expected),
            "survivor set is not the deterministic verified result")
    return value


def forced_case(k: int, shard_bytes: int, backend: str,
                control: str, candidate: str, phase: str) -> dict[str, Any]:
    return {
        "name": f"{phase}-k{k}-b{shard_bytes}-{backend}-{control[0]}v{candidate[0]}",
        "K": k, "R": k, "profile": "legacy_high_v1", "field": "gf8",
        "backend": backend, "shard_bytes": shard_bytes, "loss_count": k,
        "batch": 1, "reuse": 8, "iterations": 9, "warmup": 2,
        "seed": deterministic_seed(
            phase, k, shard_bytes, BACKENDS.index(backend),
            common.MODES.index(control), common.MODES.index(candidate)),
        "control_mode": control, "candidate_mode": candidate,
    }


def forced_backends_for_scope(scope: dict[str, Any]) -> tuple[str, ...]:
    backend = scope.get("resolved_auto_backend")
    require(backend in BACKENDS, "evidence backend has no forced-runner tier")
    return BACKENDS[:BACKENDS.index(backend) + 1]


def forced_matrix(name: str, cases: list[dict[str, Any]]) -> dict[str, Any]:
    result = {"schema": common.MATRIX_SCHEMA, "name": name, "cases": cases}
    common.normalize_matrix(result)
    return result


def rejection_cells(counts: Iterable[int]) -> list[dict[str, Any]]:
    cells = []
    for k in sorted(set(counts)):
        for shard_bytes in (4096, 65536):
            for loss in sorted({1, 4, k // 2, k - 1}):
                cells.append(main_cell(
                    f"reject-k{k}-r{k}-b{shard_bytes}-l{loss}",
                    k, k, shard_bytes, loss, "loss-rejection"))
            for r in (k - 1, k // 2):
                cells.append(main_cell(
                    f"reject-k{k}-r{r}-b{shard_bytes}-l{r}",
                    k, r, shard_bytes, r, "rate-rejection"))
    return cells


def attestation_case(identifier: str, category: str, k: int, r: int,
                     shard_bytes: int, losses: int, seed: int,
                     expects_balanced: bool) -> dict[str, Any]:
    require(re.fullmatch(r"[a-z0-9][a-z0-9-]{0,127}", identifier) is not None,
            f"invalid attestation identifier {identifier!r}")
    require(category in {
        "survivor_gate", "rejected_gate", "loss_rate_neighbor",
        "extra_aligned_confirmation", "true_tail",
    }, f"invalid attestation category {category!r}")
    require(0 < r <= k and 0 < losses <= r and shard_bytes > 0 and
            0 < seed < 1 << 64, f"invalid attestation cell {identifier}")
    cell = {
        "identifier": identifier, "category": category, "K": k, "R": r,
        "profile": "legacy_high_v1", "field": "gf8", "backend": "auto",
        "shard_bytes": shard_bytes, "loss_count": losses, "batch": 1,
        "reuse": 1, "iterations": 1, "warmup": 0, "threads": 1,
        "seed": seed,
    }
    return {
        "cell": cell,
        "benchmark_arguments": [
            "--k", str(k), "--r", str(r), "--profile", "high",
            "--field", "gf8", "--backend", "auto", "--bytes",
            str(shard_bytes), "--loss", str(losses), "--batch", "1",
            "--reuse", "1", "--iterations", "1", "--warmup", "0",
            "--threads", "1", "--seed", str(seed), "--skip-legacy",
            "--retain-samples", "--report-decode-path", "--json", "OUTPUT",
        ],
        "expected_selected_decode_path": (
            "generic" if expects_balanced else "not_generic"),
        "expected_selected_decode_rule": (
            "balanced_generic" if expects_balanced else "not_balanced_generic"),
    }


def path_attestation(survivor_signed: dict[str, Any]) -> dict[str, Any]:
    """Build the exact AUTO policy implied by the externally proven cells.

    A transform bucket, a surviving K, or an aligned byte size is not evidence
    for another cell.  Only an exact passing gate pair may select the balanced
    generic rule.  All other cells remain negative controls until a later,
    separately authenticated promotion consumes their own evidence.
    """
    survivor = unsigned(survivor_signed, SURVIVOR_SCHEMA, "survivor set")
    cases = []
    survivor_pairs = {
        (cell["K"], cell["shard_bytes"])
        for cell in survivor["survivor_cells"]
    }
    for cell in survivor["evaluated_cells"]:
        expects = (cell["K"], cell["shard_bytes"]) in survivor_pairs
        cases.append(attestation_case(
            ("survivor-" if expects else "rejected-") + cell["identifier"],
            "survivor_gate" if expects else "rejected_gate",
            cell["K"], cell["R"], cell["shard_bytes"],
            cell["loss_count"], cell["seed"], expects))

    survivor_k = sorted({cell["K"] for cell in survivor["survivor_cells"]})
    for cell in rejection_cells(survivor_k):
        cases.append(attestation_case(
            "neighbor-" + cell["identifier"], "loss_rate_neighbor",
            cell["K"], cell["R"], cell["shard_bytes"],
            cell["loss_count"], cell["seed"], False))
    for k in survivor_k:
        for shard_bytes in ALIGNED_CONFIRM_BYTES:
            if shard_bytes in EXACT_GATE_BYTES:
                continue
            cases.append(attestation_case(
                f"extra-aligned-k{k}-b{shard_bytes}",
                "extra_aligned_confirmation", k, k, shard_bytes, k,
                deterministic_seed("auto-extra-aligned", k, shard_bytes), False))
        for shard_bytes in TRUE_TAIL_BYTES:
            cases.append(attestation_case(
                f"true-tail-k{k}-b{shard_bytes}", "true_tail",
                k, k, shard_bytes, k,
                deterministic_seed("auto-true-tail", k, shard_bytes), False))

    cases.sort(key=lambda item: (
        item["cell"]["K"], item["cell"]["R"],
        item["cell"]["shard_bytes"], item["cell"]["loss_count"],
        item["cell"]["seed"], item["cell"]["identifier"]))
    identifiers = [item["cell"]["identifier"] for item in cases]
    identities = [tuple(item["cell"][key] for key in (
        "K", "R", "shard_bytes", "loss_count", "seed")) for item in cases]
    require(len(set(identifiers)) == len(identifiers),
            "attestation identifiers are duplicated")
    require(len(set(identities)) == len(identities),
            "attestation cells are duplicated")
    return signed({
        "schema": ATTESTATION_SCHEMA,
        "name": "balanced-auto-exact-cell-selected-path-attestation",
        "purpose": (
            "Authenticate AUTO path selection for exact externally proven cells "
            "and fail closed for every rejected gate, loss/rate neighbor, "
            "unproven aligned size, and true tail."),
        "candidate_commit": survivor["candidate_commit"],
        "survivor_content_sha256": survivor_signed["content_sha256"],
        "evidence_scope_sha256": survivor["evidence_scope_sha256"],
        "expected_resolved_backend": survivor["evidence_scope"][
            "resolved_auto_backend"],
        "cases": cases,
    })


def materialize_stage(plan_root: Path, survivor_signed: dict[str, Any],
                      root: Path) -> dict[str, Any]:
    require(not root.exists(), f"refusing to replace existing output {root}")
    survivor = unsigned(survivor_signed, SURVIVOR_SCHEMA, "survivor set")
    require(not survivor["required_refinement_cells"],
            "K transition refinement remains; run the listed runner_cell values first")
    root.mkdir(parents=True)
    write_json(root / "survivors.json", survivor_signed)
    survivor_cells = survivor["survivor_cells"]
    survivor_k = sorted({cell["K"] for cell in survivor_cells})
    forced_backends = forced_backends_for_scope(survivor["evidence_scope"])
    artifacts = []

    by_side: dict[int, list[dict[str, Any]]] = {}
    for cell in survivor_cells:
        by_side.setdefault(1 << (cell["K"] - 1).bit_length(), []).append(cell)
    for backend in forced_backends:
        for side, cells in sorted(by_side.items()):
            cases = [
                forced_case(cell["K"], cell["shard_bytes"], backend,
                            control, candidate, "survivor")
                for cell in cells for control, candidate in MODE_PAIRS
            ]
            path = root / "forced-survivors" / f"{backend}-t{side}.json"
            write_json(path, forced_matrix(
                f"balanced-survivors-{backend}-t{side}", cases))
            artifacts.append(artifact(path, root, len(cases), len(cases) * 12,
                                      "forced_surviving_cells"))

    for mode in CONFIRM_MODES:
        cells = [
            main_cell(f"confirm-{mode}-k{k}-b{shard_bytes}",
                      k, k, shard_bytes, k, f"aligned-confirm-{mode}")
            for k in survivor_k for shard_bytes in ALIGNED_CONFIRM_BYTES
        ]
        if cells:
            path = root / "exact-main-confirm" / f"{mode}.json"
            write_json(path, main_document(
                f"balanced-aligned-exact-main-confirm-{mode}", mode, cells,
                "Selected-mode exact-main confirmation at all aligned boundaries.",
                "Run only after the generic gate and forced survivor comparisons pass."))
            artifacts.append(artifact(path, root, len(cells), len(cells) * 12,
                                      "exact_main_aligned_confirmation"))

    rejects = rejection_cells(survivor_k)
    if rejects:
        path = root / "exact-main-auto-rejection-timing.json"
        write_json(path, main_document(
            "balanced-auto-loss-rate-rejection-timing", "auto", rejects,
            "Performance/correctness controls only; this does not attest selected path.",
            "Run after a candidate rule is frozen, alongside path-attestation.json."))
        artifacts.append(artifact(path, root, len(rejects), len(rejects) * 12,
                                  "exact_main_rejection_timing"))

    attestation = path_attestation(survivor_signed)
    attest = root / "path-attestation.json"
    write_json(attest, attestation)
    artifacts.append(artifact(
        attest, root, len(attestation["cases"]), 0,
        "same_binary_selected_path_attestation"))

    for backend in forced_backends:
        for k in survivor_k:
            cases = [
                forced_case(k, shard_bytes, backend, control, candidate, "tail")
                for shard_bytes in TRUE_TAIL_BYTES
                for control, candidate in MODE_PAIRS
            ]
            path = root / "forced-true-tails" / f"{backend}-k{k}.json"
            write_json(path, forced_matrix(
                f"balanced-true-tails-{backend}-k{k}", cases))
            artifacts.append(artifact(path, root, len(cases), len(cases) * 12,
                                      "forced_non_aligned_tail"))

    stage = signed({
        "schema": STAGE_SCHEMA,
        "name": "balanced-survivor-confirmation-stage",
        "plan_content_sha256": (load_json(plan_root / "plan.json"))["content_sha256"],
        "survivor_content_sha256": survivor_signed["content_sha256"],
        "candidate_commit": survivor["candidate_commit"],
        "evidence_scope_sha256": survivor["evidence_scope_sha256"],
        "expected_resolved_backend": survivor["evidence_scope"][
            "resolved_auto_backend"],
        "forced_confirmation_backends": list(forced_backends),
        "excluded_backends": survivor["evidence_scope"]["excluded_backends"],
        "survivor_K": survivor_k,
        "surviving_gate_cell_count": len(survivor_cells),
        "true_tail_bytes": list(TRUE_TAIL_BYTES),
        "aligned_exact_main_bytes": list(ALIGNED_CONFIRM_BYTES),
        "artifacts": artifacts,
        "timed_case_count": sum(item["case_count"] for item in artifacts
                                if item["timed_child_count"] != 0),
        "timed_child_count": sum(item["timed_child_count"] for item in artifacts),
        "path_attestation_case_count": len(attestation["cases"]),
        "path_attestation_survivor_case_count": sum(
            item["expected_selected_decode_path"] == "generic"
            for item in attestation["cases"]),
        "promotion_requires_path_attestation": bool(attestation["cases"]),
    })
    write_json(root / "stage.json", stage)
    return stage


def validate_stage(plan_root: Path, root: Path,
                   manifests: list[dict[str, Any]] | None = None,
                   scopes: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    retained_survivor = load_json(root / "survivors.json")
    survivor_path = root / "survivors.json"
    survivor = validate_survivors(plan_root, survivor_path, manifests, scopes)
    stage = unsigned(load_json(root / "stage.json"), STAGE_SCHEMA, "stage")
    require(stage["candidate_commit"] == survivor["candidate_commit"],
            "stage candidate commit differs")
    with tempfile.TemporaryDirectory(prefix="leopard2-stage-expected-") as temporary:
        expected = Path(temporary) / "stage"
        materialize_stage(plan_root, retained_survivor, expected)
        compare_trees(root, expected, "stage")
    return stage


def stable_file_identity(value: dict[str, Any], label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} identity is not an object")
    relative = value.get("relative_path")
    require(relative is None or isinstance(relative, str),
            f"{label} relative path is invalid")
    digest = value.get("sha256")
    require(isinstance(digest, str) and
            re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
            f"{label} digest is invalid")
    require(type(value.get("size")) is int and value["size"] >= 0 and
            type(value.get("mode")) is int and value["mode"] >= 0,
            f"{label} size/mode is invalid")
    return {
        "relative_path": relative, "sha256": digest,
        "size": value["size"], "mode": value["mode"],
    }


def stable_source_identity(value: dict[str, Any], candidate_commit: str) -> dict[str, Any]:
    require(isinstance(value, dict) and
            value.get("head") == candidate_commit and
            value.get("expected_commit") == candidate_commit and
            value.get("status") == "clean" and
            value.get("status_sha256") == common.EMPTY_SHA256 and
            re.fullmatch(r"[0-9a-f]{40}", str(value.get("tree"))) is not None,
            "source identity is not the clean candidate commit")
    return {
        "head": candidate_commit, "tree": value["tree"],
        "status": "clean", "status_sha256": common.EMPTY_SHA256,
    }


def stable_build_identity(value: dict[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict), "build identity is not an object")
    sources = value.get("sources")
    objects = value.get("objects")
    require(isinstance(sources, dict) and set(sources) == {
        "benchmark", "decoder", "dispatch",
    }, "build source identity set differs")
    require(isinstance(objects, dict) and set(objects) == {
        "benchmark", "decoder",
    }, "build object identity set differs")
    return {
        "cache": stable_file_identity(value.get("cache"), "CMake cache"),
        "graph": stable_file_identity(value.get("graph"), "build graph"),
        "sources": {
            key: stable_file_identity(sources[key], f"{key} source")
            for key in sorted(sources)
        },
        "objects": {
            key: stable_file_identity(objects[key], f"{key} object")
            for key in sorted(objects)
        },
        "archive": stable_file_identity(value.get("archive"), "archive"),
        "binary": stable_file_identity(value.get("binary"), "benchmark binary"),
    }


def load_exact_main_runner() -> Any:
    path = Path(__file__).resolve().parents[1] / "main_compare" / "run_abba.py"
    name = "leopard2_exact_main_build_validator"
    module = sys.modules.get(name)
    if module is not None:
        return module
    module_spec = importlib.util.spec_from_file_location(name, path)
    require(module_spec is not None and module_spec.loader is not None,
            "cannot load exact-main build validator")
    module = importlib.util.module_from_spec(module_spec)
    sys.modules[name] = module
    module_spec.loader.exec_module(module)
    return module


def canonical_candidate_build_identity(
    source_root: Path, candidate_commit: str, binary_relative: str
) -> dict[str, Any]:
    """Run the current exact-main compile/archive/link semantic validator."""
    source_root = source_root.resolve(strict=True)
    relative = Path(binary_relative)
    require(not relative.is_absolute() and ".." not in relative.parts,
            "canonical benchmark path must be source-root-relative")
    binary = (source_root / relative).resolve(strict=True)
    build_root = common.find_build_root(binary).resolve(strict=True)
    archive = (build_root / "libleopard.a").resolve(strict=True)
    runner = load_exact_main_runner()
    specification = {
        "candidate_build_dir": str(build_root),
        "candidate_executable": str(binary),
        "candidate_archive": str(archive),
        "candidate_source_root": str(source_root),
        "candidate_commit": candidate_commit,
        # Candidate validation resolves this root but does not consume baseline
        # sources.  Supplying the same clean root avoids inventing a second tree.
        "baseline_source_root": str(source_root),
    }
    try:
        provenance = runner.build_provenance(
            "candidate", specification, runner.RAW_SCHEMA_V16)
    except Exception as error:
        raise PlanError(
            f"canonical Release/AUTO build validation failed: {error}") from error
    cache = provenance.get("validated_cache")
    semantics = provenance.get("validated_compile_commands")
    require(isinstance(cache, dict) and all(
        cache.get(key) == expected
        for key, expected in REQUIRED_CANDIDATE_CACHE.items()) and
            isinstance(semantics, dict) and
            semantics.get("schema") == COMPILE_COMMANDS_SCHEMA and
            semantics.get("implementation") == "candidate" and
            semantics.get("profile") == CANDIDATE_COMPILE_PROFILE and
            semantics.get("validated_optimization") == "-O3" and
            semantics.get("validated_openmp") is True and
            isinstance(semantics.get("required_entries"), list) and
            semantics["required_entries"] and
            isinstance(provenance.get("archive_link_recipe_content"), dict),
            "canonical build provenance does not bind Release AUTO semantics")
    normalized = _normalize_bound_paths(
        provenance,
        ((str(source_root), "$SOURCE"), (str(build_root), "$BUILD")),
    )
    require(isinstance(normalized, dict), "canonical build scope is invalid")
    return {
        "schema": CANONICAL_PRODUCTION_BUILD_SCHEMA,
        "validator": CANONICAL_BUILD_VALIDATOR,
        "provenance": normalized,
        "provenance_sha256": canonical_sha256(normalized),
    }


def xorshift64(value: int) -> int:
    mask = (1 << 64) - 1
    value ^= (value << 13) & mask
    value ^= value >> 7
    value ^= (value << 17) & mask
    return value & mask


def expected_missing_indices(k: int, losses: int, seed: int) -> list[int]:
    order = list(range(k))
    state = seed ^ 0xd1b54a32d192ed03
    if state == 0:
        state = 0x9e3779b97f4a7c15
    for remaining in range(len(order), 1, -1):
        state = xorshift64(state)
        selected = state % remaining
        order[remaining - 1], order[selected] = (
            order[selected], order[remaining - 1])
    return sorted(order[:losses])


def ceil_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def coherent_production_decode_pair(path: object, rule: object) -> bool:
    return type(path) is str and type(rule) is str and \
        (path, rule) in PRODUCTION_DECODE_PATH_RULE_PAIRS


def validate_selector_pair_contract(source_root: Path) -> None:
    """Prove the Python whitelist still matches the production C++ emitters."""
    dispatch = (source_root / "Leopard2Dispatch.h").read_text(encoding="utf-8")
    decoder = (source_root / "leopard2.cpp").read_text(encoding="utf-8")
    path_names = dict(re.findall(
        r'case kDecodePath([A-Za-z0-9_]+): return "([a-z0-9_]+)";', dispatch))
    rule_names = dict(re.findall(
        r'case kDecodeRule([A-Za-z0-9_]+): return "([a-z0-9_]+)";', dispatch))
    token_pairs = set(re.findall(
        r'selection\.path = kDecodePath([A-Za-z0-9_]+);\s*'
        r'selection\.rule = kDecodeRule([A-Za-z0-9_]+);', dispatch))
    token_pairs.update(re.findall(
        r'FillTerminalDecodePathInfo\(kDecodePath([A-Za-z0-9_]+),\s*'
        r'kDecodeRule([A-Za-z0-9_]+),', decoder))
    require(token_pairs and all(path in path_names and rule in rule_names
                                for path, rule in token_pairs),
            "cannot derive production decode path/rule names")
    derived = {(path_names[path], rule_names[rule]) for path, rule in token_pairs}
    require(derived == PRODUCTION_DECODE_PATH_RULE_PAIRS,
            "production decode path/rule pair contract changed")


def expected_auto_decode_pair(case: dict[str, Any], backend: str) -> tuple[str, str]:
    """Mirror terminal-plan and AUTO selector ordering for one canonical cell.

    These attestation cells are legacy-high GF8, one-item, plan-known calls with
    no force flags and at least one missing original.  The only terminal plan is
    bounded direct repair.  The remaining prediction follows the ordered
    balanced, measured-materialized, and workspace rules in
    Leopard2Dispatch.h.  A passing gate is the *only* input allowed to satisfy
    the balanced predicate; this is intentionally narrower than an unverified
    current dispatcher rule.
    """
    cell = case["cell"]
    require(backend in AUTO_BACKENDS, "AUTO attestation backend is unknown")
    k = cell["K"]
    r = cell["R"]
    losses = cell["loss_count"]
    if 2 <= k <= DIRECT_MAX_ORIGINALS and losses <= DIRECT_MAX_LOSSES:
        return "direct", "direct"
    if case["expected_selected_decode_path"] == "generic":
        return "generic", "balanced_generic"

    padded = ceil_power_of_two(r)
    parent = ceil_power_of_two(k + padded)
    rounded = (cell["shard_bytes"] + 63) & ~63
    measured_materialized = (
        k == 224 and r == 32 and padded == 32 and parent == 256 and
        0 < losses <= 8 and rounded <= 64 * 1024 and (
            (backend in {"avx2", "avx512"} and rounded >= 24 * 1024) or
            (backend == "ssse3" and rounded >= 32 * 1024)))
    if measured_materialized:
        return "materialized", "measured_materialized"
    tiled_work_slots = 2 * padded + losses
    if tiled_work_slots < parent:
        return "tiled", "workspace_tiled"
    return "materialized", "workspace_materialized"


def expected_required_work_slots(case: dict[str, Any], path: str) -> int:
    """Mirror GetDecodePlanPathInfo workspace geometry for production paths."""
    cell = case["cell"]
    if cell["profile"] == "legacy_high_v1":
        padded = ceil_power_of_two(cell["R"])
        parent = ceil_power_of_two(cell["K"] + padded)
    else:
        require(cell["profile"] == "low_v1", "unknown decode profile")
        padded = ceil_power_of_two(cell["K"])
        parent = ceil_power_of_two(padded + cell["R"])
    if path in {"no_op", "direct"}:
        return 0
    if path in {"generic", "materialized"}:
        return parent
    require(path == "tiled", "cannot derive workspace for unknown decode path")
    # All campaign cells use legacy_high_v1.  Low-profile tiled execution uses
    # 2*P, while legacy high additionally retains one locator/output slot for
    # each requested missing original.
    if cell["profile"] == "legacy_high_v1":
        return 2 * padded + cell["loss_count"]
    return 2 * padded


def validate_attestation_output(document: object, case: dict[str, Any],
                                expected_backend: str | None = None) -> dict[str, Any]:
    require(isinstance(document, dict) and
            document.get("schema") == common.BENCHMARK_SCHEMA and
            set(document) == {
                "schema", "build", "parameters", "resolved", "correctness",
                "workload_digests", "memory", "metrics", "legacy",
            }, "attestation benchmark top-level shape differs")
    parameters = document.get("parameters")
    resolved = document.get("resolved")
    correctness = document.get("correctness")
    build_output = document.get("build")
    require(isinstance(build_output, dict) and set(build_output) == {
        "compiler", "compiler_version", "cplusplus",
    } and isinstance(build_output["compiler"], str) and
            isinstance(build_output["compiler_version"], str) and
            type(build_output["cplusplus"]) is int and build_output["cplusplus"] > 0,
            "attestation benchmark build record differs")
    require(isinstance(parameters, dict) and isinstance(resolved, dict) and
            isinstance(correctness, dict), "attestation benchmark is incomplete")
    require(set(parameters) == {
        "K", "R", "requested_profile", "requested_field", "requested_backend",
        "force_generic_decode", "force_specialized_decode", "force_tiled_decode",
        "force_materialized_decode", "skip_legacy", "retain_samples",
        "report_decode_path", "shard_bytes", "loss_count",
        "missing_original_indices", "batch", "reuse", "iterations", "warmup",
        "thread_count", "seed",
    }, "attestation benchmark parameter shape differs")
    cell = case["cell"]
    expected_parameters = {
        "K": cell["K"], "R": cell["R"],
        "requested_profile": "legacy_high_v1", "requested_field": "gf8",
        "requested_backend": "auto", "force_generic_decode": False,
        "force_specialized_decode": False, "force_tiled_decode": False,
        "force_materialized_decode": False, "skip_legacy": True,
        "retain_samples": True, "report_decode_path": True,
        "shard_bytes": cell["shard_bytes"], "loss_count": cell["loss_count"],
        "missing_original_indices": expected_missing_indices(
            cell["K"], cell["loss_count"], cell["seed"]),
        "batch": 1, "reuse": 1, "iterations": 1, "warmup": 0,
        "thread_count": 1, "seed": cell["seed"],
    }
    for key, expected in expected_parameters.items():
        actual = parameters.get(key)
        if key == "missing_original_indices":
            require(isinstance(actual, list) and
                    all(type(item) is int for item in actual) and actual == expected,
                    f"attestation benchmark {key} differs for {cell['identifier']}")
        else:
            require(type(actual) is type(expected) and actual == expected,
                    f"attestation benchmark {key} differs for {cell['identifier']}")
    transform = ceil_power_of_two(cell["R"])
    parent = ceil_power_of_two(cell["K"] + transform)
    require(set(resolved) == {
        "profile", "field", "backend", "thread_count", "parent_count",
        "padded_side", "selected_decode_path", "selected_decode_rule",
        "decode_required_work_slots", "decode_aligned_prefix_bytes",
        "decode_tail_bytes", "decode_rounded_bytes", "decode_multi_item_batch",
    } and resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and
            type(resolved.get("backend")) is str and
            resolved.get("backend") in AUTO_BACKENDS and
            (expected_backend is None or resolved.get("backend") == expected_backend) and
            type(resolved.get("thread_count")) is int and
            resolved.get("thread_count") == 1 and
            type(resolved.get("parent_count")) is int and
            resolved.get("parent_count") == parent and
            type(resolved.get("padded_side")) is int and
            resolved.get("padded_side") == transform and
            type(resolved.get("decode_aligned_prefix_bytes")) is int and
            resolved.get("decode_aligned_prefix_bytes") ==
                cell["shard_bytes"] & ~63 and
            type(resolved.get("decode_tail_bytes")) is int and
            resolved.get("decode_tail_bytes") == cell["shard_bytes"] & 63 and
            type(resolved.get("decode_rounded_bytes")) is int and
            resolved.get("decode_rounded_bytes") ==
                (cell["shard_bytes"] + 63) & ~63 and
            resolved.get("decode_multi_item_batch") is False,
            f"attestation resolved codec identity differs for {cell['identifier']}")
    selected_path = resolved["selected_decode_path"]
    selected_rule = resolved["selected_decode_rule"]
    require(coherent_production_decode_pair(selected_path, selected_rule),
            f"attestation selected path/rule pair is not a production pair: "
            f"{cell['identifier']}")
    expected_pair = expected_auto_decode_pair(case, resolved["backend"])
    require((selected_path, selected_rule) == expected_pair,
            f"attestation selected path/rule differs for {cell['identifier']}: "
            f"expected {expected_pair!r}, got {(selected_path, selected_rule)!r}")
    required_slots = expected_required_work_slots(case, selected_path)
    require(type(resolved.get("decode_required_work_slots")) is int and
            resolved["decode_required_work_slots"] == required_slots,
            f"attestation required work slots differ for {cell['identifier']}: "
            f"expected {required_slots}, got "
            f"{resolved.get('decode_required_work_slots')!r}")
    require(isinstance(correctness, dict) and set(correctness) == {
        "leopard2_round_trip", "legacy_comparison",
    } and correctness.get("leopard2_round_trip") is True and
            correctness.get("legacy_comparison") is None,
            f"attestation round trip differs for {cell['identifier']}")
    digests = document.get("workload_digests")
    require(isinstance(digests, dict) and set(digests) == {
        "algorithm", "original_data", "transmitted_parity", "recovered_originals",
    } and digests.get("algorithm") == "fnv1a64" and all(
        type(digests.get(key)) is str and
        re.fullmatch(r"[0-9a-f]{16}", digests[key]) is not None
        for key in ("original_data", "transmitted_parity", "recovered_originals")),
        f"attestation workload digests differ for {cell['identifier']}")
    memory = document.get("memory")
    require(isinstance(memory, dict) and set(memory) == {
        "scratch_alignment", "encode_scratch_bytes_per_stripe",
        "decode_scratch_bytes_per_stripe", "encode_scratch_bytes_batch",
        "decode_scratch_bytes_batch",
    } and all(type(value) is int and value >= 0 for value in memory.values()) and
            memory["scratch_alignment"] > 0 and
            memory["encode_scratch_bytes_per_stripe"] ==
                memory["encode_scratch_bytes_batch"] and
            memory["decode_scratch_bytes_per_stripe"] ==
                memory["decode_scratch_bytes_batch"],
            f"attestation memory record differs for {cell['identifier']}")
    metrics = document.get("metrics")
    require(isinstance(metrics, dict) and set(metrics) == {
        "codec_setup", "encode_execution", "decode_plan_setup",
        "decode_execution", "decode_amortized_at_reuse", "rate_semantics",
    }, f"attestation metric shape differs for {cell['identifier']}")
    for label, key in (("codec setup", "codec_setup"),
                       ("decode plan setup", "decode_plan_setup")):
        summary = metrics[key]
        require(isinstance(summary, dict) and
                isinstance(summary.get("samples_us"), list) and
                len(summary["samples_us"]) == 1 and
                isinstance(summary["samples_us"][0], (int, float)) and
                not isinstance(summary["samples_us"][0], bool) and
                summary["samples_us"][0] >= 0,
                f"attestation {label} sample differs for {cell['identifier']}")
    for label, key in (("encode", "encode_execution"),
                       ("decode", "decode_execution")):
        summary = metrics[key]
        require(isinstance(summary, dict) and
                isinstance(summary.get("samples_us_per_batch_call"), list) and
                len(summary["samples_us_per_batch_call"]) == 1 and
                isinstance(summary["samples_us_per_batch_call"][0], (int, float)) and
                not isinstance(summary["samples_us_per_batch_call"][0], bool) and
                summary["samples_us_per_batch_call"][0] > 0,
                f"attestation {label} sample differs for {cell['identifier']}")
    amortized = metrics["decode_amortized_at_reuse"]
    require(isinstance(amortized, dict) and
            type(amortized.get("reuse_count")) is int and
            amortized.get("reuse_count") == 1 and
            isinstance(metrics["rate_semantics"], str),
            f"attestation amortization record differs for {cell['identifier']}")
    require(document.get("legacy") == {
        "available": False, "unavailable_reason": "disabled by --skip-legacy",
        "codec_setup": None, "decode_timing_includes_setup": True,
        "encode_execution": None, "decode_including_setup": None,
    }, f"attestation legacy-skip record differs for {cell['identifier']}")
    return {
        "resolved_backend": resolved["backend"],
        "selected_decode_path": selected_path,
        "selected_decode_rule": selected_rule,
        "workload_digests_sha256": canonical_sha256(digests),
        "benchmark_payload_sha256": canonical_sha256(document),
    }


def attestation_identities(source_root: Path, candidate_commit: str,
                           binary_relative: str) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    validate_selector_pair_contract(source_root)
    source = stable_source_identity(
        common.git_identity(source_root, candidate_commit), candidate_commit)
    build = {
        "artifact_closure": stable_build_identity(
            common.build_identity(source_root, binary_relative)),
        "canonical_production": canonical_candidate_build_identity(
            source_root, candidate_commit, binary_relative),
    }
    collector = stable_file_identity(
        common.file_identity(Path(__file__), source_root, "attestation collector"),
        "attestation collector")
    require(collector["relative_path"] ==
            "experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py",
            "attestation collector path is not canonical")
    return source, build, collector


def load_balanced_analyzer() -> Any:
    path = Path(__file__).resolve().parents[3] / BALANCED_ANALYZER_RELATIVE_PATH
    specification = importlib.util.spec_from_file_location(
        "leopard2_balanced_promotion_analyzer", path)
    require(specification is not None and specification.loader is not None,
            "cannot load balanced evidence analyzer")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def _timing_manifest_reference(
    path: Path, snapshot: bytes, payload_digest: str,
) -> dict[str, Any]:
    require(re.fullmatch(r"[0-9a-f]{64}", payload_digest) is not None,
            "timing evidence payload digest is invalid")
    return {
        "path": str(path), "size": len(snapshot),
        "sha256": hashlib.sha256(snapshot).hexdigest(),
        "payload_digest": payload_digest,
    }


def _exact_campaign_matches(
    document: dict[str, Any], campaign: Mapping[str, Any],
) -> bool:
    expected_cells = [{
        "identifier": cell["identifier"], "k": cell["K"], "r": cell["R"],
        "shard_bytes": cell["shard_bytes"], "losses": cell["loss_count"],
        "seed": cell["seed"],
    } for cell in document["cells"]]
    return (
        campaign.get("candidate_mode") == document["candidate_mode"] and
        campaign.get("batch") == 1 and
        campaign.get("reuse") == document["reuse"] and
        campaign.get("iterations") == document["iterations"] and
        campaign.get("warmup") == document["warmup"] and
        campaign.get("threads") == 1 and campaign.get("rounds") == 3 and
        campaign.get("order") == list(EXACT_ORDER) and
        campaign.get("cells") == expected_cells
    )


def derive_promotion_timing_evidence(
    stage_root: Path,
    forced_manifest_paths: list[Path],
    exact_manifest_paths: list[Path],
    source: Mapping[str, Any],
    build: Mapping[str, Any],
) -> dict[str, Any]:
    """Authenticate every timing artifact required by one promotion stage."""
    stage = unsigned(load_json(stage_root / "stage.json"), STAGE_SCHEMA, "stage")
    survivor = unsigned(
        load_json(stage_root / "survivors.json"), SURVIVOR_SCHEMA,
        "survivor set")
    require(source.get("head") == stage["candidate_commit"] and
            source.get("tree") == survivor["evidence_scope"]["sources"][
                "candidate"]["tree"],
            "promotion timing source differs from the attested candidate")
    artifact_closure = build.get("artifact_closure")
    require(isinstance(artifact_closure, Mapping),
            "promotion timing build closure is absent")

    forced_expected: dict[str, dict[str, Any]] = {}
    exact_expected: list[tuple[dict[str, Any], dict[str, Any]]] = []
    for artifact_value in stage["artifacts"]:
        kind = artifact_value["kind"]
        if kind in {"forced_surviving_cells", "forced_non_aligned_tail"}:
            require(artifact_value["sha256"] not in forced_expected,
                    "forced timing stage artifact digest is duplicated")
            forced_expected[artifact_value["sha256"]] = artifact_value
        elif kind in {
            "exact_main_aligned_confirmation",
            "exact_main_rejection_timing",
        }:
            document = validate_main_document(
                load_json(stage_root / artifact_value["path"]),
                artifact_value["path"])
            exact_expected.append((artifact_value, document))

    require(len(forced_manifest_paths) == len(forced_expected) and
            len(exact_manifest_paths) == len(exact_expected),
            "promotion timing manifest count does not cover the stage")
    analyzer = load_balanced_analyzer()
    forced_records: list[dict[str, Any]] = []
    seen_forced: set[str] = set()
    for supplied in forced_manifest_paths:
        path = supplied.absolute()
        before = stable_regular_file_bytes(
            path, MAX_JSON_DOCUMENT_BYTES, f"forced timing manifest {path}")
        try:
            summary = analyzer.analyze(path)
        except Exception as error:
            raise PlanError(
                f"forced timing analyzer rejected {path}: {error}") from error
        after = stable_regular_file_bytes(
            path, MAX_JSON_DOCUMENT_BYTES, f"forced timing manifest {path}")
        require(before == after and isinstance(summary, dict) and
                summary.get("schema") == common.SUMMARY_SCHEMA,
                f"forced timing manifest changed or has wrong schema: {path}")
        summary_unsigned = dict(summary)
        summary_digest = summary_unsigned.pop("content_sha256", None)
        require(isinstance(summary_digest, str) and
                canonical_sha256(summary_unsigned) == summary_digest and
                summary.get("source_manifest_sha256") ==
                    hashlib.sha256(before).hexdigest(),
                f"forced timing summary identity differs: {path}")
        matrix = summary.get("matrix")
        require(isinstance(matrix, dict) and
                type(matrix.get("size")) is int and
                isinstance(matrix.get("sha256"), str),
                f"forced timing matrix identity differs: {path}")
        expected = forced_expected.get(matrix["sha256"])
        require(expected is not None and
                matrix["size"] == expected["size"] and
                matrix["sha256"] not in seen_forced,
                f"forced timing matrix is missing, duplicated, or extra: {path}")
        seen_forced.add(matrix["sha256"])
        summary_source = summary.get("source")
        summary_build = summary.get("build")
        require(isinstance(summary_source, Mapping) and
                summary_source.get("head") == source["head"] and
                summary_source.get("tree") == source["tree"] and
                summary_source.get("status") == "clean" and
                isinstance(summary_build, Mapping) and all(
                    isinstance(summary_build.get(name), Mapping) and
                    summary_build[name].get("sha256") ==
                        artifact_closure[name]["sha256"]
                    for name in ("archive", "binary")),
                f"forced timing source/build differs from attestation: {path}")
        cells = summary.get("cells")
        require(isinstance(cells, list) and
                len(cells) == expected["case_count"],
                f"forced timing cell set differs: {path}")
        maximum_upper = -math.inf
        for cell_index, cell in enumerate(cells):
            require(isinstance(cell, Mapping) and
                    isinstance(cell.get("decode_execution"), Mapping),
                    f"forced timing cell {cell_index} differs: {path}")
            decode = cell["decode_execution"]
            finite_number(
                decode.get("geometric_control_over_candidate"),
                f"forced timing {path} cell {cell_index} decode ratio",
                positive=True)
            upper_percent = finite_real(
                decode.get("ci95_upper_percent"),
                f"forced timing {path} cell {cell_index} upper CI percent")
            upper = 1.0 + upper_percent / 100.0
            require(math.isfinite(upper) and upper > 0,
                    f"forced timing {path} cell {cell_index} upper CI is invalid")
            maximum_upper = max(maximum_upper, upper)
        if expected["kind"] == "forced_surviving_cells":
            require(maximum_upper <= FORCED_SURVIVOR_MAXIMUM_RATIO,
                    f"forced survivor generic path exceeds the 2% "
                    f"same-binary regression limit: {path}")
        forced_records.append({
            "stage_artifact": {
                "path": expected["path"], "kind": expected["kind"],
                "size": expected["size"], "sha256": expected["sha256"],
            },
            "manifest": _timing_manifest_reference(
                path, before, hashlib.sha256(before).hexdigest()),
            "summary_content_sha256": summary_digest,
            "maximum_decode_ci95_upper": maximum_upper,
        })
    require(seen_forced == set(forced_expected),
            "forced timing evidence does not cover every stage artifact")

    exact_records: list[dict[str, Any]] = []
    seen_exact: set[str] = set()
    minimum_lower = math.inf
    for supplied in exact_manifest_paths:
        path = supplied.absolute()
        snapshot_bytes = stable_regular_file_bytes(
            path, MAX_JSON_DOCUMENT_BYTES, f"exact timing manifest {path}")
        manifest, scope, snapshot = verify_exact_manifest(path)
        require(snapshot == {
                    "size": len(snapshot_bytes),
                    "sha256": hashlib.sha256(snapshot_bytes).hexdigest(),
                } and scope == survivor["evidence_scope"],
                f"exact timing scope/snapshot differs: {path}")
        campaign = manifest.get("campaign")
        require(isinstance(campaign, Mapping),
                f"exact timing campaign is absent: {path}")
        matches = [
            (artifact_value, document)
            for artifact_value, document in exact_expected
            if artifact_value["sha256"] not in seen_exact and
               _exact_campaign_matches(document, campaign)
        ]
        require(len(matches) == 1,
                f"exact timing manifest is missing, duplicated, or extra: {path}")
        artifact_value, document = matches[0]
        seen_exact.add(artifact_value["sha256"])
        identities = manifest.get("identities")
        require(isinstance(identities, Mapping) and
                identities.get("candidate_source", {}).get("head") ==
                    source["head"] and all(
                    isinstance(identities.get(name), Mapping) and
                    identities[name].get("sha256") ==
                        artifact_closure[closure_name]["sha256"]
                    for name, closure_name in (
                        ("candidate_archive", "archive"),
                        ("candidate_executable", "binary"))),
                f"exact timing candidate source/build differs: {path}")
        analysis = manifest.get("analysis")
        require(isinstance(analysis, Mapping) and
                set(analysis) == {
                    cell["identifier"] for cell in document["cells"]},
                f"exact timing analysis cell set differs: {path}")
        artifact_minimum = math.inf
        for cell in document["cells"]:
            metrics = analysis[cell["identifier"]]
            require(isinstance(metrics, Mapping),
                    f"exact timing metrics are absent: {path}")
            for metric in ("decode_first_use", "decode_reuse_amortized"):
                value = metrics.get(metric)
                require(isinstance(value, Mapping) and
                        value.get("ratio_orientation") ==
                            "baseline_time_over_candidate_time",
                        f"exact timing metric orientation differs: {path}")
                lower = finite_number(
                    value.get("ci95_lower"),
                    f"exact timing {cell['identifier']} {metric} CI lower",
                    positive=True)
                artifact_minimum = min(artifact_minimum, lower)
        require(artifact_minimum >= 1.0 - NEIGHBOR_MAXIMUM_REGRESSION,
                f"exact timing exceeds the 2% neighbor regression limit: {path}")
        minimum_lower = min(minimum_lower, artifact_minimum)
        payload_digest = manifest.get("digest")
        require(isinstance(payload_digest, str),
                f"exact timing payload digest is absent: {path}")
        exact_records.append({
            "stage_artifact": {
                "path": artifact_value["path"], "kind": artifact_value["kind"],
                "size": artifact_value["size"],
                "sha256": artifact_value["sha256"],
            },
            "manifest": _timing_manifest_reference(
                path, snapshot_bytes, payload_digest),
            "minimum_decode_ci95_lower": artifact_minimum,
        })
    require(seen_exact == {
                artifact_value["sha256"]
                for artifact_value, _document in exact_expected},
            "exact timing evidence does not cover every stage artifact")
    survivor_forced_uppers = [
        record["maximum_decode_ci95_upper"] for record in forced_records
        if record["stage_artifact"]["kind"] == "forced_surviving_cells"
    ]
    require(survivor_forced_uppers,
            "promotion stage has no same-binary survivor timing evidence")
    return signed({
        "schema": PROMOTION_TIMING_SCHEMA,
        "candidate_commit": stage["candidate_commit"],
        "stage_content_sha256": (
            load_json(stage_root / "stage.json"))["content_sha256"],
        "neighbor_maximum_regression": NEIGHBOR_MAXIMUM_REGRESSION,
        "maximum_forced_survivor_decode_ci95_upper":
            max(survivor_forced_uppers),
        "minimum_exact_decode_ci95_lower": minimum_lower,
        "forced_results": sorted(
            forced_records, key=lambda item: item["stage_artifact"]["path"]),
        "exact_main_results": sorted(
            exact_records, key=lambda item: item["stage_artifact"]["path"]),
    })


def validate_promotion_timing_shape(
    stage_root: Path,
    retained: object,
) -> dict[str, Any]:
    value = unsigned(retained, PROMOTION_TIMING_SCHEMA, "promotion timing")
    require(set(value) == {
        "schema", "candidate_commit", "stage_content_sha256",
        "neighbor_maximum_regression",
        "maximum_forced_survivor_decode_ci95_upper",
        "minimum_exact_decode_ci95_lower",
        "forced_results", "exact_main_results",
    }, "promotion timing fields differ")
    stage_signed = load_json(stage_root / "stage.json")
    stage = unsigned(stage_signed, STAGE_SCHEMA, "stage")
    require(value["candidate_commit"] == stage["candidate_commit"] and
            value["stage_content_sha256"] == stage_signed["content_sha256"] and
            type(value["neighbor_maximum_regression"]) is float and
            value["neighbor_maximum_regression"] ==
                NEIGHBOR_MAXIMUM_REGRESSION,
            "promotion timing stage/threshold identity differs")
    forced = value.get("forced_results")
    exact = value.get("exact_main_results")
    require(isinstance(forced, list) and isinstance(exact, list),
            "promotion timing result lists are invalid")
    expected = {
        (item["path"], item["kind"], item["size"], item["sha256"])
        for item in stage["artifacts"]
        if item["kind"] in {
            "forced_surviving_cells", "forced_non_aligned_tail",
            "exact_main_aligned_confirmation",
            "exact_main_rejection_timing",
        }
    }
    observed: set[tuple[object, ...]] = set()
    exact_lowers: list[float] = []
    forced_survivor_uppers: list[float] = []
    for label, records, exact_records in (
        ("forced", forced, False), ("exact", exact, True),
    ):
        for index, item in enumerate(records):
            required = {
                "stage_artifact", "manifest",
                ("minimum_decode_ci95_lower" if exact_records else
                 "summary_content_sha256"),
            }
            if not exact_records:
                required.add("maximum_decode_ci95_upper")
            require(isinstance(item, Mapping) and set(item) == required and
                    isinstance(item.get("stage_artifact"), Mapping) and
                    set(item["stage_artifact"]) ==
                        {"path", "kind", "size", "sha256"} and
                    isinstance(item.get("manifest"), Mapping) and
                    set(item["manifest"]) ==
                        {"path", "size", "sha256", "payload_digest"},
                    f"promotion {label} timing record {index} differs")
            artifact_key = tuple(item["stage_artifact"][name] for name in (
                "path", "kind", "size", "sha256"))
            require(artifact_key in expected and artifact_key not in observed,
                    f"promotion {label} timing artifact is extra or duplicated")
            expected_kinds = (
                {"exact_main_aligned_confirmation",
                 "exact_main_rejection_timing"} if exact_records else
                {"forced_surviving_cells", "forced_non_aligned_tail"})
            require(item["stage_artifact"]["kind"] in expected_kinds,
                    f"promotion {label} timing artifact has the wrong kind")
            observed.add(artifact_key)
            reference = item["manifest"]
            require(isinstance(reference["path"], str) and
                    Path(reference["path"]).is_absolute() and
                    type(reference["size"]) is int and reference["size"] > 0 and
                    all(isinstance(reference[name], str) and
                        re.fullmatch(r"[0-9a-f]{64}", reference[name]) is not None
                        for name in ("sha256", "payload_digest")),
                    f"promotion {label} timing manifest reference differs")
            if exact_records:
                lower = finite_number(
                    item["minimum_decode_ci95_lower"],
                    f"promotion exact timing record {index} minimum",
                    positive=True)
                require(lower >= 0.98,
                        "promotion exact timing exceeds the 2% regression limit")
                exact_lowers.append(lower)
            else:
                require(isinstance(item["summary_content_sha256"], str) and
                        re.fullmatch(
                            r"[0-9a-f]{64}",
                            item["summary_content_sha256"]) is not None,
                        "promotion forced timing summary digest differs")
                upper = finite_number(
                    item["maximum_decode_ci95_upper"],
                    f"promotion forced timing record {index} maximum",
                    positive=True)
                if item["stage_artifact"]["kind"] == \
                        "forced_surviving_cells":
                    require(upper <= FORCED_SURVIVOR_MAXIMUM_RATIO,
                            "promotion forced survivor timing exceeds the 2% "
                            "same-binary regression limit")
                    forced_survivor_uppers.append(upper)
    require(observed == expected and exact_lowers and
            forced_survivor_uppers and
            finite_number(
                value["maximum_forced_survivor_decode_ci95_upper"],
                "promotion maximum forced survivor timing", positive=True) ==
                max(forced_survivor_uppers) and
            finite_number(
                value["minimum_exact_decode_ci95_lower"],
                "promotion minimum exact timing", positive=True) ==
                min(exact_lowers),
            "promotion timing evidence does not exactly cover the stage")
    return value


def validate_promotion_timing_evidence(
    stage_root: Path,
    retained: object,
    source: Mapping[str, Any],
    build: Mapping[str, Any],
) -> dict[str, Any]:
    value = validate_promotion_timing_shape(stage_root, retained)
    forced = value["forced_results"]
    exact = value["exact_main_results"]
    forced_paths = [
        Path(item["manifest"]["path"]) for item in forced
        if isinstance(item, Mapping) and isinstance(item.get("manifest"), Mapping)
    ]
    exact_paths = [
        Path(item["manifest"]["path"]) for item in exact
        if isinstance(item, Mapping) and isinstance(item.get("manifest"), Mapping)
    ]
    require(len(forced_paths) == len(forced) and len(exact_paths) == len(exact),
            "promotion timing manifest references are malformed")
    expected = derive_promotion_timing_evidence(
        stage_root, forced_paths, exact_paths, source, build)
    require(canonical_bytes(retained) == canonical_bytes(expected),
            "promotion timing evidence is not the deterministic verified result")
    return value


def _attestation_result_build_contract(
    result_schema: str,
) -> tuple[str, str, str, str, Mapping[str, str], bool]:
    """Select one exact canonical-build generation for an outer result."""
    require(result_schema in {
                ATTESTATION_RESULT_SCHEMA_V5,
                ATTESTATION_RESULT_SCHEMA_V6,
                ATTESTATION_RESULT_SCHEMA_V7,
                ATTESTATION_RESULT_SCHEMA_V8,
                ATTESTATION_RESULT_SCHEMA_V9,
                ATTESTATION_RESULT_SCHEMA_V10,
                ATTESTATION_RESULT_SCHEMA_V11,
                ATTESTATION_RESULT_SCHEMA_V12,
                ATTESTATION_RESULT_SCHEMA_V13,
                ATTESTATION_RESULT_SCHEMA},
            "attestation result schema is unsupported")
    if result_schema == ATTESTATION_RESULT_SCHEMA:
        return (
            CANONICAL_PRODUCTION_BUILD_SCHEMA,
            CANONICAL_BUILD_VALIDATOR,
            COMPILE_COMMANDS_SCHEMA,
            CANDIDATE_COMPILE_PROFILE,
            REQUIRED_CANDIDATE_CACHE,
            True,
        )
    if result_schema == ATTESTATION_RESULT_SCHEMA_V13:
        return (
            CANONICAL_PRODUCTION_BUILD_SCHEMA_V9,
            CANONICAL_BUILD_VALIDATOR_V9,
            COMPILE_COMMANDS_SCHEMA_V11,
            CANDIDATE_COMPILE_PROFILE_V7,
            REQUIRED_CANDIDATE_CACHE_V9,
            True,
        )
    if result_schema == ATTESTATION_RESULT_SCHEMA_V12:
        return (
            CANONICAL_PRODUCTION_BUILD_SCHEMA_V8,
            CANONICAL_BUILD_VALIDATOR_V8,
            COMPILE_COMMANDS_SCHEMA_V10,
            CANDIDATE_COMPILE_PROFILE_V6,
            REQUIRED_CANDIDATE_CACHE_V8,
            True,
        )
    if result_schema == ATTESTATION_RESULT_SCHEMA_V11:
        return (
            CANONICAL_PRODUCTION_BUILD_SCHEMA_V7,
            CANONICAL_BUILD_VALIDATOR_V7,
            COMPILE_COMMANDS_SCHEMA_V9,
            CANDIDATE_COMPILE_PROFILE_V5,
            REQUIRED_CANDIDATE_CACHE_V6,
            True,
        )
    if result_schema == ATTESTATION_RESULT_SCHEMA_V10:
        return (
            CANONICAL_PRODUCTION_BUILD_SCHEMA_V6,
            CANONICAL_BUILD_VALIDATOR_V6,
            COMPILE_COMMANDS_SCHEMA_V8,
            CANDIDATE_COMPILE_PROFILE_V4,
            REQUIRED_CANDIDATE_CACHE_V6,
            True,
        )
    if result_schema == ATTESTATION_RESULT_SCHEMA_V9:
        return (
            CANONICAL_PRODUCTION_BUILD_SCHEMA_V5,
            CANONICAL_BUILD_VALIDATOR_V5,
            COMPILE_COMMANDS_SCHEMA_V7,
            CANDIDATE_COMPILE_PROFILE_V3,
            REQUIRED_CANDIDATE_CACHE_V5,
            True,
        )
    if result_schema == ATTESTATION_RESULT_SCHEMA_V8:
        return (
            CANONICAL_PRODUCTION_BUILD_SCHEMA_V4,
            CANONICAL_BUILD_VALIDATOR_V4,
            COMPILE_COMMANDS_SCHEMA_V6,
            CANDIDATE_COMPILE_PROFILE_V3,
            REQUIRED_CANDIDATE_CACHE_V4,
            True,
        )
    if result_schema == ATTESTATION_RESULT_SCHEMA_V7:
        return (
            CANONICAL_PRODUCTION_BUILD_SCHEMA_V3,
            CANONICAL_BUILD_VALIDATOR_V3,
            COMPILE_COMMANDS_SCHEMA_V5,
            CANDIDATE_COMPILE_PROFILE_V2,
            REQUIRED_CANDIDATE_CACHE_V3,
            True,
        )
    return (
        CANONICAL_PRODUCTION_BUILD_SCHEMA_V2,
        CANONICAL_BUILD_VALIDATOR_V2,
        COMPILE_COMMANDS_SCHEMA_V4,
        CANDIDATE_COMPILE_PROFILE_V2,
        REQUIRED_CANDIDATE_CACHE_V2,
        result_schema == ATTESTATION_RESULT_SCHEMA_V6,
    )


def derive_attestation_result(
    stage_root: Path, source: dict[str, Any], build: dict[str, Any],
    collector: dict[str, Any], raw_documents: dict[str, object],
    raw_artifacts: dict[str, dict[str, Any]],
    timing_evidence: dict[str, Any] | None,
    *, result_schema: str = ATTESTATION_RESULT_SCHEMA,
) -> dict[str, Any]:
    stage_signed = load_json(stage_root / "stage.json")
    stage = unsigned(stage_signed, STAGE_SCHEMA, "stage")
    spec_signed = load_json(stage_root / "path-attestation.json")
    spec = unsigned(spec_signed, ATTESTATION_SCHEMA, "path attestation")
    survivor_signed = load_json(stage_root / "survivors.json")
    survivor = unsigned(survivor_signed, SURVIVOR_SCHEMA, "survivor set")
    require(spec["candidate_commit"] == stage["candidate_commit"] and
            spec["survivor_content_sha256"] == stage["survivor_content_sha256"] and
            survivor_signed["content_sha256"] == stage["survivor_content_sha256"] and
            spec["evidence_scope_sha256"] == survivor["evidence_scope_sha256"] and
            stage["evidence_scope_sha256"] == survivor["evidence_scope_sha256"] and
            spec["expected_resolved_backend"] ==
                survivor["evidence_scope"]["resolved_auto_backend"] and
            stage["expected_resolved_backend"] == spec["expected_resolved_backend"] and
            spec["expected_resolved_backend"] in CAMPAIGN_BACKENDS,
            "attestation specification is stale")
    require(isinstance(source, dict) and set(source) == {
        "head", "tree", "status", "status_sha256",
    } and source.get("head") == stage["candidate_commit"] and
            re.fullmatch(r"[0-9a-f]{40}", str(source.get("tree"))) is not None and
            source.get("status") == "clean" and
            source.get("status_sha256") == common.EMPTY_SHA256,
            "attestation source is not the exact clean candidate commit")
    (canonical_schema, canonical_validator, compile_schema, compile_profile,
     required_candidate_cache, require_timing) = \
        _attestation_result_build_contract(result_schema)
    if require_timing:
        require(timing_evidence is not None,
                "attestation result requires promotion timing evidence")
        validate_promotion_timing_shape(stage_root, timing_evidence)
    else:
        require(timing_evidence is None,
                "historical attestation result cannot contain timing evidence")
    require(isinstance(build, dict) and set(build) == {
        "artifact_closure", "canonical_production",
    }, "attestation benchmark build identity shape differs")
    artifacts = build.get("artifact_closure")
    canonical = build.get("canonical_production")
    require(isinstance(artifacts, dict) and set(artifacts) == {
        "cache", "graph", "sources", "objects", "archive", "binary",
    } and isinstance(artifacts.get("sources"), dict) and
            set(artifacts["sources"]) == {"benchmark", "decoder", "dispatch"} and
            isinstance(artifacts.get("objects"), dict) and
            set(artifacts["objects"]) == {"benchmark", "decoder"} and
            isinstance(canonical, dict) and set(canonical) == {
                "schema", "validator", "provenance", "provenance_sha256",
            } and canonical["schema"] == canonical_schema and
            canonical["validator"] == canonical_validator and
            isinstance(canonical["provenance"], dict) and
            canonical_sha256(canonical["provenance"]) ==
                canonical["provenance_sha256"] and
            isinstance(canonical["provenance"].get("validated_cache"), dict) and
            all(canonical["provenance"]["validated_cache"].get(key) == expected
                for key, expected in required_candidate_cache.items()) and
            isinstance(canonical["provenance"].get(
                "validated_compile_commands"), dict) and
            canonical["provenance"]["validated_compile_commands"].get(
                "schema") == compile_schema and
            canonical["provenance"]["validated_compile_commands"].get(
                "implementation") == "candidate" and
            canonical["provenance"]["validated_compile_commands"].get(
                "profile") == compile_profile and
            canonical["provenance"]["validated_compile_commands"].get(
                "validated_optimization") == "-O3" and
            canonical["provenance"]["validated_compile_commands"].get(
                "validated_openmp") is True and
            isinstance(canonical["provenance"][
                "validated_compile_commands"].get("required_entries"), list) and
            canonical["provenance"]["validated_compile_commands"][
                "required_entries"] and
            isinstance(canonical["provenance"].get(
                "archive_link_recipe_content"), dict),
            "attestation benchmark build identity shape differs")
    _validate_normalized_benchmark_attestation(
        canonical["provenance"]["validated_compile_commands"],
        source["head"], source["tree"], "$BUILD", "$SOURCE")
    _validate_normalized_build_configuration(
        canonical["provenance"]["validated_compile_commands"],
        canonical["provenance"]["validated_cache"], "$BUILD", "$SOURCE")
    for label, identity in [
        ("CMake cache", artifacts["cache"]),
        ("build graph", artifacts["graph"]),
        ("archive", artifacts["archive"]),
        ("benchmark binary", artifacts["binary"]),
        *[(f"{key} source", value) for key, value in artifacts["sources"].items()],
        *[(f"{key} object", value) for key, value in artifacts["objects"].items()],
        ("collector", collector),
    ]:
        require(isinstance(identity, dict) and set(identity) == {
            "relative_path", "sha256", "size", "mode",
        } and (identity["relative_path"] is None or
               isinstance(identity["relative_path"], str)) and
                isinstance(identity["sha256"], str) and
                re.fullmatch(r"[0-9a-f]{64}", identity["sha256"]) is not None and
                type(identity["size"]) is int and identity["size"] >= 0 and
                type(identity["mode"]) is int and identity["mode"] >= 0,
                f"attestation {label} identity differs")
    require(collector["relative_path"] ==
            "experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py",
            "attestation collector identity is not canonical")
    expected_ids = [case["cell"]["identifier"] for case in spec["cases"]]
    require(len(expected_ids) == len(set(expected_ids)) and
            set(raw_documents) == set(expected_ids) and
            set(raw_artifacts) == set(expected_ids),
            "attestation result set is missing, duplicated, or extra")
    records = []
    resolved_backend = spec["expected_resolved_backend"]
    for case in spec["cases"]:
        identifier = case["cell"]["identifier"]
        observed = validate_attestation_output(
            raw_documents[identifier], case, resolved_backend)
        artifact_value = raw_artifacts[identifier]
        require(isinstance(artifact_value, dict) and set(artifact_value) == {
            "path", "size", "sha256",
        } and isinstance(artifact_value["path"], str) and
                type(artifact_value["size"]) is int and artifact_value["size"] > 0 and
                isinstance(artifact_value["sha256"], str) and
                re.fullmatch(r"[0-9a-f]{64}",
                             artifact_value["sha256"]) is not None,
                f"attestation raw artifact differs for {identifier}")
        records.append({
            "identifier": identifier, "raw": artifact_value, **observed,
        })
    payload = {
        "schema": result_schema, "status": "complete", "valid": True,
        "plan_content_sha256": stage["plan_content_sha256"],
        "stage_content_sha256": stage_signed["content_sha256"],
        "survivor_content_sha256": stage["survivor_content_sha256"],
        "attestation_content_sha256": spec_signed["content_sha256"],
        "candidate_commit": stage["candidate_commit"],
        "resolved_backend": resolved_backend,
        "source_identity": source, "build_identity": build,
        "collector_identity": collector,
        "records": records,
    }
    if require_timing:
        payload["promotion_timing_evidence"] = timing_evidence
    return signed(payload)


def output_is_ignored(path: Path, source_root: Path) -> None:
    try:
        path.resolve().relative_to(source_root.resolve())
    except ValueError:
        return
    try:
        completed = load_exact_main_runner().run_process_bounded(
            ["git", "-C", str(source_root), "check-ignore", "-q", str(path)],
            timeout=30.0, max_stdout=64 * 1024, max_stderr=64 * 1024)
    except Exception as error:
        raise PlanError(
            f"cannot verify whether attestation output is ignored: {error}") \
            from error
    require(completed.returncode == 0,
            "attestation output inside source root must be ignored by Git")


def collect_attestation(plan_root: Path, stage_root: Path, source_root: Path,
                        binary_relative: str, output: Path,
                        forced_manifest_paths: list[Path],
                        exact_manifest_paths: list[Path]) -> dict[str, Any]:
    stage = validate_stage(plan_root, stage_root)
    candidate_commit = stage["candidate_commit"]
    relative = Path(binary_relative)
    require(not relative.is_absolute() and ".." not in relative.parts,
            "attestation benchmark path must be source-root-relative")
    stable_source_identity(
        common.git_identity(source_root, candidate_commit), candidate_commit)
    require(not output.exists(), f"refusing to replace attestation output {output}")
    output_is_ignored(output, source_root)
    binary = source_root / relative
    try:
        build_root = common.find_build_root(binary)
        refresh = load_exact_main_runner().run_process_bounded(
            ["cmake", "--build", str(build_root), "--target",
             "bench_leopard2", "-j", "1"],
            timeout=600.0, max_stdout=8 * 1024 * 1024,
            max_stderr=8 * 1024 * 1024)
    except Exception as error:
        raise PlanError(f"attestation benchmark refresh failed: {error}") from error
    require(refresh.returncode == 0,
            "attestation benchmark refresh returned failure")
    source, build, collector = attestation_identities(
        source_root, candidate_commit, binary_relative)
    timing_evidence = derive_promotion_timing_evidence(
        stage_root, forced_manifest_paths, exact_manifest_paths, source, build)
    output.mkdir(parents=True)
    raw_root = output / "raw"
    raw_root.mkdir()
    spec = unsigned(load_json(stage_root / "path-attestation.json"),
                    ATTESTATION_SCHEMA, "path attestation")
    raw_documents = {}
    raw_artifacts = {}
    for index, case in enumerate(spec["cases"]):
        identifier = case["cell"]["identifier"]
        result_path = raw_root / f"{index:04d}-{identifier}.json"
        arguments = [str(binary)] + [
            str(result_path) if item == "OUTPUT" else item
            for item in case["benchmark_arguments"]
        ]
        try:
            completed = load_exact_main_runner().run_process_bounded(
                arguments, timeout=300.0, max_stdout=1024 * 1024,
                max_stderr=1024 * 1024)
        except Exception as error:
            raise PlanError(
                f"attestation benchmark supervision failed for "
                f"{identifier}: {error}") from error
        require(completed.returncode == 0 and not completed.stdout and
                not completed.stderr,
                f"attestation benchmark failed for {identifier}")
        snapshot = stable_regular_file_bytes(
            result_path, MAX_JSON_DOCUMENT_BYTES,
            f"attestation benchmark result {identifier}")
        document = decode_json_bytes(
            snapshot, f"attestation benchmark result {identifier}")
        validate_attestation_output(
            document, case, spec["expected_resolved_backend"])
        raw_documents[identifier] = document
        raw_artifacts[identifier] = {
            "path": str(result_path.relative_to(output)),
            "size": len(snapshot),
            "sha256": hashlib.sha256(snapshot).hexdigest(),
        }
    final_source, final_build, final_collector = attestation_identities(
        source_root, candidate_commit, binary_relative)
    require((source, build, collector) ==
            (final_source, final_build, final_collector),
            "attestation source or benchmark identity changed during collection")
    result = derive_attestation_result(
        stage_root, source, build, collector, raw_documents, raw_artifacts,
        timing_evidence)
    write_json(output / "manifest.json", result)
    return result


def validate_attestation_result_files(
    stage_root: Path, manifest_path: Path, source: dict[str, Any],
    build: dict[str, Any], collector: dict[str, Any],
    *, verified_timing_evidence: dict[str, Any] | None = None,
) -> dict[str, Any]:
    stage = unsigned(load_json(stage_root / "stage.json"), STAGE_SCHEMA, "stage")
    manifest_path = manifest_path.absolute()
    retained = load_json(manifest_path)
    require(isinstance(retained, dict), "attestation result is not an object")
    result_schema = retained.get("schema")
    require(isinstance(result_schema, str) and result_schema in {
                ATTESTATION_RESULT_SCHEMA_V5,
                ATTESTATION_RESULT_SCHEMA_V6,
                ATTESTATION_RESULT_SCHEMA_V7,
                ATTESTATION_RESULT_SCHEMA_V8,
                ATTESTATION_RESULT_SCHEMA_V9,
                ATTESTATION_RESULT_SCHEMA_V10,
                ATTESTATION_RESULT_SCHEMA_V11,
                ATTESTATION_RESULT_SCHEMA_V12,
                ATTESTATION_RESULT_SCHEMA_V13,
                ATTESTATION_RESULT_SCHEMA},
            "attestation result schema differs")
    unused_canonical_schema, unused_validator, unused_compile_schema, \
        unused_compile_profile, \
        unused_cache, require_timing = \
            _attestation_result_build_contract(result_schema)
    del (unused_canonical_schema, unused_validator, unused_compile_schema,
         unused_compile_profile, unused_cache)
    value = unsigned(retained, result_schema, "attestation result")
    expected_fields = {
        "schema", "status", "valid", "plan_content_sha256",
        "stage_content_sha256", "survivor_content_sha256",
        "attestation_content_sha256", "candidate_commit", "resolved_backend",
        "source_identity", "build_identity", "collector_identity",
        "records",
    }
    if require_timing:
        expected_fields.add("promotion_timing_evidence")
    require(set(value) == expected_fields and
            value["status"] == "complete" and value["valid"] is True and
            value["candidate_commit"] == stage["candidate_commit"],
            "attestation result identity differs")
    require(value["source_identity"] == source and
            value["build_identity"] == build and
            value["collector_identity"] == collector,
            "attestation live source or benchmark identity differs")
    timing_evidence: dict[str, Any] | None = None
    if require_timing:
        if verified_timing_evidence is None:
            validate_promotion_timing_evidence(
                stage_root, value["promotion_timing_evidence"], source, build)
        else:
            validate_promotion_timing_shape(
                stage_root, verified_timing_evidence)
            require(canonical_bytes(value["promotion_timing_evidence"]) ==
                    canonical_bytes(verified_timing_evidence),
                    "attestation timing evidence differs from verified fixture")
        timing_evidence = value["promotion_timing_evidence"]
    else:
        require(verified_timing_evidence is None,
                "non-timing attestation cannot use timing evidence")

    root = manifest_path.parent
    require(root.resolve(strict=True) == root and
            stat.S_ISDIR(os.lstat(root).st_mode),
            "attestation output root traverses a symbolic link")
    records = value["records"]
    require(isinstance(records, list), "attestation records are not a list")
    raw_documents = {}
    raw_artifacts = {}
    expected_files = {Path("manifest.json")}
    for record in records:
        require(isinstance(record, dict) and isinstance(record.get("identifier"), str) and
                isinstance(record.get("raw"), dict),
                "attestation record shape differs")
        identifier = record["identifier"]
        require(identifier not in raw_documents,
                "attestation result identifier is duplicated")
        artifact_value = record["raw"]
        require(set(artifact_value) == {"path", "size", "sha256"} and
                isinstance(artifact_value["path"], str) and
                type(artifact_value["size"]) is int and artifact_value["size"] > 0 and
                isinstance(artifact_value["sha256"], str) and
                re.fullmatch(r"[0-9a-f]{64}",
                             artifact_value["sha256"]) is not None,
                f"attestation artifact shape differs for {identifier}")
        relative = Path(artifact_value["path"])
        require(not relative.is_absolute() and ".." not in relative.parts,
                f"attestation artifact path escapes output for {identifier}")
        path = root / relative
        snapshot = stable_regular_file_bytes(
            path, MAX_JSON_DOCUMENT_BYTES,
            f"attestation artifact {identifier}")
        require(len(snapshot) == artifact_value["size"] and
                hashlib.sha256(snapshot).hexdigest() == artifact_value["sha256"],
                f"attestation artifact bytes changed for {identifier}")
        raw_documents[identifier] = decode_json_bytes(
            snapshot, f"attestation artifact {identifier}")
        raw_artifacts[identifier] = artifact_value
        expected_files.add(relative)
    actual_files = set()
    for path in root.rglob("*"):
        metadata = os.lstat(path)
        require(not stat.S_ISLNK(metadata.st_mode) and
                (stat.S_ISDIR(metadata.st_mode) or
                 stat.S_ISREG(metadata.st_mode)),
                f"attestation output contains an unsafe filesystem object: {path}")
        if stat.S_ISREG(metadata.st_mode):
            actual_files.add(path.relative_to(root))
    require(actual_files == expected_files,
            "attestation output contains a missing or extra file")
    expected = derive_attestation_result(
        stage_root, source, build, collector, raw_documents, raw_artifacts,
        timing_evidence, result_schema=result_schema)
    require(canonical_bytes(retained) == canonical_bytes(expected),
            "attestation result is not the deterministic authenticated result")
    return value


def validate_attestation_result(plan_root: Path, stage_root: Path,
                                manifest_path: Path, source_root: Path,
                                binary_relative: str) -> dict[str, Any]:
    stage = validate_stage(plan_root, stage_root)
    source, build, collector = attestation_identities(
        source_root, stage["candidate_commit"], binary_relative)
    return validate_attestation_result_files(
        stage_root, manifest_path, source, build, collector)


def adversarial_resign(path: Path, mutator) -> None:
    value = load_json(path)
    require(isinstance(value, dict), "mutation target is not an object")
    value.pop("content_sha256", None)
    mutator(value)
    value["content_sha256"] = canonical_sha256(value)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def fake_evidence_scope(
    backend: str = "avx2",
    compile_schema: str = COMPILE_COMMANDS_SCHEMA,
) -> dict[str, Any]:
    require(backend in CAMPAIGN_BACKENDS, "fixture backend is outside campaign")
    require(compile_schema in {
                COMPILE_COMMANDS_SCHEMA_V4, COMPILE_COMMANDS_SCHEMA_V5,
                COMPILE_COMMANDS_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V7,
                COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
                COMPILE_COMMANDS_SCHEMA_V10,
                COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA},
            "fixture compile schema is unsupported")
    current = compile_schema == COMPILE_COMMANDS_SCHEMA
    effective_avx2 = compile_schema in (
        COMPILE_COMMANDS_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V7,
        COMPILE_COMMANDS_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V9,
        COMPILE_COMMANDS_SCHEMA_V10,
        COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA)
    configuration_file_schema = (
        BUILD_CONFIGURATION_FILE_SCHEMA if current else
        BUILD_CONFIGURATION_FILE_SCHEMA
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V11 else
        BUILD_CONFIGURATION_FILE_SCHEMA_V8
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V10 else
        BUILD_CONFIGURATION_FILE_SCHEMA_V7
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V9 else
        BUILD_CONFIGURATION_FILE_SCHEMA_V6
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
        BUILD_CONFIGURATION_FILE_SCHEMA_V5
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V7 else
        BUILD_CONFIGURATION_FILE_SCHEMA_V4
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V6 else
        BUILD_CONFIGURATION_FILE_SCHEMA_V3
        if compile_schema == COMPILE_COMMANDS_SCHEMA_V5 else
        BUILD_CONFIGURATION_FILE_SCHEMA_V2)
    configuration_record_schema, configuration_variables, candidate_cache = \
        _build_configuration_contract_for_compile_schema(
            compile_schema, configuration_file_schema)
    def fixture_artifact(path: str, kind: str, character: str,
                         mode: int = 0o644) -> dict[str, Any]:
        return {"path": path, "kind": kind, "sha256": character * 64,
                "size": 1, "mode": mode}

    def fixture_external_record(
        operand: str, role: str, kind: str,
    ) -> dict[str, Any]:
        metadata, digest, _, resolved, symlink_chain = \
            common.current_external_file_snapshot(
                Path(operand), f"fixture {role}")
        return {
            "operand": operand, "role": role,
            "lexical_symlink_chain": symlink_chain,
            "artifact": {
                "path": str(resolved), "kind": kind,
                "sha256": digest, "size": metadata.st_size,
                "mode": stat.S_IMODE(metadata.st_mode),
            },
        }

    def fixture_text(text: str) -> dict[str, Any]:
        encoded = text.encode("utf-8")
        return {"encoding": "utf-8", "text": text, "size": len(encoded),
                "sha256": hashlib.sha256(encoded).hexdigest()}

    def fixture_commit_object(raw: bytes) -> tuple[str, dict[str, Any]]:
        head = hashlib.sha1(
            f"commit {len(raw)}\0".encode("ascii") + raw).hexdigest()
        return head, {
            "encoding": "base64", "size": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(), "object_id": head,
            "base64": base64.b64encode(raw).decode("ascii"),
        }

    def fixture_build_configuration(
        root: str, source_root: str,
    ) -> dict[str, Any]:
        entries = {
            variable: "" for variable in configuration_variables
        }
        entries.update(_required_build_configuration_entries(
            configuration_file_schema, compile_schema))
        entries.update({
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_GENERATOR": "Unix Makefiles",
            "CMAKE_CONFIGURATION_TYPES": "",
            "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
        })
        material = _build_configuration_material(
            entries, configuration_variables)
        digest = hashlib.sha256(material).hexdigest()
        text = (
            f"schema={configuration_file_schema}\n"
            f"sha256={digest}\n").encode("ascii") + material
        content = fixture_text(text.decode("utf-8"))
        artifact = fixture_artifact(
            f"{root}/{BUILD_CONFIGURATION_RELATIVE_PATH}",
            "generated_build_configuration", "8")
        artifact.update({
            "size": content["size"], "sha256": content["sha256"],
        })
        return {
            "schema": configuration_record_schema,
            "artifact": artifact,
            "content": content,
            "configuration_schema": configuration_file_schema,
            "configuration_sha256": digest,
            "entries": entries,
            "embedded_build_type": "Release",
            "helper_source": fixture_artifact(
                f"{source_root}/{BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH}",
                "source_file", "9"),
        }

    compiler = fixture_artifact(
        "/usr/bin/compiler", "compiler", "a", 0o755)
    compiler_text = "fixture compiler 1.0\n"
    compiler_version = {
        "sha256": hashlib.sha256(compiler_text.encode("utf-8")).hexdigest(),
        "text": compiler_text,
    }

    def fixture_build(role: str, character: str) -> dict[str, Any]:
        baseline = role == "baseline"
        root = f"${role.upper()}_BUILD"
        source_root = "$BASELINE_SOURCE" if baseline else "$CANDIDATE_SOURCE"
        target = "leopard_main_exact.dir" if baseline else "leopard.dir"
        executable_name = "leopard_main_benchmark" if baseline else "bench_leopard2"
        archive_name = "libleopard_main_exact.a" if baseline else "libleopard.a"
        library_names = (
            BASELINE_LIBRARY_SOURCES if baseline else
            CANDIDATE_LIBRARY_SOURCES
            if compile_schema in (
                COMPILE_COMMANDS_SCHEMA_V9, COMPILE_COMMANDS_SCHEMA_V10,
                COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
            CANDIDATE_LIBRARY_SOURCES_V12
            if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
            CANDIDATE_LIBRARY_SOURCES_V11
            if compile_schema in (COMPILE_COMMANDS_SCHEMA_V6,
                                  COMPILE_COMMANDS_SCHEMA_V7) else
            CANDIDATE_LIBRARY_SOURCES_V9)
        library_pairs = []
        for name in library_names:
            source_path = f"{source_root}/{name}"
            object_relative = _normalized_compile_output(role, source_path)
            library_pairs.append({
                "source": fixture_artifact(
                    source_path, "source_file", character),
                "object": fixture_artifact(
                    f"{root}/{object_relative}",
                    "object_file", character),
            })
        benchmark_source = fixture_artifact(
            ("$CANDIDATE_SOURCE/experiments/leopard2/main_compare/"
             "legacy_main_benchmark.cpp" if baseline else
             "$CANDIDATE_SOURCE/bench/leopard2/benchmark.cpp"),
            "source_file", character)
        benchmark_relative = _normalized_compile_output(
            role, benchmark_source["path"])
        benchmark_object = fixture_artifact(
            f"{root}/{benchmark_relative}",
            "object_file", character)
        pairs = sorted([*library_pairs,
            {"source": benchmark_source, "object": benchmark_object},
        ], key=lambda pair: pair["source"]["path"])
        build_configuration = (
            None if baseline else
            fixture_build_configuration(root, source_root))
        cache = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_CONFIGURATION_TYPES": "",
            "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
            "CMAKE_CXX_FLAGS_RELEASE": (
                "-O3 -DNDEBUG" if baseline else "-O3 -DNDEBUG -O3"),
            "CMAKE_GENERATOR": "Unix Makefiles",
            **((REQUIRED_PURE_AVX2_BASELINE_CACHE if effective_avx2 else
                REQUIRED_BASELINE_CACHE) if baseline else candidate_cache),
        }
        if not baseline:
            cache.update({
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    configuration_file_schema,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
                    build_configuration["configuration_sha256"],
            })
        benchmark_relative = benchmark_object["path"][len(root) + 1:]
        objects_by_member = {
            pair["object"]["path"].rsplit("/", 1)[-1]:
                pair["object"]["path"][len(root) + 1:]
            for pair in library_pairs
        }
        members = [Path(name).name + ".o" for name in library_names]
        archive_text = fixture_text(
            f"/usr/bin/ar qc {archive_name} " +
            " ".join(objects_by_member[member] for member in members) + "\n"
            f"/usr/bin/ranlib {archive_name}\n")
        executable_text = fixture_text(
            f"/usr/bin/compiler -O3 {archive_name} {benchmark_relative} "
            f"-o {executable_name} "
            "/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so "
            "/usr/lib/x86_64-linux-gnu/libpthread.a\n")
        archive_recipe = fixture_artifact(
            f"{root}/CMakeFiles/{target}/link.txt", "build_metadata", character)
        archive_recipe.update({
            "size": archive_text["size"], "sha256": archive_text["sha256"]})
        executable_recipe = fixture_artifact(
            f"{root}/CMakeFiles/{executable_name}.dir/link.txt",
            "build_metadata", character)
        executable_recipe.update({
            "size": executable_text["size"],
            "sha256": executable_text["sha256"]})
        return {
            "build_dir": root,
            "validated_cache": cache,
            "compiler": dict(compiler),
            "compiler_invocation": {
                "invocation": "/usr/bin/compiler",
                "resolved_path": "/usr/bin/compiler",
            },
            "compiler_version_stdout": dict(compiler_version),
            "multi_config_build_tool": None,
            "multi_config_build_tool_version_stdout": None,
            "multi_config_ninja_graph": None,
            "cmake_cache": fixture_artifact(
                f"{root}/CMakeCache.txt", "build_metadata", character),
            "compile_commands": fixture_artifact(
                f"{root}/compile_commands.json", "build_metadata", character),
            "executable_link_recipe": executable_recipe,
            "archive_link_recipe": archive_recipe,
            "archive_link_recipe_content": archive_text,
            "executable_link_recipe_content": executable_text,
            "archiver": fixture_artifact(
                "/usr/bin/ar", "archiver", "b", 0o755),
            "ranlib": fixture_artifact(
                "/usr/bin/ranlib", "ranlib", "c", 0o755),
            "archive_link_tool_invocations": {
                "archiver": {"invocation": "/usr/bin/ar",
                             "resolved_path": "/usr/bin/ar"},
                "ranlib": {"invocation": "/usr/bin/ranlib",
                           "resolved_path": "/usr/bin/ranlib"},
            },
            "validated_external_link_inputs": [
                fixture_external_record(
                    "/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so",
                    "openmp_runtime_shared", "shared_library"),
                fixture_external_record(
                    "/usr/lib/x86_64-linux-gnu/libpthread.a",
                    "pthread_support_archive", "archive"),
            ],
            "validated_archive": fixture_artifact(
                f"{root}/{archive_name}", "archive", character),
            "validated_executable": fixture_artifact(
                f"{root}/{executable_name}", "executable", character, 0o755),
            "validated_archive_members": members,
            "validated_compile_commands": {
                "schema": compile_schema,
                "implementation": role,
                "profile": ((BASELINE_PURE_AVX2_COMPILE_PROFILE
                             if effective_avx2 else
                            BASELINE_COMPILE_PROFILE) if baseline else
                            CANDIDATE_COMPILE_PROFILE
                            if compile_schema == COMPILE_COMMANDS_SCHEMA else
                            CANDIDATE_COMPILE_PROFILE_V7
                            if compile_schema == COMPILE_COMMANDS_SCHEMA_V11 else
                            CANDIDATE_COMPILE_PROFILE_V6
                            if compile_schema == COMPILE_COMMANDS_SCHEMA_V10 else
                            CANDIDATE_COMPILE_PROFILE_V5
                            if compile_schema == COMPILE_COMMANDS_SCHEMA_V9 else
                            CANDIDATE_COMPILE_PROFILE_V4
                            if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
                            CANDIDATE_COMPILE_PROFILE_V3
                            if effective_avx2 else
                            CANDIDATE_COMPILE_PROFILE_V2),
                "entry_count": (
                    BASELINE_EXPECTED_COMPILE_COMMAND_COUNT if baseline else
                    CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT
                    if compile_schema in (
                        COMPILE_COMMANDS_SCHEMA_V9,
                        COMPILE_COMMANDS_SCHEMA_V10,
                        COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
                    CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT_V12
                    if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
                    CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT_V11),
                "required_sources": sorted(pair["source"]["path"] for pair in pairs),
                "validated_optimization": "-O3", "validated_openmp": True,
                "required_source_object_pairs": pairs,
                "required_entries": sorted([{
                    "directory": root,
                    "file": pair["source"]["path"],
                    "output": _normalized_compile_output(
                        role, pair["source"]["path"]),
                    "arguments": _normalized_compile_argv(
                        role, pair["source"]["path"], "/usr/bin/compiler",
                        compile_schema,
                        build_configuration=build_configuration),
                } for pair in pairs], key=lambda entry: entry["file"]),
                "isa_policy": (
                    (BASELINE_PURE_AVX2_ISA_POLICY if effective_avx2 else
                     BASELINE_NATIVE_ISA_POLICY) if baseline else
                    CANDIDATE_ISA_POLICY
                    if compile_schema in (
                        COMPILE_COMMANDS_SCHEMA_V9,
                        COMPILE_COMMANDS_SCHEMA_V10,
                        COMPILE_COMMANDS_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA) else
                    CANDIDATE_ISA_POLICY_V12
                    if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
                    CANDIDATE_ISA_POLICY_V11 if effective_avx2 else
                    CANDIDATE_ISA_POLICY_V9),
                "generated_attestation_header": None,
                "effective_build_configuration": build_configuration,
            },
        }

    def fixture_runtime(executable: str, character: str) -> dict[str, Any]:
        raw_text = (
            "linux-vdso.so.1 (0x0000000000000000)\n"
            "libc.so.6 => /lib/libc.so.6 (0x0000000000000000)\n"
            "libm.so.6 => /lib/libm.so.6 (0x0000000000000000)\n"
            "/lib/ld-linux-x86-64.so.2 (0x0000000000000000)\n")
        canonical_text = re.sub(
            r"0x[0-9A-Fa-f]+", CANONICAL_LDD_ADDRESS, raw_text)
        canonical_record = fixture_text(canonical_text)
        return {
            "executable": executable,
            "canonical_ldd_output": {
                "schema": CANONICAL_LDD_SCHEMA,
                "normalization": CANONICAL_LDD_NORMALIZATION,
                **canonical_record,
            },
            "dependencies": [
                {
                    "soname": "ld-linux-x86-64.so.2",
                    "loader_path": "/lib/ld-linux-x86-64.so.2",
                    "file": fixture_artifact(
                        "/lib/ld-linux-x86-64.so.2",
                        "dynamic_loader", character),
                },
                {
                    "soname": "libc.so.6", "loader_path": "/lib/libc.so.6",
                    "file": fixture_artifact(
                        "/lib/libc.so.6",
                        "shared_library", character),
                },
                {
                    "soname": "libm.so.6", "loader_path": "/lib/libm.so.6",
                    "file": fixture_artifact(
                        "/lib/libm.so.6", "shared_library",
                        "a" if character != "a" else "b"),
                },
                {"soname": "linux-vdso.so.1", "virtual": True},
            ],
        }

    def fixture_cpu(cpu: int, sibling: int) -> dict[str, Any]:
        return {
            "cpu": cpu, "online": "1",
            "cpuinfo": {"processor": str(cpu), "model name": "fixture"},
            "topology": {
                "core_id": "2", "physical_package_id": "0", "die_id": "0",
                "cluster_id": None, "thread_siblings_list": "2,6",
                "core_siblings_list": "0-7",
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
                    "coherency_line_size": "64", "number_of_sets": "1",
                    "ways_of_associativity": "1",
                    "physical_line_partition": "1",
                    "shared_cpu_list": "2,6", "shared_cpu_map": "44",
                    "allocation_policy": None, "write_policy": "WriteBack",
                },
                {
                    "index": 1, "level": "3", "type": "Unified", "size": "8M",
                    "coherency_line_size": "64", "number_of_sets": "1",
                    "ways_of_associativity": "1",
                    "physical_line_partition": "1",
                    "shared_cpu_list": "0-7", "shared_cpu_map": "ff",
                    "allocation_policy": None, "write_policy": "WriteBack",
                },
            ],
            "cache_index_inventory": ["index0", "index1"],
            "cache_directory_inventory_text": fixture_text("index0\nindex1\n"),
            "numa_nodes": [0],
            "numa_node_inventory": ["node0"],
            "core_class": {"core_type": None, "cpu_capacity": None},
        }

    baseline_build = fixture_build("baseline", "d")
    candidate_build = fixture_build("candidate", "e")
    baseline_archive = baseline_build["validated_archive"]
    baseline_executable = baseline_build["validated_executable"]
    candidate_archive = candidate_build["validated_archive"]
    candidate_executable = candidate_build["validated_executable"]
    baseline_raw = base64.b64decode(EXACT_MAIN_COMMIT_BASE64, validate=True)
    baseline_head, baseline_commit_object = fixture_commit_object(baseline_raw)
    require(baseline_head == EXACT_MAIN_COMMIT,
            "embedded exact-main commit fixture differs")
    candidate_tree = "c" * 40
    candidate_raw = (
        f"tree {candidate_tree}\nauthor Fixture <fixture@example.com> 1 +0000\n"
        "committer Fixture <fixture@example.com> 1 +0000\n\nfixture\n"
    ).encode("utf-8")
    candidate_head, candidate_commit_object = fixture_commit_object(candidate_raw)
    attestation_text = _benchmark_attestation_text(
        candidate_head, candidate_tree)
    attestation_content = fixture_text(attestation_text)
    attestation_artifact = fixture_artifact(
        "$CANDIDATE_BUILD/generated/leopard2-benchmark-attestation/"
        "leopard2_benchmark_source_attestation.h",
        "generated_compile_input", "f")
    attestation_artifact.update({
        "size": attestation_content["size"],
        "sha256": attestation_content["sha256"],
    })
    candidate_build["validated_compile_commands"][
        "generated_attestation_header"] = {
            "artifact": attestation_artifact,
            "content": attestation_content,
            "source_commit": candidate_head,
            "source_tree": candidate_tree,
            "source_tracked_dirty": False,
        }
    return {
        "schema": (
            EVIDENCE_SCOPE_SCHEMA if current else
            EVIDENCE_SCOPE_SCHEMA_V11
            if compile_schema == COMPILE_COMMANDS_SCHEMA_V11 else
            EVIDENCE_SCOPE_SCHEMA_V10
            if compile_schema == COMPILE_COMMANDS_SCHEMA_V10 else
            EVIDENCE_SCOPE_SCHEMA_V9
            if compile_schema == COMPILE_COMMANDS_SCHEMA_V9 else
            EVIDENCE_SCOPE_SCHEMA_V8
            if compile_schema == COMPILE_COMMANDS_SCHEMA_V8 else
            EVIDENCE_SCOPE_SCHEMA_V7
            if compile_schema == COMPILE_COMMANDS_SCHEMA_V7 else
            EVIDENCE_SCOPE_SCHEMA_V6
            if compile_schema == COMPILE_COMMANDS_SCHEMA_V6 else
            EVIDENCE_SCOPE_SCHEMA_V5
            if compile_schema == COMPILE_COMMANDS_SCHEMA_V5 else
            EVIDENCE_SCOPE_SCHEMA_V4),
        "host": {
            "system": {
                "system": "Linux", "node": "fixture", "release": "fixture",
                "version": "fixture", "machine": "x86_64", "python": "3.11",
                "page_size": 4096,
            },
            "allowed_cpu_set_at_launch": list(range(8)),
            "online_cpu_set": list(range(8)),
            "online_cpu_list_text": fixture_text("0-7\n"),
            "online_node_list_text": fixture_text("0\n"),
            "benchmark_cpu": fixture_cpu(2, 6),
            "reserved_sibling": fixture_cpu(6, 2),
            "turbo_and_pstate": {
                "/sys/devices/system/cpu/intel_pstate/no_turbo": "0",
                "/sys/devices/system/cpu/amd_pstate/status": None,
                "/sys/devices/system/cpu/cpufreq/boost": None,
            },
        },
        "sources": {
            "baseline": {"path": "$BASELINE_SOURCE", "head": EXACT_MAIN_COMMIT,
                         "tree": EXACT_MAIN_TREE, "detached": True,
                         "tracked_tree_listing_sha256": "b" * 64,
                         "tracked_status": "clean",
                         "commit_object": baseline_commit_object},
            "candidate": {"path": "$CANDIDATE_SOURCE", "head": candidate_head,
                          "tree": candidate_tree, "detached": False,
                          "tracked_tree_listing_sha256": "d" * 64,
                          "tracked_status": "clean",
                          "commit_object": candidate_commit_object},
        },
        "builds": {
            "baseline": baseline_build,
            "candidate": candidate_build,
        },
        "artifacts": {
            "baseline_archive": baseline_archive,
            "baseline_executable": baseline_executable,
            "candidate_archive": candidate_archive,
            "candidate_executable": candidate_executable,
        },
        "runtime_closures": {
            "baseline": fixture_runtime(baseline_executable["path"], "1"),
            "candidate": fixture_runtime(candidate_executable["path"], "2"),
        },
        "tools": {
            "runner": fixture_artifact("/fixture/run.py", "file", "3"),
            "taskset": fixture_artifact(
                "/usr/bin/taskset", "executable", "4", 0o755),
            "ldd": fixture_artifact(
                "/usr/bin/ldd", "executable", "5", 0o755),
            "evidence_helper": fixture_artifact(
                f"$CANDIDATE_SOURCE/{EVIDENCE_HELPER_RELATIVE_PATH}",
                "file", "6"),
        },
        "resolved_auto_backend": backend,
        "forced_confirmation_backends": list(
            BACKENDS[:BACKENDS.index(backend) + 1]),
        "excluded_backends": dict(EXCLUDED_CAMPAIGN_BACKENDS),
    }


def fake_gate_manifests(
    plan_root: Path, backend: str = "avx2"
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    plan = validate_plan(plan_root)
    manifests = []
    references = []
    scopes = []
    for item in plan["artifacts"]:
        scope = fake_evidence_scope(backend)
        gate = validate_main_document(load_json(plan_root / item["path"]), item["path"])
        campaign_cells = []
        analysis = {}
        for cell in gate["cells"]:
            campaign_cells.append({
                "identifier": cell["identifier"], "k": cell["K"], "r": cell["R"],
                "shard_bytes": cell["shard_bytes"], "losses": cell["loss_count"],
                "seed": cell["seed"],
            })
            passes = (cell["K"], cell["shard_bytes"]) in {
                (127, 65536), (128, 4096),
            }
            lower = 1.06 if passes else 0.99
            analysis[cell["identifier"]] = {
                metric: {"ci95_lower": lower,
                         "ratio_orientation": "baseline_time_over_candidate_time"}
                for metric in ("decode_first_use", "decode_reuse_amortized")
            }
        manifests.append({
            "schema": EXACT_MANIFEST_SCHEMA, "valid": True,
            "campaign": {
                "candidate_mode": "generic", "batch": 1, "reuse": 8,
                "iterations": 9, "warmup": 2, "threads": 1, "rounds": 3,
                "order": list(EXACT_ORDER), "cells": campaign_cells,
            },
            "identities": {
                "baseline_source": {"head": EXACT_MAIN_COMMIT},
                "candidate_source": {
                    "head": scope["sources"]["candidate"]["head"]},
            },
            "analysis": analysis,
        })
        references.append({
            "path": "/fixture/" + Path(item["path"]).name,
            "size": 1, "sha256": "2" * 64, "payload_digest": "3" * 64,
        })
        scopes.append(scope)
    return manifests, references, scopes


def fake_attestation_output(case: dict[str, Any], selected_path: str | None = None,
                            selected_rule: str | None = None,
                            backend: str = "avx2") -> dict[str, Any]:
    cell = case["cell"]
    expected_path, expected_rule = expected_auto_decode_pair(case, backend)
    if selected_path is None:
        selected_path = expected_path
    if selected_rule is None:
        selected_rule = expected_rule
    transform = ceil_power_of_two(cell["R"])
    parent = ceil_power_of_two(cell["K"] + transform)
    required_slots = expected_required_work_slots(case, selected_path)
    return {
        "schema": common.BENCHMARK_SCHEMA,
        "build": {
            "compiler": "fixture", "compiler_version": "1", "cplusplus": 201103,
        },
        "parameters": {
            "K": cell["K"], "R": cell["R"],
            "requested_profile": "legacy_high_v1", "requested_field": "gf8",
            "requested_backend": "auto", "force_generic_decode": False,
            "force_specialized_decode": False, "force_tiled_decode": False,
            "force_materialized_decode": False, "skip_legacy": True,
            "retain_samples": True, "report_decode_path": True,
            "shard_bytes": cell["shard_bytes"],
            "loss_count": cell["loss_count"],
            "missing_original_indices": expected_missing_indices(
                cell["K"], cell["loss_count"], cell["seed"]),
            "batch": 1, "reuse": 1, "iterations": 1, "warmup": 0,
            "thread_count": 1, "seed": cell["seed"],
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8", "backend": backend,
            "thread_count": 1, "parent_count": parent,
            "padded_side": transform, "selected_decode_path": selected_path,
            "selected_decode_rule": selected_rule,
            "decode_required_work_slots": required_slots,
            "decode_aligned_prefix_bytes": cell["shard_bytes"] & ~63,
            "decode_tail_bytes": cell["shard_bytes"] & 63,
            "decode_rounded_bytes": (cell["shard_bytes"] + 63) & ~63,
            "decode_multi_item_batch": False,
        },
        "correctness": {
            "leopard2_round_trip": True, "legacy_comparison": None,
        },
        "workload_digests": {
            "algorithm": "fnv1a64", "original_data": "1" * 16,
            "transmitted_parity": "2" * 16, "recovered_originals": "3" * 16,
        },
        "memory": {
            "scratch_alignment": 64, "encode_scratch_bytes_per_stripe": 64,
            "decode_scratch_bytes_per_stripe": 128,
            "encode_scratch_bytes_batch": 64, "decode_scratch_bytes_batch": 128,
        },
        "metrics": {
            "codec_setup": {"samples_us": [1.0]},
            "encode_execution": {"samples_us_per_batch_call": [1.0]},
            "decode_plan_setup": {"samples_us": [1.0]},
            "decode_execution": {"samples_us_per_batch_call": [1.0]},
            "decode_amortized_at_reuse": {"reuse_count": 1},
            "rate_semantics": "fixture",
        },
        "legacy": {
            "available": False, "unavailable_reason": "disabled by --skip-legacy",
            "codec_setup": None, "decode_timing_includes_setup": True,
            "encode_execution": None, "decode_including_setup": None,
        },
    }


def self_test() -> None:
    def retained_text(text: str) -> dict[str, Any]:
        encoded = text.encode("utf-8")
        return {
            "encoding": "utf-8", "text": text, "size": len(encoded),
            "sha256": hashlib.sha256(encoded).hexdigest(),
        }

    def runtime_text(closure: dict[str, Any]) -> str:
        key = ("canonical_ldd_output" if "canonical_ldd_output" in closure else
               "raw_ldd_output")
        return closure[key]["text"]

    def replace_runtime_text(closure: dict[str, Any], text: str) -> None:
        if "canonical_ldd_output" in closure:
            closure["canonical_ldd_output"] = {
                "schema": CANONICAL_LDD_SCHEMA,
                "normalization": CANONICAL_LDD_NORMALIZATION,
                **retained_text(text),
            }
        else:
            closure["raw_ldd_output"] = retained_text(text)

    exact_runner = load_exact_main_runner()
    require(tuple(exact_runner.BASELINE_LIBRARY_SOURCES) ==
                BASELINE_LIBRARY_SOURCES and
            tuple(exact_runner.CANDIDATE_LIBRARY_SOURCES_V9) ==
                CANDIDATE_LIBRARY_SOURCES_V9 and
            tuple(exact_runner.CANDIDATE_LIBRARY_SOURCES_V11) ==
                CANDIDATE_LIBRARY_SOURCES_V11 and
            tuple(exact_runner.CANDIDATE_LIBRARY_SOURCES_V12) ==
                CANDIDATE_LIBRARY_SOURCES_V12 and
            tuple(exact_runner.CANDIDATE_LIBRARY_SOURCES_V16) ==
                CANDIDATE_LIBRARY_SOURCES and
            exact_runner.BASELINE_EXPECTED_COMPILE_COMMAND_COUNT ==
                BASELINE_EXPECTED_COMPILE_COMMAND_COUNT and
            exact_runner.CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT_V16 ==
                CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT and
            len(exact_runner.CANDIDATE_LIBRARY_SOURCES_V11) == 13 and
            len(exact_runner.CANDIDATE_LIBRARY_SOURCES_V11) +
                len(exact_runner.CANDIDATE_NON_LIBRARY_COMPILE_ACTIONS) ==
                CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT_V11 and
            len(exact_runner.CANDIDATE_LIBRARY_SOURCES_V12) +
                len(exact_runner.CANDIDATE_NON_LIBRARY_COMPILE_ACTIONS) ==
                CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT_V12 and
            exact_runner.COMPILE_COMMANDS_SCHEMA_V5 ==
                COMPILE_COMMANDS_SCHEMA_V5 and
            exact_runner.COMPILE_COMMANDS_SCHEMA_V6 ==
                COMPILE_COMMANDS_SCHEMA_V6 and
            exact_runner.COMPILE_COMMANDS_SCHEMA_V7 ==
                COMPILE_COMMANDS_SCHEMA_V7 and
            exact_runner.COMPILE_COMMANDS_SCHEMA_V8 ==
                COMPILE_COMMANDS_SCHEMA_V8 and
            exact_runner.COMPILE_COMMANDS_SCHEMA_V9 ==
                COMPILE_COMMANDS_SCHEMA_V9 and
            exact_runner.COMPILE_COMMANDS_SCHEMA_V10 ==
                COMPILE_COMMANDS_SCHEMA_V10 and
            exact_runner.COMPILE_COMMANDS_SCHEMA_V11 ==
                COMPILE_COMMANDS_SCHEMA_V11 and
            exact_runner.COMPILE_COMMANDS_SCHEMA_V12 ==
                COMPILE_COMMANDS_SCHEMA and
            exact_runner.BASELINE_PURE_AVX2_COMPILE_PROFILE ==
                BASELINE_PURE_AVX2_COMPILE_PROFILE and
            exact_runner.CANDIDATE_COMPILE_PROFILE_V2 ==
                CANDIDATE_COMPILE_PROFILE_V2 and
            exact_runner.CANDIDATE_COMPILE_PROFILE_V3 ==
                CANDIDATE_COMPILE_PROFILE_V3 and
            exact_runner.CANDIDATE_COMPILE_PROFILE_V4 ==
                CANDIDATE_COMPILE_PROFILE_V4 and
            exact_runner.CANDIDATE_COMPILE_PROFILE_V5 ==
                CANDIDATE_COMPILE_PROFILE_V5 and
            exact_runner.CANDIDATE_COMPILE_PROFILE_V6 ==
                CANDIDATE_COMPILE_PROFILE_V6 and
            exact_runner.CANDIDATE_COMPILE_PROFILE_V7 ==
                CANDIDATE_COMPILE_PROFILE_V7 and
            exact_runner.CANDIDATE_COMPILE_PROFILE_V8 ==
                CANDIDATE_COMPILE_PROFILE and
            exact_runner.candidate_isa_policy(exact_runner.RAW_SCHEMA_V9) ==
                CANDIDATE_ISA_POLICY_V9 and
            exact_runner.candidate_isa_policy(exact_runner.RAW_SCHEMA_V11) ==
                CANDIDATE_ISA_POLICY_V11 and
            exact_runner.candidate_isa_policy(exact_runner.RAW_SCHEMA_V12) ==
                CANDIDATE_ISA_POLICY_V12 and
            exact_runner.candidate_isa_policy(exact_runner.RAW_SCHEMA_V13) ==
                CANDIDATE_ISA_POLICY and
            exact_runner.candidate_isa_policy(exact_runner.RAW_SCHEMA_V14) ==
                CANDIDATE_ISA_POLICY and
            exact_runner.candidate_isa_policy(exact_runner.RAW_SCHEMA_V15) ==
                CANDIDATE_ISA_POLICY and
            exact_runner.candidate_isa_policy(exact_runner.RAW_SCHEMA_V16) ==
                CANDIDATE_ISA_POLICY and
            tuple(exact_runner.BUILD_CONFIGURATION_VARIABLES_V9) ==
                BUILD_CONFIGURATION_VARIABLES and
            tuple(exact_runner.BUILD_CONFIGURATION_VARIABLES_V6) ==
                BUILD_CONFIGURATION_VARIABLES_V6 and
            tuple(exact_runner.BUILD_CONFIGURATION_VARIABLES_V7) ==
                BUILD_CONFIGURATION_VARIABLES_V7 and
            tuple(exact_runner.BUILD_CONFIGURATION_VARIABLES_V8) ==
                BUILD_CONFIGURATION_VARIABLES_V8 and
            tuple(exact_runner.BUILD_CONFIGURATION_VARIABLES_V5) ==
                BUILD_CONFIGURATION_VARIABLES_V5 and
            tuple(exact_runner.BUILD_CONFIGURATION_VARIABLES_V4) ==
                BUILD_CONFIGURATION_VARIABLES_V4 and
            tuple(exact_runner.BUILD_CONFIGURATION_VARIABLES_V3) ==
                BUILD_CONFIGURATION_VARIABLES_V3 and
            exact_runner.BUILD_CONFIGURATION_FILE_SCHEMA_V9 ==
                BUILD_CONFIGURATION_FILE_SCHEMA and
            exact_runner.BUILD_CONFIGURATION_FILE_SCHEMA_V6 ==
                BUILD_CONFIGURATION_FILE_SCHEMA_V6 and
            exact_runner.BUILD_CONFIGURATION_FILE_SCHEMA_V7 ==
                BUILD_CONFIGURATION_FILE_SCHEMA_V7 and
            exact_runner.BUILD_CONFIGURATION_FILE_SCHEMA_V8 ==
                BUILD_CONFIGURATION_FILE_SCHEMA_V8 and
            exact_runner.BUILD_CONFIGURATION_FILE_SCHEMA_V5 ==
                BUILD_CONFIGURATION_FILE_SCHEMA_V5 and
            exact_runner.BUILD_CONFIGURATION_FILE_SCHEMA_V4 ==
                BUILD_CONFIGURATION_FILE_SCHEMA_V4 and
            exact_runner.BUILD_CONFIGURATION_FILE_SCHEMA_V3 ==
                BUILD_CONFIGURATION_FILE_SCHEMA_V3 and
            exact_runner.BUILD_CONFIGURATION_RECORD_SCHEMA_V9 ==
                BUILD_CONFIGURATION_RECORD_SCHEMA and
            exact_runner.BUILD_CONFIGURATION_RECORD_SCHEMA_V6 ==
                BUILD_CONFIGURATION_RECORD_SCHEMA_V6 and
            exact_runner.BUILD_CONFIGURATION_RECORD_SCHEMA_V7 ==
                BUILD_CONFIGURATION_RECORD_SCHEMA_V7 and
            exact_runner.BUILD_CONFIGURATION_RECORD_SCHEMA_V8 ==
                BUILD_CONFIGURATION_RECORD_SCHEMA_V8 and
            exact_runner.BUILD_CONFIGURATION_RECORD_SCHEMA_V5 ==
                BUILD_CONFIGURATION_RECORD_SCHEMA_V5 and
            exact_runner.BUILD_CONFIGURATION_RECORD_SCHEMA_V4 ==
                BUILD_CONFIGURATION_RECORD_SCHEMA_V4 and
            exact_runner.BUILD_CONFIGURATION_RECORD_SCHEMA_V3 ==
                BUILD_CONFIGURATION_RECORD_SCHEMA_V3 and
            exact_runner.BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH ==
                BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH and
            exact_runner.EVIDENCE_HELPER_RELATIVE_PATH ==
                EVIDENCE_HELPER_RELATIVE_PATH,
            "balanced scope translation-unit contract drifted from its producer")
    require(exact_runner.MANIFEST_SCHEMA_V16 == EXACT_MANIFEST_SCHEMA and
            exact_runner.MANIFEST_SCHEMA_V15 == EXACT_MANIFEST_SCHEMA_V15 and
            exact_runner.MANIFEST_SCHEMA_V14 == EXACT_MANIFEST_SCHEMA_V14 and
            exact_runner.MANIFEST_SCHEMA_V13 == EXACT_MANIFEST_SCHEMA_V13 and
            exact_runner.MANIFEST_SCHEMA_V12 == EXACT_MANIFEST_SCHEMA_V12 and
            exact_runner.MANIFEST_SCHEMA_V11 == EXACT_MANIFEST_SCHEMA_V11 and
            exact_runner.MANIFEST_SCHEMA_V10 == EXACT_MANIFEST_SCHEMA_V10 and
            exact_runner.RAW_SCHEMA_V15 == EXACT_RAW_SCHEMA_V15 and
            exact_runner.RAW_SCHEMA_V14 == EXACT_RAW_SCHEMA_V14 and
            exact_runner.RAW_SCHEMA_V13 == EXACT_RAW_SCHEMA_V13 and
            exact_runner.RAW_SCHEMA_V12 == EXACT_RAW_SCHEMA_V12 and
            exact_runner.RAW_SCHEMA_V11 == EXACT_RAW_SCHEMA_V11 and
            exact_runner.RAW_SCHEMA_V10 == EXACT_RAW_SCHEMA_V10 and
            exact_runner.RAW_SCHEMA_V16 == EXACT_RAW_SCHEMA and
            exact_runner.MANIFEST_SCHEMA != EXACT_MANIFEST_SCHEMA and
            exact_runner.RAW_SCHEMA != EXACT_RAW_SCHEMA,
            "balanced exact-main v16 schema freeze drifted from its producer")
    _validate_exact_schema_pair(
        {"schema": EXACT_MANIFEST_SCHEMA},
        {"schema": EXACT_RAW_SCHEMA})
    for manifest_schema, raw_schema in (
        (exact_runner.MANIFEST_SCHEMA, exact_runner.RAW_SCHEMA),
        (exact_runner.MANIFEST_SCHEMA, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, exact_runner.RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA_V15, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA_V15),
        (EXACT_MANIFEST_SCHEMA_V14, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA_V14),
        (EXACT_MANIFEST_SCHEMA_V13, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA_V13),
        (EXACT_MANIFEST_SCHEMA_V12, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA_V12),
        (EXACT_MANIFEST_SCHEMA_V11, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA_V11),
        (EXACT_MANIFEST_SCHEMA_V10, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA_V10),
        (EXACT_MANIFEST_SCHEMA_V9, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA_V9),
        (EXACT_MANIFEST_SCHEMA_V8, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA_V8),
        (EXACT_MANIFEST_SCHEMA_V7, EXACT_RAW_SCHEMA),
        (EXACT_MANIFEST_SCHEMA, EXACT_RAW_SCHEMA_V7),
        (EXACT_MANIFEST_SCHEMA_V6, EXACT_RAW_SCHEMA_V5),
    ):
        try:
            _validate_exact_schema_pair(
                {"schema": manifest_schema}, {"schema": raw_schema})
        except PlanError:
            pass
        else:
            raise PlanError(
                "cross-version exact-main schema downgrade was accepted")
    with tempfile.TemporaryDirectory(prefix="leopard2-balanced-plan-") as temporary:
        root = Path(temporary)
        first = root / "first"
        second = root / "second"
        plan_a = generate_plan(first)
        validate_plan(first)
        plan_b = generate_plan(second)
        validate_plan(second)
        require(canonical_bytes(plan_a) == canonical_bytes(plan_b),
                "two canonical plans differ")

        exact_type_gate = load_json(
            first / "exact-main-gate" / "t128.json")
        require(isinstance(exact_type_gate, dict),
                "self-test exact-main gate is not an object")
        exact_type_gate["reuse"] = 8.0
        exact_type_gate["content_sha256"] = canonical_sha256({
            key: value for key, value in exact_type_gate.items()
            if key != "content_sha256"
        })
        try:
            validate_main_document(exact_type_gate, "type-coerced gate")
        except PlanError:
            pass
        else:
            raise PlanError(
                "balanced plan accepted float substitution for integer reuse")

        for name, payload in (
            ("duplicate.json", b'{"value":1,"value":2}'),
            ("overflow.json", b'{"value":1e9999}'),
            ("constant.json", b'{"value":NaN}'),
            ("huge-integer.json", (
                b'{"value":' + b"9" * 5000 + b"}")),
        ):
            path = root / name
            path.write_bytes(payload)
            try:
                load_json(path)
            except PlanError:
                continue
            raise PlanError(
                f"balanced plan accepted ambiguous JSON fixture {name}")

        ordinary_json = root / "ordinary.json"
        ordinary_json.write_bytes(b'{"value":1}\n')
        symbolic_json = root / "symbolic.json"
        symbolic_json.symlink_to(ordinary_json)
        try:
            load_json(symbolic_json)
        except PlanError:
            pass
        else:
            raise PlanError("balanced plan followed a JSON symbolic link")
        hard_link_json = root / "hard-link.json"
        os.link(ordinary_json, hard_link_json)
        try:
            load_json(hard_link_json)
        except PlanError:
            pass
        else:
            raise PlanError("balanced plan accepted a hard-linked JSON file")
        real_parent = root / "real-json-parent"
        real_parent.mkdir()
        nested_json = real_parent / "nested.json"
        nested_json.write_bytes(b'{"value":1}\n')
        parent_link = root / "json-parent-link"
        parent_link.symlink_to(real_parent, target_is_directory=True)
        try:
            load_json(parent_link / "nested.json")
        except PlanError:
            pass
        else:
            raise PlanError(
                "balanced plan followed a symbolic-link JSON parent")

        secure_output = root / "secure-output" / "nested" / "result.json"
        secure_document = {"schema": "secure-write-fixture/v1", "value": 7}
        write_json(secure_output, secure_document)
        require(load_json(secure_output) == secure_document,
                "secure JSON publication changed its payload")
        require(stat.S_IMODE(os.lstat(secure_output).st_mode) == 0o600 and
                stat.S_IMODE(os.lstat(secure_output.parent).st_mode) == 0o700,
                "secure JSON publication used permissive new-file modes")
        original_secure_bytes = secure_output.read_bytes()
        try:
            write_json(secure_output, {"value": 8})
        except PlanError:
            pass
        else:
            raise PlanError("secure JSON publication replaced existing output")
        require(secure_output.read_bytes() == original_secure_bytes,
                "failed duplicate JSON publication changed existing output")

        symlink_victim = root / "json-symlink-victim"
        symlink_victim.write_bytes(b"victim\n")
        symlink_output = root / "json-symlink-output"
        symlink_output.symlink_to(symlink_victim)
        try:
            write_json(symlink_output, {"value": 9})
        except PlanError:
            pass
        else:
            raise PlanError("JSON publication replaced a final symlink")
        require(symlink_victim.read_bytes() == b"victim\n",
                "JSON publication changed a symlink victim")

        linked_output = root / "linked-json-parent" / "result.json"
        linked_parent = linked_output.parent
        linked_parent.symlink_to(real_parent, target_is_directory=True)
        try:
            write_json(linked_output, {"value": 10})
        except PlanError:
            pass
        else:
            raise PlanError("JSON publication followed a parent symlink")
        require(not (real_parent / "result.json").exists(),
                "parent-symlink JSON publication changed its target")

        trusted_parent = root / "trusted-json-parent"
        trusted_parent.mkdir()
        saved_parent = root / "saved-json-parent"
        hostile_parent = root / "hostile-json-parent"
        hostile_parent.mkdir()
        swapped_output = trusted_parent / "result.json"
        original_link = os.link
        swap_performed = False

        def swap_parent_before_link(
                source, destination, *args, **kwargs):
            nonlocal swap_performed
            if not swap_performed:
                trusted_parent.rename(saved_parent)
                trusted_parent.symlink_to(
                    hostile_parent, target_is_directory=True)
                swap_performed = True
            return original_link(source, destination, *args, **kwargs)

        os.link = swap_parent_before_link
        try:
            try:
                write_json(swapped_output, {"value": 11})
            except PlanError:
                pass
            else:
                raise PlanError(
                    "JSON publication accepted a parent-directory swap")
        finally:
            os.link = original_link
        require(swap_performed and
                not (hostile_parent / "result.json").exists() and
                not (saved_parent / "result.json").exists() and
                not tuple(saved_parent.glob(".result.json.*.tmp")),
                "parent-swap JSON publication changed a victim or leaked state")

        forged = root / "forged"
        shutil.copytree(first, forged)
        gate_path = forged / "exact-main-gate" / "t128.json"
        adversarial_resign(gate_path, lambda value: value.update(
            {"candidate_mode": "auto"}))
        plan_path = forged / "plan.json"
        def reseal_plan(value):
            target = next(item for item in value["artifacts"]
                          if item["path"] == "exact-main-gate/t128.json")
            target["size"] = gate_path.stat().st_size
            target["sha256"] = file_sha256(gate_path)
        adversarial_resign(plan_path, reseal_plan)
        try:
            validate_plan(forged)
        except PlanError:
            pass
        else:
            raise PlanError("re-signed mode forgery passed canonical validation")

        forged_gate = root / "forged-gate"
        shutil.copytree(first, forged_gate)
        adversarial_resign(forged_gate / "plan.json", lambda value:
                           value["promotion_gate"].update(
                               {"minimum_ci95_lower": 0.01}))
        try:
            validate_plan(forged_gate)
        except PlanError:
            pass
        else:
            raise PlanError("re-signed gate forgery passed canonical validation")

        manifests, references, scopes = fake_gate_manifests(first)
        survivor_signed = derive_survivors(first, manifests, references, scopes)
        require(not survivor_signed["required_refinement_cells"] and
                [(cell["K"], cell["shard_bytes"])
                 for cell in survivor_signed["survivor_cells"]] == [
                    (128, 4096), (127, 65536),
                 ],
                "fixture survivor selection differs")
        require(forced_backends_for_scope(fake_evidence_scope("scalar")) ==
                    ("scalar",) and
                forced_backends_for_scope(fake_evidence_scope("ssse3")) ==
                    ("scalar", "ssse3") and
                forced_backends_for_scope(fake_evidence_scope("avx2")) == BACKENDS,
                "forced backend prefix does not match resolved AUTO tier")
        reordered_scope = json.loads(json.dumps(fake_evidence_scope()))
        baseline_semantics = reordered_scope["builds"]["baseline"][
            "validated_compile_commands"]
        baseline_pairs = baseline_semantics["required_source_object_pairs"]
        benchmark_pair = next(pair for pair in baseline_pairs
                              if pair["source"]["path"].endswith(
                                  "/legacy_main_benchmark.cpp"))
        baseline_pairs.remove(benchmark_pair)
        baseline_pairs.insert(0, benchmark_pair)
        baseline_semantics["required_sources"] = [
            pair["source"]["path"] for pair in baseline_pairs]
        validate_evidence_scope(reordered_scope)
        embedded_entries = reordered_scope["builds"]["candidate"][
            "validated_compile_commands"]["effective_build_configuration"][
                "entries"]
        _validate_embedded_build_type(embedded_entries, "Release")
        multi_entries = dict(embedded_entries)
        multi_entries.update({
            "CMAKE_BUILD_TYPE": "",
            "CMAKE_GENERATOR": "Ninja Multi-Config",
            "CMAKE_CONFIGURATION_TYPES": "Debug;Release;RelWithDebInfo",
        })
        _validate_embedded_build_type(multi_entries, "Release")
        for encoded_types, embedded in (
            ("Debug;RelWithDebInfo", "Release"),
            ("Debug;Release;Release", "Release"),
            ("Debug;Release", "Debug"),
        ):
            invalid = dict(multi_entries)
            invalid["CMAKE_CONFIGURATION_TYPES"] = encoded_types
            try:
                _validate_embedded_build_type(invalid, embedded)
            except PlanError:
                pass
            else:
                raise PlanError(
                    "invalid normalized multi-config build type was accepted")

        def reject_scope_mutation(label: str, mutate) -> None:
            forged_scopes = json.loads(json.dumps(scopes))
            mutate(forged_scopes)
            try:
                derive_survivors(first, manifests, references, forged_scopes)
            except PlanError:
                return
            raise PlanError(f"mixed/invalid evidence scope accepted {label}")

        reject_scope_mutation("host", lambda values:
                              values[-1]["host"]["system"].update({
                                  "machine": "aarch64"}))
        def move_to_other_cache_class(values) -> None:
            pair = values[-1]["host"]
            pair["benchmark_cpu"].update({"cpu": 3})
            pair["benchmark_cpu"]["cpuinfo"].update({
                "processor": "3", "cache domain": "32MiB"})
            pair["benchmark_cpu"]["topology"].update({
                "core_id": "3", "thread_siblings_list": "3,7"})
            pair["reserved_sibling"].update({"cpu": 7})
            pair["reserved_sibling"]["cpuinfo"].update({
                "processor": "7", "cache domain": "32MiB"})
            pair["reserved_sibling"]["topology"].update({
                "core_id": "3", "thread_siblings_list": "3,7"})
        reject_scope_mutation("equal-cardinality mixed cache class",
                              move_to_other_cache_class)
        reject_scope_mutation("candidate compiler", lambda values:
                              values[-1]["builds"]["candidate"]["compiler"].update({
                                  "sha256": "d" * 64}))
        reject_scope_mutation("baseline source", lambda values:
                              values[-1]["sources"]["baseline"].update({
                                  "tree": "d" * 40}))
        reject_scope_mutation("candidate CMake", lambda values:
                              values[-1]["builds"]["candidate"][
                                  "validated_cache"].update({
                                      "LEO2_BACKEND_VARIANT": "scalar"}))
        reject_scope_mutation("uniform noncanonical candidate CMake", lambda values:
                              [value["builds"]["candidate"][
                                  "validated_cache"].update({
                                  "LEO2_BACKEND_VARIANT": "scalar"})
                               for value in values])
        reject_scope_mutation(
            "candidate evidence helper path",
            lambda values: values[-1]["tools"]["evidence_helper"].update({
                "path": "$CANDIDATE_SOURCE/experiments/leopard2/"
                        "main_compare/run_abba.py"}))
        def downgrade_scope_and_drop_helper(values) -> None:
            values[-1]["schema"] = EVIDENCE_SCOPE_SCHEMA_V3
            values[-1]["tools"].pop("evidence_helper")
        reject_scope_mutation(
            "scope downgrade drops evidence helper",
            downgrade_scope_and_drop_helper)
        reject_scope_mutation(
            "candidate configuration helper path",
            lambda values: values[-1]["builds"]["candidate"][
                "validated_compile_commands"][
                    "effective_build_configuration"]["helper_source"].update({
                        "path": "$CANDIDATE_SOURCE/cmake/foreign.cmake"}))
        reject_scope_mutation(
            "candidate embedded build type",
            lambda values: values[-1]["builds"]["candidate"][
                "validated_compile_commands"][
                    "effective_build_configuration"].update({
                        "embedded_build_type": "Debug"}))

        current_v16_scope = fake_evidence_scope()
        validate_evidence_scope(current_v16_scope)
        current_v16_semantics = current_v16_scope["builds"]["candidate"][
            "validated_compile_commands"]
        current_v16_configuration = current_v16_semantics[
            "effective_build_configuration"]
        current_v16_leopard2 = next(
            entry for entry in current_v16_semantics["required_entries"]
            if entry["file"].endswith("/leopard2.cpp"))
        current_v16_t32 = next(
            entry for entry in current_v16_semantics["required_entries"]
            if entry["file"].endswith("/Leopard2BackendAVX2T32B256.cpp"))
        require(
            current_v16_semantics["entry_count"] == 27 and
            len(current_v16_semantics["required_sources"]) == 19 and
            current_v16_semantics["schema"] == COMPILE_COMMANDS_SCHEMA and
            current_v16_semantics["profile"] == CANDIDATE_COMPILE_PROFILE and
            current_v16_configuration["entries"].get(
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT") == "ON" and
            current_v16_configuration["entries"].get(
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS") == "ON" and
            current_v16_configuration["entries"].get(
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK") == "OFF" and
            current_v16_configuration["entries"].get(
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED") == "ON" and
            current_v16_leopard2["arguments"].count(
                "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1") == 1 and
            current_v16_leopard2["arguments"].count(
                "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=1") == 1 and
            current_v16_leopard2["arguments"].count(
                "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=0") == 1 and
            current_v16_leopard2["arguments"].count(
                "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1") == 1 and
            current_v16_t32["arguments"].count(
                "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1") == 1 and
            any(path.endswith("/Leopard2BackendAVX2T8K8B1024.cpp")
                for path in current_v16_semantics["required_sources"]),
            "current v16 source/compile closure fixture differs")

        historical_v8_scope = fake_evidence_scope(
            compile_schema=COMPILE_COMMANDS_SCHEMA_V4)
        validate_evidence_scope(historical_v8_scope)
        historical_v8_configuration = historical_v8_scope[
            "builds"]["candidate"]["validated_compile_commands"][
                "effective_build_configuration"]
        parsed_historical_v8 = _parse_build_configuration_bytes(
            historical_v8_configuration["content"]["text"].encode("utf-8"))
        historical_record_schema, historical_variables, historical_cache = \
            _build_configuration_contract_for_compile_schema(
                COMPILE_COMMANDS_SCHEMA_V4,
                parsed_historical_v8["configuration_schema"])
        require(
            historical_record_schema == BUILD_CONFIGURATION_RECORD_SCHEMA_V2 and
            historical_variables == BUILD_CONFIGURATION_VARIABLES_V2 and
            historical_cache == REQUIRED_CANDIDATE_CACHE_V2 and
            parsed_historical_v8["entries"] ==
                historical_v8_configuration["entries"],
            "historical v8 effective-configuration fixture no longer replays")

        historical_v9_scope = fake_evidence_scope(
            compile_schema=COMPILE_COMMANDS_SCHEMA_V5)
        validate_evidence_scope(historical_v9_scope)
        historical_v9_configuration = historical_v9_scope[
            "builds"]["candidate"]["validated_compile_commands"][
                "effective_build_configuration"]
        parsed_historical_v9 = _parse_build_configuration_bytes(
            historical_v9_configuration["content"]["text"].encode("utf-8"))
        v9_record_schema, v9_variables, v9_cache = \
            _build_configuration_contract_for_compile_schema(
                COMPILE_COMMANDS_SCHEMA_V5,
                parsed_historical_v9["configuration_schema"])
        require(
            v9_record_schema == BUILD_CONFIGURATION_RECORD_SCHEMA_V3 and
            v9_variables == BUILD_CONFIGURATION_VARIABLES_V3 and
            v9_cache == REQUIRED_CANDIDATE_CACHE_V3 and
            parsed_historical_v9["entries"] ==
                historical_v9_configuration["entries"],
            "historical v9 effective-configuration fixture no longer replays")

        historical_v10_scope = fake_evidence_scope(
            compile_schema=COMPILE_COMMANDS_SCHEMA_V6)
        validate_evidence_scope(historical_v10_scope)
        historical_v10_configuration = historical_v10_scope[
            "builds"]["candidate"]["validated_compile_commands"][
                "effective_build_configuration"]
        parsed_historical_v10 = _parse_build_configuration_bytes(
            historical_v10_configuration["content"]["text"].encode("utf-8"))
        v10_record_schema, v10_variables, v10_cache = \
            _build_configuration_contract_for_compile_schema(
                COMPILE_COMMANDS_SCHEMA_V6,
                parsed_historical_v10["configuration_schema"])
        require(
            v10_record_schema == BUILD_CONFIGURATION_RECORD_SCHEMA_V4 and
            v10_variables == BUILD_CONFIGURATION_VARIABLES_V4 and
            v10_cache == REQUIRED_CANDIDATE_CACHE_V4 and
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED" not in
                parsed_historical_v10["entries"] and
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED" not in
                parsed_historical_v10["entries"],
            "historical v10 effective-configuration fixture no longer replays")

        historical_v11_scope = fake_evidence_scope(
            compile_schema=COMPILE_COMMANDS_SCHEMA_V7)
        validate_evidence_scope(historical_v11_scope)
        historical_v11_configuration = historical_v11_scope[
            "builds"]["candidate"]["validated_compile_commands"][
                "effective_build_configuration"]
        parsed_historical_v11 = _parse_build_configuration_bytes(
            historical_v11_configuration["content"]["text"].encode("utf-8"))
        v11_record_schema, v11_variables, v11_cache = \
            _build_configuration_contract_for_compile_schema(
                COMPILE_COMMANDS_SCHEMA_V7,
                parsed_historical_v11["configuration_schema"])
        require(
            v11_record_schema == BUILD_CONFIGURATION_RECORD_SCHEMA_V5 and
            v11_variables == BUILD_CONFIGURATION_VARIABLES_V5 and
            v11_cache == REQUIRED_CANDIDATE_CACHE_V5 and
            len(historical_v11_scope["builds"]["candidate"][
                "validated_compile_commands"]["required_sources"]) == 14,
            "historical v11 effective-configuration/source fixture no longer "
            "replays")

        historical_v12_scope = fake_evidence_scope(
            compile_schema=COMPILE_COMMANDS_SCHEMA_V8)
        validate_evidence_scope(historical_v12_scope)
        historical_v12_configuration = historical_v12_scope[
            "builds"]["candidate"]["validated_compile_commands"][
                "effective_build_configuration"]
        parsed_historical_v12 = _parse_build_configuration_bytes(
            historical_v12_configuration["content"]["text"].encode("utf-8"))
        v12_record_schema, v12_variables, v12_cache = \
            _build_configuration_contract_for_compile_schema(
                COMPILE_COMMANDS_SCHEMA_V8,
                parsed_historical_v12["configuration_schema"])
        historical_v12_sources = historical_v12_scope["builds"]["candidate"][
            "validated_compile_commands"]["required_sources"]
        require(
            v12_record_schema == BUILD_CONFIGURATION_RECORD_SCHEMA_V6 and
            v12_variables == BUILD_CONFIGURATION_VARIABLES_V6 and
            v12_cache == REQUIRED_CANDIDATE_CACHE_V6 and
            len(historical_v12_sources) == 18 and
            not any(path.endswith("/Leopard2BackendAVX2T8K8B1024.cpp")
                    for path in historical_v12_sources),
            "historical v12 effective-configuration/source fixture no longer "
            "replays")

        historical_v13_scope = fake_evidence_scope(
            compile_schema=COMPILE_COMMANDS_SCHEMA_V9)
        validate_evidence_scope(historical_v13_scope)
        historical_v13_semantics = historical_v13_scope[
            "builds"]["candidate"]["validated_compile_commands"]
        historical_v13_configuration = historical_v13_semantics[
            "effective_build_configuration"]
        parsed_historical_v13 = _parse_build_configuration_bytes(
            historical_v13_configuration["content"]["text"].encode("utf-8"))
        v13_record_schema, v13_variables, v13_cache = \
            _build_configuration_contract_for_compile_schema(
                COMPILE_COMMANDS_SCHEMA_V9,
                parsed_historical_v13["configuration_schema"])
        historical_v13_leopard2 = next(
            entry for entry in historical_v13_semantics["required_entries"]
            if entry["file"].endswith("/leopard2.cpp"))
        require(
            v13_record_schema == BUILD_CONFIGURATION_RECORD_SCHEMA_V7 and
            v13_variables == BUILD_CONFIGURATION_VARIABLES_V6 and
            v13_cache == REQUIRED_CANDIDATE_CACHE_V6 and
            historical_v13_semantics["profile"] ==
                CANDIDATE_COMPILE_PROFILE_V5 and
            len(historical_v13_semantics["required_sources"]) == 19 and
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT" not in
                parsed_historical_v13["entries"] and
            "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1" not in
                historical_v13_leopard2["arguments"],
            "historical v13 effective-configuration/source fixture no longer "
            "replays")

        historical_v14_scope = fake_evidence_scope(
            compile_schema=COMPILE_COMMANDS_SCHEMA_V10)
        validate_evidence_scope(historical_v14_scope)
        historical_v14_semantics = historical_v14_scope[
            "builds"]["candidate"]["validated_compile_commands"]
        historical_v14_configuration = historical_v14_semantics[
            "effective_build_configuration"]
        parsed_historical_v14 = _parse_build_configuration_bytes(
            historical_v14_configuration["content"]["text"].encode("utf-8"))
        v14_record_schema, v14_variables, v14_cache = \
            _build_configuration_contract_for_compile_schema(
                COMPILE_COMMANDS_SCHEMA_V10,
                parsed_historical_v14["configuration_schema"])
        historical_v14_leopard2 = next(
            entry for entry in historical_v14_semantics["required_entries"]
            if entry["file"].endswith("/leopard2.cpp"))
        require(
            v14_record_schema == BUILD_CONFIGURATION_RECORD_SCHEMA_V8 and
            v14_variables == BUILD_CONFIGURATION_VARIABLES_V8 and
            v14_cache == REQUIRED_CANDIDATE_CACHE_V8 and
            historical_v14_semantics["profile"] ==
                CANDIDATE_COMPILE_PROFILE_V6 and
            parsed_historical_v14["entries"].get(
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT") == "ON" and
            "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS" not in
                parsed_historical_v14["entries"] and
            "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK" not in
                parsed_historical_v14["entries"] and
            "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1" in
                historical_v14_leopard2["arguments"] and
            "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=1" not in
                historical_v14_leopard2["arguments"] and
            "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=0" not in
                historical_v14_leopard2["arguments"],
            "historical v14 effective-configuration/source fixture no longer "
            "replays")

        historical_v15_scope = fake_evidence_scope(
            compile_schema=COMPILE_COMMANDS_SCHEMA_V11)
        validate_evidence_scope(historical_v15_scope)
        historical_v15_semantics = historical_v15_scope[
            "builds"]["candidate"]["validated_compile_commands"]
        historical_v15_configuration = historical_v15_semantics[
            "effective_build_configuration"]
        parsed_historical_v15 = _parse_build_configuration_bytes(
            historical_v15_configuration["content"]["text"].encode("utf-8"))
        v15_record_schema, v15_variables, v15_cache = \
            _build_configuration_contract_for_compile_schema(
                COMPILE_COMMANDS_SCHEMA_V11,
                parsed_historical_v15["configuration_schema"])
        historical_v15_leopard2 = next(
            entry for entry in historical_v15_semantics["required_entries"]
            if entry["file"].endswith("/leopard2.cpp"))
        historical_v15_t32 = next(
            entry for entry in historical_v15_semantics["required_entries"]
            if entry["file"].endswith("/Leopard2BackendAVX2T32B256.cpp"))
        require(
            v15_record_schema == BUILD_CONFIGURATION_RECORD_SCHEMA and
            v15_variables == BUILD_CONFIGURATION_VARIABLES and
            v15_cache == REQUIRED_CANDIDATE_CACHE_V9 and
            historical_v15_semantics["profile"] ==
                CANDIDATE_COMPILE_PROFILE_V7 and
            parsed_historical_v15["entries"].get(
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED") == "OFF" and
            "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1" not in
                historical_v15_leopard2["arguments"] and
            "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1" not in
                historical_v15_t32["arguments"],
            "historical v15 selector/compile fixture no longer replays")

        for label, relabeled_scope, scope_schema, compile_version in (
            ("v16 body as v15", current_v16_scope,
             EVIDENCE_SCOPE_SCHEMA_V11, COMPILE_COMMANDS_SCHEMA_V11),
            ("v15 body as v16", historical_v15_scope,
             EVIDENCE_SCOPE_SCHEMA, COMPILE_COMMANDS_SCHEMA),
            ("v16 body as v14", current_v16_scope,
             EVIDENCE_SCOPE_SCHEMA_V10, COMPILE_COMMANDS_SCHEMA_V10),
            ("v14 body as v16", historical_v14_scope,
             EVIDENCE_SCOPE_SCHEMA, COMPILE_COMMANDS_SCHEMA),
            ("v16 body as v13", current_v16_scope,
             EVIDENCE_SCOPE_SCHEMA_V9, COMPILE_COMMANDS_SCHEMA_V9),
            ("v13 body as v16", historical_v13_scope,
             EVIDENCE_SCOPE_SCHEMA, COMPILE_COMMANDS_SCHEMA),
            ("v16 body as v12", current_v16_scope,
             EVIDENCE_SCOPE_SCHEMA_V8, COMPILE_COMMANDS_SCHEMA_V8),
            ("v12 body as v16", historical_v12_scope,
             EVIDENCE_SCOPE_SCHEMA, COMPILE_COMMANDS_SCHEMA),
            ("v16 body as v11", current_v16_scope,
             EVIDENCE_SCOPE_SCHEMA_V7, COMPILE_COMMANDS_SCHEMA_V7),
            ("v11 body as v16", historical_v11_scope,
             EVIDENCE_SCOPE_SCHEMA, COMPILE_COMMANDS_SCHEMA),
            ("v16 body as v10", current_v16_scope,
             EVIDENCE_SCOPE_SCHEMA_V6, COMPILE_COMMANDS_SCHEMA_V6),
            ("v10 body as v16", historical_v10_scope,
             EVIDENCE_SCOPE_SCHEMA, COMPILE_COMMANDS_SCHEMA),
            ("v16 body as v9", current_v16_scope,
             EVIDENCE_SCOPE_SCHEMA_V5, COMPILE_COMMANDS_SCHEMA_V5),
            ("v9 body as v16", historical_v9_scope,
             EVIDENCE_SCOPE_SCHEMA, COMPILE_COMMANDS_SCHEMA),
        ):
            relabeled = json.loads(json.dumps(relabeled_scope))
            relabeled["schema"] = scope_schema
            for relabeled_build in relabeled["builds"].values():
                relabeled_build["validated_compile_commands"]["schema"] = \
                    compile_version
            try:
                validate_evidence_scope(relabeled)
            except PlanError:
                pass
            else:
                raise PlanError(
                    f"evidence scope accepted coherent {label} relabel")

        def coherently_enable_small_direct(values) -> None:
            value = values[-1]
            build = value["builds"]["candidate"]
            semantics = build["validated_compile_commands"]
            configuration = semantics["effective_build_configuration"]
            configuration["entries"][
                "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"] = "1"
            material = _build_configuration_material(configuration["entries"])
            digest = hashlib.sha256(material).hexdigest()
            text = (
                f"schema={BUILD_CONFIGURATION_FILE_SCHEMA}\n"
                f"sha256={digest}\n").encode("ascii") + material
            content = retained_text(text.decode("utf-8"))
            configuration["content"] = content
            configuration["artifact"].update({
                "size": content["size"], "sha256": content["sha256"],
            })
            configuration["configuration_sha256"] = digest
            build["validated_cache"][
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"] = digest
            benchmark = next(
                entry for entry in semantics["required_entries"]
                if entry["file"].endswith("/bench/leopard2/benchmark.cpp"))
            benchmark["arguments"] = [
                (f'-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256="{digest}"'
                 if token.startswith(
                     "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=") else
                 token)
                for token in benchmark["arguments"]
            ]

        reject_scope_mutation(
            "coherent noncanonical small-direct configuration",
            coherently_enable_small_direct)

        def coherently_disable_small_dual_direct(values) -> None:
            value = values[-1]
            build = value["builds"]["candidate"]
            semantics = build["validated_compile_commands"]
            configuration = semantics["effective_build_configuration"]
            configuration["entries"][
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT"] = "OFF"
            build["validated_cache"][
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT"] = "OFF"
            material = _build_configuration_material(configuration["entries"])
            digest = hashlib.sha256(material).hexdigest()
            text = (
                f"schema={BUILD_CONFIGURATION_FILE_SCHEMA}\n"
                f"sha256={digest}\n").encode("ascii") + material
            content = retained_text(text.decode("utf-8"))
            configuration["content"] = content
            configuration["artifact"].update({
                "size": content["size"], "sha256": content["sha256"],
            })
            configuration["configuration_sha256"] = digest
            build["validated_cache"][
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"] = digest
            for entry in semantics["required_entries"]:
                entry["arguments"] = [
                    ("-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=0"
                     if token ==
                        "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1" else token)
                    for token in entry["arguments"]
                ]
                if entry["file"].endswith("/bench/leopard2/benchmark.cpp"):
                    entry["arguments"] = [
                        (f'-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256="{digest}"'
                         if token.startswith(
                             "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=")
                         else token)
                        for token in entry["arguments"]
                    ]

        reject_scope_mutation(
            "coherent disabled small-dual-direct configuration",
            coherently_disable_small_dual_direct)
        for label, flags in (
            ("final optimization downgrade", "-O3 -O0"),
            ("bare optimization downgrade", "-O3 -O"),
            ("unknown optimization spelling", "-O4 -O3"),
            ("sanitizer after optimization", "-O3 -fsanitize=address"),
            ("profile after optimization", "-O3 -fprofile-generate"),
            ("LTO after optimization", "-O3 -flto"),
            ("instrumentation after optimization",
             "-O3 -finstrument-functions"),
            ("vector disable after optimization",
             "-O3 -fno-tree-vectorize"),
            ("coverage after optimization", "-O3 --coverage"),
            ("long optimize alias", "-O3 --optimize=0"),
            ("long sanitizer alias", "-O3 --sanitize=address"),
            ("long instrumentation alias", "-O3 --instrument-functions"),
            ("long LTO alias", "-O3 --lto"),
            ("long profile alias", "-O3 --profile"),
            ("inline disable", "-O3 -fno-inline"),
            ("GCC pass disable", "-O3 -fdisable-tree-vect"),
            ("Release response file", "-O3 @evil.rsp"),
        ):
            reject_scope_mutation(
                f"uniform CMake {label}",
                lambda values, item=flags: [
                    value["builds"][role]["validated_cache"].update({
                        "CMAKE_CXX_FLAGS_RELEASE": item})
                    for value in values
                    for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform mislabeled candidate compile profile",
            lambda values: [
                value["builds"]["candidate"][
                    "validated_compile_commands"].update({
                        "profile": BASELINE_COMPILE_PROFILE})
                for value in values])
        reject_scope_mutation(
            "uniform mislabeled baseline compile implementation",
            lambda values: [
                value["builds"]["baseline"][
                    "validated_compile_commands"].update({
                        "implementation": "candidate"})
                for value in values])
        reject_scope_mutation(
            "uniform obsolete compile-command schema",
            lambda values: [
                value["builds"][role][
                    "validated_compile_commands"].update({
                        "schema": "leopard2-main-compare-compile-commands/v1"})
                for value in values for role in ("baseline", "candidate")])

        def mutate_all_compiler_argv(
            values, role: str, suffix: str, mutation: str,
        ) -> None:
            for value in values:
                semantics = value["builds"][role]["validated_compile_commands"]
                entry = next(item for item in semantics["required_entries"]
                             if item["file"].endswith(suffix))
                arguments = entry["arguments"]
                source_index = arguments.index(entry["file"])
                output_index = arguments.index("-o")
                compile_index = arguments.index("-c")
                other = next(item for item in semantics["required_entries"]
                             if item["file"] != entry["file"])
                if mutation == "response":
                    arguments.insert(source_index, "@evil.rsp")
                elif mutation == "extra_source":
                    arguments.insert(source_index, other["file"])
                elif mutation == "stdin_source":
                    arguments.insert(source_index, "-")
                elif mutation == "wrong_source":
                    arguments[source_index] = other["file"]
                elif mutation == "missing_source":
                    arguments.pop(source_index)
                elif mutation == "duplicate_source":
                    arguments.append(entry["file"])
                elif mutation == "wrong_output":
                    arguments[output_index + 1] = other["output"]
                elif mutation == "coherent_output":
                    entry["output"] = other["output"]
                    arguments[output_index + 1] = other["output"]
                elif mutation == "duplicate_output":
                    arguments[compile_index:compile_index] = [
                        "-o", entry["output"]]
                elif mutation == "missing_compile":
                    arguments.pop(compile_index)
                elif mutation == "duplicate_compile":
                    arguments.insert(compile_index, "-c")
                elif mutation == "extra_definition":
                    arguments.insert(1, "-DEVIL=1")
                elif mutation == "reorder_definitions":
                    definitions = [
                        index for index, token in enumerate(arguments)
                        if token.startswith("-D")]
                    arguments[definitions[0]], arguments[definitions[1]] = \
                        arguments[definitions[1]], arguments[definitions[0]]
                elif mutation == "extra_include":
                    include = next(
                        index for index, token in enumerate(arguments)
                        if token.startswith("-I"))
                    arguments.insert(include + 1, "-I/tmp/evil")
                elif mutation == "different_include":
                    include = next(
                        index for index, token in enumerate(arguments)
                        if token.startswith("-I"))
                    arguments[include] = "-I/tmp/evil"
                elif mutation == "language_last":
                    arguments.insert(
                        arguments.index("-std=gnu++11") + 1,
                        "-std=gnu++20")
                elif mutation == "optimization_last":
                    arguments.insert(output_index, "-O0")
                elif mutation == "avx2_last":
                    arguments.insert(output_index, "-mno-avx2")
                elif mutation == "missing_avx2":
                    arguments.remove("-mavx2")
                elif mutation == "baseline_isa_last":
                    arguments.insert(output_index, "-march=x86-64")
                elif mutation == "compiler_wrapper":
                    arguments.insert(1, "--compiler-wrapper=/tmp/evil")
                else:
                    raise PlanError(
                        f"unknown compiler self-test mutation: {mutation}")

        for label, role, suffix, mutation in (
            ("compiler response file", "candidate", "/leopard.cpp", "response"),
            ("compiler second positional source", "candidate", "/leopard.cpp", "extra_source"),
            ("compiler stdin positional source", "candidate", "/leopard.cpp", "stdin_source"),
            ("compiler different positional source", "candidate", "/leopard.cpp", "wrong_source"),
            ("compiler missing positional source", "candidate", "/leopard.cpp", "missing_source"),
            ("compiler duplicate positional source", "candidate", "/leopard.cpp", "duplicate_source"),
            ("compiler wrong output", "candidate", "/leopard.cpp", "wrong_output"),
            ("compiler coherent wrong output", "candidate", "/leopard.cpp", "coherent_output"),
            ("compiler duplicate output", "candidate", "/leopard.cpp", "duplicate_output"),
            ("compiler missing compile option", "candidate", "/leopard.cpp", "missing_compile"),
            ("compiler duplicate compile option", "candidate", "/leopard.cpp", "duplicate_compile"),
            ("compiler extra definition", "candidate", "/leopard.cpp", "extra_definition"),
            ("compiler reordered definitions", "candidate", "/leopard.cpp", "reorder_definitions"),
            ("compiler extra include", "candidate", "/leopard.cpp", "extra_include"),
            ("compiler different include", "candidate", "/leopard.cpp", "different_include"),
            ("compiler last language option", "candidate", "/leopard.cpp", "language_last"),
            ("compiler last optimization option", "candidate", "/leopard.cpp", "optimization_last"),
            ("compiler last AVX2 option", "candidate", "/Leopard2BackendAVX2.cpp", "avx2_last"),
            ("compiler missing AVX2 option", "candidate", "/Leopard2BackendAVX2.cpp", "missing_avx2"),
            ("compiler last baseline ISA option", "baseline", "/leopard.cpp", "baseline_isa_last"),
            ("compiler wrapper control", "candidate", "/leopard.cpp", "compiler_wrapper"),
        ):
            reject_scope_mutation(
                "uniform " + label,
                lambda values, item_role=role, item_suffix=suffix,
                       item_mutation=mutation:
                    mutate_all_compiler_argv(
                        values, item_role, item_suffix, item_mutation))
        reject_scope_mutation("uniform incomplete sources", lambda values:
                              [value["sources"][role].pop("tree")
                               for value in values
                               for role in ("baseline", "candidate")])
        reject_scope_mutation("uniform incomplete runtime files", lambda values:
                              [dependency.update({
                                   "file": {"path": dependency["file"]["path"]}})
                               for value in values
                               for role in ("baseline", "candidate")
                               for dependency in value["runtime_closures"][role][
                                   "dependencies"]
                               if "file" in dependency])
        reject_scope_mutation("uniform path-only compile pairs", lambda values:
                              [pair.update({
                                   "source": {"path": pair["source"]["path"]},
                                   "object": {"path": pair["object"]["path"]}})
                               for value in values
                               for role in ("baseline", "candidate")
                               for pair in value["builds"][role][
                                   "validated_compile_commands"][
                                       "required_source_object_pairs"]])
        reject_scope_mutation("uniform empty compiler/CMake/link records",
                              lambda values:
                              [value["builds"][role].update({
                                   name: {} for name in (
                                       "compiler", "cmake_cache", "compile_commands",
                                       "executable_link_recipe",
                                       "archive_link_recipe")})
                               for value in values
                               for role in ("baseline", "candidate")])
        def reduce_all_outputs(values) -> None:
            for value in values:
                for role in ("baseline", "candidate"):
                    for output in ("archive", "executable"):
                        key = f"{role}_{output}"
                        retained = value["artifacts"][key]
                        reduced = {name: retained[name]
                                   for name in ("path", "kind", "sha256")}
                        value["artifacts"][key] = dict(reduced)
                        value["builds"][role][f"validated_{output}"] = dict(reduced)
        reject_scope_mutation("uniform reduced archive/executable",
                              reduce_all_outputs)
        def reduce_all_hosts_to_topology(values) -> None:
            for value in values:
                host = value["host"]
                for name in ("benchmark_cpu", "reserved_sibling"):
                    retained = host[name]
                    host[name] = {"cpu": retained["cpu"],
                                  "topology": retained["topology"]}
        reject_scope_mutation("uniform topology-only host",
                              reduce_all_hosts_to_topology)
        def truncate_all_translation_unit_closures(values) -> None:
            for value in values:
                for role, names in (
                    ("baseline", BASELINE_LIBRARY_SOURCES),
                    ("candidate", CANDIDATE_LIBRARY_SOURCES),
                ):
                    build = value["builds"][role]
                    semantics = build["validated_compile_commands"]
                    pair = next(item for item in
                                semantics["required_source_object_pairs"]
                                if item["source"]["path"].endswith(
                                    "/" + names[-1]))
                    semantics["required_source_object_pairs"].remove(pair)
                    semantics["required_sources"].remove(pair["source"]["path"])
                    semantics["entry_count"] -= 1
                    member = Path(pair["object"]["path"]).name
                    build["validated_archive_members"].remove(member)
                    relative = pair["object"]["path"][
                        len(build["build_dir"]) + 1:]
                    text = build["archive_link_recipe_content"]["text"].replace(
                        " " + relative, "", 1)
                    retained = {
                        "encoding": "utf-8", "text": text,
                        "size": len(text.encode("utf-8")),
                        "sha256": hashlib.sha256(
                            text.encode("utf-8")).hexdigest(),
                    }
                    build["archive_link_recipe_content"] = retained
                    build["archive_link_recipe"]["size"] = retained["size"]
                    build["archive_link_recipe"]["sha256"] = retained["sha256"]
        reject_scope_mutation(
            "uniform coherent translation-unit truncation",
            truncate_all_translation_unit_closures)
        def redirect_all_declared_archives(values) -> None:
            for value in values:
                for role, archive_name in (
                    ("baseline", "libleopard_main_exact.a"),
                    ("candidate", "libleopard.a"),
                ):
                    build = value["builds"][role]
                    text = build[
                        "executable_link_recipe_content"]["text"].replace(
                            f" {archive_name} ",
                            f" /tmp/{archive_name} ", 1)
                    retained = retained_text(text)
                    build["executable_link_recipe_content"] = retained
                    build["executable_link_recipe"]["size"] = retained["size"]
                    build["executable_link_recipe"]["sha256"] = retained["sha256"]
        reject_scope_mutation(
            "uniform same-basename foreign archive operands",
            redirect_all_declared_archives)
        def add_all_foreign_archive_operands(values) -> None:
            for value in values:
                for role, archive_name in (
                    ("baseline", "libleopard_main_exact.a"),
                    ("candidate", "libleopard.a"),
                ):
                    build = value["builds"][role]
                    text = build[
                        "executable_link_recipe_content"]["text"].replace(
                            f" {archive_name} ",
                            f" {archive_name} /tmp/evil.a ", 1)
                    retained = retained_text(text)
                    build["executable_link_recipe_content"] = retained
                    build["executable_link_recipe"]["size"] = retained["size"]
                    build["executable_link_recipe"]["sha256"] = retained["sha256"]
        reject_scope_mutation(
            "uniform foreign archive operands",
            add_all_foreign_archive_operands)
        def add_all_link_controls(values, control: str) -> None:
            for value in values:
                for role, archive_name in (
                    ("baseline", "libleopard_main_exact.a"),
                    ("candidate", "libleopard.a"),
                ):
                    build = value["builds"][role]
                    text = build[
                        "executable_link_recipe_content"]["text"].replace(
                            f" {archive_name} ",
                            f" {archive_name} {control} ", 1)
                    retained = retained_text(text)
                    build["executable_link_recipe_content"] = retained
                    build["executable_link_recipe"]["size"] = retained["size"]
                    build["executable_link_recipe"]["sha256"] = retained["sha256"]
        for label, control in (
            ("fused linker script", "-Wl,--script=/tmp/evil.ld"),
            ("fused search/library", "-Wl,-L,/tmp,-levil"),
            ("driver response", "@evil.rsp"),
            ("compiler specs", "-specs=/tmp/evil.specs"),
            ("alternate tool root", "-B/tmp/toolchain"),
            ("compiler plugin", "-fplugin=/tmp/evil.so"),
            ("linker plugin", "-Wl,--plugin,/tmp/evil.so"),
            ("alternate linker", "-fuse-ld=/tmp/evil-ld"),
            ("alternate runtime loader",
             "-Wl,--dynamic-linker,/tmp/evil-ld.so"),
        ):
            reject_scope_mutation(
                f"uniform {label}",
                lambda values, item=control:
                    add_all_link_controls(values, item))
        def append_all_executable_commands(values) -> None:
            for value in values:
                for role in ("baseline", "candidate"):
                    build = value["builds"][role]
                    text = build[
                        "executable_link_recipe_content"]["text"] + "\n-O3\n"
                    retained = retained_text(text)
                    build["executable_link_recipe_content"] = retained
                    build["executable_link_recipe"]["size"] = retained["size"]
                    build["executable_link_recipe"]["sha256"] = retained["sha256"]
        reject_scope_mutation(
            "uniform multiline executable recipes",
            append_all_executable_commands)
        def redirect_all_pthread_inputs(values) -> None:
            old = "/usr/lib/x86_64-linux-gnu/libpthread.a"
            new = "/tmp/libpthread.a"
            for value in values:
                for role in ("baseline", "candidate"):
                    build = value["builds"][role]
                    record = build["validated_external_link_inputs"][1]
                    record["operand"] = new
                    record["artifact"]["path"] = new
                    text = build[
                        "executable_link_recipe_content"]["text"].replace(
                            old, new, 1)
                    retained = retained_text(text)
                    build["executable_link_recipe_content"] = retained
                    build["executable_link_recipe"]["size"] = retained["size"]
                    build["executable_link_recipe"]["sha256"] = retained["sha256"]
        reject_scope_mutation(
            "uniform alternate pthread root", redirect_all_pthread_inputs)
        reject_scope_mutation(
            "uniform missing external static identity",
            lambda values: [
                value["builds"][role]["validated_external_link_inputs"].pop()
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "candidate external identity differs from baseline",
            lambda values: [
                value["builds"]["candidate"][
                    "validated_external_link_inputs"][1]["artifact"].update({
                        "sha256": "e" * 64})
                for value in values])
        reject_scope_mutation(
            "uniform external resolved path escape",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][1]["artifact"].update({
                        "path": "/tmp/libpthread.a"})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform nonexistent versioned OpenMP target",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0]["artifact"].update({
                        "path": "/usr/lib/x86_64-linux-gnu/"
                                "libgomp.so.999.999"})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform differing external roles",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][1].update({
                        "role": "openmp_runtime_shared"})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform external file mode drift",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0]["artifact"].update({
                        "mode": value["builds"][role][
                            "validated_external_link_inputs"][0]["artifact"][
                                "mode"] ^ 0o100})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform external lexical-link target drift",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0][
                        "lexical_symlink_chain"][0].update({
                            "target": "../../../x86_64-linux-gnu/"
                                      "libgomp.so.999"})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform external lexical-link mode drift",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0][
                        "lexical_symlink_chain"][0].update({"mode": 0o700})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform truncated external lexical-link chain",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0][
                        "lexical_symlink_chain"].pop()
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform fabricated pthread lexical-link chain",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][1].update({
                        "lexical_symlink_chain": [{
                            "path": "/usr/lib/x86_64-linux-gnu/libpthread.a",
                            "target": "libpthread-evil.a", "mode": 0o777,
                        }]})
                for value in values for role in ("baseline", "candidate")])
        for label, control in (
            ("recipe optimization downgrade", "-O0"),
            ("recipe sanitizer", "-fsanitize=address"),
        ):
            reject_scope_mutation(
                f"uniform {label}",
                lambda values, item=control:
                    add_all_link_controls(values, item))
        def remove_all_dynamic_loader_records(values) -> None:
            for value in values:
                for role in ("baseline", "candidate"):
                    closure = value["runtime_closures"][role]
                    closure["dependencies"] = [
                        item for item in closure["dependencies"]
                        if item["soname"] != "ld-linux-x86-64.so.2"]
                    text = "\n".join(
                        line for line in runtime_text(closure).splitlines()
                        if not line.lstrip().startswith(
                            "/lib/ld-linux-x86-64.so.2")) + "\n"
                    replace_runtime_text(closure, text)
        reject_scope_mutation(
            "uniform missing dynamic-loader records",
            remove_all_dynamic_loader_records)
        coherently_rewritten_runtime_scopes = json.loads(json.dumps(scopes))
        for value in coherently_rewritten_runtime_scopes:
            for role in ("baseline", "candidate"):
                closure = value["runtime_closures"][role]
                closure["dependencies"] = [
                    dependency for dependency in closure["dependencies"]
                    if dependency["soname"] != "libm.so.6"]
                text = "\n".join(
                    line for line in runtime_text(closure).splitlines()
                    if not line.startswith("libm.so.6 =>")) + "\n"
                replace_runtime_text(closure, text)
        derive_survivors(
            first, manifests, references, coherently_rewritten_runtime_scopes)
        def swap_all_runtime_file_records(values) -> None:
            for value in values:
                for role in ("baseline", "candidate"):
                    dependencies = value["runtime_closures"][role]["dependencies"]
                    libc = next(item for item in dependencies
                                if item["soname"] == "libc.so.6")
                    libm = next(item for item in dependencies
                                if item["soname"] == "libm.so.6")
                    libc["file"], libm["file"] = libm["file"], libc["file"]
        reject_scope_mutation(
            "uniform swapped runtime file records",
            swap_all_runtime_file_records)
        def make_all_runtime_loader_paths_noncanonical(values) -> None:
            for value in values:
                for role in ("baseline", "candidate"):
                    closure = value["runtime_closures"][role]
                    libc = next(item for item in closure["dependencies"]
                                if item["soname"] == "libc.so.6")
                    libc["loader_path"] = "/lib/./libc.so.6"
                    libc["file"]["path"] = libc["loader_path"]
                    text = runtime_text(closure).replace(
                        "/lib/libc.so.6", libc["loader_path"])
                    replace_runtime_text(closure, text)
        reject_scope_mutation(
            "uniform noncanonical runtime loader paths",
            make_all_runtime_loader_paths_noncanonical)
        def truncate_all_cache_records(values) -> None:
            for value in values:
                for name in ("benchmark_cpu", "reserved_sibling"):
                    value["host"][name]["cache_hierarchy"].pop()
                    value["host"][name]["cache_index_inventory"].pop()
        reject_scope_mutation(
            "uniform cache summary/listing mismatch", truncate_all_cache_records)
        coherently_rewritten_cache_scopes = json.loads(json.dumps(scopes))
        for value in coherently_rewritten_cache_scopes:
            for name in ("benchmark_cpu", "reserved_sibling"):
                record = value["host"][name]
                record["cache_hierarchy"].pop()
                record["cache_index_inventory"].pop()
                record["cache_directory_inventory_text"] = retained_text(
                    "".join(f"{entry}\n"
                            for entry in record["cache_index_inventory"]))
        derive_survivors(
            first, manifests, references, coherently_rewritten_cache_scopes)
        def empty_all_numa_records(values) -> None:
            for value in values:
                for name in ("benchmark_cpu", "reserved_sibling"):
                    value["host"][name]["numa_nodes"] = []
                    value["host"][name]["numa_node_inventory"] = []
        reject_scope_mutation(
            "uniform empty NUMA records", empty_all_numa_records)
        reject_scope_mutation(
            "uniform source tree/commit mismatch",
            lambda values: [
                value["sources"][role].update({"tree": "f" * 40})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation("baseline archive", lambda values:
                              values[-1]["artifacts"]["baseline_archive"].update({
                                  "sha256": "6" * 64}))
        reject_scope_mutation("baseline executable", lambda values:
                              values[-1]["artifacts"]["baseline_executable"].update({
                                  "sha256": "7" * 64}))
        reject_scope_mutation("baseline CMake cache", lambda values:
                              values[-1]["builds"]["baseline"][
                                  "cmake_cache"].update({"sha256": "8" * 64}))
        reject_scope_mutation("baseline archive recipe", lambda values:
                              values[-1]["builds"]["baseline"][
                                  "archive_link_recipe_content"].update({
                                      "sha256": "9" * 64}))
        reject_scope_mutation("baseline compile closure", lambda values:
                              values[-1]["builds"]["baseline"][
                                  "validated_compile_commands"][
                                      "required_source_object_pairs"][0][
                                          "object"].update({"sha256": "0" * 64}))
        reject_scope_mutation("baseline runtime dependency", lambda values:
                              values[-1]["runtime_closures"]["baseline"][
                                  "dependencies"][0]["file"].update({
                                      "sha256": "6" * 64}))
        reject_scope_mutation("candidate runtime dependency", lambda values:
                              values[-1]["runtime_closures"]["candidate"][
                                  "dependencies"][0]["file"].update({
                                      "sha256": "7" * 64}))
        reject_scope_mutation("uniform noncanonical baseline CMake", lambda values:
                              [value["builds"]["baseline"][
                                  "validated_cache"].update({
                                      "LEO_MAIN_HAS_MARCH_NATIVE": "0"})
                               for value in values])
        reject_scope_mutation("resolved backend", lambda values:
                              values[-1].update({
                                  "resolved_auto_backend": "ssse3"}))
        for excluded_backend in ("avx512", "neon"):
            reject_scope_mutation(
                "mixed " + excluded_backend,
                lambda values, backend=excluded_backend:
                    values[-1].update({"resolved_auto_backend": backend}))
            reject_scope_mutation(excluded_backend, lambda values, backend=excluded_backend:
                                  [value.update({
                                      "resolved_auto_backend": backend})
                                   for value in values])

        normalized_nested = _normalize_bound_paths(
            {"paths": ["/tree/build/object.o", "/tree/source.cpp"]},
            (("/tree", "$SOURCE"), ("/tree/build", "$BUILD")))
        require(normalized_nested == {
            "paths": ["$BUILD/object.o", "$SOURCE/source.cpp"]},
            "scope path normalization did not prefer the longest role root")
        normalized_collisions = _normalize_bound_paths({
            "exact": "/tree",
            "child": "/tree/child",
            "siblings": [
                "/tree-build/object.o", "/treehouse/object.o",
                "/tree.build/object.o", "/tree!object.o", "/tree#object.o",
                "/tree?object.o", "/other/tree/object.o",
            ],
            "text": ("'/tree/file' /tree-build/file /treehouse/file "
                     "/tree.build/file /tree!file /tree#file /tree?file "
                     "/other/tree/file root=/tree (/tree),/tree:/next"),
        }, (("/tree", "$ROOT"),))
        require(normalized_collisions == {
            "exact": "$ROOT",
            "child": "$ROOT/child",
            "siblings": [
                "/tree-build/object.o", "/treehouse/object.o",
                "/tree.build/object.o", "/tree!object.o", "/tree#object.o",
                "/tree?object.o", "/other/tree/object.o",
            ],
            "text": ("'$ROOT/file' /tree-build/file /treehouse/file "
                     "/tree.build/file /tree!file /tree#file /tree?file "
                     "/other/tree/file root=$ROOT ($ROOT),$ROOT:/next"),
        }, "scope path normalization rewrote a textual sibling prefix")

        def normalize_cmake_external_object(
            source_root: str, build_root: str,
        ) -> dict[str, Any]:
            object_path = (
                f"{build_root}/CMakeFiles/leopard.dir"
                f"{source_root}/LeopardFF8.cpp.o")
            recipe_text = (
                "ar qc libleopard.a CMakeFiles/leopard.dir"
                f"{source_root}/LeopardFF8.cpp.o\n")
            return _normalize_bound_paths({
                "object_path": object_path,
                "archive_link_recipe_content": retained_text(recipe_text),
            }, ((source_root, "$SOURCE"), (build_root, "$BUILD")))

        normalized_home_source = normalize_cmake_external_object(
            "/home/a", "/opt/b")
        normalized_opt_source = normalize_cmake_external_object(
            "/opt/b", "/home/a")
        expected_external_recipe = retained_text(
            "ar qc libleopard.a CMakeFiles/leopard.dir/"
            "$SOURCE/LeopardFF8.cpp.o\n")
        expected_external_object = {
            "archive_link_recipe_content": expected_external_recipe,
            "object_path": (
                "$BUILD/CMakeFiles/leopard.dir/"
                "$SOURCE/LeopardFF8.cpp.o"),
        }
        require(normalized_home_source == expected_external_object and
                normalized_opt_source == expected_external_object,
                "CMake external-source object paths are not root-independent")
        normalized_external_sibling = _normalize_bound_paths({
            "object_path": (
                "/opt/b/CMakeFiles/leopard.dir/"
                "home/a-copy/LeopardFF8.cpp.o"),
        }, (("/home/a", "$SOURCE"), ("/opt/b", "$BUILD")))
        require(normalized_external_sibling == {
            "object_path": (
                "$BUILD/CMakeFiles/leopard.dir/"
                "home/a-copy/LeopardFF8.cpp.o"),
        }, "CMake external-source handling rewrote a sibling path prefix")
        recipe_text = "/tree/build/ar qc /tree/source.cpp\n"
        normalized_recipe = _normalize_bound_paths({
            "archive_link_recipe": {
                "path": "/tree/build/link.txt", "kind": "build_metadata",
                "size": len(recipe_text.encode("utf-8")), "mode": 0o644,
                "sha256": hashlib.sha256(recipe_text.encode("utf-8")).hexdigest(),
                "mtime_ns": 1,
            },
            "archive_link_recipe_content": {
                "encoding": "utf-8", "text": recipe_text,
                "size": len(recipe_text.encode("utf-8")),
                "sha256": hashlib.sha256(recipe_text.encode("utf-8")).hexdigest(),
            },
        }, (("/tree", "$SOURCE"), ("/tree/build", "$BUILD")))
        normalized_content = normalized_recipe["archive_link_recipe_content"]
        require(normalized_content["text"] ==
                    "$BUILD/ar qc $SOURCE/source.cpp\n" and
                normalized_content["sha256"] == hashlib.sha256(
                    normalized_content["text"].encode("utf-8")).hexdigest() and
                normalized_recipe["archive_link_recipe"]["sha256"] ==
                    normalized_content["sha256"] and
                normalized_recipe["archive_link_recipe"]["size"] ==
                    normalized_content["size"],
                "normalized recipe content retained a stale byte identity")

        unsolicited_manifests = json.loads(json.dumps(manifests))
        unsolicited = main_cell(
            "gate-k80-b65536", 80, 80, 65536, 80,
            "exact-main-generic-gate")
        unsolicited_manifests[-1]["campaign"]["cells"].append({
            "identifier": unsolicited["identifier"], "k": unsolicited["K"],
            "r": unsolicited["R"], "shard_bytes": unsolicited["shard_bytes"],
            "losses": unsolicited["loss_count"], "seed": unsolicited["seed"],
        })
        unsolicited_manifests[-1]["analysis"][unsolicited["identifier"]] = {
            metric: {
                "ci95_lower": 1.06,
                "ratio_orientation": "baseline_time_over_candidate_time",
            }
            for metric in ("decode_first_use", "decode_reuse_amortized")
        }
        try:
            derive_survivors(first, unsolicited_manifests, references, scopes)
        except PlanError:
            pass
        else:
            raise PlanError("an unrequested supplemental cell entered selection")

        refinement_manifests = json.loads(json.dumps(manifests))
        for manifest in refinement_manifests:
            for cell in manifest["campaign"]["cells"]:
                if cell["k"] == 112 and cell["shard_bytes"] == 65536:
                    for metric in ("decode_first_use", "decode_reuse_amortized"):
                        manifest["analysis"][cell["identifier"]][metric][
                            "ci95_lower"] = 1.06
        refinement = derive_survivors(
            first, refinement_manifests, references, scopes)
        expected_refinement_k = list(range(97, 112)) + list(range(113, 120))
        require([cell["K"] for cell in
                 refinement["required_refinement_cells"]] == expected_refinement_k,
                "transition refinement cells differ")
        try:
            materialize_stage(first, refinement, root / "premature-stage")
        except PlanError:
            pass
        else:
            raise PlanError("stage advanced before K transition refinement")

        survivor_path = root / "survivors.json"
        write_json(survivor_path, survivor_signed)
        validate_survivors(first, survivor_path, manifests, scopes)

        forged_survivor = root / "forged-survivor.json"
        shutil.copy2(survivor_path, forged_survivor)
        def add_rejected(value):
            rejected = value["rejected_cells"][0]
            value["survivor_cells"].append(rejected)
        adversarial_resign(forged_survivor, add_rejected)
        try:
            validate_survivors(first, forged_survivor, manifests, scopes)
        except PlanError:
            pass
        else:
            raise PlanError("re-signed non-survivor entered the survivor set")

        stage_root = root / "stage"
        materialize_stage(first, survivor_signed, stage_root)
        stage = validate_stage(first, stage_root, manifests, scopes)
        require(stage["survivor_K"] == [127, 128] and
                stage["promotion_requires_path_attestation"] is True and
                stage["expected_resolved_backend"] == "avx2" and
                stage["forced_confirmation_backends"] == list(BACKENDS) and
                set(stage["excluded_backends"]) == {"avx512", "neon"},
                "stage survivor/path-attestation semantics differ")
        survivor_pairs = {(128, 4096), (127, 65536)}
        for path in (stage_root / "forced-survivors").glob("*.json"):
            _, cases = common.normalize_matrix(load_json(path))
            require(all((case["K"], case["shard_bytes"]) in survivor_pairs
                        for case in cases) and
                    {(case["K"], case["shard_bytes"]) for case in cases}
                        .issubset(survivor_pairs),
                    "a non-surviving K/byte entered forced comparison")
        for path in (stage_root / "forced-true-tails").glob("*.json"):
            _, cases = common.normalize_matrix(load_json(path))
            require(all(case["K"] in {127, 128} and
                        case["shard_bytes"] in TRUE_TAIL_BYTES and
                        case["shard_bytes"] % 64 != 0 for case in cases),
                    "aligned or non-surviving cell entered true-tail stage")
        rejection = validate_main_document(
            load_json(stage_root / "exact-main-auto-rejection-timing.json"),
            "rejection timing")
        attestation = unsigned(load_json(stage_root / "path-attestation.json"),
                               ATTESTATION_SCHEMA, "path attestation")
        positive_pairs = {
            (item["cell"]["K"], item["cell"]["shard_bytes"])
            for item in attestation["cases"]
            if item["expected_selected_decode_path"] == "generic"
        }
        require(positive_pairs == survivor_pairs and all(
            item["expected_selected_decode_rule"] == "balanced_generic"
            for item in attestation["cases"]
            if item["expected_selected_decode_path"] == "generic"),
            "attestation broadened the exact surviving gate cells")
        negative = [item for item in attestation["cases"]
                    if item["expected_selected_decode_path"] == "not_generic"]
        require(any(item["cell"]["category"] == "rejected_gate" and
                    item["cell"]["K"] == 128 and
                    item["cell"]["shard_bytes"] == 65536 and
                    item["cell"]["loss_count"] == 128 for item in negative) and
                any(item["cell"]["category"] == "rejected_gate" and
                    item["cell"]["K"] == 127 and
                    item["cell"]["shard_bytes"] == 4096 and
                    item["cell"]["loss_count"] == 127 for item in negative) and
                all(item["expected_selected_decode_rule"] ==
                    "not_balanced_generic" for item in negative),
                "rejected mixed-K full-loss gates were not fail-closed")
        require({tuple((item["cell"][key] for key in (
                    "K", "R", "shard_bytes", "loss_count", "seed")))
                 for item in attestation["cases"]
                 if item["cell"]["category"] == "loss_rate_neighbor"} ==
                {tuple((cell[key] for key in (
                    "K", "R", "shard_bytes", "loss_count", "seed")))
                 for cell in rejection["cells"]},
                "attestation loss/rate neighbors differ from rejection timing")
        require(all(item["expected_selected_decode_path"] == "not_generic"
                    for item in attestation["cases"]
                    if item["cell"]["category"] in {
                        "extra_aligned_confirmation", "true_tail",
                    }), "unproven aligned/tail cells entered AUTO")

        path_names = {pair[0] for pair in PRODUCTION_DECODE_PATH_RULE_PAIRS}
        rule_names = {pair[1] for pair in PRODUCTION_DECODE_PATH_RULE_PAIRS}
        validate_selector_pair_contract(Path(__file__).resolve().parents[3])
        require(all(coherent_production_decode_pair(*pair)
                    for pair in PRODUCTION_DECODE_PATH_RULE_PAIRS),
                "a declared production path/rule pair is not coherent")
        for path_name in path_names | {"banana"}:
            for rule_name in rule_names | {"banana"}:
                require(coherent_production_decode_pair(path_name, rule_name) ==
                        ((path_name, rule_name) in
                         PRODUCTION_DECODE_PATH_RULE_PAIRS),
                        "unknown or cross-paired selector names were accepted")
        for invalid_pair in (
            ("banana", "banana"), ("materialized", "direct"),
            ("direct", "workspace_tiled"),
        ):
            require(not coherent_production_decode_pair(*invalid_pair),
                    f"explicit invalid selector pair was accepted: {invalid_pair!r}")

        materialized_case = next(
            item for item in attestation["cases"]
            if expected_auto_decode_pair(item, "scalar") ==
                ("materialized", "workspace_materialized") and
               item["expected_selected_decode_path"] == "not_generic")
        expected_materialized = expected_auto_decode_pair(
            materialized_case, "scalar")
        for pair in sorted(PRODUCTION_DECODE_PATH_RULE_PAIRS):
            candidate_output = fake_attestation_output(
                materialized_case, pair[0], pair[1], "scalar")
            accepted = True
            try:
                validate_attestation_output(candidate_output, materialized_case)
            except PlanError:
                accepted = False
            require(accepted == (pair == expected_materialized),
                    f"negative cell pair constraint differs for {pair!r}")
        for pair in (
            ("banana", "banana"), ("materialized", "direct"),
            ("direct", "workspace_tiled"),
        ):
            try:
                validate_attestation_output(fake_attestation_output(
                    materialized_case, pair[0], pair[1], "scalar"),
                    materialized_case)
            except PlanError:
                pass
            else:
                raise PlanError(f"negative cell accepted invalid pair {pair!r}")

        direct_case = attestation_case(
            "fixture-direct", "loss_rate_neighbor", 8, 8, 4096, 4,
            deterministic_seed("fixture-direct", 8, 4), False)
        tiled_case = attestation_case(
            "fixture-tiled", "loss_rate_neighbor", 127, 63, 4096, 63,
            deterministic_seed("fixture-tiled", 127, 63), False)
        require(expected_auto_decode_pair(direct_case, "scalar") ==
                ("direct", "direct") and
                expected_auto_decode_pair(tiled_case, "scalar") ==
                ("tiled", "workspace_tiled"),
                "terminal/workspace fixture predictions differ")
        low_workspace_case = {"cell": {
            "K": 33, "R": 7, "profile": "low_v1", "loss_count": 7,
        }}
        require(expected_required_work_slots(direct_case, "no_op") == 0 and
                expected_required_work_slots(direct_case, "direct") == 0 and
                expected_required_work_slots(tiled_case, "tiled") ==
                    2 * ceil_power_of_two(63) + 63 and
                expected_required_work_slots(tiled_case, "materialized") == 256 and
                expected_required_work_slots(low_workspace_case, "tiled") == 128 and
                expected_required_work_slots(low_workspace_case, "generic") == 128,
                "production workspace formulas differ")
        validate_attestation_output(fake_attestation_output(direct_case), direct_case)
        validate_attestation_output(fake_attestation_output(tiled_case), tiled_case)

        source = {
            "head": survivor_signed["candidate_commit"],
            "tree": survivor_signed["evidence_scope"]["sources"][
                "candidate"]["tree"],
            "status": "clean",
            "status_sha256": common.EMPTY_SHA256,
        }
        identity = {
            "relative_path": "fixture", "sha256": "5" * 64,
            "size": 1, "mode": 0o644,
        }
        artifact_closure = {
            "cache": identity, "graph": identity,
            "sources": {key: identity for key in (
                "benchmark", "decoder", "dispatch")},
            "objects": {key: identity for key in ("benchmark", "decoder")},
            "archive": identity, "binary": identity,
        }
        canonical_header_path = (
            "$BUILD/generated/leopard2-benchmark-attestation/"
            "leopard2_benchmark_source_attestation.h")
        canonical_header_text = _benchmark_attestation_text(
            source["head"], source["tree"])
        canonical_header_bytes = canonical_header_text.encode("utf-8")
        canonical_header_content = {
            "encoding": "utf-8",
            "size": len(canonical_header_bytes),
            "sha256": hashlib.sha256(canonical_header_bytes).hexdigest(),
            "text": canonical_header_text,
        }
        canonical_header_artifact = {
            "path": canonical_header_path,
            "kind": "generated_compile_input",
            "size": canonical_header_content["size"],
            "mode": 0o644,
            "sha256": canonical_header_content["sha256"],
        }
        canonical_configuration_entries = {
            variable: "" for variable in BUILD_CONFIGURATION_VARIABLES
        }
        canonical_configuration_entries.update(
            REQUIRED_BUILD_CONFIGURATION_ENTRIES)
        canonical_configuration_entries.update({
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_GENERATOR": "Unix Makefiles",
            "CMAKE_CONFIGURATION_TYPES": "",
            "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
        })
        canonical_configuration_material = _build_configuration_material(
            canonical_configuration_entries)
        canonical_configuration_digest = hashlib.sha256(
            canonical_configuration_material).hexdigest()
        canonical_configuration_text = (
            f"schema={BUILD_CONFIGURATION_FILE_SCHEMA}\n"
            f"sha256={canonical_configuration_digest}\n"
        ).encode("ascii") + canonical_configuration_material
        canonical_configuration_content = retained_text(
            canonical_configuration_text.decode("utf-8"))
        canonical_configuration = {
            "schema": BUILD_CONFIGURATION_RECORD_SCHEMA,
            "artifact": {
                "path": f"$BUILD/{BUILD_CONFIGURATION_RELATIVE_PATH}",
                "kind": "generated_build_configuration",
                "size": canonical_configuration_content["size"],
                "mode": 0o644,
                "sha256": canonical_configuration_content["sha256"],
            },
            "content": canonical_configuration_content,
            "configuration_schema": BUILD_CONFIGURATION_FILE_SCHEMA,
            "configuration_sha256": canonical_configuration_digest,
            "entries": canonical_configuration_entries,
            "embedded_build_type": "Release",
            "helper_source": {
                "path":
                    f"$SOURCE/{BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH}",
                "kind": "source_file", "size": 1, "mode": 0o644,
                "sha256": "7" * 64,
            },
        }
        canonical_benchmark_arguments = [
            "/usr/bin/compiler",
            "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256="
            f'"{canonical_configuration_digest}"',
            '-DLEO2_BENCHMARK_BUILD_TYPE="Release"',
            "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1",
            "-DLEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER="
            f'"{canonical_header_path}"',
            "-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=1",
            "-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=1",
            "-I$SOURCE",
            "-Wall", "-Wextra", "-fopenmp", "-O3", "-DNDEBUG", "-O3",
            "-std=gnu++11", "-fopenmp",
            "-o",
            "CMakeFiles/bench_leopard2.dir/bench/leopard2/benchmark.cpp.o",
            "-c", "$SOURCE/bench/leopard2/benchmark.cpp",
        ]
        canonical_provenance = {
            "validated_cache": {
                **REQUIRED_CANDIDATE_CACHE,
                "CMAKE_BUILD_TYPE": "Release",
                "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
                "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG -O3",
                "CMAKE_CONFIGURATION_TYPES": "",
                "CMAKE_GENERATOR": "Unix Makefiles",
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    BUILD_CONFIGURATION_FILE_SCHEMA,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
                    canonical_configuration_digest,
            },
            "validated_compile_commands": {
                "schema": COMPILE_COMMANDS_SCHEMA,
                "implementation": "candidate",
                "profile": CANDIDATE_COMPILE_PROFILE,
                "validated_optimization": "-O3",
                "validated_openmp": True,
                "required_entries": [{
                    "directory": "$BUILD",
                    "file": "$SOURCE/bench/leopard2/benchmark.cpp",
                    "output":
                        "CMakeFiles/bench_leopard2.dir/"
                        "bench/leopard2/benchmark.cpp.o",
                    "arguments": canonical_benchmark_arguments,
                }],
                "generated_attestation_header": {
                    "artifact": canonical_header_artifact,
                    "content": canonical_header_content,
                    "source_commit": source["head"],
                    "source_tree": source["tree"],
                    "source_tracked_dirty": False,
                },
                "effective_build_configuration": canonical_configuration,
            },
            "archive_link_recipe_content": {"sha256": "6" * 64},
        }
        build = {
            "artifact_closure": artifact_closure,
            "canonical_production": {
                "schema": CANONICAL_PRODUCTION_BUILD_SCHEMA,
                "validator": CANONICAL_BUILD_VALIDATOR,
                "provenance": canonical_provenance,
                "provenance_sha256": canonical_sha256(canonical_provenance),
            },
        }
        collector = {
            **identity, "relative_path":
                "experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py",
        }

        stage_signed = load_json(stage_root / "stage.json")
        stage_unsigned = unsigned(stage_signed, STAGE_SCHEMA, "fixture stage")
        forced_timing = []
        exact_timing = []
        for index, stage_artifact in enumerate(stage_unsigned["artifacts"]):
            kind = stage_artifact["kind"]
            if kind not in {
                    "forced_surviving_cells", "forced_non_aligned_tail",
                    "exact_main_aligned_confirmation",
                    "exact_main_rejection_timing"}:
                continue
            manifest_reference = {
                "path": str(
                    (root / "fixture-timing" /
                     f"{index:04d}-{stage_artifact['sha256']}.json").absolute()),
                "size": 1,
                "sha256": hashlib.sha256(
                    f"fixture-{index}".encode("ascii")).hexdigest(),
                "payload_digest": hashlib.sha256(
                    f"payload-{index}".encode("ascii")).hexdigest(),
            }
            record = {
                "stage_artifact": {
                    key: stage_artifact[key]
                    for key in ("path", "kind", "size", "sha256")
                },
                "manifest": manifest_reference,
            }
            if kind.startswith("forced_"):
                record["summary_content_sha256"] = hashlib.sha256(
                    f"summary-{index}".encode("ascii")).hexdigest()
                record["maximum_decode_ci95_upper"] = 1.0
                forced_timing.append(record)
            else:
                record["minimum_decode_ci95_lower"] = 1.0
                exact_timing.append(record)
        timing_evidence = signed({
            "schema": PROMOTION_TIMING_SCHEMA,
            "candidate_commit": stage_unsigned["candidate_commit"],
            "stage_content_sha256": stage_signed["content_sha256"],
            "neighbor_maximum_regression": NEIGHBOR_MAXIMUM_REGRESSION,
            "maximum_forced_survivor_decode_ci95_upper": 1.0,
            "minimum_exact_decode_ci95_lower": 1.0,
            "forced_results": sorted(
                forced_timing,
                key=lambda item: item["stage_artifact"]["path"]),
            "exact_main_results": sorted(
                exact_timing,
                key=lambda item: item["stage_artifact"]["path"]),
        })
        validate_promotion_timing_shape(stage_root, timing_evidence)

        def reject_timing_mutation(label: str, mutate) -> None:
            forged_timing = json.loads(json.dumps(timing_evidence))
            forged_timing.pop("content_sha256")
            mutate(forged_timing)
            forged_timing = signed(forged_timing)
            try:
                validate_promotion_timing_shape(stage_root, forged_timing)
            except PlanError:
                return
            raise PlanError(
                f"adversarial promotion timing {label} was accepted")

        reject_timing_mutation(
            "missing-forced",
            lambda value: value["forced_results"].pop())
        reject_timing_mutation(
            "missing-exact",
            lambda value: value["exact_main_results"].pop())
        reject_timing_mutation(
            "regression",
            lambda value: (
                value["exact_main_results"][0].update({
                    "minimum_decode_ci95_lower": 0.979}),
                value.update({"minimum_exact_decode_ci95_lower": 0.979})))
        reject_timing_mutation(
            "same-binary-regression",
            lambda value: (
                next(
                    item for item in value["forced_results"]
                    if item["stage_artifact"]["kind"] ==
                        "forced_surviving_cells").update({
                            "maximum_decode_ci95_upper": 1.021}),
                value.update({
                    "maximum_forced_survivor_decode_ci95_upper": 1.021})))
        reject_timing_mutation(
            "duplicate",
            lambda value: value["exact_main_results"].append(
                value["exact_main_results"][0]))
        reject_timing_mutation(
            "wrong-stage-artifact",
            lambda value: value["exact_main_results"][0][
                "stage_artifact"].update({"sha256": "0" * 64}))

        result_root = root / "attestation-result"
        raw_root = result_root / "raw"
        raw_root.mkdir(parents=True)
        raw_documents = {}
        raw_artifacts = {}
        for index, case in enumerate(attestation["cases"]):
            identifier = case["cell"]["identifier"]
            document = fake_attestation_output(case)
            path = raw_root / f"{index:04d}-{identifier}.json"
            write_json(path, document)
            raw_documents[identifier] = document
            raw_artifacts[identifier] = {
                "path": str(path.relative_to(result_root)),
                "size": path.stat().st_size, "sha256": file_sha256(path),
            }
        result = derive_attestation_result(
            stage_root, source, build, collector, raw_documents, raw_artifacts,
            timing_evidence)
        header_mutations = {}

        missing_header = json.loads(json.dumps(build))
        missing_header["canonical_production"]["provenance"][
            "validated_compile_commands"].pop(
                "generated_attestation_header")
        header_mutations["missing"] = missing_header

        for label, field, replacement in (
            ("commit", "source_commit", "a" * 40),
            ("tree", "source_tree", "b" * 40),
            ("dirty", "source_tracked_dirty", True),
        ):
            mutated = json.loads(json.dumps(build))
            mutated["canonical_production"]["provenance"][
                "validated_compile_commands"][
                    "generated_attestation_header"][field] = replacement
            header_mutations[label] = mutated

        content_drift = json.loads(json.dumps(build))
        content_record = content_drift["canonical_production"]["provenance"][
            "validated_compile_commands"]["generated_attestation_header"]
        drifted_text = content_record["content"]["text"].replace(
            "#endif\n", "#define LEO2_UNTRUSTED_SHADOW 1\n#endif\n")
        drifted_bytes = drifted_text.encode("utf-8")
        content_record["content"].update({
            "text": drifted_text,
            "size": len(drifted_bytes),
            "sha256": hashlib.sha256(drifted_bytes).hexdigest(),
        })
        content_record["artifact"].update({
            "size": len(drifted_bytes),
            "sha256": hashlib.sha256(drifted_bytes).hexdigest(),
        })
        header_mutations["content"] = content_drift

        hash_drift = json.loads(json.dumps(build))
        hash_drift["canonical_production"]["provenance"][
            "validated_compile_commands"]["generated_attestation_header"][
                "artifact"]["sha256"] = "0" * 64
        header_mutations["hash"] = hash_drift

        shadow_input = json.loads(json.dumps(build))
        shadow_arguments = shadow_input["canonical_production"]["provenance"][
            "validated_compile_commands"]["required_entries"][0]["arguments"]
        shadow_arguments[:] = [
            token.replace(
                canonical_header_path,
                "$SOURCE/bench/leopard2/"
                "leopard2_benchmark_source_attestation.h")
            for token in shadow_arguments]
        header_mutations["source-local-compile-input"] = shadow_input

        duplicate_override = json.loads(json.dumps(build))
        duplicate_arguments = duplicate_override["canonical_production"][
            "provenance"]["validated_compile_commands"][
                "required_entries"][0]["arguments"]
        duplicate_arguments.insert(
            duplicate_arguments.index("-I$SOURCE"),
            "-DLEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER="
            '"$SOURCE/bench/leopard2/'
            'leopard2_benchmark_source_attestation.h"')
        header_mutations["duplicate-compile-input-override"] = \
            duplicate_override

        missing_enable = json.loads(json.dumps(build))
        missing_enable_arguments = missing_enable["canonical_production"][
            "provenance"]["validated_compile_commands"][
                "required_entries"][0]["arguments"]
        missing_enable_arguments.remove(
            "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1")
        header_mutations["missing-enable-definition"] = missing_enable

        undef_override = json.loads(json.dumps(build))
        undef_arguments = undef_override["canonical_production"][
            "provenance"]["validated_compile_commands"][
                "required_entries"][0]["arguments"]
        undef_arguments.insert(
            undef_arguments.index("-I$SOURCE"),
            "-ULEO2_BENCHMARK_SOURCE_ATTESTATION")
        header_mutations["undef-enable-override"] = undef_override

        undef_header = json.loads(json.dumps(build))
        undef_header_arguments = undef_header["canonical_production"][
            "provenance"]["validated_compile_commands"][
                "required_entries"][0]["arguments"]
        undef_header_arguments.insert(
            undef_header_arguments.index("-I$SOURCE"),
            "-ULEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER")
        header_mutations["undef-header-override"] = undef_header

        for label, mutated in header_mutations.items():
            provenance = mutated["canonical_production"]["provenance"]
            mutated["canonical_production"]["provenance_sha256"] = \
                canonical_sha256(provenance)
            try:
                derive_attestation_result(
                    stage_root, source, mutated, collector,
                    raw_documents, raw_artifacts, timing_evidence)
            except PlanError:
                pass
            else:
                raise PlanError(
                    "attestation accepted generated-header " + label + " drift")

        noncanonical_build = json.loads(json.dumps(build))
        provenance = noncanonical_build["canonical_production"]["provenance"]
        provenance["validated_cache"]["LEO2_BACKEND_VARIANT"] = "scalar"
        noncanonical_build["canonical_production"][
            "provenance_sha256"] = canonical_sha256(provenance)
        try:
            derive_attestation_result(
                stage_root, source, noncanonical_build, collector,
                raw_documents, raw_artifacts, timing_evidence)
        except PlanError:
            pass
        else:
            raise PlanError("noncanonical scalar-only build attested as AUTO")
        manifest_path = result_root / "manifest.json"
        write_json(manifest_path, result)
        validate_attestation_result_files(
            stage_root, manifest_path, source, build, collector,
            verified_timing_evidence=timing_evidence)
        def historical_canonical_build_fixture(
            compile_schema: str,
        ) -> dict[str, Any]:
            require(compile_schema in {
                        COMPILE_COMMANDS_SCHEMA_V4,
                        COMPILE_COMMANDS_SCHEMA_V5,
                        COMPILE_COMMANDS_SCHEMA_V6,
                        COMPILE_COMMANDS_SCHEMA_V7,
                        COMPILE_COMMANDS_SCHEMA_V8,
                        COMPILE_COMMANDS_SCHEMA_V9,
                        COMPILE_COMMANDS_SCHEMA_V10,
                        COMPILE_COMMANDS_SCHEMA_V11},
                    "historical canonical fixture schema differs")
            if compile_schema == COMPILE_COMMANDS_SCHEMA_V11:
                file_schema = BUILD_CONFIGURATION_FILE_SCHEMA
                canonical_schema = CANONICAL_PRODUCTION_BUILD_SCHEMA_V9
                validator = CANONICAL_BUILD_VALIDATOR_V9
                profile = CANDIDATE_COMPILE_PROFILE_V7
            elif compile_schema == COMPILE_COMMANDS_SCHEMA_V10:
                file_schema = BUILD_CONFIGURATION_FILE_SCHEMA_V8
                canonical_schema = CANONICAL_PRODUCTION_BUILD_SCHEMA_V8
                validator = CANONICAL_BUILD_VALIDATOR_V8
                profile = CANDIDATE_COMPILE_PROFILE_V6
            elif compile_schema == COMPILE_COMMANDS_SCHEMA_V9:
                file_schema = BUILD_CONFIGURATION_FILE_SCHEMA_V7
                canonical_schema = CANONICAL_PRODUCTION_BUILD_SCHEMA_V7
                validator = CANONICAL_BUILD_VALIDATOR_V7
                profile = CANDIDATE_COMPILE_PROFILE_V5
            elif compile_schema == COMPILE_COMMANDS_SCHEMA_V8:
                file_schema = BUILD_CONFIGURATION_FILE_SCHEMA_V6
                canonical_schema = CANONICAL_PRODUCTION_BUILD_SCHEMA_V6
                validator = CANONICAL_BUILD_VALIDATOR_V6
                profile = CANDIDATE_COMPILE_PROFILE_V4
            elif compile_schema == COMPILE_COMMANDS_SCHEMA_V7:
                file_schema = BUILD_CONFIGURATION_FILE_SCHEMA_V5
                canonical_schema = CANONICAL_PRODUCTION_BUILD_SCHEMA_V5
                validator = CANONICAL_BUILD_VALIDATOR_V5
                profile = CANDIDATE_COMPILE_PROFILE_V3
            elif compile_schema == COMPILE_COMMANDS_SCHEMA_V6:
                file_schema = BUILD_CONFIGURATION_FILE_SCHEMA_V4
                canonical_schema = CANONICAL_PRODUCTION_BUILD_SCHEMA_V4
                validator = CANONICAL_BUILD_VALIDATOR_V4
                profile = CANDIDATE_COMPILE_PROFILE_V3
            elif compile_schema == COMPILE_COMMANDS_SCHEMA_V5:
                file_schema = BUILD_CONFIGURATION_FILE_SCHEMA_V3
                canonical_schema = CANONICAL_PRODUCTION_BUILD_SCHEMA_V3
                validator = CANONICAL_BUILD_VALIDATOR_V3
                profile = CANDIDATE_COMPILE_PROFILE_V2
            else:
                file_schema = BUILD_CONFIGURATION_FILE_SCHEMA_V2
                canonical_schema = CANONICAL_PRODUCTION_BUILD_SCHEMA_V2
                validator = CANONICAL_BUILD_VALIDATOR_V2
                profile = CANDIDATE_COMPILE_PROFILE_V2
            record_schema, variables, required_cache = \
                _build_configuration_contract_for_compile_schema(
                    compile_schema, file_schema)
            entries = {variable: "" for variable in variables}
            entries.update(_required_build_configuration_entries(
                file_schema, compile_schema))
            entries.update({
                "CMAKE_BUILD_TYPE": "Release",
                "CMAKE_GENERATOR": "Unix Makefiles",
                "CMAKE_CONFIGURATION_TYPES": "",
                "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
            })
            material = _build_configuration_material(entries, variables)
            digest = hashlib.sha256(material).hexdigest()
            configuration_text = (
                f"schema={file_schema}\n"
                f"sha256={digest}\n"
            ).encode("ascii") + material
            configuration_content = retained_text(
                configuration_text.decode("utf-8"))
            configuration = json.loads(json.dumps(canonical_configuration))
            configuration.update({
                "schema": record_schema,
                "content": configuration_content,
                "configuration_schema": file_schema,
                "configuration_sha256": digest,
                "entries": entries,
            })
            configuration["artifact"].update({
                "size": configuration_content["size"],
                "sha256": configuration_content["sha256"],
            })
            historical_build = json.loads(json.dumps(build))
            historical_canonical = historical_build["canonical_production"]
            historical_canonical.update({
                "schema": canonical_schema,
                "validator": validator,
            })
            provenance = historical_canonical["provenance"]
            provenance["validated_cache"] = {
                **required_cache,
                "CMAKE_BUILD_TYPE": "Release",
                "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
                "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG -O3",
                "CMAKE_CONFIGURATION_TYPES": "",
                "CMAKE_GENERATOR": "Unix Makefiles",
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA": file_schema,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": digest,
            }
            semantics = provenance["validated_compile_commands"]
            semantics.update({
                "schema": compile_schema,
                "profile": profile,
                "effective_build_configuration": configuration,
            })
            arguments = _normalized_compile_argv(
                "candidate",
                "$CANDIDATE_SOURCE/bench/leopard2/benchmark.cpp",
                "/usr/bin/compiler", compile_schema,
                source["head"], source["tree"], configuration)
            semantics["required_entries"][0]["arguments"] = [
                token.replace("$CANDIDATE_SOURCE", "$SOURCE").replace(
                    "$CANDIDATE_BUILD", "$BUILD")
                for token in arguments
            ]
            historical_canonical["provenance_sha256"] = canonical_sha256(
                provenance)
            return historical_build

        historical_build = historical_canonical_build_fixture(
            COMPILE_COMMANDS_SCHEMA_V4)
        historical_v9_build = historical_canonical_build_fixture(
            COMPILE_COMMANDS_SCHEMA_V5)
        historical_v10_build = historical_canonical_build_fixture(
            COMPILE_COMMANDS_SCHEMA_V6)
        historical_v11_build = historical_canonical_build_fixture(
            COMPILE_COMMANDS_SCHEMA_V7)
        historical_v12_build = historical_canonical_build_fixture(
            COMPILE_COMMANDS_SCHEMA_V8)
        historical_v13_build = historical_canonical_build_fixture(
            COMPILE_COMMANDS_SCHEMA_V9)
        historical_v14_build = historical_canonical_build_fixture(
            COMPILE_COMMANDS_SCHEMA_V10)
        historical_root = root / "historical-attestation-result"
        shutil.copytree(raw_root, historical_root / "raw")
        historical = derive_attestation_result(
            stage_root, source, historical_build, collector,
            raw_documents, raw_artifacts,
            None, result_schema=ATTESTATION_RESULT_SCHEMA_V5)
        historical_manifest = historical_root / "manifest.json"
        write_json(historical_manifest, historical)
        validated_historical = validate_attestation_result_files(
            stage_root, historical_manifest, source,
            historical_build, collector)
        require(validated_historical["schema"] ==
                ATTESTATION_RESULT_SCHEMA_V5,
                "historical attestation schema did not replay strictly")

        historical_timing_root = root / "historical-timing-attestation-result"
        shutil.copytree(raw_root, historical_timing_root / "raw")
        historical_timing = derive_attestation_result(
            stage_root, source, historical_build, collector,
            raw_documents, raw_artifacts, timing_evidence,
            result_schema=ATTESTATION_RESULT_SCHEMA_V6)
        historical_timing_manifest = historical_timing_root / "manifest.json"
        write_json(historical_timing_manifest, historical_timing)
        validated_historical_timing = validate_attestation_result_files(
            stage_root, historical_timing_manifest, source,
            historical_build, collector,
            verified_timing_evidence=timing_evidence)
        require(validated_historical_timing["schema"] ==
                ATTESTATION_RESULT_SCHEMA_V6,
                "historical timing attestation schema did not replay strictly")

        historical_v9_root = root / "historical-v9-attestation-result"
        shutil.copytree(raw_root, historical_v9_root / "raw")
        historical_v9 = derive_attestation_result(
            stage_root, source, historical_v9_build, collector,
            raw_documents, raw_artifacts, timing_evidence,
            result_schema=ATTESTATION_RESULT_SCHEMA_V7)
        historical_v9_manifest = historical_v9_root / "manifest.json"
        write_json(historical_v9_manifest, historical_v9)
        validated_historical_v9 = validate_attestation_result_files(
            stage_root, historical_v9_manifest, source,
            historical_v9_build, collector,
            verified_timing_evidence=timing_evidence)
        require(validated_historical_v9["schema"] ==
                ATTESTATION_RESULT_SCHEMA_V7,
                "historical v9 attestation schema did not replay strictly")

        historical_v10_root = root / "historical-v10-attestation-result"
        shutil.copytree(raw_root, historical_v10_root / "raw")
        historical_v10 = derive_attestation_result(
            stage_root, source, historical_v10_build, collector,
            raw_documents, raw_artifacts, timing_evidence,
            result_schema=ATTESTATION_RESULT_SCHEMA_V8)
        historical_v10_manifest = historical_v10_root / "manifest.json"
        write_json(historical_v10_manifest, historical_v10)
        validated_historical_v10 = validate_attestation_result_files(
            stage_root, historical_v10_manifest, source,
            historical_v10_build, collector,
            verified_timing_evidence=timing_evidence)
        require(validated_historical_v10["schema"] ==
                ATTESTATION_RESULT_SCHEMA_V8,
                "historical v10 attestation schema did not replay strictly")

        historical_v11_root = root / "historical-v11-attestation-result"
        shutil.copytree(raw_root, historical_v11_root / "raw")
        historical_v11 = derive_attestation_result(
            stage_root, source, historical_v11_build, collector,
            raw_documents, raw_artifacts, timing_evidence,
            result_schema=ATTESTATION_RESULT_SCHEMA_V9)
        historical_v11_manifest = historical_v11_root / "manifest.json"
        write_json(historical_v11_manifest, historical_v11)
        validated_historical_v11 = validate_attestation_result_files(
            stage_root, historical_v11_manifest, source,
            historical_v11_build, collector,
            verified_timing_evidence=timing_evidence)
        require(validated_historical_v11["schema"] ==
                ATTESTATION_RESULT_SCHEMA_V9,
                "historical v11 attestation schema did not replay strictly")

        historical_v12_root = root / "historical-v12-attestation-result"
        shutil.copytree(raw_root, historical_v12_root / "raw")
        historical_v12 = derive_attestation_result(
            stage_root, source, historical_v12_build, collector,
            raw_documents, raw_artifacts, timing_evidence,
            result_schema=ATTESTATION_RESULT_SCHEMA_V10)
        historical_v12_manifest = historical_v12_root / "manifest.json"
        write_json(historical_v12_manifest, historical_v12)
        validated_historical_v12 = validate_attestation_result_files(
            stage_root, historical_v12_manifest, source,
            historical_v12_build, collector,
            verified_timing_evidence=timing_evidence)
        require(validated_historical_v12["schema"] ==
                ATTESTATION_RESULT_SCHEMA_V10,
                "historical v12 attestation schema did not replay strictly")

        historical_v13_root = root / "historical-v13-attestation-result"
        shutil.copytree(raw_root, historical_v13_root / "raw")
        historical_v13 = derive_attestation_result(
            stage_root, source, historical_v13_build, collector,
            raw_documents, raw_artifacts, timing_evidence,
            result_schema=ATTESTATION_RESULT_SCHEMA_V11)
        historical_v13_manifest = historical_v13_root / "manifest.json"
        write_json(historical_v13_manifest, historical_v13)
        validated_historical_v13 = validate_attestation_result_files(
            stage_root, historical_v13_manifest, source,
            historical_v13_build, collector,
            verified_timing_evidence=timing_evidence)
        require(validated_historical_v13["schema"] ==
                ATTESTATION_RESULT_SCHEMA_V11,
                "historical v13 attestation schema did not replay strictly")

        historical_v14_root = root / "historical-v14-attestation-result"
        shutil.copytree(raw_root, historical_v14_root / "raw")
        historical_v14 = derive_attestation_result(
            stage_root, source, historical_v14_build, collector,
            raw_documents, raw_artifacts, timing_evidence,
            result_schema=ATTESTATION_RESULT_SCHEMA_V12)
        historical_v14_manifest = historical_v14_root / "manifest.json"
        write_json(historical_v14_manifest, historical_v14)
        validated_historical_v14 = validate_attestation_result_files(
            stage_root, historical_v14_manifest, source,
            historical_v14_build, collector,
            verified_timing_evidence=timing_evidence)
        require(validated_historical_v14["schema"] ==
                ATTESTATION_RESULT_SCHEMA_V12,
                "historical v14 attestation schema did not replay strictly")

        for label, result_schema, candidate_build, candidate_timing in (
            ("current body under historical result",
             ATTESTATION_RESULT_SCHEMA_V6, build, timing_evidence),
            ("historical body under current result",
             ATTESTATION_RESULT_SCHEMA, historical_build, timing_evidence),
            ("v15 body under v14 result",
             ATTESTATION_RESULT_SCHEMA_V12, build, timing_evidence),
            ("v14 body under v15 result",
             ATTESTATION_RESULT_SCHEMA, historical_v14_build, timing_evidence),
            ("v15 body under v13 result",
             ATTESTATION_RESULT_SCHEMA_V11, build, timing_evidence),
            ("v13 body under v15 result",
             ATTESTATION_RESULT_SCHEMA, historical_v13_build, timing_evidence),
            ("v14 body under v12 result",
             ATTESTATION_RESULT_SCHEMA_V10, build, timing_evidence),
            ("v12 body under v14 result",
             ATTESTATION_RESULT_SCHEMA, historical_v12_build, timing_evidence),
            ("v14 body under v11 result",
             ATTESTATION_RESULT_SCHEMA_V9, build, timing_evidence),
            ("v11 body under v14 result",
             ATTESTATION_RESULT_SCHEMA, historical_v11_build, timing_evidence),
            ("v14 body under v10 result",
             ATTESTATION_RESULT_SCHEMA_V8, build, timing_evidence),
            ("v10 body under v14 result",
             ATTESTATION_RESULT_SCHEMA, historical_v10_build, timing_evidence),
            ("v11 body under v9 result",
             ATTESTATION_RESULT_SCHEMA_V7, build, timing_evidence),
            ("v9 body under v11 result",
             ATTESTATION_RESULT_SCHEMA, historical_v9_build, timing_evidence),
        ):
            try:
                derive_attestation_result(
                    stage_root, source, candidate_build, collector,
                    raw_documents, raw_artifacts, candidate_timing,
                    result_schema=result_schema)
            except PlanError:
                pass
            else:
                raise PlanError(
                    "attestation accepted coherent schema relabel: " + label)

        for label, result_schema, candidate_build, candidate_timing, \
                extension_name in (
            (
                "historical cache-only current-selector extension",
                ATTESTATION_RESULT_SCHEMA_V6, historical_build,
                timing_evidence,
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
            ),
            (
                "current cache-only unknown extension",
                ATTESTATION_RESULT_SCHEMA, build, timing_evidence,
                "LEO2_UNVERSIONED_SELECTOR",
            ),
        ):
            extended_build = json.loads(json.dumps(candidate_build))
            extended_canonical = extended_build["canonical_production"]
            extended_provenance = extended_canonical["provenance"]
            extended_provenance["validated_cache"][extension_name] = "OFF"
            extended_canonical["provenance_sha256"] = canonical_sha256(
                extended_provenance)
            try:
                derive_attestation_result(
                    stage_root, source, extended_build, collector,
                    raw_documents, raw_artifacts, candidate_timing,
                    result_schema=result_schema)
            except PlanError:
                pass
            else:
                raise PlanError(
                    "attestation accepted build-cache extension: " + label)

        probe_case = attestation["cases"][0]
        probe = fake_attestation_output(probe_case)
        raw_mutations = {
            "K": lambda value: value["parameters"].update({
                "K": value["parameters"]["K"] + 1}),
            "seed": lambda value: value["parameters"].update({
                "seed": value["parameters"]["seed"] + 1}),
            "profile": lambda value: value["parameters"].update({
                "requested_profile": "low_v1"}),
            "field": lambda value: value["parameters"].update({
                "requested_field": "gf16"}),
            "requested-backend": lambda value: value["parameters"].update({
                "requested_backend": "scalar"}),
            "resolved-backend": lambda value: value["resolved"].update({
                "backend": "auto"}),
            "required-work-slots": lambda value: value["resolved"].update({
                "decode_required_work_slots":
                    value["resolved"]["decode_required_work_slots"] + 1}),
            "bytes": lambda value: value["parameters"].update({
                "shard_bytes": value["parameters"]["shard_bytes"] + 1}),
            "loss": lambda value: value["parameters"].update({
                "loss_count": (1 if value["parameters"]["loss_count"] != 1
                               else 2)}),
        }
        for label, mutate in raw_mutations.items():
            forged_output = json.loads(json.dumps(probe))
            mutate(forged_output)
            try:
                validate_attestation_output(forged_output, probe_case)
            except PlanError:
                continue
            raise PlanError(f"attestation accepted wrong {label}")

        numeric_digest_output = json.loads(json.dumps(probe))
        numeric_digest_output["workload_digests"]["original_data"] = \
            int("1" * 16)
        try:
            validate_attestation_output(numeric_digest_output, probe_case)
        except PlanError:
            pass
        else:
            raise PlanError(
                "attestation accepted a numeric 16-digit FNV digest")

        for excluded_backend in ("avx512", "neon"):
            forged_output = fake_attestation_output(
                probe_case, backend="scalar")
            forged_output["resolved"]["backend"] = excluded_backend
            try:
                validate_attestation_output(forged_output, probe_case)
            except PlanError:
                pass
            else:
                raise PlanError(
                    f"attestation accepted excluded backend {excluded_backend}")

        def reject_output_mutation(label: str, mutate) -> None:
            fixture = root / ("bad-" + label)
            shutil.copytree(result_root, fixture)
            mutate(fixture)
            try:
                validate_attestation_result_files(
                    stage_root, fixture / "manifest.json", source, build, collector,
                    verified_timing_evidence=timing_evidence)
            except PlanError:
                return
            raise PlanError(f"adversarial attestation {label} was accepted")

        first_record = result["records"][0]
        first_raw = Path(first_record["raw"]["path"])
        first_case = next(item for item in attestation["cases"]
                          if item["cell"]["identifier"] ==
                          first_record["identifier"])
        def wrong_raw(path: Path, key: str, value: object) -> None:
            document = load_json(path / first_raw)
            document["resolved"][key] = value
            target = path / first_raw
            target.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n",
                              encoding="utf-8")
            def reseal(value_doc):
                record = next(item for item in value_doc["records"]
                              if item["identifier"] == first_record["identifier"])
                record["raw"]["size"] = target.stat().st_size
                record["raw"]["sha256"] = file_sha256(target)
                record["selected_decode_path"] = document["resolved"][
                    "selected_decode_path"]
                record["selected_decode_rule"] = document["resolved"][
                    "selected_decode_rule"]
                record["benchmark_payload_sha256"] = canonical_sha256(document)
            adversarial_resign(path / "manifest.json", reseal)

        wrong_path = "materialized" if first_case[
            "expected_selected_decode_path"] == "generic" else "generic"
        reject_output_mutation("wrong-path", lambda path:
                               wrong_raw(path, "selected_decode_path", wrong_path))
        wrong_rule = "workspace_materialized" if first_case[
            "expected_selected_decode_rule"] == "balanced_generic" else \
            "balanced_generic"
        reject_output_mutation("wrong-rule", lambda path:
                               wrong_raw(path, "selected_decode_rule", wrong_rule))
        reject_output_mutation("missing", lambda path:
                               (path / first_raw).unlink())
        def symbolic_raw(path: Path) -> None:
            outside = root / "outside-raw.json"
            if not outside.exists():
                shutil.copy2(path / first_raw, outside)
            (path / first_raw).unlink()
            (path / first_raw).symlink_to(outside)
        reject_output_mutation("symbolic-raw", symbolic_raw)
        def duplicate(path: Path) -> None:
            def mutate(value):
                value["records"].append(value["records"][0])
            adversarial_resign(path / "manifest.json", mutate)
        reject_output_mutation("duplicate", duplicate)
        def missing_timing(path: Path) -> None:
            adversarial_resign(
                path / "manifest.json",
                lambda value: value.pop("promotion_timing_evidence"))
        reject_output_mutation("missing-timing", missing_timing)
        reject_output_mutation("extra", lambda path:
                               (path / "extra.json").write_text("{}\n", encoding="utf-8"))
        def wrong_commit(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value.update({"candidate_commit": "6" * 40}))
        reject_output_mutation("commit", wrong_commit)
        def wrong_source_identity(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value["source_identity"].update({
                                   "head": "8" * 40,
                               }))
        reject_output_mutation("source-identity", wrong_source_identity)
        def wrong_benchmark_identity(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value["build_identity"]["artifact_closure"][
                                   "binary"].update({
                                   "sha256": "9" * 64,
                               }))
        reject_output_mutation("benchmark-identity", wrong_benchmark_identity)
        def wrong_stage(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value.update({"stage_content_sha256": "7" * 64}))
        reject_output_mutation("stage", wrong_stage)
        def resign_record(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value["records"][0].update({
                                   "selected_decode_path": "generic",
                               }))
        reject_output_mutation("resigned-output", resign_record)

        broadened = root / "broadened-stage"
        shutil.copytree(stage_root, broadened)
        spec_path = broadened / "path-attestation.json"
        def broaden(value):
            target = next(item for item in value["cases"]
                          if item["expected_selected_decode_path"] == "not_generic")
            target["expected_selected_decode_path"] = "generic"
            target["expected_selected_decode_rule"] = "balanced_generic"
        adversarial_resign(spec_path, broaden)
        stage_path = broadened / "stage.json"
        def reseal_stage(value):
            target = next(item for item in value["artifacts"]
                          if item["path"] == "path-attestation.json")
            target["size"] = spec_path.stat().st_size
            target["sha256"] = file_sha256(spec_path)
        adversarial_resign(stage_path, reseal_stage)
        try:
            validate_stage(first, broadened, manifests, scopes)
        except PlanError:
            pass
        else:
            raise PlanError("re-signed broadened AUTO dispatch stage was accepted")
    print("balanced promotion plan self-test passed: exact-cell AUTO attestation + adversarial results")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    generate_parser = commands.add_parser("generate")
    generate_parser.add_argument("--output", required=True, type=Path)
    validate_parser = commands.add_parser("validate")
    validate_parser.add_argument("--input", required=True, type=Path)
    validate_parser.add_argument("--plan", type=Path)
    select_parser = commands.add_parser("select")
    select_parser.add_argument("--plan", required=True, type=Path)
    select_parser.add_argument("--gate-manifest", action="append", required=True,
                               type=Path)
    select_parser.add_argument("--output", required=True, type=Path)
    advance_parser = commands.add_parser("advance")
    advance_parser.add_argument("--plan", required=True, type=Path)
    advance_parser.add_argument("--survivors", required=True, type=Path)
    advance_parser.add_argument("--output", required=True, type=Path)
    collect_parser = commands.add_parser("collect-attestation")
    collect_parser.add_argument("--plan", required=True, type=Path)
    collect_parser.add_argument("--stage", required=True, type=Path)
    collect_parser.add_argument("--source-root", required=True, type=Path)
    collect_parser.add_argument("--binary", default="build-release/bench_leopard2")
    collect_parser.add_argument(
        "--forced-manifest", action="append", required=True, type=Path,
        help="complete same-binary forced-run manifest; repeat for every "
             "forced stage artifact")
    collect_parser.add_argument(
        "--exact-manifest", action="append", required=True, type=Path,
        help="complete exact-main confirmation/rejection manifest; repeat "
             "for every exact-main stage artifact")
    collect_parser.add_argument("--output", required=True, type=Path)
    verify_parser = commands.add_parser("verify-attestation")
    verify_parser.add_argument("--plan", required=True, type=Path)
    verify_parser.add_argument("--stage", required=True, type=Path)
    verify_parser.add_argument("--source-root", required=True, type=Path)
    verify_parser.add_argument("--binary", default="build-release/bench_leopard2")
    verify_parser.add_argument("--manifest", required=True, type=Path)
    commands.add_parser("self-test")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "generate":
        plan = generate_plan(args.output.absolute())
        result = {"content_sha256": plan["content_sha256"],
                  "initial_case_count": plan["initial_case_count"],
                  "initial_timed_child_count": plan["initial_timed_child_count"]}
    elif args.command == "validate":
        if args.plan is None:
            plan = validate_plan(args.input.absolute())
            result = {"content_sha256": (load_json(
                args.input.absolute() / "plan.json"))["content_sha256"],
                "initial_case_count": plan["initial_case_count"]}
        else:
            stage = validate_stage(args.plan.absolute(), args.input.absolute())
            result = {"content_sha256": (load_json(
                args.input.absolute() / "stage.json"))["content_sha256"],
                "timed_case_count": stage["timed_case_count"]}
    elif args.command == "select":
        survivor = select_survivors(args.plan.absolute(), args.gate_manifest,
                                    args.output.absolute())
        result = {"content_sha256": survivor["content_sha256"],
                  "survivor_count": len(survivor["survivor_cells"]),
                  "refinement_count": len(survivor["required_refinement_cells"])}
    elif args.command == "advance":
        survivor_path = args.survivors.resolve(strict=True)
        validate_survivors(args.plan.absolute(), survivor_path)
        stage = materialize_stage(args.plan.absolute(), load_json(survivor_path),
                                  args.output.absolute())
        result = {"content_sha256": stage["content_sha256"],
                  "timed_case_count": stage["timed_case_count"],
                  "timed_child_count": stage["timed_child_count"]}
    elif args.command == "collect-attestation":
        attestation = collect_attestation(
            args.plan.absolute(), args.stage.absolute(),
            args.source_root.absolute(), args.binary, args.output.absolute(),
            args.forced_manifest, args.exact_manifest)
        result = {
            "content_sha256": attestation["content_sha256"],
            "record_count": len(attestation["records"]),
            "resolved_backend": attestation["resolved_backend"],
        }
    elif args.command == "verify-attestation":
        attestation = validate_attestation_result(
            args.plan.absolute(), args.stage.absolute(), args.manifest.absolute(),
            args.source_root.absolute(), args.binary)
        result = {
            "content_sha256": (load_json(args.manifest.absolute()))[
                "content_sha256"],
            "record_count": len(attestation["records"]),
            "resolved_backend": attestation["resolved_backend"],
            "promotion_timing_complete":
                attestation["schema"] == ATTESTATION_RESULT_SCHEMA,
        }
    else:
        self_test()
        return 0
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, json.JSONDecodeError, PlanError,
            common.EvidenceError) as error:
        raise SystemExit(f"balanced promotion plan failed: {error}") from error
