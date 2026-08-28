#!/usr/bin/env python3

"""Independent retained-only replay of the v17 GFNI exact-main campaign.

The auditor intentionally imports neither ``run_abba`` nor any of its policy
helpers.  It never executes a benchmark.  It reopens the signed manifest, raw
bundle, and retained child streams, then independently reconstructs the v17
workload, production route, CPU-isolation result, timer floors, and clustered
ABBA inference.
"""

from __future__ import annotations

import argparse
import base64
import binascii
import datetime as dt
import hashlib
import json
import math
import os
import re
import shlex
import stat
import statistics
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence


AUDIT_SCHEMA = "leopard2-main-compare-v17-independent-audit/v1"
MANIFEST_SCHEMA = "leopard2-main-compare-manifest/v17"
RAW_SCHEMA = "leopard2-main-compare-raw/v17"
BASELINE_SCHEMA = "leopard-main-benchmark-v1"
CANDIDATE_SCHEMA = "leopard2-benchmark-v35"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
BASELINE_COMPILE_PROFILE = \
    "gnu-compatible-cxx11-native-x86_64-release/v1"
BASELINE_ISA_POLICY = "whole-build -march=native"
CANDIDATE_COMPILE_PROFILE = \
    "gnu-compatible-cxx11-runtime-dispatch-production-auto-gfni-encode-x86_64-release/v9"
CANDIDATE_ISA_POLICY = (
    "portable core with ISA flags only on SSSE3, AVX2, isolated "
    "generated and fixed-shape AVX2, and GFNI translation units; "
    "AVX-512 probes disabled; AUTO context resolved AVX2 with the "
    "bounded production AUTO GF16 GFNI encode route"
)
SEALED_EXECUTABLE_PROTOCOL = "linux-sealed-executable-memfd/v1"
SEALED_EXECUTABLE_COMMAND = {
    "baseline": "<sealed-baseline-executable>",
    "candidate": "<sealed-candidate-executable>",
}
REQUIRED_EXECUTABLE_SEALS = 0x0F
MAX_LINUX_SEAL_MASK = 0x7FFF_FFFF
NONEMPTY_ARTIFACT_KINDS = {
    "archive", "archiver", "build_metadata", "build_tool", "compiler",
    "dynamic_loader", "executable", "file", "generated_build_configuration",
    "generated_compile_input", "object_file", "ranlib", "shared_library",
    "source_file",
}
EXECUTABLE_ARTIFACT_KINDS = {
    "archiver", "build_tool", "compiler", "executable", "ranlib",
}
PAIR_LEASE_SCHEMA = "leopard2-cpu-pair-lease/v1"
RESERVATION_SCHEMA = "leopard2-cpu-reservation/v1"
ISOLATION_SCHEMA = "leopard2-main-compare-isolation/v1"
SUPERVISION_SCHEMA = "leopard2-main-supervision/v1"
GIT_CAPTURE_SCHEMA = "leopard2-git-source-capture/v2"
GIT_EXECUTABLE_PROTOCOL = "linux-sealed-git-executable-memfd/v1"
GIT_METADATA_GUARD_POLICY = \
    "retained-gitdir-commondir-recursive-inotify/v1"
GIT_WORKTREE_GUARD_POLICY = \
    "retained-root-and-tracked-path-components-inotify/v1"
CANONICAL_LDD_SCHEMA = "leopard2-main-compare-canonical-ldd/v1"
CANONICAL_LDD_NORMALIZATION = "terminal-aslr-load-address/v1"
CANONICAL_LDD_ADDRESS = "<ASLR_LOAD_ADDRESS>"
COMPILE_COMMANDS_SCHEMA = "leopard2-main-compare-compile-commands/v13"
BUILD_CONFIGURATION_RECORD_SCHEMA = \
    "leopard2-main-compare-build-configuration/v10"
BUILD_CONFIGURATION_FILE_SCHEMA = \
    "leopard2-benchmark-build-configuration/v13"
NINJA_GRAPH_CLOSURE_SCHEMA = "leopard2-ninja-graph-closure/v1"
CANONICAL_NINJA_PATH = "/usr/bin/ninja"
CANONICAL_TASKSET_PATH = "/usr/bin/taskset"
CANONICAL_LDD_PATH = "/usr/bin/ldd"
CANONICAL_GIT_PATH = "/usr/bin/git"
CANONICAL_ARCHIVER_PATH = "/usr/bin/ar"
CANONICAL_RANLIB_PATH = "/usr/bin/ranlib"
CANONICAL_ARCHIVER_RESOLVED_PATH = "/usr/bin/x86_64-linux-gnu-ar"
CANONICAL_RANLIB_RESOLVED_PATH = "/usr/bin/x86_64-linux-gnu-ranlib"
BUILD_CONFIGURATION_RELATIVE_PATH = (
    "generated/leopard2-benchmark-attestation/"
    "leopard2_benchmark_build_configuration.txt"
)
BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH = \
    "cmake/Leopard2BenchmarkAttestation.cmake"
EVIDENCE_HELPER_RELATIVE_PATH = \
    "experiments/leopard2/decoder_dispatch/balanced_evidence_common.py"
ROUNDS = 3
ORDER = ("baseline", "candidate", "candidate", "baseline")
ITERATIONS = 9
WARMUP = 2
REUSE = 8
CPU = 52
SIBLING = 116
TIMER_FLOOR_US = 20_000.0
T_CRITICAL_DF2 = 4.302652729911275
CELL = {
    "identifier": "gf16-high-full",
    "k": 1000,
    "r": 200,
    "shard_bytes": 65536,
    "losses": 200,
    "seed": 131,
}
ROUTE_CONTRACT = (
    "LEGACY_HIGH_V1,GF16,AUTO,resolved_AVX2,"
    "K=1000,R=200,T=256,B=65536,native_layout,"
    "native_cantor_affine,dense_full_parity,encode_only,"
    "runtime_GFNI,startup_KAT,calibrated_AMD_1A_08,"
    "one_shot_and_one_item_batch,codec_setup_descriptive_only"
)
ORDINARY_ENCODE_API = "leo2_encode_batch:item_count=1:no_preflight_scratch"
ONE_SHOT_ENCODE_API = "leo2_encode"
CHILD_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_THREAD_LIMIT": "1",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}
STATISTICS_POLICY = {
    "method": "one mean log contrast per independent ABBA round",
    "confidence": 0.95,
    "critical_distribution": "Student-t",
    "independent_round_count_per_metric": 3,
    "constituent_pair_count_per_metric": 6,
    "degrees_of_freedom": 2,
    "child_estimator": "median of retained per-invocation samples",
    "ordinary_encode_semantics": (
        "exact Leopard1 leo_encode divided by Leopard2 ordinary "
        "leo2_encode_batch(item_count=1,no preflight scratch); "
        "codec setup excluded"
    ),
    "one_shot_encode_semantics": (
        "the same exact Leopard1 leo_encode observations divided by "
        "Leopard2 leo2_encode; codec setup excluded"
    ),
    "ratio_orientation": "baseline_time_over_candidate_time",
    "ratios_are_separate_correlated_and_must_not_be_multiplied": True,
    "same_binary_gfni_over_avx2_is_a_separate_campaign": True,
    "retained_timer_window_floor_us": 20_000.0,
}
BASELINE_LIBRARY_SOURCES = (
    "leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp",
)
CANDIDATE_LIBRARY_SOURCES = (
    "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp",
    "LeopardFF16.cpp", "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp", "Leopard2BackendAVX2Xor.cpp",
    "Leopard2BackendAVX2T16Q2.cpp", "Leopard2BackendAVX2T8K62.cpp",
    "Leopard2BackendAVX2T16K66.cpp", "Leopard2BackendAVX2T32B256.cpp",
    "Leopard2BackendAVX2T16B64.cpp", "Leopard2LowP32B64AVX2.cpp",
    "Leopard2BackendAVX2T2K4.cpp", "Leopard2BackendAVX2T8K8B1024.cpp",
    "Leopard2BackendGFNI.cpp",
)
CANDIDATE_NON_LIBRARY_ACTIONS = (
    ("tests/benchmark.cpp", "bench_leopard.dir"),
    ("bench/leopard2/benchmark.cpp", "bench_leopard2.dir"),
    ("bench/leopard2/benchmark.cpp", "bench_leopard2_prevalidated_batch.dir"),
    ("bench/leopard2/batch_preflight_benchmark.cpp",
     "bench_leopard2_batch_preflight.dir"),
    ("bench/leopard2/locator_benchmark.cpp", "bench_leopard2_locator.dir"),
    ("bench/leopard2/sparse_encode_benchmark.cpp",
     "leopard2_sparse_encode_benchmark_object.dir"),
    ("tests/leopard2/direct_oracle.cpp",
     "leopard2_sparse_encode_oracle_object.dir"),
    ("experiments/leopard2/gf16_high_encode/attribution.cpp",
     "bench_leopard2_gf16_high_attribution.dir"),
    ("tests/experiments.cpp", "experiment_leopard.dir"),
)
CANDIDATE_TARGETS = {
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
    "Leopard2BackendAVX2T16Q2.cpp": "leopard2_backend_avx2.dir",
    "Leopard2BackendAVX2T8K62.cpp": "leopard2_backend_avx2.dir",
    "Leopard2BackendAVX2T16K66.cpp": "leopard2_backend_avx2_t16_k66.dir",
    "Leopard2BackendGFNI.cpp": "leopard2_backend_gfni.dir",
}
CANDIDATE_ISA_FLAGS = {
    "Leopard2BackendSSSE3.cpp": ("-mssse3", "-mno-avx"),
    "Leopard2BackendAVX2.cpp":
        ("-mavx2", "-mno-avx512f", "-falign-functions=64"),
    "Leopard2BackendAVX2Xor.cpp":
        ("-mavx2", "-mno-avx512f", "-falign-functions=64"),
    "Leopard2BackendAVX2T32B256.cpp":
        ("-mavx2", "-mno-avx512f", "-falign-functions=64"),
    "Leopard2BackendAVX2T16B64.cpp":
        ("-mavx2", "-mno-avx512f", "-falign-functions=64"),
    "Leopard2LowP32B64AVX2.cpp":
        ("-mavx2", "-mno-avx512f", "-falign-functions=64"),
    "Leopard2BackendAVX2T2K4.cpp":
        ("-mavx2", "-mno-avx512f", "-falign-functions=64"),
    "Leopard2BackendAVX2T8K8B1024.cpp": (
        "-mavx2", "-mno-avx512f", "-flive-range-shrinkage",
        "-falign-functions=64"),
    "Leopard2BackendAVX2T16Q2.cpp":
        ("-mavx2", "-mno-avx512f", "-falign-functions=64"),
    "Leopard2BackendAVX2T8K62.cpp":
        ("-mavx2", "-mno-avx512f", "-falign-functions=64"),
    "Leopard2BackendAVX2T16K66.cpp":
        ("-mavx2", "-mno-avx512f", "-falign-functions=64"),
    "Leopard2BackendGFNI.cpp":
        ("-mavx2", "-mgfni", "-mno-avx512f", "-falign-functions=64"),
}
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
BUILD_CONFIGURATION_VARIABLES_V9 = (
    *BUILD_CONFIGURATION_VARIABLES_V6,
    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK",
)
BUILD_CONFIGURATION_VARIABLES = (
    *BUILD_CONFIGURATION_VARIABLES_V9[:17],
    "LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED",
    *BUILD_CONFIGURATION_VARIABLES_V9[17:],
)
MANIFEST_KEYS = {
    "schema", "created_utc", "valid", "validity_is_independent_of_speed",
    "raw", "campaign", "host", "isolation", "reservation", "supervision",
    "identities", "executable_snapshots", "analysis", "digest",
}
RAW_KEYS = {
    "schema", "created_utc", "validity_is_independent_of_speed", "campaign",
    "host_initial", "isolation", "reservation", "supervision",
    "input_specification", "identities_initial", "executable_snapshots",
    "invocations", "identities_final", "host_final", "analysis", "digest",
}
CAMPAIGN_KEYS = {
    "rounds", "order", "cells", "candidate_mode", "batch", "reuse",
    "iterations", "warmup", "threads", "child_environment", "benchmark_cpu",
    "reserved_sibling", "timeout_seconds", "statistics",
    "allowed_cpu_set_at_launch",
}
INVOCATION_KEYS = {
    "cell_id", "round", "slot", "implementation", "command",
    "execution_protocol", "executable_snapshot", "environment", "pinned_cpu",
    "started_utc", "duration_ns", "returncode", "stdout", "stderr", "result",
    "normalized", "identity_before", "identity_after", "reservation_before",
    "reservation_after",
}
INPUT_KEYS = {
    "runner", "taskset", "ldd", "baseline_executable", "candidate_executable",
    "baseline_archive", "candidate_archive", "baseline_build_dir",
    "candidate_build_dir", "baseline_source_root", "candidate_source_root",
    "candidate_commit", "baseline_pure_avx2",
}
INPUT_SNAPSHOT_KEYS = {
    "runner", "taskset", "ldd", "evidence_helper",
    "baseline_executable", "candidate_executable",
    "baseline_archive", "candidate_archive",
    "baseline_build", "candidate_build",
    "baseline_runtime_closure", "candidate_runtime_closure",
    "baseline_source", "candidate_source",
}
INPUT_PATH_KEYS = INPUT_KEYS - {"candidate_commit", "baseline_pure_avx2"}
BUILD_KEYS = {
    "build_dir", "cmake_cache", "compile_commands",
    "executable_link_recipe", "archive_link_recipe", "compiler",
    "compiler_version_stdout", "archiver", "ranlib",
    "validated_archive_members", "validated_executable",
    "validated_archive", "validated_cache", "validated_compile_commands",
    "archive_link_recipe_content", "executable_link_recipe_content",
    "archive_link_tool_invocations", "compiler_invocation",
    "validated_external_link_inputs", "multi_config_build_tool",
    "multi_config_build_tool_version_stdout", "multi_config_ninja_graph",
}
COMPILE_RECORD_KEYS = {
    "entry_count", "required_sources", "validated_optimization",
    "validated_openmp", "required_source_object_pairs", "isa_policy",
    "schema", "implementation", "profile", "required_entries",
    "generated_attestation_header", "effective_build_configuration",
}
ARTIFACT_KEYS = {"path", "kind", "size", "mode", "mtime_ns", "sha256"}
STREAM_KEYS = {"path", "size", "sha256"}
MANIFEST_RAW_KEYS = {"path", "size", "sha256", "payload_digest"}
CPU_STAT_FIELDS = (
    "user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal",
)
CPU_IDLE_FIELDS = {"idle", "iowait"}
HEX40 = re.compile(r"^[0-9a-f]{40}$")
HEX64 = re.compile(r"^[0-9a-f]{16}$")
HEX256 = re.compile(r"^[0-9a-f]{64}$")
GNU_CXX = re.compile(
    r"(?:g\+\+|(?:[A-Za-z0-9_.+]+-)+(?:gnu|gnueabi(?:hf)?|eabi|elf|musl|"
    r"mingw32)-g\+\+)(?:-[0-9]+(?:\.[0-9]+)*)?")
UTC = re.compile(
    r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{6})?Z$")
MAX_MANIFEST_BYTES = 256 * 1024 * 1024
MAX_RAW_BYTES = 256 * 1024 * 1024
MAX_STDOUT_BYTES = 8 * 1024 * 1024
MAX_STDERR_BYTES = 1024 * 1024
MAX_CPU_ID = 1_048_575
MAX_CPU_LIST_ENTRIES = 4096
MAX_GIT_TREE_OBJECTS = 65_536
MAX_GIT_TREE_TOTAL_BYTES = 64 * 1024 * 1024
MAX_TRACKED_SOURCE_FILES = 16_384
MAX_TRACKED_SOURCE_FILE_BYTES = 64 * 1024 * 1024
MAX_SUBMODULE_DEPTH = 16
MAX_NINJA_GRAPH_FILES = 64
MAX_NINJA_GRAPH_TOTAL_BYTES = 8 * 1024 * 1024
MASK64 = (1 << 64) - 1


class AuditError(RuntimeError):
    """Retained v17 evidence failed one independent audit gate."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def exact_keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    require(type(value) is dict and set(value) == expected,
            f"{label} keys changed")
    return value


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def reject_constant(value: str) -> Any:
    raise AuditError(f"non-finite JSON constant: {value}")


def finite_float_token(value: str) -> float:
    result = float(value)
    require(math.isfinite(result), f"non-finite JSON float: {value}")
    return result


def require_finite_tree(value: Any, label: str) -> None:
    if type(value) is float:
        require(math.isfinite(value), f"{label} contains a non-finite number")
    elif type(value) is list:
        for item in value:
            require_finite_tree(item, label)
    elif type(value) is dict:
        for item in value.values():
            require_finite_tree(item, label)
    else:
        require(type(value) in (type(None), bool, int, str),
                f"{label} contains an unsupported JSON type")


def strict_json_bytes(data: bytes, label: str) -> Any:
    require(data, f"{label} is empty")
    try:
        text = data.decode("utf-8", errors="strict")
        result = json.loads(
            text, object_pairs_hook=strict_object,
            parse_constant=reject_constant, parse_float=finite_float_token)
        require_finite_tree(result, label)
        return result
    except AuditError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError,
            RecursionError, OverflowError) as error:
        raise AuditError(f"{label} JSON is invalid: {error}") from error


def canonical_bytes(value: Any) -> bytes:
    try:
        rendered = json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False,
            ensure_ascii=False).encode("utf-8", errors="strict")
    except (TypeError, ValueError, UnicodeEncodeError) as error:
        raise AuditError(f"value is not canonical JSON: {error}") from error
    return rendered.encode("utf-8")


def git_canonical_bytes(value: Any) -> bytes:
    """Mirror the Git-capture v2 identity spelling, including Unicode."""
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
            allow_nan=False).encode("utf-8", errors="strict")
    except (TypeError, ValueError, UnicodeEncodeError, RecursionError) as error:
        raise AuditError(f"Git capture identity is not canonical JSON: {error}") \
            from error


def exact_json_equal(left: Any, right: Any) -> bool:
    return canonical_bytes(left) == canonical_bytes(right)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha1_object(kind: str, content: bytes) -> str:
    header = f"{kind} {len(content)}\0".encode("ascii")
    return hashlib.sha1(header + content, usedforsecurity=False).hexdigest()


def is_canonical_absolute_path(value: Any) -> bool:
    return (type(value) is str and bool(value) and
            not any(character in value for character in "\0\r\n") and
            Path(value).is_absolute() and os.path.normpath(value) == value)


def validate_artifact(value: Any, label: str,
                      expected_kind: str | None = None) -> dict[str, Any]:
    value = exact_keys(value, ARTIFACT_KEYS, label)
    require(is_canonical_absolute_path(value["path"]) and
            type(value["kind"]) is str and value["kind"] and
            (expected_kind is None or value["kind"] == expected_kind) and
            type(value["size"]) is int and value["size"] >= 0 and
            type(value["mode"]) is int and 0 <= value["mode"] <= 0o7777 and
            type(value["mtime_ns"]) is int and value["mtime_ns"] >= 0 and
            type(value["sha256"]) is str and
            HEX256.fullmatch(value["sha256"]) is not None,
            f"{label} file identity is invalid")
    require(value["kind"] not in NONEMPTY_ARTIFACT_KINDS or value["size"] > 0,
            f"{label} artifact is unexpectedly empty")
    if value["kind"] in EXECUTABLE_ARTIFACT_KINDS:
        require(bool(value["mode"] & 0o111), f"{label} is not executable")
    return value


def text_identity(text: str) -> dict[str, Any]:
    try:
        data = text.encode("utf-8", errors="strict")
    except UnicodeEncodeError as error:
        raise AuditError("retained text is not strict UTF-8") from error
    return {
        "encoding": "utf-8", "size": len(data),
        "sha256": sha256_bytes(data), "text": text,
    }


def validate_text(value: Any, label: str) -> dict[str, Any]:
    value = exact_keys(value, {"encoding", "size", "sha256", "text"}, label)
    require(value["encoding"] == "utf-8" and type(value["text"]) is str and
            exact_json_equal(value, text_identity(value["text"])),
            f"{label} retained text identity differs")
    return value


def validate_version_text(value: Any, label: str) -> dict[str, Any]:
    value = exact_keys(value, {"sha256", "text"}, label)
    require(type(value["text"]) is str and value["text"] and
            "\0" not in value["text"] and "\r" not in value["text"],
            f"{label} text is invalid")
    data = value["text"].encode("utf-8", errors="strict")
    require(len(data) <= 64 * 1024 and value["sha256"] == sha256_bytes(data),
            f"{label} digest differs")
    return value


def parse_cpu_list(value: str, label: str) -> set[int]:
    require(type(value) is str and value and len(value.encode("ascii")) <= 65_536,
            f"{label} is not a bounded ASCII CPU list")
    text = value.strip()
    require(text and all(character in "0123456789,-" for character in text),
            f"{label} contains invalid characters")
    result: set[int] = set()
    for component in text.split(","):
        require(component and component.count("-") <= 1,
                f"{label} contains an invalid component")
        if "-" in component:
            first_text, last_text = component.split("-", 1)
            require(first_text.isdigit() and last_text.isdigit(),
                    f"{label} contains an invalid range")
            first, last = int(first_text), int(last_text)
            require(0 <= first <= last <= MAX_CPU_ID and
                    last - first + 1 <= MAX_CPU_LIST_ENTRIES,
                    f"{label} range is out of bounds")
            result.update(range(first, last + 1))
        else:
            require(component.isdigit(), f"{label} contains an invalid CPU")
            cpu = int(component)
            require(0 <= cpu <= MAX_CPU_ID, f"{label} CPU is out of bounds")
            result.add(cpu)
        require(len(result) <= MAX_CPU_LIST_ENTRIES,
                f"{label} contains too many CPUs")
    return result


def validate_byte_identity(value: Any, label: str) -> dict[str, Any]:
    value = exact_keys(value, {"size", "sha256"}, label)
    require(type(value["size"]) is int and value["size"] >= 0 and
            type(value["sha256"]) is str and
            HEX256.fullmatch(value["sha256"]) is not None,
            f"{label} byte identity is invalid")
    return value


def decode_git_object(value: Any, kind: str, expected_id: str,
                      label: str) -> bytes:
    value = exact_keys(value, {
        "encoding", "size", "sha256", "object_id", "base64"}, label)
    require(value["encoding"] == "base64" and
            value["object_id"] == expected_id and
            type(value["base64"]) is str,
            f"{label} object identity differs")
    try:
        content = base64.b64decode(value["base64"], validate=True)
    except (ValueError, binascii.Error) as error:
        raise AuditError(f"{label} is not canonical base64") from error
    require(base64.b64encode(content).decode("ascii") == value["base64"] and
            type(value["size"]) is int and value["size"] == len(content) and
            len(content) <= (16 * 1024 * 1024 if kind == "commit"
                             else 64 * 1024 * 1024) and
            value["sha256"] == sha256_bytes(content) and
            sha1_object(kind, content) == expected_id,
            f"{label} retained bytes differ")
    return content


def safe_git_path(value: Any, label: str) -> str:
    require(type(value) is str and value and not value.startswith("/") and
            "\0" not in value and "\r" not in value and "\n" not in value and
            all(component not in ("", ".", "..")
                for component in value.split("/")),
            f"{label} is not a safe relative Git path")
    require(len(value.encode("utf-8", errors="strict")) <= 65_536,
            f"{label} exceeds its byte bound")
    return value


def parse_git_tree(content: bytes, label: str) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    offset = 0
    while offset < len(content):
        space = content.find(b" ", offset)
        nul = content.find(b"\0", space + 1 if space >= 0 else offset)
        require(space > offset and nul > space + 1 and nul + 21 <= len(content),
                f"{label} contains a truncated tree entry")
        mode_raw, name_raw = content[offset:space], content[space + 1:nul]
        object_id = content[nul + 1:nul + 21].hex()
        offset = nul + 21
        mode = mode_raw.decode("ascii", errors="strict")
        kind = ("tree" if mode == "40000" else
                "blob" if mode in ("100644", "100755", "120000") else
                "commit" if mode == "160000" else None)
        require(kind is not None and b"/" not in name_raw,
                f"{label} contains an unsupported entry")
        name = safe_git_path(
            name_raw.decode("utf-8", errors="strict"), f"{label} entry")
        entries.append({
            "name": name, "git_mode": mode, "git_type": kind,
            "object_id": object_id,
        })
        require(len(entries) <= MAX_TRACKED_SOURCE_FILES + MAX_GIT_TREE_OBJECTS,
                f"{label} has too many entries")
    sort_keys = [entry["name"].encode("utf-8") +
                 (b"/" if entry["git_type"] == "tree" else b"")
                 for entry in entries]
    require(sort_keys == sorted(set(sort_keys)),
            f"{label} entries are duplicated or unordered")
    return entries


def flatten_git_trees(root_tree: str,
                      objects: Mapping[str, bytes]) -> list[dict[str, str]]:
    require(root_tree in objects, "Git tree closure omits its root")
    result: list[dict[str, str]] = []
    reachable: set[str] = set()
    stack: list[tuple[str, str, str, frozenset[str]]] = [
        ("tree", root_tree, "", frozenset())]
    while stack:
        kind, object_id, prefix, active = stack.pop()
        require(object_id not in active,
                "Git tree closure is cyclic")
        if kind != "tree":
            raise AuditError("internal Git tree traversal variant differs")
        require(object_id in objects, "Git tree closure is incomplete")
        reachable.add(object_id)
        next_active = active | {object_id}
        for entry in reversed(parse_git_tree(
                objects[object_id], f"Git tree object {object_id}")):
            path = f"{prefix}/{entry['name']}" if prefix else entry["name"]
            safe_git_path(path, "recursive Git tree path")
            if entry["git_type"] == "tree":
                stack.append(("tree", entry["object_id"], path, next_active))
            else:
                result.append({
                    "path": path, "git_mode": entry["git_mode"],
                    "git_type": entry["git_type"],
                    "object_id": entry["object_id"],
                })
        require(len(result) <= MAX_TRACKED_SOURCE_FILES,
                "Git tree closure has too many leaves")
    result.sort(key=lambda item: item["path"])
    require(reachable == set(objects) and
            len({item["path"] for item in result}) == len(result),
            "Git tree closure is non-canonical")
    return result


def validate_git_capture(value: Any, expected_path: str, expected_head: str,
                         require_detached: bool, label: str,
                         depth: int = 0) -> dict[str, Any]:
    require(depth <= MAX_SUBMODULE_DEPTH,
            f"{label} submodule nesting exceeds its bound")
    keys = {
        "schema", "path", "head", "tree", "detached", "head_ref",
        "superproject_worktree", "tracked_tree_listing_sha256",
        "tracked_status", "commit_object", "tree_objects", "git_executable",
        "git_metadata", "worktree_guard_policy", "config", "index",
        "tracked_files", "tracked_files_sha256", "submodules",
    }
    value = exact_keys(value, keys, f"{label} Git capture")
    require(value["schema"] == GIT_CAPTURE_SCHEMA and
            is_canonical_absolute_path(expected_path) and
            value["path"] == expected_path and
            HEX40.fullmatch(expected_head) is not None and
            value["head"] == expected_head and type(value["tree"]) is str and
            HEX40.fullmatch(value["tree"]) is not None and
            type(value["detached"]) is bool and
            (not require_detached or value["detached"] is True) and
            value["tracked_status"] == "clean" and
            type(value["tracked_tree_listing_sha256"]) is str and
            HEX256.fullmatch(value["tracked_tree_listing_sha256"]) is not None,
            f"{label} Git source identity differs")
    head_ref = value["head_ref"]
    require((value["detached"] and head_ref is None) or
            (not value["detached"] and type(head_ref) is str and
             head_ref.startswith("refs/")),
            f"{label} symbolic HEAD differs")
    superproject = value["superproject_worktree"]
    require(superproject is None or
            is_canonical_absolute_path(superproject),
            f"{label} superproject identity differs")

    commit = decode_git_object(
        value["commit_object"], "commit", expected_head, f"{label} commit")
    require(b"\n\n" in commit, f"{label} commit has no message boundary")
    headers = commit.split(b"\n\n", 1)[0].split(b"\n")
    tree_line = b"tree " + value["tree"].encode("ascii")
    require(headers and headers[0] == tree_line and
            [line for line in headers if line.startswith(b"tree ")] == [tree_line],
            f"{label} commit names another tree")
    raw_trees = value["tree_objects"]
    require(type(raw_trees) is list and 0 < len(raw_trees) <= MAX_GIT_TREE_OBJECTS,
            f"{label} tree-object closure is invalid")
    trees: dict[str, bytes] = {}
    tree_total = 0
    previous = ""
    for record in raw_trees:
        object_id = record.get("object_id") if type(record) is dict else None
        require(type(object_id) is str and HEX40.fullmatch(object_id) is not None and
                previous < object_id and object_id not in trees,
                f"{label} tree objects are unordered or duplicated")
        content = decode_git_object(
            record, "tree", object_id, f"{label} tree {object_id}")
        tree_total += len(content)
        require(tree_total <= MAX_GIT_TREE_TOTAL_BYTES,
                f"{label} tree closure exceeds its byte bound")
        trees[object_id] = content
        previous = object_id
    tree_entries = flatten_git_trees(value["tree"], trees)

    git_executable = exact_keys(
        value["git_executable"], {"source", "sealed"},
        f"{label} Git executable")
    git_source = exact_keys(git_executable["source"], {
        "path", "size", "mode", "sha256"}, f"{label} Git source executable")
    require(type(git_source["path"]) is str and
            git_source["path"] == CANONICAL_GIT_PATH and
            type(git_source["size"]) is int and git_source["size"] > 0 and
            type(git_source["mode"]) is int and
            0 <= git_source["mode"] <= 0o7777 and
            git_source["mode"] & 0o111 and
            type(git_source["sha256"]) is str and
            HEX256.fullmatch(git_source["sha256"]) is not None,
            f"{label} Git executable source differs")
    git_sealed = exact_keys(git_executable["sealed"], {
        "protocol", "size", "mode", "sha256", "seals", "source_sha256"},
        f"{label} sealed Git executable")
    require(git_sealed["protocol"] == GIT_EXECUTABLE_PROTOCOL and
            git_sealed["size"] == git_source["size"] and
            git_sealed["sha256"] == git_source["sha256"] and
            git_sealed["source_sha256"] == git_source["sha256"] and
            type(git_sealed["mode"]) is int and
            0 <= git_sealed["mode"] <= 0o7777 and
            git_sealed["mode"] == git_source["mode"] and
            git_sealed["mode"] & 0o111 and
            type(git_sealed["seals"]) is int and
            0 <= git_sealed["seals"] <= MAX_LINUX_SEAL_MASK and
            git_sealed["seals"] & REQUIRED_EXECUTABLE_SEALS ==
                REQUIRED_EXECUTABLE_SEALS,
            f"{label} sealed Git executable differs")

    metadata = exact_keys(value["git_metadata"], {
        "layout", "gitdir", "commondir", "guarded_components",
        "guard_policy", "guarded_file_count", "guarded_files_sha256"},
        f"{label} Git metadata")
    require(metadata["layout"] in ("ordinary", "linked-worktree") and
            all(is_canonical_absolute_path(metadata[name])
                for name in ("gitdir", "commondir")) and
            type(metadata["guarded_components"]) is list and
            metadata["guarded_components"] ==
                sorted(set(metadata["guarded_components"])) and
            set(metadata["guarded_components"]) == {
                metadata["gitdir"], metadata["commondir"]} and
            metadata["guard_policy"] == GIT_METADATA_GUARD_POLICY and
            type(metadata["guarded_file_count"]) is int and
            metadata["guarded_file_count"] >= 0 and
            type(metadata["guarded_files_sha256"]) is str and
            HEX256.fullmatch(metadata["guarded_files_sha256"]) is not None and
            value["worktree_guard_policy"] == GIT_WORKTREE_GUARD_POLICY,
            f"{label} Git metadata closure differs")
    validate_byte_identity(value["config"], f"{label} Git config")
    index = exact_keys(value["index"], {
        "entry_count", "stage", "flags_v", "flags_f"}, f"{label} Git index")
    require(type(index["entry_count"]) is int and
            0 <= index["entry_count"] <= MAX_TRACKED_SOURCE_FILES,
            f"{label} Git index count differs")
    for name in ("stage", "flags_v", "flags_f"):
        validate_byte_identity(index[name], f"{label} Git index {name}")

    records = value["tracked_files"]
    require(type(records) is list and len(records) == index["entry_count"] and
            value["tracked_files_sha256"] ==
                sha256_bytes(git_canonical_bytes(records)),
            f"{label} tracked-file inventory differs")
    paths: list[str] = []
    submodule_records: dict[str, Mapping[str, Any]] = {}
    for record in records:
        require(type(record) is dict and record.get("kind") in {
                    "regular", "symlink", "submodule"},
                f"{label} tracked-file record differs")
        kind = record["kind"]
        common = {"path", "git_mode", "git_type", "object_id", "kind"}
        expected = (common if kind == "regular" else
                    common | {"size", "sha256", "target_encoding",
                              "target_base64"} if kind == "symlink" else
                    common | {"identity_sha256"})
        require(set(record) == expected and
                type(record["object_id"]) is str and
                HEX40.fullmatch(record["object_id"]) is not None,
                f"{label} tracked-file shape differs")
        path = safe_git_path(record["path"], f"{label} tracked path")
        paths.append(path)
        if kind == "regular":
            require(record["git_mode"] in ("100644", "100755") and
                    record["git_type"] == "blob",
                    f"{label} regular tracked-file identity differs")
        elif kind == "symlink":
            require(record["git_mode"] == "120000" and
                    record["git_type"] == "blob" and
                    record["target_encoding"] == "base64" and
                    type(record["target_base64"]) is str,
                    f"{label} symlink identity differs")
            try:
                target = base64.b64decode(record["target_base64"], validate=True)
            except (ValueError, binascii.Error) as error:
                raise AuditError(f"{label} symlink target is invalid") from error
            require(base64.b64encode(target).decode("ascii") ==
                        record["target_base64"] and
                    type(record["size"]) is int and record["size"] == len(target) and
                    record["size"] <= MAX_TRACKED_SOURCE_FILE_BYTES and
                    record["sha256"] == sha256_bytes(target) and
                    sha1_object("blob", target) == record["object_id"],
                    f"{label} symlink bytes differ")
        else:
            require(record["git_mode"] == "160000" and
                    record["git_type"] == "commit" and
                    type(record["identity_sha256"]) is str and
                    HEX256.fullmatch(record["identity_sha256"]) is not None,
                    f"{label} submodule identity differs")
            submodule_records[path] = record
    require(paths == sorted(set(paths)),
            f"{label} tracked files are unordered or duplicated")
    projected = [{name: record[name] for name in (
        "path", "git_mode", "git_type", "object_id")} for record in records]
    require(projected == tree_entries,
            f"{label} tracked files differ from retained tree objects")
    listing = b"".join((
        f"{record['git_mode']} {record['git_type']} {record['object_id']}\t"
        f"{record['path']}\0").encode("utf-8") for record in records)
    stage = b"".join((
        f"{record['git_mode']} {record['object_id']} 0\t{record['path']}\0"
    ).encode("utf-8") for record in records)
    flags = b"".join(
        f"H {record['path']}\0".encode("utf-8") for record in records)
    expected_bytes = lambda data: {"size": len(data), "sha256": sha256_bytes(data)}
    require(value["tracked_tree_listing_sha256"] == sha256_bytes(listing) and
            index["stage"] == expected_bytes(stage) and
            index["flags_v"] == expected_bytes(flags) and
            index["flags_f"] == expected_bytes(flags),
            f"{label} Git tree/index transcripts differ")
    submodules = value["submodules"]
    require(type(submodules) is list and len(submodules) == len(submodule_records),
            f"{label} submodule inventory differs")
    observed: list[str] = []
    for record in submodules:
        record = exact_keys(record, {
            "path", "object_id", "identity_sha256", "identity"},
            f"{label} submodule")
        path = record["path"]
        tracked = submodule_records.get(path)
        require(tracked is not None and record["object_id"] == tracked["object_id"] and
                record["identity_sha256"] == tracked["identity_sha256"] and
                record["identity_sha256"] ==
                    sha256_bytes(git_canonical_bytes(record["identity"])),
                f"{label} submodule cross-binding differs")
        nested_path = os.path.abspath(os.path.join(value["path"], path))
        validate_git_capture(record["identity"], nested_path,
                             record["object_id"], False,
                             f"{label} submodule {path}", depth + 1)
        observed.append(path)
    require(observed == sorted(set(observed)),
            f"{label} submodule order differs")
    return value


def parse_canonical_ldd(text: str, label: str) -> list[dict[str, Any]]:
    address = re.escape(CANONICAL_LDD_ADDRESS)
    shared = re.compile(
        rf"(?P<soname>[^\s=()]+)\s+=>\s+(?P<path>/\S+)\s+"
        rf"\((?P<address>{address})\)")
    direct = re.compile(rf"(?P<path>/\S+)\s+\((?P<address>{address})\)")
    virtual = re.compile(
        rf"(?P<soname>linux-vdso\.so\.1)\s+\((?P<address>{address})\)")
    entries: list[dict[str, Any]] = []
    rendered: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        require(line, f"{label} contains an empty line")
        match = shared.fullmatch(line)
        if match is not None:
            path = match.group("path")
            require(os.path.normpath(path) == path,
                    f"{label} shared-library path is noncanonical")
            entry = {
                "soname": match.group("soname"), "loader_path": path,
                "file_kind": "shared_library",
            }
            rendered.append(
                f"{entry['soname']} => {path} ({CANONICAL_LDD_ADDRESS})")
        else:
            match = direct.fullmatch(line)
            if match is not None:
                path = match.group("path")
                require(os.path.normpath(path) == path,
                        f"{label} dynamic-loader path is noncanonical")
                entry = {
                    "soname": Path(path).name, "loader_path": path,
                    "file_kind": "dynamic_loader",
                }
                rendered.append(f"{path} ({CANONICAL_LDD_ADDRESS})")
            else:
                match = virtual.fullmatch(line)
                require(match is not None, f"{label} has an unrecognized line")
                entry = {"soname": match.group("soname"), "virtual": True}
                rendered.append(
                    f"{entry['soname']} ({CANONICAL_LDD_ADDRESS})")
        entries.append(entry)
    require(entries and text == "\n".join(rendered) + "\n" and
            len({item["soname"] for item in entries}) == len(entries) and
            sum(item.get("file_kind") == "dynamic_loader"
                for item in entries) == 1 and
            any(item.get("file_kind") == "shared_library" for item in entries) and
            sum(item.get("virtual") is True for item in entries) <= 1,
            f"{label} canonical runtime closure differs")
    return sorted(entries, key=lambda item: item["soname"])


def validate_runtime_closure(value: Any, label: str,
                             expected_executable: str) -> dict[str, Any]:
    value = exact_keys(value, {
        "executable", "dependencies", "canonical_ldd_output"}, label)
    require(value["executable"] == expected_executable and
            type(value["dependencies"]) is list and value["dependencies"],
            f"{label} shape differs")
    output = exact_keys(value["canonical_ldd_output"], {
        "schema", "normalization", "encoding", "size", "sha256", "text"},
        f"{label} canonical ldd output")
    require(output["schema"] == CANONICAL_LDD_SCHEMA and
            output["normalization"] == CANONICAL_LDD_NORMALIZATION,
            f"{label} canonical ldd contract differs")
    validate_text({name: output[name] for name in (
        "encoding", "size", "sha256", "text")},
        f"{label} canonical ldd text")
    parsed = parse_canonical_ldd(output["text"], label)
    sonames: list[str] = []
    normalized: list[dict[str, Any]] = []
    loader_paths: list[str] = []
    for index, dependency in enumerate(value["dependencies"]):
        require(type(dependency) is dict and type(dependency.get("soname")) is str and
                dependency["soname"],
                f"{label} dependency {index} differs")
        sonames.append(dependency["soname"])
        if set(dependency) == {"soname", "virtual"}:
            require(dependency == {"soname": "linux-vdso.so.1", "virtual": True},
                    f"{label} virtual dependency differs")
            normalized.append(dict(dependency))
            continue
        require(set(dependency) == {"soname", "loader_path", "file"} and
                type(dependency["loader_path"]) is str and
                Path(dependency["loader_path"]).is_absolute(),
                f"{label} dependency variant differs")
        artifact = validate_artifact(
            dependency["file"], f"{label} dependency {index}")
        require(artifact["kind"] in {"shared_library", "dynamic_loader"} and
                artifact["path"] == dependency["loader_path"],
                f"{label} dependency file differs")
        loader_paths.append(dependency["loader_path"])
        normalized.append({
            "soname": dependency["soname"],
            "loader_path": dependency["loader_path"],
            "file_kind": artifact["kind"],
        })
    require(sonames == sorted(set(sonames)) and
            loader_paths == list(dict.fromkeys(loader_paths)) and
            normalized == parsed,
            f"{label} dependencies differ from retained ldd bytes")
    return value


def candidate_required_cache() -> dict[str, str]:
    return {
        "ENABLE_OPENMP": "ON",
        "LEO2_BACKEND_VARIANT": "auto",
        "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_BUILD_TESTS": "OFF",
        "LEO2_ENABLE_CUDA": "OFF",
        "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
        "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
        "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
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
        "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
        "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
        "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
        "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
        "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
        "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
        "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": "ON",
        "LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED": "ON",
        "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
        "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": "ON",
    }


def validate_release_flags(value: str, role: str) -> None:
    try:
        tokens = shlex.split(value, posix=True)
    except ValueError as error:
        raise AuditError(f"{role} Release flags are invalid: {error}") from error
    # This is the CMake cache value, before the campaign appends its final
    # optimization override to the generated sidecar and compiler argv.
    expected = (["-g", "-O0", "-O3"] if role == "baseline" else
                ["-O3", "-DNDEBUG"])
    require(tokens == expected,
            f"{role} global Release flags differ from the frozen build")


def validate_build_layout(cache: Mapping[str, Any], role: str) -> str | None:
    generator = cache.get("CMAKE_GENERATOR")
    build_type = cache.get("CMAKE_BUILD_TYPE")
    configurations = cache.get("CMAKE_CONFIGURATION_TYPES")
    require(type(generator) is str and type(build_type) is str and
            type(configurations) is str,
            f"{role} CMake layout fields differ")
    # The independent auditor deliberately excludes multi-config builds.  The
    # retained Ninja graph does not itself retain Ninja's expanded command
    # transcript, so accepting separately reported link recipes would leave a
    # second editable truth.  Authoritative v17 evidence uses this exact
    # single-config layout.
    require(generator == "Unix Makefiles" and
            build_type == "Release" and
            "CMAKE_MAKE_PROGRAM" not in cache,
            f"{role} is not the frozen Unix Makefiles Release layout")
    return None


def parse_build_configuration(text: str) -> dict[str, Any]:
    data = text.encode("utf-8", errors="strict")
    require(0 < len(data) <= 64 * 1024 and b"\0" not in data and
            b"\r" not in data and data.endswith(b"\n"),
            "candidate effective build-configuration bytes are invalid")
    lines = text[:-1].split("\n")
    require(len(lines) == len(BUILD_CONFIGURATION_VARIABLES) + 2 and
            lines[0] == f"schema={BUILD_CONFIGURATION_FILE_SCHEMA}" and
            lines[1].startswith("sha256="),
            "candidate effective build-configuration header differs")
    digest = lines[1].removeprefix("sha256=")
    require(HEX256.fullmatch(digest) is not None,
            "candidate effective build-configuration digest is invalid")
    entries: dict[str, str] = {}
    for expected, line in zip(BUILD_CONFIGURATION_VARIABLES, lines[2:]):
        name, separator, value = line.partition("=")
        require(separator == "=" and name == expected and name not in entries,
                "candidate effective build-configuration ordering differs")
        entries[name] = value
    material = "".join(
        f"{name}={entries[name]}\n" for name in BUILD_CONFIGURATION_VARIABLES
    ).encode("utf-8")
    require(sha256_bytes(material) == digest,
            "candidate effective build-configuration material digest differs")
    return {"configuration_sha256": digest, "entries": entries}


def validate_effective_configuration(value: Any, build_dir: str,
                                     source_root: str,
                                     benchmark_object: Mapping[str, Any],
                                     cache: Mapping[str, Any]) -> dict[str, Any]:
    value = exact_keys(value, {
        "schema", "artifact", "content", "configuration_schema",
        "configuration_sha256", "entries", "embedded_build_type",
        "helper_source"}, "candidate effective build configuration")
    require(value["schema"] == BUILD_CONFIGURATION_RECORD_SCHEMA and
            value["configuration_schema"] == BUILD_CONFIGURATION_FILE_SCHEMA and
            value["embedded_build_type"] == "Release",
            "candidate effective build-configuration contract differs")
    artifact = validate_artifact(
        value["artifact"], "candidate effective build configuration",
        "generated_build_configuration")
    require(artifact["path"] == str(Path(build_dir) /
                                     BUILD_CONFIGURATION_RELATIVE_PATH),
            "candidate effective build-configuration path differs")
    content = validate_text(
        value["content"], "candidate effective build-configuration content")
    parsed = parse_build_configuration(content["text"])
    require(content["size"] == artifact["size"] and
            content["sha256"] == artifact["sha256"] and
            value["configuration_sha256"] == parsed["configuration_sha256"] and
            exact_json_equal(value["entries"], parsed["entries"]),
            "candidate effective build-configuration bytes differ")
    entries = parsed["entries"]
    fixed = {
        "CMAKE_CXX_FLAGS": " -Wall -Wextra -fopenmp",
        "CMAKE_CXX_FLAGS_DEBUG": "-g -O0",
        "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG -O3",
        "CMAKE_CXX_FLAGS_RELWITHDEBINFO": "-O2 -g -DNDEBUG",
        "CMAKE_CXX_FLAGS_MINSIZEREL": "-Os -DNDEBUG",
        "ENABLE_OPENMP": "ON", "LEOPARD_ENABLE_GF8": "ON",
        "LEOPARD_ENABLE_GF16": "ON", "LEO2_BACKEND_VARIANT": "auto",
        "LEO2_BENCHMARK_GIT_EXECUTABLE": "/usr/bin/git",
        "LEO2_BUILD_BENCHMARKS": "ON", "LEO2_BUILD_TESTS": "OFF",
        **{name: item for name, item in candidate_required_cache().items()
           if name in BUILD_CONFIGURATION_VARIABLES},
    }
    require(set(entries) == set(BUILD_CONFIGURATION_VARIABLES) and
            all(entries.get(name) == item for name, item in fixed.items()) and
            entries["CMAKE_CXX_COMPILER"] == cache["CMAKE_CXX_COMPILER"] and
            entries["CMAKE_GENERATOR"] == cache["CMAKE_GENERATOR"] and
            entries["CMAKE_CONFIGURATION_TYPES"] ==
                cache["CMAKE_CONFIGURATION_TYPES"] and
            entries["CMAKE_BUILD_TYPE"] == cache["CMAKE_BUILD_TYPE"],
            "candidate effective build-configuration entries differ")
    helper = validate_artifact(
        value["helper_source"], "candidate attestation helper", "source_file")
    require(helper["path"] == str(Path(source_root) /
                                    BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH) and
            benchmark_object["mtime_ns"] >= artifact["mtime_ns"] and
            benchmark_object["mtime_ns"] >= helper["mtime_ns"],
            "candidate configuration/helper freshness differs")
    return value


def attestation_header_bytes(commit: str, tree: str) -> bytes:
    return (
        "#ifndef LEOPARD2_BENCHMARK_SOURCE_ATTESTATION_GENERATED_H\n"
        "#define LEOPARD2_BENCHMARK_SOURCE_ATTESTATION_GENERATED_H\n\n"
        "#undef LEO2_BENCHMARK_SOURCE_COMMIT\n"
        "#undef LEO2_BENCHMARK_SOURCE_TREE\n"
        "#undef LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY\n"
        f"#define LEO2_BENCHMARK_SOURCE_COMMIT \"{commit}\"\n"
        f"#define LEO2_BENCHMARK_SOURCE_TREE \"{tree}\"\n"
        "#define LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY 0\n\n"
        "#endif\n"
    ).encode("ascii")


def validate_generated_attestation(value: Any, build_dir: str,
                                   benchmark_object: Mapping[str, Any],
                                   commit: str, tree: str) -> dict[str, Any]:
    value = exact_keys(value, {
        "artifact", "content", "source_commit", "source_tree",
        "source_tracked_dirty"}, "candidate generated attestation")
    require(value["source_commit"] == commit and value["source_tree"] == tree and
            value["source_tracked_dirty"] is False,
            "candidate generated-attestation source differs")
    artifact = validate_artifact(
        value["artifact"], "candidate generated attestation",
        "generated_compile_input")
    expected_path = str(
        Path(build_dir) / "generated/leopard2-benchmark-attestation/"
        "leopard2_benchmark_source_attestation.h")
    require(artifact["path"] == expected_path,
            "candidate generated-attestation path differs")
    content = validate_text(value["content"],
                            "candidate generated-attestation content")
    expected = attestation_header_bytes(commit, tree)
    require(content["text"].encode("utf-8") == expected and
            content["size"] == artifact["size"] and
            content["sha256"] == artifact["sha256"] and
            benchmark_object["mtime_ns"] >= artifact["mtime_ns"],
            "candidate generated-attestation bytes/freshness differ")
    return value


def candidate_target(relative: str) -> str:
    return CANDIDATE_TARGETS.get(relative, "leopard.dir")


def expected_compile_output(role: str, source: str,
                            specification: Mapping[str, Any],
                            selected: str | None) -> str:
    configuration = f"{selected}/" if selected else ""
    path = Path(source)
    if role == "baseline":
        adapter = Path(specification["candidate_source_root"]) / \
            "experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        if path == adapter:
            return ("CMakeFiles/leopard_main_benchmark.dir/" + configuration +
                    "legacy_main_benchmark.cpp.o")
        require(path.is_relative_to(Path(specification["baseline_source_root"])),
                "baseline compile source escapes exact-main root")
        return ("CMakeFiles/leopard_main_exact.dir/" + configuration +
                source.removeprefix("/") + ".o")
    root = Path(specification["candidate_source_root"])
    require(path.is_relative_to(root),
            "candidate compile source escapes source root")
    relative = path.relative_to(root).as_posix()
    target = ("bench_leopard2.dir" if relative == "bench/leopard2/benchmark.cpp"
              else candidate_target(relative))
    return f"CMakeFiles/{target}/{configuration}{relative}.o"


def expected_candidate_compile_argv(
        source: str, specification: Mapping[str, Any], compiler: str,
        selected: str | None, configuration: Mapping[str, Any],
        compiler_path: str) -> list[str]:
    root = Path(specification["candidate_source_root"])
    relative = Path(source).relative_to(root).as_posix()
    target = ("bench_leopard2.dir" if relative == "bench/leopard2/benchmark.cpp"
              else candidate_target(relative))
    output = expected_compile_output("candidate", source, specification, selected)
    propagated: list[str] = []
    includes = [f"-I{root}"]
    if relative == "bench/leopard2/benchmark.cpp":
        digest = configuration["configuration_sha256"]
        header = Path(specification["candidate_build_dir"]) / \
            "generated/leopard2-benchmark-attestation/" \
            "leopard2_benchmark_source_attestation.h"
        definitions = [
            f'-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256="{digest}"',
            '-DLEO2_BENCHMARK_BUILD_TYPE="Release"',
            "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1",
            f'-DLEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER="{header}"',
            "-DLEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED=1",
        ]
        propagated = ["-fopenmp"]
    elif relative == "Leopard2BackendGFNI.cpp":
        definitions = [
            "-DLEO2_HAVE_AVX2_BACKEND=1", "-DLEO2_HAVE_GFNI_BACKEND=1"]
    elif relative == "Leopard2BackendAVX2T32B256.cpp":
        definitions = [
            "-DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=0",
            "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1",
            "-DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK=0",
            "-DLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
        ]
    elif relative == "Leopard2BackendAVX2T16B64.cpp":
        definitions = [
            "-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
        ]
    elif relative == "Leopard2LowP32B64AVX2.cpp":
        definitions = [
            "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
        ]
    elif relative == "Leopard2BackendAVX2T2K4.cpp":
        definitions = ["-DLEO2_HAVE_AVX2_BACKEND=1"]
    elif relative == "Leopard2BackendAVX2T8K8B1024.cpp":
        definitions = [
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1",
        ]
    elif relative == "Leopard2BackendAVX2T16K66.cpp":
        definitions = [
            "-DLEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
        ]
    elif relative in CANDIDATE_ISA_FLAGS:
        definitions = []
    elif relative in CANDIDATE_LIBRARY_SOURCES:
        definitions = [
            "-DLEO2_DISABLE_AVX2_CODEGEN=1",
            "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_HAVE_GFNI_BACKEND=1",
            "-DLEO2_HAVE_SSSE3_BACKEND=1",
        ]
        propagated = ["-fopenmp"]
    else:
        raise AuditError(f"unexpected candidate compile source {relative}")
    definitions.extend([
        "-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=1",
        "-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=1",
    ])
    if relative in CANDIDATE_LIBRARY_SOURCES and target == "leopard.dir":
        definitions.extend([
            "-DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1",
            "-DLEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1",
            "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1",
            "-DLEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED=1",
        ])
    elif relative in {
            "Leopard2BackendAVX2.cpp", "Leopard2BackendAVX2Xor.cpp",
            "Leopard2BackendAVX2T16Q2.cpp", "Leopard2BackendAVX2T8K62.cpp"}:
        definitions.extend([
            "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1",
            "-DLEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED=1",
        ])
    definitions.sort()
    source_definitions: list[str] = []
    if relative in {"leopard2.cpp", "Leopard2Backend.cpp",
                    "Leopard2BackendAVX2.cpp"}:
        source_definitions.append(
            "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1")
    if relative == "leopard2.cpp":
        source_definitions.extend([
            "-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=1",
            "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1",
            "-DLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=1",
            "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=1",
            "-DLEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS=1",
            "-DLEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK=1",
        ])
    source_definitions.sort()
    configuration_definition = (
        [f'-DCMAKE_INTDIR="{selected}"'] if selected else [])
    isa = list(CANDIDATE_ISA_FLAGS.get(relative, ()))
    if relative == "Leopard2BackendAVX2T8K8B1024.cpp":
        require(GNU_CXX.fullmatch(Path(compiler_path).name) is not None,
                "v17 candidate compiler is not the GNU driver required by argv")
    return [
        compiler, *definitions, *configuration_definition,
        *source_definitions, *includes, "-Wall", "-Wextra", "-fopenmp",
        "-O3", "-DNDEBUG", "-O3", "-std=gnu++11", *isa, *propagated,
        "-o", output, "-c", source,
    ]


def expected_compile_argv(role: str, source: str,
                          specification: Mapping[str, Any], compiler: str,
                          selected: str | None,
                          configuration: Mapping[str, Any] | None,
                          compiler_path: str) -> list[str]:
    if role == "candidate":
        require(configuration is not None,
                "candidate effective configuration is absent")
        return expected_candidate_compile_argv(
            source, specification, compiler, selected, configuration,
            compiler_path)
    adapter = str(Path(specification["candidate_source_root"]) /
                  "experiments/leopard2/main_compare/legacy_main_benchmark.cpp")
    definitions = ([] if source != adapter else [
        f'-DLEOPARD_MAIN_SOURCE_COMMIT="{MAIN_COMMIT}"',
        "-DLEO_MAIN_PURE_AVX2_PROFILE=0",
    ])
    configuration_definition = (
        [f'-DCMAKE_INTDIR="{selected}"'] if selected else [])
    output = expected_compile_output("baseline", source, specification, selected)
    return [
        compiler, *definitions, *configuration_definition,
        f"-I{specification['baseline_source_root']}",
        "-g", "-O0", "-O3", "-std=gnu++11", "-march=native",
        "-Wall", "-Wextra", "-fopenmp", "-o", output, "-c", source,
    ]


def validate_external_link_inputs(value: Any, label: str) -> list[dict[str, Any]]:
    roles = (
        ("openmp_runtime_shared", "libgomp.so", "shared_library"),
        ("pthread_support_archive", "libpthread.a", "archive"),
    )
    require(type(value) is list and len(value) == len(roles),
            f"{label} external-link closure differs")
    operands: list[str] = []
    resolved_paths: list[str] = []
    for index, (role, basename, kind) in enumerate(roles):
        record = exact_keys(value[index], {
            "operand", "role", "artifact", "lexical_symlink_chain"},
            f"{label} external input {index}")
        require(record["role"] == role and
                is_canonical_absolute_path(record["operand"]) and
                Path(record["operand"]).name == basename and
                not any(character.isspace() for character in record["operand"]) and
                "@" not in record["operand"] and "\\" not in record["operand"],
                f"{label} external operand {index} differs")
        if role == "openmp_runtime_shared":
            require(re.fullmatch(
                        r"/usr/lib/gcc/x86_64-linux-gnu/"
                        r"[0-9]+(?:\.[0-9]+)*/libgomp\.so",
                        record["operand"]) is not None,
                    f"{label} OpenMP operand is outside canonical GCC root")
        else:
            require(re.fullmatch(
                        r"/(?:usr/)?lib(?:64|/x86_64-linux-gnu)/libpthread\.a",
                        record["operand"]) is not None,
                    f"{label} pthread operand is outside canonical system root")
        artifact = validate_artifact(
            record["artifact"], f"{label} external artifact {index}", kind)
        require(artifact["size"] > 0,
                f"{label} external artifact {index} is empty")
        chain = record["lexical_symlink_chain"]
        require(type(chain) is list and len(chain) <= 40 and all(
                    type(link) is dict and set(link) == {"path", "target", "mode"} and
                    is_canonical_absolute_path(link["path"]) and
                    type(link["target"]) is str and link["target"] and
                    not any(character in link["target"] for character in "\0\r\n") and
                    type(link["mode"]) is int and 0 <= link["mode"] <= 0o7777
                    for link in chain),
                f"{label} external symlink chain {index} differs")
        require(len({link["path"] for link in chain}) == len(chain),
                f"{label} external symlink chain {index} repeats a path")
        if role == "pthread_support_archive":
            require(chain == [] and artifact["path"] == record["operand"],
                    f"{label} pthread archive lexical identity differs")
        else:
            require(chain and chain[0]["path"] == record["operand"],
                    f"{label} OpenMP symlink chain does not start at operand")
            current = Path(record["operand"])
            for link_index, link in enumerate(chain):
                require(str(current) == link["path"],
                        f"{label} OpenMP symlink chain {link_index} is discontinuous")
                target = Path(link["target"])
                current = (target if target.is_absolute()
                           else current.parent / target)
                current = Path(os.path.normpath(str(current)))
            require(str(current) == artifact["path"] and
                    re.fullmatch(
                        r"/(?:usr/)?lib(?:64|/x86_64-linux-gnu|/gcc/"
                        r"x86_64-linux-gnu/[0-9]+(?:\.[0-9]+)*)/"
                        r"libgomp\.so\.[0-9]+(?:\.[0-9]+)*",
                        artifact["path"]) is not None,
                    f"{label} OpenMP symlink chain/resolved artifact differs")
        operands.append(record["operand"])
        resolved_paths.append(artifact["path"])
    require(len(set(operands)) == len(operands) and
            len(set(resolved_paths)) == len(resolved_paths),
            f"{label} external operands/resolved files are not unique")
    return value


def ninja_relative(value: Any, label: str) -> str:
    require(type(value) is str and value and
            re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._/+@=-]*", value) is not None,
            f"{label} is not a literal Ninja path")
    path = Path(value)
    require(not path.is_absolute() and path.as_posix() == value and
            all(part not in ("", ".", "..") for part in path.parts),
            f"{label} is not a canonical relative Ninja path")
    return value


def ninja_includes(text: str, label: str) -> tuple[str, ...]:
    require(type(text) is str and text.endswith("\n") and
            "\0" not in text and "\r" not in text,
            f"{label} Ninja text differs")
    result: list[str] = []
    for line_number, line in enumerate(text[:-1].split("\n"), 1):
        if re.match(r"^(include|subninja)(?:[ \t]|$)", line) is None:
            continue
        match = re.fullmatch(r"(?:include|subninja)[ \t]+([^ \t]+)[ \t]*", line)
        require(match is not None, f"{label}:{line_number} include is not literal")
        result.append(ninja_relative(
            match.group(1), f"{label}:{line_number} include"))
    return tuple(result)


def validate_ninja_graph(value: Any, build_dir: str,
                         selected: str, label: str) -> dict[str, Any]:
    value = exact_keys(value, {"schema", "entrypoint", "files"}, label)
    require(value["schema"] == NINJA_GRAPH_CLOSURE_SCHEMA and
            value["entrypoint"] == f"build-{selected}.ninja" and
            type(value["files"]) is list and
            0 < len(value["files"]) <= MAX_NINJA_GRAPH_FILES,
            f"{label} shape differs")
    records: dict[str, dict[str, Any]] = {}
    total = 0
    for record in value["files"]:
        record = exact_keys(record, {"relative_path", "artifact", "content"},
                            f"{label} record")
        relative = ninja_relative(record["relative_path"], f"{label} path")
        require(relative not in records, f"{label} path is duplicated")
        artifact = validate_artifact(
            record["artifact"], f"{label} {relative}", "ninja_graph_input")
        content = validate_text(record["content"], f"{label} {relative} content")
        total += content["size"]
        require(artifact["path"] == str(Path(build_dir).joinpath(
                    *Path(relative).parts)) and
                artifact["size"] == content["size"] and
                artifact["sha256"] == content["sha256"] and
                total <= MAX_NINJA_GRAPH_TOTAL_BYTES,
                f"{label} {relative} bytes differ")
        records[relative] = record
    require(list(records) == sorted(records), f"{label} records are unordered")
    pending = [value["entrypoint"]]
    visited: set[str] = set()
    while pending:
        relative = pending.pop()
        require(relative in records,
                f"{label} references an unretained input {relative}")
        if relative in visited:
            continue
        visited.add(relative)
        pending.extend(ninja_includes(
            records[relative]["content"]["text"], f"{label} {relative}"))
    require(visited == set(records), f"{label} contains unreachable inputs")
    return value


def validate_archive_recipe(content: Any, identity: Mapping[str, Any] | None,
                            archive_name: str, target_dir: str,
                            objects: Sequence[str], archiver: str,
                            ranlib: str, label: str) -> None:
    content = validate_text(content, label)
    if identity is not None:
        require(content["size"] == identity["size"] and
                content["sha256"] == identity["sha256"],
                f"{label} differs from its artifact identity")
    commands: list[list[str]] = []
    for line in content["text"].splitlines():
        if not line.strip():
            continue
        try:
            command = shlex.split(line, posix=True)
        except ValueError as error:
            raise AuditError(f"cannot parse {label}: {error}") from error
        require(command, f"{label} contains an empty command")
        commands.append(command)
    require(len(commands) == 2, f"{label} command count differs")
    ar_command, ranlib_command = commands
    require(len(ar_command) >= 4 and ar_command[0] == archiver and
            ar_command[1] in {"qc", "rc", "rcs"} and
            ar_command[2] == archive_name and
            ar_command[3:] == list(objects) and
            all(item.endswith(".o") and not item.startswith("/") and
                "@" not in item and "\\" not in item for item in objects) and
            ranlib_command == [ranlib, archive_name],
            f"{label} tool/output/object closure differs")
    ordinary_prefix = f"CMakeFiles/{target_dir}/"
    ordinary = [item for item in objects if
                "/leopard2_backend_" not in f"/{item}" and
                "/leopard2_low_p32_b64_avx2.dir/" not in f"/{item}"]
    require(ordinary and all(ordinary_prefix in item for item in ordinary),
            f"{label} ordinary object target differs")


def validate_executable_recipe(content: Any,
                               identity: Mapping[str, Any] | None,
                               compiler: str, archive: str, executable: str,
                               benchmark_object: str,
                               external: Sequence[Mapping[str, Any]],
                               label: str) -> None:
    content = validate_text(content, label)
    if identity is not None:
        require(content["size"] == identity["size"] and
                content["sha256"] == identity["sha256"],
                f"{label} differs from its artifact identity")
    lines = content["text"].splitlines()
    require(len(lines) == 1 and lines[0], f"{label} command count differs")
    try:
        tokens = shlex.split(lines[0], posix=True)
    except ValueError as error:
        raise AuditError(f"cannot parse {label}: {error}") from error
    require(tokens and tokens[0] == compiler and tokens.count("-o") == 1,
            f"{label} compiler/output grammar differs")
    output_index = tokens.index("-o")
    require(output_index + 1 < len(tokens) and
            tokens[output_index + 1] == executable and
            tokens.count(archive) == tokens.count(benchmark_object) == 1,
            f"{label} project input/output closure differs")
    external_operands = [record["operand"] for record in external]
    require(all(tokens.count(operand) == 1 for operand in external_operands),
            f"{label} external input closure differs")
    consumed = {0, output_index, output_index + 1}
    for required in (archive, benchmark_object, *external_operands):
        consumed.add(tokens.index(required))
    allowed_flags = {
        "-Wall", "-Wextra", "-Wpedantic", "-fopenmp", "-g", "-DNDEBUG",
        "-O0", "-O1", "-O2", "-O3", "-Os", "-Oz",
    }
    require(all(index in consumed or token in allowed_flags
                for index, token in enumerate(tokens)) and
            not any("@" in token or token.startswith("-l") or ",-l" in token
                    for token in tokens),
            f"{label} contains an undeclared option or file input")


def validate_source_object_pair(value: Any, label: str,
                                build_dir: str) -> dict[str, Any]:
    value = exact_keys(value, {"source", "object"}, label)
    source = validate_artifact(value["source"], f"{label} source", "source_file")
    object_artifact = validate_artifact(
        value["object"], f"{label} object", "object_file")
    require(Path(object_artifact["path"]).is_relative_to(Path(build_dir)) and
            object_artifact["mtime_ns"] >= source["mtime_ns"],
            f"{label} object escapes the build or predates its source")
    return value


def validate_archive_tools(archiver_value: Any, ranlib_value: Any,
                           invocation_value: Any, label: str) \
        -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    archiver = validate_artifact(
        archiver_value, f"{label} archiver", "archiver")
    ranlib = validate_artifact(ranlib_value, f"{label} ranlib", "ranlib")
    require(archiver["path"] == CANONICAL_ARCHIVER_RESOLVED_PATH and
            ranlib["path"] == CANONICAL_RANLIB_RESOLVED_PATH,
            f"{label} resolved archive tools differ from the frozen host")
    tools = exact_keys(invocation_value, {"archiver", "ranlib"},
                       f"{label} archive tools")
    for name, artifact, canonical in (
            ("archiver", archiver, CANONICAL_ARCHIVER_PATH),
            ("ranlib", ranlib, CANONICAL_RANLIB_PATH)):
        item = exact_keys(tools[name], {"invocation", "resolved_path"},
                          f"{label} {name} invocation")
        # CMake retains the lexical /usr/bin invocation while the artifact is
        # the resolved executable (for example x86_64-linux-gnu-ar).
        require(item["invocation"] == canonical and
                item["resolved_path"] == artifact["path"],
                f"{label} {name} invocation differs")
    return archiver, ranlib, tools


def validate_build(build: Any, role: str, specification: Mapping[str, Any],
                   candidate_tree: str | None) -> dict[str, Any]:
    build = exact_keys(build, BUILD_KEYS, f"{role} build")
    build_dir = build["build_dir"]
    require(type(build_dir) is str and Path(build_dir).is_absolute() and
            build_dir == specification[f"{role}_build_dir"],
            f"{role} build directory differs")
    cache = build["validated_cache"]
    require(type(cache) is dict, f"{role} validated CMake cache is absent")
    selected = validate_build_layout(cache, role)
    configuration_count = 2 if selected else 1
    output_root = Path(build_dir) / selected if selected else Path(build_dir)
    if role == "baseline":
        executable_name, archive_name = (
            "leopard_main_benchmark", "libleopard_main_exact.a")
        archive_target = "leopard_main_exact.dir"
        profile, isa_policy = BASELINE_COMPILE_PROFILE, BASELINE_ISA_POLICY
        library_sources = BASELINE_LIBRARY_SOURCES
        benchmark_relative = \
            "experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        benchmark_source = str(Path(specification["candidate_source_root"]) /
                               benchmark_relative)
        required_cache: dict[str, Any] = {
            "CMAKE_BUILD_TYPE": cache.get("CMAKE_BUILD_TYPE"),
            "CMAKE_CXX_COMPILER": cache.get("CMAKE_CXX_COMPILER"),
            "CMAKE_CXX_FLAGS_RELEASE": "-g -O0 -O3",
            "LEO_MAIN_HAS_MARCH_NATIVE": "1",
            "LEO_MAIN_PURE_AVX2": "OFF",
            "CMAKE_CONFIGURATION_TYPES": cache.get("CMAKE_CONFIGURATION_TYPES"),
            "CMAKE_GENERATOR": cache.get("CMAKE_GENERATOR"),
        }
        base_entry_count = 5
    else:
        executable_name, archive_name = "bench_leopard2", "libleopard.a"
        archive_target = "leopard.dir"
        profile, isa_policy = CANDIDATE_COMPILE_PROFILE, CANDIDATE_ISA_POLICY
        library_sources = CANDIDATE_LIBRARY_SOURCES
        benchmark_relative = "bench/leopard2/benchmark.cpp"
        benchmark_source = str(Path(specification["candidate_source_root"]) /
                               benchmark_relative)
        compile_record_for_cache = build.get("validated_compile_commands")
        configuration_for_cache = (
            compile_record_for_cache.get("effective_build_configuration")
            if type(compile_record_for_cache) is dict else None)
        configuration_digest = (
            configuration_for_cache.get("configuration_sha256")
            if type(configuration_for_cache) is dict else None)
        required_cache = {
            "CMAKE_BUILD_TYPE": cache.get("CMAKE_BUILD_TYPE"),
            "CMAKE_CXX_COMPILER": cache.get("CMAKE_CXX_COMPILER"),
            "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
            **candidate_required_cache(),
            "CMAKE_CONFIGURATION_TYPES": cache.get("CMAKE_CONFIGURATION_TYPES"),
            "CMAKE_GENERATOR": cache.get("CMAKE_GENERATOR"),
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                BUILD_CONFIGURATION_FILE_SCHEMA,
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
                configuration_digest,
        }
        base_entry_count = len(CANDIDATE_LIBRARY_SOURCES) + \
            len(CANDIDATE_NON_LIBRARY_ACTIONS)
    if selected:
        required_cache["CMAKE_MAKE_PROGRAM"] = CANONICAL_NINJA_PATH
    require(set(cache) == set(required_cache) and
            all(cache.get(name) == expected
                for name, expected in required_cache.items()) and
            type(cache.get("CMAKE_CXX_COMPILER")) is str and
            cache["CMAKE_CXX_COMPILER"] and
            (role == "baseline" or
             HEX256.fullmatch(str(configuration_digest)) is not None),
            f"{role} validated CMake cache shape/semantics differ")
    validate_release_flags(cache["CMAKE_CXX_FLAGS_RELEASE"], role)

    metadata_paths = {
        "cmake_cache": str(Path(build_dir) / "CMakeCache.txt"),
        "compile_commands": str(Path(build_dir) / "compile_commands.json"),
        "executable_link_recipe": str(Path(build_dir) / (
            f"CMakeFiles/impl-{selected}.ninja" if selected else
            f"CMakeFiles/{executable_name}.dir/link.txt")),
        "archive_link_recipe": str(Path(build_dir) / (
            f"CMakeFiles/impl-{selected}.ninja" if selected else
            f"CMakeFiles/{archive_target}/link.txt")),
    }
    for name, path in metadata_paths.items():
        artifact = validate_artifact(
            build[name], f"{role} {name}", "build_metadata")
        require(artifact["path"] == path, f"{role} {name} path differs")
    compiler = validate_artifact(build["compiler"], f"{role} compiler", "compiler")
    archiver, ranlib, tools = validate_archive_tools(
        build["archiver"], build["ranlib"],
        build["archive_link_tool_invocations"], role)
    executable = validate_artifact(
        build["validated_executable"], f"{role} executable", "executable")
    archive = validate_artifact(
        build["validated_archive"], f"{role} archive", "archive")
    require(executable["path"] == str(output_root / executable_name) and
            archive["path"] == str(output_root / archive_name),
            f"{role} validated output paths differ")
    validate_version_text(build["compiler_version_stdout"],
                          f"{role} compiler version")
    compiler_invocation = exact_keys(build["compiler_invocation"], {
        "invocation", "resolved_path"}, f"{role} compiler invocation")
    require(compiler_invocation["invocation"] == cache["CMAKE_CXX_COMPILER"] and
            compiler_invocation["resolved_path"] == compiler["path"],
            f"{role} compiler invocation differs")
    external = validate_external_link_inputs(
        build["validated_external_link_inputs"], f"{role} executable recipe")
    require(all(executable["mtime_ns"] >= item["artifact"]["mtime_ns"]
                for item in external),
            f"{role} executable predates an external input")

    if selected:
        tool = validate_artifact(
            build["multi_config_build_tool"], f"{role} build tool", "build_tool")
        version = validate_version_text(
            build["multi_config_build_tool_version_stdout"],
            f"{role} build-tool version")
        require(tool["path"] == CANONICAL_NINJA_PATH and version is not None,
                f"{role} Ninja identity differs")
        graph = validate_ninja_graph(
            build["multi_config_ninja_graph"], build_dir, selected,
            f"{role} Ninja graph")
        link_record = next((record for record in graph["files"]
                            if record["relative_path"] ==
                            f"CMakeFiles/impl-{selected}.ninja"), None)
        require(type(link_record) is dict and all(
                    link_record["artifact"][name] ==
                        build["archive_link_recipe"][name] ==
                        build["executable_link_recipe"][name]
                    for name in ("path", "size", "mode", "mtime_ns", "sha256")) and
                all(archive["mtime_ns"] >= item["artifact"]["mtime_ns"] and
                    executable["mtime_ns"] >= item["artifact"]["mtime_ns"]
                    for item in graph["files"]),
                f"{role} retained Ninja graph/link freshness differs")
    else:
        require(build["multi_config_build_tool"] is None and
                build["multi_config_build_tool_version_stdout"] is None and
                build["multi_config_ninja_graph"] is None,
                f"{role} single-config build retained Ninja metadata")

    compile_record = exact_keys(
        build["validated_compile_commands"], COMPILE_RECORD_KEYS,
        f"{role} compile-command record")
    pairs = compile_record["required_source_object_pairs"]
    sources = compile_record["required_sources"]
    entries = compile_record["required_entries"]
    require(compile_record["schema"] == COMPILE_COMMANDS_SCHEMA and
            compile_record["implementation"] == role and
            compile_record["profile"] == profile and
            compile_record["isa_policy"] == isa_policy and
            type(compile_record["entry_count"]) is int and
            compile_record["entry_count"] == base_entry_count * configuration_count and
            compile_record["validated_optimization"] == "-O3" and
            compile_record["validated_openmp"] is True and
            type(sources) is list and sources == sorted(set(sources)) and sources and
            type(pairs) is list and type(entries) is list and
            len(pairs) == len(entries) and pairs,
            f"{role} compile-command contract differs")
    validated_pairs: list[dict[str, Any]] = []
    for index, pair in enumerate(pairs):
        validated_pairs.append(validate_source_object_pair(
            pair, f"{role} source/object pair {index}", build_dir))
    expected_sources = sorted([
        *(str(Path(specification[f"{role}_source_root"]) / name)
          for name in library_sources),
        benchmark_source,
    ])
    require([pair["source"]["path"] for pair in validated_pairs] ==
                expected_sources == sources and
            len({pair["object"]["path"] for pair in validated_pairs}) == len(pairs),
            f"{role} source/object closure differs")
    by_source = {pair["source"]["path"]: pair for pair in validated_pairs}
    benchmark_pair = by_source.get(benchmark_source)
    require(type(benchmark_pair) is dict,
            f"{role} benchmark object is not unique")
    if role == "candidate":
        require(candidate_tree is not None,
                "candidate source tree is absent")
        validate_generated_attestation(
            compile_record["generated_attestation_header"], build_dir,
            benchmark_pair["object"], specification["candidate_commit"],
            candidate_tree)
        configuration = validate_effective_configuration(
            compile_record["effective_build_configuration"], build_dir,
            specification["candidate_source_root"], benchmark_pair["object"], cache)
        require(cache["LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"] ==
                    configuration["configuration_sha256"],
                "candidate cache/configuration digest differs")
    else:
        require(compile_record["generated_attestation_header"] is None and
                compile_record["effective_build_configuration"] is None,
                "baseline unexpectedly retains candidate generated inputs")
        configuration = None
    require(type(entries) is list and
            [entry.get("file") if type(entry) is dict else None
             for entry in entries] == expected_sources,
            f"{role} retained compiler entry order differs")
    for index, entry in enumerate(entries):
        entry = exact_keys(entry, {"directory", "file", "output", "arguments"},
                           f"{role} compile entry {index}")
        pair = by_source[entry["file"]]
        output = expected_compile_output(role, entry["file"], specification, selected)
        arguments = expected_compile_argv(
            role, entry["file"], specification, cache["CMAKE_CXX_COMPILER"],
            selected, configuration, compiler["path"])
        require(entry == {
                    "directory": build_dir, "file": entry["file"],
                    "output": output, "arguments": arguments} and
                pair["object"]["path"] == str(Path(build_dir) / output),
                f"{role} compile entry {index} argv/output differs")

    archive_pairs = [pair for pair in validated_pairs if pair is not benchmark_pair]
    require(archive_pairs and
            all(Path(pair["source"]["path"]).is_relative_to(
                    Path(specification[f"{role}_source_root"]))
                for pair in archive_pairs) and
            archive["mtime_ns"] >= max(pair["object"]["mtime_ns"]
                                        for pair in archive_pairs) and
            executable["mtime_ns"] >= archive["mtime_ns"] and
            executable["mtime_ns"] >= benchmark_pair["object"]["mtime_ns"],
            f"{role} output freshness/source-root closure differs")
    source_order: Sequence[str] = library_sources
    if selected and role == "candidate":
        object_order = tuple(name for name in (
            "Leopard2BackendSSSE3.cpp", "Leopard2BackendAVX2.cpp",
            "Leopard2BackendAVX2Xor.cpp", "Leopard2BackendAVX2T16Q2.cpp",
            "Leopard2BackendAVX2T8K62.cpp", "Leopard2BackendAVX2T16K66.cpp",
            "Leopard2BackendAVX2T32B256.cpp", "Leopard2BackendAVX2T16B64.cpp",
            "Leopard2LowP32B64AVX2.cpp", "Leopard2BackendAVX2T2K4.cpp",
            "Leopard2BackendAVX2T8K8B1024.cpp", "Leopard2BackendGFNI.cpp",
        ) if name in library_sources)
        source_order = (*object_order,
                        *(name for name in library_sources
                          if name not in object_order))
    members = [Path(name).name + ".o" for name in source_order]
    require(build["validated_archive_members"] == members and
            len(members) == len(set(members)),
            f"{role} archive member closure differs")
    objects_by_member = {
        Path(pair["object"]["path"]).name:
            Path(pair["object"]["path"]).relative_to(Path(build_dir)).as_posix()
        for pair in archive_pairs
    }
    require(set(objects_by_member) == set(members),
            f"{role} archive members differ from objects")
    prefix = f"{selected}/" if selected else ""
    archive_operand = prefix + archive_name
    executable_operand = prefix + executable_name
    validate_archive_recipe(
        build["archive_link_recipe_content"],
        None if selected else build["archive_link_recipe"],
        archive_operand, archive_target,
        [objects_by_member[name] for name in members],
        tools["archiver"]["invocation"], tools["ranlib"]["invocation"],
        f"{role} archive link recipe")
    benchmark_object = Path(benchmark_pair["object"]["path"]).relative_to(
        Path(build_dir)).as_posix()
    validate_executable_recipe(
        build["executable_link_recipe_content"],
        None if selected else build["executable_link_recipe"],
        cache["CMAKE_CXX_COMPILER"], archive_operand, executable_operand,
        benchmark_object, external, f"{role} executable link recipe")
    return build


def validate_timestamp(value: Any, label: str) -> None:
    require(type(value) is str and UTC.fullmatch(value) is not None,
            f"{label} is not a canonical UTC timestamp")
    try:
        parsed = dt.datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise AuditError(f"{label} is not a valid timestamp: {error}") from error
    require(parsed.tzinfo == dt.timezone.utc and
            parsed.isoformat().replace("+00:00", "Z") == value,
            f"{label} is not canonical UTC")


def metadata_tuple(value: os.stat_result) -> tuple[int, ...]:
    return (value.st_mode, value.st_nlink, value.st_size, value.st_mtime_ns,
            value.st_ctime_ns, value.st_dev, value.st_ino)


def canonical_existing_file(path: Path) -> Path:
    require(path.is_absolute(), f"path is not absolute: {path}")
    resolved = path.resolve(strict=True)
    require(resolved == path and stat.S_ISREG(path.lstat().st_mode) and
            not path.is_symlink(), f"path is aliased or not regular: {path}")
    return path


def stable_file(path: Path, maximum: int) -> tuple[bytes, dict[str, Any]]:
    canonical_existing_file(path)
    before = path.lstat()
    require(before.st_nlink == 1 and 0 <= before.st_size <= maximum,
            f"unsafe link count or size: {path}")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    descriptor = os.open(path, flags)
    try:
        opened = os.fstat(descriptor)
        require(metadata_tuple(opened) == metadata_tuple(before),
                f"file changed before read: {path}")
        blocks: list[bytes] = []
        total = 0
        while True:
            block = os.read(descriptor, min(1024 * 1024, maximum + 1 - total))
            if not block:
                break
            blocks.append(block)
            total += len(block)
            require(total <= maximum, f"file exceeds its bound: {path}")
        closed = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after = path.lstat()
    require(metadata_tuple(before) == metadata_tuple(closed) ==
            metadata_tuple(after) and total == before.st_size,
            f"file changed while reading: {path}")
    data = b"".join(blocks)
    return data, {"sha256": sha256_bytes(data), "size": len(data)}


def safe_retained_file(root: Path, relative: Any, maximum: int,
                       label: str) -> tuple[bytes, dict[str, Any]]:
    require(type(relative) is str and relative and not os.path.isabs(relative),
            f"{label} path is not one relative path")
    parts = Path(relative).parts
    require(parts and all(part not in ("", ".", "..") for part in parts),
            f"{label} path traverses outside the evidence root")
    require(Path(relative).as_posix() == relative,
            f"{label} path is not canonical POSIX-relative text")
    path = root.joinpath(*parts)
    require(path.resolve(strict=True).is_relative_to(root) and
            path.resolve(strict=True) == path,
            f"{label} path escaped or is aliased")
    return stable_file(path, maximum)


def verify_signature(value: Any, label: str) -> dict[str, Any]:
    require(type(value) is dict, f"{label} is not an object")
    digest = value.get("digest")
    require(type(digest) is str and HEX256.fullmatch(digest) is not None,
            f"{label} has no canonical digest")
    payload = {key: item for key, item in value.items() if key != "digest"}
    require(sha256_bytes(canonical_bytes(payload)) == digest,
            f"{label} canonical digest differs")
    return value


def finite_positive(value: Any, label: str) -> float:
    require(type(value) in (int, float), f"{label} is not numeric")
    result = float(value)
    require(math.isfinite(result) and result > 0,
            f"{label} is not finite and positive")
    return result


def close_enough(actual: float, expected: float) -> bool:
    return abs(actual - expected) <= max(0.000002, abs(expected) * 0.000002)


def validate_summary(value: Any, label: str, *, setup: bool = False) \
        -> list[float]:
    require(type(value) is dict, f"{label} is not a timing summary")
    suffix = "" if setup else "_per_batch_call"
    sample_name = "samples_us" if setup else "samples_us_per_batch_call"
    names = {
        "median": f"median_us{suffix}",
        "mad": f"mad_us{suffix}",
        "minimum": f"minimum_us{suffix}",
        "maximum": f"maximum_us{suffix}",
    }
    samples_value = value.get(sample_name)
    require(type(samples_value) is list and len(samples_value) == ITERATIONS,
            f"{label} does not retain exactly {ITERATIONS} samples")
    samples = [finite_positive(item, f"{label} sample")
               for item in samples_value]
    observed_median = statistics.median(samples)
    deviations = [abs(item - observed_median) for item in samples]
    derived = {
        names["median"]: observed_median,
        names["mad"]: statistics.median(deviations),
        names["minimum"]: min(samples),
        names["maximum"]: max(samples),
    }
    for name, expected in derived.items():
        actual = value.get(name)
        require(type(actual) in (int, float) and not isinstance(actual, bool) and
                math.isfinite(float(actual)) and
                (float(actual) >= 0 if name == names["mad"]
                 else float(actual) > 0) and
                close_enough(float(actual), expected),
                f"{label} {name} is not raw-derived")
    return samples


class XorShift64:
    def __init__(self, seed: int):
        self.state = seed if seed else 0x9E3779B97F4A7C15

    def next(self) -> int:
        value = self.state
        value ^= (value << 13) & MASK64
        value ^= value >> 7
        value ^= (value << 17) & MASK64
        self.state = value & MASK64
        return self.state


def expected_losses() -> list[int]:
    order = list(range(CELL["k"]))
    random = XorShift64(CELL["seed"] ^ 0xD1B54A32D192ED03)
    for remaining in range(len(order), 1, -1):
        selected = random.next() % remaining
        order[remaining - 1], order[selected] = \
            order[selected], order[remaining - 1]
    return sorted(order[:CELL["losses"]])


def validate_digests(value: Any, label: str) -> dict[str, str]:
    value = exact_keys(value, {
        "algorithm", "original_data", "transmitted_parity",
        "recovered_originals"}, label)
    require(value["algorithm"] == "fnv1a64", f"{label} algorithm changed")
    for name in ("original_data", "transmitted_parity", "recovered_originals"):
        require(type(value[name]) is str and HEX64.fullmatch(value[name]) is not None,
                f"{label} {name} is not lowercase FNV-1a hex")
    return {name: value[name] for name in (
        "original_data", "transmitted_parity", "recovered_originals")}


def validate_wall_duration(duration_ns: Any, normalized: Mapping[str, Any],
                           role: str, label: str) -> None:
    require(type(duration_ns) is int and duration_ns > 0 and
            role in {"baseline", "candidate"},
            f"{label} wall duration is invalid")
    metric_totals = [sum(normalized["encode"]) * REUSE]
    if role == "candidate":
        metric_totals.extend((
            sum(normalized["one_shot_encode"]) * REUSE,
            sum(normalized["codec_setup"]),
        ))
    require(duration_ns / 1000.0 >= sum(metric_totals),
            f"{label} wall duration cannot enclose validated encode/setup "
            "timer work")


def expected_parameters(role: str) -> dict[str, Any]:
    common = {
        "K": CELL["k"], "R": CELL["r"],
        "shard_bytes": CELL["shard_bytes"],
        "loss_count": CELL["losses"],
        "missing_original_indices": expected_losses(),
        "batch": 1, "reuse": REUSE, "iterations": ITERATIONS,
        "warmup": WARMUP, "thread_count": 1, "seed": CELL["seed"],
    }
    if role == "baseline":
        return {**common, "logical_shard_bytes": CELL["shard_bytes"]}
    return {
        "K": CELL["k"], "R": CELL["r"],
        "requested_profile": "legacy_high_v1",
        "requested_field": "gf16", "requested_backend": "auto",
        "force_generic_decode": False, "force_specialized_decode": False,
        "force_tiled_decode": False, "force_materialized_decode": False,
        "skip_legacy": True, "retain_samples": True,
        "attest_source": True, "measure_one_shot_encode": True,
        "auto_gf16_gfni_encode_mode": 1,
        **{key: common[key] for key in (
            "shard_bytes", "loss_count", "missing_original_indices", "batch",
            "reuse", "iterations", "warmup", "thread_count", "seed")},
    }


def validate_child_result(role: str, result: Any, campaign: Mapping[str, Any],
                          identities: Mapping[str, Any]) -> dict[str, Any]:
    require(role in ("baseline", "candidate") and type(result) is dict,
            f"{role} result is not an object")
    require(result.get("schema") ==
            (BASELINE_SCHEMA if role == "baseline" else CANDIDATE_SCHEMA),
            f"{role} child schema differs")
    parameters = result.get("parameters")
    expected = expected_parameters(role)
    require(type(parameters) is dict and set(parameters) == set(expected) and
            exact_json_equal(parameters, expected),
            f"{role} child parameters differ from frozen v17")
    resolved = result.get("resolved")
    common_resolved = {
        "profile": "legacy_high_v1", "field": "gf16",
        "parent_count": 2048, "padded_side": 256, "thread_count": 1,
    }
    require(type(resolved) is dict and all(
                type(resolved.get(name)) is type(item) and resolved[name] == item
                for name, item in common_resolved.items()),
            f"{role} resolved geometry differs")
    metrics = result.get("metrics")
    require(type(metrics) is dict, f"{role} metrics are absent")
    encode = validate_summary(
        metrics.get("encode_execution"), f"{role} encode")
    digests = validate_digests(
        result.get("workload_digests"), f"{role} workload digests")
    if role == "baseline":
        require(set(resolved) == set(common_resolved) | {
                    "padded_application_bytes", "padding_policy"} and
                resolved["padded_application_bytes"] is False and
                resolved["padding_policy"] == "zero suffix per shard",
                "baseline padding identity differs")
        build = result.get("build")
        require(type(build) is dict and
                build.get("main_source_commit") == MAIN_COMMIT and
                build.get("pure_avx2") is False and
                type(build.get("cplusplus")) is int,
                "baseline child is not the native exact-main build")
        correctness = result.get("correctness")
        require(type(correctness) is dict and
                correctness.get("round_trip") is True and
                correctness.get("logical_prefix_fingerprinted") is True,
                "baseline correctness identity differs")
        require(metrics.get("codec_setup") is None and
                metrics.get("decode_timing_includes_setup") is True,
                "baseline setup semantics differ")
        window = min(encode) * REUSE
        require(window >= TIMER_FLOOR_US,
                "baseline retained encode timer window is below 20 ms")
        return {
            "digests": digests,
            "backend": "exact_main_native",
            "encode": encode,
            "retained_timer_window_us": {"encode": window},
        }

    require(set(resolved) == set(common_resolved) | {"backend", "encode_backend"} and
            resolved["backend"] == "avx2" and
            resolved["encode_backend"] == "avx2-gfni",
            "candidate context/encode backend identity differs")
    require("measure_one_shot_decode" not in parameters and
            "one_shot_decode_including_setup" not in metrics,
            "candidate v17 contains a forbidden one-shot decode measurement")
    route = {
        "auto_gf16_gfni_encode_diagnostic_mode": 1,
        "auto_gf16_gfni_encode_diagnostic_disabled": False,
        "auto_gf16_gfni_encode_mode_latched": 1,
        "auto_gf16_gfni_encode_kernel_available": True,
        "auto_gf16_gfni_encode_kernel_qualified": True,
        "auto_gf16_gfni_encode_selector_expected_selected": True,
        "auto_gf16_gfni_encode_selector_selected": True,
        "auto_gf16_gfni_encode_observed_call_count": 2,
        "auto_gf16_gfni_encode_selector_contract": ROUTE_CONTRACT,
        "auto_gf16_gfni_encode_timed_ordinary_encode_api": ORDINARY_ENCODE_API,
        "auto_gf16_gfni_encode_timed_one_shot_encode_api": ONE_SHOT_ENCODE_API,
    }
    build = result.get("build")
    require(type(build) is dict, "candidate build identity is absent")
    for name, expected_value in route.items():
        actual = build.get(name)
        require(type(actual) is type(expected_value) and actual == expected_value,
                f"candidate route field {name} differs")
    candidate_source = identities.get("candidate_source")
    require(type(candidate_source) is dict and
            build.get("source_commit") == candidate_source.get("head") and
            build.get("source_tree") == candidate_source.get("tree") and
            build.get("source_tracked_dirty") is False,
            "candidate child source attestation differs from retained source")
    correctness = result.get("correctness")
    require(type(correctness) is dict and
            correctness.get("leopard2_round_trip") is True,
            "candidate round trip failed")
    legacy = result.get("legacy")
    require(type(legacy) is dict and legacy.get("available") is False and
            legacy.get("unavailable_reason") == "disabled by --skip-legacy",
            "candidate silently used the in-tree legacy path")
    codec_setup = validate_summary(
        metrics.get("codec_setup"), "candidate codec setup", setup=True)
    one_shot = validate_summary(
        metrics.get("one_shot_encode"), "candidate one-shot encode")
    encode_window = min(encode) * REUSE
    one_shot_window = min(one_shot) * REUSE
    require(encode_window >= TIMER_FLOOR_US and
            one_shot_window >= TIMER_FLOOR_US,
            "candidate retained encode timer window is below 20 ms")
    return {
        "digests": digests,
        "backend": "avx2-gfni",
        "encode": encode,
        "one_shot_encode": one_shot,
        "codec_setup": codec_setup,
        "retained_timer_window_us": {
            "encode": encode_window,
            "one_shot_encode": one_shot_window,
        },
    }


def expected_command(role: str, taskset: str) -> list[str]:
    command = [
        taskset, "-c", str(CPU), SEALED_EXECUTABLE_COMMAND[role],
        "--k", str(CELL["k"]), "--r", str(CELL["r"]),
        "--bytes", str(CELL["shard_bytes"]),
        "--loss", str(CELL["losses"]), "--batch", "1",
        "--reuse", str(REUSE), "--iterations", str(ITERATIONS),
        "--warmup", str(WARMUP), "--threads", "1",
        "--seed", str(CELL["seed"]),
    ]
    if role == "candidate":
        command.extend((
            "--profile", "high", "--field", "gf16", "--backend", "auto",
            "--skip-legacy", "--retain-samples", "--measure-one-shot-encode",
            "--attest-source", "--auto-gf16-gfni-encode-mode", "1",
        ))
    command.extend(("--json", "-"))
    return command


def validate_acquisition_paths(input_spec: Mapping[str, Any]) -> None:
    expected_runner = str(
        Path(input_spec["candidate_source_root"]) /
        "experiments/leopard2/main_compare/run_abba.py")
    require(input_spec["runner"] == expected_runner and
            input_spec["taskset"] == CANONICAL_TASKSET_PATH and
            input_spec["ldd"] == CANONICAL_LDD_PATH,
            "v17 runner/taskset/ldd acquisition paths differ")


def validate_input_specification(input_spec: Any) -> dict[str, Any]:
    input_spec = exact_keys(input_spec, INPUT_KEYS, "input specification")
    require(all(type(input_spec[name]) is str and input_spec[name]
                for name in INPUT_KEYS - {"baseline_pure_avx2"}) and
            input_spec["baseline_pure_avx2"] is False and
            HEX40.fullmatch(input_spec["candidate_commit"]) is not None,
            "input specification does not bind the native v17 profiles")
    require(all(is_canonical_absolute_path(input_spec[name])
                for name in INPUT_PATH_KEYS),
            "input specification contains a noncanonical path")
    validate_acquisition_paths(input_spec)
    return input_spec


def validate_native_build(input_spec: Any, identities: Any) -> dict[str, Any]:
    input_spec = validate_input_specification(input_spec)
    identities = exact_keys(identities, INPUT_SNAPSHOT_KEYS,
                            "complete input snapshot")
    for name, kind in (("runner", "file"), ("taskset", "executable"),
                       ("ldd", "executable")):
        artifact = validate_artifact(identities[name], name, kind)
        require(artifact["path"] == input_spec[name],
                f"{name} identity differs from input specification")
    helper = validate_artifact(
        identities["evidence_helper"], "evidence helper", "file")
    require(helper["path"] == str(Path(input_spec["candidate_source_root"]) /
                                    EVIDENCE_HELPER_RELATIVE_PATH),
            "evidence-helper identity differs from candidate source")
    baseline_source = validate_git_capture(
        identities["baseline_source"], input_spec["baseline_source_root"],
        MAIN_COMMIT, True, "baseline")
    candidate_source = validate_git_capture(
        identities["candidate_source"], input_spec["candidate_source_root"],
        input_spec["candidate_commit"], False, "candidate")
    require(exact_json_equal(baseline_source["git_executable"],
                             candidate_source["git_executable"]),
            "baseline/candidate source captures used different Git executables")
    selector_records = [
        record for record in candidate_source["tracked_files"]
        if record.get("path") == "leopard2.cpp"]
    require(len(selector_records) == 1 and
            selector_records[0].get("kind") == "regular" and
            HEX40.fullmatch(str(selector_records[0].get("object_id"))) is not None,
            "candidate production selector source is not unique and regular")

    baseline_build = validate_build(
        identities["baseline_build"], "baseline", input_spec, None)
    candidate_build = validate_build(
        identities["candidate_build"], "candidate", input_spec,
        candidate_source["tree"])
    for role, build in (("baseline", baseline_build),
                        ("candidate", candidate_build)):
        executable = validate_artifact(
            identities[f"{role}_executable"], f"{role} top-level executable",
            "executable")
        archive = validate_artifact(
            identities[f"{role}_archive"], f"{role} top-level archive", "archive")
        require(executable["path"] == input_spec[f"{role}_executable"] and
                archive["path"] == input_spec[f"{role}_archive"] and
                exact_json_equal(executable, build["validated_executable"]) and
                exact_json_equal(archive, build["validated_archive"]),
                f"{role} top-level/build output identity differs")
        validate_runtime_closure(
            identities[f"{role}_runtime_closure"], f"{role} runtime closure",
            executable["path"])
    require(exact_json_equal(baseline_build["compiler"],
                             candidate_build["compiler"]) and
            exact_json_equal(baseline_build["compiler_version_stdout"],
                             candidate_build["compiler_version_stdout"]) and
            exact_json_equal(baseline_build["multi_config_build_tool"],
                             candidate_build["multi_config_build_tool"]) and
            exact_json_equal(
                baseline_build["multi_config_build_tool_version_stdout"],
                candidate_build["multi_config_build_tool_version_stdout"]) and
            exact_json_equal(baseline_build["archiver"],
                             candidate_build["archiver"]) and
            exact_json_equal(baseline_build["ranlib"],
                             candidate_build["ranlib"]) and
            exact_json_equal(
                baseline_build["archive_link_tool_invocations"],
                candidate_build["archive_link_tool_invocations"]) and
            exact_json_equal(baseline_build["validated_external_link_inputs"],
                             candidate_build["validated_external_link_inputs"]),
            "baseline/candidate compiler, archive tools, build tool, or "
            "external link inputs differ")
    attestation = candidate_build["validated_compile_commands"][
        "generated_attestation_header"]
    require(attestation["source_commit"] == candidate_source["head"] and
            attestation["source_tree"] == candidate_source["tree"] and
            attestation["source_tracked_dirty"] is False,
            "candidate generated attestation differs from Git provenance")
    return identities


def validate_cpu_policy(record: Any, label: str, cpu: int,
                        sibling: int) -> dict[str, Any]:
    record = exact_keys(record, {
        "cpu", "online", "cpuinfo", "topology", "frequency_policy",
        "cache_hierarchy", "cache_index_inventory",
        "cache_directory_inventory_text", "numa_nodes",
        "numa_node_inventory", "core_class"}, f"host {label}")
    require(type(record["cpu"]) is int and record["cpu"] == cpu and
            (record["online"] is None or type(record["online"]) is str),
            f"host {label} CPU identity differs")
    cpuinfo = record["cpuinfo"]
    allowed_cpuinfo = {
        "processor", "vendor_id", "cpu family", "model", "model name",
        "stepping", "microcode", "flags", "Features", "CPU implementer",
        "CPU architecture", "CPU variant", "CPU part", "CPU revision",
    }
    require(type(cpuinfo) is dict and cpuinfo and
            set(cpuinfo).issubset(allowed_cpuinfo) and
            all(type(item) is str for item in cpuinfo.values()) and
            cpuinfo.get("processor") == str(cpu) and
            cpuinfo.get("vendor_id") == "AuthenticAMD" and
            cpuinfo.get("cpu family") == "26" and
            cpuinfo.get("model") == "8" and
            type(cpuinfo.get("flags")) is str and
            {"avx2", "gfni"}.issubset(set(cpuinfo["flags"].split())) and
            any(name in cpuinfo for name in ("model name", "CPU part")),
            f"host {label} is not the calibrated AMD 1Ah/08h CPU")
    topology = exact_keys(record["topology"], {
        "core_id", "physical_package_id", "die_id", "cluster_id",
        "thread_siblings_list", "core_siblings_list"},
        f"host {label} topology")
    require(all(item is None or type(item) is str for item in topology.values()) and
            type(topology["thread_siblings_list"]) is str and
            parse_cpu_list(topology["thread_siblings_list"],
                           f"host {label} SMT siblings") == {cpu, sibling} and
            type(topology["core_siblings_list"]) is str and
            {cpu, sibling}.issubset(parse_cpu_list(
                topology["core_siblings_list"], f"host {label} core siblings")),
            f"host {label} topology differs from reserved SMT pair")
    frequency = exact_keys(record["frequency_policy"], {
        "scaling_driver", "scaling_governor", "energy_performance_preference",
        "scaling_min_freq", "scaling_max_freq", "cpuinfo_min_freq",
        "cpuinfo_max_freq"}, f"host {label} frequency policy")
    require(all(item is None or type(item) is str for item in frequency.values()),
            f"host {label} frequency policy differs")
    caches = record["cache_hierarchy"]
    cache_keys = {
        "index", "level", "type", "size", "coherency_line_size",
        "number_of_sets", "ways_of_associativity", "physical_line_partition",
        "shared_cpu_list", "shared_cpu_map", "allocation_policy", "write_policy",
    }
    require(type(caches) is list and caches,
            f"host {label} cache hierarchy is absent")
    indices: list[int] = []
    for cache in caches:
        cache = exact_keys(cache, cache_keys, f"host {label} cache")
        require(type(cache["index"]) is int and cache["index"] >= 0 and
                all(cache[name] is None or type(cache[name]) is str
                    for name in cache_keys - {"index"}) and
                all(type(cache[name]) is str and cache[name]
                    for name in ("level", "type", "size", "coherency_line_size",
                                 "shared_cpu_list", "shared_cpu_map")) and
                cpu in parse_cpu_list(cache["shared_cpu_list"],
                                      f"host {label} cache sharing"),
                f"host {label} cache identity differs")
        indices.append(cache["index"])
    require(indices == sorted(set(indices)),
            f"host {label} cache indices differ")
    inventory = record["cache_index_inventory"]
    inventory_text = validate_text(
        record["cache_directory_inventory_text"],
        f"host {label} cache inventory text")
    expected_inventory = [f"index{index}" for index in indices]
    require(inventory == expected_inventory and
            inventory_text["text"] ==
                "".join(name + "\n" for name in expected_inventory),
            f"host {label} cache inventory differs")
    nodes = record["numa_nodes"]
    require(type(nodes) is list and nodes == sorted(set(nodes)) and nodes and
            all(type(node) is int and node >= 0 for node in nodes) and
            record["numa_node_inventory"] == [f"node{node}" for node in nodes],
            f"host {label} NUMA identity differs")
    core_class = exact_keys(record["core_class"], {
        "core_type", "cpu_capacity"}, f"host {label} core class")
    require(all(item is None or type(item) is str for item in core_class.values()),
            f"host {label} core-class identity differs")
    return record


def validate_host(value: Any, label: str,
                  allowed: Sequence[int]) -> dict[str, Any]:
    value = exact_keys(value, {
        "system", "allowed_cpu_set_at_launch", "online_cpu_set",
        "online_cpu_list_text", "online_node_list_text", "benchmark_cpu",
        "reserved_sibling", "turbo_and_pstate"}, f"{label} host")
    system = exact_keys(value["system"], {
        "system", "node", "release", "version", "machine", "python",
        "page_size"}, f"{label} host system")
    require(system["system"] == "Linux" and system["machine"] == "x86_64" and
            all(type(system[name]) is str and system[name]
                for name in ("node", "release", "version", "python")) and
            type(system["page_size"]) is int and system["page_size"] > 0,
            f"{label} host system differs")
    launch = value["allowed_cpu_set_at_launch"]
    online = value["online_cpu_set"]
    require(type(launch) is list and launch == list(allowed) and
            launch == sorted(set(launch)) and
            3 <= len(launch) <= MAX_CPU_LIST_ENTRIES and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in launch) and
            type(online) is list and online == sorted(set(online)) and
            1 <= len(online) <= MAX_CPU_LIST_ENTRIES and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in online) and
            CPU in online and SIBLING in online,
            f"{label} host launch/online CPU sets differ")
    online_text = validate_text(
        value["online_cpu_list_text"], f"{label} online CPU list")
    node_text = validate_text(
        value["online_node_list_text"], f"{label} online NUMA-node list")
    online_nodes = parse_cpu_list(node_text["text"],
                                  f"{label} online NUMA nodes")
    require(parse_cpu_list(online_text["text"], f"{label} online CPUs") ==
                set(online) and online_nodes,
            f"{label} online CPU/node text differs")
    benchmark = validate_cpu_policy(
        value["benchmark_cpu"], f"{label} benchmark", CPU, SIBLING)
    sibling = validate_cpu_policy(
        value["reserved_sibling"], f"{label} sibling", SIBLING, CPU)
    require(exact_json_equal(benchmark["cache_hierarchy"],
                             sibling["cache_hierarchy"]) and
            benchmark["numa_nodes"] == sibling["numa_nodes"] and
            set(benchmark["numa_nodes"]).issubset(online_nodes) and
            exact_json_equal(benchmark["core_class"], sibling["core_class"]),
            f"{label} reserved pair cache/NUMA/core-class differs")
    turbo = exact_keys(value["turbo_and_pstate"], {
        "/sys/devices/system/cpu/intel_pstate/no_turbo",
        "/sys/devices/system/cpu/amd_pstate/status",
        "/sys/devices/system/cpu/cpufreq/boost"},
        f"{label} turbo/pstate")
    require(all(item is None or type(item) is str for item in turbo.values()),
            f"{label} turbo/pstate identity differs")
    return value


def validate_sealed_snapshot(snapshot: Any, source: Mapping[str, Any],
                             label: str) -> dict[str, Any]:
    snapshot = exact_keys(snapshot, {
        "protocol", "size", "mode", "sha256", "seals", "elf"}, label)
    require(snapshot["protocol"] == SEALED_EXECUTABLE_PROTOCOL and
            type(snapshot["size"]) is int and snapshot["size"] > 0 and
            type(snapshot["mode"]) is int and
            0 <= snapshot["mode"] <= 0o7777 and
            bool(snapshot["mode"] & 0o111) and
            type(snapshot["sha256"]) is str and
            HEX256.fullmatch(snapshot["sha256"]) is not None and
            type(snapshot["seals"]) is int and
            0 <= snapshot["seals"] <= MAX_LINUX_SEAL_MASK and
            snapshot["seals"] & REQUIRED_EXECUTABLE_SEALS ==
                REQUIRED_EXECUTABLE_SEALS and
            snapshot["elf"] is True and
            snapshot["size"] == source.get("size") and
            snapshot["mode"] == source.get("mode") and
            snapshot["sha256"] == source.get("sha256"),
            f"{label} differs from source bytes")
    return snapshot


def validate_sealed_executables(value: Any, identities: Mapping[str, Any]) -> None:
    value = exact_keys(value, {"baseline", "candidate"},
                       "sealed executable set")
    for role in ("baseline", "candidate"):
        record = exact_keys(value[role], {"source", "snapshot", "runtime_closure"},
                            f"{role} sealed executable")
        source = validate_artifact(
            record["source"], f"{role} sealed executable source", "executable")
        validate_sealed_snapshot(
            record["snapshot"], source, f"{role} sealed snapshot")
        require(exact_json_equal(source, identities.get(f"{role}_executable")) and
                type(identities.get(f"{role}_build")) is dict and
                exact_json_equal(source,
                                 identities[f"{role}_build"].get(
                                     "validated_executable")),
                f"{role} sealed source differs from build identity")
        closure = validate_runtime_closure(
            record["runtime_closure"], f"{role} sealed executable",
            SEALED_EXECUTABLE_COMMAND[role])
        source_closure = identities.get(f"{role}_runtime_closure")
        require(type(source_closure) is dict and
                exact_json_equal(closure["dependencies"],
                                 source_closure["dependencies"]) and
                exact_json_equal(closure["canonical_ldd_output"],
                                 source_closure["canonical_ldd_output"]),
                f"{role} sealed runtime closure differs")


def cpu_stat(value: Any, cpu: int, label: str) -> dict[str, Any]:
    value = exact_keys(value, {"cpu", "fields", "total_jiffies"}, label)
    fields = exact_keys(value["fields"], set(CPU_STAT_FIELDS), f"{label} fields")
    require(type(value["cpu"]) is int and value["cpu"] == cpu and
            all(type(fields[name]) is int and fields[name] >= 0
                for name in CPU_STAT_FIELDS) and
            type(value["total_jiffies"]) is int and
            value["total_jiffies"] == sum(fields.values()),
            f"{label} CPU counters differ")
    return value


def cpu_delta(before: Any, after: Any, cpu: int, label: str) -> dict[str, Any]:
    before = cpu_stat(before, cpu, f"{label} before")
    after = cpu_stat(after, cpu, f"{label} after")
    fields: dict[str, int] = {}
    for name in CPU_STAT_FIELDS:
        require(after["fields"][name] >= before["fields"][name],
                f"{label} {name} moved backwards")
        fields[name] = after["fields"][name] - before["fields"][name]
    idle = sum(fields[name] for name in CPU_IDLE_FIELDS)
    nonidle = sum(fields[name] for name in CPU_STAT_FIELDS
                  if name not in CPU_IDLE_FIELDS)
    return {
        "cpu": cpu, "fields": fields, "idle_jiffies": idle,
        "nonidle_jiffies": nonidle, "total_jiffies": idle + nonidle,
    }


def validate_reservation(value: Any) -> dict[str, Any]:
    value = exact_keys(value, {"path", "sha256", "payload", "lock"},
                       "CPU reservation")
    require(is_canonical_absolute_path(value["path"]) and
            type(value["sha256"]) is str and
            HEX256.fullmatch(value["sha256"]) is not None and
            value["lock"] == "exclusive_nonblocking",
            "CPU reservation identity differs")
    payload = exact_keys(value["payload"], {
        "benchmark_cpu", "nonce", "owner", "reserved_sibling", "schema",
        "status"}, "CPU reservation payload")
    require(payload["schema"] == RESERVATION_SCHEMA and
            payload["status"] == "held" and
            type(payload["benchmark_cpu"]) is int and
            payload["benchmark_cpu"] == CPU and
            type(payload["reserved_sibling"]) is int and
            payload["reserved_sibling"] == SIBLING and
            type(payload["owner"]) is str and payload["owner"].strip() and
            type(payload["nonce"]) is str and
            8 <= len(payload["nonce"]) <= 256 and
            value["sha256"] == sha256_bytes(canonical_bytes(payload)),
            "CPU reservation payload/digest differs")
    return value


def validate_pair_lease(value: Any) -> dict[str, Any]:
    value = exact_keys(value, {
        "device", "directory_device", "directory_inode", "inode", "lock",
        "path", "payload", "sha256"}, "CPU pair lease")
    payload = exact_keys(value["payload"], {"cpus", "schema", "uid"},
                         "CPU pair lease payload")
    require(payload["schema"] == PAIR_LEASE_SCHEMA and
            exact_json_equal(payload["cpus"], [CPU, SIBLING]) and
            all(type(cpu) is int for cpu in payload["cpus"]) and
            type(payload["uid"]) is int and payload["uid"] >= 0,
            "CPU pair lease payload differs")
    expected_path = (
        Path("/run/user") / str(payload["uid"]) / "leopard2-cpu-leases" /
        f"leopard2-cpu-pair-{payload['uid']}-{CPU}-{SIBLING}.lock")
    require(value["path"] == str(expected_path) and
            all(type(value[name]) is int and value[name] >= 0
                for name in ("device", "directory_device", "directory_inode",
                             "inode")) and
            value["lock"] == "exclusive_nonblocking_pair_wide" and
            value["sha256"] == sha256_bytes(canonical_bytes(payload)),
            "CPU pair lease path/filesystem/digest differs")
    return value


def validate_isolation(value: Any, campaign: Mapping[str, Any],
                       durations_ns: int) -> dict[str, Any]:
    value = exact_keys(value, {
        "accepted", "after", "before", "benchmark_cpu", "delta", "pair_lease",
        "policy", "reserved_sibling", "schema"}, "isolation")
    require(value["schema"] == ISOLATION_SCHEMA and
            type(value["benchmark_cpu"]) is int and
            type(value["reserved_sibling"]) is int and
            value["benchmark_cpu"] == CPU and value["reserved_sibling"] == SIBLING,
            "isolation CPU identity differs")
    expected_policy = {
        "counter_source": "/proc/stat",
        "idle_fields": ["idle", "iowait"],
        "nonidle_fields": [
            "user", "nice", "system", "irq", "softirq", "steal"],
        "reserved_sibling_max_nonidle_jiffies": 0,
    }
    require(exact_json_equal(value["policy"], expected_policy),
            "isolation policy differs")
    before = exact_keys(value["before"], {
        "benchmark_cpu", "monotonic_ns", "reserved_sibling"},
        "isolation before")
    after = exact_keys(value["after"], {
        "benchmark_cpu", "monotonic_ns", "reserved_sibling"},
        "isolation after")
    require(type(before["monotonic_ns"]) is int and
            type(after["monotonic_ns"]) is int and
            after["monotonic_ns"] > before["monotonic_ns"] and
            after["monotonic_ns"] - before["monotonic_ns"] >= durations_ns,
            "isolation interval does not cover all child durations")
    benchmark = cpu_delta(
        before["benchmark_cpu"], after["benchmark_cpu"], CPU,
        "benchmark CPU")
    sibling = cpu_delta(
        before["reserved_sibling"], after["reserved_sibling"], SIBLING,
        "reserved sibling")
    require(exact_json_equal(value["delta"], {
                "benchmark_cpu": benchmark, "reserved_sibling": sibling}) and
            benchmark["nonidle_jiffies"] > 0 and
            sibling["total_jiffies"] > 0 and
            sibling["nonidle_jiffies"] == 0 and value["accepted"] is True,
            "recomputed sibling-isolation gate failed")
    validate_pair_lease(value["pair_lease"])
    return value


def validate_supervision(value: Any, campaign: Mapping[str, Any],
                         reservation: Mapping[str, Any],
                         isolation: Mapping[str, Any]) -> dict[str, Any]:
    value = exact_keys(value, {
        "campaign_sha256", "execution_nonce",
        "isolation_after_monotonic_ns", "isolation_before_monotonic_ns",
        "launch_cpus", "reservation_nonce", "reservation_sha256",
        "reserved_cpus", "runner_finished_monotonic_ns", "runner_pid",
        "runner_started_monotonic_ns", "schema"}, "supervision")
    require(value["schema"] == SUPERVISION_SCHEMA and
            type(value["execution_nonce"]) is str and
            HEX256.fullmatch(value["execution_nonce"]) is not None and
            type(value["runner_pid"]) is int and value["runner_pid"] > 0 and
            type(value["runner_started_monotonic_ns"]) is int and
            type(value["runner_finished_monotonic_ns"]) is int,
            "supervision identity differs")
    before = isolation["before"]["monotonic_ns"]
    after = isolation["after"]["monotonic_ns"]
    payload = reservation["payload"]
    require(exact_json_equal(value["launch_cpus"],
                             campaign["allowed_cpu_set_at_launch"]) and
            exact_json_equal(value["reserved_cpus"], [CPU, SIBLING]) and
            type(value["launch_cpus"]) is list and
            all(type(cpu) is int for cpu in value["launch_cpus"]) and
            type(value["reserved_cpus"]) is list and
            all(type(cpu) is int for cpu in value["reserved_cpus"]) and
            value["campaign_sha256"] == sha256_bytes(canonical_bytes(campaign)) and
            value["reservation_sha256"] == reservation["sha256"] and
            value["reservation_nonce"] == payload["nonce"] and
            value["isolation_before_monotonic_ns"] == before and
            value["isolation_after_monotonic_ns"] == after and
            0 <= value["runner_started_monotonic_ns"] <= before <= after <=
                value["runner_finished_monotonic_ns"],
            "supervision handshake/interval does not cross-bind the campaign")
    return value


def confidence_interval(log_ratios: Sequence[float]) -> dict[str, Any]:
    require(len(log_ratios) == ROUNDS and
            all(math.isfinite(item) for item in log_ratios),
            "ABBA inference lacks three finite round contrasts")
    mean = statistics.fmean(log_ratios)
    deviation = statistics.stdev(log_ratios)
    half_width = T_CRITICAL_DF2 * deviation / math.sqrt(ROUNDS)
    lower = math.exp(mean - half_width)
    upper = math.exp(mean + half_width)
    return {
        "ratio_orientation": "baseline_time_over_candidate_time",
        "independent_round_count": ROUNDS,
        "degrees_of_freedom": ROUNDS - 1,
        "independent_round_log_contrasts": list(log_ratios),
        "geometric_speedup": math.exp(mean),
        "ci95_lower": lower,
        "ci95_upper": upper,
        "faster_ci_excludes_one": lower > 1.0,
        "promotion_lower_bound_at_least_1_05": lower >= 1.05,
        "performance_result_does_not_affect_evidence_validity": True,
    }


def analyze(invocations: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    require(len(invocations) == ROUNDS * len(ORDER),
            "analysis requires exactly twelve invocations")
    observations = {"encode": [], "one_shot_encode": []}
    within = {"encode": [], "one_shot_encode": []}
    for round_index in range(ROUNDS):
        group = invocations[round_index * 4:(round_index + 1) * 4]
        require(tuple(item["implementation"] for item in group) == ORDER,
                f"round {round_index} is not ABBA")
        for metric in ("encode", "one_shot_encode"):
            contrasts: list[float] = []
            for baseline_index, candidate_index in ((0, 1), (3, 2)):
                baseline = statistics.median(
                    group[baseline_index]["normalized"]["encode"])
                candidate_samples = group[candidate_index]["normalized"][
                    "encode" if metric == "encode" else "one_shot_encode"]
                candidate = statistics.median(candidate_samples)
                contrasts.append(math.log(baseline / candidate))
            within[metric].append(contrasts)
            observations[metric].append(statistics.fmean(contrasts))
    result: dict[str, Any] = {}
    for metric in ("encode", "one_shot_encode"):
        result[metric] = confidence_interval(observations[metric])
        result[metric]["within_round_log_contrasts"] = within[metric]
        result[metric]["constituent_pair_count"] = 2 * ROUNDS
    return {CELL["identifier"]: result}


def validate_campaign(value: Any) -> dict[str, Any]:
    value = exact_keys(value, CAMPAIGN_KEYS, "campaign")
    timeout = value["timeout_seconds"]
    allowed = value["allowed_cpu_set_at_launch"]
    require(type(timeout) in (int, float) and not isinstance(timeout, bool) and
            math.isfinite(float(timeout)) and 0 < float(timeout) <= 3600.0 and
            type(allowed) is list and allowed == sorted(set(allowed)) and
            3 <= len(allowed) <= MAX_CPU_LIST_ENTRIES and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in allowed) and
            CPU in allowed and SIBLING in allowed and len(set(allowed) - {CPU, SIBLING}) > 0,
            "campaign launch topology or timeout differs")
    expected = {
        "rounds": ROUNDS, "order": list(ORDER), "cells": [dict(CELL)],
        "candidate_mode": "auto", "batch": 1, "reuse": REUSE,
        "iterations": ITERATIONS, "warmup": WARMUP, "threads": 1,
        "child_environment": CHILD_ENVIRONMENT, "benchmark_cpu": CPU,
        "reserved_sibling": SIBLING, "statistics": STATISTICS_POLICY,
    }
    for name, expected_value in expected.items():
        require(exact_json_equal(value.get(name), expected_value),
                f"campaign {name} differs from frozen v17")
    return value


def replay(manifest_path: Path) -> dict[str, Any]:
    manifest_path = canonical_existing_file(manifest_path)
    root = manifest_path.parent.resolve(strict=True)
    require(manifest_path.parent == root, "manifest parent is aliased")
    manifest_bytes, manifest_file = stable_file(
        manifest_path, MAX_MANIFEST_BYTES)
    manifest = strict_json_bytes(manifest_bytes, "manifest")
    manifest = exact_keys(manifest, MANIFEST_KEYS, "manifest")
    verify_signature(manifest, "manifest")
    validate_timestamp(manifest["created_utc"], "manifest created_utc")
    require(manifest["schema"] == MANIFEST_SCHEMA and
            manifest["valid"] is True and
            manifest["validity_is_independent_of_speed"] is True,
            "manifest is not valid v17 evidence")
    raw_info = exact_keys(manifest["raw"], MANIFEST_RAW_KEYS,
                          "manifest raw identity")
    require(raw_info["path"] == "raw.json", "manifest raw path changed")
    raw_bytes, raw_file = safe_retained_file(
        root, raw_info["path"], MAX_RAW_BYTES, "raw bundle")
    require(type(raw_info["size"]) is int and raw_info["size"] == len(raw_bytes) and
            raw_info["sha256"] == raw_file["sha256"] and
            type(raw_info["payload_digest"]) is str and
            HEX256.fullmatch(raw_info["payload_digest"]) is not None,
            "manifest/raw file identity differs")
    raw = strict_json_bytes(raw_bytes, "raw bundle")
    raw = exact_keys(raw, RAW_KEYS, "raw bundle")
    verify_signature(raw, "raw bundle")
    validate_timestamp(raw["created_utc"], "raw created_utc")
    require(raw["schema"] == RAW_SCHEMA and
            raw["validity_is_independent_of_speed"] is True and
            raw_info["payload_digest"] == raw["digest"],
            "raw bundle is not signed v17 evidence")
    campaign = validate_campaign(raw["campaign"])
    require(type(raw["identities_initial"]) is dict and
            exact_json_equal(raw["identities_initial"], raw["identities_final"]),
            "input identities changed during the campaign")
    identities = validate_native_build(
        raw["input_specification"], raw["identities_initial"])
    reservation = validate_reservation(raw["reservation"])
    validate_sealed_executables(raw["executable_snapshots"], identities)
    input_spec = raw["input_specification"]
    require(input_spec["taskset"] == identities.get("taskset", {}).get("path"),
            "taskset path differs from retained identity")
    invocations = raw["invocations"]
    require(type(invocations) is list and len(invocations) == 12,
            "v17 campaign does not contain twelve invocations")
    expected_sequence = [
        (round_index, slot, role)
        for round_index in range(ROUNDS)
        for slot, role in enumerate(ORDER)
    ]
    accepted: list[dict[str, Any]] = []
    common_digests: dict[str, str] | None = None
    total_duration = 0
    raw_files: list[dict[str, Any]] = []
    timer_minima = {
        "baseline_encode": math.inf,
        "candidate_ordinary_encode": math.inf,
        "candidate_one_shot_encode": math.inf,
    }
    for invocation, (round_index, slot, role) in zip(
            invocations, expected_sequence):
        invocation = exact_keys(invocation, INVOCATION_KEYS,
                                f"invocation {round_index}/{slot}")
        require(invocation["cell_id"] == CELL["identifier"] and
                type(invocation["round"]) is int and
                invocation["round"] == round_index and
                type(invocation["slot"]) is int and invocation["slot"] == slot and
                invocation["implementation"] == role and
                invocation["command"] == expected_command(role, input_spec["taskset"]) and
                invocation["execution_protocol"] == SEALED_EXECUTABLE_PROTOCOL and
                exact_json_equal(invocation["executable_snapshot"],
                                 raw["executable_snapshots"][role]) and
                exact_json_equal(invocation["environment"], CHILD_ENVIRONMENT) and
                type(invocation["pinned_cpu"]) is int and
                invocation["pinned_cpu"] == CPU and
                type(invocation["duration_ns"]) is int and
                invocation["duration_ns"] > 0 and
                type(invocation["returncode"]) is int and
                invocation["returncode"] == 0,
                f"invocation {round_index}/{slot} contract differs")
        validate_timestamp(
            invocation["started_utc"], f"invocation {round_index}/{slot} start")
        require(exact_json_equal(invocation["identity_before"], identities) and
                exact_json_equal(invocation["identity_after"], identities) and
                exact_json_equal(invocation["reservation_before"], reservation) and
                exact_json_equal(invocation["reservation_after"], reservation),
                f"invocation {round_index}/{slot} provenance changed")
        total_duration += invocation["duration_ns"]
        for stream_name, maximum in (
                ("stdout", MAX_STDOUT_BYTES), ("stderr", MAX_STDERR_BYTES)):
            stream = exact_keys(invocation[stream_name], STREAM_KEYS,
                                f"invocation {round_index}/{slot} {stream_name}")
            expected_stream_path = (
                f"invocations/{CELL['identifier']}/round-{round_index}/"
                f"slot-{slot}-{role}.{stream_name}"
            )
            require(stream["path"] == expected_stream_path,
                    f"invocation {round_index}/{slot} {stream_name} path differs")
            data, observed = safe_retained_file(
                root, stream["path"], maximum,
                f"invocation {round_index}/{slot} {stream_name}")
            require(type(stream["size"]) is int and stream["size"] == len(data) and
                    type(stream["sha256"]) is str and
                    HEX256.fullmatch(stream["sha256"]) is not None and
                    stream["sha256"] == observed["sha256"],
                    f"invocation {round_index}/{slot} {stream_name} bytes differ")
            raw_files.append({
                "path": stream["path"], "size": len(data),
                "sha256": observed["sha256"],
            })
            if stream_name == "stdout":
                parsed = strict_json_bytes(
                    data, f"invocation {round_index}/{slot} stdout")
                require(exact_json_equal(parsed, invocation["result"]),
                        f"invocation {round_index}/{slot} parsed stdout differs")
        normalized = validate_child_result(
            role, invocation["result"], campaign, identities)
        require(exact_json_equal(normalized, invocation["normalized"]),
                f"invocation {round_index}/{slot} normalized result was edited")
        validate_wall_duration(
            invocation["duration_ns"], normalized, role,
            f"invocation {round_index}/{slot}")
        digests = normalized["digests"]
        if common_digests is None:
            common_digests = dict(digests)
        else:
            require(exact_json_equal(common_digests, digests),
                    "baseline/candidate workload digests differ")
        timer_minima[
            "baseline_encode" if role == "baseline"
            else "candidate_ordinary_encode"] = min(
                timer_minima[
                    "baseline_encode" if role == "baseline"
                    else "candidate_ordinary_encode"],
                normalized["retained_timer_window_us"]["encode"])
        if role == "candidate":
            timer_minima["candidate_one_shot_encode"] = min(
                timer_minima["candidate_one_shot_encode"],
                normalized["retained_timer_window_us"]["one_shot_encode"])
        accepted.append(dict(invocation, normalized=normalized))
    isolation = validate_isolation(raw["isolation"], campaign, total_duration)
    supervision = validate_supervision(
        raw["supervision"], campaign, reservation, isolation)
    host_initial = validate_host(
        raw["host_initial"], "initial", campaign["allowed_cpu_set_at_launch"])
    host_final = validate_host(
        raw["host_final"], "final", campaign["allowed_cpu_set_at_launch"])
    require(exact_json_equal(host_initial, host_final),
            "host identity changed during the campaign")
    analysis = analyze(accepted)
    require(exact_json_equal(analysis, raw["analysis"]),
            "raw clustered ABBA analysis differs from retained samples")
    for name, expected_value in (
            ("campaign", campaign), ("host", host_initial),
            ("isolation", isolation), ("reservation", reservation),
            ("supervision", supervision), ("identities", identities),
            ("executable_snapshots", raw["executable_snapshots"]),
            ("analysis", analysis)):
        require(exact_json_equal(manifest[name], expected_value),
                f"manifest {name} differs from raw evidence")
    require(all(math.isfinite(item) and item >= TIMER_FLOOR_US
                for item in timer_minima.values()),
            "global retained timer floor failed")
    raw_manifest_hash = sha256_bytes(canonical_bytes(sorted(
        raw_files, key=lambda item: item["path"])))
    auditor_bytes, auditor_file = stable_file(
        Path(__file__).resolve(strict=True), 8 * 1024 * 1024)
    del auditor_bytes
    return {
        "schema": AUDIT_SCHEMA,
        "status": "complete",
        "audit_passed": True,
        "timing_performed": False,
        "benchmark_executed": False,
        "auditor": auditor_file,
        "manifest": manifest_file,
        "raw": {
            **raw_file,
            "payload_digest": raw["digest"],
        },
        "contract": {
            "manifest_schema": MANIFEST_SCHEMA,
            "raw_schema": RAW_SCHEMA,
            "baseline_schema": BASELINE_SCHEMA,
            "candidate_schema": CANDIDATE_SCHEMA,
            "baseline_main_commit": MAIN_COMMIT,
            "cell": dict(CELL),
            "rounds": ROUNDS,
            "order": list(ORDER),
            "reuse": REUSE,
            "iterations": ITERATIONS,
            "warmup": WARMUP,
            "cpu": CPU,
            "sibling": SIBLING,
            "timer_floor_us": TIMER_FLOOR_US,
            "route_contract": ROUTE_CONTRACT,
            "build_layout_contract":
                "single-config Unix Makefiles Release only",
            "ordinary_encode_api": ORDINARY_ENCODE_API,
            "one_shot_encode_api": ONE_SHOT_ENCODE_API,
        },
        "replay": {
            "invocation_count": len(accepted),
            "candidate_route_attestation_count": sum(
                item["implementation"] == "candidate" for item in accepted),
            "baseline_required_source_object_count":
                len(identities["baseline_build"]
                    ["validated_compile_commands"]
                    ["required_source_object_pairs"]),
            "candidate_required_source_object_count":
                len(identities["candidate_build"]
                    ["validated_compile_commands"]
                    ["required_source_object_pairs"]),
            "baseline_configured_compile_action_count":
                identities["baseline_build"]["validated_compile_commands"]
                ["entry_count"],
            "candidate_configured_compile_action_count":
                identities["candidate_build"]["validated_compile_commands"]
                ["entry_count"],
            "workload_digests": common_digests,
            "minimum_retained_timer_window_us": timer_minima,
            "analysis": analysis,
            "retained_stream_count": len(raw_files),
            "retained_stream_manifest_sha256": raw_manifest_hash,
        },
        "gate_results": {
            "manifest_raw_digest_chain": True,
            "v17_schema_and_contract": True,
            "native_exact_main_identity": True,
            "full_git_source_object_closure": True,
            "full_build_source_object_archive_link_closure": True,
            "multi_config_builds_excluded_fail_closed": True,
            "generated_attestation_and_configuration_bytes": True,
            "canonical_runtime_dependency_closure": True,
            "reservation_pair_lease_supervision_closure": True,
            "complete_host_topology_identity": True,
            "candidate_gfni_route_all_six_children": True,
            "sealed_executable_identity": True,
            "zero_reserved_sibling_nonidle_jiffies": True,
            "workload_digest_identity": True,
            "retained_timer_window_floor": True,
            "wall_clock_encloses_validated_encode_and_setup_timer_work": True,
            "clustered_abba_rederived": True,
        },
        "reporting_policy": {
            "ratio_orientation": "baseline_time_over_candidate_time",
            "ordinary_and_one_shot_ratios_emitted_separately": True,
            "combined_or_stacked_ratio_emitted": False,
            "same_binary_gfni_over_avx2_is_a_separate_campaign": True,
            "performance_result_does_not_affect_evidence_validity": True,
            "production_selector_source_bytes_offline_rederived": False,
            "production_selector_bound_by_candidate_commit_tree_object_and_"
            "runtime_attestation": True,
            "multi_config_ninja_graphs_are_outside_this_auditor_contract": True,
        },
    }


def write_canonical(path: Path, value: Mapping[str, Any]) -> None:
    require(path.is_absolute(), "audit output path is not absolute")
    parent = path.parent.resolve(strict=True)
    require(parent == path.parent and not path.exists(),
            "audit output parent is aliased or output already exists")
    data = canonical_bytes(value) + b"\n"
    descriptor = os.open(
        path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0),
        0o600)
    try:
        written = 0
        while written < len(data):
            count = os.write(descriptor, data[written:])
            require(count > 0, "short audit output write")
            written += count
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    directory = os.open(parent, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def expect_failure(function: Any, label: str) -> None:
    try:
        function()
    except (AuditError, UnicodeError, json.JSONDecodeError, OSError):
        return
    raise AuditError(f"self-test accepted {label}")


def self_test() -> None:
    def clone(value: Any) -> Any:
        return strict_json_bytes(canonical_bytes(value), "self-test clone")

    def artifact(path: str, kind: str, data: bytes = b"fixture",
                 *, mode: int = 0o644, mtime_ns: int = 100) -> dict[str, Any]:
        return {
            "path": path, "kind": kind, "size": len(data), "mode": mode,
            "mtime_ns": mtime_ns, "sha256": sha256_bytes(data),
        }

    require(
        STATISTICS_POLICY.get(
            "ratios_are_separate_correlated_and_must_not_be_multiplied") is True
        and "ratios_are_independent_and_must_not_be_multiplied"
        not in STATISTICS_POLICY,
        "correlated non-composable ratio policy changed",
    )
    for kind in (
            "archive", "object_file", "source_file", "executable",
            "shared_library", "dynamic_loader"):
        empty = artifact(f"/fixture/empty-{kind}", kind, b"",
                         mode=0o755 if kind == "executable" else 0o644)
        expect_failure(lambda empty=empty, kind=kind: validate_artifact(
            empty, f"empty {kind}", kind), f"empty required {kind}")
    for kind in ("archiver", "ranlib", "compiler", "build_tool"):
        nonexecutable = artifact(
            f"/fixture/nonexec-{kind}", kind, mode=0o644)
        expect_failure(
            lambda nonexecutable=nonexecutable, kind=kind: validate_artifact(
                nonexecutable, f"nonexec {kind}", kind),
            f"non-executable {kind}")
    aliased_artifact = artifact(
        "/fixture/build/../evil/object.o", "object_file")
    expect_failure(lambda: validate_artifact(
        aliased_artifact, "aliased artifact", "object_file"),
                   "artifact path containing parent traversal")
    make_layout = {
        "CMAKE_GENERATOR": "Unix Makefiles", "CMAKE_BUILD_TYPE": "Release",
        # A single-config cache may retain the configured multi-config list;
        # it is not selected by this generator.
        "CMAKE_CONFIGURATION_TYPES": "Debug;Release",
    }
    require(validate_build_layout(make_layout, "self-test") is None,
            "Unix Makefiles Release layout was not selected")
    ninja_layout = {
        "CMAKE_GENERATOR": "Ninja Multi-Config", "CMAKE_BUILD_TYPE": "",
        "CMAKE_CONFIGURATION_TYPES": "Debug;Release",
        "CMAKE_MAKE_PROGRAM": CANONICAL_NINJA_PATH,
    }
    expect_failure(lambda: validate_build_layout(ninja_layout, "self-test"),
                   "multi-config graph without independently expanded recipes")
    acquisition = {
        "candidate_source_root": "/fixture/candidate",
        "runner": ("/fixture/candidate/experiments/leopard2/main_compare/"
                   "run_abba.py"),
        "taskset": CANONICAL_TASKSET_PATH, "ldd": CANONICAL_LDD_PATH,
    }
    validate_acquisition_paths(acquisition)
    bad_acquisition = clone(acquisition)
    bad_acquisition["taskset"] = "/tmp/taskset"
    expect_failure(lambda: validate_acquisition_paths(bad_acquisition),
                   "substituted taskset executable")
    bad_acquisition = clone(acquisition)
    bad_acquisition["runner"] = "/tmp/run_abba.py"
    expect_failure(lambda: validate_acquisition_paths(bad_acquisition),
                   "runner outside retained candidate source")
    input_fixture = {
        "runner": acquisition["runner"], "taskset": CANONICAL_TASKSET_PATH,
        "ldd": CANONICAL_LDD_PATH,
        "baseline_executable": "/fixture/baseline-build/benchmark",
        "candidate_executable": "/fixture/candidate-build/benchmark",
        "baseline_archive": "/fixture/baseline-build/library.a",
        "candidate_archive": "/fixture/candidate-build/library.a",
        "baseline_build_dir": "/fixture/baseline-build",
        "candidate_build_dir": "/fixture/candidate-build",
        "baseline_source_root": "/fixture/baseline",
        "candidate_source_root": "/fixture/candidate",
        "candidate_commit": "a" * 40, "baseline_pure_avx2": False,
    }
    validate_input_specification(input_fixture)
    bad_input = clone(input_fixture)
    bad_input["baseline_build_dir"] = "/fixture/build/../evil"
    expect_failure(lambda: validate_input_specification(bad_input),
                   "input build path containing parent traversal")
    bad_input = clone(input_fixture)
    bad_input["candidate_archive"] = "/fixture//candidate-build/library.a"
    expect_failure(lambda: validate_input_specification(bad_input),
                   "input path containing duplicate separator")
    archive_tool = artifact(
        CANONICAL_ARCHIVER_RESOLVED_PATH, "archiver", mode=0o755)
    ranlib_tool = artifact(
        CANONICAL_RANLIB_RESOLVED_PATH, "ranlib", mode=0o755)
    archive_invocations = {
        "archiver": {
            "invocation": CANONICAL_ARCHIVER_PATH,
            "resolved_path": CANONICAL_ARCHIVER_RESOLVED_PATH,
        },
        "ranlib": {
            "invocation": CANONICAL_RANLIB_PATH,
            "resolved_path": CANONICAL_RANLIB_RESOLVED_PATH,
        },
    }
    validate_archive_tools(
        archive_tool, ranlib_tool, archive_invocations, "self-test")
    bad_archive_tool = clone(archive_tool)
    bad_archive_tool["path"] = "/tmp/evil-ar"
    bad_archive_invocations = clone(archive_invocations)
    bad_archive_invocations["archiver"] = {
        "invocation": CANONICAL_ARCHIVER_PATH,
        "resolved_path": "/tmp/evil-ar"}
    expect_failure(lambda: validate_archive_tools(
        bad_archive_tool, ranlib_tool, bad_archive_invocations, "self-test"),
                   "candidate-only substituted archiver")
    fresh_pair = {
        "source": artifact(
            "/fixture/source.cpp", "source_file", mtime_ns=10),
        "object": artifact(
            "/fixture/build/source.cpp.o", "object_file", mtime_ns=20),
    }
    validate_source_object_pair(fresh_pair, "self-test pair", "/fixture/build")
    stale_pair = clone(fresh_pair)
    stale_pair["source"]["mtime_ns"] = 21
    expect_failure(lambda: validate_source_object_pair(
        stale_pair, "self-test stale pair", "/fixture/build"),
                   "object predating required source")
    for payload in (
            b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":1e309}', b'\xff'):
        expect_failure(lambda payload=payload: strict_json_bytes(
            payload, "adversarial JSON"), "adversarial JSON")
    expect_failure(lambda: strict_json_bytes(
        b"[" * 5000 + b"0" + b"]" * 5000, "recursive JSON"),
        "recursive JSON")
    expect_failure(lambda: strict_json_bytes(
        b'{"x":' + b"9" * 100_000 + b"}", "huge integer JSON"),
        "huge integer JSON")
    require(not exact_json_equal(True, 1) and
            not exact_json_equal(1, 1.0) and
            not exact_json_equal({"x": [1]}, {"x": [1.0]}),
            "type-aware JSON equality accepted numeric coercion")
    losses = expected_losses()
    require(len(losses) == CELL["losses"] and losses == sorted(set(losses)) and
            losses[0] >= 0 and losses[-1] < CELL["k"],
            "independent loss reconstruction failed")
    candidate_command = expected_command("candidate", "/usr/bin/taskset")
    require("--measure-one-shot-decode" not in candidate_command and
            candidate_command[-8:] == [
                "--skip-legacy", "--retain-samples",
                "--measure-one-shot-encode", "--attest-source",
                "--auto-gf16-gfni-encode-mode", "1", "--json", "-"],
            "candidate command reconstruction failed")
    boundary = {
        "median_us_per_batch_call": 2500.0,
        "mad_us_per_batch_call": 0.0,
        "minimum_us_per_batch_call": 2500.0,
        "maximum_us_per_batch_call": 2500.0,
        "samples_us_per_batch_call": [2500.0] * ITERATIONS,
    }
    samples = validate_summary(boundary, "timer boundary")
    require(min(samples) * REUSE == TIMER_FLOOR_US,
            "timer floor boundary failed")
    negative_mad = dict(boundary)
    negative_mad["mad_us_per_batch_call"] = -0.000001
    expect_failure(lambda: validate_summary(negative_mad, "negative MAD"),
                   "negative MAD hidden within numeric tolerance")
    timer_work = {
        "encode": [2500.0] * ITERATIONS,
        "one_shot_encode": [2000.0] * ITERATIONS,
        "codec_setup": [100.0] * ITERATIONS,
    }
    ordinary_work_us = sum(timer_work["encode"]) * REUSE
    one_shot_work_us = sum(timer_work["one_shot_encode"]) * REUSE
    setup_work_us = sum(timer_work["codec_setup"])
    validate_wall_duration(
        int((ordinary_work_us + one_shot_work_us + setup_work_us) * 1000),
        timer_work, "candidate", "self-test wall timer")
    expect_failure(lambda: validate_wall_duration(
        int(max(ordinary_work_us, one_shot_work_us, setup_work_us) * 1000),
        timer_work, "candidate", "self-test wall timer"),
                   "wall duration fitting each sequential loop but not their sum")
    expect_failure(lambda: validate_wall_duration(
        1, timer_work, "candidate", "self-test wall timer"),
                   "one-nanosecond wall duration for retained timer work")
    too_short = dict(boundary)
    too_short["samples_us_per_batch_call"] = [2499.0] + [2500.0] * 8
    too_short["minimum_us_per_batch_call"] = 2499.0
    validate_summary(too_short, "short timer")
    expect_failure(lambda: require(
        min(too_short["samples_us_per_batch_call"]) * REUSE >= TIMER_FLOOR_US,
        "timer floor"), "short retained timer")
    fake: list[dict[str, Any]] = []
    for round_index in range(ROUNDS):
        for slot, role in enumerate(ORDER):
            normalized = {
                "encode": [10.0 if role == "baseline" else 5.0] * ITERATIONS,
            }
            if role == "candidate":
                normalized["one_shot_encode"] = [4.0] * ITERATIONS
            fake.append({
                "round": round_index, "slot": slot,
                "implementation": role, "normalized": normalized,
            })
    inferred = analyze(fake)[CELL["identifier"]]
    require(close_enough(inferred["encode"]["geometric_speedup"], 2.0) and
            close_enough(
                inferred["one_shot_encode"]["geometric_speedup"], 2.5),
            "independent ABBA inference failed")
    relabeled = {"schema": "leopard2-main-compare-raw/v16"}
    expect_failure(lambda: require(relabeled["schema"] == RAW_SCHEMA,
                                   "schema relabel"),
                   "v16 schema relabeled as v17")
    route = {
        "auto_gf16_gfni_encode_kernel_available": True,
        "auto_gf16_gfni_encode_observed_call_count": 2,
    }
    route["auto_gf16_gfni_encode_kernel_available"] = 1
    expect_failure(lambda: require(
        type(route["auto_gf16_gfni_encode_kernel_available"]) is bool,
        "route bool type"), "route bool/int coercion")
    digest_fixture = {
        "algorithm": "fnv1a64", "original_data": "0" * 16,
        "transmitted_parity": "1" * 16, "recovered_originals": "2" * 16,
    }
    normalized_digests = validate_digests(digest_fixture, "self-test digests")
    require(set(normalized_digests) == {
                "original_data", "transmitted_parity", "recovered_originals"} and
            "algorithm" not in normalized_digests,
            "digest normalization retained the transport algorithm field")
    before = {"cpu": SIBLING, "fields": {
        name: 10 for name in CPU_STAT_FIELDS}, "total_jiffies": 80}
    after = {"cpu": SIBLING, "fields": {
        **before["fields"], "idle": 11}, "total_jiffies": 81}
    clean = cpu_delta(before, after, SIBLING, "self-test sibling")
    require(clean["nonidle_jiffies"] == 0 and clean["total_jiffies"] == 1,
            "clean sibling delta failed")
    contaminated = {"cpu": SIBLING, "fields": {
        **after["fields"], "system": 11}, "total_jiffies": 82}
    dirty = cpu_delta(before, contaminated, SIBLING, "self-test sibling")
    expect_failure(lambda: require(dirty["nonidle_jiffies"] == 0,
                                   "sibling contamination"),
                   "reserved sibling contamination")

    frozen_campaign = {
        "rounds": ROUNDS, "order": list(ORDER), "cells": [dict(CELL)],
        "candidate_mode": "auto", "batch": 1, "reuse": REUSE,
        "iterations": ITERATIONS, "warmup": WARMUP, "threads": 1,
        "child_environment": CHILD_ENVIRONMENT, "benchmark_cpu": CPU,
        "reserved_sibling": SIBLING, "timeout_seconds": 3600,
        "statistics": STATISTICS_POLICY,
        "allowed_cpu_set_at_launch": [0, CPU, SIBLING],
    }
    validate_campaign(frozen_campaign)
    bad_campaign = clone(frozen_campaign)
    bad_campaign["timeout_seconds"] = 3600.0001
    expect_failure(lambda: validate_campaign(bad_campaign),
                   "campaign timeout above producer maximum")
    bad_campaign = clone(frozen_campaign)
    bad_campaign["allowed_cpu_set_at_launch"][1] = float(CPU)
    expect_failure(lambda: validate_campaign(bad_campaign),
                   "floating campaign launch CPU")

    reservation_payload = {
        "benchmark_cpu": CPU, "nonce": "0123456789abcdef",
        "owner": "self-test", "reserved_sibling": SIBLING,
        "schema": RESERVATION_SCHEMA, "status": "held",
    }
    reservation = {
        "path": "/run/leopard2-self-test-reservation.json",
        "sha256": sha256_bytes(canonical_bytes(reservation_payload)),
        "payload": reservation_payload, "lock": "exclusive_nonblocking",
    }
    validate_reservation(reservation)
    bad_reservation = clone(reservation)
    bad_reservation["payload"]["benchmark_cpu"] = float(CPU)
    expect_failure(lambda: validate_reservation(bad_reservation),
                   "floating reservation CPU")
    bad_reservation = clone(reservation)
    bad_reservation["payload"]["owner"] = ""
    bad_reservation["sha256"] = sha256_bytes(canonical_bytes(
        bad_reservation["payload"]))
    expect_failure(lambda: validate_reservation(bad_reservation),
                   "coherently resigned ownerless reservation")

    lease_payload = {
        "cpus": [CPU, SIBLING], "schema": PAIR_LEASE_SCHEMA, "uid": 1000,
    }
    pair_lease = {
        "device": 1, "directory_device": 1, "directory_inode": 2,
        "inode": 3, "lock": "exclusive_nonblocking_pair_wide",
        "path": ("/run/user/1000/leopard2-cpu-leases/"
                 "leopard2-cpu-pair-1000-52-116.lock"),
        "payload": lease_payload,
        "sha256": sha256_bytes(canonical_bytes(lease_payload)),
    }
    validate_pair_lease(pair_lease)
    bad_lease = clone(pair_lease)
    bad_lease["payload"]["cpus"] = [float(CPU), float(SIBLING)]
    bad_lease["sha256"] = sha256_bytes(canonical_bytes(bad_lease["payload"]))
    expect_failure(lambda: validate_pair_lease(bad_lease),
                   "coherently resigned floating pair lease CPUs")

    zero_fields = {name: 0 for name in CPU_STAT_FIELDS}
    benchmark_after_fields = {**zero_fields, "user": 1}
    sibling_after_fields = {**zero_fields, "idle": 1}
    isolation = {
        "accepted": True,
        "before": {
            "benchmark_cpu": {
                "cpu": CPU, "fields": zero_fields, "total_jiffies": 0},
            "monotonic_ns": 20,
            "reserved_sibling": {
                "cpu": SIBLING, "fields": zero_fields, "total_jiffies": 0},
        },
        "after": {
            "benchmark_cpu": {
                "cpu": CPU, "fields": benchmark_after_fields,
                "total_jiffies": 1},
            "monotonic_ns": 30,
            "reserved_sibling": {
                "cpu": SIBLING, "fields": sibling_after_fields,
                "total_jiffies": 1},
        },
        "benchmark_cpu": CPU,
        "delta": {
            "benchmark_cpu": {
                "cpu": CPU, "fields": benchmark_after_fields,
                "idle_jiffies": 0, "nonidle_jiffies": 1,
                "total_jiffies": 1},
            "reserved_sibling": {
                "cpu": SIBLING, "fields": sibling_after_fields,
                "idle_jiffies": 1, "nonidle_jiffies": 0,
                "total_jiffies": 1},
        },
        "pair_lease": pair_lease,
        "policy": {
            "counter_source": "/proc/stat",
            "idle_fields": ["idle", "iowait"],
            "nonidle_fields": [
                "user", "nice", "system", "irq", "softirq", "steal"],
            "reserved_sibling_max_nonidle_jiffies": 0,
        },
        "reserved_sibling": SIBLING, "schema": ISOLATION_SCHEMA,
    }
    validate_isolation(isolation, frozen_campaign, 5)
    bad_isolation = clone(isolation)
    bad_isolation["benchmark_cpu"] = float(CPU)
    expect_failure(lambda: validate_isolation(
        bad_isolation, frozen_campaign, 5), "floating isolation CPU")

    supervision_campaign = frozen_campaign
    supervision_isolation = isolation
    supervision = {
        "campaign_sha256": sha256_bytes(canonical_bytes(supervision_campaign)),
        "execution_nonce": "a" * 64,
        "isolation_after_monotonic_ns": 30,
        "isolation_before_monotonic_ns": 20,
        "launch_cpus": [0, CPU, SIBLING],
        "reservation_nonce": reservation_payload["nonce"],
        "reservation_sha256": reservation["sha256"],
        "reserved_cpus": [CPU, SIBLING],
        "runner_finished_monotonic_ns": 40, "runner_pid": 123,
        "runner_started_monotonic_ns": 10,
        "schema": SUPERVISION_SCHEMA,
    }
    validate_supervision(supervision, supervision_campaign, reservation,
                         supervision_isolation)
    bad_supervision = clone(supervision)
    del bad_supervision["runner_pid"]
    expect_failure(lambda: validate_supervision(
        bad_supervision, supervision_campaign, reservation,
        supervision_isolation), "supervision missing runner PID")
    bad_supervision = clone(supervision)
    bad_supervision["reserved_cpus"] = [float(CPU), float(SIBLING)]
    expect_failure(lambda: validate_supervision(
        bad_supervision, supervision_campaign, reservation,
        supervision_isolation), "floating supervision CPUs")
    bad_supervision = clone(supervision)
    bad_supervision["runner_finished_monotonic_ns"] = 29
    expect_failure(lambda: validate_supervision(
        bad_supervision, supervision_campaign, reservation,
        supervision_isolation), "supervision interval not enclosing isolation")

    external_inputs = [
        {
            "operand": "/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so",
            "role": "openmp_runtime_shared",
            "artifact": artifact(
                "/usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0",
                "shared_library"),
            "lexical_symlink_chain": [
                {
                    "path": "/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so",
                    "target": "../../../x86_64-linux-gnu/libgomp.so.1",
                    "mode": 0o777,
                },
                {
                    "path": "/usr/lib/x86_64-linux-gnu/libgomp.so.1",
                    "target": "libgomp.so.1.0.0", "mode": 0o777,
                },
            ],
        },
        {
            "operand": "/usr/lib/x86_64-linux-gnu/libpthread.a",
            "role": "pthread_support_archive",
            "artifact": artifact(
                "/usr/lib/x86_64-linux-gnu/libpthread.a", "archive"),
            "lexical_symlink_chain": [],
        },
    ]
    validate_external_link_inputs(external_inputs, "self-test link")
    bad_external = clone(external_inputs)
    bad_external[0]["lexical_symlink_chain"][1]["target"] = \
        "libgomp.so.1.0.1"
    expect_failure(lambda: validate_external_link_inputs(
        bad_external, "self-test link"), "discontinuous external link target")
    bad_external = clone(external_inputs)
    bad_external[0]["lexical_symlink_chain"].insert(1, {
        "path": "/usr/lib/x86_64-linux-gnu/libgomp.so.1",
        "target": "libgomp.so.1", "mode": 0o777,
    })
    expect_failure(lambda: validate_external_link_inputs(
        bad_external, "self-test link"), "repeated symlink path cycle")
    bad_external = clone(external_inputs)
    bad_external[0]["lexical_symlink_chain"][0]["target"] += "\0suffix"
    expect_failure(lambda: validate_external_link_inputs(
        bad_external, "self-test link"), "NUL-bearing symlink target")

    ldd_text = (
        "/usr/lib64/ld-linux-x86-64.so.2 (<ASLR_LOAD_ADDRESS>)\n"
        "libc.so.6 => /usr/lib/x86_64-linux-gnu/libc.so.6 "
        "(<ASLR_LOAD_ADDRESS>)\n"
        "linux-vdso.so.1 (<ASLR_LOAD_ADDRESS>)\n"
    )
    ldd_identity = text_identity(ldd_text)
    runtime_closure = {
        "executable": "/fixture/executable",
        "dependencies": [
            {
                "soname": "ld-linux-x86-64.so.2",
                "loader_path": "/usr/lib64/ld-linux-x86-64.so.2",
                "file": artifact(
                    "/usr/lib64/ld-linux-x86-64.so.2", "dynamic_loader"),
            },
            {
                "soname": "libc.so.6",
                "loader_path": "/usr/lib/x86_64-linux-gnu/libc.so.6",
                "file": artifact(
                    "/usr/lib/x86_64-linux-gnu/libc.so.6", "shared_library"),
            },
            {"soname": "linux-vdso.so.1", "virtual": True},
        ],
        "canonical_ldd_output": {
            "schema": CANONICAL_LDD_SCHEMA,
            "normalization": CANONICAL_LDD_NORMALIZATION,
            **ldd_identity,
        },
    }
    validate_runtime_closure(runtime_closure, "self-test runtime",
                             "/fixture/executable")
    bad_runtime = clone(runtime_closure)
    bad_runtime["dependencies"] = []
    expect_failure(lambda: validate_runtime_closure(
        bad_runtime, "self-test runtime", "/fixture/executable"),
                   "empty runtime dependency closure")

    compile_specification = {
        "candidate_source_root": "/fixture/candidate",
        "candidate_build_dir": "/fixture/candidate-build",
        "baseline_source_root": "/fixture/baseline",
    }
    compile_configuration = {"configuration_sha256": "b" * 64}
    gfni_source = "/fixture/candidate/Leopard2BackendGFNI.cpp"
    gfni_argv = expected_candidate_compile_argv(
        gfni_source, compile_specification, "/usr/bin/g++", None,
        compile_configuration, "/usr/bin/g++")
    require(gfni_argv.count("-mgfni") == 1 and
            gfni_argv.count("-mno-avx512f") == 1 and
            not any(token.startswith("-mavx512") for token in gfni_argv),
            "GFNI translation-unit ISA partition differs")
    bad_argv = [*gfni_argv, "-mavx512f", "-include", "/tmp/evil.h"]
    expect_failure(lambda: require(
        exact_json_equal(bad_argv, gfni_argv), "compiler argv differs"),
                   "appended compiler flags and injected input")
    portable_argv = expected_candidate_compile_argv(
        "/fixture/candidate/leopard2.cpp", compile_specification,
        "/usr/bin/g++", None, compile_configuration, "/usr/bin/g++")
    require("-mgfni" not in portable_argv and
            not any(token.startswith("-mavx") for token in portable_argv),
            "portable candidate translation unit gained ISA flags")

    source_commit, source_tree = "c" * 40, "d" * 40
    header = attestation_header_bytes(source_commit, source_tree)
    attestation = {
        "artifact": artifact(
            "/fixture/candidate-build/generated/"
            "leopard2-benchmark-attestation/"
            "leopard2_benchmark_source_attestation.h",
            "generated_compile_input", header, mtime_ns=50),
        "content": text_identity(header.decode("ascii")),
        "source_commit": source_commit, "source_tree": source_tree,
        "source_tracked_dirty": False,
    }
    benchmark_object = artifact(
        "/fixture/candidate-build/CMakeFiles/bench.o", "object_file",
        mtime_ns=100)
    validate_generated_attestation(
        attestation, "/fixture/candidate-build", benchmark_object,
        source_commit, source_tree)
    bad_attestation = clone(attestation)
    bad_attestation["source_tracked_dirty"] = 0
    expect_failure(lambda: validate_generated_attestation(
        bad_attestation, "/fixture/candidate-build", benchmark_object,
        source_commit, source_tree), "integer source-dirty attestation")

    blob = b"production selector fixture\n"
    blob_id = sha1_object("blob", blob)
    tree_bytes = b"100644 leopard2.cpp\0" + bytes.fromhex(blob_id)
    tree_id = sha1_object("tree", tree_bytes)
    commit_bytes = (
        f"tree {tree_id}\nauthor Fixture <fixture@example.invalid> 0 +0000\n"
        f"committer Fixture <fixture@example.invalid> 0 +0000\n\nfixture\n"
    ).encode("ascii")
    commit_id = sha1_object("commit", commit_bytes)
    tracked = [{
        "path": "leopard2.cpp", "git_mode": "100644",
        "git_type": "blob", "object_id": blob_id, "kind": "regular",
    }]
    listing = (
        f"100644 blob {blob_id}\tleopard2.cpp\0").encode("utf-8")
    stage = (
        f"100644 {blob_id} 0\tleopard2.cpp\0").encode("utf-8")
    flags = b"H leopard2.cpp\0"
    byte_identity = lambda data: {
        "size": len(data), "sha256": sha256_bytes(data)}
    object_identity = lambda data, object_id: {
        "encoding": "base64", "size": len(data),
        "sha256": sha256_bytes(data), "object_id": object_id,
        "base64": base64.b64encode(data).decode("ascii"),
    }
    git_source = {
        "path": "/usr/bin/git", "size": 7, "mode": 0o755,
        "sha256": sha256_bytes(b"git-bin"),
    }
    git_capture = {
        "schema": GIT_CAPTURE_SCHEMA, "path": "/fixture/source",
        "head": commit_id, "tree": tree_id, "detached": True,
        "head_ref": None, "superproject_worktree": None,
        "tracked_tree_listing_sha256": sha256_bytes(listing),
        "tracked_status": "clean",
        "commit_object": object_identity(commit_bytes, commit_id),
        "tree_objects": [object_identity(tree_bytes, tree_id)],
        "git_executable": {
            "source": git_source,
            "sealed": {
                "protocol": GIT_EXECUTABLE_PROTOCOL, "size": 7,
                "mode": 0o755, "sha256": git_source["sha256"],
                "seals": REQUIRED_EXECUTABLE_SEALS,
                "source_sha256": git_source["sha256"],
            },
        },
        "git_metadata": {
            "layout": "ordinary", "gitdir": "/fixture/source/.git",
            "commondir": "/fixture/source/.git",
            "guarded_components": ["/fixture/source/.git"],
            "guard_policy": GIT_METADATA_GUARD_POLICY,
            "guarded_file_count": 1,
            "guarded_files_sha256": sha256_bytes(b"guard"),
        },
        "worktree_guard_policy": GIT_WORKTREE_GUARD_POLICY,
        "config": byte_identity(b""),
        "index": {
            "entry_count": 1, "stage": byte_identity(stage),
            "flags_v": byte_identity(flags), "flags_f": byte_identity(flags),
        },
        "tracked_files": tracked,
        "tracked_files_sha256": sha256_bytes(canonical_bytes(tracked)),
        "submodules": [],
    }
    validate_git_capture(git_capture, "/fixture/source", commit_id, True,
                         "self-test")
    bad_git = clone(git_capture)
    bad_git["tree_objects"][0]["base64"] = \
        base64.b64encode(tree_bytes + b"x").decode("ascii")
    expect_failure(lambda: validate_git_capture(
        bad_git, "/fixture/source", commit_id, True, "self-test"),
                   "altered Git tree-object bytes")
    bad_git = clone(git_capture)
    bad_git["git_executable"]["source"]["mode"] = -1
    bad_git["git_executable"]["sealed"]["mode"] = -1
    expect_failure(lambda: validate_git_capture(
        bad_git, "/fixture/source", commit_id, True, "self-test"),
                   "negative Git executable modes")
    bad_git = clone(git_capture)
    bad_git["git_executable"]["sealed"]["mode"] = 0o751
    expect_failure(lambda: validate_git_capture(
        bad_git, "/fixture/source", commit_id, True, "self-test"),
                   "Git sealed/source mode mismatch")
    bad_git = clone(git_capture)
    bad_git["git_executable"]["sealed"]["seals"] = -1
    expect_failure(lambda: validate_git_capture(
        bad_git, "/fixture/source", commit_id, True, "self-test"),
                   "negative Git executable seal mask")
    bad_git = clone(git_capture)
    bad_git["git_executable"]["source"]["path"] = "/tmp/git"
    expect_failure(lambda: validate_git_capture(
        bad_git, "/fixture/source", commit_id, True, "self-test"),
                   "substituted Git acquisition executable")

    sealed_source = artifact(
        "/fixture/executable", "executable", b"\x7fELFfixture",
        mode=0o755)
    sealed_snapshot = {
        "protocol": SEALED_EXECUTABLE_PROTOCOL,
        "size": sealed_source["size"], "mode": sealed_source["mode"],
        "sha256": sealed_source["sha256"],
        "seals": REQUIRED_EXECUTABLE_SEALS, "elf": True,
    }
    validate_sealed_snapshot(sealed_snapshot, sealed_source,
                             "self-test timed executable")
    bad_snapshot = clone(sealed_snapshot)
    bad_snapshot["seals"] = -1
    expect_failure(lambda: validate_sealed_snapshot(
        bad_snapshot, sealed_source, "self-test timed executable"),
                   "negative timed executable seal mask")
    bad_snapshot = clone(sealed_snapshot)
    bad_snapshot["mode"] = 0o751
    expect_failure(lambda: validate_sealed_snapshot(
        bad_snapshot, sealed_source, "self-test timed executable"),
                   "timed sealed/source mode mismatch")

    def cpu_policy(cpu: int, sibling: int) -> dict[str, Any]:
        cache = {
            "index": 0, "level": "1", "type": "Data", "size": "32K",
            "coherency_line_size": "64", "number_of_sets": "64",
            "ways_of_associativity": "8", "physical_line_partition": "1",
            "shared_cpu_list": "52,116", "shared_cpu_map": "1",
            "allocation_policy": "WriteAllocate", "write_policy": "WriteBack",
        }
        return {
            "cpu": cpu, "online": "1",
            "cpuinfo": {
                "processor": str(cpu), "vendor_id": "AuthenticAMD",
                "cpu family": "26", "model": "8", "model name": "fixture",
                "flags": "fpu avx2 gfni",
            },
            "topology": {
                "core_id": "1", "physical_package_id": "0", "die_id": "0",
                "cluster_id": None, "thread_siblings_list": "52,116",
                "core_siblings_list": "0,52,116",
            },
            "frequency_policy": {
                "scaling_driver": None, "scaling_governor": None,
                "energy_performance_preference": None,
                "scaling_min_freq": None, "scaling_max_freq": None,
                "cpuinfo_min_freq": None, "cpuinfo_max_freq": None,
            },
            "cache_hierarchy": [cache], "cache_index_inventory": ["index0"],
            "cache_directory_inventory_text": text_identity("index0\n"),
            "numa_nodes": [0], "numa_node_inventory": ["node0"],
            "core_class": {"core_type": None, "cpu_capacity": None},
        }

    host = {
        "system": {
            "system": "Linux", "node": "fixture", "release": "1",
            "version": "1", "machine": "x86_64", "python": "3.13",
            "page_size": 4096,
        },
        "allowed_cpu_set_at_launch": [0, CPU, SIBLING],
        "online_cpu_set": [0, CPU, SIBLING],
        "online_cpu_list_text": text_identity("0,52,116\n"),
        "online_node_list_text": text_identity("0\n"),
        "benchmark_cpu": cpu_policy(CPU, SIBLING),
        "reserved_sibling": cpu_policy(SIBLING, CPU),
        "turbo_and_pstate": {
            "/sys/devices/system/cpu/intel_pstate/no_turbo": None,
            "/sys/devices/system/cpu/amd_pstate/status": None,
            "/sys/devices/system/cpu/cpufreq/boost": None,
        },
    }
    validate_host(host, "self-test", [0, CPU, SIBLING])
    bad_host = clone(host)
    bad_host["benchmark_cpu"]["topology"]["thread_siblings_list"] = "52"
    expect_failure(lambda: validate_host(
        bad_host, "self-test", [0, CPU, SIBLING]),
                   "host missing reserved SMT sibling")
    bad_host = clone(host)
    bad_host["online_cpu_set"][1] = float(CPU)
    expect_failure(lambda: validate_host(
        bad_host, "self-test", [0, CPU, SIBLING]),
                   "floating host online CPU")
    bad_host = clone(host)
    bad_host["benchmark_cpu"]["cpuinfo"]["flags"] = "fpu avx2"
    expect_failure(lambda: validate_host(
        bad_host, "self-test", [0, CPU, SIBLING]),
                   "host without retained GFNI feature flag")
    bad_host = clone(host)
    bad_host["system"]["machine"] = "aarch64"
    expect_failure(lambda: validate_host(
        bad_host, "self-test", [0, CPU, SIBLING]),
                   "non-x86 host for x86 GFNI campaign")

    with tempfile.TemporaryDirectory(prefix="leopard-v17-audit-self-test.") \
            as directory:
        root = Path(directory).resolve()
        source = root / "source"
        source.write_bytes(b"payload")
        hardlink = root / "hardlink"
        os.link(source, hardlink)
        symlink = root / "symlink"
        symlink.symlink_to(source)
        expect_failure(lambda: stable_file(source, 1), "oversized file")
        expect_failure(lambda: stable_file(hardlink, 64), "hard-linked file")
        expect_failure(lambda: stable_file(symlink, 64), "symbolic link")
        output = root / "audit.json"
        write_canonical(output, {"b": 2, "a": 1})
        require(output.read_bytes() == b'{"a":1,"b":2}\n',
                "canonical output writer failed")
        expect_failure(lambda: write_canonical(output, {}),
                       "audit output overwrite")
    print("v17 GFNI exact-main independent auditor self-test passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--manifest", type=Path)
    source.add_argument("--report", type=Path,
                        help="compatibility alias for --manifest")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    options = parser.parse_args()
    manifest = options.manifest if options.manifest is not None else options.report
    if options.self_test:
        require(manifest is None and options.output is None,
                "--self-test cannot be combined with replay paths")
        self_test()
        return 0
    require(manifest is not None and options.output is not None,
            "--manifest/--report and --output are required")
    require(manifest.is_absolute() and options.output.is_absolute(),
            "audit paths must be absolute")
    require(manifest != options.output,
            "audit output aliases the manifest")
    result = replay(manifest)
    write_canonical(options.output, result)
    print(json.dumps({
        "audit_passed": True,
        "benchmark_executed": False,
        "manifest_sha256": result["manifest"]["sha256"],
        "output": str(options.output),
        "timing_performed": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AuditError, UnicodeError, json.JSONDecodeError, OSError) as error:
        print(f"v17 independent audit failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
