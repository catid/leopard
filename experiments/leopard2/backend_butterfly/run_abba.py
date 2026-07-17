#!/usr/bin/env python3
"""Run and replay portable, fail-closed Leopard2 butterfly ABBA evidence."""

from __future__ import print_function

import argparse
import base64
import binascii
import copy
import hashlib
import json
import math
import os
import platform
import re
import shlex
import shutil
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


SCHEMA = "leopard2-backend-butterfly-abba/v5"
RAW_SCHEMA = "leopard2-backend-butterfly-raw/v1"
RESERVATION_SCHEMA = "leopard2-cpu-reservation/v1"
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
    cell("high-gf8", 240, 16, "high", "gf8", 65536, 4, "target"),
    cell("low-gf8-64b", 32, 224, "low", "gf8", 64, 16, "neighbor"),
    cell("low-gf8-tail", 32, 224, "low", "gf8", 65, 16, "neighbor"),
    cell("low-gf8-1k", 32, 224, "low", "gf8", 1024, 16, "neighbor"),
    cell("low-gf8", 32, 224, "low", "gf8", 65536, 16, "target"),
    cell("high-gf16-64b", 1000, 200, "high", "gf16", 64, 8, "neighbor"),
    cell("high-gf16-tail", 1000, 200, "high", "gf16", 66, 8, "neighbor"),
    cell("high-gf16-1k", 1000, 200, "high", "gf16", 1024, 8, "neighbor"),
    cell("high-gf16", 1000, 200, "high", "gf16", 16384, 8, "target"),
    cell("low-gf16-64b", 128, 1024, "low", "gf16", 64, 16, "neighbor"),
    cell("low-gf16-tail", 128, 1024, "low", "gf16", 66, 16, "neighbor"),
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

BUILD_TRANSLATION_UNITS = (
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

MATRIX_SOURCE_FILES = (
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
    "tests/leopard2/legacy_golden_vectors.h",
    "tests/leopard2/test_api.cpp",
    "tests/leopard2/test_public_api_contract.cpp",
    "tests/leopard2/test_random.cpp",
    "tests/leopard2/test_locator.cpp",
    "tests/leopard2/test_boundaries.cpp",
    "tests/leopard2/test_active_lch.cpp",
    "tests/leopard2/test_gf16_tails.cpp",
    "tests/leopard2/test_gf16_padded_odd.cpp",
    "tests/leopard2/test_encoder_gf16_legacy_matrix.cpp",
    "tests/leopard2/test_low_gf16_direct_rows.cpp",
    "tests/leopard2/test_decode_high_acceptance.cpp",
    "tests/leopard2/test_decode_low_acceptance.cpp",
    "tests/leopard2/test_decode_plan_schedule.cpp",
    "tests/leopard2/test_direct_encode.cpp",
    "tests/leopard2/test_arbitrary_counts_acceptance.cpp",
    "tests/leopard2/test_max_counts.cpp",
    "tests/leopard2/test_encode_concurrency.cpp",
    "tests/leopard2/test_codec_options_abi.c",
    "tests/leopard2/test_transform_differential.cpp",
    "tests/leopard2/direct_oracle.cpp",
    "tests/leopard2/direct_oracle.h",
    "tests/leopard2/test_direct_oracle.cpp",
    "tests/leopard2/direct_repair.cpp",
    "tests/leopard2/direct_repair.h",
    "tests/leopard2/test_direct_repair.cpp",
    "tests/leopard2/fuzz_api.cpp",
    "tests/leopard2/fuzz_replay.cpp",
    "tests/cmake/test_cuda_optional.cmake",
    "cmake/leopardConfig.cmake.in",
    "tools/check_leopard2_portable_isa.sh",
    "tools/leopard2_backend_matrix.py",
)

MATRIX_COMPARE_TESTS = (
    "direct_oracle", "backend_ops", "context_backends", "r1_xor",
    "legacy_golden", "api",
    "public_api_contract", "random", "locator", "active_lch", "gf16_tails",
    "gf16_padded_odd", "gf16_legacy_encoder_matrix",
    "low_gf16_direct_rows", "decode_high_acceptance",
    "decode_low_acceptance", "decode_plan_schedule", "direct_encode",
    "arbitrary_counts_acceptance", "max_counts", "encode_concurrency",
    "codec_options_abi", "direct_repair", "boundaries",
    "transform_differential", "fuzz_smoke",
)

MATRIX_BACKEND_FAILURE_TESTS = (
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

MATRIX_TEST_SPECS = {
    "direct_oracle": ("leopard2_direct_oracle_test", []),
    "backend_ops": ("leopard2_backend_ops_test", []),
    "context_backends": ("leopard2_context_backends_test", []),
    "r1_xor": ("leopard2_r1_xor_test", []),
    "legacy_golden": ("leopard2_legacy_golden_test", []),
    "api": ("leopard2_api_test", []),
    "public_api_contract": ("leopard2_public_api_contract_test", []),
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
    "fuzz_smoke": ("leopard2_fuzz_smoke", []),
}

MATRIX_BUILD_TARGETS = tuple(MATRIX_TEST_SPECS[name][0]
                             for name in MATRIX_COMPARE_TESTS)
MATRIX_BUILD_TARGETS = (
    MATRIX_BUILD_TARGETS[:2] + ("leopard2_backend_failures_test",) +
    MATRIX_BUILD_TARGETS[2:])

MATRIX_BUILD_CACHE_KEYS = (
    "CMAKE_BUILD_TYPE", "CMAKE_GENERATOR", "CMAKE_C_FLAGS",
    "CMAKE_C_FLAGS_RELEASE", "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE",
    "CMAKE_EXE_LINKER_FLAGS", "CMAKE_EXE_LINKER_FLAGS_RELEASE",
    "CMAKE_STATIC_LINKER_FLAGS", "CMAKE_STATIC_LINKER_FLAGS_RELEASE",
    "ENABLE_OPENMP", "LEO2_BACKEND_VARIANT", "LEO2_BUILD_TESTS",
    "LEO2_BUILD_BENCHMARKS", "LEO2_BUILD_FUZZERS", "LEO2_ENABLE_CUDA",
    "CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER",
)

MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS = {
    "Leopard2Backend.cpp": 1, "Leopard2BackendAVX2.cpp": 1,
    "Leopard2BackendSSSE3.cpp": 1, "Leopard2BackendScalar.cpp": 1,
    "Leopard2CpuFeatures.cpp": 1, "Leopard2Plan.cpp": 1,
    "LeopardCommon.cpp": 1, "LeopardFF16.cpp": 1, "LeopardFF8.cpp": 1,
    "leopard.cpp": 1, "leopard2.cpp": 1,
    "tests/leopard2/direct_oracle.cpp": 13,
    "tests/leopard2/direct_repair.cpp": 1,
    "tests/leopard2/fuzz_api.cpp": 1, "tests/leopard2/fuzz_replay.cpp": 1,
    "tests/leopard2/test_active_lch.cpp": 1,
    "tests/leopard2/test_api.cpp": 1,
    "tests/leopard2/test_arbitrary_counts_acceptance.cpp": 1,
    "tests/leopard2/test_backend_failures.cpp": 1,
    "tests/leopard2/test_backend_ops.cpp": 1,
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
    "tests/leopard2/test_gf16_padded_odd.cpp": 1,
    "tests/leopard2/test_gf16_tails.cpp": 1,
    "tests/leopard2/test_legacy_golden.cpp": 1,
    "tests/leopard2/test_locator.cpp": 1,
    "tests/leopard2/test_low_gf16_direct_rows.cpp": 1,
    "tests/leopard2/test_max_counts.cpp": 1,
    "tests/leopard2/test_public_api_contract.cpp": 1,
    "tests/leopard2/test_random.cpp": 1,
    "tests/leopard2/test_transform_differential.cpp": 1,
}

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

CONFIGURED_TRANSLATION_UNITS = BUILD_TRANSLATION_UNITS + (
    "tests/benchmark.cpp",
    "tests/experiments.cpp",
)

RELEVANT_TARGETS = (
    "libleopard",
    "leopard2_backend_ssse3",
    "leopard2_backend_avx2",
    "bench_leopard2",
)

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


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False).encode("ascii")


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256_file(path):
    return sha256_bytes(Path(path).read_bytes())


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


def benchmark_arguments(item):
    return [
        "--k", str(item["K"]), "--r", str(item["R"]),
        "--profile", item["profile"], "--field", item["field"],
        "--bytes", str(item["bytes"]), "--loss", str(item["loss"]),
        "--batch", "1", "--reuse", "8", "--iterations", "7",
        "--warmup", "3", "--threads", "1", "--backend", "auto",
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


def check_raw(raw, item, label, missing_indices=None):
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
    require(set(parameters) == required_parameter_keys,
            label + " parameter key set")
    expected = {
        "K": item["K"], "R": item["R"],
        "requested_profile": item["resolved_profile"],
        "requested_field": item["field"], "requested_backend": "auto",
        "force_generic_decode": False, "force_specialized_decode": False,
        "shard_bytes": item["bytes"], "loss_count": item["loss"],
        "batch": 1, "reuse": 8, "iterations": 7, "warmup": 3,
        "thread_count": 1, "seed": 42,
    }
    for key, expected_value in expected.items():
        require(parameters[key] == expected_value,
                "{} parameter {} expected {!r}, got {!r}".format(
                    label, key, expected_value, parameters[key]))
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
        "backend": "avx2", "thread_count": 1,
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

    The benchmark does not retain the seven underlying samples in v5, so this
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


def summarize(entries, raw_by_name):
    summary = []
    for item in CELLS:
        cell_result = {
            "name": item["name"],
            "kind": item["kind"],
            "minimum_speedup_percent": item["minimum_speedup_percent"],
            "parameters": {
                "K": item["K"], "R": item["R"],
                "profile": item["resolved_profile"], "field": item["field"],
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


def file_api_targets(source_root, build_root):
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
    require(set(RELEVANT_TARGETS) <= set(documents),
            "CMake target graph omits evidence target")

    configured_units = set()
    for document in documents.values():
        for source in document.get("sources", []):
            if "compileGroupIndex" not in source:
                continue
            path = Path(source["path"])
            if not path.is_absolute():
                path = source_root / path
            require(relative_to(path, source_root) is not None,
                    "configured translation unit escapes source root: " + str(path))
            configured_units.add(path.resolve().relative_to(source_root).as_posix())
    require(configured_units == set(CONFIGURED_TRANSLATION_UNITS),
            "configured translation-unit set mismatch: missing={} extra={}".format(
                sorted(set(CONFIGURED_TRANSLATION_UNITS) - configured_units),
                sorted(configured_units - set(CONFIGURED_TRANSLATION_UNITS))))

    targets = {}
    artifact_paths = {}
    relevant_units = set()
    for name in RELEVANT_TARGETS:
        document = documents[name]
        artifacts = document.get("artifacts", [])
        require(isinstance(artifacts, list) and len(artifacts) == 1,
                "evidence target must have exactly one artifact: " + name)
        artifact = Path(artifacts[0].get("path", ""))
        if not artifact.is_absolute():
            artifact = build_root / artifact
        artifact = artifact.resolve()
        require(relative_to(artifact, build_root) is not None,
                "target artifact escapes fresh build: " + name)
        artifact_paths[name] = artifact
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
            "artifact": tagged_path(artifact, source_root, build_root),
            "dependencies": dependencies,
            "sources": sorted(sources, key=lambda value: (
                value["path"], value["compiled"])),
            "link": normalized_link,
        }
    require(relevant_units == set(BUILD_TRANSLATION_UNITS),
            "evidence target translation-unit closure mismatch")
    require(targets["bench_leopard2"]["type"] == "EXECUTABLE" and
            targets["bench_leopard2"]["dependencies"] == ["libleopard"],
            "benchmark target dependency identity")
    require(targets["libleopard"]["type"] == "STATIC_LIBRARY" and
            targets["libleopard"]["dependencies"] ==
            ["leopard2_backend_avx2", "leopard2_backend_ssse3"],
            "library target dependency identity")
    require(all(targets[name]["type"] == "OBJECT_LIBRARY" and
                targets[name]["dependencies"] == []
                for name in ("leopard2_backend_ssse3", "leopard2_backend_avx2")),
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
    artifacts, target_graph = file_api_targets(source_root, build_root)
    argv = [str(cmake), "--build", str(build_root), "--parallel", str(jobs),
            "--target", "libleopard", "bench_leopard2"]
    logical_argv = ["@tool/cmake", "--build", "@build", "--parallel",
                    str(jobs), "--target", "libleopard", "bench_leopard2"]
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
        "library": artifacts["libleopard"],
        "binary": artifacts["bench_leopard2"],
        "target_graph": target_graph,
        "rebuild": rebuild,
    }


def self_test_rebuild_record(jobs, cache_path):
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
    build_argv = ["@tool/cmake", "--build", "@build", "--parallel",
                  str(jobs), "--target", "libleopard", "bench_leopard2"]
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
                         library, binary):
    library_recipe_path = build_root / "CMakeFiles/libleopard.dir/link.txt"
    benchmark_recipe_path = build_root / "CMakeFiles/bench_leopard2.dir/link.txt"
    require(library_recipe_path.is_file() and benchmark_recipe_path.is_file(),
            "evidence collection requires CMake literal link.txt recipes")
    library_lines = [shlex.split(value) for value in
                     library_recipe_path.read_text(encoding="utf-8").splitlines()
                     if value.strip()]
    benchmark_lines = [shlex.split(value) for value in
                       benchmark_recipe_path.read_text(encoding="utf-8").splitlines()
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
    recipe_objects = {resolve_recipe_path(value, build_root)
                      for value in archive[3:]}
    require(len(recipe_objects) == len(archive[3:]) and
            recipe_objects == expected_objects,
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
            external_inputs.append({"path": tagged_path(path, source_root, build_root),
                                    "sha256": sha256_file(path)})
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
        "library_recipe_sha256": sha256_file(library_recipe_path),
        "benchmark_recipe_sha256": sha256_file(benchmark_recipe_path),
    }
    recipes["digest"] = sha256_bytes(canonical_bytes(recipes))
    return recipes, expected_objects


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
    require(set(listing) == set(expected_by_name),
            "archive member set differs from target object recipe")
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


def build_record(source_root, source_identity, compile_commands, cmake_cache,
                 library, binary, rebuild, target_graph):
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
    require(set(all_by_file) == set(CONFIGURED_TRANSLATION_UNITS),
            "configured compile-command set mismatch: missing={} extra={}".format(
                sorted(set(CONFIGURED_TRANSLATION_UNITS) - set(all_by_file)),
                sorted(set(all_by_file) - set(CONFIGURED_TRANSLATION_UNITS))))
    by_file = {relative: all_by_file[relative]
               for relative in BUILD_TRANSLATION_UNITS}
    require(set(by_file) == set(BUILD_TRANSLATION_UNITS),
            "compile command translation-unit set mismatch: missing={} extra={}".format(
                sorted(set(BUILD_TRANSLATION_UNITS) - set(by_file)),
                sorted(set(by_file) - set(BUILD_TRANSLATION_UNITS))))
    require(all(value[0]["output_sha256"] is not None for value in by_file.values()),
            "evidence target compile output is missing")
    commands = {relative: by_file[relative][0]
                for relative in sorted(by_file)}
    recipes, expected_objects = literal_link_recipes(
        source_root, build_root, tool_paths, by_file, library, binary)
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
        "configured_translation_units": sorted(all_by_file),
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


def validate_build_record(record, repo):
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
        require(set(tool) == {"logical_name", "cache_key", "basename",
                             "binary_sha256", "version", "version_sha256"},
                "tool identity keys: " + logical_name)
        require(tool["logical_name"] == logical_name and
                tool["cache_key"] == cache_key and
                isinstance(tool["basename"], str) and tool["basename"],
                "tool logical identity: " + logical_name)
        require(configuration[cache_key] == "@tool/" + logical_name,
                "configuration/tool binding: " + logical_name)
        require(sha256_bytes(tool["version"].encode("utf-8")) ==
                tool["version_sha256"], "tool version digest: " + logical_name)
        require(re.fullmatch(r"[0-9a-f]{64}", tool["binary_sha256"] or "") is not None,
                "tool binary digest: " + logical_name)
    units = record["translation_units"]
    require(set(units) == set(BUILD_TRANSLATION_UNITS),
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
            sorted(CONFIGURED_TRANSLATION_UNITS),
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
    require(set(by_unit) == set(BUILD_TRANSLATION_UNITS),
            "dependency translation-unit set")
    for relative, dependencies in by_unit.items():
        require(isinstance(dependencies, list) and dependencies == sorted(set(dependencies)) and
                "@source/" + relative in dependencies and
                all(path in manifest_paths for path in dependencies),
                "dependency list identity: " + relative)
    require(isinstance(closure["file_count"], int) and
            closure["file_count"] >= len(BUILD_TRANSLATION_UNITS),
            "dependency-closure file count")
    require(isinstance(closure["source_file_count"], int) and
            len(BUILD_TRANSLATION_UNITS) <= closure["source_file_count"] <=
            closure["file_count"], "dependency-closure source count")

    target_graph = record["target_graph"]
    require(set(target_graph) == {"configured_translation_units", "targets",
                                 "index_sha256", "codemodel_sha256", "digest"},
            "target graph key set")
    graph_payload = dict(target_graph)
    graph_digest = graph_payload.pop("digest")
    require(sha256_bytes(canonical_bytes(graph_payload)) == graph_digest and
            target_graph["configured_translation_units"] ==
            sorted(CONFIGURED_TRANSLATION_UNITS) and
            set(target_graph["targets"]) == set(RELEVANT_TARGETS),
            "target graph identity")
    targets = target_graph["targets"]
    require(targets["bench_leopard2"]["type"] == "EXECUTABLE" and
            targets["bench_leopard2"]["artifact"] == "@build/bench_leopard2" and
            targets["bench_leopard2"]["dependencies"] == ["libleopard"] and
            targets["libleopard"]["type"] == "STATIC_LIBRARY" and
            targets["libleopard"]["artifact"] == "@build/liblibleopard.a" and
            targets["libleopard"]["dependencies"] ==
            ["leopard2_backend_avx2", "leopard2_backend_ssse3"],
            "target dependency graph")
    compiled_target_sources = set()
    for name, target in targets.items():
        require(set(target) == {"type", "artifact", "dependencies", "sources",
                               "link"} and
                target["artifact"].startswith("@build/") and
                isinstance(target["sources"], list),
                "target record: " + name)
        for source in target["sources"]:
            require(set(source) == {"path", "compiled"} and
                    isinstance(source["compiled"], bool),
                    "target source record: " + name)
            if source["compiled"]:
                require(source["path"].startswith("@source/"),
                        "compiled target source is not Git-backed")
                compiled_target_sources.add(source["path"][len("@source/"):])
    require(compiled_target_sources == set(BUILD_TRANSLATION_UNITS),
            "target compiled-source closure")

    recipes = record["link_recipes"]
    require(set(recipes) == {"library", "benchmark", "external_link_inputs",
                            "library_recipe_sha256", "benchmark_recipe_sha256",
                            "digest"}, "link recipe key set")
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
            library_command[2] == "@build/liblibleopard.a",
            "literal archive output recipe")
    expected_library_objects = {
        entry["output"] for relative, entry in units.items()
        if relative != "bench/leopard2/benchmark.cpp"}
    require(len(library_command[3:]) == len(set(library_command[3:])) and
            set(library_command[3:]) == expected_library_objects and
            recipes["library"][1] ==
            ["@tool/ranlib", "@build/liblibleopard.a"],
            "literal archive object/index closure")
    benchmark_command = recipes["benchmark"][0]
    expected_benchmark_object = units[
        "bench/leopard2/benchmark.cpp"]["output"]
    for external in recipes["external_link_inputs"]:
        require(set(external) == {"path", "sha256"} and
                external["path"].startswith("@external/") and
                re.fullmatch(r"[0-9a-f]{64}", external["sha256"]) is not None,
                "external link input identity")
    file_api_libraries = [
        value["fragment"]
        for value in targets["bench_leopard2"]["link"]["fragments"]
        if value["role"] == "libraries"]
    require(file_api_libraries.count("@build/liblibleopard.a") == 1,
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

    archive = record["archive"]
    require(set(archive) == {"members", "member_count", "digest"},
            "archive manifest key set")
    archive_payload = dict(archive)
    archive_digest = archive_payload.pop("digest")
    require(sha256_bytes(canonical_bytes(archive_payload)) == archive_digest and
            archive["member_count"] == len(archive["members"]) ==
            len(BUILD_TRANSLATION_UNITS) - 1,
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
            set(member_objects) == set(expected_archive_members),
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
            build_argv[5:] == ["--target", "libleopard", "bench_leopard2"] and
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


def validate_matrix_document(document, repo, candidate_commit):
    require(set(document) == {
        "c_compiler", "compiler", "generator", "jobs", "jobs_per_variant",
        "machine", "mismatches",
        "schema", "source_changed_during_run", "source_fingerprint",
        "status", "variant_workers", "variants"},
        "matrix top-level key set")
    require(document.get("schema") == "leopard2-backend-matrix/v1",
            "matrix schema")
    require(document.get("status") == "passed", "matrix status")
    require(document.get("source_changed_during_run") is False,
            "matrix source changed during run")
    fingerprint = document.get("source_fingerprint", {})
    files = fingerprint.get("files")
    require(isinstance(files, dict), "matrix source files")
    require(set(files) == set(MATRIX_SOURCE_FILES),
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
        require(value.get("schema") == "leopard2-backend-matrix/v1" and
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
        require(set(cache) == set(MATRIX_BUILD_CACHE_KEYS) and
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
                    command["file"] in MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS and
                    command["language"] ==
                    ("C" if command["file"].endswith(".c") else "CXX") and
                    isinstance(command["argv"], list) and command["argv"] and
                    command["argv"][0] ==
                    ("@tool/cc" if command["language"] == "C" else "@tool/cxx") and
                    command["argv"].count("@source/" + command["file"]) == 1,
                    "matrix compile command identity")
            counts[command["file"]] = counts.get(command["file"], 0) + 1
        require(counts == MATRIX_EXPECTED_COMPILE_SOURCE_COUNTS,
                "matrix compile source multiset")

        commands = value["commands"]
        require(len(commands) == 2 + len(MATRIX_COMPARE_TESTS) + 2 +
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
                list(MATRIX_BUILD_TARGETS),
                "matrix build argv: " + variant)
        require(value.get("source_fingerprint") == fingerprint["digest"],
                "matrix variant source mismatch")
        tests = value.get("tests")
        required_tests = set(MATRIX_COMPARE_TESTS) | {
            "backend_failures", "portable_isa"}
        if value.get("variant") == "auto":
            required_tests.add("cuda_optional")
        require(isinstance(tests, dict) and set(tests) == required_tests,
                "matrix test set: " + str(value.get("variant")))
        executables = identity["test_executables"]
        require(set(executables) == set(MATRIX_COMPARE_TESTS),
                "matrix executable set: " + variant)
        command_index = 2
        for test_name in MATRIX_COMPARE_TESTS:
            test = tests[test_name]
            validate_matrix_command(test, "test_" + test_name,
                                    {"executable_sha256"})
            require(canonical_bytes(test) == canonical_bytes(commands[command_index]),
                    "matrix test/command crosslink: " + test_name)
            command_index += 1
            executable = executables[test_name]
            target, arguments = MATRIX_TEST_SPECS[test_name]
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
                sorted(MATRIX_BACKEND_FAILURE_TESTS) and
            canonical_bytes(failures) ==
                canonical_bytes(commands[command_index]) and
            failures["cwd"] == source_root and len(failures["argv"]) >= 4 and
            Path(failures["argv"][0]).name == "taskset" and
            failures["argv"][1:3] == ["-c", str(value["pin_cpu"])] and
            Path(failures["argv"][3]).name == tools["ctest"]["basename"] and
            failures["argv"][4:] == [
                "--test-dir", build_root, "-C", "Release", "-R",
                "^leopard2_backend_failure_", "--output-on-failure"],
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
                    "-R", "^leopard2_portable_isa$", "--output-on-failure"],
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
                        "-C", "Release", "-R", "^leopard2_cuda_optional$",
                        "--output-on-failure"],
                    "matrix CUDA-optional CTest command")
            command_index += 1
        require(command_index == len(commands), "matrix command closure: " + variant)

    recomputed_mismatches = []
    auto_tests = by_variant["auto"]["tests"]
    for variant in sorted(set(by_variant) - {"auto"}):
        for test_name in MATRIX_COMPARE_TESTS:
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


def validate_declared_template_paths(args, build):
    cache = Path(getattr(args, build + "_cmake_cache")).resolve()
    root = cache.parent
    expected = {
        build + "_compile_commands": root / "compile_commands.json",
        build + "_library": root / "liblibleopard.a",
        build: root / "bench_leopard2",
    }
    for attribute, expected_path in expected.items():
        supplied = Path(getattr(args, attribute)).resolve()
        require(supplied == expected_path.resolve(),
                "{} path substitution: expected canonical target path {}".format(
                    attribute, expected_path))
        require(supplied.is_file(), attribute + " template artifact is missing")


def matrix_record(path, repo, candidate_commit):
    raw = Path(path).read_bytes()
    require(raw, "matrix artifact is empty")
    document = parse_json_bytes(raw, "backend matrix")
    validate_matrix_document(document, repo, candidate_commit)
    return raw, {
        "sha256": sha256_bytes(raw),
        "source_fingerprint": document["source_fingerprint"]["digest"],
        "variant_count": len(document["variants"]),
    }


def reservation_record(path, cpu, sibling):
    require(fcntl is not None,
            "evidence collection requires Linux/POSIX fcntl file locking")
    path = Path(path).resolve()
    require(path.is_file(), "reservation file missing")
    handle = path.open("r+")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as error:
        handle.close()
        raise EvidenceError("cannot acquire exclusive CPU reservation: {}".format(error))
    try:
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
        return handle, {"sha256": sha256_bytes(canonical), "document": document}
    except BaseException:
        handle.close()
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
    require(manifest.get("schema") == SCHEMA, "unsupported manifest schema")
    require(manifest.get("status") == "passed", "campaign status")
    campaign = manifest.get("campaign", {})
    require(set(campaign) == {
        "order_per_round", "rounds", "cell_count", "entry_count",
        "samples_per_invocation", "warmups_per_invocation", "reuse_per_sample",
        "batch", "threads", "seed", "target_threshold_percent",
        "neighbor_floor_percent", "self_test", "started_unix",
        "finished_unix"}, "campaign key set")
    is_self_test = campaign.get("self_test") is True
    require(not is_self_test or allow_self_test,
            "self-test manifest is not production evidence")
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
        validate_build_record(provenance["builds"][build], repo)
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
        matrix_document, repo, provenance["git"]["candidate"]["commit"])
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
            benchmark_arguments(item)
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
        parsed = check_raw(raw, item, name, missing_by_cell.get(item["name"]))
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
    recomputed = summarize(entries, raw_by_name)
    require(canonical_bytes(recomputed) == canonical_bytes(manifest["summary"]),
            "summary does not replay from raw evidence")
    for item in recomputed:
        require(item["encode"]["accepted"] and item["decode"]["accepted"],
                item["name"] + " violates its promotion/neighbor floor")
    if binaries is not None:
        for build in ("baseline", "candidate"):
            path = Path(binaries[build]).resolve()
            require(path.is_file(), build + " binary missing")
            require(sha256_file(path) == binary_hashes[build],
                    build + " supplied binary hash mismatch")
    return manifest


def run_campaign(args, repo, allow_dirty=False, self_test=False):
    require(hasattr(os, "sched_getaffinity") and hasattr(os, "sched_setaffinity"),
            "Linux scheduler affinity APIs are required")
    require(isinstance(args.build_jobs, int) and 1 <= args.build_jobs <= 128,
            "build jobs must be in [1,128]")
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
        validate_declared_template_paths(args, build)
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
                    args.build_jobs, getattr(args, build + "_cmake_cache")),
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
            isolated["baseline"]["target_graph"]),
        "candidate": build_record(
            source_roots["candidate"], git_records["candidate"],
            isolated["candidate"]["compile_commands"],
            isolated["candidate"]["cmake_cache"], isolated["candidate"]["library"],
            isolated["candidate"]["binary"], isolated["candidate"]["rebuild"],
            isolated["candidate"]["target_graph"]),
    }
    binaries = {build: isolated[build]["binary"].resolve()
                for build in ("baseline", "candidate")}
    matrix_raw, matrix_info = matrix_record(
        args.matrix, repo, commits["candidate"])
    initial_affinity = sorted(os.sched_getaffinity(0))
    require(args.cpu in initial_affinity and args.reserved_sibling in initial_affinity,
            "isolated core pair is not initially allowed")
    topology = sibling_topology(args.cpu)
    require(args.reserved_sibling in topology["cpus"] and
            args.reserved_sibling != args.cpu,
            "reserved CPU is not a distinct topology sibling")
    reservation_handle, reservation = reservation_record(
        args.reservation_file, args.cpu, args.reserved_sibling)
    try:
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
            actual_argv = [str(binaries[build])] + benchmark_arguments(item)
            logical_argv = [logical_executable(
                build, builds[build]["artifacts"]["benchmark_sha256"])] + \
                benchmark_arguments(item)
            command_record = {"affinity": enforced_affinity,
                              "argv": logical_argv,
                              "environment": child_environment}
            print("[{}/{}] {}".format(index, len(jobs), name),
                  file=sys.stderr, flush=True)
            try:
                completed = subprocess.run(
                    actual_argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    env=child_environment, timeout=args.timeout, check=False)
            except subprocess.TimeoutExpired as error:
                raise EvidenceError("benchmark timeout {}: {}".format(name, error))
            require(completed.returncode == 0,
                    "benchmark failed {} rc={}".format(name, completed.returncode))
            raw = parse_json_bytes(completed.stdout, "benchmark output " + name)
            parsed = check_raw(raw, item, name, missing_by_cell.get(item["name"]))
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
        campaign_end = time.time()
        require(reported_builds.get("baseline") == reported_builds.get("candidate"),
                "baseline/candidate reported compiler identity mismatch")
        require(sorted(os.sched_getaffinity(0)) == enforced_affinity,
                "runner affinity changed during campaign")
        post_frequency = current_frequency(args.cpu)
        require(power_state(args.cpu) == static_power,
                "CPU power/governor state changed during campaign")
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
        summary = summarize(entries, raw_by_name)
        manifest = {
            "schema": SCHEMA,
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
        reservation_handle.close()
        os.sched_setaffinity(0, set(initial_affinity))
        shadow_context.cleanup()


def git_file_hashes(repo, commit, relatives):
    files = {relative: sha256_bytes(git_blob(repo, commit, relative))
             for relative in relatives}
    return {"digest": sha256_bytes(canonical_bytes(files)), "files": files}


def write_mock(path, factor):
    source = """#!/usr/bin/env python3
import json,sys
a=sys.argv[1:]
def v(name): return a[a.index(name)+1]
k=int(v('--k')); r=int(v('--r')); profile=v('--profile'); field=v('--field')
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
out={'schema':'leopard2-benchmark-v1','build':{'compiler':'mock','compiler_version':'1','cplusplus':201103},'parameters':{'K':k,'R':r,'requested_profile':resolved_profile,'requested_field':field,'requested_backend':'auto','force_generic_decode':False,'force_specialized_decode':False,'shard_bytes':byte_count,'loss_count':loss,'missing_original_indices':missing,'batch':int(v('--batch')),'reuse':int(v('--reuse')),'iterations':int(v('--iterations')),'warmup':int(v('--warmup')),'thread_count':int(v('--threads')),'seed':int(v('--seed'))},'resolved':{'profile':resolved_profile,'field':field,'backend':'avx2','thread_count':1,'parent_count':parent,'padded_side':padded},'correctness':{'leopard2_round_trip':True,'legacy_comparison':legacy_comparison},'memory':{'scratch_alignment':64,'encode_scratch_bytes_per_stripe':64,'decode_scratch_bytes_per_stripe':128,'encode_scratch_bytes_batch':64,'decode_scratch_bytes_batch':128},'metrics':{'codec_setup':setup(1.0),'encode_execution':encode,'decode_plan_setup':setup(2.0),'decode_execution':decode,'decode_amortized_at_reuse':{'reuse_count':8,'derived_median_us_per_batch_call':base*2*factor+.25,'offered_received_GB_per_s':1.0/(base*2*factor+.25),'repaired_output_GB_per_s':2.0/(base*2*factor+.25)},'rate_semantics':'offered_received counts all non-null shard pointers supplied; a plan may read a deterministic subset'},'legacy':{'available':legacy_comparison is not None,'unavailable_reason':legacy_reason,'codec_setup':None,'decode_timing_includes_setup':True,'encode_execution':rated(base*1.1,'input_GB_per_s','parity_output_GB_per_s') if legacy_comparison else None,'decode_including_setup':rated(base*2.2,'offered_received_GB_per_s','repaired_output_GB_per_s') if legacy_comparison else None}}
print(json.dumps(out,sort_keys=True,allow_nan=False))
""" % factor
    path.write_text(source, encoding="utf-8")
    path.chmod(0o755)


def write_self_test_build_files(root, source_root, binary):
    compiler = Path(shutil.which("c++") or "/usr/bin/c++").resolve()
    cc = Path(shutil.which("cc") or "/usr/bin/cc").resolve()
    ar = Path(shutil.which("ar") or "/usr/bin/ar").resolve()
    ranlib = Path(shutil.which("ranlib") or "/usr/bin/ranlib").resolve()
    compile_commands = []
    outputs = {}
    for index, relative in enumerate(CONFIGURED_TRANSLATION_UNITS):
        output = root / ("unit-{}.o".format(index))
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
    library = root / "liblibleopard.a"
    archive_objects = [outputs[value] for value in BUILD_TRANSLATION_UNITS
                       if value != "bench/leopard2/benchmark.cpp"]
    command_output([str(ar), "rcs", str(library)] +
                   [str(value) for value in archive_objects], root,
                   "self-test archive")
    library_link = root / "CMakeFiles/libleopard.dir/link.txt"
    benchmark_link = root / "CMakeFiles/bench_leopard2.dir/link.txt"
    library_link.parent.mkdir(parents=True)
    benchmark_link.parent.mkdir(parents=True)
    library_link.write_text(
        "{} qc {} {}\n{} {}\n".format(
            shlex.quote(str(ar)), shlex.quote(str(library)),
            " ".join(shlex.quote(str(value)) for value in archive_objects),
            shlex.quote(str(ranlib)), shlex.quote(str(library))),
        encoding="utf-8")
    benchmark_link.write_text(
        "{} {} -o {} {}\n".format(
            shlex.quote(str(compiler)),
            shlex.quote(str(outputs["bench/leopard2/benchmark.cpp"])),
            shlex.quote(str(binary)), shlex.quote(str(library))),
        encoding="utf-8")

    core = set(BUILD_TRANSLATION_UNITS) - {
        "Leopard2BackendSSSE3.cpp", "Leopard2BackendAVX2.cpp",
        "bench/leopard2/benchmark.cpp"}
    def target(target_type, artifact, dependencies, compiled):
        return {"type": target_type, "artifact": "@build/" + artifact.name,
                "dependencies": dependencies,
                "sources": [{"path": "@source/" + value, "compiled": True}
                            for value in sorted(compiled)],
                "link": None}
    targets = {
        "libleopard": target("STATIC_LIBRARY", library,
                             ["leopard2_backend_avx2", "leopard2_backend_ssse3"], core),
        "leopard2_backend_ssse3": target(
            "OBJECT_LIBRARY", outputs["Leopard2BackendSSSE3.cpp"], [],
            {"Leopard2BackendSSSE3.cpp"}),
        "leopard2_backend_avx2": target(
            "OBJECT_LIBRARY", outputs["Leopard2BackendAVX2.cpp"], [],
            {"Leopard2BackendAVX2.cpp"}),
        "bench_leopard2": target("EXECUTABLE", binary, ["libleopard"],
                                  {"bench/leopard2/benchmark.cpp"}),
    }
    targets["bench_leopard2"]["link"] = {
        "language": "CXX",
        "fragments": [{"role": "libraries",
                       "fragment": "@build/liblibleopard.a"}]}
    graph = {"configured_translation_units": sorted(CONFIGURED_TRANSLATION_UNITS),
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


def recompute_self_test_summary(manifest, bundle):
    raw_by_name = {}
    entries = {value["name"]: value for value in manifest["entries"]}
    missing_by_cell = {}
    for name, item, _, _, _ in expected_jobs():
        stdout, _ = decode_raw_record(bundle["raw"][name], name)
        parsed = check_raw(
            parse_json_bytes(stdout, "self-test summary " + name), item, name,
            missing_by_cell.get(item["name"]))
        missing_by_cell.setdefault(item["name"], parsed[2])
        raw_by_name[name] = parsed[:2]
    manifest["summary"] = summarize(list(entries.values()), raw_by_name)


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
                for name in MATRIX_COMPARE_TESTS
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
            for name in MATRIX_COMPARE_TESTS:
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
                 "^leopard2_backend_failure_", "--output-on-failure"],
                {"ctest_executed": True,
                 "ctest_executed_tests":
                     sorted(MATRIX_BACKEND_FAILURE_TESTS)})
            tests["backend_failures"] = failures
            commands.append(failures)
            portable = matrix_command(
                "test_portable_isa",
                [matrix_ctest, "--test-dir", matrix_build_root, "-C", "Release",
                 "-R", "^leopard2_portable_isa$", "--output-on-failure"],
                {"ctest_executed": True})
            tests["portable_isa"] = portable
            commands.append(portable)
            if value == "auto":
                cuda = matrix_command(
                    "test_cuda_optional",
                    [matrix_ctest, "--test-dir", matrix_build_root, "-C", "Release",
                     "-R", "^leopard2_cuda_optional$", "--output-on-failure"],
                    {"ctest_executed": True})
                tests["cuda_optional"] = cuda
                commands.append(cuda)
            variants.append({
                "configuration_id": configuration_id, "resumed": False,
                "variant": value, "status": "passed", "reason": "",
                "schema": "leopard2-backend-matrix/v1",
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
        atomic_json(matrix, {
            "c_compiler": matrix_c_compiler,
            "compiler": matrix_cxx_compiler,
            "generator": matrix_generator,
            "jobs": 1, "jobs_per_variant": 1,
            "machine": matrix_machine,
            "schema": "leopard2-backend-matrix/v1", "status": "passed",
            "source_changed_during_run": False, "mismatches": [],
            "source_fingerprint": fingerprint, "variant_workers": 1,
            "variants": variants,
        })
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
        )
        args.baseline_self_test_artifacts = baseline_build[3]
        args.candidate_self_test_artifacts = candidate_build[3]
        run_campaign(args, repo, allow_dirty=True, self_test=True)

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
        manifest_path = output / "abba_manifest.json"
        bundle_path = output / "abba_raw.json"
        manifest = read_json(manifest_path, "self-test manifest")
        bundle = read_json(bundle_path, "self-test raw bundle")
        validate = lambda mp, bp: validate_manifest(
            mp, repo, bp, None, None, allow_self_test=True)
        validate(manifest_path, bundle_path)

        mutations = []
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

        def mutate_matrix_test(m, b):
            def callback(document):
                document["variants"][0]["tests"]["backend_ops"]["returncode"] = 1
            mutate_embedded_matrix(m, b, callback)
        mutations.append(("coordinated matrix test", mutate_matrix_test))

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
        injected_archive = candidate_root / "injected-liblibleopard.a"
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
        mutation_count = len(mutations) + 6
    print("butterfly ABBA v5 self-test passed: path-independent replay + {} adversarial mutations".format(
        mutation_count))


def add_run_arguments(parser):
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
            print("butterfly ABBA v5 campaign passed: cells={} entries={}".format(
                len(CELLS), len(expected_jobs())))
        elif args.command == "verify":
            supplied = None
            if args.baseline is not None or args.candidate is not None:
                require(args.baseline is not None and args.candidate is not None,
                        "supply both binaries or neither")
                supplied = {"baseline": args.baseline, "candidate": args.candidate}
            validate_manifest(args.manifest, repo, args.raw_bundle,
                              supplied, args.matrix)
            print("butterfly ABBA v5 path-independent evidence replay passed")
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
