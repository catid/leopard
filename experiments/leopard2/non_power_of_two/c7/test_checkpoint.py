#!/usr/bin/env python3
"""Fail-closed validation for the retained C7 checkpoint."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import pathlib
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import unittest


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RESULTS = HERE / "results"
BACKENDS = ("scalar", "ssse3", "avx2", "auto")
BUILD_NAMES = (*BACKENDS, "asan-ubsan")
EXPECTED_RUNTIME = {
    "scalar": "scalar", "ssse3": "ssse3", "avx2": "avx2", "auto": "avx2"
}
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
GIT_SHA = re.compile(r"[0-9a-f]{40}\Z")
ARCHIVE_MEMBERS = (
    "leopard.cpp.o", "leopard2.cpp.o", "Leopard2Backend.cpp.o",
    "Leopard2BackendScalar.cpp.o", "Leopard2CpuFeatures.cpp.o",
    "Leopard2Plan.cpp.o", "LeopardCommon.cpp.o", "LeopardFF16.cpp.o",
    "LeopardFF8.cpp.o", "Leopard2BackendSSSE3.cpp.o",
    "Leopard2BackendAVX2.cpp.o",
)
SANITIZED_MEMBER_COUNTS = {
    "leopard.cpp.o": (13, 7),
    "leopard2.cpp.o": (140, 15),
    "Leopard2Backend.cpp.o": (21, 8),
    "Leopard2BackendScalar.cpp.o": (11, 6),
    "Leopard2CpuFeatures.cpp.o": (9, 5),
    "Leopard2Plan.cpp.o": (7, 5),
    "LeopardCommon.cpp.o": (13, 6),
    "LeopardFF16.cpp.o": (31, 10),
    "LeopardFF8.cpp.o": (31, 9),
    "Leopard2BackendSSSE3.cpp.o": (14, 8),
    "Leopard2BackendAVX2.cpp.o": (16, 7),
}
CORE_SOURCES = {
    "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp", "LeopardCommon.cpp", "LeopardFF16.cpp",
    "LeopardFF8.cpp", "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
}
EXPECTED_CORRECTNESS = {
    "gf8_cases": 9,
    "gf16_cases": 5,
    "coefficients": 118717,
    "gf16_vandermonde_coefficients": 1500,
    "encode_executions": 117,
    "encode_symbol_comparisons": 1030423,
    "subset_encode_executions": 117,
    "decode_plans": 403,
    "decode_executions": 403,
    "decode_symbol_comparisons": 272487,
    "maximum_loss_plans": 117,
    "unavailable_parity_plans": 175,
    "no_loss_null_calls": 14,
    "parity_rebuilds": 403,
    "odd_gf16_rejections": 10,
    "overlap_rejections": 59,
    "parity_output_overlap_rejections": 13,
    "restored_output_overlap_rejections": 12,
    "restored_input_overlap_rejections": 20,
    "selected_parity_null_rejections": 14,
    "survivor_null_rejections": 6,
    "atomic_rejection_bytes_checked": 61570,
    "read_only_input_alias_calls": 13,
    "read_only_input_alias_symbol_comparisons": 2139,
    "hot_path_allocations": 0,
    "digest_fnv64": "0xec4179e9f2776a58",
}


def digest(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def resolve_record_path(text: str) -> pathlib.Path:
    path = pathlib.Path(text)
    return path if path.is_absolute() else ROOT / path


def validate_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or not SHA256.fullmatch(value):
        raise ValueError(f"{label} is not canonical SHA-256")
    return value


def validate_git_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or not GIT_SHA.fullmatch(value):
        raise ValueError(f"{label} is not canonical git SHA")
    return value


def validate_artifact(
    record: object, label: str, *, require_retained: bool,
) -> pathlib.Path:
    if not isinstance(record, dict) or set(record) != {"bytes", "path", "sha256"}:
        raise ValueError(f"{label} artifact schema changed")
    path_text = record["path"]
    if not isinstance(path_text, str) or not path_text:
        raise ValueError(f"{label} path is invalid")
    if type(record["bytes"]) is not int or record["bytes"] < 0:
        raise ValueError(f"{label} byte count is invalid")
    expected = validate_sha(record["sha256"], f"{label} hash")
    path = resolve_record_path(path_text)
    if path.is_file():
        if path.stat().st_size != record["bytes"] or digest(path) != expected:
            raise ValueError(f"{label} retained bytes disagree with manifest")
    elif require_retained:
        raise ValueError(f"{label} retained artifact is missing")
    return path


def expected_program(names: tuple[str, ...], *, resolve: bool) -> pathlib.Path:
    for name in names:
        found = shutil.which(name)
        if found:
            path = pathlib.Path(found).absolute()
            return path.resolve() if resolve else path
    raise ValueError(f"required validation program is absent: {names}")


def compiler_program(compiler: pathlib.Path, name: str) -> pathlib.Path:
    completed = subprocess.run(
        [str(compiler), f"-print-prog-name={name}"], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    value = completed.stdout.strip()
    if not value:
        raise ValueError(f"compiler did not resolve {name}")
    path = pathlib.Path(value)
    if path.is_absolute():
        return path.absolute()
    found = shutil.which(value)
    if not found:
        raise ValueError(f"compiler-selected tool is absent: {value}")
    return pathlib.Path(found).absolute()


def validate_program(
    record: object, label: str, expected_path: pathlib.Path,
) -> pathlib.Path:
    if not isinstance(record, dict) or set(record) != {
            "path", "sha256", "version"}:
        raise ValueError(f"{label} program schema changed")
    path_text = record["path"]
    if not isinstance(path_text, str) or not pathlib.Path(path_text).is_absolute():
        raise ValueError(f"{label} path is not absolute")
    path = pathlib.Path(path_text)
    if path != expected_path or not path.is_file():
        raise ValueError(f"{label} path does not match its fixed role")
    expected_hash = validate_sha(record["sha256"], f"{label} hash")
    if digest(path) != expected_hash:
        raise ValueError(f"{label} current executable bytes changed")
    version = subprocess.run(
        [str(path), "--version"], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    ).stdout
    if record["version"] != version or not version.strip():
        raise ValueError(f"{label} version output changed")
    return path


def validate_git_artifact(
    core_sha: str, record: dict, label: str,
) -> None:
    path = resolve_record_path(record["path"])
    try:
        relative_path = path.resolve().relative_to(ROOT).as_posix()
    except ValueError as error:
        raise ValueError(f"{label} escaped the repository") from error
    committed = subprocess.run(
        ["git", "show", f"{core_sha}:{relative_path}"], cwd=ROOT,
        check=False, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    if committed.returncode != 0 or hashlib.sha256(
            committed.stdout).hexdigest() != record["sha256"]:
        raise ValueError(f"{label} is not bound to the core commit")


def parse_cmake_cache(text: str) -> dict[str, str]:
    if "FORGED" in text or not text.endswith("\n"):
        raise ValueError("CMake cache is not retained raw output")
    values: dict[str, str] = {}
    for line in text.splitlines():
        match = re.fullmatch(r"([^#/:=][^:=]*):([^=]+)=(.*)", line)
        if not match:
            continue
        key, value = match.group(1), match.group(3)
        if key in values:
            raise ValueError(f"duplicate CMake cache key: {key}")
        values[key] = value
    return values


def validate_cmake_cache_text(
    text: str, *, sanitizer: bool, backend: str, build_dir: str,
    cmake: pathlib.Path, c_compiler: pathlib.Path,
    cxx_compiler: pathlib.Path, gmake: pathlib.Path, ar: pathlib.Path,
    ranlib: pathlib.Path, cmake_linker: pathlib.Path,
) -> None:
    values = parse_cmake_cache(text)
    flags = "-fsanitize=address,undefined -fno-omit-frame-pointer" if sanitizer else ""
    expected = {
        "CMAKE_AR": str(ar),
        "CMAKE_BUILD_TYPE": "Debug" if sanitizer else "Release",
        "CMAKE_COMMAND": str(cmake),
        "CMAKE_CXX_COMPILER": str(cxx_compiler),
        "CMAKE_CXX_FLAGS": flags,
        "CMAKE_C_COMPILER": str(c_compiler),
        "CMAKE_C_FLAGS": flags,
        "CMAKE_EXE_LINKER_FLAGS": flags,
        "CMAKE_GENERATOR": "Unix Makefiles",
        "CMAKE_HOME_DIRECTORY": str(ROOT),
        "CMAKE_LINKER": str(cmake_linker),
        "CMAKE_MAKE_PROGRAM": str(gmake),
        "CMAKE_RANLIB": str(ranlib),
        "ENABLE_OPENMP": "OFF" if sanitizer else "ON",
        "LEO2_BACKEND_VARIANT": backend,
        "LEO2_BUILD_BENCHMARKS": "OFF",
        "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_BUILD_TESTS": "OFF",
        "LEO2_ENABLE_CUDA": "OFF",
        "leopard_BINARY_DIR": str((ROOT / build_dir).resolve()),
    }
    for key, value in expected.items():
        if values.get(key) != value:
            raise ValueError(f"CMake cache role changed: {key}")


def validate_configure_log_text(
    text: str, *, sanitizer: bool, c_compiler: pathlib.Path,
    cxx_compiler: pathlib.Path, build_dir: str,
) -> None:
    if not text.endswith("\n") or "FORGED" in text or "CMake Error" in text or (
            "CMake Warning" in text):
        raise ValueError("configure log is not a clean raw success")
    lines = text.splitlines()
    if len(lines) < 35 or any(not line.startswith("-- ") for line in lines):
        raise ValueError("configure log structure changed")
    family = "Clang" if sanitizer else "GNU"
    c_version = subprocess.run(
        [str(c_compiler), "-dumpfullversion", "-dumpversion"], check=True,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout.strip()
    cxx_version = subprocess.run(
        [str(cxx_compiler), "-dumpfullversion", "-dumpversion"], check=True,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout.strip()
    required = {
        f"-- The C compiler identification is {family} {c_version}",
        f"-- The CXX compiler identification is {family} {cxx_version}",
        f"-- Check for working C compiler: {c_compiler} - skipped",
        f"-- Check for working CXX compiler: {cxx_compiler} - skipped",
        "-- Detecting C compile features - done",
        "-- Detecting CXX compile features - done",
        "-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success",
        "-- Found Threads: TRUE  ",
        "-- Performing Test LEO2_FLAG_MSSSE3 - Success",
        "-- Performing Test LEO2_FLAG_MAVX2 - Success",
        f"-- Build files have been written to: {(ROOT / build_dir).resolve()}",
    }
    if not required.issubset(set(lines)):
        raise ValueError("configure log lacks required semantic events")
    if sum(line.startswith("-- Configuring done (") for line in lines) != 1 or (
            sum(line.startswith("-- Generating done (") for line in lines) != 1):
        raise ValueError("configure/generate completion is ambiguous")
    if sum(line.startswith("-- Found OpenMP: TRUE") for line in lines) != 1:
        raise ValueError("OpenMP discovery record changed")


def validate_build_log_text(
    text: str, *, sanitizer: bool, cmake: pathlib.Path,
    cxx_compiler: pathlib.Path, gmake: pathlib.Path, ar: pathlib.Path,
    ranlib: pathlib.Path,
) -> None:
    lowered = text.lower()
    if (not text.endswith("\n") or "forged" in lowered or "error:" in lowered or
            "undefined reference" in lowered or "failed" in lowered):
        raise ValueError("core build log is not a clean raw success")
    built_sources = set(re.findall(
        r"Building CXX object [^\n ]*/([^/ ]+\.cpp)\.o", text))
    if built_sources != CORE_SOURCES:
        raise ValueError("core build log source set changed")
    for required in (
        str(cmake), str(cxx_compiler), str(gmake), str(ar), str(ranlib),
        "Built target leopard2_backend_ssse3",
        "Built target leopard2_backend_avx2",
        "Linking CXX static library liblibleopard.a",
        "Built target libleopard",
    ):
        if required not in text:
            raise ValueError(f"core build log did not execute {required}")
    has_sanitizer_flags = "-fsanitize=address,undefined" in text
    if has_sanitizer_flags is not sanitizer:
        raise ValueError("core build sanitizer command evidence changed")


def validate_compile_log_text(
    stdout: str, stderr: str, *, sanitizer: bool,
    compiler: pathlib.Path, linker: pathlib.Path,
) -> None:
    linker_version = subprocess.run(
        [str(linker), "--version"], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    ).stdout.splitlines()[0]
    if (stdout != f"{linker_version}\n" or not stderr.endswith("\n") or
            "FORGED" in stderr):
        raise ValueError("standalone compiler/linker log framing changed")
    lowered = stderr.lower()
    if "undefined reference" in lowered or "error:" in lowered:
        raise ValueError("standalone compiler log contains a failure")
    for required in (
        "c7_exact_low.cpp", "liblibleopard.a", "-std=c++11", str(linker),
    ):
        if required not in stderr:
            raise ValueError("standalone verbose compile closure is incomplete")
    if sanitizer:
        if ("clang version" not in lowered or "libclang_rt.asan" not in stderr or
                "-fsanitize=address" not in stderr):
            raise ValueError("sanitizer link-driver evidence changed")
    elif (f"COLLECT_GCC={compiler}" not in stderr or
          "collect2" not in lowered):
        raise ValueError("GCC link-driver evidence changed")


def validate_symbol_scan_text(
    text: str, *, target: str, archive: bool,
    expected_counts: dict[str, int],
    expected_members: dict[str, dict[str, int]],
) -> None:
    if text and not text.endswith("\n"):
        raise ValueError("symbol scan is not canonical line output")
    counts = {"asan_lines": 0, "ubsan_lines": 0}
    members = {member: {"asan_lines": 0, "ubsan_lines": 0}
               for member in expected_members}
    line_pattern = re.compile(
        r"(?:(?:[0-9a-f]{16}) |(?: {17}))([A-Za-z?]) (\S+)\Z")
    for line in text.splitlines():
        prefix = f"{target}:"
        if not line.startswith(prefix):
            raise ValueError("symbol scan target prefix changed")
        body = line[len(prefix):]
        member = None
        if archive:
            if ":" not in body:
                raise ValueError("archive symbol line omits its member")
            member, body = body.split(":", 1)
            if member not in members:
                raise ValueError("archive symbol line names an unknown member")
        match = line_pattern.fullmatch(body)
        if not match:
            raise ValueError("symbol scan line is not nm canonical format")
        symbol = match.group(2)
        asan = "__asan_" in symbol
        ubsan = "__ubsan_" in symbol
        if not asan and not ubsan:
            raise ValueError("symbol scan contains a non-sanitizer symbol")
        counts["asan_lines"] += asan
        counts["ubsan_lines"] += ubsan
        if member is not None:
            members[member]["asan_lines"] += asan
            members[member]["ubsan_lines"] += ubsan
    if counts != expected_counts or members != expected_members:
        raise ValueError("sanitizer symbol counts or member attribution changed")
    if expected_counts["asan_lines"]:
        for required in (
            "__asan_init", "__asan_report_load1", "__asan_report_store1",
            "__ubsan_handle_pointer_overflow",
            "__ubsan_handle_type_mismatch_v1",
        ):
            if required not in text:
                raise ValueError("sanitizer symbol family is incomplete")


def validate_run_log_text(stdout: str, stderr: str, *, smoke: bool) -> None:
    if stdout != "" or stderr != ("C7 benchmark 1/1\n" if smoke else ""):
        raise ValueError("run log semantics changed")


def validate_profile(data: dict) -> None:
    if data != {
        "family": 3, "version": 1, "coordinate_map": 1,
        "systematic": "0..K-1", "parity": "K..K+R-1",
        "production_enabled": False,
    }:
        raise ValueError("C7 profile identity changed")


def validate_algebra(data: dict) -> None:
    if set(data) != {
        "coordinate_comparison", "derivation", "gf4_exhaustive",
        "large_field", "profile", "schema", "source_sha256", "status"
    }:
        raise ValueError("unexpected algebra keys")
    if data["schema"] != "leopard2-c7-algebra/v1" or data["status"] != "pass":
        raise ValueError("bad algebra status")
    if data["profile"] != {
        "identity_family": 3,
        "profile_version": 1,
        "coordinate_map_version": 1,
        "systematic_coordinates": "omega_0 .. omega_(K-1)",
        "parity_coordinates": "omega_K .. omega_(K+R-1)",
        "parent_count": "K+R", "padded_side": "K",
        "shortening": "none", "puncturing": "none",
        "future_exact_high_identity_family": 4,
    }:
        raise ValueError("algebra profile changed")
    source_hash = digest(HERE / "algebra.py")
    if data["source_sha256"] != source_hash:
        raise ValueError("algebra source hash mismatch")
    small = data["gf4_exhaustive"]
    for key, value in {
        "geometries": 120, "mds_subsets": 131038,
        "dense_coefficients": 3060, "vandermonde_oracle_coefficients": 3060,
        "affine_maps": 28800, "affine_coefficients": 734400,
        "searched_full_field_partitions": 65534,
        "transform_search_candidates": 65519, "transform_search_wins": 57,
    }.items():
        if small.get(key) != value:
            raise ValueError(f"algebra counter changed: {key}")
    if small.get("aligned_union") != {
        "available_geometries": 85, "globally_affine_geometries": 49,
        "nonaffine_changed_coefficients": 758, "dense_coefficients": 1812,
        "prefix_one_coefficients": 297, "aligned_one_coefficients": 198,
        "prefix_symbolic_fragments": 358, "aligned_symbolic_fragments": 328,
    }:
        raise ValueError("aligned-union comparison changed")
    if len(small.get("search", [])) != 15 or len(
            small.get("transform_search_focus", [])) != 5:
        raise ValueError("coordinate searches are incomplete")
    large = data["large_field"]
    if (large.get("cases"), large.get("dense_coefficients"),
            large.get("affine_coefficients")) != (10, 84781, 177888):
        raise ValueError("large-field counters changed")
    records = large.get("records")
    aligned = large.get("aligned_union_records")
    if not isinstance(records, list) or len(records) != 10 or not isinstance(
            aligned, list) or len(aligned) != 10:
        raise ValueError("large-field records are incomplete")
    declared_gf16 = records[6]
    if declared_gf16 != {
        "field": "legacy-gf16-v1", "K": 3, "R": 500,
        "dense_coefficients": 1500, "vandermonde_rows": 500,
        "coordinate_digest_sha256":
            "78f1b8baa0d6a413c44690b7a5ca9bab90d623eb2b6d19c615690e59f665ccc9",
    }:
        raise ValueError("declared GF16 Vandermonde oracle is missing")
    if aligned[0].get("available") is not False:
        raise ValueError("aligned map crosses GF8 boundary")
    gf16_nonaffine = aligned[7]
    if (gf16_nonaffine.get("globally_affine") is not False or
            gf16_nonaffine.get("changed_coefficients") != 12900):
        raise ValueError("non-affine aligned GF16 witness changed")
    if data["coordinate_comparison"].get("decision") != (
            "freeze the simple prefix map; do not add search tables to V1"):
        raise ValueError("coordinate-map decision changed")


def validate_cpp_result(
    data: dict, backend: str, timing_scope: str,
    expected_affinity: list[int], sanitizer: bool,
) -> None:
    required = {
        "affinity", "allocation_tracking", "benchmarks", "core_git_sha",
        "correctness", "library_sha256", "omp_dynamic", "omp_num_threads",
        "production_constructor_rejected", "profile", "requested_backend",
        "runtime_backend", "sanitizer", "sanitizer_features", "schema",
        "source_sha256", "status", "timing_scope",
    }
    if set(data) != required:
        raise ValueError("unexpected C++ result keys")
    if data["schema"] != "leopard2-c7-exact-low/v1" or data["status"] != "pass":
        raise ValueError("bad C++ status")
    validate_profile(data["profile"])
    if data["production_constructor_rejected"] is not True:
        raise ValueError("production unexpectedly enabled C7")
    if data["timing_scope"] != timing_scope:
        raise ValueError("timing scope changed")
    if data["affinity"] != expected_affinity:
        raise ValueError("child affinity changed")
    if data["omp_num_threads"] != "1" or data["omp_dynamic"] != "FALSE":
        raise ValueError("thread environment changed")
    validate_sha(data["source_sha256"], "C++ source hash")
    validate_sha(data["library_sha256"], "C++ library hash")
    validate_git_sha(data["core_git_sha"], "C++ core commit")
    if data["correctness"] != EXPECTED_CORRECTNESS:
        raise ValueError("C++ correctness counters changed")
    if not isinstance(data["benchmarks"], list):
        raise ValueError("C++ benchmark collection is not an array")
    if sanitizer:
        if (data["requested_backend"] != "auto" or
                data["runtime_backend"] != EXPECTED_RUNTIME["auto"] or
                data["sanitizer"] != "asan-ubsan" or
                data["allocation_tracking"] != "disabled-for-sanitizer" or
                data["sanitizer_features"] != {
                    "address": True, "undefined": True}):
            raise ValueError("sanitizer feature proof is absent")
    else:
        if (data["requested_backend"] != backend or
                data["runtime_backend"] != EXPECTED_RUNTIME[backend] or
                data["sanitizer"] != "none" or
                data["allocation_tracking"] != "global-new" or
                data["sanitizer_features"] != {
                    "address": False, "undefined": False}):
            raise ValueError("normal backend provenance changed")
    if timing_scope == "none-correctness-only" and data["benchmarks"] != []:
        raise ValueError("correctness artifact contains timing")


def validate_smoke(data: dict, expected_affinity: list[int]) -> None:
    validate_cpp_result(
        data, "auto", "non-authoritative-smoke", expected_affinity, False)
    cells = data["benchmarks"]
    if not isinstance(cells, list) or len(cells) != 1:
        raise ValueError("smoke must contain exactly one cell")
    cell = cells[0]
    expected_keys = {
        "K", "R", "batch", "bytes", "exact_coefficients", "exact_decode",
        "exact_decode_samples_us", "exact_decode_setup",
        "exact_decode_setup_samples_us", "exact_encode",
        "exact_encode_samples_us", "exact_field", "exact_setup",
        "exact_setup_samples_us", "losses", "padded_decode",
        "padded_decode_samples_us", "padded_decode_scratch",
        "padded_decode_setup", "padded_decode_setup_samples_us",
        "padded_encode", "padded_encode_samples_us", "padded_encode_scratch",
        "padded_field", "padded_setup", "padded_setup_samples_us",
        "exact_decode_terms",
    }
    if set(cell) != expected_keys:
        raise ValueError("smoke cell schema changed")
    if [cell[key] for key in ("K", "R", "bytes", "batch", "losses")] != [
            3, 253, 64, 8, 3]:
        raise ValueError("smoke geometry changed")
    if cell["exact_field"] != 1 or cell["padded_field"] != 2:
        raise ValueError("smoke field selection changed")
    for key in ("padded_encode_scratch", "padded_decode_scratch"):
        if type(cell[key]) is not int or cell[key] < 0:
            raise ValueError("smoke scratch accounting is invalid")
    sample_keys = (
        "exact_setup_samples_us", "padded_setup_samples_us",
        "exact_decode_setup_samples_us", "padded_decode_setup_samples_us",
        "exact_encode_samples_us", "padded_encode_samples_us",
        "exact_decode_samples_us", "padded_decode_samples_us",
    )
    for key in sample_keys:
        values = cell[key]
        if not isinstance(values, list) or len(values) != 7 or not all(
                type(value) in (int, float) and math.isfinite(value) and value > 0
                for value in values):
            raise ValueError(f"smoke raw samples are incomplete: {key}")
    for summary_key, sample_key in (
        ("exact_setup", "exact_setup_samples_us"),
        ("padded_setup", "padded_setup_samples_us"),
        ("exact_decode_setup", "exact_decode_setup_samples_us"),
        ("padded_decode_setup", "padded_decode_setup_samples_us"),
        ("exact_encode", "exact_encode_samples_us"),
        ("padded_encode", "padded_encode_samples_us"),
        ("exact_decode", "exact_decode_samples_us"),
        ("padded_decode", "padded_decode_samples_us"),
    ):
        summary = cell[summary_key]
        if not isinstance(summary, dict) or set(summary) != {"median_us", "mad_us"}:
            raise ValueError("smoke summary schema changed")
        samples = cell[sample_key]
        median = statistics.median(samples)
        mad = statistics.median(abs(value - median) for value in samples)
        if summary != {"median_us": median, "mad_us": mad}:
            raise ValueError("smoke summary does not match raw samples")
    if cell["exact_coefficients"] != 759 or cell["exact_decode_terms"] != 9:
        raise ValueError("smoke exact accounting changed")


def validate_manifest(data: dict) -> None:
    if set(data) != {
        "builds", "core_git_sha", "runner", "runs", "schema", "scope",
        "source", "status", "taskset",
    }:
        raise ValueError("manifest keys changed")
    if (data["schema"] != "leopard2-c7-build-run-manifest/v2" or
            data["status"] != "pass" or data["scope"] !=
            "correctness plus CPU0 non-authoritative harness smoke; no promotion timing"):
        raise ValueError("manifest scope/status changed")
    core_sha = validate_git_sha(data["core_git_sha"], "manifest core commit")
    source_path = validate_artifact(data["source"], "C7 source", require_retained=True)
    runner_path = validate_artifact(data["runner"], "C7 runner", require_retained=True)
    if source_path != (HERE / "c7_exact_low.cpp").resolve() or runner_path != (
            HERE / "run_matrix.py").resolve():
        raise ValueError("manifest source paths changed")
    validate_git_artifact(core_sha, data["source"], "C7 source")
    validate_git_artifact(core_sha, data["runner"], "C7 runner")
    taskset = validate_program(
        data["taskset"], "taskset",
        expected_program(("taskset",), resolve=True))
    builds = data["builds"]
    if not isinstance(builds, list) or [item.get("name") for item in builds] != list(
            BUILD_NAMES):
        raise ValueError("manifest build matrix changed")
    by_name: dict[str, dict] = {}
    for build in builds:
        required = {
            "ar", "backend", "build_argv", "build_dir", "build_stderr",
            "build_stdout", "c_compiler", "cmake", "cmake_cache",
            "cmake_linker", "compile_argv", "compile_stderr",
            "compile_stdout", "compiler", "configure_argv",
            "configure_stderr", "configure_stdout", "cxx_compiler",
            "executable", "gmake", "instrumentation", "library",
            "link_driver", "name", "nm", "ranlib", "sanitizer",
            "source_closure", "standalone_linker",
        }
        if set(build) != required:
            raise ValueError("build record schema changed")
        name = build["name"]
        sanitizer = name == "asan-ubsan"
        if build["sanitizer"] is not sanitizer or build["backend"] != (
                "auto" if sanitizer else name):
            raise ValueError("build backend/sanitizer label changed")
        for label in (
            "configure_stdout", "configure_stderr", "build_stdout",
            "build_stderr", "compile_stdout", "compile_stderr", "cmake_cache",
        ):
            validate_artifact(build[label], f"{name} {label}", require_retained=True)
        library = validate_artifact(
            build["library"], f"{name} archive", require_retained=False)
        executable = validate_artifact(
            build["executable"], f"{name} executable", require_retained=False)
        expected_cmake = expected_program(("cmake",), resolve=True)
        expected_nm = expected_program(("nm",), resolve=True)
        expected_gmake = expected_program(("gmake", "make"), resolve=False)
        if sanitizer:
            expected_c = expected_program(
                ("clang", "clang-18", "clang-17", "clang-16"), resolve=False)
            expected_cxx = expected_program(
                ("clang++", "clang++-18", "clang++-17", "clang++-16"),
                resolve=False)
            expected_ar = expected_program(
                ("llvm-ar-18", "llvm-ar-17", "llvm-ar-16", "llvm-ar"),
                resolve=False)
            expected_ranlib = expected_program(
                ("llvm-ranlib-18", "llvm-ranlib-17", "llvm-ranlib-16",
                 "llvm-ranlib"), resolve=False)
            expected_cmake_linker = expected_program(("ld",), resolve=False)
        else:
            expected_c = expected_program(("gcc",), resolve=True)
            expected_cxx = expected_program(("g++",), resolve=True)
            expected_ar = expected_program(
                ("x86_64-linux-gnu-ar", "ar"), resolve=False)
            expected_ranlib = expected_program(
                ("x86_64-linux-gnu-ranlib", "ranlib"), resolve=False)
            expected_cmake_linker = expected_program(
                ("x86_64-linux-gnu-ld", "ld"), resolve=False)
        expected_standalone_linker = compiler_program(expected_cxx, "ld")
        cmake = validate_program(build["cmake"], f"{name} cmake", expected_cmake)
        c_compiler = validate_program(
            build["c_compiler"], f"{name} C compiler", expected_c)
        cxx_compiler = validate_program(
            build["cxx_compiler"], f"{name} C++ compiler", expected_cxx)
        compiler_path = validate_program(
            build["compiler"], f"{name} standalone compiler", expected_cxx)
        link_driver = validate_program(
            build["link_driver"], f"{name} link driver", expected_cxx)
        standalone_linker = validate_program(
            build["standalone_linker"], f"{name} standalone linker",
            expected_standalone_linker)
        nm = validate_program(build["nm"], f"{name} nm", expected_nm)
        gmake = validate_program(
            build["gmake"], f"{name} make", expected_gmake)
        ar = validate_program(build["ar"], f"{name} ar", expected_ar)
        ranlib = validate_program(
            build["ranlib"], f"{name} ranlib", expected_ranlib)
        cmake_linker = validate_program(
            build["cmake_linker"], f"{name} CMake linker",
            expected_cmake_linker)
        if (build["compiler"] != build["cxx_compiler"] or
                build["link_driver"] != build["cxx_compiler"] or
                link_driver != compiler_path):
            raise ValueError("standalone/core/link-driver roles differ")
        configure = build["configure_argv"]
        build_argv = build["build_argv"]
        compile_argv = build["compile_argv"]
        if not all(isinstance(argv, list) and all(isinstance(item, str) for item in argv)
                   for argv in (configure, build_argv, compile_argv)):
            raise ValueError("build argv is not an exact string array")
        expected_build_dir = f"build/c7-matrix/core-{name}"
        expected_library = f"{expected_build_dir}/liblibleopard.a"
        expected_executable = f"build/c7-matrix/c7-{name}"
        if (build["build_dir"] != expected_build_dir or
                build["library"]["path"] != expected_library or
                build["executable"]["path"] != expected_executable or
                library != (ROOT / expected_library).resolve() or
                executable != (ROOT / expected_executable).resolve()):
            raise ValueError("build output paths changed")
        expected_logs = {
            "configure_stdout": f"experiments/leopard2/non_power_of_two/c7/results/logs/{name}-configure.stdout.txt",
            "configure_stderr": f"experiments/leopard2/non_power_of_two/c7/results/logs/{name}-configure.stderr.txt",
            "build_stdout": f"experiments/leopard2/non_power_of_two/c7/results/logs/{name}-core-build.stdout.txt",
            "build_stderr": f"experiments/leopard2/non_power_of_two/c7/results/logs/{name}-core-build.stderr.txt",
            "compile_stdout": f"experiments/leopard2/non_power_of_two/c7/results/logs/{name}-compile.stdout.txt",
            "compile_stderr": f"experiments/leopard2/non_power_of_two/c7/results/logs/{name}-compile.stderr.txt",
            "cmake_cache": f"experiments/leopard2/non_power_of_two/c7/results/logs/{name}-CMakeCache.txt",
        }
        if any(build[label]["path"] != path for label, path in expected_logs.items()):
            raise ValueError("build log path changed")
        expected_configure = [
            str(cmake), "-S", ".", "-B", expected_build_dir,
            "-G", "Unix Makefiles",
            f"-DCMAKE_BUILD_TYPE={'Debug' if sanitizer else 'Release'}",
            f"-DCMAKE_C_COMPILER={c_compiler}",
            f"-DCMAKE_CXX_COMPILER={cxx_compiler}",
            f"-DLEO2_BACKEND_VARIANT={build['backend']}",
            "-DLEO2_BUILD_TESTS=OFF", "-DLEO2_BUILD_BENCHMARKS=OFF",
            "-DLEO2_BUILD_FUZZERS=OFF", "-DLEO2_ENABLE_CUDA=OFF",
            f"-DENABLE_OPENMP={'OFF' if sanitizer else 'ON'}",
        ]
        if sanitizer:
            core_flags = "-fsanitize=address,undefined -fno-omit-frame-pointer"
            expected_configure += [
                f"-DCMAKE_C_FLAGS={core_flags}",
                f"-DCMAKE_CXX_FLAGS={core_flags}",
                f"-DCMAKE_EXE_LINKER_FLAGS={core_flags}",
            ]
        expected_build = [
            str(cmake), "--build", expected_build_dir, "--target",
            "libleopard", "--verbose", "--", "-j4",
        ]
        expected_compile = [
            str(compiler_path), "-v", "-Wl,-v", "-std=c++11", "-g", "-Wall",
            "-Wextra",
            "-Wpedantic", "-Werror", "-I.",
            f'-DLEO2_C7_SOURCE_SHA256="{data["source"]["sha256"]}"',
            f'-DLEO2_C7_CORE_GIT_SHA="{core_sha}"',
            f'-DLEO2_C7_LIBRARY_SHA256="{build["library"]["sha256"]}"',
        ]
        if sanitizer:
            expected_compile += [
                "-O1", "-fsanitize=address,undefined",
                "-fno-omit-frame-pointer",
                "-DLEO2_C7_DISABLE_GLOBAL_NEW_TRACKING=1",
                '-DLEO2_C7_SANITIZER_MODE="asan-ubsan"',
                "-DLEO2_C7_REQUIRE_ASAN_UBSAN=1",
            ]
        else:
            expected_compile += ["-O2"]
        expected_compile += [
            data["source"]["path"], expected_library, "-pthread",
        ]
        if not sanitizer:
            expected_compile += ["-fopenmp"]
        expected_compile += ["-o", expected_executable]
        if configure != expected_configure or build_argv != expected_build or (
                compile_argv != expected_compile):
            raise ValueError("exact configure/build/compile command changed")
        cache_text = resolve_record_path(
            build["cmake_cache"]["path"]).read_text(encoding="utf-8")
        validate_cmake_cache_text(
            cache_text, sanitizer=sanitizer, backend=build["backend"],
            build_dir=expected_build_dir, cmake=cmake,
            c_compiler=c_compiler, cxx_compiler=cxx_compiler, gmake=gmake,
            ar=ar, ranlib=ranlib, cmake_linker=cmake_linker)
        configure_stdout_text = resolve_record_path(
            build["configure_stdout"]["path"]).read_text(encoding="utf-8")
        configure_stderr_text = resolve_record_path(
            build["configure_stderr"]["path"]).read_text(encoding="utf-8")
        build_stdout_text = resolve_record_path(
            build["build_stdout"]["path"]).read_text(encoding="utf-8")
        build_stderr_text = resolve_record_path(
            build["build_stderr"]["path"]).read_text(encoding="utf-8")
        compile_stdout_text = resolve_record_path(
            build["compile_stdout"]["path"]).read_text(encoding="utf-8")
        compile_stderr_text = resolve_record_path(
            build["compile_stderr"]["path"]).read_text(encoding="utf-8")
        if configure_stderr_text != "" or build_stderr_text != "":
            raise ValueError("configure/core-build stderr is not empty")
        validate_configure_log_text(
            configure_stdout_text, sanitizer=sanitizer,
            c_compiler=c_compiler, cxx_compiler=cxx_compiler,
            build_dir=expected_build_dir)
        validate_build_log_text(
            build_stdout_text, sanitizer=sanitizer, cmake=cmake,
            cxx_compiler=cxx_compiler, gmake=gmake, ar=ar, ranlib=ranlib)
        validate_compile_log_text(
            compile_stdout_text, compile_stderr_text, sanitizer=sanitizer,
            compiler=compiler_path, linker=standalone_linker)
        instrumentation = build["instrumentation"]
        if set(instrumentation) != {
            "archive_members", "core_archive_counts",
            "core_archive_member_counts", "core_archive_symbol_scan",
            "executable_counts", "executable_symbol_scan",
            "required_compile_macro",
        }:
            raise ValueError("instrumentation schema changed")
        executable_scan = validate_artifact(
            instrumentation["executable_symbol_scan"],
            f"{name} executable nm scan", require_retained=True)
        archive_scan = validate_artifact(
            instrumentation["core_archive_symbol_scan"],
            f"{name} core archive nm scan", require_retained=True)
        scan_prefix = "experiments/leopard2/non_power_of_two/c7/results/logs"
        if (instrumentation["executable_symbol_scan"]["path"] !=
                f"{scan_prefix}/{name}-nm-executable-sanitizers.txt" or
                instrumentation["core_archive_symbol_scan"]["path"] !=
                f"{scan_prefix}/{name}-nm-core-archive-sanitizers.txt"):
            raise ValueError("nm scan path changed")
        expected_executable_counts = (
            {"asan_lines": 320, "ubsan_lines": 54} if sanitizer else
            {"asan_lines": 0, "ubsan_lines": 0})
        expected_archive_counts = (
            {"asan_lines": 306, "ubsan_lines": 86} if sanitizer else
            {"asan_lines": 0, "ubsan_lines": 0})
        expected_member_counts = {
            member: {
                "asan_lines": SANITIZED_MEMBER_COUNTS[member][0]
                if sanitizer else 0,
                "ubsan_lines": SANITIZED_MEMBER_COUNTS[member][1]
                if sanitizer else 0,
            }
            for member in ARCHIVE_MEMBERS
        }
        if (instrumentation["required_compile_macro"] is not sanitizer or
                instrumentation["archive_members"] != list(ARCHIVE_MEMBERS) or
                instrumentation["executable_counts"] !=
                expected_executable_counts or
                instrumentation["core_archive_counts"] !=
                expected_archive_counts or
                instrumentation["core_archive_member_counts"] !=
                expected_member_counts):
            raise ValueError("sanitizer instrumentation summary changed")
        validate_symbol_scan_text(
            executable_scan.read_text(encoding="utf-8"),
            target=expected_executable, archive=False,
            expected_counts=expected_executable_counts, expected_members={})
        validate_symbol_scan_text(
            archive_scan.read_text(encoding="utf-8"),
            target=expected_library, archive=True,
            expected_counts=expected_archive_counts,
            expected_members=expected_member_counts)
        closure = build["source_closure"]
        if not isinstance(closure, list) or len(closure) < 10:
            raise ValueError("core source closure is incomplete")
        closure_paths = set()
        closure_order = []
        for index, entry in enumerate(closure):
            path = validate_artifact(
                entry, f"{name} source closure {index}", require_retained=True)
            try:
                relative_path = path.relative_to(ROOT).as_posix()
            except ValueError as error:
                raise ValueError("source closure escaped repository") from error
            closure_paths.add(relative_path)
            closure_order.append(relative_path)
            git_bytes = subprocess.run(
                ["git", "show", f"{core_sha}:{relative_path}"], cwd=ROOT,
                check=False, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            )
            if git_bytes.returncode != 0 or hashlib.sha256(
                    git_bytes.stdout).hexdigest() != entry["sha256"]:
                raise ValueError("source closure is not bound to core commit")
        if len(closure_paths) != len(closure_order) or closure_order != sorted(
                closure_order):
            raise ValueError("source closure is duplicate or noncanonical")
        if not {"CMakeLists.txt", "leopard2.cpp", "LeopardFF8.cpp",
                "LeopardFF16.cpp"}.issubset(closure_paths):
            raise ValueError("core source closure lacks required units")
        by_name[name] = build

    runs = data["runs"]
    expected_runs = [*BUILD_NAMES, "smoke-nonauthoritative"]
    if not isinstance(runs, list) or [run.get("name") for run in runs] != expected_runs:
        raise ValueError("manifest run matrix changed")
    seen_cpus = []
    for run in runs:
        if set(run) != {
            "argv", "build", "environment", "kind", "name", "observed_affinity",
            "requested_cpu", "result", "stderr", "stdout",
        }:
            raise ValueError("run record schema changed")
        name = run["name"]
        smoke = name == "smoke-nonauthoritative"
        build_name = "auto" if smoke else name
        if run["build"] != build_name or run["kind"] != (
                "non-authoritative-smoke" if smoke else "correctness"):
            raise ValueError("run kind/build changed")
        cpu = run["requested_cpu"]
        if type(cpu) is not int or cpu in (15, 31) or run["observed_affinity"] != [cpu]:
            raise ValueError("run affinity is invalid or used a reserved CPU")
        if not smoke:
            seen_cpus.append(cpu)
        expected_environment = {
            "LC_ALL": "C", "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
        }
        if build_name == "asan-ubsan":
            expected_environment.update({
                "ASAN_OPTIONS": "detect_leaks=1:halt_on_error=1",
                "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1",
            })
        if run["environment"] != expected_environment:
            raise ValueError("run environment changed")
        result_path = validate_artifact(
            run["result"], f"{name} result", require_retained=True)
        stdout_path = validate_artifact(
            run["stdout"], f"{name} stdout", require_retained=True)
        stderr_path = validate_artifact(
            run["stderr"], f"{name} stderr", require_retained=True)
        result_prefix = "experiments/leopard2/non_power_of_two/c7/results"
        if (run["result"]["path"] != f"{result_prefix}/{name}.json" or
                run["stdout"]["path"] !=
                f"{result_prefix}/logs/{name}-run.stdout.txt" or
                run["stderr"]["path"] !=
                f"{result_prefix}/logs/{name}-run.stderr.txt"):
            raise ValueError("run result/log path changed")
        argv = run["argv"]
        expected_argv = [
            str(taskset), "-c", str(cpu), by_name[build_name]["executable"]["path"],
            "--backend", "auto" if build_name == "asan-ubsan" else
            by_name[build_name]["backend"], run["result"]["path"],
            "--benchmark-smoke" if smoke else "--correctness-only",
        ]
        if argv != expected_argv:
            raise ValueError("run argv changed")
        validate_run_log_text(
            stdout_path.read_text(encoding="utf-8"),
            stderr_path.read_text(encoding="utf-8"), smoke=smoke)
        build = by_name[build_name]
        child = json.loads(result_path.read_text(encoding="utf-8"))
        if (child["source_sha256"] != data["source"]["sha256"] or
                child["core_git_sha"] != core_sha or
                child["library_sha256"] != build["library"]["sha256"]):
            raise ValueError("run result is not bound to build inputs")
        if smoke:
            if cpu != 0:
                raise ValueError("retained non-authoritative smoke is not CPU0")
            validate_smoke(child, [cpu])
        else:
            validate_cpp_result(
                child, "auto" if build_name == "asan-ubsan" else build_name,
                "none-correctness-only", [cpu], build_name == "asan-ubsan")
    if len(seen_cpus) != 5 or len(set(seen_cpus)) != 5:
        raise ValueError("correctness runs were not independently pinned")


class CheckpointTests(unittest.TestCase):
    def test_algebra_canonical_regeneration(self) -> None:
        retained = RESULTS / "algebra.json"
        validate_algebra(json.loads(retained.read_text(encoding="utf-8")))
        with tempfile.TemporaryDirectory(prefix="leopard2-c7-algebra-") as directory:
            regenerated = pathlib.Path(directory) / "algebra.json"
            subprocess.run([
                sys.executable, "-X", "dev", str(HERE / "algebra.py"),
                "--output", str(regenerated),
            ], cwd=ROOT, check=True, env={
                **dict(os.environ),
                "PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0",
            })
            self.assertEqual(regenerated.read_bytes(), retained.read_bytes())

    def test_manifest_and_full_closure(self) -> None:
        validate_manifest(load("build-run-manifest.json"))

    def test_algebra_mutations_rejected(self) -> None:
        original = load("algebra.json")
        for path, value in (
            (("profile", "identity_family"), 4),
            (("gf4_exhaustive", "mds_subsets"), 1),
            (("gf4_exhaustive", "aligned_union", "globally_affine_geometries"), 85),
            (("large_field", "dense_coefficients"), 0),
            (("large_field", "records", 6, "vandermonde_rows"), 8),
            (("coordinate_comparison", "decision"), "aligned"),
            (("source_sha256",), "0" * 64),
        ):
            candidate = copy.deepcopy(original)
            cursor = candidate
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = value
            with self.assertRaises(ValueError):
                validate_algebra(candidate)

    def test_manifest_mutations_rejected(self) -> None:
        original = load("build-run-manifest.json")
        mutations: list[dict] = []
        for transform in (
            lambda item: item["source"].update(sha256="g" * 64),
            lambda item: item["builds"][0]["compile_argv"].remove("-Werror"),
            lambda item: item["builds"][4]["instrumentation"].update(
                executable_counts={"asan_lines": 320, "ubsan_lines": 53}),
            lambda item: item["runs"][0].update(requested_cpu=15),
            lambda item: item["runs"][5].update(kind="authoritative"),
            lambda item: item["runs"][5]["result"].update(sha256="0" * 64),
        ):
            candidate = copy.deepcopy(original)
            transform(candidate)
            mutations.append(candidate)
        for candidate in mutations:
            with self.assertRaises((ValueError, OSError)):
                validate_manifest(candidate)

    def test_coordinated_program_forgeries_rejected(self) -> None:
        original = load("build-run-manifest.json")

        taskset_forgery = copy.deepcopy(original)
        forged_taskset = copy.deepcopy(
            taskset_forgery["builds"][0]["cxx_compiler"])
        taskset_forgery["taskset"] = forged_taskset
        for run in taskset_forgery["runs"]:
            run["argv"][0] = forged_taskset["path"]
        with self.assertRaises(ValueError):
            validate_manifest(taskset_forgery)

        compiler_forgery = copy.deepcopy(original)
        build = compiler_forgery["builds"][0]
        forged_cxx = copy.deepcopy(
            compiler_forgery["builds"][4]["cxx_compiler"])
        for role in ("cxx_compiler", "compiler", "link_driver"):
            build[role] = copy.deepcopy(forged_cxx)
        for index, argument in enumerate(build["configure_argv"]):
            if argument.startswith("-DCMAKE_CXX_COMPILER="):
                build["configure_argv"][index] = (
                    f"-DCMAKE_CXX_COMPILER={forged_cxx['path']}")
        build["compile_argv"][0] = forged_cxx["path"]
        with self.assertRaises(ValueError):
            validate_manifest(compiler_forgery)

    def test_rehashed_semantic_log_forgeries_rejected(self) -> None:
        original = load("build-run-manifest.json")

        def require_rejection(
            record: dict, candidate: dict, forged: bytes,
        ) -> None:
            path = resolve_record_path(record["path"])
            retained = path.read_bytes()
            try:
                path.write_bytes(forged)
                record.update(
                    bytes=len(forged),
                    sha256=hashlib.sha256(forged).hexdigest())
                with self.assertRaises(ValueError):
                    validate_manifest(candidate)
            finally:
                path.write_bytes(retained)

        candidate = copy.deepcopy(original)
        record = candidate["builds"][0]["configure_stdout"]
        require_rejection(
            record, candidate, b"FORGED CONFIGURE SUCCESS\n")

        candidate = copy.deepcopy(original)
        record = candidate["builds"][0]["compile_stderr"]
        require_rejection(record, candidate, b"FORGED COMPILE SUCCESS\n")

        candidate = copy.deepcopy(original)
        record = candidate["runs"][0]["stdout"]
        require_rejection(record, candidate, b"FORGED RUN SUCCESS\n")

        candidate = copy.deepcopy(original)
        instrumentation = candidate["builds"][4]["instrumentation"]
        record = instrumentation["executable_symbol_scan"]
        target = candidate["builds"][4]["executable"]["path"]
        asan_line = f"{target}:{' ' * 17}U __asan_init\n"
        ubsan_line = (
            f"{target}:{' ' * 17}U __ubsan_handle_pointer_overflow\n")
        forged_scan = (asan_line * 320 + ubsan_line * 54).encode("utf-8")
        require_rejection(record, candidate, forged_scan)


if __name__ == "__main__":
    unittest.main(verbosity=2)
