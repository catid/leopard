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
import subprocess
import sys
import tempfile
import unittest
import statistics


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


def validate_program(record: object, label: str) -> pathlib.Path:
    if not isinstance(record, dict) or set(record) != {
            "path", "sha256", "version"}:
        raise ValueError(f"{label} program schema changed")
    path_text = record["path"]
    if not isinstance(path_text, str) or not pathlib.Path(path_text).is_absolute():
        raise ValueError(f"{label} path is not absolute")
    validate_sha(record["sha256"], f"{label} hash")
    if not isinstance(record["version"], str) or not record["version"].strip():
        raise ValueError(f"{label} version is empty")
    return pathlib.Path(path_text)


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
    if (data["schema"] != "leopard2-c7-build-run-manifest/v1" or
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
    taskset = validate_program(data["taskset"], "taskset")
    if taskset.name != "taskset":
        raise ValueError("taskset program changed")
    builds = data["builds"]
    if not isinstance(builds, list) or [item.get("name") for item in builds] != list(
            BUILD_NAMES):
        raise ValueError("manifest build matrix changed")
    by_name: dict[str, dict] = {}
    for build in builds:
        required = {
            "backend", "build_argv", "build_dir", "build_stderr",
            "build_stdout", "c_compiler", "cmake", "compile_argv",
            "compile_stderr", "compile_stdout", "compiler", "configure_argv",
            "configure_stderr", "configure_stdout", "cxx_compiler",
            "executable", "instrumentation", "library", "name", "nm",
            "sanitizer", "source_closure",
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
            "build_stderr", "compile_stdout", "compile_stderr",
        ):
            validate_artifact(build[label], f"{name} {label}", require_retained=True)
        library = validate_artifact(
            build["library"], f"{name} archive", require_retained=False)
        executable = validate_artifact(
            build["executable"], f"{name} executable", require_retained=False)
        cmake = validate_program(build["cmake"], f"{name} cmake")
        c_compiler = validate_program(build["c_compiler"], f"{name} C compiler")
        cxx_compiler = validate_program(
            build["cxx_compiler"], f"{name} C++ compiler")
        compiler_path = validate_program(
            build["compiler"], f"{name} standalone compiler")
        nm = validate_program(build["nm"], f"{name} nm")
        if build["compiler"] != build["cxx_compiler"]:
            raise ValueError("standalone and core C++ compilers differ")
        if cmake.name != "cmake" or nm.name not in {"nm", "x86_64-linux-gnu-nm"}:
            raise ValueError("build program identity changed")
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
            "libleopard", "--", "-j4",
        ]
        expected_compile = [
            str(compiler_path), "-std=c++11", "-g", "-Wall", "-Wextra",
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
        instrumentation = build["instrumentation"]
        if set(instrumentation) != {
            "has_asan_symbols", "has_ubsan_symbols", "required_compile_macro",
            "undefined_symbol_scan",
        }:
            raise ValueError("instrumentation schema changed")
        scan_path = validate_artifact(
            instrumentation["undefined_symbol_scan"], f"{name} nm scan",
            require_retained=True)
        if instrumentation["undefined_symbol_scan"]["path"] != (
                f"experiments/leopard2/non_power_of_two/c7/results/logs/"
                f"{name}-nm-undefined.txt"):
            raise ValueError("nm scan path changed")
        scan = scan_path.read_text(encoding="utf-8")
        observed_asan = "__asan_" in scan
        observed_ubsan = "__ubsan_" in scan
        if (instrumentation["required_compile_macro"] is not sanitizer or
                instrumentation["has_asan_symbols"] is not observed_asan or
                instrumentation["has_ubsan_symbols"] is not observed_ubsan or
                observed_asan is not sanitizer or observed_ubsan is not sanitizer):
            raise ValueError("sanitizer instrumentation proof changed")
        if sanitizer and "clang" not in build["compiler"]["version"].lower():
            raise ValueError("sanitizer compiler lacks __has_feature contract")
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
        validate_artifact(run["stdout"], f"{name} stdout", require_retained=True)
        validate_artifact(run["stderr"], f"{name} stderr", require_retained=True)
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
                has_ubsan_symbols=False),
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
