#!/usr/bin/env python3
"""Portable and optional live validation for a C7 matrix v3 manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re
import shutil
import subprocess
from typing import Any

from run_matrix import (
    BACKENDS,
    BUILD_NAMES,
    EXPECTED_ARCHIVE_MEMBER_COUNTS,
    EXPECTED_ARCHIVE_SANITIZER_COUNTS,
    EXPECTED_EXECUTABLE_SANITIZER_COUNTS,
    NORMALIZATION_SCHEMA,
    NORMALIZATION_TOKEN,
    PREFIX_MAP_TARGET,
)


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
GIT_SHA_RE = re.compile(r"[0-9a-f]{40}\Z")
CORE_SOURCES = {
    "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp", "LeopardCommon.cpp", "LeopardFF16.cpp",
    "LeopardFF8.cpp", "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
}
ARCHIVE_MEMBERS = tuple(EXPECTED_ARCHIVE_MEMBER_COUNTS)
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
    "decode_read_only_input_alias_calls": 117,
    "decode_read_only_input_alias_symbol_comparisons": 6025,
    "hot_path_allocations": 0,
    "digest_fnv64": "0xec4179e9f2776a58",
}
ABSOLUTE_PROJECT_PATH = re.compile(
    r"(?:[A-Za-z]:[\\/]|/)[^\s\"']*(?:[\\/])(?:"
    r"CMakeLists\.txt|cmake[\\/]leopardConfig\.cmake\.in|"
    r"experiments[\\/]leopard2[\\/]|leopard2?\.cpp|"
    r"Leopard(?:2|Common|FF)[A-Za-z0-9_]*\.(?:cpp|h))")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
        raise ValueError(f"{label} is not canonical SHA-256")
    return value


def validate_git_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or not GIT_SHA_RE.fullmatch(value):
        raise ValueError(f"{label} is not a canonical git SHA")
    return value


def resolve_path(path_text: str, source_root: pathlib.Path) -> pathlib.Path:
    path = pathlib.Path(path_text)
    return path if path.is_absolute() else source_root / path


def validate_artifact(
    record: object, label: str, source_root: pathlib.Path, *, required: bool,
    check_if_present: bool = True,
) -> pathlib.Path:
    if not isinstance(record, dict) or set(record) != {
            "bytes", "path", "sha256"}:
        raise ValueError(f"{label} artifact schema changed")
    path_text = record["path"]
    if not isinstance(path_text, str) or not path_text:
        raise ValueError(f"{label} path is invalid")
    if type(record["bytes"]) is not int or record["bytes"] < 0:
        raise ValueError(f"{label} byte count is invalid")
    expected = validate_sha(record["sha256"], f"{label} hash")
    path = resolve_path(path_text, source_root)
    if required and not path.is_file():
        raise ValueError(f"{label} retained artifact is missing")
    if path.is_file() and check_if_present:
        if path.stat().st_size != record["bytes"] or sha256(path) != expected:
            raise ValueError(f"{label} bytes disagree with the manifest")
    return path


def validate_normalized_text(
    record: object, label: str, source_root: pathlib.Path, *,
    require_token: bool,
) -> pathlib.Path:
    if not isinstance(record, dict) or set(record) != {
            "bytes", "path", "sha256", "source_root_tokens"}:
        raise ValueError(f"{label} normalized-artifact schema changed")
    token_count = record["source_root_tokens"]
    if type(token_count) is not int or token_count < 0:
        raise ValueError(f"{label} token count is invalid")
    if pathlib.Path(record["path"]).is_absolute():
        raise ValueError(f"{label} retained path is not checkout-relative")
    generic = {key: record[key] for key in ("bytes", "path", "sha256")}
    path = validate_artifact(generic, label, source_root, required=True)
    text = path.read_text(encoding="utf-8", errors="strict")
    if text.count(NORMALIZATION_TOKEN) != token_count:
        raise ValueError(f"{label} normalization token count changed")
    if require_token and token_count == 0:
        raise ValueError(f"{label} lacks its required normalization token")
    if str(source_root.resolve()) in text or ABSOLUTE_PROJECT_PATH.search(text):
        raise ValueError(f"{label} leaked an absolute checkout path")
    return path


def validate_program_record(record: object, label: str) -> pathlib.Path:
    if not isinstance(record, dict) or set(record) != {
            "path", "sha256", "version"}:
        raise ValueError(f"{label} program schema changed")
    path_text = record["path"]
    if not isinstance(path_text, str) or not pathlib.Path(path_text).is_absolute():
        raise ValueError(f"{label} path is not absolute")
    validate_sha(record["sha256"], f"{label} hash")
    if not isinstance(record["version"], str) or not record["version"].strip() or (
            not record["version"].endswith("\n")):
        raise ValueError(f"{label} version output is not canonical")
    return pathlib.Path(path_text)


def validate_program_live(record: dict, label: str) -> pathlib.Path:
    path = validate_program_record(record, label)
    if not path.is_file() or sha256(path) != record["sha256"]:
        raise ValueError(f"{label} exact executable is unavailable or changed")
    completed = subprocess.run(
        [str(path), "--version"], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if completed.stdout != record["version"]:
        raise ValueError(f"{label} exact version output changed")
    return path


def validate_argv(
    argv: object, token_count: object, label: str, *, require_token: bool,
) -> list[str]:
    if not isinstance(argv, list) or not all(isinstance(item, str) for item in argv):
        raise ValueError(f"{label} argv is not an exact string array")
    if type(token_count) is not int or token_count < 0:
        raise ValueError(f"{label} argv token count is invalid")
    text = "\n".join(argv)
    if text.count(NORMALIZATION_TOKEN) != token_count:
        raise ValueError(f"{label} argv token count changed")
    if require_token and token_count == 0:
        raise ValueError(f"{label} argv lacks its normalization token")
    if str(ROOT) in text or ABSOLUTE_PROJECT_PATH.search(text):
        raise ValueError(f"{label} argv leaked an absolute checkout path")
    return argv


def validate_git_artifact(
    commit: str, record: dict, expected_relative: str, label: str,
) -> None:
    if record["path"] != expected_relative:
        raise ValueError(f"{label} repository path changed")
    completed = subprocess.run(
        ["git", "show", f"{commit}:{expected_relative}"], cwd=ROOT,
        check=False, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if completed.returncode != 0 or hashlib.sha256(
            completed.stdout).hexdigest() != record["sha256"]:
        raise ValueError(f"{label} is not bound to the core commit")


def parse_cache(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in text.splitlines():
        match = re.fullmatch(r"([^#/:=][^:=]*):([^=]+)=(.*)", line)
        if match:
            key, value = match.group(1), match.group(3)
            if key in values:
                raise ValueError(f"duplicate CMake cache key: {key}")
            values[key] = value
    return values


def compiler_identity(record: dict, family: str) -> tuple[str, str]:
    first = record["version"].splitlines()[0]
    if family == "GNU":
        match = re.search(r"\b([0-9]+(?:\.[0-9]+){1,3})\Z", first)
    else:
        match = re.search(r"clang version ([0-9]+(?:\.[0-9]+){1,3})", first)
    if not match:
        raise ValueError("recorded compiler version is not understood")
    return family, match.group(1)


def validate_symbol_scan(
    text: str, target: str, *, archive: bool,
    expected_counts: dict[str, int],
    expected_members: dict[str, dict[str, int]],
) -> None:
    if text and not text.endswith("\n"):
        raise ValueError("symbol scan is not canonical line output")
    counts = {"asan_lines": 0, "ubsan_lines": 0}
    members = {name: {"asan_lines": 0, "ubsan_lines": 0}
               for name in expected_members}
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
                raise ValueError("archive scan omitted its member")
            member, body = body.split(":", 1)
            if member not in members:
                raise ValueError("archive scan named an unknown member")
        match = line_pattern.fullmatch(body)
        if not match:
            raise ValueError("symbol scan line is not canonical nm output")
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
        raise ValueError("sanitizer counts or archive attribution changed")
    if expected_counts["asan_lines"]:
        for symbol in (
            "__asan_init", "__asan_report_load1", "__asan_report_store1",
            "__ubsan_handle_pointer_overflow",
            "__ubsan_handle_type_mismatch_v1",
        ):
            if symbol not in text:
                raise ValueError("sanitizer symbol family is incomplete")


def live_nm(record: dict, target: pathlib.Path, retained: pathlib.Path) -> None:
    nm = validate_program_live(record, "nm")
    completed = subprocess.run(
        [str(nm), "--print-file-name", target.name], cwd=target.parent,
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if completed.stderr:
        raise ValueError("live nm emitted stderr")
    lines = [line for line in completed.stdout.splitlines()
             if b"__asan_" in line or b"__ubsan_" in line]
    actual = b"\n".join(lines) + (b"\n" if lines else b"")
    if actual != retained.read_bytes():
        raise ValueError("live nm output differs from retained bytes")


def validate_child(data: dict, backend: str, affinity: list[int], sanitizer: bool) -> None:
    if (data.get("schema") != "leopard2-c7-exact-low/v1" or
            data.get("status") != "pass" or
            data.get("production_constructor_rejected") is not True or
            data.get("correctness") != EXPECTED_CORRECTNESS or
            data.get("affinity") != affinity or data.get("benchmarks") != []):
        raise ValueError("C7 child correctness result changed")
    if data.get("requested_backend") != backend:
        raise ValueError("C7 child requested backend changed")
    if sanitizer:
        if (data.get("sanitizer") != "asan-ubsan" or
                data.get("sanitizer_features") != {
                    "address": True, "undefined": True}):
            raise ValueError("sanitizer child lacks feature proof")
    elif data.get("sanitizer") != "none":
        raise ValueError("ordinary child unexpectedly claims sanitizers")


def validate_manifest(
    data: dict[str, Any], *, source_root: pathlib.Path = ROOT,
    live: bool = False,
) -> None:
    required_top = {
        "builds", "core_git_sha", "normalization", "reproducibility",
        "runner", "runs", "schema", "scope", "source", "status",
        "taskset", "validator",
    }
    if set(data) != required_top:
        raise ValueError("manifest keys changed")
    if (data["schema"] != "leopard2-c7-build-run-manifest/v3" or
            data["status"] != "pass" or data["scope"] !=
            "correctness plus one affinity-selected non-authoritative harness "
            "smoke; no promotion timing"):
        raise ValueError("manifest status or scope changed")
    if data["normalization"] != {
        "schema": NORMALIZATION_SCHEMA,
        "token": NORMALIZATION_TOKEN,
        "operation": "replace exact source-root prefix only",
    }:
        raise ValueError("normalization identity changed")
    serialized = json.dumps(data, sort_keys=True)
    if str(source_root.resolve()) in serialized or ABSOLUTE_PROJECT_PATH.search(
            serialized):
        raise ValueError("manifest leaked an absolute checkout path")
    core_sha = validate_git_sha(data["core_git_sha"], "manifest core commit")
    source_rel = "experiments/leopard2/non_power_of_two/c7/c7_exact_low.cpp"
    runner_rel = "experiments/leopard2/non_power_of_two/c7/run_matrix.py"
    validator_rel = "experiments/leopard2/non_power_of_two/c7/validate_evidence.py"
    for key, relative, label in (
        ("source", source_rel, "C7 source"),
        ("runner", runner_rel, "C7 runner"),
        ("validator", validator_rel, "C7 validator"),
    ):
        validate_artifact(data[key], label, source_root, required=True)
        validate_git_artifact(core_sha, data[key], relative, label)
    taskset = validate_program_record(data["taskset"], "taskset")
    if live:
        validate_program_live(data["taskset"], "taskset")

    builds = data["builds"]
    if not isinstance(builds, list) or [item.get("name") for item in builds] != list(
            BUILD_NAMES):
        raise ValueError("build matrix changed")
    by_name: dict[str, dict] = {}
    build_root: str | None = None
    program_roles = (
        "ar", "c_compiler", "cmake", "cmake_linker", "compiler",
        "cxx_compiler", "gmake", "link_driver", "nm", "ranlib",
        "standalone_linker",
    )
    for build in builds:
        required = {
            "ar", "argv_source_root_tokens", "backend", "build_argv",
            "build_dir", "build_stderr", "build_stdout", "c_compiler",
            "cmake", "cmake_cache", "cmake_linker", "compile_argv",
            "compile_stderr", "compile_stdout", "compiler",
            "configure_argv", "configure_stderr", "configure_stdout",
            "cxx_compiler", "executable", "gmake", "instrumentation",
            "library", "link_driver", "name", "nm", "prefix_map_flags",
            "ranlib", "sanitizer", "source_closure", "standalone_linker",
        }
        if set(build) != required:
            raise ValueError("build record schema changed")
        name = build["name"]
        sanitizer = name == "asan-ubsan"
        if build["sanitizer"] is not sanitizer or build["backend"] != (
                "auto" if sanitizer else name):
            raise ValueError("backend/sanitizer label changed")
        records = {role: validate_program_record(build[role], f"{name} {role}")
                   for role in program_roles}
        if live:
            for role in program_roles:
                validate_program_live(build[role], f"{name} {role}")
            completed = subprocess.run(
                [str(records["compiler"]), "-print-prog-name=ld"],
                check=True, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            selected = pathlib.Path(completed.stdout.strip())
            if not selected.is_absolute():
                found = shutil.which(str(selected))
                if not found:
                    raise ValueError("live compiler-selected linker is unavailable")
                selected = pathlib.Path(found)
            if selected.resolve() != records["standalone_linker"].resolve():
                raise ValueError("live compiler-selected linker role changed")
        if (build["compiler"] != build["cxx_compiler"] or
                build["link_driver"] != build["cxx_compiler"]):
            raise ValueError("standalone compiler/link-driver role changed")
        token_counts = build["argv_source_root_tokens"]
        if not isinstance(token_counts, dict) or set(token_counts) != {
                "configure", "build", "compile"}:
            raise ValueError("build argv normalization record changed")
        configure = validate_argv(
            build["configure_argv"], token_counts["configure"],
            f"{name} configure", require_token=True)
        build_argv = validate_argv(
            build["build_argv"], token_counts["build"],
            f"{name} build", require_token=True)
        compile_argv = validate_argv(
            build["compile_argv"], token_counts["compile"],
            f"{name} compile", require_token=True)
        flags = build["prefix_map_flags"]
        mandatory = [
            f"-ffile-prefix-map={NORMALIZATION_TOKEN}={PREFIX_MAP_TARGET}",
            f"-fdebug-prefix-map={NORMALIZATION_TOKEN}={PREFIX_MAP_TARGET}",
        ]
        optional = f"-fmacro-prefix-map={NORMALIZATION_TOKEN}={PREFIX_MAP_TARGET}"
        if flags not in (mandatory, [*mandatory, optional]):
            raise ValueError("prefix-map flag set changed")
        build_dir = build["build_dir"]
        if (not isinstance(build_dir, str) or pathlib.Path(build_dir).is_absolute() or
                not build_dir.endswith(f"/core-{name}")):
            raise ValueError("build directory is not portable or canonical")
        candidate_root = build_dir.rsplit("/", 1)[0]
        if build_root is None:
            build_root = candidate_root
        elif candidate_root != build_root:
            raise ValueError("builds do not share one output root")
        expected_library = f"{build_dir}/liblibleopard.a"
        expected_executable = f"{build_root}/c7-{name}"
        if (build["library"]["path"] != expected_library or
                build["executable"]["path"] != expected_executable):
            raise ValueError("build output path changed")
        library = validate_artifact(
            build["library"], f"{name} archive", source_root, required=live,
            check_if_present=live)
        executable = validate_artifact(
            build["executable"], f"{name} executable", source_root, required=live,
            check_if_present=live)
        log_paths: dict[str, pathlib.Path] = {}
        for label in (
            "configure_stdout", "configure_stderr", "build_stdout",
            "build_stderr", "compile_stdout", "compile_stderr", "cmake_cache",
        ):
            log_paths[label] = validate_normalized_text(
                build[label], f"{name} {label}", source_root,
                require_token=label in {
                    "configure_stdout", "build_stdout", "compile_stderr",
                    "cmake_cache"})

        prefix_text = " ".join(flags)
        core_flags = prefix_text + (
            " -fsanitize=address,undefined -fno-omit-frame-pointer"
            if sanitizer else "")
        linker_flags = (
            "-fsanitize=address,undefined -fno-omit-frame-pointer"
            if sanitizer else "")
        expected_configure = [
            str(records["cmake"]), "-S", NORMALIZATION_TOKEN, "-B",
            f"{NORMALIZATION_TOKEN}/{build_dir}", "-G", "Unix Makefiles",
            f"-DCMAKE_BUILD_TYPE={'Debug' if sanitizer else 'Release'}",
            f"-DCMAKE_C_COMPILER={records['c_compiler']}",
            f"-DCMAKE_CXX_COMPILER={records['cxx_compiler']}",
            f"-DLEO2_BACKEND_VARIANT={build['backend']}",
            "-DLEO2_BUILD_TESTS=OFF", "-DLEO2_BUILD_BENCHMARKS=OFF",
            "-DLEO2_BUILD_FUZZERS=OFF", "-DLEO2_ENABLE_CUDA=OFF",
            f"-DENABLE_OPENMP={'OFF' if sanitizer else 'ON'}",
            f"-DCMAKE_C_FLAGS={core_flags}",
            f"-DCMAKE_CXX_FLAGS={core_flags}",
            f"-DCMAKE_EXE_LINKER_FLAGS={linker_flags}",
        ]
        expected_build = [
            str(records["cmake"]), "--build",
            f"{NORMALIZATION_TOKEN}/{build_dir}", "--target", "libleopard",
            "--verbose", "--", build_argv[-1],
        ]
        expected_compile_prefix = [
            str(records["compiler"]), "-v", "-Wl,-v", "-std=c++11", "-g",
            "-Wall", "-Wextra", "-Wpedantic", "-Werror", *flags,
            f"-I{NORMALIZATION_TOKEN}",
            f'-DLEO2_C7_SOURCE_SHA256="{data["source"]["sha256"]}"',
            f'-DLEO2_C7_CORE_GIT_SHA="{core_sha}"',
            f'-DLEO2_C7_LIBRARY_SHA256="{build["library"]["sha256"]}"',
        ]
        if sanitizer:
            expected_compile_prefix += [
                "-O1", "-fsanitize=address,undefined",
                "-fno-omit-frame-pointer",
                "-DLEO2_C7_DISABLE_GLOBAL_NEW_TRACKING=1",
                '-DLEO2_C7_SANITIZER_MODE="asan-ubsan"',
                "-DLEO2_C7_REQUIRE_ASAN_UBSAN=1",
            ]
        else:
            expected_compile_prefix += ["-O2"]
        expected_compile_prefix += [
            f"{NORMALIZATION_TOKEN}/{data['source']['path']}",
            f"{NORMALIZATION_TOKEN}/{expected_library}", "-pthread",
        ]
        if not sanitizer:
            expected_compile_prefix += ["-fopenmp"]
        expected_compile_prefix += [
            "-o", f"{NORMALIZATION_TOKEN}/{expected_executable}"]
        if (configure != expected_configure or build_argv != expected_build or
                compile_argv != expected_compile_prefix):
            raise ValueError("exact normalized build command changed")

        cache = parse_cache(log_paths["cmake_cache"].read_text(encoding="utf-8"))
        expected_cache = {
            "CMAKE_AR": str(records["ar"]),
            "CMAKE_BUILD_TYPE": "Debug" if sanitizer else "Release",
            "CMAKE_COMMAND": str(records["cmake"]),
            "CMAKE_CXX_COMPILER": str(records["cxx_compiler"]),
            "CMAKE_CXX_FLAGS": core_flags,
            "CMAKE_C_COMPILER": str(records["c_compiler"]),
            "CMAKE_C_FLAGS": core_flags,
            "CMAKE_EXE_LINKER_FLAGS": linker_flags,
            "CMAKE_GENERATOR": "Unix Makefiles",
            "CMAKE_HOME_DIRECTORY": NORMALIZATION_TOKEN,
            "CMAKE_LINKER": str(records["cmake_linker"]),
            "CMAKE_MAKE_PROGRAM": str(records["gmake"]),
            "CMAKE_RANLIB": str(records["ranlib"]),
            "ENABLE_OPENMP": "OFF" if sanitizer else "ON",
            "LEO2_BACKEND_VARIANT": build["backend"],
            "LEO2_BUILD_BENCHMARKS": "OFF",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_BUILD_TESTS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
            "leopard_BINARY_DIR": f"{NORMALIZATION_TOKEN}/{build_dir}",
        }
        for key, value in expected_cache.items():
            if cache.get(key) != value:
                raise ValueError(f"CMake cache role changed: {key}")
        configure_text = log_paths["configure_stdout"].read_text(encoding="utf-8")
        family = "Clang" if sanitizer else "GNU"
        _, c_version = compiler_identity(build["c_compiler"], family)
        _, cxx_version = compiler_identity(build["cxx_compiler"], family)
        required_configure = {
            f"-- The C compiler identification is {family} {c_version}",
            f"-- The CXX compiler identification is {family} {cxx_version}",
            f"-- Check for working C compiler: {records['c_compiler']} - skipped",
            f"-- Check for working CXX compiler: {records['cxx_compiler']} - skipped",
            "-- Detecting C compile features - done",
            "-- Detecting CXX compile features - done",
            "-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success",
            "-- Found Threads: TRUE  ",
            "-- Performing Test LEO2_FLAG_MSSSE3 - Success",
            "-- Performing Test LEO2_FLAG_MAVX2 - Success",
            f"-- Build files have been written to: {NORMALIZATION_TOKEN}/{build_dir}",
        }
        if not required_configure.issubset(set(configure_text.splitlines())):
            raise ValueError("configure log lacks required semantic events")
        if log_paths["configure_stderr"].read_bytes() or log_paths[
                "build_stderr"].read_bytes():
            raise ValueError("configure or core-build stderr is nonempty")
        build_text = log_paths["build_stdout"].read_text(encoding="utf-8")
        built_sources = set(re.findall(
            r"Building CXX object [^\n ]*/([^/ ]+\.cpp)\.o", build_text))
        if built_sources != CORE_SOURCES:
            raise ValueError("core build source set changed")
        for required in (
            str(records["cmake"]), str(records["cxx_compiler"]),
            str(records["gmake"]), str(records["ar"]), str(records["ranlib"]),
            *flags, "Built target leopard2_backend_ssse3",
            "Built target leopard2_backend_avx2",
            "Linking CXX static library liblibleopard.a", "Built target libleopard",
        ):
            if required not in build_text:
                raise ValueError("core build log lost an exact role or flag")
        compile_stdout = log_paths["compile_stdout"].read_text(encoding="utf-8")
        compile_stderr = log_paths["compile_stderr"].read_text(encoding="utf-8")
        linker_first = build["standalone_linker"]["version"].splitlines()[0]
        if compile_stdout != f"{linker_first}\n":
            raise ValueError("standalone linker output changed")
        for required in (
            "c7_exact_low.cpp", "liblibleopard.a", "-std=c++11",
            str(records["standalone_linker"]), *flags,
        ):
            if required not in compile_stderr:
                raise ValueError("standalone compile closure changed")

        instrumentation = build["instrumentation"]
        required_instrumentation = {
            "archive_members", "core_archive_counts",
            "core_archive_member_counts", "core_archive_symbol_scan",
            "executable_counts", "executable_symbol_scan",
            "required_compile_macro",
        }
        if set(instrumentation) != required_instrumentation:
            raise ValueError("instrumentation schema changed")
        zero_counts = {"asan_lines": 0, "ubsan_lines": 0}
        expected_executable_counts = (
            EXPECTED_EXECUTABLE_SANITIZER_COUNTS if sanitizer else zero_counts)
        expected_archive_counts = (
            EXPECTED_ARCHIVE_SANITIZER_COUNTS if sanitizer else zero_counts)
        expected_members = (
            EXPECTED_ARCHIVE_MEMBER_COUNTS if sanitizer else
            {member: dict(zero_counts) for member in ARCHIVE_MEMBERS})
        if (instrumentation["required_compile_macro"] is not sanitizer or
                instrumentation["archive_members"] != list(ARCHIVE_MEMBERS) or
                instrumentation["executable_counts"] != expected_executable_counts or
                instrumentation["core_archive_counts"] != expected_archive_counts or
                instrumentation["core_archive_member_counts"] != expected_members):
            raise ValueError("sanitizer summary or member attribution changed")
        executable_scan = validate_normalized_text(
            instrumentation["executable_symbol_scan"],
            f"{name} executable scan", source_root, require_token=False)
        archive_scan = validate_normalized_text(
            instrumentation["core_archive_symbol_scan"],
            f"{name} archive scan", source_root, require_token=False)
        validate_symbol_scan(
            executable_scan.read_text(encoding="utf-8"), executable.name,
            archive=False, expected_counts=expected_executable_counts,
            expected_members={})
        validate_symbol_scan(
            archive_scan.read_text(encoding="utf-8"), library.name,
            archive=True, expected_counts=expected_archive_counts,
            expected_members=expected_members)
        if live:
            live_nm(build["nm"], executable, executable_scan)
            live_nm(build["nm"], library, archive_scan)

        closure = build["source_closure"]
        if not isinstance(closure, list) or len(closure) < 12:
            raise ValueError("core source closure is incomplete")
        closure_paths: list[str] = []
        for index, entry in enumerate(closure):
            validate_artifact(
                entry, f"{name} source closure {index}", source_root, required=True)
            relative = entry["path"]
            if pathlib.Path(relative).is_absolute():
                raise ValueError("source closure path is not portable")
            validate_git_artifact(core_sha, entry, relative, "source closure")
            closure_paths.append(relative)
        if closure_paths != sorted(set(closure_paths)):
            raise ValueError("source closure order or uniqueness changed")
        if not {
            "CMakeLists.txt", "cmake/leopardConfig.cmake.in", "leopard2.cpp",
            "LeopardFF8.cpp", "LeopardFF16.cpp",
        }.issubset(closure_paths):
            raise ValueError("source closure lacks required configure/core inputs")
        by_name[name] = build

    reproducibility = data["reproducibility"]
    if not isinstance(reproducibility, dict) or set(reproducibility) != {
            "fingerprints", "prefix_map_target"} or reproducibility[
            "prefix_map_target"] != PREFIX_MAP_TARGET:
        raise ValueError("reproducibility record changed")
    expected_fingerprints = {
        name: {
            "library_sha256": by_name[name]["library"]["sha256"],
            "executable_sha256": by_name[name]["executable"]["sha256"],
        }
        for name in BUILD_NAMES
    }
    if reproducibility["fingerprints"] != expected_fingerprints:
        raise ValueError("reproducibility fingerprints changed")

    runs = data["runs"]
    expected_runs = [*BUILD_NAMES, "smoke-nonauthoritative"]
    if not isinstance(runs, list) or [run.get("name") for run in runs] != expected_runs:
        raise ValueError("run matrix changed")
    correctness_cpus: list[int] = []
    for run in runs:
        if set(run) != {
            "argv", "argv_source_root_tokens", "build", "environment", "kind",
            "name", "observed_affinity", "requested_cpu", "result", "stderr",
            "stdout",
        }:
            raise ValueError("run record schema changed")
        name = run["name"]
        smoke = name == "smoke-nonauthoritative"
        build_name = "auto" if smoke else name
        cpu = run["requested_cpu"]
        if type(cpu) is not int or cpu < 0 or run["observed_affinity"] != [cpu]:
            raise ValueError("run CPU or observed affinity changed")
        if not smoke:
            correctness_cpus.append(cpu)
        argv = validate_argv(
            run["argv"], run["argv_source_root_tokens"], f"{name} run",
            require_token=True)
        result = validate_normalized_text(
            run["result"], f"{name} result", source_root, require_token=False)
        stdout = validate_normalized_text(
            run["stdout"], f"{name} stdout", source_root, require_token=False)
        stderr = validate_normalized_text(
            run["stderr"], f"{name} stderr", source_root, require_token=False)
        expected_argv = [
            str(taskset), "-c", str(cpu),
            f"{NORMALIZATION_TOKEN}/{by_name[build_name]['executable']['path']}",
            "--backend", "auto" if build_name == "asan-ubsan" else
            by_name[build_name]["backend"],
            (run["result"]["path"] if pathlib.Path(run["result"]["path"]).is_absolute()
             else f"{NORMALIZATION_TOKEN}/{run['result']['path']}"),
            "--benchmark-smoke" if smoke else "--correctness-only",
        ]
        if argv != expected_argv:
            raise ValueError("normalized run command changed")
        if stdout.read_bytes() or stderr.read_text(encoding="utf-8") != (
                "C7 benchmark 1/1\n" if smoke else ""):
            raise ValueError("run log semantics changed")
        child = json.loads(result.read_text(encoding="utf-8"))
        build = by_name[build_name]
        if (child.get("source_sha256") != data["source"]["sha256"] or
                child.get("core_git_sha") != core_sha or
                child.get("library_sha256") != build["library"]["sha256"]):
            raise ValueError("run result is not bound to build inputs")
        if smoke:
            if (child.get("timing_scope") != "non-authoritative-smoke" or len(
                    child.get("benchmarks", [])) != 1 or
                    child.get("correctness") != EXPECTED_CORRECTNESS or
                    child.get("production_constructor_rejected") is not True or
                    child.get("requested_backend") != "auto"):
                raise ValueError("smoke result changed")
        else:
            validate_child(
                child, "auto" if build_name == "asan-ubsan" else build_name,
                [cpu], build_name == "asan-ubsan")
    if len(correctness_cpus) != 5 or len(set(correctness_cpus)) != 5:
        raise ValueError("correctness runs did not use five distinct CPUs")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=pathlib.Path)
    parser.add_argument(
        "--live", action="store_true",
        help="require exact recorded tools/build outputs and replay nm scans")
    arguments = parser.parse_args()
    data = json.loads(arguments.manifest.read_text(encoding="utf-8"))
    validate_manifest(data, live=arguments.live)
    print("C7 evidence validation passed (live)" if arguments.live else
          "C7 evidence validation passed (portable)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
