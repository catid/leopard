#!/usr/bin/env python3
"""Rebuild and run the hash-bound C7 correctness/smoke matrix.

The retained checkpoint deliberately uses ordinary CPUs and labels its only
timed cell non-authoritative.  A later authoritative runner may reuse the C++
harness, but must produce a separate manifest against the integrated core.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import pathlib
import re
import shlex
import shutil
import subprocess
import sys
from typing import Any, Iterable


ROOT = pathlib.Path(__file__).resolve().parents[4]
HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "c7_exact_low.cpp"
VALIDATOR = HERE / "validate_evidence.py"
BACKENDS = ("scalar", "ssse3", "avx2", "auto")
BUILD_NAMES = (*BACKENDS, "asan-ubsan")
PROGRAM_ROLES = (
    "ar", "c_compiler", "cmake", "cmake_linker", "compiler",
    "cxx_compiler", "gmake", "launcher_python", "link_driver", "nm", "ranlib",
    "standalone_linker",
)
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
GIT_SHA_RE = re.compile(r"[0-9a-f]{40}\Z")
NORMALIZATION_SCHEMA = "leopard2-source-root-prefix/v1"
NORMALIZATION_TOKEN = "${LEO2_SOURCE_ROOT}"
PREFIX_MAP_TARGET = "LEO2_SOURCE_ROOT"
PEER_ATTESTATION_SCHEMA = "leopard2-c7-peer-reproducibility/v1"
PREFIX_MAP_OPTIONS = (
    "-ffile-prefix-map", "-fdebug-prefix-map", "-fmacro-prefix-map")
EXPECTED_EXECUTABLE_SANITIZER_COUNTS = {
    "asan_lines": 320, "ubsan_lines": 54}
EXPECTED_ARCHIVE_SANITIZER_COUNTS = {
    "asan_lines": 329, "ubsan_lines": 87}
EXPECTED_ARCHIVE_MEMBER_COUNTS = {
    "leopard.cpp.o": {"asan_lines": 13, "ubsan_lines": 7},
    "leopard2.cpp.o": {"asan_lines": 141, "ubsan_lines": 15},
    "Leopard2Backend.cpp.o": {"asan_lines": 35, "ubsan_lines": 9},
    "Leopard2BackendScalar.cpp.o": {"asan_lines": 11, "ubsan_lines": 6},
    "Leopard2CpuFeatures.cpp.o": {"asan_lines": 9, "ubsan_lines": 5},
    "Leopard2Plan.cpp.o": {"asan_lines": 7, "ubsan_lines": 5},
    "LeopardCommon.cpp.o": {"asan_lines": 13, "ubsan_lines": 6},
    "LeopardFF16.cpp.o": {"asan_lines": 31, "ubsan_lines": 10},
    "LeopardFF8.cpp.o": {"asan_lines": 31, "ubsan_lines": 9},
    "Leopard2BackendSSSE3.cpp.o": {"asan_lines": 18, "ubsan_lines": 8},
    "Leopard2BackendAVX2.cpp.o": {"asan_lines": 20, "ubsan_lines": 7},
}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def relative(path: pathlib.Path) -> str:
    path = path.resolve()
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def resolve_program(name: str) -> pathlib.Path:
    found = shutil.which(name)
    if not found:
        raise RuntimeError(f"required program is unavailable: {name}")
    return pathlib.Path(found).resolve()


def resolve_first(names: Iterable[str]) -> pathlib.Path:
    for name in names:
        found = shutil.which(name)
        if found:
            # Preserve a C++ driver basename such as clang++-18.  Resolving
            # that symlink to the shared `clang` binary changes driver mode at
            # link time and silently drops the C++ runtime.
            return pathlib.Path(found).absolute()
    raise RuntimeError(
        f"none of the required programs are available: {', '.join(names)}")


def artifact(path: pathlib.Path) -> dict[str, Any]:
    return {
        "path": relative(path),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


def normalize_source_root(text: str) -> tuple[str, int]:
    """Replace only ROOT itself or ROOT followed by a path separator."""
    pattern = re.compile(re.escape(str(ROOT)) + r"(?=\Z|[/\\=\s\"'])")
    return pattern.subn(lambda _match: NORMALIZATION_TOKEN, text)


def normalize_argv(argv: list[str]) -> tuple[list[str], int]:
    normalized: list[str] = []
    count = 0
    for argument in argv:
        value, replacements = normalize_source_root(argument)
        normalized.append(value)
        count += replacements
    return normalized, count


def portable_compiler_argv(
    argv: list[str], working_directory: pathlib.Path,
) -> list[str]:
    """Spell checkout inputs relative to CMake's stable build directory.

    Clang's ASan global metadata does not honor prefix maps for an absolute
    source argv.  Relative source and include spellings are reproducible and
    still leave the prefix-map options in place for debug and macro paths.
    """
    root = str(ROOT)
    prefix = root + os.sep
    relative_root = os.path.relpath(ROOT, working_directory)
    rewritten: list[str] = []
    for argument in argv:
        if argument == root:
            rewritten.append(relative_root)
        elif argument.startswith(prefix):
            rewritten.append(os.path.relpath(argument, working_directory))
        elif argument == f"-I{root}":
            rewritten.append(f"-I{relative_root}")
        elif argument.startswith(f"-I{prefix}"):
            rewritten.append(
                "-I" + os.path.relpath(argument[2:], working_directory))
        else:
            rewritten.append(argument)
    return rewritten


def compiler_launch(argv: list[str]) -> int:
    if not argv:
        raise SystemExit("--compiler-launch requires the compiler argv")
    completed = subprocess.run(
        portable_compiler_argv(argv, pathlib.Path.cwd()), check=False)
    return completed.returncode


def normalized_text_artifact(path: pathlib.Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="strict")
    normalized, count = normalize_source_root(text)
    path.write_text(normalized, encoding="utf-8")
    record = artifact(path)
    record["source_root_tokens"] = count
    return record


def prefix_map_flags(
    c_compiler: pathlib.Path, cxx_compiler: pathlib.Path,
) -> list[str]:
    result: list[str] = []
    for index, option in enumerate(PREFIX_MAP_OPTIONS):
        flag = f"{option}={ROOT}={PREFIX_MAP_TARGET}"
        supported = []
        for compiler, language in ((c_compiler, "c"), (cxx_compiler, "c++")):
            completed = subprocess.run(
                [str(compiler), "-Werror", flag, "-x", language, "-c", "-",
                 "-o", os.devnull],
                input=b"", stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False,
            )
            supported.append(completed.returncode == 0)
        if index < 2 and not all(supported):
            raise RuntimeError(
                f"required reproducibility flag is unsupported: {flag}")
        if all(supported):
            result.append(flag)
    if len(result) < 2:
        raise AssertionError("required prefix-map flags were not selected")
    return result


def committed_sha256(commit: str, path: pathlib.Path) -> str:
    relative_path = path.resolve().relative_to(ROOT).as_posix()
    completed = subprocess.run(
        ["git", "show", f"{commit}:{relative_path}"], cwd=ROOT,
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"cannot read {relative_path} from {commit}: "
            f"{completed.stderr.decode('utf-8', errors='replace').strip()}")
    return hashlib.sha256(completed.stdout).hexdigest()


def require_committed_artifact(commit: str, record: dict[str, Any]) -> None:
    path = (ROOT / record["path"]).resolve()
    if committed_sha256(commit, path) != record["sha256"]:
        raise RuntimeError(
            f"working-tree artifact is not bound to {commit}: {record['path']}")


def run_logged(
    argv: list[str], stdout_path: pathlib.Path, stderr_path: pathlib.Path,
    *, env_additions: dict[str, str] | None = None,
) -> None:
    environment = os.environ.copy()
    if env_additions:
        environment.update(env_additions)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        completed = subprocess.run(
            argv, cwd=ROOT, env=environment, stdout=stdout, stderr=stderr,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {shlex.join(argv)}; "
            f"see {stderr_path}")


def program_record(path: pathlib.Path) -> dict[str, Any]:
    version = subprocess.run(
        [str(path), "--version"], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    ).stdout
    return {
        "path": str(path),
        "sha256": sha256(path),
        "version": version,
    }


def cmake_cache_value(path: pathlib.Path, key: str) -> str:
    prefix = f"{key}:"
    matches = []
    for line in path.read_text(encoding="utf-8", errors="strict").splitlines():
        if line.startswith(prefix) and "=" in line:
            matches.append(line.split("=", 1)[1])
    if len(matches) != 1 or not matches[0]:
        raise RuntimeError(f"CMake cache key is absent or ambiguous: {key}")
    return matches[0]


def compiler_program(compiler: pathlib.Path, name: str) -> pathlib.Path:
    output = subprocess.run(
        [str(compiler), f"-print-prog-name={name}"], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout.strip()
    if not output:
        raise RuntimeError(f"compiler did not resolve {name}: {compiler}")
    path = pathlib.Path(output)
    if path.is_absolute():
        return path.absolute()
    found = shutil.which(output)
    if not found:
        raise RuntimeError(f"compiler-selected program is unavailable: {output}")
    return pathlib.Path(found).absolute()


def filtered_symbol_scan(
    nm: pathlib.Path, target: pathlib.Path, output: pathlib.Path,
    archive_members: Iterable[str],
) -> tuple[dict[str, int], dict[str, dict[str, int]]]:
    text = subprocess.run(
        [str(nm), "--print-file-name", target.name], cwd=target.parent,
        check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout
    lines = [line for line in text.splitlines()
             if "__asan_" in line or "__ubsan_" in line]
    counts = {
        "asan_lines": sum("__asan_" in line for line in lines),
        "ubsan_lines": sum("__ubsan_" in line for line in lines),
    }
    members = {
        member: {"asan_lines": 0, "ubsan_lines": 0}
        for member in archive_members
    }
    for line in lines:
        for member in members:
            if f":{member}:" in line:
                members[member]["asan_lines"] += "__asan_" in line
                members[member]["ubsan_lines"] += "__ubsan_" in line
                break
    normalized = "\n".join(lines) + ("\n" if lines else "")
    if str(ROOT) in normalized:
        raise RuntimeError("nm output leaked the source checkout root")
    output.write_text(normalized, encoding="utf-8")
    return counts, members


def dependency_closure(build_dir: pathlib.Path) -> list[dict[str, Any]]:
    paths: set[pathlib.Path] = {
        ROOT / "CMakeLists.txt", ROOT / "cmake/leopardConfig.cmake.in"}
    for dependency in build_dir.rglob("*.o.d"):
        text = dependency.read_text(encoding="utf-8", errors="strict")
        flattened = text.replace("\\\n", " ")
        if ":" not in flattened:
            raise RuntimeError(f"malformed dependency file: {dependency}")
        for token in shlex.split(flattened.split(":", 1)[1]):
            candidate = pathlib.Path(token)
            if not candidate.is_absolute():
                candidate = (build_dir / candidate).resolve()
            else:
                candidate = candidate.resolve()
            try:
                candidate.relative_to(ROOT)
            except ValueError:
                continue
            if candidate.is_file() and "build" not in candidate.relative_to(ROOT).parts:
                paths.add(candidate)
    if len(paths) < 10:
        raise RuntimeError("core dependency closure is unexpectedly small")
    return [artifact(path) for path in sorted(paths)]


def cmake_build(
    backend: str, build_root: pathlib.Path, results: pathlib.Path,
    jobs: int, c_compiler: pathlib.Path, cxx_compiler: pathlib.Path,
    sanitizer: bool,
) -> dict[str, Any]:
    name = "asan-ubsan" if sanitizer else backend
    build_dir = build_root / f"core-{name}"
    if build_dir.exists():
        shutil.rmtree(build_dir)
    configure_stdout = results / "logs" / f"{name}-configure.stdout.txt"
    configure_stderr = results / "logs" / f"{name}-configure.stderr.txt"
    build_stdout = results / "logs" / f"{name}-core-build.stdout.txt"
    build_stderr = results / "logs" / f"{name}-core-build.stderr.txt"
    cache_output = results / "logs" / f"{name}-CMakeCache.txt"
    reproducibility_flags = prefix_map_flags(c_compiler, cxx_compiler)
    compile_flags = list(reproducibility_flags)
    if sanitizer:
        compile_flags.extend([
            "-fsanitize=address,undefined", "-fno-omit-frame-pointer"])
    flag_text = " ".join(compile_flags)
    linker_flags = (
        "-fsanitize=address,undefined -fno-omit-frame-pointer"
        if sanitizer else "")
    launcher_python = pathlib.Path(sys.executable).resolve()
    launcher = f"{launcher_python};{pathlib.Path(__file__).resolve()};--compiler-launch"
    configure_argv = [
        str(resolve_program("cmake")), "-S", str(ROOT), "-B", str(build_dir),
        "-G", "Unix Makefiles",
        f"-DCMAKE_BUILD_TYPE={'Debug' if sanitizer else 'Release'}",
        f"-DCMAKE_C_COMPILER={c_compiler}",
        f"-DCMAKE_CXX_COMPILER={cxx_compiler}",
        f"-DLEO2_BACKEND_VARIANT={'auto' if sanitizer else backend}",
        "-DLEO2_BUILD_TESTS=OFF", "-DLEO2_BUILD_BENCHMARKS=OFF",
        "-DLEO2_BUILD_FUZZERS=OFF", "-DLEO2_ENABLE_CUDA=OFF",
        f"-DENABLE_OPENMP={'OFF' if sanitizer else 'ON'}",
        f"-DCMAKE_C_FLAGS={flag_text}",
        f"-DCMAKE_CXX_FLAGS={flag_text}",
        f"-DCMAKE_EXE_LINKER_FLAGS={linker_flags}",
        f"-DCMAKE_CXX_COMPILER_LAUNCHER={launcher}",
    ]
    run_logged(configure_argv, configure_stdout, configure_stderr)
    cache_path = build_dir / "CMakeCache.txt"
    shutil.copyfile(cache_path, cache_output)
    cache_compilers = {
        "CMAKE_C_COMPILER": str(c_compiler),
        "CMAKE_CXX_COMPILER": str(cxx_compiler),
    }
    for key, expected in cache_compilers.items():
        if cmake_cache_value(cache_path, key) != expected:
            raise RuntimeError(f"CMake selected an unexpected {key}")
    make_program = pathlib.Path(
        cmake_cache_value(cache_path, "CMAKE_MAKE_PROGRAM")).absolute()
    archive_program = pathlib.Path(
        cmake_cache_value(cache_path, "CMAKE_AR")).absolute()
    ranlib_program = pathlib.Path(
        cmake_cache_value(cache_path, "CMAKE_RANLIB")).absolute()
    cmake_linker = pathlib.Path(
        cmake_cache_value(cache_path, "CMAKE_LINKER")).absolute()
    build_argv = [
        str(resolve_program("cmake")), "--build", str(build_dir),
        "--target", "libleopard", "--verbose", "--", f"-j{jobs}",
    ]
    run_logged(build_argv, build_stdout, build_stderr)
    library = build_dir / "liblibleopard.a"
    if not library.is_file():
        raise RuntimeError(f"missing library output: {library}")
    normalized_configure, configure_tokens = normalize_argv(configure_argv)
    normalized_build, build_tokens = normalize_argv(build_argv)
    return {
        "name": name,
        "backend": "auto" if sanitizer else backend,
        "sanitizer": sanitizer,
        "configure_argv": normalized_configure,
        "build_argv": normalized_build,
        "argv_source_root_tokens": {
            "configure": configure_tokens, "build": build_tokens,
        },
        "configure_stdout": normalized_text_artifact(configure_stdout),
        "configure_stderr": normalized_text_artifact(configure_stderr),
        "build_stdout": normalized_text_artifact(build_stdout),
        "build_stderr": normalized_text_artifact(build_stderr),
        "cmake_cache": normalized_text_artifact(cache_output),
        "cmake": program_record(resolve_program("cmake")),
        "c_compiler": program_record(c_compiler),
        "cxx_compiler": program_record(cxx_compiler),
        "gmake": program_record(make_program),
        "ar": program_record(archive_program),
        "ranlib": program_record(ranlib_program),
        "cmake_linker": program_record(cmake_linker),
        "launcher_python": program_record(launcher_python),
        "library": artifact(library),
        "source_closure": dependency_closure(build_dir),
        "build_dir": relative(build_dir),
        "prefix_map_flags": [normalize_source_root(flag)[0]
                             for flag in reproducibility_flags],
    }


def compile_experiment(
    build: dict[str, Any], build_root: pathlib.Path, results: pathlib.Path,
    compiler: pathlib.Path, core_git_sha: str,
) -> dict[str, Any]:
    name = build["name"]
    executable = build_root / f"c7-{name}"
    stdout_path = results / "logs" / f"{name}-compile.stdout.txt"
    stderr_path = results / "logs" / f"{name}-compile.stderr.txt"
    source_hash = sha256(SOURCE)
    library_path = ROOT / build["library"]["path"]
    reproducibility_flags = [
        flag.replace(NORMALIZATION_TOKEN, str(ROOT), 1)
        for flag in build["prefix_map_flags"]
    ]
    argv = [
        str(compiler), "-v", "-Wl,-v", "-std=c++11", "-g", "-Wall", "-Wextra",
        "-Wpedantic", "-Werror", *reproducibility_flags, f"-I{ROOT}",
        f'-DLEO2_C7_SOURCE_SHA256="{source_hash}"',
        f'-DLEO2_C7_CORE_GIT_SHA="{core_git_sha}"',
        f'-DLEO2_C7_LIBRARY_SHA256="{build["library"]["sha256"]}"',
    ]
    if build["sanitizer"]:
        argv.extend([
            "-O1", "-fsanitize=address,undefined", "-fno-omit-frame-pointer",
            "-DLEO2_C7_DISABLE_GLOBAL_NEW_TRACKING=1",
            '-DLEO2_C7_SANITIZER_MODE="asan-ubsan"',
            "-DLEO2_C7_REQUIRE_ASAN_UBSAN=1",
        ])
    else:
        argv.append("-O2")
    argv.extend([relative(SOURCE), relative(library_path), "-pthread"])
    if not build["sanitizer"]:
        argv.append("-fopenmp")
    argv.extend(["-o", relative(executable)])
    run_logged(argv, stdout_path, stderr_path)
    nm = resolve_program("nm")
    archive_program = pathlib.Path(build["ar"]["path"])
    library_members = subprocess.run(
        [str(archive_program), "t", relative(library_path)], cwd=ROOT,
        check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout.splitlines()
    if len(library_members) != len(set(library_members)) or not library_members:
        raise RuntimeError("static archive member list is invalid")
    executable_scan = results / "logs" / f"{name}-nm-executable-sanitizers.txt"
    archive_scan = results / "logs" / f"{name}-nm-core-archive-sanitizers.txt"
    executable_counts, executable_members = filtered_symbol_scan(
        nm, executable, executable_scan, ())
    archive_counts, archive_member_counts = filtered_symbol_scan(
        nm, library_path, archive_scan, library_members)
    if executable_members:
        raise RuntimeError("executable symbol scan unexpectedly has archive members")
    expected_executable = (EXPECTED_EXECUTABLE_SANITIZER_COUNTS
                           if build["sanitizer"] else
                           {"asan_lines": 0, "ubsan_lines": 0})
    expected_archive = (EXPECTED_ARCHIVE_SANITIZER_COUNTS
                        if build["sanitizer"] else
                        {"asan_lines": 0, "ubsan_lines": 0})
    if executable_counts != expected_executable or archive_counts != expected_archive:
        raise RuntimeError("sanitizer symbol family/count proof changed")
    expected_members = (
        EXPECTED_ARCHIVE_MEMBER_COUNTS if build["sanitizer"] else
        {member: {"asan_lines": 0, "ubsan_lines": 0}
         for member in library_members})
    if archive_member_counts != expected_members:
        raise RuntimeError("sanitizer archive member attribution changed")
    instrumentation = {
        "required_compile_macro": bool(build["sanitizer"]),
        "archive_members": library_members,
        "core_archive_symbol_scan": normalized_text_artifact(archive_scan),
        "core_archive_counts": archive_counts,
        "core_archive_member_counts": archive_member_counts,
        "executable_symbol_scan": normalized_text_artifact(executable_scan),
        "executable_counts": executable_counts,
    }
    standalone_linker = compiler_program(compiler, "ld")
    result = dict(build)
    normalized_compile, compile_tokens = normalize_argv(argv)
    token_counts = dict(build["argv_source_root_tokens"])
    token_counts["compile"] = compile_tokens
    result.update({
        "compiler": program_record(compiler),
        "link_driver": program_record(compiler),
        "standalone_linker": program_record(standalone_linker),
        "compile_argv": normalized_compile,
        "argv_source_root_tokens": token_counts,
        "compile_stdout": normalized_text_artifact(stdout_path),
        "compile_stderr": normalized_text_artifact(stderr_path),
        "executable": artifact(executable),
        "instrumentation": instrumentation,
        "nm": program_record(nm),
    })
    return result


def run_one(
    build: dict[str, Any], results: pathlib.Path, cpu: int, smoke: bool,
) -> dict[str, Any]:
    name = "smoke-nonauthoritative" if smoke else build["name"]
    result_path = results / f"{name}.json"
    stdout_path = results / "logs" / f"{name}-run.stdout.txt"
    stderr_path = results / "logs" / f"{name}-run.stderr.txt"
    executable = ROOT / build["executable"]["path"]
    backend = "auto" if build["sanitizer"] else build["backend"]
    argv = [
        str(resolve_program("taskset")), "-c", str(cpu),
        str(executable), "--backend", backend, str(result_path),
        "--benchmark-smoke" if smoke else "--correctness-only",
    ]
    environment = {
        "LC_ALL": "C", "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
    }
    if build["sanitizer"]:
        environment.update({
            "ASAN_OPTIONS": "detect_leaks=1:halt_on_error=1",
            "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1",
        })
    run_logged(argv, stdout_path, stderr_path, env_additions=environment)
    data = json.loads(result_path.read_text(encoding="utf-8"))
    if data.get("affinity") != [cpu]:
        raise RuntimeError(f"child affinity mismatch for {name}")
    if data.get("source_sha256") != sha256(SOURCE):
        raise RuntimeError(f"child source fingerprint mismatch for {name}")
    if data.get("library_sha256") != build["library"]["sha256"]:
        raise RuntimeError(f"child library fingerprint mismatch for {name}")
    normalized_run, run_tokens = normalize_argv(argv)
    return {
        "name": name,
        "build": build["name"],
        "kind": "non-authoritative-smoke" if smoke else "correctness",
        "requested_cpu": cpu,
        "argv": normalized_run,
        "argv_source_root_tokens": run_tokens,
        "environment": environment,
        "result": normalized_text_artifact(result_path),
        "stdout": normalized_text_artifact(stdout_path),
        "stderr": normalized_text_artifact(stderr_path),
        "observed_affinity": data["affinity"],
    }


def parse_cpus(text: str) -> list[int]:
    result = [int(item) for item in text.split(",")]
    if len(result) != 5 or len(set(result)) != 5 or any(cpu < 0 for cpu in result):
        raise argparse.ArgumentTypeError("--cpus requires five distinct CPUs")
    return result


def default_cpus(allowed: Iterable[int]) -> list[int]:
    selected = sorted(set(allowed))[:5]
    if len(selected) != 5:
        raise ValueError("C7 matrix requires at least five allowed CPUs")
    return selected


def validate_cpus(
    cpus: list[int], smoke_cpu: int, allowed: Iterable[int],
) -> None:
    allowed_set = set(allowed)
    if len(cpus) != 5 or len(set(cpus)) != 5 or any(
            type(cpu) is not int or cpu < 0 for cpu in cpus):
        raise ValueError("C7 correctness matrix requires five distinct CPUs")
    if any(cpu not in allowed_set for cpu in [*cpus, smoke_cpu]):
        raise ValueError("requested CPU is outside process affinity")


def reproducibility_fingerprints(builds: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        build["name"]: {
            "library_sha256": build["library"]["sha256"],
            "executable_sha256": build["executable"]["sha256"],
        }
        for build in builds
    }


def manifest_program_records(manifest: dict[str, Any]) -> dict[str, Any]:
    return {
        "taskset": manifest["taskset"],
        "builds": {
            build["name"]: {role: build[role] for role in PROGRAM_ROLES}
            for build in manifest["builds"]
        },
    }


def require_reproducible_peer(current: dict[str, Any], peer: dict[str, Any]) -> None:
    for key in ("schema", "core_git_sha"):
        if peer.get(key) != current.get(key):
            raise RuntimeError(f"A/B manifest {key} differs")
    if (peer.get("source") != current.get("source") or
            peer.get("runner") != current.get("runner") or
            peer.get("validator") != current.get("validator")):
        raise RuntimeError("A/B manifests do not bind the same committed tooling")
    try:
        if manifest_program_records(peer) != manifest_program_records(current):
            raise RuntimeError("A/B manifests do not bind the same exact programs")
    except (KeyError, TypeError) as error:
        raise RuntimeError("A/B peer program records are incomplete") from error
    left = current.get("reproducibility", {}).get("fingerprints")
    right = peer.get("reproducibility", {}).get("fingerprints")
    if left != right or set(left or ()) != set(BUILD_NAMES):
        raise RuntimeError("A/B archive or executable hashes differ")


def trusted_validate_manifest(
    manifest_path: pathlib.Path, source_root: pathlib.Path, *, live: bool,
) -> None:
    argv = [
        sys.executable, str(VALIDATOR), "--source-root", str(source_root),
        "--require-checkout-head",
    ]
    if live:
        argv.append("--live")
    argv.append(str(manifest_path))
    completed = subprocess.run(
        argv, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, check=False, env={
            **os.environ, "PYTHONDONTWRITEBYTECODE": "1",
        })
    if completed.returncode != 0:
        detail = completed.stderr.strip().splitlines()
        suffix = detail[-1] if detail else "no validator diagnostic"
        raise RuntimeError(
            f"A/B peer failed {'live' if live else 'portable'} validation: "
            f"{suffix}")


def authenticate_peer(
    current: dict[str, Any], peer_manifest_path: pathlib.Path,
    peer_root: pathlib.Path,
) -> dict[str, Any]:
    trusted_validate_manifest(peer_manifest_path, peer_root, live=False)
    peer = json.loads(peer_manifest_path.read_text(encoding="utf-8"))
    require_reproducible_peer(current, peer)
    # Portable validation and exact-program equality happen before live replay
    # so an untrusted manifest cannot redirect a tool invocation.
    trusted_validate_manifest(peer_manifest_path, peer_root, live=True)
    return peer


def checkout_artifact_path(
    record: dict[str, Any], source_root: pathlib.Path, label: str,
) -> pathlib.Path:
    if not isinstance(record, dict) or not isinstance(record.get("path"), str):
        raise RuntimeError(f"{label} artifact record is malformed")
    pure = pathlib.PurePosixPath(record["path"])
    if (pure.is_absolute() or pure.as_posix() != record["path"] or
            "\\" in record["path"] or ":" in record["path"] or
            any(part in ("", ".", "..") for part in pure.parts)):
        raise RuntimeError(f"{label} artifact path is not checkout-relative")
    root = source_root.resolve()
    path = root.joinpath(*pure.parts)
    try:
        path.resolve(strict=False).relative_to(root)
    except ValueError as error:
        raise RuntimeError(f"{label} artifact path escapes its checkout") from error
    if (not path.is_file() or type(record.get("bytes")) is not int or
            path.stat().st_size != record["bytes"] or
            not isinstance(record.get("sha256"), str) or
            sha256(path) != record["sha256"]):
        raise RuntimeError(f"{label} artifact bytes differ from its record")
    return path


def normalized_text_records(value: Any) -> Iterable[dict[str, Any]]:
    if isinstance(value, dict):
        if set(value) == {"bytes", "path", "sha256", "source_root_tokens"}:
            yield value
        for child in value.values():
            yield from normalized_text_records(child)
    elif isinstance(value, list):
        for child in value:
            yield from normalized_text_records(child)


def require_no_root_bytes(
    manifest: dict[str, Any], source_root: pathlib.Path,
    forbidden_roots: Iterable[pathlib.Path],
) -> dict[str, int]:
    needles = [str(path.resolve()).encode("utf-8") for path in forbidden_roots]
    serialized = json.dumps(manifest, sort_keys=True).encode("utf-8")
    if any(needle in serialized for needle in needles):
        raise RuntimeError("normalized manifest leaked an A/B checkout root")
    text_records = list(normalized_text_records(manifest))
    for index, record in enumerate(text_records):
        path = checkout_artifact_path(
            record, source_root, f"normalized text {index}")
        contents = path.read_bytes()
        if any(needle in contents for needle in needles):
            raise RuntimeError("normalized retained text leaked an A/B checkout root")
    by_name = {build["name"]: build for build in manifest["builds"]}
    for name in BUILD_NAMES:
        for key in ("library", "executable"):
            path = checkout_artifact_path(
                by_name[name][key], source_root, f"{name} {key}")
            contents = path.read_bytes()
            if any(needle in contents for needle in needles):
                raise RuntimeError("A/B binary leaked a checkout root")
    return {
        "normalized_text_records": len(text_records),
        "archives": len(BUILD_NAMES),
        "executables": len(BUILD_NAMES),
    }


def write_peer_attestation(
    peer: dict[str, Any], peer_manifest_path: pathlib.Path,
    peer_scan: dict[str, int], output: pathlib.Path,
    forbidden_roots: Iterable[pathlib.Path],
) -> dict[str, Any]:
    by_name = {build["name"]: build for build in peer["builds"]}
    normalized = sorted(
        (dict(record) for record in normalized_text_records(peer)),
        key=lambda record: record["path"])
    report = {
        "schema": PEER_ATTESTATION_SCHEMA,
        "status": "pass",
        "core_git_sha": peer["core_git_sha"],
        "peer_manifest": {
            "bytes": peer_manifest_path.stat().st_size,
            "sha256": sha256(peer_manifest_path),
        },
        "tooling": {
            key: peer[key] for key in ("source", "runner", "validator")
        },
        "fingerprints": peer["reproducibility"]["fingerprints"],
        "binary_artifacts": {
            name: {
                key: by_name[name][key] for key in ("library", "executable")
            }
            for name in BUILD_NAMES
        },
        "program_records_sha256": canonical_json_sha256(
            manifest_program_records(peer)),
        "source_closures_sha256": canonical_json_sha256({
            name: by_name[name]["source_closure"] for name in BUILD_NAMES
        }),
        "normalized_text_records_sha256": canonical_json_sha256(normalized),
        "runs_sha256": canonical_json_sha256(peer["runs"]),
        "root_scan": peer_scan,
        "checks": {
            "portable_semantics": "pass",
            "live_tools_and_outputs": "pass",
            "git_toplevel_and_head": "pass",
            "program_identity_match": "pass",
            "binary_and_text_root_scan": "pass",
        },
    }
    serialized = json.dumps(report, indent=2, sort_keys=True) + "\n"
    needles = [str(path.resolve()) for path in forbidden_roots]
    if any(needle in serialized for needle in needles):
        raise RuntimeError("peer attestation leaked an A/B checkout root")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized, encoding="utf-8")
    return artifact(output)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=pathlib.Path,
                        default=ROOT / "build/c7-matrix")
    parser.add_argument("--results-dir", type=pathlib.Path,
                        default=HERE / "results")
    parser.add_argument("--manifest", type=pathlib.Path,
                        default=HERE / "results/build-run-manifest.json")
    parser.add_argument("--jobs-per-build", type=int, default=4)
    parser.add_argument("--cpus", type=parse_cpus)
    parser.add_argument("--smoke-cpu", type=int)
    parser.add_argument("--core-git-sha")
    parser.add_argument("--compare-reproducibility-manifest", type=pathlib.Path)
    parser.add_argument("--compare-reproducibility-root", type=pathlib.Path)
    arguments = parser.parse_args()
    allowed = os.sched_getaffinity(0)
    arguments.cpus = arguments.cpus or default_cpus(allowed)
    if arguments.smoke_cpu is None:
        arguments.smoke_cpu = arguments.cpus[0]
    try:
        validate_cpus(arguments.cpus, arguments.smoke_cpu, allowed)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    if arguments.jobs_per_build < 1 or arguments.jobs_per_build > 8:
        raise SystemExit("jobs per build must be in 1..8")
    core_git_sha = arguments.core_git_sha or subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip()
    if not GIT_SHA_RE.fullmatch(core_git_sha):
        raise SystemExit("core git SHA must be 40 lowercase hexadecimal digits")
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip()
    if head != core_git_sha:
        raise SystemExit("--core-git-sha must equal the checked-out HEAD")
    for committed_source in (SOURCE, pathlib.Path(__file__), VALIDATOR):
        if committed_sha256(core_git_sha, committed_source) != sha256(
                committed_source):
            raise SystemExit(
                f"source differs from {core_git_sha}: {relative(committed_source)}")
    arguments.build_dir = arguments.build_dir.resolve()
    arguments.results_dir = arguments.results_dir.resolve()
    arguments.manifest = arguments.manifest.resolve()
    for label, path in (
        ("--build-dir", arguments.build_dir),
        ("--results-dir", arguments.results_dir),
        ("--manifest", arguments.manifest),
    ):
        try:
            path.relative_to(ROOT)
        except ValueError as error:
            raise SystemExit(f"{label} must be inside the source checkout") from error
    if bool(arguments.compare_reproducibility_manifest) != bool(
            arguments.compare_reproducibility_root):
        raise SystemExit(
            "A/B comparison requires both peer manifest and peer source root")
    arguments.build_dir.mkdir(parents=True, exist_ok=True)
    arguments.results_dir.mkdir(parents=True, exist_ok=True)

    gcc = resolve_program("gcc")
    gxx = resolve_program("g++")
    clang = resolve_first(("clang", "clang-18", "clang-17", "clang-16"))
    clangxx = resolve_first(("clang++", "clang++-18", "clang++-17", "clang++-16"))
    build_specs = [(backend, False, gcc, gxx) for backend in BACKENDS]
    build_specs.append(("auto", True, clang, clangxx))
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(
            cmake_build, backend, arguments.build_dir, arguments.results_dir,
            arguments.jobs_per_build, cc, cxx, sanitizer,
        ) for backend, sanitizer, cc, cxx in build_specs]
        core_builds = [future.result() for future in futures]
    for build in core_builds:
        for source in build["source_closure"]:
            require_committed_artifact(core_git_sha, source)
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(
            compile_experiment, build, arguments.build_dir,
            arguments.results_dir, clangxx if build["sanitizer"] else gxx,
            core_git_sha,
        ) for build in core_builds]
        builds = [future.result() for future in futures]
    by_name = {build["name"]: build for build in builds}

    correctness_specs = [
        (by_name[name], arguments.cpus[index], False)
        for index, name in enumerate((*BACKENDS, "asan-ubsan"))
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(
            run_one, build, arguments.results_dir, cpu, smoke,
        ) for build, cpu, smoke in correctness_specs]
        runs = [future.result() for future in futures]
    runs.append(run_one(
        by_name["auto"], arguments.results_dir, arguments.smoke_cpu, True))

    manifest: dict[str, Any] = {
        "schema": "leopard2-c7-build-run-manifest/v3",
        "status": "pass",
        "scope": (
            "correctness plus one affinity-selected non-authoritative harness "
            "smoke; no promotion timing"
        ),
        "normalization": {
            "schema": NORMALIZATION_SCHEMA,
            "token": NORMALIZATION_TOKEN,
            "operation": "replace exact source-root prefix only",
        },
        "core_git_sha": core_git_sha,
        "source": artifact(SOURCE),
        "runner": artifact(pathlib.Path(__file__)),
        "validator": artifact(VALIDATOR),
        "taskset": program_record(resolve_program("taskset")),
        "builds": builds,
        "runs": runs,
        "reproducibility": {
            "prefix_map_target": PREFIX_MAP_TARGET,
            "fingerprints": reproducibility_fingerprints(builds),
            "comparison": {"status": "not-run"},
        },
    }
    forbidden_roots: tuple[pathlib.Path, ...] = (ROOT,)
    if arguments.compare_reproducibility_manifest:
        peer_manifest_path = arguments.compare_reproducibility_manifest.resolve()
        peer_root = arguments.compare_reproducibility_root.resolve()
        if peer_root == ROOT:
            raise RuntimeError("A/B comparison requires a distinct source checkout")
        try:
            peer_manifest_path.relative_to(peer_root)
        except ValueError as error:
            raise RuntimeError("peer manifest is outside its source root") from error
        peer = authenticate_peer(manifest, peer_manifest_path, peer_root)
        forbidden = (ROOT, peer_root)
        current_scan = require_no_root_bytes(manifest, ROOT, forbidden)
        peer_scan = require_no_root_bytes(peer, peer_root, forbidden)
        attestation = write_peer_attestation(
            peer, peer_manifest_path, peer_scan,
            arguments.results_dir / "peer-reproducibility-attestation.json",
            forbidden)
        fingerprints = manifest["reproducibility"]["fingerprints"]
        manifest["reproducibility"]["comparison"] = {
            "status": "pass",
            "peer_manifest_sha256": sha256(peer_manifest_path),
            "fingerprints_sha256": canonical_json_sha256(fingerprints),
            "build_names": list(BUILD_NAMES),
            "checkout_roots_scanned": 2,
            "current_scan": current_scan,
            "peer_scan": peer_scan,
            "peer_attestation": attestation,
        }
        forbidden_roots = forbidden
    serialized = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    if any(str(path.resolve()) in serialized for path in forbidden_roots):
        raise RuntimeError("manifest leaked an absolute A/B source checkout root")
    arguments.manifest.parent.mkdir(parents=True, exist_ok=True)
    arguments.manifest.write_text(serialized, encoding="utf-8")
    if not SHA256_RE.fullmatch(sha256(arguments.manifest)):
        raise AssertionError("manifest hashing failed")
    # Never publish a success exit for evidence that the trusted current
    # validator cannot replay from the just-written bytes.
    trusted_validate_manifest(arguments.manifest, ROOT, live=False)
    trusted_validate_manifest(arguments.manifest, ROOT, live=True)
    return 0


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "--compiler-launch":
        raise SystemExit(compiler_launch(sys.argv[2:]))
    raise SystemExit(main())
