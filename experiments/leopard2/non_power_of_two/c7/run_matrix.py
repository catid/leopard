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
from typing import Any, Iterable


ROOT = pathlib.Path(__file__).resolve().parents[4]
HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "c7_exact_low.cpp"
BACKENDS = ("scalar", "ssse3", "avx2", "auto")
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
GIT_SHA_RE = re.compile(r"[0-9a-f]{40}\Z")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def dependency_closure(build_dir: pathlib.Path) -> list[dict[str, Any]]:
    paths: set[pathlib.Path] = {ROOT / "CMakeLists.txt"}
    for dependency in build_dir.rglob("*.o.d"):
        text = dependency.read_text(encoding="utf-8", errors="strict")
        flattened = text.replace("\\\n", " ")
        if ":" not in flattened:
            raise RuntimeError(f"malformed dependency file: {dependency}")
        for token in shlex.split(flattened.split(":", 1)[1]):
            candidate = pathlib.Path(token)
            if not candidate.is_absolute():
                candidate = (ROOT / candidate).resolve()
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
    configure_argv = [
        str(resolve_program("cmake")), "-S", ".", "-B", relative(build_dir),
        "-G", "Unix Makefiles",
        f"-DCMAKE_BUILD_TYPE={'Debug' if sanitizer else 'Release'}",
        f"-DCMAKE_C_COMPILER={c_compiler}",
        f"-DCMAKE_CXX_COMPILER={cxx_compiler}",
        f"-DLEO2_BACKEND_VARIANT={'auto' if sanitizer else backend}",
        "-DLEO2_BUILD_TESTS=OFF", "-DLEO2_BUILD_BENCHMARKS=OFF",
        "-DLEO2_BUILD_FUZZERS=OFF", "-DLEO2_ENABLE_CUDA=OFF",
        f"-DENABLE_OPENMP={'OFF' if sanitizer else 'ON'}",
    ]
    if sanitizer:
        flags = "-fsanitize=address,undefined -fno-omit-frame-pointer"
        configure_argv.extend([
            f"-DCMAKE_C_FLAGS={flags}", f"-DCMAKE_CXX_FLAGS={flags}",
            f"-DCMAKE_EXE_LINKER_FLAGS={flags}",
        ])
    run_logged(configure_argv, configure_stdout, configure_stderr)
    build_argv = [
        str(resolve_program("cmake")), "--build", relative(build_dir),
        "--target", "libleopard", "--", f"-j{jobs}",
    ]
    run_logged(build_argv, build_stdout, build_stderr)
    library = build_dir / "liblibleopard.a"
    if not library.is_file():
        raise RuntimeError(f"missing library output: {library}")
    return {
        "name": name,
        "backend": "auto" if sanitizer else backend,
        "sanitizer": sanitizer,
        "configure_argv": configure_argv,
        "build_argv": build_argv,
        "configure_stdout": artifact(configure_stdout),
        "configure_stderr": artifact(configure_stderr),
        "build_stdout": artifact(build_stdout),
        "build_stderr": artifact(build_stderr),
        "cmake": program_record(resolve_program("cmake")),
        "c_compiler": program_record(c_compiler),
        "cxx_compiler": program_record(cxx_compiler),
        "library": artifact(library),
        "source_closure": dependency_closure(build_dir),
        "build_dir": relative(build_dir),
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
    argv = [
        str(compiler), "-std=c++11", "-g", "-Wall", "-Wextra",
        "-Wpedantic", "-Werror", "-I.",
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
    nm_output = subprocess.run(
        [str(nm), "-u", str(executable)], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout
    nm_path = results / "logs" / f"{name}-nm-undefined.txt"
    nm_path.write_text(nm_output, encoding="utf-8")
    instrumentation = {
        "required_compile_macro": bool(build["sanitizer"]),
        "undefined_symbol_scan": artifact(nm_path),
        "has_asan_symbols": "__asan_" in nm_output,
        "has_ubsan_symbols": "__ubsan_" in nm_output,
    }
    if build["sanitizer"] and not (
            instrumentation["has_asan_symbols"] and
            instrumentation["has_ubsan_symbols"]):
        raise RuntimeError("sanitizer executable lacks ASan/UBSan references")
    if not build["sanitizer"] and (
            instrumentation["has_asan_symbols"] or
            instrumentation["has_ubsan_symbols"]):
        raise RuntimeError("normal executable unexpectedly references sanitizers")
    result = dict(build)
    result.update({
        "compiler": program_record(compiler),
        "compile_argv": argv,
        "compile_stdout": artifact(stdout_path),
        "compile_stderr": artifact(stderr_path),
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
        relative(executable), "--backend", backend, relative(result_path),
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
    return {
        "name": name,
        "build": build["name"],
        "kind": "non-authoritative-smoke" if smoke else "correctness",
        "requested_cpu": cpu,
        "argv": argv,
        "environment": environment,
        "result": artifact(result_path),
        "stdout": artifact(stdout_path),
        "stderr": artifact(stderr_path),
        "observed_affinity": data["affinity"],
    }


def parse_cpus(text: str) -> list[int]:
    result = [int(item) for item in text.split(",")]
    if len(result) != 5 or len(set(result)) != 5 or any(cpu < 0 for cpu in result):
        raise argparse.ArgumentTypeError("--cpus requires five distinct CPUs")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=pathlib.Path,
                        default=ROOT / "build/c7-matrix")
    parser.add_argument("--results-dir", type=pathlib.Path,
                        default=HERE / "results")
    parser.add_argument("--manifest", type=pathlib.Path,
                        default=HERE / "results/build-run-manifest.json")
    parser.add_argument("--jobs-per-build", type=int, default=4)
    parser.add_argument("--cpus", type=parse_cpus, default=parse_cpus("0,1,2,3,4"))
    parser.add_argument("--smoke-cpu", type=int, default=0)
    parser.add_argument("--core-git-sha")
    arguments = parser.parse_args()
    allowed = os.sched_getaffinity(0)
    if any(cpu not in allowed for cpu in arguments.cpus + [arguments.smoke_cpu]):
        raise SystemExit("requested CPU is outside process affinity")
    if arguments.jobs_per_build < 1 or arguments.jobs_per_build > 8:
        raise SystemExit("jobs per build must be in 1..8")
    if any(cpu in (15, 31) for cpu in arguments.cpus):
        raise SystemExit("CPUs 15 and 31 are reserved from C7 checkpoint work")
    if arguments.smoke_cpu != 0:
        raise SystemExit("the retained non-authoritative smoke must use CPU 0")
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
    for committed_source in (SOURCE, pathlib.Path(__file__)):
        if committed_sha256(core_git_sha, committed_source) != sha256(
                committed_source):
            raise SystemExit(
                f"source differs from {core_git_sha}: {relative(committed_source)}")
    arguments.build_dir = arguments.build_dir.resolve()
    arguments.results_dir = arguments.results_dir.resolve()
    arguments.manifest = arguments.manifest.resolve()
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
        "schema": "leopard2-c7-build-run-manifest/v1",
        "status": "pass",
        "scope": (
            "correctness plus CPU0 non-authoritative harness smoke; no promotion timing"
        ),
        "core_git_sha": core_git_sha,
        "source": artifact(SOURCE),
        "runner": artifact(pathlib.Path(__file__)),
        "taskset": program_record(resolve_program("taskset")),
        "builds": builds,
        "runs": runs,
    }
    arguments.manifest.parent.mkdir(parents=True, exist_ok=True)
    arguments.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not SHA256_RE.fullmatch(sha256(arguments.manifest)):
        raise AssertionError("manifest hashing failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
