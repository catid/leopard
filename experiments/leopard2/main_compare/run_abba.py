#!/usr/bin/env python3
"""Fail-closed Leopard2 versus exact Leopard-main ABBA comparison runner.

The runner deliberately does not build either implementation.  It executes two
already-built, independently linked benchmark processes on one pinned CPU,
retains their byte-for-byte output, and refuses to analyze evidence unless the
workload and build identities match the signed request.
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
import platform
import re
import shlex
import statistics
import subprocess
import sys
import time
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
RAW_SCHEMA = "leopard2-main-compare-raw/v1"
MANIFEST_SCHEMA = "leopard2-main-compare-manifest/v1"
RESERVATION_SCHEMA = "leopard2-cpu-reservation/v1"
ROUNDS = 3
ORDER = ("baseline", "candidate", "candidate", "baseline")
FNV_OFFSET = 14695981039346656037
FNV_PRIME = 1099511628211
MASK64 = (1 << 64) - 1
HEX64 = re.compile(r"^[0-9a-f]{16}$")
HEX256 = re.compile(r"^[0-9a-f]{64}$")
CHILD_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}


class EvidenceError(RuntimeError):
    """The requested run or retained evidence is not authoritative."""


@dataclass(frozen=True)
class Cell:
    identifier: str
    k: int
    r: int
    shard_bytes: int
    losses: int
    seed: int


REPRESENTATIVE_CELLS = (
    Cell("xor-gf8", 129, 1, 65536, 1, 101),
    Cell("gf8-high-one", 240, 16, 65536, 1, 103),
    Cell("gf8-high-full", 240, 16, 65536, 16, 107),
    Cell("gf8-balanced-full", 128, 128, 65536, 128, 109),
    Cell("gf16-inflation-eight", 200, 50, 65536, 8, 113),
    Cell("gf16-high-one", 1000, 200, 65536, 1, 127),
    Cell("gf16-high-full", 1000, 200, 65536, 200, 131),
    Cell("gf16-large-eight", 4096, 512, 4096, 8, 137),
)
SMOKE_CELLS = (Cell("smoke", 8, 4, 64, 1, 1),)


def statistics_policy() -> dict[str, Any]:
    return {
        "method": "one mean log contrast per independent ABBA round",
        "confidence": 0.95,
        "critical_distribution": "Student-t",
        "independent_round_count_per_metric": ROUNDS,
        "constituent_pair_count_per_metric": 2 * ROUNDS,
        "degrees_of_freedom": ROUNDS - 1,
        "child_estimator": "median of retained per-invocation samples",
        "decode_first_use_semantics": (
            "derived median plan-create plus median one-execution with codec already "
            "created; separate timing loops; excludes codec setup"),
        "decode_reuse_amortized_semantics": (
            "derived median execution plus median plan-create divided by reuse; "
            "separate timing loops; excludes codec setup"),
    }


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise EvidenceError(f"value is not canonical JSON: {error}") from error


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def signed(value: Mapping[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(dict(value))
    require("digest" not in result, "digest field already exists")
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def verify_signature(value: object, what: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{what} is not an object")
    payload = copy.deepcopy(value)
    digest = payload.pop("digest", None)
    require(isinstance(digest, str) and HEX256.fullmatch(digest) is not None,
            f"{what} has no canonical SHA-256 digest")
    require(sha256_bytes(canonical_bytes(payload)) == digest,
            f"{what} digest mismatch")
    return value


def write_json_exclusive(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = canonical_bytes(value) + b"\n"
    try:
        with path.open("xb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise EvidenceError(f"refusing to replace evidence file {path}") from error


def write_bytes_exclusive(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("xb") as stream:
            stream.write(value)
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise EvidenceError(f"refusing to replace evidence file {path}") from error


def run_checked(
    arguments: Sequence[str],
    cwd: Path | None = None,
    environment: Mapping[str, str] | None = None,
) -> bytes:
    completed = subprocess.run(
        list(arguments), cwd=None if cwd is None else str(cwd),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        timeout=30, env=None if environment is None else dict(environment))
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace").strip()
        raise EvidenceError(
            f"command failed ({completed.returncode}): {' '.join(arguments)}: {detail}")
    return completed.stdout


def git_identity(root: Path, expected_commit: str, detached: bool) -> dict[str, Any]:
    root = root.resolve(strict=True)
    head = run_checked(("git", "-C", str(root), "rev-parse", "HEAD")).decode().strip()
    require(head == expected_commit,
            f"source {root} is {head}, expected {expected_commit}")
    tree = run_checked(("git", "-C", str(root), "rev-parse", "HEAD^{tree}")) \
        .decode().strip()
    status = run_checked((
        "git", "-C", str(root), "status", "--porcelain=v1",
        "--untracked-files=normal")).decode()
    require(status == "", f"source {root} is not clean: {status!r}")
    symbolic = subprocess.run(
        ("git", "-C", str(root), "symbolic-ref", "-q", "HEAD"),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, timeout=30)
    is_detached = symbolic.returncode != 0
    if detached:
        require(is_detached, f"exact-main source {root} is not detached")
    tracked = run_checked(("git", "-C", str(root), "ls-tree", "-r", "-z", "HEAD"))
    return {
        "path": str(root),
        "head": head,
        "tree": tree,
        "detached": is_detached,
        "tracked_tree_listing_sha256": sha256_bytes(tracked),
        "tracked_status": "clean",
    }


def artifact_identity(path: Path, kind: str) -> dict[str, Any]:
    path = path.resolve(strict=True)
    require(path.is_file(), f"{kind} is not a regular file: {path}")
    stat = path.stat()
    if kind == "executable":
        require(os.access(path, os.X_OK), f"benchmark is not executable: {path}")
    if kind == "archive":
        with path.open("rb") as stream:
            require(stream.read(8) == b"!<arch>\n", f"not an ar archive: {path}")
    return {
        "path": str(path),
        "kind": kind,
        "size": stat.st_size,
        "mode": stat.st_mode & 0o7777,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256_file(path),
    }


def parse_cmake_cache(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        name_and_type, value = line.split("=", 1)
        if ":" not in name_and_type:
            continue
        name, _ = name_and_type.split(":", 1)
        values[name] = value
    return values


def compile_command_tokens(entry: Mapping[str, Any]) -> list[str]:
    arguments = entry.get("arguments")
    if isinstance(arguments, list) and all(isinstance(item, str) for item in arguments):
        return list(arguments)
    command = entry.get("command")
    require(isinstance(command, str), "compile command has neither arguments nor command")
    try:
        return shlex.split(command)
    except ValueError as error:
        raise EvidenceError(f"cannot parse compile command: {error}") from error


def validate_effective_flags(tokens: Sequence[str], what: str) -> None:
    optimizations = [
        token for token in tokens
        if re.fullmatch(r"-O(?:0|1|2|3|g|s|z|fast)", token) is not None
    ]
    require(optimizations and optimizations[-1] == "-O3",
            f"{what} last optimization flag is not -O3: {optimizations}")
    forbidden_prefixes = (
        "-fsanitize", "-fno-sanitize", "-fprofile", "-flto", "-fno-lto",
        "-finstrument-functions", "-fno-tree-vectorize", "-fno-vectorize",
        "-fno-slp-vectorize", "--coverage",
    )
    forbidden_exact = {"-pg", "-coverage"}
    rejected = [
        token for token in tokens
        if token in forbidden_exact or token.startswith(forbidden_prefixes)
    ]
    require(not rejected, f"{what} contains instrumentation/noncanonical flags: {rejected}")


def command_output_path(entry: Mapping[str, Any], tokens: Sequence[str]) -> Path:
    positions = [index for index, token in enumerate(tokens) if token == "-o"]
    require(len(positions) == 1 and positions[0] + 1 < len(tokens),
            "compile command does not have exactly one -o output")
    directory = entry.get("directory")
    require(isinstance(directory, str), "compile command has no working directory")
    output = Path(tokens[positions[0] + 1])
    if not output.is_absolute():
        output = Path(directory) / output
    return output.resolve(strict=True)


def validate_compile_commands(
    path: Path,
    implementation: str,
    specification: Mapping[str, Any],
    compiler: Path,
) -> dict[str, Any]:
    try:
        entries = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise EvidenceError(f"invalid compile_commands.json: {error}") from error
    require(isinstance(entries, list) and entries,
            f"{implementation} compile_commands.json is empty")
    by_source: dict[Path, tuple[list[str], Mapping[str, Any]]] = {}
    for entry in entries:
        require(isinstance(entry, dict) and isinstance(entry.get("file"), str),
                "compile command entry is malformed")
        source = Path(entry["file"]).resolve(strict=True)
        require(source not in by_source, f"duplicate compile command for {source}")
        tokens = compile_command_tokens(entry)
        require(tokens and Path(tokens[0]).resolve(strict=True) == compiler,
                f"compile command for {source} uses the wrong compiler")
        by_source[source] = (tokens, entry)

    baseline_root = Path(specification["baseline_source_root"]).resolve(strict=True)
    candidate_root = Path(specification["candidate_source_root"]).resolve(strict=True)
    if implementation == "baseline":
        required = {
            *(baseline_root / name for name in (
                "leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp")),
            candidate_root / "experiments/leopard2/main_compare/legacy_main_benchmark.cpp",
        }
    else:
        required = {
            *(candidate_root / name for name in (
                "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
                "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp",
                "Leopard2Plan.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp",
                "LeopardFF16.cpp", "Leopard2BackendSSSE3.cpp",
                "Leopard2BackendAVX2.cpp")),
            candidate_root / "bench/leopard2/benchmark.cpp",
        }
    required = {path.resolve(strict=True) for path in required}
    missing = sorted(str(path) for path in required - set(by_source))
    require(not missing, f"{implementation} compile commands miss sources: {missing}")
    object_records: list[dict[str, Any]] = []
    for source in required:
        tokens, entry = by_source[source]
        output = command_output_path(entry, tokens)
        validate_effective_flags(tokens, f"compile command for {source}")
        require("-fopenmp" in tokens, f"{source} was not compiled with OpenMP")
        if implementation == "baseline":
            require("-march=native" in tokens,
                    f"exact-main source {source} lacks -march=native")
        else:
            require(not any(token == "-march=native" or
                            token.startswith("-DLEO2_ENABLE_TEST_HOOKS")
                            for token in tokens),
                    f"candidate production source {source} has native/test-hook flags")
            if source.name not in {"Leopard2BackendAVX2.cpp", "Leopard2BackendSSSE3.cpp"}:
                require("-mavx2" not in tokens and "-mssse3" not in tokens,
                        f"candidate portable source {source} has an ISA-specific flag")
        source_identity = artifact_identity(source, "source_file")
        object_identity = artifact_identity(output, "object_file")
        require(object_identity["mtime_ns"] >= source_identity["mtime_ns"],
                f"object {output} predates source {source}")
        object_records.append({
            "source": source_identity,
            "object": object_identity,
        })
    if implementation == "baseline":
        adapter = (candidate_root /
                   "experiments/leopard2/main_compare/legacy_main_benchmark.cpp").resolve()
        require(any(MAIN_COMMIT in token for token in by_source[adapter][0]),
                "baseline adapter was not compiled with the exact-main commit attestation")
    else:
        avx2 = (candidate_root / "Leopard2BackendAVX2.cpp").resolve()
        ssse3 = (candidate_root / "Leopard2BackendSSSE3.cpp").resolve()
        require("-mavx2" in by_source[avx2][0] and
                "-mno-avx512f" in by_source[avx2][0],
                "candidate AVX2 backend lacks its canonical ISA isolation flags")
        require("-mssse3" in by_source[ssse3][0] and
                "-mno-avx" in by_source[ssse3][0],
                "candidate SSSE3 backend lacks its canonical ISA isolation flags")
    return {
        "entry_count": len(entries),
        "required_sources": sorted(str(path) for path in required),
        "validated_optimization": "-O3",
        "validated_openmp": True,
        "required_source_object_pairs": sorted(
            object_records, key=lambda record: record["source"]["path"]),
        "isa_policy": (
            "whole-build -march=native" if implementation == "baseline" else
            "portable core with ISA flags only on SSSE3/AVX2 translation units"),
    }


def build_provenance(
    implementation: str, specification: Mapping[str, Any]
) -> dict[str, Any]:
    build = Path(specification[f"{implementation}_build_dir"]).resolve(strict=True)
    require(build.is_dir(), f"{implementation} build path is not a directory: {build}")
    names = ({
        "executable": "leopard_main_benchmark",
        "archive": "libleopard_main_exact.a",
        "executable_link": "CMakeFiles/leopard_main_benchmark.dir/link.txt",
        "archive_link": "CMakeFiles/leopard_main_exact.dir/link.txt",
    } if implementation == "baseline" else {
        "executable": "bench_leopard2",
        "archive": "liblibleopard.a",
        "executable_link": "CMakeFiles/bench_leopard2.dir/link.txt",
        "archive_link": "CMakeFiles/libleopard.dir/link.txt",
    })
    expected_executable = (build / names["executable"]).resolve(strict=True)
    expected_archive = (build / names["archive"]).resolve(strict=True)
    require(expected_executable ==
            Path(specification[f"{implementation}_executable"]).resolve(strict=True),
            f"{implementation} executable is not the declared build artifact")
    require(expected_archive ==
            Path(specification[f"{implementation}_archive"]).resolve(strict=True),
            f"{implementation} archive is not the declared build artifact")

    cache_path = build / "CMakeCache.txt"
    commands_path = build / "compile_commands.json"
    executable_link_path = build / names["executable_link"]
    archive_link_path = build / names["archive_link"]
    cache = parse_cmake_cache(cache_path)
    require(cache.get("CMAKE_BUILD_TYPE") == "Release",
            f"{implementation} build is not CMake Release")
    validate_effective_flags(
        shlex.split(cache.get("CMAKE_CXX_FLAGS_RELEASE", "")),
        f"{implementation} CMake Release flags")
    compiler = Path(cache.get("CMAKE_CXX_COMPILER", "")).resolve(strict=True)
    candidate_root = Path(specification["candidate_source_root"]).resolve(strict=True)
    baseline_root = Path(specification["baseline_source_root"]).resolve(strict=True)
    if implementation == "baseline":
        expected_home = candidate_root / "experiments/leopard2/main_compare"
        require(Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve() ==
                expected_home.resolve(), "baseline CMake source is not the adapter directory")
        require(Path(cache.get("LEOPARD_MAIN_SOURCE_DIR", "")).resolve() == baseline_root,
                "baseline CMake cache points at the wrong exact-main source")
        required_cache = {"LEO_MAIN_HAS_MARCH_NATIVE": "1"}
        expected_archive_name = "libleopard_main_exact.a"
    else:
        require(Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve() == candidate_root,
                "candidate CMake source differs from candidate source root")
        required_cache = {
            "ENABLE_OPENMP": "ON",
            "LEO2_BACKEND_VARIANT": "auto",
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_BUILD_TESTS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
        }
        expected_archive_name = "liblibleopard.a"
    for name, expected in required_cache.items():
        require(cache.get(name) == expected,
                f"{implementation} CMake cache {name} is {cache.get(name)!r}, "
                f"expected {expected!r}")
    executable_link = executable_link_path.read_text(encoding="utf-8")
    archive_link = archive_link_path.read_text(encoding="utf-8")
    require(expected_archive_name in executable_link,
            f"{implementation} benchmark link recipe omits its declared archive")
    require(names["archive"] in archive_link,
            f"{implementation} archive recipe does not produce its declared archive")
    executable_link_tokens = shlex.split(executable_link)
    require(executable_link_tokens and
            Path(executable_link_tokens[0]).resolve(strict=True) == compiler,
            f"{implementation} link recipe uses a different compiler")
    validate_effective_flags(
        executable_link_tokens, f"{implementation} executable link recipe")
    semantics = validate_compile_commands(
        commands_path, implementation, specification, compiler)
    records = semantics["required_source_object_pairs"]
    benchmark_suffix = (
        "/experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        if implementation == "baseline" else "/bench/leopard2/benchmark.cpp")
    benchmark_records = [
        record for record in records
        if record["source"]["path"].endswith(benchmark_suffix)
    ]
    require(len(benchmark_records) == 1,
            f"{implementation} has no unique benchmark object")
    archive_records = [record for record in records if record not in benchmark_records]

    def linked_objects(recipe: str) -> list[Path]:
        result: list[Path] = []
        for line in recipe.splitlines():
            for token in shlex.split(line):
                if token.endswith(".o"):
                    item = Path(token)
                    if not item.is_absolute():
                        item = build / item
                    result.append(item.resolve(strict=True))
        return result

    archive_recipe_objects = linked_objects(archive_link)
    executable_recipe_objects = linked_objects(executable_link)
    expected_archive_objects = [
        Path(record["object"]["path"]) for record in archive_records
    ]
    expected_executable_objects = [
        Path(record["object"]["path"]) for record in benchmark_records
    ]
    require(len(archive_recipe_objects) == len(set(archive_recipe_objects)) and
            set(archive_recipe_objects) == set(expected_archive_objects),
            f"{implementation} archive object closure differs from compile commands")
    require(executable_recipe_objects == expected_executable_objects,
            f"{implementation} executable object closure differs from compile commands")
    archive_identity = artifact_identity(expected_archive, "archive")
    executable_identity = artifact_identity(expected_executable, "executable")
    require(all(archive_identity["mtime_ns"] >= record["object"]["mtime_ns"]
                for record in archive_records),
            f"{implementation} archive predates one of its object files")
    require(executable_identity["mtime_ns"] >= archive_identity["mtime_ns"] and
            all(executable_identity["mtime_ns"] >= record["object"]["mtime_ns"]
                for record in benchmark_records),
            f"{implementation} executable predates its link inputs")
    ar = Path("/usr/bin/ar").resolve(strict=True)
    members = run_checked((str(ar), "t", str(expected_archive)),
                          environment=CHILD_ENVIRONMENT).decode().splitlines()
    require(members == [path.name for path in archive_recipe_objects],
            f"{implementation} archive members differ from its link recipe")
    compiler_version = run_checked(
        (str(compiler), "--version"), environment=CHILD_ENVIRONMENT)
    return {
        "build_dir": str(build),
        "cmake_cache": artifact_identity(cache_path, "build_metadata"),
        "compile_commands": artifact_identity(commands_path, "build_metadata"),
        "executable_link_recipe": artifact_identity(
            executable_link_path, "build_metadata"),
        "archive_link_recipe": artifact_identity(archive_link_path, "build_metadata"),
        "compiler": artifact_identity(compiler, "compiler"),
        "compiler_version_stdout": {
            "sha256": sha256_bytes(compiler_version),
            "text": compiler_version.decode("utf-8", errors="strict"),
        },
        "archiver": artifact_identity(ar, "archiver"),
        "validated_archive_members": members,
        "validated_executable": executable_identity,
        "validated_archive": archive_identity,
        "validated_cache": {
            "CMAKE_BUILD_TYPE": cache["CMAKE_BUILD_TYPE"],
            "CMAKE_CXX_COMPILER": cache.get("CMAKE_CXX_COMPILER"),
            "CMAKE_CXX_FLAGS_RELEASE": cache["CMAKE_CXX_FLAGS_RELEASE"],
            **required_cache,
        },
        "validated_compile_commands": semantics,
    }


def runtime_closure(ldd: Path, executable: Path) -> dict[str, Any]:
    ldd = ldd.resolve(strict=True)
    executable = executable.resolve(strict=True)
    output = run_checked(
        (str(ldd), str(executable)), environment=CHILD_ENVIRONMENT
    ).decode("utf-8", errors="strict")
    entries: list[dict[str, Any]] = []
    for raw_line in output.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if "=>" in line:
            soname, target_and_address = (part.strip() for part in line.split("=>", 1))
            require(not target_and_address.startswith("not found"),
                    f"runtime dependency is missing: {line}")
            target = target_and_address.split(" (", 1)[0]
            require(target.startswith("/"), f"unresolved runtime dependency: {line}")
            entries.append({
                "soname": soname,
                "loader_path": target,
                "file": artifact_identity(Path(target), "shared_library"),
            })
            continue
        token = line.split(" (", 1)[0]
        if token.startswith("/"):
            entries.append({
                "soname": Path(token).name,
                "loader_path": token,
                "file": artifact_identity(Path(token), "dynamic_loader"),
            })
        elif token == "linux-vdso.so.1":
            entries.append({"soname": token, "virtual": True})
        else:
            raise EvidenceError(f"unrecognized ldd output: {line}")
    require(entries, f"ldd returned no runtime closure for {executable}")
    require(len({entry["soname"] for entry in entries}) == len(entries),
            f"duplicate runtime dependency in ldd output for {executable}")
    return {
        "executable": str(executable),
        "dependencies": sorted(entries, key=lambda item: item["soname"]),
    }


def input_snapshot(specification: Mapping[str, Any]) -> dict[str, Any]:
    ldd = Path(specification["ldd"])
    baseline_build = build_provenance("baseline", specification)
    candidate_build = build_provenance("candidate", specification)
    require(baseline_build["compiler"] == candidate_build["compiler"] and
            baseline_build["compiler_version_stdout"] ==
            candidate_build["compiler_version_stdout"],
            "baseline and candidate use different compiler binaries or versions")
    return {
        "runner": artifact_identity(Path(specification["runner"]), "file"),
        "taskset": artifact_identity(Path(specification["taskset"]), "executable"),
        "ldd": artifact_identity(ldd, "executable"),
        "baseline_executable": artifact_identity(
            Path(specification["baseline_executable"]), "executable"),
        "candidate_executable": artifact_identity(
            Path(specification["candidate_executable"]), "executable"),
        "baseline_archive": artifact_identity(
            Path(specification["baseline_archive"]), "archive"),
        "candidate_archive": artifact_identity(
            Path(specification["candidate_archive"]), "archive"),
        "baseline_build": baseline_build,
        "candidate_build": candidate_build,
        "baseline_runtime_closure": runtime_closure(
            ldd, Path(specification["baseline_executable"])),
        "candidate_runtime_closure": runtime_closure(
            ldd, Path(specification["candidate_executable"])),
        "baseline_source": git_identity(
            Path(specification["baseline_source_root"]), MAIN_COMMIT, True),
        "candidate_source": git_identity(
            Path(specification["candidate_source_root"]),
            str(specification["candidate_commit"]), False),
    }


def parse_cpu_list(text: str) -> set[int]:
    result: set[int] = set()
    for component in text.strip().split(","):
        if not component:
            continue
        if "-" in component:
            first, last = (int(item) for item in component.split("-", 1))
            require(first <= last, f"invalid CPU range {component!r}")
            result.update(range(first, last + 1))
        else:
            result.add(int(component))
    return result


def validate_topology(cpu: int, sibling: int) -> tuple[set[int], set[int]]:
    require(cpu >= 0 and sibling >= 0 and cpu != sibling,
            "benchmark CPU and reserved sibling must be distinct non-negative CPUs")
    require(hasattr(os, "sched_getaffinity") and hasattr(os, "sched_setaffinity"),
            "Linux scheduling affinity is required")
    allowed = set(os.sched_getaffinity(0))
    require(cpu in allowed and sibling in allowed,
            f"CPU pair {cpu}/{sibling} is not in initial affinity {sorted(allowed)}; "
            "launch the runner with an affinity containing both CPUs")
    topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
    require(topology.is_file(), f"missing topology file {topology}")
    siblings = parse_cpu_list(topology.read_text(encoding="ascii"))
    require(siblings == {cpu, sibling},
            f"physical core must have exactly the reserved pair {cpu}/{sibling}; "
            f"topology reports {sorted(siblings)}")
    housekeeping = allowed - {cpu, sibling}
    require(housekeeping, "no housekeeping CPU remains after reserving the physical core")
    return allowed, housekeeping


def read_text_optional(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        return None


def cpuinfo_identity(cpu: int) -> dict[str, str]:
    blocks = Path("/proc/cpuinfo").read_text(encoding="utf-8").split("\n\n")
    selected: dict[str, str] | None = None
    for block in blocks:
        values: dict[str, str] = {}
        for line in block.splitlines():
            if ":" in line:
                name, value = line.split(":", 1)
                values[name.strip()] = value.strip()
        if values.get("processor") == str(cpu):
            selected = values
            break
    require(selected is not None, f"CPU {cpu} is absent from /proc/cpuinfo")
    retained = {
        name: selected[name] for name in (
            "processor", "vendor_id", "cpu family", "model", "model name",
            "stepping", "microcode", "flags", "Features", "CPU implementer",
            "CPU architecture", "CPU variant", "CPU part", "CPU revision")
        if name in selected
    }
    require(any(name in retained for name in ("model name", "CPU part")),
            f"CPU {cpu} has no retained model identity")
    return retained


def cpu_policy_identity(cpu: int) -> dict[str, Any]:
    root = Path(f"/sys/devices/system/cpu/cpu{cpu}")
    topology_root = root / "topology"
    topology = {
        name: read_text_optional(topology_root / name) for name in (
            "core_id", "physical_package_id", "die_id", "cluster_id",
            "thread_siblings_list", "core_siblings_list")
    }
    require(topology["thread_siblings_list"] is not None,
            f"CPU {cpu} has no thread sibling topology")
    frequency_root = root / "cpufreq"
    frequency = {
        name: read_text_optional(frequency_root / name) for name in (
            "scaling_driver", "scaling_governor", "energy_performance_preference",
            "scaling_min_freq", "scaling_max_freq", "cpuinfo_min_freq",
            "cpuinfo_max_freq")
    }
    return {
        "cpu": cpu,
        "online": read_text_optional(root / "online"),
        "cpuinfo": cpuinfo_identity(cpu),
        "topology": topology,
        "frequency_policy": frequency,
    }


def host_identity(cpu: int, sibling: int, allowed_at_launch: set[int]) -> dict[str, Any]:
    require(cpu in allowed_at_launch and sibling in allowed_at_launch,
            "host identity launch CPU set does not contain the reserved pair")
    turbo_paths = (
        Path("/sys/devices/system/cpu/intel_pstate/no_turbo"),
        Path("/sys/devices/system/cpu/amd_pstate/status"),
        Path("/sys/devices/system/cpu/cpufreq/boost"),
    )
    uname = platform.uname()
    online = read_text_optional(Path("/sys/devices/system/cpu/online"))
    require(online is not None, "cannot read online CPU set")
    return {
        "system": {
            "system": uname.system,
            "node": uname.node,
            "release": uname.release,
            "version": uname.version,
            "machine": uname.machine,
            "python": platform.python_version(),
            "page_size": os.sysconf("SC_PAGE_SIZE"),
        },
        "allowed_cpu_set_at_launch": sorted(allowed_at_launch),
        "online_cpu_set": sorted(parse_cpu_list(online)),
        "benchmark_cpu": cpu_policy_identity(cpu),
        "reserved_sibling": cpu_policy_identity(sibling),
        "turbo_and_pstate": {
            str(path): read_text_optional(path) for path in turbo_paths
        },
    }


class Reservation:
    """Hold the coordinator-created canonical reservation for the whole run."""

    def __init__(self, path: Path, cpu: int, sibling: int):
        self.path = path.resolve(strict=True)
        self.cpu = cpu
        self.sibling = sibling
        self.descriptor: int | None = None
        self.identity: dict[str, Any] | None = None

    def __enter__(self) -> dict[str, Any]:
        flags = os.O_RDONLY
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        self.descriptor = os.open(self.path, flags)
        try:
            fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            os.close(self.descriptor)
            self.descriptor = None
            raise EvidenceError(f"CPU reservation is already locked: {self.path}") from error
        try:
            raw = os.read(self.descriptor, 1024 * 1024)
            payload = parse_reservation(raw, self.cpu, self.sibling)
            self.identity = {
                "path": str(self.path),
                "sha256": sha256_bytes(raw),
                "payload": payload,
                "lock": "exclusive_nonblocking",
            }
            return self.identity
        except Exception:
            fcntl.flock(self.descriptor, fcntl.LOCK_UN)
            os.close(self.descriptor)
            self.descriptor = None
            raise

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        if self.descriptor is not None:
            fcntl.flock(self.descriptor, fcntl.LOCK_UN)
            os.close(self.descriptor)
            self.descriptor = None


def parse_reservation(raw: bytes, cpu: int, sibling: int) -> dict[str, Any]:
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"invalid CPU reservation JSON: {error}") from error
    expected_keys = {
        "benchmark_cpu", "nonce", "owner", "reserved_sibling", "schema", "status",
    }
    require(isinstance(payload, dict) and set(payload) == expected_keys,
            "CPU reservation has unexpected or missing fields")
    require(raw == canonical_bytes(payload),
            "CPU reservation is not canonical JSON without a trailing newline")
    require(payload["schema"] == RESERVATION_SCHEMA and payload["status"] == "held",
            "CPU reservation is not a held v1 reservation")
    require(payload["benchmark_cpu"] == cpu and payload["reserved_sibling"] == sibling,
            "CPU reservation pair does not match the run request")
    require(isinstance(payload["owner"], str) and payload["owner"].strip(),
            "CPU reservation owner is empty")
    require(isinstance(payload["nonce"], str) and 8 <= len(payload["nonce"]) <= 256,
            "CPU reservation nonce is invalid")
    return payload


def validate_reservation_current(identity: Mapping[str, Any]) -> None:
    path = Path(str(identity.get("path", "")))
    require(path.is_file(), "CPU reservation disappeared during the campaign")
    raw = path.read_bytes()
    payload = identity.get("payload")
    require(isinstance(payload, dict), "retained CPU reservation payload is invalid")
    parsed = parse_reservation(
        raw, int(payload["benchmark_cpu"]), int(payload["reserved_sibling"]))
    require(parsed == payload and sha256_bytes(raw) == identity.get("sha256"),
            "CPU reservation changed during the campaign")


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


def expected_losses(cell: Cell) -> list[int]:
    order = list(range(cell.k))
    random = XorShift64(cell.seed ^ 0xD1B54A32D192ED03)
    for remaining in range(len(order), 1, -1):
        selected = random.next() % remaining
        order[remaining - 1], order[selected] = order[selected], order[remaining - 1]
    return sorted(order[:cell.losses])


def ceil_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def validate_cell(cell: Cell) -> None:
    require(re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", cell.identifier) is not None,
            f"invalid cell identifier {cell.identifier!r}")
    require(cell.k > 0 and cell.r > 0 and cell.r <= cell.k,
            f"cell {cell.identifier} is outside exact-main R <= K")
    require(cell.shard_bytes > 0 and cell.shard_bytes % 64 == 0,
            f"cell {cell.identifier} shard bytes are not a positive multiple of 64")
    require(0 <= cell.losses <= cell.r,
            f"cell {cell.identifier} has invalid loss count")
    parent = ceil_power_of_two(cell.k + ceil_power_of_two(cell.r))
    require(parent <= 65536, f"cell {cell.identifier} exceeds the legacy field")
    require(0 <= cell.seed <= MASK64, f"cell {cell.identifier} seed exceeds uint64")


def finite_number(value: object, what: str, positive: bool = True) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{what} is not numeric")
    result = float(value)
    require(math.isfinite(result), f"{what} is not finite")
    if positive:
        require(result > 0, f"{what} is not positive")
    else:
        require(result >= 0, f"{what} is negative")
    return result


def close_enough(actual: float, expected: float) -> bool:
    return abs(actual - expected) <= max(0.000002, abs(expected) * 0.000002)


def median(values: Sequence[float]) -> float:
    return statistics.median(values)


def validate_summary(
    summary: object, iterations: int, setup: bool = False
) -> list[float]:
    require(isinstance(summary, dict), "timing summary is not an object")
    prefix = "" if setup else "_per_batch_call"
    names = {
        "median": f"median_us{prefix}",
        "mad": f"mad_us{prefix}",
        "minimum": f"minimum_us{prefix}",
        "maximum": f"maximum_us{prefix}",
        "samples": "samples_us" if setup else "samples_us_per_batch_call",
    }
    samples_value = summary.get(names["samples"])
    require(isinstance(samples_value, list) and len(samples_value) == iterations,
            f"{names['samples']} does not contain exactly {iterations} samples")
    samples = [finite_number(value, names["samples"]) for value in samples_value]
    derived_median = median(samples)
    deviations = [abs(value - derived_median) for value in samples]
    expected = {
        names["median"]: derived_median,
        names["mad"]: median(deviations),
        names["minimum"]: min(samples),
        names["maximum"]: max(samples),
    }
    for name, derived in expected.items():
        actual = finite_number(summary.get(name), name, positive=name != names["mad"])
        require(close_enough(actual, derived), f"{name} is not derived from raw samples")
    return samples


def validate_digest_object(value: object) -> dict[str, str]:
    require(isinstance(value, dict), "workload_digests is not an object")
    require(value.get("algorithm") == "fnv1a64", "wrong workload digest algorithm")
    result: dict[str, str] = {}
    for name in ("original_data", "transmitted_parity", "recovered_originals"):
        digest = value.get(name)
        require(isinstance(digest, str) and HEX64.fullmatch(digest) is not None,
                f"invalid FNV-1a digest {name}")
        result[name] = digest
    return result


def expected_parameters(cell: Cell, campaign: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "K": cell.k,
        "R": cell.r,
        "shard_bytes": cell.shard_bytes,
        "loss_count": cell.losses,
        "missing_original_indices": expected_losses(cell),
        "batch": 1,
        "reuse": campaign["reuse"],
        "iterations": campaign["iterations"],
        "warmup": campaign["warmup"],
        "thread_count": 1,
        "seed": cell.seed,
    }


def validate_result(
    implementation: str,
    value: object,
    cell: Cell,
    campaign: Mapping[str, Any],
) -> dict[str, Any]:
    require(isinstance(value, dict), "benchmark output is not a JSON object")
    expected_schema = (
        "leopard-main-benchmark-v1" if implementation == "baseline"
        else "leopard2-benchmark-v2")
    require(value.get("schema") == expected_schema,
            f"{implementation} returned wrong schema")
    parameters = value.get("parameters")
    require(isinstance(parameters, dict), "parameters is not an object")
    for name, expected in expected_parameters(cell, campaign).items():
        require(parameters.get(name) == expected,
                f"{implementation} parameter {name} differs: "
                f"{parameters.get(name)!r} != {expected!r}")
    resolved = value.get("resolved")
    require(isinstance(resolved, dict), "resolved identity is not an object")
    padded = ceil_power_of_two(cell.r)
    parent = ceil_power_of_two(cell.k + padded)
    field = "gf8" if parent <= 256 else "gf16"
    for name, expected in {
        "profile": "legacy_high_v1", "field": field,
        "parent_count": parent, "padded_side": padded,
    }.items():
        require(resolved.get(name) == expected,
                f"{implementation} resolved {name} differs")
    if implementation == "baseline":
        require(resolved.get("thread_count") == 1,
                "baseline resolved more than one thread")
        require(value.get("build", {}).get("main_source_commit") == MAIN_COMMIT,
                "baseline did not attest the exact main commit")
        require(value.get("correctness", {}).get("round_trip") is True,
                "baseline round trip failed")
    else:
        required_candidate = {
            "requested_profile": "legacy_high_v1",
            "requested_field": "auto",
            "requested_backend": "auto",
            "force_generic_decode": False,
            "force_specialized_decode": False,
            "skip_legacy": True,
            "retain_samples": True,
        }
        for name, expected in required_candidate.items():
            require(parameters.get(name) == expected,
                    f"candidate option {name} is not comparison-safe")
        require(resolved.get("thread_count") == 1,
                "candidate resolved more than one thread")
        require(resolved.get("backend") in {"scalar", "ssse3", "avx2", "neon"},
                "candidate did not resolve a production backend")
        require(value.get("correctness", {}).get("leopard2_round_trip") is True,
                "candidate round trip failed")
        legacy = value.get("legacy")
        require(isinstance(legacy, dict) and legacy.get("available") is False and
                legacy.get("unavailable_reason") == "disabled by --skip-legacy",
                "candidate silently ran the in-tree legacy comparison")
    digests = validate_digest_object(value.get("workload_digests"))
    metrics = value.get("metrics")
    require(isinstance(metrics, dict), "metrics is not an object")
    iterations = int(campaign["iterations"])
    encode = validate_summary(metrics.get("encode_execution"), iterations)
    if implementation == "baseline":
        require(metrics.get("codec_setup") is None and
                metrics.get("decode_timing_includes_setup") is True,
                "baseline decode setup semantics are ambiguous")
        decode = validate_summary(metrics.get("decode_including_setup"), iterations)
        return {
            "digests": digests,
            "backend": "exact_main_native",
            "encode": encode,
            "decode": decode,
        }
    codec_setup = validate_summary(metrics.get("codec_setup"), iterations, setup=True)
    plan_setup = validate_summary(metrics.get("decode_plan_setup"), iterations, setup=True)
    decode = validate_summary(metrics.get("decode_execution"), iterations)
    amortized = metrics.get("decode_amortized_at_reuse")
    require(isinstance(amortized, dict) and
            amortized.get("reuse_count") == campaign["reuse"],
            "candidate amortized decode has wrong reuse")
    expected_amortized = median(decode) + median(plan_setup) / campaign["reuse"]
    actual_amortized = finite_number(
        amortized.get("derived_median_us_per_batch_call"), "amortized decode median")
    require(close_enough(actual_amortized, expected_amortized),
            "candidate amortized decode is not derived from plan and execution")
    return {
        "digests": digests,
        "backend": resolved["backend"],
        "encode": encode,
        "codec_setup": codec_setup,
        "plan_setup": plan_setup,
        "decode": decode,
    }


def child_time(normalized: Mapping[str, Any], metric: str, reuse: int) -> float:
    if metric == "encode":
        return median(normalized["encode"])
    if "plan_setup" not in normalized:
        return median(normalized["decode"])
    divisor = 1 if metric == "decode_first_use" else reuse
    # Plan creation and execution are measured in separate loops.  Combine their
    # retained medians rather than pretending sample i from each loop was one
    # jointly timed first-use call.  Codec creation is intentionally excluded.
    return median(normalized["decode"]) + median(normalized["plan_setup"]) / divisor


def confidence_interval(log_ratios: Sequence[float]) -> dict[str, Any]:
    # The two contrasts within one ABBA round share drift and are not independent.
    # Each independent observation is therefore the mean log contrast for a whole
    # ABBA round; three rounds use Student-t with df=2.
    require(len(log_ratios) == ROUNDS,
            f"paired analysis requires {ROUNDS} independent round contrasts")
    mean = statistics.fmean(log_ratios)
    deviation = statistics.stdev(log_ratios)
    # Student-t 97.5th percentile with df=2.
    critical = 4.302652729911275
    half_width = critical * deviation / math.sqrt(len(log_ratios))
    lower = math.exp(mean - half_width)
    upper = math.exp(mean + half_width)
    speedup = math.exp(mean)
    return {
        "ratio_orientation": "baseline_time_over_candidate_time",
        "independent_round_count": len(log_ratios),
        "degrees_of_freedom": len(log_ratios) - 1,
        "independent_round_log_contrasts": list(log_ratios),
        "geometric_speedup": speedup,
        "ci95_lower": lower,
        "ci95_upper": upper,
        "faster_ci_excludes_one": lower > 1.0,
        "promotion_lower_bound_at_least_1_05": lower >= 1.05,
        "performance_result_does_not_affect_evidence_validity": True,
    }


def analyze_cell(records: Sequence[Mapping[str, Any]], reuse: int) -> dict[str, Any]:
    require(len(records) == ROUNDS * len(ORDER), "cell has wrong invocation count")
    metrics = ("encode", "decode_first_use", "decode_reuse_amortized")
    observations: dict[str, list[float]] = {name: [] for name in metrics}
    raw_pairs: dict[str, list[list[float]]] = {name: [] for name in metrics}
    for round_index in range(ROUNDS):
        group = records[round_index * 4:(round_index + 1) * 4]
        require(tuple(item["implementation"] for item in group) == ORDER,
                f"round {round_index} is not ABBA")
        round_pairs: dict[str, list[float]] = {name: [] for name in metrics}
        for baseline_index, candidate_index in ((0, 1), (3, 2)):
            baseline = group[baseline_index]["normalized"]
            candidate = group[candidate_index]["normalized"]
            for metric in metrics:
                baseline_time = child_time(baseline, metric, reuse)
                candidate_time = child_time(candidate, metric, reuse)
                round_pairs[metric].append(math.log(baseline_time / candidate_time))
        for metric in metrics:
            require(len(round_pairs[metric]) == 2,
                    f"round {round_index} does not contain two ABBA contrasts")
            raw_pairs[metric].append(round_pairs[metric])
            observations[metric].append(statistics.fmean(round_pairs[metric]))
    result: dict[str, Any] = {}
    for name, values in observations.items():
        result[name] = confidence_interval(values)
        result[name]["within_round_log_contrasts"] = raw_pairs[name]
        result[name]["constituent_pair_count"] = 2 * ROUNDS
    return result


def analyze(invocations: Sequence[Mapping[str, Any]], campaign: Mapping[str, Any]) -> dict[str, Any]:
    by_cell: dict[str, list[Mapping[str, Any]]] = {
        cell["identifier"]: [] for cell in campaign["cells"]
    }
    for invocation in invocations:
        identifier = invocation.get("cell_id")
        require(identifier in by_cell, f"unknown invocation cell {identifier!r}")
        by_cell[identifier].append(invocation)
    return {
        identifier: analyze_cell(records, int(campaign["reuse"]))
        for identifier, records in by_cell.items()
    }


def safe_evidence_path(root: Path, relative: object) -> Path:
    require(isinstance(relative, str) and relative and not os.path.isabs(relative),
            "evidence path is not relative")
    path = (root / relative).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as error:
        raise EvidenceError(f"evidence path escapes output directory: {relative}") from error
    return path


def validate_host_record(
    value: object, cpu: int, sibling: int, allowed: Sequence[int]
) -> None:
    require(isinstance(value, dict), "host identity is not an object")
    require(value.get("allowed_cpu_set_at_launch") == list(allowed),
            "host identity has the wrong launch affinity")
    require(isinstance(value.get("online_cpu_set"), list) and
            cpu in value["online_cpu_set"] and sibling in value["online_cpu_set"],
            "reserved CPUs were not retained as online")
    for name, expected_cpu, expected_sibling in (
        ("benchmark_cpu", cpu, sibling), ("reserved_sibling", sibling, cpu)):
        record = value.get(name)
        require(isinstance(record, dict) and record.get("cpu") == expected_cpu,
                f"host {name} identity is invalid")
        cpuinfo = record.get("cpuinfo")
        topology = record.get("topology")
        policy = record.get("frequency_policy")
        require(isinstance(cpuinfo, dict) and
                any(key in cpuinfo for key in ("model name", "CPU part")),
                f"host {name} model identity is missing")
        require(isinstance(topology, dict) and
                isinstance(topology.get("thread_siblings_list"), str) and
                parse_cpu_list(topology["thread_siblings_list"]) ==
                {expected_cpu, expected_sibling},
                f"host {name} topology is not exactly the reserved SMT pair")
        require(isinstance(policy, dict) and
                all(key in policy for key in (
                    "scaling_driver", "scaling_governor",
                    "energy_performance_preference")),
                f"host {name} frequency policy is missing")
    require(isinstance(value.get("turbo_and_pstate"), dict),
            "host turbo/pstate identity is missing")
    require(isinstance(value.get("system"), dict) and
            isinstance(value["system"].get("release"), str),
            "host system identity is missing")


def validate_raw(
    raw: object,
    output: Path | None,
    check_files: bool,
    check_current_inputs: bool,
) -> dict[str, Any]:
    raw = verify_signature(raw, "raw bundle")
    require(raw.get("schema") == RAW_SCHEMA, "wrong raw bundle schema")
    campaign = raw.get("campaign")
    require(isinstance(campaign, dict), "campaign is not an object")
    require(campaign.get("rounds") == ROUNDS and
            tuple(campaign.get("order", ())) == ORDER,
            "campaign does not contain exactly three ABBA rounds")
    require(campaign.get("batch") == 1 and campaign.get("threads") == 1,
            "campaign is not a one-stripe, one-thread comparison")
    require(campaign.get("child_environment") == CHILD_ENVIRONMENT,
            "campaign child environment is not the strict comparison environment")
    cpu = campaign.get("benchmark_cpu")
    sibling = campaign.get("reserved_sibling")
    require(isinstance(cpu, int) and isinstance(sibling, int) and
            cpu >= 0 and sibling >= 0 and cpu != sibling,
            "campaign has no valid reserved CPU pair")
    require(isinstance(campaign.get("iterations"), int) and
            campaign["iterations"] >= 3,
            "campaign has too few timing samples")
    require(isinstance(campaign.get("reuse"), int) and campaign["reuse"] >= 1,
            "campaign reuse is invalid")
    require(isinstance(campaign.get("warmup"), int) and campaign["warmup"] >= 1,
            "campaign warmup is invalid")
    timeout = campaign.get("timeout_seconds")
    require(isinstance(timeout, (int, float)) and not isinstance(timeout, bool) and
            math.isfinite(timeout) and timeout > 0,
            "campaign timeout is invalid")
    require(campaign.get("statistics") == statistics_policy(),
            "campaign statistics policy is not the authoritative clustered ABBA policy")
    allowed = campaign.get("allowed_cpu_set_at_launch")
    require(isinstance(allowed, list) and allowed == sorted(set(allowed)) and
            all(isinstance(item, int) and item >= 0 for item in allowed) and
            cpu in allowed and sibling in allowed and len(set(allowed) - {cpu, sibling}) > 0,
            "campaign launch affinity is invalid")
    host_initial = raw.get("host_initial")
    host_final = raw.get("host_final")
    require(host_initial == host_final, "host policy/topology changed during campaign")
    validate_host_record(host_initial, cpu, sibling, allowed)
    cells_value = campaign.get("cells")
    require(isinstance(cells_value, list) and cells_value, "campaign has no cells")
    cells = [Cell(**value) for value in cells_value]
    require(len({cell.identifier for cell in cells}) == len(cells),
            "campaign cell identifiers are not unique")
    for cell in cells:
        validate_cell(cell)
    input_spec = raw.get("input_specification")
    initial = raw.get("identities_initial")
    final = raw.get("identities_final")
    require(isinstance(input_spec, dict) and isinstance(initial, dict) and
            isinstance(final, dict), "input identity is missing")
    expected_specification_keys = {
        "runner", "taskset", "ldd", "baseline_executable", "candidate_executable",
        "baseline_archive", "candidate_archive", "baseline_build_dir",
        "candidate_build_dir", "baseline_source_root", "candidate_source_root",
        "candidate_commit",
    }
    require(set(input_spec) == expected_specification_keys and
            all(isinstance(value, str) and value for value in input_spec.values()),
            "input specification is incomplete or has unexpected fields")
    require(re.fullmatch(r"[0-9a-f]{40}", input_spec["candidate_commit"]) is not None,
            "candidate commit is not a full lowercase SHA-1")
    require(initial == final, "input identities changed during the campaign")
    reservation = raw.get("reservation")
    require(isinstance(reservation, dict) and
            reservation.get("lock") == "exclusive_nonblocking" and
            isinstance(reservation.get("path"), str) and
            isinstance(reservation.get("sha256"), str) and
            HEX256.fullmatch(reservation["sha256"]) is not None,
            "retained CPU reservation identity is invalid")
    reservation_payload = reservation.get("payload")
    require(isinstance(reservation_payload, dict),
            "retained CPU reservation payload is missing")
    require(parse_reservation(canonical_bytes(reservation_payload), cpu, sibling) ==
            reservation_payload,
            "retained CPU reservation semantics differ from the campaign")
    require(reservation["sha256"] ==
            sha256_bytes(canonical_bytes(reservation_payload)),
            "retained CPU reservation hash does not match its canonical payload")
    if check_current_inputs:
        require(input_snapshot(input_spec) == initial,
                "current executable/archive/source identity differs from retained evidence")
        validate_reservation_current(reservation)
    invocations = raw.get("invocations")
    expected_count = len(cells) * ROUNDS * len(ORDER)
    require(isinstance(invocations, list) and len(invocations) == expected_count,
            f"campaign has {0 if not isinstance(invocations, list) else len(invocations)} "
            f"invocations, expected {expected_count}")
    digest_by_cell: dict[str, dict[str, str]] = {}
    candidate_backend_by_cell: dict[str, str] = {}
    cell_by_id = {cell.identifier: cell for cell in cells}
    expected_sequence = [
        (cell.identifier, round_index, slot, implementation)
        for cell in cells for round_index in range(ROUNDS)
        for slot, implementation in enumerate(ORDER)
    ]
    for invocation, expected in zip(invocations, expected_sequence):
        require(isinstance(invocation, dict), "invocation is not an object")
        observed = (
            invocation.get("cell_id"), invocation.get("round"),
            invocation.get("slot"), invocation.get("implementation"))
        require(observed == expected,
                f"invocation order/relabel mismatch: {observed!r} != {expected!r}")
        require(invocation.get("identity_before") == initial and
                invocation.get("identity_after") == initial,
                "an invocation observed a changed input identity")
        require(invocation.get("reservation_before") == reservation and
                invocation.get("reservation_after") == reservation,
                "an invocation observed a changed CPU reservation")
        require(invocation.get("returncode") == 0,
                "benchmark child did not exit successfully")
        cell = cell_by_id[expected[0]]
        expected_command = [
            input_spec["taskset"], "-c", str(cpu),
            *benchmark_arguments(expected[3], Path(input_spec[
                f"{expected[3]}_executable"]), cell, campaign),
        ]
        require(invocation.get("command") == expected_command and
                invocation.get("pinned_cpu") == cpu,
                "benchmark command or CPU pinning was edited")
        require(invocation.get("environment") == CHILD_ENVIRONMENT,
                "benchmark invocation inherited or retained an unsafe environment")
        result = invocation.get("result")
        if check_files:
            require(output is not None, "output root is required for file verification")
            for stream_name in ("stdout", "stderr"):
                stream = invocation.get(stream_name)
                require(isinstance(stream, dict), f"missing {stream_name} evidence")
                path = safe_evidence_path(output, stream.get("path"))
                require(path.is_file(), f"missing retained {stream_name}: {path}")
                data = path.read_bytes()
                require(stream.get("size") == len(data) and
                        stream.get("sha256") == sha256_bytes(data),
                        f"retained {stream_name} identity mismatch")
                if stream_name == "stdout":
                    try:
                        parsed = json.loads(data.decode("utf-8"))
                    except (UnicodeDecodeError, json.JSONDecodeError) as error:
                        raise EvidenceError(f"retained stdout is not JSON: {error}") from error
                    require(parsed == result, "parsed retained stdout differs from raw result")
        normalized = validate_result(
            expected[3], result, cell, campaign)
        require(invocation.get("normalized") == normalized,
                "retained normalized benchmark data was edited")
        digests = normalized["digests"]
        if expected[0] in digest_by_cell:
            require(digests == digest_by_cell[expected[0]],
                    f"workload FNV digests differ in cell {expected[0]}")
        else:
            digest_by_cell[expected[0]] = digests
        if expected[3] == "candidate":
            backend = normalized["backend"]
            if expected[0] in candidate_backend_by_cell:
                require(backend == candidate_backend_by_cell[expected[0]],
                        f"candidate backend changed within cell {expected[0]}")
            else:
                candidate_backend_by_cell[expected[0]] = backend
    calculated = analyze(invocations, campaign)
    require(raw.get("analysis") == calculated, "paired-log analysis was edited")
    return calculated


def benchmark_arguments(
    implementation: str, executable: Path, cell: Cell, campaign: Mapping[str, Any]
) -> list[str]:
    arguments = [
        str(executable), "--k", str(cell.k), "--r", str(cell.r),
        "--bytes", str(cell.shard_bytes), "--loss", str(cell.losses),
        "--batch", "1", "--reuse", str(campaign["reuse"]),
        "--iterations", str(campaign["iterations"]),
        "--warmup", str(campaign["warmup"]), "--threads", "1",
        "--seed", str(cell.seed),
    ]
    if implementation == "candidate":
        arguments.extend((
            "--profile", "high", "--field", "auto", "--backend", "auto",
            "--skip-legacy", "--retain-samples",
        ))
    arguments.extend(("--json", "-"))
    return arguments


def run_child(
    implementation: str,
    cell: Cell,
    round_index: int,
    slot: int,
    campaign: Mapping[str, Any],
    specification: Mapping[str, Any],
    initial_identity: Mapping[str, Any],
    reservation: Mapping[str, Any],
    output: Path,
    cpu: int,
    timeout: float,
) -> dict[str, Any]:
    executable = Path(specification[f"{implementation}_executable"])
    child_arguments = benchmark_arguments(implementation, executable, cell, campaign)
    command = [specification["taskset"], "-c", str(cpu), *child_arguments]
    before = input_snapshot(specification)
    require(before == initial_identity, "input identity changed before benchmark launch")
    validate_reservation_current(reservation)
    environment = dict(CHILD_ENVIRONMENT)
    start_utc = utc_now()
    start = time.monotonic_ns()
    try:
        completed = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=environment, check=False, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        stdout = error.stdout or b""
        stderr = error.stderr or b""
        completed = subprocess.CompletedProcess(command, 124, stdout, stderr)
    duration_ns = time.monotonic_ns() - start
    stem = f"invocations/{cell.identifier}/round-{round_index}/slot-{slot}-{implementation}"
    stdout_path = output / f"{stem}.stdout"
    stderr_path = output / f"{stem}.stderr"
    write_bytes_exclusive(stdout_path, completed.stdout)
    write_bytes_exclusive(stderr_path, completed.stderr)
    after = input_snapshot(specification)
    require(after == initial_identity, "input identity changed after benchmark launch")
    validate_reservation_current(reservation)
    require(completed.returncode == 0,
            f"{implementation} exited {completed.returncode}; see {stderr_path}")
    try:
        result = json.loads(completed.stdout.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"{implementation} stdout is not one JSON value: {error}") from error
    normalized = validate_result(implementation, result, cell, campaign)
    return {
        "cell_id": cell.identifier,
        "round": round_index,
        "slot": slot,
        "implementation": implementation,
        "command": command,
        "environment": environment,
        "pinned_cpu": cpu,
        "started_utc": start_utc,
        "duration_ns": duration_ns,
        "returncode": completed.returncode,
        "stdout": {
            "path": str(stdout_path.relative_to(output)),
            "size": len(completed.stdout),
            "sha256": sha256_bytes(completed.stdout),
        },
        "stderr": {
            "path": str(stderr_path.relative_to(output)),
            "size": len(completed.stderr),
            "sha256": sha256_bytes(completed.stderr),
        },
        "result": result,
        "normalized": normalized,
        "identity_before": before,
        "identity_after": after,
        "reservation_before": reservation,
        "reservation_after": reservation,
    }


def parse_cell(text: str) -> Cell:
    components = text.split(":")
    require(len(components) == 6,
            "--cell must be ID:K:R:BYTES:LOSSES:SEED")
    try:
        cell = Cell(
            components[0], *(int(component, 10) for component in components[1:]))
    except ValueError as error:
        raise EvidenceError(f"invalid --cell {text!r}: {error}") from error
    validate_cell(cell)
    return cell


def cells_from_options(options: argparse.Namespace) -> tuple[Cell, ...]:
    if options.cell:
        cells = tuple(parse_cell(value) for value in options.cell)
    else:
        cells = REPRESENTATIVE_CELLS if options.preset == "representative" else SMOKE_CELLS
    require(len({cell.identifier for cell in cells}) == len(cells),
            "cell identifiers must be unique")
    return cells


def run_campaign(options: argparse.Namespace) -> int:
    output = options.output.resolve()
    require(not output.exists(), f"output path already exists: {output}")
    output.mkdir(parents=True)
    cells = cells_from_options(options)
    for cell in cells:
        validate_cell(cell)
    taskset = Path(options.taskset).resolve(strict=True)
    ldd = Path(options.ldd).resolve(strict=True)
    require(taskset == Path("/usr/bin/taskset").resolve(strict=True),
            "authoritative comparison requires /usr/bin/taskset")
    require(ldd == Path("/usr/bin/ldd").resolve(strict=True),
            "authoritative comparison requires /usr/bin/ldd")
    specification = {
        "runner": str(Path(__file__).resolve()),
        "taskset": str(taskset),
        "ldd": str(ldd),
        "baseline_executable": str(options.baseline.resolve(strict=True)),
        "candidate_executable": str(options.candidate.resolve(strict=True)),
        "baseline_archive": str(options.baseline_archive.resolve(strict=True)),
        "candidate_archive": str(options.candidate_archive.resolve(strict=True)),
        "baseline_build_dir": str(options.baseline_build_dir.resolve(strict=True)),
        "candidate_build_dir": str(options.candidate_build_dir.resolve(strict=True)),
        "baseline_source_root": str(options.baseline_source_root.resolve(strict=True)),
        "candidate_source_root": str(options.candidate_source_root.resolve(strict=True)),
        "candidate_commit": options.candidate_commit,
    }
    campaign = {
        "rounds": ROUNDS,
        "order": list(ORDER),
        "cells": [asdict(cell) for cell in cells],
        "batch": 1,
        "reuse": options.reuse,
        "iterations": options.iterations,
        "warmup": options.warmup,
        "threads": 1,
        "child_environment": dict(CHILD_ENVIRONMENT),
        "benchmark_cpu": options.cpu,
        "reserved_sibling": options.reserved_sibling,
        "timeout_seconds": options.timeout,
        "statistics": statistics_policy(),
    }
    require(options.iterations >= 3, "--iterations must be at least 3")
    require(options.reuse >= 1 and options.warmup >= 1,
            "--reuse and --warmup must be positive")
    require(math.isfinite(options.timeout) and options.timeout > 0,
            "--timeout must be positive and finite")
    try:
        allowed_at_launch, housekeeping = validate_topology(
            options.cpu, options.reserved_sibling)
        campaign["allowed_cpu_set_at_launch"] = sorted(allowed_at_launch)
        host_initial = host_identity(
            options.cpu, options.reserved_sibling, allowed_at_launch)
        with Reservation(
            options.reservation_file, options.cpu, options.reserved_sibling
        ) as reservation:
            os.sched_setaffinity(0, housekeeping)
            initial = input_snapshot(specification)
            invocations: list[dict[str, Any]] = []
            for cell in cells:
                for round_index in range(ROUNDS):
                    for slot, implementation in enumerate(ORDER):
                        invocations.append(run_child(
                            implementation, cell, round_index, slot, campaign,
                            specification, initial, reservation, output,
                            options.cpu, options.timeout))
            final = input_snapshot(specification)
            require(final == initial, "input identity changed during campaign")
            host_final = host_identity(
                options.cpu, options.reserved_sibling, allowed_at_launch)
            require(host_final == host_initial,
                    "host topology/frequency policy changed during campaign")
            analysis = analyze(invocations, campaign)
            raw = signed({
                "schema": RAW_SCHEMA,
                "created_utc": utc_now(),
                "validity_is_independent_of_speed": True,
                "campaign": campaign,
                "host_initial": host_initial,
                "reservation": reservation,
                "input_specification": specification,
                "identities_initial": initial,
                "invocations": invocations,
                "identities_final": final,
                "host_final": host_final,
                "analysis": analysis,
            })
            validate_raw(raw, output, check_files=True, check_current_inputs=True)
            raw_path = output / "raw.json"
            write_json_exclusive(raw_path, raw)
            manifest = signed({
                "schema": MANIFEST_SCHEMA,
                "created_utc": utc_now(),
                "valid": True,
                "validity_is_independent_of_speed": True,
                "raw": {
                    "path": "raw.json",
                    "size": raw_path.stat().st_size,
                    "sha256": sha256_file(raw_path),
                    "payload_digest": raw["digest"],
                },
                "campaign": campaign,
                "host": host_initial,
                "reservation": reservation,
                "identities": initial,
                "analysis": analysis,
            })
            write_json_exclusive(output / "manifest.json", manifest)
    except Exception as error:
        failure = signed({
            "schema": "leopard2-main-compare-failure/v1",
            "created_utc": utc_now(),
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(),
        })
        failure_path = output / "failure.json"
        if not failure_path.exists():
            write_json_exclusive(failure_path, failure)
        raise
    print(output / "manifest.json")
    return 0


def verify_campaign(options: argparse.Namespace) -> int:
    manifest_path = options.manifest.resolve(strict=True)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    verify_signature(manifest, "manifest")
    require(manifest.get("schema") == MANIFEST_SCHEMA and manifest.get("valid") is True,
            "manifest is not valid main-comparison evidence")
    output = manifest_path.parent
    raw_info = manifest.get("raw")
    require(isinstance(raw_info, dict), "manifest has no raw bundle identity")
    raw_path = safe_evidence_path(output, raw_info.get("path"))
    require(raw_path.is_file(), "retained raw bundle is missing")
    require(raw_info.get("size") == raw_path.stat().st_size and
            raw_info.get("sha256") == sha256_file(raw_path),
            "raw bundle file identity mismatch")
    raw = json.loads(raw_path.read_text(encoding="utf-8"))
    analysis = validate_raw(
        raw, output, check_files=True,
        check_current_inputs=not options.no_current_input_check)
    require(raw_info.get("payload_digest") == raw.get("digest"),
            "manifest/raw payload identity mismatch")
    for name in ("campaign", "host", "reservation", "identities", "analysis"):
        if name == "identities":
            expected = raw["identities_initial"]
        elif name == "host":
            expected = raw["host_initial"]
        else:
            expected = raw[name]
        require(manifest.get(name) == expected,
                f"manifest {name} differs from retained raw bundle")
    require(manifest.get("analysis") == analysis, "manifest analysis was edited")
    if options.no_current_input_check:
        print("exact-main ABBA bundle structurally verified only; current build/source "
              "closure was not revalidated")
    else:
        print("exact-main ABBA evidence and current build/source closure verified")
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="run a new, isolated three-round ABBA campaign")
    run.add_argument("--baseline", required=True, type=Path)
    run.add_argument("--candidate", required=True, type=Path)
    run.add_argument("--baseline-archive", required=True, type=Path)
    run.add_argument("--candidate-archive", required=True, type=Path)
    run.add_argument("--baseline-build-dir", required=True, type=Path)
    run.add_argument("--candidate-build-dir", required=True, type=Path)
    run.add_argument("--baseline-source-root", required=True, type=Path)
    run.add_argument("--candidate-source-root", required=True, type=Path)
    run.add_argument("--candidate-commit", required=True)
    run.add_argument("--reservation-file", required=True, type=Path)
    run.add_argument("--output", required=True, type=Path)
    run.add_argument("--cpu", required=True, type=int)
    run.add_argument("--reserved-sibling", required=True, type=int)
    run.add_argument("--taskset", default="/usr/bin/taskset", type=Path)
    run.add_argument("--ldd", default="/usr/bin/ldd", type=Path)
    run.add_argument("--preset", choices=("representative", "smoke"),
                     default="representative")
    run.add_argument("--cell", action="append",
                     help="override preset with ID:K:R:BYTES:LOSSES:SEED")
    run.add_argument("--reuse", type=int, default=8)
    run.add_argument("--iterations", type=int, default=9)
    run.add_argument("--warmup", type=int, default=2)
    run.add_argument("--timeout", type=float, default=120.0)
    run.set_defaults(function=run_campaign)
    verify = commands.add_parser("verify", help="verify retained evidence and raw outputs")
    verify.add_argument("--manifest", required=True, type=Path)
    verify.add_argument("--no-current-input-check", action="store_true",
                        help="structural-only replay without revalidating original build paths")
    verify.set_defaults(function=verify_campaign)
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    options = parser().parse_args(arguments)
    try:
        return int(options.function(options))
    except EvidenceError as error:
        print(f"main comparison evidence error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
