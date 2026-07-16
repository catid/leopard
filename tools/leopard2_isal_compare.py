#!/usr/bin/env python3
"""Build, run, and validate the default-off Intel ISA-L comparison.

The runner deliberately keeps third-party code under the ignored research
cache.  It verifies exact source/tool identities, executes single-thread cells
in an ABBA provider order on one requested CPU, retains every child result, and
derives the compact summary from those results.  ISA-L and Leopard2 parity are
not wire-compatible; only identical public payload and erasure workloads are
compared.
"""

from __future__ import annotations

import argparse
import copy
import fcntl
import hashlib
import json
import math
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-isal-comparison-checkpoint/v2"
CORRECTNESS_SCHEMA = "leopard2-isal-correctness/v2"
ISA_SCHEMA = "leopard2-isal-benchmark/v2"
LEOPARD_SCHEMA = "leopard2-benchmark-v1"
ISA_URL = "https://github.com/intel/isa-l"
ISA_COMMIT = "e8cc5e87fc64b4da434f32bc1fa18184622a4998"
ISA_TREE = "e56f9556f55549c39d90e2abfca2961c6426702e"
ISA_VERSION = "2.32.1"
ISA_LICENSE_SHA256 = "bc8fd4a3d031e65e05e9c9e2add2c3f336ce527fa85c1e31031c808b58216217"
ISA_HEADER_SHA256 = "976e0cd7f3fc09ac385a75559451bca01209a12abe1f4af214acb3fc94ea0fd9"
NASM_VERSION = "2.16.03"
NASM_URL = (
    "https://www.nasm.us/pub/nasm/releasebuilds/2.16.03/"
    "nasm-2.16.03.tar.xz")
NASM_ARCHIVE_SHA256 = (
    "1412a1c760bbd05db026b6c0d1657affd6631cd0a63cddb6f73cc6d4aa616148")
ABBA = (("isa-l", "leopard2"), ("leopard2", "isa-l"),
        ("leopard2", "isa-l"), ("isa-l", "leopard2"))
SOURCE_BUNDLE_SCHEMA = "leopard2-isal-build-input-closure/v2"
TOOLCHAIN_SCHEMA = "leopard2-isal-build-toolchain/v2"
LEASE_SCHEMA = "leopard2-advisory-cpu-lease/v1"
GENERATED_EVIDENCE_PATHS = (
    "experiments/leopard2/isal_compare/checkpoint_result.json",
    "experiments/leopard2/isal_compare/correctness_result.json",
)
REQUIRED_BUILD_INPUTS = {
    "CMakeLists.txt", ".gitmodules", "sse2neon",
    "leopard2.cpp", "leopard2.h", "Leopard2Backend.cpp", "Leopard2Backend.h",
    "Leopard2BackendScalar.cpp", "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp", "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp", "Leopard2Plan.h", "LeopardCommon.cpp", "LeopardCommon.h",
    "LeopardFF8.cpp", "LeopardFF8.h", "LeopardFF16.cpp", "LeopardFF16.h",
    "bench/leopard2/benchmark.cpp",
    "experiments/leopard2/isal_compare/CMakeLists.txt",
    "experiments/leopard2/isal_compare/isal_benchmark.cpp",
    "experiments/leopard2/isal_compare/NOTICE",
    "tools/leopard2_isal_compare.py", "tools/leopard2_external_comparison.py",
}

# Bounded checkpoint only.  It covers low/balanced/high rates, tiny-batch
# packing, and two cases where public K+R still fits GF8 but Leopard's V1
# dyadic parent selects GF16.  The required 7,134-cell matrix remains open.
CHECKPOINT_CELLS = (
    {"K": 240, "R": 16, "profile": "high", "shard_bytes": 65536,
     "loss_count": 1, "batch": 1, "reuse": 8, "seed": 1101},
    {"K": 240, "R": 16, "profile": "high", "shard_bytes": 4096,
     "loss_count": 4, "batch": 8, "reuse": 8, "seed": 1102},
    {"K": 128, "R": 128, "profile": "high", "shard_bytes": 65536,
     "loss_count": 8, "batch": 1, "reuse": 8, "seed": 1201},
    {"K": 64, "R": 192, "profile": "low", "shard_bytes": 65536,
     "loss_count": 8, "batch": 1, "reuse": 8, "seed": 1301},
    {"K": 129, "R": 100, "profile": "high", "shard_bytes": 65536,
     "loss_count": 4, "batch": 1, "reuse": 8, "seed": 1401},
    {"K": 225, "R": 30, "profile": "high", "shard_bytes": 65536,
     "loss_count": 2, "batch": 1, "reuse": 8, "seed": 1402},
)


class ComparisonError(ValueError):
    """Invalid configuration, child result, or checkpoint evidence."""


def require_exact_keys(value: Mapping[str, Any], expected: set[str], label: str) -> None:
    if not isinstance(value, Mapping) or set(value) != expected:
        raise ComparisonError(
            f"{label} keys changed: got {sorted(value) if isinstance(value, Mapping) else value!r}, "
            f"expected {sorted(expected)}")


def require_hex(value: Any, length: int, label: str) -> str:
    if (not isinstance(value, str) or len(value) != length or
            any(character not in "0123456789abcdef" for character in value)):
        raise ComparisonError(f"{label} is not {length}-character lowercase hex")
    return value


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def default_cache() -> Path:
    return repo_root() / ".research" / "leopard2"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_digest(document: Mapping[str, Any]) -> str:
    payload = copy.deepcopy(dict(document))
    payload.pop("artifact_sha256", None)
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def keyed_digest(document: Mapping[str, Any], digest_key: str) -> str:
    payload = copy.deepcopy(dict(document))
    payload.pop(digest_key, None)
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def run_checked(command: Sequence[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        list(command), text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False, **kwargs)
    if completed.returncode != 0:
        raise ComparisonError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}")
    return completed


def git_identity(checkout: Path) -> tuple[str, str, str]:
    head = run_checked(["git", "-C", str(checkout), "rev-parse", "HEAD"]).stdout.strip()
    tree = run_checked(
        ["git", "-C", str(checkout), "rev-parse", "HEAD^{tree}"]).stdout.strip()
    dirty = full_git_status(checkout)
    return head, tree, dirty


def full_git_status(checkout: Path) -> str:
    return run_checked([
        "git", "-C", str(checkout), "status", "--porcelain",
        "--untracked-files=normal"]).stdout.strip()


def leopard_build_status(checkout: Path) -> str:
    command = [
        "git", "-C", str(checkout), "status", "--porcelain",
        "--untracked-files=normal", "--", ".",
    ]
    command.extend(f":(exclude){path}" for path in GENERATED_EVIDENCE_PATHS)
    return run_checked(command).stdout.strip()


def _is_build_input(relative: str) -> bool:
    path = Path(relative)
    if relative in ("CMakeLists.txt", ".gitmodules", "sse2neon"):
        return True
    if len(path.parts) == 1 and path.suffix in (".c", ".cc", ".cpp", ".h", ".hpp", ".in"):
        return True
    if relative.startswith("cmake/"):
        return True
    if relative.startswith("bench/leopard2/") and path.suffix in (".c", ".cc", ".cpp", ".h", ".hpp"):
        return True
    if relative.startswith("experiments/leopard2/isal_compare/"):
        return path.name not in ("checkpoint_result.json", "correctness_result.json")
    return relative in (
        "tools/leopard2_external_comparison.py",
        "tools/leopard2_isal_compare.py",
    )


def _committed_build_inputs(
        root: Path, commit: str) -> dict[str, tuple[str, str, str]]:
    raw = run_checked([
        "git", "-C", str(root), "ls-tree", "-r", "-z", "--full-tree", commit,
    ]).stdout
    selected: dict[str, tuple[str, str, str]] = {}
    for record in raw.split("\0"):
        if not record:
            continue
        metadata, relative = record.split("\t", 1)
        mode, object_type, object_id = metadata.split()
        if _is_build_input(relative):
            selected[relative] = (mode, object_type, object_id)
    if not REQUIRED_BUILD_INPUTS.issubset(selected):
        raise ComparisonError("build-input closure selection is incomplete")
    return selected


def _git_object_bytes(root: Path, object_id: str) -> bytes:
    completed = subprocess.run(
        ["git", "-C", str(root), "cat-file", "blob", object_id],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if completed.returncode != 0:
        raise ComparisonError(
            f"cannot read committed Git blob {object_id}: "
            f"{completed.stderr.decode('utf-8', errors='replace')}")
    return completed.stdout


def source_bundle(root: Path, commit: str = "HEAD") -> dict[str, Any]:
    source_commit = run_checked([
        "git", "-C", str(root), "rev-parse", f"{commit}^{{commit}}",
    ]).stdout.strip()
    source_tree = run_checked([
        "git", "-C", str(root), "rev-parse", f"{source_commit}^{{tree}}",
    ]).stdout.strip()
    entries: dict[str, dict[str, str]] = {}
    for relative, (mode, object_type, object_id) in sorted(
            _committed_build_inputs(root, source_commit).items()):
        if mode == "160000":
            if object_type != "commit":
                raise ComparisonError(f"submodule build input is not a commit: {relative}")
            digest = hashlib.sha256(object_id.encode("ascii")).hexdigest()
            entries[relative] = {
                "mode": mode, "git_object": object_id, "sha256": digest,
                "submodule_commit": object_id,
            }
        else:
            if object_type != "blob":
                raise ComparisonError(f"build input is not a blob: {relative}")
            digest = hashlib.sha256(_git_object_bytes(root, object_id)).hexdigest()
            entries[relative] = {
                "mode": mode, "git_object": object_id, "sha256": digest}
    encoded = json.dumps(entries, sort_keys=True, separators=(",", ":"))
    return {
        "schema": SOURCE_BUNDLE_SCHEMA,
        "source_commit": source_commit,
        "source_tree": source_tree,
        "entries": entries,
        "entry_count": len(entries),
        "bundle_sha256": hashlib.sha256(encoded.encode("utf-8")).hexdigest(),
    }


def validate_source_bundle(bundle: Mapping[str, Any], verify_local: bool = False) -> None:
    require_exact_keys(bundle, {
        "schema", "source_commit", "source_tree", "entries", "entry_count",
        "bundle_sha256"},
                       "source bundle")
    if bundle.get("schema") != SOURCE_BUNDLE_SCHEMA:
        raise ComparisonError("wrong ISA-L source-bundle schema")
    require_hex(bundle.get("source_commit"), 40, "source-bundle commit")
    require_hex(bundle.get("source_tree"), 40, "source-bundle tree")
    entries = bundle.get("entries")
    if not isinstance(entries, Mapping) or not entries:
        raise ComparisonError("ISA-L source-bundle entries are missing")
    if not REQUIRED_BUILD_INPUTS.issubset(entries):
        raise ComparisonError("ISA-L source-bundle omitted required production inputs")
    for relative, entry in entries.items():
        if not isinstance(relative, str) or not isinstance(entry, Mapping):
            raise ComparisonError("invalid ISA-L source-bundle entry")
        mode = entry.get("mode")
        expected_keys = ({"mode", "git_object", "sha256", "submodule_commit"}
                         if mode == "160000" else
                         {"mode", "git_object", "sha256"})
        require_exact_keys(entry, expected_keys, f"source bundle entry {relative}")
        require_hex(entry.get("sha256"), 64, f"source bundle entry {relative} digest")
        git_object = require_hex(
            entry.get("git_object"), 40, f"source bundle entry {relative} Git object")
        if mode not in ("100644", "100755", "120000", "160000"):
            raise ComparisonError(f"invalid source-bundle mode for {relative}")
        if mode == "160000":
            submodule = require_hex(entry.get("submodule_commit"), 40,
                                    f"source bundle submodule {relative}")
            if submodule != git_object:
                raise ComparisonError(f"source-bundle submodule identity differs: {relative}")
    if bundle.get("entry_count") != len(entries):
        raise ComparisonError("ISA-L source-bundle entry count is not derived")
    encoded = json.dumps(dict(entries), sort_keys=True, separators=(",", ":"))
    expected = hashlib.sha256(encoded.encode("utf-8")).hexdigest()
    if bundle.get("bundle_sha256") != expected:
        raise ComparisonError("ISA-L source-bundle digest is not file-derived")
    reconstructed = source_bundle(repo_root(), str(bundle.get("source_commit")))
    if reconstructed != dict(bundle):
        raise ComparisonError(
            "retained source bundle does not match its recorded Git commit/tree/blobs")
    if verify_local:
        local_commit = run_checked(
            ["git", "-C", str(repo_root()), "rev-parse", "HEAD"]).stdout.strip()
        local_tree = run_checked(
            ["git", "-C", str(repo_root()), "rev-parse", "HEAD^{tree}"]).stdout.strip()
        local_dirty = leopard_build_status(repo_root())
        if local_dirty:
            raise ComparisonError("local source checkout is not completely clean")
        if (local_commit != bundle.get("source_commit") or
                local_tree != bundle.get("source_tree")):
            raise ComparisonError("local source checkout differs from benchmark-source commit")


def build_paths(cache: Path) -> dict[str, Path]:
    return {
        "isa_source": cache / "isa-l",
        "nasm_archive": cache / "downloads" / f"nasm-{NASM_VERSION}.tar.xz",
        "nasm_source": cache / "toolchains" / f"nasm-{NASM_VERSION}-src",
        "nasm_install": cache / "toolchains" / f"nasm-{NASM_VERSION}-install",
        "isa_build": cache / "build" / f"isa-l-{ISA_COMMIT[:8]}",
        "isa_install": cache / "toolchains" / f"isa-l-{ISA_COMMIT[:8]}-install",
        "adapter_build": cache / "build" / "leopard2-isal-adapter",
        "leopard_build": cache / "build" / "leopard2-isal-leopard",
        "provenance": cache / "toolchains" / "isal-comparison-provenance.json",
    }


def clone_isa_l(source: Path) -> None:
    source.parent.mkdir(parents=True, exist_ok=True)
    if not source.exists():
        run_checked(["git", "clone", "--filter=blob:none", ISA_URL, str(source)],
                    timeout=300)
    run_checked(["git", "-C", str(source), "fetch", "origin", ISA_COMMIT],
                timeout=300)
    run_checked(["git", "-C", str(source), "checkout", "--detach", ISA_COMMIT],
                timeout=60)


def ensure_nasm(paths: Mapping[str, Path], jobs: int) -> dict[str, Any]:
    archive = paths["nasm_archive"]
    archive.parent.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        run_checked([
            "curl", "-fL", "--retry", "3", "--connect-timeout", "20",
            "-o", str(archive), NASM_URL], timeout=300)
    archive_hash = sha256_file(archive)
    if archive_hash != NASM_ARCHIVE_SHA256:
        raise ComparisonError(
            f"NASM archive hash mismatch: {archive_hash} != {NASM_ARCHIVE_SHA256}")

    install = paths["nasm_install"]
    executable = install / "bin" / "nasm"
    version_ok = False
    if executable.exists():
        version = run_checked([str(executable), "-v"]).stdout.strip()
        version_ok = f"NASM version {NASM_VERSION}" in version
    if not version_ok:
        source = paths["nasm_source"]
        if source.exists():
            shutil.rmtree(source)
        if install.exists():
            shutil.rmtree(install)
        source.mkdir(parents=True)
        run_checked([
            "tar", "-xJf", str(archive), "--strip-components=1", "-C", str(source)],
            timeout=120)
        run_checked(["./configure", f"--prefix={install}"], cwd=source, timeout=120)
        run_checked(["make", f"-j{jobs}"], cwd=source, timeout=600)
        run_checked(["make", "install"], cwd=source, timeout=120)
    version = run_checked([str(executable), "-v"]).stdout.strip()
    if f"NASM version {NASM_VERSION}" not in version:
        raise ComparisonError(f"unexpected locally built NASM version: {version}")
    return {
        "version": NASM_VERSION,
        "url": NASM_URL,
        "archive_sha256": archive_hash,
        "executable_sha256": sha256_file(executable),
        "reported_version": version,
        "path": str(executable),
    }


def configure_build_install_isa(
        paths: Mapping[str, Path], nasm: Mapping[str, Any], jobs: int) -> dict[str, Any]:
    source = paths["isa_source"]
    if not source.exists():
        clone_isa_l(source)
    head, tree, dirty = git_identity(source)
    if head != ISA_COMMIT or tree != ISA_TREE or dirty:
        raise ComparisonError(
            f"ISA-L checkout must be completely clean at {ISA_COMMIT}/{ISA_TREE}; "
            f"head={head}, tree={tree}, dirty={dirty!r}")
    build = paths["isa_build"]
    install = paths["isa_install"]
    run_checked([
        "cmake", "-S", str(source), "-B", str(build),
        "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_INSTALL_PREFIX={install}",
        f"-DCMAKE_ASM_NASM_COMPILER={nasm['path']}",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DBUILD_SHARED_LIBS=OFF", "-DISAL_BUILD_TESTS=OFF",
        "-DISAL_BUILD_PERF_TESTS=OFF", "-DISAL_BUILD_ISAL_SHIM=OFF",
        "-DISAL_BUILD_IGZIP_CLI=OFF"], timeout=120)
    run_checked(["cmake", "--build", str(build), "-j", str(jobs)], timeout=900)
    run_checked(["cmake", "--install", str(build)], timeout=180)
    library = install / "lib" / "libisal.a"
    header = install / "include" / "isa-l" / "erasure_code.h"
    source_header = source / "include" / "erasure_code.h"
    license_path = source / "LICENSE"
    if (not library.is_file() or not header.is_file() or
            not source_header.is_file() or not license_path.is_file()):
        raise ComparisonError("private ISA-L installation is incomplete")
    header_hash = sha256_file(header)
    if header_hash != ISA_HEADER_SHA256 or sha256_file(source_header) != header_hash:
        raise ComparisonError("installed ISA-L erasure-code header identity mismatch")
    return {
        "name": "Intel ISA-L",
        "version": ISA_VERSION,
        "url": ISA_URL,
        "commit": head,
        "tree": tree,
        "license": "BSD-3-Clause",
        "license_sha256": sha256_file(license_path),
        "header_sha256": header_hash,
        "library_sha256": sha256_file(library),
        "library": str(library),
        "install": str(install),
    }


def _cmake_cache_value(build: Path, name: str) -> str:
    prefix = name + ":"
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if line.startswith(prefix) and "=" in line:
            return line.split("=", 1)[1]
    raise ComparisonError(f"CMake cache variable is missing: {name} in {build}")


def _tool_identity(path: Path, version_arguments: Sequence[str]) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise ComparisonError(f"tool executable is missing: {path}")
    completed = run_checked([str(resolved), *version_arguments], timeout=30)
    reported = (completed.stdout + completed.stderr).strip()
    if not reported:
        raise ComparisonError(f"tool emitted no version identity: {path}")
    return {"sha256": sha256_file(resolved), "reported_version": reported}


def _replace_paths(value: str, replacements: Mapping[str, str]) -> str:
    result = value
    for source, token in sorted(replacements.items(), key=lambda item: -len(item[0])):
        result = result.replace(source, token)
    return result


def _compile_manifest(path: Path, replacements: Mapping[str, str]) -> dict[str, Any]:
    try:
        raw = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ComparisonError(f"invalid compile-command database {path}: {error}") from error
    if not isinstance(raw, list) or not raw:
        raise ComparisonError(f"compile-command database is empty: {path}")
    entries = []
    for item in raw:
        if not isinstance(item, Mapping) or "file" not in item or "directory" not in item:
            raise ComparisonError(f"malformed compile-command entry in {path}")
        command = item.get("command")
        if command is None and isinstance(item.get("arguments"), list):
            command = json.dumps(item["arguments"], separators=(",", ":"))
        if not isinstance(command, str):
            raise ComparisonError(f"compile-command text is missing in {path}")
        entry = {
            "directory": _replace_paths(str(item["directory"]), replacements),
            "file": _replace_paths(str(item["file"]), replacements),
            "command": _replace_paths(command, replacements),
        }
        if "output" in item:
            entry["output"] = _replace_paths(str(item["output"]), replacements)
        entries.append(entry)
    entries.sort(key=lambda entry: (
        entry["file"], entry.get("output", ""), entry["command"]))
    encoded = json.dumps(entries, sort_keys=True, separators=(",", ":"))
    return {
        "schema": "leopard2-normalized-compile-commands/v1",
        "entries": entries,
        "entry_count": len(entries),
        "sha256": hashlib.sha256(encoded.encode("utf-8")).hexdigest(),
    }


def build_identity(paths: Mapping[str, Path], nasm: Mapping[str, Any]) -> dict[str, Any]:
    root = repo_root().resolve()
    c_compiler = Path(_cmake_cache_value(paths["isa_build"], "CMAKE_C_COMPILER"))
    adapter_cxx = Path(_cmake_cache_value(paths["adapter_build"], "CMAKE_CXX_COMPILER"))
    leopard_cxx = Path(_cmake_cache_value(paths["leopard_build"], "CMAKE_CXX_COMPILER"))
    if adapter_cxx.resolve() != leopard_cxx.resolve():
        raise ComparisonError("adapter and Leopard benchmarks used different C++ compilers")
    ar = Path(_cmake_cache_value(paths["isa_build"], "CMAKE_AR"))
    ranlib = Path(_cmake_cache_value(paths["isa_build"], "CMAKE_RANLIB"))
    linker = Path(_cmake_cache_value(paths["isa_build"], "CMAKE_LINKER"))
    build_program = Path(_cmake_cache_value(paths["isa_build"], "CMAKE_MAKE_PROGRAM"))
    for build in (paths["adapter_build"], paths["leopard_build"]):
        for cache_name, expected in (
                ("CMAKE_AR", ar), ("CMAKE_RANLIB", ranlib),
                ("CMAKE_LINKER", linker), ("CMAKE_MAKE_PROGRAM", build_program)):
            observed = Path(_cmake_cache_value(build, cache_name))
            if observed.resolve() != expected.resolve():
                raise ComparisonError(
                    f"comparison builds used different {cache_name} tools")
    replacements = {
        str(root): "${LEOPARD_ROOT}",
        str(paths["isa_source"].resolve()): "${ISA_L_SOURCE}",
        str(paths["isa_install"].resolve()): "${ISA_L_INSTALL}",
        str(paths["isa_build"].resolve()): "${ISA_L_BUILD}",
        str(paths["adapter_build"].resolve()): "${ADAPTER_BUILD}",
        str(paths["leopard_build"].resolve()): "${LEOPARD_BUILD}",
        str(c_compiler.resolve()): "${CC}",
        str(adapter_cxx.resolve()): "${CXX}",
        str(Path(str(nasm["path"])).resolve()): "${NASM}",
        str(ar.resolve()): "${AR}",
        str(ranlib.resolve()): "${RANLIB}",
        str(linker.resolve()): "${LINKER}",
        str(build_program.resolve()): "${BUILD_PROGRAM}",
    }
    recipe = {
        "isa_l": {
            "configure_definitions": [
                "CMAKE_BUILD_TYPE=Release", "BUILD_SHARED_LIBS=OFF",
                "ISAL_BUILD_TESTS=OFF", "ISAL_BUILD_PERF_TESTS=OFF",
                "ISAL_BUILD_ISAL_SHIM=OFF", "ISAL_BUILD_IGZIP_CLI=OFF",
                "CMAKE_EXPORT_COMPILE_COMMANDS=ON",
            ],
            "generator": _cmake_cache_value(paths["isa_build"], "CMAKE_GENERATOR"),
            "build_target": "all", "install_target": "install",
        },
        "adapter": {
            "configure_definitions": [
                "CMAKE_BUILD_TYPE=Release", "CMAKE_EXPORT_COMPILE_COMMANDS=ON",
                "LEO2_ISAL_EXPECTED_COMMIT=" + ISA_COMMIT,
                "LEO2_ISAL_EXPECTED_LICENSE_SHA256=" + ISA_LICENSE_SHA256,
                "LEO2_ISAL_EXPECTED_HEADER_SHA256=" + ISA_HEADER_SHA256,
                "LEO2_NASM_EXECUTABLE_SHA256=" + str(nasm["executable_sha256"]),
                "LEO2_NASM_ARCHIVE_SHA256=" + NASM_ARCHIVE_SHA256,
            ],
            "generator": _cmake_cache_value(paths["adapter_build"], "CMAKE_GENERATOR"),
            "build_target": "leopard2_isal_benchmark",
        },
        "leopard2": {
            "configure_definitions": [
                "CMAKE_BUILD_TYPE=Release", "LEO2_BUILD_TESTS=OFF",
                "LEO2_BUILD_BENCHMARKS=ON", "LEO2_BUILD_FUZZERS=OFF",
                "LEO2_ENABLE_CUDA=OFF", "ENABLE_OPENMP=ON",
                "CMAKE_EXPORT_COMPILE_COMMANDS=ON",
            ],
            "generator": _cmake_cache_value(paths["leopard_build"], "CMAKE_GENERATOR"),
            "build_target": "bench_leopard2",
        },
    }
    identity = {
        "schema": TOOLCHAIN_SCHEMA,
        "recipe": recipe,
        "tools": {
            "cmake": _tool_identity(Path(shutil.which("cmake") or ""), ["--version"]),
            "cc": _tool_identity(c_compiler, ["--version"]),
            "cxx": _tool_identity(adapter_cxx, ["--version"]),
            "ar": _tool_identity(ar, ["--version"]),
            "ranlib": _tool_identity(ranlib, ["--version"]),
            "linker": _tool_identity(linker, ["--version"]),
            "build_program": _tool_identity(build_program, ["--version"]),
            "nasm": {
                "sha256": str(nasm["executable_sha256"]),
                "reported_version": str(nasm["reported_version"]),
            },
        },
        "compile_commands": {
            "isa_l": _compile_manifest(
                paths["isa_build"] / "compile_commands.json", replacements),
            "adapter": _compile_manifest(
                paths["adapter_build"] / "compile_commands.json", replacements),
            "leopard2": _compile_manifest(
                paths["leopard_build"] / "compile_commands.json", replacements),
        },
    }
    identity["identity_sha256"] = keyed_digest(identity, "identity_sha256")
    return identity


def build_benchmarks(
        paths: Mapping[str, Path], jobs: int,
        nasm: Mapping[str, Any]) -> dict[str, Any]:
    root = repo_root()
    head, tree, tracked_dirty = git_identity(root)
    status = full_git_status(root)
    if tracked_dirty or status:
        raise ComparisonError(
            "benchmark executables must be built from a completely clean commit; "
            f"HEAD={head}, status={status!r}")
    bundle = source_bundle(root)
    isa_install = paths["isa_install"]
    adapter_build = paths["adapter_build"]
    run_checked([
        "cmake", "-S", str(root / "experiments/leopard2/isal_compare"),
        "-B", str(adapter_build), "-DCMAKE_BUILD_TYPE=Release",
        f"-DLEO2_ISAL_SOURCE_ROOT={paths['isa_source']}",
        f"-DLEO2_ISAL_INSTALL_ROOT={isa_install}",
        f"-DLEO2_ISAL_EXPECTED_HEADER_SHA256={ISA_HEADER_SHA256}",
        f"-DLEO2_NASM_EXECUTABLE_SHA256={nasm['executable_sha256']}",
        f"-DLEO2_NASM_ARCHIVE_SHA256={nasm['archive_sha256']}",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"], timeout=120)
    run_checked(["cmake", "--build", str(adapter_build), "-j", str(jobs)],
                timeout=300)

    leopard_build = paths["leopard_build"]
    run_checked([
        "cmake", "-S", str(root), "-B", str(leopard_build),
        "-DCMAKE_BUILD_TYPE=Release", "-DLEO2_BUILD_TESTS=OFF",
        "-DLEO2_BUILD_BENCHMARKS=ON", "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_ENABLE_CUDA=OFF", "-DENABLE_OPENMP=ON",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"], timeout=120)
    run_checked([
        "cmake", "--build", str(leopard_build), "-j", str(jobs),
        "--target", "bench_leopard2"], timeout=600)

    isa_exe = adapter_build / "leopard2_isal_benchmark"
    leopard_exe = leopard_build / "bench_leopard2"
    if not isa_exe.is_file() or not leopard_exe.is_file():
        raise ComparisonError("comparison executables were not built")
    return {
        "isa-l": {"path": str(isa_exe), "sha256": sha256_file(isa_exe)},
        "leopard2": {"path": str(leopard_exe), "sha256": sha256_file(leopard_exe)},
        "leopard_source": {
            "commit": head, "tree": tree, "dirty": False, "clean_at_build": True,
        },
        "source_bundle": bundle,
        "build_identity": build_identity(paths, nasm),
    }


def _validate_compile_manifest(value: Mapping[str, Any], label: str) -> None:
    require_exact_keys(value, {"schema", "entries", "entry_count", "sha256"}, label)
    if value.get("schema") != "leopard2-normalized-compile-commands/v1":
        raise ComparisonError(f"wrong {label} schema")
    entries = value.get("entries")
    if not isinstance(entries, list) or not entries:
        raise ComparisonError(f"{label} entries are missing")
    for entry in entries:
        if not isinstance(entry, Mapping) or set(entry) not in (
                {"directory", "file", "command"},
                {"directory", "file", "command", "output"}):
            raise ComparisonError(f"{label} entry shape changed")
        if any(not isinstance(item, str) or not item for item in entry.values()):
            raise ComparisonError(f"{label} entry contains an invalid value")
        if not (entry["file"].startswith("${LEOPARD_ROOT}/") or
                entry["file"].startswith("${ISA_L_SOURCE}/")):
            raise ComparisonError(f"{label} source path is not portable")
        if any(prefix in text for text in entry.values()
               for prefix in ("/home/", "/Users/")):
            raise ComparisonError(f"{label} retained a checkout-specific path")
    if value.get("entry_count") != len(entries):
        raise ComparisonError(f"{label} count is not entry-derived")
    encoded = json.dumps(entries, sort_keys=True, separators=(",", ":"))
    if value.get("sha256") != hashlib.sha256(encoded.encode("utf-8")).hexdigest():
        raise ComparisonError(f"{label} digest is not entry-derived")


def validate_build_identity(
        value: Mapping[str, Any], nasm: Mapping[str, Any]) -> None:
    require_exact_keys(value, {
        "schema", "recipe", "tools", "compile_commands", "identity_sha256"},
                       "build identity")
    if value.get("schema") != TOOLCHAIN_SCHEMA:
        raise ComparisonError("wrong build-identity schema")
    if value.get("identity_sha256") != keyed_digest(value, "identity_sha256"):
        raise ComparisonError("build identity digest is not content-derived")
    recipe = value.get("recipe", {})
    require_exact_keys(recipe, {"isa_l", "adapter", "leopard2"}, "build recipe")
    for name, expected_target in (
            ("isa_l", "all"), ("adapter", "leopard2_isal_benchmark"),
            ("leopard2", "bench_leopard2")):
        record = recipe.get(name, {})
        expected_keys = ({"configure_definitions", "generator", "build_target", "install_target"}
                         if name == "isa_l" else
                         {"configure_definitions", "generator", "build_target"})
        require_exact_keys(record, expected_keys, f"{name} build recipe")
        definitions = record.get("configure_definitions")
        if (not isinstance(definitions, list) or not definitions or
                any(not isinstance(item, str) or not item for item in definitions)):
            raise ComparisonError(f"{name} configure recipe is invalid")
        if not isinstance(record.get("generator"), str) or not record["generator"]:
            raise ComparisonError(f"{name} CMake generator is missing")
        if record.get("build_target") != expected_target:
            raise ComparisonError(f"{name} build target changed")
    expected_definitions = {
        "isa_l": [
            "CMAKE_BUILD_TYPE=Release", "BUILD_SHARED_LIBS=OFF",
            "ISAL_BUILD_TESTS=OFF", "ISAL_BUILD_PERF_TESTS=OFF",
            "ISAL_BUILD_ISAL_SHIM=OFF", "ISAL_BUILD_IGZIP_CLI=OFF",
            "CMAKE_EXPORT_COMPILE_COMMANDS=ON"],
        "adapter": [
            "CMAKE_BUILD_TYPE=Release", "CMAKE_EXPORT_COMPILE_COMMANDS=ON",
            "LEO2_ISAL_EXPECTED_COMMIT=" + ISA_COMMIT,
            "LEO2_ISAL_EXPECTED_LICENSE_SHA256=" + ISA_LICENSE_SHA256,
            "LEO2_ISAL_EXPECTED_HEADER_SHA256=" + ISA_HEADER_SHA256,
            "LEO2_NASM_EXECUTABLE_SHA256=" + str(nasm.get("executable_sha256")),
            "LEO2_NASM_ARCHIVE_SHA256=" + NASM_ARCHIVE_SHA256],
        "leopard2": [
            "CMAKE_BUILD_TYPE=Release", "LEO2_BUILD_TESTS=OFF",
            "LEO2_BUILD_BENCHMARKS=ON", "LEO2_BUILD_FUZZERS=OFF",
            "LEO2_ENABLE_CUDA=OFF", "ENABLE_OPENMP=ON",
            "CMAKE_EXPORT_COMPILE_COMMANDS=ON"],
    }
    for name, expected in expected_definitions.items():
        if recipe[name]["configure_definitions"] != expected:
            raise ComparisonError(f"{name} build definitions changed")
    if recipe["isa_l"].get("install_target") != "install":
        raise ComparisonError("ISA-L install target changed")
    tools = value.get("tools", {})
    require_exact_keys(tools, {
        "cmake", "cc", "cxx", "nasm", "ar", "ranlib", "linker",
        "build_program"}, "build tools")
    for name, tool in tools.items():
        require_exact_keys(tool, {"sha256", "reported_version"}, f"{name} tool")
        require_hex(tool.get("sha256"), 64, f"{name} tool digest")
        if not isinstance(tool.get("reported_version"), str) or not tool["reported_version"]:
            raise ComparisonError(f"{name} tool version is missing")
    if (tools["nasm"]["sha256"] != nasm.get("executable_sha256") or
            tools["nasm"]["reported_version"] != nasm.get("reported_version")):
        raise ComparisonError("NASM tool identity is not provenance-derived")
    commands = value.get("compile_commands", {})
    require_exact_keys(commands, {"isa_l", "adapter", "leopard2"},
                       "compile-command manifests")
    for name, manifest in commands.items():
        _validate_compile_manifest(manifest, f"{name} compile commands")


def bootstrap(cache: Path, jobs: int) -> dict[str, Any]:
    if jobs < 1 or jobs > 8:
        raise ComparisonError("--jobs must be in 1..8")
    paths = build_paths(cache)
    cache.mkdir(parents=True, exist_ok=True)
    nasm = ensure_nasm(paths, jobs)
    for name in ("isa_build", "isa_install", "adapter_build", "leopard_build"):
        if paths[name].exists():
            shutil.rmtree(paths[name])
    isa = configure_build_install_isa(paths, nasm, jobs)
    executables = build_benchmarks(paths, jobs, nasm)
    provenance = {
        "schema": "leopard2-isal-toolchain/v2",
        "nasm": nasm,
        "isa_l": isa,
        "executables": executables,
        "constraints": {
            "third_party_location": "ignored .research cache",
            "production_dependency": False,
            "maximum_build_jobs": 8,
            "clean_rebuild": True,
        },
    }
    paths["provenance"].parent.mkdir(parents=True, exist_ok=True)
    paths["provenance"].write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    return provenance


def load_provenance(cache: Path) -> dict[str, Any]:
    path = build_paths(cache)["provenance"]
    try:
        document = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ComparisonError(f"run bootstrap first: {error}") from error
    require_exact_keys(document, {"schema", "nasm", "isa_l", "executables", "constraints"},
                       "toolchain provenance")
    if document.get("schema") != "leopard2-isal-toolchain/v2":
        raise ComparisonError("invalid toolchain provenance schema")
    constraints = document.get("constraints", {})
    require_exact_keys(constraints, {
        "third_party_location", "production_dependency", "maximum_build_jobs",
        "clean_rebuild"},
        "toolchain constraints")
    if constraints != {
            "third_party_location": "ignored .research cache",
            "production_dependency": False, "maximum_build_jobs": 8,
            "clean_rebuild": True}:
        raise ComparisonError("toolchain constraints changed")
    nasm = document.get("nasm", {})
    require_exact_keys(nasm, {"version", "url", "archive_sha256", "executable_sha256",
                              "reported_version", "path"}, "NASM provenance")
    if (nasm.get("version") != NASM_VERSION or nasm.get("url") != NASM_URL or
            nasm.get("archive_sha256") != NASM_ARCHIVE_SHA256 or
            f"NASM version {NASM_VERSION}" not in str(nasm.get("reported_version", ""))):
        raise ComparisonError("toolchain provenance has the wrong NASM identity")
    require_hex(nasm.get("executable_sha256"), 64, "NASM executable digest")
    nasm_path = Path(str(nasm.get("path", "")))
    if (not nasm_path.is_file() or
            sha256_file(nasm_path) != nasm.get("executable_sha256")):
        raise ComparisonError("NASM executable missing or changed after bootstrap")
    isa = document.get("isa_l", {})
    require_exact_keys(isa, {"name", "version", "url", "commit", "tree", "license",
                             "license_sha256", "header_sha256", "library_sha256",
                             "library", "install"},
                       "ISA-L provenance")
    if (isa.get("name") != "Intel ISA-L" or isa.get("version") != ISA_VERSION or
            isa.get("url") != ISA_URL or isa.get("commit") != ISA_COMMIT or
            isa.get("tree") != ISA_TREE or isa.get("license") != "BSD-3-Clause" or
            isa.get("license_sha256") != ISA_LICENSE_SHA256 or
            isa.get("header_sha256") != ISA_HEADER_SHA256):
        raise ComparisonError("toolchain provenance has the wrong ISA-L identity")
    require_hex(isa.get("library_sha256"), 64, "ISA-L library digest")
    isa_library = Path(str(isa.get("library", "")))
    if not isa_library.is_file() or sha256_file(isa_library) != isa.get("library_sha256"):
        raise ComparisonError("ISA-L library missing or changed after bootstrap")
    isa_header = Path(str(isa.get("install", ""))) / "include" / "isa-l" / "erasure_code.h"
    if not isa_header.is_file() or sha256_file(isa_header) != isa.get("header_sha256"):
        raise ComparisonError("ISA-L installed header missing or changed after bootstrap")
    isa_head, isa_tree, isa_dirty = git_identity(build_paths(cache)["isa_source"])
    if isa_head != ISA_COMMIT or isa_tree != ISA_TREE or isa_dirty:
        raise ComparisonError("ISA-L source checkout changed after bootstrap")
    executables = document.get("executables", {})
    require_exact_keys(executables, {
        "isa-l", "leopard2", "leopard_source", "source_bundle", "build_identity"},
        "benchmark executable provenance")
    leopard_source = executables.get("leopard_source", {})
    require_exact_keys(leopard_source, {"commit", "tree", "dirty", "clean_at_build"},
                       "Leopard source provenance")
    require_hex(leopard_source.get("commit"), 40, "Leopard source commit")
    require_hex(leopard_source.get("tree"), 40, "Leopard source tree")
    current_head = run_checked(
        ["git", "-C", str(repo_root()), "rev-parse", "HEAD"]).stdout.strip()
    current_tree = run_checked(
        ["git", "-C", str(repo_root()), "rev-parse", "HEAD^{tree}"]).stdout.strip()
    current_dirty = leopard_build_status(repo_root())
    if (leopard_source.get("dirty") is not False or
            leopard_source.get("clean_at_build") is not True or current_dirty or
            current_head != leopard_source.get("commit") or
            current_tree != leopard_source.get("tree")):
        raise ComparisonError("current benchmark source checkout is not clean")
    recorded_bundle = document.get("executables", {}).get("source_bundle", {})
    validate_source_bundle(recorded_bundle)
    if source_bundle(repo_root(), current_head) != recorded_bundle:
        raise ComparisonError("current build-input closure changed after bootstrap")
    validate_build_identity(executables.get("build_identity", {}), nasm)
    if build_identity(build_paths(cache), nasm) != executables.get("build_identity"):
        raise ComparisonError("toolchain or normalized compile commands changed after bootstrap")
    for provider in ("isa-l", "leopard2"):
        record = document.get("executables", {}).get(provider, {})
        require_exact_keys(record, {"path", "sha256"}, f"{provider} executable provenance")
        require_hex(record.get("sha256"), 64, f"{provider} executable digest")
        path_value = Path(str(record.get("path", "")))
        if not path_value.is_file() or sha256_file(path_value) != record.get("sha256"):
            raise ComparisonError(f"{provider} executable missing or changed after bootstrap")
    return document


def ceil_pow2(value: int) -> int:
    return 1 << (value - 1).bit_length()


def resolved_identity(cell: Mapping[str, Any]) -> tuple[str, str, int]:
    k, r = int(cell["K"]), int(cell["R"])
    requested = str(cell["profile"])
    profile = ("low_v1" if requested in ("low", "low_v1") or
               requested == "auto" and r > k else "legacy_high_v1")
    parent = (ceil_pow2(ceil_pow2(k) + r) if profile == "low_v1"
              else ceil_pow2(k + ceil_pow2(r)))
    return profile, "gf8" if parent <= 256 else "gf16", parent


def finite_positive(value: Any, label: str, allow_zero: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ComparisonError(f"{label} is not numeric")
    converted = float(value)
    if not math.isfinite(converted) or (converted < 0 if allow_zero else converted <= 0):
        raise ComparisonError(f"{label} is not finite and {'nonnegative' if allow_zero else 'positive'}")
    return converted


def close_enough(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=2e-6, abs_tol=2e-6)


def median(values: Sequence[float]) -> float:
    return float(statistics.median(values))


def mad(values: Sequence[float]) -> float:
    center = median(values)
    return median([abs(value - center) for value in values])


def expected_parameters(cell: Mapping[str, Any], iterations: int, warmup: int) -> dict[str, Any]:
    profile, _, _ = resolved_identity(cell)
    return {
        "K": int(cell["K"]), "R": int(cell["R"]),
        "requested_profile": profile, "requested_field": "auto",
        "requested_backend": "auto", "shard_bytes": int(cell["shard_bytes"]),
        "loss_count": int(cell["loss_count"]), "batch": int(cell["batch"]),
        "reuse": int(cell["reuse"]), "iterations": iterations,
        "warmup": warmup, "thread_count": 1, "seed": int(cell["seed"]),
    }


def validate_summary(
        summary: Mapping[str, Any], prefix: str, iterations: int,
        execution: bool, input_bytes: int = 0, output_bytes: int = 0,
        input_rate: str = "", output_rate: str = "", require_samples: bool = False) -> None:
    suffix = "_per_batch_call" if execution else ""
    expected_keys = {
        "median_us" + suffix, "mad_us" + suffix,
        "minimum_us" + suffix, "maximum_us" + suffix}
    if require_samples:
        expected_keys.add("samples_us" + suffix)
    if execution:
        expected_keys.update((input_rate, output_rate))
    require_exact_keys(summary, expected_keys, prefix)
    values = {
        name: finite_positive(summary.get(name + suffix), f"{prefix}.{name + suffix}",
                              allow_zero=(name == "mad_us"))
        for name in ("median_us", "mad_us", "minimum_us", "maximum_us")
    }
    if not (values["minimum_us"] <= values["median_us"] <= values["maximum_us"]):
        raise ComparisonError(f"{prefix} min/median/max order is invalid")
    sample_name = "samples_us" + suffix
    if require_samples:
        samples_raw = summary.get(sample_name)
        if not isinstance(samples_raw, list) or len(samples_raw) != iterations:
            raise ComparisonError(f"{prefix}.{sample_name} cardinality changed")
        samples = [finite_positive(value, f"{prefix}.{sample_name}")
                   for value in samples_raw]
        derived = {
            "median_us": median(samples), "mad_us": mad(samples),
            "minimum_us": min(samples), "maximum_us": max(samples),
        }
        for name, expected in derived.items():
            if not close_enough(values[name], expected):
                raise ComparisonError(f"{prefix}.{name + suffix} is not raw-derived")
    if execution:
        expected_input_rate = input_bytes / (values["median_us"] * 1000.0)
        observed_input = finite_positive(summary.get(input_rate), f"{prefix}.{input_rate}")
        if not close_enough(observed_input, expected_input_rate):
            raise ComparisonError(f"{prefix} throughput is not byte/time-derived")
        if output_bytes == 0:
            if summary.get(output_rate) is not None:
                raise ComparisonError(f"{prefix}.{output_rate} must be null for zero output")
        else:
            expected_output_rate = output_bytes / (values["median_us"] * 1000.0)
            observed_output = finite_positive(summary.get(output_rate),
                                              f"{prefix}.{output_rate}")
            if not close_enough(observed_output, expected_output_rate):
                raise ComparisonError(f"{prefix} throughput is not byte/time-derived")


def validate_common_parameters(
        document: Mapping[str, Any], cell: Mapping[str, Any],
        iterations: int, warmup: int,
        extra_expected: Mapping[str, Any] | None = None) -> list[int]:
    parameters = document.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ComparisonError("result parameters missing")
    expected = expected_parameters(cell, iterations, warmup)
    if extra_expected:
        expected.update(extra_expected)
    require_exact_keys(parameters, set(expected) | {"missing_original_indices"},
                       "result parameters")
    for name, value in expected.items():
        if parameters.get(name) != value:
            raise ComparisonError(
                f"result parameter {name}={parameters.get(name)!r}, expected {value!r}")
    missing = parameters.get("missing_original_indices")
    if (not isinstance(missing, list) or len(missing) != expected["loss_count"] or
            missing != sorted(set(missing)) or
            any(isinstance(index, bool) or not isinstance(index, int) or
                index < 0 or index >= expected["K"] for index in missing)):
        raise ComparisonError("missing-original indices are invalid")
    return missing


def validate_isal_result(
        document: Mapping[str, Any], cell: Mapping[str, Any],
        iterations: int, warmup: int,
        expected_library_sha256: str | None = None,
        expected_nasm_sha256: str | None = None,
        expected_oracle_mode: str = "full") -> list[int]:
    if expected_oracle_mode not in ("full", "projection"):
        raise ComparisonError("invalid expected ISA-L oracle mode")
    require_exact_keys(document, {
        "schema", "provider", "parameters", "comparison_identity",
        "correctness", "memory", "metrics"}, "ISA-L child result")
    if document.get("schema") != ISA_SCHEMA:
        raise ComparisonError("wrong ISA-L result schema")
    provider = document.get("provider", {})
    require_exact_keys(provider, {
        "name", "source_commit", "source_tree", "library_sha256", "header_sha256", "license",
        "license_sha256", "nasm_executable_sha256", "nasm_archive_sha256",
        "field", "generator", "wire_compatible"}, "ISA-L child provider")
    required_provider = {
        "name": "Intel ISA-L", "source_commit": ISA_COMMIT,
        "source_tree": ISA_TREE,
        "license": "BSD-3-Clause", "field": "gf8",
        "license_sha256": ISA_LICENSE_SHA256,
        "header_sha256": ISA_HEADER_SHA256,
        "nasm_archive_sha256": NASM_ARCHIVE_SHA256,
        "generator": "gf_gen_cauchy1_matrix", "wire_compatible": False,
    }
    if any(provider.get(name) != value for name, value in required_provider.items()):
        raise ComparisonError("ISA-L provider identity/generator/wire claim changed")
    library_hash = provider.get("library_sha256")
    require_hex(library_hash, 64, "ISA-L child library digest")
    nasm_hash = require_hex(
        provider.get("nasm_executable_sha256"), 64, "ISA-L child NASM digest")
    if expected_library_sha256 is not None and library_hash != expected_library_sha256:
        raise ComparisonError("ISA-L child library is not provenance-derived")
    if expected_nasm_sha256 is not None and nasm_hash != expected_nasm_sha256:
        raise ComparisonError("ISA-L child NASM identity is not provenance-derived")
    missing = validate_common_parameters(document, cell, iterations, warmup)
    profile, field, parent = resolved_identity(cell)
    identity = document.get("comparison_identity", {})
    require_exact_keys(identity, {
        "leopard2_profile", "leopard2_field", "leopard2_parent",
        "isa_l_field_advantage_from_padding", "scope"},
        "ISA-L comparison identity")
    expected_advantage = field == "gf16" and int(cell["K"]) + int(cell["R"]) <= 256
    if (identity.get("leopard2_profile") != profile or
            identity.get("leopard2_field") != field or
            identity.get("leopard2_parent") != parent or
            identity.get("isa_l_field_advantage_from_padding") is not expected_advantage or
            identity.get("scope") !=
            "public payload and repaired-output throughput only; parity bytes and "
            "generator matrices differ"):
        raise ComparisonError("ISA-L comparison identity/qualification changed")
    correctness = document.get("correctness", {})
    require_exact_keys(correctness, {
        "direct_source_round_trip", "systematic_generator_prefix",
        "systematic_sources_immutable",
        "independent_scalar_cauchy_coefficients_checked",
        "independent_scalar_parity_mode",
        "independent_scalar_parity_checked_bytes_per_validation",
        "independent_scalar_parity_total_bytes_per_validation",
        "independent_scalar_parity_validation_passes"},
        "ISA-L correctness evidence")
    k, r, shard_bytes = int(cell["K"]), int(cell["R"]), int(cell["shard_bytes"])
    batch = int(cell["batch"])
    total_parity_bytes = r * shard_bytes * batch
    checked_parity_bytes = (
        total_parity_bytes if expected_oracle_mode == "full" else
        r * min(shard_bytes, 64) * batch)
    if (correctness.get("direct_source_round_trip") is not True or
            correctness.get("systematic_generator_prefix") is not True or
            correctness.get("systematic_sources_immutable") is not True or
            correctness.get("independent_scalar_cauchy_coefficients_checked") != k * r or
            correctness.get("independent_scalar_parity_mode") != expected_oracle_mode or
            correctness.get("independent_scalar_parity_checked_bytes_per_validation") !=
            checked_parity_bytes or
            correctness.get("independent_scalar_parity_total_bytes_per_validation") !=
            total_parity_bytes or
            correctness.get("independent_scalar_parity_validation_passes") != 2):
        raise ComparisonError("ISA-L direct-source correctness failed")
    memory = document.get("memory", {})
    losses = int(cell["loss_count"])
    expected_memory = {
        "alignment_bytes": 64, "direct_application_buffers": True,
        "staging_copy_bytes_per_execution": 0,
        "encode_input_bytes_per_batch_call": k * shard_bytes * batch,
        "encode_output_bytes_per_batch_call": r * shard_bytes * batch,
        "decode_offered_bytes_per_batch_call": (k - losses + r) * shard_bytes * batch,
        "decode_selected_bytes_per_batch_call": k * shard_bytes * batch,
        "decode_output_bytes_per_batch_call": losses * shard_bytes * batch,
        "matrix_bytes": (k + r) * k,
        "encode_table_bytes": 32 * k * r,
        "decode_table_bytes": 32 * k * losses,
    }
    require_exact_keys(memory, set(expected_memory), "ISA-L memory evidence")
    if any(memory.get(name) != value for name, value in expected_memory.items()):
        raise ComparisonError("ISA-L memory/transfer semantics changed")
    metrics = document.get("metrics", {})
    require_exact_keys(metrics, {
        "codec_setup", "encode_execution", "decode_plan_setup",
        "decode_execution", "decode_amortized_at_reuse", "rate_semantics"},
        "ISA-L metrics")
    if metrics.get("rate_semantics") != (
            "offered_received counts all non-erased public shards; ISA-L reads "
            "the recorded deterministic K-row subset"):
        raise ComparisonError("ISA-L rate semantics changed")
    validate_summary(metrics.get("codec_setup", {}), "codec_setup", iterations,
                     False, require_samples=True)
    validate_summary(
        metrics.get("encode_execution", {}), "encode_execution", iterations, True,
        expected_memory["encode_input_bytes_per_batch_call"],
        expected_memory["encode_output_bytes_per_batch_call"],
        "input_GB_per_s", "parity_output_GB_per_s", True)
    validate_summary(metrics.get("decode_plan_setup", {}), "decode_plan_setup",
                     iterations, False, require_samples=True)
    validate_summary(
        metrics.get("decode_execution", {}), "decode_execution", iterations, True,
        expected_memory["decode_offered_bytes_per_batch_call"],
        expected_memory["decode_output_bytes_per_batch_call"],
        "offered_received_GB_per_s", "repaired_output_GB_per_s", True)
    setup = float(metrics["decode_plan_setup"]["median_us"])
    execution = float(metrics["decode_execution"]["median_us_per_batch_call"])
    expected_amortized = execution + setup / int(cell["reuse"])
    amortized = metrics.get("decode_amortized_at_reuse", {})
    require_exact_keys(amortized, {
        "reuse_count", "derived_median_us_per_batch_call",
        "offered_received_GB_per_s", "repaired_output_GB_per_s"},
        "ISA-L amortized decode")
    if (amortized.get("reuse_count") != int(cell["reuse"]) or
            not close_enough(finite_positive(
                amortized.get("derived_median_us_per_batch_call"),
                "decode amortized time"), expected_amortized)):
        raise ComparisonError("ISA-L amortized decode is not setup-derived")
    expected_input_rate = expected_memory["decode_offered_bytes_per_batch_call"] / (
        expected_amortized * 1000.0)
    observed_input_rate = finite_positive(
        amortized.get("offered_received_GB_per_s"), "ISA-L amortized input rate")
    if not close_enough(observed_input_rate, expected_input_rate):
        raise ComparisonError("ISA-L amortized input rate is not byte/time-derived")
    if expected_memory["decode_output_bytes_per_batch_call"] == 0:
        if amortized.get("repaired_output_GB_per_s") is not None:
            raise ComparisonError("ISA-L zero-output amortized rate must be null")
    else:
        expected_output_rate = expected_memory["decode_output_bytes_per_batch_call"] / (
            expected_amortized * 1000.0)
        observed_output_rate = finite_positive(
            amortized.get("repaired_output_GB_per_s"),
            "ISA-L amortized output rate")
        if not close_enough(observed_output_rate, expected_output_rate):
            raise ComparisonError("ISA-L amortized output rate is not byte/time-derived")
    return missing


def validate_leopard_result(
        document: Mapping[str, Any], cell: Mapping[str, Any],
        iterations: int, warmup: int) -> list[int]:
    if document.get("schema") != LEOPARD_SCHEMA:
        raise ComparisonError("wrong Leopard2 result schema")
    require_exact_keys(document, {
        "schema", "build", "parameters", "resolved", "correctness", "memory",
        "metrics", "legacy"}, "Leopard2 child result")
    missing = validate_common_parameters(
        document, cell, iterations, warmup,
        {"force_generic_decode": False, "force_specialized_decode": False})
    build = document.get("build", {})
    require_exact_keys(build, {"compiler", "compiler_version", "cplusplus"},
                       "Leopard2 child build")
    if (not isinstance(build.get("compiler"), str) or not build["compiler"] or
            not isinstance(build.get("compiler_version"), str) or
            not build["compiler_version"] or
            isinstance(build.get("cplusplus"), bool) or
            not isinstance(build.get("cplusplus"), int)):
        raise ComparisonError("Leopard2 child build identity is invalid")
    profile, field, parent = resolved_identity(cell)
    resolved = document.get("resolved", {})
    require_exact_keys(resolved, {
        "profile", "field", "backend", "thread_count", "parent_count",
        "padded_side"}, "Leopard2 resolved identity")
    if (resolved.get("profile") != profile or resolved.get("field") != field or
            resolved.get("parent_count") != parent or
            resolved.get("thread_count") != 1):
        raise ComparisonError("Leopard2 resolved identity changed")
    correctness = document.get("correctness", {})
    require_exact_keys(correctness, {"leopard2_round_trip", "legacy_comparison"},
                       "Leopard2 correctness")
    if correctness.get("leopard2_round_trip") is not True:
        raise ComparisonError("Leopard2 direct-source correctness failed")
    memory = document.get("memory", {})
    require_exact_keys(memory, {
        "scratch_alignment", "encode_scratch_bytes_per_stripe",
        "decode_scratch_bytes_per_stripe", "encode_scratch_bytes_batch",
        "decode_scratch_bytes_batch"}, "Leopard2 memory evidence")
    if any(isinstance(value, bool) or not isinstance(value, int) or value < 0
           for value in memory.values()):
        raise ComparisonError("Leopard2 memory evidence is invalid")
    legacy = document.get("legacy", {})
    require_exact_keys(legacy, {
        "available", "unavailable_reason", "codec_setup",
        "decode_timing_includes_setup", "encode_execution",
        "decode_including_setup"}, "Leopard2 legacy evidence")
    if (not isinstance(legacy.get("available"), bool) or
            legacy.get("decode_timing_includes_setup") is not True or
            legacy.get("codec_setup") is not None):
        raise ComparisonError("Leopard2 legacy timing contract changed")
    metrics = document.get("metrics", {})
    require_exact_keys(metrics, {
        "codec_setup", "encode_execution", "decode_plan_setup",
        "decode_execution", "decode_amortized_at_reuse", "rate_semantics"},
        "Leopard2 metrics")
    if metrics.get("rate_semantics") != (
            "offered_received counts all non-null shard pointers supplied; "
            "a plan may read a deterministic subset"):
        raise ComparisonError("Leopard2 rate semantics changed")
    k, r, shard_bytes = int(cell["K"]), int(cell["R"]), int(cell["shard_bytes"])
    losses, batch = int(cell["loss_count"]), int(cell["batch"])
    validate_summary(metrics.get("codec_setup", {}), "leopard codec setup", iterations,
                     False, require_samples=True)
    validate_summary(
        metrics.get("encode_execution", {}), "leopard encode", iterations, True,
        k * shard_bytes * batch, r * shard_bytes * batch,
        "input_GB_per_s", "parity_output_GB_per_s", True)
    validate_summary(metrics.get("decode_plan_setup", {}), "leopard plan setup",
                     iterations, False, require_samples=True)
    validate_summary(
        metrics.get("decode_execution", {}), "leopard decode", iterations, True,
        (k - losses + r) * shard_bytes * batch, losses * shard_bytes * batch,
        "offered_received_GB_per_s", "repaired_output_GB_per_s", True)
    setup = float(metrics["decode_plan_setup"]["median_us"])
    execution = float(metrics["decode_execution"]["median_us_per_batch_call"])
    expected_amortized = execution + setup / int(cell["reuse"])
    amortized = metrics.get("decode_amortized_at_reuse", {})
    require_exact_keys(amortized, {
        "reuse_count", "derived_median_us_per_batch_call",
        "offered_received_GB_per_s", "repaired_output_GB_per_s"},
        "Leopard2 amortized decode")
    if (amortized.get("reuse_count") != int(cell["reuse"]) or
            not close_enough(finite_positive(
                amortized.get("derived_median_us_per_batch_call"),
                "Leopard2 amortized decode"), expected_amortized)):
        raise ComparisonError("Leopard2 amortized decode is not setup-derived")
    expected_input_rate = (k - losses + r) * shard_bytes * batch / (
        expected_amortized * 1000.0)
    if not close_enough(finite_positive(
            amortized.get("offered_received_GB_per_s"),
            "Leopard2 amortized input rate"), expected_input_rate):
        raise ComparisonError("Leopard2 amortized input rate is not byte/time-derived")
    if losses == 0:
        if amortized.get("repaired_output_GB_per_s") is not None:
            raise ComparisonError("Leopard2 zero-output amortized rate must be null")
    else:
        expected_output_rate = losses * shard_bytes * batch / (
            expected_amortized * 1000.0)
        if not close_enough(finite_positive(
                amortized.get("repaired_output_GB_per_s"),
                "Leopard2 amortized output rate"), expected_output_rate):
            raise ComparisonError(
                "Leopard2 amortized output rate is not byte/time-derived")
    return missing


def aggregate_results(results: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    aggregate: list[dict[str, Any]] = []
    for cell_index, cell in enumerate(CHECKPOINT_CELLS):
        entry: dict[str, Any] = {"cell_index": cell_index, "cell": dict(cell), "providers": {}}
        for provider in ("isa-l", "leopard2"):
            matches = [result["document"] for result in results
                       if result["cell_index"] == cell_index and result["provider"] == provider]
            if len(matches) != 4:
                raise ComparisonError("aggregate requires four provider repetitions per cell")
            provider_summary: dict[str, Any] = {}
            for metric, key in (
                    ("codec_setup", "median_us"),
                    ("encode_execution", "median_us_per_batch_call"),
                    ("decode_plan_setup", "median_us"),
                    ("decode_execution", "median_us_per_batch_call"),
                    ("decode_amortized_at_reuse",
                     "derived_median_us_per_batch_call")):
                values = [float(document["metrics"][metric][key]) for document in matches]
                provider_summary[metric] = {
                    "median_of_run_medians_us": median(values),
                    "mad_of_run_medians_us": mad(values),
                }
                if metric in ("encode_execution", "decode_execution",
                              "decode_amortized_at_reuse"):
                    if metric == "encode_execution":
                        input_bytes = int(cell["K"]) * int(cell["shard_bytes"]) * int(cell["batch"])
                        selected_input_bytes = input_bytes
                        output_bytes = int(cell["R"]) * int(cell["shard_bytes"]) * int(cell["batch"])
                    else:
                        input_bytes = (int(cell["K"]) - int(cell["loss_count"]) +
                                       int(cell["R"])) * int(cell["shard_bytes"]) * int(cell["batch"])
                        selected_input_bytes = (int(cell["K"]) *
                                                int(cell["shard_bytes"]) *
                                                int(cell["batch"]))
                        output_bytes = int(cell["loss_count"]) * int(cell["shard_bytes"]) * int(cell["batch"])
                    elapsed = provider_summary[metric]["median_of_run_medians_us"]
                    if metric == "encode_execution":
                        provider_summary[metric]["input_GB_per_s"] = (
                            input_bytes / (elapsed * 1000.0))
                    else:
                        provider_summary[metric]["offered_input_GB_per_s"] = (
                            input_bytes / (elapsed * 1000.0))
                        provider_summary[metric]["selected_input_GB_per_s"] = (
                            selected_input_bytes / (elapsed * 1000.0))
                    provider_summary[metric]["output_GB_per_s"] = (
                        output_bytes / (elapsed * 1000.0) if output_bytes else None)
            entry["providers"][provider] = provider_summary
        entry["isa_l_speedup_vs_leopard2"] = {
            metric: (entry["providers"]["leopard2"][metric]["median_of_run_medians_us"] /
                     entry["providers"]["isa-l"][metric]["median_of_run_medians_us"])
            for metric in ("encode_execution", "decode_execution",
                           "decode_amortized_at_reuse")
        }
        aggregate.append(entry)
    return aggregate


def allowed_cpus() -> set[int]:
    if hasattr(os, "sched_getaffinity"):
        return set(os.sched_getaffinity(0))
    raise ComparisonError("sched_getaffinity is required for pinned evidence")


@contextmanager
def advisory_cpu_lease(cache: Path, cpu: int, sibling: int):
    lease_directory = cache / "leases"
    lease_directory.mkdir(parents=True, exist_ok=True)
    handles = []
    lock_names = [f"cpu-{value}.lock" for value in sorted((cpu, sibling))]
    try:
        for name in lock_names:
            handle = (lease_directory / name).open("a+", encoding="utf-8")
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                handle.close()
                raise ComparisonError(
                    f"coordinator CPU lease is already held: {name}") from error
            handle.seek(0)
            handle.truncate()
            handle.write(f"pid={os.getpid()}\n")
            handle.flush()
            handles.append(handle)
        yield {
            "schema": LEASE_SCHEMA,
            "mechanism": "fcntl-flock-exclusive-advisory",
            "scope": "cooperating Leopard2 lab jobs only; not an OS-exclusive CPU reservation",
            "cpus": sorted((cpu, sibling)),
            "lock_names": lock_names,
            "coordinator_pid": os.getpid(),
            "acquired_before_affinity": True,
            "held_through_measurement": True,
            "held_through_post_timing_integrity": True,
            "os_exclusive": False,
        }
    finally:
        for handle in reversed(handles):
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            handle.close()


def parse_cpu_list(text: str) -> list[int]:
    cpus: set[int] = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start_text, end_text = part.split("-", 1)
            start, end = int(start_text), int(end_text)
            if start < 0 or end < start:
                raise ComparisonError(f"invalid CPU range {part!r}")
            cpus.update(range(start, end + 1))
        else:
            value = int(part)
            if value < 0:
                raise ComparisonError(f"invalid CPU {part!r}")
            cpus.add(value)
    return sorted(cpus)


def optional_text(path: Path) -> str | None:
    try:
        return path.read_text().strip()
    except OSError:
        return None


def optional_khz(path: Path) -> int | None:
    value = optional_text(path)
    if value is None:
        return None
    try:
        parsed = int(value)
    except ValueError as error:
        raise ComparisonError(f"invalid kHz value in {path}: {value!r}") from error
    if parsed <= 0:
        raise ComparisonError(f"nonpositive kHz value in {path}: {parsed}")
    return parsed


def current_frequency_snapshot(cpu: int, phase: str) -> dict[str, Any]:
    if phase not in ("pre", "post"):
        raise ComparisonError("frequency snapshot phase must be pre or post")
    frequency_root = Path(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq")
    return {
        "captured_phase": phase,
        "scaling_cur_freq_khz": optional_khz(frequency_root / "scaling_cur_freq"),
        "cpuinfo_cur_freq_khz": optional_khz(frequency_root / "cpuinfo_cur_freq"),
    }


def cpu_metadata(
        cpu: int, original_allowed: set[int], reserved_idle_cpu: int,
        child_affinity_preflight: list[int],
        lease: Mapping[str, Any], pre_frequency: Mapping[str, Any],
        post_frequency: Mapping[str, Any]) -> dict[str, Any]:
    sibling_path = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
    frequency_root = Path(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq")
    cpu_info: dict[str, str] = {}
    try:
        for block in Path("/proc/cpuinfo").read_text().split("\n\n"):
            candidate: dict[str, str] = {}
            for line in block.splitlines():
                if ":" in line:
                    name, value = line.split(":", 1)
                    candidate[name.strip()] = value.strip()
            try:
                processor = int(candidate.get("processor", ""))
            except ValueError:
                continue
            if processor == cpu:
                cpu_info = candidate
                break
    except OSError:
        pass
    if not cpu_info:
        raise ComparisonError(f"/proc/cpuinfo has no metadata block for CPU {cpu}")
    siblings_text = optional_text(sibling_path)
    siblings = parse_cpu_list(siblings_text) if siblings_text else [cpu]
    return {
        "requested_cpu": cpu,
        "allowed_cpu_count": len(original_allowed),
        "allowed_cpus": sorted(original_allowed),
        "process_affinity_during_run": sorted(allowed_cpus()),
        "child_affinity_preflight": child_affinity_preflight,
        "thread_siblings_list": siblings_text,
        "thread_siblings": siblings,
        "reserved_idle_sibling_cpu": reserved_idle_cpu,
        "coordinator_lease": dict(lease),
        "scaling_driver": optional_text(frequency_root / "scaling_driver"),
        "scaling_governor": optional_text(frequency_root / "scaling_governor"),
        "energy_performance_preference": optional_text(
            frequency_root / "energy_performance_preference"),
        "cpuinfo_min_freq_khz": optional_khz(frequency_root / "cpuinfo_min_freq"),
        "cpuinfo_max_freq_khz": optional_khz(frequency_root / "cpuinfo_max_freq"),
        "scaling_min_freq_khz": optional_khz(frequency_root / "scaling_min_freq"),
        "scaling_max_freq_khz": optional_khz(frequency_root / "scaling_max_freq"),
        "current_frequency_pre": dict(pre_frequency),
        "current_frequency_post": dict(post_frequency),
        "amd_pstate_status": optional_text(
            Path("/sys/devices/system/cpu/amd_pstate/status")),
        "boost": optional_text(Path("/sys/devices/system/cpu/cpufreq/boost")),
        "no_turbo": optional_text(
            Path("/sys/devices/system/cpu/intel_pstate/no_turbo")),
        "cpuinfo_processor": int(cpu_info["processor"]),
        "cpu_vendor": cpu_info.get("vendor_id"),
        "cpu_model_name": cpu_info.get("model name"),
        "cpu_family": cpu_info.get("cpu family"),
        "cpu_model": cpu_info.get("model"),
        "cpu_stepping": cpu_info.get("stepping"),
        "microcode": cpu_info.get("microcode"),
        "platform": platform.platform(),
        "uname": list(platform.uname()),
        "python": sys.version.splitlines()[0],
        "pinning": "runner sched_setaffinity singleton; all children inherit",
        "parallelism_note": (
            "authoritative cache-sensitive measurements intentionally use one physical "
            "core; bootstrap compilation is capped at eight jobs"),
    }


def child_command(
        provider: str, executable: Path, cell: Mapping[str, Any],
        iterations: int, warmup: int, output: Path,
        oracle_mode: str | None = None) -> list[str]:
    command = [
        str(executable), "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--profile", str(cell["profile"]), "--bytes", str(cell["shard_bytes"]),
        "--loss", str(cell["loss_count"]), "--batch", str(cell["batch"]),
        "--reuse", str(cell["reuse"]), "--iterations", str(iterations),
        "--warmup", str(warmup), "--threads", "1", "--seed", str(cell["seed"]),
        "--json", str(output)]
    if provider == "leopard2":
        command[command.index("--bytes"):command.index("--bytes")] = [
            "--field", "auto", "--backend", "auto"]
        command.insert(-2, "--retain-samples")
    else:
        if oracle_mode not in ("full", "projection"):
            raise ComparisonError("ISA-L child requires an explicit oracle mode")
        command[-2:-2] = ["--oracle", oracle_mode]
    return command


class XorShift64:
    def __init__(self, seed: int):
        self.state = seed or 0x9E3779B97F4A7C15

    def next(self) -> int:
        value = self.state
        value ^= (value << 13) & 0xFFFFFFFFFFFFFFFF
        value ^= value >> 7
        value ^= (value << 17) & 0xFFFFFFFFFFFFFFFF
        self.state = value & 0xFFFFFFFFFFFFFFFF
        return self.state


def correctness_cells(count: int) -> list[dict[str, Any]]:
    if count < 16 or count > 512:
        raise ComparisonError("--cases must be in 16..512")
    boundaries = [
        (1, 1), (1, 255), (255, 1), (2, 254), (8, 248), (16, 240),
        (64, 192), (127, 129), (128, 128), (129, 100), (191, 64),
        (192, 64), (200, 50), (224, 32), (225, 30), (240, 16),
    ]
    byte_sizes = (1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64,
                  65, 127, 128, 129, 255, 256, 257, 1024, 4096)
    cells: list[dict[str, Any]] = []
    for index, (k, r) in enumerate(boundaries):
        maximum = min(k, r)
        loss = (0 if index % 4 == 0 else 1 if index % 4 == 1 else
                maximum if index % 4 == 2 else maximum // 2 or 1)
        cells.append({
            "K": k, "R": r, "profile": "low" if r > k else "high",
            "shard_bytes": byte_sizes[index % len(byte_sizes)],
            "loss_count": loss, "batch": 1 + index % 3,
            "reuse": 1 + index % 4, "seed": 0x1A2B0000 + index,
        })
    random = XorShift64(0x4C454F324953414C)
    while len(cells) < count:
        k = 1 + random.next() % 255
        r = 1 + random.next() % (256 - k)
        maximum = min(k, r)
        cells.append({
            "K": k, "R": r,
            "profile": ("auto" if len(cells) % 3 == 0 else
                        "low" if r > k else "high"),
            "shard_bytes": byte_sizes[random.next() % len(byte_sizes)],
            "loss_count": random.next() % (maximum + 1),
            "batch": 1 + random.next() % 3,
            "reuse": 1 + random.next() % 4,
            "seed": random.next(),
        })
    return cells


def portable_sources(provenance: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "isa_l": {name: provenance["isa_l"][name] for name in (
            "name", "version", "url", "commit", "tree", "license",
            "license_sha256", "header_sha256", "library_sha256")},
        "nasm": {name: provenance["nasm"][name] for name in (
            "version", "url", "archive_sha256", "executable_sha256",
            "reported_version")},
        "leopard": dict(provenance["executables"]["leopard_source"]),
        "source_bundle": provenance["executables"]["source_bundle"],
        "build_identity": provenance["executables"]["build_identity"],
    }


def validate_portable_sources(sources: Mapping[str, Any], verify_local: bool) -> None:
    require_exact_keys(sources, {
        "isa_l", "nasm", "leopard", "source_bundle", "build_identity"},
        "portable benchmark sources")
    isa = sources.get("isa_l", {})
    require_exact_keys(isa, {"name", "version", "url", "commit", "tree", "license",
                             "license_sha256", "header_sha256", "library_sha256"},
                       "portable ISA-L source")
    if (isa.get("name") != "Intel ISA-L" or isa.get("version") != ISA_VERSION or
            isa.get("url") != ISA_URL or isa.get("commit") != ISA_COMMIT or
            isa.get("tree") != ISA_TREE or isa.get("license") != "BSD-3-Clause" or
            isa.get("license_sha256") != ISA_LICENSE_SHA256 or
            isa.get("header_sha256") != ISA_HEADER_SHA256):
        raise ComparisonError("portable ISA-L identity changed")
    require_hex(isa.get("library_sha256"), 64, "portable ISA-L library digest")
    nasm = sources.get("nasm", {})
    require_exact_keys(nasm, {"version", "url", "archive_sha256", "executable_sha256",
                              "reported_version"}, "portable NASM source")
    if (nasm.get("version") != NASM_VERSION or nasm.get("url") != NASM_URL or
            nasm.get("archive_sha256") != NASM_ARCHIVE_SHA256 or
            f"NASM version {NASM_VERSION}" not in str(nasm.get("reported_version", ""))):
        raise ComparisonError("portable NASM identity changed")
    require_hex(nasm.get("executable_sha256"), 64, "portable NASM executable digest")
    leopard = sources.get("leopard", {})
    require_exact_keys(leopard, {"commit", "tree", "dirty", "clean_at_build"},
                       "portable Leopard source")
    require_hex(leopard.get("commit"), 40, "portable Leopard commit")
    require_hex(leopard.get("tree"), 40, "portable Leopard tree")
    if leopard.get("dirty") is not False or leopard.get("clean_at_build") is not True:
        raise ComparisonError("portable Leopard source was not clean")
    bundle = sources.get("source_bundle", {})
    validate_source_bundle(bundle, verify_local=verify_local)
    if (leopard.get("commit") != bundle.get("source_commit") or
            leopard.get("tree") != bundle.get("source_tree")):
        raise ComparisonError("portable Leopard identity is not source-bundle-derived")
    validate_build_identity(sources.get("build_identity", {}), nasm)
    commands = sources["build_identity"]["compile_commands"]
    for manifest_name in ("adapter", "leopard2"):
        for entry in commands[manifest_name]["entries"]:
            if entry["file"].startswith("${LEOPARD_ROOT}/"):
                relative = entry["file"][len("${LEOPARD_ROOT}/"):]
                if relative not in bundle["entries"]:
                    raise ComparisonError(
                        "Leopard compile command source is absent from source bundle")


def run_correctness(cache: Path, count: int, output: Path) -> dict[str, Any]:
    provenance = load_provenance(cache)
    executable = Path(provenance["executables"]["isa-l"]["path"])
    cells = correctness_cells(count)
    results: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="leo2-isal-correctness-", dir=str(cache)) as temporary:
        temporary_path = Path(temporary)
        for index, cell in enumerate(cells):
            result_path = temporary_path / f"case{index}.json"
            command = child_command(
                "isa-l", executable, cell, iterations=3, warmup=0,
                output=result_path, oracle_mode="full")
            completed = run_checked(command, timeout=120)
            if completed.stdout or completed.stderr:
                raise ComparisonError(f"correctness child {index} emitted output")
            document = json.loads(result_path.read_text())
            validate_isal_result(
                document, cell, 3, 0,
                provenance["isa_l"]["library_sha256"],
                provenance["nasm"]["executable_sha256"])
            results.append({
                "case_index": index, "cell": cell,
                "executable_sha256": provenance["executables"]["isa-l"]["sha256"],
                "document": document})
    artifact: dict[str, Any] = {
        "schema": CORRECTNESS_SCHEMA,
        "sources": portable_sources(provenance),
        "executables": {
            "isa-l": {"sha256": provenance["executables"]["isa-l"]["sha256"]}},
        "case_generation_seed": "0x4c454f324953414c",
        "case_count": count,
        "no_loss_cases": sum(result["cell"]["loss_count"] == 0 for result in results),
        "maximum_loss_cases": sum(
            result["cell"]["loss_count"] == min(
                result["cell"]["K"], result["cell"]["R"])
            for result in results),
        "padding_boundary_gf16_comparisons": sum(
            resolved_identity(result["cell"])[1] == "gf16" for result in results),
        "results": results,
    }
    artifact["artifact_sha256"] = canonical_digest(artifact)
    validate_correctness(artifact)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    return artifact


def validate_correctness(
        document: Mapping[str, Any], verify_local: bool = False) -> None:
    require_exact_keys(document, {
        "schema", "sources", "executables", "case_generation_seed", "case_count",
        "no_loss_cases", "maximum_loss_cases", "padding_boundary_gf16_comparisons",
        "results", "artifact_sha256"}, "correctness artifact")
    if document.get("schema") != CORRECTNESS_SCHEMA:
        raise ComparisonError("wrong correctness artifact identity")
    if document.get("case_generation_seed") != "0x4c454f324953414c":
        raise ComparisonError("correctness case-generation seed changed")
    sources = document.get("sources", {})
    validate_portable_sources(sources, verify_local=verify_local)
    executables = document.get("executables", {})
    require_exact_keys(executables, {"isa-l"}, "correctness executables")
    require_exact_keys(executables.get("isa-l", {}), {"sha256"},
                       "correctness ISA-L executable")
    require_hex(executables["isa-l"].get("sha256"), 64,
                "correctness ISA-L executable digest")
    count = document.get("case_count")
    if isinstance(count, bool) or not isinstance(count, int):
        raise ComparisonError("invalid correctness case count")
    expected_cells = correctness_cells(count)
    results = document.get("results")
    if not isinstance(results, list) or len(results) != count:
        raise ComparisonError("correctness result cardinality changed")
    for index, (result, cell) in enumerate(zip(results, expected_cells)):
        require_exact_keys(result, {
            "case_index", "cell", "executable_sha256", "document"},
            "correctness result")
        if (result.get("case_index") != index or result.get("cell") != cell or
                result.get("executable_sha256") != executables["isa-l"]["sha256"]):
            raise ComparisonError("correctness result order/cell changed")
        validate_isal_result(
            result.get("document", {}), cell, 3, 0,
            sources["isa_l"]["library_sha256"],
            sources["nasm"]["executable_sha256"])
    derived_counts = {
        "no_loss_cases": sum(cell["loss_count"] == 0 for cell in expected_cells),
        "maximum_loss_cases": sum(
            cell["loss_count"] == min(cell["K"], cell["R"])
            for cell in expected_cells),
        "padding_boundary_gf16_comparisons": sum(
            resolved_identity(cell)[1] == "gf16" for cell in expected_cells),
    }
    if any(document.get(name) != value for name, value in derived_counts.items()):
        raise ComparisonError("correctness coverage counts are not case-derived")
    if document.get("artifact_sha256") != canonical_digest(document):
        raise ComparisonError("correctness artifact digest mismatch")


def run_checkpoint(
        cache: Path, cpu: int, reserved_idle_cpu: int, output: Path,
        iterations: int, warmup: int) -> dict[str, Any]:
    original_affinity = allowed_cpus()
    if cpu not in original_affinity:
        raise ComparisonError(
            f"CPU {cpu} is outside allowed affinity {sorted(original_affinity)}")
    sibling_path = Path(
        f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
    try:
        siblings = parse_cpu_list(sibling_path.read_text().strip())
    except (OSError, ValueError) as error:
        raise ComparisonError(f"cannot verify CPU {cpu} sibling topology: {error}") from error
    if (cpu not in siblings or reserved_idle_cpu == cpu or
            reserved_idle_cpu not in siblings or reserved_idle_cpu not in original_affinity):
        raise ComparisonError(
            f"reserved idle CPU {reserved_idle_cpu} is not an allowed SMT sibling "
            f"of CPU {cpu}: {siblings}")
    if iterations < 5 or warmup < 1:
        raise ComparisonError("checkpoint requires at least 5 samples and 1 warmup")
    provenance = load_provenance(cache)
    executables = {
        provider: Path(provenance["executables"][provider]["path"])
        for provider in ("isa-l", "leopard2")}
    results: list[dict[str, Any]] = []
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
        "GOMP_CPU_AFFINITY": str(cpu), "LC_ALL": "C",
    })
    with advisory_cpu_lease(cache, cpu, reserved_idle_cpu) as lease:
        os.sched_setaffinity(0, {cpu})
        try:
            if allowed_cpus() != {cpu}:
                raise ComparisonError("runner failed to establish singleton affinity")
            preflight = run_checked([
                sys.executable, "-c",
                "import json,os;print(json.dumps(sorted(os.sched_getaffinity(0))))"],
                env=environment, timeout=30)
            if preflight.stderr:
                raise ComparisonError(
                    f"child-affinity preflight emitted stderr: {preflight.stderr}")
            child_affinity = json.loads(preflight.stdout)
            if child_affinity != [cpu]:
                raise ComparisonError(
                    f"child did not inherit singleton affinity: {child_affinity}")
            pre_frequency = current_frequency_snapshot(cpu, "pre")
            with tempfile.TemporaryDirectory(
                    prefix="leo2-isal-", dir=str(cache)) as temporary:
                temporary_path = Path(temporary)
                for cell_index, cell in enumerate(CHECKPOINT_CELLS):
                    for repetition, order in enumerate(ABBA):
                        missing_by_provider: dict[str, list[int]] = {}
                        for order_index, provider in enumerate(order):
                            result_path = temporary_path / (
                                f"cell{cell_index}.rep{repetition}.order{order_index}.{provider}.json")
                            command = child_command(
                                provider, executables[provider], cell, iterations, warmup,
                                result_path,
                                oracle_mode="projection" if provider == "isa-l" else None)
                            completed = run_checked(
                                command, env=environment, timeout=600)
                            if allowed_cpus() != {cpu}:
                                raise ComparisonError(
                                    "runner affinity changed during the measurement window")
                            if completed.stdout or completed.stderr:
                                raise ComparisonError(
                                    f"{provider} emitted unexpected output: "
                                    f"stdout={completed.stdout!r}, stderr={completed.stderr!r}")
                            try:
                                document = json.loads(result_path.read_text())
                            except (OSError, json.JSONDecodeError) as error:
                                raise ComparisonError(
                                    f"invalid {provider} JSON: {error}") from error
                            missing = (validate_isal_result(
                                document, cell, iterations, warmup,
                                provenance["isa_l"]["library_sha256"],
                                provenance["nasm"]["executable_sha256"],
                                "projection")
                                if provider == "isa-l" else
                                validate_leopard_result(
                                    document, cell, iterations, warmup))
                            missing_by_provider[provider] = missing
                            results.append({
                                "cell_index": cell_index, "repetition": repetition,
                                "order_index": order_index, "provider": provider,
                                "executable_sha256": provenance[
                                    "executables"][provider]["sha256"],
                                "document": document,
                            })
                        if missing_by_provider["isa-l"] != missing_by_provider["leopard2"]:
                            raise ComparisonError(
                                "providers did not use the same erasure pattern")
            post_timing_provenance = load_provenance(cache)
            if post_timing_provenance != provenance:
                raise ComparisonError(
                    "source, toolchain, library, or executable provenance changed "
                    "during the measurement window")
            post_frequency = current_frequency_snapshot(cpu, "post")
            host = cpu_metadata(
                cpu, original_affinity, reserved_idle_cpu, child_affinity, lease,
                pre_frequency, post_frequency)
        finally:
            os.sched_setaffinity(0, original_affinity)
    if allowed_cpus() != original_affinity:
        raise ComparisonError("runner failed to restore its original affinity")
    checkpoint: dict[str, Any] = {
        "schema": SCHEMA,
        "method": {
            "provider_order": "ABBA",
            "provider_order_by_repetition": [list(order) for order in ABBA],
            "repetitions": 4, "iterations_per_child": iterations,
            "warmup_per_child": warmup, "thread_count": 1,
            "codec_setup_separate": True, "decode_plan_setup_separate": True,
            "execution_reuses_setup": True,
            "post_timing_integrity_verified": True,
            "isa_l_timing_oracle": (
                "all K*R Cauchy coefficients plus deterministic projection of at most "
                "64 byte positions per parity shard and stripe, before and after timing"),
            "full_correctness_oracle": (
                "separate correctness artifact recomputes every ISA-L parity byte from "
                "independently derived scalar Cauchy coefficients"),
            "wire_compatible": False,
            "comparison_scope": (
                "same public K,R, shard bytes, deterministic source/erasure pattern, "
                "batch, reuse, thread count, and offered/generated/repaired bytes"),
            "decode_input_rate_semantics": (
                "offered counts every non-erased public shard pointer; selected counts "
                "the deterministic K-row subset actually consumed"),
        },
        "sources": portable_sources(provenance),
        "executables": {
            provider: {"sha256": provenance["executables"][provider]["sha256"]}
            for provider in ("isa-l", "leopard2")},
        "host": host,
        "cells": [dict(cell) for cell in CHECKPOINT_CELLS],
        "results": results,
    }
    checkpoint["aggregate"] = aggregate_results(results)
    checkpoint["remaining_gates"] = [
        "run the signed required 7,134-cell matrix where ISA-L is eligible",
        "implement persistent-pool batch scheduling before comparing thread counts above one",
        "repeat on additional x86 microarchitectures and report ISA dispatch",
        "collect available PMU/memory-bandwidth counters without weakening raw evidence checks",
        "compare 64/128-core aggregate batches on hardware exposing those CPUs",
    ]
    checkpoint["artifact_sha256"] = canonical_digest(checkpoint)
    validate_checkpoint(checkpoint)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")
    return checkpoint


def validate_advisory_lease(
        lease: Mapping[str, Any], requested_cpu: int, reserved_cpu: int) -> None:
    require_exact_keys(lease, {
        "schema", "mechanism", "scope", "cpus", "lock_names", "coordinator_pid",
        "acquired_before_affinity", "held_through_measurement",
        "held_through_post_timing_integrity", "os_exclusive"},
        "coordinator advisory lease")
    expected_cpus = sorted((requested_cpu, reserved_cpu))
    if (lease.get("schema") != LEASE_SCHEMA or
            lease.get("mechanism") != "fcntl-flock-exclusive-advisory" or
            lease.get("scope") !=
            "cooperating Leopard2 lab jobs only; not an OS-exclusive CPU reservation" or
            lease.get("cpus") != expected_cpus or
            lease.get("lock_names") != [f"cpu-{cpu}.lock" for cpu in expected_cpus] or
            isinstance(lease.get("coordinator_pid"), bool) or
            not isinstance(lease.get("coordinator_pid"), int) or
            lease.get("coordinator_pid") <= 0 or
            lease.get("acquired_before_affinity") is not True or
            lease.get("held_through_measurement") is not True or
            lease.get("held_through_post_timing_integrity") is not True or
            lease.get("os_exclusive") is not False):
        raise ComparisonError("coordinator advisory lease is invalid")


def validate_frequency_snapshot(value: Mapping[str, Any], phase: str) -> None:
    require_exact_keys(value, {
        "captured_phase", "scaling_cur_freq_khz", "cpuinfo_cur_freq_khz"},
        f"{phase} frequency snapshot")
    if value.get("captured_phase") != phase:
        raise ComparisonError(f"{phase} frequency snapshot phase changed")
    for name in ("scaling_cur_freq_khz", "cpuinfo_cur_freq_khz"):
        frequency = value.get(name)
        if (frequency is not None and
                (isinstance(frequency, bool) or not isinstance(frequency, int) or
                 frequency <= 0)):
            raise ComparisonError(f"{phase} frequency snapshot {name} is not positive kHz")


def validate_checkpoint(
        document: Mapping[str, Any], verify_local: bool = False) -> None:
    require_exact_keys(document, {
        "schema", "method", "sources", "executables", "host", "cells", "results",
        "aggregate", "remaining_gates", "artifact_sha256"}, "checkpoint")
    if document.get("schema") != SCHEMA:
        raise ComparisonError("wrong checkpoint schema")
    method = document.get("method", {})
    require_exact_keys(method, {
        "provider_order", "provider_order_by_repetition", "repetitions",
        "iterations_per_child", "warmup_per_child", "thread_count",
        "codec_setup_separate", "decode_plan_setup_separate",
        "execution_reuses_setup", "post_timing_integrity_verified",
        "isa_l_timing_oracle", "full_correctness_oracle",
        "wire_compatible", "comparison_scope",
        "decode_input_rate_semantics"}, "checkpoint method")
    if (method.get("provider_order") != "ABBA" or
            method.get("provider_order_by_repetition") != [list(order) for order in ABBA] or
            method.get("repetitions") != 4 or method.get("thread_count") != 1 or
            method.get("wire_compatible") is not False or
            method.get("codec_setup_separate") is not True or
            method.get("decode_plan_setup_separate") is not True or
            method.get("execution_reuses_setup") is not True or
            method.get("post_timing_integrity_verified") is not True or
            method.get("isa_l_timing_oracle") !=
            "all K*R Cauchy coefficients plus deterministic projection of at most "
            "64 byte positions per parity shard and stripe, before and after timing" or
            method.get("full_correctness_oracle") !=
            "separate correctness artifact recomputes every ISA-L parity byte from "
            "independently derived scalar Cauchy coefficients" or
            method.get("comparison_scope") !=
            "same public K,R, shard bytes, deterministic source/erasure pattern, "
            "batch, reuse, thread count, and offered/generated/repaired bytes" or
            method.get("decode_input_rate_semantics") !=
            "offered counts every non-erased public shard pointer; selected counts "
            "the deterministic K-row subset actually consumed"):
        raise ComparisonError("checkpoint method/fairness contract changed")
    iterations = method.get("iterations_per_child")
    warmup = method.get("warmup_per_child")
    if (isinstance(iterations, bool) or not isinstance(iterations, int) or iterations < 5 or
            isinstance(warmup, bool) or not isinstance(warmup, int) or warmup < 1):
        raise ComparisonError("checkpoint sample/warmup contract changed")
    if document.get("cells") != [dict(cell) for cell in CHECKPOINT_CELLS]:
        raise ComparisonError("checkpoint cell matrix changed")
    sources = document.get("sources", {})
    validate_portable_sources(sources, verify_local=verify_local)
    executables = document.get("executables", {})
    require_exact_keys(executables, {"isa-l", "leopard2"}, "checkpoint executables")
    for provider in ("isa-l", "leopard2"):
        executable = executables.get(provider, {})
        require_exact_keys(executable, {"sha256"}, f"checkpoint {provider} executable")
        require_hex(executable.get("sha256"), 64,
                    f"checkpoint {provider} executable digest")
    host = document.get("host", {})
    require_exact_keys(host, {
        "requested_cpu", "allowed_cpu_count", "allowed_cpus",
        "process_affinity_during_run", "child_affinity_preflight",
        "thread_siblings_list", "thread_siblings", "reserved_idle_sibling_cpu",
        "coordinator_lease", "scaling_driver", "scaling_governor",
        "energy_performance_preference", "cpuinfo_min_freq_khz",
        "cpuinfo_max_freq_khz", "scaling_min_freq_khz", "scaling_max_freq_khz",
        "current_frequency_pre", "current_frequency_post", "amd_pstate_status",
        "boost", "no_turbo", "cpuinfo_processor", "cpu_vendor",
        "cpu_model_name", "cpu_family", "cpu_model", "cpu_stepping", "microcode",
        "platform", "uname", "python", "pinning", "parallelism_note"},
        "checkpoint host")
    requested_cpu = host.get("requested_cpu")
    reserved_cpu = host.get("reserved_idle_sibling_cpu")
    siblings = host.get("thread_siblings")
    allowed = host.get("allowed_cpus")
    if (isinstance(requested_cpu, bool) or not isinstance(requested_cpu, int) or
            isinstance(reserved_cpu, bool) or not isinstance(reserved_cpu, int) or
            requested_cpu == reserved_cpu or
            not isinstance(allowed, list) or
            any(isinstance(cpu, bool) or not isinstance(cpu, int) for cpu in allowed) or
            allowed != sorted(set(allowed)) or
            host.get("allowed_cpu_count") != len(allowed) or
            requested_cpu not in allowed or reserved_cpu not in allowed or
            not isinstance(siblings, list) or requested_cpu not in siblings or
            reserved_cpu not in siblings or
            not isinstance(host.get("thread_siblings_list"), str) or
            parse_cpu_list(host["thread_siblings_list"]) != siblings or
            host.get("process_affinity_during_run") != [requested_cpu] or
            host.get("child_affinity_preflight") != [requested_cpu] or
            host.get("cpuinfo_processor") != requested_cpu or
            host.get("pinning") !=
            "runner sched_setaffinity singleton; all children inherit"):
        raise ComparisonError("checkpoint affinity evidence is invalid")
    validate_advisory_lease(host.get("coordinator_lease", {}), requested_cpu, reserved_cpu)
    for name in (
            "scaling_driver", "scaling_governor", "energy_performance_preference",
            "amd_pstate_status", "boost", "no_turbo"):
        value = host.get(name)
        if value is not None and (not isinstance(value, str) or not value):
            raise ComparisonError(f"checkpoint host {name} is neither text nor null")
    for name in (
            "cpuinfo_min_freq_khz", "cpuinfo_max_freq_khz",
            "scaling_min_freq_khz", "scaling_max_freq_khz"):
        value = host.get(name)
        if (value is not None and
                (isinstance(value, bool) or not isinstance(value, int) or value <= 0)):
            raise ComparisonError(f"checkpoint host {name} is not positive kHz or null")
    for minimum, maximum in (
            (host.get("cpuinfo_min_freq_khz"), host.get("cpuinfo_max_freq_khz")),
            (host.get("scaling_min_freq_khz"), host.get("scaling_max_freq_khz"))):
        if minimum is not None and maximum is not None and minimum > maximum:
            raise ComparisonError("checkpoint host frequency range is inverted")
    validate_frequency_snapshot(host.get("current_frequency_pre", {}), "pre")
    validate_frequency_snapshot(host.get("current_frequency_post", {}), "post")
    results = document.get("results")
    if not isinstance(results, list) or len(results) != len(CHECKPOINT_CELLS) * 8:
        raise ComparisonError("checkpoint result cardinality changed")
    seen: set[tuple[int, int, int, str]] = set()
    for result in results:
        if not isinstance(result, Mapping):
            raise ComparisonError("checkpoint result is not an object")
        require_exact_keys(result, {
            "cell_index", "repetition", "order_index", "provider",
            "executable_sha256", "document"},
            "checkpoint result")
        cell_index = result.get("cell_index")
        repetition = result.get("repetition")
        order_index = result.get("order_index")
        provider = result.get("provider")
        if (isinstance(cell_index, bool) or not isinstance(cell_index, int) or
                cell_index < 0 or cell_index >= len(CHECKPOINT_CELLS) or
                isinstance(repetition, bool) or not isinstance(repetition, int) or
                repetition < 0 or repetition >= 4 or
                isinstance(order_index, bool) or not isinstance(order_index, int) or
                order_index not in (0, 1) or
                provider not in ("isa-l", "leopard2") or
                ABBA[repetition][order_index] != provider):
            raise ComparisonError("checkpoint result schedule is invalid")
        if result.get("executable_sha256") != executables[provider]["sha256"]:
            raise ComparisonError("checkpoint child executable identity differs")
        identity = (cell_index, repetition, order_index, provider)
        if identity in seen:
            raise ComparisonError("duplicate checkpoint result")
        seen.add(identity)
        cell = CHECKPOINT_CELLS[cell_index]
        if provider == "isa-l":
            validate_isal_result(
                result.get("document", {}), cell, iterations, warmup,
                sources["isa_l"]["library_sha256"],
                sources["nasm"]["executable_sha256"], "projection")
        else:
            validate_leopard_result(result.get("document", {}), cell, iterations, warmup)
    for cell_index in range(len(CHECKPOINT_CELLS)):
        for repetition, order in enumerate(ABBA):
            documents = {result["provider"]: result["document"] for result in results
                         if result["cell_index"] == cell_index and
                         result["repetition"] == repetition}
            if (documents["isa-l"]["parameters"]["missing_original_indices"] !=
                    documents["leopard2"]["parameters"]["missing_original_indices"]):
                raise ComparisonError("paired erasure patterns differ")
    derived_aggregate = aggregate_results(results)
    if document.get("aggregate") != derived_aggregate:
        raise ComparisonError("checkpoint aggregate is not raw-result-derived")
    remaining = document.get("remaining_gates")
    if not isinstance(remaining, list) or len(remaining) < 4:
        raise ComparisonError("checkpoint must retain full-matrix limitations")
    if document.get("artifact_sha256") != canonical_digest(document):
        raise ComparisonError("checkpoint artifact digest mismatch")


def fake_summary(execution: bool, input_bytes: int = 0, output_bytes: int = 0,
                 input_rate: str = "", output_rate: str = "",
                 iterations: int = 5) -> dict[str, Any]:
    samples = [10.0 + index for index in range(iterations)]
    center = median(samples)
    suffix = "_per_batch_call" if execution else ""
    result = {
        "median_us" + suffix: center, "mad_us" + suffix: mad(samples),
        "minimum_us" + suffix: min(samples), "maximum_us" + suffix: max(samples),
        "samples_us" + suffix: samples,
    }
    if execution:
        result[input_rate] = input_bytes / (center * 1000.0)
        result[output_rate] = output_bytes / (center * 1000.0) if output_bytes else None
    return result


def fake_isal_result(
        cell: Mapping[str, Any], iterations: int = 5, warmup: int = 1,
        oracle_mode: str = "full") -> dict[str, Any]:
    params = expected_parameters(cell, iterations, warmup)
    params["missing_original_indices"] = list(range(params["loss_count"]))
    k, r, b, losses, batch = (params["K"], params["R"], params["shard_bytes"],
                              params["loss_count"], params["batch"])
    profile, field, parent = resolved_identity(cell)
    encode_input, encode_output = k * b * batch, r * b * batch
    decode_offered = (k - losses + r) * b * batch
    decode_output = losses * b * batch
    setup = fake_summary(False, iterations=iterations)
    encode = fake_summary(True, encode_input, encode_output,
                          "input_GB_per_s", "parity_output_GB_per_s", iterations)
    decode_setup = fake_summary(False, iterations=iterations)
    decode = fake_summary(True, decode_offered, decode_output,
                          "offered_received_GB_per_s", "repaired_output_GB_per_s",
                          iterations)
    center = float(setup["median_us"])
    amortized = center + center / params["reuse"]
    return {
        "schema": ISA_SCHEMA,
        "provider": {"name": "Intel ISA-L", "source_commit": ISA_COMMIT,
                     "source_tree": ISA_TREE,
                     "library_sha256": "1" * 64,
                     "header_sha256": ISA_HEADER_SHA256,
                     "license": "BSD-3-Clause",
                     "license_sha256": ISA_LICENSE_SHA256,
                     "nasm_executable_sha256": "4" * 64,
                     "nasm_archive_sha256": NASM_ARCHIVE_SHA256,
                     "field": "gf8", "generator": "gf_gen_cauchy1_matrix",
                     "wire_compatible": False},
        "parameters": params,
        "comparison_identity": {
            "leopard2_profile": profile, "leopard2_field": field,
            "leopard2_parent": parent,
            "isa_l_field_advantage_from_padding": field == "gf16",
            "scope": (
                "public payload and repaired-output throughput only; parity bytes and "
                "generator matrices differ")},
        "correctness": {"direct_source_round_trip": True,
                        "systematic_generator_prefix": True,
                        "systematic_sources_immutable": True,
                        "independent_scalar_cauchy_coefficients_checked": k * r,
                        "independent_scalar_parity_mode": oracle_mode,
                        "independent_scalar_parity_checked_bytes_per_validation": (
                            r * b * batch if oracle_mode == "full" else
                            r * min(b, 64) * batch),
                        "independent_scalar_parity_total_bytes_per_validation":
                            r * b * batch,
                        "independent_scalar_parity_validation_passes": 2},
        "memory": {"alignment_bytes": 64, "direct_application_buffers": True,
                   "staging_copy_bytes_per_execution": 0,
                   "encode_input_bytes_per_batch_call": encode_input,
                   "encode_output_bytes_per_batch_call": encode_output,
                   "decode_offered_bytes_per_batch_call": decode_offered,
                   "decode_selected_bytes_per_batch_call": k * b * batch,
                   "decode_output_bytes_per_batch_call": decode_output,
                   "matrix_bytes": (k + r) * k,
                   "encode_table_bytes": 32 * k * r,
                   "decode_table_bytes": 32 * k * losses},
        "metrics": {"codec_setup": setup, "encode_execution": encode,
                    "decode_plan_setup": decode_setup, "decode_execution": decode,
                    "decode_amortized_at_reuse": {
                        "reuse_count": params["reuse"],
                        "derived_median_us_per_batch_call": amortized,
                        "offered_received_GB_per_s": decode_offered / (amortized * 1000.0),
                        "repaired_output_GB_per_s": (
                            decode_output / (amortized * 1000.0)
                            if decode_output else None)},
                    "rate_semantics": (
                        "offered_received counts all non-erased public shards; ISA-L reads "
                        "the recorded deterministic K-row subset")},
    }


def fake_leopard_result(
        cell: Mapping[str, Any], iterations: int = 5,
        warmup: int = 1) -> dict[str, Any]:
    params = expected_parameters(cell, iterations, warmup)
    params["force_generic_decode"] = False
    params["force_specialized_decode"] = False
    params["missing_original_indices"] = list(range(params["loss_count"]))
    k, r, b, losses, batch = (params["K"], params["R"], params["shard_bytes"],
                              params["loss_count"], params["batch"])
    profile, field, parent = resolved_identity(cell)
    encode_input, encode_output = k * b * batch, r * b * batch
    decode_offered = (k - losses + r) * b * batch
    decode_output = losses * b * batch
    center = median([10.0 + index for index in range(iterations)])
    amortized = center + center / params["reuse"]
    return {
        "schema": LEOPARD_SCHEMA,
        "build": {"compiler": "test", "compiler_version": "test",
                  "cplusplus": 201103},
        "parameters": params,
        "resolved": {"profile": profile, "field": field, "parent_count": parent,
                     "thread_count": 1, "backend": "avx2", "padded_side": 1},
        "correctness": {"leopard2_round_trip": True,
                        "legacy_comparison": None},
        "memory": {"scratch_alignment": 64,
                   "encode_scratch_bytes_per_stripe": 0,
                   "decode_scratch_bytes_per_stripe": 0,
                   "encode_scratch_bytes_batch": 0,
                   "decode_scratch_bytes_batch": 0},
        "metrics": {
            "codec_setup": fake_summary(False, iterations=iterations),
            "encode_execution": fake_summary(
                True, encode_input, encode_output,
                "input_GB_per_s", "parity_output_GB_per_s", iterations),
            "decode_plan_setup": fake_summary(False, iterations=iterations),
            "decode_execution": fake_summary(
                True, decode_offered, decode_output,
                "offered_received_GB_per_s", "repaired_output_GB_per_s", iterations),
            "decode_amortized_at_reuse": {
                "reuse_count": params["reuse"],
                "derived_median_us_per_batch_call": amortized,
                "offered_received_GB_per_s": decode_offered / (amortized * 1000.0),
                "repaired_output_GB_per_s": (
                    decode_output / (amortized * 1000.0)
                    if decode_output else None),
            },
            "rate_semantics": (
                "offered_received counts all non-null shard pointers supplied; "
                "a plan may read a deterministic subset"),
        },
        "legacy": {"available": False, "unavailable_reason": "test",
                   "codec_setup": None, "decode_timing_includes_setup": True,
                   "encode_execution": None, "decode_including_setup": None},
    }


def fake_source_bundle() -> dict[str, Any]:
    return source_bundle(repo_root())


def fake_compile_manifest(name: str) -> dict[str, Any]:
    file_name = ("${ISA_L_SOURCE}/erasure_code/ec_base.c" if name == "isa_l" else
                 "${LEOPARD_ROOT}/leopard2.cpp")
    entries = [{
        "directory": "${" + name.upper() + "_BUILD}",
        "file": file_name,
        "command": "${CXX} -c " + file_name,
    }]
    encoded = json.dumps(entries, sort_keys=True, separators=(",", ":"))
    return {
        "schema": "leopard2-normalized-compile-commands/v1",
        "entries": entries, "entry_count": 1,
        "sha256": hashlib.sha256(encoded.encode("utf-8")).hexdigest(),
    }


def fake_build_identity() -> dict[str, Any]:
    identity = {
        "schema": TOOLCHAIN_SCHEMA,
        "recipe": {
            "isa_l": {
                "configure_definitions": [
                    "CMAKE_BUILD_TYPE=Release", "BUILD_SHARED_LIBS=OFF",
                    "ISAL_BUILD_TESTS=OFF", "ISAL_BUILD_PERF_TESTS=OFF",
                    "ISAL_BUILD_ISAL_SHIM=OFF", "ISAL_BUILD_IGZIP_CLI=OFF",
                    "CMAKE_EXPORT_COMPILE_COMMANDS=ON"],
                "generator": "Unix Makefiles",
                "build_target": "all", "install_target": "install"},
            "adapter": {
                "configure_definitions": [
                    "CMAKE_BUILD_TYPE=Release", "CMAKE_EXPORT_COMPILE_COMMANDS=ON",
                    "LEO2_ISAL_EXPECTED_COMMIT=" + ISA_COMMIT,
                    "LEO2_ISAL_EXPECTED_LICENSE_SHA256=" + ISA_LICENSE_SHA256,
                    "LEO2_ISAL_EXPECTED_HEADER_SHA256=" + ISA_HEADER_SHA256,
                    "LEO2_NASM_EXECUTABLE_SHA256=" + "4" * 64,
                    "LEO2_NASM_ARCHIVE_SHA256=" + NASM_ARCHIVE_SHA256],
                "generator": "Unix Makefiles",
                "build_target": "leopard2_isal_benchmark"},
            "leopard2": {
                "configure_definitions": [
                    "CMAKE_BUILD_TYPE=Release", "LEO2_BUILD_TESTS=OFF",
                    "LEO2_BUILD_BENCHMARKS=ON", "LEO2_BUILD_FUZZERS=OFF",
                    "LEO2_ENABLE_CUDA=OFF", "ENABLE_OPENMP=ON",
                    "CMAKE_EXPORT_COMPILE_COMMANDS=ON"],
                "generator": "Unix Makefiles",
                "build_target": "bench_leopard2"},
        },
        "tools": {
            "cmake": {"sha256": "5" * 64, "reported_version": "cmake test"},
            "cc": {"sha256": "6" * 64, "reported_version": "cc test"},
            "cxx": {"sha256": "7" * 64, "reported_version": "cxx test"},
            "ar": {"sha256": "8" * 64, "reported_version": "ar test"},
            "ranlib": {"sha256": "9" * 64, "reported_version": "ranlib test"},
            "linker": {"sha256": "a" * 64, "reported_version": "ld test"},
            "build_program": {"sha256": "b" * 64, "reported_version": "make test"},
            "nasm": {"sha256": "4" * 64,
                     "reported_version": f"NASM version {NASM_VERSION} test"},
        },
        "compile_commands": {
            "isa_l": fake_compile_manifest("isa_l"),
            "adapter": fake_compile_manifest("adapter"),
            "leopard2": fake_compile_manifest("leopard2"),
        },
    }
    identity["identity_sha256"] = keyed_digest(identity, "identity_sha256")
    return identity


def fake_portable_sources() -> dict[str, Any]:
    bundle = fake_source_bundle()
    return {
        "isa_l": {
            "name": "Intel ISA-L", "version": ISA_VERSION, "url": ISA_URL,
            "commit": ISA_COMMIT, "tree": ISA_TREE, "license": "BSD-3-Clause",
            "license_sha256": ISA_LICENSE_SHA256,
            "header_sha256": ISA_HEADER_SHA256,
            "library_sha256": "1" * 64,
        },
        "nasm": {
            "version": NASM_VERSION, "url": NASM_URL,
            "archive_sha256": NASM_ARCHIVE_SHA256,
            "executable_sha256": "4" * 64,
            "reported_version": f"NASM version {NASM_VERSION} test",
        },
        "leopard": {"commit": bundle["source_commit"],
                    "tree": bundle["source_tree"],
                    "dirty": False, "clean_at_build": True},
        "source_bundle": bundle,
        "build_identity": fake_build_identity(),
    }


def fake_checkpoint() -> dict[str, Any]:
    results: list[dict[str, Any]] = []
    for cell_index, cell in enumerate(CHECKPOINT_CELLS):
        for repetition, order in enumerate(ABBA):
            for order_index, provider in enumerate(order):
                document = (fake_isal_result(cell, oracle_mode="projection")
                            if provider == "isa-l" else
                            fake_leopard_result(cell))
                results.append({
                    "cell_index": cell_index, "repetition": repetition,
                    "order_index": order_index, "provider": provider,
                    "executable_sha256": (
                        "2" * 64 if provider == "isa-l" else "3" * 64),
                    "document": document,
                })
    checkpoint: dict[str, Any] = {
        "schema": SCHEMA,
        "method": {
            "provider_order": "ABBA",
            "provider_order_by_repetition": [list(order) for order in ABBA],
            "repetitions": 4, "iterations_per_child": 5,
            "warmup_per_child": 1, "thread_count": 1,
            "codec_setup_separate": True, "decode_plan_setup_separate": True,
            "execution_reuses_setup": True, "wire_compatible": False,
            "post_timing_integrity_verified": True,
            "isa_l_timing_oracle": (
                "all K*R Cauchy coefficients plus deterministic projection of at most "
                "64 byte positions per parity shard and stripe, before and after timing"),
            "full_correctness_oracle": (
                "separate correctness artifact recomputes every ISA-L parity byte from "
                "independently derived scalar Cauchy coefficients"),
            "comparison_scope": (
                "same public K,R, shard bytes, deterministic source/erasure pattern, "
                "batch, reuse, thread count, and offered/generated/repaired bytes"),
            "decode_input_rate_semantics": (
                "offered counts every non-erased public shard pointer; selected counts "
                "the deterministic K-row subset actually consumed"),
        },
        "sources": fake_portable_sources(),
        "executables": {
            "isa-l": {"sha256": "2" * 64},
            "leopard2": {"sha256": "3" * 64},
        },
        "host": {
            "requested_cpu": 7, "allowed_cpus": [7, 8],
            "allowed_cpu_count": 2,
            "process_affinity_during_run": [7],
            "child_affinity_preflight": [7],
            "thread_siblings_list": "7-8",
            "thread_siblings": [7, 8],
            "reserved_idle_sibling_cpu": 8,
            "coordinator_lease": {
                "schema": LEASE_SCHEMA,
                "mechanism": "fcntl-flock-exclusive-advisory",
                "scope": "cooperating Leopard2 lab jobs only; not an OS-exclusive CPU reservation",
                "cpus": [7, 8], "lock_names": ["cpu-7.lock", "cpu-8.lock"],
                "coordinator_pid": 123, "acquired_before_affinity": True,
                "held_through_measurement": True,
                "held_through_post_timing_integrity": True,
                "os_exclusive": False,
            },
            "scaling_driver": "test", "scaling_governor": "test",
            "energy_performance_preference": None,
            "cpuinfo_min_freq_khz": 1000000, "cpuinfo_max_freq_khz": 5000000,
            "scaling_min_freq_khz": 1000000, "scaling_max_freq_khz": 5000000,
            "current_frequency_pre": {
                "captured_phase": "pre", "scaling_cur_freq_khz": 4000000,
                "cpuinfo_cur_freq_khz": None},
            "current_frequency_post": {
                "captured_phase": "post", "scaling_cur_freq_khz": 4200000,
                "cpuinfo_cur_freq_khz": None},
            "amd_pstate_status": None, "boost": None, "no_turbo": None,
            "cpuinfo_processor": 7, "cpu_vendor": "test",
            "cpu_model_name": "test", "cpu_family": "1", "cpu_model": "1",
            "cpu_stepping": "1", "microcode": "1", "platform": "test",
            "uname": ["test"], "python": "test",
            "pinning": "runner sched_setaffinity singleton; all children inherit",
            "parallelism_note": (
                "authoritative cache-sensitive measurements intentionally use one physical "
                "core; bootstrap compilation is capped at eight jobs"),
        },
        "cells": [dict(cell) for cell in CHECKPOINT_CELLS],
        "results": results,
        "remaining_gates": ["full matrix", "multicore", "other CPUs", "PMU"],
    }
    checkpoint["aggregate"] = aggregate_results(results)
    checkpoint["artifact_sha256"] = canonical_digest(checkpoint)
    return checkpoint


def fake_correctness_artifact() -> dict[str, Any]:
    cells = correctness_cells(16)
    results = [{
        "case_index": index, "cell": cell, "executable_sha256": "2" * 64,
        "document": fake_isal_result(cell, iterations=3, warmup=0),
    } for index, cell in enumerate(cells)]
    artifact: dict[str, Any] = {
        "schema": CORRECTNESS_SCHEMA,
        "sources": fake_portable_sources(),
        "executables": {"isa-l": {"sha256": "2" * 64}},
        "case_generation_seed": "0x4c454f324953414c",
        "case_count": len(cells),
        "no_loss_cases": sum(cell["loss_count"] == 0 for cell in cells),
        "maximum_loss_cases": sum(
            cell["loss_count"] == min(cell["K"], cell["R"]) for cell in cells),
        "padding_boundary_gf16_comparisons": sum(
            resolved_identity(cell)[1] == "gf16" for cell in cells),
        "results": results,
    }
    artifact["artifact_sha256"] = canonical_digest(artifact)
    return artifact


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="leo2-isal-lease-test-") as temporary:
        with advisory_cpu_lease(Path(temporary), 7, 8) as lease:
            validate_advisory_lease(lease, 7, 8)
            try:
                with advisory_cpu_lease(Path(temporary), 7, 8):
                    raise ComparisonError("duplicate advisory CPU lease was acquired")
            except ComparisonError as error:
                if "already held" not in str(error):
                    raise
    for cell in CHECKPOINT_CELLS:
        validate_isal_result(
            fake_isal_result(cell), cell, 5, 1, "1" * 64, "4" * 64)
    base = fake_isal_result(CHECKPOINT_CELLS[4])
    mutations = []
    for path, value in (
            (("provider", "generator"), "gf_gen_rs_matrix"),
            (("provider", "wire_compatible"), True),
            (("provider", "source_commit"), "0" * 40),
            (("provider", "source_tree"), "0" * 40),
            (("provider", "library_sha256"), "0" * 64),
            (("provider", "header_sha256"), "0" * 64),
            (("provider", "license_sha256"), "0" * 64),
            (("provider", "nasm_executable_sha256"), "0" * 64),
            (("correctness", "direct_source_round_trip"), False),
            (("correctness", "systematic_sources_immutable"), False),
            (("correctness", "independent_scalar_cauchy_coefficients_checked"), 0),
            (("correctness", "independent_scalar_parity_mode"), "projection"),
            (("correctness", "independent_scalar_parity_checked_bytes_per_validation"), 1),
            (("parameters", "K"), 128),
            (("comparison_identity", "isa_l_field_advantage_from_padding"), False),
            (("memory", "staging_copy_bytes_per_execution"), 1),
            (("metrics", "encode_execution", "median_us_per_batch_call"), 13.0)):
        changed = copy.deepcopy(base)
        target: Any = changed
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        mutations.append(changed)
    for changed in mutations:
        try:
            validate_isal_result(
                changed, CHECKPOINT_CELLS[4], 5, 1, "1" * 64, "4" * 64)
        except ComparisonError:
            pass
        else:
            raise ComparisonError("ISA-L result mutation was accepted")
    checkpoint = fake_checkpoint()
    validate_checkpoint(checkpoint, verify_local=False)
    checkpoint_mutations = []
    changed = copy.deepcopy(checkpoint)
    changed["results"].pop()
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["sources"]["source_bundle"]["source_tree"] = "0" * 40
    changed["sources"]["leopard"]["tree"] = "0" * 40
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["results"][0]["provider"] = "leopard2"
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["aggregate"][0]["providers"]["isa-l"]["codec_setup"][
        "median_of_run_medians_us"] = 99.0
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["results"][0]["document"]["provider"]["generator"] = "gf_gen_rs_matrix"
    changed["aggregate"] = aggregate_results(changed["results"])
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["artifact_sha256"] = "0" * 64
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["sources"]["leopard"]["dirty"] = True
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    first_input = sorted(changed["sources"]["source_bundle"]["entries"])[0]
    changed["sources"]["source_bundle"]["entries"][first_input]["sha256"] = "0" * 64
    entries = changed["sources"]["source_bundle"]["entries"]
    encoded = json.dumps(entries, sort_keys=True, separators=(",", ":"))
    changed["sources"]["source_bundle"]["bundle_sha256"] = hashlib.sha256(
        encoded.encode("utf-8")).hexdigest()
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["host"]["child_affinity_preflight"] = [7, 8]
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    # Coordinated summary/rate/aggregate rewriting must fail because the raw
    # Leopard samples remain the authoritative inputs.
    changed = copy.deepcopy(checkpoint)
    leopard = next(result for result in changed["results"]
                   if result["provider"] == "leopard2")
    timing = leopard["document"]["metrics"]["encode_execution"]
    for name in ("median_us_per_batch_call", "mad_us_per_batch_call",
                 "minimum_us_per_batch_call", "maximum_us_per_batch_call"):
        timing[name] *= 1.75
    timing["input_GB_per_s"] /= 1.75
    timing["parity_output_GB_per_s"] /= 1.75
    changed["aggregate"] = aggregate_results(changed["results"])
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["sources"]["nasm"].pop("executable_sha256")
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["sources"]["isa_l"]["library_sha256"] = "a" * 64
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["sources"]["isa_l"]["license_sha256"] = "a" * 64
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["sources"]["isa_l"]["header_sha256"] = "a" * 64
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["sources"]["build_identity"]["tools"]["cxx"]["sha256"] = "a" * 64
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["method"]["post_timing_integrity_verified"] = False
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["sources"]["leopard"]["commit"] = "not-a-commit"
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["results"][0]["executable_sha256"] = "f" * 64
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["host"]["allowed_cpu_count"] = 0
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["host"]["current_frequency_post"]["scaling_cur_freq_khz"] = -1
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["host"]["coordinator_lease"]["os_exclusive"] = True
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["host"]["coordinator_lease"]["lock_names"] = ["forged"]
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    for changed in checkpoint_mutations:
        try:
            validate_checkpoint(changed, verify_local=False)
        except ComparisonError:
            pass
        else:
            raise ComparisonError("checkpoint mutation was accepted")
    correctness = fake_correctness_artifact()
    validate_correctness(correctness, verify_local=False)
    correctness_mutations = []
    changed = copy.deepcopy(correctness)
    changed["executables"]["isa-l"].pop("sha256")
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    changed = copy.deepcopy(correctness)
    changed["sources"]["leopard"]["commit"] = "not-a-commit"
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    changed = copy.deepcopy(correctness)
    changed["results"][0]["executable_sha256"] = "f" * 64
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    changed = copy.deepcopy(correctness)
    changed["results"][0]["document"]["correctness"][
        "independent_scalar_cauchy_coefficients_checked"] = 0
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    changed = copy.deepcopy(correctness)
    changed["results"][0]["document"]["correctness"][
        "independent_scalar_parity_mode"] = "projection"
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    for changed in correctness_mutations:
        try:
            validate_correctness(changed, verify_local=False)
        except ComparisonError:
            pass
        else:
            raise ComparisonError("correctness mutation was accepted")
    print(json.dumps({
        "cells": len(CHECKPOINT_CELLS),
        "result_mutations_rejected": len(mutations),
        "checkpoint_mutations_rejected": len(checkpoint_mutations),
        "correctness_mutations_rejected": len(correctness_mutations),
        "nasm_archive_sha256": NASM_ARCHIVE_SHA256, "status": "PASS",
    }, sort_keys=True))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    bootstrap_parser = subparsers.add_parser("bootstrap")
    bootstrap_parser.add_argument("--cache", type=Path, default=default_cache())
    bootstrap_parser.add_argument("--jobs", type=int, default=min(os.cpu_count() or 1, 8))
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--cache", type=Path, default=default_cache())
    run_parser.add_argument("--cpu", type=int, required=True)
    run_parser.add_argument("--reserved-idle-cpu", type=int, required=True)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--iterations", type=int, default=9)
    run_parser.add_argument("--warmup", type=int, default=2)
    correctness_parser = subparsers.add_parser("correctness")
    correctness_parser.add_argument("--cache", type=Path, default=default_cache())
    correctness_parser.add_argument("--cases", type=int, default=128)
    correctness_parser.add_argument("--output", type=Path, required=True)
    correctness_validate_parser = subparsers.add_parser("validate-correctness")
    correctness_validate_parser.add_argument("artifact", type=Path)
    correctness_validate_parser.add_argument(
        "--require-local-build-match", action="store_true")
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("checkpoint", type=Path)
    validate_parser.add_argument("--require-local-build-match", action="store_true")
    subparsers.add_parser("self-test")
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "bootstrap":
            output = bootstrap(arguments.cache.resolve(), arguments.jobs)
            print(json.dumps(output, indent=2, sort_keys=True))
        elif arguments.command == "run":
            output = run_checkpoint(
                arguments.cache.resolve(), arguments.cpu,
                arguments.reserved_idle_cpu, arguments.output.resolve(),
                arguments.iterations, arguments.warmup)
            print(json.dumps({
                "artifact_sha256": output["artifact_sha256"],
                "cells": len(output["cells"]), "results": len(output["results"]),
                "status": "PASS"}, sort_keys=True))
        elif arguments.command == "correctness":
            output = run_correctness(
                arguments.cache.resolve(), arguments.cases,
                arguments.output.resolve())
            print(json.dumps({
                "artifact_sha256": output["artifact_sha256"],
                "cases": output["case_count"],
                "maximum_loss_cases": output["maximum_loss_cases"],
                "no_loss_cases": output["no_loss_cases"], "status": "PASS",
            }, sort_keys=True))
        elif arguments.command == "validate-correctness":
            document = json.loads(arguments.artifact.read_text())
            validate_correctness(
                document, verify_local=arguments.require_local_build_match)
            print(json.dumps({
                "artifact_sha256": document["artifact_sha256"],
                "cases": document["case_count"], "status": "PASS",
            }, sort_keys=True))
        elif arguments.command == "validate":
            document = json.loads(arguments.checkpoint.read_text())
            validate_checkpoint(
                document, verify_local=arguments.require_local_build_match)
            print(json.dumps({
                "artifact_sha256": document["artifact_sha256"],
                "cells": len(document["cells"]), "results": len(document["results"]),
                "status": "PASS"}, sort_keys=True))
        else:
            self_test()
    except (ComparisonError, OSError, json.JSONDecodeError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
