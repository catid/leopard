#!/usr/bin/env python3
"""Build, validate, and run the default-off Jerasure comparison.

The production build never discovers Jerasure or GF-Complete.  ``bootstrap``
materializes the exact committed Leopard source, verifies two clean pinned
third-party trees, creates private static archives, and records compile, link,
static-input, runtime-library, tool, and executable closure.  ``correctness``
runs independent full-byte cases in parallel.  ``run`` is deliberately a
separate pinned ABBA timing path that requires the authoritative correctness
artifact and an idle SMT sibling lease.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import copy
import hashlib
import json
import math
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
import leopard2_isal_compare as common  # noqa: E402
sys.path.pop(0)


SCHEMA = "leopard2-jerasure-comparison-checkpoint/v1"
CORRECTNESS_SCHEMA = "leopard2-jerasure-correctness/v1"
CHILD_SCHEMA = "leopard2-jerasure-benchmark/v1"
BUILD_SCHEMA = "leopard2-jerasure-build-closure/v1"
PROVENANCE_SCHEMA = "leopard2-jerasure-toolchain/v1"
RUN_MANIFEST_SCHEMA = "leopard2-jerasure-run-manifest/v1"
CHILD_ARTIFACT_SCHEMA = "leopard2-jerasure-child-artifact/v1"
GROUP_ARTIFACT_SCHEMA = "leopard2-jerasure-abba-pair/v1"
HOST_GROUPS_SCHEMA = "leopard2-jerasure-host-groups/v1"

JERASURE_URL = "https://github.com/ceph/jerasure"
JERASURE_COMMIT = "de1739cc8483696506829b52e7fda4f6bb195e6a"
JERASURE_TREE = "fb98f85c548038a5ff294141f89603dda70dd423"
JERASURE_LICENSE_SHA256 = (
    "83b6b3ff237848fbccfa889bb52cfb13c331c1f83544b907617c7f8f31eb1769")
JERASURE_HEADER_SHA256 = (
    "a6332be8149e5094b7bdafe50935d8d45922ff52188a6ced6c5b1ae4225041c9")
REED_SOL_HEADER_SHA256 = (
    "6b6decdc41bf2e4849de0f12cbd492b1beccf8bb6494cee7f04a9e05a46726ec")

GF_COMPLETE_URL = "https://github.com/ceph/gf-complete"
GF_COMPLETE_COMMIT = "a6862d10c9db467148f20eef2c6445ac9afd94d8"
GF_COMPLETE_TREE = "5a13169b93b6e517184fbdf39033098b329d68a6"
GF_COMPLETE_LICENSE_SHA256 = (
    "cb9790699b4a3d56a43bba1dd859f7f41361cd224e8745a24eef933ea134a280")
GF_COMPLETE_HEADER_SHA256 = (
    "f2c3a9651262c053565b395296489b76815d5b5921c099400a1bc3daf1ba0e8d")

AUTHORITATIVE_CORRECTNESS_CASES = 128
EVIDENCE_ITERATIONS = 9
EVIDENCE_WARMUP = 2
MAX_K = 4096
MAX_R = 4096
MAX_MATRIX_COEFFICIENTS = 1 << 24
MAX_APPLICATION_BYTES = 8 << 30
MAX_PARENT = 65536
GF_COMPLETE_KERNEL_POLICY = "portable-compiler-defaults"
ABBA = (("jerasure", "leopard2"), ("leopard2", "jerasure"),
        ("leopard2", "jerasure"), ("jerasure", "leopard2"))
CHECKPOINT_CELLS = (
    {"K": 240, "R": 16, "profile": "high", "shard_bytes": 65536,
     "loss_count": 1, "batch": 1, "reuse": 8, "seed": 0x4A450101},
    {"K": 240, "R": 16, "profile": "high", "shard_bytes": 4096,
     "loss_count": 4, "batch": 8, "reuse": 8, "seed": 0x4A450102},
    {"K": 128, "R": 128, "profile": "high", "shard_bytes": 65536,
     "loss_count": 8, "batch": 1, "reuse": 8, "seed": 0x4A450103},
    {"K": 64, "R": 192, "profile": "low", "shard_bytes": 65536,
     "loss_count": 8, "batch": 1, "reuse": 8, "seed": 0x4A450104},
    {"K": 129, "R": 100, "profile": "high", "shard_bytes": 65536,
     "loss_count": 4, "batch": 1, "reuse": 8, "seed": 0x4A450105},
    {"K": 200, "R": 100, "profile": "high", "shard_bytes": 65536,
     "loss_count": 4, "batch": 1, "reuse": 8, "seed": 0x4A450106},
)


class ComparisonError(ValueError):
    pass


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
    payload = dict(document)
    payload.pop("artifact_sha256", None)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def mapping_digest(document: Mapping[str, Any]) -> str:
    encoded = json.dumps(document, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def atomic_write_json(path: Path, document: Mapping[str, Any]) -> None:
    """Durably replace a JSON artifact without exposing a partial document."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8") as output:
            json.dump(document, output, indent=2, sort_keys=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
        try:
            directory = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        except OSError:
            return
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def require_exact_keys(value: Mapping[str, Any], expected: set[str], label: str) -> None:
    if not isinstance(value, Mapping) or set(value) != expected:
        raise ComparisonError(
            f"{label} keys changed: expected {sorted(expected)}, "
            f"got {sorted(value) if isinstance(value, Mapping) else type(value).__name__}")


def require_hex(value: Any, length: int, label: str) -> str:
    if (not isinstance(value, str) or len(value) != length or
            any(character not in "0123456789abcdef" for character in value)):
        raise ComparisonError(f"{label} is not lowercase {length}-digit hex")
    return value


def run_checked(command: Sequence[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        list(command), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False, **kwargs)
    if completed.returncode != 0:
        raise ComparisonError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}")
    return completed


def reject_contaminating_environment() -> None:
    common.reject_contaminating_environment()
    names = sorted(name for name in os.environ if name.startswith("GF_COMPLETE_"))
    if names:
        raise ComparisonError(
            "ambient GF-Complete dispatch controls are forbidden: " + ", ".join(names))


def build_paths(cache: Path) -> dict[str, Path]:
    return {
        "jerasure_source": cache / "jerasure",
        "gf_source": cache / "gf-complete",
        "adapter_build": cache / "build" / "leopard2-jerasure-adapter",
        "leopard_build": cache / "build" / "leopard2-jerasure-leopard",
        "leopard_source": cache / "sources" / "leopard-jerasure-detached",
        "provenance": cache / "toolchains" / "jerasure-comparison-provenance.json",
    }


def git_identity(checkout: Path) -> tuple[str, str, str]:
    head = run_checked(["git", "-C", str(checkout), "rev-parse", "HEAD"]).stdout.strip()
    tree = run_checked(
        ["git", "-C", str(checkout), "rev-parse", "HEAD^{tree}"]).stdout.strip()
    dirty = run_checked([
        "git", "-C", str(checkout), "status", "--porcelain",
        "--untracked-files=all"]).stdout.strip()
    return head, tree, dirty


def clone_pinned(path: Path, url: str, commit: str) -> None:
    if not (path / ".git").is_dir():
        path.parent.mkdir(parents=True, exist_ok=True)
        run_checked(["git", "clone", "--filter=blob:none", "--no-checkout", url, str(path)],
                    timeout=180)
    run_checked(["git", "-C", str(path), "fetch", "origin", commit], timeout=180)
    run_checked(["git", "-C", str(path), "checkout", "--detach", commit], timeout=60)


def verify_source(
        path: Path, name: str, url: str, commit: str, tree: str,
        license_hash: str, headers: Mapping[str, str]) -> dict[str, Any]:
    head, actual_tree, dirty = git_identity(path)
    if head != commit or actual_tree != tree or dirty:
        raise ComparisonError(
            f"{name} must be clean at {commit}/{tree}; "
            f"got {head}/{actual_tree}, status={dirty!r}")
    license_path = path / "License.txt"
    if sha256_file(license_path) != license_hash:
        raise ComparisonError(f"{name} license digest changed")
    header_hashes = {}
    for relative, expected in headers.items():
        actual = sha256_file(path / relative)
        if actual != expected:
            raise ComparisonError(f"{name} header digest changed: {relative}")
        header_hashes[relative] = actual
    tracked = run_checked([
        "git", "-C", str(path), "ls-tree", "-r", "--full-tree", "HEAD"]).stdout.splitlines()
    tracked_sha256 = hashlib.sha256(
        ("\n".join(tracked) + "\n").encode("utf-8")).hexdigest()
    return {
        "name": name, "url": url, "commit": commit, "tree": tree,
        "license": "BSD-3-Clause", "license_sha256": license_hash,
        "headers": header_hashes, "tracked_tree_listing_sha256": tracked_sha256,
        "tracked_entry_count": len(tracked), "clean_at_build": True,
    }


def leopard_status(root: Path) -> str:
    return run_checked([
        "git", "-C", str(root), "status", "--porcelain", "--untracked-files=normal",
        "--", ".", ":(exclude)experiments/leopard2/jerasure_compare/correctness_result.json",
        ":(exclude)experiments/leopard2/jerasure_compare/checkpoint_result.json",
    ]).stdout.strip()


def leopard_tree_listing(root: Path, commit: str) -> tuple[str, int]:
    listing = run_checked([
        "git", "-C", str(root), "ls-tree", "-r", "--full-tree", commit]).stdout
    lines = listing.splitlines()
    return hashlib.sha256(listing.encode("utf-8")).hexdigest(), len(lines)


def verify_detached_leopard(
        source: Path, expected: Mapping[str, Any] | None = None) -> dict[str, Any]:
    head, tree, dirty = git_identity(source)
    symbolic = subprocess.run(
        ["git", "-C", str(source), "symbolic-ref", "-q", "HEAD"],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    listing_hash, listing_count = leopard_tree_listing(source, "HEAD")
    record = {
        "commit": head, "tree": tree, "clean_at_build": True,
        "materialization": "detached Git worktree", "detached_head": True,
        "tracked_tree_listing_sha256": listing_hash,
        "tracked_entry_count": listing_count,
    }
    if dirty or symbolic.returncode == 0:
        raise ComparisonError("Leopard benchmark source is not clean and detached")
    if expected is not None and record != expected:
        raise ComparisonError("detached Leopard source identity changed after bootstrap")
    return record


def current_source_bindings(paths: Mapping[str, Path]) -> dict[str, Any]:
    leopard = verify_detached_leopard(paths["leopard_source"])
    bindings: dict[str, Any] = {
        "leopard": {key: leopard[key] for key in (
            "commit", "tree", "tracked_tree_listing_sha256", "tracked_entry_count")}}
    for name, path in (("jerasure", paths["jerasure_source"]),
                       ("gf_complete", paths["gf_source"])):
        head, tree, dirty = git_identity(path)
        listing_hash, listing_count = leopard_tree_listing(path, "HEAD")
        if dirty:
            raise ComparisonError(f"{name} source changed after bootstrap")
        bindings[name] = {
            "commit": head, "tree": tree,
            "tracked_tree_listing_sha256": listing_hash,
            "tracked_entry_count": listing_count}
    return bindings


def materialize_leopard(paths: Mapping[str, Path]) -> dict[str, Any]:
    root = repo_root().resolve()
    status = leopard_status(root)
    if status:
        raise ComparisonError(
            "benchmark source must be committed before bootstrap; status=" + repr(status))
    head, tree, _ = git_identity(root)
    listing_hash, listing_count = leopard_tree_listing(root, head)
    destination = paths["leopard_source"].resolve()
    if destination.exists():
        completed = subprocess.run([
            "git", "-C", str(root), "worktree", "remove", "--force", str(destination)],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if completed.returncode != 0:
            shutil.rmtree(destination)
    run_checked(["git", "-C", str(root), "worktree", "prune"], timeout=30)
    destination.parent.mkdir(parents=True, exist_ok=True)
    run_checked([
        "git", "-C", str(root), "worktree", "add", "--detach", str(destination), head],
        timeout=120)
    record = verify_detached_leopard(destination)
    if (record["commit"] != head or record["tree"] != tree or
            record["tracked_tree_listing_sha256"] != listing_hash or
            record["tracked_entry_count"] != listing_count):
        raise ComparisonError("detached Leopard source differs from benchmark commit")
    return record


def cmake_cache(build: Path, key: str) -> str:
    prefix = key + ":"
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if line.startswith(prefix):
            return line.split("=", 1)[1]
    raise ComparisonError(f"CMake cache omits {key}: {build}")


def normalize_manifest_paths(paths: Mapping[str, Path]) -> dict[str, str]:
    return {
        str(paths["jerasure_source"].resolve()): "${JERASURE_SOURCE}",
        str(paths["gf_source"].resolve()): "${GF_COMPLETE_SOURCE}",
        str(paths["leopard_source"].resolve()): "${LEOPARD_SOURCE}",
        str(paths["adapter_build"].resolve()): "${ADAPTER_BUILD}",
        str(paths["leopard_build"].resolve()): "${LEOPARD_BUILD}",
    }


def build_identity(paths: Mapping[str, Path]) -> dict[str, Any]:
    adapter = paths["adapter_build"]
    leopard = paths["leopard_build"]
    cc = Path(cmake_cache(adapter, "CMAKE_C_COMPILER")).resolve()
    cxx = Path(cmake_cache(adapter, "CMAKE_CXX_COMPILER")).resolve()
    if cxx != Path(cmake_cache(leopard, "CMAKE_CXX_COMPILER")).resolve():
        raise ComparisonError("adapter and Leopard used different C++ compilers")
    ar = Path(cmake_cache(adapter, "CMAKE_AR")).resolve()
    ranlib = Path(cmake_cache(adapter, "CMAKE_RANLIB")).resolve()
    linker = Path(cmake_cache(adapter, "CMAKE_LINKER")).resolve()
    build_program = Path(cmake_cache(adapter, "CMAKE_MAKE_PROGRAM")).resolve()
    for key, expected in (("CMAKE_AR", ar), ("CMAKE_RANLIB", ranlib),
                          ("CMAKE_LINKER", linker),
                          ("CMAKE_MAKE_PROGRAM", build_program)):
        if Path(cmake_cache(leopard, key)).resolve() != expected:
            raise ComparisonError(f"builds used different {key}")
    ldd = Path(shutil.which("ldd") or "")
    readelf = Path(shutil.which("readelf") or "")
    replacements = normalize_manifest_paths(paths)
    for path, token in ((cc, "${CC}"), (cxx, "${CXX}"), (ar, "${AR}"),
                        (ranlib, "${RANLIB}"), (linker, "${LINKER}"),
                        (build_program, "${BUILD_PROGRAM}")):
        replacements[str(path)] = token

    link_files = {
        "gf_complete_util_archive": (
            adapter / "CMakeFiles/leopard2_gf_complete_util.dir/link.txt", adapter),
        "gf_complete_archive": (
            adapter / "CMakeFiles/leopard2_gf_complete_external.dir/link.txt", adapter),
        "jerasure_archive": (
            adapter / "CMakeFiles/leopard2_jerasure_external.dir/link.txt", adapter),
        "adapter_executable": (
            adapter / "CMakeFiles/leopard2_jerasure_benchmark.dir/link.txt", adapter),
        "leopard2_library": (
            leopard / "CMakeFiles/libleopard.dir/link.txt", leopard),
        "leopard2_executable": (
            leopard / "CMakeFiles/bench_leopard2.dir/link.txt", leopard),
    }
    adapter_exe = (adapter / "leopard2_jerasure_benchmark").resolve()
    leopard_exe = (leopard / "bench_leopard2").resolve()
    static_paths = {
        "gf_complete_util_archive": adapter / "libleopard2_gf_complete_util.a",
        "gf_complete_archive": adapter / "libleopard2_gf_complete_external.a",
        "jerasure_archive": adapter / "libleopard2_jerasure_external.a",
    }
    identity: dict[str, Any] = {
        "schema": BUILD_SCHEMA,
        "source_bindings": current_source_bindings(paths),
        "recipe": {
            "adapter": {
                "generator": cmake_cache(adapter, "CMAKE_GENERATOR"),
                "build_target": "leopard2_jerasure_benchmark",
                "definitions": ["CMAKE_BUILD_TYPE=Release",
                    "CMAKE_EXPORT_COMPILE_COMMANDS=ON",
                    "LEO2_JERASURE_EXPECTED_COMMIT=" + JERASURE_COMMIT,
                    "LEO2_GF_COMPLETE_EXPECTED_COMMIT=" + GF_COMPLETE_COMMIT],
            },
            "leopard2": {
                "generator": cmake_cache(leopard, "CMAKE_GENERATOR"),
                "build_target": "bench_leopard2",
                "definitions": ["CMAKE_BUILD_TYPE=Release", "LEO2_BUILD_TESTS=OFF",
                    "LEO2_BUILD_BENCHMARKS=ON", "LEO2_BUILD_FUZZERS=OFF",
                    "LEO2_ENABLE_CUDA=OFF", "ENABLE_OPENMP=ON",
                    "CMAKE_EXPORT_COMPILE_COMMANDS=ON"],
            },
        },
        "tools": {
            "cmake": common._tool_identity(Path(shutil.which("cmake") or ""), ["--version"]),
            "cc": common._tool_identity(cc, ["--version"]),
            "cxx": common._tool_identity(cxx, ["--version"]),
            "ar": common._tool_identity(ar, ["--version"]),
            "ranlib": common._tool_identity(ranlib, ["--version"]),
            "linker": common._tool_identity(linker, ["--version"]),
            "build_program": common._tool_identity(build_program, ["--version"]),
            "ldd": common._tool_identity(ldd, ["--version"]),
            "readelf": common._tool_identity(readelf, ["--version"]),
        },
        "compile_commands": {
            "adapter": common._compile_manifest(
                adapter / "compile_commands.json", replacements),
            "leopard2": common._compile_manifest(
                leopard / "compile_commands.json", replacements),
        },
        "link_commands": {
            name: common._link_manifest(path, replacements)
            for name, (path, _) in link_files.items()},
        "link_inputs": {
            name: common._link_input_manifest(path, directory, replacements)
            for name, (path, directory) in link_files.items()},
        "static_inputs": {
            name: {"normalized_path": "${ADAPTER_BUILD}/" + path.name,
                   "sha256": sha256_file(path),
                   "required_by_adapter_link": True}
            for name, path in static_paths.items()},
        "runtime_linkage": {
            "jerasure": common._runtime_link_manifest(adapter_exe, ldd, readelf),
            "leopard2": common._runtime_link_manifest(leopard_exe, ldd, readelf),
        },
    }
    identity["identity_sha256"] = mapping_digest(identity)
    return identity


def configure_and_build(paths: Mapping[str, Path], jobs: int) -> dict[str, Any]:
    if jobs < 1 or jobs > 10:
        raise ComparisonError("bootstrap --jobs must be in 1..10")
    leopard_source = materialize_leopard(paths)
    root = paths["leopard_source"]
    for build in (paths["adapter_build"], paths["leopard_build"]):
        shutil.rmtree(build, ignore_errors=True)
    run_checked([
        "cmake", "-S", str(root / "experiments/leopard2/jerasure_compare"),
        "-B", str(paths["adapter_build"]), "-DCMAKE_BUILD_TYPE=Release",
        f"-DLEO2_JERASURE_SOURCE_ROOT={paths['jerasure_source']}",
        f"-DLEO2_GF_COMPLETE_SOURCE_ROOT={paths['gf_source']}",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"], timeout=120)
    run_checked([
        "cmake", "--build", str(paths["adapter_build"]), "--target",
        "leopard2_jerasure_benchmark", "-j", str(jobs)], timeout=600)
    run_checked([
        "cmake", "-S", str(root), "-B", str(paths["leopard_build"]),
        "-DCMAKE_BUILD_TYPE=Release", "-DLEO2_BUILD_TESTS=OFF",
        "-DLEO2_BUILD_BENCHMARKS=ON", "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_ENABLE_CUDA=OFF", "-DENABLE_OPENMP=ON",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"], timeout=120)
    run_checked([
        "cmake", "--build", str(paths["leopard_build"]), "--target",
        "bench_leopard2", "-j", str(jobs)], timeout=600)
    verify_detached_leopard(root, leopard_source)
    adapter_exe = paths["adapter_build"] / "leopard2_jerasure_benchmark"
    leopard_exe = paths["leopard_build"] / "bench_leopard2"
    if not adapter_exe.is_file() or not leopard_exe.is_file():
        raise ComparisonError("comparison executables are missing")
    return {
        "leopard_source": leopard_source,
        "build_identity": build_identity(paths),
        "executables": {
            "jerasure": {"path": str(adapter_exe.resolve()),
                         "sha256": sha256_file(adapter_exe)},
            "leopard2": {"path": str(leopard_exe.resolve()),
                         "sha256": sha256_file(leopard_exe)},
        },
    }


def portable_sources(provenance: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "jerasure": provenance["jerasure"],
        "gf_complete": provenance["gf_complete"],
        "leopard": provenance["leopard_source"],
        "build_identity": provenance["build_identity"],
    }


def bootstrap(cache: Path, jobs: int) -> dict[str, Any]:
    reject_contaminating_environment()
    paths = build_paths(cache)
    clone_pinned(paths["jerasure_source"], JERASURE_URL, JERASURE_COMMIT)
    clone_pinned(paths["gf_source"], GF_COMPLETE_URL, GF_COMPLETE_COMMIT)
    jerasure = verify_source(
        paths["jerasure_source"], "Jerasure 2.0", JERASURE_URL,
        JERASURE_COMMIT, JERASURE_TREE, JERASURE_LICENSE_SHA256,
        {"include/jerasure.h": JERASURE_HEADER_SHA256,
         "include/reed_sol.h": REED_SOL_HEADER_SHA256})
    gf_complete = verify_source(
        paths["gf_source"], "GF-Complete", GF_COMPLETE_URL,
        GF_COMPLETE_COMMIT, GF_COMPLETE_TREE, GF_COMPLETE_LICENSE_SHA256,
        {"include/gf_complete.h": GF_COMPLETE_HEADER_SHA256})
    built = configure_and_build(paths, jobs)
    provenance: dict[str, Any] = {
        "schema": PROVENANCE_SCHEMA,
        "jerasure": jerasure,
        "gf_complete": gf_complete,
        "leopard_source": built["leopard_source"],
        "build_identity": built["build_identity"],
        "executables": built["executables"],
        "constraints": {
            "production_dependency": False, "default_off": True,
            "execution_threads": 1, "region_multiple_bytes": 8,
            "application_alignment_bytes": 64,
            "maximum_bootstrap_jobs": 10,
            "maximum_original_count": MAX_K,
            "maximum_recovery_count": MAX_R,
            "maximum_matrix_coefficients": MAX_MATRIX_COEFFICIENTS,
            "maximum_application_bytes": MAX_APPLICATION_BYTES,
            "maximum_leopard_parent": MAX_PARENT,
            "gf_complete_kernel_policy": GF_COMPLETE_KERNEL_POLICY,
            "wire_compatible": False,
        },
    }
    provenance["provenance_sha256"] = mapping_digest(provenance)
    atomic_write_json(paths["provenance"], provenance)
    load_provenance(cache)
    return provenance


def validate_source_record(
        value: Mapping[str, Any], name: str, url: str, commit: str, tree: str,
        license_hash: str, headers: Mapping[str, str]) -> None:
    require_exact_keys(value, {
        "name", "url", "commit", "tree", "license", "license_sha256",
        "headers", "tracked_tree_listing_sha256", "tracked_entry_count",
        "clean_at_build"}, name)
    expected = {"name": name, "url": url, "commit": commit, "tree": tree,
                "license": "BSD-3-Clause", "license_sha256": license_hash,
                "headers": dict(headers), "clean_at_build": True}
    if any(value.get(key) != expected_value for key, expected_value in expected.items()):
        raise ComparisonError(f"{name} source identity changed")
    require_hex(value.get("tracked_tree_listing_sha256"), 64, f"{name} tree listing")
    if (isinstance(value.get("tracked_entry_count"), bool) or
            not isinstance(value.get("tracked_entry_count"), int) or
            value["tracked_entry_count"] <= 0):
        raise ComparisonError(f"{name} tracked entry count is invalid")


def validate_compile_manifest(value: Mapping[str, Any], label: str) -> None:
    require_exact_keys(value, {"schema", "entries", "entry_count", "sha256"}, label)
    if value.get("schema") != "leopard2-normalized-compile-commands/v1":
        raise ComparisonError(f"wrong {label} schema")
    entries = value.get("entries")
    if not isinstance(entries, list) or not entries:
        raise ComparisonError(f"{label} entries missing")
    prefixes = ("${LEOPARD_SOURCE}/", "${JERASURE_SOURCE}/",
                "${GF_COMPLETE_SOURCE}/")
    for entry in entries:
        if (not isinstance(entry, Mapping) or set(entry) not in (
                {"directory", "file", "command"},
                {"directory", "file", "command", "output"})):
            raise ComparisonError(f"{label} entry shape changed")
        if not any(str(entry["file"]).startswith(prefix) for prefix in prefixes):
            raise ComparisonError(f"{label} source path is not closed")
        if any(prefix in str(text) for text in entry.values()
               for prefix in ("/home/", "/Users/")):
            raise ComparisonError(f"{label} retained checkout-specific path")
    if value.get("entry_count") != len(entries):
        raise ComparisonError(f"{label} entry count changed")
    expected = hashlib.sha256(json.dumps(
        entries, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    if value.get("sha256") != expected:
        raise ComparisonError(f"{label} digest is not entry-derived")


def validate_build_identity(value: Mapping[str, Any]) -> None:
    require_exact_keys(value, {
        "schema", "source_bindings", "recipe", "tools", "compile_commands", "link_commands",
        "link_inputs", "static_inputs", "runtime_linkage", "identity_sha256"},
        "build identity")
    if value.get("schema") != BUILD_SCHEMA:
        raise ComparisonError("wrong build identity schema")
    if value.get("identity_sha256") != mapping_digest(
            {key: item for key, item in value.items() if key != "identity_sha256"}):
        raise ComparisonError("build identity digest changed")
    bindings = value.get("source_bindings", {})
    require_exact_keys(bindings, {"leopard", "jerasure", "gf_complete"},
                       "build source bindings")
    for name, binding in bindings.items():
        require_exact_keys(binding, {
            "commit", "tree", "tracked_tree_listing_sha256", "tracked_entry_count"},
            f"{name} build source binding")
        require_hex(binding.get("commit"), 40, f"{name} build commit")
        require_hex(binding.get("tree"), 40, f"{name} build tree")
        require_hex(binding.get("tracked_tree_listing_sha256"), 64,
                    f"{name} build tree listing")
        if (isinstance(binding.get("tracked_entry_count"), bool) or
                not isinstance(binding.get("tracked_entry_count"), int) or
                binding["tracked_entry_count"] <= 0):
            raise ComparisonError(f"{name} build source entry count is invalid")
    recipe = value.get("recipe", {})
    require_exact_keys(recipe, {"adapter", "leopard2"}, "build recipe")
    tools = value.get("tools", {})
    require_exact_keys(tools, {
        "cmake", "cc", "cxx", "ar", "ranlib", "linker", "build_program",
        "ldd", "readelf"}, "build tools")
    for name, tool in tools.items():
        require_exact_keys(tool, {"path", "sha256", "reported_version"}, f"tool {name}")
        require_hex(tool.get("sha256"), 64, f"tool {name} digest")
        if not isinstance(tool.get("reported_version"), str) or not tool["reported_version"]:
            raise ComparisonError(f"tool {name} version missing")
    compile_commands = value.get("compile_commands", {})
    require_exact_keys(compile_commands, {"adapter", "leopard2"}, "compile commands")
    for name, manifest in compile_commands.items():
        validate_compile_manifest(manifest, f"{name} compile commands")
    for entry in compile_commands["adapter"]["entries"]:
        if str(entry["file"]).startswith("${GF_COMPLETE_SOURCE}/"):
            flags = shlex.split(str(entry["command"]))
            if any(flag.startswith("-m") and len(flag) > 2 for flag in flags):
                raise ComparisonError(
                    "GF-Complete compile closure contains a target/ISA -m flag")
    adapter_files = [entry["file"] for entry in compile_commands["adapter"]["entries"]]
    if (not any(path.startswith("${JERASURE_SOURCE}/") for path in adapter_files) or
            not any(path.startswith("${GF_COMPLETE_SOURCE}/") for path in adapter_files) or
            not any(path.startswith("${LEOPARD_SOURCE}/experiments/leopard2/jerasure_compare/")
                    for path in adapter_files)):
        raise ComparisonError("adapter compile closure omits a source family")
    expected_links = {"gf_complete_util_archive", "gf_complete_archive",
                      "jerasure_archive", "adapter_executable",
                      "leopard2_library", "leopard2_executable"}
    links = value.get("link_commands", {})
    inputs = value.get("link_inputs", {})
    require_exact_keys(links, expected_links, "link commands")
    require_exact_keys(inputs, expected_links, "link input closures")
    for name in sorted(expected_links):
        common._validate_link_manifest(links[name], f"{name} link command")
        common._validate_link_input_manifest(inputs[name], f"{name} link inputs")
    static = value.get("static_inputs", {})
    expected_static = {"gf_complete_util_archive", "gf_complete_archive",
                       "jerasure_archive"}
    require_exact_keys(static, expected_static, "static provider inputs")
    adapter_inputs = {item["normalized_path"]: item
                      for item in inputs["adapter_executable"]["inputs"]}
    for name, record in static.items():
        require_exact_keys(record, {
            "normalized_path", "sha256", "required_by_adapter_link"},
            f"static input {name}")
        require_hex(record.get("sha256"), 64, f"static input {name}")
        if record.get("required_by_adapter_link") is not True:
            raise ComparisonError(f"static input {name} not required by adapter")
        linked = adapter_inputs.get(record.get("normalized_path"))
        if linked is None or linked.get("sha256") != record.get("sha256"):
            raise ComparisonError(f"static input {name} differs from actual adapter link")
    runtime = value.get("runtime_linkage", {})
    require_exact_keys(runtime, {"jerasure", "leopard2"}, "runtime linkage")
    for name, manifest in runtime.items():
        common._validate_runtime_link_manifest(manifest, f"{name} runtime linkage")


def validate_portable_sources(value: Mapping[str, Any]) -> None:
    require_exact_keys(value, {"jerasure", "gf_complete", "leopard", "build_identity"},
                       "portable sources")
    validate_source_record(
        value["jerasure"], "Jerasure 2.0", JERASURE_URL, JERASURE_COMMIT,
        JERASURE_TREE, JERASURE_LICENSE_SHA256,
        {"include/jerasure.h": JERASURE_HEADER_SHA256,
         "include/reed_sol.h": REED_SOL_HEADER_SHA256})
    validate_source_record(
        value["gf_complete"], "GF-Complete", GF_COMPLETE_URL, GF_COMPLETE_COMMIT,
        GF_COMPLETE_TREE, GF_COMPLETE_LICENSE_SHA256,
        {"include/gf_complete.h": GF_COMPLETE_HEADER_SHA256})
    leopard = value.get("leopard", {})
    require_exact_keys(leopard, {
        "commit", "tree", "clean_at_build", "materialization", "detached_head",
        "tracked_tree_listing_sha256", "tracked_entry_count"}, "Leopard source")
    require_hex(leopard.get("commit"), 40, "Leopard commit")
    require_hex(leopard.get("tree"), 40, "Leopard tree")
    require_hex(leopard.get("tracked_tree_listing_sha256"), 64, "Leopard tree listing")
    if (leopard.get("clean_at_build") is not True or
            leopard.get("materialization") != "detached Git worktree" or
            leopard.get("detached_head") is not True):
        raise ComparisonError("Leopard build source was not clean and detached")
    validate_build_identity(value["build_identity"])
    source_records = {
        "leopard": leopard, "jerasure": value["jerasure"],
        "gf_complete": value["gf_complete"]}
    for name, binding in value["build_identity"]["source_bindings"].items():
        if binding != {key: source_records[name][key] for key in binding}:
            raise ComparisonError(f"{name} build objects are not source-bound")


def load_provenance(cache: Path) -> dict[str, Any]:
    paths = build_paths(cache)
    try:
        document = json.loads(paths["provenance"].read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ComparisonError(f"run bootstrap first: {error}") from error
    require_exact_keys(document, {
        "schema", "jerasure", "gf_complete", "leopard_source",
        "build_identity", "executables", "constraints", "provenance_sha256"},
        "provenance")
    if document.get("schema") != PROVENANCE_SCHEMA:
        raise ComparisonError("wrong provenance schema")
    if document.get("provenance_sha256") != mapping_digest(
            {key: value for key, value in document.items() if key != "provenance_sha256"}):
        raise ComparisonError("provenance digest changed")
    validate_portable_sources(portable_sources(document))
    expected_constraints = {
        "production_dependency": False, "default_off": True,
        "execution_threads": 1, "region_multiple_bytes": 8,
        "application_alignment_bytes": 64, "maximum_bootstrap_jobs": 10,
        "maximum_original_count": MAX_K,
        "maximum_recovery_count": MAX_R,
        "maximum_matrix_coefficients": MAX_MATRIX_COEFFICIENTS,
        "maximum_application_bytes": MAX_APPLICATION_BYTES,
        "maximum_leopard_parent": MAX_PARENT,
        "gf_complete_kernel_policy": GF_COMPLETE_KERNEL_POLICY,
        "wire_compatible": False}
    if document.get("constraints") != expected_constraints:
        raise ComparisonError("provenance constraints changed")
    executables = document.get("executables", {})
    require_exact_keys(executables, {"jerasure", "leopard2"}, "executables")
    for provider, record in executables.items():
        require_exact_keys(record, {"path", "sha256"}, f"{provider} executable")
        require_hex(record.get("sha256"), 64, f"{provider} executable digest")
        path = Path(str(record.get("path", "")))
        if not path.is_file() or sha256_file(path) != record["sha256"]:
            raise ComparisonError(f"{provider} executable missing or changed")
        runtime = document["build_identity"]["runtime_linkage"][provider]
        if runtime.get("executable_sha256") != record["sha256"]:
            raise ComparisonError(f"{provider} executable differs from runtime closure")
    verify_detached_leopard(paths["leopard_source"], document["leopard_source"])
    current = build_identity(paths)
    if current != document["build_identity"]:
        raise ComparisonError("current build/link/runtime closure changed after bootstrap")
    verify_source(
        paths["jerasure_source"], "Jerasure 2.0", JERASURE_URL,
        JERASURE_COMMIT, JERASURE_TREE, JERASURE_LICENSE_SHA256,
        {"include/jerasure.h": JERASURE_HEADER_SHA256,
         "include/reed_sol.h": REED_SOL_HEADER_SHA256})
    verify_source(
        paths["gf_source"], "GF-Complete", GF_COMPLETE_URL,
        GF_COMPLETE_COMMIT, GF_COMPLETE_TREE, GF_COMPLETE_LICENSE_SHA256,
        {"include/gf_complete.h": GF_COMPLETE_HEADER_SHA256})
    return document


def ceil_pow2(value: int) -> int:
    return 1 << (value - 1).bit_length()


def validate_cell_domain(cell: Mapping[str, Any]) -> None:
    k, r = int(cell["K"]), int(cell["R"])
    losses = int(cell.get("loss_count", 0))
    batch = int(cell.get("batch", 1))
    shard_bytes = int(cell.get("shard_bytes", 8))
    if k < 1 or r < 1 or k > MAX_K or r > MAX_R:
        raise ComparisonError(f"bounded adapter requires K,R in 1..{MAX_K}")
    if losses < 0 or losses > min(k, r):
        raise ComparisonError("bounded adapter loss count is invalid")
    if batch < 1 or shard_bytes < 1 or shard_bytes > 0x7FFFFFFF:
        raise ComparisonError("bounded adapter batch/shard size is invalid")
    if shard_bytes % 8:
        raise ComparisonError("bounded adapter requires an 8-byte shard multiple")
    if k * r > MAX_MATRIX_COEFFICIENTS or (
            losses and k * k > MAX_MATRIX_COEFFICIENTS):
        raise ComparisonError("bounded adapter matrix coefficient domain exceeded")
    application_bytes = (2 * k + r + losses) * shard_bytes * batch
    if application_bytes > MAX_APPLICATION_BYTES:
        raise ComparisonError("bounded adapter application-buffer domain exceeded")


def resolved_identity(cell: Mapping[str, Any]) -> tuple[str, str, int, str]:
    validate_cell_domain(cell)
    k, r = int(cell["K"]), int(cell["R"])
    requested = str(cell["profile"])
    profile = ("low_v1" if requested in ("low", "low_v1") or
               requested == "auto" and r > k else "legacy_high_v1")
    parent = (ceil_pow2(ceil_pow2(k) + r) if profile == "low_v1"
              else ceil_pow2(k + ceil_pow2(r)))
    if parent > MAX_PARENT:
        raise ComparisonError("paired Leopard2 parent exceeds 65,536 coordinates")
    leopard_field = "gf8" if parent <= 256 else "gf16"
    jerasure_field = "gf8" if k + r <= 256 else "gf16"
    return profile, leopard_field, parent, jerasure_field


def expected_parameters(
        cell: Mapping[str, Any], iterations: int, warmup: int) -> dict[str, Any]:
    profile, _, _, _ = resolved_identity(cell)
    return {
        "K": int(cell["K"]), "R": int(cell["R"]),
        "requested_profile": profile, "requested_field": "auto",
        "requested_backend": "auto", "shard_bytes": int(cell["shard_bytes"]),
        "loss_count": int(cell["loss_count"]), "batch": int(cell["batch"]),
        "reuse": int(cell["reuse"]), "iterations": iterations,
        "warmup": warmup, "thread_count": 1, "seed": int(cell["seed"]),
    }


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


def expected_missing_indices(k: int, losses: int, seed: int) -> list[int]:
    order = list(range(k))
    random = XorShift64(seed ^ 0xD1B54A32D192ED03)
    for remaining in range(len(order), 1, -1):
        selected = random.next() % remaining
        order[remaining - 1], order[selected] = order[selected], order[remaining - 1]
    return sorted(order[:losses])


def fnv1a64_update(digest: int, values: Sequence[int]) -> int:
    for value in values:
        digest ^= value
        digest = (digest * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return digest


def expected_digests(cell: Mapping[str, Any]) -> tuple[str, str]:
    k, shard_bytes = int(cell["K"]), int(cell["shard_bytes"])
    seed, batch = int(cell["seed"]), int(cell["batch"])
    losses = expected_missing_indices(k, int(cell["loss_count"]), seed)
    original_digest = 14695981039346656037
    recovered_digest = 14695981039346656037
    for stripe in range(batch):
        random = XorShift64(seed ^ ((0x9E3779B97F4A7C15 * (stripe + 1)) &
                                    0xFFFFFFFFFFFFFFFF))
        originals = bytearray(k * shard_bytes)
        for index in range(len(originals)):
            originals[index] = random.next() >> 56
        original_digest = fnv1a64_update(original_digest, originals)
        for loss in losses:
            begin = loss * shard_bytes
            recovered_digest = fnv1a64_update(
                recovered_digest, originals[begin:begin + shard_bytes])
    return f"{original_digest:016x}", f"{recovered_digest:016x}"


def validate_summary(
        summary: Mapping[str, Any], label: str, iterations: int,
        execution: bool, input_bytes: int = 0, output_bytes: int = 0,
        input_rate: str = "", output_rate: str = "") -> None:
    try:
        common.validate_summary(
            summary, label, iterations, execution, input_bytes, output_bytes,
            input_rate, output_rate, True)
    except common.ComparisonError as error:
        raise ComparisonError(str(error)) from error


def validate_execution_summary(
        summary: Mapping[str, Any], label: str, iterations: int,
        input_bytes: int, output_bytes: int,
        input_rate: str, output_rate: str) -> None:
    suffix = "_per_batch_call"
    require_exact_keys(summary, {
        "median_us" + suffix, "mad_us" + suffix,
        "minimum_us" + suffix, "maximum_us" + suffix,
        "samples_us" + suffix, input_rate, output_rate}, label)
    try:
        values = {
            name: common.finite_positive(
                summary.get(name + suffix), f"{label}.{name + suffix}",
                allow_zero=(name == "mad_us"))
            for name in ("median_us", "mad_us", "minimum_us", "maximum_us")}
        samples_raw = summary.get("samples_us" + suffix)
        if not isinstance(samples_raw, list) or len(samples_raw) != iterations:
            raise ComparisonError(f"{label} sample cardinality changed")
        samples = [common.finite_positive(value, f"{label}.sample")
                   for value in samples_raw]
    except common.ComparisonError as error:
        raise ComparisonError(str(error)) from error
    if not (values["minimum_us"] <= values["median_us"] <= values["maximum_us"]):
        raise ComparisonError(f"{label} min/median/max order is invalid")
    derived = {
        "median_us": common.median(samples), "mad_us": common.mad(samples),
        "minimum_us": min(samples), "maximum_us": max(samples)}
    if any(not common.close_enough(values[name], expected)
           for name, expected in derived.items()):
        raise ComparisonError(f"{label} summary is not raw-sample-derived")
    for byte_count, rate_name in ((input_bytes, input_rate),
                                  (output_bytes, output_rate)):
        if byte_count == 0:
            if summary.get(rate_name) is not None:
                raise ComparisonError(f"{label}.{rate_name} must be null for zero bytes")
            continue
        expected_rate = byte_count / (values["median_us"] * 1000.0)
        try:
            observed = common.finite_positive(
                summary.get(rate_name), f"{label}.{rate_name}")
        except common.ComparisonError as error:
            raise ComparisonError(str(error)) from error
        if not common.close_enough(observed, expected_rate):
            raise ComparisonError(f"{label}.{rate_name} is not byte/time-derived")


def validate_child_result(
        document: Mapping[str, Any], cell: Mapping[str, Any],
        iterations: int, warmup: int, oracle: str = "full") -> list[int]:
    require_exact_keys(document, {
        "schema", "provider", "parameters", "comparison_identity",
        "correctness", "workload_digests", "memory", "metrics"},
        "Jerasure child")
    if document.get("schema") != CHILD_SCHEMA:
        raise ComparisonError("wrong Jerasure child schema")
    provider = document.get("provider", {})
    require_exact_keys(provider, {
        "name", "source_commit", "source_tree", "header_sha256",
        "reed_sol_header_sha256", "license", "license_sha256",
        "gf_complete_commit", "gf_complete_tree", "gf_complete_header_sha256",
        "gf_complete_license_sha256", "gf_complete_simd_flags", "field",
        "generator", "wire_compatible"}, "Jerasure provider")
    _, leopard_field, parent, jerasure_field = resolved_identity(cell)
    required_provider = {
        "name": "Jerasure 2.0", "source_commit": JERASURE_COMMIT,
        "source_tree": JERASURE_TREE, "header_sha256": JERASURE_HEADER_SHA256,
        "reed_sol_header_sha256": REED_SOL_HEADER_SHA256,
        "license": "BSD-3-Clause", "license_sha256": JERASURE_LICENSE_SHA256,
        "gf_complete_commit": GF_COMPLETE_COMMIT,
        "gf_complete_tree": GF_COMPLETE_TREE,
        "gf_complete_header_sha256": GF_COMPLETE_HEADER_SHA256,
        "gf_complete_license_sha256": GF_COMPLETE_LICENSE_SHA256,
        "field": jerasure_field,
        "generator": "reed_sol_vandermonde_coding_matrix",
        "wire_compatible": False,
    }
    if any(provider.get(key) != value for key, value in required_provider.items()):
        raise ComparisonError("Jerasure provider identity changed")
    if provider.get("gf_complete_simd_flags") != GF_COMPLETE_KERNEL_POLICY:
        raise ComparisonError("GF-Complete portable kernel policy changed")
    parameters = document.get("parameters", {})
    expected = expected_parameters(cell, iterations, warmup)
    require_exact_keys(parameters, set(expected) | {"missing_original_indices"},
                       "Jerasure parameters")
    for key, value in expected.items():
        if parameters.get(key) != value:
            raise ComparisonError(f"Jerasure parameter {key} changed")
    missing = parameters.get("missing_original_indices")
    expected_missing = expected_missing_indices(
        expected["K"], expected["loss_count"], expected["seed"])
    if missing != expected_missing:
        raise ComparisonError("Jerasure losses are not seed-derived")
    profile = expected["requested_profile"]
    identity = document.get("comparison_identity", {})
    require_exact_keys(identity, {
        "leopard2_profile", "leopard2_parent", "leopard2_field",
        "jerasure_field_advantage_from_padding", "scope"}, "comparison identity")
    scope = ("public payload and repaired-output throughput only; field/basis "
             "representation, coordinates, generator matrices, and parity bytes differ")
    if (identity.get("leopard2_profile") != profile or
            identity.get("leopard2_parent") != parent or
            identity.get("leopard2_field") != leopard_field or
            identity.get("jerasure_field_advantage_from_padding") is not
            (jerasure_field == "gf8" and leopard_field == "gf16") or
            identity.get("scope") != scope):
        raise ComparisonError("Jerasure comparison qualification changed")
    k, r = expected["K"], expected["R"]
    shard_bytes, batch = expected["shard_bytes"], expected["batch"]
    checked_bytes = r * (shard_bytes if oracle == "full" else min(shard_bytes, 64)) * batch
    total_bytes = r * shard_bytes * batch
    correctness = document.get("correctness", {})
    require_exact_keys(correctness, {
        "direct_source_round_trip", "systematic_sources_immutable",
        "independent_systematic_vandermonde_coefficients_checked",
        "independent_scalar_parity_mode",
        "independent_scalar_parity_checked_bytes_per_validation",
        "independent_scalar_parity_total_bytes_per_validation",
        "independent_scalar_parity_validation_passes"}, "Jerasure correctness")
    if (correctness.get("direct_source_round_trip") is not True or
            correctness.get("systematic_sources_immutable") is not True or
            correctness.get("independent_systematic_vandermonde_coefficients_checked") != k * r or
            correctness.get("independent_scalar_parity_mode") != oracle or
            correctness.get("independent_scalar_parity_checked_bytes_per_validation") != checked_bytes or
            correctness.get("independent_scalar_parity_total_bytes_per_validation") != total_bytes or
            correctness.get("independent_scalar_parity_validation_passes") != 2):
        raise ComparisonError("Jerasure independent correctness evidence changed")
    digests = document.get("workload_digests", {})
    require_exact_keys(digests, {
        "algorithm", "original_data", "transmitted_parity", "recovered_originals"},
        "Jerasure workload digests")
    if digests.get("algorithm") != "fnv1a64":
        raise ComparisonError("Jerasure digest algorithm changed")
    for key in ("original_data", "transmitted_parity", "recovered_originals"):
        require_hex(digests.get(key), 16, f"Jerasure {key} digest")
    expected_original, expected_recovered = expected_digests(cell)
    if (digests.get("original_data") != expected_original or
            digests.get("recovered_originals") != expected_recovered):
        raise ComparisonError("Jerasure source/recovery digest is not workload-derived")
    losses = expected["loss_count"]
    memory = document.get("memory", {})
    no_loss = losses == 0
    expected_memory = {
        "alignment_bytes": 64, "region_multiple_bytes": 8,
        "direct_application_buffers": True, "staging_copy_bytes_per_execution": 0,
        "encode_input_bytes_per_batch_call": k * shard_bytes * batch,
        "encode_output_bytes_per_batch_call": r * shard_bytes * batch,
        "decode_offered_bytes_per_batch_call": (
            0 if no_loss else (k - losses + r) * shard_bytes * batch),
        "decode_selected_bytes_per_batch_call": 0 if no_loss else k * shard_bytes * batch,
        "decode_output_bytes_per_batch_call": losses * shard_bytes * batch,
        "encode_matrix_bytes": k * r * 4,
        "decode_matrix_bytes": 0 if no_loss else k * k * 4,
        "decode_id_bytes": 0 if no_loss else k * 4,
    }
    require_exact_keys(memory, set(expected_memory), "Jerasure memory")
    if any(memory.get(key) != value for key, value in expected_memory.items()):
        raise ComparisonError("Jerasure memory/traffic semantics changed")
    metrics = document.get("metrics", {})
    require_exact_keys(metrics, {
        "codec_setup", "encode_execution", "decode_plan_setup",
        "decode_execution", "decode_amortized_at_reuse", "rate_semantics"},
        "Jerasure metrics")
    if metrics.get("rate_semantics") != (
            "offered_received counts all non-erased public shards for repair; "
            "Jerasure reads the recorded deterministic K-row subset; no-loss decode "
            "is a true no-op with null throughput"):
        raise ComparisonError("Jerasure rate semantics changed")
    validate_summary(metrics["codec_setup"], "codec setup", iterations, False)
    validate_execution_summary(
        metrics["encode_execution"], "encode", iterations,
        expected_memory["encode_input_bytes_per_batch_call"],
        expected_memory["encode_output_bytes_per_batch_call"],
        "input_GB_per_s", "parity_output_GB_per_s")
    validate_summary(metrics["decode_plan_setup"], "decode setup", iterations, False)
    validate_execution_summary(
        metrics["decode_execution"], "decode", iterations,
        expected_memory["decode_offered_bytes_per_batch_call"],
        expected_memory["decode_output_bytes_per_batch_call"],
        "offered_received_GB_per_s", "repaired_output_GB_per_s")
    amortized = metrics.get("decode_amortized_at_reuse", {})
    require_exact_keys(amortized, {
        "reuse_count", "derived_median_us_per_batch_call",
        "offered_received_GB_per_s", "repaired_output_GB_per_s"},
        "Jerasure amortized decode")
    expected_us = (float(metrics["decode_execution"]["median_us_per_batch_call"]) +
                   float(metrics["decode_plan_setup"]["median_us"]) / expected["reuse"])
    if (amortized.get("reuse_count") != expected["reuse"] or
            not math.isclose(float(amortized.get("derived_median_us_per_batch_call", -1)),
                             expected_us, rel_tol=2e-6, abs_tol=2e-6)):
        raise ComparisonError("Jerasure amortized decode is not setup-derived")
    for byte_count, name in (
            (expected_memory["decode_offered_bytes_per_batch_call"],
             "offered_received_GB_per_s"),
            (expected_memory["decode_output_bytes_per_batch_call"],
             "repaired_output_GB_per_s")):
        if byte_count == 0:
            if amortized.get(name) is not None:
                raise ComparisonError(f"Jerasure zero-byte amortized {name} must be null")
        else:
            expected_rate = byte_count / (expected_us * 1000.0)
            try:
                observed = common.finite_positive(
                    amortized.get(name), f"Jerasure amortized {name}")
            except common.ComparisonError as error:
                raise ComparisonError(str(error)) from error
            if not common.close_enough(observed, expected_rate):
                raise ComparisonError(f"Jerasure amortized {name} is not derived")
    return list(missing)


def child_command(
        provider: str, executable: Path, cell: Mapping[str, Any], iterations: int,
        warmup: int, output: Path, oracle: str | None = None) -> list[str]:
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
        command.insert(-2, "--skip-legacy")
    else:
        if oracle not in ("full", "projection"):
            raise ComparisonError("Jerasure child requires explicit oracle mode")
        command[-2:-2] = ["--oracle", oracle]
    return command


def correctness_cells(count: int) -> list[dict[str, Any]]:
    if count < 16 or count > 512:
        raise ComparisonError("--cases must be in 16..512")
    boundaries = [
        (1, 1), (1, 255), (255, 1), (2, 254), (8, 248), (16, 240),
        (64, 192), (127, 129), (128, 128), (129, 100), (191, 64),
        (192, 64), (200, 50), (224, 32), (225, 30), (240, 16),
        (129, 128), (200, 100), (256, 64), (64, 256), (160, 160),
    ]
    byte_sizes = (8, 16, 24, 32, 64, 128, 256, 1024, 4096)
    cells: list[dict[str, Any]] = []
    for index, (k, r) in enumerate(boundaries):
        maximum = min(k, r)
        loss = (0 if index % 4 == 0 else 1 if index % 4 == 1 else
                maximum if index % 4 == 2 else maximum // 2 or 1)
        cells.append({
            "K": k, "R": r, "profile": "low" if r > k else "high",
            "shard_bytes": byte_sizes[index % len(byte_sizes)],
            "loss_count": loss, "batch": 1 + index % 3,
            "reuse": 1 + index % 4, "seed": 0x4A450000 + index})
    random = XorShift64(0x4C454F324A455241)
    while len(cells) < count:
        total = 2 + random.next() % 319
        k = 1 + random.next() % (total - 1)
        r = total - k
        maximum = min(k, r)
        cells.append({
            "K": k, "R": r,
            "profile": ("auto" if len(cells) % 3 == 0 else
                        "low" if r > k else "high"),
            "shard_bytes": byte_sizes[random.next() % len(byte_sizes)],
            "loss_count": random.next() % (maximum + 1),
            "batch": 1 + random.next() % 3, "reuse": 1 + random.next() % 4,
            "seed": random.next()})
    return cells[:count]


def run_correctness(cache: Path, count: int, workers: int, output: Path) -> dict[str, Any]:
    reject_contaminating_environment()
    if workers < 1 or workers > 10:
        raise ComparisonError("correctness --workers must be in 1..10")
    provenance = load_provenance(cache)
    executable = Path(provenance["executables"]["jerasure"]["path"])
    cells = correctness_cells(count)
    environment = common.controlled_child_environment()
    with tempfile.TemporaryDirectory(prefix="leo2-jerasure-correctness-", dir=str(cache)) as temp:
        temporary = Path(temp)

        def execute(item: tuple[int, dict[str, Any]]) -> dict[str, Any]:
            index, cell = item
            result_path = temporary / f"case{index}.json"
            command = child_command(
                "jerasure", executable, cell, 3, 0, result_path, oracle="full")
            completed = run_checked(command, env=environment, timeout=300)
            if completed.stdout or completed.stderr:
                raise ComparisonError(f"correctness child {index} emitted output")
            document = json.loads(result_path.read_text())
            validate_child_result(document, cell, 3, 0, "full")
            return {"case_index": index, "cell": cell,
                    "executable_sha256": provenance["executables"]["jerasure"]["sha256"],
                    "document": document}

        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(execute, enumerate(cells)))
    results.sort(key=lambda result: result["case_index"])
    artifact: dict[str, Any] = {
        "schema": CORRECTNESS_SCHEMA,
        "sources": portable_sources(provenance),
        "executables": {"jerasure": {
            "sha256": provenance["executables"]["jerasure"]["sha256"]}},
        "child_environment": common.child_environment_record(),
        "case_generation_seed": "0x4c454f324a455241",
        "case_count": count, "worker_count": workers,
        "no_loss_cases": sum(result["cell"]["loss_count"] == 0 for result in results),
        "maximum_loss_cases": sum(
            result["cell"]["loss_count"] == min(result["cell"]["K"], result["cell"]["R"])
            for result in results),
        "gf16_cases": sum(resolved_identity(result["cell"])[3] == "gf16"
                          for result in results),
        "results": results,
    }
    artifact["artifact_sha256"] = canonical_digest(artifact)
    validate_correctness(artifact)
    atomic_write_json(output, artifact)
    return artifact


def validate_correctness(
        document: Mapping[str, Any], trusted_provenance: Mapping[str, Any] | None = None) -> None:
    require_exact_keys(document, {
        "schema", "sources", "executables", "child_environment",
        "case_generation_seed", "case_count", "worker_count", "no_loss_cases",
        "maximum_loss_cases", "gf16_cases", "results", "artifact_sha256"},
        "correctness artifact")
    if (document.get("schema") != CORRECTNESS_SCHEMA or
            document.get("case_generation_seed") != "0x4c454f324a455241"):
        raise ComparisonError("wrong correctness artifact identity")
    validate_portable_sources(document["sources"])
    try:
        common.validate_child_environment_record(document["child_environment"])
    except common.ComparisonError as error:
        raise ComparisonError(str(error)) from error
    executables = document.get("executables", {})
    require_exact_keys(executables, {"jerasure"}, "correctness executables")
    require_exact_keys(executables["jerasure"], {"sha256"}, "Jerasure executable")
    executable_hash = require_hex(
        executables["jerasure"].get("sha256"), 64, "Jerasure executable digest")
    if document["sources"]["build_identity"]["runtime_linkage"]["jerasure"].get(
            "executable_sha256") != executable_hash:
        raise ComparisonError("correctness executable differs from runtime closure")
    count = document.get("case_count")
    workers = document.get("worker_count")
    if (isinstance(count, bool) or not isinstance(count, int) or count < 16 or count > 512 or
            isinstance(workers, bool) or not isinstance(workers, int) or
            workers < 1 or workers > 10):
        raise ComparisonError("correctness campaign dimensions invalid")
    cells = correctness_cells(count)
    results = document.get("results")
    if not isinstance(results, list) or len(results) != count:
        raise ComparisonError("correctness result cardinality changed")
    for index, (cell, result) in enumerate(zip(cells, results)):
        require_exact_keys(result, {
            "case_index", "cell", "executable_sha256", "document"},
            "correctness result")
        if (result.get("case_index") != index or result.get("cell") != cell or
                result.get("executable_sha256") != executable_hash):
            raise ComparisonError("correctness result order/identity changed")
        validate_child_result(result["document"], cell, 3, 0, "full")
    derived = {
        "no_loss_cases": sum(cell["loss_count"] == 0 for cell in cells),
        "maximum_loss_cases": sum(
            cell["loss_count"] == min(cell["K"], cell["R"]) for cell in cells),
        "gf16_cases": sum(resolved_identity(cell)[3] == "gf16" for cell in cells),
    }
    if any(document.get(key) != value for key, value in derived.items()):
        raise ComparisonError("correctness coverage counts are not case-derived")
    if document.get("artifact_sha256") != canonical_digest(document):
        raise ComparisonError("correctness artifact digest changed")
    if trusted_provenance is not None:
        if document["sources"] != portable_sources(trusted_provenance):
            raise ComparisonError("correctness sources differ from trusted cache")
        if executable_hash != trusted_provenance["executables"]["jerasure"]["sha256"]:
            raise ComparisonError("correctness executable differs from trusted cache")


def validate_leopard_result(
        document: Mapping[str, Any], cell: Mapping[str, Any],
        iterations: int, warmup: int) -> list[int]:
    try:
        return common.validate_leopard_result(document, cell, iterations, warmup)
    except common.ComparisonError as error:
        raise ComparisonError(str(error)) from error


def digest_triplet(document: Mapping[str, Any]) -> Mapping[str, Any]:
    value = document.get("workload_digests")
    if not isinstance(value, Mapping):
        raise ComparisonError("provider omitted workload digests")
    return value


def aggregate_results(results: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    metric_names = ("codec_setup", "encode_execution", "decode_plan_setup",
                    "decode_execution", "decode_amortized_at_reuse")
    aggregate = []
    for cell_index, cell in enumerate(CHECKPOINT_CELLS):
        entry: dict[str, Any] = {"cell_index": cell_index, "cell": dict(cell),
                                "providers": {}}
        for provider in ("jerasure", "leopard2"):
            matching = [result for result in results
                        if result["cell_index"] == cell_index and result["provider"] == provider]
            if len(matching) != 4:
                raise ComparisonError("aggregate requires four repetitions per provider")
            summaries = {}
            for metric in metric_names:
                values = []
                for result in matching:
                    summary = result["document"]["metrics"][metric]
                    key = ("derived_median_us_per_batch_call"
                           if metric == "decode_amortized_at_reuse" else
                           "median_us_per_batch_call"
                           if metric in ("encode_execution", "decode_execution") else
                           "median_us")
                    values.append(float(summary[key]))
                summaries[metric] = {
                    "median_of_run_medians_us": float(common.median(values)),
                    "mad_of_run_medians_us": float(common.mad(values)),
                }
            entry["providers"][provider] = summaries
        entry["jerasure_over_leopard2_elapsed_ratio"] = {
            metric: (entry["providers"]["jerasure"][metric]["median_of_run_medians_us"] /
                     entry["providers"]["leopard2"][metric]["median_of_run_medians_us"])
            for metric in metric_names}
        aggregate.append(entry)
    return aggregate


def correctness_binding(document: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema": "leopard2-jerasure-correctness-binding/v1",
        "artifact_sha256": document["artifact_sha256"],
        "source_identity_sha256": mapping_digest(document["sources"]),
        "jerasure_executable_sha256": document["executables"]["jerasure"]["sha256"],
        "case_count": document["case_count"],
    }


def checkpoint_method() -> dict[str, Any]:
    return {
        "provider_order": "ABBA",
        "provider_order_by_repetition": [list(order) for order in ABBA],
        "repetitions": len(ABBA),
        "iterations_per_child": EVIDENCE_ITERATIONS,
        "warmups_per_child": EVIDENCE_WARMUP,
        "raw_samples_retained_per_metric": EVIDENCE_ITERATIONS,
        "setup_and_execution_separate": True,
        "single_thread": True,
        "region_contract_bytes": 8,
        "alignment_bytes": 64,
        "parity_wire_compatible": False,
        "affinity_rechecked_after_each_child": True,
        "durable_resume_unit": (
            "one complete validated two-provider ABBA repetition; partial units discarded"),
        "post_frequency_capture": (
            "after the second child affinity check and before child JSON validation"),
        "timing_oracle": "projection with exact FNV workload digests",
        "full_correctness_oracle": (
            "separate 128-case artifact checks every generated parity byte against "
            "independent scalar systematic-Vandermonde algebra"),
        "no_loss_decode_semantics": (
            "latency/status only; zero offered/selected/output bytes and null throughput"),
    }


def checkpoint_limitations() -> list[str]:
    return [
        "Jerasure and Leopard2 parity bytes/generator matrices differ",
        "GF8-vs-GF16 padding-boundary cells are not kernel-only comparisons",
        "single-machine bounded cells do not establish the full required matrix",
        "GF-Complete is compiled at portable compiler defaults without broad ISA flags",
        "the adapter is bounded to K,R<=4096, 2^24 dense coefficients, 8 GiB "
        "of application buffers, and a Leopard2 parent no larger than 65,536",
    ]


def validate_host_metadata(
        host: Mapping[str, Any], requested_cpu: int, reserved_cpu: int,
        original_allowed: Sequence[int]) -> None:
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
        "measurement host")
    allowed = sorted(set(int(cpu) for cpu in original_allowed))
    siblings = host.get("thread_siblings")
    if (host.get("requested_cpu") != requested_cpu or
            host.get("reserved_idle_sibling_cpu") != reserved_cpu or
            host.get("allowed_cpus") != allowed or
            host.get("allowed_cpu_count") != len(allowed) or
            requested_cpu not in allowed or reserved_cpu not in allowed or
            not isinstance(siblings, list) or requested_cpu not in siblings or
            reserved_cpu not in siblings or
            not isinstance(host.get("thread_siblings_list"), str) or
            common.parse_cpu_list(host["thread_siblings_list"]) != siblings or
            host.get("process_affinity_during_run") != [requested_cpu] or
            host.get("child_affinity_preflight") != [requested_cpu] or
            host.get("cpuinfo_processor") != requested_cpu or
            host.get("pinning") !=
            "runner sched_setaffinity singleton; all children inherit" or
            host.get("parallelism_note") !=
            "authoritative timings use one physical core; Jerasure bootstrap and "
            "correctness are capped at ten workers"):
        raise ComparisonError("measurement host affinity evidence is invalid")
    try:
        common.validate_advisory_lease(
            host.get("coordinator_lease", {}), requested_cpu, reserved_cpu)
        common.validate_frequency_snapshot(host.get("current_frequency_pre", {}), "pre")
        common.validate_frequency_snapshot(host.get("current_frequency_post", {}), "post")
    except common.ComparisonError as error:
        raise ComparisonError(str(error)) from error
    for name in ("scaling_driver", "scaling_governor",
                 "energy_performance_preference", "amd_pstate_status",
                 "boost", "no_turbo", "cpu_vendor", "cpu_model_name",
                 "cpu_family", "cpu_model", "cpu_stepping", "microcode"):
        value = host.get(name)
        if value is not None and (not isinstance(value, str) or not value):
            raise ComparisonError(f"measurement host {name} is not text or null")
    for name in ("cpuinfo_min_freq_khz", "cpuinfo_max_freq_khz",
                 "scaling_min_freq_khz", "scaling_max_freq_khz"):
        value = host.get(name)
        if (value is not None and (isinstance(value, bool) or
                                   not isinstance(value, int) or value <= 0)):
            raise ComparisonError(f"measurement host {name} is not positive kHz or null")
    for minimum, maximum in ((host.get("cpuinfo_min_freq_khz"),
                              host.get("cpuinfo_max_freq_khz")),
                             (host.get("scaling_min_freq_khz"),
                              host.get("scaling_max_freq_khz"))):
        if minimum is not None and maximum is not None and minimum > maximum:
            raise ComparisonError("measurement host frequency range is inverted")
    if (not isinstance(host.get("platform"), str) or not host["platform"] or
            not isinstance(host.get("uname"), list) or not host["uname"] or
            any(not isinstance(value, str) for value in host["uname"]) or
            not isinstance(host.get("python"), str) or not host["python"]):
        raise ComparisonError("measurement host platform identity is invalid")


def expected_run_manifest(
        provenance: Mapping[str, Any], correctness: Mapping[str, Any],
        output: Path, correctness_path: Path, cpu: int, reserved_cpu: int,
        original_allowed: Sequence[int]) -> dict[str, Any]:
    manifest: dict[str, Any] = {
        "schema": RUN_MANIFEST_SCHEMA,
        "output_path": str(output.resolve()),
        "correctness_path": str(correctness_path.resolve()),
        "cpu": cpu,
        "reserved_idle_cpu": reserved_cpu,
        "original_allowed_cpus": sorted(set(int(value) for value in original_allowed)),
        "method": checkpoint_method(),
        "cells": [dict(cell) for cell in CHECKPOINT_CELLS],
        "sources_sha256": mapping_digest(portable_sources(provenance)),
        "provenance_sha256": provenance["provenance_sha256"],
        "executables": {
            provider: record["sha256"]
            for provider, record in provenance["executables"].items()},
        "correctness_binding": correctness_binding(correctness),
        "child_environment": common.child_environment_record(cpu),
        "bounds": {
            "maximum_original_count": MAX_K,
            "maximum_recovery_count": MAX_R,
            "maximum_matrix_coefficients": MAX_MATRIX_COEFFICIENTS,
            "maximum_application_bytes": MAX_APPLICATION_BYTES,
            "maximum_leopard_parent": MAX_PARENT,
            "gf_complete_kernel_policy": GF_COMPLETE_KERNEL_POLICY,
        },
    }
    manifest["artifact_sha256"] = canonical_digest(manifest)
    return manifest


def validate_run_manifest(
        document: Mapping[str, Any], expected: Mapping[str, Any]) -> None:
    require_exact_keys(document, set(expected), "run manifest")
    if document.get("schema") != RUN_MANIFEST_SCHEMA:
        raise ComparisonError("wrong run-manifest schema")
    if document.get("artifact_sha256") != canonical_digest(document):
        raise ComparisonError("run-manifest digest changed")
    if document != expected:
        raise ComparisonError("durable run state is bound to a different exact run")


def state_directory(output: Path) -> Path:
    return output.parent / f".{output.name}.state"


def establish_run_manifest(state: Path, expected: Mapping[str, Any]) -> dict[str, Any]:
    path = state / "manifest.json"
    if path.exists():
        try:
            document = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as error:
            raise ComparisonError(f"invalid durable run manifest: {error}") from error
        validate_run_manifest(document, expected)
        return dict(document)
    if state.exists() and any(state.iterdir()):
        raise ComparisonError("durable run state has artifacts but no manifest")
    atomic_write_json(path, expected)
    return dict(expected)


def make_child_artifact(
        manifest_digest: str, cell_index: int, repetition: int, order_index: int,
        provider: str, executable_sha256: str,
        document: Mapping[str, Any]) -> dict[str, Any]:
    artifact: dict[str, Any] = {
        "schema": CHILD_ARTIFACT_SCHEMA,
        "manifest_sha256": manifest_digest,
        "cell_index": cell_index,
        "repetition": repetition,
        "order_index": order_index,
        "provider": provider,
        "executable_sha256": executable_sha256,
        "document": dict(document),
    }
    artifact["artifact_sha256"] = canonical_digest(artifact)
    return artifact


def validate_child_artifact(
        artifact: Mapping[str, Any], manifest_digest: str, cell_index: int,
        repetition: int, order_index: int, provider: str,
        executable_sha256: str, validate_payload: bool = True) -> list[int]:
    require_exact_keys(artifact, {
        "schema", "manifest_sha256", "cell_index", "repetition", "order_index",
        "provider", "executable_sha256", "document", "artifact_sha256"},
        "durable child artifact")
    if (artifact.get("schema") != CHILD_ARTIFACT_SCHEMA or
            artifact.get("manifest_sha256") != manifest_digest or
            artifact.get("cell_index") != cell_index or
            artifact.get("repetition") != repetition or
            artifact.get("order_index") != order_index or
            artifact.get("provider") != provider or
            artifact.get("executable_sha256") != executable_sha256 or
            artifact.get("artifact_sha256") != canonical_digest(artifact)):
        raise ComparisonError("durable child artifact identity changed")
    if not validate_payload:
        return []
    cell = CHECKPOINT_CELLS[cell_index]
    if provider == "jerasure":
        return validate_child_result(
            artifact["document"], cell, EVIDENCE_ITERATIONS,
            EVIDENCE_WARMUP, "projection")
    return validate_leopard_result(
        artifact["document"], cell, EVIDENCE_ITERATIONS, EVIDENCE_WARMUP)


def group_directory(state: Path, cell_index: int, repetition: int) -> Path:
    return state / "groups" / f"cell{cell_index}.rep{repetition}"


def discard_partial_group(directory: Path) -> None:
    if directory.exists():
        shutil.rmtree(directory)


def validate_group_artifacts(
        directory: Path, manifest: Mapping[str, Any], cell_index: int,
        repetition: int, original_allowed: Sequence[int],
        validate_payload: bool = True) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    order = ABBA[repetition]
    child_paths = [directory / f"order{index}.{provider}.json"
                   for index, provider in enumerate(order)]
    group_path = directory / "group.json"
    required = child_paths + [group_path]
    present = [path.exists() for path in required]
    if not all(present):
        if any(present) or (directory.exists() and any(directory.iterdir())):
            discard_partial_group(directory)
        raise FileNotFoundError("ABBA repetition is incomplete")
    try:
        children = [json.loads(path.read_text()) for path in child_paths]
        group = json.loads(group_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ComparisonError(f"durable ABBA repetition is invalid: {error}") from error
    require_exact_keys(group, {
        "schema", "manifest_sha256", "cell_index", "repetition",
        "provider_order", "child_artifact_sha256", "host", "artifact_sha256"},
        "durable ABBA repetition")
    manifest_digest = str(manifest["artifact_sha256"])
    if (group.get("schema") != GROUP_ARTIFACT_SCHEMA or
            group.get("manifest_sha256") != manifest_digest or
            group.get("cell_index") != cell_index or
            group.get("repetition") != repetition or
            group.get("provider_order") != list(order) or
            group.get("child_artifact_sha256") != [
                child.get("artifact_sha256") for child in children] or
            group.get("artifact_sha256") != canonical_digest(group)):
        raise ComparisonError("durable ABBA repetition identity changed")
    missing: dict[str, list[int]] = {}
    digests: dict[str, Mapping[str, Any]] = {}
    results: list[dict[str, Any]] = []
    for order_index, (provider, child) in enumerate(zip(order, children)):
        executable_sha = str(manifest["executables"][provider])
        missing[provider] = validate_child_artifact(
            child, manifest_digest, cell_index, repetition, order_index,
            provider, executable_sha, validate_payload)
        if validate_payload:
            digests[provider] = digest_triplet(child["document"])
        results.append({key: child[key] for key in (
            "cell_index", "repetition", "order_index", "provider",
            "executable_sha256", "document")})
    if validate_payload:
        if missing["jerasure"] != missing["leopard2"]:
            raise ComparisonError("providers used different erasure patterns")
        for key in ("algorithm", "original_data", "recovered_originals"):
            if digests["jerasure"].get(key) != digests["leopard2"].get(key):
                raise ComparisonError(f"provider workload digest mismatch: {key}")
    validate_host_metadata(
        group["host"], int(manifest["cpu"]), int(manifest["reserved_idle_cpu"]),
        original_allowed)
    return results, dict(group["host"])


def load_completed_group(
        state: Path, manifest: Mapping[str, Any], cell_index: int,
        repetition: int, original_allowed: Sequence[int],
        validate_payload: bool = True) -> tuple[list[dict[str, Any]], dict[str, Any]] | None:
    directory = group_directory(state, cell_index, repetition)
    try:
        return validate_group_artifacts(
            directory, manifest, cell_index, repetition, original_allowed,
            validate_payload)
    except FileNotFoundError:
        return None


def make_host_groups(
        cpu: int, reserved_cpu: int, original_allowed: Sequence[int],
        groups: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    return {
        "schema": HOST_GROUPS_SCHEMA,
        "requested_cpu": cpu,
        "reserved_idle_sibling_cpu": reserved_cpu,
        "allowed_cpus": sorted(set(int(value) for value in original_allowed)),
        "allowed_cpu_count": len(set(int(value) for value in original_allowed)),
        "group_count": len(groups),
        "groups": [dict(group) for group in groups],
    }


def validate_host_groups(value: Mapping[str, Any]) -> None:
    require_exact_keys(value, {
        "schema", "requested_cpu", "reserved_idle_sibling_cpu", "allowed_cpus",
        "allowed_cpu_count", "group_count", "groups"}, "checkpoint host groups")
    cpu = value.get("requested_cpu")
    reserved = value.get("reserved_idle_sibling_cpu")
    allowed = value.get("allowed_cpus")
    groups = value.get("groups")
    if (value.get("schema") != HOST_GROUPS_SCHEMA or
            isinstance(cpu, bool) or not isinstance(cpu, int) or
            isinstance(reserved, bool) or not isinstance(reserved, int) or
            cpu == reserved or not isinstance(allowed, list) or
            allowed != sorted(set(allowed)) or cpu not in allowed or reserved not in allowed or
            value.get("allowed_cpu_count") != len(allowed) or
            not isinstance(groups, list) or
            value.get("group_count") != len(CHECKPOINT_CELLS) * len(ABBA) or
            len(groups) != value.get("group_count")):
        raise ComparisonError("checkpoint host-group identity is invalid")
    seen = set()
    observed_order = []
    for group in groups:
        require_exact_keys(group, {"cell_index", "repetition", "host"},
                           "checkpoint host group")
        identity = (group.get("cell_index"), group.get("repetition"))
        if (identity in seen or identity[0] not in range(len(CHECKPOINT_CELLS)) or
                identity[1] not in range(len(ABBA))):
            raise ComparisonError("checkpoint host-group order/identity changed")
        seen.add(identity)
        observed_order.append(identity)
        validate_host_metadata(group["host"], cpu, reserved, allowed)
    expected_order = [(cell_index, repetition)
                      for cell_index in range(len(CHECKPOINT_CELLS))
                      for repetition in range(len(ABBA))]
    if observed_order != expected_order:
        raise ComparisonError("checkpoint host groups are not deterministically ordered")


def run_checkpoint(
        cache: Path, cpu: int, reserved_idle_cpu: int, output: Path,
        correctness_path: Path, iterations: int, warmup: int) -> dict[str, Any]:
    if sys.platform != "linux":
        raise ComparisonError("pinned Jerasure timing is Linux-only")
    reject_contaminating_environment()
    if iterations != EVIDENCE_ITERATIONS or warmup != EVIDENCE_WARMUP:
        raise ComparisonError("authoritative timing requires 9 iterations and 2 warmups")
    provenance = load_provenance(cache)
    correctness = json.loads(correctness_path.read_text())
    validate_correctness(correctness, trusted_provenance=provenance)
    if correctness.get("case_count") != AUTHORITATIVE_CORRECTNESS_CASES:
        raise ComparisonError("timing requires the authoritative 128-case correctness artifact")
    if output.resolve() == correctness_path.resolve():
        raise ComparisonError("timing output cannot replace correctness artifact")
    allowed_before = common.allowed_cpus()
    if cpu not in allowed_before or reserved_idle_cpu not in allowed_before or cpu == reserved_idle_cpu:
        raise ComparisonError("timing CPU and idle sibling must be distinct allowed CPUs")
    sibling_text = common.optional_text(
        Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"))
    siblings = common.parse_cpu_list(sibling_text) if sibling_text else [cpu]
    if reserved_idle_cpu not in siblings:
        raise ComparisonError("reserved idle CPU is not the timed CPU's SMT sibling")
    for cell in CHECKPOINT_CELLS:
        resolved_identity(cell)
    executables = {name: Path(record["path"])
                   for name, record in provenance["executables"].items()}
    environment = common.controlled_child_environment(cpu)
    manifest_expected = expected_run_manifest(
        provenance, correctness, output, correctness_path, cpu,
        reserved_idle_cpu, allowed_before)
    durable_state = state_directory(output)
    manifest = establish_run_manifest(durable_state, manifest_expected)
    if output.exists():
        try:
            existing = json.loads(output.read_text())
        except (OSError, json.JSONDecodeError) as error:
            raise ComparisonError(f"existing checkpoint is invalid: {error}") from error
        validate_checkpoint(existing, correctness, trusted_provenance=provenance)
        host = existing["host"]
        if (host["requested_cpu"] != cpu or
                host["reserved_idle_sibling_cpu"] != reserved_idle_cpu or
                host["allowed_cpus"] != sorted(allowed_before)):
            raise ComparisonError("existing checkpoint is bound to another CPU reservation")
        return existing
    results: list[dict[str, Any]] = []
    host_groups: list[dict[str, Any]] = []
    original_affinity = set(allowed_before)
    with common.advisory_cpu_lease(cache, cpu, reserved_idle_cpu) as lease:
        os.sched_setaffinity(0, {cpu})
        try:
            if common.allowed_cpus() != {cpu}:
                raise ComparisonError("runner failed to establish singleton affinity")
            probe = run_checked([
                sys.executable, "-c", "import json,os;print(json.dumps(sorted(os.sched_getaffinity(0))))"],
                env=environment, timeout=30)
            if probe.stderr:
                raise ComparisonError("child-affinity preflight emitted stderr")
            inherited = json.loads(probe.stdout)
            if inherited != [cpu]:
                raise ComparisonError("timing child did not inherit singleton affinity")
            for cell_index, cell in enumerate(CHECKPOINT_CELLS):
                for repetition, order in enumerate(ABBA):
                    completed_group = load_completed_group(
                        durable_state, manifest, cell_index, repetition,
                        original_affinity)
                    if completed_group is None:
                        directory = group_directory(
                            durable_state, cell_index, repetition)
                        discard_partial_group(directory)
                        directory.mkdir(parents=True, exist_ok=True)
                        pre_frequency = common.current_frequency_snapshot(cpu, "pre")
                        child_artifacts: list[dict[str, Any]] = []
                        missing_by_provider: dict[str, list[int]] = {}
                        digests_by_provider: dict[str, Mapping[str, Any]] = {}
                        post_frequency: Mapping[str, Any] | None = None
                        for order_index, provider in enumerate(order):
                            raw_path = directory / f".pending.order{order_index}.{provider}.json"
                            command = child_command(
                                provider, executables[provider], cell,
                                iterations, warmup, raw_path,
                                oracle="projection" if provider == "jerasure" else None)
                            completed = run_checked(
                                command, env=environment, timeout=900)
                            if common.allowed_cpus() != {cpu}:
                                raise ComparisonError(
                                    "runner affinity changed after a timing child")
                            if order_index == len(order) - 1:
                                post_frequency = common.current_frequency_snapshot(cpu, "post")
                            if completed.stdout or completed.stderr:
                                raise ComparisonError(
                                    f"{provider} emitted unexpected output")
                            try:
                                document = json.loads(raw_path.read_text())
                            except (OSError, json.JSONDecodeError) as error:
                                raise ComparisonError(
                                    f"invalid {provider} timing JSON: {error}") from error
                            missing = (validate_child_result(
                                document, cell, iterations, warmup, "projection")
                                if provider == "jerasure" else
                                validate_leopard_result(
                                    document, cell, iterations, warmup))
                            missing_by_provider[provider] = missing
                            digests_by_provider[provider] = digest_triplet(document)
                            child = make_child_artifact(
                                str(manifest["artifact_sha256"]), cell_index,
                                repetition, order_index, provider,
                                provenance["executables"][provider]["sha256"], document)
                            atomic_write_json(
                                directory / f"order{order_index}.{provider}.json", child)
                            child_artifacts.append(child)
                            raw_path.unlink()
                        if post_frequency is None:
                            raise ComparisonError("post-pair frequency was not captured")
                        if missing_by_provider["jerasure"] != missing_by_provider["leopard2"]:
                            raise ComparisonError("providers used different erasure patterns")
                        for key in ("algorithm", "original_data", "recovered_originals"):
                            if (digests_by_provider["jerasure"].get(key) !=
                                    digests_by_provider["leopard2"].get(key)):
                                raise ComparisonError(
                                    f"provider workload digest mismatch: {key}")
                        host = common.cpu_metadata(
                            cpu=cpu, original_allowed=original_affinity,
                            reserved_idle_cpu=reserved_idle_cpu,
                            child_affinity_preflight=inherited, lease=lease,
                            pre_frequency=pre_frequency,
                            post_frequency=post_frequency)
                        host["parallelism_note"] = (
                            "authoritative timings use one physical core; Jerasure "
                            "bootstrap and correctness are capped at ten workers")
                        validate_host_metadata(
                            host, cpu, reserved_idle_cpu, original_affinity)
                        group: dict[str, Any] = {
                            "schema": GROUP_ARTIFACT_SCHEMA,
                            "manifest_sha256": manifest["artifact_sha256"],
                            "cell_index": cell_index,
                            "repetition": repetition,
                            "provider_order": list(order),
                            "child_artifact_sha256": [
                                child["artifact_sha256"] for child in child_artifacts],
                            "host": host,
                        }
                        group["artifact_sha256"] = canonical_digest(group)
                        atomic_write_json(directory / "group.json", group)
                        for raw_path in directory.glob(".pending.*.json"):
                            raw_path.unlink()
                        completed_group = validate_group_artifacts(
                            directory, manifest, cell_index, repetition,
                            original_affinity)
                    group_results, group_host = completed_group
                    results.extend(group_results)
                    host_groups.append({
                        "cell_index": cell_index, "repetition": repetition,
                        "host": group_host})
            if common.allowed_cpus() != {cpu}:
                raise ComparisonError("runner affinity changed during timing campaign")
            post_provenance = load_provenance(cache)
            if post_provenance != provenance:
                raise ComparisonError("source/build/executable provenance changed during timing")
            post_correctness = json.loads(correctness_path.read_text())
            if post_correctness != correctness:
                raise ComparisonError("correctness artifact changed during timing")
        finally:
            os.sched_setaffinity(0, original_affinity)
    if common.allowed_cpus() != original_affinity:
        raise ComparisonError("runner failed to restore its original affinity")
    host = make_host_groups(
        cpu, reserved_idle_cpu, original_affinity, host_groups)
    validate_host_groups(host)
    checkpoint: dict[str, Any] = {
        "schema": SCHEMA, "method": checkpoint_method(),
        "sources": portable_sources(provenance),
        "executables": {name: {"sha256": record["sha256"]}
                        for name, record in provenance["executables"].items()},
        "child_environment": common.child_environment_record(cpu),
        "correctness_binding": correctness_binding(correctness),
        "host": host, "cells": [dict(cell) for cell in CHECKPOINT_CELLS],
        "results": results, "aggregate": aggregate_results(results),
        "limitations": checkpoint_limitations(),
    }
    checkpoint["artifact_sha256"] = canonical_digest(checkpoint)
    validate_checkpoint(
        checkpoint, correctness, trusted_provenance=post_provenance)
    atomic_write_json(output, checkpoint)
    return checkpoint


def validate_checkpoint(
        document: Mapping[str, Any], correctness: Mapping[str, Any],
        trusted_provenance: Mapping[str, Any] | None = None) -> None:
    validate_correctness(correctness, trusted_provenance=trusted_provenance)
    if correctness.get("case_count") != AUTHORITATIVE_CORRECTNESS_CASES:
        raise ComparisonError("checkpoint requires authoritative correctness campaign")
    require_exact_keys(document, {
        "schema", "method", "sources", "executables", "child_environment",
        "correctness_binding", "host", "cells", "results", "aggregate",
        "limitations", "artifact_sha256"}, "checkpoint")
    if document.get("schema") != SCHEMA:
        raise ComparisonError("wrong checkpoint schema")
    validate_portable_sources(document["sources"])
    if document["sources"] != correctness["sources"]:
        raise ComparisonError("checkpoint/correctness source identities differ")
    if document.get("correctness_binding") != correctness_binding(correctness):
        raise ComparisonError("checkpoint is not bound to supplied correctness artifact")
    cells = document.get("cells")
    if cells != [dict(cell) for cell in CHECKPOINT_CELLS]:
        raise ComparisonError("checkpoint cells changed")
    method = document.get("method", {})
    require_exact_keys(method, set(checkpoint_method()), "checkpoint method")
    if method != checkpoint_method():
        raise ComparisonError("checkpoint method changed")
    if document.get("limitations") != checkpoint_limitations():
        raise ComparisonError("checkpoint limitations changed")
    validate_host_groups(document.get("host", {}))
    try:
        common.validate_child_environment_record(
            document["child_environment"], int(document["host"]["requested_cpu"]))
    except (common.ComparisonError, KeyError, TypeError, ValueError) as error:
        raise ComparisonError(f"invalid timing environment/host binding: {error}") from error
    executables = document.get("executables", {})
    require_exact_keys(executables, {"jerasure", "leopard2"}, "checkpoint executables")
    for provider in executables:
        require_exact_keys(executables[provider], {"sha256"}, f"{provider} executable")
        require_hex(executables[provider].get("sha256"), 64, f"{provider} executable")
        if document["sources"]["build_identity"]["runtime_linkage"][provider].get(
                "executable_sha256") != executables[provider]["sha256"]:
            raise ComparisonError(f"{provider} executable differs from runtime closure")
    if correctness["executables"]["jerasure"]["sha256"] != executables["jerasure"]["sha256"]:
        raise ComparisonError("timing Jerasure executable differs from correctness gate")
    results = document.get("results")
    expected_count = len(CHECKPOINT_CELLS) * 8
    if not isinstance(results, list) or len(results) != expected_count:
        raise ComparisonError("checkpoint result cardinality changed")
    seen = set()
    pair_evidence: dict[tuple[int, int], dict[str, tuple[list[int], Mapping[str, Any]]]] = {}
    for result in results:
        require_exact_keys(result, {
            "cell_index", "repetition", "order_index", "provider",
            "executable_sha256", "document"}, "timing result")
        cell_index, repetition, order_index = (
            result["cell_index"], result["repetition"], result["order_index"])
        provider = result["provider"]
        if (not isinstance(cell_index, int) or cell_index < 0 or
                cell_index >= len(CHECKPOINT_CELLS) or repetition not in range(4) or
                order_index not in range(2) or ABBA[repetition][order_index] != provider):
            raise ComparisonError("timing result ABBA identity changed")
        identity = (cell_index, repetition, order_index, provider)
        if identity in seen:
            raise ComparisonError("duplicate timing result")
        seen.add(identity)
        if result["executable_sha256"] != executables[provider]["sha256"]:
            raise ComparisonError("timing result executable identity changed")
        cell = CHECKPOINT_CELLS[cell_index]
        if provider == "jerasure":
            missing = validate_child_result(
                result["document"], cell, EVIDENCE_ITERATIONS, EVIDENCE_WARMUP,
                "projection")
        else:
            missing = validate_leopard_result(
                result["document"], cell, EVIDENCE_ITERATIONS, EVIDENCE_WARMUP)
        pair_evidence.setdefault((cell_index, repetition), {})[provider] = (
            missing, digest_triplet(result["document"]))
    for identity, providers in pair_evidence.items():
        if set(providers) != {"jerasure", "leopard2"}:
            raise ComparisonError(f"timing pair {identity} is incomplete")
        if providers["jerasure"][0] != providers["leopard2"][0]:
            raise ComparisonError("paired erasure patterns differ")
        for key in ("algorithm", "original_data", "recovered_originals"):
            if providers["jerasure"][1].get(key) != providers["leopard2"][1].get(key):
                raise ComparisonError(f"paired workload digest differs: {key}")
    expected_aggregate = aggregate_results(results)
    if document.get("aggregate") != expected_aggregate:
        raise ComparisonError("checkpoint aggregate is not result-derived")
    if document.get("artifact_sha256") != canonical_digest(document):
        raise ComparisonError("checkpoint artifact digest changed")
    if trusted_provenance is not None:
        if document["sources"] != portable_sources(trusted_provenance):
            raise ComparisonError("checkpoint sources differ from trusted cache")
        for provider in executables:
            if executables[provider]["sha256"] != trusted_provenance["executables"][provider]["sha256"]:
                raise ComparisonError("checkpoint executable differs from trusted cache")


def fake_host_metadata(cpu: int, reserved_cpu: int, allowed: Sequence[int]) -> dict[str, Any]:
    return {
        "requested_cpu": cpu,
        "allowed_cpu_count": len(set(allowed)),
        "allowed_cpus": sorted(set(allowed)),
        "process_affinity_during_run": [cpu],
        "child_affinity_preflight": [cpu],
        "thread_siblings_list": f"{cpu},{reserved_cpu}",
        "thread_siblings": sorted((cpu, reserved_cpu)),
        "reserved_idle_sibling_cpu": reserved_cpu,
        "coordinator_lease": {
            "schema": common.LEASE_SCHEMA,
            "mechanism": "fcntl-flock-exclusive-advisory",
            "scope": "cooperating Leopard2 lab jobs only; not an OS-exclusive CPU reservation",
            "cpus": sorted((cpu, reserved_cpu)),
            "lock_names": [f"cpu-{value}.lock" for value in sorted((cpu, reserved_cpu))],
            "coordinator_pid": 123,
            "acquired_before_affinity": True,
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
        "cpuinfo_processor": cpu, "cpu_vendor": "test",
        "cpu_model_name": "test", "cpu_family": "1", "cpu_model": "1",
        "cpu_stepping": "1", "microcode": "1", "platform": "test",
        "uname": ["test"], "python": "test",
        "pinning": "runner sched_setaffinity singleton; all children inherit",
        "parallelism_note": (
            "authoritative timings use one physical core; Jerasure bootstrap and "
            "correctness are capped at ten workers"),
    }


def run_state_self_test() -> dict[str, int]:
    cpu, reserved, allowed = 7, 8, [7, 8]
    host = fake_host_metadata(cpu, reserved, allowed)
    validate_host_metadata(host, cpu, reserved, allowed)
    manifest: dict[str, Any] = {
        "schema": RUN_MANIFEST_SCHEMA,
        "cpu": cpu, "reserved_idle_cpu": reserved,
        "executables": {"jerasure": "1" * 64, "leopard2": "2" * 64},
    }
    manifest["artifact_sha256"] = canonical_digest(manifest)
    with tempfile.TemporaryDirectory(prefix="leo2-jerasure-state-self-test-") as temp:
        state = Path(temp) / "state"
        establish_run_manifest(state, manifest)
        validate_run_manifest(json.loads((state / "manifest.json").read_text()), manifest)
        directory = group_directory(state, 0, 0)
        directory.mkdir(parents=True)
        first = make_child_artifact(
            manifest["artifact_sha256"], 0, 0, 0, ABBA[0][0],
            manifest["executables"][ABBA[0][0]], {})
        atomic_write_json(directory / f"order0.{ABBA[0][0]}.json", first)
        if load_completed_group(
                state, manifest, 0, 0, allowed, validate_payload=False) is not None:
            raise ComparisonError("partial ABBA repetition was resumed")
        if directory.exists():
            raise ComparisonError("partial ABBA repetition was not discarded")
        directory.mkdir(parents=True)
        children = []
        for order_index, provider in enumerate(ABBA[0]):
            child = make_child_artifact(
                manifest["artifact_sha256"], 0, 0, order_index, provider,
                manifest["executables"][provider], {})
            atomic_write_json(directory / f"order{order_index}.{provider}.json", child)
            children.append(child)
        group: dict[str, Any] = {
            "schema": GROUP_ARTIFACT_SCHEMA,
            "manifest_sha256": manifest["artifact_sha256"],
            "cell_index": 0, "repetition": 0,
            "provider_order": list(ABBA[0]),
            "child_artifact_sha256": [child["artifact_sha256"] for child in children],
            "host": host,
        }
        group["artifact_sha256"] = canonical_digest(group)
        atomic_write_json(directory / "group.json", group)
        if load_completed_group(
                state, manifest, 0, 0, allowed, validate_payload=False) is None:
            raise ComparisonError("complete ABBA repetition was not resumed")
        tampered = json.loads((directory / f"order0.{ABBA[0][0]}.json").read_text())
        tampered["document"] = {"tampered": True}
        atomic_write_json(directory / f"order0.{ABBA[0][0]}.json", tampered)
        try:
            load_completed_group(
                state, manifest, 0, 0, allowed, validate_payload=False)
        except ComparisonError:
            pass
        else:
            raise ComparisonError("tampered durable child artifact was accepted")
    return {"partial_groups_discarded": 1, "complete_groups_resumed": 1,
            "state_mutations_rejected": 1}


def fake_jerasure_result(
        cell: Mapping[str, Any], iterations: int, warmup: int,
        oracle: str = "projection") -> dict[str, Any]:
    params = expected_parameters(cell, iterations, warmup)
    params["missing_original_indices"] = expected_missing_indices(
        params["K"], params["loss_count"], params["seed"])
    k, r = params["K"], params["R"]
    shard_bytes, batch = params["shard_bytes"], params["batch"]
    losses = params["loss_count"]
    profile, leopard_field, parent, jerasure_field = resolved_identity(cell)
    encode_input = k * shard_bytes * batch
    encode_output = r * shard_bytes * batch
    decode_offered = 0 if losses == 0 else (k - losses + r) * shard_bytes * batch
    decode_selected = 0 if losses == 0 else k * shard_bytes * batch
    decode_output = losses * shard_bytes * batch
    setup = common.fake_summary(False, iterations=iterations)
    encode = common.fake_summary(
        True, encode_input, encode_output, "input_GB_per_s",
        "parity_output_GB_per_s", iterations)
    plan = common.fake_summary(False, iterations=iterations)
    decode = common.fake_summary(
        True, decode_offered, decode_output, "offered_received_GB_per_s",
        "repaired_output_GB_per_s", iterations)
    if decode_offered == 0:
        decode["offered_received_GB_per_s"] = None
    amortized_us = (float(decode["median_us_per_batch_call"]) +
                    float(plan["median_us"]) / params["reuse"])
    original, recovered = expected_digests(cell)
    return {
        "schema": CHILD_SCHEMA,
        "provider": {
            "name": "Jerasure 2.0", "source_commit": JERASURE_COMMIT,
            "source_tree": JERASURE_TREE, "header_sha256": JERASURE_HEADER_SHA256,
            "reed_sol_header_sha256": REED_SOL_HEADER_SHA256,
            "license": "BSD-3-Clause", "license_sha256": JERASURE_LICENSE_SHA256,
            "gf_complete_commit": GF_COMPLETE_COMMIT,
            "gf_complete_tree": GF_COMPLETE_TREE,
            "gf_complete_header_sha256": GF_COMPLETE_HEADER_SHA256,
            "gf_complete_license_sha256": GF_COMPLETE_LICENSE_SHA256,
            "gf_complete_simd_flags": GF_COMPLETE_KERNEL_POLICY,
            "field": jerasure_field,
            "generator": "reed_sol_vandermonde_coding_matrix",
            "wire_compatible": False,
        },
        "parameters": params,
        "comparison_identity": {
            "leopard2_profile": profile, "leopard2_parent": parent,
            "leopard2_field": leopard_field,
            "jerasure_field_advantage_from_padding": (
                jerasure_field == "gf8" and leopard_field == "gf16"),
            "scope": (
                "public payload and repaired-output throughput only; field/basis "
                "representation, coordinates, generator matrices, and parity bytes differ"),
        },
        "correctness": {
            "direct_source_round_trip": True,
            "systematic_sources_immutable": True,
            "independent_systematic_vandermonde_coefficients_checked": k * r,
            "independent_scalar_parity_mode": oracle,
            "independent_scalar_parity_checked_bytes_per_validation": (
                r * (shard_bytes if oracle == "full" else min(shard_bytes, 64)) * batch),
            "independent_scalar_parity_total_bytes_per_validation": (
                r * shard_bytes * batch),
            "independent_scalar_parity_validation_passes": 2,
        },
        "workload_digests": {
            "algorithm": "fnv1a64", "original_data": original,
            "transmitted_parity": "4" * 16, "recovered_originals": recovered},
        "memory": {
            "alignment_bytes": 64, "region_multiple_bytes": 8,
            "direct_application_buffers": True,
            "staging_copy_bytes_per_execution": 0,
            "encode_input_bytes_per_batch_call": encode_input,
            "encode_output_bytes_per_batch_call": encode_output,
            "decode_offered_bytes_per_batch_call": decode_offered,
            "decode_selected_bytes_per_batch_call": decode_selected,
            "decode_output_bytes_per_batch_call": decode_output,
            "encode_matrix_bytes": k * r * 4,
            "decode_matrix_bytes": 0 if losses == 0 else k * k * 4,
            "decode_id_bytes": 0 if losses == 0 else k * 4,
        },
        "metrics": {
            "codec_setup": setup, "encode_execution": encode,
            "decode_plan_setup": plan, "decode_execution": decode,
            "decode_amortized_at_reuse": {
                "reuse_count": params["reuse"],
                "derived_median_us_per_batch_call": amortized_us,
                "offered_received_GB_per_s": (
                    decode_offered / (amortized_us * 1000.0)
                    if decode_offered else None),
                "repaired_output_GB_per_s": (
                    decode_output / (amortized_us * 1000.0)
                    if decode_output else None),
            },
            "rate_semantics": (
                "offered_received counts all non-erased public shards for repair; "
                "Jerasure reads the recorded deterministic K-row subset; no-loss "
                "decode is a true no-op with null throughput"),
        },
    }


def fake_checkpoint(correctness: Mapping[str, Any]) -> dict[str, Any]:
    jerasure_sha = correctness["executables"]["jerasure"]["sha256"]
    leopard_sha = correctness["sources"]["build_identity"]["runtime_linkage"][
        "leopard2"]["executable_sha256"]
    results = []
    hosts = []
    for cell_index, cell in enumerate(CHECKPOINT_CELLS):
        original, recovered = expected_digests(cell)
        for repetition, order in enumerate(ABBA):
            hosts.append({
                "cell_index": cell_index, "repetition": repetition,
                "host": fake_host_metadata(7, 8, [7, 8])})
            for order_index, provider in enumerate(order):
                if provider == "jerasure":
                    document = fake_jerasure_result(
                        cell, EVIDENCE_ITERATIONS, EVIDENCE_WARMUP, "projection")
                    executable_sha = jerasure_sha
                else:
                    document = common.fake_leopard_result(
                        cell, EVIDENCE_ITERATIONS, EVIDENCE_WARMUP,
                        common.LEOPARD_SCHEMA_V2)
                    document["workload_digests"].update(
                        original_data=original, recovered_originals=recovered)
                    executable_sha = leopard_sha
                results.append({
                    "cell_index": cell_index, "repetition": repetition,
                    "order_index": order_index, "provider": provider,
                    "executable_sha256": executable_sha, "document": document})
    checkpoint: dict[str, Any] = {
        "schema": SCHEMA, "method": checkpoint_method(),
        "sources": copy.deepcopy(correctness["sources"]),
        "executables": {
            "jerasure": {"sha256": jerasure_sha},
            "leopard2": {"sha256": leopard_sha}},
        "child_environment": common.child_environment_record(7),
        "correctness_binding": correctness_binding(correctness),
        "host": make_host_groups(7, 8, [7, 8], hosts),
        "cells": [dict(cell) for cell in CHECKPOINT_CELLS],
        "results": results,
        "limitations": checkpoint_limitations(),
    }
    checkpoint["aggregate"] = aggregate_results(results)
    checkpoint["artifact_sha256"] = canonical_digest(checkpoint)
    return checkpoint


def run_mutation_tests(correctness_path: Path) -> dict[str, int]:
    try:
        correctness = json.loads(correctness_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ComparisonError(f"mutation test needs correctness artifact: {error}") from error
    validate_correctness(correctness)
    checkpoint = fake_checkpoint(correctness)
    validate_checkpoint(checkpoint, correctness)
    checkpoint_mutations = []
    changed = copy.deepcopy(checkpoint)
    changed["method"]["unvalidated_claim"] = True
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["method"]["no_loss_decode_semantics"] = "claims throughput"
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["limitations"].pop()
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["host"]["groups"][0]["host"]["process_affinity_during_run"] = [7, 8]
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["host"]["groups"][0]["host"]["current_frequency_post"][
        "scaling_cur_freq_khz"] = -1
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["host"]["groups"][0]["host"]["coordinator_lease"][
        "lock_names"] = ["forged"]
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    jerasure = next(result for result in changed["results"]
                    if result["provider"] == "jerasure")
    jerasure["document"]["provider"]["gf_complete_simd_flags"] = "-march=native"
    changed["aggregate"] = aggregate_results(changed["results"])
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["results"][0]["executable_sha256"] = "f" * 64
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["aggregate"][0]["providers"]["jerasure"]["codec_setup"][
        "median_of_run_medians_us"] = 999.0
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["child_environment"]["variables"]["LD_PRELOAD"] = "/tmp/inject.so"
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    for changed in checkpoint_mutations:
        try:
            validate_checkpoint(changed, correctness)
        except ComparisonError:
            pass
        else:
            raise ComparisonError("adversarial checkpoint mutation was accepted")
    correctness_mutations = []
    changed = copy.deepcopy(correctness)
    changed["worker_count"] = 11
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    changed = copy.deepcopy(correctness)
    changed["sources"]["leopard"]["commit"] = "0" * 40
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    changed = copy.deepcopy(correctness)
    changed["child_environment"]["variables"]["LD_PRELOAD"] = "/tmp/inject.so"
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    changed = copy.deepcopy(correctness)
    changed["results"][0]["document"]["provider"][
        "gf_complete_simd_flags"] = "-mpclmul"
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    changed = copy.deepcopy(correctness)
    compile_manifest = changed["sources"]["build_identity"]["compile_commands"][
        "adapter"]
    gf_entry = next(entry for entry in compile_manifest["entries"]
                    if entry["file"].startswith("${GF_COMPLETE_SOURCE}/"))
    gf_entry["command"] += " -march=x86-64-v3"
    compile_manifest["sha256"] = hashlib.sha256(json.dumps(
        compile_manifest["entries"], sort_keys=True,
        separators=(",", ":")).encode("utf-8")).hexdigest()
    build_identity = changed["sources"]["build_identity"]
    build_identity["identity_sha256"] = mapping_digest({
        key: value for key, value in build_identity.items()
        if key != "identity_sha256"})
    changed["artifact_sha256"] = canonical_digest(changed)
    correctness_mutations.append(changed)
    changed = copy.deepcopy(correctness)
    changed["artifact_sha256"] = "0" * 64
    correctness_mutations.append(changed)
    for changed in correctness_mutations:
        try:
            validate_correctness(changed)
        except ComparisonError:
            pass
        else:
            raise ComparisonError("adversarial correctness mutation was accepted")
    no_loss = fake_jerasure_result({
        "K": 8, "R": 8, "profile": "high", "shard_bytes": 64,
        "loss_count": 0, "batch": 1, "reuse": 8, "seed": 7}, 5, 0, "full")
    validate_child_result(no_loss, {
        "K": 8, "R": 8, "profile": "high", "shard_bytes": 64,
        "loss_count": 0, "batch": 1, "reuse": 8, "seed": 7}, 5, 0, "full")
    no_loss["metrics"]["decode_execution"]["offered_received_GB_per_s"] = 1.0
    try:
        validate_child_result(no_loss, {
            "K": 8, "R": 8, "profile": "high", "shard_bytes": 64,
            "loss_count": 0, "batch": 1, "reuse": 8, "seed": 7}, 5, 0, "full")
    except ComparisonError:
        pass
    else:
        raise ComparisonError("no-loss throughput claim was accepted")
    return {"checkpoint_mutations_rejected": len(checkpoint_mutations),
            "correctness_mutations_rejected": len(correctness_mutations),
            "no_loss_mutations_rejected": 1}


def self_test() -> None:
    # Normal configuration must remain entirely independent of third parties.
    production = (repo_root() / "CMakeLists.txt").read_text()
    forbidden = (
        "add_subdirectory(experiments/leopard2/jerasure_compare",
        "leo2_jerasure_source_root", "leo2_gf_complete_source_root",
        "target_link_libraries(libleopard jerasure",
    )
    if any(pattern in production.lower() for pattern in forbidden):
        raise ComparisonError("production CMake unexpectedly references Jerasure")
    sys.path.insert(0, str(repo_root() / "tools"))
    try:
        import leopard2_external_comparison as audit
    finally:
        sys.path.pop(0)
    aligned = audit.classify("jerasure", {
        "K": 129, "R": 100, "requested_profile": "legacy_high_v1",
        "shard_bytes": 64, "thread_count": 1})
    if aligned.get("status") not in ("adapter-required", "adapter-available-unmeasured"):
        raise ComparisonError("aligned Jerasure audit classification changed")
    unaligned = audit.classify("jerasure", {
        "K": 8, "R": 4, "requested_profile": "legacy_high_v1",
        "shard_bytes": 65, "thread_count": 1})
    if unaligned.get("status") != "excluded":
        raise ComparisonError("unaligned Jerasure workload did not fail closed")
    outside = audit.classify("jerasure", {
        "K": MAX_K + 1, "R": 1, "requested_profile": "legacy_high_v1",
        "shard_bytes": 64, "loss_count": 1, "batch": 1, "thread_count": 1})
    inflated = audit.classify("jerasure", {
        "K": 40000, "R": 20000, "requested_profile": "legacy_high_v1",
        "shard_bytes": 64, "loss_count": 1, "batch": 1, "thread_count": 1})
    memory_bound = audit.classify("jerasure", {
        "K": MAX_K, "R": MAX_R, "requested_profile": "legacy_high_v1",
        "shard_bytes": 1048576, "loss_count": 1, "batch": 1,
        "thread_count": 1})
    if (outside.get("status") != "excluded" or
            inflated.get("status") != "excluded" or
            not any("parent exceeds" in reason for reason in inflated.get("reasons", [])) or
            memory_bound.get("status") != "excluded"):
        raise ComparisonError("bounded Jerasure classifications did not fail closed")
    cells = correctness_cells(32)
    if len(cells) != 32 or not any(resolved_identity(cell)[3] == "gf16" for cell in cells):
        raise ComparisonError("correctness generator lost GF16 coverage")
    for cell in cells:
        original, recovered = expected_digests(cell)
        require_hex(original, 16, "synthetic original digest")
        require_hex(recovered, 16, "synthetic recovered digest")
    state_counts = run_state_self_test()
    print(json.dumps({
        "status": "PASS", "normal_build_optional": True,
        "aligned_contract_bytes": 8, "synthetic_cells": len(cells),
        "host_metadata_binding": True,
        "bounded_domain_rejections": 3,
        **state_counts}, sort_keys=True))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    bootstrap_parser = subparsers.add_parser("bootstrap")
    bootstrap_parser.add_argument("--cache", type=Path, default=default_cache())
    bootstrap_parser.add_argument("--jobs", type=int, default=min(os.cpu_count() or 1, 10))
    correctness_parser = subparsers.add_parser("correctness")
    correctness_parser.add_argument("--cache", type=Path, default=default_cache())
    correctness_parser.add_argument("--cases", type=int,
                                    default=AUTHORITATIVE_CORRECTNESS_CASES)
    correctness_parser.add_argument("--workers", type=int,
                                    default=min(os.cpu_count() or 1, 10))
    correctness_parser.add_argument("--output", type=Path, required=True)
    validate_correctness_parser = subparsers.add_parser("validate-correctness")
    validate_correctness_parser.add_argument("artifact", type=Path)
    validate_correctness_parser.add_argument("--cache", type=Path, default=default_cache())
    validate_correctness_parser.add_argument("--require-local-build-match", action="store_true")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--cache", type=Path, default=default_cache())
    run_parser.add_argument("--cpu", type=int, required=True)
    run_parser.add_argument("--reserved-idle-cpu", type=int, required=True)
    run_parser.add_argument("--correctness-artifact", type=Path, required=True)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--iterations", type=int, default=EVIDENCE_ITERATIONS)
    run_parser.add_argument("--warmup", type=int, default=EVIDENCE_WARMUP)
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("checkpoint", type=Path)
    validate_parser.add_argument("--correctness-artifact", type=Path, required=True)
    validate_parser.add_argument("--cache", type=Path, default=default_cache())
    validate_parser.add_argument("--require-local-build-match", action="store_true")
    mutation_parser = subparsers.add_parser("mutation-test")
    mutation_parser.add_argument("correctness_artifact", type=Path)
    subparsers.add_parser("self-test")
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "bootstrap":
            document = bootstrap(arguments.cache.resolve(), arguments.jobs)
            print(json.dumps({"status": "PASS",
                              "provenance_sha256": document["provenance_sha256"]},
                             sort_keys=True))
        elif arguments.command == "correctness":
            document = run_correctness(
                arguments.cache.resolve(), arguments.cases, arguments.workers,
                arguments.output.resolve())
            print(json.dumps({"status": "PASS", "cases": document["case_count"],
                              "artifact_sha256": document["artifact_sha256"]},
                             sort_keys=True))
        elif arguments.command == "validate-correctness":
            document = json.loads(arguments.artifact.read_text())
            provenance = (load_provenance(arguments.cache.resolve())
                          if arguments.require_local_build_match else None)
            validate_correctness(document, trusted_provenance=provenance)
            print(json.dumps({"status": "PASS", "cases": document["case_count"],
                              "artifact_sha256": document["artifact_sha256"]},
                             sort_keys=True))
        elif arguments.command == "run":
            document = run_checkpoint(
                arguments.cache.resolve(), arguments.cpu, arguments.reserved_idle_cpu,
                arguments.output.resolve(), arguments.correctness_artifact.resolve(),
                arguments.iterations, arguments.warmup)
            print(json.dumps({"status": "PASS", "results": len(document["results"]),
                              "artifact_sha256": document["artifact_sha256"]},
                             sort_keys=True))
        elif arguments.command == "validate":
            checkpoint = json.loads(arguments.checkpoint.read_text())
            correctness = json.loads(arguments.correctness_artifact.read_text())
            provenance = (load_provenance(arguments.cache.resolve())
                          if arguments.require_local_build_match else None)
            validate_checkpoint(checkpoint, correctness, trusted_provenance=provenance)
            print(json.dumps({"status": "PASS", "results": len(checkpoint["results"]),
                              "artifact_sha256": checkpoint["artifact_sha256"]},
                             sort_keys=True))
        elif arguments.command == "mutation-test":
            counts = run_mutation_tests(arguments.correctness_artifact.resolve())
            print(json.dumps({"status": "PASS", **counts}, sort_keys=True))
        else:
            self_test()
    except (ComparisonError, common.ComparisonError, OSError,
            json.JSONDecodeError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
