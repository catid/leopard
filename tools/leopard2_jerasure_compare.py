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
    actual_head, actual_tree, dirty = git_identity(destination)
    symbolic = subprocess.run([
        "git", "-C", str(destination), "symbolic-ref", "-q", "HEAD"],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    actual_listing, actual_count = leopard_tree_listing(destination, "HEAD")
    if (actual_head != head or actual_tree != tree or dirty or symbolic.returncode == 0 or
            actual_listing != listing_hash or actual_count != listing_count):
        raise ComparisonError("detached Leopard source differs from benchmark commit")
    return {
        "commit": head, "tree": tree, "clean_at_build": True,
        "materialization": "detached Git worktree", "detached_head": True,
        "tracked_tree_listing_sha256": listing_hash,
        "tracked_entry_count": listing_count,
    }


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
    actual_head, actual_tree, dirty = git_identity(root)
    if (actual_head != leopard_source["commit"] or
            actual_tree != leopard_source["tree"] or dirty):
        raise ComparisonError("detached Leopard source changed while building")
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
            "wire_compatible": False,
        },
    }
    provenance["provenance_sha256"] = mapping_digest(provenance)
    paths["provenance"].parent.mkdir(parents=True, exist_ok=True)
    paths["provenance"].write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
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
        "schema", "recipe", "tools", "compile_commands", "link_commands",
        "link_inputs", "static_inputs", "runtime_linkage", "identity_sha256"},
        "build identity")
    if value.get("schema") != BUILD_SCHEMA:
        raise ComparisonError("wrong build identity schema")
    if value.get("identity_sha256") != mapping_digest(
            {key: item for key, item in value.items() if key != "identity_sha256"}):
        raise ComparisonError("build identity digest changed")
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


def resolved_identity(cell: Mapping[str, Any]) -> tuple[str, str, int, str]:
    k, r = int(cell["K"]), int(cell["R"])
    requested = str(cell["profile"])
    profile = ("low_v1" if requested in ("low", "low_v1") or
               requested == "auto" and r > k else "legacy_high_v1")
    parent = (ceil_pow2(ceil_pow2(k) + r) if profile == "low_v1"
              else ceil_pow2(k + ceil_pow2(r)))
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
    if not isinstance(provider.get("gf_complete_simd_flags"), str):
        raise ComparisonError("GF-Complete SIMD build identity is missing")
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
    expected_memory = {
        "alignment_bytes": 64, "region_multiple_bytes": 8,
        "direct_application_buffers": True, "staging_copy_bytes_per_execution": 0,
        "encode_input_bytes_per_batch_call": k * shard_bytes * batch,
        "encode_output_bytes_per_batch_call": r * shard_bytes * batch,
        "decode_offered_bytes_per_batch_call": (k - losses + r) * shard_bytes * batch,
        "decode_selected_bytes_per_batch_call": k * shard_bytes * batch,
        "decode_output_bytes_per_batch_call": losses * shard_bytes * batch,
        "encode_matrix_bytes": k * r * 4,
        "decode_matrix_bytes": k * k * 4,
        "decode_id_bytes": k * 4,
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
            "offered_received counts all non-erased public shards; Jerasure reads "
            "the recorded deterministic K-row subset"):
        raise ComparisonError("Jerasure rate semantics changed")
    validate_summary(metrics["codec_setup"], "codec setup", iterations, False)
    validate_summary(
        metrics["encode_execution"], "encode", iterations, True,
        expected_memory["encode_input_bytes_per_batch_call"],
        expected_memory["encode_output_bytes_per_batch_call"],
        "input_GB_per_s", "parity_output_GB_per_s")
    validate_summary(metrics["decode_plan_setup"], "decode setup", iterations, False)
    validate_summary(
        metrics["decode_execution"], "decode", iterations, True,
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
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
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
    executables = {name: Path(record["path"])
                   for name, record in provenance["executables"].items()}
    environment = common.controlled_child_environment(cpu)
    results: list[dict[str, Any]] = []
    pre_frequency: Mapping[str, Any]
    post_frequency: Mapping[str, Any]
    original_affinity = os.sched_getaffinity(0)
    with common.advisory_cpu_lease(cache, cpu, reserved_idle_cpu) as lease:
        os.sched_setaffinity(0, {cpu})
        try:
            probe = run_checked([
                sys.executable, "-c", "import json,os;print(json.dumps(sorted(os.sched_getaffinity(0))))"],
                env=environment, timeout=30)
            inherited = json.loads(probe.stdout)
            if inherited != [cpu]:
                raise ComparisonError("timing child did not inherit singleton affinity")
            pre_frequency = common.current_frequency_snapshot(cpu, "pre")
            with tempfile.TemporaryDirectory(prefix="leo2-jerasure-timing-", dir=str(cache)) as temp:
                temporary = Path(temp)
                for cell_index, cell in enumerate(CHECKPOINT_CELLS):
                    for repetition, order in enumerate(ABBA):
                        pair_digests: dict[str, Mapping[str, Any]] = {}
                        pair_missing: dict[str, list[int]] = {}
                        for order_index, provider in enumerate(order):
                            result_path = temporary / (
                                f"cell{cell_index}.rep{repetition}.order{order_index}.{provider}.json")
                            command = child_command(
                                provider, executables[provider], cell, iterations, warmup,
                                result_path, oracle="projection" if provider == "jerasure" else None)
                            completed = run_checked(command, env=environment, timeout=900)
                            if completed.stdout or completed.stderr:
                                raise ComparisonError(f"{provider} emitted unexpected output")
                            document = json.loads(result_path.read_text())
                            missing = (validate_child_result(
                                document, cell, iterations, warmup, "projection")
                                if provider == "jerasure" else
                                validate_leopard_result(document, cell, iterations, warmup))
                            pair_digests[provider] = digest_triplet(document)
                            pair_missing[provider] = missing
                            results.append({
                                "cell_index": cell_index, "repetition": repetition,
                                "order_index": order_index, "provider": provider,
                                "executable_sha256": provenance["executables"][provider]["sha256"],
                                "document": document})
                        if pair_missing["jerasure"] != pair_missing["leopard2"]:
                            raise ComparisonError("providers used different erasure patterns")
                        for key in ("algorithm", "original_data", "recovered_originals"):
                            if pair_digests["jerasure"].get(key) != pair_digests["leopard2"].get(key):
                                raise ComparisonError(f"provider workload digest mismatch: {key}")
            post_frequency = common.current_frequency_snapshot(cpu, "post")
            post_provenance = load_provenance(cache)
            if post_provenance != provenance:
                raise ComparisonError("source/build/executable provenance changed during timing")
            post_correctness = json.loads(correctness_path.read_text())
            if post_correctness != correctness:
                raise ComparisonError("correctness artifact changed during timing")
        finally:
            os.sched_setaffinity(0, original_affinity)
    method = {
        "provider_order": "ABBA", "provider_order_by_repetition": [list(x) for x in ABBA],
        "repetitions": 4, "iterations_per_child": iterations,
        "warmups_per_child": warmup, "retain_all_samples": True,
        "setup_and_execution_separate": True, "single_thread": True,
        "region_contract_bytes": 8, "alignment_bytes": 64,
        "parity_wire_compatible": False,
        "timing_oracle": "projection with exact FNV workload digests",
        "full_correctness_oracle": (
            "separate 128-case artifact checks every generated parity byte against "
            "independent scalar systematic-Vandermonde algebra"),
    }
    host = common.cpu_metadata(
        cpu, reserved_idle_cpu, original_affinity, inherited, lease,
        pre_frequency, post_frequency)
    checkpoint: dict[str, Any] = {
        "schema": SCHEMA, "method": method,
        "sources": portable_sources(provenance),
        "executables": {name: {"sha256": record["sha256"]}
                        for name, record in provenance["executables"].items()},
        "child_environment": common.child_environment_record(cpu),
        "correctness_binding": correctness_binding(correctness),
        "host": host, "cells": [dict(cell) for cell in CHECKPOINT_CELLS],
        "results": results, "aggregate": aggregate_results(results),
        "limitations": [
            "Jerasure and Leopard2 parity bytes/generator matrices differ",
            "GF8-vs-GF16 padding-boundary cells are not kernel-only comparisons",
            "single-machine bounded cells do not establish the full required matrix"],
    }
    checkpoint["artifact_sha256"] = canonical_digest(checkpoint)
    validate_checkpoint(checkpoint, correctness)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")
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
    if (method.get("provider_order_by_repetition") != [list(x) for x in ABBA] or
            method.get("repetitions") != 4 or
            method.get("iterations_per_child") != EVIDENCE_ITERATIONS or
            method.get("warmups_per_child") != EVIDENCE_WARMUP or
            method.get("single_thread") is not True or
            method.get("region_contract_bytes") != 8 or
            method.get("alignment_bytes") != 64 or
            method.get("parity_wire_compatible") is not False):
        raise ComparisonError("checkpoint method changed")
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
            validate_child_result(
                result["document"], cell, EVIDENCE_ITERATIONS, EVIDENCE_WARMUP,
                "projection")
        else:
            validate_leopard_result(
                result["document"], cell, EVIDENCE_ITERATIONS, EVIDENCE_WARMUP)
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
    cells = correctness_cells(32)
    if len(cells) != 32 or not any(resolved_identity(cell)[3] == "gf16" for cell in cells):
        raise ComparisonError("correctness generator lost GF16 coverage")
    for cell in cells:
        original, recovered = expected_digests(cell)
        require_hex(original, 16, "synthetic original digest")
        require_hex(recovered, 16, "synthetic recovered digest")
    print(json.dumps({
        "status": "PASS", "normal_build_optional": True,
        "aligned_contract_bytes": 8, "synthetic_cells": len(cells)}, sort_keys=True))


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
        else:
            self_test()
    except (ComparisonError, common.ComparisonError, OSError,
            json.JSONDecodeError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
