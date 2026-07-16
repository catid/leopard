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
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-isal-comparison-checkpoint/v1"
CORRECTNESS_SCHEMA = "leopard2-isal-correctness/v1"
ISA_SCHEMA = "leopard2-isal-benchmark/v1"
LEOPARD_SCHEMA = "leopard2-benchmark-v1"
ISA_URL = "https://github.com/intel/isa-l"
ISA_COMMIT = "e8cc5e87fc64b4da434f32bc1fa18184622a4998"
ISA_VERSION = "2.32.1"
NASM_VERSION = "2.16.03"
NASM_URL = (
    "https://www.nasm.us/pub/nasm/releasebuilds/2.16.03/"
    "nasm-2.16.03.tar.xz")
NASM_ARCHIVE_SHA256 = (
    "1412a1c760bbd05db026b6c0d1657affd6631cd0a63cddb6f73cc6d4aa616148")
ABBA = (("isa-l", "leopard2"), ("leopard2", "isa-l"),
        ("leopard2", "isa-l"), ("isa-l", "leopard2"))
SOURCE_BUNDLE_PATHS = (
    "CMakeLists.txt",
    "bench/leopard2/benchmark.cpp",
    "experiments/leopard2/isal_compare/CMakeLists.txt",
    "experiments/leopard2/isal_compare/NOTICE",
    "experiments/leopard2/isal_compare/isal_benchmark.cpp",
    "tools/leopard2_external_comparison.py",
    "tools/leopard2_isal_compare.py",
)

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


def run_checked(command: Sequence[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        list(command), text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False, **kwargs)
    if completed.returncode != 0:
        raise ComparisonError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}")
    return completed


def git_identity(checkout: Path) -> tuple[str, str]:
    head = run_checked(["git", "-C", str(checkout), "rev-parse", "HEAD"]).stdout.strip()
    dirty = run_checked(
        ["git", "-C", str(checkout), "status", "--porcelain",
         "--untracked-files=no"]).stdout.strip()
    return head, dirty


def full_git_status(checkout: Path) -> str:
    return run_checked([
        "git", "-C", str(checkout), "status", "--porcelain",
        "--untracked-files=normal"]).stdout.strip()


def source_bundle(root: Path) -> dict[str, Any]:
    files: dict[str, str] = {}
    for relative in SOURCE_BUNDLE_PATHS:
        path = root / relative
        if not path.is_file():
            raise ComparisonError(f"source-bundle file is missing: {relative}")
        files[relative] = sha256_file(path)
    encoded = json.dumps(files, sort_keys=True, separators=(",", ":"))
    return {
        "schema": "leopard2-isal-source-bundle/v1",
        "files": files,
        "bundle_sha256": hashlib.sha256(encoded.encode("utf-8")).hexdigest(),
    }


def validate_source_bundle(bundle: Mapping[str, Any], verify_local: bool = False) -> None:
    if bundle.get("schema") != "leopard2-isal-source-bundle/v1":
        raise ComparisonError("wrong ISA-L source-bundle schema")
    files = bundle.get("files")
    if not isinstance(files, Mapping) or set(files) != set(SOURCE_BUNDLE_PATHS):
        raise ComparisonError("ISA-L source-bundle file set changed")
    for relative, digest in files.items():
        if (not isinstance(relative, str) or not isinstance(digest, str) or
                len(digest) != 64 or
                any(character not in "0123456789abcdef" for character in digest)):
            raise ComparisonError("invalid ISA-L source-bundle entry")
    encoded = json.dumps(dict(files), sort_keys=True, separators=(",", ":"))
    expected = hashlib.sha256(encoded.encode("utf-8")).hexdigest()
    if bundle.get("bundle_sha256") != expected:
        raise ComparisonError("ISA-L source-bundle digest is not file-derived")
    if verify_local and source_bundle(repo_root()) != dict(bundle):
        raise ComparisonError("retained source bundle differs from local source files")


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
    head, dirty = git_identity(source)
    if head != ISA_COMMIT or dirty:
        raise ComparisonError(
            f"ISA-L checkout must be clean at {ISA_COMMIT}; head={head}, dirty={dirty!r}")
    build = paths["isa_build"]
    install = paths["isa_install"]
    run_checked([
        "cmake", "-S", str(source), "-B", str(build),
        "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_INSTALL_PREFIX={install}",
        f"-DCMAKE_ASM_NASM_COMPILER={nasm['path']}",
        "-DBUILD_SHARED_LIBS=OFF", "-DISAL_BUILD_TESTS=OFF",
        "-DISAL_BUILD_PERF_TESTS=OFF", "-DISAL_BUILD_ISAL_SHIM=OFF",
        "-DISAL_BUILD_IGZIP_CLI=OFF"], timeout=120)
    run_checked(["cmake", "--build", str(build), "-j", str(jobs)], timeout=900)
    run_checked(["cmake", "--install", str(build)], timeout=180)
    library = install / "lib" / "libisal.a"
    header = install / "include" / "isa-l" / "erasure_code.h"
    license_path = source / "LICENSE"
    if not library.is_file() or not header.is_file() or not license_path.is_file():
        raise ComparisonError("private ISA-L installation is incomplete")
    return {
        "name": "Intel ISA-L",
        "version": ISA_VERSION,
        "url": ISA_URL,
        "commit": head,
        "license": "BSD-3-Clause",
        "license_sha256": sha256_file(license_path),
        "library_sha256": sha256_file(library),
        "library": str(library),
        "install": str(install),
    }


def build_benchmarks(paths: Mapping[str, Path], jobs: int) -> dict[str, Any]:
    root = repo_root()
    head, tracked_dirty = git_identity(root)
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
        f"-DLEO2_ISAL_INSTALL_ROOT={isa_install}"], timeout=120)
    run_checked(["cmake", "--build", str(adapter_build), "-j", str(jobs)],
                timeout=300)

    leopard_build = paths["leopard_build"]
    run_checked([
        "cmake", "-S", str(root), "-B", str(leopard_build),
        "-DCMAKE_BUILD_TYPE=Release", "-DLEO2_BUILD_TESTS=OFF",
        "-DLEO2_BUILD_BENCHMARKS=ON", "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_ENABLE_CUDA=OFF", "-DENABLE_OPENMP=ON"], timeout=120)
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
            "commit": head, "dirty_tracked": False, "clean_at_build": True,
        },
        "source_bundle": bundle,
    }


def bootstrap(cache: Path, jobs: int) -> dict[str, Any]:
    if jobs < 1 or jobs > 8:
        raise ComparisonError("--jobs must be in 1..8")
    paths = build_paths(cache)
    cache.mkdir(parents=True, exist_ok=True)
    nasm = ensure_nasm(paths, jobs)
    isa = configure_build_install_isa(paths, nasm, jobs)
    executables = build_benchmarks(paths, jobs)
    provenance = {
        "schema": "leopard2-isal-toolchain/v1",
        "nasm": nasm,
        "isa_l": isa,
        "executables": executables,
        "constraints": {
            "third_party_location": "ignored .research cache",
            "production_dependency": False,
            "maximum_build_jobs": 8,
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
    if document.get("schema") != "leopard2-isal-toolchain/v1":
        raise ComparisonError("invalid toolchain provenance schema")
    if document.get("nasm", {}).get("archive_sha256") != NASM_ARCHIVE_SHA256:
        raise ComparisonError("toolchain provenance has the wrong NASM archive")
    if document.get("isa_l", {}).get("commit") != ISA_COMMIT:
        raise ComparisonError("toolchain provenance has the wrong ISA-L commit")
    leopard_source = document.get("executables", {}).get("leopard_source", {})
    head, tracked_dirty = git_identity(repo_root())
    if (leopard_source.get("commit") != head or
            leopard_source.get("dirty_tracked") is not False or
            leopard_source.get("clean_at_build") is not True or tracked_dirty):
        raise ComparisonError(
            "current tracked sources differ from the clean benchmark build commit")
    recorded_bundle = document.get("executables", {}).get("source_bundle", {})
    validate_source_bundle(recorded_bundle)
    if source_bundle(repo_root()) != recorded_bundle:
        raise ComparisonError("current adapter/runner/CMake sources changed after bootstrap")
    for provider in ("isa-l", "leopard2"):
        record = document.get("executables", {}).get(provider, {})
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
        iterations: int, warmup: int) -> list[int]:
    parameters = document.get("parameters")
    if not isinstance(parameters, Mapping):
        raise ComparisonError("result parameters missing")
    expected = expected_parameters(cell, iterations, warmup)
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
        iterations: int, warmup: int) -> list[int]:
    if document.get("schema") != ISA_SCHEMA:
        raise ComparisonError("wrong ISA-L result schema")
    provider = document.get("provider", {})
    required_provider = {
        "name": "Intel ISA-L", "source_commit": ISA_COMMIT,
        "license": "BSD-3-Clause", "field": "gf8",
        "generator": "gf_gen_cauchy1_matrix", "wire_compatible": False,
    }
    if any(provider.get(name) != value for name, value in required_provider.items()):
        raise ComparisonError("ISA-L provider identity/generator/wire claim changed")
    library_hash = provider.get("library_sha256")
    if (not isinstance(library_hash, str) or len(library_hash) != 64 or
            any(character not in "0123456789abcdef" for character in library_hash)):
        raise ComparisonError("ISA-L library hash is invalid")
    missing = validate_common_parameters(document, cell, iterations, warmup)
    profile, field, parent = resolved_identity(cell)
    identity = document.get("comparison_identity", {})
    expected_advantage = field == "gf16" and int(cell["K"]) + int(cell["R"]) <= 256
    if (identity.get("leopard2_profile") != profile or
            identity.get("leopard2_field") != field or
            identity.get("leopard2_parent") != parent or
            identity.get("isa_l_field_advantage_from_padding") is not expected_advantage or
            "parity bytes and generator matrices differ" not in identity.get("scope", "")):
        raise ComparisonError("ISA-L comparison identity/qualification changed")
    correctness = document.get("correctness", {})
    if (correctness.get("direct_source_round_trip") is not True or
            correctness.get("systematic_generator_prefix") is not True):
        raise ComparisonError("ISA-L direct-source correctness failed")
    memory = document.get("memory", {})
    k, r, shard_bytes = int(cell["K"]), int(cell["R"]), int(cell["shard_bytes"])
    losses, batch = int(cell["loss_count"]), int(cell["batch"])
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
    if any(memory.get(name) != value for name, value in expected_memory.items()):
        raise ComparisonError("ISA-L memory/transfer semantics changed")
    metrics = document.get("metrics", {})
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
    missing = validate_common_parameters(document, cell, iterations, warmup)
    profile, field, parent = resolved_identity(cell)
    resolved = document.get("resolved", {})
    if (resolved.get("profile") != profile or resolved.get("field") != field or
            resolved.get("parent_count") != parent or
            resolved.get("thread_count") != 1):
        raise ComparisonError("Leopard2 resolved identity changed")
    if document.get("correctness", {}).get("leopard2_round_trip") is not True:
        raise ComparisonError("Leopard2 direct-source correctness failed")
    metrics = document.get("metrics", {})
    k, r, shard_bytes = int(cell["K"]), int(cell["R"]), int(cell["shard_bytes"])
    losses, batch = int(cell["loss_count"]), int(cell["batch"])
    validate_summary(metrics.get("codec_setup", {}), "leopard codec setup", iterations,
                     False)
    validate_summary(
        metrics.get("encode_execution", {}), "leopard encode", iterations, True,
        k * shard_bytes * batch, r * shard_bytes * batch,
        "input_GB_per_s", "parity_output_GB_per_s")
    validate_summary(metrics.get("decode_plan_setup", {}), "leopard plan setup",
                     iterations, False)
    validate_summary(
        metrics.get("decode_execution", {}), "leopard decode", iterations, True,
        (k - losses + r) * shard_bytes * batch, losses * shard_bytes * batch,
        "offered_received_GB_per_s", "repaired_output_GB_per_s")
    setup = float(metrics["decode_plan_setup"]["median_us"])
    execution = float(metrics["decode_execution"]["median_us_per_batch_call"])
    expected_amortized = execution + setup / int(cell["reuse"])
    amortized = metrics.get("decode_amortized_at_reuse", {})
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


def cpu_metadata(
        cpu: int, original_allowed: set[int], reserved_idle_cpu: int,
        child_affinity_preflight: list[int]) -> dict[str, Any]:
    sibling_path = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
    governor_path = Path(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_governor")
    def optional_text(path: Path) -> str | None:
        try:
            return path.read_text().strip()
        except OSError:
            return None
    cpu_info: dict[str, str] = {}
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if not line.strip():
                break
            if ":" in line:
                name, value = line.split(":", 1)
                cpu_info[name.strip()] = value.strip()
    except OSError:
        pass
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
        "sibling_reservation_asserted_by_caller": True,
        "scaling_governor": optional_text(governor_path),
        "boost": optional_text(Path("/sys/devices/system/cpu/cpufreq/boost")),
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
        "sibling_reservation": (
            "caller explicitly reserved the recorded SMT sibling idle for the "
            "measurement window"),
        "parallelism_note": (
            "authoritative cache-sensitive measurements intentionally use one physical "
            "core; bootstrap compilation is capped at eight jobs"),
    }


def child_command(
        provider: str, executable: Path, cell: Mapping[str, Any],
        iterations: int, warmup: int, output: Path) -> list[str]:
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
                output=result_path)
            completed = run_checked(command, timeout=120)
            if completed.stdout or completed.stderr:
                raise ComparisonError(f"correctness child {index} emitted output")
            document = json.loads(result_path.read_text())
            validate_isal_result(document, cell, 3, 0)
            results.append({"case_index": index, "cell": cell, "document": document})
    artifact: dict[str, Any] = {
        "schema": CORRECTNESS_SCHEMA,
        "source_commit": ISA_COMMIT,
        "executable_sha256": provenance["executables"]["isa-l"]["sha256"],
        "source_bundle": provenance["executables"]["source_bundle"],
        "leopard_source": provenance["executables"]["leopard_source"],
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


def validate_correctness(document: Mapping[str, Any]) -> None:
    if (document.get("schema") != CORRECTNESS_SCHEMA or
            document.get("source_commit") != ISA_COMMIT):
        raise ComparisonError("wrong correctness artifact identity")
    validate_source_bundle(document.get("source_bundle", {}), verify_local=True)
    leopard_source = document.get("leopard_source", {})
    if (leopard_source.get("dirty_tracked") is not False or
            leopard_source.get("clean_at_build") is not True):
        raise ComparisonError("correctness artifact was not built from clean sources")
    count = document.get("case_count")
    if isinstance(count, bool) or not isinstance(count, int):
        raise ComparisonError("invalid correctness case count")
    expected_cells = correctness_cells(count)
    results = document.get("results")
    if not isinstance(results, list) or len(results) != count:
        raise ComparisonError("correctness result cardinality changed")
    for index, (result, cell) in enumerate(zip(results, expected_cells)):
        if (result.get("case_index") != index or result.get("cell") != cell):
            raise ComparisonError("correctness result order/cell changed")
        validate_isal_result(result.get("document", {}), cell, 3, 0)
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
    os.sched_setaffinity(0, {cpu})
    try:
        if allowed_cpus() != {cpu}:
            raise ComparisonError("runner failed to establish singleton affinity")
        preflight = run_checked([
            sys.executable, "-c",
            "import json,os;print(json.dumps(sorted(os.sched_getaffinity(0))))"],
            env=environment, timeout=30)
        if preflight.stderr:
            raise ComparisonError(f"child-affinity preflight emitted stderr: {preflight.stderr}")
        child_affinity = json.loads(preflight.stdout)
        if child_affinity != [cpu]:
            raise ComparisonError(
                f"child did not inherit singleton affinity: {child_affinity}")
        with tempfile.TemporaryDirectory(prefix="leo2-isal-", dir=str(cache)) as temporary:
            temporary_path = Path(temporary)
            for cell_index, cell in enumerate(CHECKPOINT_CELLS):
                for repetition, order in enumerate(ABBA):
                    missing_by_provider: dict[str, list[int]] = {}
                    for order_index, provider in enumerate(order):
                        result_path = temporary_path / (
                            f"cell{cell_index}.rep{repetition}.order{order_index}.{provider}.json")
                        command = child_command(
                            provider, executables[provider], cell, iterations, warmup,
                            result_path)
                        completed = run_checked(command, env=environment, timeout=600)
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
                            raise ComparisonError(f"invalid {provider} JSON: {error}") from error
                        missing = (validate_isal_result(document, cell, iterations, warmup)
                                   if provider == "isa-l" else
                                   validate_leopard_result(document, cell, iterations, warmup))
                        missing_by_provider[provider] = missing
                        results.append({
                            "cell_index": cell_index, "repetition": repetition,
                            "order_index": order_index, "provider": provider,
                            "document": document,
                        })
                    if missing_by_provider["isa-l"] != missing_by_provider["leopard2"]:
                        raise ComparisonError("providers did not use the same erasure pattern")
        host = cpu_metadata(
            cpu, original_affinity, reserved_idle_cpu, child_affinity)
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
            "wire_compatible": False,
            "comparison_scope": (
                "same public K,R, shard bytes, deterministic source/erasure pattern, "
                "batch, reuse, thread count, and offered/generated/repaired bytes"),
            "decode_input_rate_semantics": (
                "offered counts every non-erased public shard pointer; selected counts "
                "the deterministic K-row subset actually consumed"),
        },
        "sources": {
            "isa_l": {name: provenance["isa_l"][name] for name in (
                "name", "version", "url", "commit", "license",
                "license_sha256", "library_sha256")},
            "nasm": {name: provenance["nasm"][name] for name in (
                "version", "url", "archive_sha256", "executable_sha256",
                "reported_version")},
            "leopard": provenance["executables"]["leopard_source"],
            "source_bundle": provenance["executables"]["source_bundle"],
        },
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


def validate_checkpoint(document: Mapping[str, Any]) -> None:
    if document.get("schema") != SCHEMA:
        raise ComparisonError("wrong checkpoint schema")
    method = document.get("method", {})
    if (method.get("provider_order") != "ABBA" or
            method.get("provider_order_by_repetition") != [list(order) for order in ABBA] or
            method.get("repetitions") != 4 or method.get("thread_count") != 1 or
            method.get("wire_compatible") is not False or
            method.get("codec_setup_separate") is not True or
            method.get("decode_plan_setup_separate") is not True or
            method.get("execution_reuses_setup") is not True or
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
    if (sources.get("isa_l", {}).get("commit") != ISA_COMMIT or
            sources.get("isa_l", {}).get("license") != "BSD-3-Clause" or
            sources.get("nasm", {}).get("version") != NASM_VERSION or
            sources.get("nasm", {}).get("url") != NASM_URL or
            sources.get("nasm", {}).get("archive_sha256") != NASM_ARCHIVE_SHA256):
        raise ComparisonError("checkpoint source/tool identity changed")
    leopard_source = sources.get("leopard", {})
    if (leopard_source.get("dirty_tracked") is not False or
            leopard_source.get("clean_at_build") is not True):
        raise ComparisonError("checkpoint executables were not built from clean sources")
    validate_source_bundle(sources.get("source_bundle", {}), verify_local=True)
    for provider in ("isa-l", "leopard2"):
        executable = document.get("executables", {}).get(provider, {})
        digest = executable.get("sha256")
        if (not isinstance(digest, str) or len(digest) != 64 or
                any(character not in "0123456789abcdef" for character in digest) or
                "path" in executable):
            raise ComparisonError("checkpoint executable identity is not portable")
    host = document.get("host", {})
    requested_cpu = host.get("requested_cpu")
    reserved_cpu = host.get("reserved_idle_sibling_cpu")
    siblings = host.get("thread_siblings")
    if (isinstance(requested_cpu, bool) or not isinstance(requested_cpu, int) or
            isinstance(reserved_cpu, bool) or not isinstance(reserved_cpu, int) or
            requested_cpu == reserved_cpu or
            requested_cpu not in host.get("allowed_cpus", []) or
            reserved_cpu not in host.get("allowed_cpus", []) or
            not isinstance(siblings, list) or requested_cpu not in siblings or
            reserved_cpu not in siblings or
            host.get("process_affinity_during_run") != [requested_cpu] or
            host.get("child_affinity_preflight") != [requested_cpu] or
            host.get("sibling_reservation_asserted_by_caller") is not True or
            host.get("pinning") !=
            "runner sched_setaffinity singleton; all children inherit"):
        raise ComparisonError("checkpoint affinity evidence is invalid")
    results = document.get("results")
    if not isinstance(results, list) or len(results) != len(CHECKPOINT_CELLS) * 8:
        raise ComparisonError("checkpoint result cardinality changed")
    seen: set[tuple[int, int, int, str]] = set()
    for result in results:
        if not isinstance(result, Mapping):
            raise ComparisonError("checkpoint result is not an object")
        cell_index = result.get("cell_index")
        repetition = result.get("repetition")
        order_index = result.get("order_index")
        provider = result.get("provider")
        if (isinstance(cell_index, bool) or not isinstance(cell_index, int) or
                cell_index < 0 or cell_index >= len(CHECKPOINT_CELLS) or
                isinstance(repetition, bool) or not isinstance(repetition, int) or
                repetition < 0 or repetition >= 4 or order_index not in (0, 1) or
                provider not in ("isa-l", "leopard2") or
                ABBA[repetition][order_index] != provider):
            raise ComparisonError("checkpoint result schedule is invalid")
        identity = (cell_index, repetition, order_index, provider)
        if identity in seen:
            raise ComparisonError("duplicate checkpoint result")
        seen.add(identity)
        cell = CHECKPOINT_CELLS[cell_index]
        if provider == "isa-l":
            validate_isal_result(result.get("document", {}), cell, iterations, warmup)
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
                 input_rate: str = "", output_rate: str = "") -> dict[str, Any]:
    samples = [10.0, 11.0, 12.0, 13.0, 14.0]
    suffix = "_per_batch_call" if execution else ""
    result = {
        "median_us" + suffix: 12.0, "mad_us" + suffix: 1.0,
        "minimum_us" + suffix: 10.0, "maximum_us" + suffix: 14.0,
    }
    if execution:
        result[input_rate] = input_bytes / 12000.0
        result[output_rate] = output_bytes / 12000.0 if output_bytes else None
    else:
        result["samples_us"] = samples
    return result


def fake_isal_result(cell: Mapping[str, Any]) -> dict[str, Any]:
    params = expected_parameters(cell, 5, 1)
    params["missing_original_indices"] = list(range(params["loss_count"]))
    k, r, b, losses, batch = (params["K"], params["R"], params["shard_bytes"],
                              params["loss_count"], params["batch"])
    profile, field, parent = resolved_identity(cell)
    encode_input, encode_output = k * b * batch, r * b * batch
    decode_offered = (k - losses + r) * b * batch
    decode_output = losses * b * batch
    setup = fake_summary(False)
    encode = fake_summary(True, encode_input, encode_output,
                          "input_GB_per_s", "parity_output_GB_per_s")
    decode_setup = fake_summary(False)
    decode = fake_summary(True, decode_offered, decode_output,
                          "offered_received_GB_per_s", "repaired_output_GB_per_s")
    encode["samples_us_per_batch_call"] = [10.0, 11.0, 12.0, 13.0, 14.0]
    decode["samples_us_per_batch_call"] = [10.0, 11.0, 12.0, 13.0, 14.0]
    amortized = 12.0 + 12.0 / params["reuse"]
    return {
        "schema": ISA_SCHEMA,
        "provider": {"name": "Intel ISA-L", "source_commit": ISA_COMMIT,
                     "library_sha256": "1" * 64, "license": "BSD-3-Clause",
                     "field": "gf8", "generator": "gf_gen_cauchy1_matrix",
                     "wire_compatible": False},
        "parameters": params,
        "comparison_identity": {
            "leopard2_profile": profile, "leopard2_field": field,
            "leopard2_parent": parent,
            "isa_l_field_advantage_from_padding": field == "gf16",
            "scope": "public only; parity bytes and generator matrices differ"},
        "correctness": {"direct_source_round_trip": True,
                        "systematic_generator_prefix": True},
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
                            if decode_output else None)}},
    }


def fake_leopard_result(cell: Mapping[str, Any]) -> dict[str, Any]:
    params = expected_parameters(cell, 5, 1)
    params["missing_original_indices"] = list(range(params["loss_count"]))
    k, r, b, losses, batch = (params["K"], params["R"], params["shard_bytes"],
                              params["loss_count"], params["batch"])
    profile, field, parent = resolved_identity(cell)
    encode_input, encode_output = k * b * batch, r * b * batch
    decode_offered = (k - losses + r) * b * batch
    decode_output = losses * b * batch
    amortized = 12.0 + 12.0 / params["reuse"]
    return {
        "schema": LEOPARD_SCHEMA,
        "parameters": params,
        "resolved": {"profile": profile, "field": field, "parent_count": parent,
                     "thread_count": 1, "backend": "avx2", "padded_side": 1},
        "correctness": {"leopard2_round_trip": True},
        "metrics": {
            "codec_setup": fake_summary(False),
            "encode_execution": fake_summary(
                True, encode_input, encode_output,
                "input_GB_per_s", "parity_output_GB_per_s"),
            "decode_plan_setup": fake_summary(False),
            "decode_execution": fake_summary(
                True, decode_offered, decode_output,
                "offered_received_GB_per_s", "repaired_output_GB_per_s"),
            "decode_amortized_at_reuse": {
                "reuse_count": params["reuse"],
                "derived_median_us_per_batch_call": amortized,
                "offered_received_GB_per_s": decode_offered / (amortized * 1000.0),
                "repaired_output_GB_per_s": (
                    decode_output / (amortized * 1000.0)
                    if decode_output else None),
            },
        },
    }


def fake_source_bundle() -> dict[str, Any]:
    return source_bundle(repo_root())


def fake_checkpoint() -> dict[str, Any]:
    results: list[dict[str, Any]] = []
    for cell_index, cell in enumerate(CHECKPOINT_CELLS):
        for repetition, order in enumerate(ABBA):
            for order_index, provider in enumerate(order):
                document = (fake_isal_result(cell) if provider == "isa-l" else
                            fake_leopard_result(cell))
                results.append({
                    "cell_index": cell_index, "repetition": repetition,
                    "order_index": order_index, "provider": provider,
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
            "comparison_scope": "same public workload",
            "decode_input_rate_semantics": (
                "offered counts every non-erased public shard pointer; selected counts "
                "the deterministic K-row subset actually consumed"),
        },
        "sources": {
            "isa_l": {"commit": ISA_COMMIT, "license": "BSD-3-Clause"},
            "nasm": {"version": NASM_VERSION, "url": NASM_URL,
                     "archive_sha256": NASM_ARCHIVE_SHA256},
            "leopard": {"commit": "1" * 40, "dirty_tracked": False,
                        "clean_at_build": True},
            "source_bundle": fake_source_bundle(),
        },
        "executables": {
            "isa-l": {"sha256": "2" * 64},
            "leopard2": {"sha256": "3" * 64},
        },
        "host": {
            "requested_cpu": 7, "allowed_cpus": [7, 8],
            "allowed_cpu_count": 2,
            "process_affinity_during_run": [7],
            "child_affinity_preflight": [7],
            "thread_siblings": [7, 8],
            "reserved_idle_sibling_cpu": 8,
            "sibling_reservation_asserted_by_caller": True,
            "pinning": "runner sched_setaffinity singleton; all children inherit",
        },
        "cells": [dict(cell) for cell in CHECKPOINT_CELLS],
        "results": results,
        "remaining_gates": ["full matrix", "multicore", "other CPUs", "PMU"],
    }
    checkpoint["aggregate"] = aggregate_results(results)
    checkpoint["artifact_sha256"] = canonical_digest(checkpoint)
    return checkpoint


def self_test() -> None:
    for cell in CHECKPOINT_CELLS:
        validate_isal_result(fake_isal_result(cell), cell, 5, 1)
    base = fake_isal_result(CHECKPOINT_CELLS[4])
    mutations = []
    for path, value in (
            (("provider", "generator"), "gf_gen_rs_matrix"),
            (("provider", "wire_compatible"), True),
            (("provider", "source_commit"), "0" * 40),
            (("correctness", "direct_source_round_trip"), False),
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
            validate_isal_result(changed, CHECKPOINT_CELLS[4], 5, 1)
        except ComparisonError:
            pass
        else:
            raise ComparisonError("ISA-L result mutation was accepted")
    checkpoint = fake_checkpoint()
    validate_checkpoint(checkpoint)
    checkpoint_mutations = []
    changed = copy.deepcopy(checkpoint)
    changed["results"].pop()
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
    changed["sources"]["leopard"]["dirty_tracked"] = True
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["sources"]["source_bundle"]["files"][SOURCE_BUNDLE_PATHS[0]] = "0" * 64
    files = changed["sources"]["source_bundle"]["files"]
    encoded = json.dumps(files, sort_keys=True, separators=(",", ":"))
    changed["sources"]["source_bundle"]["bundle_sha256"] = hashlib.sha256(
        encoded.encode("utf-8")).hexdigest()
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    changed = copy.deepcopy(checkpoint)
    changed["host"]["child_affinity_preflight"] = [7, 8]
    changed["artifact_sha256"] = canonical_digest(changed)
    checkpoint_mutations.append(changed)
    for changed in checkpoint_mutations:
        try:
            validate_checkpoint(changed)
        except ComparisonError:
            pass
        else:
            raise ComparisonError("checkpoint mutation was accepted")
    print(json.dumps({
        "cells": len(CHECKPOINT_CELLS),
        "result_mutations_rejected": len(mutations),
        "checkpoint_mutations_rejected": len(checkpoint_mutations),
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
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("checkpoint", type=Path)
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
            validate_correctness(document)
            print(json.dumps({
                "artifact_sha256": document["artifact_sha256"],
                "cases": document["case_count"], "status": "PASS",
            }, sort_keys=True))
        elif arguments.command == "validate":
            document = json.loads(arguments.checkpoint.read_text())
            validate_checkpoint(document)
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
