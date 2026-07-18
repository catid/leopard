#!/usr/bin/env python3
"""Run resumable exact-mask versus mature-prefix encoder experiments.

The C++ cell benchmark contains both paths in one binary and interleaves their
samples.  This runner adds deterministic manifests, source/binary identity,
resumable per-cell artifacts, CPU pinning, and fail-closed authority rules.
It never configures or builds Leopard: a pinned measurement must consume a
clean, already-built executable whose embedded Git SHA matches the source tree.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


SCHEMA = "leopard2-sparse-encode-crossover/v1"
JOB_SCHEMA = "leopard2-sparse-encode-crossover-job/v1"
BENCHMARK_SCHEMA = "leopard2-sparse-encode-benchmark-v2"
ATTESTATION_SCHEMA = "leopard2-benchmark-isolation-attestation/v1"
SOURCE_FILES = (
    "CMakeLists.txt",
    "LeopardCommon.cpp",
    "LeopardCommon.h",
    "Leopard2Backend.cpp",
    "Leopard2Backend.h",
    "Leopard2BackendAVX2.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendScalar.cpp",
    "Leopard2CpuFeatures.cpp",
    "Leopard2Direct.h",
    "Leopard2Dispatch.h",
    "Leopard2Plan.cpp",
    "Leopard2Plan.h",
    "LeopardFF8.cpp",
    "LeopardFF8.h",
    "LeopardFF16.cpp",
    "LeopardFF16.h",
    "leopard.cpp",
    "leopard.h",
    "leopard2.cpp",
    "leopard2.h",
    "bench/leopard2/sparse_encode_benchmark.cpp",
    "tools/leopard2_sparse_encode_benchmark_json_test.py",
    "tools/leopard2_sparse_encode_crossover.py",
)
KNOWN_BACKENDS = ("scalar", "ssse3", "avx2", "auto")
CELL_SCHEMA = "leopard2-sparse-encode-cells/v1"


class CrossoverError(RuntimeError):
    pass


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest_value(value: Any) -> str:
    return digest_bytes(canonical_bytes(value))


def atomic_write(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(value)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write(
        path,
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True).encode("utf-8")
        + b"\n",
    )


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CrossoverError(f"cannot read {path}: {error}") from error


def parse_csv_unsigned(text: str, name: str, maximum: int) -> list[int]:
    result: list[int] = []
    for raw in text.split(","):
        item = raw.strip()
        if not item or not item.isdigit():
            raise CrossoverError(f"{name} must contain unsigned decimal integers")
        value = int(item)
        if value <= 0 or value > maximum:
            raise CrossoverError(f"{name} values must be in 1..{maximum}")
        if value in result:
            raise CrossoverError(f"{name} contains duplicate {value}")
        result.append(value)
    if not result:
        raise CrossoverError(f"{name} is empty")
    return sorted(result)


def parse_backends(text: str) -> list[str]:
    result = [item.strip().lower() for item in text.split(",") if item.strip()]
    if not result or len(result) != len(set(result)):
        raise CrossoverError("backends must be a nonempty unique list")
    invalid = sorted(set(result) - set(KNOWN_BACKENDS))
    if invalid:
        raise CrossoverError(f"unsupported backends: {','.join(invalid)}")
    return result


def allowed_cpus() -> list[int]:
    if hasattr(os, "sched_getaffinity"):
        try:
            cpus = sorted(os.sched_getaffinity(0))
            if cpus:
                return cpus
        except OSError:
            pass
    return list(range(os.cpu_count() or 1))


def compact_cpu_list(cpus: list[int]) -> str:
    values = sorted(set(cpus))
    if not values:
        return ""
    ranges: list[str] = []
    first = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(first) if first == previous else f"{first}-{previous}")
        first = previous = value
    ranges.append(str(first) if first == previous else f"{first}-{previous}")
    return ",".join(ranges)


def read_optional(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError):
        return None


def machine_identity(cpus: list[int]) -> dict[str, Any]:
    uname = platform.uname()
    model = None
    cpuinfo = read_optional(Path("/proc/cpuinfo"))
    if cpuinfo:
        for line in cpuinfo.splitlines():
            if line.lower().startswith("model name") and ":" in line:
                model = line.split(":", 1)[1].strip()
                break
    return {
        "allowed_cpu_list": compact_cpu_list(cpus),
        "architecture": platform.machine(),
        "cpu_model": model,
        "logical_cpus_allowed": len(cpus),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "uname": {
            "machine": uname.machine,
            "node": uname.node,
            "release": uname.release,
            "system": uname.system,
            "version": uname.version,
        },
    }


def build_metadata(executable: Path) -> dict[str, Any]:
    result: dict[str, Any] = {"executable_sha256": digest_bytes(executable.read_bytes())}
    cache = executable.parent / "CMakeCache.txt"
    if not cache.is_file():
        result["cmake_cache"] = None
        return result
    data = cache.read_bytes()
    fields: dict[str, str] = {}
    for raw in data.decode("utf-8", errors="replace").splitlines():
        if raw.startswith("//") or raw.startswith("#") or "=" not in raw:
            continue
        key_type, value = raw.split("=", 1)
        key = key_type.split(":", 1)[0]
        if key in (
            "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS",
            "CMAKE_CXX_FLAGS_RELEASE", "LEO2_BACKEND_VARIANT",
            "LEOPARD_ENABLE_GF8", "LEOPARD_ENABLE_GF16",
        ):
            fields[key] = value
    result["cmake_cache"] = {
        "path": str(cache), "sha256": digest_bytes(data), "selected_fields": fields,
    }
    return result


def source_fingerprint(source: Path) -> dict[str, Any]:
    files: dict[str, str] = {}
    missing: list[str] = []
    for relative in SOURCE_FILES:
        path = source / relative
        if path.is_file():
            files[relative] = digest_bytes(path.read_bytes())
        else:
            missing.append(relative)
    if missing:
        raise CrossoverError(f"required source files are missing: {','.join(missing)}")
    return {"digest": digest_value(files), "files": files}


def git_state(source: Path) -> dict[str, Any]:
    def git(*arguments: str) -> str:
        completed = subprocess.run(
            ["git", *arguments], cwd=source, check=False,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=30,
        )
        if completed.returncode != 0:
            raise CrossoverError(
                f"git {' '.join(arguments)} failed: {completed.stderr.strip()}"
            )
        return completed.stdout.strip()

    sha = git("rev-parse", "HEAD")
    dirty_text = git("status", "--porcelain", "--untracked-files=normal")
    if len(sha) != 40 or any(character not in "0123456789abcdef" for character in sha):
        raise CrossoverError("Git HEAD is not a full lowercase SHA-1")
    return {"dirty": bool(dirty_text), "sha": sha}


def base_cells() -> list[dict[str, Any]]:
    return [
        {"profile": "high", "field": "gf8", "K": 192, "R": 64,
         "mask_name": "edge_sparse", "requested_parity": [0, 31, 63]},
        {"profile": "high", "field": "gf8", "K": 192, "R": 64,
         "mask_name": "scattered_sparse", "requested_parity": [3, 15, 39, 63]},
        {"profile": "low", "field": "gf8", "K": 32, "R": 192,
         "mask_name": "edge_sparse", "requested_parity": [0, 31, 32, 95, 191]},
        {"profile": "low", "field": "gf8", "K": 32, "R": 192,
         "mask_name": "scattered_sparse", "requested_parity": [3, 23, 47, 111, 191]},
        {"profile": "high", "field": "gf16", "K": 1000, "R": 200,
         "mask_name": "edge_sparse", "requested_parity": [0, 63, 127, 199]},
        {"profile": "high", "field": "gf16", "K": 1000, "R": 200,
         "mask_name": "scattered_sparse", "requested_parity": [3, 31, 95, 159, 199]},
        {"profile": "low", "field": "gf16", "K": 128, "R": 896,
         "mask_name": "edge_sparse", "requested_parity": [0, 127, 128, 383, 895]},
        {"profile": "low", "field": "gf16", "K": 128, "R": 896,
         "mask_name": "scattered_sparse", "requested_parity": [7, 63, 135, 255, 519, 895]},
    ]


def make_cells(backends: list[str], shard_bytes: list[int]) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for base in base_cells():
        for backend in sorted(backends):
            for byte_count in sorted(shard_bytes):
                cell = dict(base)
                cell.update(backend=backend, shard_bytes=byte_count)
                cells.append(cell)
    return sorted(
        cells,
        key=lambda cell: (
            cell["field"], cell["profile"], cell["K"], cell["R"],
            cell["mask_name"], cell["backend"], cell["shard_bytes"],
        ),
    )


def validate_cells(value: Any) -> list[dict[str, Any]]:
    if not isinstance(value, dict) or value.get("schema") != CELL_SCHEMA:
        raise CrossoverError("cell manifest has an unknown schema")
    raw_cells = value.get("cells")
    if not isinstance(raw_cells, list) or not raw_cells:
        raise CrossoverError("cell manifest contains no cells")
    cells: list[dict[str, Any]] = []
    for index, raw in enumerate(raw_cells):
        if not isinstance(raw, dict):
            raise CrossoverError(f"cell {index} is not an object")
        expected_keys = {
            "profile", "field", "K", "R", "mask_name",
            "requested_parity", "backend", "shard_bytes",
        }
        if set(raw) != expected_keys:
            raise CrossoverError(f"cell {index} has unknown or missing fields")
        profile = raw["profile"]
        field = raw["field"]
        backend = raw["backend"]
        k = raw["K"]
        r = raw["R"]
        byte_count = raw["shard_bytes"]
        requested_raw = raw["requested_parity"]
        if profile not in ("high", "low") or field not in ("gf8", "gf16"):
            raise CrossoverError(f"cell {index} has an invalid profile or field")
        if backend not in KNOWN_BACKENDS:
            raise CrossoverError(f"cell {index} has an invalid backend")
        if not isinstance(k, int) or isinstance(k, bool) or k <= 0 or (
            not isinstance(r, int) or isinstance(r, bool) or r <= 0
        ):
            raise CrossoverError(f"cell {index} has invalid K or R")
        if not isinstance(byte_count, int) or isinstance(byte_count, bool) or (
            byte_count <= 0 or byte_count > 1 << 30
        ):
            raise CrossoverError(f"cell {index} has invalid shard bytes")
        if field == "gf16" and byte_count % 2 != 0:
            raise CrossoverError(f"cell {index} has an odd GF16 byte count")
        if not isinstance(raw["mask_name"], str) or not raw["mask_name"]:
            raise CrossoverError(f"cell {index} has no mask name")
        requested = list(range(r)) if requested_raw == "all" else requested_raw
        if not isinstance(requested, list) or not requested or any(
            not isinstance(item, int) or isinstance(item, bool)
            or item < 0 or item >= r for item in requested
        ) or requested != sorted(set(requested)):
            raise CrossoverError(f"cell {index} has invalid requested parity")
        side_value = r if profile == "high" else k
        side = 1
        while side < side_value:
            side <<= 1
        order = 256 if field == "gf8" else 65536
        span = k + side if profile == "high" else side + r
        if side < 2 or span > order or (profile == "low" and side > order // 2):
            raise CrossoverError(f"cell {index} does not fit its field")
        cell = dict(raw)
        cell["requested_parity"] = requested
        cells.append(cell)
    cells.sort(key=lambda cell: (
        cell["field"], cell["profile"], cell["K"], cell["R"],
        cell["mask_name"], cell["backend"], cell["shard_bytes"],
    ))
    identities = [canonical_bytes(cell) for cell in cells]
    if len(identities) != len(set(identities)):
        raise CrossoverError("cell manifest contains duplicate cells")
    return cells


def load_cells(path: Path) -> list[dict[str, Any]]:
    return validate_cells(read_json(path))


def mask_text(indices: list[int]) -> str:
    return ",".join(str(value) for value in indices)


def benchmark_command(
    executable: Path, cell: dict[str, Any], settings: dict[str, Any]
) -> list[str]:
    return [
        str(executable),
        "--profile", cell["profile"],
        "--field", cell["field"],
        "--k", str(cell["K"]),
        "--r", str(cell["R"]),
        "--bytes", str(cell["shard_bytes"]),
        "--requested-parity", mask_text(cell["requested_parity"]),
        "--backend", cell["backend"],
        "--iterations", str(settings["iterations"]),
        "--samples", str(settings["samples"]),
        "--warmups", str(settings["warmups"]),
        "--setup-iterations", str(settings["setup_iterations"]),
        "--reuse", mask_text(settings["reuse"]),
        "--memory-mib", str(settings["memory_mib"]),
        "--seed", str(settings["seed"]),
    ]


def validate_metric(metric: Any, samples: int, name: str) -> None:
    if not isinstance(metric, dict) or not isinstance(metric.get("samples"), list):
        raise CrossoverError(f"missing metric {name}")
    values = metric["samples"]
    if len(values) != samples or not all(
        isinstance(value, (int, float)) and value > 0 for value in values
    ):
        raise CrossoverError(f"invalid samples for {name}")
    if abs(float(metric.get("median", -1)) - statistics.median(values)) > 0.002:
        raise CrossoverError(f"wrong median for {name}")
    if abs(float(metric.get("minimum", -1)) - min(values)) > 0.002 or (
        abs(float(metric.get("maximum", -1)) - max(values)) > 0.002
    ):
        raise CrossoverError(f"wrong extrema for {name}")
    if not isinstance(metric.get("mad"), (int, float)) or metric["mad"] < 0:
        raise CrossoverError(f"invalid MAD for {name}")


def validate_benchmark_result(
    result: Any,
    cell: dict[str, Any],
    settings: dict[str, Any],
    expected_git: dict[str, Any],
) -> None:
    if not isinstance(result, dict) or result.get("schema") != BENCHMARK_SCHEMA:
        raise CrossoverError("benchmark emitted an unknown schema")
    if result.get("authoritative") is not False:
        raise CrossoverError("raw benchmark improperly claimed authority")
    build = result.get("build")
    if not isinstance(build, dict):
        raise CrossoverError("benchmark omitted build identity")
    if build.get("source_git_sha") != expected_git["sha"]:
        raise CrossoverError("benchmark binary source SHA differs from runner source")
    if build.get("source_dirty") != int(expected_git["dirty"]):
        raise CrossoverError("benchmark binary dirty marker differs from runner source")
    if build.get("library_test_hooks") is not False:
        raise CrossoverError(
            "benchmark library contains test-hook instrumentation; "
            "configure authoritative builds with LEO2_BUILD_TESTS=OFF"
        )
    if not isinstance(build.get("compiler"), str) or not build["compiler"]:
        raise CrossoverError("benchmark binary omitted compiler identity")
    if not isinstance(build.get("compiler_version"), str) or not build["compiler_version"]:
        raise CrossoverError("benchmark binary omitted compiler version")
    if not isinstance(build.get("cplusplus"), int):
        raise CrossoverError("benchmark binary omitted C++ language identity")
    parameters = result.get("parameters")
    expected_parameters = {
        "profile": cell["profile"], "field": cell["field"],
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": cell["shard_bytes"],
        "requested_parity": cell["requested_parity"],
        "requested_backend": cell["backend"],
        "iterations": settings["iterations"], "samples": settings["samples"],
        "warmups": settings["warmups"],
        "setup_iterations": settings["setup_iterations"],
        "reuse": settings["reuse"], "memory_mib": settings["memory_mib"],
        "seed": settings["seed"],
    }
    if not isinstance(parameters, dict):
        raise CrossoverError("benchmark omitted parameters")
    for name, value in expected_parameters.items():
        if parameters.get(name) != value:
            raise CrossoverError(f"benchmark parameter differs: {name}")
    resolved = result.get("resolved")
    if not isinstance(resolved, dict) or not isinstance(resolved.get("backend"), str):
        raise CrossoverError("benchmark omitted resolved backend")
    if cell["backend"] != "auto" and resolved["backend"] != cell["backend"]:
        raise CrossoverError("benchmark resolved a different backend")
    plan = result.get("plan")
    if not isinstance(plan, dict):
        raise CrossoverError("benchmark omitted plan accounting")
    if not (
        isinstance(plan.get("retained_butterflies"), int)
        and isinstance(plan.get("prefix_butterflies"), int)
        and 0 < plan["retained_butterflies"] <= plan["prefix_butterflies"]
        <= plan.get("full_butterflies", -1)
    ):
        raise CrossoverError("cell emitted invalid forward-plan accounting")
    if plan["retained_butterflies"] == plan["prefix_butterflies"] and (
        cell["requested_parity"] != list(range(cell["R"]))
    ):
        raise CrossoverError("sparse requested-output cell retained prefix work")
    available = plan.get("inverse_candidate_available")
    inverse_operations = plan.get("inverse_operations")
    inverse_groups = plan.get("inverse_source_groups")
    inverse_prefix = plan.get("inverse_active_prefix")
    mature_zero_fill = plan.get("mature_zero_fill_shards")
    if not isinstance(available, bool) or not all(
        isinstance(value, int) and value >= 0
        for value in (
            inverse_operations, inverse_groups, inverse_prefix,
            mature_zero_fill,
        )
    ) or plan.get("pruned_zero_fill_shards") != 0:
        raise CrossoverError("benchmark emitted invalid inverse-plan accounting")
    padded_side = resolved.get("padded_side")
    if not isinstance(padded_side, int) or padded_side < 2:
        raise CrossoverError("benchmark omitted the padded transform side")
    if available:
        if not (
            0 < inverse_prefix < padded_side
            and inverse_operations > 0
            and inverse_groups == inverse_prefix // 4
            and mature_zero_fill == padded_side - inverse_prefix
        ):
            raise CrossoverError("benchmark emitted inconsistent inverse plan")
    elif any(
        value != 0
        for value in (
            inverse_operations, inverse_groups, inverse_prefix,
            mature_zero_fill,
        )
    ):
        raise CrossoverError("unavailable inverse plan retained accounting")
    correctness = result.get("correctness")
    if not isinstance(correctness, dict) or (
        correctness.get("exact_prefix_parity_match") is not True
    ):
        raise CrossoverError("exact parity differs from mature prefix parity")
    if correctness.get("inverse_pruned_parity_match") is not True:
        raise CrossoverError("inverse-pruned parity differs from mature parity")
    digest = correctness.get("digest_fnv1a64")
    if not isinstance(digest, str) or len(digest) != 18 or not digest.startswith("0x"):
        raise CrossoverError("benchmark emitted an invalid parity digest")
    metrics = result.get("metrics")
    if not isinstance(metrics, dict):
        raise CrossoverError("benchmark omitted metrics")
    names = (
        "schedule_setup_ns", "inverse_schedule_setup_ns",
        "prefix_execution_ns", "exact_prepared_execution_ns",
        "exact_call_local_total_ns", "prefix_pruned_inverse_execution_ns",
        "exact_pruned_inverse_execution_ns",
    )
    for name in names:
        validate_metric(metrics.get(name), settings["samples"], name)
    amortized = metrics.get("amortized_exact")
    if not isinstance(amortized, list) or (
        [row.get("reuse") for row in amortized] != settings["reuse"]
    ):
        raise CrossoverError("benchmark emitted wrong reuse rows")
    setup = float(metrics["schedule_setup_ns"]["median"])
    execution = float(metrics["exact_prepared_execution_ns"]["median"])
    prefix = float(metrics["prefix_execution_ns"]["median"])
    call_local = float(metrics["exact_call_local_total_ns"]["median"])
    prefix_inverse = float(metrics["prefix_pruned_inverse_execution_ns"]["median"])
    exact_inverse = float(metrics["exact_pruned_inverse_execution_ns"]["median"])
    if abs(float(metrics.get("prefix_over_exact_prepared", -1)) - prefix / execution) > 0.002:
        raise CrossoverError("benchmark emitted wrong prepared-exact ratio")
    if abs(float(metrics.get("prefix_over_exact_call_local", -1)) - prefix / call_local) > 0.002:
        raise CrossoverError("benchmark emitted wrong call-local ratio")
    if abs(float(metrics.get("mature_over_pruned_inverse_prefix", -1)) - prefix / prefix_inverse) > 0.002:
        raise CrossoverError("benchmark emitted wrong inverse-prefix ratio")
    if abs(float(metrics.get("mature_over_pruned_inverse_exact", -1)) - execution / exact_inverse) > 0.002:
        raise CrossoverError("benchmark emitted wrong inverse-exact ratio")
    for row in amortized:
        modeled = execution + setup / row["reuse"]
        if abs(float(row.get("modeled_ns", -1)) - modeled) > 0.002:
            raise CrossoverError("benchmark emitted wrong amortized total")
        if abs(float(row.get("prefix_over_modeled_exact", -1)) - prefix / modeled) > 0.002:
            raise CrossoverError("benchmark emitted wrong amortized ratio")


def load_attestation(path: Path, cpu: int) -> dict[str, Any]:
    value = read_json(path)
    if not isinstance(value, dict) or value.get("schema") != ATTESTATION_SCHEMA:
        raise CrossoverError("isolation attestation has an unknown schema")
    if value.get("cpu") != cpu:
        raise CrossoverError("isolation attestation names a different CPU")
    if value.get("smt_sibling_idle") is not True or (
        value.get("competing_work_idle") is not True
    ):
        raise CrossoverError("isolation attestation does not prove an idle host")
    if not isinstance(value.get("operator"), str) or not value["operator"].strip():
        raise CrossoverError("isolation attestation omits the operator")
    if not isinstance(value.get("timestamp_utc"), str) or not value["timestamp_utc"].strip():
        raise CrossoverError("isolation attestation omits the timestamp")
    return value


def make_manifest(
    source: Path,
    executable: Path,
    mode: str,
    cells: list[dict[str, Any]],
    settings: dict[str, Any],
    pin_cpu: int | None,
    attestation: dict[str, Any] | None,
) -> dict[str, Any]:
    cpus = allowed_cpus()
    fingerprint = source_fingerprint(source)
    state = git_state(source)
    executable_sha = digest_bytes(executable.read_bytes())
    build = build_metadata(executable)
    machine = machine_identity(cpus)
    identity = {
        "cells": cells,
        "executable_sha256": executable_sha,
        "git": state,
        "machine": machine,
        "mode": mode,
        "pin_cpu": pin_cpu,
        "settings": settings,
        "source_fingerprint": fingerprint,
        "attestation": attestation,
        "build_metadata": build,
    }
    configuration_id = digest_value(identity)
    jobs = []
    for cell in cells:
        job_identity = {
            "cell": cell,
            "configuration_id": configuration_id,
            "executable_sha256": executable_sha,
        }
        jobs.append({
            "cell": cell,
            "configuration_id": configuration_id,
            "job_id": digest_value(job_identity)[:20],
        })
    return {
        "attestation": attestation,
        "authoritative_requested": mode == "pinned",
        "build_metadata": build,
        "configuration_id": configuration_id,
        "executable": str(executable),
        "executable_sha256": executable_sha,
        "git": state,
        "jobs": jobs,
        "machine": machine,
        "mode": mode,
        "pin_cpu": pin_cpu,
        "schema": SCHEMA,
        "settings": settings,
        "source": str(source),
        "source_fingerprint": fingerprint,
    }


def validate_job_artifacts(
    result_dir: Path,
    job: dict[str, Any],
    expected: dict[str, Any],
    manifest: dict[str, Any],
) -> None:
    if job.get("schema") != JOB_SCHEMA or job.get("job_id") != expected["job_id"]:
        raise CrossoverError("job artifact has stale identity")
    if job.get("configuration_id") != manifest["configuration_id"]:
        raise CrossoverError("job artifact belongs to another manifest")
    for stream in ("stdout", "stderr"):
        relative = job.get(f"{stream}_path")
        expected_digest = job.get(f"{stream}_sha256")
        if not isinstance(relative, str) or not isinstance(expected_digest, str):
            raise CrossoverError("job artifact omits retained streams")
        path = result_dir / relative
        if not path.is_file() or digest_bytes(path.read_bytes()) != expected_digest:
            raise CrossoverError(f"job {stream} artifact changed")
    if job.get("status") == "passed":
        validate_benchmark_result(
            job.get("benchmark"), expected["cell"], manifest["settings"], manifest["git"]
        )


def run_job(
    job: dict[str, Any], manifest: dict[str, Any], result_dir: Path,
    executable: Path, taskset: str | None, timeout: int, resume: bool,
) -> dict[str, Any]:
    job_dir = result_dir / "jobs"
    result_path = job_dir / f"{job['job_id']}.json"
    if resume and result_path.is_file():
        existing = read_json(result_path)
        validate_job_artifacts(result_dir, existing, job, manifest)
        return existing
    command = benchmark_command(executable, job["cell"], manifest["settings"])
    if taskset is not None:
        command = [taskset, "-c", str(manifest["pin_cpu"]), *command]
    environment = os.environ.copy()
    environment.update(OMP_DYNAMIC="FALSE", OMP_NUM_THREADS="1")
    completed = subprocess.run(
        command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=timeout, env=environment,
    )
    stdout_relative = f"jobs/{job['job_id']}.stdout"
    stderr_relative = f"jobs/{job['job_id']}.stderr"
    atomic_write(result_dir / stdout_relative, completed.stdout)
    atomic_write(result_dir / stderr_relative, completed.stderr)
    benchmark = None
    error = None
    status = "failed"
    if completed.returncode == 0:
        try:
            benchmark = json.loads(completed.stdout.decode("utf-8"))
            validate_benchmark_result(
                benchmark, job["cell"], manifest["settings"], manifest["git"]
            )
            status = "passed"
        except (UnicodeError, ValueError, CrossoverError) as caught:
            error = str(caught)
    else:
        error = f"benchmark exited {completed.returncode}"
    result = {
        "benchmark": benchmark,
        "cell": job["cell"],
        "command": command,
        "configuration_id": manifest["configuration_id"],
        "error": error,
        "job_id": job["job_id"],
        "returncode": completed.returncode,
        "schema": JOB_SCHEMA,
        "status": status,
        "stderr_path": stderr_relative,
        "stderr_sha256": digest_bytes(completed.stderr),
        "stdout_path": stdout_relative,
        "stdout_sha256": digest_bytes(completed.stdout),
    }
    atomic_write_json(result_path, result)
    return result


def summarize_jobs(jobs: list[dict[str, Any]], reuse: list[int]) -> dict[str, Any]:
    passed = [job for job in jobs if job.get("status") == "passed"]
    prepared: list[float] = []
    call_local: list[float] = []
    amortized: dict[int, list[float]] = {value: [] for value in reuse}
    inverse_prefix: list[float] = []
    inverse_exact: list[float] = []
    inverse_prefix_amortized: dict[int, list[float]] = {
        value: [] for value in reuse
    }
    inverse_exact_amortized: dict[int, list[float]] = {
        value: [] for value in reuse
    }
    for job in passed:
        metrics = job["benchmark"]["metrics"]
        prefix = float(metrics["prefix_execution_ns"]["median"])
        exact = float(metrics["exact_prepared_execution_ns"]["median"])
        setup = float(metrics["schedule_setup_ns"]["median"])
        call = float(metrics["exact_call_local_total_ns"]["median"])
        prepared.append((prefix / exact - 1.0) * 100.0)
        call_local.append((prefix / call - 1.0) * 100.0)
        for reuse_count in reuse:
            modeled = exact + setup / reuse_count
            amortized[reuse_count].append((prefix / modeled - 1.0) * 100.0)
        if job["benchmark"]["plan"]["inverse_candidate_available"]:
            inverse_setup = float(
                metrics["inverse_schedule_setup_ns"]["median"]
            )
            pruned_prefix = float(
                metrics["prefix_pruned_inverse_execution_ns"]["median"]
            )
            pruned_exact = float(
                metrics["exact_pruned_inverse_execution_ns"]["median"]
            )
            inverse_prefix.append((prefix / pruned_prefix - 1.0) * 100.0)
            inverse_exact.append((exact / pruned_exact - 1.0) * 100.0)
            for reuse_count in reuse:
                inverse_prefix_modeled = (
                    pruned_prefix + inverse_setup / reuse_count
                )
                inverse_exact_modeled = pruned_exact + inverse_setup / reuse_count
                inverse_prefix_amortized[reuse_count].append(
                    (prefix / inverse_prefix_modeled - 1.0) * 100.0
                )
                inverse_exact_amortized[reuse_count].append(
                    (exact / inverse_exact_modeled - 1.0) * 100.0
                )

    def summary(values: list[float]) -> dict[str, Any]:
        return {
            "cells": len(values),
            "gain_max_percent": max(values) if values else None,
            "gain_median_percent": statistics.median(values) if values else None,
            "gain_min_percent": min(values) if values else None,
            "regressions": sum(value < 0 for value in values),
            "promotions_at_5_percent": sum(value >= 5.0 for value in values),
            "severe_regressions_below_minus_2_percent":
                sum(value < -2.0 for value in values),
        }

    return {
        "call_local_total": summary(call_local),
        "exact_prepared_execution": summary(prepared),
        "inverse_pruned_prefix_execution": summary(inverse_prefix),
        "inverse_pruned_exact_execution": summary(inverse_exact),
        "jobs_failed": len(jobs) - len(passed),
        "jobs_passed": len(passed),
        "jobs_total": len(jobs),
        "modeled_amortized": {
            str(value): summary(amortized[value]) for value in reuse
        },
        "modeled_amortized_inverse_prefix": {
            str(value): summary(inverse_prefix_amortized[value])
            for value in reuse
        },
        "modeled_amortized_inverse_exact": {
            str(value): summary(inverse_exact_amortized[value])
            for value in reuse
        },
    }


def load_manifest(result_dir: Path) -> dict[str, Any]:
    manifest = read_json(result_dir / "manifest.json")
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise CrossoverError("manifest has an unknown schema")
    if not isinstance(manifest.get("jobs"), list) or not manifest.get("configuration_id"):
        raise CrossoverError("manifest is incomplete")
    return manifest


def load_jobs(result_dir: Path, manifest: dict[str, Any]) -> list[dict[str, Any]]:
    expected = {job["job_id"]: job for job in manifest["jobs"]}
    paths = {path.stem: path for path in sorted((result_dir / "jobs").glob("*.json"))}
    if set(paths) != set(expected):
        raise CrossoverError("job file set differs from manifest")
    jobs = []
    for job_id in sorted(expected):
        value = read_json(paths[job_id])
        validate_job_artifacts(result_dir, value, expected[job_id], manifest)
        jobs.append(value)
    return jobs


def require_compatible_directory(result_dir: Path, manifest: dict[str, Any]) -> None:
    path = result_dir / "manifest.json"
    if not path.exists():
        if (result_dir / "jobs").is_dir() and any((result_dir / "jobs").iterdir()):
            raise CrossoverError("result directory has jobs but no manifest")
        return
    previous = load_manifest(result_dir)
    if previous["configuration_id"] != manifest["configuration_id"]:
        raise CrossoverError("result directory belongs to a different configuration")


def run_matrix(arguments: argparse.Namespace) -> int:
    source = Path(arguments.source).resolve()
    executable = Path(arguments.executable).resolve()
    result_dir = Path(arguments.result_dir).resolve()
    if not executable.is_file():
        raise CrossoverError(f"benchmark executable is missing: {executable}")
    cpus = allowed_cpus()
    pin_cpu: int | None = None
    taskset: str | None = None
    attestation = None
    state = git_state(source)
    if arguments.command == "pinned":
        if state["dirty"]:
            raise CrossoverError("pinned evidence requires a clean source tree")
        pin_cpu = arguments.cpu
        if pin_cpu is None or pin_cpu not in cpus:
            raise CrossoverError(
                f"pinned CPU must be in allowed affinity {compact_cpu_list(cpus)}"
            )
        if arguments.workers != 1:
            raise CrossoverError("pinned measurements require exactly one worker")
        taskset = shutil.which(arguments.taskset)
        if taskset is None:
            raise CrossoverError("pinned measurements require taskset")
        if not arguments.isolation_attestation:
            raise CrossoverError("pinned measurements require --isolation-attestation")
        attestation = load_attestation(
            Path(arguments.isolation_attestation).resolve(), pin_cpu
        )
    shard_bytes = parse_csv_unsigned(arguments.bytes, "bytes", 1 << 30)
    reuse = parse_csv_unsigned(arguments.reuse, "reuse", 1000000)
    backends = parse_backends(arguments.backends)
    if arguments.cell_manifest:
        cell_manifest = Path(arguments.cell_manifest).resolve()
        cells = load_cells(cell_manifest)
        settings_cell_manifest = {
            "path": str(cell_manifest),
            "sha256": digest_bytes(cell_manifest.read_bytes()),
        }
    else:
        cells = make_cells(backends, shard_bytes)
        settings_cell_manifest = None
    if arguments.workers is None:
        arguments.workers = min(128, len(cpus), len(cells))
    settings = {
        "iterations": arguments.iterations,
        "memory_mib": arguments.memory_mib,
        "placement_policy": "taskset single CPU" if pin_cpu is not None
            else "inherited allowed affinity with independent worker processes",
        "reuse": reuse,
        "samples": arguments.samples,
        "seed": arguments.seed,
        "setup_iterations": arguments.setup_iterations,
        "timeout_seconds": arguments.timeout,
        "warmups": arguments.warmups,
        "workers": arguments.workers,
        "cell_manifest": settings_cell_manifest,
    }
    manifest = make_manifest(
        source, executable, arguments.command, cells, settings,
        pin_cpu, attestation,
    )
    require_compatible_directory(result_dir, manifest)
    atomic_write_json(result_dir / "manifest.json", manifest)
    expected_job_ids = {job["job_id"] for job in manifest["jobs"]}
    existing_job_ids = {
        path.stem for path in (result_dir / "jobs").glob("*.json")
    } if (result_dir / "jobs").is_dir() else set()
    extra_job_ids = sorted(existing_job_ids - expected_job_ids)
    if extra_job_ids:
        raise CrossoverError(
            "result directory contains stale job artifacts: "
            + ",".join(extra_job_ids)
        )
    print(
        f"sparse encode crossover: mode={arguments.command} cells={len(cells)} "
        f"allowed={manifest['machine']['allowed_cpu_list']}"
        + (f" pinned={pin_cpu}" if pin_cpu is not None else ""),
        flush=True,
    )
    jobs: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=arguments.workers) as pool:
        futures = {
            pool.submit(
                run_job, job, manifest, result_dir, executable,
                taskset, arguments.timeout, not arguments.no_resume,
            ): job
            for job in manifest["jobs"]
        }
        for future in concurrent.futures.as_completed(futures):
            job = futures[future]
            result = future.result()
            jobs.append(result)
            print(f"{job['job_id']} {result['status']}", flush=True)
    jobs.sort(key=lambda value: value["job_id"])
    fingerprint_after = source_fingerprint(source)
    git_after = git_state(source)
    executable_after = digest_bytes(executable.read_bytes())
    source_changed = (
        fingerprint_after["digest"] != manifest["source_fingerprint"]["digest"]
        or git_after != manifest["git"]
    )
    executable_changed = executable_after != manifest["executable_sha256"]
    failed = any(job["status"] != "passed" for job in jobs)
    authoritative = (
        arguments.command == "pinned" and not failed and not source_changed
        and not executable_changed and attestation is not None
    )
    status = "passed" if not failed and not source_changed and not executable_changed else "failed"
    summary = summarize_jobs(jobs, reuse)
    matrix = {
        "attestation": attestation,
        "authoritative": authoritative,
        "executable_changed_during_run": executable_changed,
        "jobs": jobs,
        "manifest_configuration_id": manifest["configuration_id"],
        "schema": SCHEMA,
        "source_changed_during_run": source_changed,
        "status": status,
        "summary": summary,
    }
    atomic_write_json(result_dir / "matrix.json", matrix)
    atomic_write_json(result_dir / "summary.json", summary)
    print(
        f"matrix: {status}; passed={summary['jobs_passed']}/{summary['jobs_total']}; "
        f"authoritative={str(authoritative).lower()} ({result_dir / 'matrix.json'})"
    )
    return 0 if status == "passed" else 1


def analyze(arguments: argparse.Namespace) -> int:
    result_dir = Path(arguments.result_dir).resolve()
    manifest = load_manifest(result_dir)
    source = Path(manifest.get("source", "")).resolve()
    executable = Path(manifest.get("executable", "")).resolve()
    if source_fingerprint(source) != manifest.get("source_fingerprint"):
        raise CrossoverError("current source fingerprint differs from the manifest")
    if git_state(source) != manifest.get("git"):
        raise CrossoverError("current Git identity differs from the manifest")
    if not executable.is_file() or (
        digest_bytes(executable.read_bytes()) != manifest.get("executable_sha256")
    ):
        raise CrossoverError("current executable identity differs from the manifest")
    jobs = load_jobs(result_dir, manifest)
    matrix = read_json(result_dir / "matrix.json")
    if not isinstance(matrix, dict) or matrix.get("schema") != SCHEMA:
        raise CrossoverError("matrix has an unknown schema")
    if matrix.get("manifest_configuration_id") != manifest["configuration_id"]:
        raise CrossoverError("matrix belongs to another manifest")
    if canonical_bytes(matrix.get("jobs")) != canonical_bytes(jobs):
        raise CrossoverError("matrix job snapshot differs from retained jobs")
    summary = summarize_jobs(jobs, manifest["settings"]["reuse"])
    if canonical_bytes(summary) != canonical_bytes(matrix.get("summary")):
        raise CrossoverError("matrix summary differs from validated jobs")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if matrix.get("status") == "passed" else 1


def self_test() -> int:
    cells_a = make_cells(["avx2", "scalar"], [1024, 64])
    cells_b = make_cells(["scalar", "avx2"], [64, 1024])
    if canonical_bytes(cells_a) != canonical_bytes(cells_b) or len(cells_a) != 32:
        raise CrossoverError("cell generation is not deterministic")
    if parse_csv_unsigned("64,1,8", "reuse", 100) != [1, 8, 64]:
        raise CrossoverError("numeric list normalization failed")
    custom_cells = validate_cells({
        "schema": CELL_SCHEMA,
        "cells": [{
            "profile": "low", "field": "gf8", "K": 17, "R": 18,
            "mask_name": "all", "requested_parity": list(range(18)),
            "backend": "scalar", "shard_bytes": 65536,
        }],
    })
    if len(custom_cells) != 1 or custom_cells[0]["K"] != 17:
        raise CrossoverError("custom cell manifest validation failed")
    malformed_cells = dict(custom_cells[0])
    malformed_cells["requested_parity"] = [18]
    try:
        validate_cells({"schema": CELL_SCHEMA, "cells": [malformed_cells]})
    except CrossoverError:
        pass
    else:
        raise CrossoverError("out-of-range custom cell was accepted")
    with tempfile.TemporaryDirectory(prefix="leopard2-sparse-crossover-") as directory:
        root = Path(directory)
        for relative in SOURCE_FILES:
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(relative + "\n", encoding="utf-8")
        before = source_fingerprint(root)
        (root / SOURCE_FILES[0]).write_text("changed\n", encoding="utf-8")
        after = source_fingerprint(root)
        if before["digest"] == after["digest"]:
            raise CrossoverError("source fingerprint ignored a mutation")
        attestation = {
            "schema": ATTESTATION_SCHEMA, "cpu": 7,
            "smt_sibling_idle": True, "competing_work_idle": True,
            "operator": "self-test", "timestamp_utc": "2000-01-01T00:00:00Z",
        }
        path = root / "attestation.json"
        atomic_write_json(path, attestation)
        if load_attestation(path, 7) != attestation:
            raise CrossoverError("attestation round trip failed")
        attestation["smt_sibling_idle"] = False
        atomic_write_json(path, attestation)
        try:
            load_attestation(path, 7)
        except CrossoverError:
            pass
        else:
            raise CrossoverError("unsafe isolation attestation was accepted")

    cell = dict(base_cells()[0])
    cell.update(backend="scalar", shard_bytes=64)
    settings = {
        "iterations": 2, "memory_mib": 16, "reuse": [1, 8, 64],
        "samples": 3, "seed": 7, "setup_iterations": 2, "warmups": 1,
    }
    state = {"dirty": False, "sha": "a" * 40}
    setup = 10.0
    prefix = 30.0
    exact = 20.0
    benchmark = {
        "schema": BENCHMARK_SCHEMA,
        "authoritative": False,
        "build": {
            "source_git_sha": state["sha"], "source_dirty": 0,
            "library_test_hooks": False,
            "compiler": "self-test", "compiler_version": "1", "cplusplus": 201103,
        },
        "parameters": {
            "profile": cell["profile"], "field": cell["field"],
            "K": cell["K"], "R": cell["R"], "shard_bytes": 64,
            "requested_parity": cell["requested_parity"],
            "requested_backend": "scalar", "iterations": 2, "samples": 3,
            "warmups": 1, "setup_iterations": 2, "reuse": [1, 8, 64],
            "memory_mib": 16, "seed": 7,
        },
        "resolved": {"backend": "scalar", "padded_side": 64},
        "plan": {
            "full_butterflies": 100, "prefix_butterflies": 80,
            "retained_butterflies": 50,
            "inverse_candidate_available": False,
            "inverse_operations": 0, "inverse_source_groups": 0,
            "inverse_active_prefix": 0, "mature_zero_fill_shards": 0,
            "pruned_zero_fill_shards": 0,
        },
        "correctness": {
            "exact_prefix_parity_match": True,
            "inverse_pruned_parity_match": True,
            "digest_fnv1a64": "0x0123456789abcdef",
        },
        "metrics": {
            "schedule_setup_ns": {
                "median": setup, "mad": 0, "minimum": 10, "maximum": 10,
                "samples": [10, 10, 10],
            },
            "inverse_schedule_setup_ns": {
                "median": setup, "mad": 0, "minimum": 10, "maximum": 10,
                "samples": [10, 10, 10],
            },
            "prefix_execution_ns": {
                "median": prefix, "mad": 0, "minimum": 30, "maximum": 30,
                "samples": [30, 30, 30],
            },
            "exact_prepared_execution_ns": {
                "median": exact, "mad": 0, "minimum": 20, "maximum": 20,
                "samples": [20, 20, 20],
            },
            "exact_call_local_total_ns": {
                "median": 25, "mad": 0, "minimum": 25, "maximum": 25,
                "samples": [25, 25, 25],
            },
            "prefix_pruned_inverse_execution_ns": {
                "median": 24, "mad": 0, "minimum": 24, "maximum": 24,
                "samples": [24, 24, 24],
            },
            "exact_pruned_inverse_execution_ns": {
                "median": 18, "mad": 0, "minimum": 18, "maximum": 18,
                "samples": [18, 18, 18],
            },
            "prefix_over_exact_prepared": prefix / exact,
            "prefix_over_exact_call_local": prefix / 25,
            "mature_over_pruned_inverse_prefix": prefix / 24,
            "mature_over_pruned_inverse_exact": exact / 18,
            "amortized_exact": [
                {"reuse": reuse, "modeled_ns": exact + setup / reuse,
                 "prefix_over_modeled_exact": prefix / (exact + setup / reuse)}
                for reuse in (1, 8, 64)
            ],
        },
    }
    validate_benchmark_result(benchmark, cell, settings, state)
    unavailable_summary = summarize_jobs(
        [{"status": "passed", "benchmark": benchmark}], [1, 8, 64]
    )
    if unavailable_summary["inverse_pruned_prefix_execution"]["cells"] != 0:
        raise CrossoverError("summary included an unavailable inverse candidate")
    candidate_benchmark = json.loads(json.dumps(benchmark))
    candidate_benchmark["plan"].update(
        inverse_candidate_available=True,
        inverse_operations=63,
        inverse_source_groups=4,
        inverse_active_prefix=17,
        mature_zero_fill_shards=47,
    )
    candidate_summary = summarize_jobs(
        [{"status": "passed", "benchmark": candidate_benchmark}], [1, 8, 64]
    )
    inverse_gain = candidate_summary[
        "inverse_pruned_prefix_execution"
    ]["gain_median_percent"]
    if inverse_gain is None or abs(inverse_gain - 25.0) > 0.001:
        raise CrossoverError("summary computed the wrong inverse-prefix gain")
    benchmark["build"]["source_git_sha"] = "b" * 40
    try:
        validate_benchmark_result(benchmark, cell, settings, state)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("mismatched binary source SHA was accepted")
    benchmark["build"]["source_git_sha"] = state["sha"]
    benchmark["build"]["source_dirty"] = 1
    try:
        validate_benchmark_result(benchmark, cell, settings, state)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("mismatched binary dirty marker was accepted")
    benchmark["build"]["source_dirty"] = 0
    benchmark["build"]["library_test_hooks"] = True
    try:
        validate_benchmark_result(benchmark, cell, settings, state)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("instrumented benchmark archive was accepted")
    del benchmark["build"]["library_test_hooks"]
    try:
        validate_benchmark_result(benchmark, cell, settings, state)
    except CrossoverError:
        pass
    else:
        raise CrossoverError("benchmark without instrumentation identity was accepted")
    print("PASS sparse encode crossover self-test cells=32 identity_mutations=6")
    return 0


def positive(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("screen", "pinned"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--source", default=".")
        sub.add_argument("--executable", required=True)
        sub.add_argument("--result-dir", required=True)
        sub.add_argument("--cell-manifest")
        sub.add_argument("--backends", default="auto")
        sub.add_argument(
            "--bytes", default="64,1024" if command == "screen"
            else "64,1024,65536,262144"
        )
        sub.add_argument("--reuse", default="1,8,64")
        sub.add_argument("--iterations", type=positive, default=8 if command == "screen" else 64)
        sub.add_argument("--samples", type=positive, default=3 if command == "screen" else 9)
        sub.add_argument("--warmups", type=positive, default=1 if command == "screen" else 4)
        sub.add_argument(
            "--setup-iterations", type=positive,
            default=8 if command == "screen" else 64,
        )
        sub.add_argument("--memory-mib", type=positive, default=768)
        sub.add_argument("--seed", type=positive, default=0x535041525345454E)
        sub.add_argument("--workers", type=positive, default=None if command == "screen" else 1)
        sub.add_argument("--timeout", type=positive, default=180)
        sub.add_argument("--no-resume", action="store_true")
        sub.add_argument("--cpu", type=int)
        sub.add_argument("--taskset", default="taskset")
        sub.add_argument("--isolation-attestation")
    analyze_parser = subparsers.add_parser("analyze")
    analyze_parser.add_argument("--result-dir", required=True)
    subparsers.add_parser("self-test")
    return parser


def main() -> int:
    arguments = build_parser().parse_args()
    if arguments.command in ("screen", "pinned"):
        if arguments.samples < 3 or arguments.samples > 101 or arguments.samples % 2 == 0:
            raise CrossoverError("samples must be odd and in 3..101")
        if arguments.command == "screen" and (
            arguments.cpu is not None or arguments.isolation_attestation is not None
        ):
            raise CrossoverError("screen mode does not accept pinned-only options")
        return run_matrix(arguments)
    if arguments.command == "analyze":
        return analyze(arguments)
    return self_test()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (CrossoverError, subprocess.TimeoutExpired) as error:
        print(f"sparse encode crossover: {error}", file=sys.stderr)
        raise SystemExit(1)
