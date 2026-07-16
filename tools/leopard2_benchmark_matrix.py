#!/usr/bin/env python3
"""Generate and collect deterministic Leopard2 benchmark matrices.

The generated job specification is consumed by tools/leopard2_lab.py.  Timing
jobs write benchmark JSON to stdout so the lab runner can retain stdout,
stderr, affinity, timeout, and exit status independently for every cell.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


SCHEMA = "leopard2-benchmark-matrix/v2"
SPEC_SCHEMA = "leopard2-benchmark-spec/v2"
LAB_MANIFEST_SCHEMA = "leopard2-lab-manifest/v2"
MODE_AUTOMATIC = "automatic"
MODE_FORCED_SPECIALIZED = "forced-specialized"
MODE_FORCED_GENERIC = "forced-generic"
EXPECTED_CELL_FIELDS = (
    "K",
    "R",
    "requested_profile",
    "requested_field",
    "requested_backend",
    "force_generic_decode",
    "force_specialized_decode",
    "shard_bytes",
    "loss_count",
    "batch",
    "reuse",
    "iterations",
    "warmup",
    "thread_count",
    "seed",
)

RATE_CASES = {
    "low": (
        (8, 248, "low"),
        (16, 240, "low"),
        (32, 224, "low"),
        (64, 192, "low"),
        (100, 156, "low"),
        (127, 129, "low"),
    ),
    "balanced": ((128, 128, "high"), (128, 128, "low"), (256, 256, "high")),
    "high": (
        (192, 64, "high"),
        (224, 32, "high"),
        (240, 16, "high"),
        (248, 8, "high"),
        (1000, 200, "high"),
        (4096, 512, "high"),
    ),
    "padding": (
        (129, 1, "high"),
        (129, 100, "high"),
        (200, 50, "high"),
        (225, 30, "high"),
    ),
}

FULL_SHARD_BYTES = (
    64,
    256,
    1024,
    4096,
    16384,
    65536,
    262144,
    1048576,
    4194304,
    16777216,
)


class MatrixError(ValueError):
    """Invalid matrix input or result."""


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def digest(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_identity(path: Path) -> dict:
    try:
        data = path.read_bytes()
    except OSError as error:
        raise MatrixError(f"cannot read {path}: {error}") from error
    return {"sha256": hashlib.sha256(data).hexdigest(), "size_bytes": len(data)}


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True, ensure_ascii=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def load_json(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except OSError as error:
        raise MatrixError(f"cannot read {path}: {error}") from error
    except ValueError as error:
        raise MatrixError(f"invalid JSON in {path}: {error}") from error


def loss_counts(k: int, r: int) -> Tuple[int, ...]:
    maximum = min(k, r)
    candidates = (0, 1, 2, 4, 8, r // 4, r // 2, maximum)
    return tuple(sorted({value for value in candidates if 0 <= value <= maximum}))


def memory_limit_mb(k: int, r: int, shard_bytes: int, batch: int) -> int:
    # The benchmark retains originals, parity, received copies, legacy work,
    # and codec scratch. Eight public-codeword equivalents plus 2 GiB leaves
    # headroom for the generic parent workspace without disabling the runner's
    # per-process safety limit.
    public_bytes = (k + r) * shard_bytes * batch
    estimate = public_bytes * 8 + 2 * 1024 * 1024 * 1024
    # Do not cap this at 64 GiB: large required-matrix cells can need hundreds
    # of GiB and an artificially low RLIMIT_AS guarantees a meaningless
    # allocation failure.  leopard2_lab records cells that exceed 80% of host
    # memory as explicitly unavailable before launching them.
    return max(2048, (estimate + (1 << 20) - 1) >> 20)


def _job(
    benchmark: str,
    category: str,
    k: int,
    r: int,
    profile: str,
    shard_bytes: int,
    losses: int,
    batch: int,
    reuse: int,
    threads: int,
    mode: str,
    iterations: int,
    warmup: int,
) -> Dict[str, object]:
    if mode not in (MODE_AUTOMATIC, MODE_FORCED_SPECIALIZED, MODE_FORCED_GENERIC):
        raise MatrixError(f"invalid decoder request mode: {mode}")
    job_id = (
        f"{category}.k{k}.r{r}.{profile}.b{shard_bytes}.l{losses}."
        f"batch{batch}.reuse{reuse}.t{threads}.{mode}"
    )
    pair_identity = (
        f"{category}|{k}|{r}|{profile}|{shard_bytes}|{losses}|"
        f"{batch}|{reuse}|{threads}"
    ).encode("ascii")
    pair_seed = int.from_bytes(hashlib.sha256(pair_identity).digest()[:4], "big")
    cpu_group = "pair-" + hashlib.sha256(pair_identity).hexdigest()
    command = [
        benchmark,
        "--k",
        str(k),
        "--r",
        str(r),
        "--profile",
        profile,
        "--field",
        "auto",
        "--backend",
        "auto",
        "--bytes",
        str(shard_bytes),
        "--loss",
        str(losses),
        "--batch",
        str(batch),
        "--reuse",
        str(reuse),
        "--iterations",
        str(iterations),
        "--warmup",
        str(warmup),
        "--threads",
        str(threads),
        "--seed",
        str(pair_seed),
        "--json",
        "-",
    ]
    if mode == MODE_FORCED_GENERIC:
        command.append("--force-generic")
    elif mode == MODE_FORCED_SPECIALIZED:
        command.append("--force-specialized")
    requested_profile = {
        "high": "legacy_high_v1",
        "low": "low_v1",
        "exact": "exact_experimental_v1",
        "auto": "auto",
    }[profile]
    benchmark_cell = {
        "K": k,
        "R": r,
        "requested_profile": requested_profile,
        "requested_field": "auto",
        "requested_backend": "auto",
        "force_generic_decode": mode == MODE_FORCED_GENERIC,
        "force_specialized_decode": mode == MODE_FORCED_SPECIALIZED,
        "shard_bytes": shard_bytes,
        "loss_count": losses,
        "batch": batch,
        "reuse": reuse,
        "iterations": iterations,
        "warmup": warmup,
        "thread_count": threads,
        "seed": pair_seed,
    }
    memory_mb = memory_limit_mb(k, r, shard_bytes, batch)
    return {
        "id": job_id,
        "command": command,
        "cpu_count": threads,
        "cpu_group": cpu_group,
        "benchmark_cell": benchmark_cell,
        "memory_mb": memory_mb,
        "minimum_memory_mb": memory_mb,
        "timeout_seconds": 1800,
    }


def _checkpoint_cells() -> Iterable[Tuple[str, int, int, str]]:
    yield "low", 16, 240, "low"
    yield "balanced", 128, 128, "high"
    yield "high", 240, 16, "high"


def make_spec(
    benchmark: str,
    preset: str,
    workers: int,
    iterations: int,
    warmup: int,
    pinned_cpu: int | None = None,
) -> dict:
    if workers <= 0 or workers > 128:
        raise MatrixError("workers must be in [1, 128]")
    jobs: List[Dict[str, object]] = []
    if preset == "smoke":
        cells = (
            ("low", 16, 240, "low", 4096, 1),
            ("high", 240, 16, "high", 4096, 1),
        )
        for category, k, r, profile, shard_bytes, losses in cells:
            for mode in (MODE_AUTOMATIC, MODE_FORCED_SPECIALIZED, MODE_FORCED_GENERIC):
                jobs.append(_job(
                    benchmark, category, k, r, profile, shard_bytes, losses,
                    1, 2, 1, mode, iterations, warmup))
    elif preset == "checkpoint":
        for category, k, r, profile in _checkpoint_cells():
            for shard_bytes in (4096, 65536, 1048576):
                for losses in (0, 1, min(8, k, r)):
                    modes = ((MODE_AUTOMATIC,) if losses == 0 else
                             (MODE_AUTOMATIC, MODE_FORCED_SPECIALIZED,
                              MODE_FORCED_GENERIC))
                    for mode in modes:
                        jobs.append(_job(
                            benchmark, category, k, r, profile, shard_bytes,
                            losses, 1, 8, 1, mode, iterations, warmup))
    elif preset == "balanced-crossover":
        for shard_bytes in (256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536):
            for losses in (1, 2, 4, 8, 16, 32, 64, 128):
                for mode in (MODE_AUTOMATIC, MODE_FORCED_SPECIALIZED,
                             MODE_FORCED_GENERIC):
                    jobs.append(_job(
                        benchmark, "balanced-crossover", 128, 128, "high",
                        shard_bytes, losses, 1, 16, 1, mode,
                        iterations, warmup))
    elif preset == "required":
        for category in sorted(RATE_CASES):
            for k, r, profile in RATE_CASES[category]:
                for shard_bytes in FULL_SHARD_BYTES:
                    for losses in loss_counts(k, r):
                        modes = ((MODE_FORCED_SPECIALIZED, MODE_FORCED_GENERIC)
                                 if losses else (MODE_AUTOMATIC,))
                        for mode in modes:
                            jobs.append(_job(
                                benchmark, category, k, r, profile, shard_bytes,
                                losses, 1, 8, 1, mode, iterations, warmup))
        # Reuse/batch and thread scaling are separate from the count/size/loss
        # Cartesian product to avoid impossible multi-terabyte allocations.
        for category, k, r, profile in _checkpoint_cells():
            for batch, reuse in ((8, 8), (64, 64), (1024, 1024)):
                shard_bytes = 65536 if batch <= 8 else (4096 if batch <= 64 else 256)
                jobs.append(_job(
                    benchmark, category + "-reuse", k, r, profile,
                    shard_bytes, min(8, k, r), batch, reuse, 1,
                    MODE_AUTOMATIC,
                    iterations, warmup))
            for threads in (1, 2, 4, 8, 16, 32, 64, 128):
                jobs.append(_job(
                    benchmark, category + "-scaling", k, r, profile,
                    4096, min(8, k, r), 128, 8, threads,
                    MODE_AUTOMATIC, iterations, warmup))
    else:
        raise MatrixError(f"unknown preset: {preset}")

    if pinned_cpu is not None:
        if pinned_cpu < 0:
            raise MatrixError("pinned CPU must be non-negative")
        if any(int(job["cpu_count"]) != 1 for job in jobs):
            raise MatrixError("--pinned-cpu is only valid for all-single-thread presets")
        for job in jobs:
            job.pop("cpu_count")
            job["cpu_set"] = [pinned_cpu]

    ids = [str(job["id"]) for job in jobs]
    if len(ids) != len(set(ids)):
        duplicates = sorted(key for key, count in Counter(ids).items() if count > 1)
        raise MatrixError("duplicate job ids: " + ", ".join(duplicates))
    spec = {
        "schema": SPEC_SCHEMA,
        "root": str(Path.cwd().resolve()),
        "base_seed": 0x4C454F32,
        "workers": workers,
        "defaults": {
            "cpu_policy": "physical-first",
            "timeout_seconds": 1800,
            "memory_mb": 4096,
        },
        "metadata": {
            "preset": preset,
            "serial_timing_jobs": workers == 1,
            "isolation_status": "not-established-by-generator",
            "benchmark": benchmark,
            "pinned_cpu": pinned_cpu,
        },
        "jobs": sorted(jobs, key=lambda item: str(item["id"])),
    }
    spec["spec_digest"] = digest(spec)
    return spec


def _positive_finite(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise MatrixError(f"{label} must be numeric")
    converted = float(value)
    if not math.isfinite(converted) or converted <= 0.0:
        raise MatrixError(f"{label} must be finite and positive")
    return converted


def _integer(value: object, label: str, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise MatrixError(f"{label} must be an integer")
    if value < 0 or (value == 0 and not allow_zero):
        qualifier = "non-negative" if allow_zero else "positive"
        raise MatrixError(f"{label} must be {qualifier}")
    return value


def _validate_benchmark(value: object, job_id: str) -> dict:
    if not isinstance(value, dict) or value.get("schema") != "leopard2-benchmark-v1":
        raise MatrixError(f"job {job_id} did not emit leopard2-benchmark-v1 JSON")
    for key in ("build", "parameters", "resolved", "correctness", "memory", "metrics", "legacy"):
        if not isinstance(value.get(key), dict):
            raise MatrixError(f"job {job_id} benchmark JSON is missing {key}")
    if value["correctness"].get("leopard2_round_trip") is not True:
        raise MatrixError(f"job {job_id} did not report a successful round trip")
    parameters = value["parameters"]
    for name in ("K", "R", "shard_bytes", "batch", "reuse", "iterations",
                 "thread_count"):
        _integer(parameters.get(name), f"job {job_id} parameter {name}")
    for name in ("loss_count", "warmup", "seed"):
        _integer(
            parameters.get(name), f"job {job_id} parameter {name}",
            allow_zero=True)
    for name in ("requested_profile", "requested_field", "requested_backend"):
        if not isinstance(parameters.get(name), str) or not parameters[name]:
            raise MatrixError(f"job {job_id} parameter {name} is missing")
    if not isinstance(parameters.get("force_generic_decode"), bool):
        raise MatrixError(f"job {job_id} has an invalid force_generic_decode value")
    if not isinstance(parameters.get("force_specialized_decode"), bool):
        raise MatrixError(f"job {job_id} has an invalid force_specialized_decode value")
    if (parameters["force_generic_decode"] and
            parameters["force_specialized_decode"]):
        raise MatrixError(f"job {job_id} forces two decoder policies")
    missing = parameters.get("missing_original_indices")
    if (not isinstance(missing, list) or
            not all(isinstance(index, int) and not isinstance(index, bool) and index >= 0
                    for index in missing) or
            len(missing) != len(set(missing)) or missing != sorted(missing)):
        raise MatrixError(f"job {job_id} has invalid missing original indices")
    if parameters.get("loss_count") != len(missing):
        raise MatrixError(f"job {job_id} loss count does not match missing indices")
    if (parameters["loss_count"] > min(parameters["K"], parameters["R"]) or
            any(index >= parameters["K"] for index in missing)):
        raise MatrixError(f"job {job_id} missing original indices are out of range")
    resolved = value["resolved"]
    for name in ("profile", "field", "backend"):
        if not isinstance(resolved.get(name), str) or not resolved[name]:
            raise MatrixError(f"job {job_id} resolved {name} is missing")
    for name in ("thread_count", "parent_count", "padded_side"):
        _integer(resolved.get(name), f"job {job_id} resolved {name}")
    metrics = value["metrics"].get("decode_execution")
    if not isinstance(metrics, dict):
        raise MatrixError(f"job {job_id} decode execution metrics are missing")
    _positive_finite(
        metrics.get("median_us_per_batch_call"),
        f"job {job_id} decode median")
    return value


def _validate_manifest(value: object) -> dict:
    if not isinstance(value, dict) or value.get("schema") != LAB_MANIFEST_SCHEMA:
        raise MatrixError("unsupported lab manifest")
    expected_digest = value.get("manifest_digest")
    unsigned = dict(value)
    unsigned.pop("manifest_digest", None)
    if expected_digest != digest(unsigned):
        raise MatrixError("lab manifest digest does not match its contents")
    source_spec = value.get("source_spec")
    if (not isinstance(source_spec, dict) or
            source_spec.get("schema") != SPEC_SCHEMA or
            not isinstance(source_spec.get("digest"), str) or
            len(source_spec["digest"]) != 64 or
            any(character not in "0123456789abcdef"
                for character in source_spec["digest"]) or
            not isinstance(source_spec.get("metadata"), dict)):
        raise MatrixError("lab manifest is missing source specification metadata")
    jobs = value.get("jobs")
    if not isinstance(jobs, list):
        raise MatrixError("lab manifest jobs are missing")
    ids = [job.get("id") for job in jobs if isinstance(job, dict)]
    if (len(ids) != len(jobs) or not all(isinstance(job_id, str) for job_id in ids) or
            ids != sorted(ids) or len(ids) != len(set(ids))):
        raise MatrixError("lab manifest jobs are invalid or duplicated")
    for job in jobs:
        unsigned_job = dict(job)
        expected_job_digest = unsigned_job.pop("job_digest", None)
        if expected_job_digest != digest(unsigned_job):
            raise MatrixError(f"job {job['id']} digest does not match its contents")
        if not isinstance(job.get("cpu_group"), str) or not job["cpu_group"]:
            raise MatrixError(f"job {job['id']} is missing its scheduled pair group")
        expected_cell = job.get("benchmark_cell")
        if (not isinstance(expected_cell, dict) or
                set(expected_cell) != set(EXPECTED_CELL_FIELDS)):
            raise MatrixError(
                f"job {job['id']} expected benchmark cell has the wrong fields")
    return value


def collect(manifest_path: Path, results_dir: Path) -> dict:
    manifest = _validate_manifest(load_json(manifest_path))
    records = []
    outcomes: Counter[str] = Counter()
    for job in manifest.get("jobs", []):
        if not isinstance(job, dict) or not isinstance(job.get("id"), str):
            raise MatrixError("manifest contains an invalid job")
        executable = job.get("executable")
        if (not isinstance(executable, dict) or
                not isinstance(executable.get("sha256"), str) or
                len(executable["sha256"]) != 64 or
                any(character not in "0123456789abcdef"
                    for character in executable["sha256"])):
            raise MatrixError(f"job {job['id']} is missing executable identity")
        job_id = job["id"]
        job_dir = results_dir / "jobs" / job_id
        lab_result = load_json(job_dir / "result.json")
        if (not isinstance(lab_result, dict) or
                lab_result.get("schema") != "leopard2-lab-result/v1" or
                lab_result.get("state") != "complete"):
            raise MatrixError(f"job {job_id} has an invalid lab result")
        unsigned_result = dict(lab_result)
        expected_result_digest = unsigned_result.pop("result_digest", None)
        if expected_result_digest != digest(unsigned_result):
            raise MatrixError(f"job {job_id} terminal result digest is invalid")
        outputs = lab_result.get("outputs")
        if (not isinstance(outputs, dict) or
                outputs.get("stdout") != file_identity(job_dir / "stdout.txt") or
                outputs.get("stderr") != file_identity(job_dir / "stderr.txt")):
            raise MatrixError(f"job {job_id} output content identity is invalid")
        if lab_result.get("job_digest") != job.get("job_digest"):
            raise MatrixError(f"job {job_id} result is stale")
        if lab_result.get("job_id") != job_id:
            raise MatrixError(f"job {job_id} result identifies a different job")
        result_cpu_set = lab_result.get("cpu_set")
        if result_cpu_set != job.get("cpu_set"):
            raise MatrixError(f"job {job_id} result CPU set differs from its manifest")
        outcome = str(lab_result.get("outcome"))
        outcomes[outcome] += 1
        record = {
            "job_id": job_id,
            "outcome": outcome,
            "exit_code": lab_result.get("exit_code"),
            "cpu_set": lab_result.get("cpu_set"),
            "duration_seconds": lab_result.get("duration_seconds"),
            "executable": executable,
            "cpu_group": job["cpu_group"],
        }
        if outcome == "success":
            benchmark = _validate_benchmark(load_json(job_dir / "stdout.txt"), job_id)
            expected_cell = job.get("benchmark_cell")
            actual_cell = {
                name: benchmark["parameters"].get(name)
                for name in EXPECTED_CELL_FIELDS}
            if actual_cell != expected_cell:
                mismatches = sorted(
                    name for name in EXPECTED_CELL_FIELDS
                    if actual_cell.get(name) != expected_cell.get(name))
                raise MatrixError(
                    f"job {job_id} benchmark parameters differ from its manifest: "
                    + ", ".join(mismatches))
            record["benchmark"] = benchmark
            record["build_digest"] = digest(benchmark["build"])
        else:
            record["stderr"] = (job_dir / "stderr.txt").read_text(
                encoding="utf-8", errors="replace")[-4096:]
        records.append(record)

    pair_map: Dict[str, Dict[str, dict]] = {}
    pair_parameters: Dict[str, dict] = {}
    for record in records:
        benchmark = record.get("benchmark")
        if not isinstance(benchmark, dict):
            continue
        parameters = benchmark["parameters"]
        resolved = benchmark["resolved"]
        identity = {
            "K": parameters.get("K"),
            "R": parameters.get("R"),
            "requested_profile": parameters.get("requested_profile"),
            "requested_field": parameters.get("requested_field"),
            "requested_backend": parameters.get("requested_backend"),
            "shard_bytes": parameters.get("shard_bytes"),
            "loss_count": parameters.get("loss_count"),
            "missing_original_indices": parameters.get("missing_original_indices"),
            "batch": parameters.get("batch"),
            "reuse": parameters.get("reuse"),
            "iterations": parameters.get("iterations"),
            "warmup": parameters.get("warmup"),
            "thread_count": parameters.get("thread_count"),
            "seed": parameters.get("seed"),
            "resolved_profile": resolved.get("profile"),
            "resolved_field": resolved.get("field"),
            "resolved_backend": resolved.get("backend"),
            "resolved_thread_count": resolved.get("thread_count"),
            "parent_count": resolved.get("parent_count"),
            "padded_side": resolved.get("padded_side"),
            "binary_sha256": record["executable"]["sha256"],
            "build_digest": record["build_digest"],
            "cpu_set": record["cpu_set"],
        }
        key = record["cpu_group"]
        if parameters["force_generic_decode"]:
            mode = MODE_FORCED_GENERIC
        elif parameters["force_specialized_decode"]:
            mode = MODE_FORCED_SPECIALIZED
        else:
            mode = MODE_AUTOMATIC
        pair = pair_map.setdefault(key, {})
        previous_identity = pair_parameters.get(key)
        if previous_identity is not None and identity != previous_identity:
            mismatches = sorted(
                name for name in identity
                if identity.get(name) != previous_identity.get(name))
            raise MatrixError(
                f"scheduled pair group {key} emitted different identities: "
                + ", ".join(mismatches))
        if mode in pair:
            raise MatrixError(
                f"duplicate {mode} pair member for jobs "
                f"{pair[mode]['job_id']} and {record['job_id']}")
        pair[mode] = record
        pair_parameters[key] = identity
    comparisons = []
    dispatcher_checks = []
    for key, pair in sorted(pair_map.items()):
        if MODE_FORCED_SPECIALIZED not in pair or MODE_FORCED_GENERIC not in pair:
            continue
        forced_specialized = pair[MODE_FORCED_SPECIALIZED][
            "benchmark"]["metrics"]["decode_execution"]
        forced_generic = pair[MODE_FORCED_GENERIC][
            "benchmark"]["metrics"]["decode_execution"]
        forced_specialized_us = _positive_finite(
            forced_specialized.get("median_us_per_batch_call"),
            f"job {pair[MODE_FORCED_SPECIALIZED]['job_id']} decode median")
        forced_generic_us = _positive_finite(
            forced_generic.get("median_us_per_batch_call"),
            f"job {pair[MODE_FORCED_GENERIC]['job_id']} decode median")
        comparisons.append({
            "parameters": pair_parameters[key],
            "forced_specialized_job": pair[MODE_FORCED_SPECIALIZED]["job_id"],
            "forced_generic_job": pair[MODE_FORCED_GENERIC]["job_id"],
            "forced_specialized_median_us": forced_specialized_us,
            "forced_generic_median_us": forced_generic_us,
            "forced_specialized_speedup_vs_forced_generic": round(
                forced_generic_us / forced_specialized_us, 9),
        })
        if MODE_AUTOMATIC in pair:
            automatic = pair[MODE_AUTOMATIC][
                "benchmark"]["metrics"]["decode_execution"]
            automatic_us = _positive_finite(
                automatic.get("median_us_per_batch_call"),
                f"job {pair[MODE_AUTOMATIC]['job_id']} decode median")
            best_forced_mode = (MODE_FORCED_SPECIALIZED
                                if forced_specialized_us <= forced_generic_us
                                else MODE_FORCED_GENERIC)
            best_forced_us = min(forced_specialized_us, forced_generic_us)
            dispatcher_checks.append({
                "parameters": pair_parameters[key],
                "automatic_job": pair[MODE_AUTOMATIC]["job_id"],
                "forced_specialized_job": pair[MODE_FORCED_SPECIALIZED]["job_id"],
                "forced_generic_job": pair[MODE_FORCED_GENERIC]["job_id"],
                "automatic_median_us": automatic_us,
                "best_forced_mode_by_median": best_forced_mode,
                "best_forced_median_us": best_forced_us,
                "automatic_ratio_to_best_forced": round(
                    automatic_us / best_forced_us, 9),
            })
    output = {
        "schema": SCHEMA,
        "manifest_digest": manifest.get("manifest_digest"),
        "source_spec": manifest.get("source_spec"),
        "summary": dict(sorted(outcomes.items())),
        "record_count": len(records),
        "comparison_count": len(comparisons),
        "dispatcher_check_count": len(dispatcher_checks),
        "comparisons": comparisons,
        "dispatcher_checks": dispatcher_checks,
        "records": records,
    }
    output["content_digest"] = digest(output)
    return output


def self_test() -> None:
    smoke = make_spec("/tmp/bench_leopard2", "smoke", 1, 3, 1)
    smoke_groups = Counter(str(job.get("cpu_group")) for job in smoke["jobs"])
    if sorted(smoke_groups.values()) != [3, 3]:
        raise MatrixError("paired jobs do not share stable CPU-assignment groups")
    first = make_spec("/tmp/bench_leopard2", "checkpoint", 1, 3, 1)
    second = make_spec("/tmp/bench_leopard2", "checkpoint", 1, 3, 1)
    if canonical_bytes(first) != canonical_bytes(second):
        raise MatrixError("spec generation is not deterministic")
    if len(first["jobs"]) != 63:
        raise MatrixError("checkpoint job-count invariant failed")
    pinned = make_spec("/tmp/bench_leopard2", "checkpoint", 1, 3, 1, 7)
    if any(job.get("cpu_set") != [7] or "cpu_count" in job for job in pinned["jobs"]):
        raise MatrixError("pinned checkpoint invariant failed")
    required = make_spec("/tmp/bench_leopard2", "required", 1, 3, 1)
    if len(required["jobs"]) != 2483:
        raise MatrixError("required preset job-count invariant failed")
    crossover = make_spec("/tmp/bench_leopard2", "balanced-crossover", 1, 3, 1, 7)
    if len(crossover["jobs"]) != 216:
        raise MatrixError("balanced crossover job-count invariant failed")
    for k, r in ((1, 1), (8, 248), (240, 16)):
        losses = loss_counts(k, r)
        if losses != tuple(sorted(set(losses))) or losses[0] != 0 or losses[-1] != min(k, r):
            raise MatrixError("loss-count boundary invariant failed")
    if memory_limit_mb(16, 240, 65536, 1) < 2048:
        raise MatrixError("memory-limit floor invariant failed")
    if memory_limit_mb(4096, 512, 16777216, 1) <= 65536:
        raise MatrixError("large-cell memory limit is still incorrectly capped")
    scaling_jobs = [job for job in required["jobs"] if "-scaling." in str(job["id"])]
    if len(scaling_jobs) != 24:
        raise MatrixError("thread-scaling job-count invariant failed")
    scaling_work = set()
    for job in scaling_jobs:
        command = list(job["command"])
        scaling_work.add((
            command[command.index("--bytes") + 1],
            command[command.index("--loss") + 1],
            command[command.index("--batch") + 1],
            command[command.index("--reuse") + 1],
        ))
    if scaling_work != {("4096", "8", "128", "8")}:
        raise MatrixError("thread scaling changes the amount of measured work")

    def expect_error(label: str, callback) -> None:
        try:
            callback()
        except MatrixError:
            return
        raise MatrixError(f"negative collector test did not reject {label}")

    with tempfile.TemporaryDirectory(prefix="leopard2-benchmark-matrix-") as directory:
        root = Path(directory)
        executable = {"path": "/tmp/bench_leopard2", "sha256": "a" * 64,
                      "size_bytes": 12345}

        def write_fixture(name: str, mutation=None, duplicate=False):
            case = root / name
            jobs = []
            lab_results = {}
            benchmarks = {}
            modes = [
                (MODE_AUTOMATIC, 12.0),
                (MODE_FORCED_SPECIALIZED, 10.0),
                (MODE_FORCED_GENERIC, 20.0),
            ]
            if duplicate:
                modes.append((MODE_FORCED_SPECIALIZED, 11.0))
            for ordinal, (mode, median) in enumerate(modes):
                forced_generic = mode == MODE_FORCED_GENERIC
                forced_specialized = mode == MODE_FORCED_SPECIALIZED
                job_id = f"pair.{mode}.{ordinal}"
                job = {
                    "id": job_id,
                    "cpu_set": [7],
                    "cpu_group": "self-test-pair",
                    "executable": executable,
                    "command": (["/tmp/bench_leopard2", "--force-generic"]
                                if forced_generic else
                                (["/tmp/bench_leopard2", "--force-specialized"]
                                 if forced_specialized else ["/tmp/bench_leopard2"])),
                }
                job["job_digest"] = digest(job)
                jobs.append(job)
                lab_results[job_id] = {
                    "schema": "leopard2-lab-result/v1",
                    "state": "complete",
                    "job_id": job_id,
                    "job_digest": job["job_digest"],
                    "outcome": "success",
                    "exit_code": 0,
                    "cpu_set": [7],
                    "duration_seconds": 0.1,
                }
                benchmarks[job_id] = {
                    "schema": "leopard2-benchmark-v1",
                    "build": {"compiler": "self-test", "cplusplus": 201103},
                    "parameters": {
                        "K": 240, "R": 16,
                        "requested_profile": "legacy_high_v1",
                        "requested_field": "auto", "requested_backend": "auto",
                        "force_generic_decode": forced_generic,
                        "force_specialized_decode": forced_specialized,
                        "shard_bytes": 4096, "loss_count": 2,
                        "missing_original_indices": [3, 11],
                        "batch": 1, "reuse": 8, "iterations": 3,
                        "warmup": 1, "thread_count": 1, "seed": 17,
                    },
                    "resolved": {
                        "profile": "legacy_high_v1", "field": "gf8",
                        "backend": "avx2", "thread_count": 1,
                        "parent_count": 256, "padded_side": 16,
                    },
                    "correctness": {"leopard2_round_trip": True},
                    "memory": {},
                    "metrics": {
                        "decode_execution": {"median_us_per_batch_call": median}},
                    "legacy": {},
                }
                job["benchmark_cell"] = {
                    name: value for name, value in
                    benchmarks[job_id]["parameters"].items()
                    if name != "missing_original_indices"}
                job["job_digest"] = digest({
                    key: value for key, value in job.items()
                    if key != "job_digest"})
                lab_results[job_id]["job_digest"] = job["job_digest"]
            jobs.sort(key=lambda job: str(job["id"]))
            manifest = {
                "schema": LAB_MANIFEST_SCHEMA,
                "root": "/tmp",
                "source_spec": {
                    "schema": SPEC_SCHEMA, "digest": "b" * 64,
                    "metadata": {"preset": "self-test"}},
                "jobs": jobs,
            }
            if mutation:
                mutation(manifest, lab_results, benchmarks)
            manifest["manifest_digest"] = digest(manifest)
            write_json(case / "manifest.json", manifest)
            for job in jobs:
                job_id = str(job["id"])
                job_dir = case / "results" / "jobs" / job_id
                write_json(job_dir / "stdout.txt", benchmarks[job_id])
                (job_dir / "stderr.txt").write_text("", encoding="utf-8")
                lab_results[job_id]["outputs"] = {
                    "stdout": file_identity(job_dir / "stdout.txt"),
                    "stderr": file_identity(job_dir / "stderr.txt"),
                }
                lab_results[job_id]["result_digest"] = digest(lab_results[job_id])
                write_json(job_dir / "result.json", lab_results[job_id])
            return case / "manifest.json", case / "results"

        manifest_path, results_path = write_fixture("valid")
        collected = collect(manifest_path, results_path)
        if (collected["record_count"] != 3 or collected["comparison_count"] != 1 or
                collected["dispatcher_check_count"] != 1 or
                collected["comparisons"][0][
                    "forced_specialized_speedup_vs_forced_generic"] != 2.0 or
                collected["dispatcher_checks"][0][
                    "automatic_ratio_to_best_forced"] != 1.2 or
                collected["source_spec"]["metadata"]["preset"] != "self-test"):
            raise MatrixError("collector pairing invariant failed")

        bad_schema = write_fixture(
            "bad-schema",
            lambda manifest, results, benchmarks: benchmarks[
                sorted(benchmarks)[0]].update(schema="not-a-benchmark"))
        expect_error("benchmark schema", lambda: collect(*bad_schema))

        bad_spec_schema = write_fixture(
            "bad-spec-schema",
            lambda manifest, results, benchmarks: manifest[
                "source_spec"].update(schema="leopard2-benchmark-spec/v999"))
        expect_error("source specification schema", lambda: collect(*bad_spec_schema))

        stale_result = write_fixture(
            "stale-result",
            lambda manifest, results, benchmarks: results[
                sorted(results)[0]].update(job_digest="stale"))
        expect_error("stale result digest", lambda: collect(*stale_result))

        failed_round_trip = write_fixture(
            "failed-round-trip",
            lambda manifest, results, benchmarks: benchmarks[
                sorted(benchmarks)[0]]["correctness"].update(
                    leopard2_round_trip=False))
        expect_error("failed round trip", lambda: collect(*failed_round_trip))

        policy_mismatch = write_fixture(
            "policy-mismatch",
            lambda manifest, results, benchmarks: benchmarks[
                sorted(benchmarks)[0]]["parameters"].update(
                    force_generic_decode=True))
        expect_error("command/result policy mismatch", lambda: collect(*policy_mismatch))

        missing_pattern_mismatch = write_fixture(
            "missing-pattern-mismatch",
            lambda manifest, results, benchmarks: benchmarks[
                sorted(benchmarks)[0]]["parameters"].update(
                    missing_original_indices=[4, 12]))
        expect_error(
            "scheduled pair missing-index mismatch",
            lambda: collect(*missing_pattern_mismatch))

        resolved_identity_mismatch = write_fixture(
            "resolved-identity-mismatch",
            lambda manifest, results, benchmarks: benchmarks[
                sorted(benchmarks)[0]]["resolved"].update(backend="scalar"))
        expect_error(
            "scheduled pair resolved-identity mismatch",
            lambda: collect(*resolved_identity_mismatch))

        def remove_expected_cell_field(manifest, results, benchmarks):
            job = manifest["jobs"][0]
            job["benchmark_cell"].pop("seed")
            unsigned_job = dict(job)
            unsigned_job.pop("job_digest", None)
            job["job_digest"] = digest(unsigned_job)
            results[job["id"]]["job_digest"] = job["job_digest"]

        incomplete_expected_cell = write_fixture(
            "incomplete-expected-cell", remove_expected_cell_field)
        expect_error(
            "incomplete expected benchmark cell",
            lambda: collect(*incomplete_expected_cell))

        for parameter_name, replacement in (("K", 239), ("seed", 18),
                                             ("shard_bytes", 8192)):
            parameter_mismatch = write_fixture(
                "parameter-mismatch-" + parameter_name,
                lambda manifest, results, benchmarks, name=parameter_name,
                       value=replacement: benchmarks[sorted(benchmarks)[0]][
                           "parameters"].update({name: value}))
            expect_error(
                parameter_name + " parameter mismatch",
                lambda fixture=parameter_mismatch: collect(*fixture))

        zero_timing = write_fixture(
            "zero-timing",
            lambda manifest, results, benchmarks: benchmarks[
                sorted(benchmarks)[0]]["metrics"]["decode_execution"].update(
                    median_us_per_batch_call=0.0))
        expect_error("non-positive timing", lambda: collect(*zero_timing))

        nonfinite_timing = write_fixture(
            "nonfinite-timing",
            lambda manifest, results, benchmarks: benchmarks[
                sorted(benchmarks)[0]]["metrics"]["decode_execution"].update(
                    median_us_per_batch_call=float("nan")))
        expect_error("non-finite timing", lambda: collect(*nonfinite_timing))

        duplicate_pair = write_fixture("duplicate-pair", duplicate=True)
        expect_error("duplicate pair member", lambda: collect(*duplicate_pair))

        corrupted_output = write_fixture("corrupted-output")
        corrupted_job = sorted(
            (corrupted_output[1] / "jobs").iterdir(), key=lambda path: path.name)[0]
        with (corrupted_job / "stdout.txt").open("a", encoding="utf-8") as output:
            output.write("corruption\n")
        expect_error("post-run output corruption", lambda: collect(*corrupted_output))

        manifest_path, results_path = write_fixture("stale-manifest")
        stale_manifest = load_json(manifest_path)
        stale_manifest["source_spec"]["metadata"]["preset"] = "changed"
        write_json(manifest_path, stale_manifest)
        expect_error("stale manifest digest", lambda: collect(manifest_path, results_path))
    print(json.dumps({
        "balanced_crossover_jobs": len(crossover["jobs"]),
        "checkpoint_jobs": len(first["jobs"]),
        "required_jobs": len(required["jobs"]),
        "status": "PASS",
    }, sort_keys=True))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate", help="write a leopard2_lab job specification")
    generate.add_argument("--benchmark", required=True)
    generate.add_argument(
        "--preset",
        choices=("smoke", "checkpoint", "balanced-crossover", "required"),
        default="checkpoint")
    generate.add_argument("--workers", type=int, default=1)
    generate.add_argument("--iterations", type=int, default=9)
    generate.add_argument("--warmup", type=int, default=2)
    generate.add_argument("--pinned-cpu", type=int)
    generate.add_argument("--output", required=True)
    collect_parser = subparsers.add_parser("collect", help="collect benchmark JSON from lab results")
    collect_parser.add_argument("--manifest", required=True)
    collect_parser.add_argument("--results-dir", required=True)
    collect_parser.add_argument("--output", required=True)
    subparsers.add_parser("self-test")
    args = parser.parse_args(argv)
    try:
        if args.command == "generate":
            spec = make_spec(
                args.benchmark, args.preset, args.workers, args.iterations,
                args.warmup, args.pinned_cpu)
            write_json(Path(args.output), spec)
            print(json.dumps({"jobs": len(spec["jobs"]), "output": args.output,
                              "spec_digest": spec["spec_digest"]}, sort_keys=True))
        elif args.command == "collect":
            output = collect(Path(args.manifest), Path(args.results_dir))
            write_json(Path(args.output), output)
            print(json.dumps({"records": output["record_count"],
                              "comparisons": output["comparison_count"],
                              "dispatcher_checks": output["dispatcher_check_count"],
                              "output": args.output}, sort_keys=True))
        else:
            self_test()
    except MatrixError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
