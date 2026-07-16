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
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


SCHEMA = "leopard2-benchmark-matrix/v1"
SPEC_SCHEMA = "leopard2-benchmark-spec/v1"

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
    return max(2048, min(65536, (estimate + (1 << 20) - 1) >> 20))


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
    force_generic: bool,
    iterations: int,
    warmup: int,
) -> Dict[str, object]:
    mode = "generic" if force_generic else "specialized"
    job_id = (
        f"{category}.k{k}.r{r}.{profile}.b{shard_bytes}.l{losses}."
        f"batch{batch}.reuse{reuse}.t{threads}.{mode}"
    )
    pair_identity = (
        f"{category}|{k}|{r}|{profile}|{shard_bytes}|{losses}|"
        f"{batch}|{reuse}|{threads}"
    ).encode("ascii")
    pair_seed = int.from_bytes(hashlib.sha256(pair_identity).digest()[:4], "big")
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
    if force_generic:
        command.append("--force-generic")
    return {
        "id": job_id,
        "command": command,
        "cpu_count": threads,
        "memory_mb": memory_limit_mb(k, r, shard_bytes, batch),
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
            for force_generic in (False, True):
                jobs.append(_job(
                    benchmark, category, k, r, profile, shard_bytes, losses,
                    1, 2, 1, force_generic, iterations, warmup))
    elif preset == "checkpoint":
        for category, k, r, profile in _checkpoint_cells():
            for shard_bytes in (4096, 65536, 1048576):
                for losses in (0, 1, min(8, k, r)):
                    jobs.append(_job(
                        benchmark, category, k, r, profile, shard_bytes,
                        losses, 1, 8, 1, False, iterations, warmup))
                    if losses:
                        jobs.append(_job(
                            benchmark, category, k, r, profile, shard_bytes,
                            losses, 1, 8, 1, True, iterations, warmup))
    elif preset == "balanced-crossover":
        for shard_bytes in (256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536):
            for losses in (1, 2, 4, 8, 16, 32, 64, 128):
                for force_generic in (False, True):
                    jobs.append(_job(
                        benchmark, "balanced-crossover", 128, 128, "high",
                        shard_bytes, losses, 1, 16, 1, force_generic,
                        iterations, warmup))
    elif preset == "required":
        for category in sorted(RATE_CASES):
            for k, r, profile in RATE_CASES[category]:
                for shard_bytes in FULL_SHARD_BYTES:
                    for losses in loss_counts(k, r):
                        for force_generic in ((False, True) if losses else (False,)):
                            jobs.append(_job(
                                benchmark, category, k, r, profile, shard_bytes,
                                losses, 1, 8, 1, force_generic, iterations, warmup))
        # Reuse/batch and thread scaling are separate from the count/size/loss
        # Cartesian product to avoid impossible multi-terabyte allocations.
        for category, k, r, profile in _checkpoint_cells():
            for batch, reuse in ((8, 8), (64, 64), (1024, 1024)):
                shard_bytes = 65536 if batch <= 8 else (4096 if batch <= 64 else 256)
                jobs.append(_job(
                    benchmark, category + "-reuse", k, r, profile,
                    shard_bytes, min(8, k, r), batch, reuse, 1, False,
                    iterations, warmup))
            for threads in (1, 2, 4, 8, 16, 32, 64, 128):
                jobs.append(_job(
                    benchmark, category + "-scaling", k, r, profile,
                    4096, min(8, k, r), max(threads, 8), 8, threads,
                    False, iterations, warmup))
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
            "authoritative_isolation": workers == 1,
            "benchmark": benchmark,
            "pinned_cpu": pinned_cpu,
        },
        "jobs": sorted(jobs, key=lambda item: str(item["id"])),
    }
    spec["spec_digest"] = digest(spec)
    return spec


def _validate_benchmark(value: object, job_id: str) -> dict:
    if not isinstance(value, dict) or value.get("schema") != "leopard2-benchmark-v1":
        raise MatrixError(f"job {job_id} did not emit leopard2-benchmark-v1 JSON")
    for key in ("parameters", "resolved", "correctness", "memory", "metrics", "legacy"):
        if not isinstance(value.get(key), dict):
            raise MatrixError(f"job {job_id} benchmark JSON is missing {key}")
    if value["correctness"].get("leopard2_round_trip") is not True:
        raise MatrixError(f"job {job_id} did not report a successful round trip")
    return value


def collect(manifest_path: Path, results_dir: Path) -> dict:
    manifest = load_json(manifest_path)
    if not isinstance(manifest, dict) or manifest.get("schema") != "leopard2-lab-manifest/v1":
        raise MatrixError("unsupported lab manifest")
    records = []
    outcomes: Counter[str] = Counter()
    for job in manifest.get("jobs", []):
        if not isinstance(job, dict) or not isinstance(job.get("id"), str):
            raise MatrixError("manifest contains an invalid job")
        job_id = job["id"]
        job_dir = results_dir / "jobs" / job_id
        lab_result = load_json(job_dir / "result.json")
        if not isinstance(lab_result, dict) or lab_result.get("schema") != "leopard2-lab-result/v1":
            raise MatrixError(f"job {job_id} has an invalid lab result")
        if lab_result.get("job_digest") != job.get("job_digest"):
            raise MatrixError(f"job {job_id} result is stale")
        outcome = str(lab_result.get("outcome"))
        outcomes[outcome] += 1
        record = {
            "job_id": job_id,
            "outcome": outcome,
            "exit_code": lab_result.get("exit_code"),
            "cpu_set": lab_result.get("cpu_set"),
            "duration_seconds": lab_result.get("duration_seconds"),
        }
        if outcome == "success":
            benchmark = _validate_benchmark(load_json(job_dir / "stdout.txt"), job_id)
            record["benchmark"] = benchmark
        else:
            record["stderr"] = (job_dir / "stderr.txt").read_text(
                encoding="utf-8", errors="replace")[-4096:]
        records.append(record)

    pair_map: Dict[Tuple[object, ...], Dict[bool, dict]] = {}
    for record in records:
        benchmark = record.get("benchmark")
        if not isinstance(benchmark, dict):
            continue
        parameters = benchmark["parameters"]
        key = tuple(parameters.get(name) for name in (
            "K", "R", "requested_profile", "requested_field",
            "requested_backend", "shard_bytes", "loss_count", "batch",
            "reuse", "thread_count", "seed"))
        pair_map.setdefault(key, {})[bool(parameters.get("force_generic_decode"))] = record
    comparisons = []
    for key, pair in sorted(pair_map.items(), key=lambda item: repr(item[0])):
        if False not in pair or True not in pair:
            continue
        specialized = pair[False]["benchmark"]["metrics"]["decode_execution"]
        generic = pair[True]["benchmark"]["metrics"]["decode_execution"]
        specialized_us = float(specialized["median_us_per_batch_call"])
        generic_us = float(generic["median_us_per_batch_call"])
        comparisons.append({
            "parameters": list(key),
            "specialized_job": pair[False]["job_id"],
            "generic_job": pair[True]["job_id"],
            "specialized_median_us": specialized_us,
            "generic_median_us": generic_us,
            "specialized_speedup_vs_generic": round(generic_us / specialized_us, 9),
        })
    output = {
        "schema": SCHEMA,
        "manifest_digest": manifest.get("manifest_digest"),
        "summary": dict(sorted(outcomes.items())),
        "record_count": len(records),
        "comparison_count": len(comparisons),
        "comparisons": comparisons,
        "records": records,
    }
    output["content_digest"] = digest(output)
    return output


def self_test() -> None:
    first = make_spec("/tmp/bench_leopard2", "checkpoint", 1, 3, 1)
    second = make_spec("/tmp/bench_leopard2", "checkpoint", 1, 3, 1)
    if canonical_bytes(first) != canonical_bytes(second):
        raise MatrixError("spec generation is not deterministic")
    if len(first["jobs"]) != 45:
        raise MatrixError("checkpoint job-count invariant failed")
    pinned = make_spec("/tmp/bench_leopard2", "checkpoint", 1, 3, 1, 7)
    if any(job.get("cpu_set") != [7] or "cpu_count" in job for job in pinned["jobs"]):
        raise MatrixError("pinned checkpoint invariant failed")
    required = make_spec("/tmp/bench_leopard2", "required", 1, 3, 1)
    if len(required["jobs"]) < 1000:
        raise MatrixError("required preset unexpectedly omitted the broad matrix")
    crossover = make_spec("/tmp/bench_leopard2", "balanced-crossover", 1, 3, 1, 7)
    if len(crossover["jobs"]) != 144:
        raise MatrixError("balanced crossover job-count invariant failed")
    for k, r in ((1, 1), (8, 248), (240, 16)):
        losses = loss_counts(k, r)
        if losses != tuple(sorted(set(losses))) or losses[0] != 0 or losses[-1] != min(k, r):
            raise MatrixError("loss-count boundary invariant failed")
    if memory_limit_mb(16, 240, 65536, 1) < 2048:
        raise MatrixError("memory-limit floor invariant failed")
    with tempfile.TemporaryDirectory(prefix="leopard2-benchmark-matrix-") as directory:
        root = Path(directory)
        jobs = []
        for generic, median in ((False, 10.0), (True, 20.0)):
            job_id = "pair.generic" if generic else "pair.specialized"
            job = {"id": job_id, "job_digest": "digest-" + job_id}
            jobs.append(job)
            job_dir = root / "results" / "jobs" / job_id
            write_json(job_dir / "result.json", {
                "schema": "leopard2-lab-result/v1",
                "job_digest": job["job_digest"],
                "outcome": "success",
                "exit_code": 0,
                "cpu_set": [0],
                "duration_seconds": 0.1,
            })
            write_json(job_dir / "stdout.txt", {
                "schema": "leopard2-benchmark-v1",
                "parameters": {
                    "K": 240, "R": 16, "requested_profile": "legacy_high_v1",
                    "requested_field": "auto", "requested_backend": "auto",
                    "force_generic_decode": generic, "shard_bytes": 4096,
                    "loss_count": 2, "batch": 1, "reuse": 8,
                    "thread_count": 1, "seed": 17,
                },
                "resolved": {},
                "correctness": {"leopard2_round_trip": True},
                "memory": {},
                "metrics": {"decode_execution": {"median_us_per_batch_call": median}},
                "legacy": {},
            })
        write_json(root / "manifest.json", {
            "schema": "leopard2-lab-manifest/v1",
            "manifest_digest": "self-test",
            "jobs": jobs,
        })
        collected = collect(root / "manifest.json", root / "results")
        if (collected["record_count"] != 2 or collected["comparison_count"] != 1 or
                collected["comparisons"][0]["specialized_speedup_vs_generic"] != 2.0):
            raise MatrixError("collector pairing invariant failed")
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
                              "output": args.output}, sort_keys=True))
        else:
            self_test()
    except MatrixError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
