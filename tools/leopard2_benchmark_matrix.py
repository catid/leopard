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
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

_TOOLS_DIRECTORY = str(Path(__file__).resolve().parent)
if _TOOLS_DIRECTORY not in sys.path:
    sys.path.insert(0, _TOOLS_DIRECTORY)

from leopard2_perf_evidence import (  # noqa: E402
    perf_probe_command,
    probe_command_matches_request,
    read_perf_stat_evidence,
)


SCHEMA = "leopard2-benchmark-matrix/v3"
SPEC_SCHEMA = "leopard2-benchmark-spec/v3"
LAB_MANIFEST_SCHEMA = "leopard2-lab-manifest/v3"
MODE_AUTOMATIC = "automatic"
MODE_FORCED_SPECIALIZED = "forced-specialized"
MODE_FORCED_GENERIC = "forced-generic"
ORDER_SINGLE = "single"
ORDER_AB = "ab"
ORDER_BA = "ba"
ORDER_LAYOUT = {
    ORDER_SINGLE: (MODE_AUTOMATIC,),
    # A is the specialized decoder and B is the generic decoder.  Automatic
    # follows the first pair so every nonzero-loss cell has one dispatcher
    # observation without changing the forced-path counterbalance.
    ORDER_AB: (MODE_FORCED_SPECIALIZED, MODE_FORCED_GENERIC, MODE_AUTOMATIC),
    ORDER_BA: (MODE_FORCED_GENERIC, MODE_FORCED_SPECIALIZED),
}
EXPECTED_CELL_FIELDS = (
    "K",
    "R",
    "requested_profile",
    "requested_field",
    "requested_backend",
    "force_generic_decode",
    "force_specialized_decode",
    "force_tiled_decode",
    "force_materialized_decode",
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
    # Larger analogues of the 1:3, 100:156, and 127:129 GF8 low-rate cells.
    # Their low-profile parents are 2048 or 8192 coordinates, so AUTO must
    # resolve them to GF16 rather than accidentally reusing a GF8-only case.
    "low-gf16": (
        (512, 1536, "low"),
        (1600, 2496, "low"),
        (2032, 2064, "low"),
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
GF16_LOW_RESOLUTION = {
    (512, 1536): (512, 2048),
    (1600, 2496): (2048, 8192),
    (2032, 2064): (2048, 8192),
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
DEFAULT_PERF_EVENTS = (
    "cycles",
    "instructions",
    "cache-references",
    "cache-misses",
    "branches",
    "branch-misses",
    "page-faults",
    "context-switches",
    "cpu-migrations",
    "dTLB-loads",
    "dTLB-load-misses",
)
PERF_EVENT_RE = re.compile(r"^[A-Za-z0-9_.:/=-]+$")


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


def _ordered_runs(losses: int) -> Tuple[Tuple[str, int, str], ...]:
    if losses == 0:
        return ((ORDER_SINGLE, 0, MODE_AUTOMATIC),)
    return tuple(
        (trial, sequence, mode)
        for trial in (ORDER_AB, ORDER_BA)
        for sequence, mode in enumerate(ORDER_LAYOUT[trial]))


def _parse_ordered_job_id(job_id: str) -> Tuple[str, int, str]:
    match = re.search(
        r"\.order-(single|ab|ba)\.slot([0-9]{2})-([a-z-]+)$", job_id)
    if match is None:
        raise MatrixError(f"job {job_id} has no valid counterbalance suffix")
    trial, sequence_text, mode = match.groups()
    sequence = int(sequence_text)
    layout = ORDER_LAYOUT[trial]
    if sequence >= len(layout) or layout[sequence] != mode:
        raise MatrixError(
            f"job {job_id} has an invalid {trial.upper()} mode sequence")
    return trial, sequence, mode


def _mode_from_cell(cell: Mapping[str, object]) -> str:
    force_generic = cell.get("force_generic_decode")
    force_specialized = cell.get("force_specialized_decode")
    if not isinstance(force_generic, bool) or not isinstance(force_specialized, bool):
        raise MatrixError("expected cell has invalid decoder policy flags")
    if force_generic and force_specialized:
        raise MatrixError("expected cell forces two decoder policies")
    if force_generic:
        return MODE_FORCED_GENERIC
    if force_specialized:
        return MODE_FORCED_SPECIALIZED
    return MODE_AUTOMATIC


def _validate_profile_resolution(
    job_id: str, cell: Mapping[str, object], benchmark: Mapping[str, object],
) -> None:
    if not job_id.startswith("low-gf16."):
        return
    counts = (cell.get("K"), cell.get("R"))
    expectation = GF16_LOW_RESOLUTION.get(counts)
    if expectation is None:
        raise MatrixError(f"job {job_id} is not a declared GF16 low-rate analogue")
    padded_side, parent_count = expectation
    resolved = benchmark["resolved"]
    if (resolved.get("profile") != "low_v1" or
            resolved.get("field") != "gf16" or
            resolved.get("padded_side") != padded_side or
            resolved.get("parent_count") != parent_count):
        raise MatrixError(
            f"job {job_id} did not resolve to its declared GF16 low parent")


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
    order_trial: str,
    order_sequence: int,
    iterations: int,
    warmup: int,
) -> Dict[str, object]:
    if mode not in (MODE_AUTOMATIC, MODE_FORCED_SPECIALIZED, MODE_FORCED_GENERIC):
        raise MatrixError(f"invalid decoder request mode: {mode}")
    layout = ORDER_LAYOUT.get(order_trial)
    if (layout is None or order_sequence < 0 or order_sequence >= len(layout) or
            layout[order_sequence] != mode):
        raise MatrixError("decoder mode does not match its counterbalance position")
    job_id = (
        f"{category}.k{k}.r{r}.{profile}.b{shard_bytes}.l{losses}."
        f"batch{batch}.reuse{reuse}.t{threads}.order-{order_trial}."
        f"slot{order_sequence:02d}-{mode}"
    )
    # All AB/BA repetitions of one cell intentionally share both the random
    # missing pattern and lab CPU-assignment group.  Collection splits the
    # shared scheduling group by the order suffix, so repeated modes do not
    # collide while both repetitions still run on the same assigned CPU set.
    cell_identity = (
        f"{category}|{k}|{r}|{profile}|{shard_bytes}|{losses}|"
        f"{batch}|{reuse}|{threads}"
    ).encode("ascii")
    pair_seed = int.from_bytes(hashlib.sha256(cell_identity).digest()[:4], "big")
    cpu_group = "pair-" + hashlib.sha256(cell_identity).hexdigest()
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
        "force_tiled_decode": False,
        "force_materialized_decode": False,
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
        # Counterbalanced observations are one atomic evidence unit.  The lab
        # runner resumes the whole logical cell only when all members were
        # completed by the same runner invocation.
        "resume_group": cpu_group,
        "benchmark_cell": benchmark_cell,
        "memory_mb": memory_mb,
        "minimum_memory_mb": memory_mb,
        "timeout_seconds": 1800,
    }


def _checkpoint_cells() -> Iterable[Tuple[str, int, int, str]]:
    yield "low", 16, 240, "low"
    yield "balanced", 128, 128, "high"
    yield "high", 240, 16, "high"


def _dimension_signature(job: Mapping[str, object]) -> Tuple[object, ...]:
    job_id = job.get("id")
    cell = job.get("benchmark_cell")
    if not isinstance(job_id, str) or ".k" not in job_id:
        raise MatrixError("matrix job is missing its category prefix")
    if not isinstance(cell, dict) or set(cell) != set(EXPECTED_CELL_FIELDS):
        raise MatrixError(f"job {job_id} has an incomplete expected cell")
    trial, sequence, declared_mode = _parse_ordered_job_id(job_id)
    mode = _mode_from_cell(cell)
    if mode != declared_mode:
        raise MatrixError(f"job {job_id} suffix disagrees with its expected cell")
    return (
        job_id.split(".k", 1)[0],
        cell["K"], cell["R"], cell["requested_profile"],
        cell["requested_field"], cell["requested_backend"],
        cell["shard_bytes"], cell["loss_count"], cell["batch"],
        cell["reuse"], cell["iterations"], cell["warmup"],
        cell["thread_count"], trial, sequence, mode,
    )


def _expected_required_dimensions(
    iterations: int, warmup: int,
) -> Counter[Tuple[object, ...]]:
    """Independently enumerate the required preset's exact dimension set.

    This intentionally does not call RATE_CASES, loss_counts, _ordered_runs,
    _job, or the production shard-size selector.  It is a compact oracle that
    catches missing dimensions and same-cardinality mode/order swaps.
    """
    expected_rates = {
        "balanced": ((128, 128, "high"), (128, 128, "low"),
                     (256, 256, "high")),
        "high": ((192, 64, "high"), (224, 32, "high"),
                 (240, 16, "high"), (248, 8, "high"),
                 (1000, 200, "high"), (4096, 512, "high")),
        "low": ((8, 248, "low"), (16, 240, "low"),
                (32, 224, "low"), (64, 192, "low"),
                (100, 156, "low"), (127, 129, "low")),
        "low-gf16": ((512, 1536, "low"), (1600, 2496, "low"),
                     (2032, 2064, "low")),
        "padding": ((129, 1, "high"), (129, 100, "high"),
                    (200, 50, "high"), (225, 30, "high")),
    }
    shard_sizes = (64, 256, 1024, 4096, 16384, 65536, 262144,
                   1048576, 4194304, 16777216)
    nonzero_layout = (
        (ORDER_AB, 0, MODE_FORCED_SPECIALIZED),
        (ORDER_AB, 1, MODE_FORCED_GENERIC),
        (ORDER_AB, 2, MODE_AUTOMATIC),
        (ORDER_BA, 0, MODE_FORCED_GENERIC),
        (ORDER_BA, 1, MODE_FORCED_SPECIALIZED),
    )
    profile_name = {
        "high": "legacy_high_v1",
        "low": "low_v1",
    }
    rows: List[Tuple[object, ...]] = []

    def add(
        category: str, k: int, r: int, profile: str, shard_bytes: int,
        losses: int, batch: int, reuse: int, threads: int,
        trial: str, sequence: int, mode: str,
    ) -> None:
        rows.append((
            category, k, r, profile_name[profile], "auto", "auto",
            shard_bytes, losses, batch, reuse, iterations, warmup, threads,
            trial, sequence, mode,
        ))

    for category in sorted(expected_rates):
        for k, r, profile in expected_rates[category]:
            maximum = min(k, r)
            losses = tuple(sorted({
                value for value in (0, 1, 2, 4, 8, r // 4, r // 2, maximum)
                if 0 <= value <= maximum
            }))
            for shard_bytes in shard_sizes:
                for loss_count in losses:
                    layout = ((ORDER_SINGLE, 0, MODE_AUTOMATIC),
                              ) if loss_count == 0 else nonzero_layout
                    for trial, sequence, mode in layout:
                        add(category, k, r, profile, shard_bytes, loss_count,
                            1, 8, 1, trial, sequence, mode)

    checkpoint = (
        ("low", 16, 240, "low"),
        ("balanced", 128, 128, "high"),
        ("high", 240, 16, "high"),
    )
    batch_sizes = ((1, 65536), (8, 65536), (64, 4096), (1024, 256))
    for category, k, r, profile in checkpoint:
        for batch, shard_bytes in batch_sizes:
            for reuse in (1, 8, 64, 1024):
                for trial, sequence, mode in nonzero_layout:
                    add(category + "-reuse", k, r, profile, shard_bytes, 8,
                        batch, reuse, 1, trial, sequence, mode)
        for threads in (1, 2, 4, 8, 16, 32, 64, 128):
            add(category + "-scaling", k, r, profile, 4096, 8, 128, 8,
                threads, ORDER_SINGLE, 0, MODE_AUTOMATIC)
    return Counter(rows)


def _validate_required_dimensions(
    jobs: Sequence[Mapping[str, object]], iterations: int, warmup: int,
) -> None:
    actual = Counter(_dimension_signature(job) for job in jobs)
    expected = _expected_required_dimensions(iterations, warmup)
    if actual == expected:
        return
    missing = list((expected - actual).elements())[:3]
    extra = list((actual - expected).elements())[:3]
    raise MatrixError(
        "required dimension/mode set differs from its independent oracle; "
        f"missing={missing!r}, extra={extra!r}")


def make_spec(
    benchmark: str,
    preset: str,
    workers: int,
    iterations: int,
    warmup: int,
    pinned_cpu: int | None = None,
    perf_stat: str | None = None,
    perf_events: Sequence[str] = DEFAULT_PERF_EVENTS,
    require_perf_counters: bool = False,
) -> dict:
    if workers <= 0 or workers > 128:
        raise MatrixError("workers must be in [1, 128]")
    if require_perf_counters and perf_stat is None:
        raise MatrixError("required performance counters need --perf-stat")
    if perf_stat is not None:
        if not perf_stat:
            raise MatrixError("perf-stat command must not be empty")
        if (not perf_events or len(perf_events) != len(set(perf_events)) or
                not all(isinstance(event, str) and PERF_EVENT_RE.match(event)
                        for event in perf_events)):
            raise MatrixError(
                "performance counter events must be unique simple perf event names")
    jobs: List[Dict[str, object]] = []
    if preset == "smoke":
        cells = (
            ("low", 16, 240, "low", 4096, 1),
            ("high", 240, 16, "high", 4096, 1),
        )
        for category, k, r, profile, shard_bytes, losses in cells:
            for order_trial, order_sequence, mode in _ordered_runs(losses):
                jobs.append(_job(
                    benchmark, category, k, r, profile, shard_bytes, losses,
                    1, 2, 1, mode, order_trial, order_sequence,
                    iterations, warmup))
    elif preset == "checkpoint":
        for category, k, r, profile in _checkpoint_cells():
            for shard_bytes in (4096, 65536, 1048576):
                for losses in (0, 1, min(8, k, r)):
                    for order_trial, order_sequence, mode in _ordered_runs(losses):
                        jobs.append(_job(
                            benchmark, category, k, r, profile, shard_bytes,
                            losses, 1, 8, 1, mode, order_trial,
                            order_sequence, iterations, warmup))
    elif preset == "balanced-crossover":
        for shard_bytes in (256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536):
            for losses in (1, 2, 4, 8, 16, 32, 64, 128):
                for order_trial, order_sequence, mode in _ordered_runs(losses):
                    jobs.append(_job(
                        benchmark, "balanced-crossover", 128, 128, "high",
                        shard_bytes, losses, 1, 16, 1, mode,
                        order_trial, order_sequence, iterations, warmup))
    elif preset == "required":
        for category in sorted(RATE_CASES):
            for k, r, profile in RATE_CASES[category]:
                for shard_bytes in FULL_SHARD_BYTES:
                    for losses in loss_counts(k, r):
                        for order_trial, order_sequence, mode in _ordered_runs(losses):
                            jobs.append(_job(
                                benchmark, category, k, r, profile, shard_bytes,
                                losses, 1, 8, 1, mode, order_trial,
                                order_sequence, iterations, warmup))
        # Reuse/batch and thread scaling are separate from the count/size/loss
        # Cartesian product to avoid impossible multi-terabyte allocations.
        for category, k, r, profile in _checkpoint_cells():
            for batch in (1, 8, 64, 1024):
                # Bound each batch to roughly the same public byte footprint;
                # reuse then varies independently without changing allocation.
                shard_bytes = (65536 if batch <= 8 else
                               (4096 if batch == 64 else 256))
                for reuse in (1, 8, 64, 1024):
                    for order_trial, order_sequence, mode in _ordered_runs(
                            min(8, k, r)):
                        jobs.append(_job(
                            benchmark, category + "-reuse", k, r, profile,
                            shard_bytes, min(8, k, r), batch, reuse, 1,
                            mode, order_trial, order_sequence,
                            iterations, warmup))
            for threads in (1, 2, 4, 8, 16, 32, 64, 128):
                jobs.append(_job(
                    benchmark, category + "-scaling", k, r, profile,
                    4096, min(8, k, r), 128, 8, threads,
                    MODE_AUTOMATIC, ORDER_SINGLE, 0, iterations, warmup))
    else:
        raise MatrixError(f"unknown preset: {preset}")

    if pinned_cpu is not None:
        if pinned_cpu < 0:
            raise MatrixError("pinned CPU must be non-negative")
        for job in jobs:
            # Pin every single-thread comparison cell.  Mixed presets retain
            # normal topology-aware CPU groups for their multi-thread scaling
            # rows, which cannot be represented by one pinned CPU.
            if int(job["cpu_count"]) == 1:
                job.pop("cpu_count")
                job["cpu_set"] = [pinned_cpu]

    ids = [str(job["id"]) for job in jobs]
    if len(ids) != len(set(ids)):
        duplicates = sorted(key for key, count in Counter(ids).items() if count > 1)
        raise MatrixError("duplicate job ids: " + ", ".join(duplicates))
    if preset == "required":
        _validate_required_dimensions(jobs, iterations, warmup)
    defaults = {
        "cpu_policy": "physical-first",
        "timeout_seconds": 1800,
        "memory_mb": 4096,
    }
    if perf_stat is not None:
        defaults["performance_counters"] = {
            "provider": "linux-perf-stat",
            "command": perf_stat,
            "events": list(perf_events),
            "optional": not require_perf_counters,
        }
    spec = {
        "schema": SPEC_SCHEMA,
        "root": str(Path.cwd().resolve()),
        "base_seed": 0x4C454F32,
        "workers": workers,
        "defaults": defaults,
        "metadata": {
            "preset": preset,
            "serial_timing_jobs": workers == 1,
            "isolation_status": "not-established-by-generator",
            "counterbalance": (
                "A=forced-specialized,B=forced-generic; "
                "AB=S/G/automatic then BA=G/S"),
            "counterbalance_temporal_order_requires_workers": 1,
            "benchmark": benchmark,
            "pinned_cpu": pinned_cpu,
            "performance_counters": {
                "requested": perf_stat is not None,
                "provider": "linux-perf-stat" if perf_stat is not None else None,
                "events": list(perf_events) if perf_stat is not None else [],
                "optional": (
                    not require_perf_counters if perf_stat is not None else None),
            },
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
    if not isinstance(parameters.get("force_tiled_decode"), bool):
        raise MatrixError(f"job {job_id} has an invalid force_tiled_decode value")
    if not isinstance(parameters.get("force_materialized_decode"), bool):
        raise MatrixError(
            f"job {job_id} has an invalid force_materialized_decode value")
    if (parameters["force_generic_decode"] and
            parameters["force_specialized_decode"]):
        raise MatrixError(f"job {job_id} forces two decoder policies")
    if (parameters["force_tiled_decode"] and
            parameters["force_materialized_decode"]):
        raise MatrixError(f"job {job_id} forces two workspace kernels")
    if (parameters["force_generic_decode"] and
            (parameters["force_tiled_decode"] or
             parameters["force_materialized_decode"])):
        raise MatrixError(
            f"job {job_id} combines generic and specialized workspace policies")
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
        if (not isinstance(job.get("resume_group"), str) or
                job["resume_group"] != job["cpu_group"]):
            raise MatrixError(
                f"job {job['id']} is missing its atomic resume group")
        expected_cell = job.get("benchmark_cell")
        if (not isinstance(expected_cell, dict) or
                set(expected_cell) != set(EXPECTED_CELL_FIELDS)):
            raise MatrixError(
                f"job {job['id']} expected benchmark cell has the wrong fields")
        _, _, declared_mode = _parse_ordered_job_id(job["id"])
        if declared_mode != _mode_from_cell(expected_cell):
            raise MatrixError(
                f"job {job['id']} counterbalance suffix disagrees with its cell")
        counters = job.get("performance_counters")
        if counters is not None:
            counter_executable = counters.get("executable") if isinstance(
                counters, dict) else None
            if (not isinstance(counters, dict) or
                    set(counters) != {
                        "provider", "events", "optional", "executable",
                        "probe_command"} or
                    counters.get("provider") != "linux-perf-stat" or
                    not isinstance(counters.get("events"), list) or
                    not counters["events"] or
                    len(counters["events"]) != len(set(counters["events"])) or
                    not all(isinstance(event, str) and PERF_EVENT_RE.match(event)
                            for event in counters["events"] or []) or
                    not isinstance(counters.get("optional"), bool) or
                    not isinstance(counter_executable, dict) or
                    not isinstance(counter_executable.get("path"), str) or
                    not isinstance(counter_executable.get("sha256"), str) or
                    len(counter_executable["sha256"]) != 64 or
                    any(character not in "0123456789abcdef"
                        for character in counter_executable["sha256"]) or
                    not isinstance(counter_executable.get("size_bytes"), int) or
                    counter_executable["size_bytes"] < 0 or
                    not probe_command_matches_request(
                        counters, counters.get("probe_command"))):
                raise MatrixError(
                    f"job {job['id']} has an invalid performance counter request")
    return value


def _validate_performance_evidence(
    job: Mapping[str, object], lab_result: Mapping[str, object], job_dir: Path,
) -> dict | None:
    request = job.get("performance_counters")
    evidence = lab_result.get("performance_counters")
    outputs = lab_result.get("outputs")
    if not isinstance(outputs, dict):
        raise MatrixError(f"job {job['id']} output identities are missing")
    if request is None:
        if evidence is not None or "performance_counters" in outputs:
            raise MatrixError(
                f"job {job['id']} has unexpected performance counter evidence")
        return None
    if (not isinstance(request, dict) or not isinstance(evidence, dict) or
            not set(evidence).issubset({
                "provider", "events", "optional", "executable", "probe",
                "status", "raw_output", "measurements", "detail"}) or
            evidence.get("provider") != request.get("provider") or
            evidence.get("events") != request.get("events") or
            evidence.get("optional") != request.get("optional") or
            evidence.get("executable") != request.get("executable") or
            evidence.get("status") not in
            ("available", "partial", "unavailable")):
        raise MatrixError(
            f"job {job['id']} performance counter evidence differs from its request")
    probe = evidence.get("probe")
    if (not isinstance(probe, dict) or
            not {"status", "command", "cpu_set", "exit_code",
                 "duration_seconds"}.issubset(probe) or
            not set(probe).issubset({
                "status", "command", "cpu_set", "exit_code",
                "duration_seconds", "stderr_tail", "detail"}) or
            probe.get("status") not in ("available", "unavailable") or
            probe.get("cpu_set") != job.get("cpu_set") or
            not probe_command_matches_request(request, probe.get("command")) or
            (probe.get("exit_code") is not None and
             (isinstance(probe.get("exit_code"), bool) or
              not isinstance(probe.get("exit_code"), int))) or
            (probe.get("status") == "available" and
             (probe.get("exit_code") != 0 or "detail" in probe)) or
            (probe.get("status") == "unavailable" and
             (not isinstance(probe.get("detail"), str) or
              not probe.get("detail"))) or
            isinstance(probe.get("duration_seconds"), bool) or
            not isinstance(probe.get("duration_seconds"), (int, float)) or
            not math.isfinite(float(probe["duration_seconds"])) or
            float(probe["duration_seconds"]) < 0.0 or
            ("stderr_tail" in probe and
             not isinstance(probe["stderr_tail"], str))):
        raise MatrixError(f"job {job['id']} has invalid counter preflight evidence")
    raw_output = evidence.get("raw_output")
    if raw_output is None:
        if "performance_counters" in outputs:
            raise MatrixError(
                f"job {job['id']} has an unsigned counter output identity")
    elif (raw_output != "perf-stat.txt" or
          not isinstance(outputs.get("performance_counters"), dict)):
        raise MatrixError(
            f"job {job['id']} performance counter output identity is invalid")
    measurements = evidence.get("measurements", [])
    if not isinstance(measurements, list):
        raise MatrixError(f"job {job['id']} counter measurements are invalid")
    status = evidence["status"]
    if raw_output is None and measurements:
        raise MatrixError(
            f"job {job['id']} has counter measurements without raw output")

    if raw_output is None:
        expected_detail = probe.get(
            "detail", "performance counters were unavailable during preflight")
        if status != "unavailable" or evidence.get("detail") != expected_detail:
            raise MatrixError(
                f"job {job['id']} counter status disagrees with its preflight")
        expected_keys = {
            "provider", "events", "optional", "executable", "probe",
            "status", "detail"}
        if set(evidence) != expected_keys:
            raise MatrixError(
                f"job {job['id']} has noncanonical unavailable counter evidence")
        return dict(evidence)

    if probe["status"] != "available":
        raise MatrixError(
            f"job {job['id']} has raw counters from an unavailable preflight")
    try:
        (raw_identity, derived_measurements, derived_status,
         derived_detail) = read_perf_stat_evidence(
             job_dir / "perf-stat.txt", request["events"])
    except OSError as error:
        raise MatrixError(
            f"cannot read retained counters for job {job['id']}: {error}") from error
    if raw_identity != outputs["performance_counters"]:
        raise MatrixError(
            f"job {job['id']} performance counter output identity is invalid")
    if measurements != derived_measurements or status != derived_status:
        raise MatrixError(
            f"job {job['id']} counter JSON differs from retained raw output")
    if ((derived_detail is None and "detail" in evidence) or
            (derived_detail is not None and
             evidence.get("detail") != derived_detail)):
        raise MatrixError(
            f"job {job['id']} counter detail differs from retained raw output")
    expected_keys = {
        "provider", "events", "optional", "executable", "probe", "status",
        "raw_output", "measurements"}
    if derived_detail is not None:
        expected_keys.add("detail")
    if set(evidence) != expected_keys:
        raise MatrixError(
            f"job {job['id']} has noncanonical raw counter evidence")
    return dict(evidence)


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
        run_epoch = lab_result.get("run_epoch")
        if (not isinstance(run_epoch, str) or len(run_epoch) != 64 or
                any(character not in "0123456789abcdef"
                    for character in run_epoch)):
            raise MatrixError(f"job {job_id} has no valid run epoch")
        order_trial, order_sequence, _ = _parse_ordered_job_id(job_id)
        performance_counters = _validate_performance_evidence(
            job, lab_result, job_dir)
        record = {
            "job_id": job_id,
            "outcome": outcome,
            "exit_code": lab_result.get("exit_code"),
            "cpu_set": lab_result.get("cpu_set"),
            "duration_seconds": lab_result.get("duration_seconds"),
            "executable": executable,
            "cpu_group": job["cpu_group"],
            "resume_group": job["resume_group"],
            "run_epoch": run_epoch,
            "order_trial": order_trial,
            "order_sequence": order_sequence,
            "performance_counter_request": job.get("performance_counters"),
            "performance_counters": performance_counters,
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
            _validate_profile_resolution(job_id, expected_cell, benchmark)
            record["benchmark"] = benchmark
            record["build_digest"] = digest(benchmark["build"])
        else:
            record["stderr"] = (job_dir / "stderr.txt").read_text(
                encoding="utf-8", errors="replace")[-4096:]
        records.append(record)

    group_epochs: Dict[str, set[str]] = {}
    for record in records:
        group_epochs.setdefault(record["resume_group"], set()).add(
            record["run_epoch"])
    mixed_groups = sorted(
        group for group, epochs in group_epochs.items() if len(epochs) != 1)
    if mixed_groups:
        raise MatrixError(
            "logical resume groups contain mixed run epochs: " +
            ", ".join(mixed_groups[:3]))

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
            "resume_group": record["resume_group"],
            "run_epoch": record["run_epoch"],
            "order_trial": record["order_trial"],
            "performance_counter_request": record[
                "performance_counter_request"],
            "performance_counter_status": (
                record["performance_counters"]["status"]
                if record["performance_counters"] is not None else
                "not-requested"),
        }
        # AB and BA deliberately share one cpu_group so the lab assigns the
        # same CPU set.  They remain distinct comparison repetitions here.
        key = record["cpu_group"] + "|order-" + record["order_trial"]
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
            "order_trial": pair_parameters[key]["order_trial"],
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
                "order_trial": pair_parameters[key]["order_trial"],
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
        "performance_counter_summary": dict(sorted(Counter(
            (record["performance_counters"]["status"]
             if record["performance_counters"] is not None else
             "not-requested")
            for record in records).items())),
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
    if sorted(smoke_groups.values()) != [5, 5]:
        raise MatrixError("paired jobs do not share stable CPU-assignment groups")
    counterbalanced = (
        (ORDER_AB, 0, MODE_FORCED_SPECIALIZED),
        (ORDER_AB, 1, MODE_FORCED_GENERIC),
        (ORDER_AB, 2, MODE_AUTOMATIC),
        (ORDER_BA, 0, MODE_FORCED_GENERIC),
        (ORDER_BA, 1, MODE_FORCED_SPECIALIZED),
    )
    for cpu_group in smoke_groups:
        members = [job for job in smoke["jobs"] if job["cpu_group"] == cpu_group]
        if tuple(_parse_ordered_job_id(str(job["id"])) for job in members) != counterbalanced:
            raise MatrixError("sorted smoke jobs do not execute adjacent S/G/A then G/S")
        positions = [smoke["jobs"].index(job) for job in members]
        if positions != list(range(min(positions), max(positions) + 1)):
            raise MatrixError("counterbalanced members are not lexicographically adjacent")
        if len({job["benchmark_cell"]["seed"] for job in members}) != 1:
            raise MatrixError("AB/BA members do not share a deterministic seed")
    first = make_spec("/tmp/bench_leopard2", "checkpoint", 1, 3, 1)
    second = make_spec("/tmp/bench_leopard2", "checkpoint", 1, 3, 1)
    if canonical_bytes(first) != canonical_bytes(second):
        raise MatrixError("spec generation is not deterministic")
    if len(first["jobs"]) != 99:
        raise MatrixError("checkpoint job-count invariant failed")
    pinned = make_spec("/tmp/bench_leopard2", "checkpoint", 1, 3, 1, 7)
    if any(job.get("cpu_set") != [7] or "cpu_count" in job for job in pinned["jobs"]):
        raise MatrixError("pinned checkpoint invariant failed")
    instrumented = make_spec(
        "/tmp/bench_leopard2", "smoke", 1, 3, 1, 7,
        "/tmp/perf", ("cycles", "instructions"), False)
    instrumented_request = instrumented["defaults"].get(
        "performance_counters", {})
    if (instrumented_request.get("provider") != "linux-perf-stat" or
            instrumented_request.get("events") != ["cycles", "instructions"] or
            instrumented_request.get("optional") is not True or
            instrumented["metadata"]["performance_counters"]["requested"] is
            not True):
        raise MatrixError("performance counter generation invariant failed")
    try:
        make_spec(
            "/tmp/bench_leopard2", "smoke", 1, 3, 1, None,
            "/tmp/perf", ("cycles", "cycles"), False)
    except MatrixError:
        pass
    else:
        raise MatrixError("duplicate performance events were accepted")
    required = make_spec("/tmp/bench_leopard2", "required", 1, 3, 1)
    if len(required["jobs"]) != 7134:
        raise MatrixError("required preset job-count invariant failed")
    pinned_required = make_spec(
        "/tmp/bench_leopard2", "required", 1, 3, 1, 7)
    for job in pinned_required["jobs"]:
        threads = int(job["benchmark_cell"]["thread_count"])
        if threads == 1:
            if job.get("cpu_set") != [7] or "cpu_count" in job:
                raise MatrixError("single-thread required comparison was not pinned")
        elif job.get("cpu_count") != threads or "cpu_set" in job:
            raise MatrixError("multi-thread required scaling row was incorrectly pinned")
    del pinned_required
    crossover = make_spec("/tmp/bench_leopard2", "balanced-crossover", 1, 3, 1, 7)
    if len(crossover["jobs"]) != 360:
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

    reuse_jobs = [job for job in required["jobs"] if "-reuse.k" in str(job["id"])]
    if len(reuse_jobs) != 240:
        raise MatrixError("reuse/batch job-count invariant failed")
    for category in ("low-reuse", "balanced-reuse", "high-reuse"):
        logical = {
            (job["benchmark_cell"]["batch"], job["benchmark_cell"]["reuse"],
             job["benchmark_cell"]["shard_bytes"])
            for job in reuse_jobs if str(job["id"]).startswith(category + ".")
        }
        expected_logical = {
            (batch, reuse, shard_bytes)
            for batch, shard_bytes in ((1, 65536), (8, 65536),
                                       (64, 4096), (1024, 256))
            for reuse in (1, 8, 64, 1024)
        }
        if logical != expected_logical:
            raise MatrixError(f"{category} does not independently vary reuse and batch")

    gf16_cases = {
        (job["benchmark_cell"]["K"], job["benchmark_cell"]["R"])
        for job in required["jobs"]
        if str(job["id"]).startswith("low-gf16.")
    }
    if gf16_cases != {(512, 1536), (1600, 2496), (2032, 2064)}:
        raise MatrixError("larger GF16 low-rate analogues are incomplete")
    expected_resolutions = {
        (512, 1536): (512, 2048),
        (1600, 2496): (2048, 8192),
        (2032, 2064): (2048, 8192),
    }
    if GF16_LOW_RESOLUTION != expected_resolutions:
        raise MatrixError("GF16 low-rate resolved-parent table changed unexpectedly")
    for (k, r), (expected_padded, expected_parent) in expected_resolutions.items():
        padded = 1 << (k - 1).bit_length()
        parent = 1 << (padded + r - 1).bit_length()
        if (padded != expected_padded or parent != expected_parent or
                parent <= 256):
            raise MatrixError("GF16 analogue parent-size expectation is wrong")
        _validate_profile_resolution(
            f"low-gf16.k{k}.r{r}.fixture", {"K": k, "R": r},
            {"resolved": {
                "profile": "low_v1", "field": "gf16",
                "padded_side": padded, "parent_count": parent}})

    # Every comparable nonzero cell has exactly one automatic observation and
    # two counterbalanced forced repetitions, all contiguous in sorted order.
    required_groups: Dict[str, List[Dict[str, object]]] = {}
    for job in required["jobs"]:
        if "-scaling.k" in str(job["id"]):
            continue
        required_groups.setdefault(str(job["cpu_group"]), []).append(job)
    for members in required_groups.values():
        losses = int(members[0]["benchmark_cell"]["loss_count"])
        parsed = tuple(_parse_ordered_job_id(str(job["id"])) for job in members)
        expected_order = ((ORDER_SINGLE, 0, MODE_AUTOMATIC),
                          ) if losses == 0 else counterbalanced
        if parsed != expected_order:
            raise MatrixError("required cell has incomplete or swapped mode/order coverage")
        positions = [required["jobs"].index(job) for job in members]
        if positions != list(range(min(positions), max(positions) + 1)):
            raise MatrixError("required counterbalance group is not adjacent")
        if len({job["benchmark_cell"]["seed"] for job in members}) != 1:
            raise MatrixError("required AB/BA trial seeds differ")

    def expect_error(label: str, callback) -> None:
        try:
            callback()
        except MatrixError:
            return
        raise MatrixError(f"negative collector test did not reject {label}")

    wrong_dimension = json.loads(json.dumps(required["jobs"]))
    wrong_dimension[0]["benchmark_cell"]["shard_bytes"] = 65
    expect_error(
        "same-cardinality required dimension substitution",
        lambda: _validate_required_dimensions(wrong_dimension, 3, 1))
    expect_error(
        "swapped AB slot",
        lambda: _parse_ordered_job_id(
            "cell.k1.r1.high.b64.l1.batch1.reuse1.t1."
            "order-ab.slot00-forced-generic"))
    expect_error(
        "missing counterbalance suffix",
        lambda: _parse_ordered_job_id("cell.k1.r1.high.forced-generic"))
    expect_error(
        "wrong GF16 resolved parent",
        lambda: _validate_profile_resolution(
            "low-gf16.k512.r1536.fixture", {"K": 512, "R": 1536},
            {"resolved": {
                "profile": "low_v1", "field": "gf16",
                "padded_side": 512, "parent_count": 4096}}))

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
                (MODE_FORCED_SPECIALIZED, 10.0, ORDER_AB, 0, "pair"),
                (MODE_FORCED_GENERIC, 20.0, ORDER_AB, 1, "pair"),
                (MODE_AUTOMATIC, 12.0, ORDER_AB, 2, "pair"),
                (MODE_FORCED_GENERIC, 22.0, ORDER_BA, 0, "pair"),
                (MODE_FORCED_SPECIALIZED, 11.0, ORDER_BA, 1, "pair"),
            ]
            if duplicate:
                modes.append((
                    MODE_FORCED_SPECIALIZED, 11.0, ORDER_AB, 0,
                    "duplicate"))
            for mode, median, trial, sequence, id_prefix in modes:
                forced_generic = mode == MODE_FORCED_GENERIC
                forced_specialized = mode == MODE_FORCED_SPECIALIZED
                job_id = (
                    f"{id_prefix}.order-{trial}.slot{sequence:02d}-{mode}")
                job = {
                    "id": job_id,
                    "cpu_set": [7],
                    "cpu_group": "self-test-pair",
                    "resume_group": "self-test-pair",
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
                    "run_epoch": "a" * 64,
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
                        "force_tiled_decode": False,
                        "force_materialized_decode": False,
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
        if (collected["record_count"] != 5 or collected["comparison_count"] != 2 or
                collected["dispatcher_check_count"] != 1 or
                [comparison["order_trial"]
                 for comparison in collected["comparisons"]] != [ORDER_AB, ORDER_BA] or
                any(comparison["forced_specialized_speedup_vs_forced_generic"] != 2.0
                    for comparison in collected["comparisons"]) or
                collected["dispatcher_checks"][0][
                    "automatic_ratio_to_best_forced"] != 1.2 or
                collected["source_spec"]["metadata"]["preset"] != "self-test"):
            raise MatrixError("collector pairing invariant failed")

        def add_unavailable_counters(manifest, results, benchmarks):
            del benchmarks
            counter_executable = {
                "path": "/tmp/perf", "sha256": "c" * 64,
                "size_bytes": 67890,
            }
            for job in manifest["jobs"]:
                request = {
                    "provider": "linux-perf-stat",
                    "events": ["cycles", "instructions"],
                    "optional": True,
                    "executable": counter_executable,
                }
                request["probe_command"] = perf_probe_command(
                    counter_executable["path"], request["events"],
                    "/usr/bin/python3")
                job["performance_counters"] = request
                unsigned_job = dict(job)
                unsigned_job.pop("job_digest", None)
                job["job_digest"] = digest(unsigned_job)
                result = results[job["id"]]
                result["job_digest"] = job["job_digest"]
                result["performance_counters"] = {
                    "provider": request["provider"],
                    "events": request["events"],
                    "optional": request["optional"],
                    "executable": counter_executable,
                    "status": "unavailable",
                    "probe": {
                        "status": "unavailable",
                        "command": request["probe_command"],
                        "cpu_set": job["cpu_set"],
                        "exit_code": 255,
                        "duration_seconds": 0.001,
                        "detail": "permission denied fixture",
                    },
                    "detail": "permission denied fixture",
                }

        counter_fixture = write_fixture(
            "unavailable-counters", add_unavailable_counters)
        counter_collected = collect(*counter_fixture)
        if (counter_collected["performance_counter_summary"] != {
                "unavailable": 5} or
                any(record["performance_counters"]["status"] != "unavailable"
                    for record in counter_collected["records"])):
            raise MatrixError(
                "collector did not preserve unavailable counter evidence")

        def promote_available_counters(fixture):
            manifest_path, results_path = fixture
            manifest = load_json(manifest_path)
            for job in manifest["jobs"]:
                job_dir = results_path / "jobs" / job["id"]
                raw_path = job_dir / "perf-stat.txt"
                raw_path.write_text(
                    "1000;;cycles;1.000;100.00;\n"
                    "2000;;instructions;1.000;100.00;\n",
                    encoding="utf-8")
                result_path = job_dir / "result.json"
                result = load_json(result_path)
                evidence = result["performance_counters"]
                evidence.update({
                    "status": "available",
                    "raw_output": "perf-stat.txt",
                    "measurements": [
                        {
                            "event": "cycles", "reported_event": "cycles",
                            "raw_value": "1000", "unit": "",
                            "status": "counted", "value": 1000.0,
                            "runtime": "1.000", "running_percentage": 100.0,
                        },
                        {
                            "event": "instructions",
                            "reported_event": "instructions",
                            "raw_value": "2000", "unit": "",
                            "status": "counted", "value": 2000.0,
                            "runtime": "1.000", "running_percentage": 100.0,
                        },
                    ],
                })
                evidence.pop("detail", None)
                evidence["probe"].update(status="available", exit_code=0)
                evidence["probe"].pop("detail", None)
                result["outputs"]["performance_counters"] = file_identity(raw_path)
                result.pop("result_digest", None)
                result["result_digest"] = digest(result)
                write_json(result_path, result)
            return manifest_path, results_path

        available_fixture = promote_available_counters(write_fixture(
            "available-counters", add_unavailable_counters))
        available_collected = collect(*available_fixture)
        if available_collected["performance_counter_summary"] != {
                "available": 5}:
            raise MatrixError(
                "collector did not preserve available counter evidence")

        def mutate_first_counter_result(fixture, mutation):
            manifest_path, results_path = fixture
            manifest = load_json(manifest_path)
            job = manifest["jobs"][0]
            result_path = results_path / "jobs" / job["id"] / "result.json"
            result = load_json(result_path)
            result.pop("result_digest", None)
            mutation(result)
            result["result_digest"] = digest(result)
            write_json(result_path, result)
            return manifest_path, results_path

        missing_raw_fixture = promote_available_counters(write_fixture(
            "available-without-raw", add_unavailable_counters))
        missing_raw_fixture = mutate_first_counter_result(
            missing_raw_fixture,
            lambda result: (
                result["performance_counters"].pop("raw_output"),
                result["outputs"].pop("performance_counters")))
        expect_error(
            "available counter evidence without retained raw output",
            lambda: collect(*missing_raw_fixture))

        wrong_reported_fixture = promote_available_counters(write_fixture(
            "wrong-reported-counter", add_unavailable_counters))
        wrong_reported_fixture = mutate_first_counter_result(
            wrong_reported_fixture,
            lambda result: result["performance_counters"]["measurements"][
                0].update(reported_event="instructions"))
        expect_error(
            "counter evidence relabeled from a different reported event",
            lambda: collect(*wrong_reported_fixture))

        wrong_value_fixture = promote_available_counters(write_fixture(
            "counter-value-disagrees-with-raw", add_unavailable_counters))
        wrong_value_fixture = mutate_first_counter_result(
            wrong_value_fixture,
            lambda result: result["performance_counters"]["measurements"][
                0].update(raw_value="999999999999", value=999999999999.0))
        expect_error(
            "counter JSON value contradicting retained raw output",
            lambda: collect(*wrong_value_fixture))

        bogus_probe_fixture = promote_available_counters(write_fixture(
            "bogus-counter-probe", add_unavailable_counters))
        bogus_probe_fixture = mutate_first_counter_result(
            bogus_probe_fixture,
            lambda result: result["performance_counters"]["probe"].update(
                command=["bogus-probe"], exit_code=99))
        expect_error(
            "forged probe command and successful exit",
            lambda: collect(*bogus_probe_fixture))

        mixed_epoch = write_fixture(
            "mixed-run-epoch",
            lambda manifest, results, benchmarks: results[
                sorted(results)[0]].update(run_epoch="b" * 64))
        expect_error(
            "mixed logical-group run epochs",
            lambda: collect(*mixed_epoch))

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
        "required_dimension_signatures": len(
            _expected_required_dimensions(3, 1)),
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
    generate.add_argument(
        "--perf-stat",
        help=("perf executable to content-address and use for optional Linux "
              "hardware counters"))
    generate.add_argument(
        "--perf-events",
        default=",".join(DEFAULT_PERF_EVENTS),
        help="comma-separated perf events (default: portable CPU/cache/TLB set)")
    generate.add_argument(
        "--require-perf-counters", action="store_true",
        help="record jobs unavailable instead of running when counters cannot be read")
    generate.add_argument("--output", required=True)
    collect_parser = subparsers.add_parser("collect", help="collect benchmark JSON from lab results")
    collect_parser.add_argument("--manifest", required=True)
    collect_parser.add_argument("--results-dir", required=True)
    collect_parser.add_argument("--output", required=True)
    subparsers.add_parser("self-test")
    args = parser.parse_args(argv)
    try:
        if args.command == "generate":
            perf_events = tuple(
                event.strip() for event in args.perf_events.split(",")
                if event.strip())
            spec = make_spec(
                args.benchmark, args.preset, args.workers, args.iterations,
                args.warmup, args.pinned_cpu, args.perf_stat, perf_events,
                args.require_perf_counters)
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
