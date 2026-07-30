#!/usr/bin/env python3
"""Measure the two-block T=8 encoder candidate against control and main.

This runner is deliberately narrow.  It consumes immutable, already-built
executables, runs the reusable prevalidated-batch API on one logical CPU, and
rejects a round if its reserved SMT sibling performs non-idle work.  The
candidate and control must report the same clean source commit and tree; exact
Leopard main must report the pinned historical commit.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-gf8-t8-two-block-abba/v1"
SUMMARY_SCHEMA = "leopard2-gf8-t8-two-block-summary/v1"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
T95_DF2 = 4.302652729911275
TARGET_CONTROL_FLOOR = 1.05
TARGET_MAIN_FLOOR = 1.0
NEIGHBOR_FLOOR = 1.0 / 1.02
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
TARGET_ORDER = (
    ("main", "candidate", "control", "control", "candidate", "main"),
    ("control", "main", "candidate", "candidate", "main", "control"),
    ("candidate", "control", "main", "main", "control", "candidate"),
)
NEIGHBOR_ORDER = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
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
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def load_support() -> Any:
    path = Path(__file__).resolve().parents[1] / "main_compare" / "run_abba.py"
    specification = importlib.util.spec_from_file_location(
        "t8_two_block_main_compare_support", path)
    require(specification is not None and specification.loader is not None,
            f"cannot load evidence support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_exclusive(path: Path, value: object) -> None:
    payload = json.dumps(
        value, indent=2, sort_keys=True, allow_nan=False
    ).encode("utf-8") + b"\n"
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        require(os.write(descriptor, payload) == len(payload),
                f"short write: {path}")
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    status = resolved.stat()
    require(status.st_size > 0 and os.access(resolved, os.X_OK),
            f"benchmark is not executable: {resolved}")
    return {
        "path": str(resolved),
        "size": status.st_size,
        "mode": status.st_mode & 0o777,
        "sha256": sha256(resolved),
    }


def target_cells(target_bytes: int = 64) -> list[dict[str, Any]]:
    cells = []
    index = 0
    for k in range(9, 17):
        for r in range(5, 9):
            cells.append({
                "id": f"target-k{k}-r{r}-b{target_bytes}",
                "K": k,
                "R": r,
                "bytes": target_bytes,
                "role": "target",
                "seed": 0x142E000 + index,
            })
            index += 1
    return cells


def neighbor_cells(target_bytes: int = 64) -> list[dict[str, Any]]:
    cells = []
    byte_neighbors = (32, 33, 63, 65) if target_bytes == 64 \
        else (64, 65, 255, 257, 1024)
    for k, r in ((9, 5), (16, 8)):
        for shard_bytes in byte_neighbors:
            cells.append({
                "id": f"bytes-k{k}-r{r}-b{shard_bytes}",
                "K": k,
                "R": r,
                "bytes": shard_bytes,
                "role": "byte_neighbor",
            })
    for k, r in (
        (8, 5), (8, 8), (17, 5), (17, 8),
        (9, 4), (16, 4), (9, 9), (16, 9),
    ):
        cells.append({
            "id": f"shape-k{k}-r{r}-b{target_bytes}",
            "K": k,
            "R": r,
            "bytes": target_bytes,
            "role": "shape_neighbor",
        })
    for index, cell in enumerate(cells):
        cell["seed"] = 0x142F000 + index
    return cells


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    common = [
        "/usr/bin/prlimit", "--as=201326592",
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", "1",
        "--batch", "64", "--reuse", "64",
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]), "--json", "-",
    ]
    if implementation == "main":
        return common
    return common[:-2] + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
        "--json", "-",
    ]


def positive_metric(result: Mapping[str, Any], name: str) -> float:
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "benchmark metrics are absent")
    value = metrics.get(name)
    require(isinstance(value, dict), f"missing benchmark metric: {name}")
    median = value.get("median_us_per_batch_call")
    require(isinstance(median, (int, float)) and
            not isinstance(median, bool) and
            math.isfinite(float(median)) and float(median) > 0,
            f"invalid benchmark metric: {name}")
    return float(median)


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    require(isinstance(result, dict), "benchmark output is not an object")
    expected_parameters = {
        "K": cell["K"], "R": cell["R"], "shard_bytes": cell["bytes"],
        "loss_count": 1, "batch": 64, "reuse": 64,
        "iterations": iterations, "warmup": warmup,
        "thread_count": 1, "seed": cell["seed"],
    }
    parameters = result.get("parameters")
    require(isinstance(parameters, dict) and
            all(parameters.get(name) == value
                for name, value in expected_parameters.items()),
            "benchmark parameters differ from the frozen cell")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(isinstance(resolved, dict) and
            resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and
            resolved.get("thread_count") == 1,
            "benchmark resolved a different profile, field, or thread count")
    require(isinstance(digests, dict) and digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and len(digests[name]) == 16
                for name in (
                    "original_data", "transmitted_parity",
                    "recovered_originals")),
            "benchmark digests are incomplete")
    if implementation == "main":
        require(result.get("schema") == "leopard-main-benchmark-v1" and
                isinstance(correctness, dict) and
                correctness.get("round_trip") is True,
                "exact-main benchmark identity or round trip failed")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("main_source_commit") == MAIN_COMMIT,
                "exact-main source identity changed")
    else:
        require(result.get("schema") == "leopard2-benchmark-v5" and
                resolved.get("backend") == "avx2" and
                isinstance(correctness, dict) and
                correctness.get("leopard2_round_trip") is True,
                "Leopard2 benchmark identity or round trip failed")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("prevalidated_batch_experiment") is True and
                build.get("source_commit") == source_commit and
                build.get("source_tree") == source_tree and
                build.get("source_tracked_dirty") is False,
                "Leopard2 embedded source identity changed")
    return {
        "encode_us": positive_metric(result, "encode_execution"),
        "digests": dict(digests),
    }


def run_one(
    implementation: str,
    identity: Mapping[str, Any],
    cell: Mapping[str, Any],
    cpu: int,
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    executable = Path(str(identity["path"]))
    require(sha256(executable) == identity["sha256"],
            f"{implementation} binary changed before execution")
    command = benchmark_command(
        implementation, executable, cell, cpu, iterations, warmup)
    start = time.monotonic_ns()
    completed = subprocess.run(
        command, env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=60.0, check=False)
    elapsed = time.monotonic_ns() - start
    require(completed.returncode == 0,
            f"{implementation} failed: " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:])
    require(sha256(executable) == identity["sha256"],
            f"{implementation} binary changed after execution")
    try:
        result = json.loads(completed.stdout.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(
            f"{implementation} output is not one JSON object: {error}") from error
    return {
        "implementation": implementation,
        "elapsed_ns": elapsed,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "normalized": validate_result(
            implementation, result, cell, source_commit, source_tree,
            iterations, warmup),
        "result": result,
    }


def confidence_interval(round_log_ratios: Sequence[float]) -> dict[str, Any]:
    require(len(round_log_ratios) == 3,
            "three independent round contrasts are required")
    center = statistics.mean(round_log_ratios)
    half = T95_DF2 * statistics.stdev(round_log_ratios) / math.sqrt(3)
    return {
        "speedup": math.exp(center),
        "ci95": [math.exp(center - half), math.exp(center + half)],
        "round_log_ratios": list(round_log_ratios),
    }


def analyze(
    cell: Mapping[str, Any],
    rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    reference = rounds[0]["invocations"][0]["normalized"]["digests"]
    labels = ("control", "main") if cell["role"] == "target" else ("control",)
    contrasts: dict[str, list[float]] = {label: [] for label in labels}
    for round_value in rounds:
        require(round_value["isolation"]["accepted"] is True,
                "contaminated round cannot be analyzed")
        invocations = round_value["invocations"]
        require(all(item["normalized"]["digests"] == reference
                    for item in invocations),
                "candidate/control/main workload digests differ")
        candidate = [
            item["normalized"]["encode_us"] for item in invocations
            if item["implementation"] == "candidate"
        ]
        require(len(candidate) == 2, "round lacks two candidate observations")
        candidate_log = statistics.mean(math.log(value) for value in candidate)
        for label in labels:
            baseline = [
                item["normalized"]["encode_us"] for item in invocations
                if item["implementation"] == label
            ]
            require(len(baseline) == 2,
                    f"round lacks two {label} observations")
            contrasts[label].append(
                statistics.mean(math.log(value) for value in baseline) -
                candidate_log)
    output = {
        "cell": dict(cell),
        "digests": reference,
        "control_over_candidate": confidence_interval(contrasts["control"]),
    }
    if "main" in contrasts:
        output["main_over_candidate"] = confidence_interval(contrasts["main"])
    return output


def acquire_global_lock() -> int:
    descriptor = os.open(LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        os.close(descriptor)
        raise EvidenceError(f"authoritative lock is busy: {LOCK_PATH}")
    return descriptor


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--iterations", type=int, default=9)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument(
        "--target-bytes", type=int, choices=(64, 256), default=64,
        help="dense target shard size; default preserves the original campaign")
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    require(options.iterations >= 3 and options.warmup >= 1,
            "insufficient benchmark repetitions")
    require(not options.output.exists(), "output path already exists")
    options.output.mkdir(parents=True)
    lock_descriptor = acquire_global_lock()
    identities = {
        "candidate": file_identity(options.candidate),
        "control": file_identity(options.control),
        "main": file_identity(options.main),
    }
    require(len({identity["sha256"] for identity in identities.values()}) == 3,
            "candidate, control, and main binaries are not distinct")
    cells = target_cells(options.target_bytes) + \
        neighbor_cells(options.target_bytes)
    raw: dict[str, Any] = {
        "schema": SCHEMA,
        "created_utc": SUPPORT.utc_now(),
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "main_commit": MAIN_COMMIT,
        "cpu": options.cpu,
        "reserved_sibling": options.sibling,
        "iterations": options.iterations,
        "warmup": options.warmup,
        "target_bytes": options.target_bytes,
        "batch": 64,
        "reuse": 64,
        "identities": identities,
        "cells": [],
    }
    try:
        require(set(os.sched_getaffinity(0)) == {options.cpu},
                "runner must be singleton-pinned to the benchmark CPU")
        sibling_text = Path(
            f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
            "thread_siblings_list").read_text(encoding="ascii")
        require(SUPPORT.parse_cpu_list(sibling_text) ==
                {options.cpu, options.sibling},
                "requested CPUs are not one SMT pair")
        with SUPPORT.StableLeaseAnchor(), \
                SUPPORT.PairLease(options.cpu, options.sibling) as pair_lease:
            raw["pair_lease"] = pair_lease
            before_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = SUPPORT.cpu_stat_snapshot(options.sibling)
            before_ns = time.monotonic_ns()
            time.sleep(5.0)
            presample = SUPPORT.isolation_record(
                options.cpu, options.sibling, pair_lease,
                before_ns, time.monotonic_ns(),
                before_cpu, SUPPORT.cpu_stat_snapshot(options.cpu),
                before_sibling, SUPPORT.cpu_stat_snapshot(options.sibling))
            raw["presample"] = presample
            require(
                presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] == 0 and
                presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                "CPU pair was not quiet during the presample")
            for cell_index, cell in enumerate(cells):
                orders = TARGET_ORDER if cell["role"] == "target" \
                    else NEIGHBOR_ORDER
                cell_raw = {"cell": dict(cell), "rounds": []}
                for round_index, order in enumerate(orders):
                    before_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
                    before_sibling = SUPPORT.cpu_stat_snapshot(options.sibling)
                    before_ns = time.monotonic_ns()
                    invocations = [
                        run_one(
                            label, identities[label], cell, options.cpu,
                            options.source_commit, options.source_tree,
                            options.iterations, options.warmup)
                        for label in order
                    ]
                    isolation = SUPPORT.isolation_record(
                        options.cpu, options.sibling, pair_lease,
                        before_ns, time.monotonic_ns(),
                        before_cpu, SUPPORT.cpu_stat_snapshot(options.cpu),
                        before_sibling,
                        SUPPORT.cpu_stat_snapshot(options.sibling))
                    cell_raw["rounds"].append({
                        "round": round_index,
                        "order": list(order),
                        "invocations": invocations,
                        "isolation": isolation,
                    })
                    require(isolation["accepted"],
                            f"contaminated {cell['id']} round {round_index}")
                raw["cells"].append(cell_raw)
                print(
                    f"{cell_index + 1}/{len(cells)} {cell['id']}",
                    file=sys.stderr, flush=True)
        analyses = [
            analyze(item["cell"], item["rounds"]) for item in raw["cells"]
        ]
        target_failures = []
        neighbor_failures = []
        for result in analyses:
            if result["cell"]["role"] == "target":
                if (result["control_over_candidate"]["ci95"][0] <
                        TARGET_CONTROL_FLOOR or
                        result["main_over_candidate"]["ci95"][0] <
                        TARGET_MAIN_FLOOR):
                    target_failures.append(result["cell"]["id"])
            elif result["control_over_candidate"]["ci95"][1] < NEIGHBOR_FLOOR:
                neighbor_failures.append(result["cell"]["id"])
        raw["completed_utc"] = SUPPORT.utc_now()
        write_exclusive(options.output / "raw.json", raw)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "accepted" if not target_failures and
                not neighbor_failures else "rejected",
            "source_commit": options.source_commit,
            "source_tree": options.source_tree,
            "main_commit": MAIN_COMMIT,
            "target_bytes": options.target_bytes,
            "cell_count": len(analyses),
            "target_count": len(target_cells(options.target_bytes)),
            "neighbor_count": len(neighbor_cells(options.target_bytes)),
            "process_count": sum(
                len(round_value["invocations"])
                for item in raw["cells"] for round_value in item["rounds"]),
            "all_digests_matched": True,
            "all_rounds_zero_sibling_nonidle": True,
            "target_failures": target_failures,
            "neighbor_failures": neighbor_failures,
            "cells": analyses,
            "binary_sha256": {
                name: identity["sha256"]
                for name, identity in identities.items()
            },
            "raw_sha256": sha256(options.output / "raw.json"),
        }
        write_exclusive(options.output / "summary.json", summary)
        print(json.dumps({
            "status": summary["status"],
            "cells": summary["cell_count"],
            "processes": summary["process_count"],
            "target_failures": target_failures,
            "neighbor_failures": neighbor_failures,
        }, sort_keys=True))
        return 0 if summary["status"] == "accepted" else 2
    except Exception as error:
        raw["failed_utc"] = SUPPORT.utc_now()
        raw["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
        }
        write_exclusive(options.output / "failure.json", raw)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    finally:
        os.close(lock_descriptor)


if __name__ == "__main__":
    raise SystemExit(main())
