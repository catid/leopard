#!/usr/bin/env python3
"""Measure the amortized GF8 T=4 AVX2 encoder against control and main.

The runner consumes immutable, already-built executables.  It pins every
benchmark process to one logical CPU, reserves its SMT sibling, validates
candidate/control executable-section identity, checks embedded source
attestation and route selection, and retains three balanced ABBA-style round
contrasts per cell.
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


SCHEMA = "leopard2-gf8-t4-batch-abba/v1"
SUMMARY_SCHEMA = "leopard2-gf8-t4-batch-summary/v1"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
T95 = {
    3: 4.302652729911275,
    5: 2.7764451051977987,
    9: 2.306004135204166,
}
TARGET_CONTROL_FLOOR = 1.05
KEY_MAIN_FLOOR = 1.05
NEIGHBOR_FLOOR = 1.0 / 1.02
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
ORDERS = (
    ("main", "candidate", "control", "control", "candidate", "main"),
    ("control", "main", "candidate", "candidate", "main", "control"),
    ("candidate", "control", "main", "main", "control", "candidate"),
)
CONTROL_ONLY_ORDERS = (
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
    path = Path(__file__).resolve().with_name("run_t8_two_block_abba.py")
    specification = importlib.util.spec_from_file_location(
        "t4_batch_t8_evidence_support", path)
    require(specification is not None and specification.loader is not None,
            f"cannot load evidence support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


T8_SUPPORT = load_support()
MAIN_SUPPORT = T8_SUPPORT.SUPPORT


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def selected_shape(k: int, r: int) -> bool:
    return r in (3, 4) and (
        3 <= k <= 7 or 9 <= k <= 11
    )


def main_supported(cell: Mapping[str, Any]) -> bool:
    return int(cell["R"]) <= int(cell["K"])


def production_selected(k: int, r: int, shard_bytes: int) -> bool:
    return selected_shape(k, r) and \
        32 <= shard_bytes <= 2048 and shard_bytes % 32 == 0


def checkpoint_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    all_shapes = [
        (k, r)
        for k in (3, 4, 5, 6, 7, 9, 10, 11)
        for r in (3, 4)
    ]
    representative_shapes = ((3, 3), (4, 4), (7, 3), (11, 4))

    for k, r in all_shapes:
        for shard_bytes in (64, 512, 1024, 1984, 2048):
            cells.append({
                "id": f"target-k{k}-r{r}-b{shard_bytes}-q64",
                "K": k,
                "R": r,
                "bytes": shard_bytes,
                "batch": 64,
                "role": "target",
            })
    for k, r in representative_shapes:
        for shard_bytes in (64, 1024, 2048):
            for batch in (1, 8):
                cells.append({
                    "id": f"scale-k{k}-r{r}-b{shard_bytes}-q{batch}",
                    "K": k,
                    "R": r,
                    "bytes": shard_bytes,
                    "batch": batch,
                    "role": "target",
                })
    for k, r in representative_shapes:
        for batch in (1, 8, 64):
            cells.append({
                "id": f"neighbor-k{k}-r{r}-b2112-q{batch}",
                "K": k,
                "R": r,
                "bytes": 2112,
                "batch": batch,
                "role": "boundary_neighbor",
            })

    for index, cell in enumerate(cells):
        cell["seed"] = 0x74340000 + index
        cell["reuse"] = max(128, 8192 // int(cell["batch"]))
    return cells


def smoke_cells() -> list[dict[str, Any]]:
    cells = [
        {
            "id": "target-k4-r4-b64-q64",
            "K": 4, "R": 4, "bytes": 64, "batch": 64,
            "role": "target",
        },
        {
            "id": "target-k4-r4-b1024-q8",
            "K": 4, "R": 4, "bytes": 1024, "batch": 8,
            "role": "target",
        },
        {
            "id": "target-k3-r4-b64-q64-no-main",
            "K": 3, "R": 4, "bytes": 64, "batch": 64,
            "role": "target",
        },
        {
            "id": "neighbor-k4-r4-b2112-q64",
            "K": 4, "R": 4, "bytes": 2112, "batch": 64,
            "role": "boundary_neighbor",
        },
    ]
    for index, cell in enumerate(cells):
        cell["seed"] = 0x7434F000 + index
        cell["reuse"] = max(128, 8192 // int(cell["batch"]))
    return cells


def holdout_cells() -> list[dict[str, Any]]:
    cells = [
        {
            "id": "holdout-k11-r4-b1984-q64",
            "K": 11, "R": 4, "bytes": 1984, "batch": 64,
            "role": "target",
        },
        {
            "id": "holdout-k4-r4-b2048-q1",
            "K": 4, "R": 4, "bytes": 2048, "batch": 1,
            "role": "target",
        },
        {
            "id": "holdout-k11-r4-b2048-q1",
            "K": 11, "R": 4, "bytes": 2048, "batch": 1,
            "role": "target",
        },
    ]
    for index, cell in enumerate(cells):
        cell["seed"] = 0x7434E000 + index
        cell["reuse"] = max(128, 8192 // int(cell["batch"]))
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
        "/usr/bin/prlimit", "--as=268435456",
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", "1",
        "--batch", str(cell["batch"]), "--reuse", str(cell["reuse"]),
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


def positive_encode_metric(result: Mapping[str, Any]) -> float:
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "benchmark metrics are absent")
    encode = metrics.get("encode_execution")
    require(isinstance(encode, dict), "encode metric is absent")
    median = encode.get("median_us_per_batch_call")
    require(isinstance(median, (int, float)) and
            not isinstance(median, bool) and
            math.isfinite(float(median)) and float(median) > 0,
            "encode median is not finite and positive")
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
        "K": cell["K"],
        "R": cell["R"],
        "shard_bytes": cell["bytes"],
        "loss_count": 1,
        "batch": cell["batch"],
        "reuse": cell["reuse"],
        "iterations": iterations,
        "warmup": warmup,
        "thread_count": 1,
        "seed": cell["seed"],
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
    require(isinstance(digests, dict) and
            digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and
                len(digests[name]) == 16
                for name in (
                    "original_data", "transmitted_parity",
                    "recovered_originals")),
            "benchmark workload digests are incomplete")

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
        expected_selected = (
            implementation == "candidate" and
            production_selected(
                int(cell["K"]), int(cell["R"]), int(cell["bytes"]))
        )
        require(isinstance(build, dict) and
                build.get("prevalidated_batch_experiment") is True and
                build.get("high_t4_batch_selected") is expected_selected and
                build.get("source_commit") == source_commit and
                build.get("source_tree") == source_tree and
                build.get("source_tracked_dirty") is False,
                "Leopard2 route marker or source identity changed")
    return {
        "encode_us": positive_encode_metric(result),
        "digests": dict(digests),
        "schema": result["schema"],
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
    failure_output: Path,
) -> dict[str, Any]:
    executable = Path(str(identity["path"]))
    require(sha256(executable) == identity["sha256"],
            f"{implementation} binary changed before execution")
    command = benchmark_command(
        implementation, executable, cell, cpu, iterations, warmup)
    started_ns = time.monotonic_ns()
    failure_prefix = failure_output / (
        f"failure-{implementation}-{started_ns}")

    def persist_failure(stdout: bytes, stderr: bytes) -> None:
        T8_SUPPORT.write_bytes_exclusive(
            failure_prefix.with_suffix(".stdout"), stdout)
        T8_SUPPORT.write_bytes_exclusive(
            failure_prefix.with_suffix(".stderr"), stderr)

    try:
        completed = subprocess.run(
            command, env=CHILD_ENVIRONMENT,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=60.0, check=False)
    except subprocess.TimeoutExpired as error:
        persist_failure(error.stdout or b"", error.stderr or b"")
        raise EvidenceError(
            f"{implementation} timed out after 60 seconds") from error
    elapsed_ns = time.monotonic_ns() - started_ns
    if completed.returncode != 0:
        persist_failure(completed.stdout, completed.stderr)
        raise EvidenceError(
            f"{implementation} failed: " +
            completed.stderr.decode(
                "utf-8", errors="replace")[-1000:])
    require(sha256(executable) == identity["sha256"],
            f"{implementation} binary changed after execution")
    try:
        result = json.loads(completed.stdout.decode("utf-8"))
        normalized = validate_result(
            implementation, result, cell, source_commit, source_tree,
            iterations, warmup)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        persist_failure(completed.stdout, completed.stderr)
        raise EvidenceError(
            f"{implementation} output is not one JSON object: {error}") \
            from error
    except Exception:
        persist_failure(completed.stdout, completed.stderr)
        raise
    return {
        "implementation": implementation,
        "elapsed_ns": elapsed_ns,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "normalized": normalized,
        "result": result,
    }


def confidence_interval(round_log_ratios: Sequence[float]) -> dict[str, Any]:
    require(len(round_log_ratios) in T95,
            "round count has no predeclared Student-t threshold")
    center = statistics.mean(round_log_ratios)
    half_width = \
        T95[len(round_log_ratios)] * statistics.stdev(round_log_ratios) / \
        math.sqrt(len(round_log_ratios))
    return {
        "speedup": math.exp(center),
        "ci95": [
            math.exp(center - half_width),
            math.exp(center + half_width),
        ],
        "round_log_ratios": list(round_log_ratios),
    }


def analyze(
    cell: Mapping[str, Any],
    rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    reference = rounds[0]["invocations"][0]["normalized"]["digests"]
    labels = ("control", "main") if main_supported(cell) else ("control",)
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
        require(len(candidate) == 2,
                "round lacks two candidate observations")
        candidate_log = statistics.mean(math.log(value) for value in candidate)
        for label in contrasts:
            baseline = [
                item["normalized"]["encode_us"] for item in invocations
                if item["implementation"] == label
            ]
            require(len(baseline) == 2,
                    f"round lacks two {label} observations")
            contrasts[label].append(
                statistics.mean(math.log(value) for value in baseline) -
                candidate_log)
    result = {
        "cell": dict(cell),
        "digests": reference,
        "control_over_candidate":
            confidence_interval(contrasts["control"]),
    }
    if "main" in contrasts:
        result["main_over_candidate"] = \
            confidence_interval(contrasts["main"])
    return result


def key_main_cell(cell: Mapping[str, Any]) -> bool:
    return cell["K"] == 4 and cell["R"] == 4 and \
        cell["bytes"] in (64, 1024) and cell["role"] == "target"


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
    parser.add_argument(
        "--preset", choices=("smoke", "checkpoint", "holdout"),
        default="checkpoint")
    parser.add_argument("--rounds", type=int, choices=tuple(T95), default=3)
    parser.add_argument("--iterations", type=int, default=9)
    parser.add_argument("--warmup", type=int, default=4)
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    require(options.iterations >= 3 and options.warmup >= 1,
            "insufficient benchmark repetitions")
    require(not options.output.exists(), "output path already exists")
    options.output.mkdir(parents=True)
    cells = (
        smoke_cells() if options.preset == "smoke"
        else holdout_cells() if options.preset == "holdout"
        else checkpoint_cells()
    )
    round_orders = tuple(
        ORDERS[index % len(ORDERS)] for index in range(options.rounds))
    control_only_round_orders = tuple(
        CONTROL_ONLY_ORDERS[index % len(CONTROL_ONLY_ORDERS)]
        for index in range(options.rounds))
    raw: dict[str, Any] = {
        "schema": SCHEMA,
        "created_utc": MAIN_SUPPORT.utc_now(),
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "main_commit": MAIN_COMMIT,
        "preset": options.preset,
        "cpu": options.cpu,
        "reserved_sibling": options.sibling,
        "iterations": options.iterations,
        "warmup": options.warmup,
        "rounds": options.rounds,
        "cells": [],
    }
    lock_descriptor: int | None = None
    try:
        lock_descriptor = acquire_global_lock()
        identities = {
            "candidate": T8_SUPPORT.file_identity(options.candidate),
            "control": T8_SUPPORT.file_identity(options.control),
            "main": T8_SUPPORT.file_identity(options.main),
        }
        raw["identities"] = identities
        require(
            len({identity["sha256"] for identity in identities.values()}) == 3,
            "candidate, control, and main binaries are not distinct")
        executable_sections = {
            name: T8_SUPPORT.executable_sections_identity(
                Path(str(identities[name]["path"])))
            for name in ("candidate", "control")
        }
        raw["executable_sections"] = executable_sections
        require(
            executable_sections["candidate"]["sections"] ==
                executable_sections["control"]["sections"],
            "candidate and control executable instruction sections differ")
        require(set(os.sched_getaffinity(0)) == {options.cpu},
                "runner must be singleton-pinned to the benchmark CPU")
        sibling_text = Path(
            f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
            "thread_siblings_list").read_text(encoding="ascii")
        require(MAIN_SUPPORT.parse_cpu_list(sibling_text) ==
                {options.cpu, options.sibling},
                "requested CPUs are not one SMT pair")
        raw["host"] = {
            "runner_affinity": sorted(os.sched_getaffinity(0)),
            "benchmark_cpu": MAIN_SUPPORT.cpu_policy_identity(options.cpu),
            "reserved_sibling":
                MAIN_SUPPORT.cpu_policy_identity(options.sibling),
        }
        with MAIN_SUPPORT.StableLeaseAnchor(), \
                MAIN_SUPPORT.PairLease(
                    options.cpu, options.sibling) as pair_lease:
            raw["pair_lease"] = pair_lease
            before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
            before_ns = time.monotonic_ns()
            time.sleep(5.0)
            presample = MAIN_SUPPORT.isolation_record(
                options.cpu, options.sibling, pair_lease,
                before_ns, time.monotonic_ns(),
                before_cpu,
                MAIN_SUPPORT.cpu_stat_snapshot(options.cpu),
                before_sibling,
                MAIN_SUPPORT.cpu_stat_snapshot(options.sibling))
            raw["presample"] = presample
            require(
                presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] == 0 and
                presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                "CPU pair was not quiet during the presample")

            for cell_index, cell in enumerate(cells):
                cell_raw = {"cell": dict(cell), "rounds": []}
                orders = round_orders if main_supported(cell) \
                    else control_only_round_orders
                for round_index, order in enumerate(orders):
                    before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
                    before_sibling = \
                        MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
                    before_ns = time.monotonic_ns()
                    invocations = []
                    for label in order:
                        invocations.append(run_one(
                            label, identities[label], cell, options.cpu,
                            options.source_commit, options.source_tree,
                            options.iterations, options.warmup,
                            options.output))
                    isolation = MAIN_SUPPORT.isolation_record(
                        options.cpu, options.sibling, pair_lease,
                        before_ns, time.monotonic_ns(),
                        before_cpu,
                        MAIN_SUPPORT.cpu_stat_snapshot(options.cpu),
                        before_sibling,
                        MAIN_SUPPORT.cpu_stat_snapshot(options.sibling))
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
        target_control_failures = [
            result["cell"]["id"] for result in analyses
            if result["cell"]["role"] == "target" and
            result["control_over_candidate"]["ci95"][0] <
                TARGET_CONTROL_FLOOR
        ]
        key_main_failures = [
            result["cell"]["id"] for result in analyses
            if key_main_cell(result["cell"]) and
            result["main_over_candidate"]["ci95"][0] < KEY_MAIN_FLOOR
        ]
        credible_neighbor_regressions = [
            result["cell"]["id"] for result in analyses
            if result["cell"]["role"] == "boundary_neighbor" and
            result["control_over_candidate"]["ci95"][1] < NEIGHBOR_FLOOR
        ]
        main_slower_cells = [
            result["cell"]["id"] for result in analyses
            if "main_over_candidate" in result and
            result["main_over_candidate"]["speedup"] < 1.0
        ]
        raw["completed_utc"] = MAIN_SUPPORT.utc_now()
        T8_SUPPORT.write_exclusive(options.output / "raw.json", raw)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": (
                "accepted"
                if not target_control_failures and
                not key_main_failures and
                not credible_neighbor_regressions
                else "rejected"
            ),
            "source_commit": options.source_commit,
            "source_tree": options.source_tree,
            "main_commit": MAIN_COMMIT,
            "preset": options.preset,
            "cell_count": len(analyses),
            "process_count": sum(
                len(round_value["invocations"])
                for cell_value in raw["cells"]
                for round_value in cell_value["rounds"]),
            "all_digests_matched": True,
            "all_rounds_zero_sibling_nonidle": True,
            "target_control_failures": target_control_failures,
            "key_main_failures": key_main_failures,
            "credible_neighbor_regressions":
                credible_neighbor_regressions,
            "main_slower_cells": main_slower_cells,
            "cells": analyses,
            "binary_sha256": {
                name: identity["sha256"]
                for name, identity in identities.items()
            },
            "candidate_control_executable_sections_sha256":
                executable_sections["candidate"]["combined_sha256"],
            "raw_sha256": sha256(options.output / "raw.json"),
        }
        T8_SUPPORT.write_exclusive(
            options.output / "summary.json", summary)
        print(json.dumps({
            "status": summary["status"],
            "cells": summary["cell_count"],
            "processes": summary["process_count"],
            "target_control_failures": target_control_failures,
            "key_main_failures": key_main_failures,
            "credible_neighbor_regressions":
                credible_neighbor_regressions,
            "main_slower_cells": main_slower_cells,
        }, sort_keys=True))
        return 0 if summary["status"] == "accepted" else 2
    except Exception as error:
        raw["failed_utc"] = MAIN_SUPPORT.utc_now()
        raw["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
        }
        T8_SUPPORT.write_exclusive(
            options.output / "failure.json", raw)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    finally:
        if lock_descriptor is not None:
            os.close(lock_descriptor)


if __name__ == "__main__":
    raise SystemExit(main())
