#!/usr/bin/env python3
"""Run fail-closed same-source GF16 AVX2 neighboring-regime evidence.

This runner compares the exact parent/candidate pair a725bf4/7921998.  It is
intentionally narrow: one pinned CPU, one reserved idle SMT sibling, explicit
GF16 and AVX2, and three independent paired-log ABBA rounds per cell.  A
contaminated round aborts the campaign and is never summarized as timing.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Mapping


CONTROL_COMMIT = "a725bf46846714c4f1e438e5a6a43c2e29b37613"
CANDIDATE_COMMIT = "7921998a599994d063617faf831624551fd63255"
RAW_SCHEMA = "leopard2-gf16-avx2-neighbor-abba/v1"
SUMMARY_SCHEMA = "leopard2-gf16-avx2-neighbor-summary/v1"
MANIFEST_SCHEMA = "leopard2-gf16-avx2-neighbor-manifest/v1"
ORDER = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
REGRESSION_FLOOR = 1.0 / 1.02
TWO_SIDED_T95_DF2 = 4.302652729911275
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


def cell(identifier: str, profile: str, k: int, r: int, shard_bytes: int,
         reuse: int, role: str) -> dict[str, Any]:
    return {
        "id": identifier,
        "profile": profile,
        "K": k,
        "R": r,
        "shard_bytes": shard_bytes,
        "loss_count": min(k, r),
        "reuse": reuse,
        "iterations": 9,
        "warmup": 2,
        "seed": 0x792100 + len(CELLS_BUILD),
        "role": role,
    }


# A list is used while constructing stable per-cell seeds, then frozen below.
CELLS_BUILD: list[dict[str, Any]] = []


def add(*args: Any) -> None:
    CELLS_BUILD.append(cell(*args))


# High-rate neighboring K/R and shard-size regimes.  The 4096/4098 pair is an
# immediate aligned/ragged *physical* GF16 control around the accepted target.
# Native GF16 stores complete two-byte symbols; the versioned padded-odd wire
# layout maps a 4097-byte payload onto this same 4098-byte transform geometry.
add("high-1000-200-4k", "high", 1000, 200, 4096, 4, "high_neighbor")
add("high-1000-200-64k", "high", 1000, 200, 65536, 1, "high_neighbor")
add("high-4096-512-4k", "high", 4096, 512, 4096, 2, "target_recheck")
add("high-4096-512-4098", "high", 4096, 512, 4098, 2, "ragged_neighbor")

# Low-profile GF16 cells exercise the same unconditional non-fused range
# kernel through both encoding and maximum-loss Algorithm 4 decoding.
add("low-512-1536-4k", "low", 512, 1536, 4096, 4, "low_neighbor")
add("low-1600-2496-4k", "low", 1600, 2496, 4096, 2, "low_neighbor")
add("low-2032-2064-4k", "low", 2032, 2064, 4096, 2, "low_neighbor")

# Fused/split physical boundaries.  Native GF16 requires complete two-byte
# symbols; padded-odd payloads 65 and 129 map to physical sizes 66 and 130.
add("high-tail-64", "high", 1000, 200, 64, 64, "tail_control")
add("high-tail-66", "high", 1000, 200, 66, 64, "tail_split")
add("high-tail-128", "high", 1000, 200, 128, 64, "tail_control")
add("high-tail-130", "high", 1000, 200, 130, 64, "tail_split")
add("low-tail-66", "low", 512, 1536, 66, 32, "tail_split")
add("low-tail-130", "low", 512, 1536, 130, 32, "tail_split")
CELLS = tuple(CELLS_BUILD)


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def canonical(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_exclusive(path: Path, value: object) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(path, flags, 0o644)
    try:
        payload = json.dumps(value, indent=2, sort_keys=True,
                             allow_nan=False).encode("utf-8") + b"\n"
        written = os.write(descriptor, payload)
        require(written == len(payload), f"short write: {path}")
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def load_support() -> Any:
    path = Path(__file__).resolve().parents[1] / "main_compare" / "run_abba.py"
    specification = importlib.util.spec_from_file_location(
        "gf16_neighbor_main_compare_support", path)
    require(specification is not None and specification.loader is not None,
            f"cannot load evidence support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()


def git(root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(root), *arguments], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return completed.stdout.strip()


def source_identity(root: Path, commit: str) -> dict[str, Any]:
    require(git(root, "rev-parse", "HEAD") == commit,
            f"unexpected source commit: {root}")
    require(git(root, "status", "--porcelain=v1", "--untracked-files=no") == "",
            f"tracked source changes present: {root}")
    build = root / "build" / "neighbor-release"
    executable = (build / "bench_leopard2").resolve(strict=True)
    archive = (build / "libleopard.a").resolve(strict=True)
    avx2_object = (build / "CMakeFiles" / "leopard2_backend_avx2.dir" /
                   "Leopard2BackendAVX2.cpp.o").resolve(strict=True)
    link_recipe = (build / "CMakeFiles" / "bench_leopard2.dir" /
                   "link.txt").resolve(strict=True)
    return {
        "commit": commit,
        "root": str(root.resolve()),
        "executable": str(executable),
        "executable_sha256": sha256(executable),
        "archive_sha256": sha256(archive),
        "avx2_object_sha256": sha256(avx2_object),
        "link_recipe_sha256": sha256(link_recipe),
        "cmake_cache_sha256": sha256((build / "CMakeCache.txt").resolve(strict=True)),
    }


def validate_topology(cpu: int, sibling: int) -> dict[str, Any]:
    allowed = set(os.sched_getaffinity(0))
    require(allowed == {cpu},
            f"runner must start singleton-pinned to CPU {cpu}; got {sorted(allowed)}")
    records = []
    for logical in (cpu, sibling):
        root = Path(f"/sys/devices/system/cpu/cpu{logical}/topology")
        siblings = SUPPORT.parse_cpu_list(
            (root / "thread_siblings_list").read_text(encoding="ascii").strip())
        records.append({
            "cpu": logical,
            "thread_siblings": sorted(siblings),
            "core_id": int((root / "core_id").read_text(encoding="ascii")),
            "package_id": int((root / "physical_package_id").read_text(encoding="ascii")),
        })
    require(all(record["thread_siblings"] == sorted((cpu, sibling))
                for record in records), "requested CPUs are not mutual SMT siblings")
    require(records[0]["core_id"] == records[1]["core_id"] and
            records[0]["package_id"] == records[1]["package_id"],
            "requested CPUs do not share one physical core")
    return {"runner_affinity": sorted(allowed), "records": records}


def result_metric(result: Mapping[str, Any], name: str) -> float:
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "benchmark metrics are absent")
    metric = metrics.get(name)
    require(isinstance(metric, dict), f"benchmark metric is absent: {name}")
    value = metric.get("median_us_per_batch_call")
    require(isinstance(value, (int, float)) and not isinstance(value, bool) and
            math.isfinite(float(value)) and float(value) > 0,
            f"benchmark metric is invalid: {name}")
    samples = metric.get("samples_us_per_batch_call")
    require(isinstance(samples, list) and len(samples) == 9 and
            all(isinstance(item, (int, float)) and not isinstance(item, bool) and
                math.isfinite(float(item)) and float(item) > 0 for item in samples),
            f"retained samples are invalid: {name}")
    return float(value)


def validate_result(result: object, cell_value: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(result, dict) and result.get("schema") == "leopard2-benchmark-v3",
            "benchmark did not emit schema v3")
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(isinstance(parameters, dict) and isinstance(resolved, dict) and
            isinstance(correctness, dict) and isinstance(digests, dict),
            "benchmark identity/correctness fields are incomplete")
    expected_profile = "legacy_high_v1" if cell_value["profile"] == "high" else "low_v1"
    require(parameters.get("K") == cell_value["K"] and
            parameters.get("R") == cell_value["R"] and
            parameters.get("shard_bytes") == cell_value["shard_bytes"] and
            parameters.get("loss_count") == cell_value["loss_count"] and
            parameters.get("reuse") == cell_value["reuse"] and
            parameters.get("iterations") == cell_value["iterations"] and
            parameters.get("warmup") == cell_value["warmup"] and
            parameters.get("seed") == cell_value["seed"] and
            parameters.get("batch") == 1 and
            parameters.get("thread_count") == 1 and
            parameters.get("requested_profile") == expected_profile and
            parameters.get("requested_field") == "gf16" and
            parameters.get("requested_backend") == "avx2",
            "benchmark parameters differ from the predeclared cell")
    require(resolved.get("profile") == expected_profile and
            resolved.get("field") == "gf16" and
            resolved.get("backend") == "avx2" and
            resolved.get("thread_count") == 1 and
            resolved.get("selected_decode_path") == "tiled" and
            resolved.get("selected_decode_rule") == "workspace_tiled",
            "benchmark resolved a different profile/field/backend/decode path/thread count")
    require(correctness.get("leopard2_round_trip") is True,
            "benchmark round trip failed")
    require(digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and len(digests[name]) == 16
                for name in ("original_data", "transmitted_parity", "recovered_originals")),
            "benchmark workload digests are incomplete")
    return {
        "encode_us": result_metric(result, "encode_execution"),
        "decode_us": result_metric(result, "decode_execution"),
        "digests": dict(digests),
        "resolved": dict(resolved),
    }


def benchmark_command(executable: str, cell_value: Mapping[str, Any], cpu: int) -> list[str]:
    return [
        "/usr/bin/taskset", "-c", str(cpu), executable,
        "--k", str(cell_value["K"]), "--r", str(cell_value["R"]),
        "--profile", str(cell_value["profile"]), "--field", "gf16",
        "--backend", "avx2", "--bytes", str(cell_value["shard_bytes"]),
        "--loss", str(cell_value["loss_count"]), "--batch", "1",
        "--reuse", str(cell_value["reuse"]),
        "--iterations", str(cell_value["iterations"]),
        "--warmup", str(cell_value["warmup"]), "--threads", "1",
        "--seed", str(cell_value["seed"]), "--skip-legacy",
        "--retain-samples", "--report-decode-path", "--json", "-",
    ]


def run_one(label: str, identity: Mapping[str, Any], cell_value: Mapping[str, Any],
            cpu: int) -> dict[str, Any]:
    command = benchmark_command(str(identity["executable"]), cell_value, cpu)
    before_hash = sha256(Path(str(identity["executable"])))
    require(before_hash == identity["executable_sha256"], "binary changed before launch")
    start = time.monotonic_ns()
    completed = SUPPORT.run_process_bounded(
        command, environment=CHILD_ENVIRONMENT, timeout=600.0,
        max_stdout=16 * 1024 * 1024, max_stderr=2 * 1024 * 1024)
    duration_ns = time.monotonic_ns() - start
    require(completed.returncode == 0,
            f"{label} benchmark failed ({completed.returncode}): " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:])
    require(sha256(Path(str(identity["executable"]))) == before_hash,
            "binary changed after launch")
    try:
        parsed = json.loads(completed.stdout.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"benchmark stdout is not one JSON object: {error}") from error
    normalized = validate_result(parsed, cell_value)
    return {
        "implementation": label,
        "command": command,
        "environment": dict(CHILD_ENVIRONMENT),
        "runner_affinity": sorted(os.sched_getaffinity(0)),
        "duration_ns": duration_ns,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "stderr": completed.stderr.decode("utf-8", errors="replace"),
        "normalized": normalized,
        "result": parsed,
    }


def ratio_summary(round_values: list[float]) -> dict[str, Any]:
    require(len(round_values) == 3 and all(value > 0 for value in round_values),
            "three positive round ratios are required")
    logs = [math.log(value) for value in round_values]
    center = statistics.mean(logs)
    half = TWO_SIDED_T95_DF2 * statistics.stdev(logs) / math.sqrt(len(logs))
    speedup = math.exp(center)
    return {
        "round_control_over_candidate": round_values,
        "speedup_control_over_candidate": speedup,
        "confidence_interval_95": [math.exp(center - half), math.exp(center + half)],
        "candidate_slowdown_fraction": 1.0 / speedup - 1.0,
        "regression_floor_control_over_candidate": REGRESSION_FLOOR,
        "point_estimate_regression_over_2_percent": speedup < REGRESSION_FLOOR,
        "credible_regression_over_2_percent": math.exp(center + half) < REGRESSION_FLOOR,
    }


def summarize_cell(cell_value: Mapping[str, Any], rounds: list[dict[str, Any]]) -> dict[str, Any]:
    require(len(rounds) == 3, "cell did not complete three rounds")
    reference_digests = rounds[0]["invocations"][0]["normalized"]["digests"]
    for round_value in rounds:
        require(round_value["isolation"]["accepted"] is True,
                "contaminated round cannot be summarized")
        for invocation in round_value["invocations"]:
            require(invocation["normalized"]["digests"] == reference_digests,
                    "control/candidate workload digests differ")
    output: dict[str, Any] = {
        "cell": dict(cell_value),
        "digests": reference_digests,
        "encode": {},
        "decode": {},
    }
    for output_name, metric_name in (("encode", "encode_us"), ("decode", "decode_us")):
        ratios = []
        for round_value in rounds:
            control = [item["normalized"][metric_name]
                       for item in round_value["invocations"]
                       if item["implementation"] == "control"]
            candidate = [item["normalized"][metric_name]
                         for item in round_value["invocations"]
                         if item["implementation"] == "candidate"]
            require(len(control) == 2 and len(candidate) == 2,
                    "ABBA round is incomplete")
            ratios.append(math.exp(statistics.mean(map(math.log, control)) -
                                   statistics.mean(map(math.log, candidate))))
        output[output_name] = ratio_summary(ratios)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--control-root", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--sibling", type=int, required=True)
    return parser.parse_args()


def main() -> int:
    options = parse_args()
    require(options.cpu >= 0 and options.sibling >= 0 and options.cpu != options.sibling,
            "CPU and sibling must be distinct non-negative values")
    topology = validate_topology(options.cpu, options.sibling)
    control = source_identity(options.control_root.resolve(strict=True), CONTROL_COMMIT)
    candidate = source_identity(options.candidate_root.resolve(strict=True), CANDIDATE_COMMIT)
    options.output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "created_utc": SUPPORT.utc_now(),
        "cells": list(CELLS),
        "order_by_round": [list(value) for value in ORDER],
        "metrics": ["encode_execution", "decode_execution"],
        "field": "gf16",
        "backend": "avx2",
        "threads": 1,
        "batch": 1,
        "promotion_policy": {
            "regression_floor_control_over_candidate": REGRESSION_FLOOR,
            "regression_definition": "candidate time exceeds control by more than 2 percent",
            "statistical_method": "one mean log contrast per independent ABBA round; two-sided Student-t 95 percent CI with df=2",
            "metrics_judged_separately": True,
            "contamination_policy": "abort campaign; do not summarize timings",
        },
        "source": {"control": control, "candidate": candidate},
        "topology": topology,
    }
    write_exclusive(options.output / "manifest.json", manifest)
    raw: dict[str, Any] = {
        "schema": RAW_SCHEMA,
        "manifest": manifest,
        "host": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python": platform.python_version(),
            "cpu": options.cpu,
            "reserved_sibling": options.sibling,
        },
        "cells": [],
    }
    identities = {"control": control, "candidate": candidate}
    try:
        with SUPPORT.StableLeaseAnchor(), \
                SUPPORT.PairLease(options.cpu, options.sibling) as pair_lease:
            raw["pair_lease"] = pair_lease
            # A quiet five-second presample prevents knowingly starting during
            # an already-active sibling workload.
            pre_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
            pre_sibling = SUPPORT.cpu_stat_snapshot(options.sibling)
            pre_start = time.monotonic_ns()
            time.sleep(5.0)
            pre_end = time.monotonic_ns()
            presample = SUPPORT.isolation_record(
                options.cpu, options.sibling, pair_lease, pre_start, pre_end,
                pre_cpu, SUPPORT.cpu_stat_snapshot(options.cpu), pre_sibling,
                SUPPORT.cpu_stat_snapshot(options.sibling))
            raw["presample"] = presample
            require(presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] == 0 and
                    presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                    "CPU pair was not quiet during the five-second presample")
            for cell_value in CELLS:
                cell_raw = {"cell": dict(cell_value), "rounds": []}
                for round_index, order in enumerate(ORDER):
                    before_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
                    before_sibling = SUPPORT.cpu_stat_snapshot(options.sibling)
                    before_ns = time.monotonic_ns()
                    invocations = [
                        run_one(label, identities[label], cell_value, options.cpu)
                        for label in order
                    ]
                    after_ns = time.monotonic_ns()
                    isolation = SUPPORT.isolation_record(
                        options.cpu, options.sibling, pair_lease,
                        before_ns, after_ns, before_cpu,
                        SUPPORT.cpu_stat_snapshot(options.cpu), before_sibling,
                        SUPPORT.cpu_stat_snapshot(options.sibling))
                    round_value = {
                        "round": round_index,
                        "order": list(order),
                        "invocations": invocations,
                        "isolation": isolation,
                    }
                    cell_raw["rounds"].append(round_value)
                    if not isolation["accepted"]:
                        raw["cells"].append(cell_raw)
                        raise EvidenceError(
                            f"contamination in {cell_value['id']} round {round_index}: "
                            f"sibling nonidle="
                            f"{isolation['delta']['reserved_sibling']['nonidle_jiffies']}")
                raw["cells"].append(cell_raw)
        raw["completed_utc"] = SUPPORT.utc_now()
        write_exclusive(options.output / "raw.json", raw)
        cell_summaries = [
            summarize_cell(item["cell"], item["rounds"])
            for item in raw["cells"]
        ]
        regressions = []
        for item in cell_summaries:
            for metric in ("encode", "decode"):
                if item[metric]["point_estimate_regression_over_2_percent"]:
                    regressions.append({"cell": item["cell"]["id"], "metric": metric,
                                        "result": item[metric]})
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "accepted" if not regressions else "neighbor_regression_detected",
            "source": {"control_commit": CONTROL_COMMIT,
                       "candidate_commit": CANDIDATE_COMMIT},
            "cell_count": len(cell_summaries),
            "process_count": len(cell_summaries) * len(ORDER) * 4,
            "all_digests_matched": True,
            "all_rounds_zero_sibling_nonidle": True,
            "regressions_over_2_percent": regressions,
            "cells": cell_summaries,
            "raw_sha256": sha256(options.output / "raw.json"),
            "manifest_sha256": sha256(options.output / "manifest.json"),
        }
        write_exclusive(options.output / "summary.json", summary)
        print(json.dumps({
            "status": summary["status"],
            "cells": len(cell_summaries),
            "processes": summary["process_count"],
            "regressions": len(regressions),
            "output": str(options.output),
        }, sort_keys=True))
        return 0 if not regressions else 2
    except Exception as error:
        raw["failed_utc"] = SUPPORT.utc_now()
        raw["failure"] = {"type": type(error).__name__, "message": str(error)}
        write_exclusive(options.output / "failure.json", raw)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
