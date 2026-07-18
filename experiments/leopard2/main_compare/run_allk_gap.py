#!/usr/bin/env python3
"""Broad diagnostic Leopard-main versus Leopard2 all-K gap map.

This is intentionally not promotion evidence.  It saturates all allowed CPUs
with independent single-thread cells to find algorithmic regions worth fixing.
Each cell uses an ABBA process order and retains exact JSON from independently
linked exact-main and Leopard2 executables.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence


MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
CURRENT_COMMIT = "2fce390c9855b6c86b7e20fa86625db500757859"
ORDER = ("main", "leopard2", "leopard2", "main")
CHILD_ENV = {
    "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1", "OMP_PROC_BIND": "TRUE",
    "OMP_PLACES": "cores", "PATH": "/usr/bin:/bin", "TZ": "UTC",
}


@dataclasses.dataclass(frozen=True)
class Cell:
    identifier: str
    family: str
    k: int
    r: int
    shard_bytes: int
    losses: int
    redundancy_band: str
    loss_band: str
    seed: int
    iterations: int
    reuse: int
    warmup: int


def ceil_pow2(value: int) -> int:
    return 1 << (value - 1).bit_length()


def parent_for(k: int, r: int) -> tuple[int, int]:
    padded = ceil_pow2(r)
    return ceil_pow2(k + padded), padded


def gf8_r_values(k: int) -> list[tuple[str, int]]:
    max_padded = 1 << ((256 - k).bit_length() - 1)
    maximum = min(k, max_padded)
    values = (("low-R", 1), ("mid-R", max(1, (maximum + 1) // 2)),
              ("max-GF8-R", maximum))
    result: list[tuple[str, int]] = []
    seen: set[int] = set()
    for name, value in values:
        if value not in seen:
            seen.add(value)
            result.append((name, value))
    return result


GF16_K = (
    129, 130, 191, 192, 193, 200, 224, 225, 239, 240, 241, 248,
    249, 255, 256, 257, 300, 511, 512, 513, 1000, 2048, 4096,
)


def gf16_r_values(k: int) -> list[tuple[str, int]]:
    if k <= 255:
        gap = 256 - k
        first_forced = (1 << (gap.bit_length() - 1)) + 1
        values = (("first-GF16-R", first_forced),
                  ("mid-R", min(k, max(first_forced, (k + 3) // 4))),
                  ("high-R", min(k, 512)))
    else:
        values = (("low-R", 1), ("mid-R", min(k, max(2, k // 8))),
                  ("high-R", min(k, 512)))
    result: list[tuple[str, int]] = []
    seen: set[int] = set()
    for name, value in values:
        parent, _ = parent_for(k, value)
        if value <= k and parent > 256 and value not in seen:
            seen.add(value)
            result.append((name, value))
    return result


def make_cells() -> list[Cell]:
    cells: list[Cell] = []
    seed = 0x38000000
    for k in range(1, 256):
        for redundancy_band, r in gf8_r_values(k):
            parent, _ = parent_for(k, r)
            assert parent <= 256
            for shard_bytes in (4096, 65536):
                for loss_band, losses in (("one-loss", 1), ("max-loss", r)):
                    if loss_band == "max-loss" and losses == 1:
                        continue
                    seed += 1
                    cells.append(Cell(
                        f"gf8-k{k}-r{r}-b{shard_bytes}-l{losses}", "gf8-all-k",
                        k, r, shard_bytes, losses, redundancy_band, loss_band,
                        seed, 5, 16 if shard_bytes == 4096 else 4, 1))
    for k in GF16_K:
        for redundancy_band, r in gf16_r_values(k):
            parent, _ = parent_for(k, r)
            assert 256 < parent <= 65536
            for shard_bytes in (512, 4096):
                for loss_band, losses in (("one-loss", 1), ("max-loss", r)):
                    if loss_band == "max-loss" and losses == 1:
                        continue
                    seed += 1
                    cells.append(Cell(
                        f"gf16-k{k}-r{r}-b{shard_bytes}-l{losses}",
                        "gf16-representative", k, r, shard_bytes, losses,
                        redundancy_band, loss_band, seed, 5,
                        16 if shard_bytes == 512 else 8, 1))
    assert len({cell.identifier for cell in cells}) == len(cells)
    return cells


def direct_locator_cutoff(field: str, parent: int) -> int:
    if field == "gf8":
        if parent <= 8: return parent
        if parent == 16: return 8
        if parent in (32, 128): return 9
        if parent == 64: return 8
        return 7
    if parent <= 32: return parent
    if parent == 64: return 34
    if parent == 128: return 24
    if parent in (256, 512): return 16
    return 14


def classify_paths(cell: Cell, result: Mapping[str, Any]) -> dict[str, Any]:
    resolved = result["resolved"]
    profile = resolved["profile"]
    field = resolved["field"]
    backend = resolved["backend"]
    parent = int(resolved["parent_count"])
    padded = int(resolved["padded_side"])
    assert profile == "legacy_high_v1"
    if padded == 1:
        encode = "direct-xor-single-parity"
    else:
        encode = "specialized-high-transform"
    if padded == 1 and cell.losses == 1:
        decode = "direct-xor-one-loss"
        workspace = "direct-bounded"
        locator = "none"
    elif 2 <= cell.k <= 16 and cell.losses <= 4:
        decode = "direct-small-loss-matrix"
        workspace = "direct-bounded"
        locator = "none"
    elif (field == "gf8" and cell.k == 128 and cell.r == 128 and
          cell.losses == 128 and 256 <= cell.shard_bytes <= 1024 * 1024):
        decode = "generic-full-parent-transform"
        workspace = "materialized-N"
        locator = "active-parent-" + (
            "sparse-direct" if padded <= direct_locator_cutoff(field, parent)
            else "walsh")
    else:
        decode = "specialized-high-transform"
        calibrated_materialized = (
            field == "gf8" and cell.k == 224 and cell.r == 32 and
            1 <= cell.losses <= 8 and cell.shard_bytes <= 65536 and
            backend in ("avx2", "ssse3") and
            cell.shard_bytes >= (24576 if backend == "avx2" else 32768))
        tiled_slots = 2 * padded + cell.losses
        if calibrated_materialized or tiled_slots >= parent:
            workspace = "materialized-N"
        else:
            workspace = "tiled-2T-plus-losses"
        permanent_cached = padded > cell.r and \
            cell.r <= direct_locator_cutoff(field, parent)
        effective = cell.r if permanent_cached else padded
        locator = "active-parent-" + (
            "sparse-direct" if effective <= direct_locator_cutoff(field, parent)
            else "walsh")
    return {
        "profile": profile, "field": field, "backend": backend,
        "parent_count": parent, "padded_side": padded,
        "parent_inflation": parent / float(cell.k + cell.r),
        "encode_path": encode, "decode_path": decode,
        "decode_workspace": workspace, "locator_setup": locator,
    }


def command(role: str, executable: Path, cell: Cell, cpu: int,
            with_current_legacy: bool) -> list[str]:
    common = [
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell.k), "--r", str(cell.r),
        "--bytes", str(cell.shard_bytes), "--loss", str(cell.losses),
        "--batch", "1", "--reuse", str(cell.reuse),
        "--iterations", str(cell.iterations), "--warmup", str(cell.warmup),
        "--threads", "1", "--seed", str(cell.seed),
    ]
    if role == "leopard2":
        common.extend(("--profile", "high", "--field", "auto",
                       "--backend", "auto", "--retain-samples"))
        if not with_current_legacy:
            common.append("--skip-legacy")
    common.extend(("--json", "-"))
    return common


def run_one(command_value: Sequence[str], timeout: float) -> dict[str, Any]:
    started = time.time_ns()
    completed = subprocess.run(
        list(command_value), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env=CHILD_ENV, timeout=timeout)
    finished = time.time_ns()
    record: dict[str, Any] = {
        "command": list(command_value), "returncode": completed.returncode,
        "duration_ns": finished - started,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr": completed.stderr.decode("utf-8", errors="replace"),
    }
    if completed.returncode == 0:
        record["result"] = json.loads(completed.stdout.decode("utf-8"))
    else:
        record["stdout"] = completed.stdout.decode("utf-8", errors="replace")
    return record


def metric(record: Mapping[str, Any], path: Sequence[str]) -> float:
    value: Any = record["result"]
    for key in path:
        value = value[key]
    return float(value)


def geometric(values: Sequence[float]) -> float:
    assert values and all(value > 0 and math.isfinite(value) for value in values)
    return math.exp(statistics.fmean(math.log(value) for value in values))


def gap_tags(cell: Cell, paths: Mapping[str, Any], encode_speedup: float,
             decode_first_speedup: float, decode_reuse_speedup: float) -> list[str]:
    tags: list[str] = []
    if encode_speedup < 1.05:
        if paths["encode_path"] == "direct-xor-single-parity":
            tags.append("encode:R1-xor/API-overhead")
        else:
            tags.append("encode:legacy-high-transform-not-faster")
        if cell.k % paths["padded_side"] != 0:
            tags.append("encode:partial-final-message-block")
    if decode_first_speedup < 1.05 or decode_reuse_speedup < 1.05:
        decode_path = paths["decode_path"]
        if decode_path.startswith("direct-"):
            tags.append("decode:direct-path-overhead-or-kernel")
        elif decode_path.startswith("generic-"):
            tags.append("decode:generic-fallback-crossover")
        else:
            tags.append("decode:specialized-high-crossover")
        if paths["decode_workspace"] == "materialized-N":
            tags.append("decode:materialized-N-workspace")
        if paths["locator_setup"].endswith("walsh"):
            tags.append("decode:active-parent-Walsh-setup")
    if paths["parent_inflation"] >= 1.5:
        tags.append("common:dyadic-parent-inflation>=1.5x")
    if cell.shard_bytes <= 4096:
        tags.append("common:small-shard-fixed-cost")
    return sorted(set(tags))


def analyze_cell(cell: Cell, invocations: Sequence[Mapping[str, Any]], cpu: int) \
        -> dict[str, Any]:
    failures = [record for record in invocations if record["returncode"] != 0]
    result: dict[str, Any] = {
        "schema": "leopard2-all-k-gap-cell/v1", "cell": dataclasses.asdict(cell),
        "cpu": cpu, "order": list(ORDER), "invocations": list(invocations),
        "valid": not failures, "diagnostic_not_promotion_evidence": True,
    }
    if failures:
        result["failures"] = failures
        return result
    main = [record for record in invocations if
            Path(record["command"][3]).name == "leopard_main_benchmark"]
    current = [record for record in invocations if
               Path(record["command"][3]).name == "bench_leopard2"]
    assert len(main) == len(current) == 2
    paths = classify_paths(cell, current[0]["result"])
    assert all(classify_paths(cell, record["result"]) == paths for record in current)
    main_encode = geometric([metric(record, ("metrics", "encode_execution",
        "median_us_per_batch_call")) for record in main])
    current_encode = geometric([metric(record, ("metrics", "encode_execution",
        "median_us_per_batch_call")) for record in current])
    main_decode = geometric([metric(record, ("metrics", "decode_including_setup",
        "median_us_per_batch_call")) for record in main])
    current_plan = geometric([metric(record, ("metrics", "decode_plan_setup",
        "median_us")) for record in current])
    current_decode = geometric([metric(record, ("metrics", "decode_execution",
        "median_us_per_batch_call")) for record in current])
    current_codec = geometric([metric(record, ("metrics", "codec_setup",
        "median_us")) for record in current])
    first = current_plan + current_decode
    amortized = current_decode + current_plan / cell.reuse
    encode_speedup = main_encode / current_encode
    first_speedup = main_decode / first
    reuse_speedup = main_decode / amortized
    paths["gap_tags"] = gap_tags(
        cell, paths, encode_speedup, first_speedup, reuse_speedup)
    result.update({
        "selected": paths,
        "timing_us": {
            "main_encode_execution": main_encode,
            "leopard2_codec_setup": current_codec,
            "leopard2_encode_execution": current_encode,
            "main_decode_including_setup": main_decode,
            "leopard2_decode_plan_setup": current_plan,
            "leopard2_decode_execution": current_decode,
            "leopard2_decode_first_use": first,
            "leopard2_decode_amortized_at_reuse": amortized,
        },
        "speedup_main_over_leopard2": {
            "encode": encode_speedup,
            "decode_first_use": first_speedup,
            "decode_at_reuse": reuse_speedup,
            "decode_execution_only_optimistic": main_decode / current_decode,
        },
        "significantly_beats_main_1_05": {
            "encode": encode_speedup >= 1.05,
            "decode_first_use": first_speedup >= 1.05,
            "decode_at_reuse": reuse_speedup >= 1.05,
        },
        "memory": current[0]["result"]["memory"],
    })
    if all(record["result"]["legacy"]["available"] for record in current):
        current_legacy_encode = geometric([metric(
            record, ("legacy", "encode_execution", "median_us_per_batch_call"))
            for record in current])
        current_legacy_decode = geometric([metric(
            record, ("legacy", "decode_including_setup",
                     "median_us_per_batch_call"))
            for record in current])
        result["timing_us"].update({
            "current_tree_legacy_encode_execution": current_legacy_encode,
            "current_tree_legacy_decode_including_setup": current_legacy_decode,
        })
        result["attribution_speedup"] = {
            "exact_main_over_current_tree_legacy": {
                "encode": main_encode / current_legacy_encode,
                "decode_including_setup": main_decode / current_legacy_decode,
            },
            "current_tree_legacy_over_leopard2": {
                "encode": current_legacy_encode / current_encode,
                "decode_first_use": current_legacy_decode / first,
                "decode_at_reuse": current_legacy_decode / amortized,
                "decode_execution_only_optimistic":
                    current_legacy_decode / current_decode,
            },
        }
    return result


def run_cell(cell: Cell, index: int, cpus: Sequence[int], main: Path,
             current: Path, output: Path, timeout: float,
             with_current_legacy: bool) -> dict[str, Any]:
    path = output / "cells" / (cell.identifier + ".json")
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    cpu = cpus[index % len(cpus)]
    invocations = []
    for role in ORDER:
        executable = main if role == "main" else current
        try:
            invocations.append(run_one(command(
                role, executable, cell, cpu, with_current_legacy), timeout))
        except subprocess.TimeoutExpired as error:
            invocations.append({
                "command": list(error.cmd), "returncode": -999,
                "duration_ns": int(timeout * 1e9), "stderr": "timeout",
            })
            break
    result = analyze_cell(cell, invocations, cpu)
    temporary = path.with_suffix(".tmp")
    temporary.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_text(json.dumps(result, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)
    return result


def summarize(results: Sequence[Mapping[str, Any]], metadata: Mapping[str, Any]) \
        -> dict[str, Any]:
    valid = [result for result in results if result.get("valid") is True]
    failed = [result for result in results if result.get("valid") is not True]
    metrics = ("encode", "decode_first_use", "decode_at_reuse")
    summary: dict[str, Any] = {
        "schema": "leopard2-all-k-gap-summary/v1", "metadata": dict(metadata),
        "cell_count": len(results), "valid_cell_count": len(valid),
        "failed_cell_count": len(failed),
        "diagnostic_not_promotion_evidence": True,
        "threshold": "Leopard2 significant when main_time/Leopard2_time >= 1.05",
        "metrics": {}, "attribution_metrics": {}, "gap_tags": {},
        "worst_cells": {},
    }
    for name in metrics:
        values = [result["speedup_main_over_leopard2"][name] for result in valid]
        slow = [value for value in values if value < 1.05]
        summary["metrics"][name] = {
            "count": len(values), "gap_count": len(slow),
            "gap_fraction": len(slow) / len(values) if values else None,
            "median_speedup": statistics.median(values) if values else None,
            "p10_speedup": sorted(values)[max(0, len(values) // 10 - 1)]
            if values else None,
        }
        ranked = sorted(valid, key=lambda result:
                        result["speedup_main_over_leopard2"][name])[:30]
        summary["worst_cells"][name] = [{
            "cell": result["cell"], "selected": result["selected"],
            "speedup": result["speedup_main_over_leopard2"][name],
            "timing_us": result["timing_us"],
        } for result in ranked]
    attribution_paths = {
        "exact_main_over_current_tree_legacy_encode":
            ("exact_main_over_current_tree_legacy", "encode"),
        "exact_main_over_current_tree_legacy_decode":
            ("exact_main_over_current_tree_legacy", "decode_including_setup"),
        "current_tree_legacy_over_leopard2_encode":
            ("current_tree_legacy_over_leopard2", "encode"),
        "current_tree_legacy_over_leopard2_decode_first_use":
            ("current_tree_legacy_over_leopard2", "decode_first_use"),
        "current_tree_legacy_over_leopard2_decode_at_reuse":
            ("current_tree_legacy_over_leopard2", "decode_at_reuse"),
    }
    if valid and all("attribution_speedup" in result for result in valid):
        for name, (leg, metric_name) in attribution_paths.items():
            values = [result["attribution_speedup"][leg][metric_name]
                      for result in valid]
            ordered = sorted(values)
            summary["attribution_metrics"][name] = {
                "count": len(values),
                "median_ratio": statistics.median(values),
                "p10_ratio": ordered[max(0, len(ordered) // 10 - 1)],
                "p90_ratio": ordered[min(len(ordered) - 1,
                                         9 * len(ordered) // 10)],
                "count_below_1_0": sum(value < 1.0 for value in values),
                "count_below_1_05": sum(value < 1.05 for value in values),
            }
    tags: dict[str, list[Mapping[str, Any]]] = {}
    for result in valid:
        for tag in result["selected"]["gap_tags"]:
            tags.setdefault(tag, []).append(result)
    summary["gap_tags"] = {
        tag: {
            "cell_count": len(items),
            "median_encode_speedup": statistics.median(
                item["speedup_main_over_leopard2"]["encode"] for item in items),
            "median_decode_first_speedup": statistics.median(
                item["speedup_main_over_leopard2"]["decode_first_use"]
                for item in items),
            "median_decode_reuse_speedup": statistics.median(
                item["speedup_main_over_leopard2"]["decode_at_reuse"]
                for item in items),
        } for tag, items in sorted(tags.items(), key=lambda pair: -len(pair[1]))
    }
    summary["failed_cells"] = [result["cell"] for result in failed]
    return summary


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main", type=Path, required=True)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=min(30, os.cpu_count() or 1))
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--with-current-legacy", action="store_true",
                        help="also time the current tree's retained legacy API")
    options = parser.parse_args(arguments)
    main_exe = options.main.resolve(strict=True)
    current_exe = options.current.resolve(strict=True)
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cpus = sorted(os.sched_getaffinity(0))
    workers = min(options.workers, len(cpus), 30)
    assert workers > 0
    cpus = cpus[:workers]
    cells = make_cells()
    metadata = {
        "main_commit": MAIN_COMMIT, "current_commit": CURRENT_COMMIT,
        "main_executable": str(main_exe), "current_executable": str(current_exe),
        "main_sha256": hashlib.sha256(main_exe.read_bytes()).hexdigest(),
        "current_sha256": hashlib.sha256(current_exe.read_bytes()).hexdigest(),
        "allowed_cpus": sorted(os.sched_getaffinity(0)), "used_cpus": cpus,
        "workers": workers, "order": list(ORDER),
        "with_current_legacy": options.with_current_legacy,
        "matrix": {"gf8_K": [1, 255], "gf8_shard_bytes": [4096, 65536],
                   "gf16_K": list(GF16_K), "gf16_shard_bytes": [512, 4096]},
        "measurement_note": "all CPUs saturated; diagnostic crossover map, not isolated promotion evidence",
    }
    (output / "manifest.json").write_text(
        json.dumps({"metadata": metadata,
                    "cells": [dataclasses.asdict(cell) for cell in cells]},
                   sort_keys=True) + "\n", encoding="utf-8")
    results: list[Mapping[str, Any]] = [None] * len(cells)  # type: ignore[list-item]
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        future_map = {
            executor.submit(run_cell, cell, index, cpus, main_exe, current_exe,
                            output, options.timeout,
                            options.with_current_legacy): index
            for index, cell in enumerate(cells)
        }
        for future in concurrent.futures.as_completed(future_map):
            index = future_map[future]
            results[index] = future.result()
            completed += 1
            if completed % 50 == 0 or completed == len(cells):
                print(f"{completed}/{len(cells)} cells", flush=True)
    with (output / "cells.jsonl").open("w", encoding="utf-8") as stream:
        for result in results:
            stream.write(json.dumps(result, sort_keys=True) + "\n")
    gaps = [result for result in results if result.get("valid") is True and
            not all(result["significantly_beats_main_1_05"].values())]
    with (output / "gap_cells.jsonl").open("w", encoding="utf-8") as stream:
        for result in gaps:
            stream.write(json.dumps(result, sort_keys=True) + "\n")
    summary = summarize(results, metadata)
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output / "summary.json")
    return 0 if not summary["failed_cell_count"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
