#!/usr/bin/env python3
"""Qualify the extended GF8/AVX2 one-block T=8 binding.

The candidate and control must come from the same clean source tree.  Their
core and benchmark text is expected to be byte-identical; a nonzero volatile
data marker disables only the extended byte-count selector in the control.
Exact Leopard main is run only for public shapes satisfying its R <= K API
contract.  Every process is serial, address-space bounded, and pinned to one
CPU while its SMT sibling is reserved and checked for non-idle work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

import run_t8_two_block_abba as base


SCHEMA = "leopard2-gf8-t8-one-block-extended-abba/v1"
SUMMARY_SCHEMA = "leopard2-gf8-t8-one-block-extended-summary/v1"
TARGET_BYTES = (128, 192, 256, 320, 384, 448, 512)
BYTE_NEIGHBORS = (
    64, 65, 127, 129, 191, 193, 255, 257, 319,
    321, 383, 385, 447, 449, 511, 513, 576, 640,
)
BEYOND_SCHEMA = "leopard2-gf8-t8-one-block-beyond512-abba/v1"
BEYOND_SUMMARY_SCHEMA = \
    "leopard2-gf8-t8-one-block-beyond512-summary/v1"
BEYOND_TARGET_BYTES = (576, 640, 704, 768, 832, 896, 960, 1024)
BEYOND_BYTE_NEIGHBORS = (
    512, 513, 575, 577, 639, 641, 703, 705, 767,
    769, 831, 833, 895, 897, 959, 961, 1023, 1025, 1088,
)
BEYOND_PRODUCTION_MASKS = {
    576: 0xFFFF,
    640: 0xFFFF,
    704: 0x4FFD,
    768: 0x4FDD,
    832: 0x4FDD,
    896: 0x0FD4,
    960: 0x4EFD,
    1024: 0x4FCC,
}
TARGET_CONTROL_FLOOR = 1.05
TARGET_MAIN_FLOOR = 1.0
NEIGHBOR_FLOOR = 1.0 / 1.02


def beyond_production_selected(
    k: int,
    r: int,
    shard_bytes: int,
) -> bool:
    if not (5 <= k <= 8 and 5 <= r <= 8):
        return False
    if shard_bytes == 64 or (
            128 <= shard_bytes <= 512 and shard_bytes % 64 == 0):
        return True
    mask = BEYOND_PRODUCTION_MASKS.get(shard_bytes)
    if mask is None:
        return False
    bit = 4 * (k - 5) + (r - 5)
    return (mask & (1 << bit)) != 0


def executable_sections_identity(executable: Path) -> dict[str, Any]:
    """Hash every ELF section carrying executable instructions."""
    readelf_name = shutil.which("readelf")
    objcopy_name = shutil.which("objcopy")
    base.require(readelf_name is not None and objcopy_name is not None,
                 "readelf and objcopy are required for executable provenance")
    readelf = Path(readelf_name).resolve(strict=True)
    objcopy = Path(objcopy_name).resolve(strict=True)
    section_table = subprocess.run(
        [str(readelf), "-W", "-S", str(executable)],
        env=base.CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=10.0, check=False)
    base.require(section_table.returncode == 0,
                 "readelf failed: " +
                 section_table.stderr.decode(
                     "utf-8", errors="replace")[-1000:])
    section_pattern = re.compile(
        r"^\s*\[\s*\d+\]\s+(\.\S+)\s+\S+\s+"
        r"[0-9A-Fa-f]+\s+[0-9A-Fa-f]+\s+[0-9A-Fa-f]+\s+"
        r"\S+\s+([A-Z]*)\s+")
    section_names: list[str] = []
    for line in section_table.stdout.decode(
            "utf-8", errors="strict").splitlines():
        match = section_pattern.match(line)
        if match and "X" in match.group(2):
            section_names.append(match.group(1))
    base.require(".text" in section_names and section_names,
                 "ELF executable sections are incomplete")

    records: list[dict[str, Any]] = []
    combined = hashlib.sha256()
    with tempfile.TemporaryDirectory(
            prefix="leopard-t8-executable-sections-") as directory:
        root = Path(directory)
        for index, section_name in enumerate(section_names):
            output = root / f"section-{index}.bin"
            copied_elf = root / f"copy-{index}.elf"
            completed = subprocess.run(
                [str(objcopy), "--dump-section",
                 f"{section_name}={output}", str(executable),
                 str(copied_elf)],
                env=base.CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, timeout=10.0, check=False)
            base.require(completed.returncode == 0 and output.is_file() and
                         copied_elf.is_file(),
                         f"objcopy failed for {section_name}: " +
                         completed.stderr.decode(
                             "utf-8", errors="replace")[-1000:])
            payload = output.read_bytes()
            digest = hashlib.sha256(payload).hexdigest()
            records.append({
                "name": section_name,
                "size": len(payload),
                "sha256": digest,
            })
            combined.update(section_name.encode("utf-8"))
            combined.update(b"\0")
            combined.update(payload)
    return {
        "sections": records,
        "combined_sha256": combined.hexdigest(),
        "readelf": base.file_identity(readelf),
        "objcopy": base.file_identity(objcopy),
    }


def target_cells(
    target_bytes: Sequence[int] = TARGET_BYTES,
    seed_base: int = 0x1540000,
) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    index = 0
    for k in range(5, 9):
        for r in range(5, 9):
            role = "target_main" if r <= k else "target_control"
            for shard_bytes in target_bytes:
                cells.append({
                    "id": f"{role}-k{k}-r{r}-b{shard_bytes}",
                    "K": k,
                    "R": r,
                    "bytes": shard_bytes,
                    "role": role,
                    "seed": seed_base + index,
                })
                index += 1
    return cells


def neighbor_cells(
    target_bytes: Sequence[int] = TARGET_BYTES,
    byte_neighbors: Sequence[int] = BYTE_NEIGHBORS,
    seed_base: int = 0x1550000,
) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for k, r in ((5, 5), (8, 8)):
        for shard_bytes in byte_neighbors:
            cells.append({
                "id": f"bytes-k{k}-r{r}-b{shard_bytes}",
                "K": k,
                "R": r,
                "bytes": shard_bytes,
                "role": "neighbor",
            })
    for shard_bytes in target_bytes:
        for k, r in ((4, 4), (9, 5), (8, 4), (8, 9)):
            cells.append({
                "id": f"shape-k{k}-r{r}-b{shard_bytes}",
                "K": k,
                "R": r,
                "bytes": shard_bytes,
                "role": "neighbor",
            })
    for index, cell in enumerate(cells):
        cell["seed"] = seed_base + index
    return cells


def confidence_interval(values: Sequence[float]) -> dict[str, Any]:
    base.require(len(values) == 3,
                 "three independent round contrasts are required")
    center = statistics.mean(values)
    half = base.T95_DF2 * statistics.stdev(values) / math.sqrt(3)
    return {
        "speedup": math.exp(center),
        "ci95": [math.exp(center - half), math.exp(center + half)],
        "round_log_ratios": list(values),
    }


def analyze(
    cell: Mapping[str, Any],
    rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    reference = rounds[0]["invocations"][0]["normalized"]["digests"]
    labels = ["control"]
    if cell["role"] == "target_main":
        labels.append("main")
    contrasts: dict[str, list[float]] = {label: [] for label in labels}
    for round_value in rounds:
        base.require(round_value["isolation"]["accepted"] is True,
                     "contaminated round cannot be analyzed")
        invocations = round_value["invocations"]
        base.require(
            all(item["normalized"]["digests"] == reference
                for item in invocations),
            "candidate/control/main workload digests differ")
        candidate = [
            item["normalized"]["encode_us"] for item in invocations
            if item["implementation"] == "candidate"
        ]
        base.require(len(candidate) == 2,
                     "round lacks two candidate observations")
        candidate_log = statistics.mean(
            math.log(value) for value in candidate)
        for label in labels:
            baseline = [
                item["normalized"]["encode_us"] for item in invocations
                if item["implementation"] == label
            ]
            base.require(len(baseline) == 2,
                         f"round lacks two {label} observations")
            contrasts[label].append(
                statistics.mean(math.log(value) for value in baseline) -
                candidate_log)
    output = {
        "cell": dict(cell),
        "digests": reference,
        "control_over_candidate":
            confidence_interval(contrasts["control"]),
    }
    if "main" in contrasts:
        output["main_over_candidate"] = confidence_interval(
            contrasts["main"])
    return output


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
        "--beyond-512", action="store_true",
        help="qualify the separate 576-through-1024-byte selector")
    parser.add_argument(
        "--holdout", action="store_true",
        help="use an independent seed family for confirmation")
    parser.add_argument(
        "--final-selector", action="store_true",
        help=(
            "validate the discovery/holdout intersection; excluded target "
            "shapes become same-source regression neighbors"
        ))
    return parser.parse_args()


def run_cell(
    cell: Mapping[str, Any],
    identities: Mapping[str, Mapping[str, Any]],
    options: argparse.Namespace,
    pair_lease: Mapping[str, Any],
    output: Path,
) -> dict[str, Any]:
    if cell["role"] == "target_main":
        orders = base.TARGET_ORDER
    else:
        orders = base.NEIGHBOR_ORDER
    cell_raw: dict[str, Any] = {"cell": dict(cell), "rounds": []}
    for round_index, order in enumerate(orders):
        before_cpu = base.SUPPORT.cpu_stat_snapshot(options.cpu)
        before_sibling = base.SUPPORT.cpu_stat_snapshot(options.sibling)
        before_ns = time.monotonic_ns()
        invocations = []
        for slot_index, label in enumerate(order):
            invocation = base.run_one(
                label, identities[label], cell, options.cpu,
                options.source_commit, options.source_tree,
                options.iterations, options.warmup, 0,
                campaign=options.campaign,
                final_selector=options.final_selector,
                failure_output=output)
            invocations.append(invocation)
            base.write_exclusive(
                output /
                f"round-{round_index}-slot-{slot_index}.json",
                invocation)
        isolation = base.SUPPORT.isolation_record(
            options.cpu, options.sibling, pair_lease,
            before_ns, time.monotonic_ns(),
            before_cpu, base.SUPPORT.cpu_stat_snapshot(options.cpu),
            before_sibling,
            base.SUPPORT.cpu_stat_snapshot(options.sibling))
        round_record = {
            "round": round_index,
            "order": list(order),
            "invocations": invocations,
            "isolation": isolation,
        }
        cell_raw["rounds"].append(round_record)
        base.write_exclusive(
            output / f"round-{round_index}.json", round_record)
        base.require(isolation["accepted"],
                     f"contaminated {cell['id']} round {round_index}")
    return cell_raw


def main() -> int:
    options = parse_arguments()
    base.require(options.iterations >= 3 and options.warmup >= 1,
                 "insufficient benchmark repetitions")
    base.require(not options.output.exists(),
                 "output path already exists")
    options.output.mkdir(parents=True)
    if options.beyond_512:
        options.campaign = "one-block-beyond512"
        schema = BEYOND_SCHEMA
        summary_schema = BEYOND_SUMMARY_SCHEMA
        target_bytes = BEYOND_TARGET_BYTES
        byte_neighbors = BEYOND_BYTE_NEIGHBORS
        target_seed = 0x1640000
        neighbor_seed = 0x1650000
    else:
        options.campaign = "one-block-extended"
        schema = SCHEMA
        summary_schema = SUMMARY_SCHEMA
        target_bytes = TARGET_BYTES
        byte_neighbors = BYTE_NEIGHBORS
        target_seed = 0x1540000
        neighbor_seed = 0x1550000
    if options.holdout:
        target_seed += 0x10000
        neighbor_seed += 0x10000
    base.require(not options.final_selector or options.beyond_512,
                 "--final-selector requires --beyond-512")
    cells = target_cells(target_bytes, target_seed) + \
        neighbor_cells(target_bytes, byte_neighbors, neighbor_seed)
    if options.beyond_512:
        for cell in cells:
            selected = beyond_production_selected(
                int(cell["K"]), int(cell["R"]), int(cell["bytes"]))
            cell["candidate_selected"] = selected
            if options.final_selector and \
                    cell["role"].startswith("target") and not selected:
                cell["role"] = "excluded_neighbor"
    raw: dict[str, Any] = {
        "schema": schema,
        "created_utc": base.SUPPORT.utc_now(),
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "main_commit": base.MAIN_COMMIT,
        "cpu": options.cpu,
        "reserved_sibling": options.sibling,
        "iterations": options.iterations,
        "warmup": options.warmup,
        "target_bytes": list(target_bytes),
        "campaign": options.campaign,
        "holdout": options.holdout,
        "final_selector": options.final_selector,
        "batch": 64,
        "reuse": 64,
        "cells": [],
    }
    lock_descriptor: int | None = None
    try:
        lock_descriptor = base.acquire_global_lock()
        identities = {
            "candidate": base.file_identity(options.candidate),
            "control": base.file_identity(options.control),
            "main": base.file_identity(options.main),
        }
        raw["identities"] = identities
        base.require(
            len({identity["sha256"] for identity in identities.values()}) == 3,
            "candidate, control, and main binaries are not distinct")
        executable_sections = {
            name: executable_sections_identity(Path(str(identity["path"])))
            for name, identity in identities.items()
            if name in ("candidate", "control")
        }
        raw["executable_sections"] = executable_sections
        base.require(
            executable_sections["candidate"]["sections"] ==
                executable_sections["control"]["sections"],
            "candidate and control executable instruction sections differ")
        base.require(set(os.sched_getaffinity(0)) == {options.cpu},
                     "runner must be singleton-pinned to the benchmark CPU")
        sibling_text = Path(
            f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
            "thread_siblings_list").read_text(encoding="ascii")
        base.require(base.SUPPORT.parse_cpu_list(sibling_text) ==
                     {options.cpu, options.sibling},
                     "requested CPUs are not one SMT pair")
        raw["host"] = {
            "runner_affinity": sorted(os.sched_getaffinity(0)),
            "benchmark_cpu":
                base.SUPPORT.cpu_policy_identity(options.cpu),
            "reserved_sibling":
                base.SUPPORT.cpu_policy_identity(options.sibling),
        }
        cells_root = options.output / "cells"
        cells_root.mkdir()
        base.write_exclusive(options.output / "campaign.json", raw)
        with base.SUPPORT.StableLeaseAnchor(), \
                base.SUPPORT.PairLease(
                    options.cpu, options.sibling) as pair_lease:
            raw["pair_lease"] = pair_lease
            before_cpu = base.SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = base.SUPPORT.cpu_stat_snapshot(options.sibling)
            before_ns = time.monotonic_ns()
            time.sleep(5.0)
            presample = base.SUPPORT.isolation_record(
                options.cpu, options.sibling, pair_lease,
                before_ns, time.monotonic_ns(),
                before_cpu, base.SUPPORT.cpu_stat_snapshot(options.cpu),
                before_sibling,
                base.SUPPORT.cpu_stat_snapshot(options.sibling))
            raw["presample"] = presample
            base.require(
                presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] == 0 and
                presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                "CPU pair was not quiet during the presample")
            for cell_index, cell in enumerate(cells):
                cell_output = cells_root / str(cell["id"])
                cell_output.mkdir()
                raw["active_cell"] = dict(cell)
                cell_raw = run_cell(
                    cell, identities, options, pair_lease, cell_output)
                raw["cells"].append(cell_raw)
                raw.pop("active_cell", None)
                base.write_exclusive(
                    cell_output / "cell.json", cell_raw)
                print(
                    f"{cell_index + 1}/{len(cells)} {cell['id']}",
                    file=sys.stderr, flush=True)

        analyses = [
            analyze(item["cell"], item["rounds"]) for item in raw["cells"]
        ]
        target_failures: list[str] = []
        neighbor_failures: list[str] = []
        for result in analyses:
            role = result["cell"]["role"]
            if role.startswith("target"):
                failed = (
                    result["control_over_candidate"]["ci95"][0] <
                    TARGET_CONTROL_FLOOR
                )
                if role == "target_main":
                    failed = failed or (
                        result["main_over_candidate"]["ci95"][0] <
                        TARGET_MAIN_FLOOR
                    )
                if failed:
                    target_failures.append(result["cell"]["id"])
            elif result["control_over_candidate"]["ci95"][1] < NEIGHBOR_FLOOR:
                neighbor_failures.append(result["cell"]["id"])

        raw["completed_utc"] = base.SUPPORT.utc_now()
        base.write_exclusive(options.output / "raw.json", raw)
        summary = {
            "schema": summary_schema,
            "status": "accepted" if not target_failures and
                not neighbor_failures else "rejected",
            "source_commit": options.source_commit,
            "source_tree": options.source_tree,
            "main_commit": base.MAIN_COMMIT,
            "target_bytes": list(target_bytes),
            "campaign": options.campaign,
            "holdout": options.holdout,
            "final_selector": options.final_selector,
            "cell_count": len(analyses),
            "target_main_count": sum(
                item["cell"]["role"] == "target_main"
                for item in analyses),
            "target_control_count": sum(
                item["cell"]["role"] == "target_control"
                for item in analyses),
            "neighbor_count": sum(
                not item["cell"]["role"].startswith("target")
                for item in analyses),
            "process_count": sum(
                len(round_value["invocations"])
                for item in raw["cells"]
                for round_value in item["rounds"]),
            "all_digests_matched": True,
            "all_rounds_zero_sibling_nonidle": True,
            "target_failures": target_failures,
            "neighbor_failures": neighbor_failures,
            "cells": analyses,
            "binary_sha256": {
                name: identity["sha256"]
                for name, identity in identities.items()
            },
            "candidate_control_executable_sections_sha256":
                executable_sections["candidate"]["combined_sha256"],
            "raw_sha256": base.sha256(options.output / "raw.json"),
        }
        base.write_exclusive(options.output / "summary.json", summary)
        print(json.dumps({
            "status": summary["status"],
            "cells": summary["cell_count"],
            "processes": summary["process_count"],
            "target_failures": target_failures,
            "neighbor_failures": neighbor_failures,
        }, sort_keys=True))
        return 0 if summary["status"] == "accepted" else 2
    except Exception as error:
        raw["failed_utc"] = base.SUPPORT.utc_now()
        raw["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
        }
        base.write_exclusive(options.output / "failure.json", raw)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    finally:
        if lock_descriptor is not None:
            os.close(lock_descriptor)


if __name__ == "__main__":
    raise SystemExit(main())
