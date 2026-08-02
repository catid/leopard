#!/usr/bin/env python3
"""Static checks for the final T32/B256 immutable-binary campaign."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path
from typing import Any


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_runner() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_t32_b256_generated_abba.py")
    specification = importlib.util.spec_from_file_location(
        "leopard2_t32_b256_generated_runner_test_target", path)
    require(specification is not None and specification.loader is not None,
            "cannot load T32/B256 runner")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()
BASE = RUNNER.BASE


def synthetic_round(order: tuple[str, ...]) -> dict[str, Any]:
    digests = {
        "original_data": "1" * 16,
        "transmitted_parity": "2" * 16,
        "recovered_originals": "1" * 16,
    }
    times = {"candidate": 1.0, "control": 1.25, "main": 1.5}
    return {
        "order": list(order),
        "invocations": [{
            "implementation": label,
            "normalized": {
                "encode_us": times[label],
                "digests": dict(digests),
            },
        } for label in order],
        "isolation": {"accepted": True},
    }


def benchmark_result(
    cell: dict[str, Any],
    source_commit: str,
    source_tree: str,
) -> dict[str, Any]:
    return {
        "schema": BASE.CANDIDATE_SCHEMA,
        "parameters": {
            "K": cell["K"],
            "R": cell["R"],
            "shard_bytes": cell["bytes"],
            "loss_count": cell["loss"],
            "batch": cell["batch"],
            "reuse": cell["reuse"],
            "iterations": 15,
            "warmup": 64,
            "thread_count": 1,
            "seed": cell["seed"],
        },
        "resolved": {
            "profile": "legacy_high_v1",
            "field": "gf8",
            "backend": "avx2",
            "thread_count": 1,
        },
        "correctness": {"leopard2_round_trip": True},
        "workload_digests": {
            "algorithm": "fnv1a64",
            "original_data": "1" * 16,
            "transmitted_parity": "2" * 16,
            "recovered_originals": "1" * 16,
        },
        "build": {
            "source_commit": source_commit,
            "source_tree": source_tree,
            "source_tracked_dirty": False,
        },
        "metrics": {
            "encode_execution": {"median_us_per_batch_call": 1.0},
        },
    }


def main_benchmark_result(cell: dict[str, Any]) -> dict[str, Any]:
    physical_bytes = (cell["bytes"] + 63) & ~63
    padded = physical_bytes != cell["bytes"]
    return {
        "schema": "leopard-main-benchmark-v2" if padded else
            "leopard-main-benchmark-v1",
        "parameters": {
            "K": cell["K"],
            "R": cell["R"],
            "shard_bytes": physical_bytes,
            "logical_shard_bytes": cell["bytes"],
            "loss_count": cell["loss"],
            "batch": cell["batch"],
            "reuse": cell["reuse"],
            "iterations": 15,
            "warmup": 64,
            "thread_count": 1,
            "seed": cell["seed"],
        },
        "resolved": {
            "profile": "legacy_high_v1",
            "field": "gf8",
            "thread_count": 1,
            "padded_application_bytes": padded,
            "padding_policy": "zero suffix per shard",
        },
        "correctness": {
            "round_trip": True,
            "logical_prefix_fingerprinted": True,
        },
        "workload_digests": {
            "algorithm": "fnv1a64",
            "original_data": "1" * 16,
            "transmitted_parity": "2" * 16,
            "recovered_originals": "1" * 16,
        },
        "build": {"main_source_commit": BASE.MAIN_COMMIT},
        "metrics": {
            "encode_execution": {"median_us_per_batch_call": 1.0},
        },
    }


def main() -> int:
    cells = RUNNER.cells()
    targets = [cell for cell in cells if cell["role"] == "target"]
    neighbors = [cell for cell in cells if cell["role"] == "neighbor"]
    require(len(cells) == 19 and len(targets) == 6 and len(neighbors) == 13,
            "final T32/B256 matrix cardinality changed")
    require(
        {(cell["bytes"], cell["loss"], cell["batch"])
         for cell in targets} == {
            (256, 1, 1), (256, 2, 1), (256, 4, 1),
            (256, 8, 1), (256, 16, 1), (256, 32, 1),
        },
        "target byte/loss/batch coverage changed")
    require(any(cell["bytes"] == 256 and cell["batch"] == 8
                for cell in neighbors),
            "inert prevalidated-batch control changed")
    require({cell["bytes"] for cell in neighbors} >= {64, 255, 256, 257, 1024},
            "byte-neighbor coverage changed")
    require({(cell["K"], cell["R"]) for cell in neighbors} >= {
                (16, 8), (16, 16), (31, 31), (31, 32),
                (32, 31), (32, 33), (33, 32), (64, 64)},
            "shape-neighbor coverage changed")
    require({cell["id"] for cell in cells if not cell["compare_main"]} == {
                "shape-neighbor-k31-r32-b256-l1-q1",
                "shape-neighbor-k32-r33-b256-l1-q1",
            } and
            all(cell["compare_main"] is True
                for cell in cells if cell["R"] <= cell["K"]),
            "exact-main compatibility classification changed")
    require(len(BASE.TARGET_ORDER) == 9 and
            len(BASE.select_round_orders(BASE.TARGET_ORDER, 9)) == 9,
            "nine-round mirrored schedule changed")

    loss32 = next(cell for cell in targets if cell["loss"] == 32)
    command = BASE.benchmark_command(
        "candidate", Path("/frozen/candidate"), loss32, 13, 15, 64)
    require(command[command.index("--loss") + 1] == "32" and
            command[command.index("--bytes") + 1] == "256" and
            command[command.index("--backend") + 1] == "avx2",
            "loss-aware candidate command changed")

    source_commit = "a" * 40
    source_tree = "b" * 40
    BASE.validate_result(
        "candidate", benchmark_result(loss32, source_commit, source_tree),
        loss32, source_commit, source_tree, 15, 64)

    for logical_bytes, physical_bytes in ((255, 256), (257, 320)):
        byte_neighbor = next(
            cell for cell in neighbors if cell["bytes"] == logical_bytes)
        main_command = BASE.benchmark_command(
            "main", Path("/frozen/main"), byte_neighbor, 13, 15, 64)
        candidate_command = BASE.benchmark_command(
            "candidate", Path("/frozen/candidate"), byte_neighbor,
            13, 15, 64)
        require(
            main_command[main_command.index("--bytes") + 1] ==
                str(physical_bytes) and
            main_command[main_command.index("--logical-bytes") + 1] ==
                str(logical_bytes) and
            candidate_command[candidate_command.index("--bytes") + 1] ==
                str(logical_bytes) and
            "--logical-bytes" not in candidate_command,
            "exact-main logical-byte adapter changed")
        BASE.validate_result(
            "main", main_benchmark_result(byte_neighbor), byte_neighbor,
            source_commit, source_tree, 15, 64)

    neighbor = neighbors[0]
    rounds = [synthetic_round(order) for order in BASE.TARGET_ORDER]
    analysis = BASE.analyze(neighbor, rounds)
    require(math.isclose(
                analysis["control_over_candidate"]["speedup"], 1.25) and
            math.isclose(
                analysis["main_over_candidate"]["speedup"], 1.5),
            "exact-main neighbor contrast changed")
    require(BASE.ALLOW_MULTIPLE_TARGETS is True and
            BASE.TARGET_CONTROL_FLOOR == 1.05 and
            BASE.TARGET_MAIN_FLOOR == 1.05 and
            math.isclose(BASE.NEIGHBOR_FLOOR, 1.0 / 1.02),
            "promotion thresholds changed")
    print("T32/B256 final ABBA runner checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
