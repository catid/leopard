#!/usr/bin/env python3
"""Validate the extended GF8 T=4 AVX2 binding thresholds.

This wrapper deliberately leaves the original sub-2-KiB campaign unchanged.
It reuses that runner's frozen-binary, source-attestation, parity-digest,
CPU-isolation, and balanced-ABBA machinery while replacing only the production
selector and cell matrix for the per-(K,R) extended thresholds.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any, Mapping


def load_base() -> Any:
    path = Path(__file__).resolve().with_name("run_t4_batch_abba.py")
    specification = importlib.util.spec_from_file_location(
        "t4_extended_abba_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load T=4 ABBA runner: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.SCHEMA = "leopard2-gf8-t4-extended-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-t4-extended-summary/v1"
BASE_BENCHMARK_COMMAND = BASE.benchmark_command
BASE_REQUIRE = BASE.require
BASE_VALIDATE_RESULT = BASE.validate_result


MAXIMUM_BYTES = {
    (3, 3): 16 * 1024,
    (3, 4): 2 * 1024,
    (4, 3): 12 * 1024,
    (4, 4): 6 * 1024,
    (5, 3): 6 * 1024,
    (5, 4): 4 * 1024,
    (6, 3): 16 * 1024,
    (6, 4): 8 * 1024,
    (7, 3): 16 * 1024,
    (7, 4): 4 * 1024,
    (9, 3): 8 * 1024,
    (9, 4): 6 * 1024,
    (10, 3): 8 * 1024,
    (10, 4): 6 * 1024,
    (11, 3): 6 * 1024,
    (11, 4): 4 * 1024,
}
EXTENDED_GRID = (3072, 4096, 6144, 8192, 12288, 16384)


def production_selected(k: int, r: int, shard_bytes: int) -> bool:
    maximum = MAXIMUM_BYTES.get((k, r), 0)
    return 32 <= shard_bytes <= maximum and shard_bytes % 32 == 0


def main_supported(cell: Mapping[str, Any]) -> bool:
    return int(cell["R"]) <= int(cell["K"]) and \
        int(cell["bytes"]) % 64 == 0


def shared_candidate_control() -> bool:
    try:
        candidate_index = sys.argv.index("--candidate") + 1
        control_index = sys.argv.index("--control") + 1
        return Path(sys.argv[candidate_index]).samefile(
            Path(sys.argv[control_index]))
    except (OSError, ValueError, IndexError):
        return False


def require(condition: bool, message: str) -> None:
    if not condition and \
            message == "candidate, control, and main binaries are not distinct" \
            and shared_candidate_control():
        return
    BASE_REQUIRE(condition, message)


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    command = BASE_BENCHMARK_COMMAND(
        implementation, executable, cell, cpu, iterations, warmup)
    if implementation == "control":
        command[-2:-2] = ["--disable-high-t4-binding"]
    return command


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    normalized = BASE_VALIDATE_RESULT(
        implementation, result, cell, source_commit, source_tree,
        iterations, warmup)
    if implementation != "main":
        if not isinstance(result, dict):
            raise BASE.EvidenceError("benchmark output is not an object")
        build = result.get("build")
        BASE_REQUIRE(
            isinstance(build, dict) and
            build.get("high_t4_batch_diagnostic_disabled") is
                (implementation == "control"),
            "same-executable T=4 diagnostic mode changed")
    return normalized


def make_cell(
    identifier: str,
    k: int,
    r: int,
    shard_bytes: int,
    batch: int,
    role: str,
) -> dict[str, Any]:
    return {
        "id": identifier,
        "K": k,
        "R": r,
        "bytes": shard_bytes,
        "batch": batch,
        "role": role,
    }


def finalize_cells(
    cells: list[dict[str, Any]],
    seed_base: int,
) -> list[dict[str, Any]]:
    identifiers = [str(cell["id"]) for cell in cells]
    if len(identifiers) != len(set(identifiers)):
        raise RuntimeError("extended T=4 campaign contains duplicate cell IDs")
    for index, cell in enumerate(cells):
        cell["seed"] = seed_base + index
        cell["reuse"] = max(128, 8192 // int(cell["batch"]))
    return cells


def checkpoint_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    extended_shapes = [
        shape for shape, maximum in sorted(MAXIMUM_BYTES.items())
        if maximum > 2048
    ]

    # Sweep every discovery grid point that the production table selects.
    for k, r in extended_shapes:
        maximum = MAXIMUM_BYTES[(k, r)]
        for shard_bytes in EXTENDED_GRID:
            if shard_bytes > maximum:
                continue
            cells.append(make_cell(
                f"target-k{k}-r{r}-b{shard_bytes}-q64",
                k, r, shard_bytes, 64, "target"))

    # The shortcut is legal even for one binding item.  Its least-amortized
    # points are explicit neighboring-regime regression gates rather than
    # being inferred from the batch-64 promotion targets.
    for k, r in extended_shapes:
        maximum = MAXIMUM_BYTES[(k, r)]
        for batch in (1, 8):
            cells.append(make_cell(
                f"scale-k{k}-r{r}-b{maximum}-q{batch}",
                k, r, maximum, batch, "boundary_neighbor"))

    # Immediately beyond each threshold the frozen candidate must select the
    # same old path as the diagnostic control, without a credible regression.
    for k, r in sorted(MAXIMUM_BYTES):
        maximum = MAXIMUM_BYTES[(k, r)]
        cells.append(make_cell(
            f"neighbor-k{k}-r{r}-b{maximum + 32}-q64",
            k, r, maximum + 32, 64, "boundary_neighbor"))

    return finalize_cells(cells, 0x74350000)


def smoke_cells() -> list[dict[str, Any]]:
    cells = [
        make_cell(
            "smoke-k3-r3-b16384-q64",
            3, 3, 16384, 64, "target"),
        make_cell(
            "smoke-k4-r4-b6144-q1",
            4, 4, 6144, 1, "boundary_neighbor"),
        make_cell(
            "smoke-k5-r3-b6144-q64",
            5, 3, 6144, 64, "target"),
        make_cell(
            "smoke-k11-r4-b4096-q8",
            11, 4, 4096, 8, "boundary_neighbor"),
        make_cell(
            "smoke-k3-r4-b2080-q64-no-main",
            3, 4, 2080, 64, "boundary_neighbor"),
    ]
    return finalize_cells(cells, 0x7435F000)


def holdout_cells() -> list[dict[str, Any]]:
    cells = [
        make_cell(
            "holdout-k4-r3-b12288-q64",
            4, 3, 12288, 64, "target"),
        make_cell(
            "holdout-k6-r3-b8192-q64",
            6, 3, 8192, 64, "target"),
        make_cell(
            "holdout-k6-r4-b4096-q64",
            6, 4, 4096, 64, "target"),
        make_cell(
            "holdout-k6-r4-b6144-q64",
            6, 4, 6144, 64, "target"),
        make_cell(
            "holdout-k7-r4-b4096-q64",
            7, 4, 4096, 64, "target"),
        make_cell(
            "holdout-k9-r3-b6144-q64",
            9, 3, 6144, 64, "target"),
        make_cell(
            "holdout-k9-r4-b6144-q64",
            9, 4, 6144, 64, "target"),
        make_cell(
            "holdout-k10-r3-b6144-q64",
            10, 3, 6144, 64, "target"),
        make_cell(
            "holdout-k10-r4-b3072-q64",
            10, 4, 3072, 64, "target"),
        make_cell(
            "holdout-k10-r4-b6144-q64",
            10, 4, 6144, 64, "target"),
        make_cell(
            "holdout-k11-r3-b6144-q64",
            11, 3, 6144, 64, "target"),
        make_cell(
            "holdout-k11-r4-b3072-q64",
            11, 4, 3072, 64, "target"),
        make_cell(
            "holdout-k11-r4-b4096-q64",
            11, 4, 4096, 64, "target"),
    ]
    return finalize_cells(cells, 0x7435E000)


def key_main_cell(cell: Mapping[str, Any]) -> bool:
    return production_selected(
        int(cell["K"]), int(cell["R"]), int(cell["bytes"])) and \
        int(cell["R"]) <= int(cell["K"])


BASE.production_selected = production_selected
BASE.main_supported = main_supported
BASE.require = require
BASE.benchmark_command = benchmark_command
BASE.validate_result = validate_result
BASE.checkpoint_cells = checkpoint_cells
BASE.smoke_cells = smoke_cells
BASE.holdout_cells = holdout_cells
BASE.key_main_cell = key_main_cell


if __name__ == "__main__":
    raise SystemExit(BASE.main())
