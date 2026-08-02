#!/usr/bin/env python3
"""Measure the integrated GF8 T=4 B512 and K8 terminal promotion.

This configured specialization of the hardened T=4 ABBA runner covers all
26 previously promoted cells as regression gates, the eight K=4..7 B=512
cells, the six K=8 B=64/128/256 cells, and two unselected K=8 B=512 controls.
It compares an ordinary candidate with a same-source binary whose complete
T=4 terminal family is disabled by one retained data word, plus the exact
Leopard main codec.  The shared runner proves executable-section identity and
full-file equality after masking only that selector and the GNU build ID.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any


BENCHMARK_CPU = 14
RESERVED_SIBLING = 30


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_t4_packed_terminal_family_abba.py")
    specification = importlib.util.spec_from_file_location(
        "t4_b512_k8_hardened_abba_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load hardened ABBA base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-t4-b512-k8-terminal-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-t4-b512-k8-terminal-summary/v1"
BASE.BENCHMARK_CPU = BENCHMARK_CPU
BASE.RESERVED_SIBLING = RESERVED_SIBLING
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.git_capture.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.git_capture._build_provenance.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.link_common.__file__).resolve(),
)
BASE.EXPECTED_BINARY_SHA256 = {
    "candidate":
        "c28333b52978d7cc2ca8a4e09375bc8bbf04c434f759747250811363a9bfcee6",
    "control":
        "5990cdf02e192fa21759e039e96c08255dcbb362e530d073c6c8d24349c80cf9",
    "main":
        "a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910",
}
BASE.REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE = True


def campaign_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []

    def append(
        original_count: int,
        recovery_count: int,
        shard_bytes: int,
        role: str,
        selected: bool,
    ) -> None:
        index = len(cells)
        cells.append({
            "id": (
                f"{role}-k{original_count}-r{recovery_count}"
                f"-b{shard_bytes}-q1"
            ),
            "K": original_count,
            "R": recovery_count,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": BASE.CELL_REUSE,
            "role": role,
            "candidate_selected": selected,
            "control_selected": False,
            "seed": 0x74370000 + index,
        })

    for original_count in range(4, 8):
        for recovery_count in (3, 4):
            for shard_bytes in (64, 128, 256):
                append(original_count, recovery_count, shard_bytes,
                       "preexisting_regression", True)
    for recovery_count in (3, 4):
        append(4, recovery_count, 1024, "preexisting_regression", True)

    for original_count in range(4, 8):
        for recovery_count in (3, 4):
            append(original_count, recovery_count, 512,
                   "b512_target", True)

    for recovery_count in (3, 4):
        for shard_bytes in (64, 128, 256):
            append(8, recovery_count, shard_bytes, "k8_target", True)
        append(8, recovery_count, 512, "k8_b512_control", False)

    BASE.require(
        len(cells) == 42 and
        len({cell["id"] for cell in cells}) == len(cells) and
        sum(cell["role"] == "preexisting_regression" for cell in cells) ==
            26 and
        sum(cell["role"] == "b512_target" for cell in cells) == 8 and
        sum(cell["role"] == "k8_target" for cell in cells) == 6 and
        sum(cell["role"] == "k8_b512_control" for cell in cells) == 2 and
        sum(cell["candidate_selected"] for cell in cells) == 40,
        "integrated T=4 acceptance matrix is incomplete")
    return cells


BASE.campaign_cells = campaign_cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
