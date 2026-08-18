#!/usr/bin/env python3
"""Qualify the exact GF8 K=62/R=8/64-byte fused AVX2 encoder.

The target is compared with exact Leopard main and with the same frozen
Leopard2 executable after disabling only the K62/R8/B64 arithmetic leaf.
Adjacent K, R, and byte-count cells must prove selector equivalence inside a
two-percent band.  In particular, the dedicated selector leaves the R=9 T16
terminal enabled, avoiding the confounded family-wide control used during
exploration.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any


SUPPORT_PATH = Path(__file__).resolve().with_name(
    "run_k9r5_b1024_terminal_abba.py")


def load_support() -> Any:
    specification = importlib.util.spec_from_file_location(
        "k62r8_b64_fused_evidence_support", SUPPORT_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load qualification support: {SUPPORT_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()
BASE = SUPPORT.BASE
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-k62r8-b64-fused-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-k62r8-b64-fused-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L22g_k62r8_b64_fused_modeE"
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.CANDIDATE_SCHEMA = "leopard2-benchmark-v5"
BASE.CONTROL_SCHEMA = "leopard2-benchmark-v29"
BASE.CONTROL_EXTRA_ARGUMENTS = ("--k62r8-b64-fused-mode", "0")
BASE.CONTROL_BUILD_MARKER = "k62r8_b64_fused_diagnostic_disabled"
BASE.TARGET_CONTROL_FLOOR = 1.05
BASE.TARGET_MAIN_FLOOR = 1.05
BASE.NEIGHBOR_EQUIVALENCE_LOWER = 1.0 / 1.02
BASE.NEIGHBOR_EQUIVALENCE_UPPER = 1.02
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    SUPPORT_PATH.resolve(),
    SUPPORT.BASE_PATH.resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
    SUPPORT.PROVENANCE_PATH.resolve(),
)


def cells() -> list[dict[str, Any]]:
    values = (
        ("target-k62-r8-b64-q1", 62, 8, 64, "target", True),
        ("k-control-k61-r8-b64-q1", 61, 8, 64, "neighbor", False),
        ("k-control-k63-r8-b64-q1", 63, 8, 64, "neighbor", False),
        ("r-control-k62-r7-b64-q1", 62, 7, 64, "neighbor", False),
        ("r-control-k62-r9-b64-q1", 62, 9, 64, "neighbor", False),
        ("byte-control-k62-r8-b63-q1", 62, 8, 63, "neighbor", False),
        ("byte-control-k62-r8-b65-q1", 62, 8, 65, "neighbor", False),
    )
    result = []
    for index, (name, k, r, shard_bytes, role, compare_main) in enumerate(
            values):
        result.append({
            "id": name,
            "K": k,
            "R": r,
            "bytes": shard_bytes,
            "loss": 1,
            "batch": 1,
            "reuse": 8192,
            "role": role,
            "compare_main": compare_main,
            "measure_one_shot": True,
            "seed": 0x3E080040 + index,
        })
    BASE.require(
        len(result) == 7 and
        sum(cell["role"] == "target" for cell in result) == 1 and
        len({cell["id"] for cell in result}) == len(result) and
        len({cell["seed"] for cell in result}) == len(result) and
        {cell["id"] for cell in result if cell["compare_main"]} == {
            "target-k62-r8-b64-q1",
        } and
        all(cell["loss"] == 1 and cell["batch"] == 1 and
            cell["reuse"] == 8192 and cell["measure_one_shot"] is True
            for cell in result),
        "K62/R8/B64 qualification matrix is incomplete")
    return result


BASE.cells = cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
