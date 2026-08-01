#!/usr/bin/env python3
"""Qualify the ordinary K=9/R=5/256-byte AVX2 encode terminal.

This narrow campaign reuses the packed-terminal evidence engine: one immutable
same-source binary, an initialized selector word, singleton CPU pinning, an
SMT-pair lease, exact Leopard main on the target, and candidate/control
regression neighbors everywhere the terminal must remain inert.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any


def load_support() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_k5r5_b64_terminal_abba.py")
    specification = importlib.util.spec_from_file_location(
        "k9r5_b256_packed_terminal_support", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load packed-terminal support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()
SUPPORT.SCHEMA = "leopard2-gf8-k9r5-b256-terminal-abba/v1"
SUPPORT.SUMMARY_SCHEMA = "leopard2-gf8-k9r5-b256-terminal-summary/v1"
SUPPORT.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L25g_k9r5_b256_terminal_modeE"
SUPPORT.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
SUPPORT.CONTROL_EXTRA_ARGUMENTS = ("--disable-k9r5-b256-terminal",)
SUPPORT.CONTROL_BUILD_MARKER = \
    "k9r5_b256_terminal_diagnostic_disabled"
SUPPORT.NEIGHBOR_ORDER = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
)
SUPPORT.__doc__ = __doc__
SUPPORT.__file__ = __file__


def cells() -> list[dict[str, Any]]:
    values = [
        ("target-k9-r5-b256-q1", 9, 5, 256, 1, "target"),
        ("byte-neighbor-k9-r5-b255-q1", 9, 5, 255, 1, "neighbor"),
        ("byte-neighbor-k9-r5-b257-q1", 9, 5, 257, 1, "neighbor"),
        ("byte-neighbor-k9-r5-b512-q1", 9, 5, 512, 1, "neighbor"),
        ("shape-neighbor-k8-r5-b256-q1", 8, 5, 256, 1, "neighbor"),
        ("shape-neighbor-k9-r4-b256-q1", 9, 4, 256, 1, "neighbor"),
        ("shape-neighbor-k10-r5-b256-q1", 10, 5, 256, 1, "neighbor"),
        ("shape-neighbor-k9-r6-b256-q1", 9, 6, 256, 1, "neighbor"),
        ("existing-terminal-k5-r5-b64-q1", 5, 5, 64, 1, "neighbor"),
        ("existing-terminal-k16-r8-b256-q1", 16, 8, 256, 1, "neighbor"),
        ("batch-neighbor-k9-r5-b256-q2", 9, 5, 256, 2, "neighbor"),
        ("batch-neighbor-k9-r5-b256-q8", 9, 5, 256, 8, "neighbor"),
    ]
    result = []
    for index, (name, k, r, shard_bytes, batch, role) in enumerate(values):
        result.append({
            "id": name,
            "K": k,
            "R": r,
            "bytes": shard_bytes,
            "batch": batch,
            "reuse": max(8192, 65536 // batch),
            "role": role,
            "seed": 0x09525600 + index,
        })
    return result


SUPPORT.cells = cells


if __name__ == "__main__":
    raise SystemExit(SUPPORT.main())
