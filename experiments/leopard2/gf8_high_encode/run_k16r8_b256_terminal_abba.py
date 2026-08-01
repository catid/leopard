#!/usr/bin/env python3
"""Qualify the ordinary K=16/R=8/256-byte AVX2 encode terminal.

This narrow campaign reuses the packed-terminal evidence engine: immutable
same-source binaries, identical executable sections, an explicit initialized
selector word, singleton CPU pinning, an SMT-pair lease, exact Leopard main on
the target, and candidate/control regression neighbors everywhere the terminal
must remain inert.
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
        "k16r8_b256_packed_terminal_support", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load packed-terminal support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()
SUPPORT.SCHEMA = "leopard2-gf8-k16r8-b256-terminal-abba/v1"
SUPPORT.SUMMARY_SCHEMA = "leopard2-gf8-k16r8-b256-terminal-summary/v1"
SUPPORT.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L26g_k16r8_b256_terminal_modeE"
SUPPORT.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
SUPPORT.CONTROL_EXTRA_ARGUMENTS = ("--disable-k16r8-b256-terminal",)
SUPPORT.CONTROL_BUILD_MARKER = \
    "k16r8_b256_terminal_diagnostic_disabled"
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
        ("target-k16-r8-b256-q1", 16, 8, 256, 1, "target"),
        ("byte-neighbor-k16-r8-b255-q1", 16, 8, 255, 1, "neighbor"),
        ("byte-neighbor-k16-r8-b257-q1", 16, 8, 257, 1, "neighbor"),
        ("byte-neighbor-k16-r8-b512-q1", 16, 8, 512, 1, "neighbor"),
        ("shape-neighbor-k15-r8-b256-q1", 15, 8, 256, 1, "neighbor"),
        ("shape-neighbor-k16-r7-b256-q1", 16, 7, 256, 1, "neighbor"),
        ("shape-neighbor-k15-r7-b256-q1", 15, 7, 256, 1, "neighbor"),
        ("existing-terminal-k5-r5-b64-q1", 5, 5, 64, 1, "neighbor"),
        ("terminal-neighbor-k1-r1-b256-q1", 1, 1, 256, 1, "neighbor"),
        ("terminal-neighbor-k2-r1-b256-q1", 2, 1, 256, 1, "neighbor"),
        ("batch-neighbor-k16-r8-b256-q2", 16, 8, 256, 2, "neighbor"),
        ("batch-neighbor-k16-r8-b256-q8", 16, 8, 256, 8, "neighbor"),
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
            "seed": 0x16825600 + index,
        })
    return result


SUPPORT.cells = cells


if __name__ == "__main__":
    raise SystemExit(SUPPORT.main())
