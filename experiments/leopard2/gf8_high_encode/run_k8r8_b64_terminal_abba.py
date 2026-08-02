#!/usr/bin/env python3
"""Qualify the ordinary K=8/R=8/64-byte AVX2 encode terminal.

This is a configured specialization of the hardened fixed-terminal ABBA
runner.  It compares one frozen candidate with an executable-section-identical
same-source compile-time control and exact Leopard main.  Byte, shape, and
larger-transform neighbors use the same-source contrast to detect code-layout
or added-dispatch regressions without conflating them with historical main.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_k5r5_b64_terminal_abba.py")
    specification = importlib.util.spec_from_file_location(
        "k8r8_b64_terminal_evidence_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load evidence base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-k8r8-b64-terminal-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-k8r8-b64-terminal-summary/v1"
BASE.MODE_SYMBOL = "_ZN12_GLOBAL__N_1L24g_k8r8_b64_terminal_modeE"
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
)


def cells() -> list[dict[str, Any]]:
    values = (
        ("target-k8-r8-b64-q1", 8, 8, 64, "target"),
        ("byte-neighbor-k8-r8-b63-q1", 8, 8, 63, "neighbor"),
        ("byte-neighbor-k8-r8-b65-q1", 8, 8, 65, "neighbor"),
        ("byte-neighbor-k8-r8-b128-q1", 8, 8, 128, "neighbor"),
        ("shape-neighbor-k8-r7-b64-q1", 8, 7, 64, "neighbor"),
        ("shape-neighbor-k7-r8-b64-q1", 7, 8, 64, "neighbor"),
        ("shape-neighbor-k9-r8-b64-q1", 9, 8, 64, "neighbor"),
        ("scale-neighbor-k64-r32-b64-q1", 64, 32, 64, "neighbor"),
        ("scale-neighbor-k100-r28-b64-q1", 100, 28, 64, "neighbor"),
    )
    return [
        {
            "id": name,
            "K": original_count,
            "R": recovery_count,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": 8192,
            "role": role,
            "seed": 0x4B385200 + index,
        }
        for index, (
            name, original_count, recovery_count, shard_bytes, role
        ) in enumerate(values)
    ]


BASE.cells = cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
