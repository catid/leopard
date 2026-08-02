#!/usr/bin/env python3
"""Qualify balanced T=32/64/128 GF8 AVX2 64-byte encode terminals.

This configured specialization of the hardened fixed-terminal ABBA runner
compares a frozen candidate with an executable-section-identical same-source
control and exact Leopard main.  It retains byte and immediate K/R neighbors,
the rejected T=16 shape, and unrelated prior tiny terminals as same-source
layout controls.  T=128 additionally requires a 1.05 lower confidence bound
against exact main before its selector can remain enabled.
"""

from __future__ import annotations

import importlib.util
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Sequence


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_k5r5_b64_terminal_abba.py")
    specification = importlib.util.spec_from_file_location(
        "balanced_b64_terminal_evidence_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load evidence base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-balanced-b64-terminal-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-balanced-b64-terminal-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L28g_balanced_b64_terminal_modeE"
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.TARGET_MAIN_FLOOR = 1.05
BASE.TARGET_ORDER = BASE.TARGET_ORDER * 3
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
)


def confidence_interval(values: Sequence[float]) -> dict[str, Any]:
    """Student interval over three, four, or nine independent ABBA rounds."""
    BASE.require(len(values) in (3, 4, 9),
                 "three, four, or nine round contrasts are required")
    critical = {
        3: BASE.T95_DF2,
        4: BASE.T95_DF3,
        9: 2.306004135204166,
    }[len(values)]
    center = statistics.mean(values)
    half_width = critical * statistics.stdev(values) / math.sqrt(len(values))
    return {
        "speedup": math.exp(center),
        "ci95": [
            math.exp(center - half_width),
            math.exp(center + half_width),
        ],
        "round_log_ratios": list(values),
    }


BASE.confidence_interval = confidence_interval


def cells() -> list[dict[str, Any]]:
    values = (
        ("target-k32-r32-b64-q1", 32, 32, 64, "target"),
        ("target-k64-r64-b64-q1", 64, 64, 64, "target"),
        ("target-k128-r128-b64-q1", 128, 128, 64, "target"),
        ("byte-k32-r32-b63-q1", 32, 32, 63, "neighbor"),
        ("byte-k32-r32-b65-q1", 32, 32, 65, "neighbor"),
        ("byte-k32-r32-b256-q1", 32, 32, 256, "neighbor"),
        ("byte-k64-r64-b63-q1", 64, 64, 63, "neighbor"),
        ("byte-k64-r64-b65-q1", 64, 64, 65, "neighbor"),
        ("byte-k64-r64-b256-q1", 64, 64, 256, "neighbor"),
        ("byte-k128-r128-b63-q1", 128, 128, 63, "neighbor"),
        ("byte-k128-r128-b65-q1", 128, 128, 65, "neighbor"),
        ("byte-k128-r128-b256-q1", 128, 128, 256, "neighbor"),
        ("shape-k31-r32-b64-q1", 31, 32, 64, "neighbor"),
        ("shape-k32-r31-b64-q1", 32, 31, 64, "neighbor"),
        ("shape-k63-r64-b64-q1", 63, 64, 64, "neighbor"),
        ("shape-k64-r63-b64-q1", 64, 63, 64, "neighbor"),
        ("shape-k127-r128-b64-q1", 127, 128, 64, "neighbor"),
        ("shape-k128-r127-b64-q1", 128, 127, 64, "neighbor"),
        ("rejected-k16-r16-b64-q1", 16, 16, 64, "neighbor"),
        ("terminal-k8-r8-b64-q1", 8, 8, 64, "neighbor"),
        ("terminal-k5-r5-b64-q1", 5, 5, 64, "neighbor"),
        ("terminal-k16-r8-b256-q1", 16, 8, 256, "neighbor"),
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
            "seed": 0x42363400 + index,
        }
        for index, (
            name, original_count, recovery_count, shard_bytes, role
        ) in enumerate(values)
    ]


BASE.cells = cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
