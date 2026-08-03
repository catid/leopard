#!/usr/bin/env python3
"""Final immutable-binary qualification for the T=32/B=256 encoder.

This specialization compares a generated-kernel build with an
executable-section-identical same-source control and exact Leopard main.  It
times the exact one-stripe target across practical erasure counts, then proves
byte, shape, and batch neighbors do not regress.  Main-compatible cells retain
exact-main observations; R>K shapes remain same-source controls because the
legacy API rejects them.  Only exact K=R=32/B=256/q1 cells are held to the
five-percent candidate promotion floor.
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
        "t32_b256_generated_evidence_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load evidence base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-t32-b256-generated-abba/v4"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-t32-b256-generated-summary/v4"
BASE.MODE_SYMBOL = \
    "_ZN7leopard7backend12_GLOBAL__N_1L25g_t32_b256_generated_modeE"
BASE.TARGET_CONTROL_FLOOR = 1.05
BASE.TARGET_MAIN_FLOOR = 1.05
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.REQUIRE_EXPECTED_IDENTITIES = True
BASE.REQUIRE_BUILD_CLOSURE = True
BASE.REQUIRE_FULL_ELF_IDENTITY = True
BASE.TARGET_ORDER = BASE.TARGET_ORDER * 3
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
)


def cells() -> list[dict[str, Any]]:
    values = []
    for loss in (1, 2, 4, 8, 16, 32):
        values.append((
            f"target-k32-r32-b256-l{loss}-q1",
            32, 32, 256, loss, 1, "target"))
    values.append((
        "batch-neighbor-k32-r32-b256-l1-q8",
        32, 32, 256, 1, 8, "neighbor"))
    values.extend((
        ("byte-neighbor-k32-r32-b255-l1-q1",
         32, 32, 255, 1, 1, "neighbor"),
        ("byte-neighbor-k32-r32-b257-l1-q1",
         32, 32, 257, 1, 1, "neighbor"),
        ("byte-neighbor-k32-r32-b64-l1-q1",
         32, 32, 64, 1, 1, "neighbor"),
        ("byte-neighbor-k32-r32-b1024-l1-q1",
         32, 32, 1024, 1, 1, "neighbor"),
        ("shape-neighbor-k31-r32-b256-l1-q1",
         31, 32, 256, 1, 1, "neighbor"),
        ("shape-neighbor-k32-r31-b256-l1-q1",
         32, 31, 256, 1, 1, "neighbor"),
        ("shape-neighbor-k31-r31-b256-l1-q1",
         31, 31, 256, 1, 1, "neighbor"),
        ("shape-neighbor-k33-r32-b256-l1-q1",
         33, 32, 256, 1, 1, "neighbor"),
        ("shape-neighbor-k32-r33-b256-l1-q1",
         32, 33, 256, 1, 1, "neighbor"),
        ("side-neighbor-k16-r16-b256-l1-q1",
         16, 16, 256, 1, 1, "neighbor"),
        ("side-neighbor-k64-r64-b256-l1-q1",
         64, 64, 256, 1, 1, "neighbor"),
        ("practical-neighbor-k16-r8-b256-l4-q1",
         16, 8, 256, 4, 1, "neighbor"),
    ))
    result = []
    for index, (name, k, r, shard_bytes, loss, batch, role) in enumerate(
            values):
        result.append({
            "id": name,
            "K": k,
            "R": r,
            "bytes": shard_bytes,
            "loss": loss,
            "batch": batch,
            "measure_one_shot": batch == 1,
            "reuse": 8192,
            "role": role,
            "compare_main": r <= k,
            "seed": 0x54324200 + index,
        })
    BASE.require(
        len(result) == 19 and
        sum(cell["role"] == "target" for cell in result) == 6 and
        len({cell["id"] for cell in result}) == len(result) and
        len({cell["seed"] for cell in result}) == len(result) and
        all(1 <= cell["loss"] <= cell["R"] for cell in result) and
        all(cell["measure_one_shot"] == (cell["batch"] == 1)
            for cell in result) and
        {cell["id"] for cell in result if not cell["compare_main"]} == {
            "shape-neighbor-k31-r32-b256-l1-q1",
            "shape-neighbor-k32-r33-b256-l1-q1",
        },
        "T32/B256 final qualification matrix is incomplete")
    return result


BASE.cells = cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
