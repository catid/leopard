#!/usr/bin/env python3
"""Early immutable-binary screen for the generated T=32/B=256 encoder.

This specialization compares a generated-kernel build with an
executable-section-identical same-source control and exact Leopard main.  It
is intentionally target-only: neighboring/layout and structural build gates
remain mandatory before promotion if this arithmetic screen succeeds.
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
BASE.SCHEMA = "leopard2-gf8-t32-b256-generated-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-t32-b256-generated-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN7leopard7backend12_GLOBAL__N_1L25g_t32_b256_generated_modeE"
BASE.TARGET_CONTROL_FLOOR = 1.05
BASE.TARGET_MAIN_FLOOR = 1.05
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
)


def cells() -> list[dict[str, Any]]:
    return [{
        "id": "target-k32-r32-b256-q1",
        "K": 32,
        "R": 32,
        "bytes": 256,
        "batch": 1,
        "reuse": 8192,
        "role": "target",
        "seed": 0x54324232,
    }]


BASE.cells = cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
