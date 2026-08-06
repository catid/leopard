#!/usr/bin/env python3
"""Qualify the packed K=9..16/R=5..8/B=64 T=8 encode terminal.

The campaign compares one immutable Leopard2 executable with its packed
terminal enabled by default and disabled through a process-local diagnostic
selector.  Every promoted shape is also compared with the exact Leopard main
codec.  Byte and K/R boundary cells prove that the selector is inert outside
the intended region.
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
        "t8_two_block_b64_terminal_evidence_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load evidence base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-t8-two-block-b64-terminal-abba/v1"
BASE.SUMMARY_SCHEMA = \
    "leopard2-gf8-t8-two-block-b64-terminal-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L35g_high_t8_two_block_b64_packed_modeE"
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.CONTROL_EXTRA_ARGUMENTS = (
    "--high-t8-two-block-b64-terminal-mode", "0")
BASE.CONTROL_BUILD_MARKER = \
    "high_t8_two_block_b64_terminal_diagnostic_disabled"
BASE.CONTROL_SCHEMA = "leopard2-benchmark-v25"
BASE.TARGET_CONTROL_FLOOR = 1.05
BASE.TARGET_MAIN_FLOOR = 1.05
BASE.REQUIRE_EXPECTED_IDENTITIES = True
BASE.CHILD_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_THREAD_LIMIT": "1",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
)


def cells() -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for original_count in range(9, 17):
        for recovery_count in range(5, 9):
            index = len(result)
            result.append({
                "id": (
                    f"target-k{original_count}-r{recovery_count}"
                    "-b64-q1"),
                "K": original_count,
                "R": recovery_count,
                "bytes": 64,
                "batch": 1,
                "reuse": 8192,
                "role": "target",
                "measure_one_shot": True,
                "seed": 0x28364000 + index,
            })

    boundaries = (
        ("byte-k9-r5-b63", 9, 5, 63),
        ("byte-k9-r5-b65", 9, 5, 65),
        ("byte-k16-r8-b63", 16, 8, 63),
        ("byte-k16-r8-b65", 16, 8, 65),
        ("shape-k8-r8-b64", 8, 8, 64),
        ("shape-k9-r4-b64", 9, 4, 64),
        ("shape-k16-r4-b64", 16, 4, 64),
        ("shape-k17-r8-b64", 17, 8, 64),
        ("size-k9-r5-b128", 9, 5, 128),
        ("size-k16-r8-b256", 16, 8, 256),
    )
    for name, original_count, recovery_count, shard_bytes in boundaries:
        index = len(result)
        result.append({
            "id": f"neighbor-{name}-q1",
            "K": original_count,
            "R": recovery_count,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": 8192,
            "role": "neighbor",
            "measure_one_shot": True,
            "seed": 0x28364000 + index,
        })
    return result


BASE.cells = cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
