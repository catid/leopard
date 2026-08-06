#!/usr/bin/env python3
"""Qualify the packed B256 T=8 tail/two-block encode terminal.

The campaign compares one immutable Leopard2 executable with its terminal
enabled and disabled through a process-local selector.  Every promoted shape
is also compared with the exact Leopard main codec.  Excluded K/R cells and
byte/count boundaries prove that the selector is inert outside the measured
ownership region.
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
        "t8_two_block_b256_terminal_evidence_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load evidence base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-t8-two-block-b256-terminal-abba/v1"
BASE.SUMMARY_SCHEMA = \
    "leopard2-gf8-t8-two-block-b256-terminal-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L36g_high_t8_two_block_b256_packed_modeE"
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.CONTROL_EXTRA_ARGUMENTS = (
    "--high-t8-two-block-b256-terminal-mode", "0")
BASE.CONTROL_BUILD_MARKER = \
    "high_t8_two_block_b256_terminal_diagnostic_disabled"
BASE.CONTROL_SCHEMA = "leopard2-benchmark-v26"
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


def promoted(original_count: int, recovery_count: int) -> bool:
    return (
        original_count == 10 or
        (original_count == 11 and recovery_count == 5) or
        (original_count >= 13 and
         not (original_count == 16 and recovery_count == 8))
    )


def cells() -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for original_count in range(10, 17):
        for recovery_count in range(5, 9):
            if not promoted(original_count, recovery_count):
                continue
            index = len(result)
            result.append({
                "id": (
                    f"target-k{original_count}-r{recovery_count}"
                    "-b256-q1"),
                "K": original_count,
                "R": recovery_count,
                "bytes": 256,
                "batch": 1,
                "reuse": 8192,
                "role": "target",
                "measure_one_shot": True,
                "seed": 0x28374000 + index,
            })

    boundaries = (
        ("excluded-k11-r6", 11, 6, 256),
        ("excluded-k11-r7", 11, 7, 256),
        ("excluded-k11-r8", 11, 8, 256),
        ("excluded-k12-r5", 12, 5, 256),
        ("excluded-k12-r6", 12, 6, 256),
        ("excluded-k12-r7", 12, 7, 256),
        ("excluded-k12-r8", 12, 8, 256),
        ("existing-k9-r5", 9, 5, 256),
        ("existing-k16-r8", 16, 8, 256),
        ("byte-k10-r5-b255", 10, 5, 255),
        ("byte-k10-r5-b257", 10, 5, 257),
        ("byte-k16-r7-b255", 16, 7, 255),
        ("byte-k16-r7-b257", 16, 7, 257),
        ("shape-k10-r4", 10, 4, 256),
        ("shape-k17-r8", 17, 8, 256),
        ("prior-b64-k10-r5", 10, 5, 64),
        ("prior-b64-k16-r7", 16, 7, 64),
        ("size-k13-r5-b128", 13, 5, 128),
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
            "seed": 0x28374000 + index,
        })
    return result


BASE.cells = cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
