#!/usr/bin/env python3
"""Measure the exact-B64 GF8 T=32 final-IFFT range candidate.

This is a configured specialization of the hardened descriptive ABBA runner.
It consumes three already-built immutable executables: the ordinary candidate,
a same-source control built with
LEO2_DIAGNOSTIC_DISABLE_T32_FINAL_IFFT2_RANGE=1, and exact Leopard main.  The
shared runner verifies that candidate and control differ only in the initialized
selector object and GNU build-id note, validates embedded source identities and
full-file SHA-256 before and after execution, pins every process to CPU14 while
reserving its SMT sibling CPU30, and retains every normalized process result in
raw.json before writing a deterministic summary.json.

All 72 K/R/byte cells run candidate, control, and exact main.  The B=64 cells
are the selected candidate region; B=256 and B=1024 are inert controls.  The
two promotion-primary cells are K=64/R=32/B=64 and K=127/R=28/B=64.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any


BENCHMARK_CPU = 14
RESERVED_SIBLING = 30
ORIGINAL_COUNTS = (63, 64, 65, 96, 100, 126, 127, 128)
RECOVERY_COUNTS = (28, 31, 32)
BYTE_COUNTS = (64, 256, 1024)
PRIMARY_TARGETS = frozenset(((64, 32, 64), (127, 28, 64)))


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_t4_packed_terminal_family_abba.py")
    specification = importlib.util.spec_from_file_location(
        "t32_final_ifft_hardened_abba_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load hardened ABBA base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-t32-final-ifft-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-t32-final-ifft-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN7leopard3ff8L29g_high_final_ifft2_range_modeE"
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
        "6737262ea4b206690678c01a742c6e8e56f99cfdec2ead968f0c26aa9db69d42",
    "control":
        "b88bd5c7c8a57915ebb4cfd09e74e39ee9dcd1691f7fea1d4d16d614be691a64",
    "main":
        "a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910",
}
BASE.REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE = True


def campaign_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    index = 0
    for original_count in ORIGINAL_COUNTS:
        for recovery_count in RECOVERY_COUNTS:
            for shard_bytes in BYTE_COUNTS:
                shape = (original_count, recovery_count, shard_bytes)
                primary = shape in PRIMARY_TARGETS
                selected = shard_bytes == 64
                role = "primary_target" if primary else \
                    "selected_target" if selected else "unselected_control"
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
                    "primary_target": primary,
                    "candidate_selected": selected,
                    "control_selected": False,
                    "seed": 0x74320000 + index,
                })
                index += 1

    expected = {
        (original_count, recovery_count, shard_bytes)
        for original_count in ORIGINAL_COUNTS
        for recovery_count in RECOVERY_COUNTS
        for shard_bytes in BYTE_COUNTS
    }
    actual = {(cell["K"], cell["R"], cell["bytes"]) for cell in cells}
    BASE.require(
        len(cells) == 72 and actual == expected and
        len({cell["id"] for cell in cells}) == len(cells) and
        sum(cell["candidate_selected"] for cell in cells) == 24 and
        sum(cell["primary_target"] for cell in cells) == 2 and
        {shape for shape in actual if shape in PRIMARY_TARGETS} ==
            set(PRIMARY_TARGETS),
        "T=32 final-IFFT acceptance matrix is incomplete")
    return cells


BASE.campaign_cells = campaign_cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
