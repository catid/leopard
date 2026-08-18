#!/usr/bin/env python3
"""Qualify the exact GF8 K=66/R=16/64-byte second-tail AVX2 encoder.

The target is compared with exact Leopard main and with the same frozen
Leopard2 executable after disabling only the coordinate-81 tail update.
K65/K67, R15/R17, and 63/65-byte cells are selector-inert controls whose
complete confidence intervals must stay inside the two-percent band.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path
from typing import Any, Mapping


SUPPORT_PATH = Path(__file__).resolve().with_name(
    "run_k9r5_b1024_terminal_abba.py")


def load_support() -> Any:
    specification = importlib.util.spec_from_file_location(
        "k66r16_b64_tail_evidence_support", SUPPORT_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load qualification support: {SUPPORT_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()
BASE = SUPPORT.BASE
BASE.__doc__ = __doc__
SUPPORT.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-k66r16-b64-tail-abba/v2"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-k66r16-b64-tail-summary/v2"
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L22g_k66r16_b64_tail_modeE"
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.CANDIDATE_SCHEMA = "leopard2-benchmark-v5"
BASE.CONTROL_SCHEMA = "leopard2-benchmark-v30"
BASE.CONTROL_EXTRA_ARGUMENTS = ("--k66r16-b64-tail-mode", "0")
BASE.CONTROL_BUILD_MARKER = "k66r16_b64_tail_diagnostic_disabled"
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


# The candidate and diagnostic control are the same executable and differ only
# by a command-line selector.  Running separate byte-for-byte copies can put
# their text on different physical cache pages and create a stable false label
# effect for these sub-microsecond kernels.  Freeze the input once and execute
# both labels through that one immutable path.  Both caller inputs and hashes
# are still checked independently and must identify the same source file.
_BASE_FREEZE_EXECUTABLE = BASE.freeze_executable
_BASE_RUN_ONE = BASE.run_one
_BASE_BUILD_CLOSURE_IDENTITY = BASE.build_closure_identity
_SHARED_FROZEN_CANDIDATE: Mapping[str, Any] | None = None
_SHARED_FROZEN_PHYSICAL_IDENTITY: Mapping[str, Any] | None = None


def physical_file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    status = resolved.stat()
    BASE.require(resolved.is_file() and status.st_size > 0,
                 f"physical identity is not a nonempty file: {resolved}")
    return {
        "path": str(resolved),
        "device": status.st_dev,
        "inode": status.st_ino,
    }


def require_shared_frozen_physical_identity() -> None:
    BASE.require(_SHARED_FROZEN_PHYSICAL_IDENTITY is not None,
                 "shared frozen physical identity was not established")
    observed = physical_file_identity(Path(str(
        _SHARED_FROZEN_PHYSICAL_IDENTITY["path"])))
    BASE.require(observed == _SHARED_FROZEN_PHYSICAL_IDENTITY,
                 "shared frozen executable inode changed")


def freeze_shared_candidate_control(
    source: Path,
    expected_sha256: str | None,
    destination: Path,
) -> dict[str, Any]:
    global _SHARED_FROZEN_CANDIDATE, _SHARED_FROZEN_PHYSICAL_IDENTITY
    if destination.name == "candidate":
        BASE.require(_SHARED_FROZEN_CANDIDATE is None and
                     _SHARED_FROZEN_PHYSICAL_IDENTITY is None,
                     "shared candidate was frozen more than once")
        record = _BASE_FREEZE_EXECUTABLE(
            source, expected_sha256, destination)
        record = dict(record)
        _SHARED_FROZEN_PHYSICAL_IDENTITY = physical_file_identity(
            Path(str(record["frozen"]["path"])))
        record["frozen_physical_identity"] = dict(
            _SHARED_FROZEN_PHYSICAL_IDENTITY)
        _SHARED_FROZEN_CANDIDATE = record
        return record
    if destination.name == "control":
        BASE.require(_SHARED_FROZEN_CANDIDATE is not None,
                     "control was frozen before the shared candidate")
        input_identity = BASE.T8_SUPPORT.file_identity(source)
        BASE.require(
            expected_sha256 is None or
            (re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is not None and
             input_identity["sha256"] == expected_sha256),
            f"control input binary SHA-256 changed: {source}")
        BASE.require(
            input_identity == _SHARED_FROZEN_CANDIDATE["input"],
            "candidate/control inputs do not name one physical executable")
        frozen_identity = dict(_SHARED_FROZEN_CANDIDATE["frozen"])
        BASE.require(
            Path(str(frozen_identity["path"])).parent ==
                destination.parent.resolve(strict=True) and
            not destination.exists(),
            "shared control destination is ambiguous")
        return {
            "input": input_identity,
            "frozen": frozen_identity,
            "shared_physical_executable": "candidate",
            "frozen_physical_identity": dict(
                _SHARED_FROZEN_PHYSICAL_IDENTITY),
        }
    BASE.require(destination.name == "main",
                 f"unexpected executable freeze label: {destination.name}")
    return _BASE_FREEZE_EXECUTABLE(source, expected_sha256, destination)


def run_one_from_shared_frozen_executable(
    implementation: str,
    identity: Mapping[str, Any],
    cell: Mapping[str, Any],
    cpu: int,
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
    failure_output: Path,
) -> dict[str, Any]:
    if implementation in {"candidate", "control"}:
        require_shared_frozen_physical_identity()
        BASE.require(
            identity.get("path") ==
                _SHARED_FROZEN_PHYSICAL_IDENTITY["path"],
            "candidate/control invocation escaped the shared executable")
    result = _BASE_RUN_ONE(
        implementation, identity, cell, cpu, source_commit, source_tree,
        iterations, warmup, failure_output)
    if implementation in {"candidate", "control"}:
        require_shared_frozen_physical_identity()
    return result


def build_closure_with_shared_frozen_executable(options: Any) -> Any:
    require_shared_frozen_physical_identity()
    result = _BASE_BUILD_CLOSURE_IDENTITY(options)
    require_shared_frozen_physical_identity()
    return result


BASE.freeze_executable = freeze_shared_candidate_control
BASE.run_one = run_one_from_shared_frozen_executable
BASE.build_closure_identity = build_closure_with_shared_frozen_executable


def cells() -> list[dict[str, Any]]:
    values = (
        ("target-k66-r16-b64-q5", 66, 16, 64, "target", True),
        ("k-control-k65-r16-b64-q5", 65, 16, 64, "neighbor", False),
        ("k-control-k67-r16-b64-q5", 67, 16, 64, "neighbor", False),
        ("r-control-k66-r15-b64-q5", 66, 15, 64, "neighbor", False),
        ("r-control-k66-r17-b64-q3", 66, 17, 64, "neighbor", False),
        ("byte-control-k66-r16-b63-q5", 66, 16, 63, "neighbor", False),
        ("byte-control-k66-r16-b65-q5", 66, 16, 65, "neighbor", False),
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
            "seed": 0x42100040 + index,
        })
    BASE.require(
        len(result) == 7 and
        sum(cell["role"] == "target" for cell in result) == 1 and
        len({cell["id"] for cell in result}) == len(result) and
        len({cell["seed"] for cell in result}) == len(result) and
        {cell["id"] for cell in result if cell["compare_main"]} == {
            "target-k66-r16-b64-q5",
        } and
        all(cell["loss"] == 1 and cell["batch"] == 1 and
            cell["reuse"] == 8192 and cell["measure_one_shot"] is True
            for cell in result),
        "K66/R16/B64 qualification matrix is incomplete")
    return result


BASE.cells = cells


if __name__ == "__main__":
    raise SystemExit(BASE.main())
