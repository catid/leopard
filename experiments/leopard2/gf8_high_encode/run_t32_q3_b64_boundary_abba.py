#!/usr/bin/env python3
"""Qualify the GF8 T=R=32, K=65..96, 64-byte packed boundary.

This is an intentionally narrow production-promotion campaign.  It qualifies
only public one-shot encode and a one-item batch call at the 64-byte boundary;
it makes no claim about multi-stripe batching or other shard sizes.  K=77,
78, and 79 are the regression focus, while K=65, 76, 80, 81, 95, and 96
prove that the aggregate route is useful across its complete-block and ragged
tail boundaries.  Every active cell must clear a five-percent lower confidence
bound against the same frozen production executable with the route disabled,
and must not regress against the canonical exact Leopard-main executable.

Candidate and control are two command-line modes of one immutable inode.  The
five inactive K/R/byte neighbors must prove that both execution contrasts have
their complete confidence intervals inside the two-percent equivalence band.
All invocations retain and compare original, parity, and recovery digests.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


PARENT_PATH = Path(__file__).resolve().with_name(
    "run_k66r16_b64_tail_abba.py")


def load_parent() -> Any:
    specification = importlib.util.spec_from_file_location(
        "t32_q3_b64_boundary_evidence_parent", PARENT_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load qualification support: {PARENT_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


PARENT = load_parent()
BASE = PARENT.BASE
SUPPORT_DEPENDENCIES = tuple(BASE.RUNNER_DEPENDENCIES)
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-t32-q3-b64-boundary-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-t32-q3-b64-boundary-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L28g_balanced_b64_terminal_modeE"
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.CANDIDATE_SCHEMA = "leopard2-benchmark-v22"
BASE.CONTROL_SCHEMA = "leopard2-benchmark-v22"
BASE.CONTROL_EXTRA_ARGUMENTS = (
    "--balanced-b64-terminal-mode", "0",
)
BASE.CONTROL_BUILD_MARKER = None
BASE.TARGET_CONTROL_FLOOR = 1.05
BASE.TARGET_MAIN_FLOOR = 1.0
BASE.NEIGHBOR_EQUIVALENCE_LOWER = 1.0 / 1.02
BASE.NEIGHBOR_EQUIVALENCE_UPPER = 1.02
BASE.CANONICAL_MAIN_SHA256 = \
    "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93"
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    *SUPPORT_DEPENDENCIES,
)


def cells() -> list[dict[str, Any]]:
    values = (
        ("target-k77-r32-b64-q3", 77, "focus", "target", True),
        ("target-k78-r32-b64-q3", 78, "focus", "target", True),
        ("target-k79-r32-b64-q3", 79, "focus", "target", True),
        ("context-k65-r32-b64-q3", 65, "context", "target", True),
        ("context-k76-r32-b64-q3", 76, "context", "target", True),
        ("context-k80-r32-b64-q3", 80, "context", "target", True),
        ("context-k81-r32-b64-q3", 81, "context", "target", True),
        ("context-k95-r32-b64-q3", 95, "context", "target", True),
        ("context-k96-r32-b64-q3", 96, "context", "target", True),
        ("r-neighbor-k79-r31-b64-q3", 79, "inactive", "neighbor", False),
        ("r-neighbor-k79-r33-b64-q2", 79, "inactive", "neighbor", False),
        ("byte-neighbor-k79-r32-b63-q3", 79, "inactive", "neighbor", False),
        ("byte-neighbor-k79-r32-b65-q3", 79, "inactive", "neighbor", False),
        ("k-neighbor-k97-r32-b64-q4", 97, "inactive", "neighbor", False),
    )
    result = []
    for index, (name, k, qualification, role, compare_main) in enumerate(
            values):
        recovery_count = 31 if "r31" in name else \
            33 if "r33" in name else 32
        shard_bytes = 63 if "b63" in name else \
            65 if "b65" in name else 64
        result.append({
            "id": name,
            "K": k,
            "R": recovery_count,
            "bytes": shard_bytes,
            "loss": 1,
            "batch": 1,
            "reuse": 8192,
            "role": role,
            "qualification": qualification,
            "compare_main": compare_main,
            "measure_one_shot": True,
            "seed": 0x54334200 + index,
        })
    active_k = {65, 76, 77, 78, 79, 80, 81, 95, 96}
    inactive_shapes = {
        (79, 31, 64),
        (79, 33, 64),
        (79, 32, 63),
        (79, 32, 65),
        (97, 32, 64),
    }
    BASE.require(
        len(result) == 14 and
        {cell["K"] for cell in result if cell["role"] == "target"} ==
            active_k and
        {(cell["K"], cell["R"], cell["bytes"]) for cell in result
         if cell["role"] == "neighbor"} == inactive_shapes and
        {cell["K"] for cell in result
         if cell["qualification"] == "focus"} == {77, 78, 79} and
        all(cell["compare_main"] == (cell["role"] == "target")
            for cell in result) and
        len({cell["id"] for cell in result}) == len(result) and
        len({cell["seed"] for cell in result}) == len(result) and
        all(cell["loss"] == 1 and cell["batch"] == 1 and
            cell["reuse"] == 8192 and
            cell["measure_one_shot"] is True for cell in result),
        "T32/Q3/B64 qualification matrix is incomplete")
    return result


BASE.cells = cells


_BASE_BENCHMARK_COMMAND = BASE.benchmark_command


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    """Select both same-inode labels explicitly; never rely on defaults."""
    command = _BASE_BENCHMARK_COMMAND(
        implementation, executable, cell, cpu, iterations, warmup)
    if implementation == "candidate":
        json_index = command.index("--json")
        command[json_index:json_index] = [
            "--balanced-b64-terminal-mode", "1"]
    if implementation in {"candidate", "control"}:
        option_indices = [
            index for index, value in enumerate(command)
            if value == "--balanced-b64-terminal-mode"
        ]
        expected = "1" if implementation == "candidate" else "0"
        BASE.require(
            len(option_indices) == 1 and
            option_indices[0] + 1 < len(command) and
            command[option_indices[0] + 1] == expected,
            f"{implementation} did not select its exact diagnostic mode")
    return command


BASE.benchmark_command = benchmark_command


_BASE_VALIDATE_RESULT = BASE.validate_result


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    normalized = _BASE_VALIDATE_RESULT(
        implementation, result, cell, source_commit, source_tree,
        iterations, warmup)
    if implementation in {"candidate", "control"}:
        BASE.require(isinstance(result, Mapping),
                     "Leopard2 diagnostic result is malformed")
        build = result.get("build")
        expected_mode = 1 if implementation == "candidate" else 0
        BASE.require(
            isinstance(build, Mapping) and
            build.get("balanced_b64_terminal_diagnostic_mode") ==
                expected_mode and
            build.get("balanced_b64_terminal_enabled") is
                (expected_mode == 1),
            f"{implementation} output does not attest its exact balanced "
            "B64 diagnostic mode")
    return normalized


BASE.validate_result = validate_result


_BASE_ANALYZE = BASE.analyze


def _require_ratio_lower_bound(
    result: Mapping[str, Any], key: str, floor: float, cell_id: str,
) -> float:
    ratio = result.get(key)
    BASE.require(isinstance(ratio, Mapping) and
                 isinstance(ratio.get("ci95"), list) and
                 len(ratio["ci95"]) == 2,
                 f"{cell_id} lacks the {key} confidence interval")
    lower = ratio["ci95"][0]
    BASE.require(isinstance(lower, (int, float)) and lower >= floor,
                 f"{cell_id} {key} lower CI {lower!r} < {floor}")
    return float(lower)


def analyze(
    cell: Mapping[str, Any], rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    result = _BASE_ANALYZE(cell, rounds)
    cell_id = str(cell["id"])
    if cell.get("role") == "target":
        lower_bounds = {
            "encode_control": _require_ratio_lower_bound(
                result, "control_over_candidate",
                BASE.TARGET_CONTROL_FLOOR, cell_id),
            "one_shot_control": _require_ratio_lower_bound(
                result, "one_shot_control_over_candidate",
                BASE.TARGET_CONTROL_FLOOR, cell_id),
            "encode_main": _require_ratio_lower_bound(
                result, "main_over_candidate",
                BASE.TARGET_MAIN_FLOOR, cell_id),
            "one_shot_main": _require_ratio_lower_bound(
                result, "one_shot_main_over_candidate",
                BASE.TARGET_MAIN_FLOOR, cell_id),
        }
        result["active_boundary_validation"] = {
            "same_binary_floor": BASE.TARGET_CONTROL_FLOOR,
            "exact_main_floor": BASE.TARGET_MAIN_FLOOR,
            "lower_ci95": lower_bounds,
            "scope": "one-shot and one-item batch at exactly 64 bytes",
            "accepted": True,
        }
        return result

    confidence_intervals = {}
    for name, key in (
        ("encode_execution", "control_over_candidate"),
        ("one_shot_encode", "one_shot_control_over_candidate"),
    ):
        ratio = result.get(key)
        BASE.require(isinstance(ratio, Mapping) and
                     isinstance(ratio.get("ci95"), list) and
                     len(ratio["ci95"]) == 2,
                     f"inactive {cell_id} lacks the {name} contrast")
        lower, upper = ratio["ci95"]
        BASE.require(
            isinstance(lower, (int, float)) and
            isinstance(upper, (int, float)) and
            lower >= BASE.NEIGHBOR_EQUIVALENCE_LOWER and
            upper <= BASE.NEIGHBOR_EQUIVALENCE_UPPER,
            f"inactive {cell_id} {name} does not prove selector equivalence "
            f"inside [{BASE.NEIGHBOR_EQUIVALENCE_LOWER}, "
            f"{BASE.NEIGHBOR_EQUIVALENCE_UPPER}]: CI [{lower}, {upper}]")
        confidence_intervals[name] = [lower, upper]
    result["inactive_boundary_validation"] = {
        "equivalence_band": [
            BASE.NEIGHBOR_EQUIVALENCE_LOWER,
            BASE.NEIGHBOR_EQUIVALENCE_UPPER,
        ],
        "ci95": confidence_intervals,
        "accepted": True,
    }
    return result


BASE.analyze = analyze


if __name__ == "__main__":
    raise SystemExit(BASE.main())
