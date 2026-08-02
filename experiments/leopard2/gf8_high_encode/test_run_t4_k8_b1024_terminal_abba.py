#!/usr/bin/env python3
"""Deterministic fail-closed tests for the K8/B1024 T=4 ABBA runner."""

from __future__ import annotations

import copy
import importlib.util
import os
import sys
from pathlib import Path
from typing import Any, Callable


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_module(filename: str, module_name: str) -> Any:
    path = Path(__file__).resolve().with_name(filename)
    specification = importlib.util.spec_from_file_location(module_name, path)
    require(specification is not None and specification.loader is not None,
            f"cannot load {filename}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_module(
    "run_t4_k8_b1024_terminal_abba.py",
    "leopard2_t4_k8_b1024_runner_test_target")
PRIOR = load_module(
    "run_t4_k8_b512_terminal_abba.py",
    "leopard2_t4_k8_b512_runner_test_reference")
BASE = RUNNER.BASE


def expect_failure(action: Callable[[], object], message: str) -> None:
    try:
        action()
    except Exception:
        return
    raise RuntimeError(message)


def passing_analyses(cells: list[dict[str, Any]]) -> list[dict[str, Any]]:
    analyses: list[dict[str, Any]] = []
    target_roles = {
        "k8_regression", "k8_b512_target", "k8_b1024_target"}
    for cell in cells:
        if cell["role"] in target_roles:
            control_speedup = 1.20
            main_speedup = 1.18
            control_ci = [1.10, 1.30]
            main_ci = [1.08, 1.28]
        else:
            control_speedup = 1.00
            main_speedup = 1.10
            control_ci = [0.99, 1.01]
            main_ci = [1.08, 1.12]
        analyses.append({
            "cell": dict(cell),
            "candidate_vs_control": {
                "speedup": control_speedup,
                "ci95": control_ci,
            },
            "candidate_vs_main": {
                "speedup": main_speedup,
                "ci95": main_ci,
            },
        })
    return analyses


def one_failure(
    analyses: list[dict[str, Any]],
    role: str,
    contrast: str,
    key: str,
    value: Any,
) -> dict[str, Any]:
    changed = copy.deepcopy(analyses)
    item = next(entry for entry in changed if entry["cell"]["role"] == role)
    item[contrast][key] = value
    result = RUNNER.evaluate_promotion(changed)
    require(result["passed"] is False and result["failure_count"] == 1,
            f"{role} synthetic failure was not isolated")
    return result["failures"][0]


def main() -> int:
    cells = RUNNER.campaign_cells()
    prior = PRIOR.campaign_cells()
    require(len(cells) == 46, "K8/B1024 matrix size changed")
    require(cells[:42] == prior[:42],
            "first 42 B512 campaign cells changed")
    require(
        [(cell["K"], cell["R"], cell["bytes"], cell["seed"])
         for cell in cells[42:44]] ==
        [(cell["K"], cell["R"], cell["bytes"], cell["seed"])
         for cell in prior[42:44]] and
        all(cell["role"] == "k8_b1024_target" and
            cell["candidate_selected"] is True and
            cell["control_selected"] is False
            for cell in cells[42:44]),
        "B1024 controls were not deterministically re-roled")
    require(
        [(cell["K"], cell["R"], cell["bytes"], cell["seed"])
         for cell in cells[44:]] == [
            (8, 3, 2048, 0x7437002c),
            (8, 4, 2048, 0x7437002d),
        ] and
        all(cell["role"] == "k8_b2048_control" and
            cell["candidate_selected"] is False and
            cell["control_selected"] is False
            for cell in cells[44:]),
        "B2048 inert controls changed")
    expected_roles = {
        "preexisting_regression": 26,
        "b512_regression": 8,
        "k8_regression": 6,
        "k8_b512_target": 2,
        "k8_b1024_target": 2,
        "k8_b2048_control": 2,
    }
    require(
        {role: sum(cell["role"] == role for cell in cells)
         for role in expected_roles} == expected_roles and
        len({cell["id"] for cell in cells}) == 46 and
        len({cell["seed"] for cell in cells}) == 46 and
        sum(cell["candidate_selected"] for cell in cells) == 44 and
        sum(cell["control_selected"] for cell in cells) == 34,
        "K8/B1024 role or selector accounting changed")

    runner_path = Path(RUNNER.__file__).resolve()
    require(BASE.BENCHMARK_CPU == 14 and BASE.RESERVED_SIBLING == 30,
            "authoritative CPU pair changed")
    require(BASE.MODE_SYMBOL ==
            "_ZN12_GLOBAL__N_1L25g_k8r3r4_t4_terminal_modeE" and
            BASE.AUXILIARY_MODE_EXPECTATIONS == {
                "_ZN12_GLOBAL__N_1L28g_high_t4_batch_binding_modeE": {
                    "candidate": 1,
                    "control": 1,
                },
            },
            "selector evidence policy changed")
    require(BASE.SCHEMA ==
            "leopard2-gf8-t4-k8-b1024-terminal-abba/v1" and
            BASE.SUMMARY_SCHEMA ==
            "leopard2-gf8-t4-k8-b1024-terminal-summary/v1",
            "K8/B1024 evidence schema changed")
    require(BASE.RUNNER_PATH == runner_path and os.access(runner_path, os.X_OK),
            "K8/B1024 runner is not its executable provenance root")
    expected_dependencies = (
        runner_path,
        Path(BASE.__file__).resolve(),
        Path(BASE.T8_SUPPORT.__file__).resolve(),
        Path(BASE.MAIN_SUPPORT.__file__).resolve(),
        Path(BASE.MAIN_SUPPORT.git_capture.__file__).resolve(),
        Path(BASE.MAIN_SUPPORT.git_capture._build_provenance.__file__).resolve(),
        Path(BASE.MAIN_SUPPORT.link_common.__file__).resolve(),
    )
    require(BASE.RUNNER_DEPENDENCIES == expected_dependencies and
            len(set(expected_dependencies)) == 7,
            "K8/B1024 runner dependency chain changed")
    require(BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL is True and
            BASE.CONTROL_EXTRA_ARGUMENTS ==
                ("--disable-k8r3r4-t4-terminal",) and
            BASE.CONTROL_BUILD_MARKER ==
                "k8r3r4_t4_terminal_diagnostic_disabled" and
            BASE.CANDIDATE_SCHEMA == "leopard2-benchmark-v5" and
            BASE.CONTROL_SCHEMA == "leopard2-benchmark-v11" and
            BASE.REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE is False and
            BASE.REQUIRE_EQUAL_EXECUTABLE_PATH_LENGTHS is True and
            BASE.EXPECTED_BINARY_SHA256 == {
                "candidate":
                    "b52e01f87fd2095e4a590ef5df1fc1fb0b8dc303eb65c2a0defa73d49696ad43",
                "control":
                    "b52e01f87fd2095e4a590ef5df1fc1fb0b8dc303eb65c2a0defa73d49696ad43",
                "main":
                    "a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910",
            },
            "fail-closed frozen binary policy changed")
    sample_cell = cells[42]
    candidate_command = BASE.benchmark_command(
        "candidate", Path("/frozen/a/benchmark"), sample_cell, 14, 15, 6)
    control_command = BASE.benchmark_command(
        "control", Path("/frozen/b/benchmark"), sample_cell, 14, 15, 6)
    require("--disable-k8r3r4-t4-terminal" not in candidate_command and
            control_command.count("--disable-k8r3r4-t4-terminal") == 1,
            "runtime control argument attribution changed")
    equal_paths = {
        "candidate": {"path": "/frozen/a/benchmark"},
        "control": {"path": "/frozen/b/benchmark"},
        "main": {"path": "/frozen/m/benchmark"},
    }
    require(BASE.executable_path_lengths(equal_paths) == {
                "candidate": 19, "control": 19, "main": 19},
            "equal executable paths were not accepted")
    unequal_paths = copy.deepcopy(equal_paths)
    unequal_paths["control"]["path"] = "/frozen/control/benchmark"
    expect_failure(lambda: BASE.executable_path_lengths(unequal_paths),
                   "unequal executable paths were accepted")

    analyses = passing_analyses(cells)
    passed = RUNNER.evaluate_promotion(analyses)
    require(passed["passed"] is True and passed["failure_count"] == 0 and
            passed["role_counts"] == expected_roles,
            "valid synthetic campaign failed")

    failure = one_failure(
        analyses, "k8_b1024_target", "candidate_vs_control", "ci95",
        [1.049, 1.20])
    require(failure["contrast"] == "candidate_vs_control" and
            failure["statistic"] == "ci95_lower" and
            failure["required_relation"] == ">=",
            "B1024 same-source target gate changed")
    failure = one_failure(
        analyses, "k8_regression", "candidate_vs_main", "ci95",
        [1.049, 1.20])
    require(failure["contrast"] == "candidate_vs_main",
            "prior K8 exact-main target gate changed")
    failure = one_failure(
        analyses, "preexisting_regression", "candidate_vs_control",
        "speedup", 1.021)
    require(failure["required_relation"] == "<=",
            "K4..7 upper same-path gate changed")
    failure = one_failure(
        analyses, "b512_regression", "candidate_vs_control", "speedup",
        0.979)
    require(failure["required_relation"] == ">=",
            "K4..7 lower same-path gate changed")
    failure = one_failure(
        analyses, "preexisting_regression", "candidate_vs_main", "speedup",
        0.979)
    require(failure["contrast"] == "candidate_vs_main",
            "K4..7 exact-main regression gate changed")
    failure = one_failure(
        analyses, "k8_b2048_control", "candidate_vs_control", "speedup",
        1.021)
    require(failure["required_relation"] == "<=",
            "B2048 same-path gate changed")

    descriptive = copy.deepcopy(analyses)
    b2048 = next(item for item in descriptive
                 if item["cell"]["role"] == "k8_b2048_control")
    b2048["candidate_vs_main"]["speedup"] = 0.01
    b2048["candidate_vs_main"]["ci95"] = [0.005, 0.02]
    require(RUNNER.evaluate_promotion(descriptive)["passed"] is True,
            "B2048 exact-main descriptive result became a gate")
    expect_failure(lambda: RUNNER.evaluate_promotion(analyses[:-1]),
                   "incomplete promotion matrix was accepted")

    print("K8/B1024 T4 ABBA runner self-test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
