#!/usr/bin/env python3
"""Isolated fail-closed checks for the final97 75-round confirmation."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import math
import statistics
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any, Mapping


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_runner() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_current_atlas_final97_confirmation_abba.py")
    specification = importlib.util.spec_from_file_location(
        "current_atlas_final97_confirmation_test_target", path)
    require(specification is not None and specification.loader is not None,
            "cannot load final97 confirmation runner")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()
BASE = RUNNER.BASE


def require_rejected(callback: Any, message: str) -> None:
    try:
        callback()
    except BASE.EvidenceError:
        return
    raise RuntimeError(message)


def synthetic_round(
    order: tuple[str, ...], round_index: int,
    missing: list[int], *, control_primary: float = 1.0,
    control_companion: float = 1.0, main_primary: float = 1.05,
    main_companion: float = 1.05,
) -> dict[str, Any]:
    values = {
        "candidate": (1.0, 1.0),
        "control": (control_primary, control_companion),
        "main": (main_primary, main_companion),
    }
    digests = {
        "algorithm": "fnv1a64",
        "original_data": "1" * 16,
        "transmitted_parity": "2" * 16,
        "recovered_originals": "3" * 16,
    }
    return {
        "attempt": 0,
        "round": round_index,
        "order": list(order),
        "invocations": [{
            "implementation": label,
            "normalized": {
                "encode_us": values[label][0],
                "one_shot_encode_us": values[label][1],
                "digests": digests,
                "missing_original_indices": missing,
            },
        } for label in order],
        "isolation": {"accepted": True},
        "discarded_attempts": [],
    }


def synthetic_campaign() -> tuple[dict[str, Any], dict[str, Any]]:
    orders = RUNNER.select_round_orders(
        BASE.TARGET_ORDER, RUNNER.CONFIRMATORY_ROUNDS)
    raw_cells = []
    for cell in RUNNER.cells():
        missing = list(cell["missing_original_indices"])
        rounds = [
            synthetic_round(order, index, missing)
            for index, order in enumerate(orders)
        ]
        raw_cells.append({"cell": cell, "rounds": rounds})
    identities = {
        "candidate": {"path": "/frozen/shared", "sha256": "a" * 64},
        "control": {"path": "/frozen/shared", "sha256": "a" * 64},
        "main": {"path": "/frozen/main", "sha256": "b" * 64},
    }
    raw = {
        "schema": RUNNER.RAW_SCHEMA,
        "source_commit": "c" * 40,
        "source_tree": "d" * 40,
        "main_commit": BASE.MAIN_COMMIT,
        "requested_rounds": RUNNER.CONFIRMATORY_ROUNDS,
        "cpu": 30,
        "iterations": 31,
        "warmup": 64,
        "identities": identities,
        "cells": raw_cells,
    }
    analyses = [
        RUNNER.SUPPORT.analyze(raw_cell["cell"], raw_cell["rounds"])
        for raw_cell in raw_cells
    ]
    summary = {
        "schema": RUNNER.PRELIMINARY_SUMMARY_SCHEMA,
        "status": "accepted",
        "source_commit": raw["source_commit"],
        "source_tree": raw["source_tree"],
        "main_commit": raw["main_commit"],
        "cell_count": RUNNER.EXPECTED_CELL_COUNT,
        "process_count": RUNNER.EXPECTED_ACCEPTED_PROCESS_COUNT,
        "discarded_process_count": 0,
        "discarded_round_attempts": 0,
        "binary_sha256": {
            name: identity["sha256"]
            for name, identity in identities.items()
        },
        "cells": analyses,
    }
    return raw, summary


def exercise_schemas_and_dependencies() -> None:
    require(
        Path(RUNNER.__file__).stat().st_mode & 0o111 != 0 and
        Path(__file__).stat().st_mode & 0o111 != 0 and
        RUNNER.RAW_SCHEMA.endswith("/v3") and
        RUNNER.PRELIMINARY_SUMMARY_SCHEMA.endswith("/v3") and
        RUNNER.FINAL_SUMMARY_SCHEMA.endswith("/v3") and
        RUNNER.ROUND_PAYLOAD_CLOSURE_SCHEMA.endswith("/v3") and
        len({RUNNER.RAW_SCHEMA, RUNNER.PRELIMINARY_SUMMARY_SCHEMA,
             RUNNER.FINAL_SUMMARY_SCHEMA,
             RUNNER.ROUND_PAYLOAD_CLOSURE_SCHEMA}) == 4 and
        BASE.SCHEMA == RUNNER.RAW_SCHEMA and
        BASE.SUMMARY_SCHEMA == RUNNER.PRELIMINARY_SUMMARY_SCHEMA,
        "confirmation evidence schemas are not distinct v3 contracts")
    expected_dependencies = (
        Path(RUNNER.__file__).resolve(), *RUNNER.SUPPORT_DEPENDENCIES)
    require(
        BASE.RUNNER_PATH == Path(RUNNER.__file__).resolve() and
        BASE.RUNNER_DEPENDENCIES == expected_dependencies and
        len(set(expected_dependencies)) == len(expected_dependencies) and
        BASE.MAX_ISOLATION_ATTEMPTS == 8 and
        Path(RUNNER.SUPPORT.__file__).resolve() in expected_dependencies and
        RUNNER.SUPPORT.ATLAS_MANIFEST_PATH.resolve() in expected_dependencies
        and RUNNER.SUPPORT.ATLAS_SUMMARY_PATH.resolve() in expected_dependencies
        and RUNNER.SUPPORT.ATLAS_GENERATOR_PATH.resolve() in
            expected_dependencies,
        "confirmation omitted part of the inherited provenance closure")
    require(
        RUNNER.REJECTED_PILOT_SUMMARY_SHA256 ==
            "a4df8880ca7c70d892f6b0da0647792daed163e0316612beeb6e9f9bb529c29c"
        and RUNNER.REJECTED_PILOT_RAW_SHA256 ==
            "5823b85c41b2f3cff081dfc1ee776b139a61354427afdfa294703081be3d4d1c",
        "rejected pilot evidence binding changed")


def exercise_fixed_matrix_and_orders() -> None:
    cells = RUNNER.cells()
    require(
        len(cells) == 97 and cells == RUNNER.SUPPORT.cells() and
        RUNNER._canonical_sha256(cells) == RUNNER.MATRIX_SHA256 and
        len({cell["id"] for cell in cells}) == 97 and
        all(cell["compare_main"] is True and
            cell["measure_one_shot"] is True for cell in cells),
        "confirmation does not contain the exact complete final97 matrix")
    orders = RUNNER.select_round_orders(BASE.TARGET_ORDER, 75)
    require(
        len(orders) == 75 and
        all(orders[index] == BASE.TARGET_ORDER[index % 3]
            for index in range(75)) and
        all(orders.count(order) == 25 for order in BASE.TARGET_ORDER),
        "75-round schedule is not 25 complete balanced cycles")
    for position in range(6):
        require(
            {label: sum(order[position] == label for order in orders)
             for label in ("candidate", "control", "main")} ==
                {"candidate": 25, "control": 25, "main": 25},
            "a process position is not balanced across labels")
    require_rejected(
        lambda: RUNNER.select_round_orders(BASE.TARGET_ORDER, 25),
        "the rejected pilot round count was accepted")
    require_rejected(
        lambda: RUNNER.select_round_orders(BASE.TARGET_ORDER[:2], 75),
        "a partial ABBA cycle was accepted")


def exercise_confidence_interval() -> None:
    values = [((index % 9) - 4) / 10_000.0 for index in range(75)]
    observed = RUNNER.confidence_interval(values)
    center = statistics.mean(values)
    half_width = RUNNER.T95_DF74 * statistics.stdev(values) / math.sqrt(75)
    require(
        observed["degrees_of_freedom"] == 74 and
        observed["t_critical"] == 1.992543495180924 and
        observed["confidence_level"] == 0.95 and
        observed["round_log_ratios"] == values and
        math.isclose(observed["speedup"], math.exp(center),
                     rel_tol=0.0, abs_tol=1e-15) and
        all(math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-15)
            for actual, expected in zip(
                observed["ci95"],
                (math.exp(center - half_width),
                 math.exp(center + half_width)))),
        "75-round Student interval changed")
    require_rejected(lambda: RUNNER.confidence_interval(values[:-1]),
                     "74 contrasts were accepted")
    require_rejected(lambda: RUNNER.confidence_interval(values + [0.0]),
                     "76 contrasts were accepted")
    malformed = list(values)
    malformed[0] = float("nan")
    require_rejected(lambda: RUNNER.confidence_interval(malformed),
                     "a non-finite contrast was accepted")


def exercise_argument_and_command_contract() -> None:
    digest = "a" * 64
    arguments = [
        str(BASE.RUNNER_PATH), "--candidate", "/build/bench_leopard2",
        "--candidate-sha256", digest,
        "--control", "/build/bench_leopard2",
        "--control-sha256", digest,
        "--main", "/frozen/main",
        "--main-sha256", BASE.CANONICAL_MAIN_SHA256,
        "--source-commit", "c" * 40, "--source-tree", "d" * 40,
        "--output", "/results/final97-v3", "--cpu", "30",
        "--sibling", "31",
    ]
    saved_argv = sys.argv
    try:
        sys.argv = arguments
        parsed = RUNNER.parse_arguments()
        require(parsed.rounds == 75 and parsed.iterations == 31 and
                parsed.warmup == 64 and
                parsed.candidate_archive is None and
                parsed.candidate_compile_commands is None,
                "authoritative 75-round CLI defaults changed")
        for option, value in (
            ("--rounds", "25"),
            ("--iterations", "30"),
            ("--warmup", "63"),
        ):
            sys.argv = arguments + [option, value]
            with contextlib.redirect_stderr(io.StringIO()):
                try:
                    RUNNER.parse_arguments()
                except SystemExit as error:
                    require(error.code != 0,
                            f"invalid {option} exited successfully")
                else:
                    raise RuntimeError(f"CLI accepted invalid {option}")
    finally:
        sys.argv = saved_argv

    representative = RUNNER.cells()[0]
    candidate = BASE.benchmark_command(
        "candidate", Path("/frozen/shared"), representative, 30, 31, 64)
    control = BASE.benchmark_command(
        "control", Path("/frozen/shared"), representative, 30, 31, 64)
    main = BASE.benchmark_command(
        "main", Path("/frozen/main"), representative, 30, 31, 64)
    require(
        candidate == control and
        "--measure-one-shot-decode" in candidate and
        "--measure-one-shot-decode" not in main and
        BASE.CANONICAL_MAIN_SHA256 ==
            "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93",
        "same-binary route or canonical exact-main identity changed")


def exercise_full_matrix_validation_and_gates() -> None:
    raw, preliminary = synthetic_campaign()
    counts = RUNNER.validate_confirmation_matrix(raw, preliminary)
    require(
        counts["accepted_process_count"] == 43_650 and
        counts["accepted_round_count"] == 7_275 and
        counts["discarded_process_count"] == 0 and
        counts["discarded_attempt_count"] == 0,
        "complete synthetic campaign counts changed")

    first = raw["cells"][0]
    saved_cell = first["cell"]
    changed_cell = deepcopy(saved_cell)
    changed_cell["missing_original_indices"] = []
    first["cell"] = changed_cell
    require_rejected(
        lambda: RUNNER.validate_confirmation_matrix(raw, preliminary),
        "mutated erasure metadata was accepted")
    first["cell"] = saved_cell

    raw["cells"][0], raw["cells"][1] = raw["cells"][1], raw["cells"][0]
    require_rejected(
        lambda: RUNNER.validate_confirmation_matrix(raw, preliminary),
        "reordered cells were accepted")
    raw["cells"][0], raw["cells"][1] = raw["cells"][1], raw["cells"][0]

    removed = first["rounds"].pop()
    require_rejected(
        lambda: RUNNER.validate_confirmation_matrix(raw, preliminary),
        "a 74-round cell was accepted")
    first["rounds"].append(removed)
    first["rounds"].append(removed)
    require_rejected(
        lambda: RUNNER.validate_confirmation_matrix(raw, preliminary),
        "a 76-round cell was accepted")
    first["rounds"].pop()

    saved_order = first["rounds"][0]["order"]
    changed_order = list(saved_order)
    changed_order[0], changed_order[1] = changed_order[1], changed_order[0]
    first["rounds"][0]["order"] = changed_order
    require_rejected(
        lambda: RUNNER.validate_confirmation_matrix(raw, preliminary),
        "an altered ABBA order was accepted")
    first["rounds"][0]["order"] = saved_order

    first["rounds"][0]["isolation"]["accepted"] = False
    require_rejected(
        lambda: RUNNER.validate_confirmation_matrix(raw, preliminary),
        "a contaminated accepted round was accepted")
    first["rounds"][0]["isolation"]["accepted"] = True

    saved_attempt = first["rounds"][0]["attempt"]
    saved_discards = first["rounds"][0]["discarded_attempts"]
    first["rounds"][0]["discarded_attempts"] = [{
        "attempt": index,
        "order": list(first["rounds"][0]["order"]),
        "invocations": first["rounds"][0]["invocations"],
        "isolation": {"accepted": False},
    } for index in range(8)]
    first["rounds"][0]["attempt"] = 8
    require_rejected(
        lambda: RUNNER.validate_confirmation_matrix(raw, preliminary),
        "more than seven discarded isolation attempts were accepted")
    first["rounds"][0]["discarded_attempts"] = saved_discards
    first["rounds"][0]["attempt"] = saved_attempt

    preliminary["cells"][0]["control_over_candidate"]["speedup"] = 9.0
    require_rejected(
        lambda: RUNNER.validate_confirmation_matrix(raw, preliminary),
        "a rewritten preliminary statistic was accepted")
    preliminary["cells"][0]["control_over_candidate"]["speedup"] = 1.0
    RUNNER.validate_confirmation_matrix(raw, preliminary)

    # A structurally complete compact campaign cannot pass closure without the
    # full journal; this catches accidental fail-open artifact handling without
    # creating 7,275 fixture files.
    require_rejected(
        lambda: RUNNER.validate_round_payload_artifacts(
            raw, preliminary, Path("/nonexistent/final97-v3")),
        "missing full round payload artifacts were accepted")

    accepted = RUNNER.apply_final_acceptance(deepcopy(preliminary))
    contract = accepted["acceptance_contract"]
    require(
        accepted["status"] == "accepted" and
        accepted["schema"] == RUNNER.FINAL_SUMMARY_SCHEMA and
        not accepted["exact_main_regressions"] and
        not accepted["same_binary_equivalence_failures"] and
        contract["accepted_rounds_per_cell"] == 75 and
        contract["balanced_order_cycles"] == 25 and
        contract["iterations_per_process"] == 31 and
        contract["warmup_iterations_per_process"] == 64 and
        contract["maximum_isolation_attempts_per_round"] == 8 and
        contract["discarded_isolation_attempts_used_in_inference"] is False and
        contract["pooling"] == contract["trimming"] == "none" and
        contract["outcome_selected_cells"] is False and
        contract["statistical_campaign_retries"] == 0 and
        contract["rejected_v2_pilot"]["used_for_acceptance"] is False,
        "valid standalone confirmation was rejected or mixed with the pilot")

    def rejected_gate(
        ratio_name: str, interval: list[float], failure_name: str,
    ) -> None:
        changed = deepcopy(preliminary)
        changed["cells"][0][ratio_name]["ci95"] = interval
        RUNNER.apply_final_acceptance(changed)
        require(changed["status"] == "rejected" and changed[failure_name],
                f"gate accepted {ratio_name} interval {interval}")

    rejected_gate(
        "control_over_candidate", [0.97, 1.0],
        "same_binary_equivalence_failures")
    rejected_gate(
        "one_shot_control_over_candidate", [1.0, 1.021],
        "same_binary_equivalence_failures")
    rejected_gate(
        "main_over_candidate", [0.999, 1.01],
        "exact_main_regressions")
    rejected_gate(
        "one_shot_main_over_candidate", [0.999, 1.01],
        "exact_main_regressions")

    substituted = deepcopy(preliminary)
    substituted["cells"][-1] = deepcopy(substituted["cells"][0])
    require_rejected(lambda: RUNNER.apply_final_acceptance(substituted),
                     "a selected/substituted final matrix was accepted")


def main() -> int:
    exercise_schemas_and_dependencies()
    exercise_fixed_matrix_and_orders()
    exercise_confidence_interval()
    exercise_argument_and_command_contract()
    exercise_full_matrix_validation_and_gates()
    print("Final97 standalone 75-round confirmation checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
