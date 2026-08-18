#!/usr/bin/env python3
"""Standalone 75-round confirmation of the complete final97 matrix.

This wrapper deliberately reuses the hardened final97 collector, immutable
binary handling, command validation, round journal, and build-provenance
closure.  It changes only the predeclared statistical design: all 97 cells are
measured for 75 accepted rounds, comprising 25 complete repetitions of the
three balanced six-process ABBA orders.  The rejected 25-round campaign is a
pilot and contributes no observations to this confirmation.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import math
import os
import stat
import statistics
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any, Mapping, Sequence


SUPPORT_PATH = Path(__file__).resolve().with_name(
    "run_current_atlas_regression_screen_abba.py")


def load_support() -> Any:
    specification = importlib.util.spec_from_file_location(
        "current_atlas_final97_confirmation_support", SUPPORT_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load final97 support: {SUPPORT_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()
BASE = SUPPORT.BASE
SUPPORT_DEPENDENCIES = tuple(BASE.RUNNER_DEPENDENCIES)

RAW_SCHEMA = "leopard2-current-atlas-final97-confirmation-abba/v3"
PRELIMINARY_SUMMARY_SCHEMA = \
    "leopard2-current-atlas-final97-confirmation-preliminary-summary/v3"
FINAL_SUMMARY_SCHEMA = \
    "leopard2-current-atlas-final97-confirmation-summary/v3"
ROUND_PAYLOAD_CLOSURE_SCHEMA = \
    "leopard2-current-atlas-final97-confirmation-round-payload-closure/v3"
ROUND_PAYLOAD_SCHEMA = "leopard2-current-atlas-round-payload/v1"
CONFIRMATORY_ROUNDS = 75
BALANCED_ORDER_CYCLES = 25
CONFIDENCE_LEVEL = 0.95
T95_DF74 = 1.992543495180924
EXPECTED_CELL_COUNT = 97
EXPECTED_ACCEPTED_PROCESS_COUNT = 43_650
EXPECTED_ACCEPTED_ARTIFACT_COUNT = 7_275
MATRIX_SHA256 = \
    "fe21e525ff084564af37a4e7b43f17df8e971f96a252199125e26d50e445da8b"
REJECTED_PILOT_SUMMARY_SHA256 = \
    "a4df8880ca7c70d892f6b0da0647792daed163e0316612beeb6e9f9bb529c29c"
REJECTED_PILOT_RAW_SHA256 = \
    "5823b85c41b2f3cff081dfc1ee776b139a61354427afdfa294703081be3d4d1c"

BASE.__doc__ = __doc__
SUPPORT.__doc__ = __doc__
BASE.SCHEMA = RAW_SCHEMA
BASE.SUMMARY_SCHEMA = PRELIMINARY_SUMMARY_SCHEMA
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    *SUPPORT_DEPENDENCIES,
)
BASE.require(
    len(set(BASE.RUNNER_DEPENDENCIES)) == len(BASE.RUNNER_DEPENDENCIES) and
    BASE.MAX_ISOLATION_ATTEMPTS == 8,
    "confirmation dependency closure or isolation-attempt cap changed")


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"),
        ensure_ascii=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def cells() -> list[dict[str, Any]]:
    """Return an independently copied, hash-pinned complete final97 matrix."""
    result = deepcopy(SUPPORT.cells())
    BASE.require(
        len(result) == EXPECTED_CELL_COUNT and
        len({cell.get("id") for cell in result}) == EXPECTED_CELL_COUNT and
        all(cell.get("compare_main") is True and
            cell.get("measure_one_shot") is True for cell in result) and
        _canonical_sha256(result) == MATRIX_SHA256,
        "confirmation matrix differs from the predeclared full 97 cells")
    return result


BASE.cells = cells


def confidence_interval(values: Sequence[float]) -> dict[str, Any]:
    """Two-sided Student interval over exactly 75 round-level contrasts."""
    BASE.require(
        len(values) == CONFIRMATORY_ROUNDS and
        all(isinstance(value, (int, float)) and
            not isinstance(value, bool) and math.isfinite(float(value))
            for value in values),
        "exactly 75 finite numeric round contrasts are required")
    numeric = [float(value) for value in values]
    center = statistics.mean(numeric)
    half_width = T95_DF74 * statistics.stdev(numeric) / math.sqrt(
        CONFIRMATORY_ROUNDS)
    return {
        "speedup": math.exp(center),
        "ci95": [
            math.exp(center - half_width),
            math.exp(center + half_width),
        ],
        "round_log_ratios": numeric,
        "confidence_level": CONFIDENCE_LEVEL,
        "degrees_of_freedom": CONFIRMATORY_ROUNDS - 1,
        "t_critical": T95_DF74,
    }


BASE.confidence_interval = confidence_interval


def select_round_orders(
    orders: Sequence[Sequence[str]], requested_rounds: int | None,
) -> tuple[tuple[str, ...], ...]:
    """Return exactly 25 repetitions of the three balanced ABBA orders."""
    normalized = tuple(tuple(order) for order in orders)
    BASE.require(
        normalized == tuple(BASE.TARGET_ORDER) and
        len(normalized) == 3 and requested_rounds == CONFIRMATORY_ROUNDS,
        "confirmation requires its exact three-order, 75-round schedule")
    return tuple(normalized[index % 3]
                 for index in range(CONFIRMATORY_ROUNDS))


BASE.select_round_orders = select_round_orders


def parse_arguments() -> argparse.Namespace:
    """Expose one fixed confirmation design and no statistical retry knob."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--candidate-sha256")
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--control-sha256")
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument("--main-sha256", default=BASE.CANONICAL_MAIN_SHA256)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--iterations", type=int, choices=(31,), default=31)
    parser.add_argument("--warmup", type=int, choices=(64,), default=64)
    parser.add_argument(
        "--rounds", type=int, choices=(CONFIRMATORY_ROUNDS,),
        default=CONFIRMATORY_ROUNDS)
    parser.set_defaults(
        candidate_archive=None,
        candidate_archive_sha256=None,
        control_archive=None,
        control_archive_sha256=None,
        candidate_compile_commands=None,
        candidate_compile_commands_sha256=None,
        control_compile_commands=None,
        control_compile_commands_sha256=None,
    )
    options = parser.parse_args()
    BASE.require(
        options.main_sha256 == BASE.CANONICAL_MAIN_SHA256,
        "exact-main executable SHA-256 differs from the canonical AVX2 "
        "Leopard-main benchmark")
    return options


BASE.parse_arguments = parse_arguments


def _read_json(path: Path) -> Mapping[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as source:
            value = json.load(source)
    except (OSError, json.JSONDecodeError) as error:
        raise BASE.EvidenceError(f"cannot read evidence JSON {path}: {error}") \
            from error
    BASE.require(isinstance(value, Mapping),
                 f"evidence JSON is not an object: {path}")
    return value


def _ratio_names() -> tuple[str, ...]:
    return (
        "control_over_candidate", "one_shot_control_over_candidate",
        "main_over_candidate", "one_shot_main_over_candidate",
    )


def validate_confirmation_matrix(
    raw: Mapping[str, Any], summary: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate the exact 97-by-75 design and recompute every analysis."""
    raw_cells = raw.get("cells")
    expected_cells = cells()
    expected_orders = select_round_orders(
        BASE.TARGET_ORDER, CONFIRMATORY_ROUNDS)
    BASE.require(
        raw.get("schema") == RAW_SCHEMA and
        raw.get("requested_rounds") == CONFIRMATORY_ROUNDS and
        raw.get("iterations") == 31 and raw.get("warmup") == 64 and
        isinstance(raw_cells, list) and
        len(raw_cells) == EXPECTED_CELL_COUNT and
        summary.get("schema") == PRELIMINARY_SUMMARY_SCHEMA and
        summary.get("cell_count") == EXPECTED_CELL_COUNT and
        summary.get("source_commit") == raw.get("source_commit") and
        summary.get("source_tree") == raw.get("source_tree") and
        summary.get("main_commit") == raw.get("main_commit"),
        "confirmation raw or preliminary summary identity is malformed")

    accepted_processes = 0
    discarded_processes = 0
    discarded_attempt_count = 0
    for cell_index, raw_cell in enumerate(raw_cells):
        cell = raw_cell.get("cell")
        rounds = raw_cell.get("rounds")
        BASE.require(
            cell == expected_cells[cell_index] and
            isinstance(rounds, list) and
            len(rounds) == CONFIRMATORY_ROUNDS,
            "confirmation cell matrix or round count changed")
        position_counts = [
            {label: 0 for label in ("candidate", "control", "main")}
            for _ in range(6)
        ]
        label_counts = {
            label: 0 for label in ("candidate", "control", "main")}
        for round_index, accepted_round in enumerate(rounds):
            expected_order = list(expected_orders[round_index])
            discarded = accepted_round.get("discarded_attempts")
            BASE.require(
                accepted_round.get("round") == round_index and
                accepted_round.get("order") == expected_order and
                isinstance(discarded, list) and
                len(discarded) < BASE.MAX_ISOLATION_ATTEMPTS,
                "accepted round order, index, or discard list changed")
            attempts = [(True, attempt) for attempt in discarded]
            attempts.append((False, accepted_round))
            for discarded_index, (is_discarded, attempt) in enumerate(
                    attempts):
                invocations = attempt.get("invocations")
                isolation = attempt.get("isolation")
                expected_attempt_index = discarded_index
                BASE.require(
                    attempt.get("attempt") == expected_attempt_index and
                    attempt.get("order") == expected_order and
                    isinstance(invocations, list) and
                    len(invocations) == 6 and
                    [item.get("implementation")
                     for item in invocations] == expected_order and
                    isinstance(isolation, Mapping) and
                    isolation.get("accepted") is (not is_discarded),
                    "confirmation round attempt is malformed")
                if is_discarded:
                    discarded_processes += 6
                    discarded_attempt_count += 1
                else:
                    accepted_processes += 6
                    for position, label in enumerate(expected_order):
                        position_counts[position][label] += 1
                        label_counts[label] += 1
        BASE.require(
            label_counts == {"candidate": 150, "control": 150, "main": 150}
            and all(counts == {"candidate": 25, "control": 25, "main": 25}
                    for counts in position_counts),
            "confirmation order cycle is not label/position balanced")

    BASE.require(
        accepted_processes == EXPECTED_ACCEPTED_PROCESS_COUNT and
        summary.get("process_count") == accepted_processes and
        summary.get("discarded_process_count") == discarded_processes and
        summary.get("discarded_round_attempts") == discarded_attempt_count and
        discarded_processes == 6 * discarded_attempt_count,
        "confirmation accepted or discarded process count changed")
    recomputed = [
        SUPPORT.analyze(raw_cell["cell"], raw_cell["rounds"])
        for raw_cell in raw_cells
    ]
    BASE.require(summary.get("cells") == recomputed,
                 "preliminary analyses do not recompute from all 75 rounds")
    for analysis in recomputed:
        for name in _ratio_names():
            ratio = analysis.get(name)
            BASE.require(
                isinstance(ratio, Mapping) and
                ratio.get("degrees_of_freedom") == 74 and
                ratio.get("confidence_level") == CONFIDENCE_LEVEL and
                ratio.get("t_critical") == T95_DF74 and
                isinstance(ratio.get("round_log_ratios"), list) and
                len(ratio["round_log_ratios"]) == CONFIRMATORY_ROUNDS,
                "analysis is not the predeclared 75-round Student interval")
    identities = raw.get("identities")
    BASE.require(
        isinstance(identities, Mapping) and
        set(identities) == {"candidate", "control", "main"} and
        summary.get("binary_sha256") == {
            name: identities[name].get("sha256")
            for name in ("candidate", "control", "main")
        },
        "summary binary identities do not match compact raw evidence")
    return {
        "accepted_process_count": accepted_processes,
        "discarded_process_count": discarded_processes,
        "discarded_attempt_count": discarded_attempt_count,
        "accepted_round_count": EXPECTED_ACCEPTED_ARTIFACT_COUNT,
        "recomputed_analyses": recomputed,
    }


def validate_round_payload_artifacts(
    raw: Mapping[str, Any], summary: Mapping[str, Any], output: Path,
) -> dict[str, Any]:
    """Reopen and semantically validate every accepted/discarded artifact."""
    matrix = validate_confirmation_matrix(raw, summary)
    raw_cells = raw["cells"]
    expected_orders = select_round_orders(
        BASE.TARGET_ORDER, CONFIRMATORY_ROUNDS)
    identities = raw["identities"]
    expected_sequence = 0
    accepted_artifacts = 0
    discarded_artifacts = 0
    artifact_hashes: list[str] = []
    seen_paths: set[str] = set()
    for cell_index, raw_cell in enumerate(raw_cells):
        cell = raw_cell["cell"]
        for round_index, accepted_round in enumerate(raw_cell["rounds"]):
            expected_order = list(expected_orders[round_index])
            attempts = [
                (True, attempt)
                for attempt in accepted_round["discarded_attempts"]
            ]
            attempts.append((False, accepted_round))
            for is_discarded, attempt in attempts:
                invocations = attempt["invocations"]
                isolation = attempt["isolation"]
                references = [
                    invocation.get("round_payload_artifact")
                    for invocation in invocations
                ]
                BASE.require(
                    all(isinstance(reference, Mapping)
                        for reference in references),
                    "compact invocation lacks its full round payload")
                artifact_identity = {
                    name: references[0].get(name)
                    for name in ("path", "size", "sha256")
                }
                expected_path = (
                    output / "round_payloads" /
                    f"round-{expected_sequence:06d}.json").resolve()
                BASE.require(
                    all({name: reference.get(name) for name in
                         ("path", "size", "sha256")} == artifact_identity
                        for reference in references) and
                    [reference.get("invocation_index")
                     for reference in references] == list(range(6)) and
                    artifact_identity["path"] == str(expected_path) and
                    artifact_identity["path"] not in seen_paths,
                    "round payload references are duplicated or out of order")
                BASE.require(
                    BASE.support_file_identity(expected_path) ==
                        artifact_identity,
                    "round payload identity changed")
                payload = _read_json(expected_path)
                full_invocations = payload.get("invocations")
                BASE.require(
                    payload.get("schema") == ROUND_PAYLOAD_SCHEMA and
                    payload.get("sequence") == expected_sequence and
                    payload.get("isolation") == isolation and
                    isinstance(full_invocations, list) and
                    len(full_invocations) == 6 and
                    [item.get("implementation")
                     for item in full_invocations] == expected_order,
                    "full round payload differs from compact evidence")
                for compact, full in zip(invocations, full_invocations):
                    BASE.require(
                        isinstance(full.get("result"), Mapping) and
                        isinstance(full.get("command"), list),
                        "full child result or command is absent")
                    implementation = full.get("implementation")
                    identity = identities.get(implementation)
                    BASE.require(
                        implementation in {"candidate", "control", "main"}
                        and isinstance(identity, Mapping),
                        "full child implementation identity is absent")
                    expected_command = SUPPORT.benchmark_command(
                        implementation, Path(str(identity["path"])), cell,
                        raw["cpu"], raw["iterations"], raw["warmup"])
                    expected_normalized = SUPPORT.validate_result(
                        implementation, full["result"], cell,
                        raw["source_commit"], raw["source_tree"],
                        raw["iterations"], raw["warmup"])
                    BASE.require(
                        full["command"] == expected_command and
                        full.get("normalized") == expected_normalized,
                        "full child command or result is not reproducible")
                    expected_compact = {
                        name: value for name, value in full.items()
                        if name not in {"result", "command"}
                    }
                    compact_normalized = dict(
                        expected_compact["normalized"])
                    compact_normalized.pop("encode_samples_us", None)
                    compact_normalized.pop(
                        "one_shot_encode_samples_us", None)
                    expected_compact["normalized"] = compact_normalized
                    observed_compact = {
                        name: value for name, value in compact.items()
                        if name != "round_payload_artifact"
                    }
                    BASE.require(observed_compact == expected_compact,
                                 "compact invocation differs from full payload")
                BASE.require(
                    BASE.support_file_identity(expected_path) ==
                        artifact_identity,
                    "round payload changed during semantic validation")
                seen_paths.add(str(expected_path))
                artifact_hashes.append(str(artifact_identity["sha256"]))
                if is_discarded:
                    discarded_artifacts += 1
                else:
                    accepted_artifacts += 1
                expected_sequence += 1

    payload_directory = output / "round_payloads"
    try:
        directory_entries = list(os.scandir(payload_directory))
    except OSError as error:
        raise BASE.EvidenceError(
            f"cannot enumerate round payload directory: {error}") from error
    actual_paths: set[str] = set()
    for entry in directory_entries:
        try:
            mode = entry.stat(follow_symlinks=False).st_mode
        except OSError as error:
            raise BASE.EvidenceError(
                f"cannot stat round payload entry {entry.path}: {error}") \
                from error
        BASE.require(stat.S_ISREG(mode) and not entry.is_symlink(),
                     "round payload directory contains a non-regular entry")
        actual_paths.add(str(Path(entry.path).resolve()))
    BASE.require(
        accepted_artifacts == EXPECTED_ACCEPTED_ARTIFACT_COUNT and
        discarded_artifacts == matrix["discarded_attempt_count"] and
        expected_sequence == accepted_artifacts + discarded_artifacts and
        len(directory_entries) == expected_sequence and
        len(actual_paths) == expected_sequence and actual_paths == seen_paths,
        "round payload directory is not the exact campaign journal")
    ordered_hash = hashlib.sha256(
        ("\n".join(artifact_hashes) + "\n").encode("ascii")).hexdigest()
    return {
        "schema": ROUND_PAYLOAD_CLOSURE_SCHEMA,
        "matrix_sha256": MATRIX_SHA256,
        "rounds_per_cell": CONFIRMATORY_ROUNDS,
        "balanced_order_cycles": BALANCED_ORDER_CYCLES,
        "artifact_count": expected_sequence,
        "accepted_artifact_count": accepted_artifacts,
        "discarded_artifact_count": discarded_artifacts,
        "accepted_process_count": matrix["accepted_process_count"],
        "discarded_process_count": matrix["discarded_process_count"],
        "ordered_artifact_sha256": ordered_hash,
    }


def apply_final_acceptance(summary: dict[str, Any]) -> dict[str, Any]:
    """Gate all 97 cells and both operations using this campaign alone."""
    SUPPORT.apply_final_acceptance(summary)
    contract = summary["acceptance_contract"]
    contract.update({
        "accepted_rounds_per_cell": CONFIRMATORY_ROUNDS,
        "balanced_order_cycles": BALANCED_ORDER_CYCLES,
        "inference_unit": "one predeclared six-process ABBA round contrast",
        "student_t_degrees_of_freedom": 74,
        "student_t_two_sided_critical": T95_DF74,
        "iterations_per_process": 31,
        "warmup_iterations_per_process": 64,
        "maximum_isolation_attempts_per_round": 8,
        "isolation_attempt_policy":
            "discard contaminated attempts; retain exactly one accepted round",
        "discarded_isolation_attempts_used_in_inference": False,
        "pooling": "none",
        "trimming": "none",
        "matrix_selection": "all_97_predeclared_final97_cells",
        "matrix_sha256": MATRIX_SHA256,
        "outcome_selected_cells": False,
        "statistical_campaign_retries": 0,
        "decision_evidence": "this standalone 75-round campaign only",
        "rejected_v2_pilot": {
            "summary_sha256": REJECTED_PILOT_SUMMARY_SHA256,
            "raw_sha256": REJECTED_PILOT_RAW_SHA256,
            "used_for_acceptance": False,
        },
    })
    summary["schema"] = FINAL_SUMMARY_SCHEMA
    return summary


def main() -> int:
    """Collect one full confirmation and then close all evidence under lock."""
    options = BASE.parse_arguments()
    captured_stdout = io.StringIO()
    with contextlib.redirect_stdout(captured_stdout):
        preliminary_status = BASE.main()
    try:
        finalization_lock, finalization_cpu = \
            SUPPORT._begin_isolated_finalization(options.cpu, options.sibling)
    except Exception as error:
        SUPPORT._bind_failure_journal(options.output)
        print(f"evidence rejected: cannot isolate finalization: {error}",
              file=sys.stderr)
        return 1
    try:
        summary_path = options.output / "summary.json"
        if preliminary_status not in (0, 2):
            SUPPORT._bind_failure_journal(options.output)
            output = captured_stdout.getvalue()
            if output:
                print(output, end="")
            return preliminary_status
        if not summary_path.exists():
            SUPPORT._bind_failure_journal(options.output)
            print("evidence rejected: preliminary summary is absent",
                  file=sys.stderr)
            return 1
        summary = dict(_read_json(summary_path))
        summary["preliminary_generic_status"] = summary.get("status")
        summary["finalization_cpu"] = finalization_cpu
        finalization_failed = False
        try:
            raw_path = options.output / "raw.json"
            observed_raw_sha256 = BASE.sha256(raw_path)
            BASE.require(
                summary.get("raw_sha256") == observed_raw_sha256,
                "preliminary summary does not bind compact raw evidence")
            raw = _read_json(raw_path)
            BASE.require(BASE.sha256(raw_path) == observed_raw_sha256,
                         "compact raw evidence changed while being read")
            SUPPORT.validate_final_live_identities(raw, options)
            summary["round_payload_closure"] = \
                validate_round_payload_artifacts(raw, summary, options.output)
            SUPPORT.validate_final_live_identities(raw, options)
            apply_final_acceptance(summary)
        except Exception as error:
            finalization_failed = True
            summary["schema"] = FINAL_SUMMARY_SCHEMA
            summary["status"] = "rejected"
            summary["finalization_failure"] = {
                "type": type(error).__name__, "message": str(error),
            }
        SUPPORT._replace_json(summary_path, summary)
        print(json.dumps({
            "status": summary.get("status"),
            "cells": summary.get("cell_count"),
            "processes": summary.get("process_count"),
            "discarded_processes": summary.get("discarded_process_count"),
            "exact_main_regressions": len(
                summary.get("exact_main_regressions", [])),
            "same_binary_equivalence_failures": len(
                summary.get("same_binary_equivalence_failures", [])),
        }, sort_keys=True))
        if finalization_failed:
            return 1
        return 0 if summary.get("status") == "accepted" else 2
    finally:
        os.close(finalization_lock)


if __name__ == "__main__":
    raise SystemExit(main())
