#!/usr/bin/env python3
"""Adversarial replay checks for the K65/R65/B64 v2 evidence contract."""

from __future__ import annotations

import argparse
import contextlib
import copy
import importlib.util
import io
import os
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path
from typing import Any, Mapping


def load_runner() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_k65r65_b64_packed_terminal_abba.py")
    specification = importlib.util.spec_from_file_location(
        "leopard2_k65r65_b64_v2_replay_target", path)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load K65/R65/B64 v2 runner")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()
BASE = RUNNER.BASE
SEALED_AUTHORITY = Path(
    "/home/catid/leopard-exact-main-authority-6d4f690-a1/Aauthority")


def authority_binding() -> dict[str, Any]:
    return RUNNER.validate_exact_main_authority_binding({
        "schema": RUNNER.AUTHORITY_BINDING_SCHEMA,
        "lane_root": "/frozen/exact-main-authority",
        "authority_schema": RUNNER.EXACT_MAIN_AUTHORITY_SCHEMA,
        "record_sha256": RUNNER.EXACT_MAIN_AUTHORITY_RECORD_SHA256,
        "ledger_sha256": RUNNER.EXACT_MAIN_AUTHORITY_LEDGER_SHA256,
        "verifier_schema": RUNNER.AUTHORITY_VERIFIER.VERIFIER_SCHEMA,
        "verifier_verdict_sha256": "7" * 64,
        "executable_sha256": RUNNER.EXACT_MAIN_EXECUTABLE_SHA256,
        "archive_sha256": RUNNER.EXACT_MAIN_ARCHIVE_SHA256,
        "combined_sha256": RUNNER.EXACT_MAIN_COMBINED_SHA256,
        "source_commit": RUNNER.EXACT_MAIN_SOURCE_COMMIT,
        "source_tree": RUNNER.EXACT_MAIN_SOURCE_TREE,
        "adapter_commit": RUNNER.EXACT_MAIN_ADAPTER_COMMIT,
        "pure_avx2": True,
        "historical_non_authority": list(RUNNER.HISTORICAL_NON_AUTHORITY),
    })


def terminal(attempt: int, outcome: str = "rejected") -> dict[str, Any]:
    completed = outcome in {"accepted", "rejected"}
    return {
        "schema": RUNNER.ATTEMPT_TERMINAL_SCHEMA,
        "generation": RUNNER.GENERATION,
        "attempt": attempt,
        "attempt_budget": RUNNER.ATTEMPT_BUDGET,
        "outcome": outcome,
        "promotable": outcome == "accepted",
        "output_root": f"/evidence/k65-v2-a{attempt}",
        "summary_schema": RUNNER.FINAL_SUMMARY_SCHEMA_V2
            if completed else "absent",
        "summary_sha256": str(attempt) * 64 if completed else None,
        "raw_sha256": chr(ord("a") + attempt) * 64 if completed else None,
        "failure_sha256": None if completed else "f" * 64,
        "authority_record_sha256":
            RUNNER.EXACT_MAIN_AUTHORITY_RECORD_SHA256,
    }


def invocation(label: str, value: float) -> dict[str, Any]:
    return {
        "implementation": label,
        "normalized": {
            "encode_us": value,
            "one_shot_encode_us": value,
            "digests": {
                "original_data": "1" * 16,
                "transmitted_parity": "2" * 16,
                "recovered_originals": "1" * 16,
            },
        },
    }


def round_record(
    index: int, order: tuple[str, ...], *, control: float, main: float,
) -> dict[str, Any]:
    values = {"candidate": 1.0, "control": control, "main": main}
    return {
        "round": index,
        "attempt": 0,
        "order": list(order),
        "invocations": [invocation(label, values[label]) for label in order],
        "isolation": {
            "accepted": True,
            "delta": {"reserved_sibling": {"nonidle_jiffies": 0}},
        },
        "discarded_attempts": [],
    }


def complete_fixture(
    *, target_control: float = 1.06, target_main: float = 1.06,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    authority = authority_binding()
    preregistration = RUNNER.validate_attempt_preregistration(
        1, [], verify_files=False)
    raw_cells = []
    analyses = []
    process_count = 0
    for cell in RUNNER.cells():
        cycle = BASE.TARGET_ORDER if cell["compare_main"] else BASE.NEIGHBOR_ORDER
        orders = BASE.select_round_orders(cycle, RUNNER.CONFIRMATORY_ROUNDS)
        if cell["role"] == "target":
            control, main = target_control, target_main
        elif cell["main_floor"] is not None:
            control, main = 1.0, 0.99
        elif cell["compare_main"]:
            control, main = 1.0, 0.50
        else:
            control, main = 1.0, 1.0
        rounds = [
            round_record(index, order, control=control, main=main)
            for index, order in enumerate(orders)
        ]
        raw_cells.append({"cell": copy.deepcopy(cell), "rounds": rounds})
        analyses.append(BASE.analyze(cell, rounds))
        process_count += sum(len(item["invocations"]) for item in rounds)
    raw = {
        "schema": RUNNER.RAW_SCHEMA_V2,
        "requested_rounds": RUNNER.CONFIRMATORY_ROUNDS,
        "iterations": RUNNER.FIXED_ITERATIONS,
        "warmup": RUNNER.FIXED_WARMUP,
        "cells": raw_cells,
        "exact_main_authority": copy.deepcopy(authority),
        "attempt_preregistration": copy.deepcopy(preregistration),
    }
    preliminary = {
        "schema": RUNNER.PRELIMINARY_SUMMARY_SCHEMA_V2,
        "status": "accepted",
        "source_commit": "e" * 40,
        "source_tree": "f" * 40,
        "main_commit": RUNNER.EXACT_MAIN_SOURCE_COMMIT,
        "generation": RUNNER.GENERATION,
        "cell_count": len(analyses),
        "process_count": process_count,
        "discarded_process_count": 0,
        "discarded_round_attempts": 0,
        "all_digests_matched": True,
        "all_rounds_zero_sibling_nonidle": True,
        "binary_sha256": {
            "candidate": "c" * 64,
            "control": "c" * 64,
            "main": RUNNER.EXACT_MAIN_EXECUTABLE_SHA256,
        },
        "candidate_control_executable_sections_sha256": "8" * 64,
        "mode_words": {"shared_binary_default": {"value": 1}},
        "cells": analyses,
        "exact_main_authority": copy.deepcopy(authority),
        "attempt_preregistration": copy.deepcopy(preregistration),
        "raw_sha256": "d" * 64,
        "preliminary_generic_status": "accepted",
        "finalization_cpu": 1,
        "preliminary_summary_sha256": "9" * 64,
    }
    return raw, preliminary, authority, preregistration


def main_result(cell: Mapping[str, Any]) -> dict[str, Any]:
    samples = [1.0] * RUNNER.FIXED_ITERATIONS
    metric = {
        "median_us_per_batch_call": 1.0,
        "mad_us_per_batch_call": 0.0,
        "minimum_us_per_batch_call": 1.0,
        "maximum_us_per_batch_call": 1.0,
        "samples_us_per_batch_call": samples,
    }
    return {
        "schema": "leopard-main-benchmark-v1",
        "parameters": {
            "K": cell["K"], "R": cell["R"],
            "shard_bytes": cell["bytes"],
            "logical_shard_bytes": cell["bytes"],
            "loss_count": cell["loss"], "batch": cell["batch"],
            "reuse": cell["reuse"],
            "iterations": RUNNER.FIXED_ITERATIONS,
            "warmup": RUNNER.FIXED_WARMUP,
            "thread_count": 1, "seed": cell["seed"],
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "thread_count": 1, "padded_application_bytes": False,
            "padding_policy": "zero suffix per shard",
        },
        "correctness": {
            "round_trip": True, "logical_prefix_fingerprinted": True,
        },
        "workload_digests": {
            "algorithm": "fnv1a64", "original_data": "1" * 16,
            "transmitted_parity": "2" * 16,
            "recovered_originals": "1" * 16,
        },
        "build": {
            "main_source_commit": RUNNER.EXACT_MAIN_SOURCE_COMMIT,
            "pure_avx2": True, "cplusplus": 201103,
        },
        "metrics": {"encode_execution": metric},
    }


def mocked_main_artifacts(
    *, promotion_error: BaseException | None = None,
    finalization_entry_error: BaseException | None = None,
) -> tuple[int, dict[str, dict[str, Any]], str, str]:
    with tempfile.TemporaryDirectory(
            prefix="leopard-k65-v2-finalizer-") as directory:
        output = Path(directory) / "attempt-1"
        raw, preliminary, authority, preregistration = complete_fixture()
        for name in (
            "preliminary_generic_status", "finalization_cpu",
            "preliminary_summary_sha256"):
            preliminary.pop(name)
        options = argparse.Namespace(
            output=output, cpu=7, sibling=3, attempt=1,
            exact_main_authority=copy.deepcopy(authority),
            exact_main_authority_lane=Path(authority["lane_root"]),
            attempt_preregistration=copy.deepcopy(preregistration),
        )

        def collect() -> int:
            output.mkdir()
            BASE.T8_SUPPORT.write_exclusive(output / "raw.json", raw)
            preliminary["raw_sha256"] = RUNNER._sha256_file(
                output / "raw.json")
            BASE.T8_SUPPORT.write_exclusive(
                output / "summary.json", preliminary)
            return 0

        captured_stdout = io.StringIO()
        captured_stderr = io.StringIO()
        descriptor = os.open(os.devnull, os.O_RDONLY) \
            if finalization_entry_error is None else None
        with contextlib.ExitStack() as stack:
            stack.enter_context(mock.patch.object(
                BASE, "parse_arguments", return_value=options))
            stack.enter_context(mock.patch.object(
                BASE, "main", side_effect=collect))
            if finalization_entry_error is None:
                stack.enter_context(mock.patch.object(
                    RUNNER, "_begin_isolated_finalization",
                    return_value=(descriptor, 11)))
            else:
                stack.enter_context(mock.patch.object(
                    RUNNER, "_begin_isolated_finalization",
                    side_effect=finalization_entry_error))
            stack.enter_context(mock.patch.object(
                RUNNER, "validate_final_live_identities", return_value=None))
            stack.enter_context(mock.patch.object(
                RUNNER, "bind_exact_main_authority",
                return_value=copy.deepcopy(authority)))
            if promotion_error is not None:
                stack.enter_context(mock.patch.object(
                    RUNNER, "validate_promotable_final_summary",
                    side_effect=promotion_error))
            stack.enter_context(contextlib.redirect_stdout(captured_stdout))
            stack.enter_context(contextlib.redirect_stderr(captured_stderr))
            status = RUNNER.main()
        artifacts = {
            name: RUNNER._read_json_object(output / name, f"test {name}")
            for name in (
                "summary.json", "preliminary-summary.json", "failure.json",
                "attempt-terminal.json")
            if (output / name).is_file()
        }
        return (
            status, artifacts, captured_stdout.getvalue(),
            captured_stderr.getvalue())


class K65V2ReplayTests(unittest.TestCase):
    def assert_rejected(self, callback: Any) -> None:
        with self.assertRaises(BASE.EvidenceError):
            callback()

    def test_v1_projection_and_matrix_are_frozen(self) -> None:
        self.assertEqual(
            RUNNER._canonical_sha256(RUNNER.cells()), RUNNER.MATRIX_SHA256)
        self.assertEqual(
            RUNNER._canonical_sha256(RUNNER.generation_projection(1)),
            RUNNER.V1_SEMANTIC_PROJECTION_SHA256)
        self.assertEqual(
            RUNNER.generation_projection(1)["canonical_main_sha256"],
            RUNNER.HISTORICAL_EXECUTABLE_SHA256)
        self.assertEqual(
            RUNNER.generation_projection(2)["canonical_main_sha256"],
            RUNNER.EXACT_MAIN_EXECUTABLE_SHA256)
        self.assert_rejected(lambda: RUNNER.generation_projection(3))

    @unittest.skipUnless(
        (SEALED_AUTHORITY / "baseline-authority.json").is_file(),
        "sealed exact-main authority lane is not retained on this host")
    def test_retained_exact_main_authority_replays_to_the_pins(self) -> None:
        binding = RUNNER.bind_exact_main_authority(SEALED_AUTHORITY)
        self.assertEqual(
            binding["record_sha256"],
            RUNNER.EXACT_MAIN_AUTHORITY_RECORD_SHA256)
        self.assertEqual(
            binding["executable_sha256"],
            RUNNER.EXACT_MAIN_EXECUTABLE_SHA256)
        self.assertEqual(
            binding["combined_sha256"],
            RUNNER.EXACT_MAIN_COMBINED_SHA256)
        sealed_executable = SEALED_AUTHORITY / (
            "artifacts/canonical-first/leopard_main_benchmark")
        self.assertEqual(sealed_executable.stat().st_mode & 0o777, 0o400)
        self.assert_rejected(lambda:
            RUNNER.validate_exact_main_launch_copy(
                sealed_executable, binding))

    def test_main_launch_copy_is_executable_and_outside_sealed_lane(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-v2-main-launch-") as directory:
            root = Path(directory)
            lane = root / "authority"
            lane.mkdir()
            sealed = lane / "leopard_main_benchmark"
            sealed.write_bytes(b"synthetic exact-main authority bytes\n")
            sealed.chmod(0o400)
            launch_copy = root / "launch-copy"
            launch_copy.write_bytes(sealed.read_bytes())
            launch_copy.chmod(0o555)
            digest = RUNNER._sha256_file(launch_copy)
            binding = authority_binding()
            binding["lane_root"] = str(lane.resolve())
            binding["executable_sha256"] = digest
            with mock.patch.object(
                    RUNNER, "EXACT_MAIN_EXECUTABLE_SHA256", digest):
                self.assert_rejected(lambda:
                    RUNNER.validate_exact_main_launch_copy(sealed, binding))
                self.assertEqual(
                    RUNNER.validate_exact_main_launch_copy(
                        launch_copy, binding),
                    launch_copy.resolve())
                launch_copy.chmod(0o400)
                self.assert_rejected(lambda:
                    RUNNER.validate_exact_main_launch_copy(
                        launch_copy, binding))

    def test_main_child_requires_closed_pure_avx2_authority_build(self) -> None:
        cell = RUNNER.cells()[0]
        pristine = main_result(cell)
        normalized = BASE.validate_result(
            "main", pristine, cell, "e" * 40, "f" * 40,
            RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP)
        self.assertEqual(normalized["exact_main_authority_attribution"], {
            "pure_avx2": True,
            "main_source_commit": RUNNER.EXACT_MAIN_SOURCE_COMMIT,
            "authority_record_sha256":
                RUNNER.EXACT_MAIN_AUTHORITY_RECORD_SHA256,
            "executable_sha256": RUNNER.EXACT_MAIN_EXECUTABLE_SHA256,
        })
        for mutation in (
            lambda value: value["build"].pop("pure_avx2"),
            lambda value: value["build"].update({"pure_avx2": False}),
            lambda value: value["build"].update({"extra": False}),
            lambda value: value["build"].update(
                {"main_source_commit": "0" * 40}),
            lambda value: value["build"].update({"cplusplus": True}),
            lambda value: value.update({"schema": "leopard-main-benchmark-v2"}),
        ):
            malformed = copy.deepcopy(pristine)
            mutation(malformed)
            self.assert_rejected(lambda malformed=malformed:
                BASE.validate_result(
                    "main", malformed, cell, "e" * 40, "f" * 40,
                    RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP))

    def test_attempt_budget_is_gapless_and_outcome_independent(self) -> None:
        first = terminal(1)
        second = terminal(2, "failed")
        record = RUNNER.validate_attempt_preregistration(
            3, [first, second], verify_files=False)
        self.assertEqual(record["attempts_remaining"], 0)
        self.assertEqual([item["attempt"] for item in record["lineage"]], [1, 2])
        cases = (
            (4, [first, second, terminal(3)]),
            (3, [first]),
            (3, [second, first]),
            (2, [terminal(1, "accepted")]),
        )
        for attempt, lineage in cases:
            self.assert_rejected(lambda attempt=attempt, lineage=lineage:
                RUNNER.validate_attempt_preregistration(
                    attempt, lineage, verify_files=False))
        preliminary = terminal(1)
        preliminary["summary_schema"] = RUNNER.PRELIMINARY_SUMMARY_SCHEMA_V2
        self.assert_rejected(lambda: RUNNER.validate_attempt_preregistration(
            2, [preliminary], verify_files=False))

    def test_attempt_lineage_replays_retained_artifact_hashes(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-v2-lineage-") as directory:
            root = Path(directory) / "attempt-1"
            root.mkdir()
            summary = {
                "schema": RUNNER.FINAL_SUMMARY_SCHEMA_V2,
                "status": "rejected",
                "promotable": False,
            }
            raw = {"schema": RUNNER.RAW_SCHEMA_V2}
            BASE.T8_SUPPORT.write_exclusive(root / "summary.json", summary)
            BASE.T8_SUPPORT.write_exclusive(root / "raw.json", raw)
            retained = terminal(1)
            retained["output_root"] = str(root.resolve())
            retained["summary_sha256"] = RUNNER._sha256_file(
                root / "summary.json")
            retained["raw_sha256"] = RUNNER._sha256_file(root / "raw.json")
            BASE.T8_SUPPORT.write_exclusive(
                root / "attempt-terminal.json", retained)
            validated = RUNNER.validate_attempt_preregistration(
                2, [retained], verify_files=True)
            self.assertEqual(validated["lineage"], [retained])
            (root / "raw.json").write_bytes(b"{}\n")
            self.assert_rejected(lambda:
                RUNNER.validate_attempt_preregistration(
                    2, [retained], verify_files=True))

    def test_authority_binding_rejects_every_pinned_drift(self) -> None:
        pristine = authority_binding()
        mutations = (
            ("record_sha256", "0" * 64),
            ("ledger_sha256", "0" * 64),
            ("executable_sha256", RUNNER.HISTORICAL_EXECUTABLE_SHA256),
            ("archive_sha256", RUNNER.HISTORICAL_ARCHIVE_SHA256),
            ("combined_sha256", "0" * 64),
            ("source_commit", "0" * 40),
            ("source_tree", "0" * 40),
            ("adapter_commit", "0" * 40),
            ("pure_avx2", False),
        )
        for name, replacement in mutations:
            malformed = copy.deepcopy(pristine)
            malformed[name] = replacement
            self.assert_rejected(lambda malformed=malformed:
                RUNNER.validate_exact_main_authority_binding(malformed))
        relabeled = copy.deepcopy(pristine)
        relabeled["historical_non_authority"][0]["authority"] = True
        self.assert_rejected(lambda:
            RUNNER.validate_exact_main_authority_binding(relabeled))
        for lane_root in (
            "/frozen/../exact-main-authority",
            "/frozen/exact-main-authority/",
            "/",
        ):
            malformed = copy.deepcopy(pristine)
            malformed["lane_root"] = lane_root
            self.assert_rejected(lambda malformed=malformed:
                RUNNER.validate_exact_main_authority_binding(malformed))

        malformed_terminal = terminal(1)
        malformed_terminal["output_root"] = "/evidence/../k65-v2-a1"
        self.assert_rejected(lambda: RUNNER.validate_attempt_terminal(
            malformed_terminal, verify_files=False))

        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-v2-authority-link-") as directory:
            target = Path(directory) / "authority"
            target.mkdir()
            link = Path(directory) / "authority-link"
            link.symlink_to(target, target_is_directory=True)
            self.assert_rejected(lambda:
                RUNNER.bind_exact_main_authority(link))

    def test_only_complete_final_summary_is_promotable(self) -> None:
        raw, preliminary, authority, preregistration = complete_fixture()
        final = RUNNER.apply_final_acceptance(
            preliminary, raw, authority, preregistration)
        self.assertEqual(final["schema"], RUNNER.FINAL_SUMMARY_SCHEMA_V2)
        self.assertEqual(final["status"], "accepted")
        self.assertIs(final["promotable"], True)
        self.assertEqual(
            RUNNER.validate_promotable_final_summary(final), final)
        self.assert_rejected(lambda:
            RUNNER.validate_promotable_final_summary(preliminary))
        relabeled = copy.deepcopy(preliminary)
        relabeled["schema"] = RUNNER.FINAL_SUMMARY_SCHEMA_V2
        self.assert_rejected(lambda:
            RUNNER.validate_promotable_final_summary(relabeled))
        self.assert_rejected(lambda:
            RUNNER.apply_final_acceptance(
                final, raw, authority, preregistration))

        for mutate in (
            lambda value: value.update({"process_count": 0}),
            lambda value: value.update({"all_digests_matched": False}),
            lambda value: value.update(
                {"all_rounds_zero_sibling_nonidle": False}),
            lambda value: value["binary_sha256"].update(
                {"main": "0" * 64}),
            lambda value: value["binary_sha256"].update(
                {"candidate": RUNNER.EXACT_MAIN_EXECUTABLE_SHA256,
                 "control": RUNNER.EXACT_MAIN_EXECUTABLE_SHA256}),
            lambda value: value["mode_words"][
                "shared_binary_default"].update({"value": 2}),
        ):
            malformed = copy.deepcopy(final)
            mutate(malformed)
            self.assert_rejected(lambda malformed=malformed:
                RUNNER.validate_promotable_final_summary(malformed))

    def test_partial_or_reordered_campaign_cannot_finalize(self) -> None:
        mutators = (
            lambda raw: raw["cells"][0]["rounds"].pop(),
            lambda raw: raw["cells"].pop(),
            lambda raw: raw["cells"][0].update({"failed_round": {}}),
            lambda raw: raw["cells"][0]["rounds"][0]["isolation"].update(
                {"accepted": False}),
            lambda raw: raw["cells"][0]["rounds"][1].update({"round": 0}),
        )
        for mutate in mutators:
            raw, preliminary, authority, preregistration = complete_fixture()
            mutate(raw)
            self.assert_rejected(lambda raw=raw, preliminary=preliminary,
                                        authority=authority,
                                        preregistration=preregistration:
                RUNNER.apply_final_acceptance(
                    preliminary, raw, authority, preregistration))

    def test_discarded_full_orders_remain_valid_but_are_not_inference(self) -> None:
        cases = (
            (0, len(BASE.TARGET_ORDER[0])),
            (1, len(BASE.NEIGHBOR_ORDER[0])),
        )
        for cell_index, discarded_processes in cases:
            with self.subTest(cell_index=cell_index):
                raw, preliminary, authority, preregistration = \
                    complete_fixture()
                retained_round = raw["cells"][cell_index]["rounds"][0]
                discarded = copy.deepcopy(retained_round)
                discarded.pop("round")
                discarded.pop("discarded_attempts")
                discarded["isolation"]["accepted"] = False
                retained_round["attempt"] = 1
                retained_round["discarded_attempts"] = [discarded]
                preliminary["discarded_round_attempts"] = 1
                preliminary["discarded_process_count"] = discarded_processes
                final = RUNNER.apply_final_acceptance(
                    preliminary, raw, authority, preregistration)
                self.assertEqual(final["status"], "accepted")
                self.assertEqual(
                    RUNNER.validate_promotable_final_summary(final), final)
                for impossible_count in (3, 5):
                    malformed = copy.deepcopy(final)
                    malformed["discarded_process_count"] = impossible_count
                    self.assert_rejected(lambda malformed=malformed:
                        RUNNER.validate_promotable_final_summary(malformed))

    def test_threshold_failure_is_final_but_not_promotable(self) -> None:
        raw, preliminary, authority, preregistration = complete_fixture(
            target_control=1.04)
        final = RUNNER.apply_final_acceptance(
            preliminary, raw, authority, preregistration)
        self.assertEqual(final["status"], "rejected")
        self.assertIs(final["promotable"], False)
        self.assertGreater(len(final["gate_failures"]), 0)
        self.assert_rejected(lambda:
            RUNNER.validate_promotable_final_summary(final))

    def test_historical_hash_cannot_escape_non_authority_disclosure(self) -> None:
        raw, preliminary, authority, preregistration = complete_fixture()
        final = RUNNER.apply_final_acceptance(
            preliminary, raw, authority, preregistration)
        final["escaped_authority_sha256"] = RUNNER.HISTORICAL_EXECUTABLE_SHA256
        self.assert_rejected(lambda:
            RUNNER.validate_promotable_final_summary(final))

    def test_final_live_identity_revalidation_fails_closed(self) -> None:
        runner_identity = BASE.T8_SUPPORT.file_identity(BASE.RUNNER_PATH)
        dependency_identities = [
            BASE.support_file_identity(path)
            for path in BASE.RUNNER_DEPENDENCIES
        ]
        closure = {"schema": "test-build-closure"}
        raw = {
            "runner_after": runner_identity,
            "runner_dependencies_after": dependency_identities,
            "identities_after": {},
            "input_identities_after": {},
            "build_closure_after": closure,
        }
        options = argparse.Namespace()
        with mock.patch.object(
                BASE, "build_closure_identity", return_value=closure):
            RUNNER.validate_final_live_identities(raw, options)
            malformed = copy.deepcopy(raw)
            malformed["runner_after"]["sha256"] = "0" * 64
            self.assert_rejected(lambda:
                RUNNER.validate_final_live_identities(malformed, options))

    def test_main_atomically_replaces_preliminary_with_final_terminal(self) -> None:
        status, artifacts, _, _ = mocked_main_artifacts()
        self.assertEqual(status, 0)
        final = artifacts["summary.json"]
        self.assertEqual(
            RUNNER.validate_promotable_final_summary(final), final)
        retained = artifacts["preliminary-summary.json"]
        self.assertEqual(
            retained["schema"], RUNNER.PRELIMINARY_SUMMARY_SCHEMA_V2)
        self.assertNotIn("promotable", retained)
        terminal_record = artifacts["attempt-terminal.json"]
        self.assertEqual(terminal_record["outcome"], "accepted")
        self.assertIs(terminal_record["promotable"], True)
        self.assertNotIn("failure.json", artifacts)

    def test_promotion_gate_failure_is_published_only_as_failed(self) -> None:
        error = BASE.EvidenceError("synthetic final promotion rejection")
        status, artifacts, _, _ = mocked_main_artifacts(
            promotion_error=error)
        self.assertEqual(status, 1)
        self.assertEqual(artifacts["summary.json"]["status"], "rejected")
        self.assertIs(artifacts["summary.json"]["promotable"], False)
        self.assertEqual(
            artifacts["attempt-terminal.json"]["outcome"], "failed")
        self.assertIs(
            artifacts["attempt-terminal.json"]["promotable"], False)
        self.assertIn("failure.json", artifacts)

    def test_finalization_entry_failure_is_recorded_in_attempt_ledger(self) -> None:
        error = BASE.EvidenceError("synthetic finalization isolation failure")
        status, artifacts, _, _ = mocked_main_artifacts(
            finalization_entry_error=error)
        self.assertEqual(status, 1)
        self.assertEqual(artifacts["summary.json"]["status"], "rejected")
        self.assertIs(artifacts["summary.json"]["promotable"], False)
        self.assertEqual(
            artifacts["attempt-terminal.json"]["outcome"], "failed")
        self.assertIn("failure.json", artifacts)

    def test_preliminary_failure_without_journal_is_recorded(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-k65-v2-preliminary-failure-") as directory:
            output = Path(directory) / "attempt-1"
            authority = authority_binding()
            preregistration = RUNNER.validate_attempt_preregistration(
                1, [], verify_files=False)
            options = argparse.Namespace(
                output=output, cpu=7, sibling=3, attempt=1,
                exact_main_authority=authority,
                exact_main_authority_lane=Path(authority["lane_root"]),
                attempt_preregistration=preregistration,
            )

            def fail_without_journal() -> int:
                output.mkdir()
                return 1

            with mock.patch.object(
                    BASE, "parse_arguments", return_value=options), \
                    mock.patch.object(
                        BASE, "main", side_effect=fail_without_journal), \
                    mock.patch.object(
                        RUNNER, "_begin_isolated_finalization") as begin, \
                    contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                status = RUNNER.main()
            self.assertEqual(status, 1)
            begin.assert_not_called()
            failure = RUNNER._read_json_object(
                output / "failure.json", "test preliminary failure")
            self.assertEqual(failure["generation"], RUNNER.GENERATION)
            terminal_record = RUNNER._read_json_object(
                output / "attempt-terminal.json", "test attempt terminal")
            self.assertEqual(terminal_record["outcome"], "failed")
            self.assertEqual(terminal_record["summary_schema"], "absent")


if __name__ == "__main__":
    unittest.main()
