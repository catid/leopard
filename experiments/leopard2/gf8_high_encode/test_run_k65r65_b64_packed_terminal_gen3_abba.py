#!/usr/bin/env python3
"""Plan-only and state-machine tests for the K65 generation-3 runner shell."""

from __future__ import annotations

import copy
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
RUNNER_PATH = HERE / "run_k65r65_b64_packed_terminal_gen3_abba.py"
V2_PATH = HERE / "run_k65r65_b64_packed_terminal_abba.py"
TEMPLATE_PATH = HERE / \
    "k65r65_b64_packed_terminal_gen3_preregistration.template.json"
PLAN_SHA256 = \
    "415415cf72be1756b8cb8c3ec69d41fae605e0bd7fab66eaa41bf23a5c49c399"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


runner = load_module("k65_gen3_plan_runner_under_test", RUNNER_PATH)


def final_preregistration() -> dict:
    return runner.prereg.preregistration_record(
        authority="unit-test-fixture-not-an-authorization",
        authorized_utc="2026-08-30T00:00:00Z",
        clock_ticks_per_second=100,
        candidate_primary_cpus=[8, 9, 10, 11, 12, 13, 14, 15],
        excluded_pairs=[],
        track_b_permitted=False,
        setup_invalid_budget=5,
        environment_rejected_budget=8,
        evidence_attempt_budget=3,
        scan_window_count=60,
        scan_nominal_window_ns=1_000_000_000,
        bridge_minimum_window_count=5,
        bridge_maximum_window_count=120,
        bridge_nominal_window_ns=1_000_000_000,
        maximum_handoff_elapsed_ns=120_000_000_000,
        freeze_point="armed",
        candidate_executable_mode="frozen-sha256",
        candidate_executable_sha256="a" * 64,
        candidate_executable_size=1_234_567,
        candidate_build_provenance_sha256="4a" * 32,
        candidate_reproducible_build_core_sha256="5b" * 32,
        candidate_authority_record_sha256="6c" * 32,
        candidate_authority_ledger_sha256="7d" * 32,
        candidate_source_commit="1" * 40,
        candidate_source_tree="2" * 40,
        host_machine_id_sha256="b" * 64,
        host_name="qualification-host.example",
        host_architecture="x86_64",
        host_cpu_model="Unit Test CPU",
        output_lane_binding={
            "schema": runner.prereg.OUTPUT_LANE_SCHEMA,
            "path": "/unit/k65-generation-3-output-lane",
            "device": 1,
            "inode": 1,
            "uid": os.geteuid(),
            "mode": 0o500,
            "link_count": 6,
            "initial_mtime_ns": 1,
            "initial_ctime_ns": 1,
            "file_handle": {
                "schema": runner.prereg.OUTPUT_LANE_FILE_HANDLE_SCHEMA,
                "handle_type": 1,
                "handle_hex": "00",
            },
            "lane_manifest": {
                "schema":
                    runner.prereg.OUTPUT_LANE_MANIFEST_BINDING_SCHEMA,
                "name": runner.prereg.OUTPUT_LANE_MANIFEST_FILE,
                "device": 1,
                "inode": 2,
                "uid": os.geteuid(),
                "mode": 0o400,
                "link_count": 1,
                "initial_mtime_ns": 1,
                "initial_ctime_ns": 1,
                "sha256": "0" * 64,
                "size": 1,
                "file_handle": {
                    "schema": runner.prereg.OUTPUT_LANE_FILE_HANDLE_SCHEMA,
                    "handle_type": 1,
                    "handle_hex": "00",
                },
            },
        },
        child_launch_context=
            runner.prereg.recommended_launch_context_record(),
        controller_bindings=[
            {"path": path, "sha256": f"{index + 1:064x}"}
            for index, path in enumerate(runner.prereg.REQUIRED_CONTROLLER_PATHS)
        ],
        verify_files=False,
    )


def template() -> dict:
    return runner.prereg.load_preregistration_template(
        TEMPLATE_PATH.read_bytes(), verify_files=False)


PAIR = {"benchmark_cpu": 8, "reserved_sibling": 72}
class K65Generation3PlanRunnerTests(unittest.TestCase):
    def assertRejected(self, function, pattern: str | None = None) -> None:
        with self.assertRaises((
            runner.PlanError,
            runner.prereg.PreregistrationError,
            runner.contract.QualificationError,
        )) as raised:
            function()
        if pattern is not None:
            self.assertIn(pattern, str(raised.exception))

    def test_plan_contains_all_1650_exact_children_without_execution(self) -> None:
        plan = runner.campaign_plan_record(
            preregistration=template())
        self.assertFalse(plan["safe_to_execute"])
        self.assertFalse(plan["candidate_timing_performed"])
        self.assertFalse(plan["preregistration_ratified"])
        self.assertNotIn("selected_pair", plan)
        self.assertNotIn("artifacts", plan)
        self.assertEqual(plan["child_process_count"], 1650)
        self.assertEqual(len(plan["child_plans"]), 1650)
        self.assertEqual(
            [sum(child["cell_index"] == index for child in plan["child_plans"])
             for index in range(13)],
            [150, 100, 100, 100, 100, 100, 100, 150, 150, 150, 150, 150, 150])
        self.assertEqual(
            runner.contract.canonical_sha256(plan), PLAN_SHA256)
        self.assertEqual(
            runner.validate_campaign_plan(plan, template()), plan)

    def test_every_child_command_and_round_order_match_frozen_v2(self) -> None:
        v2 = load_module("k65_v2_for_gen3_plan_comparison", V2_PATH)
        self.assertEqual(runner.cells(), v2.cells())
        plans = runner.child_plans()
        index = 0
        for cell, expected_cell in zip(runner.cells(), v2.cells()):
            self.assertEqual(cell, expected_cell)
            base_orders = v2.BASE.TARGET_ORDER if cell["compare_main"] \
                else v2.BASE.NEIGHBOR_ORDER
            orders = v2.select_round_orders(base_orders, 25)
            for round_index, order in enumerate(orders):
                for slot, implementation in enumerate(order):
                    item = plans[index]
                    self.assertEqual(
                        (item["cell_id"], item["round"], item["slot"],
                         item["implementation"]),
                        (cell["id"], round_index, slot, implementation))
                    self.assertEqual(
                        item["argv_tail"],
                        v2.benchmark_command(
                            implementation, Path("/logical/executable"),
                            cell, 8, 31, 64)[6:])
                    self.assertEqual(
                        item["timeout_budget"],
                        v2.child_timeout_budget(implementation, cell, 31, 64))
                    index += 1
        self.assertEqual(index, 1650)

    def test_preregistration_campaign_must_match_frozen_plan_constants(self) -> None:
        base = template()
        mutations = {
            "generation": 4,
            "matrix_sha256": "0" * 64,
            "cell_count": 14,
            "rounds_per_cell": 26,
            "iterations_per_child": 41,
            "warmup_per_child": 65,
            "reuse_per_child": 8193,
            "expected_child_process_count": 1651,
        }
        for key, value in mutations.items():
            with self.subTest(key=key):
                campaign = copy.deepcopy(base["campaign"])
                campaign[key] = value
                changed = copy.deepcopy(base)
                changed["campaign"] = campaign
                with mock.patch.object(
                    runner.prereg, "campaign_contract_record",
                    return_value=copy.deepcopy(campaign),
                ):
                    self.assertRejected(lambda changed=changed: (
                        runner.campaign_plan_record(
                            preregistration=changed)), "frozen plan constants")

    def test_state_machine_accepts_only_one_way_paths(self) -> None:
        success = runner.attempt_state_history_record(
            lane_class="evidence", lane_index=1,
            states=(
                "INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
                "BRIDGING", "BRIDGED", "ARMING", "PRESAMPLING", "ARMED", "TIMING",
                "FINALIZING", "ACCEPTED",
            ),
            selected_pair=PAIR)
        self.assertEqual(runner.validate_attempt_state_history(success), success)
        self.assertTrue(success["evidence_attempt_committed"])
        self.assertEqual(success["budget_committed"], "evidence_attempt")
        self.assertEqual(success["terminal_state"], "ACCEPTED")
        self.assertEqual(
            [event["selected_pair"] for event in success["events"][:8]],
            [None] * 8)
        self.assertTrue(all(
            event["selected_pair"] == PAIR for event in success["events"][8:]))

        invalid = (
            ("INIT", "TIMING"),
            ("INIT", "PREREGISTERED", "ACCEPTED"),
            ("INIT", "PREREGISTERED", "QUALIFYING", "ENV_REJECTED", "TIMING"),
            ("INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
             "BRIDGING", "BRIDGED", "ARMING", "ARMED", "ACCEPTED"),
            ("INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
             "BRIDGING", "BRIDGED", "ARMING", "ARMED", "TIMING", "ACCEPTED"),
        )
        for states in invalid:
            with self.subTest(states=states):
                self.assertRejected(lambda states=states: (
                    runner.attempt_state_history_record(
                        lane_class="evidence", lane_index=1, states=states,
                        selected_pair=PAIR)))

    def test_budget_partition_and_crash_after_timing_are_structural(self) -> None:
        registration = final_preregistration()
        setup = runner.attempt_state_history_record(
            lane_class="setup", lane_index=1,
            states=(
                "INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
                "BRIDGING", "BRIDGED", "ARMING", "SETUP_INVALID",
            ), selected_pair=None)
        environment = runner.attempt_state_history_record(
            lane_class="environment", lane_index=1,
            states=("INIT", "PREREGISTERED", "QUALIFYING", "ENV_REJECTED"),
            selected_pair=None)
        killed_after_timing = runner.attempt_state_history_record(
            lane_class="evidence", lane_index=1,
            states=(
                "INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
                "BRIDGING", "BRIDGED", "ARMING", "PRESAMPLING", "ARMED", "TIMING",
            ), selected_pair=PAIR)
        accepted = runner.attempt_state_history_record(
            lane_class="evidence", lane_index=2,
            states=(
                "INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
                "BRIDGING", "BRIDGED", "ARMING", "PRESAMPLING", "ARMED", "TIMING",
                "FINALIZING", "ACCEPTED",
            ), selected_pair=PAIR)
        ledger = runner.budget_ledger_record(
            registration, [setup, environment, killed_after_timing, accepted])
        self.assertEqual(ledger["setup_invalid_used"], 1)
        self.assertEqual(ledger["environment_rejected_used"], 1)
        self.assertEqual(ledger["evidence_attempts_used"], 2)
        self.assertTrue(killed_after_timing["evidence_attempt_committed"])
        self.assertIsNone(killed_after_timing["terminal_state"])
        self.assertEqual(runner.validate_budget_ledger(
            ledger, registration), ledger)

        abandoned_qualification = runner.attempt_state_history_record(
            lane_class="environment", lane_index=2,
            states=("INIT", "PREREGISTERED", "QUALIFYING"),
            selected_pair=None)
        abandoned_armed = runner.attempt_state_history_record(
            lane_class="evidence", lane_index=3,
            states=(
                "INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
                "BRIDGING", "BRIDGED", "ARMING", "PRESAMPLING", "ARMED",
            ), selected_pair=PAIR)
        abandoned_presampling = runner.attempt_state_history_record(
            lane_class="environment", lane_index=3,
            states=(
                "INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
                "BRIDGING", "BRIDGED", "ARMING", "PRESAMPLING",
            ), selected_pair=None)
        self.assertEqual(abandoned_qualification["budget_committed"],
                         "environment_rejected")
        self.assertEqual(abandoned_presampling["budget_committed"],
                         "environment_rejected")
        self.assertEqual(abandoned_armed["budget_committed"], "evidence_attempt")
        self.assertFalse(abandoned_qualification["evidence_attempt_committed"])
        self.assertTrue(abandoned_armed["evidence_attempt_committed"])
        self.assertTrue(abandoned_armed["events"][-1][
            "evidence_attempt_committed"])

    def test_plan_has_no_artifact_or_pair_injection_surface(self) -> None:
        plan = runner.campaign_plan_record(preregistration=template())
        self.assertNotIn("selected_pair", plan)
        self.assertNotIn("artifacts", plan)
        self.assertTrue(all("command" not in child for child in plan["child_plans"]))
        self.assertTrue(all(child["argv_tail"][0] == "--k"
                            for child in plan["child_plans"]))

    def test_frozen_pair_is_write_once_and_budget_limits_fail_closed(self) -> None:
        registration = final_preregistration()
        histories = []
        for index in range(1, 4):
            histories.append(runner.attempt_state_history_record(
                lane_class="evidence", lane_index=index,
                states=(
                    "INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
                    "BRIDGING", "BRIDGED", "ARMING", "PRESAMPLING", "ARMED", "TIMING",
                    "REJECTED",
                ), selected_pair=PAIR))
        ledger = runner.budget_ledger_record(registration, histories)
        self.assertEqual(ledger["evidence_attempts_used"], 3)

        fourth = runner.attempt_state_history_record(
            lane_class="evidence", lane_index=4,
            states=(
                "INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED",
                "BRIDGING", "BRIDGED", "ARMING", "PRESAMPLING", "ARMED", "TIMING",
            ), selected_pair=PAIR)
        self.assertRejected(lambda: runner.budget_ledger_record(
            registration, [*histories, fourth]), "budget is exhausted")

        changed = copy.deepcopy(histories[1])
        other = {"benchmark_cpu": 9, "reserved_sibling": 73}
        changed = runner.attempt_state_history_record(
            lane_class="evidence", lane_index=2,
            states=[event["to_state"] for event in changed["events"]],
            selected_pair=other)
        self.assertRejected(lambda: runner.budget_ledger_record(
            registration, [histories[0], changed]), "frozen pair changed")

    def test_cli_is_canonical_plan_only_and_never_ratifies_template(self) -> None:
        command = [
            sys.executable, "-I", "-S", "-B", str(RUNNER_PATH),
            "--plan-only", "--preregistration", str(TEMPLATE_PATH),
        ]
        completed = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False, timeout=30)
        self.assertEqual(completed.returncode, 0, completed.stderr.decode())
        self.assertEqual(completed.stderr, b"")
        value = runner.contract.strict_json_loads(completed.stdout, "CLI plan")
        self.assertEqual(completed.stdout, runner.contract.canonical_json_bytes(value))
        self.assertFalse(value["preregistration_ratified"])
        self.assertFalse(value["safe_to_execute"])
        self.assertEqual(value["child_process_count"], 1650)

    def test_plan_runner_has_no_acquisition_or_execution_primitives(self) -> None:
        source = RUNNER_PATH.read_text("utf-8")
        for forbidden in (
            "import os", "import time", "sched_setaffinity", "sched_getaffinity",
            '"/proc/', '"/sys/', "subprocess", "Popen", "check_output",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, source)
        self.assertIn("The only CLI mode is ``--plan-only``", source)

    def test_loading_plan_runner_does_not_mutate_frozen_v2(self) -> None:
        v2 = load_module("k65_v2_before_second_gen3_load", V2_PATH)
        before = {
            "generation": v2.GENERATION,
            "attempt_budget": v2.ATTEMPT_BUDGET,
            "raw_schema": v2.RAW_SCHEMA_V2,
            "summary_schema": v2.FINAL_SUMMARY_SCHEMA_V2,
            "matrix": v2.MATRIX_SHA256,
            "projection": v2._canonical_sha256(v2.generation_projection(2)),
            "base_id": id(v2.BASE),
            "cells": copy.deepcopy(v2.cells()),
        }
        second = load_module("k65_gen3_plan_runner_second_load", RUNNER_PATH)
        after = {
            "generation": v2.GENERATION,
            "attempt_budget": v2.ATTEMPT_BUDGET,
            "raw_schema": v2.RAW_SCHEMA_V2,
            "summary_schema": v2.FINAL_SUMMARY_SCHEMA_V2,
            "matrix": v2.MATRIX_SHA256,
            "projection": v2._canonical_sha256(v2.generation_projection(2)),
            "base_id": id(v2.BASE),
            "cells": copy.deepcopy(v2.cells()),
        }
        self.assertEqual(after, before)
        self.assertFalse(hasattr(second, "BASE"))
        with self.assertRaises(v2.BASE.EvidenceError):
            v2.generation_projection(3)


if __name__ == "__main__":
    unittest.main()
