#!/usr/bin/env python3
"""Host-independent adversarial tests for K65 generation-3 execution authority."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import os
from pathlib import Path
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "k65_gen3_execution_contract.py"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


execution = load_module("k65_gen3_execution_contract_under_test", MODULE_PATH)
prereg = execution.prereg
contract = execution.contract
bridge_contract = execution.bridge_contract

CANDIDATE_SHA256 = "a" * 64
CANDIDATE_SIZE = 2_000_000
BUILD_PROVENANCE_SHA256 = "4a" * 32
REPRODUCIBLE_CORE_SHA256 = "5b" * 32
CANDIDATE_AUTHORITY_RECORD_SHA256 = "6c" * 32
CANDIDATE_AUTHORITY_LEDGER_SHA256 = "7d" * 32


def preregistration(
    candidate_sha256: str = CANDIDATE_SHA256, *,
    candidate_size: int = CANDIDATE_SIZE,
    candidate_build_provenance_sha256: str = BUILD_PROVENANCE_SHA256,
    candidate_reproducible_build_core_sha256: str = REPRODUCIBLE_CORE_SHA256,
    candidate_authority_record_sha256: str =
        CANDIDATE_AUTHORITY_RECORD_SHA256,
    candidate_authority_ledger_sha256: str =
        CANDIDATE_AUTHORITY_LEDGER_SHA256,
    controller_bindings: list[dict] | None = None,
    output_lane_binding: dict | None = None,
) -> dict:
    controllers = controller_bindings if controller_bindings is not None else [
        {
            "path": path,
            "sha256": hashlib.sha256(path.encode("ascii")).hexdigest(),
        }
        for path in prereg.REQUIRED_CONTROLLER_PATHS
    ]
    handle = {
        "schema": prereg.OUTPUT_LANE_FILE_HANDLE_SCHEMA,
        "handle_type": 1,
        "handle_hex": "00",
    }
    lane_binding = output_lane_binding if output_lane_binding is not None else {
        "schema": prereg.OUTPUT_LANE_SCHEMA,
        "path": "/unit/k65-generation-3-output-lane",
        "device": 1,
        "inode": 1,
        "uid": os.geteuid(),
        "mode": 0o500,
        "link_count": 6,
        "initial_mtime_ns": 1,
        "initial_ctime_ns": 1,
        "file_handle": copy.deepcopy(handle),
        "lane_manifest": {
            "schema": prereg.OUTPUT_LANE_MANIFEST_BINDING_SCHEMA,
            "name": prereg.OUTPUT_LANE_MANIFEST_FILE,
            "device": 1,
            "inode": 2,
            "uid": os.geteuid(),
            "mode": 0o400,
            "link_count": 1,
            "initial_mtime_ns": 1,
            "initial_ctime_ns": 1,
            "sha256": "0" * 64,
            "size": 1,
            "file_handle": copy.deepcopy(handle),
        },
    }
    return prereg.preregistration_record(
        authority="host-independent unit-test fixture",
        authorized_utc="2026-08-31T00:00:00Z",
        clock_ticks_per_second=100,
        candidate_primary_cpus=[2], excluded_pairs=[],
        track_b_permitted=False,
        setup_invalid_budget=5, environment_rejected_budget=8,
        evidence_attempt_budget=3,
        scan_window_count=2, scan_nominal_window_ns=250_000_000,
        bridge_minimum_window_count=2, bridge_maximum_window_count=3,
        bridge_nominal_window_ns=250_000_000,
        maximum_handoff_elapsed_ns=2_000_000_000,
        freeze_point="armed", candidate_executable_mode="frozen-sha256",
        candidate_executable_sha256=candidate_sha256,
        candidate_executable_size=candidate_size,
        candidate_build_provenance_sha256=candidate_build_provenance_sha256,
        candidate_reproducible_build_core_sha256=
            candidate_reproducible_build_core_sha256,
        candidate_authority_record_sha256=candidate_authority_record_sha256,
        candidate_authority_ledger_sha256=candidate_authority_ledger_sha256,
        candidate_source_commit="1" * 40,
        candidate_source_tree="2" * 40,
        host_machine_id_sha256="3" * 64,
        host_name="foureyes.lan", host_architecture="x86_64",
        host_cpu_model="AMD Ryzen Threadripper PRO 9985WX 64-Cores",
        output_lane_binding=lane_binding,
        child_launch_context=prereg.recommended_launch_context_record(),
        controller_bindings=controllers, verify_files=False,
    )


def qualification_records(
    registration: dict, expected_frozen_pair: dict | None = None,
) -> tuple[dict, dict, dict]:
    """Retarget the verifier's pure snapshot fixture to preregistered Track A."""
    acquisition_data, bridge_data, unused_policy, geometry = \
        execution.pair_verifier._self_test_fixture()
    old_acquisition = contract.strict_json_loads(acquisition_data)
    old_bridge = contract.strict_json_loads(bridge_data)
    policy = registration["qualification"]["policies"][0]
    old_scan = old_acquisition["scan"]
    scan_endpoints = [old_scan["windows"][0]["before"]] + [
        window["after"] for window in old_scan["windows"]
    ]
    scan = contract.pair_qualification_scan_record(
        policy,
        allowed_cpu_set_at_launch=old_scan["allowed_cpu_set_at_launch"],
        topology_before=old_scan["topology_before"],
        topology_after=old_scan["topology_after"],
        snapshots=scan_endpoints,
        frozen_pair_from_prior_attempt=expected_frozen_pair,
    )
    acquisition = {
        "schema": bridge_contract.ACQUISITION_SCHEMA,
        "acquisition_method": bridge_contract.ACQUISITION_METHOD,
        "sources": dict(bridge_contract._SOURCE_ITEMS),
        "policy": copy.deepcopy(policy),
        "policy_sha256": contract.canonical_sha256(policy),
        "requested_window_count": 2,
        "nominal_window_ns": 250_000_000,
        "frozen_pair_from_prior_attempt": copy.deepcopy(expected_frozen_pair),
        "allowed_cpu_set_at_launch": [2, 6],
        "allowed_cpu_set_after_scan": [2, 6],
        "clock_ticks_per_second_at_launch": 100,
        "clock_ticks_per_second_after_scan": 100,
        "topology_before_sha256": contract.canonical_sha256(
            scan["topology_before"]),
        "topology_after_sha256": contract.canonical_sha256(
            scan["topology_after"]),
        "scan": scan,
        "scan_sha256": contract.canonical_sha256(scan),
        "host_mutation_performed": False,
        "candidate_timing_performed": False,
    }
    bridge_endpoints = [old_bridge["windows"][0]["before"]] + [
        window["after"] for window in old_bridge["windows"]
    ]
    bridge = bridge_contract.pair_qualification_bridge_record(
        acquisition, expected_policy=policy,
        expected_frozen_pair=expected_frozen_pair,
        expected_acquisition_window_count=2,
        expected_acquisition_nominal_window_ns=250_000_000,
        expected_bridge_geometry=geometry, snapshots=bridge_endpoints,
    )
    verdict = execution.pair_verifier.require_accepted_pair_qualification_bundle(
        contract.canonical_json_bytes(acquisition),
        contract.canonical_json_bytes(bridge),
        expected_policy=policy,
        expected_policy_sha256=contract.canonical_sha256(policy),
        expected_frozen_pair=expected_frozen_pair,
        expected_acquisition_window_count=2,
        expected_acquisition_nominal_window_ns=250_000_000,
        expected_bridge_geometry=geometry,
    )
    return acquisition, bridge, verdict


def armed_closure_kwargs(bridge: dict) -> dict:
    return {
        "authority_bundle_sha256": "8a" * 32,
        "attempt_manifest_sha256": "9b" * 32,
        "lane_binding_sha256": "ac" * 32,
        "armed_monotonic_ns": bridge["bridge_finished_monotonic_ns"],
    }


def fixture(
    candidate_sha256: str = CANDIDATE_SHA256, *,
    candidate_size: int = CANDIDATE_SIZE,
    candidate_build_provenance_sha256: str = BUILD_PROVENANCE_SHA256,
    candidate_reproducible_build_core_sha256: str = REPRODUCIBLE_CORE_SHA256,
    candidate_authority_record_sha256: str =
        CANDIDATE_AUTHORITY_RECORD_SHA256,
    candidate_authority_ledger_sha256: str =
        CANDIDATE_AUTHORITY_LEDGER_SHA256,
    controller_bindings: list[dict] | None = None,
    output_lane_binding: dict | None = None,
    expected_frozen_pair: dict | None = None,
) -> dict:
    registration = preregistration(
        candidate_sha256,
        candidate_size=candidate_size,
        candidate_build_provenance_sha256=
            candidate_build_provenance_sha256,
        candidate_reproducible_build_core_sha256=
            candidate_reproducible_build_core_sha256,
        candidate_authority_record_sha256=
            candidate_authority_record_sha256,
        candidate_authority_ledger_sha256=
            candidate_authority_ledger_sha256,
        controller_bindings=controller_bindings,
        output_lane_binding=output_lane_binding)
    acquisition, bridge, verdict = qualification_records(
        registration, expected_frozen_pair)
    authority = execution.host_authority_record(
        machine_id_sha256="3" * 64, hostname="foureyes.lan",
        architecture="x86_64",
        cpu_model="AMD Ryzen Threadripper PRO 9985WX 64-Cores",
    )
    host = execution.host_instance_record(
        authority=authority,
        boot_id="12345678-1234-1234-1234-123456789abc",
        kernel_release="6.8.0-unit-test", online_cpus=[2, 6],
        allowed_cpus=[2, 6], clock_ticks_per_second=100,
        topology_sha256=acquisition["topology_before_sha256"],
    )
    selection = prereg.qualification_track_selection_record(
        registration, policy_a_scan=acquisition["scan"],
        policy_b_scan=None, expected_frozen_pair=expected_frozen_pair)
    evidence = {
        "policy_a_scan": acquisition["scan"],
        "policy_b_scan": None,
        "track_selection": selection,
        "expected_frozen_pair": copy.deepcopy(expected_frozen_pair),
        "acquisition_data": contract.canonical_json_bytes(acquisition),
        "bridge_data": contract.canonical_json_bytes(bridge),
        "independent_verdict_data": contract.canonical_json_bytes(verdict),
    }
    qualification = execution.qualification_binding_record(
        registration, host, evidence=evidence)
    roles = {
        "candidate": {
            "role": "candidate", "raw_sha256": candidate_sha256,
            "size": CANDIDATE_SIZE,
            "launch_protocol": execution.CANDIDATE_LAUNCH_PROTOCOL,
            "handle_id": "4" * 64,
            "handle_device": 10, "handle_inode": 20,
        },
        "control": {
            "role": "control", "raw_sha256": candidate_sha256,
            "size": CANDIDATE_SIZE,
            "launch_protocol": execution.CONTROL_LAUNCH_PROTOCOL,
            "handle_id": "4" * 64,
            "handle_device": 10, "handle_inode": 20,
        },
        "main": {
            "role": "main",
            "raw_sha256": execution.EXACT_MAIN_PATH_VARIANT_RAW_SHA256,
            "size": execution.EXACT_MAIN_PATH_VARIANT_SIZE,
            "launch_protocol": execution.EXACT_MAIN_LAUNCH_PROTOCOL,
            "handle_id": "5" * 64,
            "handle_device": 10, "handle_inode": 21,
        },
    }
    exact_main = execution.exact_main_authority_record(
        verifier_verdict_sha256="6" * 64)
    source_authority = execution.candidate_source_authority_record(
        registration,
        build_provenance_sha256=candidate_build_provenance_sha256,
        reproducible_build_core_sha256=
            candidate_reproducible_build_core_sha256,
        authority_record_sha256=candidate_authority_record_sha256,
        authority_ledger_sha256=candidate_authority_ledger_sha256,
        host_authority=host["authority"])
    artifacts = execution.artifact_bundle_record(
        registration, roles=roles, exact_main_authority=exact_main,
        candidate_source_authority=source_authority,
        host_authority=host["authority"])
    plan = execution.plan_contract.campaign_plan_record(
        preregistration=registration)
    armed = execution.armed_record(
        registration, plan, host, artifacts, qualification,
        qualification_evidence=evidence, evidence_attempt=1,
        prior_armed_chain_sha256=None,
        **armed_closure_kwargs(bridge))
    return {
        "preregistration": registration,
        "acquisition": acquisition,
        "bridge": bridge,
        "verdict": verdict,
        "host": host,
        "selection": selection,
        "evidence": evidence,
        "qualification": qualification,
        "roles": roles,
        "exact_main": exact_main,
        "source_authority": source_authority,
        "artifacts": artifacts,
        "plan": plan,
        "armed": armed,
    }


class K65Generation3ExecutionContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.base = fixture()

    def fresh(self) -> dict:
        return copy.deepcopy(self.base)

    def assertRejected(self, function, pattern: str | None = None) -> None:
        with self.assertRaises(ValueError) as raised:
            function()
        if pattern is not None:
            self.assertIn(pattern, str(raised.exception))

    def test_host_authority_and_instance_are_exact_fixed_points(self) -> None:
        data = self.fresh()
        host = data["host"]
        self.assertEqual(execution.validate_host_instance(host), host)
        self.assertEqual(
            execution.validate_host_authority(host["authority"]),
            host["authority"])
        self.assertEqual(host["authority"]["architecture"], "x86_64")
        self.assertEqual(host["allowed_cpus"], [2, 6])

        mutations = []
        extra = copy.deepcopy(host)
        extra["unexpected"] = False
        mutations.append(extra)
        bool_cpu = copy.deepcopy(host)
        bool_cpu["online_cpus"][0] = False
        mutations.append(bool_cpu)
        unsorted = copy.deepcopy(host)
        unsorted["allowed_cpus"] = [6, 2]
        mutations.append(unsorted)
        outside = copy.deepcopy(host)
        outside["allowed_cpus"] = [2, 7]
        mutations.append(outside)
        bad_boot = copy.deepcopy(host)
        bad_boot["boot_id"] = "NOT-A-UUID"
        mutations.append(bad_boot)
        wrong_arch = copy.deepcopy(host)
        wrong_arch["authority"]["architecture"] = "aarch64"
        mutations.append(wrong_arch)
        uppercase_host = copy.deepcopy(host)
        uppercase_host["authority"]["hostname"] = "FourEyes.lan"
        mutations.append(uppercase_host)
        for index, mutation in enumerate(mutations):
            with self.subTest(index=index):
                self.assertRejected(
                    lambda mutation=mutation:
                        execution.validate_host_instance(mutation))

    def test_artifact_bundle_binds_one_candidate_control_handle(self) -> None:
        data = self.fresh()
        bundle = data["artifacts"]
        registration = data["preregistration"]
        self.assertEqual(
            execution.validate_artifact_bundle(
                bundle, registration,
                host_authority=data["host"]["authority"]), bundle)
        roles = bundle["roles"]
        self.assertTrue(bundle["candidate_control_same_handle"])
        self.assertEqual(roles["candidate"]["handle_id"],
                         roles["control"]["handle_id"])
        self.assertEqual(roles["candidate"]["raw_sha256"],
                         registration["candidate_executable"]["sha256"])
        self.assertNotEqual(roles["main"]["handle_id"],
                            roles["candidate"]["handle_id"])

        mutations = []
        split_handle = copy.deepcopy(bundle)
        split_handle["roles"]["control"]["handle_id"] = "7" * 64
        mutations.append(split_handle)
        split_bytes = copy.deepcopy(bundle)
        split_bytes["roles"]["control"]["raw_sha256"] = "8" * 64
        mutations.append(split_bytes)
        split_size = copy.deepcopy(bundle)
        split_size["roles"]["control"]["size"] += 1
        mutations.append(split_size)
        reratified_size = copy.deepcopy(bundle)
        reratified_size["roles"]["candidate"]["size"] += 1
        reratified_size["roles"]["control"]["size"] += 1
        mutations.append(reratified_size)
        same_main = copy.deepcopy(bundle)
        same_main["roles"]["main"]["handle_id"] = "4" * 64
        mutations.append(same_main)
        canonical_main = copy.deepcopy(bundle)
        canonical_main["roles"]["main"]["raw_sha256"] = \
            prereg.EXACT_MAIN_EXECUTABLE_SHA256
        canonical_main["exact_main_authority"]["raw_sha256"] = \
            prereg.EXACT_MAIN_EXECUTABLE_SHA256
        mutations.append(canonical_main)
        normalized_drift = copy.deepcopy(bundle)
        normalized_drift["exact_main_authority"][
            "normalized_combined_sha256"] = "9" * 64
        mutations.append(normalized_drift)
        source_drift = copy.deepcopy(bundle)
        source_drift["candidate_source"]["tree"] = "3" * 40
        mutations.append(source_drift)
        controller_drift = copy.deepcopy(bundle)
        controller_drift["controller_bindings_sha256"] = "b" * 64
        mutations.append(controller_drift)
        bool_claim = copy.deepcopy(bundle)
        bool_claim["candidate_control_same_handle"] = 1
        mutations.append(bool_claim)
        for index, mutation in enumerate(mutations):
            with self.subTest(index=index):
                self.assertRejected(lambda mutation=mutation:
                    execution.validate_artifact_bundle(
                        mutation, registration,
                        host_authority=data["host"]["authority"]))

    def test_every_role_is_checked_against_the_generation2_denylist(self) -> None:
        self.assertRejected(
            lambda: fixture(prereg.EXACT_MAIN_EXECUTABLE_SHA256),
            "generation-2 artifact")

    def test_qualification_replays_verifier_and_all_exact_hash_joins(self) -> None:
        data = self.fresh()
        binding = data["qualification"]
        observed = execution.validate_qualification_binding(
            binding, data["preregistration"], data["host"],
            evidence=data["evidence"])
        self.assertEqual(observed, binding)
        self.assertTrue(binding["bridge_accepted"])
        self.assertFalse(binding["candidate_timing_performed"])
        self.assertEqual(binding["selected_pair"],
                         {"benchmark_cpu": 2, "reserved_sibling": 6})
        self.assertEqual(
            binding["campaign_presample_before_sha256"],
            data["bridge"]["campaign_presample_before_sha256"])

        mutations = []
        noncanonical_acquisition = copy.deepcopy(data["evidence"])
        noncanonical_acquisition["acquisition_data"] += b" "
        mutations.append((data["host"], noncanonical_acquisition))
        false_verdict = copy.deepcopy(data["evidence"])
        verdict = copy.deepcopy(data["verdict"])
        verdict["bridge_accepted"] = False
        false_verdict["independent_verdict_data"] = \
            contract.canonical_json_bytes(verdict)
        mutations.append((data["host"], false_verdict))
        different_scan = copy.deepcopy(data["evidence"])
        different_scan["policy_a_scan"]["scan_finished_monotonic_ns"] += 1
        mutations.append((data["host"], different_scan))
        wrong_track_selection = copy.deepcopy(data["evidence"])
        wrong_track_selection["track_selection"]["selected_pair"] = {
            "benchmark_cpu": 6, "reserved_sibling": 2}
        mutations.append((data["host"], wrong_track_selection))
        wrong_host_topology = copy.deepcopy(data["host"])
        wrong_host_topology["topology_sha256"] = "7" * 64
        mutations.append((wrong_host_topology, data["evidence"]))
        wrong_host_allowed = copy.deepcopy(data["host"])
        wrong_host_allowed["online_cpus"] = [2, 6, 7]
        wrong_host_allowed["allowed_cpus"] = [2, 6, 7]
        mutations.append((wrong_host_allowed, data["evidence"]))
        wrong_host_authority = copy.deepcopy(data["host"])
        wrong_host_authority["authority"]["machine_id_sha256"] = "8" * 64
        mutations.append((wrong_host_authority, data["evidence"]))
        for index, (host, evidence) in enumerate(mutations):
            with self.subTest(index=index):
                self.assertRejected(lambda host=host, evidence=evidence:
                    execution.validate_qualification_binding(
                        binding, data["preregistration"], host,
                        evidence=evidence))

    def test_armed_record_replays_plan_artifacts_and_qualification(self) -> None:
        data = self.fresh()
        armed = data["armed"]
        observed = execution.validate_armed_record(
            armed, data["preregistration"], data["plan"], data["host"],
            data["artifacts"], data["qualification"],
            qualification_evidence=data["evidence"])
        self.assertEqual(observed, armed)
        self.assertEqual(armed["state"], "ARMED")
        self.assertEqual(armed["evidence_attempt_limit"], 3)
        self.assertEqual(armed["publication_rule"],
                         "atomic-no-replace-before-child/v1")
        self.assertTrue(armed["crash_after_publication_consumes"])
        self.assertEqual(execution.validate_armed_record_shape(armed), armed)

        # Shape-only recovery authenticates syntax and chaining, not sidecars.
        shape_only = copy.deepcopy(armed)
        shape_only["artifact_bundle_sha256"] = "7" * 64
        self.assertEqual(
            execution.validate_armed_record_shape(shape_only), shape_only)
        self.assertRejected(lambda: execution.validate_armed_record(
            shape_only, data["preregistration"], data["plan"], data["host"],
            data["artifacts"], data["qualification"],
            qualification_evidence=data["evidence"]))

        second = execution.armed_record(
            data["preregistration"], data["plan"], data["host"],
            data["artifacts"], data["qualification"],
            qualification_evidence=data["evidence"], evidence_attempt=2,
            prior_armed_chain_sha256=contract.canonical_sha256(armed),
            **armed_closure_kwargs(data["bridge"]))
        self.assertEqual(second["evidence_attempt"], 2)
        self.assertEqual(second["prior_armed_chain_sha256"],
                         contract.canonical_sha256(armed))
        self.assertRejected(lambda: execution.armed_record(
            data["preregistration"], data["plan"], data["host"],
            data["artifacts"], data["qualification"],
            qualification_evidence=data["evidence"], evidence_attempt=2,
            prior_armed_chain_sha256=None,
            **armed_closure_kwargs(data["bridge"])), "prior ARMED")
        self.assertRejected(lambda: execution.armed_record(
            data["preregistration"], data["plan"], data["host"],
            data["artifacts"], data["qualification"],
            qualification_evidence=data["evidence"], evidence_attempt=1,
            prior_armed_chain_sha256="c" * 64,
            **armed_closure_kwargs(data["bridge"])), "first evidence")
        self.assertRejected(lambda: execution.armed_record(
            data["preregistration"], data["plan"], data["host"],
            data["artifacts"], data["qualification"],
            qualification_evidence=data["evidence"], evidence_attempt=4,
            prior_armed_chain_sha256="c" * 64,
            **armed_closure_kwargs(data["bridge"])))

    def test_armed_record_rejects_resigned_or_stale_inputs(self) -> None:
        data = self.fresh()
        mutations = []
        changed_plan = copy.deepcopy(data["plan"])
        changed_plan["child_plans"][0]["argv_tail"][-1] = "not-json-output"
        mutations.append((changed_plan, data["artifacts"],
                          data["qualification"], data["evidence"]))
        changed_artifacts = copy.deepcopy(data["artifacts"])
        changed_artifacts["roles"]["candidate"]["handle_id"] = "d" * 64
        mutations.append((data["plan"], changed_artifacts,
                          data["qualification"], data["evidence"]))
        changed_qualification = copy.deepcopy(data["qualification"])
        changed_qualification["bridge_sha256"] = "e" * 64
        mutations.append((data["plan"], data["artifacts"],
                          changed_qualification, data["evidence"]))
        changed_evidence = copy.deepcopy(data["evidence"])
        verdict = copy.deepcopy(data["verdict"])
        verdict["candidate_timing_performed"] = True
        changed_evidence["independent_verdict_data"] = \
            contract.canonical_json_bytes(verdict)
        mutations.append((data["plan"], data["artifacts"],
                          data["qualification"], changed_evidence))
        for index, (plan, artifacts, qualification, evidence) in \
                enumerate(mutations):
            with self.subTest(index=index):
                self.assertRejected(lambda plan=plan, artifacts=artifacts,
                                           qualification=qualification,
                                           evidence=evidence:
                    execution.armed_record(
                        data["preregistration"], plan, data["host"], artifacts,
                        qualification, qualification_evidence=evidence,
                        evidence_attempt=1, prior_armed_chain_sha256=None,
                        **armed_closure_kwargs(data["bridge"])))

    def test_armed_wire_shape_rejects_claim_and_hash_mutations(self) -> None:
        data = self.fresh()
        armed = data["armed"]
        mutations = []
        extra = copy.deepcopy(armed)
        extra["safe_to_execute"] = True
        mutations.append(extra)
        unarmed = copy.deepcopy(armed)
        unarmed["state"] = "TIMING"
        mutations.append(unarmed)
        nonboolean = copy.deepcopy(armed)
        nonboolean["crash_after_publication_consumes"] = 1
        mutations.append(nonboolean)
        wrong_pair = copy.deepcopy(armed)
        wrong_pair["selected_pair"] = {
            "benchmark_cpu": 6, "reserved_sibling": 2}
        mutations.append(wrong_pair)
        wrong_source = copy.deepcopy(armed)
        wrong_source["candidate_source"]["commit"] = "3" * 40
        mutations.append(wrong_source)
        for index, mutation in enumerate(mutations):
            with self.subTest(index=index):
                self.assertRejected(lambda mutation=mutation:
                    execution.validate_armed_record(
                        mutation, data["preregistration"], data["plan"],
                        data["host"], data["artifacts"],
                        data["qualification"],
                        qualification_evidence=data["evidence"]))

    def test_all_post_import_validation_is_pure(self) -> None:
        data = self.fresh()
        with mock.patch("builtins.open", side_effect=AssertionError(
                "execution contract attempted live file I/O")):
            self.assertEqual(execution.validate_armed_record(
                data["armed"], data["preregistration"], data["plan"],
                data["host"], data["artifacts"], data["qualification"],
                qualification_evidence=data["evidence"]), data["armed"])


if __name__ == "__main__":
    unittest.main()
