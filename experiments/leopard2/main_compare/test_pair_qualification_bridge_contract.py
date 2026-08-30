#!/usr/bin/env python3
"""Deterministic, host-free tests for the gapless bridge contract."""

from __future__ import annotations

import ast
import copy
import importlib.util
from pathlib import Path
import sys
import unittest


HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


fixtures = load_module(
    "leopard2_pair_qualification_acquire_bridge_fixtures",
    HERE / "test_pair_qualification_acquire.py")
bridge = load_module(
    "leopard2_pair_qualification_bridge_contract_test",
    HERE / "pair_qualification_bridge_contract.py")
contract = bridge.contract


def geometry_fixture(
    *, minimum: int = 2, maximum: int = 4,
    nominal_ns: int = 250_000_000, maximum_elapsed_ns: int = 2_000_000_000,
) -> dict:
    return bridge.bridge_geometry_record(
        minimum_window_count=minimum,
        maximum_window_count=maximum,
        nominal_window_ns=nominal_ns,
        maximum_handoff_elapsed_ns=maximum_elapsed_ns)


def bridge_snapshots(
    acquisition: dict, *, window_count: int = 2,
    nonidle_cpu: int | None = None, stalled_cpu: int | None = None,
) -> list[dict]:
    endpoints = [copy.deepcopy(acquisition["scan"]["windows"][-1]["after"])]
    for unused_index in range(window_count):
        before = endpoints[-1]
        counters = {}
        for key, value in before["cpus"].items():
            cpu = int(key)
            fields = dict(value["fields"])
            increment = 1 if cpu == stalled_cpu else 25
            if cpu == nonidle_cpu:
                fields["user"] += 1
                fields["idle"] += increment - 1
            else:
                fields["idle"] += increment
            counters[cpu] = fields
        endpoints.append(contract.shared_snapshot_record(
            read_started_monotonic_ns=
            before["read_finished_monotonic_ns"] + 250_000_000,
            read_finished_monotonic_ns=
            before["read_finished_monotonic_ns"] + 251_000_000,
            counters=counters))
    return endpoints


def bridge_fixture(
    *, window_count: int = 2, nonidle_cpu: int | None = None,
    stalled_cpu: int | None = None,
) -> tuple[dict, dict, dict, dict]:
    acquisition, policy, unused_reader = fixtures.acquisition_fixture()
    geometry = geometry_fixture()
    record = bridge.pair_qualification_bridge_record(
        acquisition, expected_policy=policy, expected_frozen_pair=None,
        expected_acquisition_window_count=2,
        expected_acquisition_nominal_window_ns=250_000_000,
        expected_bridge_geometry=geometry,
        snapshots=bridge_snapshots(
            acquisition, window_count=window_count,
            nonidle_cpu=nonidle_cpu, stalled_cpu=stalled_cpu))
    return record, acquisition, policy, geometry


class PairQualificationBridgeContractTest(unittest.TestCase):
    def assertBridgeFailure(self, function, *args, **kwargs) -> None:  # noqa: N802
        with self.assertRaises(contract.QualificationError):
            function(*args, **kwargs)

    def validate(self, record, acquisition, policy, geometry):
        return bridge.validate_pair_qualification_bridge(
            record, acquisition, expected_policy=policy,
            expected_frozen_pair=None,
            expected_acquisition_window_count=2,
            expected_acquisition_nominal_window_ns=250_000_000,
            expected_bridge_geometry=geometry)

    def test_golden_bridge_is_gapless_canonical_and_conditioned(self) -> None:
        record, acquisition, policy, geometry = bridge_fixture()
        self.assertEqual(self.validate(record, acquisition, policy, geometry), record)
        self.assertEqual(record["bridge_head_sha256"], record["scan_tail_sha256"])
        self.assertEqual(
            record["windows"][0]["before"],
            acquisition["scan"]["windows"][-1]["after"])
        self.assertEqual(
            record["windows"][-1]["after"],
            record["campaign_presample_before"])
        self.assertEqual(record["selected_pair"], {
            "benchmark_cpu": 2, "reserved_sibling": 6,
        })
        self.assertEqual(record["guarded_cpus"], [2, 3, 6, 7])
        self.assertTrue(record["bridge_accepted"])
        self.assertEqual(record["nonidle_guarded_cpus"], [])
        self.assertEqual(record["not_live_guarded_cpus"], [])
        self.assertEqual(record["shared_host_claim_ceiling"], {
            "promotion_eligible": False,
            "host_exclusivity_proved": False,
            "whole_campaign_interval_observed": False,
            "causal_performance_claim_allowed": False,
        })
        self.assertEqual(
            contract.canonical_sha256(record),
            "8bedcecacc32dcfe01396c8d6ccfc394258d9e03960a8300e468601382c98c08")

    def test_spliced_scan_tail_is_rejected_before_bridge_construction(self) -> None:
        acquisition, policy, unused_reader = fixtures.acquisition_fixture()
        snapshots = bridge_snapshots(acquisition)
        snapshots[0]["cpus"]["2"]["fields"]["idle"] += 1
        snapshots[0]["cpus"]["2"]["total_jiffies"] += 1
        self.assertBridgeFailure(
            bridge.pair_qualification_bridge_record,
            acquisition, expected_policy=policy, expected_frozen_pair=None,
            expected_acquisition_window_count=2,
            expected_acquisition_nominal_window_ns=250_000_000,
            expected_bridge_geometry=geometry_fixture(), snapshots=snapshots)

    def test_missing_middle_window_cannot_be_relabelled_as_gapless(self) -> None:
        record, acquisition, policy, geometry = bridge_fixture(window_count=3)
        first = record["windows"][0]
        third = record["windows"][2]
        relabelled_third = contract.pair_observation_window_record(
            policy, phase="bridge", index=1, before=third["before"],
            after=third["after"], observed_cpus=record["observed_cpus"])
        record["windows"] = [first, relabelled_third]
        self.assertBridgeFailure(
            self.validate, record, acquisition, policy, geometry)

    def test_bridge_deadline_and_window_geometry_are_external_bounds(self) -> None:
        acquisition, policy, unused_reader = fixtures.acquisition_fixture()
        snapshots = bridge_snapshots(acquisition)
        too_short = geometry_fixture(maximum_elapsed_ns=500_000_000)
        self.assertBridgeFailure(
            bridge.pair_qualification_bridge_record,
            acquisition, expected_policy=policy, expected_frozen_pair=None,
            expected_acquisition_window_count=2,
            expected_acquisition_nominal_window_ns=250_000_000,
            expected_bridge_geometry=too_short, snapshots=snapshots)
        impossible = {
            "minimum_window_count": 5,
            "maximum_window_count": 5,
            "nominal_window_ns": 250_000_000,
            "maximum_handoff_elapsed_ns": 1_000_000_000,
        }
        self.assertBridgeFailure(bridge.validate_bridge_geometry, impossible)

    def test_policy_pair_geometry_and_acquisition_relabels_are_rejected(self) -> None:
        record, acquisition, policy, geometry = bridge_fixture()
        mutations = []
        for key, value in (
                ("policy_sha256", "0" * 64),
                ("selected_pair", {"benchmark_cpu": 3, "reserved_sibling": 7}),
                ("scan_tail_sha256", "0" * 64),
                ("bridge_geometry_sha256", "0" * 64),
                ("acquisition_sha256", "0" * 64),
                ("candidate_timing_performed", True)):
            mutation = copy.deepcopy(record)
            mutation[key] = value
            mutations.append(mutation)
        claim = copy.deepcopy(record)
        claim["shared_host_claim_ceiling"]["promotion_eligible"] = True
        mutations.append(claim)
        for mutation in mutations:
            with self.subTest(field_differences={
                    key for key in record
                    if not contract.exact_json_equal(record[key], mutation[key])}):
                self.assertBridgeFailure(
                    self.validate, mutation, acquisition, policy, geometry)

        malformed_acquisition = copy.deepcopy(acquisition)
        malformed_acquisition["scan_sha256"] = "0" * 64
        self.assertBridgeFailure(
            self.validate, record, malformed_acquisition, policy, geometry)
        changed_geometry = geometry_fixture(maximum=5)
        self.assertBridgeFailure(
            bridge.validate_pair_qualification_bridge,
            record, acquisition, expected_policy=policy,
            expected_frozen_pair=None,
            expected_acquisition_window_count=2,
            expected_acquisition_nominal_window_ns=250_000_000,
            expected_bridge_geometry=changed_geometry)

    def test_nonidle_and_not_live_bridges_are_valid_rejected_diagnostics(self) -> None:
        nonidle, acquisition, policy, geometry = bridge_fixture(nonidle_cpu=2)
        self.assertEqual(self.validate(nonidle, acquisition, policy, geometry), nonidle)
        self.assertFalse(nonidle["bridge_accepted"])
        self.assertEqual(nonidle["nonidle_guarded_cpus"], [2])

        stalled, acquisition, policy, geometry = bridge_fixture(stalled_cpu=2)
        self.assertEqual(self.validate(stalled, acquisition, policy, geometry), stalled)
        self.assertFalse(stalled["bridge_accepted"])
        self.assertEqual(stalled["not_live_guarded_cpus"], [2])

    def test_producer_acceptance_attestation_cannot_override_recomputation(self) -> None:
        record, acquisition, policy, geometry = bridge_fixture(nonidle_cpu=2)
        self.assertFalse(record["bridge_accepted"])
        record["bridge_accepted"] = True
        record["nonidle_guarded_cpus"] = []
        self.assertBridgeFailure(
            self.validate, record, acquisition, policy, geometry)

    def test_strict_loader_rejects_duplicate_multivalue_and_nonfinite_json(self) -> None:
        record, acquisition, policy, geometry = bridge_fixture()
        canonical = contract.canonical_json_bytes(record)
        self.assertEqual(bridge.load_pair_qualification_bridge(
            canonical, acquisition, expected_policy=policy,
            expected_frozen_pair=None,
            expected_acquisition_window_count=2,
            expected_acquisition_nominal_window_ns=250_000_000,
            expected_bridge_geometry=geometry), record)
        for invalid in (
                b'{"schema":1,"schema":2}\n', canonical + canonical,
                b'{"value":NaN}\n', b'{"value":1e309}\n'):
            with self.subTest(invalid=invalid[:32]):
                self.assertBridgeFailure(
                    bridge.load_pair_qualification_bridge,
                    invalid, acquisition, expected_policy=policy,
                    expected_frozen_pair=None,
                    expected_acquisition_window_count=2,
                    expected_acquisition_nominal_window_ns=250_000_000,
                    expected_bridge_geometry=geometry)

    def test_bridge_module_has_no_host_io_or_process_surface(self) -> None:
        tree = ast.parse(
            (HERE / "pair_qualification_bridge_contract.py").read_text("utf-8"))
        imports = {
            alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
            for alias in node.names
        }
        from_imports = {
            node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
        }
        called_names = {
            node.func.id for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }
        self.assertEqual(imports, {"copy", "importlib.util", "sys"})
        self.assertEqual(from_imports, {"__future__", "pathlib", "typing"})
        self.assertTrue({"open", "sleep", "system", "fork", "exec", "kill"}.
                        isdisjoint(called_names))
        text = (HERE / "pair_qualification_bridge_contract.py").read_text("utf-8")
        self.assertNotIn("sched_setaffinity", text)
        self.assertNotIn("subprocess", text)


if __name__ == "__main__":
    unittest.main()
