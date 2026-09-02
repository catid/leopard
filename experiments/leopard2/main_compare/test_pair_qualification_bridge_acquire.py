#!/usr/bin/env python3
"""Host-free adversarial tests for live v19 bridge acquisition."""

from __future__ import annotations

import ast
import copy
import importlib.util
import os
from pathlib import Path
import sys
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "pair_qualification_bridge_acquire.py"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


fixtures = load_module(
    "leopard2_pair_qualification_v19_bridge_acquire_fixtures",
    HERE / "test_pair_qualification_acquire.py")
live = load_module(
    "leopard2_pair_qualification_v19_bridge_acquire_test", MODULE_PATH)
contract = live.contract


def policy_fixture() -> dict:
    return contract.qualification_policy_record(
        clock_ticks_per_second=100, candidate_primary_cpus=[2, 3],
        excluded_pairs=[], domain_mode="pair-only")


def proc_stat_with_user(snapshot_index: int, cpu: int, user: int) -> bytes:
    data = fixtures.proc_stat(snapshot_index)
    idle = 1000 + cpu + snapshot_index * 25
    original = f"cpu{cpu} 0 0 0 {idle} 0 0 0 0 0 0".encode("ascii")
    changed = f"cpu{cpu} {user} 0 0 {idle - user} 0 0 0 0 0 0".encode(
        "ascii")
    if original not in data:
        raise RuntimeError("test /proc/stat row is missing")
    return data.replace(original, changed)


def bridge_fixture(
    *, proc_records=None, topology_variants=None, allowed_values=None,
    tick_values=None,
) -> tuple[dict, dict, dict, fixtures.FakeHostReader]:
    records = proc_records or [
        fixtures.proc_stat(0), fixtures.proc_stat(1), fixtures.proc_stat(2),
        fixtures.proc_stat(6), fixtures.proc_stat(10),
    ]
    reader = fixtures.FakeHostReader(
        proc_records=records, topology_variants=topology_variants,
        allowed_values=allowed_values, tick_values=tick_values)
    policy = policy_fixture()
    acquisition = live.acquire.acquire_pair_qualification(
        reader, policy_value=policy, window_count=2,
        nominal_window_ns=250_000_000,
        frozen_pair_from_prior_attempt=None)
    return acquisition, policy, live.v19_bridge_geometry_record(), reader


def acquire_bridge(
    acquisition: dict, policy: dict, geometry: dict,
    reader: fixtures.FakeHostReader, *, scan_tail=None,
) -> dict:
    tail = copy.deepcopy(
        acquisition["scan"]["windows"][-1]["after"]
        if scan_tail is None else scan_tail)
    return live.acquire_v19_pair_qualification_bridge(
        reader, acquisition, tail, expected_policy=policy,
        expected_policy_sha256=contract.canonical_sha256(policy),
        expected_frozen_pair=None, expected_acquisition_window_count=2,
        expected_acquisition_nominal_window_ns=250_000_000,
        expected_bridge_geometry=geometry)


class PairQualificationBridgeAcquireTest(unittest.TestCase):
    def assertAcquireFailure(self, function, *args, **kwargs) -> None:  # noqa: N802
        with self.assertRaises(live.BridgeAcquisitionError):
            function(*args, **kwargs)

    def test_accepted_chain_is_exact_fixed_and_independently_replayed(self) -> None:
        acquisition, policy, geometry, reader = bridge_fixture()
        record = acquire_bridge(acquisition, policy, geometry, reader)
        tail = acquisition["scan"]["windows"][-1]["after"]
        self.assertEqual(record["windows"][0]["before"], tail)
        self.assertEqual(record["bridge_head_sha256"],
                         record["scan_tail_sha256"])
        self.assertEqual(record["bridge_window_count"], 2)
        self.assertEqual(record["bridge_geometry"], geometry)
        self.assertEqual(record["guarded_cpus"], [2, 6])
        self.assertTrue(record["bridge_accepted"])
        self.assertEqual(reader.sleep_calls, [
            250_000_000, 250_000_000, 1_000_000_000, 1_000_000_000,
        ])
        self.assertEqual(reader.allowed_calls, 4)
        self.assertEqual(reader.tick_calls, 4)
        self.assertEqual(reader.proc_calls, 5)
        verdict = live.verifier.require_accepted_pair_qualification_bundle(
            contract.canonical_json_bytes(acquisition),
            contract.canonical_json_bytes(record), expected_policy=policy,
            expected_policy_sha256=contract.canonical_sha256(policy),
            expected_frozen_pair=None,
            expected_acquisition_window_count=2,
            expected_acquisition_nominal_window_ns=250_000_000,
            expected_bridge_geometry=geometry)
        self.assertEqual(verdict["bridge_sha256"],
                         contract.canonical_sha256(record))
        self.assertEqual(record["shared_host_claim_ceiling"], {
            "promotion_eligible": False,
            "host_exclusivity_proved": False,
            "whole_campaign_interval_observed": False,
            "causal_performance_claim_allowed": False,
        })

    def test_nonidle_and_not_live_bridges_fail_independent_acceptance(self) -> None:
        variants = (
            [
                fixtures.proc_stat(0), fixtures.proc_stat(1),
                fixtures.proc_stat(2), proc_stat_with_user(6, 2, 1),
                proc_stat_with_user(10, 2, 1),
            ],
            [
                fixtures.proc_stat(0), fixtures.proc_stat(1),
                fixtures.proc_stat(2), fixtures.proc_stat(2),
                fixtures.proc_stat(2),
            ],
        )
        for records in variants:
            with self.subTest(last=records[-1]):
                acquisition, policy, geometry, reader = bridge_fixture(
                    proc_records=records)
                self.assertAcquireFailure(
                    acquire_bridge, acquisition, policy, geometry, reader)

    def test_short_window_and_deadline_overrun_fail_closed(self) -> None:
        for scale in (0, 3):
            acquisition, policy, geometry, reader = bridge_fixture()
            reader.sleep_scale = scale
            with self.subTest(scale=scale):
                self.assertAcquireFailure(
                    acquire_bridge, acquisition, policy, geometry, reader)

        acquisition, policy, geometry, reader = bridge_fixture()
        reader.sleep_scale = 4
        self.assertAcquireFailure(
            acquire_bridge, acquisition, policy, geometry, reader)
        self.assertEqual(reader.sleep_calls, [
            250_000_000, 250_000_000, 1_000_000_000,
        ])

    def test_monotonic_domain_and_initial_deadline_fail_before_sleep(self) -> None:
        for now in (0, None):
            acquisition, policy, geometry, reader = bridge_fixture()
            if now is None:
                tail = acquisition["scan"]["windows"][-1]["after"]
                reader.now = tail["read_finished_monotonic_ns"] + 4_000_000_001
            else:
                reader.now = now
            with self.subTest(now=reader.now):
                self.assertAcquireFailure(
                    acquire_bridge, acquisition, policy, geometry, reader)
                self.assertEqual(reader.sleep_calls, [
                    250_000_000, 250_000_000,
                ])

    def test_start_and_finish_host_projection_drift_fail_closed(self) -> None:
        base = fixtures.topology_files()
        changed = fixtures.topology_files(core_3=9)
        allowed = [2, 3, 6, 7]

        acquisition, policy, geometry, reader = bridge_fixture(
            topology_variants=[base, base, changed])
        self.assertAcquireFailure(
            acquire_bridge, acquisition, policy, geometry, reader)
        self.assertEqual(reader.proc_calls, 3)

        acquisition, policy, geometry, reader = bridge_fixture(
            allowed_values=[allowed, allowed, [2, 6, 7]])
        self.assertAcquireFailure(
            acquire_bridge, acquisition, policy, geometry, reader)
        self.assertEqual(reader.proc_calls, 3)

        acquisition, policy, geometry, reader = bridge_fixture(
            topology_variants=[base, base, base, changed])
        self.assertAcquireFailure(
            acquire_bridge, acquisition, policy, geometry, reader)

        acquisition, policy, geometry, reader = bridge_fixture(
            allowed_values=[allowed, allowed, allowed, [2, 6, 7]])
        self.assertAcquireFailure(
            acquire_bridge, acquisition, policy, geometry, reader)

        for tick_values in ([100, 100, 250], [100, 100, 100, 250]):
            acquisition, policy, geometry, reader = bridge_fixture(
                tick_values=tick_values)
            with self.subTest(tick_values=tick_values):
                self.assertAcquireFailure(
                    acquire_bridge, acquisition, policy, geometry, reader)

    def test_wrong_scan_tail_splice_is_rejected_before_host_access(self) -> None:
        acquisition, policy, geometry, reader = bridge_fixture()
        tail = copy.deepcopy(acquisition["scan"]["windows"][-1]["after"])
        tail["cpus"]["2"]["fields"]["idle"] += 1
        tail["cpus"]["2"]["total_jiffies"] += 1
        allowed_calls = reader.allowed_calls
        self.assertAcquireFailure(
            acquire_bridge, acquisition, policy, geometry, reader,
            scan_tail=tail)
        self.assertEqual(reader.allowed_calls, allowed_calls)

    def test_reader_failure_is_translated_but_interrupt_propagates(self) -> None:
        acquisition, policy, geometry, reader = bridge_fixture(
            proc_records=[
                fixtures.proc_stat(0), fixtures.proc_stat(1),
                fixtures.proc_stat(2),
            ])
        self.assertAcquireFailure(
            acquire_bridge, acquisition, policy, geometry, reader)

        acquisition, policy, geometry, reader = bridge_fixture()

        def interrupted(unused_duration_ns: int) -> None:
            raise KeyboardInterrupt("injected interruption")

        reader.sleep_ns = interrupted
        with self.assertRaises(KeyboardInterrupt):
            acquire_bridge(acquisition, policy, geometry, reader)

    def test_geometry_and_pair_only_policy_are_closed(self) -> None:
        geometry = live.v19_bridge_geometry_record()
        self.assertEqual(live.validate_v19_bridge_geometry(geometry), geometry)
        for key, value in (
                ("minimum_window_count", 1),
                ("maximum_window_count", 3),
                ("nominal_window_ns", 999_999_999),
                ("maximum_handoff_elapsed_ns", 5_000_000_001)):
            mutation = dict(geometry)
            mutation[key] = value
            with self.subTest(key=key):
                self.assertAcquireFailure(
                    live.validate_v19_bridge_geometry, mutation)

        acquisition, unused_policy, geometry, reader = bridge_fixture()
        changed_policy = contract.qualification_policy_record(
            clock_ticks_per_second=100, candidate_primary_cpus=[2, 3],
            excluded_pairs=[], domain_mode="pair-and-domain")
        self.assertAcquireFailure(
            live.acquire_v19_pair_qualification_bridge,
            reader, acquisition,
            acquisition["scan"]["windows"][-1]["after"],
            expected_policy=changed_policy,
            expected_policy_sha256=contract.canonical_sha256(changed_policy),
            expected_frozen_pair=None,
            expected_acquisition_window_count=2,
            expected_acquisition_nominal_window_ns=250_000_000,
            expected_bridge_geometry=geometry)

        acquisition, policy, geometry, reader = bridge_fixture()
        allowed_calls = reader.allowed_calls
        self.assertAcquireFailure(
            live.acquire_v19_pair_qualification_bridge,
            reader, acquisition,
            acquisition["scan"]["windows"][-1]["after"],
            expected_policy=policy, expected_policy_sha256="0" * 64,
            expected_frozen_pair=None,
            expected_acquisition_window_count=2,
            expected_acquisition_nominal_window_ns=250_000_000,
            expected_bridge_geometry=geometry)
        self.assertEqual(reader.allowed_calls, allowed_calls)

    def test_fake_path_has_no_mutator_launcher_writer_or_runner_import(self) -> None:
        acquisition, policy, geometry, reader = bridge_fixture()
        with mock.patch.object(os, "system", side_effect=AssertionError), \
                mock.patch.object(os, "fork", side_effect=AssertionError), \
                mock.patch.object(os, "execve", side_effect=AssertionError), \
                mock.patch.object(
                    os, "sched_setaffinity", side_effect=AssertionError,
                    create=True), \
                mock.patch.object(os, "kill", side_effect=AssertionError), \
                mock.patch.object(
                    os, "setpriority", side_effect=AssertionError,
                    create=True), \
                mock.patch.object(os, "open", side_effect=AssertionError), \
                mock.patch.object(os, "write", side_effect=AssertionError), \
                mock.patch.object(os, "unlink", side_effect=AssertionError), \
                mock.patch.object(os, "fsync", side_effect=AssertionError):
            self.assertTrue(
                acquire_bridge(acquisition, policy, geometry, reader)
                ["bridge_accepted"])

        tree = ast.parse(MODULE_PATH.read_text("utf-8"))
        imports = {
            alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
            for alias in node.names
        }
        from_imports = {
            node.module for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        }
        self.assertEqual(imports, {"copy", "importlib.util", "sys"})
        self.assertEqual(from_imports, {"__future__", "pathlib", "typing"})
        text = MODULE_PATH.read_text("utf-8")
        for forbidden in (
                "run_abba", "run_k65", "subprocess", "sched_setaffinity",
                "os.write", "os.open", "fork(", "execve(",
                "write_acquisition_record_exclusive", "SystemHostReader"):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, text)


if __name__ == "__main__":
    unittest.main()
