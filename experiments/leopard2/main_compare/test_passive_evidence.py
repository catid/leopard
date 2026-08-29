#!/usr/bin/python3

"""Deterministic tests for the passive shared-host evidence boundary."""

from __future__ import annotations

import copy
import importlib.util
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[3]
WRAPPER = ROOT / "experiments/leopard2/main_compare" / \
    "run_authoritative_v17_gfni_main_compare.sh"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


audit = load_module(
    "passive_evidence_auditor",
    ROOT / "experiments/leopard2/main_compare" /
    "audit_v17_gfni_main_compare.py",
)
census = load_module(
    "passive_environment_census",
    ROOT / "experiments/leopard2/main_compare" /
    "passive_environment_census.py",
)
runner_fixtures = load_module(
    "passive_evidence_runner_fixtures",
    ROOT / "experiments/leopard2/main_compare" / "test_run_abba.py",
)


def cpu_record(cpu: int, *, user: int = 0, idle: int = 0,
               system: int = 0, iowait: int = 0) -> dict:
    fields = {name: 0 for name in census.CPU_FIELDS}
    fields.update(user=user, idle=idle, system=system, iowait=iowait)
    return {
        "cpu": cpu,
        "fields": fields,
        "total_jiffies": sum(fields[name] for name in census.CPU_FIELDS[:8]),
        "nonidle_jiffies": sum(fields[name]
                                for name in census.NONIDLE_FIELDS),
    }


def raw_cpu_record(cpu: int, *, user: int = 0, idle: int = 0,
                   system: int = 0) -> dict:
    fields = {name: 0 for name in census.CPU_FIELDS[:8]}
    fields.update(user=user, idle=idle, system=system)
    return {"cpu": cpu, "fields": fields,
            "total_jiffies": sum(fields.values())}


def controller() -> dict:
    return {
        "schema": "leopard2-v17-passive-controller-affinity/v1",
        "acquisition_generation": "passive-v1",
        "wrapper_pid": 100,
        "before_allowed_cpus": [0, 1, 52, 116],
        "after_allowed_cpus": [0, 1],
        "runner_launch_allowed_cpus": [0, 1, 52, 116],
        "benchmark_cpu": 52,
        "reserved_sibling": 116,
        "affinity_mutation_scope":
            "wrapper-process-and-owned-descendants-only",
        "active_affinity_supervisor_executed": False,
    }


def v18_controller() -> dict:
    return {
        **controller(),
        "schema": "leopard2-v18-passive-controller-affinity/v1",
        "acquisition_generation": "passive-v2",
        "before_allowed_cpus": [0, 52, 116],
        "after_allowed_cpus": [0],
        "runner_launch_allowed_cpus": [0, 52, 116],
    }


def census_cpu_from_raw(value: dict, *, extra_user: int = 0) -> dict:
    fields = {name: 0 for name in census.CPU_FIELDS}
    fields.update(value["fields"])
    fields["user"] += extra_user
    return {
        "cpu": value["cpu"],
        "fields": fields,
        "total_jiffies": sum(fields[name] for name in census.CPU_FIELDS[:8]),
        "nonidle_jiffies": sum(fields[name] for name in census.NONIDLE_FIELDS),
    }


def v18_snapshot(raw: dict, phase: str, *, cpu52_extra: int = 0,
                 cpu116_extra: int = 0,
                 include_pair_eligible: bool = True) -> dict:
    value = snapshot(phase)
    isolation = raw["isolation"]
    if phase == "pre":
        outer = isolation["before"]
        boundary = outer["monotonic_ns"] - 1
        value["scan_started_monotonic_ns"] = boundary - 2
        value["scan_finished_monotonic_ns"] = boundary - 1
    else:
        outer = isolation["after"]
        boundary = outer["monotonic_ns"] + 1
        value["scan_started_monotonic_ns"] = boundary + 1
        value["scan_finished_monotonic_ns"] = boundary + 2
    value["activity_boundary_monotonic_ns"] = boundary
    value["proc_stat"] = {
        "52": census_cpu_from_raw(
            outer["benchmark_cpu"], extra_user=cpu52_extra),
        "116": census_cpu_from_raw(
            outer["reserved_sibling"], extra_user=cpu116_extra),
    }
    wrapper = value["same_uid_processes"]["entries"][0]
    wrapper["cpus_allowed"] = [0]
    wrapper["reserved_pair_intersection"] = []
    if include_pair_eligible:
        pair_eligible = copy.deepcopy(wrapper)
        pair_eligible.update({
            "pid": 101,
            "tid": 101,
            "process_identity": {
                "device": 1, "inode": 5, "starttime_ticks": 4},
            "thread_identity": {
                "device": 1, "inode": 6, "starttime_ticks": 4},
            "comm_sha256": "b" * 64,
            "cpus_allowed": [0, 52],
            "reserved_pair_intersection": [52],
        })
        value["same_uid_processes"]["entries"].append(pair_eligible)
        value["same_uid_processes"]["summary"].update({
            "retained_process_count": 2,
            "retained_thread_count": 2,
            "pair_eligible_process_count": 1,
            "pair_eligible_thread_count": 1,
            "confined_to_reserved_subset_process_count": 0,
            "confined_to_reserved_subset_thread_count": 0,
        })
    reseal_snapshot(value)
    return value


def raw_evidence() -> dict:
    return {
        "schema": census.RAW_SCHEMA_V17,
        "supervision": None,
        "campaign": {"allowed_cpu_set_at_launch": [0, 1, 52, 116]},
        "invocations": [
            {"duration_ns": 1}
            for _ in range(census.EXPECTED_INVOCATION_COUNT)
        ],
        "isolation": {
            "schema": "leopard2-main-compare-isolation/v1",
            "accepted": True,
            "benchmark_cpu": 52,
            "reserved_sibling": 116,
            "before": {
                "monotonic_ns": 4,
                "benchmark_cpu": raw_cpu_record(52),
                "reserved_sibling": raw_cpu_record(116),
            },
            "after": {
                "monotonic_ns": 16,
                "benchmark_cpu": raw_cpu_record(52, user=1),
                "reserved_sibling": raw_cpu_record(116, idle=1),
            },
            "delta": {}, "pair_lease": {}, "policy": {},
        },
    }


def reseal_snapshot(value: dict) -> None:
    payload = dict(value)
    payload.pop("digest", None)
    value["digest"] = census.digest_payload(payload)


def snapshot(phase: str) -> dict:
    entry = {
        "pid": 100,
        "tid": 100,
        "uid": 1000,
        "process_identity": {
            "device": 1, "inode": 2, "starttime_ticks": 3},
        "thread_identity": {
            "device": 1, "inode": 4, "starttime_ticks": 3},
        "comm_size": 4,
        "comm_sha256": "a" * 64,
        "cpus_allowed": [0, 1],
        "reserved_pair_intersection": [],
    }
    if phase == "pre":
        started, finished, boundary = 1, 2, 3
        proc_stat = {
            "52": cpu_record(52), "116": cpu_record(116)}
    else:
        boundary, started, finished = 20, 21, 22
        proc_stat = {
            "52": cpu_record(52, user=1),
            "116": cpu_record(116, idle=1),
        }
    payload = {
        "schema": census.SCHEMA,
        "phase": phase,
        "semantics": census.SEMANTICS,
        "mutation_policy": census.MUTATION_POLICY,
        "uid": 1000,
        "boot_id": "01234567-89ab-cdef-0123-456789abcdef",
        "namespaces": {
            "pid": {"device": 1, "inode": 10},
            "user": {"device": 1, "inode": 11},
        },
        "clock_ticks_per_second": 100,
        "reserved_cpus": [52, 116],
        "scan_started_monotonic_ns": started,
        "scan_finished_monotonic_ns": finished,
        "activity_boundary_monotonic_ns": boundary,
        "collector": {
            "pid": 200, "allowed_cpus": [0, 1],
            "reserved_pair_excluded": True,
        },
        "proc_stat": proc_stat,
        "same_uid_processes": {
            "entries": [entry],
            "vanished_during_scan": [],
            "summary": {
                "retained_process_count": 1,
                "retained_thread_count": 1,
                "pair_eligible_process_count": 0,
                "pair_eligible_thread_count": 0,
                "nonuniform_process_count": 0,
                "confined_to_reserved_subset_process_count": 0,
                "confined_to_reserved_subset_thread_count": 0,
                "vanished_record_count": 0,
            },
        },
        "capabilities": {
            "atomic_snapshot": False,
            "interval_complete": False,
            "records_execution_history": False,
            "establishes_cpu_exclusivity": False,
        },
    }
    return {**payload, "digest": census.digest_payload(payload)}


class PassiveEvidenceTests(unittest.TestCase):
    def test_passive_policy_is_permanently_nonpromotion(self) -> None:
        policy = census.compare(
            snapshot("pre"), snapshot("post"), raw_evidence(), controller())
        self.assertFalse(policy["cpu_pair_exclusive"])
        self.assertFalse(policy["interval_complete_task_observation"])
        self.assertFalse(policy["foreign_cpu52_work_attributable"])
        self.assertFalse(policy["promotion_eligible"])
        self.assertFalse(policy["promotion_passed"])

    def test_v18_policy_gates_only_retained_windows_and_discloses_endpoints(
            self) -> None:
        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V18)
        pre = v18_snapshot(raw, "pre")
        post = v18_snapshot(
            raw, "post", cpu52_extra=5, cpu116_extra=7)
        policy = census.compare(pre, post, raw, v18_controller())
        self.assertEqual(policy["schema"], census.POLICY_SCHEMA_V2)
        self.assertEqual(policy["acquisition_generation"], "passive-v2")
        self.assertEqual(
            policy["windowed_contamination"]["schema"],
            "leopard2-v18-windowed-census-policy-observation/v1")
        self.assertTrue(policy["windowed_contamination"][
            "all_benchmark_cpu_excess_zero"])
        self.assertTrue(policy["windowed_contamination"][
            "all_reserved_sibling_nonidle_zero"])
        self.assertFalse(policy["outer_disclosure"]["gated"])
        self.assertGreater(
            policy["outer_disclosure"]["outer_cpu_activity"]["52"]
            ["nonidle_jiffies"], 5)
        self.assertGreaterEqual(
            policy["outer_disclosure"]["outer_cpu_activity"]["116"]
            ["nonidle_jiffies"], 7)
        self.assertEqual(policy["shared_host_exposure"]["pre"]
                         ["confined_to_reserved_subset_thread_count"], 0)
        self.assertEqual(policy["shared_host_exposure"]["post"]
                         ["confined_to_reserved_subset_thread_count"], 0)
        self.assertEqual(policy["shared_host_exposure"]["pre"]
                         ["pair_eligible_thread_count"], 1)
        self.assertEqual(policy["shared_host_exposure"]["post"]
                         ["pair_eligible_thread_count"], 1)
        self.assertFalse(policy["promotion_eligible"])
        self.assertFalse(policy["promotion_passed"])
        self.assertFalse(policy["causal_performance_claim_eligible"])
        self.assertFalse(policy["cpu_pair_exclusive"])

    def test_v18_policy_rejects_confined_thread_and_persistent_affinity_change(
            self) -> None:
        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V18)
        pre = v18_snapshot(raw, "pre")
        post = v18_snapshot(raw, "post")
        for value in (pre, post):
            entry = value["same_uid_processes"]["entries"][1]
            entry["cpus_allowed"] = [52]
            entry["reserved_pair_intersection"] = [52]
            value["same_uid_processes"]["summary"].update({
                "confined_to_reserved_subset_process_count": 1,
                "confined_to_reserved_subset_thread_count": 1,
            })
            reseal_snapshot(value)
        with self.assertRaises(census.CensusError):
            census.compare(pre, post, raw, v18_controller())

        pre = v18_snapshot(raw, "pre")
        post = v18_snapshot(raw, "post")
        post_entry = post["same_uid_processes"]["entries"][1]
        post_entry["cpus_allowed"] = [0, 116]
        post_entry["reserved_pair_intersection"] = [116]
        reseal_snapshot(post)
        with self.assertRaises(census.CensusError):
            census.compare(pre, post, raw, v18_controller())

    def test_v18_census_rejects_bool_as_integer_in_window_closure(self) -> None:
        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V18)
        raw["invocations"][0]["cpu_window"]["delta"]["benchmark_cpu"] \
            ["fields"]["nice"] = False
        raw["isolation"]["invocation_windows"][0]["delta"] \
            ["benchmark_cpu"]["fields"]["nice"] = False
        with self.assertRaises(census.CensusError):
            census.compare(
                v18_snapshot(raw, "pre", include_pair_eligible=False),
                v18_snapshot(raw, "post", include_pair_eligible=False),
                raw, v18_controller())

    def test_v18_counter_epochs_must_fit_campaign_endpoints_everywhere(
            self) -> None:
        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V18)
        for invocation, retained in zip(
                raw["invocations"],
                raw["isolation"]["invocation_windows"]):
            for window in (invocation["cpu_window"], retained):
                for phase in ("before", "after"):
                    for endpoint in ("benchmark_cpu", "reserved_sibling"):
                        window[phase][endpoint]["fields"]["user"] += 1000
                        window[phase][endpoint]["total_jiffies"] += 1000
        raw = runner_fixtures.resign(raw)
        windows = [invocation["cpu_window"]
                   for invocation in raw["invocations"]]
        with self.assertRaises(runner_fixtures.runner.EvidenceError):
            runner_fixtures.runner.validate_raw(
                raw, None, check_files=False, check_current_inputs=False)
        with self.assertRaises(audit.AuditError):
            audit.validate_isolation_v2(raw["isolation"], windows)
        with self.assertRaises(census.CensusError):
            census.compare(
                v18_snapshot(raw, "pre", include_pair_eligible=False),
                v18_snapshot(raw, "post", include_pair_eligible=False),
                raw, v18_controller())

    def test_v18_zero_thresholds_do_not_reuse_legacy_allowance(self) -> None:
        runner_policy = runner_fixtures.runner.cpu_window_policy()
        audit_policy = audit.cpu_window_policy()
        census_policy = census.window_policy()
        for policy in (runner_policy, audit_policy, census_policy):
            self.assertEqual(
                policy["benchmark_cpu_max_nonidle_excess_jiffies"], 0)
            self.assertEqual(
                policy["reserved_sibling_max_nonidle_jiffies"], 0)
        self.assertEqual(
            census.
            LEGACY_V17_CPU52_AGGREGATE_EXCESS_LIMIT_JIFFIES,
            16)
        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V18)
        arguments = (
            v18_snapshot(raw, "pre", include_pair_eligible=False),
            v18_snapshot(raw, "post", include_pair_eligible=False),
            raw, v18_controller())
        expected = census.compare(*arguments)
        with mock.patch.object(
                census,
                "LEGACY_V17_CPU52_AGGREGATE_EXCESS_LIMIT_JIFFIES",
                10_000):
            self.assertEqual(census.compare(*arguments), expected)
        forbidden_legacy_literal = (
            "MAX_NONIDLE_EXCESS_OVER_CHILD_WALL_CEILING_JIFFIES")
        for path in (
                ROOT / "experiments/leopard2/main_compare" /
                "passive_environment_census.py",
                ROOT / "experiments/leopard2/main_compare" /
                "audit_v17_gfni_main_compare.py"):
            with self.subTest(path=path):
                self.assertNotIn(
                    forbidden_legacy_literal,
                    path.read_text(encoding="utf-8"))

    def test_v18_auditor_self_test_covers_window_boundaries(self) -> None:
        completed = subprocess.run(
            [sys.executable, "-I", "-S", "-B", str(
                ROOT / "experiments/leopard2/main_compare" /
                "audit_v17_gfni_main_compare.py"), "--self-test"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(
            completed.stdout,
            "v17 GFNI exact-main independent auditor self-test passed\n")

    def test_v18_audit_claim_preserves_planned_window_fields(self) -> None:
        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V18)
        observation = audit.windowed_observation(raw["isolation"])
        fields = audit.audit_mode_result_fields("windowed", observation)
        self.assertEqual(
            [item["index"] for item in observation["per_invocation"]],
            list(range(12)))
        claim = fields["isolation_claim"]
        self.assertEqual(
            claim["windowed_screen"],
            "per retained benchmark invocation")
        self.assertIs(claim["out_of_window_activity_gated"], False)
        self.assertIs(
            claim[
                "reserved_sibling_zero_nonidle_jiffies_in_every_retained_window"],
            True)

    def test_persistent_affinity_change_rejects(self) -> None:
        pre, post = snapshot("pre"), snapshot("post")
        post["same_uid_processes"]["entries"][0]["cpus_allowed"] = [0]
        post["same_uid_processes"]["entries"][0][
            "reserved_pair_intersection"] = []
        post["same_uid_processes"]["summary"].update(
            pair_eligible_process_count=0, pair_eligible_thread_count=0)
        reseal_snapshot(post)
        with self.assertRaises(census.CensusError):
            census.compare(pre, post, raw_evidence(), controller())

    def test_persistent_nonwrapper_affinity_change_rejects(self) -> None:
        pre, post = snapshot("pre"), snapshot("post")
        for value in (pre, post):
            entry = copy.deepcopy(
                value["same_uid_processes"]["entries"][0])
            entry.update(
                pid=101,
                tid=101,
                process_identity={
                    "device": 1, "inode": 5, "starttime_ticks": 4},
                thread_identity={
                    "device": 1, "inode": 6, "starttime_ticks": 4},
                comm_sha256="b" * 64,
            )
            value["same_uid_processes"]["entries"].append(entry)
            value["same_uid_processes"]["summary"].update(
                retained_process_count=2, retained_thread_count=2)
        post["same_uid_processes"]["entries"][1]["cpus_allowed"] = [0]
        reseal_snapshot(pre)
        reseal_snapshot(post)
        with self.assertRaises(census.CensusError):
            census.compare(pre, post, raw_evidence(), controller())

    def test_outer_sibling_activity_rejects(self) -> None:
        post = snapshot("post")
        post["proc_stat"]["116"] = cpu_record(116, idle=1, system=1)
        reseal_snapshot(post)
        with self.assertRaises(census.CensusError):
            census.compare(
                snapshot("pre"), post, raw_evidence(), controller())

    def test_raw_counters_must_fit_outer_interval(self) -> None:
        raw = raw_evidence()
        raw["isolation"]["before"]["benchmark_cpu"]["fields"]["user"] = 2
        raw["isolation"]["before"]["benchmark_cpu"]["total_jiffies"] = 2
        raw["isolation"]["after"]["benchmark_cpu"]["fields"]["user"] = 3
        raw["isolation"]["after"]["benchmark_cpu"]["total_jiffies"] = 3
        with self.assertRaises(census.CensusError):
            census.compare(snapshot("pre"), snapshot("post"), raw, controller())

    def test_outer_benchmark_cpu_excess_boundary(self) -> None:
        post = snapshot("post")
        post["proc_stat"]["52"] = cpu_record(52, user=17)
        reseal_snapshot(post)
        policy = census.compare(
            snapshot("pre"), post, raw_evidence(), controller())
        self.assertEqual(
            policy["outer_contamination"]
            ["benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies"],
            16)

        post["proc_stat"]["52"] = cpu_record(52, user=18)
        reseal_snapshot(post)
        with self.assertRaises(census.CensusError):
            census.compare(
                snapshot("pre"), post, raw_evidence(), controller())

    def test_outer_auxiliary_classes_are_disclosed_not_zero_gated(self) -> None:
        post = snapshot("post")
        post["proc_stat"]["52"] = cpu_record(52, user=1, iowait=3)
        reseal_snapshot(post)
        policy = census.compare(
            snapshot("pre"), post, raw_evidence(), controller())
        self.assertEqual(
            policy["outer_contamination"]
            ["benchmark_cpu_auxiliary_class_jiffies"]["iowait"], 3)
        self.assertIn(
            "no auxiliary class has a separate zero-tolerance gate",
            policy["outer_contamination"]["auxiliary_class_policy"])

    def test_single_confined_thread_in_nonuniform_process_rejects(self) -> None:
        value = snapshot("pre")
        second = copy.deepcopy(value["same_uid_processes"]["entries"][0])
        second["tid"] = 101
        second["thread_identity"] = {
            "device": 1, "inode": 5, "starttime_ticks": 3}
        second["cpus_allowed"] = [52]
        second["reserved_pair_intersection"] = [52]
        value["same_uid_processes"]["entries"].append(second)
        value["same_uid_processes"]["summary"].update(
            retained_thread_count=2,
            pair_eligible_process_count=1,
            pair_eligible_thread_count=1,
            nonuniform_process_count=1,
            confined_to_reserved_subset_process_count=1,
            confined_to_reserved_subset_thread_count=1,
        )
        reseal_snapshot(value)
        with self.assertRaises(census.CensusError):
            census.validate_snapshot(value, "pre")

    def test_summary_boolean_is_not_an_integer(self) -> None:
        value = snapshot("pre")
        value["same_uid_processes"]["summary"][
            "confined_to_reserved_subset_thread_count"] = False
        reseal_snapshot(value)
        with self.assertRaises(census.CensusError):
            census.validate_snapshot(value, "pre")

    def test_snapshot_numeric_coercions_reject(self) -> None:
        mutations = (
            lambda value: value.__setitem__(
                "clock_ticks_per_second", 100.0),
            lambda value: value["capabilities"].__setitem__(
                "atomic_snapshot", 0),
            lambda value: value["same_uid_processes"]["entries"][0]
                .__setitem__("uid", 1000.0),
            lambda value: value.__setitem__(
                "reserved_cpus", [52.0, 116.0]),
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                value = snapshot("pre")
                mutation(value)
                reseal_snapshot(value)
                with self.assertRaises(census.CensusError):
                    census.validate_snapshot(value, "pre")

        value = snapshot("pre")
        entry = value["same_uid_processes"]["entries"][0]
        entry["cpus_allowed"] = [0, 1, 52]
        entry["reserved_pair_intersection"] = [52.0]
        value["same_uid_processes"]["summary"].update(
            pair_eligible_process_count=1, pair_eligible_thread_count=1)
        reseal_snapshot(value)
        with self.assertRaises(census.CensusError):
            census.validate_snapshot(value, "pre")

        controller_value = controller()
        controller_value["runner_launch_allowed_cpus"] = [
            0.0, 1.0, 52.0, 116.0]
        with self.assertRaises(census.CensusError):
            census.validate_controller(
                controller_value, census.RAW_SCHEMA_V17)

    def test_namespace_identity_follows_the_target(self) -> None:
        with tempfile.TemporaryDirectory(prefix="leopard-passive-ns-test.") \
                as directory:
            root = Path(directory)
            target = root / "target"
            target.write_bytes(b"namespace fixture")
            link = root / "namespace"
            link.symlink_to(target)
            identity = census.namespace_identity(link)
            metadata = target.stat()
            self.assertEqual(identity, {
                "device": metadata.st_dev, "inode": metadata.st_ino})

    def test_controller_launch_mask_is_cross_bound(self) -> None:
        raw = raw_evidence()
        raw["campaign"]["allowed_cpu_set_at_launch"] = [0, 1, 52]
        with self.assertRaises(census.CensusError):
            census.compare(
                snapshot("pre"), snapshot("post"), raw, controller())

        raw = raw_evidence()
        raw["campaign"]["allowed_cpu_set_at_launch"] = [
            0.0, 1.0, 52.0, 116.0]
        with self.assertRaises(census.CensusError):
            census.compare(
                snapshot("pre"), snapshot("post"), raw, controller())

    def test_null_supervision_must_be_present(self) -> None:
        raw = raw_evidence()
        del raw["supervision"]
        with self.assertRaises(census.CensusError):
            census.compare(
                snapshot("pre"), snapshot("post"), raw, controller())

    def test_child_wall_durations_must_fit_isolation_window(self) -> None:
        raw = raw_evidence()
        raw["invocations"][0]["duration_ns"] = 2
        with self.assertRaises(census.CensusError):
            census.compare(
                snapshot("pre"), snapshot("post"), raw, controller())

    def test_full_proc_identity_device_is_persistent(self) -> None:
        post = snapshot("post")
        post["same_uid_processes"]["entries"][0]["process_identity"][
            "device"] = 9
        reseal_snapshot(post)
        with self.assertRaises(census.CensusError):
            census.compare(
                snapshot("pre"), post, raw_evidence(), controller())

    def test_auditor_modes_are_anti_relabeling(self) -> None:
        campaign, reservation, isolation = {}, {}, {}
        self.assertIsNone(audit.supervision_for_mode(
            None, None, campaign, reservation, isolation, "absent"))
        with self.assertRaises(audit.AuditError):
            audit.supervision_for_mode(
                {}, None, campaign, reservation, isolation, "absent")
        with self.assertRaises(audit.AuditError):
            audit.supervision_for_mode(
                None, None, campaign, reservation, isolation, "required")

    def test_passive_foreign_jiffy_boundary(self) -> None:
        isolation = {
            "before": {"monotonic_ns": 1},
            "after": {"monotonic_ns": 2},
            "delta": {
                "benchmark_cpu": {"nonidle_jiffies": 17, "idle_jiffies": 0},
                "reserved_sibling": {"nonidle_jiffies": 0},
            },
        }
        result = audit.passive_contamination(isolation, 1)
        self.assertEqual(
            result[
                "benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies"],
            16)
        over = copy.deepcopy(isolation)
        over["delta"]["benchmark_cpu"]["nonidle_jiffies"] = 18
        with self.assertRaises(audit.AuditError):
            audit.passive_contamination(over, 1)

    def test_passive_auditor_contract_is_permanently_nonpromotion(self) -> None:
        isolation = {
            "before": {"monotonic_ns": 1},
            "after": {"monotonic_ns": 2},
            "delta": {
                "benchmark_cpu": {
                    "nonidle_jiffies": 1, "idle_jiffies": 0},
                "reserved_sibling": {"nonidle_jiffies": 0},
            },
        }
        contamination = audit.passive_contamination(isolation, 1)
        fields = audit.audit_mode_result_fields("absent", contamination)
        self.assertEqual(fields["schema"], audit.PASSIVE_AUDIT_SCHEMA)
        self.assertEqual(
            fields["audit_mode"], "passive-shared-host-observation")
        self.assertEqual(
            fields["evidence_class"],
            "passive-shared-host-observation/v1")
        for name in (
                "affinity_supervisor_binding_verified", "promotion_eligible",
                "promotion_passed", "causal_performance_claim_eligible",
                "cpu_pair_exclusive"):
            self.assertIs(fields[name], False)
        claim = fields["isolation_claim"]
        self.assertIsNone(claim["campaign_supervision"])
        for name in (
                "foreign_process_affinity_mutation_claimed",
                "foreign_process_signalling_claimed",
                "benchmark_cpu_exclusive_ownership_claimed",
                "same_uid_pair_exclusion_certificate",
                "cgroup_or_os_exclusive_certificate",
                "interval_complete_task_observation",
                "benchmark_cpu_foreign_work_attributable",
                "causal_performance_claim_eligible", "promotion_eligible",
                "promotion_passed"):
            self.assertIs(claim[name], False)
        favorable = audit.confidence_interval(
            [math.log(1.20)] * audit.ROUNDS)
        self.assertTrue(favorable[
            "promotion_lower_bound_at_least_1_05"])
        self.assertIs(fields["promotion_eligible"], False)
        self.assertIs(fields["promotion_passed"], False)

    def test_auditor_self_test_rejects_irrelevant_mode_option(self) -> None:
        completed = subprocess.run(
            [sys.executable, "-I", "-S", "-B", str(
                ROOT / "experiments/leopard2/main_compare" /
                "audit_v17_gfni_main_compare.py"),
             "--self-test", "--supervision", "absent"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            text=True,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("replay options", completed.stderr)

    def test_wrapper_separates_passive_from_active_contract(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        self.assertIn("--passive-shared-host", text)
        self.assertIn("--attempt N --attempt-budget 3", text)
        self.assertIn("passive-windowed-shared-host-observation/v1", text)
        self.assertIn("v18-gfni-main-passive-not-promoted-envelope/v1", text)
        self.assertIn("v18-gfni-main-failed-envelope/v1", text)
        self.assertIn("leopard2-passive-shared-host-policy/v2", text)
        self.assertIn("windowed_contamination", text)
        self.assertIn("passive_windowed_contamination_gate", text)
        self.assertIn(
            "windowed_benchmark_cpu_nonidle_excess_jiffies", text)
        self.assertIn(
            "out_of_window_reserved_sibling_nonidle_jiffies", text)
        self.assertIn("attempt_statement", text)
        self.assertIn("sealed_attempt_envelopes", text)
        self.assertIn('if [[ "$passive_mode" == true ||', text)
        self.assertIn("active_affinity_supervisor_executed:false", text)
        self.assertIn("--supervision windowed", text)
        self.assertIn("--supervision absent", text)
        self.assertIn("fresh active-v17 acquisition is disabled", text)
        self.assertIn("wrapper-process-and-owned-descendants-only", text)

    def test_wrapper_passive_contract_self_test_covers_float_timeout(
            self) -> None:
        wrapper_text = WRAPPER.read_text(encoding="utf-8")
        self.assertIn(
            '{"campaign":{"timeout_seconds":600.0}}', wrapper_text)
        self.assertIn(
            '"--timeout", $timeout_argument', wrapper_text)
        completed = subprocess.run(
            [str(WRAPPER), "--self-test-passive-contract"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(
            completed.stdout,
            "v18 passive wrapper contract self-test passed\n")


if __name__ == "__main__":
    unittest.main()
