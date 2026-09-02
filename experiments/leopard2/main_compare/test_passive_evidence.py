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


def v19_controller(raw: dict) -> dict:
    qualification = raw["pair_qualification"]
    selected = qualification["selected_pair"]
    benchmark_cpu = selected["benchmark_cpu"]
    reserved_sibling = selected["reserved_sibling"]
    launch = list(qualification["acquisition"]["allowed_cpu_set_at_launch"])
    verified_ns = qualification["bridge"]["bridge_finished_monotonic_ns"]
    return {
        "schema": census.CONTROLLER_SCHEMA_V19,
        "acquisition_generation": "passive-v3",
        "wrapper_pid": 100,
        "before_allowed_cpus": list(launch),
        "after_allowed_cpus": [
            cpu for cpu in launch
            if cpu not in (benchmark_cpu, reserved_sibling)
        ],
        "runner_launch_allowed_cpus": list(launch),
        "benchmark_cpu": benchmark_cpu,
        "reserved_sibling": reserved_sibling,
        "affinity_mutation_scope":
            "wrapper-process-and-owned-descendants-only",
        "active_affinity_supervisor_executed": False,
        "verified_acquisition_sha256": qualification["acquisition_sha256"],
        "verified_bridge_sha256": qualification["bridge_sha256"],
        "pair_verification_completed_monotonic_ns": verified_ns,
        "affinity_narrowing_started_monotonic_ns": verified_ns,
        "affinity_narrowing_finished_monotonic_ns": verified_ns,
    }


def v19_snapshot(raw: dict, phase: str) -> dict:
    qualification = raw["pair_qualification"]
    selected = qualification["selected_pair"]
    benchmark_cpu = selected["benchmark_cpu"]
    reserved_sibling = selected["reserved_sibling"]
    controller_value = v19_controller(raw)
    value = snapshot(phase)
    lease_uid = raw["isolation"]["pair_lease"]["payload"]["uid"]
    value["uid"] = lease_uid
    value["reserved_cpus"] = sorted((benchmark_cpu, reserved_sibling))
    value["collector"]["allowed_cpus"] = [
        controller_value["after_allowed_cpus"][0]]
    wrapper = value["same_uid_processes"]["entries"][0]
    wrapper["uid"] = lease_uid
    wrapper["cpus_allowed"] = list(controller_value["after_allowed_cpus"])
    wrapper["reserved_pair_intersection"] = []
    if phase == "pre":
        tail = qualification["bridge"]["campaign_presample_before"]
        started = tail["read_finished_monotonic_ns"] + 1
        finished = started + 1
        boundary = finished + 1
        first_started = qualification["first_window_handoff"] \
            ["first_window_before"]["read_started_monotonic_ns"]
        if boundary > first_started:
            raise RuntimeError("v19 fixture has no pre-census handoff gap")
        sources = {
            str(benchmark_cpu): tail["cpus"][str(benchmark_cpu)],
            str(reserved_sibling): tail["cpus"][str(reserved_sibling)],
        }
    else:
        outer = raw["isolation"]["after"]
        boundary = outer["monotonic_ns"] + 1
        started = boundary + 1
        finished = started + 1
        sources = {
            str(benchmark_cpu): outer["benchmark_cpu"],
            str(reserved_sibling): outer["reserved_sibling"],
        }
    value.update({
        "scan_started_monotonic_ns": started,
        "scan_finished_monotonic_ns": finished,
        "activity_boundary_monotonic_ns": boundary,
        "proc_stat": {
            key: census_cpu_from_raw(record)
            for key, record in sources.items()
        },
    })
    reseal_snapshot(value)
    return value


def rebuild_v19_isolation(raw: dict, *, before_ns: int | None = None) -> None:
    pair = raw["pair_qualification"]["selected_pair"]
    benchmark_cpu = pair["benchmark_cpu"]
    reserved_sibling = pair["reserved_sibling"]
    windows = [invocation["cpu_window"] for invocation in raw["invocations"]]
    current = raw["isolation"]
    raw["isolation"] = runner_fixtures.runner.isolation_record_v2(
        benchmark_cpu, reserved_sibling, current["pair_lease"],
        current["before"]["monotonic_ns"] if before_ns is None else before_ns,
        current["after"]["monotonic_ns"],
        current["before"]["benchmark_cpu"],
        current["after"]["benchmark_cpu"],
        current["before"]["reserved_sibling"],
        current["after"]["reserved_sibling"], windows)


def raw_evidence() -> dict:
    return {
        "schema": census.RAW_SCHEMA_V17,
        "supervision": None,
        "campaign": {
            "allowed_cpu_set_at_launch": [0, 1, 52, 116],
            "benchmark_cpu": 52,
            "reserved_sibling": 116,
        },
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


def relabel_windowed_raw_pair(
    raw: dict,
    cpu_pair: tuple[int, int],
) -> dict:
    """Relabel only pair-oriented v18 fields for internal-validator tests."""
    benchmark_cpu, reserved_sibling = cpu_pair
    raw = copy.deepcopy(raw)
    legacy_pair = set(audit.LEGACY_CPU_PAIR)
    raw["campaign"]["allowed_cpu_set_at_launch"] = sorted(
        (set(raw["campaign"]["allowed_cpu_set_at_launch"]) - legacy_pair) |
        set(cpu_pair))
    raw["campaign"].update(
        benchmark_cpu=benchmark_cpu,
        reserved_sibling=reserved_sibling,
    )
    reservation_payload = raw["reservation"]["payload"]
    reservation_payload.update(
        benchmark_cpu=benchmark_cpu,
        reserved_sibling=reserved_sibling,
    )
    raw["reservation"]["sha256"] = audit.sha256_bytes(
        audit.canonical_bytes(reservation_payload))

    isolation = raw["isolation"]
    isolation.update(
        benchmark_cpu=benchmark_cpu,
        reserved_sibling=reserved_sibling,
    )
    pair_lease = isolation["pair_lease"]
    pair_lease["payload"]["cpus"] = sorted(cpu_pair)
    pair_lease["path"] = (
        f"/run/user/{pair_lease['payload']['uid']}/leopard2-cpu-leases/"
        f"leopard2-cpu-pair-{pair_lease['payload']['uid']}-"
        f"{min(cpu_pair)}-{max(cpu_pair)}.lock"
    )
    pair_lease["sha256"] = audit.sha256_bytes(
        audit.canonical_bytes(pair_lease["payload"]))
    for endpoint in (isolation["before"], isolation["after"]):
        endpoint["benchmark_cpu"]["cpu"] = benchmark_cpu
        endpoint["reserved_sibling"]["cpu"] = reserved_sibling
    isolation["delta"]["benchmark_cpu"]["cpu"] = benchmark_cpu
    isolation["delta"]["reserved_sibling"]["cpu"] = reserved_sibling

    for invocation, retained in zip(
            raw["invocations"], isolation["invocation_windows"]):
        invocation["pinned_cpu"] = benchmark_cpu
        invocation["command"][2] = str(benchmark_cpu)
        for reservation_name in ("reservation_before", "reservation_after"):
            invocation[reservation_name] = copy.deepcopy(raw["reservation"])
        for window in (invocation["cpu_window"], retained):
            window.update(
                benchmark_cpu=benchmark_cpu,
                reserved_sibling=reserved_sibling,
            )
            for endpoint in (window["before"], window["after"]):
                endpoint["benchmark_cpu"]["cpu"] = benchmark_cpu
                endpoint["reserved_sibling"]["cpu"] = reserved_sibling
            window["delta"]["benchmark_cpu"]["cpu"] = benchmark_cpu
            window["delta"]["reserved_sibling"]["cpu"] = reserved_sibling
    return raw


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
    def test_legacy_census_policy_outputs_are_frozen(self) -> None:
        v17 = census.compare(
            snapshot("pre"), snapshot("post"), raw_evidence(), controller())
        self.assertEqual(
            census.digest_payload(v17),
            "d54e6f35278ae28e5f28d8a5597c422bf3a8f4cf106574c45e277eaa4d4220cb")

        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V18)
        v18 = census.compare(
            v18_snapshot(raw, "pre"),
            v18_snapshot(raw, "post", cpu52_extra=5, cpu116_extra=7),
            raw, v18_controller())
        self.assertEqual(
            census.digest_payload(v18),
            "f9a19b78a2d0fc005fcd60ddbb8b0be83aa989f607dff5dd71dcb27cfd759495")

    def test_v19_conditioned_policy_replays_dynamic_pair(self) -> None:
        with mock.patch.object(
                runner_fixtures.runner.os, "getuid", return_value=4242):
            raw = runner_fixtures.synthetic_raw(
                raw_schema=runner_fixtures.runner.RAW_SCHEMA_V19)
        self.assertEqual(v19_snapshot(raw, "pre")["uid"], 4242)
        policy = census.compare(
            v19_snapshot(raw, "pre"), v19_snapshot(raw, "post"), raw,
            v19_controller(raw),
            expected_v19_attempt=runner_fixtures.expected_v19_attempt_one())
        self.assertEqual(policy["schema"], census.POLICY_SCHEMA_V3)
        self.assertEqual(policy["acquisition_generation"], "passive-v3")
        self.assertEqual(policy["pair_qualification"]["selected_pair"], {
            "benchmark_cpu": 1, "reserved_sibling": 65})
        self.assertEqual(
            policy["pair_qualification"]["terminal"], "NOT_PROMOTED")
        self.assertTrue(policy["pair_qualification"]
                        ["fixed_point_replayed_independently"])
        self.assertTrue(policy["handoff_census"]["accepted"])
        self.assertFalse(policy["cpu_pair_exclusive"])
        self.assertFalse(policy["causal_performance_claim_eligible"])
        self.assertFalse(policy["promotion_eligible"])
        self.assertFalse(policy["promotion_passed"])
        self.assertEqual(policy["shared_host_claim_ceiling"], {
            "promotion_eligible": False,
            "host_exclusivity_proved": False,
            "whole_campaign_interval_observed": False,
            "causal_performance_claim_allowed": False,
        })

    def test_v19_compare_is_pure_and_performs_no_live_host_observation(
            self) -> None:
        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V19)
        with (
                mock.patch.object(
                    census, "bounded_read",
                    side_effect=AssertionError("unexpected host read")),
                mock.patch.object(
                    census.os, "sched_getaffinity",
                    side_effect=AssertionError("unexpected affinity read")),
                mock.patch.object(
                    census.time, "clock_gettime_ns",
                    side_effect=AssertionError("unexpected clock read"))):
            policy = census.compare(
                v19_snapshot(raw, "pre"), v19_snapshot(raw, "post"), raw,
                v19_controller(raw),
                expected_v19_attempt=
                    runner_fixtures.expected_v19_attempt_one())
        self.assertEqual(policy["schema"], census.POLICY_SCHEMA_V3)

    def test_v19_census_rejects_conditioning_and_join_splices(self) -> None:
        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V19)
        pre = v19_snapshot(raw, "pre")
        post = v19_snapshot(raw, "post")
        controller_value = v19_controller(raw)
        attempt = runner_fixtures.expected_v19_attempt_one()

        def reject(label: str, expected_error: str, mutation) -> None:
            candidate_pre = copy.deepcopy(pre)
            candidate_post = copy.deepcopy(post)
            candidate_raw = copy.deepcopy(raw)
            candidate_controller = copy.deepcopy(controller_value)
            candidate_attempt = copy.deepcopy(attempt)
            mutation(candidate_pre, candidate_post, candidate_raw,
                     candidate_controller, candidate_attempt)
            with self.subTest(splice=label):
                with self.assertRaises(census.CensusError) as caught:
                    census.compare(
                        candidate_pre, candidate_post, candidate_raw,
                        candidate_controller,
                        expected_v19_attempt=candidate_attempt)
                self.assertIn(expected_error, str(caught.exception))

        def snapshot_pair(pre_value, _post, _raw, _controller, _attempt):
            pre_value["reserved_cpus"] = [1, 66]
            reseal_snapshot(pre_value)

        def snapshot_time(pre_value, _post, raw_value, _controller, _attempt):
            first_ns = raw_value["pair_qualification"] \
                ["first_window_handoff"]["first_window_before"] \
                ["read_started_monotonic_ns"]
            pre_value["activity_boundary_monotonic_ns"] = first_ns + 1
            reseal_snapshot(pre_value)

        def snapshot_counter(
                pre_value, _post, _raw, _controller, _attempt):
            record = pre_value["proc_stat"]["1"]
            record["fields"]["idle"] -= 1
            record["total_jiffies"] -= 1
            reseal_snapshot(pre_value)

        def policy_splice(_pre, _post, raw_value, _controller, _attempt):
            raw_value["pair_qualification"]["policy"] \
                ["candidate_primary_cpus"].pop()

        def attempt_splice(_pre, _post, _raw, _controller, attempt_value):
            attempt_value.clear()
            attempt_value.update(census.v19_attempt_record(
                attempt=2,
                prior_attempt_failure_sha256="0" * 64,
                prior_attempt_acquisition_sha256="1" * 64,
                prior_attempt_selection_status=
                    "no-candidate-pair-qualified"))

        def pair_splice(_pre, _post, raw_value, _controller, _attempt):
            raw_value["pair_qualification"]["selected_pair"] \
                ["reserved_sibling"] = 66

        def bridge_splice(_pre, _post, raw_value, _controller, _attempt):
            raw_value["pair_qualification"]["bridge"] \
                ["bridge_tail_sha256"] = "0" * 64

        def handoff_splice(_pre, _post, raw_value, _controller, _attempt):
            raw_value["pair_qualification"]["first_window_handoff"] \
                ["handoff_elapsed_ns"] = 0

        def lineage_splice(_pre, _post, raw_value, _controller, _attempt):
            raw_value["pair_qualification"]["v18_failure_lineage"] \
                ["attempts"][0]["envelope_sha256sums_sha256"] = "0" * 64

        def host_splice(_pre, _post, raw_value, _controller, _attempt):
            for name in ("host_initial", "host_final"):
                raw_value[name]["system"]["release"] = "changed-fixture"

        def claim_splice(_pre, _post, raw_value, _controller, _attempt):
            raw_value["pair_qualification"]["shared_host_claim_ceiling"] \
                ["promotion_eligible"] = True

        def controller_splice(
                _pre, _post, _raw, controller_record, _attempt):
            controller_record["runner_launch_allowed_cpus"].pop()

        def controller_bridge_splice(
                _pre, _post, _raw, controller_record, _attempt):
            controller_record["verified_bridge_sha256"] = "0" * 64

        def controller_time_splice(
                _pre, _post, _raw, controller_record, _attempt):
            controller_record[
                "pair_verification_completed_monotonic_ns"] -= 1

        def controller_order_splice(
                _pre, _post, _raw, controller_record, _attempt):
            controller_record[
                "pair_verification_completed_monotonic_ns"] += 1

        def acquisition_campaign_splice(
                _pre, _post, raw_value, controller_record, _attempt):
            raw_value["campaign"]["allowed_cpu_set_at_launch"].append(128)
            controller_record["before_allowed_cpus"].append(128)
            controller_record["after_allowed_cpus"].append(128)
            controller_record["runner_launch_allowed_cpus"].append(128)

        def isolation_bridge_splice(
                _pre, _post, raw_value, _controller, _attempt):
            rebuild_v19_isolation(
                raw_value,
                before_ns=raw_value["isolation"]["before"]
                    ["monotonic_ns"] + 1)

        def first_window_splice(
                _pre, _post, raw_value, _controller, _attempt):
            invocation = raw_value["invocations"][0]
            old_window = invocation["cpu_window"]
            changed_before = copy.deepcopy(old_window["before"])
            for name in ("monotonic_ns", "read_started_monotonic_ns",
                         "read_finished_monotonic_ns"):
                changed_before[name] += 1
            invocation["cpu_window"] = runner_fixtures.runner.cpu_window_record(
                1, 65, invocation["cell_id"], invocation["round"],
                invocation["slot"], invocation["implementation"],
                changed_before, old_window["after"],
                old_window["child_started_monotonic_ns"],
                invocation["duration_ns"])
            rebuild_v19_isolation(raw_value)

        def lease_splice(_pre, _post, raw_value, _controller, _attempt):
            lease = raw_value["isolation"]["pair_lease"]
            lease["payload"]["cpus"] = [1, 66]
            uid = lease["payload"]["uid"]
            lease["path"] = (
                f"/run/user/{uid}/leopard2-cpu-leases/"
                f"leopard2-cpu-pair-{uid}-1-66.lock")
            lease["sha256"] = census.digest_payload(lease["payload"])

        def lease_owner_splice(
                _pre, _post, raw_value, _controller, _attempt):
            lease = raw_value["isolation"]["pair_lease"]
            lease["payload"]["uid"] += 1
            uid = lease["payload"]["uid"]
            lease["path"] = (
                f"/run/user/{uid}/leopard2-cpu-leases/"
                f"leopard2-cpu-pair-{uid}-1-65.lock")
            lease["sha256"] = census.digest_payload(lease["payload"])

        def status_splice(_pre, _post, raw_value, _controller, _attempt):
            raw_value["pair_qualification"]["record_status"] = "failed"

        mutations = (
            ("controller affinity",
             "passive controller affinity contract differs",
             controller_splice),
            ("controller bridge binding",
             "passive-v3 bridge verification and affinity narrowing",
             controller_bridge_splice),
            ("controller verification ordering",
             "passive-v3 bridge verification and affinity narrowing",
             controller_time_splice),
            ("controller verification before narrowing",
             "passive-v3 controller verification record differs",
             controller_order_splice),
            ("snapshot selected pair", "pre census host identity differs",
             snapshot_pair),
            ("snapshot handoff timestamp",
             "passive-v3 census is not nested", snapshot_time),
            ("snapshot handoff counter", "CPU 1 counters escape",
             snapshot_counter),
            ("qualification policy", "v19 qualification policy differs",
             policy_splice),
            ("external attempt authority", "v19 attempt differs",
             attempt_splice),
            ("selected pair",
             "v19 selected pair, acquisition hash, or status differs",
             pair_splice),
            ("qualification bridge", "v19 bridge differs",
             bridge_splice),
            ("first-window handoff", "v19 handoff differs",
             handoff_splice),
            ("v18 failure lineage", "v19 failure lineage differs",
             lineage_splice),
            ("host identity hash", "v19 qualification host identity differs",
             host_splice),
            ("claim ceiling", "v19 claim ceiling", claim_splice),
            ("qualification status",
             "requires complete NOT_PROMOTED v19 evidence", status_splice),
            ("acquisition/campaign affinity",
             "passive-v3 controller or campaign affinity differs",
             acquisition_campaign_splice),
            ("dynamic pair lease", "passive-v3 CPU pair lease identity differs",
             lease_splice),
            ("pair lease census owner",
             "passive-v3 pair lease owner differs from the census",
             lease_owner_splice),
            ("bridge/isolation join",
             "passive-v3 bridge tail or first-window handoff differs",
             isolation_bridge_splice),
            ("handoff/first-window join",
             "passive-v3 bridge tail or first-window handoff differs",
             first_window_splice),
        )
        for label, expected_error, mutation in mutations:
            reject(label, expected_error, mutation)
        if sys.flags.optimize == 0:
            completed = subprocess.run([
                sys.executable, "-O", "-I", "-S", "-B", str(Path(__file__)),
                "PassiveEvidenceTests."
                "test_v19_census_rejects_conditioning_and_join_splices",
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_v19_census_cli_accepts_external_attempt_and_dynamic_verify(
            self) -> None:
        raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V19)
        values = {
            "pre.json": v19_snapshot(raw, "pre"),
            "post.json": v19_snapshot(raw, "post"),
            "raw.json": raw,
            "controller.json": v19_controller(raw),
        }
        script = ROOT / "experiments/leopard2/main_compare" / \
            "passive_environment_census.py"
        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-census-cli-") as temporary:
            root = Path(temporary)
            for name, value in values.items():
                (root / name).write_bytes(census.canonical_bytes(value) + b"\n")
            for optimized in (False, True):
                output = root / f"policy-{'opt' if optimized else 'normal'}.json"
                command = [sys.executable]
                if optimized:
                    command.append("-O")
                command.extend((
                    "-I", "-S", "-B", str(script), "compare",
                    "--pre", str(root / "pre.json"),
                    "--post", str(root / "post.json"),
                    "--raw", str(root / "raw.json"),
                    "--controller", str(root / "controller.json"),
                    "--output", str(output), "--v19-attempt", "1"))
                completed = subprocess.run(
                    command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    check=False, text=True)
                with self.subTest(optimized=optimized):
                    self.assertEqual(
                        completed.returncode, 0, completed.stderr)
                    self.assertEqual(
                        census.load_json(output)["schema"],
                        census.POLICY_SCHEMA_V3)

            verified = subprocess.run([
                sys.executable, "-I", "-S", "-B", str(script), "verify",
                "--input", str(root / "pre.json"), "--phase", "pre",
                "--generation", "passive-v3", "--reserved-cpus", "1,65",
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                text=True)
            self.assertEqual(verified.returncode, 0, verified.stderr)

            missing_authority = subprocess.run([
                sys.executable, "-I", "-S", "-B", str(script), "compare",
                "--pre", str(root / "pre.json"),
                "--post", str(root / "post.json"),
                "--raw", str(root / "raw.json"),
                "--controller", str(root / "controller.json"),
                "--output", str(root / "missing-authority.json"),
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                text=True)
            self.assertNotEqual(missing_authority.returncode, 0)
            self.assertIn(
                "passive-v3 compare requires --v19-attempt",
                missing_authority.stderr)

            missing_pair = subprocess.run([
                sys.executable, "-I", "-S", "-B", str(script), "verify",
                "--input", str(root / "pre.json"), "--phase", "pre",
                "--generation", "passive-v3",
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                text=True)
            self.assertNotEqual(missing_pair.returncode, 0)
            self.assertIn(
                "passive-v3 verify requires --reserved-cpus",
                missing_pair.stderr)

            legacy_values = {
                "legacy-pre.json": snapshot("pre"),
                "legacy-post.json": snapshot("post"),
                "legacy-raw.json": raw_evidence(),
                "legacy-controller.json": controller(),
            }
            for name, value in legacy_values.items():
                (root / name).write_bytes(census.canonical_bytes(value) + b"\n")
            legacy_flag = subprocess.run([
                sys.executable, "-I", "-S", "-B", str(script), "compare",
                "--pre", str(root / "legacy-pre.json"),
                "--post", str(root / "legacy-post.json"),
                "--raw", str(root / "legacy-raw.json"),
                "--controller", str(root / "legacy-controller.json"),
                "--output", str(root / "legacy-policy.json"),
                "--v19-attempt", "1",
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                text=True)
            self.assertNotEqual(legacy_flag.returncode, 0)
            self.assertIn(
                "v19 attempt flags require a passive-v3 raw record",
                legacy_flag.stderr)

            legacy_pair = subprocess.run([
                sys.executable, "-I", "-S", "-B", str(script), "verify",
                "--input", str(root / "legacy-pre.json"), "--phase", "pre",
                "--generation", "passive-v1", "--reserved-cpus", "52,116",
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                text=True)
            self.assertNotEqual(legacy_pair.returncode, 0)
            self.assertIn(
                "--reserved-cpus is only valid for passive-v3 verify",
                legacy_pair.stderr)

    def test_census_self_test_stdout_is_frozen_normal_and_optimized(self) -> None:
        script = ROOT / "experiments/leopard2/main_compare" / \
            "passive_environment_census.py"
        for optimized in (False, True):
            command = [sys.executable]
            if optimized:
                command.append("-O")
            command.extend(("-I", "-S", "-B", str(script), "self-test"))
            completed = subprocess.run(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, text=True)
            with self.subTest(optimized=optimized):
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertEqual(
                    completed.stdout,
                    "passive environment census self-test passed\n")

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

    def test_pair_generic_internal_validators_preserve_oriented_roles(
            self) -> None:
        cpu_pair = (7, 3)
        raw = relabel_windowed_raw_pair(
            runner_fixtures.synthetic_raw(
                raw_schema=runner_fixtures.runner.RAW_SCHEMA_V18),
            cpu_pair,
        )
        windows = [invocation["cpu_window"]
                   for invocation in raw["invocations"]]

        self.assertEqual(
            audit.expected_command(
                "baseline", "/usr/bin/taskset", cpu_pair=cpu_pair)[2],
            "7")
        audit.validate_campaign(raw["campaign"], cpu_pair=cpu_pair)
        audit.validate_reservation(raw["reservation"], cpu_pair=cpu_pair)
        lease = audit.validate_pair_lease(
            raw["isolation"]["pair_lease"], cpu_pair=cpu_pair)
        self.assertEqual(lease["payload"]["cpus"], [3, 7])
        self.assertTrue(lease["path"].endswith("-3-7.lock"))
        audit.validate_cpu_window(
            windows[0], (0, 0, "baseline"),
            raw["invocations"][0]["duration_ns"], cpu_pair=cpu_pair)
        audit.validate_isolation_v2(
            raw["isolation"], windows, cpu_pair=cpu_pair)
        observation = audit.windowed_observation(
            raw["isolation"], cpu_pair=cpu_pair)
        self.assertEqual(
            (observation["benchmark_cpu"], observation["reserved_sibling"]),
            cpu_pair)

        census.validate_v18_isolation(
            raw["isolation"], raw, cpu_pair=cpu_pair)
        controller_value = v18_controller()
        controller_value.update(
            before_allowed_cpus=[0, 3, 7],
            after_allowed_cpus=[0],
            runner_launch_allowed_cpus=[0, 3, 7],
            benchmark_cpu=cpu_pair[0],
            reserved_sibling=cpu_pair[1],
        )
        census.validate_controller(
            controller_value, census.RAW_SCHEMA_V18, cpu_pair=cpu_pair)
        generic_snapshot = snapshot("pre")
        generic_snapshot["reserved_cpus"] = [3, 7]
        generic_snapshot["proc_stat"] = {
            "3": cpu_record(3), "7": cpu_record(7)}
        reseal_snapshot(generic_snapshot)
        census.validate_snapshot(
            generic_snapshot, "pre", cpu_pair=cpu_pair)

        isolation = raw["isolation"]
        supervision = {
            "schema": audit.SUPERVISION_SCHEMA,
            "execution_nonce": "a" * 64,
            "runner_pid": 123,
            "runner_started_monotonic_ns": 0,
            "runner_finished_monotonic_ns":
                isolation["after"]["monotonic_ns"] + 1,
            "launch_cpus": raw["campaign"]["allowed_cpu_set_at_launch"],
            "reserved_cpus": [3, 7],
            "campaign_sha256": audit.sha256_bytes(
                audit.canonical_bytes(raw["campaign"])),
            "reservation_sha256": raw["reservation"]["sha256"],
            "reservation_nonce": raw["reservation"]["payload"]["nonce"],
            "isolation_before_monotonic_ns":
                isolation["before"]["monotonic_ns"],
            "isolation_after_monotonic_ns":
                isolation["after"]["monotonic_ns"],
        }
        audit.validate_supervision(
            supervision, raw["campaign"], raw["reservation"], isolation,
            cpu_pair=cpu_pair)
        unsorted_supervision = copy.deepcopy(supervision)
        unsorted_supervision["reserved_cpus"] = [7, 3]
        with self.assertRaises(audit.AuditError):
            audit.validate_supervision(
                unsorted_supervision, raw["campaign"], raw["reservation"],
                isolation, cpu_pair=cpu_pair)

        legacy_raw = runner_fixtures.synthetic_raw(
            raw_schema=runner_fixtures.runner.RAW_SCHEMA_V18)
        reversed_host = copy.deepcopy(legacy_raw["host_initial"])
        reversed_host["benchmark_cpu"], reversed_host["reserved_sibling"] = (
            reversed_host["reserved_sibling"],
            reversed_host["benchmark_cpu"],
        )
        for record in (
                reversed_host["benchmark_cpu"],
                reversed_host["reserved_sibling"]):
            record["cpuinfo"]["flags"] = "avx2 gfni"
        audit.validate_host(
            reversed_host, "reversed", [0, 52, 116],
            cpu_pair=(116, 52))

        self.assertEqual(
            audit.cpu_pair_for_raw_schema(audit.RAW_SCHEMA_V17),
            audit.LEGACY_CPU_PAIR)
        self.assertEqual(
            census.cpu_pair_for_generation("passive-v2"),
            census.LEGACY_CPU_PAIR)
        with self.assertRaises(audit.AuditError):
            audit.cpu_pair_for_raw_schema("leopard2-main-compare-raw/v19")
        with self.assertRaises(census.CensusError):
            census.cpu_pair_for_raw_schema(
                "leopard2-main-compare-raw/v19")
        with self.assertRaises(census.CensusError):
            census.cpu_pair_for_generation("passive-v3")
        with self.assertRaises(census.CensusError):
            census.validate_controller(
                controller_value, "leopard2-main-compare-raw/v19",
                cpu_pair=cpu_pair)

    def test_census_uses_oriented_pair_not_sorted_pair_for_v17_gate(
            self) -> None:
        cpu_pair = (116, 52)
        raw = raw_evidence()
        raw["campaign"].update(
            benchmark_cpu=cpu_pair[0], reserved_sibling=cpu_pair[1])
        isolation = raw["isolation"]
        isolation.update(
            benchmark_cpu=cpu_pair[0], reserved_sibling=cpu_pair[1])
        for endpoint in (isolation["before"], isolation["after"]):
            endpoint["benchmark_cpu"]["cpu"] = cpu_pair[0]
            endpoint["reserved_sibling"]["cpu"] = cpu_pair[1]

        pre, post = snapshot("pre"), snapshot("post")
        pre["proc_stat"] = {
            "52": cpu_record(52), "116": cpu_record(116)}
        post["proc_stat"] = {
            "52": cpu_record(52, idle=1),
            "116": cpu_record(116, user=17),
        }
        reseal_snapshot(pre)
        reseal_snapshot(post)
        controller_value = controller()
        controller_value.update(
            benchmark_cpu=cpu_pair[0], reserved_sibling=cpu_pair[1])

        with mock.patch.object(
                census, "cpu_pair_for_raw_schema", return_value=cpu_pair):
            policy = census.compare(pre, post, raw, controller_value)
        self.assertEqual(
            policy["outer_contamination"]
            ["benchmark_cpu_nonidle_jiffies"],
            17)
        self.assertEqual(
            policy["outer_contamination"]
            ["reserved_sibling_nonidle_jiffies"],
            0)

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
