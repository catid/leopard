#!/usr/bin/python3

"""Capture and verify read-only endpoint censuses for passive v17-v19 evidence.

This controller never changes affinity and never signals another process.  Its
snapshots disclose affinity eligibility at two non-atomic endpoints; they do
not claim interval completeness or CPU exclusivity.
"""

from __future__ import annotations

import argparse
import errno
import hashlib
import json
import os
import re
import stat
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-passive-same-uid-endpoint-census/v1"
POLICY_SCHEMA_V1 = "leopard2-passive-shared-host-policy/v1"
POLICY_SCHEMA_V2 = "leopard2-passive-shared-host-policy/v2"
POLICY_SCHEMA_V3 = "leopard2-passive-shared-host-policy/v3"
POLICY_SCHEMA = POLICY_SCHEMA_V1
RAW_SCHEMA_V17 = "leopard2-main-compare-raw/v17"
RAW_SCHEMA_V18 = "leopard2-main-compare-raw/v18"
RAW_SCHEMA_V19 = "leopard2-main-compare-raw/v19"
ISOLATION_SCHEMA_V1 = "leopard2-main-compare-isolation/v1"
ISOLATION_SCHEMA_V2 = "leopard2-main-compare-isolation/v2"
CPU_WINDOW_SCHEMA = "leopard2-main-compare-invocation-cpu-window/v1"
PAIR_LEASE_SCHEMA = "leopard2-cpu-pair-lease/v1"
CONTROLLER_SCHEMA_V17 = "leopard2-v17-passive-controller-affinity/v1"
CONTROLLER_SCHEMA_V18 = "leopard2-v18-passive-controller-affinity/v1"
CONTROLLER_SCHEMA_V19 = "leopard2-v19-passive-controller-affinity/v1"
SEMANTICS = "non-atomic-endpoint-affinity-eligibility"
MUTATION_POLICY = "read-only; no affinity mutation, signal, or procfs write"
MAX_CPU_ID = 1_048_575
MAX_CPU_LIST_ENTRIES = 4096
MAX_STATUS_BYTES = 1 << 20
MAX_STAT_BYTES = 1 << 20
MAX_TASKS = 1_000_000
CLOCK_TICKS_PER_SECOND = 100
LEGACY_V17_CPU52_AGGREGATE_EXCESS_LIMIT_JIFFIES = 16
EXPECTED_INVOCATION_COUNT = 12
LEGACY_BENCHMARK_CPU = 52
LEGACY_RESERVED_SIBLING = 116
LEGACY_CPU_PAIR = (LEGACY_BENCHMARK_CPU, LEGACY_RESERVED_SIBLING)
HEX256 = re.compile(r"^[0-9a-f]{64}$")
BOOT_ID = re.compile(
    r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$")
CPU_FIELDS = (
    "user", "nice", "system", "idle", "iowait", "irq", "softirq",
    "steal", "guest", "guest_nice",
)
NONIDLE_FIELDS = ("user", "nice", "system", "irq", "softirq", "steal")
V19_COUNTER_FIELDS = CPU_FIELDS[:8]
V19_PAIR_KEYS = {"benchmark_cpu", "reserved_sibling"}
V19_POLICY_KEYS = {
    "schema", "clock_ticks_per_second", "candidate_primary_cpus",
    "excluded_pairs", "domain_mode", "candidate_rule", "selection_rule",
    "frozen_pair_rule", "observation",
}
V19_OBSERVATION_KEYS = {
    "counter_source", "counter_fields", "idle_fields", "nonidle_fields",
    "tick_cap_formula", "liveness_lower_slack_jiffies",
    "nonidle_limit_jiffies", "liveness_rule",
}
V19_TOPOLOGY_KEYS = {"schema", "cpus"}
V19_TOPOLOGY_CPU_KEYS = {
    "cpu", "online", "physical_package_id", "die_id", "core_id",
    "thread_siblings", "domain_cpus",
}
V19_SNAPSHOT_KEYS = {
    "cpus", "read_started_monotonic_ns", "read_finished_monotonic_ns",
}
V19_COUNTER_KEYS = {"cpu", "fields", "total_jiffies"}
V19_DELTA_KEYS = {
    "cpu", "fields", "idle_jiffies", "nonidle_jiffies", "total_jiffies",
}
V19_WINDOW_KEYS = {
    "schema", "phase", "index", "before", "after", "window_ns",
    "wide_window_ns", "elapsed_tick_cap", "wide_elapsed_tick_cap",
    "liveness_lower_bound_jiffies", "delta", "nonidle_cpus",
    "not_live_cpus", "clock_ticks_per_second", "policy",
}
V19_AGGREGATE_KEYS = {
    *V19_DELTA_KEYS, "nonidle_window_count", "not_live_window_count",
    "overall_live",
}
V19_CANDIDATE_KEYS = {
    "benchmark_cpu", "reserved_sibling", "domain_cpus", "pair_idle",
    "pair_live", "domain_idle", "domain_live", "qualified",
}
V19_SCAN_KEYS = {
    "schema", "policy", "allowed_cpu_set_at_launch", "topology_before",
    "topology_after", "topology_reread_identical", "windows",
    "cpu_aggregates", "candidate_pairs", "eligible_pairs", "selected",
    "frozen_pair_from_prior_attempt", "selection_status",
    "scan_started_monotonic_ns", "scan_finished_monotonic_ns",
    "realized_scan_total_ns", "wide_scan_total_ns",
    "overall_elapsed_tick_cap", "overall_wide_elapsed_tick_cap",
    "overall_liveness_lower_bound_jiffies", "scan_window_count",
    "max_proc_stat_read_ns", "excluded_pair_count", "candidate_pair_count",
    "eligible_pair_count", "candidate_timing_performed",
}
V19_ACQUISITION_KEYS = {
    "schema", "acquisition_method", "sources", "policy", "policy_sha256",
    "requested_window_count", "nominal_window_ns",
    "frozen_pair_from_prior_attempt", "allowed_cpu_set_at_launch",
    "allowed_cpu_set_after_scan", "clock_ticks_per_second_at_launch",
    "clock_ticks_per_second_after_scan", "topology_before_sha256",
    "topology_after_sha256", "scan", "scan_sha256",
    "host_mutation_performed", "candidate_timing_performed",
}
V19_BRIDGE_GEOMETRY_KEYS = {
    "minimum_window_count", "maximum_window_count", "nominal_window_ns",
    "maximum_handoff_elapsed_ns",
}
V19_QUALIFICATION_GEOMETRY_KEYS = {"window_count", "nominal_window_ns"}
V19_BRIDGE_KEYS = {
    "schema", "handoff_rule", "acceptance_rule", "acquisition_sha256",
    "scan_sha256", "policy_sha256", "selected_pair",
    "frozen_pair_from_prior_attempt", "acquisition_geometry",
    "acquisition_geometry_sha256", "bridge_geometry",
    "bridge_geometry_sha256", "observed_cpus", "guarded_cpus",
    "scan_tail_sha256", "windows", "cpu_aggregates",
    "bridge_head_sha256", "bridge_tail_sha256",
    "campaign_presample_before", "campaign_presample_before_sha256",
    "bridge_started_monotonic_ns", "bridge_finished_monotonic_ns",
    "bridge_deadline_monotonic_ns", "realized_bridge_total_ns",
    "wide_bridge_total_ns", "overall_elapsed_tick_cap",
    "overall_wide_elapsed_tick_cap", "overall_liveness_lower_bound_jiffies",
    "bridge_window_count", "nonidle_guarded_cpus", "not_live_guarded_cpus",
    "bridge_accepted", "host_mutation_performed",
    "candidate_timing_performed", "shared_host_claim_ceiling",
}
V19_ATTEMPT_KEYS = {
    "schema", "attempt", "attempt_budget", "frozen_pair_from_prior_attempt",
    "prior_attempt_failure_sha256", "prior_attempt_acquisition_sha256",
    "prior_attempt_selection_status", "prior_attempt_selected_pair",
    "fresh_selection_permitted", "selection_rule", "frozen_pair_rule",
}
V19_HANDOFF_KEYS = {
    "schema", "selected_pair", "bridge_tail_sha256",
    "bridge_tail_read_finished_monotonic_ns", "first_window_before",
    "first_window_before_sha256", "first_window_before_read_started_monotonic_ns",
    "handoff_elapsed_ns", "maximum_handoff_elapsed_ns",
    "selected_pair_nonidle_delta", "accepted", "failure_terminal",
}
V19_FIRST_WINDOW_KEYS = {
    "benchmark_cpu", "monotonic_ns", "read_finished_monotonic_ns",
    "read_started_monotonic_ns", "reserved_sibling",
}
V19_LINEAGE_KEYS = {"schema", "source_commit", "source_tree", "attempts"}
V19_LINEAGE_ENTRY_KEYS = {
    "attempt", "envelope", "envelope_sha256sums_sha256",
}
V19_CLAIM_KEYS = {
    "promotion_eligible", "host_exclusivity_proved",
    "whole_campaign_interval_observed", "causal_performance_claim_allowed",
}
V19_RECORD_KEYS = {
    "schema", "stage", "record_status", "attempt", "policy",
    "policy_sha256", "host_identity_sha256", "qualification_geometry",
    "qualification_geometry_sha256", "acquisition", "acquisition_sha256",
    "selected_pair", "selection_status", "bridge", "bridge_sha256",
    "bridge_geometry_sha256", "first_window_handoff",
    "v18_failure_lineage", "v18_failure_lineage_sha256",
    "shared_host_claim_ceiling", "terminal", "host_mutation_performed",
    "candidate_timing_performed",
}


V19_MAX_COUNTER = (1 << 64) - 1
V19_MAX_MONOTONIC_NS = (1 << 63) - 1
V19_MAX_WINDOWS = 4096
V19_MAX_JSON_BYTES = 64 * 1024 * 1024
V19_POLICY_SCHEMA = "leopard2-pair-qualification-policy/v1"
V19_TOPOLOGY_SCHEMA = "leopard2-pair-qualification-topology/v1"
V19_WINDOW_SCHEMA = "leopard2-pair-qualification-window/v1"
V19_SCAN_SCHEMA = "leopard2-pair-qualification-scan/v1"
V19_ACQUISITION_SCHEMA = "leopard2-pair-qualification-acquisition/v1"
V19_ACQUISITION_METHOD = \
    "read-only-linux-sysfs-proc-stat-shared-snapshot-chain/v1"
V19_BRIDGE_SCHEMA = "leopard2-pair-qualification-bridge/v1"
V19_HANDOFF_RULE = \
    "exact-scan-tail-shared-n-plus-one-presample-before/v1"
V19_BRIDGE_ACCEPTANCE_RULE = \
    "selected-pair-or-domain-zero-nonidle-live-fixed-point/v1"
V19_RECORD_SCHEMA = "leopard2-main-compare-pair-qualification/v1"
V19_ATTEMPT_SCHEMA = "leopard2-main-compare-pair-qualified-attempt/v1"
V19_HANDOFF_SCHEMA = "leopard2-main-compare-first-window-handoff/v1"
V19_LINEAGE_SCHEMA = "leopard2-main-compare-v18-failure-lineage/v1"
V19_CANDIDATE_RULE = "explicit-primary-symmetric-smt2-online-allowed/v1"
V19_SELECTION_RULE = "lowest-primary-unless-frozen/v1"
V19_FROZEN_PAIR_RULE = "require-frozen-when-present-no-fallback/v1"
V19_TICK_CAP_FORMULA = \
    "(elapsed_ns * clock_ticks_per_second) // 1000000000 + 1"
V19_LIVENESS_LOWER_SLACK_JIFFIES = 2
V19_NONIDLE_LIMIT_JIFFIES = 0
V19_LIVENESS_RULE = (
    "narrow_tick_cap - 2 >= 1 and narrow_tick_cap - 2 <= "
    "total_jiffies_delta <= wide_tick_cap"
)
V19_COUNTER_SOURCE = "/proc/stat"
V19_IDLE_FIELDS = ("idle", "iowait")
V19_NONIDLE_FIELDS = (
    "user", "nice", "system", "irq", "softirq", "steal",
)
V19_CANDIDATE_PRIMARY_CPUS = tuple(range(1, 64))
V19_ATTEMPT_BUDGET = 2
V19_QUALIFICATION_WINDOW_COUNT = 12
V19_QUALIFICATION_NOMINAL_WINDOW_NS = 7_000_000_000
V19_BRIDGE_WINDOW_COUNT = 2
V19_BRIDGE_NOMINAL_WINDOW_NS = 1_000_000_000
V19_MAXIMUM_HANDOFF_ELAPSED_NS = 5_000_000_000
V19_SUCCESS_TERMINAL = "NOT_PROMOTED"
V19_SOURCE_ITEMS = (
    ("allowed_cpu_set_at_launch", "sched_getaffinity(0)"),
    ("clock_ticks_per_second", "sysconf(SC_CLK_TCK)"),
    ("counter_source", "/proc/stat"),
    ("monotonic_clock", "time.monotonic_ns"),
    ("online_cpus", "/sys/devices/system/cpu/online"),
    ("physical_package_id",
     "/sys/devices/system/cpu/cpu{cpu}/topology/physical_package_id"),
    ("die_id", "/sys/devices/system/cpu/cpu{cpu}/topology/die_id"),
    ("core_id", "/sys/devices/system/cpu/cpu{cpu}/topology/core_id"),
    ("thread_siblings",
     "/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"),
    ("domain_cpus",
     "/sys/devices/system/cpu/cpu{cpu}/cache/index3/shared_cpu_list"),
)
V19_V18_SOURCE_COMMIT = "c8f825d0a033d31d220b0ebce9cc8871e8c2fc6d"
V19_V18_SOURCE_TREE = "2c17a0a7bcea20274d2593cb204442c4c817e464"
V19_V18_FAILURES = (
    (1, ".research/leopard-79h/c8f825d-v18-passive-main-a1",
     "ce65c3a49ef1c1d89ba51ea03d0af4742d6790e6f2ea2662917d9ef9a9d945d7"),
    (2, ".research/leopard-79h/c8f825d-v18-passive-main-a2",
     "a1bf0eda157c251f33f7260ebd76931d88054d460bd07a97bcba2811384b2c10"),
    (3, ".research/leopard-79h/c8f825d-v18-passive-main-a3",
     "fe5b40cc98753cbd794ee019cb0e2643d0ccee0aca4c5fd7b2e0b27df8a86139"),
)


class CensusError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CensusError(message)


def cpu_pair_for_raw_schema(raw_schema: str) -> tuple[int, int]:
    require(raw_schema in (RAW_SCHEMA_V17, RAW_SCHEMA_V18),
            "passive campaign raw schema has no frozen CPU-pair contract")
    return LEGACY_CPU_PAIR


def cpu_pair_for_generation(generation: str) -> tuple[int, int]:
    require(generation in ("passive-v1", "passive-v2"),
            "passive census generation has no frozen CPU-pair contract")
    return LEGACY_CPU_PAIR


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def exact_json_equal(left: Any, right: Any) -> bool:
    """Compare JSON values without Python's bool/int coercion."""
    return canonical_bytes(left) == canonical_bytes(right)


def digest_payload(value: Mapping[str, Any]) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def v19_exact_json_equal(left: Any, right: Any) -> bool:
    """Compare v19 JSON without Python bool/int or int/float coercion."""
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return set(left) == set(right) and all(
            v19_exact_json_equal(left[key], right[key]) for key in left)
    if type(left) is list:
        return len(left) == len(right) and all(
            v19_exact_json_equal(a, b) for a, b in zip(left, right))
    return left == right


def v19_canonical_bytes(value: Any) -> bytes:
    """Mirror the v19 contract spelling, including Unicode and final LF."""
    try:
        return (json.dumps(
            value, allow_nan=False, ensure_ascii=False,
            separators=(",", ":"), sort_keys=True) + "\n").encode(
                "utf-8", errors="strict")
    except (TypeError, ValueError, UnicodeError, RecursionError) as error:
        raise CensusError(f"v19 value is not canonical finite JSON: {error}") \
            from error


def v19_sha256(value: Any) -> str:
    return hashlib.sha256(v19_canonical_bytes(value)).hexdigest()


def v19_bounded_int(
    value: Any, minimum: int, maximum: int, label: str,
) -> int:
    require(type(value) is int and minimum <= value <= maximum,
            f"{label} is not an exact bounded integer")
    return value


def v19_hex64(value: Any, label: str) -> str:
    require(type(value) is str and HEX256.fullmatch(value) is not None,
            f"{label} is not exact lowercase SHA-256")
    return value


def v19_cpu_list(
    value: Any, label: str, *, allow_empty: bool = False,
) -> list[int]:
    require(type(value) is list and (allow_empty or value) and
            len(value) <= MAX_CPU_LIST_ENTRIES,
            f"{label} is not a bounded CPU list")
    result = [
        v19_bounded_int(cpu, 0, MAX_CPU_ID, f"{label} entry")
        for cpu in value
    ]
    require(result == sorted(set(result)), f"{label} is not canonical")
    return result


def v19_pair(value: Any, label: str, *, allow_none: bool = False
             ) -> dict[str, int] | None:
    if value is None:
        require(allow_none, f"{label} is absent")
        return None
    pair = exact_keys(value, V19_PAIR_KEYS, label)
    benchmark = v19_bounded_int(
        pair["benchmark_cpu"], 0, MAX_CPU_ID, f"{label} benchmark CPU")
    sibling = v19_bounded_int(
        pair["reserved_sibling"], 0, MAX_CPU_ID, f"{label} reserved sibling")
    require(benchmark != sibling, f"{label} repeats one logical CPU")
    return {"benchmark_cpu": benchmark, "reserved_sibling": sibling}


def v19_claim_ceiling() -> dict[str, bool]:
    return {
        "promotion_eligible": False,
        "host_exclusivity_proved": False,
        "whole_campaign_interval_observed": False,
        "causal_performance_claim_allowed": False,
    }


def v19_observation_policy() -> dict[str, Any]:
    return {
        "counter_source": V19_COUNTER_SOURCE,
        "counter_fields": list(V19_COUNTER_FIELDS),
        "idle_fields": list(V19_IDLE_FIELDS),
        "nonidle_fields": list(V19_NONIDLE_FIELDS),
        "tick_cap_formula": V19_TICK_CAP_FORMULA,
        "liveness_lower_slack_jiffies":
            V19_LIVENESS_LOWER_SLACK_JIFFIES,
        "nonidle_limit_jiffies": V19_NONIDLE_LIMIT_JIFFIES,
        "liveness_rule": V19_LIVENESS_RULE,
    }


def v19_policy_record() -> dict[str, Any]:
    return {
        "schema": V19_POLICY_SCHEMA,
        "clock_ticks_per_second": CLOCK_TICKS_PER_SECOND,
        "candidate_primary_cpus": list(V19_CANDIDATE_PRIMARY_CPUS),
        "excluded_pairs": [],
        "domain_mode": "pair-only",
        "candidate_rule": V19_CANDIDATE_RULE,
        "selection_rule": V19_SELECTION_RULE,
        "frozen_pair_rule": V19_FROZEN_PAIR_RULE,
        "observation": v19_observation_policy(),
    }


def validate_v19_policy(value: Any) -> dict[str, Any]:
    policy = exact_keys(value, V19_POLICY_KEYS, "v19 qualification policy")
    exact_keys(
        policy["observation"], V19_OBSERVATION_KEYS,
        "v19 qualification observation policy")
    expected = v19_policy_record()
    require(v19_exact_json_equal(policy, expected),
            "v19 qualification policy differs from preregistration")
    return expected


def v19_attempt_record(
    *, attempt: int, prior_attempt_failure_sha256: Any = None,
    prior_attempt_acquisition_sha256: Any = None,
    prior_attempt_selection_status: Any = None,
    prior_attempt_selected_pair: Any = None,
) -> dict[str, Any]:
    number = v19_bounded_int(
        attempt, 1, V19_ATTEMPT_BUDGET, "v19 attempt")
    selected = v19_pair(
        prior_attempt_selected_pair, "v19 prior selected pair",
        allow_none=True)
    if number == 1:
        require(prior_attempt_failure_sha256 is None and
                prior_attempt_acquisition_sha256 is None and
                prior_attempt_selection_status is None and selected is None,
                "v19 attempt 1 carries prior-attempt state")
        failure_hash = acquisition_hash = selection_status = None
        frozen = None
        fresh = True
    else:
        failure_hash = v19_hex64(
            prior_attempt_failure_sha256, "v19 prior attempt failure hash")
        require(prior_attempt_selection_status in (
                    "not-acquired", "no-candidate-pair-qualified",
                    "selected-lowest-primary"),
                "v19 prior selection status differs")
        selection_status = prior_attempt_selection_status
        if selection_status == "not-acquired":
            require(prior_attempt_acquisition_sha256 is None and
                    selected is None,
                    "v19 not-acquired prior attempt carries acquisition")
            acquisition_hash = None
        else:
            acquisition_hash = v19_hex64(
                prior_attempt_acquisition_sha256,
                "v19 prior attempt acquisition hash")
        if selection_status == "selected-lowest-primary":
            require(selected is not None,
                    "v19 selected prior attempt lacks its pair")
            frozen = dict(selected)
            fresh = False
        else:
            require(selected is None,
                    "v19 unselected prior attempt carries a pair")
            frozen = None
            fresh = True
    return {
        "schema": V19_ATTEMPT_SCHEMA,
        "attempt": number,
        "attempt_budget": V19_ATTEMPT_BUDGET,
        "frozen_pair_from_prior_attempt": frozen,
        "prior_attempt_failure_sha256": failure_hash,
        "prior_attempt_acquisition_sha256": acquisition_hash,
        "prior_attempt_selection_status": selection_status,
        "prior_attempt_selected_pair": selected,
        "fresh_selection_permitted": fresh,
        "selection_rule": V19_SELECTION_RULE,
        "frozen_pair_rule": V19_FROZEN_PAIR_RULE,
    }


def validate_v19_attempt(value: Any, expected_value: Any) -> dict[str, Any]:
    attempt = exact_keys(value, V19_ATTEMPT_KEYS, "v19 attempt record")
    expected = exact_keys(
        expected_value, V19_ATTEMPT_KEYS, "expected v19 attempt record")
    rebuilt = v19_attempt_record(
        attempt=expected["attempt"],
        prior_attempt_failure_sha256=
            expected["prior_attempt_failure_sha256"],
        prior_attempt_acquisition_sha256=
            expected["prior_attempt_acquisition_sha256"],
        prior_attempt_selection_status=
            expected["prior_attempt_selection_status"],
        prior_attempt_selected_pair=expected["prior_attempt_selected_pair"])
    require(v19_exact_json_equal(expected, rebuilt) and
            v19_exact_json_equal(attempt, rebuilt),
            "v19 attempt differs from external authority")
    return rebuilt


def v19_lineage_record() -> dict[str, Any]:
    return {
        "schema": V19_LINEAGE_SCHEMA,
        "source_commit": V19_V18_SOURCE_COMMIT,
        "source_tree": V19_V18_SOURCE_TREE,
        "attempts": [{
            "attempt": attempt,
            "envelope": envelope,
            "envelope_sha256sums_sha256": digest,
        } for attempt, envelope, digest in V19_V18_FAILURES],
    }


def validate_v19_lineage(value: Any) -> dict[str, Any]:
    lineage = exact_keys(value, V19_LINEAGE_KEYS, "v19 failure lineage")
    require(type(lineage["attempts"]) is list and
            len(lineage["attempts"]) == len(V19_V18_FAILURES),
            "v19 failure lineage attempt count differs")
    for entry in lineage["attempts"]:
        exact_keys(entry, V19_LINEAGE_ENTRY_KEYS,
                   "v19 failure lineage entry")
    expected = v19_lineage_record()
    require(v19_exact_json_equal(lineage, expected),
            "v19 failure lineage differs from frozen v18 failures")
    return expected


def v19_counter_fields(value: Any, label: str) -> dict[str, int]:
    require(type(value) is dict and set(value) == set(V19_COUNTER_FIELDS),
            f"{label} counter fields differ")
    return {
        field: v19_bounded_int(
            value[field], 0, V19_MAX_COUNTER, f"{label} {field}")
        for field in V19_COUNTER_FIELDS
    }


def v19_counter(value: Any, expected_cpu: int, label: str
                ) -> dict[str, Any]:
    record = exact_keys(value, V19_COUNTER_KEYS, label)
    require(type(record["cpu"]) is int and record["cpu"] == expected_cpu,
            f"{label} CPU differs")
    fields = v19_counter_fields(record["fields"], label)
    total = v19_bounded_int(
        record["total_jiffies"], 0, V19_MAX_COUNTER, f"{label} total")
    expected = {"cpu": expected_cpu, "fields": fields,
                "total_jiffies": total}
    require(total == sum(fields.values()) and
            v19_exact_json_equal(record, expected),
            f"{label} total or canonical shape differs")
    return expected


def v19_snapshot(
    value: Any, expected_cpus: Sequence[int], label: str,
) -> dict[str, Any]:
    snapshot = exact_keys(value, V19_SNAPSHOT_KEYS, label)
    cpus = list(expected_cpus)
    require(cpus == sorted(set(cpus)) and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in cpus),
            f"{label} expected CPU universe differs")
    started = v19_bounded_int(
        snapshot["read_started_monotonic_ns"], 0, V19_MAX_MONOTONIC_NS,
        f"{label} read start")
    finished = v19_bounded_int(
        snapshot["read_finished_monotonic_ns"], started,
        V19_MAX_MONOTONIC_NS, f"{label} read finish")
    table = snapshot["cpus"]
    require(type(table) is dict and set(table) == {str(cpu) for cpu in cpus},
            f"{label} CPU table differs")
    expected = {
        "cpus": {
            str(cpu): v19_counter(
                table[str(cpu)], cpu, f"{label} CPU {cpu}")
            for cpu in cpus
        },
        "read_started_monotonic_ns": started,
        "read_finished_monotonic_ns": finished,
    }
    require(v19_exact_json_equal(snapshot, expected),
            f"{label} is not canonical")
    return expected


def v19_counter_delta(
    before: Mapping[str, Any], after: Mapping[str, Any], label: str,
) -> dict[str, Any]:
    cpu = before["cpu"]
    require(type(cpu) is int and after.get("cpu") == cpu and
            type(after.get("cpu")) is int,
            f"{label} CPU changed")
    fields: dict[str, int] = {}
    for field in V19_COUNTER_FIELDS:
        first = before["fields"][field]
        last = after["fields"][field]
        require(last >= first, f"{label} {field} counter regressed")
        fields[field] = last - first
    idle = sum(fields[field] for field in V19_IDLE_FIELDS)
    nonidle = sum(fields[field] for field in V19_NONIDLE_FIELDS)
    total = sum(fields.values())
    v19_bounded_int(total, 0, V19_MAX_COUNTER, f"{label} total")
    return {
        "cpu": cpu, "fields": fields, "idle_jiffies": idle,
        "nonidle_jiffies": nonidle, "total_jiffies": total,
    }


def v19_topology(value: Any) -> dict[str, Any]:
    topology = exact_keys(value, V19_TOPOLOGY_KEYS, "v19 topology")
    raw_entries = topology["cpus"]
    require(topology["schema"] == V19_TOPOLOGY_SCHEMA and
            type(raw_entries) is list and
            0 < len(raw_entries) <= MAX_CPU_LIST_ENTRIES,
            "v19 topology shape differs")
    entries: list[dict[str, Any]] = []
    for index, raw in enumerate(raw_entries):
        entry = exact_keys(
            raw, V19_TOPOLOGY_CPU_KEYS, f"v19 topology CPU {index}")
        cpu = v19_bounded_int(
            entry["cpu"], 0, MAX_CPU_ID, "v19 topology CPU")
        require(type(entry["online"]) is bool,
                "v19 topology online flag differs")
        siblings = v19_cpu_list(
            entry["thread_siblings"], "v19 thread siblings")
        domain = v19_cpu_list(entry["domain_cpus"], "v19 domain CPUs")
        require(len(siblings) in (1, 2) and cpu in siblings and cpu in domain,
                "v19 topology membership differs")
        entries.append({
            "cpu": cpu, "online": entry["online"],
            "physical_package_id": v19_bounded_int(
                entry["physical_package_id"], 0, (1 << 31) - 1,
                "v19 package ID"),
            "die_id": v19_bounded_int(
                entry["die_id"], 0, (1 << 31) - 1, "v19 die ID"),
            "core_id": v19_bounded_int(
                entry["core_id"], 0, (1 << 31) - 1, "v19 core ID"),
            "thread_siblings": siblings, "domain_cpus": domain,
        })
    require([entry["cpu"] for entry in entries] ==
            sorted({entry["cpu"] for entry in entries}),
            "v19 topology CPUs are unordered or duplicated")
    by_cpu = {entry["cpu"]: entry for entry in entries}
    domain_ids: dict[tuple[int, ...], int] = {}
    domain_id_by_cpu: dict[int, int] = {}
    for entry in entries:
        signature = tuple(entry["domain_cpus"])
        domain_id_by_cpu[entry["cpu"]] = domain_ids.setdefault(
            signature, len(domain_ids))
    physical_groups: dict[tuple[int, int, int], list[int]] = {}
    for entry in entries:
        physical = (entry["physical_package_id"], entry["die_id"],
                    entry["core_id"])
        prior = physical_groups.setdefault(
            physical, entry["thread_siblings"])
        require(prior == entry["thread_siblings"],
                "v19 physical core names multiple sibling groups")
        for sibling_cpu in entry["thread_siblings"]:
            require(sibling_cpu in by_cpu,
                    "v19 topology names a missing sibling")
            sibling = by_cpu[sibling_cpu]
            require(sibling["thread_siblings"] == entry["thread_siblings"] and
                    (sibling["physical_package_id"], sibling["die_id"],
                     sibling["core_id"]) == physical and
                    domain_id_by_cpu[sibling_cpu] ==
                        domain_id_by_cpu[entry["cpu"]],
                    "v19 topology sibling symmetry differs")
        for domain_cpu in entry["domain_cpus"]:
            require(domain_cpu in by_cpu,
                    "v19 topology names a missing domain CPU")
            member = by_cpu[domain_cpu]
            require(domain_id_by_cpu[domain_cpu] ==
                        domain_id_by_cpu[entry["cpu"]] and
                    member["physical_package_id"] ==
                        entry["physical_package_id"] and
                    member["die_id"] == entry["die_id"],
                    "v19 topology domain symmetry differs")
    expected = {"schema": V19_TOPOLOGY_SCHEMA, "cpus": entries}
    require(v19_exact_json_equal(topology, expected),
            "v19 topology is not canonical")
    return expected


def v19_candidate_pairs(
    policy: Mapping[str, Any], topology: Mapping[str, Any],
    allowed: Sequence[int],
) -> list[dict[str, Any]]:
    by_cpu = {entry["cpu"]: entry for entry in topology["cpus"]}
    allowed_set = set(allowed)
    for cpu in allowed:
        require(cpu in by_cpu and by_cpu[cpu]["online"],
                "v19 launch affinity contains a missing or offline CPU")
    candidates: list[dict[str, Any]] = []
    physical_pairs: set[frozenset[int]] = set()
    excluded_by_primary = {
        exclusion["benchmark_cpu"]: exclusion
        for exclusion in policy["excluded_pairs"]
    }
    for primary in policy["candidate_primary_cpus"]:
        require(primary in by_cpu,
                "v19 candidate primary is missing from topology")
        entry = by_cpu[primary]
        require(len(entry["thread_siblings"]) == 2,
                "v19 candidate is not exact SMT2")
        sibling = next(
            cpu for cpu in entry["thread_siblings"] if cpu != primary)
        physical = frozenset((primary, sibling))
        require(physical not in physical_pairs,
                "v19 policy repeats a physical core")
        physical_pairs.add(physical)
        exclusion = excluded_by_primary.get(primary)
        if exclusion is not None:
            require(exclusion["reserved_sibling"] == sibling,
                    "v19 excluded pair differs from topology")
            continue
        if not (entry["online"] and by_cpu[sibling]["online"] and
                primary in allowed_set and sibling in allowed_set):
            continue
        candidates.append({
            "benchmark_cpu": primary, "reserved_sibling": sibling,
            "domain_cpus": list(entry["domain_cpus"]),
        })
    return candidates


def v19_tick_cap(elapsed_ns: int, ticks: int) -> int:
    return (elapsed_ns * ticks) // 1_000_000_000 + 1


def v19_window(
    value: Any, *, policy: Mapping[str, Any], phase: str, index: int,
    observed_cpus: Sequence[int],
) -> dict[str, Any]:
    window = exact_keys(value, V19_WINDOW_KEYS, "v19 observation window")
    require(phase in ("qualification", "bridge") and
            type(index) is int and 0 <= index < V19_MAX_WINDOWS,
            "v19 observation window identity differs")
    cpus = list(observed_cpus)
    before = v19_snapshot(window["before"], cpus, "v19 window before")
    after = v19_snapshot(window["after"], cpus, "v19 window after")
    narrow = (after["read_started_monotonic_ns"] -
              before["read_finished_monotonic_ns"])
    wide = (after["read_finished_monotonic_ns"] -
            before["read_started_monotonic_ns"])
    require(narrow > 0 and wide >= narrow,
            "v19 observation interval differs")
    ticks = policy["clock_ticks_per_second"]
    narrow_cap = v19_tick_cap(narrow, ticks)
    wide_cap = v19_tick_cap(wide, ticks)
    lower = narrow_cap - V19_LIVENESS_LOWER_SLACK_JIFFIES
    deltas: dict[str, Any] = {}
    nonidle: list[int] = []
    not_live: list[int] = []
    for cpu in cpus:
        delta = v19_counter_delta(
            before["cpus"][str(cpu)], after["cpus"][str(cpu)],
            f"v19 window CPU {cpu}")
        deltas[str(cpu)] = delta
        if delta["nonidle_jiffies"] > V19_NONIDLE_LIMIT_JIFFIES:
            nonidle.append(cpu)
        if not (lower >= 1 and lower <= delta["total_jiffies"] <= wide_cap):
            not_live.append(cpu)
    expected = {
        "schema": V19_WINDOW_SCHEMA, "phase": phase, "index": index,
        "before": before, "after": after, "window_ns": narrow,
        "wide_window_ns": wide, "elapsed_tick_cap": narrow_cap,
        "wide_elapsed_tick_cap": wide_cap,
        "liveness_lower_bound_jiffies": lower, "delta": deltas,
        "nonidle_cpus": nonidle, "not_live_cpus": not_live,
        "clock_ticks_per_second": ticks,
        "policy": v19_observation_policy(),
    }
    require(v19_exact_json_equal(window, expected),
            "v19 observation window differs from fixed-point replay")
    return expected


def v19_scan(
    value: Any, policy: Mapping[str, Any], frozen_value: Any,
) -> dict[str, Any]:
    scan = exact_keys(value, V19_SCAN_KEYS, "v19 qualification scan")
    require(scan["schema"] == V19_SCAN_SCHEMA and
            v19_exact_json_equal(scan["policy"], policy),
            "v19 scan policy or schema differs")
    frozen = v19_pair(
        frozen_value, "v19 frozen pair", allow_none=True)
    require(v19_exact_json_equal(
                scan["frozen_pair_from_prior_attempt"], frozen),
            "v19 scan frozen pair differs")
    topology_before = v19_topology(scan["topology_before"])
    topology_after = v19_topology(scan["topology_after"])
    require(v19_exact_json_equal(topology_before, topology_after),
            "v19 topology changed during qualification")
    allowed = v19_cpu_list(
        scan["allowed_cpu_set_at_launch"], "v19 launch affinity")
    candidates = v19_candidate_pairs(policy, topology_before, allowed)
    observed = sorted({
        cpu for candidate in candidates for cpu in candidate["domain_cpus"]
    })
    raw_windows = scan["windows"]
    require(type(raw_windows) is list and
            0 < len(raw_windows) <= V19_MAX_WINDOWS,
            "v19 qualification window list differs")
    windows: list[dict[str, Any]] = []
    for index, raw_window in enumerate(raw_windows):
        window = v19_window(
            raw_window, policy=policy, phase="qualification", index=index,
            observed_cpus=observed)
        if windows:
            require(v19_exact_json_equal(
                        windows[-1]["after"], window["before"]),
                    "v19 qualification windows do not share endpoints")
        windows.append(window)
    first, last = windows[0]["before"], windows[-1]["after"]
    narrow = (last["read_started_monotonic_ns"] -
              first["read_finished_monotonic_ns"])
    wide = (last["read_finished_monotonic_ns"] -
            first["read_started_monotonic_ns"])
    require(narrow > 0 and wide >= narrow,
            "v19 qualification aggregate interval differs")
    ticks = policy["clock_ticks_per_second"]
    overall_cap = v19_tick_cap(narrow, ticks)
    overall_wide_cap = v19_tick_cap(wide, ticks)
    overall_lower = overall_cap - V19_LIVENESS_LOWER_SLACK_JIFFIES
    nonidle_sets = [set(window["nonidle_cpus"]) for window in windows]
    not_live_sets = [set(window["not_live_cpus"]) for window in windows]
    aggregates: dict[str, Any] = {}
    for cpu in observed:
        aggregate = v19_counter_delta(
            first["cpus"][str(cpu)], last["cpus"][str(cpu)],
            f"v19 scan aggregate CPU {cpu}")
        summed = {
            field: sum(window["delta"][str(cpu)]["fields"][field]
                       for window in windows)
            for field in V19_COUNTER_FIELDS
        }
        require(summed == aggregate["fields"],
                "v19 qualification windows do not telescope")
        aggregates[str(cpu)] = {
            **aggregate,
            "nonidle_window_count": sum(
                cpu in values for values in nonidle_sets),
            "not_live_window_count": sum(
                cpu in values for values in not_live_sets),
            "overall_live": (
                overall_lower >= 1 and
                overall_lower <= aggregate["total_jiffies"] <=
                    overall_wide_cap),
        }
    candidate_results: list[dict[str, Any]] = []
    eligible: list[dict[str, int]] = []
    for candidate in candidates:
        pair_cpus = (
            candidate["benchmark_cpu"], candidate["reserved_sibling"])
        domain = candidate["domain_cpus"]
        pair_idle = all(
            aggregates[str(cpu)]["nonidle_jiffies"] == 0 and
            aggregates[str(cpu)]["nonidle_window_count"] == 0
            for cpu in pair_cpus)
        pair_live = all(
            aggregates[str(cpu)]["overall_live"] and
            aggregates[str(cpu)]["not_live_window_count"] == 0
            for cpu in pair_cpus)
        domain_idle = all(
            aggregates[str(cpu)]["nonidle_jiffies"] == 0 and
            aggregates[str(cpu)]["nonidle_window_count"] == 0
            for cpu in domain)
        domain_live = all(
            aggregates[str(cpu)]["overall_live"] and
            aggregates[str(cpu)]["not_live_window_count"] == 0
            for cpu in domain)
        qualified = pair_idle and pair_live and (
            policy["domain_mode"] == "pair-only" or
            (domain_idle and domain_live))
        result = {
            **candidate, "pair_idle": pair_idle, "pair_live": pair_live,
            "domain_idle": domain_idle, "domain_live": domain_live,
            "qualified": qualified,
        }
        candidate_results.append(result)
        if qualified:
            eligible.append({
                "benchmark_cpu": candidate["benchmark_cpu"],
                "reserved_sibling": candidate["reserved_sibling"],
            })
    if frozen is not None:
        by_cpu = {entry["cpu"]: entry for entry in topology_before["cpus"]}
        static_pairs = []
        excluded_primaries = {
            exclusion["benchmark_cpu"]
            for exclusion in policy["excluded_pairs"]
        }
        for primary in policy["candidate_primary_cpus"]:
            entry = by_cpu[primary]
            sibling = next(
                cpu for cpu in entry["thread_siblings"] if cpu != primary)
            if primary not in excluded_primaries:
                static_pairs.append({
                    "benchmark_cpu": primary,
                    "reserved_sibling": sibling})
        require(frozen in static_pairs,
                "v19 frozen pair is not a static candidate")
        if frozen in eligible:
            selected = dict(frozen)
            selection_status = "selected-frozen-pair"
        else:
            selected = None
            selection_status = "frozen-pair-did-not-requalify"
    elif eligible:
        selected = dict(eligible[0])
        selection_status = "selected-lowest-primary"
    else:
        selected = None
        selection_status = "no-candidate-pair-qualified"
    expected = {
        "schema": V19_SCAN_SCHEMA, "policy": dict(policy),
        "allowed_cpu_set_at_launch": allowed,
        "topology_before": topology_before, "topology_after": topology_after,
        "topology_reread_identical": True, "windows": windows,
        "cpu_aggregates": aggregates, "candidate_pairs": candidate_results,
        "eligible_pairs": eligible, "selected": selected,
        "frozen_pair_from_prior_attempt": frozen,
        "selection_status": selection_status,
        "scan_started_monotonic_ns": first["read_started_monotonic_ns"],
        "scan_finished_monotonic_ns": last["read_finished_monotonic_ns"],
        "realized_scan_total_ns": narrow, "wide_scan_total_ns": wide,
        "overall_elapsed_tick_cap": overall_cap,
        "overall_wide_elapsed_tick_cap": overall_wide_cap,
        "overall_liveness_lower_bound_jiffies": overall_lower,
        "scan_window_count": len(windows),
        "max_proc_stat_read_ns": max(
            snapshot["read_finished_monotonic_ns"] -
            snapshot["read_started_monotonic_ns"]
            for snapshot in (
                [windows[0]["before"]] +
                [window["after"] for window in windows])),
        "excluded_pair_count": len(policy["excluded_pairs"]),
        "candidate_pair_count": len(candidate_results),
        "eligible_pair_count": len(eligible),
        "candidate_timing_performed": False,
    }
    require(v19_exact_json_equal(scan, expected),
            "v19 qualification scan differs from fixed-point replay")
    return expected


def v19_acquisition_geometry() -> dict[str, int]:
    return {
        "requested_window_count": V19_QUALIFICATION_WINDOW_COUNT,
        "nominal_window_ns": V19_QUALIFICATION_NOMINAL_WINDOW_NS,
    }


def v19_qualification_geometry() -> dict[str, int]:
    return {
        "window_count": V19_QUALIFICATION_WINDOW_COUNT,
        "nominal_window_ns": V19_QUALIFICATION_NOMINAL_WINDOW_NS,
    }


def v19_bridge_geometry() -> dict[str, int]:
    return {
        "minimum_window_count": V19_BRIDGE_WINDOW_COUNT,
        "maximum_window_count": V19_BRIDGE_WINDOW_COUNT,
        "nominal_window_ns": V19_BRIDGE_NOMINAL_WINDOW_NS,
        "maximum_handoff_elapsed_ns": V19_MAXIMUM_HANDOFF_ELAPSED_NS,
    }


def validate_v19_bridge_geometry(value: Any) -> dict[str, int]:
    geometry = exact_keys(
        value, V19_BRIDGE_GEOMETRY_KEYS, "v19 bridge geometry")
    minimum = v19_bounded_int(
        geometry["minimum_window_count"], 1, V19_MAX_WINDOWS,
        "v19 minimum bridge window count")
    maximum = v19_bounded_int(
        geometry["maximum_window_count"], minimum, V19_MAX_WINDOWS,
        "v19 maximum bridge window count")
    nominal = v19_bounded_int(
        geometry["nominal_window_ns"], 1, 86_400_000_000_000,
        "v19 nominal bridge duration")
    handoff = v19_bounded_int(
        geometry["maximum_handoff_elapsed_ns"], 1,
        V19_MAX_MONOTONIC_NS, "v19 maximum bridge handoff duration")
    require(minimum * nominal <= handoff,
            "v19 bridge geometry cannot contain its minimum windows")
    expected = v19_bridge_geometry()
    require(v19_exact_json_equal(geometry, expected),
            "v19 bridge geometry differs from preregistration")
    return expected


def v19_acquisition(
    value: Any, policy: Mapping[str, Any], frozen_value: Any,
) -> dict[str, Any]:
    acquisition = exact_keys(
        value, V19_ACQUISITION_KEYS, "v19 qualification acquisition")
    require(len(v19_canonical_bytes(acquisition)) <= V19_MAX_JSON_BYTES,
            "v19 acquisition exceeds its canonical byte bound")
    frozen = v19_pair(
        frozen_value, "v19 acquisition frozen pair", allow_none=True)
    require(acquisition["schema"] == V19_ACQUISITION_SCHEMA and
            acquisition["acquisition_method"] == V19_ACQUISITION_METHOD and
            type(acquisition["sources"]) is dict and
            v19_exact_json_equal(
                acquisition["sources"], dict(V19_SOURCE_ITEMS)) and
            v19_exact_json_equal(acquisition["policy"], policy) and
            acquisition["policy_sha256"] == v19_sha256(policy) and
            type(acquisition["requested_window_count"]) is int and
            acquisition["requested_window_count"] ==
                V19_QUALIFICATION_WINDOW_COUNT and
            type(acquisition["nominal_window_ns"]) is int and
            acquisition["nominal_window_ns"] ==
                V19_QUALIFICATION_NOMINAL_WINDOW_NS and
            v19_exact_json_equal(
                acquisition["frozen_pair_from_prior_attempt"], frozen),
            "v19 acquisition identity differs")
    allowed_launch = v19_cpu_list(
        acquisition["allowed_cpu_set_at_launch"],
        "v19 acquisition launch affinity")
    allowed_after = v19_cpu_list(
        acquisition["allowed_cpu_set_after_scan"],
        "v19 acquisition final affinity")
    require(allowed_launch == allowed_after,
            "v19 launch affinity changed during acquisition")
    ticks_launch = v19_bounded_int(
        acquisition["clock_ticks_per_second_at_launch"], 1, 1_000_000,
        "v19 acquisition launch clock rate")
    ticks_after = v19_bounded_int(
        acquisition["clock_ticks_per_second_after_scan"], 1, 1_000_000,
        "v19 acquisition final clock rate")
    require(ticks_launch == ticks_after == policy["clock_ticks_per_second"],
            "v19 acquisition clock rate changed or differs from policy")
    scan = v19_scan(acquisition["scan"], policy, frozen)
    require(scan["scan_window_count"] == V19_QUALIFICATION_WINDOW_COUNT and
            scan["allowed_cpu_set_at_launch"] == allowed_launch and
            all(window["window_ns"] >=
                V19_QUALIFICATION_NOMINAL_WINDOW_NS
                for window in scan["windows"]),
            "v19 acquisition scan geometry or affinity differs")
    before_hash = v19_sha256(scan["topology_before"])
    after_hash = v19_sha256(scan["topology_after"])
    scan_hash = v19_sha256(scan)
    require(acquisition["topology_before_sha256"] == before_hash and
            acquisition["topology_after_sha256"] == after_hash and
            acquisition["scan_sha256"] == scan_hash and
            type(acquisition["host_mutation_performed"]) is bool and
            acquisition["host_mutation_performed"] is False and
            type(acquisition["candidate_timing_performed"]) is bool and
            acquisition["candidate_timing_performed"] is False and
            scan["candidate_timing_performed"] is False,
            "v19 acquisition hashes, mutation, or timing flags differ")
    expected = {
        "schema": V19_ACQUISITION_SCHEMA,
        "acquisition_method": V19_ACQUISITION_METHOD,
        "sources": dict(V19_SOURCE_ITEMS), "policy": dict(policy),
        "policy_sha256": v19_sha256(policy),
        "requested_window_count": V19_QUALIFICATION_WINDOW_COUNT,
        "nominal_window_ns": V19_QUALIFICATION_NOMINAL_WINDOW_NS,
        "frozen_pair_from_prior_attempt": frozen,
        "allowed_cpu_set_at_launch": allowed_launch,
        "allowed_cpu_set_after_scan": allowed_after,
        "clock_ticks_per_second_at_launch": ticks_launch,
        "clock_ticks_per_second_after_scan": ticks_after,
        "topology_before_sha256": before_hash,
        "topology_after_sha256": after_hash, "scan": scan,
        "scan_sha256": scan_hash, "host_mutation_performed": False,
        "candidate_timing_performed": False,
    }
    require(v19_exact_json_equal(acquisition, expected),
            "v19 acquisition differs from fixed-point replay")
    return expected


def v19_bridge_from_validated_acquisition(
    value: Any, acquisition: Mapping[str, Any],
    policy: Mapping[str, Any], frozen_value: Any,
) -> dict[str, Any]:
    bridge = exact_keys(value, V19_BRIDGE_KEYS, "v19 qualification bridge")
    require(len(v19_canonical_bytes(bridge)) <= V19_MAX_JSON_BYTES,
            "v19 bridge exceeds its canonical byte bound")
    frozen = v19_pair(
        frozen_value, "v19 bridge frozen pair", allow_none=True)
    geometry = validate_v19_bridge_geometry(bridge["bridge_geometry"])
    require(bridge["schema"] == V19_BRIDGE_SCHEMA and
            bridge["handoff_rule"] == V19_HANDOFF_RULE and
            bridge["acceptance_rule"] == V19_BRIDGE_ACCEPTANCE_RULE and
            v19_exact_json_equal(
                bridge["frozen_pair_from_prior_attempt"], frozen),
            "v19 bridge identity differs")
    scan = acquisition["scan"]
    selected = v19_pair(scan["selected"], "v19 bridge selected pair")
    require(selected is not None and scan["selection_status"] in (
                "selected-lowest-primary", "selected-frozen-pair"),
            "v19 bridge lacks a selected acquisition")
    scan_tail = scan["windows"][-1]["after"]
    observed = sorted(int(cpu) for cpu in scan_tail["cpus"])
    require(observed and len(observed) <= MAX_CPU_LIST_ENTRIES,
            "v19 bridge observed CPU universe differs")
    raw_windows = bridge["windows"]
    require(type(raw_windows) is list and
            geometry["minimum_window_count"] <= len(raw_windows) <=
                geometry["maximum_window_count"],
            "v19 bridge window count differs")
    windows: list[dict[str, Any]] = []
    for index, raw_window in enumerate(raw_windows):
        window = v19_window(
            raw_window, policy=policy, phase="bridge", index=index,
            observed_cpus=observed)
        if index == 0:
            require(v19_exact_json_equal(window["before"], scan_tail),
                    "v19 bridge head differs from qualification tail")
        else:
            require(v19_exact_json_equal(
                        windows[-1]["after"], window["before"]),
                    "v19 bridge windows do not share endpoints")
        windows.append(window)
    require(all(window["window_ns"] >= V19_BRIDGE_NOMINAL_WINDOW_NS
                for window in windows),
            "v19 bridge contains a short observation window")
    first, last = windows[0]["before"], windows[-1]["after"]
    started = first["read_finished_monotonic_ns"]
    finished = last["read_finished_monotonic_ns"]
    require(started <= V19_MAX_MONOTONIC_NS -
                V19_MAXIMUM_HANDOFF_ELAPSED_NS,
            "v19 bridge deadline overflows monotonic time")
    deadline = started + V19_MAXIMUM_HANDOFF_ELAPSED_NS
    require(finished <= deadline,
            "v19 bridge closes after its fixed deadline")
    narrow = last["read_started_monotonic_ns"] - started
    wide = (last["read_finished_monotonic_ns"] -
            first["read_started_monotonic_ns"])
    require(narrow > 0 and wide >= narrow,
            "v19 bridge aggregate interval differs")
    ticks = policy["clock_ticks_per_second"]
    overall_cap = v19_tick_cap(narrow, ticks)
    overall_wide_cap = v19_tick_cap(wide, ticks)
    overall_lower = overall_cap - V19_LIVENESS_LOWER_SLACK_JIFFIES
    nonidle_sets = [set(window["nonidle_cpus"]) for window in windows]
    not_live_sets = [set(window["not_live_cpus"]) for window in windows]
    aggregates: dict[str, Any] = {}
    for cpu in observed:
        aggregate = v19_counter_delta(
            first["cpus"][str(cpu)], last["cpus"][str(cpu)],
            f"v19 bridge aggregate CPU {cpu}")
        summed = {
            field: sum(window["delta"][str(cpu)]["fields"][field]
                       for window in windows)
            for field in V19_COUNTER_FIELDS
        }
        require(summed == aggregate["fields"],
                "v19 bridge windows do not telescope")
        aggregates[str(cpu)] = {
            **aggregate,
            "nonidle_window_count": sum(
                cpu in values for values in nonidle_sets),
            "not_live_window_count": sum(
                cpu in values for values in not_live_sets),
            "overall_live": (
                overall_lower >= 1 and
                overall_lower <= aggregate["total_jiffies"] <=
                    overall_wide_cap),
        }
    candidates = [
        candidate for candidate in scan["candidate_pairs"]
        if candidate["benchmark_cpu"] == selected["benchmark_cpu"] and
        candidate["reserved_sibling"] == selected["reserved_sibling"]
    ]
    require(len(candidates) == 1 and candidates[0]["qualified"] is True,
            "v19 bridge pair was not independently qualified")
    guarded = (
        list(candidates[0]["domain_cpus"])
        if policy["domain_mode"] == "pair-and-domain"
        else sorted(selected.values()))
    require(set(guarded) <= set(observed),
            "v19 bridge guard escapes observed CPUs")
    nonidle_guarded = [
        cpu for cpu in guarded
        if aggregates[str(cpu)]["nonidle_jiffies"] != 0 or
        aggregates[str(cpu)]["nonidle_window_count"] != 0
    ]
    not_live_guarded = [
        cpu for cpu in guarded
        if not aggregates[str(cpu)]["overall_live"] or
        aggregates[str(cpu)]["not_live_window_count"] != 0
    ]
    acquisition_geometry = v19_acquisition_geometry()
    scan_tail_hash = v19_sha256(scan_tail)
    bridge_tail_hash = v19_sha256(last)
    expected = {
        "schema": V19_BRIDGE_SCHEMA, "handoff_rule": V19_HANDOFF_RULE,
        "acceptance_rule": V19_BRIDGE_ACCEPTANCE_RULE,
        "acquisition_sha256": v19_sha256(acquisition),
        "scan_sha256": v19_sha256(scan),
        "policy_sha256": v19_sha256(policy),
        "selected_pair": selected,
        "frozen_pair_from_prior_attempt": frozen,
        "acquisition_geometry": acquisition_geometry,
        "acquisition_geometry_sha256": v19_sha256(acquisition_geometry),
        "bridge_geometry": geometry,
        "bridge_geometry_sha256": v19_sha256(geometry),
        "observed_cpus": observed, "guarded_cpus": guarded,
        "scan_tail_sha256": scan_tail_hash, "windows": windows,
        "cpu_aggregates": aggregates,
        "bridge_head_sha256": v19_sha256(first),
        "bridge_tail_sha256": bridge_tail_hash,
        "campaign_presample_before": last,
        "campaign_presample_before_sha256": bridge_tail_hash,
        "bridge_started_monotonic_ns": started,
        "bridge_finished_monotonic_ns": finished,
        "bridge_deadline_monotonic_ns": deadline,
        "realized_bridge_total_ns": narrow,
        "wide_bridge_total_ns": wide,
        "overall_elapsed_tick_cap": overall_cap,
        "overall_wide_elapsed_tick_cap": overall_wide_cap,
        "overall_liveness_lower_bound_jiffies": overall_lower,
        "bridge_window_count": len(windows),
        "nonidle_guarded_cpus": nonidle_guarded,
        "not_live_guarded_cpus": not_live_guarded,
        "bridge_accepted": not nonidle_guarded and not not_live_guarded,
        "host_mutation_performed": False,
        "candidate_timing_performed": False,
        "shared_host_claim_ceiling": v19_claim_ceiling(),
    }
    require(v19_exact_json_equal(bridge, expected),
            "v19 bridge differs from fixed-point replay")
    return expected


def v19_first_window_before(
    value: Any, selected: Mapping[str, int],
) -> dict[str, Any]:
    before = exact_keys(
        value, V19_FIRST_WINDOW_KEYS, "v19 first-window before snapshot")
    started = v19_bounded_int(
        before["read_started_monotonic_ns"], 0, V19_MAX_MONOTONIC_NS,
        "v19 first-window read start")
    finished = v19_bounded_int(
        before["read_finished_monotonic_ns"], started,
        V19_MAX_MONOTONIC_NS, "v19 first-window read finish")
    require(type(before["monotonic_ns"]) is int and
            before["monotonic_ns"] == finished,
            "v19 first-window monotonic boundary differs")
    expected = {
        "benchmark_cpu": v19_counter(
            before["benchmark_cpu"], selected["benchmark_cpu"],
            "v19 first-window benchmark CPU"),
        "monotonic_ns": finished,
        "read_finished_monotonic_ns": finished,
        "read_started_monotonic_ns": started,
        "reserved_sibling": v19_counter(
            before["reserved_sibling"], selected["reserved_sibling"],
            "v19 first-window reserved sibling"),
    }
    require(v19_exact_json_equal(before, expected),
            "v19 first-window before snapshot differs")
    return expected


def v19_nonidle_delta(
    before: Mapping[str, Any], after: Mapping[str, Any], label: str,
) -> int:
    total = 0
    for field in V19_COUNTER_FIELDS:
        first = before["fields"][field]
        last = after["fields"][field]
        require(last >= first, f"{label} {field} counter regressed")
        if field in V19_NONIDLE_FIELDS:
            total += last - first
    return v19_bounded_int(total, 0, V19_MAX_COUNTER, f"{label} nonidle")


def v19_handoff(
    value: Any, bridge: Mapping[str, Any], selected: Mapping[str, int],
) -> dict[str, Any]:
    handoff = exact_keys(value, V19_HANDOFF_KEYS, "v19 first-window handoff")
    require(bridge["bridge_accepted"] is True and
            v19_exact_json_equal(bridge["selected_pair"], selected) and
            v19_exact_json_equal(handoff["selected_pair"], selected),
            "v19 handoff lacks its accepted selected bridge")
    tail = bridge["campaign_presample_before"]
    tail_hash = v19_sha256(tail)
    require(bridge["bridge_tail_sha256"] == tail_hash and
            bridge["campaign_presample_before_sha256"] == tail_hash,
            "v19 handoff bridge-tail binding differs")
    before = v19_first_window_before(
        handoff["first_window_before"], selected)
    tail_finished = v19_bounded_int(
        tail["read_finished_monotonic_ns"], 0, V19_MAX_MONOTONIC_NS,
        "v19 bridge-tail read finish")
    first_started = before["read_started_monotonic_ns"]
    require(first_started >= tail_finished,
            "v19 first-window starts before the bridge tail")
    elapsed = first_started - tail_finished
    deltas: dict[str, int] = {}
    for role in ("benchmark_cpu", "reserved_sibling"):
        cpu = selected[role]
        require(str(cpu) in tail["cpus"],
                "v19 bridge tail omits the selected pair")
        tail_counter = v19_counter(
            tail["cpus"][str(cpu)], cpu, f"v19 bridge-tail {role}")
        deltas[role] = v19_nonidle_delta(
            tail_counter, before[role], f"v19 first-window {role}")
    late = elapsed > V19_MAXIMUM_HANDOFF_ELAPSED_NS
    nonidle = any(delta != 0 for delta in deltas.values())
    failure = (
        "first-window-handoff-late" if late else
        "first-window-handoff-selected-pair-nonidle" if nonidle else None)
    expected = {
        "schema": V19_HANDOFF_SCHEMA, "selected_pair": dict(selected),
        "bridge_tail_sha256": tail_hash,
        "bridge_tail_read_finished_monotonic_ns": tail_finished,
        "first_window_before": before,
        "first_window_before_sha256": v19_sha256(before),
        "first_window_before_read_started_monotonic_ns": first_started,
        "handoff_elapsed_ns": elapsed,
        "maximum_handoff_elapsed_ns": V19_MAXIMUM_HANDOFF_ELAPSED_NS,
        "selected_pair_nonidle_delta": deltas,
        "accepted": not late and not nonidle,
        "failure_terminal": failure,
    }
    require(v19_exact_json_equal(handoff, expected),
            "v19 handoff differs from fixed-point replay")
    return expected


def validate_v19_qualification(
    value: Any, expected_attempt: Any, *,
    expected_host_identity_sha256: Any,
) -> tuple[dict[str, Any], tuple[int, int]]:
    record = exact_keys(
        value, V19_RECORD_KEYS, "pair-qualified v19 record")
    require(len(v19_canonical_bytes(record)) <= V19_MAX_JSON_BYTES,
            "pair-qualified v19 record exceeds its canonical byte bound")
    require(record["schema"] == V19_RECORD_SCHEMA and
            record["stage"] == "complete" and
            record["record_status"] == "complete" and
            record["terminal"] == V19_SUCCESS_TERMINAL,
            "auditor requires complete NOT_PROMOTED v19 evidence")
    policy = validate_v19_policy(record["policy"])
    policy_hash = v19_sha256(policy)
    require(record["policy_sha256"] == policy_hash,
            "v19 qualification policy hash differs")
    host_hash = v19_hex64(
        record["host_identity_sha256"], "v19 host identity hash")
    require(host_hash == v19_hex64(
                expected_host_identity_sha256,
                "expected v19 host identity hash"),
            "v19 qualification host identity differs")
    attempt = validate_v19_attempt(record["attempt"], expected_attempt)
    geometry = exact_keys(
        record["qualification_geometry"],
        V19_QUALIFICATION_GEOMETRY_KEYS, "v19 qualification geometry")
    expected_geometry = v19_qualification_geometry()
    require(v19_exact_json_equal(geometry, expected_geometry) and
            record["qualification_geometry_sha256"] ==
                v19_sha256(expected_geometry),
            "v19 qualification geometry differs")
    acquisition = v19_acquisition(
        record["acquisition"], policy,
        attempt["frozen_pair_from_prior_attempt"])
    selected = v19_pair(
        acquisition["scan"]["selected"], "v19 selected pair")
    require(selected is not None and
            v19_exact_json_equal(record["selected_pair"], selected) and
            record["selection_status"] ==
                acquisition["scan"]["selection_status"] and
            record["acquisition_sha256"] == v19_sha256(acquisition),
            "v19 selected pair, acquisition hash, or status differs")
    bridge = v19_bridge_from_validated_acquisition(
        record["bridge"], acquisition, policy,
        attempt["frozen_pair_from_prior_attempt"])
    require(bridge["bridge_accepted"] is True and
            record["bridge_sha256"] == v19_sha256(bridge) and
            record["bridge_geometry_sha256"] ==
                v19_sha256(v19_bridge_geometry()),
            "v19 record lacks an accepted bridge fixed point")
    handoff = v19_handoff(
        record["first_window_handoff"], bridge, selected)
    require(handoff["accepted"] is True,
            "v19 record lacks an accepted first-window handoff")
    lineage = validate_v19_lineage(record["v18_failure_lineage"])
    claims = exact_keys(
        record["shared_host_claim_ceiling"], V19_CLAIM_KEYS,
        "v19 shared-host claim ceiling")
    require(v19_exact_json_equal(claims, v19_claim_ceiling()) and
            record["v18_failure_lineage_sha256"] == v19_sha256(lineage) and
            type(record["host_mutation_performed"]) is bool and
            record["host_mutation_performed"] is False and
            type(record["candidate_timing_performed"]) is bool and
            record["candidate_timing_performed"] is False,
            "v19 claim ceiling, lineage hash, mutation, or timing differs")
    expected = {
        "schema": V19_RECORD_SCHEMA, "stage": "complete",
        "record_status": "complete", "attempt": attempt,
        "policy": policy, "policy_sha256": policy_hash,
        "host_identity_sha256": host_hash,
        "qualification_geometry": expected_geometry,
        "qualification_geometry_sha256": v19_sha256(expected_geometry),
        "acquisition": acquisition,
        "acquisition_sha256": v19_sha256(acquisition),
        "selected_pair": selected,
        "selection_status": acquisition["scan"]["selection_status"],
        "bridge": bridge, "bridge_sha256": v19_sha256(bridge),
        "bridge_geometry_sha256": v19_sha256(v19_bridge_geometry()),
        "first_window_handoff": handoff,
        "v18_failure_lineage": lineage,
        "v18_failure_lineage_sha256": v19_sha256(lineage),
        "shared_host_claim_ceiling": v19_claim_ceiling(),
        "terminal": V19_SUCCESS_TERMINAL,
        "host_mutation_performed": False,
        "candidate_timing_performed": False,
    }
    require(v19_exact_json_equal(record, expected),
            "v19 record differs from full fixed-point replay")
    return expected, (
        selected["benchmark_cpu"], selected["reserved_sibling"])


def bounded_read(path: Path, maximum: int) -> bytes:
    with path.open("rb", buffering=0) as source:
        metadata = os.fstat(source.fileno())
        if stat.S_ISREG(metadata.st_mode) and metadata.st_size > 0:
            require(metadata.st_size <= maximum,
                    f"oversized record: {path}")
            data = source.read(metadata.st_size + 1)
            require(len(data) == metadata.st_size and source.read(1) == b"",
                    f"record changed while reading: {path}")
            return data
        chunks: list[bytes] = []
        total = 0
        while total <= maximum:
            chunk = source.read(min(1 << 20, maximum + 1 - total))
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
    require(total <= maximum, f"oversized record: {path}")
    return b"".join(chunks)


def parse_cpu_list(text: str) -> list[int]:
    require(text != "" and text == text.strip(), "invalid empty/padded CPU list")
    values: list[int] = []
    previous = -1
    for token in text.split(","):
        require(token != "", "empty CPU-list token")
        if "-" in token:
            require(token.count("-") == 1, "invalid CPU range")
            first_text, last_text = token.split("-", 1)
        else:
            first_text = last_text = token
        require(first_text.isascii() and last_text.isascii() and
                first_text.isdecimal() and last_text.isdecimal() and
                len(first_text) <= 7 and len(last_text) <= 7,
                "non-decimal CPU-list token")
        require((first_text == "0" or not first_text.startswith("0")) and
                (last_text == "0" or not last_text.startswith("0")),
                "noncanonical CPU-list integer")
        first, last = int(first_text), int(last_text)
        require(0 <= first <= last <= MAX_CPU_ID and first > previous,
                "overlapping, unordered, or out-of-range CPU list")
        count = last - first + 1
        require(count <= MAX_CPU_LIST_ENTRIES and
                len(values) + count <= MAX_CPU_LIST_ENTRIES,
                "CPU list has too many entries")
        values.extend(range(first, last + 1))
        previous = last
    require(values and format_cpu_list(values) == text,
            "CPU list is empty or noncanonical")
    return values


def parse_oriented_cpu_pair(text: str, label: str) -> dict[str, int]:
    """Parse a bounded canonical BENCHMARK,SIBLING command-line pair."""
    require(type(text) is str, f"{label} must be BENCHMARK,SIBLING")
    components = text.split(",")
    require(len(components) == 2 and
            all(component and component.isascii() and
                component.isdecimal() and len(component) <= 7 and
                (component == "0" or not component.startswith("0"))
                for component in components),
            f"{label} must be BENCHMARK,SIBLING")
    return v19_pair({
        "benchmark_cpu": int(components[0], 10),
        "reserved_sibling": int(components[1], 10),
    }, label)


def format_cpu_list(values: Sequence[int]) -> str:
    ordered = sorted(set(values))
    require(ordered and all(type(value) is int and 0 <= value <= MAX_CPU_ID
                            for value in ordered), "invalid CPU values")
    ranges: list[str] = []
    start = previous = ordered[0]
    for value in ordered[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(ranges)


def parse_status(data: bytes) -> dict[str, Any]:
    text = data.decode("utf-8", "strict")
    require("\x00" not in text, "NUL in proc status")
    fields: dict[str, str] = {}
    for line in text.splitlines():
        if ":" not in line:
            continue
        name, value = line.split(":", 1)
        if name in ("Uid", "Cpus_allowed_list", "Name"):
            require(name not in fields, f"duplicate proc status field: {name}")
            fields[name] = value.strip()
    require(set(fields) == {"Uid", "Cpus_allowed_list", "Name"},
            "proc status lacks required fields")
    uid_tokens = fields["Uid"].split()
    require(len(uid_tokens) == 4 and all(token.isdecimal()
                                        for token in uid_tokens),
            "invalid proc status Uid")
    cpus = parse_cpu_list(fields["Cpus_allowed_list"])
    require(format_cpu_list(cpus) == fields["Cpus_allowed_list"],
            "noncanonical Cpus_allowed_list")
    name_bytes = fields["Name"].encode("utf-8")
    return {
        "uid": int(uid_tokens[0]),
        "cpus_allowed": cpus,
        "comm_size": len(name_bytes),
        "comm_sha256": hashlib.sha256(name_bytes).hexdigest(),
    }


def parse_starttime(data: bytes) -> int:
    require(b"\x00" not in data and data.endswith(b"\n"),
            "invalid proc stat framing")
    close = data.rfind(b")")
    require(close > 1 and data[0:data.find(b" ")].isdigit(),
            "invalid proc stat identity")
    tail = data[close + 2:].split()
    # tail[0] is field 3 (state), so field 22 (starttime) is tail[19].
    require(len(tail) >= 20 and tail[19].isdigit(),
            "proc stat lacks starttime")
    return int(tail[19])


def identity(path: Path, stat_path: Path) -> dict[str, int]:
    node = path.stat(follow_symlinks=False)
    require(stat.S_ISDIR(node.st_mode), f"proc task is not a directory: {path}")
    return {
        "device": node.st_dev,
        "inode": node.st_ino,
        "starttime_ticks": parse_starttime(bounded_read(stat_path, MAX_STAT_BYTES)),
    }


def namespace_identity(path: Path) -> dict[str, int]:
    link = path.stat(follow_symlinks=False)
    require(stat.S_ISLNK(link.st_mode),
            f"namespace path is not a symlink: {path}")
    # The procfs symlink belongs to the short-lived collector process.  Follow
    # it so endpoint snapshots bind the namespace object shared by both
    # collectors, rather than two necessarily different procfs dentries.
    node = path.stat(follow_symlinks=True)
    require(node.st_dev >= 0 and node.st_ino > 0,
            f"namespace target identity is invalid: {path}")
    return {"device": node.st_dev, "inode": node.st_ino}


def read_cpu_stat(cpus: Sequence[int]) -> dict[str, Any]:
    wanted = set(cpus)
    observed: dict[int, dict[str, Any]] = {}
    data = bounded_read(Path("/proc/stat"), 8 << 20)
    for line in data.splitlines():
        parts = line.split()
        if not parts or not parts[0].startswith(b"cpu") or parts[0] == b"cpu":
            continue
        suffix = parts[0][3:]
        if not suffix.isdigit() or int(suffix) not in wanted:
            continue
        cpu = int(suffix)
        require(cpu not in observed and 8 <= len(parts) - 1 <= len(CPU_FIELDS),
                "invalid or duplicate /proc/stat CPU record")
        values = []
        for token in parts[1:]:
            require(token.isdigit(), "non-integer /proc/stat counter")
            values.append(int(token))
        values.extend([0] * (len(CPU_FIELDS) - len(values)))
        fields = dict(zip(CPU_FIELDS, values))
        observed[cpu] = {
            "cpu": cpu,
            "fields": fields,
            "total_jiffies": sum(values[:8]),
            "nonidle_jiffies": sum(fields[name] for name in NONIDLE_FIELDS),
        }
    require(set(observed) == wanted, "reserved CPU missing from /proc/stat")
    return {str(cpu): observed[cpu] for cpu in sorted(observed)}


def vanished_record(pid: int, tid: int | None, stage: str) -> dict[str, Any]:
    return {"pid": pid, "tid": tid, "stage": stage}


def thread_entry(pid: int, tid: int, uid: int,
                 reserved: Sequence[int]) -> dict[str, Any]:
    process = Path(f"/proc/{pid}")
    task = process / "task" / str(tid)
    process_before = identity(process, process / "stat")
    task_before = identity(task, task / "stat")
    parsed = parse_status(bounded_read(task / "status", MAX_STATUS_BYTES))
    require(parsed["uid"] == uid, "same-UID task changed ownership")
    affinity = sorted(os.sched_getaffinity(tid))
    require(affinity == parsed["cpus_allowed"],
            "procfs and sched_getaffinity disagree")
    task_after = identity(task, task / "stat")
    process_after = identity(process, process / "stat")
    require(process_before == process_after and task_before == task_after,
            "proc task identity changed during scan")
    intersection = sorted(set(affinity).intersection(reserved))
    return {
        "pid": pid,
        "tid": tid,
        "uid": uid,
        "process_identity": process_before,
        "thread_identity": task_before,
        "comm_size": parsed["comm_size"],
        "comm_sha256": parsed["comm_sha256"],
        "cpus_allowed": affinity,
        "reserved_pair_intersection": intersection,
    }


def scan_same_uid(uid: int, reserved: Sequence[int]) -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    vanished: list[dict[str, Any]] = []
    numeric = sorted(
        (int(path.name) for path in Path("/proc").iterdir()
         if path.name.isdecimal()),
    )
    require(len(numeric) <= MAX_TASKS, "proc process census is unbounded")
    for pid in numeric:
        process = Path(f"/proc/{pid}")
        try:
            parsed = parse_status(bounded_read(process / "status", MAX_STATUS_BYTES))
        except FileNotFoundError:
            vanished.append(vanished_record(pid, None, "process-status"))
            continue
        if parsed["uid"] != uid:
            continue
        try:
            tids = sorted(int(path.name) for path in (process / "task").iterdir()
                          if path.name.isdecimal())
        except FileNotFoundError:
            vanished.append(vanished_record(pid, None, "task-list"))
            continue
        require(len(entries) + len(tids) <= MAX_TASKS,
                "proc thread census is unbounded")
        for tid in tids:
            try:
                entries.append(thread_entry(pid, tid, uid, reserved))
            except OSError as error:
                if error.errno in (errno.ENOENT, errno.ESRCH):
                    vanished.append(vanished_record(pid, tid, "thread-read"))
                    continue
                raise
    entries.sort(key=lambda item: (item["pid"], item["tid"]))
    require(len({(item["pid"], item["tid"]) for item in entries}) == len(entries),
            "duplicate census thread identity")
    by_process: dict[int, list[dict[str, Any]]] = {}
    for entry in entries:
        by_process.setdefault(entry["pid"], []).append(entry)
    pair_processes = sum(any(item["reserved_pair_intersection"] for item in group)
                         for group in by_process.values())
    nonuniform = sum(len({tuple(item["cpus_allowed"]) for item in group}) > 1
                     for group in by_process.values())
    confined_threads = sum(
        set(item["cpus_allowed"]).issubset(set(reserved)) for item in entries)
    confined_processes = sum(
        any(set(item["cpus_allowed"]).issubset(set(reserved)) for item in group)
        for group in by_process.values())
    return {
        "entries": entries,
        "vanished_during_scan": vanished,
        "summary": {
            "retained_process_count": len(by_process),
            "retained_thread_count": len(entries),
            "pair_eligible_process_count": pair_processes,
            "pair_eligible_thread_count": sum(
                bool(item["reserved_pair_intersection"]) for item in entries),
            "nonuniform_process_count": nonuniform,
            "confined_to_reserved_subset_process_count": confined_processes,
            "confined_to_reserved_subset_thread_count": confined_threads,
            "vanished_record_count": len(vanished),
        },
    }


def capture(phase: str, reserved: Sequence[int]) -> dict[str, Any]:
    require(phase in ("pre", "post"), "invalid census phase")
    require(len(reserved) == 2 and len(set(reserved)) == 2 and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in reserved),
            "invalid reserved CPU pair")
    reserved = sorted(reserved)
    uid = os.getuid()
    collector_affinity = sorted(os.sched_getaffinity(0))
    require(not set(collector_affinity).intersection(reserved),
            "census collector affinity includes a reserved CPU")
    clock_ticks = os.sysconf("SC_CLK_TCK")
    require(type(clock_ticks) is int and
            clock_ticks == CLOCK_TICKS_PER_SECOND,
            "passive contract requires SC_CLK_TCK=100")
    if phase == "post":
        boundary_ns = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
        cpu_stat = read_cpu_stat(reserved)
    scan_started = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    census = scan_same_uid(uid, reserved)
    scan_finished = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    if phase == "pre":
        cpu_stat = read_cpu_stat(reserved)
        boundary_ns = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    require(scan_finished >= scan_started, "census monotonic interval regressed")
    summary = census["summary"]
    require(summary["retained_process_count"] > 0 and
            summary["retained_thread_count"] > 0,
            "same-UID census retained no tasks")
    require(summary["confined_to_reserved_subset_thread_count"] == 0,
            "same-UID thread is confined to the reserved pair")
    payload = {
        "schema": SCHEMA,
        "phase": phase,
        "semantics": SEMANTICS,
        "mutation_policy": MUTATION_POLICY,
        "uid": uid,
        "boot_id": bounded_read(Path("/proc/sys/kernel/random/boot_id"), 128)
            .decode("ascii", "strict").strip(),
        "namespaces": {
            "pid": namespace_identity(Path("/proc/self/ns/pid")),
            "user": namespace_identity(Path("/proc/self/ns/user")),
        },
        "clock_ticks_per_second": clock_ticks,
        "reserved_cpus": reserved,
        "scan_started_monotonic_ns": scan_started,
        "scan_finished_monotonic_ns": scan_finished,
        "activity_boundary_monotonic_ns": boundary_ns,
        "collector": {
            "pid": os.getpid(),
            "allowed_cpus": collector_affinity,
            "reserved_pair_excluded": True,
        },
        "proc_stat": cpu_stat,
        "same_uid_processes": census,
        "capabilities": {
            "atomic_snapshot": False,
            "interval_complete": False,
            "records_execution_history": False,
            "establishes_cpu_exclusivity": False,
        },
    }
    return {**payload, "digest": digest_payload(payload)}


def exact_keys(value: Any, keys: set[str], label: str) -> Mapping[str, Any]:
    require(type(value) is dict and set(value) == keys, f"{label} keys differ")
    return value


def validate_snapshot(
    value: Any,
    phase: str,
    *,
    require_no_confined_threads: bool = True,
    cpu_pair: tuple[int, int] = LEGACY_CPU_PAIR,
) -> Mapping[str, Any]:
    canonical_pair = sorted(cpu_pair)
    pair_set = set(cpu_pair)
    value = exact_keys(value, {
        "schema", "phase", "semantics", "mutation_policy", "uid", "boot_id",
        "namespaces", "clock_ticks_per_second", "reserved_cpus",
        "scan_started_monotonic_ns", "scan_finished_monotonic_ns",
        "activity_boundary_monotonic_ns", "collector", "proc_stat",
        "same_uid_processes", "capabilities", "digest",
    }, f"{phase} census")
    payload = dict(value)
    digest = payload.pop("digest")
    require(type(digest) is str and HEX256.fullmatch(digest) is not None and
            digest == digest_payload(payload), f"{phase} census digest differs")
    capabilities = exact_keys(value["capabilities"], {
        "atomic_snapshot", "interval_complete", "records_execution_history",
        "establishes_cpu_exclusivity",
    }, f"{phase} capabilities")
    require(value["schema"] == SCHEMA and value["phase"] == phase and
            value["semantics"] == SEMANTICS and
            value["mutation_policy"] == MUTATION_POLICY and
            type(value["clock_ticks_per_second"]) is int and
            value["clock_ticks_per_second"] == CLOCK_TICKS_PER_SECOND and
            all(capabilities[name] is False for name in capabilities),
            f"{phase} census contract differs")
    reserved_cpus = value["reserved_cpus"]
    require(type(value["uid"]) is int and value["uid"] >= 0 and
            type(value["boot_id"]) is str and BOOT_ID.fullmatch(value["boot_id"])
            is not None and type(reserved_cpus) is list and
            all(type(cpu) is int for cpu in reserved_cpus) and
            reserved_cpus == canonical_pair,
            f"{phase} census host identity differs")
    namespaces = exact_keys(value["namespaces"], {"pid", "user"},
                            f"{phase} namespaces")
    for name in ("pid", "user"):
        record = exact_keys(namespaces[name], {"device", "inode"},
                            f"{phase} {name} namespace")
        require(all(type(record[field]) is int and record[field] > 0
                    for field in ("device", "inode")),
                f"{phase} namespace identity is invalid")
    started = value["scan_started_monotonic_ns"]
    finished = value["scan_finished_monotonic_ns"]
    boundary = value["activity_boundary_monotonic_ns"]
    require(all(type(item) is int and item >= 0
                for item in (started, finished, boundary)) and
            started <= finished and
            ((phase == "pre" and finished <= boundary) or
             (phase == "post" and boundary <= started)),
            f"{phase} census interval is invalid")
    collector = exact_keys(value["collector"], {
        "pid", "allowed_cpus", "reserved_pair_excluded"},
        f"{phase} collector")
    collector_cpus = collector["allowed_cpus"]
    require(type(collector["pid"]) is int and collector["pid"] > 0 and
            collector["reserved_pair_excluded"] is True and
            type(collector_cpus) is list and collector_cpus ==
                sorted(set(collector_cpus)) and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID
                for cpu in collector_cpus) and collector_cpus and
            not set(collector_cpus).intersection(value["reserved_cpus"]),
            f"{phase} collector used a reserved CPU")
    proc_stat = exact_keys(value["proc_stat"], {
        str(cpu) for cpu in canonical_pair},
                           f"{phase} proc stat")
    for cpu in canonical_pair:
        record = exact_keys(proc_stat[str(cpu)], {
            "cpu", "fields", "total_jiffies", "nonidle_jiffies"},
            f"{phase} CPU {cpu}")
        fields = exact_keys(record["fields"], set(CPU_FIELDS),
                            f"{phase} CPU {cpu} fields")
        require(type(record["cpu"]) is int and record["cpu"] == cpu and
                all(type(fields[name]) is int and fields[name] >= 0
                    for name in CPU_FIELDS) and
                type(record["total_jiffies"]) is int and
                type(record["nonidle_jiffies"]) is int and
                record["total_jiffies"] == sum(fields[name]
                                                for name in CPU_FIELDS[:8]) and
                record["nonidle_jiffies"] == sum(fields[name]
                                                  for name in NONIDLE_FIELDS),
                f"{phase} CPU counters differ")
    census = exact_keys(value["same_uid_processes"], {
        "entries", "vanished_during_scan", "summary"},
        f"{phase} same-UID census")
    entries = census["entries"]
    require(type(entries) is list and 0 < len(entries) <= MAX_TASKS,
            f"{phase} retained thread census is invalid")
    identities: list[tuple[int, int]] = []
    by_process: dict[int, list[Mapping[str, Any]]] = {}
    for index, entry_value in enumerate(entries):
        entry = exact_keys(entry_value, {
            "pid", "tid", "uid", "process_identity", "thread_identity",
            "comm_size", "comm_sha256", "cpus_allowed",
            "reserved_pair_intersection"}, f"{phase} entry {index}")
        require(all(type(entry[name]) is int and entry[name] > 0
                    for name in ("pid", "tid")) and
                type(entry["uid"]) is int and
                entry["uid"] == value["uid"] and
                type(entry["comm_size"]) is int and entry["comm_size"] >= 0 and
                type(entry["comm_sha256"]) is str and
                HEX256.fullmatch(entry["comm_sha256"]) is not None,
                f"{phase} entry identity differs")
        for identity_name in ("process_identity", "thread_identity"):
            identity_value = exact_keys(entry[identity_name], {
                "device", "inode", "starttime_ticks"},
                f"{phase} entry {identity_name}")
            require(all(type(identity_value[name]) is int and
                        identity_value[name] >= 0 for name in identity_value) and
                    identity_value["inode"] > 0,
                    f"{phase} proc identity differs")
        cpus = entry["cpus_allowed"]
        intersection = entry["reserved_pair_intersection"]
        require(type(cpus) is list and cpus == sorted(set(cpus)) and cpus and
                all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in cpus)
                and type(intersection) is list and
                all(type(cpu) is int for cpu in intersection) and
                intersection == sorted(set(cpus).intersection(pair_set)),
                f"{phase} entry affinity differs")
        identities.append((entry["pid"], entry["tid"]))
        by_process.setdefault(entry["pid"], []).append(entry)
    require(identities == sorted(identities) and
            len(set(identities)) == len(identities),
            f"{phase} census entries are not sorted and unique")
    for group in by_process.values():
        require(len({canonical_bytes(item["process_identity"])
                     for item in group}) == 1,
                f"{phase} process identity differs across threads")
    vanished = census["vanished_during_scan"]
    require(type(vanished) is list and len(vanished) <= MAX_TASKS,
            f"{phase} vanished census is invalid")
    for index, item_value in enumerate(vanished):
        item = exact_keys(item_value, {"pid", "tid", "stage"},
                          f"{phase} vanished {index}")
        require(type(item["pid"]) is int and item["pid"] > 0 and
                (item["tid"] is None or
                 (type(item["tid"]) is int and item["tid"] > 0)) and
                item["stage"] in ("process-status", "task-list", "thread-read"),
                f"{phase} vanished record differs")
    summary = exact_keys(census["summary"], {
        "retained_process_count", "retained_thread_count",
        "pair_eligible_process_count", "pair_eligible_thread_count",
        "nonuniform_process_count",
        "confined_to_reserved_subset_process_count",
        "confined_to_reserved_subset_thread_count", "vanished_record_count"},
        f"{phase} summary")
    expected_summary = {
        "retained_process_count": len(by_process),
        "retained_thread_count": len(entries),
        "pair_eligible_process_count": sum(
            any(item["reserved_pair_intersection"] for item in group)
            for group in by_process.values()),
        "pair_eligible_thread_count": sum(
            bool(item["reserved_pair_intersection"]) for item in entries),
        "nonuniform_process_count": sum(
            len({tuple(item["cpus_allowed"]) for item in group}) > 1
            for group in by_process.values()),
        "confined_to_reserved_subset_process_count": sum(
            any(set(item["cpus_allowed"]).issubset(pair_set)
                for item in group) for group in by_process.values()),
        "confined_to_reserved_subset_thread_count": sum(
            set(item["cpus_allowed"]).issubset(pair_set)
            for item in entries),
        "vanished_record_count": len(vanished),
    }
    require(all(type(item) is int and item >= 0 for item in summary.values()) and
            exact_json_equal(summary, expected_summary) and
            (not require_no_confined_threads or
             summary["confined_to_reserved_subset_thread_count"] == 0),
            f"{phase} census summary differs")
    return value


def cpu_delta(before: Mapping[str, Any], after: Mapping[str, Any]) -> dict[str, Any]:
    require(before["cpu"] == after["cpu"], "CPU identity changed")
    fields = {}
    for name in CPU_FIELDS:
        first = before["fields"][name]
        last = after["fields"][name]
        require(type(first) is int and type(last) is int and 0 <= first <= last,
                "CPU counter regressed")
        fields[name] = last - first
    return {
        "cpu": before["cpu"],
        "fields": fields,
        "total_jiffies": sum(fields[name] for name in CPU_FIELDS[:8]),
        "nonidle_jiffies": sum(fields[name] for name in NONIDLE_FIELDS),
    }


def validate_controller(
    value: Any,
    raw_schema: str,
    *,
    cpu_pair: tuple[int, int] = LEGACY_CPU_PAIR,
) -> Mapping[str, Any]:
    benchmark_cpu, reserved_sibling = cpu_pair
    contracts = {
        RAW_SCHEMA_V17: (CONTROLLER_SCHEMA_V17, "passive-v1"),
        RAW_SCHEMA_V18: (CONTROLLER_SCHEMA_V18, "passive-v2"),
        RAW_SCHEMA_V19: (CONTROLLER_SCHEMA_V19, "passive-v3"),
    }
    require(raw_schema in contracts,
            "passive controller raw schema differs")
    keys = {
        "schema", "acquisition_generation", "wrapper_pid",
        "before_allowed_cpus", "after_allowed_cpus",
        "runner_launch_allowed_cpus", "benchmark_cpu", "reserved_sibling",
        "affinity_mutation_scope", "active_affinity_supervisor_executed",
    }
    if raw_schema == RAW_SCHEMA_V19:
        keys.update({
            "verified_acquisition_sha256", "verified_bridge_sha256",
            "pair_verification_completed_monotonic_ns",
            "affinity_narrowing_started_monotonic_ns",
            "affinity_narrowing_finished_monotonic_ns",
        })
    value = exact_keys(value, keys, "passive controller affinity")
    before = value["before_allowed_cpus"]
    after = value["after_allowed_cpus"]
    launch = value["runner_launch_allowed_cpus"]
    expected_schema, expected_generation = contracts[raw_schema]
    require(value["schema"] == expected_schema and
            value["acquisition_generation"] == expected_generation and
            type(value["wrapper_pid"]) is int and value["wrapper_pid"] > 0 and
            type(before) is list and before == sorted(set(before)) and
            type(after) is list and after == sorted(set(after)) and after and
            type(launch) is list and launch == sorted(set(launch)) and
            len(before) <= MAX_CPU_LIST_ENTRIES and
            len(after) <= MAX_CPU_LIST_ENTRIES and
            len(launch) <= MAX_CPU_LIST_ENTRIES and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID
                for cpu in before + after + launch) and
            benchmark_cpu in before and reserved_sibling in before and
            after == [cpu for cpu in before if cpu not in set(cpu_pair)] and
            launch == before and
            type(value["benchmark_cpu"]) is int and
            value["benchmark_cpu"] == benchmark_cpu and
            type(value["reserved_sibling"]) is int and
            value["reserved_sibling"] == reserved_sibling and
            value["affinity_mutation_scope"] ==
                "wrapper-process-and-owned-descendants-only" and
            value["active_affinity_supervisor_executed"] is False,
            "passive controller affinity contract differs")
    if raw_schema == RAW_SCHEMA_V19:
        times = (
            value["pair_verification_completed_monotonic_ns"],
            value["affinity_narrowing_started_monotonic_ns"],
            value["affinity_narrowing_finished_monotonic_ns"],
        )
        require(all(type(value[name]) is str and
                    HEX256.fullmatch(value[name]) is not None
                    for name in ("verified_acquisition_sha256",
                                 "verified_bridge_sha256")) and
                all(type(item) is int and
                    0 <= item <= V19_MAX_MONOTONIC_NS for item in times) and
                times[0] <= times[1] <= times[2],
                "passive-v3 controller verification record differs")
    return value


def validate_v19_pair_lease(
    value: Any,
    *,
    cpu_pair: tuple[int, int],
) -> Mapping[str, Any]:
    canonical_pair = sorted(cpu_pair)
    value = exact_keys(value, {
        "device", "directory_device", "directory_inode", "inode", "lock",
        "path", "payload", "sha256",
    }, "passive-v3 CPU pair lease")
    payload = exact_keys(
        value["payload"], {"cpus", "schema", "uid"},
        "passive-v3 CPU pair lease payload")
    uid = payload["uid"]
    expected_path = (
        f"/run/user/{uid}/leopard2-cpu-leases/"
        f"leopard2-cpu-pair-{uid}-{canonical_pair[0]}-"
        f"{canonical_pair[1]}.lock"
    )
    require(payload["schema"] == PAIR_LEASE_SCHEMA and
            type(payload["cpus"]) is list and
            payload["cpus"] == canonical_pair and
            all(type(cpu) is int for cpu in payload["cpus"]) and
            type(uid) is int and uid >= 0 and
            type(value["path"]) is str and value["path"] == expected_path and
            all(type(value[name]) is int and value[name] >= 0
                for name in ("device", "directory_device", "directory_inode",
                             "inode")) and
            value["lock"] == "exclusive_nonblocking_pair_wide" and
            type(value["sha256"]) is str and
            HEX256.fullmatch(value["sha256"]) is not None and
            value["sha256"] == digest_payload(payload),
            "passive-v3 CPU pair lease identity differs")
    return value


def raw_cpu_stat(value: Any, cpu: int, label: str) -> Mapping[str, Any]:
    value = exact_keys(value, {"cpu", "fields", "total_jiffies"}, label)
    fields = exact_keys(value["fields"], set(CPU_FIELDS[:8]),
                        f"{label} fields")
    require(type(value["cpu"]) is int and value["cpu"] == cpu and
            all(type(fields[name]) is int and fields[name] >= 0
                for name in CPU_FIELDS[:8]) and
            type(value["total_jiffies"]) is int and
            value["total_jiffies"] == sum(fields.values()),
            f"{label} counters differ")
    return value


def raw_isolation_delta(before: Any, after: Any, cpu: int,
                        label: str) -> dict[str, Any]:
    before = raw_cpu_stat(before, cpu, f"{label} before")
    after = raw_cpu_stat(after, cpu, f"{label} after")
    fields: dict[str, int] = {}
    for name in CPU_FIELDS[:8]:
        require(after["fields"][name] >= before["fields"][name],
                f"{label} {name} counter regressed")
        fields[name] = after["fields"][name] - before["fields"][name]
    idle = fields["idle"] + fields["iowait"]
    nonidle = sum(fields[name] for name in NONIDLE_FIELDS)
    return {
        "cpu": cpu,
        "fields": fields,
        "idle_jiffies": idle,
        "nonidle_jiffies": nonidle,
        "total_jiffies": idle + nonidle,
    }


def window_policy() -> dict[str, Any]:
    return {
        "benchmark_cpu_max_nonidle_excess_jiffies": 0,
        "benchmark_cpu_tick_ceiling_formula":
            "floor(child_wall_duration_ns * clock_ticks_per_second / "
            "1000000000) + 1",
        "counter_source": "/proc/stat",
        "idle_fields": ["idle", "iowait"],
        "interpretation":
            "sampled 100 Hz rejection screen over the retained benchmark "
            "window only; not process ownership attribution, not an "
            "interference upper bound, and not proof of CPU exclusivity",
        "nonidle_fields": list(NONIDLE_FIELDS),
        "reserved_sibling_max_nonidle_jiffies": 0,
    }


def validate_window_endpoint(
    value: Any,
    label: str,
    *,
    before: bool,
    cpu_pair: tuple[int, int] = LEGACY_CPU_PAIR,
) -> Mapping[str, Any]:
    benchmark_cpu, reserved_sibling = cpu_pair
    value = exact_keys(value, {
        "benchmark_cpu", "monotonic_ns", "read_finished_monotonic_ns",
        "read_started_monotonic_ns", "reserved_sibling"}, label)
    started = value["read_started_monotonic_ns"]
    finished = value["read_finished_monotonic_ns"]
    boundary = value["monotonic_ns"]
    require(type(started) is int and type(finished) is int and
            type(boundary) is int and 0 <= started <= finished and
            boundary == (finished if before else started),
            f"{label} read boundaries differ")
    raw_cpu_stat(value["benchmark_cpu"], benchmark_cpu,
                 f"{label} benchmark CPU")
    raw_cpu_stat(value["reserved_sibling"], reserved_sibling,
                 f"{label} reserved sibling")
    return value


def validate_v18_window(value: Any, round_index: int, slot: int, role: str,
                        duration_ns: int, *,
                        cpu_pair: tuple[int, int] = LEGACY_CPU_PAIR
                        ) -> Mapping[str, Any]:
    benchmark_cpu, reserved_sibling = cpu_pair
    label = f"v18 invocation {round_index}/{slot} CPU window"
    value = exact_keys(value, {
        "accepted", "after", "before", "benchmark_cpu",
        "benchmark_cpu_nonidle_excess_jiffies",
        "benchmark_cpu_tick_ceiling_jiffies", "cell_id",
        "child_started_monotonic_ns", "child_wall_duration_ns",
        "clock_ticks_per_second", "delta", "implementation", "policy",
        "reserved_sibling", "reserved_sibling_nonidle_jiffies", "round",
        "schema", "slot", "window_ns"}, label)
    require(value["schema"] == CPU_WINDOW_SCHEMA and
            value["cell_id"] == "gf16-high-full" and
            type(value["round"]) is int and value["round"] == round_index and
            type(value["slot"]) is int and value["slot"] == slot and
            value["implementation"] == role and
            type(value["benchmark_cpu"]) is int and
            value["benchmark_cpu"] == benchmark_cpu and
            type(value["reserved_sibling"]) is int and
            value["reserved_sibling"] == reserved_sibling and
            type(value["clock_ticks_per_second"]) is int and
            value["clock_ticks_per_second"] == CLOCK_TICKS_PER_SECOND and
            type(value["child_started_monotonic_ns"]) is int and
            value["child_started_monotonic_ns"] >= 0 and
            type(value["child_wall_duration_ns"]) is int and
            value["child_wall_duration_ns"] == duration_ns and duration_ns > 0,
            f"{label} identity differs")
    before = validate_window_endpoint(
        value["before"], f"{label} before", before=True,
        cpu_pair=cpu_pair)
    after = validate_window_endpoint(
        value["after"], f"{label} after", before=False,
        cpu_pair=cpu_pair)
    started = value["child_started_monotonic_ns"]
    finished = started + duration_ns
    before_ns = before["monotonic_ns"]
    after_ns = after["monotonic_ns"]
    benchmark = raw_isolation_delta(
        before["benchmark_cpu"], after["benchmark_cpu"], benchmark_cpu,
        f"{label} benchmark CPU")
    sibling = raw_isolation_delta(
        before["reserved_sibling"], after["reserved_sibling"],
        reserved_sibling,
        f"{label} reserved sibling")
    ceiling = duration_ns * CLOCK_TICKS_PER_SECOND // 1_000_000_000 + 1
    excess = max(0, benchmark["nonidle_jiffies"] - ceiling)
    require(before_ns <= started < finished <= after_ns and
            before["read_finished_monotonic_ns"] <= started and
            finished <= after["read_started_monotonic_ns"] and
            type(value["window_ns"]) is int and
            value["window_ns"] == after_ns - before_ns and
            value["window_ns"] >= duration_ns and
            exact_json_equal(value["delta"], {
                "benchmark_cpu": benchmark, "reserved_sibling": sibling}) and
            exact_json_equal(value["policy"], window_policy()) and
            type(value["benchmark_cpu_tick_ceiling_jiffies"]) is int and
            value["benchmark_cpu_tick_ceiling_jiffies"] == ceiling and
            type(value["benchmark_cpu_nonidle_excess_jiffies"]) is int and
            value["benchmark_cpu_nonidle_excess_jiffies"] == excess and
            type(value["reserved_sibling_nonidle_jiffies"]) is int and
            value["reserved_sibling_nonidle_jiffies"] ==
                sibling["nonidle_jiffies"] and
            value["accepted"] is True and excess == 0 and
            sibling["nonidle_jiffies"] == 0,
            f"{label} rejection screen failed")
    return value


def validate_v18_counter_chain(
    windows: Sequence[Mapping[str, Any]],
    outer_before: Mapping[str, Any],
    outer_after: Mapping[str, Any],
    *,
    cpu_pair: tuple[int, int] = LEGACY_CPU_PAIR,
) -> None:
    """Independently nest every retained counter epoch in the outer pair."""
    benchmark_cpu, reserved_sibling = cpu_pair
    for endpoint, cpu in (("benchmark_cpu", benchmark_cpu),
                          ("reserved_sibling", reserved_sibling)):
        previous = raw_cpu_stat(
            outer_before[endpoint], cpu, f"outer {endpoint} before")
        for index, window in enumerate(windows):
            current_before = raw_cpu_stat(
                window["before"][endpoint], cpu,
                f"window {index} {endpoint} before")
            current_after = raw_cpu_stat(
                window["after"][endpoint], cpu,
                f"window {index} {endpoint} after")
            require(all(
                previous["fields"][name] <=
                current_before["fields"][name] <=
                current_after["fields"][name]
                for name in CPU_FIELDS[:8]),
                "v18 per-invocation CPU counter epochs are not ordered")
            previous = current_after
        final = raw_cpu_stat(
            outer_after[endpoint], cpu, f"outer {endpoint} after")
        require(all(previous["fields"][name] <= final["fields"][name]
                    for name in CPU_FIELDS[:8]),
                "v18 per-invocation CPU counters escape outer endpoints")


def validate_v18_isolation(
    value: Any,
    raw: Mapping[str, Any],
    *,
    cpu_pair: tuple[int, int] = LEGACY_CPU_PAIR,
) -> Mapping[str, Any]:
    benchmark_cpu, reserved_sibling = cpu_pair
    value = exact_keys(value, {
        "accepted", "after", "before", "benchmark_cpu", "delta",
        "invocation_windows", "out_of_window", "pair_lease", "policy",
        "reserved_sibling", "retained_window_count", "schema", "windowed",
        "windows_schema"}, "v18 windowed isolation")
    if raw.get("schema") == RAW_SCHEMA_V19:
        validate_v19_pair_lease(value["pair_lease"], cpu_pair=cpu_pair)
    invocations = raw.get("invocations")
    order = ("baseline", "candidate", "candidate", "baseline")
    require(type(invocations) is list and len(invocations) == 12,
            "v18 campaign invocation census differs")
    windows: list[Mapping[str, Any]] = []
    for index, invocation in enumerate(invocations):
        round_index, slot = divmod(index, len(order))
        role = order[slot]
        require(type(invocation) is dict and
                invocation.get("cell_id") == "gf16-high-full" and
                type(invocation.get("round")) is int and
                invocation["round"] == round_index and
                type(invocation.get("slot")) is int and
                invocation["slot"] == slot and
                invocation.get("implementation") == role and
                type(invocation.get("duration_ns")) is int and
                invocation["duration_ns"] > 0,
                f"v18 invocation {index} identity differs")
        windows.append(validate_v18_window(
            invocation.get("cpu_window"), round_index, slot, role,
            invocation["duration_ns"], cpu_pair=cpu_pair))
    require(value["schema"] == ISOLATION_SCHEMA_V2 and
            value["windows_schema"] == CPU_WINDOW_SCHEMA and
            type(value["benchmark_cpu"]) is int and
            value["benchmark_cpu"] == benchmark_cpu and
            type(value["reserved_sibling"]) is int and
            value["reserved_sibling"] == reserved_sibling and
            type(value["retained_window_count"]) is int and
            value["retained_window_count"] == 12 and
            exact_json_equal(value["invocation_windows"], windows) and
            all(windows[index]["after"]["read_finished_monotonic_ns"] <=
                windows[index + 1]["before"]["read_started_monotonic_ns"]
                for index in range(11)),
            "v18 isolation invocation-window closure differs")
    before = exact_keys(value["before"], {
        "benchmark_cpu", "monotonic_ns", "reserved_sibling"},
        "v18 isolation before")
    after = exact_keys(value["after"], {
        "benchmark_cpu", "monotonic_ns", "reserved_sibling"},
        "v18 isolation after")
    require(type(before["monotonic_ns"]) is int and
            type(after["monotonic_ns"]) is int and
            0 <= before["monotonic_ns"] < after["monotonic_ns"] and
            before["monotonic_ns"] <=
                windows[0]["before"]["read_started_monotonic_ns"] and
            windows[-1]["after"]["read_finished_monotonic_ns"] <=
                after["monotonic_ns"],
            "v18 isolation interval does not enclose retained windows")
    validate_v18_counter_chain(
        windows, before, after, cpu_pair=cpu_pair)
    benchmark = raw_isolation_delta(
        before["benchmark_cpu"], after["benchmark_cpu"], benchmark_cpu,
        "v18 global benchmark CPU")
    sibling = raw_isolation_delta(
        before["reserved_sibling"], after["reserved_sibling"],
        reserved_sibling,
        "v18 global reserved sibling")
    benchmark_nonidle = sum(
        window["delta"]["benchmark_cpu"]["nonidle_jiffies"]
        for window in windows)
    benchmark_ceiling = sum(
        window["benchmark_cpu_tick_ceiling_jiffies"] for window in windows)
    benchmark_excess = sum(
        window["benchmark_cpu_nonidle_excess_jiffies"] for window in windows)
    sibling_nonidle = sum(
        window["reserved_sibling_nonidle_jiffies"] for window in windows)
    child_duration = sum(window["child_wall_duration_ns"] for window in windows)
    window_duration = sum(window["window_ns"] for window in windows)
    expected_windowed = {
        "benchmark_cpu_auxiliary_class_jiffies": {
            name: sum(window["delta"]["benchmark_cpu"]["fields"][name]
                      for window in windows)
            for name in ("iowait", "irq", "nice", "softirq", "steal")},
        "benchmark_cpu_nonidle_excess_jiffies": benchmark_excess,
        "benchmark_cpu_nonidle_jiffies": benchmark_nonidle,
        "benchmark_cpu_tick_ceiling_jiffies": benchmark_ceiling,
        "child_wall_duration_total_ns": child_duration,
        "clock_ticks_per_second": CLOCK_TICKS_PER_SECOND,
        "reserved_sibling_nonidle_jiffies": sibling_nonidle,
        "window_duration_total_ns": window_duration,
    }
    expected_out_of_window = {
        "benchmark_cpu_nonidle_jiffies":
            benchmark["nonidle_jiffies"] - benchmark_nonidle,
        "duration_ns":
            after["monotonic_ns"] - before["monotonic_ns"] - window_duration,
        "gated": False,
        "interpretation":
            "runner overhead outside every retained benchmark window; "
            "disclosed so the shared-host exposure is visible, deliberately "
            "not gated because it contains no retained measurement",
        "reserved_sibling_nonidle_jiffies":
            sibling["nonidle_jiffies"] - sibling_nonidle,
    }
    require(exact_json_equal(value["delta"], {
                "benchmark_cpu": benchmark, "reserved_sibling": sibling}) and
            all(type(item) is int and item >= 0 for item in (
                expected_out_of_window["benchmark_cpu_nonidle_jiffies"],
                expected_out_of_window["duration_ns"],
                expected_out_of_window["reserved_sibling_nonidle_jiffies"])) and
            all(
                benchmark["fields"][name] >= sum(
                    window["delta"]["benchmark_cpu"]["fields"][name]
                    for window in windows) and
                sibling["fields"][name] >= sum(
                    window["delta"]["reserved_sibling"]["fields"][name]
                    for window in windows)
                for name in CPU_FIELDS[:8]) and
            exact_json_equal(value["windowed"], expected_windowed) and
            exact_json_equal(value["out_of_window"], expected_out_of_window) and
            exact_json_equal(value["policy"], {
                **window_policy(),
                "benchmark_cpu_max_nonidle_excess_jiffies": 0,
                "reserved_sibling_campaign_nonidle_gated": False,
                "windowed_gate": "per retained invocation"}) and
            benchmark["nonidle_jiffies"] > 0 and
            sibling["total_jiffies"] > 0 and
            benchmark_excess == 0 and sibling_nonidle == 0 and
            value["accepted"] is True,
            "v18 isolation aggregate or rejection policy differs")
    return value


def endpoint_identity_key(entry: Mapping[str, Any]) -> tuple[int, ...]:
    process = entry["process_identity"]
    thread = entry["thread_identity"]
    return (
        entry["pid"], entry["tid"],
        process["device"], process["inode"], process["starttime_ticks"],
        thread["device"], thread["inode"], thread["starttime_ticks"],
    )


def legacy_v17_child_wall_reference(raw: Mapping[str, Any]) -> dict[str, int]:
    """Return the descriptive wall-time reference for the CPU52 screen.

    This arithmetic does not attribute /proc/stat ticks to a process.  It only
    defines a conservative rejection reference that is independently
    recomputed from the twelve retained child wall durations.
    """
    invocations = raw.get("invocations")
    require(type(invocations) is list and
            len(invocations) == EXPECTED_INVOCATION_COUNT,
            "passive campaign invocation census differs")
    durations: list[int] = []
    for index, invocation in enumerate(invocations):
        require(type(invocation) is dict,
                f"passive invocation {index} is not an object")
        duration = invocation.get("duration_ns")
        require(type(duration) is int and duration > 0,
                f"passive invocation {index} duration differs")
        durations.append(duration)
    total = sum(durations)
    numerator = total * CLOCK_TICKS_PER_SECOND
    ceiling = (numerator + 1_000_000_000 - 1) // 1_000_000_000
    return {
        "child_wall_duration_total_ns": total,
        "child_wall_time_ceiling_jiffies": ceiling,
    }


def load_json(path: Path) -> Any:
    data = bounded_read(path, 256 << 20)
    return json.loads(data.decode("utf-8", "strict"),
                      object_pairs_hook=lambda pairs: _unique_object(pairs))


def _unique_object(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        require(key not in value, f"duplicate JSON key: {key}")
        value[key] = item
    return value


def compare(
    pre: Mapping[str, Any], post: Mapping[str, Any],
    raw: Mapping[str, Any], controller: Mapping[str, Any], *,
    expected_v19_attempt: Any = None,
) -> dict[str, Any]:
    require(type(raw) is dict and raw.get("schema") in (
                RAW_SCHEMA_V17, RAW_SCHEMA_V18, RAW_SCHEMA_V19),
            "passive campaign raw schema differs")
    raw_schema = raw["schema"]
    v19_mode = raw_schema == RAW_SCHEMA_V19
    require((expected_v19_attempt is not None) == v19_mode,
            "external v19 attempt authority is required only for passive-v3")
    qualification: dict[str, Any] | None = None
    if v19_mode:
        host_initial = raw.get("host_initial")
        host_final = raw.get("host_final")
        require(type(host_initial) is dict and type(host_final) is dict and
                v19_exact_json_equal(host_initial, host_final),
                "v19 host identity changed during the campaign")
        qualification, cpu_pair = validate_v19_qualification(
            raw.get("pair_qualification"), expected_v19_attempt,
            expected_host_identity_sha256=v19_sha256(host_initial))
    else:
        cpu_pair = cpu_pair_for_raw_schema(raw_schema)
    benchmark_cpu, reserved_sibling = cpu_pair
    pre = validate_snapshot(pre, "pre", cpu_pair=cpu_pair)
    post = validate_snapshot(post, "post", cpu_pair=cpu_pair)
    controller = validate_controller(
        controller, raw_schema, cpu_pair=cpu_pair)
    require(pre["uid"] == post["uid"] and
            pre["boot_id"] == post["boot_id"] and
            pre["namespaces"] == post["namespaces"] and
            pre["reserved_cpus"] == post["reserved_cpus"],
            "passive census identity changed")
    require("supervision" in raw and
            raw["supervision"] is None,
            "passive campaign supervision is not null")
    campaign = raw.get("campaign")
    raw_launch_cpus = (campaign.get("allowed_cpu_set_at_launch")
                       if type(campaign) is dict else None)
    require(type(campaign) is dict and type(raw_launch_cpus) is list and
            raw_launch_cpus == sorted(set(raw_launch_cpus)) and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID
                for cpu in raw_launch_cpus) and
            raw_launch_cpus == controller["runner_launch_allowed_cpus"] and
            type(campaign.get("benchmark_cpu")) is int and
            campaign["benchmark_cpu"] == benchmark_cpu and
            type(campaign.get("reserved_sibling")) is int and
            campaign["reserved_sibling"] == reserved_sibling,
            "passive runner launch affinity differs from the controller")
    if v19_mode:
        require(qualification is not None and v19_exact_json_equal(
                    raw_launch_cpus,
                    qualification["acquisition"]
                        ["allowed_cpu_set_at_launch"]) and
                v19_exact_json_equal(
                    controller["before_allowed_cpus"], raw_launch_cpus),
                "passive-v3 controller or campaign affinity differs from "
                "qualification acquisition")
        require(controller["verified_acquisition_sha256"] ==
                    qualification["acquisition_sha256"] and
                controller["verified_bridge_sha256"] ==
                    qualification["bridge_sha256"] and
                qualification["bridge"]["bridge_finished_monotonic_ns"] <=
                    controller[
                        "pair_verification_completed_monotonic_ns"] <=
                    controller[
                        "affinity_narrowing_started_monotonic_ns"] <=
                    controller[
                        "affinity_narrowing_finished_monotonic_ns"] <=
                    pre["scan_started_monotonic_ns"],
                "passive-v3 bridge verification and affinity narrowing "
                "sequence differs")
    isolation = raw.get("isolation")
    if raw_schema in (RAW_SCHEMA_V18, RAW_SCHEMA_V19):
        isolation = validate_v18_isolation(
            isolation, raw, cpu_pair=cpu_pair)
    else:
        require(type(isolation) is dict and
                set(isolation) == {
                    "accepted", "after", "before", "benchmark_cpu", "delta",
                    "pair_lease", "policy", "reserved_sibling", "schema"} and
                isolation.get("schema") == ISOLATION_SCHEMA_V1 and
                type(isolation.get("benchmark_cpu")) is int and
                isolation.get("benchmark_cpu") == benchmark_cpu and
                type(isolation.get("reserved_sibling")) is int and
                isolation.get("reserved_sibling") == reserved_sibling and
                isolation.get("accepted") is True,
                "passive v17 campaign isolation was not accepted")
    if v19_mode:
        require(isolation["pair_lease"]["payload"]["uid"] == pre["uid"],
                "passive-v3 pair lease owner differs from the census")
    isolation_before = exact_keys(isolation["before"], {
        "benchmark_cpu", "monotonic_ns", "reserved_sibling"},
        "raw isolation before")
    isolation_after = exact_keys(isolation["after"], {
        "benchmark_cpu", "monotonic_ns", "reserved_sibling"},
        "raw isolation after")
    before_ns = isolation_before["monotonic_ns"]
    after_ns = isolation_after["monotonic_ns"]
    first_retained_before: Mapping[str, Any] | None = None
    if v19_mode:
        require(qualification is not None,
                "passive-v3 isolation lacks pair qualification")
        bridge_tail = qualification["bridge"]["campaign_presample_before"]
        expected_isolation_before = {
            "benchmark_cpu": bridge_tail["cpus"][str(benchmark_cpu)],
            "monotonic_ns": bridge_tail["read_finished_monotonic_ns"],
            "reserved_sibling": bridge_tail["cpus"][str(reserved_sibling)],
        }
        first_retained_before = isolation["invocation_windows"][0]["before"]
        require(v19_exact_json_equal(
                    isolation_before, expected_isolation_before) and
                v19_exact_json_equal(
                    first_retained_before,
                    qualification["first_window_handoff"]
                        ["first_window_before"]),
                "passive-v3 bridge tail or first-window handoff differs")
        require(type(before_ns) is int and type(after_ns) is int and
                before_ns <= pre["scan_started_monotonic_ns"] <=
                    pre["scan_finished_monotonic_ns"] <=
                    pre["activity_boundary_monotonic_ns"] <=
                    first_retained_before[
                        "read_started_monotonic_ns"] <= after_ns <=
                    post["activity_boundary_monotonic_ns"] <=
                    post["scan_started_monotonic_ns"],
                "passive-v3 census is not nested in the qualified handoff "
                "and campaign interval")
    else:
        require(type(before_ns) is int and type(after_ns) is int and
                pre["scan_finished_monotonic_ns"] <=
                    pre["activity_boundary_monotonic_ns"] <=
                    before_ns <= after_ns <=
                    post["activity_boundary_monotonic_ns"] <=
                    post["scan_started_monotonic_ns"],
                "passive census does not enclose campaign isolation")
    for cpu, raw_name in ((benchmark_cpu, "benchmark_cpu"),
                          (reserved_sibling, "reserved_sibling")):
        raw_before = raw_cpu_stat(
            isolation_before[raw_name], cpu, f"raw CPU {cpu} before")
        raw_after = raw_cpu_stat(
            isolation_after[raw_name], cpu, f"raw CPU {cpu} after")
        outer_before = pre["proc_stat"][str(cpu)]["fields"]
        outer_after = post["proc_stat"][str(cpu)]["fields"]
        if v19_mode:
            require(first_retained_before is not None,
                    "passive-v3 first retained counter is absent")
            first_counter = raw_cpu_stat(
                first_retained_before[raw_name], cpu,
                f"passive-v3 first-window CPU {cpu}")
            require(all(
                raw_before["fields"][name] <= outer_before[name] <=
                    first_counter["fields"][name] <=
                    raw_after["fields"][name] <= outer_after[name]
                for name in CPU_FIELDS[:8]),
                f"CPU {cpu} counters escape the passive-v3 nested interval")
        else:
            require(all(
                outer_before[name] <= raw_before["fields"][name] <=
                    raw_after["fields"][name] <= outer_after[name]
                for name in CPU_FIELDS[:8]),
                f"raw CPU {cpu} counters escape the passive outer interval")
    cpus = pre["reserved_cpus"]
    outer = {
        str(cpu): cpu_delta(pre["proc_stat"][str(cpu)],
                            post["proc_stat"][str(cpu)]) for cpu in cpus
    }
    benchmark = outer[str(benchmark_cpu)]
    sibling = outer[str(reserved_sibling)]
    outer_contamination: dict[str, Any] | None = None
    if raw_schema == RAW_SCHEMA_V17:
        require(sibling["nonidle_jiffies"] == 0,
                "outer passive sibling activity gate failed")
        wall_reference = legacy_v17_child_wall_reference(raw)
        require(after_ns > before_ns and
                wall_reference["child_wall_duration_total_ns"] <=
                    after_ns - before_ns,
                "passive child wall durations escape the isolation interval")
        observed = benchmark["nonidle_jiffies"]
        excess = max(
            0, observed - wall_reference["child_wall_time_ceiling_jiffies"])
        require(
            excess <=
                LEGACY_V17_CPU52_AGGREGATE_EXCESS_LIMIT_JIFFIES,
            "outer benchmark CPU nonidle excess rejection threshold exceeded")
        outer_contamination = {
            "clock_ticks_per_second": CLOCK_TICKS_PER_SECOND,
            **wall_reference,
            "benchmark_cpu_nonidle_jiffies": observed,
            "benchmark_cpu_auxiliary_class_jiffies": {
                name: benchmark["fields"][name]
                for name in ("nice", "iowait", "irq", "softirq", "steal")
            },
            "benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies":
                excess,
            "policy_max_nonidle_excess_over_child_wall_ceiling_jiffies":
                LEGACY_V17_CPU52_AGGREGATE_EXCESS_LIMIT_JIFFIES,
            "reserved_sibling_nonidle_jiffies": sibling["nonidle_jiffies"],
            "auxiliary_class_policy": (
                "nice, irq, softirq, and steal are included in aggregate "
                "nonidle; iowait is disclosed as diagnostic idle accounting; "
                "no auxiliary class has a separate zero-tolerance gate"
            ),
            "interpretation": (
                "descriptive one-sided rejection screen only; the excess is "
                "not process ownership attribution or an upper bound on "
                "interference"
            ),
        }
    pre_entries = {
        endpoint_identity_key(item): item
        for item in pre["same_uid_processes"]["entries"]
    }
    post_entries = {
        endpoint_identity_key(item): item
        for item in post["same_uid_processes"]["entries"]
    }
    persistent = sorted(set(pre_entries).intersection(post_entries))
    require(all(pre_entries[key]["cpus_allowed"] ==
                post_entries[key]["cpus_allowed"] for key in persistent),
            "persistent same-UID task changed affinity between endpoints")
    wrapper_pid = controller["wrapper_pid"]
    wrapper_keys = [
        key for key in persistent if key[0] == wrapper_pid and key[1] == wrapper_pid
    ]
    require(len(wrapper_keys) == 1 and
            pre_entries[wrapper_keys[0]]["cpus_allowed"] ==
                controller["after_allowed_cpus"] and
            post_entries[wrapper_keys[0]]["cpus_allowed"] ==
                controller["after_allowed_cpus"],
            "wrapper endpoint identity/affinity is not persistent")
    wrapper = pre_entries[wrapper_keys[0]]
    if v19_mode:
        require(qualification is not None and
                first_retained_before is not None,
                "passive-v3 result lacks its conditioning chain")
        return {
            "schema": POLICY_SCHEMA_V3,
            "status": "complete",
            "acquisition_generation": "passive-v3",
            "evidence_class":
                "conditioned-passive-windowed-shared-host-observation/v1",
            "policy_evaluation_complete": True,
            "cpu_pair_exclusive": False,
            "interval_complete_task_observation": False,
            "benchmark_cpu_foreign_work_attributable": False,
            "causal_performance_claim_eligible": False,
            "promotion_eligible": False,
            "promotion_passed": False,
            "census_collector_foreign_process_affinity_mutation_performed":
                False,
            "census_collector_foreign_process_signal_sent": False,
            "campaign_supervision": None,
            "pair_qualification": {
                "schema": qualification["schema"],
                "attempt": dict(qualification["attempt"]),
                "selected_pair": dict(qualification["selected_pair"]),
                "policy_sha256": qualification["policy_sha256"],
                "acquisition_sha256": qualification["acquisition_sha256"],
                "bridge_sha256": qualification["bridge_sha256"],
                "first_window_before_sha256": qualification[
                    "first_window_handoff"]["first_window_before_sha256"],
                "v18_failure_lineage_sha256": qualification[
                    "v18_failure_lineage_sha256"],
                "shared_host_claim_ceiling": dict(
                    qualification["shared_host_claim_ceiling"]),
                "terminal": qualification["terminal"],
                "fixed_point_replayed_independently": True,
            },
            "controller_affinity": {
                "schema": controller["schema"],
                "runner_launch_allowed_cpus": list(
                    controller["runner_launch_allowed_cpus"]),
                "narrowed_wrapper_allowed_cpus": list(
                    controller["after_allowed_cpus"]),
                "launch_matches_qualification_acquisition": True,
                "verified_acquisition_sha256":
                    controller["verified_acquisition_sha256"],
                "verified_bridge_sha256":
                    controller["verified_bridge_sha256"],
                "pair_verification_completed_monotonic_ns": controller[
                    "pair_verification_completed_monotonic_ns"],
                "affinity_narrowing_started_monotonic_ns": controller[
                    "affinity_narrowing_started_monotonic_ns"],
                "affinity_narrowing_finished_monotonic_ns": controller[
                    "affinity_narrowing_finished_monotonic_ns"],
                "selected_pair_removed_only_after_verified_bridge": True,
                "active_affinity_supervisor_executed": False,
            },
            "handoff_census": {
                "schema":
                    "leopard2-v19-conditioned-pre-census-handoff/v1",
                "bridge_tail_sha256": qualification["bridge"]
                    ["bridge_tail_sha256"],
                "bridge_tail_read_finished_monotonic_ns": before_ns,
                "pre_census_digest": pre["digest"],
                "pre_scan_started_monotonic_ns":
                    pre["scan_started_monotonic_ns"],
                "pre_scan_finished_monotonic_ns":
                    pre["scan_finished_monotonic_ns"],
                "pre_activity_boundary_monotonic_ns":
                    pre["activity_boundary_monotonic_ns"],
                "first_window_before_sha256": qualification[
                    "first_window_handoff"]["first_window_before_sha256"],
                "first_window_before_read_started_monotonic_ns":
                    first_retained_before["read_started_monotonic_ns"],
                "pair_counters_nested_between_bridge_and_first_window": True,
                "accepted": True,
            },
            "windowed_contamination": {
                "schema":
                    "leopard2-v19-windowed-census-policy-observation/v1",
                "gated": True,
                "retained_window_count": isolation["retained_window_count"],
                "windowed": isolation["windowed"],
                "all_benchmark_cpu_excess_zero":
                    isolation["windowed"]
                        ["benchmark_cpu_nonidle_excess_jiffies"] == 0,
                "all_reserved_sibling_nonidle_zero":
                    isolation["windowed"]
                        ["reserved_sibling_nonidle_jiffies"] == 0,
                "policy": window_policy(),
            },
            "outer_disclosure": {
                "gated": False,
                "outer_cpu_activity": outer,
                "isolation_out_of_window": isolation["out_of_window"],
                "interpretation": (
                    "pre-handoff endpoint and post-campaign counters are "
                    "retained for shared-host disclosure only and do not "
                    "establish exclusivity or causality"
                ),
            },
            "shared_host_exposure": {
                phase: {
                    name: census["same_uid_processes"]["summary"][name]
                    for name in (
                        "retained_process_count", "retained_thread_count",
                        "pair_eligible_process_count",
                        "pair_eligible_thread_count",
                        "nonuniform_process_count",
                        "confined_to_reserved_subset_process_count",
                        "confined_to_reserved_subset_thread_count",
                        "vanished_record_count")
                }
                for phase, census in (("pre", pre), ("post", post))
            },
            "persistent_same_uid_thread_count": len(persistent),
            "wrapper_endpoint": {
                "pid": wrapper_pid,
                "process_identity": wrapper["process_identity"],
                "thread_identity": wrapper["thread_identity"],
                "cpus_allowed": wrapper["cpus_allowed"],
            },
            "shared_host_claim_ceiling": v19_claim_ceiling(),
            "interpretation": (
                "the selected pair is qualified and bridged into a nested "
                "pre-census and twelve retained rejection screens; shared-"
                "host observations remain nonexclusive, noncausal, and never "
                "promotion evidence"
            ),
        }
    if raw_schema == RAW_SCHEMA_V18:
        return {
            "schema": POLICY_SCHEMA_V2,
            "status": "complete",
            "acquisition_generation": "passive-v2",
            "evidence_class": "passive-windowed-shared-host-observation/v1",
            "policy_evaluation_complete": True,
            "cpu_pair_exclusive": False,
            "interval_complete_task_observation": False,
            "foreign_cpu52_work_attributable": False,
            "causal_performance_claim_eligible": False,
            "promotion_eligible": False,
            "promotion_passed": False,
            "census_collector_foreign_process_affinity_mutation_performed": False,
            "census_collector_foreign_process_signal_sent": False,
            "campaign_supervision": None,
            "windowed_contamination": {
                "schema":
                    "leopard2-v18-windowed-census-policy-observation/v1",
                "gated": True,
                "retained_window_count": isolation["retained_window_count"],
                "windowed": isolation["windowed"],
                "all_benchmark_cpu_excess_zero":
                    isolation["windowed"][
                        "benchmark_cpu_nonidle_excess_jiffies"] == 0,
                "all_reserved_sibling_nonidle_zero":
                    isolation["windowed"][
                        "reserved_sibling_nonidle_jiffies"] == 0,
                "policy": window_policy(),
            },
            "outer_disclosure": {
                "gated": False,
                "outer_cpu_activity": outer,
                "isolation_out_of_window": isolation["out_of_window"],
                "interpretation": (
                    "endpoint and whole-campaign CPU counters are retained "
                    "for shared-host disclosure only and do not affect v18 "
                    "acceptance"
                ),
            },
            "shared_host_exposure": {
                phase: {
                    name: census["same_uid_processes"]["summary"][name]
                    for name in (
                        "retained_process_count", "retained_thread_count",
                        "pair_eligible_process_count",
                        "pair_eligible_thread_count",
                        "nonuniform_process_count",
                        "confined_to_reserved_subset_process_count",
                        "confined_to_reserved_subset_thread_count",
                        "vanished_record_count")
                }
                for phase, census in (("pre", pre), ("post", post))
            },
            "persistent_same_uid_thread_count": len(persistent),
            "wrapper_endpoint": {
                "pid": wrapper_pid,
                "process_identity": wrapper["process_identity"],
                "thread_identity": wrapper["thread_identity"],
                "cpus_allowed": wrapper["cpus_allowed"],
            },
            "interpretation": (
                "only the twelve retained launch-through-reap windows are "
                "screened; endpoint affinity eligibility and outer counters "
                "do not establish CPU exclusivity or freedom from transient "
                "interference"
            ),
        }
    require(outer_contamination is not None,
            "v17 outer contamination policy was not evaluated")
    return {
        "schema": POLICY_SCHEMA_V1,
        "status": "complete",
        "policy_evaluation_complete": True,
        "cpu_pair_exclusive": False,
        "interval_complete_task_observation": False,
        "foreign_cpu52_work_attributable": False,
        "causal_performance_claim_eligible": False,
        "promotion_eligible": False,
        "promotion_passed": False,
        "census_collector_foreign_process_affinity_mutation_performed": False,
        "census_collector_foreign_process_signal_sent": False,
        "campaign_supervision": None,
        "outer_cpu_activity": outer,
        "outer_contamination": outer_contamination,
        "persistent_same_uid_thread_count": len(persistent),
        "wrapper_endpoint": {
            "pid": wrapper_pid,
            "process_identity": wrapper["process_identity"],
            "thread_identity": wrapper["thread_identity"],
            "cpus_allowed": wrapper["cpus_allowed"],
        },
        "interpretation": (
            "endpoint affinity eligibility and jiffy counters do not establish "
            "benchmark-CPU exclusivity or freedom from transient interference"
        ),
    }


def write_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    require(path.is_absolute() and path.parent.resolve(strict=True) == path.parent and
            not path.exists(), "output path is not a new canonical absolute file")
    data = canonical_bytes(value) + b"\n"
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                         getattr(os, "O_CLOEXEC", 0), 0o600)
    try:
        view = memoryview(data)
        while view:
            count = os.write(descriptor, view)
            require(count > 0, "short census output write")
            view = view[count:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def self_test() -> None:
    require(parse_cpu_list("0-1,4,52,116") == [0, 1, 4, 52, 116],
            "CPU-list parser failed")
    require(format_cpu_list([0, 1, 4, 52, 116]) == "0-1,4,52,116",
            "CPU-list formatter failed")
    for bad in (
            "", "01", "1-0", "1-1", "1,1", "1,,2", "1-2,2-3",
            "1-2,3", " 1", "١"):
        try:
            parse_cpu_list(bad)
        except CensusError:
            pass
        else:
            raise CensusError(f"self-test accepted CPU list: {bad!r}")
    require(len(parse_cpu_list("0-4095")) == MAX_CPU_LIST_ENTRIES,
            "maximum CPU-list size failed")
    try:
        parse_cpu_list("0-4096")
    except CensusError:
        pass
    else:
        raise CensusError("self-test accepted oversized CPU list")
    for bad in ("01,2", "1", "1,1", "1,1048576"):
        try:
            parse_oriented_cpu_pair(bad, "self-test pair")
        except CensusError:
            pass
        else:
            raise CensusError(f"self-test accepted bounded input: {bad!r}")
    require(parse_oriented_cpu_pair("1,65", "self-test pair") == {
                "benchmark_cpu": 1, "reserved_sibling": 65},
            "oriented CPU-pair parser failed")
    require(v19_exact_json_equal(
                v19_attempt_record(attempt=1), {
                    "schema": V19_ATTEMPT_SCHEMA,
                    "attempt": 1,
                    "attempt_budget": V19_ATTEMPT_BUDGET,
                    "frozen_pair_from_prior_attempt": None,
                    "prior_attempt_failure_sha256": None,
                    "prior_attempt_acquisition_sha256": None,
                    "prior_attempt_selection_status": None,
                    "prior_attempt_selected_pair": None,
                    "fresh_selection_permitted": True,
                    "selection_rule": V19_SELECTION_RULE,
                    "frozen_pair_rule": V19_FROZEN_PAIR_RULE,
                }), "v19 attempt self-test failed")
    base = {
        "cpu": 52,
        "fields": {name: 10 for name in CPU_FIELDS},
        "total_jiffies": 80,
        "nonidle_jiffies": 60,
    }
    after = json.loads(json.dumps(base))
    after["fields"]["user"] += 1
    delta = cpu_delta(base, after)
    require(delta["fields"]["user"] == 1 and delta["nonidle_jiffies"] == 1,
            "CPU delta self-test failed")
    print("passive environment census self-test passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    capture_parser = subparsers.add_parser("capture")
    capture_parser.add_argument("--phase", choices=("pre", "post"), required=True)
    capture_parser.add_argument("--reserved-cpus", default="52,116")
    capture_parser.add_argument("--output", type=Path, required=True)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--input", type=Path, required=True)
    verify_parser.add_argument("--phase", choices=("pre", "post"), required=True)
    verify_parser.add_argument(
        "--generation", choices=("passive-v1", "passive-v2", "passive-v3"),
        default="passive-v1")
    verify_parser.add_argument(
        "--reserved-cpus",
        help=("canonical dynamic pair required for passive-v3 snapshot-only "
              "verification"))
    compare_parser = subparsers.add_parser("compare")
    compare_parser.add_argument("--pre", type=Path, required=True)
    compare_parser.add_argument("--post", type=Path, required=True)
    compare_parser.add_argument("--raw", type=Path, required=True)
    compare_parser.add_argument("--controller", type=Path, required=True)
    compare_parser.add_argument("--output", type=Path, required=True)
    compare_parser.add_argument("--v19-attempt", type=int, choices=(1, 2))
    compare_parser.add_argument("--v19-prior-failure-sha256")
    compare_parser.add_argument("--v19-prior-acquisition-sha256")
    compare_parser.add_argument(
        "--v19-prior-selection-status",
        choices=("not-acquired", "no-candidate-pair-qualified",
                 "selected-lowest-primary"))
    compare_parser.add_argument(
        "--v19-prior-selected-pair", metavar="BENCHMARK,SIBLING")
    subparsers.add_parser("self-test")
    options = parser.parse_args()
    if options.command == "self-test":
        self_test()
        return 0
    if options.command == "capture":
        reserved = parse_cpu_list(options.reserved_cpus)
        value = capture(options.phase, reserved)
        write_exclusive(options.output, value)
        return 0
    if options.command == "verify":
        if options.generation == "passive-v3":
            require(options.reserved_cpus is not None,
                    "passive-v3 verify requires --reserved-cpus")
            reserved = parse_cpu_list(options.reserved_cpus)
            require(len(reserved) == 2,
                    "passive-v3 --reserved-cpus must contain exactly two CPUs")
            cpu_pair = (reserved[0], reserved[1])
        else:
            require(options.reserved_cpus is None,
                    "--reserved-cpus is only valid for passive-v3 verify")
            cpu_pair = cpu_pair_for_generation(options.generation)
        validate_snapshot(
            load_json(options.input), options.phase,
            cpu_pair=cpu_pair)
        return 0
    raw = load_json(options.raw)
    v19_options = (
        options.v19_attempt, options.v19_prior_failure_sha256,
        options.v19_prior_acquisition_sha256,
        options.v19_prior_selection_status,
        options.v19_prior_selected_pair,
    )
    expected_v19_attempt = None
    if type(raw) is dict and raw.get("schema") == RAW_SCHEMA_V19:
        require(options.v19_attempt is not None,
                "passive-v3 compare requires --v19-attempt")
        prior_pair = None
        if options.v19_prior_selected_pair is not None:
            prior_pair = parse_oriented_cpu_pair(
                options.v19_prior_selected_pair,
                "--v19-prior-selected-pair")
        expected_v19_attempt = v19_attempt_record(
            attempt=options.v19_attempt,
            prior_attempt_failure_sha256=
                options.v19_prior_failure_sha256,
            prior_attempt_acquisition_sha256=
                options.v19_prior_acquisition_sha256,
            prior_attempt_selection_status=
                options.v19_prior_selection_status,
            prior_attempt_selected_pair=prior_pair)
    else:
        require(all(value is None for value in v19_options),
                "v19 attempt flags require a passive-v3 raw record")
    result = compare(
        load_json(options.pre), load_json(options.post), raw,
        load_json(options.controller),
        expected_v19_attempt=expected_v19_attempt)
    write_exclusive(options.output, result)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (CensusError, UnicodeError, json.JSONDecodeError, OSError) as error:
        print(f"passive environment census failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
