#!/usr/bin/env python3
"""Pure, fixed-point contract for a gapless qualification-to-campaign bridge.

The bridge starts at the exact final endpoint of a validated qualification
acquisition and ends at the snapshot that a later campaign must reuse as its
presample ``before`` endpoint.  This module performs no file I/O, host discovery,
sleeping, process launch, affinity management, or timing acquisition.  All policy
and geometry values are supplied by an external preregistration.
"""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import sys
from typing import Any, Mapping, NoReturn, Sequence


def _load_contract() -> Any:
    module_name = "leopard2_pair_qualification_contract_for_bridge"
    expected = Path(__file__).resolve().with_name("pair_qualification_contract.py")
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        if Path(getattr(loaded, "__file__", "")).resolve() != expected:
            raise RuntimeError("pair qualification contract came from another path")
        return loaded
    specification = importlib.util.spec_from_file_location(module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load the pair qualification contract")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    if Path(getattr(module, "__file__", "")).resolve() != expected:
        raise RuntimeError("pair qualification contract resolved to another path")
    return module


contract = _load_contract()

ACQUISITION_SCHEMA = "leopard2-pair-qualification-acquisition/v1"
ACQUISITION_METHOD = \
    "read-only-linux-sysfs-proc-stat-shared-snapshot-chain/v1"
BRIDGE_SCHEMA = "leopard2-pair-qualification-bridge/v1"
HANDOFF_RULE = "exact-scan-tail-shared-n-plus-one-presample-before/v1"
ACCEPTANCE_RULE = "selected-pair-or-domain-zero-nonidle-live-fixed-point/v1"
MAX_NOMINAL_WINDOW_NS = 24 * 60 * 60 * 1_000_000_000

_SOURCE_ITEMS = (
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

ACQUISITION_KEYS = frozenset((
    "schema", "acquisition_method", "sources", "policy", "policy_sha256",
    "requested_window_count", "nominal_window_ns",
    "frozen_pair_from_prior_attempt", "allowed_cpu_set_at_launch",
    "allowed_cpu_set_after_scan", "clock_ticks_per_second_at_launch",
    "clock_ticks_per_second_after_scan", "topology_before_sha256",
    "topology_after_sha256", "scan", "scan_sha256",
    "host_mutation_performed", "candidate_timing_performed",
))
PAIR_KEYS = frozenset(("benchmark_cpu", "reserved_sibling"))
ACQUISITION_GEOMETRY_KEYS = frozenset((
    "requested_window_count", "nominal_window_ns",
))
BRIDGE_GEOMETRY_KEYS = frozenset((
    "minimum_window_count", "maximum_window_count", "nominal_window_ns",
    "maximum_handoff_elapsed_ns",
))
CPU_AGGREGATE_KEYS = frozenset((
    "cpu", "fields", "idle_jiffies", "nonidle_jiffies", "total_jiffies",
    "nonidle_window_count", "not_live_window_count", "overall_live",
))
CLAIM_CEILING_KEYS = frozenset((
    "promotion_eligible", "host_exclusivity_proved",
    "whole_campaign_interval_observed", "causal_performance_claim_allowed",
))
BRIDGE_KEYS = frozenset((
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
    "overall_wide_elapsed_tick_cap",
    "overall_liveness_lower_bound_jiffies", "bridge_window_count",
    "nonidle_guarded_cpus", "not_live_guarded_cpus", "bridge_accepted",
    "host_mutation_performed", "candidate_timing_performed",
    "shared_host_claim_ceiling",
))


class BridgeError(contract.QualificationError):
    """A retained bridge or one of its bindings violated the closed contract."""


def _fail(message: str) -> NoReturn:
    raise BridgeError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _bounded_int(value: Any, minimum: int, maximum: int, label: str) -> int:
    _require(type(value) is int, f"{label} is not an exact integer")
    _require(minimum <= value <= maximum, f"{label} is outside its structural bound")
    return value


def _pair(value: Any, label: str) -> dict[str, int] | None:
    if value is None:
        return None
    _require(type(value) is dict and set(value) == PAIR_KEYS,
             f"{label} is not an exact pair object")
    benchmark = _bounded_int(
        value["benchmark_cpu"], 0, contract.MAX_CPU_ID,
        f"{label} benchmark CPU")
    sibling = _bounded_int(
        value["reserved_sibling"], 0, contract.MAX_CPU_ID,
        f"{label} reserved sibling")
    _require(benchmark != sibling, f"{label} repeats one logical CPU")
    return {"benchmark_cpu": benchmark, "reserved_sibling": sibling}


def _canonical_cpu_list(value: Any, label: str) -> list[int]:
    _require(type(value) is list and 0 < len(value) <= contract.MAX_CPU_COUNT,
             f"{label} is not a bounded non-empty list")
    result = [
        _bounded_int(cpu, 0, contract.MAX_CPU_ID, f"{label} entry")
        for cpu in value
    ]
    _require(result == sorted(set(result)), f"{label} is not canonical")
    return result


def _claim_ceiling() -> dict[str, bool]:
    return {
        "promotion_eligible": False,
        "host_exclusivity_proved": False,
        "whole_campaign_interval_observed": False,
        "causal_performance_claim_allowed": False,
    }


def _acquisition_geometry(
    window_count: int, nominal_window_ns: int,
) -> dict[str, int]:
    return {
        "requested_window_count": window_count,
        "nominal_window_ns": nominal_window_ns,
    }


def bridge_geometry_record(
    *, minimum_window_count: int, maximum_window_count: int,
    nominal_window_ns: int, maximum_handoff_elapsed_ns: int,
) -> dict[str, int]:
    minimum = _bounded_int(
        minimum_window_count, 1, contract.MAX_WINDOWS,
        "minimum bridge window count")
    maximum = _bounded_int(
        maximum_window_count, minimum, contract.MAX_WINDOWS,
        "maximum bridge window count")
    nominal = _bounded_int(
        nominal_window_ns, 1, MAX_NOMINAL_WINDOW_NS,
        "nominal bridge window duration")
    elapsed = _bounded_int(
        maximum_handoff_elapsed_ns, 1, contract.MAX_MONOTONIC_NS,
        "maximum bridge handoff duration")
    _require(minimum * nominal <= elapsed,
             "bridge geometry cannot contain its minimum windows")
    return {
        "minimum_window_count": minimum,
        "maximum_window_count": maximum,
        "nominal_window_ns": nominal,
        "maximum_handoff_elapsed_ns": elapsed,
    }


def validate_bridge_geometry(value: Any) -> dict[str, int]:
    _require(type(value) is dict and set(value) == BRIDGE_GEOMETRY_KEYS,
             "bridge geometry has an unexpected key set")
    expected = bridge_geometry_record(
        minimum_window_count=value["minimum_window_count"],
        maximum_window_count=value["maximum_window_count"],
        nominal_window_ns=value["nominal_window_ns"],
        maximum_handoff_elapsed_ns=value["maximum_handoff_elapsed_ns"],
    )
    _require(contract.exact_json_equal(value, expected),
             "bridge geometry is not canonical")
    return copy.deepcopy(expected)


def validate_acquisition_for_bridge(
    value: Any, *, expected_policy: Any, expected_frozen_pair: Any,
    expected_window_count: int, expected_nominal_window_ns: int,
) -> dict[str, Any]:
    """Independently recompute the acquisition fixed point without its producer."""
    _require(type(value) is dict and set(value) == ACQUISITION_KEYS,
             "bridge acquisition has an unexpected key set")
    _require(value["schema"] == ACQUISITION_SCHEMA,
             "bridge acquisition schema differs")
    _require(value["acquisition_method"] == ACQUISITION_METHOD,
             "bridge acquisition method differs")
    sources = dict(_SOURCE_ITEMS)
    _require(type(value["sources"]) is dict and
             contract.exact_json_equal(value["sources"], sources),
             "bridge acquisition sources differ")
    policy = contract.validate_qualification_policy(expected_policy)
    _require(contract.exact_json_equal(value["policy"], policy),
             "bridge acquisition policy differs")
    _require(value["policy_sha256"] == contract.canonical_sha256(policy),
             "bridge acquisition policy hash differs")
    window_count = _bounded_int(
        expected_window_count, 1, contract.MAX_WINDOWS,
        "expected acquisition window count")
    nominal = _bounded_int(
        expected_nominal_window_ns, 1, MAX_NOMINAL_WINDOW_NS,
        "expected acquisition nominal duration")
    _require(type(value["requested_window_count"]) is int and
             value["requested_window_count"] == window_count,
             "bridge acquisition window count differs")
    _require(type(value["nominal_window_ns"]) is int and
             value["nominal_window_ns"] == nominal,
             "bridge acquisition nominal duration differs")
    frozen = _pair(expected_frozen_pair, "expected frozen pair")
    _require(contract.exact_json_equal(
        value["frozen_pair_from_prior_attempt"], frozen),
        "bridge acquisition frozen pair differs")
    allowed_launch = _canonical_cpu_list(
        value["allowed_cpu_set_at_launch"], "launch allowed CPU set")
    allowed_after = _canonical_cpu_list(
        value["allowed_cpu_set_after_scan"], "final allowed CPU set")
    _require(allowed_launch == allowed_after,
             "allowed CPU set changed across bridge acquisition")
    ticks_launch = _bounded_int(
        value["clock_ticks_per_second_at_launch"], 1,
        contract.MAX_CLOCK_TICKS_PER_SECOND, "launch clock tick rate")
    ticks_after = _bounded_int(
        value["clock_ticks_per_second_after_scan"], 1,
        contract.MAX_CLOCK_TICKS_PER_SECOND, "final clock tick rate")
    _require(ticks_launch == ticks_after == policy["clock_ticks_per_second"],
             "clock tick rate changed or differs from bridge policy")
    scan = contract.validate_pair_qualification_scan(
        value["scan"], expected_policy=policy, expected_frozen_pair=frozen)
    _require(scan["scan_window_count"] == window_count,
             "bridge acquisition scan has the wrong window count")
    _require(scan["allowed_cpu_set_at_launch"] == allowed_launch,
             "bridge acquisition scan allowed CPU set differs")
    _require(all(window["window_ns"] >= nominal for window in scan["windows"]),
             "bridge acquisition scan contains a short window")
    before_hash = contract.canonical_sha256(scan["topology_before"])
    after_hash = contract.canonical_sha256(scan["topology_after"])
    scan_hash = contract.canonical_sha256(scan)
    _require(value["topology_before_sha256"] == before_hash and
             value["topology_after_sha256"] == after_hash,
             "bridge acquisition topology hash differs")
    _require(value["scan_sha256"] == scan_hash,
             "bridge acquisition scan hash differs")
    _require(type(value["host_mutation_performed"]) is bool and
             value["host_mutation_performed"] is False,
             "bridge acquisition reports host mutation")
    _require(type(value["candidate_timing_performed"]) is bool and
             value["candidate_timing_performed"] is False and
             scan["candidate_timing_performed"] is False,
             "bridge acquisition reports candidate timing")
    expected = {
        "schema": ACQUISITION_SCHEMA,
        "acquisition_method": ACQUISITION_METHOD,
        "sources": sources,
        "policy": copy.deepcopy(policy),
        "policy_sha256": contract.canonical_sha256(policy),
        "requested_window_count": window_count,
        "nominal_window_ns": nominal,
        "frozen_pair_from_prior_attempt": copy.deepcopy(frozen),
        "allowed_cpu_set_at_launch": allowed_launch,
        "allowed_cpu_set_after_scan": allowed_after,
        "clock_ticks_per_second_at_launch": ticks_launch,
        "clock_ticks_per_second_after_scan": ticks_after,
        "topology_before_sha256": before_hash,
        "topology_after_sha256": after_hash,
        "scan": copy.deepcopy(scan),
        "scan_sha256": scan_hash,
        "host_mutation_performed": False,
        "candidate_timing_performed": False,
    }
    _require(contract.exact_json_equal(value, expected),
             "bridge acquisition is not canonical")
    return copy.deepcopy(expected)


def _tick_cap(elapsed_ns: int, ticks_per_second: int) -> int:
    return (elapsed_ns * ticks_per_second) // 1_000_000_000 + 1


def _counter_delta(
    before: Mapping[str, Any], after: Mapping[str, Any], cpu: int,
) -> dict[str, Any]:
    fields: dict[str, int] = {}
    for field in contract.COUNTER_FIELDS:
        earlier = before["fields"][field]
        later = after["fields"][field]
        _require(later >= earlier, f"bridge CPU {cpu} counter regressed")
        fields[field] = later - earlier
    idle = sum(fields[field] for field in contract.IDLE_FIELDS)
    nonidle = sum(fields[field] for field in contract.NONIDLE_FIELDS)
    total = _bounded_int(
        sum(fields.values()), 0, contract.MAX_COUNTER,
        f"bridge CPU {cpu} aggregate total")
    return {
        "cpu": cpu,
        "fields": fields,
        "idle_jiffies": idle,
        "nonidle_jiffies": nonidle,
        "total_jiffies": total,
    }


def _bridge_from_validated_acquisition(
    acquisition: dict[str, Any], *, policy: dict[str, Any],
    frozen_pair: dict[str, int] | None, acquisition_window_count: int,
    acquisition_nominal_window_ns: int, geometry: dict[str, int],
    snapshots: Sequence[Any],
) -> dict[str, Any]:
    scan = acquisition["scan"]
    selected = contract.selected_pair_from_scan(
        scan, expected_policy=policy, expected_frozen_pair=frozen_pair)
    scan_tail = scan["windows"][-1]["after"]
    observed_cpus = sorted(int(cpu) for cpu in scan_tail["cpus"])
    _require(type(snapshots) in (list, tuple),
             "bridge snapshots are not a sequence")
    window_count = len(snapshots) - 1
    _require(geometry["minimum_window_count"] <= window_count <=
             geometry["maximum_window_count"],
             "bridge snapshot chain has the wrong bounded length")
    _require(contract.exact_json_equal(snapshots[0], scan_tail),
             "bridge head differs from the qualification scan tail")
    windows = [
        contract.pair_observation_window_record(
            policy, phase="bridge", index=index,
            before=snapshots[index], after=snapshots[index + 1],
            observed_cpus=observed_cpus)
        for index in range(window_count)
    ]
    _require(all(
        window["window_ns"] >= geometry["nominal_window_ns"]
        for window in windows), "bridge contains a short observation window")
    first = windows[0]["before"]
    last = windows[-1]["after"]
    started = first["read_finished_monotonic_ns"]
    finished = last["read_finished_monotonic_ns"]
    _require(started <= contract.MAX_MONOTONIC_NS -
             geometry["maximum_handoff_elapsed_ns"],
             "bridge deadline overflows the monotonic clock bound")
    deadline = started + geometry["maximum_handoff_elapsed_ns"]
    _require(finished <= deadline, "bridge closes after its fixed deadline")
    narrow_ns = last["read_started_monotonic_ns"] - started
    wide_ns = (
        last["read_finished_monotonic_ns"] -
        first["read_started_monotonic_ns"]
    )
    _require(narrow_ns > 0 and wide_ns >= narrow_ns,
             "bridge aggregate interval is invalid")
    ticks = policy["clock_ticks_per_second"]
    overall_cap = _tick_cap(narrow_ns, ticks)
    overall_wide_cap = _tick_cap(wide_ns, ticks)
    overall_lower = overall_cap - contract.LIVENESS_LOWER_SLACK_JIFFIES
    nonidle_sets = [set(window["nonidle_cpus"]) for window in windows]
    not_live_sets = [set(window["not_live_cpus"]) for window in windows]
    aggregates: dict[str, Any] = {}
    for cpu in observed_cpus:
        aggregate = _counter_delta(
            first["cpus"][str(cpu)], last["cpus"][str(cpu)], cpu)
        summed_fields = {
            field: sum(
                window["delta"][str(cpu)]["fields"][field]
                for window in windows)
            for field in contract.COUNTER_FIELDS
        }
        _require(summed_fields == aggregate["fields"],
                 "bridge window chain does not telescope")
        aggregates[str(cpu)] = {
            **aggregate,
            "nonidle_window_count": sum(
                cpu in nonidle for nonidle in nonidle_sets),
            "not_live_window_count": sum(
                cpu in not_live for not_live in not_live_sets),
            "overall_live": (
                overall_lower >= 1 and
                overall_lower <= aggregate["total_jiffies"] <= overall_wide_cap
            ),
        }
    candidates = [
        candidate for candidate in scan["candidate_pairs"]
        if candidate["benchmark_cpu"] == selected["benchmark_cpu"] and
        candidate["reserved_sibling"] == selected["reserved_sibling"]
    ]
    _require(len(candidates) == 1 and candidates[0]["qualified"] is True,
             "selected bridge pair is not one qualified candidate")
    guarded_cpus = (
        list(candidates[0]["domain_cpus"])
        if policy["domain_mode"] == "pair-and-domain"
        else sorted(selected.values())
    )
    _require(set(guarded_cpus) <= set(observed_cpus),
             "bridge guard CPUs escape the observed CPU universe")
    nonidle_guarded = [
        cpu for cpu in guarded_cpus
        if aggregates[str(cpu)]["nonidle_jiffies"] != 0 or
        aggregates[str(cpu)]["nonidle_window_count"] != 0
    ]
    not_live_guarded = [
        cpu for cpu in guarded_cpus
        if not aggregates[str(cpu)]["overall_live"] or
        aggregates[str(cpu)]["not_live_window_count"] != 0
    ]
    acquisition_geometry = _acquisition_geometry(
        acquisition_window_count, acquisition_nominal_window_ns)
    scan_tail_hash = contract.canonical_sha256(scan_tail)
    bridge_tail_hash = contract.canonical_sha256(last)
    return {
        "schema": BRIDGE_SCHEMA,
        "handoff_rule": HANDOFF_RULE,
        "acceptance_rule": ACCEPTANCE_RULE,
        "acquisition_sha256": contract.canonical_sha256(acquisition),
        "scan_sha256": contract.canonical_sha256(scan),
        "policy_sha256": contract.canonical_sha256(policy),
        "selected_pair": copy.deepcopy(selected),
        "frozen_pair_from_prior_attempt": copy.deepcopy(frozen_pair),
        "acquisition_geometry": acquisition_geometry,
        "acquisition_geometry_sha256": contract.canonical_sha256(
            acquisition_geometry),
        "bridge_geometry": copy.deepcopy(geometry),
        "bridge_geometry_sha256": contract.canonical_sha256(geometry),
        "observed_cpus": observed_cpus,
        "guarded_cpus": guarded_cpus,
        "scan_tail_sha256": scan_tail_hash,
        "windows": windows,
        "cpu_aggregates": aggregates,
        "bridge_head_sha256": contract.canonical_sha256(first),
        "bridge_tail_sha256": bridge_tail_hash,
        "campaign_presample_before": copy.deepcopy(last),
        "campaign_presample_before_sha256": bridge_tail_hash,
        "bridge_started_monotonic_ns": started,
        "bridge_finished_monotonic_ns": finished,
        "bridge_deadline_monotonic_ns": deadline,
        "realized_bridge_total_ns": narrow_ns,
        "wide_bridge_total_ns": wide_ns,
        "overall_elapsed_tick_cap": overall_cap,
        "overall_wide_elapsed_tick_cap": overall_wide_cap,
        "overall_liveness_lower_bound_jiffies": overall_lower,
        "bridge_window_count": window_count,
        "nonidle_guarded_cpus": nonidle_guarded,
        "not_live_guarded_cpus": not_live_guarded,
        "bridge_accepted": not nonidle_guarded and not not_live_guarded,
        "host_mutation_performed": False,
        "candidate_timing_performed": False,
        "shared_host_claim_ceiling": _claim_ceiling(),
    }


def pair_qualification_bridge_record(
    acquisition_value: Any, *, expected_policy: Any,
    expected_frozen_pair: Any, expected_acquisition_window_count: int,
    expected_acquisition_nominal_window_ns: int,
    expected_bridge_geometry: Any, snapshots: Sequence[Any],
) -> dict[str, Any]:
    policy = contract.validate_qualification_policy(expected_policy)
    frozen = _pair(expected_frozen_pair, "expected frozen pair")
    acquisition = validate_acquisition_for_bridge(
        acquisition_value, expected_policy=policy,
        expected_frozen_pair=frozen,
        expected_window_count=expected_acquisition_window_count,
        expected_nominal_window_ns=expected_acquisition_nominal_window_ns)
    geometry = validate_bridge_geometry(expected_bridge_geometry)
    return _bridge_from_validated_acquisition(
        acquisition, policy=policy, frozen_pair=frozen,
        acquisition_window_count=expected_acquisition_window_count,
        acquisition_nominal_window_ns=expected_acquisition_nominal_window_ns,
        geometry=geometry, snapshots=snapshots)


def validate_pair_qualification_bridge(
    value: Any, acquisition_value: Any, *, expected_policy: Any,
    expected_frozen_pair: Any, expected_acquisition_window_count: int,
    expected_acquisition_nominal_window_ns: int,
    expected_bridge_geometry: Any,
) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == BRIDGE_KEYS,
             "pair qualification bridge has an unexpected key set")
    _require(value["schema"] == BRIDGE_SCHEMA,
             "pair qualification bridge schema differs")
    _require(value["handoff_rule"] == HANDOFF_RULE and
             value["acceptance_rule"] == ACCEPTANCE_RULE,
             "pair qualification bridge rule differs")
    policy = contract.validate_qualification_policy(expected_policy)
    frozen = _pair(expected_frozen_pair, "expected frozen pair")
    acquisition = validate_acquisition_for_bridge(
        acquisition_value, expected_policy=policy,
        expected_frozen_pair=frozen,
        expected_window_count=expected_acquisition_window_count,
        expected_nominal_window_ns=expected_acquisition_nominal_window_ns)
    geometry = validate_bridge_geometry(expected_bridge_geometry)
    scan_tail = acquisition["scan"]["windows"][-1]["after"]
    observed_cpus = sorted(int(cpu) for cpu in scan_tail["cpus"])
    _require(type(value["windows"]) is list and
             geometry["minimum_window_count"] <= len(value["windows"]) <=
             geometry["maximum_window_count"],
             "pair qualification bridge windows have the wrong length")
    windows: list[dict[str, Any]] = []
    endpoints: list[dict[str, Any]] = []
    for index, raw_window in enumerate(value["windows"]):
        window = contract.validate_pair_observation_window(
            raw_window, expected_policy=policy, expected_phase="bridge",
            expected_index=index, expected_cpus=observed_cpus)
        if index == 0:
            _require(contract.exact_json_equal(window["before"], scan_tail),
                     "bridge head differs from the qualification scan tail")
            endpoints.append(window["before"])
        else:
            _require(contract.exact_json_equal(
                windows[-1]["after"], window["before"]),
                "bridge windows do not share an exact N+1 endpoint")
        endpoints.append(window["after"])
        windows.append(window)
    expected = _bridge_from_validated_acquisition(
        acquisition, policy=policy, frozen_pair=frozen,
        acquisition_window_count=expected_acquisition_window_count,
        acquisition_nominal_window_ns=expected_acquisition_nominal_window_ns,
        geometry=geometry, snapshots=endpoints)
    _require(contract.exact_json_equal(value, expected),
             "pair qualification bridge differs from fixed-point recomputation")
    return copy.deepcopy(expected)


def load_pair_qualification_bridge(
    data: bytes, acquisition_value: Any, *, expected_policy: Any,
    expected_frozen_pair: Any, expected_acquisition_window_count: int,
    expected_acquisition_nominal_window_ns: int,
    expected_bridge_geometry: Any,
) -> dict[str, Any]:
    return validate_pair_qualification_bridge(
        contract.strict_json_loads(data, "pair qualification bridge"),
        acquisition_value, expected_policy=expected_policy,
        expected_frozen_pair=expected_frozen_pair,
        expected_acquisition_window_count=expected_acquisition_window_count,
        expected_acquisition_nominal_window_ns=
        expected_acquisition_nominal_window_ns,
        expected_bridge_geometry=expected_bridge_geometry)


__all__ = (
    "ACCEPTANCE_RULE", "BRIDGE_SCHEMA", "BridgeError", "HANDOFF_RULE",
    "bridge_geometry_record", "load_pair_qualification_bridge",
    "pair_qualification_bridge_record", "validate_acquisition_for_bridge",
    "validate_bridge_geometry", "validate_pair_qualification_bridge",
)
