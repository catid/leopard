#!/usr/bin/env python3
"""Pure, policy-parameterized CPU-pair qualification records.

This module deliberately performs no host discovery, file I/O, sleeping, affinity
management, or timing.  It consumes already-retained topology and counter snapshots,
constructs closed JSON records, and independently recomputes those records during
validation.  Acquisition and independent-auditor wiring belong to a later versioned
milestone.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from typing import Any, NoReturn


POLICY_SCHEMA = "leopard2-pair-qualification-policy/v1"
TOPOLOGY_SCHEMA = "leopard2-pair-qualification-topology/v1"
WINDOW_SCHEMA = "leopard2-pair-qualification-window/v1"
SCAN_SCHEMA = "leopard2-pair-qualification-scan/v1"

COUNTER_SOURCE = "/proc/stat"
COUNTER_FIELDS = (
    "user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal",
)
IDLE_FIELDS = ("idle", "iowait")
NONIDLE_FIELDS = ("user", "nice", "system", "irq", "softirq", "steal")
TICK_CAP_FORMULA = "(elapsed_ns * clock_ticks_per_second) // 1000000000 + 1"
LIVENESS_LOWER_SLACK_JIFFIES = 2
LIVENESS_RULE = (
    f"narrow_tick_cap - {LIVENESS_LOWER_SLACK_JIFFIES} >= 1 and "
    f"narrow_tick_cap - {LIVENESS_LOWER_SLACK_JIFFIES} <= "
    "total_jiffies_delta <= wide_tick_cap"
)
CANDIDATE_RULE = "explicit-primary-symmetric-smt2-online-allowed/v1"
SELECTION_RULE = "lowest-primary-unless-frozen/v1"
FROZEN_PAIR_RULE = "require-frozen-when-present-no-fallback/v1"

DOMAIN_MODES = frozenset(("pair-only", "pair-and-domain"))
WINDOW_PHASES = frozenset(("qualification", "bridge"))
NONIDLE_LIMIT_JIFFIES = 0

# These are parser/resource bounds, not evidence thresholds.
MAX_JSON_BYTES = 64 << 20
MAX_CPU_COUNT = 4096
MAX_CPU_ID = (1 << 20) - 1
MAX_WINDOWS = 4096
MAX_COUNTER = (1 << 64) - 1
MAX_MONOTONIC_NS = (1 << 63) - 1
MAX_CLOCK_TICKS_PER_SECOND = 1_000_000
MAX_REASON_LENGTH = 512
DISPLAY_EXCLUDED_CODE_POINTS = frozenset((
    0x200B, 0x200C, 0x200D, 0x200E, 0x200F,
    0x2028, 0x2029, 0x202A, 0x202B, 0x202C, 0x202D, 0x202E,
    0x2066, 0x2067, 0x2068, 0x2069, 0xFEFF,
))

POLICY_KEYS = frozenset((
    "schema", "clock_ticks_per_second", "candidate_primary_cpus",
    "excluded_pairs", "domain_mode", "candidate_rule", "selection_rule",
    "frozen_pair_rule", "observation",
))
OBSERVATION_POLICY_KEYS = frozenset((
    "counter_source", "counter_fields", "idle_fields", "nonidle_fields",
    "tick_cap_formula", "liveness_lower_slack_jiffies",
    "nonidle_limit_jiffies", "liveness_rule",
))
PAIR_KEYS = frozenset(("benchmark_cpu", "reserved_sibling"))
EXCLUSION_KEYS = frozenset(("benchmark_cpu", "reserved_sibling", "reason"))
TOPOLOGY_KEYS = frozenset(("schema", "cpus"))
TOPOLOGY_CPU_KEYS = frozenset((
    "cpu", "online", "physical_package_id", "die_id", "core_id",
    "thread_siblings", "domain_cpus",
))
SNAPSHOT_KEYS = frozenset((
    "cpus", "read_started_monotonic_ns", "read_finished_monotonic_ns",
))
CPU_COUNTER_KEYS = frozenset(("cpu", "fields", "total_jiffies"))
CPU_DELTA_KEYS = frozenset((
    "cpu", "fields", "idle_jiffies", "nonidle_jiffies", "total_jiffies",
))
WINDOW_POLICY_KEYS = OBSERVATION_POLICY_KEYS
WINDOW_KEYS = frozenset((
    "schema", "phase", "index", "before", "after", "window_ns",
    "wide_window_ns", "elapsed_tick_cap", "wide_elapsed_tick_cap",
    "liveness_lower_bound_jiffies", "delta", "nonidle_cpus",
    "not_live_cpus", "clock_ticks_per_second", "policy",
))
CPU_AGGREGATE_KEYS = frozenset((
    "cpu", "fields", "idle_jiffies", "nonidle_jiffies", "total_jiffies",
    "nonidle_window_count", "not_live_window_count", "overall_live",
))
CANDIDATE_RESULT_KEYS = frozenset((
    "benchmark_cpu", "reserved_sibling", "domain_cpus", "pair_idle",
    "pair_live", "domain_idle", "domain_live", "qualified",
))
SCAN_KEYS = frozenset((
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
))


class QualificationError(ValueError):
    """A retained qualification object violates the pure contract."""


def _fail(message: str) -> NoReturn:
    raise QualificationError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _require_dict(value: Any, expected_keys: frozenset[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label} is not an exact JSON object")
    _require(set(value) == expected_keys, f"{label} has an unexpected key set")
    return value


def _bounded_int(value: Any, minimum: int, maximum: int, label: str) -> int:
    _require(type(value) is int, f"{label} is not an exact integer")
    _require(minimum <= value <= maximum, f"{label} is outside its structural bound")
    return value


def _cpu_id(value: Any, label: str) -> int:
    return _bounded_int(value, 0, MAX_CPU_ID, label)


def _bounded_text(value: Any, label: str) -> str:
    _require(type(value) is str and 0 < len(value) <= MAX_REASON_LENGTH,
             f"{label} is not a bounded non-empty string")
    _require(all(
        ord(character) >= 0x20 and
        not (0x7F <= ord(character) <= 0x9F) and
        not (0xD800 <= ord(character) <= 0xDFFF) and
        ord(character) not in DISPLAY_EXCLUDED_CODE_POINTS
        for character in value
    ), f"{label} contains a non-display-safe character")
    return value


def _canonical_int_list(
    values: Any,
    *,
    label: str,
    maximum_length: int = MAX_CPU_COUNT,
    allow_empty: bool = True,
) -> list[int]:
    _require(type(values) in (list, tuple), f"{label} is not an integer sequence")
    _require((allow_empty or values) and len(values) <= maximum_length,
             f"{label} has an invalid length")
    result = [_cpu_id(value, f"{label} entry") for value in values]
    _require(len(set(result)) == len(result), f"{label} contains a duplicate")
    return sorted(result)


def exact_json_equal(left: Any, right: Any) -> bool:
    """Compare JSON values without Python's bool/int or int/float coercions."""
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return set(left) == set(right) and all(
            exact_json_equal(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            exact_json_equal(a, b) for a, b in zip(left, right)
        )
    return left == right


def strict_json_loads(data: bytes, label: str = "qualification JSON") -> Any:
    """Load exactly one finite JSON value with duplicate-key rejection."""
    _require(type(data) is bytes, f"{label} is not a byte string")
    _require(0 < len(data) <= MAX_JSON_BYTES, f"{label} has an invalid byte length")

    def object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise QualificationError(f"{label} contains duplicate object key {key!r}")
            result[key] = value
        return result

    def reject_constant(token: str) -> NoReturn:
        raise QualificationError(f"{label} contains non-finite number {token}")

    def finite_float(token: str) -> float:
        value = float(token)
        if not math.isfinite(value):
            raise QualificationError(f"{label} contains a non-finite float")
        return value

    try:
        text = data.decode("utf-8")
        return json.loads(
            text,
            object_pairs_hook=object_pairs,
            parse_constant=reject_constant,
            parse_float=finite_float,
        )
    except QualificationError:
        raise
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError, RecursionError) as error:
        raise QualificationError(f"{label} is not strict singleton JSON") from error


def canonical_json_bytes(value: Any) -> bytes:
    try:
        return (json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        ) + "\n").encode("utf-8")
    except (TypeError, ValueError, UnicodeError, RecursionError) as error:
        raise QualificationError("value is not canonical finite JSON") from error


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _observation_policy() -> dict[str, Any]:
    return {
        "counter_source": COUNTER_SOURCE,
        "counter_fields": list(COUNTER_FIELDS),
        "idle_fields": list(IDLE_FIELDS),
        "nonidle_fields": list(NONIDLE_FIELDS),
        "tick_cap_formula": TICK_CAP_FORMULA,
        "liveness_lower_slack_jiffies": LIVENESS_LOWER_SLACK_JIFFIES,
        "nonidle_limit_jiffies": NONIDLE_LIMIT_JIFFIES,
        "liveness_rule": LIVENESS_RULE,
    }


def _pair_record(benchmark_cpu: Any, reserved_sibling: Any, label: str) -> dict[str, int]:
    benchmark = _cpu_id(benchmark_cpu, f"{label} benchmark CPU")
    sibling = _cpu_id(reserved_sibling, f"{label} reserved sibling")
    _require(benchmark != sibling, f"{label} names the same logical CPU twice")
    return {"benchmark_cpu": benchmark, "reserved_sibling": sibling}


def _validate_pair(value: Any, label: str) -> dict[str, int]:
    pair = _require_dict(value, PAIR_KEYS, label)
    expected = _pair_record(pair["benchmark_cpu"], pair["reserved_sibling"], label)
    _require(exact_json_equal(pair, expected), f"{label} is not canonical")
    return expected


def qualification_policy_record(
    *,
    clock_ticks_per_second: int,
    candidate_primary_cpus: Sequence[int],
    excluded_pairs: Sequence[Mapping[str, Any]],
    domain_mode: str,
) -> dict[str, Any]:
    ticks = _bounded_int(
        clock_ticks_per_second, 1, MAX_CLOCK_TICKS_PER_SECOND,
        "qualification clock tick rate",
    )
    primaries = _canonical_int_list(
        candidate_primary_cpus, label="candidate primary CPUs",
    )
    primary_set = set(primaries)
    _require(type(excluded_pairs) in (list, tuple),
             "excluded pairs are not a sequence")
    _require(len(excluded_pairs) <= len(primaries), "too many excluded pairs")
    exclusions: list[dict[str, Any]] = []
    seen_physical_pairs: set[frozenset[int]] = set()
    seen_excluded_primaries: set[int] = set()
    for index, value in enumerate(excluded_pairs):
        exclusion = _require_dict(value, EXCLUSION_KEYS, f"excluded pair {index}")
        pair = _pair_record(
            exclusion["benchmark_cpu"], exclusion["reserved_sibling"],
            f"excluded pair {index}",
        )
        _require(pair["benchmark_cpu"] in primary_set,
                 "excluded pair primary is outside the candidate policy")
        _require(pair["benchmark_cpu"] not in seen_excluded_primaries,
                 "duplicate excluded candidate primary")
        seen_excluded_primaries.add(pair["benchmark_cpu"])
        reason = _bounded_text(exclusion["reason"], f"excluded pair {index} reason")
        physical = frozenset(pair.values())
        _require(physical not in seen_physical_pairs, "duplicate excluded physical pair")
        seen_physical_pairs.add(physical)
        exclusions.append({**pair, "reason": reason})
    exclusions.sort(key=lambda item: (item["benchmark_cpu"], item["reserved_sibling"]))
    _require(type(domain_mode) is str and domain_mode in DOMAIN_MODES,
             "domain mode is not recognized")
    return {
        "schema": POLICY_SCHEMA,
        "clock_ticks_per_second": ticks,
        "candidate_primary_cpus": primaries,
        "excluded_pairs": exclusions,
        "domain_mode": domain_mode,
        "candidate_rule": CANDIDATE_RULE,
        "selection_rule": SELECTION_RULE,
        "frozen_pair_rule": FROZEN_PAIR_RULE,
        "observation": _observation_policy(),
    }


def validate_qualification_policy(value: Any) -> dict[str, Any]:
    policy = _require_dict(value, POLICY_KEYS, "qualification policy")
    _require(policy["schema"] == POLICY_SCHEMA, "qualification policy schema differs")
    _require_dict(policy["observation"], OBSERVATION_POLICY_KEYS,
                  "qualification observation policy")
    expected = qualification_policy_record(
        clock_ticks_per_second=policy["clock_ticks_per_second"],
        candidate_primary_cpus=policy["candidate_primary_cpus"],
        excluded_pairs=policy["excluded_pairs"],
        domain_mode=policy["domain_mode"],
    )
    _require(exact_json_equal(policy, expected), "qualification policy is not canonical")
    return copy.deepcopy(expected)


def _topology_cpu_record(value: Any, label: str) -> dict[str, Any]:
    entry = _require_dict(value, TOPOLOGY_CPU_KEYS, label)
    cpu = _cpu_id(entry["cpu"], f"{label} CPU")
    _require(type(entry["online"]) is bool, f"{label} online flag is not boolean")
    package = _bounded_int(
        entry["physical_package_id"], 0, (1 << 31) - 1, f"{label} package ID",
    )
    die = _bounded_int(entry["die_id"], 0, (1 << 31) - 1, f"{label} die ID")
    core = _bounded_int(entry["core_id"], 0, (1 << 31) - 1, f"{label} core ID")
    siblings = _canonical_int_list(
        entry["thread_siblings"], label=f"{label} thread siblings", allow_empty=False,
    )
    _require(len(siblings) in (1, 2), f"{label} is not an SMT1/SMT2 topology")
    _require(cpu in siblings, f"{label} does not contain itself in thread siblings")
    domain = _canonical_int_list(
        entry["domain_cpus"], label=f"{label} domain CPUs", allow_empty=False,
    )
    _require(cpu in domain, f"{label} does not contain itself in its domain")
    return {
        "cpu": cpu,
        "online": entry["online"],
        "physical_package_id": package,
        "die_id": die,
        "core_id": core,
        "thread_siblings": siblings,
        "domain_cpus": domain,
    }


def topology_record(cpus: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    _require(type(cpus) in (list, tuple), "topology CPUs are not a sequence")
    _require(0 < len(cpus) <= MAX_CPU_COUNT, "topology CPU count is invalid")
    entries = [
        _topology_cpu_record(value, f"topology CPU {index}")
        for index, value in enumerate(cpus)
    ]
    entries.sort(key=lambda item: item["cpu"])
    _require(len({item["cpu"] for item in entries}) == len(entries),
             "topology contains a duplicate logical CPU")
    by_cpu = {item["cpu"]: item for item in entries}
    domain_signature_ids: dict[tuple[int, ...], int] = {}
    domain_id_by_cpu: dict[int, int] = {}
    for entry in entries:
        signature = tuple(entry["domain_cpus"])
        domain_id_by_cpu[entry["cpu"]] = domain_signature_ids.setdefault(
            signature, len(domain_signature_ids),
        )
    threads_by_physical_core: dict[tuple[int, int, int], list[int]] = {}
    for entry in entries:
        physical_core = (
            entry["physical_package_id"], entry["die_id"], entry["core_id"],
        )
        prior_threads = threads_by_physical_core.setdefault(
            physical_core, entry["thread_siblings"],
        )
        _require(prior_threads == entry["thread_siblings"],
                 "one physical core identity names multiple sibling groups")
        for sibling_cpu in entry["thread_siblings"]:
            _require(sibling_cpu in by_cpu, "topology names a missing thread sibling")
            sibling = by_cpu[sibling_cpu]
            _require(sibling["thread_siblings"] == entry["thread_siblings"],
                     "topology thread siblings are asymmetric")
            _require(
                sibling["physical_package_id"] == entry["physical_package_id"] and
                sibling["die_id"] == entry["die_id"] and
                sibling["core_id"] == entry["core_id"],
                "topology thread siblings disagree on physical core identity",
            )
            _require(domain_id_by_cpu[sibling_cpu] == domain_id_by_cpu[entry["cpu"]],
                     "topology thread siblings disagree on domain membership")
        for domain_cpu in entry["domain_cpus"]:
            _require(domain_cpu in by_cpu, "topology names a missing domain CPU")
            domain_member = by_cpu[domain_cpu]
            _require(domain_id_by_cpu[domain_cpu] == domain_id_by_cpu[entry["cpu"]],
                     "topology domain membership is asymmetric")
            _require(
                domain_member["physical_package_id"] == entry["physical_package_id"] and
                domain_member["die_id"] == entry["die_id"],
                "topology domain crosses a package or die",
            )
    return {"schema": TOPOLOGY_SCHEMA, "cpus": entries}


def validate_topology(value: Any) -> dict[str, Any]:
    topology = _require_dict(value, TOPOLOGY_KEYS, "qualification topology")
    _require(topology["schema"] == TOPOLOGY_SCHEMA, "qualification topology schema differs")
    expected = topology_record(topology["cpus"])
    _require(exact_json_equal(topology, expected), "qualification topology is not canonical")
    return copy.deepcopy(expected)


def _allowed_cpu_list(value: Any) -> list[int]:
    allowed = _canonical_int_list(value, label="allowed CPU set at launch")
    _require(type(value) is list and exact_json_equal(value, allowed),
             "allowed CPU set at launch is not canonical")
    return allowed


def derive_candidate_pairs(
    policy_value: Any,
    topology_value: Any,
    allowed_cpu_set_at_launch: Sequence[int],
) -> list[dict[str, Any]]:
    policy = validate_qualification_policy(policy_value)
    topology = validate_topology(topology_value)
    allowed = _canonical_int_list(
        allowed_cpu_set_at_launch, label="allowed CPU set at launch",
    )
    by_cpu = {entry["cpu"]: entry for entry in topology["cpus"]}
    for cpu in allowed:
        _require(cpu in by_cpu and by_cpu[cpu]["online"],
                 "allowed CPU set contains a missing or offline CPU")
    excluded_by_primary = {
        item["benchmark_cpu"]: item for item in policy["excluded_pairs"]
    }
    candidates: list[dict[str, Any]] = []
    seen_physical_pairs: set[frozenset[int]] = set()
    for primary in policy["candidate_primary_cpus"]:
        _require(primary in by_cpu, "candidate primary is missing from topology")
        primary_entry = by_cpu[primary]
        _require(len(primary_entry["thread_siblings"]) == 2,
                 "candidate primary is not in an exact SMT2 core")
        sibling = next(
            cpu for cpu in primary_entry["thread_siblings"] if cpu != primary
        )
        physical = frozenset((primary, sibling))
        _require(physical not in seen_physical_pairs,
                 "candidate primary policy names one physical core more than once")
        seen_physical_pairs.add(physical)
        exclusion = excluded_by_primary.get(primary)
        if exclusion is not None:
            _require(exclusion["reserved_sibling"] == sibling,
                     "excluded pair disagrees with retained topology")
            continue
        if not (primary_entry["online"] and by_cpu[sibling]["online"] and
                primary in allowed and sibling in allowed):
            continue
        candidates.append({
            "benchmark_cpu": primary,
            "reserved_sibling": sibling,
            "domain_cpus": list(primary_entry["domain_cpus"]),
        })
    return copy.deepcopy(candidates)


def _counter_fields_record(value: Any, label: str) -> dict[str, int]:
    _require(type(value) is dict and set(value) == set(COUNTER_FIELDS),
             f"{label} has an unexpected counter-field set")
    return {
        field: _bounded_int(value[field], 0, MAX_COUNTER, f"{label} {field}")
        for field in COUNTER_FIELDS
    }


def shared_snapshot_record(
    *,
    read_started_monotonic_ns: int,
    read_finished_monotonic_ns: int,
    counters: Mapping[int, Mapping[str, int]],
) -> dict[str, Any]:
    started = _bounded_int(
        read_started_monotonic_ns, 0, MAX_MONOTONIC_NS,
        "snapshot read-start timestamp",
    )
    finished = _bounded_int(
        read_finished_monotonic_ns, 0, MAX_MONOTONIC_NS,
        "snapshot read-finish timestamp",
    )
    _require(started <= finished, "snapshot read timestamps are reversed")
    _require(type(counters) is dict and len(counters) <= MAX_CPU_COUNT,
             "snapshot counters are not a bounded mapping")
    parsed: list[tuple[int, dict[str, int]]] = []
    for raw_cpu, raw_fields in counters.items():
        cpu = _cpu_id(raw_cpu, "snapshot counter CPU")
        fields = _counter_fields_record(raw_fields, f"snapshot CPU {cpu}")
        parsed.append((cpu, fields))
    parsed.sort(key=lambda item: item[0])
    _require(len({cpu for cpu, _ in parsed}) == len(parsed),
             "snapshot contains a duplicate logical CPU")
    return {
        "cpus": {
            str(cpu): {
                "cpu": cpu,
                "fields": fields,
                "total_jiffies": _bounded_int(
                    sum(fields.values()), 0, MAX_COUNTER,
                    f"snapshot CPU {cpu} total",
                ),
            }
            for cpu, fields in parsed
        },
        "read_started_monotonic_ns": started,
        "read_finished_monotonic_ns": finished,
    }


def _validate_shared_snapshot(value: Any, expected_cpus: Sequence[int], label: str) -> dict[str, Any]:
    snapshot = _require_dict(value, SNAPSHOT_KEYS, label)
    expected_cpu_list = _canonical_int_list(expected_cpus, label=f"{label} expected CPUs")
    _require(type(snapshot["cpus"]) is dict, f"{label} CPU table is not an object")
    _require(set(snapshot["cpus"]) == {str(cpu) for cpu in expected_cpu_list},
             f"{label} CPU table differs from the observed CPU universe")
    counters: dict[int, dict[str, int]] = {}
    for cpu in expected_cpu_list:
        key = str(cpu)
        record = _require_dict(
            snapshot["cpus"][key], CPU_COUNTER_KEYS, f"{label} CPU {cpu}",
        )
        _require(record["cpu"] == cpu and type(record["cpu"]) is int,
                 f"{label} CPU identity differs")
        fields = _counter_fields_record(record["fields"], f"{label} CPU {cpu}")
        total = _bounded_int(
            record["total_jiffies"], 0, MAX_COUNTER, f"{label} CPU {cpu} total",
        )
        _require(total == sum(fields.values()), f"{label} CPU total is inconsistent")
        counters[cpu] = fields
    expected = shared_snapshot_record(
        read_started_monotonic_ns=snapshot["read_started_monotonic_ns"],
        read_finished_monotonic_ns=snapshot["read_finished_monotonic_ns"],
        counters=counters,
    )
    _require(exact_json_equal(snapshot, expected), f"{label} is not canonical")
    return expected


def _counter_delta(before: Mapping[str, Any], after: Mapping[str, Any], label: str) -> dict[str, Any]:
    _require_dict(before, CPU_COUNTER_KEYS, f"{label} before counter")
    _require_dict(after, CPU_COUNTER_KEYS, f"{label} after counter")
    cpu = _cpu_id(before["cpu"], f"{label} CPU")
    _require(type(after["cpu"]) is int and after["cpu"] == cpu,
             f"{label} CPU identity changed")
    before_fields = _counter_fields_record(before["fields"], f"{label} before")
    after_fields = _counter_fields_record(after["fields"], f"{label} after")
    deltas: dict[str, int] = {}
    for field in COUNTER_FIELDS:
        _require(after_fields[field] >= before_fields[field],
                 f"counter_epoch_regressed: {label} {field}")
        deltas[field] = after_fields[field] - before_fields[field]
    idle = sum(deltas[field] for field in IDLE_FIELDS)
    nonidle = sum(deltas[field] for field in NONIDLE_FIELDS)
    total = sum(deltas.values())
    _bounded_int(idle, 0, MAX_COUNTER, f"{label} idle delta")
    _bounded_int(nonidle, 0, MAX_COUNTER, f"{label} nonidle delta")
    _bounded_int(total, 0, MAX_COUNTER, f"{label} total delta")
    return {
        "cpu": cpu,
        "fields": deltas,
        "idle_jiffies": idle,
        "nonidle_jiffies": nonidle,
        "total_jiffies": total,
    }


def _tick_cap(elapsed_ns: int, ticks_per_second: int) -> int:
    return (elapsed_ns * ticks_per_second) // 1_000_000_000 + 1


def _pair_observation_window_record(
    policy: dict[str, Any],
    *,
    phase: str,
    index: int,
    before: Any,
    after: Any,
    observed_cpus: Sequence[int],
) -> dict[str, Any]:
    _require(type(phase) is str and phase in WINDOW_PHASES,
             "qualification window phase differs")
    window_index = _bounded_int(index, 0, MAX_WINDOWS - 1, "qualification window index")
    cpus = _canonical_int_list(observed_cpus, label="observed CPUs")
    before_record = _validate_shared_snapshot(before, cpus, "window before snapshot")
    after_record = _validate_shared_snapshot(after, cpus, "window after snapshot")
    narrow_ns = (
        after_record["read_started_monotonic_ns"] -
        before_record["read_finished_monotonic_ns"]
    )
    wide_ns = (
        after_record["read_finished_monotonic_ns"] -
        before_record["read_started_monotonic_ns"]
    )
    _require(narrow_ns > 0, "qualification window has no positive narrowed interval")
    _require(wide_ns >= narrow_ns, "qualification window read intervals overlap")
    ticks = policy["clock_ticks_per_second"]
    narrow_cap = _tick_cap(narrow_ns, ticks)
    wide_cap = _tick_cap(wide_ns, ticks)
    lower = narrow_cap - LIVENESS_LOWER_SLACK_JIFFIES
    delta: dict[str, Any] = {}
    nonidle_cpus: list[int] = []
    not_live_cpus: list[int] = []
    for cpu in cpus:
        cpu_delta = _counter_delta(
            before_record["cpus"][str(cpu)], after_record["cpus"][str(cpu)],
            f"qualification window CPU {cpu}",
        )
        delta[str(cpu)] = cpu_delta
        if cpu_delta["nonidle_jiffies"] > NONIDLE_LIMIT_JIFFIES:
            nonidle_cpus.append(cpu)
        if not (
            lower >= 1 and lower <= cpu_delta["total_jiffies"] <= wide_cap
        ):
            not_live_cpus.append(cpu)
    return {
        "schema": WINDOW_SCHEMA,
        "phase": phase,
        "index": window_index,
        "before": copy.deepcopy(before_record),
        "after": copy.deepcopy(after_record),
        "window_ns": narrow_ns,
        "wide_window_ns": wide_ns,
        "elapsed_tick_cap": narrow_cap,
        "wide_elapsed_tick_cap": wide_cap,
        "liveness_lower_bound_jiffies": lower,
        "delta": delta,
        "nonidle_cpus": nonidle_cpus,
        "not_live_cpus": not_live_cpus,
        "clock_ticks_per_second": ticks,
        "policy": _observation_policy(),
    }


def pair_observation_window_record(
    policy_value: Any,
    *,
    phase: str,
    index: int,
    before: Any,
    after: Any,
    observed_cpus: Sequence[int],
) -> dict[str, Any]:
    return _pair_observation_window_record(
        validate_qualification_policy(policy_value),
        phase=phase,
        index=index,
        before=before,
        after=after,
        observed_cpus=observed_cpus,
    )


def _validate_pair_observation_window(
    value: Any,
    *,
    policy: dict[str, Any],
    expected_phase: str,
    expected_index: int,
    expected_cpus: Sequence[int],
) -> dict[str, Any]:
    window = _require_dict(value, WINDOW_KEYS, "pair observation window")
    _require(window["schema"] == WINDOW_SCHEMA, "pair observation window schema differs")
    _require_dict(window["policy"], WINDOW_POLICY_KEYS, "window observation policy")
    expected = _pair_observation_window_record(
        policy,
        phase=expected_phase,
        index=expected_index,
        before=window["before"],
        after=window["after"],
        observed_cpus=expected_cpus,
    )
    _require(exact_json_equal(window, expected),
             "pair observation window differs from fixed-point recomputation")
    return copy.deepcopy(expected)


def validate_pair_observation_window(
    value: Any,
    *,
    expected_policy: Any,
    expected_phase: str,
    expected_index: int,
    expected_cpus: Sequence[int],
) -> dict[str, Any]:
    return _validate_pair_observation_window(
        value,
        policy=validate_qualification_policy(expected_policy),
        expected_phase=expected_phase,
        expected_index=expected_index,
        expected_cpus=expected_cpus,
    )


def _canonical_frozen_pair(value: Any, label: str) -> dict[str, int] | None:
    if value is None:
        return None
    return _validate_pair(value, label)


def _candidate_pair_identity(value: Mapping[str, Any]) -> dict[str, int]:
    return {
        "benchmark_cpu": value["benchmark_cpu"],
        "reserved_sibling": value["reserved_sibling"],
    }


def _static_candidate_identities(
    policy: Mapping[str, Any], topology: Mapping[str, Any],
) -> list[dict[str, int]]:
    by_cpu = {entry["cpu"]: entry for entry in topology["cpus"]}
    excluded = {entry["benchmark_cpu"] for entry in policy["excluded_pairs"]}
    identities: list[dict[str, int]] = []
    seen_physical_pairs: set[frozenset[int]] = set()
    for primary in policy["candidate_primary_cpus"]:
        _require(primary in by_cpu, "candidate primary is missing from topology")
        entry = by_cpu[primary]
        _require(len(entry["thread_siblings"]) == 2,
                 "candidate primary is not in an exact SMT2 core")
        sibling = next(cpu for cpu in entry["thread_siblings"] if cpu != primary)
        physical = frozenset((primary, sibling))
        _require(physical not in seen_physical_pairs,
                 "candidate primary policy names one physical core more than once")
        seen_physical_pairs.add(physical)
        if primary not in excluded:
            identities.append({
                "benchmark_cpu": primary,
                "reserved_sibling": sibling,
            })
    return identities


def pair_qualification_scan_record(
    policy_value: Any,
    *,
    allowed_cpu_set_at_launch: Sequence[int],
    topology_before: Any,
    topology_after: Any,
    snapshots: Sequence[Any],
    frozen_pair_from_prior_attempt: Any,
) -> dict[str, Any]:
    policy = validate_qualification_policy(policy_value)
    before_topology = validate_topology(topology_before)
    after_topology = validate_topology(topology_after)
    _require(exact_json_equal(before_topology, after_topology),
             "qualification topology changed across the scan")
    allowed = _canonical_int_list(
        allowed_cpu_set_at_launch, label="allowed CPU set at launch",
    )
    candidates = derive_candidate_pairs(policy, before_topology, allowed)
    observed_cpus = sorted({
        cpu for candidate in candidates for cpu in candidate["domain_cpus"]
    })
    _require(type(snapshots) in (list, tuple), "scan snapshots are not a sequence")
    _require(2 <= len(snapshots) <= MAX_WINDOWS + 1,
             "scan does not contain a bounded N+1 snapshot chain")
    endpoint_records = [
        _validate_shared_snapshot(value, observed_cpus, f"scan snapshot {index}")
        for index, value in enumerate(snapshots)
    ]
    windows = [
        _pair_observation_window_record(
            policy,
            phase="qualification",
            index=index,
            before=endpoint_records[index],
            after=endpoint_records[index + 1],
            observed_cpus=observed_cpus,
        )
        for index in range(len(endpoint_records) - 1)
    ]
    first = endpoint_records[0]
    last = endpoint_records[-1]
    narrow_scan_ns = (
        last["read_started_monotonic_ns"] - first["read_finished_monotonic_ns"]
    )
    wide_scan_ns = (
        last["read_finished_monotonic_ns"] - first["read_started_monotonic_ns"]
    )
    _require(narrow_scan_ns > 0 and wide_scan_ns >= narrow_scan_ns,
             "qualification scan interval is invalid")
    ticks = policy["clock_ticks_per_second"]
    overall_cap = _tick_cap(narrow_scan_ns, ticks)
    overall_wide_cap = _tick_cap(wide_scan_ns, ticks)
    overall_lower = overall_cap - LIVENESS_LOWER_SLACK_JIFFIES
    aggregates: dict[str, Any] = {}
    nonidle_sets = [set(window["nonidle_cpus"]) for window in windows]
    not_live_sets = [set(window["not_live_cpus"]) for window in windows]
    for cpu in observed_cpus:
        aggregate = _counter_delta(
            first["cpus"][str(cpu)], last["cpus"][str(cpu)],
            f"qualification aggregate CPU {cpu}",
        )
        summed_fields = {
            field: sum(window["delta"][str(cpu)]["fields"][field] for window in windows)
            for field in COUNTER_FIELDS
        }
        _require(summed_fields == aggregate["fields"],
                 "qualification window chain does not telescope")
        nonidle_window_count = sum(
            cpu in nonidle for nonidle in nonidle_sets
        )
        not_live_window_count = sum(
            cpu in not_live for not_live in not_live_sets
        )
        aggregates[str(cpu)] = {
            **aggregate,
            "nonidle_window_count": nonidle_window_count,
            "not_live_window_count": not_live_window_count,
            "overall_live": (
                overall_lower >= 1 and
                overall_lower <= aggregate["total_jiffies"] <= overall_wide_cap
            ),
        }
    candidate_results: list[dict[str, Any]] = []
    eligible_pairs: list[dict[str, int]] = []
    for candidate in candidates:
        pair_cpus = (candidate["benchmark_cpu"], candidate["reserved_sibling"])
        domain_cpus = candidate["domain_cpus"]
        pair_idle = all(
            aggregates[str(cpu)]["nonidle_jiffies"] == 0 and
            aggregates[str(cpu)]["nonidle_window_count"] == 0
            for cpu in pair_cpus
        )
        pair_live = all(
            aggregates[str(cpu)]["overall_live"] and
            aggregates[str(cpu)]["not_live_window_count"] == 0
            for cpu in pair_cpus
        )
        domain_idle = all(
            aggregates[str(cpu)]["nonidle_jiffies"] == 0 and
            aggregates[str(cpu)]["nonidle_window_count"] == 0
            for cpu in domain_cpus
        )
        domain_live = all(
            aggregates[str(cpu)]["overall_live"] and
            aggregates[str(cpu)]["not_live_window_count"] == 0
            for cpu in domain_cpus
        )
        qualified = pair_idle and pair_live and (
            policy["domain_mode"] == "pair-only" or (domain_idle and domain_live)
        )
        result = {
            **candidate,
            "pair_idle": pair_idle,
            "pair_live": pair_live,
            "domain_idle": domain_idle,
            "domain_live": domain_live,
            "qualified": qualified,
        }
        candidate_results.append(result)
        if qualified:
            eligible_pairs.append(_candidate_pair_identity(candidate))
    frozen = _canonical_frozen_pair(
        frozen_pair_from_prior_attempt, "frozen pair from prior attempt",
    )
    if frozen is not None:
        _require(frozen in _static_candidate_identities(policy, before_topology),
                 "frozen pair is not a static policy candidate")
        if frozen in eligible_pairs:
            selected = copy.deepcopy(frozen)
            selection_status = "selected-frozen-pair"
        else:
            selected = None
            selection_status = "frozen-pair-did-not-requalify"
    elif eligible_pairs:
        selected = copy.deepcopy(eligible_pairs[0])
        selection_status = "selected-lowest-primary"
    else:
        selected = None
        selection_status = "no-candidate-pair-qualified"
    max_read_ns = max(
        snapshot["read_finished_monotonic_ns"] -
        snapshot["read_started_monotonic_ns"]
        for snapshot in endpoint_records
    )
    return {
        "schema": SCAN_SCHEMA,
        "policy": copy.deepcopy(policy),
        "allowed_cpu_set_at_launch": allowed,
        "topology_before": copy.deepcopy(before_topology),
        "topology_after": copy.deepcopy(after_topology),
        "topology_reread_identical": True,
        "windows": windows,
        "cpu_aggregates": aggregates,
        "candidate_pairs": candidate_results,
        "eligible_pairs": eligible_pairs,
        "selected": selected,
        "frozen_pair_from_prior_attempt": frozen,
        "selection_status": selection_status,
        "scan_started_monotonic_ns": first["read_started_monotonic_ns"],
        "scan_finished_monotonic_ns": last["read_finished_monotonic_ns"],
        "realized_scan_total_ns": narrow_scan_ns,
        "wide_scan_total_ns": wide_scan_ns,
        "overall_elapsed_tick_cap": overall_cap,
        "overall_wide_elapsed_tick_cap": overall_wide_cap,
        "overall_liveness_lower_bound_jiffies": overall_lower,
        "scan_window_count": len(windows),
        "max_proc_stat_read_ns": max_read_ns,
        "excluded_pair_count": len(policy["excluded_pairs"]),
        "candidate_pair_count": len(candidate_results),
        "eligible_pair_count": len(eligible_pairs),
        "candidate_timing_performed": False,
    }


def validate_pair_qualification_scan(
    value: Any,
    *,
    expected_policy: Any,
    expected_frozen_pair: Any,
) -> dict[str, Any]:
    scan = _require_dict(value, SCAN_KEYS, "pair qualification scan")
    _require(scan["schema"] == SCAN_SCHEMA, "pair qualification scan schema differs")
    policy = validate_qualification_policy(expected_policy)
    _require(exact_json_equal(scan["policy"], policy),
             "pair qualification scan policy differs from preregistration")
    frozen = _canonical_frozen_pair(expected_frozen_pair, "expected frozen pair")
    _require(exact_json_equal(scan["frozen_pair_from_prior_attempt"], frozen),
             "pair qualification scan erased or changed the frozen pair")
    before_topology = validate_topology(scan["topology_before"])
    after_topology = validate_topology(scan["topology_after"])
    allowed = _allowed_cpu_list(scan["allowed_cpu_set_at_launch"])
    candidates = derive_candidate_pairs(policy, before_topology, allowed)
    observed_cpus = sorted({
        cpu for candidate in candidates for cpu in candidate["domain_cpus"]
    })
    _require(type(scan["windows"]) is list and 0 < len(scan["windows"]) <= MAX_WINDOWS,
             "pair qualification windows are not a bounded non-empty list")
    windows: list[dict[str, Any]] = []
    endpoints: list[dict[str, Any]] = []
    for index, raw_window in enumerate(scan["windows"]):
        window = _validate_pair_observation_window(
            raw_window,
            policy=policy,
            expected_phase="qualification",
            expected_index=index,
            expected_cpus=observed_cpus,
        )
        if index == 0:
            endpoints.append(window["before"])
        else:
            _require(exact_json_equal(windows[-1]["after"], window["before"]),
                     "qualification windows do not share an exact N+1 endpoint")
        endpoints.append(window["after"])
        windows.append(window)
    expected = pair_qualification_scan_record(
        policy,
        allowed_cpu_set_at_launch=allowed,
        topology_before=before_topology,
        topology_after=after_topology,
        snapshots=endpoints,
        frozen_pair_from_prior_attempt=frozen,
    )
    _require(exact_json_equal(scan, expected),
             "pair qualification scan differs from fixed-point recomputation")
    return copy.deepcopy(expected)


def load_pair_qualification_scan(
    data: bytes,
    *,
    expected_policy: Any,
    expected_frozen_pair: Any,
) -> dict[str, Any]:
    return validate_pair_qualification_scan(
        strict_json_loads(data, "pair qualification scan"),
        expected_policy=expected_policy,
        expected_frozen_pair=expected_frozen_pair,
    )


def selected_pair_from_scan(
    value: Any,
    *,
    expected_policy: Any,
    expected_frozen_pair: Any,
) -> dict[str, int]:
    validated = validate_pair_qualification_scan(
        value,
        expected_policy=expected_policy,
        expected_frozen_pair=expected_frozen_pair,
    )
    status = validated["selection_status"]
    if status == "no-candidate-pair-qualified":
        _fail("no_candidate_pair_qualified")
    if status == "frozen-pair-did-not-requalify":
        _fail("frozen_pair_did_not_requalify")
    _require(status in ("selected-lowest-primary", "selected-frozen-pair"),
             "selected-pair status differs")
    return copy.deepcopy(_validate_pair(validated["selected"], "selected pair"))


__all__ = (
    "CANDIDATE_RULE", "COUNTER_FIELDS", "COUNTER_SOURCE", "DOMAIN_MODES",
    "FROZEN_PAIR_RULE", "IDLE_FIELDS", "LIVENESS_RULE", "NONIDLE_FIELDS",
    "POLICY_SCHEMA", "QualificationError", "SCAN_SCHEMA", "SELECTION_RULE",
    "TOPOLOGY_SCHEMA", "WINDOW_SCHEMA", "canonical_json_bytes", "canonical_sha256",
    "derive_candidate_pairs", "exact_json_equal", "load_pair_qualification_scan",
    "pair_observation_window_record", "pair_qualification_scan_record",
    "qualification_policy_record", "selected_pair_from_scan",
    "shared_snapshot_record", "strict_json_loads", "topology_record",
    "validate_pair_observation_window", "validate_pair_qualification_scan",
    "validate_qualification_policy", "validate_topology",
)
