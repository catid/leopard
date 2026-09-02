#!/usr/bin/env python3
"""Pure fixed-point contract for conditioned-passive v19 evidence.

This module binds the prospectively frozen qualification and bridge records to
the two-attempt state machine, the first retained invocation handoff, the three
exhausted v18 failure envelopes, the shared-host claim ceiling, and the
``NOT_PROMOTED`` terminal.  It performs no host access, file I/O, process
launch, affinity change, benchmark timing, or evidence publication.
"""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import sys
from typing import Any, Mapping, NoReturn


def _load_local_module(module_name: str, filename: str) -> Any:
    expected = Path(__file__).resolve().with_name(filename)
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        if Path(getattr(loaded, "__file__", "")).resolve() != expected:
            raise RuntimeError(f"{filename} came from another path")
        return loaded
    specification = importlib.util.spec_from_file_location(module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    if Path(getattr(module, "__file__", "")).resolve() != expected:
        raise RuntimeError(f"{filename} resolved to another path")
    return module


bridge = _load_local_module(
    "leopard2_pair_qualification_bridge_for_v19_contract",
    "pair_qualification_bridge_contract.py")
verifier = _load_local_module(
    "leopard2_pair_qualification_verifier_for_v19_contract",
    "pair_qualification_verify.py")
contract = bridge.contract

SCHEMA = "leopard2-main-compare-pair-qualification/v1"
ATTEMPT_SCHEMA = "leopard2-main-compare-pair-qualified-attempt/v1"
HANDOFF_SCHEMA = "leopard2-main-compare-first-window-handoff/v1"
V18_FAILURE_LINEAGE_SCHEMA = \
    "leopard2-main-compare-v18-failure-lineage/v1"
ATTEMPT_BUDGET = 2
QUALIFICATION_WINDOW_COUNT = 12
QUALIFICATION_NOMINAL_WINDOW_NS = 7_000_000_000
BRIDGE_WINDOW_COUNT = 2
BRIDGE_NOMINAL_WINDOW_NS = 1_000_000_000
MAXIMUM_HANDOFF_ELAPSED_NS = 5_000_000_000
SUCCESS_TERMINAL = "NOT_PROMOTED"

STAGES = (
    "identity", "acquired", "bridged", "handoff", "campaign", "complete",
)
RECORD_STATUSES = ("in-progress", "failed", "complete")
FAILURE_TERMINALS = frozenset((
    "identity-failed",
    "no-candidate-pair-qualified",
    "frozen-pair-did-not-requalify",
    "bridge-not-accepted",
    "first-window-handoff-unavailable",
    "first-window-handoff-late",
    "first-window-handoff-selected-pair-nonidle",
    "campaign-rejected",
))

ATTEMPT_KEYS = frozenset((
    "schema", "attempt", "attempt_budget", "frozen_pair_from_prior_attempt",
    "prior_attempt_failure_sha256", "prior_attempt_acquisition_sha256",
    "prior_attempt_selection_status", "prior_attempt_selected_pair",
    "fresh_selection_permitted", "selection_rule", "frozen_pair_rule",
))
QUALIFICATION_GEOMETRY_KEYS = frozenset((
    "window_count", "nominal_window_ns",
))
HANDOFF_KEYS = frozenset((
    "schema", "selected_pair", "bridge_tail_sha256",
    "bridge_tail_read_finished_monotonic_ns", "first_window_before",
    "first_window_before_sha256", "first_window_before_read_started_monotonic_ns",
    "handoff_elapsed_ns", "maximum_handoff_elapsed_ns",
    "selected_pair_nonidle_delta", "accepted", "failure_terminal",
))
PAIR_WINDOW_BEFORE_KEYS = frozenset((
    "benchmark_cpu", "monotonic_ns", "read_finished_monotonic_ns",
    "read_started_monotonic_ns", "reserved_sibling",
))
CPU_COUNTER_KEYS = frozenset(("cpu", "fields", "total_jiffies"))
NONIDLE_DELTA_KEYS = frozenset(("benchmark_cpu", "reserved_sibling"))
LINEAGE_ENTRY_KEYS = frozenset((
    "attempt", "envelope", "envelope_sha256sums_sha256",
))
LINEAGE_KEYS = frozenset((
    "schema", "source_commit", "source_tree", "attempts",
))
CLAIM_CEILING_KEYS = frozenset((
    "promotion_eligible", "host_exclusivity_proved",
    "whole_campaign_interval_observed", "causal_performance_claim_allowed",
))
RECORD_KEYS = frozenset((
    "schema", "stage", "record_status", "attempt", "policy",
    "policy_sha256", "host_identity_sha256", "qualification_geometry",
    "qualification_geometry_sha256", "acquisition", "acquisition_sha256",
    "selected_pair", "selection_status", "bridge", "bridge_sha256",
    "bridge_geometry_sha256", "first_window_handoff",
    "v18_failure_lineage", "v18_failure_lineage_sha256",
    "shared_host_claim_ceiling", "terminal", "host_mutation_performed",
    "candidate_timing_performed",
))

# The v18 ancestry is content-addressed historical authority, not a live-path
# dependency.  Its three envelope paths are disclosure labels bound to the
# exact retained SHA256SUMS digests below.  Separate v18 compatibility gates
# replay the sealed copies; v19 replay intentionally remains possible if the
# archival paths are offline after those bytes have been pinned.
_V18_SOURCE_COMMIT = "c8f825d0a033d31d220b0ebce9cc8871e8c2fc6d"
_V18_SOURCE_TREE = "2c17a0a7bcea20274d2593cb204442c4c817e464"
_V18_FAILURES = (
    (1, ".research/leopard-79h/c8f825d-v18-passive-main-a1",
     "ce65c3a49ef1c1d89ba51ea03d0af4742d6790e6f2ea2662917d9ef9a9d945d7"),
    (2, ".research/leopard-79h/c8f825d-v18-passive-main-a2",
     "a1bf0eda157c251f33f7260ebd76931d88054d460bd07a97bcba2811384b2c10"),
    (3, ".research/leopard-79h/c8f825d-v18-passive-main-a3",
     "fe5b40cc98753cbd794ee019cb0e2643d0ccee0aca4c5fd7b2e0b27df8a86139"),
)


class PairQualifiedV19Error(bridge.BridgeError):
    """The conditioned-passive v19 fixed point was not established."""


def _fail(message: str) -> NoReturn:
    raise PairQualifiedV19Error(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _hex64(value: Any, label: str) -> str:
    _require(type(value) is str and len(value) == 64 and
             all(character in "0123456789abcdef" for character in value),
             f"{label} is not exact lowercase SHA-256")
    return value


def _bounded_int(value: Any, minimum: int, maximum: int, label: str) -> int:
    _require(type(value) is int and minimum <= value <= maximum,
             f"{label} is not an exact bounded integer")
    return value


def _pair(value: Any, label: str) -> dict[str, int] | None:
    if value is None:
        return None
    _require(type(value) is dict and set(value) == {
        "benchmark_cpu", "reserved_sibling",
    }, f"{label} is not an exact pair object")
    benchmark = _bounded_int(
        value["benchmark_cpu"], 0, contract.MAX_CPU_ID,
        f"{label} benchmark CPU")
    sibling = _bounded_int(
        value["reserved_sibling"], 0, contract.MAX_CPU_ID,
        f"{label} reserved sibling")
    _require(benchmark != sibling, f"{label} repeats one logical CPU")
    return {"benchmark_cpu": benchmark, "reserved_sibling": sibling}


def _claim_ceiling() -> dict[str, bool]:
    return {
        "promotion_eligible": False,
        "host_exclusivity_proved": False,
        "whole_campaign_interval_observed": False,
        "causal_performance_claim_allowed": False,
    }


def qualification_geometry_record() -> dict[str, int]:
    return {
        "window_count": QUALIFICATION_WINDOW_COUNT,
        "nominal_window_ns": QUALIFICATION_NOMINAL_WINDOW_NS,
    }


def bridge_geometry_record() -> dict[str, int]:
    return bridge.bridge_geometry_record(
        minimum_window_count=BRIDGE_WINDOW_COUNT,
        maximum_window_count=BRIDGE_WINDOW_COUNT,
        nominal_window_ns=BRIDGE_NOMINAL_WINDOW_NS,
        maximum_handoff_elapsed_ns=MAXIMUM_HANDOFF_ELAPSED_NS)


def v18_failure_lineage_record() -> dict[str, Any]:
    return {
        "schema": V18_FAILURE_LINEAGE_SCHEMA,
        "source_commit": _V18_SOURCE_COMMIT,
        "source_tree": _V18_SOURCE_TREE,
        "attempts": [
            {
                "attempt": attempt,
                "envelope": envelope,
                "envelope_sha256sums_sha256": digest,
            }
            for attempt, envelope, digest in _V18_FAILURES
        ],
    }


def validate_v18_failure_lineage(value: Any) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == LINEAGE_KEYS,
             "v18 failure lineage has an unexpected key set")
    _require(value["schema"] == V18_FAILURE_LINEAGE_SCHEMA and
             value["source_commit"] == _V18_SOURCE_COMMIT and
             value["source_tree"] == _V18_SOURCE_TREE and
             type(value["attempts"]) is list and
             len(value["attempts"]) == len(_V18_FAILURES),
             "v18 failure lineage identity differs")
    for entry in value["attempts"]:
        _require(type(entry) is dict and set(entry) == LINEAGE_ENTRY_KEYS,
                 "v18 failure lineage entry has an unexpected key set")
        _bounded_int(entry["attempt"], 1, 3, "v18 failure attempt")
        _hex64(entry["envelope_sha256sums_sha256"],
               "v18 failure envelope hash")
    expected = v18_failure_lineage_record()
    _require(contract.exact_json_equal(value, expected),
             "v18 failure lineage differs from the exhausted envelopes")
    return copy.deepcopy(expected)


def pair_qualified_attempt_record(
    *, attempt: int, prior_attempt_failure_sha256: Any = None,
    prior_attempt_acquisition_sha256: Any = None,
    prior_attempt_selection_status: Any = None,
    prior_attempt_selected_pair: Any = None,
) -> dict[str, Any]:
    number = _bounded_int(attempt, 1, ATTEMPT_BUDGET, "v19 attempt")
    selected = _pair(prior_attempt_selected_pair, "prior selected pair")
    if number == 1:
        _require(prior_attempt_failure_sha256 is None and
                 prior_attempt_acquisition_sha256 is None and
                 prior_attempt_selection_status is None and selected is None,
                 "attempt 1 cannot carry prior-attempt state")
        failure_hash = acquisition_hash = selection_status = None
        frozen = None
        fresh = True
    else:
        failure_hash = _hex64(
            prior_attempt_failure_sha256, "prior attempt failure hash")
        _require(prior_attempt_selection_status in (
            "not-acquired", "no-candidate-pair-qualified",
            "selected-lowest-primary",
        ), "attempt 2 prior selection status differs")
        selection_status = prior_attempt_selection_status
        if selection_status == "not-acquired":
            _require(prior_attempt_acquisition_sha256 is None and selected is None,
                     "not-acquired prior attempt carries acquisition state")
            acquisition_hash = None
        else:
            acquisition_hash = _hex64(
                prior_attempt_acquisition_sha256,
                "prior attempt acquisition hash")
        if selection_status == "selected-lowest-primary":
            _require(selected is not None,
                     "selected prior attempt lacks its frozen pair")
            frozen = copy.deepcopy(selected)
            fresh = False
        else:
            _require(selected is None,
                     "unselected prior attempt carries a pair")
            frozen = None
            fresh = True
    return {
        "schema": ATTEMPT_SCHEMA,
        "attempt": number,
        "attempt_budget": ATTEMPT_BUDGET,
        "frozen_pair_from_prior_attempt": frozen,
        "prior_attempt_failure_sha256": failure_hash,
        "prior_attempt_acquisition_sha256": acquisition_hash,
        "prior_attempt_selection_status": selection_status,
        "prior_attempt_selected_pair": copy.deepcopy(selected),
        "fresh_selection_permitted": fresh,
        "selection_rule": contract.SELECTION_RULE,
        "frozen_pair_rule": contract.FROZEN_PAIR_RULE,
    }


def validate_pair_qualified_attempt(value: Any) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == ATTEMPT_KEYS,
             "v19 attempt record has an unexpected key set")
    expected = pair_qualified_attempt_record(
        attempt=value["attempt"],
        prior_attempt_failure_sha256=value["prior_attempt_failure_sha256"],
        prior_attempt_acquisition_sha256=
        value["prior_attempt_acquisition_sha256"],
        prior_attempt_selection_status=value["prior_attempt_selection_status"],
        prior_attempt_selected_pair=value["prior_attempt_selected_pair"])
    _require(contract.exact_json_equal(value, expected),
             "v19 attempt record is not canonical")
    return copy.deepcopy(expected)


def _cpu_counter(value: Any, expected_cpu: int, label: str) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == CPU_COUNTER_KEYS and
             type(value["cpu"]) is int and value["cpu"] == expected_cpu and
             type(value["fields"]) is dict and
             set(value["fields"]) == set(contract.COUNTER_FIELDS),
             f"{label} has an unexpected CPU-counter shape")
    fields: dict[str, int] = {}
    for field in contract.COUNTER_FIELDS:
        fields[field] = _bounded_int(
            value["fields"][field], 0, contract.MAX_COUNTER,
            f"{label} {field}")
    total = _bounded_int(
        value["total_jiffies"], 0, contract.MAX_COUNTER, f"{label} total")
    _require(total == sum(fields.values()), f"{label} total differs")
    expected = {"cpu": expected_cpu, "fields": fields, "total_jiffies": total}
    _require(contract.exact_json_equal(value, expected),
             f"{label} is not canonical")
    return expected


def _first_window_before(
    value: Any, selected_pair: dict[str, int],
) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == PAIR_WINDOW_BEFORE_KEYS,
             "first-window before snapshot has an unexpected key set")
    started = _bounded_int(
        value["read_started_monotonic_ns"], 0, contract.MAX_MONOTONIC_NS,
        "first-window read-start timestamp")
    finished = _bounded_int(
        value["read_finished_monotonic_ns"], started,
        contract.MAX_MONOTONIC_NS, "first-window read-finish timestamp")
    _require(type(value["monotonic_ns"]) is int and
             value["monotonic_ns"] == finished,
             "first-window monotonic boundary differs from read finish")
    expected = {
        "benchmark_cpu": _cpu_counter(
            value["benchmark_cpu"], selected_pair["benchmark_cpu"],
            "first-window benchmark CPU"),
        "monotonic_ns": finished,
        "read_finished_monotonic_ns": finished,
        "read_started_monotonic_ns": started,
        "reserved_sibling": _cpu_counter(
            value["reserved_sibling"], selected_pair["reserved_sibling"],
            "first-window reserved sibling"),
    }
    _require(contract.exact_json_equal(value, expected),
             "first-window before snapshot is not canonical")
    return copy.deepcopy(expected)


def _nonidle_delta(
    before: Mapping[str, Any], after: Mapping[str, Any], label: str,
) -> int:
    total = 0
    for field in contract.COUNTER_FIELDS:
        first = before["fields"][field]
        last = after["fields"][field]
        _require(last >= first, f"{label} {field} counter regressed")
        if field in contract.NONIDLE_FIELDS:
            total += last - first
    return _bounded_int(total, 0, contract.MAX_COUNTER, f"{label} nonidle delta")


def first_window_handoff_record(
    bridge_value: Any, selected_pair_value: Any,
    first_window_before_value: Any,
) -> dict[str, Any]:
    """Build a handoff from an already validated, accepted bridge record.

    This producer helper is not an independent bridge validator.  The public
    v19 record validator always runs the pure bridge fixed point and independent
    verifier before calling it.
    """
    selected = _pair(selected_pair_value, "handoff selected pair")
    _require(selected is not None and type(bridge_value) is dict and
             set(bridge_value) == bridge.BRIDGE_KEYS and
             bridge_value.get("schema") == bridge.BRIDGE_SCHEMA and
             bridge_value.get("bridge_accepted") is True and
             contract.exact_json_equal(bridge_value.get("selected_pair"), selected),
             "handoff bridge is not accepted for the selected pair")
    tail = bridge_value.get("campaign_presample_before")
    _require(type(tail) is dict and type(tail.get("cpus")) is dict,
             "handoff bridge tail is malformed")
    tail_hash = contract.canonical_sha256(tail)
    _require(bridge_value.get("bridge_tail_sha256") == tail_hash and
             bridge_value.get("campaign_presample_before_sha256") == tail_hash,
             "handoff bridge-tail binding differs")
    before = _first_window_before(first_window_before_value, selected)
    tail_finished = _bounded_int(
        tail.get("read_finished_monotonic_ns"), 0, contract.MAX_MONOTONIC_NS,
        "bridge-tail read-finish timestamp")
    first_started = before["read_started_monotonic_ns"]
    _require(first_started >= tail_finished,
             "first-window snapshot starts before the bridge tail")
    elapsed = first_started - tail_finished
    deltas: dict[str, int] = {}
    for role in ("benchmark_cpu", "reserved_sibling"):
        cpu = selected[role]
        key = str(cpu)
        _require(key in tail["cpus"],
                 "selected pair is absent from the bridge tail")
        tail_counter = _cpu_counter(
            tail["cpus"][key], cpu, f"bridge-tail {role}")
        deltas[role] = _nonidle_delta(
            tail_counter, before[role], f"first-window {role}")
    late = elapsed > MAXIMUM_HANDOFF_ELAPSED_NS
    nonidle = any(value != 0 for value in deltas.values())
    failure_terminal = (
        "first-window-handoff-late" if late else
        "first-window-handoff-selected-pair-nonidle" if nonidle else None
    )
    return {
        "schema": HANDOFF_SCHEMA,
        "selected_pair": copy.deepcopy(selected),
        "bridge_tail_sha256": tail_hash,
        "bridge_tail_read_finished_monotonic_ns": tail_finished,
        "first_window_before": copy.deepcopy(before),
        "first_window_before_sha256": contract.canonical_sha256(before),
        "first_window_before_read_started_monotonic_ns": first_started,
        "handoff_elapsed_ns": elapsed,
        "maximum_handoff_elapsed_ns": MAXIMUM_HANDOFF_ELAPSED_NS,
        "selected_pair_nonidle_delta": deltas,
        "accepted": not late and not nonidle,
        "failure_terminal": failure_terminal,
    }


def _validated_components(
    *, policy_value: Any, expected_policy_sha256: Any,
    attempt_value: Any, acquisition_value: Any, bridge_value: Any,
    first_window_before_value: Any,
) -> tuple[
    dict[str, Any], str, dict[str, Any], dict[str, Any] | None,
    dict[str, int] | None, str | None, dict[str, Any] | None,
    dict[str, Any] | None,
]:
    policy = contract.validate_qualification_policy(policy_value)
    _require(policy["domain_mode"] == "pair-only",
             "v19 evidence requires the pair-only qualification policy")
    policy_hash = _hex64(expected_policy_sha256, "expected policy hash")
    _require(policy_hash == contract.canonical_sha256(policy),
             "expected policy hash differs from the qualification policy")
    attempt = validate_pair_qualified_attempt(attempt_value)
    frozen = attempt["frozen_pair_from_prior_attempt"]
    acquisition_record = selected = selection_status = None
    if acquisition_value is not None:
        acquisition_record = bridge.validate_acquisition_for_bridge(
            acquisition_value, expected_policy=policy,
            expected_frozen_pair=frozen,
            expected_window_count=QUALIFICATION_WINDOW_COUNT,
            expected_nominal_window_ns=QUALIFICATION_NOMINAL_WINDOW_NS)
        selection_status = acquisition_record["scan"]["selection_status"]
        selected = _pair(acquisition_record["scan"]["selected"],
                         "acquisition selected pair")
        if frozen is None:
            _require(selection_status in (
                "selected-lowest-primary", "no-candidate-pair-qualified",
            ), "fresh v19 acquisition selection status differs")
        else:
            _require(selection_status in (
                "selected-frozen-pair", "frozen-pair-did-not-requalify",
            ), "frozen v19 acquisition selection status differs")
    bridge_record = handoff = None
    if bridge_value is not None:
        _require(acquisition_record is not None and selected is not None,
                 "v19 bridge lacks a selected acquisition")
        geometry = bridge_geometry_record()
        bridge_record = bridge.validate_pair_qualification_bridge(
            bridge_value, acquisition_record, expected_policy=policy,
            expected_frozen_pair=frozen,
            expected_acquisition_window_count=QUALIFICATION_WINDOW_COUNT,
            expected_acquisition_nominal_window_ns=
            QUALIFICATION_NOMINAL_WINDOW_NS,
            expected_bridge_geometry=geometry)
        _require(bridge_record["bridge_accepted"] is True,
                 "v19 retains only an accepted qualification bridge")
        try:
            verdict = verifier.require_accepted_pair_qualification_bundle(
                contract.canonical_json_bytes(acquisition_record),
                contract.canonical_json_bytes(bridge_record),
                expected_policy=policy, expected_policy_sha256=policy_hash,
                expected_frozen_pair=frozen,
                expected_acquisition_window_count=QUALIFICATION_WINDOW_COUNT,
                expected_acquisition_nominal_window_ns=
                QUALIFICATION_NOMINAL_WINDOW_NS,
                expected_bridge_geometry=geometry)
        except Exception as error:
            raise PairQualifiedV19Error(
                f"v19 bridge independent verification failed: {error}") from error
        _require(contract.exact_json_equal(verdict["selected_pair"], selected),
                 "independent bridge selected another pair")
    if first_window_before_value is not None:
        _require(bridge_record is not None and selected is not None,
                 "first-window handoff lacks an accepted bridge")
        handoff = first_window_handoff_record(
            bridge_record, selected, first_window_before_value)
    return (
        policy, policy_hash, attempt, acquisition_record, selected,
        selection_status, bridge_record, handoff,
    )


def _validate_stage_status(
    *, stage: Any, record_status: Any, terminal: Any,
    acquisition: dict[str, Any] | None, selected: dict[str, int] | None,
    bridge_record: dict[str, Any] | None, handoff: dict[str, Any] | None,
) -> tuple[str, str, str | None]:
    _require(type(stage) is str and stage in STAGES, "v19 stage differs")
    _require(type(record_status) is str and record_status in RECORD_STATUSES,
             "v19 record status differs")
    required_rank = STAGES.index(stage)
    _require((acquisition is not None) == (required_rank >= 1) and
             (bridge_record is not None) == (required_rank >= 2) and
             (handoff is not None) == (required_rank >= 3),
             "v19 stage does not match its retained prefix")
    if acquisition is not None and selected is None:
        _require(stage == "acquired" and record_status == "failed",
                 "unselected acquisition escaped its failure terminal")
    if handoff is not None and not handoff["accepted"]:
        _require(stage == "handoff" and record_status == "failed",
                 "rejected handoff escaped its failure terminal")
    if record_status == "complete":
        _require(stage == "complete" and terminal == SUCCESS_TERMINAL and
                 handoff is not None and handoff["accepted"],
                 "complete v19 evidence lacks its NOT_PROMOTED terminal")
        return stage, record_status, SUCCESS_TERMINAL
    if record_status == "in-progress":
        _require(stage != "complete" and terminal is None and
                 (handoff is None or handoff["accepted"]),
                 "in-progress v19 evidence carries a terminal")
        return stage, record_status, None
    _require(type(terminal) is str and terminal in FAILURE_TERMINALS,
             "failed v19 evidence has an unknown terminal")
    expected_terminal: str
    if stage == "identity":
        expected_terminal = "identity-failed"
    elif stage == "acquired":
        if selected is None:
            selection_status = acquisition["scan"]["selection_status"]
            expected_terminal = (
                "no-candidate-pair-qualified"
                if selection_status == "no-candidate-pair-qualified"
                else "frozen-pair-did-not-requalify"
            )
        else:
            expected_terminal = "bridge-not-accepted"
    elif stage == "bridged":
        expected_terminal = "first-window-handoff-unavailable"
    elif stage == "handoff":
        _require(handoff is not None and not handoff["accepted"],
                 "failed handoff terminal lacks a rejected handoff")
        expected_terminal = handoff["failure_terminal"]
    elif stage == "campaign":
        expected_terminal = "campaign-rejected"
    else:
        _fail("complete stage cannot be a failed v19 record")
    _require(terminal == expected_terminal,
             "v19 failure terminal differs from its retained prefix")
    return stage, record_status, terminal


def pair_qualified_v19_record(
    *, stage: str, record_status: str, terminal: Any,
    policy_value: Any, expected_policy_sha256: Any,
    host_identity_sha256: Any, attempt_value: Any,
    acquisition_value: Any = None, bridge_value: Any = None,
    first_window_before_value: Any = None,
) -> dict[str, Any]:
    host_hash = _hex64(host_identity_sha256, "host identity hash")
    (
        policy, policy_hash, attempt, acquisition_record, selected,
        selection_status, bridge_record, handoff,
    ) = _validated_components(
        policy_value=policy_value,
        expected_policy_sha256=expected_policy_sha256,
        attempt_value=attempt_value, acquisition_value=acquisition_value,
        bridge_value=bridge_value,
        first_window_before_value=first_window_before_value)
    stage, record_status, terminal = _validate_stage_status(
        stage=stage, record_status=record_status, terminal=terminal,
        acquisition=acquisition_record, selected=selected,
        bridge_record=bridge_record, handoff=handoff)
    qualification_geometry = qualification_geometry_record()
    bridge_geometry = bridge_geometry_record()
    lineage = v18_failure_lineage_record()
    return {
        "schema": SCHEMA,
        "stage": stage,
        "record_status": record_status,
        "attempt": copy.deepcopy(attempt),
        "policy": copy.deepcopy(policy),
        "policy_sha256": policy_hash,
        "host_identity_sha256": host_hash,
        "qualification_geometry": qualification_geometry,
        "qualification_geometry_sha256":
        contract.canonical_sha256(qualification_geometry),
        "acquisition": copy.deepcopy(acquisition_record),
        "acquisition_sha256": (
            contract.canonical_sha256(acquisition_record)
            if acquisition_record is not None else None),
        "selected_pair": copy.deepcopy(selected),
        "selection_status": selection_status,
        "bridge": copy.deepcopy(bridge_record),
        "bridge_sha256": (
            contract.canonical_sha256(bridge_record)
            if bridge_record is not None else None),
        "bridge_geometry_sha256": contract.canonical_sha256(bridge_geometry),
        "first_window_handoff": copy.deepcopy(handoff),
        "v18_failure_lineage": lineage,
        "v18_failure_lineage_sha256": contract.canonical_sha256(lineage),
        "shared_host_claim_ceiling": _claim_ceiling(),
        "terminal": terminal,
        "host_mutation_performed": False,
        "candidate_timing_performed": False,
    }


def validate_pair_qualified_v19_record(
    value: Any, *, expected_policy: Any, expected_policy_sha256: Any,
    expected_host_identity_sha256: Any, expected_attempt: Any,
) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == RECORD_KEYS,
             "pair-qualified v19 record has an unexpected key set")
    _require(value["schema"] == SCHEMA,
             "pair-qualified v19 schema differs")
    _require(contract.exact_json_equal(value["policy"], expected_policy),
             "pair-qualified v19 policy differs from preregistration")
    _require(value["policy_sha256"] == expected_policy_sha256,
             "pair-qualified v19 policy hash differs from preregistration")
    _require(value["host_identity_sha256"] == expected_host_identity_sha256,
             "pair-qualified v19 host identity differs")
    _require(contract.exact_json_equal(value["attempt"], expected_attempt),
             "pair-qualified v19 attempt differs from expectation")
    _require(type(value["qualification_geometry"]) is dict and
             set(value["qualification_geometry"]) ==
             QUALIFICATION_GEOMETRY_KEYS and
             contract.exact_json_equal(
                 value["qualification_geometry"],
                 qualification_geometry_record()) and
             value["qualification_geometry_sha256"] ==
             contract.canonical_sha256(qualification_geometry_record()),
             "pair-qualified v19 qualification geometry differs")
    validate_v18_failure_lineage(value["v18_failure_lineage"])
    _require(value["v18_failure_lineage_sha256"] ==
             contract.canonical_sha256(v18_failure_lineage_record()),
             "pair-qualified v19 lineage hash differs")
    _require(type(value["shared_host_claim_ceiling"]) is dict and
             set(value["shared_host_claim_ceiling"]) == CLAIM_CEILING_KEYS and
             contract.exact_json_equal(
                 value["shared_host_claim_ceiling"], _claim_ceiling()),
             "pair-qualified v19 claim ceiling differs")
    _require(type(value["host_mutation_performed"]) is bool and
             value["host_mutation_performed"] is False and
             type(value["candidate_timing_performed"]) is bool and
             value["candidate_timing_performed"] is False,
             "pair-qualified v19 record reports mutation or timing")
    first_before = None
    if value["first_window_handoff"] is not None:
        _require(type(value["first_window_handoff"]) is dict and
                 set(value["first_window_handoff"]) == HANDOFF_KEYS,
                 "pair-qualified v19 handoff has an unexpected key set")
        first_before = value["first_window_handoff"]["first_window_before"]
    expected = pair_qualified_v19_record(
        stage=value["stage"], record_status=value["record_status"],
        terminal=value["terminal"], policy_value=expected_policy,
        expected_policy_sha256=expected_policy_sha256,
        host_identity_sha256=expected_host_identity_sha256,
        attempt_value=expected_attempt,
        acquisition_value=value["acquisition"],
        bridge_value=value["bridge"],
        first_window_before_value=first_before)
    _require(contract.exact_json_equal(value, expected),
             "pair-qualified v19 record differs from fixed-point recomputation")
    return copy.deepcopy(expected)


def load_pair_qualified_v19_record(
    data: bytes, *, expected_policy: Any, expected_policy_sha256: Any,
    expected_host_identity_sha256: Any, expected_attempt: Any,
) -> dict[str, Any]:
    validated = validate_pair_qualified_v19_record(
        contract.strict_json_loads(data, "pair-qualified v19 record"),
        expected_policy=expected_policy,
        expected_policy_sha256=expected_policy_sha256,
        expected_host_identity_sha256=expected_host_identity_sha256,
        expected_attempt=expected_attempt)
    _require(data == contract.canonical_json_bytes(validated),
             "pair-qualified v19 wire bytes are not canonical")
    return validated


__all__ = (
    "ATTEMPT_BUDGET", "ATTEMPT_SCHEMA", "BRIDGE_NOMINAL_WINDOW_NS",
    "BRIDGE_WINDOW_COUNT", "FAILURE_TERMINALS", "HANDOFF_SCHEMA",
    "MAXIMUM_HANDOFF_ELAPSED_NS", "PairQualifiedV19Error",
    "QUALIFICATION_NOMINAL_WINDOW_NS", "QUALIFICATION_WINDOW_COUNT", "SCHEMA",
    "SUCCESS_TERMINAL", "V18_FAILURE_LINEAGE_SCHEMA", "bridge_geometry_record",
    "first_window_handoff_record", "load_pair_qualified_v19_record",
    "pair_qualified_attempt_record", "pair_qualified_v19_record",
    "qualification_geometry_record",
    "validate_pair_qualified_attempt", "validate_pair_qualified_v19_record",
    "validate_v18_failure_lineage", "v18_failure_lineage_record",
)
