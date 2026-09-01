#!/usr/bin/env python3
"""Generation-3 K65 logical plan/state contract; this file cannot run benchmarks.

The only CLI mode is ``--plan-only``.  It emits every one of the 1,650 frozen
logical child argument tails without accepting executable paths or a CPU pair,
launching a child, reading host topology, acquiring a lease, or touching an
evidence lane.  The separately hashed acquisition runner may materialize these
logical roles only from an atomically durable ARMED record and sealed handles.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
from pathlib import Path, PurePosixPath
import re
import sys
from typing import Any, Mapping, NoReturn, Sequence


HERE = Path(__file__).resolve().parent
PREREGISTRATION_PATH = HERE / "k65_gen3_preregistration.py"


def _load_preregistration_contract() -> Any:
    module_name = "leopard2_k65_gen3_preregistration_for_plan_runner"
    expected = PREREGISTRATION_PATH.resolve(strict=True)
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        if Path(str(getattr(loaded, "__file__", ""))).resolve() != expected:
            raise RuntimeError("generation-3 preregistration came from another path")
        return loaded
    specification = importlib.util.spec_from_file_location(module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load generation-3 preregistration contract")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    if Path(str(getattr(module, "__file__", ""))).resolve() != expected:
        raise RuntimeError("generation-3 preregistration resolved to another path")
    return module


prereg = _load_preregistration_contract()
contract = prereg.contract

GENERATION = 3
PLAN_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-campaign-plan/v2"
STATE_HISTORY_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-state-history/v2"
STATE_EVENT_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-state-event/v2"
BUDGET_LEDGER_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-budget-ledger/v2"

ROUNDS = 25
ITERATIONS = 31
WARMUP = 64
REUSE = 8192
EXPECTED_CHILD_PROCESS_COUNT = 1650
MATRIX_SHA256 = prereg.V2_MATRIX_SHA256
SELECTOR_ARGUMENT = "--k65r65-b64-packed-terminal-mode"
CHILD_TIMEOUT_SETUP_SECONDS = 120
CHILD_TIMEOUT_LOGICAL_BYTES_PER_SECOND = 1 << 30

TARGET_ORDERS = (
    ("main", "candidate", "control", "control", "candidate", "main"),
    ("control", "main", "candidate", "candidate", "main", "control"),
    ("candidate", "control", "main", "main", "control", "candidate"),
)
NEIGHBOR_ORDERS = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
)

STATES = (
    "INIT", "PREREGISTERED", "QUALIFYING", "QUALIFIED", "BRIDGING",
    "BRIDGED", "ARMING", "PRESAMPLING", "ARMED", "TIMING", "FINALIZING",
    "SETUP_INVALID", "ENV_REJECTED", "REJECTED", "ACCEPTED",
)
TERMINAL_STATES = frozenset((
    "SETUP_INVALID", "ENV_REJECTED", "REJECTED", "ACCEPTED",
))
ALLOWED_TRANSITIONS = {
    "INIT": frozenset(("PREREGISTERED",)),
    "PREREGISTERED": frozenset(("QUALIFYING",)),
    "QUALIFYING": frozenset(("ENV_REJECTED", "QUALIFIED")),
    "QUALIFIED": frozenset(("BRIDGING",)),
    "BRIDGING": frozenset(("ENV_REJECTED", "BRIDGED")),
    "BRIDGED": frozenset(("ARMING",)),
    "ARMING": frozenset(("SETUP_INVALID", "PRESAMPLING")),
    "PRESAMPLING": frozenset(("ENV_REJECTED", "ARMED")),
    "ARMED": frozenset(("TIMING",)),
    "TIMING": frozenset(("REJECTED", "FINALIZING")),
    "FINALIZING": frozenset(("REJECTED", "ACCEPTED")),
    "SETUP_INVALID": frozenset(),
    "ENV_REJECTED": frozenset(),
    "REJECTED": frozenset(),
    "ACCEPTED": frozenset(),
}

PAIR_KEYS = frozenset(("benchmark_cpu", "reserved_sibling"))
STATE_EVENT_KEYS = frozenset((
    "schema", "index", "from_state", "to_state", "selected_pair",
    "budget_committed", "evidence_attempt_committed", "terminal",
))
STATE_HISTORY_KEYS = frozenset((
    "schema", "generation", "lane_class", "lane_index", "events",
    "highest_state", "terminal_state", "budget_committed",
    "evidence_attempt_committed", "frozen_pair",
))
BUDGET_LEDGER_KEYS = frozenset((
    "schema", "generation", "setup_invalid_used", "setup_invalid_limit",
    "environment_rejected_used", "environment_rejected_limit",
    "evidence_attempts_used", "evidence_attempts_limit", "frozen_pair",
    "histories",
))
CHILD_PLAN_KEYS = frozenset((
    "index", "cell_index", "cell_id", "round", "slot", "implementation",
    "argv_tail", "timeout_budget",
))
TIMEOUT_KEYS = frozenset((
    "timeout_seconds", "setup_seconds", "logical_bytes_per_second_floor",
    "measured_metric_count", "calls_per_metric", "logical_bytes_per_call",
    "logical_byte_visits",
))
PLAN_KEYS = frozenset((
    "schema", "generation", "mode", "safe_to_execute",
    "candidate_timing_performed", "preregistration_schema",
    "preregistration_sha256", "preregistration_ratified",
    "campaign", "state_graph",
    "cells", "child_plans", "child_process_count",
))


class PlanError(ValueError):
    """A generation-3 plan, state history, or budget ledger is invalid."""


def _fail(message: str) -> NoReturn:
    raise PlanError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _exact_object(value: Any, keys: frozenset[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == keys,
             f"{label} is not an exact object")
    return value


def _bounded_int(value: Any, minimum: int, maximum: int, label: str) -> int:
    _require(type(value) is int and minimum <= value <= maximum,
             f"{label} is outside its exact integer bounds")
    return value


def _pair(value: Any, label: str) -> dict[str, int]:
    record = _exact_object(value, PAIR_KEYS, label)
    benchmark = _bounded_int(
        record["benchmark_cpu"], 0, contract.MAX_CPU_ID,
        f"{label} benchmark CPU")
    sibling = _bounded_int(
        record["reserved_sibling"], 0, contract.MAX_CPU_ID,
        f"{label} reserved sibling")
    _require(benchmark != sibling, f"{label} repeats one logical CPU")
    return {"benchmark_cpu": benchmark, "reserved_sibling": sibling}


def _absolute_path(value: Any, label: str) -> str:
    _require(type(value) is str and value != "" and "\x00" not in value,
             f"{label} is not a path")
    path = PurePosixPath(value)
    _require(path.is_absolute() and ".." not in path.parts and "." not in path.parts,
             f"{label} is not a canonical absolute POSIX path")
    return str(path)


def cells() -> list[dict[str, Any]]:
    values = (
        ("target-k65-r65-b64-q1", 65, 65, 64, "target", True, 1.05),
        ("k-control-k64-r65-b64-q1", 64, 65, 64, "neighbor", False, None),
        ("k-control-k66-r65-b64-q1", 66, 65, 64, "neighbor", False, None),
        ("r-control-k65-r64-b64-q1", 65, 64, 64, "neighbor", False, None),
        ("r-control-k65-r66-b64-q1", 65, 66, 64, "neighbor", False, None),
        ("byte-control-k65-r65-b63-q1", 65, 65, 63, "neighbor", False, None),
        ("byte-control-k65-r65-b65-q1", 65, 65, 65, "neighbor", False, None),
        ("retained-k65-r65-b8192-q1", 65, 65, 8192, "neighbor", True, 0.98),
        ("balanced-layout-k64-r64-b64-q1", 64, 64, 64, "neighbor", True, 0.98),
        ("balanced-layout-k128-r128-b64-q1", 128, 128, 64, "neighbor", True, 0.98),
        ("layout-context-k79-r32-b64-q1", 79, 32, 64, "neighbor", True, None),
        ("layout-context-k62-r8-b64-q1", 62, 8, 64, "neighbor", True, None),
        ("layout-context-k66-r16-b64-q1", 66, 16, 64, "neighbor", True, None),
    )
    result = []
    for index, (name, k, r, size, role, compare_main, floor) in enumerate(values):
        result.append({
            "id": name,
            "K": k,
            "R": r,
            "bytes": size,
            "loss": 1,
            "batch": 1,
            "reuse": REUSE,
            "role": role,
            "compare_main": compare_main,
            "main_floor": floor,
            "measure_one_shot": True,
            "seed": 0x6565B640 + index,
        })
    # Generation 2's frozen hash predates the shared canonical helper and omits
    # its trailing LF.  Preserve that exact historical byte projection here.
    legacy_bytes = contract.canonical_json_bytes(result)[:-1]
    _require(hashlib.sha256(legacy_bytes).hexdigest() == MATRIX_SHA256,
             "generation-3 matrix differs from frozen generation 2")
    return result


def state_graph_record() -> dict[str, list[str]]:
    return {state: sorted(ALLOWED_TRANSITIONS[state]) for state in STATES}


def _budget_for_state(state: str, *, evidence_committed: bool) -> str | None:
    if evidence_committed:
        return "evidence_attempt"
    if state == "SETUP_INVALID":
        return "setup_invalid"
    if state == "ENV_REJECTED":
        return "environment_rejected"
    return None


def _abandoned_budget_for_state(state: str) -> str:
    if state in {
        "QUALIFYING", "QUALIFIED", "BRIDGING", "BRIDGED", "PRESAMPLING",
    }:
        return "environment_rejected"
    _require(state in {"INIT", "PREREGISTERED", "ARMING"},
             "abandoned pre-timing state cannot be classified")
    return "setup_invalid"


def attempt_state_history_record(
    *, lane_class: str, lane_index: int, states: Sequence[str],
    selected_pair: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Construct one exact one-way history and classify its consumed budget."""
    _require(lane_class in {"setup", "environment", "evidence"},
             "state-history lane class differs")
    index = _bounded_int(lane_index, 1, 1000, "state-history lane index")
    _require(type(states) in (list, tuple) and len(states) >= 1 and
             all(type(state) is str and state in STATES for state in states),
             "state sequence is invalid")
    sequence = list(states)
    _require(sequence[0] == "INIT", "state sequence does not start at INIT")
    _require(len(sequence) == len(set(sequence)),
             "state sequence revisits a state")
    pair = _pair(selected_pair, "state-history selected pair") \
        if selected_pair is not None else None
    events: list[dict[str, Any]] = []
    terminal = sequence[-1] if sequence[-1] in TERMINAL_STATES else None
    final_evidence_committed = "ARMED" in sequence
    final_budget = _budget_for_state(
        sequence[-1], evidence_committed=final_evidence_committed)
    if final_budget is None:
        final_budget = _abandoned_budget_for_state(sequence[-1])
    evidence_committed = False
    armed_reached = False
    for event_index, state in enumerate(sequence):
        previous = sequence[event_index - 1] if event_index else None
        if previous is not None:
            _require(state in ALLOWED_TRANSITIONS[previous],
                     f"illegal state transition {previous}->{state}")
        if state == "ARMED":
            armed_reached = True
            evidence_committed = True
            _require(pair is not None, "ARMED state does not bind a pair")
        if state == "TIMING":
            _require(armed_reached, "TIMING was reached without ARMED")
        event_pair = copy.deepcopy(pair) if armed_reached else None
        budget = "evidence_attempt" if evidence_committed else None
        if event_index == len(sequence) - 1:
            budget = final_budget
        events.append({
            "schema": STATE_EVENT_SCHEMA,
            "index": event_index,
            "from_state": previous,
            "to_state": state,
            "selected_pair": event_pair,
            "budget_committed": budget,
            "evidence_attempt_committed": evidence_committed,
            "terminal": state in TERMINAL_STATES,
        })
    budget = final_budget
    expected_lane_class = {
        "setup_invalid": "setup",
        "environment_rejected": "environment",
        "evidence_attempt": "evidence",
        None: lane_class,
    }[budget]
    _require(lane_class == expected_lane_class,
             "state history is stored in the wrong budget lane")
    _require(not (terminal is not None and sequence[-1] not in TERMINAL_STATES),
             "state-history terminal differs")
    return {
        "schema": STATE_HISTORY_SCHEMA,
        "generation": GENERATION,
        "lane_class": lane_class,
        "lane_index": index,
        "events": events,
        "highest_state": sequence[-1],
        "terminal_state": terminal,
        "budget_committed": budget,
        "evidence_attempt_committed": evidence_committed,
        "frozen_pair": copy.deepcopy(pair) if armed_reached else None,
    }


def validate_attempt_state_history(value: Any) -> dict[str, Any]:
    record = _exact_object(value, STATE_HISTORY_KEYS, "state history")
    _require(type(record["events"]) is list and record["events"],
             "state history has no events")
    states: list[str] = []
    pair = None
    for index, event in enumerate(record["events"]):
        item = _exact_object(event, STATE_EVENT_KEYS, f"state event {index}")
        _require(item["schema"] == STATE_EVENT_SCHEMA and item["index"] == index,
                 f"state event {index} metadata differs")
        states.append(item["to_state"])
        if item["selected_pair"] is not None:
            observed_pair = _pair(item["selected_pair"], f"state event {index} pair")
            _require(pair is None or pair == observed_pair,
                     "state history changes its frozen pair")
            pair = observed_pair
    rebuilt = attempt_state_history_record(
        lane_class=record["lane_class"], lane_index=record["lane_index"],
        states=states, selected_pair=pair)
    _require(contract.exact_json_equal(record, rebuilt),
             "state history differs from its fixed point")
    return copy.deepcopy(rebuilt)


def budget_ledger_record(
    preregistration: Mapping[str, Any], histories: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    registration = prereg.validate_preregistration(
        preregistration, verify_files=False)
    _require(type(histories) in (list, tuple), "histories are not a sequence")
    validated = [validate_attempt_state_history(item) for item in histories]
    for lane_class in ("setup", "environment", "evidence"):
        indexes = [item["lane_index"] for item in validated
                   if item["lane_class"] == lane_class]
        _require(indexes == list(range(1, len(indexes) + 1)),
                 f"{lane_class} lane indexes are not a gapless prefix")
    pairs = [item["frozen_pair"] for item in validated
             if item["frozen_pair"] is not None]
    _require(not pairs or all(pair == pairs[0] for pair in pairs),
             "generation-3 frozen pair changed across attempts")
    setup_used = sum(item["budget_committed"] == "setup_invalid"
                     for item in validated)
    environment_used = sum(item["budget_committed"] == "environment_rejected"
                           for item in validated)
    evidence_used = sum(item["evidence_attempt_committed"] for item in validated)
    budgets = registration["budgets"]
    _require(setup_used <= budgets["setup_invalid"],
             "setup-invalid budget is exhausted")
    _require(environment_used <= budgets["environment_rejected"],
             "environment-rejected budget is exhausted")
    _require(evidence_used <= budgets["evidence_attempts"],
             "evidence-attempt budget is exhausted")
    return {
        "schema": BUDGET_LEDGER_SCHEMA,
        "generation": GENERATION,
        "setup_invalid_used": setup_used,
        "setup_invalid_limit": budgets["setup_invalid"],
        "environment_rejected_used": environment_used,
        "environment_rejected_limit": budgets["environment_rejected"],
        "evidence_attempts_used": evidence_used,
        "evidence_attempts_limit": budgets["evidence_attempts"],
        "frozen_pair": copy.deepcopy(pairs[0]) if pairs else None,
        "histories": validated,
    }


def validate_budget_ledger(
    value: Any, preregistration: Mapping[str, Any],
) -> dict[str, Any]:
    record = _exact_object(value, BUDGET_LEDGER_KEYS, "budget ledger")
    rebuilt = budget_ledger_record(preregistration, record["histories"])
    _require(contract.exact_json_equal(record, rebuilt),
             "budget ledger differs from its fixed point")
    return copy.deepcopy(rebuilt)


def child_timeout_budget(
    implementation: str, cell: Mapping[str, Any],
) -> dict[str, int]:
    _require(implementation in {"candidate", "control", "main"},
             "child implementation differs")
    metric_count = 2
    if implementation != "main" and cell["measure_one_shot"] is True:
        metric_count += 1
    calls_per_metric = ITERATIONS * int(cell["reuse"]) + WARMUP
    logical_bytes_per_call = (
        (int(cell["K"]) + int(cell["R"])) * int(cell["bytes"]) *
        int(cell["batch"])
    )
    visits = metric_count * calls_per_metric * logical_bytes_per_call
    timeout = CHILD_TIMEOUT_SETUP_SECONDS + (
        visits + CHILD_TIMEOUT_LOGICAL_BYTES_PER_SECOND - 1
    ) // CHILD_TIMEOUT_LOGICAL_BYTES_PER_SECOND
    return {
        "timeout_seconds": timeout,
        "setup_seconds": CHILD_TIMEOUT_SETUP_SECONDS,
        "logical_bytes_per_second_floor":
            CHILD_TIMEOUT_LOGICAL_BYTES_PER_SECOND,
        "measured_metric_count": metric_count,
        "calls_per_metric": calls_per_metric,
        "logical_bytes_per_call": logical_bytes_per_call,
        "logical_byte_visits": visits,
    }


def benchmark_argv_tail(
    implementation: str, cell: Mapping[str, Any],
) -> list[str]:
    """Return only arguments after the executable; paths/CPUs are post-ARMED."""
    _require(implementation in {"candidate", "control", "main"},
             "benchmark implementation differs")
    arguments = [
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", str(cell["loss"]),
        "--batch", str(cell["batch"]), "--reuse", str(cell["reuse"]),
        "--iterations", str(ITERATIONS), "--warmup", str(WARMUP),
        "--threads", "1", "--seed", str(cell["seed"]),
    ]
    if implementation != "main":
        arguments.extend((
            "--profile", "high", "--field", "gf8", "--backend", "avx2",
            "--skip-legacy", "--retain-samples", "--attest-source",
            "--measure-one-shot-encode", SELECTOR_ARGUMENT,
            "1" if implementation == "candidate" else "0",
        ))
    arguments.extend(("--json", "-"))
    return arguments


def child_plans() -> list[dict[str, Any]]:
    """Return the immutable logical schedule without launch authority."""
    result: list[dict[str, Any]] = []
    for cell_index, cell in enumerate(cells()):
        orders = TARGET_ORDERS if cell["compare_main"] else NEIGHBOR_ORDERS
        for round_index in range(ROUNDS):
            order = orders[round_index % len(orders)]
            for slot, implementation in enumerate(order):
                result.append({
                    "index": len(result),
                    "cell_index": cell_index,
                    "cell_id": cell["id"],
                    "round": round_index,
                    "slot": slot,
                    "implementation": implementation,
                    "argv_tail": benchmark_argv_tail(implementation, cell),
                    "timeout_budget": child_timeout_budget(implementation, cell),
                })
    _require(len(result) == EXPECTED_CHILD_PROCESS_COUNT,
             "generation-3 child process count differs")
    return result


def _validate_preregistration_for_plan(value: Any) -> tuple[dict[str, Any], bool]:
    if type(value) is not dict:
        _fail("plan preregistration is not an object")
    if value.get("schema") == prereg.TEMPLATE_SCHEMA:
        return prereg.validate_preregistration_template(
            value, verify_files=False), False
    if value.get("schema") == prereg.PREREGISTRATION_SCHEMA:
        return prereg.validate_preregistration(value, verify_files=False), True
    _fail("plan preregistration schema differs")


def campaign_plan_record(
    *, preregistration: Mapping[str, Any],
) -> dict[str, Any]:
    registration, ratified = _validate_preregistration_for_plan(preregistration)
    campaign = registration["campaign"]
    frozen_campaign_fields = {
        "generation": GENERATION,
        "matrix_sha256": MATRIX_SHA256,
        "cell_count": len(cells()),
        "rounds_per_cell": ROUNDS,
        "iterations_per_child": ITERATIONS,
        "warmup_per_child": WARMUP,
        "reuse_per_child": REUSE,
        "expected_child_process_count": EXPECTED_CHILD_PROCESS_COUNT,
    }
    _require(
        all(campaign[key] == expected
            for key, expected in frozen_campaign_fields.items()),
        "preregistration campaign differs from the frozen plan constants")
    plans = child_plans()
    _require(
        campaign["expected_child_process_count"] == len(plans),
        "preregistration child-process count differs from the generated plan")
    return {
        "schema": PLAN_SCHEMA,
        "generation": GENERATION,
        "mode": "plan-only",
        "safe_to_execute": False,
        "candidate_timing_performed": False,
        "preregistration_schema": registration["schema"],
        "preregistration_sha256": contract.canonical_sha256(registration),
        "preregistration_ratified": ratified,
        "campaign": copy.deepcopy(campaign),
        "state_graph": state_graph_record(),
        "cells": cells(),
        "child_plans": plans,
        "child_process_count": len(plans),
    }


def validate_campaign_plan(value: Any, preregistration: Mapping[str, Any]) -> dict[str, Any]:
    record = _exact_object(value, PLAN_KEYS, "campaign plan")
    _require(type(record["child_plans"]) is list and
             len(record["child_plans"]) == EXPECTED_CHILD_PROCESS_COUNT,
             "campaign plan child list differs")
    for index, child in enumerate(record["child_plans"]):
        item = _exact_object(child, CHILD_PLAN_KEYS, f"child plan {index}")
        _exact_object(item["timeout_budget"], TIMEOUT_KEYS,
                      f"child plan {index} timeout")
        _require(item["index"] == index, f"child plan {index} index differs")
    rebuilt = campaign_plan_record(preregistration=preregistration)
    _require(contract.exact_json_equal(record, rebuilt),
             "campaign plan differs from its fixed point")
    return copy.deepcopy(rebuilt)


def load_campaign_plan(data: bytes, preregistration: Mapping[str, Any]) -> dict[str, Any]:
    return validate_campaign_plan(
        contract.strict_json_loads(data, "generation-3 campaign plan"),
        preregistration)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-only", action="store_true", required=True)
    parser.add_argument("--preregistration", type=Path, required=True)
    options = parser.parse_args()
    try:
        data = options.preregistration.read_bytes()
        value = contract.strict_json_loads(data, "plan preregistration")
        if type(value) is dict and value.get("schema") == prereg.TEMPLATE_SCHEMA:
            registration = prereg.load_preregistration_template(
                data, verify_files=False)
        elif (type(value) is dict and
              value.get("schema") == prereg.PREREGISTRATION_SCHEMA):
            registration = prereg.load_preregistration(
                data, verify_files=True)
        else:
            _fail("plan preregistration schema differs")
        plan = campaign_plan_record(preregistration=registration)
        sys.stdout.buffer.write(contract.canonical_json_bytes(plan))
    except (OSError, PlanError, prereg.PreregistrationError,
            contract.QualificationError) as error:
        print(f"K65 generation-3 plan rejected: {error}", file=sys.stderr)
        return 1
    return 0


__all__ = (
    "ALLOWED_TRANSITIONS", "BUDGET_LEDGER_SCHEMA", "EXPECTED_CHILD_PROCESS_COUNT",
    "GENERATION", "PLAN_SCHEMA", "PlanError", "STATE_HISTORY_SCHEMA",
    "attempt_state_history_record", "benchmark_argv_tail", "budget_ledger_record",
    "campaign_plan_record", "cells", "child_plans", "child_timeout_budget",
    "contract", "load_campaign_plan", "state_graph_record",
    "validate_attempt_state_history", "validate_budget_ledger",
    "validate_campaign_plan",
)


if __name__ == "__main__":
    raise SystemExit(main())
