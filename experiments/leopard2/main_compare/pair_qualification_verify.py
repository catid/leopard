#!/usr/bin/env python3
"""Independent offline verifier for pair qualification acquisition and bridge.

Only evidence files explicitly named on the command line are read.  The verifier
does not inspect procfs/sysfs, acquire clocks, launch children, mutate affinity, or
perform candidate timing.  It independently binds canonical wire bytes to the
preregistered policy and geometry and then replays the pure fixed points.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
from pathlib import Path
import sys
from typing import Any, NoReturn, Sequence


def _load_bridge_contract() -> Any:
    module_name = "leopard2_pair_qualification_bridge_for_verifier"
    expected = Path(__file__).resolve().with_name(
        "pair_qualification_bridge_contract.py")
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        if Path(getattr(loaded, "__file__", "")).resolve() != expected:
            raise RuntimeError("pair qualification bridge came from another path")
        return loaded
    specification = importlib.util.spec_from_file_location(module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load the pair qualification bridge contract")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    if Path(getattr(module, "__file__", "")).resolve() != expected:
        raise RuntimeError("pair qualification bridge resolved to another path")
    return module


bridge = _load_bridge_contract()
contract = bridge.contract

VERDICT_SCHEMA = "leopard2-pair-qualification-independent-verdict/v1"
VERIFICATION_METHOD = "canonical-wire-independent-fixed-point-replay/v1"
SELF_TEST_OUTPUT = \
    "leopard2 pair qualification independent verifier self-test passed"
PAIR_KEYS = frozenset(("benchmark_cpu", "reserved_sibling"))
VERDICT_KEYS = frozenset((
    "schema", "verification_method", "policy_sha256", "acquisition_sha256",
    "bridge_sha256", "selected_pair", "bridge_accepted",
    "fixed_point_verified", "verifier_host_access_performed",
    "candidate_timing_performed", "shared_host_claim_ceiling",
))


class VerificationError(bridge.BridgeError):
    """Independent replay did not establish the requested verdict."""


def _fail(message: str) -> NoReturn:
    raise VerificationError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _hex64(value: Any, label: str) -> str:
    _require(type(value) is str and len(value) == 64 and
             all(character in "0123456789abcdef" for character in value),
             f"{label} is not exact lowercase SHA-256")
    return value


def _pair(value: Any, label: str) -> dict[str, int] | None:
    if value is None:
        return None
    _require(type(value) is dict and set(value) == PAIR_KEYS,
             f"{label} is not an exact pair object")
    benchmark = value["benchmark_cpu"]
    sibling = value["reserved_sibling"]
    _require(type(benchmark) is int and type(sibling) is int and
             0 <= benchmark <= contract.MAX_CPU_ID and
             0 <= sibling <= contract.MAX_CPU_ID and benchmark != sibling,
             f"{label} is outside the CPU-pair bound")
    return {"benchmark_cpu": benchmark, "reserved_sibling": sibling}


def _independent_cross_check(
    acquisition: dict[str, Any], bridge_record: dict[str, Any], *,
    policy: dict[str, Any], frozen_pair: dict[str, int] | None,
) -> None:
    scan = acquisition["scan"]
    selected = contract.selected_pair_from_scan(
        scan, expected_policy=policy, expected_frozen_pair=frozen_pair)
    _require(contract.exact_json_equal(bridge_record["selected_pair"], selected),
             "independent selected-pair replay differs")
    _require(bridge_record["policy_sha256"] == contract.canonical_sha256(policy),
             "independent policy hash differs")
    _require(bridge_record["scan_sha256"] == contract.canonical_sha256(scan),
             "independent scan hash differs")
    _require(bridge_record["acquisition_sha256"] ==
             contract.canonical_sha256(acquisition),
             "independent acquisition hash differs")
    scan_tail = scan["windows"][-1]["after"]
    _require(contract.exact_json_equal(
        bridge_record["windows"][0]["before"], scan_tail),
        "independent scan-to-bridge handoff differs")
    _require(contract.exact_json_equal(
        bridge_record["windows"][-1]["after"],
        bridge_record["campaign_presample_before"]),
        "independent bridge-to-presample handoff differs")
    nonidle = [
        cpu for cpu in bridge_record["guarded_cpus"]
        if bridge_record["cpu_aggregates"][str(cpu)]["nonidle_jiffies"] != 0 or
        bridge_record["cpu_aggregates"][str(cpu)]["nonidle_window_count"] != 0
    ]
    not_live = [
        cpu for cpu in bridge_record["guarded_cpus"]
        if not bridge_record["cpu_aggregates"][str(cpu)]["overall_live"] or
        bridge_record["cpu_aggregates"][str(cpu)]["not_live_window_count"] != 0
    ]
    _require(nonidle == bridge_record["nonidle_guarded_cpus"] and
             not_live == bridge_record["not_live_guarded_cpus"] and
             bridge_record["bridge_accepted"] == (not nonidle and not not_live),
             "independent bridge acceptance replay differs")
    expected_ceiling = {
        "promotion_eligible": False,
        "host_exclusivity_proved": False,
        "whole_campaign_interval_observed": False,
        "causal_performance_claim_allowed": False,
    }
    _require(contract.exact_json_equal(
        bridge_record["shared_host_claim_ceiling"], expected_ceiling),
        "independent shared-host claim ceiling differs")
    _require(bridge_record["host_mutation_performed"] is False and
             bridge_record["candidate_timing_performed"] is False,
             "independent bridge replay found mutation or candidate timing")


def verify_pair_qualification_bundle(
    acquisition_data: bytes, bridge_data: bytes, *, expected_policy: Any,
    expected_policy_sha256: Any, expected_frozen_pair: Any,
    expected_acquisition_window_count: int,
    expected_acquisition_nominal_window_ns: int,
    expected_bridge_geometry: Any,
) -> dict[str, Any]:
    policy = contract.validate_qualification_policy(expected_policy)
    policy_hash = _hex64(expected_policy_sha256, "expected policy hash")
    _require(policy_hash == contract.canonical_sha256(policy),
             "preregistered policy hash differs from policy bytes")
    frozen = _pair(expected_frozen_pair, "expected frozen pair")
    acquisition_raw = contract.strict_json_loads(
        acquisition_data, "independent acquisition")
    acquisition = bridge.validate_acquisition_for_bridge(
        acquisition_raw, expected_policy=policy,
        expected_frozen_pair=frozen,
        expected_window_count=expected_acquisition_window_count,
        expected_nominal_window_ns=expected_acquisition_nominal_window_ns)
    _require(acquisition_data == contract.canonical_json_bytes(acquisition),
             "independent acquisition wire bytes are not canonical")
    bridge_raw = contract.strict_json_loads(
        bridge_data, "independent qualification bridge")
    bridge_record = bridge.validate_pair_qualification_bridge(
        bridge_raw, acquisition, expected_policy=policy,
        expected_frozen_pair=frozen,
        expected_acquisition_window_count=expected_acquisition_window_count,
        expected_acquisition_nominal_window_ns=
        expected_acquisition_nominal_window_ns,
        expected_bridge_geometry=expected_bridge_geometry)
    _require(bridge_data == contract.canonical_json_bytes(bridge_record),
             "independent bridge wire bytes are not canonical")
    _independent_cross_check(
        acquisition, bridge_record, policy=policy, frozen_pair=frozen)
    verdict = {
        "schema": VERDICT_SCHEMA,
        "verification_method": VERIFICATION_METHOD,
        "policy_sha256": policy_hash,
        "acquisition_sha256": contract.canonical_sha256(acquisition),
        "bridge_sha256": contract.canonical_sha256(bridge_record),
        "selected_pair": copy.deepcopy(bridge_record["selected_pair"]),
        "bridge_accepted": bridge_record["bridge_accepted"],
        "fixed_point_verified": True,
        "verifier_host_access_performed": False,
        "candidate_timing_performed": False,
        "shared_host_claim_ceiling": copy.deepcopy(
            bridge_record["shared_host_claim_ceiling"]),
    }
    _require(type(verdict) is dict and set(verdict) == VERDICT_KEYS,
             "independent verdict key set differs")
    return verdict


def require_accepted_pair_qualification_bundle(
    acquisition_data: bytes, bridge_data: bytes, **keywords: Any,
) -> dict[str, Any]:
    verdict = verify_pair_qualification_bundle(
        acquisition_data, bridge_data, **keywords)
    _require(verdict["bridge_accepted"] is True,
             "pair qualification bridge was valid but not accepted")
    return verdict


def _snapshot(started: int, idle: int) -> dict[str, Any]:
    return contract.shared_snapshot_record(
        read_started_monotonic_ns=started,
        read_finished_monotonic_ns=started + 1_000_000,
        counters={
            2: {field: idle if field == "idle" else 0
                for field in contract.COUNTER_FIELDS},
            6: {field: idle if field == "idle" else 0
                for field in contract.COUNTER_FIELDS},
        })


def _self_test_fixture() -> tuple[bytes, bytes, dict[str, Any], dict[str, int]]:
    policy = contract.qualification_policy_record(
        clock_ticks_per_second=100, candidate_primary_cpus=[2],
        excluded_pairs=[], domain_mode="pair-only")
    topology = contract.topology_record([
        {
            "cpu": cpu, "online": True, "physical_package_id": 0,
            "die_id": 0, "core_id": 0, "thread_siblings": [2, 6],
            "domain_cpus": [2, 6],
        }
        for cpu in (2, 6)
    ])
    qualification_snapshots = [
        _snapshot(1_000_000_000, 1_000),
        _snapshot(1_251_000_000, 1_025),
        _snapshot(1_502_000_000, 1_050),
    ]
    scan = contract.pair_qualification_scan_record(
        policy, allowed_cpu_set_at_launch=[2, 6],
        topology_before=topology, topology_after=topology,
        snapshots=qualification_snapshots,
        frozen_pair_from_prior_attempt=None)
    acquisition = {
        "schema": bridge.ACQUISITION_SCHEMA,
        "acquisition_method": bridge.ACQUISITION_METHOD,
        "sources": dict(bridge._SOURCE_ITEMS),
        "policy": copy.deepcopy(policy),
        "policy_sha256": contract.canonical_sha256(policy),
        "requested_window_count": 2,
        "nominal_window_ns": 250_000_000,
        "frozen_pair_from_prior_attempt": None,
        "allowed_cpu_set_at_launch": [2, 6],
        "allowed_cpu_set_after_scan": [2, 6],
        "clock_ticks_per_second_at_launch": 100,
        "clock_ticks_per_second_after_scan": 100,
        "topology_before_sha256": contract.canonical_sha256(topology),
        "topology_after_sha256": contract.canonical_sha256(topology),
        "scan": scan,
        "scan_sha256": contract.canonical_sha256(scan),
        "host_mutation_performed": False,
        "candidate_timing_performed": False,
    }
    bridge_snapshots = [
        qualification_snapshots[-1],
        _snapshot(1_753_000_000, 1_075),
        _snapshot(2_004_000_000, 1_100),
    ]
    geometry = bridge.bridge_geometry_record(
        minimum_window_count=2, maximum_window_count=3,
        nominal_window_ns=250_000_000,
        maximum_handoff_elapsed_ns=2_000_000_000)
    bridge_record = bridge.pair_qualification_bridge_record(
        acquisition, expected_policy=policy, expected_frozen_pair=None,
        expected_acquisition_window_count=2,
        expected_acquisition_nominal_window_ns=250_000_000,
        expected_bridge_geometry=geometry, snapshots=bridge_snapshots)
    return (
        contract.canonical_json_bytes(acquisition),
        contract.canonical_json_bytes(bridge_record), policy, geometry,
    )


def self_test() -> None:
    acquisition_data, bridge_data, policy, geometry = _self_test_fixture()
    keywords = {
        "expected_policy": policy,
        "expected_policy_sha256": contract.canonical_sha256(policy),
        "expected_frozen_pair": None,
        "expected_acquisition_window_count": 2,
        "expected_acquisition_nominal_window_ns": 250_000_000,
        "expected_bridge_geometry": geometry,
    }
    verdict = require_accepted_pair_qualification_bundle(
        acquisition_data, bridge_data, **keywords)
    _require(verdict["selected_pair"] == {
        "benchmark_cpu": 2, "reserved_sibling": 6,
    }, "self-test selected pair differs")
    mutated = contract.strict_json_loads(bridge_data, "self-test bridge")
    mutated["bridge_accepted"] = False
    try:
        verify_pair_qualification_bundle(
            acquisition_data, contract.canonical_json_bytes(mutated), **keywords)
    except contract.QualificationError:
        pass
    else:
        _fail("self-test accepted a producer-only bridge attestation")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="independently verify a pair qualification bridge")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("self-test")
    verify = subparsers.add_parser("verify")
    verify.add_argument("--policy", required=True)
    verify.add_argument("--acquisition", required=True)
    verify.add_argument("--bridge", required=True)
    verify.add_argument("--expected-policy-sha256", required=True)
    verify.add_argument("--expected-acquisition-window-count", required=True, type=int)
    verify.add_argument("--expected-acquisition-nominal-window-ns", required=True,
                        type=int)
    verify.add_argument("--minimum-bridge-window-count", required=True, type=int)
    verify.add_argument("--maximum-bridge-window-count", required=True, type=int)
    verify.add_argument("--bridge-nominal-window-ns", required=True, type=int)
    verify.add_argument("--maximum-bridge-handoff-elapsed-ns", required=True,
                        type=int)
    verify.add_argument("--frozen-benchmark-cpu", type=int)
    verify.add_argument("--frozen-reserved-sibling", type=int)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    options = _parser().parse_args(argv)
    try:
        if options.command == "self-test":
            self_test()
            print(SELF_TEST_OUTPUT)
            return 0
        policy_data = Path(options.policy).read_bytes()
        policy = contract.strict_json_loads(policy_data, "preregistered policy")
        _require(policy_data == contract.canonical_json_bytes(policy),
                 "preregistered policy wire bytes are not canonical")
        neither_frozen = (
            options.frozen_benchmark_cpu is None and
            options.frozen_reserved_sibling is None
        )
        both_frozen = (
            options.frozen_benchmark_cpu is not None and
            options.frozen_reserved_sibling is not None
        )
        _require(neither_frozen or both_frozen,
                 "frozen pair CLI requires both CPU roles or neither")
        frozen = (
            {
                "benchmark_cpu": options.frozen_benchmark_cpu,
                "reserved_sibling": options.frozen_reserved_sibling,
            }
            if both_frozen else None
        )
        geometry = bridge.bridge_geometry_record(
            minimum_window_count=options.minimum_bridge_window_count,
            maximum_window_count=options.maximum_bridge_window_count,
            nominal_window_ns=options.bridge_nominal_window_ns,
            maximum_handoff_elapsed_ns=
            options.maximum_bridge_handoff_elapsed_ns)
        verdict = require_accepted_pair_qualification_bundle(
            Path(options.acquisition).read_bytes(),
            Path(options.bridge).read_bytes(),
            expected_policy=policy,
            expected_policy_sha256=options.expected_policy_sha256,
            expected_frozen_pair=frozen,
            expected_acquisition_window_count=
            options.expected_acquisition_window_count,
            expected_acquisition_nominal_window_ns=
            options.expected_acquisition_nominal_window_ns,
            expected_bridge_geometry=geometry)
        sys.stdout.buffer.write(contract.canonical_json_bytes(verdict))
        return 0
    except (OSError, contract.QualificationError) as error:
        print(f"pair qualification verification failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = (
    "SELF_TEST_OUTPUT", "VERDICT_SCHEMA", "VerificationError",
    "require_accepted_pair_qualification_bundle", "self_test",
    "verify_pair_qualification_bundle",
)
