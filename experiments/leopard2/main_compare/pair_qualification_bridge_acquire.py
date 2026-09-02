#!/usr/bin/env python3
"""Read-only live acquisition for the prospectively frozen v19 bridge.

The caller supplies a validated qualification acquisition and its exact final
shared snapshot.  This module keeps that snapshot as the bridge head, acquires
exactly two further shared ``/proc/stat`` endpoints through an injected reader,
and returns the existing pure bridge record only after independent fixed-point
verification accepts it.  It has no CLI, benchmark import, process-launch
surface, host mutator, or evidence writer.
"""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import sys
from typing import Any, Callable, NoReturn, Sequence


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


acquire = _load_local_module(
    "leopard2_pair_qualification_acquire_for_v19_bridge",
    "pair_qualification_acquire.py")
bridge = _load_local_module(
    "leopard2_pair_qualification_bridge_contract_for_v19_acquisition",
    "pair_qualification_bridge_contract.py")
verifier = _load_local_module(
    "leopard2_pair_qualification_verifier_for_v19_bridge",
    "pair_qualification_verify.py")
contract = bridge.contract

V19_BRIDGE_WINDOW_COUNT = 2
V19_BRIDGE_NOMINAL_WINDOW_NS = 1_000_000_000
V19_MAXIMUM_HANDOFF_ELAPSED_NS = 5_000_000_000


class BridgeAcquisitionError(bridge.BridgeError):
    """Live host state or evidence violated the frozen v19 bridge boundary."""


def _fail(message: str) -> NoReturn:
    raise BridgeAcquisitionError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _reader_call(label: str, function: Callable[..., Any], *args: Any) -> Any:
    """Translate normal injected-reader failures without swallowing interrupts."""
    try:
        return function(*args)
    except BridgeAcquisitionError:
        raise
    except Exception as error:
        raise BridgeAcquisitionError(f"{label} failed: {error}") from error


def _canonical_cpu_list(value: Any, label: str) -> list[int]:
    _require(type(value) in (list, tuple) and
             0 < len(value) <= contract.MAX_CPU_COUNT,
             f"{label} is not a bounded CPU sequence")
    cpus: list[int] = []
    for cpu in value:
        _require(type(cpu) is int and 0 <= cpu <= contract.MAX_CPU_ID,
                 f"{label} contains an invalid CPU")
        cpus.append(cpu)
    _require(cpus == sorted(set(cpus)), f"{label} is not canonical")
    return cpus


def _policy_sha256(value: Any, policy: dict[str, Any]) -> str:
    _require(type(value) is str and len(value) == 64 and
             all(character in "0123456789abcdef" for character in value),
             "expected policy hash is not exact lowercase SHA-256")
    expected = contract.canonical_sha256(policy)
    _require(value == expected,
             "expected policy hash differs from the qualification policy")
    return expected


def v19_bridge_geometry_record() -> dict[str, int]:
    """Return a fresh copy of the prospectively frozen bridge geometry."""
    return bridge.bridge_geometry_record(
        minimum_window_count=V19_BRIDGE_WINDOW_COUNT,
        maximum_window_count=V19_BRIDGE_WINDOW_COUNT,
        nominal_window_ns=V19_BRIDGE_NOMINAL_WINDOW_NS,
        maximum_handoff_elapsed_ns=V19_MAXIMUM_HANDOFF_ELAPSED_NS)


def validate_v19_bridge_geometry(value: Any) -> dict[str, int]:
    geometry = bridge.validate_bridge_geometry(value)
    expected = v19_bridge_geometry_record()
    _require(contract.exact_json_equal(geometry, expected),
             "live bridge geometry differs from the frozen v19 geometry")
    return copy.deepcopy(expected)


def _check_live_host_projection(
    reader: acquire.HostReader, *, policy: dict[str, Any],
    expected_allowed_cpus: Sequence[int], expected_topology: dict[str, Any],
    phase: str,
) -> None:
    allowed = _canonical_cpu_list(
        _reader_call(f"{phase} affinity read", reader.allowed_cpus),
        f"{phase} allowed CPU set")
    _require(allowed == list(expected_allowed_cpus),
             f"allowed CPU set changed at {phase}")
    ticks = _reader_call(
        f"{phase} clock-tick read", reader.clock_ticks_per_second)
    _require(type(ticks) is int and
             ticks == policy["clock_ticks_per_second"],
             f"clock tick rate changed at {phase}")
    topology = _reader_call(
        f"{phase} topology acquisition", acquire._acquire_topology,
        reader, allowed, policy["candidate_primary_cpus"])
    _require(contract.exact_json_equal(topology, expected_topology),
             f"qualification topology changed at {phase}")


def _capture_snapshot(
    reader: acquire.HostReader, observed_cpus: Sequence[int], phase: str,
) -> dict[str, Any]:
    return _reader_call(
        f"{phase} shared snapshot acquisition", acquire._capture_snapshot,
        reader, observed_cpus)


def acquire_v19_pair_qualification_bridge(
    reader: acquire.HostReader, acquisition_value: Any,
    scan_tail_value: Any, *, expected_policy: Any,
    expected_policy_sha256: Any, expected_frozen_pair: Any,
    expected_acquisition_window_count: int,
    expected_acquisition_nominal_window_ns: int,
    expected_bridge_geometry: Any,
) -> dict[str, Any]:
    """Acquire, independently verify, and return one accepted v19 bridge.

    ``scan_tail_value`` is deliberately explicit: the live caller must hand the
    exact in-memory endpoint produced by qualification across this boundary.
    Normal reader failures are translated to :class:`BridgeAcquisitionError`;
    ``KeyboardInterrupt`` and other ``BaseException`` subclasses propagate.
    """
    policy = contract.validate_qualification_policy(expected_policy)
    _require(policy["domain_mode"] == "pair-only",
             "live v19 bridge requires the pair-only policy")
    policy_hash = _policy_sha256(expected_policy_sha256, policy)
    geometry = validate_v19_bridge_geometry(expected_bridge_geometry)
    acquisition_record = bridge.validate_acquisition_for_bridge(
        acquisition_value, expected_policy=policy,
        expected_frozen_pair=expected_frozen_pair,
        expected_window_count=expected_acquisition_window_count,
        expected_nominal_window_ns=expected_acquisition_nominal_window_ns)
    scan = acquisition_record["scan"]
    scan_tail = scan["windows"][-1]["after"]
    _require(contract.exact_json_equal(scan_tail_value, scan_tail),
             "supplied bridge head differs from the exact qualification scan tail")
    allowed_cpus = acquisition_record["allowed_cpu_set_after_scan"]
    expected_topology = scan["topology_after"]
    observed_cpus = sorted(int(cpu) for cpu in scan_tail["cpus"])
    _check_live_host_projection(
        reader, policy=policy, expected_allowed_cpus=allowed_cpus,
        expected_topology=expected_topology, phase="bridge start")

    bridge_started = scan_tail["read_finished_monotonic_ns"]
    _require(bridge_started <= contract.MAX_MONOTONIC_NS -
             geometry["maximum_handoff_elapsed_ns"],
             "v19 bridge deadline overflows the monotonic clock bound")
    deadline = bridge_started + geometry["maximum_handoff_elapsed_ns"]
    now = _reader_call(
        "bridge-start monotonic read", reader.monotonic_ns)
    _require(type(now) is int and
             bridge_started <= now <= deadline,
             "live reader is outside the qualification monotonic handoff")
    _require(now <= deadline -
             V19_BRIDGE_WINDOW_COUNT * V19_BRIDGE_NOMINAL_WINDOW_NS,
             "v19 bridge cannot fit its fixed windows before the deadline")

    snapshots = [copy.deepcopy(scan_tail_value)]
    for index in range(V19_BRIDGE_WINDOW_COUNT):
        remaining_windows = V19_BRIDGE_WINDOW_COUNT - index
        now = _reader_call(
            f"bridge window {index} monotonic read", reader.monotonic_ns)
        _require(type(now) is int and
                 bridge_started <= now <= deadline,
                 f"bridge window {index} starts outside the handoff deadline")
        _require(now <= deadline -
                 remaining_windows * V19_BRIDGE_NOMINAL_WINDOW_NS,
                 f"bridge window {index} cannot fit before the deadline")
        _reader_call(
            f"bridge window {index} sleep", reader.sleep_ns,
            V19_BRIDGE_NOMINAL_WINDOW_NS)
        snapshots.append(_capture_snapshot(
            reader, observed_cpus, f"bridge window {index}"))

    _check_live_host_projection(
        reader, policy=policy, expected_allowed_cpus=allowed_cpus,
        expected_topology=expected_topology, phase="bridge finish")
    try:
        record = bridge.pair_qualification_bridge_record(
            acquisition_record, expected_policy=policy,
            expected_frozen_pair=expected_frozen_pair,
            expected_acquisition_window_count=
            expected_acquisition_window_count,
            expected_acquisition_nominal_window_ns=
            expected_acquisition_nominal_window_ns,
            expected_bridge_geometry=geometry, snapshots=snapshots)
        verdict = verifier.require_accepted_pair_qualification_bundle(
            contract.canonical_json_bytes(acquisition_record),
            contract.canonical_json_bytes(record), expected_policy=policy,
            expected_policy_sha256=policy_hash,
            expected_frozen_pair=expected_frozen_pair,
            expected_acquisition_window_count=
            expected_acquisition_window_count,
            expected_acquisition_nominal_window_ns=
            expected_acquisition_nominal_window_ns,
            expected_bridge_geometry=geometry)
    except Exception as error:
        raise BridgeAcquisitionError(
            f"v19 bridge fixed-point verification failed: {error}") from error
    _require(verdict["bridge_sha256"] == contract.canonical_sha256(record),
             "independent verifier returned a different bridge binding")
    return copy.deepcopy(record)


__all__ = (
    "BridgeAcquisitionError", "V19_BRIDGE_NOMINAL_WINDOW_NS",
    "V19_BRIDGE_WINDOW_COUNT", "V19_MAXIMUM_HANDOFF_ELAPSED_NS",
    "acquire_v19_pair_qualification_bridge", "validate_v19_bridge_geometry",
    "v19_bridge_geometry_record",
)
