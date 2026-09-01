#!/usr/bin/env python3
"""Pure generation-3 K65 execution-authority records.

This module is deliberately incapable of discovering a host, opening an
executable, acquiring qualification evidence, publishing a lane, or launching a
child.  It validates already-retained records and canonical wire bytes.  The
acquisition runner remains responsible for obtaining stable executable handles
and for durably publishing the returned ARMED record before any child launch.

The final ARMED fixed point replays the qualification acquisition, bridge, and
independent verdict.  A caller cannot turn opaque hashes into launch authority.
"""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import re
import sys
from typing import Any, Mapping, NoReturn


HERE = Path(__file__).resolve().parent
MAIN_COMPARE = HERE.parent / "main_compare"
PREREGISTRATION_PATH = HERE / "k65_gen3_preregistration.py"
PLAN_PATH = HERE / "run_k65r65_b64_packed_terminal_gen3_abba.py"
PAIR_VERIFIER_PATH = MAIN_COMPARE / "pair_qualification_verify.py"


def _load_path_module(name: str, path: Path) -> Any:
    expected = path.resolve(strict=True)
    loaded = sys.modules.get(name)
    if loaded is not None:
        if Path(str(getattr(loaded, "__file__", ""))).resolve() != expected:
            raise RuntimeError(f"{name} came from another path")
        return loaded
    specification = importlib.util.spec_from_file_location(name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {expected}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    if Path(str(getattr(module, "__file__", ""))).resolve() != expected:
        raise RuntimeError(f"{name} resolved to another path")
    return module


prereg = _load_path_module(
    "leopard2_k65_gen3_preregistration_for_execution_contract",
    PREREGISTRATION_PATH,
)
pair_verifier = _load_path_module(
    "leopard2_pair_qualification_verifier_for_k65_execution_contract",
    PAIR_VERIFIER_PATH,
)
plan_contract = _load_path_module(
    "leopard2_k65_gen3_plan_for_execution_contract", PLAN_PATH)
contract = prereg.contract
bridge_contract = pair_verifier.bridge


GENERATION = 3
PREREGISTRATION_SCHEMA_EXPECTATION = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-preregistration/v4"
HOST_AUTHORITY_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-host-authority/v1"
HOST_INSTANCE_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-host-instance/v1"
ARTIFACT_BUNDLE_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-artifact-bundle/v1"
CANDIDATE_SOURCE_AUTHORITY_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-candidate-source-authority/v1"
QUALIFICATION_BINDING_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-qualification-binding/v1"
ARMED_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-armed/v2"

EXACT_MAIN_AUTHORITY_RECORD_SHA256 = \
    "ce3d6e647fd098f558ee01855523c0643b8c175578315d02a168a87992e5ae01"
EXACT_MAIN_VERIFIER_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-baseline-verification/v1"
EXACT_MAIN_PATH_VARIANT_RAW_SHA256 = \
    "0baae845bbf30d2b3b213c02501a31c2d15fd125965d898f7a100fa6d0ede46d"
EXACT_MAIN_PATH_VARIANT_SIZE = 1_175_456
EXACT_MAIN_NORMALIZED_IDENTITY_SCHEMA = \
    "leopard2-gf8-exact-main-normalized-code-identity/v1"
EXACT_MAIN_NORMALIZED_COMBINED_SHA256 = \
    "ddfef166af6c1dafd989019f87694526693623fa4ea2aa9e4d74f97c012fa093"
EXACT_MAIN_SOURCE_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
EXACT_MAIN_SOURCE_TREE = "b7c8830d96a978f6ec14fe747095f066e351ae72"
EXACT_MAIN_ADAPTER_COMMIT = "6d4f690fe5ba7cf08f2f8f80a263765266b462e0"

CANDIDATE_LAUNCH_PROTOCOL = "leopard2-k65-mode-1-json-stdin/v1"
CONTROL_LAUNCH_PROTOCOL = "leopard2-k65-mode-0-json-stdin/v1"
EXACT_MAIN_LAUNCH_PROTOCOL = "leopard1-exact-main-json-stdin/v1"
PUBLICATION_RULE = "atomic-no-replace-before-child/v1"

HOST_AUTHORITY_KEYS = frozenset((
    "schema", "machine_id_sha256", "hostname", "architecture", "cpu_model",
))
HOST_INSTANCE_KEYS = frozenset((
    "schema", "authority", "boot_id", "kernel_release", "online_cpus",
    "allowed_cpus", "clock_ticks_per_second", "topology_sha256",
))
SOURCE_KEYS = frozenset(("commit", "tree"))
ROLE_KEYS = frozenset((
    "role", "raw_sha256", "size", "launch_protocol", "handle_id",
    "handle_device", "handle_inode",
))
CANDIDATE_SOURCE_AUTHORITY_KEYS = frozenset((
    "schema", "commit", "tree", "build_provenance_sha256",
    "reproducible_build_core_sha256", "authority_record_sha256",
    "authority_ledger_sha256", "verification_method",
))
ROLES_KEYS = frozenset(("candidate", "control", "main"))
EXACT_MAIN_AUTHORITY_KEYS = frozenset((
    "record_sha256", "verifier_schema", "verifier_verdict_sha256", "role",
    "raw_sha256", "normalized_identity_schema",
    "normalized_combined_sha256", "source_commit", "source_tree",
    "adapter_commit",
))
ARTIFACT_BUNDLE_KEYS = frozenset((
    "schema", "preregistration_sha256", "candidate_source",
    "candidate_source_authority", "controller_bindings_sha256", "roles",
    "candidate_control_same_handle", "exact_main_authority",
))
PAIR_KEYS = frozenset(("benchmark_cpu", "reserved_sibling"))
QUALIFICATION_BINDING_KEYS = frozenset((
    "schema", "preregistration_sha256", "host_instance_sha256",
    "selected_track", "selected_pair", "policy_sha256",
    "track_selection_sha256", "acquisition_sha256", "bridge_sha256",
    "independent_verdict_sha256", "campaign_presample_before_sha256",
    "bridge_accepted", "candidate_timing_performed",
))
QUALIFICATION_EVIDENCE_KEYS = frozenset((
    "policy_a_scan", "policy_b_scan", "track_selection",
    "expected_frozen_pair", "acquisition_data", "bridge_data",
    "independent_verdict_data",
))
ARMED_KEYS = frozenset((
    "schema", "generation", "state", "evidence_attempt",
    "evidence_attempt_limit", "prior_armed_chain_sha256",
    "preregistration_sha256", "plan_sha256", "host_instance_sha256",
    "artifact_bundle_sha256", "qualification_binding_sha256",
    "authority_bundle_sha256", "attempt_manifest_sha256",
    "lane_binding_sha256",
    "controller_bindings_sha256", "candidate_source", "selected_pair",
    "armed_monotonic_ns", "bridge_deadline_monotonic_ns",
    "publication_rule", "crash_after_publication_consumes",
))


class ExecutionContractError(prereg.PreregistrationError):
    """A purported generation-3 execution-authority record is invalid."""


def _fail(message: str) -> NoReturn:
    raise ExecutionContractError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _exact_object(value: Any, keys: frozenset[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == keys,
             f"{label} is not an exact object")
    return value


def _sha256(value: Any, label: str) -> str:
    _require(type(value) is str and
             re.fullmatch(r"[0-9a-f]{64}", value) is not None,
             f"{label} is not a lowercase SHA-256")
    return value


def _git_oid(value: Any, label: str) -> str:
    _require(type(value) is str and
             re.fullmatch(r"[0-9a-f]{40}", value) is not None,
             f"{label} is not a lowercase 40-hex object ID")
    return value


def _bounded_int(value: Any, minimum: int, maximum: int, label: str) -> int:
    _require(type(value) is int and minimum <= value <= maximum,
             f"{label} is outside its exact integer bounds")
    return value


def _bounded_ascii(value: Any, maximum: int, label: str) -> str:
    _require(type(value) is str and value.strip() == value and
             0 < len(value) <= maximum and
             all(0x20 <= ord(character) <= 0x7e for character in value),
             f"{label} is not bounded canonical ASCII")
    return value


def _canonical_cpus(value: Any, label: str) -> list[int]:
    _require(type(value) is list and 0 < len(value) <= contract.MAX_CPU_COUNT,
             f"{label} is not a bounded non-empty list")
    result = [
        _bounded_int(cpu, 0, contract.MAX_CPU_ID, f"{label} entry")
        for cpu in value
    ]
    _require(result == sorted(set(result)), f"{label} is not canonical")
    return result


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


def _validated_preregistration(
    value: Any, *, host_authority: Mapping[str, Any],
) -> dict[str, Any]:
    registration = prereg.validate_preregistration(value, verify_files=False)
    supplied_host = validate_host_authority(host_authority)
    registration_host = {
        key: supplied_host[key]
        for key in ("machine_id_sha256", "hostname", "architecture", "cpu_model")
    }
    _require(registration["schema"] == PREREGISTRATION_SCHEMA_EXPECTATION and
             registration["safe_to_execute"] is False and
             registration["candidate_executable"]["mode"] == "frozen-sha256" and
             contract.exact_json_equal(
                 registration["host_authority"], registration_host),
             "execution authority requires the frozen non-executable v2 preregistration")
    return registration


_require(
    prereg.PREREGISTRATION_SCHEMA == PREREGISTRATION_SCHEMA_EXPECTATION and
    prereg.EXACT_MAIN_AUTHORITY_RECORD_SHA256 ==
        EXACT_MAIN_AUTHORITY_RECORD_SHA256 and
    prereg.EXACT_MAIN_SOURCE_COMMIT == EXACT_MAIN_SOURCE_COMMIT and
    prereg.EXACT_MAIN_SOURCE_TREE == EXACT_MAIN_SOURCE_TREE and
    prereg.EXACT_MAIN_ADAPTER_COMMIT == EXACT_MAIN_ADAPTER_COMMIT and
    pair_verifier.VERDICT_SCHEMA ==
        "leopard2-pair-qualification-independent-verdict/v1",
    "generation-3 execution dependencies moved from their reviewed fixed points",
)


def host_authority_record(
    *, machine_id_sha256: str, hostname: str, architecture: str,
    cpu_model: str,
) -> dict[str, Any]:
    machine = _sha256(machine_id_sha256, "host machine-id hash")
    host = _bounded_ascii(hostname, 253, "host name")
    _require(re.fullmatch(
        r"(?=.{1,253}\Z)[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?"
        r"(?:\.[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?)*",
        host) is not None, "host name is not canonical lower-case DNS syntax")
    arch = _bounded_ascii(architecture, 64, "host architecture")
    _require(arch == "x86_64", "host architecture is not the frozen x86_64 target")
    model = _bounded_ascii(cpu_model, 512, "host CPU model")
    return {
        "schema": HOST_AUTHORITY_SCHEMA,
        "machine_id_sha256": machine,
        "hostname": host,
        "architecture": arch,
        "cpu_model": model,
    }


def validate_host_authority(value: Any) -> dict[str, Any]:
    record = _exact_object(value, HOST_AUTHORITY_KEYS, "host authority")
    expected = host_authority_record(
        machine_id_sha256=record["machine_id_sha256"],
        hostname=record["hostname"], architecture=record["architecture"],
        cpu_model=record["cpu_model"])
    _require(contract.exact_json_equal(record, expected),
             "host authority differs from its fixed point")
    return copy.deepcopy(expected)


def host_instance_record(
    *, authority: Mapping[str, Any], boot_id: str, kernel_release: str,
    online_cpus: list[int], allowed_cpus: list[int],
    clock_ticks_per_second: int, topology_sha256: str,
) -> dict[str, Any]:
    host = validate_host_authority(authority)
    boot = _bounded_ascii(boot_id, 36, "host boot ID")
    _require(re.fullmatch(
        r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
        boot) is not None, "host boot ID is not a canonical UUID")
    kernel = _bounded_ascii(kernel_release, 256, "kernel release")
    online = _canonical_cpus(online_cpus, "online CPU list")
    allowed = _canonical_cpus(allowed_cpus, "allowed CPU list")
    _require(set(allowed) <= set(online),
             "allowed CPU list is not a subset of online CPUs")
    ticks = _bounded_int(
        clock_ticks_per_second, 1, contract.MAX_CLOCK_TICKS_PER_SECOND,
        "host clock tick rate")
    topology = _sha256(topology_sha256, "host topology hash")
    return {
        "schema": HOST_INSTANCE_SCHEMA,
        "authority": host,
        "boot_id": boot,
        "kernel_release": kernel,
        "online_cpus": online,
        "allowed_cpus": allowed,
        "clock_ticks_per_second": ticks,
        "topology_sha256": topology,
    }


def validate_host_instance(value: Any) -> dict[str, Any]:
    record = _exact_object(value, HOST_INSTANCE_KEYS, "host instance")
    expected = host_instance_record(
        authority=record["authority"], boot_id=record["boot_id"],
        kernel_release=record["kernel_release"],
        online_cpus=record["online_cpus"], allowed_cpus=record["allowed_cpus"],
        clock_ticks_per_second=record["clock_ticks_per_second"],
        topology_sha256=record["topology_sha256"])
    _require(contract.exact_json_equal(record, expected),
             "host instance differs from its fixed point")
    return copy.deepcopy(expected)


def exact_main_authority_record(
    *, verifier_verdict_sha256: str,
) -> dict[str, Any]:
    return {
        "record_sha256": EXACT_MAIN_AUTHORITY_RECORD_SHA256,
        "verifier_schema": EXACT_MAIN_VERIFIER_SCHEMA,
        "verifier_verdict_sha256": _sha256(
            verifier_verdict_sha256, "exact-main verifier verdict hash"),
        "role": "path_variant",
        "raw_sha256": EXACT_MAIN_PATH_VARIANT_RAW_SHA256,
        "normalized_identity_schema": EXACT_MAIN_NORMALIZED_IDENTITY_SCHEMA,
        "normalized_combined_sha256": EXACT_MAIN_NORMALIZED_COMBINED_SHA256,
        "source_commit": EXACT_MAIN_SOURCE_COMMIT,
        "source_tree": EXACT_MAIN_SOURCE_TREE,
        "adapter_commit": EXACT_MAIN_ADAPTER_COMMIT,
    }


def validate_exact_main_authority(value: Any) -> dict[str, Any]:
    record = _exact_object(
        value, EXACT_MAIN_AUTHORITY_KEYS, "exact-main path-variant authority")
    expected = exact_main_authority_record(
        verifier_verdict_sha256=record["verifier_verdict_sha256"])
    _require(contract.exact_json_equal(record, expected),
             "exact-main path-variant authority differs")
    return copy.deepcopy(expected)


def _artifact_role(value: Any, expected_role: str) -> dict[str, Any]:
    record = _exact_object(value, ROLE_KEYS, f"{expected_role} artifact role")
    _require(record["role"] == expected_role,
             f"{expected_role} artifact label differs")
    raw_hash = _sha256(record["raw_sha256"], f"{expected_role} raw hash")
    size = _bounded_int(
        record["size"], 1, (1 << 63) - 1, f"{expected_role} artifact size")
    handle = _sha256(record["handle_id"], f"{expected_role} handle ID")
    device = _bounded_int(
        record["handle_device"], 0, (1 << 64) - 1,
        f"{expected_role} handle device")
    inode = _bounded_int(
        record["handle_inode"], 1, (1 << 64) - 1,
        f"{expected_role} handle inode")
    protocols = {
        "candidate": CANDIDATE_LAUNCH_PROTOCOL,
        "control": CONTROL_LAUNCH_PROTOCOL,
        "main": EXACT_MAIN_LAUNCH_PROTOCOL,
    }
    _require(record["launch_protocol"] == protocols[expected_role],
             f"{expected_role} launch protocol differs")
    return {
        "role": expected_role,
        "raw_sha256": raw_hash,
        "size": size,
        "launch_protocol": protocols[expected_role],
        "handle_id": handle,
        "handle_device": device,
        "handle_inode": inode,
    }


def candidate_source_authority_record(
    preregistration: Mapping[str, Any], *,
    build_provenance_sha256: str,
    reproducible_build_core_sha256: str,
    authority_record_sha256: str,
    authority_ledger_sha256: str,
    host_authority: Mapping[str, Any],
) -> dict[str, Any]:
    """Bind the live-validated clean build closure to the ratified source."""
    registration = _validated_preregistration(
        preregistration, host_authority=host_authority)
    source = registration["candidate_source"]
    candidate = registration["candidate_executable"]
    provenance_sha = _sha256(
        build_provenance_sha256, "candidate build-provenance hash")
    reproducible_sha = _sha256(
        reproducible_build_core_sha256,
        "candidate reproducible-build core hash")
    authority_record_sha = _sha256(
        authority_record_sha256, "candidate authority-record file hash")
    authority_ledger_sha = _sha256(
        authority_ledger_sha256, "candidate authority-ledger hash")
    _require(
        provenance_sha == candidate["build_provenance_sha256"] and
        reproducible_sha == candidate["reproducible_build_core_sha256"] and
        authority_record_sha == candidate["authority_record_sha256"] and
        authority_ledger_sha == candidate["authority_ledger_sha256"],
        "candidate build closure differs from preregistration")
    return {
        "schema": CANDIDATE_SOURCE_AUTHORITY_SCHEMA,
        "commit": _git_oid(source["commit"], "candidate source commit"),
        "tree": _git_oid(source["tree"], "candidate source tree"),
        "build_provenance_sha256": provenance_sha,
        "reproducible_build_core_sha256": reproducible_sha,
        "authority_record_sha256": authority_record_sha,
        "authority_ledger_sha256": authority_ledger_sha,
        "verification_method":
            "tracked-git-tree-object-compile-archive-link-clean-replay/v1",
    }


def validate_candidate_source_authority(
    value: Any, preregistration: Mapping[str, Any], *,
    host_authority: Mapping[str, Any],
) -> dict[str, Any]:
    record = _exact_object(
        value, CANDIDATE_SOURCE_AUTHORITY_KEYS,
        "candidate source authority")
    expected = candidate_source_authority_record(
        preregistration,
        build_provenance_sha256=record["build_provenance_sha256"],
        reproducible_build_core_sha256=
            record["reproducible_build_core_sha256"],
        authority_record_sha256=record["authority_record_sha256"],
        authority_ledger_sha256=record["authority_ledger_sha256"],
        host_authority=host_authority)
    _require(contract.exact_json_equal(record, expected),
             "candidate source authority differs from its fixed point")
    return copy.deepcopy(expected)


def artifact_bundle_record(
    preregistration: Mapping[str, Any], *, roles: Mapping[str, Any],
    exact_main_authority: Mapping[str, Any],
    candidate_source_authority: Mapping[str, Any],
    host_authority: Mapping[str, Any],
) -> dict[str, Any]:
    registration = _validated_preregistration(
        preregistration, host_authority=host_authority)
    role_object = _exact_object(roles, ROLES_KEYS, "artifact roles")
    normalized_roles = {
        role: _artifact_role(role_object[role], role)
        for role in ("candidate", "control", "main")
    }
    candidate = normalized_roles["candidate"]
    control = normalized_roles["control"]
    main = normalized_roles["main"]
    _require(
        candidate["raw_sha256"] == control["raw_sha256"] ==
            registration["candidate_executable"]["sha256"] and
        candidate["size"] == control["size"] ==
            registration["candidate_executable"]["size"] and
        candidate["handle_id"] == control["handle_id"] and
        candidate["handle_device"] == control["handle_device"] and
        candidate["handle_inode"] == control["handle_inode"],
        "candidate and control are not two protocols over one frozen handle",
    )
    _require(main["raw_sha256"] == EXACT_MAIN_PATH_VARIANT_RAW_SHA256 and
             main["size"] == EXACT_MAIN_PATH_VARIANT_SIZE and
             (main["handle_device"], main["handle_inode"]) !=
                 (candidate["handle_device"], candidate["handle_inode"]) and
             main["handle_id"] != candidate["handle_id"],
             "main artifact is not the distinct frozen path-variant handle")
    authority = validate_exact_main_authority(exact_main_authority)
    _require(main["raw_sha256"] == authority["raw_sha256"],
             "main role and exact-main authority raw hashes differ")
    for role in normalized_roles.values():
        prereg.reject_denylisted_evidence_hash(
            role["raw_sha256"], registration)
    source = _exact_object(
        registration["candidate_source"], SOURCE_KEYS, "candidate source")
    candidate_source = {
        "commit": _git_oid(source["commit"], "candidate source commit"),
        "tree": _git_oid(source["tree"], "candidate source tree"),
    }
    source_authority = validate_candidate_source_authority(
        candidate_source_authority, registration,
        host_authority=host_authority)
    return {
        "schema": ARTIFACT_BUNDLE_SCHEMA,
        "preregistration_sha256": contract.canonical_sha256(registration),
        "candidate_source": candidate_source,
        "candidate_source_authority": source_authority,
        "controller_bindings_sha256": contract.canonical_sha256(
            registration["controller_bindings"]),
        "roles": normalized_roles,
        "candidate_control_same_handle": True,
        "exact_main_authority": authority,
    }


def validate_artifact_bundle(
    value: Any, preregistration: Mapping[str, Any], *,
    host_authority: Mapping[str, Any],
) -> dict[str, Any]:
    record = _exact_object(value, ARTIFACT_BUNDLE_KEYS, "artifact bundle")
    expected = artifact_bundle_record(
        preregistration, roles=record["roles"],
        exact_main_authority=record["exact_main_authority"],
        candidate_source_authority=record["candidate_source_authority"],
        host_authority=host_authority)
    _require(contract.exact_json_equal(record, expected),
             "artifact bundle differs from its fixed point")
    return copy.deepcopy(expected)


def _qualification_evidence(value: Any) -> dict[str, Any]:
    evidence = _exact_object(
        value, QUALIFICATION_EVIDENCE_KEYS, "qualification evidence inputs")
    for field in ("acquisition_data", "bridge_data", "independent_verdict_data"):
        _require(type(evidence[field]) is bytes,
                 f"qualification {field} is not canonical wire bytes")
    return evidence


def qualification_binding_record(
    preregistration: Mapping[str, Any], host_instance: Mapping[str, Any], *,
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    """Replay and bind one accepted, untimed qualification-to-presample chain."""
    host = validate_host_instance(host_instance)
    registration = _validated_preregistration(
        preregistration, host_authority=host["authority"])
    inputs = _qualification_evidence(evidence)
    selection = prereg.validate_qualification_track_selection(
        inputs["track_selection"], registration,
        policy_a_scan=inputs["policy_a_scan"],
        policy_b_scan=inputs["policy_b_scan"],
        expected_frozen_pair=inputs["expected_frozen_pair"])
    selected_pair = selection["selected_pair"]
    _require(selection["selected_track"] in {"pair-and-domain", "pair-only"} and
             selected_pair is not None and
             selection["selection_status"] in {
                 "selected-lowest-primary", "selected-frozen-pair"} and
             selection["candidate_timing_performed"] is False,
             "qualification track selection did not select one untimed pair")
    pair = _pair(selected_pair, "qualification selected pair")
    policy_index = 0 if selection["selected_track"] == "pair-and-domain" else 1
    policies = registration["qualification"]["policies"]
    _require(policy_index < len(policies), "selected qualification track is absent")
    policy = contract.validate_qualification_policy(policies[policy_index])
    policy_sha = contract.canonical_sha256(policy)
    _require(registration["qualification"]["policy_sha256s"][policy_index] ==
             policy_sha and
             host["clock_ticks_per_second"] == policy["clock_ticks_per_second"],
             "host tick rate or selected policy hash differs")
    geometry = registration["qualification"]["geometry"]
    bridge_geometry = {
        "minimum_window_count": geometry["bridge_minimum_window_count"],
        "maximum_window_count": geometry["bridge_maximum_window_count"],
        "nominal_window_ns": geometry["bridge_nominal_window_ns"],
        "maximum_handoff_elapsed_ns": geometry["maximum_handoff_elapsed_ns"],
    }
    expected_verdict = pair_verifier.require_accepted_pair_qualification_bundle(
        inputs["acquisition_data"], inputs["bridge_data"],
        expected_policy=policy, expected_policy_sha256=policy_sha,
        expected_frozen_pair=inputs["expected_frozen_pair"],
        expected_acquisition_window_count=geometry["scan_window_count"],
        expected_acquisition_nominal_window_ns=geometry["scan_nominal_window_ns"],
        expected_bridge_geometry=bridge_geometry)
    verdict = contract.strict_json_loads(
        inputs["independent_verdict_data"], "independent qualification verdict")
    _require(inputs["independent_verdict_data"] ==
             contract.canonical_json_bytes(verdict),
             "independent qualification verdict bytes are not canonical")
    _require(contract.exact_json_equal(verdict, expected_verdict),
             "independent qualification verdict differs from replay")
    acquisition = contract.strict_json_loads(
        inputs["acquisition_data"], "qualification acquisition")
    bridge = contract.strict_json_loads(
        inputs["bridge_data"], "qualification bridge")
    selected_scan = (inputs["policy_a_scan"] if policy_index == 0
                     else inputs["policy_b_scan"])
    _require(contract.exact_json_equal(acquisition["scan"], selected_scan),
             "selected track scan and acquisition scan differ")
    topology = contract.validate_topology(acquisition["scan"]["topology_before"])
    topology_cpus = {item["cpu"] for item in topology["cpus"]}
    _require(
        acquisition["topology_before_sha256"] == host["topology_sha256"] ==
            acquisition["topology_after_sha256"] and
        acquisition["allowed_cpu_set_at_launch"] == host["allowed_cpus"] ==
            acquisition["allowed_cpu_set_after_scan"] and
        topology_cpus <= set(host["online_cpus"]) and
        set(pair.values()) <= set(host["allowed_cpus"]),
        "qualification acquisition differs from the frozen host instance",
    )
    _require(
        contract.exact_json_equal(bridge["selected_pair"], pair) and
        contract.exact_json_equal(verdict["selected_pair"], pair) and
        bridge["policy_sha256"] == verdict["policy_sha256"] == policy_sha and
        bridge["acquisition_sha256"] == verdict["acquisition_sha256"] and
        verdict["bridge_sha256"] == contract.canonical_sha256(bridge) and
        bridge["bridge_accepted"] is True and
        verdict["bridge_accepted"] is True and
        bridge["candidate_timing_performed"] is False and
        verdict["candidate_timing_performed"] is False,
        "qualification pair, hashes, acceptance, or timing flags differ",
    )
    return {
        "schema": QUALIFICATION_BINDING_SCHEMA,
        "preregistration_sha256": contract.canonical_sha256(registration),
        "host_instance_sha256": contract.canonical_sha256(host),
        "selected_track": selection["selected_track"],
        "selected_pair": pair,
        "policy_sha256": policy_sha,
        "track_selection_sha256": contract.canonical_sha256(selection),
        "acquisition_sha256": contract.canonical_sha256(acquisition),
        "bridge_sha256": contract.canonical_sha256(bridge),
        "independent_verdict_sha256": contract.canonical_sha256(verdict),
        "campaign_presample_before_sha256":
            bridge["campaign_presample_before_sha256"],
        "bridge_accepted": True,
        "candidate_timing_performed": False,
    }


def validate_qualification_binding(
    value: Any, preregistration: Mapping[str, Any],
    host_instance: Mapping[str, Any], *, evidence: Mapping[str, Any],
) -> dict[str, Any]:
    record = _exact_object(
        value, QUALIFICATION_BINDING_KEYS, "qualification binding")
    expected = qualification_binding_record(
        preregistration, host_instance, evidence=evidence)
    _require(contract.exact_json_equal(record, expected),
             "qualification binding differs from its replayed fixed point")
    return copy.deepcopy(expected)


def validate_armed_record_shape(value: Any) -> dict[str, Any]:
    """Validate ledger recovery syntax without conferring launch authority.

    Historical ARMED ledger entries may be used only to recover their gapless
    attempt numbers and canonical hash chain.  Creation or use of a new entry
    must go through :func:`validate_armed_record`, which fully replays every
    qualification and artifact sidecar.
    """
    record = _exact_object(value, ARMED_KEYS, "ARMED record shape")
    _require(record["schema"] == ARMED_SCHEMA and
             record["generation"] == GENERATION and
             record["state"] == "ARMED",
             "ARMED record shape metadata differs")
    limit = _bounded_int(
        record["evidence_attempt_limit"], 1, 1000,
        "ARMED evidence-attempt limit")
    _require(limit == 3, "ARMED evidence-attempt limit differs")
    attempt = _bounded_int(
        record["evidence_attempt"], 1, limit, "ARMED evidence attempt")
    prior_value = record["prior_armed_chain_sha256"]
    if attempt == 1:
        _require(prior_value is None,
                 "first ARMED record unexpectedly has a prior chain")
        prior = None
    else:
        prior = _sha256(prior_value, "prior ARMED chain hash")
    hashes = {
        name: _sha256(record[name], f"ARMED {name}")
        for name in (
            "preregistration_sha256", "plan_sha256", "host_instance_sha256",
            "artifact_bundle_sha256", "qualification_binding_sha256",
            "authority_bundle_sha256", "attempt_manifest_sha256",
            "lane_binding_sha256",
            "controller_bindings_sha256",
        )
    }
    source = _exact_object(
        record["candidate_source"], SOURCE_KEYS, "ARMED candidate source")
    candidate_source = {
        "commit": _git_oid(source["commit"], "ARMED candidate source commit"),
        "tree": _git_oid(source["tree"], "ARMED candidate source tree"),
    }
    pair = _pair(record["selected_pair"], "ARMED selected pair")
    armed_ns = _bounded_int(
        record["armed_monotonic_ns"], 1, (1 << 63) - 1,
        "ARMED monotonic timestamp")
    deadline_ns = _bounded_int(
        record["bridge_deadline_monotonic_ns"], 1, (1 << 63) - 1,
        "ARMED bridge deadline")
    _require(armed_ns <= deadline_ns,
             "ARMED timestamp exceeds its bridge deadline")
    _require(record["publication_rule"] == PUBLICATION_RULE and
             record["crash_after_publication_consumes"] is True,
             "ARMED publication or crash-consumption rule differs")
    return {
        "schema": ARMED_SCHEMA,
        "generation": GENERATION,
        "state": "ARMED",
        "evidence_attempt": attempt,
        "evidence_attempt_limit": limit,
        "prior_armed_chain_sha256": prior,
        **hashes,
        "candidate_source": candidate_source,
        "selected_pair": pair,
        "armed_monotonic_ns": armed_ns,
        "bridge_deadline_monotonic_ns": deadline_ns,
        "publication_rule": PUBLICATION_RULE,
        "crash_after_publication_consumes": True,
    }


def armed_record(
    preregistration: Mapping[str, Any], plan: Mapping[str, Any],
    host_instance: Mapping[str, Any], artifact_bundle: Mapping[str, Any],
    qualification_binding: Mapping[str, Any], *,
    qualification_evidence: Mapping[str, Any], evidence_attempt: int,
    prior_armed_chain_sha256: str | None,
    authority_bundle_sha256: str, attempt_manifest_sha256: str,
    lane_binding_sha256: str,
    armed_monotonic_ns: int,
) -> dict[str, Any]:
    """Recompute every authority input and return the record to publish durably."""
    host = validate_host_instance(host_instance)
    registration = _validated_preregistration(
        preregistration, host_authority=host["authority"])
    logical_plan = plan_contract.validate_campaign_plan(plan, registration)
    artifacts = validate_artifact_bundle(
        artifact_bundle, registration, host_authority=host["authority"])
    qualification = validate_qualification_binding(
        qualification_binding, registration, host,
        evidence=qualification_evidence)
    limit = registration["budgets"]["evidence_attempts"]
    attempt = _bounded_int(evidence_attempt, 1, limit, "evidence attempt")
    if attempt == 1:
        _require(prior_armed_chain_sha256 is None,
                 "first evidence attempt unexpectedly has a prior ARMED chain")
        prior = None
    else:
        prior = _sha256(
            prior_armed_chain_sha256, "prior ARMED chain hash")
    registration_sha = contract.canonical_sha256(registration)
    host_sha = contract.canonical_sha256(host)
    controller_sha = contract.canonical_sha256(
        registration["controller_bindings"])
    authority_bundle_sha = _sha256(
        authority_bundle_sha256, "ARMED authority-bundle hash")
    attempt_manifest_sha = _sha256(
        attempt_manifest_sha256, "ARMED attempt-manifest hash")
    lane_binding_sha = _sha256(
        lane_binding_sha256, "ARMED lane-binding hash")
    bridge = contract.strict_json_loads(
        qualification_evidence["bridge_data"], "ARMED qualification bridge")
    bridge_finished_ns = _bounded_int(
        bridge.get("bridge_finished_monotonic_ns"), 1, (1 << 63) - 1,
        "qualification bridge finish timestamp")
    bridge_deadline_ns = _bounded_int(
        bridge.get("bridge_deadline_monotonic_ns"), 1, (1 << 63) - 1,
        "qualification bridge deadline")
    armed_ns = _bounded_int(
        armed_monotonic_ns, 1, (1 << 63) - 1,
        "ARMED monotonic timestamp")
    _require(bridge_finished_ns <= armed_ns <= bridge_deadline_ns,
             "ARMED timestamp is outside the qualification handoff")
    _require(
        logical_plan["preregistration_sha256"] == registration_sha and
        artifacts["preregistration_sha256"] == registration_sha ==
            qualification["preregistration_sha256"] and
        artifacts["controller_bindings_sha256"] == controller_sha and
        qualification["host_instance_sha256"] == host_sha,
        "ARMED inputs do not share the preregistration, controller, or host joins",
    )
    return {
        "schema": ARMED_SCHEMA,
        "generation": GENERATION,
        "state": "ARMED",
        "evidence_attempt": attempt,
        "evidence_attempt_limit": limit,
        "prior_armed_chain_sha256": prior,
        "preregistration_sha256": registration_sha,
        "plan_sha256": contract.canonical_sha256(logical_plan),
        "host_instance_sha256": host_sha,
        "artifact_bundle_sha256": contract.canonical_sha256(artifacts),
        "qualification_binding_sha256": contract.canonical_sha256(qualification),
        "authority_bundle_sha256": authority_bundle_sha,
        "attempt_manifest_sha256": attempt_manifest_sha,
        "lane_binding_sha256": lane_binding_sha,
        "controller_bindings_sha256": controller_sha,
        "candidate_source": copy.deepcopy(artifacts["candidate_source"]),
        "selected_pair": copy.deepcopy(qualification["selected_pair"]),
        "armed_monotonic_ns": armed_ns,
        "bridge_deadline_monotonic_ns": bridge_deadline_ns,
        "publication_rule": PUBLICATION_RULE,
        "crash_after_publication_consumes": True,
    }


def validate_armed_record(
    value: Any, preregistration: Mapping[str, Any], plan: Mapping[str, Any],
    host_instance: Mapping[str, Any], artifact_bundle: Mapping[str, Any],
    qualification_binding: Mapping[str, Any], *,
    qualification_evidence: Mapping[str, Any],
) -> dict[str, Any]:
    record = validate_armed_record_shape(value)
    expected = armed_record(
        preregistration, plan, host_instance, artifact_bundle,
        qualification_binding,
        qualification_evidence=qualification_evidence,
        evidence_attempt=record["evidence_attempt"],
        prior_armed_chain_sha256=record["prior_armed_chain_sha256"],
        authority_bundle_sha256=record["authority_bundle_sha256"],
        attempt_manifest_sha256=record["attempt_manifest_sha256"],
        lane_binding_sha256=record["lane_binding_sha256"],
        armed_monotonic_ns=record["armed_monotonic_ns"])
    _require(contract.exact_json_equal(record, expected),
             "ARMED record differs from its fully replayed fixed point")
    return copy.deepcopy(expected)


__all__ = (
    "ARMED_SCHEMA", "ARTIFACT_BUNDLE_SCHEMA", "ExecutionContractError",
    "CANDIDATE_SOURCE_AUTHORITY_SCHEMA",
    "EXACT_MAIN_NORMALIZED_COMBINED_SHA256",
    "EXACT_MAIN_PATH_VARIANT_RAW_SHA256", "HOST_AUTHORITY_SCHEMA",
    "HOST_INSTANCE_SCHEMA", "PREREGISTRATION_SCHEMA_EXPECTATION",
    "PUBLICATION_RULE", "QUALIFICATION_BINDING_SCHEMA", "armed_record",
    "artifact_bundle_record", "candidate_source_authority_record",
    "exact_main_authority_record",
    "host_authority_record", "host_instance_record",
    "qualification_binding_record", "validate_armed_record",
    "validate_armed_record_shape",
    "validate_artifact_bundle", "validate_candidate_source_authority",
    "validate_exact_main_authority",
    "validate_host_authority", "validate_host_instance",
    "validate_qualification_binding",
)
