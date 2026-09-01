#!/usr/bin/env python3
"""Pure generation-3 preregistration and closed-generation lineage contract.

The module deliberately does not discover topology, read CPU counters, acquire a
lease, build a binary, or launch a benchmark.  It separates the mechanically
verifiable generation-2 history from the policy choices that still require explicit
authority.  A template is therefore never accepted as an executable
preregistration.
"""

from __future__ import annotations

import argparse
import copy
import ctypes
import errno
import hashlib
import importlib.util
import inspect
import os
from pathlib import Path, PurePosixPath
import re
import secrets
import stat
import sys
from typing import Any, Mapping, NoReturn, Sequence


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
MAIN_COMPARE = HERE.parent / "main_compare"
PAIR_CONTRACT_PATH = MAIN_COMPARE / "pair_qualification_contract.py"
PAIR_BRIDGE_PATH = MAIN_COMPARE / "pair_qualification_bridge_contract.py"


def _load_pair_contract() -> Any:
    module_name = "leopard2_pair_qualification_contract_for_k65_gen3_prereg"
    expected = PAIR_CONTRACT_PATH.resolve(strict=True)
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        if Path(str(getattr(loaded, "__file__", ""))).resolve() != expected:
            raise RuntimeError("pair qualification contract came from another path")
        return loaded
    specification = importlib.util.spec_from_file_location(module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load pair qualification contract")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    if Path(str(getattr(module, "__file__", ""))).resolve() != expected:
        raise RuntimeError("pair qualification contract resolved to another path")
    return module


contract = _load_pair_contract()


def _load_bridge_contract() -> Any:
    module_name = "leopard2_bridge_contract_for_k65_gen3_prereg"
    expected = PAIR_BRIDGE_PATH.resolve(strict=True)
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        if Path(str(getattr(loaded, "__file__", ""))).resolve() != expected:
            raise RuntimeError("pair qualification bridge came from another path")
        return loaded
    specification = importlib.util.spec_from_file_location(module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load pair qualification bridge")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    if Path(str(getattr(module, "__file__", ""))).resolve() != expected:
        raise RuntimeError("pair qualification bridge resolved to another path")
    return module


bridge_contract = _load_bridge_contract()

PREREGISTRATION_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-preregistration/v4"
TEMPLATE_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-preregistration-template/v4"
OUTPUT_LANE_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-output-lane/v4"
OUTPUT_LANE_FILE_HANDLE_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-output-lane-file-handle/v1"
OUTPUT_LANE_MANIFEST_BINDING_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-lane-manifest-binding/v2"
OUTPUT_LANE_TOPOLOGY_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-lane-binding/v4"
OUTPUT_LANE_PUBLICATION_RULE = \
    "preallocated-exact-payload-checkpoints-no-root-mutation/v2"
OUTPUT_LANE_MANIFEST_FILE = "lane.json"
OUTPUT_LANE_LOCK_FILE = ".arming.lock"
OUTPUT_LANE_ATTEMPTS_DIRECTORY = "attempts"
OUTPUT_LANE_CHILD_PROCESS_COUNT = 1650
OUTPUT_LANE_PREARMED_BOUNDARY_STATES = (
    "PREREGISTERED", "QUALIFYING", "QUALIFIED", "BRIDGING", "BRIDGED",
    "ARMING", "PRESAMPLING",
)
CLOSED_GENERATION_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-closed-generation/v1"
REPORTING_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-conditioned-reporting/v1"
NO_POOLING_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-no-pooling/v1"
TRACK_SELECTION_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-qualification-track-selection/v1"
LAUNCH_CONTEXT_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-gen3-launch-context/v1"

GENERATION = 3
V2_GENERATION = 2
V2_ATTEMPT_BUDGET = 3
V2_RUNNER_RELATIVE = \
    "experiments/leopard2/gf8_high_encode/run_k65r65_b64_packed_terminal_abba.py"
V2_RUNNER_SHA256 = \
    "2a5489ba7c1866135e5fc1577c3f4290e851bb41837d0e46e1118ac7699397ca"
V2_SOURCE_COMMIT = "c1b9fcfd6cc798c92103dce437038d5a22e56d25"
V2_SOURCE_TREE = "d46f1fb3688666410c4ca2eb98fa2c522eeab47f"
V2_TERMINAL_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-attempt-terminal/v2"
EXACT_MAIN_AUTHORITY_RECORD_SHA256 = \
    "ce3d6e647fd098f558ee01855523c0643b8c175578315d02a168a87992e5ae01"
EXACT_MAIN_SOURCE_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
EXACT_MAIN_SOURCE_TREE = "b7c8830d96a978f6ec14fe747095f066e351ae72"
EXACT_MAIN_ADAPTER_COMMIT = "6d4f690fe5ba7cf08f2f8f80a263765266b462e0"
V2_CANDIDATE_EXECUTABLE_SHA256 = \
    "07f364a68df55610cafe23ce22b22d5e7e5f1a1d303d357b8f1373bbabc93bc3"
EXACT_MAIN_EXECUTABLE_SHA256 = \
    "e1d849056a7f061b127e14c0d8f71165ceaa9065f2d28d7469cae33e7e19eadf"
V2_MATRIX_SHA256 = \
    "ab9572c4101b2af5eda4b7cfab17e979239f698c4c1196e660cf6f5e3f4af27c"
REQUIRED_CONTROLLER_PATHS = (
    "experiments/leopard2/main_compare/pair_qualification_contract.py",
    "experiments/leopard2/main_compare/pair_qualification_acquire.py",
    "experiments/leopard2/main_compare/pair_qualification_bridge_contract.py",
    "experiments/leopard2/main_compare/pair_qualification_verify.py",
    "experiments/leopard2/main_compare/git_capture.py",
    "experiments/leopard2/gf8_high_encode/k65_gen3_preregistration.py",
    "experiments/leopard2/gf8_high_encode/k65_gen3_exact_source_bootstrap.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k65r65_b64_packed_terminal_gen3_abba.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k65r65_b64_packed_terminal_abba.py",
    "experiments/leopard2/gf8_high_encode/k65_gen3_execution_contract.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k65r65_b64_packed_terminal_gen3_acquire.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k66r16_b64_tail_abba.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k9r5_b1024_terminal_abba.py",
    "experiments/leopard2/gf8_high_encode/"
    "run_k5r5_b64_terminal_abba.py",
    "experiments/leopard2/gf8_high_encode/run_t8_two_block_abba.py",
    "experiments/leopard2/main_compare/run_abba.py",
    "experiments/leopard2/decoder_dispatch/balanced_evidence_common.py",
    "tools/leopard2_build_provenance.py",
    "tools/leopard2_exact_main_baseline_acquire.py",
    "tools/leopard2_exact_main_baseline.py",
    "tools/leopard2_exact_main_baseline_record.py",
    "tools/leopard2_exact_main_baseline_verifier.py",
)

_ATTEMPTS = (
    {
        "attempt": 1,
        "failure_class": "environment_rejected",
        "failure_message": "CPU pair was not quiet during the presample",
        "terminal_sha256":
            "8b9728b41619ff3f7c7ea57dfe7d2a97057df2639807a665eda3b0236e7e5a48",
        "failure_sha256":
            "0a0dfa4d8ef81699a14586d2f1a232910c15aef87407373ae26771e371b44561",
    },
    {
        "attempt": 2,
        "failure_class": "setup_invalid",
        "failure_message":
            "production source is dirty or differs from the requested commit/tree",
        "terminal_sha256":
            "4cae64c7214d15efa7eb71dc794e96fc09cd493fa4de3375655cd62e1b9b70e1",
        "failure_sha256":
            "ee5d211777dd9b973efc70a3b2e00759af1668712fb0f6ebef7b89253e25b287",
    },
    {
        "attempt": 3,
        "failure_class": "contaminated_in_campaign",
        "failure_message":
            "contaminated target-k65-r65-b64-q1 round 4 for 3 attempts",
        "terminal_sha256":
            "48dc4e968ba11a9defcf6d1b87f35253cf5e168a56e28cc454e59bc39c2eb206",
        "failure_sha256":
            "086f8e063d579f8422b4e286e9b02b22bbf4c6db60d462f279f60c378d8d5dc1",
    },
)
_LINEAGE_FILES = (
    {
        "attempt": 2,
        "path":
            ".research/leopard-79h/c1b9fcf-k65r65-b64-v2-attempt2-lineage.json",
        "sha256":
            "a33a13a33caa549640c579945f9494bd03d0dd4c7e9478a356fa1818e4b9aaae",
    },
    {
        "attempt": 3,
        "path":
            ".research/leopard-79h/c1b9fcf-k65r65-b64-v2-attempt3-lineage.json",
        "sha256":
            "d3a3b0ee8fde0827d0c56e3e9dbbfa2833e2364afc50b0f4e01d54a0d7ed3dbc",
    },
)

ATTEMPT_TERMINAL_KEYS = frozenset((
    "schema", "generation", "attempt", "attempt_budget", "outcome",
    "promotable", "output_root", "summary_schema", "summary_sha256",
    "raw_sha256", "failure_sha256", "authority_record_sha256",
))
ATTEMPT_LINEAGE_KEYS = frozenset((
    "attempt", "failure_class", "lane", "terminal_path", "terminal_sha256",
    "failure_path", "failure_sha256",
))
LINEAGE_FILE_KEYS = frozenset(("attempt", "path", "sha256"))
CLOSED_GENERATION_KEYS = frozenset((
    "schema", "generation", "source_commit", "source_tree", "runner_path",
    "runner_sha256", "attempt_budget", "attempts", "lineage_files",
    "exact_main_authority", "generation_exhausted", "valid_ratio_count",
    "conclusion",
))
AUTHORITY_KEYS = frozenset((
    "record_sha256", "source_commit", "source_tree", "adapter_commit",
))
DENYLIST_ENTRY_KEYS = frozenset(("label", "sha256"))
NO_POOLING_KEYS = frozenset((
    "schema", "rule", "prior_generation_role", "prior_timing_reused",
    "prior_timing_resummarized", "prior_timing_pooled", "denylisted_artifacts",
))
REPORTING_KEYS = frozenset((
    "schema", "evidence_category", "universal_statements", "correlation_warning",
    "scope_warning", "track_a_claim_ceiling", "track_b_claim_ceiling",
    "track_a_ratio_field", "track_b_ratio_field", "prohibited_unqualified_words",
))
CLAIM_KEYS = frozenset((
    "promotion_eligible", "host_exclusivity_proved",
    "whole_campaign_interval_observed", "causal_performance_claim_allowed",
))
CAMPAIGN_KEYS = frozenset((
    "generation", "matrix_sha256", "cell_count", "rounds_per_cell",
    "iterations_per_child", "warmup_per_child", "reuse_per_child",
    "expected_child_process_count", "inference_unit", "pooling", "trimming",
    "prior_generation_timing_used", "statistical_campaign_retries",
    "target_control_ci95_lower_floor", "target_main_ci95_lower_floor",
    "neighbor_selector_ci95_equivalence_band",
    "retained_exact_main_ci95_lower_floor",
    "maximum_isolation_attempts_per_round", "confidence_level",
    "student_t_degrees_of_freedom", "student_t_two_sided_critical",
))
BUDGET_KEYS = frozenset((
    "setup_invalid", "environment_rejected", "evidence_attempts",
    "evidence_budget_commit_state", "child_launch_requires_prior_budget_commit",
    "crash_after_evidence_commit_consumes",
))
GEOMETRY_KEYS = frozenset((
    "scan_window_count", "scan_nominal_window_ns",
    "bridge_minimum_window_count", "bridge_maximum_window_count",
    "bridge_nominal_window_ns", "maximum_handoff_elapsed_ns",
))
QUALIFICATION_KEYS = frozenset((
    "policy_evaluation_order", "policies", "policy_sha256s", "track_b_permitted",
    "geometry", "freeze_point", "selection_rule", "frozen_pair_rule",
    "pre_campaign_only", "candidate_timing_performed",
))
CANDIDATE_BINDING_KEYS = frozenset((
    "mode", "sha256", "size", "build_provenance_sha256",
    "reproducible_build_core_sha256", "authority_record_sha256",
    "authority_ledger_sha256",
))
SOURCE_KEYS = frozenset(("commit", "tree"))
HOST_AUTHORITY_KEYS = frozenset((
    "machine_id_sha256", "hostname", "architecture", "cpu_model",
))
CONTROLLER_BINDING_KEYS = frozenset(("path", "sha256"))
OUTPUT_LANE_KEYS = frozenset((
    "schema", "path", "device", "inode", "uid", "mode", "link_count",
    "initial_mtime_ns", "initial_ctime_ns", "file_handle", "lane_manifest",
))
OUTPUT_LANE_FILE_HANDLE_KEYS = frozenset((
    "schema", "handle_type", "handle_hex",
))
OUTPUT_LANE_MANIFEST_BINDING_KEYS = frozenset((
    "schema", "name", "device", "inode", "uid", "mode", "link_count",
    "initial_mtime_ns", "initial_ctime_ns", "sha256", "size", "file_handle",
))
LAUNCH_CONTEXT_KEYS = frozenset((
    "schema", "nice", "scheduler_policy", "scheduler_priority", "umask",
    "signal_mask", "signal_dispositions", "cpu_affinity_rule", "rlimits",
    "launcher_rlimit_overrides",
))
RESOURCE_LIMIT_KEYS = frozenset(("name", "soft", "hard"))
RESOURCE_LIMIT_NAMES = (
    "RLIMIT_AS", "RLIMIT_CORE", "RLIMIT_CPU", "RLIMIT_DATA", "RLIMIT_FSIZE",
    "RLIMIT_MEMLOCK", "RLIMIT_MSGQUEUE", "RLIMIT_NICE", "RLIMIT_NOFILE",
    "RLIMIT_NPROC", "RLIMIT_RSS", "RLIMIT_RTPRIO", "RLIMIT_RTTIME",
    "RLIMIT_SIGPENDING", "RLIMIT_STACK",
)
RATIFICATION_KEYS = frozenset(("status", "authority", "authorized_utc"))
PREREGISTRATION_KEYS = frozenset((
    "schema", "generation", "ratification", "closed_generation_2",
    "qualification", "budgets", "campaign", "candidate_executable",
    "candidate_source", "host_authority", "launch_context",
    "controller_bindings", "output_lane", "no_pooling", "reporting",
    "safe_to_execute",
))
TRACK_SELECTION_KEYS = frozenset((
    "schema", "generation", "policy_a_scan_sha256", "policy_b_scan_sha256",
    "policy_b_retained", "policy_b_consulted", "selected_track",
    "selected_pair", "selection_status", "claim_ceiling",
    "candidate_timing_performed",
))
TEMPLATE_KEYS = frozenset((
    "schema", "generation", "ratification_status", "safe_to_execute",
    "closed_generation_2", "campaign", "no_pooling", "reporting",
    "required_decisions",
))
DECISION_KEYS = frozenset(("name", "recommended_value", "authorized_value"))


class PreregistrationError(ValueError):
    """A generation-3 preregistration or lineage binding is invalid."""


def _fail(message: str) -> NoReturn:
    raise PreregistrationError(message)


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


def _sha256(value: Any, label: str) -> str:
    _require(type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value) is not None,
             f"{label} is not a lowercase SHA-256")
    return value


def _git_oid(value: Any, label: str) -> str:
    _require(type(value) is str and re.fullmatch(r"[0-9a-f]{40}", value) is not None,
             f"{label} is not a lowercase 40-hex object ID")
    return value


def _bounded_ascii(value: Any, maximum: int, label: str) -> str:
    _require(type(value) is str and value.strip() == value and
             0 < len(value) <= maximum and
             all(0x20 <= ord(character) <= 0x7E for character in value),
             f"{label} is not bounded printable ASCII")
    return value


class _LinuxFileHandleHeader(ctypes.Structure):
    _fields_ = (
        ("handle_bytes", ctypes.c_uint),
        ("handle_type", ctypes.c_int),
    )


def _capture_output_lane_file_handle(
    directory_fd: int, name: bytes, label: str,
) -> tuple[dict[str, Any], int]:
    """Capture Linux's opaque inode-generation handle for one directory."""
    _require(sys.platform.startswith("linux") and
             type(directory_fd) is int and type(name) is bytes and name and
             b"\0" not in name,
             f"{label} file-handle request is invalid")
    libc = ctypes.CDLL(None, use_errno=True)
    operation = getattr(libc, "name_to_handle_at", None)
    _require(operation is not None,
             f"{label} filesystem has no durable file-handle API")
    operation.argtypes = (
        ctypes.c_int, ctypes.c_char_p, ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    )
    operation.restype = ctypes.c_int
    probe = _LinuxFileHandleHeader(0, 0)
    mount_id = ctypes.c_int(-1)
    ctypes.set_errno(0)
    result = operation(
        directory_fd, name, ctypes.byref(probe), ctypes.byref(mount_id), 0)
    error_number = ctypes.get_errno()
    _require(result == -1 and error_number == errno.EOVERFLOW and
             1 <= probe.handle_bytes <= 128,
             f"{label} filesystem cannot export a bounded durable handle")
    header_size = ctypes.sizeof(_LinuxFileHandleHeader)
    storage = ctypes.create_string_buffer(header_size + probe.handle_bytes)
    header = _LinuxFileHandleHeader.from_buffer(storage)
    header.handle_bytes = probe.handle_bytes
    ctypes.set_errno(0)
    result = operation(
        directory_fd, name, ctypes.byref(storage), ctypes.byref(mount_id), 0)
    error_number = ctypes.get_errno()
    _require(result == 0 and error_number == 0 and
             header.handle_bytes == probe.handle_bytes and
             1 <= mount_id.value <= (1 << 31) - 1 and
             0 <= header.handle_type <= (1 << 31) - 1,
             f"{label} durable file handle could not be captured")
    handle = bytes(storage)[header_size:header_size + header.handle_bytes]
    _require(len(handle) == header.handle_bytes,
             f"{label} durable file handle is truncated")
    return {
        "schema": OUTPUT_LANE_FILE_HANDLE_SCHEMA,
        "handle_type": header.handle_type,
        "handle_hex": handle.hex(),
    }, mount_id.value


def _validate_output_lane_file_handle(value: Any) -> dict[str, Any]:
    record = _exact_object(
        value, OUTPUT_LANE_FILE_HANDLE_KEYS,
        "generation-3 output-lane file handle")
    handle = record["handle_hex"]
    _require(
        record["schema"] == OUTPUT_LANE_FILE_HANDLE_SCHEMA and
        type(record["handle_type"]) is int and
        0 <= record["handle_type"] <= (1 << 31) - 1 and
        type(handle) is str and 2 <= len(handle) <= 256 and
        len(handle) % 2 == 0 and
        re.fullmatch(r"[0-9a-f]+", handle) is not None,
        "generation-3 output-lane file handle differs")
    return copy.deepcopy(record)


def _validate_output_lane_manifest_binding(value: Any) -> dict[str, Any]:
    record = copy.deepcopy(_exact_object(
        value, OUTPUT_LANE_MANIFEST_BINDING_KEYS,
        "generation-3 output-lane manifest binding"))
    _require(
        record["schema"] == OUTPUT_LANE_MANIFEST_BINDING_SCHEMA and
        record["name"] == OUTPUT_LANE_MANIFEST_FILE and
        type(record["device"]) is int and record["device"] >= 0 and
        type(record["inode"]) is int and record["inode"] > 0 and
        type(record["uid"]) is int and record["uid"] == os.geteuid() and
        type(record["mode"]) is int and record["mode"] == 0o400 and
        type(record["link_count"]) is int and record["link_count"] == 1 and
        type(record["initial_mtime_ns"]) is int and
        record["initial_mtime_ns"] > 0 and
        type(record["initial_ctime_ns"]) is int and
        record["initial_ctime_ns"] > 0 and
        type(record["size"]) is int and 1 <= record["size"] <= 16 * 1024 * 1024,
        "generation-3 output-lane manifest binding differs")
    record["sha256"] = _sha256(
        record["sha256"], "generation-3 output-lane manifest hash")
    record["file_handle"] = _validate_output_lane_file_handle(
        record["file_handle"])
    return copy.deepcopy(record)


def _validate_output_lane_binding(value: Any) -> dict[str, Any]:
    record = copy.deepcopy(_exact_object(
        value, OUTPUT_LANE_KEYS, "generation-3 output-lane binding"))
    path_text = _bounded_ascii(
        record["path"], 4096, "generation-3 output-lane path")
    path = Path(path_text)
    _require(
        record["schema"] == OUTPUT_LANE_SCHEMA and path.is_absolute() and
        str(path) == path_text and ".." not in path.parts and
        type(record["device"]) is int and
        0 <= record["device"] <= (1 << 64) - 1 and
        type(record["inode"]) is int and
        0 < record["inode"] <= (1 << 64) - 1 and
        type(record["uid"]) is int and 0 <= record["uid"] <= (1 << 32) - 1 and
        type(record["mode"]) is int and record["mode"] == 0o500 and
        type(record["link_count"]) is int and
        3 <= record["link_count"] <= 1003 and
        type(record["initial_mtime_ns"]) is int and
        record["initial_mtime_ns"] > 0 and
        type(record["initial_ctime_ns"]) is int and
        record["initial_ctime_ns"] > 0,
        "generation-3 output-lane binding differs")
    record["file_handle"] = _validate_output_lane_file_handle(
        record["file_handle"])
    record["lane_manifest"] = _validate_output_lane_manifest_binding(
        record["lane_manifest"])
    return copy.deepcopy(record)


def _lane_regular_slot_binding(
    metadata: os.stat_result, name: str,
) -> dict[str, Any]:
    _require(stat.S_ISREG(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and metadata.st_nlink == 1 and
             stat.S_IMODE(metadata.st_mode) == 0o400 and metadata.st_size == 0,
             f"generation-3 lane slot {name} is unsafe")
    return {
        "name": name,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o400,
        "link_count": 1,
        "initial_mtime_ns": metadata.st_mtime_ns,
        "initial_ctime_ns": metadata.st_ctime_ns,
    }


def _create_lane_regular_slot(lane_fd: int, name: str) -> dict[str, Any]:
    descriptor = -1
    try:
        descriptor = os.open(
            name, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            os.O_NOFOLLOW, 0o600, dir_fd=lane_fd)
        os.fchmod(descriptor, 0o400)
        # One post-chmod fsync persists both the newly created empty inode and
        # its final sealed mode; the containing directory is fsynced after the
        # complete slot inventory is linked.
        os.fsync(descriptor)
        return _lane_regular_slot_binding(os.fstat(descriptor), name)
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _lane_lock_binding(metadata: os.stat_result) -> dict[str, Any]:
    _require(stat.S_ISREG(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and metadata.st_nlink == 1 and
             stat.S_IMODE(metadata.st_mode) == 0o600 and metadata.st_size == 0,
             "generation-3 preallocated lane lock is unsafe")
    return {
        "name": OUTPUT_LANE_LOCK_FILE,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o600,
        "link_count": 1,
        "initial_mtime_ns": metadata.st_mtime_ns,
        "initial_ctime_ns": metadata.st_ctime_ns,
    }


def _lane_attempts_binding(
    metadata: os.stat_result, file_handle: Mapping[str, Any],
) -> dict[str, Any]:
    _require(stat.S_ISDIR(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and metadata.st_nlink == 2 and
             stat.S_IMODE(metadata.st_mode) == 0o500,
             "generation-3 preallocated attempts directory is unsafe")
    return {
        "name": OUTPUT_LANE_ATTEMPTS_DIRECTORY,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o500,
        "initial_link_count": 2,
        "file_handle": _validate_output_lane_file_handle(file_handle),
    }


def _journal_checkpoint_name(sequence: int) -> str:
    _require(type(sequence) is int and
             0 <= sequence <= 2 * OUTPUT_LANE_CHILD_PROCESS_COUNT + 1,
             "generation-3 journal checkpoint sequence is invalid")
    if sequence == 0:
        target = "armed"
    elif sequence == 2 * OUTPUT_LANE_CHILD_PROCESS_COUNT + 1:
        target = "complete"
    else:
        child_index = (sequence - 1) // 2
        kind = "intent" if sequence % 2 == 1 else "result"
        target = f"{kind}-{child_index:04d}"
    return f"checkpoint-{sequence:04d}-{target}.marker"


def _lane_progress_directory_binding(
    metadata: os.stat_result, *, name: str,
    file_handle: Mapping[str, Any], markers: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    _require(stat.S_ISDIR(metadata.st_mode) and
             metadata.st_uid == os.geteuid() and metadata.st_nlink == 2 and
             stat.S_IMODE(metadata.st_mode) == 0o500,
             "generation-3 journal checkpoint directory is unsafe")
    return {
        "name": name,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "uid": metadata.st_uid,
        "mode": 0o500,
        "link_count": 2,
        "initial_mtime_ns": metadata.st_mtime_ns,
        "initial_ctime_ns": metadata.st_ctime_ns,
        "file_handle": _validate_output_lane_file_handle(file_handle),
        "markers": copy.deepcopy(list(markers)),
    }


def _lane_topology_record(
    lane: Path, *, lane_nonce: str, setup_invalid_budget: int,
    environment_rejected_budget: int, evidence_attempt_budget: int,
    arming_lock: Mapping[str, Any], attempts_directory: Mapping[str, Any],
    campaign_markers: Sequence[Mapping[str, Any]],
    campaign_progress_directories: Sequence[Mapping[str, Any]],
    prearmed_history_markers: Sequence[Mapping[str, Any]],
    prearmed_boundary_markers: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    return {
        "schema": OUTPUT_LANE_TOPOLOGY_SCHEMA,
        "generation": GENERATION,
        "lane_nonce": lane_nonce,
        "absolute_path": str(lane),
        "budgets": {
            "setup_invalid": setup_invalid_budget,
            "environment_rejected": environment_rejected_budget,
            "evidence_attempts": evidence_attempt_budget,
        },
        "arming_lock": copy.deepcopy(arming_lock),
        "attempts_directory": copy.deepcopy(attempts_directory),
        "campaign_markers": copy.deepcopy(list(campaign_markers)),
        "campaign_progress_directories": copy.deepcopy(
            list(campaign_progress_directories)),
        "prearmed_history_markers": copy.deepcopy(
            list(prearmed_history_markers)),
        "prearmed_boundary_markers": copy.deepcopy(
            list(prearmed_boundary_markers)),
        "publication_rule": OUTPUT_LANE_PUBLICATION_RULE,
    }


def output_lane_binding_record(
    output_lane: Path, *, setup_invalid_budget: int,
    environment_rejected_budget: int, evidence_attempt_budget: int,
) -> dict[str, Any]:
    """Preallocate and seal the complete lane topology before ratification."""
    _require(
        isinstance(output_lane, Path) and output_lane.is_absolute() and
        not output_lane.is_symlink(),
        "generation-3 output lane is not a canonical directory")
    setup_budget = _bounded_int(
        setup_invalid_budget, 1, 1000, "setup-invalid budget")
    environment_budget = _bounded_int(
        environment_rejected_budget, 1, 1000,
        "environment-rejected budget")
    evidence_budget = _bounded_int(
        evidence_attempt_budget, 1, 1000, "evidence-attempt budget")
    _require(evidence_budget == 3,
             "evidence-attempt budget must preserve the frozen authority budget")
    prearmed_limit = setup_budget + environment_budget + evidence_budget
    descriptor = -1
    lock_fd = -1
    attempts_fd = -1
    manifest_fd = -1
    try:
        resolved = output_lane.resolve(strict=True)
        _require(output_lane == resolved,
                 "generation-3 output lane is not a canonical directory")
        descriptor = os.open(
            output_lane,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW)
        retained = os.fstat(descriptor)
        pathname = os.stat(output_lane, follow_symlinks=False)
        _require(
            stat.S_ISDIR(retained.st_mode) and
            stat.S_ISDIR(pathname.st_mode) and
            (retained.st_dev, retained.st_ino) ==
                (pathname.st_dev, pathname.st_ino) and
            retained.st_uid == pathname.st_uid == os.geteuid() and
            retained.st_mtime_ns == pathname.st_mtime_ns and
            retained.st_ctime_ns == pathname.st_ctime_ns and
            stat.S_IMODE(retained.st_mode) == 0o700 and
            retained.st_nlink == pathname.st_nlink == 2 and
            not os.listdir(descriptor),
            "generation-3 output lane is not an empty owner-private directory")

        lock_fd = os.open(
            OUTPUT_LANE_LOCK_FILE,
            os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            os.O_NOFOLLOW, 0o600, dir_fd=descriptor)
        os.fchmod(lock_fd, 0o600)
        arming_lock = _lane_lock_binding(os.fstat(lock_fd))
        os.mkdir(OUTPUT_LANE_ATTEMPTS_DIRECTORY, 0o700, dir_fd=descriptor)
        attempts_fd = os.open(
            OUTPUT_LANE_ATTEMPTS_DIRECTORY,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=descriptor)
        os.fchmod(attempts_fd, 0o500)
        attempts_handle, attempts_mount_id = _capture_output_lane_file_handle(
            descriptor, os.fsencode(OUTPUT_LANE_ATTEMPTS_DIRECTORY),
            "preallocated generation-3 attempts directory")
        attempts = _lane_attempts_binding(
            os.fstat(attempts_fd), attempts_handle)

        campaign_markers = [
            _create_lane_regular_slot(
                descriptor, f"attempt-{index:03d}-transcript.jsonl")
            for index in range(1, evidence_budget + 1)
        ]
        campaign_progress_directories: list[dict[str, Any]] = []
        for index in range(1, evidence_budget + 1):
            progress_name = f"attempt-{index:03d}-journal-checkpoints"
            os.mkdir(progress_name, 0o700, dir_fd=descriptor)
            progress_fd = os.open(
                progress_name,
                os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW,
                dir_fd=descriptor)
            try:
                checkpoint_markers = [
                    _create_lane_regular_slot(
                        progress_fd, _journal_checkpoint_name(sequence))
                    for sequence in range(
                        0, 2 * OUTPUT_LANE_CHILD_PROCESS_COUNT + 2)
                ]
                os.fsync(progress_fd)
                os.fchmod(progress_fd, 0o500)
                os.fsync(progress_fd)
                progress_handle, progress_mount_id = \
                    _capture_output_lane_file_handle(
                        descriptor, os.fsencode(progress_name),
                        f"generation-3 journal checkpoints {index}")
                _require(progress_mount_id == attempts_mount_id,
                         "generation-3 checkpoint directory crossed a mount")
                campaign_progress_directories.append(
                    _lane_progress_directory_binding(
                        os.fstat(progress_fd), name=progress_name,
                        file_handle=progress_handle,
                        markers=checkpoint_markers))
            finally:
                os.close(progress_fd)
        history_markers = [
            _create_lane_regular_slot(
                descriptor, f"prearmed-{index:04d}-history.jsonl")
            for index in range(1, prearmed_limit + 1)
        ]
        boundary_markers = [
            _create_lane_regular_slot(
                descriptor,
                f"prearmed-{index:04d}-{state.lower()}-reached.marker")
            for index in range(1, prearmed_limit + 1)
            for state in OUTPUT_LANE_PREARMED_BOUNDARY_STATES
        ]
        topology = _lane_topology_record(
            resolved, lane_nonce=secrets.token_hex(32),
            setup_invalid_budget=setup_budget,
            environment_rejected_budget=environment_budget,
            evidence_attempt_budget=evidence_budget,
            arming_lock=arming_lock, attempts_directory=attempts,
            campaign_markers=campaign_markers,
            campaign_progress_directories=campaign_progress_directories,
            prearmed_history_markers=history_markers,
            prearmed_boundary_markers=boundary_markers)
        manifest_data = contract.canonical_json_bytes(topology)
        manifest_fd = os.open(
            OUTPUT_LANE_MANIFEST_FILE,
            os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            os.O_NOFOLLOW, 0o600, dir_fd=descriptor)
        written = 0
        while written < len(manifest_data):
            count = os.write(manifest_fd, manifest_data[written:])
            _require(count > 0,
                     "generation-3 lane manifest write made no progress")
            written += count
        os.fsync(manifest_fd)
        os.fchmod(manifest_fd, 0o400)
        os.fsync(manifest_fd)
        os.fsync(lock_fd)
        os.fsync(attempts_fd)
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o500)
        os.fsync(descriptor)

        retained = os.fstat(descriptor)
        pathname = os.stat(resolved, follow_symlinks=False)
        manifest_metadata = os.fstat(manifest_fd)
        manifest_path_metadata = os.stat(
            OUTPUT_LANE_MANIFEST_FILE, dir_fd=descriptor,
            follow_symlinks=False)
        expected_names = {
            OUTPUT_LANE_LOCK_FILE, OUTPUT_LANE_ATTEMPTS_DIRECTORY,
            OUTPUT_LANE_MANIFEST_FILE,
            *(item["name"] for item in campaign_markers),
            *(item["name"] for item in campaign_progress_directories),
            *(item["name"] for item in history_markers),
            *(item["name"] for item in boundary_markers),
        }
        _require(
            stat.S_ISDIR(retained.st_mode) and
            stat.S_IMODE(retained.st_mode) == 0o500 and
            retained.st_uid == pathname.st_uid == os.geteuid() and
            retained.st_nlink == pathname.st_nlink ==
                3 + evidence_budget and
            retained.st_mtime_ns == pathname.st_mtime_ns and
            retained.st_ctime_ns == pathname.st_ctime_ns and
            (retained.st_dev, retained.st_ino) ==
                (pathname.st_dev, pathname.st_ino) and
            set(os.listdir(descriptor)) == expected_names and
            stat.S_ISREG(manifest_metadata.st_mode) and
            stat.S_IMODE(manifest_metadata.st_mode) == 0o400 and
            manifest_metadata.st_uid == os.geteuid() and
            manifest_metadata.st_nlink == 1 and
            manifest_metadata.st_size == len(manifest_data) and
            (manifest_metadata.st_dev, manifest_metadata.st_ino) ==
                (manifest_path_metadata.st_dev, manifest_path_metadata.st_ino),
            "generation-3 output lane did not reach its sealed topology")
        retained_handle, retained_mount_id = _capture_output_lane_file_handle(
            descriptor, b".", "retained generation-3 output lane")
        pathname_handle, pathname_mount_id = _capture_output_lane_file_handle(
            -100, os.fsencode(resolved), "path generation-3 output lane")
        _require(retained_mount_id == pathname_mount_id and
                 contract.exact_json_equal(
                     retained_handle, pathname_handle),
            "generation-3 output lane file handles differ during capture")
        _require(retained_mount_id == attempts_mount_id,
                 "generation-3 attempts directory crossed a mount boundary")
        manifest_handle, manifest_mount_id = _capture_output_lane_file_handle(
            descriptor, os.fsencode(OUTPUT_LANE_MANIFEST_FILE),
            "retained generation-3 lane manifest")
        manifest_path_handle, manifest_path_mount_id = \
            _capture_output_lane_file_handle(
                -100, os.fsencode(resolved / OUTPUT_LANE_MANIFEST_FILE),
                "path generation-3 lane manifest")
        _require(
            retained_mount_id == manifest_mount_id == manifest_path_mount_id and
            contract.exact_json_equal(
                manifest_handle, manifest_path_handle),
            "generation-3 lane manifest file handles differ during capture")
        return {
            "schema": OUTPUT_LANE_SCHEMA,
            "path": str(resolved),
            "device": retained.st_dev,
            "inode": retained.st_ino,
            "uid": retained.st_uid,
            "mode": 0o500,
            "link_count": retained.st_nlink,
            "initial_mtime_ns": retained.st_mtime_ns,
            "initial_ctime_ns": retained.st_ctime_ns,
            "file_handle": retained_handle,
            "lane_manifest": {
                "schema": OUTPUT_LANE_MANIFEST_BINDING_SCHEMA,
                "name": OUTPUT_LANE_MANIFEST_FILE,
                "device": manifest_metadata.st_dev,
                "inode": manifest_metadata.st_ino,
                "uid": manifest_metadata.st_uid,
                "mode": 0o400,
                "link_count": 1,
                "initial_mtime_ns": manifest_metadata.st_mtime_ns,
                "initial_ctime_ns": manifest_metadata.st_ctime_ns,
                "sha256": hashlib.sha256(manifest_data).hexdigest(),
                "size": len(manifest_data),
                "file_handle": manifest_handle,
            },
        }
    except OSError as error:
        raise PreregistrationError(
            "generation-3 output lane cannot be captured") from error
    finally:
        if manifest_fd >= 0:
            os.close(manifest_fd)
        if attempts_fd >= 0:
            os.close(attempts_fd)
        if lock_fd >= 0:
            os.close(lock_fd)
        if descriptor >= 0:
            os.close(descriptor)


def recommended_launch_context_record() -> dict[str, Any]:
    """Return the exact child process context proposed for ratification."""
    limits = (
        ("RLIMIT_AS", 201_326_592, 201_326_592),
        ("RLIMIT_CORE", 0, -1),
        ("RLIMIT_CPU", -1, -1),
        ("RLIMIT_DATA", -1, -1),
        ("RLIMIT_FSIZE", -1, -1),
        ("RLIMIT_MEMLOCK", 67_514_826_752, 67_514_826_752),
        ("RLIMIT_MSGQUEUE", 819_200, 819_200),
        ("RLIMIT_NICE", 0, 0),
        ("RLIMIT_NOFILE", 1_048_576, 1_048_576),
        ("RLIMIT_NPROC", 2_059_714, 2_059_714),
        ("RLIMIT_RSS", -1, -1),
        ("RLIMIT_RTPRIO", 0, 0),
        ("RLIMIT_RTTIME", -1, -1),
        ("RLIMIT_SIGPENDING", 2_059_714, 2_059_714),
        ("RLIMIT_STACK", 8_388_608, -1),
    )
    return {
        "schema": LAUNCH_CONTEXT_SCHEMA,
        "nice": 0,
        "scheduler_policy": 0,
        "scheduler_priority": 0,
        "umask": 0o022,
        "signal_mask": "empty",
        "signal_dispositions": "default-except-sigkill-sigstop",
        "cpu_affinity_rule": "selected-benchmark-cpu-singleton-at-exec",
        "rlimits": [
            {"name": name, "soft": soft, "hard": hard}
            for name, soft, hard in limits
        ],
        "launcher_rlimit_overrides": [{
            "name": "RLIMIT_AS", "soft": 201_326_592,
            "hard": 201_326_592,
        }],
    }


def validate_launch_context_contract(value: Any) -> dict[str, Any]:
    record = _exact_object(
        value, LAUNCH_CONTEXT_KEYS, "generation-3 launch context")
    _require(all(type(record[name]) is int for name in (
                 "nice", "scheduler_policy", "scheduler_priority",
                 "umask")) and
             record["schema"] == LAUNCH_CONTEXT_SCHEMA and
             record["nice"] == 0 and
             record["scheduler_policy"] == 0 and
             record["scheduler_priority"] == 0 and
             record["umask"] == 0o022 and
             record["signal_mask"] == "empty" and
             record["signal_dispositions"] ==
                "default-except-sigkill-sigstop" and
             record["cpu_affinity_rule"] ==
                "selected-benchmark-cpu-singleton-at-exec",
             "generation-3 launch context normative policy differs")
    _require(type(record["rlimits"]) is list and
             len(record["rlimits"]) == len(RESOURCE_LIMIT_NAMES),
             "generation-3 launch-context resource limits differ")
    limits: list[dict[str, Any]] = []
    for expected_name, raw in zip(RESOURCE_LIMIT_NAMES, record["rlimits"]):
        limit = _exact_object(
            raw, RESOURCE_LIMIT_KEYS,
            f"generation-3 resource limit {expected_name}")
        soft, hard = limit["soft"], limit["hard"]
        _require(limit["name"] == expected_name and
                 type(soft) is int and type(hard) is int and
                 (soft == -1 or soft >= 0) and
                 (hard == -1 or hard >= 0) and
                 (hard == -1 or (soft != -1 and soft <= hard)),
                 f"generation-3 resource limit {expected_name} is invalid")
        limits.append({"name": expected_name, "soft": soft, "hard": hard})
    overrides = record["launcher_rlimit_overrides"]
    _require(type(overrides) is list and len(overrides) == 1,
             "generation-3 launcher resource-limit override differs")
    override = _exact_object(
        overrides[0], RESOURCE_LIMIT_KEYS,
        "generation-3 launcher resource-limit override")
    expected_override = {
        "name": "RLIMIT_AS", "soft": 201_326_592, "hard": 201_326_592,
    }
    _require(contract.exact_json_equal(override, expected_override) and
             contract.exact_json_equal(limits[0], expected_override),
             "generation-3 AS limit does not match its launcher override")
    return {
        "schema": LAUNCH_CONTEXT_SCHEMA,
        "nice": 0,
        "scheduler_policy": 0,
        "scheduler_priority": 0,
        "umask": 0o022,
        "signal_mask": "empty",
        "signal_dispositions": "default-except-sigkill-sigstop",
        "cpu_affinity_rule": "selected-benchmark-cpu-singleton-at-exec",
        "rlimits": limits,
        "launcher_rlimit_overrides": [expected_override],
    }


def _host_authority_record(
    *, machine_id_sha256: Any, hostname: Any, architecture: Any,
    cpu_model: Any,
) -> dict[str, str]:
    machine = _sha256(machine_id_sha256, "host machine-id hash")
    host = _bounded_ascii(hostname, 253, "host name")
    _require(re.fullmatch(
        r"(?=.{1,253}\Z)[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?"
        r"(?:\.[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?)*", host) is not None,
        "host name is not canonical lower-case DNS syntax")
    arch = _bounded_ascii(architecture, 64, "host architecture")
    _require(arch == "x86_64", "host architecture is not x86_64")
    model = _bounded_ascii(cpu_model, 512, "host CPU model")
    return {
        "machine_id_sha256": machine,
        "hostname": host,
        "architecture": arch,
        "cpu_model": model,
    }


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _strict_json_path(path: Path, label: str) -> Any:
    try:
        data = path.read_bytes()
    except OSError as error:
        raise PreregistrationError(f"cannot read {label}: {path}") from error
    return contract.strict_json_loads(data, label)


def _repo_relative(value: Any, label: str) -> str:
    _require(type(value) is str and value != "", f"{label} is not a path")
    path = PurePosixPath(value)
    _require(not path.is_absolute() and ".." not in path.parts and "." not in path.parts,
             f"{label} is not a canonical repository-relative path")
    return str(path)


def _lane_relative(attempt: int) -> str:
    return f".research/leopard-79h/c1b9fcf-k65r65-b64-v2-a{attempt}"


def current_controller_bindings(
    repo_root: Path = REPO_ROOT,
) -> list[dict[str, str]]:
    """Hash the exact reviewed controller set without inferring authority."""
    root = repo_root.resolve(strict=True)
    result = []
    for relative in sorted(REQUIRED_CONTROLLER_PATHS):
        path = root / relative
        _require(path.is_file() and not path.is_symlink(),
                 f"controller is not one regular file: {relative}")
        result.append({"path": relative, "sha256": _sha256_file(path)})
    return result


def _terminal_projection(value: Any, attempt: int, expected: Mapping[str, Any]) -> dict[str, Any]:
    terminal = _exact_object(value, ATTEMPT_TERMINAL_KEYS,
                             f"generation-2 attempt {attempt} terminal")
    _require(
        terminal["schema"] == V2_TERMINAL_SCHEMA and
        terminal["generation"] == V2_GENERATION and
        terminal["attempt"] == attempt and
        terminal["attempt_budget"] == V2_ATTEMPT_BUDGET and
        terminal["outcome"] == "failed" and
        terminal["promotable"] is False and
        terminal["summary_schema"] == "absent" and
        terminal["summary_sha256"] is None and
        terminal["raw_sha256"] is None and
        terminal["failure_sha256"] == expected["failure_sha256"] and
        terminal["authority_record_sha256"] == EXACT_MAIN_AUTHORITY_RECORD_SHA256,
        f"generation-2 attempt {attempt} terminal semantics differ")
    return copy.deepcopy(terminal)


def closed_generation2_lineage_record(
    repo_root: Path = REPO_ROOT, *, verify_files: bool = True,
) -> dict[str, Any]:
    """Return the exact, exhausted three-attempt generation-2 history."""
    root = repo_root.resolve(strict=True)
    attempts: list[dict[str, Any]] = []
    terminals: list[dict[str, Any]] = []
    for expected in _ATTEMPTS:
        attempt = int(expected["attempt"])
        lane = _lane_relative(attempt)
        terminal_rel = f"{lane}/attempt-terminal.json"
        failure_rel = f"{lane}/failure.json"
        if verify_files:
            lane_path = root / lane
            _require(lane_path.is_dir() and not lane_path.is_symlink(),
                     f"generation-2 lane a{attempt} is not retained")
            terminal_path = root / terminal_rel
            failure_path = root / failure_rel
            _require(terminal_path.is_file() and not terminal_path.is_symlink(),
                     f"generation-2 terminal a{attempt} is not regular")
            _require(failure_path.is_file() and not failure_path.is_symlink(),
                     f"generation-2 failure a{attempt} is not regular")
            _require(_sha256_file(terminal_path) == expected["terminal_sha256"],
                     f"generation-2 terminal a{attempt} hash differs")
            _require(_sha256_file(failure_path) == expected["failure_sha256"],
                     f"generation-2 failure a{attempt} hash differs")
            terminal = _terminal_projection(
                _strict_json_path(terminal_path, f"generation-2 terminal a{attempt}"),
                attempt, expected)
            failure = _strict_json_path(
                failure_path, f"generation-2 failure a{attempt}")
            _require(
                type(failure) is dict and
                failure.get("schema") ==
                    "leopard2-gf8-k65r65-b64-packed-terminal-abba/v2" and
                failure.get("generation") == V2_GENERATION and
                failure.get("source_commit") == V2_SOURCE_COMMIT and
                failure.get("source_tree") == V2_SOURCE_TREE and
                type(failure.get("failure")) is dict and
                failure["failure"].get("type") == "EvidenceError" and
                failure["failure"].get("message") == expected["failure_message"],
                f"generation-2 failure a{attempt} semantics differ")
            terminals.append(terminal)
        attempts.append({
            "attempt": attempt,
            "failure_class": expected["failure_class"],
            "lane": lane,
            "terminal_path": terminal_rel,
            "terminal_sha256": expected["terminal_sha256"],
            "failure_path": failure_rel,
            "failure_sha256": expected["failure_sha256"],
        })
    lineage_files = copy.deepcopy(list(_LINEAGE_FILES))
    if verify_files:
        for lineage in lineage_files:
            path = root / lineage["path"]
            _require(path.is_file() and not path.is_symlink(),
                     f"generation-2 lineage attempt {lineage['attempt']} is absent")
            _require(_sha256_file(path) == lineage["sha256"],
                     f"generation-2 lineage attempt {lineage['attempt']} hash differs")
            observed = _strict_json_path(path, "generation-2 attempt lineage")
            _require(contract.exact_json_equal(
                observed, terminals[: int(lineage["attempt"]) - 1]),
                f"generation-2 lineage attempt {lineage['attempt']} is not exact")
        runner = root / V2_RUNNER_RELATIVE
        _require(runner.is_file() and not runner.is_symlink() and
                 _sha256_file(runner) == V2_RUNNER_SHA256,
                 "frozen generation-2 runner hash differs")
    return {
        "schema": CLOSED_GENERATION_SCHEMA,
        "generation": V2_GENERATION,
        "source_commit": V2_SOURCE_COMMIT,
        "source_tree": V2_SOURCE_TREE,
        "runner_path": V2_RUNNER_RELATIVE,
        "runner_sha256": V2_RUNNER_SHA256,
        "attempt_budget": V2_ATTEMPT_BUDGET,
        "attempts": attempts,
        "lineage_files": lineage_files,
        "exact_main_authority": {
            "record_sha256": EXACT_MAIN_AUTHORITY_RECORD_SHA256,
            "source_commit": EXACT_MAIN_SOURCE_COMMIT,
            "source_tree": EXACT_MAIN_SOURCE_TREE,
            "adapter_commit": EXACT_MAIN_ADAPTER_COMMIT,
        },
        "generation_exhausted": True,
        "valid_ratio_count": 0,
        "conclusion":
            "generation 2 exhausted attempts 1-3 without valid ratios",
    }


def validate_closed_generation2_lineage(
    value: Any, repo_root: Path = REPO_ROOT, *, verify_files: bool = True,
) -> dict[str, Any]:
    _exact_object(value, CLOSED_GENERATION_KEYS, "closed generation-2 lineage")
    expected = closed_generation2_lineage_record(
        repo_root, verify_files=verify_files)
    _require(contract.exact_json_equal(value, expected),
             "closed generation-2 lineage differs")
    return copy.deepcopy(expected)


def no_pooling_record(closed_generation: Mapping[str, Any]) -> dict[str, Any]:
    validated = validate_closed_generation2_lineage(
        closed_generation, verify_files=False)
    artifacts = [
        {"label": "generation-2 candidate/control executable",
         "sha256": V2_CANDIDATE_EXECUTABLE_SHA256},
        {"label": "generation-2 exact-main executable",
         "sha256": EXACT_MAIN_EXECUTABLE_SHA256},
    ]
    for attempt in validated["attempts"]:
        artifacts.extend((
            {"label": f"generation-2 a{attempt['attempt']} terminal",
             "sha256": attempt["terminal_sha256"]},
            {"label": f"generation-2 a{attempt['attempt']} failure",
             "sha256": attempt["failure_sha256"]},
        ))
    for lineage in validated["lineage_files"]:
        artifacts.append({
            "label": f"generation-2 attempt-{lineage['attempt']} lineage",
            "sha256": lineage["sha256"],
        })
    artifacts.sort(key=lambda item: item["label"])
    _require(len({item["sha256"] for item in artifacts}) == len(artifacts),
             "generation-2 denylist contains a hash collision")
    return {
        "schema": NO_POOLING_SCHEMA,
        "rule": "generation-2-artifacts-are-lineage-only-never-evidence/v1",
        "prior_generation_role": "closed-failure-lineage-only",
        "prior_timing_reused": False,
        "prior_timing_resummarized": False,
        "prior_timing_pooled": False,
        "denylisted_artifacts": artifacts,
    }


def validate_no_pooling(value: Any, closed_generation: Mapping[str, Any]) -> dict[str, Any]:
    record = _exact_object(value, NO_POOLING_KEYS, "no-pooling contract")
    _require(type(record["denylisted_artifacts"]) is list,
             "no-pooling denylist is not a list")
    for index, item in enumerate(record["denylisted_artifacts"]):
        _exact_object(item, DENYLIST_ENTRY_KEYS, f"denylist entry {index}")
        _sha256(item["sha256"], f"denylist entry {index} hash")
    expected = no_pooling_record(closed_generation)
    _require(contract.exact_json_equal(record, expected),
             "no-pooling contract differs")
    return copy.deepcopy(expected)


def _shared_host_claim_ceiling() -> dict[str, bool]:
    """No pre-campaign scan establishes campaign-wide shared-host validity."""
    return {
        "promotion_eligible": False,
        "host_exclusivity_proved": False,
        "whole_campaign_interval_observed": False,
        "causal_performance_claim_allowed": False,
    }


def reporting_contract_record() -> dict[str, Any]:
    return {
        "schema": REPORTING_SCHEMA,
        "evidence_category":
            "conditioned passive shared-host observation",
        "universal_statements": [
            "generation 3, attempt N of at most 3",
            "generation 2 produced zero valid ratios",
            "generation-2 timings were not reused, re-summarised, or pooled",
        ],
        "correlation_warning":
            "batch-of-one and one-shot ratios share Leopard1 numerator observations; "
            "they are correlated estimates and must not be multiplied or stacked",
        "scope_warning":
            "every number is host-, compiler-, API-, and workload-specific",
        "track_a_claim_ceiling": _shared_host_claim_ceiling(),
        "track_b_claim_ceiling": _shared_host_claim_ceiling(),
        "track_a_ratio_field": "conditioned_nominal_ratio",
        "track_b_ratio_field": "conditioned_nominal_ratio",
        "prohibited_unqualified_words": ["faster", "speedup"],
    }


def validate_reporting_contract(value: Any) -> dict[str, Any]:
    record = _exact_object(value, REPORTING_KEYS, "reporting contract")
    _exact_object(record["track_a_claim_ceiling"], CLAIM_KEYS,
                  "Track-A claim ceiling")
    _exact_object(record["track_b_claim_ceiling"], CLAIM_KEYS,
                  "Track-B claim ceiling")
    expected = reporting_contract_record()
    _require(contract.exact_json_equal(record, expected),
             "reporting contract differs")
    return copy.deepcopy(expected)


def qualification_track_selection_record(
    preregistration: Mapping[str, Any], *, policy_a_scan: Any,
    policy_b_scan: Any, expected_frozen_pair: Any,
) -> dict[str, Any]:
    """Apply the frozen A-then-B order without using any performance outcome."""
    registration = validate_preregistration(
        preregistration, verify_files=False)
    qualification = registration["qualification"]
    policies = qualification["policies"]
    scan_a = contract.validate_pair_qualification_scan(
        policy_a_scan, expected_policy=policies[0],
        expected_frozen_pair=expected_frozen_pair)
    _require(
        scan_a["scan_window_count"] ==
            qualification["geometry"]["scan_window_count"],
        "Policy A scan window count differs from the preregistration")
    scan_b = None
    if qualification["track_b_permitted"]:
        _require(len(policies) == 2 and policy_b_scan is not None,
                 "permitted Policy B scan is absent")
        scan_b = contract.validate_pair_qualification_scan(
            policy_b_scan, expected_policy=policies[1],
            expected_frozen_pair=expected_frozen_pair)
        _require(
            scan_b["scan_window_count"] ==
                qualification["geometry"]["scan_window_count"],
            "Policy B scan window count differs from the preregistration")
        _require(
            contract.exact_json_equal(
                scan_a["allowed_cpu_set_at_launch"],
                scan_b["allowed_cpu_set_at_launch"]) and
            contract.exact_json_equal(
                scan_a["topology_before"], scan_b["topology_before"]) and
            contract.exact_json_equal(
                scan_a["topology_after"], scan_b["topology_after"]) and
            len(scan_a["windows"]) == len(scan_b["windows"]) and
            all(contract.exact_json_equal(a["before"], b["before"]) and
                contract.exact_json_equal(a["after"], b["after"])
                for a, b in zip(scan_a["windows"], scan_b["windows"])) and
            scan_a["scan_started_monotonic_ns"] ==
                scan_b["scan_started_monotonic_ns"] and
            scan_a["scan_finished_monotonic_ns"] ==
                scan_b["scan_finished_monotonic_ns"],
            "Policy A and Policy B do not share one retained snapshot chain")
    else:
        _require(len(policies) == 1 and policy_b_scan is None,
                 "forbidden Policy B scan was supplied")

    status_a = scan_a["selection_status"]
    selected_track = None
    selected_pair = None
    policy_b_consulted = False
    selection_status = status_a
    if status_a in {"selected-lowest-primary", "selected-frozen-pair"}:
        selected_track = "pair-and-domain"
        selected_pair = copy.deepcopy(scan_a["selected"])
    elif status_a == "no-candidate-pair-qualified" and scan_b is not None:
        policy_b_consulted = True
        status_b = scan_b["selection_status"]
        selection_status = status_b
        if status_b in {"selected-lowest-primary", "selected-frozen-pair"}:
            selected_track = "pair-only"
            selected_pair = copy.deepcopy(scan_b["selected"])
        else:
            _require(status_b in {
                "no-candidate-pair-qualified",
                "frozen-pair-did-not-requalify",
            }, "Policy B selection status differs")
    else:
        _require(status_a in {
            "no-candidate-pair-qualified",
            "frozen-pair-did-not-requalify",
        }, "Policy A selection status differs")

    return {
        "schema": TRACK_SELECTION_SCHEMA,
        "generation": GENERATION,
        "policy_a_scan_sha256": contract.canonical_sha256(scan_a),
        "policy_b_scan_sha256": contract.canonical_sha256(scan_b)
            if scan_b is not None else None,
        "policy_b_retained": scan_b is not None,
        "policy_b_consulted": policy_b_consulted,
        "selected_track": selected_track,
        "selected_pair": selected_pair,
        "selection_status": selection_status,
        "claim_ceiling": _shared_host_claim_ceiling(),
        "candidate_timing_performed": False,
    }


def validate_qualification_track_selection(
    value: Any, preregistration: Mapping[str, Any], *, policy_a_scan: Any,
    policy_b_scan: Any, expected_frozen_pair: Any,
) -> dict[str, Any]:
    record = _exact_object(value, TRACK_SELECTION_KEYS,
                           "qualification track selection")
    _exact_object(record["claim_ceiling"], CLAIM_KEYS,
                  "qualification track claim ceiling")
    expected = qualification_track_selection_record(
        preregistration, policy_a_scan=policy_a_scan,
        policy_b_scan=policy_b_scan,
        expected_frozen_pair=expected_frozen_pair)
    _require(contract.exact_json_equal(record, expected),
             "qualification track selection differs")
    return copy.deepcopy(expected)


def campaign_contract_record() -> dict[str, Any]:
    return {
        "generation": GENERATION,
        "matrix_sha256": V2_MATRIX_SHA256,
        "cell_count": 13,
        "rounds_per_cell": 25,
        "iterations_per_child": 31,
        "warmup_per_child": 64,
        "reuse_per_child": 8192,
        "expected_child_process_count": 1650,
        "inference_unit": "within-round log contrast",
        "pooling": "none",
        "trimming": "none",
        "prior_generation_timing_used": False,
        "statistical_campaign_retries": 0,
        "target_control_ci95_lower_floor": 1.05,
        "target_main_ci95_lower_floor": 1.05,
        "neighbor_selector_ci95_equivalence_band": [1.0 / 1.02, 1.02],
        "retained_exact_main_ci95_lower_floor": 0.98,
        "maximum_isolation_attempts_per_round": 3,
        "confidence_level": 0.95,
        "student_t_degrees_of_freedom": 24,
        "student_t_two_sided_critical": 2.0638985616280205,
    }


def validate_campaign_contract(value: Any) -> dict[str, Any]:
    record = _exact_object(value, CAMPAIGN_KEYS, "campaign contract")
    expected = campaign_contract_record()
    _require(contract.exact_json_equal(record, expected),
             "campaign contract differs")
    return copy.deepcopy(expected)


_PREREGISTRATION_PARAMETER_DECISIONS = {
    "authority": "generation_authorized",
    "authorized_utc": "generation_authorized",
    "clock_ticks_per_second": "clock_ticks_per_second",
    "candidate_primary_cpus": "candidate_primary_cpus",
    "excluded_pairs": "excluded_pairs",
    "track_b_permitted": "track_b_permitted",
    "setup_invalid_budget": "setup_invalid_budget",
    "environment_rejected_budget": "environment_rejected_budget",
    "evidence_attempt_budget": "evidence_attempt_budget",
    "scan_window_count": "scan_window_count",
    "scan_nominal_window_ns": "scan_nominal_window_ns",
    "bridge_minimum_window_count": "bridge_geometry",
    "bridge_maximum_window_count": "bridge_geometry",
    "bridge_nominal_window_ns": "bridge_geometry",
    "maximum_handoff_elapsed_ns": "bridge_geometry",
    "freeze_point": "pair_freeze_point",
    "candidate_executable_mode": "candidate_executable_binding",
    "candidate_executable_sha256": "candidate_executable_binding",
    "candidate_executable_size": "candidate_executable_binding",
    "candidate_build_provenance_sha256": "candidate_build_closure",
    "candidate_reproducible_build_core_sha256": "candidate_build_closure",
    "candidate_authority_record_sha256": "candidate_build_closure",
    "candidate_authority_ledger_sha256": "candidate_build_closure",
    "candidate_source_commit": "candidate_source_commit_and_tree",
    "candidate_source_tree": "candidate_source_commit_and_tree",
    "host_machine_id_sha256": "host_authority",
    "host_name": "host_authority",
    "host_architecture": "host_authority",
    "host_cpu_model": "host_authority",
    "output_lane_binding": "output_lane_binding",
    "child_launch_context": "child_launch_context",
    "controller_bindings": "controller_file_bindings",
}
_PREREGISTRATION_CONTEXT_PARAMETERS = frozenset(("repo_root", "verify_files"))


def _required_decisions() -> list[dict[str, Any]]:
    decisions = [
        {"name": "generation_authorized", "recommended_value": True,
         "authorized_value": None},
        {"name": "clock_ticks_per_second",
         "recommended_value": "freeze-observed-SC_CLK_TCK-before-acquisition",
         "authorized_value": None},
        {"name": "candidate_primary_cpus",
         "recommended_value": list(range(8, 16)), "authorized_value": None},
        {"name": "excluded_pairs", "recommended_value": [],
         "authorized_value": None},
        {"name": "track_b_permitted", "recommended_value": False,
         "authorized_value": None},
        {"name": "setup_invalid_budget", "recommended_value": 5,
         "authorized_value": None},
        {"name": "environment_rejected_budget", "recommended_value": 8,
         "authorized_value": None},
        {"name": "evidence_attempt_budget", "recommended_value": 3,
         "authorized_value": None},
        {"name": "scan_window_count", "recommended_value": 60,
         "authorized_value": None},
        {"name": "scan_nominal_window_ns", "recommended_value": 1_000_000_000,
         "authorized_value": None},
        {"name": "bridge_geometry", "recommended_value": {
            "minimum_window_count": 5,
            "maximum_window_count": 120,
            "nominal_window_ns": 1_000_000_000,
            "maximum_handoff_elapsed_ns": 120_000_000_000,
        }, "authorized_value": None},
        {"name": "pair_freeze_point", "recommended_value": "armed",
         "authorized_value": None},
        {"name": "candidate_executable_binding",
         "recommended_value": "freeze-sha256-before-acquisition",
         "authorized_value": None},
        {"name": "candidate_source_commit_and_tree",
         "recommended_value": "freeze-before-acquisition",
         "authorized_value": None},
        {"name": "candidate_build_closure",
         "recommended_value":
             "validated-provenance-plus-clean-reproducible-build-core",
         "authorized_value": None},
        {"name": "host_authority",
         "recommended_value": "freeze-stable-host-before-acquisition",
         "authorized_value": None},
        {"name": "output_lane_binding",
         "recommended_value":
             "freeze-canonical-path-and-linux-file-handle-before-ratification",
         "authorized_value": None},
        {"name": "child_launch_context",
         "recommended_value": recommended_launch_context_record(),
         "authorized_value": None},
        {"name": "controller_file_bindings",
         "recommended_value": "hash-reviewed-bytes-at-ratification",
         "authorized_value": None},
    ]
    parameters = set(inspect.signature(preregistration_record).parameters)
    _require(
        parameters ==
            set(_PREREGISTRATION_PARAMETER_DECISIONS) |
            _PREREGISTRATION_CONTEXT_PARAMETERS,
        "preregistration constructor decision coverage is stale")
    _require(
        {decision["name"] for decision in decisions} ==
            set(_PREREGISTRATION_PARAMETER_DECISIONS.values()),
        "preregistration decision list does not cover every policy parameter")
    return decisions


def preregistration_template_record(
    repo_root: Path = REPO_ROOT, *, verify_files: bool = True,
) -> dict[str, Any]:
    closed = closed_generation2_lineage_record(
        repo_root, verify_files=verify_files)
    return {
        "schema": TEMPLATE_SCHEMA,
        "generation": GENERATION,
        "ratification_status": "requires-explicit-user-or-external-authority",
        "safe_to_execute": False,
        "closed_generation_2": closed,
        "campaign": campaign_contract_record(),
        "no_pooling": no_pooling_record(closed),
        "reporting": reporting_contract_record(),
        "required_decisions": _required_decisions(),
    }


def validate_preregistration_template(
    value: Any, repo_root: Path = REPO_ROOT, *, verify_files: bool = True,
) -> dict[str, Any]:
    record = _exact_object(value, TEMPLATE_KEYS, "preregistration template")
    _require(record["schema"] == TEMPLATE_SCHEMA and
             record["generation"] == GENERATION and
             record["safe_to_execute"] is False,
             "preregistration template metadata differs")
    _require(type(record["required_decisions"]) is list and
             len(record["required_decisions"]) == len(_required_decisions()),
             "preregistration template decision list differs")
    for index, decision in enumerate(record["required_decisions"]):
        _exact_object(decision, DECISION_KEYS, f"required decision {index}")
        _require(decision["authorized_value"] is None,
                 "template contains an unauthorized decision")
    expected = preregistration_template_record(
        repo_root, verify_files=verify_files)
    _require(contract.exact_json_equal(record, expected),
             "preregistration template differs")
    return copy.deepcopy(expected)


def preregistration_record(
    *, authority: str, authorized_utc: str, clock_ticks_per_second: int,
    candidate_primary_cpus: Sequence[int], excluded_pairs: Sequence[Mapping[str, Any]],
    track_b_permitted: bool, setup_invalid_budget: int,
    environment_rejected_budget: int, evidence_attempt_budget: int,
    scan_window_count: int, scan_nominal_window_ns: int,
    bridge_minimum_window_count: int, bridge_maximum_window_count: int,
    bridge_nominal_window_ns: int, maximum_handoff_elapsed_ns: int,
    freeze_point: str, candidate_executable_mode: str,
    candidate_executable_sha256: str | None,
    candidate_executable_size: int,
    candidate_build_provenance_sha256: str,
    candidate_reproducible_build_core_sha256: str,
    candidate_authority_record_sha256: str,
    candidate_authority_ledger_sha256: str,
    candidate_source_commit: str, candidate_source_tree: str,
    host_machine_id_sha256: str, host_name: str, host_architecture: str,
    host_cpu_model: str, output_lane_binding: Mapping[str, Any],
    child_launch_context: Mapping[str, Any],
    controller_bindings: Sequence[Mapping[str, Any]],
    repo_root: Path = REPO_ROOT, verify_files: bool = True,
) -> dict[str, Any]:
    """Construct a final preregistration only from explicitly authorized values."""
    _require(type(authority) is str and authority.strip() == authority and
             0 < len(authority) <= 512 and
             all(0x20 <= ord(character) <= 0x7E for character in authority),
             "ratifying authority is not a bounded identity")
    _require(type(authorized_utc) is str and re.fullmatch(
        r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z",
        authorized_utc) is not None, "ratification time is not canonical UTC")
    _require(type(track_b_permitted) is bool,
             "Track-B permission is not boolean")
    _require(freeze_point == "armed", "pair freeze point is not ARMED")
    _require(candidate_executable_mode == "frozen-sha256",
             "candidate executable must be frozen before ratification")
    candidate_sha = _sha256(candidate_executable_sha256,
                            "candidate executable hash")
    candidate_size = _bounded_int(
        candidate_executable_size, 1, 64 * 1024 * 1024,
        "candidate executable size")
    build_provenance_sha = _sha256(
        candidate_build_provenance_sha256,
        "candidate build-provenance hash")
    reproducible_core_sha = _sha256(
        candidate_reproducible_build_core_sha256,
        "candidate reproducible-build core hash")
    candidate_authority_record_sha = _sha256(
        candidate_authority_record_sha256,
        "candidate authority-record file hash")
    candidate_authority_ledger_sha = _sha256(
        candidate_authority_ledger_sha256,
        "candidate authority-ledger hash")
    candidate_source = {
        "commit": _git_oid(candidate_source_commit, "candidate source commit"),
        "tree": _git_oid(candidate_source_tree, "candidate source tree"),
    }
    host_authority = _host_authority_record(
        machine_id_sha256=host_machine_id_sha256, hostname=host_name,
        architecture=host_architecture, cpu_model=host_cpu_model)
    output_lane = _validate_output_lane_binding(output_lane_binding)
    launch_context = validate_launch_context_contract(child_launch_context)
    _require(type(controller_bindings) in (list, tuple),
             "controller bindings are not a sequence")
    controllers: list[dict[str, str]] = []
    for index, raw_binding in enumerate(controller_bindings):
        binding = _exact_object(
            raw_binding, CONTROLLER_BINDING_KEYS,
            f"controller binding {index}")
        controllers.append({
            "path": _repo_relative(binding["path"],
                                   f"controller binding {index} path"),
            "sha256": _sha256(binding["sha256"],
                              f"controller binding {index} hash"),
        })
    controllers.sort(key=lambda item: item["path"])
    _require([item["path"] for item in controllers] ==
             sorted(REQUIRED_CONTROLLER_PATHS),
             "controller bindings do not cover the exact required files")
    if verify_files:
        root = repo_root.resolve(strict=True)
        for binding in controllers:
            path = root / binding["path"]
            _require(path.is_file() and not path.is_symlink() and
                     _sha256_file(path) == binding["sha256"],
                     f"controller bytes differ: {binding['path']}")
    ticks = _bounded_int(clock_ticks_per_second, 1,
                         contract.MAX_CLOCK_TICKS_PER_SECOND, "clock tick rate")
    setup_budget = _bounded_int(setup_invalid_budget, 1, 1000,
                                "setup-invalid budget")
    environment_budget = _bounded_int(environment_rejected_budget, 1, 1000,
                                      "environment-rejected budget")
    evidence_budget = _bounded_int(evidence_attempt_budget, 1, 1000,
                                   "evidence-attempt budget")
    _require(evidence_budget == 3,
             "evidence-attempt budget must preserve the frozen authority budget")
    scan_count = _bounded_int(scan_window_count, 1, contract.MAX_WINDOWS,
                              "scan window count")
    scan_nominal = _bounded_int(scan_nominal_window_ns, 1,
                                24 * 60 * 60 * 1_000_000_000,
                                "scan nominal window duration")
    bridge_geometry = {
        "minimum_window_count": bridge_minimum_window_count,
        "maximum_window_count": bridge_maximum_window_count,
        "nominal_window_ns": bridge_nominal_window_ns,
        "maximum_handoff_elapsed_ns": maximum_handoff_elapsed_ns,
    }
    try:
        bridge_geometry = bridge_contract.validate_bridge_geometry(
            bridge_geometry)
    except PreregistrationError:
        raise
    except Exception as error:
        raise PreregistrationError("bridge geometry was rejected") from error
    policy_a = contract.qualification_policy_record(
        clock_ticks_per_second=ticks,
        candidate_primary_cpus=candidate_primary_cpus,
        excluded_pairs=excluded_pairs,
        domain_mode="pair-and-domain")
    policies = [policy_a]
    policy_order = ["pair-and-domain"]
    if track_b_permitted:
        policies.append(contract.qualification_policy_record(
            clock_ticks_per_second=ticks,
            candidate_primary_cpus=candidate_primary_cpus,
            excluded_pairs=excluded_pairs,
            domain_mode="pair-only"))
        policy_order.append("pair-only")
    closed = closed_generation2_lineage_record(
        repo_root, verify_files=verify_files)
    return {
        "schema": PREREGISTRATION_SCHEMA,
        "generation": GENERATION,
        "ratification": {
            "status": "authorized",
            "authority": authority,
            "authorized_utc": authorized_utc,
        },
        "closed_generation_2": closed,
        "qualification": {
            "policy_evaluation_order": policy_order,
            "policies": policies,
            "policy_sha256s": [contract.canonical_sha256(policy)
                               for policy in policies],
            "track_b_permitted": track_b_permitted,
            "geometry": {
                "scan_window_count": scan_count,
                "scan_nominal_window_ns": scan_nominal,
                "bridge_minimum_window_count":
                    bridge_geometry["minimum_window_count"],
                "bridge_maximum_window_count":
                    bridge_geometry["maximum_window_count"],
                "bridge_nominal_window_ns":
                    bridge_geometry["nominal_window_ns"],
                "maximum_handoff_elapsed_ns":
                    bridge_geometry["maximum_handoff_elapsed_ns"],
            },
            "freeze_point": freeze_point,
            "selection_rule": contract.SELECTION_RULE,
            "frozen_pair_rule": contract.FROZEN_PAIR_RULE,
            "pre_campaign_only": True,
            "candidate_timing_performed": False,
        },
        "budgets": {
            "setup_invalid": setup_budget,
            "environment_rejected": environment_budget,
            "evidence_attempts": evidence_budget,
            "evidence_budget_commit_state": "ARMED",
            "child_launch_requires_prior_budget_commit": True,
            "crash_after_evidence_commit_consumes": "evidence_attempt",
        },
        "campaign": campaign_contract_record(),
        "candidate_executable": {
            "mode": candidate_executable_mode,
            "sha256": candidate_sha,
            "size": candidate_size,
            "build_provenance_sha256": build_provenance_sha,
            "reproducible_build_core_sha256": reproducible_core_sha,
            "authority_record_sha256": candidate_authority_record_sha,
            "authority_ledger_sha256": candidate_authority_ledger_sha,
        },
        "candidate_source": candidate_source,
        "host_authority": host_authority,
        "output_lane": output_lane,
        "launch_context": launch_context,
        "controller_bindings": controllers,
        "no_pooling": no_pooling_record(closed),
        "reporting": reporting_contract_record(),
        # Ratification freezes policy and inputs.  Only the separately hashed
        # execution contract may turn an atomically durable ARMED record into
        # child-launch authority.
        "safe_to_execute": False,
    }


def validate_preregistration(
    value: Any, repo_root: Path = REPO_ROOT, *, verify_files: bool = True,
) -> dict[str, Any]:
    record = _exact_object(value, PREREGISTRATION_KEYS, "preregistration")
    _require(record["schema"] == PREREGISTRATION_SCHEMA and
             record["generation"] == GENERATION and
             record["safe_to_execute"] is False,
             "preregistration metadata differs")
    ratification = _exact_object(
        record["ratification"], RATIFICATION_KEYS, "ratification")
    qualification = _exact_object(
        record["qualification"], QUALIFICATION_KEYS, "qualification")
    budgets = _exact_object(record["budgets"], BUDGET_KEYS, "budgets")
    geometry = _exact_object(
        qualification["geometry"], GEOMETRY_KEYS, "qualification geometry")
    candidate = _exact_object(
        record["candidate_executable"], CANDIDATE_BINDING_KEYS,
        "candidate executable binding")
    source = _exact_object(
        record["candidate_source"], SOURCE_KEYS, "candidate source")
    host = _exact_object(
        record["host_authority"], HOST_AUTHORITY_KEYS, "host authority")
    output_lane = _validate_output_lane_binding(record["output_lane"])
    launch_context = validate_launch_context_contract(
        record["launch_context"])
    controllers = record["controller_bindings"]
    _require(type(controllers) is list,
             "controller bindings are not a list")
    for index, binding in enumerate(controllers):
        _exact_object(binding, CONTROLLER_BINDING_KEYS,
                      f"controller binding {index}")
    policies = qualification.get("policies")
    _require(type(policies) is list and len(policies) in (1, 2),
             "qualification policy list differs")
    first = contract.validate_qualification_policy(policies[0])
    _require(first["domain_mode"] == "pair-and-domain",
             "Policy A is not pair-and-domain")
    if len(policies) == 2:
        second = contract.validate_qualification_policy(policies[1])
        _require(second["domain_mode"] == "pair-only" and
                 second["clock_ticks_per_second"] ==
                    first["clock_ticks_per_second"] and
                 second["candidate_primary_cpus"] ==
                    first["candidate_primary_cpus"] and
                 second["excluded_pairs"] == first["excluded_pairs"],
                 "Policy B differs from Policy A beyond domain mode")
    rebuilt = preregistration_record(
        authority=ratification["authority"],
        authorized_utc=ratification["authorized_utc"],
        clock_ticks_per_second=first["clock_ticks_per_second"],
        candidate_primary_cpus=first["candidate_primary_cpus"],
        excluded_pairs=first["excluded_pairs"],
        track_b_permitted=qualification["track_b_permitted"],
        setup_invalid_budget=budgets["setup_invalid"],
        environment_rejected_budget=budgets["environment_rejected"],
        evidence_attempt_budget=budgets["evidence_attempts"],
        scan_window_count=geometry["scan_window_count"],
        scan_nominal_window_ns=geometry["scan_nominal_window_ns"],
        bridge_minimum_window_count=geometry["bridge_minimum_window_count"],
        bridge_maximum_window_count=geometry["bridge_maximum_window_count"],
        bridge_nominal_window_ns=geometry["bridge_nominal_window_ns"],
        maximum_handoff_elapsed_ns=geometry["maximum_handoff_elapsed_ns"],
        freeze_point=qualification["freeze_point"],
        candidate_executable_mode=candidate["mode"],
        candidate_executable_sha256=candidate["sha256"],
        candidate_executable_size=candidate["size"],
        candidate_build_provenance_sha256=
            candidate["build_provenance_sha256"],
        candidate_reproducible_build_core_sha256=
            candidate["reproducible_build_core_sha256"],
        candidate_authority_record_sha256=
            candidate["authority_record_sha256"],
        candidate_authority_ledger_sha256=
            candidate["authority_ledger_sha256"],
        candidate_source_commit=source["commit"],
        candidate_source_tree=source["tree"],
        host_machine_id_sha256=host["machine_id_sha256"],
        host_name=host["hostname"],
        host_architecture=host["architecture"],
        host_cpu_model=host["cpu_model"],
        output_lane_binding=output_lane,
        child_launch_context=launch_context,
        controller_bindings=controllers,
        repo_root=repo_root, verify_files=verify_files)
    _require(contract.exact_json_equal(record, rebuilt),
             "preregistration differs from its fixed point")
    return copy.deepcopy(rebuilt)


def _validate_output_lane_path_and_lineage(
    output_lane: Path, preregistration: Mapping[str, Any], repo_root: Path,
) -> tuple[Path, dict[str, Any]]:
    _require(
        isinstance(output_lane, Path) and output_lane.is_absolute() and
        not output_lane.is_symlink(),
        "generation-3 output lane is not a canonical directory")
    try:
        record = validate_preregistration(
            preregistration, repo_root, verify_files=False)
        root = repo_root.resolve(strict=True)
        lane = output_lane.resolve(strict=True)
    except OSError as error:
        raise PreregistrationError(
            "generation-3 output lane cannot be resolved") from error
    binding = _validate_output_lane_binding(record["output_lane"])
    _require(
        output_lane == lane == Path(binding["path"]) and lane.is_dir(),
        "generation-3 output lane differs from its ratified canonical path")
    for item in record["closed_generation_2"]["attempts"]:
        prior_lexical = root.joinpath(*PurePosixPath(item["lane"]).parts)
        try:
            prior_resolved = prior_lexical.resolve(strict=False)
        except OSError as error:
            raise PreregistrationError(
                "generation-2 lane cannot be resolved") from error
        for prior in {prior_lexical, prior_resolved}:
            _require(
                lane != prior and prior not in lane.parents and
                lane not in prior.parents,
                "generation-3 output lane overlaps a generation-2 lane")
    return lane, binding


def validate_output_lane_descriptor_identity_for_preregistration(
    output_lane: Path, lane_fd: int, preregistration: Mapping[str, Any],
    repo_root: Path = REPO_ROOT,
) -> Path:
    """Join a retained directory capability to the one ratified lane inode."""
    lane, binding = _validate_output_lane_path_and_lineage(
        output_lane, preregistration, repo_root)
    try:
        retained = os.fstat(lane_fd)
        pathname = os.lstat(lane)
    except OSError as error:
        raise PreregistrationError(
            "generation-3 output lane identity cannot be revalidated") \
            from error
    retained_handle, retained_mount_id = _capture_output_lane_file_handle(
        lane_fd, b".", "retained generation-3 output lane")
    pathname_handle, pathname_mount_id = _capture_output_lane_file_handle(
        -100, os.fsencode(lane), "path generation-3 output lane")
    _require(
        stat.S_ISDIR(retained.st_mode) and
        stat.S_ISDIR(pathname.st_mode) and
        retained.st_uid == pathname.st_uid == binding["uid"] == os.geteuid() and
        stat.S_IMODE(retained.st_mode) == binding["mode"] == 0o500 and
        stat.S_IMODE(pathname.st_mode) == binding["mode"] and
        retained.st_nlink == pathname.st_nlink == binding["link_count"] and
        retained.st_mtime_ns == pathname.st_mtime_ns ==
            binding["initial_mtime_ns"] and
        retained.st_ctime_ns == pathname.st_ctime_ns ==
            binding["initial_ctime_ns"] and
        (retained.st_dev, retained.st_ino) ==
            (pathname.st_dev, pathname.st_ino) ==
            (binding["device"], binding["inode"]),
        "generation-3 output lane differs from its ratified inode")
    _require(
        retained_mount_id == pathname_mount_id and
        contract.exact_json_equal(
            retained_handle, binding["file_handle"]) and
        contract.exact_json_equal(pathname_handle, binding["file_handle"]),
        "generation-3 output lane differs from its ratified lifetime handle")
    manifest_binding = _validate_output_lane_manifest_binding(
        binding["lane_manifest"])
    manifest_fd = -1
    try:
        manifest_fd = os.open(
            manifest_binding["name"],
            os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW, dir_fd=lane_fd)
        manifest_retained = os.fstat(manifest_fd)
        manifest_pathname = os.stat(
            manifest_binding["name"], dir_fd=lane_fd,
            follow_symlinks=False)
        _require(
            stat.S_ISREG(manifest_retained.st_mode) and
            manifest_retained.st_uid == manifest_binding["uid"] and
            manifest_retained.st_nlink == manifest_binding["link_count"] and
            stat.S_IMODE(manifest_retained.st_mode) ==
                manifest_binding["mode"] and
            manifest_retained.st_size == manifest_binding["size"] and
            manifest_retained.st_mtime_ns ==
                manifest_binding["initial_mtime_ns"] and
            manifest_retained.st_ctime_ns ==
                manifest_binding["initial_ctime_ns"] and
            (manifest_retained.st_dev, manifest_retained.st_ino) ==
                (manifest_binding["device"], manifest_binding["inode"]) ==
                (manifest_pathname.st_dev, manifest_pathname.st_ino),
            "generation-3 lane manifest differs from its ratified inode")
        manifest_data = bytearray()
        while len(manifest_data) < manifest_binding["size"]:
            block = os.read(
                manifest_fd, manifest_binding["size"] - len(manifest_data))
            _require(block, "generation-3 lane manifest is truncated")
            manifest_data.extend(block)
        _require(os.read(manifest_fd, 1) == b"" and
                 hashlib.sha256(manifest_data).hexdigest() ==
                    manifest_binding["sha256"],
                 "generation-3 lane manifest bytes differ")
        retained_manifest_handle, retained_manifest_mount_id = \
            _capture_output_lane_file_handle(
                lane_fd, os.fsencode(manifest_binding["name"]),
                "retained generation-3 lane manifest")
        pathname_manifest_handle, pathname_manifest_mount_id = \
            _capture_output_lane_file_handle(
                -100, os.fsencode(lane / manifest_binding["name"]),
                "path generation-3 lane manifest")
        _require(
            retained_manifest_mount_id == pathname_manifest_mount_id ==
                retained_mount_id and
            contract.exact_json_equal(
                retained_manifest_handle, manifest_binding["file_handle"]) and
            contract.exact_json_equal(
                pathname_manifest_handle, manifest_binding["file_handle"]),
            "generation-3 lane manifest differs from its ratified lifetime "
            "handle")
    finally:
        if manifest_fd >= 0:
            os.close(manifest_fd)
    return lane


def open_output_lane_identity_for_preregistration(
    output_lane: Path, preregistration: Mapping[str, Any],
    repo_root: Path = REPO_ROOT,
) -> tuple[Path, int]:
    """Open and retain the one pre-ratified output-lane capability."""
    lane, unused_binding = _validate_output_lane_path_and_lineage(
        output_lane, preregistration, repo_root)
    descriptor = -1
    try:
        descriptor = os.open(
            lane,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW)
        validate_output_lane_descriptor_identity_for_preregistration(
            lane, descriptor, preregistration, repo_root)
        retained = descriptor
        descriptor = -1
        return lane, retained
    except OSError as error:
        raise PreregistrationError(
            "generation-3 output lane cannot be opened") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def validate_output_lane_identity_for_preregistration(
    output_lane: Path, preregistration: Mapping[str, Any],
    repo_root: Path = REPO_ROOT,
) -> Path:
    """Reject every path or inode except the one frozen at ratification."""
    lane, descriptor = open_output_lane_identity_for_preregistration(
        output_lane, preregistration, repo_root)
    os.close(descriptor)
    return lane


def validate_output_lane_for_preregistration(
    output_lane: Path, preregistration: Mapping[str, Any],
    repo_root: Path = REPO_ROOT,
) -> Path:
    """Require the exact pre-ratified topology and reject v2 overlap."""
    return validate_output_lane_identity_for_preregistration(
        output_lane, preregistration, repo_root)


def reject_denylisted_evidence_hash(
    sha256: Any, preregistration: Mapping[str, Any],
) -> None:
    observed = _sha256(sha256, "ingested evidence hash")
    registration = validate_preregistration(
        preregistration, verify_files=False)
    denylisted = {
        item["sha256"]
        for item in registration["no_pooling"]["denylisted_artifacts"]
    }
    _require(observed not in denylisted,
             "generation-2 artifact was ingested as generation-3 evidence")


def load_preregistration(
    data: bytes, repo_root: Path = REPO_ROOT, *, verify_files: bool = True,
) -> dict[str, Any]:
    validated = validate_preregistration(
        contract.strict_json_loads(data, "generation-3 preregistration"),
        repo_root, verify_files=verify_files)
    _require(data == contract.canonical_json_bytes(validated),
             "generation-3 preregistration bytes are not canonical")
    return validated


def load_preregistration_template(
    data: bytes, repo_root: Path = REPO_ROOT, *, verify_files: bool = True,
) -> dict[str, Any]:
    validated = validate_preregistration_template(
        contract.strict_json_loads(data, "generation-3 preregistration template"),
        repo_root, verify_files=verify_files)
    _require(data == contract.canonical_json_bytes(validated),
             "generation-3 preregistration template bytes are not canonical")
    return validated


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--verify", type=Path)
    action.add_argument("--verify-template", type=Path)
    options = parser.parse_args()
    path = options.verify if options.verify is not None else options.verify_template
    try:
        data = path.read_bytes()
        if options.verify is not None:
            load_preregistration(data, REPO_ROOT, verify_files=True)
            print("K65 generation-3 preregistration verified")
        else:
            # A template is permanently non-executable, so structural validation
            # is portable.  Final ratification still requires retained v2 bytes.
            load_preregistration_template(data, REPO_ROOT, verify_files=False)
            print("K65 generation-3 unratified template verified; execution forbidden")
    except (OSError, PreregistrationError, contract.QualificationError) as error:
        print(f"K65 generation-3 preregistration rejected: {error}", file=sys.stderr)
        return 1
    return 0


__all__ = (
    "CLOSED_GENERATION_SCHEMA", "GENERATION", "NO_POOLING_SCHEMA",
    "OUTPUT_LANE_FILE_HANDLE_SCHEMA", "OUTPUT_LANE_SCHEMA",
    "PREREGISTRATION_SCHEMA", "PreregistrationError",
    "REPORTING_SCHEMA",
    "TEMPLATE_SCHEMA", "TRACK_SELECTION_SCHEMA", "campaign_contract_record",
    "closed_generation2_lineage_record", "current_controller_bindings",
    "contract", "load_preregistration", "load_preregistration_template",
    "no_pooling_record", "preregistration_record", "preregistration_template_record",
    "open_output_lane_identity_for_preregistration",
    "output_lane_binding_record",
    "qualification_track_selection_record",
    "reject_denylisted_evidence_hash", "reporting_contract_record",
    "validate_closed_generation2_lineage", "validate_no_pooling",
    "validate_output_lane_identity_for_preregistration",
    "validate_output_lane_descriptor_identity_for_preregistration",
    "validate_output_lane_for_preregistration", "validate_preregistration",
    "validate_preregistration_template", "validate_qualification_track_selection",
    "validate_reporting_contract",
)


if __name__ == "__main__":
    raise SystemExit(main())
