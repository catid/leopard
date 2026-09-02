#!/usr/bin/env python3
"""Replay a sealed conditioned-passive v19 wrapper envelope.

This verifier is intentionally retained-only.  It never captures host state,
changes affinity, signals a process, builds code, or executes a benchmark.  It
reconstructs owner-only copies of the retained campaign and controller sources
and asks the frozen runner, independent auditor, passive census, and pair
qualification verifier to replay the evidence under one external attempt
authority.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, NoReturn, Sequence


SUCCESS_TERMINAL_SCHEMA = \
    "leopard2-v19-gfni-main-conditioned-passive-not-promoted-envelope/v1"
FAILURE_TERMINAL_SCHEMA = \
    "leopard2-v19-gfni-main-conditioned-passive-failed-envelope/v1"
COMPLETE_CORE_SCHEMA = \
    "leopard2-v19-gfni-main-conditioned-passive-core-manifest/v1"
FAILED_CORE_SCHEMA = \
    "leopard2-v19-gfni-main-conditioned-passive-failed-core-manifest/v1"
ATTEMPT_LINEAGE_SCHEMA = \
    "leopard2-v19-gfni-main-conditioned-attempt-lineage/v1"
FAILURE_EXIT_STATUS_SCHEMA = \
    "leopard2-v19-conditioned-campaign-exit-status/v1"
GENERATION = "passive-v3"
EVIDENCE_CLASS = \
    "conditioned-passive-windowed-shared-host-observation/v1"
PREREGISTRATION_SCHEMA = \
    "leopard2-v19-conditioned-main-wrapper-preregistration/v1"
PREREGISTRATION_SHA256 = \
    "27c1b7d76a0ecdbe194d6e6b62c01e48b1c7d10fc8ef99ebad4d76238669f0c1"
CONTROLLER_CLOSURE_SCHEMA = "leopard2-v19-controller-closure/v1"
POLICY_SHA256 = \
    "e2d605446ab3eb9b726d90ab5c62a9a4c5f4f39ce8bb91e7246061017f5f4494"
QUALIFICATION_WINDOW_COUNT = 12
QUALIFICATION_NOMINAL_WINDOW_NS = 7_000_000_000
BRIDGE_WINDOW_COUNT = 2
BRIDGE_NOMINAL_WINDOW_NS = 1_000_000_000
MAXIMUM_HANDOFF_ELAPSED_NS = 5_000_000_000
HEX40 = re.compile(r"^[0-9a-f]{40}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")

CLAIM_CEILING = {
    "promotion_eligible": False,
    "host_exclusivity_proved": False,
    "whole_campaign_interval_observed": False,
    "causal_performance_claim_allowed": False,
}

CONTROLLER_BINDINGS = (
    ("pair_qualification_acquire.py",
     "experiments/leopard2/main_compare/pair_qualification_acquire.py"),
    ("pair_qualification_bridge_acquire.py",
     "experiments/leopard2/main_compare/"
     "pair_qualification_bridge_acquire.py"),
    ("pair_qualification_bridge_contract.py",
     "experiments/leopard2/main_compare/"
     "pair_qualification_bridge_contract.py"),
    ("pair_qualification_contract.py",
     "experiments/leopard2/main_compare/pair_qualification_contract.py"),
    ("pair_qualification_verify.py",
     "experiments/leopard2/main_compare/pair_qualification_verify.py"),
    ("pair_qualified_v19_contract.py",
     "experiments/leopard2/main_compare/pair_qualified_v19_contract.py"),
)

TRUSTED_SOURCE_BINDINGS = (
    ("run-authoritative.sh",
     "experiments/leopard2/main_compare/"
     "run_authoritative_v17_gfni_main_compare.sh"),
    ("conditioned_v19_wrapper_replay.py",
     "experiments/leopard2/main_compare/conditioned_v19_wrapper_replay.py"),
    ("run_abba.py", "experiments/leopard2/main_compare/run_abba.py"),
    ("git_capture.py", "experiments/leopard2/main_compare/git_capture.py"),
    ("leopard2_build_provenance.py", "tools/leopard2_build_provenance.py"),
    ("balanced_evidence_common.py",
     "experiments/leopard2/decoder_dispatch/balanced_evidence_common.py"),
    ("audit_v17_gfni_main_compare.py",
     "experiments/leopard2/main_compare/audit_v17_gfni_main_compare.py"),
    ("passive_environment_census.py",
     "experiments/leopard2/main_compare/passive_environment_census.py"),
    *CONTROLLER_BINDINGS,
)

CONTROLLER_SOURCE_FILES = (
    "run_abba.py",
    "git_capture.py",
    "leopard2_build_provenance.py",
    "balanced_evidence_common.py",
    "audit_v17_gfni_main_compare.py",
    "passive_environment_census.py",
    *(lane for lane, _source in CONTROLLER_BINDINGS),
)

COMMON_CORE_ENTRIES = frozenset((
    "SHA256SUMS",
    "manifest.json",
    "attempt-lineage.json",
    "conditioned-v19-preregistration.json",
    "conditioned-v19-controller-closure.json",
    "conditioned_v19_wrapper_replay.py",
    "run-authoritative.sh",
    *CONTROLLER_SOURCE_FILES,
    "campaign",
))
SUCCESS_CORE_ENTRIES = COMMON_CORE_ENTRIES | frozenset((
    "pair-qualification-policy.json",
    "pair-qualification-acquisition.json",
    "pair-qualification-bridge.json",
    "pair-qualification-verdict.json",
    "controller-affinity.json",
    "environment-census-pre.json",
    "environment-census-post.json",
    "passive-environment-policy.json",
    "audit.json",
))
FAILURE_CORE_ENTRIES = COMMON_CORE_ENTRIES | frozenset((
    "campaign-exit-status.json",
))

SUCCESS_TERMINAL_KEYS = frozenset((
    "schema", "status", "acquisition_generation", "evidence_class",
    "attempt", "attempt_budget", "attempt_lineage_sha256",
    "source_commit", "source_tree", "controller_source_commit",
    "controller_source_tree", "selected_pair",
    "pair_qualification_terminal", "promotion_eligible",
    "promotion_passed", "campaign_exit_status", "core_manifest_sha256",
    "core_sha256sums_sha256", "preregistration_sha256",
    "controller_closure_sha256", "shared_host_claim_ceiling",
))
FAILURE_TERMINAL_KEYS = frozenset((
    "schema", "status", "acquisition_generation", "evidence_class",
    "attempt", "attempt_budget", "attempt_lineage_sha256",
    "source_commit", "source_tree", "controller_source_commit",
    "controller_source_tree", "selected_pair",
    "pair_qualification_stage", "pair_qualification_terminal",
    "promotion_eligible", "promotion_passed", "campaign_exit_status",
    "campaign_exit_status_sha256", "failure_sha256",
    "core_manifest_sha256", "core_sha256sums_sha256",
    "preregistration_sha256", "controller_closure_sha256",
    "shared_host_claim_ceiling",
))
COMMON_CORE_KEYS = frozenset((
    "schema", "status", "acquisition_generation", "evidence_class",
    "attempt", "attempt_budget", "attempt_lineage_sha256",
    "source_commit", "source_tree", "controller_source_commit",
    "controller_source_tree", "selected_pair",
    "preregistration_sha256", "controller_closure_sha256",
    "shared_host_claim_ceiling", "resource_envelope",
    "promotion_eligible", "promotion_passed", "campaign_exit_status",
))
COMPLETE_CORE_KEYS = COMMON_CORE_KEYS | frozenset((
    "campaign_manifest_sha256", "campaign_raw_sha256", "audit_sha256",
    "postseal_audit_sha256", "controller_affinity_sha256",
    "environment_census_pre_sha256", "environment_census_post_sha256",
    "passive_environment_policy_sha256",
    "pair_qualification_policy_sha256",
    "pair_qualification_acquisition_sha256",
    "pair_qualification_bridge_sha256",
    "pair_qualification_verdict_sha256",
    "independent_auditor_supervision_mode",
    "active_affinity_supervisor_executed",
))
FAILED_CORE_KEYS = COMMON_CORE_KEYS | frozenset((
    "campaign_failure_sha256", "failure_verify_status", "failure_verified",
    "pair_qualification_stage", "pair_qualification_terminal",
    "campaign_exit_status_sha256",
))
LINEAGE_KEYS = frozenset((
    "schema", "acquisition_generation", "attempt", "attempt_budget",
    "source_commit", "source_tree", "controller_source_commit",
    "controller_source_tree", "prior_attempts", "external_attempt",
    "v18_failure_lineage_sha256", "preregistration_sha256",
    "controller_closure_sha256",
))
PRIOR_ATTEMPT_KEYS = frozenset((
    "attempt", "envelope", "terminal", "terminal_schema",
    "envelope_sha256sums_sha256", "terminal_sha256",
))


class ReplayError(RuntimeError):
    """A retained v19 wrapper envelope failed closed."""


def fail(message: str) -> NoReturn:
    raise ReplayError(message)


def require(condition: Any, message: str) -> None:
    if not condition:
        fail(message)


def exact_json_equal(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return (set(left) == set(right) and
                all(exact_json_equal(left[key], right[key]) for key in left))
    if type(left) is list:
        return (len(left) == len(right) and
                all(exact_json_equal(a, b) for a, b in zip(left, right)))
    return left == right


def canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) +
        "\n").encode("utf-8")


def unique_object(pairs: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        require(key not in result, "wrapper replay JSON has a duplicate key")
        result[key] = value
    return result


def read_bytes(path: Path, *, maximum: int = 256 << 20) -> bytes:
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
    require(hasattr(os, "O_NOFOLLOW"), "O_NOFOLLOW is unavailable")
    flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    try:
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
                before.st_uid == os.geteuid() and
                before.st_gid == os.getegid() and
                not stat.S_IMODE(before.st_mode) & 0o222,
                f"retained file metadata is unsafe: {path}")
        require(0 <= before.st_size <= maximum,
                f"retained file exceeds replay bound: {path}")

        def read_once() -> bytes:
            chunks = []
            remaining = before.st_size
            while remaining:
                chunk = os.read(descriptor, min(1 << 20, remaining))
                require(chunk != b"", f"retained file shortened: {path}")
                chunks.append(chunk)
                remaining -= len(chunk)
            require(os.read(descriptor, 1) == b"",
                    f"retained file grew while read: {path}")
            return b"".join(chunks)

        first = read_once()
        os.lseek(descriptor, 0, os.SEEK_SET)
        second = read_once()
        after = os.fstat(descriptor)
        path_after = path.lstat()
        stable_fields = (
            "st_dev", "st_ino", "st_mode", "st_nlink", "st_uid", "st_gid",
            "st_size", "st_mtime_ns", "st_ctime_ns",
        )
        require(first == second and all(
                    getattr(before, field) == getattr(after, field) ==
                    getattr(path_after, field) for field in stable_fields),
                f"retained file changed while read: {path}")
        return first
    finally:
        os.close(descriptor)


def load_json(path: Path, *, canonical: bool = True) -> dict[str, Any]:
    data = read_bytes(path)
    value = json.loads(
        data.decode("utf-8"), object_pairs_hook=unique_object,
        parse_constant=lambda token: fail(f"non-finite JSON value: {token}"))
    require(type(value) is dict, f"retained JSON is not an object: {path}")
    if canonical:
        require(data == canonical_json_bytes(value),
                f"retained JSON is not canonical: {path}")
    return value


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(read_bytes(path))


def hex40(value: Any, label: str) -> str:
    require(type(value) is str and HEX40.fullmatch(value) is not None,
            f"{label} is not exact lowercase Git identity")
    return value


def hex64(value: Any, label: str) -> str:
    require(type(value) is str and HEX64.fullmatch(value) is not None,
            f"{label} is not exact lowercase SHA-256")
    return value


def pair(value: Any, label: str, *, nullable: bool = False) -> dict[str, int] | None:
    if nullable and value is None:
        return None
    require(type(value) is dict and set(value) == {
        "benchmark_cpu", "reserved_sibling"}, f"{label} shape differs")
    benchmark = value["benchmark_cpu"]
    sibling = value["reserved_sibling"]
    require(type(benchmark) is int and type(sibling) is int and
            0 <= benchmark <= 1_048_575 and 0 <= sibling <= 1_048_575 and
            benchmark != sibling, f"{label} value differs")
    return {"benchmark_cpu": benchmark, "reserved_sibling": sibling}


def verify_tree_safety(root: Path) -> None:
    root_metadata = root.lstat()
    require(root.is_absolute() and stat.S_ISDIR(root_metadata.st_mode) and
            not root.is_symlink() and root_metadata.st_uid == os.geteuid() and
            not stat.S_IMODE(root_metadata.st_mode) & 0o222,
            "sealed v19 replay root metadata differs")
    for directory, directory_names, filenames in os.walk(root):
        directory_path = Path(directory)
        metadata = directory_path.lstat()
        require(stat.S_ISDIR(metadata.st_mode) and
                not directory_path.is_symlink() and
                metadata.st_uid == os.geteuid() and
                not stat.S_IMODE(metadata.st_mode) & 0o222,
                f"unsafe sealed directory: {directory_path}")
        for name in (*directory_names, *filenames):
            path = directory_path / name
            item = path.lstat()
            require(not stat.S_ISLNK(item.st_mode) and
                    item.st_uid == os.geteuid() and
                    not stat.S_IMODE(item.st_mode) & 0o222,
                    f"unsafe sealed node: {path}")
            if stat.S_ISREG(item.st_mode):
                require(item.st_nlink == 1,
                        f"hard-linked sealed file: {path}")
            else:
                require(stat.S_ISDIR(item.st_mode),
                        f"special sealed node: {path}")


def verify_tree_metadata_manifest(root: Path) -> None:
    document = load_json(root / "TREE-METADATA.json")
    require(set(document) == {
                "entries", "excluded_paths", "final_mode_policy", "root",
                "schema", "self_policy", "uid_gid_policy"} and
            document["schema"] ==
                "leopard2-authoritative-tree-metadata/v1" and
            document["root"] == "." and
            document["excluded_paths"] == ["TREE-METADATA.json"] and
            document["final_mode_policy"] ==
                "observed mode with all write bits removed",
            "v19 tree metadata contract differs")
    root_stat = root.lstat()
    ownership = {
        "gid": os.getegid(),
        "rule": "every retained node has the invoking effective uid and gid",
        "uid": os.geteuid(),
    }
    self_policy = {
        "gid": os.getegid(),
        "mode": "0400",
        "nlink": 1,
        "sha256_binding":
            "exactly one ./TREE-METADATA.json checksum entry",
        "type": "file",
        "uid": os.geteuid(),
    }
    manifest_stat = (root / "TREE-METADATA.json").lstat()
    require(exact_json_equal(document["uid_gid_policy"], ownership) and
            exact_json_equal(document["self_policy"], self_policy) and
            root_stat.st_uid == os.geteuid() and
            root_stat.st_gid == os.getegid() and
            stat.S_ISREG(manifest_stat.st_mode) and
            not stat.S_ISLNK(manifest_stat.st_mode) and
            manifest_stat.st_uid == os.geteuid() and
            manifest_stat.st_gid == os.getegid() and
            stat.S_IMODE(manifest_stat.st_mode) == 0o400 and
            manifest_stat.st_nlink == 1,
            "v19 tree metadata ownership/self policy differs")
    actual = []
    for path in (root, *sorted(root.rglob("*"), key=lambda item: item.as_posix())):
        relative = "." if path == root else path.relative_to(root).as_posix()
        if relative == "TREE-METADATA.json":
            continue
        metadata = path.lstat()
        if stat.S_ISDIR(metadata.st_mode):
            kind = "directory"
        elif stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1:
            kind = "file"
        else:
            fail(f"v19 tree metadata node differs: {relative}")
        actual.append({
            "gid": metadata.st_gid,
            "mode": format(stat.S_IMODE(metadata.st_mode), "04o"),
            "nlink": metadata.st_nlink,
            "path": relative,
            "type": kind,
            "uid": metadata.st_uid,
        })
    actual.sort(key=lambda record: record["path"])
    require(exact_json_equal(document["entries"], actual),
            "v19 tree metadata has extra, missing, or changed nodes")


def regular_files(root: Path, *, omit: frozenset[str]) -> list[str]:
    result = []
    for path in root.rglob("*"):
        if path.is_file() and not path.is_symlink():
            relative = path.relative_to(root).as_posix()
            if relative not in omit:
                result.append(relative)
    return sorted(result, key=lambda value: value.encode("utf-8"))


def verify_exact_file_tree(root: Path, expected_files: Iterable[str]) -> None:
    expected = set()
    expected_directories = set()
    for value in expected_files:
        require(type(value) is str and value != "", "retained path is empty")
        relative = PurePosixPath(value)
        require(not relative.is_absolute() and
                all(part not in ("", ".", "..") for part in relative.parts),
                f"retained path is not canonical and relative: {value}")
        canonical = relative.as_posix()
        require(canonical not in expected, f"duplicate retained path: {value}")
        expected.add(canonical)
        for parent in relative.parents:
            if parent != PurePosixPath("."):
                expected_directories.add(parent.as_posix())
    actual_files = set(regular_files(root, omit=frozenset()))
    actual_directories = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*") if path.is_dir() and not path.is_symlink()
    }
    require(actual_files == expected and
            actual_directories == expected_directories,
            "retained campaign has extra, missing, or mistyped nodes")


def verify_checksum_file(
    root: Path, checksum_path: Path, expected_paths: Sequence[str],
) -> None:
    data = read_bytes(checksum_path, maximum=4 << 20)
    require(data.endswith(b"\n") and b"\r" not in data,
            f"checksum file wire format differs: {checksum_path}")
    lines = data.decode("ascii").splitlines()
    expected_lines = []
    for relative in expected_paths:
        require(relative and not relative.startswith("/") and
                ".." not in Path(relative).parts,
                "checksum path escapes its root")
        path = root / relative
        expected_lines.append(f"{sha256_file(path)}  ./{relative}")
    require(lines == expected_lines,
            f"checksum file set, order, or digest differs: {checksum_path}")


def load_local_module(name: str, path: Path) -> Any:
    specification = importlib.util.spec_from_file_location(name, path)
    require(specification is not None and specification.loader is not None,
            f"cannot load retained module: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    require(Path(module.__file__).resolve() == path.resolve(),
            f"retained module resolved outside its core: {path}")
    return module


def validate_attempt_lineage(
    lineage: dict[str, Any], *, trusted_source_root: Path,
    terminal: dict[str, Any], envelope: Path, authority_envelope: Path,
    core: Path,
    require_trusted_source_match: bool,
) -> dict[str, Any]:
    require(set(lineage) == LINEAGE_KEYS and
            lineage["schema"] == ATTEMPT_LINEAGE_SCHEMA and
            lineage["acquisition_generation"] == GENERATION,
            "v19 attempt lineage identity differs")
    attempt_number = lineage["attempt"]
    require(type(attempt_number) is int and 1 <= attempt_number <= 2 and
            lineage["attempt_budget"] == 2,
            "v19 attempt lineage budget differs")
    hex40(lineage["source_commit"], "lineage source commit")
    hex40(lineage["source_tree"], "lineage source tree")
    hex40(lineage["controller_source_commit"],
          "lineage controller source commit")
    hex40(lineage["controller_source_tree"],
          "lineage controller source tree")
    require(lineage["source_commit"] == terminal["source_commit"] and
            lineage["source_tree"] == terminal["source_tree"] and
            lineage["controller_source_commit"] ==
                terminal["controller_source_commit"] and
            lineage["controller_source_tree"] ==
                terminal["controller_source_tree"] and
            attempt_number == terminal["attempt"] and
            lineage["attempt_budget"] == terminal["attempt_budget"],
            "v19 attempt lineage differs from terminal")
    prior = lineage["prior_attempts"]
    require(type(prior) is list and len(prior) == attempt_number - 1,
            "v19 prior-attempt lineage length differs")
    prior_qualification: dict[str, Any] | None = None
    for index, entry in enumerate(prior, start=1):
        require(type(entry) is dict and set(entry) == PRIOR_ATTEMPT_KEYS and
                type(entry["attempt"]) is int and
                entry["attempt"] == index,
                "v19 prior-attempt lineage entry differs")
        expected_envelope = authority_envelope.parent / (
            f"{lineage['source_commit'][:7]}-v19-conditioned-main-a{index}")
        retained_envelope = envelope.parent / expected_envelope.name
        require(type(entry["envelope"]) is str and
                entry["envelope"] == str(expected_envelope) and
                entry["terminal"] == "FAILED.json" and
                entry["terminal_schema"] == FAILURE_TERMINAL_SCHEMA,
                "v19 prior-attempt envelope identity differs")
        hex64(entry["envelope_sha256sums_sha256"],
              "prior envelope checksum hash")
        hex64(entry["terminal_sha256"], "prior terminal hash")
        prior_terminal = load_json(retained_envelope / "FAILED.json")
        require(prior_terminal["schema"] == FAILURE_TERMINAL_SCHEMA and
                prior_terminal["status"] == "failed" and
                prior_terminal["acquisition_generation"] == GENERATION and
                type(prior_terminal["attempt"]) is int and
                prior_terminal["attempt"] == index and
                prior_terminal["source_commit"] == lineage["source_commit"] and
                prior_terminal["source_tree"] == lineage["source_tree"] and
                prior_terminal["controller_source_commit"] ==
                    lineage["controller_source_commit"] and
                prior_terminal["controller_source_tree"] ==
                    lineage["controller_source_tree"],
                "v19 prior-attempt terminal lineage differs")
        for core_name, _relative in TRUSTED_SOURCE_BINDINGS:
            require(read_bytes(core / core_name) ==
                    read_bytes(retained_envelope / "core" / core_name),
                    "v19 prior/current retained source closure differs")
        prior_result = verify_envelope(
            retained_envelope, trusted_source_root,
            require_trusted_source_match=require_trusted_source_match,
            authority_envelope=expected_envelope)
        require(prior_result["attempt"] == index and
                prior_result["status"] == "failed" and
                sha256_file(retained_envelope / "SHA256SUMS") ==
                    entry["envelope_sha256sums_sha256"] and
                sha256_file(retained_envelope / "FAILED.json") ==
                    entry["terminal_sha256"],
                "v19 prior-attempt retained envelope binding differs")
        prior_failure = load_json(
            retained_envelope / "core/campaign/failure.json")
        prior_qualification = prior_failure.get("pair_qualification")
        require(type(prior_qualification) is dict,
                "v19 prior failure lacks pair qualification")
    v19 = load_local_module(
        "leopard2_conditioned_wrapper_pair_v19",
        trusted_source_root /
        "experiments/leopard2/main_compare/"
        "pair_qualified_v19_contract.py")
    try:
        authority = v19.validate_pair_qualified_attempt(
            lineage["external_attempt"])
    except v19.PairQualifiedV19Error as error:
        raise ReplayError(f"external v19 attempt authority is invalid: {error}") \
            from error
    require(authority["attempt"] == attempt_number and
            authority["attempt_budget"] == 2,
            "external v19 attempt authority differs from lineage")
    if attempt_number == 1:
        require(prior == [], "attempt one retained prior attempts")
    else:
        require(authority["prior_attempt_failure_sha256"] ==
                prior[-1]["terminal_sha256"],
                "prior failure authority differs from lineage terminal")
        require(prior_qualification is not None,
                "attempt two lacks verified prior qualification")
        prior_acquisition = prior_qualification.get("acquisition")
        prior_selected = prior_qualification.get("selected_pair")
        if prior_acquisition is None:
            prior_status = "not-acquired"
            prior_acquisition_hash = None
        else:
            prior_acquisition_hash = prior_qualification.get(
                "acquisition_sha256")
            prior_status = (
                "selected-lowest-primary" if prior_selected is not None else
                "no-candidate-pair-qualified")
        require(authority["prior_attempt_acquisition_sha256"] ==
                    prior_acquisition_hash and
                authority["prior_attempt_selection_status"] == prior_status and
                exact_json_equal(authority["prior_attempt_selected_pair"],
                                 prior_selected),
                "prior acquisition/selection authority differs from envelope")
    hex64(lineage["v18_failure_lineage_sha256"],
          "v18 failure lineage hash")
    expected_v18_lineage_sha256 = sha256_bytes(
        canonical_json_bytes(v19.v18_failure_lineage_record()))
    # This joins the trusted controller's content-addressed historical record.
    # Physical v18 archive availability is deliberately covered by the
    # separate exact-v18 replay gate rather than made a v19 replay dependency.
    require(lineage["v18_failure_lineage_sha256"] ==
            expected_v18_lineage_sha256,
            "v18 failure lineage differs from trusted controller constants")
    require(lineage["preregistration_sha256"] == PREREGISTRATION_SHA256,
            "lineage preregistration hash differs")
    hex64(lineage["controller_closure_sha256"],
          "lineage controller closure hash")
    return authority


def attempt_cli_arguments(authority: dict[str, Any]) -> list[str]:
    result = ["--v19-attempt", str(authority["attempt"])]
    mappings = (
        ("prior_attempt_failure_sha256", "--v19-prior-failure-sha256"),
        ("prior_attempt_acquisition_sha256",
         "--v19-prior-acquisition-sha256"),
        ("prior_attempt_selection_status", "--v19-prior-selection-status"),
    )
    for field, option in mappings:
        value = authority[field]
        if value is not None:
            result.extend((option, str(value)))
    selected = authority["prior_attempt_selected_pair"]
    if selected is not None:
        result.extend((
            "--v19-prior-selected-pair",
            f"{selected['benchmark_cpu']},{selected['reserved_sibling']}",
        ))
    return result


def validate_controller_closure(core: Path) -> tuple[dict[str, Any], str]:
    record_path = core / "conditioned-v19-controller-closure.json"
    record = load_json(record_path)
    require(set(record) == {"schema", "files"} and
            record["schema"] == CONTROLLER_CLOSURE_SCHEMA and
            type(record["files"]) is list and
            len(record["files"]) == len(CONTROLLER_BINDINGS),
            "conditioned v19 controller closure shape differs")
    expected = []
    for lane_path, source_path in CONTROLLER_BINDINGS:
        retained = core / lane_path
        expected.append({
            "lane_path": lane_path,
            "source_path": source_path,
            "sha256": sha256_file(retained),
        })
    require(exact_json_equal(record["files"], expected),
            "conditioned v19 controller closure differs")
    return record, sha256_file(record_path)


def validate_trusted_source_closure(
    core: Path, trusted_source_root: Path,
) -> None:
    trusted_source_root = trusted_source_root.resolve(strict=True)
    trusted_metadata = trusted_source_root.lstat()
    require(stat.S_ISDIR(trusted_metadata.st_mode) and
            not trusted_source_root.is_symlink() and
            trusted_metadata.st_uid == os.geteuid() and
            trusted_metadata.st_gid == os.getegid() and
            not stat.S_IMODE(trusted_metadata.st_mode) & 0o222,
            "trusted v19 replay source root is unsafe")
    for core_name, relative in TRUSTED_SOURCE_BINDINGS:
        trusted_path = trusted_source_root / relative
        require(read_bytes(core / core_name) == read_bytes(trusted_path),
                f"retained v19 source differs from trusted source: {core_name}")


def validate_preflight_identity(
    raw: dict[str, Any], campaign_manifest: dict[str, Any],
    preregistration: dict[str, Any],
) -> None:
    preflight = preregistration.get("build_preflight")
    require(type(preflight) is dict,
            "conditioned v19 build preflight is absent")
    expected = {
        "candidate_source": {
            "head": preflight.get("source_commit"),
            "tree": preflight.get("source_tree"),
        },
        "candidate_executable": preflight.get("candidate_binary_sha256"),
        "candidate_archive": preflight.get("candidate_archive_sha256"),
        "baseline_source": {
            "head": preflight.get("baseline_commit"),
            "tree": preflight.get("baseline_tree"),
        },
        "baseline_executable": preflight.get("baseline_binary_sha256"),
        "baseline_archive": preflight.get("baseline_archive_sha256"),
    }
    for value in expected.values():
        if type(value) is dict:
            hex40(value.get("head"), "preflight source commit")
            hex40(value.get("tree"), "preflight source tree")
        else:
            hex64(value, "preflight artifact hash")

    identities = [
        raw.get("identities_initial"), raw.get("identities_final"),
        campaign_manifest.get("identities"),
    ]
    invocations = raw.get("invocations")
    require(type(invocations) is list,
            "v19 preflight identity invocation set differs")
    for invocation in invocations:
        require(type(invocation) is dict,
                "v19 preflight identity invocation differs")
        identities.extend((
            invocation.get("identity_before"),
            invocation.get("identity_after"),
        ))
    for identity in identities:
        require(type(identity) is dict,
                "v19 retained build identity is absent")
        for role in ("candidate", "baseline"):
            source = identity.get(f"{role}_source")
            executable = identity.get(f"{role}_executable")
            archive = identity.get(f"{role}_archive")
            require(type(source) is dict and
                    source.get("head") ==
                        expected[f"{role}_source"]["head"] and
                    source.get("tree") ==
                        expected[f"{role}_source"]["tree"] and
                    type(executable) is dict and
                    executable.get("sha256") ==
                        expected[f"{role}_executable"] and
                    type(archive) is dict and
                    archive.get("sha256") == expected[f"{role}_archive"],
                    f"v19 {role} preflight identity differs")


def validate_failure_preflight_identity(
    failure: dict[str, Any], preregistration: dict[str, Any],
) -> None:
    initial = failure.get("identities_initial")
    invocations = failure.get("invocations")
    require(type(initial) is dict and type(invocations) is list,
            "v19 failure lacks its preflight identity boundary")
    # A failure has no final campaign manifest/identity.  Reuse the exact
    # success identity validator with the retained initial identity as both
    # endpoints, while still validating every retained invocation endpoint.
    validate_preflight_identity({
        "identities_initial": initial,
        "identities_final": initial,
        "invocations": invocations,
    }, {"identities": initial}, preregistration)


def validate_common_records(
    envelope: Path, terminal_path: Path, *, success: bool,
    trusted_source_root: Path, authority_envelope: Path,
    require_trusted_source_match: bool,
) -> tuple[
    Path, dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]
]:
    core = envelope / "core"
    require(set(path.name for path in envelope.iterdir()) ==
            ({"NOT_PROMOTED.json", "postseal-audit.json", "SHA256SUMS",
              "TREE-METADATA.json", "core"} if success else
             {"FAILED.json", "SHA256SUMS", "TREE-METADATA.json", "core"}),
            "sealed v19 envelope file set differs")
    expected_core_entries = SUCCESS_CORE_ENTRIES if success else \
        FAILURE_CORE_ENTRIES
    require(set(path.name for path in core.iterdir()) == expected_core_entries,
            "sealed v19 core file set differs")
    core_checksum_paths = regular_files(
        core, omit=frozenset(("SHA256SUMS",)))
    verify_checksum_file(
        core, core / "SHA256SUMS", core_checksum_paths)
    outer_checksum_paths = [
        "NOT_PROMOTED.json" if success else "FAILED.json",
    ]
    if success:
        outer_checksum_paths.append("postseal-audit.json")
    outer_checksum_paths.extend(("core/SHA256SUMS", "TREE-METADATA.json"))
    verify_checksum_file(
        envelope, envelope / "SHA256SUMS", outer_checksum_paths)

    terminal = load_json(terminal_path)
    manifest = load_json(core / "manifest.json")
    lineage_path = core / "attempt-lineage.json"
    lineage = load_json(lineage_path)
    if require_trusted_source_match:
        validate_trusted_source_closure(core, trusted_source_root)
    expected_terminal_keys = (
        SUCCESS_TERMINAL_KEYS if success else FAILURE_TERMINAL_KEYS)
    expected_manifest_keys = COMPLETE_CORE_KEYS if success else FAILED_CORE_KEYS
    require(set(terminal) == expected_terminal_keys and
            set(manifest) == expected_manifest_keys,
            "v19 terminal or core manifest key set differs")
    require(terminal["schema"] == (
                SUCCESS_TERMINAL_SCHEMA if success else
                FAILURE_TERMINAL_SCHEMA) and
            manifest["schema"] == (
                COMPLETE_CORE_SCHEMA if success else FAILED_CORE_SCHEMA) and
            terminal["status"] == ("complete" if success else "failed") and
            manifest["status"] == terminal["status"] and
            terminal["acquisition_generation"] == GENERATION and
            manifest["acquisition_generation"] == GENERATION and
            terminal["evidence_class"] == EVIDENCE_CLASS and
            manifest["evidence_class"] == EVIDENCE_CLASS and
            terminal["attempt_budget"] == 2 and
            manifest["attempt_budget"] == 2 and
            terminal["attempt"] == manifest["attempt"] and
            terminal["attempt_lineage_sha256"] ==
                manifest["attempt_lineage_sha256"] == sha256_file(lineage_path),
            "v19 terminal/core identity differs")
    for record in (terminal, manifest):
        hex40(record["source_commit"], "v19 source commit")
        hex40(record["source_tree"], "v19 source tree")
        hex40(record["controller_source_commit"],
              "v19 controller source commit")
        hex40(record["controller_source_tree"],
              "v19 controller source tree")
        require(record["promotion_eligible"] is False and
                record["promotion_passed"] is False and
                exact_json_equal(record["shared_host_claim_ceiling"],
                                 CLAIM_CEILING),
                "v19 promotion or claim ceiling differs")
    require(terminal["source_commit"] == manifest["source_commit"] and
            terminal["source_tree"] == manifest["source_tree"] and
            terminal["controller_source_commit"] ==
                manifest["controller_source_commit"] and
            terminal["controller_source_tree"] ==
                manifest["controller_source_tree"] and
            terminal["core_manifest_sha256"] ==
                sha256_file(core / "manifest.json") and
            terminal["core_sha256sums_sha256"] ==
                sha256_file(core / "SHA256SUMS"),
            "v19 terminal core binding differs")

    preregistration_path = core / "conditioned-v19-preregistration.json"
    preregistration = load_json(preregistration_path)
    require(preregistration.get("schema") == PREREGISTRATION_SCHEMA and
            sha256_file(preregistration_path) == PREREGISTRATION_SHA256 and
            terminal["preregistration_sha256"] == PREREGISTRATION_SHA256 and
            manifest["preregistration_sha256"] == PREREGISTRATION_SHA256 and
            exact_json_equal(manifest["resource_envelope"],
                             preregistration.get("resource_envelope")),
            "v19 preregistration/resource binding differs")
    preflight = preregistration.get("build_preflight")
    require(type(preflight) is dict and
            terminal["source_commit"] == manifest["source_commit"] ==
                lineage["source_commit"] == preflight.get("source_commit") and
            terminal["source_tree"] == manifest["source_tree"] ==
                lineage["source_tree"] == preflight.get("source_tree"),
            "v19 candidate source authority differs from preregistration")
    require(preregistration.get("attempt_contract", {}).get(
                "failure_lineage_sha256") ==
            lineage["v18_failure_lineage_sha256"],
            "v19 preregistration/v18 lineage binding differs")
    _closure, closure_hash = validate_controller_closure(core)
    require(terminal["controller_closure_sha256"] == closure_hash and
            manifest["controller_closure_sha256"] == closure_hash and
            lineage["controller_closure_sha256"] == closure_hash,
            "v19 controller closure hash binding differs")
    authority = validate_attempt_lineage(
        lineage, trusted_source_root=trusted_source_root, terminal=terminal,
        envelope=envelope, authority_envelope=authority_envelope, core=core,
        require_trusted_source_match=require_trusted_source_match)
    require(lineage["source_commit"] == manifest["source_commit"] and
            lineage["source_tree"] == manifest["source_tree"] and
            lineage["controller_source_commit"] ==
                manifest["controller_source_commit"] and
            lineage["controller_source_tree"] ==
                manifest["controller_source_tree"] and
            lineage["preregistration_sha256"] == PREREGISTRATION_SHA256,
            "v19 lineage source/preregistration binding differs")
    return core, terminal, manifest, lineage, authority


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    data = read_bytes(source)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(destination, flags, 0o600)
    try:
        os.fchmod(descriptor, 0o600)
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            require(written > 0, f"short replay copy: {destination}")
            view = view[written:]
    finally:
        os.close(descriptor)


def copy_tree_owner_only(source: Path, destination: Path) -> None:
    require(not destination.exists(), "replay destination already exists")
    destination.mkdir(mode=0o700)
    for root, directory_names, filenames in os.walk(source):
        root_path = Path(root)
        relative = root_path.relative_to(source)
        target_root = destination / relative
        target_root.chmod(0o700)
        for name in sorted(directory_names):
            target = target_root / name
            target.mkdir(mode=0o700)
            target.chmod(0o700)
        for name in sorted(filenames):
            copy_file(root_path / name, target_root / name)


def controller_tree(trusted_source_root: Path, root: Path) -> Path:
    main_compare = root / "experiments/leopard2/main_compare"
    decoder = root / "experiments/leopard2/decoder_dispatch"
    tools = root / "tools"
    main_compare.mkdir(mode=0o700, parents=True)
    decoder.mkdir(mode=0o700, parents=True)
    tools.mkdir(mode=0o700, parents=True)
    for filename in CONTROLLER_SOURCE_FILES:
        if filename == "leopard2_build_provenance.py":
            destination = tools / filename
            source = trusted_source_root / "tools" / filename
        elif filename == "balanced_evidence_common.py":
            destination = decoder / filename
            source = trusted_source_root / \
                "experiments/leopard2/decoder_dispatch" / filename
        else:
            destination = main_compare / filename
            source = trusted_source_root / \
                "experiments/leopard2/main_compare" / filename
        copy_file(source, destination)
    return main_compare


def python_command(script: Path, *arguments: str) -> list[str]:
    optimization = ["-O"] if not __debug__ else []
    return ["/usr/bin/python3", "-I", "-S", "-B", *optimization,
            str(script), *arguments]


def run_checked(arguments: Sequence[str]) -> subprocess.CompletedProcess[bytes]:
    completed = subprocess.run(
        list(arguments), stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False, timeout=120)
    require(completed.returncode == 0,
            "retained consumer replay failed: " +
            completed.stderr.decode("utf-8", errors="replace"))
    return completed


def validate_success(
    envelope: Path, core: Path, terminal: dict[str, Any],
    manifest: dict[str, Any], lineage: dict[str, Any],
    authority: dict[str, Any], trusted_source_root: Path,
    require_preflight_identity: bool,
) -> dict[str, Any]:
    selected = pair(terminal["selected_pair"], "success selected pair")
    require(selected is not None and
            exact_json_equal(manifest["selected_pair"], selected) and
            terminal["pair_qualification_terminal"] == "NOT_PROMOTED" and
            type(terminal["campaign_exit_status"]) is int and
            terminal["campaign_exit_status"] == 0 and
            type(manifest["campaign_exit_status"]) is int and
            manifest["campaign_exit_status"] == 0 and
            manifest["independent_auditor_supervision_mode"] ==
                "conditioned" and
            manifest["active_affinity_supervisor_executed"] is False,
            "v19 success terminal/core semantics differ")
    bindings = {
        "campaign_manifest_sha256": core / "campaign/manifest.json",
        "campaign_raw_sha256": core / "campaign/raw.json",
        "audit_sha256": core / "audit.json",
        "postseal_audit_sha256": envelope / "postseal-audit.json",
        "controller_affinity_sha256": core / "controller-affinity.json",
        "environment_census_pre_sha256":
            core / "environment-census-pre.json",
        "environment_census_post_sha256":
            core / "environment-census-post.json",
        "passive_environment_policy_sha256":
            core / "passive-environment-policy.json",
        "pair_qualification_policy_sha256":
            core / "pair-qualification-policy.json",
        "pair_qualification_acquisition_sha256":
            core / "pair-qualification-acquisition.json",
        "pair_qualification_bridge_sha256":
            core / "pair-qualification-bridge.json",
        "pair_qualification_verdict_sha256":
            core / "pair-qualification-verdict.json",
    }
    for field, path in bindings.items():
        require(manifest[field] == sha256_file(path),
                f"v19 core hash binding differs: {field}")
    require(read_bytes(core / "audit.json") ==
            read_bytes(envelope / "postseal-audit.json"),
            "v19 postseal audit is not byte-identical")

    raw = load_json(core / "campaign/raw.json")
    expected_campaign_files = ["manifest.json", "raw.json"]
    invocations = raw.get("invocations")
    require(type(invocations) is list,
            "v19 campaign invocation inventory differs")
    for invocation in invocations:
        require(type(invocation) is dict,
                "v19 campaign invocation record differs")
        for stream_name in ("stdout", "stderr"):
            stream = invocation.get(stream_name)
            require(type(stream) is dict and type(stream.get("path")) is str,
                    "v19 campaign stream inventory differs")
            expected_campaign_files.append(stream["path"])
    verify_exact_file_tree(core / "campaign", expected_campaign_files)
    campaign_manifest = load_json(core / "campaign/manifest.json")
    if require_preflight_identity:
        validate_preflight_identity(
            raw, campaign_manifest,
            load_json(core / "conditioned-v19-preregistration.json"))
    controller = load_json(core / "controller-affinity.json")
    pre = load_json(core / "environment-census-pre.json")
    post = load_json(core / "environment-census-post.json")
    audit = load_json(core / "audit.json")
    policy = load_json(core / "passive-environment-policy.json")
    qualification = raw.get("pair_qualification")
    require(type(qualification) is dict and
            qualification.get("record_status") == "complete" and
            qualification.get("stage") == "complete" and
            qualification.get("terminal") == "NOT_PROMOTED" and
            exact_json_equal(qualification.get("attempt"), authority) and
            exact_json_equal(qualification.get("selected_pair"), selected),
            "v19 complete qualification differs from authority/pair")
    require(qualification.get("v18_failure_lineage_sha256") ==
            lineage["v18_failure_lineage_sha256"],
            "v19 wrapper/runner v18 lineage hash differs")
    require(exact_json_equal(campaign_manifest.get("pair_qualification"),
                             qualification),
            "v19 campaign manifest qualification differs")
    require(raw.get("supervision") is None and
            raw.get("campaign", {}).get("benchmark_cpu") ==
                selected["benchmark_cpu"] and
            raw.get("campaign", {}).get("reserved_sibling") ==
                selected["reserved_sibling"],
            "v19 campaign dynamic pair/supervision differs")
    require(controller.get("schema") ==
                "leopard2-v19-passive-controller-affinity/v1" and
            controller.get("acquisition_generation") == GENERATION and
            controller.get("benchmark_cpu") == selected["benchmark_cpu"] and
            controller.get("reserved_sibling") ==
                selected["reserved_sibling"] and
            controller.get("active_affinity_supervisor_executed") is False and
            controller.get("verified_acquisition_sha256") ==
                qualification.get("acquisition_sha256") and
            controller.get("verified_bridge_sha256") ==
                qualification.get("bridge_sha256"),
            "v19 controller pair/hash binding differs")
    require(exact_json_equal(audit.get("pair_qualification"),
                             policy.get("pair_qualification")) and
            audit.get("schema") ==
                "leopard2-main-compare-v19-conditioned-passive-"
                "independent-audit/v1" and
            audit.get("acquisition_generation") == GENERATION and
            audit.get("evidence_class") == EVIDENCE_CLASS and
            policy.get("schema") ==
                "leopard2-passive-shared-host-policy/v3" and
            policy.get("acquisition_generation") == GENERATION and
            policy.get("evidence_class") == EVIDENCE_CLASS,
            "v19 retained audit/census identity differs")
    for value in (audit, policy):
        require(value.get("promotion_eligible") is False and
                value.get("promotion_passed") is False and
                value.get("cpu_pair_exclusive") is False and
                value.get("causal_performance_claim_eligible") is False,
                "v19 retained consumer overclaims evidence")

    pair_policy = load_json(core / "pair-qualification-policy.json")
    acquisition = load_json(core / "pair-qualification-acquisition.json")
    bridge = load_json(core / "pair-qualification-bridge.json")
    require(sha256_bytes(canonical_json_bytes(pair_policy)) == POLICY_SHA256 and
            exact_json_equal(qualification.get("policy"), pair_policy) and
            exact_json_equal(qualification.get("acquisition"), acquisition) and
            exact_json_equal(qualification.get("bridge"), bridge),
            "v19 standalone pair records differ from campaign evidence")

    attempt_arguments = attempt_cli_arguments(authority)
    with tempfile.TemporaryDirectory(
            prefix="leopard-v19-wrapper-consumers-") as temporary:
        root = Path(temporary)
        retained_campaign = root / "campaign"
        copy_tree_owner_only(core / "campaign", retained_campaign)
        main_compare = controller_tree(
            trusted_source_root, root / "controller")

        runner_result = run_checked(python_command(
            main_compare / "run_abba.py", "verify", "--manifest",
            str(retained_campaign / "manifest.json"),
            "--no-current-input-check", *attempt_arguments))
        require(runner_result.stderr == b"",
                "v19 runner replay wrote stderr")

        verdict = run_checked(python_command(
            main_compare / "pair_qualification_verify.py", "verify",
            "--policy", str(core / "pair-qualification-policy.json"),
            "--acquisition",
            str(core / "pair-qualification-acquisition.json"),
            "--bridge", str(core / "pair-qualification-bridge.json"),
            "--expected-policy-sha256", POLICY_SHA256,
            "--expected-acquisition-window-count",
            str(QUALIFICATION_WINDOW_COUNT),
            "--expected-acquisition-nominal-window-ns",
            str(QUALIFICATION_NOMINAL_WINDOW_NS),
            "--minimum-bridge-window-count", str(BRIDGE_WINDOW_COUNT),
            "--maximum-bridge-window-count", str(BRIDGE_WINDOW_COUNT),
            "--bridge-nominal-window-ns", str(BRIDGE_NOMINAL_WINDOW_NS),
            "--maximum-bridge-handoff-elapsed-ns",
            str(MAXIMUM_HANDOFF_ELAPSED_NS),
            *(() if authority["frozen_pair_from_prior_attempt"] is None else (
                "--frozen-benchmark-cpu",
                str(authority["frozen_pair_from_prior_attempt"]
                    ["benchmark_cpu"]),
                "--frozen-reserved-sibling",
                str(authority["frozen_pair_from_prior_attempt"]
                    ["reserved_sibling"]),
            ))))
        require(verdict.stdout ==
                read_bytes(core / "pair-qualification-verdict.json") and
                verdict.stderr == b"",
                "independent v19 pair verdict differs")

        audit_output = root / "audit.json"
        audit_result = run_checked(python_command(
            main_compare / "audit_v17_gfni_main_compare.py",
            "--manifest", str(retained_campaign / "manifest.json"),
            "--output", str(audit_output), "--supervision", "conditioned",
            *attempt_arguments))
        require(audit_result.stderr == b"" and
                audit_output.read_bytes() == read_bytes(core / "audit.json"),
                "independent conditioned audit replay differs")

        reserved = \
            f"{selected['benchmark_cpu']},{selected['reserved_sibling']}"
        for phase, path in (("pre", core / "environment-census-pre.json"),
                            ("post", core / "environment-census-post.json")):
            verified = run_checked(python_command(
                main_compare / "passive_environment_census.py", "verify",
                "--input", str(path), "--phase", phase,
                "--generation", GENERATION, "--reserved-cpus", reserved))
            require(verified.stdout == b"" and verified.stderr == b"",
                    "passive-v3 census snapshot replay wrote output")
        policy_output = root / "policy.json"
        compared = run_checked(python_command(
            main_compare / "passive_environment_census.py", "compare",
            "--pre", str(core / "environment-census-pre.json"),
            "--post", str(core / "environment-census-post.json"),
            "--raw", str(retained_campaign / "raw.json"),
            "--controller", str(core / "controller-affinity.json"),
            "--output", str(policy_output), *attempt_arguments))
        require(compared.stdout == b"" and compared.stderr == b"" and
                policy_output.read_bytes() ==
                read_bytes(core / "passive-environment-policy.json"),
                "passive-v3 census comparison differs")
    return {
        "attempt": authority["attempt"],
        "campaign_exit_status": terminal["campaign_exit_status"],
        "controller_source_commit": terminal["controller_source_commit"],
        "controller_source_tree": terminal["controller_source_tree"],
        "outer_sha256sums_sha256": sha256_file(envelope / "SHA256SUMS"),
        "selected_pair": selected,
        "source_commit": terminal["source_commit"],
        "source_tree": terminal["source_tree"],
        "status": "complete",
        "terminal": "NOT_PROMOTED",
        "terminal_file": "NOT_PROMOTED.json",
        "terminal_sha256": sha256_file(envelope / "NOT_PROMOTED.json"),
    }


def validate_failure(
    core: Path, terminal: dict[str, Any], manifest: dict[str, Any],
    lineage: dict[str, Any], authority: dict[str, Any],
    trusted_source_root: Path, require_preflight_identity: bool,
) -> dict[str, Any]:
    terminal_pair = pair(
        terminal["selected_pair"], "failure terminal selected pair",
        nullable=True)
    manifest_pair = pair(
        manifest["selected_pair"], "failure manifest selected pair",
        nullable=True)
    exit_status_path = core / "campaign-exit-status.json"
    # run_abba's failure record does not contain its caller's shell status.
    # The wrapper therefore retains this separately produced, core-hashed
    # record as the replay source of truth for campaign_exit_status.
    exit_status_record = load_json(exit_status_path)
    require(set(exit_status_record) == {"campaign_exit_status", "schema"} and
            exit_status_record["schema"] == FAILURE_EXIT_STATUS_SCHEMA and
            type(exit_status_record["campaign_exit_status"]) is int and
            1 <= exit_status_record["campaign_exit_status"] <= 255,
            "v19 retained campaign exit status differs")
    retained_exit_status = exit_status_record["campaign_exit_status"]
    exit_status_sha256 = sha256_file(exit_status_path)
    require(exact_json_equal(terminal_pair, manifest_pair) and
            type(terminal["campaign_exit_status"]) is int and
            terminal["campaign_exit_status"] == retained_exit_status and
            type(manifest["campaign_exit_status"]) is int and
            manifest["campaign_exit_status"] ==
                terminal["campaign_exit_status"] and
            terminal["campaign_exit_status_sha256"] ==
                manifest["campaign_exit_status_sha256"] ==
                exit_status_sha256 and
            type(manifest["failure_verify_status"]) is int and
            manifest["failure_verify_status"] == 0 and
            manifest["failure_verified"] is True,
            "v19 failed terminal/core status differs")
    failure_path = core / "campaign/failure.json"
    require(terminal["failure_sha256"] ==
            manifest["campaign_failure_sha256"] == sha256_file(failure_path),
            "v19 failure hash binding differs")
    failure = load_json(failure_path)
    if require_preflight_identity:
        validate_failure_preflight_identity(
            failure,
            load_json(core / "conditioned-v19-preregistration.json"))
    retained = failure.get("retained_files")
    require(type(retained) is list,
            "v19 failure retained-file inventory differs")
    expected_campaign_files = ["failure.json"]
    for record in retained:
        require(type(record) is dict and type(record.get("path")) is str,
                "v19 failure retained-file record differs")
        expected_campaign_files.append(record["path"])
    verify_exact_file_tree(core / "campaign", expected_campaign_files)
    qualification = failure.get("pair_qualification")
    qualification_pair = pair(
        qualification.get("selected_pair") if type(qualification) is dict
        else None,
        "failure qualification selected pair", nullable=True)
    require(type(qualification) is dict and
            exact_json_equal(qualification.get("attempt"), authority) and
            exact_json_equal(qualification_pair, terminal_pair) and
            qualification.get("stage") ==
                terminal["pair_qualification_stage"] ==
                manifest["pair_qualification_stage"] and
            qualification.get("terminal") ==
                terminal["pair_qualification_terminal"] ==
                manifest["pair_qualification_terminal"],
            "v19 failure qualification/terminal binding differs")
    require(qualification.get("v18_failure_lineage_sha256") ==
            lineage["v18_failure_lineage_sha256"],
            "v19 failed wrapper/runner v18 lineage hash differs")
    attempt_arguments = attempt_cli_arguments(authority)
    with tempfile.TemporaryDirectory(
            prefix="leopard-v19-wrapper-failure-") as temporary:
        root = Path(temporary)
        campaign = root / "campaign"
        copy_tree_owner_only(core / "campaign", campaign)
        main_compare = controller_tree(
            trusted_source_root, root / "controller")
        replayed = run_checked(python_command(
            main_compare / "run_abba.py", "verify-failure", "--failure",
            str(campaign / "failure.json"), *attempt_arguments))
        require(replayed.stderr == b"",
                "v19 failure runner replay wrote stderr")
    return {
        "attempt": authority["attempt"],
        "campaign_exit_status": terminal["campaign_exit_status"],
        "controller_source_commit": terminal["controller_source_commit"],
        "controller_source_tree": terminal["controller_source_tree"],
        "outer_sha256sums_sha256": sha256_file(core.parent / "SHA256SUMS"),
        "selected_pair": terminal_pair,
        "source_commit": terminal["source_commit"],
        "source_tree": terminal["source_tree"],
        "status": "failed",
        "terminal": qualification["terminal"],
        "terminal_file": "FAILED.json",
        "terminal_sha256": sha256_file(core.parent / "FAILED.json"),
    }


def verify_envelope(
    envelope: Path, trusted_source_root: Path, *,
    require_trusted_source_match: bool = False,
    require_preflight_identity: bool = False,
    authority_envelope: Path | None = None,
) -> dict[str, Any]:
    require(envelope.is_absolute(), "v19 envelope path is not absolute")
    canonical = Path(os.path.abspath(os.fspath(envelope)))
    require(canonical == envelope and not stat.S_ISLNK(envelope.lstat().st_mode),
            "v19 envelope path is not canonical or is a symbolic link")
    resolved = envelope.resolve(strict=True)
    require(resolved == envelope,
            "v19 envelope path resolves to another location")
    if authority_envelope is None:
        authority_envelope = envelope
    require(authority_envelope.is_absolute() and
            Path(os.path.abspath(os.fspath(authority_envelope))) ==
                authority_envelope and
            authority_envelope.name == envelope.name,
            "v19 authority envelope label differs")
    verify_tree_safety(envelope)
    verify_tree_metadata_manifest(envelope)
    success_path = envelope / "NOT_PROMOTED.json"
    failure_path = envelope / "FAILED.json"
    require(success_path.exists() != failure_path.exists(),
            "v19 envelope must retain exactly one terminal")
    success = success_path.exists()
    terminal_path = success_path if success else failure_path
    core, terminal, manifest, lineage, authority = validate_common_records(
        envelope, terminal_path, success=success,
        trusted_source_root=trusted_source_root,
        authority_envelope=authority_envelope,
        require_trusted_source_match=require_trusted_source_match)
    if success:
        return validate_success(
            envelope, core, terminal, manifest, lineage, authority,
            trusted_source_root, require_preflight_identity)
    return validate_failure(
        core, terminal, manifest, lineage, authority, trusted_source_root,
        require_preflight_identity)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--envelope", required=True, type=Path)
    parser.add_argument("--trusted-source-root", required=True, type=Path)
    parser.add_argument("--authority-envelope", type=Path)
    parser.add_argument("--require-trusted-source-match", action="store_true")
    parser.add_argument("--require-preflight-identity", action="store_true")
    options = parser.parse_args(argv)
    try:
        result = verify_envelope(
            options.envelope, options.trusted_source_root,
            require_trusted_source_match=options.require_trusted_source_match,
            require_preflight_identity=options.require_preflight_identity,
            authority_envelope=options.authority_envelope)
        sys.stdout.buffer.write(canonical_json_bytes(result))
        return 0
    except (ReplayError, OSError, UnicodeError, json.JSONDecodeError,
            subprocess.SubprocessError) as error:
        print(f"conditioned v19 wrapper replay failed: {error}",
              file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = (
    "ATTEMPT_LINEAGE_SCHEMA", "CLAIM_CEILING", "COMPLETE_CORE_KEYS",
    "COMPLETE_CORE_SCHEMA", "CONTROLLER_BINDINGS", "EVIDENCE_CLASS",
    "FAILED_CORE_KEYS", "FAILED_CORE_SCHEMA", "FAILURE_TERMINAL_KEYS",
    "FAILURE_TERMINAL_SCHEMA", "GENERATION", "PREREGISTRATION_SHA256",
    "ReplayError", "SUCCESS_TERMINAL_KEYS", "SUCCESS_TERMINAL_SCHEMA",
    "canonical_json_bytes", "verify_envelope",
)
