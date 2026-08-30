#!/usr/bin/env python3
"""Qualify the exact GF8 K=65/R=65/64-byte packed AVX2 encoder.

The candidate and mature control are two runtime modes of one frozen physical
executable.  Twenty-five balanced ABBA rounds compare ordinary batch-of-one
and one-shot encoding against both that control and canonical Leopard main.
Immediate selector boundaries, the retained 8-KiB path, and hot-code layout
canaries must prove that the selector is inert.  Exact-main comparisons for
the retained path and the two established balanced terminals are additionally
held to a two-percent regression floor; the remaining layout comparisons are
retained as context rather than used as unestablished release gates.

The sealed exact-main authority tree remains immutable and non-executable.
``--main`` names a byte-identical executable copy installed outside that tree;
the pinned SHA-256, sealed authority replay, and child payload jointly bind it.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import fcntl
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


SUPPORT_PATH = Path(__file__).resolve().with_name(
    "run_k66r16_b64_tail_abba.py")


def load_support() -> Any:
    specification = importlib.util.spec_from_file_location(
        "k65r65_b64_packed_terminal_evidence_support", SUPPORT_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load qualification support: {SUPPORT_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()
BASE = SUPPORT.BASE
PRODUCTION_SUPPORT = SUPPORT.SUPPORT
BASE.__doc__ = __doc__

RAW_SCHEMA_V1 = "leopard2-gf8-k65r65-b64-packed-terminal-abba/v1"
SUMMARY_SCHEMA_V1 = \
    "leopard2-gf8-k65r65-b64-packed-terminal-summary/v1"
RAW_SCHEMA_V2 = "leopard2-gf8-k65r65-b64-packed-terminal-abba/v2"
PRELIMINARY_SUMMARY_SCHEMA_V2 = \
    "leopard2-gf8-k65r65-b64-packed-terminal-preliminary-summary/v2"
FINAL_SUMMARY_SCHEMA_V2 = \
    "leopard2-gf8-k65r65-b64-packed-terminal-summary/v2"
AUTHORITY_BINDING_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-exact-main-authority-binding/v2"
ATTEMPT_PREREGISTRATION_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-attempt-preregistration/v2"
ATTEMPT_TERMINAL_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-attempt-terminal/v2"
GENERATION = 2

BASE.SCHEMA = RAW_SCHEMA_V2
BASE.SUMMARY_SCHEMA = PRELIMINARY_SUMMARY_SCHEMA_V2
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L33g_k65r65_b64_packed_terminal_modeE"
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.CANDIDATE_SCHEMA = "leopard2-benchmark-v31"
BASE.CONTROL_SCHEMA = "leopard2-benchmark-v31"
BASE.CONTROL_EXTRA_ARGUMENTS = ()
BASE.CONTROL_BUILD_MARKER = \
    "k65r65_b64_packed_terminal_diagnostic_disabled"
BASE.REQUIRE_EXPECTED_IDENTITIES = True
BASE.REQUIRE_BUILD_CLOSURE = True
BASE.REQUIRE_FULL_ELF_IDENTITY = True
BASE.TARGET_CONTROL_FLOOR = 1.05
BASE.TARGET_MAIN_FLOOR = 1.05
BASE.NEIGHBOR_EQUIVALENCE_LOWER = 1.0 / 1.02
BASE.NEIGHBOR_EQUIVALENCE_UPPER = 1.02
BASE.RETAINED_MAIN_FLOOR = 0.98
BASE.MAX_ISOLATION_ATTEMPTS = 3
BASE.RUNNER_PATH = Path(__file__).resolve()
_INHERITED_RUNNER_DEPENDENCIES = tuple(BASE.RUNNER_DEPENDENCIES)

SOURCE_ROOT = Path(__file__).resolve().parents[3]
TOOLS_ROOT = SOURCE_ROOT / "tools"
IDENTITY_CONTRACT_PATH = TOOLS_ROOT / "leopard2_exact_main_baseline.py"
AUTHORITY_RECORD_PATH = \
    TOOLS_ROOT / "leopard2_exact_main_baseline_record.py"
AUTHORITY_VERIFIER_PATH = \
    TOOLS_ROOT / "leopard2_exact_main_baseline_verifier.py"


def _load_authority_module(
    module_name: str, path: Path,
) -> Any:
    resolved = path.resolve(strict=True)
    existing = sys.modules.get(module_name)
    if existing is not None and Path(
            str(getattr(existing, "__file__", ""))).resolve() == resolved:
        return existing
    specification = importlib.util.spec_from_file_location(
        module_name, resolved)
    BASE.require(
        specification is not None and specification.loader is not None,
        f"cannot load exact-main authority dependency: {resolved}")
    module = importlib.util.module_from_spec(specification)
    previous = sys.modules.get(module_name)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        if previous is None:
            sys.modules.pop(module_name, None)
        else:
            sys.modules[module_name] = previous
        raise
    return module


IDENTITY_CONTRACT = _load_authority_module(
    "leopard2_exact_main_baseline", IDENTITY_CONTRACT_PATH)
AUTHORITY_RECORD = _load_authority_module(
    "leopard2_exact_main_baseline_record", AUTHORITY_RECORD_PATH)
AUTHORITY_VERIFIER = _load_authority_module(
    "leopard2_exact_main_baseline_verifier", AUTHORITY_VERIFIER_PATH)

BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    *_INHERITED_RUNNER_DEPENDENCIES,
    IDENTITY_CONTRACT_PATH.resolve(),
    AUTHORITY_RECORD_PATH.resolve(),
    AUTHORITY_VERIFIER_PATH.resolve(),
)

EXACT_MAIN_AUTHORITY_SCHEMA = \
    "leopard2-gf8-exact-main-pure-avx2-baseline/v1"
EXACT_MAIN_AUTHORITY_RECORD_SHA256 = \
    "ce3d6e647fd098f558ee01855523c0643b8c175578315d02a168a87992e5ae01"
EXACT_MAIN_AUTHORITY_LEDGER_SHA256 = \
    "20ca2823aed2330474cba186b0029e1cfd39dfced3b6c76fa966dcac1ea3fae2"
EXACT_MAIN_EXECUTABLE_SHA256 = \
    "e1d849056a7f061b127e14c0d8f71165ceaa9065f2d28d7469cae33e7e19eadf"
EXACT_MAIN_ARCHIVE_SHA256 = \
    "9a3e3501f96da63c19d6c0b46fd4040967be4ec584c20ad5102cf1072e9b32a3"
EXACT_MAIN_COMBINED_SHA256 = \
    "ddfef166af6c1dafd989019f87694526693623fa4ea2aa9e4d74f97c012fa093"
EXACT_MAIN_SOURCE_COMMIT = \
    "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
EXACT_MAIN_SOURCE_TREE = "b7c8830d96a978f6ec14fe747095f066e351ae72"
EXACT_MAIN_ADAPTER_COMMIT = \
    "6d4f690fe5ba7cf08f2f8f80a263765266b462e0"
HISTORICAL_EXECUTABLE_SHA256 = \
    "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93"
HISTORICAL_ARCHIVE_SHA256 = \
    "49250fb7528df955e57e96f84b309cf0c26f203122d112f3f501336c0eb131d0"
HISTORICAL_NON_AUTHORITY = (
    {
        "authority": False,
        "label": "historical_archive_sha256",
        "value": HISTORICAL_ARCHIVE_SHA256,
    },
    {
        "authority": False,
        "label": "historical_executable_sha256",
        "value": HISTORICAL_EXECUTABLE_SHA256,
    },
)
ATTEMPT_BUDGET = 3

BASE.CANONICAL_MAIN_SHA256 = EXACT_MAIN_EXECUTABLE_SHA256

BASE.require(
    len({
        RAW_SCHEMA_V1, SUMMARY_SCHEMA_V1, RAW_SCHEMA_V2,
        PRELIMINARY_SUMMARY_SCHEMA_V2, FINAL_SUMMARY_SCHEMA_V2,
        AUTHORITY_BINDING_SCHEMA, ATTEMPT_PREREGISTRATION_SCHEMA,
        ATTEMPT_TERMINAL_SCHEMA,
    }) == 8 and
    AUTHORITY_RECORD.AUTHORITY_SCHEMA == EXACT_MAIN_AUTHORITY_SCHEMA and
    AUTHORITY_RECORD.BENCHMARK_SCHEMA == "leopard-main-benchmark-v1" and
    AUTHORITY_RECORD.SEAL_PROTOCOL == "owner-only-tree-sha256sums/v1" and
    all(item["authority"] is False for item in HISTORICAL_NON_AUTHORITY) and
    {item["value"] for item in HISTORICAL_NON_AUTHORITY}.isdisjoint({
        EXACT_MAIN_AUTHORITY_RECORD_SHA256,
        EXACT_MAIN_AUTHORITY_LEDGER_SHA256,
        EXACT_MAIN_EXECUTABLE_SHA256,
        EXACT_MAIN_ARCHIVE_SHA256,
        EXACT_MAIN_COMBINED_SHA256,
    }),
    "K65 v2 schema or exact-main authority constants are ambiguous")

AUTHORITY_BINDING_KEYS = frozenset((
    "schema", "lane_root", "authority_schema", "record_sha256",
    "ledger_sha256", "verifier_schema", "verifier_verdict_sha256",
    "executable_sha256", "archive_sha256", "combined_sha256",
    "source_commit", "source_tree", "adapter_commit", "pure_avx2",
    "historical_non_authority",
))
ATTEMPT_TERMINAL_KEYS = frozenset((
    "schema", "generation", "attempt", "attempt_budget", "outcome",
    "promotable", "output_root", "summary_schema", "summary_sha256",
    "raw_sha256", "failure_sha256", "authority_record_sha256",
))
ATTEMPT_PREREGISTRATION_KEYS = frozenset((
    "schema", "generation", "attempt", "attempt_budget",
    "attempts_remaining", "lineage", "decision_evidence",
    "statistical_campaign_retries", "outcome_selected_cells",
))
FINAL_SUMMARY_KEYS = frozenset((
    "schema", "status", "source_commit", "source_tree", "main_commit",
    "cell_count", "process_count", "discarded_process_count",
    "discarded_round_attempts", "all_digests_matched",
    "all_rounds_zero_sibling_nonidle", "cells", "binary_sha256",
    "candidate_control_executable_sections_sha256", "mode_words",
    "raw_sha256", "generation", "exact_main_authority",
    "attempt_preregistration", "preliminary_generic_status",
    "finalization_cpu", "preliminary_summary_sha256",
    "acceptance_contract", "gate_failures", "promotable",
))


def _canonical_json_bytes(value: Any) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"),
            ensure_ascii=True, allow_nan=False).encode("ascii")
    except (TypeError, ValueError, UnicodeError) as error:
        raise BASE.EvidenceError("value is not canonical JSON") from error


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _strict_json_path(path: Path, label: str) -> Any:
    BASE.require(path.is_file() and not path.is_symlink(),
                 f"{label} is not one regular non-symlink file")
    try:
        return IDENTITY_CONTRACT.strict_json_loads(path.read_bytes(), label)
    except Exception as error:
        raise BASE.EvidenceError(f"{label} is not strict JSON") from error


def _lower_sha256_or_none(value: Any, label: str) -> str | None:
    BASE.require(
        value is None or
        (type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value)),
        f"{label} is not null or one lowercase SHA-256")
    return value


def _canonical_absolute_posix_path(value: Any, label: str) -> str:
    BASE.require(
        type(value) is str and 1 < len(value) <= 4096 and
        value.startswith("/") and "\x00" not in value and
        Path(value).as_posix() == value and
        all(part not in {"", ".", ".."} for part in value.split("/")[1:]),
        f"{label} is not absolute canonical POSIX")
    return value


def validate_exact_main_authority_binding(value: Any) -> dict[str, Any]:
    """Validate and detach the closed v2 authority citation."""
    BASE.require(type(value) is dict and set(value) == AUTHORITY_BINDING_KEYS,
                 "exact-main authority binding keys changed")
    expected = {
        "schema": AUTHORITY_BINDING_SCHEMA,
        "authority_schema": EXACT_MAIN_AUTHORITY_SCHEMA,
        "record_sha256": EXACT_MAIN_AUTHORITY_RECORD_SHA256,
        "ledger_sha256": EXACT_MAIN_AUTHORITY_LEDGER_SHA256,
        "verifier_schema": AUTHORITY_VERIFIER.VERIFIER_SCHEMA,
        "executable_sha256": EXACT_MAIN_EXECUTABLE_SHA256,
        "archive_sha256": EXACT_MAIN_ARCHIVE_SHA256,
        "combined_sha256": EXACT_MAIN_COMBINED_SHA256,
        "source_commit": EXACT_MAIN_SOURCE_COMMIT,
        "source_tree": EXACT_MAIN_SOURCE_TREE,
        "adapter_commit": EXACT_MAIN_ADAPTER_COMMIT,
        "pure_avx2": True,
        "historical_non_authority": list(HISTORICAL_NON_AUTHORITY),
    }
    for name, expected_value in expected.items():
        BASE.require(
            type(value.get(name)) is type(expected_value) and
            value[name] == expected_value,
            f"exact-main authority binding {name} changed")
    _canonical_absolute_posix_path(
        value.get("lane_root"), "exact-main authority binding lane root")
    verdict_sha = value.get("verifier_verdict_sha256")
    BASE.require(
        type(verdict_sha) is str and
        re.fullmatch(r"[0-9a-f]{64}", verdict_sha) is not None,
        "exact-main verifier verdict SHA-256 is malformed")
    return copy.deepcopy(value)


def bind_exact_main_authority(lane_root: Path) -> dict[str, Any]:
    """Independently replay and bind the promoted exact-main authority lane."""
    BASE.require(not lane_root.is_symlink(),
                 "exact-main authority lane is a symbolic link")
    resolved = lane_root.resolve(strict=True)
    BASE.require(resolved.is_dir(),
                 "exact-main authority lane is not one directory")
    verdict = AUTHORITY_VERIFIER.verify_sealed_lane(resolved)
    BASE.require(
        type(verdict) is dict and
        verdict.get("schema") == AUTHORITY_VERIFIER.VERIFIER_SCHEMA and
        verdict.get("outcome") == "promoted_authority" and
        verdict.get("promoted") is True and
        type(verdict.get("record")) is dict and
        verdict["record"].get("schema") == EXACT_MAIN_AUTHORITY_SCHEMA and
        verdict["record"].get("record_sha256") ==
            EXACT_MAIN_AUTHORITY_RECORD_SHA256 and
        verdict["record"].get("canonical_bytes_identical") is True and
        type(verdict.get("seal")) is dict and
        verdict["seal"].get("protocol") == AUTHORITY_RECORD.SEAL_PROTOCOL and
        verdict["seal"].get("sha256sums_sha256") ==
            EXACT_MAIN_AUTHORITY_LEDGER_SHA256,
        "sealed exact-main lane did not replay as the pinned authority")
    record_path = resolved / "baseline-authority.json"
    record = AUTHORITY_RECORD.load_baseline_authority_record(
        record_path.read_bytes())
    builds = record.get("builds")
    identity = record.get("identity")
    source = record.get("source")
    BASE.require(
        record.get("schema") == EXACT_MAIN_AUTHORITY_SCHEMA and
        record.get("record_sha256") == EXACT_MAIN_AUTHORITY_RECORD_SHA256 and
        type(record.get("lane")) is dict and
        record["lane"].get("root") == str(resolved) and
        record["lane"].get("seal_protocol") ==
            AUTHORITY_RECORD.SEAL_PROTOCOL and
        record["lane"].get("attempt_budget") == ATTEMPT_BUDGET and
        type(builds) is dict and type(identity) is dict and
        type(source) is dict,
        "exact-main authority record identity changed")
    for role in ("canonical_first", "canonical_second"):
        BASE.require(
            type(builds.get(role)) is dict and
            builds[role].get("executable", {}).get("sha256") ==
                EXACT_MAIN_EXECUTABLE_SHA256 and
            builds[role].get("archive", {}).get("sha256") ==
                EXACT_MAIN_ARCHIVE_SHA256 and
            identity.get(role, {}).get("combined_sha256") ==
                EXACT_MAIN_COMBINED_SHA256,
            f"exact-main {role} authority identity changed")
    BASE.require(
        type(builds.get("path_variant")) is dict and
        identity.get("path_variant", {}).get("combined_sha256") ==
            EXACT_MAIN_COMBINED_SHA256 and
        identity.get("combined_sha256") == EXACT_MAIN_COMBINED_SHA256 and
        identity.get("normalized_match") is True and
        identity.get("canonical_raw_bytes_identical") is True and
        identity.get("path_variant_raw_bytes_differ") is True and
        source.get("baseline", {}).get("commit") == EXACT_MAIN_SOURCE_COMMIT and
        source.get("baseline", {}).get("tree") == EXACT_MAIN_SOURCE_TREE and
        source.get("adapter_repository", {}).get("commit") ==
            EXACT_MAIN_ADAPTER_COMMIT and
        record.get("configure", {}).get("pure_avx2") is True and
        record.get("superseded_references") ==
            list(HISTORICAL_NON_AUTHORITY) and
        type(record.get("promotion")) is dict and
        all(field is True for name, field in record["promotion"].items()
            if name != "gate") and
        all(item.get("authority") is False
            for item in record["superseded_references"]) and
        _sha256_file(resolved / "SHA256SUMS") ==
            EXACT_MAIN_AUTHORITY_LEDGER_SHA256,
        "exact-main authority provenance or promotion claims changed")
    binding = {
        "schema": AUTHORITY_BINDING_SCHEMA,
        "lane_root": str(resolved),
        "authority_schema": EXACT_MAIN_AUTHORITY_SCHEMA,
        "record_sha256": EXACT_MAIN_AUTHORITY_RECORD_SHA256,
        "ledger_sha256": EXACT_MAIN_AUTHORITY_LEDGER_SHA256,
        "verifier_schema": verdict["schema"],
        "verifier_verdict_sha256": verdict["verdict_sha256"],
        "executable_sha256": EXACT_MAIN_EXECUTABLE_SHA256,
        "archive_sha256": EXACT_MAIN_ARCHIVE_SHA256,
        "combined_sha256": EXACT_MAIN_COMBINED_SHA256,
        "source_commit": EXACT_MAIN_SOURCE_COMMIT,
        "source_tree": EXACT_MAIN_SOURCE_TREE,
        "adapter_commit": EXACT_MAIN_ADAPTER_COMMIT,
        "pure_avx2": True,
        "historical_non_authority": list(HISTORICAL_NON_AUTHORITY),
    }
    return validate_exact_main_authority_binding(binding)


def validate_exact_main_launch_copy(
    path: Path, authority: Mapping[str, Any],
) -> Path:
    """Require an executable copy without mutating the sealed authority."""
    binding = validate_exact_main_authority_binding(authority)
    try:
        resolved = path.resolve(strict=True)
        lane_root = Path(binding["lane_root"]).resolve(strict=True)
    except OSError as error:
        raise BASE.EvidenceError(
            "exact-main launch copy or authority lane is absent") from error
    BASE.require(
        not resolved.is_relative_to(lane_root),
        "--main must be an executable copy outside the sealed authority lane")
    BASE.require(
        resolved.is_file() and os.access(resolved, os.X_OK),
        "exact-main launch copy is not executable")
    BASE.require(
        _sha256_file(resolved) == EXACT_MAIN_EXECUTABLE_SHA256,
        "exact-main launch copy differs from the sealed authority executable")
    return resolved


def validate_attempt_terminal(
    value: Any, *, verify_files: bool, allow_accepted: bool = False,
) -> dict[str, Any]:
    """Validate one prior, non-promotable v2 campaign terminal record."""
    BASE.require(type(value) is dict and set(value) == ATTEMPT_TERMINAL_KEYS,
                 "attempt terminal record keys changed")
    BASE.require(
        value.get("schema") == ATTEMPT_TERMINAL_SCHEMA and
        type(value.get("generation")) is int and
        value["generation"] == GENERATION and
        type(value.get("attempt")) is int and
        1 <= value["attempt"] <= ATTEMPT_BUDGET and
        type(value.get("attempt_budget")) is int and
        value["attempt_budget"] == ATTEMPT_BUDGET and
        value.get("outcome") in (
            {"accepted", "rejected", "failed"} if allow_accepted else
            {"rejected", "failed"}) and
        type(value.get("promotable")) is bool and
        value["promotable"] is (value.get("outcome") == "accepted") and
        value.get("authority_record_sha256") ==
            EXACT_MAIN_AUTHORITY_RECORD_SHA256,
        "attempt terminal metadata is invalid")
    output_root = _canonical_absolute_posix_path(
        value.get("output_root"), "attempt terminal output root")
    for name in ("summary_sha256", "raw_sha256", "failure_sha256"):
        _lower_sha256_or_none(value.get(name), f"attempt terminal {name}")
    if value["outcome"] in {"accepted", "rejected"}:
        BASE.require(
            value.get("summary_schema") == FINAL_SUMMARY_SCHEMA_V2 and
            value.get("summary_sha256") is not None and
            value.get("raw_sha256") is not None and
            value.get("failure_sha256") is None,
            "completed attempt terminal does not bind final evidence")
    else:
        BASE.require(
            value.get("summary_schema") in {
                "absent", FINAL_SUMMARY_SCHEMA_V2} and
            ((value["summary_schema"] == "absent") ==
             (value.get("summary_sha256") is None)) and
            value.get("failure_sha256") is not None,
            "failed attempt terminal artifact bindings are inconsistent")
    if verify_files:
        root = Path(output_root)
        BASE.require(
            root.resolve(strict=True) == root and root.is_dir() and
            not root.is_symlink(),
            "prior attempt output root is not retained canonically")
        terminal_path = root / "attempt-terminal.json"
        observed_terminal = _strict_json_path(
            terminal_path, "prior attempt terminal")
        BASE.require(observed_terminal == value,
                     "prior attempt terminal bytes differ from lineage")
        bindings = (
            ("summary.json", value["summary_sha256"]),
            ("raw.json", value["raw_sha256"]),
            ("failure.json", value["failure_sha256"]),
        )
        for name, expected_sha in bindings:
            path = root / name
            if expected_sha is None:
                BASE.require(not path.exists() and not path.is_symlink(),
                             f"prior attempt unexpectedly contains {name}")
            else:
                BASE.require(
                    path.is_file() and not path.is_symlink() and
                    _sha256_file(path) == expected_sha,
                    f"prior attempt {name} differs from lineage")
        if value["summary_schema"] == FINAL_SUMMARY_SCHEMA_V2:
            summary = _strict_json_path(root / "summary.json",
                                        "prior final summary")
            BASE.require(
                type(summary) is dict and
                summary.get("schema") == FINAL_SUMMARY_SCHEMA_V2 and
                summary.get("status") == (
                    "accepted" if value["outcome"] == "accepted" else
                    "rejected") and
                summary.get("promotable") is value["promotable"],
                "attempt final summary differs from its terminal")
    return copy.deepcopy(value)


def validate_attempt_preregistration(
    attempt: int, lineage: Any, *, verify_files: bool = True,
) -> dict[str, Any]:
    """Freeze one gapless campaign attempt without outcome selection."""
    BASE.require(type(attempt) is int and 1 <= attempt <= ATTEMPT_BUDGET,
                 "campaign attempt is outside the preregistered budget")
    BASE.require(type(lineage) is list and len(lineage) == attempt - 1,
                 "attempt lineage is not the exact prior-attempt prefix")
    validated = [
        validate_attempt_terminal(item, verify_files=verify_files)
        for item in lineage
    ]
    BASE.require(
        [item["attempt"] for item in validated] == list(range(1, attempt)) and
        len({item["output_root"] for item in validated}) == len(validated) and
        all(item["outcome"] in {"rejected", "failed"} and
            item["promotable"] is False for item in validated),
        "attempt lineage is reordered, duplicated, selected, or incomplete")
    value = {
        "schema": ATTEMPT_PREREGISTRATION_SCHEMA,
        "generation": GENERATION,
        "attempt": attempt,
        "attempt_budget": ATTEMPT_BUDGET,
        "attempts_remaining": ATTEMPT_BUDGET - attempt,
        "lineage": validated,
        "decision_evidence": "this standalone 25-round campaign only",
        "statistical_campaign_retries": 0,
        "outcome_selected_cells": False,
    }
    BASE.require(set(value) == ATTEMPT_PREREGISTRATION_KEYS,
                 "attempt preregistration keys changed")
    return copy.deepcopy(value)


def _load_attempt_lineage(path: Path | None, attempt: int) -> list[Any]:
    if attempt == 1:
        BASE.require(path is None,
                     "attempt 1 must not cite an attempt-lineage file")
        return []
    BASE.require(path is not None,
                 "later attempts require the retained attempt lineage")
    value = _strict_json_path(path, "attempt lineage")
    BASE.require(type(value) is list, "attempt lineage JSON is not an array")
    return value

CONFIRMATORY_ROUNDS = 25
FIXED_ITERATIONS = 31
FIXED_WARMUP = 64
FIXED_REUSE = 8192
SELECTOR_ARGUMENT = "--k65r65-b64-packed-terminal-mode"
SELECTOR_CONTRACT = (
    "LEGACY_HIGH_V1,GF8,AVX2,K=65,R=65,T=128,B=64,"
    "native_layout,auto_encode,one_shot_and_one_item_batch"
)
ORDINARY_API_MARKER = "leo2_encode_batch:item_count=1:no_preflight_scratch"
ONE_SHOT_API_MARKER = "leo2_encode"
CHILD_TIMEOUT_SETUP_SECONDS = 120
# The retained B8192 child takes about 36 seconds on the qualification host.
# This deliberately conservative 1-GiB/s logical-traffic floor leaves more
# than a twenty-fold wall-time margin without turning a hung child into a
# multi-hour wait.
CHILD_TIMEOUT_LOGICAL_BYTES_PER_SECOND = 1024 * 1024 * 1024


def parse_arguments() -> argparse.Namespace:
    """Expose no timing-policy escape hatch for confirmatory evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--candidate-sha256", required=True)
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--control-sha256", required=True)
    parser.add_argument(
        "--main", required=True, type=Path,
        help=(
            "executable copy outside the sealed authority lane; install the "
            "canonical-first authority artifact byte-for-byte with mode 0555"
        ))
    parser.add_argument(
        "--main-sha256", default=BASE.CANONICAL_MAIN_SHA256)
    parser.add_argument(
        "--exact-main-authority-lane", required=True, type=Path)
    parser.add_argument("--generation", type=int, choices=(GENERATION,),
                        default=GENERATION)
    parser.add_argument("--attempt", type=int,
                        choices=range(1, ATTEMPT_BUDGET + 1), required=True)
    parser.add_argument("--attempt-lineage", type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.set_defaults(
        iterations=FIXED_ITERATIONS,
        warmup=FIXED_WARMUP,
        rounds=CONFIRMATORY_ROUNDS,
        candidate_archive=None,
        candidate_archive_sha256=None,
        control_archive=None,
        control_archive_sha256=None,
        candidate_compile_commands=None,
        candidate_compile_commands_sha256=None,
        control_compile_commands=None,
        control_compile_commands_sha256=None,
    )
    options = parser.parse_args()
    BASE.require(
        options.main_sha256 == BASE.CANONICAL_MAIN_SHA256,
        "exact-main executable SHA-256 differs from the canonical AVX2 "
        "Leopard-main benchmark")
    BASE.require(
        re.fullmatch(r"[0-9a-f]{64}", options.candidate_sha256) is not None
        and options.control_sha256 == options.candidate_sha256,
        "candidate/control must provide one identical lowercase SHA-256")
    options.exact_main_authority = bind_exact_main_authority(
        options.exact_main_authority_lane)
    options.main = validate_exact_main_launch_copy(
        options.main, options.exact_main_authority)
    lineage = _load_attempt_lineage(options.attempt_lineage, options.attempt)
    options.attempt_preregistration = validate_attempt_preregistration(
        options.attempt, lineage, verify_files=True)
    return options


BASE.parse_arguments = parse_arguments


_BASE_SELECT_ROUND_ORDERS = BASE.select_round_orders


def select_round_orders(
    orders: Sequence[Sequence[str]], requested_rounds: int | None,
) -> tuple[tuple[str, ...], ...]:
    BASE.require(requested_rounds == CONFIRMATORY_ROUNDS,
                 "K65/R65/B64 evidence requires exactly 25 rounds")
    return _BASE_SELECT_ROUND_ORDERS(orders, requested_rounds)


BASE.select_round_orders = select_round_orders


_BASE_CONFIDENCE_INTERVAL = BASE.confidence_interval


def confidence_interval(values: Sequence[float]) -> dict[str, Any]:
    BASE.require(len(values) == CONFIRMATORY_ROUNDS,
                 "K65/R65/B64 evidence requires 25 independent contrasts")
    result = _BASE_CONFIDENCE_INTERVAL(values)
    BASE.require(
        result.get("confidence_level") == 0.95 and
        result.get("degrees_of_freedom") == 24 and
        result.get("t_critical") == PRODUCTION_SUPPORT.T95_DF24,
        "K65/R65/B64 confidence interval is not the fixed df=24 interval")
    return result


BASE.confidence_interval = confidence_interval


def cells() -> list[dict[str, Any]]:
    # main_floor is deliberately present in the attested cell payload.  Only
    # cells with a prior exact-main baseline become release gates.
    values = (
        ("target-k65-r65-b64-q1", 65, 65, 64,
         "target", True, 1.05),
        ("k-control-k64-r65-b64-q1", 64, 65, 64,
         "neighbor", False, None),
        ("k-control-k66-r65-b64-q1", 66, 65, 64,
         "neighbor", False, None),
        ("r-control-k65-r64-b64-q1", 65, 64, 64,
         "neighbor", False, None),
        ("r-control-k65-r66-b64-q1", 65, 66, 64,
         "neighbor", False, None),
        ("byte-control-k65-r65-b63-q1", 65, 65, 63,
         "neighbor", False, None),
        ("byte-control-k65-r65-b65-q1", 65, 65, 65,
         "neighbor", False, None),
        ("retained-k65-r65-b8192-q1", 65, 65, 8192,
         "neighbor", True, 0.98),
        ("balanced-layout-k64-r64-b64-q1", 64, 64, 64,
         "neighbor", True, 0.98),
        ("balanced-layout-k128-r128-b64-q1", 128, 128, 64,
         "neighbor", True, 0.98),
        ("layout-context-k79-r32-b64-q1", 79, 32, 64,
         "neighbor", True, None),
        ("layout-context-k62-r8-b64-q1", 62, 8, 64,
         "neighbor", True, None),
        ("layout-context-k66-r16-b64-q1", 66, 16, 64,
         "neighbor", True, None),
    )
    result = []
    for index, (name, k, r, shard_bytes, role, compare_main,
                main_floor) in enumerate(values):
        result.append({
            "id": name,
            "K": k,
            "R": r,
            "bytes": shard_bytes,
            "loss": 1,
            "batch": 1,
            "reuse": FIXED_REUSE,
            "role": role,
            "compare_main": compare_main,
            "main_floor": main_floor,
            "measure_one_shot": True,
            "seed": 0x6565B640 + index,
        })
    BASE.require(
        len(result) == 13 and
        sum(cell["role"] == "target" for cell in result) == 1 and
        len({cell["id"] for cell in result}) == len(result) and
        len({cell["seed"] for cell in result}) == len(result) and
        {cell["id"] for cell in result if cell["compare_main"]} == {
            "target-k65-r65-b64-q1",
            "retained-k65-r65-b8192-q1",
            "balanced-layout-k64-r64-b64-q1",
            "balanced-layout-k128-r128-b64-q1",
            "layout-context-k79-r32-b64-q1",
            "layout-context-k62-r8-b64-q1",
            "layout-context-k66-r16-b64-q1",
        } and
        {cell["id"] for cell in result
         if cell["role"] == "neighbor" and cell["main_floor"] == 0.98} == {
            "retained-k65-r65-b8192-q1",
            "balanced-layout-k64-r64-b64-q1",
            "balanced-layout-k128-r128-b64-q1",
        } and
        all(cell["loss"] == 1 and cell["batch"] == 1 and
            cell["reuse"] == FIXED_REUSE and
            cell["measure_one_shot"] is True
            for cell in result),
        "K65/R65/B64 qualification matrix is incomplete")
    return result


BASE.cells = cells

MATRIX_SHA256 = \
    "ab9572c4101b2af5eda4b7cfab17e979239f698c4c1196e660cf6f5e3f4af27c"
V1_SEMANTIC_PROJECTION_SHA256 = \
    "1f8bbf9e0f9daab2b2e4616e5895dca3639b1b3b1212383c9f729b5292b0854f"


def generation_projection(generation: int) -> dict[str, Any]:
    """Return the frozen replay surface for K65 evidence generations."""
    BASE.require(type(generation) is int and generation in (1, GENERATION),
                 "unknown K65 evidence generation")
    projection = {
        "raw_schema": RAW_SCHEMA_V1 if generation == 1 else RAW_SCHEMA_V2,
        "summary_schema": SUMMARY_SCHEMA_V1 if generation == 1 else
            FINAL_SUMMARY_SCHEMA_V2,
        "mode_symbol": BASE.MODE_SYMBOL,
        "candidate_schema": BASE.CANDIDATE_SCHEMA,
        "control_schema": BASE.CONTROL_SCHEMA,
        "control_build_marker": BASE.CONTROL_BUILD_MARKER,
        "target_control_floor": BASE.TARGET_CONTROL_FLOOR,
        "target_main_floor": BASE.TARGET_MAIN_FLOOR,
        "neighbor_band": [
            BASE.NEIGHBOR_EQUIVALENCE_LOWER,
            BASE.NEIGHBOR_EQUIVALENCE_UPPER,
        ],
        "retained_main_floor": BASE.RETAINED_MAIN_FLOOR,
        "canonical_main_sha256": HISTORICAL_EXECUTABLE_SHA256
            if generation == 1 else EXACT_MAIN_EXECUTABLE_SHA256,
        "max_isolation_attempts": BASE.MAX_ISOLATION_ATTEMPTS,
        "rounds": CONFIRMATORY_ROUNDS,
        "iterations": FIXED_ITERATIONS,
        "warmup": FIXED_WARMUP,
        "reuse": FIXED_REUSE,
        "selector_argument": SELECTOR_ARGUMENT,
        "selector_contract": SELECTOR_CONTRACT,
        "ordinary_api": ORDINARY_API_MARKER,
        "one_shot_api": ONE_SHOT_API_MARKER,
        "cells": cells(),
    }
    if generation == GENERATION:
        projection["exact_main_authority"] = {
            "record_sha256": EXACT_MAIN_AUTHORITY_RECORD_SHA256,
            "ledger_sha256": EXACT_MAIN_AUTHORITY_LEDGER_SHA256,
            "executable_sha256": EXACT_MAIN_EXECUTABLE_SHA256,
            "archive_sha256": EXACT_MAIN_ARCHIVE_SHA256,
            "combined_sha256": EXACT_MAIN_COMBINED_SHA256,
            "source_commit": EXACT_MAIN_SOURCE_COMMIT,
            "source_tree": EXACT_MAIN_SOURCE_TREE,
            "adapter_commit": EXACT_MAIN_ADAPTER_COMMIT,
        }
    return copy.deepcopy(projection)


BASE.require(
    _canonical_sha256(cells()) == MATRIX_SHA256 and
    _canonical_sha256(generation_projection(1)) ==
        V1_SEMANTIC_PROJECTION_SHA256,
    "K65 v1 semantic projection or 13-cell matrix changed")


_BASE_BENCHMARK_COMMAND = BASE.benchmark_command


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    command = _BASE_BENCHMARK_COMMAND(
        implementation, executable, cell, cpu, iterations, warmup)
    if implementation == "main":
        return command
    BASE.require(SELECTOR_ARGUMENT not in command,
                 "K65/R65 selector was already present in the command")
    BASE.require(command[-2:] == ["--json", "-"],
                 "benchmark JSON command suffix changed")
    mode = "0" if implementation == "control" else "1"
    BASE.require(implementation in {"candidate", "control"},
                 f"unexpected Leopard2 implementation label: {implementation}")
    return command[:-2] + [SELECTOR_ARGUMENT, mode, "--json", "-"]


BASE.benchmark_command = benchmark_command


def child_timeout_budget(
    implementation: str,
    cell: Mapping[str, Any],
    iterations: int,
    warmup: int,
) -> dict[str, int]:
    """Budget long fixed-reuse children from their timed logical traffic."""
    BASE.require(implementation in {"candidate", "control", "main"},
                 f"unexpected timeout implementation: {implementation}")
    values = (
        cell.get("K"), cell.get("R"), cell.get("bytes"),
        cell.get("batch"), cell.get("reuse"), iterations, warmup)
    BASE.require(all(type(value) is int and value > 0 for value in values),
                 "timeout workload dimensions are invalid")
    measured_metric_count = 2
    if implementation != "main" and cell.get("measure_one_shot") is True:
        measured_metric_count += 1
    calls_per_metric = iterations * int(cell["reuse"]) + warmup
    logical_bytes_per_call = (
        (int(cell["K"]) + int(cell["R"])) * int(cell["bytes"]) *
        int(cell["batch"]))
    logical_byte_visits = (
        measured_metric_count * calls_per_metric * logical_bytes_per_call)
    timeout_seconds = CHILD_TIMEOUT_SETUP_SECONDS + math.ceil(
        logical_byte_visits / CHILD_TIMEOUT_LOGICAL_BYTES_PER_SECOND)
    return {
        "timeout_seconds": timeout_seconds,
        "setup_seconds": CHILD_TIMEOUT_SETUP_SECONDS,
        "logical_bytes_per_second_floor":
            CHILD_TIMEOUT_LOGICAL_BYTES_PER_SECOND,
        "measured_metric_count": measured_metric_count,
        "calls_per_metric": calls_per_metric,
        "logical_bytes_per_call": logical_bytes_per_call,
        "logical_byte_visits": logical_byte_visits,
    }


def run_one_with_workload_timeout(
    implementation: str,
    identity: Mapping[str, Any],
    cell: Mapping[str, Any],
    cpu: int,
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
    failure_output: Path,
) -> dict[str, Any]:
    """Run one child with the K65 campaign's attested workload budget."""
    executable = Path(str(identity["path"]))
    BASE.require(BASE.sha256(executable) == identity["sha256"],
                 f"{implementation} binary changed before execution")
    command = BASE.benchmark_command(
        implementation, executable, cell, cpu, iterations, warmup)
    timeout_budget = child_timeout_budget(
        implementation, cell, iterations, warmup)
    started_ns = time.monotonic_ns()
    try:
        completed = subprocess.run(
            command, env=BASE.CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout_budget["timeout_seconds"], check=False)
    except subprocess.TimeoutExpired as error:
        retained = BASE.retain_failure_streams(
            failure_output, implementation, started_ns,
            error.stdout, error.stderr)
        raise BASE.EvidenceError(
            f"{implementation} timed out after "
            f"{timeout_budget['timeout_seconds']} workload-sized seconds; "
            f"retained output: {retained}") from error
    elapsed_ns = time.monotonic_ns() - started_ns
    if completed.returncode != 0:
        retained = BASE.retain_failure_streams(
            failure_output, implementation, started_ns,
            completed.stdout, completed.stderr)
        raise BASE.EvidenceError(
            f"{implementation} failed: " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:] +
            f"; retained output: {retained}")
    try:
        BASE.require(BASE.sha256(executable) == identity["sha256"],
                     f"{implementation} binary changed after execution")
        result = json.loads(completed.stdout.decode("utf-8"))
        normalized = BASE.validate_result(
            implementation, result, cell, source_commit, source_tree,
            iterations, warmup)
    except Exception as error:
        retained = BASE.retain_failure_streams(
            failure_output, implementation, started_ns,
            completed.stdout, completed.stderr)
        raise BASE.EvidenceError(
            f"{implementation} output was rejected: {error}; "
            f"retained output: {retained}") from error
    return {
        "implementation": implementation,
        "command": command,
        "timeout_budget": timeout_budget,
        "elapsed_ns": elapsed_ns,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "normalized": normalized,
        "result": result,
    }


# Preserve K66's before/after same-inode checks and replace only its inherited
# fixed-60-second inner launcher.
SUPPORT._BASE_RUN_ONE = run_one_with_workload_timeout


_BASE_VALIDATE_RESULT = BASE.validate_result


def _require_exact_keys(
    value: Mapping[str, Any], expected: Mapping[str, Any], label: str,
) -> None:
    for key, expected_value in expected.items():
        BASE.require(key in value and type(value[key]) is type(expected_value)
                     and value[key] == expected_value,
                     f"{label} {key} differs from the v31 contract")


def _require_v31_payload_key_sets(result: Mapping[str, Any]) -> None:
    BASE.require(set(result) == {
        "schema", "build", "parameters", "resolved", "correctness",
        "memory", "metrics", "legacy", "workload_digests",
    }, "v31 top-level payload keys changed")
    sections = {
        name: result.get(name)
        for name in (
            "build", "parameters", "resolved", "correctness", "memory",
            "metrics", "legacy", "workload_digests")
    }
    BASE.require(all(isinstance(value, Mapping)
                     for value in sections.values()),
                 "v31 payload section is not an object")

    expected_build = {
        "compiler", "compiler_version", "cplusplus",
        "k8r3r4_t4_terminal_diagnostic_disabled",
        "t8_full_parity_terminal_diagnostic_disabled",
        "k16r8_b256_terminal_diagnostic_disabled",
        "k9r5_b256_terminal_diagnostic_disabled",
        "k9r5_b1024_terminal_diagnostic_disabled",
        "k65r65_b64_packed_terminal_diagnostic_mode",
        "k65r65_b64_packed_terminal_diagnostic_disabled",
        "k65r65_b64_packed_terminal_mode_latched",
        "k65r65_b64_packed_terminal_selector_expected_selected",
        "k65r65_b64_packed_terminal_selector_selected",
        "k65r65_b64_packed_terminal_selector_contract",
        "k65r65_b64_packed_terminal_timed_ordinary_encode_api",
        "k65r65_b64_packed_terminal_timed_one_shot_encode_api",
        "source_commit", "source_tree", "source_tracked_dirty",
    }
    optional_build = "high_t16_q2_b64_fused_diagnostic_disabled"
    if optional_build in sections["build"]:
        expected_build.add(optional_build)
    BASE.require(set(sections["build"]) == expected_build,
                 "v31 build payload keys changed")
    BASE.require(set(sections["parameters"]) == {
        "K", "R", "requested_profile", "requested_field",
        "requested_backend", "force_generic_decode",
        "force_specialized_decode", "force_tiled_decode",
        "force_materialized_decode", "skip_legacy", "retain_samples",
        "attest_source", "measure_one_shot_encode",
        "k65r65_b64_packed_terminal_mode", "shard_bytes", "loss_count",
        "missing_original_indices", "batch", "reuse", "iterations",
        "warmup", "thread_count", "seed",
    }, "v31 parameter payload keys changed")
    BASE.require(set(sections["resolved"]) == {
        "profile", "field", "backend", "thread_count", "parent_count",
        "padded_side",
    }, "v31 resolved payload keys changed")
    BASE.require(set(sections["correctness"]) == {
        "leopard2_round_trip", "legacy_comparison",
    }, "v31 correctness payload keys changed")
    BASE.require(set(sections["memory"]) == {
        "scratch_alignment", "encode_scratch_bytes_per_stripe",
        "decode_scratch_bytes_per_stripe", "encode_scratch_bytes_batch",
        "decode_scratch_bytes_batch",
    }, "v31 memory payload keys changed")
    BASE.require(set(sections["metrics"]) == {
        "codec_setup", "encode_execution", "one_shot_encode",
        "decode_plan_setup", "decode_execution",
        "decode_amortized_at_reuse", "rate_semantics",
    }, "v31 metrics payload keys changed")
    BASE.require(set(sections["legacy"]) == {
        "available", "unavailable_reason", "codec_setup",
        "decode_timing_includes_setup", "encode_execution",
        "decode_including_setup",
    }, "v31 legacy payload keys changed")
    BASE.require(set(sections["workload_digests"]) == {
        "algorithm", "original_data", "transmitted_parity",
        "recovered_originals",
    }, "v31 workload-digest payload keys changed")


def _ceil_power_of_two(value: int) -> int:
    BASE.require(type(value) is int and value > 0,
                 "codec dimension is not positive")
    return 1 << (value - 1).bit_length()


def expected_codec_geometry(cell: Mapping[str, Any]) -> tuple[int, int]:
    recovery_side = _ceil_power_of_two(int(cell["R"]))
    parent_count = _ceil_power_of_two(recovery_side + int(cell["K"]))
    return parent_count, recovery_side


def expected_missing_original_indices(
    cell: Mapping[str, Any],
) -> list[int]:
    """Reproduce benchmark.cpp SelectLosses and its uint64 xorshift."""
    original_count = int(cell["K"])
    loss_count = int(cell["loss"])
    mask = (1 << 64) - 1
    state = (int(cell["seed"]) ^ 0xD1B54A32D192ED03) & mask
    if state == 0:
        state = 0x9E3779B97F4A7C15
    order = list(range(original_count))
    for remaining in range(original_count, 1, -1):
        value = state
        value ^= (value << 13) & mask
        value ^= value >> 7
        value ^= (value << 17) & mask
        state = value & mask
        selected = state % remaining
        order[remaining - 1], order[selected] = (
            order[selected], order[remaining - 1])
    return sorted(order[:loss_count])


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    normalized = _BASE_VALIDATE_RESULT(
        implementation, result, cell, source_commit, source_tree,
        iterations, warmup)
    if implementation == "main":
        BASE.require(type(result) is dict,
                     "exact-main child payload is not an object")
        build = result.get("build")
        BASE.require(
            result.get("schema") == AUTHORITY_RECORD.BENCHMARK_SCHEMA and
            type(build) is dict and set(build) == {
                "main_source_commit", "pure_avx2", "cplusplus"} and
            build.get("main_source_commit") == EXACT_MAIN_SOURCE_COMMIT and
            build.get("pure_avx2") is True and
            type(build.get("cplusplus")) is int and build["cplusplus"] > 0,
            "exact-main child is not the sealed pure-AVX2 authority build")
        normalized["exact_main_authority_attribution"] = {
            "pure_avx2": True,
            "main_source_commit": EXACT_MAIN_SOURCE_COMMIT,
            "authority_record_sha256": EXACT_MAIN_AUTHORITY_RECORD_SHA256,
            "executable_sha256": EXACT_MAIN_EXECUTABLE_SHA256,
        }
        return normalized
    BASE.require(isinstance(result, Mapping),
                 "Leopard2 v31 result is not an object")
    _require_v31_payload_key_sets(result)
    mode = 0 if implementation == "control" else 1
    exact_target_shape = (
        cell.get("K") == 65 and cell.get("R") == 65 and
        cell.get("bytes") == 64)
    expected_selected = mode == 1 and exact_target_shape
    parameters = result.get("parameters")
    build = result.get("build")
    resolved = result.get("resolved")
    BASE.require(isinstance(parameters, Mapping) and
                 isinstance(build, Mapping) and
                 isinstance(resolved, Mapping),
                 "Leopard2 v31 parameters, resolved geometry, or build "
                 "markers are absent")
    expected_missing = expected_missing_original_indices(cell)
    observed_missing = parameters.get("missing_original_indices")
    BASE.require(
        type(observed_missing) is list and
        all(type(index) is int for index in observed_missing) and
        observed_missing == expected_missing,
        "parameter missing_original_indices differs in type or value from "
        "the deterministic selection")
    _require_exact_keys(parameters, {
        "K": int(cell["K"]),
        "R": int(cell["R"]),
        "requested_profile": "legacy_high_v1",
        "requested_field": "gf8",
        "requested_backend": "avx2",
        "force_generic_decode": False,
        "force_specialized_decode": False,
        "force_tiled_decode": False,
        "force_materialized_decode": False,
        "skip_legacy": True,
        "retain_samples": True,
        "attest_source": True,
        "measure_one_shot_encode": True,
        "k65r65_b64_packed_terminal_mode": mode,
        "shard_bytes": int(cell["bytes"]),
        "loss_count": int(cell["loss"]),
        "missing_original_indices": expected_missing,
        "batch": int(cell["batch"]),
        "reuse": int(cell["reuse"]),
        "iterations": iterations,
        "warmup": warmup,
        "thread_count": 1,
        "seed": int(cell["seed"]),
    }, "parameter")
    _require_exact_keys(build, {
        "k8r3r4_t4_terminal_diagnostic_disabled": False,
        "t8_full_parity_terminal_diagnostic_disabled": False,
        "k16r8_b256_terminal_diagnostic_disabled": False,
        "k9r5_b256_terminal_diagnostic_disabled": False,
        "k9r5_b1024_terminal_diagnostic_disabled": False,
    }, "unrelated build marker")
    BASE.require(
        type(build.get("compiler")) is str and bool(build["compiler"]) and
        type(build.get("compiler_version")) is str and
        bool(build["compiler_version"]) and
        type(build.get("cplusplus")) is int and build["cplusplus"] > 0 and
        ("high_t16_q2_b64_fused_diagnostic_disabled" not in build or
         type(build["high_t16_q2_b64_fused_diagnostic_disabled"]) is bool),
        "v31 compiler or optional build marker is malformed")
    if "high_t16_q2_b64_fused_diagnostic_disabled" in build:
        BASE.require(
            build["high_t16_q2_b64_fused_diagnostic_disabled"] is False,
            "unrelated high-T16 fused disable marker is active")
    parent_count, padded_side = expected_codec_geometry(cell)
    _require_exact_keys(resolved, {
        "profile": "legacy_high_v1",
        "field": "gf8",
        "backend": "avx2",
        "thread_count": 1,
        "parent_count": parent_count,
        "padded_side": padded_side,
    }, "resolved geometry")
    expected_build = {
        "k65r65_b64_packed_terminal_diagnostic_mode": mode,
        "k65r65_b64_packed_terminal_diagnostic_disabled": mode == 0,
        "k65r65_b64_packed_terminal_mode_latched": mode,
        "k65r65_b64_packed_terminal_selector_expected_selected":
            expected_selected,
        "k65r65_b64_packed_terminal_selector_selected": expected_selected,
        "k65r65_b64_packed_terminal_selector_contract": SELECTOR_CONTRACT,
        "k65r65_b64_packed_terminal_timed_ordinary_encode_api":
            ORDINARY_API_MARKER,
    }
    if cell.get("measure_one_shot") is True:
        expected_build[
            "k65r65_b64_packed_terminal_timed_one_shot_encode_api"] = \
                ONE_SHOT_API_MARKER
    _require_exact_keys(build, expected_build, "build marker")
    normalized["k65r65_b64_packed_terminal_attribution"] = {
        "requested_mode": mode,
        "latched_mode": build[
            "k65r65_b64_packed_terminal_mode_latched"],
        "selector_expected_selected": expected_selected,
        "selector_selected": build[
            "k65r65_b64_packed_terminal_selector_selected"],
        "ordinary_api": build[
            "k65r65_b64_packed_terminal_timed_ordinary_encode_api"],
        "one_shot_api": build.get(
            "k65r65_b64_packed_terminal_timed_one_shot_encode_api"),
    }
    return normalized


BASE.validate_result = validate_result


# Bypass the inherited K9 generation's fail-fast neighbor gate.  Generation 2
# retains every structurally valid 25-round result and applies its frozen K65
# gates only in the independent final summary below.
_BASE_ANALYZE = PRODUCTION_SUPPORT._BASE_ANALYZE
BASE.require(callable(_BASE_ANALYZE),
             "unversioned structural analysis helper is unavailable")


def _validated_ci(
    ratio: Any, cell_id: str, metric_name: str,
) -> list[float]:
    BASE.require(isinstance(ratio, Mapping) and
                 isinstance(ratio.get("ci95"), list) and
                 len(ratio["ci95"]) == 2,
                 f"{cell_id} {metric_name} contrast is missing")
    lower, upper = ratio["ci95"]
    BASE.require(
        isinstance(lower, (int, float)) and not isinstance(lower, bool) and
        isinstance(upper, (int, float)) and not isinstance(upper, bool) and
        math.isfinite(float(lower)) and math.isfinite(float(upper)) and
        lower <= upper,
        f"{cell_id} {metric_name} confidence interval is malformed")
    return [float(lower), float(upper)]


def analyze(
    cell: Mapping[str, Any], rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    result = _BASE_ANALYZE(cell, rounds)
    if cell.get("role") == "target":
        ratios = {
            "encode_execution_control": (
                result.get("control_over_candidate"),
                BASE.TARGET_CONTROL_FLOOR),
            "one_shot_encode_control": (
                result.get("one_shot_control_over_candidate"),
                BASE.TARGET_CONTROL_FLOOR),
            "encode_execution_main": (
                result.get("main_over_candidate"), BASE.TARGET_MAIN_FLOOR),
            "one_shot_encode_main": (
                result.get("one_shot_main_over_candidate"),
                BASE.TARGET_MAIN_FLOOR),
        }
        lower_bounds = {}
        accepted = True
        for name, (ratio, floor) in ratios.items():
            lower = _validated_ci(ratio, "target", name)[0]
            lower_bounds[name] = lower
            accepted = accepted and lower >= floor
        result["target_acceptance_validation"] = {
            "control_floor": BASE.TARGET_CONTROL_FLOOR,
            "main_floor": BASE.TARGET_MAIN_FLOOR,
            "lower_ci95": lower_bounds,
            "accepted": accepted,
        }
        return result
    if cell.get("role") != "neighbor":
        return result

    cell_id = str(cell.get("id"))
    selector_cis = {
        "encode_execution": _validated_ci(
            result.get("control_over_candidate"), cell_id,
            "encode_execution"),
        "one_shot_encode": _validated_ci(
            result.get("one_shot_control_over_candidate"), cell_id,
            "one_shot_encode"),
    }
    selector_accepted = all(
        BASE.NEIGHBOR_EQUIVALENCE_LOWER <= interval[0] and
        interval[1] <= BASE.NEIGHBOR_EQUIVALENCE_UPPER
        for interval in selector_cis.values())
    result["neighbor_selector_inertness_validation"] = {
        "equivalence_band": [
            BASE.NEIGHBOR_EQUIVALENCE_LOWER,
            BASE.NEIGHBOR_EQUIVALENCE_UPPER,
        ],
        "ci95": selector_cis,
        "accepted": selector_accepted,
    }

    if not cell.get("compare_main"):
        return result
    ratios = {
        "encode_execution": result.get("main_over_candidate"),
        "one_shot_encode": result.get("one_shot_main_over_candidate"),
    }
    context: dict[str, Any] = {}
    for name, ratio in ratios.items():
        interval = _validated_ci(
            ratio, cell_id, f"{name} exact-main")
        context[name] = {
            "speedup": ratio.get("speedup"),
            "ci95": interval,
        }
    floor = cell.get("main_floor")
    accepted = True
    if floor is not None:
        BASE.require(type(floor) is float and floor == BASE.RETAINED_MAIN_FLOOR,
                     f"{cell_id} exact-main floor is not established")
        accepted = all(item["ci95"][0] >= floor for item in context.values())
    result["exact_main_context"] = {
        "policy": "gated" if floor is not None else "context_only",
        "floor": floor,
        "metrics": context,
        "accepted": accepted,
    }
    return result


BASE.analyze = analyze


def _validated_ratio_record(value: Any, label: str) -> dict[str, Any]:
    expected_keys = {
        "speedup", "ci95", "round_log_ratios", "confidence_level",
        "degrees_of_freedom", "t_critical",
    }
    BASE.require(type(value) is dict and set(value) == expected_keys,
                 f"{label} confidence record keys changed")
    samples = value.get("round_log_ratios")
    BASE.require(
        type(samples) is list and len(samples) == CONFIRMATORY_ROUNDS and
        all(isinstance(item, (int, float)) and not isinstance(item, bool) and
            math.isfinite(float(item)) for item in samples),
        f"{label} does not retain 25 finite round contrasts")
    recomputed = confidence_interval(samples)
    BASE.require(value == recomputed,
                 f"{label} differs from its retained round contrasts")
    return copy.deepcopy(value)


def _acceptance_failures(
    analyses: Any,
) -> list[dict[str, Any]]:
    expected_cells = cells()
    BASE.require(
        type(analyses) is list and len(analyses) == len(expected_cells) and
        [item.get("cell") if type(item) is dict else None
         for item in analyses] == expected_cells,
        "final summary does not contain the frozen 13-cell matrix")
    failures: list[dict[str, Any]] = []
    for analysis, cell in zip(analyses, expected_cells):
        cell_id = cell["id"]
        control = _validated_ratio_record(
            analysis.get("control_over_candidate"),
            f"{cell_id} batch control")
        one_control = _validated_ratio_record(
            analysis.get("one_shot_control_over_candidate"),
            f"{cell_id} one-shot control")
        if cell["role"] == "target":
            main = _validated_ratio_record(
                analysis.get("main_over_candidate"),
                f"{cell_id} batch main")
            one_main = _validated_ratio_record(
                analysis.get("one_shot_main_over_candidate"),
                f"{cell_id} one-shot main")
            lower_ci95 = {
                "encode_execution_control": control["ci95"][0],
                "one_shot_encode_control": one_control["ci95"][0],
                "encode_execution_main": main["ci95"][0],
                "one_shot_encode_main": one_main["ci95"][0],
            }
            accepted = (
                lower_ci95["encode_execution_control"] >=
                    BASE.TARGET_CONTROL_FLOOR and
                lower_ci95["one_shot_encode_control"] >=
                    BASE.TARGET_CONTROL_FLOOR and
                lower_ci95["encode_execution_main"] >=
                    BASE.TARGET_MAIN_FLOOR and
                lower_ci95["one_shot_encode_main"] >=
                    BASE.TARGET_MAIN_FLOOR)
            BASE.require(
                analysis.get("target_acceptance_validation") == {
                    "control_floor": BASE.TARGET_CONTROL_FLOOR,
                    "main_floor": BASE.TARGET_MAIN_FLOOR,
                    "lower_ci95": lower_ci95,
                    "accepted": accepted,
                },
                "target acceptance projection changed")
            if not accepted:
                failures.append({
                    "cell_id": cell_id,
                    "gate": "target-control-and-main-ci95-lower-floors",
                    "observed_lower_ci95": lower_ci95,
                })
            continue

        selector_cis = {
            "encode_execution": list(control["ci95"]),
            "one_shot_encode": list(one_control["ci95"]),
        }
        selector_accepted = all(
            BASE.NEIGHBOR_EQUIVALENCE_LOWER <= interval[0] and
            interval[1] <= BASE.NEIGHBOR_EQUIVALENCE_UPPER
            for interval in selector_cis.values())
        BASE.require(
            analysis.get("neighbor_selector_inertness_validation") == {
                "equivalence_band": [
                    BASE.NEIGHBOR_EQUIVALENCE_LOWER,
                    BASE.NEIGHBOR_EQUIVALENCE_UPPER,
                ],
                "ci95": selector_cis,
                "accepted": selector_accepted,
            },
            f"{cell_id} selector-inertness projection changed")
        if not selector_accepted:
            failures.append({
                "cell_id": cell_id,
                "gate": "neighbor-selector-ci95-equivalence-band",
                "observed_ci95": selector_cis,
            })
        if not cell["compare_main"]:
            BASE.require(
                "main_over_candidate" not in analysis and
                "one_shot_main_over_candidate" not in analysis and
                "exact_main_context" not in analysis,
                f"{cell_id} acquired an undeclared exact-main contrast")
            continue
        main = _validated_ratio_record(
            analysis.get("main_over_candidate"),
            f"{cell_id} batch main")
        one_main = _validated_ratio_record(
            analysis.get("one_shot_main_over_candidate"),
            f"{cell_id} one-shot main")
        floor = cell["main_floor"]
        main_accepted = floor is None or (
            main["ci95"][0] >= floor and one_main["ci95"][0] >= floor)
        expected_context = {
            "policy": "gated" if floor is not None else "context_only",
            "floor": floor,
            "metrics": {
                "encode_execution": {
                    "speedup": main["speedup"],
                    "ci95": list(main["ci95"]),
                },
                "one_shot_encode": {
                    "speedup": one_main["speedup"],
                    "ci95": list(one_main["ci95"]),
                },
            },
            "accepted": main_accepted,
        }
        BASE.require(analysis.get("exact_main_context") == expected_context,
                     f"{cell_id} exact-main policy projection changed")
        if not main_accepted:
            failures.append({
                "cell_id": cell_id,
                "gate": "retained-exact-main-ci95-lower-floor",
                "floor": floor,
                "observed_ci95": {
                    "encode_execution": list(main["ci95"]),
                    "one_shot_encode": list(one_main["ci95"]),
                },
            })
    return failures


def _validate_complete_campaign(
    raw: Mapping[str, Any], summary: Mapping[str, Any],
    authority: Mapping[str, Any], preregistration: Mapping[str, Any],
) -> list[dict[str, Any]]:
    expected_cells = cells()
    BASE.require(
        raw.get("schema") == RAW_SCHEMA_V2 and
        raw.get("requested_rounds") == CONFIRMATORY_ROUNDS and
        raw.get("iterations") == FIXED_ITERATIONS and
        raw.get("warmup") == FIXED_WARMUP and
        raw.get("exact_main_authority") == authority and
        raw.get("attempt_preregistration") == preregistration and
        summary.get("schema") == PRELIMINARY_SUMMARY_SCHEMA_V2 and
        summary.get("generation") == GENERATION and
        summary.get("exact_main_authority") == authority and
        summary.get("attempt_preregistration") == preregistration and
        summary.get("cell_count") == len(expected_cells) and
        _canonical_sha256(expected_cells) == MATRIX_SHA256,
        "preliminary campaign policy or authority binding changed")
    raw_cells = raw.get("cells")
    analyses = summary.get("cells")
    BASE.require(
        type(raw_cells) is list and len(raw_cells) == len(expected_cells) and
        type(analyses) is list and len(analyses) == len(expected_cells),
        "campaign is missing one or more frozen cells")
    process_count = 0
    discarded_process_count = 0
    discarded_attempt_count = 0
    for index, (raw_cell, expected_cell, analysis) in enumerate(zip(
            raw_cells, expected_cells, analyses)):
        BASE.require(
            type(raw_cell) is dict and
            raw_cell.get("cell") == expected_cell and
            "failed_round" not in raw_cell and
            type(raw_cell.get("rounds")) is list and
            len(raw_cell["rounds"]) == CONFIRMATORY_ROUNDS,
            f"cell {index} is not one complete 25-round campaign")
        for round_index, round_value in enumerate(raw_cell["rounds"]):
            discarded = round_value.get("discarded_attempts") \
                if type(round_value) is dict else None
            BASE.require(
                type(round_value) is dict and
                round_value.get("round") == round_index and
                type(round_value.get("isolation")) is dict and
                round_value["isolation"].get("accepted") is True and
                type(round_value.get("invocations")) is list and
                type(discarded) is list and
                len(discarded) <= BASE.MAX_ISOLATION_ATTEMPTS - 1 and
                all(type(item) is dict and
                    type(item.get("isolation")) is dict and
                    item["isolation"].get("accepted") is False
                    for item in discarded),
                f"cell {index} round {round_index} is incomplete or contaminated")
            process_count += len(round_value["invocations"])
            discarded_process_count += sum(
                len(item.get("invocations", [])) for item in discarded)
            discarded_attempt_count += len(discarded)
        recomputed = analyze(expected_cell, raw_cell["rounds"])
        BASE.require(analysis == recomputed,
                     f"cell {index} summary differs from retained rounds")
    BASE.require(
        summary.get("process_count") == process_count and
        summary.get("discarded_process_count") == discarded_process_count and
        summary.get("discarded_round_attempts") == discarded_attempt_count and
        summary.get("all_digests_matched") is True and
        summary.get("all_rounds_zero_sibling_nonidle") is True and
        type(summary.get("binary_sha256")) is dict and
        summary["binary_sha256"].get("main") ==
            EXACT_MAIN_EXECUTABLE_SHA256 and
        summary["binary_sha256"].get("candidate") ==
            summary["binary_sha256"].get("control") and
        summary["binary_sha256"].get("candidate") !=
            summary["binary_sha256"].get("main"),
        "preliminary summary counts, isolation, or binary identities changed")
    return _acceptance_failures(analyses)


def _acceptance_contract(
    authority: Mapping[str, Any], preregistration: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": "leopard2-gf8-k65r65-b64-packed-terminal-acceptance/v2",
        "accepted_rounds_per_cell": CONFIRMATORY_ROUNDS,
        "inference_unit": "within-round log contrast",
        "student_t_degrees_of_freedom": CONFIRMATORY_ROUNDS - 1,
        "student_t_two_sided_critical": PRODUCTION_SUPPORT.T95_DF24,
        "iterations_per_process": FIXED_ITERATIONS,
        "warmup_iterations_per_process": FIXED_WARMUP,
        "reuse_per_process": FIXED_REUSE,
        "maximum_isolation_attempts_per_round":
            BASE.MAX_ISOLATION_ATTEMPTS,
        "campaign_attempt_budget": ATTEMPT_BUDGET,
        "discarded_isolation_attempts_used_in_inference": False,
        "pooling": "none",
        "trimming": "none",
        "matrix_selection": "all_13_predeclared_k65r65_b64_cells",
        "matrix_sha256": MATRIX_SHA256,
        "target_control_ci95_lower_floor": BASE.TARGET_CONTROL_FLOOR,
        "target_main_ci95_lower_floor": BASE.TARGET_MAIN_FLOOR,
        "neighbor_selector_ci95_equivalence_band": [
            BASE.NEIGHBOR_EQUIVALENCE_LOWER,
            BASE.NEIGHBOR_EQUIVALENCE_UPPER,
        ],
        "retained_exact_main_ci95_lower_floor": BASE.RETAINED_MAIN_FLOOR,
        "outcome_selected_cells": False,
        "statistical_campaign_retries": 0,
        "decision_evidence": "this standalone 25-round campaign only",
        "attempt": preregistration["attempt"],
        "exact_main_authority": copy.deepcopy(authority),
        "historical_references_non_authoritative": True,
    }


def apply_final_acceptance(
    preliminary: Mapping[str, Any], raw: Mapping[str, Any],
    authority: Mapping[str, Any], preregistration: Mapping[str, Any],
) -> dict[str, Any]:
    """Create the only promotable K65 v2 summary from complete evidence."""
    BASE.require(
        type(preliminary) is dict and
        preliminary.get("schema") == PRELIMINARY_SUMMARY_SCHEMA_V2 and
        "promotable" not in preliminary and
        "acceptance_contract" not in preliminary and
        "finalization_failure" not in preliminary,
        "input is not one unfinalized K65 v2 preliminary summary")
    binding = validate_exact_main_authority_binding(authority)
    attempt = validate_attempt_preregistration(
        preregistration.get("attempt") if type(preregistration) is dict else -1,
        preregistration.get("lineage") if type(preregistration) is dict else None,
        verify_files=False)
    BASE.require(attempt == preregistration,
                 "attempt preregistration block changed")
    failures = _validate_complete_campaign(
        raw, preliminary, binding, attempt)
    final = copy.deepcopy(preliminary)
    for name in (
        "target_control_failure", "target_main_failure",
        "target_control_failure_by_metric", "target_main_failure_by_metric",
        "credible_neighbor_regressions",
        "credible_neighbor_one_shot_regressions",
    ):
        final.pop(name, None)
    final["generation"] = GENERATION
    final["exact_main_authority"] = binding
    final["attempt_preregistration"] = attempt
    final["acceptance_contract"] = _acceptance_contract(binding, attempt)
    final["gate_failures"] = failures
    final["status"] = "accepted" if not failures else "rejected"
    final["promotable"] = not failures
    final["schema"] = FINAL_SUMMARY_SCHEMA_V2
    return final


def _walk_historical_hash_paths(
    value: Any, path: tuple[Any, ...] = (),
) -> list[tuple[Any, ...]]:
    found: list[tuple[Any, ...]] = []
    if type(value) is dict:
        for name, child in value.items():
            found.extend(_walk_historical_hash_paths(child, path + (name,)))
    elif type(value) is list:
        for index, child in enumerate(value):
            found.extend(_walk_historical_hash_paths(child, path + (index,)))
    elif value in {HISTORICAL_EXECUTABLE_SHA256, HISTORICAL_ARCHIVE_SHA256}:
        found.append(path)
    return found


def validate_promotable_final_summary(value: Any) -> dict[str, Any]:
    """Accept only a complete standalone accepted v2 final summary."""
    expected_process_count = sum(
        len((BASE.TARGET_ORDER if cell["compare_main"] else
             BASE.NEIGHBOR_ORDER)[0]) * CONFIRMATORY_ROUNDS
        for cell in cells())
    binary_sha256 = value.get("binary_sha256") \
        if type(value) is dict else None
    BASE.require(
        type(value) is dict and
        set(value) == FINAL_SUMMARY_KEYS and
        value.get("schema") == FINAL_SUMMARY_SCHEMA_V2 and
        value.get("status") == "accepted" and
        value.get("promotable") is True and
        value.get("generation") == GENERATION and
        value.get("cell_count") == len(cells()) and
        value.get("process_count") == expected_process_count and
        type(value.get("discarded_process_count")) is int and
        value["discarded_process_count"] >= 0 and
        type(value.get("discarded_round_attempts")) is int and
        0 <= value["discarded_round_attempts"] <=
            len(cells()) * CONFIRMATORY_ROUNDS *
                (BASE.MAX_ISOLATION_ATTEMPTS - 1) and
        len(BASE.NEIGHBOR_ORDER[0]) *
            value["discarded_round_attempts"] <=
            value["discarded_process_count"] <=
            len(BASE.TARGET_ORDER[0]) *
                value["discarded_round_attempts"] and
        (value["discarded_process_count"] -
         len(BASE.NEIGHBOR_ORDER[0]) *
             value["discarded_round_attempts"]) %
            (len(BASE.TARGET_ORDER[0]) -
             len(BASE.NEIGHBOR_ORDER[0])) == 0 and
        value.get("all_digests_matched") is True and
        value.get("all_rounds_zero_sibling_nonidle") is True and
        type(binary_sha256) is dict and
        set(binary_sha256) == {"candidate", "control", "main"} and
        all(type(item) is str and
            re.fullmatch(r"[0-9a-f]{64}", item) is not None
            for item in binary_sha256.values()) and
        binary_sha256["candidate"] == binary_sha256["control"] and
        binary_sha256["candidate"] != binary_sha256["main"] and
        binary_sha256["main"] == EXACT_MAIN_EXECUTABLE_SHA256 and
        type(value.get("candidate_control_executable_sections_sha256"))
            is str and
        re.fullmatch(
            r"[0-9a-f]{64}",
            value["candidate_control_executable_sections_sha256"])
            is not None and
        type(value.get("source_commit")) is str and
        re.fullmatch(r"[0-9a-f]{40}", value["source_commit"])
            is not None and
        type(value.get("source_tree")) is str and
        re.fullmatch(r"[0-9a-f]{40}", value["source_tree"])
            is not None and
        value.get("main_commit") == EXACT_MAIN_SOURCE_COMMIT and
        type(value.get("raw_sha256")) is str and
        re.fullmatch(r"[0-9a-f]{64}", value["raw_sha256"])
            is not None and
        type(value.get("preliminary_summary_sha256")) is str and
        re.fullmatch(
            r"[0-9a-f]{64}", value["preliminary_summary_sha256"])
            is not None and
        type(value.get("finalization_cpu")) is int and
        value["finalization_cpu"] >= 0 and
        value.get("preliminary_generic_status") in {"accepted", "rejected"} and
        type(value.get("mode_words")) is dict and
        value["mode_words"].get("shared_binary_default", {}).get("value")
            == 1 and
        value.get("gate_failures") == [] and
        "finalization_failure" not in value,
        "summary is not an accepted promotable K65 v2 final")
    authority = validate_exact_main_authority_binding(
        value.get("exact_main_authority"))
    preregistration = value.get("attempt_preregistration")
    attempt = validate_attempt_preregistration(
        preregistration.get("attempt")
            if type(preregistration) is dict else -1,
        preregistration.get("lineage")
            if type(preregistration) is dict else None,
        verify_files=False)
    BASE.require(
        preregistration == attempt and
        value.get("acceptance_contract") ==
            _acceptance_contract(authority, attempt) and
        not _acceptance_failures(value.get("cells")),
        "promotable final summary policy projection changed")
    historical_paths = _walk_historical_hash_paths(value)
    BASE.require(
        historical_paths and all(
            len(path) >= 3 and path[-3] == "historical_non_authority" and
            path[-1] == "value" for path in historical_paths),
        "historical exact-main hashes escaped their non-authority disclosure")
    return copy.deepcopy(value)


def _read_json_object(path: Path, label: str) -> dict[str, Any]:
    value = _strict_json_path(path, label)
    BASE.require(type(value) is dict, f"{label} is not one JSON object")
    return value


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _replace_json(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(path.name + ".finalizing")
    BASE.require(not temporary.exists() and not temporary.is_symlink(),
                 f"stale finalization file exists: {temporary}")
    BASE.T8_SUPPORT.write_exclusive(temporary, value)
    os.replace(temporary, path)
    _fsync_directory(path.parent)


def _begin_isolated_finalization(
    benchmark_cpu: int, sibling_cpu: int,
) -> tuple[int, int]:
    finalization_cpu = None
    for candidate in range(os.cpu_count() or 0):
        if candidate in {benchmark_cpu, sibling_cpu}:
            continue
        try:
            os.sched_setaffinity(0, {candidate})
        except OSError:
            continue
        if os.sched_getaffinity(0) == {candidate}:
            finalization_cpu = candidate
            break
    BASE.require(finalization_cpu is not None,
                 "no CPU outside the timed pair is available for finalization")
    descriptor = os.open(BASE.LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
    except Exception:
        os.close(descriptor)
        raise
    return descriptor, int(finalization_cpu)


def validate_final_live_identities(
    raw: Mapping[str, Any], options: argparse.Namespace,
) -> None:
    BASE.require(
        BASE.T8_SUPPORT.file_identity(BASE.RUNNER_PATH) ==
            raw.get("runner_after") and
        [BASE.support_file_identity(path)
         for path in BASE.RUNNER_DEPENDENCIES] ==
            raw.get("runner_dependencies_after"),
        "runner or support identity changed before finalization")
    for collection_name in ("identities_after", "input_identities_after"):
        collection = raw.get(collection_name)
        BASE.require(type(collection) is dict,
                     f"{collection_name} is absent")
        for identity in collection.values():
            BASE.require(
                type(identity) is dict and
                BASE.T8_SUPPORT.file_identity(Path(str(identity["path"]))) ==
                    identity,
                f"{collection_name} changed before finalization")
    BASE.require(
        BASE.build_closure_identity(options) == raw.get("build_closure_after"),
        "production build closure changed before finalization")


def _bind_failure_journal(
    output: Path, authority: Mapping[str, Any],
    preregistration: Mapping[str, Any],
) -> bool:
    failure_path = output / "failure.json"
    if not failure_path.exists():
        return False
    failure = _read_json_object(failure_path, "K65 campaign failure")
    failure["generation"] = GENERATION
    failure["exact_main_authority"] = copy.deepcopy(authority)
    failure["attempt_preregistration"] = copy.deepcopy(preregistration)
    _replace_json(failure_path, failure)
    return True


def _write_finalization_failure(
    output: Path, error: BaseException,
    authority: Mapping[str, Any], preregistration: Mapping[str, Any],
) -> None:
    failure_path = output / "failure.json"
    BASE.require(not failure_path.exists() and not failure_path.is_symlink(),
                 "finalization failure collided with another failure record")
    raw_path = output / "raw.json"
    value = {
        "schema": "leopard2-gf8-k65r65-b64-packed-terminal-finalization-failure/v2",
        "generation": GENERATION,
        "failed_utc": BASE.MAIN_SUPPORT.utc_now(),
        "failure": {
            "type": type(error).__name__,
            "message": str(error)[:4096],
        },
        "raw_sha256": _sha256_file(raw_path) if raw_path.is_file() else None,
        "exact_main_authority": copy.deepcopy(authority),
        "attempt_preregistration": copy.deepcopy(preregistration),
    }
    BASE.T8_SUPPORT.write_exclusive(failure_path, value)
    _fsync_directory(output)


def _write_attempt_terminal(
    options: argparse.Namespace, outcome: str,
) -> dict[str, Any]:
    root = options.output.resolve(strict=True)
    summary_path = root / "summary.json"
    raw_path = root / "raw.json"
    failure_path = root / "failure.json"
    summary = _read_json_object(summary_path, "K65 final summary") \
        if summary_path.is_file() else None
    value = {
        "schema": ATTEMPT_TERMINAL_SCHEMA,
        "generation": GENERATION,
        "attempt": options.attempt,
        "attempt_budget": ATTEMPT_BUDGET,
        "outcome": outcome,
        "promotable": outcome == "accepted",
        "output_root": str(root),
        "summary_schema": summary["schema"] if summary is not None else "absent",
        "summary_sha256": _sha256_file(summary_path)
            if summary is not None else None,
        "raw_sha256": _sha256_file(raw_path) if raw_path.is_file() else None,
        "failure_sha256": _sha256_file(failure_path)
            if failure_path.is_file() else None,
        "authority_record_sha256": EXACT_MAIN_AUTHORITY_RECORD_SHA256,
    }
    validate_attempt_terminal(
        value, verify_files=False, allow_accepted=True)
    terminal_path = root / "attempt-terminal.json"
    BASE.require(not terminal_path.exists() and not terminal_path.is_symlink(),
                 "attempt terminal already exists")
    BASE.T8_SUPPORT.write_exclusive(terminal_path, value)
    _fsync_directory(root)
    return value


def main() -> int:
    """Collect preliminary evidence, then independently finalize generation 2."""
    options = BASE.parse_arguments()
    captured_stdout = io.StringIO()
    with contextlib.redirect_stdout(captured_stdout):
        preliminary_status = BASE.main()
    if not options.output.exists():
        output = captured_stdout.getvalue()
        if output:
            print(output, end="")
        return preliminary_status if preliminary_status != 0 else 1

    authority = options.exact_main_authority
    preregistration = options.attempt_preregistration
    if preliminary_status not in (0, 2):
        if not _bind_failure_journal(
                options.output, authority, preregistration):
            error = BASE.EvidenceError(
                "preliminary collection failed without a failure journal")
            _write_finalization_failure(
                options.output, error, authority, preregistration)
        _write_attempt_terminal(options, "failed")
        output = captured_stdout.getvalue()
        if output:
            print(output, end="")
        return preliminary_status

    try:
        finalization_lock, finalization_cpu = _begin_isolated_finalization(
            options.cpu, options.sibling)
    except Exception as error:
        summary_path = options.output / "summary.json"
        if summary_path.is_file():
            failed_summary = _read_json_object(
                summary_path, "K65 preliminary summary")
            failed_summary["schema"] = FINAL_SUMMARY_SCHEMA_V2
            failed_summary["generation"] = GENERATION
            failed_summary["status"] = "rejected"
            failed_summary["promotable"] = False
            failed_summary["exact_main_authority"] = copy.deepcopy(authority)
            failed_summary["attempt_preregistration"] = copy.deepcopy(
                preregistration)
            failed_summary["finalization_failure"] = {
                "type": type(error).__name__,
                "message": str(error)[:4096],
            }
            _replace_json(summary_path, failed_summary)
        _write_finalization_failure(
            options.output, error, authority, preregistration)
        _write_attempt_terminal(options, "failed")
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    try:
        summary_path = options.output / "summary.json"
        raw_path = options.output / "raw.json"
        if not summary_path.is_file() or not raw_path.is_file():
            error = BASE.EvidenceError(
                "preliminary raw evidence or summary is absent")
            _write_finalization_failure(
                options.output, error, authority, preregistration)
            _write_attempt_terminal(options, "failed")
            print(f"evidence rejected: {error}", file=sys.stderr)
            return 1

        preliminary = _read_json_object(
            summary_path, "K65 preliminary summary")
        raw = _read_json_object(raw_path, "K65 raw evidence")
        original_raw_sha256 = _sha256_file(raw_path)
        finalization_failed = False
        try:
            BASE.require(
                preliminary.get("schema") ==
                    PRELIMINARY_SUMMARY_SCHEMA_V2 and
                preliminary.get("raw_sha256") == original_raw_sha256,
                "preliminary summary does not bind raw evidence")
            raw["generation"] = GENERATION
            raw["exact_main_authority"] = copy.deepcopy(authority)
            raw["attempt_preregistration"] = copy.deepcopy(preregistration)
            _replace_json(raw_path, raw)
            preliminary["generation"] = GENERATION
            preliminary["exact_main_authority"] = copy.deepcopy(authority)
            preliminary["attempt_preregistration"] = copy.deepcopy(
                preregistration)
            preliminary["raw_sha256"] = _sha256_file(raw_path)
            preliminary["preliminary_generic_status"] = preliminary.get(
                "status")
            preliminary["finalization_cpu"] = finalization_cpu
            _replace_json(summary_path, preliminary)
            retained_preliminary = options.output / "preliminary-summary.json"
            BASE.require(
                not retained_preliminary.exists() and
                not retained_preliminary.is_symlink(),
                "retained preliminary summary already exists")
            BASE.T8_SUPPORT.write_exclusive(retained_preliminary, preliminary)
            _fsync_directory(options.output)

            validate_final_live_identities(raw, options)
            rebound = bind_exact_main_authority(
                options.exact_main_authority_lane)
            BASE.require(rebound == authority,
                         "exact-main authority changed before finalization")
            final = copy.deepcopy(preliminary)
            final["preliminary_summary_sha256"] = _sha256_file(
                retained_preliminary)
            final = apply_final_acceptance(
                final, raw, rebound, preregistration)
            if final["status"] == "accepted":
                validate_promotable_final_summary(final)
            _replace_json(summary_path, final)
        except Exception as error:
            finalization_failed = True
            failed_summary = copy.deepcopy(preliminary)
            failed_summary["schema"] = FINAL_SUMMARY_SCHEMA_V2
            failed_summary["generation"] = GENERATION
            failed_summary["status"] = "rejected"
            failed_summary["promotable"] = False
            failed_summary["exact_main_authority"] = copy.deepcopy(authority)
            failed_summary["attempt_preregistration"] = copy.deepcopy(
                preregistration)
            failed_summary["finalization_failure"] = {
                "type": type(error).__name__,
                "message": str(error)[:4096],
            }
            _replace_json(summary_path, failed_summary)
            _write_finalization_failure(
                options.output, error, authority, preregistration)

        final_summary = _read_json_object(summary_path, "K65 final summary")
        outcome = "failed" if finalization_failed else final_summary["status"]
        _write_attempt_terminal(options, outcome)
        print(json.dumps({
            "status": final_summary["status"],
            "promotable": final_summary.get("promotable", False),
            "generation": GENERATION,
            "attempt": options.attempt,
            "cells": final_summary.get("cell_count"),
            "processes": final_summary.get("process_count"),
            "discarded_processes": final_summary.get(
                "discarded_process_count"),
            "gate_failures": len(final_summary.get("gate_failures", [])),
        }, sort_keys=True))
        if finalization_failed:
            return 1
        if final_summary["status"] == "accepted":
            return 0
        return 2
    finally:
        os.close(finalization_lock)


if __name__ == "__main__":
    raise SystemExit(main())
