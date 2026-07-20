#!/usr/bin/env python3
"""Generate and advance the exact-main-first balanced decoder matrix."""

from __future__ import annotations

import argparse
import base64
import binascii
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
from typing import Any, Iterable

import balanced_evidence_common as common


PLAN_SCHEMA = "leopard2-balanced-promotion-plan/v2"
MAIN_CELL_SCHEMA = "leopard2-main-compare-cell-list/v2"
SURVIVOR_SCHEMA = "leopard2-balanced-promotion-survivors/v4"
STAGE_SCHEMA = "leopard2-balanced-promotion-stage/v5"
ATTESTATION_SCHEMA = "leopard2-balanced-auto-path-attestation/v5"
ATTESTATION_RESULT_SCHEMA = "leopard2-balanced-auto-path-result/v4"
EXACT_MANIFEST_SCHEMA = "leopard2-main-compare-manifest/v5"
EXACT_MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
EXACT_MAIN_TREE = "b7c8830d96a978f6ec14fe747095f066e351ae72"
# Exact commit payload used only by deterministic scope self-tests.  Retaining
# the signed commit bytes makes the fixture exercise the same object-ID/tree
# proof as real schema-v5 evidence rather than substituting a plausible hash.
EXACT_MAIN_COMMIT_BASE64 = (
    "dHJlZSBiN2M4ODMwZDk2YTk3OGY2ZWMxNGZlNzQ3MDk1ZjA2NmUzNTFhZTcyCnBhcmVudCAy"
    "MmRkYzc4MDQ5OThkMzFjOGYxYTI2MTdlZTcyMGUwNjNiMWZhNmNkCnBhcmVudCAzNjQyN2Rk"
    "MjViZjY2NWY0MTA1MjVhMDNhMWYwYTBlYTkyNzUxNTBiCmF1dGhvciBDaHJpcyBUYXlsb3Ig"
    "PG1yY2F0aWRAZ21haWwuY29tPiAxNzExMTkyODE2IC0wNzAwCmNvbW1pdHRlciBHaXRIdWIg"
    "PG5vcmVwbHlAZ2l0aHViLmNvbT4gMTcxMTE5MjgxNiAtMDcwMApncGdzaWcgLS0tLS1CRUdJ"
    "TiBQR1AgU0lHTkFUVVJFLS0tLS0KIAogd3NGY0JBQUJDQUFRQlFKbC9ycndDUkMxYVE3dXU1"
    "VWhsQUFBWUVVUUFKQkh1ZjBtUWRxQW9mclZTR3ZHT2ZRZwogSllQQVNZTmh5azhoUEs0N3E3"
    "REhPWVI0UUN3eHZySC9VaVpYcGMyZndSK3FlNUxSRkhqUjFaNE43K3AvK0J6cgogMnd6cXBZ"
    "NmR1SUt5NzVvdjFNR0RGa0FocmpPdlBXbDFCR2RlcTNCYnRQL0NzQ0JVSVEyOCtNZG5OVGRq"
    "dGprVwogMkJMWUFYN3VnU2pKWGQ1OEE0eTJ2MHJvUDZKeVRWbHR3b0NRcUtlYUlTVkNDa0Y0"
    "amcvcmJTZHRUck9WeE1oZAogenlPdmtTMzlzcGVOUk9UaUlsTTB4N3VESkR6ZFZja2draEV2"
    "N2dDN3Z6VWdzRFRHYUh6WXkwTU5ISDlmaHhBaAogcEdNQXEzV3p4ZGgxdjJLQTF5cExvZ1VV"
    "ekwwSzNFRFpvRjREL1RZMG96SDdESm1TUytzV1BVRWpUSHNXQzg5eAogdFR2bk8vSzhLNjZj"
    "SkZDQ1N1LzBvSm4yVFdhQnZUdldnNXJ5WFN2RlNibmVjNGRzQlAzczBZVS9odGJRbDN4Tgog"
    "TTkvMmw3M21ra2NJNHIxOVRHQU0zV0JzYTVxKzBrUHk5QnE2SHloeTM5ampOV2pkTnVzZFFG"
    "a3p0dXNSVFJsNwogYUx5SlNRTXB4cWVzeHRacWVEMjRmQmZYUmFRMlhhYkxFUytsTThBd0ly"
    "dytMUFFScWcrQ20rb3Rxb21BeG9QaQogemVmcUU3OGxXVkJRenVESkpLZnNHaERGY1BFaWJJ"
    "UXYyN0hOT3BsaTBtWkFJbGFnckJiOEF4cjRkbytVS0d6cAogMkhYYlprTjZFTWVGYVRKNzY4"
    "MUhLZHlQRFdiWUM3STVPRFBJNktkUW1kaWFJNmZldUJ4QTgybjU2R3laQ2w1OQogclFHZnFM"
    "dGU3VFVvRnVMTE53c1YKID00MTBSCiAtLS0tLUVORCBQR1AgU0lHTkFUVVJFLS0tLS0KIAoK"
    "TWVyZ2UgcHVsbCByZXF1ZXN0ICMyMyBmcm9tIGdibGV0cjQyL21hc3RlcgoKSGFuZGxlIHRo"
    "ZSBjYXNlIHdoZW4gbm8gb3JpZ2luYWwgZGF0YSB3YXMgbG9zdA==")
PROMOTION_GATE = 1.05
BOUNDARY_K = (
    5, 7, 8, 9, 14, 15, 16, 17, 29, 30, 31, 32, 33, 62, 63, 64,
    65, 96, 112, 120, 124, 125, 126, 127, 128,
)
EXACT_GATE_BYTES = (256, 4096, 65536)
ALIGNED_CONFIRM_BYTES = (192, 4032, 4096, 65536, 1024 * 1024)
TRUE_TAIL_BYTES = (193, 4033, 4097, 1024 * 1024 + 1)
BACKENDS = ("scalar", "ssse3", "avx2")
MODE_PAIRS = (("generic", "materialized"), ("generic", "tiled"))
CONFIRM_MODES = ("generic", "auto", "materialized", "tiled")
EXACT_ORDER = ("baseline", "candidate", "candidate", "baseline")
FORCED_ORDERS = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
# This promotion campaign has same-binary forced confirmation for exactly these
# backends.  AVX-512 and NEON are intentionally rejected rather than borrowing
# generic speed evidence collected under a kernel that the confirmation runner
# cannot force.  A later schema can add them with their own forced matrices.
CAMPAIGN_BACKENDS = ("scalar", "ssse3", "avx2")
AUTO_BACKENDS = CAMPAIGN_BACKENDS
EXACT_RAW_SCHEMA = "leopard2-main-compare-raw/v5"
MAX_GIT_COMMIT_BYTES = 1024 * 1024
MAX_RETAINED_TEXT_BYTES = 1024 * 1024
MAX_CPU_ID = 1_048_575
MAX_CPU_LIST_ENTRIES = 4096
BASELINE_LIBRARY_SOURCES = (
    "leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp",
)
CANDIDATE_LIBRARY_SOURCES = (
    "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp", "Leopard2Plan.cpp",
    "LeopardCommon.cpp", "LeopardFF8.cpp", "LeopardFF16.cpp",
    "Leopard2BackendSSSE3.cpp", "Leopard2BackendAVX2.cpp",
    "Leopard2BackendAVX512.cpp",
)
BASELINE_EXPECTED_COMPILE_COMMAND_COUNT = 5
CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT = 20
COMPILE_COMMANDS_SCHEMA = "leopard2-main-compare-compile-commands/v2"
BASELINE_COMPILE_PROFILE = \
    "gnu-compatible-cxx11-native-x86_64-release/v1"
CANDIDATE_COMPILE_PROFILE = \
    "gnu-compatible-cxx11-runtime-dispatch-x86_64-release/v1"
REQUIRED_CANDIDATE_CACHE = {
    "CMAKE_BUILD_TYPE": "Release",
    "ENABLE_OPENMP": "ON",
    "LEO2_BACKEND_VARIANT": "auto",
    "LEO2_BUILD_BENCHMARKS": "ON",
    "LEO2_BUILD_FUZZERS": "OFF",
    "LEO2_BUILD_TESTS": "OFF",
    "LEO2_ENABLE_CUDA": "OFF",
}
REQUIRED_BASELINE_CACHE = {
    "CMAKE_BUILD_TYPE": "Release",
    "LEO_MAIN_HAS_MARCH_NATIVE": "1",
}
EXCLUDED_CAMPAIGN_BACKENDS = {
    "avx512": "no forced AVX-512 runner coverage in this campaign schema",
    "neon": "no native NEON forced-runner coverage in this campaign schema",
}
# Exact string pairs emitted by DecodePathName/DecodePathRuleName in
# Leopard2Dispatch.h.  This list is deliberately complete: it prevents a
# negative predicate ("not generic") from accepting unknown or cross-paired
# strings while the exact AUTO prediction below narrows each attestation cell
# to the one pair its plan and selector geometry can produce.
PRODUCTION_DECODE_PATH_RULE_PAIRS = frozenset({
    ("no_op", "no_op"),
    ("direct", "direct"),
    ("generic", "forced_generic"),
    ("materialized", "forced_materialized"),
    ("tiled", "forced_tiled"),
    ("generic", "balanced_generic"),
    ("tiled", "measured_batch_tiled"),
    ("materialized", "measured_materialized"),
    ("tiled", "workspace_tiled"),
    ("materialized", "workspace_materialized"),
    ("materialized", "unsupported_profile"),
})
# Production direct-plan limits mirrored from leopard2.cpp.  The attestation
# source/build identity binds that file, so changing these values requires an
# explicit evidence-protocol update rather than silently broadening acceptance.
DIRECT_MAX_ORIGINALS = 16
DIRECT_MAX_LOSSES = 4


class PlanError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PlanError(message)


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            block = stream.read(1024 * 1024)
            if not block:
                return digest.hexdigest()
            digest.update(block)


def deterministic_seed(namespace: str, *values: int) -> int:
    seed = int.from_bytes(hashlib.sha256(canonical_bytes({
        "namespace": namespace, "values": values,
    })).digest()[:8], "big")
    return seed or 1


def signed(value: dict[str, Any]) -> dict[str, Any]:
    result = dict(value)
    require("content_sha256" not in result, "document is already signed")
    result["content_sha256"] = canonical_sha256(result)
    return result


def unsigned(value: object, schema: str, label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    result = dict(value)
    digest = result.pop("content_sha256", None)
    require(result.get("schema") == schema, f"{label} schema differs")
    require(isinstance(digest, str) and len(digest) == 64 and
            canonical_sha256(result) == digest, f"{label} digest differs")
    return result


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as output:
        json.dump(value, output, indent=2, sort_keys=True, allow_nan=False)
        output.write("\n")


def load_json(path: Path) -> object:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def main_cell(identifier: str, k: int, r: int, shard_bytes: int,
              losses: int, namespace: str) -> dict[str, Any]:
    require(re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", identifier) is not None,
            f"invalid exact-main identifier {identifier!r}")
    require(0 < r <= k and 0 < losses <= r and
            shard_bytes > 0 and shard_bytes % 64 == 0,
            f"invalid exact-main cell {identifier}")
    seed = deterministic_seed(namespace, k, r, shard_bytes, losses)
    return {
        "identifier": identifier, "K": k, "R": r,
        "shard_bytes": shard_bytes, "loss_count": losses, "seed": seed,
        "runner_cell": f"{identifier}:{k}:{r}:{shard_bytes}:{losses}:{seed}",
    }


def validate_main_cell(cell: object, label: str) -> dict[str, Any]:
    require(isinstance(cell, dict) and set(cell) == {
        "identifier", "K", "R", "shard_bytes", "loss_count", "seed",
        "runner_cell",
    }, f"{label} shape differs")
    identifier = cell["identifier"]
    require(isinstance(identifier, str) and
            re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", identifier) is not None,
            f"{label} identifier is invalid")
    for key in ("K", "R", "shard_bytes", "loss_count", "seed"):
        require(type(cell[key]) is int and cell[key] > 0,
                f"{label} {key} is invalid")
    require(cell["R"] <= cell["K"] and cell["loss_count"] <= cell["R"] and
            cell["shard_bytes"] % 64 == 0 and cell["seed"] < 1 << 64,
            f"{label} is outside exact-main limits")
    require(cell["runner_cell"] == (
        f"{identifier}:{cell['K']}:{cell['R']}:{cell['shard_bytes']}:"
        f"{cell['loss_count']}:{cell['seed']}"),
        f"{label} runner encoding differs")
    return cell


def main_document(name: str, mode: str, cells: list[dict[str, Any]],
                  purpose: str, run_condition: str) -> dict[str, Any]:
    require(mode in ("auto", "generic", "materialized", "tiled"),
            "invalid exact-main mode")
    return signed({
        "schema": MAIN_CELL_SCHEMA, "name": name, "candidate_mode": mode,
        "exact_main_commit": EXACT_MAIN_COMMIT, "reuse": 8,
        "iterations": 9, "warmup": 2, "purpose": purpose,
        "run_condition": run_condition, "cells": cells,
    })


def validate_main_document(value: object, label: str) -> dict[str, Any]:
    result = unsigned(value, MAIN_CELL_SCHEMA, label)
    require(set(result) == {
        "schema", "name", "candidate_mode", "exact_main_commit", "reuse",
        "iterations", "warmup", "purpose", "run_condition", "cells",
    }, f"{label} fields differ")
    require(result["candidate_mode"] in CONFIRM_MODES and
            result["exact_main_commit"] == EXACT_MAIN_COMMIT and
            result["reuse"] == 8 and result["iterations"] == 9 and
            result["warmup"] == 2 and isinstance(result["name"], str) and
            isinstance(result["purpose"], str) and
            isinstance(result["run_condition"], str),
            f"{label} semantics differ")
    cells = result["cells"]
    require(isinstance(cells, list) and cells, f"{label} has no cells")
    for index, cell in enumerate(cells):
        validate_main_cell(cell, f"{label} cell {index}")
    require(len({cell["identifier"] for cell in cells}) == len(cells),
            f"{label} identifiers are duplicated")
    return result


def transform_groups(counts: Iterable[int]) -> dict[int, list[int]]:
    result: dict[int, list[int]] = {}
    for k in sorted(set(counts)):
        result.setdefault(1 << (k - 1).bit_length(), []).append(k)
    return result


def artifact(path: Path, root: Path, case_count: int, child_count: int,
             kind: str) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(root)), "kind": kind,
        "case_count": case_count, "timed_child_count": child_count,
        "size": path.stat().st_size, "sha256": file_sha256(path),
    }


def generate_plan(root: Path) -> dict[str, Any]:
    require(not root.exists(), f"refusing to replace existing output {root}")
    root.mkdir(parents=True)
    artifacts = []
    for side, counts in sorted(transform_groups(BOUNDARY_K).items()):
        cells = [
            main_cell(f"gate-k{k}-b{shard_bytes}", k, k, shard_bytes, k,
                      "exact-main-generic-gate")
            for k in counts for shard_bytes in EXACT_GATE_BYTES
        ]
        path = root / "exact-main-gate" / f"t{side}.json"
        write_json(path, main_document(
            f"balanced-generic-exact-main-gate-t{side}", "generic", cells,
            "External gate for balanced full-loss generic decode.",
            "Run before any internal forced-path promotion comparison."))
        artifacts.append(artifact(path, root, len(cells), len(cells) * 12,
                                  "exact_main_gate"))
    plan = signed({
        "schema": PLAN_SCHEMA,
        "name": "balanced-current-source-exact-main-first",
        "exact_main_commit": EXACT_MAIN_COMMIT,
        "promotion_gate": {
            "metrics": ["decode_first_use", "decode_reuse_amortized"],
            "minimum_ci95_lower": PROMOTION_GATE,
            "neighbor_maximum_regression": 0.02,
        },
        "runner_ordering": {
            "exact_main": {
                "round_count": 3, "order_in_every_round": list(EXACT_ORDER),
                "labels": ["ABBA", "ABBA", "ABBA"],
            },
            "forced_same_binary": {
                "round_count": 3,
                "orders": [list(order) for order in FORCED_ORDERS],
                "labels": ["ABBA", "BAAB", "ABBA"],
            },
        },
        "workflow": {
            "select_command": (
                "select verifies every schema-v5 exact-main gate manifest and "
                "derives an evidence-bound survivor set"),
            "refinement": (
                "opposite outcomes separated by an unmeasured K emit exact "
                "runner_cell refinements; advance refuses until none remain"),
            "advance_command": (
                "advance emits only surviving K/byte forced comparisons, "
                "aligned exact-main confirmations, true tails, and rejection controls"),
            "promotion": (
                "AUTO selected-path attestation and exact-main confirmation are "
                "required before a dispatcher change"),
        },
        "artifacts": artifacts,
        "initial_case_count": sum(item["case_count"] for item in artifacts),
        "initial_timed_child_count": sum(
            item["timed_child_count"] for item in artifacts),
    })
    write_json(root / "plan.json", plan)
    return plan


def json_tree(root: Path) -> list[Path]:
    return sorted(path.relative_to(root) for path in root.rglob("*.json"))


def compare_trees(actual: Path, expected: Path, label: str) -> None:
    actual_files = json_tree(actual)
    expected_files = json_tree(expected)
    require(actual_files == expected_files, f"{label} file set differs")
    for relative in expected_files:
        require((actual / relative).read_bytes() == (expected / relative).read_bytes(),
                f"{label} canonical bytes differ: {relative}")


def structural_plan(root: Path) -> dict[str, Any]:
    plan = unsigned(load_json(root / "plan.json"), PLAN_SCHEMA, "plan")
    require(set(plan) == {
        "schema", "name", "exact_main_commit", "promotion_gate",
        "runner_ordering", "workflow", "artifacts", "initial_case_count",
        "initial_timed_child_count",
    }, "plan fields differ")
    require(plan["exact_main_commit"] == EXACT_MAIN_COMMIT and
            plan["promotion_gate"] == {
                "metrics": ["decode_first_use", "decode_reuse_amortized"],
                "minimum_ci95_lower": PROMOTION_GATE,
                "neighbor_maximum_regression": 0.02,
            }, "plan gate semantics differ")
    total_cases = 0
    total_children = 0
    for item in plan["artifacts"]:
        require(isinstance(item, dict) and set(item) == {
            "path", "kind", "case_count", "timed_child_count", "size", "sha256",
        } and item["kind"] == "exact_main_gate", "plan artifact differs")
        path = root / item["path"]
        require(path.is_file() and path.stat().st_size == item["size"] and
                file_sha256(path) == item["sha256"], "plan artifact bytes differ")
        document = validate_main_document(load_json(path), str(path))
        require(document["candidate_mode"] == "generic" and
                len(document["cells"]) == item["case_count"] and
                item["timed_child_count"] == item["case_count"] * 12,
                "gate artifact accounting differs")
        total_cases += item["case_count"]
        total_children += item["timed_child_count"]
    require(total_cases == plan["initial_case_count"] and
            total_children == plan["initial_timed_child_count"],
            "plan totals differ")
    return plan


def validate_plan(root: Path) -> dict[str, Any]:
    plan = structural_plan(root)
    with tempfile.TemporaryDirectory(prefix="leopard2-plan-expected-") as temporary:
        expected = Path(temporary) / "plan"
        generate_plan(expected)
        compare_trees(root, expected, "plan")
    return plan


def planned_gate_cells(root: Path, plan: dict[str, Any]) -> dict[tuple[int, int], dict[str, Any]]:
    result = {}
    for item in plan["artifacts"]:
        document = validate_main_document(load_json(root / item["path"]), item["path"])
        for cell in document["cells"]:
            key = (cell["K"], cell["shard_bytes"])
            require(key not in result, "planned gate cell is duplicated")
            result[key] = cell
    return result


def _normalize_bound_paths(value: object, replacements: tuple[tuple[str, str], ...]) -> object:
    """Remove volatile filesystem metadata while retaining exact byte identity."""
    require(replacements and all(isinstance(original, str) and original and
                                 isinstance(marker, str) and marker
                                 for original, marker in replacements) and
            all(Path(original).is_absolute() and original != "/" and
                os.path.normpath(original) == original
                for original, _ in replacements) and
            len({original for original, _ in replacements}) == len(replacements) and
            len({marker for _, marker in replacements}) == len(replacements),
            "scope path replacements are invalid or ambiguous")
    # A baseline build may live below the candidate source tree.  Replacing the
    # longest roots first prevents a broad source marker from swallowing the
    # more specific, role-distinct build identity.
    ordered = tuple(sorted(replacements, key=lambda item: len(item[0]), reverse=True))

    token_prefix_delimiters = frozenset(" \t\r\n\"'`=(:,;[{")
    token_suffix_delimiters = frozenset(" \t\r\n\"'`,:;)]}")
    cmake_external_object_suffix = re.compile(
        r"CMakeFiles/[^/\s\"'`]+\.dir$")

    def replace_root_tokens(text: str, original: str, marker: str) -> str:
        """Replace an absolute root only at a standalone path-token boundary."""
        pieces: list[str] = []
        cursor = 0
        while True:
            index = text.find(original, cursor)
            if index < 0:
                pieces.append(text[cursor:])
                return "".join(pieces)
            end = index + len(original)
            prefix = text[:index]
            cmake_match = cmake_external_object_suffix.search(prefix)
            cmake_external_source = cmake_match is not None and (
                cmake_match.start() == 0 or
                prefix[cmake_match.start() - 1] == "/" or
                prefix[cmake_match.start() - 1] in token_prefix_delimiters)
            fused_include = index >= 2 and text[index - 2:index] == "-I" and (
                index == 2 or text[index - 3] in token_prefix_delimiters)
            before_ok = index == 0 or \
                text[index - 1] in token_prefix_delimiters or \
                cmake_external_source or fused_include
            after_ok = end == len(text) or text[end] == "/" or \
                text[end] in token_suffix_delimiters
            pieces.append(text[cursor:index])
            if before_ok and after_ok:
                # CMake encodes an absolute external source below
                # CMakeFiles/<target>.dir by stripping no leading slash:
                # target.dir/home/user/source.cpp.o.  Preserve a separator
                # before the abstract root marker in that one typed context.
                pieces.append(("/" if cmake_external_source else "") + marker)
            else:
                pieces.append(original)
            cursor = end

    def visit(item: object) -> object:
        if isinstance(item, dict):
            result = {
                key: visit(child) for key, child in sorted(item.items())
                if key not in {"device", "inode", "mtime_ns"}
            }
            # Retained recipe text legitimately contains source/build roots.
            # Replacing those roots changes the normalized bytes, so recompute
            # their normalized content identity instead of leaving a stale
            # size/hash pair or retaining machine-specific absolute paths.
            if set(result) == {"encoding", "size", "sha256", "text"} and \
                    result.get("encoding") == "utf-8" and \
                    isinstance(result.get("text"), str):
                encoded = result["text"].encode("utf-8")
                result["size"] = len(encoded)
                result["sha256"] = hashlib.sha256(encoded).hexdigest()
            for prefix in ("archive", "executable"):
                content = result.get(f"{prefix}_link_recipe_content")
                identity = result.get(f"{prefix}_link_recipe")
                if isinstance(content, dict) and isinstance(identity, dict) and \
                        type(content.get("size")) is int and \
                        isinstance(content.get("sha256"), str):
                    identity["size"] = content["size"]
                    identity["sha256"] = content["sha256"]
            return result
        if isinstance(item, list):
            return [visit(child) for child in item]
        if isinstance(item, str):
            result = item
            for original, marker in ordered:
                result = replace_root_tokens(result, original, marker)
            return result
        return item

    return visit(value)


def _parse_scope_cpu_list(value: object, label: str) -> set[int]:
    require(isinstance(value, str) and value and
            len(value.encode("utf-8")) <= 65_536,
            f"{label} is not a bounded CPU-list string")
    result: set[int] = set()
    for component in value.split(","):
        component = component.strip()
        require(re.fullmatch(r"[0-9]+(?:-[0-9]+)?", component) is not None,
                f"{label} is malformed")
        bounds = [int(item) for item in component.split("-", 1)]
        first, last = bounds[0], bounds[-1]
        require(0 <= first <= last <= MAX_CPU_ID,
                f"{label} range is reversed or outside the CPU-ID bound")
        result.update(range(first, last + 1))
        require(len(result) <= MAX_CPU_LIST_ENTRIES,
                f"{label} expands beyond the CPU-count bound")
    require(result, f"{label} is empty")
    return result


def _validate_scope_artifact(
    value: object, label: str, expected_kind: str | None = None,
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "path", "kind", "size", "mode", "sha256"},
            f"{label} normalized artifact shape differs")
    require(isinstance(value.get("path"), str) and value["path"] and
            isinstance(value.get("kind"), str) and value["kind"] and
            (expected_kind is None or value["kind"] == expected_kind) and
            type(value.get("size")) is int and value["size"] >= 0 and
            type(value.get("mode")) is int and 0 <= value["mode"] <= 0o7777 and
            isinstance(value.get("sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is not None,
            f"{label} normalized artifact is invalid")
    return value


def _validate_scope_text(value: object, label: str) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "encoding", "size", "sha256", "text"} and
            value.get("encoding") == "utf-8" and
            isinstance(value.get("text"), str),
            f"{label} normalized retained text shape differs")
    encoded = value["text"].encode("utf-8")
    require(type(value.get("size")) is int and
            0 < value["size"] <= MAX_RETAINED_TEXT_BYTES and
            value["size"] == len(encoded) and "\x00" not in value["text"] and
            isinstance(value.get("sha256"), str) and
            value["sha256"] == hashlib.sha256(encoded).hexdigest(),
            f"{label} normalized retained text identity differs")
    return value


def _validate_scope_numeric_directory_inventory(
    value: object, prefix: str, label: str,
) -> list[str]:
    retained = _validate_scope_text(value, label)
    text = retained["text"]
    names = text.splitlines()
    require(names and text == "".join(name + "\n" for name in names) and
            all(name.startswith(prefix) and
                name.removeprefix(prefix).isdigit() for name in names) and
            names == sorted(set(names),
                            key=lambda name: int(name.removeprefix(prefix))),
            f"{label} is not a canonical numeric directory inventory")
    return names


def _validate_scope_commit_object(
    value: object, expected_commit: str, expected_tree: str, label: str,
) -> dict[str, Any]:
    require(isinstance(expected_commit, str) and
            re.fullmatch(r"[0-9a-f]{40}", expected_commit) is not None and
            isinstance(expected_tree, str) and
            re.fullmatch(r"[0-9a-f]{40}", expected_tree) is not None,
            f"{label} normalized expected Git identity is invalid")
    require(isinstance(value, dict) and set(value) == {
                "encoding", "size", "sha256", "object_id", "base64"} and
            value.get("encoding") == "base64" and
            isinstance(value.get("base64"), str),
            f"{label} normalized Git commit-object shape differs")
    try:
        raw = base64.b64decode(value["base64"], validate=True)
    except (ValueError, binascii.Error) as error:
        raise PlanError(
            f"{label} normalized Git commit object is not canonical base64") from error
    canonical = base64.b64encode(raw).decode("ascii")
    object_id = hashlib.sha1(
        f"commit {len(raw)}\0".encode("ascii") + raw).hexdigest()
    require(value["base64"] == canonical and
            type(value.get("size")) is int and value["size"] == len(raw) and
            0 < len(raw) <= MAX_GIT_COMMIT_BYTES and
            value.get("sha256") == hashlib.sha256(raw).hexdigest() and
            value.get("object_id") == expected_commit and
            object_id == expected_commit,
            f"{label} normalized Git commit-object byte identity differs")
    require(b"\n\n" in raw,
            f"{label} normalized Git commit has no header/message boundary")
    header_lines = raw.split(b"\n\n", 1)[0].splitlines()
    expected_tree_line = b"tree " + expected_tree.encode("ascii")
    trees = [line for line in header_lines if line.startswith(b"tree ")]
    require(header_lines and header_lines[0] == expected_tree_line and
            trees == [expected_tree_line],
            f"{label} normalized Git commit names another tree")
    return value


def _validate_scope_source(value: object, role: str) -> dict[str, Any]:
    require(role in {"baseline", "candidate"} and isinstance(value, dict) and
            set(value) == {
                "path", "head", "tree", "detached",
                "tracked_tree_listing_sha256", "tracked_status",
                "commit_object"},
            f"{role} normalized source identity shape differs")
    expected_path = "$BASELINE_SOURCE" if role == "baseline" else "$CANDIDATE_SOURCE"
    expected_head = EXACT_MAIN_COMMIT if role == "baseline" else None
    require(value.get("path") == expected_path and
            isinstance(value.get("head"), str) and
            re.fullmatch(r"[0-9a-f]{40}", value["head"]) is not None and
            (expected_head is None or value["head"] == expected_head) and
            isinstance(value.get("tree"), str) and
            re.fullmatch(r"[0-9a-f]{40}", value["tree"]) is not None and
            type(value.get("detached")) is bool and
            (role != "baseline" or value["detached"] is True) and
            isinstance(value.get("tracked_tree_listing_sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}",
                         value["tracked_tree_listing_sha256"]) is not None and
            value.get("tracked_status") == "clean",
            f"{role} normalized source identity differs")
    _validate_scope_commit_object(
        value.get("commit_object"), value["head"], value["tree"], role)
    return value


def _normalized_compile_output(role: str, source: str) -> str:
    baseline = role == "baseline"
    if baseline:
        if source == ("$CANDIDATE_SOURCE/experiments/leopard2/main_compare/"
                      "legacy_main_benchmark.cpp"):
            return ("CMakeFiles/leopard_main_benchmark.dir/"
                    "legacy_main_benchmark.cpp.o")
        require(source.startswith("$BASELINE_SOURCE/"),
                "baseline normalized compiler source escapes its root")
        relative = source[len("$BASELINE_SOURCE/"):]
        return ("CMakeFiles/leopard_main_exact.dir/$BASELINE_SOURCE/" +
                relative + ".o")
    require(role == "candidate" and
            source.startswith("$CANDIDATE_SOURCE/"),
            "candidate normalized compiler source escapes its root")
    relative = source[len("$CANDIDATE_SOURCE/"):]
    if relative == "bench/leopard2/benchmark.cpp":
        return "CMakeFiles/bench_leopard2.dir/bench/leopard2/benchmark.cpp.o"
    backend_targets = {
        "Leopard2BackendSSSE3.cpp": "leopard2_backend_ssse3.dir",
        "Leopard2BackendAVX2.cpp": "leopard2_backend_avx2.dir",
        "Leopard2BackendAVX512.cpp": "leopard2_backend_avx512.dir",
    }
    return f"CMakeFiles/{backend_targets.get(relative, 'leopard.dir')}/{relative}.o"


def _normalized_compile_argv(
    role: str, source: str, compiler_invocation: str,
) -> list[str]:
    output = _normalized_compile_output(role, source)
    if role == "baseline":
        adapter = ("$CANDIDATE_SOURCE/experiments/leopard2/main_compare/"
                   "legacy_main_benchmark.cpp")
        definitions = ([] if source != adapter else [
            f'-DLEOPARD_MAIN_SOURCE_COMMIT="{EXACT_MAIN_COMMIT}"',
        ])
        return [
            compiler_invocation, *definitions, "-I$BASELINE_SOURCE",
            "-g", "-O0", "-O3", "-std=gnu++11", "-march=native",
            "-Wall", "-Wextra", "-fopenmp",
            "-o", output, "-c", source,
        ]

    relative = source[len("$CANDIDATE_SOURCE/"):]
    isolated_flags = {
        "Leopard2BackendSSSE3.cpp": ["-mssse3", "-mno-avx"],
        "Leopard2BackendAVX2.cpp": [
            "-mavx2", "-mno-avx512f", "-falign-functions=64"],
        "Leopard2BackendAVX512.cpp": [
            "-mavx2", "-mavx512f", "-mavx512bw", "-mavx512vl",
            "-mprefer-vector-width=256", "-falign-functions=64"],
    }
    if relative == "Leopard2BackendAVX512.cpp":
        definitions = ["-DLEO2_HAVE_AVX2_BACKEND=1"]
    elif relative in isolated_flags or relative == "bench/leopard2/benchmark.cpp":
        definitions = []
    else:
        definitions = [
            "-DLEO2_DISABLE_AVX2_CODEGEN=1",
            "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
            "-DLEO2_HAVE_AVX2_BACKEND=1",
            "-DLEO2_HAVE_AVX512_BACKEND=1",
            "-DLEO2_HAVE_SSSE3_BACKEND=1",
        ]
    return [
        compiler_invocation, *definitions, "-I$CANDIDATE_SOURCE",
        "-Wall", "-Wextra", "-fopenmp", "-O3", "-DNDEBUG", "-O3",
        "-std=gnu++11", *isolated_flags.get(relative, []),
        *([] if relative in isolated_flags else ["-fopenmp"]),
        "-o", output, "-c", source,
    ]


def _validate_scope_build(build: object, role: str) -> dict[str, Any]:
    require(role in {"baseline", "candidate"} and isinstance(build, dict) and
            set(build) == {
                "build_dir", "cmake_cache", "compile_commands",
                "executable_link_recipe", "archive_link_recipe", "compiler",
                "compiler_version_stdout", "archiver", "ranlib",
                "validated_archive_members", "validated_executable",
                "validated_archive", "validated_cache",
                "validated_compile_commands", "archive_link_recipe_content",
                "executable_link_recipe_content",
                "archive_link_tool_invocations", "compiler_invocation",
                "validated_external_link_inputs"},
            f"{role} normalized build identity shape differs")
    expected_root = "$BASELINE_BUILD" if role == "baseline" else "$CANDIDATE_BUILD"
    require(build.get("build_dir") == expected_root,
            f"{role} normalized build directory differs")
    baseline = role == "baseline"
    archive_name = "libleopard_main_exact.a" if baseline else "libleopard.a"
    executable_name = "leopard_main_benchmark" if baseline else "bench_leopard2"
    archive_target = "leopard_main_exact.dir" if baseline else "leopard.dir"
    library_sources = (
        BASELINE_LIBRARY_SOURCES if baseline else CANDIDATE_LIBRARY_SOURCES)
    expected_entry_count = (
        BASELINE_EXPECTED_COMPILE_COMMAND_COUNT if baseline else
        CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT)
    expected_compile_profile = (
        BASELINE_COMPILE_PROFILE if baseline else CANDIDATE_COMPILE_PROFILE)
    metadata_paths = {
        "cmake_cache": f"{expected_root}/CMakeCache.txt",
        "compile_commands": f"{expected_root}/compile_commands.json",
        "executable_link_recipe": (
            f"{expected_root}/CMakeFiles/{executable_name}.dir/link.txt"),
        "archive_link_recipe": (
            f"{expected_root}/CMakeFiles/{archive_target}/link.txt"),
    }
    for name in ("cmake_cache", "compile_commands", "executable_link_recipe",
                 "archive_link_recipe"):
        record = _validate_scope_artifact(
            build.get(name), f"{role} {name}", "build_metadata")
        require(record["path"] == metadata_paths[name],
                f"{role} normalized {name} path differs")
    compiler = _validate_scope_artifact(
        build.get("compiler"), f"{role} compiler", "compiler")
    archiver = _validate_scope_artifact(
        build.get("archiver"), f"{role} archiver", "archiver")
    ranlib = _validate_scope_artifact(
        build.get("ranlib"), f"{role} ranlib", "ranlib")
    archive = _validate_scope_artifact(
        build.get("validated_archive"), f"{role} archive", "archive")
    executable = _validate_scope_artifact(
        build.get("validated_executable"), f"{role} executable", "executable")
    require(archive["path"] == f"{expected_root}/{archive_name}" and
            executable["path"] == f"{expected_root}/{executable_name}",
            f"{role} normalized output paths differ")
    version = build.get("compiler_version_stdout")
    require(isinstance(version, dict) and set(version) == {"sha256", "text"} and
            isinstance(version.get("text"), str) and version["text"] and
            isinstance(version.get("sha256"), str) and
            version["sha256"] == hashlib.sha256(
                version["text"].encode("utf-8")).hexdigest(),
            f"{role} compiler version identity differs")
    cache = build.get("validated_cache")
    required_cache = (REQUIRED_BASELINE_CACHE if role == "baseline"
                      else REQUIRED_CANDIDATE_CACHE)
    expected_cache_keys = {
        "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS_RELEASE",
        *required_cache,
    }
    require(isinstance(cache, dict) and set(cache) == expected_cache_keys and
            isinstance(cache.get("CMAKE_CXX_COMPILER"), str) and
            cache["CMAKE_CXX_COMPILER"] and
            isinstance(cache.get("CMAKE_CXX_FLAGS_RELEASE"), str) and
            all(cache.get(key) == expected
                for key, expected in required_cache.items()),
            f"{role} validated CMake cache differs")
    try:
        release_flags = shlex.split(
            cache["CMAKE_CXX_FLAGS_RELEASE"], posix=True)
    except ValueError as error:
        raise PlanError(
            f"cannot parse {role} normalized CMake Release flags: {error}") \
            from error
    try:
        common.validate_effective_flags(
            release_flags, f"{role} normalized CMake Release flags",
            "release")
    except common.EvidenceError as error:
        raise PlanError(str(error)) from error
    compiler_invocation = build.get("compiler_invocation")
    require(isinstance(compiler_invocation, dict) and
            set(compiler_invocation) == {"invocation", "resolved_path"} and
            compiler_invocation.get("invocation") ==
                cache["CMAKE_CXX_COMPILER"] and
            compiler_invocation.get("resolved_path") == compiler["path"],
            f"{role} normalized compiler invocation differs")
    semantics = build.get("validated_compile_commands")
    require(isinstance(semantics, dict) and set(semantics) == {
                "entry_count", "required_sources", "validated_optimization",
                "validated_openmp", "required_source_object_pairs", "isa_policy",
                "schema", "implementation", "profile", "required_entries"} and
            type(semantics.get("entry_count")) is int and
            semantics["entry_count"] == expected_entry_count and
            semantics.get("schema") == COMPILE_COMMANDS_SCHEMA and
            semantics.get("implementation") == role and
            semantics.get("profile") == expected_compile_profile and
            semantics.get("validated_optimization") == "-O3" and
            semantics.get("validated_openmp") is True and
            isinstance(semantics.get("isa_policy"), str) and
            semantics["isa_policy"] and
            isinstance(semantics.get("required_sources"), list) and
            semantics["required_sources"] and
            len(semantics["required_sources"]) == len(set(
                semantics["required_sources"])) and
            isinstance(semantics.get("required_source_object_pairs"), list) and
            semantics["required_source_object_pairs"] and
            isinstance(semantics.get("required_entries"), list) and
            len(semantics["required_entries"]) ==
                len(semantics["required_source_object_pairs"]),
            f"{role} normalized compile-command identity differs")
    entries_by_source: dict[str, dict[str, Any]] = {}
    for index, entry in enumerate(semantics["required_entries"]):
        require(isinstance(entry, dict) and set(entry) == {
                    "directory", "file", "output", "arguments"} and
                entry.get("directory") == expected_root and
                isinstance(entry.get("file"), str) and entry["file"] and
                isinstance(entry.get("output"), str) and entry["output"] and
                isinstance(entry.get("arguments"), list) and
                all(isinstance(token, str) and token and "@" not in token
                    for token in entry["arguments"]) and
                entry["file"] not in entries_by_source,
                f"{role} normalized compiler entry {index} differs")
        entries_by_source[entry["file"]] = entry
    sources: list[str] = []
    objects: list[str] = []
    for index, pair in enumerate(semantics["required_source_object_pairs"]):
        require(isinstance(pair, dict) and set(pair) == {"source", "object"},
                f"{role} normalized compile pair {index} shape differs")
        source = _validate_scope_artifact(
            pair["source"], f"{role} source {index}", "source_file")
        obj = _validate_scope_artifact(
            pair["object"], f"{role} object {index}", "object_file")
        entry = entries_by_source.get(source["path"])
        expected_output = _normalized_compile_output(role, source["path"])
        require(entry == {
                    "directory": expected_root,
                    "file": source["path"],
                    "output": expected_output,
                    "arguments": _normalized_compile_argv(
                        role, source["path"],
                        compiler_invocation["invocation"]),
                } and
                obj["path"] == f"{expected_root}/{expected_output}",
                f"{role} normalized compiler argv/output {index} differs")
        sources.append(source["path"])
        objects.append(obj["path"])
    require(sources == semantics["required_sources"] and
            len(sources) == len(set(sources)) and
            len(objects) == len(set(objects)) and
            set(entries_by_source) == set(sources) and
            semantics["entry_count"] >= len(sources),
            f"{role} normalized compile source/object closure differs")
    benchmark_suffix = (
        "/experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        if baseline else "/bench/leopard2/benchmark.cpp")
    implementation_source = "$BASELINE_SOURCE" if baseline else "$CANDIDATE_SOURCE"
    expected_sources = sorted([
        *(f"{implementation_source}/{name}" for name in library_sources),
        "$CANDIDATE_SOURCE" + benchmark_suffix,
    ])
    require(len(sources) == len(expected_sources) and
            set(sources) == set(expected_sources),
            f"{role} normalized translation-unit set differs from the producer")
    benchmark_pairs = [
        pair for pair in semantics["required_source_object_pairs"]
        if pair["source"]["path"].endswith(benchmark_suffix)]
    require(len(benchmark_pairs) == 1 and
            all(path.startswith(expected_root + "/") for path in objects),
            f"{role} normalized benchmark/build object closure differs")
    archive_pairs = [pair for pair in semantics["required_source_object_pairs"]
                     if pair not in benchmark_pairs]
    require(archive_pairs and
            all(pair["source"]["path"].startswith("$CANDIDATE_SOURCE/")
                for pair in benchmark_pairs) and
            all(pair["source"]["path"].startswith(implementation_source + "/")
                for pair in archive_pairs),
            f"{role} normalized compile sources escape their source roots")
    members = build.get("validated_archive_members")
    expected_members = [Path(name).name + ".o" for name in library_sources]
    require(isinstance(members, list) and members == expected_members and
            len(members) == len(set(members)) and
            all(isinstance(member, str) and member and "/" not in member
                for member in members),
            f"{role} normalized archive members differ")
    objects_by_member = {
        pair["object"]["path"].rsplit("/", 1)[-1]: pair["object"]["path"]
        for pair in archive_pairs
    }
    require(len(objects_by_member) == len(archive_pairs) and
            set(objects_by_member) == set(members),
            f"{role} normalized archive member/object closure differs")
    archive_text = _validate_scope_text(
        build.get("archive_link_recipe_content"),
        f"{role} archive link recipe")
    executable_text = _validate_scope_text(
        build.get("executable_link_recipe_content"),
        f"{role} executable link recipe")
    require(archive_text["size"] == build["archive_link_recipe"]["size"] and
            archive_text["sha256"] == build["archive_link_recipe"]["sha256"] and
            executable_text["size"] == build["executable_link_recipe"]["size"] and
            executable_text["sha256"] ==
                build["executable_link_recipe"]["sha256"],
            f"{role} normalized retained link bytes differ from file identity")
    tools = build.get("archive_link_tool_invocations")
    require(isinstance(tools, dict) and set(tools) == {"archiver", "ranlib"},
            f"{role} archive tool invocation shape differs")
    for name, artifact in (("archiver", archiver), ("ranlib", ranlib)):
        invocation = tools.get(name)
        require(isinstance(invocation, dict) and set(invocation) == {
                    "invocation", "resolved_path"} and
                isinstance(invocation.get("invocation"), str) and
                invocation["invocation"] and
                invocation.get("resolved_path") == artifact["path"],
                f"{role} normalized {name} invocation differs")
    archive_commands = [
        shlex.split(line, posix=True)
        for line in archive_text["text"].splitlines() if line.strip()
    ]
    expected_archive_objects = [
        objects_by_member[member][len(expected_root) + 1:]
        for member in members
    ]
    require(len(archive_commands) == 2 and
            len(archive_commands[0]) >= 4 and
            archive_commands[0][0] == tools["archiver"]["invocation"] and
            archive_commands[0][1] in {"qc", "rc", "rcs"} and
            archive_commands[0][2] == archive_name and
            archive_commands[0][3:] == expected_archive_objects and
            archive_commands[1] == [tools["ranlib"]["invocation"], archive_name],
            f"{role} normalized archive recipe semantics differ")
    executable_tokens = shlex.split(executable_text["text"], posix=True)
    expected_benchmark_object = benchmark_pairs[0]["object"]["path"][
        len(expected_root) + 1:]
    external_link_inputs = build.get("validated_external_link_inputs")
    require(isinstance(external_link_inputs, list),
            f"{role} normalized external-link identity is absent")
    for index, record in enumerate(external_link_inputs):
        require(isinstance(record, dict) and
                isinstance(record.get("artifact"), dict),
                f"{role} normalized external-link identity {index} differs")
        _validate_scope_artifact(
            record["artifact"],
            f"{role} normalized external-link artifact {index}")
    try:
        common.validate_executable_link_semantics(
            executable_tokens,
            compiler_invocation=compiler_invocation["invocation"],
            archive_name=archive_name, executable_name=executable_name,
            benchmark_object=expected_benchmark_object,
            external_link_inputs=external_link_inputs,
            label=f"{role} normalized executable link recipe")
    except common.EvidenceError as error:
        raise PlanError(
            f"{role} normalized executable recipe semantics differ: {error}") \
            from error
    require(compiler["path"], f"{role} normalized compiler path is empty")
    return build


def _parse_scope_ldd_output(value: str, label: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for raw_line in value.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if "=>" in line:
            soname, target_and_address = (
                part.strip() for part in line.split("=>", 1))
            require(not target_and_address.startswith("not found"),
                    f"{label} retains a missing dependency")
            target = target_and_address.split(" (", 1)[0]
            require(target.startswith("/") and os.path.normpath(target) == target,
                    f"{label} retains an unresolved dependency")
            entries.append({
                "soname": soname, "loader_path": target,
                "file_kind": "shared_library",
            })
            continue
        token = line.split(" (", 1)[0]
        if token.startswith("/"):
            require(os.path.normpath(token) == token,
                    f"{label} retains a non-canonical dynamic-loader path")
            entries.append({
                "soname": Path(token).name, "loader_path": token,
                "file_kind": "dynamic_loader",
            })
        elif token == "linux-vdso.so.1":
            entries.append({"soname": token, "virtual": True})
        else:
            raise PlanError(f"{label} contains unrecognized ldd output: {line}")
    require(entries and
            len({entry["soname"] for entry in entries}) == len(entries) and
            sum(entry.get("file_kind") == "dynamic_loader"
                for entry in entries) == 1 and
            any(entry.get("file_kind") == "shared_library"
                for entry in entries) and
            sum(entry.get("virtual") is True for entry in entries) <= 1,
            f"{label} is not an internally consistent, unique loader record")
    return sorted(entries, key=lambda item: item["soname"])


def _validate_runtime_closure(value: object, label: str) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "executable", "dependencies", "raw_ldd_output"} and
            isinstance(value.get("executable"), str) and value["executable"] and
            isinstance(value.get("dependencies"), list) and value["dependencies"],
            f"{label} normalized runtime closure shape differs")
    raw = _validate_scope_text(
        value.get("raw_ldd_output"), f"{label} normalized raw ldd output")
    parsed = _parse_scope_ldd_output(
        raw["text"], f"{label} normalized raw ldd output")
    sonames: list[str] = []
    loader_paths: list[str] = []
    file_paths: list[str] = []
    for dependency in value["dependencies"]:
        require(isinstance(dependency, dict) and
                isinstance(dependency.get("soname"), str) and
                dependency["soname"],
                f"{label} normalized runtime dependency is invalid")
        sonames.append(dependency["soname"])
        if set(dependency) == {"soname", "virtual"}:
            require(dependency["soname"] == "linux-vdso.so.1" and
                    dependency.get("virtual") is True,
                    f"{label} normalized virtual dependency differs")
            continue
        require(set(dependency) == {"soname", "loader_path", "file"} and
                isinstance(dependency.get("loader_path"), str) and
                dependency["loader_path"],
                f"{label} normalized runtime dependency variant differs")
        file_record = _validate_scope_artifact(
            dependency["file"], f"{label} runtime {dependency['soname']}")
        require(file_record["kind"] in {"shared_library", "dynamic_loader"} and
                file_record["path"] == dependency["loader_path"],
                f"{label} normalized runtime dependency kind/loader path differs")
        loader_paths.append(dependency["loader_path"])
        file_paths.append(file_record["path"])
    require(sonames == sorted(sonames) and len(sonames) == len(set(sonames)) and
            len(loader_paths) == len(set(loader_paths)) and
            len(file_paths) == len(set(file_paths)),
            f"{label} normalized runtime closure is not sorted and unique")
    normalized_dependencies = []
    for dependency in value["dependencies"]:
        if dependency.get("virtual") is True:
            normalized_dependencies.append(dict(dependency))
        else:
            normalized_dependencies.append({
                "soname": dependency["soname"],
                "loader_path": dependency["loader_path"],
                "file_kind": dependency["file"]["kind"],
            })
    require(normalized_dependencies == parsed,
            f"{label} normalized dependency summary differs from raw ldd bytes")
    return value


def _validate_scope_cpu_policy(
    value: object, label: str, expected_sibling: int,
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "cpu", "online", "cpuinfo", "topology", "frequency_policy",
                "cache_hierarchy", "cache_index_inventory",
                "cache_directory_inventory_text", "numa_nodes",
                "numa_node_inventory", "core_class"} and
            type(value.get("cpu")) is int and
            0 <= value["cpu"] <= MAX_CPU_ID and
            (value.get("online") is None or isinstance(value["online"], str)),
            f"{label} normalized CPU policy shape differs")
    cpu = value["cpu"]
    cpuinfo = value.get("cpuinfo")
    allowed_cpuinfo = {
        "processor", "vendor_id", "cpu family", "model", "model name",
        "stepping", "microcode", "flags", "Features", "CPU implementer",
        "CPU architecture", "CPU variant", "CPU part", "CPU revision",
    }
    require(isinstance(cpuinfo, dict) and cpuinfo and
            set(cpuinfo).issubset(allowed_cpuinfo) and
            cpuinfo.get("processor") == str(cpu) and
            any(name in cpuinfo for name in ("model name", "CPU part")) and
            all(isinstance(item, str) for item in cpuinfo.values()),
            f"{label} normalized CPU model identity differs")
    topology = value.get("topology")
    require(isinstance(topology, dict) and set(topology) == {
                "core_id", "physical_package_id", "die_id", "cluster_id",
                "thread_siblings_list", "core_siblings_list"} and
            all(item is None or isinstance(item, str)
                for item in topology.values()) and
            _parse_scope_cpu_list(topology.get("thread_siblings_list"),
                                  f"{label} thread siblings") == {
                cpu, expected_sibling} and
            {cpu, expected_sibling}.issubset(_parse_scope_cpu_list(
                topology.get("core_siblings_list"),
                f"{label} core siblings")),
            f"{label} normalized CPU topology differs")
    frequency = value.get("frequency_policy")
    require(isinstance(frequency, dict) and set(frequency) == {
                "scaling_driver", "scaling_governor",
                "energy_performance_preference", "scaling_min_freq",
                "scaling_max_freq", "cpuinfo_min_freq", "cpuinfo_max_freq"} and
            all(item is None or isinstance(item, str)
                for item in frequency.values()),
            f"{label} normalized frequency identity differs")
    caches = value.get("cache_hierarchy")
    require(isinstance(caches, list) and caches,
            f"{label} normalized cache hierarchy is absent")
    indices = []
    for cache in caches:
        require(isinstance(cache, dict) and set(cache) == {
                    "index", "level", "type", "size", "coherency_line_size",
                    "number_of_sets", "ways_of_associativity",
                    "physical_line_partition", "shared_cpu_list", "shared_cpu_map",
                    "allocation_policy", "write_policy"} and
                type(cache.get("index")) is int and cache["index"] >= 0 and
                all(cache.get(name) is None or isinstance(cache[name], str)
                    for name in set(cache) - {"index"}) and
                all(isinstance(cache.get(name), str) and cache[name]
                    for name in ("level", "type", "size", "coherency_line_size",
                                 "shared_cpu_list", "shared_cpu_map")) and
                cpu in _parse_scope_cpu_list(
                    cache["shared_cpu_list"], f"{label} cache shared CPUs"),
                f"{label} normalized cache identity differs")
        indices.append(cache["index"])
    require(indices == sorted(indices) and len(indices) == len(set(indices)),
            f"{label} normalized cache indices differ")
    cache_inventory = value.get("cache_index_inventory")
    raw_cache_inventory = _validate_scope_numeric_directory_inventory(
        value.get("cache_directory_inventory_text"), "index",
        f"{label} raw cache-directory inventory")
    require(isinstance(cache_inventory, list) and
            cache_inventory == [f"index{index}" for index in indices] and
            cache_inventory == raw_cache_inventory,
            f"{label} normalized cache inventory differs")
    numa = value.get("numa_nodes")
    node_inventory = value.get("numa_node_inventory")
    core_class = value.get("core_class")
    require(isinstance(numa, list) and numa and numa == sorted(set(numa)) and
            all(type(node) is int and 0 <= node <= MAX_CPU_ID for node in numa) and
            isinstance(node_inventory, list) and
            node_inventory == [f"node{node}" for node in numa] and
            isinstance(core_class, dict) and set(core_class) == {
                "core_type", "cpu_capacity"} and
            all(item is None or isinstance(item, str)
                for item in core_class.values()),
            f"{label} normalized NUMA/core-class identity differs")
    return value


def _validate_scope_host(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "system", "allowed_cpu_set_at_launch", "online_cpu_set",
                "online_cpu_list_text", "online_node_list_text",
                "benchmark_cpu", "reserved_sibling", "turbo_and_pstate"},
            "gate normalized host identity shape differs")
    system = value.get("system")
    require(isinstance(system, dict) and set(system) == {
                "system", "node", "release", "version", "machine", "python",
                "page_size"} and
            all(isinstance(system.get(name), str) and system[name]
                for name in ("system", "node", "release", "version", "machine",
                             "python")) and
            type(system.get("page_size")) is int and system["page_size"] > 0,
            "gate normalized host system identity differs")
    allowed = value.get("allowed_cpu_set_at_launch")
    online = value.get("online_cpu_set")
    require(isinstance(allowed, list) and allowed == sorted(set(allowed)) and
            3 <= len(allowed) <= MAX_CPU_LIST_ENTRIES and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in allowed) and
            isinstance(online, list) and online == sorted(set(online)) and
            1 <= len(online) <= MAX_CPU_LIST_ENTRIES and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in online),
            "gate normalized host CPU sets differ")
    raw_online = _validate_scope_text(
        value.get("online_cpu_list_text"), "gate normalized online CPU list")
    raw_nodes = _validate_scope_text(
        value.get("online_node_list_text"), "gate normalized online NUMA list")
    online_nodes = _parse_scope_cpu_list(
        raw_nodes["text"], "gate normalized online NUMA list")
    require(_parse_scope_cpu_list(
                raw_online["text"], "gate normalized online CPU list") ==
                set(online),
            "gate normalized online CPU summary differs from sysfs bytes")
    benchmark = value.get("benchmark_cpu")
    sibling = value.get("reserved_sibling")
    require(isinstance(benchmark, dict) and isinstance(sibling, dict) and
            type(benchmark.get("cpu")) is int and
            type(sibling.get("cpu")) is int and
            benchmark["cpu"] != sibling["cpu"],
            "gate normalized host reserved pair differs")
    _validate_scope_cpu_policy(
        benchmark, "benchmark CPU", sibling["cpu"])
    _validate_scope_cpu_policy(
        sibling, "reserved sibling", benchmark["cpu"])
    pair = {benchmark["cpu"], sibling["cpu"]}
    require(pair.issubset(set(allowed)) and pair.issubset(set(online)) and
            benchmark["cache_hierarchy"] == sibling["cache_hierarchy"] and
            benchmark["cache_index_inventory"] ==
                sibling["cache_index_inventory"] and
            benchmark["numa_nodes"] == sibling["numa_nodes"] and
            benchmark["numa_node_inventory"] ==
                sibling["numa_node_inventory"] and
            set(benchmark["numa_nodes"]).issubset(online_nodes) and
            benchmark["core_class"] == sibling["core_class"],
            "gate normalized pair is outside CPU sets or differs in cache/class/domain")
    turbo = value.get("turbo_and_pstate")
    require(isinstance(turbo, dict) and set(turbo) == {
                "/sys/devices/system/cpu/intel_pstate/no_turbo",
                "/sys/devices/system/cpu/amd_pstate/status",
                "/sys/devices/system/cpu/cpufreq/boost"} and
            all(item is None or isinstance(item, str) for item in turbo.values()),
            "gate normalized host turbo/pstate identity differs")
    return value


def selection_scope_from_verified_bundle(
    manifest: dict[str, Any], raw: dict[str, Any]
) -> dict[str, Any]:
    """Derive the one environment in which every gate cell is meaningful."""
    require(raw.get("schema") == EXACT_RAW_SCHEMA,
            "exact-main raw bundle is not complete-identity schema v5")
    identities = manifest.get("identities")
    host = manifest.get("host")
    require(isinstance(identities, dict) and isinstance(host, dict),
            "exact-main scope identity is incomplete")
    baseline_build = identities.get("baseline_build")
    candidate_build = identities.get("candidate_build")
    baseline_source = identities.get("baseline_source")
    candidate_source = identities.get("candidate_source")
    require(all(isinstance(item, dict) for item in (
        baseline_build, candidate_build, baseline_source, candidate_source)),
        "exact-main baseline/candidate build or source scope is absent")
    baseline_source_root = baseline_source.get("path")
    candidate_source_root = candidate_source.get("path")
    baseline_build_root = baseline_build.get("build_dir")
    candidate_build_root = candidate_build.get("build_dir")
    require(all(isinstance(item, str) and item for item in (
        baseline_source_root, candidate_source_root,
        baseline_build_root, candidate_build_root)),
        "exact-main role-specific source/build roots are absent")

    backends = {
        invocation.get("normalized", {}).get("backend")
        for invocation in raw.get("invocations", [])
        if isinstance(invocation, dict) and
           invocation.get("implementation") == "candidate" and
           isinstance(invocation.get("normalized"), dict)
    }
    require(len(backends) == 1, "exact-main candidate resolved multiple/no backends")
    resolved_backend = next(iter(backends))
    require(resolved_backend in CAMPAIGN_BACKENDS,
            "resolved AUTO backend lacks forced confirmation coverage; "
            "AVX-512 and NEON are outside this campaign")

    replacements = (
        (baseline_source_root, "$BASELINE_SOURCE"),
        (candidate_source_root, "$CANDIDATE_SOURCE"),
        (baseline_build_root, "$BASELINE_BUILD"),
        (candidate_build_root, "$CANDIDATE_BUILD"),
    )
    scope = {
        "schema": "leopard2-balanced-evidence-scope/v3",
        # Retain the exact pinned CPU pair.  This is intentionally stricter than
        # topology cardinality: the current evidence format lacks cache/NUMA and
        # heterogeneous-core descriptors needed to equate different pairs.
        "host": _normalize_bound_paths(host, replacements),
        "sources": {
            "baseline": _normalize_bound_paths(baseline_source, replacements),
            "candidate": _normalize_bound_paths(candidate_source, replacements),
        },
        "builds": {
            "baseline": _normalize_bound_paths(baseline_build, replacements),
            "candidate": _normalize_bound_paths(candidate_build, replacements),
        },
        "artifacts": {
            key: _normalize_bound_paths(identities.get(key), replacements)
            for key in (
                "baseline_archive", "baseline_executable",
                "candidate_archive", "candidate_executable")
        },
        "runtime_closures": {
            "baseline": _normalize_bound_paths(
                identities.get("baseline_runtime_closure"), replacements),
            "candidate": _normalize_bound_paths(
                identities.get("candidate_runtime_closure"), replacements),
        },
        "tools": {
            key: _normalize_bound_paths(identities.get(key), replacements)
            for key in ("runner", "taskset", "ldd")
        },
        "resolved_auto_backend": resolved_backend,
        "forced_confirmation_backends": list(
            BACKENDS[:BACKENDS.index(resolved_backend) + 1]),
        "excluded_backends": dict(EXCLUDED_CAMPAIGN_BACKENDS),
    }
    return scope


def validate_evidence_scope(scope: object) -> dict[str, Any]:
    require(isinstance(scope, dict) and set(scope) == {
        "schema", "host", "sources", "builds", "artifacts",
        "runtime_closures", "tools",
        "resolved_auto_backend", "forced_confirmation_backends",
        "excluded_backends",
    } and scope.get("schema") == "leopard2-balanced-evidence-scope/v3",
            "gate evidence scope shape differs")
    backend = scope.get("resolved_auto_backend")
    require(backend in CAMPAIGN_BACKENDS and
            scope.get("forced_confirmation_backends") ==
                list(BACKENDS[:BACKENDS.index(backend) + 1]) and
            scope.get("excluded_backends") == EXCLUDED_CAMPAIGN_BACKENDS,
            "gate evidence backend coverage declaration differs")
    host = scope.get("host")
    sources = scope.get("sources")
    builds = scope.get("builds")
    artifacts = scope.get("artifacts")
    closures = scope.get("runtime_closures")
    tools = scope.get("tools")
    _validate_scope_host(host)
    require(isinstance(sources, dict) and set(sources) == {"baseline", "candidate"} and
            isinstance(builds, dict) and set(builds) == {"baseline", "candidate"},
            "gate evidence role-specific source/build scope differs")
    _validate_scope_source(sources["baseline"], "baseline")
    _validate_scope_source(sources["candidate"], "candidate")
    _validate_scope_build(builds["baseline"], "baseline")
    _validate_scope_build(builds["candidate"], "candidate")
    require(builds["baseline"]["compiler"] == builds["candidate"]["compiler"] and
            builds["baseline"]["compiler_version_stdout"] ==
                builds["candidate"]["compiler_version_stdout"],
            "gate evidence baseline/candidate compiler identity differs")
    require(builds["baseline"]["validated_external_link_inputs"] ==
                builds["candidate"]["validated_external_link_inputs"],
            "gate evidence baseline/candidate external link inputs differ")
    require(isinstance(artifacts, dict) and set(artifacts) == {
        "baseline_archive", "baseline_executable",
        "candidate_archive", "candidate_executable",
    }, "gate evidence artifact closure shape differs")
    for key, artifact_value in artifacts.items():
        _validate_scope_artifact(
            artifact_value, key,
            "archive" if key.endswith("archive") else "executable")
    require(artifacts["baseline_archive"] ==
                builds["baseline"]["validated_archive"] and
            artifacts["baseline_executable"] ==
                builds["baseline"]["validated_executable"] and
            artifacts["candidate_archive"] ==
                builds["candidate"]["validated_archive"] and
            artifacts["candidate_executable"] ==
                builds["candidate"]["validated_executable"],
            "gate evidence top-level/build artifact closure differs")
    require(isinstance(closures, dict) and set(closures) == {
        "baseline", "candidate",
    }, "gate evidence runtime closure shape differs")
    _validate_runtime_closure(closures["baseline"], "baseline")
    _validate_runtime_closure(closures["candidate"], "candidate")
    require(closures["baseline"]["executable"] ==
                artifacts["baseline_executable"].get("path") and
            closures["candidate"]["executable"] ==
                artifacts["candidate_executable"].get("path"),
            "gate evidence runtime closure names a different executable")
    require(isinstance(tools, dict) and set(tools) == {"runner", "taskset", "ldd"},
            "gate evidence tool scope shape differs")
    for key, tool in tools.items():
        _validate_scope_artifact(
            tool, key, "file" if key == "runner" else "executable")
    return scope


def verify_exact_manifest(
    path: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    try:
        document, raw, _, snapshot = \
            load_exact_main_runner().verified_campaign_bundle(
                path, no_current_input_check=True)
    except Exception as error:
        raise PlanError(f"exact-main verifier rejected {path}: {error}") from error
    require(isinstance(document, dict) and
            document.get("schema") == EXACT_MANIFEST_SCHEMA and
            document.get("valid") is True,
            f"exact-main manifest is not valid schema v5: {path}")
    require(isinstance(snapshot, dict) and set(snapshot) == {"size", "sha256"} and
            type(snapshot.get("size")) is int and snapshot["size"] > 0 and
            isinstance(snapshot.get("sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", snapshot["sha256"]) is not None,
            f"exact-main manifest snapshot identity is invalid: {path}")
    return (document, validate_evidence_scope(
        selection_scope_from_verified_bundle(document, raw)), snapshot)


def manifest_cell(cell: object) -> dict[str, Any]:
    require(isinstance(cell, dict) and set(cell) == {
        "identifier", "k", "r", "shard_bytes", "losses", "seed",
    }, "exact-main campaign cell shape differs")
    require(isinstance(cell["identifier"], str),
            "exact-main campaign identifier is invalid")
    for key in ("k", "r", "shard_bytes", "losses", "seed"):
        require(type(cell[key]) is int and cell[key] > 0,
                f"exact-main campaign {key} is invalid")
    expected = main_cell(
        f"gate-k{cell['k']}-b{cell['shard_bytes']}",
        cell["k"], cell["r"], cell["shard_bytes"], cell["losses"],
                         "exact-main-generic-gate")
    require(cell == {
        "identifier": expected["identifier"], "k": expected["K"],
        "r": expected["R"], "shard_bytes": expected["shard_bytes"],
        "losses": expected["loss_count"], "seed": expected["seed"],
    }, "exact-main campaign cell differs from its canonical gate cell")
    return expected


def refinement_cells(evaluated: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: dict[tuple[int, int], dict[str, Any]] = {}
    grouped: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for cell in evaluated:
        side = 1 << (cell["K"] - 1).bit_length()
        grouped.setdefault((side, cell["shard_bytes"]), []).append(cell)
    for (_, shard_bytes), cells in grouped.items():
        cells.sort(key=lambda item: item["K"])
        for first, second in zip(cells, cells[1:]):
            if first["passes_gate"] == second["passes_gate"] or \
                    second["K"] == first["K"] + 1:
                continue
            for k in range(first["K"] + 1, second["K"]):
                candidate = main_cell(
                    f"gate-k{k}-b{shard_bytes}", k, k, shard_bytes, k,
                    "exact-main-generic-gate")
                result[(k, shard_bytes)] = candidate
    return [result[key] for key in sorted(result)]


def derive_survivors(plan_root: Path, manifests: list[dict[str, Any]],
                     references: list[dict[str, Any]],
                     scopes: list[dict[str, Any]]) -> dict[str, Any]:
    plan = validate_plan(plan_root)
    require(len(manifests) == len(references) == len(scopes) and manifests,
            "gate manifests, references, and evidence scopes differ")
    scope = validate_evidence_scope(scopes[0])
    for candidate_scope in scopes[1:]:
        validate_evidence_scope(candidate_scope)
    require(all(candidate == scope for candidate in scopes),
            "gate manifests do not share one normalized evidence environment")
    planned = planned_gate_cells(plan_root, plan)
    evaluated = []
    seen = set()
    candidate_commit = None
    for manifest in manifests:
        campaign = manifest.get("campaign")
        identities = manifest.get("identities")
        analysis = manifest.get("analysis")
        require(isinstance(campaign, dict) and isinstance(identities, dict) and
                isinstance(analysis, dict), "exact-main manifest is incomplete")
        require(campaign.get("candidate_mode") == "generic" and
                campaign.get("batch") == 1 and campaign.get("reuse") == 8 and
                campaign.get("iterations") == 9 and campaign.get("warmup") == 2 and
                campaign.get("threads") == 1 and campaign.get("rounds") == 3 and
                campaign.get("order") == list(EXACT_ORDER),
                "exact-main gate campaign semantics differ")
        baseline = identities.get("baseline_source", {})
        candidate = identities.get("candidate_source", {})
        require(baseline.get("head") == EXACT_MAIN_COMMIT and
                isinstance(candidate.get("head"), str),
                "exact-main source identity differs")
        if candidate_commit is None:
            candidate_commit = candidate["head"]
        require(candidate["head"] == candidate_commit,
                "gate manifests use different candidate commits")
        require(scope["sources"]["candidate"]["head"] == candidate_commit and
                scope["sources"]["baseline"]["head"] == EXACT_MAIN_COMMIT,
                "normalized evidence scope names a different candidate commit")
        raw_cells = campaign.get("cells")
        require(isinstance(raw_cells, list), "gate campaign cells are absent")
        for raw_cell in raw_cells:
            cell = manifest_cell(raw_cell)
            key = (cell["K"], cell["shard_bytes"])
            require(key not in seen and 5 <= cell["K"] <= 128 and
                    cell["R"] == cell["K"] and
                    cell["loss_count"] == cell["K"] and
                    cell["shard_bytes"] in EXACT_GATE_BYTES,
                    "gate manifest cell is duplicated or outside the refinement domain")
            seen.add(key)
            metrics = analysis.get(cell["identifier"])
            require(isinstance(metrics, dict), "gate analysis cell is absent")
            lower = {}
            for metric in ("decode_first_use", "decode_reuse_amortized"):
                result = metrics.get(metric)
                require(isinstance(result, dict) and
                        isinstance(result.get("ci95_lower"), (int, float)) and
                        not isinstance(result.get("ci95_lower"), bool) and
                        result.get("ratio_orientation") ==
                            "baseline_time_over_candidate_time",
                        "gate analysis metric is invalid")
                lower[metric] = float(result["ci95_lower"])
            evaluated.append({
                **cell, "ci95_lower": lower,
                "passes_gate": all(value >= PROMOTION_GATE for value in lower.values()),
            })
    require(set(planned).issubset(seen), "verified manifests omit planned gate cells")
    planned_evaluated = [
        cell for cell in evaluated
        if (cell["K"], cell["shard_bytes"]) in planned
    ]
    requested_refinements = {
        (cell["K"], cell["shard_bytes"])
        for cell in refinement_cells(planned_evaluated)
    }
    supplemental = seen - set(planned)
    require(supplemental.issubset(requested_refinements),
            "verified manifests contain an unrequested supplemental gate cell")
    evaluated.sort(key=lambda item: (item["shard_bytes"], item["K"]))
    survivors = [item for item in evaluated if item["passes_gate"]]
    rejected = [item for item in evaluated if not item["passes_gate"]]
    return signed({
        "schema": SURVIVOR_SCHEMA,
        "plan_content_sha256": (load_json(plan_root / "plan.json"))["content_sha256"],
        "exact_main_commit": EXACT_MAIN_COMMIT,
        "candidate_commit": candidate_commit,
        "promotion_minimum_ci95_lower": PROMOTION_GATE,
        "evidence_scope": scope,
        "evidence_scope_sha256": canonical_sha256(scope),
        "gate_manifests": references,
        "evaluated_cells": evaluated,
        "survivor_cells": survivors,
        "rejected_cells": rejected,
        "required_refinement_cells": refinement_cells(evaluated),
    })


def select_survivors(plan_root: Path, paths: list[Path], output: Path) -> dict[str, Any]:
    require(not output.exists(), f"refusing to replace {output}")
    manifests = []
    scopes = []
    references = []
    for path in paths:
        resolved = path.resolve(strict=True)
        manifest, scope, snapshot = verify_exact_manifest(resolved)
        manifests.append(manifest)
        scopes.append(scope)
        references.append({
            "path": str(resolved), **snapshot,
            "payload_digest": manifest["digest"],
        })
    result = derive_survivors(plan_root, manifests, references, scopes)
    write_json(output, result)
    return result


def validate_survivors(plan_root: Path, path: Path,
                       manifests: list[dict[str, Any]] | None = None,
                       scopes: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    retained = load_json(path)
    value = unsigned(retained, SURVIVOR_SCHEMA, "survivor set")
    require(set(value) == {
        "schema", "plan_content_sha256", "exact_main_commit", "candidate_commit",
        "promotion_minimum_ci95_lower", "gate_manifests", "evaluated_cells",
        "survivor_cells", "rejected_cells", "required_refinement_cells",
        "evidence_scope", "evidence_scope_sha256",
    }, "survivor fields differ")
    references = value["gate_manifests"]
    if manifests is None:
        manifests = []
        scopes = []
        require(isinstance(references, list) and references,
                "survivor gate references are absent")
        for reference in references:
            require(isinstance(reference, dict) and set(reference) == {
                "path", "size", "sha256", "payload_digest",
            }, "survivor gate reference differs")
            manifest_path = Path(reference["path"])
            manifest, scope, snapshot = verify_exact_manifest(manifest_path)
            require(snapshot == {
                        "size": reference["size"],
                        "sha256": reference["sha256"],
                    }, "survivor gate manifest bytes changed")
            require(manifest["digest"] == reference["payload_digest"],
                    "survivor gate payload changed")
            manifests.append(manifest)
            scopes.append(scope)
    require(scopes is not None, "verified evidence scopes are absent")
    expected = derive_survivors(plan_root, manifests, references, scopes)
    require(canonical_bytes(retained) == canonical_bytes(expected),
            "survivor set is not the deterministic verified result")
    return value


def forced_case(k: int, shard_bytes: int, backend: str,
                control: str, candidate: str, phase: str) -> dict[str, Any]:
    return {
        "name": f"{phase}-k{k}-b{shard_bytes}-{backend}-{control[0]}v{candidate[0]}",
        "K": k, "R": k, "profile": "legacy_high_v1", "field": "gf8",
        "backend": backend, "shard_bytes": shard_bytes, "loss_count": k,
        "batch": 1, "reuse": 8, "iterations": 9, "warmup": 2,
        "seed": deterministic_seed(
            phase, k, shard_bytes, BACKENDS.index(backend),
            common.MODES.index(control), common.MODES.index(candidate)),
        "control_mode": control, "candidate_mode": candidate,
    }


def forced_backends_for_scope(scope: dict[str, Any]) -> tuple[str, ...]:
    backend = scope.get("resolved_auto_backend")
    require(backend in BACKENDS, "evidence backend has no forced-runner tier")
    return BACKENDS[:BACKENDS.index(backend) + 1]


def forced_matrix(name: str, cases: list[dict[str, Any]]) -> dict[str, Any]:
    result = {"schema": common.MATRIX_SCHEMA, "name": name, "cases": cases}
    common.normalize_matrix(result)
    return result


def rejection_cells(counts: Iterable[int]) -> list[dict[str, Any]]:
    cells = []
    for k in sorted(set(counts)):
        for shard_bytes in (4096, 65536):
            for loss in sorted({1, 4, k // 2, k - 1}):
                cells.append(main_cell(
                    f"reject-k{k}-r{k}-b{shard_bytes}-l{loss}",
                    k, k, shard_bytes, loss, "loss-rejection"))
            for r in (k - 1, k // 2):
                cells.append(main_cell(
                    f"reject-k{k}-r{r}-b{shard_bytes}-l{r}",
                    k, r, shard_bytes, r, "rate-rejection"))
    return cells


def attestation_case(identifier: str, category: str, k: int, r: int,
                     shard_bytes: int, losses: int, seed: int,
                     expects_balanced: bool) -> dict[str, Any]:
    require(re.fullmatch(r"[a-z0-9][a-z0-9-]{0,127}", identifier) is not None,
            f"invalid attestation identifier {identifier!r}")
    require(category in {
        "survivor_gate", "rejected_gate", "loss_rate_neighbor",
        "extra_aligned_confirmation", "true_tail",
    }, f"invalid attestation category {category!r}")
    require(0 < r <= k and 0 < losses <= r and shard_bytes > 0 and
            0 < seed < 1 << 64, f"invalid attestation cell {identifier}")
    cell = {
        "identifier": identifier, "category": category, "K": k, "R": r,
        "profile": "legacy_high_v1", "field": "gf8", "backend": "auto",
        "shard_bytes": shard_bytes, "loss_count": losses, "batch": 1,
        "reuse": 1, "iterations": 1, "warmup": 0, "threads": 1,
        "seed": seed,
    }
    return {
        "cell": cell,
        "benchmark_arguments": [
            "--k", str(k), "--r", str(r), "--profile", "high",
            "--field", "gf8", "--backend", "auto", "--bytes",
            str(shard_bytes), "--loss", str(losses), "--batch", "1",
            "--reuse", "1", "--iterations", "1", "--warmup", "0",
            "--threads", "1", "--seed", str(seed), "--skip-legacy",
            "--retain-samples", "--report-decode-path", "--json", "OUTPUT",
        ],
        "expected_selected_decode_path": (
            "generic" if expects_balanced else "not_generic"),
        "expected_selected_decode_rule": (
            "balanced_generic" if expects_balanced else "not_balanced_generic"),
    }


def path_attestation(survivor_signed: dict[str, Any]) -> dict[str, Any]:
    """Build the exact AUTO policy implied by the externally proven cells.

    A transform bucket, a surviving K, or an aligned byte size is not evidence
    for another cell.  Only an exact passing gate pair may select the balanced
    generic rule.  All other cells remain negative controls until a later,
    separately authenticated promotion consumes their own evidence.
    """
    survivor = unsigned(survivor_signed, SURVIVOR_SCHEMA, "survivor set")
    cases = []
    survivor_pairs = {
        (cell["K"], cell["shard_bytes"])
        for cell in survivor["survivor_cells"]
    }
    for cell in survivor["evaluated_cells"]:
        expects = (cell["K"], cell["shard_bytes"]) in survivor_pairs
        cases.append(attestation_case(
            ("survivor-" if expects else "rejected-") + cell["identifier"],
            "survivor_gate" if expects else "rejected_gate",
            cell["K"], cell["R"], cell["shard_bytes"],
            cell["loss_count"], cell["seed"], expects))

    survivor_k = sorted({cell["K"] for cell in survivor["survivor_cells"]})
    for cell in rejection_cells(survivor_k):
        cases.append(attestation_case(
            "neighbor-" + cell["identifier"], "loss_rate_neighbor",
            cell["K"], cell["R"], cell["shard_bytes"],
            cell["loss_count"], cell["seed"], False))
    for k in survivor_k:
        for shard_bytes in ALIGNED_CONFIRM_BYTES:
            if shard_bytes in EXACT_GATE_BYTES:
                continue
            cases.append(attestation_case(
                f"extra-aligned-k{k}-b{shard_bytes}",
                "extra_aligned_confirmation", k, k, shard_bytes, k,
                deterministic_seed("auto-extra-aligned", k, shard_bytes), False))
        for shard_bytes in TRUE_TAIL_BYTES:
            cases.append(attestation_case(
                f"true-tail-k{k}-b{shard_bytes}", "true_tail",
                k, k, shard_bytes, k,
                deterministic_seed("auto-true-tail", k, shard_bytes), False))

    cases.sort(key=lambda item: (
        item["cell"]["K"], item["cell"]["R"],
        item["cell"]["shard_bytes"], item["cell"]["loss_count"],
        item["cell"]["seed"], item["cell"]["identifier"]))
    identifiers = [item["cell"]["identifier"] for item in cases]
    identities = [tuple(item["cell"][key] for key in (
        "K", "R", "shard_bytes", "loss_count", "seed")) for item in cases]
    require(len(set(identifiers)) == len(identifiers),
            "attestation identifiers are duplicated")
    require(len(set(identities)) == len(identities),
            "attestation cells are duplicated")
    return signed({
        "schema": ATTESTATION_SCHEMA,
        "name": "balanced-auto-exact-cell-selected-path-attestation",
        "purpose": (
            "Authenticate AUTO path selection for exact externally proven cells "
            "and fail closed for every rejected gate, loss/rate neighbor, "
            "unproven aligned size, and true tail."),
        "candidate_commit": survivor["candidate_commit"],
        "survivor_content_sha256": survivor_signed["content_sha256"],
        "evidence_scope_sha256": survivor["evidence_scope_sha256"],
        "expected_resolved_backend": survivor["evidence_scope"][
            "resolved_auto_backend"],
        "cases": cases,
    })


def materialize_stage(plan_root: Path, survivor_signed: dict[str, Any],
                      root: Path) -> dict[str, Any]:
    require(not root.exists(), f"refusing to replace existing output {root}")
    survivor = unsigned(survivor_signed, SURVIVOR_SCHEMA, "survivor set")
    require(not survivor["required_refinement_cells"],
            "K transition refinement remains; run the listed runner_cell values first")
    root.mkdir(parents=True)
    write_json(root / "survivors.json", survivor_signed)
    survivor_cells = survivor["survivor_cells"]
    survivor_k = sorted({cell["K"] for cell in survivor_cells})
    forced_backends = forced_backends_for_scope(survivor["evidence_scope"])
    artifacts = []

    by_side: dict[int, list[dict[str, Any]]] = {}
    for cell in survivor_cells:
        by_side.setdefault(1 << (cell["K"] - 1).bit_length(), []).append(cell)
    for backend in forced_backends:
        for side, cells in sorted(by_side.items()):
            cases = [
                forced_case(cell["K"], cell["shard_bytes"], backend,
                            control, candidate, "survivor")
                for cell in cells for control, candidate in MODE_PAIRS
            ]
            path = root / "forced-survivors" / f"{backend}-t{side}.json"
            write_json(path, forced_matrix(
                f"balanced-survivors-{backend}-t{side}", cases))
            artifacts.append(artifact(path, root, len(cases), len(cases) * 12,
                                      "forced_surviving_cells"))

    for mode in CONFIRM_MODES:
        cells = [
            main_cell(f"confirm-{mode}-k{k}-b{shard_bytes}",
                      k, k, shard_bytes, k, f"aligned-confirm-{mode}")
            for k in survivor_k for shard_bytes in ALIGNED_CONFIRM_BYTES
        ]
        if cells:
            path = root / "exact-main-confirm" / f"{mode}.json"
            write_json(path, main_document(
                f"balanced-aligned-exact-main-confirm-{mode}", mode, cells,
                "Selected-mode exact-main confirmation at all aligned boundaries.",
                "Run only after the generic gate and forced survivor comparisons pass."))
            artifacts.append(artifact(path, root, len(cells), len(cells) * 12,
                                      "exact_main_aligned_confirmation"))

    rejects = rejection_cells(survivor_k)
    if rejects:
        path = root / "exact-main-auto-rejection-timing.json"
        write_json(path, main_document(
            "balanced-auto-loss-rate-rejection-timing", "auto", rejects,
            "Performance/correctness controls only; this does not attest selected path.",
            "Run after a candidate rule is frozen, alongside path-attestation.json."))
        artifacts.append(artifact(path, root, len(rejects), len(rejects) * 12,
                                  "exact_main_rejection_timing"))

    attestation = path_attestation(survivor_signed)
    attest = root / "path-attestation.json"
    write_json(attest, attestation)
    artifacts.append(artifact(
        attest, root, len(attestation["cases"]), 0,
        "same_binary_selected_path_attestation"))

    for backend in forced_backends:
        for k in survivor_k:
            cases = [
                forced_case(k, shard_bytes, backend, control, candidate, "tail")
                for shard_bytes in TRUE_TAIL_BYTES
                for control, candidate in MODE_PAIRS
            ]
            path = root / "forced-true-tails" / f"{backend}-k{k}.json"
            write_json(path, forced_matrix(
                f"balanced-true-tails-{backend}-k{k}", cases))
            artifacts.append(artifact(path, root, len(cases), len(cases) * 12,
                                      "forced_non_aligned_tail"))

    stage = signed({
        "schema": STAGE_SCHEMA,
        "name": "balanced-survivor-confirmation-stage",
        "plan_content_sha256": (load_json(plan_root / "plan.json"))["content_sha256"],
        "survivor_content_sha256": survivor_signed["content_sha256"],
        "candidate_commit": survivor["candidate_commit"],
        "evidence_scope_sha256": survivor["evidence_scope_sha256"],
        "expected_resolved_backend": survivor["evidence_scope"][
            "resolved_auto_backend"],
        "forced_confirmation_backends": list(forced_backends),
        "excluded_backends": survivor["evidence_scope"]["excluded_backends"],
        "survivor_K": survivor_k,
        "surviving_gate_cell_count": len(survivor_cells),
        "true_tail_bytes": list(TRUE_TAIL_BYTES),
        "aligned_exact_main_bytes": list(ALIGNED_CONFIRM_BYTES),
        "artifacts": artifacts,
        "timed_case_count": sum(item["case_count"] for item in artifacts
                                if item["timed_child_count"] != 0),
        "timed_child_count": sum(item["timed_child_count"] for item in artifacts),
        "path_attestation_case_count": len(attestation["cases"]),
        "path_attestation_survivor_case_count": sum(
            item["expected_selected_decode_path"] == "generic"
            for item in attestation["cases"]),
        "promotion_requires_path_attestation": bool(attestation["cases"]),
    })
    write_json(root / "stage.json", stage)
    return stage


def validate_stage(plan_root: Path, root: Path,
                   manifests: list[dict[str, Any]] | None = None,
                   scopes: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    retained_survivor = load_json(root / "survivors.json")
    survivor_path = root / "survivors.json"
    survivor = validate_survivors(plan_root, survivor_path, manifests, scopes)
    stage = unsigned(load_json(root / "stage.json"), STAGE_SCHEMA, "stage")
    require(stage["candidate_commit"] == survivor["candidate_commit"],
            "stage candidate commit differs")
    with tempfile.TemporaryDirectory(prefix="leopard2-stage-expected-") as temporary:
        expected = Path(temporary) / "stage"
        materialize_stage(plan_root, retained_survivor, expected)
        compare_trees(root, expected, "stage")
    return stage


def stable_file_identity(value: dict[str, Any], label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} identity is not an object")
    relative = value.get("relative_path")
    require(relative is None or isinstance(relative, str),
            f"{label} relative path is invalid")
    digest = value.get("sha256")
    require(isinstance(digest, str) and
            re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
            f"{label} digest is invalid")
    require(type(value.get("size")) is int and value["size"] >= 0 and
            type(value.get("mode")) is int and value["mode"] >= 0,
            f"{label} size/mode is invalid")
    return {
        "relative_path": relative, "sha256": digest,
        "size": value["size"], "mode": value["mode"],
    }


def stable_source_identity(value: dict[str, Any], candidate_commit: str) -> dict[str, Any]:
    require(isinstance(value, dict) and
            value.get("head") == candidate_commit and
            value.get("expected_commit") == candidate_commit and
            value.get("status") == "clean" and
            value.get("status_sha256") == common.EMPTY_SHA256 and
            re.fullmatch(r"[0-9a-f]{40}", str(value.get("tree"))) is not None,
            "source identity is not the clean candidate commit")
    return {
        "head": candidate_commit, "tree": value["tree"],
        "status": "clean", "status_sha256": common.EMPTY_SHA256,
    }


def stable_build_identity(value: dict[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict), "build identity is not an object")
    sources = value.get("sources")
    objects = value.get("objects")
    require(isinstance(sources, dict) and set(sources) == {
        "benchmark", "decoder", "dispatch",
    }, "build source identity set differs")
    require(isinstance(objects, dict) and set(objects) == {
        "benchmark", "decoder",
    }, "build object identity set differs")
    return {
        "cache": stable_file_identity(value.get("cache"), "CMake cache"),
        "graph": stable_file_identity(value.get("graph"), "build graph"),
        "sources": {
            key: stable_file_identity(sources[key], f"{key} source")
            for key in sorted(sources)
        },
        "objects": {
            key: stable_file_identity(objects[key], f"{key} object")
            for key in sorted(objects)
        },
        "archive": stable_file_identity(value.get("archive"), "archive"),
        "binary": stable_file_identity(value.get("binary"), "benchmark binary"),
    }


def load_exact_main_runner() -> Any:
    path = Path(__file__).resolve().parents[1] / "main_compare" / "run_abba.py"
    name = "leopard2_exact_main_build_validator"
    module = sys.modules.get(name)
    if module is not None:
        return module
    module_spec = importlib.util.spec_from_file_location(name, path)
    require(module_spec is not None and module_spec.loader is not None,
            "cannot load exact-main build validator")
    module = importlib.util.module_from_spec(module_spec)
    sys.modules[name] = module
    module_spec.loader.exec_module(module)
    return module


def canonical_candidate_build_identity(
    source_root: Path, binary_relative: str
) -> dict[str, Any]:
    """Run the exact-main schema-v5 compile/archive/link semantic validator."""
    source_root = source_root.resolve(strict=True)
    relative = Path(binary_relative)
    require(not relative.is_absolute() and ".." not in relative.parts,
            "canonical benchmark path must be source-root-relative")
    binary = (source_root / relative).resolve(strict=True)
    build_root = common.find_build_root(binary).resolve(strict=True)
    archive = (build_root / "libleopard.a").resolve(strict=True)
    runner = load_exact_main_runner()
    specification = {
        "candidate_build_dir": str(build_root),
        "candidate_executable": str(binary),
        "candidate_archive": str(archive),
        "candidate_source_root": str(source_root),
        # Candidate validation resolves this root but does not consume baseline
        # sources.  Supplying the same clean root avoids inventing a second tree.
        "baseline_source_root": str(source_root),
    }
    try:
        provenance = runner.build_provenance(
            "candidate", specification, runner.RAW_SCHEMA)
    except Exception as error:
        raise PlanError(
            f"canonical Release/AUTO build validation failed: {error}") from error
    cache = provenance.get("validated_cache")
    semantics = provenance.get("validated_compile_commands")
    require(isinstance(cache, dict) and all(
        cache.get(key) == expected
        for key, expected in REQUIRED_CANDIDATE_CACHE.items()) and
            isinstance(semantics, dict) and
            semantics.get("schema") == COMPILE_COMMANDS_SCHEMA and
            semantics.get("implementation") == "candidate" and
            semantics.get("profile") == CANDIDATE_COMPILE_PROFILE and
            semantics.get("validated_optimization") == "-O3" and
            semantics.get("validated_openmp") is True and
            isinstance(semantics.get("required_entries"), list) and
            semantics["required_entries"] and
            isinstance(provenance.get("archive_link_recipe_content"), dict),
            "canonical build provenance does not bind Release AUTO semantics")
    normalized = _normalize_bound_paths(
        provenance,
        ((str(source_root), "$SOURCE"), (str(build_root), "$BUILD")),
    )
    require(isinstance(normalized, dict), "canonical build scope is invalid")
    return {
        "schema": "leopard2-canonical-production-build/v1",
        "validator": "exact-main/run_abba.py build_provenance schema v5",
        "provenance": normalized,
        "provenance_sha256": canonical_sha256(normalized),
    }


def xorshift64(value: int) -> int:
    mask = (1 << 64) - 1
    value ^= (value << 13) & mask
    value ^= value >> 7
    value ^= (value << 17) & mask
    return value & mask


def expected_missing_indices(k: int, losses: int, seed: int) -> list[int]:
    order = list(range(k))
    state = seed ^ 0xd1b54a32d192ed03
    if state == 0:
        state = 0x9e3779b97f4a7c15
    for remaining in range(len(order), 1, -1):
        state = xorshift64(state)
        selected = state % remaining
        order[remaining - 1], order[selected] = (
            order[selected], order[remaining - 1])
    return sorted(order[:losses])


def ceil_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def coherent_production_decode_pair(path: object, rule: object) -> bool:
    return type(path) is str and type(rule) is str and \
        (path, rule) in PRODUCTION_DECODE_PATH_RULE_PAIRS


def validate_selector_pair_contract(source_root: Path) -> None:
    """Prove the Python whitelist still matches the production C++ emitters."""
    dispatch = (source_root / "Leopard2Dispatch.h").read_text(encoding="utf-8")
    decoder = (source_root / "leopard2.cpp").read_text(encoding="utf-8")
    path_names = dict(re.findall(
        r'case kDecodePath([A-Za-z0-9_]+): return "([a-z0-9_]+)";', dispatch))
    rule_names = dict(re.findall(
        r'case kDecodeRule([A-Za-z0-9_]+): return "([a-z0-9_]+)";', dispatch))
    token_pairs = set(re.findall(
        r'selection\.path = kDecodePath([A-Za-z0-9_]+);\s*'
        r'selection\.rule = kDecodeRule([A-Za-z0-9_]+);', dispatch))
    token_pairs.update(re.findall(
        r'FillTerminalDecodePathInfo\(kDecodePath([A-Za-z0-9_]+),\s*'
        r'kDecodeRule([A-Za-z0-9_]+),', decoder))
    require(token_pairs and all(path in path_names and rule in rule_names
                                for path, rule in token_pairs),
            "cannot derive production decode path/rule names")
    derived = {(path_names[path], rule_names[rule]) for path, rule in token_pairs}
    require(derived == PRODUCTION_DECODE_PATH_RULE_PAIRS,
            "production decode path/rule pair contract changed")


def expected_auto_decode_pair(case: dict[str, Any], backend: str) -> tuple[str, str]:
    """Mirror terminal-plan and AUTO selector ordering for one canonical cell.

    These attestation cells are legacy-high GF8, one-item, plan-known calls with
    no force flags and at least one missing original.  The only terminal plan is
    bounded direct repair.  The remaining prediction follows the ordered
    balanced, measured-materialized, and workspace rules in
    Leopard2Dispatch.h.  A passing gate is the *only* input allowed to satisfy
    the balanced predicate; this is intentionally narrower than an unverified
    current dispatcher rule.
    """
    cell = case["cell"]
    require(backend in AUTO_BACKENDS, "AUTO attestation backend is unknown")
    k = cell["K"]
    r = cell["R"]
    losses = cell["loss_count"]
    if 2 <= k <= DIRECT_MAX_ORIGINALS and losses <= DIRECT_MAX_LOSSES:
        return "direct", "direct"
    if case["expected_selected_decode_path"] == "generic":
        return "generic", "balanced_generic"

    padded = ceil_power_of_two(r)
    parent = ceil_power_of_two(k + padded)
    rounded = (cell["shard_bytes"] + 63) & ~63
    measured_materialized = (
        k == 224 and r == 32 and padded == 32 and parent == 256 and
        0 < losses <= 8 and rounded <= 64 * 1024 and (
            (backend in {"avx2", "avx512"} and rounded >= 24 * 1024) or
            (backend == "ssse3" and rounded >= 32 * 1024)))
    if measured_materialized:
        return "materialized", "measured_materialized"
    tiled_work_slots = 2 * padded + losses
    if tiled_work_slots < parent:
        return "tiled", "workspace_tiled"
    return "materialized", "workspace_materialized"


def expected_required_work_slots(case: dict[str, Any], path: str) -> int:
    """Mirror GetDecodePlanPathInfo workspace geometry for production paths."""
    cell = case["cell"]
    if cell["profile"] == "legacy_high_v1":
        padded = ceil_power_of_two(cell["R"])
        parent = ceil_power_of_two(cell["K"] + padded)
    else:
        require(cell["profile"] == "low_v1", "unknown decode profile")
        padded = ceil_power_of_two(cell["K"])
        parent = ceil_power_of_two(padded + cell["R"])
    if path in {"no_op", "direct"}:
        return 0
    if path in {"generic", "materialized"}:
        return parent
    require(path == "tiled", "cannot derive workspace for unknown decode path")
    # All campaign cells use legacy_high_v1.  Low-profile tiled execution uses
    # 2*P, while legacy high additionally retains one locator/output slot for
    # each requested missing original.
    if cell["profile"] == "legacy_high_v1":
        return 2 * padded + cell["loss_count"]
    return 2 * padded


def validate_attestation_output(document: object, case: dict[str, Any],
                                expected_backend: str | None = None) -> dict[str, Any]:
    require(isinstance(document, dict) and
            document.get("schema") == common.BENCHMARK_SCHEMA and
            set(document) == {
                "schema", "build", "parameters", "resolved", "correctness",
                "workload_digests", "memory", "metrics", "legacy",
            }, "attestation benchmark top-level shape differs")
    parameters = document.get("parameters")
    resolved = document.get("resolved")
    correctness = document.get("correctness")
    build_output = document.get("build")
    require(isinstance(build_output, dict) and set(build_output) == {
        "compiler", "compiler_version", "cplusplus",
    } and isinstance(build_output["compiler"], str) and
            isinstance(build_output["compiler_version"], str) and
            type(build_output["cplusplus"]) is int and build_output["cplusplus"] > 0,
            "attestation benchmark build record differs")
    require(isinstance(parameters, dict) and isinstance(resolved, dict) and
            isinstance(correctness, dict), "attestation benchmark is incomplete")
    require(set(parameters) == {
        "K", "R", "requested_profile", "requested_field", "requested_backend",
        "force_generic_decode", "force_specialized_decode", "force_tiled_decode",
        "force_materialized_decode", "skip_legacy", "retain_samples",
        "report_decode_path", "shard_bytes", "loss_count",
        "missing_original_indices", "batch", "reuse", "iterations", "warmup",
        "thread_count", "seed",
    }, "attestation benchmark parameter shape differs")
    cell = case["cell"]
    expected_parameters = {
        "K": cell["K"], "R": cell["R"],
        "requested_profile": "legacy_high_v1", "requested_field": "gf8",
        "requested_backend": "auto", "force_generic_decode": False,
        "force_specialized_decode": False, "force_tiled_decode": False,
        "force_materialized_decode": False, "skip_legacy": True,
        "retain_samples": True, "report_decode_path": True,
        "shard_bytes": cell["shard_bytes"], "loss_count": cell["loss_count"],
        "missing_original_indices": expected_missing_indices(
            cell["K"], cell["loss_count"], cell["seed"]),
        "batch": 1, "reuse": 1, "iterations": 1, "warmup": 0,
        "thread_count": 1, "seed": cell["seed"],
    }
    for key, expected in expected_parameters.items():
        actual = parameters.get(key)
        if key == "missing_original_indices":
            require(isinstance(actual, list) and
                    all(type(item) is int for item in actual) and actual == expected,
                    f"attestation benchmark {key} differs for {cell['identifier']}")
        else:
            require(type(actual) is type(expected) and actual == expected,
                    f"attestation benchmark {key} differs for {cell['identifier']}")
    transform = ceil_power_of_two(cell["R"])
    parent = ceil_power_of_two(cell["K"] + transform)
    require(set(resolved) == {
        "profile", "field", "backend", "thread_count", "parent_count",
        "padded_side", "selected_decode_path", "selected_decode_rule",
        "decode_required_work_slots", "decode_aligned_prefix_bytes",
        "decode_tail_bytes", "decode_rounded_bytes", "decode_multi_item_batch",
    } and resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and
            type(resolved.get("backend")) is str and
            resolved.get("backend") in AUTO_BACKENDS and
            (expected_backend is None or resolved.get("backend") == expected_backend) and
            type(resolved.get("thread_count")) is int and
            resolved.get("thread_count") == 1 and
            type(resolved.get("parent_count")) is int and
            resolved.get("parent_count") == parent and
            type(resolved.get("padded_side")) is int and
            resolved.get("padded_side") == transform and
            type(resolved.get("decode_aligned_prefix_bytes")) is int and
            resolved.get("decode_aligned_prefix_bytes") ==
                cell["shard_bytes"] & ~63 and
            type(resolved.get("decode_tail_bytes")) is int and
            resolved.get("decode_tail_bytes") == cell["shard_bytes"] & 63 and
            type(resolved.get("decode_rounded_bytes")) is int and
            resolved.get("decode_rounded_bytes") ==
                (cell["shard_bytes"] + 63) & ~63 and
            resolved.get("decode_multi_item_batch") is False,
            f"attestation resolved codec identity differs for {cell['identifier']}")
    selected_path = resolved["selected_decode_path"]
    selected_rule = resolved["selected_decode_rule"]
    require(coherent_production_decode_pair(selected_path, selected_rule),
            f"attestation selected path/rule pair is not a production pair: "
            f"{cell['identifier']}")
    expected_pair = expected_auto_decode_pair(case, resolved["backend"])
    require((selected_path, selected_rule) == expected_pair,
            f"attestation selected path/rule differs for {cell['identifier']}: "
            f"expected {expected_pair!r}, got {(selected_path, selected_rule)!r}")
    required_slots = expected_required_work_slots(case, selected_path)
    require(type(resolved.get("decode_required_work_slots")) is int and
            resolved["decode_required_work_slots"] == required_slots,
            f"attestation required work slots differ for {cell['identifier']}: "
            f"expected {required_slots}, got "
            f"{resolved.get('decode_required_work_slots')!r}")
    require(isinstance(correctness, dict) and set(correctness) == {
        "leopard2_round_trip", "legacy_comparison",
    } and correctness.get("leopard2_round_trip") is True and
            correctness.get("legacy_comparison") is None,
            f"attestation round trip differs for {cell['identifier']}")
    digests = document.get("workload_digests")
    require(isinstance(digests, dict) and set(digests) == {
        "algorithm", "original_data", "transmitted_parity", "recovered_originals",
    } and digests.get("algorithm") == "fnv1a64" and all(
        re.fullmatch(r"[0-9a-f]{16}", str(digests.get(key))) is not None
        for key in ("original_data", "transmitted_parity", "recovered_originals")),
        f"attestation workload digests differ for {cell['identifier']}")
    memory = document.get("memory")
    require(isinstance(memory, dict) and set(memory) == {
        "scratch_alignment", "encode_scratch_bytes_per_stripe",
        "decode_scratch_bytes_per_stripe", "encode_scratch_bytes_batch",
        "decode_scratch_bytes_batch",
    } and all(type(value) is int and value >= 0 for value in memory.values()) and
            memory["scratch_alignment"] > 0 and
            memory["encode_scratch_bytes_per_stripe"] ==
                memory["encode_scratch_bytes_batch"] and
            memory["decode_scratch_bytes_per_stripe"] ==
                memory["decode_scratch_bytes_batch"],
            f"attestation memory record differs for {cell['identifier']}")
    metrics = document.get("metrics")
    require(isinstance(metrics, dict) and set(metrics) == {
        "codec_setup", "encode_execution", "decode_plan_setup",
        "decode_execution", "decode_amortized_at_reuse", "rate_semantics",
    }, f"attestation metric shape differs for {cell['identifier']}")
    for label, key in (("codec setup", "codec_setup"),
                       ("decode plan setup", "decode_plan_setup")):
        summary = metrics[key]
        require(isinstance(summary, dict) and
                isinstance(summary.get("samples_us"), list) and
                len(summary["samples_us"]) == 1 and
                isinstance(summary["samples_us"][0], (int, float)) and
                not isinstance(summary["samples_us"][0], bool) and
                summary["samples_us"][0] >= 0,
                f"attestation {label} sample differs for {cell['identifier']}")
    for label, key in (("encode", "encode_execution"),
                       ("decode", "decode_execution")):
        summary = metrics[key]
        require(isinstance(summary, dict) and
                isinstance(summary.get("samples_us_per_batch_call"), list) and
                len(summary["samples_us_per_batch_call"]) == 1 and
                isinstance(summary["samples_us_per_batch_call"][0], (int, float)) and
                not isinstance(summary["samples_us_per_batch_call"][0], bool) and
                summary["samples_us_per_batch_call"][0] > 0,
                f"attestation {label} sample differs for {cell['identifier']}")
    amortized = metrics["decode_amortized_at_reuse"]
    require(isinstance(amortized, dict) and
            type(amortized.get("reuse_count")) is int and
            amortized.get("reuse_count") == 1 and
            isinstance(metrics["rate_semantics"], str),
            f"attestation amortization record differs for {cell['identifier']}")
    require(document.get("legacy") == {
        "available": False, "unavailable_reason": "disabled by --skip-legacy",
        "codec_setup": None, "decode_timing_includes_setup": True,
        "encode_execution": None, "decode_including_setup": None,
    }, f"attestation legacy-skip record differs for {cell['identifier']}")
    return {
        "resolved_backend": resolved["backend"],
        "selected_decode_path": selected_path,
        "selected_decode_rule": selected_rule,
        "workload_digests_sha256": canonical_sha256(digests),
        "benchmark_payload_sha256": canonical_sha256(document),
    }


def attestation_identities(source_root: Path, candidate_commit: str,
                           binary_relative: str) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    validate_selector_pair_contract(source_root)
    source = stable_source_identity(
        common.git_identity(source_root, candidate_commit), candidate_commit)
    build = {
        "artifact_closure": stable_build_identity(
            common.build_identity(source_root, binary_relative)),
        "canonical_production": canonical_candidate_build_identity(
            source_root, binary_relative),
    }
    collector = stable_file_identity(
        common.file_identity(Path(__file__), source_root, "attestation collector"),
        "attestation collector")
    require(collector["relative_path"] ==
            "experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py",
            "attestation collector path is not canonical")
    return source, build, collector


def derive_attestation_result(
    stage_root: Path, source: dict[str, Any], build: dict[str, Any],
    collector: dict[str, Any], raw_documents: dict[str, object],
    raw_artifacts: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    stage_signed = load_json(stage_root / "stage.json")
    stage = unsigned(stage_signed, STAGE_SCHEMA, "stage")
    spec_signed = load_json(stage_root / "path-attestation.json")
    spec = unsigned(spec_signed, ATTESTATION_SCHEMA, "path attestation")
    survivor_signed = load_json(stage_root / "survivors.json")
    survivor = unsigned(survivor_signed, SURVIVOR_SCHEMA, "survivor set")
    require(spec["candidate_commit"] == stage["candidate_commit"] and
            spec["survivor_content_sha256"] == stage["survivor_content_sha256"] and
            survivor_signed["content_sha256"] == stage["survivor_content_sha256"] and
            spec["evidence_scope_sha256"] == survivor["evidence_scope_sha256"] and
            stage["evidence_scope_sha256"] == survivor["evidence_scope_sha256"] and
            spec["expected_resolved_backend"] ==
                survivor["evidence_scope"]["resolved_auto_backend"] and
            stage["expected_resolved_backend"] == spec["expected_resolved_backend"] and
            spec["expected_resolved_backend"] in CAMPAIGN_BACKENDS,
            "attestation specification is stale")
    require(isinstance(source, dict) and set(source) == {
        "head", "tree", "status", "status_sha256",
    } and source.get("head") == stage["candidate_commit"] and
            re.fullmatch(r"[0-9a-f]{40}", str(source.get("tree"))) is not None and
            source.get("status") == "clean" and
            source.get("status_sha256") == common.EMPTY_SHA256,
            "attestation source is not the exact clean candidate commit")
    require(isinstance(build, dict) and set(build) == {
        "artifact_closure", "canonical_production",
    }, "attestation benchmark build identity shape differs")
    artifacts = build.get("artifact_closure")
    canonical = build.get("canonical_production")
    require(isinstance(artifacts, dict) and set(artifacts) == {
        "cache", "graph", "sources", "objects", "archive", "binary",
    } and isinstance(artifacts.get("sources"), dict) and
            set(artifacts["sources"]) == {"benchmark", "decoder", "dispatch"} and
            isinstance(artifacts.get("objects"), dict) and
            set(artifacts["objects"]) == {"benchmark", "decoder"} and
            isinstance(canonical, dict) and set(canonical) == {
                "schema", "validator", "provenance", "provenance_sha256",
            } and canonical["schema"] == "leopard2-canonical-production-build/v1" and
            canonical["validator"] ==
                "exact-main/run_abba.py build_provenance schema v5" and
            isinstance(canonical["provenance"], dict) and
            canonical_sha256(canonical["provenance"]) ==
                canonical["provenance_sha256"] and
            isinstance(canonical["provenance"].get("validated_cache"), dict) and
            all(canonical["provenance"]["validated_cache"].get(key) == expected
                for key, expected in REQUIRED_CANDIDATE_CACHE.items()) and
            isinstance(canonical["provenance"].get(
                "validated_compile_commands"), dict) and
            canonical["provenance"]["validated_compile_commands"].get(
                "schema") == COMPILE_COMMANDS_SCHEMA and
            canonical["provenance"]["validated_compile_commands"].get(
                "implementation") == "candidate" and
            canonical["provenance"]["validated_compile_commands"].get(
                "profile") == CANDIDATE_COMPILE_PROFILE and
            canonical["provenance"]["validated_compile_commands"].get(
                "validated_optimization") == "-O3" and
            canonical["provenance"]["validated_compile_commands"].get(
                "validated_openmp") is True and
            isinstance(canonical["provenance"][
                "validated_compile_commands"].get("required_entries"), list) and
            canonical["provenance"]["validated_compile_commands"][
                "required_entries"] and
            isinstance(canonical["provenance"].get(
                "archive_link_recipe_content"), dict),
            "attestation benchmark build identity shape differs")
    for label, identity in [
        ("CMake cache", artifacts["cache"]),
        ("build graph", artifacts["graph"]),
        ("archive", artifacts["archive"]),
        ("benchmark binary", artifacts["binary"]),
        *[(f"{key} source", value) for key, value in artifacts["sources"].items()],
        *[(f"{key} object", value) for key, value in artifacts["objects"].items()],
        ("collector", collector),
    ]:
        require(isinstance(identity, dict) and set(identity) == {
            "relative_path", "sha256", "size", "mode",
        } and (identity["relative_path"] is None or
               isinstance(identity["relative_path"], str)) and
                re.fullmatch(r"[0-9a-f]{64}", str(identity["sha256"])) is not None and
                type(identity["size"]) is int and identity["size"] >= 0 and
                type(identity["mode"]) is int and identity["mode"] >= 0,
                f"attestation {label} identity differs")
    require(collector["relative_path"] ==
            "experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py",
            "attestation collector identity is not canonical")
    expected_ids = [case["cell"]["identifier"] for case in spec["cases"]]
    require(len(expected_ids) == len(set(expected_ids)) and
            set(raw_documents) == set(expected_ids) and
            set(raw_artifacts) == set(expected_ids),
            "attestation result set is missing, duplicated, or extra")
    records = []
    resolved_backend = spec["expected_resolved_backend"]
    for case in spec["cases"]:
        identifier = case["cell"]["identifier"]
        observed = validate_attestation_output(
            raw_documents[identifier], case, resolved_backend)
        artifact_value = raw_artifacts[identifier]
        require(isinstance(artifact_value, dict) and set(artifact_value) == {
            "path", "size", "sha256",
        } and isinstance(artifact_value["path"], str) and
                type(artifact_value["size"]) is int and artifact_value["size"] > 0 and
                re.fullmatch(r"[0-9a-f]{64}",
                             str(artifact_value["sha256"])) is not None,
                f"attestation raw artifact differs for {identifier}")
        records.append({
            "identifier": identifier, "raw": artifact_value, **observed,
        })
    return signed({
        "schema": ATTESTATION_RESULT_SCHEMA, "status": "complete", "valid": True,
        "plan_content_sha256": stage["plan_content_sha256"],
        "stage_content_sha256": stage_signed["content_sha256"],
        "survivor_content_sha256": stage["survivor_content_sha256"],
        "attestation_content_sha256": spec_signed["content_sha256"],
        "candidate_commit": stage["candidate_commit"],
        "resolved_backend": resolved_backend,
        "source_identity": source, "build_identity": build,
        "collector_identity": collector, "records": records,
    })


def output_is_ignored(path: Path, source_root: Path) -> None:
    try:
        path.resolve().relative_to(source_root.resolve())
    except ValueError:
        return
    completed = subprocess.run(
        ["git", "-C", str(source_root), "check-ignore", "-q", str(path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    require(completed.returncode == 0,
            "attestation output inside source root must be ignored by Git")


def collect_attestation(plan_root: Path, stage_root: Path, source_root: Path,
                        binary_relative: str, output: Path) -> dict[str, Any]:
    stage = validate_stage(plan_root, stage_root)
    candidate_commit = stage["candidate_commit"]
    relative = Path(binary_relative)
    require(not relative.is_absolute() and ".." not in relative.parts,
            "attestation benchmark path must be source-root-relative")
    stable_source_identity(
        common.git_identity(source_root, candidate_commit), candidate_commit)
    require(not output.exists(), f"refusing to replace attestation output {output}")
    output_is_ignored(output, source_root)
    binary = source_root / relative
    common.refresh_build(binary)
    source, build, collector = attestation_identities(
        source_root, candidate_commit, binary_relative)
    output.mkdir(parents=True)
    raw_root = output / "raw"
    raw_root.mkdir()
    spec = unsigned(load_json(stage_root / "path-attestation.json"),
                    ATTESTATION_SCHEMA, "path attestation")
    raw_documents = {}
    raw_artifacts = {}
    for index, case in enumerate(spec["cases"]):
        identifier = case["cell"]["identifier"]
        result_path = raw_root / f"{index:04d}-{identifier}.json"
        arguments = [str(binary)] + [
            str(result_path) if item == "OUTPUT" else item
            for item in case["benchmark_arguments"]
        ]
        completed = subprocess.run(
            arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(completed.returncode == 0 and not completed.stdout and
                not completed.stderr and result_path.is_file(),
                f"attestation benchmark failed for {identifier}")
        document = load_json(result_path)
        validate_attestation_output(
            document, case, spec["expected_resolved_backend"])
        raw_documents[identifier] = document
        raw_artifacts[identifier] = {
            "path": str(result_path.relative_to(output)),
            "size": result_path.stat().st_size, "sha256": file_sha256(result_path),
        }
    final_source, final_build, final_collector = attestation_identities(
        source_root, candidate_commit, binary_relative)
    require((source, build, collector) ==
            (final_source, final_build, final_collector),
            "attestation source or benchmark identity changed during collection")
    result = derive_attestation_result(
        stage_root, source, build, collector, raw_documents, raw_artifacts)
    write_json(output / "manifest.json", result)
    return result


def validate_attestation_result_files(
    stage_root: Path, manifest_path: Path, source: dict[str, Any],
    build: dict[str, Any], collector: dict[str, Any],
) -> dict[str, Any]:
    stage = unsigned(load_json(stage_root / "stage.json"), STAGE_SCHEMA, "stage")
    root = manifest_path.resolve(strict=True).parent
    retained = load_json(manifest_path)
    value = unsigned(retained, ATTESTATION_RESULT_SCHEMA, "attestation result")
    require(set(value) == {
        "schema", "status", "valid", "plan_content_sha256",
        "stage_content_sha256", "survivor_content_sha256",
        "attestation_content_sha256", "candidate_commit", "resolved_backend",
        "source_identity", "build_identity", "collector_identity", "records",
    } and value["status"] == "complete" and value["valid"] is True and
            value["candidate_commit"] == stage["candidate_commit"],
            "attestation result identity differs")
    require(value["source_identity"] == source and
            value["build_identity"] == build and
            value["collector_identity"] == collector,
            "attestation live source or benchmark identity differs")
    records = value["records"]
    require(isinstance(records, list), "attestation records are not a list")
    raw_documents = {}
    raw_artifacts = {}
    expected_files = {Path("manifest.json")}
    for record in records:
        require(isinstance(record, dict) and isinstance(record.get("identifier"), str) and
                isinstance(record.get("raw"), dict),
                "attestation record shape differs")
        identifier = record["identifier"]
        require(identifier not in raw_documents,
                "attestation result identifier is duplicated")
        artifact_value = record["raw"]
        require(set(artifact_value) == {"path", "size", "sha256"} and
                isinstance(artifact_value["path"], str) and
                type(artifact_value["size"]) is int and artifact_value["size"] > 0 and
                re.fullmatch(r"[0-9a-f]{64}",
                             str(artifact_value["sha256"])) is not None,
                f"attestation artifact shape differs for {identifier}")
        relative = Path(artifact_value["path"])
        require(not relative.is_absolute() and ".." not in relative.parts,
                f"attestation artifact path escapes output for {identifier}")
        path = root / relative
        require(path.is_file() and path.stat().st_size == artifact_value["size"] and
                file_sha256(path) == artifact_value["sha256"],
                f"attestation artifact bytes changed for {identifier}")
        raw_documents[identifier] = load_json(path)
        raw_artifacts[identifier] = artifact_value
        expected_files.add(relative)
    actual_files = {
        path.relative_to(root) for path in root.rglob("*") if path.is_file()
    }
    require(actual_files == expected_files,
            "attestation output contains a missing or extra file")
    expected = derive_attestation_result(
        stage_root, source, build, collector, raw_documents, raw_artifacts)
    require(canonical_bytes(retained) == canonical_bytes(expected),
            "attestation result is not the deterministic authenticated result")
    return value


def validate_attestation_result(plan_root: Path, stage_root: Path,
                                manifest_path: Path, source_root: Path,
                                binary_relative: str) -> dict[str, Any]:
    stage = validate_stage(plan_root, stage_root)
    source, build, collector = attestation_identities(
        source_root, stage["candidate_commit"], binary_relative)
    return validate_attestation_result_files(
        stage_root, manifest_path, source, build, collector)


def adversarial_resign(path: Path, mutator) -> None:
    value = load_json(path)
    require(isinstance(value, dict), "mutation target is not an object")
    value.pop("content_sha256", None)
    mutator(value)
    value["content_sha256"] = canonical_sha256(value)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def fake_evidence_scope(backend: str = "avx2") -> dict[str, Any]:
    require(backend in CAMPAIGN_BACKENDS, "fixture backend is outside campaign")
    def fixture_artifact(path: str, kind: str, character: str,
                         mode: int = 0o644) -> dict[str, Any]:
        return {"path": path, "kind": kind, "sha256": character * 64,
                "size": 1, "mode": mode}

    def fixture_external_record(
        operand: str, role: str, kind: str,
    ) -> dict[str, Any]:
        metadata, digest, _, resolved, symlink_chain = \
            common.current_external_file_snapshot(
                Path(operand), f"fixture {role}")
        return {
            "operand": operand, "role": role,
            "lexical_symlink_chain": symlink_chain,
            "artifact": {
                "path": str(resolved), "kind": kind,
                "sha256": digest, "size": metadata.st_size,
                "mode": stat.S_IMODE(metadata.st_mode),
            },
        }

    def fixture_text(text: str) -> dict[str, Any]:
        encoded = text.encode("utf-8")
        return {"encoding": "utf-8", "text": text, "size": len(encoded),
                "sha256": hashlib.sha256(encoded).hexdigest()}

    def fixture_commit_object(raw: bytes) -> tuple[str, dict[str, Any]]:
        head = hashlib.sha1(
            f"commit {len(raw)}\0".encode("ascii") + raw).hexdigest()
        return head, {
            "encoding": "base64", "size": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(), "object_id": head,
            "base64": base64.b64encode(raw).decode("ascii"),
        }

    compiler = fixture_artifact(
        "/usr/bin/compiler", "compiler", "a", 0o755)
    compiler_text = "fixture compiler 1.0\n"
    compiler_version = {
        "sha256": hashlib.sha256(compiler_text.encode("utf-8")).hexdigest(),
        "text": compiler_text,
    }

    def fixture_build(role: str, character: str) -> dict[str, Any]:
        baseline = role == "baseline"
        root = f"${role.upper()}_BUILD"
        source_root = "$BASELINE_SOURCE" if baseline else "$CANDIDATE_SOURCE"
        target = "leopard_main_exact.dir" if baseline else "leopard.dir"
        executable_name = "leopard_main_benchmark" if baseline else "bench_leopard2"
        archive_name = "libleopard_main_exact.a" if baseline else "libleopard.a"
        library_names = (
            BASELINE_LIBRARY_SOURCES if baseline else CANDIDATE_LIBRARY_SOURCES)
        library_pairs = []
        for name in library_names:
            source_path = f"{source_root}/{name}"
            object_relative = _normalized_compile_output(role, source_path)
            library_pairs.append({
                "source": fixture_artifact(
                    source_path, "source_file", character),
                "object": fixture_artifact(
                    f"{root}/{object_relative}",
                    "object_file", character),
            })
        benchmark_source = fixture_artifact(
            ("$CANDIDATE_SOURCE/experiments/leopard2/main_compare/"
             "legacy_main_benchmark.cpp" if baseline else
             "$CANDIDATE_SOURCE/bench/leopard2/benchmark.cpp"),
            "source_file", character)
        benchmark_relative = _normalized_compile_output(
            role, benchmark_source["path"])
        benchmark_object = fixture_artifact(
            f"{root}/{benchmark_relative}",
            "object_file", character)
        pairs = sorted([*library_pairs,
            {"source": benchmark_source, "object": benchmark_object},
        ], key=lambda pair: pair["source"]["path"])
        cache = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
            "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
            **(REQUIRED_BASELINE_CACHE if baseline else REQUIRED_CANDIDATE_CACHE),
        }
        benchmark_relative = benchmark_object["path"][len(root) + 1:]
        objects_by_member = {
            pair["object"]["path"].rsplit("/", 1)[-1]:
                pair["object"]["path"][len(root) + 1:]
            for pair in library_pairs
        }
        members = [Path(name).name + ".o" for name in library_names]
        archive_text = fixture_text(
            f"/usr/bin/ar qc {archive_name} " +
            " ".join(objects_by_member[member] for member in members) + "\n"
            f"/usr/bin/ranlib {archive_name}\n")
        executable_text = fixture_text(
            f"/usr/bin/compiler -O3 {archive_name} {benchmark_relative} "
            f"-o {executable_name} "
            "/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so "
            "/usr/lib/x86_64-linux-gnu/libpthread.a\n")
        archive_recipe = fixture_artifact(
            f"{root}/CMakeFiles/{target}/link.txt", "build_metadata", character)
        archive_recipe.update({
            "size": archive_text["size"], "sha256": archive_text["sha256"]})
        executable_recipe = fixture_artifact(
            f"{root}/CMakeFiles/{executable_name}.dir/link.txt",
            "build_metadata", character)
        executable_recipe.update({
            "size": executable_text["size"],
            "sha256": executable_text["sha256"]})
        return {
            "build_dir": root,
            "validated_cache": cache,
            "compiler": dict(compiler),
            "compiler_invocation": {
                "invocation": "/usr/bin/compiler",
                "resolved_path": "/usr/bin/compiler",
            },
            "compiler_version_stdout": dict(compiler_version),
            "cmake_cache": fixture_artifact(
                f"{root}/CMakeCache.txt", "build_metadata", character),
            "compile_commands": fixture_artifact(
                f"{root}/compile_commands.json", "build_metadata", character),
            "executable_link_recipe": executable_recipe,
            "archive_link_recipe": archive_recipe,
            "archive_link_recipe_content": archive_text,
            "executable_link_recipe_content": executable_text,
            "archiver": fixture_artifact(
                "/usr/bin/ar", "archiver", "b", 0o755),
            "ranlib": fixture_artifact(
                "/usr/bin/ranlib", "ranlib", "c", 0o755),
            "archive_link_tool_invocations": {
                "archiver": {"invocation": "/usr/bin/ar",
                             "resolved_path": "/usr/bin/ar"},
                "ranlib": {"invocation": "/usr/bin/ranlib",
                           "resolved_path": "/usr/bin/ranlib"},
            },
            "validated_external_link_inputs": [
                fixture_external_record(
                    "/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so",
                    "openmp_runtime_shared", "shared_library"),
                fixture_external_record(
                    "/usr/lib/x86_64-linux-gnu/libpthread.a",
                    "pthread_support_archive", "archive"),
            ],
            "validated_archive": fixture_artifact(
                f"{root}/{archive_name}", "archive", character),
            "validated_executable": fixture_artifact(
                f"{root}/{executable_name}", "executable", character, 0o755),
            "validated_archive_members": members,
            "validated_compile_commands": {
                "schema": COMPILE_COMMANDS_SCHEMA,
                "implementation": role,
                "profile": (BASELINE_COMPILE_PROFILE if baseline else
                            CANDIDATE_COMPILE_PROFILE),
                "entry_count": (
                    BASELINE_EXPECTED_COMPILE_COMMAND_COUNT if baseline else
                    CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT),
                "required_sources": sorted(pair["source"]["path"] for pair in pairs),
                "validated_optimization": "-O3", "validated_openmp": True,
                "required_source_object_pairs": pairs,
                "required_entries": sorted([{
                    "directory": root,
                    "file": pair["source"]["path"],
                    "output": _normalized_compile_output(
                        role, pair["source"]["path"]),
                    "arguments": _normalized_compile_argv(
                        role, pair["source"]["path"], "/usr/bin/compiler"),
                } for pair in pairs], key=lambda entry: entry["file"]),
                "isa_policy": (
                    "whole-build -march=native" if baseline else
                    "portable core with ISA flags only on SSSE3, AVX2, and "
                    "AVX-512VL translation units"),
            },
        }

    def fixture_runtime(executable: str, character: str) -> dict[str, Any]:
        raw_text = (
            "linux-vdso.so.1 (0x0000000000000000)\n"
            "libc.so.6 => /lib/libc.so.6 (0x0000000000000000)\n"
            "libm.so.6 => /lib/libm.so.6 (0x0000000000000000)\n"
            "/lib/ld-linux-x86-64.so.2 (0x0000000000000000)\n")
        return {
            "executable": executable,
            "raw_ldd_output": fixture_text(raw_text),
            "dependencies": [
                {
                    "soname": "ld-linux-x86-64.so.2",
                    "loader_path": "/lib/ld-linux-x86-64.so.2",
                    "file": fixture_artifact(
                        "/lib/ld-linux-x86-64.so.2",
                        "dynamic_loader", character),
                },
                {
                    "soname": "libc.so.6", "loader_path": "/lib/libc.so.6",
                    "file": fixture_artifact(
                        "/lib/libc.so.6",
                        "shared_library", character),
                },
                {
                    "soname": "libm.so.6", "loader_path": "/lib/libm.so.6",
                    "file": fixture_artifact(
                        "/lib/libm.so.6", "shared_library",
                        "a" if character != "a" else "b"),
                },
                {"soname": "linux-vdso.so.1", "virtual": True},
            ],
        }

    def fixture_cpu(cpu: int, sibling: int) -> dict[str, Any]:
        return {
            "cpu": cpu, "online": "1",
            "cpuinfo": {"processor": str(cpu), "model name": "fixture"},
            "topology": {
                "core_id": "2", "physical_package_id": "0", "die_id": "0",
                "cluster_id": None, "thread_siblings_list": "2,6",
                "core_siblings_list": "0-7",
            },
            "frequency_policy": {
                "scaling_driver": "fixture", "scaling_governor": "performance",
                "energy_performance_preference": "performance",
                "scaling_min_freq": "1", "scaling_max_freq": "2",
                "cpuinfo_min_freq": "1", "cpuinfo_max_freq": "2",
            },
            "cache_hierarchy": [
                {
                    "index": 0, "level": "1", "type": "Data", "size": "32K",
                    "coherency_line_size": "64", "number_of_sets": "1",
                    "ways_of_associativity": "1",
                    "physical_line_partition": "1",
                    "shared_cpu_list": "2,6", "shared_cpu_map": "44",
                    "allocation_policy": None, "write_policy": "WriteBack",
                },
                {
                    "index": 1, "level": "3", "type": "Unified", "size": "8M",
                    "coherency_line_size": "64", "number_of_sets": "1",
                    "ways_of_associativity": "1",
                    "physical_line_partition": "1",
                    "shared_cpu_list": "0-7", "shared_cpu_map": "ff",
                    "allocation_policy": None, "write_policy": "WriteBack",
                },
            ],
            "cache_index_inventory": ["index0", "index1"],
            "cache_directory_inventory_text": fixture_text("index0\nindex1\n"),
            "numa_nodes": [0],
            "numa_node_inventory": ["node0"],
            "core_class": {"core_type": None, "cpu_capacity": None},
        }

    baseline_build = fixture_build("baseline", "d")
    candidate_build = fixture_build("candidate", "e")
    baseline_archive = baseline_build["validated_archive"]
    baseline_executable = baseline_build["validated_executable"]
    candidate_archive = candidate_build["validated_archive"]
    candidate_executable = candidate_build["validated_executable"]
    baseline_raw = base64.b64decode(EXACT_MAIN_COMMIT_BASE64, validate=True)
    baseline_head, baseline_commit_object = fixture_commit_object(baseline_raw)
    require(baseline_head == EXACT_MAIN_COMMIT,
            "embedded exact-main commit fixture differs")
    candidate_tree = "c" * 40
    candidate_raw = (
        f"tree {candidate_tree}\nauthor Fixture <fixture@example.com> 1 +0000\n"
        "committer Fixture <fixture@example.com> 1 +0000\n\nfixture\n"
    ).encode("utf-8")
    candidate_head, candidate_commit_object = fixture_commit_object(candidate_raw)
    return {
        "schema": "leopard2-balanced-evidence-scope/v3",
        "host": {
            "system": {
                "system": "Linux", "node": "fixture", "release": "fixture",
                "version": "fixture", "machine": "x86_64", "python": "3.11",
                "page_size": 4096,
            },
            "allowed_cpu_set_at_launch": list(range(8)),
            "online_cpu_set": list(range(8)),
            "online_cpu_list_text": fixture_text("0-7\n"),
            "online_node_list_text": fixture_text("0\n"),
            "benchmark_cpu": fixture_cpu(2, 6),
            "reserved_sibling": fixture_cpu(6, 2),
            "turbo_and_pstate": {
                "/sys/devices/system/cpu/intel_pstate/no_turbo": "0",
                "/sys/devices/system/cpu/amd_pstate/status": None,
                "/sys/devices/system/cpu/cpufreq/boost": None,
            },
        },
        "sources": {
            "baseline": {"path": "$BASELINE_SOURCE", "head": EXACT_MAIN_COMMIT,
                         "tree": EXACT_MAIN_TREE, "detached": True,
                         "tracked_tree_listing_sha256": "b" * 64,
                         "tracked_status": "clean",
                         "commit_object": baseline_commit_object},
            "candidate": {"path": "$CANDIDATE_SOURCE", "head": candidate_head,
                          "tree": candidate_tree, "detached": False,
                          "tracked_tree_listing_sha256": "d" * 64,
                          "tracked_status": "clean",
                          "commit_object": candidate_commit_object},
        },
        "builds": {
            "baseline": baseline_build,
            "candidate": candidate_build,
        },
        "artifacts": {
            "baseline_archive": baseline_archive,
            "baseline_executable": baseline_executable,
            "candidate_archive": candidate_archive,
            "candidate_executable": candidate_executable,
        },
        "runtime_closures": {
            "baseline": fixture_runtime(baseline_executable["path"], "1"),
            "candidate": fixture_runtime(candidate_executable["path"], "2"),
        },
        "tools": {
            "runner": fixture_artifact("/fixture/run.py", "file", "3"),
            "taskset": fixture_artifact(
                "/usr/bin/taskset", "executable", "4", 0o755),
            "ldd": fixture_artifact(
                "/usr/bin/ldd", "executable", "5", 0o755),
        },
        "resolved_auto_backend": backend,
        "forced_confirmation_backends": list(
            BACKENDS[:BACKENDS.index(backend) + 1]),
        "excluded_backends": dict(EXCLUDED_CAMPAIGN_BACKENDS),
    }


def fake_gate_manifests(
    plan_root: Path, backend: str = "avx2"
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    plan = validate_plan(plan_root)
    manifests = []
    references = []
    scopes = []
    for item in plan["artifacts"]:
        scope = fake_evidence_scope(backend)
        gate = validate_main_document(load_json(plan_root / item["path"]), item["path"])
        campaign_cells = []
        analysis = {}
        for cell in gate["cells"]:
            campaign_cells.append({
                "identifier": cell["identifier"], "k": cell["K"], "r": cell["R"],
                "shard_bytes": cell["shard_bytes"], "losses": cell["loss_count"],
                "seed": cell["seed"],
            })
            passes = (cell["K"], cell["shard_bytes"]) in {
                (127, 65536), (128, 4096),
            }
            lower = 1.06 if passes else 0.99
            analysis[cell["identifier"]] = {
                metric: {"ci95_lower": lower,
                         "ratio_orientation": "baseline_time_over_candidate_time"}
                for metric in ("decode_first_use", "decode_reuse_amortized")
            }
        manifests.append({
            "schema": EXACT_MANIFEST_SCHEMA, "valid": True,
            "campaign": {
                "candidate_mode": "generic", "batch": 1, "reuse": 8,
                "iterations": 9, "warmup": 2, "threads": 1, "rounds": 3,
                "order": list(EXACT_ORDER), "cells": campaign_cells,
            },
            "identities": {
                "baseline_source": {"head": EXACT_MAIN_COMMIT},
                "candidate_source": {
                    "head": scope["sources"]["candidate"]["head"]},
            },
            "analysis": analysis,
        })
        references.append({
            "path": "/fixture/" + Path(item["path"]).name,
            "size": 1, "sha256": "2" * 64, "payload_digest": "3" * 64,
        })
        scopes.append(scope)
    return manifests, references, scopes


def fake_attestation_output(case: dict[str, Any], selected_path: str | None = None,
                            selected_rule: str | None = None,
                            backend: str = "avx2") -> dict[str, Any]:
    cell = case["cell"]
    expected_path, expected_rule = expected_auto_decode_pair(case, backend)
    if selected_path is None:
        selected_path = expected_path
    if selected_rule is None:
        selected_rule = expected_rule
    transform = ceil_power_of_two(cell["R"])
    parent = ceil_power_of_two(cell["K"] + transform)
    required_slots = expected_required_work_slots(case, selected_path)
    return {
        "schema": common.BENCHMARK_SCHEMA,
        "build": {
            "compiler": "fixture", "compiler_version": "1", "cplusplus": 201103,
        },
        "parameters": {
            "K": cell["K"], "R": cell["R"],
            "requested_profile": "legacy_high_v1", "requested_field": "gf8",
            "requested_backend": "auto", "force_generic_decode": False,
            "force_specialized_decode": False, "force_tiled_decode": False,
            "force_materialized_decode": False, "skip_legacy": True,
            "retain_samples": True, "report_decode_path": True,
            "shard_bytes": cell["shard_bytes"],
            "loss_count": cell["loss_count"],
            "missing_original_indices": expected_missing_indices(
                cell["K"], cell["loss_count"], cell["seed"]),
            "batch": 1, "reuse": 1, "iterations": 1, "warmup": 0,
            "thread_count": 1, "seed": cell["seed"],
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8", "backend": backend,
            "thread_count": 1, "parent_count": parent,
            "padded_side": transform, "selected_decode_path": selected_path,
            "selected_decode_rule": selected_rule,
            "decode_required_work_slots": required_slots,
            "decode_aligned_prefix_bytes": cell["shard_bytes"] & ~63,
            "decode_tail_bytes": cell["shard_bytes"] & 63,
            "decode_rounded_bytes": (cell["shard_bytes"] + 63) & ~63,
            "decode_multi_item_batch": False,
        },
        "correctness": {
            "leopard2_round_trip": True, "legacy_comparison": None,
        },
        "workload_digests": {
            "algorithm": "fnv1a64", "original_data": "1" * 16,
            "transmitted_parity": "2" * 16, "recovered_originals": "3" * 16,
        },
        "memory": {
            "scratch_alignment": 64, "encode_scratch_bytes_per_stripe": 64,
            "decode_scratch_bytes_per_stripe": 128,
            "encode_scratch_bytes_batch": 64, "decode_scratch_bytes_batch": 128,
        },
        "metrics": {
            "codec_setup": {"samples_us": [1.0]},
            "encode_execution": {"samples_us_per_batch_call": [1.0]},
            "decode_plan_setup": {"samples_us": [1.0]},
            "decode_execution": {"samples_us_per_batch_call": [1.0]},
            "decode_amortized_at_reuse": {"reuse_count": 1},
            "rate_semantics": "fixture",
        },
        "legacy": {
            "available": False, "unavailable_reason": "disabled by --skip-legacy",
            "codec_setup": None, "decode_timing_includes_setup": True,
            "encode_execution": None, "decode_including_setup": None,
        },
    }


def self_test() -> None:
    def retained_text(text: str) -> dict[str, Any]:
        encoded = text.encode("utf-8")
        return {
            "encoding": "utf-8", "text": text, "size": len(encoded),
            "sha256": hashlib.sha256(encoded).hexdigest(),
        }

    exact_runner = load_exact_main_runner()
    require(tuple(exact_runner.BASELINE_LIBRARY_SOURCES) ==
                BASELINE_LIBRARY_SOURCES and
            tuple(exact_runner.CANDIDATE_LIBRARY_SOURCES) ==
                CANDIDATE_LIBRARY_SOURCES and
            exact_runner.BASELINE_EXPECTED_COMPILE_COMMAND_COUNT ==
                BASELINE_EXPECTED_COMPILE_COMMAND_COUNT and
            exact_runner.CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT ==
                CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT,
            "balanced scope translation-unit contract drifted from its producer")
    with tempfile.TemporaryDirectory(prefix="leopard2-balanced-plan-") as temporary:
        root = Path(temporary)
        first = root / "first"
        second = root / "second"
        plan_a = generate_plan(first)
        validate_plan(first)
        plan_b = generate_plan(second)
        validate_plan(second)
        require(canonical_bytes(plan_a) == canonical_bytes(plan_b),
                "two canonical plans differ")

        forged = root / "forged"
        shutil.copytree(first, forged)
        gate_path = forged / "exact-main-gate" / "t128.json"
        adversarial_resign(gate_path, lambda value: value.update(
            {"candidate_mode": "auto"}))
        plan_path = forged / "plan.json"
        def reseal_plan(value):
            target = next(item for item in value["artifacts"]
                          if item["path"] == "exact-main-gate/t128.json")
            target["size"] = gate_path.stat().st_size
            target["sha256"] = file_sha256(gate_path)
        adversarial_resign(plan_path, reseal_plan)
        try:
            validate_plan(forged)
        except PlanError:
            pass
        else:
            raise PlanError("re-signed mode forgery passed canonical validation")

        forged_gate = root / "forged-gate"
        shutil.copytree(first, forged_gate)
        adversarial_resign(forged_gate / "plan.json", lambda value:
                           value["promotion_gate"].update(
                               {"minimum_ci95_lower": 0.01}))
        try:
            validate_plan(forged_gate)
        except PlanError:
            pass
        else:
            raise PlanError("re-signed gate forgery passed canonical validation")

        manifests, references, scopes = fake_gate_manifests(first)
        survivor_signed = derive_survivors(first, manifests, references, scopes)
        require(not survivor_signed["required_refinement_cells"] and
                [(cell["K"], cell["shard_bytes"])
                 for cell in survivor_signed["survivor_cells"]] == [
                    (128, 4096), (127, 65536),
                 ],
                "fixture survivor selection differs")
        require(forced_backends_for_scope(fake_evidence_scope("scalar")) ==
                    ("scalar",) and
                forced_backends_for_scope(fake_evidence_scope("ssse3")) ==
                    ("scalar", "ssse3") and
                forced_backends_for_scope(fake_evidence_scope("avx2")) == BACKENDS,
                "forced backend prefix does not match resolved AUTO tier")
        reordered_scope = json.loads(json.dumps(fake_evidence_scope()))
        baseline_semantics = reordered_scope["builds"]["baseline"][
            "validated_compile_commands"]
        baseline_pairs = baseline_semantics["required_source_object_pairs"]
        benchmark_pair = next(pair for pair in baseline_pairs
                              if pair["source"]["path"].endswith(
                                  "/legacy_main_benchmark.cpp"))
        baseline_pairs.remove(benchmark_pair)
        baseline_pairs.insert(0, benchmark_pair)
        baseline_semantics["required_sources"] = [
            pair["source"]["path"] for pair in baseline_pairs]
        validate_evidence_scope(reordered_scope)

        def reject_scope_mutation(label: str, mutate) -> None:
            forged_scopes = json.loads(json.dumps(scopes))
            mutate(forged_scopes)
            try:
                derive_survivors(first, manifests, references, forged_scopes)
            except PlanError:
                return
            raise PlanError(f"mixed/invalid evidence scope accepted {label}")

        reject_scope_mutation("host", lambda values:
                              values[-1]["host"]["system"].update({
                                  "machine": "aarch64"}))
        def move_to_other_cache_class(values) -> None:
            pair = values[-1]["host"]
            pair["benchmark_cpu"].update({"cpu": 3})
            pair["benchmark_cpu"]["cpuinfo"].update({
                "processor": "3", "cache domain": "32MiB"})
            pair["benchmark_cpu"]["topology"].update({
                "core_id": "3", "thread_siblings_list": "3,7"})
            pair["reserved_sibling"].update({"cpu": 7})
            pair["reserved_sibling"]["cpuinfo"].update({
                "processor": "7", "cache domain": "32MiB"})
            pair["reserved_sibling"]["topology"].update({
                "core_id": "3", "thread_siblings_list": "3,7"})
        reject_scope_mutation("equal-cardinality mixed cache class",
                              move_to_other_cache_class)
        reject_scope_mutation("candidate compiler", lambda values:
                              values[-1]["builds"]["candidate"]["compiler"].update({
                                  "sha256": "d" * 64}))
        reject_scope_mutation("baseline source", lambda values:
                              values[-1]["sources"]["baseline"].update({
                                  "tree": "d" * 40}))
        reject_scope_mutation("candidate CMake", lambda values:
                              values[-1]["builds"]["candidate"][
                                  "validated_cache"].update({
                                      "LEO2_BACKEND_VARIANT": "scalar"}))
        reject_scope_mutation("uniform noncanonical candidate CMake", lambda values:
                              [value["builds"]["candidate"][
                                  "validated_cache"].update({
                                      "LEO2_BACKEND_VARIANT": "scalar"})
                               for value in values])
        for label, flags in (
            ("final optimization downgrade", "-O3 -O0"),
            ("bare optimization downgrade", "-O3 -O"),
            ("unknown optimization spelling", "-O4 -O3"),
            ("sanitizer after optimization", "-O3 -fsanitize=address"),
            ("profile after optimization", "-O3 -fprofile-generate"),
            ("LTO after optimization", "-O3 -flto"),
            ("instrumentation after optimization",
             "-O3 -finstrument-functions"),
            ("vector disable after optimization",
             "-O3 -fno-tree-vectorize"),
            ("coverage after optimization", "-O3 --coverage"),
            ("long optimize alias", "-O3 --optimize=0"),
            ("long sanitizer alias", "-O3 --sanitize=address"),
            ("long instrumentation alias", "-O3 --instrument-functions"),
            ("long LTO alias", "-O3 --lto"),
            ("long profile alias", "-O3 --profile"),
            ("inline disable", "-O3 -fno-inline"),
            ("GCC pass disable", "-O3 -fdisable-tree-vect"),
            ("Release response file", "-O3 @evil.rsp"),
        ):
            reject_scope_mutation(
                f"uniform CMake {label}",
                lambda values, item=flags: [
                    value["builds"][role]["validated_cache"].update({
                        "CMAKE_CXX_FLAGS_RELEASE": item})
                    for value in values
                    for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform mislabeled candidate compile profile",
            lambda values: [
                value["builds"]["candidate"][
                    "validated_compile_commands"].update({
                        "profile": BASELINE_COMPILE_PROFILE})
                for value in values])
        reject_scope_mutation(
            "uniform mislabeled baseline compile implementation",
            lambda values: [
                value["builds"]["baseline"][
                    "validated_compile_commands"].update({
                        "implementation": "candidate"})
                for value in values])
        reject_scope_mutation(
            "uniform obsolete compile-command schema",
            lambda values: [
                value["builds"][role][
                    "validated_compile_commands"].update({
                        "schema": "leopard2-main-compare-compile-commands/v1"})
                for value in values for role in ("baseline", "candidate")])

        def mutate_all_compiler_argv(
            values, role: str, suffix: str, mutation: str,
        ) -> None:
            for value in values:
                semantics = value["builds"][role]["validated_compile_commands"]
                entry = next(item for item in semantics["required_entries"]
                             if item["file"].endswith(suffix))
                arguments = entry["arguments"]
                source_index = arguments.index(entry["file"])
                output_index = arguments.index("-o")
                compile_index = arguments.index("-c")
                other = next(item for item in semantics["required_entries"]
                             if item["file"] != entry["file"])
                if mutation == "response":
                    arguments.insert(source_index, "@evil.rsp")
                elif mutation == "extra_source":
                    arguments.insert(source_index, other["file"])
                elif mutation == "stdin_source":
                    arguments.insert(source_index, "-")
                elif mutation == "wrong_source":
                    arguments[source_index] = other["file"]
                elif mutation == "missing_source":
                    arguments.pop(source_index)
                elif mutation == "duplicate_source":
                    arguments.append(entry["file"])
                elif mutation == "wrong_output":
                    arguments[output_index + 1] = other["output"]
                elif mutation == "coherent_output":
                    entry["output"] = other["output"]
                    arguments[output_index + 1] = other["output"]
                elif mutation == "duplicate_output":
                    arguments[compile_index:compile_index] = [
                        "-o", entry["output"]]
                elif mutation == "missing_compile":
                    arguments.pop(compile_index)
                elif mutation == "duplicate_compile":
                    arguments.insert(compile_index, "-c")
                elif mutation == "extra_definition":
                    arguments.insert(1, "-DEVIL=1")
                elif mutation == "reorder_definitions":
                    definitions = [
                        index for index, token in enumerate(arguments)
                        if token.startswith("-D")]
                    arguments[definitions[0]], arguments[definitions[1]] = \
                        arguments[definitions[1]], arguments[definitions[0]]
                elif mutation == "extra_include":
                    include = next(
                        index for index, token in enumerate(arguments)
                        if token.startswith("-I"))
                    arguments.insert(include + 1, "-I/tmp/evil")
                elif mutation == "different_include":
                    include = next(
                        index for index, token in enumerate(arguments)
                        if token.startswith("-I"))
                    arguments[include] = "-I/tmp/evil"
                elif mutation == "language_last":
                    arguments.insert(
                        arguments.index("-std=gnu++11") + 1,
                        "-std=gnu++20")
                elif mutation == "optimization_last":
                    arguments.insert(output_index, "-O0")
                elif mutation == "avx2_last":
                    arguments.insert(output_index, "-mno-avx2")
                elif mutation == "missing_avx2":
                    arguments.remove("-mavx2")
                elif mutation == "baseline_isa_last":
                    arguments.insert(output_index, "-march=x86-64")
                elif mutation == "compiler_wrapper":
                    arguments.insert(1, "--compiler-wrapper=/tmp/evil")
                else:
                    raise PlanError(
                        f"unknown compiler self-test mutation: {mutation}")

        for label, role, suffix, mutation in (
            ("compiler response file", "candidate", "/leopard.cpp", "response"),
            ("compiler second positional source", "candidate", "/leopard.cpp", "extra_source"),
            ("compiler stdin positional source", "candidate", "/leopard.cpp", "stdin_source"),
            ("compiler different positional source", "candidate", "/leopard.cpp", "wrong_source"),
            ("compiler missing positional source", "candidate", "/leopard.cpp", "missing_source"),
            ("compiler duplicate positional source", "candidate", "/leopard.cpp", "duplicate_source"),
            ("compiler wrong output", "candidate", "/leopard.cpp", "wrong_output"),
            ("compiler coherent wrong output", "candidate", "/leopard.cpp", "coherent_output"),
            ("compiler duplicate output", "candidate", "/leopard.cpp", "duplicate_output"),
            ("compiler missing compile option", "candidate", "/leopard.cpp", "missing_compile"),
            ("compiler duplicate compile option", "candidate", "/leopard.cpp", "duplicate_compile"),
            ("compiler extra definition", "candidate", "/leopard.cpp", "extra_definition"),
            ("compiler reordered definitions", "candidate", "/leopard.cpp", "reorder_definitions"),
            ("compiler extra include", "candidate", "/leopard.cpp", "extra_include"),
            ("compiler different include", "candidate", "/leopard.cpp", "different_include"),
            ("compiler last language option", "candidate", "/leopard.cpp", "language_last"),
            ("compiler last optimization option", "candidate", "/leopard.cpp", "optimization_last"),
            ("compiler last AVX2 option", "candidate", "/Leopard2BackendAVX2.cpp", "avx2_last"),
            ("compiler missing AVX2 option", "candidate", "/Leopard2BackendAVX2.cpp", "missing_avx2"),
            ("compiler last baseline ISA option", "baseline", "/leopard.cpp", "baseline_isa_last"),
            ("compiler wrapper control", "candidate", "/leopard.cpp", "compiler_wrapper"),
        ):
            reject_scope_mutation(
                "uniform " + label,
                lambda values, item_role=role, item_suffix=suffix,
                       item_mutation=mutation:
                    mutate_all_compiler_argv(
                        values, item_role, item_suffix, item_mutation))
        reject_scope_mutation("uniform incomplete sources", lambda values:
                              [value["sources"][role].pop("tree")
                               for value in values
                               for role in ("baseline", "candidate")])
        reject_scope_mutation("uniform incomplete runtime files", lambda values:
                              [dependency.update({
                                   "file": {"path": dependency["file"]["path"]}})
                               for value in values
                               for role in ("baseline", "candidate")
                               for dependency in value["runtime_closures"][role][
                                   "dependencies"]
                               if "file" in dependency])
        reject_scope_mutation("uniform path-only compile pairs", lambda values:
                              [pair.update({
                                   "source": {"path": pair["source"]["path"]},
                                   "object": {"path": pair["object"]["path"]}})
                               for value in values
                               for role in ("baseline", "candidate")
                               for pair in value["builds"][role][
                                   "validated_compile_commands"][
                                       "required_source_object_pairs"]])
        reject_scope_mutation("uniform empty compiler/CMake/link records",
                              lambda values:
                              [value["builds"][role].update({
                                   name: {} for name in (
                                       "compiler", "cmake_cache", "compile_commands",
                                       "executable_link_recipe",
                                       "archive_link_recipe")})
                               for value in values
                               for role in ("baseline", "candidate")])
        def reduce_all_outputs(values) -> None:
            for value in values:
                for role in ("baseline", "candidate"):
                    for output in ("archive", "executable"):
                        key = f"{role}_{output}"
                        retained = value["artifacts"][key]
                        reduced = {name: retained[name]
                                   for name in ("path", "kind", "sha256")}
                        value["artifacts"][key] = dict(reduced)
                        value["builds"][role][f"validated_{output}"] = dict(reduced)
        reject_scope_mutation("uniform reduced archive/executable",
                              reduce_all_outputs)
        def reduce_all_hosts_to_topology(values) -> None:
            for value in values:
                host = value["host"]
                for name in ("benchmark_cpu", "reserved_sibling"):
                    retained = host[name]
                    host[name] = {"cpu": retained["cpu"],
                                  "topology": retained["topology"]}
        reject_scope_mutation("uniform topology-only host",
                              reduce_all_hosts_to_topology)
        def truncate_all_translation_unit_closures(values) -> None:
            for value in values:
                for role, names in (
                    ("baseline", BASELINE_LIBRARY_SOURCES),
                    ("candidate", CANDIDATE_LIBRARY_SOURCES),
                ):
                    build = value["builds"][role]
                    semantics = build["validated_compile_commands"]
                    pair = next(item for item in
                                semantics["required_source_object_pairs"]
                                if item["source"]["path"].endswith(
                                    "/" + names[-1]))
                    semantics["required_source_object_pairs"].remove(pair)
                    semantics["required_sources"].remove(pair["source"]["path"])
                    semantics["entry_count"] -= 1
                    member = Path(pair["object"]["path"]).name
                    build["validated_archive_members"].remove(member)
                    relative = pair["object"]["path"][
                        len(build["build_dir"]) + 1:]
                    text = build["archive_link_recipe_content"]["text"].replace(
                        " " + relative, "", 1)
                    retained = {
                        "encoding": "utf-8", "text": text,
                        "size": len(text.encode("utf-8")),
                        "sha256": hashlib.sha256(
                            text.encode("utf-8")).hexdigest(),
                    }
                    build["archive_link_recipe_content"] = retained
                    build["archive_link_recipe"]["size"] = retained["size"]
                    build["archive_link_recipe"]["sha256"] = retained["sha256"]
        reject_scope_mutation(
            "uniform coherent translation-unit truncation",
            truncate_all_translation_unit_closures)
        def redirect_all_declared_archives(values) -> None:
            for value in values:
                for role, archive_name in (
                    ("baseline", "libleopard_main_exact.a"),
                    ("candidate", "libleopard.a"),
                ):
                    build = value["builds"][role]
                    text = build[
                        "executable_link_recipe_content"]["text"].replace(
                            f" {archive_name} ",
                            f" /tmp/{archive_name} ", 1)
                    retained = retained_text(text)
                    build["executable_link_recipe_content"] = retained
                    build["executable_link_recipe"]["size"] = retained["size"]
                    build["executable_link_recipe"]["sha256"] = retained["sha256"]
        reject_scope_mutation(
            "uniform same-basename foreign archive operands",
            redirect_all_declared_archives)
        def add_all_foreign_archive_operands(values) -> None:
            for value in values:
                for role, archive_name in (
                    ("baseline", "libleopard_main_exact.a"),
                    ("candidate", "libleopard.a"),
                ):
                    build = value["builds"][role]
                    text = build[
                        "executable_link_recipe_content"]["text"].replace(
                            f" {archive_name} ",
                            f" {archive_name} /tmp/evil.a ", 1)
                    retained = retained_text(text)
                    build["executable_link_recipe_content"] = retained
                    build["executable_link_recipe"]["size"] = retained["size"]
                    build["executable_link_recipe"]["sha256"] = retained["sha256"]
        reject_scope_mutation(
            "uniform foreign archive operands",
            add_all_foreign_archive_operands)
        def add_all_link_controls(values, control: str) -> None:
            for value in values:
                for role, archive_name in (
                    ("baseline", "libleopard_main_exact.a"),
                    ("candidate", "libleopard.a"),
                ):
                    build = value["builds"][role]
                    text = build[
                        "executable_link_recipe_content"]["text"].replace(
                            f" {archive_name} ",
                            f" {archive_name} {control} ", 1)
                    retained = retained_text(text)
                    build["executable_link_recipe_content"] = retained
                    build["executable_link_recipe"]["size"] = retained["size"]
                    build["executable_link_recipe"]["sha256"] = retained["sha256"]
        for label, control in (
            ("fused linker script", "-Wl,--script=/tmp/evil.ld"),
            ("fused search/library", "-Wl,-L,/tmp,-levil"),
            ("driver response", "@evil.rsp"),
            ("compiler specs", "-specs=/tmp/evil.specs"),
            ("alternate tool root", "-B/tmp/toolchain"),
            ("compiler plugin", "-fplugin=/tmp/evil.so"),
            ("linker plugin", "-Wl,--plugin,/tmp/evil.so"),
            ("alternate linker", "-fuse-ld=/tmp/evil-ld"),
            ("alternate runtime loader",
             "-Wl,--dynamic-linker,/tmp/evil-ld.so"),
        ):
            reject_scope_mutation(
                f"uniform {label}",
                lambda values, item=control:
                    add_all_link_controls(values, item))
        def redirect_all_pthread_inputs(values) -> None:
            old = "/usr/lib/x86_64-linux-gnu/libpthread.a"
            new = "/tmp/libpthread.a"
            for value in values:
                for role in ("baseline", "candidate"):
                    build = value["builds"][role]
                    record = build["validated_external_link_inputs"][1]
                    record["operand"] = new
                    record["artifact"]["path"] = new
                    text = build[
                        "executable_link_recipe_content"]["text"].replace(
                            old, new, 1)
                    retained = retained_text(text)
                    build["executable_link_recipe_content"] = retained
                    build["executable_link_recipe"]["size"] = retained["size"]
                    build["executable_link_recipe"]["sha256"] = retained["sha256"]
        reject_scope_mutation(
            "uniform alternate pthread root", redirect_all_pthread_inputs)
        reject_scope_mutation(
            "uniform missing external static identity",
            lambda values: [
                value["builds"][role]["validated_external_link_inputs"].pop()
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "candidate external identity differs from baseline",
            lambda values: [
                value["builds"]["candidate"][
                    "validated_external_link_inputs"][1]["artifact"].update({
                        "sha256": "e" * 64})
                for value in values])
        reject_scope_mutation(
            "uniform external resolved path escape",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][1]["artifact"].update({
                        "path": "/tmp/libpthread.a"})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform nonexistent versioned OpenMP target",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0]["artifact"].update({
                        "path": "/usr/lib/x86_64-linux-gnu/"
                                "libgomp.so.999.999"})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform differing external roles",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][1].update({
                        "role": "openmp_runtime_shared"})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform external file mode drift",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0]["artifact"].update({
                        "mode": value["builds"][role][
                            "validated_external_link_inputs"][0]["artifact"][
                                "mode"] ^ 0o100})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform external lexical-link target drift",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0][
                        "lexical_symlink_chain"][0].update({
                            "target": "../../../x86_64-linux-gnu/"
                                      "libgomp.so.999"})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform external lexical-link mode drift",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0][
                        "lexical_symlink_chain"][0].update({"mode": 0o700})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform truncated external lexical-link chain",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][0][
                        "lexical_symlink_chain"].pop()
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation(
            "uniform fabricated pthread lexical-link chain",
            lambda values: [
                value["builds"][role][
                    "validated_external_link_inputs"][1].update({
                        "lexical_symlink_chain": [{
                            "path": "/usr/lib/x86_64-linux-gnu/libpthread.a",
                            "target": "libpthread-evil.a", "mode": 0o777,
                        }]})
                for value in values for role in ("baseline", "candidate")])
        for label, control in (
            ("recipe optimization downgrade", "-O0"),
            ("recipe sanitizer", "-fsanitize=address"),
        ):
            reject_scope_mutation(
                f"uniform {label}",
                lambda values, item=control:
                    add_all_link_controls(values, item))
        def remove_all_dynamic_loader_records(values) -> None:
            for value in values:
                for role in ("baseline", "candidate"):
                    closure = value["runtime_closures"][role]
                    closure["dependencies"] = [
                        item for item in closure["dependencies"]
                        if item["soname"] != "ld-linux-x86-64.so.2"]
                    text = "\n".join(
                        line for line in
                        closure["raw_ldd_output"]["text"].splitlines()
                        if not line.lstrip().startswith(
                            "/lib/ld-linux-x86-64.so.2")) + "\n"
                    closure["raw_ldd_output"] = {
                        "encoding": "utf-8", "text": text,
                        "size": len(text.encode("utf-8")),
                        "sha256": hashlib.sha256(
                            text.encode("utf-8")).hexdigest(),
                    }
        reject_scope_mutation(
            "uniform missing dynamic-loader records",
            remove_all_dynamic_loader_records)
        coherently_rewritten_runtime_scopes = json.loads(json.dumps(scopes))
        for value in coherently_rewritten_runtime_scopes:
            for role in ("baseline", "candidate"):
                closure = value["runtime_closures"][role]
                closure["dependencies"] = [
                    dependency for dependency in closure["dependencies"]
                    if dependency["soname"] != "libm.so.6"]
                text = "\n".join(
                    line for line in closure["raw_ldd_output"]["text"].splitlines()
                    if not line.startswith("libm.so.6 =>")) + "\n"
                closure["raw_ldd_output"] = retained_text(text)
        derive_survivors(
            first, manifests, references, coherently_rewritten_runtime_scopes)
        def swap_all_runtime_file_records(values) -> None:
            for value in values:
                for role in ("baseline", "candidate"):
                    dependencies = value["runtime_closures"][role]["dependencies"]
                    libc = next(item for item in dependencies
                                if item["soname"] == "libc.so.6")
                    libm = next(item for item in dependencies
                                if item["soname"] == "libm.so.6")
                    libc["file"], libm["file"] = libm["file"], libc["file"]
        reject_scope_mutation(
            "uniform swapped runtime file records",
            swap_all_runtime_file_records)
        def make_all_runtime_loader_paths_noncanonical(values) -> None:
            for value in values:
                for role in ("baseline", "candidate"):
                    closure = value["runtime_closures"][role]
                    libc = next(item for item in closure["dependencies"]
                                if item["soname"] == "libc.so.6")
                    libc["loader_path"] = "/lib/./libc.so.6"
                    libc["file"]["path"] = libc["loader_path"]
                    text = closure["raw_ldd_output"]["text"].replace(
                        "/lib/libc.so.6", libc["loader_path"])
                    encoded = text.encode("utf-8")
                    closure["raw_ldd_output"] = {
                        "encoding": "utf-8", "text": text,
                        "size": len(encoded),
                        "sha256": hashlib.sha256(encoded).hexdigest(),
                    }
        reject_scope_mutation(
            "uniform noncanonical runtime loader paths",
            make_all_runtime_loader_paths_noncanonical)
        def truncate_all_cache_records(values) -> None:
            for value in values:
                for name in ("benchmark_cpu", "reserved_sibling"):
                    value["host"][name]["cache_hierarchy"].pop()
                    value["host"][name]["cache_index_inventory"].pop()
        reject_scope_mutation(
            "uniform cache summary/listing mismatch", truncate_all_cache_records)
        coherently_rewritten_cache_scopes = json.loads(json.dumps(scopes))
        for value in coherently_rewritten_cache_scopes:
            for name in ("benchmark_cpu", "reserved_sibling"):
                record = value["host"][name]
                record["cache_hierarchy"].pop()
                record["cache_index_inventory"].pop()
                record["cache_directory_inventory_text"] = retained_text(
                    "".join(f"{entry}\n"
                            for entry in record["cache_index_inventory"]))
        derive_survivors(
            first, manifests, references, coherently_rewritten_cache_scopes)
        def empty_all_numa_records(values) -> None:
            for value in values:
                for name in ("benchmark_cpu", "reserved_sibling"):
                    value["host"][name]["numa_nodes"] = []
                    value["host"][name]["numa_node_inventory"] = []
        reject_scope_mutation(
            "uniform empty NUMA records", empty_all_numa_records)
        reject_scope_mutation(
            "uniform source tree/commit mismatch",
            lambda values: [
                value["sources"][role].update({"tree": "f" * 40})
                for value in values for role in ("baseline", "candidate")])
        reject_scope_mutation("baseline archive", lambda values:
                              values[-1]["artifacts"]["baseline_archive"].update({
                                  "sha256": "6" * 64}))
        reject_scope_mutation("baseline executable", lambda values:
                              values[-1]["artifacts"]["baseline_executable"].update({
                                  "sha256": "7" * 64}))
        reject_scope_mutation("baseline CMake cache", lambda values:
                              values[-1]["builds"]["baseline"][
                                  "cmake_cache"].update({"sha256": "8" * 64}))
        reject_scope_mutation("baseline archive recipe", lambda values:
                              values[-1]["builds"]["baseline"][
                                  "archive_link_recipe_content"].update({
                                      "sha256": "9" * 64}))
        reject_scope_mutation("baseline compile closure", lambda values:
                              values[-1]["builds"]["baseline"][
                                  "validated_compile_commands"][
                                      "required_source_object_pairs"][0][
                                          "object"].update({"sha256": "0" * 64}))
        reject_scope_mutation("baseline runtime dependency", lambda values:
                              values[-1]["runtime_closures"]["baseline"][
                                  "dependencies"][0]["file"].update({
                                      "sha256": "6" * 64}))
        reject_scope_mutation("candidate runtime dependency", lambda values:
                              values[-1]["runtime_closures"]["candidate"][
                                  "dependencies"][0]["file"].update({
                                      "sha256": "7" * 64}))
        reject_scope_mutation("uniform noncanonical baseline CMake", lambda values:
                              [value["builds"]["baseline"][
                                  "validated_cache"].update({
                                      "LEO_MAIN_HAS_MARCH_NATIVE": "0"})
                               for value in values])
        reject_scope_mutation("resolved backend", lambda values:
                              values[-1].update({
                                  "resolved_auto_backend": "ssse3"}))
        for excluded_backend in ("avx512", "neon"):
            reject_scope_mutation(
                "mixed " + excluded_backend,
                lambda values, backend=excluded_backend:
                    values[-1].update({"resolved_auto_backend": backend}))
            reject_scope_mutation(excluded_backend, lambda values, backend=excluded_backend:
                                  [value.update({
                                      "resolved_auto_backend": backend})
                                   for value in values])

        normalized_nested = _normalize_bound_paths(
            {"paths": ["/tree/build/object.o", "/tree/source.cpp"]},
            (("/tree", "$SOURCE"), ("/tree/build", "$BUILD")))
        require(normalized_nested == {
            "paths": ["$BUILD/object.o", "$SOURCE/source.cpp"]},
            "scope path normalization did not prefer the longest role root")
        normalized_collisions = _normalize_bound_paths({
            "exact": "/tree",
            "child": "/tree/child",
            "siblings": [
                "/tree-build/object.o", "/treehouse/object.o",
                "/tree.build/object.o", "/tree!object.o", "/tree#object.o",
                "/tree?object.o", "/other/tree/object.o",
            ],
            "text": ("'/tree/file' /tree-build/file /treehouse/file "
                     "/tree.build/file /tree!file /tree#file /tree?file "
                     "/other/tree/file root=/tree (/tree),/tree:/next"),
        }, (("/tree", "$ROOT"),))
        require(normalized_collisions == {
            "exact": "$ROOT",
            "child": "$ROOT/child",
            "siblings": [
                "/tree-build/object.o", "/treehouse/object.o",
                "/tree.build/object.o", "/tree!object.o", "/tree#object.o",
                "/tree?object.o", "/other/tree/object.o",
            ],
            "text": ("'$ROOT/file' /tree-build/file /treehouse/file "
                     "/tree.build/file /tree!file /tree#file /tree?file "
                     "/other/tree/file root=$ROOT ($ROOT),$ROOT:/next"),
        }, "scope path normalization rewrote a textual sibling prefix")

        def normalize_cmake_external_object(
            source_root: str, build_root: str,
        ) -> dict[str, Any]:
            object_path = (
                f"{build_root}/CMakeFiles/leopard.dir"
                f"{source_root}/LeopardFF8.cpp.o")
            recipe_text = (
                "ar qc libleopard.a CMakeFiles/leopard.dir"
                f"{source_root}/LeopardFF8.cpp.o\n")
            return _normalize_bound_paths({
                "object_path": object_path,
                "archive_link_recipe_content": retained_text(recipe_text),
            }, ((source_root, "$SOURCE"), (build_root, "$BUILD")))

        normalized_home_source = normalize_cmake_external_object(
            "/home/a", "/opt/b")
        normalized_opt_source = normalize_cmake_external_object(
            "/opt/b", "/home/a")
        expected_external_recipe = retained_text(
            "ar qc libleopard.a CMakeFiles/leopard.dir/"
            "$SOURCE/LeopardFF8.cpp.o\n")
        expected_external_object = {
            "archive_link_recipe_content": expected_external_recipe,
            "object_path": (
                "$BUILD/CMakeFiles/leopard.dir/"
                "$SOURCE/LeopardFF8.cpp.o"),
        }
        require(normalized_home_source == expected_external_object and
                normalized_opt_source == expected_external_object,
                "CMake external-source object paths are not root-independent")
        normalized_external_sibling = _normalize_bound_paths({
            "object_path": (
                "/opt/b/CMakeFiles/leopard.dir/"
                "home/a-copy/LeopardFF8.cpp.o"),
        }, (("/home/a", "$SOURCE"), ("/opt/b", "$BUILD")))
        require(normalized_external_sibling == {
            "object_path": (
                "$BUILD/CMakeFiles/leopard.dir/"
                "home/a-copy/LeopardFF8.cpp.o"),
        }, "CMake external-source handling rewrote a sibling path prefix")
        recipe_text = "/tree/build/ar qc /tree/source.cpp\n"
        normalized_recipe = _normalize_bound_paths({
            "archive_link_recipe": {
                "path": "/tree/build/link.txt", "kind": "build_metadata",
                "size": len(recipe_text.encode("utf-8")), "mode": 0o644,
                "sha256": hashlib.sha256(recipe_text.encode("utf-8")).hexdigest(),
                "mtime_ns": 1,
            },
            "archive_link_recipe_content": {
                "encoding": "utf-8", "text": recipe_text,
                "size": len(recipe_text.encode("utf-8")),
                "sha256": hashlib.sha256(recipe_text.encode("utf-8")).hexdigest(),
            },
        }, (("/tree", "$SOURCE"), ("/tree/build", "$BUILD")))
        normalized_content = normalized_recipe["archive_link_recipe_content"]
        require(normalized_content["text"] ==
                    "$BUILD/ar qc $SOURCE/source.cpp\n" and
                normalized_content["sha256"] == hashlib.sha256(
                    normalized_content["text"].encode("utf-8")).hexdigest() and
                normalized_recipe["archive_link_recipe"]["sha256"] ==
                    normalized_content["sha256"] and
                normalized_recipe["archive_link_recipe"]["size"] ==
                    normalized_content["size"],
                "normalized recipe content retained a stale byte identity")

        unsolicited_manifests = json.loads(json.dumps(manifests))
        unsolicited = main_cell(
            "gate-k80-b65536", 80, 80, 65536, 80,
            "exact-main-generic-gate")
        unsolicited_manifests[-1]["campaign"]["cells"].append({
            "identifier": unsolicited["identifier"], "k": unsolicited["K"],
            "r": unsolicited["R"], "shard_bytes": unsolicited["shard_bytes"],
            "losses": unsolicited["loss_count"], "seed": unsolicited["seed"],
        })
        unsolicited_manifests[-1]["analysis"][unsolicited["identifier"]] = {
            metric: {
                "ci95_lower": 1.06,
                "ratio_orientation": "baseline_time_over_candidate_time",
            }
            for metric in ("decode_first_use", "decode_reuse_amortized")
        }
        try:
            derive_survivors(first, unsolicited_manifests, references, scopes)
        except PlanError:
            pass
        else:
            raise PlanError("an unrequested supplemental cell entered selection")

        refinement_manifests = json.loads(json.dumps(manifests))
        for manifest in refinement_manifests:
            for cell in manifest["campaign"]["cells"]:
                if cell["k"] == 112 and cell["shard_bytes"] == 65536:
                    for metric in ("decode_first_use", "decode_reuse_amortized"):
                        manifest["analysis"][cell["identifier"]][metric][
                            "ci95_lower"] = 1.06
        refinement = derive_survivors(
            first, refinement_manifests, references, scopes)
        expected_refinement_k = list(range(97, 112)) + list(range(113, 120))
        require([cell["K"] for cell in
                 refinement["required_refinement_cells"]] == expected_refinement_k,
                "transition refinement cells differ")
        try:
            materialize_stage(first, refinement, root / "premature-stage")
        except PlanError:
            pass
        else:
            raise PlanError("stage advanced before K transition refinement")

        survivor_path = root / "survivors.json"
        write_json(survivor_path, survivor_signed)
        validate_survivors(first, survivor_path, manifests, scopes)

        forged_survivor = root / "forged-survivor.json"
        shutil.copy2(survivor_path, forged_survivor)
        def add_rejected(value):
            rejected = value["rejected_cells"][0]
            value["survivor_cells"].append(rejected)
        adversarial_resign(forged_survivor, add_rejected)
        try:
            validate_survivors(first, forged_survivor, manifests, scopes)
        except PlanError:
            pass
        else:
            raise PlanError("re-signed non-survivor entered the survivor set")

        stage_root = root / "stage"
        materialize_stage(first, survivor_signed, stage_root)
        stage = validate_stage(first, stage_root, manifests, scopes)
        require(stage["survivor_K"] == [127, 128] and
                stage["promotion_requires_path_attestation"] is True and
                stage["expected_resolved_backend"] == "avx2" and
                stage["forced_confirmation_backends"] == list(BACKENDS) and
                set(stage["excluded_backends"]) == {"avx512", "neon"},
                "stage survivor/path-attestation semantics differ")
        survivor_pairs = {(128, 4096), (127, 65536)}
        for path in (stage_root / "forced-survivors").glob("*.json"):
            _, cases = common.normalize_matrix(load_json(path))
            require(all((case["K"], case["shard_bytes"]) in survivor_pairs
                        for case in cases) and
                    {(case["K"], case["shard_bytes"]) for case in cases}
                        .issubset(survivor_pairs),
                    "a non-surviving K/byte entered forced comparison")
        for path in (stage_root / "forced-true-tails").glob("*.json"):
            _, cases = common.normalize_matrix(load_json(path))
            require(all(case["K"] in {127, 128} and
                        case["shard_bytes"] in TRUE_TAIL_BYTES and
                        case["shard_bytes"] % 64 != 0 for case in cases),
                    "aligned or non-surviving cell entered true-tail stage")
        rejection = validate_main_document(
            load_json(stage_root / "exact-main-auto-rejection-timing.json"),
            "rejection timing")
        attestation = unsigned(load_json(stage_root / "path-attestation.json"),
                               ATTESTATION_SCHEMA, "path attestation")
        positive_pairs = {
            (item["cell"]["K"], item["cell"]["shard_bytes"])
            for item in attestation["cases"]
            if item["expected_selected_decode_path"] == "generic"
        }
        require(positive_pairs == survivor_pairs and all(
            item["expected_selected_decode_rule"] == "balanced_generic"
            for item in attestation["cases"]
            if item["expected_selected_decode_path"] == "generic"),
            "attestation broadened the exact surviving gate cells")
        negative = [item for item in attestation["cases"]
                    if item["expected_selected_decode_path"] == "not_generic"]
        require(any(item["cell"]["category"] == "rejected_gate" and
                    item["cell"]["K"] == 128 and
                    item["cell"]["shard_bytes"] == 65536 and
                    item["cell"]["loss_count"] == 128 for item in negative) and
                any(item["cell"]["category"] == "rejected_gate" and
                    item["cell"]["K"] == 127 and
                    item["cell"]["shard_bytes"] == 4096 and
                    item["cell"]["loss_count"] == 127 for item in negative) and
                all(item["expected_selected_decode_rule"] ==
                    "not_balanced_generic" for item in negative),
                "rejected mixed-K full-loss gates were not fail-closed")
        require({tuple((item["cell"][key] for key in (
                    "K", "R", "shard_bytes", "loss_count", "seed")))
                 for item in attestation["cases"]
                 if item["cell"]["category"] == "loss_rate_neighbor"} ==
                {tuple((cell[key] for key in (
                    "K", "R", "shard_bytes", "loss_count", "seed")))
                 for cell in rejection["cells"]},
                "attestation loss/rate neighbors differ from rejection timing")
        require(all(item["expected_selected_decode_path"] == "not_generic"
                    for item in attestation["cases"]
                    if item["cell"]["category"] in {
                        "extra_aligned_confirmation", "true_tail",
                    }), "unproven aligned/tail cells entered AUTO")

        path_names = {pair[0] for pair in PRODUCTION_DECODE_PATH_RULE_PAIRS}
        rule_names = {pair[1] for pair in PRODUCTION_DECODE_PATH_RULE_PAIRS}
        validate_selector_pair_contract(Path(__file__).resolve().parents[3])
        require(all(coherent_production_decode_pair(*pair)
                    for pair in PRODUCTION_DECODE_PATH_RULE_PAIRS),
                "a declared production path/rule pair is not coherent")
        for path_name in path_names | {"banana"}:
            for rule_name in rule_names | {"banana"}:
                require(coherent_production_decode_pair(path_name, rule_name) ==
                        ((path_name, rule_name) in
                         PRODUCTION_DECODE_PATH_RULE_PAIRS),
                        "unknown or cross-paired selector names were accepted")
        for invalid_pair in (
            ("banana", "banana"), ("materialized", "direct"),
            ("direct", "workspace_tiled"),
        ):
            require(not coherent_production_decode_pair(*invalid_pair),
                    f"explicit invalid selector pair was accepted: {invalid_pair!r}")

        materialized_case = next(
            item for item in attestation["cases"]
            if expected_auto_decode_pair(item, "scalar") ==
                ("materialized", "workspace_materialized") and
               item["expected_selected_decode_path"] == "not_generic")
        expected_materialized = expected_auto_decode_pair(
            materialized_case, "scalar")
        for pair in sorted(PRODUCTION_DECODE_PATH_RULE_PAIRS):
            candidate_output = fake_attestation_output(
                materialized_case, pair[0], pair[1], "scalar")
            accepted = True
            try:
                validate_attestation_output(candidate_output, materialized_case)
            except PlanError:
                accepted = False
            require(accepted == (pair == expected_materialized),
                    f"negative cell pair constraint differs for {pair!r}")
        for pair in (
            ("banana", "banana"), ("materialized", "direct"),
            ("direct", "workspace_tiled"),
        ):
            try:
                validate_attestation_output(fake_attestation_output(
                    materialized_case, pair[0], pair[1], "scalar"),
                    materialized_case)
            except PlanError:
                pass
            else:
                raise PlanError(f"negative cell accepted invalid pair {pair!r}")

        direct_case = attestation_case(
            "fixture-direct", "loss_rate_neighbor", 8, 8, 4096, 4,
            deterministic_seed("fixture-direct", 8, 4), False)
        tiled_case = attestation_case(
            "fixture-tiled", "loss_rate_neighbor", 127, 63, 4096, 63,
            deterministic_seed("fixture-tiled", 127, 63), False)
        require(expected_auto_decode_pair(direct_case, "scalar") ==
                ("direct", "direct") and
                expected_auto_decode_pair(tiled_case, "scalar") ==
                ("tiled", "workspace_tiled"),
                "terminal/workspace fixture predictions differ")
        low_workspace_case = {"cell": {
            "K": 33, "R": 7, "profile": "low_v1", "loss_count": 7,
        }}
        require(expected_required_work_slots(direct_case, "no_op") == 0 and
                expected_required_work_slots(direct_case, "direct") == 0 and
                expected_required_work_slots(tiled_case, "tiled") ==
                    2 * ceil_power_of_two(63) + 63 and
                expected_required_work_slots(tiled_case, "materialized") == 256 and
                expected_required_work_slots(low_workspace_case, "tiled") == 128 and
                expected_required_work_slots(low_workspace_case, "generic") == 128,
                "production workspace formulas differ")
        validate_attestation_output(fake_attestation_output(direct_case), direct_case)
        validate_attestation_output(fake_attestation_output(tiled_case), tiled_case)

        source = {
            "head": survivor_signed["candidate_commit"],
            "tree": "4" * 40, "status": "clean",
            "status_sha256": common.EMPTY_SHA256,
        }
        identity = {
            "relative_path": "fixture", "sha256": "5" * 64,
            "size": 1, "mode": 0o644,
        }
        artifact_closure = {
            "cache": identity, "graph": identity,
            "sources": {key: identity for key in (
                "benchmark", "decoder", "dispatch")},
            "objects": {key: identity for key in ("benchmark", "decoder")},
            "archive": identity, "binary": identity,
        }
        canonical_provenance = {
            "validated_cache": dict(REQUIRED_CANDIDATE_CACHE),
            "validated_compile_commands": {
                "schema": COMPILE_COMMANDS_SCHEMA,
                "implementation": "candidate",
                "profile": CANDIDATE_COMPILE_PROFILE,
                "validated_optimization": "-O3",
                "validated_openmp": True,
                "required_entries": [{"fixture": True}],
            },
            "archive_link_recipe_content": {"sha256": "6" * 64},
        }
        build = {
            "artifact_closure": artifact_closure,
            "canonical_production": {
                "schema": "leopard2-canonical-production-build/v1",
                "validator":
                    "exact-main/run_abba.py build_provenance schema v5",
                "provenance": canonical_provenance,
                "provenance_sha256": canonical_sha256(canonical_provenance),
            },
        }
        collector = {
            **identity, "relative_path":
                "experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py",
        }
        result_root = root / "attestation-result"
        raw_root = result_root / "raw"
        raw_root.mkdir(parents=True)
        raw_documents = {}
        raw_artifacts = {}
        for index, case in enumerate(attestation["cases"]):
            identifier = case["cell"]["identifier"]
            document = fake_attestation_output(case)
            path = raw_root / f"{index:04d}-{identifier}.json"
            write_json(path, document)
            raw_documents[identifier] = document
            raw_artifacts[identifier] = {
                "path": str(path.relative_to(result_root)),
                "size": path.stat().st_size, "sha256": file_sha256(path),
            }
        result = derive_attestation_result(
            stage_root, source, build, collector, raw_documents, raw_artifacts)
        noncanonical_build = json.loads(json.dumps(build))
        provenance = noncanonical_build["canonical_production"]["provenance"]
        provenance["validated_cache"]["LEO2_BACKEND_VARIANT"] = "scalar"
        noncanonical_build["canonical_production"][
            "provenance_sha256"] = canonical_sha256(provenance)
        try:
            derive_attestation_result(
                stage_root, source, noncanonical_build, collector,
                raw_documents, raw_artifacts)
        except PlanError:
            pass
        else:
            raise PlanError("noncanonical scalar-only build attested as AUTO")
        manifest_path = result_root / "manifest.json"
        write_json(manifest_path, result)
        validate_attestation_result_files(
            stage_root, manifest_path, source, build, collector)

        probe_case = attestation["cases"][0]
        probe = fake_attestation_output(probe_case)
        raw_mutations = {
            "K": lambda value: value["parameters"].update({
                "K": value["parameters"]["K"] + 1}),
            "seed": lambda value: value["parameters"].update({
                "seed": value["parameters"]["seed"] + 1}),
            "profile": lambda value: value["parameters"].update({
                "requested_profile": "low_v1"}),
            "field": lambda value: value["parameters"].update({
                "requested_field": "gf16"}),
            "requested-backend": lambda value: value["parameters"].update({
                "requested_backend": "scalar"}),
            "resolved-backend": lambda value: value["resolved"].update({
                "backend": "auto"}),
            "required-work-slots": lambda value: value["resolved"].update({
                "decode_required_work_slots":
                    value["resolved"]["decode_required_work_slots"] + 1}),
            "bytes": lambda value: value["parameters"].update({
                "shard_bytes": value["parameters"]["shard_bytes"] + 1}),
            "loss": lambda value: value["parameters"].update({
                "loss_count": (1 if value["parameters"]["loss_count"] != 1
                               else 2)}),
        }
        for label, mutate in raw_mutations.items():
            forged_output = json.loads(json.dumps(probe))
            mutate(forged_output)
            try:
                validate_attestation_output(forged_output, probe_case)
            except PlanError:
                continue
            raise PlanError(f"attestation accepted wrong {label}")

        for excluded_backend in ("avx512", "neon"):
            forged_output = fake_attestation_output(
                probe_case, backend="scalar")
            forged_output["resolved"]["backend"] = excluded_backend
            try:
                validate_attestation_output(forged_output, probe_case)
            except PlanError:
                pass
            else:
                raise PlanError(
                    f"attestation accepted excluded backend {excluded_backend}")

        def reject_output_mutation(label: str, mutate) -> None:
            fixture = root / ("bad-" + label)
            shutil.copytree(result_root, fixture)
            mutate(fixture)
            try:
                validate_attestation_result_files(
                    stage_root, fixture / "manifest.json", source, build, collector)
            except PlanError:
                return
            raise PlanError(f"adversarial attestation {label} was accepted")

        first_record = result["records"][0]
        first_raw = Path(first_record["raw"]["path"])
        first_case = next(item for item in attestation["cases"]
                          if item["cell"]["identifier"] ==
                          first_record["identifier"])
        def wrong_raw(path: Path, key: str, value: object) -> None:
            document = load_json(path / first_raw)
            document["resolved"][key] = value
            target = path / first_raw
            target.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n",
                              encoding="utf-8")
            def reseal(value_doc):
                record = next(item for item in value_doc["records"]
                              if item["identifier"] == first_record["identifier"])
                record["raw"]["size"] = target.stat().st_size
                record["raw"]["sha256"] = file_sha256(target)
                record["selected_decode_path"] = document["resolved"][
                    "selected_decode_path"]
                record["selected_decode_rule"] = document["resolved"][
                    "selected_decode_rule"]
                record["benchmark_payload_sha256"] = canonical_sha256(document)
            adversarial_resign(path / "manifest.json", reseal)

        wrong_path = "materialized" if first_case[
            "expected_selected_decode_path"] == "generic" else "generic"
        reject_output_mutation("wrong-path", lambda path:
                               wrong_raw(path, "selected_decode_path", wrong_path))
        wrong_rule = "workspace_materialized" if first_case[
            "expected_selected_decode_rule"] == "balanced_generic" else \
            "balanced_generic"
        reject_output_mutation("wrong-rule", lambda path:
                               wrong_raw(path, "selected_decode_rule", wrong_rule))
        reject_output_mutation("missing", lambda path:
                               (path / first_raw).unlink())
        def duplicate(path: Path) -> None:
            def mutate(value):
                value["records"].append(value["records"][0])
            adversarial_resign(path / "manifest.json", mutate)
        reject_output_mutation("duplicate", duplicate)
        reject_output_mutation("extra", lambda path:
                               (path / "extra.json").write_text("{}\n", encoding="utf-8"))
        def wrong_commit(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value.update({"candidate_commit": "6" * 40}))
        reject_output_mutation("commit", wrong_commit)
        def wrong_source_identity(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value["source_identity"].update({
                                   "head": "8" * 40,
                               }))
        reject_output_mutation("source-identity", wrong_source_identity)
        def wrong_benchmark_identity(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value["build_identity"]["artifact_closure"][
                                   "binary"].update({
                                   "sha256": "9" * 64,
                               }))
        reject_output_mutation("benchmark-identity", wrong_benchmark_identity)
        def wrong_stage(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value.update({"stage_content_sha256": "7" * 64}))
        reject_output_mutation("stage", wrong_stage)
        def resign_record(path: Path) -> None:
            adversarial_resign(path / "manifest.json", lambda value:
                               value["records"][0].update({
                                   "selected_decode_path": "generic",
                               }))
        reject_output_mutation("resigned-output", resign_record)

        broadened = root / "broadened-stage"
        shutil.copytree(stage_root, broadened)
        spec_path = broadened / "path-attestation.json"
        def broaden(value):
            target = next(item for item in value["cases"]
                          if item["expected_selected_decode_path"] == "not_generic")
            target["expected_selected_decode_path"] = "generic"
            target["expected_selected_decode_rule"] = "balanced_generic"
        adversarial_resign(spec_path, broaden)
        stage_path = broadened / "stage.json"
        def reseal_stage(value):
            target = next(item for item in value["artifacts"]
                          if item["path"] == "path-attestation.json")
            target["size"] = spec_path.stat().st_size
            target["sha256"] = file_sha256(spec_path)
        adversarial_resign(stage_path, reseal_stage)
        try:
            validate_stage(first, broadened, manifests, scopes)
        except PlanError:
            pass
        else:
            raise PlanError("re-signed broadened AUTO dispatch stage was accepted")
    print("balanced promotion plan self-test passed: exact-cell AUTO attestation + adversarial results")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    generate_parser = commands.add_parser("generate")
    generate_parser.add_argument("--output", required=True, type=Path)
    validate_parser = commands.add_parser("validate")
    validate_parser.add_argument("--input", required=True, type=Path)
    validate_parser.add_argument("--plan", type=Path)
    select_parser = commands.add_parser("select")
    select_parser.add_argument("--plan", required=True, type=Path)
    select_parser.add_argument("--gate-manifest", action="append", required=True,
                               type=Path)
    select_parser.add_argument("--output", required=True, type=Path)
    advance_parser = commands.add_parser("advance")
    advance_parser.add_argument("--plan", required=True, type=Path)
    advance_parser.add_argument("--survivors", required=True, type=Path)
    advance_parser.add_argument("--output", required=True, type=Path)
    collect_parser = commands.add_parser("collect-attestation")
    collect_parser.add_argument("--plan", required=True, type=Path)
    collect_parser.add_argument("--stage", required=True, type=Path)
    collect_parser.add_argument("--source-root", required=True, type=Path)
    collect_parser.add_argument("--binary", default="build-release/bench_leopard2")
    collect_parser.add_argument("--output", required=True, type=Path)
    verify_parser = commands.add_parser("verify-attestation")
    verify_parser.add_argument("--plan", required=True, type=Path)
    verify_parser.add_argument("--stage", required=True, type=Path)
    verify_parser.add_argument("--source-root", required=True, type=Path)
    verify_parser.add_argument("--binary", default="build-release/bench_leopard2")
    verify_parser.add_argument("--manifest", required=True, type=Path)
    commands.add_parser("self-test")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "generate":
        plan = generate_plan(args.output.absolute())
        result = {"content_sha256": plan["content_sha256"],
                  "initial_case_count": plan["initial_case_count"],
                  "initial_timed_child_count": plan["initial_timed_child_count"]}
    elif args.command == "validate":
        if args.plan is None:
            plan = validate_plan(args.input.absolute())
            result = {"content_sha256": (load_json(
                args.input.absolute() / "plan.json"))["content_sha256"],
                "initial_case_count": plan["initial_case_count"]}
        else:
            stage = validate_stage(args.plan.absolute(), args.input.absolute())
            result = {"content_sha256": (load_json(
                args.input.absolute() / "stage.json"))["content_sha256"],
                "timed_case_count": stage["timed_case_count"]}
    elif args.command == "select":
        survivor = select_survivors(args.plan.absolute(), args.gate_manifest,
                                    args.output.absolute())
        result = {"content_sha256": survivor["content_sha256"],
                  "survivor_count": len(survivor["survivor_cells"]),
                  "refinement_count": len(survivor["required_refinement_cells"])}
    elif args.command == "advance":
        survivor_path = args.survivors.resolve(strict=True)
        validate_survivors(args.plan.absolute(), survivor_path)
        stage = materialize_stage(args.plan.absolute(), load_json(survivor_path),
                                  args.output.absolute())
        result = {"content_sha256": stage["content_sha256"],
                  "timed_case_count": stage["timed_case_count"],
                  "timed_child_count": stage["timed_child_count"]}
    elif args.command == "collect-attestation":
        attestation = collect_attestation(
            args.plan.absolute(), args.stage.absolute(),
            args.source_root.absolute(), args.binary, args.output.absolute())
        result = {
            "content_sha256": attestation["content_sha256"],
            "record_count": len(attestation["records"]),
            "resolved_backend": attestation["resolved_backend"],
        }
    elif args.command == "verify-attestation":
        attestation = validate_attestation_result(
            args.plan.absolute(), args.stage.absolute(), args.manifest.absolute(),
            args.source_root.absolute(), args.binary)
        result = {
            "content_sha256": (load_json(args.manifest.absolute()))[
                "content_sha256"],
            "record_count": len(attestation["records"]),
            "resolved_backend": attestation["resolved_backend"],
        }
    else:
        self_test()
        return 0
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, json.JSONDecodeError, PlanError,
            common.EvidenceError) as error:
        raise SystemExit(f"balanced promotion plan failed: {error}") from error
