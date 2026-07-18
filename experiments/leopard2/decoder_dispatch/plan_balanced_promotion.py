#!/usr/bin/env python3
"""Generate and advance the exact-main-first balanced decoder matrix."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Iterable

import balanced_evidence_common as common


PLAN_SCHEMA = "leopard2-balanced-promotion-plan/v2"
MAIN_CELL_SCHEMA = "leopard2-main-compare-cell-list/v2"
SURVIVOR_SCHEMA = "leopard2-balanced-promotion-survivors/v2"
STAGE_SCHEMA = "leopard2-balanced-promotion-stage/v3"
ATTESTATION_SCHEMA = "leopard2-balanced-auto-path-attestation/v3"
ATTESTATION_RESULT_SCHEMA = "leopard2-balanced-auto-path-result/v2"
EXACT_MANIFEST_SCHEMA = "leopard2-main-compare-manifest/v4"
EXACT_MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
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
EXACT_RAW_SCHEMA = "leopard2-main-compare-raw/v4"
REQUIRED_CANDIDATE_CACHE = {
    "CMAKE_BUILD_TYPE": "Release",
    "ENABLE_OPENMP": "ON",
    "LEO2_BACKEND_VARIANT": "auto",
    "LEO2_BUILD_BENCHMARKS": "ON",
    "LEO2_BUILD_FUZZERS": "OFF",
    "LEO2_BUILD_TESTS": "OFF",
    "LEO2_ENABLE_CUDA": "OFF",
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
                "select verifies every schema-v4 exact-main gate manifest and "
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


def _cpu_list_count(value: object, label: str) -> int:
    require(isinstance(value, str), f"{label} is not a CPU-list string")
    result = set()
    for component in value.split(","):
        component = component.strip()
        require(component and re.fullmatch(r"[0-9]+(?:-[0-9]+)?", component),
                f"{label} is malformed")
        bounds = [int(item) for item in component.split("-", 1)]
        first, last = (bounds[0], bounds[-1])
        require(first <= last, f"{label} range is reversed")
        result.update(range(first, last + 1))
    return len(result)


def _normalized_cpu_policy(value: object, label: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{label} is not an object")
    cpuinfo = value.get("cpuinfo")
    topology = value.get("topology")
    frequency = value.get("frequency_policy")
    require(isinstance(cpuinfo, dict) and isinstance(topology, dict) and
            isinstance(frequency, dict), f"{label} policy is incomplete")
    # Logical CPU and core identifiers vary between the independently pinned
    # shards.  Retain the CPU signature, sibling cardinality,
    # frequency policy, and online state instead of those incidental IDs.
    normalized_cpuinfo = {
        key: item for key, item in cpuinfo.items() if key != "processor"
    }
    thread_list = topology.get("thread_siblings_list")
    core_list = topology.get("core_siblings_list")
    return {
        "online": value.get("online"),
        "cpuinfo": normalized_cpuinfo,
        "topology": {
            "thread_sibling_count": _cpu_list_count(
                thread_list, f"{label} thread siblings"),
            "core_sibling_count": _cpu_list_count(
                core_list, f"{label} core siblings"),
        },
        "frequency_policy": frequency,
    }


def _normalize_bound_paths(value: object, replacements: tuple[tuple[str, str], ...]) -> object:
    """Remove volatile filesystem metadata while retaining exact byte identity."""
    if isinstance(value, dict):
        return {
            key: _normalize_bound_paths(item, replacements)
            for key, item in sorted(value.items())
            if key not in {"device", "inode", "mtime_ns"}
        }
    if isinstance(value, list):
        return [_normalize_bound_paths(item, replacements) for item in value]
    if isinstance(value, str):
        result = value
        for original, marker in replacements:
            result = result.replace(original, marker)
        return result
    return value


def selection_scope_from_verified_bundle(
    manifest: dict[str, Any], raw: dict[str, Any]
) -> dict[str, Any]:
    """Derive the one environment in which every gate cell is meaningful."""
    require(raw.get("schema") == EXACT_RAW_SCHEMA,
            "exact-main raw bundle is not hardened schema v4")
    identities = manifest.get("identities")
    host = manifest.get("host")
    require(isinstance(identities, dict) and isinstance(host, dict),
            "exact-main scope identity is incomplete")
    candidate_build = identities.get("candidate_build")
    candidate_source = identities.get("candidate_source")
    require(isinstance(candidate_build, dict) and isinstance(candidate_source, dict),
            "exact-main candidate build/source scope is absent")
    cache = candidate_build.get("validated_cache")
    compile_semantics = candidate_build.get("validated_compile_commands")
    require(isinstance(cache, dict) and all(
        cache.get(key) == expected
        for key, expected in REQUIRED_CANDIDATE_CACHE.items()),
        "exact-main candidate is not the canonical Release AUTO build")
    require(isinstance(compile_semantics, dict) and
            compile_semantics.get("validated_optimization") == "-O3" and
            compile_semantics.get("validated_openmp") is True and
            isinstance(compile_semantics.get("required_source_object_pairs"), list) and
            compile_semantics["required_source_object_pairs"] and
            isinstance(candidate_build.get("archive_link_recipe_content"), dict) and
            isinstance(candidate_build.get("archive_link_tool_invocations"), dict) and
            isinstance(candidate_build.get("ranlib"), dict),
            "exact-main candidate lacks hardened compile/link provenance")
    source_root = candidate_source.get("path")
    build_root = candidate_build.get("build_dir")
    require(isinstance(source_root, str) and isinstance(build_root, str),
            "exact-main candidate path roots are absent")

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

    allowed = host.get("allowed_cpu_set_at_launch")
    online = host.get("online_cpu_set")
    require(isinstance(allowed, list) and allowed and
            all(type(item) is int and item >= 0 for item in allowed) and
            isinstance(online, list) and online and
            all(type(item) is int and item >= 0 for item in online),
            "exact-main CPU-set scope is invalid")
    replacements = ((source_root, "$SOURCE"), (build_root, "$BUILD"))
    build_scope = _normalize_bound_paths(candidate_build, replacements)
    require(isinstance(build_scope, dict), "normalized build scope is invalid")
    scope = {
        "schema": "leopard2-balanced-evidence-scope/v1",
        "host": {
            "system": host.get("system"),
            "allowed_cpu_set_at_launch": allowed,
            "online_cpu_set": online,
            "benchmark_cpu_class": _normalized_cpu_policy(
                host.get("benchmark_cpu"), "benchmark CPU"),
            "reserved_sibling_class": _normalized_cpu_policy(
                host.get("reserved_sibling"), "reserved sibling"),
            "turbo_and_pstate": host.get("turbo_and_pstate"),
        },
        "compiler_and_build": build_scope,
        "candidate_source": _normalize_bound_paths(candidate_source, replacements),
        "resolved_auto_backend": resolved_backend,
        "forced_confirmation_backends": list(
            BACKENDS[:BACKENDS.index(resolved_backend) + 1]),
        "excluded_backends": dict(EXCLUDED_CAMPAIGN_BACKENDS),
    }
    return scope


def validate_evidence_scope(scope: object) -> dict[str, Any]:
    require(isinstance(scope, dict) and set(scope) == {
        "schema", "host", "compiler_and_build", "candidate_source",
        "resolved_auto_backend", "forced_confirmation_backends",
        "excluded_backends",
    } and scope.get("schema") == "leopard2-balanced-evidence-scope/v1",
            "gate evidence scope shape differs")
    backend = scope.get("resolved_auto_backend")
    require(backend in CAMPAIGN_BACKENDS and
            scope.get("forced_confirmation_backends") ==
                list(BACKENDS[:BACKENDS.index(backend) + 1]) and
            scope.get("excluded_backends") == EXCLUDED_CAMPAIGN_BACKENDS,
            "gate evidence backend coverage declaration differs")
    host = scope.get("host")
    build = scope.get("compiler_and_build")
    source = scope.get("candidate_source")
    require(isinstance(host, dict) and set(host) == {
        "system", "allowed_cpu_set_at_launch", "online_cpu_set",
        "benchmark_cpu_class", "reserved_sibling_class", "turbo_and_pstate",
    } and isinstance(host["system"], dict) and
            isinstance(host["benchmark_cpu_class"], dict) and
            isinstance(host["reserved_sibling_class"], dict) and
            isinstance(host["allowed_cpu_set_at_launch"], list) and
            host["allowed_cpu_set_at_launch"] and
            isinstance(host["online_cpu_set"], list) and host["online_cpu_set"],
            "gate evidence host/topology scope differs")
    require(isinstance(build, dict) and isinstance(source, dict) and
            isinstance(build.get("validated_cache"), dict) and
            all(build["validated_cache"].get(key) == expected
                for key, expected in REQUIRED_CANDIDATE_CACHE.items()) and
            isinstance(build.get("compiler"), dict) and
            isinstance(build.get("validated_compile_commands"), dict) and
            build["validated_compile_commands"].get(
                "validated_optimization") == "-O3" and
            build["validated_compile_commands"].get("validated_openmp") is True and
            isinstance(build.get("archive_link_recipe_content"), dict) and
            isinstance(source.get("head"), str) and
            re.fullmatch(r"[0-9a-f]{40}", source["head"]) is not None,
            "gate evidence compiler/CMake/build scope differs")
    return scope


def verify_exact_manifest(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    runner = Path(__file__).resolve().parents[1] / "main_compare" / "run_abba.py"
    command = [sys.executable, str(runner), "verify", "--manifest", str(path),
               "--no-current-input-check"]
    completed = subprocess.run(command, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True)
    require(completed.returncode == 0 and not completed.stderr,
            f"exact-main verifier rejected {path}: {completed.stderr.strip()}")
    document = load_json(path)
    require(isinstance(document, dict) and
            document.get("schema") == EXACT_MANIFEST_SCHEMA and
            document.get("valid") is True,
            f"exact-main manifest is not valid schema v4: {path}")
    raw_record = document.get("raw")
    require(isinstance(raw_record, dict) and
            isinstance(raw_record.get("path"), str),
            f"exact-main manifest raw identity is incomplete: {path}")
    relative = Path(raw_record["path"])
    require(not relative.is_absolute() and ".." not in relative.parts,
            "exact-main raw path escapes its evidence directory")
    raw_path = path.parent / relative
    raw = load_json(raw_path)
    require(isinstance(raw, dict), "exact-main raw bundle is not an object")
    return document, validate_evidence_scope(
        selection_scope_from_verified_bundle(document, raw))


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
        require(scope["candidate_source"]["head"] == candidate_commit,
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
        manifest, scope = verify_exact_manifest(resolved)
        manifests.append(manifest)
        scopes.append(scope)
        references.append({
            "path": str(resolved), "size": resolved.stat().st_size,
            "sha256": file_sha256(resolved), "payload_digest": manifest["digest"],
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
            require(manifest_path.is_file() and
                    manifest_path.stat().st_size == reference["size"] and
                    file_sha256(manifest_path) == reference["sha256"],
                    "survivor gate manifest bytes changed")
            manifest, scope = verify_exact_manifest(manifest_path)
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
    """Run the exact-main schema-v4 compile/archive/link semantic validator."""
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
            semantics.get("validated_optimization") == "-O3" and
            semantics.get("validated_openmp") is True and
            isinstance(provenance.get("archive_link_recipe_content"), dict),
            "canonical build provenance does not bind Release AUTO semantics")
    normalized = _normalize_bound_paths(
        provenance,
        ((str(source_root), "$SOURCE"), (str(build_root), "$BUILD")),
    )
    require(isinstance(normalized, dict), "canonical build scope is invalid")
    return {
        "schema": "leopard2-canonical-production-build/v1",
        "validator": "exact-main/run_abba.py build_provenance schema v4",
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
                "exact-main/run_abba.py build_provenance schema v4" and
            isinstance(canonical["provenance"], dict) and
            canonical_sha256(canonical["provenance"]) ==
                canonical["provenance_sha256"] and
            isinstance(canonical["provenance"].get("validated_cache"), dict) and
            all(canonical["provenance"]["validated_cache"].get(key) == expected
                for key, expected in REQUIRED_CANDIDATE_CACHE.items()) and
            isinstance(canonical["provenance"].get(
                "validated_compile_commands"), dict) and
            canonical["provenance"]["validated_compile_commands"].get(
                "validated_optimization") == "-O3" and
            canonical["provenance"]["validated_compile_commands"].get(
                "validated_openmp") is True and
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
    return {
        "schema": "leopard2-balanced-evidence-scope/v1",
        "host": {
            "system": {"system": "Linux", "machine": "x86_64"},
            "allowed_cpu_set_at_launch": list(range(8)),
            "online_cpu_set": list(range(8)),
            "benchmark_cpu_class": {"cpuinfo": {"model": "fixture"}},
            "reserved_sibling_class": {"cpuinfo": {"model": "fixture"}},
            "turbo_and_pstate": {"fixture": "stable"},
        },
        "compiler_and_build": {
            "validated_cache": dict(REQUIRED_CANDIDATE_CACHE),
            "compiler": {"sha256": "a" * 64},
            "validated_compile_commands": {
                "validated_optimization": "-O3", "validated_openmp": True,
            },
            "archive_link_recipe_content": {"sha256": "b" * 64},
        },
        "candidate_source": {"head": "1" * 40, "tree": "c" * 40},
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
                "candidate_source": {"head": "1" * 40},
            },
            "analysis": analysis,
        })
        references.append({
            "path": "/fixture/" + Path(item["path"]).name,
            "size": 1, "sha256": "2" * 64, "payload_digest": "3" * 64,
        })
        scopes.append(fake_evidence_scope(backend))
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
        reject_scope_mutation("compiler", lambda values:
                              values[-1]["compiler_and_build"]["compiler"].update({
                                  "sha256": "d" * 64}))
        reject_scope_mutation("CMake", lambda values:
                              values[-1]["compiler_and_build"][
                                  "validated_cache"].update({
                                      "LEO2_BACKEND_VARIANT": "scalar"}))
        reject_scope_mutation("uniform noncanonical CMake", lambda values:
                              [value["compiler_and_build"][
                                  "validated_cache"].update({
                                      "LEO2_BACKEND_VARIANT": "scalar"})
                               for value in values])
        reject_scope_mutation("resolved backend", lambda values:
                              values[-1].update({
                                  "resolved_auto_backend": "ssse3"}))
        for excluded_backend in ("avx512", "neon"):
            reject_scope_mutation(excluded_backend, lambda values, backend=excluded_backend:
                                  [value.update({
                                      "resolved_auto_backend": backend})
                                   for value in values])

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
            "head": "1" * 40, "tree": "4" * 40, "status": "clean",
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
                "validated_optimization": "-O3",
                "validated_openmp": True,
            },
            "archive_link_recipe_content": {"sha256": "6" * 64},
        }
        build = {
            "artifact_closure": artifact_closure,
            "canonical_production": {
                "schema": "leopard2-canonical-production-build/v1",
                "validator":
                    "exact-main/run_abba.py build_provenance schema v4",
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
