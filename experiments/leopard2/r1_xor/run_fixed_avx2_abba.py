#!/usr/bin/env python3
"""Run the evidence-grade fixed-size AVX2 GF8 R=1 A/B campaign.

This runner never builds code.  Candidate and control must be equal-length
hard links to one lane-owned frozen ``bench_leopard2`` inode.  Their only
runtime difference is ``--r1-fixed-avx2-mode 1`` versus ``0``; the benchmark
forces the older small-reduction experiment off before codec construction.
Exact Leopard main is measured in the same alternating MAAM/AMMA segments as
the same-source BAAB/ABBA comparison.  Each completed cell is an atomic resume
unit and every child is address-space limited.

The implementation deliberately reuses the reviewed provenance, isolation,
statistics, atomic-artifact, and resume machinery from
``run_small_reduction_abba.py``.  This file owns the fixed-candidate matrix,
runtime command, schema-v23 diagnostic checks, selector identity, and campaign
contract; importing the old runner does not modify or execute its campaign.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-r1-fixed-avx2-abba/v1"
MANIFEST_SCHEMA = "leopard2-r1-fixed-avx2-manifest/v1"
CELL_SCHEMA = "leopard2-r1-fixed-avx2-cell/v1"
SUMMARY_SCHEMA = "leopard2-r1-fixed-avx2-summary/v1"
FAILURE_SCHEMA = "leopard2-r1-fixed-avx2-failure/v1"
BENCHMARK_SCHEMA = "leopard2-benchmark-v23"
MODE_SYMBOL = "_ZN12_GLOBAL__N_1L24g_r1_fixed_avx2_xor_modeE"
BENCHMARK_CPU = 4
RESERVED_SIBLING = 20
TARGET_K = (3, 4, 5, 7, 8, 24, 31, 32, 63, 64, 127, 128, 224, 255)
TARGET_BYTES = (64, 256)
CONTROL_BYTES = (128, 512)
POSITIONS = ("first", "middle", "last")
SELECTOR_CONTRACT = (
    "LEGACY_HIGH_V1,GF8,AVX2,R=1,K=3..255,B=64|256,native_layout,"
    "auto_encode_decode")


def load_module(path: Path, name: str) -> Any:
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load campaign dependency: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER_PATH = Path(__file__).resolve()
REPOSITORY = RUNNER_PATH.parents[3]
BASE_PATH = REPOSITORY / (
    "experiments/leopard2/r1_xor/run_small_reduction_abba.py")
BASE = load_module(BASE_PATH, "r1_fixed_avx2_evidence_base")
EvidenceError = BASE.EvidenceError
require = BASE.require
BASE_VALIDATE_RESULT = BASE.validate_result
BASE_RUNNER_DEPENDENCIES = tuple(BASE.RUNNER_DEPENDENCIES)


ORDINARY_BATCH_ROUTE_PROOF = {
    **BASE.ORDINARY_BATCH_ROUTE_PROOF,
    "required_schema": BENCHMARK_SCHEMA,
    "path_attestation_scope": (
        "Schema v23 identifies the ordinary zero-preflight one-item batch "
        "APIs and the actual leo2_encode/leo2_decode one-shot lanes. The "
        "fixed diagnostic additionally proves that the older small-reduction "
        "policy was disabled before immutable codec classification."),
}


def validate_candidate_identity(
    commit: str,
    tree: str,
    executable_sha256: str,
) -> None:
    """Require explicit immutable identities; build provenance binds them."""
    require(re.fullmatch(r"[0-9a-f]{40}", commit or "") is not None and
            re.fullmatch(r"[0-9a-f]{40}", tree or "") is not None and
            re.fullmatch(r"[0-9a-f]{64}", executable_sha256 or "") is not
                None,
            "candidate identities must be full lowercase hashes")
    require(set(commit) != {"0"} and set(tree) != {"0"} and
            set(executable_sha256) != {"0"},
            "candidate identities cannot be all-zero sentinels")


def expected_reduction_path(k: int, shard_bytes: int, mode: int) -> str:
    require(mode in (0, 1), "fixed AVX2 diagnostic mode is invalid")
    if k == 2 and shard_bytes in (64, 256):
        return "k2_terminal"
    if mode == 1 and 3 <= k <= 255 and shard_bytes in TARGET_BYTES:
        return "fixed_avx2"
    return "pairwise"


def campaign_cells() -> list[dict[str, Any]]:
    specifications: list[tuple[str, int, int, str]] = []
    for k in TARGET_K:
        for shard_bytes in TARGET_BYTES:
            for position in POSITIONS:
                specifications.append(("target_fixed", k, shard_bytes,
                                       position))

    # K=2 exercises the adjacent mature terminal.  First and last are its two
    # distinct missing positions; "middle" would duplicate last.
    for shard_bytes in TARGET_BYTES:
        for position in ("first", "last"):
            specifications.append(("control_k2", 2, shard_bytes, position))

    # Exercise both byte boundaries for every target K while rotating loss
    # position deterministically.  Across each three consecutive K values all
    # exclusion shapes are represented without tripling the negative matrix.
    for k_index, k in enumerate(TARGET_K):
        for byte_index, shard_bytes in enumerate(CONTROL_BYTES):
            specifications.append((
                "control_bytes", k, shard_bytes,
                POSITIONS[(k_index + byte_index) % len(POSITIONS)]))

    cells: list[dict[str, Any]] = []
    for index, (role, k, shard_bytes, position) in enumerate(specifications):
        seed = BASE.seed_for_loss(
            k, position, 0x463158000000 + index * 0x1000)
        candidate_path = expected_reduction_path(k, shard_bytes, 1)
        control_path = expected_reduction_path(k, shard_bytes, 0)
        selected = candidate_path != control_path
        require(selected == (role == "target_fixed"),
                "cell role disagrees with the fixed AVX2 selector")
        cells.append({
            "id": (f"{role}-k{k}-r1-b{shard_bytes}-"
                   f"m{position}-q1"),
            "K": k,
            "R": 1,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": BASE.reuse_for_cell(k, shard_bytes),
            "loss": 1,
            "missing_position": position,
            "expected_missing_original_indices": [
                BASE.desired_loss(k, position)],
            "role": role,
            "candidate_selected": selected,
            "candidate_reduction_path": candidate_path,
            "control_reduction_path": control_path,
            "seed": seed,
        })

    role_counts = {
        role: sum(cell["role"] == role for cell in cells)
        for role in ("target_fixed", "control_k2", "control_bytes")
    }
    require(len(cells) == 116 and role_counts == {
                "target_fixed": 84,
                "control_k2": 4,
                "control_bytes": 28,
            } and
            sum(bool(cell["candidate_selected"]) for cell in cells) == 84 and
            len({cell["id"] for cell in cells}) == len(cells) and
            len({cell["seed"] for cell in cells}) == len(cells) and
            all(BASE.selected_loss(cell["K"], cell["seed"]) ==
                cell["expected_missing_original_indices"][0]
                for cell in cells),
            "fixed AVX2 campaign matrix is incomplete")
    return cells


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    common = [
        "/usr/bin/prlimit", f"--as={BASE.ADDRESS_SPACE_BYTES}",
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", "1",
        "--bytes", str(cell["bytes"]), "--loss", "1",
        "--batch", "1", "--reuse", str(cell["reuse"]),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]), "--json", "-",
    ]
    if implementation == "main":
        return common
    require(implementation in {"candidate", "control"},
            "unknown Leopard2 implementation lane")
    mode = "1" if implementation == "candidate" else "0"
    return common[:-2] + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--measure-one-shot-decode",
        "--r1-fixed-avx2-mode", mode, "--json", "-",
    ]


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    if implementation == "main":
        return BASE_VALIDATE_RESULT(
            implementation, result, cell, iterations, warmup)
    require(implementation in {"candidate", "control"} and
            isinstance(result, Mapping),
            "fixed AVX2 result lane is invalid")
    mode = 1 if implementation == "candidate" else 0
    expected_path = str(cell[
        "candidate_reduction_path" if mode == 1 else
        "control_reduction_path"])
    parameters = result.get("parameters")
    build = result.get("build")
    require(result.get("schema") == BENCHMARK_SCHEMA and
            isinstance(parameters, Mapping) and
            parameters.get("r1_fixed_avx2_mode") == mode and
            isinstance(build, Mapping) and
            build.get("r1_fixed_avx2_diagnostic_mode") == mode and
            build.get("r1_fixed_avx2_candidate_enabled") is True and
            build.get("r1_small_reduction_codec_enabled") is False and
            build.get("r1_encode_reduction_path") == expected_path and
            build.get("r1_decode_reduction_path") == expected_path and
            build.get("r1_fixed_avx2_selector_contract") ==
                SELECTOR_CONTRACT,
            "fixed AVX2 mode, disabled legacy experiment, selector contract, "
            "or execution route differs")

    # Reuse the mature schema-v12 API/metrics validator after validating all
    # schema-v23 fields above.  The synthetic small-mode fields are local to
    # this adapter and only let the shared routine apply its existing metric,
    # API-scope, digest, and amortization checks.
    adapted = json.loads(json.dumps(result))
    adapted["schema"] = "leopard2-benchmark-v12"
    adapted["build"]["r1_small_reduction_diagnostic_mode"] = mode
    adapted["build"]["r1_small_reduction_codec_enabled"] = mode == 1
    normalized = BASE_VALIDATE_RESULT(
        implementation, adapted, cell, iterations, warmup)
    normalized["mode"] = mode
    normalized["reduction_path"] = expected_path
    normalized["fixed_avx2_candidate_enabled"] = True
    normalized["small_reduction_codec_enabled"] = False
    return normalized


def make_contract(
    options: argparse.Namespace,
    cells: Sequence[Mapping[str, Any]],
    identities: Mapping[str, Mapping[str, Any]],
    shared: Mapping[str, Any],
    provenance: Mapping[str, Any],
    isa: Mapping[str, Any],
    dependencies: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    lane_owned_copy_proof = BASE.validate_lane_owned_copy_identities(
        identities, provenance)
    mode_word = BASE.mode_word_identity(
        Path(str(identities["candidate"]["path"])))
    require(mode_word["value"] == 1,
            "shared production binary does not default fixed AVX2 mode to one")
    candidate_example = benchmark_command(
        "candidate", Path(str(identities["candidate"]["path"])),
        cells[0], options.cpu, options.iterations, options.warmup)
    control_example = benchmark_command(
        "control", Path(str(identities["control"]["path"])),
        cells[0], options.cpu, options.iterations, options.warmup)
    require(len(candidate_example) == len(control_example) and
            all(len(left) == len(right)
                for left, right in zip(candidate_example, control_example)),
            "candidate/control argv allocation shape differs")
    return {
        "schema": SCHEMA,
        "main_commit": BASE.MAIN_COMMIT,
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "binary_identities": dict(identities),
        "shared_inode_and_equal_path_proof": dict(shared),
        "lane_owned_copy_proof": lane_owned_copy_proof,
        "source_archive_build_provenance": dict(provenance),
        "pure_avx2_no_evex": dict(isa),
        "mode_word_default": mode_word,
        "runtime_attribution": {
            "candidate_arguments": ["--r1-fixed-avx2-mode", "1"],
            "control_arguments": ["--r1-fixed-avx2-mode", "0"],
            "small_reduction_forced_off_before_codec_creation": True,
            "only_runtime_difference": True,
            "candidate_control_equal_argument_lengths": True,
        },
        "api_lane_contract": BASE.API_LANES,
        "ordinary_batch_route_proof": ORDINARY_BATCH_ROUTE_PROOF,
        "benchmark_api_scope": (
            "Schema v23 times ordinary leo2_encode_batch(item_count=1), "
            "leo2_decode_plan_execute_batch(item_count=1), and distinct "
            "leo2_encode/leo2_decode one-shot calls without prevalidated "
            "bindings or preflight scratch."),
        "isolation": {
            "canonical_lock": str(BASE.LOCK_PATH),
            "benchmark_cpu": options.cpu,
            "reserved_sibling": options.sibling,
            "child_address_space_bytes": BASE.ADDRESS_SPACE_BYTES,
            "child_environment": BASE.CHILD_ENVIRONMENT,
        },
        "timing": {
            "rounds": options.rounds,
            "iterations": options.iterations,
            "warmup": options.warmup,
            "reuse_policy": {
                "formula": (
                    "clamp(nearest_integer(16777216 / (K * shard_bytes)), "
                    "64, 262144)"),
                "target_input_bytes_per_metric_sample":
                    BASE.TARGET_INPUT_BYTES_PER_METRIC_SAMPLE,
                "minimum": BASE.MIN_REUSE,
                "maximum": BASE.MAX_REUSE,
                "serialized_per_cell": True,
                "same_value_for_main_candidate_control": True,
            },
            "reuse_values": sorted({cell["reuse"] for cell in cells}),
            "round_segments": [
                [list(segment) for segment in BASE.ROUND_SEGMENTS[index % 2]]
                for index in range(options.rounds)],
            "method": (
                "alternating exact-main MAAM/AMMA and same-source "
                "BAAB/ABBA geometric contrasts; Student-t 95% interval over "
                "independent round log contrasts"),
        },
        "matrix": [dict(cell) for cell in cells],
        "matrix_role_semantics": (
            "target means fixed mode resolves fixed_avx2 while mode zero "
            "resolves the mature pairwise path. Controls require identical "
            "routes in both modes."),
        "runner_dependencies": list(dependencies),
    }


def configure_base(main_sha256: str | None = None) -> None:
    """Install this campaign's policy into the shared evidence engine."""
    BASE.SCHEMA = SCHEMA
    BASE.MANIFEST_SCHEMA = MANIFEST_SCHEMA
    BASE.CELL_SCHEMA = CELL_SCHEMA
    BASE.SUMMARY_SCHEMA = SUMMARY_SCHEMA
    BASE.FAILURE_SCHEMA = FAILURE_SCHEMA
    BASE.MODE_SYMBOL = MODE_SYMBOL
    BASE.BENCHMARK_CPU = BENCHMARK_CPU
    BASE.RESERVED_SIBLING = RESERVED_SIBLING
    BASE.RUNNER_PATH = RUNNER_PATH
    BASE.RUNNER_DEPENDENCIES = (RUNNER_PATH,) + BASE_RUNNER_DEPENDENCIES
    BASE.ORDINARY_BATCH_ROUTE_PROOF = ORDINARY_BATCH_ROUTE_PROOF
    BASE.validate_frozen_candidate_identity = validate_candidate_identity
    BASE.expected_reduction_path = expected_reduction_path
    BASE.campaign_cells = campaign_cells
    BASE.benchmark_command = benchmark_command
    BASE.validate_result = validate_result
    BASE.make_contract = make_contract
    if main_sha256 is not None:
        BASE.MAIN_SHA256 = main_sha256


def fixture_result(
    cell: Mapping[str, Any], mode: int,
    iterations: int = 5, warmup: int = 2,
) -> dict[str, Any]:
    def setup(value: float) -> dict[str, Any]:
        return {"median_us": value, "samples_us": [value] * iterations}

    def execution(value: float) -> dict[str, Any]:
        return {
            "median_us_per_batch_call": value,
            "samples_us_per_batch_call": [value] * iterations,
        }

    path = str(cell[
        "candidate_reduction_path" if mode == 1 else
        "control_reduction_path"])
    return {
        "schema": BENCHMARK_SCHEMA,
        "parameters": {
            "K": cell["K"], "R": 1, "shard_bytes": cell["bytes"],
            "loss_count": 1, "missing_original_indices":
                cell["expected_missing_original_indices"],
            "batch": 1, "reuse": cell["reuse"],
            "iterations": iterations, "warmup": warmup,
            "thread_count": 1, "seed": cell["seed"],
            "measure_one_shot_decode": True, "skip_legacy": True,
            "retain_samples": True,
            "requested_profile": "legacy_high_v1",
            "requested_field": "gf8", "requested_backend": "avx2",
            "force_generic_decode": False,
            "force_specialized_decode": False,
            "force_tiled_decode": False,
            "force_materialized_decode": False,
            "r1_fixed_avx2_mode": mode,
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "backend": "avx2", "thread_count": 1, "padded_side": 1,
        },
        "correctness": {"leopard2_round_trip": True},
        "workload_digests": {
            "algorithm": "fnv1a64", "original_data": "0" * 16,
            "transmitted_parity": "1" * 16,
            "recovered_originals": "2" * 16,
        },
        "build": {
            "r1_fixed_avx2_diagnostic_mode": mode,
            "r1_fixed_avx2_candidate_enabled": True,
            "r1_small_reduction_codec_enabled": False,
            "r1_encode_reduction_path": path,
            "r1_decode_reduction_path": path,
            "r1_fixed_avx2_selector_contract": SELECTOR_CONTRACT,
            **BASE.ORDINARY_BATCH_ROUTE_PROOF["required_build_fields"],
        },
        "memory": dict(BASE.ORDINARY_BATCH_ROUTE_PROOF[
            "required_memory_fields"]),
        "metrics": {
            "codec_setup": setup(1.0),
            "encode_execution": execution(2.0),
            "one_shot_encode": execution(2.5),
            "decode_plan_setup": setup(3.0),
            "decode_execution": execution(4.0),
            "decode_amortized_at_reuse": {
                "reuse_count": cell["reuse"],
                "derived_median_us_per_batch_call":
                    4.0 + 3.0 / float(cell["reuse"]),
            },
            "one_shot_decode_including_setup": execution(7.0),
        },
    }


def self_test() -> int:
    configure_base()
    validate_candidate_identity("1" * 40, "2" * 40, "3" * 64)
    for identity in (
            ("0" * 40, "2" * 40, "3" * 64),
            ("1" * 40, "bad", "3" * 64),
            ("1" * 40, "2" * 40, "3" * 63)):
        try:
            validate_candidate_identity(*identity)
        except EvidenceError:
            pass
        else:
            raise EvidenceError("self-test accepted an invalid identity")

    cells = campaign_cells()
    require(len(cells) == 116 and
            sum(cell["candidate_selected"] for cell in cells) == 84 and
            {cell["K"] for cell in cells if cell["role"] == "target_fixed"} ==
                set(TARGET_K) and
            {cell["bytes"] for cell in cells
             if cell["role"] == "target_fixed"} == set(TARGET_BYTES) and
            all({cell["missing_position"] for cell in cells
                 if cell["role"] == "target_fixed" and
                    cell["K"] == k and cell["bytes"] == shard_bytes} ==
                set(POSITIONS)
                for k in TARGET_K for shard_bytes in TARGET_BYTES),
            "self-test matrix coverage failed")
    require(expected_reduction_path(3, 64, 1) == "fixed_avx2" and
            expected_reduction_path(255, 256, 1) == "fixed_avx2" and
            expected_reduction_path(3, 64, 0) == "pairwise" and
            expected_reduction_path(3, 128, 1) == "pairwise" and
            expected_reduction_path(2, 64, 1) == "k2_terminal",
            "self-test selector boundary failed")

    sample = cells[0]
    candidate = benchmark_command(
        "candidate", Path("/frozen/a/benchmark"), sample,
        BENCHMARK_CPU, 5, 2)
    control = benchmark_command(
        "control", Path("/frozen/b/benchmark"), sample,
        BENCHMARK_CPU, 5, 2)
    require(candidate[-4:] == [
                "--r1-fixed-avx2-mode", "1", "--json", "-"] and
            control[-4:] == [
                "--r1-fixed-avx2-mode", "0", "--json", "-"] and
            "--r1-small-reduction-mode" not in candidate and
            "--r1-small-reduction-mode" not in control and
            len(candidate) == len(control) and
            all(len(left) == len(right)
                for left, right in zip(candidate, control)),
            "self-test runtime attribution argv failed")

    fixture = fixture_result(sample, 1)
    normalized = validate_result("candidate", fixture, sample, 5, 2)
    require(normalized["reduction_path"] == "fixed_avx2" and
            normalized["fixed_avx2_candidate_enabled"] is True and
            normalized["small_reduction_codec_enabled"] is False,
            "self-test schema-v23 fixture was mislabeled")
    for mutation in (
            "schema", "mode", "candidate_enabled", "small_enabled",
            "encode_route", "decode_route", "selector_contract"):
        changed = json.loads(json.dumps(fixture))
        if mutation == "schema":
            changed["schema"] = "leopard2-benchmark-v22"
        elif mutation == "mode":
            changed["build"]["r1_fixed_avx2_diagnostic_mode"] = 0
        elif mutation == "candidate_enabled":
            changed["build"]["r1_fixed_avx2_candidate_enabled"] = False
        elif mutation == "small_enabled":
            changed["build"]["r1_small_reduction_codec_enabled"] = True
        elif mutation == "encode_route":
            changed["build"]["r1_encode_reduction_path"] = "pairwise"
        elif mutation == "decode_route":
            changed["build"]["r1_decode_reduction_path"] = "pairwise"
        else:
            changed["build"]["r1_fixed_avx2_selector_contract"] += "!"
        try:
            validate_result("candidate", changed, sample, 5, 2)
        except EvidenceError:
            pass
        else:
            raise EvidenceError(
                f"self-test accepted diagnostic mutation: {mutation}")

    require(BASE.inspect_isa_disassembly(
        "  10: c5 fd ef c0 vpxor ymm0,ymm0,ymm0") == {
            "evex_prefixed_instruction_count": 0,
            "ymm_operand_instruction_count": 1,
        }, "self-test ISA parser changed")
    print(json.dumps({
        "status": "self-test-passed",
        "cells": len(cells),
        "targets": 84,
        "controls": 32,
    }, sort_keys=True))
    return 0


def parse_arguments(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--control", type=Path)
    parser.add_argument("--main", type=Path)
    parser.add_argument("--candidate-sha256")
    parser.add_argument("--main-sha256")
    parser.add_argument("--build-dir", type=Path,
                        help="provenance root containing bench_leopard2")
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--source-commit")
    parser.add_argument("--source-tree")
    parser.add_argument("--main-build-dir", type=Path,
                        help="provenance root containing leopard_main_benchmark")
    parser.add_argument("--main-source-root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cpu", type=int, default=BENCHMARK_CPU)
    parser.add_argument("--sibling", type=int, default=RESERVED_SIBLING)
    parser.add_argument("--rounds", type=int,
                        choices=tuple(BASE.T95), default=3)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--warmup", type=int, default=6)
    parser.add_argument("--timeout", type=float, default=60.0)
    options = parser.parse_args(arguments)
    if not options.self_test:
        required = {
            "candidate": options.candidate,
            "control": options.control,
            "main": options.main,
            "candidate-sha256": options.candidate_sha256,
            "main-sha256": options.main_sha256,
            "build-dir": options.build_dir,
            "source-root": options.source_root,
            "source-commit": options.source_commit,
            "source-tree": options.source_tree,
            "main-build-dir": options.main_build_dir,
            "main-source-root": options.main_source_root,
            "output": options.output,
        }
        missing = sorted(name for name, value in required.items()
                         if value is None)
        require(not missing,
                "missing required arguments: " + ", ".join(missing))
        require(re.fullmatch(r"[0-9a-f]{64}",
                             options.main_sha256 or "") is not None,
                "exact-main identity must be a full lowercase hash")
    return options


def main(arguments: Sequence[str] | None = None) -> int:
    options = parse_arguments(arguments)
    if options.self_test:
        return self_test()
    configure_base(options.main_sha256)
    return BASE.run_campaign(options)


if __name__ == "__main__":
    raise SystemExit(main())
