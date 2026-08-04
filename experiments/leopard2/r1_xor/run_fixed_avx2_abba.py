#!/usr/bin/env python3
"""Run the evidence-grade fixed-size AVX2 GF8 R=1 A/B campaign.

This runner never builds code.  Candidate and control must be equal-length
hard links to one lane-owned frozen benchmark inode.  By default that inode is
``bench_leopard2`` and the existing ordinary-public-API campaign is unchanged.
``--prevalidated-binding`` instead requires the separately compiled
``bench_leopard2_prevalidated_batch`` target and reports binding setup apart
from binding execution.  In either campaign the only candidate/control runtime
difference is ``--r1-fixed-avx2-mode 1`` versus ``0``; the benchmark forces the
older small-reduction experiment off before codec construction.
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
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-r1-fixed-avx2-abba/v1"
MANIFEST_SCHEMA = "leopard2-r1-fixed-avx2-manifest/v1"
CELL_SCHEMA = "leopard2-r1-fixed-avx2-cell/v1"
SUMMARY_SCHEMA = "leopard2-r1-fixed-avx2-summary/v1"
FAILURE_SCHEMA = "leopard2-r1-fixed-avx2-failure/v1"
BINDING_SCHEMA = "leopard2-r1-fixed-avx2-binding-abba/v1"
BINDING_MANIFEST_SCHEMA = "leopard2-r1-fixed-avx2-binding-manifest/v1"
BINDING_CELL_SCHEMA = "leopard2-r1-fixed-avx2-binding-cell/v1"
BINDING_SUMMARY_SCHEMA = "leopard2-r1-fixed-avx2-binding-summary/v1"
BINDING_FAILURE_SCHEMA = "leopard2-r1-fixed-avx2-binding-failure/v1"
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
PREVALIDATED_ENCODE_API = "leo2_encode_batch_binding_execute"
PREVALIDATED_DECODE_API = "leo2_decode_batch_binding_execute"
PREVALIDATED_ENCODE_SETUP_API = "leo2_encode_batch_binding_create"
PREVALIDATED_DECODE_SETUP_API = "leo2_decode_batch_binding_create"
PREVALIDATED_BUILD_FIELDS = (
    "prevalidated_batch_experiment",
    "high_t4_batch_diagnostic_disabled",
    "high_t4_batch_selected",
    "high_t8_one_block_extended_enabled",
    "high_t8_one_block_beyond_512_enabled",
    "high_t8_one_kilobyte_extension_enabled",
    "high_t8_tiny_binding_enabled",
    "high_t8_ragged_binding_enabled",
    "high_t8_one_block_selected",
    "high_t8_two_block_128_192_enabled",
    "high_t8_two_block_320_enabled",
    "high_t8_two_block_extended_enabled",
    "high_t8_two_block_selected",
)
PREVALIDATED_ATTRIBUTION_BUILD_FIELDS = (
    "r1_prevalidated_binding_attribution",
    "r1_encode_binding_setup_api",
    "r1_decode_binding_setup_api",
    "r1_binding_setup_reported_separately",
    "r1_binding_item_count",
)
_PREVALIDATED_BINDING = False


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
BASE_CAPTURE_SOURCES_AND_BUILDS = BASE.capture_sources_and_builds
BASE_ANALYZE_CELL = BASE.analyze_cell
BASE_RUNNER_DEPENDENCIES = tuple(BASE.RUNNER_DEPENDENCIES)
BASE_API_LANES = BASE.API_LANES


ORDINARY_BATCH_ROUTE_PROOF = {
    **BASE.ORDINARY_BATCH_ROUTE_PROOF,
    "required_schema": BENCHMARK_SCHEMA,
    "path_attestation_scope": (
        "Schema v23 identifies the ordinary zero-preflight one-item batch "
        "APIs and the actual leo2_encode/leo2_decode one-shot lanes. The "
        "fixed diagnostic additionally proves that the older small-reduction "
        "policy was disabled before immutable codec classification."),
}

PREVALIDATED_BINDING_ROUTE_PROOF = {
    "required_schema": BENCHMARK_SCHEMA,
    "required_batch_item_count": 1,
    "required_loss_count": 1,
    "required_build_fields": {
        "prevalidated_batch_experiment": True,
        "r1_prevalidated_binding_attribution": True,
        "r1_encode_binding_setup_api": PREVALIDATED_ENCODE_SETUP_API,
        "r1_decode_binding_setup_api": PREVALIDATED_DECODE_SETUP_API,
        "r1_binding_setup_reported_separately": True,
        "r1_binding_item_count": 1,
        "r1_timed_encode_api": PREVALIDATED_ENCODE_API,
        "r1_timed_reused_decode_api": PREVALIDATED_DECODE_API,
    },
    "required_parameter_fields": {"prevalidated_binding": True},
    "required_metric_fields": [
        "encode_binding_setup", "decode_binding_setup"],
    "path_attestation_scope": (
        "Schema v23 plus the compile-time prevalidated marker, explicit CLI "
        "attribution marker, exact binding create/execute API names, and "
        "separate setup metrics prove that encode_execution and "
        "decode_execution time immutable one-item binding executors."),
}

PREVALIDATED_API_LANES = {
    **BASE_API_LANES,
    "main_encode": BASE_API_LANES["main_encode"],
    "main_decode_one_call": BASE_API_LANES["main_decode_one_call"],
    "binding_encode_setup": {
        "implementation": "leopard2",
        "api": PREVALIDATED_ENCODE_SETUP_API,
        "item_count": 1,
        "metric": "encode_binding_setup",
    },
    "binding_encode_execution": {
        "implementation": "leopard2",
        "api": PREVALIDATED_ENCODE_API,
        "item_count": 1,
        "setup_excluded": True,
        "metric": "encode_execution",
    },
    "decode_plan_setup": BASE_API_LANES["leopard2_decode_plan_setup"],
    "binding_decode_setup": {
        "implementation": "leopard2",
        "api": PREVALIDATED_DECODE_SETUP_API,
        "item_count": 1,
        "metric": "decode_binding_setup",
    },
    "binding_decode_execution": {
        "implementation": "leopard2",
        "api": PREVALIDATED_DECODE_API,
        "item_count": 1,
        "setup_excluded": True,
        "execution_scope": "one_loss_direct_xor",
        "metric": "decode_execution",
    },
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
    command = common[:-2] + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--measure-one-shot-decode",
    ]
    if _PREVALIDATED_BINDING:
        command.append("--prevalidated-binding")
    return command + [
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
    metrics = result.get("metrics")
    require(result.get("schema") == BENCHMARK_SCHEMA and
            isinstance(parameters, Mapping) and
            parameters.get("r1_fixed_avx2_mode") == mode and
            isinstance(build, Mapping) and
            isinstance(metrics, Mapping) and
            build.get("r1_fixed_avx2_diagnostic_mode") == mode and
            build.get("r1_fixed_avx2_candidate_enabled") is True and
            build.get("r1_small_reduction_codec_enabled") is False and
            build.get("r1_encode_reduction_path") == expected_path and
            build.get("r1_decode_reduction_path") == expected_path and
            build.get("r1_fixed_avx2_selector_contract") ==
                SELECTOR_CONTRACT,
            "fixed AVX2 mode, disabled legacy experiment, selector contract, "
            "or execution route differs")

    if _PREVALIDATED_BINDING:
        require(parameters.get("prevalidated_binding") is True and
                all(name in build for name in PREVALIDATED_BUILD_FIELDS) and
                build.get("prevalidated_batch_experiment") is True and
                all(type(build.get(name)) is bool
                    for name in PREVALIDATED_BUILD_FIELDS) and
                all(build.get(name) == value for name, value in
                    PREVALIDATED_BINDING_ROUTE_PROOF[
                        "required_build_fields"].items()) and
                all(name in metrics for name in
                    PREVALIDATED_BINDING_ROUTE_PROOF[
                        "required_metric_fields"]),
                "prevalidated build identity, explicit binding opt-in, exact "
                "timed APIs, or separate binding setup metrics differ")
    else:
        require("prevalidated_binding" not in parameters and
                all(name not in build for name in PREVALIDATED_BUILD_FIELDS) and
                all(name not in build
                    for name in PREVALIDATED_ATTRIBUTION_BUILD_FIELDS) and
                "encode_binding_setup" not in metrics and
                "decode_binding_setup" not in metrics,
                "ordinary campaign unexpectedly selected a prevalidated "
                "binding benchmark")

    # Reuse the mature schema-v12 API/metrics validator after validating all
    # schema-v23 fields above.  The synthetic small-mode fields are local to
    # this adapter and only let the shared routine apply its existing metric,
    # API-scope, digest, and amortization checks.
    adapted = json.loads(json.dumps(result))
    encode_binding_setup = 0.0
    decode_binding_setup = 0.0
    if _PREVALIDATED_BINDING:
        encode_binding_setup = BASE.metric_median(
            metrics, "encode_binding_setup", "median_us", iterations,
            "samples_us")
        decode_binding_setup = BASE.metric_median(
            metrics, "decode_binding_setup", "median_us", iterations,
            "samples_us")
        for name in PREVALIDATED_BUILD_FIELDS:
            del adapted["build"][name]
        for name in PREVALIDATED_ATTRIBUTION_BUILD_FIELDS:
            del adapted["build"][name]
        del adapted["parameters"]["prevalidated_binding"]
        del adapted["metrics"]["encode_binding_setup"]
        del adapted["metrics"]["decode_binding_setup"]
        adapted["build"].update(
            ORDINARY_BATCH_ROUTE_PROOF["required_build_fields"])
    adapted["schema"] = "leopard2-benchmark-v12"
    adapted["build"]["r1_small_reduction_diagnostic_mode"] = mode
    adapted["build"]["r1_small_reduction_codec_enabled"] = mode == 1
    normalized = BASE_VALIDATE_RESULT(
        implementation, adapted, cell, iterations, warmup)
    normalized["mode"] = mode
    normalized["reduction_path"] = expected_path
    normalized["fixed_avx2_candidate_enabled"] = True
    normalized["small_reduction_codec_enabled"] = False
    if _PREVALIDATED_BINDING:
        ordinary_metrics = normalized["metrics_us"]
        encode_execution = ordinary_metrics["batch_encode"]
        decode_execution = ordinary_metrics["decode_execution"]
        codec_setup = ordinary_metrics["codec_setup"]
        plan_setup = ordinary_metrics["decode_plan_setup"]
        reuse = float(cell["reuse"])
        normalized["metrics_us"] = {
            "codec_setup": codec_setup,
            "binding_encode_setup": encode_binding_setup,
            "binding_encode_execution": encode_execution,
            "binding_encode_first_use":
                encode_binding_setup + encode_execution,
            "binding_encode_reuse_amortized":
                encode_execution + encode_binding_setup / reuse,
            "decode_plan_setup": plan_setup,
            "binding_decode_setup": decode_binding_setup,
            "binding_decode_execution": decode_execution,
            "binding_decode_first_use":
                plan_setup + decode_binding_setup + decode_execution,
            "binding_decode_reuse_amortized":
                decode_execution +
                (plan_setup + decode_binding_setup) / reuse,
        }
        normalized["api_lanes"] = PREVALIDATED_API_LANES
        normalized["prevalidated_binding"] = True
    return normalized


def candidate_target_name() -> str:
    return (
        "bench_leopard2_prevalidated_batch" if _PREVALIDATED_BINDING else
        "bench_leopard2")


def capture_sources_and_builds(
    options: argparse.Namespace,
) -> dict[str, Any]:
    """Bind provenance to the selected ordinary or prevalidated target."""
    if not _PREVALIDATED_BINDING:
        return BASE_CAPTURE_SOURCES_AND_BUILDS(options)

    candidate_source = options.source_root.resolve(strict=True)
    main_source = options.main_source_root.resolve(strict=True)
    try:
        sources = BASE.MAIN_SUPPORT.git_capture.capture_git_identities((
            (candidate_source, options.source_commit, False),
            (main_source, BASE.MAIN_COMMIT, True),
        ))
    except Exception as error:
        raise EvidenceError(f"source capture failed: {error}") from error
    require(sources[0].get("tree") == options.source_tree and
            sources[0].get("tracked_status") == "clean" and
            sources[1].get("tree") is not None and
            sources[1].get("tracked_status") == "clean",
            "source tree or tracked cleanliness differs")

    candidate_build = options.build_dir.resolve(strict=True)
    target = candidate_target_name()
    candidate_executable = candidate_build / target
    try:
        candidate_provenance = \
            BASE.BUILD_PROVENANCE.candidate_build_provenance(
                candidate_build, candidate_source, candidate_executable,
                target)
    except Exception as error:
        raise EvidenceError(f"candidate build provenance failed: {error}") \
            from error

    main_build = options.main_build_dir.resolve(strict=True)
    main_build_executable = main_build / "leopard_main_benchmark"
    main_archive = main_build / "libleopard_main_exact.a"
    main_specification = {
        "baseline_build_dir": str(main_build),
        "baseline_executable": str(main_build_executable),
        "baseline_archive": str(main_archive),
        "baseline_source_root": str(main_source),
        "candidate_source_root": str(candidate_source),
        "baseline_pure_avx2": True,
    }
    try:
        main_provenance = BASE.MAIN_SUPPORT.build_provenance(
            "baseline", main_specification)
    except Exception as error:
        raise EvidenceError(f"exact-main build provenance failed: {error}") \
            from error
    return {
        "candidate_source": sources[0],
        "main_source": sources[1],
        "candidate_build": candidate_provenance,
        "main_build": main_provenance,
    }


def analyze_binding_cell(
    cell: Mapping[str, Any],
    rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Analyze binding execution without naming it as an ordinary API call."""
    require(bool(rounds), "cell has no completed rounds")
    reference = rounds[0]["segments"][0]["invocations"][0][
        "normalized"]["digests"]
    same_source_names = (
        "codec_setup", "binding_encode_setup",
        "binding_encode_execution", "binding_encode_first_use",
        "binding_encode_reuse_amortized", "decode_plan_setup",
        "binding_decode_setup", "binding_decode_execution",
        "binding_decode_first_use", "binding_decode_reuse_amortized",
    )
    exact_main_comparisons = {
        "exact_main_leo_encode_vs_binding_encode_execution_"
        "semantically_different":
            ("legacy_encode", "binding_encode_execution"),
        "exact_main_leo_encode_vs_binding_encode_first_use_"
        "semantically_different":
            ("legacy_encode", "binding_encode_first_use"),
        "exact_main_leo_decode_vs_binding_decode_execution_"
        "semantically_different":
            ("legacy_decode", "binding_decode_execution"),
        "exact_main_leo_decode_vs_binding_decode_first_use_"
        "semantically_different":
            ("legacy_decode", "binding_decode_first_use"),
    }
    logs = {
        "candidate_vs_control": {name: [] for name in same_source_names},
        "candidate_vs_main": {
            name: [] for name in exact_main_comparisons},
    }
    for round_index, round_value in enumerate(rounds):
        require(round_value.get("round") == round_index,
                "round sequence is not contiguous and zero-based")
        isolation = round_value.get("isolation")
        require(isinstance(isolation, Mapping) and
                isolation.get("accepted") is True and
                isolation.get("delta", {}).get(
                    "reserved_sibling", {}).get("nonidle_jiffies") == 0,
                "contaminated round cannot be analyzed")
        segments = round_value.get("segments")
        require(isinstance(segments, list) and len(segments) == 2,
                "round does not contain two balanced segments")
        expected_segments = BASE.ROUND_SEGMENTS[round_index % 2]
        for segment_index, segment in enumerate(segments):
            require(isinstance(segment, Mapping),
                    "round segment is not an object")
            expected_order = expected_segments[segment_index]
            invocations = segment.get("invocations")
            require(isinstance(invocations, list) and
                    segment.get("order") == list(expected_order) and
                    [item.get("implementation")
                     if isinstance(item, Mapping) else None
                     for item in invocations] == list(expected_order) and
                    len(invocations) == 4 and
                    all(item["normalized"]["digests"] == reference
                        for item in invocations),
                    "balanced segment order or workload digests differ")
            baseline = segment.get("baseline")
            expected_baseline = (
                "main" if "main" in expected_order else "control")
            require(baseline == expected_baseline,
                    "segment baseline differs from its declared order")
            candidates = [item for item in invocations
                          if item["implementation"] == "candidate"]
            baselines = [item for item in invocations
                         if item["implementation"] == baseline]
            require(len(candidates) == len(baselines) == 2,
                    "balanced segment lacks two observations per lane")
            contrast = (
                "candidate_vs_main" if baseline == "main" else
                "candidate_vs_control")
            comparisons = (
                exact_main_comparisons if baseline == "main" else
                {name: (name, name) for name in same_source_names})
            for name, (baseline_name, candidate_name) in comparisons.items():
                candidate_log = statistics.mean(math.log(
                    item["normalized"]["metrics_us"][candidate_name])
                    for item in candidates)
                baseline_log = statistics.mean(math.log(
                    item["normalized"]["metrics_us"][baseline_name])
                    for item in baselines)
                logs[contrast][name].append(baseline_log - candidate_log)
    return {
        "cell": dict(cell),
        "digests": dict(reference),
        "candidate_vs_control": {
            name: BASE.confidence_interval(values)
            for name, values in logs["candidate_vs_control"].items()
        },
        "candidate_vs_main": {
            name: BASE.confidence_interval(values)
            for name, values in logs["candidate_vs_main"].items()
        },
        "exact_main_comparison_semantics": {
            "equivalent_public_api_semantics": False,
            "label": "cross_api_structural_reference_only",
            "main_calls": ["leo_encode", "leo_decode"],
            "leopard2_calls": [
                PREVALIDATED_ENCODE_API, PREVALIDATED_DECODE_API],
            "reason": (
                "Exact Leopard main has no reusable prevalidated-binding API; "
                "these ratios compare distinct call/setup contracts and are "
                "not ordinary-public-API speedup claims."),
        },
    }


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
        "schema": BASE.SCHEMA,
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
            "shared_prevalidated_binding_argument":
                _PREVALIDATED_BINDING,
            "small_reduction_forced_off_before_codec_creation": True,
            "only_runtime_difference": True,
            "candidate_control_equal_argument_lengths": True,
        },
        "api_lane_contract": BASE.API_LANES,
        ("prevalidated_binding_route_proof" if _PREVALIDATED_BINDING else
         "ordinary_batch_route_proof"): (
            PREVALIDATED_BINDING_ROUTE_PROOF if _PREVALIDATED_BINDING else
            ORDINARY_BATCH_ROUTE_PROOF),
        "benchmark_api_scope": (
            "Schema v23 times leo2_encode_batch_binding_execute and "
            "leo2_decode_batch_binding_execute while reporting both binding "
            "creation costs separately. Exact-main leo_encode/leo_decode "
            "ratios are cross-API structural references, not ordinary API "
            "speedup claims."
            if _PREVALIDATED_BINDING else
            "Schema v23 times ordinary leo2_encode_batch(item_count=1), "
            "leo2_decode_plan_execute_batch(item_count=1), and distinct "
            "leo2_encode/leo2_decode one-shot calls without prevalidated "
            "bindings or preflight scratch."),
        "exact_main_comparison_semantics": (
            {
                "equivalent_public_api_semantics": False,
                "label": "cross_api_structural_reference_only",
                "main_calls": ["leo_encode", "leo_decode"],
                "leopard2_calls": [
                    PREVALIDATED_ENCODE_API, PREVALIDATED_DECODE_API],
            } if _PREVALIDATED_BINDING else {
                "equivalent_public_api_semantics": True,
                "label": "ordinary_public_api_comparison",
            }),
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


def configure_base(
    main_sha256: str | None = None,
    prevalidated_binding: bool = False,
) -> None:
    """Install this campaign's policy into the shared evidence engine."""
    global _PREVALIDATED_BINDING
    _PREVALIDATED_BINDING = prevalidated_binding
    BASE.SCHEMA = BINDING_SCHEMA if prevalidated_binding else SCHEMA
    BASE.MANIFEST_SCHEMA = (
        BINDING_MANIFEST_SCHEMA if prevalidated_binding else MANIFEST_SCHEMA)
    BASE.CELL_SCHEMA = (
        BINDING_CELL_SCHEMA if prevalidated_binding else CELL_SCHEMA)
    BASE.SUMMARY_SCHEMA = (
        BINDING_SUMMARY_SCHEMA if prevalidated_binding else SUMMARY_SCHEMA)
    BASE.FAILURE_SCHEMA = (
        BINDING_FAILURE_SCHEMA if prevalidated_binding else FAILURE_SCHEMA)
    BASE.MODE_SYMBOL = MODE_SYMBOL
    BASE.BENCHMARK_CPU = BENCHMARK_CPU
    BASE.RESERVED_SIBLING = RESERVED_SIBLING
    BASE.RUNNER_PATH = RUNNER_PATH
    BASE.RUNNER_DEPENDENCIES = (RUNNER_PATH,) + BASE_RUNNER_DEPENDENCIES
    BASE.ORDINARY_BATCH_ROUTE_PROOF = ORDINARY_BATCH_ROUTE_PROOF
    BASE.API_LANES = (
        PREVALIDATED_API_LANES if prevalidated_binding else BASE_API_LANES)
    BASE.validate_frozen_candidate_identity = validate_candidate_identity
    BASE.expected_reduction_path = expected_reduction_path
    BASE.campaign_cells = campaign_cells
    BASE.benchmark_command = benchmark_command
    BASE.validate_result = validate_result
    BASE.make_contract = make_contract
    BASE.capture_sources_and_builds = (
        capture_sources_and_builds if prevalidated_binding else
        BASE_CAPTURE_SOURCES_AND_BUILDS)
    BASE.analyze_cell = (
        analyze_binding_cell if prevalidated_binding else BASE_ANALYZE_CELL)
    if main_sha256 is not None:
        BASE.MAIN_SHA256 = main_sha256


def fixture_result(
    cell: Mapping[str, Any], mode: int,
    iterations: int = 5, warmup: int = 2,
    prevalidated_binding: bool = False,
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
    result = {
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
    if prevalidated_binding:
        result["parameters"]["prevalidated_binding"] = True
        result["build"].update({
            name: False for name in PREVALIDATED_BUILD_FIELDS})
        result["build"]["prevalidated_batch_experiment"] = True
        result["build"].update({
            "r1_prevalidated_binding_attribution": True,
            "r1_encode_binding_setup_api": PREVALIDATED_ENCODE_SETUP_API,
            "r1_decode_binding_setup_api": PREVALIDATED_DECODE_SETUP_API,
            "r1_binding_setup_reported_separately": True,
            "r1_binding_item_count": 1,
            "r1_timed_encode_api": PREVALIDATED_ENCODE_API,
            "r1_timed_reused_decode_api": PREVALIDATED_DECODE_API,
        })
        result["metrics"]["encode_binding_setup"] = setup(0.5)
        result["metrics"]["decode_binding_setup"] = setup(0.75)
    return result


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

    configure_base(prevalidated_binding=True)
    require(candidate_target_name() == "bench_leopard2_prevalidated_batch",
            "self-test selected the wrong prevalidated provenance target")
    binding_candidate = benchmark_command(
        "candidate", Path("/frozen/a/benchmark"), sample,
        BENCHMARK_CPU, 5, 2)
    binding_control = benchmark_command(
        "control", Path("/frozen/b/benchmark"), sample,
        BENCHMARK_CPU, 5, 2)
    require("--prevalidated-binding" in binding_candidate and
            "--prevalidated-binding" in binding_control and
            binding_candidate[-4:] == [
                "--r1-fixed-avx2-mode", "1", "--json", "-"] and
            binding_control[-4:] == [
                "--r1-fixed-avx2-mode", "0", "--json", "-"] and
            len(binding_candidate) == len(binding_control) and
            all(len(left) == len(right) for left, right in
                zip(binding_candidate, binding_control)),
            "self-test prevalidated binding argv attribution failed")

    binding_fixture = fixture_result(
        sample, 1, prevalidated_binding=True)
    binding_normalized = validate_result(
        "candidate", binding_fixture, sample, 5, 2)
    require(binding_normalized.get("prevalidated_binding") is True and
            binding_normalized["api_lanes"][
                "binding_encode_execution"]["api"] ==
                PREVALIDATED_ENCODE_API and
            binding_normalized["api_lanes"][
                "binding_decode_execution"]["api"] ==
                PREVALIDATED_DECODE_API and
            math.isclose(binding_normalized["metrics_us"][
                "binding_encode_first_use"], 2.5) and
            math.isclose(binding_normalized["metrics_us"][
                "binding_decode_first_use"], 7.75),
            "self-test binding setup/execution normalization changed")
    for mutation in (
            "compile_marker", "explicit_marker", "parameter_marker",
            "encode_api", "decode_api", "encode_setup_api",
            "decode_setup_api", "setup_separate", "item_count",
            "encode_setup_metric", "decode_setup_metric"):
        changed = json.loads(json.dumps(binding_fixture))
        if mutation == "compile_marker":
            changed["build"]["prevalidated_batch_experiment"] = False
        elif mutation == "explicit_marker":
            changed["build"]["r1_prevalidated_binding_attribution"] = False
        elif mutation == "parameter_marker":
            changed["parameters"]["prevalidated_binding"] = False
        elif mutation == "encode_api":
            changed["build"]["r1_timed_encode_api"] = "leo2_encode_batch"
        elif mutation == "decode_api":
            changed["build"]["r1_timed_reused_decode_api"] = \
                "leo2_decode_plan_execute_batch"
        elif mutation == "encode_setup_api":
            changed["build"]["r1_encode_binding_setup_api"] += "!"
        elif mutation == "decode_setup_api":
            changed["build"]["r1_decode_binding_setup_api"] += "!"
        elif mutation == "setup_separate":
            changed["build"]["r1_binding_setup_reported_separately"] = False
        elif mutation == "item_count":
            changed["build"]["r1_binding_item_count"] = 2
        elif mutation == "encode_setup_metric":
            del changed["metrics"]["encode_binding_setup"]
        else:
            del changed["metrics"]["decode_binding_setup"]
        try:
            validate_result("candidate", changed, sample, 5, 2)
        except EvidenceError:
            pass
        else:
            raise EvidenceError(
                "self-test accepted prevalidated binding mutation: " +
                mutation)

    configure_base()
    require(candidate_target_name() == "bench_leopard2" and
            "--prevalidated-binding" not in benchmark_command(
                "candidate", Path("/frozen/a/benchmark"), sample,
                BENCHMARK_CPU, 5, 2),
            "self-test failed to restore ordinary campaign behavior")

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
        "prevalidated_binding_contract": "passed",
    }, sort_keys=True))
    return 0


def parse_arguments(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--prevalidated-binding", action="store_true",
        help=("benchmark the explicit one-item prevalidated binding "
              "create/execute contract instead of ordinary batch APIs"))
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--control", type=Path)
    parser.add_argument("--main", type=Path)
    parser.add_argument("--candidate-sha256")
    parser.add_argument("--main-sha256")
    parser.add_argument("--build-dir", type=Path,
                        help=("provenance root containing bench_leopard2 or "
                              "bench_leopard2_prevalidated_batch"))
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
    configure_base(options.main_sha256, options.prevalidated_binding)
    return BASE.run_campaign(options)


if __name__ == "__main__":
    raise SystemExit(main())
