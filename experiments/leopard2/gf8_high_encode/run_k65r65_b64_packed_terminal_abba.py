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
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
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

BASE.SCHEMA = "leopard2-gf8-k65r65-b64-packed-terminal-abba/v1"
BASE.SUMMARY_SCHEMA = \
    "leopard2-gf8-k65r65-b64-packed-terminal-summary/v1"
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
BASE.CANONICAL_MAIN_SHA256 = \
    "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93"
BASE.MAX_ISOLATION_ATTEMPTS = 3
BASE.RUNNER_PATH = Path(__file__).resolve()
_INHERITED_RUNNER_DEPENDENCIES = tuple(BASE.RUNNER_DEPENDENCIES)
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    *_INHERITED_RUNNER_DEPENDENCIES,
)

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
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument(
        "--main-sha256", default=BASE.CANONICAL_MAIN_SHA256)
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


_BASE_ANALYZE = BASE.analyze


def _require_ci_inside_equivalence_band(
    ratio: Any, cell_id: str, metric_name: str,
) -> list[float]:
    BASE.require(isinstance(ratio, Mapping) and
                 isinstance(ratio.get("ci95"), list) and
                 len(ratio["ci95"]) == 2,
                 f"{cell_id} {metric_name} selector contrast is missing")
    lower, upper = ratio["ci95"]
    BASE.require(
        isinstance(lower, (int, float)) and not isinstance(lower, bool) and
        isinstance(upper, (int, float)) and not isinstance(upper, bool) and
        lower >= BASE.NEIGHBOR_EQUIVALENCE_LOWER and
        upper <= BASE.NEIGHBOR_EQUIVALENCE_UPPER,
        f"{cell_id} {metric_name} does not prove selector equivalence "
        f"inside [{BASE.NEIGHBOR_EQUIVALENCE_LOWER}, "
        f"{BASE.NEIGHBOR_EQUIVALENCE_UPPER}]: CI [{lower}, {upper}]")
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
        for name, (ratio, floor) in ratios.items():
            BASE.require(isinstance(ratio, Mapping) and
                         isinstance(ratio.get("ci95"), list) and
                         len(ratio["ci95"]) == 2,
                         f"target {name} contrast is missing")
            lower = ratio["ci95"][0]
            BASE.require(isinstance(lower, (int, float)) and
                         not isinstance(lower, bool) and lower >= floor,
                         f"target {name} lower CI {lower!r} < {floor}")
            lower_bounds[name] = lower
        result["target_acceptance_validation"] = {
            "control_floor": BASE.TARGET_CONTROL_FLOOR,
            "main_floor": BASE.TARGET_MAIN_FLOOR,
            "lower_ci95": lower_bounds,
            "accepted": True,
        }
        return result
    if cell.get("role") != "neighbor":
        return result

    cell_id = str(cell.get("id"))
    selector_cis = {
        "encode_execution": _require_ci_inside_equivalence_band(
            result.get("control_over_candidate"), cell_id,
            "encode_execution"),
        "one_shot_encode": _require_ci_inside_equivalence_band(
            result.get("one_shot_control_over_candidate"), cell_id,
            "one_shot_encode"),
    }
    result["neighbor_selector_inertness_validation"] = {
        "equivalence_band": [
            BASE.NEIGHBOR_EQUIVALENCE_LOWER,
            BASE.NEIGHBOR_EQUIVALENCE_UPPER,
        ],
        "ci95": selector_cis,
        "accepted": True,
    }

    if not cell.get("compare_main"):
        return result
    ratios = {
        "encode_execution": result.get("main_over_candidate"),
        "one_shot_encode": result.get("one_shot_main_over_candidate"),
    }
    context: dict[str, Any] = {}
    for name, ratio in ratios.items():
        BASE.require(isinstance(ratio, Mapping) and
                     isinstance(ratio.get("ci95"), list) and
                     len(ratio["ci95"]) == 2,
                     f"{cell_id} {name} exact-main contrast is missing")
        context[name] = {
            "speedup": ratio.get("speedup"),
            "ci95": list(ratio["ci95"]),
        }
    floor = cell.get("main_floor")
    if floor is not None:
        BASE.require(type(floor) is float and floor == BASE.RETAINED_MAIN_FLOOR,
                     f"{cell_id} exact-main floor is not established")
        for name, ratio in ratios.items():
            lower = ratio["ci95"][0]
            BASE.require(isinstance(lower, (int, float)) and
                         not isinstance(lower, bool) and lower >= floor,
                         f"{cell_id} {name} regressed against exact main: "
                         f"lower CI {lower!r} < {floor}")
    result["exact_main_context"] = {
        "policy": "gated" if floor is not None else "context_only",
        "floor": floor,
        "metrics": context,
        "accepted": True,
    }
    return result


BASE.analyze = analyze


if __name__ == "__main__":
    raise SystemExit(BASE.main())
