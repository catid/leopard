#!/usr/bin/env python3
"""Qualify the fixed T=16/B=64 three-pass AVX2 encode terminal.

The target is measured through both the direct ``leo2_encode`` entry point and
the one-item ``leo2_encode_batch`` wrapper.  A single immutable Leopard2
binary supplies candidate and control through the diagnostic selector; exact
Leopard main supplies the historical old-API reference.  Non-target cells
exercise the same binary without touching the selector and therefore measure
layout/noise rather than a different executable.
"""

from __future__ import annotations

import importlib.util
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_k5r5_b64_terminal_abba.py")
    specification = importlib.util.spec_from_file_location(
        "t16_b64_three_pass_evidence_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load evidence base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-t16-b64-three-pass-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-t16-b64-three-pass-summary/v1"
BASE.MODE_SYMBOL = "_ZN7leopard3ff8L30g_high_t16_b64_three_pass_modeE"
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.SHARED_MODE_DEFAULT = 2
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.TARGET_MAIN_FLOOR = 1.05
BASE.TARGET_ORDER = BASE.TARGET_ORDER * 3
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
)
BASE_BENCHMARK_COMMAND = BASE.benchmark_command


def confidence_interval(values: Sequence[float]) -> dict[str, Any]:
    BASE.require(len(values) in (3, 9),
                 "three or nine independent ABBA contrasts are required")
    critical = {3: BASE.T95_DF2, 9: 2.306004135204166}[len(values)]
    center = statistics.mean(values)
    half_width = critical * statistics.stdev(values) / math.sqrt(len(values))
    return {
        "speedup": math.exp(center),
        "ci95": [
            math.exp(center - half_width),
            math.exp(center + half_width),
        ],
        "round_log_ratios": list(values),
    }


BASE.confidence_interval = confidence_interval


def cells() -> list[dict[str, Any]]:
    values = (
        ("target-direct-k16-r16-b64", 16, 16, 64, "direct", "target"),
        ("target-q1-k16-r16-b64", 16, 16, 64, "q1", "target"),
        ("byte-k16-r16-b63-q1", 16, 16, 63, "q1", "neighbor"),
        ("byte-k16-r16-b65-q1", 16, 16, 65, "q1", "neighbor"),
        ("byte-k16-r16-b128-q1", 16, 16, 128, "q1", "neighbor"),
        ("shape-k15-r16-b64-q1", 15, 16, 64, "q1", "neighbor"),
        ("shape-k16-r15-b64-q1", 16, 15, 64, "q1", "neighbor"),
        ("transform-t8-k8-r8-b64-q1", 8, 8, 64, "q1", "neighbor"),
        ("transform-t32-k32-r32-b64-q1", 32, 32, 64, "q1", "neighbor"),
    )
    return [
        {
            "id": name,
            "K": original_count,
            "R": recovery_count,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": 8192,
            "metric": metric,
            "role": role,
            "seed": 0x54313600 + index,
        }
        for index, (
            name, original_count, recovery_count, shard_bytes, metric, role
        ) in enumerate(values)
    ]


BASE.cells = cells


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    command = BASE_BENCHMARK_COMMAND(
        implementation, executable, cell, cpu, iterations, warmup)
    if implementation == "main":
        return command
    command = [value for value in command if value != "--attest-source"]
    if cell["role"] == "target":
        command[-2:-2] = [
            "--t16-b64-three-pass-mode",
            "1" if implementation == "candidate" else "2",
        ]
    return command


BASE.benchmark_command = benchmark_command


def metric(result: Mapping[str, Any], name: str) -> float:
    metrics = result.get("metrics")
    BASE.require(isinstance(metrics, dict), "benchmark metrics are absent")
    value = metrics.get(name)
    BASE.require(isinstance(value, dict), f"{name} metric is absent")
    median = value.get("median_us_per_batch_call")
    BASE.require(isinstance(median, (int, float)) and
                 not isinstance(median, bool) and
                 math.isfinite(float(median)) and float(median) > 0,
                 f"{name} median is not finite and positive")
    return float(median)


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    del source_commit, source_tree
    BASE.require(isinstance(result, dict), "benchmark output is not an object")
    expected = {
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": cell["bytes"], "loss_count": 1,
        "batch": 1, "reuse": cell["reuse"],
        "iterations": iterations, "warmup": warmup,
        "thread_count": 1, "seed": cell["seed"],
    }
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    BASE.require(isinstance(parameters, dict) and
                 all(parameters.get(key) == value
                     for key, value in expected.items()),
                 "benchmark parameters differ from the frozen cell")
    BASE.require(isinstance(resolved, dict) and
                 resolved.get("profile") == "legacy_high_v1" and
                 resolved.get("field") == "gf8" and
                 resolved.get("thread_count") == 1,
                 "benchmark resolved a different profile/field/thread count")
    BASE.require(isinstance(digests, dict) and
                 digests.get("algorithm") == "fnv1a64" and
                 all(isinstance(digests.get(name), str) and
                     len(digests[name]) == 16 for name in (
                         "original_data", "transmitted_parity",
                         "recovered_originals")),
                 "benchmark workload digests are incomplete")
    if implementation == "main":
        BASE.require(result.get("schema") == "leopard-main-benchmark-v1" and
                     isinstance(correctness, dict) and
                     correctness.get("round_trip") is True,
                     "exact-main identity or round trip failed")
        build = result.get("build")
        BASE.require(isinstance(build, dict) and
                     build.get("main_source_commit") == BASE.MAIN_COMMIT,
                     "exact-main source identity changed")
        encode_us = metric(result, "encode_execution")
    else:
        target = cell["role"] == "target"
        BASE.require(result.get("schema") == (
                         "leopard2-benchmark-v17" if target
                         else "leopard2-benchmark-v2") and
                     resolved.get("backend") == "avx2" and
                     isinstance(correctness, dict) and
                     correctness.get("leopard2_round_trip") is True,
                     "Leopard2 schema/backend/round trip failed")
        build = result.get("build")
        BASE.require(isinstance(build, dict), "Leopard2 build identity absent")
        if target:
            expected_mode = 1 if implementation == "candidate" else 2
            BASE.require(
                build.get("t16_b64_three_pass_diagnostic_mode") ==
                    expected_mode and
                build.get("t16_b64_three_pass_enabled") is
                    (implementation == "candidate") and
                build.get("t16_b64_timed_q1_encode_api") ==
                    "leo2_encode_batch:item_count=1:no_preflight_scratch" and
                build.get("t16_b64_timed_direct_encode_api") == "leo2_encode",
                "runtime selector or timed API identity differs from label")
        encode_us = metric(result,
                           "one_shot_encode" if cell["metric"] == "direct"
                           else "encode_execution")
    return {
        "encode_us": encode_us,
        "metric": "old_leo_encode" if implementation == "main"
                  else cell["metric"],
        "digests": dict(digests),
        "schema": result["schema"],
    }


BASE.validate_result = validate_result


if __name__ == "__main__":
    raise SystemExit(BASE.main())
