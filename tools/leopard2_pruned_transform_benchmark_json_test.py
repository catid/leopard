#!/usr/bin/env python3
"""Smoke and fail-closed schema checks for the C1 transform benchmark."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys


class SmokeError(RuntimeError):
    pass


def run(executable: pathlib.Path, arguments: list[str], expect_success: bool) -> str:
    completed = subprocess.run(
        [str(executable), *arguments],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=60,
    )
    if expect_success and completed.returncode != 0:
        raise SmokeError(
            f"benchmark failed ({completed.returncode}): {completed.stderr.strip()}"
        )
    if not expect_success and completed.returncode == 0:
        raise SmokeError("benchmark accepted a malformed invocation")
    return completed.stdout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SmokeError(message)


def validate_result(raw: str, expected: dict[str, object]) -> None:
    try:
        result = json.loads(raw)
    except json.JSONDecodeError as error:
        raise SmokeError(f"benchmark output is not JSON: {error}") from error

    require(
        result.get("schema") == "leopard2-pruned-transform-benchmark-v1",
        "unexpected schema",
    )
    require(result.get("authoritative") is False, "smoke result claimed authority")
    require(
        isinstance(result.get("authority_note"), str)
        and "isolated runner" in result["authority_note"],
        "missing authority caveat",
    )
    build = result.get("build")
    require(isinstance(build, dict), "missing build identity")
    require(isinstance(build.get("source_git_sha"), str), "missing source SHA")
    require(build.get("source_dirty") in (0, 1), "invalid dirty marker")

    parameters = result.get("parameters")
    require(isinstance(parameters, dict), "missing parameters")
    for name, value in expected.items():
        require(parameters.get(name) == value, f"parameter {name} differs")

    resolved = result.get("resolved")
    require(isinstance(resolved, dict), "missing resolved section")
    require(resolved.get("backend") == "scalar", "wrong runtime backend")
    require(
        isinstance(resolved.get("workspace_instances"), int)
        and 1 <= resolved["workspace_instances"] <= parameters["iterations"],
        "invalid workspace instance count",
    )

    plan = result.get("plan")
    require(isinstance(plan, dict), "missing plan section")
    operations = plan.get("operations")
    fused = plan.get("fused_four_descriptors")
    require(
        isinstance(operations, int)
        and isinstance(fused, int)
        and 0 <= fused * 4 <= operations,
        "invalid operation/fusion counts",
    )
    require(
        plan.get("effective_execution_steps") == operations - fused * 3,
        "effective step count differs",
    )
    require(
        0 < operations <= plan.get("full_butterflies", -1),
        "pruned plan exceeds the full graph",
    )
    require(
        plan.get("flat_storage_bytes", 0) > 0
        and plan.get("fused_storage_bytes", 0) >= plan["flat_storage_bytes"],
        "invalid plan storage accounting",
    )

    correctness = result.get("correctness")
    require(
        isinstance(correctness, dict)
        and correctness.get("matches_full") is True,
        "full-transform comparison failed",
    )
    digest = correctness.get("digest_fnv1a64")
    require(
        isinstance(digest, str)
        and len(digest) == 18
        and digest.startswith("0x"),
        "invalid correctness digest",
    )

    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "missing metrics")
    for name in ("setup_ns", "full_ns", "flat_ns", "fused_ns"):
        metric = metrics.get(name)
        require(isinstance(metric, dict), f"missing {name}")
        samples = metric.get("samples")
        require(
            isinstance(samples, list) and len(samples) == parameters["samples"],
            f"wrong {name} sample count",
        )
        require(
            all(isinstance(value, (int, float)) and value > 0 for value in samples),
            f"invalid {name} samples",
        )
        ordered = sorted(float(value) for value in samples)
        expected_median = ordered[len(ordered) // 2]
        require(
            math.isclose(float(metric.get("median", -1)), expected_median, abs_tol=0.001),
            f"wrong {name} median",
        )
    for name in ("full_over_flat", "full_over_fused"):
        require(
            isinstance(metrics.get(name), (int, float))
            and math.isfinite(metrics[name])
            and metrics[name] > 0,
            f"invalid {name}",
        )


def main() -> int:
    if len(sys.argv) != 2:
        raise SmokeError("usage: leopard2_pruned_transform_benchmark_json_test.py BENCHMARK")
    executable = pathlib.Path(sys.argv[1]).resolve()
    require(executable.is_file(), "benchmark executable is missing")

    common = [
        "--backend", "scalar",
        "--iterations", "4",
        "--samples", "3",
        "--warmups", "1",
        "--setup-iterations", "2",
        "--memory-mib", "16",
        "--size", "16",
        "--input-active", "9",
        "--output-requested", "7",
    ]
    cases = [
        (
            ["--field", "gf8", "--direction", "forward", "--bytes", "65", *common],
            {"field": "gf8", "direction": "forward", "shard_bytes": 65},
        ),
        (
            [
                "--field", "gf16", "--direction", "inverse", "--bytes", "66",
                "--shift", "16", "--input-shape", "holey",
                "--output-shape", "holey", *common,
            ],
            {
                "field": "gf16", "direction": "inverse", "shard_bytes": 66,
                "shift": 16, "input_shape": "holey", "output_shape": "holey",
                "full_input_prefix": 14, "full_output_prefix": 14,
            },
        ),
    ]
    for arguments, expected in cases:
        expected.update(
            size=16,
            input_active=9,
            output_requested=7,
            requested_backend="scalar",
            iterations=4,
            samples=3,
            warmups=1,
            setup_iterations=2,
            memory_mib=16,
        )
        expected.setdefault("input_shape", "prefix")
        expected.setdefault("output_shape", "prefix")
        expected.setdefault("full_input_prefix", 9)
        expected.setdefault("full_output_prefix", 7)
        validate_result(run(executable, arguments, True), expected)

    run(executable, ["--field", "gf16", "--bytes", "65"], False)
    run(executable, ["--size", "3"], False)
    run(executable, ["--input-active", "0"], False)
    run(executable, ["--size", "4294967360"], False)
    run(executable, ["--samples", "4"], False)
    print("PASS pruned transform benchmark JSON smoke cases=2 rejects=5")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (SmokeError, subprocess.TimeoutExpired) as error:
        print(f"FAIL pruned transform benchmark JSON smoke: {error}", file=sys.stderr)
        raise SystemExit(1)
