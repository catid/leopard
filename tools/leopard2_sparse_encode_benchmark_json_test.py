#!/usr/bin/env python3
"""Smoke and fail-closed schema checks for sparse encoder crossover cells."""

from __future__ import annotations

import json
import math
import os
import pathlib
import statistics
import subprocess
import sys


class SmokeError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SmokeError(message)


def run(executable: pathlib.Path, arguments: list[str], success: bool) -> str:
    environment = os.environ.copy()
    environment.update(OMP_DYNAMIC="FALSE", OMP_NUM_THREADS="1")
    completed = subprocess.run(
        [str(executable), *arguments],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=60,
        env=environment,
    )
    if success and completed.returncode != 0:
        raise SmokeError(
            f"benchmark failed ({completed.returncode}): {completed.stderr.strip()}"
        )
    if not success and completed.returncode == 0:
        raise SmokeError("benchmark accepted a malformed invocation")
    return completed.stdout


def validate_metric(metric: object, samples: int, name: str) -> None:
    require(isinstance(metric, dict), f"missing metric {name}")
    values = metric.get("samples")
    require(
        isinstance(values, list) and len(values) == samples,
        f"wrong sample count for {name}",
    )
    require(
        all(isinstance(value, (int, float)) and value > 0 for value in values),
        f"invalid samples for {name}",
    )
    expected_median = statistics.median(values)
    median = metric.get("median")
    require(
        isinstance(median, (int, float))
        and math.isclose(float(median), expected_median, abs_tol=0.002),
        f"wrong median for {name}",
    )
    expected_mad = statistics.median(
        abs(float(value) - expected_median) for value in values
    )
    mad = metric.get("mad")
    require(
        isinstance(mad, (int, float))
        and math.isclose(float(mad), expected_mad, abs_tol=0.002),
        f"wrong MAD for {name}",
    )
    require(metric.get("minimum") == min(values), f"wrong minimum for {name}")
    require(metric.get("maximum") == max(values), f"wrong maximum for {name}")


def validate(raw: str, expected: dict[str, object]) -> None:
    try:
        result = json.loads(raw)
    except json.JSONDecodeError as error:
        raise SmokeError(f"benchmark output is not JSON: {error}") from error
    require(
        result.get("schema") == "leopard2-sparse-encode-benchmark-v3",
        "unexpected schema",
    )
    require(result.get("authoritative") is False, "cell claimed authority")
    require(
        "fail-closed pinned runner" in str(result.get("authority_note")),
        "missing authority caveat",
    )
    build = result.get("build")
    require(isinstance(build, dict), "missing build identity")
    require(
        isinstance(build.get("source_git_sha"), str)
        and len(build["source_git_sha"]) >= 7,
        "missing source SHA",
    )
    require(build.get("source_dirty") in (0, 1), "invalid dirty marker")
    require(
        build.get("library_test_hooks") is False,
        "benchmark linked test-hook instrumentation",
    )
    require(
        isinstance(build.get("compiler"), str) and build["compiler"],
        "missing compiler identity",
    )
    require(
        isinstance(build.get("compiler_version"), str) and build["compiler_version"],
        "missing compiler version",
    )
    require(isinstance(build.get("cplusplus"), int), "missing C++ language identity")

    parameters = result.get("parameters")
    require(isinstance(parameters, dict), "missing parameters")
    for name, value in expected.items():
        require(parameters.get(name) == value, f"parameter {name} differs")
    require(parameters.get("reuse") == [1, 8, 64], "reuse points differ")

    resolved = result.get("resolved")
    require(isinstance(resolved, dict), "missing resolved section")
    require(resolved.get("backend") == "scalar", "wrong runtime backend")
    require(
        isinstance(resolved.get("padded_side"), int)
        and resolved["padded_side"] >= 2,
        "invalid padded side",
    )

    plan = result.get("plan")
    require(isinstance(plan, dict), "missing plan accounting")
    prefix = plan.get("prefix_butterflies")
    retained = plan.get("retained_butterflies")
    require(
        isinstance(prefix, int)
        and isinstance(retained, int)
        and 0 < retained < prefix <= plan.get("full_butterflies", -1),
        "sparse test mask did not reduce the mature prefix graph",
    )
    require(
        isinstance(plan.get("schedule_bytes"), int)
        and plan["schedule_bytes"] > 0
        and isinstance(plan.get("dependency_workspace_bytes"), int)
        and plan["dependency_workspace_bytes"] > 0,
        "invalid schedule storage",
    )

    correctness = result.get("correctness")
    require(
        isinstance(correctness, dict)
        and correctness.get("exact_prefix_parity_match") is True,
        "parity comparison failed",
    )
    digest = correctness.get("digest_fnv1a64")
    require(
        isinstance(digest, str) and len(digest) == 18 and digest.startswith("0x"),
        "invalid parity digest",
    )

    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "missing metrics")
    names = (
        "schedule_setup_ns",
        "prefix_execution_ns",
        "exact_prepared_execution_ns",
        "exact_call_local_total_ns",
    )
    for name in names:
        validate_metric(metrics.get(name), int(parameters["samples"]), name)
    for name in ("prefix_over_exact_prepared", "prefix_over_exact_call_local"):
        require(
            isinstance(metrics.get(name), (int, float))
            and math.isfinite(metrics[name])
            and metrics[name] > 0,
            f"invalid ratio {name}",
        )
    amortized = metrics.get("amortized_exact")
    require(
        isinstance(amortized, list)
        and [row.get("reuse") for row in amortized] == [1, 8, 64],
        "invalid amortization rows",
    )
    setup = float(metrics["schedule_setup_ns"]["median"])
    execution = float(metrics["exact_prepared_execution_ns"]["median"])
    mature = float(metrics["prefix_execution_ns"]["median"])
    for row in amortized:
        modeled = execution + setup / row["reuse"]
        require(
            math.isclose(float(row.get("modeled_ns", -1)), modeled, abs_tol=0.002),
            "wrong amortized exact model",
        )
        require(
            math.isclose(
                float(row.get("prefix_over_modeled_exact", -1)),
                mature / modeled,
                abs_tol=0.002,
            ),
            "wrong amortized ratio",
        )


def main() -> int:
    if len(sys.argv) != 2:
        raise SmokeError("usage: leopard2_sparse_encode_benchmark_json_test.py BENCHMARK")
    executable = pathlib.Path(sys.argv[1]).resolve()
    require(executable.is_file(), "benchmark executable is missing")
    common = [
        "--backend", "scalar",
        "--bytes", "64",
        "--iterations", "2",
        "--samples", "3",
        "--warmups", "1",
        "--setup-iterations", "2",
        "--reuse", "1,8,64",
        "--memory-mib", "16",
    ]
    cases = (
        (
            ["--profile", "high", "--field", "gf8", "--k", "48", "--r", "16",
             "--requested-parity", "0,7,15", *common],
            {"profile": "high", "field": "gf8", "K": 48, "R": 16,
             "requested_parity": [0, 7, 15]},
        ),
        (
            ["--profile", "low", "--field", "gf8", "--k", "16", "--r", "48",
             "--requested-parity", "0,15,16,31,47", *common],
            {"profile": "low", "field": "gf8", "K": 16, "R": 48,
             "requested_parity": [0, 15, 16, 31, 47]},
        ),
        (
            ["--profile", "high", "--field", "gf16", "--k", "48", "--r", "16",
             "--requested-parity", "0,7,15", *common],
            {"profile": "high", "field": "gf16", "K": 48, "R": 16,
             "requested_parity": [0, 7, 15]},
        ),
        (
            ["--profile", "low", "--field", "gf16", "--k", "16", "--r", "48",
             "--requested-parity", "0,15,16,31,47", *common],
            {"profile": "low", "field": "gf16", "K": 16, "R": 48,
             "requested_parity": [0, 15, 16, 31, 47]},
        ),
    )
    for arguments, expected in cases:
        expected.update(
            shard_bytes=64,
            requested_backend="scalar",
            iterations=2,
            samples=3,
            warmups=1,
            setup_iterations=2,
            memory_mib=16,
            seed=0x535041525345454E,
        )
        validate(run(executable, arguments, True), expected)

    rejects = (
        ["--field", "gf16", "--bytes", "65"],
        ["--profile", "low", "--field", "gf8", "--k", "129", "--r", "1"],
        ["--requested-parity", ""],
        ["--r", "8", "--requested-parity", "8"],
        ["--samples", "4"],
        ["--reuse", "1,1"],
        ["--memory-mib", "0"],
    )
    for arguments in rejects:
        run(executable, arguments, False)
    print("PASS sparse encode benchmark JSON smoke cases=4 rejects=7")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (SmokeError, subprocess.TimeoutExpired) as error:
        print(f"FAIL sparse encode benchmark JSON smoke: {error}", file=sys.stderr)
        raise SystemExit(1)
