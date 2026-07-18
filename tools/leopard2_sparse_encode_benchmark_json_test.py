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
    median = metric.get("median")
    require(
        isinstance(median, (int, float))
        and math.isclose(float(median), statistics.median(values), abs_tol=0.001),
        f"wrong median for {name}",
    )
    require(metric.get("minimum") == min(values), f"wrong minimum for {name}")
    require(metric.get("maximum") == max(values), f"wrong maximum for {name}")


def validate(raw: str, expected: dict[str, object]) -> None:
    try:
        result = json.loads(raw)
    except json.JSONDecodeError as error:
        raise SmokeError(f"benchmark output is not JSON: {error}") from error
    require(
        result.get("schema") == "leopard2-sparse-encode-benchmark-v2",
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
        isinstance(build.get("library_test_hooks"), bool),
        "missing instrumentation marker",
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
    require(
        plan.get("inverse_candidate_available") is True
        and isinstance(plan.get("inverse_operations"), int)
        and plan["inverse_operations"] > 0
        and isinstance(plan.get("inverse_source_groups"), int)
        and isinstance(plan.get("inverse_active_prefix"), int)
        and 0 < plan["inverse_active_prefix"] < resolved["padded_side"],
        "strict-prefix smoke cell did not compile an inverse source plan",
    )
    require(
        plan["inverse_source_groups"] == plan["inverse_active_prefix"] // 4,
        "inverse source plan has wrong complete four-source group count",
    )
    require(
        isinstance(plan.get("mature_zero_fill_shards"), int)
        and plan["mature_zero_fill_shards"] > 0
        and plan.get("pruned_zero_fill_shards") == 0,
        "invalid inverse-prefix zero-fill accounting",
    )

    correctness = result.get("correctness")
    require(
        isinstance(correctness, dict)
        and correctness.get("exact_prefix_parity_match") is True,
        "parity comparison failed",
    )
    require(
        correctness.get("inverse_pruned_parity_match") is True,
        "inverse-pruned parity comparison failed",
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
        "inverse_schedule_setup_ns",
        "prefix_execution_ns",
        "exact_prepared_execution_ns",
        "exact_call_local_total_ns",
        "prefix_pruned_inverse_execution_ns",
        "exact_pruned_inverse_execution_ns",
    )
    for name in names:
        validate_metric(metrics.get(name), int(parameters["samples"]), name)
    for name in (
        "prefix_over_exact_prepared",
        "prefix_over_exact_call_local",
        "mature_over_pruned_inverse_prefix",
        "mature_over_pruned_inverse_exact",
    ):
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
    prefix_inverse = float(
        metrics["prefix_pruned_inverse_execution_ns"]["median"]
    )
    exact_inverse = float(
        metrics["exact_pruned_inverse_execution_ns"]["median"]
    )
    require(
        math.isclose(
            float(metrics["mature_over_pruned_inverse_prefix"]),
            mature / prefix_inverse,
            abs_tol=0.002,
        ),
        "wrong inverse-prefix ratio",
    )
    require(
        math.isclose(
            float(metrics["mature_over_pruned_inverse_exact"]),
            execution / exact_inverse,
            abs_tol=0.002,
        ),
        "wrong inverse-exact ratio",
    )
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
            ["--profile", "high", "--field", "gf8", "--k", "49", "--r", "16",
             "--requested-parity", "0,7,15", *common],
            {"profile": "high", "field": "gf8", "K": 49, "R": 16,
             "requested_parity": [0, 7, 15]},
        ),
        (
            ["--profile", "low", "--field", "gf8", "--k", "17", "--r", "48",
             "--requested-parity", "0,15,16,31,47", *common],
            {"profile": "low", "field": "gf8", "K": 17, "R": 48,
             "requested_parity": [0, 15, 16, 31, 47]},
        ),
        (
            ["--profile", "high", "--field", "gf16", "--k", "49", "--r", "16",
             "--requested-parity", "0,7,15", *common],
            {"profile": "high", "field": "gf16", "K": 49, "R": 16,
             "requested_parity": [0, 7, 15]},
        ),
        (
            ["--profile", "low", "--field", "gf16", "--k", "17", "--r", "48",
             "--requested-parity", "0,15,16,31,47", *common],
            {"profile": "low", "field": "gf16", "K": 17, "R": 48,
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
