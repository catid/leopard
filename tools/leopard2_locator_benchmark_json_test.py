#!/usr/bin/env python3
"""Validate the locator microbenchmark's machine-readable evidence contract."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def exact_keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    require(set(value) == expected, f"{label} keys changed: {sorted(value)}")


def require_rejected(executable: Path, arguments: list[str], label: str) -> None:
    completed = subprocess.run(
        [str(executable), *arguments], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True, check=False, timeout=5)
    require(completed.returncode == 2, f"{label} was not rejected")


def main() -> int:
    require(len(sys.argv) == 2,
            "usage: leopard2_locator_benchmark_json_test.py EXECUTABLE")
    executable = Path(sys.argv[1]).resolve()
    require_rejected(executable, ["--warmup", "-1"], "negative integer")
    require_rejected(
        executable, ["--calls", "42949672960"], "overflowing integer")
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    command = [
        str(executable), "--field", "gf8", "--n", "256", "--erasures", "8",
        "--calls", "4", "--iterations", "5", "--warmup", "1", "--seed", "7",
    ]
    completed = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        check=False, env=environment)
    require(completed.returncode == 0,
            f"benchmark failed: {completed.stdout}\n{completed.stderr}")
    require(not completed.stderr, f"benchmark stderr: {completed.stderr}")
    document = json.loads(completed.stdout)

    exact_keys(document, {
        "schema", "build", "runtime", "field", "parent_n", "erasure_count",
        "calls_per_sample", "iterations", "warmup_calls", "seed",
        "measurement_order", "dispatcher", "direct", "active_walsh",
    }, "document")
    require(document["schema"] == "leopard2-locator-benchmark-v2",
            "locator benchmark schema changed")
    exact_keys(document["build"], {
        "source_git_sha", "source_dirty_at_configure", "compiler",
        "compiler_version"}, "build")
    exact_keys(document["runtime"], {
        "allowed_cpus", "openmp_macro", "omp_num_threads", "omp_dynamic"},
        "runtime")
    if sys.platform.startswith("linux"):
        require(document["runtime"]["allowed_cpus"],
                "Linux affinity identity is empty")
    require(document["runtime"]["omp_num_threads"] == "1" and
            document["runtime"]["omp_dynamic"] == "FALSE",
            "controlled OpenMP identity missing")
    require(document["measurement_order"] == "ABBA",
            "measurement order changed")
    require(document["iterations"] == 5 and document["calls_per_sample"] == 4,
            "benchmark dimensions changed")
    require(document["dispatcher"] == "active_walsh",
            "production dispatcher boundary changed")
    for name in ("direct", "active_walsh"):
        exact_keys(document[name], {"median_us", "mad_us", "samples_us"}, name)
        samples = document[name]["samples_us"]
        require(isinstance(samples, list) and len(samples) == 5 and
                all(isinstance(value, (int, float)) and value > 0
                    for value in samples),
                f"{name} raw samples invalid")
        require(document[name]["median_us"] == sorted(samples)[2],
                f"{name} median does not match raw samples")
    print("leopard2 locator benchmark JSON test passed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"locator benchmark JSON test failed: {error}", file=sys.stderr)
        raise SystemExit(1)
