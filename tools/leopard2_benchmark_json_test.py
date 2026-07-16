#!/usr/bin/env python3
"""Regression-check bench_leopard2's default and external-evidence JSON shapes."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run(executable: Path, skip_legacy: bool) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="leo2-benchmark-json-") as temporary:
        output = Path(temporary) / "result.json"
        command = [
            str(executable), "--k", "3", "--r", "2", "--profile", "high",
            "--field", "gf8", "--backend", "auto", "--bytes", "64",
            "--loss", "1", "--batch", "1", "--reuse", "1",
            "--iterations", "1", "--warmup", "0", "--threads", "1",
            "--seed", "7", "--json", str(output),
        ]
        if skip_legacy:
            command.insert(-2, "--skip-legacy")
        completed = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=False)
        require(completed.returncode == 0,
                f"benchmark failed: {completed.stdout}\n{completed.stderr}")
        require(not completed.stdout and not completed.stderr,
                "benchmark emitted output while writing JSON")
        return json.loads(output.read_text())


def validate_common(document: dict[str, Any]) -> None:
    require(set(document) == {
        "schema", "build", "parameters", "resolved", "correctness",
        "memory", "metrics", "legacy"}, "top-level JSON keys changed")
    require(document["schema"] == "leopard2-benchmark-v1", "schema changed")
    require(set(document["build"]) == {
        "compiler", "compiler_version", "cplusplus"}, "build keys changed")
    require(set(document["resolved"]) == {
        "profile", "field", "backend", "thread_count", "parent_count",
        "padded_side"}, "resolved keys changed")
    require(set(document["correctness"]) == {
        "leopard2_round_trip", "legacy_comparison"}, "correctness keys changed")
    require(set(document["memory"]) == {
        "scratch_alignment", "encode_scratch_bytes_per_stripe",
        "decode_scratch_bytes_per_stripe", "encode_scratch_bytes_batch",
        "decode_scratch_bytes_batch"}, "memory keys changed")
    require(set(document["metrics"]) == {
        "codec_setup", "encode_execution", "decode_plan_setup",
        "decode_execution", "decode_amortized_at_reuse", "rate_semantics"},
        "metrics keys changed")
    require(set(document["legacy"]) == {
        "available", "unavailable_reason", "codec_setup",
        "decode_timing_includes_setup", "encode_execution",
        "decode_including_setup"}, "legacy keys changed")
    for metric in ("codec_setup", "decode_plan_setup"):
        require(set(document["metrics"][metric]) == {
            "median_us", "mad_us", "minimum_us", "maximum_us"},
            f"default {metric} unexpectedly retained raw samples")
    for metric in ("encode_execution", "decode_execution"):
        require(not any(key.startswith("samples_us")
                        for key in document["metrics"][metric]),
                f"default {metric} unexpectedly retained raw samples")


def main() -> int:
    if len(sys.argv) != 2:
        raise RuntimeError("usage: leopard2_benchmark_json_test.py BENCH_LEOPARD2")
    executable = Path(sys.argv[1]).resolve()
    default = run(executable, False)
    validate_common(default)
    require(set(default["parameters"]) == {
        "K", "R", "requested_profile", "requested_field",
        "requested_backend", "force_generic_decode",
        "force_specialized_decode", "shard_bytes", "loss_count",
        "missing_original_indices", "batch", "reuse", "iterations",
        "warmup", "thread_count", "seed"},
        "default parameter structure changed")
    require(default["legacy"]["available"] is True and
            default["correctness"]["legacy_comparison"] == "matched",
            "default behavior no longer executes the available legacy oracle")

    external = run(executable, True)
    validate_common(external)
    require(set(external["parameters"]) == set(default["parameters"]) | {"skip_legacy"},
            "external-evidence parameter structure changed")
    require(external["parameters"]["skip_legacy"] is True,
            "external-evidence mode was not recorded")
    require(external["legacy"] == {
        "available": False,
        "unavailable_reason":
            "disabled by --skip-legacy for symmetric external comparison",
        "codec_setup": None,
        "decode_timing_includes_setup": True,
        "encode_execution": None,
        "decode_including_setup": None,
    }, "external-evidence mode did not completely skip legacy work")
    require(external["correctness"]["legacy_comparison"] is None,
            "external-evidence mode claimed a legacy comparison")
    print("leopard2 benchmark JSON regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
