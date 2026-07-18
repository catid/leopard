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


def run(
    executable: Path,
    external_evidence: bool,
    report_decode_path: bool = False,
) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="leo2-benchmark-json-") as temporary:
        output = Path(temporary) / "result.json"
        command = [
            str(executable), "--k", "3", "--r", "2", "--profile", "high",
            "--field", "auto", "--backend", "auto", "--bytes", "64",
            "--loss", "1", "--batch", "1", "--reuse", "1",
            "--iterations", "1", "--warmup", "0", "--threads", "1",
            "--seed", "7",
        ]
        if external_evidence:
            command.extend(("--skip-legacy", "--retain-samples"))
        if report_decode_path:
            command.append("--report-decode-path")
        command.extend(("--json", str(output)))
        completed = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, check=False)
        require(completed.returncode == 0,
                f"benchmark failed: {completed.stdout}\n{completed.stderr}")
        require(not completed.stdout and not completed.stderr,
                "benchmark emitted output while writing JSON")
        return json.loads(output.read_text())


def validate_common(document: dict[str, Any], retain_samples: bool) -> None:
    expected_top = {
        "schema", "build", "parameters", "resolved", "correctness",
        "memory", "metrics", "legacy"}
    if document["schema"] in {
        "leopard2-benchmark-v2", "leopard2-benchmark-v3"
    }:
        expected_top.add("workload_digests")
    require(set(document) == expected_top, "top-level JSON keys changed")
    require(set(document["build"]) == {
        "compiler", "compiler_version", "cplusplus"}, "build keys changed")
    expected_resolved = {
        "profile", "field", "backend", "thread_count", "parent_count",
        "padded_side"}
    if document["schema"] == "leopard2-benchmark-v3":
        expected_resolved.update({
            "selected_decode_path", "selected_decode_rule",
            "decode_required_work_slots", "decode_aligned_prefix_bytes",
            "decode_tail_bytes", "decode_rounded_bytes",
            "decode_multi_item_batch",
        })
    require(set(document["resolved"]) == expected_resolved,
            "resolved keys changed")
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
        expected = {
            "median_us", "mad_us", "minimum_us", "maximum_us"}
        if retain_samples:
            expected.add("samples_us")
        require(set(document["metrics"][metric]) == expected,
                f"{metric} raw-sample structure changed")
        samples = document["metrics"][metric].get("samples_us")
        require((isinstance(samples, list) and len(samples) == 1) == retain_samples,
                f"{metric} raw-sample cardinality changed")
    for metric in ("encode_execution", "decode_execution"):
        expected = {
            "median_us_per_batch_call", "mad_us_per_batch_call",
            "minimum_us_per_batch_call", "maximum_us_per_batch_call",
            ("input_GB_per_s" if metric == "encode_execution" else
             "offered_received_GB_per_s"),
            ("parity_output_GB_per_s" if metric == "encode_execution" else
             "repaired_output_GB_per_s")}
        if retain_samples:
            expected.add("samples_us_per_batch_call")
        require(set(document["metrics"][metric]) == expected,
                f"{metric} structure changed")
        samples = document["metrics"][metric].get("samples_us_per_batch_call")
        require((isinstance(samples, list) and len(samples) == 1) == retain_samples,
                f"{metric} raw-sample cardinality changed")


def validate_workload_digests(document: dict[str, Any]) -> None:
    digests = document.get("workload_digests")
    require(isinstance(digests, dict) and set(digests) == {
        "algorithm", "original_data", "transmitted_parity",
        "recovered_originals"}, "workload digest structure changed")
    require(digests["algorithm"] == "fnv1a64", "workload digest algorithm changed")
    for name in ("original_data", "transmitted_parity", "recovered_originals"):
        value = digests[name]
        require(isinstance(value, str) and len(value) == 16 and
                all(character in "0123456789abcdef" for character in value),
                f"workload digest {name} is not lowercase FNV-1a hex")


def validate_isal_comparison_contract(document: dict[str, Any]) -> None:
    # Exercise the exact parser used by future ISA-L collection, rather than
    # letting this executable-shape regression and the retained-artifact
    # validator drift independently.
    import leopard2_isal_compare as comparison

    cell = {
        "K": 3, "R": 2, "profile": "high", "shard_bytes": 64,
        "loss_count": 1, "batch": 1, "reuse": 1, "seed": 7,
    }
    comparison.validate_leopard_result(document, cell, 1, 0)


def main() -> int:
    if len(sys.argv) != 2:
        raise RuntimeError("usage: leopard2_benchmark_json_test.py BENCH_LEOPARD2")
    executable = Path(sys.argv[1]).resolve()
    default = run(executable, False)
    require(default["schema"] == "leopard2-benchmark-v1",
            "default benchmark schema changed")
    validate_common(default, False)
    require(set(default["parameters"]) == {
        "K", "R", "requested_profile", "requested_field",
        "requested_backend", "force_generic_decode",
        "force_specialized_decode", "force_tiled_decode",
        "force_materialized_decode", "shard_bytes", "loss_count",
        "missing_original_indices", "batch", "reuse", "iterations",
        "warmup", "thread_count", "seed"},
        "default parameter structure changed")
    require(default["parameters"]["force_tiled_decode"] is False and
            default["parameters"]["force_materialized_decode"] is False,
            "default benchmark forced a specialized workspace kernel")
    require(default["legacy"]["available"] is True and
            default["correctness"]["legacy_comparison"] == "matched",
            "default behavior no longer executes the available legacy oracle")

    external = run(executable, True)
    require(external["schema"] in {
        "leopard2-benchmark-v1", "leopard2-benchmark-v2"},
        "external-evidence benchmark schema changed")
    validate_common(external, True)
    expected_external_parameters = set(default["parameters"]) | {"skip_legacy"}
    if external["schema"] == "leopard2-benchmark-v2":
        expected_external_parameters.add("retain_samples")
        validate_workload_digests(external)
    require(set(external["parameters"]) == expected_external_parameters,
            "external-evidence parameter structure changed")
    require(external["parameters"]["skip_legacy"] is True,
            "external-evidence mode was not recorded")
    if external["schema"] == "leopard2-benchmark-v2":
        require(external["parameters"]["retain_samples"] is True,
                "external-evidence raw-sample mode was not recorded")
    unavailable_reason = (
        "disabled by --skip-legacy" if
        external["schema"] == "leopard2-benchmark-v2" else
        "disabled by --skip-legacy for symmetric external comparison")
    require(external["legacy"] == {
        "available": False,
        "unavailable_reason": unavailable_reason,
        "codec_setup": None,
        "decode_timing_includes_setup": True,
        "encode_execution": None,
        "decode_including_setup": None,
    }, "external-evidence mode did not completely skip legacy work")
    require(external["correctness"]["legacy_comparison"] is None,
            "external-evidence mode claimed a legacy comparison")
    validate_isal_comparison_contract(external)

    path_report = run(executable, True, True)
    require(path_report["schema"] == "leopard2-benchmark-v3",
            "path-report benchmark schema changed")
    validate_common(path_report, True)
    validate_workload_digests(path_report)
    require(set(path_report["parameters"]) ==
            (expected_external_parameters | {"report_decode_path"}),
            "path-report parameter structure changed")
    require(path_report["parameters"]["report_decode_path"] is True,
            "path-report opt-in was not recorded")
    resolved = path_report["resolved"]
    require(resolved["selected_decode_path"] == "direct" and
            resolved["selected_decode_rule"] == "direct" and
            resolved["decode_required_work_slots"] == 0,
            "path-report did not observe the production direct repair")
    require(resolved["decode_aligned_prefix_bytes"] == 64 and
            resolved["decode_tail_bytes"] == 0 and
            resolved["decode_rounded_bytes"] == 64 and
            resolved["decode_multi_item_batch"] is False,
            "path-report byte/batch geometry differs")
    print("leopard2 benchmark JSON regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
