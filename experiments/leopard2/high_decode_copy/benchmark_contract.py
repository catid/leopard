#!/usr/bin/env python3
"""Validate the private Algorithm 5 copy/no-copy attribution benchmark."""

from __future__ import annotations

import argparse
import copy
import json
import math
import re
import statistics
import subprocess
import tempfile
from pathlib import Path
from typing import Any


HEX64 = re.compile(r"^[0-9a-f]{16}$")


class ContractError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def finite_number(value: object, name: str, *, positive: bool = True) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{name} is not numeric")
    result = float(value)
    require(math.isfinite(result) and (result > 0 if positive else result >= 0),
            f"{name} is not finite and bounded")
    return result


def close_enough(actual: float, expected: float) -> bool:
    return abs(actual - expected) <= max(0.000002, abs(expected) * 0.000002)


def validate_summary(summary: object, iterations: int, *, setup: bool) -> float:
    require(isinstance(summary, dict), "timing summary is not an object")
    suffix = "" if setup else "_per_batch_call"
    sample_name = "samples_us" if setup else "samples_us_per_batch_call"
    timing_keys = {
        f"median_us{suffix}", f"mad_us{suffix}", f"minimum_us{suffix}",
        f"maximum_us{suffix}", sample_name}
    if setup:
        require(set(summary) == timing_keys,
                "setup timing summary contains an unbound claim")
    samples_value = summary.get(sample_name)
    require(type(samples_value) is list and len(samples_value) == iterations,
            "timing summary has the wrong retained sample count")
    samples = [finite_number(value, sample_name) for value in samples_value]
    median = statistics.median(samples)
    expected = {
        f"median_us{suffix}": median,
        f"mad_us{suffix}": statistics.median(
            abs(value - median) for value in samples),
        f"minimum_us{suffix}": min(samples),
        f"maximum_us{suffix}": max(samples),
    }
    for name, derived in expected.items():
        actual = finite_number(summary.get(name), name, positive="mad_us" not in name)
        require(close_enough(actual, derived),
                f"{name} is not derived from retained samples")
    return median


def validate_rate(value: object, byte_count: int, median_us: float, name: str) -> None:
    actual = finite_number(value, name)
    expected = byte_count / (median_us * 1000.0)
    require(close_enough(actual, expected), f"{name} is not derived from its median")


def validate_metrics(result: dict[str, Any], parameters: dict[str, Any]) -> None:
    names = ("K", "R", "shard_bytes", "loss_count", "batch", "reuse",
             "iterations", "warmup", "thread_count", "seed")
    require(all(type(parameters.get(name)) is int for name in names) and
            parameters["K"] > 0 and parameters["R"] > 0 and
            parameters["shard_bytes"] > 0 and
            0 < parameters["loss_count"] <= min(parameters["K"], parameters["R"]) and
            parameters["batch"] == 1 and parameters["reuse"] >= 1 and
            parameters["iterations"] >= 1 and parameters["warmup"] >= 0 and
            parameters["thread_count"] == 1 and parameters["seed"] >= 0,
            "benchmark numeric controls are not exact bounded integers")
    missing = parameters.get("missing_original_indices")
    require(type(missing) is list and len(missing) == parameters["loss_count"] and
            all(type(index) is int and 0 <= index < parameters["K"]
                for index in missing) and missing == sorted(set(missing)),
            "benchmark loss set is not a sorted unique integer array")
    metrics = result.get("metrics")
    require(isinstance(metrics, dict) and set(metrics) == {
                "codec_setup", "encode_execution", "decode_plan_setup",
                "decode_execution", "decode_amortized_at_reuse", "rate_semantics"},
            "benchmark timing metrics are incomplete")
    require(set(metrics["encode_execution"]) == {
                "median_us_per_batch_call", "mad_us_per_batch_call",
                "minimum_us_per_batch_call", "maximum_us_per_batch_call",
                "samples_us_per_batch_call", "input_GB_per_s",
                "parity_output_GB_per_s"} and
            set(metrics["decode_execution"]) == {
                "median_us_per_batch_call", "mad_us_per_batch_call",
                "minimum_us_per_batch_call", "maximum_us_per_batch_call",
                "samples_us_per_batch_call", "offered_received_GB_per_s",
                "repaired_output_GB_per_s"},
            "execution timing summary contains an unbound claim")
    iterations = parameters["iterations"]
    validate_summary(metrics["codec_setup"], iterations, setup=True)
    encode_median = validate_summary(
        metrics["encode_execution"], iterations, setup=False)
    plan_median = validate_summary(metrics["decode_plan_setup"], iterations, setup=True)
    decode_median = validate_summary(
        metrics["decode_execution"], iterations, setup=False)
    k = parameters["K"]
    r = parameters["R"]
    shard_bytes = parameters["shard_bytes"]
    losses = parameters["loss_count"]
    validate_rate(metrics["encode_execution"].get("input_GB_per_s"),
                  k * shard_bytes, encode_median, "encode input rate")
    validate_rate(metrics["encode_execution"].get("parity_output_GB_per_s"),
                  r * shard_bytes, encode_median, "encode parity rate")
    decode_input = (k - losses + r) * shard_bytes
    decode_output = losses * shard_bytes
    validate_rate(metrics["decode_execution"].get("offered_received_GB_per_s"),
                  decode_input, decode_median, "decode input rate")
    validate_rate(metrics["decode_execution"].get("repaired_output_GB_per_s"),
                  decode_output, decode_median, "decode output rate")
    amortized = metrics["decode_amortized_at_reuse"]
    require(isinstance(amortized, dict) and set(amortized) == {
                "reuse_count", "derived_median_us_per_batch_call",
                "offered_received_GB_per_s", "repaired_output_GB_per_s"} and
            type(amortized.get("reuse_count")) is int and
            amortized["reuse_count"] == parameters["reuse"],
            "amortized decode summary has the wrong shape/reuse")
    derived = decode_median + plan_median / parameters["reuse"]
    actual = finite_number(
        amortized.get("derived_median_us_per_batch_call"), "amortized median")
    require(close_enough(actual, derived),
            "amortized median is not plan/reuse plus execution")
    validate_rate(amortized.get("offered_received_GB_per_s"),
                  decode_input, derived, "amortized input rate")
    validate_rate(amortized.get("repaired_output_GB_per_s"),
                  decode_output, derived, "amortized output rate")
    require(metrics["rate_semantics"] ==
            "offered_received counts all non-null shard pointers supplied; "
            "a plan may read a deterministic subset",
            "rate semantics changed")


def validate_document(
    document: object,
    *,
    mode: str,
    workspace: str,
    field: str,
    tail_bytes: int,
    evaluator: str | None = None,
) -> dict[str, Any]:
    require(type(tail_bytes) is int and 0 <= tail_bytes < 64,
            "expected decode tail is not an exact bounded integer")
    require(isinstance(document, dict), "benchmark result is not an object")
    result = document
    require(result.get("schema") == "leopard2-benchmark-v4",
            "attribution benchmark schema changed")
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    attribution = result.get("high_evaluator_attribution")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(isinstance(parameters, dict) and isinstance(resolved, dict) and
            isinstance(attribution, dict) and isinstance(correctness, dict) and
            isinstance(digests, dict), "benchmark attestation is incomplete")
    require(parameters.get("high_evaluator_mode") == mode and
            parameters.get("requested_profile") == "legacy_high_v1" and
            parameters.get("requested_field") == field and
            parameters.get("force_specialized_decode") is True and
            parameters.get("force_generic_decode") is False and
            parameters.get("force_materialized_decode") is
                (workspace == "materialized") and
            parameters.get("force_tiled_decode") is (workspace == "tiled") and
            parameters.get("skip_legacy") is True and
            parameters.get("retain_samples") is True and
            parameters.get("report_decode_path") is True,
            "benchmark request attestation differs from the signed role")
    resolved_tail = resolved.get("decode_tail_bytes")
    require(type(resolved_tail) is int and 0 <= resolved_tail < 64,
            "resolved decode tail is not an exact bounded integer")
    require(resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == field and
            resolved.get("selected_decode_path") == workspace and
            resolved.get("selected_decode_rule") == "forced_" + workspace and
            resolved.get("high_evaluator_mode") == mode and
            resolved_tail == tail_bytes,
            "resolved field/path/tail/mode attestation differs")
    require(correctness.get("leopard2_round_trip") is True and
            correctness.get("legacy_comparison") is None,
            "benchmark did not attest a private round trip")
    require(digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and
                HEX64.fullmatch(digests[name]) is not None
                for name in ("original_data", "transmitted_parity",
                             "recovered_originals")),
            "workload digest attestation is malformed")
    output_blocks = attribution.get("output_blocks")
    two_way = attribution.get("fft_butterfly2_out_of_place")
    four_way = attribution.get("fft_butterfly4_out_of_place")
    fallbacks = attribution.get("compatibility_copy_fallbacks")
    pruned_blocks = attribution.get("pruned_output_blocks")
    mature_blocks = attribution.get("mature_output_blocks")
    require(attribution.get("mode") == mode and
            attribution.get("invariant_passed") is True and
            all(type(value) is int and value >= 0 for value in
                (output_blocks, two_way, four_way, fallbacks,
                 pruned_blocks, mature_blocks)) and
            output_blocks > 0 and
            pruned_blocks + mature_blocks == output_blocks,
            "high-evaluator counter attestation is malformed")
    if mode == "no-copy":
        require(fallbacks == 0 and two_way + four_way > 0,
                "no-copy role entered a copy fallback or no out-of-place layer")
    elif mode == "copy-fallback":
        require(fallbacks == output_blocks and two_way + four_way == 0,
                "copy role did not force exactly one fallback per output block")
    else:
        raise ContractError("unknown signed high-evaluator mode")
    if evaluator is not None:
        require(evaluator == "mature" and pruned_blocks == 0 and
                mature_blocks == output_blocks,
                "full-block role did not exclusively use the mature evaluator")
    validate_metrics(result, parameters)
    return result


def validate_pair(
    no_copy: object,
    copy_fallback: object,
    *,
    workspace: str,
    field: str,
    tail_bytes: int,
    evaluator: str | None = None,
) -> None:
    first = validate_document(
        no_copy, mode="no-copy", workspace=workspace, field=field,
        tail_bytes=tail_bytes, evaluator=evaluator)
    second = validate_document(
        copy_fallback, mode="copy-fallback", workspace=workspace, field=field,
        tail_bytes=tail_bytes, evaluator=evaluator)
    require(first["workload_digests"] == second["workload_digests"] and
            first["parameters"]["missing_original_indices"] ==
                second["parameters"]["missing_original_indices"] and
            first["memory"] == second["memory"],
            "copy/no-copy roles do not describe one identical workload")
    first_parameters = dict(first["parameters"])
    second_parameters = dict(second["parameters"])
    first_parameters.pop("high_evaluator_mode")
    second_parameters.pop("high_evaluator_mode")
    require(first_parameters == second_parameters,
            "copy/no-copy parameters differ outside the mode selector")
    first_resolved = dict(first["resolved"])
    second_resolved = dict(second["resolved"])
    first_resolved.pop("high_evaluator_mode")
    second_resolved.pop("high_evaluator_mode")
    require(first_resolved == second_resolved,
            "copy/no-copy resolution differs outside the mode selector")


def benchmark_command(
    executable: Path,
    output: Path,
    *,
    k: int,
    r: int,
    shard_bytes: int,
    losses: int,
    seed: int,
    field: str,
    workspace: str,
    mode: str,
) -> list[str]:
    return [
        str(executable), "--k", str(k), "--r", str(r),
        "--profile", "high", "--field", field, "--backend", "scalar",
        "--bytes", str(shard_bytes), "--loss", str(losses), "--batch", "1",
        "--reuse", "1", "--iterations", "1", "--warmup", "0",
        "--threads", "1", "--seed", str(seed), "--force-specialized",
        "--force-" + workspace, "--skip-legacy", "--retain-samples",
        "--report-decode-path", "--high-evaluator-mode", mode,
        "--json", str(output),
    ]


def run_checked(arguments: list[str], *, expect_success: bool = True) -> None:
    completed = subprocess.run(
        arguments, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=30, check=False,
        env={
            "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
            "OMP_NUM_THREADS": "1", "OMP_PLACES": "cores",
            "OMP_PROC_BIND": "TRUE", "PATH": "/usr/bin:/bin", "TZ": "UTC",
        })
    if expect_success:
        require(completed.returncode == 0 and not completed.stdout and
                not completed.stderr,
                "attribution benchmark failed or emitted an unexpected stream")
    else:
        require(completed.returncode != 0,
                "attribution benchmark accepted a malformed selector contract")


def smoke(executable: Path) -> None:
    executable = executable.resolve(strict=True)
    cases = (
        ("gf8", "materialized", 192, 64, 64, 4, 101, 0, None),
        ("gf8", "materialized", 192, 64, 65, 4, 103, 1, None),
        ("gf8", "tiled", 240, 16, 256, 8, 107, 0, None),
        ("gf8", "tiled", 240, 16, 257, 8, 109, 1, None),
        ("gf16", "materialized", 257, 63, 64, 4, 113, 0, None),
        ("gf16", "materialized", 257, 63, 66, 4, 127, 2, None),
        ("gf16", "tiled", 1000, 200, 256, 8, 131, 0, None),
        ("gf16", "tiled", 1000, 200, 258, 8, 137, 2, None),
        ("gf8", "materialized", 128, 128, 64, 128, 139, 0, "mature"),
        ("gf8", "tiled", 128, 128, 64, 128, 149, 0, "mature"),
        ("gf16", "materialized", 256, 256, 64, 256, 151, 0, "mature"),
        ("gf16", "tiled", 256, 256, 64, 256, 157, 0, "mature"),
    )
    with tempfile.TemporaryDirectory(prefix="leo2-high-copy-contract-") as root_text:
        root = Path(root_text)
        for index, (field, workspace, k, r, byte_count, losses, seed,
                    tail_bytes, evaluator) in enumerate(cases):
            documents: dict[str, object] = {}
            for mode in ("no-copy", "copy-fallback"):
                output = root / f"{index}-{mode}.json"
                run_checked(benchmark_command(
                    executable, output, k=k, r=r, shard_bytes=byte_count,
                    losses=losses, seed=seed, field=field,
                    workspace=workspace, mode=mode))
                documents[mode] = json.loads(output.read_text(encoding="utf-8"))
            validate_pair(
                documents["no-copy"], documents["copy-fallback"],
                workspace=workspace, field=field, tail_bytes=tail_bytes,
                evaluator=evaluator)

        base = benchmark_command(
            executable, root / "invalid.json", k=192, r=64, shard_bytes=64,
            losses=4, seed=149, field="gf8", workspace="materialized",
            mode="no-copy")
        missing_mode = base[:]
        mode_index = missing_mode.index("--high-evaluator-mode")
        del missing_mode[mode_index:mode_index + 2]
        run_checked(missing_mode, expect_success=False)
        invalid_mode = base[:]
        invalid_mode[invalid_mode.index("no-copy")] = "unknown"
        run_checked(invalid_mode, expect_success=False)
        missing_path = [item for item in base if item != "--force-materialized"]
        run_checked(missing_path, expect_success=False)
    print("high-decode copy attribution benchmark contract passed: 12 paired cells")


def synthetic_document(mode: str) -> dict[str, Any]:
    copy_mode = mode == "copy-fallback"
    document = {
        "schema": "leopard2-benchmark-v4",
        "parameters": {
            "requested_profile": "legacy_high_v1", "requested_field": "gf8",
            "force_specialized_decode": True, "force_generic_decode": False,
            "force_materialized_decode": True, "force_tiled_decode": False,
            "skip_legacy": True, "retain_samples": True,
            "report_decode_path": True, "high_evaluator_mode": mode,
            "K": 8, "R": 8, "shard_bytes": 1, "loss_count": 1,
            "missing_original_indices": [7],
            "batch": 1, "reuse": 1, "iterations": 1, "warmup": 0,
            "thread_count": 1, "seed": 1,
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "selected_decode_path": "materialized",
            "selected_decode_rule": "forced_materialized",
            "decode_tail_bytes": 1, "high_evaluator_mode": mode,
        },
        "high_evaluator_attribution": {
            "mode": mode, "output_blocks": 2,
            "fft_butterfly2_out_of_place": 0 if copy_mode else 16,
            "fft_butterfly4_out_of_place": 0,
            "compatibility_copy_fallbacks": 2 if copy_mode else 0,
            "pruned_output_blocks": 0,
            "mature_output_blocks": 2,
            "invariant_passed": True,
        },
        "correctness": {"leopard2_round_trip": True,
                        "legacy_comparison": None},
        "workload_digests": {
            "algorithm": "fnv1a64", "original_data": "0" * 16,
            "transmitted_parity": "1" * 16,
            "recovered_originals": "2" * 16,
        },
        "memory": {"decode_scratch_bytes_per_stripe": 64},
    }
    document["metrics"] = {
        "codec_setup": {"median_us": 2.0, "mad_us": 0.0,
                        "minimum_us": 2.0, "maximum_us": 2.0,
                        "samples_us": [2.0]},
        "encode_execution": {
            "median_us_per_batch_call": 1.0, "mad_us_per_batch_call": 0.0,
            "minimum_us_per_batch_call": 1.0, "maximum_us_per_batch_call": 1.0,
            "samples_us_per_batch_call": [1.0],
            "input_GB_per_s": 0.008, "parity_output_GB_per_s": 0.008},
        "decode_plan_setup": {"median_us": 3.0, "mad_us": 0.0,
                              "minimum_us": 3.0, "maximum_us": 3.0,
                              "samples_us": [3.0]},
        "decode_execution": {
            "median_us_per_batch_call": 4.0, "mad_us_per_batch_call": 0.0,
            "minimum_us_per_batch_call": 4.0, "maximum_us_per_batch_call": 4.0,
            "samples_us_per_batch_call": [4.0],
            "offered_received_GB_per_s": 0.00375,
            "repaired_output_GB_per_s": 0.00025},
        "decode_amortized_at_reuse": {
            "reuse_count": 1, "derived_median_us_per_batch_call": 7.0,
            "offered_received_GB_per_s": 15.0 / 7000.0,
            "repaired_output_GB_per_s": 1.0 / 7000.0},
        "rate_semantics":
            "offered_received counts all non-null shard pointers supplied; "
            "a plan may read a deterministic subset",
    }
    return document


def self_test() -> None:
    no_copy = synthetic_document("no-copy")
    fallback = synthetic_document("copy-fallback")
    validate_pair(no_copy, fallback, workspace="materialized", field="gf8",
                  tail_bytes=1, evaluator="mature")
    mutations = []
    wrong_mode = copy.deepcopy(no_copy)
    wrong_mode["resolved"]["high_evaluator_mode"] = "copy-fallback"
    mutations.append(wrong_mode)
    wrong_counter = copy.deepcopy(no_copy)
    wrong_counter["high_evaluator_attribution"]["compatibility_copy_fallbacks"] = 1
    mutations.append(wrong_counter)
    wrong_path = copy.deepcopy(no_copy)
    wrong_path["resolved"]["selected_decode_path"] = "tiled"
    mutations.append(wrong_path)
    wrong_digest = copy.deepcopy(no_copy)
    wrong_digest["workload_digests"]["recovered_originals"] = "g" * 16
    mutations.append(wrong_digest)
    wrong_evaluator = copy.deepcopy(no_copy)
    wrong_evaluator["high_evaluator_attribution"]["pruned_output_blocks"] = 1
    wrong_evaluator["high_evaluator_attribution"]["mature_output_blocks"] = 1
    mutations.append(wrong_evaluator)
    wrong_timing = copy.deepcopy(no_copy)
    wrong_timing["metrics"]["decode_execution"]["median_us_per_batch_call"] = 5.0
    mutations.append(wrong_timing)
    bool_control = copy.deepcopy(no_copy)
    bool_control["parameters"]["reuse"] = True
    mutations.append(bool_control)
    bool_loss = copy.deepcopy(no_copy)
    bool_loss["parameters"]["missing_original_indices"] = [True]
    mutations.append(bool_loss)
    bool_tail = copy.deepcopy(no_copy)
    bool_tail["resolved"]["decode_tail_bytes"] = True
    mutations.append(bool_tail)
    for mutation in mutations:
        try:
            validate_pair(mutation, fallback, workspace="materialized", field="gf8",
                          tail_bytes=1, evaluator="mature")
        except ContractError:
            continue
        raise ContractError(
            "adversarial mode/counter/path/digest/evaluator mutation passed")
    print("high-decode copy attribution contract self-test passed: 9 mutations rejected")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("self-test")
    smoke_parser = subparsers.add_parser("smoke")
    smoke_parser.add_argument("--benchmark", required=True, type=Path)
    options = parser.parse_args()
    if options.command == "self-test":
        self_test()
    else:
        smoke(options.benchmark)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ContractError, OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"benchmark_contract.py: {error}", file=__import__("sys").stderr)
        raise SystemExit(1)
