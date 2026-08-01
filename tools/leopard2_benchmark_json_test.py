#!/usr/bin/env python3
"""Regression-check bench_leopard2's default and external-evidence JSON shapes."""

from __future__ import annotations

import copy
import importlib.util
import json
import math
import os
import re
import signal
import stat
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

MAX_JSON_BYTES = 64 * 1024 * 1024
PRODUCTION_DECODE_PATH_RULE_PAIRS = frozenset({
    ("no_op", "no_op"),
    ("direct", "direct"),
    ("generic", "forced_generic"),
    ("materialized", "forced_materialized"),
    ("tiled", "forced_tiled"),
    ("generic", "balanced_generic"),
    ("tiled", "measured_batch_tiled"),
    ("materialized", "measured_materialized"),
    ("tiled", "workspace_tiled"),
    ("materialized", "workspace_materialized"),
    ("materialized", "unsupported_profile"),
})
_PROCESS_RUNNER: Any = None


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def process_runner() -> Any:
    global _PROCESS_RUNNER
    if _PROCESS_RUNNER is None:
        path = (
            Path(__file__).resolve().parents[1] / "experiments" / "leopard2" /
            "main_compare" / "run_abba.py")
        spec = importlib.util.spec_from_file_location(
            "leopard2_benchmark_json_process_runner", path)
        require(spec is not None and spec.loader is not None,
                "cannot load bounded benchmark process runner")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        _PROCESS_RUNNER = module
    return _PROCESS_RUNNER


def run_process(command: list[str], timeout: float = 120.0) -> Any:
    if not sys.platform.startswith("linux"):
        process = subprocess.Popen(
            command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=(os.name == "posix"),
            creationflags=(
                subprocess.CREATE_NEW_PROCESS_GROUP
                if os.name == "nt" else 0))
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired as error:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
            process.communicate()
            raise RuntimeError(
                f"bounded benchmark command exceeded {timeout:.3f} seconds") \
                from error
        require(len(stdout) <= 1024 * 1024 and len(stderr) <= 1024 * 1024,
                "bounded benchmark command exceeded its output limit")
        return subprocess.CompletedProcess(
            command, process.returncode, stdout, stderr)
    try:
        return process_runner().run_process_bounded(
            command, timeout=timeout, max_stdout=1024 * 1024,
            max_stderr=1024 * 1024)
    except Exception as error:
        raise RuntimeError(f"bounded benchmark command failed: {error}") from error


def load_json(path: Path) -> dict[str, Any]:
    def pairs(items: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in items:
            require(key not in result, f"duplicate JSON key {key!r}")
            result[key] = value
        return result

    def constant(value: str) -> object:
        raise RuntimeError(f"non-standard JSON constant {value!r}")

    def floating(value: str) -> float:
        result = float(value)
        require(math.isfinite(result), f"non-finite JSON number {value!r}")
        return result

    path = path.absolute()
    require(path.parent.resolve(strict=True) == path.parent and
            stat.S_ISDIR(os.lstat(path.parent).st_mode),
            "benchmark JSON parent traverses a symbolic link")
    metadata = os.lstat(path)
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
            0 < metadata.st_size <= MAX_JSON_BYTES,
            "benchmark JSON is not a bounded single-link regular file")
    descriptor = os.open(
        path, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))
    try:
        initial = os.fstat(descriptor)
        require((initial.st_mode, initial.st_nlink, initial.st_size,
                 initial.st_mtime_ns, initial.st_ctime_ns, initial.st_dev,
                 initial.st_ino) ==
                (metadata.st_mode, metadata.st_nlink, metadata.st_size,
                 metadata.st_mtime_ns, metadata.st_ctime_ns, metadata.st_dev,
                 metadata.st_ino),
                "benchmark JSON changed before its read")
        payload = bytearray()
        while True:
            block = os.read(descriptor, min(
                1024 * 1024, MAX_JSON_BYTES + 1 - len(payload)))
            if not block:
                break
            payload.extend(block)
            require(len(payload) <= MAX_JSON_BYTES,
                    "benchmark JSON exceeds its byte bound")
        final = os.fstat(descriptor)
        path_final = os.lstat(path)
        identity = lambda value: (
            value.st_mode, value.st_nlink, value.st_size, value.st_mtime_ns,
            value.st_ctime_ns, value.st_dev, value.st_ino)
        require(identity(initial) == identity(final) == identity(path_final) and
                len(payload) == final.st_size,
                "benchmark JSON changed during its read")
    finally:
        os.close(descriptor)
    try:
        value = json.loads(
            bytes(payload).decode("utf-8", errors="strict"),
            object_pairs_hook=pairs, parse_constant=constant,
            parse_float=floating)
    except (UnicodeDecodeError, ValueError, OverflowError,
            RecursionError) as error:
        raise RuntimeError(f"invalid strict benchmark JSON: {error}") from error
    require(isinstance(value, dict), "benchmark JSON is not an object")
    return value


def require_finite_numbers(value: object, label: str) -> None:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return
    if isinstance(value, (int, float)):
        try:
            converted = float(value)
        except (OverflowError, ValueError) as error:
            raise RuntimeError(f"{label} is outside finite range") from error
        require(math.isfinite(converted), f"{label} is not finite")
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            require_finite_numbers(item, f"{label}[{index}]")
        return
    if isinstance(value, dict):
        for key, item in value.items():
            require_finite_numbers(item, f"{label}.{key}")
        return
    raise RuntimeError(f"{label} has an unsupported JSON value")


def finite_metric(
    value: object,
    label: str,
    *,
    positive: bool = False,
) -> float:
    require(type(value) in {int, float}, f"{label} is not numeric")
    result = float(value)
    require(math.isfinite(result) and
            (result > 0.0 if positive else result >= 0.0),
            f"{label} is not "
            f"{'positive' if positive else 'nonnegative and finite'}")
    return result


def median(values: list[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) * 0.5


def approximately_equal(left: float, right: float) -> bool:
    return abs(left - right) <= max(
        0.000003, 0.0000005 * max(abs(left), abs(right)))


def validate_timing_summary(
    summary: object,
    label: str,
    *,
    retain_samples: bool,
    iterations: int,
    execution: bool,
    input_rate_name: str | None = None,
    output_rate_name: str | None = None,
    input_bytes: int = 0,
    output_bytes: int = 0,
) -> None:
    require(isinstance(summary, dict), f"{label} is not an object")
    suffix = "_us_per_batch_call" if execution else "_us"
    fields = {
        f"median{suffix}", f"mad{suffix}",
        f"minimum{suffix}", f"maximum{suffix}",
    }
    sample_name = (
        "samples_us_per_batch_call" if execution else "samples_us")
    if retain_samples:
        fields.add(sample_name)
    if execution:
        require(input_rate_name is not None and output_rate_name is not None,
                f"{label} rate names are absent")
        fields.update({input_rate_name, output_rate_name})
    require(set(summary) == fields, f"{label} fields changed")

    minimum = finite_metric(
        summary[f"minimum{suffix}"], f"{label} minimum",
        positive=execution)
    maximum = finite_metric(
        summary[f"maximum{suffix}"], f"{label} maximum",
        positive=execution)
    central = finite_metric(
        summary[f"median{suffix}"], f"{label} median",
        positive=execution)
    deviation = finite_metric(summary[f"mad{suffix}"], f"{label} MAD")
    require(minimum <= central <= maximum and
            deviation <= maximum - minimum + 0.000003,
            f"{label} summary ordering is inconsistent")

    samples = summary.get(sample_name)
    require((isinstance(samples, list) and
             len(samples) == iterations) == retain_samples,
            f"{label} raw-sample cardinality changed")
    if retain_samples:
        observed = [
            finite_metric(item, f"{label} sample {index}",
                          positive=execution)
            for index, item in enumerate(samples)
        ]
        observed_median = median(observed)
        observed_mad = median([
            abs(item - observed_median) for item in observed])
        require(
            approximately_equal(minimum, min(observed)) and
            approximately_equal(maximum, max(observed)) and
            approximately_equal(central, observed_median) and
            approximately_equal(deviation, observed_mad),
            f"{label} summary is not derived from retained samples")

    if execution:
        for rate_name, byte_count in (
                (input_rate_name, input_bytes),
                (output_rate_name, output_bytes)):
            rate = summary[rate_name]
            if byte_count == 0:
                require(rate is None,
                        f"{label} {rate_name} is non-null for zero bytes")
            else:
                finite_metric(rate, f"{label} {rate_name}", positive=True)


def require_schema_modes_rejected(
    executable: Path,
    arguments: tuple[str, ...],
    expected_diagnostic: str,
) -> None:
    completed = run_process([str(executable), *arguments], timeout=30.0)
    require(completed.returncode == 1,
            f"benchmark accepted incompatible schema modes {arguments!r}")
    stdout = completed.stdout.decode("utf-8", errors="strict")
    stderr = completed.stderr.decode("utf-8", errors="strict")
    require(stdout == "",
            "schema-mode rejection unexpectedly wrote stdout")
    require(stderr ==
            f"leopard2 benchmark: {expected_diagnostic}\n",
            "schema-mode rejection diagnostic changed: "
            f"{stderr!r}")


def run(
    executable: Path,
    external_evidence: bool,
    report_decode_path: bool = False,
    attest_source: bool = False,
    report_direct_executor: bool = False,
    measure_one_shot_decode: bool = False,
    disable_k16r8_b256_terminal: bool = False,
    disable_k9r5_b256_terminal: bool = False,
    disable_k9r6r8_b256_terminal: bool = False,
    *,
    k: int = 3,
    r: int = 2,
    losses: int = 1,
) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="leo2-benchmark-json-") as temporary:
        output = Path(temporary) / "result.json"
        command = [
            str(executable), "--k", str(k), "--r", str(r), "--profile", "high",
            "--field", "auto", "--backend", "auto", "--bytes", "64",
            "--loss", str(losses), "--batch", "1", "--reuse", "1",
            "--iterations", "1", "--warmup", "0", "--threads", "1",
            "--seed", "7",
        ]
        if external_evidence:
            command.extend(("--skip-legacy", "--retain-samples"))
        if report_decode_path:
            command.append("--report-decode-path")
        if attest_source:
            command.append("--attest-source")
        if report_direct_executor:
            command.append("--report-direct-executor")
        if measure_one_shot_decode:
            command.append("--measure-one-shot-decode")
        if disable_k16r8_b256_terminal:
            command.append("--disable-k16r8-b256-terminal")
        if disable_k9r5_b256_terminal:
            command.append("--disable-k9r5-b256-terminal")
        if disable_k9r6r8_b256_terminal:
            command.append("--disable-k9r6r8-b256-terminal")
        command.extend(("--json", str(output)))
        completed = run_process(command)
        require(completed.returncode == 0,
                "benchmark failed: " +
                completed.stdout.decode("utf-8", errors="replace") + "\n" +
                completed.stderr.decode("utf-8", errors="replace"))
        require(not completed.stdout and not completed.stderr,
                "benchmark emitted output while writing JSON")
        return load_json(output)


def validate_common(document: dict[str, Any], retain_samples: bool) -> None:
    require(isinstance(document, dict) and type(retain_samples) is bool,
            "benchmark document/retention mode is invalid")
    require_finite_numbers(document, "benchmark")
    require(type(document.get("schema")) is str and
            document["schema"] in {
                "leopard2-benchmark-v1", "leopard2-benchmark-v2",
                "leopard2-benchmark-v3", "leopard2-benchmark-v5",
                "leopard2-benchmark-v6", "leopard2-benchmark-v7",
                "leopard2-benchmark-v8", "leopard2-benchmark-v9",
                "leopard2-benchmark-v10",
            }, "benchmark schema is unsupported")
    expected_top = {
        "schema", "build", "parameters", "resolved", "correctness",
        "memory", "metrics", "legacy"}
    if document["schema"] in {
        "leopard2-benchmark-v2", "leopard2-benchmark-v3",
        "leopard2-benchmark-v5", "leopard2-benchmark-v6",
        "leopard2-benchmark-v7", "leopard2-benchmark-v8",
        "leopard2-benchmark-v9", "leopard2-benchmark-v10",
    }:
        expected_top.add("workload_digests")
    require(set(document) == expected_top, "top-level JSON keys changed")
    for section in (
            "build", "parameters", "resolved", "correctness", "memory",
            "metrics", "legacy"):
        require(isinstance(document[section], dict),
                f"benchmark {section} is not an object")
    expected_build = {
        "compiler", "compiler_version", "cplusplus",
        "k16r8_b256_terminal_diagnostic_disabled",
        "k9r5_b256_terminal_diagnostic_disabled",
    }
    if document["schema"] == "leopard2-benchmark-v10":
        expected_build.add("k9r6r8_b256_terminal_diagnostic_disabled")
    if document["schema"] in {
            "leopard2-benchmark-v5", "leopard2-benchmark-v7"}:
        expected_build.update({
            "source_commit", "source_tree", "source_tracked_dirty"})
    if document["schema"] in {
            "leopard2-benchmark-v6", "leopard2-benchmark-v7"}:
        expected_build.add("equal_rounded_multi_loss_enabled")
    if document["schema"] in {
            "leopard2-benchmark-v8", "leopard2-benchmark-v9"}:
        expected_build.add("one_shot_equal_rounded_direct_enabled")
    if document["schema"] == "leopard2-benchmark-v9":
        expected_build.update({
            "one_shot_equal_rounded_direct_enabled",
            "cauchy_log_reuse_enabled",
        })
    require(set(document["build"]) == expected_build, "build keys changed")
    require(type(document["build"]["compiler"]) is str and
            type(document["build"]["compiler_version"]) is str and
            type(document["build"]["cplusplus"]) is int and
            type(document["build"][
                "k16r8_b256_terminal_diagnostic_disabled"]) is bool and
            type(document["build"][
                "k9r5_b256_terminal_diagnostic_disabled"]) is bool and
            document["build"]["cplusplus"] > 0,
            "benchmark build value types changed")
    if document["schema"] == "leopard2-benchmark-v10":
        require(type(document["build"][
                    "k9r6r8_b256_terminal_diagnostic_disabled"]) is bool,
                "K9/R=6..8 terminal selector is not Boolean")
    if document["schema"] in {
            "leopard2-benchmark-v6", "leopard2-benchmark-v7"}:
        require(type(document["build"][
                    "equal_rounded_multi_loss_enabled"]) is bool,
                "equal-rounded build selector is not Boolean")
    if document["schema"] in {
            "leopard2-benchmark-v8", "leopard2-benchmark-v9"}:
        require(type(document["build"][
                    "one_shot_equal_rounded_direct_enabled"]) is bool,
                "one-shot build selector is not Boolean")
    if document["schema"] == "leopard2-benchmark-v9":
        require(type(document["build"][
                    "cauchy_log_reuse_enabled"]) is bool,
                "Cauchy-log-reuse build selector is not Boolean")
    if document["schema"] in {
            "leopard2-benchmark-v5", "leopard2-benchmark-v7"}:
        for name in ("source_commit", "source_tree"):
            value = document["build"][name]
            require(value == "unknown" or
                    (isinstance(value, str) and len(value) == 40 and
                     all(character in "0123456789abcdef" for character in value)),
                    f"benchmark {name} is not a lowercase Git identity")
        require(isinstance(document["build"]["source_tracked_dirty"], bool),
                "benchmark dirty state is not Boolean")
    expected_resolved = {
        "profile", "field", "backend", "thread_count", "parent_count",
        "padded_side"}
    if document["schema"] in {
        "leopard2-benchmark-v3", "leopard2-benchmark-v6",
        "leopard2-benchmark-v7",
    }:
        expected_resolved.update({
            "selected_decode_path", "selected_decode_rule",
            "decode_required_work_slots", "decode_aligned_prefix_bytes",
            "decode_tail_bytes", "decode_rounded_bytes",
            "decode_multi_item_batch",
        })
    if document["schema"] in {
            "leopard2-benchmark-v6", "leopard2-benchmark-v7"}:
        expected_resolved.add("selected_direct_executor")
    require(set(document["resolved"]) == expected_resolved,
            "resolved keys changed")
    parameters = document["parameters"]
    for name in ("K", "R", "shard_bytes", "batch", "reuse", "iterations",
                 "thread_count"):
        require(type(parameters.get(name)) is int and parameters[name] > 0,
                f"benchmark parameter {name} is invalid")
    require(type(parameters.get("loss_count")) is int and
            0 <= parameters["loss_count"] <=
                min(parameters["K"], parameters["R"]) and
            isinstance(parameters.get("missing_original_indices"), list) and
            len(parameters["missing_original_indices"]) ==
                parameters["loss_count"] and
            all(type(index) is int and 0 <= index < parameters["K"]
                for index in parameters["missing_original_indices"]) and
            len(set(parameters["missing_original_indices"])) ==
                parameters["loss_count"],
            "benchmark loss parameters are invalid")
    require(type(document["resolved"]["thread_count"]) is int and
            document["resolved"]["thread_count"] > 0 and
            type(document["resolved"]["parent_count"]) is int and
            document["resolved"]["parent_count"] > 0 and
            type(document["resolved"]["padded_side"]) is int and
            document["resolved"]["padded_side"] > 0,
            "resolved numeric codec identity is invalid")
    if document["schema"] in {
        "leopard2-benchmark-v3", "leopard2-benchmark-v6",
        "leopard2-benchmark-v7",
    }:
        pair = (
            document["resolved"]["selected_decode_path"],
            document["resolved"]["selected_decode_rule"])
        require(all(type(item) is str for item in pair) and
                pair in PRODUCTION_DECODE_PATH_RULE_PAIRS,
                "resolved decode path/rule pair is not a production pair")
    if document["schema"] in {
            "leopard2-benchmark-v6", "leopard2-benchmark-v7"}:
        executor = document["resolved"]["selected_direct_executor"]
        require(type(executor) is str and
                executor in {"none", "output_major", "source_major"} and
                (executor == "none" or pair == ("direct", "direct")),
                "resolved direct executor is inconsistent with decode path")
    require(set(document["correctness"]) == {
        "leopard2_round_trip", "legacy_comparison"}, "correctness keys changed")
    require(document["correctness"]["leopard2_round_trip"] is True and
            (document["correctness"]["legacy_comparison"] is None or
             document["correctness"]["legacy_comparison"] == "matched"),
            "correctness results are invalid")
    expected_memory = {
        "scratch_alignment", "encode_scratch_bytes_per_stripe",
        "decode_scratch_bytes_per_stripe", "encode_scratch_bytes_batch",
        "decode_scratch_bytes_batch"}
    if document["schema"] in {
            "leopard2-benchmark-v8", "leopard2-benchmark-v9"}:
        expected_memory.update({
            "one_shot_decode_scratch_bytes_per_stripe",
            "one_shot_decode_scratch_bytes_batch",
        })
    require(set(document["memory"]) == expected_memory,
            "memory keys changed")
    require(all(type(value) is int and value >= 0
                for value in document["memory"].values()) and
            document["memory"]["scratch_alignment"] > 0,
            "memory values are invalid")
    expected_metrics = {
        "codec_setup", "encode_execution", "decode_plan_setup",
        "decode_execution", "decode_amortized_at_reuse", "rate_semantics"}
    if document["schema"] in {
            "leopard2-benchmark-v8", "leopard2-benchmark-v9"}:
        expected_metrics.add("one_shot_decode_including_setup")
    require(set(document["metrics"]) == expected_metrics,
        "metrics keys changed")
    require(set(document["legacy"]) == {
        "available", "unavailable_reason", "codec_setup",
        "decode_timing_includes_setup", "encode_execution",
        "decode_including_setup"}, "legacy keys changed")
    iterations = parameters["iterations"]
    batch = parameters["batch"]
    k = parameters["K"]
    r = parameters["R"]
    shard_bytes = parameters["shard_bytes"]
    losses = parameters["loss_count"]
    encode_input_bytes = k * shard_bytes * batch
    encode_output_bytes = r * shard_bytes * batch
    decode_input_bytes = (k - losses + r) * shard_bytes * batch
    decode_output_bytes = losses * shard_bytes * batch

    validate_timing_summary(
        document["metrics"]["codec_setup"], "codec_setup",
        retain_samples=retain_samples, iterations=iterations, execution=False)
    validate_timing_summary(
        document["metrics"]["decode_plan_setup"], "decode_plan_setup",
        retain_samples=retain_samples, iterations=iterations, execution=False)
    validate_timing_summary(
        document["metrics"]["encode_execution"], "encode_execution",
        retain_samples=retain_samples, iterations=iterations, execution=True,
        input_rate_name="input_GB_per_s",
        output_rate_name="parity_output_GB_per_s",
        input_bytes=encode_input_bytes, output_bytes=encode_output_bytes)
    validate_timing_summary(
        document["metrics"]["decode_execution"], "decode_execution",
        retain_samples=retain_samples, iterations=iterations, execution=True,
        input_rate_name="offered_received_GB_per_s",
        output_rate_name="repaired_output_GB_per_s",
        input_bytes=decode_input_bytes, output_bytes=decode_output_bytes)
    if document["schema"] in {
            "leopard2-benchmark-v8", "leopard2-benchmark-v9"}:
        require(parameters.get("measure_one_shot_decode") is True,
                "one-shot decode benchmark opt-in was not recorded")
        validate_timing_summary(
            document["metrics"]["one_shot_decode_including_setup"],
            "one_shot_decode_including_setup",
            retain_samples=retain_samples, iterations=iterations,
            execution=True,
            input_rate_name="offered_received_GB_per_s",
            output_rate_name="repaired_output_GB_per_s",
            input_bytes=decode_input_bytes,
            output_bytes=decode_output_bytes)

    amortized = document["metrics"]["decode_amortized_at_reuse"]
    require(isinstance(amortized, dict) and set(amortized) == {
        "reuse_count", "derived_median_us_per_batch_call",
        "offered_received_GB_per_s", "repaired_output_GB_per_s",
    } and type(amortized["reuse_count"]) is int and
            amortized["reuse_count"] == parameters["reuse"],
            "decode amortization structure changed")
    derived = finite_metric(
        amortized["derived_median_us_per_batch_call"],
        "decode amortized median", positive=True)
    expected_derived = (
        float(document["metrics"]["decode_execution"][
            "median_us_per_batch_call"]) +
        float(document["metrics"]["decode_plan_setup"]["median_us"]) /
            parameters["reuse"])
    require(approximately_equal(derived, expected_derived),
            "decode amortization is inconsistent with setup/execution medians")
    for name, byte_count in (
            ("offered_received_GB_per_s", decode_input_bytes),
            ("repaired_output_GB_per_s", decode_output_bytes)):
        if byte_count == 0:
            require(amortized[name] is None,
                    f"decode amortized {name} is non-null for zero bytes")
        else:
            finite_metric(
                amortized[name], f"decode amortized {name}", positive=True)
    require(document["metrics"]["rate_semantics"] ==
            "offered_received counts all non-null shard pointers supplied; "
            "a plan may read a deterministic subset",
            "benchmark rate semantics changed")

    legacy = document["legacy"]
    require(type(legacy["available"]) is bool and
            legacy["codec_setup"] is None and
            legacy["decode_timing_includes_setup"] is True,
            "legacy timing contract changed")
    if legacy["available"]:
        require(legacy["unavailable_reason"] is None and
                document["correctness"]["legacy_comparison"] == "matched",
                "available legacy comparison is inconsistent")
        validate_timing_summary(
            legacy["encode_execution"], "legacy encode_execution",
            retain_samples=retain_samples, iterations=iterations,
            execution=True, input_rate_name="input_GB_per_s",
            output_rate_name="parity_output_GB_per_s",
            input_bytes=encode_input_bytes, output_bytes=encode_output_bytes)
        validate_timing_summary(
            legacy["decode_including_setup"], "legacy decode_including_setup",
            retain_samples=retain_samples, iterations=iterations,
            execution=True, input_rate_name="offered_received_GB_per_s",
            output_rate_name="repaired_output_GB_per_s",
            input_bytes=decode_input_bytes, output_bytes=decode_output_bytes)
    else:
        require(type(legacy["unavailable_reason"]) is str and
                bool(legacy["unavailable_reason"]) and
                legacy["encode_execution"] is None and
                legacy["decode_including_setup"] is None and
                document["correctness"]["legacy_comparison"] is None,
                "unavailable legacy comparison is inconsistent")


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


def validate_direct_executor_report(
    document: dict[str, Any],
    expected_executor: str,
    *,
    expect_direct_path: bool,
    expected_schema: str = "leopard2-benchmark-v6",
) -> None:
    require(document.get("schema") == expected_schema,
            "direct-executor benchmark schema changed")
    validate_common(document, True)
    validate_workload_digests(document)
    parameters = document["parameters"]
    require(parameters.get("report_decode_path") is True and
            parameters.get("report_direct_executor") is True,
            "direct-executor opt-in was not completely recorded")
    resolved = document["resolved"]
    executor = resolved.get("selected_direct_executor")
    require(executor in {"none", "output_major", "source_major"},
            "direct-executor metadata is malformed")
    require(executor == expected_executor,
            "direct-executor metadata differs from the selected production path")
    if expect_direct_path:
        require(resolved.get("selected_decode_path") == "direct" and
                resolved.get("selected_decode_rule") == "direct",
                "direct-executor report did not observe direct repair")
    else:
        require(resolved.get("selected_decode_path") in {
                    "generic", "materialized", "tiled"} and
                resolved.get("selected_decode_rule") not in {
                    "no_op", "direct", "unsupported_profile"} and
                executor == "none",
                "direct-executor report did not observe transform repair")


def require_direct_executor_rejected(
    document: dict[str, Any],
    expected_executor: str,
    *,
    expect_direct_path: bool,
    label: str,
) -> None:
    try:
        validate_direct_executor_report(
            document, expected_executor,
            expect_direct_path=expect_direct_path)
    except RuntimeError:
        return
    raise RuntimeError(f"direct-executor validator accepted {label}")


def require_common_rejected(document: dict[str, Any], label: str) -> None:
    try:
        validate_common(document, True)
    except (KeyError, RuntimeError, TypeError, ValueError):
        return
    raise RuntimeError(f"benchmark validator accepted {label}")


def validate_prevalidated_batch_variant(
    document: dict[str, Any],
    retain_samples: bool,
) -> None:
    build_fields = {
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
    }
    metric_fields = {"encode_binding_setup", "decode_binding_setup"}
    require(isinstance(document.get("build"), dict) and
            build_fields <= set(document["build"]),
            "prevalidated benchmark build markers are incomplete")
    require(isinstance(document.get("metrics"), dict) and
            metric_fields <= set(document["metrics"]),
            "prevalidated benchmark setup metrics are incomplete")

    normalized = copy.deepcopy(document)
    for name in build_fields:
        del normalized["build"][name]
    for name in metric_fields:
        del normalized["metrics"][name]
    validate_common(normalized, retain_samples)

    build = document["build"]
    require(build["prevalidated_batch_experiment"] is True and
            all(type(build[name]) is bool for name in build_fields),
            "prevalidated benchmark build markers are not Boolean")
    require(build["high_t4_batch_selected"] is False and
            build["high_t8_one_block_selected"] is False and
            build["high_t8_two_block_selected"] is False,
            "K=1/R=1 unexpectedly selected a T=4/T=8 encode binding")

    iterations = document["parameters"]["iterations"]
    for name in sorted(metric_fields):
        validate_timing_summary(
            document["metrics"][name], name,
            retain_samples=retain_samples, iterations=iterations,
            execution=False)
        finite_metric(
            document["metrics"][name]["median_us"],
            f"{name} median", positive=True)


def require_prevalidated_rejected(
    document: dict[str, Any], label: str,
) -> None:
    try:
        validate_prevalidated_batch_variant(document, True)
    except (KeyError, RuntimeError, TypeError, ValueError):
        return
    raise RuntimeError(
        f"prevalidated benchmark validator accepted {label}")


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
    if len(sys.argv) not in (2, 3):
        raise RuntimeError(
            "usage: leopard2_benchmark_json_test.py BENCH_LEOPARD2 "
            "[BENCH_LEOPARD2_ALLK]")
    executable = Path(sys.argv[1]).resolve()
    default = run(executable, False)
    require(default["schema"] == "leopard2-benchmark-v1",
            "default benchmark schema changed")
    validate_common(default, False)
    require(default["build"][
                "k16r8_b256_terminal_diagnostic_disabled"] is False,
            "default benchmark disabled the K16/R8 terminal")
    require(default["build"][
                "k9r5_b256_terminal_diagnostic_disabled"] is False,
            "default benchmark disabled the K9/R5 terminal")
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
    require(external["build"][
                "k16r8_b256_terminal_diagnostic_disabled"] is False,
            "external-evidence benchmark disabled the K16/R8 terminal")
    require(external["build"][
                "k9r5_b256_terminal_diagnostic_disabled"] is False,
            "external-evidence benchmark disabled the K9/R5 terminal")
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

    disabled_k16_terminal = run(
        executable, True, disable_k16r8_b256_terminal=True)
    validate_common(disabled_k16_terminal, True)
    validate_workload_digests(disabled_k16_terminal)
    require(disabled_k16_terminal["build"][
                "k16r8_b256_terminal_diagnostic_disabled"] is True,
            "K16/R8 attribution option was not recorded")
    require(disabled_k16_terminal["workload_digests"] ==
                external["workload_digests"],
            "inert K16/R8 attribution option changed the workload")

    disabled_k9_terminal = run(
        executable, True, disable_k9r5_b256_terminal=True)
    validate_common(disabled_k9_terminal, True)
    validate_workload_digests(disabled_k9_terminal)
    require(disabled_k9_terminal["build"][
                "k9r5_b256_terminal_diagnostic_disabled"] is True,
            "K9/R5 attribution option was not recorded")
    require(disabled_k9_terminal["workload_digests"] ==
                external["workload_digests"],
            "inert K9/R5 attribution option changed the workload")

    disabled_k9r6r8_terminal = run(
        executable, True, disable_k9r6r8_b256_terminal=True)
    require(disabled_k9r6r8_terminal["schema"] ==
                "leopard2-benchmark-v10",
            "K9/R=6..8 attribution did not version its build marker")
    validate_common(disabled_k9r6r8_terminal, True)
    validate_workload_digests(disabled_k9r6r8_terminal)
    require(disabled_k9r6r8_terminal["build"][
                "k9r6r8_b256_terminal_diagnostic_disabled"] is True,
            "K9/R=6..8 attribution option was not recorded")
    require(disabled_k9r6r8_terminal["workload_digests"] ==
                external["workload_digests"],
            "inert K9/R=6..8 attribution option changed the workload")

    prevalidated_executable = executable.with_name(
        "bench_leopard2_prevalidated_batch" + executable.suffix)
    require(prevalidated_executable.is_file(),
            "prevalidated benchmark target is absent")
    prevalidated = run(
        prevalidated_executable, True, k=1, r=1, losses=1)
    require(prevalidated["schema"] == "leopard2-benchmark-v2",
            "prevalidated benchmark external schema changed")
    validate_prevalidated_batch_variant(prevalidated, True)
    validate_workload_digests(prevalidated)
    require(prevalidated["resolved"]["profile"] == "legacy_high_v1" and
            prevalidated["resolved"]["field"] == "gf8" and
            prevalidated["correctness"]["leopard2_round_trip"] is True,
            "prevalidated K=1/R=1 binding round trip changed")

    missing_binding_metric = copy.deepcopy(prevalidated)
    del missing_binding_metric["metrics"]["decode_binding_setup"]
    require_prevalidated_rejected(
        missing_binding_metric, "missing decode binding setup")
    malformed_binding_marker = copy.deepcopy(prevalidated)
    malformed_binding_marker["build"][
        "prevalidated_batch_experiment"] = 1
    require_prevalidated_rejected(
        malformed_binding_marker, "numeric experiment marker")

    one_shot = run(
        executable, True, measure_one_shot_decode=True,
        k=17, r=17, losses=8)
    require(one_shot["schema"] == "leopard2-benchmark-v9",
            "one-shot decode benchmark schema changed")
    validate_common(one_shot, True)
    validate_workload_digests(one_shot)
    require(set(one_shot["parameters"]) ==
            (expected_external_parameters | {"measure_one_shot_decode"}),
            "one-shot decode parameter structure changed")
    require(one_shot["parameters"]["measure_one_shot_decode"] is True,
            "one-shot decode opt-in was not recorded")
    require(one_shot["memory"][
                "one_shot_decode_scratch_bytes_per_stripe"] > 0 and
            one_shot["memory"]["one_shot_decode_scratch_bytes_batch"] ==
                one_shot["memory"][
                    "one_shot_decode_scratch_bytes_per_stripe"],
            "one-shot decode scratch accounting changed")

    path_diagnostic = (
        "--attest-source and --report-decode-path use distinct JSON schemas")
    for arguments in (
            ("--attest-source", "--report-decode-path"),
            ("--report-decode-path", "--attest-source")):
        require_schema_modes_rejected(
            executable, arguments, path_diagnostic)
    one_shot_diagnostic = (
        "--measure-one-shot-decode currently uses a standalone schema and "
        "cannot be combined with path/source attestation")
    for incompatible in (
            "--attest-source", "--report-decode-path",
            "--report-direct-executor"):
        for arguments in (
                ("--measure-one-shot-decode", incompatible),
                (incompatible, "--measure-one-shot-decode")):
            require_schema_modes_rejected(
                executable, arguments, one_shot_diagnostic)

    standard_attested = run(executable, True, attest_source=True)
    require(standard_attested["schema"] == "leopard2-benchmark-v5",
            "standard benchmark source-attested schema changed")
    validate_common(standard_attested, True)
    validate_workload_digests(standard_attested)
    require(set(standard_attested["parameters"]) ==
            (expected_external_parameters | {"attest_source"}),
            "standard source-attested parameter structure changed")
    require(standard_attested["parameters"]["attest_source"] is True,
            "standard benchmark did not record source-attestation opt-in")
    expected_header = os.environ.get(
        "LEO2_EXPECTED_SOURCE_ATTESTATION_HEADER")
    expected_commit = os.environ.get("LEO2_EXPECTED_SOURCE_COMMIT")
    expected_tree = os.environ.get("LEO2_EXPECTED_SOURCE_TREE")
    expected_dirty = os.environ.get("LEO2_EXPECTED_SOURCE_DIRTY")
    if expected_header is not None:
        header = Path(expected_header).read_text(encoding="utf-8")

        def macro(name: str) -> str:
            match = re.search(
                rf"^#define {re.escape(name)} (.+)$", header, re.MULTILINE)
            require(match is not None,
                    f"generated source attestation omits {name}")
            return match.group(1)

        expected_commit = macro(
            "LEO2_BENCHMARK_SOURCE_COMMIT").strip('"')
        expected_tree = macro("LEO2_BENCHMARK_SOURCE_TREE").strip('"')
        expected_dirty = macro(
            "LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY")
    if expected_commit is not None:
        require(standard_attested["build"]["source_commit"] == expected_commit,
                "standard benchmark source commit differs from CMake identity")
        require(standard_attested["build"]["source_tree"] == expected_tree,
                "standard benchmark source tree differs from CMake identity")
        require(standard_attested["build"]["source_tracked_dirty"] is
                (expected_dirty != "0"),
                "standard benchmark dirty flag differs from CMake identity")

    attested_direct = run(
        executable, True, attest_source=True,
        report_direct_executor=True, k=3, r=2, losses=2)
    validate_direct_executor_report(
        attested_direct, "output_major", expect_direct_path=True,
        expected_schema="leopard2-benchmark-v7")
    require(set(attested_direct["parameters"]) ==
            (expected_external_parameters | {
                "attest_source", "report_decode_path",
                "report_direct_executor"}),
            "source-attested direct-executor parameters changed")
    require(all(attested_direct["build"][name] ==
                standard_attested["build"][name]
                for name in (
                    "compiler", "compiler_version", "cplusplus",
                    "source_commit", "source_tree",
                    "source_tracked_dirty")),
            "source-attested direct executor changed source identity")

    if len(sys.argv) == 3:
        attested_executable = Path(sys.argv[2]).resolve()
        attested_default = run(attested_executable, False)
        require(attested_default["schema"] == "leopard2-benchmark-v1",
                "attested target changed its default schema without opt-in")
        validate_common(attested_default, False)
        require(attested_default["build"] == default["build"],
                "attested target exposed source identity without opt-in")
        attested = run(attested_executable, True, attest_source=True)
        require(attested["schema"] == "leopard2-benchmark-v5",
                "all-K source-attested benchmark schema changed")
        validate_common(attested, True)
        validate_workload_digests(attested)
        require(attested["build"] == standard_attested["build"],
                "standard and all-K source identities differ")

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

    direct_executor = run(
        executable, True, report_direct_executor=True,
        k=3, r=2, losses=2)
    validate_direct_executor_report(
        direct_executor, "output_major", expect_direct_path=True)
    require(set(direct_executor["parameters"]) ==
            (expected_external_parameters |
             {"report_decode_path", "report_direct_executor"}),
            "direct-executor parameter structure changed")

    direct_xor_executor = run(
        executable, True, report_direct_executor=True,
        k=1, r=1, losses=1)
    validate_direct_executor_report(
        direct_xor_executor, "none", expect_direct_path=True)
    require(set(direct_xor_executor["parameters"]) ==
            set(direct_executor["parameters"]),
            "direct-XOR executor parameter structure changed")

    transform_executor = run(
        executable, True, report_direct_executor=True,
        # K=5 is intentionally eligible in small-direct experiment modes.
        # Keep this production smoke outside every current direct-repair
        # family so modes 0, 1, and 2 all exercise a transform executor.
        k=17, r=5, losses=5)
    validate_direct_executor_report(
        transform_executor, "none", expect_direct_path=False)
    require(set(transform_executor["parameters"]) ==
            set(direct_executor["parameters"]),
            "transform direct-executor parameter structure changed")

    missing_executor = json.loads(json.dumps(direct_executor))
    del missing_executor["resolved"]["selected_direct_executor"]
    require_direct_executor_rejected(
        missing_executor, "output_major", expect_direct_path=True,
        label="missing selected_direct_executor")

    malformed_executor = json.loads(json.dumps(direct_executor))
    malformed_executor["resolved"]["selected_direct_executor"] = "output-major"
    require_direct_executor_rejected(
        malformed_executor, "output_major", expect_direct_path=True,
        label="malformed selected_direct_executor")

    inconsistent_executor = json.loads(json.dumps(transform_executor))
    inconsistent_executor["resolved"]["selected_direct_executor"] = \
        "output_major"
    require_direct_executor_rejected(
        inconsistent_executor, "none", expect_direct_path=False,
        label="transform path with a direct executor")

    unknown_pair = json.loads(json.dumps(path_report))
    unknown_pair["resolved"]["selected_decode_rule"] = "totally_bogus"
    require_common_rejected(unknown_pair, "unknown decode path/rule pair")

    cross_pair = json.loads(json.dumps(path_report))
    cross_pair["resolved"]["selected_decode_path"] = "generic"
    cross_pair["resolved"]["selected_decode_rule"] = "workspace_tiled"
    require_common_rejected(cross_pair, "cross-paired decode path/rule")

    nonfinite_metric = json.loads(json.dumps(external))
    nonfinite_metric["metrics"]["decode_execution"][
        "median_us_per_batch_call"] = float("nan")
    require_common_rejected(nonfinite_metric, "NaN decode timing")

    boolean_metric = json.loads(json.dumps(external))
    boolean_metric["metrics"]["decode_execution"][
        "minimum_us_per_batch_call"] = True
    require_common_rejected(boolean_metric, "Boolean decode timing")

    inconsistent_summary = json.loads(json.dumps(external))
    inconsistent_summary["metrics"]["decode_execution"][
        "median_us_per_batch_call"] += 1.0
    require_common_rejected(
        inconsistent_summary, "timing summary inconsistent with samples")

    numeric_digest = json.loads(json.dumps(external))
    numeric_digest["workload_digests"]["original_data"] = int("1" * 16)
    try:
        validate_workload_digests(numeric_digest)
    except RuntimeError:
        pass
    else:
        raise RuntimeError("benchmark validator accepted a numeric FNV digest")

    with tempfile.TemporaryDirectory(
            prefix="leo2-benchmark-json-adversarial-") as temporary:
        root = Path(temporary)
        ordinary = root / "ordinary.json"
        ordinary.write_text('{"value":1}\n', encoding="utf-8")
        require(load_json(ordinary) == {"value": 1},
                "strict benchmark loader rejected an ordinary document")

        symbolic = root / "symbolic.json"
        symbolic.symlink_to(ordinary)
        try:
            load_json(symbolic)
        except (OSError, RuntimeError):
            pass
        else:
            raise RuntimeError("strict benchmark loader followed a symlink")

        hard_link = root / "hard-link.json"
        os.link(ordinary, hard_link)
        try:
            load_json(hard_link)
        except (OSError, RuntimeError):
            pass
        else:
            raise RuntimeError("strict benchmark loader accepted a hard link")

        nonstandard = root / "nonstandard.json"
        nonstandard.write_text('{"value":NaN}\n', encoding="utf-8")
        try:
            load_json(nonstandard)
        except (OSError, RuntimeError):
            pass
        else:
            raise RuntimeError(
                "strict benchmark loader accepted non-standard NaN")

    if sys.platform.startswith("linux"):
        descendant_code = (
            "import subprocess,sys,time;"
            "subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)'],"
            "start_new_session=True);time.sleep(60)")
        try:
            run_process(
                [sys.executable, "-c", descendant_code], timeout=0.2)
        except RuntimeError as error:
            require("exceeded" in str(error),
                    f"bounded timeout failed for an unexpected reason: {error}")
        else:
            raise RuntimeError(
                "bounded benchmark runner accepted an over-time descendant tree")

    print("leopard2 benchmark JSON regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
