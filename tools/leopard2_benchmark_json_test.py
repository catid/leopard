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
    ("tiled", "translated_low"),
    ("materialized", "translated_low"),
    ("materialized", "unsupported_profile"),
})
R1_REDUCTION_PATHS = frozenset({
    "not_applicable", "k1_copy", "k2_terminal", "pairwise", "dense",
    "coarse", "fused_final", "group4", "fixed_avx2",
})
R1_TIMED_ENCODE_API = (
    "leo2_encode_batch:item_count=1:no_preflight_scratch")
R1_TIMED_ONE_SHOT_ENCODE_API = "leo2_encode"
R1_TIMED_REUSED_DECODE_API = (
    "leo2_decode_plan_execute_batch:item_count=1:"
    "no_preflight_scratch:one_loss_direct_xor")
R1_TIMED_ONE_SHOT_DECODE_API = "leo2_decode:one_loss"
R1_PREVALIDATED_TIMED_ENCODE_API = "leo2_encode_batch_binding_execute"
R1_PREVALIDATED_TIMED_REUSED_DECODE_API = (
    "leo2_decode_batch_binding_execute")
R1_ENCODE_BINDING_SETUP_API = "leo2_encode_batch_binding_create"
R1_DECODE_BINDING_SETUP_API = "leo2_decode_batch_binding_create"
R1_FIXED_AVX2_SELECTOR_CONTRACT = (
    "LEGACY_HIGH_V1,GF8,AVX2,R=1,K=3..255,B=64|256,"
    "native_layout,auto_encode_decode")
K65R65_B64_SELECTOR_CONTRACT = (
    "LEGACY_HIGH_V1,GF8,AVX2,K=65,R=65,T=128,B=64,native_layout,"
    "auto_encode,one_shot_and_one_item_batch")
K65R65_B64_TIMED_ORDINARY_ENCODE_API = (
    "leo2_encode_batch:item_count=1:no_preflight_scratch")
K65R65_B64_TIMED_ONE_SHOT_ENCODE_API = "leo2_encode"
K65R65_B64_AVX512_GFNI_SELECTOR_CONTRACT = (
    "LEGACY_HIGH_V1,GF8,AUTO,K=65,R=65,T=128,B=64,native_layout,"
    "packed_terminal,runtime_AVX512F_BW_VL_GFNI,startup_KAT,"
    "calibrated_AMD_1A_08,one_shot_and_one_item_batch")
ONE_SHOT_SCHEMAS = frozenset({
    "leopard2-benchmark-v8", "leopard2-benchmark-v9",
    "leopard2-benchmark-v12", "leopard2-benchmark-v15",
    "leopard2-benchmark-v16", "leopard2-benchmark-v23",
})
R1_SCHEMAS = frozenset({
    "leopard2-benchmark-v12", "leopard2-benchmark-v23",
})
TRANSIENT_PLAN_SCHEMAS = frozenset({
    "leopard2-benchmark-v15", "leopard2-benchmark-v16",
})
TRANSIENT_SETUP_POLICIES = {
    0: "reusable_metadata_superset",
    1: "exact_byte_selected_metadata",
    2: "exact_byte_tiny_high_regular_metadata",
    3: "exact_byte_tiny_all_regular_metadata",
}
TRANSIENT_EXECUTION_EXCLUDES = [
    "decode_plan_create_destroy",
    "public_one_shot_no_loss_and_direct_dispatch",
]
TRANSIENT_EXECUTION_GUARDS = [
    "transform_route_and_exact_shard_bytes",
    "presence_metadata_overlap",
]
TRANSIENT_SETUP_EXCLUDES = [
    "r1_terminal_probe",
    "equal_rounded_direct_repair_probe",
]
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
    disable_k9r5_b1024_terminal: bool = False,
    disable_k9r6r8_b256_terminal: bool = False,
    r1_small_reduction_mode: int | None = None,
    r1_fixed_avx2_mode: int | None = None,
    prevalidated_binding: bool = False,
    one_shot_plan_setup_mode: int | None = None,
    gf8_avx2_walsh_locator_mode: int | None = None,
    k65r65_b64_packed_terminal_mode: int | None = None,
    k65r65_b64_avx512_gfni_mode: int | None = None,
    measure_one_shot_encode: bool = False,
    *,
    k: int = 3,
    r: int = 2,
    losses: int = 1,
    shard_bytes: int = 64,
) -> dict[str, Any]:
    require(r1_small_reduction_mode is None or
            (type(r1_small_reduction_mode) is int and
             r1_small_reduction_mode in {0, 1}),
            "R=1 small-reduction mode must be absent, zero, or one")
    require(r1_fixed_avx2_mode is None or
            (type(r1_fixed_avx2_mode) is int and
             r1_fixed_avx2_mode in {0, 1}),
            "fixed AVX2 R=1 mode must be absent, zero, or one")
    require(r1_small_reduction_mode is None or
            r1_fixed_avx2_mode is None,
            "R=1 diagnostic modes must be mutually exclusive")
    require(type(prevalidated_binding) is bool and
            (not prevalidated_binding or r1_fixed_avx2_mode is not None),
            "prevalidated binding requires a fixed AVX2 R=1 mode")
    require(one_shot_plan_setup_mode is None or
            (type(one_shot_plan_setup_mode) is int and
             one_shot_plan_setup_mode in {0, 1, 2, 3}),
            "one-shot setup mode must be absent or an integer from zero to three")
    require(gf8_avx2_walsh_locator_mode is None or
            (type(gf8_avx2_walsh_locator_mode) is int and
             gf8_avx2_walsh_locator_mode in {0, 1}),
            "GF8 AVX2 Walsh-locator mode must be absent, zero, or one")
    require(k65r65_b64_packed_terminal_mode is None or
            (type(k65r65_b64_packed_terminal_mode) is int and
             k65r65_b64_packed_terminal_mode in {0, 1}),
            "K65/R65/B64 packed-terminal mode must be absent, zero, or one")
    require(k65r65_b64_avx512_gfni_mode is None or
            (type(k65r65_b64_avx512_gfni_mode) is int and
             k65r65_b64_avx512_gfni_mode in {0, 1}),
            "K65/R65/B64 AVX-512/GFNI mode must be absent, zero, or one")
    require(k65r65_b64_packed_terminal_mode is None or
            k65r65_b64_avx512_gfni_mode is None,
            "K65/R65/B64 diagnostic modes must be mutually exclusive")
    require(type(measure_one_shot_encode) is bool,
            "one-shot encode measurement mode must be Boolean")
    require(type(shard_bytes) is int and shard_bytes > 0,
            "benchmark shard byte count must be a positive integer")
    with tempfile.TemporaryDirectory(prefix="leo2-benchmark-json-") as temporary:
        output = Path(temporary) / "result.json"
        diagnostic_auto_gfni = k65r65_b64_avx512_gfni_mode is not None
        diagnostic_avx2 = (
            r1_small_reduction_mode is not None or
            r1_fixed_avx2_mode is not None or
            one_shot_plan_setup_mode is not None or
            gf8_avx2_walsh_locator_mode is not None or
            k65r65_b64_packed_terminal_mode is not None)
        field = "gf8" if diagnostic_avx2 or diagnostic_auto_gfni else "auto"
        backend = "avx2" if diagnostic_avx2 else "auto"
        command = [
            str(executable), "--k", str(k), "--r", str(r), "--profile", "high",
            "--field", field, "--backend", backend,
            "--bytes", str(shard_bytes),
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
        if measure_one_shot_encode:
            command.append("--measure-one-shot-encode")
        if disable_k16r8_b256_terminal:
            command.append("--disable-k16r8-b256-terminal")
        if disable_k9r5_b256_terminal:
            command.append("--disable-k9r5-b256-terminal")
        if disable_k9r5_b1024_terminal:
            command.append("--disable-k9r5-b1024-terminal")
        if disable_k9r6r8_b256_terminal:
            command.append("--disable-k9r6r8-b256-terminal")
        if r1_small_reduction_mode is not None:
            command.extend((
                "--r1-small-reduction-mode",
                str(r1_small_reduction_mode)))
        if r1_fixed_avx2_mode is not None:
            command.extend((
                "--r1-fixed-avx2-mode",
                str(r1_fixed_avx2_mode)))
        if prevalidated_binding:
            command.append("--prevalidated-binding")
        if one_shot_plan_setup_mode is not None:
            command.extend((
                "--one-shot-plan-setup-mode",
                str(one_shot_plan_setup_mode)))
        if gf8_avx2_walsh_locator_mode is not None:
            command.extend((
                "--gf8-avx2-walsh-locator-mode",
                str(gf8_avx2_walsh_locator_mode)))
        if k65r65_b64_packed_terminal_mode is not None:
            command.extend((
                "--k65r65-b64-packed-terminal-mode",
                str(k65r65_b64_packed_terminal_mode)))
        if k65r65_b64_avx512_gfni_mode is not None:
            command.extend((
                "--k65r65-b64-avx512-gfni-mode",
                str(k65r65_b64_avx512_gfni_mode)))
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
                "leopard2-benchmark-v10", "leopard2-benchmark-v11",
                "leopard2-benchmark-v12", "leopard2-benchmark-v15",
                "leopard2-benchmark-v16", "leopard2-benchmark-v23",
                "leopard2-benchmark-v31", "leopard2-benchmark-v32",
            }, "benchmark schema is unsupported")
    expected_top = {
        "schema", "build", "parameters", "resolved", "correctness",
        "memory", "metrics", "legacy"}
    if document["schema"] in {
        "leopard2-benchmark-v2", "leopard2-benchmark-v3",
        "leopard2-benchmark-v5", "leopard2-benchmark-v6",
        "leopard2-benchmark-v7", "leopard2-benchmark-v8",
        "leopard2-benchmark-v9", "leopard2-benchmark-v10",
        "leopard2-benchmark-v11", "leopard2-benchmark-v12",
        "leopard2-benchmark-v15", "leopard2-benchmark-v16",
        "leopard2-benchmark-v23", "leopard2-benchmark-v31",
        "leopard2-benchmark-v32",
    }:
        expected_top.add("workload_digests")
    require(set(document) == expected_top, "top-level JSON keys changed")
    for section in (
            "build", "parameters", "resolved", "correctness", "memory",
            "metrics", "legacy"):
        require(isinstance(document[section], dict),
                f"benchmark {section} is not an object")
    is_k65r65_v31 = document["schema"] == "leopard2-benchmark-v31"
    k65r65_v31_attested = (
        is_k65r65_v31 and
        document["parameters"].get("attest_source") is True)
    k65r65_v31_one_shot = (
        is_k65r65_v31 and
        document["parameters"].get("measure_one_shot_encode") is True)
    is_k65r65_v32 = document["schema"] == "leopard2-benchmark-v32"
    k65r65_v32_attested = (
        is_k65r65_v32 and
        document["parameters"].get("attest_source") is True)
    expected_build = {
        "compiler", "compiler_version", "cplusplus",
        "k8r3r4_t4_terminal_diagnostic_disabled",
        "t8_full_parity_terminal_diagnostic_disabled",
        "k16r8_b256_terminal_diagnostic_disabled",
        "k9r5_b256_terminal_diagnostic_disabled",
        "k9r5_b1024_terminal_diagnostic_disabled",
    }
    if "high_t16_q2_b64_fused_diagnostic_disabled" in document["build"]:
        expected_build.add("high_t16_q2_b64_fused_diagnostic_disabled")
    if document["schema"] == "leopard2-benchmark-v10":
        expected_build.add("k9r6r8_b256_terminal_diagnostic_disabled")
    if document["schema"] == "leopard2-benchmark-v12":
        expected_build.update({
            "r1_small_reduction_diagnostic_mode",
            "r1_small_reduction_codec_enabled",
            "r1_encode_reduction_path",
            "r1_decode_reduction_path",
            "r1_timed_encode_api",
            "r1_timed_one_shot_encode_api",
            "r1_timed_reused_decode_api",
            "r1_timed_one_shot_decode_api",
        })
    if document["schema"] == "leopard2-benchmark-v23":
        expected_build.update({
            "r1_fixed_avx2_diagnostic_mode",
            "r1_fixed_avx2_candidate_enabled",
            "r1_small_reduction_codec_enabled",
            "r1_fixed_avx2_selector_contract",
            "r1_encode_reduction_path",
            "r1_decode_reduction_path",
            "r1_timed_encode_api",
            "r1_timed_one_shot_encode_api",
            "r1_timed_reused_decode_api",
            "r1_timed_one_shot_decode_api",
        })
    if document["schema"] == "leopard2-benchmark-v31":
        expected_build.update({
            "k65r65_b64_packed_terminal_diagnostic_mode",
            "k65r65_b64_packed_terminal_diagnostic_disabled",
            "k65r65_b64_packed_terminal_mode_latched",
            "k65r65_b64_packed_terminal_selector_expected_selected",
            "k65r65_b64_packed_terminal_selector_selected",
            "k65r65_b64_packed_terminal_selector_contract",
            "k65r65_b64_packed_terminal_timed_ordinary_encode_api",
        })
        if k65r65_v31_one_shot:
            expected_build.add(
                "k65r65_b64_packed_terminal_timed_one_shot_encode_api")
    if is_k65r65_v32:
        expected_build.update({
            "k65r65_b64_avx512_gfni_diagnostic_mode",
            "k65r65_b64_avx512_gfni_diagnostic_disabled",
            "k65r65_b64_avx512_gfni_mode_latched",
            "k65r65_b64_avx512_gfni_kernel_qualified",
            "k65r65_b64_avx512_gfni_selector_expected_selected",
            "k65r65_b64_avx512_gfni_selector_selected",
            "k65r65_b64_avx512_gfni_observed_call_count",
            "k65r65_b64_avx512_gfni_selector_contract",
            "k65r65_b64_avx512_gfni_timed_ordinary_encode_api",
            "k65r65_b64_avx512_gfni_timed_one_shot_encode_api",
        })
    if document["schema"] in TRANSIENT_PLAN_SCHEMAS:
        expected_build.update({
            "one_shot_plan_setup_diagnostic_mode",
            "one_shot_plan_setup_policy",
            "one_shot_transient_execution_scope",
            "one_shot_transient_execution_excludes",
            "one_shot_transient_execution_guards",
            "one_shot_transient_execution_scratch_policy",
        })
    if document["schema"] == "leopard2-benchmark-v16":
        expected_build.update({
            "gf8_avx2_walsh_locator_diagnostic_mode",
            "gf8_avx2_walsh_locator_enabled",
        })
    if document["schema"] in {
            "leopard2-benchmark-v5", "leopard2-benchmark-v7"} or \
            k65r65_v31_attested or k65r65_v32_attested:
        expected_build.update({
            "source_commit", "source_tree", "source_tracked_dirty"})
    if document["schema"] in {
            "leopard2-benchmark-v6", "leopard2-benchmark-v7"}:
        expected_build.add("equal_rounded_multi_loss_enabled")
    if document["schema"] in ONE_SHOT_SCHEMAS:
        expected_build.add("one_shot_equal_rounded_direct_enabled")
    if document["schema"] in {
            "leopard2-benchmark-v9", "leopard2-benchmark-v12",
            "leopard2-benchmark-v15", "leopard2-benchmark-v16",
            "leopard2-benchmark-v23"}:
        expected_build.update({
            "one_shot_equal_rounded_direct_enabled",
            "cauchy_log_reuse_enabled",
        })
    require(set(document["build"]) == expected_build, "build keys changed")
    require(type(document["build"]["compiler"]) is str and
            type(document["build"]["compiler_version"]) is str and
            type(document["build"]["cplusplus"]) is int and
            type(document["build"][
                "k8r3r4_t4_terminal_diagnostic_disabled"]) is bool and
            type(document["build"][
                "t8_full_parity_terminal_diagnostic_disabled"]) is bool and
            type(document["build"][
                "k16r8_b256_terminal_diagnostic_disabled"]) is bool and
            type(document["build"][
                "k9r5_b256_terminal_diagnostic_disabled"]) is bool and
            type(document["build"][
                "k9r5_b1024_terminal_diagnostic_disabled"]) is bool and
            document["build"]["cplusplus"] > 0,
            "benchmark build value types changed")
    if "high_t16_q2_b64_fused_diagnostic_disabled" in document["build"]:
        require(type(document["build"][
                    "high_t16_q2_b64_fused_diagnostic_disabled"]) is bool,
                "high-T16 fused terminal selector is not Boolean")
    if document["schema"] == "leopard2-benchmark-v12":
        build = document["build"]
        mode = build["r1_small_reduction_diagnostic_mode"]
        require(type(mode) is int and mode in {0, 1} and
                type(build["r1_small_reduction_codec_enabled"]) is bool and
                build["r1_small_reduction_codec_enabled"] is (mode == 1),
                "R=1 small-reduction mode metadata is invalid")
    if document["schema"] == "leopard2-benchmark-v23":
        build = document["build"]
        mode = build["r1_fixed_avx2_diagnostic_mode"]
        require(type(mode) is int and mode in {0, 1} and
                build["r1_fixed_avx2_candidate_enabled"] is True and
                build["r1_small_reduction_codec_enabled"] is False and
                build["r1_fixed_avx2_selector_contract"] ==
                    R1_FIXED_AVX2_SELECTOR_CONTRACT,
                "fixed AVX2 R=1 mode metadata is invalid")
    if document["schema"] == "leopard2-benchmark-v31":
        build = document["build"]
        mode = build["k65r65_b64_packed_terminal_diagnostic_mode"]
        selector_expected = mode == 1 and \
            document["parameters"]["K"] == 65 and \
            document["parameters"]["R"] == 65 and \
            document["parameters"]["shard_bytes"] == 64
        require(type(mode) is int and mode in {0, 1} and
                type(build[
                    "k65r65_b64_packed_terminal_diagnostic_disabled"]) is
                    bool and
                build[
                    "k65r65_b64_packed_terminal_diagnostic_disabled"] is
                    (mode == 0) and
                type(build[
                    "k65r65_b64_packed_terminal_mode_latched"]) is int and
                build["k65r65_b64_packed_terminal_mode_latched"] == mode and
                type(build[
                    "k65r65_b64_packed_terminal_selector_expected_selected"])
                    is bool and
                type(build[
                    "k65r65_b64_packed_terminal_selector_selected"]) is bool and
                build[
                    "k65r65_b64_packed_terminal_selector_expected_selected"]
                    is selector_expected and
                build["k65r65_b64_packed_terminal_selector_selected"] is
                    selector_expected and
                build["k65r65_b64_packed_terminal_selector_contract"] ==
                    K65R65_B64_SELECTOR_CONTRACT and
                build[
                    "k65r65_b64_packed_terminal_timed_ordinary_encode_api"] ==
                    K65R65_B64_TIMED_ORDINARY_ENCODE_API and
                (not k65r65_v31_one_shot or
                 build[
                    "k65r65_b64_packed_terminal_timed_one_shot_encode_api"] ==
                    K65R65_B64_TIMED_ONE_SHOT_ENCODE_API),
                "K65/R65/B64 packed-terminal mode metadata is invalid")
    if is_k65r65_v32:
        build = document["build"]
        mode = build["k65r65_b64_avx512_gfni_diagnostic_mode"]
        qualified = build["k65r65_b64_avx512_gfni_kernel_qualified"]
        selector_expected = (
            mode == 1 and qualified and
            document["parameters"]["K"] == 65 and
            document["parameters"]["R"] == 65 and
            document["parameters"]["shard_bytes"] == 64 and
            document["resolved"]["backend"] == "avx2")
        require(type(mode) is int and mode in {0, 1} and
                type(build[
                    "k65r65_b64_avx512_gfni_diagnostic_disabled"]) is bool and
                build[
                    "k65r65_b64_avx512_gfni_diagnostic_disabled"] is
                    (mode == 0) and
                type(build["k65r65_b64_avx512_gfni_mode_latched"]) is int and
                build["k65r65_b64_avx512_gfni_mode_latched"] == mode and
                type(qualified) is bool and
                type(build[
                    "k65r65_b64_avx512_gfni_selector_expected_selected"])
                    is bool and
                type(build[
                    "k65r65_b64_avx512_gfni_selector_selected"]) is bool and
                build[
                    "k65r65_b64_avx512_gfni_selector_expected_selected"] is
                    selector_expected and
                build["k65r65_b64_avx512_gfni_selector_selected"] is
                    selector_expected and
                type(build[
                    "k65r65_b64_avx512_gfni_observed_call_count"]) is int and
                build["k65r65_b64_avx512_gfni_observed_call_count"] ==
                    (2 if selector_expected else 0) and
                build["k65r65_b64_avx512_gfni_selector_contract"] ==
                    K65R65_B64_AVX512_GFNI_SELECTOR_CONTRACT and
                build[
                    "k65r65_b64_avx512_gfni_timed_ordinary_encode_api"] ==
                    K65R65_B64_TIMED_ORDINARY_ENCODE_API and
                build[
                    "k65r65_b64_avx512_gfni_timed_one_shot_encode_api"] ==
                    K65R65_B64_TIMED_ONE_SHOT_ENCODE_API,
                "K65/R65/B64 AVX-512/GFNI mode metadata is invalid")
    if document["schema"] in R1_SCHEMAS:
        build = document["build"]
        for name in (
                "r1_encode_reduction_path", "r1_decode_reduction_path"):
            require(type(build[name]) is str and
                    build[name] in R1_REDUCTION_PATHS,
                    f"R=1 reduction path {name} is invalid")
        require(build["r1_timed_encode_api"] == R1_TIMED_ENCODE_API and
                build["r1_timed_one_shot_encode_api"] ==
                    R1_TIMED_ONE_SHOT_ENCODE_API and
                build["r1_timed_reused_decode_api"] ==
                    R1_TIMED_REUSED_DECODE_API and
                build["r1_timed_one_shot_decode_api"] ==
                    R1_TIMED_ONE_SHOT_DECODE_API,
                "R=1 timed API scope changed")
    if document["schema"] in TRANSIENT_PLAN_SCHEMAS:
        build = document["build"]
        mode = build["one_shot_plan_setup_diagnostic_mode"]
        require(type(mode) is int and mode in TRANSIENT_SETUP_POLICIES and
                build["one_shot_plan_setup_policy"] ==
                    TRANSIENT_SETUP_POLICIES[mode] and
                build["one_shot_transient_execution_scope"] ==
                    "diagnostic_transform_plan_execution_only" and
                build["one_shot_transient_execution_excludes"] ==
                    TRANSIENT_EXECUTION_EXCLUDES and
                build["one_shot_transient_execution_guards"] ==
                    TRANSIENT_EXECUTION_GUARDS and
                build["one_shot_transient_execution_scratch_policy"] ==
                    "exact_plan_query",
                "one-shot transient attribution metadata is invalid")
    if document["schema"] == "leopard2-benchmark-v16":
        build = document["build"]
        mode = build["gf8_avx2_walsh_locator_diagnostic_mode"]
        require(type(mode) is int and mode in {0, 1} and
                type(build["gf8_avx2_walsh_locator_enabled"]) is bool and
                build["gf8_avx2_walsh_locator_enabled"] is (mode == 1),
                "GF8 AVX2 Walsh-locator attribution metadata is invalid")
    if document["schema"] == "leopard2-benchmark-v10":
        require(type(document["build"][
                    "k9r6r8_b256_terminal_diagnostic_disabled"]) is bool,
                "K9/R=6..8 terminal selector is not Boolean")
    if document["schema"] in {
            "leopard2-benchmark-v6", "leopard2-benchmark-v7"}:
        require(type(document["build"][
                    "equal_rounded_multi_loss_enabled"]) is bool,
                "equal-rounded build selector is not Boolean")
    if document["schema"] in ONE_SHOT_SCHEMAS:
        require(type(document["build"][
                    "one_shot_equal_rounded_direct_enabled"]) is bool,
                "one-shot build selector is not Boolean")
    if document["schema"] in {
            "leopard2-benchmark-v9", "leopard2-benchmark-v12",
            "leopard2-benchmark-v15", "leopard2-benchmark-v16",
            "leopard2-benchmark-v23"}:
        require(type(document["build"][
                    "cauchy_log_reuse_enabled"]) is bool,
                "Cauchy-log-reuse build selector is not Boolean")
    if document["schema"] in {
            "leopard2-benchmark-v5", "leopard2-benchmark-v7"} or \
            k65r65_v31_attested or k65r65_v32_attested:
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
    if document["schema"] in TRANSIENT_PLAN_SCHEMAS:
        expected_resolved.update({
            "one_shot_transient_decode_path",
            "one_shot_transient_decode_rule",
            "one_shot_transient_low_input_pruned_plans",
            "one_shot_transient_low_output_pruned_plans",
            "one_shot_transient_high_input_pruned_plans",
            "one_shot_transient_high_output_pruned_plans",
            "one_shot_transient_setup_scope",
            "one_shot_transient_setup_excludes",
            "one_shot_transient_setup_includes_no_loss_prefix_probe",
        })
    if document["schema"] == "leopard2-benchmark-v16":
        expected_resolved.add("gf8_avx2_walsh_locator_enabled")
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
    if document["schema"] in TRANSIENT_PLAN_SCHEMAS:
        mode = parameters.get("one_shot_plan_setup_mode")
        require(type(mode) is int and mode in TRANSIENT_SETUP_POLICIES and
                mode == document["build"][
                    "one_shot_plan_setup_diagnostic_mode"],
                "one-shot setup parameter is invalid or inconsistent")
    if document["schema"] == "leopard2-benchmark-v16":
        mode = parameters.get("gf8_avx2_walsh_locator_mode")
        require(type(mode) is int and mode in {0, 1} and
                mode == document["build"][
                    "gf8_avx2_walsh_locator_diagnostic_mode"],
                "GF8 AVX2 Walsh-locator parameter is invalid or inconsistent")
    if document["schema"] == "leopard2-benchmark-v23":
        mode = parameters.get("r1_fixed_avx2_mode")
        require(type(mode) is int and mode in {0, 1} and
                mode == document["build"]["r1_fixed_avx2_diagnostic_mode"],
                "fixed AVX2 R=1 parameter is invalid or inconsistent")
    if document["schema"] == "leopard2-benchmark-v31":
        expected_parameters = {
            "K", "R", "requested_profile", "requested_field",
            "requested_backend", "force_generic_decode",
            "force_specialized_decode", "force_tiled_decode",
            "force_materialized_decode", "skip_legacy", "retain_samples",
            "k65r65_b64_packed_terminal_mode", "shard_bytes", "loss_count",
            "missing_original_indices", "batch", "reuse", "iterations",
            "warmup", "thread_count", "seed",
        }
        if k65r65_v31_attested:
            expected_parameters.add("attest_source")
        if k65r65_v31_one_shot:
            expected_parameters.add("measure_one_shot_encode")
        mode = parameters.get("k65r65_b64_packed_terminal_mode")
        require(set(parameters) == expected_parameters and
                type(mode) is int and mode in {0, 1} and
                mode == document["build"][
                    "k65r65_b64_packed_terminal_diagnostic_mode"] and
                parameters["requested_profile"] == "legacy_high_v1" and
                parameters["requested_field"] == "gf8" and
                parameters["requested_backend"] == "avx2" and
                parameters["batch"] == 1 and
                parameters["thread_count"] == 1 and
                parameters["skip_legacy"] is True and
                parameters["retain_samples"] is True and
                all(parameters[name] is False for name in (
                    "force_generic_decode", "force_specialized_decode",
                    "force_tiled_decode", "force_materialized_decode")) and
                type(parameters["warmup"]) is int and
                parameters["warmup"] >= 0 and
                type(parameters["reuse"]) is int and
                parameters["reuse"] > 0 and
                type(parameters["iterations"]) is int and
                parameters["iterations"] > 0 and
                type(parameters["seed"]) is int and
                0 <= parameters["seed"] <= (1 << 64) - 1 and
                (not k65r65_v31_attested or
                 parameters["attest_source"] is True) and
                (not k65r65_v31_one_shot or
                 parameters["measure_one_shot_encode"] is True),
                "K65/R65/B64 packed-terminal parameters are invalid or "
                "ambiguous")
    if is_k65r65_v32:
        expected_parameters = {
            "K", "R", "requested_profile", "requested_field",
            "requested_backend", "force_generic_decode",
            "force_specialized_decode", "force_tiled_decode",
            "force_materialized_decode", "skip_legacy", "retain_samples",
            "measure_one_shot_encode",
            "k65r65_b64_avx512_gfni_mode", "shard_bytes", "loss_count",
            "missing_original_indices", "batch", "reuse", "iterations",
            "warmup", "thread_count", "seed",
        }
        if k65r65_v32_attested:
            expected_parameters.add("attest_source")
        mode = parameters.get("k65r65_b64_avx512_gfni_mode")
        require(set(parameters) == expected_parameters and
                type(mode) is int and mode in {0, 1} and
                mode == document["build"][
                    "k65r65_b64_avx512_gfni_diagnostic_mode"] and
                parameters["requested_profile"] == "legacy_high_v1" and
                parameters["requested_field"] == "gf8" and
                parameters["requested_backend"] == "auto" and
                parameters["batch"] == 1 and
                parameters["thread_count"] == 1 and
                parameters["skip_legacy"] is True and
                parameters["retain_samples"] is True and
                parameters["measure_one_shot_encode"] is True and
                all(parameters[name] is False for name in (
                    "force_generic_decode", "force_specialized_decode",
                    "force_tiled_decode", "force_materialized_decode")) and
                type(parameters["warmup"]) is int and
                parameters["warmup"] >= 0 and
                type(parameters["reuse"]) is int and
                parameters["reuse"] > 0 and
                type(parameters["iterations"]) is int and
                parameters["iterations"] > 0 and
                type(parameters["seed"]) is int and
                0 <= parameters["seed"] <= (1 << 64) - 1 and
                (not k65r65_v32_attested or
                 parameters["attest_source"] is True),
                "K65/R65/B64 AVX-512/GFNI parameters are invalid or "
                "ambiguous")
    require(type(document["resolved"]["thread_count"]) is int and
            document["resolved"]["thread_count"] > 0 and
            type(document["resolved"]["parent_count"]) is int and
            document["resolved"]["parent_count"] > 0 and
            type(document["resolved"]["padded_side"]) is int and
            document["resolved"]["padded_side"] > 0,
            "resolved numeric codec identity is invalid")
    if document["schema"] == "leopard2-benchmark-v31":
        resolved = document["resolved"]
        require(resolved["profile"] == "legacy_high_v1" and
                resolved["field"] == "gf8" and
                resolved["backend"] == "avx2" and
                type(resolved["thread_count"]) is int and
                resolved["thread_count"] == 1,
                "K65/R65/B64 resolved codec identity is invalid")
    if is_k65r65_v32:
        resolved = document["resolved"]
        require(resolved["profile"] == "legacy_high_v1" and
                resolved["field"] == "gf8" and
                resolved["backend"] in {
                    "scalar", "ssse3", "avx2", "neon", "avx512",
                    "avx2-gfni"} and
                type(resolved["thread_count"]) is int and
                resolved["thread_count"] == 1,
                "K65/R65/B64 AVX-512/GFNI resolved codec identity is "
                "invalid")
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
    if document["schema"] in TRANSIENT_PLAN_SCHEMAS:
        resolved = document["resolved"]
        transient_pair = (
            resolved["one_shot_transient_decode_path"],
            resolved["one_shot_transient_decode_rule"])
        require(all(type(item) is str for item in transient_pair) and
                transient_pair in PRODUCTION_DECODE_PATH_RULE_PAIRS and
                all(type(resolved[name]) is int and resolved[name] >= 0
                    for name in (
                        "one_shot_transient_low_input_pruned_plans",
                        "one_shot_transient_low_output_pruned_plans",
                        "one_shot_transient_high_input_pruned_plans",
                        "one_shot_transient_high_output_pruned_plans")) and
                resolved["one_shot_transient_setup_scope"] ==
                    "diagnostic_pattern_setup" and
                resolved["one_shot_transient_setup_excludes"] ==
                    TRANSIENT_SETUP_EXCLUDES and
                type(resolved[
                    "one_shot_transient_setup_includes_no_loss_prefix_probe"])
                    is bool,
                "one-shot transient resolved attribution is invalid")
    if document["schema"] == "leopard2-benchmark-v16":
        require(type(document["resolved"][
                    "gf8_avx2_walsh_locator_enabled"]) is bool and
                document["resolved"][
                    "gf8_avx2_walsh_locator_enabled"] is
                    document["build"][
                        "gf8_avx2_walsh_locator_enabled"],
                "resolved GF8 AVX2 Walsh-locator mode is inconsistent")
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
    if document["schema"] in ONE_SHOT_SCHEMAS:
        expected_memory.update({
            "one_shot_decode_scratch_bytes_per_stripe",
            "one_shot_decode_scratch_bytes_batch",
        })
    if document["schema"] in TRANSIENT_PLAN_SCHEMAS:
        expected_memory.add("one_shot_transient_scratch_bytes_per_stripe")
    if document["schema"] in R1_SCHEMAS:
        expected_memory.update({
            "encode_batch_preflight_scratch_bytes",
            "decode_batch_preflight_scratch_bytes",
        })
    require(set(document["memory"]) == expected_memory,
            "memory keys changed")
    require(all(type(value) is int and value >= 0
                for value in document["memory"].values()) and
            document["memory"]["scratch_alignment"] > 0,
            "memory values are invalid")
    if document["schema"] in TRANSIENT_PLAN_SCHEMAS:
        require(document["memory"][
                    "one_shot_transient_scratch_bytes_per_stripe"] <=
                document["memory"][
                    "one_shot_decode_scratch_bytes_per_stripe"],
                "transient plan scratch exceeds the public one-shot bound")
    expected_metrics = {
        "codec_setup", "encode_execution", "decode_plan_setup",
        "decode_execution", "decode_amortized_at_reuse", "rate_semantics"}
    if document["schema"] in ONE_SHOT_SCHEMAS:
        expected_metrics.add("one_shot_decode_including_setup")
    if document["schema"] in R1_SCHEMAS:
        expected_metrics.add("one_shot_encode")
    if k65r65_v31_one_shot or is_k65r65_v32:
        expected_metrics.add("one_shot_encode")
    if document["schema"] in TRANSIENT_PLAN_SCHEMAS:
        expected_metrics.update({
            "one_shot_transient_plan_setup",
            "one_shot_transient_execution",
            "one_shot_transient_amortized_at_reuse",
        })
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
    if (document["schema"] in R1_SCHEMAS or k65r65_v31_one_shot or
            is_k65r65_v32):
        validate_timing_summary(
            document["metrics"]["one_shot_encode"], "one_shot_encode",
            retain_samples=retain_samples, iterations=iterations,
            execution=True, input_rate_name="input_GB_per_s",
            output_rate_name="parity_output_GB_per_s",
            input_bytes=encode_input_bytes, output_bytes=encode_output_bytes)
    validate_timing_summary(
        document["metrics"]["decode_execution"], "decode_execution",
        retain_samples=retain_samples, iterations=iterations, execution=True,
        input_rate_name="offered_received_GB_per_s",
        output_rate_name="repaired_output_GB_per_s",
        input_bytes=decode_input_bytes, output_bytes=decode_output_bytes)
    if document["schema"] in ONE_SHOT_SCHEMAS:
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
    if document["schema"] in TRANSIENT_PLAN_SCHEMAS:
        validate_timing_summary(
            document["metrics"]["one_shot_transient_plan_setup"],
            "one_shot_transient_plan_setup",
            retain_samples=retain_samples, iterations=iterations,
            execution=False)
        validate_timing_summary(
            document["metrics"]["one_shot_transient_execution"],
            "one_shot_transient_execution",
            retain_samples=retain_samples, iterations=iterations,
            execution=True,
            input_rate_name="offered_received_GB_per_s",
            output_rate_name="repaired_output_GB_per_s",
            input_bytes=decode_input_bytes,
            output_bytes=decode_output_bytes)
        transient_amortized = document["metrics"][
            "one_shot_transient_amortized_at_reuse"]
        require(isinstance(transient_amortized, dict) and
                set(transient_amortized) == {
                    "reuse_count", "derived_median_us_per_batch_call",
                    "offered_received_GB_per_s", "repaired_output_GB_per_s",
                } and type(transient_amortized["reuse_count"]) is int and
                transient_amortized["reuse_count"] == parameters["reuse"],
                "transient plan amortization structure changed")
        transient_derived = finite_metric(
            transient_amortized["derived_median_us_per_batch_call"],
            "transient plan amortized median", positive=True)
        transient_expected = (
            float(document["metrics"]["one_shot_transient_execution"][
                "median_us_per_batch_call"]) +
            float(document["metrics"]["one_shot_transient_plan_setup"][
                "median_us"]) / parameters["reuse"])
        require(approximately_equal(transient_derived, transient_expected),
                "transient plan amortization is inconsistent with its medians")
        for name, byte_count in (
                ("offered_received_GB_per_s", decode_input_bytes),
                ("repaired_output_GB_per_s", decode_output_bytes)):
            if byte_count == 0:
                require(transient_amortized[name] is None,
                        f"transient amortized {name} is non-null for zero bytes")
            else:
                finite_metric(
                    transient_amortized[name],
                    f"transient amortized {name}", positive=True)

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
    expected = {
        "algorithm", "original_data", "transmitted_parity",
        "recovered_originals"}
    if document.get("schema") in TRANSIENT_PLAN_SCHEMAS:
        expected.add("recovered_originals_provenance")
    require(isinstance(digests, dict) and set(digests) == expected,
            "workload digest structure changed")
    require(digests["algorithm"] == "fnv1a64", "workload digest algorithm changed")
    for name in ("original_data", "transmitted_parity", "recovered_originals"):
        value = digests[name]
        require(isinstance(value, str) and len(value) == 16 and
                all(character in "0123456789abcdef" for character in value),
                f"workload digest {name} is not lowercase FNV-1a hex")
    if document.get("schema") in TRANSIENT_PLAN_SCHEMAS:
        require(digests["recovered_originals_provenance"] ==
                "diagnostic_one_shot_transform_execution_after_poison",
                "transient recovered digest provenance changed")


def validate_r1_small_reduction_report(
    document: dict[str, Any],
    expected_mode: int,
    expected_parameters: set[str],
    expected_k: int,
    expected_path: str,
) -> None:
    require(type(expected_mode) is int and expected_mode in {0, 1},
            "invalid expected R=1 small-reduction mode")
    require(type(expected_k) is int and expected_k >= 2 and
            expected_path in {"pairwise", "k2_terminal", "fused_final"},
            "invalid expected R=1 small-reduction path")
    require(document.get("schema") == "leopard2-benchmark-v12",
            "R=1 small-reduction benchmark schema changed")
    validate_common(document, True)
    validate_workload_digests(document)

    build = document["build"]
    require(build["r1_small_reduction_diagnostic_mode"] == expected_mode and
            build["r1_small_reduction_codec_enabled"] is
                (expected_mode == 1) and
            build["r1_encode_reduction_path"] == expected_path and
            build["r1_decode_reduction_path"] == expected_path,
            "R=1 reduction attribution selected an unexpected path")

    parameters = document["parameters"]
    require(set(parameters) == expected_parameters and
            parameters["K"] == expected_k and parameters["R"] == 1 and
            parameters["requested_profile"] == "legacy_high_v1" and
            parameters["requested_field"] == "gf8" and
            parameters["requested_backend"] == "avx2" and
            parameters["shard_bytes"] == 64 and
            parameters["loss_count"] == 1 and
            parameters["batch"] == 1 and parameters["reuse"] == 1 and
            parameters["iterations"] == 1 and
            parameters["warmup"] == 0 and
            parameters["thread_count"] == 1 and
            parameters["skip_legacy"] is True and
            parameters["retain_samples"] is True and
            parameters["measure_one_shot_decode"] is True,
            "R=1 small-reduction benchmark is not the ordinary one-item cell")

    require(document["resolved"] == {
        "profile": "legacy_high_v1",
        "field": "gf8",
        "backend": "avx2",
        "thread_count": 1,
        "parent_count": 4,
        "padded_side": 1,
    }, "R=1 small-reduction resolved codec identity changed")
    memory = document["memory"]
    require(memory["encode_batch_preflight_scratch_bytes"] == 0 and
            memory["decode_batch_preflight_scratch_bytes"] == 0 and
            memory["one_shot_decode_scratch_bytes_batch"] ==
                memory["one_shot_decode_scratch_bytes_per_stripe"],
            "R=1 benchmark did not retain the no-preflight one-item API scope")
    require(document["legacy"] == {
        "available": False,
        "unavailable_reason": "disabled by --skip-legacy",
        "codec_setup": None,
        "decode_timing_includes_setup": True,
        "encode_execution": None,
        "decode_including_setup": None,
    }, "R=1 attribution did not completely exclude the legacy timing lane")


def validate_r1_fixed_avx2_report(
    document: dict[str, Any],
    expected_mode: int,
    expected_parameters: set[str],
    expected_path: str,
) -> None:
    require(type(expected_mode) is int and expected_mode in {0, 1},
            "invalid expected fixed AVX2 R=1 mode")
    require(expected_path in {"pairwise", "fixed_avx2"},
            "invalid expected fixed AVX2 R=1 path")
    require(document.get("schema") == "leopard2-benchmark-v23",
            "fixed AVX2 R=1 benchmark schema changed")
    validate_common(document, True)
    validate_workload_digests(document)

    build = document["build"]
    require(build["r1_fixed_avx2_diagnostic_mode"] == expected_mode and
            build["r1_fixed_avx2_candidate_enabled"] is True and
            build["r1_small_reduction_codec_enabled"] is False and
            build["r1_encode_reduction_path"] == expected_path and
            build["r1_decode_reduction_path"] == expected_path,
            "fixed AVX2 R=1 attribution selected an unexpected path")

    parameters = document["parameters"]
    require(set(parameters) == expected_parameters and
            parameters["K"] == 3 and parameters["R"] == 1 and
            parameters["requested_profile"] == "legacy_high_v1" and
            parameters["requested_field"] == "gf8" and
            parameters["requested_backend"] == "avx2" and
            parameters["shard_bytes"] == 64 and
            parameters["loss_count"] == 1 and
            parameters["batch"] == 1 and parameters["reuse"] == 1 and
            parameters["iterations"] == 1 and
            parameters["warmup"] == 0 and
            parameters["thread_count"] == 1 and
            parameters["skip_legacy"] is True and
            parameters["retain_samples"] is True and
            parameters["measure_one_shot_decode"] is True and
            parameters["r1_fixed_avx2_mode"] == expected_mode,
            "fixed AVX2 R=1 benchmark is not the ordinary one-item cell")

    require(document["resolved"] == {
        "profile": "legacy_high_v1",
        "field": "gf8",
        "backend": "avx2",
        "thread_count": 1,
        "parent_count": 4,
        "padded_side": 1,
    }, "fixed AVX2 R=1 resolved codec identity changed")
    memory = document["memory"]
    require(memory["encode_batch_preflight_scratch_bytes"] == 0 and
            memory["decode_batch_preflight_scratch_bytes"] == 0 and
            memory["one_shot_decode_scratch_bytes_batch"] ==
                memory["one_shot_decode_scratch_bytes_per_stripe"],
            "fixed AVX2 R=1 report changed the one-item API scope")
    require(document["legacy"] == {
        "available": False,
        "unavailable_reason": "disabled by --skip-legacy",
        "codec_setup": None,
        "decode_timing_includes_setup": True,
        "encode_execution": None,
        "decode_including_setup": None,
    }, "fixed AVX2 R=1 attribution did not exclude the legacy lane")


def validate_prevalidated_r1_fixed_avx2_report(
    document: dict[str, Any],
    expected_mode: int,
    expected_parameters: set[str],
    expected_path: str,
) -> None:
    require(type(expected_mode) is int and expected_mode in {0, 1},
            "invalid expected prevalidated fixed AVX2 R=1 mode")
    require(expected_path in {"pairwise", "fixed_avx2"},
            "invalid expected prevalidated fixed AVX2 R=1 path")
    validate_prevalidated_batch_variant(document, True)
    validate_workload_digests(document)

    build = document["build"]
    require(build["r1_fixed_avx2_diagnostic_mode"] == expected_mode and
            build["r1_fixed_avx2_candidate_enabled"] is True and
            build["r1_small_reduction_codec_enabled"] is False and
            build["r1_encode_reduction_path"] == expected_path and
            build["r1_decode_reduction_path"] == expected_path,
            "prevalidated fixed AVX2 R=1 selected an unexpected path")

    parameters = document["parameters"]
    require(set(parameters) == expected_parameters and
            parameters["K"] == 3 and parameters["R"] == 1 and
            parameters["requested_profile"] == "legacy_high_v1" and
            parameters["requested_field"] == "gf8" and
            parameters["requested_backend"] == "avx2" and
            parameters["shard_bytes"] == 64 and
            parameters["loss_count"] == 1 and
            parameters["batch"] == 1 and parameters["reuse"] == 1 and
            parameters["iterations"] == 1 and
            parameters["warmup"] == 0 and
            parameters["thread_count"] == 1 and
            parameters["skip_legacy"] is True and
            parameters["retain_samples"] is True and
            parameters["measure_one_shot_decode"] is True and
            parameters["r1_fixed_avx2_mode"] == expected_mode and
            parameters["prevalidated_binding"] is True,
            "prevalidated fixed AVX2 R=1 report changed its one-item cell")
    require(document["resolved"] == {
        "profile": "legacy_high_v1",
        "field": "gf8",
        "backend": "avx2",
        "thread_count": 1,
        "parent_count": 4,
        "padded_side": 1,
    }, "prevalidated fixed AVX2 R=1 codec identity changed")
    memory = document["memory"]
    require(memory["encode_batch_preflight_scratch_bytes"] == 0 and
            memory["decode_batch_preflight_scratch_bytes"] == 0 and
            memory["one_shot_decode_scratch_bytes_batch"] ==
                memory["one_shot_decode_scratch_bytes_per_stripe"],
            "prevalidated fixed AVX2 R=1 report changed API scratch scope")
    require(document["legacy"] == {
        "available": False,
        "unavailable_reason": "disabled by --skip-legacy",
        "codec_setup": None,
        "decode_timing_includes_setup": True,
        "encode_execution": None,
        "decode_including_setup": None,
    }, "prevalidated fixed AVX2 R=1 report included the legacy lane")


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


def require_workload_digests_rejected(
    document: dict[str, Any],
    label: str,
) -> None:
    try:
        validate_workload_digests(document)
    except (KeyError, RuntimeError, TypeError, ValueError):
        return
    raise RuntimeError(f"workload-digest validator accepted {label}")


def validate_transient_plan_report(
    document: dict[str, Any],
    expected_parameters: set[str],
    *,
    setup_mode: int,
    walsh_mode: int | None,
) -> None:
    expected_schema = (
        "leopard2-benchmark-v16" if walsh_mode is not None else
        "leopard2-benchmark-v15")
    require(document.get("schema") == expected_schema,
            "transient-plan attribution schema changed")
    validate_common(document, True)
    validate_workload_digests(document)
    parameters = document["parameters"]
    require(set(parameters) == expected_parameters and
            parameters["requested_profile"] == "legacy_high_v1" and
            parameters["requested_field"] == "gf8" and
            parameters["requested_backend"] == "avx2" and
            parameters["batch"] == 1 and parameters["thread_count"] == 1 and
            parameters["skip_legacy"] is True and
            parameters["retain_samples"] is True and
            parameters["measure_one_shot_decode"] is True and
            parameters["one_shot_plan_setup_mode"] == setup_mode,
            "transient-plan attribution parameters changed")
    require(document["resolved"]["profile"] == "legacy_high_v1" and
            document["resolved"]["field"] == "gf8" and
            document["resolved"]["backend"] == "avx2" and
            document["resolved"]["one_shot_transient_decode_path"] not in {
                "no_op", "direct"},
            "transient-plan attribution did not select an AVX2 transform")
    require(document["legacy"] == {
        "available": False,
        "unavailable_reason": "disabled by --skip-legacy",
        "codec_setup": None,
        "decode_timing_includes_setup": True,
        "encode_execution": None,
        "decode_including_setup": None,
    }, "transient-plan attribution did not exclude the legacy lane")
    if walsh_mode is None:
        require("gf8_avx2_walsh_locator_mode" not in parameters,
                "v15 unexpectedly records a Walsh-locator mode")
    else:
        require(parameters["gf8_avx2_walsh_locator_mode"] == walsh_mode and
                document["resolved"][
                    "gf8_avx2_walsh_locator_enabled"] is (walsh_mode == 1),
                "v16 did not resolve the requested Walsh-locator mode")


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
    r1_binding_fields = {
        "r1_prevalidated_binding_attribution",
        "r1_encode_binding_setup_api",
        "r1_decode_binding_setup_api",
        "r1_binding_setup_reported_separately",
        "r1_binding_item_count",
    }
    require(isinstance(document.get("build"), dict) and
            build_fields <= set(document["build"]),
            "prevalidated benchmark build markers are incomplete")
    require(isinstance(document.get("metrics"), dict) and
            metric_fields <= set(document["metrics"]),
            "prevalidated benchmark setup metrics are incomplete")
    require(isinstance(document.get("parameters"), dict),
            "prevalidated benchmark parameters are missing")
    r1_binding = "prevalidated_binding" in document["parameters"]
    if r1_binding:
        build = document["build"]
        require(document.get("schema") == "leopard2-benchmark-v23" and
                document["parameters"]["prevalidated_binding"] is True and
                r1_binding_fields <= set(build) and
                build["r1_prevalidated_binding_attribution"] is True and
                build["r1_encode_binding_setup_api"] ==
                    R1_ENCODE_BINDING_SETUP_API and
                build["r1_decode_binding_setup_api"] ==
                    R1_DECODE_BINDING_SETUP_API and
                build["r1_binding_setup_reported_separately"] is True and
                type(build["r1_binding_item_count"]) is int and
                build["r1_binding_item_count"] == 1 and
                build["r1_timed_encode_api"] ==
                    R1_PREVALIDATED_TIMED_ENCODE_API and
                build["r1_timed_one_shot_encode_api"] ==
                    R1_TIMED_ONE_SHOT_ENCODE_API and
                build["r1_timed_reused_decode_api"] ==
                    R1_PREVALIDATED_TIMED_REUSED_DECODE_API and
                build["r1_timed_one_shot_decode_api"] ==
                    R1_TIMED_ONE_SHOT_DECODE_API,
                "prevalidated R=1 binding attribution is invalid")

    normalized = copy.deepcopy(document)
    for name in build_fields:
        del normalized["build"][name]
    for name in metric_fields:
        del normalized["metrics"][name]
    if r1_binding:
        for name in r1_binding_fields:
            del normalized["build"][name]
        del normalized["parameters"]["prevalidated_binding"]
        normalized["build"]["r1_timed_encode_api"] = R1_TIMED_ENCODE_API
        normalized["build"]["r1_timed_reused_decode_api"] = (
            R1_TIMED_REUSED_DECODE_API)
    validate_common(normalized, retain_samples)

    build = document["build"]
    require(build["prevalidated_batch_experiment"] is True and
            all(type(build[name]) is bool for name in build_fields),
            "prevalidated benchmark build markers are not Boolean")
    require(build["high_t4_batch_selected"] is False and
            build["high_t8_one_block_selected"] is False and
            build["high_t8_two_block_selected"] is False,
            "R=1 unexpectedly selected a T=4/T=8 encode binding")

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
                "t8_full_parity_terminal_diagnostic_disabled"] is False,
            "default benchmark disabled the T=8 full-parity terminal")
    require(default["build"][
                "k16r8_b256_terminal_diagnostic_disabled"] is False,
            "default benchmark disabled the K16/R8 terminal")
    require(default["build"][
                "k9r5_b256_terminal_diagnostic_disabled"] is False,
            "default benchmark disabled the K9/R5 terminal")
    require(default["build"][
                "k9r5_b1024_terminal_diagnostic_disabled"] is False,
            "default benchmark disabled the K9/R5/1024 terminal")
    require(all(not name.startswith("k65r65_b64_packed_terminal_")
                for name in default["build"]) and
            "k65r65_b64_packed_terminal_mode" not in default["parameters"],
            "default benchmark exposed the opt-in K65/R65/B64 diagnostic")
    require(all(not name.startswith("k65r65_b64_avx512_gfni_")
                for name in default["build"]) and
            "k65r65_b64_avx512_gfni_mode" not in default["parameters"],
            "default benchmark exposed the opt-in K65/R65/B64 AVX-512/"
            "GFNI diagnostic")
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
    require(external["build"][
                "k9r5_b1024_terminal_diagnostic_disabled"] is False,
            "external-evidence benchmark disabled the K9/R5/1024 terminal")
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

    k65r65_mode_zero = run(
        executable, True, k65r65_b64_packed_terminal_mode=0,
        k=65, r=65, losses=1)
    k65r65_mode_one = run(
        executable, True, k65r65_b64_packed_terminal_mode=1,
        k=65, r=65, losses=1)
    for mode, document in enumerate((k65r65_mode_zero, k65r65_mode_one)):
        require(document["schema"] == "leopard2-benchmark-v31",
                "K65/R65/B64 diagnostic schema changed")
        validate_common(document, True)
        validate_workload_digests(document)
        require(document["build"][
                    "k65r65_b64_packed_terminal_diagnostic_mode"] == mode and
                document["build"][
                    "k65r65_b64_packed_terminal_diagnostic_disabled"] is
                    (mode == 0) and
                document["build"][
                    "k65r65_b64_packed_terminal_mode_latched"] == mode and
                document["build"][
                    "k65r65_b64_packed_terminal_selector_expected_selected"]
                    is (mode == 1) and
                document["build"][
                    "k65r65_b64_packed_terminal_selector_selected"] is
                    (mode == 1) and
                document["parameters"][
                    "k65r65_b64_packed_terminal_mode"] == mode,
                "K65/R65/B64 diagnostic selector was not recorded exactly")
    require(k65r65_mode_zero["workload_digests"] ==
                k65r65_mode_one["workload_digests"],
            "K65/R65/B64 packed terminal changed encoded or recovered data")

    k65r65_attested = run(
        executable, True, attest_source=True,
        k65r65_b64_packed_terminal_mode=1,
        k=65, r=65, losses=1)
    k65r65_one_shot = run(
        executable, True, measure_one_shot_encode=True,
        k65r65_b64_packed_terminal_mode=1,
        k=65, r=65, losses=1)
    k65r65_attested_one_shot_zero = run(
        executable, True, attest_source=True,
        measure_one_shot_encode=True,
        k65r65_b64_packed_terminal_mode=0,
        k=65, r=65, losses=1)
    k65r65_attested_one_shot_one = run(
        executable, True, attest_source=True,
        measure_one_shot_encode=True,
        k65r65_b64_packed_terminal_mode=1,
        k=65, r=65, losses=1)
    k65r65_variants = (
        k65r65_mode_zero, k65r65_mode_one, k65r65_attested,
        k65r65_one_shot, k65r65_attested_one_shot_zero,
        k65r65_attested_one_shot_one,
    )
    for document in k65r65_variants:
        validate_common(document, True)
        validate_workload_digests(document)
    require(all(document["workload_digests"] ==
                k65r65_mode_zero["workload_digests"]
                for document in k65r65_variants),
            "K65/R65/B64 evidence options changed the workload")
    require(set(k65r65_attested["build"]) ==
                set(k65r65_mode_one["build"]) |
                {"source_commit", "source_tree", "source_tracked_dirty"} and
            set(k65r65_attested["parameters"]) ==
                set(k65r65_mode_one["parameters"]) | {"attest_source"},
            "K65/R65/B64 source attestation keys are ambiguous")
    require(set(k65r65_one_shot["build"]) ==
                set(k65r65_mode_one["build"]) | {
                    "k65r65_b64_packed_terminal_timed_one_shot_encode_api"} and
            set(k65r65_one_shot["parameters"]) ==
                set(k65r65_mode_one["parameters"]) |
                {"measure_one_shot_encode"} and
            set(k65r65_one_shot["metrics"]) ==
                set(k65r65_mode_one["metrics"]) | {"one_shot_encode"},
            "K65/R65/B64 one-shot encode keys are ambiguous")

    k65r65_neighbor_shapes = (
        (64, 65, 64), (66, 65, 64),
        (65, 64, 64), (65, 66, 64),
        (65, 65, 63), (65, 65, 65),
    )
    for neighbor_k, neighbor_r, neighbor_bytes in k65r65_neighbor_shapes:
        neighbor_documents = tuple(
            run(executable, True,
                k65r65_b64_packed_terminal_mode=mode,
                k=neighbor_k, r=neighbor_r, losses=1,
                shard_bytes=neighbor_bytes)
            for mode in (0, 1))
        for mode, document in enumerate(neighbor_documents):
            validate_common(document, True)
            validate_workload_digests(document)
            build = document["build"]
            require(build[
                        "k65r65_b64_packed_terminal_mode_latched"] == mode and
                    build[
                        "k65r65_b64_packed_terminal_selector_expected_selected"]
                        is False and
                    build[
                        "k65r65_b64_packed_terminal_selector_selected"] is
                        False,
                    "K65/R65/B64 neighbor selected the exact terminal")
        require(neighbor_documents[0]["workload_digests"] ==
                    neighbor_documents[1]["workload_digests"],
                "K65/R65/B64 neighbor mode changed the workload")

    malformed_k65r65_build_mode = copy.deepcopy(k65r65_mode_zero)
    malformed_k65r65_build_mode["build"][
        "k65r65_b64_packed_terminal_diagnostic_mode"] = False
    require_common_rejected(
        malformed_k65r65_build_mode,
        "Boolean K65/R65/B64 packed-terminal build mode")
    malformed_k65r65_disabled = copy.deepcopy(k65r65_mode_zero)
    malformed_k65r65_disabled["build"][
        "k65r65_b64_packed_terminal_diagnostic_disabled"] = False
    require_common_rejected(
        malformed_k65r65_disabled,
        "inconsistent K65/R65/B64 packed-terminal disabled marker")
    malformed_k65r65_parameter = copy.deepcopy(k65r65_mode_zero)
    malformed_k65r65_parameter["parameters"][
        "k65r65_b64_packed_terminal_mode"] = 1
    require_common_rejected(
        malformed_k65r65_parameter,
        "inconsistent K65/R65/B64 packed-terminal parameter")
    missing_k65r65_parameter = copy.deepcopy(k65r65_mode_zero)
    del missing_k65r65_parameter["parameters"][
        "k65r65_b64_packed_terminal_mode"]
    require_common_rejected(
        missing_k65r65_parameter,
        "missing K65/R65/B64 packed-terminal parameter")

    malformed_k65r65_latched = copy.deepcopy(k65r65_mode_zero)
    malformed_k65r65_latched["build"][
        "k65r65_b64_packed_terminal_mode_latched"] = True
    require_common_rejected(
        malformed_k65r65_latched,
        "Boolean K65/R65/B64 latched mode")
    malformed_k65r65_expected = copy.deepcopy(k65r65_mode_zero)
    malformed_k65r65_expected["build"][
        "k65r65_b64_packed_terminal_selector_expected_selected"] = True
    require_common_rejected(
        malformed_k65r65_expected,
        "incorrect K65/R65/B64 expected selector")
    malformed_k65r65_selected = copy.deepcopy(k65r65_mode_zero)
    malformed_k65r65_selected["build"][
        "k65r65_b64_packed_terminal_selector_selected"] = True
    require_common_rejected(
        malformed_k65r65_selected,
        "incorrect K65/R65/B64 actual selector")
    malformed_k65r65_contract = copy.deepcopy(k65r65_mode_zero)
    malformed_k65r65_contract["build"][
        "k65r65_b64_packed_terminal_selector_contract"] = "ambiguous"
    require_common_rejected(
        malformed_k65r65_contract,
        "ambiguous K65/R65/B64 selector contract")
    malformed_k65r65_ordinary_api = copy.deepcopy(k65r65_mode_zero)
    malformed_k65r65_ordinary_api["build"][
        "k65r65_b64_packed_terminal_timed_ordinary_encode_api"] = \
        "leo2_encode"
    require_common_rejected(
        malformed_k65r65_ordinary_api,
        "incorrect K65/R65/B64 ordinary timed API")

    unexpected_k65r65_attestation = copy.deepcopy(k65r65_mode_zero)
    unexpected_k65r65_attestation["build"]["source_commit"] = "unknown"
    require_common_rejected(
        unexpected_k65r65_attestation,
        "unrequested K65/R65/B64 source attestation")
    missing_k65r65_attestation = copy.deepcopy(k65r65_attested)
    del missing_k65r65_attestation["build"]["source_commit"]
    require_common_rejected(
        missing_k65r65_attestation,
        "incomplete K65/R65/B64 source attestation")
    malformed_k65r65_attestation_parameter = copy.deepcopy(k65r65_attested)
    malformed_k65r65_attestation_parameter["parameters"][
        "attest_source"] = False
    require_common_rejected(
        malformed_k65r65_attestation_parameter,
        "false K65/R65/B64 source-attestation parameter")

    unexpected_k65r65_one_shot = copy.deepcopy(k65r65_mode_zero)
    unexpected_k65r65_one_shot["build"][
        "k65r65_b64_packed_terminal_timed_one_shot_encode_api"] = \
        K65R65_B64_TIMED_ONE_SHOT_ENCODE_API
    require_common_rejected(
        unexpected_k65r65_one_shot,
        "unrequested K65/R65/B64 one-shot API marker")
    missing_k65r65_one_shot_marker = copy.deepcopy(k65r65_one_shot)
    del missing_k65r65_one_shot_marker["build"][
        "k65r65_b64_packed_terminal_timed_one_shot_encode_api"]
    require_common_rejected(
        missing_k65r65_one_shot_marker,
        "missing K65/R65/B64 one-shot API marker")
    malformed_k65r65_one_shot_api = copy.deepcopy(k65r65_one_shot)
    malformed_k65r65_one_shot_api["build"][
        "k65r65_b64_packed_terminal_timed_one_shot_encode_api"] = \
        "leo2_encode_batch"
    require_common_rejected(
        malformed_k65r65_one_shot_api,
        "incorrect K65/R65/B64 one-shot API marker")
    missing_k65r65_one_shot_metric = copy.deepcopy(k65r65_one_shot)
    del missing_k65r65_one_shot_metric["metrics"]["one_shot_encode"]
    require_common_rejected(
        missing_k65r65_one_shot_metric,
        "missing K65/R65/B64 one-shot metric")
    malformed_k65r65_one_shot_parameter = copy.deepcopy(k65r65_one_shot)
    malformed_k65r65_one_shot_parameter["parameters"][
        "measure_one_shot_encode"] = False
    require_common_rejected(
        malformed_k65r65_one_shot_parameter,
        "false K65/R65/B64 one-shot parameter")

    k65r65_parameter_mutations = (
        ("force_generic_decode", True),
        ("force_specialized_decode", True),
        ("force_tiled_decode", True),
        ("force_materialized_decode", True),
        ("warmup", True), ("warmup", -1),
        ("reuse", True), ("reuse", 0),
        ("iterations", True), ("iterations", 0),
        ("seed", True), ("seed", -1), ("seed", 1 << 64),
    )
    for name, value in k65r65_parameter_mutations:
        malformed = copy.deepcopy(k65r65_mode_zero)
        malformed["parameters"][name] = value
        require_common_rejected(
            malformed, f"invalid K65/R65/B64 {name} parameter")
    for name, value in (
            ("profile", "low_v1"), ("field", "gf16"),
            ("backend", "scalar"), ("thread_count", 2),
            ("thread_count", True)):
        malformed = copy.deepcopy(k65r65_mode_zero)
        malformed["resolved"][name] = value
        require_common_rejected(
            malformed, f"invalid K65/R65/B64 resolved {name}")

    require_schema_modes_rejected(
        executable, ("--k65r65-b64-packed-terminal-mode", "2"),
        "--k65r65-b64-packed-terminal-mode must be exactly 0 or 1")
    require_schema_modes_rejected(
        executable, ("--k65r65-b64-packed-terminal-mode", "0"),
        "--k65r65-b64-packed-terminal-mode requires explicit "
        "high/GF8/AVX2, batch=1, one thread, --skip-legacy, and "
        "--retain-samples")
    k65r65_contract_arguments = (
        "--k", "65", "--r", "65", "--profile", "high",
        "--field", "gf8", "--backend", "avx2", "--bytes", "64",
        "--loss", "1", "--batch", "1", "--reuse", "1",
        "--iterations", "1", "--warmup", "0", "--threads", "1",
        "--skip-legacy", "--retain-samples",
        "--k65r65-b64-packed-terminal-mode", "0",
    )
    k65r65_conflict_diagnostic = (
        "--k65r65-b64-packed-terminal-mode cannot be combined with another "
        "diagnostic mode, terminal-disable control, or decode-path override")
    k65r65_conflicting_modes = (
        ("--force-generic",),
        ("--force-specialized",),
        ("--force-tiled",),
        ("--force-materialized",),
        ("--report-decode-path",),
        ("--report-direct-executor",),
        ("--measure-one-shot-decode",),
        ("--one-shot-plan-setup-mode", "0"),
        ("--low-p32-b64-terminal-mode", "0"),
        ("--low-p128-b64-terminal-mode", "0"),
        ("--low-p16-partial-direct-output-mode", "0"),
        ("--gf8-avx2-walsh-locator-mode", "0"),
        ("--small-dual-regular-fallback-mode", "0"),
        ("--r1-small-reduction-mode", "0"),
        ("--r1-fixed-avx2-mode", "0"),
        ("--k8r3r4-t4-terminal-mode", "0"),
        ("--balanced-b64-terminal-mode", "0"),
        ("--k62r8-b64-fused-mode", "0"),
        ("--k66r16-b64-tail-mode", "0"),
        ("--high-t16-prepared-terminal-mode", "0"),
        ("--high-t8-two-block-b64-terminal-mode", "0"),
        ("--high-t8-two-block-b256-terminal-mode", "0"),
        ("--high-t8-two-block-b1024-terminal-mode", "0"),
        ("--disable-t8-full-parity-terminal",),
        ("--disable-k16r8-b256-terminal",),
        ("--disable-k9r5-b256-terminal",),
        ("--disable-k9r5-b1024-terminal",),
        ("--disable-k9r6r8-b256-terminal",),
    )
    for conflicting_mode in k65r65_conflicting_modes:
        require_schema_modes_rejected(
            executable, k65r65_contract_arguments + conflicting_mode,
            k65r65_conflict_diagnostic)

    k65r65_gfni_mode_zero = run(
        executable, True, measure_one_shot_encode=True,
        k65r65_b64_avx512_gfni_mode=0,
        k=65, r=65, losses=1)
    k65r65_gfni_mode_one = run(
        executable, True, measure_one_shot_encode=True,
        k65r65_b64_avx512_gfni_mode=1,
        k=65, r=65, losses=1)
    for mode, document in enumerate((
            k65r65_gfni_mode_zero, k65r65_gfni_mode_one)):
        require(document["schema"] == "leopard2-benchmark-v32",
                "K65/R65/B64 AVX-512/GFNI diagnostic schema changed")
        validate_common(document, True)
        validate_workload_digests(document)
        build = document["build"]
        selected = build[
            "k65r65_b64_avx512_gfni_selector_selected"]
        require(build["k65r65_b64_avx512_gfni_diagnostic_mode"] == mode and
                build["k65r65_b64_avx512_gfni_mode_latched"] == mode and
                build[
                    "k65r65_b64_avx512_gfni_observed_call_count"] ==
                    (2 if selected else 0) and
                document["parameters"][
                    "k65r65_b64_avx512_gfni_mode"] == mode,
                "K65/R65/B64 AVX-512/GFNI diagnostic selector was not "
                "recorded exactly")
    require(k65r65_gfni_mode_zero["workload_digests"] ==
                k65r65_gfni_mode_one["workload_digests"] and
            k65r65_gfni_mode_zero["resolved"] ==
                k65r65_gfni_mode_one["resolved"],
            "K65/R65/B64 AVX-512/GFNI mode changed the workload or "
            "resolved codec")

    k65r65_gfni_attested = run(
        executable, True, attest_source=True,
        measure_one_shot_encode=True,
        k65r65_b64_avx512_gfni_mode=1,
        k=65, r=65, losses=1)
    validate_common(k65r65_gfni_attested, True)
    validate_workload_digests(k65r65_gfni_attested)
    require(k65r65_gfni_attested["workload_digests"] ==
                k65r65_gfni_mode_one["workload_digests"] and
            set(k65r65_gfni_attested["build"]) ==
                set(k65r65_gfni_mode_one["build"]) |
                {"source_commit", "source_tree", "source_tracked_dirty"} and
            set(k65r65_gfni_attested["parameters"]) ==
                set(k65r65_gfni_mode_one["parameters"]) | {"attest_source"},
            "K65/R65/B64 AVX-512/GFNI source attestation is ambiguous")

    for neighbor_k, neighbor_r, neighbor_bytes in k65r65_neighbor_shapes:
        neighbor_documents = tuple(
            run(executable, True, measure_one_shot_encode=True,
                k65r65_b64_avx512_gfni_mode=mode,
                k=neighbor_k, r=neighbor_r, losses=1,
                shard_bytes=neighbor_bytes)
            for mode in (0, 1))
        for mode, document in enumerate(neighbor_documents):
            validate_common(document, True)
            validate_workload_digests(document)
            build = document["build"]
            require(build[
                        "k65r65_b64_avx512_gfni_mode_latched"] == mode and
                    build[
                        "k65r65_b64_avx512_gfni_selector_expected_selected"]
                        is False and
                    build[
                        "k65r65_b64_avx512_gfni_selector_selected"] is False and
                    build[
                        "k65r65_b64_avx512_gfni_observed_call_count"] == 0,
                    "K65/R65/B64 AVX-512/GFNI neighbor selected the exact "
                    "leaf")
        require(neighbor_documents[0]["workload_digests"] ==
                    neighbor_documents[1]["workload_digests"] and
                neighbor_documents[0]["resolved"] ==
                    neighbor_documents[1]["resolved"],
                "K65/R65/B64 AVX-512/GFNI neighbor mode changed the "
                "workload or resolved codec")

    for field_name, bad_value, label in (
            ("k65r65_b64_avx512_gfni_diagnostic_mode", False,
             "Boolean mode"),
            ("k65r65_b64_avx512_gfni_mode_latched", True,
             "Boolean latched mode"),
            ("k65r65_b64_avx512_gfni_kernel_qualified", 1,
             "numeric qualification"),
            ("k65r65_b64_avx512_gfni_observed_call_count", True,
             "Boolean call count"),
            ("k65r65_b64_avx512_gfni_selector_contract", "ambiguous",
             "ambiguous contract"),
            ("k65r65_b64_avx512_gfni_timed_ordinary_encode_api",
             "leo2_encode", "incorrect ordinary API"),
            ("k65r65_b64_avx512_gfni_timed_one_shot_encode_api",
             "leo2_encode_batch", "incorrect one-shot API")):
        malformed = copy.deepcopy(k65r65_gfni_mode_zero)
        malformed["build"][field_name] = bad_value
        require_common_rejected(
            malformed, f"K65/R65/B64 AVX-512/GFNI {label}")
    malformed_gfni_parameter = copy.deepcopy(k65r65_gfni_mode_zero)
    malformed_gfni_parameter["parameters"][
        "k65r65_b64_avx512_gfni_mode"] = 1
    require_common_rejected(
        malformed_gfni_parameter,
        "inconsistent K65/R65/B64 AVX-512/GFNI parameter")
    malformed_gfni_one_shot = copy.deepcopy(k65r65_gfni_mode_zero)
    malformed_gfni_one_shot["parameters"]["measure_one_shot_encode"] = False
    require_common_rejected(
        malformed_gfni_one_shot,
        "false K65/R65/B64 AVX-512/GFNI one-shot parameter")

    require_schema_modes_rejected(
        executable, ("--k65r65-b64-avx512-gfni-mode", "2"),
        "--k65r65-b64-avx512-gfni-mode must be exactly 0 or 1")
    require_schema_modes_rejected(
        executable, ("--k65r65-b64-avx512-gfni-mode", "0"),
        "--k65r65-b64-avx512-gfni-mode requires explicit high/GF8/AUTO, "
        "batch=1, one thread, --skip-legacy, --retain-samples, and "
        "--measure-one-shot-encode")
    k65r65_gfni_contract_arguments = (
        "--k", "65", "--r", "65", "--profile", "high",
        "--field", "gf8", "--backend", "auto", "--bytes", "64",
        "--loss", "1", "--batch", "1", "--reuse", "1",
        "--iterations", "1", "--warmup", "0", "--threads", "1",
        "--skip-legacy", "--retain-samples", "--measure-one-shot-encode",
        "--k65r65-b64-avx512-gfni-mode", "0",
    )
    require_schema_modes_rejected(
        executable, k65r65_gfni_contract_arguments + (
            "--k65r65-b64-packed-terminal-mode", "0"),
        "K65/R65/B64 diagnostic modes are mutually exclusive")

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

    disabled_k9_b1024_terminal = run(
        executable, True, disable_k9r5_b1024_terminal=True)
    validate_common(disabled_k9_b1024_terminal, True)
    validate_workload_digests(disabled_k9_b1024_terminal)
    require(disabled_k9_b1024_terminal["build"][
                "k9r5_b1024_terminal_diagnostic_disabled"] is True,
            "K9/R5/1024 attribution option was not recorded")
    require(disabled_k9_b1024_terminal["workload_digests"] ==
                external["workload_digests"],
            "inert K9/R5/1024 attribution option changed the workload")

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
    require_schema_modes_rejected(
        prevalidated_executable, k65r65_contract_arguments,
        "--k65r65-b64-packed-terminal-mode requires the ordinary benchmark")
    require_schema_modes_rejected(
        prevalidated_executable, k65r65_gfni_contract_arguments,
        "--k65r65-b64-avx512-gfni-mode requires the ordinary benchmark")
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

    transient_parameters = expected_external_parameters | {
        "measure_one_shot_decode", "one_shot_plan_setup_mode"}
    transient = run(
        executable, True, measure_one_shot_decode=True,
        one_shot_plan_setup_mode=3, k=17, r=5, losses=5)
    validate_transient_plan_report(
        transient, transient_parameters, setup_mode=3, walsh_mode=None)

    walsh_parameters = transient_parameters | {
        "gf8_avx2_walsh_locator_mode"}
    walsh_disabled = run(
        executable, True, measure_one_shot_decode=True,
        one_shot_plan_setup_mode=3, gf8_avx2_walsh_locator_mode=0,
        k=17, r=5, losses=5)
    walsh_enabled = run(
        executable, True, measure_one_shot_decode=True,
        one_shot_plan_setup_mode=3, gf8_avx2_walsh_locator_mode=1,
        k=17, r=5, losses=5)
    validate_transient_plan_report(
        walsh_disabled, walsh_parameters, setup_mode=3, walsh_mode=0)
    validate_transient_plan_report(
        walsh_enabled, walsh_parameters, setup_mode=3, walsh_mode=1)
    require(transient["workload_digests"] ==
                walsh_disabled["workload_digests"] ==
                walsh_enabled["workload_digests"],
            "transient setup or Walsh attribution changed recovered data")

    malformed_transient_mode = copy.deepcopy(transient)
    malformed_transient_mode["build"][
        "one_shot_plan_setup_diagnostic_mode"] = False
    require_common_rejected(
        malformed_transient_mode, "Boolean transient setup mode")
    malformed_transient_policy = copy.deepcopy(transient)
    malformed_transient_policy["build"][
        "one_shot_plan_setup_policy"] = "unknown"
    require_common_rejected(
        malformed_transient_policy, "unknown transient setup policy")
    malformed_transient_exclusions = copy.deepcopy(transient)
    malformed_transient_exclusions["build"][
        "one_shot_transient_execution_excludes"] = []
    require_common_rejected(
        malformed_transient_exclusions, "omitted transient wrapper exclusions")
    malformed_transient_guards = copy.deepcopy(transient)
    malformed_transient_guards["build"][
        "one_shot_transient_execution_guards"] = [
            "transform_route_and_exact_shard_bytes"]
    require_common_rejected(
        malformed_transient_guards, "omitted transient wrapper guard")
    malformed_setup_exclusions = copy.deepcopy(transient)
    malformed_setup_exclusions["resolved"][
        "one_shot_transient_setup_excludes"] = []
    require_common_rejected(
        malformed_setup_exclusions, "omitted transient setup exclusions")
    malformed_walsh_mode = copy.deepcopy(walsh_enabled)
    malformed_walsh_mode["build"][
        "gf8_avx2_walsh_locator_diagnostic_mode"] = True
    require_common_rejected(malformed_walsh_mode, "Boolean Walsh mode")
    inconsistent_walsh = copy.deepcopy(walsh_enabled)
    inconsistent_walsh["resolved"][
        "gf8_avx2_walsh_locator_enabled"] = False
    require_common_rejected(inconsistent_walsh, "inconsistent Walsh mode")
    malformed_digest_provenance = copy.deepcopy(transient)
    malformed_digest_provenance["workload_digests"][
        "recovered_originals_provenance"] = "ordinary_reused_plan"
    require_workload_digests_rejected(
        malformed_digest_provenance, "incorrect recovered digest provenance")

    r1_parameters = expected_external_parameters | {"measure_one_shot_decode"}
    r1_mode_zero = run(
        executable, True, measure_one_shot_decode=True,
        r1_small_reduction_mode=0, k=2, r=1, losses=1)
    r1_mode_one = run(
        executable, True, measure_one_shot_decode=True,
        r1_small_reduction_mode=1, k=2, r=1, losses=1)
    validate_r1_small_reduction_report(
        r1_mode_zero, 0, r1_parameters, 2, "k2_terminal")
    validate_r1_small_reduction_report(
        r1_mode_one, 1, r1_parameters, 2, "k2_terminal")
    require(r1_mode_zero["workload_digests"] ==
                r1_mode_one["workload_digests"],
            "R=1 reduction attribution changed the encoded or recovered data")

    r1_differential_zero = run(
        executable, True, measure_one_shot_decode=True,
        r1_small_reduction_mode=0, k=3, r=1, losses=1)
    r1_differential_one = run(
        executable, True, measure_one_shot_decode=True,
        r1_small_reduction_mode=1, k=3, r=1, losses=1)
    validate_r1_small_reduction_report(
        r1_differential_zero, 0, r1_parameters, 3, "pairwise")
    validate_r1_small_reduction_report(
        r1_differential_one, 1, r1_parameters, 3, "fused_final")
    require(r1_differential_zero["workload_digests"] ==
                r1_differential_one["workload_digests"],
            "R=1 differential attribution changed encoded or recovered data")

    malformed_r1_mode = copy.deepcopy(r1_mode_zero)
    malformed_r1_mode["build"][
        "r1_small_reduction_diagnostic_mode"] = False
    require_common_rejected(malformed_r1_mode, "Boolean R=1 mode")
    malformed_r1_enabled = copy.deepcopy(r1_mode_zero)
    malformed_r1_enabled["build"][
        "r1_small_reduction_codec_enabled"] = 0
    require_common_rejected(
        malformed_r1_enabled, "numeric R=1 codec-enabled marker")
    malformed_r1_path = copy.deepcopy(r1_mode_zero)
    malformed_r1_path["build"]["r1_encode_reduction_path"] = "unknown"
    require_common_rejected(malformed_r1_path, "unknown R=1 reduction path")
    malformed_r1_api = copy.deepcopy(r1_mode_zero)
    malformed_r1_api["build"]["r1_timed_encode_api"] = "leo2_encode"
    require_common_rejected(malformed_r1_api, "ambiguous R=1 timed API")

    require_schema_modes_rejected(
        executable, ("--r1-small-reduction-mode", "2"),
        "--r1-small-reduction-mode must be exactly 0 or 1")
    require_schema_modes_rejected(
        executable, ("--r1-small-reduction-mode", "0"),
        "--r1-small-reduction-mode requires explicit high/GF8/AVX2, "
        "R=1, one loss, batch=1, one thread, --skip-legacy, "
        "--retain-samples, and --measure-one-shot-decode")
    r1_contract_arguments = (
        "--k", "2", "--r", "1", "--profile", "high",
        "--field", "gf8", "--backend", "avx2", "--bytes", "64",
        "--loss", "1", "--batch", "1", "--reuse", "1",
        "--iterations", "1", "--warmup", "0", "--threads", "1",
        "--skip-legacy", "--retain-samples", "--measure-one-shot-decode",
        "--r1-small-reduction-mode", "0",
    )
    require_schema_modes_rejected(
        executable,
        r1_contract_arguments + ("--disable-k9r6r8-b256-terminal",),
        "--r1-small-reduction-mode requires explicit high/GF8/AVX2, "
        "R=1, one loss, batch=1, one thread, --skip-legacy, "
        "--retain-samples, and --measure-one-shot-decode")
    require_schema_modes_rejected(
        prevalidated_executable, r1_contract_arguments,
        "--r1-small-reduction-mode requires the ordinary benchmark")

    fixed_r1_parameters = r1_parameters | {"r1_fixed_avx2_mode"}
    fixed_r1_mode_zero = run(
        executable, True, measure_one_shot_decode=True,
        r1_fixed_avx2_mode=0, k=3, r=1, losses=1)
    fixed_r1_mode_one = run(
        executable, True, measure_one_shot_decode=True,
        r1_fixed_avx2_mode=1, k=3, r=1, losses=1)
    validate_r1_fixed_avx2_report(
        fixed_r1_mode_zero, 0, fixed_r1_parameters, "pairwise")
    validate_r1_fixed_avx2_report(
        fixed_r1_mode_one, 1, fixed_r1_parameters, "fixed_avx2")
    require(fixed_r1_mode_zero["workload_digests"] ==
                fixed_r1_mode_one["workload_digests"],
            "fixed AVX2 R=1 attribution changed encoded or recovered data")

    prevalidated_fixed_r1_parameters = fixed_r1_parameters | {
        "prevalidated_binding"}
    prevalidated_fixed_r1_mode_zero = run(
        prevalidated_executable, True, measure_one_shot_decode=True,
        r1_fixed_avx2_mode=0, prevalidated_binding=True,
        k=3, r=1, losses=1)
    prevalidated_fixed_r1_mode_one = run(
        prevalidated_executable, True, measure_one_shot_decode=True,
        r1_fixed_avx2_mode=1, prevalidated_binding=True,
        k=3, r=1, losses=1)
    validate_prevalidated_r1_fixed_avx2_report(
        prevalidated_fixed_r1_mode_zero, 0,
        prevalidated_fixed_r1_parameters, "pairwise")
    validate_prevalidated_r1_fixed_avx2_report(
        prevalidated_fixed_r1_mode_one, 1,
        prevalidated_fixed_r1_parameters, "fixed_avx2")
    require(fixed_r1_mode_zero["workload_digests"] ==
                fixed_r1_mode_one["workload_digests"] ==
                prevalidated_fixed_r1_mode_zero["workload_digests"] ==
                prevalidated_fixed_r1_mode_one["workload_digests"],
            "prevalidated binding changed encoded or recovered R=1 data")

    malformed_binding_attribution = copy.deepcopy(
        prevalidated_fixed_r1_mode_zero)
    malformed_binding_attribution["build"][
        "r1_prevalidated_binding_attribution"] = False
    require_prevalidated_rejected(
        malformed_binding_attribution, "disabled R=1 binding attribution")
    malformed_binding_api = copy.deepcopy(prevalidated_fixed_r1_mode_zero)
    malformed_binding_api["build"]["r1_timed_encode_api"] = (
        R1_TIMED_ENCODE_API)
    require_prevalidated_rejected(
        malformed_binding_api, "ambiguous R=1 binding timed API")
    malformed_binding_parameter = copy.deepcopy(
        prevalidated_fixed_r1_mode_zero)
    malformed_binding_parameter["parameters"]["prevalidated_binding"] = False
    require_prevalidated_rejected(
        malformed_binding_parameter, "disabled R=1 binding parameter")
    missing_r1_binding_setup = copy.deepcopy(prevalidated_fixed_r1_mode_zero)
    del missing_r1_binding_setup["metrics"]["decode_binding_setup"]
    require_prevalidated_rejected(
        missing_r1_binding_setup, "missing R=1 decode binding setup")

    malformed_fixed_r1_mode = copy.deepcopy(fixed_r1_mode_zero)
    malformed_fixed_r1_mode["build"][
        "r1_fixed_avx2_diagnostic_mode"] = False
    require_common_rejected(
        malformed_fixed_r1_mode, "Boolean fixed AVX2 R=1 mode")
    malformed_fixed_r1_candidate = copy.deepcopy(fixed_r1_mode_zero)
    malformed_fixed_r1_candidate["build"][
        "r1_fixed_avx2_candidate_enabled"] = False
    require_common_rejected(
        malformed_fixed_r1_candidate, "disabled fixed AVX2 candidate")
    malformed_fixed_r1_small = copy.deepcopy(fixed_r1_mode_zero)
    malformed_fixed_r1_small["build"][
        "r1_small_reduction_codec_enabled"] = True
    require_common_rejected(
        malformed_fixed_r1_small, "enabled small reduction in fixed report")
    malformed_fixed_r1_contract = copy.deepcopy(fixed_r1_mode_zero)
    malformed_fixed_r1_contract["build"][
        "r1_fixed_avx2_selector_contract"] = "ambiguous"
    require_common_rejected(
        malformed_fixed_r1_contract, "ambiguous fixed selector contract")
    malformed_fixed_r1_path = copy.deepcopy(fixed_r1_mode_zero)
    malformed_fixed_r1_path["build"]["r1_decode_reduction_path"] = "unknown"
    require_common_rejected(
        malformed_fixed_r1_path, "unknown fixed AVX2 R=1 path")

    require_schema_modes_rejected(
        executable, ("--r1-fixed-avx2-mode", "2"),
        "--r1-fixed-avx2-mode must be exactly 0 or 1")
    require_schema_modes_rejected(
        executable, ("--r1-fixed-avx2-mode", "0"),
        "--r1-fixed-avx2-mode requires explicit high/GF8/AVX2, "
        "R=1, one loss, batch=1, one thread, --skip-legacy, "
        "--retain-samples, and --measure-one-shot-decode")
    fixed_r1_contract_arguments = (
        "--k", "3", "--r", "1", "--profile", "high",
        "--field", "gf8", "--backend", "avx2", "--bytes", "64",
        "--loss", "1", "--batch", "1", "--reuse", "1",
        "--iterations", "1", "--warmup", "0", "--threads", "1",
        "--skip-legacy", "--retain-samples", "--measure-one-shot-decode",
        "--r1-fixed-avx2-mode", "0",
    )
    require_schema_modes_rejected(
        executable,
        fixed_r1_contract_arguments +
            ("--disable-k9r6r8-b256-terminal",),
        "--r1-fixed-avx2-mode requires explicit high/GF8/AVX2, "
        "R=1, one loss, batch=1, one thread, --skip-legacy, "
        "--retain-samples, and --measure-one-shot-decode")
    require_schema_modes_rejected(
        executable,
        fixed_r1_contract_arguments +
            ("--r1-small-reduction-mode", "0"),
        "--r1-small-reduction-mode requires explicit high/GF8/AVX2, "
        "R=1, one loss, batch=1, one thread, --skip-legacy, "
        "--retain-samples, and --measure-one-shot-decode")
    require_schema_modes_rejected(
        prevalidated_executable, fixed_r1_contract_arguments,
        "--r1-fixed-avx2-mode in the prevalidated benchmark requires "
        "--prevalidated-binding")
    require_schema_modes_rejected(
        prevalidated_executable, ("--prevalidated-binding",),
        "--prevalidated-binding requires --r1-fixed-avx2-mode in the "
        "prevalidated benchmark")

    path_diagnostic = (
        "--attest-source and --report-decode-path use distinct JSON schemas")
    for arguments in (
            ("--attest-source", "--report-decode-path"),
            ("--report-decode-path", "--attest-source")):
        require_schema_modes_rejected(
            executable, arguments, path_diagnostic)
    one_shot_diagnostic = (
        "--measure-one-shot-decode currently uses a standalone schema and "
        "cannot be combined with path attestation or unrelated "
        "source-attestation modes")
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
