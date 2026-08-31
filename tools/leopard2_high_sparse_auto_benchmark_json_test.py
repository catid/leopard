#!/usr/bin/env python3
"""Strict smoke contract for sparse-high production-AUTO telemetry."""

from __future__ import annotations

import copy
import json
import math
import os
import pathlib
import statistics
import subprocess
import sys
from typing import Any


SCHEMA = "leopard2-high-sparse-auto-benchmark-v1"
CONFIG_SCHEMA = "leopard2-benchmark-build-configuration/v16"
AUTHORITY_NOTE = (
    "raw telemetry is non-authoritative; authority requires the pinned paired "
    "runner"
)
AVX2_UNAVAILABLE_DIAGNOSTIC = (
    "leopard2 sparse-high AUTO benchmark: effective backend must be avx2 for "
    "sparse-high telemetry"
)
EXPECTED_METHODOLOGY = {
    "codec_setup_scope":
        "codec_create only; context and policy reused; destroy excluded",
    "binding_setup_scope":
        "binding_create only against preallocated descriptors; destroy excluded; "
        "null outside binding API",
    "execution_scope":
        "one public API call including ordinary validation; calls averaged within "
        "each sample",
    "route_witness_scope":
        "one untimed public execution, or two for binding reuse proof; disarmed "
        "before warmup",
    "timing_allocation_scope":
        "all benchmark-owned shards, descriptors, scratch, bindings, and sample "
        "vectors are prepared outside execution spans",
    "amortization_formula":
        "execution + codec_setup/reuse + binding_setup/reuse when applicable; one "
        "binding per codec",
    "affinity_scope":
        "affinity is established externally and captured before context creation",
    "production_autotuning": False,
}


class SmokeError(RuntimeError):
    pass


class SkipBenchmark(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SmokeError(message)


def exact_keys(value: object, expected: set[str], name: str) -> dict[str, Any]:
    require(isinstance(value, dict), f"{name} is not an object")
    actual = set(value)
    require(
        actual == expected,
        f"{name} keys differ: missing={sorted(expected - actual)}, "
        f"extra={sorted(actual - expected)}",
    )
    return value


def is_int(value: object) -> bool:
    return type(value) is int


def is_number(value: object) -> bool:
    return type(value) in (int, float) and math.isfinite(float(value))


def is_hex(value: object, digits: int, prefix: bool = False) -> bool:
    if not isinstance(value, str):
        return False
    expected = digits + (2 if prefix else 0)
    if len(value) != expected or (prefix and not value.startswith("0x")):
        return False
    body = value[2:] if prefix else value
    return all(character in "0123456789abcdef" for character in body)


def is_sparse_high_campaign_tuple(k: object, r: object, shard_bytes: object) -> bool:
    if not (is_int(k) and is_int(r) and is_int(shard_bytes)):
        return False
    qualified_shape = k in {2, 3, 4, 8, 12, 16} and r in {2, 4, 8, 16}
    if not qualified_shape:
        return False
    if shard_bytes == 4096:
        return True
    boundary_shape = (k, r) in {(2, 16), (16, 2)}
    return boundary_shape and shard_bytes in {1024, 1088, 2048, 4032, 4160, 65536}


def next_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def load_exact_json(raw: str) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise SmokeError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    try:
        value = json.loads(raw, object_pairs_hook=reject_duplicates)
    except json.JSONDecodeError as error:
        raise SmokeError(f"benchmark output is not JSON: {error}") from error
    require(isinstance(value, dict), "benchmark output is not a JSON object")
    return value


def validate_summary(
    value: object,
    sample_count: int,
    name: str,
    require_positive: bool,
) -> dict[str, Any]:
    summary = exact_keys(
        value,
        {"samples_us", "median_us", "mad_us", "minimum_us", "maximum_us"},
        name,
    )
    samples = summary["samples_us"]
    require(
        isinstance(samples, list) and len(samples) == sample_count,
        f"{name} sample count differs",
    )
    require(
        all(
            is_number(sample)
            and (float(sample) > 0.0 if require_positive else float(sample) >= 0.0)
            for sample in samples
        ),
        f"{name} contains invalid samples",
    )
    expected_median = statistics.median(float(sample) for sample in samples)
    expected_mad = statistics.median(
        abs(float(sample) - expected_median) for sample in samples
    )
    for key, expected in (
        ("median_us", expected_median),
        ("mad_us", expected_mad),
        ("minimum_us", min(float(sample) for sample in samples)),
        ("maximum_us", max(float(sample) for sample in samples)),
    ):
        require(is_number(summary[key]), f"{name}.{key} is not finite")
        require(
            math.isclose(float(summary[key]), expected, rel_tol=0.0, abs_tol=0.000002),
            f"{name}.{key} is not derived from samples",
        )
    return summary


def validate(
    document: dict[str, Any],
    expected: dict[str, Any],
    executable: pathlib.Path,
) -> None:
    exact_keys(
        document,
        {
            "schema",
            "authoritative",
            "authority_note",
            "build",
            "runtime",
            "parameters",
            "resolved",
            "qualification",
            "correctness",
            "memory",
            "metrics",
            "methodology",
        },
        "document",
    )
    require(document["schema"] == SCHEMA, "schema differs")
    require(document["authoritative"] is False, "raw telemetry claimed authority")
    require(document["authority_note"] == AUTHORITY_NOTE,
            "authority caveat differs")

    build = exact_keys(
        document["build"],
        {
            "compiler",
            "compiler_version",
            "cplusplus",
            "backend_variant",
            "build_type",
            "build_configuration_schema",
            "build_configuration_sha256",
            "source_commit",
            "source_tree",
            "source_tracked_dirty",
            "library_test_hooks",
            "high_sparse_tables_compiled",
            "high_sparse_auto_compiled_default",
        },
        "build",
    )
    for name in ("compiler", "compiler_version", "backend_variant", "build_type"):
        require(isinstance(build[name], str) and build[name], f"build.{name} missing")
    require(
        build["backend_variant"] in {"auto", "avx2"},
        "build backend variant is outside the AVX2 telemetry contract",
    )
    require(is_int(build["cplusplus"]), "build.cplusplus is not an integer")
    require(
        build["build_configuration_schema"] == CONFIG_SCHEMA,
        "build configuration generation differs",
    )
    require(
        is_hex(build["build_configuration_sha256"], 64),
        "build configuration digest is invalid",
    )
    for name in ("source_commit", "source_tree"):
        require(
            build[name] == "unknown" or is_hex(build[name], 40),
            f"build.{name} is invalid",
        )
    require(type(build["source_tracked_dirty"]) is bool, "dirty marker is not bool")
    require(build["library_test_hooks"] is False, "benchmark used test hooks")
    require(build["high_sparse_tables_compiled"] is True, "sparse tables absent")
    require(
        type(build["high_sparse_auto_compiled_default"]) is bool,
        "compiled AUTO marker is not bool",
    )

    runtime = exact_keys(
        document["runtime"],
        {"linux_procfs_affinity_attested", "executable_path", "allowed_cpus"},
        "runtime",
    )
    require(
        type(runtime["linux_procfs_affinity_attested"]) is bool,
        "runtime affinity marker is not bool",
    )
    require(isinstance(runtime["executable_path"], str), "runtime path is not text")
    require(isinstance(runtime["allowed_cpus"], list), "runtime CPUs are not a list")
    require(
        all(is_int(cpu) and cpu >= 0 for cpu in runtime["allowed_cpus"]),
        "runtime CPU identity is invalid",
    )
    if sys.platform.startswith("linux"):
        require(runtime["linux_procfs_affinity_attested"] is True,
                "Linux affinity is not attested")
        inherited_cpus = sorted(os.sched_getaffinity(0))
        require(
            runtime["allowed_cpus"] == inherited_cpus,
            "runtime affinity differs from the parser process affinity",
        )
        require(
            pathlib.Path(runtime["executable_path"]).samefile(executable),
            "runtime executable path differs",
        )
    else:
        require(
            runtime["linux_procfs_affinity_attested"] is False
            and runtime["executable_path"] == ""
            and runtime["allowed_cpus"] == [],
            "non-Linux runtime identity must be empty and unattested",
        )

    parameters = exact_keys(
        document["parameters"],
        {
            "requested_backend",
            "policy",
            "api",
            "K",
            "R",
            "Q",
            "parity_index",
            "profile",
            "field",
            "shard_layout",
            "codec_flags",
            "shard_bytes",
            "batch",
            "reuse",
            "iterations",
            "setup_iterations",
            "calls_per_sample",
            "warmups",
            "requested_thread_count",
            "seed",
            "memory_mib",
        },
        "parameters",
    )
    for name, value in expected["parameters"].items():
        require(parameters[name] == value, f"parameters.{name} differs")
    require(parameters["Q"] == 1, "benchmark is not Q=1")
    require(parameters["profile"] == "legacy_high_v1", "profile differs")
    require(parameters["field"] == "gf8", "field differs")
    require(parameters["shard_layout"] == "native_v1", "layout differs")
    require(parameters["codec_flags"] == 0, "codec flags are nonzero")
    for name in (
        "K", "R", "Q", "parity_index", "codec_flags", "shard_bytes", "batch", "reuse",
        "iterations", "setup_iterations", "calls_per_sample", "warmups",
        "requested_thread_count", "seed", "memory_mib",
    ):
        require(is_int(parameters[name]), f"parameters.{name} is not an integer")
    require(0 <= parameters["parity_index"] < parameters["R"], "parity row invalid")
    require(
        is_sparse_high_campaign_tuple(
            parameters["K"], parameters["R"], parameters["shard_bytes"]),
        "cell is outside the 36-cell sparse-high campaign envelope",
    )

    resolved = exact_keys(
        document["resolved"],
        {
            "effective_backend",
            "thread_count",
            "parent_count",
            "padded_side",
            "direct_generator_rows",
            "auto_direct_selected",
            "selected_route",
        },
        "resolved",
    )
    require(resolved["effective_backend"] == "avx2", "effective backend differs")
    require(
        is_int(resolved["thread_count"])
        and resolved["thread_count"] == parameters["requested_thread_count"],
        "effective thread count differs",
    )
    expected_padded_side = next_power_of_two(parameters["R"])
    expected_parent_count = next_power_of_two(
        parameters["K"] + expected_padded_side
    )
    require(
        is_int(resolved["parent_count"])
        and resolved["parent_count"] == expected_parent_count,
        "parent count differs from the legacy-high layout",
    )
    require(
        is_int(resolved["padded_side"])
        and resolved["padded_side"] == expected_padded_side,
        "padded side differs from the legacy-high layout",
    )
    require(
        is_int(resolved["direct_generator_rows"])
        and resolved["direct_generator_rows"] == expected["rows"],
        "prepared row count differs",
    )
    require(
        resolved["auto_direct_selected"] is expected["direct"],
        "AUTO selector result differs",
    )
    require(
        resolved["selected_route"] == ("direct" if expected["direct"] else "transform"),
        "selected route differs",
    )

    qualification = exact_keys(
        document["qualification"],
        {
            "route_witness_armed",
            "witness_public_executions",
            "expected_item_calls",
            "direct_calls",
            "transform_calls",
            "witness_disabled_before_timing",
        },
        "qualification",
    )
    public_executions = 2 if parameters["api"] == "binding" else 1
    expected_calls = public_executions * parameters["batch"]
    require(qualification["route_witness_armed"] is True, "witness was not armed")
    require(
        qualification["witness_disabled_before_timing"] is True,
        "witness remained enabled during timing",
    )
    require(
        qualification["witness_public_executions"] == public_executions
        and qualification["expected_item_calls"] == expected_calls,
        "witness call accounting differs",
    )
    require(
        qualification["direct_calls"] == (expected_calls if expected["direct"] else 0)
        and qualification["transform_calls"]
        == (0 if expected["direct"] else expected_calls),
        "actual route witness differs",
    )
    for name in (
        "witness_public_executions", "expected_item_calls", "direct_calls",
        "transform_calls",
    ):
        require(is_int(qualification[name]), f"qualification.{name} is not integer")

    correctness = exact_keys(
        document["correctness"],
        {
            "oracle_algorithm",
            "independent_oracle_match",
            "input_immutable",
            "unrequested_outputs_untouched",
            "post_timing_recheck_match",
            "input_checksum_fnv1a64",
            "parity_checksum_fnv1a64",
        },
        "correctness",
    )
    require(
        correctness["oracle_algorithm"]
        == "legacy-gf8-direct-systematic-generator-v1",
        "oracle identity differs",
    )
    for name in (
        "independent_oracle_match",
        "input_immutable",
        "unrequested_outputs_untouched",
        "post_timing_recheck_match",
    ):
        require(correctness[name] is True, f"correctness.{name} failed")
    for name in ("input_checksum_fnv1a64", "parity_checksum_fnv1a64"):
        require(is_hex(correctness[name], 16, prefix=True), f"{name} is invalid")

    memory = exact_keys(
        document["memory"],
        {
            "scratch_alignment",
            "scratch_bytes_per_item",
            "scratch_bytes_batch",
            "logical_input_bytes_per_call",
            "requested_output_bytes_per_call",
            "estimated_benchmark_storage_bytes",
        },
        "memory",
    )
    require(
        all(is_int(value) and value >= 0 for value in memory.values()),
        "memory accounting is not nonnegative integers",
    )
    require(memory["scratch_alignment"] == 64, "scratch alignment differs")
    require(
        memory["scratch_bytes_batch"]
        == memory["scratch_bytes_per_item"] * parameters["batch"],
        "batch scratch accounting differs",
    )
    require(
        memory["logical_input_bytes_per_call"]
        == parameters["K"] * parameters["shard_bytes"] * parameters["batch"],
        "logical input byte accounting differs",
    )
    require(
        memory["requested_output_bytes_per_call"]
        == parameters["shard_bytes"] * parameters["batch"],
        "requested output byte accounting differs",
    )
    expected_storage = (
        (parameters["K"] + parameters["R"] + 1)
        * parameters["shard_bytes"]
        * parameters["batch"]
        + memory["scratch_bytes_batch"]
    )
    require(
        memory["estimated_benchmark_storage_bytes"] == expected_storage,
        "estimated benchmark storage accounting differs",
    )

    metrics = exact_keys(
        document["metrics"],
        {"codec_setup", "binding_setup", "execution", "amortized"},
        "metrics",
    )
    codec_setup = validate_summary(
        metrics["codec_setup"], parameters["setup_iterations"],
        "metrics.codec_setup", False,
    )
    if parameters["api"] == "binding":
        binding_setup = validate_summary(
            metrics["binding_setup"], parameters["setup_iterations"],
            "metrics.binding_setup", False,
        )
    else:
        require(metrics["binding_setup"] is None, "binding setup must be null")
        binding_setup = None
    execution = validate_summary(
        metrics["execution"], parameters["iterations"],
        "metrics.execution", True,
    )
    amortized = exact_keys(
        metrics["amortized"],
        {
            "reuse_count",
            "derived_median_us_per_api_call",
            "logical_input_GB_per_s",
            "requested_parity_output_GB_per_s",
        },
        "metrics.amortized",
    )
    require(
        is_int(amortized["reuse_count"])
        and amortized["reuse_count"] == parameters["reuse"],
        "reuse differs",
    )
    expected_us = float(execution["median_us"]) + float(
        codec_setup["median_us"]
    ) / parameters["reuse"]
    if binding_setup is not None:
        expected_us += float(binding_setup["median_us"]) / parameters["reuse"]
    require(
        is_number(amortized["derived_median_us_per_api_call"])
        and math.isclose(
            float(amortized["derived_median_us_per_api_call"]),
            expected_us,
            rel_tol=0.0,
            abs_tol=0.000003,
        ),
        "amortized time is not derived from setup and execution",
    )
    for key, byte_key in (
        ("logical_input_GB_per_s", "logical_input_bytes_per_call"),
        ("requested_parity_output_GB_per_s", "requested_output_bytes_per_call"),
    ):
        expected_rate = memory[byte_key] / (expected_us * 1000.0)
        require(
            is_number(amortized[key])
            and math.isclose(
                float(amortized[key]), expected_rate,
                rel_tol=0.000005,
                abs_tol=0.000003,
            ),
            f"{key} is not byte/time-derived",
        )

    methodology = exact_keys(
        document["methodology"],
        {
            "codec_setup_scope",
            "binding_setup_scope",
            "execution_scope",
            "route_witness_scope",
            "timing_allocation_scope",
            "amortization_formula",
            "affinity_scope",
            "production_autotuning",
        },
        "methodology",
    )
    require(methodology == EXPECTED_METHODOLOGY, "methodology contract differs")


def run_process(executable: pathlib.Path, arguments: list[str]) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    threads = "1"
    if "--threads" in arguments:
        threads = arguments[arguments.index("--threads") + 1]
    environment.update(
        OMP_DYNAMIC="FALSE",
        OMP_NUM_THREADS=threads,
        OMP_THREAD_LIMIT=threads,
    )
    return subprocess.run(
        [str(executable), *arguments],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=60,
        env=environment,
    )


def base_arguments() -> list[str]:
    return [
        "--k", "2", "--r", "16", "--bytes", "1024",
        "--parity-index", "7", "--iterations", "3",
        "--setup-iterations", "3", "--calls-per-sample", "1",
        "--warmups", "1", "--reuse", "8", "--memory-mib", "16",
    ]


def synthetic_document(
    executable: pathlib.Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    parameters = {
        "requested_backend": "auto",
        "policy": "tables_off_auto_off",
        "api": "one_shot",
        "K": 2,
        "R": 16,
        "Q": 1,
        "parity_index": 7,
        "profile": "legacy_high_v1",
        "field": "gf8",
        "shard_layout": "native_v1",
        "codec_flags": 0,
        "shard_bytes": 1024,
        "batch": 1,
        "reuse": 1,
        "iterations": 1,
        "setup_iterations": 1,
        "calls_per_sample": 1,
        "warmups": 0,
        "requested_thread_count": 1,
        "seed": 0x48534155544F5631,
        "memory_mib": 16,
    }
    linux_runtime = sys.platform.startswith("linux")
    if linux_runtime:
        allowed_cpus = sorted(os.sched_getaffinity(0))
    else:
        allowed_cpus = []
    summary_one = {
        "samples_us": [1.0],
        "median_us": 1.0,
        "mad_us": 0.0,
        "minimum_us": 1.0,
        "maximum_us": 1.0,
    }
    summary_three = {
        "samples_us": [3.0],
        "median_us": 3.0,
        "mad_us": 0.0,
        "minimum_us": 3.0,
        "maximum_us": 3.0,
    }
    document = {
        "schema": SCHEMA,
        "authoritative": False,
        "authority_note": AUTHORITY_NOTE,
        "build": {
            "compiler": "synthetic",
            "compiler_version": "synthetic-v1",
            "cplusplus": 201103,
            "backend_variant": "auto",
            "build_type": "Release",
            "build_configuration_schema": CONFIG_SCHEMA,
            "build_configuration_sha256": "0" * 64,
            "source_commit": "unknown",
            "source_tree": "unknown",
            "source_tracked_dirty": False,
            "library_test_hooks": False,
            "high_sparse_tables_compiled": True,
            "high_sparse_auto_compiled_default": False,
        },
        "runtime": {
            "linux_procfs_affinity_attested": linux_runtime,
            "executable_path": str(executable) if linux_runtime else "",
            "allowed_cpus": allowed_cpus,
        },
        "parameters": parameters,
        "resolved": {
            "effective_backend": "avx2",
            "thread_count": 1,
            "parent_count": 32,
            "padded_side": 16,
            "direct_generator_rows": 0,
            "auto_direct_selected": False,
            "selected_route": "transform",
        },
        "qualification": {
            "route_witness_armed": True,
            "witness_public_executions": 1,
            "expected_item_calls": 1,
            "direct_calls": 0,
            "transform_calls": 1,
            "witness_disabled_before_timing": True,
        },
        "correctness": {
            "oracle_algorithm": "legacy-gf8-direct-systematic-generator-v1",
            "independent_oracle_match": True,
            "input_immutable": True,
            "unrequested_outputs_untouched": True,
            "post_timing_recheck_match": True,
            "input_checksum_fnv1a64": "0x0123456789abcdef",
            "parity_checksum_fnv1a64": "0xfedcba9876543210",
        },
        "memory": {
            "scratch_alignment": 64,
            "scratch_bytes_per_item": 64,
            "scratch_bytes_batch": 64,
            "logical_input_bytes_per_call": 2048,
            "requested_output_bytes_per_call": 1024,
            "estimated_benchmark_storage_bytes": 19520,
        },
        "metrics": {
            "codec_setup": summary_one,
            "binding_setup": None,
            "execution": summary_three,
            "amortized": {
                "reuse_count": 1,
                "derived_median_us_per_api_call": 4.0,
                "logical_input_GB_per_s": 0.512,
                "requested_parity_output_GB_per_s": 0.256,
            },
        },
        "methodology": dict(EXPECTED_METHODOLOGY),
    }
    expected = {
        "parameters": parameters,
        "rows": 0,
        "direct": False,
    }
    return document, expected


def exercise_mutation_rejections(
    document: dict[str, Any],
    expected: dict[str, Any],
    executable: pathlib.Path,
) -> None:
    mutations: list[tuple[str, dict[str, Any]]] = []
    value = copy.deepcopy(document)
    value["extra"] = None
    mutations.append(("extra top-level key", value))
    value = copy.deepcopy(document)
    del value["build"]["library_test_hooks"]
    mutations.append(("missing build key", value))
    value = copy.deepcopy(document)
    value["build"]["library_test_hooks"] = True
    mutations.append(("hook archive claim", value))
    value = copy.deepcopy(document)
    value["authority_note"] = f"prefix {AUTHORITY_NOTE} suffix"
    mutations.append(("authority caveat drift", value))
    value = copy.deepcopy(document)
    value["runtime"]["allowed_cpus"] = [999999, 999999]
    mutations.append(("runtime CPU provenance", value))
    value = copy.deepcopy(document)
    value["resolved"]["thread_count"] = True
    mutations.append(("boolean thread count", value))
    value = copy.deepcopy(document)
    value["resolved"]["direct_generator_rows"] = False
    mutations.append(("boolean prepared row count", value))
    value = copy.deepcopy(document)
    value["resolved"]["parent_count"] = 1
    mutations.append(("parent layout count", value))
    value = copy.deepcopy(document)
    value["resolved"]["padded_side"] = 1
    mutations.append(("padded layout side", value))
    value = copy.deepcopy(document)
    value["parameters"]["codec_flags"] = False
    mutations.append(("boolean codec flags", value))
    value = copy.deepcopy(document)
    value["qualification"]["transform_calls"] += 1
    mutations.append(("route witness count", value))
    value = copy.deepcopy(document)
    value["metrics"]["execution"]["extra"] = 0
    mutations.append(("extra timing key", value))
    value = copy.deepcopy(document)
    value["metrics"]["amortized"]["derived_median_us_per_api_call"] += 1
    mutations.append(("amortized formula", value))
    value = copy.deepcopy(document)
    value["metrics"]["amortized"]["reuse_count"] = True
    mutations.append(("boolean reuse count", value))
    value = copy.deepcopy(document)
    value["build"]["backend_variant"] = "scalar"
    mutations.append(("unsupported build backend", value))
    value = copy.deepcopy(document)
    value["memory"]["estimated_benchmark_storage_bytes"] = 0
    mutations.append(("storage accounting", value))
    value = copy.deepcopy(document)
    value["memory"]["scratch_alignment"] = 3
    mutations.append(("scratch alignment", value))
    value = copy.deepcopy(document)
    value["methodology"]["execution_scope"] = "arbitrary"
    mutations.append(("methodology drift", value))
    value = copy.deepcopy(document)
    value["parameters"]["Q"] = True
    mutations.append(("boolean integer", value))
    for name, value in mutations:
        expect_validation_failure(value, expected, executable, name)


def run_success(
    executable: pathlib.Path,
    extra: list[str],
    expected: dict[str, Any],
) -> dict[str, Any]:
    completed = run_process(executable, [*base_arguments(), *extra])
    if completed.returncode != 0:
        if (completed.stderr.strip() == AVX2_UNAVAILABLE_DIAGNOSTIC
                and not completed.stdout.strip()):
            raise SkipBenchmark("host AUTO/AVX2 runtime is unavailable")
        raise SmokeError(
            f"benchmark failed ({completed.returncode}): {completed.stderr.strip()}"
        )
    document = load_exact_json(completed.stdout)
    validate(document, expected, executable)
    return document


def run_failure(executable: pathlib.Path, arguments: list[str], name: str) -> None:
    completed = run_process(executable, arguments)
    require(completed.returncode != 0, f"benchmark accepted {name}")
    require(not completed.stdout.strip(), f"failed {name} emitted JSON")


def expect_validation_failure(
    document: dict[str, Any],
    expected: dict[str, Any],
    executable: pathlib.Path,
    name: str,
) -> None:
    try:
        validate(document, expected, executable)
    except SmokeError:
        return
    raise SmokeError(f"strict parser accepted mutation: {name}")


def main() -> int:
    if len(sys.argv) != 2:
        raise SmokeError("usage: leopard2_high_sparse_auto_benchmark_json_test.py EXECUTABLE")
    executable = pathlib.Path(sys.argv[1]).resolve(strict=True)
    fixture, fixture_expected = synthetic_document(executable)
    validate(fixture, fixture_expected, executable)
    exercise_mutation_rejections(fixture, fixture_expected, executable)

    common_parameters = {
        "K": 2,
        "R": 16,
        "Q": 1,
        "parity_index": 7,
        "shard_bytes": 1024,
        "reuse": 8,
        "iterations": 3,
        "setup_iterations": 3,
        "calls_per_sample": 1,
        "warmups": 1,
        "seed": 0x48534155544F5631,
        "memory_mib": 16,
    }
    valid = [
        *base_arguments(), "--api", "batch", "--batch", "1",
        "--backend", "auto", "--threads", "1",
        "--policy", "tables-on-auto-on",
    ]
    invalid_cases = [
        ([*valid, "--force-direct"], "force-direct"),
        ([*valid, "--force-transform"], "force-transform"),
        ([*valid, "--mode", "direct"], "mode override"),
        ([*base_arguments(), "--api", "one-shot", "--batch", "4",
          "--backend", "auto", "--threads", "1", "--policy",
          "tables-on-auto-on"], "one-shot batch four"),
        ([*base_arguments(), "--api", "batch", "--batch", "2",
          "--backend", "auto", "--threads", "1", "--policy",
          "tables-on-auto-on"], "unsupported batch"),
        ([*valid, "--threads", "2"], "unsupported thread count"),
        ([*valid, "--backend", "scalar"], "unsupported backend"),
        ([*valid, "--policy", "tables-off-auto-on"], "invalid policy state"),
        ([*valid, "--reuse", "2"], "unsupported reuse"),
        ([*valid, "--parity-index", "16"], "invalid parity row"),
        ([*valid, "--iterations", "0"], "zero iterations"),
        ([*valid, "--k", "5", "--bytes", "4096"],
         "unsupported K at the standard byte count"),
        ([*valid, "--r", "3", "--bytes", "4096"],
         "unsupported R at the standard byte count"),
        ([*valid, "--k", "3", "--r", "16", "--bytes", "1024"],
         "non-boundary shape at a boundary byte count"),
        ([*valid, "--bytes", "1056"], "unsupported boundary byte count"),
        ([*base_arguments(), "--backend", "auto", "--policy",
          "tables-on-auto-on"], "missing API"),
        ([*base_arguments(), "--api", "batch", "--policy",
          "tables-on-auto-on"], "missing backend"),
        ([*base_arguments(), "--api", "batch", "--backend", "auto"],
         "missing policy"),
    ]
    for arguments, name in invalid_cases:
        run_failure(executable, arguments, name)

    cases = [
        (
            ["--api", "one-shot", "--batch", "1", "--backend", "auto",
             "--threads", "1", "--policy", "tables-off-auto-off"],
            {"requested_backend": "auto", "policy": "tables_off_auto_off",
             "api": "one_shot", "batch": 1, "requested_thread_count": 1},
            0,
            False,
        ),
        (
            ["--api", "batch", "--batch", "4", "--backend", "auto",
             "--threads", "1", "--policy", "tables-on-auto-off"],
            {"requested_backend": "auto", "policy": "tables_on_auto_off",
             "api": "batch", "batch": 4, "requested_thread_count": 1},
            16,
            False,
        ),
        (
            ["--api", "one-shot", "--batch", "1", "--backend", "auto",
             "--threads", "1", "--policy", "tables-on-auto-on"],
            {"requested_backend": "auto", "policy": "tables_on_auto_on",
             "api": "one_shot", "batch": 1, "requested_thread_count": 1},
            16,
            True,
        ),
        (
            ["--api", "batch", "--batch", "1", "--backend", "auto",
             "--threads", "1", "--policy", "tables-on-auto-on"],
            {"requested_backend": "auto", "policy": "tables_on_auto_on",
             "api": "batch", "batch": 1, "requested_thread_count": 1},
            16,
            True,
        ),
        (
            ["--api", "batch", "--batch", "4", "--backend", "auto",
             "--threads", "1", "--policy", "tables-on-auto-on"],
            {"requested_backend": "auto", "policy": "tables_on_auto_on",
             "api": "batch", "batch": 4, "requested_thread_count": 1},
            16,
            True,
        ),
        (
            ["--api", "batch", "--batch", "16", "--backend", "auto",
             "--threads", "1", "--policy", "tables-on-auto-on"],
            {"requested_backend": "auto", "policy": "tables_on_auto_on",
             "api": "batch", "batch": 16, "requested_thread_count": 1},
            16,
            True,
        ),
        (
            ["--api", "binding", "--batch", "1", "--backend", "auto",
             "--threads", "1", "--policy", "tables-on-auto-on"],
            {"requested_backend": "auto", "policy": "tables_on_auto_on",
             "api": "binding", "batch": 1, "requested_thread_count": 1},
            16,
            True,
        ),
        (
            ["--api", "binding", "--batch", "4", "--backend", "auto",
             "--threads", "1", "--policy", "tables-on-auto-on"],
            {"requested_backend": "auto", "policy": "tables_on_auto_on",
             "api": "binding", "batch": 4, "requested_thread_count": 1},
            16,
            True,
        ),
        (
            ["--api", "binding", "--batch", "16", "--backend", "auto",
             "--threads", "1", "--policy", "tables-on-auto-on"],
            {"requested_backend": "auto", "policy": "tables_on_auto_on",
             "api": "binding", "batch": 16, "requested_thread_count": 1},
            16,
            True,
        ),
        (
            ["--api", "batch", "--batch", "1", "--backend", "avx2",
             "--threads", "1", "--policy", "tables-on-auto-on"],
            {"requested_backend": "avx2", "policy": "tables_on_auto_on",
             "api": "batch", "batch": 1, "requested_thread_count": 1},
            16,
            False,
        ),
        (
            ["--api", "batch", "--batch", "4", "--backend", "auto",
             "--threads", "4", "--policy", "tables-on-auto-on"],
            {"requested_backend": "auto", "policy": "tables_on_auto_on",
             "api": "batch", "batch": 4, "requested_thread_count": 4},
            16,
            False,
        ),
    ]
    first_document: dict[str, Any] | None = None
    first_expected: dict[str, Any] | None = None
    for arguments, identity, rows, direct in cases:
        expected = {
            "parameters": {**common_parameters, **identity},
            "rows": rows,
            "direct": direct,
        }
        document = run_success(executable, arguments, expected)
        if first_document is None:
            first_document = document
            first_expected = expected

    require(first_document is not None and first_expected is not None,
            "no positive document was retained")
    exercise_mutation_rejections(first_document, first_expected, executable)

    print("sparse-high production-AUTO benchmark schema: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SkipBenchmark as error:
        print(f"sparse-high production-AUTO benchmark schema: SKIP ({error})")
        raise SystemExit(0)
    except (SmokeError, OSError, subprocess.SubprocessError) as error:
        print(f"sparse-high production-AUTO benchmark schema: FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
