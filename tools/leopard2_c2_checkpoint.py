#!/usr/bin/env python3
"""Validate and merge standalone C2 truncated-transform checkpoint results."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pathlib
import platform
import re
import subprocess
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-c2-checkpoint/v2"
RESULT_SCHEMA = "leopard2-c2-cpp-v2"
EXPECTED_BACKENDS = {
    "auto": "avx2",
    "scalar": "scalar",
    "ssse3": "ssse3",
    "avx2": "avx2",
}
EXPECTED_CASE_MATRIX = (
    ("gf8_prefix_n128_s0_forward", 8, 128, 0, "forward", 65, 33, 2112),
    ("gf8_prefix_n128_s0_inverse", 8, 128, 0, "inverse", 65, 33, 2112),
    ("gf8_last_coset_n128_s128_forward", 8, 128, 128, "forward", 65, 33, 2112),
    ("gf8_last_coset_n128_s128_inverse", 8, 128, 128, "inverse", 65, 33, 2112),
    ("gf8_benchmark_n256_s0_forward", 8, 256, 0, "forward", 129, 65, 4160),
    ("gf8_benchmark_n256_s0_inverse", 8, 256, 0, "inverse", 129, 65, 4160),
    ("gf16_small_n256_s0_forward", 16, 256, 0, "forward", 129, 65, 4160),
    ("gf16_small_n256_s0_inverse", 16, 256, 0, "inverse", 129, 65, 4160),
    ("gf16_small_last_n256_s65280_forward", 16, 256, 65280, "forward", 129, 65, 4160),
    ("gf16_small_last_n256_s65280_inverse", 16, 256, 65280, "inverse", 129, 65, 4160),
    ("gf16_medium_n1024_s0_forward", 16, 1024, 0, "forward", 513, 257, 16448),
    ("gf16_medium_n1024_s0_inverse", 16, 1024, 0, "inverse", 513, 257, 16448),
    ("gf16_medium_last_n1024_s64512_forward", 16, 1024, 64512, "forward", 513, 257, 16448),
    ("gf16_medium_last_n1024_s64512_inverse", 16, 1024, 64512, "inverse", 513, 257, 16448),
    ("gf16_large_n4096_s0_forward", 16, 4096, 0, "forward", 2049, 1025, 65600),
    ("gf16_large_n4096_s0_inverse", 16, 4096, 0, "inverse", 2049, 1025, 65600),
    ("gf16_large_last_n4096_s61440_forward", 16, 4096, 61440, "forward", 2049, 1025, 65600),
    ("gf16_large_last_n4096_s61440_inverse", 16, 4096, 61440, "inverse", 2049, 1025, 65600),
    ("gf16_deep_n8192_s57344_forward", 16, 8192, 57344, "forward", 4097, 2049, 131136),
    ("gf16_deep_n8192_s57344_inverse", 16, 8192, 57344, "inverse", 4097, 2049, 131136),
    ("gf16_irregular_forward", 16, 1024, 32768, "forward", 466, 473, 30272),
    ("gf16_irregular_inverse", 16, 1024, 32768, "inverse", 466, 473, 30272),
)
EXPECTED_BENCHMARK_MATRIX = (
    ("gf8_benchmark_n256_s0_forward", 8, 256, "forward"),
    ("gf8_benchmark_n256_s0_inverse", 8, 256, "inverse"),
    ("gf16_medium_n1024_s0_forward", 16, 1024, "forward"),
    ("gf16_medium_n1024_s0_inverse", 16, 1024, "inverse"),
    ("gf16_large_n4096_s0_forward", 16, 4096, "forward"),
    ("gf16_large_n4096_s0_inverse", 16, 4096, "inverse"),
)
CASE_PARITY_KEYS = (
    "name",
    "field_bits",
    "parent_size",
    "shift",
    "direction",
    "active",
    "requested",
    "compared_bytes",
    "pair_operations",
    "complete_operations",
    "maximum_complete_block",
    "scratch_slots",
    "peak_live_slots",
    "execution_scratch_bytes",
    "padded_scratch_bytes",
    "serialized_schedule_bytes",
    "resident_plan_bytes_excluding_allocator_headers",
)
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
GIT_SHA_PATTERN = re.compile(r"^[0-9a-f]{40}$")
DIGEST_PATTERN = re.compile(r"^0x[0-9a-f]{16}$")
EXPECTED_CORRECTNESS_DIGEST = "0xeebedd5febef7e07"
MAX_RUNTIME_RATIO = 100.0
MIN_RUNTIME_RATIO = 0.01
EXPECTED_SOURCE_RELATIVE = pathlib.Path(
    "experiments/leopard2/c2_truncated_cpp/c2_truncated.cpp"
)


class CheckpointError(RuntimeError):
    pass


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise CheckpointError(f"cannot hash {path}: {error}") from error
    return digest.hexdigest()


def repository_head(repository: pathlib.Path) -> str:
    try:
        value = subprocess.check_output(
            ["git", "-C", str(repository), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.PIPE,
        ).strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise CheckpointError(f"cannot resolve repository HEAD: {error}") from error
    if not GIT_SHA_PATTERN.fullmatch(value):
        raise CheckpointError("repository HEAD is not a full lowercase Git SHA")
    return value


def exact_int(value: Any, name: str, minimum: int, maximum: int) -> int:
    if type(value) is not int or value < minimum or value > maximum:
        raise CheckpointError(
            f"{name} must be an integer in [{minimum}, {maximum}]"
        )
    return value


def finite_number(
    value: Any, name: str, minimum: float, maximum: float
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CheckpointError(f"{name} must be numeric")
    result = float(value)
    if not math.isfinite(result) or result < minimum or result > maximum:
        raise CheckpointError(
            f"{name} must be finite and in [{minimum}, {maximum}]"
        )
    return result


def load_result(path: pathlib.Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CheckpointError(f"cannot read {path}: {error}") from error
    if not isinstance(payload, dict):
        raise CheckpointError(f"{path}: result root is not an object")
    if payload.get("schema_version") != RESULT_SCHEMA:
        raise CheckpointError(f"{path}: wrong result schema")
    if payload.get("status") != "pass":
        raise CheckpointError(f"{path}: checkpoint did not pass")
    if payload.get("wire_identity") != "existing padded dyadic parent":
        raise CheckpointError(f"{path}: result changed wire identity")
    return payload


def case_geometry(item: Mapping[str, Any]) -> tuple[Any, ...]:
    return tuple(item.get(key) for key in CASE_PARITY_KEYS[:8])


def case_parity(item: Mapping[str, Any]) -> tuple[Any, ...]:
    return tuple(item.get(key) for key in CASE_PARITY_KEYS)


def validate_case_matrix(
    payload: Mapping[str, Any], path: pathlib.Path
) -> tuple[tuple[Any, ...], ...]:
    exact_int(payload.get("context_thread_count"), f"{path}: context threads", 1, 1)
    shard_bytes = exact_int(payload.get("shard_bytes"), f"{path}: shard bytes", 64, 64)
    pointer_bytes = exact_int(payload.get("pointer_bytes"), f"{path}: pointer bytes", 4, 16)
    staging = exact_int(payload.get("pair_staging_shards"), f"{path}: staging shards", 4, 4)
    exact_int(
        payload.get("production_multiplier_checks"),
        f"{path}: multiplier checks",
        65790,
        65790,
    )
    digest = payload.get("correctness_digest_fnv1a64")
    if not isinstance(digest, str) or not DIGEST_PATTERN.fullmatch(digest):
        raise CheckpointError(f"{path}: malformed correctness digest")
    if digest != EXPECTED_CORRECTNESS_DIGEST:
        raise CheckpointError(f"{path}: correctness digest differs from the fixed matrix")
    cases = payload.get("cases")
    if not isinstance(cases, list) or len(cases) != len(EXPECTED_CASE_MATRIX):
        raise CheckpointError(f"{path}: exact case count differs")
    if not all(isinstance(item, dict) for item in cases):
        raise CheckpointError(f"{path}: a case is not an object")
    if tuple(case_geometry(item) for item in cases) != EXPECTED_CASE_MATRIX:
        raise CheckpointError(f"{path}: exact case matrix differs")

    parity: list[tuple[Any, ...]] = []
    for index, item in enumerate(cases):
        if not isinstance(item, dict):
            raise CheckpointError(f"{path}: case {index} is not an object")
        prefix = f"{path}: case {item.get('name', index)}"
        field_bits = exact_int(item.get("field_bits"), f"{prefix} field", 8, 16)
        if field_bits not in (8, 16):
            raise CheckpointError(f"{prefix}: unsupported field")
        parent = exact_int(item.get("parent_size"), f"{prefix} parent", 1, 65536)
        if parent & (parent - 1):
            raise CheckpointError(f"{prefix}: non-dyadic parent")
        shift = exact_int(item.get("shift"), f"{prefix} shift", 0, 65535)
        if shift & (parent - 1) or shift + parent > (1 << field_bits):
            raise CheckpointError(f"{prefix}: invalid aligned coset")
        requested = exact_int(item.get("requested"), f"{prefix} requested", 0, parent)
        exact_int(item.get("active"), f"{prefix} active", 0, parent)
        compared = exact_int(
            item.get("compared_bytes"), f"{prefix} compared bytes", 0, parent * 64
        )
        if compared != requested * shard_bytes:
            raise CheckpointError(f"{prefix}: compared-byte count is inconsistent")
        pairs = exact_int(
            item.get("pair_operations"), f"{prefix} pair operations", 0, parent * 32
        )
        complete = exact_int(
            item.get("complete_operations"), f"{prefix} complete operations", 0, parent
        )
        maximum_complete = exact_int(
            item.get("maximum_complete_block"),
            f"{prefix} maximum complete block",
            0,
            parent,
        )
        if maximum_complete and maximum_complete & (maximum_complete - 1):
            raise CheckpointError(f"{prefix}: complete block is not dyadic")
        if bool(complete) != bool(maximum_complete):
            raise CheckpointError(f"{prefix}: complete block accounting differs")
        if pairs == 0 and complete == 0 and (item.get("active") or requested):
            raise CheckpointError(f"{prefix}: nonempty case has no operations")
        slots = exact_int(
            item.get("scratch_slots"), f"{prefix} scratch slots", 0, 2 * parent
        )
        peak = exact_int(
            item.get("peak_live_slots"), f"{prefix} peak slots", 0, 2 * parent
        )
        if slots != peak:
            raise CheckpointError(f"{prefix}: allocator slot/peak accounting differs")
        execution_scratch = exact_int(
            item.get("execution_scratch_bytes"),
            f"{prefix} execution scratch",
            0,
            2**40,
        )
        expected_scratch = (
            slots * shard_bytes
            + staging * shard_bytes
            + maximum_complete * pointer_bytes
        )
        if execution_scratch != expected_scratch:
            raise CheckpointError(f"{prefix}: execution scratch bound is inconsistent")
        padded_scratch = exact_int(
            item.get("padded_scratch_bytes"), f"{prefix} padded scratch", 1, 2**40
        )
        if padded_scratch != parent * (shard_bytes + pointer_bytes):
            raise CheckpointError(f"{prefix}: padded scratch bound is inconsistent")
        exact_int(
            item.get("serialized_schedule_bytes"),
            f"{prefix} serialized schedule",
            1,
            parent * 4096,
        )
        exact_int(
            item.get("resident_plan_bytes_excluding_allocator_headers"),
            f"{prefix} resident plan",
            1,
            parent * 8192,
        )
        parity.append(case_parity(item))
    return tuple(parity)


def validate_binding(
    payload: Mapping[str, Any],
    path: pathlib.Path,
    expected_source_sha: str,
    expected_core_sha: str,
    expected_core_matrix_sha: str,
    expected_library_sha: str,
    expected_sanitizer: str,
) -> None:
    binding = payload.get("build_binding")
    if not isinstance(binding, dict):
        raise CheckpointError(f"{path}: missing build binding")
    source_sha = binding.get("source_sha256")
    core_sha = binding.get("core_git_sha")
    core_matrix_sha = binding.get("core_matrix_sha256")
    library_sha = binding.get("linked_library_sha256")
    sanitizer = binding.get("sanitizer_mode")
    compiler = binding.get("compiler")
    if source_sha != expected_source_sha or not SHA256_PATTERN.fullmatch(str(source_sha)):
        raise CheckpointError(f"{path}: result is not bound to the selected source")
    if core_sha != expected_core_sha or not GIT_SHA_PATTERN.fullmatch(str(core_sha)):
        raise CheckpointError(f"{path}: result is not bound to the selected core revision")
    if (
        core_matrix_sha != expected_core_matrix_sha
        or not SHA256_PATTERN.fullmatch(str(core_matrix_sha))
    ):
        raise CheckpointError(f"{path}: result is not bound to the selected core matrix")
    if library_sha != expected_library_sha or not SHA256_PATTERN.fullmatch(str(library_sha)):
        raise CheckpointError(f"{path}: result is not bound to the selected library")
    if sanitizer != expected_sanitizer:
        raise CheckpointError(f"{path}: sanitizer build binding differs")
    if not isinstance(compiler, str) or not compiler or compiler == "unbound":
        raise CheckpointError(f"{path}: compiler build binding is missing")


def validate_runtime_environment(
    payload: Mapping[str, Any],
    path: pathlib.Path,
    *,
    require_single_cpu: bool = False,
) -> tuple[str, str, str]:
    environment = payload.get("runtime_environment")
    if not isinstance(environment, dict):
        raise CheckpointError(f"{path}: missing raw runtime environment")
    if environment.get("affinity_api") != "sched_getaffinity":
        raise CheckpointError(f"{path}: sched_getaffinity evidence is unavailable")
    cpus = environment.get("process_affinity_cpus")
    if not isinstance(cpus, list) or not cpus:
        raise CheckpointError(f"{path}: process affinity is empty or malformed")
    checked_cpus = [
        exact_int(cpu, f"{path}: affinity CPU", 0, 1_048_575) for cpu in cpus
    ]
    if checked_cpus != sorted(set(checked_cpus)):
        raise CheckpointError(f"{path}: affinity CPUs are not sorted and unique")
    if require_single_cpu and len(checked_cpus) != 1:
        raise CheckpointError(f"{path}: benchmark was not actually pinned to one CPU")
    if environment.get("omp_num_threads_env") != "1":
        raise CheckpointError(f"{path}: OMP_NUM_THREADS raw evidence is not 1")
    exact_int(
        environment.get("openmp_max_threads"),
        f"{path}: OpenMP maximum threads",
        1,
        1,
    )
    identity = []
    for key in ("hostname", "machine", "cpu_model"):
        value = environment.get(key)
        if not isinstance(value, str) or not value or value == "unavailable":
            raise CheckpointError(f"{path}: raw {key} evidence is unavailable")
        identity.append(value)
    return tuple(identity)  # type: ignore[return-value]


def validate_benchmarks(
    benchmark: Mapping[str, Any],
    path: pathlib.Path,
    case_by_name: Mapping[str, Mapping[str, Any]],
) -> list[float]:
    rows = benchmark.get("benchmarks")
    if not isinstance(rows, list) or len(rows) != len(EXPECTED_BENCHMARK_MATRIX):
        raise CheckpointError(f"{path}: exact benchmark row count differs")
    geometry = tuple(
        (item.get("name"), item.get("field_bits"), item.get("parent_size"), item.get("direction"))
        if isinstance(item, dict) else ()
        for item in rows
    )
    if geometry != EXPECTED_BENCHMARK_MATRIX:
        raise CheckpointError(f"{path}: exact benchmark matrix differs")
    slowdowns: list[float] = []
    for item in rows:
        name = item["name"]
        prefix = f"{path}: benchmark {name}"
        exact_int(item.get("repetitions"), f"{prefix} repetitions", 3, 1_000_000)
        setup = finite_number(item.get("setup_median_us"), f"{prefix} setup median", 0.001, 1e9)
        setup_mad = finite_number(item.get("setup_mad_us"), f"{prefix} setup MAD", 0.0, 1e9)
        candidate = finite_number(item.get("candidate_median_us"), f"{prefix} candidate median", 0.001, 1e9)
        candidate_mad = finite_number(item.get("candidate_mad_us"), f"{prefix} candidate MAD", 0.0, 1e9)
        padded = finite_number(item.get("padded_median_us"), f"{prefix} padded median", 0.001, 1e9)
        padded_mad = finite_number(item.get("padded_mad_us"), f"{prefix} padded MAD", 0.0, 1e9)
        if setup_mad > setup * 10 or candidate_mad > candidate * 10 or padded_mad > padded * 10:
            raise CheckpointError(f"{prefix}: MAD is not credible")
        slowdown = candidate / padded
        if slowdown < MIN_RUNTIME_RATIO or slowdown > MAX_RUNTIME_RATIO:
            raise CheckpointError(f"{prefix}: runtime ratio is outside the evidence bound")
        if slowdown < 1.0:
            raise CheckpointError(f"{prefix}: retained negative-result row is not slower")
        reported = finite_number(
            item.get("padded_over_candidate"), f"{prefix} reported ratio", 0.0, 1e6
        )
        expected = padded / candidate
        if abs(reported - expected) > max(0.00051, expected * 0.02):
            raise CheckpointError(f"{prefix}: reported ratio differs from medians")
        case = case_by_name[name]
        for field in ("candidate_scratch_bytes", "padded_scratch_bytes"):
            value = exact_int(item.get(field), f"{prefix} {field}", 1, 2**40)
            case_field = (
                "execution_scratch_bytes"
                if field == "candidate_scratch_bytes"
                else "padded_scratch_bytes"
            )
            if value != case[case_field]:
                raise CheckpointError(f"{prefix}: benchmark scratch differs from case")
        slowdowns.append(slowdown)
    return slowdowns


def merge(
    backend_paths: Sequence[pathlib.Path],
    sanitizer_path: pathlib.Path,
    benchmark_path: pathlib.Path,
    source_path: pathlib.Path,
    library_paths: Mapping[str, pathlib.Path],
    sanitizer_library_path: pathlib.Path,
    core_git_sha: str,
    core_matrix_sha256: str,
    repository_path: pathlib.Path,
) -> dict[str, Any]:
    if not GIT_SHA_PATTERN.fullmatch(core_git_sha):
        raise CheckpointError("core Git SHA must be exactly 40 lowercase hex digits")
    if not SHA256_PATTERN.fullmatch(core_matrix_sha256):
        raise CheckpointError("core matrix SHA-256 must be exactly 64 lowercase hex digits")
    if set(library_paths) != set(EXPECTED_BACKENDS):
        raise CheckpointError("library bindings must name auto, scalar, SSSE3, and AVX2")
    repository = repository_path.resolve()
    expected_source = (repository / EXPECTED_SOURCE_RELATIVE).resolve()
    if source_path.resolve() != expected_source:
        raise CheckpointError("selected source is not the repository C2 checkpoint source")
    observed_head = repository_head(repository)
    if core_git_sha != observed_head:
        raise CheckpointError("core Git SHA differs from the repository HEAD")
    source_sha = sha256(expected_source)
    library_hashes = {label: sha256(path) for label, path in library_paths.items()}
    sanitizer_library_sha = sha256(sanitizer_library_path)

    if len(backend_paths) != len(EXPECTED_BACKENDS):
        raise CheckpointError("exactly auto, scalar, SSSE3, and AVX2 are required")
    backend_results: dict[str, tuple[pathlib.Path, dict[str, Any]]] = {}
    for path in backend_paths:
        payload = load_result(path)
        label = payload.get("requested_backend")
        if label in backend_results:
            raise CheckpointError(f"duplicate backend result: {label}")
        backend_results[label] = (path, payload)
    if set(backend_results) != set(EXPECTED_BACKENDS):
        raise CheckpointError(f"backend labels differ: {sorted(backend_results)}")

    reference_parity: tuple[tuple[Any, ...], ...] | None = None
    reference_host: tuple[str, str, str] | None = None
    digests: set[str] = set()
    for label, (path, payload) in backend_results.items():
        if payload.get("mode") != "correctness":
            raise CheckpointError(f"{path}: backend result is not correctness-only")
        if payload.get("runtime_backend") != EXPECTED_BACKENDS[label]:
            raise CheckpointError(
                f"{path}: runtime backend {payload.get('runtime_backend')} "
                f"does not match {EXPECTED_BACKENDS[label]}"
            )
        validate_binding(
            payload,
            path,
            source_sha,
            core_git_sha,
            core_matrix_sha256,
            library_hashes[label],
            "none",
        )
        host = validate_runtime_environment(payload, path)
        if reference_host is None:
            reference_host = host
        elif host != reference_host:
            raise CheckpointError(f"{path}: raw host/machine identity differs")
        parity = validate_case_matrix(payload, path)
        if reference_parity is None:
            reference_parity = parity
        elif parity != reference_parity:
            raise CheckpointError(f"{path}: backend scratch/operation/schedule parity differs")
        digests.add(payload["correctness_digest_fnv1a64"])
    if len(digests) != 1:
        raise CheckpointError("backend output digests differ")
    if reference_parity is None or reference_host is None:
        raise CheckpointError("backend reference evidence is missing")
    first = backend_results["auto"][1]

    sanitizer = load_result(sanitizer_path)
    if sanitizer.get("requested_backend") != "asan-ubsan":
        raise CheckpointError("sanitizer result is not labeled as ASan/UBSan")
    if sanitizer.get("mode") != "correctness":
        raise CheckpointError("sanitizer result is not correctness-only")
    if sanitizer.get("runtime_backend") != first.get("runtime_backend"):
        raise CheckpointError("sanitizer runtime backend differs from auto")
    validate_binding(
        sanitizer,
        sanitizer_path,
        source_sha,
        core_git_sha,
        core_matrix_sha256,
        sanitizer_library_sha,
        "asan-ubsan",
    )
    if validate_runtime_environment(sanitizer, sanitizer_path) != reference_host:
        raise CheckpointError("sanitizer raw host/machine identity differs")
    if validate_case_matrix(sanitizer, sanitizer_path) != reference_parity:
        raise CheckpointError("sanitizer scratch/operation/schedule parity differs")
    if sanitizer.get("correctness_digest_fnv1a64") not in digests:
        raise CheckpointError("sanitizer output digest differs")

    benchmark = load_result(benchmark_path)
    if benchmark.get("requested_backend") != "auto":
        raise CheckpointError("benchmark must use the auto backend")
    if benchmark.get("runtime_backend") != first.get("runtime_backend"):
        raise CheckpointError("benchmark runtime backend differs from auto")
    if benchmark.get("mode") != "all":
        raise CheckpointError("benchmark result is not timing mode")
    validate_binding(
        benchmark,
        benchmark_path,
        source_sha,
        core_git_sha,
        core_matrix_sha256,
        library_hashes["auto"],
        "none",
    )
    benchmark_host = validate_runtime_environment(
        benchmark, benchmark_path, require_single_cpu=True
    )
    if benchmark_host != reference_host:
        raise CheckpointError("benchmark raw host/machine identity differs")
    if validate_case_matrix(benchmark, benchmark_path) != reference_parity:
        raise CheckpointError("benchmark scratch/operation/schedule parity differs")
    if benchmark.get("correctness_digest_fnv1a64") not in digests:
        raise CheckpointError("benchmark output digest differs")
    case_by_name = {item["name"]: item for item in benchmark["cases"]}
    slowdowns = validate_benchmarks(benchmark, benchmark_path, case_by_name)

    scratch_ratios = [
        item["execution_scratch_bytes"] / item["padded_scratch_bytes"]
        for item in first["cases"]
    ]
    prefix_ratios = [
        ratio
        for ratio, item in zip(scratch_ratios, first["cases"])
        if "irregular" not in item["name"]
    ]
    irregular_ratios = [
        ratio
        for ratio, item in zip(scratch_ratios, first["cases"])
        if "irregular" in item["name"]
    ]

    return {
        "schema_version": SCHEMA,
        "status": "pass",
        "disposition": "not_promoted",
        "wire_identity": "existing padded dyadic parent",
        "negative_result": (
            "The scalar ragged-boundary executor is {:.1f}-{:.1f}x slower "
            "than the production padded SIMD transform in retained cells; "
            "compact scratch alone is not a promotion argument."
        ).format(min(slowdowns), max(slowdowns)),
        "inconclusive": (
            "A fused SIMD boundary kernel and a more compact generated "
            "schedule were not implemented, so their crossover remains open."
        ),
        "correctness": {
            "digest_fnv1a64": next(iter(digests)),
            "cases_per_backend": len(first["cases"]),
            "bytes_compared_per_backend": sum(
                item["compared_bytes"] for item in first["cases"]
            ),
            "production_multiplier_checks_per_backend": first[
                "production_multiplier_checks"
            ],
            "backends": {
                label: {
                    "runtime_backend": payload["runtime_backend"],
                    "result_sha256": sha256(path),
                    "linked_library_sha256": library_hashes[label],
                }
                for label, (path, payload) in sorted(backend_results.items())
            },
            "asan_ubsan": {
                "status": "pass",
                "runtime_backend": sanitizer["runtime_backend"],
                "result_sha256": sha256(sanitizer_path),
                "linked_library_sha256": sanitizer_library_sha,
            },
        },
        "memory": {
            "accounting": (
                "caller inputs/outputs and ordinary call frames excluded; "
                "arena slots plus four 64-byte pair-staging shards and "
                "maximum complete-block pointers"
            ),
            "prefix_scratch_ratio_min": min(prefix_ratios),
            "prefix_scratch_ratio_max": max(prefix_ratios),
            "irregular_scratch_ratio_min": min(irregular_ratios),
            "irregular_scratch_ratio_max": max(irregular_ratios),
            "largest_serialized_schedule_bytes": max(
                item["serialized_schedule_bytes"] for item in first["cases"]
            ),
            "largest_execution_scratch_bytes": max(
                item["execution_scratch_bytes"] for item in first["cases"]
            ),
            "largest_padded_scratch_bytes": max(
                item["padded_scratch_bytes"] for item in first["cases"]
            ),
        },
        "benchmark": {
            "scope": (
                "single process; raw OMP_NUM_THREADS/openmp_max_threads are 1; "
                "raw sched_getaffinity contains one CPU; medians and MAD; "
                "input packing/output scatter included"
            ),
            "runtime_backend": benchmark["runtime_backend"],
            "slowdown_min": min(slowdowns),
            "slowdown_max": max(slowdowns),
            "rows": benchmark["benchmarks"],
            "result_sha256": sha256(benchmark_path),
        },
        "source": {
            "path": str(source_path),
            "sha256": source_sha,
            "core_git_sha": core_git_sha,
            "core_matrix_sha256": core_matrix_sha256,
        },
        "machine": {
            "hostname": benchmark_host[0],
            "platform": benchmark_host[1],
            "processor": benchmark_host[2],
            "benchmark_process_affinity_cpus": benchmark[
                "runtime_environment"
            ]["process_affinity_cpus"],
            "merge_process_affinity_cpus": sorted(os.sched_getaffinity(0)),
            "python": platform.python_version(),
        },
    }


def parse_binding(value: str) -> tuple[str, pathlib.Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("library binding must be label=path")
    label, path = value.split("=", 1)
    if not label or not path:
        raise argparse.ArgumentTypeError("library binding must be label=path")
    return label, pathlib.Path(path)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", action="append", type=pathlib.Path, required=True)
    parser.add_argument("--library", action="append", type=parse_binding, required=True)
    parser.add_argument("--sanitizer", type=pathlib.Path, required=True)
    parser.add_argument("--sanitizer-library", type=pathlib.Path, required=True)
    parser.add_argument("--benchmark", type=pathlib.Path, required=True)
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--core-git-sha", required=True)
    parser.add_argument("--core-matrix-sha256", required=True)
    parser.add_argument("--repository", type=pathlib.Path, default=pathlib.Path("."))
    parser.add_argument("--output", type=pathlib.Path, required=True)
    arguments = parser.parse_args(argv)
    try:
        libraries: dict[str, pathlib.Path] = {}
        for label, path in arguments.library:
            if label in libraries:
                raise CheckpointError(f"duplicate library binding: {label}")
            libraries[label] = path
        payload = merge(
            arguments.backend,
            arguments.sanitizer,
            arguments.benchmark,
            arguments.source,
            libraries,
            arguments.sanitizer_library,
            arguments.core_git_sha,
            arguments.core_matrix_sha256,
            arguments.repository,
        )
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(
            "C2 checkpoint merged: backends=4 cases/backend={} digest={}".format(
                payload["correctness"]["cases_per_backend"],
                payload["correctness"]["digest_fnv1a64"],
            )
        )
        return 0
    except CheckpointError as error:
        print(f"C2 checkpoint error: {error}", file=os.sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
