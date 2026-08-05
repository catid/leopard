#!/usr/bin/env python3
"""Measure the production small-dual decoder against transform and Leopard1.

This deliberately narrow runner is for final policy tuning, not release-wide
benchmark provenance.  It executes caller-provided, read-only binaries through
already-open descriptors while holding the canonical GF8 campaign lock.  Every
ABBA attempt is pinned to one logical CPU and is rejected if the non-idle
counter of its SMT sibling changes.  Reusable-plan and public one-shot modes
are reported separately.  Full benchmark JSON, workload digests, binary
identities, and before/after hashes are retained for independent review.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import re
import resource
import signal
import shlex
import stat
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable


SCHEMA = "leopard2-small-dual-promotion/v1"
ORDER = ("baseline", "candidate", "candidate", "baseline")
REVERSE_ORDER = ("candidate", "baseline", "baseline", "candidate")
CANONICAL_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")
MAX_CAPTURE_BYTES = 32 << 10
MAX_CELLS = 32
MAX_BINARY_BYTES = 128 << 20
MAX_BINARY_AGGREGATE_BYTES = 256 << 20
F_ADD_SEALS = getattr(fcntl, "F_ADD_SEALS", 1033)
F_GET_SEALS = getattr(fcntl, "F_GET_SEALS", 1034)
REQUIRED_SEALS = (
    getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
    getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
    getattr(fcntl, "F_SEAL_GROW", 0x0004) |
    getattr(fcntl, "F_SEAL_WRITE", 0x0008)
)
CHILD_ENV = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PROC_BIND": "TRUE",
    "OMP_THREAD_LIMIT": "1",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def strict_json_loads(encoded: bytes, label: str) -> Any:
    try:
        text = encoded.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise EvidenceError(f"{label} is not strict UTF-8") from error

    def reject_constant(value: str) -> None:
        raise EvidenceError(f"{label} contains non-finite JSON number {value}")

    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            require(key not in result, f"{label} contains duplicate key {key}")
            result[key] = value
        return result

    try:
        return json.loads(
            text, parse_constant=reject_constant,
            object_pairs_hook=unique_object)
    except json.JSONDecodeError as error:
        raise EvidenceError(f"{label} is invalid JSON") from error


def sha256_descriptor(descriptor: int) -> str:
    digest = hashlib.sha256()
    offset = 0
    while True:
        block = os.pread(descriptor, 1 << 20, offset)
        if not block:
            break
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def descriptor_identity(
    descriptor: int, resolved: Path, expected: os.stat_result | None = None,
) -> dict[str, Any]:
    status = os.fstat(descriptor)
    require(stat.S_ISREG(status.st_mode),
            f"benchmark binary is not regular: {resolved}")
    require(0 < status.st_size <= MAX_BINARY_BYTES,
            f"benchmark binary has an unsafe size: {resolved}")
    require(status.st_mode & 0o222 == 0,
            f"benchmark binary is writable: {resolved}")
    if expected is not None:
        require((status.st_dev, status.st_ino, status.st_size,
                 status.st_mtime_ns) ==
                (expected.st_dev, expected.st_ino, expected.st_size,
                 expected.st_mtime_ns),
                f"benchmark descriptor changed: {resolved}")
    pathname_status = resolved.stat()
    require((pathname_status.st_dev, pathname_status.st_ino,
             pathname_status.st_size, pathname_status.st_mtime_ns) ==
            (status.st_dev, status.st_ino, status.st_size,
             status.st_mtime_ns),
            f"benchmark pathname no longer names its descriptor: {resolved}")
    return {
        "path": str(resolved),
        "size": status.st_size,
        "mode": oct(status.st_mode & 0o7777),
        "device": status.st_dev,
        "inode": status.st_ino,
        "mtime_ns": status.st_mtime_ns,
        "sha256": sha256_descriptor(descriptor),
    }


def open_binary(path: Path) -> tuple[int, dict[str, Any]]:
    resolved = path.resolve(strict=True)
    descriptor = os.open(
        resolved, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        return descriptor, descriptor_identity(descriptor, resolved)
    except Exception:
        os.close(descriptor)
        raise


def sealed_binary_snapshot(
    source_descriptor: int, name: str,
) -> tuple[int, dict[str, Any]]:
    require(hasattr(os, "memfd_create"), "sealed memfd is unavailable")
    descriptor = os.memfd_create(
        f"leopard2-small-dual-{name}",
        getattr(os, "MFD_CLOEXEC", 0x0001) |
        getattr(os, "MFD_ALLOW_SEALING", 0x0002))
    try:
        offset = 0
        while True:
            block = os.pread(source_descriptor, 1 << 20, offset)
            if not block:
                break
            view = memoryview(block)
            while view:
                written = os.write(descriptor, view)
                require(written > 0, "short sealed-binary write")
                view = view[written:]
            offset += len(block)
        os.fchmod(descriptor, 0o555)
        fcntl.fcntl(descriptor, F_ADD_SEALS, REQUIRED_SEALS)
        seals = fcntl.fcntl(descriptor, F_GET_SEALS)
        require(seals & REQUIRED_SEALS == REQUIRED_SEALS,
                "benchmark snapshot is not write-sealed")
        status = os.fstat(descriptor)
        require(stat.S_ISREG(status.st_mode) and status.st_size == offset,
                "sealed benchmark snapshot differs in size")
        identity = {
            "name": name,
            "size": status.st_size,
            "mode": oct(status.st_mode & 0o7777),
            "device": status.st_dev,
            "inode": status.st_ino,
            "seals": seals,
            "sha256": sha256_descriptor(descriptor),
        }
        return descriptor, identity
    except Exception:
        os.close(descriptor)
        raise


def canonical_lock_open() -> tuple[int, dict[str, Any]]:
    descriptor = os.open(
        CANONICAL_LOCK,
        os.O_RDWR | os.O_CREAT | os.O_CLOEXEC | os.O_NOFOLLOW,
        0o600)
    try:
        status = os.fstat(descriptor)
        require(stat.S_ISREG(status.st_mode) and status.st_uid == os.getuid() and
                status.st_nlink == 1 and status.st_mode & 0o777 == 0o600,
                "canonical benchmark lock metadata is invalid")
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        pathname_status = CANONICAL_LOCK.stat()
        require((pathname_status.st_dev, pathname_status.st_ino) ==
                (status.st_dev, status.st_ino),
                "canonical benchmark lock pathname changed")
        return descriptor, {
            "path": str(CANONICAL_LOCK),
            "device": status.st_dev,
            "inode": status.st_ino,
            "uid": status.st_uid,
            "mode": oct(status.st_mode & 0o7777),
            "nlink": status.st_nlink,
        }
    except Exception:
        os.close(descriptor)
        raise


def validate_lock_identity(
    descriptor: int, identity: dict[str, Any],
) -> None:
    status = os.fstat(descriptor)
    pathname_status = CANONICAL_LOCK.stat()
    expected = (identity["device"], identity["inode"])
    require((status.st_dev, status.st_ino) == expected and
            (pathname_status.st_dev, pathname_status.st_ino) == expected and
            stat.S_ISREG(status.st_mode) and status.st_uid == os.getuid() and
            status.st_nlink == 1 and status.st_mode & 0o777 == 0o600,
            "canonical benchmark lock identity changed")


def read_bounded(path: Path, limit: int, label: str) -> bytes:
    resolved = path.resolve(strict=True)
    status_before = resolved.stat()
    require(stat.S_ISREG(status_before.st_mode) and
            status_before.st_size <= limit,
            f"{label} is absent, non-regular, or too large")
    data = resolved.read_bytes()
    status_after = resolved.stat()
    require(len(data) == status_before.st_size and
            (status_before.st_dev, status_before.st_ino,
             status_before.st_size, status_before.st_mtime_ns) ==
            (status_after.st_dev, status_after.st_ino,
             status_after.st_size, status_after.st_mtime_ns),
            f"{label} changed while read")
    return data


def build_configuration_identity(
    build: Path, selector: str,
) -> dict[str, Any]:
    resolved = build.resolve(strict=True)
    require(resolved.is_dir(), "build path is not a directory")
    cache_path = resolved / "CMakeCache.txt"
    commands_path = resolved / "compile_commands.json"
    executable_path = resolved / "bench_leopard2"
    cache_bytes = read_bounded(cache_path, 4 << 20, "CMake cache")
    commands_bytes = read_bounded(
        commands_path, 32 << 20, "compile-command database")
    executable_bytes = read_bounded(
        executable_path, MAX_BINARY_BYTES, "build benchmark executable")
    try:
        cache_text = cache_bytes.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise EvidenceError("CMake cache is not strict UTF-8") from error
    expected_cache = f"LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT:BOOL={selector}"
    require(cache_text.splitlines().count(expected_cache) == 1,
            "build cache has the wrong small-dual selector")
    commands = strict_json_loads(commands_bytes, "compile-command database")
    require(isinstance(commands, list),
            "compile-command database is not an array")
    entries = [entry for entry in commands if isinstance(entry, dict) and
               isinstance(entry.get("file"), str) and
               Path(entry["file"]).name == "leopard2.cpp"]
    require(len(entries) == 1 and isinstance(entries[0].get("command"), str),
            "portable Leopard2 compile action is not unique")
    tokens = shlex.split(entries[0]["command"], posix=True)
    expected_definition = (
        "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=" +
        ("1" if selector == "ON" else "0"))
    opposite_definition = (
        "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=" +
        ("0" if selector == "ON" else "1"))
    require(tokens.count(expected_definition) == 1 and
            opposite_definition not in tokens,
            "portable Leopard2 compile action has the wrong selector")
    cache_entries: dict[str, tuple[str, str]] = {}
    for line in cache_text.splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        left, value = line.split("=", 1)
        if ":" not in left:
            continue
        key, kind = left.split(":", 1)
        require(key not in cache_entries, "CMake cache key is duplicated")
        cache_entries[key] = (kind, value)
    require(cache_entries.get("LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT") ==
            ("BOOL", selector), "parsed selector differs")
    normalized_cache = {}
    for key, (kind, value) in cache_entries.items():
        value = value.replace(str(resolved), "<BUILD>")
        if key == "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT":
            value = "<SELECTOR>"
        elif key == "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
            value = "<CONFIGURATION-SHA256>"
        normalized_cache[key] = (kind, value)
    configuration_hash = cache_entries[
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"][1]
    require(re.fullmatch(r"[0-9a-f]{64}", configuration_hash) is not None,
            "effective configuration hash is invalid")

    def normalize_compile_value(value: Any) -> Any:
        if isinstance(value, str):
            return value.replace(str(resolved), "<BUILD>").replace(
                expected_definition,
                "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=<SELECTOR>"
            ).replace(configuration_hash, "<CONFIGURATION-SHA256>")
        if isinstance(value, list):
            return [normalize_compile_value(item) for item in value]
        if isinstance(value, dict):
            return {
                key: normalize_compile_value(item)
                for key, item in value.items()
            }
        return value

    require(all(isinstance(entry, dict) for entry in commands),
            "compile-command database contains a non-object entry")
    normalized_commands = [normalize_compile_value(entry)
                           for entry in commands]
    normalized_commands.sort(key=lambda entry: json.dumps(
        entry, sort_keys=True, separators=(",", ":")))
    normalized_commands_bytes = json.dumps(
        normalized_commands, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return {
        "path": str(resolved),
        "selector": selector,
        "cache_sha256": hashlib.sha256(cache_bytes).hexdigest(),
        "compile_commands_sha256": hashlib.sha256(commands_bytes).hexdigest(),
        "benchmark_sha256": hashlib.sha256(executable_bytes).hexdigest(),
        "normalized_cache_sha256": hashlib.sha256(json.dumps(
            normalized_cache, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")).hexdigest(),
        "normalized_compile_commands_sha256": hashlib.sha256(
            normalized_commands_bytes).hexdigest(),
        "leopard2_source": str(Path(entries[0]["file"]).resolve(strict=True)),
        "leopard2_compile_definition": expected_definition,
    }


def cpu_counters(cpu: int) -> dict[str, Any]:
    prefix = f"cpu{cpu} "
    with Path("/proc/stat").open("r", encoding="ascii") as source:
        for line in source:
            if line.startswith(prefix):
                values = tuple(int(value) for value in line.split()[1:])
                require(len(values) >= 8, "short per-CPU /proc/stat record")
                # idle and iowait do not execute on the sibling.  guest time is
                # already included in user/nice, so do not count it twice.
                non_idle = sum(values[index] for index in (0, 1, 2, 5, 6, 7))
                return {"raw": values, "non_idle": non_idle}
    raise EvidenceError(f"CPU {cpu} is absent from /proc/stat")


def child_setup(cpu: int) -> None:
    os.sched_setaffinity(0, {cpu})
    resource.setrlimit(resource.RLIMIT_AS, (256 << 20, 256 << 20))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(
        resource.RLIMIT_FSIZE, (MAX_CAPTURE_BYTES, MAX_CAPTURE_BYTES))


def benchmark_arguments(
    implementation: str,
    cell: dict[str, int | str],
    iterations: int,
    warmup: int,
    decode_mode: str,
    source_attestation: bool = False,
) -> list[str]:
    common = [
        "--k", str(cell["k"]),
        "--r", str(cell["r"]),
        "--loss", str(cell["losses"]),
        "--batch", "1",
        "--reuse", str(cell["reuse"]),
        "--iterations", str(iterations),
        "--warmup", str(warmup),
        "--threads", "1",
        "--seed", str(cell["seed"]),
        "--json", "-",
    ]
    logical_bytes = int(cell["bytes"])
    if implementation == "main":
        physical_bytes = (logical_bytes + 63) & ~63
        return [
            "--bytes", str(physical_bytes),
            "--logical-bytes", str(logical_bytes),
            *common,
        ]
    if source_attestation:
        diagnostics = ["--attest-source"]
    elif decode_mode == "one_shot":
        diagnostics = ["--retain-samples", "--measure-one-shot-decode"]
    else:
        diagnostics = [
            "--retain-samples",
            "--report-decode-path",
            "--report-direct-executor",
        ]
    return [
        "--profile", "high",
        "--field", "gf8",
        "--backend", "avx2",
        "--bytes", str(logical_bytes),
        "--skip-legacy",
        *diagnostics,
        *common,
    ]


def run_benchmark(
    descriptor: int,
    lock_descriptor: int,
    implementation: str,
    cell: dict[str, int | str],
    cpu: int,
    iterations: int,
    warmup: int,
    timeout: float,
    decode_mode: str,
    source_attestation: bool = False,
) -> dict[str, Any]:
    executable = f"/proc/self/fd/{descriptor}"
    arguments = benchmark_arguments(
        implementation, cell, iterations, warmup, decode_mode,
        source_attestation)
    started = time.monotonic_ns()
    with tempfile.TemporaryFile() as stdout_file, \
            tempfile.TemporaryFile() as stderr_file:
        process = subprocess.Popen(
            [executable, *arguments],
            executable=executable,
            stdin=subprocess.DEVNULL,
            stdout=stdout_file,
            stderr=stderr_file,
            env=CHILD_ENV,
            # Keep the canonical flock's open-file description alive in every
            # timed child.  If the coordinator dies, another campaign still
            # cannot overlap the surviving executable.
            pass_fds=(descriptor, lock_descriptor),
            preexec_fn=lambda: child_setup(cpu),
            start_new_session=True,
        )
        try:
            return_code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired as error:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            raise EvidenceError(
                f"{implementation} benchmark timed out") from error
        descendant_survived = False
        try:
            os.killpg(process.pid, 0)
            descendant_survived = True
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        require(not descendant_survived,
                f"{implementation} benchmark left a descendant")
        stdout_file.seek(0, os.SEEK_END)
        stderr_file.seek(0, os.SEEK_END)
        stdout_size = stdout_file.tell()
        stderr_size = stderr_file.tell()
        require(stdout_size <= MAX_CAPTURE_BYTES and
                stderr_size <= MAX_CAPTURE_BYTES,
                f"{implementation} benchmark output exceeded its bound")
        stdout_file.seek(0)
        stderr_file.seek(0)
        stdout = stdout_file.read(MAX_CAPTURE_BYTES + 1)
        stderr = stderr_file.read(MAX_CAPTURE_BYTES + 1)
    elapsed = time.monotonic_ns() - started
    require(return_code == 0,
            f"{implementation} benchmark failed ({return_code}): "
            f"{stderr.decode('utf-8', errors='replace')}")
    require(not stderr,
            f"{implementation} benchmark emitted stderr")
    result = strict_json_loads(
        stdout, f"{implementation} benchmark output")
    require(isinstance(result, dict), "benchmark JSON is not an object")
    correctness = result.get("correctness")
    require(isinstance(correctness, dict), "benchmark correctness is absent")
    if implementation == "main":
        require(correctness.get("round_trip") is True,
                "Leopard1 round trip failed")
    else:
        require(correctness.get("leopard2_round_trip") is True,
                "Leopard2 round trip failed")
        resolved = result.get("resolved")
        require(isinstance(resolved, dict) and
                resolved.get("backend") == "avx2" and
                resolved.get("field") == "gf8" and
                resolved.get("profile") == "legacy_high_v1",
                "Leopard2 did not resolve to legacy-high GF8/AVX2")
    return {
        "implementation": implementation,
        "arguments": arguments,
        "elapsed_ns": elapsed,
        "result": result,
    }


def finite_number(value: Any, label: str, positive: bool = False) -> float:
    require(not isinstance(value, bool) and isinstance(value, (int, float)),
            f"{label} is not numeric")
    result = float(value)
    require(math.isfinite(result), f"{label} is not finite")
    require(not positive or result > 0, f"{label} is not positive")
    return result


def validate_summary_samples(
    value: Any,
    sample_key: str,
    median_key: str,
    mad_key: str,
    minimum_key: str,
    maximum_key: str,
    iterations: int,
    label: str,
) -> tuple[list[float], float]:
    require(isinstance(value, dict), f"{label} is not an object")
    samples_value = value.get(sample_key)
    require(isinstance(samples_value, list) and
            len(samples_value) == iterations,
            f"{label} sample count differs")
    samples = [finite_number(sample, f"{label} sample", positive=True)
               for sample in samples_value]
    median = statistics.median(samples)
    mad = statistics.median(abs(sample - median) for sample in samples)
    reported = {
        "median": finite_number(value.get(median_key),
                                f"{label} median", positive=True),
        "mad": finite_number(value.get(mad_key), f"{label} MAD"),
        "minimum": finite_number(value.get(minimum_key),
                                 f"{label} minimum", positive=True),
        "maximum": finite_number(value.get(maximum_key),
                                 f"{label} maximum", positive=True),
    }
    expected = {
        "median": median,
        "mad": mad,
        "minimum": min(samples),
        "maximum": max(samples),
    }
    for key in expected:
        require(math.isclose(reported[key], expected[key],
                             rel_tol=2e-6, abs_tol=2e-6),
                f"{label} {key} was not derived from retained samples")
    return samples, median


def validate_digests(value: Any, label: str) -> dict[str, str]:
    require(isinstance(value, dict) and
            set(value) == {"algorithm", "original_data",
                           "transmitted_parity", "recovered_originals"} and
            value.get("algorithm") == "fnv1a64",
            f"{label} digest object is invalid")
    for key in ("original_data", "transmitted_parity", "recovered_originals"):
        require(isinstance(value[key], str) and
                re.fullmatch(r"[0-9a-f]{16}", value[key]) is not None,
                f"{label} {key} digest is invalid")
    return value


def validate_benchmark_result(
    invocation: dict[str, Any],
    cell: dict[str, int | str],
    iterations: int,
    warmup: int,
    main_commit: str,
    decode_mode: str,
) -> None:
    implementation = invocation["implementation"]
    result = invocation["result"]
    logical_bytes = int(cell["bytes"])
    physical_bytes = (logical_bytes + 63) & ~63
    if implementation == "main":
        expected_schema = ("leopard-main-benchmark-v1" if
                           physical_bytes == logical_bytes else
                           "leopard-main-benchmark-v2")
        require(result.get("schema") == expected_schema,
                "Leopard1 benchmark schema differs")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("main_source_commit") == main_commit,
                "Leopard1 source commit differs")
    else:
        expected_schema = ("leopard2-benchmark-v9" if
                           decode_mode == "one_shot" else
                           "leopard2-benchmark-v6")
        require(result.get("schema") == expected_schema,
                f"{implementation} benchmark schema differs")
    parameters = result.get("parameters")
    require(isinstance(parameters, dict), "benchmark parameters are absent")
    expected = {
        "K": int(cell["k"]), "R": int(cell["r"]),
        "loss_count": int(cell["losses"]), "batch": 1,
        "reuse": int(cell["reuse"]), "iterations": iterations,
        "warmup": warmup, "thread_count": 1, "seed": int(cell["seed"]),
    }
    for key, value in expected.items():
        require(type(parameters.get(key)) is int and
                parameters.get(key) == value,
                f"{implementation} parameter {key} differs")
    missing = parameters.get("missing_original_indices")
    require(isinstance(missing, list) and
            len(missing) == int(cell["losses"]) and
            all(type(index) is int and 0 <= index < int(cell["k"])
                for index in missing) and len(set(missing)) == len(missing),
            f"{implementation} missing-original set is invalid")
    if implementation == "main":
        require(parameters.get("shard_bytes") == physical_bytes and
                parameters.get("logical_shard_bytes") == logical_bytes,
                "Leopard1 physical/logical byte parameters differ")
        resolved = result.get("resolved")
        padded_application_bytes = resolved.get(
            "padded_application_bytes") if isinstance(resolved, dict) else None
        require(isinstance(resolved, dict) and
                resolved.get("profile") == "legacy_high_v1" and
                resolved.get("field") == "gf8" and
                type(padded_application_bytes) is bool and
                padded_application_bytes == (physical_bytes != logical_bytes) and
                resolved.get("padding_policy") == "zero suffix per shard",
                "Leopard1 resolved profile or padding policy differs")
    else:
        require(parameters.get("shard_bytes") == logical_bytes and
                parameters.get("requested_profile") == "legacy_high_v1" and
                parameters.get("requested_field") == "gf8" and
                parameters.get("requested_backend") == "avx2" and
                parameters.get("retain_samples") is True and
                ((decode_mode == "one_shot" and
                  parameters.get("measure_one_shot_decode") is True and
                  "report_direct_executor" not in parameters) or
                 (decode_mode == "reusable" and
                  parameters.get("report_direct_executor") is True and
                  "measure_one_shot_decode" not in parameters)),
                f"{implementation} requested policy differs")
    validate_digests(result.get("workload_digests"), implementation)
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "benchmark metrics are absent")
    if implementation == "main":
        validate_summary_samples(
            metrics.get("decode_including_setup"),
            "samples_us_per_batch_call", "median_us_per_batch_call",
            "mad_us_per_batch_call", "minimum_us_per_batch_call",
            "maximum_us_per_batch_call", iterations,
            "Leopard1 decode including setup")
    else:
        _, setup_median = validate_summary_samples(
            metrics.get("decode_plan_setup"), "samples_us", "median_us",
            "mad_us", "minimum_us", "maximum_us", iterations,
            f"{implementation} decode plan setup")
        _, execution_median = validate_summary_samples(
            metrics.get("decode_execution"),
            "samples_us_per_batch_call", "median_us_per_batch_call",
            "mad_us_per_batch_call", "minimum_us_per_batch_call",
            "maximum_us_per_batch_call", iterations,
            f"{implementation} decode execution")
        amortized = metrics.get("decode_amortized_at_reuse")
        require(isinstance(amortized, dict) and
                amortized.get("reuse_count") == int(cell["reuse"]),
                f"{implementation} amortized reuse differs")
        reported_amortized = finite_number(
            amortized.get("derived_median_us_per_batch_call"),
            f"{implementation} amortized decode", positive=True)
        require(math.isclose(
                    reported_amortized,
                    execution_median + setup_median / int(cell["reuse"]),
                    rel_tol=2e-6, abs_tol=2e-6),
                f"{implementation} amortized decode is not derived")
        if decode_mode == "one_shot":
            validate_summary_samples(
                metrics.get("one_shot_decode_including_setup"),
                "samples_us_per_batch_call", "median_us_per_batch_call",
                "mad_us_per_batch_call", "minimum_us_per_batch_call",
                "maximum_us_per_batch_call", iterations,
                f"{implementation} one-shot decode including setup")


def validate_source_attestations(
    candidate: dict[str, Any], control: dict[str, Any], commit: str,
) -> dict[str, Any]:
    builds: list[dict[str, Any]] = []
    for invocation in (candidate, control):
        result = invocation["result"]
        require(result.get("schema") == "leopard2-benchmark-v5",
                "source-attestation schema differs")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("source_commit") == commit and
                build.get("source_tracked_dirty") is False and
                isinstance(build.get("source_tree"), str) and
                re.fullmatch(r"[0-9a-f]{40}", build["source_tree"]) is not None,
                "Leopard2 source attestation differs")
        builds.append(build)
    require(builds[0]["source_commit"] == builds[1]["source_commit"] and
            builds[0]["source_tree"] == builds[1]["source_tree"],
            "candidate and control source identities differ")
    return {
        "candidate": candidate,
        "control": control,
        "source_commit": builds[0]["source_commit"],
        "source_tree": builds[0]["source_tree"],
        "source_tracked_dirty": False,
    }


def decode_time(invocation: dict[str, Any], metric: str) -> float:
    result = invocation["result"]
    metrics = result["metrics"]
    implementation = invocation["implementation"]
    if implementation == "main":
        require(metric in ("amortized", "one_shot"),
                "invalid Leopard1 metric selector")
        return float(metrics["decode_including_setup"]
                     ["median_us_per_batch_call"])
    if metric == "one_shot":
        return float(metrics["one_shot_decode_including_setup"]
                     ["median_us_per_batch_call"])
    if metric == "amortized":
        return float(metrics["decode_amortized_at_reuse"]
                     ["derived_median_us_per_batch_call"])
    require(metric == "execution", "invalid Leopard2 metric selector")
    return float(metrics["decode_execution"]["median_us_per_batch_call"])


def validate_workload_equality(
    invocations: Iterable[dict[str, Any]], comparison: str,
) -> None:
    retained = list(invocations)
    require(retained, f"{comparison} has no invocations")
    expected_digests = retained[0]["result"]["workload_digests"]
    expected_missing = retained[0]["result"]["parameters"] \
        ["missing_original_indices"]
    for invocation in retained:
        require(invocation["result"]["workload_digests"] == expected_digests,
                f"{comparison} workload digests differ")
        require(invocation["result"]["parameters"]
                    ["missing_original_indices"] == expected_missing,
                f"{comparison} missing-original set differs")


def expected_candidate_direct(cell: dict[str, int | str]) -> bool:
    k = int(cell["k"])
    losses = int(cell["losses"])
    return (5 <= k <= 16 and 5 <= int(cell["r"]) <= 8 and
            5 <= losses <= 8 and int(cell["bytes"]) >= 1024 and
            not (losses == k and k >= 7))


def validate_routes(
    cell: dict[str, int | str], comparison: str,
    invocations: Iterable[dict[str, Any]],
    decode_mode: str,
) -> None:
    # The standalone public one-shot timing schema deliberately omits internal
    # route metadata.  Its exact ON/OFF build selectors and source attestation
    # bind the policy comparison; focused production tests cover route choice.
    if decode_mode == "one_shot":
        return
    for invocation in invocations:
        if invocation["implementation"] == "main":
            continue
        path = invocation["result"]["resolved"].get("selected_decode_path")
        rule = invocation["result"]["resolved"].get("selected_decode_rule")
        executor = invocation["result"]["resolved"].get(
            "selected_direct_executor")
        if int(cell["losses"]) == 4:
            require(path == "direct" and rule == "direct" and
                    executor == "output_major",
                    f"established four-loss route differs for {cell['id']}")
        elif invocation["implementation"] == "candidate":
            require((path == "direct") == expected_candidate_direct(cell),
                    f"candidate route mismatch for {cell['id']}: {path}")
            if path == "direct":
                require(rule == "direct" and executor == "source_major",
                        f"candidate direct executor differs for {cell['id']}")
            else:
                require((path, rule) in (
                            ("tiled", "workspace_tiled"),
                            ("materialized", "workspace_materialized"),
                            ("materialized", "translated_low")) and
                        executor == "none",
                        f"candidate transform route differs for {cell['id']}")
        elif comparison == "same_source":
            require((path, rule) in (
                        ("tiled", "workspace_tiled"),
                        ("materialized", "workspace_materialized"),
                        ("materialized", "translated_low")) and
                    executor == "none",
                    f"transform-only control route differs for {cell['id']}")


def run_attempt(
    comparison: str,
    descriptors: dict[str, int],
    lock_descriptor: int,
    cell: dict[str, int | str],
    cpu: int,
    sibling: int,
    rounds: int,
    iterations: int,
    warmup: int,
    timeout: float,
    main_commit: str,
    decode_mode: str,
) -> dict[str, Any]:
    baseline_implementation = "control" if comparison == "same_source" else "main"
    sibling_before = cpu_counters(sibling)
    rounds_out: list[dict[str, Any]] = []
    all_invocations: list[dict[str, Any]] = []
    for round_index in range(rounds):
        order = ORDER if round_index % 2 == 0 else REVERSE_ORDER
        invocations: list[dict[str, Any]] = []
        for role in order:
            implementation = (
                "candidate" if role == "candidate" else baseline_implementation)
            invocation = run_benchmark(
                descriptors[implementation], lock_descriptor,
                implementation, cell, cpu,
                iterations, warmup, timeout, decode_mode)
            validate_benchmark_result(
                invocation, cell, iterations, warmup, main_commit,
                decode_mode)
            invocations.append(invocation)
            all_invocations.append(invocation)
        pair_indices = ((0, 1), (3, 2)) if order == ORDER else \
            ((1, 0), (2, 3))
        if decode_mode == "one_shot":
            metric_names = ("one_shot",)
        else:
            metric_names = (("execution", "amortized") if
                            comparison == "same_source" else ("amortized",))
        metric_results: dict[str, Any] = {}
        for metric in metric_names:
            paired_log_speedups: list[float] = []
            pairs: list[dict[str, float]] = []
            for baseline_index, candidate_index in pair_indices:
                baseline_us = decode_time(
                    invocations[baseline_index], metric)
                candidate_us = decode_time(
                    invocations[candidate_index], metric)
                log_speedup = math.log(baseline_us / candidate_us)
                paired_log_speedups.append(log_speedup)
                pairs.append({
                    "baseline_us": baseline_us,
                    "candidate_us": candidate_us,
                    "speedup": math.exp(log_speedup),
                })
            metric_results[metric] = {
                "pairs": pairs,
                "round_log_speedup": statistics.mean(paired_log_speedups),
                "round_geometric_mean_speedup": math.exp(
                    statistics.mean(paired_log_speedups)),
            }
        rounds_out.append({
            "round": round_index,
            "order": order,
            "metrics": metric_results,
            "invocations": invocations,
        })
    validate_workload_equality(all_invocations, comparison)
    validate_routes(cell, comparison, all_invocations, decode_mode)
    sibling_after = cpu_counters(sibling)
    sibling_delta = sibling_after["non_idle"] - sibling_before["non_idle"]
    t_critical = {3: 4.302652729911275, 4: 3.182446305284263,
                  5: 2.7764451051977987}[rounds]
    aggregate_metrics: dict[str, Any] = {}
    for metric in rounds_out[0]["metrics"]:
        round_logs = [value["metrics"][metric]["round_log_speedup"]
                      for value in rounds_out]
        mean_log = statistics.mean(round_logs)
        standard_error = statistics.stdev(round_logs) / math.sqrt(rounds)
        aggregate_metrics[metric] = {
            "round_log_speedups": round_logs,
            "geometric_mean_speedup": math.exp(mean_log),
            "confidence_level": 0.95,
            "confidence_interval_lower": math.exp(
                mean_log - t_critical * standard_error),
            "confidence_interval_upper": math.exp(
                mean_log + t_critical * standard_error),
        }
    return {
        "comparison": comparison,
        "sibling_before": sibling_before,
        "sibling_after": sibling_after,
        "sibling_non_idle_delta": sibling_delta,
        "accepted": sibling_delta == 0,
        "metrics": aggregate_metrics,
        "rounds": rounds_out,
    }


def parse_cell(text: str) -> dict[str, int | str]:
    fields = text.split(":")
    require(len(fields) == 7,
            "cell must be ID:K:R:BYTES:LOSSES:SEED:REUSE")
    identifier = fields[0]
    require(identifier and all(character.isalnum() or character in "-_"
                               for character in identifier),
            "cell ID is invalid")
    try:
        k, r, shard_bytes, losses, seed, reuse = (
            int(value) for value in fields[1:])
    except ValueError as error:
        raise EvidenceError("cell contains a non-integer value") from error
    require(5 <= k <= 16 and 5 <= r <= 8 and r <= k,
            "cell is outside the exact-main small-dual count scope")
    require(4 <= losses <= min(k, r), "cell loss count is invalid")
    require(1 <= shard_bytes <= (1 << 20) and
            0 <= seed <= (1 << 64) - 1 and 1 <= reuse <= 4096,
            "cell byte/seed/reuse value is invalid")
    return {"id": identifier, "k": k, "r": r, "bytes": shard_bytes,
            "losses": losses, "seed": seed, "reuse": reuse}


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(value, sort_keys=True, separators=(",", ":")) +
               "\n").encode("utf-8")
    require(len(encoded) <= (128 << 20), "result JSON exceeds its bound")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=str(path.parent))
    linked = False
    committed = False
    try:
        with os.fdopen(descriptor, "wb", closefd=True) as output:
            output.write(encoded)
            output.flush()
            os.fsync(output.fileno())
        try:
            os.link(temporary_name, path, follow_symlinks=False)
            linked = True
        except FileExistsError as error:
            raise EvidenceError(f"output already exists: {path}") from error
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
        os.unlink(temporary_name)
        temporary_name = ""
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
        committed = True
    finally:
        if temporary_name:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
        if linked and not committed:
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass


def machine_record(cpu: int, sibling: int) -> dict[str, Any]:
    result: dict[str, Any] = {
        "allowed_cpus": sorted(os.sched_getaffinity(0)),
        "benchmark_cpu": cpu,
        "reserved_sibling": sibling,
        "uname": tuple(os.uname()),
    }
    for name, path in (
        ("governor", Path(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_governor")),
        ("driver", Path(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_driver")),
        ("siblings", Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")),
    ):
        try:
            result[name] = path.read_text(encoding="ascii").strip()
        except OSError:
            result[name] = None
    return result


def parse_cpu_list(text: str) -> set[int]:
    result: set[int] = set()
    for part in text.strip().split(","):
        require(part != "", "empty CPU-list component")
        if "-" in part:
            first_text, last_text = part.split("-", 1)
            first, last = int(first_text), int(last_text)
            require(0 <= first <= last, "invalid CPU-list range")
            result.update(range(first, last + 1))
        else:
            value = int(part)
            require(value >= 0, "negative CPU-list member")
            result.add(value)
    return result


def validate_supervised_pair(cpu: int, sibling: int) -> str:
    nonce = os.environ.get("LEO2_AFFINITY_EXECUTION_NONCE")
    require(isinstance(nonce, str) and
            re.fullmatch(r"[0-9a-f]{64}", nonce) is not None,
            "campaign must run under leopard2_affinity_supervisor.py")
    cpu_siblings = parse_cpu_list(Path(
        f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"
    ).read_text(encoding="ascii"))
    sibling_siblings = parse_cpu_list(Path(
        f"/sys/devices/system/cpu/cpu{sibling}/topology/thread_siblings_list"
    ).read_text(encoding="ascii"))
    require(cpu_siblings == sibling_siblings == {cpu, sibling},
            "declared CPUs are not one exact SMT pair")
    coordinator_cpus = os.sched_getaffinity(0)
    require(cpu not in coordinator_cpus and sibling not in coordinator_cpus,
            "supervisor did not exclude the reserved pair from coordinator")
    require(coordinator_cpus, "coordinator has no housekeeping CPU")
    return nonce


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument("--candidate-build-dir", required=True, type=Path)
    parser.add_argument("--control-build-dir", required=True, type=Path)
    parser.add_argument("--candidate-commit", required=True)
    parser.add_argument("--main-commit", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cell", action="append", required=True)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--decode-mode", choices=("reusable", "one_shot"),
                        default="reusable")
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--iterations", type=int, default=9)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--max-attempts", type=int, default=2)
    options = parser.parse_args()
    require(options.cpu != options.sibling, "CPU pair members must differ")
    require(3 <= options.rounds <= 4 and 5 <= options.iterations <= 31 and
            1 <= options.warmup <= 10 and 1 <= options.max_attempts <= 2 and
            math.isfinite(options.timeout) and 1 <= options.timeout <= 120,
            "campaign counts are below the evidence floor")
    cells = [parse_cell(value) for value in options.cell]
    require(1 <= len(cells) <= MAX_CELLS and
            len({cell["id"] for cell in cells}) == len(cells),
            "duplicate cell ID")
    require(re.fullmatch(r"[0-9a-f]{40}", options.candidate_commit) is not None and
            re.fullmatch(r"[0-9a-f]{40}", options.main_commit) is not None,
            "source commit is not a full hexadecimal object ID")
    execution_nonce = validate_supervised_pair(options.cpu, options.sibling)
    resource.setrlimit(resource.RLIMIT_AS, (768 << 20, 768 << 20))
    paths = {
        "candidate": options.candidate,
        "control": options.control,
        "main": options.main,
    }
    lock_descriptor, lock_identity = canonical_lock_open()
    original_descriptors: dict[str, int] = {}
    snapshot_descriptors: dict[str, int] = {}
    identities_before: dict[str, dict[str, Any]] = {}
    snapshot_identities: dict[str, dict[str, Any]] = {}
    started_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    try:
        validate_lock_identity(lock_descriptor, lock_identity)
        build_identities = {
            "candidate": build_configuration_identity(
                options.candidate_build_dir, "ON"),
            "control": build_configuration_identity(
                options.control_build_dir, "OFF"),
        }
        require(build_identities["candidate"]["normalized_cache_sha256"] ==
                    build_identities["control"]["normalized_cache_sha256"] and
                build_identities["candidate"]
                    ["normalized_compile_commands_sha256"] ==
                    build_identities["control"]
                    ["normalized_compile_commands_sha256"] and
                build_identities["candidate"]["leopard2_source"] ==
                    build_identities["control"]["leopard2_source"],
                "candidate/control normalized build configurations differ")
        for name, path in paths.items():
            descriptor, identity = open_binary(path)
            original_descriptors[name] = descriptor
            identities_before[name] = identity
        require(sum(identity["size"] for identity in
                    identities_before.values()) <=
                MAX_BINARY_AGGREGATE_BYTES,
                "benchmark binaries exceed the aggregate snapshot bound")
        for name, descriptor in original_descriptors.items():
            snapshot, snapshot_identity = sealed_binary_snapshot(
                descriptor, name)
            snapshot_descriptors[name] = snapshot
            snapshot_identities[name] = snapshot_identity
            require(identities_before[name]["sha256"] ==
                    snapshot_identity["sha256"],
                    f"{name} sealed snapshot differs")
        require(build_identities["candidate"]["benchmark_sha256"] ==
                    identities_before["candidate"]["sha256"] and
                build_identities["control"]["benchmark_sha256"] ==
                    identities_before["control"]["sha256"],
                "frozen candidate/control differs from declared build target")
        attestation_cell = dict(cells[0])
        attestation_cell["reuse"] = 1
        source_attestations = validate_source_attestations(
            run_benchmark(
                snapshot_descriptors["candidate"], lock_descriptor,
                "candidate", attestation_cell, options.cpu, 1, 1,
                options.timeout, options.decode_mode,
                source_attestation=True),
            run_benchmark(
                snapshot_descriptors["control"], lock_descriptor,
                "control", attestation_cell, options.cpu, 1, 1,
                options.timeout, options.decode_mode,
                source_attestation=True),
            options.candidate_commit)
        results: list[dict[str, Any]] = []
        for cell_index, cell in enumerate(cells, 1):
            cell_result: dict[str, Any] = {"cell": cell, "comparisons": {}}
            for comparison in ("same_source", "exact_main"):
                attempts: list[dict[str, Any]] = []
                for _ in range(options.max_attempts):
                    validate_lock_identity(lock_descriptor, lock_identity)
                    attempt = run_attempt(
                        comparison, snapshot_descriptors, lock_descriptor, cell,
                        options.cpu,
                        options.sibling, options.rounds, options.iterations,
                        options.warmup, options.timeout, options.main_commit,
                        options.decode_mode)
                    attempts.append(attempt)
                    if attempt["accepted"]:
                        break
                require(attempts[-1]["accepted"],
                        f"sibling activity persisted for {cell['id']} "
                        f"{comparison}")
                cell_result["comparisons"][comparison] = {
                    "accepted_attempt": len(attempts) - 1,
                    "attempts": attempts,
                }
            results.append(cell_result)
            print(f"[{cell_index}/{len(cells)}] {cell['id']} accepted",
                  flush=True)
        identities_after = {
            name: descriptor_identity(
                original_descriptors[name], Path(identity["path"]))
            for name, identity in identities_before.items()
        }
        require(identities_after == identities_before,
                "a benchmark binary changed during the campaign")
        for name, descriptor in snapshot_descriptors.items():
            require(sha256_descriptor(descriptor) ==
                    snapshot_identities[name]["sha256"] and
                    fcntl.fcntl(descriptor, F_GET_SEALS) & REQUIRED_SEALS ==
                    REQUIRED_SEALS,
                    f"{name} sealed snapshot changed")
        build_identities_after = {
            "candidate": build_configuration_identity(
                options.candidate_build_dir, "ON"),
            "control": build_configuration_identity(
                options.control_build_dir, "OFF"),
        }
        require(build_identities_after == build_identities,
                "candidate/control build metadata changed")
        validate_lock_identity(lock_descriptor, lock_identity)
        summary: dict[str, Any] = {}
        comparison_metrics = ((
            ("same_source", ("one_shot",)),
            ("exact_main", ("one_shot",)),
        ) if options.decode_mode == "one_shot" else (
            ("same_source", ("execution", "amortized")),
            ("exact_main", ("amortized",)),
        ))
        for comparison, metric_names in comparison_metrics:
            summary[comparison] = {}
            for metric in metric_names:
                accepted = [
                    value["comparisons"][comparison]["attempts"]
                         [value["comparisons"][comparison]["accepted_attempt"]]
                         ["metrics"][metric]
                    for value in results
                ]
                speedups = [value["geometric_mean_speedup"]
                            for value in accepted]
                lower_bounds = [value["confidence_interval_lower"]
                                for value in accepted]
                summary[comparison][metric] = {
                    "cells": len(speedups),
                    "minimum_speedup": min(speedups),
                    "median_speedup": statistics.median(speedups),
                    "geometric_mean_speedup": math.exp(statistics.mean(
                        math.log(value) for value in speedups)),
                    "minimum_95_percent_lower_bound": min(lower_bounds),
                    "cells_with_credible_5_percent_gain": sum(
                        value >= 1.05 for value in lower_bounds),
                    "point_estimate_regressions": sum(
                        value < 1.0 for value in speedups),
                }
        report = {
            "schema": SCHEMA,
            "started_utc": started_utc,
            "finished_utc": time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "candidate_commit": options.candidate_commit,
            "main_commit": options.main_commit,
            "execution_nonce": execution_nonce,
            "comparison_semantics": {
                "same_source_execution":
                    "reused decode execution only",
                "same_source_amortized":
                    "decode execution plus plan setup divided by reuse",
                "exact_main_amortized":
                    "Leopard2 reused-plan amortized decode versus Leopard1 "
                    "decode including per-call setup",
                "irregular_byte_cells":
                    "Leopard2 exact logical bytes versus Leopard1 64-byte "
                    "zero-padded physical shards",
                "one_shot":
                    "public leo2_decode including transient plan setup versus "
                    "Leopard1 decode including per-call setup; reuse is only "
                    "the number of independent calls inside each timing sample",
            },
            "parameters": {
                "decode_mode": options.decode_mode,
                "rounds": options.rounds,
                "iterations": options.iterations,
                "warmup": options.warmup,
                "timeout_seconds": options.timeout,
                "max_attempts": options.max_attempts,
            },
            "machine": machine_record(options.cpu, options.sibling),
            "canonical_lock": lock_identity,
            "builds": {
                "before": build_identities,
                "after": build_identities_after,
            },
            "binaries": {
                "sources_before": identities_before,
                "sources_after": identities_after,
                "sealed_snapshots": snapshot_identities,
            },
            "source_attestations": source_attestations,
            "summary": summary,
            "results": results,
        }
        atomic_json(options.output, report)
        print(json.dumps(summary, sort_keys=True), flush=True)
    finally:
        for descriptor in (*snapshot_descriptors.values(),
                           *original_descriptors.values(), lock_descriptor):
            try:
                os.close(descriptor)
            except OSError:
                pass
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (EvidenceError, OSError, subprocess.SubprocessError) as error:
        print(f"small-dual promotion error: {error}", file=sys.stderr)
        raise SystemExit(1)
