#!/usr/bin/env python3
"""Reproduce Leopard2 tiny direct-encoder crossover measurements.

The runner consumes already-built ``bench_leopard2_direct_encode`` binaries.
It never configures or builds the project, which keeps authoritative pinned
measurements separate from compilation and other machine-wide work.  Every
cell is a resumable JSON job with deterministic input data and retained logs.

Typical use::

    tools/leopard2_direct_encode_crossover.py screen --build-root build
    tools/leopard2_direct_encode_crossover.py pinned --build-root build --cpu 16
    tools/leopard2_direct_encode_crossover.py analyze \
        --result-dir results/leopard2/direct-encode-crossover/pinned

The default build-root lookup accepts the repository's conventional trees
``build/direct-encode-SCALAR`` and ``build/SCALAR``.  ``--executable-root`` may
instead name a different root or contain a literal ``{backend}`` placeholder.
The forced-path benchmark is intentionally bounded to ``K,R <= 16``; the
adjacent 17 fallback is a correctness/dispatch test, not a meaningful
forced-direct crossover cell.
"""

from __future__ import print_function

import argparse
import concurrent.futures
import hashlib
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path


SCHEMA = "leopard2-direct-encode-crossover/v1"
JOB_SCHEMA = "leopard2-direct-encode-crossover-job/v1"
ANALYSIS_SCHEMA = "leopard2-direct-encode-crossover-analysis/v1"
BENCHMARK_SCHEMA = "leopard2-direct-encode-benchmark-v1"
KNOWN_BACKENDS = ("scalar", "ssse3", "avx2", "avx512")
SOURCE_FILES = (
    "CMakeLists.txt",
    "LeopardCommon.cpp",
    "LeopardCommon.h",
    "LeopardFF8.cpp",
    "LeopardFF8.h",
    "LeopardFF16.cpp",
    "LeopardFF16.h",
    "Leopard2Dispatch.h",
    "Leopard2Direct.h",
    "Leopard2Plan.cpp",
    "Leopard2Plan.h",
    "leopard.cpp",
    "leopard.h",
    "leopard2.cpp",
    "leopard2.h",
    "bench/leopard2/direct_encode_benchmark.cpp",
    "tools/leopard2_direct_encode_crossover.py",
)


class CrossoverError(Exception):
    """An actionable configuration, execution, or result error."""


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def digest_value(value):
    return digest_bytes(canonical_bytes(value))


def normalized_output(value):
    return value.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def atomic_write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True, ensure_ascii=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, str(path))
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def compact_cpu_list(cpus):
    values = sorted(set(int(cpu) for cpu in cpus))
    if not values:
        return ""
    result = []
    first = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        result.append(str(first) if first == previous else "{}-{}".format(first, previous))
        first = previous = value
    result.append(str(first) if first == previous else "{}-{}".format(first, previous))
    return ",".join(result)


def allowed_cpus():
    if hasattr(os, "sched_getaffinity"):
        try:
            result = sorted(os.sched_getaffinity(0))
            if result:
                return result
        except OSError:
            pass
    return list(range(os.cpu_count() or 1))


def read_optional(path):
    try:
        return Path(path).read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError):
        return None


def cpu_description(cpus):
    model = None
    flags = []
    cpuinfo = read_optional("/proc/cpuinfo")
    if cpuinfo:
        wanted = set(cpus)
        fallback = None
        for block in cpuinfo.split("\n\n"):
            fields = {}
            for line in block.splitlines():
                if ":" in line:
                    key, value = line.split(":", 1)
                    fields[key.strip().lower()] = value.strip()
            if fallback is None:
                fallback = fields
            try:
                processor = int(fields.get("processor", "-1"))
            except ValueError:
                processor = -1
            if processor in wanted:
                model = fields.get("model name", fields.get("hardware"))
                flags = fields.get("flags", fields.get("features", "")).split()
                break
        if model is None and fallback:
            model = fallback.get("model name", fallback.get("hardware"))
            flags = fallback.get("flags", fallback.get("features", "")).split()
    return model, sorted(set(flags))


def cpu_topology(cpu):
    root = Path("/sys/devices/system/cpu/cpu{}".format(cpu))
    result = {"cpu": cpu}
    fields = {
        "core_id": root / "topology/core_id",
        "package_id": root / "topology/physical_package_id",
        "thread_siblings": root / "topology/thread_siblings_list",
        "scaling_governor": root / "cpufreq/scaling_governor",
    }
    for key, path in fields.items():
        value = read_optional(path)
        if value is not None:
            result[key] = value
    return result


def machine_identity(cpus):
    model, flags = cpu_description(cpus)
    uname = platform.uname()
    result = {
        "allowed_cpu_list": compact_cpu_list(cpus),
        "architecture": platform.machine(),
        "cpu_flags": flags,
        "cpu_model": model,
        "logical_cpus_allowed": len(cpus),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "uname": {
            "machine": uname.machine,
            "node": uname.node,
            "release": uname.release,
            "system": uname.system,
            "version": uname.version,
        },
    }
    memory = read_optional("/proc/meminfo")
    if memory:
        for line in memory.splitlines():
            if line.startswith("MemTotal:"):
                result["memory_total"] = line.split(":", 1)[1].strip()
                break
    no_turbo = read_optional("/sys/devices/system/cpu/intel_pstate/no_turbo")
    boost = read_optional("/sys/devices/system/cpu/cpufreq/boost")
    if no_turbo is not None:
        result["intel_pstate_no_turbo"] = no_turbo
    if boost is not None:
        result["cpufreq_boost"] = boost
    return result


def source_fingerprint(source):
    files = {}
    missing = []
    for relative in SOURCE_FILES:
        path = source / relative
        if path.is_file():
            files[relative] = digest_bytes(path.read_bytes())
        else:
            missing.append(relative)
    if missing:
        raise CrossoverError(
            "required crossover input is missing: {}".format(", ".join(missing))
        )
    return {"digest": digest_value(files), "files": files}


def parse_backends(text):
    result = []
    for item in text.split(","):
        backend = item.strip().lower()
        if backend and backend not in result:
            result.append(backend)
    invalid = [item for item in result if item not in KNOWN_BACKENDS]
    if not result or invalid:
        raise CrossoverError(
            "backends must be a comma-separated subset of {}".format(
                ",".join(KNOWN_BACKENDS)
            )
        )
    return result


def executable_candidates(root, backend):
    rendered = Path(str(root).format(backend=backend))
    roots = [
        rendered,
        rendered / backend,
        rendered / ("direct-encode-" + backend),
    ]
    names = ("bench_leopard2_direct_encode", "bench_leopard2_direct_encode.exe")
    candidates = []
    for directory in roots:
        for name in names:
            candidates.extend((directory / name, directory / "Release" / name))
    unique = []
    seen = set()
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            seen.add(key)
            unique.append(candidate)
    return unique


def find_executable(root, backend):
    candidates = executable_candidates(root, backend)
    for candidate in candidates:
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            return candidate.resolve()
    raise CrossoverError(
        "benchmark executable for {} was not found; checked {}".format(
            backend, ", ".join(str(item) for item in candidates)
        )
    )


def cmake_build_metadata(executable):
    """Fingerprint the CMake inputs nearest an already-built benchmark."""
    cache = None
    directory = executable.parent
    for candidate_root in (directory, directory.parent, directory.parent.parent):
        candidate = candidate_root / "CMakeCache.txt"
        if candidate.is_file():
            cache = candidate
            break
    if cache is None:
        raise CrossoverError(
            "cannot find CMakeCache.txt near benchmark executable {}".format(
                executable
            )
        )
    try:
        cache_bytes = cache.read_bytes()
        cache_lines = cache_bytes.decode("utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise CrossoverError("cannot read {}: {}".format(cache, error))
    prefixes = (
        "CMAKE_BUILD_TYPE", "CMAKE_C_COMPILER", "CMAKE_C_FLAGS",
        "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS", "CMAKE_GENERATOR",
        "LEO2_BACKEND_VARIANT", "LEO2_BUILD_BENCHMARKS",
        "LEO2_BUILD_FUZZERS", "LEO2_BUILD_TESTS", "LEO2_ENABLE_CUDA",
    )
    entries = {}
    for line in cache_lines:
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        typed_key, value = line.split("=", 1)
        key = typed_key.split(":", 1)[0]
        if any(key == prefix or key.startswith(prefix + "_") for prefix in prefixes):
            entries[key] = value
    build_root = cache.parent
    extra_files = {}
    for relative in (
            "compile_commands.json", "CMakeFiles/CMakeConfigureLog.yaml",
            "CMakeFiles/CMakeOutput.log", "build.ninja",
            "CMakeFiles/leopard.dir/flags.make",
            "CMakeFiles/bench_leopard2_direct_encode.dir/flags.make"):
        path = build_root / relative
        if path.is_file():
            extra_files[relative] = digest_bytes(path.read_bytes())
    return {
        "cmake_cache": str(cache.resolve()),
        "cmake_cache_sha256": digest_bytes(cache_bytes),
        "entries": entries,
        "extra_file_sha256": extra_files,
    }


def cell(region, backend, k, r, profile, field, shard_bytes, q, mask):
    return {
        "backend": backend,
        "field": field,
        "k": k,
        "mask": mask,
        "profile": profile,
        "q": q,
        "r": r,
        "region": region,
        "shard_bytes": shard_bytes,
    }


def threshold_k(backend):
    return 3 if backend == "scalar" else 2


def compact_grid(backends, r):
    """Cells at each promoted condition and its deliberately excluded neighbors."""
    result = []
    for backend in backends:
        minimum_k = threshold_k(backend)
        # Threshold, large-shard behavior, and the bounded edge in both fields.
        for field in ("gf8", "gf16"):
            for k, shard_bytes in (
                    (minimum_k, 1024), (minimum_k, 65536), (16, 4096)):
                result.append(cell(
                    "candidate", backend, k, r, "low", field,
                    shard_bytes, 1, "0"))
        # One adjacent complete tile and an equal-Q sparse mask expose the two
        # AUTO inputs (byte shape and actual output count) independently.
        result.append(cell(
            "candidate", backend, minimum_k, r, "low", "gf8",
            1088, 1, "0"))
        if r > 1:
            result.append(cell(
                "candidate_sparse_output", backend, minimum_k, r,
                "low", "gf8", 4096, 1, str(r - 1)))
        # These neighbors are intentionally outside AUTO and guard each clause.
        if backend == "scalar":
            for field, shard_bytes in (
                    ("gf8", 1024), ("gf8", 65536), ("gf16", 4096)):
                result.append(cell(
                    "excluded_scalar_k2", backend, 2, r, "low", field,
                    shard_bytes, 1, "0"))
        if r >= 2:
            result.append(cell(
                "excluded_q2", backend, minimum_k, r, "low", "gf8",
                4096, 2, "0-1"))
        if r >= 3:
            result.append(cell(
                "excluded_q2_holey", backend, minimum_k, r, "low", "gf8",
                4096, 2, "0,{}".format(r - 1)))
        result.append(cell(
            "excluded_high_profile", backend, minimum_k, r, "high", "gf8",
            4096, 1, "0"))
        for shard_bytes in (63, 1023):
            result.append(cell(
                "excluded_ragged_tail", backend, minimum_k, r, "low", "gf8",
                shard_bytes, 1, "0"))
    return result


def full_grid(backends, r):
    result = compact_grid(backends, r)
    seen = {digest_value(item) for item in result}

    def add(item):
        identity = digest_value(item)
        if identity not in seen:
            seen.add(identity)
            result.append(item)

    for backend in backends:
        minimum_k = threshold_k(backend)
        candidate_ks = list(range(minimum_k, 17))
        for field in ("gf8", "gf16"):
            for k in candidate_ks:
                for shard_bytes in (1024, 1088, 2048, 4096, 16384, 65536, 1048576):
                    add(cell(
                        "candidate", backend, k, r, "low", field,
                        shard_bytes, 1, "0"))
            # A one-output sparse request exercises a non-prefix parity row.
            for k in (minimum_k, 4, 8, 16):
                if k < minimum_k:
                    continue
                add(cell(
                    "candidate_sparse_output", backend, k, r, "low", field,
                    4096, 1, str(r - 1)))
        for k in (minimum_k, 4, 8, 16):
            if k < minimum_k:
                continue
            for q in (2, 4, 8):
                if q > r:
                    continue
                add(cell(
                    "excluded_q_gt_1", backend, k, r, "low", "gf8",
                    4096, q, "0-{}".format(q - 1)))
            for profile in ("high",):
                add(cell(
                    "excluded_high_profile", backend, k, r, profile, "gf8",
                    4096, 1, "0"))
        for shard_bytes in (1, 63, 960, 1023, 1025, 1087, 1089, 4095, 4097):
            add(cell(
                "excluded_byte_boundary", backend, minimum_k, r, "low", "gf8",
                shard_bytes, 1, "0"))
        if backend == "scalar":
            for field in ("gf8", "gf16"):
                for shard_bytes in (1024, 4096, 65536, 1048576):
                    add(cell(
                        "excluded_scalar_k2", backend, 2, r, "low", field,
                        shard_bytes, 1, "0"))
        for k in range(1, minimum_k):
            for field in ("gf8", "gf16"):
                add(cell(
                    "excluded_k_below_auto", backend, k, r, "low", field,
                    4096, 1, "0"))
    return result


def parse_r_values(text):
    if text.strip().lower() == "all":
        return list(range(1, 17))
    result = []
    for item in text.split(","):
        try:
            value = int(item.strip())
        except ValueError:
            raise CrossoverError("R must be 'all' or a comma-separated integer list")
        if value < 1 or value > 16:
            raise CrossoverError("every R must be in [1,16]")
        if value not in result:
            result.append(value)
    if not result:
        raise CrossoverError("at least one R value is required")
    return sorted(result)


def sorted_grid(backends, r_values, full):
    values = []
    for r in r_values:
        values.extend(full_grid(backends, r) if full else compact_grid(backends, r))
    return sorted(values, key=lambda item: (
        item["backend"], item["region"], item["profile"], item["field"],
        item["r"], item["k"], item["q"], item["shard_bytes"], item["mask"]
    ))


def stable_seed(cell_value):
    # Keep zero reserved and remain within the benchmark's uint64 parser.
    return int(digest_value({"seed_for": cell_value})[:16], 16) or 1


def invocation_order(mode, job_id, abba_rounds):
    if mode == "pinned":
        return ("direct", "transform", "transform", "direct") * abba_rounds
    # Alternating the two-cell screening order prevents one path from always
    # receiving the colder process launch while keeping the order deterministic.
    return (("direct", "transform") if int(job_id[:8], 16) % 2 == 0
            else ("transform", "direct"))


def make_jobs(cells, executables, build_metadata, source_state, machine, settings):
    jobs = []
    for cell_value in cells:
        executable = executables[cell_value["backend"]]
        executable_sha256 = digest_bytes(executable.read_bytes())
        identity = {
            "cell": cell_value,
            "build_metadata": build_metadata[cell_value["backend"]],
            "executable": str(executable),
            "executable_sha256": executable_sha256,
            "machine": machine,
            "settings": settings,
            "source_fingerprint": source_state["digest"],
        }
        job_id = digest_value(identity)[:24]
        jobs.append({
            "cell": cell_value,
            "build_metadata": build_metadata[cell_value["backend"]],
            "configuration_id": digest_value(identity),
            "executable": str(executable),
            "executable_sha256": executable_sha256,
            "invocation_order": list(invocation_order(
                settings["mode"], job_id, settings["abba_rounds"])),
            "job_id": job_id,
            "seed": stable_seed(cell_value),
        })
    return jobs


def benchmark_argv(job, timed_mode, raw_path, settings):
    item = job["cell"]
    benchmark = settings["benchmark"]
    argv = [
        job["executable"],
        "--k", str(item["k"]),
        "--r", str(item["r"]),
        "--profile", item["profile"],
        "--field", item["field"],
        "--bytes", str(item["shard_bytes"]),
        "--q", str(item["q"]),
        "--requested-parity", item["mask"],
        "--batch", str(benchmark["batch"]),
        "--reuse", str(benchmark["reuse"]),
        "--iterations", str(benchmark["iterations"]),
        "--warmups", str(benchmark["warmups"]),
        "--threads", "1",
        "--seed", str(job["seed"]),
        "--mode", timed_mode,
        "--json", str(raw_path),
    ]
    if settings["mode"] == "pinned":
        argv = [settings["taskset"], "-c", str(settings["pin_cpu"])] + argv
    return argv


def run_command(argv, cwd, stdout_path, stderr_path, timeout, environment):
    timed_out = False
    try:
        completed = subprocess.run(
            [str(item) for item in argv], cwd=str(cwd), env=environment,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout
        )
        returncode = completed.returncode
        stdout = normalized_output(completed.stdout)
        stderr = normalized_output(completed.stderr)
    except subprocess.TimeoutExpired as error:
        timed_out = True
        returncode = 124
        stdout = normalized_output(error.stdout or b"")
        stderr = normalized_output(error.stderr or b"")
    stdout_path.write_bytes(stdout)
    stderr_path.write_bytes(stderr)
    return {
        "argv": [str(item) for item in argv],
        "cwd": str(cwd),
        "returncode": returncode,
        "stderr_log": stderr_path.name,
        "stderr_sha256": digest_bytes(stderr),
        "stdout_log": stdout_path.name,
        "stdout_sha256": digest_bytes(stdout),
        "timed_out": timed_out,
    }


def requested_indices(mask, r):
    indices = set()
    for part in mask.split(","):
        bounds = part.split("-", 1)
        try:
            first = int(bounds[0])
            last = int(bounds[1]) if len(bounds) == 2 else first
        except ValueError:
            raise CrossoverError("invalid generated parity mask {!r}".format(mask))
        if first < 0 or last < first or last >= r:
            raise CrossoverError("generated parity mask is out of range: {!r}".format(mask))
        indices.update(range(first, last + 1))
    return sorted(indices)


def validate_raw(raw, job, timed_mode, settings):
    if raw.get("schema") != BENCHMARK_SCHEMA:
        raise CrossoverError("benchmark emitted an unknown schema")
    parameters = raw.get("parameters", {})
    resolved = raw.get("resolved", {})
    correctness = raw.get("correctness", {})
    item = job["cell"]
    expected = {
        "K": item["k"], "R": item["r"], "Q": item["q"],
        "forced_mode": "force_" + timed_mode,
        "shard_bytes": item["shard_bytes"],
        "seed": job["seed"],
        "requested_parity_indices": requested_indices(item["mask"], item["r"]),
        "requested_profile": ("low_v1" if item["profile"] == "low"
                              else "legacy_high_v1"),
        "requested_field": item["field"],
        "batch": settings["benchmark"]["batch"],
        "reuse": settings["benchmark"]["reuse"],
        "iterations": settings["benchmark"]["iterations"],
        "warmups": settings["benchmark"]["warmups"],
        "thread_count": 1,
    }
    for key, value in expected.items():
        if parameters.get(key) != value:
            raise CrossoverError(
                "benchmark parameter {} is {!r}, expected {!r}".format(
                    key, parameters.get(key), value
                )
            )
    if resolved.get("profile") != ("low_v1" if item["profile"] == "low"
                                    else "legacy_high_v1"):
        raise CrossoverError("benchmark resolved an unexpected profile")
    if resolved.get("field") != item["field"]:
        raise CrossoverError("benchmark resolved an unexpected field")
    if item["backend"] != "auto" and resolved.get("backend") != item["backend"]:
        raise CrossoverError(
            "benchmark resolved backend {!r}, expected {!r}".format(
                resolved.get("backend"), item["backend"]
            )
        )
    if correctness.get("direct_transform_parity_match") is not True:
        raise CrossoverError("benchmark direct/transform correctness check failed")
    if correctness.get("unrequested_outputs_untouched") is not True:
        raise CrossoverError("benchmark unrequested-output guard check failed")
    checksum = correctness.get("parity_checksum_fnv1a64")
    if not isinstance(checksum, str) or not re.fullmatch(r"0x[0-9a-fA-F]{16}", checksum):
        raise CrossoverError("benchmark omitted a well-formed parity checksum")
    if resolved.get("direct_capable") is not True or resolved.get("thread_count") != 1:
        raise CrossoverError("benchmark resolved unexpected capability or thread count")
    expected_direct = timed_mode == "direct"
    if resolved.get("timed_path_is_direct") is not expected_direct:
        raise CrossoverError("benchmark timed a path other than the forced path")
    try:
        median_us = float(raw["metrics"]["encode_execution"][
            "median_us_per_batch_call"])
        mad_us = float(raw["metrics"]["encode_execution"][
            "mad_us_per_batch_call"])
    except (KeyError, TypeError, ValueError):
        raise CrossoverError("benchmark omitted numeric execution metrics")
    if median_us <= 0 or mad_us < 0:
        raise CrossoverError("benchmark emitted invalid execution metrics")
    return median_us, mad_us, checksum.lower()


def validate_execution_inputs(job):
    executable = Path(job["executable"])
    try:
        executable_sha256 = digest_bytes(executable.read_bytes())
    except OSError as error:
        raise CrossoverError("cannot re-read executable {}: {}".format(
            executable, error
        ))
    if executable_sha256 != job["executable_sha256"]:
        raise CrossoverError("benchmark executable changed during the run")
    if canonical_bytes(cmake_build_metadata(executable)) != canonical_bytes(
            job["build_metadata"]):
        raise CrossoverError("benchmark CMake metadata changed during the run")


def artifact_path(root, relative, description):
    if not isinstance(relative, str) or not relative:
        raise CrossoverError("{} path is missing".format(description))
    root = root.resolve()
    candidate = (root / relative).resolve()
    try:
        common = os.path.commonpath((str(root), str(candidate)))
    except ValueError:
        common = ""
    if common != str(root):
        raise CrossoverError("{} escapes the result directory".format(description))
    return candidate


def validate_job_artifacts(result_dir, result, expected_job, settings):
    if result.get("status") != "passed":
        return
    for key in (
            "build_metadata", "cell", "configuration_id", "executable",
            "executable_sha256", "job_id", "seed"):
        if canonical_bytes(result.get(key)) != canonical_bytes(expected_job.get(key)):
            raise CrossoverError(
                "passed job {} has stale {} metadata".format(
                    expected_job["job_id"], key
                )
            )
    order = expected_job.get("invocation_order", [])
    commands = result.get("commands")
    measurements = result.get("measurements")
    if (not order or not isinstance(commands, list) or
            not isinstance(measurements, list) or
            len(commands) != len(order) or len(measurements) != len(order)):
        raise CrossoverError(
            "passed job {} has an incomplete command/measurement sequence".format(
                expected_job["job_id"]
            )
        )
    checksums = set()
    log_root = result_dir / "logs" / expected_job["job_id"]
    for index, timed_mode in enumerate(order):
        command = commands[index]
        measurement = measurements[index]
        if (command.get("returncode") != 0 or command.get("timed_out") is not False or
                measurement.get("sequence_index") != index or
                measurement.get("timed_mode") != timed_mode):
            raise CrossoverError(
                "passed job {} has an invalid sequence entry {}".format(
                    expected_job["job_id"], index
                )
            )
        for stream in ("stdout", "stderr"):
            name = command.get(stream + "_log")
            if not isinstance(name, str) or Path(name).name != name:
                raise CrossoverError("job log name is missing or unsafe")
            path = artifact_path(log_root, name, stream + " log")
            try:
                value = path.read_bytes()
            except OSError as error:
                raise CrossoverError("cannot read {}: {}".format(path, error))
            if digest_bytes(value) != command.get(stream + "_sha256"):
                raise CrossoverError("{} hash does not match the job record".format(path))
        raw_path = artifact_path(
            result_dir, measurement.get("benchmark_json"), "benchmark JSON"
        )
        try:
            raw_bytes = raw_path.read_bytes()
            raw = json.loads(raw_bytes.decode("utf-8"))
        except (OSError, UnicodeError, ValueError) as error:
            raise CrossoverError("cannot read {}: {}".format(raw_path, error))
        if digest_bytes(raw_bytes) != measurement.get("benchmark_json_sha256"):
            raise CrossoverError("{} hash does not match the job record".format(raw_path))
        median_us, mad_us, checksum = validate_raw(
            raw, expected_job, timed_mode, settings
        )
        try:
            recorded_median = float(measurement.get("median_us", -1))
            recorded_mad = float(measurement.get("mad_us", -1))
        except (TypeError, ValueError):
            raise CrossoverError("passed job contains non-numeric recorded metrics")
        if recorded_median != median_us or recorded_mad != mad_us:
            raise CrossoverError("raw metrics differ from the passed job summary")
        checksums.add(checksum)
    if len(checksums) != 1:
        raise CrossoverError("passed job raw artifacts contain different parity checksums")
    recomputed = summarize_measurements(measurements)
    if canonical_bytes(recomputed) != canonical_bytes(result.get("summary")):
        raise CrossoverError("passed job aggregate does not match its raw measurements")


def median(values):
    return float(statistics.median(values))


def summarize_measurements(measurements):
    by_mode = {"direct": [], "transform": []}
    for measurement in measurements:
        by_mode[measurement["timed_mode"]].append(measurement["median_us"])
    if not by_mode["direct"] or not by_mode["transform"]:
        raise CrossoverError("job did not measure both encoder paths")
    direct = median(by_mode["direct"])
    transform = median(by_mode["transform"])
    gain = (transform / direct - 1.0) * 100.0
    rounds = []
    if len(measurements) % 4 == 0 and len(measurements) >= 4:
        for offset in range(0, len(measurements), 4):
            group = measurements[offset:offset + 4]
            if [item["timed_mode"] for item in group] != [
                    "direct", "transform", "transform", "direct"]:
                rounds = []
                break
            round_direct = median([group[0]["median_us"], group[3]["median_us"]])
            round_transform = median([group[1]["median_us"], group[2]["median_us"]])
            rounds.append({
                "direct_median_us": round_direct,
                "gain_percent": (round_transform / round_direct - 1.0) * 100.0,
                "transform_median_us": round_transform,
            })
    return {
        "direct_invocation_medians_us": by_mode["direct"],
        "direct_median_us": direct,
        "gain_percent": gain,
        "rounds": rounds,
        "transform_invocation_medians_us": by_mode["transform"],
        "transform_median_us": transform,
    }


def run_job(job, context):
    result_path = context["result_dir"] / "jobs" / (job["job_id"] + ".json")
    if context["resume"] and result_path.is_file():
        try:
            previous = json.loads(result_path.read_text(encoding="utf-8"))
            if (previous.get("configuration_id") == job["configuration_id"] and
                    previous.get("status") == "passed"):
                validate_job_artifacts(
                    context["result_dir"], previous, job, context["settings"]
                )
                return previous
        except (CrossoverError, OSError, UnicodeError, ValueError):
            pass

    log_dir = context["result_dir"] / "logs" / job["job_id"]
    raw_dir = context["result_dir"] / "raw" / job["job_id"]
    log_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "build_metadata": job["build_metadata"],
        "cell": job["cell"],
        "commands": [],
        "configuration_id": job["configuration_id"],
        "executable": job["executable"],
        "executable_sha256": job["executable_sha256"],
        "job_id": job["job_id"],
        "measurements": [],
        "reason": "",
        "resumed": False,
        "schema": JOB_SCHEMA,
        "seed": job["seed"],
        "status": "failed",
    }
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["LEO2_EXPECT_BACKEND"] = job["cell"]["backend"]
    checksums = set()
    try:
        for index, timed_mode in enumerate(job["invocation_order"]):
            validate_execution_inputs(job)
            label = "{:02d}-{}".format(index, timed_mode)
            raw_path = raw_dir / (label + ".json")
            stdout_path = log_dir / (label + ".stdout.log")
            stderr_path = log_dir / (label + ".stderr.log")
            argv = benchmark_argv(job, timed_mode, raw_path, context["settings"])
            command = run_command(
                argv, context["source"], stdout_path, stderr_path,
                context["timeout"], environment
            )
            command["environment_overrides"] = {
                "LEO2_EXPECT_BACKEND": job["cell"]["backend"],
                "OMP_NUM_THREADS": "1",
            }
            result["commands"].append(command)
            validate_execution_inputs(job)
            if command["returncode"] != 0:
                raise CrossoverError(
                    "{} exited with status {}".format(label, command["returncode"])
                )
            try:
                raw_bytes = raw_path.read_bytes()
                raw = json.loads(raw_bytes.decode("utf-8"))
            except (OSError, UnicodeError, ValueError) as error:
                raise CrossoverError("cannot read {}: {}".format(raw_path, error))
            median_us, mad_us, checksum = validate_raw(
                raw, job, timed_mode, context["settings"]
            )
            checksums.add(checksum)
            result["measurements"].append({
                "benchmark_json": str(raw_path.relative_to(context["result_dir"])),
                "benchmark_json_sha256": digest_bytes(raw_bytes),
                "mad_us": mad_us,
                "median_us": median_us,
                "sequence_index": index,
                "timed_mode": timed_mode,
            })
        if len(checksums) != 1:
            raise CrossoverError("forced paths emitted different parity checksums")
        result["summary"] = summarize_measurements(result["measurements"])
        result["status"] = "passed"
    except (CrossoverError, OSError, ValueError) as error:
        result["reason"] = str(error)
    atomic_write_json(result_path, result)
    return result


def summarize_region(results, promotion_percent):
    passed = [item for item in results if item.get("status") == "passed"]
    gains = sorted(float(item["summary"]["gain_percent"]) for item in passed)
    round_gains = sorted(
        float(round_value["gain_percent"])
        for item in passed for round_value in item["summary"].get("rounds", [])
    )
    return {
        "cell_count": len(results),
        "failed_count": len(results) - len(passed),
        "gain_max_percent": max(gains) if gains else None,
        "gain_median_percent": median(gains) if gains else None,
        "gain_min_percent": min(gains) if gains else None,
        "improvement_count": sum(value > 0 for value in gains),
        "promotion_count": sum(value >= promotion_percent for value in gains),
        "regression_count": sum(value < 0 for value in gains),
        "round_gain_median_percent": median(round_gains) if round_gains else None,
        "round_gain_min_percent": min(round_gains) if round_gains else None,
        "round_regression_count": sum(value < 0 for value in round_gains),
        "severe_regression_count": sum(value < -2.0 for value in gains),
    }


def analyze_results(
        results, promotion_percent, manifest_configuration_id=None,
        run_status=None, source_changed_during_run=None,
        execution_input_error=None):
    ordered = sorted(results, key=lambda item: item.get("job_id", ""))
    regions = {}
    for item in ordered:
        region = item.get("cell", {}).get("region", "unknown")
        regions.setdefault(region, []).append(item)
    candidate = []
    excluded = []
    for region, values in regions.items():
        if region.startswith("candidate"):
            candidate.extend(values)
        else:
            excluded.extend(values)
    candidate_summary = summarize_region(candidate, promotion_percent)
    analysis = {
        "candidate": candidate_summary,
        "excluded_neighbors": summarize_region(excluded, promotion_percent),
        "execution_input_error": execution_input_error,
        "jobs_failed": sum(item.get("status") != "passed" for item in ordered),
        "jobs_passed": sum(item.get("status") == "passed" for item in ordered),
        "jobs_total": len(ordered),
        "manifest_configuration_id": manifest_configuration_id,
        "promotion_percent": promotion_percent,
        "regions": {
            region: summarize_region(values, promotion_percent)
            for region, values in sorted(regions.items())
        },
        "run_status": run_status,
        "schema": ANALYSIS_SCHEMA,
        "source_changed_during_run": source_changed_during_run,
    }
    analysis["candidate"]["all_cells_meet_promotion_threshold"] = bool(candidate) and (
        candidate_summary["failed_count"] == 0 and
        candidate_summary["promotion_count"] == candidate_summary["cell_count"]
    )
    return analysis


def load_manifest(result_dir):
    path = result_dir / "manifest.json"
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CrossoverError("cannot read {}: {}".format(path, error))
    if manifest.get("schema") != SCHEMA or not manifest.get("configuration_id"):
        raise CrossoverError("{} has an unknown or incomplete schema".format(path))
    if not isinstance(manifest.get("jobs"), list):
        raise CrossoverError("{} does not contain an expected job list".format(path))
    return manifest


def require_compatible_result_dir(result_dir, manifest):
    path = result_dir / "manifest.json"
    if not path.exists():
        job_dir = result_dir / "jobs"
        if job_dir.is_dir() and any(job_dir.glob("*.json")):
            raise CrossoverError(
                "result directory has jobs but no manifest: {}".format(result_dir)
            )
        return
    previous = load_manifest(result_dir)
    if previous.get("configuration_id") != manifest.get("configuration_id"):
        raise CrossoverError(
            "result directory belongs to configuration {}; new configuration is {}; "
            "select a new --result-dir rather than mixing jobs".format(
                previous.get("configuration_id"), manifest.get("configuration_id")
            )
        )


def load_job_results(result_dir, manifest):
    job_dir = result_dir / "jobs"
    if not job_dir.is_dir():
        raise CrossoverError("job directory does not exist: {}".format(job_dir))
    expected = {}
    for job in manifest["jobs"]:
        job_id = job.get("job_id")
        configuration_id = job.get("configuration_id")
        if not job_id or not configuration_id or job_id in expected:
            raise CrossoverError("manifest contains a duplicate or incomplete job")
        expected[job_id] = job
    actual_paths = {path.stem: path for path in sorted(job_dir.glob("*.json"))}
    missing = sorted(set(expected) - set(actual_paths))
    extra = sorted(set(actual_paths) - set(expected))
    if missing or extra:
        raise CrossoverError(
            "job set does not match manifest (missing: {}; stale/extra: {})".format(
                ",".join(missing) or "none", ",".join(extra) or "none"
            )
        )
    results = []
    for job_id in sorted(expected):
        path = actual_paths[job_id]
        try:
            item = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, ValueError) as error:
            raise CrossoverError("cannot read {}: {}".format(path, error))
        if item.get("schema") != JOB_SCHEMA:
            raise CrossoverError("{} has an unknown job schema".format(path))
        if item.get("job_id") != job_id:
            raise CrossoverError("{} contains the wrong job ID".format(path))
        if item.get("configuration_id") != expected[job_id]["configuration_id"]:
            raise CrossoverError("{} belongs to a stale configuration".format(path))
        validate_job_artifacts(
            result_dir, item, expected[job_id], manifest["settings"]
        )
        results.append(item)
    return results


def write_merged(
        result_dir, manifest, results, source_end, promotion_percent,
        execution_input_error=""):
    ordered = sorted(results, key=lambda item: item["job_id"])
    source_changed = source_end is not None and (
        source_end["digest"] != manifest["source_fingerprint"]["digest"]
    )
    failed = any(item.get("status") != "passed" for item in ordered)
    run_status = "failed" if failed or source_changed or execution_input_error else "passed"
    analysis = analyze_results(
        ordered, promotion_percent, manifest["configuration_id"],
        run_status, source_changed, execution_input_error or None
    )
    merged = {
        "analysis": analysis,
        "execution_input_error": execution_input_error or None,
        "jobs": ordered,
        "manifest_configuration_id": manifest["configuration_id"],
        "schema": SCHEMA,
        "source_changed_during_run": source_changed,
        "source_fingerprint": manifest["source_fingerprint"],
        "source_fingerprint_after": source_end,
        "status": run_status,
    }
    atomic_write_json(result_dir / "matrix.json", merged)
    atomic_write_json(result_dir / "analysis.json", analysis)
    return merged


def print_analysis(analysis):
    candidate = analysis["candidate"]
    print(
        "candidate cells: {}/{} passed; gain min={}; median={}; "
        "regressions={}; >= {:.1f}%={}".format(
            candidate["cell_count"] - candidate["failed_count"],
            candidate["cell_count"],
            ("n/a" if candidate["gain_min_percent"] is None else
             "{:.2f}%".format(candidate["gain_min_percent"])),
            ("n/a" if candidate["gain_median_percent"] is None else
             "{:.2f}%".format(candidate["gain_median_percent"])),
            candidate["regression_count"], analysis["promotion_percent"],
            candidate["promotion_count"],
        )
    )
    neighbors = analysis["excluded_neighbors"]
    print(
        "excluded neighbors: {}/{} passed; gain min={}; median={}; regressions={}".format(
            neighbors["cell_count"] - neighbors["failed_count"],
            neighbors["cell_count"],
            ("n/a" if neighbors["gain_min_percent"] is None else
             "{:.2f}%".format(neighbors["gain_min_percent"])),
            ("n/a" if neighbors["gain_median_percent"] is None else
             "{:.2f}%".format(neighbors["gain_median_percent"])),
            neighbors["regression_count"],
        )
    )


def resolve_path(source, value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (source / path).resolve()


def run_matrix(arguments):
    source = Path(arguments.source).resolve()
    executable_root = resolve_path(
        source, arguments.executable_root or arguments.build_root
    )
    result_dir = resolve_path(source, arguments.result_dir)
    backends = parse_backends(arguments.backends)
    cpus = allowed_cpus()
    pin_cpu = None
    taskset = None
    if arguments.command == "pinned":
        pin_cpu = cpus[0] if arguments.cpu is None else arguments.cpu
        if pin_cpu not in cpus:
            raise CrossoverError(
                "pinned CPU {} is outside allowed affinity {}".format(
                    pin_cpu, compact_cpu_list(cpus)
                )
            )
        taskset = shutil.which(arguments.taskset)
        if not taskset:
            raise CrossoverError("pinned mode requires the taskset executable")
        if arguments.workers != 1:
            raise CrossoverError("pinned ABBA measurements require --workers 1")
    source_state = source_fingerprint(source)
    machine = machine_identity(cpus)
    if pin_cpu is not None:
        machine["pinned_cpu_topology"] = cpu_topology(pin_cpu)
    executables = {
        backend: find_executable(executable_root, backend) for backend in backends
    }
    build_metadata = {
        backend: cmake_build_metadata(executable)
        for backend, executable in executables.items()
    }
    settings = {
        "abba_rounds": arguments.abba_rounds if arguments.command == "pinned" else 0,
        "benchmark": {
            "batch": arguments.batch,
            "iterations": arguments.iterations,
            "reuse": arguments.reuse,
            "warmups": arguments.warmups,
        },
        "mode": arguments.command,
        "pin_cpu": pin_cpu,
        "placement_policy": (
            "external taskset single-CPU pin"
            if arguments.command == "pinned" else
            "unpinned discovery jobs inherit allowed affinity; OS schedules workers"
        ),
        "taskset": str(Path(taskset).resolve()) if taskset else None,
        "timeout_seconds_per_invocation": arguments.timeout,
        "workers": arguments.workers,
    }
    if arguments.command == "pinned":
        settings["isolation_prerequisite"] = (
            "taskset pins the benchmark process only; the caller must keep its "
            "SMT sibling and competing system work idle"
        )
    r_values = parse_r_values(arguments.r)
    cells = sorted_grid(backends, r_values, arguments.full)
    jobs = make_jobs(
        cells, executables, build_metadata, source_state, machine, settings
    )
    manifest_identity = {
        "jobs": jobs,
        "machine": machine,
        "settings": settings,
        "source_fingerprint": source_state,
    }
    manifest = {
        "configuration_id": digest_value(manifest_identity),
        "executables": {
            backend: {
                "build_metadata": build_metadata[backend],
                "path": str(path), "sha256": digest_bytes(path.read_bytes())
            } for backend, path in sorted(executables.items())
        },
        "jobs": jobs,
        "machine": machine,
        "schema": SCHEMA,
        "settings": settings,
        "source_fingerprint": source_state,
    }
    require_compatible_result_dir(result_dir, manifest)
    atomic_write_json(result_dir / "manifest.json", manifest)
    print(
        "direct-encode crossover: mode={} cells={} backends={} cpus={}{}".format(
            arguments.command, len(jobs), ",".join(backends),
            machine["allowed_cpu_list"],
            " pinned={}".format(pin_cpu) if pin_cpu is not None else "",
        ), flush=True
    )
    context = {
        "result_dir": result_dir,
        "resume": not arguments.no_resume,
        "settings": settings,
        "source": source,
        "timeout": arguments.timeout,
    }
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        futures = {executor.submit(run_job, job, context): job for job in jobs}
        for future in concurrent.futures.as_completed(futures):
            job = futures[future]
            result = future.result()
            results.append(result)
            print(
                "{} {}: {}".format(
                    job["cell"]["backend"], job["job_id"], result["status"]
                ), flush=True
            )
    source_end = source_fingerprint(source)
    execution_input_error = ""
    try:
        checked = set()
        for job in jobs:
            key = (job["executable"], job["executable_sha256"])
            if key not in checked:
                checked.add(key)
                validate_execution_inputs(job)
    except CrossoverError as error:
        execution_input_error = str(error)
    merged = write_merged(
        result_dir, manifest, results, source_end, arguments.promotion_percent,
        execution_input_error
    )
    print_analysis(merged["analysis"])
    print("matrix: {} ({})".format(merged["status"], result_dir / "matrix.json"))
    return 0 if merged["status"] == "passed" else 1


def analyze_command(arguments):
    result_dir = Path(arguments.result_dir).resolve()
    manifest = load_manifest(result_dir)
    results = load_job_results(result_dir, manifest)
    matrix_path = result_dir / "matrix.json"
    try:
        matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CrossoverError("cannot read completed matrix {}: {}".format(
            matrix_path, error
        ))
    if (matrix.get("schema") != SCHEMA or
            matrix.get("manifest_configuration_id") != manifest["configuration_id"]):
        raise CrossoverError("matrix does not match the current manifest")
    if (matrix.get("source_fingerprint", {}).get("digest") !=
            manifest.get("source_fingerprint", {}).get("digest")):
        raise CrossoverError("matrix source fingerprint does not match the manifest")
    matrix_jobs = sorted(matrix.get("jobs", []), key=lambda item: item.get("job_id", ""))
    ordered_results = sorted(results, key=lambda item: item.get("job_id", ""))
    if canonical_bytes(matrix_jobs) != canonical_bytes(ordered_results):
        raise CrossoverError("matrix job snapshot differs from the validated job files")
    source_changed = matrix.get("source_changed_during_run")
    run_status = matrix.get("status")
    if not isinstance(source_changed, bool) or run_status not in ("passed", "failed"):
        raise CrossoverError("matrix omits source-change or status evidence")
    analysis = analyze_results(
        results, arguments.promotion_percent, manifest["configuration_id"],
        run_status, source_changed, matrix.get("execution_input_error")
    )
    output = Path(arguments.output).resolve() if arguments.output else result_dir / "analysis.json"
    atomic_write_json(output, analysis)
    print_analysis(analysis)
    print("analysis: {}".format(output))
    return 0 if analysis["jobs_failed"] == 0 and run_status == "passed" else 1


def self_test():
    assert compact_cpu_list([3, 2, 1, 7, 9, 8]) == "1-3,7-9"
    assert digest_value({"b": 2, "a": 1}) == digest_value({"a": 1, "b": 2})
    assert stable_seed({"cell": 1}) == stable_seed({"cell": 1})
    assert parse_backends("avx512") == ["avx512"]
    assert parser().parse_args(["screen"]).backends == "scalar,ssse3,avx2"
    grid = sorted_grid(["scalar", "ssse3", "avx2"], [16], False)
    candidates = [item for item in grid if item["region"] == "candidate"]
    assert candidates
    for item in candidates:
        assert item["profile"] == "low" and item["q"] == 1
        assert item["shard_bytes"] >= 1024 and item["shard_bytes"] % 64 == 0
        assert item["k"] >= threshold_k(item["backend"])
    regions = {item["region"] for item in grid}
    assert {
        "excluded_scalar_k2", "excluded_q2", "excluded_high_profile",
        "excluded_ragged_tail",
    }.issubset(regions)
    assert {item["field"] for item in candidates} == {"gf8", "gf16"}
    assert any(item["region"] == "candidate_sparse_output" for item in grid)
    assert any(item["region"] == "excluded_q2_holey" for item in grid)
    r1 = sorted_grid(["scalar"], [1], False)
    assert all(item["r"] == 1 and item["q"] == 1 for item in r1)
    assert parse_r_values("1,16,1") == [1, 16]
    assert parse_r_values("all") == list(range(1, 17))
    assert requested_indices("0,2-3", 4) == [0, 2, 3]
    assert invocation_order("pinned", "00000000", 2) == (
        "direct", "transform", "transform", "direct",
        "direct", "transform", "transform", "direct",
    )
    measurements = [
        {"timed_mode": "direct", "median_us": 8.0},
        {"timed_mode": "transform", "median_us": 10.0},
        {"timed_mode": "transform", "median_us": 10.0},
        {"timed_mode": "direct", "median_us": 8.0},
    ]
    summary = summarize_measurements(measurements)
    assert summary["gain_percent"] == 25.0
    assert summary["rounds"][0]["gain_percent"] == 25.0
    fake = {
        "cell": candidates[0], "job_id": "a", "status": "passed",
        "summary": summary,
    }
    analysis = analyze_results([fake], 5.0)
    assert analysis["candidate"]["gain_min_percent"] == 25.0
    assert analysis["candidate"]["regression_count"] == 0
    with tempfile.TemporaryDirectory(prefix="leo2-crossover-self-test-") as directory:
        root = Path(directory)
        path = root / "stable.json"
        atomic_write_json(path, {"z": [3, 2, 1], "a": "stable"})
        assert json.loads(path.read_text(encoding="utf-8"))["a"] == "stable"
        record = run_command(
            [sys.executable, "-c", "print('ok')"], root,
            root / "stdout.log", root / "stderr.log", 5, os.environ.copy()
        )
        assert record["returncode"] == 0
        assert (root / "stdout.log").read_bytes() == b"ok\n"
        result_root = root / "results"
        manifest = {
            "configuration_id": "manifest-a",
            "jobs": [
                {"job_id": "a", "configuration_id": "cell-a"},
                {"job_id": "b", "configuration_id": "cell-b"},
            ],
            "schema": SCHEMA,
            "settings": {},
        }
        atomic_write_json(result_root / "manifest.json", manifest)
        require_compatible_result_dir(result_root, manifest)
        incompatible = dict(manifest)
        incompatible["configuration_id"] = "manifest-b"
        try:
            require_compatible_result_dir(result_root, incompatible)
            raise AssertionError("stale manifest was accepted")
        except CrossoverError:
            pass
        atomic_write_json(result_root / "jobs" / "a.json", {
            "configuration_id": "cell-a", "job_id": "a", "schema": JOB_SCHEMA,
        })
        try:
            load_job_results(result_root, manifest)
            raise AssertionError("partial job set was accepted")
        except CrossoverError:
            pass
        atomic_write_json(result_root / "jobs" / "b.json", {
            "configuration_id": "cell-b", "job_id": "b", "schema": JOB_SCHEMA,
        })
        atomic_write_json(result_root / "jobs" / "stale.json", {
            "configuration_id": "old", "job_id": "stale", "schema": JOB_SCHEMA,
        })
        try:
            load_job_results(result_root, manifest)
            raise AssertionError("stale extra job was accepted")
        except CrossoverError:
            pass
        (result_root / "jobs" / "stale.json").unlink()
        assert len(load_job_results(result_root, manifest)) == 2
    print("leopard2 direct-encode crossover self-test passed (no codec required)")
    return 0


def add_run_arguments(parser, default_result):
    default_source = str(Path(__file__).resolve().parents[1])
    parser.add_argument("--source", default=default_source)
    parser.add_argument("--build-root", default="build")
    parser.add_argument(
        "--executable-root", default=None,
        help="executable root; may contain a literal {backend} placeholder"
    )
    parser.add_argument("--result-dir", default=default_result)
    parser.add_argument("--backends", default="scalar,ssse3,avx2")
    parser.add_argument(
        "--r", default="16",
        help="R value, comma-separated R values, or 'all' (default 16)"
    )
    parser.add_argument("--batch", type=int, default=1)
    parser.add_argument("--reuse", type=int, default=32)
    parser.add_argument("--iterations", type=int, default=7)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--promotion-percent", type=float, default=5.0)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--no-resume", action="store_true")


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command")
    screen = subparsers.add_parser(
        "screen", help="run the compact grid without authoritative pinning"
    )
    add_run_arguments(
        screen, "results/leopard2/direct-encode-crossover/screen"
    )
    screen.set_defaults(workers=min(len(allowed_cpus()), 128))
    pinned = subparsers.add_parser(
        "pinned", help="run isolated, externally pinned ABBA measurements"
    )
    add_run_arguments(
        pinned, "results/leopard2/direct-encode-crossover/pinned"
    )
    pinned.set_defaults(iterations=15, warmups=4, reuse=64, workers=1)
    pinned.add_argument("--cpu", type=int, default=None)
    pinned.add_argument("--taskset", default="taskset")
    pinned.add_argument("--abba-rounds", type=int, default=3)
    analyze = subparsers.add_parser(
        "analyze", help="deterministically reanalyze completed job JSON"
    )
    analyze.add_argument(
        "--result-dir",
        default="results/leopard2/direct-encode-crossover/pinned",
    )
    analyze.add_argument("--promotion-percent", type=float, default=5.0)
    analyze.add_argument("--output", default=None)
    subparsers.add_parser("self-test", help="test the runner without a built codec")
    return result


def main():
    arguments = parser().parse_args()
    try:
        if arguments.command == "self-test":
            return self_test()
        if arguments.command == "analyze":
            return analyze_command(arguments)
        if arguments.command in ("screen", "pinned"):
            numeric = (
                arguments.batch, arguments.reuse, arguments.iterations,
                arguments.warmups, arguments.timeout, arguments.workers,
            )
            if any(value <= 0 for value in numeric):
                raise CrossoverError(
                    "batch, reuse, iterations, warmups, timeout, and workers "
                    "must be positive"
                )
            if arguments.command == "pinned" and arguments.abba_rounds <= 0:
                raise CrossoverError("abba-rounds must be positive")
            return run_matrix(arguments)
        parser().print_help()
        return 2
    except (CrossoverError, OSError, subprocess.SubprocessError) as error:
        print("direct-encode crossover error: {}".format(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
