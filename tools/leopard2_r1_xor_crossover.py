#!/usr/bin/env python3
"""Resumable, paired Leopard2 R=1 XOR crossover measurements.

The runner compares already-built binaries from an exact baseline checkout and
an exact candidate checkout.  It never configures or builds either tree.
Pinned runs execute one job at a time on one allowed logical CPU in deterministic
baseline/candidate/candidate/baseline (ABBA) order.  Each invocation retains the
benchmark JSON and stdout/stderr, and completed jobs resume only after their
identity and artifact hashes have been revalidated.

Typical use::

    tools/leopard2_r1_xor_crossover.py run \
      --baseline ../baseline/build/release/bench_leopard2 \
      --candidate build/release/bench_leopard2 \
      --cpu 15 --result-dir .research/leopard2/r1-xor/pinned
    tools/leopard2_r1_xor_crossover.py analyze \
      --result-dir .research/leopard2/r1-xor/pinned

Compilation and unrelated workloads must be stopped during a pinned run.  The
manifest records CPU topology plus SHA-256 identities for the executable and
CMake cache files.  The recorded source commit and dirty flag are observations
of the source path named by CMake at manifest creation time: the benchmark does
not embed a commit identifier, so those fields do not cryptographically bind the
executable to a source tree.  The operator also remains responsible for
reserving the physical core's SMT sibling.
"""

from __future__ import print_function

import argparse
import hashlib
import json
import os
import platform
import random
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path


SCHEMA = "leopard2-r1-xor-crossover/v1"
JOB_SCHEMA = "leopard2-r1-xor-crossover-job/v1"
ANALYSIS_SCHEMA = "leopard2-r1-xor-crossover-analysis/v1"
BENCHMARK_SCHEMA = "leopard2-benchmark-v2"
BACKENDS = ("scalar", "ssse3", "avx2")
ACCEPTED_BACKENDS = BACKENDS + ("avx512",)


class CrossoverError(Exception):
    pass


class SelfTestChecks:
    def __init__(self):
        self.count = 0

    def require(self, condition, description):
        self.count += 1
        if not condition:
            raise CrossoverError(
                "R=1 crossover self-test failed: {}".format(description))


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True).encode("utf-8")


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def digest_value(value):
    return digest_bytes(canonical(value))


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True,
                      ensure_ascii=True)
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


def normalized(value):
    return value.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def allowed_cpus():
    if hasattr(os, "sched_getaffinity"):
        try:
            values = sorted(os.sched_getaffinity(0))
            if values:
                return values
        except OSError:
            pass
    return list(range(os.cpu_count() or 1))


def read_optional(path):
    try:
        return Path(path).read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError):
        return None


def cpu_topology(cpu):
    root = Path("/sys/devices/system/cpu/cpu{}".format(cpu))
    result = {"cpu": cpu}
    for key, relative in (
            ("core_id", "topology/core_id"),
            ("package_id", "topology/physical_package_id"),
            ("thread_siblings", "topology/thread_siblings_list"),
            ("governor", "cpufreq/scaling_governor"),
            ("energy_preference", "cpufreq/energy_performance_preference")):
        value = read_optional(root / relative)
        if value is not None:
            result[key] = value
    return result


def machine_identity(cpu):
    model = None
    cpuinfo = read_optional("/proc/cpuinfo")
    if cpuinfo:
        for block in cpuinfo.split("\n\n"):
            fields = {}
            for line in block.splitlines():
                if ":" in line:
                    key, value = line.split(":", 1)
                    fields[key.strip()] = value.strip()
            if fields.get("processor") == str(cpu):
                model = fields.get("model name", fields.get("Hardware"))
                break
    return {
        "allowed_cpus": allowed_cpus(),
        "architecture": platform.machine(),
        "cpu_model": model,
        "pin": cpu_topology(cpu),
        "platform": platform.platform(),
        "python": platform.python_version(),
    }


def parse_cache(path):
    entries = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise CrossoverError("cannot read {}: {}".format(path, error))
    for line in lines:
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        typed, value = line.split("=", 1)
        key = typed.split(":", 1)[0]
        if key in (
                "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER",
                "CMAKE_CXX_COMPILER_ID", "CMAKE_CXX_COMPILER_VERSION",
                "CMAKE_CXX_FLAGS", "CMAKE_GENERATOR",
                "CMAKE_HOME_DIRECTORY", "LEO2_BACKEND_VARIANT",
                "LEO2_BUILD_BENCHMARKS", "LEO2_BUILD_TESTS",
                "LEO2_ENABLE_CUDA") or key.startswith("CMAKE_CXX_FLAGS_"):
            entries[key] = value
    return entries


def executable_identity(path):
    executable = Path(path).resolve()
    if not executable.is_file() or not os.access(str(executable), os.X_OK):
        raise CrossoverError("benchmark is not executable: {}".format(executable))
    cache = None
    directory = executable.parent
    for root in (directory, directory.parent, directory.parent.parent):
        candidate = root / "CMakeCache.txt"
        if candidate.is_file():
            cache = candidate
            break
    if cache is None:
        raise CrossoverError(
            "CMakeCache.txt not found near {}".format(executable))
    entries = parse_cache(cache)
    source = entries.get("CMAKE_HOME_DIRECTORY")
    commit = None
    dirty = None
    if source and (Path(source) / ".git").exists():
        completed = subprocess.run(
            ["git", "-C", source, "status", "--porcelain=v1"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if completed.returncode == 0:
            dirty = bool(completed.stdout.strip())
        completed = subprocess.run(
            ["git", "-C", source, "rev-parse", "HEAD"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if completed.returncode == 0:
            commit = completed.stdout.decode("ascii").strip()
    return {
        "build_cache": str(cache.resolve()),
        "build_cache_sha256": digest_bytes(cache.read_bytes()),
        "cmake": entries,
        "executable": str(executable),
        "executable_sha256": digest_bytes(executable.read_bytes()),
        "source_commit": commit,
        "source_dirty": dirty,
    }


def cell(region, backend, k, r, shard_bytes, batch):
    return {
        "backend": backend,
        "batch": batch,
        "field": "gf8",
        "k": k,
        "losses": 1,
        "profile": "high",
        "r": r,
        "region": region,
        "shard_bytes": shard_bytes,
    }


def compact_grid(backends):
    targets = (
        (3, 1, 1), (3, 1, 17), (3, 1, 65),
        (9, 1, 1024), (9, 8, 4096), (9, 64, 4096),
        (10, 1, 65), (32, 1, 65536), (32, 1, 1048576),
        (240, 1, 65536),
    )
    neighbors_k2 = ((2, 1, 1), (2, 1, 65), (2, 1, 1048576))
    neighbors_r2 = (
        (9, 1, 65), (9, 1, 4096), (9, 8, 4096),
        (32, 1, 65536),
    )
    result = []
    for backend in backends:
        for k, batch, shard_bytes in targets:
            result.append(cell("target_r1", backend, k, 1,
                               shard_bytes, batch))
        for k, batch, shard_bytes in neighbors_k2:
            result.append(cell("neighbor_k2", backend, k, 1,
                               shard_bytes, batch))
        for k, batch, shard_bytes in neighbors_r2:
            result.append(cell("neighbor_r2", backend, k, 2,
                               shard_bytes, batch))
    return result


def full_grid(backends):
    result = compact_grid(backends)
    seen = {digest_value(item) for item in result}
    for backend in backends:
        for k in (2, 3, 4, 9, 10, 31, 32, 33, 64, 127, 128, 240):
            for shard_bytes in (
                    1, 2, 3, 7, 15, 16, 17, 31, 32, 33,
                    63, 64, 65, 127, 128, 129, 257, 1024,
                    4096, 65536, 1048576):
                item = cell("neighbor_k2" if k == 2 else "target_r1",
                            backend, k, 1, shard_bytes, 1)
                identity = digest_value(item)
                if identity not in seen:
                    seen.add(identity)
                    result.append(item)
        for batch in (8, 64):
            for k, shard_bytes in ((3, 17), (9, 1024), (32, 4096)):
                item = cell("target_r1", backend, k, 1, shard_bytes, batch)
                identity = digest_value(item)
                if identity not in seen:
                    seen.add(identity)
                    result.append(item)
    return result


def tiny_threshold_grid(backends):
    """R=1 selector boundaries and unaffected binary-layout controls.

    K=2 exercises the dense two-input reducer and K=7 the exact-arity
    fused-final reducer.  Adjacent K values and both sides of the byte
    thresholds make a compiler-layout or overly broad selector change visible
    instead of attributing every small-code movement to the target kernel.
    """
    result = []
    counts = (2, 3, 4, 6, 7, 8)
    shard_sizes = (1024, 2047, 2048, 2049, 3072, 4095, 4096)
    for backend in backends:
        for k in counts:
            for shard_bytes in shard_sizes:
                target = k in (2, 7) and 2048 <= shard_bytes < 4096
                result.append(cell(
                    "target_r1" if target else "neighbor_r1",
                    backend, k, 1, shard_bytes, 1))
    return result


def final_remainders_grid(backends):
    result = []
    target_counts = (7, 12, 13, 14, 15, 20, 21, 22, 23)
    counts = (6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
              17, 18, 19, 20, 21, 22, 23, 24, 31, 255)
    for backend in backends:
        for k in counts:
            shard_sizes = [
                1024, 1041, 4095, 4096, 4097, 4113,
                65536, 65553, 1048576, 1048593]
            if k == 7:
                shard_sizes.extend((
                    4194304, 4194321,
                    8388608, 8388625,
                    16777216, 16777233))
            for shard_bytes in shard_sizes:
                result.append(cell(
                    "target_r1" if k in target_counts else "neighbor_r1",
                    backend, k, 1, shard_bytes, 1))
    return result


def final_remainders_confirmation_grid(backends):
    """Focused repeat grid for final-remainder cache/layout attribution.

    The broad screen found a K=7 crossover between one and eight MiB plus
    bimodal one-MiB observations around K=20..23.  It also observed a
    surprising aligned K=255 regression even though that codec cannot select
    the candidate callback.  Keep those target and untouched-control cells in
    one deterministic grid so a resumed campaign cannot silently omit the
    negative evidence.
    """
    result = []
    target_counts = {7, 12, 13, 14, 15, 20, 21, 22, 23}
    k7_sizes = (
        65536, 65553,
        1048576, 1048593,
        2097152, 2097169,
        3145728, 3145745,
        4194304, 4194321,
        6291456, 6291473,
        8388608, 8388625,
        16777216, 16777233,
    )
    one_mib_counts = (19, 20, 21, 22, 23, 24)
    one_mib_sizes = (1048576, 1048593)
    positive_anchors = (
        (14, 65536), (14, 65553),
        (20, 65536), (20, 65553),
        (22, 65536), (22, 65553),
        (23, 65536), (23, 65553),
    )
    untouched_layout_controls = (
        (16, 65536), (16, 1048576),
        (254, 1048575), (254, 1048576),
        (254, 1048577), (254, 1048593),
        (255, 1048575), (255, 1048576),
        (255, 1048577), (255, 1048593),
    )
    for backend in backends:
        for shard_bytes in k7_sizes:
            result.append(cell(
                "target_r1", backend, 7, 1, shard_bytes, 1))
        for k in one_mib_counts:
            for shard_bytes in one_mib_sizes:
                result.append(cell(
                    "target_r1" if k in target_counts else "neighbor_r1",
                    backend, k, 1, shard_bytes, 1))
        for k, shard_bytes in positive_anchors:
            result.append(cell(
                "target_r1", backend, k, 1, shard_bytes, 1))
        for k, shard_bytes in untouched_layout_controls:
            result.append(cell(
                "neighbor_r1", backend, k, 1, shard_bytes, 1))
    return result


def final_remainders_production_grid(backends):
    """Final selector boundaries, arbitrary tails, and immediate controls."""
    result = []
    target_counts = (12, 13, 14, 15, 22, 23)
    target_sizes = (
        1024, 2048,
        4095, 4096, 4097, 4113,
        65536, 65553,
        1048575, 1048576, 1048577, 1048593,
    )
    k7_sizes = (
        1024, 2048,
        4095, 4096, 4097, 4113,
        65536, 65553,
        1048576, 1048593,
        2097152, 3145728,
        4194303, 4194304, 4194305, 4194321,
        8388608,
    )
    control_cases = (
        (6, 4096), (8, 4096), (11, 4096), (16, 4096),
        (20, 4096), (21, 4096), (24, 4096),
        (6, 65553), (8, 65553), (11, 65553), (16, 65553),
        (20, 65553), (21, 65553), (24, 65553),
        (6, 1048576), (8, 1048576), (11, 1048576),
        (16, 1048576), (20, 1048576), (21, 1048576),
        (24, 1048576),
    )
    for backend in backends:
        for shard_bytes in k7_sizes:
            result.append(cell(
                "target_r1" if 4096 <= shard_bytes <= 4194304
                else "neighbor_r1",
                backend, 7, 1, shard_bytes, 1))
        for k in target_counts:
            for shard_bytes in target_sizes:
                result.append(cell(
                    "target_r1" if 4096 <= shard_bytes <= 1048576
                    else "neighbor_r1",
                    backend, k, 1, shard_bytes, 1))
        for k, shard_bytes in control_cases:
            result.append(cell(
                "neighbor_r1", backend, k, 1, shard_bytes, 1))
    return result


def parse_backends(value):
    result = []
    for item in value.split(","):
        item = item.strip().lower()
        if item and item not in result:
            result.append(item)
    if not result or any(item not in ACCEPTED_BACKENDS for item in result):
        raise CrossoverError(
            "backends must be a subset of {}".format(
                ",".join(ACCEPTED_BACKENDS)
            )
        )
    return result


def stable_seed(value):
    return int(digest_value({"seed": value})[:16], 16) or 1


def choose_reuse(item, target_work_bytes, maximum_reuse):
    offered = (item["k"] + item["r"]) * item["shard_bytes"] * item["batch"]
    return max(1, min(maximum_reuse, target_work_bytes // max(1, offered)))


def make_manifest(args):
    backends = parse_backends(args.backends)
    if args.grid == "full":
        cells = full_grid(backends)
    elif args.grid == "tiny-thresholds":
        cells = tiny_threshold_grid(backends)
    elif args.grid == "final-remainders":
        cells = final_remainders_grid(backends)
    elif args.grid == "final-remainders-confirmation":
        cells = final_remainders_confirmation_grid(backends)
    elif args.grid == "final-remainders-production":
        cells = final_remainders_production_grid(backends)
    else:
        cells = compact_grid(backends)
    cells = sorted(cells, key=lambda item: (
        item["backend"], item["region"], item["r"], item["k"],
        item["shard_bytes"], item["batch"]))
    executables = {
        "baseline": executable_identity(args.baseline),
        "candidate": executable_identity(args.candidate),
    }
    settings = {
        "abba_rounds": args.abba_rounds,
        "backends": backends,
        "cpu": args.cpu,
        "grid": args.grid,
        "iterations": args.iterations,
        "maximum_reuse": args.maximum_reuse,
        "target_work_bytes": args.target_work_bytes,
        "timeout_seconds": args.timeout,
        "warmup": args.warmup,
    }
    jobs = []
    for item in cells:
        identity = {
            "cell": item,
            "executables": executables,
            "settings": settings,
        }
        job_id = digest_value(identity)[:24]
        jobs.append({
            "cell": item,
            "job_id": job_id,
            "order": ["baseline", "candidate", "candidate", "baseline"] *
                     args.abba_rounds,
            "reuse": choose_reuse(
                item, args.target_work_bytes, args.maximum_reuse),
            "seed": stable_seed(item),
        })
    manifest = {
        "schema": SCHEMA,
        "executables": executables,
        "jobs": jobs,
        "machine": machine_identity(args.cpu),
        "settings": settings,
    }
    manifest["configuration_id"] = digest_value(manifest)
    return manifest


def validate_benchmark(raw, job):
    if raw.get("schema") != BENCHMARK_SCHEMA:
        raise CrossoverError("benchmark emitted an unknown schema")
    item = job["cell"]
    parameters = raw.get("parameters", {})
    selector_names = ("force_tiled_decode", "force_materialized_decode")
    selector_presence = tuple(name in parameters for name in selector_names)
    if selector_presence == (True, True):
        for name in selector_names:
            if type(parameters[name]) is not bool or parameters[name] is not False:
                raise CrossoverError(
                    "benchmark workspace selector {} is not false Boolean".format(name))
    elif selector_presence != (False, False):
        raise CrossoverError("benchmark workspace selector pair is partial")
    expected = {
        "K": item["k"], "R": item["r"], "batch": item["batch"],
        "force_generic_decode": False,
        "force_specialized_decode": False,
        "iterations": job["iterations"], "loss_count": item["losses"],
        "requested_backend": item["backend"],
        "requested_field": item["field"],
        "requested_profile": "legacy_high_v1",
        "retain_samples": True, "reuse": job["reuse"],
        "seed": job["seed"], "shard_bytes": item["shard_bytes"],
        "skip_legacy": True, "thread_count": 1,
        "warmup": job["warmup"],
    }
    for key, value in expected.items():
        actual = parameters.get(key)
        if ((type(value) is bool and
             (type(actual) is not bool or actual is not value)) or
                (type(value) is not bool and actual != value)):
            raise CrossoverError(
                "benchmark parameter {} is {!r}, expected {!r}".format(
                    key, actual, value))
    resolved = raw.get("resolved", {})
    if resolved.get("backend") != item["backend"]:
        raise CrossoverError("benchmark resolved the wrong backend")
    if resolved.get("profile") != "legacy_high_v1" or \
            resolved.get("field") != item["field"]:
        raise CrossoverError("benchmark resolved the wrong wire profile")
    if raw.get("correctness", {}).get("leopard2_round_trip") is not True:
        raise CrossoverError("benchmark round-trip check failed")
    metrics = raw.get("metrics", {})
    values = {}
    for name in ("encode_execution", "decode_execution"):
        try:
            value = float(metrics[name]["median_us_per_batch_call"])
            mad = float(metrics[name]["mad_us_per_batch_call"])
        except (KeyError, TypeError, ValueError):
            raise CrossoverError("benchmark omitted {} timing".format(name))
        if value <= 0 or mad < 0:
            raise CrossoverError("benchmark emitted invalid timing")
        values[name] = {"median_us": value, "mad_us": mad}
    digests = raw.get("workload_digests")
    if not isinstance(digests, dict) or not all(
            isinstance(digests.get(key), str) for key in
            ("original_data", "transmitted_parity", "recovered_originals")):
        raise CrossoverError("benchmark omitted workload digests")
    return values, digests


def command_for(manifest, job, variant, raw_path):
    item = job["cell"]
    executable = manifest["executables"][variant]["executable"]
    return [
        "taskset", "-c", str(manifest["settings"]["cpu"]), executable,
        "--k", str(item["k"]), "--r", str(item["r"]),
        "--profile", item["profile"], "--field", item["field"],
        "--backend", item["backend"], "--bytes", str(item["shard_bytes"]),
        "--loss", str(item["losses"]), "--batch", str(item["batch"]),
        "--reuse", str(job["reuse"]), "--iterations", str(job["iterations"]),
        "--warmup", str(job["warmup"]), "--threads", "1",
        "--seed", str(job["seed"]), "--skip-legacy", "--retain-samples",
        "--json", str(raw_path),
    ]


def job_with_settings(manifest, source_job):
    job = dict(source_job)
    job["iterations"] = manifest["settings"]["iterations"]
    job["warmup"] = manifest["settings"]["warmup"]
    return job


def run_job(root, manifest, source_job):
    job = job_with_settings(manifest, source_job)
    log_dir = root / "logs" / job["job_id"]
    raw_dir = root / "raw" / job["job_id"]
    log_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    records = []
    reference_digests = None
    environment = os.environ.copy()
    environment["OMP_DYNAMIC"] = "FALSE"
    environment["OMP_NUM_THREADS"] = "1"
    for sequence, variant in enumerate(job["order"]):
        stem = "{:02d}-{}".format(sequence, variant)
        raw_path = raw_dir / (stem + ".json")
        stdout_path = log_dir / (stem + ".stdout")
        stderr_path = log_dir / (stem + ".stderr")
        argv = command_for(manifest, job, variant, raw_path)
        try:
            completed = subprocess.run(
                argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env=environment, timeout=manifest["settings"]["timeout_seconds"])
            returncode = completed.returncode
            stdout = normalized(completed.stdout)
            stderr = normalized(completed.stderr)
            timed_out = False
        except subprocess.TimeoutExpired as error:
            returncode = 124
            stdout = normalized(error.stdout or b"")
            stderr = normalized(error.stderr or b"")
            timed_out = True
        stdout_path.write_bytes(stdout)
        stderr_path.write_bytes(stderr)
        if returncode != 0 or timed_out:
            raise CrossoverError(
                "{} failed; see {}".format(job["job_id"], stderr_path))
        try:
            raw_bytes = raw_path.read_bytes()
            raw = json.loads(raw_bytes.decode("utf-8"))
        except (OSError, UnicodeError, ValueError) as error:
            raise CrossoverError("cannot parse {}: {}".format(raw_path, error))
        metrics, digests = validate_benchmark(raw, job)
        if reference_digests is None:
            reference_digests = digests
        elif canonical(reference_digests) != canonical(digests):
            raise CrossoverError("baseline/candidate workload digests differ")
        records.append({
            "argv": argv,
            "benchmark_json": str(raw_path.relative_to(root)),
            "benchmark_json_sha256": digest_bytes(raw_bytes),
            "metrics": metrics,
            "sequence": sequence,
            "stderr": str(stderr_path.relative_to(root)),
            "stderr_sha256": digest_bytes(stderr),
            "stdout": str(stdout_path.relative_to(root)),
            "stdout_sha256": digest_bytes(stdout),
            "variant": variant,
        })
    result = {
        "schema": JOB_SCHEMA,
        "cell": job["cell"],
        "configuration_id": manifest["configuration_id"],
        "job_id": job["job_id"],
        "measurements": records,
        "reuse": job["reuse"],
        "seed": job["seed"],
        "status": "passed",
        "workload_digests": reference_digests,
    }
    atomic_json(root / "jobs" / (job["job_id"] + ".json"), result)
    return result


def safe_artifact(root, relative):
    candidate = (root / relative).resolve()
    if os.path.commonpath((str(root.resolve()), str(candidate))) != str(root.resolve()):
        raise CrossoverError("artifact path escapes result directory")
    return candidate


def valid_existing(root, manifest, job):
    path = root / "jobs" / (job["job_id"] + ".json")
    if not path.is_file():
        return False
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        if value.get("schema") != JOB_SCHEMA or value.get("status") != "passed":
            return False
        if value.get("configuration_id") != manifest["configuration_id"] or \
                value.get("job_id") != job["job_id"] or \
                canonical(value.get("cell")) != canonical(job["cell"]):
            return False
        if value.get("reuse") != job["reuse"] or \
                value.get("seed") != job["seed"]:
            return False
        runtime_job = job_with_settings(manifest, job)
        measurements = value.get("measurements", [])
        if len(measurements) != len(job["order"]):
            return False
        reference_digests = None
        for index, measurement in enumerate(measurements):
            variant = job["order"][index]
            stem = "{:02d}-{}".format(index, variant)
            expected_paths = {
                "benchmark_json": "raw/{}/{}.json".format(job["job_id"], stem),
                "stdout": "logs/{}/{}.stdout".format(job["job_id"], stem),
                "stderr": "logs/{}/{}.stderr".format(job["job_id"], stem),
            }
            if measurement.get("sequence") != index or \
                    measurement.get("variant") != variant:
                return False
            artifacts = {}
            for key in ("benchmark_json", "stdout", "stderr"):
                if measurement.get(key) != expected_paths[key]:
                    return False
                artifact = safe_artifact(root, measurement[key])
                data = artifact.read_bytes()
                if digest_bytes(data) != measurement[key + "_sha256"]:
                    return False
                artifacts[key] = (artifact, data)
            raw_path, raw_bytes = artifacts["benchmark_json"]
            raw = json.loads(raw_bytes.decode("utf-8"))
            metrics, digests = validate_benchmark(raw, runtime_job)
            if canonical(metrics) != canonical(measurement.get("metrics")):
                return False
            expected_argv = command_for(
                manifest, runtime_job, variant, raw_path)
            if canonical(expected_argv) != canonical(measurement.get("argv")):
                return False
            if reference_digests is None:
                reference_digests = digests
            elif canonical(reference_digests) != canonical(digests):
                return False
        if canonical(reference_digests) != \
                canonical(value.get("workload_digests")):
            return False
        return True
    except (KeyError, OSError, TypeError, UnicodeError, ValueError,
            CrossoverError):
        return False


def command_run(args):
    if shutil.which("taskset") is None:
        raise CrossoverError("taskset is required for pinned measurements")
    if args.cpu not in allowed_cpus():
        raise CrossoverError("CPU {} is outside process affinity".format(args.cpu))
    if args.abba_rounds < 2 or args.iterations < 3 or args.warmup < 1:
        raise CrossoverError("use at least 2 ABBA rounds, 3 iterations, and 1 warmup")
    manifest = make_manifest(args)
    root = Path(args.result_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "manifest.json"
    if manifest_path.is_file():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if canonical(existing) != canonical(manifest):
            raise CrossoverError(
                "existing manifest differs; choose a fresh result directory")
    else:
        atomic_json(manifest_path, manifest)
    total = len(manifest["jobs"])
    for index, job in enumerate(manifest["jobs"]):
        if valid_existing(root, manifest, job):
            print("[{}/{}] resume {}".format(index + 1, total, job["job_id"]),
                  flush=True)
            continue
        print("[{}/{}] run {} {} K={} R={} B={} batch={}".format(
            index + 1, total, job["job_id"], job["cell"]["backend"],
            job["cell"]["k"], job["cell"]["r"],
            job["cell"]["shard_bytes"], job["cell"]["batch"]), flush=True)
        run_job(root, manifest, job)
    print("completed {} deterministic ABBA jobs".format(total))


def percentile(values, fraction):
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def confidence_interval(values, seed, draws=10000):
    if not values:
        raise CrossoverError("cannot bootstrap an empty sample")
    generator = random.Random(seed)
    bootstraps = []
    for _ in range(draws):
        sample = [values[generator.randrange(len(values))]
                  for _ in range(len(values))]
        bootstraps.append(statistics.median(sample))
    return percentile(bootstraps, 0.025), percentile(bootstraps, 0.975)


def analyze_metric(result, metric):
    measurements = result["measurements"]
    if len(measurements) % 4 != 0:
        raise CrossoverError("ABBA result has an incomplete round")
    improvements = []
    for offset in range(0, len(measurements), 4):
        round_values = measurements[offset:offset + 4]
        if [item["variant"] for item in round_values] != \
                ["baseline", "candidate", "candidate", "baseline"]:
            raise CrossoverError("result is not in ABBA order")
        baseline = statistics.mean(
            [round_values[0]["metrics"][metric]["median_us"],
             round_values[3]["metrics"][metric]["median_us"]])
        candidate = statistics.mean(
            [round_values[1]["metrics"][metric]["median_us"],
             round_values[2]["metrics"][metric]["median_us"]])
        improvements.append(100.0 * (baseline / candidate - 1.0))
    median = statistics.median(improvements)
    seed = int(digest_value({
        "job": result["job_id"], "metric": metric})[:16], 16)
    lower, upper = confidence_interval(improvements, seed)
    return {
        "abba_round_improvements_percent": improvements,
        "credible_improvement_at_least_5_percent": lower >= 5.0,
        "credible_regression_over_2_percent": upper < -2.0,
        "median_improvement_percent": median,
        "observed_regression_over_2_percent": median < -2.0,
        "percentile_95": [lower, upper],
    }


def evaluate_gate(cells, requested_backends):
    target_wins = {
        backend: {"decode": 0, "encode": 0}
        for backend in requested_backends
    }
    regressions = []
    for item in cells:
        backend = item["cell"]["backend"]
        if backend not in target_wins:
            raise CrossoverError(
                "analysis contains unrequested backend {}".format(backend))
        for metric in ("decode", "encode"):
            if item["cell"]["region"] == "target_r1" and \
                    item[metric]["credible_improvement_at_least_5_percent"]:
                target_wins[backend][metric] += 1
            if item[metric]["observed_regression_over_2_percent"]:
                regressions.append({
                    "cell": item["cell"], "job_id": item["job_id"],
                    "metric": metric,
                    "median_improvement_percent":
                        item[metric]["median_improvement_percent"],
                    "credible": item[metric]["credible_regression_over_2_percent"],
                })
    every_group_wins = all(
        target_wins[backend][metric] >= 1
        for backend in requested_backends
        for metric in ("decode", "encode"))
    no_regressions = not regressions
    return {
        "backend_target_wins": target_wins,
        "credible_target_win_for_every_requested_backend_and_metric":
            every_group_wins,
        "no_observed_neighbor_or_target_regression_over_2_percent":
            no_regressions,
        "observed_regressions": regressions,
        "passed": every_group_wins and no_regressions,
        "promotion_rule": (
            "at least one credible >=5% target win for every requested backend "
            "and both encode/decode; zero observed >2% target/neighbor regressions"),
        "requested_backends": list(requested_backends),
    }


def enforce_gate(gate, reporting_only):
    if not gate["passed"] and not reporting_only:
        raise CrossoverError(
            "promotion gate failed (use --reporting-only to retain a negative report)")


def self_test_resume_validation(checks):
    with tempfile.TemporaryDirectory(prefix="leopard2-r1-xor-self-test-") as name:
        root = Path(name).resolve()
        source_job = {
            "cell": cell("target_r1", "scalar", 3, 1, 17, 1),
            "job_id": "resume-validation",
            "order": ["baseline", "candidate", "candidate", "baseline"],
            "reuse": 7,
            "seed": 42,
        }
        manifest = {
            "configuration_id": "self-test",
            "executables": {
                "baseline": {"executable": str(root / "baseline")},
                "candidate": {"executable": str(root / "candidate")},
            },
            "settings": {
                "backends": ["scalar"], "cpu": 0,
                "iterations": 3, "warmup": 1,
            },
        }
        job = job_with_settings(manifest, source_job)
        records = []
        digests = {
            "algorithm": "fnv1a64",
            "original_data": "01",
            "recovered_originals": "02",
            "transmitted_parity": "03",
        }
        for sequence, variant in enumerate(job["order"]):
            stem = "{:02d}-{}".format(sequence, variant)
            raw_path = root / "raw" / job["job_id"] / (stem + ".json")
            stdout_path = root / "logs" / job["job_id"] / (stem + ".stdout")
            stderr_path = root / "logs" / job["job_id"] / (stem + ".stderr")
            raw = {
                "schema": BENCHMARK_SCHEMA,
                "parameters": {
                    "K": 3, "R": 1, "batch": 1, "iterations": 3,
                    "force_generic_decode": False,
                    "force_specialized_decode": False,
                    "loss_count": 1, "requested_backend": "scalar",
                    "requested_field": "gf8",
                    "requested_profile": "legacy_high_v1",
                    "retain_samples": True, "reuse": 7, "seed": 42,
                    "shard_bytes": 17, "skip_legacy": True,
                    "thread_count": 1, "warmup": 1,
                },
                "resolved": {
                    "backend": "scalar", "field": "gf8",
                    "profile": "legacy_high_v1",
                },
                "correctness": {"leopard2_round_trip": True},
                "metrics": {
                    "encode_execution": {
                        "median_us_per_batch_call": 10.0,
                        "mad_us_per_batch_call": 0.1,
                    },
                    "decode_execution": {
                        "median_us_per_batch_call": 11.0,
                        "mad_us_per_batch_call": 0.2,
                    },
                },
                "workload_digests": digests,
            }
            if variant == "candidate":
                raw["parameters"].update({
                    "force_tiled_decode": False,
                    "force_materialized_decode": False,
                })
            atomic_json(raw_path, raw)
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            stdout_path.write_bytes(b"")
            stderr_path.write_bytes(b"")
            raw_bytes = raw_path.read_bytes()
            metrics, parsed_digests = validate_benchmark(raw, job)
            checks.require(
                canonical(parsed_digests) == canonical(digests),
                "validated workload digests differ from fixture")
            records.append({
                "argv": command_for(manifest, job, variant, raw_path),
                "benchmark_json": str(raw_path.relative_to(root)),
                "benchmark_json_sha256": digest_bytes(raw_bytes),
                "metrics": metrics,
                "sequence": sequence,
                "stderr": str(stderr_path.relative_to(root)),
                "stderr_sha256": digest_bytes(b""),
                "stdout": str(stdout_path.relative_to(root)),
                "stdout_sha256": digest_bytes(b""),
                "variant": variant,
            })
        result_path = root / "jobs" / (job["job_id"] + ".json")
        pristine_result = {
            "schema": JOB_SCHEMA,
            "cell": job["cell"],
            "configuration_id": manifest["configuration_id"],
            "job_id": job["job_id"],
            "measurements": records,
            "reuse": job["reuse"],
            "seed": job["seed"],
            "status": "passed",
            "workload_digests": digests,
        }
        atomic_json(result_path, pristine_result)
        checks.require(
            valid_existing(root, manifest, source_job),
            "pristine resumable job was rejected")

        def rejects_job_mutation(mutator):
            mutated = json.loads(json.dumps(pristine_result))
            mutator(mutated)
            atomic_json(result_path, mutated)
            checks.require(
                not valid_existing(root, manifest, source_job),
                "mutated resumable job was accepted")
            atomic_json(result_path, pristine_result)

        rejects_job_mutation(
            lambda value: value["measurements"][0]["metrics"]
                ["encode_execution"].update({"median_us": 9.0}))
        rejects_job_mutation(
            lambda value: value["measurements"][0]["argv"].append("--tampered"))
        rejects_job_mutation(lambda value: value.update({"seed": 43}))
        rejects_job_mutation(lambda value: value.update({"reuse": 8}))
        rejects_job_mutation(
            lambda value: value["measurements"][0].update({"variant": "candidate"}))
        rejects_job_mutation(
            lambda value: value["workload_digests"].update({"original_data": "ff"}))

        raw_path = root / records[0]["benchmark_json"]
        pristine_raw = json.loads(raw_path.read_text(encoding="utf-8"))

        for force_selector in (
                "force_generic_decode", "force_specialized_decode"):
            for mutation in ("missing", "active", "integer"):
                changed = json.loads(json.dumps(pristine_raw))
                if mutation == "missing":
                    changed["parameters"].pop(force_selector)
                elif mutation == "active":
                    changed["parameters"][force_selector] = True
                else:
                    changed["parameters"][force_selector] = 0
                try:
                    validate_benchmark(changed, job)
                except CrossoverError:
                    pass
                else:
                    raise CrossoverError(
                        "{} {} mutation was accepted".format(
                            force_selector, mutation))

        selector_fixture = json.loads(json.dumps(pristine_raw))
        selector_fixture["parameters"].update({
            "force_tiled_decode": False,
            "force_materialized_decode": False,
        })
        validate_benchmark(selector_fixture, job)
        for selector in ("force_tiled_decode", "force_materialized_decode"):
            changed = json.loads(json.dumps(selector_fixture))
            changed["parameters"].pop(selector)
            try:
                validate_benchmark(changed, job)
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "partial workspace selector pair was accepted")
            changed = json.loads(json.dumps(selector_fixture))
            changed["parameters"][selector] = True
            try:
                validate_benchmark(changed, job)
            except CrossoverError:
                pass
            else:
                raise CrossoverError(
                    "active workspace selector was accepted")

        def rejects_raw_mutation(mutator):
            mutated = json.loads(json.dumps(pristine_raw))
            mutator(mutated)
            atomic_json(raw_path, mutated)
            mutated_result = json.loads(json.dumps(pristine_result))
            mutated_result["measurements"][0]["benchmark_json_sha256"] = \
                digest_bytes(raw_path.read_bytes())
            atomic_json(result_path, mutated_result)
            checks.require(
                not valid_existing(root, manifest, source_job),
                "mutated raw benchmark was accepted")
            atomic_json(raw_path, pristine_raw)
            atomic_json(result_path, pristine_result)

        # Updating the cached file hash cannot conceal invalid benchmark
        # parameters or a mismatch with cached derived metrics/digests.
        rejects_raw_mutation(
            lambda value: value["parameters"].update({"seed": 43}))
        rejects_raw_mutation(
            lambda value: value["metrics"]["encode_execution"].update(
                {"median_us_per_batch_call": 9.0}))
        rejects_raw_mutation(
            lambda value: value["workload_digests"].update(
                {"original_data": "ff"}))


def command_analyze(args):
    root = Path(args.result_dir).resolve()
    try:
        manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise CrossoverError("cannot read manifest: {}".format(error))
    cells = []
    missing = []
    for job in manifest.get("jobs", []):
        if not valid_existing(root, manifest, job):
            missing.append(job["job_id"])
            continue
        result = json.loads((root / "jobs" / (job["job_id"] + ".json"))
                            .read_text(encoding="utf-8"))
        cells.append({
            "cell": result["cell"],
            "decode": analyze_metric(result, "decode_execution"),
            "encode": analyze_metric(result, "encode_execution"),
            "job_id": result["job_id"],
        })
    if missing:
        raise CrossoverError("incomplete jobs: {}".format(", ".join(missing)))
    gate = evaluate_gate(cells, manifest["settings"]["backends"])
    analysis = {
        "schema": ANALYSIS_SCHEMA,
        "cells": cells,
        "configuration_id": manifest["configuration_id"],
        "gate": gate,
    }
    atomic_json(root / "analysis.json", analysis)
    print(json.dumps(gate, indent=2, sort_keys=True))
    enforce_gate(gate, args.reporting_only)


def command_self_test(_args):
    checks = SelfTestChecks()
    checks.require(
        parse_backends("avx512") == ["avx512"],
        "AVX-512 backend parser result")
    defaults = parser().parse_args([
        "run", "--baseline", "baseline", "--candidate", "candidate",
        "--result-dir", "results", "--cpu", "0",
    ])
    checks.require(
        defaults.backends == ",".join(BACKENDS),
        "default backend list")
    values = compact_grid(BACKENDS)
    checks.require(
        {item["backend"] for item in values} == set(BACKENDS),
        "compact grid backend coverage")
    checks.require(
        min(item["shard_bytes"] for item in values) == 1,
        "compact grid minimum shard bytes")
    checks.require(
        max(item["shard_bytes"] for item in values) == 1048576,
        "compact grid maximum shard bytes")
    checks.require(
        max(item["k"] for item in values) >= 240,
        "compact grid large-K coverage")
    checks.require(
        {item["batch"] for item in values}.issuperset({1, 8, 64}),
        "compact grid batch coverage")
    confirmation = final_remainders_confirmation_grid(["avx2"])
    checks.require(
        len(confirmation) == 46,
        "final-remainders confirmation grid size")
    checks.require(
        any(item["k"] == 7 and item["shard_bytes"] == 3145728
            for item in confirmation),
        "final-remainders K=7 confirmation cell")
    checks.require(
        any(item["k"] == 255 and item["shard_bytes"] == 1048576
            for item in confirmation),
        "final-remainders K=255 confirmation cell")
    production = final_remainders_production_grid(["avx2"])
    checks.require(
        len(production) == 110,
        "final-remainders production grid size")
    for k in (7, 12, 13, 14, 15, 22, 23):
        for shard_bytes in (1024, 2048):
            control = next(item for item in production
                           if item["k"] == k and
                           item["shard_bytes"] == shard_bytes)
            checks.require(
                control["region"] == "neighbor_r1",
                "small-shard K={} neighbor classification".format(k))
    checks.require(
        any(item["k"] == 7 and item["shard_bytes"] == 4194305
            for item in production),
        "K=7 upper-bound neighbor presence")
    checks.require(
        next(item for item in production
             if item["k"] == 7 and item["shard_bytes"] == 4194305)[
                 "region"] == "neighbor_r1",
        "K=7 upper-bound neighbor classification")
    for k in (12, 13, 14, 15, 22, 23):
        checks.require(
            any(item["k"] == k and item["shard_bytes"] == 1048577
                for item in production),
            "K={} upper-bound neighbor presence".format(k))
        checks.require(
            next(item for item in production
                 if item["k"] == k and
                 item["shard_bytes"] == 1048577)["region"] ==
                    "neighbor_r1",
            "K={} upper-bound neighbor classification".format(k))
    tiny = tiny_threshold_grid(["avx2"])
    checks.require(len(tiny) == 42, "tiny-threshold grid size")
    for k in (2, 7):
        for shard_bytes in (2048, 2049, 3072, 4095):
            checks.require(
                next(item for item in tiny
                     if item["k"] == k and
                     item["shard_bytes"] == shard_bytes)["region"] ==
                    "target_r1",
                "tiny-threshold target classification")
    for k in (3, 4, 6, 8):
        checks.require(
            all(item["region"] == "neighbor_r1"
                for item in tiny if item["k"] == k),
            "tiny-threshold control classification")
    fake = {
        "job_id": "0" * 24,
        "measurements": [],
    }
    for _ in range(4):
        for variant, timing in (
                ("baseline", 10.0), ("candidate", 8.0),
                ("candidate", 8.0), ("baseline", 10.0)):
            fake["measurements"].append({
                "variant": variant,
                "metrics": {
                    "encode_execution": {"median_us": timing},
                    "decode_execution": {"median_us": timing},
                },
            })
    result = analyze_metric(fake, "encode_execution")
    checks.require(
        result["median_improvement_percent"] == 25.0,
        "paired improvement calculation")
    checks.require(
        result["credible_improvement_at_least_5_percent"],
        "paired confidence calculation")
    positive_cells = []
    for backend in BACKENDS:
        metric = {
            "credible_improvement_at_least_5_percent": True,
            "credible_regression_over_2_percent": False,
            "median_improvement_percent": 10.0,
            "observed_regression_over_2_percent": False,
        }
        positive_cells.append({
            "cell": cell("target_r1", backend, 3, 1, 17, 1),
            "decode": dict(metric), "encode": dict(metric),
            "job_id": backend,
        })
    positive_gate = evaluate_gate(positive_cells, BACKENDS)
    checks.require(
        positive_gate["passed"],
        "positive promotion gate")

    zero_win_cells = json.loads(json.dumps(positive_cells))
    zero_win_cells[0]["encode"][
        "credible_improvement_at_least_5_percent"] = False
    zero_win_gate = evaluate_gate(zero_win_cells, BACKENDS)
    checks.require(
        not zero_win_gate["passed"],
        "zero-win promotion rejection")
    try:
        enforce_gate(zero_win_gate, False)
        raise AssertionError("zero-win gate unexpectedly passed")
    except CrossoverError:
        pass
    enforce_gate(zero_win_gate, True)

    regression_cells = json.loads(json.dumps(positive_cells))
    regression_cells[0]["decode"]["median_improvement_percent"] = -3.0
    regression_cells[0]["decode"]["observed_regression_over_2_percent"] = True
    regression_gate = evaluate_gate(regression_cells, BACKENDS)
    checks.require(
        not regression_gate["passed"],
        "regression promotion rejection")
    try:
        enforce_gate(regression_gate, False)
        raise AssertionError("regression gate unexpectedly passed")
    except CrossoverError:
        pass
    enforce_gate(regression_gate, True)
    self_test_resume_validation(checks)
    checks.require(
        choose_reuse(cell("x", "scalar", 9, 1, 4096, 1),
                     1024 * 1024, 100000) > 1,
        "reuse calculation")
    print(
        "Leopard2 R=1 XOR crossover self-test passed: jobs={} checks={}".format(
            len(values), checks.count))


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command")
    run = subparsers.add_parser("run", help="run/resume pinned ABBA jobs")
    run.add_argument("--baseline", required=True)
    run.add_argument("--candidate", required=True)
    run.add_argument("--result-dir", required=True)
    run.add_argument("--cpu", type=int, required=True)
    run.add_argument("--backends", default=",".join(BACKENDS))
    run.add_argument(
        "--grid", choices=(
            "compact", "full", "tiny-thresholds", "final-remainders",
            "final-remainders-confirmation", "final-remainders-production",
            ),
        default="compact")
    run.add_argument("--abba-rounds", type=int, default=5)
    run.add_argument("--iterations", type=int, default=9)
    run.add_argument("--warmup", type=int, default=2)
    run.add_argument("--target-work-bytes", type=int, default=268435456)
    run.add_argument("--maximum-reuse", type=int, default=100000)
    run.add_argument("--timeout", type=int, default=180)
    run.set_defaults(function=command_run)
    analyze = subparsers.add_parser("analyze", help="validate and analyze jobs")
    analyze.add_argument("--result-dir", required=True)
    analyze.add_argument(
        "--reporting-only", action="store_true",
        help="write a negative analysis without failing the command")
    analyze.set_defaults(function=command_analyze)
    self_test = subparsers.add_parser("self-test", help="test the runner")
    self_test.set_defaults(function=command_self_test)
    return result


def main(argv=None):
    arguments = parser().parse_args(argv)
    if not hasattr(arguments, "function"):
        parser().print_help()
        return 2
    try:
        arguments.function(arguments)
        return 0
    except CrossoverError as error:
        print("R=1 XOR crossover error: {}".format(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
