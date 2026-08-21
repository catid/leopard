#!/usr/bin/env python3

"""Frozen same-binary qualification for the larger K65/R65 T128 leaf."""

import argparse
import fcntl
import hashlib
import json
import math
import os
import platform
import stat
import statistics
import subprocess
import sys
import time
import traceback
import types
from pathlib import Path


TARGET_ROUNDS = 25
NEIGHBOR_ROUNDS = 9
SUSTAINED_ROUNDS = 9
ITERATIONS = 31
WARMUP = 64
BASE_REUSE_BYTES = 524288
SUSTAINED_REUSE_MULTIPLIER = 16
MAX_ROUND_ATTEMPTS = 5
CAMPAIGN_DEADLINE_SECONDS = 7200
MAX_RESULT_BYTES = 8 * 1024 * 1024
MAX_VALIDATOR_BYTES = 8 * 1024 * 1024
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
GIT = Path("/usr/bin/git")
T_CRITICAL = {8: 2.306004135204166, 24: 2.0638985616280205}
CHILD_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}
GIT_ENVIRONMENT = dict(CHILD_ENVIRONMENT, **{
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_CONFIG_SYSTEM": "/dev/null",
    "GIT_NO_REPLACE_OBJECTS": "1",
    "GIT_OPTIONAL_LOCKS": "0",
})
CONTRACT = (
    "LEGACY_HIGH_V1,GF8,AUTO,K=65,R=65,T=128,"
    "B=128|256|512|1024|2048|4096,native_layout,"
    "native_cantor_affine,dense_packed_full,64_byte_vector_tiles,"
    "max_microtile_bytes=256,runtime_AVX512F_BW_VL_GFNI,"
    "startup_KAT_B128_B320,calibrated_AMD_1A_08,"
    "one_shot_and_one_item_batch"
)

# id, K, R, bytes, rounds, role, promotion region, reuse multiplier
TARGETS = tuple(
    (f"target-k65-r65-b{shard_bytes}", 65, 65, shard_bytes,
     TARGET_ROUNDS, "target", region, 1)
    for region, byte_values in (
        ("small", (128, 256, 512)),
        ("large", (1024, 2048, 4096)),
    )
    for shard_bytes in byte_values
)
BYTE_NEIGHBORS = tuple(
    (f"neighbor-k65-r65-b{shard_bytes}", 65, 65, shard_bytes,
     NEIGHBOR_ROUNDS, "neighbor", "inactive", 1)
    for shard_bytes in (127, 129, 192, 320, 768, 4095, 4097, 8192)
)
SHAPE_NEIGHBORS = tuple(
    (f"neighbor-k{k}-r{r}-b{shard_bytes}", k, r, shard_bytes,
     NEIGHBOR_ROUNDS, "neighbor", "inactive", 1)
    for shard_bytes in (512, 4096)
    for k, r in ((64, 65), (66, 65), (65, 64), (65, 66))
)
SUSTAINED = tuple(
    (f"sustained-k65-r65-b{shard_bytes}", 65, 65, shard_bytes,
     SUSTAINED_ROUNDS, "sustained", "endpoint", SUSTAINED_REUSE_MULTIPLIER)
    for shard_bytes in (128, 4096)
)
CELLS = TARGETS + BYTE_NEIGHBORS + SHAPE_NEIGHBORS + SUSTAINED


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def metadata_matches(left, right):
    return (
        left.st_dev == right.st_dev and left.st_ino == right.st_ino and
        left.st_mode == right.st_mode and left.st_nlink == right.st_nlink and
        left.st_size == right.st_size and
        left.st_mtime_ns == right.st_mtime_ns and
        left.st_ctime_ns == right.st_ctime_ns)


def stable_file(path, *, return_bytes=False, maximum_bytes=None,
                require_single_link=False, require_read_only=False):
    before = path.lstat()
    require(stat.S_ISREG(before.st_mode), f"not a regular file: {path}")
    require(not path.is_symlink(), f"symlink is not evidence: {path}")
    if require_single_link:
        require(before.st_nlink == 1, f"file has multiple links: {path}")
    if require_read_only:
        require(before.st_mode & 0o222 == 0,
                f"frozen file remains writable: {path}")
    if maximum_bytes is not None:
        require(before.st_size <= maximum_bytes,
                f"file exceeds evidence size bound: {path}")
    descriptor = os.open(
        path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        opened = os.fstat(descriptor)
        require(metadata_matches(before, opened),
                f"file changed while opening: {path}")
        digest = hashlib.sha256()
        retained = bytearray() if return_bytes else None
        total = 0
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            total += len(block)
            if maximum_bytes is not None:
                require(total <= maximum_bytes,
                        f"file exceeds evidence size bound: {path}")
            digest.update(block)
            if retained is not None:
                retained.extend(block)
        after_descriptor = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after_path = path.lstat()
    require(metadata_matches(before, after_descriptor) and
            metadata_matches(before, after_path) and total == before.st_size,
            f"file changed while hashing: {path}")
    identity = {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size": before.st_size,
        "mode": stat.S_IMODE(before.st_mode),
        "links": before.st_nlink,
        "device": before.st_dev,
        "inode": before.st_ino,
        "mtime_ns": before.st_mtime_ns,
        "ctime_ns": before.st_ctime_ns,
    }
    if retained is None:
        return identity
    return bytes(retained), identity


def sha256(path):
    return stable_file(path)["sha256"]


def file_identity(path, require_single_link=False, require_read_only=False):
    return stable_file(
        path,
        require_single_link=require_single_link,
        require_read_only=require_read_only)


def write_bytes_exclusive(path, data):
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise


def strict_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_json_constant(value):
    raise RuntimeError(f"non-finite JSON constant: {value}")


def strict_json_bytes(data):
    text = data.decode("utf-8", errors="strict")
    return json.loads(
        text, object_pairs_hook=strict_object,
        parse_constant=reject_json_constant)


def load_validator(path, expected_identity):
    source, identity = stable_file(
        path, return_bytes=True, maximum_bytes=MAX_VALIDATOR_BYTES,
        require_single_link=True, require_read_only=True)
    require(identity == expected_identity,
            "validator changed before direct source execution")
    module_name = "leopard2_frozen_benchmark_validator"
    module = types.ModuleType(module_name)
    module.__file__ = str(path)
    sys.modules[module_name] = module
    code = compile(source, str(path), "exec", dont_inherit=True, optimize=0)
    exec(code, module.__dict__)
    require(callable(getattr(module, "validate_common", None)) and
            callable(getattr(module, "validate_workload_digests", None)),
            "frozen benchmark validator lacks required entry points")
    return module


def cpu_snapshot(cpu):
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            values = [int(value) for value in line.split()[1:]]
            require(len(values) >= 8, "CPU scheduler counters are incomplete")
            return {
                "cpu": cpu,
                "values": values,
                "idle": values[3] + values[4],
                "nonidle": sum(values[:3]) + sum(values[5:8]),
                "total": sum(values),
            }
    raise RuntimeError(f"CPU {cpu} is absent from /proc/stat")


def cpu_delta(before, after):
    require(before["cpu"] == after["cpu"] and
            len(before["values"]) == len(after["values"]),
            "CPU scheduler counter shape changed")
    values = [new - old for old, new in zip(
        before["values"], after["values"])]
    require(all(value >= 0 for value in values),
            "CPU scheduler counter moved backwards")
    return {
        "cpu": before["cpu"],
        "values": values,
        "idle": after["idle"] - before["idle"],
        "nonidle": after["nonidle"] - before["nonidle"],
        "total": after["total"] - before["total"],
    }


def acquire_lock(inherited_fd):
    if inherited_fd is None:
        descriptor = os.open(
            LOCK_PATH,
            os.O_RDWR | os.O_CREAT | os.O_CLOEXEC | os.O_NOFOLLOW,
            0o600)
        lock_mode = "runner-acquired"
    else:
        descriptor = inherited_fd
    canonical = LOCK_PATH.lstat()
    require(stat.S_ISREG(canonical.st_mode),
            "canonical benchmark lock is not a regular file")
    if inherited_fd is not None:
        inherited = os.fstat(descriptor)
        require((inherited.st_dev, inherited.st_ino) ==
                (canonical.st_dev, canonical.st_ino),
                "inherited lock descriptor is not the canonical lock")
        lock_mode = "inherited-across-build-copy-campaign"
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as error:
        if inherited_fd is None:
            os.close(descriptor)
        raise RuntimeError(f"canonical benchmark lock is busy: {LOCK_PATH}") \
            from error
    locked = os.fstat(descriptor)
    require((locked.st_dev, locked.st_ino) ==
            (canonical.st_dev, canonical.st_ino),
            "canonical benchmark lock changed during acquisition")
    return descriptor, {
        "path": str(LOCK_PATH),
        "mode": lock_mode,
        "descriptor": descriptor,
        "device": locked.st_dev,
        "inode": locked.st_ino,
    }


def verify_lock_continuity(descriptor, identity, cpu):
    canonical = LOCK_PATH.lstat()
    held = os.fstat(descriptor)
    require(stat.S_ISREG(canonical.st_mode) and canonical.st_nlink == 1,
            "canonical benchmark lock is not one regular file")
    require((canonical.st_dev, canonical.st_ino) ==
            (identity["device"], identity["inode"]) ==
            (held.st_dev, held.st_ino),
            "canonical benchmark lock continuity was lost")
    require(set(os.sched_getaffinity(0)) == {cpu},
            "campaign controller affinity changed")
    return {
        "device": held.st_dev,
        "inode": held.st_ino,
        "links": held.st_nlink,
        "affinity": sorted(os.sched_getaffinity(0)),
    }


def git_bytes(repository, arguments, label, maximum_bytes=8 * 1024 * 1024):
    command = [str(GIT), "-C", str(repository)] + list(arguments)
    completed = subprocess.run(
        command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False, timeout=120,
        env=GIT_ENVIRONMENT)
    require(completed.returncode == 0,
            f"{label} failed: " + completed.stderr.decode(
                "utf-8", errors="replace"))
    require(not completed.stderr, f"{label} emitted stderr")
    require(len(completed.stdout) <= maximum_bytes,
            f"{label} output exceeds evidence bound")
    return completed.stdout


def verify_source_binding(repository, source_commit, source_tree,
                          source_archive_identity, runner_identity,
                          validator_identity):
    require(repository.is_absolute(), "repository path is not absolute")
    repository_metadata = repository.lstat()
    require(stat.S_ISDIR(repository_metadata.st_mode) and
            not repository.is_symlink(),
            "repository path is not one real directory")
    root = git_bytes(
        repository, ("rev-parse", "--show-toplevel"),
        "Git top-level query").decode("utf-8").strip()
    require(Path(root).resolve() == repository.resolve(),
            "repository is not the Git top level")
    head = git_bytes(
        repository, ("rev-parse", "HEAD"), "Git HEAD query").decode(
            "ascii").strip()
    tree = git_bytes(
        repository, ("rev-parse", "HEAD^{tree}"),
        "Git tree query").decode("ascii").strip()
    require(head == source_commit and tree == source_tree,
            "live source does not match the declared commit/tree")
    status = git_bytes(
        repository,
        ("status", "--porcelain=v1", "--untracked-files=normal"),
        "Git status query")
    require(not status, "live source is not clean during qualification")

    expected_sources = (
        ("experiments/leopard2/aes_gfni_t128/"
         "run_larger_shards_abba.py", runner_identity),
        ("tools/leopard2_benchmark_json_test.py", validator_identity),
    )
    source_bindings = {}
    for relative, frozen_identity in expected_sources:
        source = git_bytes(
            repository, ("show", f"{source_commit}:{relative}"),
            f"frozen source query for {relative}",
            maximum_bytes=MAX_VALIDATOR_BYTES)
        digest = hashlib.sha256(source).hexdigest()
        require(digest == frozen_identity["sha256"] and
                len(source) == frozen_identity["size"],
                f"frozen file differs from committed source: {relative}")
        source_bindings[relative] = {
            "sha256": digest,
            "size": len(source),
        }

    archive_command = [
        str(GIT), "-C", str(repository), "archive", "--format=tar",
        "--prefix=source/", source_commit,
    ]
    process = subprocess.Popen(
        archive_command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, env=GIT_ENVIRONMENT)
    require(process.stdout is not None and process.stderr is not None,
            "Git archive pipes are unavailable")
    archive_digest = hashlib.sha256()
    archive_size = 0
    for block in iter(lambda: process.stdout.read(1024 * 1024), b""):
        archive_digest.update(block)
        archive_size += len(block)
    archive_stderr = process.stderr.read()
    returncode = process.wait(timeout=120)
    require(returncode == 0 and not archive_stderr,
            "canonical Git archive generation failed: " +
            archive_stderr.decode("utf-8", errors="replace"))
    require(archive_digest.hexdigest() ==
            source_archive_identity["sha256"] and
            archive_size == source_archive_identity["size"],
            "source archive does not equal canonical commit archive")
    return {
        "repository": str(repository),
        "git": file_identity(GIT),
        "head": head,
        "tree": tree,
        "status_porcelain": "",
        "committed_frozen_sources": source_bindings,
        "canonical_archive": {
            "sha256": archive_digest.hexdigest(),
            "size": archive_size,
            "format": "git archive --format=tar --prefix=source/",
        },
    }


def finite_positive(value, label):
    require(type(value) in (int, float), f"{label} is not numeric")
    value = float(value)
    require(math.isfinite(value) and value > 0, f"{label} is not positive")
    return value


def timing(document, name):
    summary = document["metrics"][name]
    samples = summary["samples_us_per_batch_call"]
    require(len(samples) == ITERATIONS, f"{name} sample count changed")
    samples = [finite_positive(value, f"{name} sample") for value in samples]
    observed = statistics.median(samples)
    reported = finite_positive(
        summary["median_us_per_batch_call"], f"{name} median")
    require(abs(observed - reported) <= 0.000003,
            f"{name} median is not derived from retained samples")
    return reported, samples


def append_event(path, event):
    encoded = json.dumps(event, sort_keys=True, separators=(",", ":")) + "\n"
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_APPEND | os.O_CREAT | os.O_CLOEXEC |
        os.O_NOFOLLOW,
        0o600)
    try:
        metadata = os.fstat(descriptor)
        require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
                "campaign journal is not one regular file")
        data = encoded.encode("utf-8")
        written = 0
        while written < len(data):
            count = os.write(descriptor, data[written:])
            require(count > 0, "campaign journal write made no progress")
            written += count
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def report_temporary_path(path):
    return path.with_name(path.name + ".tmp")


def write_report(path, report):
    temporary = report_temporary_path(path)
    data = (json.dumps(report, indent=2, sort_keys=True) + "\n").encode(
        "utf-8")
    write_bytes_exclusive(temporary, data)
    os.replace(temporary, path)
    directory = os.open(path.parent, os.O_RDONLY | os.O_CLOEXEC)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def cell_reuse(cell):
    shard_bytes = cell[3]
    multiplier = cell[7]
    return max(1, BASE_REUSE_BYTES // shard_bytes) * multiplier


def run_launch(binary, cpu, cell, mode, seed, source_commit, source_tree,
               validator, invocations, journal, label, deadline):
    cell_id, k, r, shard_bytes, _, role, _, _ = cell
    reuse = cell_reuse(cell)
    invocation = invocations / label
    invocation.mkdir(mode=0o700)
    output = invocation / "result.json"
    stdout_path = invocation / "stdout"
    stderr_path = invocation / "stderr"
    command_path = invocation / "command.json"
    command = [
        "/usr/bin/taskset", "-c", str(cpu), str(binary),
        "--k", str(k), "--r", str(r),
        "--profile", "high", "--field", "gf8",
        "--backend", "auto", "--bytes", str(shard_bytes),
        "--loss", "1", "--batch", "1", "--reuse", str(reuse),
        "--iterations", str(ITERATIONS), "--warmup", str(WARMUP),
        "--threads", "1", "--seed", str(seed),
        "--skip-legacy", "--retain-samples",
        "--measure-one-shot-encode", "--attest-source",
        "--k65r65-t128-avx512-gfni-mode", str(mode),
        "--json", str(output),
    ]
    write_bytes_exclusive(
        command_path,
        (json.dumps(command, separators=(",", ":")) + "\n").encode("utf-8"))
    append_event(journal, {
        "event": "launch_started",
        "label": label,
        "cell": cell_id,
        "mode": mode,
        "command": command,
        "invocation": str(invocation),
    })
    started = time.monotonic_ns()
    remaining = deadline - time.monotonic()
    require(remaining > 0, "campaign deadline expired before launch")
    try:
        completed = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False, timeout=min(180, remaining), env=CHILD_ENVIRONMENT)
    except subprocess.TimeoutExpired as error:
        write_bytes_exclusive(stdout_path, error.stdout or b"")
        write_bytes_exclusive(stderr_path, error.stderr or b"")
        append_event(journal, {
            "event": "launch_failed",
            "label": label,
            "reason": "timeout",
        })
        raise RuntimeError(f"benchmark timed out: {label}") from error
    elapsed_ns = time.monotonic_ns() - started
    write_bytes_exclusive(stdout_path, completed.stdout)
    write_bytes_exclusive(stderr_path, completed.stderr)
    if completed.returncode != 0:
        append_event(journal, {
            "event": "launch_failed",
            "label": label,
            "reason": "nonzero",
            "returncode": completed.returncode,
        })
        raise RuntimeError(
            "benchmark failed: " + completed.stdout.decode(
                "utf-8", errors="replace") + completed.stderr.decode(
                "utf-8", errors="replace"))
    if completed.stdout or completed.stderr:
        append_event(journal, {
            "event": "launch_rejected",
            "label": label,
            "reason": "unexpected_terminal_output",
            "stdout": file_identity(stdout_path),
            "stderr": file_identity(stderr_path),
        })
        raise RuntimeError("benchmark emitted unexpected terminal output")
    if not output.is_file():
        append_event(journal, {
            "event": "launch_rejected",
            "label": label,
            "reason": "missing_retained_json",
        })
        raise RuntimeError("benchmark omitted retained JSON output")
    raw_document, raw_result_identity = stable_file(
        output, return_bytes=True, maximum_bytes=MAX_RESULT_BYTES,
        require_single_link=True)
    try:
        document = strict_json_bytes(raw_document)
        validator.validate_common(document, True)
        validator.validate_workload_digests(document)
    except BaseException as error:
        append_event(journal, {
            "event": "launch_rejected",
            "label": label,
            "reason": type(error).__name__,
            "message": str(error),
            "raw_result": raw_result_identity,
        })
        raise

    require(document["schema"] == "leopard2-benchmark-v33",
            "benchmark schema changed")
    parameters = document["parameters"]
    require((parameters["K"], parameters["R"],
             parameters["shard_bytes"]) == (k, r, shard_bytes),
            "workload identity changed")
    require(parameters["k65r65_t128_avx512_gfni_mode"] == mode,
            "requested diagnostic mode changed")
    require(parameters["requested_backend"] == "auto" and
            parameters["requested_profile"] == "legacy_high_v1" and
            parameters["requested_field"] == "gf8" and
            parameters["thread_count"] == 1 and
            parameters["batch"] == 1 and parameters["reuse"] == reuse and
            parameters["loss_count"] == 1 and
            parameters["seed"] == seed and
            parameters["iterations"] == ITERATIONS and
            parameters["warmup"] == WARMUP and
            parameters["measure_one_shot_encode"] is True and
            parameters["retain_samples"] is True and
            parameters["attest_source"] is True,
            "benchmark contract parameters changed")
    require(document["resolved"]["backend"] == "avx2",
            "AUTO did not resolve to the measured AVX2 control table")
    build = document["build"]
    require(build["source_commit"] == source_commit and
            build["source_tree"] == source_tree and
            build["source_tracked_dirty"] is False,
            "embedded source attestation changed")
    require(build["k65r65_t128_avx512_gfni_diagnostic_mode"] == mode and
            build["k65r65_t128_avx512_gfni_mode_latched"] == mode and
            build["k65r65_t128_avx512_gfni_kernel_qualified"] is True and
            build["k65r65_t128_avx512_gfni_selector_contract"] == CONTRACT,
            "runtime qualification metadata changed")
    selectable = role in {"target", "sustained"}
    selected = build["k65r65_t128_avx512_gfni_selector_selected"]
    expected_selected = bool(selectable and mode == 1)
    require(build[
                "k65r65_t128_avx512_gfni_selector_expected_selected"] is
            expected_selected and selected is expected_selected,
            "production selector did not match the exact-cell contract")
    require(build[
                "k65r65_t128_avx512_gfni_observed_call_count"] ==
            (2 if selected else 0),
            "untimed route-probe call count changed")
    require(build[
                "k65r65_t128_avx512_gfni_observed_vector_tile_count"] ==
            (2 * (shard_bytes // 64) if selected else 0),
            "untimed route-probe vector-tile count changed")
    require(document["correctness"]["leopard2_round_trip"] is True,
            "round-trip correctness failed")
    encode_us, encode_samples = timing(document, "encode_execution")
    one_shot_us, one_shot_samples = timing(document, "one_shot_encode")
    digests = document["workload_digests"]
    require(digests["algorithm"] == "fnv1a64",
            "workload digest algorithm changed")
    record = {
        "label": label,
        "mode": mode,
        "elapsed_ns": elapsed_ns,
        "reuse": reuse,
        "encode_us": encode_us,
        "one_shot_us": one_shot_us,
        "encode_samples_us": encode_samples,
        "one_shot_samples_us": one_shot_samples,
        "workload_digests": digests,
        "selector_selected": selected,
        "observed_call_count": build[
            "k65r65_t128_avx512_gfni_observed_call_count"],
        "observed_vector_tile_count": build[
            "k65r65_t128_avx512_gfni_observed_vector_tile_count"],
        "source": {
            "source_commit": build["source_commit"],
            "source_tree": build["source_tree"],
            "source_tracked_dirty": build["source_tracked_dirty"],
        },
        "raw": {
            "command": file_identity(
                command_path, require_single_link=True),
            "stdout": file_identity(
                stdout_path, require_single_link=True),
            "stderr": file_identity(
                stderr_path, require_single_link=True),
            "result": raw_result_identity,
        },
    }
    require(file_identity(output, require_single_link=True) ==
            raw_result_identity,
            "validated retained JSON changed before record publication")
    append_event(journal, {
        "event": "launch_validated",
        "label": label,
        "record": record,
    })
    return record


def summarize(log_values):
    count = len(log_values)
    require(count > 1, "not enough independent rounds")
    mean = statistics.mean(log_values)
    sample_sd = statistics.stdev(log_values)
    t_value = T_CRITICAL[count - 1]
    radius = t_value * sample_sd / math.sqrt(count)
    return {
        "rounds": count,
        "log_mean": mean,
        "log_sample_sd": sample_sd,
        "t_critical": t_value,
        "geometric_mean_speedup": math.exp(mean),
        "ci95": [math.exp(mean - radius), math.exp(mean + radius)],
    }


def run_cell(binary, cpu, sibling, cell, source_commit, source_tree, seed,
             validator, invocations, journal, deadline):
    cell_id, k, r, shard_bytes, rounds, role, region, _ = cell
    records = []
    expected_digests = None
    encode_logs = []
    one_shot_logs = []
    rejected_attempts = []
    for round_index in range(rounds):
        for attempt in range(MAX_ROUND_ATTEMPTS):
            order = (1, 0, 0, 1) if round_index % 2 == 0 \
                else (0, 1, 1, 0)
            launches = []
            before_cpu = cpu_snapshot(cpu)
            before_sibling = cpu_snapshot(sibling)
            try:
                for launch_index, mode in enumerate(order):
                    label = (
                        f"{cell_id}-round{round_index}-attempt{attempt}-"
                        f"slot{launch_index}-mode{mode}")
                    launch = run_launch(
                        binary, cpu, cell, mode, seed,
                        source_commit, source_tree, validator,
                        invocations, journal, label, deadline)
                    if expected_digests is None:
                        expected_digests = launch["workload_digests"]
                    require(launch["workload_digests"] == expected_digests,
                            f"{cell_id} digest changed across modes or launches")
                    launches.append(launch)
            except BaseException as error:
                append_event(journal, {
                    "event": "launch_attempt_failed",
                    "cell": cell_id,
                    "round": round_index,
                    "attempt": attempt,
                    "type": type(error).__name__,
                    "message": str(error),
                })
                raise
            isolation = {
                "benchmark_cpu": cpu_delta(
                    before_cpu, cpu_snapshot(cpu)),
                "reserved_sibling": cpu_delta(
                    before_sibling, cpu_snapshot(sibling)),
            }
            attempt_record = {
                "round": round_index,
                "attempt": attempt,
                "order": list(order),
                "launches": launches,
                "isolation": isolation,
            }
            isolation_passed = (
                isolation["benchmark_cpu"]["nonidle"] > 0 and
                isolation["reserved_sibling"]["total"] > 0 and
                isolation["reserved_sibling"]["idle"] > 0 and
                isolation["reserved_sibling"]["nonidle"] == 0)
            append_event(journal, {
                "event": "round_attempt_complete",
                "cell": cell_id,
                "round": round_index,
                "attempt": attempt,
                "isolation": isolation,
                "accepted": isolation_passed,
            })
            if isolation_passed:
                break
            rejected_attempts.append(attempt_record)
        else:
            raise RuntimeError(
                f"SMT sibling remained active in {cell_id} round "
                f"{round_index} for {MAX_ROUND_ATTEMPTS} attempts")
        controls = [launch for launch in launches if launch["mode"] == 0]
        candidates = [launch for launch in launches if launch["mode"] == 1]
        require(len(controls) == 2 and len(candidates) == 2,
                "round is not balanced")
        encode_log = statistics.mean(
            math.log(item["encode_us"]) for item in controls) - \
            statistics.mean(
                math.log(item["encode_us"]) for item in candidates)
        one_shot_log = statistics.mean(
            math.log(item["one_shot_us"]) for item in controls) - \
            statistics.mean(
                math.log(item["one_shot_us"]) for item in candidates)
        encode_logs.append(encode_log)
        one_shot_logs.append(one_shot_log)
        record = {
            "round": round_index,
            "attempt": attempt,
            "order": list(order),
            "log_control_over_candidate_encode": encode_log,
            "log_control_over_candidate_one_shot": one_shot_log,
            "launches": launches,
            "isolation": isolation,
        }
        records.append(record)
        append_event(journal, {
            "event": "round_complete",
            "cell": cell_id,
            "round": round_index,
            "encode_log": encode_log,
            "one_shot_log": one_shot_log,
        })
    return {
        "id": cell_id,
        "K": k,
        "R": r,
        "shard_bytes": shard_bytes,
        "role": role,
        "promotion_region": region,
        "rounds": rounds,
        "reuse": cell_reuse(cell),
        "workload_digests": expected_digests,
        "encode_execution": summarize(encode_logs),
        "one_shot_encode": summarize(one_shot_logs),
        "records": records,
        "rejected_contaminated_attempts": rejected_attempts,
    }


def base_report(options, binary_hash, started, results, provenance,
                lock_identity, quiet_presample):
    return {
        "schema": "leopard2-k65r65-t128-larger-avx512-gfni-frozen-abba/v1",
        "claim_scope": "same-binary current-control; no exact-main claim",
        "binary": str(options.binary),
        "binary_sha256_pre": binary_hash,
        "controller_command": list(sys.argv),
        "provenance": provenance,
        "canonical_lock": lock_identity,
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "source_tracked_dirty": False,
        "cpu": options.cpu,
        "sibling": options.sibling,
        "affinity": sorted(os.sched_getaffinity(0)),
        "iterations_per_launch": ITERATIONS,
        "warmup_per_launch": WARMUP,
        "child_environment": dict(CHILD_ENVIRONMENT),
        "git_environment": dict(GIT_ENVIRONMENT),
        "base_reuse_bytes_per_sample": BASE_REUSE_BYTES,
        "launches_per_round": 4,
        "max_round_attempts": MAX_ROUND_ATTEMPTS,
        "target_rounds": TARGET_ROUNDS,
        "neighbor_rounds": NEIGHBOR_ROUNDS,
        "sustained_rounds": SUSTAINED_ROUNDS,
        "sustained_reuse_multiplier": SUSTAINED_REUSE_MULTIPLIER,
        "started_unix_ns": started,
        "platform": platform.platform(),
        "quiet_presample": quiet_presample,
        "cells": results,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--build-binary", type=Path, required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--runner-sha256", required=True)
    parser.add_argument("--validator", type=Path, required=True)
    parser.add_argument("--validator-sha256", required=True)
    parser.add_argument("--source-archive", type=Path, required=True)
    parser.add_argument("--source-archive-sha256", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--sibling", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--journal", type=Path, required=True)
    parser.add_argument("--invocations", type=Path, required=True)
    parser.add_argument("--lock-fd", type=int)
    options = parser.parse_args()

    lock_descriptor, lock_identity = acquire_lock(options.lock_fd)
    require(set(os.sched_getaffinity(0)) == {options.cpu},
            "campaign controller must be singleton-pinned to the benchmark CPU")
    artifact_root = options.output.parent.resolve()
    require(artifact_root.is_dir(), "artifact directory is absent")
    for path in (
            options.binary, options.validator, options.source_archive,
            options.output,
            options.journal, options.invocations):
        require(path.is_absolute(), f"artifact path is not absolute: {path}")
        require(path.parent.resolve() == artifact_root,
                f"artifact escaped lane directory: {path}")
    require(options.build_binary.is_absolute(),
            "mutable build benchmark path is not absolute")
    require(options.repository.is_absolute(),
            "repository path is not absolute")
    runner_path = Path(__file__).absolute()
    artifact_paths = (
        options.binary, options.validator, options.source_archive,
        options.output, report_temporary_path(options.output),
        options.journal, options.invocations, runner_path,
    )
    canonical_artifact_paths = [path.resolve() for path in artifact_paths]
    require(len(set(canonical_artifact_paths)) ==
            len(canonical_artifact_paths),
            "campaign artifact paths or report temporary path alias")
    require(not options.output.exists(),
            "refusing to overwrite campaign report")
    require(not report_temporary_path(options.output).exists(),
            "refusing to reuse campaign report temporary path")
    require(not options.journal.exists(),
            "refusing to overwrite campaign journal")
    require(not options.invocations.exists(),
            "refusing to overwrite invocation evidence")
    options.invocations.mkdir(mode=0o700)

    require(runner_path.parent.resolve() == artifact_root,
            "qualification runner is not the frozen lane copy")
    frozen_binary = file_identity(
        options.binary, require_single_link=True, require_read_only=True)
    build_binary = file_identity(options.build_binary)
    runner_identity = file_identity(
        runner_path, require_single_link=True, require_read_only=True)
    validator_identity = file_identity(
        options.validator, require_single_link=True, require_read_only=True)
    source_archive_identity = file_identity(
        options.source_archive, require_single_link=True,
        require_read_only=True)
    require(frozen_binary["sha256"] == options.binary_sha256,
            "frozen benchmark pre-run SHA-256 mismatch")
    require(build_binary["sha256"] == frozen_binary["sha256"],
            "frozen and mutable benchmark hashes differ")
    require((frozen_binary["device"], frozen_binary["inode"]) !=
            (build_binary["device"], build_binary["inode"]),
            "frozen benchmark is not an independent byte copy")
    require(runner_identity["sha256"] == options.runner_sha256,
            "frozen runner SHA-256 mismatch")
    require(validator_identity["sha256"] == options.validator_sha256,
            "frozen validator SHA-256 mismatch")
    require(source_archive_identity["sha256"] ==
            options.source_archive_sha256,
            "frozen source archive SHA-256 mismatch")
    for value, label, length in (
            (options.binary_sha256, "binary SHA-256", 64),
            (options.runner_sha256, "runner SHA-256", 64),
            (options.validator_sha256, "validator SHA-256", 64),
            (options.source_archive_sha256, "source archive SHA-256", 64),
            (options.source_commit, "source commit", 40),
            (options.source_tree, "source tree", 40)):
        require(len(value) == length and
                all(character in "0123456789abcdef" for character in value),
                f"invalid {label}")
    validator = load_validator(options.validator, validator_identity)
    source_binding = verify_source_binding(
        options.repository, options.source_commit, options.source_tree,
        source_archive_identity, runner_identity, validator_identity)
    provenance = {
        "frozen_binary": frozen_binary,
        "mutable_build_binary_pre": build_binary,
        "runner_pre": runner_identity,
        "validator_pre": validator_identity,
        "source_archive_pre": source_archive_identity,
        "source_binding_pre": source_binding,
    }
    siblings_path = Path(
        f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
        "thread_siblings_list")
    observed_siblings = siblings_path.read_text(encoding="ascii").strip()
    require(observed_siblings == f"{options.cpu},{options.sibling}",
            "requested CPU/sibling topology changed")

    started = time.time_ns()
    deadline = time.monotonic() + CAMPAIGN_DEADLINE_SECONDS
    results = []
    append_event(options.journal, {
        "event": "campaign_started",
        "binary_sha256": frozen_binary["sha256"],
        "runner_sha256": runner_identity["sha256"],
        "validator_sha256": validator_identity["sha256"],
        "source_archive_sha256": source_archive_identity["sha256"],
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "controller_command": list(sys.argv),
        "canonical_lock": lock_identity,
        "started_unix_ns": started,
    })

    quiet_presample = {}

    def retain_failure(error):
        append_event(options.journal, {
            "event": "campaign_failed",
            "type": type(error).__name__,
            "message": str(error),
        })
        failure = base_report(
            options, frozen_binary["sha256"], started, results,
            provenance, lock_identity, quiet_presample)
        failure.update({
            "status": "failed",
            "claim_passed": False,
            "completed_unix_ns": time.time_ns(),
            "completed_cell_count": len(results),
            "total_cell_count": len(CELLS),
            "journal_at_failure": file_identity(options.journal),
            "failure": {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc(),
            },
        })
        write_report(options.output, failure)

    try:
        before_quiet_cpu = cpu_snapshot(options.cpu)
        before_quiet_sibling = cpu_snapshot(options.sibling)
        time.sleep(0.25)
        quiet_presample = {
            "benchmark_cpu": cpu_delta(
                before_quiet_cpu, cpu_snapshot(options.cpu)),
            "reserved_sibling": cpu_delta(
                before_quiet_sibling, cpu_snapshot(options.sibling)),
        }
        require(quiet_presample["reserved_sibling"]["total"] > 0 and
                quiet_presample["reserved_sibling"]["idle"] > 0 and
                quiet_presample["reserved_sibling"]["nonidle"] == 0,
                "reserved SMT sibling was active during quiet presample")
        append_event(options.journal, {
            "event": "quiet_presample",
            "isolation": quiet_presample,
        })
        for index, cell in enumerate(CELLS):
            require(time.monotonic() < deadline,
                    "campaign deadline expired before cell")
            verify_lock_continuity(
                lock_descriptor, lock_identity, options.cpu)
            require(sha256(options.binary) == frozen_binary["sha256"],
                    "frozen benchmark changed before cell")
            result = run_cell(
                options.binary, options.cpu, options.sibling, cell,
                options.source_commit, options.source_tree,
                0x6565B000 + index, validator,
                options.invocations, options.journal, deadline)
            results.append(result)
            append_event(options.journal, {
                "event": "cell_complete",
                "cell": result["id"],
                "encode_execution": result["encode_execution"],
                "one_shot_encode": result["one_shot_encode"],
            })
            checkpoint = base_report(
                options, frozen_binary["sha256"], started, results,
                provenance, lock_identity, quiet_presample)
            checkpoint.update({
                "status": "running",
                "completed_cell_count": len(results),
                "total_cell_count": len(CELLS),
                "journal_checkpoint": file_identity(options.journal),
            })
            write_report(options.output, checkpoint)
        require(len(results) == len(CELLS), "campaign cell count changed")
        lock_closure = verify_lock_continuity(
            lock_descriptor, lock_identity, options.cpu)
        source_binding_post = verify_source_binding(
            options.repository, options.source_commit, options.source_tree,
            source_archive_identity, runner_identity, validator_identity)
        require(source_binding_post == source_binding,
                "source binding changed during qualification")
        post_identities = {
            "frozen_binary": file_identity(
                options.binary, require_single_link=True,
                require_read_only=True),
            "mutable_build_binary": file_identity(options.build_binary),
            "runner": file_identity(
                runner_path, require_single_link=True, require_read_only=True),
            "validator": file_identity(
                options.validator, require_single_link=True,
                require_read_only=True),
            "source_archive": file_identity(
                options.source_archive, require_single_link=True,
                require_read_only=True),
            "source_binding": source_binding_post,
            "canonical_lock": lock_closure,
        }
        require(post_identities["frozen_binary"] == frozen_binary and
                post_identities["mutable_build_binary"] == build_binary and
                post_identities["runner"] == runner_identity and
                post_identities["validator"] == validator_identity and
                post_identities["source_archive"] == source_archive_identity,
                "campaign provenance changed during execution")

        target_gates = {}
        sustained_gates = {}
        neighbor_route_gates = {}
        for item in results:
            speed_gate = (
                item["encode_execution"]["ci95"][0] >= 1.05 and
                item["one_shot_encode"]["ci95"][0] >= 1.05)
            if item["role"] == "target":
                target_gates[item["id"]] = speed_gate
            elif item["role"] == "sustained":
                sustained_gates[item["id"]] = speed_gate
            else:
                route_inert = all(
                    launch["selector_selected"] is False and
                    launch["observed_call_count"] == 0 and
                    launch["observed_vector_tile_count"] == 0
                    for record in item["records"]
                    for launch in record["launches"])
                neighbor_route_gates[item["id"]] = route_inert
        region_gates = {
            "small": all(target_gates[
                f"target-k65-r65-b{shard_bytes}"]
                for shard_bytes in (128, 256, 512)) and
                sustained_gates["sustained-k65-r65-b128"],
            "large": all(target_gates[
                f"target-k65-r65-b{shard_bytes}"]
                for shard_bytes in (1024, 2048, 4096)) and
                sustained_gates["sustained-k65-r65-b4096"],
        }
        route_safety = all(neighbor_route_gates.values())
        claim_passed = all(region_gates.values()) and route_safety
        gate_results = {
            "target_speed": target_gates,
            "promotion_regions_including_sustained_endpoint": region_gates,
            "sustained_endpoint_speed": sustained_gates,
            "neighbor_route_inert": neighbor_route_gates,
            "global_route_safety": route_safety,
        }
        append_event(options.journal, {
            "event": "campaign_complete",
            "claim_passed": claim_passed,
            "gate_results": gate_results,
            "post_identities": post_identities,
        })
        journal_closure = file_identity(options.journal)
        report = base_report(
            options, frozen_binary["sha256"], started, results,
            provenance, lock_identity, quiet_presample)
        report.update({
            "status": "complete",
            "claim_passed": claim_passed,
            "binary_sha256_post": post_identities[
                "frozen_binary"]["sha256"],
            "post_identities": post_identities,
            "journal_closure": journal_closure,
            "observed_sibling_list": observed_siblings,
            "completed_unix_ns": time.time_ns(),
            "completed_cell_count": len(results),
            "total_cell_count": len(CELLS),
            "gate_policy": {
                "targets": (
                    "both per-cell 95% lower confidence bounds >= 1.05; "
                    "small and large regions promote all-or-none with their "
                    "sustained endpoint"),
                "sustained_endpoints": (
                    "both 95% lower confidence bounds >= 1.05 at 16x reuse"),
                "neighbors": (
                    "descriptive timing only; exact selector false, observed "
                    "call/tile counts zero, and identical retained digests"),
            },
            "gate_results": gate_results,
        })
        write_report(options.output, report)
    except BaseException as error:
        retain_failure(error)
        raise

    report_identity = file_identity(options.output)
    print(json.dumps({
        "claim_passed": claim_passed,
        "binary_sha256": frozen_binary["sha256"],
        "target_speedups": {
            item["id"]: {
                "encode": item["encode_execution"],
                "one_shot": item["one_shot_encode"],
            }
            for item in results if item["role"] == "target"
        },
        "sustained_speedups": {
            item["id"]: {
                "encode": item["encode_execution"],
                "one_shot": item["one_shot_encode"],
            }
            for item in results if item["role"] == "sustained"
        },
        "gate_results": report["gate_results"],
        "report": str(options.output),
        "report_identity": report_identity,
        "journal": str(options.journal),
        "journal_identity": journal_closure,
    }, indent=2, sort_keys=True))
    require(claim_passed, "campaign promotion gates failed")
    del lock_descriptor


if __name__ == "__main__":
    main()
