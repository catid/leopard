#!/usr/bin/env python3

"""Frozen same-binary qualification for AUTO GF16 GFNI encode routing.

The exact production cell is inferential.  Every inactive cell is a
deterministic selector-boundary control: its timings are retained for audit,
but neither its timing ratio nor any comparison between neighboring cells is a
promotion gate.
"""

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
INACTIVE_ROUNDS = 2
ITERATIONS = 9
WARMUP = 2
REUSE = 8
MIN_RETAINED_TIMER_WINDOW_US = 20_000.0
MAX_ROUND_ATTEMPTS = 5
CAMPAIGN_DEADLINE_SECONDS = 7200
MAX_RESULT_BYTES = 8 * 1024 * 1024
MAX_VALIDATOR_BYTES = 8 * 1024 * 1024
EXPECTED_CPU = 52
EXPECTED_SIBLING = 116
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
GIT = Path("/usr/bin/git")
TASKSET = Path("/usr/bin/taskset")
T_CRITICAL_24 = 2.0638985616280205
WORKLOAD_SEED = 0xA016F016
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
    "LEGACY_HIGH_V1,GF16,AUTO,resolved_AVX2,K=1000,R=200,T=256,"
    "B=65536,native_layout,native_cantor_affine,dense_full_parity,"
    "encode_only,runtime_GFNI,startup_KAT,calibrated_AMD_1A_08,"
    "one_shot_and_one_item_batch,codec_setup_descriptive_only"
)
FAILURE_REPORT_WRITTEN = False
JOURNAL_IDENTITIES = {}
ACTIVE_CELL_EVIDENCE = None

# id, K, R, bytes, rounds, role, profile, field, backend, rationale
CELLS = (
    ("target-k1000-r200-b65536-high-gf16-auto", 1000, 200, 65536,
     TARGET_ROUNDS,
     "target", "high", "gf16", "auto", "exact production selector cell"),
    ("inactive-k999", 999, 200, 65536, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "auto", "lower K boundary"),
    ("inactive-k1001", 1001, 200, 65536, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "auto", "upper K boundary"),
    ("inactive-r199", 1000, 199, 65536, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "auto", "lower R boundary"),
    ("inactive-r201", 1000, 201, 65536, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "auto", "upper R boundary"),
    ("inactive-b65534", 1000, 200, 65534, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "auto", "lower byte boundary"),
    ("inactive-b65538", 1000, 200, 65538, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "auto", "upper byte boundary"),
    ("inactive-b32768", 1000, 200, 32768, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "auto", "smaller byte control"),
    ("inactive-b131072", 1000, 200, 131072, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "auto", "larger byte control"),
    ("inactive-low-gf16", 1000, 200, 65536, INACTIVE_ROUNDS,
     "inactive", "low", "gf16", "auto", "profile boundary"),
    ("inactive-explicit-avx2", 1000, 200, 65536, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "avx2", "AVX2 backend boundary"),
    ("inactive-explicit-gfni", 1000, 200, 65536, INACTIVE_ROUNDS,
     "inactive", "high", "gf16", "gfni", "GFNI backend boundary"),
    ("inactive-high-gf8", 128, 64, 65536, INACTIVE_ROUNDS,
     "inactive", "high", "gf8", "auto", "field boundary"),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def expected_route(k, r, shard_bytes, profile, field, backend, mode):
    require(mode in (0, 1), "diagnostic mode is not binary")
    context_backend = "avx2-gfni" if backend == "gfni" else "avx2"
    exact_codec_identity = (
        k == 1000 and r == 200 and profile == "high" and
        field == "gf16" and backend == "auto")
    available = bool(mode == 1 and exact_codec_identity)
    selected = bool(available and shard_bytes == 65536)
    return {
        "kernel_available": available,
        "kernel_qualified": available,
        "selector_selected": selected,
        "observed_call_count": 2 if selected else 0,
        "context_backend": context_backend,
        "encode_backend": "avx2-gfni" if selected else context_backend,
    }


def verify_static_contract():
    require((TARGET_ROUNDS, INACTIVE_ROUNDS, ITERATIONS, WARMUP, REUSE,
             MIN_RETAINED_TIMER_WINDOW_US) ==
            (25, 2, 9, 2, 8, 20_000.0),
            "frozen campaign counts changed")
    require((EXPECTED_CPU, EXPECTED_SIBLING) == (52, 116),
            "frozen CPU topology changed")
    require(WORKLOAD_SEED == 0xA016F016,
            "frozen workload seed changed")
    require(CONTRACT ==
            "LEGACY_HIGH_V1,GF16,AUTO,resolved_AVX2,K=1000,R=200,T=256,"
            "B=65536,native_layout,native_cantor_affine,dense_full_parity,"
            "encode_only,runtime_GFNI,startup_KAT,calibrated_AMD_1A_08,"
            "one_shot_and_one_item_batch,codec_setup_descriptive_only",
            "frozen selector contract changed")
    require(len(CELLS) == 13 and
            len({cell[0] for cell in CELLS}) == len(CELLS),
            "campaign cell census changed or contains duplicate IDs")
    targets = [cell for cell in CELLS if cell[5] == "target"]
    require(len(targets) == 1 and targets[0][1:9] ==
            (1000, 200, 65536, 25, "target", "high", "gf16", "auto"),
            "exact target cell changed")
    observed_inactive = {
        (cell[0], cell[1], cell[2], cell[3], cell[6], cell[7], cell[8])
        for cell in CELLS if cell[5] == "inactive"
    }
    expected_inactive = {
        ("inactive-k999", 999, 200, 65536, "high", "gf16", "auto"),
        ("inactive-k1001", 1001, 200, 65536, "high", "gf16", "auto"),
        ("inactive-r199", 1000, 199, 65536, "high", "gf16", "auto"),
        ("inactive-r201", 1000, 201, 65536, "high", "gf16", "auto"),
        ("inactive-b65534", 1000, 200, 65534, "high", "gf16", "auto"),
        ("inactive-b65538", 1000, 200, 65538, "high", "gf16", "auto"),
        ("inactive-b32768", 1000, 200, 32768, "high", "gf16", "auto"),
        ("inactive-b131072", 1000, 200, 131072, "high", "gf16", "auto"),
        ("inactive-low-gf16", 1000, 200, 65536, "low", "gf16", "auto"),
        ("inactive-explicit-avx2", 1000, 200, 65536,
         "high", "gf16", "avx2"),
        ("inactive-explicit-gfni", 1000, 200, 65536,
         "high", "gf16", "gfni"),
        ("inactive-high-gf8", 128, 64, 65536, "high", "gf8", "auto"),
    }
    require(observed_inactive == expected_inactive and
            all(cell[4] == INACTIVE_ROUNDS
                for cell in CELLS if cell[5] == "inactive"),
            "inactive selector-boundary cell census changed")


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


def retained_invocation_evidence(command_path, stdout_path, stderr_path,
                                 result_path):
    evidence = {
        "command": file_identity(command_path, require_single_link=True),
        "stdout": file_identity(stdout_path, require_single_link=True),
        "stderr": file_identity(stderr_path, require_single_link=True),
    }
    if result_path.exists():
        evidence["result"] = stable_file(
            result_path, maximum_bytes=MAX_RESULT_BYTES,
            require_single_link=True)
    else:
        evidence["result"] = None
    return evidence


def verify_retained_raw_evidence(results):
    identities = []
    accepted_launches = 0
    rejected_contaminated_launches = 0
    for cell in results:
        groups = (
            (cell["records"], "accepted"),
            (cell["rejected_contaminated_attempts"], "rejected"),
        )
        for attempts, disposition in groups:
            for attempt in attempts:
                for launch in attempt["launches"]:
                    if disposition == "accepted":
                        accepted_launches += 1
                    else:
                        rejected_contaminated_launches += 1
                    for name in ("command", "stdout", "stderr", "result"):
                        expected = launch["raw"][name]
                        observed = file_identity(
                            Path(expected["path"]), require_single_link=True)
                        require(observed == expected,
                                "retained raw evidence changed: " +
                                expected["path"])
                        identities.append({
                            "launch": launch["label"],
                            "disposition": disposition,
                            "kind": name,
                            "path": expected["path"],
                            "sha256": expected["sha256"],
                            "size": expected["size"],
                        })
    paths = [item["path"] for item in identities]
    require(len(paths) == len(set(paths)),
            "retained raw evidence paths are not unique")
    encoded = json.dumps(
        identities, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return {
        "accepted_launches": accepted_launches,
        "rejected_contaminated_launches": rejected_contaminated_launches,
        "file_count": len(identities),
        "manifest_sha256": hashlib.sha256(encoded).hexdigest(),
        "manifest_semantics": (
            "SHA-256 of canonical JSON records ordered by cell, disposition, "
            "attempt, launch, and command/stdout/stderr/result"),
    }


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
        lock_scope = "campaign-only"
    else:
        descriptor = inherited_fd
        lock_scope = "wrapper-build-copy-campaign"
    try:
        canonical = LOCK_PATH.lstat()
        opened = os.fstat(descriptor)
        require(stat.S_ISREG(canonical.st_mode) and
                stat.S_ISREG(opened.st_mode) and
                canonical.st_nlink == 1 and opened.st_nlink == 1,
                "canonical benchmark lock is not one regular file")
        require((opened.st_dev, opened.st_ino) ==
                (canonical.st_dev, canonical.st_ino),
                "lock descriptor is not the canonical lock")
        if inherited_fd is not None:
            lock_mode = "inherited-across-build-copy-campaign"
    except BaseException:
        if inherited_fd is None:
            os.close(descriptor)
        raise
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as error:
        if inherited_fd is None:
            os.close(descriptor)
        raise RuntimeError(
            f"canonical benchmark lock is busy: {LOCK_PATH}") from error
    try:
        locked = os.fstat(descriptor)
        canonical_after = LOCK_PATH.lstat()
        require(stat.S_ISREG(canonical_after.st_mode) and
                canonical_after.st_nlink == 1 and locked.st_nlink == 1 and
                (locked.st_dev, locked.st_ino) ==
                (canonical.st_dev, canonical.st_ino) ==
                (canonical_after.st_dev, canonical_after.st_ino),
                "canonical benchmark lock changed during acquisition")
    except BaseException:
        if inherited_fd is None:
            os.close(descriptor)
        raise
    return descriptor, {
        "path": str(LOCK_PATH),
        "mode": lock_mode,
        "scope": lock_scope,
        "descriptor": descriptor,
        "device": locked.st_dev,
        "inode": locked.st_ino,
    }


def verify_lock_continuity(descriptor, identity, cpu):
    canonical = LOCK_PATH.lstat()
    held = os.fstat(descriptor)
    require(stat.S_ISREG(canonical.st_mode) and
            stat.S_ISREG(held.st_mode) and
            canonical.st_nlink == 1 and held.st_nlink == 1,
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
        ("experiments/leopard2/gfni_codec/"
         "run_auto_gf16_gfni_encode_abba.py", runner_identity),
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
    retained_windows = [sample * REUSE for sample in samples]
    require(min(retained_windows) >= MIN_RETAINED_TIMER_WINDOW_US,
            f"{name} retained timer window is shorter than the fixed floor")
    return reported, samples, min(retained_windows)


def append_event(path, event):
    encoded = json.dumps(event, sort_keys=True, separators=(",", ":")) + "\n"
    key = str(path)
    expected = JOURNAL_IDENTITIES.get(key)
    flags = os.O_WRONLY | os.O_APPEND | os.O_CLOEXEC | os.O_NOFOLLOW
    if expected is None:
        flags |= os.O_CREAT | os.O_EXCL
    descriptor = os.open(
        path, flags, 0o600)
    try:
        metadata = os.fstat(descriptor)
        require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
                "campaign journal is not one regular file")
        if expected is None:
            require(metadata.st_size == 0,
                    "new campaign journal was not empty")
        else:
            require((metadata.st_dev, metadata.st_ino, metadata.st_size) ==
                    (expected["device"], expected["inode"], expected["size"]),
                    "campaign journal continuity was lost")
        data = encoded.encode("utf-8")
        written = 0
        while written < len(data):
            count = os.write(descriptor, data[written:])
            require(count > 0, "campaign journal write made no progress")
            written += count
        os.fsync(descriptor)
        after = os.fstat(descriptor)
        path_after = path.lstat()
        require(stat.S_ISREG(path_after.st_mode) and
                path_after.st_nlink == 1 and
                (after.st_dev, after.st_ino) ==
                (path_after.st_dev, path_after.st_ino) and
                after.st_size == metadata.st_size + len(data),
                "campaign journal changed during append")
        JOURNAL_IDENTITIES[key] = {
            "device": after.st_dev,
            "inode": after.st_ino,
            "size": after.st_size,
        }
    finally:
        os.close(descriptor)
    if expected is None:
        directory = os.open(path.parent, os.O_RDONLY | os.O_CLOEXEC)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)


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


def emergency_output_from_argv():
    values = [sys.argv[index + 1] for index, value in enumerate(sys.argv[:-1])
              if value == "--output"]
    if len(values) != 1:
        return None
    path = Path(values[0])
    if not path.is_absolute() or path.exists() or \
            report_temporary_path(path).exists():
        return None
    try:
        parent = path.parent.lstat()
    except OSError:
        return None
    if not stat.S_ISDIR(parent.st_mode) or path.parent.is_symlink():
        return None
    return path


def retain_bootstrap_failure(error):
    global FAILURE_REPORT_WRITTEN
    if FAILURE_REPORT_WRITTEN:
        return
    output = emergency_output_from_argv()
    if output is None:
        return
    write_report(output, {
        "schema": "leopard2-auto-gf16-gfni-encode-frozen-abba/v1",
        "status": "failed",
        "claim_passed": False,
        "failure_phase": "bootstrap_or_preflight",
        "controller_command": list(sys.argv),
        "completed_unix_ns": time.time_ns(),
        "failure": {
            "type": type(error).__name__,
            "message": str(error),
            "traceback": traceback.format_exc(),
        },
    })
    FAILURE_REPORT_WRITTEN = True


def run_launch(binary, cpu, cell, mode, seed, source_commit, source_tree,
               validator, invocations, journal, label, deadline):
    (cell_id, k, r, shard_bytes, _, role, profile, field, backend,
     _) = cell
    invocation = invocations / label
    invocation.mkdir(mode=0o700)
    output = invocation / "result.json"
    stdout_path = invocation / "stdout"
    stderr_path = invocation / "stderr"
    command_path = invocation / "command.json"
    command = [
        str(TASKSET), "-c", str(cpu), str(binary),
        "--k", str(k), "--r", str(r),
        "--profile", profile, "--field", field,
        "--backend", backend, "--bytes", str(shard_bytes),
        "--loss", "1", "--batch", "1", "--reuse", str(REUSE),
        "--iterations", str(ITERATIONS), "--warmup", str(WARMUP),
        "--threads", "1", "--seed", str(seed),
        "--skip-legacy", "--retain-samples",
        "--measure-one-shot-encode", "--attest-source",
        "--auto-gf16-gfni-encode-mode", str(mode),
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
            "raw": retained_invocation_evidence(
                command_path, stdout_path, stderr_path, output),
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
            "raw": retained_invocation_evidence(
                command_path, stdout_path, stderr_path, output),
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
            "raw": retained_invocation_evidence(
                command_path, stdout_path, stderr_path, output),
        })
        raise RuntimeError("benchmark emitted unexpected terminal output")
    if not output.is_file():
        append_event(journal, {
            "event": "launch_rejected",
            "label": label,
            "reason": "missing_retained_json",
            "raw": retained_invocation_evidence(
                command_path, stdout_path, stderr_path, output),
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
            "raw": retained_invocation_evidence(
                command_path, stdout_path, stderr_path, output),
        })
        raise

    require(document["schema"] == "leopard2-benchmark-v35",
            "benchmark schema changed")
    parameters = document["parameters"]
    require((parameters["K"], parameters["R"],
             parameters["shard_bytes"]) == (k, r, shard_bytes),
            "workload identity changed")
    require(parameters["auto_gf16_gfni_encode_mode"] == mode,
            "requested diagnostic mode changed")
    expected_profile = "legacy_high_v1" if profile == "high" else "low_v1"
    expected_requested_backend = (
        "avx2-gfni" if backend == "gfni" else backend)
    require(parameters["requested_backend"] == expected_requested_backend and
            parameters["requested_profile"] == expected_profile and
            parameters["requested_field"] == field and
            parameters["skip_legacy"] is True and
            parameters["thread_count"] == 1 and
            parameters["batch"] == 1 and parameters["reuse"] == REUSE and
            parameters["loss_count"] == 1 and
            parameters["seed"] == seed and
            parameters["iterations"] == ITERATIONS and
            parameters["warmup"] == WARMUP and
            parameters["measure_one_shot_encode"] is True and
            parameters["retain_samples"] is True and
            parameters["attest_source"] is True,
            "benchmark contract parameters changed")
    route = expected_route(
        k, r, shard_bytes, profile, field, backend, mode)
    require(document["resolved"]["backend"] == route["context_backend"] and
            document["resolved"]["profile"] == expected_profile and
            document["resolved"]["field"] == field,
            "workload did not resolve to the expected context table")
    build = document["build"]
    require(build["source_commit"] == source_commit and
            build["source_tree"] == source_tree and
            build["source_tracked_dirty"] is False,
            "embedded source attestation changed")
    require(build["auto_gf16_gfni_encode_diagnostic_mode"] == mode and
            build["auto_gf16_gfni_encode_diagnostic_disabled"] is
            (mode == 0) and
            build["auto_gf16_gfni_encode_mode_latched"] == mode and
            build["auto_gf16_gfni_encode_kernel_available"] is
            route["kernel_available"] and
            build["auto_gf16_gfni_encode_kernel_qualified"] is
            route["kernel_qualified"] and
            build["auto_gf16_gfni_encode_selector_contract"] == CONTRACT and
            build["auto_gf16_gfni_encode_timed_ordinary_encode_api"] ==
            "leo2_encode_batch:item_count=1:no_preflight_scratch" and
            build["auto_gf16_gfni_encode_timed_one_shot_encode_api"] ==
            "leo2_encode",
            "runtime qualification metadata changed")
    selected = build["auto_gf16_gfni_encode_selector_selected"]
    require(build[
                "auto_gf16_gfni_encode_selector_expected_selected"] is
            route["selector_selected"] and
            selected is route["selector_selected"],
            "production selector did not match the exact-cell contract")
    require(build[
                "auto_gf16_gfni_encode_observed_call_count"] ==
            route["observed_call_count"],
            "untimed route-probe call count changed")
    require(document["resolved"]["encode_backend"] ==
            route["encode_backend"],
            "resolved encode backend differs from the selected route")
    require(document["correctness"]["leopard2_round_trip"] is True,
            "round-trip correctness failed")
    encode_us, encode_samples, encode_min_window_us = timing(
        document, "encode_execution")
    one_shot_us, one_shot_samples, one_shot_min_window_us = timing(
        document, "one_shot_encode")
    digests = document["workload_digests"]
    require(digests["algorithm"] == "fnv1a64",
            "workload digest algorithm changed")
    record = {
        "label": label,
        "mode": mode,
        "elapsed_ns": elapsed_ns,
        "reuse": REUSE,
        "encode_us": encode_us,
        "one_shot_us": one_shot_us,
        "encode_samples_us": encode_samples,
        "one_shot_samples_us": one_shot_samples,
        "retained_timer_window_us": {
            "required_floor": MIN_RETAINED_TIMER_WINDOW_US,
            "encode_execution_min": encode_min_window_us,
            "one_shot_encode_min": one_shot_min_window_us,
        },
        "workload_digests": digests,
        "resolved": document["resolved"],
        "expected_route": route,
        "kernel_available": build[
            "auto_gf16_gfni_encode_kernel_available"],
        "kernel_qualified": build[
            "auto_gf16_gfni_encode_kernel_qualified"],
        "selector_selected": selected,
        "observed_call_count": build[
            "auto_gf16_gfni_encode_observed_call_count"],
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


def summarize(log_values, inferential):
    count = len(log_values)
    require(count > 0, "not enough completed rounds")
    mean = statistics.mean(log_values)
    result = {
        "rounds": count,
        "log_mean": mean,
        "geometric_mean_speedup": math.exp(mean),
        "inferential": inferential,
    }
    if inferential:
        require(count == TARGET_ROUNDS,
                "inferential target round count changed")
        sample_sd = statistics.stdev(log_values)
        radius = T_CRITICAL_24 * sample_sd / math.sqrt(count)
        result.update({
            "log_sample_sd": sample_sd,
            "t_critical": T_CRITICAL_24,
            "ci95": [math.exp(mean - radius), math.exp(mean + radius)],
        })
    else:
        result.update({
            "log_sample_sd": statistics.stdev(log_values)
            if count > 1 else None,
            "t_critical": None,
            "ci95": None,
        })
    return result


def run_cell(binary, cpu, sibling, cell, source_commit, source_tree, seed,
             validator, invocations, journal, deadline):
    global ACTIVE_CELL_EVIDENCE
    (cell_id, k, r, shard_bytes, rounds, role, profile, field, backend,
     rationale) = cell
    records = []
    expected_digests = None
    encode_logs = []
    one_shot_logs = []
    rejected_attempts = []
    ACTIVE_CELL_EVIDENCE = {
        "id": cell_id,
        "completed_rounds": 0,
        "rejected_contaminated_attempts": rejected_attempts,
    }
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
                    "label": label,
                    "invocation": str(invocations / label),
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
        ACTIVE_CELL_EVIDENCE["completed_rounds"] = len(records)
        append_event(journal, {
            "event": "round_complete",
            "cell": cell_id,
            "round": round_index,
            "encode_log": encode_log,
            "one_shot_log": one_shot_log,
        })
    result = {
        "id": cell_id,
        "K": k,
        "R": r,
        "shard_bytes": shard_bytes,
        "role": role,
        "profile": profile,
        "field": field,
        "backend": backend,
        "rationale": rationale,
        "rounds": rounds,
        "reuse": REUSE,
        "workload_digests": expected_digests,
        "encode_execution": summarize(encode_logs, role == "target"),
        "one_shot_encode": summarize(one_shot_logs, role == "target"),
        "records": records,
        "rejected_contaminated_attempts": rejected_attempts,
    }
    ACTIVE_CELL_EVIDENCE = None
    return result


def base_report(options, binary_hash, started, results, provenance,
                lock_identity, quiet_presample):
    return {
        "schema": "leopard2-auto-gf16-gfni-encode-frozen-abba/v1",
        "claim_scope": (
            "same-binary mode-0 control versus mode-1 candidate at the exact "
            "K1000/R200/B65536/high/GF16/AUTO cell; no exact-main claim"),
        "binary": str(options.binary),
        "binary_sha256_pre": binary_hash,
        "controller_command": list(sys.argv),
        "provenance": provenance,
        "canonical_lock": lock_identity,
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "source_tracked_dirty": False,
        "workload_seed": WORKLOAD_SEED,
        "cpu": options.cpu,
        "sibling": options.sibling,
        "affinity": sorted(os.sched_getaffinity(0)),
        "iterations_per_launch": ITERATIONS,
        "warmup_per_launch": WARMUP,
        "child_environment": dict(CHILD_ENVIRONMENT),
        "git_environment": dict(GIT_ENVIRONMENT),
        "reuse_per_sample": REUSE,
        "retained_timer_duration_floor_us":
            MIN_RETAINED_TIMER_WINDOW_US,
        "launches_per_round": 4,
        "max_round_attempts": MAX_ROUND_ATTEMPTS,
        "target_rounds": TARGET_ROUNDS,
        "inactive_control_rounds": INACTIVE_ROUNDS,
        "benchmark_schema_required": "leopard2-benchmark-v35",
        "diagnostic_option": "--auto-gf16-gfni-encode-mode",
        "co_primary_gate": {
            "ordinary_encode_execution_lcb95_min": 1.05,
            "one_shot_encode_lcb95_min": 1.05,
        },
        "started_unix_ns": started,
        "platform": platform.platform(),
        "quiet_presample": quiet_presample,
        "active_incomplete_cell": ACTIVE_CELL_EVIDENCE,
        "cells": results,
    }


def self_test():
    verify_static_contract()

    duplicate_rejected = False
    try:
        strict_json_bytes(b'{"value":1,"value":2}')
    except RuntimeError:
        duplicate_rejected = True
    require(duplicate_rejected, "duplicate-key self-test failed")

    nonfinite_rejected = False
    try:
        strict_json_bytes(b'{"value":NaN}')
    except RuntimeError:
        nonfinite_rejected = True
    require(nonfinite_rejected, "non-finite self-test failed")

    target_off = expected_route(
        1000, 200, 65536, "high", "gf16", "auto", 0)
    target_on = expected_route(
        1000, 200, 65536, "high", "gf16", "auto", 1)
    byte_neighbor_on = expected_route(
        1000, 200, 65534, "high", "gf16", "auto", 1)
    explicit_on = expected_route(
        1000, 200, 65536, "high", "gf16", "gfni", 1)
    require(not target_off["kernel_available"] and
            not target_off["selector_selected"] and
            target_on["kernel_available"] and
            target_on["selector_selected"] and
            target_on["observed_call_count"] == 2 and
            target_on["encode_backend"] == "avx2-gfni" and
            byte_neighbor_on["kernel_available"] and
            not byte_neighbor_on["selector_selected"] and
            byte_neighbor_on["observed_call_count"] == 0 and
            not explicit_on["kernel_available"] and
            explicit_on["context_backend"] == "avx2-gfni" and
            explicit_on["encode_backend"] == "avx2-gfni",
            "route-contract self-test failed")

    retained = [2500.0] * ITERATIONS
    document = {"metrics": {"probe": {
        "samples_us_per_batch_call": retained,
        "median_us_per_batch_call": 2500.0,
    }}}
    median, samples, minimum = timing(document, "probe")
    require(median == 2500.0 and samples == retained and
            minimum == MIN_RETAINED_TIMER_WINDOW_US,
            "timer-floor self-test failed")
    summary = summarize([0.0] * TARGET_ROUNDS, True)
    require(summary["geometric_mean_speedup"] == 1.0 and
            summary["ci95"] == [1.0, 1.0],
            "paired-log summary self-test failed")

    print(json.dumps({
        "schema": "leopard2-auto-gf16-gfni-encode-runner-self-test/v1",
        "passed": True,
        "cells": len(CELLS),
        "target_rounds": TARGET_ROUNDS,
        "inactive_rounds": INACTIVE_ROUNDS,
        "iterations": ITERATIONS,
        "warmup": WARMUP,
        "reuse": REUSE,
        "timer_floor_us": MIN_RETAINED_TIMER_WINDOW_US,
    }, sort_keys=True))


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

    verify_static_contract()
    require(options.cpu == EXPECTED_CPU and
            options.sibling == EXPECTED_SIBLING,
            "this frozen campaign is fixed to CPU 52 and sibling 116")
    lock_descriptor, lock_identity = acquire_lock(options.lock_fd)
    require(set(os.sched_getaffinity(0)) == {options.cpu},
            "campaign controller must be singleton-pinned to the benchmark CPU")
    artifact_root = options.output.parent.resolve()
    artifact_parent = options.output.parent.lstat()
    require(stat.S_ISDIR(artifact_parent.st_mode) and
            not options.output.parent.is_symlink() and
            artifact_root.is_dir(),
            "artifact directory is absent, aliased, or symlinked")
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
    taskset_identity = file_identity(TASKSET)
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
        "taskset_pre": taskset_identity,
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
        global FAILURE_REPORT_WRITTEN
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
        FAILURE_REPORT_WRITTEN = True

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
        for cell in CELLS:
            require(time.monotonic() < deadline,
                    "campaign deadline expired before cell")
            verify_lock_continuity(
                lock_descriptor, lock_identity, options.cpu)
            require(sha256(options.binary) == frozen_binary["sha256"],
                    "frozen benchmark changed before cell")
            result = run_cell(
                options.binary, options.cpu, options.sibling, cell,
                options.source_commit, options.source_tree,
                WORKLOAD_SEED, validator,
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
        raw_evidence_closure = verify_retained_raw_evidence(results)
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
            "taskset": file_identity(TASKSET),
            "source_binding": source_binding_post,
            "canonical_lock": lock_closure,
            "retained_raw_evidence": raw_evidence_closure,
        }
        require(post_identities["frozen_binary"] == frozen_binary and
                post_identities["mutable_build_binary"] == build_binary and
                post_identities["runner"] == runner_identity and
                post_identities["validator"] == validator_identity and
                post_identities["source_archive"] == source_archive_identity and
                post_identities["taskset"] == taskset_identity,
                "campaign provenance changed during execution")

        targets = [item for item in results if item["role"] == "target"]
        require(len(targets) == 1, "exactly one inferential target is required")
        target = targets[0]
        ordinary_speed = target["encode_execution"]["ci95"][0] >= 1.05
        one_shot_speed = target["one_shot_encode"]["ci95"][0] >= 1.05
        route_gates = {
            item["id"]: all(
                launch["expected_route"] == expected_route(
                    item["K"], item["R"], item["shard_bytes"],
                    item["profile"], item["field"], item["backend"],
                    launch["mode"]) and
                launch["kernel_available"] is
                    launch["expected_route"]["kernel_available"] and
                launch["kernel_qualified"] is
                    launch["expected_route"]["kernel_qualified"] and
                launch["selector_selected"] is
                    launch["expected_route"]["selector_selected"] and
                launch["observed_call_count"] ==
                    launch["expected_route"]["observed_call_count"] and
                launch["resolved"]["backend"] ==
                    launch["expected_route"]["context_backend"] and
                launch["resolved"]["encode_backend"] ==
                    launch["expected_route"]["encode_backend"] and
                launch["resolved"]["profile"] ==
                    ("legacy_high_v1" if item["profile"] == "high"
                     else "low_v1") and
                launch["resolved"]["field"] == item["field"]
                for record in (
                    item["records"] +
                    item["rejected_contaminated_attempts"])
                for launch in record["launches"])
            for item in results
        }
        digest_gates = {
            item["id"]: all(
                launch["workload_digests"] == item["workload_digests"]
                for record in (
                    item["records"] +
                    item["rejected_contaminated_attempts"])
                for launch in record["launches"])
            for item in results
        }
        result_by_id = {item["id"]: item for item in results}
        backend_digest_ids = (
            "target-k1000-r200-b65536-high-gf16-auto",
            "inactive-explicit-avx2",
            "inactive-explicit-gfni",
        )
        backend_digest_reference = result_by_id[
            backend_digest_ids[0]]["workload_digests"]
        cross_backend_digest_identity = all(
            result_by_id[cell_id]["workload_digests"] ==
            backend_digest_reference
            for cell_id in backend_digest_ids)
        timer_window_gates = {
            item["id"]: all(
                launch["retained_timer_window_us"][
                    "encode_execution_min"] >=
                MIN_RETAINED_TIMER_WINDOW_US and
                launch["retained_timer_window_us"][
                    "one_shot_encode_min"] >=
                MIN_RETAINED_TIMER_WINDOW_US
                for attempt in (
                    item["records"] +
                    item["rejected_contaminated_attempts"])
                for launch in attempt["launches"])
            for item in results
        }
        route_safety = all(route_gates.values())
        digest_safety = (
            all(digest_gates.values()) and cross_backend_digest_identity)
        timer_window_safety = all(timer_window_gates.values())
        claim_passed = (
            ordinary_speed and one_shot_speed and route_safety and
            digest_safety and timer_window_safety)
        gate_results = {
            "target_ordinary_encode_lcb95_at_least_1_05": ordinary_speed,
            "target_one_shot_encode_lcb95_at_least_1_05": one_shot_speed,
            "per_cell_availability_selection_count_and_resolved_route":
                route_gates,
            "per_cell_cross_mode_digest_identity": digest_gates,
            "target_explicit_backend_digest_identity":
                cross_backend_digest_identity,
            "per_cell_retained_timer_duration_floor": timer_window_gates,
            "global_route_safety": route_safety,
            "global_digest_safety": digest_safety,
            "global_retained_timer_duration_safety": timer_window_safety,
            "inactive_timing_used_as_gate": False,
        }
        append_event(options.journal, {
            "event": "campaign_complete",
            "claim_passed": claim_passed,
            "gate_results": gate_results,
            "post_identities": post_identities,
            "retained_raw_evidence": raw_evidence_closure,
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
                "target": (
                    "ordinary encode_execution and one_shot_encode are "
                    "co-primary; both 25-round paired-log t-interval 95% "
                    "lower bounds must be >= 1.05"),
                "route": (
                    "every launch must match the frozen mode-specific "
                    "availability, qualification, selection, call-count, "
                    "context-backend, and encode-backend contract; only the "
                    "target selects in mode 1, while exact-K/R byte controls "
                    "may qualify without selecting"),
                "digests": (
                    "all retained workload digests must be identical across "
                    "modes, order slots, attempts, and rounds within a cell; "
                    "the target and both explicit-backend controls must also "
                    "match under the common workload seed"),
                "timer_duration": (
                    "every retained encode_execution and one_shot_encode "
                    "sample multiplied by reuse must span at least 20 "
                    "milliseconds"),
                "inactive_controls": (
                    "timing is retained as descriptive evidence only and is "
                    "never a promotion gate"),
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
        "inactive_descriptive_speedups": {
            item["id"]: {
                "encode": item["encode_execution"],
                "one_shot": item["one_shot_encode"],
            }
            for item in results if item["role"] == "inactive"
        },
        "gate_results": report["gate_results"],
        "report": str(options.output),
        "report_identity": report_identity,
        "journal": str(options.journal),
        "journal_identity": journal_closure,
    }, indent=2, sort_keys=True))
    require(claim_passed, "campaign promotion gates failed")
    os.close(lock_descriptor)


if __name__ == "__main__":
    try:
        if sys.argv[1:] == ["--self-test"]:
            self_test()
        else:
            main()
    except BaseException as error:
        retain_bootstrap_failure(error)
        raise
