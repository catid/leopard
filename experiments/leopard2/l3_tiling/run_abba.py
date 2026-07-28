#!/usr/bin/env python3
"""Run fail-closed GF16 AVX2 L3-tiling ABBA evidence.

Two caller-named lanes are compared from frozen, campaign-owned benchmark
executables.  Each lane must either identify an existing benchmark from an
exact clean Git worktree, or request a fresh Release build from such a
worktree.  The runner:

* serializes build, validation, and timing with the canonical Leopard lock;
* verifies the benchmark's compiled-in source commit and tree;
* copies it into a read-only lane artifact and verifies that artifact before
  and after every child process;
* reserves one of the two measured physical CPU pairs;
* rejects any round with activity on the reserved SMT sibling; and
* retains all nine child samples for three independent ABBA rounds.

This program intentionally makes no online policy decision.  It produces raw
evidence and a descriptive comparison for an offline promotion decision.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import platform
import re
import shutil
import stat
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


CANONICAL_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")
ALLOWED_CPU_PAIRS = frozenset(((4, 20), (14, 30)))
RAW_SCHEMA = "leopard2-l3-tiling-abba-raw/v1"
SUMMARY_SCHEMA = "leopard2-l3-tiling-abba-summary/v1"
MANIFEST_SCHEMA = "leopard2-l3-tiling-abba-manifest/v1"
FAILURE_SCHEMA = "leopard2-l3-tiling-abba-failure/v1"
ATTESTATION_SCHEMA = "leopard2-l3-tiling-lane-attestation/v1"
ARTIFACT_SCHEMA = "leopard2-l3-tiling-frozen-artifact/v1"
ORDER_SIDES = (
    ("left", "right", "right", "left"),
    ("right", "left", "left", "right"),
    ("left", "right", "right", "left"),
)
ROUNDS = len(ORDER_SIDES)
RETAINED_SAMPLES = 9
WARMUP = 2
REUSE = 1
REGRESSION_FLOOR = 1.0 / 1.02
PROMOTION_FLOOR = 1.05
TWO_SIDED_T95_DF2 = 4.302652729911275
LABEL_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
HEX40_PATTERN = re.compile(r"^[0-9a-f]{40}$")
HEX64_PATTERN = re.compile(r"^[0-9a-f]{64}$")
MAX_BINARY_BYTES = 512 * 1024 * 1024
MAX_BUILD_OUTPUT = 8 * 1024 * 1024
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


@dataclass(frozen=True)
class Cell:
    identifier: str
    profile: str
    k: int
    r: int
    shard_bytes: int
    losses: int
    role: str
    seed: int

    def json(self) -> dict[str, Any]:
        return {
            "id": self.identifier,
            "profile": self.profile,
            "K": self.k,
            "R": self.r,
            "shard_bytes": self.shard_bytes,
            "loss_count": self.losses,
            "reuse": REUSE,
            "iterations": RETAINED_SAMPLES,
            "warmup": WARMUP,
            "seed": self.seed,
            "role": self.role,
        }


# Nine losses are above every production direct-repair limit and therefore
# make the reported path a transform traversal.  The maximum-loss boundary
# cell uses all 512 recovery shards.
CELLS = (
    Cell("high-live-96m", "high", 300, 100, 393216, 9,
         "generic_live_set_boundary", 0x463000),
    Cell("high-live-256m", "high", 300, 100, 1048576, 9,
         "generic_live_set_above_cache", 0x463001),
    Cell("low-live-128m", "low", 100, 193, 524288, 9,
         "low_profile_live_set_boundary", 0x463002),
    Cell("low-live-256m", "low", 100, 193, 1048576, 9,
         "low_profile_live_set_above_cache", 0x463003),
    Cell("high-side-256-64k", "high", 1000, 200, 65536, 9,
         "legacy_high_side_256", 0x463004),
    Cell("high-side-512-64k", "high", 2000, 500, 65536, 9,
         "legacy_high_side_512_below_large_l3", 0x463005),
    Cell("high-side-512-128k", "high", 2000, 500, 131072, 9,
         "legacy_high_side_512_above_large_l3", 0x463006),
    Cell("low-side-256-loss9", "low", 200, 800, 65536, 9,
         "low_profile_side_override_diagnostic", 0x463007),
    Cell("high-side-512-maxloss", "high", 2000, 512, 65536, 512,
         "actual_decode_rows_maximum_loss", 0x463008),
)


@dataclass(frozen=True)
class LaneRequest:
    label: str
    root: Path
    commit: str
    binary: Path | None
    build: bool


class EvidenceError(RuntimeError):
    """Evidence cannot be accepted under the declared campaign contract."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def load_support() -> Any:
    path = Path(__file__).resolve().parents[1] / "main_compare" / "run_abba.py"
    specification = importlib.util.spec_from_file_location(
        "leopard2_l3_tiling_main_compare_support", path)
    require(specification is not None and specification.loader is not None,
            f"cannot load benchmark support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()


def canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise EvidenceError(f"value is not canonical JSON: {error}") from error


def sha256(path: Path) -> str:
    resolved = path.resolve(strict=True)
    before = os.lstat(resolved)
    require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
            0 <= before.st_size <= MAX_BINARY_BYTES,
            f"file is not a bounded single-link regular file: {resolved}")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    descriptor = os.open(resolved, flags)
    try:
        initial = os.fstat(descriptor)
        require((initial.st_dev, initial.st_ino) ==
                (before.st_dev, before.st_ino) and
                stat.S_ISREG(initial.st_mode) and initial.st_nlink == 1 and
                initial.st_size <= MAX_BINARY_BYTES,
                f"file identity changed before hashing: {resolved}")
        digest = hashlib.sha256()
        retained = 0
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            retained += len(block)
            require(retained <= MAX_BINARY_BYTES,
                    f"file grew beyond the identity bound: {resolved}")
            digest.update(block)
        final = os.fstat(descriptor)
        after = os.lstat(resolved)
        require(retained == initial.st_size and
                (final.st_dev, final.st_ino, final.st_size,
                 final.st_mtime_ns, final.st_ctime_ns, final.st_nlink) ==
                (initial.st_dev, initial.st_ino, initial.st_size,
                 initial.st_mtime_ns, initial.st_ctime_ns, initial.st_nlink) and
                (after.st_dev, after.st_ino) ==
                (initial.st_dev, initial.st_ino),
                f"file changed while hashing: {resolved}")
        return digest.hexdigest()
    finally:
        os.close(descriptor)


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    metadata = os.lstat(resolved)
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
            f"identity path is not a single-link regular file: {resolved}")
    digest = sha256(resolved)
    return {
        "path": str(resolved),
        "sha256": digest,
        "size": metadata.st_size,
        "mode": stat.S_IMODE(metadata.st_mode),
        "uid": metadata.st_uid,
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "mtime_ns": metadata.st_mtime_ns,
        "ctime_ns": metadata.st_ctime_ns,
    }


def verify_file_identity(identity: Mapping[str, Any], *, executable: bool) -> None:
    path = Path(str(identity.get("path", ""))).resolve(strict=True)
    metadata = os.lstat(path)
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
            metadata.st_uid == identity.get("uid") and
            metadata.st_size == identity.get("size") and
            stat.S_IMODE(metadata.st_mode) == identity.get("mode") and
            metadata.st_dev == identity.get("device") and
            metadata.st_ino == identity.get("inode") and
            metadata.st_mtime_ns == identity.get("mtime_ns") and
            metadata.st_ctime_ns == identity.get("ctime_ns") and
            sha256(path) == identity.get("sha256"),
            f"frozen artifact identity changed: {path}")
    require(not executable or os.access(path, os.X_OK),
            f"frozen benchmark is no longer executable: {path}")


def write_exclusive(path: Path, value: object) -> None:
    payload = json.dumps(
        value, indent=2, sort_keys=True, allow_nan=False
    ).encode("utf-8") + b"\n"
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | \
        getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags, 0o600)
    try:
        written = 0
        while written != len(payload):
            count = os.write(descriptor, payload[written:])
            require(count > 0, f"short write: {path}")
            written += count
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def git(root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["/usr/bin/git", "-C", str(root), *arguments],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env={
            "GIT_NO_REPLACE_OBJECTS": "1",
            "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
        },
        check=False)
    require(completed.returncode == 0,
            f"git {' '.join(arguments)} failed for {root}: "
            f"{completed.stderr.strip()}")
    return completed.stdout.strip()


def require_commit(value: str, label: str) -> str:
    require(HEX40_PATTERN.fullmatch(value) is not None,
            f"{label} must be one lowercase 40-hex commit")
    return value


def source_identity(request: LaneRequest) -> dict[str, Any]:
    root = request.root.resolve(strict=True)
    require(root.is_dir() and
            Path(git(root, "rev-parse", "--show-toplevel")).resolve(
                strict=True) == root,
            f"lane {request.label} root is not its Git worktree top: {root}")
    require(not git(root, "for-each-ref", "--format=%(refname)", "refs/replace"),
            f"lane {request.label} repository contains replace refs")
    head = require_commit(git(root, "rev-parse", "HEAD"),
                          f"lane {request.label} HEAD")
    tree = require_commit(git(root, "rev-parse", "HEAD^{tree}"),
                          f"lane {request.label} tree")
    require(head == request.commit,
            f"lane {request.label} HEAD {head} != requested {request.commit}")
    status = git(
        root, "status", "--porcelain=v1", "--untracked-files=normal",
        "--ignore-submodules=none")
    require(status == "", f"lane {request.label} source is dirty: {status!r}")
    return {
        "label": request.label,
        "root": str(root),
        "commit": head,
        "tree": tree,
        "status": "clean",
        "status_sha256": hashlib.sha256(b"").hexdigest(),
    }


def cmake_home(cache: Path) -> Path:
    prefix = "CMAKE_HOME_DIRECTORY:INTERNAL="
    matches = [
        line[len(prefix):] for line in cache.read_text(
            encoding="utf-8", errors="strict").splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1 and matches[0],
            f"{cache}: expected one CMAKE_HOME_DIRECTORY")
    return Path(matches[0]).resolve(strict=True)


def find_build_root(binary: Path, source_root: Path) -> Path:
    resolved = binary.resolve(strict=True)
    require(resolved.is_file() and os.access(resolved, os.X_OK),
            f"benchmark is absent or not executable: {resolved}")
    for candidate in (resolved.parent, *resolved.parents):
        cache = candidate / "CMakeCache.txt"
        if cache.is_file():
            require(cmake_home(cache) == source_root.resolve(strict=True),
                    f"{cache}: CMake source does not match lane worktree")
            return candidate
        if candidate == source_root:
            break
    raise EvidenceError(f"cannot find the lane CMakeCache above {resolved}")


def build_lane(request: LaneRequest, output: Path, jobs: int) -> tuple[Path, dict[str, Any]]:
    build_root = output / "builds" / request.label
    build_root.parent.mkdir(mode=0o700, exist_ok=True)
    require(not build_root.exists(), f"lane build directory exists: {build_root}")
    configure = [
        "/usr/bin/cmake", "-S", str(request.root), "-B", str(build_root),
        "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
        "-DLEO2_BUILD_TESTS=OFF", "-DLEO2_BUILD_BENCHMARKS=ON",
        "-DLEO2_ENABLE_CUDA=OFF", "-DLEO2_BACKEND_VARIANT=auto",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
    ]
    build = [
        "/usr/bin/cmake", "--build", str(build_root),
        "--target", "bench_leopard2", "-j", str(jobs),
    ]
    records = []
    for name, command in (("configure", configure), ("build", build)):
        completed = SUPPORT.run_process_bounded(
            command, environment=CHILD_ENVIRONMENT, timeout=1800.0,
            max_stdout=MAX_BUILD_OUTPUT, max_stderr=MAX_BUILD_OUTPUT)
        records.append({
            "name": name,
            "command": command,
            "returncode": completed.returncode,
            "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
            "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
            "stdout": completed.stdout.decode("utf-8", errors="replace"),
            "stderr": completed.stderr.decode("utf-8", errors="replace"),
        })
        require(completed.returncode == 0,
                f"{request.label} {name} failed ({completed.returncode})")
    binary = (build_root / "bench_leopard2").resolve(strict=True)
    return binary, {
        "mode": "built_by_runner",
        "build_root": str(build_root.resolve()),
        "jobs": jobs,
        "steps": records,
    }


def existing_lane_binary(request: LaneRequest) -> tuple[Path, dict[str, Any]]:
    require(request.binary is not None, "existing lane binary is absent")
    binary = request.binary
    if not binary.is_absolute():
        binary = request.root / binary
    resolved = binary.resolve(strict=True)
    return resolved, {"mode": "caller_supplied"}


def source_binary_identity(
    binary: Path, source: Mapping[str, Any], build_record: Mapping[str, Any]
) -> dict[str, Any]:
    root = Path(str(source["root"]))
    binary_identity = file_identity(binary)
    try:
        build_root = find_build_root(binary, root)
    except EvidenceError:
        # A caller may already have copied the executable into an immutable
        # artifact directory that intentionally carries no mutable build tree.
        # The mandatory runtime schema-v5 probe below is the authoritative
        # compiled-source binding for this mode.
        require(build_record.get("mode") == "caller_supplied",
                "runner-built benchmark lost its CMake build identity")
        return {
            "build": dict(build_record),
            "build_metadata": "external_frozen_binary_runtime_attestation_required",
            "binary": binary_identity,
        }
    cache = build_root / "CMakeCache.txt"
    generated_candidates = [
        build_root / "generated" / "leopard2-benchmark-attestation" /
            "leopard2_benchmark_source_attestation.h",
        build_root / "source_attestation.h",
    ]
    generated_existing = [
        path.resolve() for path in generated_candidates if path.is_file()
    ]
    require(len(generated_existing) == 1,
            "benchmark source-attestation header is absent or ambiguous")
    generated = generated_existing[0]
    header = generated.read_text(encoding="ascii", errors="strict")
    require(
        f'#define LEO2_BENCHMARK_SOURCE_COMMIT "{source["commit"]}"' in header and
        f'#define LEO2_BENCHMARK_SOURCE_TREE "{source["tree"]}"' in header and
        "#define LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY 0" in header,
        "generated benchmark source attestation differs from lane source")
    closure: dict[str, Any] = {}
    archive_candidates = sorted({
        path.resolve() for name in ("libleopard.a", "leopard.lib")
        for path in build_root.rglob(name) if path.is_file()
    })
    require(len(archive_candidates) == 1,
            "CMake build must contain exactly one production Leopard archive")
    closure["production_archive"] = file_identity(archive_candidates[0])
    optional_paths: dict[str, list[Path]] = {
        "compile_commands": [build_root / "compile_commands.json"],
        "benchmark_link_recipe": list(
            build_root.rglob("CMakeFiles/bench_leopard2.dir/link.txt")) +
            [build_root / "bench_link.txt"],
        "archive_link_recipe": [
            *build_root.rglob("CMakeFiles/leopard.dir/link.txt"),
            *build_root.rglob("CMakeFiles/libleopard.dir/link.txt"),
            build_root / "archive_link.txt",
        ],
    }
    for name, candidates in optional_paths.items():
        existing = sorted({
            path.resolve() for path in candidates if path.is_file()
        })
        require(len(existing) <= 1,
                f"CMake build contains ambiguous {name} files")
        if existing:
            closure[name] = file_identity(existing[0])
    return {
        "build": dict(build_record),
        "build_root": str(build_root),
        "binary": binary_identity,
        "cmake_cache": file_identity(cache),
        "generated_source_attestation": file_identity(generated),
        "build_closure": closure,
    }


def freeze_binary(
    source_binary: Path, destination: Path, expected_sha256: str
) -> dict[str, Any]:
    require(HEX64_PATTERN.fullmatch(expected_sha256) is not None,
            "source binary SHA-256 is invalid")
    destination.parent.mkdir(mode=0o700, parents=True, exist_ok=False)
    source_path = source_binary.resolve(strict=True)
    source_before = os.lstat(source_path)
    require(stat.S_ISREG(source_before.st_mode) and source_before.st_nlink == 1 and
            source_before.st_size <= MAX_BINARY_BYTES,
            f"source benchmark is not a bounded regular file: {source_path}")
    source_flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
        getattr(os, "O_NOFOLLOW", 0)
    destination_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | \
        getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    source_fd = os.open(source_path, source_flags)
    destination_fd = os.open(destination, destination_flags, 0o500)
    digest = hashlib.sha256()
    copied = 0
    try:
        initial = os.fstat(source_fd)
        require((initial.st_dev, initial.st_ino) ==
                (source_before.st_dev, source_before.st_ino),
                "source benchmark changed before freezing")
        while True:
            block = os.read(source_fd, 1024 * 1024)
            if not block:
                break
            copied += len(block)
            require(copied <= MAX_BINARY_BYTES,
                    "source benchmark exceeds the freeze bound")
            digest.update(block)
            offset = 0
            while offset != len(block):
                written = os.write(destination_fd, block[offset:])
                require(written > 0, "short write while freezing benchmark")
                offset += written
        os.fsync(destination_fd)
        os.fchmod(destination_fd, 0o500)
        source_after = os.fstat(source_fd)
        require(
            copied == initial.st_size and
            (source_after.st_dev, source_after.st_ino, source_after.st_size,
             source_after.st_mtime_ns, source_after.st_ctime_ns,
             source_after.st_nlink) ==
            (initial.st_dev, initial.st_ino, initial.st_size,
             initial.st_mtime_ns, initial.st_ctime_ns, initial.st_nlink),
            "source benchmark changed while freezing")
    finally:
        os.close(destination_fd)
        os.close(source_fd)
    actual_sha256 = digest.hexdigest()
    require(actual_sha256 == expected_sha256,
            "frozen bytes differ from pre-freeze source identity")
    identity = file_identity(destination)
    require(identity["sha256"] == expected_sha256 and
            identity["mode"] == 0o500 and identity["uid"] == os.getuid(),
            "frozen benchmark identity is not immutable for the campaign")
    os.chmod(destination.parent, 0o500)
    verify_file_identity(identity, executable=True)
    return {
        "schema": ARTIFACT_SCHEMA,
        "source_path": str(source_path),
        "identity": identity,
    }


class CanonicalLock:
    """Own the exact campaign-wide lock for build, test, and timing."""

    def __init__(self) -> None:
        self.descriptor: int | None = None
        self.identity: dict[str, Any] | None = None

    def __enter__(self) -> dict[str, Any]:
        before = None
        try:
            before = os.lstat(CANONICAL_LOCK)
        except FileNotFoundError:
            pass
        flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0) | \
            getattr(os, "O_NOFOLLOW", 0)
        self.descriptor = os.open(CANONICAL_LOCK, flags, 0o600)
        metadata = os.fstat(self.descriptor)
        path_metadata = os.lstat(CANONICAL_LOCK)
        require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o600 and
                (metadata.st_dev, metadata.st_ino) ==
                (path_metadata.st_dev, path_metadata.st_ino) and
                (before is None or
                 (metadata.st_dev, metadata.st_ino) ==
                 (before.st_dev, before.st_ino)),
                "canonical benchmark lock has unsafe identity or permissions")
        waiting_ns = time.monotonic_ns()
        print(f"waiting for canonical lock {CANONICAL_LOCK}", flush=True)
        fcntl.flock(self.descriptor, fcntl.LOCK_EX)
        acquired_ns = time.monotonic_ns()
        self.identity = {
            "path": str(CANONICAL_LOCK),
            "device": metadata.st_dev,
            "inode": metadata.st_ino,
            "uid": metadata.st_uid,
            "mode": stat.S_IMODE(metadata.st_mode),
            "lock": "blocking_exclusive_flock",
            "waiting_monotonic_ns": waiting_ns,
            "acquired_monotonic_ns": acquired_ns,
        }
        print("acquired canonical lock", flush=True)
        return self.identity

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        descriptor = self.descriptor
        self.descriptor = None
        if descriptor is not None:
            try:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            finally:
                os.close(descriptor)


def validate_topology(cpu: int, sibling: int) -> dict[str, Any]:
    require((cpu, sibling) in ALLOWED_CPU_PAIRS,
            f"CPU pair must be one of {sorted(ALLOWED_CPU_PAIRS)}")
    initial_affinity = set(os.sched_getaffinity(0))
    require(cpu in initial_affinity,
            f"benchmark CPU {cpu} is outside launch affinity")
    records = []
    for logical in (cpu, sibling):
        topology = Path(f"/sys/devices/system/cpu/cpu{logical}/topology")
        siblings = SUPPORT.parse_cpu_list(
            (topology / "thread_siblings_list").read_text(
                encoding="ascii").strip())
        records.append({
            "cpu": logical,
            "thread_siblings": sorted(siblings),
            "core_id": int((topology / "core_id").read_text(
                encoding="ascii").strip()),
            "package_id": int((topology / "physical_package_id").read_text(
                encoding="ascii").strip()),
        })
    require(all(record["thread_siblings"] == sorted((cpu, sibling))
                for record in records),
            "requested CPUs are not the exact mutual SMT pair")
    require(records[0]["core_id"] == records[1]["core_id"] and
            records[0]["package_id"] == records[1]["package_id"],
            "requested CPUs do not identify one physical core")
    return {
        "launch_affinity": sorted(initial_affinity),
        "records": records,
    }


def pin_runner(cpu: int) -> None:
    os.sched_setaffinity(0, {cpu})
    require(set(os.sched_getaffinity(0)) == {cpu},
            f"runner could not pin itself to CPU {cpu}")


def benchmark_command(
    executable: str, cell: Cell, cpu: int, *, attest_source: bool
) -> list[str]:
    if attest_source:
        arguments = [
            "--k", "8", "--r", "4", "--profile", "high",
            "--field", "gf16", "--backend", "avx2", "--bytes", "64",
            "--loss", "4", "--batch", "1", "--reuse", "1",
            "--iterations", "1", "--warmup", "1", "--threads", "1",
            "--seed", "4603904", "--skip-legacy", "--retain-samples",
            "--attest-source", "--json", "-",
        ]
    else:
        arguments = [
            "--k", str(cell.k), "--r", str(cell.r),
            "--profile", cell.profile, "--field", "gf16",
            "--backend", "avx2", "--bytes", str(cell.shard_bytes),
            "--loss", str(cell.losses), "--batch", "1",
            "--reuse", str(REUSE), "--iterations", str(RETAINED_SAMPLES),
            "--warmup", str(WARMUP), "--threads", "1",
            "--seed", str(cell.seed), "--skip-legacy", "--retain-samples",
            "--report-decode-path", "--json", "-",
        ]
    return ["/usr/bin/taskset", "-c", str(cpu), executable, *arguments]


def result_metric(
    result: Mapping[str, Any], name: str, sample_count: int
) -> tuple[float, list[float]]:
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "benchmark metrics are absent")
    metric = metrics.get(name)
    require(isinstance(metric, dict), f"benchmark metric is absent: {name}")
    value = metric.get("median_us_per_batch_call")
    samples = metric.get("samples_us_per_batch_call")
    require(isinstance(value, (int, float)) and not isinstance(value, bool) and
            math.isfinite(float(value)) and float(value) > 0,
            f"benchmark median is invalid: {name}")
    require(isinstance(samples, list) and len(samples) == sample_count and
            all(isinstance(item, (int, float)) and not isinstance(item, bool) and
                math.isfinite(float(item)) and float(item) > 0
                for item in samples),
            f"retained benchmark samples are invalid: {name}")
    return float(value), [float(item) for item in samples]


def validate_attestation_result(
    result: object, source: Mapping[str, Any]
) -> dict[str, Any]:
    require(isinstance(result, dict) and
            result.get("schema") == "leopard2-benchmark-v5",
            "source attestation did not emit benchmark schema v5")
    build = result.get("build")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    require(isinstance(build, dict) and isinstance(resolved, dict) and
            isinstance(correctness, dict),
            "source attestation result is incomplete")
    require(build.get("source_commit") == source["commit"] and
            build.get("source_tree") == source["tree"] and
            build.get("source_tracked_dirty") is False,
            "compiled source attestation differs from requested clean source")
    require(resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf16" and
            resolved.get("backend") == "avx2" and
            resolved.get("thread_count") == 1 and
            correctness.get("leopard2_round_trip") is True,
            "source-attestation probe resolved unexpected codec semantics")
    result_metric(result, "encode_execution", 1)
    result_metric(result, "decode_execution", 1)
    return {
        "source_commit": build["source_commit"],
        "source_tree": build["source_tree"],
        "source_tracked_dirty": build["source_tracked_dirty"],
        "resolved": dict(resolved),
    }


def validate_timing_result(result: object, cell: Cell) -> dict[str, Any]:
    require(isinstance(result, dict) and
            result.get("schema") == "leopard2-benchmark-v3",
            "timed benchmark did not emit schema v3")
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(isinstance(parameters, dict) and isinstance(resolved, dict) and
            isinstance(correctness, dict) and isinstance(digests, dict),
            "timed benchmark identity/correctness fields are incomplete")
    expected_profile = "legacy_high_v1" if cell.profile == "high" else "low_v1"
    require(
        parameters.get("K") == cell.k and parameters.get("R") == cell.r and
        parameters.get("shard_bytes") == cell.shard_bytes and
        parameters.get("loss_count") == cell.losses and
        parameters.get("reuse") == REUSE and
        parameters.get("iterations") == RETAINED_SAMPLES and
        parameters.get("warmup") == WARMUP and
        parameters.get("seed") == cell.seed and
        parameters.get("batch") == 1 and
        parameters.get("thread_count") == 1 and
        parameters.get("requested_profile") == expected_profile and
        parameters.get("requested_field") == "gf16" and
        parameters.get("requested_backend") == "avx2" and
        parameters.get("skip_legacy") is True and
        parameters.get("retain_samples") is True and
        parameters.get("report_decode_path") is True,
        "benchmark parameters differ from the predeclared cell")
    require(
        resolved.get("profile") == expected_profile and
        resolved.get("field") == "gf16" and
        resolved.get("backend") == "avx2" and
        resolved.get("thread_count") == 1 and
        resolved.get("selected_decode_path") == "tiled" and
        resolved.get("selected_decode_rule") == "workspace_tiled" and
        isinstance(resolved.get("decode_required_work_slots"), int) and
        resolved["decode_required_work_slots"] > 0,
        "benchmark resolved an unexpected profile/backend/decode path")
    require(correctness.get("leopard2_round_trip") is True,
            "benchmark round trip failed")
    require(digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and len(digests[name]) == 16
                for name in (
                    "original_data", "transmitted_parity",
                    "recovered_originals")),
            "benchmark workload digests are incomplete")
    encode, encode_samples = result_metric(
        result, "encode_execution", RETAINED_SAMPLES)
    decode, decode_samples = result_metric(
        result, "decode_execution", RETAINED_SAMPLES)
    return {
        "encode_us": encode,
        "encode_samples_us": encode_samples,
        "decode_us": decode,
        "decode_samples_us": decode_samples,
        "digests": dict(digests),
        "resolved": dict(resolved),
    }


def run_one(
    lane: Mapping[str, Any], cell: Cell, cpu: int, timeout: float,
    *, attest_source: bool,
) -> dict[str, Any]:
    artifact = lane["artifact"]["identity"]
    verify_file_identity(artifact, executable=True)
    command = benchmark_command(
        str(artifact["path"]), cell, cpu, attest_source=attest_source)
    started_ns = time.monotonic_ns()
    completed = SUPPORT.run_process_bounded(
        command, environment=CHILD_ENVIRONMENT, timeout=timeout,
        max_stdout=16 * 1024 * 1024, max_stderr=2 * 1024 * 1024)
    duration_ns = time.monotonic_ns() - started_ns
    verify_file_identity(artifact, executable=True)
    require(completed.returncode == 0,
            f"{lane['label']} benchmark failed ({completed.returncode}): " +
            completed.stderr.decode("utf-8", errors="replace")[-2000:])
    try:
        parsed = json.loads(completed.stdout.decode("utf-8", errors="strict"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(
            f"{lane['label']} stdout is not one JSON object: {error}") from error
    normalized = (
        validate_attestation_result(parsed, lane["source"])
        if attest_source else validate_timing_result(parsed, cell)
    )
    return {
        "lane": lane["label"],
        "command": command,
        "environment": dict(CHILD_ENVIRONMENT),
        "runner_affinity": sorted(os.sched_getaffinity(0)),
        "duration_ns": duration_ns,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "stderr": completed.stderr.decode("utf-8", errors="replace"),
        "normalized": normalized,
        "result": parsed,
    }


def isolation_window(
    cpu: int, sibling: int, pair_lease: Mapping[str, Any],
    operation: Any,
) -> tuple[Any, dict[str, Any]]:
    before_cpu = SUPPORT.cpu_stat_snapshot(cpu)
    before_sibling = SUPPORT.cpu_stat_snapshot(sibling)
    before_ns = time.monotonic_ns()
    result = operation()
    after_ns = time.monotonic_ns()
    isolation = SUPPORT.isolation_record(
        cpu, sibling, pair_lease, before_ns, after_ns,
        before_cpu, SUPPORT.cpu_stat_snapshot(cpu),
        before_sibling, SUPPORT.cpu_stat_snapshot(sibling))
    require(isolation["accepted"] is True,
            "reserved SMT sibling was active during measured work")
    return result, isolation


def ratio_summary(values: list[float]) -> dict[str, Any]:
    require(len(values) == ROUNDS and all(value > 0 for value in values),
            "three positive ABBA round ratios are required")
    logs = [math.log(value) for value in values]
    center = statistics.mean(logs)
    half = TWO_SIDED_T95_DF2 * statistics.stdev(logs) / math.sqrt(ROUNDS)
    speedup = math.exp(center)
    lower = math.exp(center - half)
    upper = math.exp(center + half)
    return {
        "round_left_time_over_right_time": values,
        "speedup_of_right_over_left": speedup,
        "confidence_interval_95": [lower, upper],
        "right_slowdown_fraction": 1.0 / speedup - 1.0,
        "regression_floor_left_time_over_right_time": REGRESSION_FLOOR,
        "promotion_floor_left_time_over_right_time": PROMOTION_FLOOR,
        "point_estimate_right_regression_over_2_percent":
            speedup < REGRESSION_FLOOR,
        "credible_right_regression_over_2_percent":
            upper < REGRESSION_FLOOR,
        "credible_right_improvement_at_least_5_percent":
            lower >= PROMOTION_FLOOR,
    }


def summarize_cell(
    cell: Cell, rounds: list[dict[str, Any]], left: str, right: str
) -> dict[str, Any]:
    require(len(rounds) == ROUNDS, "cell did not complete three rounds")
    reference_digests = rounds[0]["invocations"][0]["normalized"]["digests"]
    paths: dict[str, list[dict[str, Any]]] = {left: [], right: []}
    for round_value in rounds:
        require(round_value["isolation"]["accepted"] is True,
                "contaminated round cannot be summarized")
        for invocation in round_value["invocations"]:
            require(invocation["normalized"]["digests"] == reference_digests,
                    "lane workload digests differ")
            lane = invocation["lane"]
            resolved = invocation["normalized"]["resolved"]
            paths[lane].append({
                "path": resolved["selected_decode_path"],
                "rule": resolved["selected_decode_rule"],
                "required_work_slots":
                    resolved["decode_required_work_slots"],
            })
    output: dict[str, Any] = {
        "cell": cell.json(),
        "digests": reference_digests,
        "decode_paths": {
            lane: {
                "unique": sorted({
                    canonical(value).decode("utf-8") for value in records
                }),
                "invocation_count": len(records),
            }
            for lane, records in paths.items()
        },
    }
    for output_name, metric_name in (
        ("encode", "encode_us"), ("decode", "decode_us")
    ):
        ratios = []
        for round_value in rounds:
            left_values = [
                item["normalized"][metric_name]
                for item in round_value["invocations"] if item["lane"] == left
            ]
            right_values = [
                item["normalized"][metric_name]
                for item in round_value["invocations"] if item["lane"] == right
            ]
            require(len(left_values) == 2 and len(right_values) == 2,
                    "ABBA round is incomplete")
            ratios.append(math.exp(
                statistics.mean(map(math.log, left_values)) -
                statistics.mean(map(math.log, right_values))))
        output[output_name] = ratio_summary(ratios)
    return output


def parse_lane(values: Sequence[str], *, build: bool) -> LaneRequest:
    expected = 3 if build else 4
    require(len(values) == expected,
            f"lane specification requires {expected} values")
    label, root_text, commit = values[:3]
    require(LABEL_PATTERN.fullmatch(label) is not None,
            f"invalid lane label: {label!r}")
    root = Path(root_text).resolve(strict=True)
    binary = None if build else Path(values[3])
    return LaneRequest(
        label=label, root=root,
        commit=require_commit(commit, f"lane {label} commit"),
        binary=binary, build=build)


def parse_args(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Authoritative L3-aware tiling ABBA comparison")
    group = parser.add_argument_group("lanes")
    group.add_argument(
        "--lane", nargs=4, action="append", default=[], metavar=(
            "LABEL", "SOURCE_ROOT", "COMMIT", "BENCH_BINARY"),
        help="use an existing CMake-built bench_leopard2")
    group.add_argument(
        "--build-lane", nargs=3, action="append", default=[], metavar=(
            "LABEL", "SOURCE_ROOT", "COMMIT"),
        help="configure/build/freeze bench_leopard2 inside the campaign")
    parser.add_argument(
        "--comparison", nargs=2, required=True, metavar=("LEFT", "RIGHT"),
        help="lane labels; results report RIGHT speedup over LEFT")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int, choices=(4, 14), required=True)
    parser.add_argument("--sibling", type=int, choices=(20, 30), required=True)
    parser.add_argument(
        "--cells", nargs="+", metavar="ID",
        help="run only these cell IDs (default: the complete canonical matrix)")
    parser.add_argument("--timeout", type=float, default=1800.0)
    options = parser.parse_args(arguments)
    requests = [
        *(parse_lane(value, build=False) for value in options.lane),
        *(parse_lane(value, build=True) for value in options.build_lane),
    ]
    require(len(requests) == 2, "exactly two --lane/--build-lane entries are required")
    labels = [request.label for request in requests]
    require(len(set(labels)) == 2, "lane labels must be unique")
    require(list(options.comparison) == labels or
            set(options.comparison) == set(labels),
            "--comparison must name the two declared lanes exactly")
    require((options.cpu, options.sibling) in ALLOWED_CPU_PAIRS,
            "--cpu/--sibling must select 4/20 or 14/30")
    require(math.isfinite(options.timeout) and
            1.0 <= options.timeout <= 3600.0,
            "--timeout must be between 1 and 3600 seconds")
    known_cells = {cell.identifier: cell for cell in CELLS}
    requested_cells = (
        [cell.identifier for cell in CELLS]
        if options.cells is None else options.cells
    )
    require(requested_cells and len(set(requested_cells)) == len(requested_cells),
            "--cells must name at least one unique cell")
    require(all(identifier in known_cells for identifier in requested_cells),
            "--cells contains an unknown cell ID")
    requested_set = set(requested_cells)
    # Canonical matrix order makes repeated campaigns byte-stable even if the
    # command line listed the selected cells in another order.
    options.selected_cells = tuple(
        cell for cell in CELLS if cell.identifier in requested_set)
    options.requests = requests
    return options


def output_is_safe(output: Path, roots: Sequence[Path]) -> Path:
    absolute = Path(os.path.abspath(output))
    require(not absolute.exists(), f"output already exists: {absolute}")
    for root in roots:
        resolved_root = root.resolve(strict=True)
        try:
            absolute.relative_to(resolved_root)
        except ValueError:
            continue
        completed = subprocess.run(
            ["/usr/bin/git", "-C", str(resolved_root), "check-ignore",
             "-q", str(absolute)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(completed.returncode == 0,
                "output inside a lane source root must be ignored by Git")
    return absolute


def prepare_lane(
    request: LaneRequest, output: Path, jobs: int
) -> dict[str, Any]:
    before = source_identity(request)
    if request.build:
        binary, build_record = build_lane(request, output, jobs)
    else:
        binary, build_record = existing_lane_binary(request)
    source_binary = source_binary_identity(binary, before, build_record)
    after_build = source_identity(request)
    require(after_build == before,
            f"lane {request.label} source changed during build capture")
    destination = output / "artifacts" / request.label / "bench_leopard2"
    artifact = freeze_binary(
        binary, destination, source_binary["binary"]["sha256"])
    return {
        "schema": ATTESTATION_SCHEMA,
        "label": request.label,
        "source": before,
        "source_binary": source_binary,
        "artifact": artifact,
    }


def verify_lane(lane: Mapping[str, Any]) -> None:
    verify_file_identity(lane["artifact"]["identity"], executable=True)
    source_binary = lane["source_binary"]
    verify_file_identity(source_binary["binary"], executable=True)
    for name in ("cmake_cache", "generated_source_attestation"):
        if name in source_binary:
            verify_file_identity(source_binary[name], executable=False)
    for identity in source_binary.get("build_closure", {}).values():
        verify_file_identity(identity, executable=False)
    source = lane["source"]
    request = LaneRequest(
        label=str(lane["label"]), root=Path(str(source["root"])),
        commit=str(source["commit"]), binary=None, build=False)
    require(source_identity(request) == source,
            f"lane {lane['label']} source changed during campaign")


def main(arguments: Sequence[str] | None = None) -> int:
    try:
        options = parse_args(arguments)
        output = output_is_safe(
            options.output, [request.root for request in options.requests])
        topology = validate_topology(options.cpu, options.sibling)
    except Exception as error:
        print(f"evidence rejected before output creation: {error}", file=sys.stderr)
        return 1

    output.mkdir(mode=0o700, parents=True, exist_ok=False)
    partial: dict[str, Any] = {
        "schema": RAW_SCHEMA,
        "started_utc": SUPPORT.utc_now(),
        "cells": [],
    }
    try:
        with CanonicalLock() as canonical_lock:
            jobs = min(128, max(1, len(os.sched_getaffinity(0))))
            lanes = {
                request.label: prepare_lane(request, output, jobs)
                for request in options.requests
            }
            left, right = options.comparison
            order = tuple(
                tuple(left if side == "left" else right for side in round_sides)
                for round_sides in ORDER_SIDES)
            manifest = {
                "schema": MANIFEST_SCHEMA,
                "created_utc": SUPPORT.utc_now(),
                "canonical_lock": canonical_lock,
                "comparison": {
                    "left": left,
                    "right": right,
                    "ratio_semantics":
                        "left time divided by right time; above one means right is faster",
                },
                "cells": [cell.json() for cell in options.selected_cells],
                "order_by_round": [list(value) for value in order],
                "rounds": ROUNDS,
                "retained_samples_per_invocation": RETAINED_SAMPLES,
                "field": "gf16",
                "backend": "avx2",
                "threads": 1,
                "batch": 1,
                "timeout_seconds": options.timeout,
                "statistics": {
                    "method":
                        "one mean-log contrast per independent ABBA round",
                    "confidence": 0.95,
                    "critical_distribution": "Student-t",
                    "degrees_of_freedom": 2,
                    "promotion_floor": PROMOTION_FLOOR,
                    "regression_floor": REGRESSION_FLOOR,
                },
                "isolation_policy": {
                    "canonical_lock": str(CANONICAL_LOCK),
                    "cpu_pair": [options.cpu, options.sibling],
                    "allowed_pairs": [list(value)
                                      for value in sorted(ALLOWED_CPU_PAIRS)],
                    "reserved_sibling_max_nonidle_jiffies": 0,
                    "contamination": "abort and emit failure; never summarize",
                },
                "lanes": lanes,
                "topology": topology,
            }
            write_exclusive(output / "manifest.json", manifest)
            partial["manifest"] = manifest
            partial["host"] = {
                "platform": platform.platform(),
                "processor": platform.processor(),
                "python": platform.python_version(),
                "cpu": options.cpu,
                "reserved_sibling": options.sibling,
            }
            partial["canonical_lock"] = canonical_lock
            partial["lanes"] = lanes

            pin_runner(options.cpu)
            with SUPPORT.StableLeaseAnchor(), \
                    SUPPORT.PairLease(
                        options.cpu, options.sibling) as pair_lease:
                partial["pair_lease"] = pair_lease

                pre_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
                pre_sibling = SUPPORT.cpu_stat_snapshot(options.sibling)
                pre_start = time.monotonic_ns()
                time.sleep(5.0)
                pre_end = time.monotonic_ns()
                presample = SUPPORT.isolation_record(
                    options.cpu, options.sibling, pair_lease,
                    pre_start, pre_end, pre_cpu,
                    SUPPORT.cpu_stat_snapshot(options.cpu), pre_sibling,
                    SUPPORT.cpu_stat_snapshot(options.sibling))
                partial["presample"] = presample
                require(
                    presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] == 0 and
                    presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                    "CPU pair was not quiet during the five-second presample")

                partial["source_attestation_probes"] = []
                for label in (left, right):
                    probe, isolation = isolation_window(
                        options.cpu, options.sibling, pair_lease,
                        lambda label=label: run_one(
                            lanes[label], options.selected_cells[0], options.cpu,
                            options.timeout, attest_source=True))
                    partial["source_attestation_probes"].append({
                        "lane": label,
                        "invocation": probe,
                        "isolation": isolation,
                    })

                for cell in options.selected_cells:
                    cell_raw: dict[str, Any] = {
                        "cell": cell.json(), "rounds": [],
                    }
                    for round_index, round_order in enumerate(order):
                        invocations, isolation = isolation_window(
                            options.cpu, options.sibling, pair_lease,
                            lambda round_order=round_order: [
                                run_one(
                                    lanes[label], cell, options.cpu,
                                    options.timeout, attest_source=False)
                                for label in round_order
                            ])
                        round_value = {
                            "round": round_index,
                            "order": list(round_order),
                            "invocations": invocations,
                            "isolation": isolation,
                        }
                        cell_raw["rounds"].append(round_value)
                        print(
                            f"completed {cell.identifier} round "
                            f"{round_index + 1}/{ROUNDS}", flush=True)
                    partial["cells"].append(cell_raw)

            for lane in lanes.values():
                verify_lane(lane)
            partial["completed_utc"] = SUPPORT.utc_now()
            write_exclusive(output / "raw.json", partial)
            summaries = [
                summarize_cell(cell, raw_cell["rounds"], left, right)
                for cell, raw_cell in zip(
                    options.selected_cells, partial["cells"])
            ]
            regressions = []
            promotions = []
            for summary in summaries:
                for metric in ("encode", "decode"):
                    if summary[metric][
                            "point_estimate_right_regression_over_2_percent"]:
                        regressions.append({
                            "cell": summary["cell"]["id"],
                            "metric": metric,
                            "result": summary[metric],
                        })
                    if summary[metric][
                            "credible_right_improvement_at_least_5_percent"]:
                        promotions.append({
                            "cell": summary["cell"]["id"],
                            "metric": metric,
                            "result": summary[metric],
                        })
            summary_document = {
                "schema": SUMMARY_SCHEMA,
                "status": (
                    "right_lane_regression_detected"
                    if regressions else "accepted"),
                "comparison": {"left": left, "right": right},
                "cell_count": len(summaries),
                "timed_process_count": len(summaries) * ROUNDS * 4,
                "attestation_process_count": 2,
                "all_digests_matched": True,
                "all_rounds_zero_sibling_nonidle": True,
                "right_lane_regressions_over_2_percent": regressions,
                "right_lane_credible_improvements_at_least_5_percent": promotions,
                "cells": summaries,
                "raw_sha256": sha256(output / "raw.json"),
                "manifest_sha256": sha256(output / "manifest.json"),
                "artifact_sha256": {
                    label: lane["artifact"]["identity"]["sha256"]
                    for label, lane in lanes.items()
                },
            }
            write_exclusive(output / "summary.json", summary_document)
            print(json.dumps({
                "status": summary_document["status"],
                "left": left,
                "right": right,
                "cells": len(summaries),
                "timed_processes":
                    summary_document["timed_process_count"],
                "regressions": len(regressions),
                "promotions": len(promotions),
                "output": str(output),
            }, sort_keys=True))
            return 0 if not regressions else 2
    except Exception as error:
        partial["failed_utc"] = SUPPORT.utc_now()
        partial["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
        }
        try:
            write_exclusive(output / "failure.json", {
                "schema": FAILURE_SCHEMA,
                "partial": partial,
            })
        except Exception as write_error:
            print(
                f"could not retain failure record: {write_error}",
                file=sys.stderr)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
