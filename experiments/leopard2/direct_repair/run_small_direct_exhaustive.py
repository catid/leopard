#!/usr/bin/env python3
"""Run, resume, and independently verify the exhaustive small-direct oracle.

The execution binary is accepted only when it is the exact CMake target linked
to a validated production Leopard archive.  Every source, object, recipe,
configuration, helper, archive, and executable identity in that closure is
copied into the campaign so that ``--verify --no-current-input-check`` remains
meaningful after the original build tree has disappeared.
"""

from __future__ import annotations

import argparse
import ctypes
import errno
import fcntl
import functools
import hashlib
import importlib.util
import itertools
import json
import math
import os
import re
import select
import selectors
import shlex
import signal
import stat
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA = "leopard2-small-direct-exhaustive-campaign/v3"
SHARD_SCHEMA = "leopard2-small-direct-exhaustive/v1"
BUILD_SCHEMA_V1 = "leopard2-small-direct-exhaustive-build/v1"
BUILD_SCHEMA = "leopard2-small-direct-exhaustive-build/v2"
BUNDLE_SCHEMA = "leopard2-small-direct-exhaustive-bundle/v1"
BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2 = (
    "leopard2-benchmark-build-configuration-attestation/v2")
BUILD_CONFIGURATION_ATTESTATION_SCHEMA = (
    "leopard2-benchmark-build-configuration-attestation/v3")
CLEAN_REBUILD_SCHEMA_V1 = "leopard2-small-direct-clean-rebuild/v1"
CLEAN_REBUILD_SCHEMA = "leopard2-small-direct-clean-rebuild/v2"
EXPECTED_PATTERNS = 1_982_812
DEFAULT_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")
CANONICAL_GIT = Path("/usr/bin/git")
CANONICAL_CMAKE = Path("/usr/bin/cmake")
CANONICAL_TASKSET = Path("/usr/bin/taskset")
MAX_FILE_BYTES = 256 * 1024 * 1024
MAX_JSON_BYTES = 16 * 1024 * 1024
MAX_SHARD_JSON_BYTES = 1 * 1024 * 1024
MAX_GUARDIAN_STATUS_BYTES = 32 * 1024
CGROUP_TEARDOWN_SECONDS = 5.0
SUPERVISOR_EOF_GRACE_SECONDS = 0.25
SUPERVISOR_GRACE_SECONDS = 3 * CGROUP_TEARDOWN_SECONDS + 5.0
SUPERVISOR_PYTHON_FLAGS = ("-I", "-S", "-B")
CURRENT_INTERPRETER_EXECUTABLE = "/proc/self/exe"
FNV_OFFSET = 14_695_981_039_346_656_037
FNV_PRIME = 1_099_511_628_211
U64_MASK = (1 << 64) - 1
MODE_CACHE_VARIABLE = "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"
PROCESS_SUPPORT_RELATIVE = "experiments/leopard2/main_compare/run_abba.py"
FORBIDDEN_EXPERIMENT_MACROS_V2 = (
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE",
)
FORBIDDEN_EXPERIMENT_MACROS = (
    *FORBIDDEN_EXPERIMENT_MACROS_V2,
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
)
MODE_DEFINITIONS = {
    1: "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=1",
    2: "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=2",
}
BUILD_CONFIGURATION_VARIABLES_V2 = (
    "CMAKE_BUILD_TYPE",
    "CMAKE_GENERATOR",
    "CMAKE_CONFIGURATION_TYPES",
    "CMAKE_CXX_COMPILER",
    "CMAKE_CXX_FLAGS",
    "CMAKE_CXX_FLAGS_DEBUG",
    "CMAKE_CXX_FLAGS_RELEASE",
    "CMAKE_CXX_FLAGS_RELWITHDEBINFO",
    "CMAKE_CXX_FLAGS_MINSIZEREL",
    "ENABLE_OPENMP",
    "LEOPARD_ENABLE_GF8",
    "LEOPARD_ENABLE_GF16",
    "LEO2_BACKEND_VARIANT",
    "LEO2_BENCHMARK_GIT_EXECUTABLE",
    "LEO2_BUILD_BENCHMARKS",
    "LEO2_BUILD_TESTS",
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE",
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE",
)
BUILD_CONFIGURATION_VARIABLES = (
    *BUILD_CONFIGURATION_VARIABLES_V2[:-1],
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
    BUILD_CONFIGURATION_VARIABLES_V2[-1],
)
CONTROLLED_CONFIGURATION_KEYS_V1 = (
    "CMAKE_BUILD_TYPE",
    "CMAKE_AR",
    "CMAKE_COMMAND",
    "CMAKE_GENERATOR",
    "CMAKE_EXPORT_COMPILE_COMMANDS",
    "CMAKE_C_COMPILER",
    "CMAKE_CXX_COMPILER",
    "CMAKE_CXX_FLAGS",
    "CMAKE_CXX_FLAGS_RELEASE",
    "CMAKE_LINKER",
    "CMAKE_MAKE_PROGRAM",
    "CMAKE_RANLIB",
    "ENABLE_OPENMP",
    "LEOPARD_ENABLE_GF8",
    "LEOPARD_ENABLE_GF16",
    "LEO2_BACKEND_VARIANT",
    "LEO2_BUILD_ALLK_DIAGNOSTIC",
    "LEO2_BUILD_BENCHMARKS",
    "LEO2_BUILD_FUZZERS",
    "LEO2_BUILD_TESTS",
    "LEO2_ENABLE_CUDA",
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE",
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE",
    "LEO2_FLAG_MAVX2",
    "LEO2_FLAG_MNO_AVX512F",
    "LEO2_FLAG_MAVX512F",
    "LEO2_FLAG_MAVX512BW",
    "LEO2_FLAG_MAVX512VL",
)
CONTROLLED_CONFIGURATION_KEYS = (
    *CONTROLLED_CONFIGURATION_KEYS_V1[:-6],
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
    *CONTROLLED_CONFIGURATION_KEYS_V1[-6:],
)
TARGET_NAME = "leopard2_small_direct_exhaustive_test"
BENCHMARK_NAME = "bench_leopard2"
BOUND_SOURCE_FILES = (
    "CMakeLists.txt",
    "leopard.cpp",
    "leopard2.cpp",
    "Leopard2Backend.cpp",
    "Leopard2Backend.h",
    "Leopard2BackendScalar.cpp",
    "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp",
    "Leopard2BackendAVX2Xor.cpp",
    "Leopard2BackendGFNI.cpp",
    "Leopard2CpuFeatures.cpp",
    "Leopard2Direct.h",
    "Leopard2Dispatch.h",
    "Leopard2Plan.cpp",
    "Leopard2Plan.h",
    "LeopardCommon.cpp",
    "LeopardCommon.h",
    "LeopardFF8.cpp",
    "LeopardFF8.h",
    "LeopardFF16.cpp",
    "LeopardFF16.h",
    "leopard.h",
    "leopard2.h",
    "bench/leopard2/benchmark.cpp",
    "cmake/GenerateBenchmarkSourceAttestation.cmake",
    "cmake/Leopard2BenchmarkAttestation.cmake",
    "tests/leopard2/direct_oracle.cpp",
    "tests/leopard2/direct_oracle.h",
    "tests/leopard2/test_small_direct_exhaustive.cpp",
    "tools/leopard2_build_provenance.py",
    "tools/leopard2_direct_encode_crossover.py",
    PROCESS_SUPPORT_RELATIVE,
    "experiments/leopard2/direct_repair/run_small_direct_exhaustive.py",
)
CHILD_ENV = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}
EXECUTION_POLICY = {
    "schema": "leopard2-small-direct-execution-policy/v3",
    "executable_binding": "linux-sealed-memfd",
    "supervisor_interpreter": (
        "linux-proc-self-exe-current-inode-base-prefix-argv0"),
    "supervisor_python_flags": list(SUPERVISOR_PYTHON_FLAGS),
    "containment": (
        "guardian-fork-server-subreaper-plus-delegated-cgroup-v2-per-shard"),
    "descendant_cleanup": (
        "shard-and-campaign-subreaper-proc-children-pidfd-kill-and-reap"),
    "output_capture": "bounded-nonblocking-pipes-before-atomic-publication",
    "parent_liveness": "inherited-lock-and-eof-pipe",
    "supervisor_eof_grace_seconds": SUPERVISOR_EOF_GRACE_SECONDS,
    "guardian_cleanup_grace_seconds": CGROUP_TEARDOWN_SECONDS,
    "supervisor_parent_death_signal": "prctl-pdeathsig-sigcont-with-ppid-race-check",
    "campaign_guardian": (
        "sealed-runner-fork-server-handshake-before-cgroup-create"),
    "orphaned_cgroup_cleanup": (
        "guardian-plus-last-validated-shard-supervisor-remove-campaign-scope"),
    "delegated_cgroup_cleanup": (
        "validated-no-link-recursive-postorder-within-exact-shard-root"),
}
LINUX_F_ADD_SEALS = getattr(fcntl, "F_ADD_SEALS", 1033)
LINUX_F_GET_SEALS = getattr(fcntl, "F_GET_SEALS", 1034)
LINUX_F_SEAL_SEAL = getattr(fcntl, "F_SEAL_SEAL", 0x0001)
LINUX_F_SEAL_SHRINK = getattr(fcntl, "F_SEAL_SHRINK", 0x0002)
LINUX_F_SEAL_GROW = getattr(fcntl, "F_SEAL_GROW", 0x0004)
LINUX_F_SEAL_WRITE = getattr(fcntl, "F_SEAL_WRITE", 0x0008)
LINUX_PR_SET_PDEATHSIG = 1
LINUX_PR_SET_CHILD_SUBREAPER = 36
LINUX_PR_GET_CHILD_SUBREAPER = 37


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def load_process_containment_support() -> tuple[Any, Path]:
    """Load the repository's audited subreaper/procfs/pidfd implementation."""
    path = Path(__file__).resolve().parents[3] / PROCESS_SUPPORT_RELATIVE
    resolved = path.resolve(strict=True)
    specification = importlib.util.spec_from_file_location(
        "leopard2_small_direct_exhaustive_process_containment", resolved)
    require(specification is not None and specification.loader is not None,
            "cannot import process-containment support")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    try:
        specification.loader.exec_module(module)
    except Exception:
        sys.modules.pop(specification.name, None)
        raise
    require(callable(getattr(module, "LinuxDescendantContainment", None)),
            "process-containment support contract changed")
    return module, resolved


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False) +
            "\n").encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def record_digest(value: Mapping[str, Any], digest_key: str = "digest") -> str:
    copy = dict(value)
    copy.pop(digest_key, None)
    return sha256_bytes(canonical_bytes(copy))


def decode_json_bytes(payload: bytes, label: str) -> Any:
    def object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            require(key not in result,
                    "JSON contains duplicate key %r" % key)
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise EvidenceError(
            "JSON contains non-finite constant %s" % value)

    def finite_float(value: str) -> float:
        parsed = float(value)
        require(math.isfinite(parsed),
                "JSON contains non-finite number %s" % value)
        return parsed

    try:
        return json.loads(
            payload.decode("utf-8", errors="strict"),
            object_pairs_hook=object_pairs,
            parse_constant=reject_constant, parse_float=finite_float)
    except (UnicodeError, ValueError, OverflowError, RecursionError) as error:
        raise EvidenceError("invalid JSON in %s: %s" % (label, error)) \
            from error


def file_snapshot(
        path: Path, maximum_bytes: int = MAX_FILE_BYTES,
        retain_content: bool = True,
) -> tuple[dict[str, Any], bytes]:
    require(type(maximum_bytes) is int and maximum_bytes >= 0,
            "file byte bound is invalid")
    resolved = Path(os.path.abspath(path))
    linked_before = resolved.lstat()
    require(resolved.resolve(strict=True) == resolved and
            stat.S_ISREG(linked_before.st_mode),
            "%s is not a canonical regular file" % resolved)
    descriptor = os.open(
        resolved, os.O_RDONLY | os.O_CLOEXEC |
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))
    try:
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode),
                "%s is not a regular file" % resolved)
        require(0 <= before.st_size <= maximum_bytes,
                "%s exceeds its retained byte bound" % resolved)
        digest = hashlib.sha256()
        size = 0
        chunks = [] if retain_content else None
        for block in iter(lambda: os.read(descriptor, 1 << 20), b""):
            require(size + len(block) <= maximum_bytes,
                    "%s grew beyond its retained byte bound" % resolved)
            digest.update(block)
            if chunks is not None:
                chunks.append(block)
            size += len(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    def stable_fields(value: os.stat_result) -> tuple[int, ...]:
        return (
            value.st_dev, value.st_ino, value.st_mode, value.st_uid,
            value.st_gid, value.st_nlink, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns,
        )

    require(stable_fields(linked_before) == stable_fields(before) ==
            stable_fields(after) and size == before.st_size,
            "%s changed while read" % resolved)
    linked_after = resolved.lstat()
    require(stable_fields(linked_after) == stable_fields(after),
            "%s pathname changed while read" % resolved)
    identity = {
        "path": str(resolved), "size": size, "sha256": digest.hexdigest(),
    }
    return identity, b"".join(chunks or ())


def load_json(path: Path, maximum_bytes: int = MAX_JSON_BYTES) -> Any:
    unused_identity, payload = file_snapshot(path, maximum_bytes)
    return decode_json_bytes(payload, str(path))


def file_identity(
        path: Path, maximum_bytes: int = MAX_FILE_BYTES) -> dict[str, Any]:
    identity, unused_payload = file_snapshot(
        path, maximum_bytes, retain_content=False)
    return identity


def _process_executable_inode_identity() -> dict[str, Any]:
    """Hash the immutable inode currently executing this Python process."""
    descriptor = os.open(
        "/proc/self/exe", os.O_RDONLY | os.O_CLOEXEC |
        getattr(os, "O_NONBLOCK", 0))
    try:
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode) and
                0 <= before.st_size <= MAX_FILE_BYTES,
                "current interpreter executable is invalid")
        digest = hashlib.sha256()
        size = 0
        for block in iter(lambda: os.read(descriptor, 1 << 20), b""):
            require(size + len(block) <= MAX_FILE_BYTES,
                    "current interpreter exceeds its retained byte bound")
            digest.update(block)
            size += len(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    def stable(value: os.stat_result) -> tuple[int, ...]:
        return (
            value.st_dev, value.st_ino, value.st_mode, value.st_size)
    require(stable(before) == stable(after) and size == before.st_size,
            "current interpreter changed while read")
    return {
        "device": after.st_dev, "inode": after.st_ino,
        "mode": after.st_mode, "size": size,
        "sha256": digest.hexdigest(),
    }


def process_executable_identity(
        retained_path: str | None = None) -> dict[str, Any]:
    """Capture one truthful path plus the currently executing Python inode."""
    current = _process_executable_inode_identity()
    if retained_path is None:
        reported = os.readlink("/proc/self/exe")
        require(reported.startswith("/") and
                not reported.endswith(" (deleted)"),
                "current interpreter has no retainable canonical path")
        retained_path = str(Path(reported).resolve(strict=True))
    require(isinstance(retained_path, str) and
            Path(retained_path).is_absolute() and
            Path(retained_path) == Path(os.path.abspath(retained_path)),
            "retained interpreter path is invalid")
    executable = Path(retained_path)
    linked_before = executable.lstat()
    require(executable.resolve(strict=True) == executable and
            stat.S_ISREG(linked_before.st_mode),
            "retained interpreter path is not a canonical regular file")
    linked_after = executable.lstat()
    linked_identity = (
        linked_before.st_dev, linked_before.st_ino,
        linked_before.st_mode, linked_before.st_size)
    linked_after_identity = (
        linked_after.st_dev, linked_after.st_ino,
        linked_after.st_mode, linked_after.st_size)
    current_identity = (
        current["device"], current["inode"],
        current["mode"], current["size"])
    require(linked_identity == linked_after_identity == current_identity,
            "retained interpreter pathname differs from its executing inode")
    return {
        "schema": "leopard2-process-executable/v1",
        "path": str(executable),
        **current,
    }


def process_executable_identity_matches(identity: Any) -> bool:
    """Revalidate the running inode without trusting its now-mutable path."""
    require(isinstance(identity, dict) and set(identity) == {
                "schema", "path", "device", "inode", "mode", "size",
                "sha256"} and
            identity.get("schema") == "leopard2-process-executable/v1" and
            isinstance(identity.get("path"), str) and
            Path(identity["path"]).is_absolute() and
            identity["path"] == str(Path(identity["path"])) and
            Path(identity["path"]) == Path(os.path.abspath(identity["path"])) and
            all(type(identity.get(name)) is int
                for name in ("device", "inode", "mode", "size")) and
            identity["device"] >= 0 and identity["inode"] >= 0 and
            stat.S_ISREG(identity["mode"]) and
            0 <= identity["size"] <= MAX_FILE_BYTES and
            isinstance(identity.get("sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", identity["sha256"]) is not None,
            "retained interpreter identity is malformed")
    current = _process_executable_inode_identity()
    return all(identity[name] == current[name] for name in current)


@functools.lru_cache(maxsize=1)
def interpreter_discovery_argv0() -> str:
    """Return a canonical argv[0] used only for Python prefix discovery.

    Internal helpers are executed through ``/proc/self/exe`` so that the
    kernel selects the exact inode already validated by
    :func:`process_executable_identity`.  CPython also uses argv[0] to locate
    its standard library, however.  Python 3.13 cannot complete isolated
    startup when ``/proc/self/exe`` names a deleted inode, as may happen after
    an executable-path replacement.  A canonical base-installation argv[0]
    preserves prefix discovery without changing the executable selected by
    the kernel.
    """
    base_prefix = Path(os.path.abspath(sys.base_prefix))
    require(base_prefix.is_absolute() and
            base_prefix.resolve(strict=True) == base_prefix and
            stat.S_ISDIR(base_prefix.stat().st_mode),
            "current interpreter base prefix is invalid")
    name = "python%d.%d%s" % (
        sys.version_info.major, sys.version_info.minor, sys.abiflags)
    candidate = base_prefix / "bin" / name
    linked = candidate.lstat()
    require(candidate.resolve(strict=True) == candidate and
            stat.S_ISREG(linked.st_mode) and os.access(candidate, os.X_OK),
            "current interpreter base executable is invalid")
    candidate_identity = file_identity(candidate)
    current_identity = _process_executable_inode_identity()
    require(candidate_identity["size"] == current_identity["size"] and
            candidate_identity["sha256"] == current_identity["sha256"],
            "current interpreter differs from its base executable")
    return str(candidate)


def launch_current_interpreter(
        command: Sequence[str], **options: Any) -> subprocess.Popen[bytes]:
    """Launch an isolated helper with the current executable inode."""
    require(isinstance(command, Sequence) and command and
            command[0] == interpreter_discovery_argv0() and
            "executable" not in options,
            "current interpreter launch command is invalid")
    return subprocess.Popen(
        command, executable=CURRENT_INTERPRETER_EXECUTABLE, **options)


def exact_file_identity(
        path: Path, expected: Mapping[str, Any], label: str) -> None:
    require(file_identity(path) == expected,
            "%s identity changed" % label)


def ensure_directory(path: Path, root: Path) -> Path:
    require(path.is_absolute() and
            (path == root or root in path.parents),
            "campaign directory escapes its root: %s" % path)
    if path != root:
        canonical_directory(path.parent, root, "campaign parent directory")
        if path.exists() or path.is_symlink():
            require(not path.is_symlink(),
                    "campaign directory is a symlink: %s" % path)
        else:
            os.mkdir(path, 0o700)
    canonical_directory(path, root, "campaign directory")
    return path


def canonical_artifact(path: Path, root: Path, label: str) -> None:
    linked = path.lstat()
    require(path.is_absolute() and
            (path == root or root in path.parents) and
            path.resolve(strict=True) == path and
            stat.S_ISREG(linked.st_mode) and linked.st_uid == os.getuid() and
            linked.st_nlink == 1,
            "%s escapes the canonical campaign tree" % label)


def canonical_directory(path: Path, root: Path, label: str) -> None:
    linked = path.lstat()
    require(path.is_absolute() and
            (path == root or root in path.parents) and
            path.resolve(strict=True) == path and
            stat.S_ISDIR(linked.st_mode) and linked.st_uid == os.getuid() and
            stat.S_IMODE(linked.st_mode) == 0o700 and
            linked.st_nlink >= 2,
            "%s escapes the canonical campaign tree" % label)


def cleanup_stale_temporaries(
        directory: Path, root: Path, patterns: Sequence[str]) -> None:
    """Remove only owned canonical staging files with an exact known name."""
    for path in tuple(directory.iterdir()):
        if not any(re.fullmatch(pattern, path.name) for pattern in patterns):
            continue
        canonical_artifact(path, root, "stale atomic staging file")
        require(path.lstat().st_size <= MAX_FILE_BYTES,
                "stale atomic staging file exceeds its byte bound")
        path.unlink()


def prepare_output_root(
        requested: Path, source_root: Path, workers: int) -> Path:
    output = Path(os.path.abspath(requested))
    require(output != Path("/") and
            not output.is_relative_to(source_root),
            "campaign output must be outside the source tree")
    if output.exists() or output.is_symlink():
        require(not output.is_symlink() and output.resolve(strict=True) ==
                    output,
                "campaign output must be a canonical directory")
    else:
        parent = output.parent
        linked_parent = parent.lstat()
        require(parent.resolve(strict=True) == parent and
                stat.S_ISDIR(linked_parent.st_mode) and
                os.access(parent, os.W_OK | os.X_OK),
                "campaign output parent is not a canonical "
                "writable directory")
        os.mkdir(output, 0o700)
    canonical_directory(output, output, "campaign output")
    cleanup_stale_temporaries(
        output, output,
        (r"(?:request|manifest)\.json\.tmp-[0-9]+",))

    allowed_top = {
        "request.json", "manifest.json", "artifacts",
        *("shard-%04d" % index for index in range(workers)),
    }
    require({path.name for path in output.iterdir()} <= allowed_top,
            "campaign contains unknown top-level artifacts")
    for path in output.iterdir():
        if path.name in ("request.json", "manifest.json"):
            canonical_artifact(path, output, "retained campaign metadata")
        elif path.name == "artifacts":
            canonical_directory(path, output, "retained artifact directory")
            cleanup_stale_temporaries(
                path, output,
                (re.escape(TARGET_NAME) + r"\.tmp-[0-9]+",))
            allowed_artifacts = {TARGET_NAME, "provenance"}
            require({item.name for item in path.iterdir()} <=
                    allowed_artifacts,
                    "artifact directory contains unknown entries")
            for item in path.iterdir():
                if item.name == "provenance":
                    canonical_directory(
                        item, output, "retained provenance directory")
                    cleanup_stale_temporaries(
                        item, output,
                        (r"file-[0-9]{4}\.tmp-[0-9]+",))
                    require(all(re.fullmatch(r"file-[0-9]{4}", child.name)
                                for child in item.iterdir()),
                            "provenance directory contains unknown entries")
                    for child in item.iterdir():
                        canonical_artifact(
                            child, output, "retained provenance entry")
                else:
                    canonical_artifact(
                        item, output, "retained frozen executable")
        else:
            canonical_directory(path, output, "retained shard directory")
            cleanup_stale_temporaries(
                path, output,
                (r"(?:stdout\.json|stderr\.txt|result\.json)"
                 r"\.tmp-[0-9]+",))
            require({item.name for item in path.iterdir()} <=
                    {"stdout.json", "stderr.txt", "result.json"},
                    "shard directory contains unknown entries")
            for item in path.iterdir():
                canonical_artifact(item, output, "retained shard entry")
    return output


def atomic_bytes(
        path: Path, payload: bytes, maximum_bytes: int,
        mode: int = 0o600) -> None:
    require(path.parent.is_dir() and
            path.parent.resolve(strict=True) == path.parent,
            "atomic publication parent is not a canonical directory")
    require(type(payload) is bytes and type(maximum_bytes) is int and
            0 <= len(payload) <= maximum_bytes <= MAX_FILE_BYTES and
            type(mode) is int and 0 <= mode <= 0o777,
            "atomic publication payload exceeds its byte bound")
    temporary = path.with_name(path.name + ".tmp-%d" % os.getpid())
    descriptor = -1
    try:
        descriptor = os.open(
            temporary,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0),
            mode)
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            require(written > 0, "short atomic publication write")
            offset += written
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = -1
        os.replace(temporary, path)
        directory = os.open(
            path.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary.exists() or temporary.is_symlink():
            canonical_artifact(
                temporary, path.parent,
                "failed atomic publication staging file")
            temporary.unlink()


def atomic_json(path: Path, value: Any) -> None:
    payload = json.dumps(
        value, sort_keys=True, indent=2,
        allow_nan=False).encode("utf-8") + b"\n"
    atomic_bytes(path, payload, MAX_JSON_BYTES)


def acquire_lock(path: Path, timeout: float):
    require(path == DEFAULT_LOCK and path.is_absolute() and
            math.isfinite(timeout) and timeout >= 0,
            "exhaustive campaign requires the canonical GF8 lock")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(
        path, os.O_RDWR | os.O_CREAT | os.O_CLOEXEC |
        getattr(os, "O_NOFOLLOW", 0),
        0o600)
    stream = os.fdopen(descriptor, "r+")
    deadline = time.monotonic() + timeout
    while True:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            os.fchmod(stream.fileno(), 0o600)
            lock_identity(stream, path)
            return stream
        except BlockingIOError:
            if time.monotonic() >= deadline:
                stream.close()
                raise EvidenceError("timed out waiting for canonical GF8 lock")
            time.sleep(1.0)
        except Exception:
            stream.close()
            raise


def close_campaign_lock(stream: Any) -> None:
    """Close locally without unlocking an inherited guardian descriptor.

    The lifecycle guardian inherits this descriptor with ``pass_fds``.  Those
    descriptors share one open-file description, so an explicit ``LOCK_UN``
    here would also release the guardian's lock while it may still be tearing
    down descendants.  Closing is sufficient: the kernel releases the flock
    after the final inherited descriptor is closed.
    """
    stream.close()


def lock_identity(stream: Any, path: Path) -> dict[str, Any]:
    require(path == DEFAULT_LOCK and path.is_absolute(),
            "noncanonical GF8 lock path")
    opened = os.fstat(stream.fileno())
    linked = path.lstat()
    require(stat.S_ISREG(opened.st_mode) and opened.st_uid == os.getuid() and
            opened.st_nlink == 1 and
            stat.S_IMODE(opened.st_mode) == 0o600 and
            opened.st_dev == linked.st_dev and
            opened.st_ino == linked.st_ino,
            "canonical GF8 lock was replaced or has unsafe metadata")
    return {
        "schema": "leopard2-global-benchmark-lock/v1",
        "path": str(path), "device": opened.st_dev,
        "inode": opened.st_ino, "uid": opened.st_uid,
        "mode": stat.S_IMODE(opened.st_mode), "nlink": opened.st_nlink,
        "mechanism": "fcntl-flock-exclusive",
    }


def fnv_mix_u64(digest: int, value: int) -> int:
    # Keep this independent from the C++ verifier.  Unrolling avoids making
    # exact descriptor reconstruction the dominant campaign cost.
    digest = ((digest ^ (value & 0xff)) * FNV_PRIME) & U64_MASK
    digest = ((digest ^ ((value >> 8) & 0xff)) * FNV_PRIME) & U64_MASK
    digest = ((digest ^ ((value >> 16) & 0xff)) * FNV_PRIME) & U64_MASK
    digest = ((digest ^ ((value >> 24) & 0xff)) * FNV_PRIME) & U64_MASK
    digest = ((digest ^ ((value >> 32) & 0xff)) * FNV_PRIME) & U64_MASK
    digest = ((digest ^ ((value >> 40) & 0xff)) * FNV_PRIME) & U64_MASK
    digest = ((digest ^ ((value >> 48) & 0xff)) * FNV_PRIME) & U64_MASK
    return ((digest ^ ((value >> 56) & 0xff)) * FNV_PRIME) & U64_MASK


@functools.lru_cache(maxsize=8)
def reconstruct_expected_shards(shard_count: int) -> tuple[dict[str, Any], ...]:
    require(type(shard_count) is int and 1 <= shard_count <= 128,
            "exact reconstruction requires 1..128 shards")
    states = [{
        "assigned_patterns": 0,
        "recovered_shards": 0,
        "verified_basis_symbols": 0,
        "digest_fnv1a64": FNV_OFFSET,
        "ordinal_digest_fnv1a64": FNV_OFFSET,
    } for _ in range(shard_count)]
    ordinal = 0
    for k in range(5, 17):
        for r in range(5, 9):
            for losses in range(5, min(k, r) + 1):
                for missing in itertools.combinations(range(k), losses):
                    for parities in itertools.combinations(range(r), losses):
                        state = states[ordinal % shard_count]
                        state["assigned_patterns"] += 1
                        state["recovered_shards"] += losses
                        state["verified_basis_symbols"] += losses * k
                        digest = state["digest_fnv1a64"]
                        digest = fnv_mix_u64(digest, k)
                        digest = fnv_mix_u64(digest, r)
                        digest = fnv_mix_u64(digest, losses)
                        for index in missing:
                            digest = fnv_mix_u64(digest, index)
                        for index in parities:
                            digest = fnv_mix_u64(digest, index)
                        state["digest_fnv1a64"] = digest
                        state["ordinal_digest_fnv1a64"] = fnv_mix_u64(
                            state["ordinal_digest_fnv1a64"], ordinal)
                        ordinal += 1
    require(ordinal == EXPECTED_PATTERNS,
            "independent descriptor enumeration changed pattern count")
    result = []
    for index, state in enumerate(states):
        result.append({
            "shard_index": index,
            "shard_count": shard_count,
            "total_patterns": EXPECTED_PATTERNS,
            "assigned_patterns": state["assigned_patterns"],
            "recovered_shards": state["recovered_shards"],
            "verified_basis_symbols": state["verified_basis_symbols"],
            "basis_seed": 0,
            "assignment": "global_ordinal_mod_shard_count",
            "digest_fnv1a64": "%016x" % state["digest_fnv1a64"],
            "ordinal_digest_fnv1a64":
                "%016x" % state["ordinal_digest_fnv1a64"],
        })
    require(sum(item["assigned_patterns"] for item in result) ==
            EXPECTED_PATTERNS,
            "independent descriptor shards do not cover the exact matrix")
    return tuple(result)


def selected_worker_cpus(workers: int) -> list[int]:
    allowed = sorted(os.sched_getaffinity(0))
    require(type(workers) is int and 1 <= workers <= 128 and
            workers <= len(allowed),
            "workers must be positive, at most 128, and no larger than "
            "the process-visible CPU set")
    return allowed[:workers]


def validate_run_limits(
        workers: int, timeout: float, lock_timeout: float) -> list[int]:
    require(type(workers) is int and
            type(timeout) in (int, float) and math.isfinite(timeout) and
            timeout > 0 and
            type(lock_timeout) in (int, float) and
            math.isfinite(lock_timeout) and lock_timeout >= 0,
            "workers/timeouts are invalid")
    return selected_worker_cpus(workers)


def require_mode_macro(
        arguments: Sequence[str], expected_mode: int | None,
        label: str) -> None:
    require(isinstance(arguments, list) and
            all(isinstance(argument, str) for argument in arguments),
            "%s has invalid compiler arguments" % label)
    references = [
        argument for argument in arguments
        if MODE_CACHE_VARIABLE in argument
    ]
    expected = (
        [] if expected_mode is None
        else [MODE_DEFINITIONS[expected_mode]]
    )
    require(references == expected,
            "%s has an alternate, conflicting, or missing diagnostic "
            "mode definition" % label)


def require_forbidden_experiment_macros_absent(
        arguments: Sequence[str], label: str) -> None:
    references = [
        argument for argument in arguments
        if any(name in argument for name in FORBIDDEN_EXPERIMENT_MACROS)
    ]
    require(not references,
            "%s enables or alters a forbidden experiment macro" % label)


def normalize_build_arguments(
        arguments: Sequence[str], build_root: Path,
        source_root: Path) -> list[str]:
    build_text = str(build_root)
    source_text = str(source_root)
    return [
        argument.replace(build_text, "<BUILD>").replace(
            source_text, "<SOURCE>")
        for argument in arguments
    ]


def require_compile_source_match(
        entry: Mapping[str, Any], arguments: Sequence[str],
        declared_source: Path) -> None:
    positions = [
        index for index, token in enumerate(arguments)
        if token == "-c"
    ]
    require(len(positions) == 1 and positions[0] + 1 < len(arguments),
            "compile command does not have one source operand")
    argv_source = Path(arguments[positions[0] + 1])
    if not argv_source.is_absolute():
        argv_source = Path(entry["directory"]) / argv_source
    require(argv_source.resolve(strict=True) == declared_source,
            "compile argv source differs from entry file")


def build_record_contract(
        build_schema: str) -> tuple[
            str, tuple[str, ...], tuple[str, ...], tuple[str, ...], str]:
    if build_schema == BUILD_SCHEMA:
        return (
            BUILD_CONFIGURATION_ATTESTATION_SCHEMA,
            BUILD_CONFIGURATION_VARIABLES,
            FORBIDDEN_EXPERIMENT_MACROS,
            CONTROLLED_CONFIGURATION_KEYS,
            CLEAN_REBUILD_SCHEMA,
        )
    if build_schema == BUILD_SCHEMA_V1:
        return (
            BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2,
            BUILD_CONFIGURATION_VARIABLES_V2,
            FORBIDDEN_EXPERIMENT_MACROS_V2,
            CONTROLLED_CONFIGURATION_KEYS_V1,
            CLEAN_REBUILD_SCHEMA_V1,
        )
    raise EvidenceError("retained exhaustive build schema is unsupported")


def build_configuration_digest(
        entries: Mapping[str, Any],
        attestation_schema: str = BUILD_CONFIGURATION_ATTESTATION_SCHEMA
) -> str:
    if attestation_schema == BUILD_CONFIGURATION_ATTESTATION_SCHEMA:
        variables = BUILD_CONFIGURATION_VARIABLES
    elif attestation_schema == BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2:
        variables = BUILD_CONFIGURATION_VARIABLES_V2
    else:
        raise EvidenceError(
            "effective configuration attestation schema is unsupported")
    require(isinstance(entries, dict) and
            set(entries) == set(variables) and
            all(isinstance(entries[name], str)
                for name in variables),
            "effective configuration entries have the wrong shape")
    material = "".join(
        "%s=%s\n" % (name, entries[name])
        for name in variables
    )
    return sha256_bytes(material.encode("utf-8"))


def validate_shard(
        value: Any, shard_index: int, shard_count: int,
        expected: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "mode", "shard_index", "shard_count", "total_patterns",
        "assigned_patterns", "recovered_shards", "verified_basis_symbols",
        "basis_seed", "assignment", "digest_fnv1a64",
        "ordinal_digest_fnv1a64",
    }, "shard result has the wrong shape")
    reconstructed = dict(expected)
    expected_value = {"schema": SHARD_SCHEMA, "mode": value.get("mode"),
                      **reconstructed}
    require(type(value.get("mode")) is int and value["mode"] in (1, 2) and
            value == expected_value and
            value["shard_index"] == shard_index and
            value["shard_count"] == shard_count,
            "shard result identity/counts/digests are invalid")
    return value


def load_module(path: Path, name: str) -> Any:
    specification = importlib.util.spec_from_file_location(name, path)
    require(specification is not None and specification.loader is not None,
            "cannot import provenance helper %s" % path)
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    try:
        specification.loader.exec_module(module)
    except Exception:
        sys.modules.pop(name, None)
        raise
    return module


def git_output(source: Path, *arguments: str) -> str:
    completed = subprocess.run(
        [str(CANONICAL_GIT), "-C", str(source), *arguments],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, env={
            "GIT_CONFIG_GLOBAL": "/dev/null", "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_NO_REPLACE_OBJECTS": "1", "GIT_OPTIONAL_LOCKS": "0",
            "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
        }, check=False)
    require(completed.returncode == 0,
            "git source identity failed: %s" %
            completed.stderr.decode(errors="replace"))
    return completed.stdout.decode("utf-8", errors="strict").strip()


def source_identity(source: Path) -> dict[str, Any]:
    status = git_output(source, "status", "--short", "--untracked-files=all")
    require(not status, "exhaustive evidence requires a clean source worktree")
    files = {}
    for relative in BOUND_SOURCE_FILES:
        path = source / relative
        require(path.is_file(), "bound source is missing %s" % relative)
        files[relative] = file_identity(path)
    result = {
        "root": str(source),
        "head": git_output(source, "rev-parse", "HEAD"),
        "head_tree": git_output(source, "rev-parse", "HEAD^{tree}"),
        "status_short": [],
        "files": files,
    }
    result["digest"] = record_digest(result)
    return result


def capture_build_closure(binary: Path) -> dict[str, Any]:
    source = Path(__file__).resolve(strict=True).parents[3]
    build = binary.resolve(strict=True).parent
    expected_binary = (build / TARGET_NAME).resolve(strict=True)
    require(binary.resolve(strict=True) == expected_binary,
            "binary is not the exact exhaustive CMake target")
    provenance_path = source / "tools/leopard2_build_provenance.py"
    configuration_path = source / "tools/leopard2_direct_encode_crossover.py"
    provenance = load_module(
        provenance_path, "leopard2_small_direct_exhaustive_provenance")
    configuration = load_module(
        configuration_path, "leopard2_small_direct_exhaustive_configuration")
    try:
        production = provenance.candidate_build_provenance(
            build, source, build / BENCHMARK_NAME, BENCHMARK_NAME)
        metadata = configuration.cmake_build_metadata(build / BENCHMARK_NAME)
        configuration.validate_build_configuration_attestation(
            metadata["effective_configuration_attestation"],
            build / configuration.BUILD_CONFIGURATION_RELATIVE_PATH)
        cache_identity, cache_bytes = provenance.file_snapshot(
            build / "CMakeCache.txt", "exhaustive CMake cache")
        cache = provenance.parse_cmake_cache(cache_bytes)
    except Exception as error:
        raise EvidenceError(
            "shared production provenance rejected exhaustive build: %s" %
            error) from error
    mode_text = cache.get(MODE_CACHE_VARIABLE)
    require(mode_text in ("1", "2"),
            "exhaustive target requires CMake small-direct mode 1 or 2")
    mode = int(mode_text)
    require(cache.get("LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN") == "OFF" and
            cache.get("LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE") == "OFF" and
            cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT") == "OFF",
            "canonical small-direct evidence requires source-plan, "
            "high-direct, and generalized one-loss experiments OFF")
    controlled_configuration = {
        name: cache.get(name) for name in CONTROLLED_CONFIGURATION_KEYS
    }
    require(all(isinstance(item, str)
                for item in controlled_configuration.values()),
            "canonical small-direct build omits a controlled cache entry")
    attested = metadata["effective_configuration_attestation"]
    require(
        isinstance(attested, dict) and
        attested.get("schema") == BUILD_CONFIGURATION_ATTESTATION_SCHEMA and
        isinstance(attested.get("entries"), dict) and
        attested["entries"].get(MODE_CACHE_VARIABLE) == mode_text and
        all(attested["entries"].get(name) == "OFF"
            for name in FORBIDDEN_EXPERIMENT_MACROS) and
        attested.get("sha256") ==
            build_configuration_digest(
                attested["entries"], BUILD_CONFIGURATION_ATTESTATION_SCHEMA),
        "current effective configuration has the wrong schema, selectors, "
        "or digest")

    commands_identity, commands_bytes = provenance.file_snapshot(
        build / "compile_commands.json",
        "exhaustive compile commands")
    try:
        compile_entries = decode_json_bytes(
            commands_bytes, "exhaustive compile database")
    except EvidenceError as error:
        raise EvidenceError("exhaustive compile database is invalid") from error
    require(isinstance(compile_entries, list) and compile_entries,
            "exhaustive compile database is empty")
    c_compiler = Path(cache["CMAKE_C_COMPILER"]).resolve(strict=True)
    cxx_compiler = Path(cache["CMAKE_CXX_COMPILER"]).resolve(strict=True)
    by_output = {}
    for entry in compile_entries:
        try:
            tokens = provenance._compile_tokens(entry)
            source_operand = Path(entry["file"]).resolve(strict=True)
            provenance._require_compile_driver(
                tokens, source_operand, c_compiler, cxx_compiler)
            require_compile_source_match(entry, tokens, source_operand)
            output = provenance._compile_output(entry, tokens, build)
        except Exception as error:
            raise EvidenceError(
                "exhaustive compile entry is invalid: %s" % error) from error
        require(output not in by_output,
                "exhaustive compile database has duplicate output")
        by_output[output] = (entry, tokens, source_operand)

    recipe_path = build / ("CMakeFiles/%s.dir/link.txt" % TARGET_NAME)
    recipe_identity, recipe_bytes = provenance.file_snapshot(
        recipe_path, "exhaustive link recipe")
    try:
        target_objects = provenance._recipe_objects(
            recipe_bytes, build, "exhaustive link recipe")
    except Exception as error:
        raise EvidenceError("exhaustive link recipe is invalid: %s" %
                            error) from error
    expected_sources = {
        (source / "tests/leopard2/direct_oracle.cpp").resolve(strict=True),
        (source /
         "tests/leopard2/test_small_direct_exhaustive.cpp").resolve(
             strict=True),
    }
    target_records = []
    found_sources = set()
    for obj in target_objects:
        require(obj in by_output,
                "exhaustive object has no exact compile command")
        entry, tokens, source_operand = by_output[obj]
        require(source_operand in expected_sources,
                "unexpected source in exhaustive target closure")
        require_forbidden_experiment_macros_absent(
            tokens, "exhaustive target source")
        if source_operand.name == "test_small_direct_exhaustive.cpp":
            require_mode_macro(
                tokens, mode, "exhaustive selector source")
        else:
            require_mode_macro(
                tokens, None, "mode-independent direct oracle")
        try:
            flag_profile = provenance._validate_compile_flags(
                tokens, source_operand)
        except Exception as error:
            raise EvidenceError(
                "exhaustive source has invalid compile flags: %s" %
                error) from error
        found_sources.add(source_operand)
        target_records.append({
            "source": provenance.file_identity(
                source_operand, "exhaustive source"),
            "object": provenance.file_identity(obj, "exhaustive object"),
            "compile_entry": {
                "directory": entry["directory"], "file": entry["file"],
                "output": entry["output"], "arguments": tokens,
            },
            "flag_profile": flag_profile,
        })
    require(found_sources == expected_sources,
            "exhaustive target source closure is incomplete")

    for record in production["source_object_compile_closure"]:
        arguments = record["compile_entry"]["arguments"]
        require_forbidden_experiment_macros_absent(
            arguments, "production source")
        selected_source = (
            record["role"] == "archive" and
            Path(record["source"]["path"]).resolve(strict=True) ==
                (source / "leopard2.cpp").resolve(strict=True))
        if selected_source:
            require_mode_macro(
                arguments, mode, "production leopard2.cpp compile")
        else:
            require_mode_macro(
                arguments, None, "unselected production source")

    try:
        recipe_text = recipe_bytes.decode("utf-8", errors="strict")
        link_tokens = []
        for line in recipe_text.splitlines():
            if line.strip():
                link_tokens.extend(shlex.split(line, posix=True))
    except (UnicodeError, ValueError) as error:
        raise EvidenceError("cannot parse exhaustive link recipe") from error
    archive = Path(production["archive"]["path"]).resolve(strict=True)
    archive_operands = []
    for token in link_tokens:
        operand = Path(token)
        if operand.name != archive.name:
            continue
        if not operand.is_absolute():
            operand = build / operand
        archive_operands.append(operand.resolve(strict=True))
    require(archive_operands == [archive],
            "exhaustive target does not link exactly the validated archive")
    exhaustive_identity = provenance.file_identity(
        expected_binary, "exhaustive executable")
    require(all(exhaustive_identity["mtime_ns"] >=
                record["object"]["mtime_ns"] for record in target_records) and
            exhaustive_identity["mtime_ns"] >=
                production["archive"]["mtime_ns"],
            "exhaustive executable predates a link input")
    sidecar_identity = file_identity(
        build / configuration.BUILD_CONFIGURATION_RELATIVE_PATH)
    result = {
        "schema": BUILD_SCHEMA,
        "mode": mode,
        "source": source_identity(source),
        "production_closure": production,
        "effective_configuration": metadata,
        "controlled_configuration": controlled_configuration,
        "cache": cache_identity,
        "compile_commands": commands_identity,
        "target_link_recipe": recipe_identity,
        "target_link_arguments": normalize_build_arguments(
            link_tokens, build, source),
        "target_compile_closure": sorted(
            target_records, key=lambda item: item["source"]["path"]),
        "mode_independent_object_sha256": next(
            record["object"]["sha256"] for record in target_records
            if Path(record["source"]["path"]).name ==
                "direct_oracle.cpp"),
        "exhaustive_executable": exhaustive_identity,
        "configuration_sidecar": sidecar_identity,
        "provenance_helper": file_identity(provenance_path),
        "configuration_reader": file_identity(configuration_path),
    }
    result["digest"] = record_digest(result)
    return result


def build_reproducible_material(
        value: Mapping[str, Any]) -> dict[str, Any]:
    source_root = Path(value["source"]["root"])
    production = value["production_closure"]
    build_root = Path(production["build_root"])

    def compact(identity: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "size": identity["size"], "sha256": identity["sha256"],
        }

    def relative_source(identity: Mapping[str, Any]) -> str:
        path = Path(identity["path"])
        require(path.is_relative_to(source_root),
                "reproducible source escapes the source root")
        return str(path.relative_to(source_root))

    production_objects = []
    for record in production["source_object_compile_closure"]:
        production_objects.append({
            "role": record["role"],
            "source": relative_source(record["source"]),
            "object": compact(record["object"]),
            "flag_profile": record["flag_profile"],
            "arguments": normalize_build_arguments(
                record["compile_entry"]["arguments"],
                build_root, source_root),
        })
    target_objects = []
    target_build_root = Path(value["cache"]["path"]).parent
    for record in value["target_compile_closure"]:
        target_objects.append({
            "source": relative_source(record["source"]),
            "object": compact(record["object"]),
            "flag_profile": record["flag_profile"],
            "arguments": normalize_build_arguments(
                record["compile_entry"]["arguments"],
                target_build_root, source_root),
        })
    return {
        "mode": value["mode"],
        "controlled_configuration": value["controlled_configuration"],
        "production_validated_cache": production["validated_cache"],
        "compiler": compact(production["compiler"]),
        "compiler_version_sha256":
            production["compiler_version_sha256"],
        "archiver": compact(production["archiver"]),
        "production_objects": sorted(
            production_objects,
            key=lambda item: (item["role"], item["source"])),
        "archive_members": production["archive_members"],
        "archive_member_identities":
            production["archive_member_identities"],
        "archive": compact(production["archive"]),
        "target_objects": sorted(
            target_objects, key=lambda item: item["source"]),
        "target_link_arguments": value["target_link_arguments"],
        "exhaustive_executable": compact(
            value["exhaustive_executable"]),
    }


def run_checked(
        command: Sequence[str], label: str, timeout: float = 1800.0,
        extra_environment: Mapping[str, str] | None = None) -> None:
    environment = {
        "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
        "TZ": "UTC",
    }
    if extra_environment is not None:
        environment.update(extra_environment)
    process: subprocess.Popen[bytes] | None = None
    process_support = None
    stderr = b""
    try:
        process_support, unused_support_path = \
            load_process_containment_support()
        with process_support.LinuxDescendantContainment() as containment:
            process = subprocess.Popen(
                list(command), stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                env=environment, start_new_session=True)
            containment.spawned_process = process
            containment.observe_spawn(process)
            containment.attach(process.pid)
            try:
                unused_stdout, stderr = process.communicate(timeout=timeout)
            finally:
                # A successful immediate child may still have daemonized work
                # into another session.  Prove the entire owned tree empty on
                # success, timeout, and every exceptional exit.
                containment.terminate_and_reap(process)
    except (OSError, subprocess.SubprocessError) as error:
        raise EvidenceError("%s failed to execute: %s" %
                            (label, error)) from error
    except Exception as error:
        if (process_support is not None and
                isinstance(error, process_support.EvidenceError)):
            raise EvidenceError("%s containment failed: %s" %
                                (label, error)) from error
        raise
    finally:
        if process is not None and process.stderr is not None:
            process.stderr.close()
    require(process is not None and process.returncode == 0,
            "%s failed with %d: %s" % (
                label, -1 if process is None else process.returncode,
                stderr[-8192:].decode(errors="replace")))


def controlled_configure_profile(
        controlled: Mapping[str, str], mode: int,
        build_schema: str = BUILD_SCHEMA) -> list[str]:
    unused_attestation, unused_variables, disabled_experiments, \
        unused_controlled, unused_rebuild = build_record_contract(build_schema)
    del unused_attestation, unused_variables, unused_controlled, unused_rebuild
    profile = [
        "-G", "Unix Makefiles",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DENABLE_OPENMP=ON",
        "-DLEOPARD_ENABLE_GF8=ON",
        "-DLEOPARD_ENABLE_GF16=ON",
        "-DLEO2_BACKEND_VARIANT=%s" %
            controlled["LEO2_BACKEND_VARIANT"],
        "-DLEO2_BUILD_ALLK_DIAGNOSTIC=OFF",
        "-DLEO2_BUILD_BENCHMARKS=ON",
        "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_BUILD_TESTS=ON",
        "-DLEO2_ENABLE_CUDA=OFF",
        "-DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF",
        "-DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF",
        "-D%s=%d" % (MODE_CACHE_VARIABLE, mode),
        "-DLEO2_FLAG_MAVX512F=FALSE",
        "-DLEO2_FLAG_MAVX512BW=FALSE",
        "-DLEO2_FLAG_MAVX512VL=FALSE",
    ]
    if "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT" in disabled_experiments:
        profile.insert(
            profile.index("-D%s=%d" % (MODE_CACHE_VARIABLE, mode)),
            "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF")
    return profile


def controlled_configure_environment(
        controlled: Mapping[str, str]) -> dict[str, str]:
    return {
        "CC": controlled["CMAKE_C_COMPILER"],
        "CXX": controlled["CMAKE_CXX_COMPILER"],
    }


def clean_rebuild_proof(
        candidate: Mapping[str, Any]) -> dict[str, Any]:
    controlled = candidate["controlled_configuration"]
    source = Path(candidate["source"]["root"])
    tool_keys = {
        "cmake": "CMAKE_COMMAND",
        "make": "CMAKE_MAKE_PROGRAM",
        "linker": "CMAKE_LINKER",
        "ranlib": "CMAKE_RANLIB",
        "archiver": "CMAKE_AR",
    }
    tools = {}
    for label, cache_key in tool_keys.items():
        path = Path(controlled[cache_key]).resolve(strict=True)
        require(path.is_file() and os.access(path, os.X_OK),
                "controlled rebuild tool is unavailable: %s" % path)
        tools[label] = {
            "configured_path": controlled[cache_key],
            "identity": file_identity(path),
        }

    configure_profile = controlled_configure_profile(
        controlled, candidate["mode"])
    configure_environment = controlled_configure_environment(controlled)
    with tempfile.TemporaryDirectory(
            prefix="leopard2-small-direct-rebuild-") as temporary:
        build = Path(temporary).resolve(strict=True)
        configure = [
            tools["cmake"]["identity"]["path"],
            "-S", str(source), "-B", str(build),
            *configure_profile,
        ]
        run_checked(
            configure, "controlled exhaustive configure",
            extra_environment=configure_environment)
        build_command = [
            tools["cmake"]["identity"]["path"],
            "--build", str(build), "--target",
            BENCHMARK_NAME, TARGET_NAME, "--parallel",
            str(min(128, len(os.sched_getaffinity(0)))),
        ]
        run_checked(
            build_command, "controlled exhaustive rebuild",
            extra_environment=configure_environment)
        rebuilt = capture_build_closure(build / TARGET_NAME)
        require(rebuilt["source"] == candidate["source"],
                "source identity changed during clean exhaustive rebuild")
        candidate_material = build_reproducible_material(candidate)
        rebuilt_material = build_reproducible_material(rebuilt)
        require(rebuilt_material == candidate_material,
                "clean exhaustive rebuild did not reproduce the "
                "archive/objects/link/executable bytes and semantics")
    result = {
        "schema": CLEAN_REBUILD_SCHEMA,
        "mode": candidate["mode"],
        "tools": tools,
        "configure_profile": configure_profile,
        "configure_environment": configure_environment,
        "candidate_material": candidate_material,
        "rebuilt_material": rebuilt_material,
    }
    result["digest"] = record_digest(result)
    return result


def capture_current_build(binary: Path) -> dict[str, Any]:
    result = capture_build_closure(binary)
    result["reproducible_rebuild"] = clean_rebuild_proof(result)
    result["digest"] = record_digest(result)
    return result


def build_without_rebuild_proof(
        value: Mapping[str, Any]) -> dict[str, Any]:
    result = dict(value)
    result.pop("reproducible_rebuild", None)
    result["digest"] = record_digest(result)
    return result


def validate_build_record(value: Any) -> None:
    build_keys = {
        "schema", "mode", "source", "production_closure",
        "effective_configuration", "controlled_configuration",
        "cache", "compile_commands",
        "target_link_recipe", "target_link_arguments",
        "target_compile_closure",
        "mode_independent_object_sha256", "exhaustive_executable",
        "configuration_sidecar", "provenance_helper",
        "configuration_reader", "reproducible_rebuild", "digest",
    }
    require(isinstance(value, dict) and set(value) == build_keys and
            value.get("schema") in (BUILD_SCHEMA_V1, BUILD_SCHEMA) and
            type(value.get("mode")) is int and
            value["mode"] in (1, 2) and
            value.get("digest") == record_digest(value) and
            isinstance(value.get("source"), dict) and
            isinstance(value.get("production_closure"), dict) and
            value["production_closure"].get("schema") ==
                "leopard2-production-build-closure/v1" and
            isinstance(value.get("target_compile_closure"), list) and
            len(value["target_compile_closure"]) == 2,
            "retained exhaustive build provenance is invalid")

    def valid_path(path: Any) -> bool:
        return (isinstance(path, str) and Path(path).is_absolute() and
                Path(path) == Path(os.path.abspath(path)))

    def valid_sha256(digest: Any) -> bool:
        return (isinstance(digest, str) and
                re.fullmatch(r"[0-9a-f]{64}", digest) is not None)

    def validate_minimal_identity(
            identity: Any, label: str) -> Mapping[str, Any]:
        require(isinstance(identity, dict) and set(identity) == {
                    "path", "size", "sha256",
                } and valid_path(identity.get("path")) and
                type(identity.get("size")) is int and
                identity["size"] >= 0 and
                valid_sha256(identity.get("sha256")),
                "%s identity is invalid" % label)
        return identity

    def validate_rich_identity(
            identity: Any, label: str) -> Mapping[str, Any]:
        require(isinstance(identity, dict) and set(identity) == {
                    "path", "sha256", "device", "inode", "mode", "size",
                    "mtime_ns", "ctime_ns",
                } and valid_path(identity.get("path")) and
                valid_sha256(identity.get("sha256")) and
                all(type(identity.get(field)) is int and
                    identity[field] >= 0 for field in (
                        "device", "inode", "mode", "size", "mtime_ns",
                        "ctime_ns")) and
                stat.S_ISREG(identity["mode"]),
                "%s identity is invalid" % label)
        return identity

    build_schema = value["schema"]
    (attestation_schema, configuration_variables, disabled_experiments,
     controlled_keys, clean_rebuild_schema) = \
        build_record_contract(build_schema)
    mode = value["mode"]
    source = value["source"]
    controlled = value["controlled_configuration"]
    fixed_configuration = {
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_GENERATOR": "Unix Makefiles",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_CXX_FLAGS": "",
        "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
        "ENABLE_OPENMP": "ON",
        "LEOPARD_ENABLE_GF8": "ON",
        "LEOPARD_ENABLE_GF16": "ON",
        "LEO2_BUILD_ALLK_DIAGNOSTIC": "OFF",
        "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_BUILD_TESTS": "ON",
        "LEO2_ENABLE_CUDA": "OFF",
        "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
        "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
        **({
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
        } if build_schema == BUILD_SCHEMA else {}),
        MODE_CACHE_VARIABLE: str(mode),
    }
    require(isinstance(controlled, dict) and
            set(controlled) == set(controlled_keys) and
            all(isinstance(item, str) for item in controlled.values()) and
            all(controlled[name] == expected
                for name, expected in fixed_configuration.items()) and
            controlled["LEO2_BACKEND_VARIANT"] in ("auto", "avx2") and
            all(Path(controlled[name]).is_absolute()
                for name in (
                    "CMAKE_AR", "CMAKE_COMMAND", "CMAKE_C_COMPILER",
                    "CMAKE_CXX_COMPILER", "CMAKE_LINKER",
                    "CMAKE_MAKE_PROGRAM", "CMAKE_RANLIB")) and
            controlled["LEO2_FLAG_MAVX2"] not in ("", "0", "FALSE", "OFF") and
            controlled["LEO2_FLAG_MNO_AVX512F"] not in
                ("", "0", "FALSE", "OFF") and
            all(controlled[name] in ("", "0", "FALSE", "OFF")
                for name in ("LEO2_FLAG_MAVX512F",
                             "LEO2_FLAG_MAVX512BW",
                             "LEO2_FLAG_MAVX512VL")),
            "retained controlled build configuration is invalid")
    require(set(source) == {
                "root", "head", "head_tree", "status_short", "files",
                "digest",
            } and
            isinstance(source.get("root"), str) and
            Path(source["root"]).is_absolute() and
            Path(source["root"]) ==
                Path(os.path.abspath(source["root"])) and
            isinstance(source.get("head"), str) and
            re.fullmatch(r"[0-9a-f]{40}", source["head"]) is not None and
            isinstance(source.get("head_tree"), str) and
            re.fullmatch(r"[0-9a-f]{40}",
                         source["head_tree"]) is not None and
            source.get("digest") == record_digest(source) and
            source.get("status_short") == [] and
            isinstance(source.get("files"), dict) and
            set(source["files"]) == set(BOUND_SOURCE_FILES),
            "retained source closure is invalid")
    source_root = Path(source["root"])
    for relative in BOUND_SOURCE_FILES:
        identity = validate_minimal_identity(
            source["files"][relative],
            "retained bound source %s" % relative)
        expected_path = source_root / relative
        require(identity["path"] == str(expected_path),
                "retained bound source identity is invalid: %s" % relative)

    rich_top_level = (
        ("cache", "retained CMake cache"),
        ("compile_commands", "retained compile commands"),
        ("target_link_recipe", "retained exhaustive link recipe"),
        ("exhaustive_executable", "retained exhaustive executable"),
    )
    for field, label in rich_top_level:
        validate_rich_identity(value[field], label)
    link_arguments = value["target_link_arguments"]
    require(isinstance(link_arguments, list) and link_arguments and
            all(isinstance(argument, str) and argument
                for argument in link_arguments) and
            sum(argument in ("libleopard.a", "<BUILD>/libleopard.a")
                for argument in link_arguments) == 1,
            "retained exhaustive target link arguments are invalid")
    minimal_top_level = (
        ("configuration_sidecar", "retained configuration sidecar"),
        ("provenance_helper", "retained provenance helper"),
        ("configuration_reader", "retained configuration reader"),
    )
    for field, label in minimal_top_level:
        validate_minimal_identity(value[field], label)

    production = value["production_closure"]
    production_cache = production.get("validated_cache")
    general = "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"
    require(isinstance(production_cache, dict) and
            production_cache.get(MODE_CACHE_VARIABLE) == str(mode) and
            all(production_cache.get(name) == "OFF"
                for name in disabled_experiments) and
            ((build_schema == BUILD_SCHEMA and
              production_cache.get(general) == "OFF") or
             (build_schema == BUILD_SCHEMA_V1 and
              general not in production_cache)),
            "retained production cache has the wrong versioned selectors")
    records = production.get("source_object_compile_closure")
    require(isinstance(records, list) and records,
            "retained production compile closure is empty")
    for record in records:
        require(isinstance(record, dict) and
                record.get("role") in ("archive", "benchmark") and
                isinstance(record.get("compile_entry"), dict) and
                isinstance(record["compile_entry"].get("arguments"), list),
                "retained production compile record is invalid")
        arguments = record["compile_entry"]["arguments"]
        require_forbidden_experiment_macros_absent(
            arguments, "retained production source")
        selected_source = (
            record["role"] == "archive" and
            record.get("source", {}).get("path") ==
                str(Path(source["root"]) / "leopard2.cpp"))
        if selected_source:
            require_mode_macro(
                arguments, mode, "retained production leopard2.cpp")
        else:
            require_mode_macro(
                arguments, None, "retained unselected production source")
    targets = value["target_compile_closure"]
    expected_names = {
        "direct_oracle.cpp", "test_small_direct_exhaustive.cpp",
    }
    require(all(isinstance(item, dict) and
                isinstance(item.get("source"), dict) and
                isinstance(item["source"].get("path"), str)
                for item in targets) and
            {Path(item["source"]["path"]).name
             for item in targets} == expected_names,
            "retained exhaustive target source closure changed")
    for record in targets:
        require(isinstance(record, dict) and set(record) == {
                    "source", "object", "compile_entry", "flag_profile",
                } and isinstance(record["flag_profile"], str) and
                bool(record["flag_profile"]),
                "retained exhaustive target compile record is invalid")
        target_source = validate_rich_identity(
            record["source"], "retained exhaustive target source")
        target_object = validate_rich_identity(
            record["object"], "retained exhaustive target object")
        compile_entry = record["compile_entry"]
        require(isinstance(compile_entry, dict) and set(compile_entry) == {
                    "directory", "file", "output", "arguments",
                } and valid_path(compile_entry.get("directory")) and
                isinstance(compile_entry.get("file"), str) and
                bool(compile_entry["file"]) and
                isinstance(compile_entry.get("output"), str) and
                bool(compile_entry["output"]) and
                isinstance(compile_entry.get("arguments"), list) and
                compile_entry["arguments"] and
                all(isinstance(argument, str) and argument
                    for argument in compile_entry["arguments"]),
                "retained exhaustive target compile entry is invalid")
        compile_file = Path(compile_entry["file"])
        compile_output = Path(compile_entry["output"])
        directory = Path(compile_entry["directory"])
        if not compile_file.is_absolute():
            compile_file = directory / compile_file
        if not compile_output.is_absolute():
            compile_output = directory / compile_output
        require(Path(os.path.abspath(compile_file)) ==
                    Path(target_source["path"]) and
                Path(os.path.abspath(compile_output)) ==
                    Path(target_object["path"]),
                "retained exhaustive target compile paths changed")
        target_source_path = Path(target_source["path"])
        require(target_source_path.is_relative_to(source_root),
                "retained exhaustive target source escapes the "
                "bound source closure")
        source_relative = target_source_path.relative_to(source_root)
        bound_source = source["files"].get(str(source_relative))
        require(bound_source is not None and
                all(target_source[field] == bound_source[field]
                    for field in ("path", "size", "sha256")),
                "retained exhaustive target source is outside the "
                "bound source closure")
        arguments = compile_entry["arguments"]
        require_forbidden_experiment_macros_absent(
            arguments, "retained exhaustive target source")
        name = Path(record["source"]["path"]).name
        if name == "test_small_direct_exhaustive.cpp":
            require_mode_macro(
                arguments, mode, "retained exhaustive selector compile")
        else:
            require_mode_macro(
                arguments, None,
                "retained mode-independent direct oracle")
            require(valid_sha256(
                        value.get("mode_independent_object_sha256")) and
                    value["mode_independent_object_sha256"] ==
                        target_object["sha256"],
                    "retained direct-oracle object identity changed")

    metadata = value["effective_configuration"]
    require(isinstance(metadata, dict) and set(metadata) == {
                "binding_scope", "build_root", "cmake_cache",
                "cmake_cache_sha256",
                "effective_configuration_attestation", "entries",
                "executable", "extra_file_sha256",
            } and isinstance(metadata["binding_scope"], str) and
            bool(metadata["binding_scope"]) and
            valid_path(metadata["build_root"]) and
            valid_path(metadata["cmake_cache"]) and
            valid_sha256(metadata["cmake_cache_sha256"]) and
            isinstance(metadata["entries"], dict) and
            all(isinstance(key, str) and isinstance(item, str)
                for key, item in metadata["entries"].items()) and
            isinstance(metadata["extra_file_sha256"], dict) and
            all(isinstance(key, str) and valid_sha256(item)
                for key, item in metadata["extra_file_sha256"].items()) and
            metadata["cmake_cache"] == value["cache"]["path"] and
            metadata["cmake_cache_sha256"] == value["cache"]["sha256"],
            "retained effective configuration metadata is invalid")
    effective_executable = metadata["executable"]
    require(isinstance(effective_executable, dict) and
            set(effective_executable) == {
                "mtime_ns", "path", "sha256", "size",
            } and type(effective_executable["mtime_ns"]) is int and
            effective_executable["mtime_ns"] >= 0 and
            valid_path(effective_executable["path"]) and
            valid_sha256(effective_executable["sha256"]) and
            type(effective_executable["size"]) is int and
            effective_executable["size"] >= 0,
            "retained effective-configuration executable is invalid")
    attestation = metadata["effective_configuration_attestation"]
    require(isinstance(attestation, dict) and set(attestation) == {
                "schema", "entries", "path", "sha256",
            } and attestation.get("schema") ==
                attestation_schema and
            isinstance(attestation.get("entries"), dict) and
            set(attestation["entries"]) == set(configuration_variables) and
            all(isinstance(attestation["entries"][name], str)
                for name in configuration_variables) and
            attestation["entries"].get(MODE_CACHE_VARIABLE) == str(mode) and
            all(attestation["entries"].get(name) == "OFF"
                for name in disabled_experiments) and
            attestation.get("sha256") ==
                build_configuration_digest(
                    attestation["entries"], attestation_schema) and
            metadata["entries"].get(MODE_CACHE_VARIABLE) == str(mode) and
            all(metadata["entries"].get(name) == "OFF"
                for name in disabled_experiments) and
            ((build_schema == BUILD_SCHEMA and
              metadata["entries"].get(general) == "OFF") or
             (build_schema == BUILD_SCHEMA_V1 and
              general not in metadata["entries"])) and
            valid_path(attestation.get("path")) and
            valid_sha256(attestation.get("sha256")) and
            attestation["path"] == value["configuration_sidecar"]["path"] and
            value["configuration_sidecar"]["size"] > 0,
            "retained effective configuration has the wrong mode/schema")

    proof = value["reproducible_rebuild"]
    require(isinstance(proof, dict) and set(proof) == {
                "schema", "mode", "tools", "configure_profile",
                "configure_environment",
                "candidate_material", "rebuilt_material", "digest",
            } and
            proof["schema"] == clean_rebuild_schema and
            type(proof["mode"]) is int and proof["mode"] == mode and
            proof["digest"] == record_digest(proof) and
            proof["configure_profile"] ==
                controlled_configure_profile(
                    controlled, mode, build_schema) and
            proof["configure_environment"] ==
                controlled_configure_environment(controlled) and
            isinstance(proof["tools"], dict) and
            set(proof["tools"]) ==
                {"cmake", "make", "linker", "ranlib", "archiver"},
            "retained clean-rebuild proof is invalid")
    tool_keys = {
        "cmake": "CMAKE_COMMAND",
        "make": "CMAKE_MAKE_PROGRAM",
        "linker": "CMAKE_LINKER",
        "ranlib": "CMAKE_RANLIB",
        "archiver": "CMAKE_AR",
    }
    for label, cache_key in tool_keys.items():
        tool = proof["tools"][label]
        require(isinstance(tool, dict) and set(tool) == {
                    "configured_path", "identity",
                } and tool["configured_path"] == controlled[cache_key],
                "retained rebuild tool configuration changed")
        validate_minimal_identity(
            tool["identity"], "retained rebuild %s" % label)
    try:
        expected_material = build_reproducible_material(value)
    except (KeyError, TypeError, ValueError) as error:
        raise EvidenceError(
            "retained reproducible build material is malformed") from error
    require(proof["candidate_material"] == expected_material and
            proof["rebuilt_material"] == expected_material,
            "retained clean rebuild does not reproduce the exact build")


def iter_path_identities(value: Any) -> Iterable[dict[str, Any]]:
    if isinstance(value, dict):
        if (isinstance(value.get("path"), str) and
                type(value.get("size")) is int and
                isinstance(value.get("sha256"), str) and
                re.fullmatch(r"[0-9a-f]{64}", value["sha256"])):
            yield {
                "path": value["path"], "size": value["size"],
                "sha256": value["sha256"],
            }
        for child in value.values():
            yield from iter_path_identities(child)
    elif isinstance(value, list):
        for child in value:
            yield from iter_path_identities(child)


def freeze_file(
        origin: Mapping[str, Any], destination: Path, mode: int,
        campaign_root: Path) -> dict[str, Any]:
    source = Path(origin["path"])
    expected = {
        "path": str(source.resolve(strict=True)),
        "size": origin["size"], "sha256": origin["sha256"],
    }
    source_identity, source_payload = file_snapshot(source)
    require(source_identity == expected, "bundle origin identity changed")
    if destination.exists():
        canonical_artifact(destination, campaign_root,
                           "retained frozen file")
        actual = file_identity(destination)
        require(actual["size"] == origin["size"] and
                actual["sha256"] == origin["sha256"] and
                stat.S_IMODE(destination.lstat().st_mode) == mode,
                "retained frozen file changed")
    else:
        atomic_bytes(destination, source_payload, MAX_FILE_BYTES, mode)
    exact_file_identity(source, expected, "bundle origin")
    frozen = file_identity(destination)
    require(frozen["size"] == origin["size"] and
            frozen["sha256"] == origin["sha256"],
            "frozen artifact copy changed bytes")
    return {
        "origin": dict(origin), "frozen": frozen, "mode": mode,
    }


def freeze_bundle(
        build_record: Mapping[str, Any], extra_identities: Sequence[
            Mapping[str, Any]], output: Path) -> dict[str, Any]:
    root = output.resolve(strict=True)
    directory = ensure_directory(output / "artifacts" / "provenance", root)
    identities = {}
    for identity in itertools.chain(
            iter_path_identities(build_record), extra_identities):
        path = str(Path(identity["path"]).resolve(strict=True))
        retained = {
            "path": path, "size": identity["size"],
            "sha256": identity["sha256"],
        }
        if path in identities:
            require(identities[path] == retained,
                    "one provenance path has conflicting identities")
        identities[path] = retained
    entries = []
    for index, path in enumerate(sorted(identities)):
        destination = directory / ("file-%04d" % index)
        entries.append(freeze_file(
            identities[path], destination, 0o444, root))
    result = {
        "schema": BUNDLE_SCHEMA,
        "entries": entries,
        "origin_count": len(entries),
    }
    result["digest"] = record_digest(result)
    return result


def validate_bundle(
        value: Any, build_record: Mapping[str, Any],
        extra_identities: Sequence[Mapping[str, Any]], campaign_root: Path
) -> None:
    require(isinstance(value, dict) and value.get("schema") == BUNDLE_SCHEMA and
            value.get("digest") == record_digest(value) and
            type(value.get("origin_count")) is int and
            isinstance(value.get("entries"), list) and
            value["origin_count"] == len(value["entries"]),
            "retained provenance bundle header is invalid")
    expected_origins = {}
    for identity in itertools.chain(
            iter_path_identities(build_record), extra_identities):
        path = str(Path(identity["path"]))
        minimal = {
            "path": path, "size": identity["size"],
            "sha256": identity["sha256"],
        }
        if path in expected_origins:
            require(expected_origins[path] == minimal,
                    "retained build has conflicting path identities")
        expected_origins[path] = minimal
    require(len(expected_origins) == len(value["entries"]),
            "provenance bundle omits or duplicates an origin")
    seen = set()
    directory = campaign_root / "artifacts" / "provenance"
    canonical_directory(directory, campaign_root,
                        "provenance bundle directory")
    for index, entry in enumerate(value["entries"]):
        require(isinstance(entry, dict) and
                set(entry) == {"origin", "frozen", "mode"} and
                entry["mode"] == 0o444,
                "provenance bundle entry is malformed")
        origin = entry["origin"]
        require(origin in expected_origins.values() and
                origin["path"] not in seen,
                "provenance bundle origin is unknown or duplicated")
        seen.add(origin["path"])
        expected_path = directory / ("file-%04d" % index)
        frozen = entry["frozen"]
        require(frozen.get("path") == str(expected_path) and
                frozen.get("size") == origin["size"] and
                frozen.get("sha256") == origin["sha256"],
                "frozen provenance entry is not bound to its origin")
        canonical_artifact(expected_path, campaign_root,
                           "frozen provenance entry")
        exact_file_identity(expected_path, frozen,
                            "frozen provenance entry")
        require(stat.S_IMODE(expected_path.lstat().st_mode) == 0o444,
                "frozen provenance entry mode changed")
    require(seen == set(expected_origins),
            "provenance bundle origin set is incomplete")
    actual_names = {path.name for path in directory.iterdir()}
    require(actual_names == {
                "file-%04d" % index for index in range(len(value["entries"]))
            },
            "provenance bundle contains extra or missing files")


def retained_shard(
        root: Path, shard_index: int, request: Mapping[str, Any],
        expected: Mapping[str, Any]) -> dict[str, Any] | None:
    shard_count = request["shard_count"]
    directory = root / ("shard-%04d" % shard_index)
    result_path = directory / "result.json"
    stdout_path = directory / "stdout.json"
    stderr_path = directory / "stderr.txt"
    if not result_path.exists():
        return None
    for path, label in ((result_path, "result"), (stdout_path, "stdout"),
                        (stderr_path, "stderr")):
        canonical_artifact(path, root, "retained shard %s" % label)
    record = load_json(result_path, MAX_SHARD_JSON_BYTES)
    require(isinstance(record, dict) and set(record) == {
        "shard", "stdout", "stderr", "binary", "cpu", "command",
        "execution",
    }, "retained shard envelope has the wrong shape")
    cpu = request["allowed_cpus"][shard_index]
    command = [
        request["taskset"]["path"], "-c", str(cpu),
        request["frozen_binary"]["path"],
        "--shard-index", str(shard_index),
        "--shard-count", str(shard_count),
    ]
    execution = {
        "policy": request["execution_policy"],
        "interpreter": request["interpreter"],
        "binary_snapshot": sealed_snapshot_identity(
            request["frozen_binary"], 0o500),
        "runner_snapshot": sealed_snapshot_identity(
            request["runner"], 0o400),
        "taskset_snapshot": sealed_snapshot_identity(
            request["taskset"], 0o500),
    }
    require(record["stdout"] == file_identity(stdout_path) and
            record["stderr"] == file_identity(stderr_path) and
            record["binary"] == request["frozen_binary"] and
            record["cpu"] == cpu and record["command"] == command and
            record["execution"] == execution,
            "retained shard stream/binary/command/execution identity changed")
    parsed = load_json(stdout_path, MAX_SHARD_JSON_BYTES)
    require(record["shard"] == parsed,
            "retained shard result differs from stdout")
    validated = validate_shard(
        parsed, shard_index, shard_count, expected)
    require(validated["mode"] == request["build_provenance"]["mode"],
            "retained shard mode differs from build provenance")
    return validated


def sealed_snapshot_identity(
        identity: Mapping[str, Any], mode: int) -> dict[str, Any]:
    return {
        "schema": "leopard2-linux-sealed-snapshot/v1",
        "sha256": identity["sha256"], "size": identity["size"],
        "mode": mode,
        "seals": (
            LINUX_F_SEAL_SEAL | LINUX_F_SEAL_SHRINK |
            LINUX_F_SEAL_GROW | LINUX_F_SEAL_WRITE),
    }


def linux_memfd_create(name: str, flags: int) -> int:
    """Call memfd_create through libc when CPython has no os wrapper."""
    wrapper = getattr(os, "memfd_create", None)
    if callable(wrapper):
        return wrapper(name, flags)
    require(sys.platform.startswith("linux"),
            "immutable shard execution requires Linux memfd_create")
    try:
        function = ctypes.CDLL(None, use_errno=True).memfd_create
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            "immutable shard execution requires Linux memfd_create") \
            from error
    function.argtypes = (ctypes.c_char_p, ctypes.c_uint)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    descriptor = function(name.encode("utf-8"), ctypes.c_uint(flags))
    if descriptor >= 0:
        return descriptor
    number = ctypes.get_errno()
    raise OSError(number or errno.EPERM,
                  os.strerror(number or errno.EPERM))


def sealed_file_snapshot(
        path: Path, expected: Mapping[str, Any],
        label: str, mode: int) -> tuple[int, dict[str, Any]]:
    """Copy exact bytes into an immutable Linux memfd used for execution."""
    require(sys.platform.startswith("linux"),
            "immutable shard execution requires Linux memfd_create")
    identity, payload = file_snapshot(path)
    require(identity == expected, "%s identity changed" % label)
    descriptor = linux_memfd_create(
        "leopard2-" + re.sub(r"[^A-Za-z0-9_.-]", "-", label),
        getattr(os, "MFD_CLOEXEC", 0x0001) |
        getattr(os, "MFD_ALLOW_SEALING", 0x0002))
    try:
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            require(written > 0, "short immutable snapshot write")
            offset += written
        os.fchmod(descriptor, mode)
        seals = (LINUX_F_SEAL_SEAL | LINUX_F_SEAL_SHRINK |
                 LINUX_F_SEAL_GROW | LINUX_F_SEAL_WRITE)
        fcntl.fcntl(descriptor, LINUX_F_ADD_SEALS, seals)
        retained = fcntl.fcntl(descriptor, LINUX_F_GET_SEALS)
        status = os.fstat(descriptor)
        require(retained & seals == seals and status.st_size == len(payload),
                "immutable snapshot seals or size changed")
        os.lseek(descriptor, 0, os.SEEK_SET)
        exact_file_identity(path, expected, label)
        snapshot = sealed_snapshot_identity(identity, mode)
        require(snapshot["sha256"] == sha256_bytes(payload) and
                snapshot["size"] == len(payload) and
                snapshot["seals"] == retained,
                "immutable snapshot record changed")
        return descriptor, snapshot
    except BaseException:
        os.close(descriptor)
        raise


def validate_sealed_descriptor(
        descriptor: int, expected_size: int, expected_sha256: str,
        expected_mode: int, label: str) -> None:
    """Validate an inherited sealed memfd without changing its file offset."""
    require(type(descriptor) is int and descriptor >= 0 and
            type(expected_size) is int and 0 <= expected_size <=
                MAX_FILE_BYTES and
            isinstance(expected_sha256, str) and
            re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is not None,
            "%s immutable descriptor arguments are invalid" % label)
    status = os.fstat(descriptor)
    seals = fcntl.fcntl(descriptor, LINUX_F_GET_SEALS)
    required_seals = (
        LINUX_F_SEAL_SEAL | LINUX_F_SEAL_SHRINK |
        LINUX_F_SEAL_GROW | LINUX_F_SEAL_WRITE)
    require(stat.S_ISREG(status.st_mode) and
            stat.S_IMODE(status.st_mode) == expected_mode and
            status.st_size == expected_size and
            seals & required_seals == required_seals,
            "%s immutable descriptor metadata changed" % label)
    digest = hashlib.sha256()
    offset = 0
    while offset < expected_size:
        block = os.pread(
            descriptor, min(1 << 20, expected_size - offset), offset)
        require(block, "%s immutable descriptor was truncated" % label)
        digest.update(block)
        offset += len(block)
    require(os.pread(descriptor, 1, expected_size) == b"" and
            digest.hexdigest() == expected_sha256,
            "%s immutable descriptor bytes changed" % label)


def child_subreaper_state() -> int:
    """Return the exact Linux child-subreaper state for this process."""
    library = ctypes.CDLL(None, use_errno=True)
    value = ctypes.c_int(-1)
    result = library.prctl(
        LINUX_PR_GET_CHILD_SUBREAPER, ctypes.byref(value), 0, 0, 0)
    if result != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    require(value.value in (0, 1),
            "Linux returned an invalid child-subreaper state")
    return value.value


def arm_parent_death_continue(expected_parent_pid: int) -> None:
    """Wake a stopped supervisor when its exact coordinator disappears."""
    require(type(expected_parent_pid) is int and expected_parent_pid > 1 and
            os.getppid() == expected_parent_pid,
            "shard supervisor coordinator identity changed before startup")
    library = ctypes.CDLL(None, use_errno=True)
    result = library.prctl(
        LINUX_PR_SET_PDEATHSIG, int(signal.SIGCONT), 0, 0, 0)
    if result != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    require(os.getppid() == expected_parent_pid,
            "shard supervisor coordinator exited during startup")


def set_child_subreaper(value: int = 1) -> None:
    """Set and verify child-subreaper ownership for exact descendant reaping."""
    require(value in (0, 1), "invalid child-subreaper state")
    library = ctypes.CDLL(None, use_errno=True)
    result = library.prctl(
        LINUX_PR_SET_CHILD_SUBREAPER, value, 0, 0, 0)
    if result != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    require(child_subreaper_state() == value,
            "Linux did not apply the requested child-subreaper state")


def adopted_child_pids() -> list[int]:
    """Return exact direct children adopted by this single-purpose process."""
    task_root = Path("/proc/%d/task" % os.getpid())
    result: set[int] = set()
    for task in tuple(task_root.iterdir()):
        if not task.name.isdigit():
            continue
        try:
            descriptor = os.open(
                task / "children",
                os.O_RDONLY | os.O_CLOEXEC |
                getattr(os, "O_NOFOLLOW", 0) |
                getattr(os, "O_NONBLOCK", 0))
        except FileNotFoundError:
            continue
        try:
            payload = os.read(descriptor, 1 << 20)
            require(not os.read(descriptor, 1),
                    "adopted descendant inventory exceeds its byte bound")
        finally:
            os.close(descriptor)
        try:
            result.update(int(item) for item in payload.split())
        except ValueError as error:
            raise EvidenceError(
                "adopted descendant inventory is malformed") from error
    require(all(pid > 0 and pid != os.getpid() for pid in result),
            "adopted descendant inventory contains an invalid PID")
    return sorted(result)


def linux_pidfd_open(pid: int) -> int:
    wrapper = getattr(os, "pidfd_open", None)
    if callable(wrapper):
        return wrapper(pid, 0)
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_open
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            "pidfd descendant containment is unavailable") from error
    function.argtypes = (ctypes.c_int, ctypes.c_uint)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    descriptor = function(ctypes.c_int(pid), ctypes.c_uint(0))
    if descriptor >= 0:
        return descriptor
    number = ctypes.get_errno()
    if number == errno.ESRCH:
        raise ProcessLookupError(number, os.strerror(number))
    raise OSError(number or errno.EPERM,
                  os.strerror(number or errno.EPERM))


def linux_pidfd_send_signal(
        descriptor: int, signal_number: int) -> None:
    wrapper = getattr(signal, "pidfd_send_signal", None)
    if callable(wrapper):
        wrapper(descriptor, signal_number, None, 0)
        return
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_send_signal
    except (AttributeError, OSError) as error:
        raise EvidenceError(
            "pidfd descendant containment is unavailable") from error
    function.argtypes = (
        ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_uint)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    result = function(
        ctypes.c_int(descriptor), ctypes.c_int(signal_number),
        None, ctypes.c_uint(0))
    if result == 0:
        return
    number = ctypes.get_errno()
    if number == errno.ESRCH:
        raise ProcessLookupError(number, os.strerror(number))
    raise OSError(number or errno.EPERM,
                  os.strerror(number or errno.EPERM))


def pidfd_signal(pid: int, signal_number: int) -> None:
    """Signal one exact process identity without a PID-reuse race."""
    require(isinstance(signal_number, int) and
            not isinstance(signal_number, bool) and signal_number >= 0,
            "pidfd signal number is invalid")
    try:
        descriptor = linux_pidfd_open(pid)
    except ProcessLookupError:
        return
    try:
        try:
            linux_pidfd_send_signal(descriptor, signal_number)
        except ProcessLookupError:
            pass
    finally:
        os.close(descriptor)


def pidfd_kill(pid: int) -> None:
    pidfd_signal(pid, signal.SIGKILL)


def validate_pidfd_support() -> None:
    descriptor = linux_pidfd_open(os.getpid())
    try:
        linux_pidfd_send_signal(descriptor, 0)
    finally:
        os.close(descriptor)


def kill_and_reap_descendants(
        process: subprocess.Popen[bytes],
        timeout: float = CGROUP_TEARDOWN_SECONDS) -> None:
    """Kill/reap the shard tree even when members escaped their cgroup."""
    deadline = time.monotonic() + timeout
    if process.poll() is None:
        pidfd_kill(process.pid)
        try:
            process.wait(timeout=max(0.001, deadline - time.monotonic()))
        except subprocess.TimeoutExpired as error:
            raise EvidenceError(
                "contained shard leader could not be killed/reaped") \
                from error
    while True:
        for pid in adopted_child_pids():
            pidfd_kill(pid)
        try:
            pid, unused_status = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            return
        if pid > 0:
            continue
        if time.monotonic() >= deadline:
            raise EvidenceError(
                "adopted shard descendants could not be killed/reaped")
        time.sleep(0.01)


def kill_and_reap_campaign_children(
        timeout: float = CGROUP_TEARDOWN_SECONDS) -> None:
    """Parent fallback for supervisors or descendants that failed cleanup.

    The campaign coordinator enables subreaper mode only after proving it has
    no pre-existing direct children.  Consequently every child visible here
    belongs to this launch, including a shard descendant that migrated out of
    its delegated cgroup before its supervisor was forcibly terminated.
    """
    deadline = time.monotonic() + timeout
    empty_scans = 0
    while True:
        children = adopted_child_pids()
        for pid in children:
            pidfd_kill(pid)
        while True:
            try:
                pid, unused_status = os.waitpid(-1, os.WNOHANG)
            except ChildProcessError:
                pid = 0
                break
            if pid <= 0:
                break
        remaining = adopted_child_pids()
        if not remaining:
            empty_scans += 1
            if empty_scans >= 2:
                return
        else:
            empty_scans = 0
        if time.monotonic() >= deadline:
            raise EvidenceError(
                "campaign descendants could not be killed/reaped: %s" %
                remaining)
        time.sleep(0.01)


def cgroup_v2_current() -> Path:
    require(sys.platform.startswith("linux"),
            "shard containment requires Linux cgroup v2")
    records = [
        line.split(":", 2) for line in
        Path("/proc/self/cgroup").read_text(encoding="ascii").splitlines()
    ]
    unified = [record[2] for record in records
               if len(record) == 3 and record[0] == "0" and record[1] == ""]
    require(len(unified) == 1 and unified[0].startswith("/"),
            "unified cgroup-v2 membership is unavailable")
    root = Path("/sys/fs/cgroup").resolve(strict=True)
    current = (root / unified[0].lstrip("/")).resolve(strict=True)
    require(current == root or root in current.parents,
            "current cgroup-v2 membership escapes its mount")
    for name in ("cgroup.procs", "cgroup.events", "cgroup.kill"):
        require((current / name).exists(),
                "cgroup-v2 containment omits %s" % name)
    require(os.access(current, os.W_OK | os.X_OK),
            "current cgroup-v2 scope is not delegated to this runner")
    return current


def cgroup_write(path: Path, payload: bytes) -> None:
    descriptor = os.open(
        path, os.O_WRONLY | os.O_CLOEXEC |
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))
    try:
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            require(written > 0, "short cgroup control write")
            offset += written
    finally:
        os.close(descriptor)


def cgroup_pids(scope: Path) -> list[int]:
    descriptor = os.open(
        scope / "cgroup.procs",
        os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0) |
        getattr(os, "O_NONBLOCK", 0))
    try:
        payload = os.read(descriptor, 1 << 20)
        require(not os.read(descriptor, 1),
                "cgroup process inventory exceeds its byte bound")
    finally:
        os.close(descriptor)
    try:
        return [int(line) for line in payload.splitlines() if line]
    except ValueError as error:
        raise EvidenceError("cgroup process inventory is malformed") from error


def signal_cgroup(scope: Path) -> None:
    if not scope.exists():
        return
    try:
        cgroup_write(scope / "cgroup.kill", b"1")
    except FileNotFoundError:
        return


def cgroup_subtree_directories(scope: Path) -> list[Path]:
    """Return a validated post-order inventory rooted at one cgroup.

    Shard code runs under the same UID as the campaign and can create child
    cgroups inside its delegated scope.  Never follow a link or walk above the
    exact supplied root while inventorying those descendants.
    """
    root = Path("/sys/fs/cgroup").resolve(strict=True)
    require(scope.is_absolute() and
            (scope == root or root in scope.parents),
            "contained cgroup subtree escapes its mount")
    try:
        linked = scope.lstat()
        resolved = scope.resolve(strict=True)
    except FileNotFoundError:
        return []
    require(stat.S_ISDIR(linked.st_mode) and resolved == scope,
            "contained cgroup subtree root changed")

    postorder: list[Path] = []
    stack: list[tuple[Path, bool]] = [(scope, False)]
    visited: set[tuple[int, int]] = set()
    while stack:
        current, expanded = stack.pop()
        try:
            metadata = current.lstat()
            resolved_current = current.resolve(strict=True)
        except FileNotFoundError:
            continue
        require(stat.S_ISDIR(metadata.st_mode) and
                resolved_current == current and
                (current == scope or scope in current.parents),
                "contained cgroup subtree entry escaped its root")
        identity = (metadata.st_dev, metadata.st_ino)
        if expanded:
            postorder.append(current)
            continue
        require(identity not in visited,
                "contained cgroup subtree contains a directory cycle")
        visited.add(identity)
        stack.append((current, True))
        try:
            entries = tuple(current.iterdir())
        except FileNotFoundError:
            continue
        children = []
        for child in entries:
            try:
                child_metadata = child.lstat()
            except FileNotFoundError:
                continue
            if not stat.S_ISDIR(child_metadata.st_mode):
                continue
            try:
                resolved_child = child.resolve(strict=True)
            except FileNotFoundError:
                continue
            require(child.parent == current and resolved_child == child and
                    scope in child.parents,
                    "contained cgroup child escaped its root")
            children.append(child)
        for child in sorted(children, reverse=True):
            stack.append((child, False))
    return postorder


def cgroup_subtree_pids(scope: Path) -> list[int]:
    """Return every process in an exact cgroup subtree."""
    result: set[int] = set()
    for current in cgroup_subtree_directories(scope):
        try:
            result.update(cgroup_pids(current))
        except FileNotFoundError:
            continue
    return sorted(result)


def remove_empty_cgroup(
        scope: Path, timeout: float = CGROUP_TEARDOWN_SECONDS) -> None:
    """Remove an empty cgroup and every empty descendant, deepest first."""
    require(type(timeout) in (int, float) and math.isfinite(timeout) and
            timeout >= 0,
            "contained cgroup teardown timeout is invalid")
    deadline = time.monotonic() + timeout
    while scope.exists():
        directories = cgroup_subtree_directories(scope)
        if not directories:
            return
        residual_pids: set[int] = set()
        residual_directories: list[str] = []
        for current in directories:
            try:
                pids = cgroup_pids(current)
            except FileNotFoundError:
                continue
            residual_pids.update(pids)
            if pids:
                residual_directories.append(str(current))
                continue
            try:
                current.rmdir()
            except FileNotFoundError:
                continue
            except OSError as error:
                if error.errno not in (
                        errno.EBUSY, errno.ENOTEMPTY, errno.ENOENT):
                    raise
                residual_directories.append(str(current))
        if not scope.exists():
            return
        if time.monotonic() >= deadline:
            raise EvidenceError(
                "contained cgroup subtree remained after teardown: "
                "pids=%s directories=%s" % (
                    sorted(residual_pids),
                    residual_directories[:16]))
        time.sleep(0.01)


def remove_empty_campaign_cgroup(
        campaign_scope: Path, parent_cgroup: Path) -> None:
    """Let the last shard supervisor reclaim an orphaned campaign scope.

    All shard scopes are created before any supervisor starts.  Therefore an
    empty campaign scope means every sibling supervisor has already completed
    its own bounded teardown.  Earlier supervisors merely observe ENOTEMPTY;
    the last one removes the validated scope.  This also works after SIGKILL
    closes the coordinator's inherited liveness pipe.
    """
    root = Path("/sys/fs/cgroup").resolve(strict=True)
    require(parent_cgroup.is_absolute() and
            parent_cgroup.resolve(strict=True) == parent_cgroup and
            (parent_cgroup == root or root in parent_cgroup.parents) and
            campaign_scope.is_absolute() and
            campaign_scope.parent == parent_cgroup and
            re.fullmatch(r"leopard2-exhaustive-[0-9]+-[0-9a-f]+",
                         campaign_scope.name) is not None,
            "internal campaign cgroup is invalid")
    try:
        linked = campaign_scope.lstat()
    except FileNotFoundError:
        return
    try:
        resolved = campaign_scope.resolve(strict=True)
    except FileNotFoundError:
        return
    require(stat.S_ISDIR(linked.st_mode) and resolved == campaign_scope,
            "internal campaign cgroup changed")
    try:
        campaign_scope.rmdir()
    except FileNotFoundError:
        return
    except OSError as error:
        if error.errno in (errno.EBUSY, errno.ENOTEMPTY):
            return
        raise


def validate_future_campaign_cgroup(
        campaign_scope: Path, parent_cgroup: Path,
        expected_shard_names: Sequence[str]) -> tuple[str, ...]:
    """Validate the exact future cgroup namespace before its first mkdir."""
    root = Path("/sys/fs/cgroup").resolve(strict=True)
    require(parent_cgroup.is_absolute() and
            parent_cgroup.resolve(strict=True) == parent_cgroup and
            (parent_cgroup == root or root in parent_cgroup.parents) and
            campaign_scope.is_absolute() and
            campaign_scope.parent == parent_cgroup and
            re.fullmatch(r"leopard2-exhaustive-[0-9]+-[0-9a-f]+",
                         campaign_scope.name) is not None,
            "campaign guardian received an invalid future cgroup")
    names = tuple(expected_shard_names)
    require(1 <= len(names) <= 128 and
            len(names) == len(set(names)) and
            list(names) == sorted(names) and
            all(re.fullmatch(r"shard-[0-9]{4}", name) is not None
                for name in names),
            "campaign guardian received invalid shard scope names")
    return names


def guardian_cleanup_campaign_cgroup(
        campaign_scope: Path, parent_cgroup: Path,
        expected_shard_names: Sequence[str]) -> None:
    """Gracefully drain, then kill/remove, the exact accepted inventory."""
    names = validate_future_campaign_cgroup(
        campaign_scope, parent_cgroup, expected_shard_names)
    try:
        linked = campaign_scope.lstat()
    except FileNotFoundError:
        return
    try:
        resolved = campaign_scope.resolve(strict=True)
    except FileNotFoundError:
        return
    require(stat.S_ISDIR(linked.st_mode) and resolved == campaign_scope,
            "campaign guardian scope changed")

    def directory_records() -> set[tuple[str, int, int]]:
        records: set[tuple[str, int, int]] = set()
        for current in cgroup_subtree_directories(campaign_scope):
            try:
                metadata = current.lstat()
                resolved_current = current.resolve(strict=True)
            except FileNotFoundError:
                continue
            require(stat.S_ISDIR(metadata.st_mode) and
                    resolved_current == current and
                    (current == campaign_scope or
                     campaign_scope in current.parents),
                    "campaign guardian subtree identity escaped its root")
            records.add((
                str(current.relative_to(campaign_scope)),
                metadata.st_dev, metadata.st_ino))
        return records

    initial_directories = directory_records()
    allowed = set(names)
    children = []
    unexpected_children = []
    try:
        entries = tuple(campaign_scope.iterdir())
    except FileNotFoundError:
        return
    for child in entries:
        try:
            metadata = child.lstat()
        except FileNotFoundError:
            continue
        if not stat.S_ISDIR(metadata.st_mode):
            continue
        try:
            resolved_child = child.resolve(strict=True)
        except FileNotFoundError:
            continue
        require(child.parent == campaign_scope and resolved_child == child,
                "campaign guardian found a child outside its campaign scope")
        # An unexpected group is evidence corruption, but it is still inside
        # this exact private campaign scope.  Retain the violation and remove
        # the subtree before failing so a dead coordinator cannot leak it.
        if child.name not in allowed:
            unexpected_children.append(child.name)
        children.append(child)

    # A shard adversary can SIGSTOP its supervisor after migrating upward.
    # Resume every exact member through a pidfd so the supervisor can observe
    # parent EOF and perform its stronger subreaper sweep.  Empty, unlaunched
    # scopes are removed immediately.  Only residual scopes are cgroup-killed
    # after the shared bounded grace period.
    def current_pids(child: Path) -> list[int] | None:
        if not child.exists():
            return None
        try:
            return cgroup_subtree_pids(child)
        except FileNotFoundError:
            return None

    def remove_if_empty(child: Path) -> None:
        try:
            remove_empty_cgroup(child)
        except FileNotFoundError:
            pass

    poller = select.poll()
    tracked: dict[int, int] = {}

    def track(pid: int) -> None:
        try:
            descriptor = linux_pidfd_open(pid)
        except ProcessLookupError:
            return
        try:
            poller.register(
                descriptor, select.POLLIN | select.POLLHUP | select.POLLERR)
        except BaseException:
            os.close(descriptor)
            raise
        tracked[descriptor] = pid

    def forget_exited() -> None:
        for descriptor, unused_events in poller.poll(0):
            if descriptor not in tracked:
                continue
            poller.unregister(descriptor)
            os.close(descriptor)
            del tracked[descriptor]

    def signal_tracked(signal_number: int) -> None:
        for descriptor in tuple(tracked):
            try:
                linux_pidfd_send_signal(descriptor, signal_number)
            except ProcessLookupError:
                pass

    try:
        for child in sorted(children):
            pids = current_pids(child)
            if pids is None:
                continue
            if not pids:
                remove_if_empty(child)
                continue
            for pid in pids:
                track(pid)
        deadline = time.monotonic() + CGROUP_TEARDOWN_SECONDS
        while True:
            forget_exited()
            signal_tracked(signal.SIGCONT)
            remaining = []
            for child in sorted(children):
                pids = current_pids(child)
                if pids is None:
                    continue
                if not pids:
                    remove_if_empty(child)
                else:
                    remaining.append(child)
            if not remaining and not tracked:
                break
            if time.monotonic() >= deadline:
                break
            time.sleep(0.01)

        # Graceful supervisors have now had a full bounded cleanup window.
        # Kill only residual cgroup members and still-live exact initial
        # identities, including a supervisor that moved itself upward.
        for child in remaining:
            signal_cgroup(child)
            remove_if_empty(child)
        signal_tracked(signal.SIGKILL)
        kill_deadline = time.monotonic() + CGROUP_TEARDOWN_SECONDS
        while tracked:
            forget_exited()
            if not tracked:
                break
            require(time.monotonic() < kill_deadline,
                    "campaign guardian retained live supervisor identities")
            time.sleep(0.01)

        # A tracked same-UID process can create or enter another cgroup after
        # the initial inventory when SIGCONT resumes it.  Once every exact
        # pidfd identity is dead, atomically kill any other process still
        # inside this private campaign root, inventory the now-stable subtree,
        # and remove it deepest-first.  The recursive walker validates every
        # directory against campaign_scope and cannot escape to parent_cgroup.
        signal_cgroup(campaign_scope)
        subtree_deadline = time.monotonic() + CGROUP_TEARDOWN_SECONDS
        while campaign_scope.exists():
            in_scope_pids = cgroup_subtree_pids(campaign_scope)
            if not in_scope_pids:
                break
            require(time.monotonic() < subtree_deadline,
                    "campaign guardian retained late in-scope processes: %s" %
                    in_scope_pids)
            signal_cgroup(campaign_scope)
            time.sleep(0.01)
        final_directories = directory_records()
        late_directories = final_directories - initial_directories
        remove_empty_cgroup(campaign_scope)
        require(not campaign_scope.exists(),
                "campaign guardian could not remove its campaign scope")
        require(not unexpected_children and not late_directories,
                "campaign guardian removed unexpected cgroups: "
                "initial=%s late=%s" % (
                    sorted(unexpected_children),
                    sorted(record[0] for record in late_directories)))
    finally:
        for descriptor in tuple(tracked):
            try:
                poller.unregister(descriptor)
            except (KeyError, OSError):
                pass
            os.close(descriptor)
            del tracked[descriptor]


def write_guardian_status(descriptor: int, payload: bytes) -> bool:
    """Write one bounded, pipe-atomic guardian status record."""
    require(type(descriptor) is int and descriptor >= 0 and
            isinstance(payload, bytes) and 1 <= len(payload) <= 64 and
            payload.endswith(b"\n") and b"\n" not in payload[:-1] and
            all(byte < 128 for byte in payload),
            "campaign guardian status record is invalid")
    try:
        written = os.write(descriptor, payload)
    except BrokenPipeError:
        return False
    require(written == len(payload),
            "campaign guardian status write was short")
    return True


def guardian_reap_owned_children(
        processes: Mapping[int, subprocess.Popen[bytes]],
        allow_eof_grace: bool) -> None:
    """Kill and reap the guardian's complete, immutable launch lineage."""
    if allow_eof_grace:
        for process in processes.values():
            if process.poll() is None:
                pidfd_signal(process.pid, signal.SIGCONT)
        deadline = time.monotonic() + SUPERVISOR_EOF_GRACE_SECONDS
        while (time.monotonic() < deadline and
               any(process.poll() is None
                   for process in processes.values())):
            time.sleep(0.01)
    for process in processes.values():
        if process.poll() is None:
            pidfd_signal(process.pid, signal.SIGCONT)
            pidfd_kill(process.pid)
    deadline = time.monotonic() + CGROUP_TEARDOWN_SECONDS
    for process in processes.values():
        if process.poll() is not None:
            continue
        try:
            process.wait(timeout=max(0.001, deadline - time.monotonic()))
        except subprocess.TimeoutExpired as error:
            raise EvidenceError(
                "campaign guardian could not reap a shard supervisor") \
                from error
    # Once every direct supervisor has been reaped, any double-forked or
    # upward-migrated shard process is a direct child of this subreaper.
    kill_and_reap_campaign_children()


def supervise_campaign_guardian(
        campaign_scope: Path, parent_cgroup: Path,
        launches: Sequence[tuple[int, int]], shard_count: int,
        output_root: Path, timeout: float,
        runner_descriptor: int, runner_size: int, runner_sha256: str,
        binary_descriptor: int, binary_size: int, binary_sha256: str,
        taskset_descriptor: int, taskset_size: int, taskset_sha256: str,
        interpreter_identity: Mapping[str, Any],
        parent_watch_descriptor: int, ready_descriptor: int,
        status_descriptor: int) -> int:
    """Fork, monitor, and reap every supervisor in one owned lineage."""
    launch_records = tuple(launches)
    require(type(shard_count) is int and 1 <= shard_count <= 128 and
            type(timeout) in (int, float) and math.isfinite(timeout) and
            timeout > 0 and
            1 <= len(launch_records) <= shard_count and
            all(type(index) is int and type(cpu) is int and
                0 <= index < shard_count and cpu >= 0
                for index, cpu in launch_records) and
            [index for index, unused_cpu in launch_records] ==
                sorted({index for index, unused_cpu in launch_records}) and
            len({cpu for unused_index, cpu in launch_records}) ==
                len(launch_records) and
            all(cpu in os.sched_getaffinity(0)
                for unused_index, cpu in launch_records),
            "campaign guardian launch records are invalid")
    names = tuple(
        "shard-%04d" % index for index, unused_cpu in launch_records)
    validate_future_campaign_cgroup(
        campaign_scope, parent_cgroup, names)
    require(not campaign_scope.exists(),
            "campaign cgroup existed before guardian readiness")
    require(cgroup_v2_current() == parent_cgroup,
            "campaign guardian started outside its parent cgroup")
    output_root = Path(os.path.abspath(output_root))
    require(output_root == output_root.resolve(strict=True),
            "campaign guardian output root is not canonical")
    canonical_directory(
        output_root, output_root, "campaign guardian output root")
    for shard_index, unused_cpu in launch_records:
        directory = output_root / ("shard-%04d" % shard_index)
        canonical_directory(
            directory, output_root, "campaign guardian shard directory")
        require(not (directory / "stdout.json").exists() and
                not (directory / "stdout.json").is_symlink() and
                not (directory / "stderr.txt").exists() and
                not (directory / "stderr.txt").is_symlink(),
                "campaign guardian found pre-existing shard output")
    validate_sealed_descriptor(
        runner_descriptor, runner_size, runner_sha256, 0o400,
        "campaign guardian runner")
    validate_sealed_descriptor(
        binary_descriptor, binary_size, binary_sha256, 0o500,
        "campaign guardian binary")
    validate_sealed_descriptor(
        taskset_descriptor, taskset_size, taskset_sha256, 0o500,
        "campaign guardian taskset")
    require(process_executable_identity_matches(interpreter_identity),
            "campaign guardian interpreter inode identity changed")
    for descriptor, label in (
            (parent_watch_descriptor, "parent control"),
            (ready_descriptor, "readiness"),
            (status_descriptor, "status")):
        descriptor_status = os.fstat(descriptor)
        require(stat.S_ISFIFO(descriptor_status.st_mode),
                "campaign guardian %s descriptor is invalid" % label)
    validate_pidfd_support()
    require(not adopted_child_pids(),
            "campaign guardian found pre-existing child processes")
    set_child_subreaper()
    arm_parent_death_continue(os.getppid())

    selector = selectors.DefaultSelector()
    processes: dict[int, subprocess.Popen[bytes]] = {}
    parent_dead = False
    campaign_failed = False
    try:
        os.set_blocking(parent_watch_descriptor, False)
        os.set_blocking(status_descriptor, True)
        selector.register(parent_watch_descriptor, selectors.EVENT_READ)
        # Readiness is the last fallible initialization step.  No cgroup may
        # exist and no supervisor may be forked until the coordinator receives
        # this byte, creates every exact scope, and replies with "L".
        written = os.write(ready_descriptor, b"R")
        require(written == 1, "campaign guardian readiness write was short")
        os.close(ready_descriptor)
        ready_descriptor = -1

        launch_authorized = False
        while not launch_authorized:
            events = selector.select()
            require(events, "campaign guardian control wait was empty")
            try:
                block = os.read(parent_watch_descriptor, 2)
            except BlockingIOError:
                continue
            if not block:
                parent_dead = True
                break
            require(block == b"L",
                    "campaign guardian received invalid launch control")
            launch_authorized = True

        if launch_authorized:
            require(campaign_scope.resolve(strict=True) == campaign_scope and
                    campaign_scope.is_dir(),
                    "campaign cgroup was not created as authorized")
            for name in names:
                scope = campaign_scope / name
                require(scope.resolve(strict=True) == scope and
                        scope.is_dir(),
                        "campaign shard cgroup was not created as authorized")

        for shard_index, cpu in launch_records if launch_authorized else ():
            events = selector.select(0)
            if events:
                block = os.read(parent_watch_descriptor, 1)
                require(not block,
                        "campaign guardian received extra control data")
                parent_dead = True
                break
            directory = output_root / ("shard-%04d" % shard_index)
            scope = campaign_scope / ("shard-%04d" % shard_index)
            command = [
                *supervisor_command_prefix(runner_descriptor),
                "--scope", str(scope),
                "--parent-cgroup", str(parent_cgroup),
                "--coordinator-pid", str(os.getpid()),
                "--binary-fd", str(binary_descriptor),
                "--binary-size", str(binary_size),
                "--binary-sha256", binary_sha256,
                "--taskset-fd", str(taskset_descriptor),
                "--taskset-size", str(taskset_size),
                "--taskset-sha256", taskset_sha256,
                "--cpu", str(cpu),
                "--shard-index", str(shard_index),
                "--shard-count", str(shard_count),
                "--stdout", str(directory / "stdout.json"),
                "--stderr", str(directory / "stderr.txt"),
                "--output-root", str(output_root),
                "--parent-watch-fd", str(parent_watch_descriptor),
                "--timeout", str(timeout),
            ]
            process = launch_current_interpreter(
                command, stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                start_new_session=True,
                pass_fds=(
                    runner_descriptor, binary_descriptor,
                    taskset_descriptor, parent_watch_descriptor),
                env=CHILD_ENV)
            processes[shard_index] = process
            if not write_guardian_status(
                    status_descriptor,
                    ("S %d %d\n" % (shard_index, process.pid)).encode(
                        "ascii")):
                parent_dead = True
                break

        active = set(processes)
        deadline = time.monotonic() + timeout + SUPERVISOR_GRACE_SECONDS
        while active and not parent_dead and not campaign_failed:
            events = selector.select(0.05)
            if events:
                block = os.read(parent_watch_descriptor, 1)
                require(not block,
                        "campaign guardian received extra control data")
                parent_dead = True
                break
            for shard_index in sorted(tuple(active)):
                process = processes[shard_index]
                return_code = process.poll()
                if return_code is None:
                    continue
                active.remove(shard_index)
                if not write_guardian_status(
                        status_descriptor,
                        ("E %d %d\n" % (
                            shard_index, return_code)).encode("ascii")):
                    parent_dead = True
                    break
                if return_code != 0:
                    campaign_failed = True
                    break
            if active and time.monotonic() >= deadline:
                shard_index = min(active)
                if write_guardian_status(
                        status_descriptor,
                        ("E %d 124\n" % shard_index).encode("ascii")):
                    campaign_failed = True
                else:
                    parent_dead = True

        if (launch_authorized and not parent_dead and not campaign_failed and
                len(processes) == len(launch_records)):
            require(not adopted_child_pids(),
                    "successful supervisors left adopted descendants")
            if write_guardian_status(status_descriptor, b"D\n"):
                while True:
                    events = selector.select()
                    require(events,
                            "campaign guardian completion wait was empty")
                    block = os.read(parent_watch_descriptor, 1)
                    if not block:
                        parent_dead = True
                        break
                    raise EvidenceError(
                        "campaign guardian received extra control data")
            else:
                parent_dead = True
    finally:
        selector.close()
        if ready_descriptor >= 0:
            try:
                os.close(ready_descriptor)
            except OSError:
                pass
        try:
            guardian_reap_owned_children(processes, parent_dead)
        finally:
            try:
                guardian_cleanup_campaign_cgroup(
                    campaign_scope, parent_cgroup, names)
            finally:
                try:
                    os.close(status_descriptor)
                except OSError:
                    pass
    return 0


def campaign_guardian_main(arguments: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--campaign-scope", type=Path, required=True)
    parser.add_argument("--parent-cgroup", type=Path, required=True)
    parser.add_argument("--launch", action="append", required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--timeout", type=float, required=True)
    parser.add_argument("--runner-fd", type=int, required=True)
    parser.add_argument("--runner-size", type=int, required=True)
    parser.add_argument("--runner-sha256", required=True)
    parser.add_argument("--binary-fd", type=int, required=True)
    parser.add_argument("--binary-size", type=int, required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--taskset-fd", type=int, required=True)
    parser.add_argument("--taskset-size", type=int, required=True)
    parser.add_argument("--taskset-sha256", required=True)
    parser.add_argument("--interpreter-path", required=True)
    parser.add_argument("--interpreter-device", type=int, required=True)
    parser.add_argument("--interpreter-inode", type=int, required=True)
    parser.add_argument("--interpreter-mode", type=int, required=True)
    parser.add_argument("--interpreter-size", type=int, required=True)
    parser.add_argument("--interpreter-sha256", required=True)
    parser.add_argument("--parent-watch-fd", type=int, required=True)
    parser.add_argument("--ready-fd", type=int, required=True)
    parser.add_argument("--status-fd", type=int, required=True)
    options = parser.parse_args(arguments)
    launches = []
    for value in options.launch:
        match = re.fullmatch(r"([0-9]+):([0-9]+)", value)
        require(match is not None,
                "campaign guardian launch argument is invalid")
        launches.append((int(match.group(1)), int(match.group(2))))
    interpreter = {
        "schema": "leopard2-process-executable/v1",
        "path": options.interpreter_path,
        "device": options.interpreter_device,
        "inode": options.interpreter_inode,
        "mode": options.interpreter_mode,
        "size": options.interpreter_size,
        "sha256": options.interpreter_sha256,
    }
    return supervise_campaign_guardian(
        options.campaign_scope, options.parent_cgroup,
        launches, options.shard_count, options.output_root, options.timeout,
        options.runner_fd, options.runner_size, options.runner_sha256,
        options.binary_fd, options.binary_size, options.binary_sha256,
        options.taskset_fd, options.taskset_size, options.taskset_sha256,
        interpreter, options.parent_watch_fd, options.ready_fd,
        options.status_fd)


def bounded_append(target: bytearray, block: bytes, limit: int) -> bool:
    remaining = max(0, limit - len(target))
    target.extend(block[:remaining])
    return len(block) > remaining


def supervise_shard(
        scope: Path, parent_cgroup: Path, coordinator_pid: int,
        binary_descriptor: int,
        binary_size: int, binary_sha256: str, taskset_descriptor: int,
        taskset_size: int, taskset_sha256: str, cpu: int,
        shard_index: int, shard_count: int,
        stdout_path: Path, stderr_path: Path, output_root: Path,
        parent_watch_descriptor: int, timeout: float) -> int:
    """Internal one-shard supervisor executed in a dedicated cgroup."""
    require(type(cpu) is int and cpu >= 0 and
            type(coordinator_pid) is int and coordinator_pid > 1 and
            type(shard_index) is int and type(shard_count) is int and
            0 <= shard_index < shard_count <= 128 and
            type(timeout) in (int, float) and math.isfinite(timeout) and
            timeout > 0,
            "internal shard execution arguments are invalid")
    arm_parent_death_continue(coordinator_pid)
    require(cgroup_v2_current() == parent_cgroup,
            "internal supervisor started outside its delegated cgroup")
    cgroup_root = Path("/sys/fs/cgroup").resolve(strict=True)
    require(parent_cgroup.is_absolute() and
            parent_cgroup.resolve(strict=True) == parent_cgroup and
            (parent_cgroup == cgroup_root or
             cgroup_root in parent_cgroup.parents) and
            scope.is_absolute() and scope.resolve(strict=True) == scope and
            scope.parent.parent == parent_cgroup and
            re.fullmatch(r"leopard2-exhaustive-[0-9]+-[0-9a-f]+",
                         scope.parent.name) is not None and
            scope.name == "shard-%04d" % shard_index and scope.is_dir(),
            "internal shard cgroup is invalid")
    output_root = Path(os.path.abspath(output_root))
    require(output_root == output_root.resolve(strict=True),
            "internal campaign root is not canonical")
    canonical_directory(
        output_root, output_root, "internal campaign root")
    directory = output_root / ("shard-%04d" % shard_index)
    canonical_directory(
        directory, output_root, "internal shard output directory")
    require(stdout_path == directory / "stdout.json" and
            stderr_path == directory / "stderr.txt" and
            not stdout_path.exists() and not stdout_path.is_symlink() and
            not stderr_path.exists() and not stderr_path.is_symlink(),
            "internal shard output paths are invalid")
    require(cpu in os.sched_getaffinity(0),
            "internal shard CPU is outside the allowed affinity")
    validate_sealed_descriptor(
        binary_descriptor, binary_size, binary_sha256, 0o500,
        "exhaustive binary")
    validate_sealed_descriptor(
        taskset_descriptor, taskset_size, taskset_sha256, 0o500,
        "taskset helper")
    watch_status = os.fstat(parent_watch_descriptor)
    require(stat.S_ISFIFO(watch_status.st_mode),
            "internal parent-liveness descriptor is invalid")
    validate_pidfd_support()
    set_child_subreaper()
    cgroup_write(scope / "cgroup.procs", str(os.getpid()).encode("ascii"))
    command = [
        "/proc/self/fd/%d" % taskset_descriptor, "-c", str(cpu),
        "/proc/self/fd/%d" % binary_descriptor,
        "--shard-index", str(shard_index),
        "--shard-count", str(shard_count),
    ]
    selector = selectors.DefaultSelector()
    process: subprocess.Popen[bytes] | None = None
    stdout = bytearray()
    stderr = bytearray()
    failure: str | None = None
    return_code = -int(signal.SIGKILL)
    try:
        process = subprocess.Popen(
            command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, start_new_session=True,
            pass_fds=(taskset_descriptor, binary_descriptor), env=CHILD_ENV)
        require(process.stdout is not None and process.stderr is not None,
                "internal supervisor cannot capture shard output")
        streams = (
            (process.stdout, stdout, MAX_SHARD_JSON_BYTES, "stdout"),
            (process.stderr, stderr, MAX_SHARD_JSON_BYTES, "stderr"),
        )
        for stream, target, limit, name in streams:
            os.set_blocking(stream.fileno(), False)
            selector.register(
                stream, selectors.EVENT_READ,
                ("output", target, limit, name))
        os.set_blocking(parent_watch_descriptor, False)
        selector.register(
            parent_watch_descriptor, selectors.EVENT_READ,
            ("parent", None, 0, "parent"))
        open_outputs = len(streams)
        deadline = time.monotonic() + timeout
        while True:
            if open_outputs == 0:
                polled = process.poll()
                if polled is not None:
                    return_code = polled
                    break
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                failure = "shard exceeded its execution timeout"
                break
            for key, unused_mask in selector.select(min(remaining, 0.05)):
                kind, target, limit, name = key.data
                stream = key.fileobj
                descriptor = (
                    stream if isinstance(stream, int) else stream.fileno())
                try:
                    block = os.read(descriptor, 65536)
                except BlockingIOError:
                    continue
                if kind == "parent":
                    failure = (
                        "campaign parent exited before shard completion"
                        if not block else
                        "campaign parent-liveness pipe received data")
                    break
                if not block:
                    selector.unregister(stream)
                    open_outputs -= 1
                    continue
                require(isinstance(target, bytearray),
                        "internal output capture state is invalid")
                if bounded_append(target, block, limit):
                    failure = "shard %s exceeded its byte limit" % name
                    break
            if failure is not None:
                break
    finally:
        selector.close()
        if process is not None:
            if process.stdout is not None:
                process.stdout.close()
            if process.stderr is not None:
                process.stderr.close()
        # Move the supervisor out before atomically killing every inherited
        # member, including descendants that called setsid() and double-forked.
        cgroup_write(
            parent_cgroup / "cgroup.procs",
            str(os.getpid()).encode("ascii"))
        cleanup_error: BaseException | None = None
        try:
            # Cgroup kill is an acceleration for the ordinary case.  It is
            # not the trust boundary: a same-UID descendant may have migrated
            # upward before cleanup.
            signal_cgroup(scope)
        except BaseException as error:
            cleanup_error = error
        if process is not None:
            try:
                kill_and_reap_descendants(process)
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
            if isinstance(process.returncode, int):
                return_code = process.returncode
        try:
            remove_empty_cgroup(scope)
        except BaseException as error:
            if cleanup_error is None:
                cleanup_error = error
        try:
            remove_empty_campaign_cgroup(scope.parent, parent_cgroup)
        except BaseException as error:
            if cleanup_error is None:
                cleanup_error = error
        if cleanup_error is not None:
            raise cleanup_error
    if failure is not None:
        diagnostic = ("\n" + failure + "\n").encode("utf-8")
        bounded_append(stderr, diagnostic, MAX_SHARD_JSON_BYTES)
    atomic_bytes(stdout_path, bytes(stdout), MAX_SHARD_JSON_BYTES)
    atomic_bytes(stderr_path, bytes(stderr), MAX_SHARD_JSON_BYTES)
    if failure is not None:
        return 124
    return 0 if return_code == 0 else 125


def supervisor_main(arguments: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--scope", type=Path, required=True)
    parser.add_argument("--parent-cgroup", type=Path, required=True)
    parser.add_argument("--coordinator-pid", type=int, required=True)
    parser.add_argument("--binary-fd", type=int, required=True)
    parser.add_argument("--binary-size", type=int, required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--taskset-fd", type=int, required=True)
    parser.add_argument("--taskset-size", type=int, required=True)
    parser.add_argument("--taskset-sha256", required=True)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--stdout", type=Path, required=True)
    parser.add_argument("--stderr", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--parent-watch-fd", type=int, required=True)
    parser.add_argument("--timeout", type=float, required=True)
    options = parser.parse_args(arguments)
    return supervise_shard(
        options.scope, options.parent_cgroup, options.coordinator_pid,
        options.binary_fd,
        options.binary_size, options.binary_sha256, options.taskset_fd,
        options.taskset_size, options.taskset_sha256, options.cpu,
        options.shard_index, options.shard_count, options.stdout,
        options.stderr, options.output_root, options.parent_watch_fd,
        options.timeout)


def terminate_supervisor(
        process: subprocess.Popen[bytes], scope: Path,
        graceful_deadline: float | None = None) -> None:
    """Reap one supervisor, allowing its exact pidfd sweep before force.

    The caller closes the shared parent-liveness writer before supplying a
    graceful deadline.  That EOF asks every supervisor to kill and reap its
    own complete descendant tree.  Cgroup/SIGKILL is only the bounded fallback
    when the supervisor cannot complete that protocol.
    """
    if process.poll() is None and graceful_deadline is not None:
        remaining = graceful_deadline - time.monotonic()
        if remaining > 0:
            try:
                process.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                pass
    if process.poll() is None:
        signal_cgroup(scope)
        pidfd_kill(process.pid)
    try:
        process.wait(timeout=CGROUP_TEARDOWN_SECONDS)
    except subprocess.TimeoutExpired as error:
        raise EvidenceError("shard supervisor could not be reaped") from error
    if scope.exists():
        signal_cgroup(scope)
        remove_empty_cgroup(scope)


def supervisor_command_prefix(runner_descriptor: int) -> list[str]:
    require(type(runner_descriptor) is int and runner_descriptor >= 0,
            "supervisor runner descriptor is invalid")
    return [
        interpreter_discovery_argv0(),
        *SUPERVISOR_PYTHON_FLAGS,
        "/proc/self/fd/%d" % runner_descriptor,
        "--internal-shard-supervisor",
    ]


def guardian_command(
        runner_descriptor: int, runner_snapshot: Mapping[str, Any],
        binary_descriptor: int, binary_snapshot: Mapping[str, Any],
        taskset_descriptor: int, taskset_snapshot: Mapping[str, Any],
        interpreter_identity: Mapping[str, Any],
        campaign_scope: Path, parent_cgroup: Path,
        launches: Sequence[tuple[int, int]], shard_count: int,
        output_root: Path, timeout: float,
        parent_watch_descriptor: int, ready_descriptor: int,
        status_descriptor: int) -> list[str]:
    launch_records = tuple(launches)
    names = tuple(
        "shard-%04d" % index for index, unused_cpu in launch_records)
    validate_future_campaign_cgroup(
        campaign_scope, parent_cgroup, names)
    for snapshot, label in (
            (runner_snapshot, "runner"),
            (binary_snapshot, "binary"),
            (taskset_snapshot, "taskset")):
        require(isinstance(snapshot, Mapping) and
                type(snapshot.get("size")) is int and
                isinstance(snapshot.get("sha256"), str),
                "campaign guardian %s snapshot is invalid" % label)
    require(type(shard_count) is int and 1 <= shard_count <= 128 and
            launch_records and
            all(type(index) is int and type(cpu) is int
                for index, cpu in launch_records) and
            type(timeout) in (int, float) and math.isfinite(timeout) and
            timeout > 0,
            "campaign guardian command arguments are invalid")
    command = [
        interpreter_discovery_argv0(),
        *SUPERVISOR_PYTHON_FLAGS,
        "/proc/self/fd/%d" % runner_descriptor,
        "--internal-campaign-guardian",
        "--campaign-scope", str(campaign_scope),
        "--parent-cgroup", str(parent_cgroup),
        "--shard-count", str(shard_count),
        "--output-root", str(output_root),
        "--timeout", str(timeout),
    ]
    for shard_index, cpu in launch_records:
        command.extend(("--launch", "%d:%d" % (shard_index, cpu)))
    command.extend((
        "--runner-fd", str(runner_descriptor),
        "--runner-size", str(runner_snapshot["size"]),
        "--runner-sha256", runner_snapshot["sha256"],
        "--binary-fd", str(binary_descriptor),
        "--binary-size", str(binary_snapshot["size"]),
        "--binary-sha256", binary_snapshot["sha256"],
        "--taskset-fd", str(taskset_descriptor),
        "--taskset-size", str(taskset_snapshot["size"]),
        "--taskset-sha256", taskset_snapshot["sha256"],
        "--interpreter-path", interpreter_identity["path"],
        "--interpreter-device", str(interpreter_identity["device"]),
        "--interpreter-inode", str(interpreter_identity["inode"]),
        "--interpreter-mode", str(interpreter_identity["mode"]),
        "--interpreter-size", str(interpreter_identity["size"]),
        "--interpreter-sha256", interpreter_identity["sha256"],
        "--parent-watch-fd", str(parent_watch_descriptor),
        "--ready-fd", str(ready_descriptor),
        "--status-fd", str(status_descriptor),
    ))
    return command


def wait_guardian_ready(
        process: subprocess.Popen[bytes], descriptor: int,
        timeout: float = CGROUP_TEARDOWN_SECONDS) -> None:
    """Require the validated guardian handshake before the first cgroup mkdir."""
    require(type(descriptor) is int and descriptor >= 0 and
            type(timeout) in (int, float) and math.isfinite(timeout) and
            timeout > 0,
            "campaign guardian readiness arguments are invalid")
    os.set_blocking(descriptor, False)
    selector = selectors.DefaultSelector()
    payload = bytearray()
    deadline = time.monotonic() + timeout
    try:
        selector.register(descriptor, selectors.EVENT_READ)
        while True:
            remaining = deadline - time.monotonic()
            require(remaining > 0,
                    "campaign guardian readiness timed out")
            events = selector.select(min(remaining, 0.05))
            if not events:
                require(process.poll() is None,
                        "campaign guardian exited before readiness")
                continue
            block = os.read(descriptor, 2)
            if block:
                payload.extend(block)
                require(len(payload) <= 1,
                        "campaign guardian readiness payload is invalid")
                continue
            break
    finally:
        selector.close()
    require(bytes(payload) == b"R" and process.poll() is None,
            "campaign guardian did not establish lifecycle ownership")


def start_campaign_guardian(
        descriptor: int, campaign_scope: Path) -> None:
    """Authorize launches only after every exact cgroup has been created."""
    require(type(descriptor) is int and descriptor >= 0 and
            campaign_scope.is_absolute() and campaign_scope.is_dir() and
            campaign_scope.resolve(strict=True) == campaign_scope,
            "campaign guardian launch authorization is invalid")
    written = os.write(descriptor, b"L")
    require(written == 1,
            "campaign guardian launch authorization was short")


def wait_guardian_completion(
        process: subprocess.Popen[bytes], descriptor: int,
        expected_indices: Sequence[int], timeout: float) -> dict[int, int]:
    """Validate the fork server's bounded status protocol through completion."""
    expected = tuple(expected_indices)
    require(type(descriptor) is int and descriptor >= 0 and
            expected and expected == tuple(sorted(set(expected))) and
            type(timeout) in (int, float) and math.isfinite(timeout) and
            timeout > 0,
            "campaign guardian completion arguments are invalid")
    os.set_blocking(descriptor, False)
    selector = selectors.DefaultSelector()
    payload = bytearray()
    launched: dict[int, int] = {}
    completed: dict[int, int] = {}
    next_launch = 0
    deadline = time.monotonic() + timeout

    def accept_record(record: bytes) -> bool:
        nonlocal next_launch
        if record == b"D":
            require(next_launch == len(expected) and
                    set(completed) == set(expected) and
                    all(return_code == 0
                        for return_code in completed.values()),
                    "campaign guardian completed an incomplete shard set")
            return True
        start = re.fullmatch(rb"S ([0-9]+) ([0-9]+)", record)
        if start is not None:
            shard_index = int(start.group(1))
            pid = int(start.group(2))
            require(next_launch < len(expected) and
                    shard_index == expected[next_launch] and pid > 1 and
                    shard_index not in launched and
                    pid not in launched.values(),
                    "campaign guardian launch status is invalid")
            launched[shard_index] = pid
            next_launch += 1
            return False
        end = re.fullmatch(rb"E ([0-9]+) (-?[0-9]+)", record)
        require(end is not None,
                "campaign guardian status record is malformed")
        shard_index = int(end.group(1))
        return_code = int(end.group(2))
        require(shard_index in launched and
                shard_index not in completed and
                -(1 << 31) <= return_code < (1 << 31),
                "campaign guardian exit status is invalid")
        completed[shard_index] = return_code
        require(return_code == 0,
                "shard %d failed with %d; see its stderr.txt" % (
                    shard_index, return_code))
        return False

    try:
        selector.register(descriptor, selectors.EVENT_READ)
        while True:
            remaining = deadline - time.monotonic()
            require(remaining > 0,
                    "campaign guardian completion timed out")
            events = selector.select(min(remaining, 0.05))
            if not events:
                require(process.poll() is None,
                        "campaign guardian exited before completion")
                continue
            block = os.read(descriptor, 4096)
            if not block:
                require(process.poll() is None,
                        "campaign guardian exited before completion")
                raise EvidenceError(
                    "campaign guardian status closed before completion")
            require(len(payload) + len(block) <=
                        MAX_GUARDIAN_STATUS_BYTES,
                    "campaign guardian status exceeded its byte bound")
            payload.extend(block)
            while b"\n" in payload:
                record, unused_separator, remainder = payload.partition(b"\n")
                payload[:] = remainder
                if accept_record(bytes(record)):
                    require(not payload and process.poll() is None,
                            "campaign guardian completion status is invalid")
                    return launched
    finally:
        selector.close()


def launch_shards(
        output: Path, request: Mapping[str, Any],
        expected_shards: Sequence[Mapping[str, Any]],
        lock_descriptor: int | None = None) -> list[dict[str, Any]]:
    workers = request["workers"]
    shard_count = request["shard_count"]
    require(workers == shard_count == len(expected_shards) and workers > 0,
            "one deterministic shard per worker is required")
    binary = Path(request["frozen_binary"]["path"])
    exact_file_identity(binary, request["frozen_binary"],
                        "frozen exhaustive executable")
    exact_file_identity(
        Path(request["taskset"]["path"]), request["taskset"],
        "canonical taskset helper")
    results: list[dict[str, Any] | None] = [None] * shard_count
    pending = []
    for shard_index in range(shard_count):
        retained = retained_shard(
            output, shard_index, request, expected_shards[shard_index])
        if retained is None:
            pending.append(shard_index)
        else:
            results[shard_index] = retained
    if not pending:
        return [value for value in results if value is not None]

    runner_path = Path(__file__).resolve(strict=True)
    exact_file_identity(
        runner_path, request["runner"], "exhaustive campaign runner")
    if "interpreter" in request:
        interpreter_identity = request["interpreter"]
        require(process_executable_identity_matches(interpreter_identity),
                "campaign interpreter inode identity changed")
    else:
        interpreter_identity = process_executable_identity()
    taskset_path = Path(request["taskset"]["path"])
    binary_descriptor = -1
    runner_descriptor = -1
    taskset_descriptor = -1
    parent_watch_read = -1
    parent_watch_write = -1
    guardian_ready_read = -1
    guardian_ready_write = -1
    guardian_status_read = -1
    guardian_status_write = -1
    guardian_process: subprocess.Popen[bytes] | None = None
    campaign_scope: Path | None = None
    created_scopes: set[Path] = set()
    previous_subreaper: int | None = None
    parent_containment_active = False
    try:
        binary_descriptor, binary_snapshot = sealed_file_snapshot(
            binary, request["frozen_binary"], "exhaustive-binary", 0o500)
        runner_descriptor, runner_snapshot = sealed_file_snapshot(
            runner_path, request["runner"], "exhaustive-runner", 0o400)
        taskset_descriptor, taskset_snapshot = sealed_file_snapshot(
            taskset_path, request["taskset"], "taskset-helper", 0o500)
        parent_watch_read, parent_watch_write = os.pipe2(
            os.O_CLOEXEC | os.O_NONBLOCK)
        guardian_ready_read, guardian_ready_write = os.pipe2(
            os.O_CLOEXEC | os.O_NONBLOCK)
        guardian_status_read, guardian_status_write = os.pipe2(
            os.O_CLOEXEC | os.O_NONBLOCK)
        require(request.get("execution_policy") != EXECUTION_POLICY or
                lock_descriptor is not None,
                "production shard launch requires the held campaign lock")
        if lock_descriptor is not None:
            require(type(lock_descriptor) is int and lock_descriptor >= 0,
                    "campaign lock descriptor is invalid")
            lock_status = os.fstat(lock_descriptor)
            retained_lock = request.get("lock", {})
            require(stat.S_ISREG(lock_status.st_mode) and
                    lock_status.st_dev == retained_lock.get("device") and
                    lock_status.st_ino == retained_lock.get("inode") and
                    lock_status.st_uid == retained_lock.get("uid") and
                    lock_status.st_nlink == retained_lock.get("nlink") and
                    stat.S_IMODE(lock_status.st_mode) ==
                        retained_lock.get("mode"),
                    "campaign lock descriptor is invalid")
        require(binary_snapshot["sha256"] ==
                    request["frozen_binary"]["sha256"] and
                binary_snapshot["size"] ==
                    request["frozen_binary"]["size"] and
                runner_snapshot["sha256"] ==
                    request["runner"]["sha256"] and
                runner_snapshot["size"] == request["runner"]["size"] and
                taskset_snapshot["sha256"] ==
                    request["taskset"]["sha256"] and
                taskset_snapshot["size"] == request["taskset"]["size"],
                "sealed execution snapshot identity changed")
        # Prepare all application-visible paths before the guardian handshake.
        # Once readiness transfers lifecycle ownership, only exact cgroup
        # creation and the one-byte launch authorization remain.
        for shard_index in pending:
            directory = ensure_directory(
                output / ("shard-%04d" % shard_index), output)
            cleanup_stale_temporaries(
                directory, output,
                (r"(?:stdout\.json|stderr\.txt|result\.json)"
                 r"\.tmp-[0-9]+",))
            for stale in (
                    directory / "stdout.json",
                    directory / "stderr.txt"):
                if stale.exists() or stale.is_symlink():
                    canonical_artifact(stale, output, "stale shard output")
                    stale.unlink()
            require(not (directory / "result.json").is_symlink(),
                    "pending shard result is a symlink")
        parent_cgroup = cgroup_v2_current()
        validate_pidfd_support()
        previous_subreaper = child_subreaper_state()
        require(not adopted_child_pids(),
                "campaign containment found pre-existing child processes")
        try:
            set_child_subreaper(1)
        except BaseException:
            if child_subreaper_state() != previous_subreaper:
                set_child_subreaper(previous_subreaper)
            raise
        parent_containment_active = True
        campaign_scope = parent_cgroup / (
            "leopard2-exhaustive-%d-%x" % (os.getpid(), time.time_ns()))
        launches = tuple(
            (shard_index, request["allowed_cpus"][shard_index])
            for shard_index in pending)
        expected_scope_names = tuple(
            "shard-%04d" % shard_index
            for shard_index, unused_cpu in launches)
        validate_future_campaign_cgroup(
            campaign_scope, parent_cgroup, expected_scope_names)
        guardian_launch = guardian_command(
            runner_descriptor, runner_snapshot,
            binary_descriptor, binary_snapshot,
            taskset_descriptor, taskset_snapshot,
            interpreter_identity, campaign_scope, parent_cgroup,
            launches, shard_count, output,
            request["timeout_seconds_per_shard"],
            parent_watch_read, guardian_ready_write,
            guardian_status_write)
        guardian_descriptors = [
            runner_descriptor, binary_descriptor, taskset_descriptor,
            parent_watch_read, guardian_ready_write, guardian_status_write,
        ]
        if lock_descriptor is not None:
            guardian_descriptors.append(lock_descriptor)
        guardian_process = launch_current_interpreter(
            guardian_launch, stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            start_new_session=True,
            pass_fds=tuple(guardian_descriptors), env=CHILD_ENV)
        os.close(guardian_ready_write)
        guardian_ready_write = -1
        os.close(guardian_status_write)
        guardian_status_write = -1
        os.close(parent_watch_read)
        parent_watch_read = -1
        wait_guardian_ready(
            guardian_process, guardian_ready_read)
        os.close(guardian_ready_read)
        guardian_ready_read = -1
        require(guardian_process.poll() is None,
                "campaign guardian exited before cgroup creation")
        os.mkdir(campaign_scope)
        for shard_index in pending:
            require(guardian_process.poll() is None,
                    "campaign guardian exited during cgroup creation")
            scope = campaign_scope / ("shard-%04d" % shard_index)
            os.mkdir(scope)
            created_scopes.add(scope)
        start_campaign_guardian(parent_watch_write, campaign_scope)
        launched = wait_guardian_completion(
            guardian_process, guardian_status_read, pending,
            request["timeout_seconds_per_shard"] +
            SUPERVISOR_GRACE_SECONDS + CGROUP_TEARDOWN_SECONDS)
        require(set(launched) == set(pending),
                "campaign guardian omitted a shard supervisor")
        completed_count = shard_count - len(pending)
        for shard_index in pending:
            directory = output / ("shard-%04d" % shard_index)
            stderr_path = directory / "stderr.txt"
            parsed = load_json(
                directory / "stdout.json", MAX_SHARD_JSON_BYTES)
            value = validate_shard(
                parsed, shard_index, shard_count,
                expected_shards[shard_index])
            require(value["mode"] == request["build_provenance"]["mode"],
                    "shard mode differs from build provenance")
            cpu = request["allowed_cpus"][shard_index]
            command = [
                request["taskset"]["path"], "-c", str(cpu), str(binary),
                "--shard-index", str(shard_index),
                "--shard-count", str(shard_count),
            ]
            envelope = {
                "shard": value,
                "stdout": file_identity(directory / "stdout.json"),
                "stderr": file_identity(stderr_path),
                "binary": request["frozen_binary"],
                "cpu": cpu,
                "command": command,
                "execution": {
                    "policy": request.get(
                        "execution_policy", EXECUTION_POLICY),
                    "interpreter": (
                        request["interpreter"]
                        if "interpreter" in request else
                        interpreter_identity),
                    "binary_snapshot": binary_snapshot,
                    "runner_snapshot": runner_snapshot,
                    "taskset_snapshot": taskset_snapshot,
                },
            }
            atomic_json(directory / "result.json", envelope)
            results[shard_index] = value
            completed_count += 1
            print("[%d/%d] shard %d verified %d matrices" % (
                completed_count, shard_count, shard_index,
                value["assigned_patterns"]), flush=True)
    finally:
        cleanup_error: BaseException | None = None
        # EOF transfers cancellation to the guardian, which is the actual
        # parent and subreaper for every supervisor and escaped descendant.
        if parent_watch_write >= 0:
            try:
                os.close(parent_watch_write)
                parent_watch_write = -1
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if guardian_process is not None:
            try:
                try:
                    guardian_return = guardian_process.wait(
                        timeout=SUPERVISOR_GRACE_SECONDS)
                except subprocess.TimeoutExpired:
                    pidfd_kill(guardian_process.pid)
                    guardian_return = guardian_process.wait(
                        timeout=CGROUP_TEARDOWN_SECONDS)
                    raise EvidenceError(
                        "campaign guardian exceeded its cleanup deadline")
                require(guardian_return == 0,
                        "campaign guardian failed with %d" %
                        guardian_return)
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
            guardian_process = None
        # If the guardian itself failed, its entire lineage reparents to this
        # second-level subreaper before these cgroup and pidfd fallbacks.
        for scope in tuple(created_scopes):
            try:
                signal_cgroup(scope)
                remove_empty_cgroup(scope)
                created_scopes.discard(scope)
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if parent_containment_active:
            try:
                kill_and_reap_campaign_children()
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if campaign_scope is not None and campaign_scope.exists():
            try:
                require(not tuple(
                            child for child in campaign_scope.iterdir()
                            if child.is_dir()),
                        "campaign cgroup contains an unknown child scope")
                campaign_scope.rmdir()
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        if parent_containment_active:
            try:
                require(previous_subreaper in (0, 1),
                        "previous campaign subreaper state was lost")
                set_child_subreaper(previous_subreaper)
                parent_containment_active = False
            except BaseException as error:
                if cleanup_error is None:
                    cleanup_error = error
        for descriptor in (
                guardian_ready_write, guardian_ready_read,
                guardian_status_write, guardian_status_read,
                parent_watch_write, parent_watch_read,
                taskset_descriptor, runner_descriptor, binary_descriptor):
            if descriptor >= 0:
                os.close(descriptor)
        if cleanup_error is not None:
            raise cleanup_error
    require(all(value is not None for value in results),
            "exhaustive shard set is incomplete")
    exact_file_identity(binary, request["frozen_binary"],
                        "frozen exhaustive executable")
    exact_file_identity(
        Path(request["taskset"]["path"]), request["taskset"],
        "canonical taskset helper")
    exact_file_identity(
        runner_path, request["runner"], "exhaustive campaign runner")
    require(process_executable_identity_matches(interpreter_identity),
            "campaign interpreter inode identity changed")
    return [value for value in results if value is not None]


def validate_request(value: Any, campaign_root: Path) -> list[dict[str, Any]]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "source_binary", "frozen_binary", "build_provenance",
        "provenance_bundle", "runner", "taskset", "interpreter",
        "workers", "shard_count",
        "expected_patterns", "timeout_seconds_per_shard", "allowed_cpus",
        "assignment", "basis_seed", "expected_shards", "child_environment",
        "execution_policy", "lock", "resume_policy", "digest",
    } and value.get("schema") == SCHEMA and
            value.get("digest") == record_digest(value),
            "exhaustive request is malformed")
    workers = value["workers"]
    require(type(workers) is int and 1 <= workers <= 128 and
            value["shard_count"] == workers and
            value["expected_patterns"] == EXPECTED_PATTERNS and
            value["assignment"] == "global_ordinal_mod_shard_count" and
            value["basis_seed"] == 0 and
            value["child_environment"] == CHILD_ENV and
            value["execution_policy"] == EXECUTION_POLICY and
            value["resume_policy"] ==
                "validated deterministic shard resume" and
            type(value["timeout_seconds_per_shard"]) in (int, float) and
            math.isfinite(value["timeout_seconds_per_shard"]) and
            value["timeout_seconds_per_shard"] > 0,
            "exhaustive request policy is invalid")
    cpus = value["allowed_cpus"]
    require(isinstance(cpus, list) and len(cpus) == workers and
            cpus == sorted(set(cpus)) and
            all(type(cpu) is int and cpu >= 0 for cpu in cpus),
            "exhaustive request CPU map is invalid")
    validate_build_record(value["build_provenance"])
    require(value["build_provenance"]["mode"] in (1, 2),
            "exhaustive request mode is invalid")
    interpreter = value["interpreter"]
    require(isinstance(interpreter, dict) and set(interpreter) == {
                "schema", "path", "device", "inode", "mode", "size",
                "sha256",
            } and
            interpreter["schema"] == "leopard2-process-executable/v1" and
            isinstance(interpreter["path"], str) and
            Path(interpreter["path"]).is_absolute() and
            Path(interpreter["path"]) ==
                Path(os.path.abspath(interpreter["path"])) and
            all(type(interpreter.get(name)) is int and
                interpreter[name] >= 0
                for name in ("device", "inode", "mode", "size")) and
            stat.S_ISREG(interpreter["mode"]) and
            isinstance(interpreter.get("sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}",
                         interpreter["sha256"]) is not None,
            "request interpreter identity is invalid")
    reconstructed = [
        dict(item) for item in reconstruct_expected_shards(workers)
    ]
    require(value["expected_shards"] == reconstructed,
            "exhaustive request descriptor expectations changed")
    frozen_path = campaign_root / "artifacts" / TARGET_NAME
    require(value["frozen_binary"].get("path") == str(frozen_path) and
            value["source_binary"].get("sha256") ==
                value["frozen_binary"].get("sha256") and
            value["source_binary"].get("size") ==
                value["frozen_binary"].get("size"),
            "exhaustive request frozen binary is invalid")
    build_executable = value["build_provenance"]["exhaustive_executable"]
    require(build_executable.get("path") ==
                value["source_binary"].get("path") and
            build_executable.get("size") ==
                value["source_binary"].get("size") and
            build_executable.get("sha256") ==
                value["source_binary"].get("sha256"),
            "request source binary differs from its build provenance")
    source_runner = value["build_provenance"]["source"]["files"][
        "experiments/leopard2/direct_repair/"
        "run_small_direct_exhaustive.py"]
    require(value["runner"] == source_runner and
            value["taskset"].get("path") ==
                str(CANONICAL_TASKSET),
            "request runner/taskset is not bound to canonical inputs")
    lock = value["lock"]
    require(isinstance(lock, dict) and set(lock) == {
        "schema", "path", "device", "inode", "uid", "mode", "nlink",
        "mechanism",
    } and lock["schema"] == "leopard2-global-benchmark-lock/v1" and
            lock["path"] == str(DEFAULT_LOCK) and
            lock["mode"] == 0o600 and lock["nlink"] == 1 and
            lock["mechanism"] == "fcntl-flock-exclusive" and
            all(type(lock[key]) is int and lock[key] >= 0
                for key in ("device", "inode", "uid")),
            "request canonical lock identity is invalid")
    canonical_artifact(frozen_path, campaign_root,
                       "frozen exhaustive executable")
    exact_file_identity(frozen_path, value["frozen_binary"],
                        "frozen exhaustive executable")
    require(stat.S_IMODE(frozen_path.lstat().st_mode) == 0o555,
            "frozen exhaustive executable mode changed")
    validate_bundle(
        value["provenance_bundle"], value["build_provenance"],
        (value["runner"], value["taskset"], value["interpreter"]),
        campaign_root)
    return reconstructed


def verify_manifest(
        manifest_path: Path, no_current_input_check: bool = False) -> int:
    resolved_manifest = manifest_path.resolve(strict=True)
    campaign_root = resolved_manifest.parent
    require(resolved_manifest == campaign_root / "manifest.json",
            "manifest must use its canonical campaign filename")
    canonical_directory(
        campaign_root, campaign_root, "exhaustive campaign root")
    canonical_artifact(
        resolved_manifest, campaign_root, "exhaustive manifest")
    manifest = load_json(resolved_manifest)
    require(isinstance(manifest, dict) and set(manifest) == {
        "schema", "request", "mode", "shard_count", "verified_patterns",
        "recovered_shards", "verified_basis_symbols", "shards", "digest",
    } and manifest.get("schema") == SCHEMA and
            manifest.get("digest") == record_digest(manifest),
            "exhaustive manifest is malformed")
    request = manifest["request"]
    expected = validate_request(request, campaign_root)
    request_path = campaign_root / "request.json"
    canonical_artifact(
        request_path, campaign_root, "exhaustive request")
    require(load_json(request_path) == request,
            "retained exhaustive request changed")
    expected_top = {
        "request.json", "manifest.json", "artifacts",
        *("shard-%04d" % index for index in range(request["shard_count"])),
    }
    require({path.name for path in campaign_root.iterdir()} == expected_top,
            "campaign contains extra or missing top-level artifacts")
    artifacts = campaign_root / "artifacts"
    canonical_directory(artifacts, campaign_root, "artifact directory")
    require({path.name for path in artifacts.iterdir()} ==
            {TARGET_NAME, "provenance"},
            "campaign artifact directory contains extra or missing files")

    shards = []
    for index, expected_shard in enumerate(expected):
        directory = campaign_root / ("shard-%04d" % index)
        canonical_directory(directory, campaign_root, "shard directory")
        require({path.name for path in directory.iterdir()} ==
                    {"stdout.json", "stderr.txt", "result.json"},
                "shard directory contains extra, missing, or linked artifacts")
        value = retained_shard(
            campaign_root, index, request, expected_shard)
        require(value is not None, "retained shard is incomplete")
        shards.append(value)
    mode_set = {value["mode"] for value in shards}
    require(mode_set == {request["build_provenance"]["mode"]} and
            manifest["mode"] == next(iter(mode_set)) and
            manifest["shard_count"] == len(shards) ==
                request["shard_count"] and
            manifest["verified_patterns"] ==
                sum(value["assigned_patterns"] for value in shards) ==
                EXPECTED_PATTERNS and
            manifest["recovered_shards"] ==
                sum(value["recovered_shards"] for value in shards) and
            manifest["verified_basis_symbols"] ==
                sum(value["verified_basis_symbols"] for value in shards) and
            manifest["shards"] == shards,
            "manifest does not reconstruct from exact retained shards")

    if not no_current_input_check:
        exact_file_identity(
            Path(request["source_binary"]["path"]),
            request["source_binary"], "current exhaustive executable")
        current = capture_build_closure(
            Path(request["source_binary"]["path"]))
        require(current == build_without_rebuild_proof(
                    request["build_provenance"]),
                "current exhaustive build provenance changed")
        exact_file_identity(
            Path(request["runner"]["path"]), request["runner"],
            "current exhaustive runner")
        exact_file_identity(
            Path(request["taskset"]["path"]), request["taskset"],
            "current taskset helper")
        require(process_executable_identity_matches(request["interpreter"]),
                "current campaign interpreter inode changed")
        linked_lock = DEFAULT_LOCK.lstat()
        lock = request["lock"]
        require(stat.S_ISREG(linked_lock.st_mode) and
                linked_lock.st_dev == lock["device"] and
                linked_lock.st_ino == lock["inode"] and
                linked_lock.st_uid == lock["uid"] and
                linked_lock.st_nlink == lock["nlink"] and
                stat.S_IMODE(linked_lock.st_mode) == lock["mode"],
                "current canonical lock identity changed")
    print("verified %d exact coefficient matrices in %d shards (mode %d)" % (
        manifest["verified_patterns"], len(shards), manifest["mode"]))
    return 0


def run(options: argparse.Namespace) -> int:
    require(options.binary is not None and options.output is not None,
            "--binary and --output are required for a campaign")
    # Reject adversarial cardinalities before provenance capture, artifact
    # allocation, or the O(EXPECTED_PATTERNS) exact reconstruction.
    selected_cpus = validate_run_limits(
        options.workers, options.timeout, options.lock_timeout)
    source_binary = options.binary.resolve(strict=True)
    require(source_binary.is_file() and os.access(source_binary, os.X_OK),
            "exhaustive verifier is not executable")
    source_root = Path(__file__).resolve(strict=True).parents[3]

    lock = acquire_lock(options.lock, options.lock_timeout)
    try:
        output = prepare_output_root(
            options.output, source_root, options.workers)
        held_lock = lock_identity(lock, options.lock)
        build_record = capture_current_build(source_binary)
        source_binary_identity = file_identity(source_binary)
        require(source_binary_identity["sha256"] ==
                    build_record["exhaustive_executable"]["sha256"] and
                source_binary_identity["size"] ==
                    build_record["exhaustive_executable"]["size"],
                "source binary differs from its build closure")
        artifacts = ensure_directory(output / "artifacts", output)
        frozen_path = artifacts / TARGET_NAME
        frozen_record = freeze_file(
            source_binary_identity, frozen_path, 0o555, output)
        frozen_binary = frozen_record["frozen"]
        runner_identity = file_identity(Path(__file__))
        require(CANONICAL_TASKSET.is_file() and
                os.access(CANONICAL_TASKSET, os.X_OK),
                "canonical /usr/bin/taskset is required")
        taskset_identity = file_identity(CANONICAL_TASKSET)
        interpreter_identity = process_executable_identity()
        bundle = freeze_bundle(
            build_record,
            (runner_identity, taskset_identity, interpreter_identity),
            output)
        expected = [
            dict(value)
            for value in reconstruct_expected_shards(options.workers)
        ]
        request = {
            "schema": SCHEMA,
            "source_binary": source_binary_identity,
            "frozen_binary": frozen_binary,
            "build_provenance": build_record,
            "provenance_bundle": bundle,
            "runner": runner_identity,
            "taskset": taskset_identity,
            "interpreter": interpreter_identity,
            "workers": options.workers,
            "shard_count": options.workers,
            "expected_patterns": EXPECTED_PATTERNS,
            "timeout_seconds_per_shard": options.timeout,
            "allowed_cpus": selected_cpus,
            "assignment": "global_ordinal_mod_shard_count",
            "basis_seed": 0,
            "expected_shards": expected,
            "child_environment": CHILD_ENV,
            "execution_policy": EXECUTION_POLICY,
            "lock": held_lock,
            "resume_policy": "validated deterministic shard resume",
        }
        request["digest"] = record_digest(request)
        request_path = output / "request.json"
        if request_path.exists():
            require(load_json(request_path) == request,
                    "resume request differs from retained campaign")
        else:
            atomic_json(request_path, request)
        manifest_path = output / "manifest.json"
        if manifest_path.exists():
            return verify_manifest(manifest_path, False)
        allowed_resume_names = {
            "request.json", "artifacts",
            *("shard-%04d" % index
              for index in range(options.workers)),
        }
        require({path.name for path in output.iterdir()} <=
                allowed_resume_names,
                "incomplete campaign contains unknown artifacts")
        shards = launch_shards(
            output, request, expected, lock.fileno())
        require(lock_identity(lock, options.lock) == held_lock,
                "canonical lock identity changed during campaign")
        require(capture_build_closure(source_binary) ==
                    build_without_rebuild_proof(build_record),
                "build/source provenance changed during campaign")
        exact_file_identity(source_binary, source_binary_identity,
                            "source exhaustive executable")
        manifest = {
            "schema": SCHEMA,
            "request": request,
            "mode": build_record["mode"],
            "shard_count": len(shards),
            "verified_patterns": sum(
                value["assigned_patterns"] for value in shards),
            "recovered_shards": sum(
                value["recovered_shards"] for value in shards),
            "verified_basis_symbols": sum(
                value["verified_basis_symbols"] for value in shards),
            "shards": shards,
        }
        manifest["digest"] = record_digest(manifest)
        atomic_json(manifest_path, manifest)
        require(lock_identity(lock, options.lock) == held_lock,
                "canonical lock identity changed at manifest commit")
        return verify_manifest(manifest_path, False)
    finally:
        close_campaign_lock(lock)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--binary", type=Path)
    result.add_argument("--output", type=Path)
    result.add_argument(
        "--workers", type=int,
        default=min(128, len(os.sched_getaffinity(0))))
    result.add_argument("--timeout", type=float, default=3600.0)
    result.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    result.add_argument("--lock-timeout", type=float, default=3600.0)
    result.add_argument("--verify", type=Path,
                        help="verify a retained manifest instead of running")
    result.add_argument("--no-current-input-check", action="store_true",
                        help="offline verify using only frozen artifacts")
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    try:
        if arguments is None:
            arguments = sys.argv[1:]
        if (arguments and
                arguments[0] == "--internal-campaign-guardian"):
            return campaign_guardian_main(arguments[1:])
        if (arguments and
                arguments[0] == "--internal-shard-supervisor"):
            return supervisor_main(arguments[1:])
        options = parser().parse_args(arguments)
        if options.verify is not None:
            require(options.binary is None and options.output is None,
                    "--verify cannot be combined with --binary/--output")
            return verify_manifest(
                options.verify, options.no_current_input_check)
        require(not options.no_current_input_check,
                "--no-current-input-check requires --verify")
        return run(options)
    except (AttributeError, EvidenceError, IndexError, KeyError, OSError,
            subprocess.SubprocessError, TypeError, ValueError) as error:
        print("error: %s" % error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
