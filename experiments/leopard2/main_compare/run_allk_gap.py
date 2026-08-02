#!/usr/bin/env python3
"""Broad diagnostic Leopard-main versus Leopard2 all-K gap map.

This is intentionally not promotion evidence.  It saturates all allowed CPUs
with independent single-thread cells to find algorithmic regions worth fixing.
Each cell uses an ABBA process order and retains exact JSON from independently
linked exact-main and Leopard2 executables.
"""

from __future__ import annotations

import argparse
import copy
import concurrent.futures
import ctypes
import dataclasses
import errno
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import secrets
import selectors
import signal
import stat
import statistics
import subprocess
import sys
import threading
import time
from typing import Any, Callable, Mapping, Sequence


_MAIN_COMPARE_DIR = Path(__file__).resolve().parent
if str(_MAIN_COMPARE_DIR) not in sys.path:
    sys.path.insert(0, str(_MAIN_COMPARE_DIR))
import git_capture
_GIT_CAPTURE_PATH = (_MAIN_COMPARE_DIR / "git_capture.py").resolve(strict=True)
if Path(git_capture.__file__).resolve(strict=True) != _GIT_CAPTURE_PATH:
    raise RuntimeError(
        "all-K Git capture helper resolved outside this source tree")

TOOLS_DIRECTORY = Path(__file__).resolve().parents[3] / "tools"
if str(TOOLS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIRECTORY))
from leopard2_build_provenance import (
    BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V2,
    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3,
    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4,
    CANONICAL_REPLAY_RECIPE_SCHEMA,
    CORE_LIBRARY_SOURCES,
    LEGACY_REPLAY_RECIPE_SCHEMA,
    PRODUCTION_BUILD_CLOSURE_SCHEMA,
    PRODUCTION_BUILD_CLOSURE_SCHEMA_V1,
    REPLAY_INVOCATION_SCHEMA,
    REPLAY_INVOCATION_SCHEMA_V1,
    REPRODUCIBLE_BUILD_PROOF_SCHEMA,
    REPRODUCIBLE_BUILD_PROOF_SCHEMA_V2,
    candidate_build_provenance,
    compare_reproducible_builds,
    validate_reproducible_build_proof,
    verify_reproducible_candidate_build,
)


MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
PURE_AVX2_MAIN_SHA256 = (
    "a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910")
ORDER = ("main", "leopard2", "leopard2", "main")
CHILD_ENV = {
    "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1", "OMP_PROC_BIND": "TRUE",
    "OMP_PLACES": "cores", "PATH": "/usr/bin:/bin", "TZ": "UTC",
}
LINUX_F_ADD_SEALS = getattr(fcntl, "F_ADD_SEALS", 1033)
LINUX_F_GET_SEALS = getattr(fcntl, "F_GET_SEALS", 1034)
LINUX_F_SEAL_SEAL = getattr(fcntl, "F_SEAL_SEAL", 0x0001)
LINUX_F_SEAL_SHRINK = getattr(fcntl, "F_SEAL_SHRINK", 0x0002)
LINUX_F_SEAL_GROW = getattr(fcntl, "F_SEAL_GROW", 0x0004)
LINUX_F_SEAL_WRITE = getattr(fcntl, "F_SEAL_WRITE", 0x0008)
LINUX_MFD_CLOEXEC = getattr(os, "MFD_CLOEXEC", 0x0001)
LINUX_MFD_ALLOW_SEALING = getattr(os, "MFD_ALLOW_SEALING", 0x0002)
AUTHORITATIVE_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")
MAX_BENCHMARK_STDOUT_BYTES = 8 * 1024 * 1024
MAX_BENCHMARK_STDERR_BYTES = 1024 * 1024
MAX_HELPER_STDOUT_BYTES = 64 * 1024 * 1024
MAX_HELPER_STDERR_BYTES = 8 * 1024 * 1024
MAX_PERSISTED_CELL_BYTES = (
    len(ORDER) *
    (2 * MAX_BENCHMARK_STDOUT_BYTES +
     6 * MAX_BENCHMARK_STDERR_BYTES) +
    16 * 1024 * 1024
)
MAX_PERSISTED_MANIFEST_BYTES = 64 * 1024 * 1024
MAX_SOURCE_ATTESTATION_PREFLIGHTS = 1024
MASK64 = (1 << 64) - 1
CHILD_REAP_TIMEOUT_SECONDS = 5.0
MAX_SUPERVISOR_CONTROL_BYTES = 64 * 1024
BOUNDED_SUPERVISOR_MODE = "--internal-bounded-process-supervisor"
RUN_CONTRACT_SCHEMA_V4 = "leopard2-all-k-gap-contract/v4"
RUN_CONTRACT_SCHEMA_V5 = "leopard2-all-k-gap-contract/v5"
RUN_CONTRACT_SCHEMA_V6 = "leopard2-all-k-gap-contract/v6"
RUN_CONTRACT_SCHEMA = "leopard2-all-k-gap-contract/v7"
MANIFEST_SCHEMA_V4 = "leopard2-all-k-gap-manifest/v4"
MANIFEST_SCHEMA_V5 = "leopard2-all-k-gap-manifest/v5"
MANIFEST_SCHEMA_V6 = "leopard2-all-k-gap-manifest/v6"
MANIFEST_SCHEMA = "leopard2-all-k-gap-manifest/v7"
ALL_K_EVIDENCE_CONTRACTS = {
    RUN_CONTRACT_SCHEMA_V4: {
        "closure": PRODUCTION_BUILD_CLOSURE_SCHEMA_V1,
        "configuration": BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V2,
        "proof": REPRODUCIBLE_BUILD_PROOF_SCHEMA_V2,
        "replay_plan": LEGACY_REPLAY_RECIPE_SCHEMA,
        "replay_invocation": REPLAY_INVOCATION_SCHEMA_V1,
    },
    RUN_CONTRACT_SCHEMA: {
        "closure": PRODUCTION_BUILD_CLOSURE_SCHEMA,
        "configuration": BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
        "proof": REPRODUCIBLE_BUILD_PROOF_SCHEMA,
        "replay_plan": CANONICAL_REPLAY_RECIPE_SCHEMA,
        "replay_invocation": REPLAY_INVOCATION_SCHEMA,
    },
    RUN_CONTRACT_SCHEMA_V5: {
        "closure": PRODUCTION_BUILD_CLOSURE_SCHEMA,
        "configuration": BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3,
        "proof": REPRODUCIBLE_BUILD_PROOF_SCHEMA,
        "replay_plan": CANONICAL_REPLAY_RECIPE_SCHEMA,
        "replay_invocation": REPLAY_INVOCATION_SCHEMA,
    },
    RUN_CONTRACT_SCHEMA_V6: {
        "closure": PRODUCTION_BUILD_CLOSURE_SCHEMA,
        "configuration": BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4,
        "proof": REPRODUCIBLE_BUILD_PROOF_SCHEMA,
        "replay_plan": CANONICAL_REPLAY_RECIPE_SCHEMA,
        "replay_invocation": REPLAY_INVOCATION_SCHEMA,
    },
}
MANIFEST_TO_CONTRACT_SCHEMA = {
    MANIFEST_SCHEMA_V4: RUN_CONTRACT_SCHEMA_V4,
    MANIFEST_SCHEMA_V5: RUN_CONTRACT_SCHEMA_V5,
    MANIFEST_SCHEMA_V6: RUN_CONTRACT_SCHEMA_V6,
    MANIFEST_SCHEMA: RUN_CONTRACT_SCHEMA,
}
RUN_CONTRACT_KEYS = frozenset((
    "schema", "main_commit", "current_commit", "expected_main_sha256",
    "current_source_initial", "current_build_initial",
    "current_reproducible_build", "main_executable_initial",
    "current_executable_initial", "main_executable_snapshot",
    "current_executable_snapshot", "current_source_attestation_identity",
    "main_isa", "current_isa", "benchmark_lock", "allowed_cpus", "used_cpus",
    "workers", "order", "timeout_seconds", "with_current_legacy", "matrix",
    "measurement_note",
))
ALL_K_BUILD_CACHE_KEYS_V2 = frozenset((
    "CMAKE_BUILD_TYPE", "CMAKE_EXPORT_COMPILE_COMMANDS", "CMAKE_GENERATOR",
    "ENABLE_OPENMP", "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA",
    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256",
    "LEO2_BUILD_BENCHMARKS", "LEO2_BUILD_FUZZERS", "LEO2_ENABLE_CUDA",
    "LEO2_BUILD_ALLK_DIAGNOSTIC", "LEO2_BUILD_TESTS",
    "LEO2_BACKEND_VARIANT", "LEO2_BENCHMARK_GIT_EXECUTABLE",
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE",
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE", "LEOPARD_ENABLE_GF8",
    "LEOPARD_ENABLE_GF16", "LEO2_FLAG_FALIGN_FUNCTIONS_64",
    "LEO2_FLAG_MAVX2", "LEO2_FLAG_MNO_AVX512F",
    "LEO2_FLAG_MAVX512F", "LEO2_FLAG_MAVX512BW",
    "LEO2_FLAG_MAVX512VL", "LEO2_LOCATOR_GIT_EXECUTABLE",
    "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE", "CMAKE_AR",
    "CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER", "CMAKE_LINKER",
    "CMAKE_MAKE_PROGRAM", "CMAKE_RANLIB",
))
ALL_K_BUILD_CACHE_KEYS_V3 = frozenset((
    *ALL_K_BUILD_CACHE_KEYS_V2,
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
))
ALL_K_BUILD_CACHE_KEYS_V4 = frozenset((
    *ALL_K_BUILD_CACHE_KEYS_V3,
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR",
    "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE",
    "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING",
    "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING",
    "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING",
    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT",
))
ALL_K_BUILD_CACHE_KEYS = frozenset((
    *ALL_K_BUILD_CACHE_KEYS_V4,
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED",
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED",
))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def inherited_pass_fds(
        descriptors: Sequence[int], label: str) -> tuple[int, ...]:
    result = tuple(sorted(set(descriptors)))
    require(all(type(descriptor) is int and descriptor >= 0
                for descriptor in result),
            f"{label} inherited descriptor set is invalid")
    return result


def linux_memfd_create(name: str) -> int:
    flags = LINUX_MFD_CLOEXEC | LINUX_MFD_ALLOW_SEALING
    if hasattr(os, "memfd_create"):
        return os.memfd_create(name, flags)
    libc = ctypes.CDLL(None, use_errno=True)
    creator = getattr(libc, "memfd_create", None)
    require(creator is not None,
            "all-K executable snapshots require Linux memfd_create support")
    creator.argtypes = (ctypes.c_char_p, ctypes.c_uint)
    creator.restype = ctypes.c_int
    descriptor = creator(name.encode("utf-8"), flags)
    if descriptor < 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), name)
    return descriptor


def require_hex(value: str, label: str) -> str:
    require(isinstance(value, str) and len(value) == 40 and
            all(character in "0123456789abcdef" for character in value),
            f"{label} must be exactly 40 lowercase hexadecimal characters")
    return value


def require_sha256(value: str, label: str) -> str:
    require(isinstance(value, str) and len(value) == 64 and
            all(character in "0123456789abcdef" for character in value),
            f"{label} must be exactly 64 lowercase hexadecimal characters")
    return value


def inspect_isa_disassembly(disassembly: str) -> dict[str, int]:
    evex = 0
    ymm = 0
    for line in disassembly.splitlines():
        fields = line.lstrip().split()
        if fields and fields[0].endswith(":") and len(fields) > 1 and \
                fields[1] == "62":
            evex += 1
        if re.search(r"\bymm[0-9]+\b", line):
            ymm += 1
    return {"evex_prefixed_instruction_count": evex,
            "ymm_operand_instruction_count": ymm}


def inspect_pure_avx2(
        snapshot_descriptor: int, label: str, *,
        inherited_descriptors: Sequence[int] = ()) -> dict[str, Any]:
    require(type(snapshot_descriptor) is int and snapshot_descriptor >= 0,
            f"{label} snapshot descriptor is invalid")
    snapshot_path = Path(f"/proc/self/fd/{snapshot_descriptor}")
    pass_fds = inherited_pass_fds(
        (*inherited_descriptors, snapshot_descriptor), label)
    completed = run_helper_bounded(
        ["/usr/bin/objdump", "-d", "-M", "intel", str(snapshot_path)],
        timeout=180.0, max_stdout=MAX_HELPER_STDOUT_BYTES,
        max_stderr=MAX_HELPER_STDERR_BYTES,
        inherited_descriptors=pass_fds)
    stdout = completed.stdout.decode("utf-8", errors="strict")
    stderr = completed.stderr.decode("utf-8", errors="replace")
    require(completed.returncode == 0,
            f"objdump failed for {label}: {stderr.strip()}")
    result = inspect_isa_disassembly(stdout)
    require(result["evex_prefixed_instruction_count"] == 0,
            f"{label} contains EVEX-prefixed instructions")
    require(result["ymm_operand_instruction_count"] > 0,
            f"{label} contains no visible AVX2 YMM instructions")
    version = run_helper_bounded(
        ["/usr/bin/objdump", "--version"], timeout=30.0,
        max_stdout=1024 * 1024, max_stderr=1024 * 1024,
        inherited_descriptors=inherited_pass_fds(
            inherited_descriptors, f"{label} objdump version"))
    require(version.returncode == 0,
            f"objdump --version failed for {label}: " +
            version.stderr.decode("utf-8", errors="replace").strip())
    version_lines = version.stdout.decode(
        "utf-8", errors="strict").splitlines()
    require(bool(version_lines),
            f"objdump --version returned no text for {label}")
    return {**result, "objdump_version": version_lines[0]}


class IdentityBoundDescriptor:
    """Own one raw descriptor without ever closing a recycled descriptor.

    Raw descriptor numbers are process-global and can be reused by an
    interruption handler immediately after close(2).  Keep the acquired inode
    identity until close is observed, and retain ownership when a close fails
    before taking effect.  A second close attempt makes one-shot interruption
    failures self-cleaning while the post-failure probe prevents that retry
    from touching a recycled number.
    """

    def __init__(self, label: str):
        self.label = label
        self._descriptor: int | None = None
        self._identity: tuple[int, int, int, int] | None = None

    @staticmethod
    def _current_identity(descriptor: int) -> tuple[int, int, int, int]:
        metadata = os.fstat(descriptor)
        return (
            metadata.st_dev, metadata.st_ino,
            stat.S_IFMT(metadata.st_mode), metadata.st_rdev,
        )

    @staticmethod
    def _fallback_identity(descriptor: int) -> tuple[int, int, int, int]:
        """Read descriptor identity without fstat for interrupted cleanup."""
        metadata = os.stat(f"/proc/self/fd/{descriptor}")
        return (
            metadata.st_dev, metadata.st_ino,
            stat.S_IFMT(metadata.st_mode), metadata.st_rdev,
        )

    def _cleanup_identity(
            self, descriptor: int,
            ) -> tuple[tuple[int, int, int, int], BaseException | None]:
        try:
            return self._current_identity(descriptor), None
        except OSError as error:
            if error.errno == errno.EBADF:
                raise
            primary: BaseException = error
        except BaseException as error:
            primary = error
        try:
            return self._fallback_identity(descriptor), primary
        except BaseException as fallback_error:
            if hasattr(primary, "add_note"):
                primary.add_note(
                    f"{self.label} fallback identity probe also failed: "
                    f"{type(fallback_error).__name__}: {fallback_error}")
            raise primary

    @property
    def descriptor(self) -> int | None:
        return self._descriptor

    @property
    def closed(self) -> bool:
        return self._descriptor is None

    def owns(self, descriptor: int) -> bool:
        return self._descriptor == descriptor

    def adopt(self, descriptor: int) -> int:
        require(type(descriptor) is int and descriptor >= 0 and
                self._descriptor is None,
                f"{self.label} descriptor adoption is invalid")
        # Ownership transfers at entry.  Publish the descriptor before the
        # first fallible probe so an interrupted identity read can still be
        # cleaned through this object.
        self._descriptor = descriptor
        try:
            self._identity = self._current_identity(descriptor)
            return descriptor
        except BaseException as error:
            close_descriptor_preserving(
                self, error, f"{self.label} failed adoption")
            raise

    def open(
            self, path: str | os.PathLike[str], flags: int,
            mode: int = 0o777, *, dir_fd: int | None = None) -> int:
        require(self._descriptor is None,
                f"{self.label} descriptor is already open")
        descriptor = os.open(path, flags, mode, dir_fd=dir_fd)
        try:
            return self.adopt(descriptor)
        except BaseException as error:
            # adopt() normally completed cleanup itself.  Retry through the
            # same identity-bound owner if its first cleanup probe was also
            # interrupted; never fall back to an unverified raw close.
            if self.owns(descriptor):
                close_descriptor_preserving(
                    self, error, f"{self.label} failed open")
            raise

    def release(self) -> int:
        descriptor = self._descriptor
        require(descriptor is not None,
                f"{self.label} descriptor is not open")
        self._descriptor = None
        self._identity = None
        return descriptor

    def refresh_identity(self) -> None:
        descriptor = self._descriptor
        require(descriptor is not None,
                f"{self.label} descriptor is not open")
        self._identity = self._current_identity(descriptor)

    def _clear_if_current(
            self, descriptor: int,
            identity: tuple[int, int, int, int]) -> None:
        if self._descriptor == descriptor and self._identity == identity:
            self._descriptor = None
            self._identity = None

    def _probe_after_close_failure(
            self, descriptor: int,
            identity: tuple[int, int, int, int]) -> str:
        probe_error: BaseException | None = None
        try:
            current = self._current_identity(descriptor)
        except OSError as error:
            if error.errno == errno.EBADF:
                self._clear_if_current(descriptor, identity)
                return "closed"
            probe_error = error
        except BaseException as error:
            probe_error = error
        if probe_error is not None:
            try:
                current = self._fallback_identity(descriptor)
            except FileNotFoundError:
                self._clear_if_current(descriptor, identity)
                return "closed"
            except BaseException as fallback_error:
                if hasattr(probe_error, "add_note"):
                    probe_error.add_note(
                        f"{self.label} fallback close probe also failed: "
                        f"{type(fallback_error).__name__}: {fallback_error}")
                raise probe_error
        if current != identity:
            # The original close completed and an interruption handler reused
            # its number.  Drop the stale claim without touching the new file.
            self._clear_if_current(descriptor, identity)
            return "recycled"
        return "open"

    def close(self) -> None:
        descriptor = self._descriptor
        identity = self._identity
        identity_error: BaseException | None = None
        if descriptor is None:
            return
        if identity is None:
            try:
                identity, identity_error = \
                    self._cleanup_identity(descriptor)
            except OSError as error:
                if error.errno == errno.EBADF:
                    self._descriptor = None
                    return
                raise
            self._identity = identity
        try:
            current, current_error = self._cleanup_identity(descriptor)
            if identity_error is None:
                identity_error = current_error
        except OSError as error:
            if error.errno == errno.EBADF:
                self._clear_if_current(descriptor, identity)
                return
            raise
        if current != identity:
            self._clear_if_current(descriptor, identity)
            raise RuntimeError(
                f"{self.label} descriptor number was recycled before close")

        try:
            os.close(descriptor)
        except BaseException as first_error:
            state = self._probe_after_close_failure(descriptor, identity)
            if state == "open":
                try:
                    os.close(descriptor)
                except BaseException as retry_error:
                    retry_state = self._probe_after_close_failure(
                        descriptor, identity)
                    if hasattr(first_error, "add_note"):
                        first_error.add_note(
                            f"{self.label} descriptor close retry failed "
                            f"({retry_state}): {type(retry_error).__name__}: "
                            f"{retry_error}")
                    raise first_error
                else:
                    self._clear_if_current(descriptor, identity)
            raise first_error
        else:
            self._clear_if_current(descriptor, identity)
            if identity_error is not None:
                raise identity_error


def close_descriptor_preserving(
        owner: IdentityBoundDescriptor,
        primary: BaseException | None, label: str) -> None:
    """Close one owner without replacing an active primary exception."""
    try:
        owner.close()
    except BaseException as cleanup_error:
        if primary is None:
            raise
        if hasattr(primary, "add_note"):
            primary.add_note(
                f"{label} cleanup also failed: "
                f"{type(cleanup_error).__name__}: {cleanup_error}")


def close_descriptors_preserving(
        owners: Sequence[tuple[IdentityBoundDescriptor, str]],
        primary: BaseException | None) -> None:
    """Close every owner, preserving the first active failure."""
    first_cleanup: BaseException | None = None
    for owner, label in owners:
        try:
            owner.close()
        except BaseException as cleanup_error:
            active = primary if primary is not None else first_cleanup
            if active is None:
                first_cleanup = cleanup_error
            elif hasattr(active, "add_note"):
                active.add_note(
                    f"{label} cleanup also failed: "
                    f"{type(cleanup_error).__name__}: {cleanup_error}")
    if primary is None and first_cleanup is not None:
        raise first_cleanup


class CanonicalBenchmarkLock:
    """Fail-closed exclusive ownership of the shared benchmark lock inode."""

    def __init__(self, path: Path = AUTHORITATIVE_LOCK):
        self.path = Path(path)
        self._descriptor_owner = IdentityBoundDescriptor(
            "canonical benchmark lock")
        self.identity: dict[str, Any] | None = None

    @property
    def descriptor(self) -> int | None:
        return self._descriptor_owner.descriptor

    def validate_current(self) -> dict[str, Any]:
        require(self.descriptor is not None and self.identity is not None,
                "canonical benchmark lock is not held")
        try:
            descriptor = os.fstat(self.descriptor)
            path = os.lstat(self.path)
        except OSError as error:
            raise RuntimeError(
                f"canonical benchmark lock revalidation failed: {error}"
            ) from error
        current = {
            "device": descriptor.st_dev,
            "inode": descriptor.st_ino,
            "lock": "exclusive",
            "mode": stat.S_IMODE(descriptor.st_mode),
            "nlink": descriptor.st_nlink,
            "path": str(self.path),
            "uid": descriptor.st_uid,
        }
        require(
            stat.S_ISREG(descriptor.st_mode) and
            descriptor.st_uid == os.getuid() and
            descriptor.st_nlink == 1 and
            stat.S_IMODE(descriptor.st_mode) == 0o600 and
            (descriptor.st_dev, descriptor.st_ino) ==
            (path.st_dev, path.st_ino) and
            current == self.identity,
            "canonical benchmark lock path was replaced or changed")
        return dict(self.identity)

    def __enter__(self) -> "CanonicalBenchmarkLock":
        require(self.path == AUTHORITATIVE_LOCK,
                f"all-K runs require canonical lock {AUTHORITATIVE_LOCK}")
        require(hasattr(os, "O_NOFOLLOW"),
                "canonical benchmark locking requires O_NOFOLLOW")
        flags = os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW
        flags |= getattr(os, "O_CLOEXEC", 0)
        try:
            self._descriptor_owner.open(self.path, flags, 0o600)
        except OSError as error:
            raise RuntimeError(
                f"cannot open canonical benchmark lock {self.path}: {error}"
            ) from error
        try:
            descriptor = os.fstat(self.descriptor)
            path = os.lstat(self.path)
            require(
                stat.S_ISREG(descriptor.st_mode) and
                descriptor.st_uid == os.getuid() and
                descriptor.st_nlink == 1 and
                stat.S_IMODE(descriptor.st_mode) == 0o600 and
                (descriptor.st_dev, descriptor.st_ino) ==
                (path.st_dev, path.st_ino),
                "canonical benchmark lock has unsafe metadata")
            print(f"waiting for exclusive benchmark lock {self.path}",
                  flush=True)
            fcntl.flock(self.descriptor, fcntl.LOCK_EX)
            locked = os.fstat(self.descriptor)
            self.identity = {
                "device": locked.st_dev,
                "inode": locked.st_ino,
                "lock": "exclusive",
                "mode": stat.S_IMODE(locked.st_mode),
                "nlink": locked.st_nlink,
                "path": str(self.path),
                "uid": locked.st_uid,
            }
            self.validate_current()
            print("acquired benchmark lock", flush=True)
            return self
        except BaseException as error:
            close_descriptor_preserving(
                self._descriptor_owner, error,
                "canonical benchmark lock acquisition")
            if self._descriptor_owner.closed:
                self.identity = None
            raise

    def __exit__(self, exc_type, exc, traceback) -> None:
        validation_error: BaseException | None = None
        try:
            self.validate_current()
        except BaseException as error:
            validation_error = error
        # Do not issue LOCK_UN: all subprocesses receive a duplicate of this
        # open-file description.  Closing our descriptor lets a surviving
        # helper or nested build child retain the campaign lock until it exits.
        # LOCK_UN here would release the shared flock out from under it.
        close_error: BaseException | None = None
        try:
            self._descriptor_owner.close()
        except BaseException as error:
            close_error = error
        if self._descriptor_owner.closed:
            self.identity = None
        if close_error is not None:
            if validation_error is not None and \
                    hasattr(close_error, "add_note"):
                close_error.add_note(
                    "canonical benchmark lock validation also failed: "
                    f"{type(validation_error).__name__}: {validation_error}")
            raise close_error
        if validation_error is not None:
            raise validation_error


def benchmark_lock(path: Path) -> CanonicalBenchmarkLock:
    return CanonicalBenchmarkLock(path)


def git_identity(
        source_root: Path, requested_commit: str, *,
        inherited_descriptors: Sequence[int] = ()) -> dict[str, Any]:
    requested = require_hex(requested_commit, "current source commit")
    try:
        return git_capture.capture_git_identity(
            source_root, requested, require_detached=False,
            inherited_descriptors=inherited_descriptors)
    except git_capture.GitCaptureError as error:
        raise RuntimeError(str(error)) from error


def file_identity(path: Path, label: str) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    owner = IdentityBoundDescriptor(f"{label} identity source")
    descriptor = owner.open(resolved, flags)
    try:
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode), f"{label} is not a regular file")
        digest = hashlib.sha256()
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
        after = os.fstat(descriptor)
    finally:
        owner.close()
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_nlink",
        "st_size", "st_mtime_ns", "st_ctime_ns",
    )
    require(all(getattr(before, name) == getattr(after, name)
                for name in stable_fields),
            f"{label} changed while it was hashed")
    path_status = resolved.stat()
    require(all(getattr(after, name) == getattr(path_status, name)
                for name in stable_fields),
            f"{label} path changed while it was hashed")
    require(os.access(resolved, os.X_OK), f"{label} is not executable")
    return {
        "path": str(resolved), "sha256": digest.hexdigest(),
        "device": after.st_dev, "inode": after.st_ino,
        "mode": after.st_mode, "uid": after.st_uid, "gid": after.st_gid,
        "links": after.st_nlink, "size": after.st_size,
        "mtime_ns": after.st_mtime_ns, "ctime_ns": after.st_ctime_ns,
    }


def sealed_snapshot_identity(descriptor: int, label: str) -> dict[str, Any]:
    status = os.fstat(descriptor)
    require(stat.S_ISREG(status.st_mode), f"{label} snapshot is not regular")
    digest = hashlib.sha256()
    offset = 0
    while offset < status.st_size:
        chunk = os.pread(descriptor, min(1024 * 1024, status.st_size - offset),
                         offset)
        require(bool(chunk), f"{label} snapshot ended before its recorded size")
        digest.update(chunk)
        offset += len(chunk)
    seals = fcntl.fcntl(descriptor, LINUX_F_GET_SEALS)
    required_seals = (LINUX_F_SEAL_SEAL | LINUX_F_SEAL_SHRINK |
                      LINUX_F_SEAL_GROW | LINUX_F_SEAL_WRITE)
    require(seals & required_seals == required_seals,
            f"{label} snapshot is not immutably sealed")
    return {
        "kind": "linux-sealed-memfd-v1",
        "sha256": digest.hexdigest(), "size": status.st_size,
        "mode": status.st_mode, "seals": seals,
    }


def snapshot_executable(
        path: Path, label: str, *,
        descriptor_owner: IdentityBoundDescriptor | None = None) \
        -> tuple[dict[str, Any], int, dict[str, Any]]:
    require(sys.platform.startswith("linux") and hasattr(os, "pread"),
            "all-K executable snapshots require Linux sealed memfd support")
    source_identity = file_identity(path, label)
    source_owner = IdentityBoundDescriptor(f"{label} snapshot source")
    snapshot_owner = (
        descriptor_owner if descriptor_owner is not None else
        IdentityBoundDescriptor(f"{label} sealed snapshot"))
    require(snapshot_owner.closed,
            f"{label} sealed snapshot owner is already in use")
    source = source_owner.open(
        source_identity["path"],
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    snapshot_descriptor = -1
    failure: BaseException | None = None
    try:
        snapshot_descriptor = linux_memfd_create(
            "leopard2-allk-" + label.replace(" ", "-"))
        snapshot_owner.adopt(snapshot_descriptor)
        digest = hashlib.sha256()
        while True:
            chunk = os.read(source, 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            view = memoryview(chunk)
            while view:
                written = os.write(snapshot_descriptor, view)
                require(written > 0, f"{label} snapshot write made no progress")
                view = view[written:]
        require(digest.hexdigest() == source_identity["sha256"],
                f"{label} changed between identity capture and snapshot copy")
        require(os.pread(snapshot_descriptor, 4, 0) == b"\x7fELF",
                f"{label} is not an ELF executable")
        require(file_identity(Path(source_identity["path"]), label) ==
                source_identity,
                f"{label} changed while its executable snapshot was created")
        os.fchmod(snapshot_descriptor, 0o500)
        snapshot_owner.refresh_identity()
        required_seals = (LINUX_F_SEAL_SEAL | LINUX_F_SEAL_SHRINK |
                          LINUX_F_SEAL_GROW | LINUX_F_SEAL_WRITE)
        fcntl.fcntl(
            snapshot_descriptor, LINUX_F_ADD_SEALS, required_seals)
        snapshot_identity = sealed_snapshot_identity(
            snapshot_descriptor, label)
        require(snapshot_identity["sha256"] == source_identity["sha256"] and
                snapshot_identity["size"] == source_identity["size"],
                f"{label} sealed snapshot does not match its source identity")
        source_owner.close()
        if descriptor_owner is None:
            snapshot_descriptor = snapshot_owner.release()
        return source_identity, snapshot_descriptor, snapshot_identity
    except BaseException as error:
        failure = error
        close_descriptor_preserving(
            snapshot_owner, error, f"{label} sealed snapshot")
        raise
    finally:
        close_descriptor_preserving(
            source_owner, failure, f"{label} snapshot source")


def require_reproduced_snapshot_binding(
        build: Mapping[str, Any], proof: Mapping[str, Any],
        snapshot_source: Mapping[str, Any],
        snapshot: Mapping[str, Any]) -> None:
    """Bind the immutable timed bytes to the validated/reproduced closure."""
    closure_executable = build.get("executable")
    require(isinstance(closure_executable, dict),
            "current build closure has no executable identity")
    try:
        validate_reproducible_build_proof(
            proof, build, label="current all-K build")
    except Exception as error:
        raise RuntimeError(
            "current build has no valid reproducible executable proof: " +
            str(error)) from error
    require(snapshot_source == closure_executable,
            "Leopard2 snapshot source differs from the reproduced build "
            "closure executable")
    require(snapshot.get("sha256") == closure_executable.get("sha256") and
            snapshot.get("size") == closure_executable.get("size"),
            "Leopard2 timed snapshot bytes differ from the reproduced build "
            "closure executable")
    require(proof.get("executable_sha256") ==
            closure_executable.get("sha256"),
            "Leopard2 reproducible proof differs from the build closure "
            "executable")


def snapshot_reproduced_executable(
        path: Path, build: Mapping[str, Any],
        proof: Mapping[str, Any], *,
        descriptor_owner: IdentityBoundDescriptor | None = None,
        ) -> tuple[dict[str, Any], int, dict[str, Any]]:
    """Snapshot only the executable whose exact closure was reproduced."""
    owner = (
        descriptor_owner if descriptor_owner is not None else
        IdentityBoundDescriptor("reproduced Leopard2 snapshot"))
    source, descriptor, snapshot = snapshot_executable(
        path, "Leopard2 benchmark", descriptor_owner=owner)
    try:
        require_reproduced_snapshot_binding(
            build, proof, source, snapshot)
    except BaseException as error:
        close_descriptor_preserving(
            owner, error, "reproduced Leopard2 snapshot")
        raise
    if descriptor_owner is None:
        descriptor = owner.release()
    return source, descriptor, snapshot


class ExecutableSnapshotOwner:
    """Close every successfully captured executable memfd on every exit."""

    def __init__(self) -> None:
        self._owners: list[IdentityBoundDescriptor] = []

    @property
    def descriptors(self) -> list[int]:
        return [
            owner.descriptor for owner in self._owners
            if owner.descriptor is not None
        ]

    def _discard_closed(self, owner: IdentityBoundDescriptor) -> None:
        if owner.closed:
            self._owners = [
                retained for retained in self._owners
                if retained is not owner
            ]

    def capture(self, path: Path, label: str) \
            -> tuple[dict[str, Any], int, dict[str, Any]]:
        owner = IdentityBoundDescriptor(f"{label} owned snapshot")
        try:
            # Publish the empty owner before any descriptor can be acquired.
            self._owners.append(owner)
            return snapshot_executable(
                path, label, descriptor_owner=owner)
        except BaseException as error:
            close_descriptor_preserving(
                owner, error, f"{label} owned snapshot")
            self._discard_closed(owner)
            raise

    def capture_reproduced(
            self, path: Path, build: Mapping[str, Any],
            proof: Mapping[str, Any]) \
            -> tuple[dict[str, Any], int, dict[str, Any]]:
        owner = IdentityBoundDescriptor("owned reproduced Leopard2 snapshot")
        try:
            self._owners.append(owner)
            return snapshot_reproduced_executable(
                path, build, proof, descriptor_owner=owner)
        except BaseException as error:
            close_descriptor_preserving(
                owner, error, "owned reproduced Leopard2 snapshot")
            self._discard_closed(owner)
            raise

    def close(self) -> None:
        errors: list[BaseException] = []
        for owner in reversed(tuple(self._owners)):
            try:
                owner.close()
            except BaseException as error:
                errors.append(error)
            self._discard_closed(owner)
        if errors:
            raise RuntimeError(
                "cannot close executable snapshot: " +
                "; ".join(
                    f"{type(error).__name__}: {error}"
                    for error in errors)
            ) from errors[0]


def canonical_digest(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"),
        allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def canonical_equal(left: Any, right: Any) -> bool:
    """Compare JSON values without Python's bool/int/float coercions."""
    return canonical_digest(left) == canonical_digest(right)


class AtomicEvidenceDirectory:
    """Lifetime owner of one canonical all-K output directory inode."""

    def __init__(
            self, path: Path, descriptor_owner: IdentityBoundDescriptor,
            identity: tuple[int, int]) -> None:
        self.path = path
        require(not descriptor_owner.closed,
                "all-K evidence root descriptor is closed")
        self._descriptor_owner = descriptor_owner
        self.identity = identity

    @property
    def descriptor(self) -> int | None:
        return self._descriptor_owner.descriptor

    @staticmethod
    def _flags() -> int:
        require(hasattr(os, "O_NOFOLLOW") and hasattr(os, "O_DIRECTORY"),
                "secure evidence writes require O_NOFOLLOW and O_DIRECTORY")
        return (os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
                getattr(os, "O_CLOEXEC", 0))

    @staticmethod
    def _absolute(path: Path) -> Path:
        try:
            return Path(os.path.abspath(os.fspath(path)))
        except (OSError, TypeError, ValueError) as error:
            raise RuntimeError(
                f"evidence directory path is invalid: {error}") from error

    @classmethod
    def _open_absolute(cls, path: Path) -> IdentityBoundDescriptor:
        require(path.is_absolute(),
                "evidence directory path is not absolute")
        current = IdentityBoundDescriptor(
            "absolute evidence traversal root")
        child: IdentityBoundDescriptor | None = None
        try:
            current.open("/", cls._flags())
            for component in path.parts[1:]:
                require(component not in {"", ".", ".."} and
                        "/" not in component,
                        f"invalid evidence directory component: {component!r}")
                child = IdentityBoundDescriptor(
                    f"absolute evidence component {component}")
                child.open(
                    component, cls._flags(), dir_fd=current.descriptor)
                metadata = os.fstat(child.descriptor)
                entry = os.stat(
                    component, dir_fd=current.descriptor,
                    follow_symlinks=False)
                require(stat.S_ISDIR(metadata.st_mode) and
                        (metadata.st_dev, metadata.st_ino) ==
                        (entry.st_dev, entry.st_ino),
                        f"evidence directory component changed: {component}")
                current.close()
                current = child
                child = None
            return current
        except BaseException as error:
            if child is not None:
                close_descriptor_preserving(
                    child, error,
                    "absolute evidence child traversal")
            close_descriptor_preserving(
                current, error, "absolute evidence traversal")
            raise

    @classmethod
    def open_or_create(cls, requested: Path) -> "AtomicEvidenceDirectory":
        path = cls._absolute(requested)
        require(path != path.parent,
                "the filesystem root cannot be an evidence directory")
        current_owner: IdentityBoundDescriptor | None = \
            IdentityBoundDescriptor("canonical evidence traversal root")
        child_owner: IdentityBoundDescriptor | None = None
        try:
            current_owner.open("/", cls._flags())
            current = Path("/")
            for component in path.parts[1:]:
                created = False
                try:
                    child_owner = IdentityBoundDescriptor(
                        f"canonical evidence component {component}")
                    child_owner.open(
                        component, cls._flags(),
                        dir_fd=current_owner.descriptor)
                except FileNotFoundError:
                    child_owner = None
                    try:
                        os.mkdir(
                            component, 0o700,
                            dir_fd=current_owner.descriptor)
                    except FileExistsError as error:
                        raise RuntimeError(
                            f"evidence directory appeared during creation: "
                            f"{current / component}") from error
                    initial = os.stat(
                        component, dir_fd=current_owner.descriptor,
                        follow_symlinks=False)
                    require(stat.S_ISDIR(initial.st_mode) and
                            initial.st_uid == os.getuid(),
                            f"new evidence directory is unsafe: "
                            f"{current / component}")
                    # mkdir honors umask.  Establish the exact contract before
                    # opening a possibly mode-000 directory.
                    os.chmod(
                        component, 0o700, dir_fd=current_owner.descriptor,
                        follow_symlinks=False)
                    child_owner = IdentityBoundDescriptor(
                        f"new canonical evidence component {component}")
                    child_owner.open(
                        component, cls._flags(),
                        dir_fd=current_owner.descriptor)
                    created = True
                metadata = os.fstat(child_owner.descriptor)
                entry = os.stat(
                    component, dir_fd=current_owner.descriptor,
                    follow_symlinks=False)
                require(stat.S_ISDIR(metadata.st_mode) and
                        (metadata.st_dev, metadata.st_ino) ==
                        (entry.st_dev, entry.st_ino),
                        f"evidence directory entry changed: "
                        f"{current / component}")
                if created:
                    require(metadata.st_uid == os.getuid() and
                            stat.S_IMODE(metadata.st_mode) == 0o700 and
                            (metadata.st_dev, metadata.st_ino) ==
                            (initial.st_dev, initial.st_ino),
                            f"new evidence directory was replaced: "
                            f"{current / component}")
                    os.fsync(current_owner.descriptor)
                current_owner.close()
                current_owner = child_owner
                child_owner = None
                current /= component
            metadata = os.fstat(current_owner.descriptor)
            require(stat.S_ISDIR(metadata.st_mode) and
                    metadata.st_uid == os.getuid() and
                    stat.S_IMODE(metadata.st_mode) & 0o022 == 0,
                    f"evidence directory is not safely owned: {path}")
            # Previous all-K runs may have mode 0755.  They were not writable
            # by another user, so normalize them to the new exact mode.
            os.fchmod(current_owner.descriptor, 0o700)
            current_owner.refresh_identity()
            metadata = os.fstat(current_owner.descriptor)
            identity = (metadata.st_dev, metadata.st_ino)
            result = cls(path, current_owner, identity)
            result.validate_current()
            current_owner = None
            return result
        except OSError as error:
            raise RuntimeError(
                f"cannot open canonical evidence directory {path}: "
                f"{error}") from error
        finally:
            primary_error = sys.exc_info()[1]
            if child_owner is not None:
                close_descriptor_preserving(
                    child_owner, primary_error,
                    "canonical evidence child traversal")
            if current_owner is not None:
                close_descriptor_preserving(
                    current_owner, primary_error,
                    "canonical evidence traversal")

    def validate_current(self) -> None:
        require(self.descriptor is not None,
                "evidence directory descriptor is closed")
        rebound: IdentityBoundDescriptor | None = None
        try:
            metadata = os.fstat(self.descriptor)
            rebound = self._open_absolute(self.path)
            current = os.fstat(rebound.descriptor)
            require(stat.S_ISDIR(metadata.st_mode) and
                    metadata.st_uid == os.getuid() and
                    stat.S_IMODE(metadata.st_mode) == 0o700 and
                    (metadata.st_dev, metadata.st_ino) == self.identity ==
                    (current.st_dev, current.st_ino),
                    "evidence output directory was replaced or changed")
        except OSError as error:
            raise RuntimeError(
                f"cannot revalidate evidence output directory: {error}"
            ) from error
        finally:
            if rebound is not None:
                close_descriptor_preserving(
                    rebound, sys.exc_info()[1],
                    "revalidated evidence directory")

    @staticmethod
    def _parts(relative: str) -> tuple[str, ...]:
        require(isinstance(relative, str) and relative and
                not os.path.isabs(relative),
                "evidence destination is not relative")
        parts = tuple(Path(relative).parts)
        require(parts and all(
            component not in {"", ".", ".."} and "/" not in component
            for component in parts),
            "evidence destination is invalid")
        return parts

    def _open_parent(
            self, relative: str, create: bool,
            ) -> tuple[IdentityBoundDescriptor, str]:
        self.validate_current()
        parts = self._parts(relative)
        require(self.descriptor is not None,
                "evidence directory descriptor is closed")
        current = IdentityBoundDescriptor("held evidence parent root")
        child: IdentityBoundDescriptor | None = None
        try:
            current.adopt(os.dup(self.descriptor))
            for component in parts[:-1]:
                created = False
                try:
                    child = IdentityBoundDescriptor(
                        f"held evidence child {component}")
                    child.open(
                        component, self._flags(),
                        dir_fd=current.descriptor)
                except FileNotFoundError:
                    child = None
                    require(create,
                            f"evidence child directory is absent: {relative}")
                    try:
                        os.mkdir(
                            component, 0o700, dir_fd=current.descriptor)
                        created = True
                    except FileExistsError:
                        # Parallel all-K lanes may race to create "cells".
                        # The common validation below accepts only the exact
                        # owner-only directory that won.
                        created = False
                    if created:
                        initial = os.stat(
                            component, dir_fd=current.descriptor,
                            follow_symlinks=False)
                        require(stat.S_ISDIR(initial.st_mode) and
                                initial.st_uid == os.getuid(),
                                f"new evidence child directory is unsafe: "
                                f"{relative}")
                        os.chmod(
                            component, 0o700, dir_fd=current.descriptor,
                            follow_symlinks=False)
                    else:
                        observed = os.stat(
                            component, dir_fd=current.descriptor,
                            follow_symlinks=False)
                        require(stat.S_ISDIR(observed.st_mode) and
                                observed.st_uid == os.getuid() and
                                stat.S_IMODE(observed.st_mode) & 0o022 == 0,
                                f"racing evidence child directory is unsafe: "
                                f"{relative}")
                        os.chmod(
                            component, 0o700, dir_fd=current.descriptor,
                            follow_symlinks=False)
                    child = IdentityBoundDescriptor(
                        f"new held evidence child {component}")
                    child.open(
                        component, self._flags(),
                        dir_fd=current.descriptor)
                metadata = os.fstat(child.descriptor)
                entry = os.stat(
                    component, dir_fd=current.descriptor,
                    follow_symlinks=False)
                require(stat.S_ISDIR(metadata.st_mode) and
                        metadata.st_uid == os.getuid() and
                        stat.S_IMODE(metadata.st_mode) == 0o700 and
                        (metadata.st_dev, metadata.st_ino) ==
                        (entry.st_dev, entry.st_ino),
                        f"evidence child directory is unsafe: {relative}")
                if created:
                    require((metadata.st_dev, metadata.st_ino) ==
                            (initial.st_dev, initial.st_ino),
                            f"new evidence child directory was replaced: "
                            f"{relative}")
                    os.fsync(current.descriptor)
                current.close()
                current = child
                child = None
            return current, parts[-1]
        except BaseException as error:
            if child is not None:
                close_descriptor_preserving(
                    child, error, "held evidence child traversal")
            close_descriptor_preserving(
                current, error, "held evidence parent traversal")
            raise

    def _revalidate_parent(self, relative: str, descriptor: int) -> None:
        rebound, unused_name = self._open_parent(relative, create=False)
        del unused_name
        try:
            expected = os.fstat(descriptor)
            current = os.fstat(rebound.descriptor)
            require((expected.st_dev, expected.st_ino) ==
                    (current.st_dev, current.st_ino),
                    f"evidence destination parent was replaced: {relative}")
        finally:
            close_descriptor_preserving(
                rebound, sys.exc_info()[1],
                "revalidated evidence parent")

    def read_optional(
            self, relative: str, maximum_bytes: int, label: str, *,
            mutation_hook: Callable[[], None] | None = None,
            ) -> bytes | None:
        """Read one exact evidence inode through the lifetime-held root.

        The parent is created when absent so a first-run lookup and a resume
        lookup use the same componentwise no-follow traversal.  Existing
        evidence must satisfy the exact modes established by ``write_atomic``;
        two complete reads plus metadata/path rebinding reject in-place and
        A-to-B-to-A pathname substitution.
        """
        require(type(maximum_bytes) is int and maximum_bytes > 0,
                f"{label} byte bound is invalid")
        directory_owner, name = self._open_parent(relative, create=True)
        directory = directory_owner.descriptor
        file_owner = IdentityBoundDescriptor(f"held evidence {label}")
        try:
            flags = (os.O_RDONLY | os.O_NOFOLLOW |
                     getattr(os, "O_CLOEXEC", 0) |
                     getattr(os, "O_NONBLOCK", 0))
            try:
                file_owner.open(name, flags, dir_fd=directory)
            except FileNotFoundError:
                self._revalidate_parent(relative, directory)
                return None
            descriptor = file_owner.descriptor
            metadata = os.fstat(descriptor)
            current = os.stat(
                name, dir_fd=directory, follow_symlinks=False)
            identity = (
                metadata.st_dev, metadata.st_ino, metadata.st_mode,
                metadata.st_uid, metadata.st_gid, metadata.st_nlink,
                metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns,
            )
            current_identity = (
                current.st_dev, current.st_ino, current.st_mode,
                current.st_uid, current.st_gid, current.st_nlink,
                current.st_size, current.st_mtime_ns, current.st_ctime_ns,
            )
            require(stat.S_ISREG(metadata.st_mode) and
                    metadata.st_uid == os.getuid() and
                    metadata.st_nlink == 1 and
                    stat.S_IMODE(metadata.st_mode) == 0o600 and
                    metadata.st_size <= maximum_bytes and
                    identity == current_identity,
                    f"{label} is not exact owner-only evidence")
            values: list[bytes] = []
            for pass_index in range(2):
                os.lseek(descriptor, 0, os.SEEK_SET)
                payload = bytearray()
                hook_called = False
                while len(payload) <= maximum_bytes:
                    block = os.read(
                        descriptor,
                        min(1024 * 1024,
                            maximum_bytes + 1 - len(payload)))
                    if not block:
                        break
                    payload.extend(block)
                    if (mutation_hook is not None and pass_index == 0 and
                            not hook_called):
                        hook_called = True
                        mutation_hook()
                require(len(payload) <= maximum_bytes,
                        f"{label} exceeds {maximum_bytes} bytes")
                values.append(bytes(payload))
            final = os.fstat(descriptor)
            final_current = os.stat(
                name, dir_fd=directory, follow_symlinks=False)
            final_identity = (
                final.st_dev, final.st_ino, final.st_mode,
                final.st_uid, final.st_gid, final.st_nlink,
                final.st_size, final.st_mtime_ns, final.st_ctime_ns,
            )
            final_current_identity = (
                final_current.st_dev, final_current.st_ino,
                final_current.st_mode, final_current.st_uid,
                final_current.st_gid, final_current.st_nlink,
                final_current.st_size, final_current.st_mtime_ns,
                final_current.st_ctime_ns,
            )
            require(values[0] == values[1] and
                    len(values[0]) == metadata.st_size and
                    identity == final_identity == final_current_identity,
                    f"{label} changed while it was being read")
            self._revalidate_parent(relative, directory)
            return values[0]
        except OSError as error:
            raise RuntimeError(
                f"cannot read held evidence {label}: {error}") from error
        finally:
            primary_error = sys.exc_info()[1]
            close_descriptors_preserving((
                (file_owner, f"held evidence {label}"),
                (directory_owner, "held evidence parent directory"),
            ), primary_error)

    @staticmethod
    def _validate_destination(
            directory: int, name: str, path: Path) -> None:
        try:
            metadata = os.stat(
                name, dir_fd=directory, follow_symlinks=False)
        except FileNotFoundError:
            return
        require(stat.S_ISREG(metadata.st_mode) and
                metadata.st_uid == os.getuid() and metadata.st_nlink == 1 and
                stat.S_IMODE(metadata.st_mode) == 0o600,
                f"refusing unsafe evidence destination {path}")

    def write_atomic(self, relative: str, encoded: bytes) -> None:
        require(isinstance(encoded, bytes),
                "evidence payload is not bytes")
        directory_owner, name = self._open_parent(relative, create=True)
        directory = directory_owner.descriptor
        destination = self.path / relative
        temporary_name: str | None = None
        temporary_owner = IdentityBoundDescriptor(
            "secure evidence temporary")
        try:
            self._validate_destination(directory, name, destination)
            flags = (os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW |
                     getattr(os, "O_CLOEXEC", 0))
            for unused_attempt in range(128):
                candidate = (
                    f".{name}.tmp-{os.getpid()}-{secrets.token_hex(16)}")
                try:
                    temporary_owner.open(
                        candidate, flags, 0o600, dir_fd=directory)
                    temporary_name = candidate
                    break
                except FileExistsError:
                    continue
            require(temporary_owner.descriptor is not None and
                    temporary_name is not None,
                    f"cannot allocate secure evidence temporary in "
                    f"{destination.parent}")
            temporary_descriptor = temporary_owner.descriptor
            os.fchmod(temporary_descriptor, 0o600)
            view = memoryview(encoded)
            written = 0
            while written < len(view):
                try:
                    count = os.write(
                        temporary_descriptor, view[written:])
                except InterruptedError:
                    continue
                require(count > 0,
                        "secure evidence temporary write made no progress")
                written += count
            os.fsync(temporary_descriptor)
            temporary = os.fstat(temporary_descriptor)
            require(stat.S_ISREG(temporary.st_mode) and
                    temporary.st_uid == os.getuid() and
                    temporary.st_nlink == 1 and
                    temporary.st_size == len(encoded) and
                    stat.S_IMODE(temporary.st_mode) == 0o600,
                    "secure evidence temporary metadata changed")
            self._validate_destination(directory, name, destination)
            self._revalidate_parent(relative, directory)
            os.replace(
                temporary_name, name,
                src_dir_fd=directory, dst_dir_fd=directory)
            temporary_name = None
            published = os.stat(
                name, dir_fd=directory, follow_symlinks=False)
            require((published.st_dev, published.st_ino) ==
                    (temporary.st_dev, temporary.st_ino) and
                    stat.S_ISREG(published.st_mode) and
                    published.st_uid == os.getuid() and
                    published.st_nlink == 1 and
                    stat.S_IMODE(published.st_mode) == 0o600,
                    f"published evidence file changed: {relative}")
            os.fsync(directory)
            self._revalidate_parent(relative, directory)
        except OSError as error:
            raise RuntimeError(
                f"cannot atomically replace evidence file {destination}: "
                f"{error}") from error
        finally:
            primary_error = sys.exc_info()[1]
            cleanup_errors: list[tuple[str, BaseException]] = []
            try:
                temporary_owner.close()
            except BaseException as error:
                cleanup_errors.append(
                    ("close secure evidence temporary", error))
            if temporary_name is not None:
                try:
                    os.unlink(temporary_name, dir_fd=directory)
                except FileNotFoundError:
                    pass
                except BaseException as error:
                    cleanup_errors.append(
                        ("remove secure evidence temporary", error))
            try:
                directory_owner.close()
            except BaseException as error:
                cleanup_errors.append(
                    ("close held evidence parent directory", error))
            if cleanup_errors:
                cleanup_detail = "; ".join(
                    f"{operation}: {type(error).__name__}: {error}"
                    for operation, error in cleanup_errors)
                if primary_error is not None:
                    raise RuntimeError(
                        f"{primary_error}; atomic evidence cleanup also "
                        f"failed: {cleanup_detail}") from primary_error
                raise RuntimeError(
                    "atomic evidence cleanup failed: " + cleanup_detail
                ) from cleanup_errors[0][1]

    def close(self) -> None:
        self._descriptor_owner.close()


def _write_bytes_atomic(
        path: Path, encoded: bytes, *,
        evidence_directory: AtomicEvidenceDirectory | None = None) -> None:
    require(isinstance(path, Path) and isinstance(encoded, bytes),
            "evidence destination or payload is invalid")
    owned = evidence_directory
    close_owned = False
    if owned is None:
        owned = AtomicEvidenceDirectory.open_or_create(path.parent)
        close_owned = True
        relative = path.name
    else:
        absolute = AtomicEvidenceDirectory._absolute(path)
        try:
            relative = absolute.relative_to(owned.path).as_posix()
        except ValueError as error:
            raise RuntimeError(
                "evidence destination is outside the held output directory"
            ) from error
    try:
        owned.write_atomic(relative, encoded)
    finally:
        if close_owned:
            owned.close()


def write_json_atomic(
        path: Path, value: Any, *, pretty: bool = False,
        evidence_directory: AtomicEvidenceDirectory | None = None) -> None:
    if pretty:
        encoded = json.dumps(
            value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    else:
        encoded = json.dumps(
            value, sort_keys=True, allow_nan=False) + "\n"
    _write_bytes_atomic(
        path, encoded.encode("utf-8", errors="strict"),
        evidence_directory=evidence_directory)


def validate_run_contract_evidence(
        contract: Mapping[str, Any]) -> Mapping[str, Any]:
    """Bind each all-K outer schema to one exact replay-evidence generation."""
    require(isinstance(contract, Mapping),
            "all-K run contract is not an object")
    require(set(contract) == RUN_CONTRACT_KEYS,
            "all-K run contract fields differ")
    contract_schema = contract.get("schema")
    expected = ALL_K_EVIDENCE_CONTRACTS.get(contract_schema)
    require(expected is not None,
            "all-K run contract schema is unsupported")
    build = contract.get("current_build_initial")
    proof = contract.get("current_reproducible_build")
    require(isinstance(build, Mapping) and isinstance(proof, Mapping),
            "all-K run contract build/replay evidence is absent")
    cache = build.get("validated_cache")
    immutable = proof.get("immutable_replay")
    transport = (
        immutable.get("recipe_transport")
        if isinstance(immutable, Mapping) else None)
    invocations = (
        immutable.get("invocations")
        if isinstance(immutable, Mapping) else None)
    expected_cache_keys = {
        RUN_CONTRACT_SCHEMA_V4: ALL_K_BUILD_CACHE_KEYS_V2,
        RUN_CONTRACT_SCHEMA_V5: ALL_K_BUILD_CACHE_KEYS_V3,
        RUN_CONTRACT_SCHEMA_V6: ALL_K_BUILD_CACHE_KEYS_V4,
        RUN_CONTRACT_SCHEMA: ALL_K_BUILD_CACHE_KEYS,
    }[contract_schema]
    require(
        build.get("schema") == expected["closure"] and
        isinstance(cache, Mapping) and
        set(cache) == expected_cache_keys and
        cache.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") ==
            expected["configuration"] and
        proof.get("schema") == expected["proof"] and
        isinstance(transport, Mapping) and
        transport.get("schema") == expected["replay_plan"] and
        isinstance(invocations, Mapping) and
        set(invocations) == {"configure", "build"} and
        all(isinstance(invocation, Mapping) and
            invocation.get("schema") == expected["replay_invocation"]
            for invocation in invocations.values()),
        "all-K run contract replay-evidence schema tuple differs")
    require(
        cache.get("LEO2_BUILD_ALLK_DIAGNOSTIC") == "ON" and
        cache.get("LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN") == "OFF" and
        cache.get("LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE") == "OFF" and
        cache.get("LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE") in {"0", "1", "2"},
        "all-K run contract selector tuple differs")
    if contract_schema == RUN_CONTRACT_SCHEMA:
        require(
            cache.get("LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR") == "OFF" and
            cache.get("LEO2_EXPERIMENT_CAUCHY_LOG_REUSE") == "ON" and
            cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT") == "ON" and
            cache.get("LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING") == "ON" and
            cache.get("LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING") == "ON" and
            cache.get("LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING") == "ON" and
            cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED") == "OFF" and
            cache.get(
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED") == "OFF" and
            cache.get(
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT") == "ON",
            "current all-K run contract does not bind the production "
            "selector tuple")
    elif contract_schema == RUN_CONTRACT_SCHEMA_V6:
        require(
            cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT") == "ON" and
            cache.get("LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT") ==
                "ON" and
            cache.get("LEO2_EXPERIMENT_CAUCHY_LOG_REUSE") == "ON" and
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED" not in cache and
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED" not in cache,
            "v6 all-K run contract selector tuple differs")
    elif contract_schema == RUN_CONTRACT_SCHEMA_V5:
        require(
            cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT") == "OFF" and
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT" not in cache and
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE" not in cache,
            "v5 all-K run contract selector tuple differs")
    else:
        require(
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT" not in cache and
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT" not in cache and
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE" not in cache,
            "v4 all-K run contract contains an unversioned selector")
    try:
        validate_reproducible_build_proof(
            proof, build, label="all-K run contract")
    except Exception as error:
        raise RuntimeError(
            "all-K run contract reproducible-build proof is invalid: " +
            str(error)) from error
    return contract


def validate_manifest(
        document: Mapping[str, Any], contract: Mapping[str, Any],
        contract_digest: str, cells: Sequence["Cell"],
        current_source: Mapping[str, Any],
        current_snapshot: Mapping[str, Any],
        attestation_identity: Mapping[str, Any]) -> None:
    require(set(document) == {
        "schema", "run_contract", "run_contract_sha256", "cells", "completion",
        "current_source_attestation_preflights",
    }, "all-K manifest keys changed")
    validate_run_contract_evidence(contract)
    manifest_schema = document.get("schema")
    require(
        MANIFEST_TO_CONTRACT_SCHEMA.get(manifest_schema) ==
            contract.get("schema"),
        "all-K manifest/run-contract schema tuple mismatch")
    validate_run_contract_evidence(document.get("run_contract"))
    require(canonical_equal(document.get("run_contract"), contract),
            "existing all-K manifest run contract does not match this request")
    require(document.get("run_contract_sha256") == contract_digest,
            "existing all-K manifest contract digest mismatch")
    require(canonical_digest(document.get("run_contract")) == contract_digest,
            "existing all-K manifest contract bytes are internally inconsistent")
    require(canonical_equal(
        document.get("cells"), [dataclasses.asdict(cell) for cell in cells]),
            "existing all-K manifest cell matrix mismatch")
    require(canonical_equal(
        contract.get("matrix"), describe_cell_matrix(cells)),
        "all-K run contract matrix description does not match its cells")
    completion = document.get("completion")
    require(completion is None or isinstance(completion, dict),
            "all-K manifest completion must be null or an object")
    preflights = document.get("current_source_attestation_preflights")
    require(isinstance(preflights, list) and
            1 <= len(preflights) <= MAX_SOURCE_ATTESTATION_PREFLIGHTS,
            "all-K manifest source-attestation preflights are malformed")
    for index, preflight in enumerate(preflights):
        require(isinstance(preflight, dict),
                f"all-K manifest preflight {index} is not an object")
        require(canonical_equal(source_attestation_identity(
            preflight, current_source, current_snapshot),
            attestation_identity),
            f"all-K manifest preflight {index} identity mismatch")


@dataclasses.dataclass(frozen=True)
class Cell:
    identifier: str
    family: str
    k: int
    r: int
    shard_bytes: int
    losses: int
    redundancy_band: str
    loss_band: str
    seed: int
    iterations: int
    reuse: int
    warmup: int


def describe_cell_matrix(cells: Sequence[Cell]) -> dict[str, Any]:
    """Return the exact compact matrix description bound into the contract.

    Keep this derived from the cell list rather than from the default campaign
    constants.  Focused campaign wrappers deliberately replace ``make_cells``;
    a hard-coded description would then sign a matrix different from the one
    actually executed.
    """
    require(bool(cells), "all-K cell matrix is empty")
    gf8 = [cell for cell in cells if cell.family == "gf8-all-k"]
    gf16 = [cell for cell in cells if cell.family == "gf16-representative"]
    require(len(gf8) + len(gf16) == len(cells),
            "all-K cell matrix contains an unknown family")
    require(bool(gf8), "all-K cell matrix has no GF8 cells")
    gf8_k = sorted({cell.k for cell in gf8})
    require(gf8_k == list(range(gf8_k[0], gf8_k[-1] + 1)),
            "all-K GF8 K range is not contiguous")
    return {
        "gf8_K": [gf8_k[0], gf8_k[-1]],
        "gf8_shard_bytes": sorted({cell.shard_bytes for cell in gf8}),
        "gf16_K": sorted({cell.k for cell in gf16}),
        "gf16_shard_bytes":
            sorted({cell.shard_bytes for cell in gf16}),
        "gf8_only": not gf16,
        "cell_count": len(cells),
    }


def ceil_pow2(value: int) -> int:
    return 1 << (value - 1).bit_length()


def parent_for(k: int, r: int) -> tuple[int, int]:
    padded = ceil_pow2(r)
    return ceil_pow2(k + padded), padded


class XorShift64:
    def __init__(self, seed: int):
        self.state = seed if seed else 0x9E3779B97F4A7C15

    def next(self) -> int:
        value = self.state
        value ^= (value << 13) & MASK64
        value ^= value >> 7
        value ^= (value << 17) & MASK64
        self.state = value & MASK64
        return self.state


def expected_losses(cell: Cell) -> list[int]:
    order = list(range(cell.k))
    random = XorShift64(cell.seed ^ 0xD1B54A32D192ED03)
    for remaining in range(len(order), 1, -1):
        selected = random.next() % remaining
        order[remaining - 1], order[selected] = \
            order[selected], order[remaining - 1]
    return sorted(order[:cell.losses])


def gf8_r_values(k: int) -> list[tuple[str, int]]:
    max_padded = 1 << ((256 - k).bit_length() - 1)
    maximum = min(k, max_padded)
    values = (("low-R", 1), ("mid-R", max(1, (maximum + 1) // 2)),
              ("max-GF8-R", maximum))
    result: list[tuple[str, int]] = []
    seen: set[int] = set()
    for name, value in values:
        if value not in seen:
            seen.add(value)
            result.append((name, value))
    return result


GF16_K = (
    129, 130, 191, 192, 193, 200, 224, 225, 239, 240, 241, 248,
    249, 255, 256, 257, 300, 511, 512, 513, 1000, 2048, 4096,
)


def gf16_r_values(k: int) -> list[tuple[str, int]]:
    if k <= 255:
        gap = 256 - k
        first_forced = (1 << (gap.bit_length() - 1)) + 1
        values = (("first-GF16-R", first_forced),
                  ("mid-R", min(k, max(first_forced, (k + 3) // 4))),
                  ("high-R", min(k, 512)))
    else:
        values = (("low-R", 1), ("mid-R", min(k, max(2, k // 8))),
                  ("high-R", min(k, 512)))
    result: list[tuple[str, int]] = []
    seen: set[int] = set()
    for name, value in values:
        parent, _ = parent_for(k, value)
        if value <= k and parent > 256 and value not in seen:
            seen.add(value)
            result.append((name, value))
    return result


def make_cells(*, gf8_only: bool = False) -> list[Cell]:
    cells: list[Cell] = []
    seed = 0x38000000
    for k in range(1, 256):
        for redundancy_band, r in gf8_r_values(k):
            parent, _ = parent_for(k, r)
            require(parent <= 256,
                    "GF8 all-K matrix generated an out-of-field parent")
            for shard_bytes in (4096, 65536):
                for loss_band, losses in (("one-loss", 1), ("max-loss", r)):
                    if loss_band == "max-loss" and losses == 1:
                        continue
                    seed += 1
                    cells.append(Cell(
                        f"gf8-k{k}-r{r}-b{shard_bytes}-l{losses}", "gf8-all-k",
                        k, r, shard_bytes, losses, redundancy_band, loss_band,
                        seed, 5, 16 if shard_bytes == 4096 else 4, 1))
    if not gf8_only:
        for k in GF16_K:
            for redundancy_band, r in gf16_r_values(k):
                parent, _ = parent_for(k, r)
                require(256 < parent <= 65536,
                        "GF16 representative matrix generated an invalid parent")
                for shard_bytes in (512, 4096):
                    for loss_band, losses in (("one-loss", 1),
                                              ("max-loss", r)):
                        if loss_band == "max-loss" and losses == 1:
                            continue
                        seed += 1
                        cells.append(Cell(
                            f"gf16-k{k}-r{r}-b{shard_bytes}-l{losses}",
                            "gf16-representative", k, r, shard_bytes, losses,
                            redundancy_band, loss_band, seed, 5,
                            16 if shard_bytes == 512 else 8, 1))
    require(len({cell.identifier for cell in cells}) == len(cells),
            "all-K cell identifiers are not unique")
    return cells


def direct_locator_cutoff(field: str, parent: int) -> int:
    if field == "gf8":
        if parent <= 8: return parent
        if parent == 16: return 8
        if parent in (32, 128): return 9
        if parent == 64: return 8
        return 7
    if parent <= 32: return parent
    if parent == 64: return 34
    if parent == 128: return 24
    if parent in (256, 512): return 16
    return 14


def classify_paths(cell: Cell, result: Mapping[str, Any]) -> dict[str, Any]:
    resolved = result.get("resolved")
    require(isinstance(resolved, dict),
            "all-K current benchmark has no resolved codec geometry")
    profile = resolved.get("profile")
    field = resolved.get("field")
    backend = resolved.get("backend")
    parent = resolved.get("parent_count")
    padded = resolved.get("padded_side")
    expected_parent, expected_padded = parent_for(cell.k, cell.r)
    require(type(parent) is int and type(padded) is int and
            (parent, padded) == (expected_parent, expected_padded),
            "all-K current benchmark reported inconsistent parent geometry")
    require(profile == "legacy_high_v1",
            "all-K current benchmark selected a non-legacy-high profile")
    expected_field = "gf8" if cell.family == "gf8-all-k" else "gf16"
    require(field == expected_field,
            "all-K current benchmark selected the wrong finite field")
    require(backend == "avx2",
            "all-K current benchmark did not select explicit AVX2")
    decode = resolved.get("selected_decode_path")
    decode_rule = resolved.get("selected_decode_rule")
    require(decode in ("direct", "generic", "materialized", "tiled"),
            "all-K current benchmark did not report its selected decode path")
    require(isinstance(decode_rule, str) and decode_rule,
            "all-K current benchmark did not report its decode rule")
    translated_low_capable = ceil_pow2(cell.k) == padded
    allowed_rules = {
        "direct": {"direct"},
        "generic": {"balanced_generic"},
        "materialized": {
            "measured_materialized", "workspace_materialized",
        },
        "tiled": {"workspace_tiled"},
    }
    if translated_low_capable:
        allowed_rules["materialized"].add("translated_low")
        allowed_rules["tiled"].add("translated_low")
    require(decode_rule in allowed_rules[decode],
            "all-K current benchmark reported an incoherent decode path/rule")
    thread_count = resolved.get("thread_count")
    aligned_prefix = resolved.get("decode_aligned_prefix_bytes")
    tail = resolved.get("decode_tail_bytes")
    rounded = resolved.get("decode_rounded_bytes")
    multi_item_batch = resolved.get("decode_multi_item_batch")
    expected_tail = cell.shard_bytes & 63
    expected_aligned = cell.shard_bytes - expected_tail
    expected_rounded = expected_aligned + (64 if expected_tail else 0)
    require(type(thread_count) is int and thread_count == 1,
            "all-K current benchmark reported the wrong resolved thread count")
    require(type(aligned_prefix) is int and
            aligned_prefix == expected_aligned and
            type(tail) is int and tail == expected_tail and
            type(rounded) is int and rounded == expected_rounded,
            "all-K current benchmark reported inconsistent byte geometry")
    require(type(multi_item_batch) is bool and not multi_item_batch,
            "all-K current benchmark reported an inconsistent batch geometry")
    work_slots = resolved.get("decode_required_work_slots")
    if decode == "direct":
        expected_work_slots = 0
    elif decode in ("generic", "materialized"):
        expected_work_slots = parent
    elif decode_rule == "translated_low":
        expected_work_slots = parent
    else:
        expected_work_slots = 2 * padded + cell.losses
    require(type(work_slots) is int and work_slots == expected_work_slots,
            "all-K current benchmark reported inconsistent decode workspace")
    if padded == 1:
        encode = "direct-xor-single-parity"
    else:
        encode = "specialized-high-transform"
    if decode == "direct":
        workspace = "direct-bounded"
        locator = "none"
    elif decode == "generic":
        workspace = "materialized-N"
        locator = "active-parent-" + (
            "sparse-direct" if padded <= direct_locator_cutoff(field, parent)
            else "walsh")
    else:
        workspace = (
            "materialized-N"
            if decode == "materialized" or decode_rule == "translated_low"
            else "tiled-side-plus-losses")
        permanent_cached = padded > cell.r and \
            cell.r <= direct_locator_cutoff(field, parent)
        effective = cell.r if permanent_cached else padded
        locator = "active-parent-" + (
            "sparse-direct" if effective <= direct_locator_cutoff(field, parent)
            else "walsh")
    return {
        "profile": profile, "field": field, "backend": backend,
        "parent_count": parent, "padded_side": padded,
        "parent_inflation": parent / float(cell.k + cell.r),
        "encode_path": encode, "decode_path": decode,
        "decode_rule": decode_rule,
        "decode_workspace": workspace, "locator_setup": locator,
        "decode_required_work_slots": work_slots,
    }


def command(role: str, executable: Path, cell: Cell, cpu: int,
            with_current_legacy: bool) -> list[str]:
    common = [
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell.k), "--r", str(cell.r),
        "--bytes", str(cell.shard_bytes), "--loss", str(cell.losses),
        "--batch", "1", "--reuse", str(cell.reuse),
        "--iterations", str(cell.iterations), "--warmup", str(cell.warmup),
        "--threads", "1", "--seed", str(cell.seed),
    ]
    if role == "leopard2":
        field = "gf8" if cell.family == "gf8-all-k" else "gf16"
        common.extend(("--profile", "high", "--field", field,
                       "--backend", "avx2", "--retain-samples",
                       "--report-decode-path"))
        if not with_current_legacy:
            common.append("--skip-legacy")
    common.extend(("--json", "-"))
    return common


def source_attestation_command(executable: Path, cpu: int) -> list[str]:
    return [
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", "3", "--r", "2", "--bytes", "64", "--loss", "1",
        "--batch", "1", "--reuse", "1", "--iterations", "1",
        "--warmup", "0", "--threads", "1", "--seed", "7",
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
        "--json", "-",
    ]


def strict_json_bytes(payload: bytes, label: str) -> Any:
    require(isinstance(payload, bytes), f"{label} is not a byte string")
    try:
        text = payload.decode("utf-8", errors="strict")

        def reject_constant(value: str) -> Any:
            raise ValueError(f"non-finite JSON constant {value}")

        def reject_duplicate_keys(
                pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
            result: dict[str, Any] = {}
            for key, value in pairs:
                if key in result:
                    raise ValueError(f"duplicate JSON key {key!r}")
                result[key] = value
            return result

        def finite_float(value: str) -> float:
            result = float(value)
            if not math.isfinite(result):
                raise ValueError(f"non-finite JSON number {value}")
            return result

        result = json.loads(
            text, parse_constant=reject_constant,
            parse_float=finite_float,
            object_pairs_hook=reject_duplicate_keys)
        pending = [result]
        while pending:
            value = pending.pop()
            if type(value) is float and not math.isfinite(value):
                raise ValueError("non-finite decoded JSON number")
            if isinstance(value, dict):
                pending.extend(value.values())
            elif isinstance(value, list):
                pending.extend(value)
        return result
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError,
            OverflowError, RecursionError) as error:
        raise RuntimeError(f"{label} is not strict finite JSON: {error}") from error


def _open_existing_evidence_directory(
        path: Path, label: str,
        ) -> tuple[Path, IdentityBoundDescriptor, tuple[int, ...]]:
    require(hasattr(os, "O_NOFOLLOW") and hasattr(os, "O_DIRECTORY"),
            f"{label} requires O_NOFOLLOW and O_DIRECTORY")
    lexical = Path(os.path.abspath(os.fspath(path)))
    try:
        owner = AtomicEvidenceDirectory._open_absolute(lexical)
    except OSError as error:
        raise RuntimeError(
            f"cannot open {label} directory {lexical}: {error}"
        ) from error
    try:
        descriptor = owner.descriptor
        metadata = os.fstat(descriptor)
        identity = (
            metadata.st_dev, metadata.st_ino, metadata.st_mode,
            metadata.st_uid, metadata.st_gid, metadata.st_nlink,
            metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns,
        )
        require(stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o700,
                f"{label} directory is not exact owner-only evidence")
        return lexical, owner, identity
    except BaseException as error:
        close_descriptor_preserving(
            owner, error, f"{label} directory")
        raise


def _validate_existing_evidence_directory(
        path: Path, descriptor: int, identity: tuple[int, ...],
        label: str) -> None:
    rebound: IdentityBoundDescriptor | None = None
    try:
        metadata = os.fstat(descriptor)
        rebound = AtomicEvidenceDirectory._open_absolute(path)
        current = os.fstat(rebound.descriptor)
    except OSError as error:
        raise RuntimeError(
            f"cannot revalidate {label} directory: {error}"
        ) from error
    finally:
        if rebound is not None:
            close_descriptor_preserving(
                rebound, sys.exc_info()[1],
                f"revalidated {label} directory")
    metadata_identity = (
        metadata.st_dev, metadata.st_ino, metadata.st_mode,
        metadata.st_uid, metadata.st_gid, metadata.st_nlink,
        metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns,
    )
    current_identity = (
        current.st_dev, current.st_ino, current.st_mode,
        current.st_uid, current.st_gid, current.st_nlink,
        current.st_size, current.st_mtime_ns, current.st_ctime_ns,
    )
    require(stat.S_ISDIR(metadata.st_mode) and
            metadata.st_uid == os.getuid() and
            stat.S_IMODE(metadata.st_mode) == 0o700 and
            metadata_identity == identity == current_identity,
            f"{label} directory changed while reading evidence")


def load_strict_json(
        path: Path, label: str, maximum_bytes: int,
        *, mutation_hook: Callable[[], None] | None = None) -> Any:
    """Load persisted evidence without JSON normalization ambiguities."""
    require(type(maximum_bytes) is int and maximum_bytes > 0,
            f"{label} byte bound is invalid")
    require(isinstance(path, Path) and path.name not in {"", ".", ".."} and
            "/" not in path.name,
            f"{label} path is invalid")
    parent, directory_owner, directory_identity = \
        _open_existing_evidence_directory(path.parent, label)
    directory_descriptor = directory_owner.descriptor
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(
            path.name, flags, dir_fd=directory_descriptor)
    except OSError as error:
        close_descriptor_preserving(
            directory_owner, error, f"{label} directory")
        raise RuntimeError(f"cannot read {label} {path}: {error}") from error
    try:
        metadata = os.fstat(descriptor)
        current = os.stat(
            path.name, dir_fd=directory_descriptor,
            follow_symlinks=False)
        require(stat.S_ISREG(metadata.st_mode) and
                metadata.st_uid == os.getuid() and metadata.st_nlink == 1 and
                stat.S_IMODE(metadata.st_mode) == 0o600 and
                (metadata.st_dev, metadata.st_ino) ==
                (current.st_dev, current.st_ino),
                f"{label} is not exact owner-only evidence")
        require(metadata.st_size <= maximum_bytes,
                f"{label} exceeds {maximum_bytes} bytes")
        initial_identity = (
            metadata.st_dev, metadata.st_ino, metadata.st_mode,
            metadata.st_uid, metadata.st_gid, metadata.st_nlink,
            metadata.st_size,
            metadata.st_mtime_ns, metadata.st_ctime_ns,
        )
        payload = bytearray()
        hook_called = False
        while len(payload) <= maximum_bytes:
            chunk = os.read(
                descriptor, min(1024 * 1024, maximum_bytes + 1 - len(payload)))
            if not chunk:
                break
            payload.extend(chunk)
            if mutation_hook is not None and not hook_called:
                hook_called = True
                mutation_hook()
        require(len(payload) <= maximum_bytes,
                f"{label} exceeds {maximum_bytes} bytes")
        os.lseek(descriptor, 0, os.SEEK_SET)
        verification = bytearray()
        while len(verification) <= maximum_bytes:
            chunk = os.read(
                descriptor,
                min(1024 * 1024,
                    maximum_bytes + 1 - len(verification)))
            if not chunk:
                break
            verification.extend(chunk)
        require(len(verification) <= maximum_bytes and
                verification == payload,
                f"{label} changed while it was being read")
        final_metadata = os.fstat(descriptor)
        final_current = os.stat(
            path.name, dir_fd=directory_descriptor,
            follow_symlinks=False)
        final_identity = (
            final_metadata.st_dev, final_metadata.st_ino,
            final_metadata.st_mode, final_metadata.st_uid,
            final_metadata.st_gid, final_metadata.st_nlink,
            final_metadata.st_size, final_metadata.st_mtime_ns,
            final_metadata.st_ctime_ns,
        )
        require(len(payload) == metadata.st_size and
                final_identity == initial_identity and
                (final_current.st_dev, final_current.st_ino,
                 final_current.st_mode, final_current.st_uid,
                 final_current.st_gid, final_current.st_nlink,
                 final_current.st_size, final_current.st_mtime_ns,
                 final_current.st_ctime_ns) == initial_identity,
                f"{label} changed while it was being read")
        _validate_existing_evidence_directory(
            parent, directory_descriptor, directory_identity, label)
    except OSError as error:
        raise RuntimeError(f"cannot read {label} {path}: {error}") from error
    finally:
        primary_error = sys.exc_info()[1]
        cleanup_errors: list[tuple[str, BaseException]] = []
        try:
            os.close(descriptor)
        except BaseException as error:
            cleanup_errors.append((f"close {label}", error))
        try:
            directory_owner.close()
        except BaseException as error:
            cleanup_errors.append((f"close {label} directory", error))
        if cleanup_errors:
            detail = "; ".join(
                f"{operation}: {type(error).__name__}: {error}"
                for operation, error in cleanup_errors)
            if primary_error is not None:
                raise RuntimeError(
                    f"{primary_error}; evidence cleanup also failed: "
                    f"{detail}") from primary_error
            raise RuntimeError(
                "evidence cleanup failed: " + detail
            ) from cleanup_errors[0][1]
    return strict_json_bytes(bytes(payload), label)


def _write_descriptor_all(descriptor: int, value: bytes) -> None:
    offset = 0
    while offset < len(value):
        try:
            written = os.write(descriptor, value[offset:])
        except InterruptedError:
            continue
        require(written > 0, "bounded-process supervisor write made no progress")
        offset += written


def _bounded_process_supervisor(arguments: Sequence[str]) -> int:
    """Single-threaded exec boundary for the audited descendant owner."""
    require(len(arguments) == 8,
            "bounded-process supervisor argument count changed")
    try:
        command_value = json.loads(arguments[0])
        timeout = float(arguments[1])
        max_stdout = int(arguments[2])
        max_stderr = int(arguments[3])
        environment_value = json.loads(arguments[4])
        inherited_value = json.loads(arguments[5])
        control_descriptor = int(arguments[6])
        launch_descriptor = int(arguments[7])
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise RuntimeError(
            f"bounded-process supervisor arguments are invalid: {error}"
        ) from error
    require(isinstance(command_value, list) and command_value and
            all(isinstance(item, str) and item for item in command_value),
            "bounded-process supervisor command is invalid")
    require(timeout > 0 and math.isfinite(timeout),
            "bounded-process supervisor timeout is invalid")
    require(type(max_stdout) is int and
            0 <= max_stdout <= MAX_HELPER_STDOUT_BYTES and
            type(max_stderr) is int and
            0 <= max_stderr <= MAX_HELPER_STDERR_BYTES,
            "bounded-process supervisor output limit is invalid")
    require(isinstance(environment_value, dict) and
            1 <= len(environment_value) <= 128 and
            all(isinstance(name, str) and name and
                isinstance(value, str) and "\0" not in name + value
                for name, value in environment_value.items()),
            "bounded-process supervisor environment is invalid")
    require(isinstance(inherited_value, list) and
            len(inherited_value) <= 128 and
            all(type(descriptor) is int and descriptor >= 0
                for descriptor in inherited_value) and
            all(type(descriptor) is int and descriptor >= 0 for descriptor in (
                control_descriptor, launch_descriptor)),
        "bounded-process supervisor descriptor is invalid")
    inherited_descriptors = inherited_pass_fds(
        inherited_value, "bounded-process supervisor")

    # The parent binds this process to an exact PID/start-time identity before
    # releasing the gate.  Thus even an identity-capture failure cannot race a
    # benchmark launch or leave an unattributable detached descendant.
    try:
        launch = os.read(launch_descriptor, 1)
    finally:
        os.close(launch_descriptor)
    require(launch == b"\x01",
            "bounded-process supervisor launch was not authorized")

    # Import only after the exec boundary.  run_abba's audited containment
    # requires a single-threaded owner; the all-K coordinator deliberately
    # uses many worker threads, so invoking it in-process would be unsound.
    module_directory = str(Path(__file__).resolve().parent)
    inserted = module_directory not in sys.path
    if inserted:
        sys.path.insert(0, module_directory)
    try:
        from run_abba import (  # pylint: disable=import-outside-toplevel
            EvidenceError as ContainmentError,
            run_process_bounded as audited_run_process_bounded,
        )
    finally:
        if inserted:
            sys.path.remove(module_directory)

    status: dict[str, Any]
    try:
        completed = audited_run_process_bounded(
            command_value, environment=environment_value, timeout=timeout,
            max_stdout=max_stdout, max_stderr=max_stderr,
            inherited_descriptors=inherited_descriptors)
        _write_descriptor_all(1, completed.stdout)
        _write_descriptor_all(2, completed.stderr)
        status = {"status": "ok", "returncode": completed.returncode}
    except ContainmentError as error:
        message = str(error)
        status = {
            "status": (
                "timeout" if message.startswith("command exceeded ")
                else "error"),
            "message": message,
        }
    except Exception as error:  # pragma: no cover - emergency protocol guard
        status = {
            "status": "error",
            "message": f"{type(error).__name__}: {error}",
        }
    encoded = json.dumps(
        status, sort_keys=True, separators=(",", ":"),
        allow_nan=False).encode("utf-8")
    require(len(encoded) <= MAX_SUPERVISOR_CONTROL_BYTES,
            "bounded-process supervisor control record is oversized")
    try:
        _write_descriptor_all(control_descriptor, encoded)
    finally:
        os.close(control_descriptor)
    return 0


def _open_supervisor_proc_directory(pid: int) -> int:
    require(type(pid) is int and pid > 0,
            "bounded-process supervisor tree PID is invalid")
    required = ("O_DIRECTORY", "O_NOFOLLOW")
    require(all(hasattr(os, name) for name in required),
            "bounded-process supervisor identity requires no-follow procfs")
    flags = (os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
             getattr(os, "O_CLOEXEC", 0))
    descriptor = os.open(f"/proc/{pid}", flags)
    try:
        metadata = os.fstat(descriptor)
        current = os.stat(f"/proc/{pid}", follow_symlinks=False)
        require(stat.S_ISDIR(metadata.st_mode) and
                (metadata.st_dev, metadata.st_ino) ==
                (current.st_dev, current.st_ino),
                f"bounded-process supervisor tree PID {pid} changed while "
                "opening procfs")
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _read_supervisor_proc_record(
        pid: int, proc_descriptor: int,
        ) -> tuple[int, int, str, int, int] | None:
    """Read one process through a held proc-directory inode."""
    stat_descriptor: int | None = None
    try:
        flags = (os.O_RDONLY | os.O_NOFOLLOW |
                 getattr(os, "O_CLOEXEC", 0) |
                 getattr(os, "O_NONBLOCK", 0))
        stat_descriptor = os.open(
            "stat", flags, dir_fd=proc_descriptor)
        chunks = bytearray()
        while len(chunks) <= 64 * 1024:
            block = os.read(
                stat_descriptor, min(4096, 64 * 1024 + 1 - len(chunks)))
            if not block:
                break
            chunks.extend(block)
        require(len(chunks) <= 64 * 1024,
                f"bounded-process supervisor tree PID {pid} stat is oversized")
        payload = bytes(chunks)
    except OSError as error:
        if error.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise RuntimeError(
            f"cannot inspect bounded-process supervisor tree PID {pid}: "
            f"{error}") from error
    finally:
        if stat_descriptor is not None:
            os.close(stat_descriptor)
    closing = payload.rfind(b")")
    require(closing > 0 and closing + 2 < len(payload),
            f"bounded-process supervisor tree PID {pid} has malformed stat")
    fields = payload[closing + 2:].split()
    require(len(fields) >= 20,
            f"bounded-process supervisor tree PID {pid} has truncated stat")
    try:
        state = fields[0].decode("ascii")
        parent = int(fields[1])
        starttime = int(fields[19])
    except (UnicodeDecodeError, ValueError) as error:
        raise RuntimeError(
            f"bounded-process supervisor tree PID {pid} has invalid stat"
        ) from error
    require(len(state) == 1 and parent >= 0 and starttime >= 0,
            f"bounded-process supervisor tree PID {pid} has invalid identity")
    try:
        metadata = os.fstat(proc_descriptor)
        current = os.stat(f"/proc/{pid}", follow_symlinks=False)
    except OSError as error:
        if error.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise RuntimeError(
            f"cannot rebind bounded-process supervisor tree PID {pid}: "
            f"{error}") from error
    if ((metadata.st_dev, metadata.st_ino) !=
            (current.st_dev, current.st_ino)):
        return None
    return parent, starttime, state, metadata.st_dev, metadata.st_ino


def _supervisor_proc_record(
        pid: int) -> tuple[int, int, str, int, int] | None:
    """Return ppid/start/state plus the exact proc-directory inode."""
    descriptor: int | None = None
    try:
        descriptor = _open_supervisor_proc_directory(pid)
        return _read_supervisor_proc_record(pid, descriptor)
    except OSError as error:
        if error.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise RuntimeError(
            f"cannot inspect bounded-process supervisor tree PID {pid}: "
            f"{error}") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _supervisor_proc_snapshot(
) -> dict[int, tuple[int, int, str, int, int]]:
    """Best-effort procfs snapshot; owned descendants remain mandatory."""
    try:
        names = os.listdir("/proc")
    except OSError as error:
        raise RuntimeError(
            f"cannot enumerate bounded-process supervisor tree: {error}"
        ) from error
    result: dict[int, tuple[int, int, str]] = {}
    for name in names:
        if not name.isascii() or not name.isdigit():
            continue
        pid = int(name)
        try:
            record = _supervisor_proc_record(pid)
        except RuntimeError:
            try:
                owner = os.stat(
                    f"/proc/{pid}", follow_symlinks=False).st_uid
            except OSError:
                continue
            if owner == os.getuid():
                raise
            continue
        if record is not None:
            result[pid] = record
    return result


def _open_pidfd(pid: int) -> int:
    require(type(pid) is int and pid > 0,
            "pidfd target is invalid")
    if hasattr(os, "pidfd_open"):
        return os.pidfd_open(pid, 0)
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_open
    except (AttributeError, OSError) as error:
        raise RuntimeError(
            f"bounded-process supervisor cleanup requires pidfd_open: {error}"
        ) from error
    ctypes.set_errno(0)
    descriptor = function(ctypes.c_int(pid), ctypes.c_uint(0))
    if descriptor < 0:
        number = ctypes.get_errno()
        raise OSError(number or errno.EPERM,
                      os.strerror(number or errno.EPERM), str(pid))
    return descriptor


def _send_pidfd_signal(descriptor: int, signal_number: int) -> None:
    if hasattr(signal, "pidfd_send_signal"):
        try:
            signal.pidfd_send_signal(descriptor, signal_number)
        except ProcessLookupError:
            pass
        return
    try:
        send_signal = ctypes.CDLL(None, use_errno=True).pidfd_send_signal
    except (AttributeError, OSError) as error:
        raise RuntimeError(
            f"bounded-process supervisor cleanup requires "
            f"pidfd_send_signal: {error}") from error
    ctypes.set_errno(0)
    status = send_signal(
        ctypes.c_int(descriptor), ctypes.c_int(signal_number), None,
        ctypes.c_uint(0))
    if status != 0:
        number = ctypes.get_errno()
        if number != errno.ESRCH:
            raise RuntimeError(
                "cannot signal bounded-process supervisor pidfd: " +
                os.strerror(number or errno.EPERM))


class RetainedPidfds:
    """Own pidfds bound to a held proc inode across pidfd_open."""

    def __init__(self) -> None:
        self.descriptors: dict[
            tuple[int, int, int, int], IdentityBoundDescriptor] = {}

    def retain(
            self, pid: int,
            expected_identity: tuple[int, int, int, int] | None = None,
            ) -> tuple[int, int, int, int] | None:
        if expected_identity is not None:
            require(isinstance(expected_identity, tuple) and
                    len(expected_identity) == 4 and
                    all(type(value) is int for value in expected_identity) and
                    expected_identity[0] == pid,
                    "expected bounded-process identity is invalid")
            existing = self.descriptors.get(expected_identity)
            if existing is not None and not existing.closed:
                return expected_identity
            if existing is not None:
                self.descriptors.pop(expected_identity, None)
        proc_owner = IdentityBoundDescriptor(
            f"bounded-process PID {pid} proc directory")
        pidfd_owner: IdentityBoundDescriptor | None = None
        identity: tuple[int, int, int, int] | None = None
        failure: BaseException | None = None
        try:
            try:
                proc_owner.adopt(
                    _open_supervisor_proc_directory(pid))
            except OSError as error:
                if error.errno in (errno.ENOENT, errno.ESRCH):
                    return None
                raise RuntimeError(
                    f"cannot retain bounded-process supervisor tree PID "
                    f"{pid}: {error}") from error
            proc_descriptor = proc_owner.descriptor
            require(proc_descriptor is not None,
                    "held proc descriptor was lost")
            before = _read_supervisor_proc_record(pid, proc_descriptor)
            if before is None:
                return None
            identity = (pid, before[1], before[3], before[4])
            if expected_identity is not None and identity != expected_identity:
                return None
            existing = self.descriptors.get(identity)
            if existing is not None and not existing.closed:
                return identity
            if existing is not None:
                self.descriptors.pop(identity, None)

            pidfd_owner = IdentityBoundDescriptor(
                f"bounded-process PID {pid} pidfd")
            try:
                pidfd_owner.adopt(_open_pidfd(pid))
            except OSError as error:
                if error.errno == errno.ESRCH:
                    return None
                raise RuntimeError(
                    f"cannot retain bounded-process supervisor tree PID "
                    f"{pid}: {error}") from error
            after = _read_supervisor_proc_record(pid, proc_descriptor)
            if after is None:
                return None
            after_identity = (
                pid, after[1], after[3], after[4])
            if after_identity != identity or (
                    expected_identity is not None and
                    after_identity != expected_identity):
                return None
            existing = self.descriptors.get(identity)
            if existing is not None and not existing.closed:
                return identity
            if existing is not None:
                self.descriptors.pop(identity, None)
            # The owner object, rather than a naked fd number, is the transfer
            # unit.  If an exception lands after insertion, the registry and
            # the local variable still name the same safe owner.
            self.descriptors[identity] = pidfd_owner
            return identity
        except BaseException as error:
            failure = error
            raise
        finally:
            cleanup_errors: list[tuple[str, BaseException]] = []
            try:
                proc_owner.close()
            except BaseException as cleanup_error:
                cleanup_errors.append(("proc directory", cleanup_error))
            if pidfd_owner is not None and (
                    identity is None or
                    self.descriptors.get(identity) is not pidfd_owner):
                try:
                    pidfd_owner.close()
                except BaseException as cleanup_error:
                    cleanup_errors.append(
                        ("untransferred pidfd", cleanup_error))
            if cleanup_errors:
                detail = "; ".join(
                    f"{role}: {type(error).__name__}: {error}"
                    for role, error in cleanup_errors)
                if failure is not None:
                    if hasattr(failure, "add_note"):
                        failure.add_note(
                            f"bounded-process PID {pid} cleanup also "
                            f"failed: {detail}")
                else:
                    raise RuntimeError(
                        f"bounded-process PID {pid} cleanup failed: "
                        f"{detail}") from cleanup_errors[0][1]

    def signal(
            self, identity: tuple[int, int, int, int],
            signal_number: int) -> None:
        owner = self.descriptors.get(identity)
        require(owner is not None and owner.descriptor is not None,
                f"no retained pidfd for process identity {identity}")
        _send_pidfd_signal(owner.descriptor, signal_number)

    def close(self) -> None:
        errors: list[BaseException] = []
        for identity, owner in tuple(self.descriptors.items()):
            try:
                owner.close()
            except BaseException as error:
                errors.append(error)
            if owner.closed and self.descriptors.get(identity) is owner:
                self.descriptors.pop(identity, None)
        if errors:
            raise RuntimeError(
                "cannot close retained pidfd: " +
                "; ".join(
                    f"{type(error).__name__}: {error}"
                    for error in errors)
            ) from errors[0]


def _pidfd_signal_exact(
        pidfds: RetainedPidfds, identity: tuple[int, int, int, int],
        signal_number: int) -> None:
    """Signal only a pidfd retained while this identity was observed."""
    pidfds.signal(identity, signal_number)


def _supervisor_descendants(
        root: tuple[int, int, int, int],
        snapshot: Mapping[int, tuple[int, int, str, int, int]],
        known: set[tuple[int, int, int, int]],
        ) -> set[tuple[int, int, int, int]]:
    """Find the exact root tree plus previously observed reparented members."""
    root_pid, root_starttime, root_proc_dev, root_proc_ino = root
    result = {
        identity for identity in known
        if (identity[0] in snapshot and
            (snapshot[identity[0]][1],
             snapshot[identity[0]][3],
             snapshot[identity[0]][4]) ==
            (identity[1], identity[2], identity[3]))
    }
    root_record = snapshot.get(root_pid)
    if (root_record is not None and
            (root_record[1], root_record[3], root_record[4]) ==
            (root_starttime, root_proc_dev, root_proc_ino)):
        result.add(root)
    changed = True
    while changed:
        changed = False
        parents = {identity[0] for identity in result}
        for pid, record in snapshot.items():
            identity = (pid, record[1], record[3], record[4])
            if record[0] in parents and identity not in result:
                result.add(identity)
                changed = True
    known.update(result)
    return result


def terminate_supervisor_tree_bounded(
        supervisor: subprocess.Popen[bytes],
        root: tuple[int, int, int, int],
        timeout: float = CHILD_REAP_TIMEOUT_SECONDS,
        pidfds: RetainedPidfds | None = None) -> None:
    owned = pidfds is None
    registry = pidfds if pidfds is not None else RetainedPidfds()
    try:
        require(registry.retain(root[0], root) == root,
                "bounded-process supervisor root identity changed before "
                "pidfd retention")
        _terminate_supervisor_tree_with_pidfds(
            supervisor, root, timeout, registry)
    finally:
        if owned:
            registry.close()


def _terminate_supervisor_tree_with_pidfds(
        supervisor: subprocess.Popen[bytes],
        root: tuple[int, int, int, int], timeout: float,
        pidfds: RetainedPidfds) -> None:
    """Freeze and kill a wedged supervisor and every separate-session child.

    The supervisor is itself the audited child subreaper.  It remains alive
    while this fallback discovers and stops its descendants, so double-forked
    or setsid children cannot lose their ancestry before they are captured.
    """
    require(supervisor.pid == root[0] and timeout > 0,
            "bounded-process supervisor cleanup identity is invalid")
    deadline = time.monotonic() + timeout
    known: set[tuple[int, int, int, int]] = {root}
    stable_scans = 0
    while stable_scans < 2:
        snapshot = _supervisor_proc_snapshot()
        targets = _supervisor_descendants(root, snapshot, known)
        # Freeze the subreaper first, then every descendant.  Repeated scans
        # close the fork race before any identity is killed or reparented.
        ordered = ([root] if root in targets else []) + sorted(
            targets - {root})
        for identity in ordered:
            retained = pidfds.retain(identity[0], identity)
            if retained is not None:
                _pidfd_signal_exact(pidfds, retained, signal.SIGSTOP)
            else:
                known.discard(identity)
        after = _supervisor_proc_snapshot()
        current = _supervisor_descendants(root, after, known)
        all_stopped = all(
            after.get(identity[0]) is not None and
            (after[identity[0]][1],
             after[identity[0]][3],
             after[identity[0]][4]) ==
            (identity[1], identity[2], identity[3]) and
            after[identity[0]][2] in {"T", "t", "Z", "X", "x"}
            for identity in current)
        if current == targets and all_stopped:
            stable_scans += 1
        else:
            stable_scans = 0
        require(time.monotonic() < deadline,
                "cannot freeze bounded-process supervisor tree")
        time.sleep(0.01)

    # Kill descendants before the root so their ancestry remains attributable
    # throughout the complete capture.  PID handles prevent reuse races.
    snapshot = _supervisor_proc_snapshot()
    targets = _supervisor_descendants(root, snapshot, known)
    for identity in sorted(targets - {root}, reverse=True):
        retained = pidfds.retain(identity[0], identity)
        if retained is not None:
            _pidfd_signal_exact(pidfds, retained, signal.SIGKILL)
    if root in targets:
        _pidfd_signal_exact(pidfds, root, signal.SIGKILL)
    remaining = max(0.001, deadline - time.monotonic())
    try:
        supervisor.wait(timeout=remaining)
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            "bounded-process supervisor could not be reaped after SIGKILL"
        ) from error

    while True:
        live = []
        for identity in known:
            record = _supervisor_proc_record(identity[0])
            if (record is not None and
                    (record[1], record[3], record[4]) ==
                    (identity[1], identity[2], identity[3]) and
                    record[2] not in {"Z", "X", "x"}):
                live.append(identity)
        if not live:
            return
        require(time.monotonic() < deadline,
                "bounded-process supervisor left live separate-session "
                "descendants")
        for identity in live:
            retained = pidfds.retain(identity[0], identity)
            if retained is not None:
                _pidfd_signal_exact(pidfds, retained, signal.SIGKILL)
        time.sleep(0.01)


def _close_supervisor_capture(
        supervisor: subprocess.Popen[bytes], control_read: int) -> None:
    """Close every coordinator-side capture descriptor, including partial use."""
    errors: list[OSError] = []
    for stream in (supervisor.stdout, supervisor.stderr):
        if stream is None:
            continue
        try:
            stream.close()
        except OSError as error:
            errors.append(error)
    try:
        os.close(control_read)
    except OSError as error:
        if error.errno != errno.EBADF:
            errors.append(error)
    if errors:
        raise RuntimeError(
            "cannot close bounded-process supervisor capture descriptors: " +
            "; ".join(str(error) for error in errors)
        ) from errors[0]


def _cleanup_after_supervisor_collection_failure(
        supervisor: subprocess.Popen[bytes], control_read: int,
        root: tuple[int, int, int, int] | None,
        pidfds: RetainedPidfds) -> None:
    """Reap a failed collector's exact tree and close every capture endpoint."""
    errors: list[BaseException] = []
    capture_closed = False
    if root is not None:
        try:
            terminate_supervisor_tree_bounded(
                supervisor, root, pidfds=pidfds)
        except BaseException as error:
            errors.append(error)
    else:
        # The launch gate was never released.  Closing capture endpoints keeps
        # an error-reporting child from blocking on a full pipe while it exits.
        try:
            _close_supervisor_capture(supervisor, control_read)
            capture_closed = True
        except BaseException as error:
            errors.append(error)
        try:
            supervisor.wait(timeout=CHILD_REAP_TIMEOUT_SECONDS)
        except BaseException as error:
            errors.append(error)
    if not capture_closed:
        try:
            _close_supervisor_capture(supervisor, control_read)
        except BaseException as error:
            errors.append(error)
    if errors:
        raise RuntimeError(
            "bounded-process supervisor collection cleanup failed: " +
            "; ".join(f"{type(error).__name__}: {error}" for error in errors)
        ) from errors[0]


def _collect_supervisor(
        supervisor: subprocess.Popen[bytes], control_read: int,
        root: tuple[int, int, int, int], pidfds: RetainedPidfds,
        deadline: float, max_stdout: int,
        max_stderr: int) -> tuple[bytes, bytes, bytes]:
    require(supervisor.stdout is not None and supervisor.stderr is not None,
            "bounded-process supervisor capture pipes are absent")
    selector = selectors.DefaultSelector()
    stdout = bytearray()
    stderr = bytearray()
    control = bytearray()
    channels = {
        supervisor.stdout.fileno(): (supervisor.stdout, stdout, max_stdout),
        supervisor.stderr.fileno(): (supervisor.stderr, stderr, max_stderr),
        control_read: (control_read, control, MAX_SUPERVISOR_CONTROL_BYTES),
    }
    failure: str | None = None
    leader_exit_deadline: float | None = None
    try:
        for stream, unused_output, unused_limit in channels.values():
            del unused_output, unused_limit
            os.set_blocking(
                stream if isinstance(stream, int) else stream.fileno(),
                False)
            selector.register(stream, selectors.EVENT_READ)
        while selector.get_map() or supervisor.poll() is None:
            observed = supervisor.poll()
            now = time.monotonic()
            if observed is not None and leader_exit_deadline is None:
                leader_exit_deadline = min(
                    deadline, now + CHILD_REAP_TIMEOUT_SECONDS)
            effective_deadline = (
                min(deadline, leader_exit_deadline)
                if leader_exit_deadline is not None else deadline)
            remaining = effective_deadline - now
            if remaining <= 0:
                failure = (
                    "bounded-process supervisor exceeded its audited "
                    "cleanup bound")
                break
            events = selector.select(min(remaining, 0.05))
            for key, unused_mask in events:
                del unused_mask
                descriptor = key.fd
                stream, output, limit = channels[descriptor]
                try:
                    block = os.read(descriptor, 65536)
                except BlockingIOError:
                    continue
                if not block:
                    selector.unregister(stream)
                    continue
                output.extend(block)
                if len(output) > limit:
                    failure = (
                        "bounded-process supervisor exceeded an output "
                        "byte limit")
                    break
            if failure is not None:
                break

        if failure is not None:
            terminate_supervisor_tree_bounded(
                supervisor, root, pidfds=pidfds)
            raise RuntimeError(failure)
        require(supervisor.returncode is not None,
                "bounded-process supervisor was not reaped")
        return bytes(stdout), bytes(stderr), bytes(control)
    finally:
        selector.close()
        _close_supervisor_capture(supervisor, control_read)


def _run_supervised_command(
        command_value: Sequence[str], timeout: float, *,
        max_stdout: int, max_stderr: int,
        inherited_descriptors: Sequence[int],
        environment: Mapping[str, str]) -> subprocess.CompletedProcess[bytes]:
    require(isinstance(command_value, Sequence) and
            not isinstance(command_value, (str, bytes, bytearray)) and
            command_value and
            len(command_value) <= 512 and
            all(isinstance(item, str) and item and "\0" not in item
                for item in command_value),
            "benchmark command is invalid")
    require(isinstance(timeout, (int, float)) and
            not isinstance(timeout, bool) and
            timeout > 0 and math.isfinite(timeout) and timeout <= 3600.0,
            "benchmark timeout must be finite and positive")
    require(type(max_stdout) is int and
            0 <= max_stdout <= MAX_HELPER_STDOUT_BYTES and
            type(max_stderr) is int and
            0 <= max_stderr <= MAX_HELPER_STDERR_BYTES,
            "supervised command output limit is invalid")
    pass_fds = inherited_pass_fds(
        inherited_descriptors, "supervised command")
    for descriptor in pass_fds:
        try:
            os.fstat(descriptor)
        except OSError as error:
            raise RuntimeError(
                f"benchmark inherited descriptor {descriptor} is invalid: "
                f"{error}") from error
    require(isinstance(environment, Mapping) and
            1 <= len(environment) <= 128 and
            all(isinstance(name, str) and name and "\0" not in name and
                isinstance(value, str) and "\0" not in value
                for name, value in environment.items()),
            "supervised command environment is invalid")
    command_json = json.dumps(
        list(command_value), ensure_ascii=True, separators=(",", ":"))
    environment_json = json.dumps(
        dict(environment), ensure_ascii=True, separators=(",", ":"),
        sort_keys=True)
    inherited_json = json.dumps(list(pass_fds), separators=(",", ":"))
    control_read, control_write = os.pipe2(getattr(os, "O_CLOEXEC", 0))
    try:
        launch_read, launch_write = os.pipe2(
            getattr(os, "O_CLOEXEC", 0))
    except BaseException:
        # A process-wide descriptor shortage can make the second allocation
        # fail.  Roll back the already-created control pipe without obscuring
        # the allocation error.
        for descriptor in (control_read, control_write):
            try:
                os.close(descriptor)
            except OSError:
                pass
        raise
    supervisor: subprocess.Popen[bytes] | None = None
    pidfds = RetainedPidfds()
    root: tuple[int, int, int, int] | None = None
    try:
        supervisor = subprocess.Popen(
            [
                sys.executable, str(Path(__file__).resolve()),
                BOUNDED_SUPERVISOR_MODE, command_json, repr(float(timeout)),
                str(max_stdout), str(max_stderr), environment_json,
                inherited_json,
                str(control_write), str(launch_read),
            ],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, env=CHILD_ENV, start_new_session=True,
            pass_fds=tuple(sorted({
                *pass_fds, control_write, launch_read,
            })))
    except BaseException:
        for descriptor in (
                control_read, control_write, launch_read, launch_write):
            os.close(descriptor)
        raise
    os.close(control_write)
    os.close(launch_read)
    try:
        root = pidfds.retain(supervisor.pid)
        require(root is not None,
                "bounded-process supervisor disappeared during launch")
        _write_descriptor_all(launch_write, b"\x01")
    except BaseException:
        # Closing the unreleased gate makes the child reject launch before it
        # imports containment code or starts the benchmark.
        os.close(launch_write)
        try:
            try:
                _collect_supervisor(
                    supervisor, control_read,
                    (root if root is not None else
                     (supervisor.pid, -1, -1, -1)),
                    pidfds, time.monotonic() + CHILD_REAP_TIMEOUT_SECONDS,
                    max_stdout, max_stderr)
            except BaseException as collection_error:
                try:
                    _cleanup_after_supervisor_collection_failure(
                        supervisor, control_read, root, pidfds)
                except BaseException as cleanup_error:
                    raise RuntimeError(
                        "bounded-process supervisor launch, collection, and "
                        "cleanup all failed"
                    ) from cleanup_error
                raise RuntimeError(
                    "bounded-process supervisor launch and collection both "
                    "failed"
                ) from collection_error
        except BaseException:
            try:
                pidfds.close()
            except BaseException as close_error:
                raise RuntimeError(
                    "bounded-process supervisor launch failed and retained "
                    "pidfd cleanup also failed"
                ) from close_error
            raise
        else:
            pidfds.close()
        raise
    os.close(launch_write)
    try:
        try:
            require(root is not None,
                    "bounded-process supervisor root identity was lost")
            stdout, stderr, control_bytes = _collect_supervisor(
                supervisor, control_read, root, pidfds,
                time.monotonic() + timeout +
                2.0 * CHILD_REAP_TIMEOUT_SECONDS + 5.0,
                max_stdout, max_stderr)
        except BaseException:
            # _collect_supervisor normally owns timeout/output-limit cleanup,
            # but selector, pipe, or close failures can interrupt it at any
            # point.  Never return through such an exception with a live tree
            # (and therefore a leaked inherited campaign lock).
            try:
                _cleanup_after_supervisor_collection_failure(
                    supervisor, control_read, root, pidfds)
            except BaseException as cleanup_error:
                # This exception is raised while handling the collection
                # exception, so Python retains that original exception as its
                # context and the explicit cause retains the cleanup failure.
                raise RuntimeError(
                    "bounded-process supervisor collection failed and exact "
                    "tree cleanup also failed"
                ) from cleanup_error
            raise
    except BaseException:
        try:
            pidfds.close()
        except BaseException as close_error:
            # As above, the active operation/cleanup exception remains the
            # context of this explicit descriptor-close failure.
            raise RuntimeError(
                "bounded-process supervisor failed and retained pidfd cleanup "
                "also failed"
            ) from close_error
        raise
    else:
        pidfds.close()
    control = strict_json_bytes(
        control_bytes, "bounded-process supervisor control")
    require(isinstance(control, dict),
            "bounded-process supervisor control is not an object")
    require(supervisor.returncode == 0,
            "bounded-process supervisor failed: " +
            stderr.decode("utf-8", errors="replace"))
    require(len(stdout) <= max_stdout and len(stderr) <= max_stderr,
            "bounded-process supervisor exceeded an output byte limit")
    status = control.get("status")
    if status == "timeout":
        require(set(control) == {"status", "message"} and
                isinstance(control.get("message"), str),
                "bounded-process supervisor timeout record is malformed")
        raise subprocess.TimeoutExpired(
            list(command_value), timeout, output=stdout, stderr=stderr)
    if status == "error":
        require(set(control) == {"status", "message"} and
                isinstance(control.get("message"), str),
                "bounded-process supervisor error record is malformed")
        raise RuntimeError(
            "bounded-process supervisor rejected the command: " +
            control["message"])
    require(status == "ok" and set(control) == {"status", "returncode"} and
            type(control.get("returncode")) is int,
            "bounded-process supervisor success record is malformed")
    return subprocess.CompletedProcess(
        list(command_value), control["returncode"], stdout, stderr)


def run_helper_bounded(
        command_value: Sequence[str], *, timeout: float = 120.0,
        max_stdout: int = MAX_HELPER_STDOUT_BYTES,
        max_stderr: int = MAX_HELPER_STDERR_BYTES,
        inherited_descriptors: Sequence[int] = (),
        env: Mapping[str, str] = CHILD_ENV,
        ) -> subprocess.CompletedProcess[bytes]:
    return _run_supervised_command(
        command_value, timeout, max_stdout=max_stdout,
        max_stderr=max_stderr,
        inherited_descriptors=inherited_descriptors,
        environment=env)


def run_process_bounded(
        command_value: Sequence[str], timeout: float,
        snapshot_descriptor: int,
        lock_descriptor: int) -> subprocess.CompletedProcess[bytes]:
    require(type(snapshot_descriptor) is int and snapshot_descriptor >= 0 and
            type(lock_descriptor) is int and lock_descriptor >= 0,
            "benchmark inherited descriptor is invalid")
    return _run_supervised_command(
        command_value, timeout,
        max_stdout=MAX_BENCHMARK_STDOUT_BYTES,
        max_stderr=MAX_BENCHMARK_STDERR_BYTES,
        inherited_descriptors=(snapshot_descriptor, lock_descriptor),
        environment=CHILD_ENV)


def validate_invocation_record(
        record: Mapping[str, Any], label: str) -> None:
    require(isinstance(record, dict), f"{label} is not an object")
    returncode = record.get("returncode")
    require(type(returncode) is int, f"{label} return code is not an integer")
    common_keys = {
        "role", "command", "executable_snapshot_sha256", "returncode",
        "duration_ns", "stdout_sha256", "stdout", "stderr",
    }
    expected_keys = common_keys | ({"result"} if returncode == 0 else set())
    require(set(record) == expected_keys,
            f"{label} keys differ from its return-code contract")
    require(record.get("role") in {"main", "leopard2"} and
            isinstance(record.get("command"), list) and record["command"] and
            all(isinstance(item, str) and item for item in record["command"]),
            f"{label} role or command is malformed")
    snapshot = record.get("executable_snapshot_sha256")
    require(isinstance(snapshot, str) and
            re.fullmatch(r"[0-9a-f]{64}", snapshot) is not None,
            f"{label} executable snapshot digest is malformed")
    duration = record.get("duration_ns")
    require(type(duration) is int and duration > 0,
            f"{label} duration is not a positive integer")
    stdout = record.get("stdout")
    stderr = record.get("stderr")
    require(isinstance(stdout, str) and isinstance(stderr, str),
            f"{label} retained streams are not text")
    try:
        stdout_bytes = stdout.encode("utf-8", errors="strict")
        stderr_bytes = stderr.encode("utf-8", errors="strict")
    except UnicodeEncodeError as error:
        raise RuntimeError(f"{label} retained streams are not strict UTF-8") \
            from error
    require(len(stdout_bytes) <= MAX_BENCHMARK_STDOUT_BYTES and
            len(stderr_bytes) <= MAX_BENCHMARK_STDERR_BYTES,
            f"{label} retained stream exceeds its byte bound")
    require(record.get("stdout_sha256") ==
            hashlib.sha256(stdout_bytes).hexdigest(),
            f"{label} stdout SHA-256 differs from retained raw stdout")
    if returncode == 0:
        require(stderr == "", f"{label} successful invocation wrote stderr")
        parsed = strict_json_bytes(stdout_bytes, f"{label} raw stdout")
        require(canonical_equal(parsed, record.get("result")),
                f"{label} parsed result differs from retained raw stdout")


def run_one(role: str, command_value: Sequence[str], timeout: float,
            snapshot_descriptor: int,
            snapshot_identity: Mapping[str, Any],
            lock_descriptor: int) -> dict[str, Any]:
    execution_command = list(command_value)
    require(len(execution_command) >= 4 and
            execution_command[0:2] == ["/usr/bin/taskset", "-c"],
            "all-K execution command has an unexpected launcher shape")
    execution_command[3] = f"/proc/self/fd/{snapshot_descriptor}"
    started = time.monotonic_ns()
    completed = run_process_bounded(
        execution_command, timeout, snapshot_descriptor, lock_descriptor)
    finished = time.monotonic_ns()
    try:
        stdout = completed.stdout.decode("utf-8", errors="strict")
        stderr = completed.stderr.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise RuntimeError("benchmark streams are not strict UTF-8") from error
    record: dict[str, Any] = {
        "role": role, "command": list(command_value),
        "executable_snapshot_sha256": snapshot_identity["sha256"],
        "returncode": completed.returncode,
        "duration_ns": finished - started,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stdout": stdout, "stderr": stderr,
    }
    if completed.returncode == 0:
        record["result"] = strict_json_bytes(
            completed.stdout, "benchmark raw stdout")
    validate_invocation_record(record, "live benchmark invocation")
    return record


def metric(record: Mapping[str, Any], path: Sequence[str]) -> float:
    value: Any = record["result"]
    for key in path:
        value = value[key]
    label = "/".join(path)
    require(type(value) in (int, float),
            f"benchmark metric {label} is not an exact JSON number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise RuntimeError(
            f"benchmark metric {label} is outside the finite float range"
        ) from error
    require(result > 0 and math.isfinite(result),
            f"benchmark metric {label} is not finite and positive")
    return result


def exact_positive_finite(value: Any, label: str) -> float:
    require(type(value) in (int, float),
            f"{label} is not an exact JSON number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise RuntimeError(f"{label} is outside the finite float range") \
            from error
    require(result > 0 and math.isfinite(result),
            f"{label} is not finite and positive")
    return result


def geometric(values: Sequence[float]) -> float:
    require(bool(values) and
            all(value > 0 and math.isfinite(value) for value in values),
            "benchmark timing samples must be finite and positive")
    try:
        result = math.exp(
            statistics.fmean(math.log(value) for value in values))
    except OverflowError as error:
        raise RuntimeError("benchmark geometric mean overflowed") from error
    require(result > 0 and math.isfinite(result),
            "benchmark geometric mean is not finite and positive")
    return result


def finite_ratio(numerator: float, denominator: float, label: str) -> float:
    require(numerator > 0 and denominator > 0 and
            math.isfinite(numerator) and math.isfinite(denominator),
            f"{label} operands are not finite and positive")
    result = numerator / denominator
    require(result > 0 and math.isfinite(result),
            f"{label} is not finite and positive")
    return result


def validate_invocation_identities(
        invocations: Sequence[Mapping[str, Any]], main_commit: str,
        current_source: Mapping[str, Any],
        main_snapshot: Mapping[str, Any],
        current_snapshot: Mapping[str, Any]) -> None:
    require(len(invocations) == len(ORDER),
            "a valid all-K cell must contain exactly four invocations")
    require(tuple(record.get("role") for record in invocations) == ORDER,
            "all-K invocation role order mismatch")
    for index, record in enumerate(invocations):
        require(record.get("returncode") == 0,
                f"all-K invocation {index} did not exit successfully")
        result = record.get("result")
        require(isinstance(result, dict),
                f"all-K invocation {index} has no JSON result")
        build = result.get("build")
        require(isinstance(build, dict),
                f"all-K invocation {index} has no build identity")
        if record["role"] == "main":
            expected_snapshot = main_snapshot
            require(build.get("main_source_commit") == main_commit,
                    "exact-main benchmark embedded commit mismatch")
        else:
            expected_snapshot = current_snapshot
            require(result.get("schema") == "leopard2-benchmark-v3",
                    "Leopard2 benchmark decode-path schema mismatch")
            require(result.get("parameters", {}).get(
                "report_decode_path") is True,
                "Leopard2 benchmark did not record decode-path reporting")
        require(record.get("executable_snapshot_sha256") ==
                expected_snapshot.get("sha256"),
                f"all-K invocation {index} executable snapshot mismatch")
        validate_invocation_record(record, f"all-K invocation {index}")


def validate_source_attestation(
        record: Mapping[str, Any], current_source: Mapping[str, Any],
        current_snapshot: Mapping[str, Any]) -> None:
    require(record.get("returncode") == 0,
            "Leopard2 source-attestation preflight failed")
    require(record.get("executable_snapshot_sha256") ==
            current_snapshot.get("sha256"),
            "Leopard2 source-attestation snapshot mismatch")
    result = record.get("result")
    require(isinstance(result, dict) and
            result.get("schema") == "leopard2-benchmark-v5",
            "Leopard2 benchmark source-attestation schema mismatch")
    require(result.get("parameters", {}).get("attest_source") is True,
            "Leopard2 benchmark did not record source attestation")
    build = result.get("build")
    require(isinstance(build, dict),
            "Leopard2 source attestation has no build identity")
    require(build.get("source_commit") == current_source.get("head"),
            "Leopard2 benchmark embedded commit mismatch")
    require(build.get("source_tree") == current_source.get("tree"),
            "Leopard2 benchmark embedded tree mismatch")
    require(build.get("source_tracked_dirty") is False,
            "Leopard2 benchmark was built from a tracked-dirty tree")
    resolved = result.get("resolved", {})
    require((resolved.get("profile"), resolved.get("field"),
             resolved.get("backend")) ==
            ("legacy_high_v1", "gf8", "avx2"),
            "Leopard2 source attestation resolved the wrong codec identity")
    require(result.get("correctness", {}).get(
        "leopard2_round_trip") is True,
        "Leopard2 source-attestation round trip failed")
    validate_invocation_record(record, "Leopard2 source-attestation invocation")


def source_attestation_identity(
        record: Mapping[str, Any], current_source: Mapping[str, Any],
        current_snapshot: Mapping[str, Any]) -> dict[str, Any]:
    """Extract only deterministic evidence from a complete raw preflight."""
    validate_source_attestation(record, current_source, current_snapshot)
    result = record["result"]
    require(isinstance(result, dict) and "metrics" in result,
            "Leopard2 source attestation has no volatile metrics section")
    stable_result = {
        key: copy.deepcopy(value)
        for key, value in result.items() if key != "metrics"
    }
    identity = {
        "schema": "leopard2-all-k-source-attestation-identity/v1",
        "role": record["role"],
        "command": copy.deepcopy(record["command"]),
        "executable_snapshot_sha256":
            record["executable_snapshot_sha256"],
        "returncode": record["returncode"],
        "result_without_metrics": stable_result,
    }
    # This is also a final guard against accidentally promoting a volatile
    # non-finite value into the persistent run contract.
    canonical_digest(identity)
    return identity


def validate_correctness(cell: Cell,
                         invocations: Sequence[Mapping[str, Any]]) -> None:
    expected_parameters = {
        "K": cell.k, "R": cell.r, "shard_bytes": cell.shard_bytes,
        "loss_count": cell.losses, "batch": 1, "reuse": cell.reuse,
        "iterations": cell.iterations, "warmup": cell.warmup,
        "thread_count": 1, "seed": cell.seed,
    }
    fingerprints: list[Mapping[str, Any]] = []
    expected_parent, expected_padded = parent_for(cell.k, cell.r)
    expected_field = "gf8" if cell.family == "gf8-all-k" else "gf16"
    for index, record in enumerate(invocations):
        result = record["result"]
        parameters = result.get("parameters")
        require(isinstance(parameters, dict),
                f"all-K invocation {index} has no parameters")
        for name, value in expected_parameters.items():
            require(type(parameters.get(name)) is int and
                    parameters.get(name) == value,
                    f"all-K invocation {index} parameter {name} mismatch")
        missing = parameters.get("missing_original_indices")
        require(isinstance(missing, list) and
                all(type(value) is int for value in missing) and
                canonical_equal(missing, expected_losses(cell)),
                f"all-K invocation {index} missing-original indices mismatch")
        resolved = result.get("resolved")
        require(isinstance(resolved, dict) and
                resolved.get("profile") == "legacy_high_v1" and
                resolved.get("field") == expected_field and
                type(resolved.get("thread_count")) is int and
                resolved.get("thread_count") == 1 and
                type(resolved.get("parent_count")) is int and
                resolved.get("parent_count") == expected_parent and
                type(resolved.get("padded_side")) is int and
                resolved.get("padded_side") == expected_padded,
                f"all-K invocation {index} resolved geometry mismatch")
        if record["role"] == "leopard2":
            expected_options: dict[str, Any] = {
                "requested_profile": "legacy_high_v1",
                "requested_field": expected_field,
                "requested_backend": "avx2",
                "force_generic_decode": False,
                "force_specialized_decode": False,
                "force_tiled_decode": False,
                "force_materialized_decode": False,
                "skip_legacy": "--skip-legacy" in record["command"],
                "retain_samples": True,
                "report_decode_path": True,
            }
            require(record["command"].count("--skip-legacy") <= 1,
                    f"all-K invocation {index} repeats --skip-legacy")
            for name, value in expected_options.items():
                require(type(parameters.get(name)) is type(value) and
                        parameters.get(name) == value,
                        f"all-K invocation {index} option {name} mismatch")
        correctness = result.get("correctness")
        require(isinstance(correctness, dict),
                f"all-K invocation {index} has no correctness record")
        if record["role"] == "main":
            require(correctness.get("round_trip") is True,
                    "exact-main benchmark round trip failed")
        else:
            require(correctness.get("leopard2_round_trip") is True,
                    "Leopard2 benchmark round trip failed")
        fingerprint = result.get("workload_digests")
        expected_digest_keys = {
            "algorithm", "original_data", "transmitted_parity",
            "recovered_originals",
        }
        require(isinstance(fingerprint, dict) and
                set(fingerprint) == expected_digest_keys,
                f"all-K invocation {index} workload digest structure mismatch")
        require(fingerprint["algorithm"] == "fnv1a64",
                f"all-K invocation {index} workload digest algorithm mismatch")
        for name in ("original_data", "transmitted_parity",
                     "recovered_originals"):
            value = fingerprint[name]
            require(isinstance(value, str) and len(value) == 16 and
                    all(character in "0123456789abcdef"
                        for character in value),
                    f"all-K invocation {index} workload digest {name} "
                    "is not lowercase FNV-1a hex")
        # original_data covers every systematic shard, whereas
        # recovered_originals covers only the sorted missing-original set.
        # They are therefore expected to differ for a partial-loss workload.
        # The independently built implementations must still agree on the
        # complete fingerprint below, and each benchmark has already proved
        # its own byte-for-byte round trip above.  When every original is
        # missing, both traversals cover the same shards in the same order and
        # equality is an additional useful consistency check.
        if cell.losses == cell.k:
            require(fingerprint["recovered_originals"] ==
                    fingerprint["original_data"],
                    f"all-K invocation {index} full-loss recovered-original "
                    "digest does not match its input")
        fingerprints.append(fingerprint)
    require(all(canonical_equal(fingerprint, fingerprints[0])
                for fingerprint in fingerprints),
            "Leopard1 and Leopard2 workload digests differ")


def gap_tags(cell: Cell, paths: Mapping[str, Any], encode_speedup: float,
             decode_first_speedup: float, decode_reuse_speedup: float) -> list[str]:
    tags: list[str] = []
    if encode_speedup < 1.05:
        if paths["encode_path"] == "direct-xor-single-parity":
            tags.append("encode:R1-xor/API-overhead")
        else:
            tags.append("encode:legacy-high-transform-not-faster")
        if cell.k % paths["padded_side"] != 0:
            tags.append("encode:partial-final-message-block")
    if decode_first_speedup < 1.05 or decode_reuse_speedup < 1.05:
        decode_path = paths["decode_path"]
        if decode_path == "direct":
            tags.append("decode:direct-path-overhead-or-kernel")
        elif decode_path == "generic":
            tags.append("decode:generic-fallback-crossover")
        else:
            tags.append("decode:specialized-high-crossover")
        if paths["decode_workspace"] == "materialized-N":
            tags.append("decode:materialized-N-workspace")
        if paths["locator_setup"].endswith("walsh"):
            tags.append("decode:active-parent-Walsh-setup")
    if paths["parent_inflation"] >= 1.5:
        tags.append("common:dyadic-parent-inflation>=1.5x")
    if cell.shard_bytes <= 4096:
        tags.append("common:small-shard-fixed-cost")
    return sorted(set(tags))


def analyze_cell(cell: Cell, invocations: Sequence[Mapping[str, Any]], cpu: int,
                 contract_digest: str, main_commit: str,
                 current_source: Mapping[str, Any],
                 main_snapshot: Mapping[str, Any],
                 current_snapshot: Mapping[str, Any]) -> dict[str, Any]:
    failures = [record for record in invocations if record["returncode"] != 0]
    result: dict[str, Any] = {
        "schema": "leopard2-all-k-gap-cell/v2",
        "run_contract_sha256": contract_digest,
        "cell": dataclasses.asdict(cell),
        "cpu": cpu, "order": list(ORDER), "invocations": list(invocations),
        "valid": not failures, "diagnostic_not_promotion_evidence": True,
    }
    if failures:
        result["failures"] = failures
        return result
    validate_invocation_identities(
        invocations, main_commit, current_source,
        main_snapshot, current_snapshot)
    validate_correctness(cell, invocations)
    main = [record for record in invocations if record["role"] == "main"]
    current = [record for record in invocations if record["role"] == "leopard2"]
    require(len(main) == len(current) == 2,
            "all-K ABBA role cardinality mismatch")
    paths = classify_paths(cell, current[0]["result"])
    require(all(classify_paths(cell, record["result"]) == paths
                for record in current),
            "Leopard2 selected-path identity changed within an all-K cell")
    main_encode = geometric([metric(record, ("metrics", "encode_execution",
        "median_us_per_batch_call")) for record in main])
    current_encode = geometric([metric(record, ("metrics", "encode_execution",
        "median_us_per_batch_call")) for record in current])
    main_decode = geometric([metric(record, ("metrics", "decode_including_setup",
        "median_us_per_batch_call")) for record in main])
    current_plan = geometric([metric(record, ("metrics", "decode_plan_setup",
        "median_us")) for record in current])
    current_decode = geometric([metric(record, ("metrics", "decode_execution",
        "median_us_per_batch_call")) for record in current])
    current_codec = geometric([metric(record, ("metrics", "codec_setup",
        "median_us")) for record in current])
    first = current_plan + current_decode
    amortized = current_decode + current_plan / cell.reuse
    require(first > 0 and amortized > 0 and
            math.isfinite(first) and math.isfinite(amortized),
            "derived Leopard2 decode timing is not finite and positive")
    encode_speedup = finite_ratio(
        main_encode, current_encode, "encode speedup")
    first_speedup = finite_ratio(
        main_decode, first, "first-use decode speedup")
    reuse_speedup = finite_ratio(
        main_decode, amortized, "reused-plan decode speedup")
    optimistic_speedup = finite_ratio(
        main_decode, current_decode, "execution-only decode speedup")
    paths["gap_tags"] = gap_tags(
        cell, paths, encode_speedup, first_speedup, reuse_speedup)
    result.update({
        "selected": paths,
        "timing_us": {
            "main_encode_execution": main_encode,
            "leopard2_codec_setup": current_codec,
            "leopard2_encode_execution": current_encode,
            "main_decode_including_setup": main_decode,
            "leopard2_decode_plan_setup": current_plan,
            "leopard2_decode_execution": current_decode,
            "leopard2_decode_first_use": first,
            "leopard2_decode_amortized_at_reuse": amortized,
        },
        "speedup_main_over_leopard2": {
            "encode": encode_speedup,
            "decode_first_use": first_speedup,
            "decode_at_reuse": reuse_speedup,
            "decode_execution_only_optimistic": optimistic_speedup,
        },
        "significantly_beats_main_1_05": {
            "encode": encode_speedup >= 1.05,
            "decode_first_use": first_speedup >= 1.05,
            "decode_at_reuse": reuse_speedup >= 1.05,
        },
        "memory": current[0]["result"]["memory"],
    })
    legacy_available: list[bool] = []
    for index, record in enumerate(current):
        legacy = record["result"].get("legacy")
        require(isinstance(legacy, dict),
                f"Leopard2 invocation {index} has no legacy result")
        available = legacy.get("available")
        require(type(available) is bool,
                f"Leopard2 invocation {index} legacy availability "
                "is not an exact boolean")
        expected_available = "--skip-legacy" not in record["command"]
        require(available is expected_available,
                f"Leopard2 invocation {index} legacy availability "
                "contradicts its command")
        legacy_available.append(available)
    require(len(set(legacy_available)) == 1,
            "Leopard2 legacy availability changed within an all-K cell")
    if legacy_available[0]:
        current_legacy_encode = geometric([metric(
            record, ("legacy", "encode_execution", "median_us_per_batch_call"))
            for record in current])
        current_legacy_decode = geometric([metric(
            record, ("legacy", "decode_including_setup",
                     "median_us_per_batch_call"))
            for record in current])
        result["timing_us"].update({
            "current_tree_legacy_encode_execution": current_legacy_encode,
            "current_tree_legacy_decode_including_setup": current_legacy_decode,
        })
        result["attribution_speedup"] = {
            "exact_main_over_current_tree_legacy": {
                "encode": finite_ratio(
                    main_encode, current_legacy_encode,
                    "exact-main/current-legacy encode ratio"),
                "decode_including_setup": finite_ratio(
                    main_decode, current_legacy_decode,
                    "exact-main/current-legacy decode ratio"),
            },
            "current_tree_legacy_over_leopard2": {
                "encode": finite_ratio(
                    current_legacy_encode, current_encode,
                    "current-legacy/Leopard2 encode ratio"),
                "decode_first_use": finite_ratio(
                    current_legacy_decode, first,
                    "current-legacy/Leopard2 first-use decode ratio"),
                "decode_at_reuse": finite_ratio(
                    current_legacy_decode, amortized,
                    "current-legacy/Leopard2 reused-plan decode ratio"),
                "decode_execution_only_optimistic":
                    finite_ratio(
                        current_legacy_decode, current_decode,
                        "current-legacy/Leopard2 execution-only decode ratio"),
            },
        }
    return result


def validate_cell_document(
        document: Mapping[str, Any], cell: Cell, cpu: int,
        contract_digest: str, main: Path, current: Path,
        with_current_legacy: bool, main_commit: str,
        current_source: Mapping[str, Any],
        main_snapshot: Mapping[str, Any],
        current_snapshot: Mapping[str, Any]) -> None:
    require(document.get("schema") == "leopard2-all-k-gap-cell/v2",
            f"stored cell {cell.identifier} schema mismatch")
    require(document.get("run_contract_sha256") == contract_digest,
            f"stored cell {cell.identifier} contract mismatch")
    require(canonical_equal(document.get("cell"), dataclasses.asdict(cell)),
            f"stored cell {cell.identifier} parameters mismatch")
    require(type(document.get("cpu")) is int and document.get("cpu") == cpu,
            f"stored cell {cell.identifier} CPU mismatch")
    require(canonical_equal(document.get("order"), list(ORDER)),
            f"stored cell {cell.identifier} order mismatch")
    invocations = document.get("invocations")
    require(isinstance(invocations, list) and 1 <= len(invocations) <= len(ORDER),
            f"stored cell {cell.identifier} invocation count is invalid")
    for slot, record in enumerate(invocations):
        require(isinstance(record, dict),
                f"stored cell {cell.identifier} invocation is not an object")
        role = ORDER[slot]
        executable = main if role == "main" else current
        require(record.get("role") == role,
                f"stored cell {cell.identifier} invocation role mismatch")
        require(canonical_equal(record.get("command"), command(
            role, executable, cell, cpu, with_current_legacy)),
            f"stored cell {cell.identifier} command mismatch")
        validate_invocation_record(
            record, f"stored cell {cell.identifier} invocation {slot}")
    recomputed = analyze_cell(
        cell, invocations, cpu, contract_digest, main_commit, current_source,
        main_snapshot, current_snapshot)
    require(set(document) == set(recomputed),
            f"stored cell {cell.identifier} keys differ from exact "
            "recomputation")
    require(document == recomputed and
            canonical_digest(document) == canonical_digest(recomputed),
            f"stored cell {cell.identifier} differs from exact recomputation "
            "of retained invocations")


def run_cell(cell: Cell, cpu: int, main: Path,
             current: Path, output: Path, timeout: float,
             with_current_legacy: bool, contract_digest: str,
             main_commit: str,
             current_source: Mapping[str, Any],
             lock_descriptor: int,
             main_snapshot_descriptor: int,
             main_snapshot: Mapping[str, Any],
             current_snapshot_descriptor: int,
             current_snapshot: Mapping[str, Any],
             evidence_directory: AtomicEvidenceDirectory | None = None,
             ) -> dict[str, Any]:
    path = output / "cells" / (cell.identifier + ".json")
    stored: Any | None = None
    if evidence_directory is not None:
        relative = path.relative_to(evidence_directory.path).as_posix()
        payload = evidence_directory.read_optional(
            relative, MAX_PERSISTED_CELL_BYTES,
            f"stored cell {cell.identifier}")
        if payload is not None:
            stored = strict_json_bytes(
                payload, f"stored cell {cell.identifier}")
    elif path.is_file():
        stored = load_strict_json(
            path, f"stored cell {cell.identifier}",
            MAX_PERSISTED_CELL_BYTES)
    if stored is not None:
        require(isinstance(stored, dict),
                f"stored cell {cell.identifier} is not a JSON object")
        validate_cell_document(
            stored, cell, cpu, contract_digest, main, current,
            with_current_legacy, main_commit, current_source,
            main_snapshot, current_snapshot)
        if stored["valid"] is True:
            return stored
    invocations = []
    for role in ORDER:
        executable = main if role == "main" else current
        snapshot_descriptor = (main_snapshot_descriptor if role == "main" else
                               current_snapshot_descriptor)
        snapshot_identity = (main_snapshot if role == "main" else
                             current_snapshot)
        logical_command = command(
            role, executable, cell, cpu, with_current_legacy)
        try:
            invocations.append(run_one(
                role, logical_command, timeout,
                snapshot_descriptor, snapshot_identity, lock_descriptor))
        except subprocess.TimeoutExpired as error:
            empty_stdout = ""
            invocations.append({
                "role": role, "command": logical_command,
                "executable_snapshot_sha256": snapshot_identity["sha256"],
                "returncode": -999,
                "duration_ns": max(1, int(timeout * 1e9)),
                "stdout_sha256": hashlib.sha256(b"").hexdigest(),
                "stdout": empty_stdout, "stderr": "timeout",
            })
            break
    result = analyze_cell(
        cell, invocations, cpu, contract_digest, main_commit, current_source,
        main_snapshot, current_snapshot)
    write_json_atomic(
        path, result, evidence_directory=evidence_directory)
    return result


def summarize(results: Sequence[Mapping[str, Any]], metadata: Mapping[str, Any]) \
        -> dict[str, Any]:
    require(isinstance(metadata, Mapping),
            "all-K summary metadata is not an object")
    run_contract = metadata.get("run_contract")
    require(isinstance(run_contract, Mapping) and
            type(run_contract.get("with_current_legacy")) is bool,
            "all-K summary metadata has no exact legacy-attribution contract")
    with_current_legacy = run_contract["with_current_legacy"]
    valid: list[Mapping[str, Any]] = []
    failed: list[Mapping[str, Any]] = []
    metrics = ("encode", "decode_first_use", "decode_at_reuse")
    attribution_presence: list[bool] = []
    for index, result in enumerate(results):
        require(isinstance(result, Mapping),
                f"all-K summary cell {index} is not an object")
        is_valid = result.get("valid")
        require(type(is_valid) is bool,
                f"all-K summary cell {index} validity is not an exact boolean")
        require(isinstance(result.get("cell"), Mapping),
                f"all-K summary cell {index} has no cell identity")
        if not is_valid:
            failed.append(result)
            continue
        valid.append(result)
        speedups = result.get("speedup_main_over_leopard2")
        require(isinstance(speedups, Mapping),
                f"all-K summary cell {index} has no speedup map")
        for name in metrics:
            exact_positive_finite(
                speedups.get(name),
                f"all-K summary cell {index} speedup {name}")
        selected = result.get("selected")
        require(isinstance(selected, Mapping),
                f"all-K summary cell {index} has no selected-path map")
        gap_tag_values = selected.get("gap_tags")
        require(isinstance(gap_tag_values, list) and
                all(isinstance(tag, str) for tag in gap_tag_values),
                f"all-K summary cell {index} gap tags are malformed")
        require(isinstance(result.get("timing_us"), Mapping),
                f"all-K summary cell {index} has no timing map")
        has_attribution = "attribution_speedup" in result
        attribution_presence.append(has_attribution)
        if has_attribution:
            attribution = result["attribution_speedup"]
            require(isinstance(attribution, Mapping),
                    f"all-K summary cell {index} attribution is not an object")
            expected_attribution = {
                "exact_main_over_current_tree_legacy": (
                    "encode", "decode_including_setup"),
                "current_tree_legacy_over_leopard2": (
                    "encode", "decode_first_use", "decode_at_reuse",
                    "decode_execution_only_optimistic"),
            }
            for leg, names in expected_attribution.items():
                values = attribution.get(leg)
                require(isinstance(values, Mapping),
                        f"all-K summary cell {index} attribution leg "
                        f"{leg} is missing")
                for name in names:
                    exact_positive_finite(
                        values.get(name),
                        f"all-K summary cell {index} attribution "
                        f"{leg}/{name}")
    if attribution_presence:
        require(all(value == attribution_presence[0]
                    for value in attribution_presence),
                "all-K valid cells have inconsistent legacy attribution")
        require(attribution_presence[0] is with_current_legacy,
                "all-K attribution presence contradicts the run contract")
    summary: dict[str, Any] = {
        "schema": "leopard2-all-k-gap-summary/v2", "metadata": dict(metadata),
        "cell_count": len(results), "valid_cell_count": len(valid),
        "failed_cell_count": len(failed),
        "diagnostic_not_promotion_evidence": True,
        "threshold": "Leopard2 significant when main_time/Leopard2_time >= 1.05",
        "metrics": {}, "attribution_metrics": {}, "gap_tags": {},
        "worst_cells": {},
    }
    for name in metrics:
        values = [result["speedup_main_over_leopard2"][name] for result in valid]
        slow = [value for value in values if value < 1.05]
        summary["metrics"][name] = {
            "count": len(values), "gap_count": len(slow),
            "gap_fraction": len(slow) / len(values) if values else None,
            "median_speedup": statistics.median(values) if values else None,
            "p10_speedup": sorted(values)[max(0, len(values) // 10 - 1)]
            if values else None,
        }
        ranked = sorted(valid, key=lambda result:
                        result["speedup_main_over_leopard2"][name])[:30]
        summary["worst_cells"][name] = [{
            "cell": result["cell"], "selected": result["selected"],
            "speedup": result["speedup_main_over_leopard2"][name],
            "timing_us": result["timing_us"],
        } for result in ranked]
    attribution_paths = {
        "exact_main_over_current_tree_legacy_encode":
            ("exact_main_over_current_tree_legacy", "encode"),
        "exact_main_over_current_tree_legacy_decode":
            ("exact_main_over_current_tree_legacy", "decode_including_setup"),
        "current_tree_legacy_over_leopard2_encode":
            ("current_tree_legacy_over_leopard2", "encode"),
        "current_tree_legacy_over_leopard2_decode_first_use":
            ("current_tree_legacy_over_leopard2", "decode_first_use"),
        "current_tree_legacy_over_leopard2_decode_at_reuse":
            ("current_tree_legacy_over_leopard2", "decode_at_reuse"),
    }
    if valid and attribution_presence[0]:
        for name, (leg, metric_name) in attribution_paths.items():
            values = [result["attribution_speedup"][leg][metric_name]
                      for result in valid]
            ordered = sorted(values)
            summary["attribution_metrics"][name] = {
                "count": len(values),
                "median_ratio": statistics.median(values),
                "p10_ratio": ordered[max(0, len(ordered) // 10 - 1)],
                "p90_ratio": ordered[min(len(ordered) - 1,
                                         9 * len(ordered) // 10)],
                "count_below_1_0": sum(value < 1.0 for value in values),
                "count_below_1_05": sum(value < 1.05 for value in values),
            }
    tags: dict[str, list[Mapping[str, Any]]] = {}
    for result in valid:
        for tag in result["selected"]["gap_tags"]:
            tags.setdefault(tag, []).append(result)
    summary["gap_tags"] = {
        tag: {
            "cell_count": len(items),
            "median_encode_speedup": statistics.median(
                item["speedup_main_over_leopard2"]["encode"] for item in items),
            "median_decode_first_speedup": statistics.median(
                item["speedup_main_over_leopard2"]["decode_first_use"]
                for item in items),
            "median_decode_reuse_speedup": statistics.median(
                item["speedup_main_over_leopard2"]["decode_at_reuse"]
                for item in items),
        } for tag, items in sorted(tags.items(), key=lambda pair: -len(pair[1]))
    }
    summary["failed_cells"] = [result["cell"] for result in failed]
    return summary


def run_serial_cpu_lanes(
        items: Sequence[Any], cpus: Sequence[int],
        execute: Callable[[Any, int], Mapping[str, Any]]) \
        -> list[Mapping[str, Any]]:
    """Run one serial queue per pinned CPU, with queues concurrent."""
    require(bool(cpus) and len(set(cpus)) == len(cpus),
            "serial CPU lanes require distinct CPUs")
    lanes: list[list[tuple[int, Any]]] = [[] for _ in cpus]
    for index, item in enumerate(items):
        lanes[index % len(cpus)].append((index, item))
    results: list[Mapping[str, Any] | None] = [None] * len(items)
    progress_lock = threading.Lock()
    completed = 0

    def run_lane(lane_index: int) -> None:
        nonlocal completed
        cpu = cpus[lane_index]
        for index, item in lanes[lane_index]:
            results[index] = execute(item, cpu)
            with progress_lock:
                completed += 1
                if completed % 50 == 0 or completed == len(items):
                    print(f"{completed}/{len(items)} cells", flush=True)

    with concurrent.futures.ThreadPoolExecutor(
            max_workers=min(len(cpus), len(items) or 1)) as executor:
        futures = [
            executor.submit(run_lane, lane_index)
            for lane_index, lane in enumerate(lanes) if lane
        ]
        for future in concurrent.futures.as_completed(futures):
            future.result()
    require(all(result is not None for result in results),
            "serial CPU lane scheduler omitted a result")
    return [result for result in results if result is not None]


def _run_with_snapshot_owner_held(
        options: argparse.Namespace, lock_guard: CanonicalBenchmarkLock,
        snapshots: ExecutableSnapshotOwner,
        evidence_directory: AtomicEvidenceDirectory) -> int:
    lock_identity = lock_guard.validate_current()
    lock_descriptor = lock_guard.descriptor
    require(lock_descriptor is not None,
            "canonical benchmark lock descriptor is unavailable")
    inherited_descriptors = (lock_descriptor,)
    main_commit = require_hex(options.main_commit, "exact-main source commit")
    expected_main_sha256 = require_sha256(
        options.main_sha256, "exact-main executable SHA-256")
    current_commit = require_hex(options.current_commit,
                                 "current source commit")
    current_source_initial = git_identity(
        options.current_source_root, current_commit,
        inherited_descriptors=inherited_descriptors)
    current_build_initial = candidate_build_provenance(
        options.current.resolve(strict=True).parent,
        options.current_source_root, options.current,
        "bench_leopard2_allk",
        inherited_descriptors=inherited_descriptors)
    current_reproducible_build = verify_reproducible_candidate_build(
        current_build_initial,
        jobs=options.reproducible_build_jobs,
        inherited_descriptors=inherited_descriptors)
    (main_executable_initial, main_snapshot_descriptor,
     main_snapshot_initial) = snapshots.capture(
        options.main, "exact-main benchmark")
    (current_executable_initial, current_snapshot_descriptor,
     current_snapshot_initial) = snapshots.capture_reproduced(
        options.current, current_build_initial, current_reproducible_build)
    lock_guard.validate_current()
    require(main_executable_initial["sha256"] == expected_main_sha256,
            "exact-main executable SHA-256 does not match the frozen "
            "pure-AVX2 baseline")
    main_exe = Path(main_executable_initial["path"])
    current_exe = Path(current_executable_initial["path"])
    main_isa = inspect_pure_avx2(
        main_snapshot_descriptor, "exact-main benchmark",
        inherited_descriptors=inherited_descriptors)
    current_isa = inspect_pure_avx2(
        current_snapshot_descriptor, "Leopard2 benchmark",
        inherited_descriptors=inherited_descriptors)
    output = evidence_directory.path
    cpus = sorted(os.sched_getaffinity(0))
    require(options.workers > 0, "workers must be positive")
    require(options.timeout > 0 and math.isfinite(options.timeout),
            "timeout must be finite and positive")
    workers = min(options.workers, len(cpus), 128)
    require(workers > 0, "the process has no allowed CPUs")
    cpus = cpus[:workers]
    current_source_attestation = run_one(
        "leopard2",
        source_attestation_command(current_exe, cpus[0]),
        options.timeout, current_snapshot_descriptor,
        current_snapshot_initial, lock_descriptor)
    validate_source_attestation(
        current_source_attestation, current_source_initial,
        current_snapshot_initial)
    current_source_attestation_identity = source_attestation_identity(
        current_source_attestation, current_source_initial,
        current_snapshot_initial)
    cells = make_cells(gf8_only=options.gf8_only)
    contract = {
        "schema": RUN_CONTRACT_SCHEMA,
        "main_commit": main_commit, "current_commit": current_commit,
        "expected_main_sha256": expected_main_sha256,
        "current_source_initial": current_source_initial,
        "current_build_initial": current_build_initial,
        "current_reproducible_build": current_reproducible_build,
        "main_executable_initial": main_executable_initial,
        "current_executable_initial": current_executable_initial,
        "main_executable_snapshot": main_snapshot_initial,
        "current_executable_snapshot": current_snapshot_initial,
        "current_source_attestation_identity":
            current_source_attestation_identity,
        "main_isa": main_isa, "current_isa": current_isa,
        "benchmark_lock": lock_identity,
        "allowed_cpus": sorted(os.sched_getaffinity(0)), "used_cpus": cpus,
        "workers": workers, "order": list(ORDER),
        "timeout_seconds": options.timeout,
        "with_current_legacy": options.with_current_legacy,
        "matrix": describe_cell_matrix(cells),
        "measurement_note": "all CPUs saturated; diagnostic crossover map, not isolated promotion evidence",
    }
    validate_run_contract_evidence(contract)
    contract_digest = canonical_digest(contract)
    manifest_path = output / "manifest.json"
    manifest: dict[str, Any] = {
        "schema": MANIFEST_SCHEMA,
        "run_contract": contract,
        "run_contract_sha256": contract_digest,
        "cells": [dataclasses.asdict(cell) for cell in cells],
        "current_source_attestation_preflights": [
            current_source_attestation
        ],
        "completion": None,
    }
    existing_manifest: Any | None = None
    manifest_payload = evidence_directory.read_optional(
        "manifest.json", MAX_PERSISTED_MANIFEST_BYTES,
        "existing all-K manifest")
    if manifest_payload is not None:
        existing_manifest = strict_json_bytes(
            manifest_payload, "existing all-K manifest")
    if existing_manifest is not None:
        require(isinstance(existing_manifest, dict),
                "existing all-K manifest is not a JSON object")
        validate_manifest(
            existing_manifest, contract, contract_digest, cells,
            current_source_initial, current_snapshot_initial,
            current_source_attestation_identity)
        manifest = existing_manifest
        require(len(manifest["current_source_attestation_preflights"]) <
                MAX_SOURCE_ATTESTATION_PREFLIGHTS,
                "all-K manifest source-attestation preflight limit reached")
        manifest["current_source_attestation_preflights"].append(
            current_source_attestation)
        validate_manifest(
            manifest, contract, contract_digest, cells,
            current_source_initial, current_snapshot_initial,
            current_source_attestation_identity)
        write_json_atomic(
            manifest_path, manifest,
            evidence_directory=evidence_directory)
    else:
        write_json_atomic(
            manifest_path, manifest,
            evidence_directory=evidence_directory)
    lock_guard.validate_current()
    results = run_serial_cpu_lanes(
        cells, cpus,
        lambda cell, cpu: run_cell(
            cell, cpu, main_exe, current_exe, output, options.timeout,
            options.with_current_legacy, contract_digest,
            main_commit, current_source_initial, lock_descriptor,
            main_snapshot_descriptor, main_snapshot_initial,
            current_snapshot_descriptor, current_snapshot_initial,
            evidence_directory=evidence_directory))
    lock_guard.validate_current()
    current_source_final = git_identity(
        options.current_source_root, current_commit,
        inherited_descriptors=inherited_descriptors)
    current_build_final = candidate_build_provenance(
        options.current.resolve(strict=True).parent,
        options.current_source_root, options.current,
        "bench_leopard2_allk",
        inherited_descriptors=inherited_descriptors)
    main_executable_final = file_identity(main_exe, "exact-main benchmark")
    current_executable_final = file_identity(
        current_exe, "Leopard2 benchmark")
    main_snapshot_final = sealed_snapshot_identity(
        main_snapshot_descriptor, "exact-main benchmark")
    current_snapshot_final = sealed_snapshot_identity(
        current_snapshot_descriptor, "Leopard2 benchmark")
    require(canonical_equal(current_source_final, current_source_initial),
            "current source identity changed during the all-K run")
    require(canonical_equal(current_build_final, current_build_initial),
            "current source/object/archive/link closure changed during the "
            "all-K run")
    require(canonical_equal(
        main_executable_final, main_executable_initial),
            "exact-main executable identity changed during the all-K run")
    require(canonical_equal(
        current_executable_final, current_executable_initial),
            "Leopard2 executable identity changed during the all-K run")
    require(canonical_equal(main_snapshot_final, main_snapshot_initial),
            "exact-main sealed executable snapshot changed")
    require(canonical_equal(current_snapshot_final, current_snapshot_initial),
            "Leopard2 sealed executable snapshot changed")
    completion = {
        "current_source_final": current_source_final,
        "current_build_final": current_build_final,
        "main_executable_final": main_executable_final,
        "current_executable_final": current_executable_final,
        "main_executable_snapshot_final": main_snapshot_final,
        "current_executable_snapshot_final": current_snapshot_final,
    }
    if manifest.get("completion") is not None:
        require(canonical_equal(manifest["completion"], completion),
                "existing all-K completion identity mismatch")
    manifest["completion"] = completion
    validate_manifest(
        manifest, contract, contract_digest, cells,
        current_source_initial, current_snapshot_initial,
        current_source_attestation_identity)
    write_json_atomic(
        manifest_path, manifest,
        evidence_directory=evidence_directory)
    cells_text = "".join(json.dumps(
        result, sort_keys=True, allow_nan=False) + "\n"
                         for result in results)
    cells_path = output / "cells.jsonl"
    _write_bytes_atomic(
        cells_path, cells_text.encode("utf-8", errors="strict"),
        evidence_directory=evidence_directory)
    gaps = [result for result in results if result.get("valid") is True and
            not all(result["significantly_beats_main_1_05"].values())]
    gap_text = "".join(json.dumps(
        result, sort_keys=True, allow_nan=False) + "\n"
                       for result in gaps)
    gap_path = output / "gap_cells.jsonl"
    _write_bytes_atomic(
        gap_path, gap_text.encode("utf-8", errors="strict"),
        evidence_directory=evidence_directory)
    metadata = {
        "run_contract": contract,
        "run_contract_sha256": contract_digest,
        "completion": completion,
    }
    summary = summarize(results, metadata)
    write_json_atomic(
        output / "summary.json", summary, pretty=True,
        evidence_directory=evidence_directory)
    print(output / "summary.json")
    return 0 if not summary["failed_cell_count"] else 1


def _run_with_snapshot_owner(
        options: argparse.Namespace, lock_guard: CanonicalBenchmarkLock,
        snapshots: ExecutableSnapshotOwner) -> int:
    evidence_directory = AtomicEvidenceDirectory.open_or_create(
        options.output)
    try:
        return _run_with_snapshot_owner_held(
            options, lock_guard, snapshots, evidence_directory)
    finally:
        evidence_directory.close()


def run(options: argparse.Namespace,
        lock_guard: CanonicalBenchmarkLock) -> int:
    snapshots = ExecutableSnapshotOwner()
    try:
        return _run_with_snapshot_owner(options, lock_guard, snapshots)
    finally:
        snapshots.close()


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main", type=Path, required=True)
    parser.add_argument("--main-commit", default=MAIN_COMMIT)
    parser.add_argument("--main-sha256", default=PURE_AVX2_MAIN_SHA256,
                        help="required frozen pure-AVX2 baseline digest")
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--current-source-root", type=Path, required=True)
    parser.add_argument("--current-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int,
                        default=min(128, os.cpu_count() or 1))
    parser.add_argument(
        "--reproducible-build-jobs", type=int, default=1,
        help="parallel jobs for the one-time clean reproducibility rebuild; "
             "kept separate from benchmark workers to bound peak memory")
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--gf8-only", action="store_true",
                        help="run the 2,522-cell GF8 K=1..255 matrix only")
    parser.add_argument("--lock", type=Path,
                        default=AUTHORITATIVE_LOCK,
                        help="canonical exclusive lock shared by every "
                             "build/test/timing lane (cannot be overridden)")
    parser.add_argument("--with-current-legacy", action="store_true",
                        help="also time the current tree's retained legacy API")
    options = parser.parse_args(arguments)
    with benchmark_lock(options.lock) as lock_guard:
        return run(options, lock_guard)


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == BOUNDED_SUPERVISOR_MODE:
        raise SystemExit(_bounded_process_supervisor(sys.argv[2:]))
    raise SystemExit(main())
