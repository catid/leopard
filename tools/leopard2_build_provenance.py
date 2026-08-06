#!/usr/bin/env python3
"""Fail-closed provenance for production Leopard2 benchmark builds.

The benchmark's embedded source attestation covers only the benchmark
translation unit.  It cannot prove which objects were archived into Leopard2.
This module binds a benchmark to its exact Release CMake cache, compiler
commands, source/object pairs, archive recipe/members, link recipe, and output
bytes.  The final-map and all-K runners both use the same closure.
"""

from __future__ import annotations

from collections import Counter
from contextlib import ExitStack, nullcontext
import ctypes
import errno
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import select
import selectors
import shlex
import shutil
import signal
import stat
import struct
import subprocess
import sys
import tempfile
import time
from typing import Any, Mapping, Sequence
import unicodedata


MAX_FILE_BYTES = 256 * 1024 * 1024
MAX_METADATA_BYTES = 16 * 1024 * 1024
MAX_METADATA_JSON_DEPTH = 64
MAX_TRACKED_SOURCE_FILES = 16 * 1024
MAX_TRACKED_SOURCE_FILE_BYTES = 64 * 1024 * 1024
MAX_TRACKED_SOURCE_TOTAL_BYTES = 512 * 1024 * 1024
MAX_TRACKED_SOURCE_PATH_BYTES = 4096
MAX_GITFILE_BYTES = 4096
MAX_REPLAY_RECIPE_FILES = 4096
MAX_REPLAY_RECIPE_TOTAL_BYTES = 128 * 1024 * 1024
CHILD_REAP_TIMEOUT_SECONDS = 5.0
PR_SET_CHILD_SUBREAPER = 36
PR_GET_CHILD_SUBREAPER = 37
PR_SET_NO_NEW_PRIVS = 38

# Linux Landlock ABI 3 adds TRUNCATE mediation.  Canonical replay commands use
# it to make every pathname write outside their runner-owned workspace fail,
# even if a same-UID adversary transiently replaces an output with a symlink.
LANDLOCK_CREATE_RULESET_VERSION = 1
LANDLOCK_RULE_PATH_BENEATH = 1
LANDLOCK_ACCESS_FS_WRITE_FILE = 1 << 1
LANDLOCK_ACCESS_FS_REMOVE_DIR = 1 << 4
LANDLOCK_ACCESS_FS_REMOVE_FILE = 1 << 5
LANDLOCK_ACCESS_FS_MAKE_CHAR = 1 << 6
LANDLOCK_ACCESS_FS_MAKE_DIR = 1 << 7
LANDLOCK_ACCESS_FS_MAKE_REG = 1 << 8
LANDLOCK_ACCESS_FS_MAKE_SOCK = 1 << 9
LANDLOCK_ACCESS_FS_MAKE_FIFO = 1 << 10
LANDLOCK_ACCESS_FS_MAKE_BLOCK = 1 << 11
LANDLOCK_ACCESS_FS_MAKE_SYM = 1 << 12
LANDLOCK_ACCESS_FS_REFER = 1 << 13
LANDLOCK_ACCESS_FS_TRUNCATE = 1 << 14
LANDLOCK_WRITE_ACCESS = (
    LANDLOCK_ACCESS_FS_WRITE_FILE |
    LANDLOCK_ACCESS_FS_REMOVE_DIR |
    LANDLOCK_ACCESS_FS_REMOVE_FILE |
    LANDLOCK_ACCESS_FS_MAKE_CHAR |
    LANDLOCK_ACCESS_FS_MAKE_DIR |
    LANDLOCK_ACCESS_FS_MAKE_REG |
    LANDLOCK_ACCESS_FS_MAKE_SOCK |
    LANDLOCK_ACCESS_FS_MAKE_FIFO |
    LANDLOCK_ACCESS_FS_MAKE_BLOCK |
    LANDLOCK_ACCESS_FS_MAKE_SYM |
    LANDLOCK_ACCESS_FS_REFER |
    LANDLOCK_ACCESS_FS_TRUNCATE
)
SYS_LANDLOCK_CREATE_RULESET = 444
SYS_LANDLOCK_ADD_RULE = 445
SYS_LANDLOCK_RESTRICT_SELF = 446
CLONE_NEWNS = 0x00020000
CLONE_NEWUSER = 0x10000000
MS_NOSUID = 0x00000002
MS_NODEV = 0x00000004
MS_REC = 0x00004000
MS_PRIVATE = 0x00040000
PRIVATE_REPLAY_TMPFS_BYTES = 512 * 1024 * 1024

# Linux inotify constants.  Provenance command containment is already
# Linux-only; inotify gives retained pathname guards an event history, which
# inode/mtime checks alone cannot provide after an A -> B -> A rename.
IN_MODIFY = 0x00000002
IN_ATTRIB = 0x00000004
IN_CLOSE_WRITE = 0x00000008
IN_MOVED_FROM = 0x00000040
IN_MOVED_TO = 0x00000080
IN_CREATE = 0x00000100
IN_DELETE = 0x00000200
IN_DELETE_SELF = 0x00000400
IN_MOVE_SELF = 0x00000800
IN_UNMOUNT = 0x00002000
IN_Q_OVERFLOW = 0x00004000
IN_IGNORED = 0x00008000
IN_ONLYDIR = 0x01000000
IN_DONT_FOLLOW = 0x02000000
IN_EXCL_UNLINK = 0x04000000
IN_MASK_ADD = 0x20000000
IN_CLOEXEC = getattr(os, "O_CLOEXEC", 0)
IN_NONBLOCK = getattr(os, "O_NONBLOCK", 0)
INOTIFY_EVENT = struct.Struct("iIII")

# A retained CMake recipe is executed by Make through a shell.  shlex splitting
# alone is therefore not a security boundary: ``object.o&&evil`` is one shlex
# token but two shell commands.  The canonical production recipes need none of
# these expansion, substitution, globbing, redirection, comment, or control
# characters, so reject them wherever they occur, including inside an option.
SHELL_RECIPE_META = frozenset(";&|<>()`$*?[]{}~#")

CMAKE_CACHE_ENTRY_TYPES = frozenset((
    "BOOL", "FILEPATH", "INTERNAL", "PATH", "STATIC", "STRING",
    "UNINITIALIZED",
))

# The production evidence contract below gives these cache entries semantic
# meaning.  Bind that meaning to the exact type emitted by the canonical
# configure rather than accepting a forged value under an unrelated CMake
# type.  Entries that this helper does not consume still receive syntax and
# known-type validation, but do not acquire a new semantic requirement here.
CMAKE_CACHE_REQUIRED_ENTRY_TYPES = {
    "CMAKE_AR": frozenset(("FILEPATH",)),
    "CMAKE_BUILD_TYPE": frozenset(("STRING",)),
    "CMAKE_C_COMPILER": frozenset(("FILEPATH",)),
    "CMAKE_C_COMPILER_ARG1": frozenset(("STRING",)),
    "CMAKE_C_COMPILER_LAUNCHER": frozenset(("STRING",)),
    "CMAKE_C_COMPILER_TARGET": frozenset(("STRING",)),
    "CMAKE_CXX_COMPILER": frozenset(("FILEPATH",)),
    "CMAKE_CXX_COMPILER_ARG1": frozenset(("STRING",)),
    "CMAKE_CXX_COMPILER_LAUNCHER": frozenset(("STRING",)),
    "CMAKE_CXX_COMPILER_TARGET": frozenset(("STRING",)),
    "CMAKE_CXX_FLAGS": frozenset(("STRING",)),
    "CMAKE_CXX_FLAGS_RELEASE": frozenset(("STRING",)),
    "CMAKE_EXPORT_COMPILE_COMMANDS": frozenset(("UNINITIALIZED",)),
    "CMAKE_GENERATOR": frozenset(("INTERNAL",)),
    "CMAKE_HOME_DIRECTORY": frozenset(("INTERNAL",)),
    "CMAKE_LINKER": frozenset(("FILEPATH",)),
    "CMAKE_MAKE_PROGRAM": frozenset(("FILEPATH",)),
    "CMAKE_POSITION_INDEPENDENT_CODE": frozenset(("BOOL",)),
    "CMAKE_RANLIB": frozenset(("FILEPATH",)),
    "CMAKE_SYSROOT": frozenset(("PATH",)),
    "CMAKE_TOOLCHAIN_FILE": frozenset(("FILEPATH",)),
    "ENABLE_OPENMP": frozenset(("BOOL",)),
    "LEOPARD_ENABLE_GF8": frozenset(("BOOL",)),
    "LEOPARD_ENABLE_GF16": frozenset(("BOOL",)),
    "LEO2_BACKEND_VARIANT": frozenset(("STRING",)),
    "LEO2_BENCHMARK_GIT_EXECUTABLE": frozenset(("FILEPATH",)),
    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
        frozenset(("INTERNAL",)),
    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
        frozenset(("INTERNAL",)),
    "LEO2_BUILD_ALLK_DIAGNOSTIC": frozenset(("BOOL",)),
    "LEO2_BUILD_BENCHMARKS": frozenset(("BOOL",)),
    "LEO2_BUILD_FUZZERS": frozenset(("BOOL",)),
    "LEO2_BUILD_TESTS": frozenset(("BOOL",)),
    "LEO2_ENABLE_CUDA": frozenset(("BOOL",)),
    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT": frozenset(("BOOL",)),
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": frozenset(("BOOL",)),
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED":
        frozenset(("BOOL",)),
    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK":
        frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT":
        frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": frozenset(("STRING",)),
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK": frozenset(("BOOL",)),
    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT":
        frozenset(("BOOL",)),
    "LEO2_FLAG_FALIGN_FUNCTIONS_64": frozenset(("INTERNAL",)),
    "LEO2_FLAG_MAVX2": frozenset(("INTERNAL",)),
    "LEO2_FLAG_MGFNI": frozenset(("INTERNAL",)),
    "LEO2_FLAG_MNO_AVX512F": frozenset(("INTERNAL",)),
    "LEO2_LOCATOR_GIT_EXECUTABLE": frozenset(("FILEPATH",)),
    # The pure-AVX2 evidence configure deliberately seeds these otherwise
    # successful probe variables to FALSE on its command line.
    "LEO2_FLAG_MAVX512F": frozenset(("UNINITIALIZED",)),
    "LEO2_FLAG_MAVX512BW": frozenset(("UNINITIALIZED",)),
    "LEO2_FLAG_MAVX512VL": frozenset(("UNINITIALIZED",)),
}

CORE_LIBRARY_SOURCES = {
    "leopard.cpp",
    "leopard2.cpp",
    "Leopard2Backend.cpp",
    "Leopard2BackendScalar.cpp",
    "Leopard2CpuFeatures.cpp",
    "Leopard2Plan.cpp",
    "LeopardCommon.cpp",
}

UNCONDITIONAL_AVX2_LIBRARY_SOURCES_V5 = frozenset((
    "Leopard2BackendAVX2.cpp",
    "Leopard2BackendAVX2Xor.cpp",
))
UNCONDITIONAL_AVX2_LIBRARY_SOURCES_V6 = frozenset((
    *UNCONDITIONAL_AVX2_LIBRARY_SOURCES_V5,
    "Leopard2BackendAVX2T2K4.cpp",
))
UNCONDITIONAL_AVX2_LIBRARY_SOURCES = frozenset((
    *UNCONDITIONAL_AVX2_LIBRARY_SOURCES_V6,
    "Leopard2BackendAVX2T8K8B1024.cpp",
))
CONDITIONAL_AVX2_LIBRARY_SOURCES = frozenset((
    "Leopard2BackendAVX2T16B64.cpp",
    "Leopard2BackendAVX2T16Q2.cpp",
    "Leopard2BackendAVX2T32B256.cpp",
))
PRODUCTION_AVX2_LIBRARY_SOURCES = (
    UNCONDITIONAL_AVX2_LIBRARY_SOURCES |
    CONDITIONAL_AVX2_LIBRARY_SOURCES)

GIT_ENVIRONMENT = {
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_NO_REPLACE_OBJECTS": "1",
    "GIT_OPTIONAL_LOCKS": "0",
    "LANG": "C",
    "LC_ALL": "C",
    "PATH": "/usr/bin:/bin",
}

PRODUCTION_BUILD_CLOSURE_SCHEMA_V1 = \
    "leopard2-production-build-closure/v1"
PRODUCTION_BUILD_CLOSURE_SCHEMA = \
    "leopard2-production-build-closure/v2"
BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V2 = \
    "leopard2-benchmark-build-configuration/v2"
BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3 = \
    "leopard2-benchmark-build-configuration/v3"
BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4 = \
    "leopard2-benchmark-build-configuration/v4"
BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5 = \
    "leopard2-benchmark-build-configuration/v5"
BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6 = \
    "leopard2-benchmark-build-configuration/v6"
BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7 = \
    "leopard2-benchmark-build-configuration/v7"
BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8 = \
    "leopard2-benchmark-build-configuration/v8"
BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9 = \
    "leopard2-benchmark-build-configuration/v9"
BENCHMARK_BUILD_CONFIGURATION_SCHEMA = \
    "leopard2-benchmark-build-configuration/v10"
REPRODUCIBLE_BUILD_PROOF_SCHEMA_V2 = \
    "leopard2-reproducible-build-proof/v2"
REPRODUCIBLE_BUILD_PROOF_SCHEMA = \
    "leopard2-reproducible-build-proof/v3"
REPLAY_INVOCATION_SCHEMA_V1 = "leopard2-replay-invocation/v1"
REPLAY_INVOCATION_SCHEMA = "leopard2-replay-invocation/v2"
LEGACY_REPLAY_RECIPE_SCHEMA = "leopard2-replay-recipe-transport/v2"
CANONICAL_REPLAY_RECIPE_SCHEMA = \
    "leopard2-canonical-replay-build-plan/v2"


class BuildProvenanceError(RuntimeError):
    """The candidate build is not valid benchmark evidence."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise BuildProvenanceError(message)


_OWNER_TERMINAL_EXCEPTIONS = (KeyboardInterrupt, SystemExit)


def _add_owner_failure_note(
    target: BaseException,
    context: str,
    relation: str,
    failure: BaseException,
) -> None:
    """Attach one cleanup diagnostic and its already-retained subfailures."""
    try:
        target.add_note(
            f"{context}: {relation} was "
            f"{type(failure).__name__}: {failure}")
        for note in getattr(failure, "__notes__", ()):
            target.add_note(f"{context}: nested cleanup detail: {note}")
    except AttributeError:
        pass


def _owner_exception_precedence(
    earlier: BaseException | None,
    later: BaseException | None,
    context: str,
) -> BaseException | None:
    """Select one authoritative owner failure without losing termination.

    KeyboardInterrupt and SystemExit are control-flow requests, not ordinary
    diagnostics.  An earlier terminal exception remains authoritative through
    all later cleanup.  A later cleanup terminal supersedes an earlier
    ordinary failure.  Ordinary-only callers retain their existing first-error
    behavior and may wrap the result at their public boundary.
    """
    require(isinstance(context, str) and context,
            "owner exception-precedence context is invalid")
    if earlier is None:
        return later
    if later is None:
        return earlier
    if earlier is later:
        # Nested owner cleanup can catch and re-present the already-selected
        # terminal exception.  Adding that exception's notes back onto itself
        # would extend the same list while iterating it and grow without bound.
        return earlier
    earlier_terminal = isinstance(earlier, _OWNER_TERMINAL_EXCEPTIONS)
    later_terminal = isinstance(later, _OWNER_TERMINAL_EXCEPTIONS)
    if earlier_terminal:
        _add_owner_failure_note(
            earlier, context, "later cleanup failure", later)
        return earlier
    if later_terminal:
        _add_owner_failure_note(
            later, context, "earlier failure", earlier)
        return later
    _add_owner_failure_note(
        earlier, context, "later cleanup failure", later)
    return earlier


def _ordinary_owner_cleanup_error(
    cleanup: BaseException, message: str,
) -> BaseException:
    """Retain terminal cleanup verbatim and wrap an ordinary close failure."""
    if isinstance(cleanup, _OWNER_TERMINAL_EXCEPTIONS):
        return cleanup
    failure = BuildProvenanceError(message)
    failure.__cause__ = cleanup
    return failure


def _raise_owner_exit_failure(
    primary: BaseException,
    cleanup: BaseException,
    context: str,
    ordinary_message: str,
) -> None:
    """Apply owner precedence when cleanup fails during exception unwinding."""
    authoritative = _owner_exception_precedence(primary, cleanup, context)
    if isinstance(authoritative, _OWNER_TERMINAL_EXCEPTIONS):
        raise authoritative.with_traceback(authoritative.__traceback__)
    failure = BuildProvenanceError(ordinary_message)
    raise failure from cleanup


class _LandlockRulesetAttr(ctypes.Structure):
    _fields_ = (("handled_access_fs", ctypes.c_uint64),)


class _LandlockPathBeneathAttr(ctypes.Structure):
    _pack_ = 1
    _fields_ = (
        ("allowed_access", ctypes.c_uint64),
        ("parent_fd", ctypes.c_int32),
    )


def _landlock_restrict_writes_to(
    root: Path | str, retained_root_descriptor: int,
    additional_write_descriptors: Sequence[int] = (),
) -> None:
    """Restrict this process and its descendants to writes below ``root``.

    This runs in the already single-threaded subprocess child immediately
    before exec.  Reads and execution remain unrestricted; every pathname
    mutation class supported by Landlock ABI 3 is handled and allowed only
    below exact O_PATH directories retained by the parent.  Linux accepts only
    directory FDs as LANDLOCK_RULE_PATH_BENEATH parents.  Additional regular
    descriptors must therefore be writable anonymous memfds that were already
    opened by the trusted parent; Landlock intentionally does not revoke
    writes through such an inherited FD.  This is the descriptor-only
    transport used for private replay artifact sinks.
    """
    requested = Path(os.path.abspath(os.fspath(root)))
    require(sys.platform.startswith("linux") and requested.is_absolute(),
            "canonical replay write sandbox requires Linux and an absolute root")
    require(type(retained_root_descriptor) is int and
            retained_root_descriptor >= 0,
            "canonical replay write sandbox descriptor is invalid")
    try:
        retained_status = os.fstat(retained_root_descriptor)
        lexical_status = os.stat(requested, follow_symlinks=False)
    except OSError as error:
        raise BuildProvenanceError(
            f"canonical replay write sandbox root is unavailable: {error}") \
            from error
    require(
        stat.S_ISDIR(retained_status.st_mode) and
        (retained_status.st_dev, retained_status.st_ino) ==
        (lexical_status.st_dev, lexical_status.st_ino),
        "canonical replay write sandbox root changed before restriction")
    library = ctypes.CDLL(None, use_errno=True)
    syscall = library.syscall
    syscall.restype = ctypes.c_long
    ctypes.set_errno(0)
    abi = syscall(
        ctypes.c_long(SYS_LANDLOCK_CREATE_RULESET),
        ctypes.c_void_p(), ctypes.c_size_t(0),
        ctypes.c_uint(LANDLOCK_CREATE_RULESET_VERSION))
    if abi < 3:
        number = ctypes.get_errno()
        raise BuildProvenanceError(
            "canonical replay requires Landlock ABI 3 or later: " +
            (os.strerror(number) if abi < 0 and number else f"ABI {abi}"))

    ruleset = -1
    try:
        attributes = _LandlockRulesetAttr(LANDLOCK_WRITE_ACCESS)
        ctypes.set_errno(0)
        ruleset = syscall(
            ctypes.c_long(SYS_LANDLOCK_CREATE_RULESET),
            ctypes.byref(attributes), ctypes.sizeof(attributes),
            ctypes.c_uint(0))
        if ruleset < 0:
            number = ctypes.get_errno()
            raise BuildProvenanceError(
                "cannot create canonical replay Landlock ruleset: " +
                os.strerror(number or errno.EPERM))
        rules = [(retained_root_descriptor, LANDLOCK_WRITE_ACCESS)]
        for descriptor in additional_write_descriptors:
            status = os.fstat(descriptor)
            require(
                stat.S_ISREG(status.st_mode) or stat.S_ISDIR(status.st_mode),
                "canonical replay Landlock output descriptor is unsafe")
            if stat.S_ISDIR(status.st_mode):
                rules.append((descriptor, LANDLOCK_WRITE_ACCESS))
                continue
            flags = fcntl.fcntl(descriptor, fcntl.F_GETFL)
            try:
                descriptor_name = os.readlink(
                    f"/proc/self/fd/{descriptor}")
                fcntl.fcntl(
                    descriptor, getattr(fcntl, "F_GET_SEALS", 1034))
            except OSError as error:
                raise BuildProvenanceError(
                    "canonical replay regular output descriptor is not an "
                    f"anonymous sealable memfd: {error}") from error
            require(
                status.st_nlink == 0 and
                descriptor_name.startswith("/memfd:") and
                descriptor_name.endswith(" (deleted)") and
                flags & os.O_ACCMODE in (os.O_WRONLY, os.O_RDWR),
                "canonical replay regular output descriptor is not a "
                "writable anonymous memfd")
        for descriptor, allowed_access in rules:
            rule = _LandlockPathBeneathAttr(allowed_access, descriptor)
            ctypes.set_errno(0)
            result = syscall(
                ctypes.c_long(SYS_LANDLOCK_ADD_RULE),
                ctypes.c_int(ruleset),
                ctypes.c_int(LANDLOCK_RULE_PATH_BENEATH),
                ctypes.byref(rule), ctypes.c_uint(0))
            if result != 0:
                number = ctypes.get_errno()
                raise BuildProvenanceError(
                    "cannot add canonical replay Landlock path rule: " +
                    os.strerror(number or errno.EPERM))
        prctl = library.prctl
        prctl.restype = ctypes.c_int
        ctypes.set_errno(0)
        if prctl(
                ctypes.c_int(PR_SET_NO_NEW_PRIVS), ctypes.c_ulong(1),
                ctypes.c_ulong(0), ctypes.c_ulong(0),
                ctypes.c_ulong(0)) != 0:
            number = ctypes.get_errno()
            raise BuildProvenanceError(
                "cannot enable no_new_privs for canonical replay: " +
                os.strerror(number or errno.EPERM))
        ctypes.set_errno(0)
        result = syscall(
            ctypes.c_long(SYS_LANDLOCK_RESTRICT_SELF),
            ctypes.c_int(ruleset), ctypes.c_uint(0))
        if result != 0:
            number = ctypes.get_errno()
            raise BuildProvenanceError(
                "cannot enter canonical replay Landlock ruleset: " +
                os.strerror(number or errno.EPERM))
    finally:
        if ruleset >= 0:
            _close_raw_descriptor_with_precedence(
                ruleset, sys.exception(),
                "canonical replay Landlock ruleset",
                "cannot close canonical replay Landlock ruleset")


def _mount_private_replay_tmpfs(
    root: Path | str, retained_root_descriptor: int,
    directory_relatives: Sequence[str], outside_uid: int, outside_gid: int,
) -> None:
    """Enter a private user/mount namespace and cover ``root`` with tmpfs."""
    requested = Path(os.path.abspath(os.fspath(root)))
    require(
        sys.platform.startswith("linux") and requested.is_absolute() and
        type(retained_root_descriptor) is int and
        retained_root_descriptor >= 0 and
        type(outside_uid) is int and outside_uid >= 0 and
        type(outside_gid) is int and outside_gid >= 0 and
        all(
            isinstance(relative, str) and
            relative not in ("", ".", "..") and
            not Path(relative).is_absolute() and
            all(component not in ("", ".", "..")
                for component in Path(relative).parts)
            for relative in directory_relatives),
        "private replay tmpfs request is malformed")
    retained_status = os.fstat(retained_root_descriptor)
    lexical_status = os.stat(requested, follow_symlinks=False)
    require(
        stat.S_ISDIR(retained_status.st_mode) and
        (retained_status.st_dev, retained_status.st_ino) ==
        (lexical_status.st_dev, lexical_status.st_ino),
        "private replay tmpfs root changed before namespace creation")
    library = ctypes.CDLL(None, use_errno=True)
    unshare = library.unshare
    unshare.argtypes = (ctypes.c_int,)
    unshare.restype = ctypes.c_int
    ctypes.set_errno(0)
    if unshare(ctypes.c_int(CLONE_NEWUSER | CLONE_NEWNS)) != 0:
        number = ctypes.get_errno()
        raise BuildProvenanceError(
            "cannot create private replay user/mount namespace: " +
            os.strerror(number or errno.EPERM))

    def write_mapping(path: str, content: bytes, *, optional: bool = False) -> None:
        descriptor = -1
        try:
            descriptor = os.open(
                path, os.O_WRONLY | getattr(os, "O_CLOEXEC", 0))
            view = memoryview(content)
            while view:
                written = os.write(descriptor, view)
                require(written > 0,
                        f"private replay namespace mapping write stalled: "
                        f"{path}")
                view = view[written:]
        except FileNotFoundError:
            if not optional:
                raise
        finally:
            if descriptor >= 0:
                _close_raw_descriptor_with_precedence(
                    descriptor, sys.exception(),
                    f"private replay namespace mapping {path}",
                    f"cannot close private replay namespace mapping {path}")

    try:
        write_mapping("/proc/self/setgroups", b"deny\n", optional=True)
        write_mapping(
            "/proc/self/uid_map", f"0 {outside_uid} 1\n".encode("ascii"))
        write_mapping(
            "/proc/self/gid_map", f"0 {outside_gid} 1\n".encode("ascii"))
        os.setresgid(0, 0, 0)
        os.setresuid(0, 0, 0)
    except OSError as error:
        raise BuildProvenanceError(
            f"cannot map private replay namespace credentials: {error}") \
            from error

    mount = library.mount
    mount.argtypes = (
        ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
        ctypes.c_ulong, ctypes.c_void_p,
    )
    mount.restype = ctypes.c_int
    ctypes.set_errno(0)
    if mount(
            ctypes.c_char_p(), ctypes.c_char_p(b"/"), ctypes.c_char_p(),
            ctypes.c_ulong(MS_REC | MS_PRIVATE), ctypes.c_void_p()) != 0:
        number = ctypes.get_errno()
        raise BuildProvenanceError(
            "cannot privatize replay mount propagation: " +
            os.strerror(number or errno.EPERM))
    # mount(2) rejects procfs magic-link targets on some kernels.  The exact
    # inode was retained and compared above; mount in this new namespace and
    # compare the covered lexical path immediately afterwards.
    target = os.fsencode(requested)
    options = (
        f"mode=0700,size={PRIVATE_REPLAY_TMPFS_BYTES},nr_inodes=8192"
    ).encode("ascii")
    ctypes.set_errno(0)
    if mount(
            ctypes.c_char_p(b"leopard2-private-replay"),
            ctypes.c_char_p(target), ctypes.c_char_p(b"tmpfs"),
            ctypes.c_ulong(MS_NOSUID | MS_NODEV),
            ctypes.c_char_p(options)) != 0:
        number = ctypes.get_errno()
        raise BuildProvenanceError(
            "cannot mount private replay tmpfs: " +
            os.strerror(number or errno.EPERM))
    try:
        mounted = os.stat(requested, follow_symlinks=False)
        require(
            stat.S_ISDIR(mounted.st_mode) and
            mounted.st_dev != retained_status.st_dev,
            "private replay tmpfs did not cover its staging root")
        mounted_descriptor = os.open(
            requested,
            getattr(os, "O_PATH", os.O_RDONLY) |
            getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_DIRECTORY", 0) |
            getattr(os, "O_NOFOLLOW", 0))
        try:
            mounted_status = os.fstat(mounted_descriptor)
            require(
                mounted_status.st_dev == mounted.st_dev and
                mounted_status.st_ino == mounted.st_ino,
                "private replay mounted-root descriptor differs")
            os.dup2(
                mounted_descriptor, retained_root_descriptor,
                inheritable=True)
        finally:
            _close_raw_descriptor_with_precedence(
                mounted_descriptor, sys.exception(),
                "private replay mounted-root descriptor",
                "cannot close private replay mounted-root descriptor")
        for relative in sorted(set(directory_relatives)):
            (Path(f"/proc/self/fd/{retained_root_descriptor}") /
             relative).mkdir(
                mode=0o700, parents=True, exist_ok=False)
    except OSError as error:
        raise BuildProvenanceError(
            f"cannot initialize private replay tmpfs: {error}") from error


def _linux_memfd_create(name: str, flags: int) -> int:
    """Create a memfd even when this CPython omitted the os wrapper."""
    wrapper = getattr(os, "memfd_create", None)
    if callable(wrapper):
        return wrapper(name, flags)
    require(sys.platform.startswith("linux"),
            "immutable executable snapshots require Linux memfd")
    try:
        function = ctypes.CDLL(None, use_errno=True).memfd_create
    except (AttributeError, OSError) as error:
        raise BuildProvenanceError(
            f"Linux memfd_create is unavailable: {error}") from error
    function.argtypes = (ctypes.c_char_p, ctypes.c_uint)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    descriptor = function(name.encode("utf-8"), ctypes.c_uint(flags))
    if descriptor >= 0:
        return descriptor
    number = ctypes.get_errno()
    raise BuildProvenanceError(
        "cannot create immutable executable memfd: " +
        os.strerror(number or errno.EPERM))


KCMP_FILE = 0
_SYS_KCMP_BY_MACHINE = {
    "x86_64": 312,
    "amd64": 312,
    "i386": 349,
    "i486": 349,
    "i586": 349,
    "i686": 349,
    "armv6l": 378,
    "armv7l": 378,
    "aarch64": 272,
    "arm64": 272,
    "riscv64": 272,
    "loongarch64": 272,
    "ppc": 354,
    "ppc64": 354,
    "ppc64le": 354,
    "s390": 343,
    "s390x": 343,
}


def _linux_same_open_file_description(
    first: int, second: int,
) -> bool:
    """Compare exact open-file descriptions, not merely inode metadata."""
    require(sys.platform.startswith("linux"),
            "interruption-safe descriptor closure requires Linux")
    try:
        os.fstat(first)
    except OSError as error:
        if error.errno == errno.EBADF:
            return False
        raise
    machine = os.uname().machine.lower()
    number = _SYS_KCMP_BY_MACHINE.get(machine)
    require(number is not None,
            f"Linux kcmp syscall number is unknown for {machine!r}")
    syscall = ctypes.CDLL(None, use_errno=True).syscall
    syscall.restype = ctypes.c_long
    ctypes.set_errno(0)
    result = syscall(
        ctypes.c_long(number),
        ctypes.c_int(os.getpid()), ctypes.c_int(os.getpid()),
        ctypes.c_int(KCMP_FILE),
        ctypes.c_ulong(first), ctypes.c_ulong(second))
    if result >= 0:
        return result == 0
    observed = ctypes.get_errno()
    if observed == errno.EBADF:
        return False
    raise BuildProvenanceError(
        "cannot compare retained open-file descriptions with Linux kcmp: " +
        os.strerror(observed or errno.EPERM))


def _linux_close_descriptor_once(descriptor: int) -> None:
    """Issue one non-retried Linux close syscall.

    Linux releases the numeric descriptor early in close(2).  Retrying after
    any reported error can therefore close an unrelated descriptor recycled
    by asynchronous code.  This helper deliberately issues exactly one raw
    syscall and never retries its numeric argument.
    """
    require(sys.platform.startswith("linux"),
            "interruption-safe descriptor closure requires Linux")
    function = ctypes.CDLL(None, use_errno=True).close
    function.argtypes = (ctypes.c_int,)
    function.restype = ctypes.c_int
    ctypes.set_errno(0)
    result = function(ctypes.c_int(descriptor))
    if result == 0:
        return
    observed = ctypes.get_errno()
    if observed == errno.EBADF:
        return
    raise OSError(
        observed or errno.EIO,
        os.strerror(observed or errno.EIO))


def _close_descriptor_with_ofd_guard(
    descriptor: int,
    context: str,
    ownership_consumed: list[bool],
) -> None:
    """Close one FD while an OFD duplicate disambiguates close-time ABA.

    An async handler can run after close(2), reopen the same inode into the
    recycled numeric slot, and then raise.  stat fields cannot distinguish
    that new open-file description from the one this owner held.  Keep a
    duplicate of the original description until the first close returns.  On
    an exception, Linux kcmp says whether the numeric slot still aliases that
    duplicate.  Only a proven alias is force-closed; a recycled slot is left
    untouched.  The duplicate is then released with one non-retried raw close.
    """
    require(type(descriptor) is int and descriptor >= 0 and
            isinstance(context, str) and context and
            isinstance(ownership_consumed, list) and
            ownership_consumed == [False],
            f"{context} descriptor-close guard arguments are invalid")
    guard = os.dup(descriptor)
    failure: BaseException | None = None
    candidate_is_owned = False
    try:
        os.close(descriptor)
    except BaseException as error:
        failure = error
        try:
            candidate_is_owned = _linux_same_open_file_description(
                descriptor, guard)
        except BaseException as comparison_error:
            failure = _owner_exception_precedence(
                failure, comparison_error,
                f"{context} open-file-description comparison")

    if candidate_is_owned:
        try:
            _linux_close_descriptor_once(descriptor)
        except BaseException as cleanup_error:
            failure = _owner_exception_precedence(
                failure, cleanup_error,
                f"{context} retained descriptor close")
    try:
        _linux_close_descriptor_once(guard)
    except BaseException as cleanup_error:
        failure = _owner_exception_precedence(
            failure, cleanup_error,
            f"{context} guard descriptor close")
    ownership_consumed[0] = True
    if failure is not None:
        raise failure.with_traceback(failure.__traceback__)


class _OwnedDescriptor:
    """Own one raw descriptor across interruption-prone Python handoffs.

    Directory descriptors cannot be wrapped by ``os.fdopen``.  Keeping their
    identity with the integer lets a retry distinguish an already-closed (or
    recycled) descriptor from the object originally acquired.
    """

    def __init__(self) -> None:
        self._descriptor = -1
        self._identity: tuple[int, int, int] | None = None

    @property
    def descriptor(self) -> int:
        return self._descriptor

    def open(self, path: Path | str, flags: int, mode: int = 0o777,
             *, dir_fd: int | None = None) -> int:
        require(self._descriptor < 0,
                "owned descriptor was opened more than once")
        if dir_fd is None:
            self._descriptor = os.open(path, flags, mode)
        else:
            self._descriptor = os.open(path, flags, mode, dir_fd=dir_fd)
        try:
            status = os.fstat(self._descriptor)
        except BaseException as primary:
            descriptor = self._descriptor
            consumed = [False]
            try:
                _close_descriptor_with_ofd_guard(
                    descriptor, "owned descriptor open cleanup", consumed)
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error, "owned descriptor open",
                    "owned descriptor open cleanup failed: "
                    f"{cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            finally:
                if consumed[0]:
                    self._descriptor = -1
                    self._identity = None
            raise
        self._identity = (status.st_dev, status.st_ino, status.st_mode)
        return self._descriptor

    def close(self) -> None:
        descriptor = self._descriptor
        if descriptor < 0:
            return
        try:
            status = os.fstat(descriptor)
        except OSError as error:
            if error.errno == errno.EBADF:
                self._descriptor = -1
                self._identity = None
                return
            raise
        current = (status.st_dev, status.st_ino, status.st_mode)
        if self._identity is not None and current != self._identity:
            # An interruption handler reused the numeric descriptor after the
            # original close completed.  It is no longer ours to close.
            self._descriptor = -1
            self._identity = None
            return
        consumed = [False]
        try:
            _close_descriptor_with_ofd_guard(
                descriptor, "owned descriptor", consumed)
        finally:
            if consumed[0]:
                self._descriptor = -1
                self._identity = None

    def __enter__(self) -> "_OwnedDescriptor":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, tb
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc is None:
                raise
            _raise_owner_exit_failure(
                exc, cleanup_error, "owned descriptor",
                f"owned descriptor cleanup failed: {cleanup_error}; "
                f"primary failure: {type(exc).__name__}: {exc}")


def _close_owned_descriptor_after_failure(
    owner: _OwnedDescriptor,
    primary: BaseException,
    context: str,
) -> None:
    """Close one partially published descriptor under owner precedence."""
    _close_owned_descriptors_after_failure(
        (owner,), primary, context)


def _close_owned_descriptors_after_failure(
    owners: Sequence[_OwnedDescriptor],
    primary: BaseException,
    context: str,
) -> None:
    """Close every partial owner before applying terminal precedence."""
    cleanup_failure: BaseException | None = None
    cleanup_details: list[str] = []
    for owner in owners:
        try:
            owner.close()
        except BaseException as error:
            cleanup_details.append(f"{type(error).__name__}: {error}")
            observed = _ordinary_owner_cleanup_error(
                error, f"{context} cleanup failed: {error}")
            cleanup_failure = _owner_exception_precedence(
                cleanup_failure, observed, f"{context} cleanup")
    if cleanup_failure is None:
        return
    details = "; ".join(cleanup_details)
    _raise_owner_exit_failure(
        primary, cleanup_failure, context,
        f"{context} cleanup failed: {details}; primary failure: "
        f"{type(primary).__name__}: {primary}")


def _close_raw_descriptor_with_precedence(
    descriptor: int,
    primary: BaseException | None,
    context: str,
    ordinary_message: str,
    *,
    ownership_consumed: list[bool] | None = None,
) -> None:
    """Close one raw descriptor without masking active terminal control flow."""
    consumed = [False]
    try:
        _close_descriptor_with_ofd_guard(
            descriptor, context, consumed)
    except BaseException as cleanup_error:
        detail = f"{ordinary_message}: {cleanup_error}"
        observed = _ordinary_owner_cleanup_error(
            cleanup_error, detail)
        if primary is not None:
            _raise_owner_exit_failure(
                primary, observed, context,
                f"{detail}; primary failure: "
                f"{type(primary).__name__}: {primary}")
        if observed is cleanup_error:
            raise
        raise observed from cleanup_error
    finally:
        if ownership_consumed is not None and consumed[0]:
            ownership_consumed[0] = True


def _require_safe_unicode(value: str, label: str) -> None:
    """Reject characters that can reframe or visually spoof evidence.

    Cc covers C0, DEL, and C1 controls.  Cs covers unpaired surrogates.
    Cf includes bidi overrides/isolates and zero-width format controls; those
    are valid Unicode but unsafe in command/path evidence.  Zl/Zp are Unicode
    line delimiters and are rejected for the same framing reason.
    """
    forbidden_categories = frozenset(("Cc", "Cs", "Cf", "Zl", "Zp"))
    require(not any(
                unicodedata.category(character) in forbidden_categories
                for character in value),
            f"{label} contains an unsafe control, format, surrogate, or "
            "line-separator character")


def _strict_json_loads(
    encoded: bytes, label: str, *, maximum_depth: int = MAX_METADATA_JSON_DEPTH,
) -> Any:
    """Decode bounded metadata without accepting ambiguous JSON extensions."""
    require(0 <= len(encoded) <= MAX_METADATA_BYTES,
            f"{label} exceeds its retained byte bound")

    def object_from_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require_safe_unicode(key, label)
            require(key not in result, f"{label} contains duplicate key {key!r}")
            result[key] = value
        return result

    def finite_float(token: str) -> float:
        value = float(token)
        require(math.isfinite(value),
                f"{label} contains a non-finite number")
        return value

    def reject_constant(token: str) -> None:
        raise BuildProvenanceError(
            f"{label} contains non-standard constant {token!r}")

    try:
        text = encoded.decode("utf-8", errors="strict")
        value = json.loads(
            text,
            object_pairs_hook=object_from_pairs,
            parse_constant=reject_constant,
            parse_float=finite_float,
        )
    except BuildProvenanceError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, OverflowError,
            RecursionError, ValueError) as error:
        raise BuildProvenanceError(f"{label} is invalid: {error}") from error

    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        require(depth <= maximum_depth,
                f"{label} exceeds its maximum nesting depth")
        if isinstance(current, dict):
            stack.extend((item, depth + 1) for item in current.values())
        elif isinstance(current, list):
            stack.extend((item, depth + 1) for item in current)
        elif isinstance(current, str):
            _require_safe_unicode(current, label)
        elif isinstance(current, float):
            require(math.isfinite(current),
                    f"{label} contains a non-finite number")
    return value


def _stable_fields(status: os.stat_result) -> tuple[int, ...]:
    return (
        status.st_dev,
        status.st_ino,
        status.st_mode,
        status.st_uid,
        status.st_gid,
        status.st_nlink,
        status.st_size,
        status.st_mtime_ns,
        status.st_ctime_ns,
    )


class _InotifyMutationGuard:
    """Retain a fail-closed event history for exact pathnames or a tree."""

    _DIRECTORY_MASK = (
        IN_CLOSE_WRITE | IN_MOVED_FROM | IN_MOVED_TO | IN_CREATE | IN_DELETE |
        IN_DELETE_SELF | IN_MOVE_SELF | IN_UNMOUNT | IN_IGNORED |
        IN_EXCL_UNLINK
    )
    _FILE_MASK = (
        IN_MODIFY | IN_ATTRIB | IN_CLOSE_WRITE | IN_DELETE_SELF |
        IN_MOVE_SELF | IN_UNMOUNT | IN_IGNORED | IN_DONT_FOLLOW
    )

    def __init__(self, label: str) -> None:
        require(sys.platform.startswith("linux"),
                f"{label} pathname mutation guard requires Linux")
        self.label = label
        self._file: Any | None = None
        # wd -> None means every event is relevant.  Otherwise only directory
        # entry events naming one of the retained byte strings are relevant;
        # self/unmount/overflow events remain unconditionally relevant.
        self.rules: dict[int, set[bytes] | None] = {}
        self.mutations: list[str] = []
        descriptor = -1
        retained_file = None
        try:
            library = ctypes.CDLL(None, use_errno=True)
            function = library.inotify_init1
            function.argtypes = (ctypes.c_int,)
            function.restype = ctypes.c_int
            ctypes.set_errno(0)
            descriptor = function(
                ctypes.c_int(IN_CLOEXEC | IN_NONBLOCK))
            if descriptor < 0:
                number = ctypes.get_errno()
                raise OSError(number or errno.EPERM,
                              os.strerror(number or errno.EPERM))
            retained_file = os.fdopen(descriptor, "rb", buffering=0)
            descriptor = retained_file.fileno()
            self._library = library
            self._file = retained_file
        except BaseException as error:
            cleanup_error: BaseException | None = None
            try:
                if retained_file is None and descriptor >= 0:
                    _close_raw_descriptor_with_precedence(
                        descriptor, error,
                        f"{label} pathname mutation guard raw descriptor",
                        f"cannot close interrupted {label} pathname "
                        "mutation guard raw descriptor")
                elif retained_file is not None:
                    retained_file.close()
            except BaseException as observed:
                cleanup_error = observed
            if self._file is retained_file:
                self._file = None
            try:
                self._close_without_verification()
            except BaseException as observed:
                cleanup_error = _owner_exception_precedence(
                    cleanup_error, observed,
                    f"{label} pathname mutation guard constructor cleanup")
            if cleanup_error is not None:
                _raise_owner_exit_failure(
                    error, cleanup_error,
                    f"{label} pathname mutation guard constructor",
                    f"cannot close interrupted {label} pathname mutation "
                    f"guard: {cleanup_error}; primary failure: "
                    f"{type(error).__name__}: {error}")
            if not isinstance(error, (AttributeError, OSError)):
                raise
            raise BuildProvenanceError(
                f"cannot initialize {label} pathname mutation guard: "
                f"{error}") from error

    @property
    def descriptor(self) -> int:
        retained = self._file
        if retained is None or retained.closed:
            return -1
        return retained.fileno()

    def _add_watch(
        self, path: Path, mask: int, names: set[bytes] | None,
    ) -> None:
        require(self.descriptor >= 0,
                f"{self.label} pathname mutation guard is closed")
        try:
            function = self._library.inotify_add_watch
            function.argtypes = (
                ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32)
            function.restype = ctypes.c_int
            encoded = os.fsencode(path)
            require(b"\0" not in encoded,
                    f"{self.label} pathname contains NUL")
            ctypes.set_errno(0)
            watch = function(
                ctypes.c_int(self.descriptor), ctypes.c_char_p(encoded),
                ctypes.c_uint32(mask | IN_MASK_ADD))
            if watch < 0:
                number = ctypes.get_errno()
                raise OSError(number or errno.EPERM,
                              os.strerror(number or errno.EPERM))
        except OSError as error:
            raise BuildProvenanceError(
                f"cannot guard {self.label} pathname {path}: {error}") \
                from error
        previous = self.rules.get(watch)
        if previous is None and watch in self.rules:
            return
        if names is None:
            self.rules[watch] = None
        elif previous is None:
            self.rules[watch] = set(names)
        else:
            previous.update(names)

    @staticmethod
    def _absolute_lexical(path: Path) -> Path:
        return Path(os.path.abspath(os.fspath(path)))

    def add_file_path(self, path: Path | str) -> None:
        """Watch every pathname component plus the final regular-file inode."""
        absolute = self._absolute_lexical(Path(path))
        parts = absolute.parts
        require(absolute.is_absolute() and len(parts) >= 2,
                f"{self.label} guarded file path is invalid")
        current = Path(parts[0])
        for component in parts[1:]:
            self._add_watch(
                current,
                self._DIRECTORY_MASK | IN_ONLYDIR,
                {os.fsencode(component)})
            current = current / component
        self._add_watch(absolute, self._FILE_MASK, None)

    def add_tree(self, root: Path | str) -> None:
        """Watch all existing directories and every mutation below root."""
        absolute = self._absolute_lexical(Path(root))
        parts = absolute.parts
        require(absolute.is_absolute() and len(parts) >= 2,
                f"{self.label} guarded tree path is invalid")
        current = Path(parts[0])
        for component in parts[1:]:
            self._add_watch(
                current,
                self._DIRECTORY_MASK | IN_ONLYDIR,
                {os.fsencode(component)})
            current = current / component
        try:
            directories = [absolute]
            for directory, child_directories, _files in os.walk(
                    absolute, topdown=True, followlinks=False):
                child_directories.sort()
                parent = Path(directory)
                for name in child_directories:
                    child = parent / name
                    require(not child.is_symlink(),
                            f"{self.label} tree contains a directory symlink")
                    directories.append(child)
            for directory in directories:
                self._add_watch(
                    directory,
                    self._DIRECTORY_MASK | IN_ONLYDIR,
                    None)
        except OSError as error:
            raise BuildProvenanceError(
                f"cannot enumerate guarded {self.label} tree: {error}") \
                from error

    def add_directory_path(self, path: Path | str) -> None:
        """Watch an exact directory and all pathname components leading to it."""
        absolute = self._absolute_lexical(Path(path))
        parts = absolute.parts
        require(absolute.is_absolute() and len(parts) >= 2,
                f"{self.label} guarded directory path is invalid")
        current = Path(parts[0])
        for component in parts[1:]:
            self._add_watch(
                current,
                self._DIRECTORY_MASK | IN_ONLYDIR,
                {os.fsencode(component)})
            current = current / component
        # Child-entry changes are checked through retained openat descriptors
        # and hashes.  Here only mutation of the directory itself is relevant.
        self._add_watch(
            absolute, self._DIRECTORY_MASK | IN_ONLYDIR, set())

    def add_exact_directory_entries(
        self, directory: Path | str, names: Sequence[str],
    ) -> None:
        """Watch only topology changes for selected child entry names."""
        absolute = self._absolute_lexical(Path(directory))
        encoded = {os.fsencode(name) for name in names}
        require(absolute.is_absolute() and encoded and
                all(name and b"/" not in name and b"\0" not in name
                    for name in encoded),
                f"{self.label} exact directory entry set is invalid")
        mask = (
            IN_MOVED_FROM | IN_MOVED_TO | IN_CREATE | IN_DELETE |
            IN_DELETE_SELF | IN_MOVE_SELF | IN_UNMOUNT | IN_IGNORED |
            IN_EXCL_UNLINK | IN_ONLYDIR)
        self._add_watch(absolute, mask, encoded)

    def accept_exact_regular_creation(self, name: str) -> None:
        """Consume exactly one runner-owned IN_CREATE event.

        All exact output names are watched before any reservation is created.
        Therefore any additional event, including an adversarial create on a
        different retained name, makes this operation fail closed.
        """
        require(isinstance(name, str) and name and "/" not in name,
                f"{self.label} accepted creation name is invalid")
        self._drain()
        encoded = os.fsencode(name)
        pattern = re.compile(
            r"^wd=[0-9]+ mask=0x100 name=" +
            re.escape(repr(encoded)) + r"$")
        require(len(self.mutations) == 1 and
                pattern.fullmatch(self.mutations[0]) is not None,
                f"{self.label} reservation creation event differs: " +
                "; ".join(self.mutations[:8]))
        self.mutations.clear()

    def _drain(self) -> None:
        require(self.descriptor >= 0,
                f"{self.label} pathname mutation guard is closed")
        while True:
            try:
                data = os.read(self.descriptor, 1 << 20)
            except BlockingIOError:
                return
            except OSError as error:
                raise BuildProvenanceError(
                    f"cannot read {self.label} pathname mutation guard: "
                    f"{error}") from error
            if not data:
                raise BuildProvenanceError(
                    f"{self.label} pathname mutation guard reached EOF")
            offset = 0
            while offset < len(data):
                require(len(data) - offset >= INOTIFY_EVENT.size,
                        f"{self.label} pathname mutation event is truncated")
                watch, mask, _cookie, name_size = INOTIFY_EVENT.unpack_from(
                    data, offset)
                offset += INOTIFY_EVENT.size
                require(name_size <= len(data) - offset,
                        f"{self.label} pathname mutation name is truncated")
                name = data[offset:offset + name_size].rstrip(b"\0")
                offset += name_size
                if mask & IN_Q_OVERFLOW:
                    self.mutations.append("inotify queue overflow")
                    continue
                rule = self.rules.get(watch)
                self_event = bool(mask & (
                    IN_DELETE_SELF | IN_MOVE_SELF | IN_UNMOUNT | IN_IGNORED))
                if (watch not in self.rules or self_event or rule is None or
                        name in rule):
                    self.mutations.append(
                        f"wd={watch} mask=0x{mask:x} name={name!r}")

    def verify(self) -> None:
        self._drain()
        require(not self.mutations,
                f"{self.label} pathname changed while retained: " +
                "; ".join(self.mutations[:8]))

    def _close_without_verification(self) -> None:
        retained = self._file
        if retained is not None:
            try:
                retained.close()
            except BaseException as error:
                observed = _ordinary_owner_cleanup_error(
                    error,
                    f"cannot force-close {self.label} pathname mutation "
                    f"guard: {error}")
                if observed is error:
                    raise
                raise observed from error
            if retained.closed and self._file is retained:
                self._file = None

    def close(self) -> None:
        retained = self._file
        if retained is None:
            return
        failure: BaseException | None = None
        if self.descriptor >= 0:
            try:
                self.verify()
            except BaseException as error:
                failure = error
        try:
            retained.close()
        except BaseException as error:
            close_failure = _ordinary_owner_cleanup_error(
                error,
                f"cannot close {self.label} pathname mutation guard: "
                f"{error}")
            failure = _owner_exception_precedence(
                failure, close_failure,
                f"{self.label} pathname mutation guard close")
        if retained.closed and self._file is retained:
            self._file = None
        if failure is not None:
            raise failure

    def __enter__(self) -> "_InotifyMutationGuard":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, tb
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc is None:
                raise
            _raise_owner_exit_failure(
                exc, cleanup_error,
                f"{self.label} pathname mutation guard owner",
                f"{self.label} pathname mutation guard cleanup failed: "
                f"{cleanup_error}; primary failure: "
                f"{type(exc).__name__}: {exc}")


class _ReplayArtifactSink:
    """One parent-owned memfd populated only by the private replay child."""

    def __init__(
        self, logical_path: Path | str, label: str, *,
        maximum_bytes: int = MAX_FILE_BYTES,
    ) -> None:
        self.logical_path = Path(os.path.abspath(os.fspath(logical_path)))
        self.label = label
        self.maximum_bytes = maximum_bytes
        self._file: Any | None = None
        self.content = b""
        self.identity: dict[str, Any] = {}
        self._fields: tuple[int, ...] | None = None
        self._seals = 0
        allow_sealing = getattr(os, "MFD_ALLOW_SEALING", 0x0002)
        close_on_exec = getattr(os, "MFD_CLOEXEC", 0x0001)
        descriptor = -1
        retained = None
        try:
            descriptor = _linux_memfd_create(
                "leo2-output-" + re.sub(
                    r"[^A-Za-z0-9_.-]", "-", label)[:80],
                allow_sealing | close_on_exec)
            retained = os.fdopen(descriptor, "r+b", buffering=0)
            os.fchmod(descriptor, 0o600)
            self._file = retained
        except BaseException as primary:
            cleanup_error: BaseException | None = None
            try:
                if retained is not None:
                    retained.close()
                elif descriptor >= 0:
                    _close_raw_descriptor_with_precedence(
                        descriptor, primary,
                        f"{label} replay output raw descriptor",
                        f"cannot close interrupted {label} replay output "
                        "raw descriptor")
            except BaseException as observed:
                cleanup_error = observed
            if cleanup_error is not None:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    f"{label} replay output sink constructor",
                    f"cannot close interrupted {label} replay output sink: "
                    f"{cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            raise

    @property
    def descriptor(self) -> int:
        retained = self._file
        if retained is None or retained.closed:
            return -1
        return retained.fileno()

    @property
    def proc_path(self) -> str:
        require(self.descriptor >= 0,
                f"{self.label} replay output sink is closed")
        return f"/proc/self/fd/{self.descriptor}"

    def seal(self, *, executable: bool = False, allow_empty: bool = False) -> None:
        require(self.descriptor >= 0 and not self.identity,
                f"{self.label} replay output sink was already sealed")
        before = os.fstat(self.descriptor)
        require(
            stat.S_ISREG(before.st_mode) and
            0 <= before.st_size <= self.maximum_bytes and
            (allow_empty or before.st_size > 0),
            f"{self.label} replay output sink is empty, malformed, or too "
            "large")
        os.fchmod(self.descriptor, 0o700 if executable else 0o600)
        content = _read_bounded_descriptor(
            self.descriptor, before.st_size, self.label, self.maximum_bytes)
        seals = (
            getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
            getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
            getattr(fcntl, "F_SEAL_GROW", 0x0004) |
            getattr(fcntl, "F_SEAL_WRITE", 0x0008))
        fcntl.fcntl(
            self.descriptor, getattr(fcntl, "F_ADD_SEALS", 1033), seals)
        after = os.fstat(self.descriptor)
        observed_seals = fcntl.fcntl(
            self.descriptor, getattr(fcntl, "F_GET_SEALS", 1034))
        require(
            after.st_size == len(content) and
            observed_seals & seals == seals,
            f"{self.label} replay output sink could not be frozen")
        self.content = content
        # A denied write to a fully sealed memfd may still advance mtime/ctime
        # on Linux.  Those timestamps are therefore evidence metadata, not an
        # immutable-property check.  Bind the exact inode and all stable
        # ownership/type/size fields, then independently recheck the content
        # digest and seal set in verify().
        self._fields = (
            after.st_dev, after.st_ino, after.st_mode,
            after.st_uid, after.st_gid, after.st_nlink, after.st_size,
        )
        self._seals = observed_seals
        self.identity = {
            "path": str(self.logical_path),
            "sha256": hashlib.sha256(content).hexdigest(),
            "device": after.st_dev,
            "inode": after.st_ino,
            "mode": after.st_mode,
            "uid": after.st_uid,
            "gid": after.st_gid,
            "links": after.st_nlink,
            "size": after.st_size,
            "mtime_ns": after.st_mtime_ns,
            "ctime_ns": after.st_ctime_ns,
        }

    def verify(self) -> None:
        require(self.descriptor >= 0 and self._fields is not None and
                bool(self.identity),
                f"{self.label} replay output sink is not sealed")
        status = os.fstat(self.descriptor)
        content = _read_bounded_descriptor(
            self.descriptor, status.st_size, self.label, self.maximum_bytes)
        observed_seals = fcntl.fcntl(
            self.descriptor, getattr(fcntl, "F_GET_SEALS", 1034))
        require(
            (
                status.st_dev, status.st_ino, status.st_mode,
                status.st_uid, status.st_gid, status.st_nlink,
                status.st_size,
            ) == self._fields and
            observed_seals == self._seals and
            hashlib.sha256(content).hexdigest() ==
                self.identity["sha256"] and
            content == self.content,
            f"{self.label} sealed replay output changed")

    def close(self) -> None:
        retained = self._file
        if retained is None:
            return
        failure: BaseException | None = None
        if self.identity:
            try:
                self.verify()
            except BaseException as error:
                failure = error
        try:
            retained.close()
        except BaseException as error:
            failure = _owner_exception_precedence(
                failure, error, f"{self.label} replay output sink close")
        if retained.closed and self._file is retained:
            self._file = None
        if failure is not None:
            raise failure

    def __enter__(self) -> "_ReplayArtifactSink":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, tb
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc is None:
                raise
            _raise_owner_exit_failure(
                exc, cleanup_error, f"sealed {self.label} owner",
                f"sealed {self.label} cleanup failed: {cleanup_error}; "
                f"primary failure: {type(exc).__name__}: {exc}")


class _RetainedFileSnapshot:
    """Hold the exact regular file whose bytes and pathname were attested.

    Keeping this descriptor open lets callers execute or reopen the attested
    object through /proc/self/fd without a pathname ABA window.  The lexical
    command remains argv[0]; _run() uses the descriptor only as execve's
    executable pathname.
    """

    def __init__(
        self, path: Path | str, label: str, *,
        maximum_bytes: int = MAX_FILE_BYTES,
    ) -> None:
        require(type(maximum_bytes) is int and maximum_bytes >= 0,
                f"{label} retained byte bound is invalid")
        self.requested = Path(os.path.abspath(os.fspath(path)))
        self.label = label
        self.maximum_bytes = maximum_bytes
        self.resolved: Path | None = None
        self._source_file: Any | None = None
        self._executable_file: Any | None = None
        self.identity: dict[str, Any] = {}
        self.executable_identity: dict[str, Any] = {}
        self.content = b""
        self._fields: tuple[int, ...] | None = None
        self._path_guard: _InotifyMutationGuard | None = None
        try:
            expected_resolved = self.requested.resolve(strict=True)
            self._path_guard = _InotifyMutationGuard(label)
            self._path_guard.add_file_path(self.requested)
            if (self._path_guard._absolute_lexical(self.requested) !=
                    expected_resolved):
                self._path_guard.add_file_path(expected_resolved)
            self._open(expected_resolved)
            self._path_guard.verify()
        except BaseException as primary:
            try:
                self._close_without_verification()
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    f"retained {label} constructor",
                    f"retained {label} constructor cleanup failed: "
                    f"{cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            raise

    def _verify_path(self, status: os.stat_result, action: str) -> None:
        require(self.resolved is not None,
                f"{self.label} retained pathname was lost")
        try:
            path_status = self.resolved.lstat()
            requested_status = self.requested.stat()
            resolved_after = self.requested.resolve(strict=True)
        except (OSError, RuntimeError) as error:
            raise BuildProvenanceError(
                f"{self.label} pathname changed {action}: {error}") from error
        require(_stable_fields(status) == _stable_fields(path_status) and
                _stable_fields(status) == _stable_fields(requested_status) and
                resolved_after == self.resolved,
                f"{self.label} pathname changed {action}")

    def _open(self, expected_resolved: Path) -> None:
        descriptor = -1
        source_file = None
        try:
            requested_before = self.requested.stat()
            require(stat.S_ISREG(requested_before.st_mode),
                    f"{self.label} is not a regular file")
            self.resolved = self.requested.resolve(strict=True)
            require(self.resolved == expected_resolved,
                    f"{self.label} pathname changed before it was opened")
            path_before = self.resolved.lstat()
            require(
                _stable_fields(requested_before) == _stable_fields(path_before),
                f"{self.label} pathname changed before it was opened")
            descriptor = os.open(
                self.resolved,
                os.O_RDONLY |
                getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_NONBLOCK", 0) |
                getattr(os, "O_NOFOLLOW", 0),
            )
            source_file = os.fdopen(descriptor, "rb", buffering=0)
            descriptor = source_file.fileno()
            before = os.fstat(descriptor)
            require(stat.S_ISREG(before.st_mode),
                    f"{self.label} is not a regular file")
            require(_stable_fields(path_before) == _stable_fields(before),
                    f"{self.label} pathname changed while it was opened")
            require(0 <= before.st_size <= self.maximum_bytes,
                    f"{self.label} exceeds its retained byte bound")
            chunks: list[bytes] = []
            remaining = before.st_size
            while remaining:
                block = os.read(
                    descriptor, min(1 << 20, remaining))
                require(bool(block),
                        f"{self.label} ended before its recorded size")
                chunks.append(block)
                remaining -= len(block)
            require(not os.read(descriptor, 1),
                    f"{self.label} grew beyond its recorded size while read")
            after = os.fstat(descriptor)
            require(_stable_fields(before) == _stable_fields(after),
                    f"{self.label} changed while read")
            self._verify_path(after, "while read")
            content = b"".join(chunks)
            fields = _stable_fields(after)
            identity = {
                "path": str(self.resolved),
                "sha256": hashlib.sha256(content).hexdigest(),
                "device": after.st_dev,
                "inode": after.st_ino,
                "mode": after.st_mode,
                "uid": after.st_uid,
                "gid": after.st_gid,
                "links": after.st_nlink,
                "size": after.st_size,
                "mtime_ns": after.st_mtime_ns,
                "ctime_ns": after.st_ctime_ns,
            }
            self.content = content
            self._fields = fields
            self.identity = identity
            self._source_file = source_file
        except BaseException as primary:
            try:
                if source_file is None and descriptor >= 0:
                    _close_raw_descriptor_with_precedence(
                        descriptor, primary,
                        f"retained {self.label} raw descriptor",
                        f"cannot close interrupted retained {self.label} "
                        "raw descriptor")
                elif source_file is not None:
                    source_file.close()
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    f"retained {self.label} open",
                    f"cannot close interrupted retained {self.label}: "
                    f"{cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            if self._source_file is source_file:
                self._source_file = None
            try:
                self._close_without_verification()
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    f"retained {self.label} open cleanup",
                    f"cannot finish interrupted retained {self.label} "
                    f"cleanup: {cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            if isinstance(primary, BuildProvenanceError):
                raise
            if not isinstance(primary, (OSError, RuntimeError)):
                raise
            raise BuildProvenanceError(
                f"cannot open stable {self.label}: {primary}") from primary

    @property
    def descriptor(self) -> int:
        source = self._source_file
        if source is None or source.closed:
            return -1
        return source.fileno()

    @property
    def proc_path(self) -> str:
        require(self.descriptor >= 0,
                f"{self.label} retained descriptor is closed")
        return f"/proc/self/fd/{self.descriptor}"

    @property
    def _executable_descriptor(self) -> int:
        executable = self._executable_file
        if executable is None or executable.closed:
            return -1
        return executable.fileno()

    @property
    def executable_descriptor(self) -> int:
        """Return a sealed immutable copy suitable for descriptor exec."""
        require(self.descriptor >= 0,
                f"{self.label} retained descriptor is closed")
        if self._executable_descriptor >= 0:
            self._verify_executable_descriptor()
            return self._executable_descriptor
        allow_sealing = getattr(os, "MFD_ALLOW_SEALING", 0x0002)
        close_on_exec = getattr(os, "MFD_CLOEXEC", 0x0001)
        descriptor = -1
        executable = None
        try:
            descriptor = _linux_memfd_create(
                "leo2-" + re.sub(
                    r"[^A-Za-z0-9_.-]", "-", self.label)[:80],
                allow_sealing | close_on_exec)
            # os.fdopen() owns the raw descriptor only after it returns
            # successfully.  Keep the raw integer under this cleanup boundary
            # so an injected BaseException from fdopen cannot leak the memfd.
            executable = os.fdopen(
                descriptor, "r+b", buffering=0)
            descriptor = executable.fileno()
            view = memoryview(self.content)
            while view:
                written = os.write(descriptor, view)
                require(written > 0,
                        f"{self.label} immutable executable write stalled")
                view = view[written:]
            os.fchmod(descriptor, stat.S_IMODE(self.identity["mode"]))
            seals = (
                getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
                getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
                getattr(fcntl, "F_SEAL_GROW", 0x0004) |
                getattr(fcntl, "F_SEAL_WRITE", 0x0008))
            fcntl.fcntl(
                descriptor, getattr(fcntl, "F_ADD_SEALS", 1033), seals)
            observed_seals = fcntl.fcntl(
                descriptor, getattr(fcntl, "F_GET_SEALS", 1034))
            require(observed_seals & seals == seals,
                    f"{self.label} immutable executable is not fully sealed")
            executable_identity = self._executable_descriptor_identity(
                descriptor)

            # Keep ownership local until every operation that can fail has
            # completed.  Publishing the descriptor earlier lets an
            # asynchronous BaseException close the local fd while leaving the
            # object claiming ownership of its number; if that number is then
            # recycled, close() can close an unrelated file.
            self.executable_identity = executable_identity
            self._executable_file = executable
            return self._executable_descriptor
        except BaseException as primary:
            try:
                if executable is None and descriptor >= 0:
                    _close_raw_descriptor_with_precedence(
                        descriptor, primary,
                        f"immutable {self.label} executable raw descriptor",
                        f"cannot close interrupted immutable {self.label} "
                        "executable raw descriptor")
                elif executable is not None:
                    executable.close()
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    f"immutable {self.label} executable constructor",
                    f"cannot close interrupted immutable {self.label} "
                    f"executable: {cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            if self._executable_file is executable:
                self._executable_file = None
            self.executable_identity = {}
            raise

    def _executable_descriptor_identity(
            self, descriptor: int) -> dict[str, Any]:
        require(type(descriptor) is int and descriptor >= 0,
                f"{self.label} immutable executable is closed")
        status = os.fstat(descriptor)
        seals = (
            getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
            getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
            getattr(fcntl, "F_SEAL_GROW", 0x0004) |
            getattr(fcntl, "F_SEAL_WRITE", 0x0008))
        observed_seals = fcntl.fcntl(
            descriptor,
            getattr(fcntl, "F_GET_SEALS", 1034))
        sealed_content = _read_bounded_descriptor(
            descriptor, status.st_size,
            f"{self.label} immutable executable", self.maximum_bytes)
        sealed_digest = hashlib.sha256(sealed_content).hexdigest()
        require(stat.S_ISREG(status.st_mode) and
                status.st_size == len(self.content) and
                observed_seals & seals == seals and
                sealed_digest == self.identity["sha256"],
                f"{self.label} immutable executable changed")
        return {
            "sha256": sealed_digest,
            "size": status.st_size,
            "seals": observed_seals,
            "source_sha256": self.identity["sha256"],
        }

    def _verify_executable_descriptor(self) -> None:
        require(self._executable_descriptor >= 0,
                f"{self.label} immutable executable is closed")
        require(
            self._executable_descriptor_identity(
                self._executable_descriptor) == self.executable_identity,
            f"{self.label} immutable executable identity changed")

    def verify(self) -> None:
        require(self.descriptor >= 0 and self._fields is not None,
                f"{self.label} retained descriptor is closed")
        try:
            current = os.fstat(self.descriptor)
        except OSError as error:
            raise BuildProvenanceError(
                f"cannot verify retained {self.label}: {error}") from error
        require(_stable_fields(current) == self._fields,
                f"{self.label} changed while retained")
        self._verify_path(current, "while retained")
        if self._executable_descriptor >= 0:
            self._verify_executable_descriptor()

    def _close_without_verification(self) -> None:
        failure: BaseException | None = None
        executable = self._executable_file
        if executable is not None:
            try:
                executable.close()
            except BaseException as error:
                observed = _ordinary_owner_cleanup_error(
                    error,
                    f"cannot force-close immutable {self.label} "
                    f"executable: {error}")
                failure = _owner_exception_precedence(
                    failure, observed,
                    f"retained {self.label} forced executable close")
            if executable.closed and self._executable_file is executable:
                self._executable_file = None
        source = self._source_file
        if source is not None:
            try:
                source.close()
            except BaseException as error:
                observed = _ordinary_owner_cleanup_error(
                    error,
                    f"cannot force-close retained {self.label}: {error}")
                failure = _owner_exception_precedence(
                    failure, observed,
                    f"retained {self.label} forced source close")
            if source.closed and self._source_file is source:
                self._source_file = None
        guard = self._path_guard
        if guard is not None:
            try:
                guard._close_without_verification()
            except BaseException as error:
                failure = _owner_exception_precedence(
                    failure, error,
                    f"retained {self.label} forced guard close")
            if guard.descriptor < 0 and self._path_guard is guard:
                self._path_guard = None
        if failure is not None:
            raise failure.with_traceback(failure.__traceback__)

    def close(self) -> None:
        if self._source_file is None and self._executable_file is None and \
                self._path_guard is None:
            return
        verification_error: BaseException | None = None
        if self.descriptor >= 0:
            try:
                self.verify()
            except BaseException as error:
                verification_error = error
        guard = self._path_guard
        if guard is not None:
            try:
                guard.close()
            except BaseException as error:
                verification_error = _owner_exception_precedence(
                    verification_error, error,
                    f"retained {self.label} guard close")
            if guard.descriptor < 0 and self._path_guard is guard:
                self._path_guard = None
        source = self._source_file
        if source is not None:
            try:
                source.close()
            except BaseException as error:
                close_failure = _ordinary_owner_cleanup_error(
                    error, f"cannot close retained {self.label}: {error}")
                verification_error = _owner_exception_precedence(
                    verification_error, close_failure,
                    f"retained {self.label} source close")
            if source.closed and self._source_file is source:
                self._source_file = None
        executable = self._executable_file
        if executable is not None:
            try:
                executable.close()
            except BaseException as error:
                close_failure = _ordinary_owner_cleanup_error(
                    error,
                    f"cannot close immutable {self.label} executable: "
                    f"{error}")
                verification_error = _owner_exception_precedence(
                    verification_error, close_failure,
                    f"retained {self.label} executable close")
            if executable.closed and self._executable_file is executable:
                self._executable_file = None
        if verification_error is not None:
            raise verification_error

    def __enter__(self) -> "_RetainedFileSnapshot":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, tb
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc is None:
                raise
            _raise_owner_exit_failure(
                exc, cleanup_error, f"retained {self.label} owner",
                f"retained {self.label} cleanup failed: {cleanup_error}; "
                f"primary failure: {type(exc).__name__}: {exc}")


class _RetainedDirectoryTree:
    """Hold one exact directory inode and every mutation below it."""

    def __init__(self, path: Path | str, label: str) -> None:
        self.requested = Path(os.path.abspath(os.fspath(path)))
        self.label = label
        self.resolved: Path | None = None
        self._owner = _OwnedDescriptor()
        self._fields: tuple[int, ...] | None = None
        self.guard: _InotifyMutationGuard | None = None
        try:
            expected = self.requested.resolve(strict=True)
            require(expected.is_dir(), f"{label} is not a directory")
            self.guard = _InotifyMutationGuard(label)
            self.guard.add_directory_path(self.requested)
            if self.guard._absolute_lexical(self.requested) != expected:
                self.guard.add_directory_path(expected)
            self.guard.add_tree(expected)
            self._owner.open(
                expected,
                os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_NOFOLLOW", 0))
            status = os.fstat(self.descriptor)
            require(stat.S_ISDIR(status.st_mode),
                    f"{label} retained object is not a directory")
            self.resolved = expected
            self._fields = _stable_fields(status)
            self.verify()
        except BaseException as primary:
            try:
                self._close_without_verification()
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    f"retained {label} directory constructor",
                    f"retained {label} directory constructor cleanup "
                    f"failed: {cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            raise

    @property
    def descriptor(self) -> int:
        return self._owner.descriptor

    def verify(self) -> None:
        require(self.descriptor >= 0 and self.resolved is not None and
                self._fields is not None and self.guard is not None,
                f"{self.label} retained directory is closed")
        self.guard.verify()
        try:
            current = os.fstat(self.descriptor)
            requested = self.requested.stat()
            resolved_status = self.resolved.lstat()
            resolved_after = self.requested.resolve(strict=True)
        except (OSError, RuntimeError) as error:
            raise BuildProvenanceError(
                f"{self.label} directory pathname changed while retained: "
                f"{error}") from error
        require(_stable_fields(current) == self._fields and
                _stable_fields(requested) == self._fields and
                _stable_fields(resolved_status) == self._fields and
                resolved_after == self.resolved,
                f"{self.label} directory pathname changed while retained")
        self.guard.verify()

    def _close_without_verification(self) -> None:
        failure: BaseException | None = None
        try:
            self._owner.close()
        except BaseException as error:
            failure = _ordinary_owner_cleanup_error(
                error,
                f"cannot force-close retained {self.label} directory: "
                f"{error}")
        guard = self.guard
        if guard is not None:
            try:
                guard._close_without_verification()
            except BaseException as error:
                failure = _owner_exception_precedence(
                    failure, error,
                    f"retained {self.label} forced directory guard close")
            if guard.descriptor < 0 and self.guard is guard:
                self.guard = None
        if failure is not None:
            raise failure.with_traceback(failure.__traceback__)

    def close(self) -> None:
        if self.descriptor < 0 and self.guard is None:
            return
        failure: BaseException | None = None
        if self.descriptor >= 0:
            try:
                self.verify()
            except BaseException as error:
                failure = error
        guard = self.guard
        if guard is not None:
            try:
                guard.close()
            except BaseException as error:
                failure = _owner_exception_precedence(
                    failure, error,
                    f"retained {self.label} directory guard close")
            if guard.descriptor < 0 and self.guard is guard:
                self.guard = None
        if self.descriptor >= 0:
            try:
                self._owner.close()
            except BaseException as error:
                close_failure = _ordinary_owner_cleanup_error(
                    error,
                    f"cannot close retained {self.label} directory: {error}")
                failure = _owner_exception_precedence(
                    failure, close_failure,
                    f"retained {self.label} directory descriptor close")
        if failure is not None:
            raise failure

    def __enter__(self) -> "_RetainedDirectoryTree":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, tb
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc is None:
                raise
            _raise_owner_exit_failure(
                exc, cleanup_error,
                f"retained {self.label} directory owner",
                f"retained {self.label} directory cleanup failed: "
                f"{cleanup_error}; primary failure: "
                f"{type(exc).__name__}: {exc}")


def _retain_exact_symlink_directory(
    path: Path, expected_targets: Mapping[str, str], label: str,
) -> _RetainedDirectoryTree:
    """Retain a private tool prefix and bind every role to its exact target.

    The directory guard is armed before this verification.  Enumerating through
    its retained descriptor then closes the pre-guard window where a same-UID
    process could replace one symlink after its initial readlink check but
    before mutation monitoring began.
    """
    require(isinstance(expected_targets, Mapping) and expected_targets and
            all(isinstance(role, str) and
                re.fullmatch(r"[A-Za-z0-9_.+-]+", role) is not None and
                isinstance(target, str) and target.startswith(
                    "/proc/self/fd/") and "\0" not in target
                for role, target in expected_targets.items()),
            f"{label} expected mappings are invalid")
    guard = _RetainedDirectoryTree(path, label)
    try:
        observed: dict[str, str] = {}
        with os.scandir(guard.descriptor) as entries:
            for entry in entries:
                role = entry.name
                require(role not in observed and entry.is_symlink(),
                        f"{label} contains a duplicate or non-symlink entry")
                target = os.readlink(role, dir_fd=guard.descriptor)
                _require_safe_unicode(target, f"{label} mapping target")
                observed[role] = target
        require(observed == dict(expected_targets),
                f"{label} mappings changed before retention")
        guard.verify()
        return guard
    except BaseException as primary:
        try:
            guard._close_without_verification()
        except BaseException as cleanup_error:
            _raise_owner_exit_failure(
                primary, cleanup_error,
                f"{label} exact symlink-directory retention",
                f"{label} exact symlink-directory cleanup failed: "
                f"{cleanup_error}; primary failure: "
                f"{type(primary).__name__}: {primary}")
        raise


def _parse_git_directory_reference(
    content: bytes, base: Path, label: str, *, prefix: bytes = b"",
) -> Path:
    require(0 < len(content) <= MAX_GITFILE_BYTES,
            f"{label} is empty or oversized")
    require(content.endswith(b"\n") and content.count(b"\n") == 1,
            f"{label} is not one canonical LF-terminated record")
    record = content[:-1]
    require(record.startswith(prefix) and len(record) > len(prefix),
            f"{label} has an invalid record prefix")
    try:
        value = record[len(prefix):].decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise BuildProvenanceError(
            f"{label} path is not strict UTF-8") from error
    _require_safe_unicode(value, f"{label} path")
    require("\0" not in value and len(os.fsencode(value)) <=
            MAX_TRACKED_SOURCE_PATH_BYTES,
            f"{label} path is invalid or oversized")
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    try:
        return path.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise BuildProvenanceError(
            f"{label} path cannot be resolved: {error}") from error


class _RetainedGitMetadata:
    """Retain the worktree .git entry, exact gitdir, and common metadata."""

    def __init__(self, source: Path) -> None:
        self.source = source
        self.stack = ExitStack()
        self.entry_file: _RetainedFileSnapshot | None = None
        self.gitdir: _RetainedDirectoryTree | None = None
        self.common: _RetainedDirectoryTree | None = None
        self.common_file: _RetainedFileSnapshot | None = None
        try:
            entry = source / ".git"
            entry_status = entry.lstat()
            if stat.S_ISDIR(entry_status.st_mode):
                gitdir_path = entry.resolve(strict=True)
            elif stat.S_ISREG(entry_status.st_mode):
                self.entry_file = self.stack.enter_context(
                    _RetainedFileSnapshot(
                        entry, "tracked source .git file",
                        maximum_bytes=MAX_GITFILE_BYTES))
                gitdir_path = _parse_git_directory_reference(
                    self.entry_file.content, entry.parent,
                    "tracked source .git file", prefix=b"gitdir: ")
            else:
                raise BuildProvenanceError(
                    "tracked source .git entry is not a directory or "
                    "regular gitfile")

            self.gitdir = self.stack.enter_context(
                _RetainedDirectoryTree(
                    gitdir_path, "tracked source Git directory"))
            commondir_path = self.gitdir.resolved / "commondir"
            try:
                commondir_status = commondir_path.lstat()
            except FileNotFoundError:
                commondir_status = None
            if commondir_status is None:
                self.common = self.gitdir
            else:
                require(stat.S_ISREG(commondir_status.st_mode),
                        "tracked source Git commondir is not a regular file")
                self.common_file = self.stack.enter_context(
                    _RetainedFileSnapshot(
                        commondir_path, "tracked source Git commondir",
                        maximum_bytes=MAX_GITFILE_BYTES))
                common_path = _parse_git_directory_reference(
                    self.common_file.content, self.gitdir.resolved,
                    "tracked source Git commondir")
                if common_path == self.gitdir.resolved:
                    self.common = self.gitdir
                else:
                    self.common = self.stack.enter_context(
                        _RetainedDirectoryTree(
                            common_path,
                            "tracked source common Git directory"))

            require(self.common is not None and
                    self.common.resolved is not None,
                    "tracked source common Git directory was lost")
            for root in {self.gitdir.resolved, self.common.resolved}:
                alternates = root / "objects" / "info" / "alternates"
                require(not alternates.exists() and not alternates.is_symlink(),
                        "tracked source Git object alternates are unsupported")
            self.verify()
        except BaseException as primary:
            try:
                self.stack.close()
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    "retained Git metadata constructor",
                    "retained Git metadata constructor cleanup failed: "
                    f"{cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            raise

    @property
    def descriptor(self) -> int:
        require(self.gitdir is not None and self.gitdir.descriptor >= 0,
                "tracked source Git directory descriptor is closed")
        return self.gitdir.descriptor

    def verify(self) -> None:
        if self.entry_file is not None:
            self.entry_file.verify()
        require(self.gitdir is not None,
                "tracked source Git directory was lost")
        self.gitdir.verify()
        if self.common_file is not None:
            self.common_file.verify()
        if self.common is not None and self.common is not self.gitdir:
            self.common.verify()

    def close(self) -> None:
        failure: BaseException | None = None
        try:
            self.verify()
        except BaseException as error:
            failure = error
        try:
            self.stack.close()
        except BaseException as error:
            failure = _owner_exception_precedence(
                failure, error, "retained Git metadata stack close")
        if failure is not None:
            raise failure

    def __enter__(self) -> "_RetainedGitMetadata":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, tb
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc is None:
                raise
            _raise_owner_exit_failure(
                exc, cleanup_error, "retained Git metadata owner",
                "retained Git metadata cleanup failed: "
                f"{cleanup_error}; primary failure: "
                f"{type(exc).__name__}: {exc}")


def file_snapshot(
    path: Path | str, label: str, *, maximum_bytes: int = MAX_FILE_BYTES,
) -> tuple[dict[str, Any], bytes]:
    """Read one stable regular-file snapshot and close its retained handle."""
    with _RetainedFileSnapshot(
            path, label, maximum_bytes=maximum_bytes) as snapshot:
        return snapshot.identity, snapshot.content


def file_identity(
    path: Path | str, label: str, *, maximum_bytes: int = MAX_FILE_BYTES,
) -> dict[str, Any]:
    return file_snapshot(path, label, maximum_bytes=maximum_bytes)[0]


def _read_bounded_descriptor(
    descriptor: int, size: int, label: str, maximum_bytes: int,
) -> bytes:
    require(0 <= size <= maximum_bytes,
            f"{label} exceeds its retained byte bound")
    try:
        os.lseek(descriptor, 0, os.SEEK_SET)
        chunks: list[bytes] = []
        remaining = size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            require(bool(block), f"{label} ended before its recorded size")
            chunks.append(block)
            remaining -= len(block)
        require(not os.read(descriptor, 1),
                f"{label} grew beyond its recorded size while read")
        return b"".join(chunks)
    except OSError as error:
        raise BuildProvenanceError(
            f"cannot read stable {label}: {error}") from error


class _OpenDirectoryTree:
    """Traverse one retained directory without following component symlinks."""

    def __init__(self, root: Path, stack: ExitStack, label: str) -> None:
        self.root = root
        self.stack = stack
        self.label = label
        self._owners: list[_OwnedDescriptor] = []
        owner = _OwnedDescriptor()
        self._owners.append(owner)
        try:
            descriptor = owner.open(
                root,
                os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_NOFOLLOW", 0))
            stack.enter_context(owner)
        except BaseException as error:
            _close_owned_descriptor_after_failure(
                owner, error, f"retained {label} root constructor")
            if not isinstance(error, OSError):
                raise
            raise BuildProvenanceError(
                f"cannot open retained {label} root: {error}") from error
        self.root_descriptor = descriptor
        self.directories: dict[tuple[str, ...], int] = {(): descriptor}

    def directory(self, parts: tuple[str, ...]) -> int:
        current_parts: tuple[str, ...] = ()
        current = self.root_descriptor
        for component in parts:
            current_parts += (component,)
            existing = self.directories.get(current_parts)
            if existing is not None:
                current = existing
                continue
            owner = _OwnedDescriptor()
            self._owners.append(owner)
            try:
                descriptor = owner.open(
                    component,
                    os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_DIRECTORY", 0) |
                    getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=current)
                self.stack.enter_context(owner)
            except BaseException as error:
                _close_owned_descriptor_after_failure(
                    owner, error,
                    f"{self.label} directory component constructor")
                if not isinstance(error, OSError):
                    raise
                raise BuildProvenanceError(
                    f"{self.label} directory component "
                    f"{'/'.join(current_parts)!r} is unstable or unsafe: "
                    f"{error}") from error
            self.directories[current_parts] = descriptor
            current = descriptor
        return current

    def open_regular(self, relative: str) -> tuple[int, os.stat_result]:
        parts = tuple(relative.split("/"))
        require(parts and all(
                    component not in ("", ".", "..")
                    for component in parts),
                f"{self.label} contains unsafe relative path {relative!r}")
        parent = self.directory(parts[:-1])
        owner = _OwnedDescriptor()
        descriptor = -1
        try:
            descriptor = owner.open(
                parts[-1],
                os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_NONBLOCK", 0) |
                getattr(os, "O_NOFOLLOW", 0),
                dir_fd=parent)
            status = os.fstat(descriptor)
        except BaseException as error:
            _close_owned_descriptor_after_failure(
                owner, error,
                f"{self.label} regular-file descriptor constructor")
            if not isinstance(error, OSError):
                raise
            raise BuildProvenanceError(
                f"cannot open stable {self.label} file {relative!r}: "
                f"{error}") from error
        if not stat.S_ISREG(status.st_mode):
            primary = BuildProvenanceError(
                f"{self.label} file {relative!r} is not regular")
            _close_owned_descriptor_after_failure(
                owner, primary,
                f"{self.label} non-regular descriptor cleanup")
            raise primary
        try:
            self.stack.enter_context(owner)
        except BaseException as primary:
            _close_owned_descriptor_after_failure(
                owner, primary,
                f"{self.label} descriptor callback registration")
            raise
        self._owners.append(owner)
        return descriptor, status


class _CreateDirectoryTree:
    """Create files below one retained private root using openat only."""

    def __init__(self, root: Path, stack: ExitStack, label: str) -> None:
        self.root = root
        self.stack = stack
        self.label = label
        self._owners: list[_OwnedDescriptor] = []
        owner = _OwnedDescriptor()
        self._owners.append(owner)
        try:
            root.mkdir(mode=0o700)
            descriptor = owner.open(
                root,
                os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_NOFOLLOW", 0))
            stack.enter_context(owner)
        except BaseException as error:
            _close_owned_descriptor_after_failure(
                owner, error, f"retained {label} root creation")
            if not isinstance(error, OSError):
                raise
            raise BuildProvenanceError(
                f"cannot create retained {label} root: {error}") from error
        self.root_descriptor = descriptor
        self.directories: dict[tuple[str, ...], int] = {(): descriptor}

    def directory(self, parts: tuple[str, ...]) -> int:
        current_parts: tuple[str, ...] = ()
        current = self.root_descriptor
        for component in parts:
            current_parts += (component,)
            existing = self.directories.get(current_parts)
            if existing is not None:
                current = existing
                continue
            try:
                os.mkdir(component, mode=0o700, dir_fd=current)
            except FileExistsError:
                pass
            except OSError as error:
                raise BuildProvenanceError(
                    f"cannot create {self.label} directory "
                    f"{'/'.join(current_parts)!r}: {error}") from error
            owner = _OwnedDescriptor()
            self._owners.append(owner)
            try:
                descriptor = owner.open(
                    component,
                    os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_DIRECTORY", 0) |
                    getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=current)
                self.stack.enter_context(owner)
            except BaseException as error:
                _close_owned_descriptor_after_failure(
                    owner, error,
                    f"{self.label} created-directory descriptor")
                if not isinstance(error, OSError):
                    raise
                raise BuildProvenanceError(
                    f"{self.label} directory component "
                    f"{'/'.join(current_parts)!r} is unstable or unsafe: "
                    f"{error}") from error
            self.directories[current_parts] = descriptor
            current = descriptor
        return current

    def write_exclusive(
        self, relative: str, content: bytes, mode: int,
    ) -> None:
        parts = tuple(relative.split("/"))
        require(parts and all(
                    component not in ("", ".", "..")
                    for component in parts),
                f"{self.label} contains unsafe relative path {relative!r}")
        parent = self.directory(parts[:-1])
        flags = (
            os.O_WRONLY | os.O_CREAT | os.O_EXCL |
            getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0))
        owner = _OwnedDescriptor()
        primary: BaseException | None = None
        try:
            descriptor = owner.open(
                parts[-1], flags, mode | 0o400, dir_fd=parent)
            os.fchmod(descriptor, mode)
            view = memoryview(content)
            while view:
                written = os.write(descriptor, view)
                require(written > 0,
                        f"{self.label} write stalled")
                view = view[written:]
            os.fsync(descriptor)
        except BaseException as error:
            primary = error
        try:
            owner.close()
        except BaseException as cleanup_error:
            if primary is not None:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    f"{self.label} exclusive-file owner",
                    f"cannot close stable {self.label} file {relative!r}: "
                    f"{cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            if not isinstance(cleanup_error, OSError):
                raise
            raise BuildProvenanceError(
                f"cannot write stable {self.label} file {relative!r}: "
                f"{cleanup_error}") from cleanup_error
        if primary is None:
            return
        if not isinstance(primary, OSError):
            raise primary.with_traceback(primary.__traceback__)
        raise BuildProvenanceError(
            f"cannot write stable {self.label} file {relative!r}: "
            f"{primary}") from primary


def _git_source_state(
    source_descriptor: int, inherited_descriptors: Sequence[int],
    git_snapshot: _RetainedFileSnapshot,
    git_metadata: _RetainedGitMetadata,
) -> tuple[list[str], dict[str, Any]]:
    require(git_snapshot.descriptor >= 0 and git_snapshot.resolved is not None,
            "source-inventory Git executable is not retained")
    git_snapshot.verify()
    git_metadata.verify()
    inherited = tuple(inherited_descriptors) + (
        source_descriptor, git_metadata.descriptor)
    source_argument = f"/proc/self/fd/{source_descriptor}"

    def invoke(
        arguments: Sequence[str], label: str,
        maximum_bytes: int = MAX_METADATA_BYTES,
    ) -> bytes:
        return _run(
            (
                str(git_snapshot.resolved),
                "-C", source_argument,
                f"--git-dir=/proc/self/fd/{git_metadata.descriptor}",
                f"--work-tree={source_argument}",
                *arguments,
            ),
            label, maximum_bytes=maximum_bytes,
            inherited_descriptors=inherited,
            executable_descriptor=git_snapshot.executable_descriptor)

    raw_paths = invoke(
        ("ls-files", "--recurse-submodules", "-z"),
        "tracked source inventory", 32 << 20)
    commit = invoke(
        ("rev-parse", "--verify", "HEAD"),
        "tracked source commit").strip()
    tree = invoke(
        ("rev-parse", "--verify", "HEAD^{tree}"),
        "tracked source tree").strip()
    status_bytes = invoke(
        ("status", "--porcelain=v1", "--untracked-files=normal",
         "--ignore-submodules=none"),
        "tracked source status", 32 << 20).rstrip()
    git_metadata.verify()
    git_snapshot.verify()

    paths: list[str] = []
    for encoded in raw_paths.split(b"\0"):
        if not encoded:
            continue
        require(len(encoded) <= MAX_TRACKED_SOURCE_PATH_BYTES,
                "tracked source path exceeds its byte bound")
        try:
            relative = encoded.decode("utf-8", errors="strict")
        except UnicodeDecodeError as error:
            raise BuildProvenanceError(
                "tracked source path is not strict UTF-8") from error
        _require_safe_unicode(relative, "tracked source path")
        require(not relative.startswith("/") and
                all(component not in ("", ".", "..")
                    for component in relative.split("/")),
                f"tracked source inventory contains unsafe path {relative!r}")
        paths.append(relative)
    require(paths and len(paths) <= MAX_TRACKED_SOURCE_FILES and
            len(paths) == len(set(paths)),
            "tracked source inventory is empty, oversized, or duplicated")
    require(paths == sorted(paths),
            "tracked source inventory is not canonical")
    try:
        commit_text = commit.decode("ascii")
        tree_text = tree.decode("ascii")
    except UnicodeDecodeError as error:
        raise BuildProvenanceError(
            "tracked source Git identity is not ASCII") from error
    require(re.fullmatch(r"[0-9a-f]{40,64}", commit_text) is not None and
            re.fullmatch(r"[0-9a-f]{40,64}", tree_text) is not None,
            "tracked source Git identity is non-canonical")
    return paths, {
        "commit": commit_text,
        "tree": tree_text,
        "dirty": bool(status_bytes),
        "status_sha256": hashlib.sha256(status_bytes).hexdigest(),
    }


def _capture_tracked_source_tree(
    source: Path, destination: Path | None = None,
    *, inherited_descriptors: Sequence[int] = (),
) -> dict[str, Any]:
    """Capture one simultaneous, bounded tracked-tree state.

    Every source descriptor remains open until every file has been read and
    revalidated by content and by root-relative pathname.  Consequently a
    replacement restored before the final check cannot create a mixed tree.
    """
    source = source.resolve(strict=True)
    require(source.is_dir(), "tracked source root is not a directory")
    path_guard = _InotifyMutationGuard("tracked source root")
    records: list[dict[str, Any]] = []
    total = 0
    try:
        path_guard.add_directory_path(source)
        with ExitStack() as stack:
            tree = _OpenDirectoryTree(source, stack, "tracked source")
            git_snapshot = stack.enter_context(
                _RetainedFileSnapshot(
                    Path("/usr/bin/git").resolve(strict=True),
                    "source-inventory Git executable"))
            git_metadata = stack.enter_context(
                _RetainedGitMetadata(source))
            paths_before, git_before = _git_source_state(
                tree.root_descriptor, inherited_descriptors,
                git_snapshot, git_metadata)
            retained: list[
                tuple[str, int, tuple[int, ...], str, int]
            ] = []
            output_tree: _CreateDirectoryTree | None = None
            if destination is not None:
                require(not destination.exists(),
                        "tracked source snapshot destination already exists")
                output_tree = _CreateDirectoryTree(
                    destination, stack, "tracked source snapshot")
            for relative in paths_before:
                descriptor, status = tree.open_regular(relative)
                require(status.st_size <= MAX_TRACKED_SOURCE_FILE_BYTES,
                        f"tracked source file {relative!r} exceeds its "
                        "individual byte bound")
                total += status.st_size
                require(total <= MAX_TRACKED_SOURCE_TOTAL_BYTES,
                        "tracked source tree exceeds its total byte bound")
                content = _read_bounded_descriptor(
                    descriptor, status.st_size,
                    f"tracked source file {relative!r}",
                    MAX_TRACKED_SOURCE_FILE_BYTES)
                after = os.fstat(descriptor)
                require(_stable_fields(status) == _stable_fields(after),
                        f"tracked source file {relative!r} changed while read")
                digest = hashlib.sha256(content).hexdigest()
                retained.append((
                    relative, descriptor, _stable_fields(after), digest,
                    stat.S_IMODE(after.st_mode)))
                records.append({
                    "path": relative,
                    "sha256": digest,
                    "size": after.st_size,
                    "mode": stat.S_IMODE(after.st_mode),
                })
                if output_tree is not None:
                    output_tree.write_exclusive(
                        relative, content, stat.S_IMODE(after.st_mode))

            paths_after, git_after = _git_source_state(
                tree.root_descriptor, inherited_descriptors,
                git_snapshot, git_metadata)
            require(paths_after == paths_before and git_after == git_before,
                    "tracked source Git state changed during capture")
            for relative, descriptor, fields, digest, _mode in retained:
                current = os.fstat(descriptor)
                require(_stable_fields(current) == fields,
                        f"tracked source file {relative!r} changed while "
                        "retained")
                content = _read_bounded_descriptor(
                    descriptor, current.st_size,
                    f"retained tracked source file {relative!r}",
                    MAX_TRACKED_SOURCE_FILE_BYTES)
                require(hashlib.sha256(content).hexdigest() == digest,
                        f"tracked source file {relative!r} changed and was "
                        "restored while retained")
                reopened, reopened_status = tree.open_regular(relative)
                del reopened
                require(_stable_fields(reopened_status) == fields,
                        f"tracked source file {relative!r} pathname changed "
                        "while retained")
            git_metadata.verify()
            git_snapshot.verify()
            path_guard.verify()
    finally:
        path_guard.close()
    return {
        "schema": "leopard2-tracked-source-tree/v1",
        "total_bytes": total,
        "files": records,
        "git": git_before,
        "git_tool": dict(git_snapshot.identity),
    }


class _RetainedPrivateSourceTree:
    """Guard a runner-owned tracked source snapshot throughout clean replay."""

    def __init__(
        self, source: Path, destination: Path,
        *, inherited_descriptors: Sequence[int] = (),
    ) -> None:
        self.destination = destination
        self.manifest = _capture_tracked_source_tree(
            source, destination,
            inherited_descriptors=inherited_descriptors)
        self.guard = _InotifyMutationGuard(
            "runner-owned tracked source snapshot")
        try:
            self.guard.add_tree(destination)
            self._verify_files()
            self.guard.verify()
        except BaseException as primary:
            try:
                self.guard._close_without_verification()
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    "runner-owned source constructor",
                    "runner-owned source constructor cleanup failed: "
                    f"{cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            raise

    def _verify_files(self) -> None:
        expected = self.manifest.get("files")
        require(isinstance(expected, list),
                "runner-owned source manifest is malformed")
        with ExitStack() as stack:
            tree = _OpenDirectoryTree(
                self.destination, stack, "runner-owned tracked source")
            observed = []
            total = 0
            for record in expected:
                require(isinstance(record, dict) and
                        isinstance(record.get("path"), str),
                        "runner-owned source manifest record is malformed")
                relative = record["path"]
                descriptor, status = tree.open_regular(relative)
                content = _read_bounded_descriptor(
                    descriptor, status.st_size,
                    f"runner-owned source file {relative!r}",
                    MAX_TRACKED_SOURCE_FILE_BYTES)
                total += len(content)
                observed.append({
                    "path": relative,
                    "sha256": hashlib.sha256(content).hexdigest(),
                    "size": len(content),
                    "mode": stat.S_IMODE(status.st_mode),
                })
            require(total == self.manifest.get("total_bytes") and
                    observed == expected,
                    "runner-owned tracked source snapshot changed")

    def verify(self) -> None:
        self.guard.verify()
        self._verify_files()
        self.guard.verify()

    def close(self) -> None:
        failure: BaseException | None = None
        try:
            self.verify()
        except BaseException as error:
            failure = error
        try:
            self.guard.close()
        except BaseException as error:
            failure = _owner_exception_precedence(
                failure, error, "runner-owned source guard close")
        if failure is not None:
            raise failure

    def __enter__(self) -> "_RetainedPrivateSourceTree":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, tb
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc is None:
                raise
            _raise_owner_exit_failure(
                exc, cleanup_error, "runner-owned source owner",
                "runner-owned source cleanup failed: "
                f"{cleanup_error}; primary failure: "
                f"{type(exc).__name__}: {exc}")


def _verify_tracked_source_manifest(
    source: Path, manifest: Mapping[str, Any],
) -> dict[str, Any]:
    require(manifest.get("schema") == "leopard2-tracked-source-tree/v1" and
            isinstance(manifest.get("files"), list) and
            isinstance(manifest.get("git"), dict) and
            isinstance(manifest.get("git_tool"), dict),
            "tracked source manifest is malformed")
    git = manifest["git"]
    git_tool = manifest["git_tool"]
    require(isinstance(git.get("commit"), str) and
            re.fullmatch(r"[0-9a-f]{40,64}", git["commit"]) is not None and
            isinstance(git.get("tree"), str) and
            re.fullmatch(r"[0-9a-f]{40,64}", git["tree"]) is not None and
            type(git.get("dirty")) is bool and
            isinstance(git.get("status_sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", git["status_sha256"]) is not None,
            "tracked source manifest Git identity is malformed")
    require(isinstance(git_tool.get("path"), str) and
            isinstance(git_tool.get("sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", git_tool["sha256"]) is not None and
            type(git_tool.get("size")) is int and git_tool["size"] >= 0 and
            type(git_tool.get("device")) is int and
            type(git_tool.get("inode")) is int,
            "tracked source manifest Git tool identity is malformed")
    expected = manifest["files"]
    require(expected and len(expected) <= MAX_TRACKED_SOURCE_FILES,
            "tracked source manifest file count is invalid")
    expected_paths = [
        record.get("path") if isinstance(record, dict) else None
        for record in expected
    ]
    require(all(isinstance(path, str) for path in expected_paths) and
            expected_paths == sorted(expected_paths) and
            len(expected_paths) == len(set(expected_paths)) and
            type(manifest.get("total_bytes")) is int and
            0 <= manifest["total_bytes"] <= MAX_TRACKED_SOURCE_TOTAL_BYTES,
            "tracked source manifest path order or total is invalid")
    with ExitStack() as stack:
        tree = _OpenDirectoryTree(source, stack, "tracked source replay")
        observed = []
        total = 0
        for record in expected:
            require(isinstance(record, dict) and
                    isinstance(record.get("path"), str) and
                    type(record.get("size")) is int and
                    type(record.get("mode")) is int and
                    isinstance(record.get("sha256"), str) and
                    re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is not
                    None and
                    0 <= record["size"] <= MAX_TRACKED_SOURCE_FILE_BYTES and
                    0 <= record["mode"] <= 0o7777,
                    "tracked source manifest file record is malformed")
            relative = record["path"]
            descriptor, status = tree.open_regular(relative)
            require(status.st_size <= MAX_TRACKED_SOURCE_FILE_BYTES,
                    f"tracked source replay file {relative!r} exceeds its "
                    "individual byte bound")
            content = _read_bounded_descriptor(
                descriptor, status.st_size,
                f"tracked source replay file {relative!r}",
                MAX_TRACKED_SOURCE_FILE_BYTES)
            total += len(content)
            require(total <= MAX_TRACKED_SOURCE_TOTAL_BYTES,
                    "tracked source replay exceeds its total byte bound")
            observed.append({
                "path": relative,
                "sha256": hashlib.sha256(content).hexdigest(),
                "size": len(content),
                "mode": stat.S_IMODE(status.st_mode),
            })
        require(observed == expected and
                total == manifest.get("total_bytes"),
                "tracked source replay differs from its retained manifest")
    return {
        "schema": manifest["schema"],
        "total_bytes": total,
        "files": [dict(record) for record in expected],
        "git": dict(git),
        "git_tool": dict(git_tool),
    }


def _linux_prctl(option: int, argument: object) -> None:
    """Invoke the Linux child-subreaper controls used by command containment."""
    require(sys.platform.startswith("linux"),
            "bounded command containment requires Linux")
    try:
        function = ctypes.CDLL(None, use_errno=True).prctl
    except (AttributeError, OSError) as error:
        raise BuildProvenanceError(
            f"Linux child-subreaper prctl is unavailable: {error}") from error
    ctypes.set_errno(0)
    result = function(
        ctypes.c_int(option), argument, ctypes.c_ulong(0),
        ctypes.c_ulong(0), ctypes.c_ulong(0))
    if result != 0:
        number = ctypes.get_errno()
        raise BuildProvenanceError(
            "Linux child-subreaper prctl failed: " +
            os.strerror(number or errno.EPERM))


def _get_child_subreaper() -> int:
    value = ctypes.c_int(-1)
    _linux_prctl(PR_GET_CHILD_SUBREAPER, ctypes.byref(value))
    require(value.value in (0, 1),
            "Linux returned an invalid child-subreaper state")
    return value.value


def _set_child_subreaper(value: int) -> None:
    require(value in (0, 1), "invalid child-subreaper state")
    _linux_prctl(PR_SET_CHILD_SUBREAPER, ctypes.c_ulong(value))
    require(_get_child_subreaper() == value,
            "Linux did not retain the requested child-subreaper state")


def _read_proc_process_record_descriptor(
    pid: int, directory_descriptor: int,
) -> tuple[int, int, int, int, str, int, int] | None:
    """Read one process through a caller-retained no-follow proc directory."""
    stat_descriptor = -1
    try:
        directory_before = os.fstat(directory_descriptor)
        stat_descriptor = os.open(
            "stat",
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_descriptor)
        data = os.read(stat_descriptor, 4096)
        require(data and len(data) < 4096 and not os.read(stat_descriptor, 1),
                f"Linux process {pid} has oversized procfs stat data")
        directory_after = os.fstat(directory_descriptor)
        require(_stable_fields(directory_before) ==
                _stable_fields(directory_after),
                f"Linux process {pid} procfs directory changed while read")
    except (FileNotFoundError, ProcessLookupError):
        return None
    except OSError as error:
        if error.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise BuildProvenanceError(
            f"cannot inspect Linux process {pid}: {error}") from error
    finally:
        primary = sys.exception()
        if stat_descriptor >= 0:
            _close_raw_descriptor_with_precedence(
                stat_descriptor, primary,
                f"Linux process {pid} procfs stat descriptor",
                f"cannot close Linux process {pid} procfs stat descriptor")
    closing = data.rfind(b")")
    require(closing > 0 and closing + 2 < len(data),
            f"Linux process {pid} has malformed procfs stat data")
    fields = data[closing + 2:].split()
    require(len(fields) >= 20,
            f"Linux process {pid} has truncated procfs stat data")
    try:
        state = fields[0].decode("ascii")
        ppid, pgrp, session = int(fields[1]), int(fields[2]), int(fields[3])
        starttime = int(fields[19])
    except (UnicodeDecodeError, ValueError) as error:
        raise BuildProvenanceError(
            f"Linux process {pid} has invalid procfs stat fields") from error
    require(len(state) == 1 and min(ppid, pgrp, session, starttime) >= 0,
            f"Linux process {pid} has invalid procfs identity")
    return (
        ppid, pgrp, session, starttime, state,
        directory_after.st_dev, directory_after.st_ino,
    )


def _open_proc_process_record(
    pid: int,
) -> tuple[tuple[int, int, int, int, str, int, int], int] | None:
    """Open and return an exact process record with its directory FD held."""
    directory_descriptor = -1
    try:
        directory_descriptor = os.open(
            f"/proc/{pid}",
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0))
        record = _read_proc_process_record_descriptor(
            pid, directory_descriptor)
        if record is None:
            _close_raw_descriptor_with_precedence(
                directory_descriptor, None,
                f"vanished Linux process {pid} procfs directory",
                f"cannot close vanished Linux process {pid} procfs "
                "directory")
            return None
        path_status = os.stat(
            f"/proc/{pid}", follow_symlinks=False)
        descriptor_status = os.fstat(directory_descriptor)
        require(_stable_fields(path_status) ==
                _stable_fields(descriptor_status) and
                record[5:] ==
                (descriptor_status.st_dev, descriptor_status.st_ino),
                f"Linux process {pid} procfs pathname changed while opened")
        return record, directory_descriptor
    except (FileNotFoundError, ProcessLookupError) as primary:
        if directory_descriptor >= 0:
            _close_raw_descriptor_with_precedence(
                directory_descriptor, primary,
                f"Linux process {pid} procfs directory",
                f"cannot close Linux process {pid} procfs directory")
        return None
    except OSError as error:
        if directory_descriptor >= 0:
            _close_raw_descriptor_with_precedence(
                directory_descriptor, error,
                f"Linux process {pid} procfs directory",
                f"cannot close Linux process {pid} procfs directory")
        if error.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise BuildProvenanceError(
            f"cannot inspect Linux process {pid}: {error}") from error
    except BaseException as primary:
        if directory_descriptor >= 0:
            _close_raw_descriptor_with_precedence(
                directory_descriptor, primary,
                f"Linux process {pid} procfs directory",
                f"cannot close Linux process {pid} procfs directory")
        raise


def _proc_process_record(
    pid: int,
) -> tuple[int, int, int, int, str, int, int] | None:
    """Return a process record while closing its temporary proc descriptor."""
    opened = _open_proc_process_record(pid)
    if opened is None:
        return None
    record, descriptor = opened
    try:
        return record
    finally:
        _close_raw_descriptor_with_precedence(
            descriptor, sys.exception(),
            f"Linux process {pid} temporary procfs directory",
            f"cannot close Linux process {pid} temporary procfs directory")


def _close_retained_mapping_descriptor(
    mapping: dict[Any, int], key: Any,
) -> None:
    """Close one mapped FD without same-inode numeric-descriptor ABA."""
    descriptor = mapping[key]
    consumed = [False]
    try:
        _close_descriptor_with_ofd_guard(
            descriptor, "retained mapping descriptor", consumed)
    finally:
        if consumed[0] and mapping.get(key) == descriptor:
            del mapping[key]


class _ProcProcessSnapshot:
    """One procfs snapshot whose task-directory FDs stay open while used."""

    def __init__(self) -> None:
        require(Path("/proc/self/stat").is_file(),
                "bounded command containment requires mounted Linux procfs")
        self.records: dict[
            int, tuple[int, int, int, int, str, int, int]] = {}
        self.descriptors: dict[int, int] = {}
        try:
            names = os.listdir("/proc")
        except OSError as error:
            raise BuildProvenanceError(
                f"cannot enumerate Linux procfs: {error}") from error
        try:
            for name in names:
                if not name.isascii() or not name.isdigit():
                    continue
                pid = int(name)
                try:
                    opened = _open_proc_process_record(pid)
                except BuildProvenanceError:
                    try:
                        owner = os.stat(
                            f"/proc/{name}", follow_symlinks=False).st_uid
                    except OSError:
                        continue
                    if owner == os.getuid():
                        raise
                    continue
                if opened is not None:
                    record, descriptor = opened
                    self.records[pid] = record
                    self.descriptors[pid] = descriptor
            require(os.getpid() in self.records,
                    "Linux procfs does not expose the provenance runner")
        except BaseException as primary:
            try:
                self.close()
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    "Linux procfs snapshot constructor",
                    "Linux procfs snapshot constructor cleanup failed: "
                    f"{cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            raise

    def close(self) -> None:
        errors = []
        failure: BaseException | None = None
        for pid in list(self.descriptors):
            try:
                _close_retained_mapping_descriptor(
                    self.descriptors, pid)
            except BaseException as error:
                errors.append(str(error))
                observed = _ordinary_owner_cleanup_error(
                    error,
                    "cannot close Linux procfs snapshot descriptor "
                    f"{pid}: {error}")
                failure = _owner_exception_precedence(
                    failure, observed,
                    "Linux procfs snapshot descriptor cleanup")
            if pid not in self.descriptors:
                self.records.pop(pid, None)
        if not self.descriptors:
            self.records.clear()
        if isinstance(failure, _OWNER_TERMINAL_EXCEPTIONS):
            raise failure.with_traceback(failure.__traceback__)
        require(not errors,
                "cannot close Linux procfs snapshot descriptors: " +
                "; ".join(errors))

    def __enter__(self) -> "_ProcProcessSnapshot":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, tb
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc is None:
                raise
            _raise_owner_exit_failure(
                exc, cleanup_error, "procfs snapshot owner",
                "procfs snapshot cleanup failed: "
                f"{cleanup_error}; primary failure: "
                f"{type(exc).__name__}: {exc}")


def _proc_process_snapshot() -> _ProcProcessSnapshot:
    """Snapshot procfs and retain every observed task directory."""
    return _ProcProcessSnapshot()


def _linux_pidfd_open(pid: int) -> int | None:
    """Open a race-free process handle with a libc fallback."""
    python_wrapper = getattr(os, "pidfd_open", None)
    if callable(python_wrapper):
        try:
            return python_wrapper(pid, 0)
        except OSError as error:
            if error.errno == errno.ESRCH:
                return None
            raise BuildProvenanceError(
                f"cannot open pidfd for process {pid}: {error}") from error
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_open
    except (AttributeError, OSError) as error:
        raise BuildProvenanceError(
            f"Linux pidfd_open is unavailable: {error}") from error
    ctypes.set_errno(0)
    descriptor = function(ctypes.c_int(pid), ctypes.c_uint(0))
    if descriptor >= 0:
        return descriptor
    number = ctypes.get_errno()
    if number == errno.ESRCH:
        return None
    raise BuildProvenanceError(
        f"cannot open pidfd for process {pid}: " +
        os.strerror(number or errno.EPERM))


def _linux_pidfd_signal(descriptor: int, signal_number: int) -> None:
    python_wrapper = getattr(signal, "pidfd_send_signal", None)
    if callable(python_wrapper):
        try:
            python_wrapper(descriptor, signal_number, None, 0)
            return
        except OSError as error:
            if error.errno == errno.ESRCH:
                return
            raise BuildProvenanceError(
                "cannot signal a contained process through pidfd: " +
                str(error)) from error
    try:
        function = ctypes.CDLL(None, use_errno=True).pidfd_send_signal
    except (AttributeError, OSError) as error:
        raise BuildProvenanceError(
            f"Linux pidfd_send_signal is unavailable: {error}") from error
    ctypes.set_errno(0)
    result = function(
        ctypes.c_int(descriptor), ctypes.c_int(signal_number), None,
        ctypes.c_uint(0))
    if result != 0:
        number = ctypes.get_errno()
        if number == errno.ESRCH:
            return
        raise BuildProvenanceError(
            "cannot signal a contained process through pidfd: " +
            os.strerror(number or errno.EPERM))


def _validate_pidfd_support() -> None:
    descriptor = _linux_pidfd_open(os.getpid())
    require(descriptor is not None,
            "Linux pidfd support cannot identify the provenance runner")
    try:
        _linux_pidfd_signal(descriptor, 0)
    finally:
        _close_raw_descriptor_with_precedence(
            descriptor, sys.exception(),
            "Linux pidfd support probe",
            "cannot close Linux pidfd support-probe descriptor")


class _LinuxDescendantContainment:
    """Own and reap the complete tree of one direct command.

    Process groups are insufficient because a command can call setsid() and
    double-fork while retaining benchmark locks or output descriptors.  A
    temporary subreaper makes every such orphan our direct child.  Procfs
    ancestry plus pidfds identifies and kills the exact tree after successful
    commands as well as after errors and timeouts.
    """

    def __init__(self) -> None:
        self.runner_pid = os.getpid()
        self.previous_subreaper: int | None = None
        self.leader: tuple[int, int, int, int] | None = None
        self.known: set[tuple[int, int, int, int]] = set()
        self.pidfds: dict[tuple[int, int, int, int], int] = {}
        self.procfds: dict[tuple[int, int, int, int], int] = {}
        self.process: subprocess.Popen[bytes] | None = None
        self.active = False
        self.proven_empty = False

    @staticmethod
    def _direct_children(
        snapshot: Mapping[
            int, tuple[int, int, int, int, str, int, int]],
        parent: int,
    ) -> set[tuple[int, int, int, int]]:
        return {
            (pid, record[3], record[5], record[6])
            for pid, record in snapshot.items()
                if record[0] == parent}

    def __enter__(self) -> "_LinuxDescendantContainment":
        require(sys.platform.startswith("linux"),
                "bounded command containment requires Linux")
        try:
            task_count = sum(
                1 for name in os.listdir("/proc/self/task")
                if name.isascii() and name.isdigit())
        except OSError as error:
            raise BuildProvenanceError(
                f"cannot inspect runner threads: {error}") from error
        require(task_count == 1,
                "bounded command containment requires a single-threaded runner")
        _validate_pidfd_support()
        self.previous_subreaper = _get_child_subreaper()
        try:
            _set_child_subreaper(1)
            with _proc_process_snapshot() as snapshot:
                require(not self._direct_children(
                            snapshot.records, self.runner_pid),
                        "bounded command containment found pre-existing "
                        "children")
            self.active = True
            return self
        except BaseException as primary:
            try:
                _set_child_subreaper(self.previous_subreaper)
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    "Linux descendant containment constructor",
                    "Linux descendant containment constructor cleanup "
                    f"failed: {cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            self.previous_subreaper = None
            raise

    def attach(self, process: subprocess.Popen[bytes]) -> None:
        require(self.active and self.process is process and
                self.leader is None and process.pid > 0,
                "invalid bounded command attachment")
        opened = _open_proc_process_record(process.pid)
        require(opened is not None and opened[0][0] == self.runner_pid,
                "spawned command is not an owned direct child")
        record, observed_descriptor = opened
        self.leader = (
            process.pid, record[3], record[5], record[6])
        try:
            require(self._retain_pidfd(
                        self.leader, record,
                        observed_descriptor=observed_descriptor),
                    "spawned command escaped before its pidfd was retained")
        finally:
            _close_raw_descriptor_with_precedence(
                observed_descriptor, sys.exception(),
                "Linux descendant attachment procfs descriptor",
                "cannot close Linux descendant attachment procfs "
                "descriptor")
        self.known.add(self.leader)

    @staticmethod
    def _same_process_record(
        first: tuple[int, int, int, int, str, int, int],
        second: tuple[int, int, int, int, str, int, int],
    ) -> bool:
        # State is allowed to change while pidfd_open runs.  Parent, process
        # group, session, start time, and the held/no-follow proc-directory
        # inode together bind the observed process.
        return first[:4] == second[:4] and first[5:] == second[5:]

    def _retain_pidfd(
        self,
        identity: tuple[int, int, int, int],
        expected: tuple[int, int, int, int, str, int, int],
        *, observed_descriptor: int | None = None,
    ) -> bool:
        """Retain a race-free handle to one exact observed descendant."""
        if identity in self.pidfds:
            return True
        pid, starttime, proc_device, proc_inode = identity
        require(expected[3] == starttime and
                expected[5:] == (proc_device, proc_inode),
                "process snapshot identity does not match its start time")
        owned_observation = -1
        retained_proc = -1
        descriptor = -1
        retained_successfully = False
        try:
            if observed_descriptor is None:
                opened = _open_proc_process_record(pid)
                if opened is None:
                    return False
                before, owned_observation = opened
                observed_descriptor = owned_observation
            else:
                before = _read_proc_process_record_descriptor(
                    pid, observed_descriptor)
            if (before is None or
                    not self._same_process_record(before, expected)):
                return False
            retained_proc = os.dup(observed_descriptor)
            # The exact proc inode is now pinned before pidfd_open.  Re-reading
            # through this same descriptor after pidfd_open detects death/PID
            # reuse even if starttime and a recycled numeric inode collide.
            opened_pidfd = _linux_pidfd_open(pid)
            if opened_pidfd is None:
                return False
            descriptor = opened_pidfd
            retained_after = _read_proc_process_record_descriptor(
                pid, retained_proc)
            opened_after = _open_proc_process_record(pid)
            if opened_after is None:
                return False
            after, after_descriptor = opened_after
            try:
                if (retained_after is None or
                        not self._same_process_record(
                            retained_after, expected) or
                        not self._same_process_record(after, expected)):
                    return False
            finally:
                _close_raw_descriptor_with_precedence(
                    after_descriptor, sys.exception(),
                    "Linux descendant post-pidfd procfs descriptor",
                    "cannot close Linux descendant post-pidfd procfs "
                    "descriptor")
            self.pidfds[identity] = descriptor
            self.procfds[identity] = retained_proc
            retained_successfully = True
        finally:
            primary = sys.exception()
            cleanup_failure: BaseException | None = None
            for cleanup_descriptor in (
                    owned_observation,
                    descriptor if not retained_successfully else -1,
                    retained_proc if not retained_successfully else -1):
                if cleanup_descriptor < 0:
                    continue
                try:
                    _close_raw_descriptor_with_precedence(
                        cleanup_descriptor, primary,
                        "Linux descendant retained raw descriptor",
                        "cannot close Linux descendant retained raw "
                        "descriptor")
                except BaseException as error:
                    cleanup_failure = _owner_exception_precedence(
                        cleanup_failure, error,
                        "Linux descendant retained-descriptor cleanup")
            if cleanup_failure is not None:
                if primary is not None:
                    _raise_owner_exit_failure(
                        primary, cleanup_failure,
                        "Linux descendant retained-descriptor owner",
                        "Linux descendant retained-descriptor cleanup "
                        f"failed: {cleanup_failure}; primary failure: "
                        f"{type(primary).__name__}: {primary}")
                raise cleanup_failure.with_traceback(
                    cleanup_failure.__traceback__)
        return retained_successfully

    def _close_pidfds(self) -> None:
        errors = []
        failure: BaseException | None = None
        for mapping in (self.pidfds, self.procfds):
            for identity in list(mapping):
                descriptor = mapping[identity]
                try:
                    _close_retained_mapping_descriptor(
                        mapping, identity)
                except BaseException as error:
                    errors.append(f"fd {descriptor}: {error}")
                    observed = _ordinary_owner_cleanup_error(
                        error,
                        f"cannot close retained process descriptor "
                        f"{descriptor}: {error}")
                    failure = _owner_exception_precedence(
                        failure, observed,
                        "Linux descendant retained process cleanup")
        if isinstance(failure, _OWNER_TERMINAL_EXCEPTIONS):
            raise failure.with_traceback(failure.__traceback__)
        require(not errors, "cannot close retained process descriptors: " +
                "; ".join(errors))

    def _signal_retained(
        self, targets: set[tuple[int, int, int, int]], signal_number: int,
    ) -> None:
        for identity in sorted(targets, reverse=True):
            descriptor = self.pidfds.get(identity)
            if descriptor is not None:
                _linux_pidfd_signal(descriptor, signal_number)

    def leader_exited(self) -> bool:
        require(self.leader is not None and self.leader in self.pidfds,
                "bounded command leader has no retained pidfd")
        descriptor = self.pidfds[self.leader]
        try:
            poller = select.poll()
            poller.register(descriptor, select.POLLIN)
            events = poller.poll(0)
        except (OSError, ValueError) as error:
            raise BuildProvenanceError(
                f"cannot poll bounded command leader pidfd: {error}") \
                from error
        return bool(events)

    def _discover(
        self,
        snapshot: _ProcProcessSnapshot | Mapping[
            int, tuple[int, int, int, int, str, int, int]],
    ) -> set[tuple[int, int, int, int]]:
        if isinstance(snapshot, _ProcProcessSnapshot):
            records = snapshot.records
            observed_descriptors = snapshot.descriptors
        else:
            records = snapshot
            observed_descriptors = {}
        targets = {
            identity for identity in self.known
            if (identity[0] in records and
                records[identity[0]][3] == identity[1] and
                records[identity[0]][5:] == identity[2:])
        }
        targets.update(self._direct_children(records, self.runner_pid))
        changed = True
        while changed:
            changed = False
            parents = {identity[0] for identity in targets}
            for pid, record in records.items():
                identity = (pid, record[3], record[5], record[6])
                if record[0] in parents and identity not in targets:
                    targets.add(identity)
                    changed = True
        for identity in sorted(targets):
            record = records.get(identity[0])
            if (record is not None and record[3] == identity[1] and
                    record[5:] == identity[2:]):
                if self._retain_pidfd(
                        identity, record,
                        observed_descriptor=observed_descriptors.get(
                            identity[0])):
                    # Only race-free retained identities survive beyond this
                    # ancestry-qualified snapshot.  A vanished observation may
                    # not later rebind an unrelated same-tick PID by number.
                    self.known.add(identity)
        return targets

    def _terminate_and_reap_owned_tree(
        self, process: subprocess.Popen[bytes] | None,
    ) -> None:
        deadline = time.monotonic() + CHILD_REAP_TIMEOUT_SECONDS
        empty_scans = 0
        while True:
            with _proc_process_snapshot() as snapshot:
                targets = self._discover(snapshot)
                self._signal_retained(targets, signal.SIGKILL)
                records = dict(snapshot.records)
            for identity in sorted(self.known):
                if process is not None and identity == self.leader:
                    continue
                record = records.get(identity[0])
                if (record is None or record[3] != identity[1] or
                        record[5:] != identity[2:]):
                    continue
                try:
                    os.waitpid(identity[0], os.WNOHANG)
                except (ChildProcessError, ProcessLookupError):
                    pass
                except OSError as error:
                    raise BuildProvenanceError(
                        f"cannot reap contained process {identity[0]}: "
                        f"{error}") from error

            with _proc_process_snapshot() as live_snapshot:
                live = self._discover(live_snapshot)
            other_live = (
                live if process is None else
                {identity for identity in live
                 if identity != self.leader})
            leader_done = process is None or self.leader_exited()
            if leader_done and not other_live:
                empty_scans += 1
                if empty_scans >= 2:
                    if process is not None:
                        try:
                            process.wait(timeout=0.1)
                        except subprocess.TimeoutExpired as error:
                            raise BuildProvenanceError(
                                "bounded command pidfd became readable "
                                "before the leader was waitable") from error
                    self.proven_empty = True
                    return
            else:
                empty_scans = 0
            remaining = deadline - time.monotonic()
            require(remaining > 0,
                    "bounded command descendants survived SIGKILL")
            time.sleep(min(0.01, remaining))

    def terminate_and_reap(self) -> None:
        process = self.process
        require(self.active and process is not None and
                self.leader is not None and self.leader[0] == process.pid,
                "bounded command containment is not attached")
        self._terminate_and_reap_owned_tree(process)

    def terminate_unattached_and_reap(self) -> None:
        require(self.active and self.process is None and self.leader is None,
                "unattached command containment has invalid state")
        self._terminate_and_reap_owned_tree(None)

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del tb
        if not self.active:
            return
        cleanup_error: BaseException | None = None
        restore_error: BaseException | None = None
        close_error: BaseException | None = None
        previous = self.previous_subreaper
        try:
            if self.process is not None and not self.proven_empty:
                try:
                    # This path deliberately runs even when the leader exited
                    # zero, proving that no detached FD holder survived.
                    if self.leader is None:
                        opened = _open_proc_process_record(self.process.pid)
                        if opened is not None:
                            record, observed_descriptor = opened
                            identity = (
                                self.process.pid, record[3],
                                record[5], record[6])
                            try:
                                if self._retain_pidfd(
                                        identity, record,
                                        observed_descriptor=
                                        observed_descriptor):
                                    self.leader = identity
                                    self.known.add(identity)
                            finally:
                                _close_raw_descriptor_with_precedence(
                                    observed_descriptor, sys.exception(),
                                    "Linux descendant exit procfs "
                                    "descriptor",
                                    "cannot close Linux descendant exit "
                                    "procfs descriptor")
                    if self.leader is None:
                        # The direct child may already have been reaped by an
                        # installed SIGCHLD handler while detached descendants
                        # remain adopted by this subreaper.  Fall back to the
                        # complete unattached-child sweep instead of leaking
                        # those descendants on the attachment error path.
                        unattached_process = self.process
                        self.process = None
                        self.terminate_unattached_and_reap()
                        unattached_process.poll()
                    else:
                        self.terminate_and_reap()
                except BaseException as error:
                    cleanup_error = error
            elif self.process is None:
                try:
                    # Popen can fail after creating a child, and an asynchronous
                    # exception can land before the returned object is attached.
                    # Prove that no such direct child or detached descendant
                    # survives before restoring the prior subreaper state.
                    self.terminate_unattached_and_reap()
                except BaseException as error:
                    cleanup_error = error
        finally:
            self.active = False
            self.previous_subreaper = None
            if previous is None:
                restore_error = BuildProvenanceError(
                    "previous child-subreaper state was lost")
            else:
                try:
                    _set_child_subreaper(previous)
                except BaseException as error:
                    restore_error = error
            try:
                self._close_pidfds()
            except BaseException as error:
                close_error = error
        if (cleanup_error is not None or restore_error is not None or
                close_error is not None):
            authoritative = (
                exc if isinstance(exc, BaseException) else None)
            for phase, error in (
                    ("descendant cleanup", cleanup_error),
                    ("subreaper restore", restore_error),
                    ("pidfd close", close_error)):
                if error is not None:
                    authoritative = _owner_exception_precedence(
                        authoritative, error,
                        f"Linux descendant containment {phase}")
            details = []
            if cleanup_error is not None:
                details.append(
                    "descendant cleanup failed: " +
                    f"{type(cleanup_error).__name__}: {cleanup_error}")
            if restore_error is not None:
                details.append(
                    "subreaper restore failed: " +
                    f"{type(restore_error).__name__}: {restore_error}")
            if close_error is not None:
                details.append(
                    "pidfd close failed: " +
                    f"{type(close_error).__name__}: {close_error}")
            if exc is not None:
                details.append(
                    f"primary failure: {type(exc).__name__}: {exc}")
            if isinstance(authoritative, _OWNER_TERMINAL_EXCEPTIONS):
                raise authoritative.with_traceback(
                    authoritative.__traceback__)
            raise BuildProvenanceError("; ".join(details)) from (
                cleanup_error or restore_error or close_error)


def _run(
    command: Sequence[str], label: str, *, maximum_bytes: int = 4 << 20,
    timeout: float = 120, inherited_descriptors: Sequence[int] = (),
    executable_descriptor: int | None = None,
    environment_overrides: Mapping[str, str] | None = None,
    write_sandbox_root: Path | str | None = None,
    write_sandbox_descriptors: Sequence[int] = (),
    private_tmpfs_root: Path | str | None = None,
    private_tmpfs_directories: Sequence[str] = (),
    private_tmpfs_descriptor: int | None = None,
) -> bytes:
    # This argv is passed directly to execve through Popen(shell=False), so
    # CR/LF and empty data arguments are not shell syntax.  Only argv[0] must
    # be nonempty, and NUL is the sole forbidden character in every element.
    require(isinstance(command, Sequence) and
            not isinstance(command, (str, bytes)) and command and
            isinstance(command[0], str) and command[0] and
            all(isinstance(item, str) and "\0" not in item
                for item in command),
            f"{label} command argv is invalid")
    require(type(maximum_bytes) is int and maximum_bytes >= 0,
            f"{label} output byte bound is invalid")
    require(isinstance(timeout, (int, float)) and
            not isinstance(timeout, bool) and
            math.isfinite(float(timeout)) and timeout > 0,
            f"{label} timeout is invalid")
    sandbox_owner: _OwnedDescriptor | None = None
    sandbox_root: Path | None = None
    tmpfs_owner: _OwnedDescriptor | None = None
    tmpfs_root: Path | None = None
    tmpfs_descriptor = -1
    outside_uid = os.getuid()
    outside_gid = os.getgid()
    descriptor_set = set(inherited_descriptors)
    require(all(type(descriptor) is int and descriptor >= 0
                for descriptor in write_sandbox_descriptors),
            f"{label} additional write descriptor set is invalid")
    descriptor_set.update(write_sandbox_descriptors)
    if write_sandbox_root is not None:
        sandbox_root = Path(os.path.abspath(os.fspath(write_sandbox_root)))
        require(sys.platform.startswith("linux") and
                sandbox_root.is_absolute(),
                f"{label} write sandbox root is invalid")
        sandbox_owner = _OwnedDescriptor()
        try:
            sandbox_descriptor = sandbox_owner.open(
                sandbox_root,
                getattr(os, "O_PATH", os.O_RDONLY) |
                getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_NOFOLLOW", 0))
            lexical_status = os.stat(
                sandbox_root, follow_symlinks=False)
            retained_status = os.fstat(sandbox_descriptor)
            require(
                stat.S_ISDIR(lexical_status.st_mode) and
                (lexical_status.st_dev, lexical_status.st_ino) ==
                (retained_status.st_dev, retained_status.st_ino),
                f"{label} write sandbox root is unstable or unsafe")
        except BaseException as primary:
            _close_owned_descriptor_after_failure(
                sandbox_owner, primary,
                f"{label} write-sandbox descriptor")
            raise
        descriptor_set.add(sandbox_descriptor)
    if private_tmpfs_root is not None:
        require(sandbox_root is not None,
                f"{label} private tmpfs requires a write sandbox")
        tmpfs_root = Path(os.path.abspath(os.fspath(private_tmpfs_root)))
        require(tmpfs_root.is_absolute() and
                tmpfs_root.is_relative_to(sandbox_root) and
                tmpfs_root != sandbox_root,
                f"{label} private tmpfs root is invalid")
        if private_tmpfs_descriptor is None:
            tmpfs_owner = _OwnedDescriptor()
            try:
                tmpfs_descriptor = tmpfs_owner.open(
                    tmpfs_root,
                    getattr(os, "O_PATH", os.O_RDONLY) |
                    getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_DIRECTORY", 0) |
                    getattr(os, "O_NOFOLLOW", 0))
            except BaseException as primary:
                owners = [tmpfs_owner]
                if sandbox_owner is not None:
                    owners.append(sandbox_owner)
                _close_owned_descriptors_after_failure(
                    owners, primary,
                    f"{label} private-tmpfs descriptor")
                raise
        else:
            require(type(private_tmpfs_descriptor) is int and
                    private_tmpfs_descriptor >= 0,
                    f"{label} private tmpfs descriptor is invalid")
            tmpfs_descriptor = private_tmpfs_descriptor
            try:
                tmpfs_status = os.fstat(tmpfs_descriptor)
                tmpfs_lexical = os.stat(
                    tmpfs_root, follow_symlinks=False)
            except OSError as error:
                if sandbox_owner is not None:
                    _close_owned_descriptor_after_failure(
                        sandbox_owner, error,
                        f"{label} private-tmpfs validation")
                raise BuildProvenanceError(
                    f"{label} private tmpfs descriptor is unavailable: "
                    f"{error}") from error
            require(
                stat.S_ISDIR(tmpfs_status.st_mode) and
                (tmpfs_status.st_dev, tmpfs_status.st_ino) ==
                (tmpfs_lexical.st_dev, tmpfs_lexical.st_ino),
                f"{label} private tmpfs descriptor differs from its root")
        descriptor_set.add(tmpfs_descriptor)
    else:
        require(not private_tmpfs_directories and
                private_tmpfs_descriptor is None,
                f"{label} private tmpfs directories lack a root")
    if executable_descriptor is not None:
        require(type(executable_descriptor) is int and
                executable_descriptor >= 0,
                f"{label} executable descriptor is invalid")
        descriptor_set.add(executable_descriptor)
    pass_fds = tuple(sorted(descriptor_set))
    require(all(type(descriptor) is int and descriptor >= 0
                for descriptor in pass_fds),
            f"{label} inherited descriptor set is invalid")
    for descriptor in pass_fds:
        try:
            os.fstat(descriptor)
        except OSError as error:
            raise BuildProvenanceError(
                f"{label} inherited descriptor {descriptor} is invalid: "
                f"{error}") from error
    environment = dict(GIT_ENVIRONMENT)
    if environment_overrides is not None:
        require(all(
                    isinstance(name, str) and
                    re.fullmatch(r"[A-Z][A-Z0-9_]*", name) is not None and
                    isinstance(value, str)
                    for name, value in environment_overrides.items()),
                f"{label} environment overrides are invalid")
        for name, value in environment_overrides.items():
            _require_safe_unicode(value, f"{label} environment {name}")
            require("\0" not in value,
                    f"{label} environment {name} contains NUL")
            environment[name] = value

    process: subprocess.Popen[bytes] | None = None
    selector = selectors.DefaultSelector()
    stdout = bytearray()
    stderr = bytearray()
    failure: str | None = None
    returncode = -int(signal.SIGKILL)

    def child_restrictions() -> None:
        if tmpfs_root is not None:
            _mount_private_replay_tmpfs(
                tmpfs_root, tmpfs_descriptor,
                private_tmpfs_directories, outside_uid, outside_gid)
        if sandbox_root is not None and sandbox_owner is not None:
            write_descriptors = tuple(dict.fromkeys((
                *write_sandbox_descriptors,
                *((tmpfs_descriptor,) if tmpfs_root is not None else ()),
            )))
            _landlock_restrict_writes_to(
                sandbox_root, sandbox_owner.descriptor,
                write_descriptors)

    try:
        with _LinuxDescendantContainment() as containment:
            process = subprocess.Popen(
                list(command), stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env=environment, pass_fds=pass_fds,
                start_new_session=True,
                preexec_fn=(
                    child_restrictions
                    if sandbox_root is not None else None),
                executable=(
                    f"/proc/self/fd/{executable_descriptor}"
                    if executable_descriptor is not None else None))
            containment.process = process
            containment.attach(process)
            require(process.stdout is not None and process.stderr is not None,
                    f"{label} output pipes are unavailable")
            streams = {
                process.stdout.fileno(): (process.stdout, stdout,
                                          maximum_bytes),
                process.stderr.fileno(): (process.stderr, stderr,
                                          maximum_bytes),
            }
            for stream, _output, _limit in streams.values():
                os.set_blocking(stream.fileno(), False)
                selector.register(stream, selectors.EVENT_READ)

            def consume(
                events: Sequence[tuple[selectors.SelectorKey, int]],
            ) -> None:
                nonlocal failure
                for key, _mask in events:
                    descriptor = key.fileobj.fileno()
                    stream, output, limit = streams[descriptor]
                    try:
                        block = os.read(descriptor, 65536)
                    except BlockingIOError:
                        continue
                    if not block:
                        selector.unregister(stream)
                        continue
                    output.extend(block)
                    if len(output) > limit and failure is None:
                        failure = (
                            f"{label} output exceeds its retained byte bound")

            deadline = time.monotonic() + float(timeout)
            while True:
                if containment.leader_exited():
                    break
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    failure = (
                        f"{label} exceeded {float(timeout):.3f} seconds")
                    break
                events = selector.select(min(remaining, 0.1))
                consume(events)
                if failure is not None:
                    break

            # Required on every path, including rc=0, before any inherited lock
            # can leave this call's lifetime.  Do this as soon as the direct
            # leader exits: a detached descendant may intentionally retain the
            # stdout/stderr write ends, so waiting for pipe EOF first would
            # misclassify the successful direct command as a timeout.
            containment.terminate_and_reap()
            if isinstance(process.returncode, int):
                returncode = process.returncode

            # Every contained writer is now gone.  Drain bytes already buffered
            # in the nonblocking pipes and require prompt EOF; an open pipe with
            # no owned writer would mean a descriptor escaped containment.
            drain_deadline = (
                time.monotonic() + CHILD_REAP_TIMEOUT_SECONDS)
            while selector.get_map():
                remaining = drain_deadline - time.monotonic()
                require(remaining > 0,
                        f"{label} output pipes remained open after "
                        "descendant cleanup")
                events = selector.select(min(remaining, 0.05))
                if not events:
                    continue
                consume(events)
                if failure is not None and (
                        len(stdout) > maximum_bytes or
                        len(stderr) > maximum_bytes):
                    break
    finally:
        primary = sys.exception()
        cleanup_failure: BaseException | None = None
        cleanup_details: list[str] = []

        def close_resource(resource: Any, description: str) -> None:
            nonlocal cleanup_failure
            try:
                resource.close()
            except BaseException as error:
                cleanup_details.append(
                    f"{description}: {type(error).__name__}: {error}")
                observed = _ordinary_owner_cleanup_error(
                    error, f"cannot close {description}: {error}")
                cleanup_failure = _owner_exception_precedence(
                    cleanup_failure, observed,
                    f"{label} command resource cleanup")

        close_resource(selector, f"{label} output selector")
        if process is not None:
            if process.stdout is not None:
                close_resource(process.stdout, f"{label} stdout pipe")
            if process.stderr is not None:
                close_resource(process.stderr, f"{label} stderr pipe")
        if sandbox_owner is not None:
            close_resource(
                sandbox_owner, f"{label} write-sandbox descriptor")
        if tmpfs_owner is not None:
            close_resource(
                tmpfs_owner, f"{label} private-tmpfs descriptor")
        if cleanup_failure is not None:
            details = "; ".join(cleanup_details)
            if primary is not None:
                _raise_owner_exit_failure(
                    primary, cleanup_failure,
                    f"{label} command owner",
                    f"{label} command cleanup failed: {details}; "
                    f"primary failure: {type(primary).__name__}: {primary}")
            if isinstance(
                    cleanup_failure, _OWNER_TERMINAL_EXCEPTIONS):
                raise cleanup_failure.with_traceback(
                    cleanup_failure.__traceback__)
            raise BuildProvenanceError(
                f"{label} command cleanup failed: {details}") from (
                    cleanup_failure)

    if failure is not None:
        raise BuildProvenanceError(failure)
    require(returncode == 0,
            f"{label} failed with rc={returncode}: "
            f"{stderr.decode(errors='replace').strip()}")
    return bytes(stdout)


def _replay_invocation_record(
    command: Sequence[str], *,
    executable_label: str,
    executable: _RetainedFileSnapshot,
    inherited_descriptors: Sequence[int],
    environment_overrides: Mapping[str, str] | None,
    maximum_bytes: int,
    timeout: float,
    write_sandbox_root: Path | str,
    write_sandbox_descriptors: Sequence[int] = (),
    private_tmpfs_root: Path | str | None = None,
    private_tmpfs_directories: Sequence[str] = (),
    private_tmpfs_descriptor: int | None = None,
) -> dict[str, Any]:
    """Publish the exact argv, environment, and immutable exec binding."""
    require(isinstance(executable_label, str) and executable_label,
            "replay invocation executable label is invalid")
    descriptor = executable.executable_descriptor
    pass_fds = sorted({descriptor, *inherited_descriptors})
    environment = dict(GIT_ENVIRONMENT)
    if environment_overrides is not None:
        environment.update(environment_overrides)
    sandbox_root = Path(os.path.abspath(os.fspath(write_sandbox_root)))
    require(sandbox_root.is_absolute(),
            "replay invocation write sandbox root is invalid")
    require(all(type(descriptor) is int and descriptor >= 0
                for descriptor in write_sandbox_descriptors),
            "replay invocation write descriptor set is invalid")
    private_tmpfs = None
    if private_tmpfs_root is not None:
        tmpfs_root = Path(os.path.abspath(os.fspath(private_tmpfs_root)))
        require(
            tmpfs_root.is_absolute() and
            tmpfs_root.is_relative_to(sandbox_root) and
            type(private_tmpfs_descriptor) is int and
            private_tmpfs_descriptor >= 0 and
            isinstance(private_tmpfs_directories, Sequence) and
            not isinstance(private_tmpfs_directories, (str, bytes)),
            "replay invocation private tmpfs is invalid")
        private_tmpfs = {
            "root": str(tmpfs_root),
            "descriptor": private_tmpfs_descriptor,
            "bytes": PRIVATE_REPLAY_TMPFS_BYTES,
            "directories": sorted(set(private_tmpfs_directories)),
        }
    else:
        require(not private_tmpfs_directories and
                private_tmpfs_descriptor is None,
                "replay invocation tmpfs directories lack a root")
    return {
        "schema": REPLAY_INVOCATION_SCHEMA,
        "argv": list(command),
        "environment": {
            name: environment[name] for name in sorted(environment)
        },
        "pass_fds": pass_fds,
        "maximum_output_bytes": maximum_bytes,
        "timeout_seconds": timeout,
        "write_sandbox_root": str(sandbox_root),
        "write_sandbox_descriptors":
            sorted(set(write_sandbox_descriptors)),
        "private_tmpfs": private_tmpfs,
        "executable": {
            "label": executable_label,
            "logical_argv0": command[0],
            "sealed_descriptor": descriptor,
            "sealed_sha256": executable.executable_identity["sha256"],
            "sealed_size": executable.executable_identity["size"],
            "sealed_seals": executable.executable_identity["seals"],
        },
    }


def parse_cmake_cache(content: bytes) -> dict[str, str]:
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise BuildProvenanceError("CMakeCache.txt is not strict UTF-8") from error
    require("\0" not in text and "\r" not in text,
            "CMakeCache.txt contains a forbidden delimiter")
    result: dict[str, str] = {}
    # CMakeCache.txt is LF-framed.  splitlines() would incorrectly interpret
    # Unicode line separators as records.  They are unsafe evidence content
    # and are rejected below together with controls and bidi format marks.
    for line in text.split("\n"):
        _require_safe_unicode(line, "CMakeCache.txt")
        if not line or line.startswith(("#", "//")):
            continue
        require("=" in line,
                "CMakeCache.txt contains an unframed cache record")
        name_and_type, value = line.split("=", 1)
        require(name_and_type.count(":") == 1,
                "CMakeCache.txt contains a malformed typed key")
        name, entry_type = name_and_type.split(":", 1)
        require(bool(name) and
                re.fullmatch(r"[A-Za-z0-9_.+-]+", name) is not None and
                entry_type in CMAKE_CACHE_ENTRY_TYPES and
                name not in result,
                f"CMakeCache.txt contains a malformed or duplicate key "
                f"{name!r}")
        allowed_types = CMAKE_CACHE_REQUIRED_ENTRY_TYPES.get(name)
        require(allowed_types is None or entry_type in allowed_types,
                f"CMakeCache.txt key {name} has type {entry_type}, "
                f"expected one of {sorted(allowed_types or ())}")
        result[name] = value
    return result


def _cmake_true(value: str | None) -> bool:
    return (value or "").upper() in {"1", "ON", "TRUE", "YES", "Y"}


def _cmake_false(value: str | None) -> bool:
    return (value or "").upper() in {
        "0", "OFF", "FALSE", "NO", "N", "IGNORE", "NOTFOUND", "",
    } or (value or "").upper().endswith("-NOTFOUND")


def _resolve_build_operand(build: Path, directory: Path, value: str) -> Path:
    operand = Path(value)
    if not operand.is_absolute():
        operand = directory / operand
    resolved = operand.resolve(strict=False)
    require(resolved.is_relative_to(build),
            f"build operand escapes the candidate build: {value}")
    return resolved


def _lexical_build_operand(build: Path, directory: Path, value: str) -> Path:
    """Resolve a build operand lexically without consulting filesystem links."""
    operand = Path(value)
    lexical = Path(os.path.abspath(os.fspath(
        operand if operand.is_absolute() else directory / operand)))
    require(lexical.is_relative_to(build),
            f"build operand escapes the candidate build: {value}")
    return lexical


def _require_shell_literal_tokens(
    tokens: Sequence[str], label: str,
) -> None:
    """Reject shell syntax even when shlex leaves it fused into one token."""
    def response_indirection(token: str) -> bool:
        return token.startswith("@") or token.startswith("-Wl,@")

    require(tokens and all(
                token and not response_indirection(token)
                for token in tokens),
            f"{label} contains response-file indirection")
    require(tokens and all(
                token and not any(character in token
                                  for character in SHELL_RECIPE_META) and
                not any(character in token for character in "\0\r\n")
                for token in tokens),
            f"{label} contains shell control, substitution, or expansion "
            "syntax")


def _recipe_commands(content: bytes, label: str) -> list[list[str]]:
    """Parse a CMake Unix-Makefiles link recipe without shell delegation."""
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise BuildProvenanceError(f"{label} is not strict UTF-8") from error
    require("\0" not in text and "\r" not in text,
            f"{label} contains a forbidden delimiter")
    commands: list[list[str]] = []
    for line in text.split("\n"):
        if not line.strip():
            continue
        try:
            tokens = shlex.split(line, posix=True)
        except ValueError as error:
            raise BuildProvenanceError(
                f"cannot parse {label}: {error}") from error
        require(tokens and all(
                    token and not token.startswith(("@", "-Wl,@"))
                    for token in tokens),
                f"{label} contains an empty command or response file")
        _require_shell_literal_tokens(tokens, label)
        commands.append(tokens)
    require(commands, f"{label} contains no command")
    return commands


def _require_exact_driver(
    token: str, expected: Path, label: str,
) -> None:
    """Require CMake's retained driver spelling, not merely a symlink alias."""
    require(expected.is_absolute() and token == str(expected),
            f"{label} uses another command driver")
    try:
        resolved = Path(token).resolve(strict=True)
        expected_resolved = expected.resolve(strict=True)
    except (OSError, ValueError) as error:
        raise BuildProvenanceError(
            f"{label} has an invalid command driver") from error
    require(resolved == expected_resolved,
            f"{label} uses another command driver")


def _normalize_root_token(
    token: str, root_path: Path, replacement: str,
) -> str:
    """Replace exact root-path occurrences, never arbitrary substrings.

    CMake emits the root as a standalone token, after path-valued option
    prefixes such as -I, or after an '='.  Replacing an occurrence followed by
    punctuation (for example ``/tmp/build-suffix``) would hide a semantic argv
    difference between the candidate and clean rebuild.
    """
    root = str(root_path)
    require(root_path.is_absolute() and root and "\0" not in token and
            replacement,
            "compile/link build-root normalization input is invalid")
    attached_path_options = (
        "-I", "-L", "-F", "-isystem", "-iquote", "-idirafter",
        "-include", "-imacros",
    )
    result: list[str] = []
    cursor = 0
    while True:
        position = token.find(root, cursor)
        if position < 0:
            result.append(token[cursor:])
            break
        end = position + len(root)
        prefix = token[:position]
        left_boundary = (
            position == 0 or token[position - 1] in " \t'\"=,:"
        )
        if not left_boundary:
            word_start = max(
                token.rfind(" ", 0, position),
                token.rfind("\t", 0, position),
                token.rfind("'", 0, position),
                token.rfind('"', 0, position),
            ) + 1
            option_prefix = token[word_start:position]
            left_boundary = option_prefix in attached_path_options
        right_boundary = (
            end == len(token) or token[end] in "/ \t'\""
        )
        result.append(token[cursor:position])
        if left_boundary and right_boundary:
            result.append(replacement)
        else:
            result.append(root)
        cursor = end
    return "".join(result)


def _normalize_build_token(token: str, build: Path) -> str:
    return _normalize_root_token(token, build, "${BUILD_ROOT}")


def _normalize_build_argv(
    tokens: Sequence[str], build: Path, source: Path | None = None,
) -> list[str]:
    result = [_normalize_build_token(token, build) for token in tokens]
    if source is not None:
        result = [
            _normalize_root_token(token, source, "${SOURCE_ROOT}")
            for token in result
        ]
    return result


def _canonical_compile_argv(
    entry: Mapping[str, Any],
    tokens: Sequence[str],
    source: Path,
    build: Path,
    source_root: Path | None = None,
) -> list[str]:
    """Validate CMake's direct compile shape and retain its exact semantics."""
    require(tokens and all(
                token and not any(character in token for character in "\0\r\n")
                for token in tokens),
            f"compile command contains an invalid token for {source.name}")
    _require_shell_literal_tokens(
        tokens, f"compile command for {source.name}")
    compile_positions = [
        index for index, token in enumerate(tokens) if token == "-c"
    ]
    require(compile_positions == [len(tokens) - 2],
            f"compile command has a non-canonical source operand for "
            f"{source.name}")
    retained_source = entry.get("file")
    require(isinstance(retained_source, str) and
            tokens[-1] == retained_source and
            Path(tokens[-1]).resolve(strict=True) == source,
            "candidate compile command source metadata differs from argv")
    _compile_output(entry, tokens, build)
    return _normalize_build_argv(tokens, build, source_root)


def _archive_recipe_semantics(
    content: bytes,
    build: Path,
    archive: Path,
    archiver: Path,
    ranlib: Path,
    *,
    lexical_build_operands: bool = False,
) -> tuple[list[Path], list[list[str]]]:
    """Validate the exact two-command CMake static-archive recipe shape."""
    label = "candidate archive link recipe"
    commands = _recipe_commands(content, label)
    require(len(commands) == 2,
            f"{label} must contain exactly archive and index commands")
    archive_command, index_command = commands
    require(len(archive_command) >= 4,
            f"{label} archive command is incomplete")
    _require_exact_driver(archive_command[0], archiver, label)
    require(archive_command[1] == "qc",
            f"{label} uses non-canonical archive flags")
    resolve_operand = (
        _lexical_build_operand
        if lexical_build_operands else _resolve_build_operand)
    archive_operand = resolve_operand(
        build, build, archive_command[2])
    require(archive_operand == archive,
            f"{label} writes another archive")
    object_tokens = archive_command[3:]
    require(all(not token.startswith("-") and
                token.endswith((".o", ".obj"))
                for token in object_tokens),
            f"{label} contains a non-object or option-shaped archive operand")
    objects = [
        resolve_operand(build, build, token)
        for token in object_tokens
    ]
    require(objects and len(objects) == len(set(objects)),
            f"{label} object closure is empty or contains duplicates")

    require(len(index_command) == 2,
            f"{label} index command has extra operands")
    _require_exact_driver(index_command[0], ranlib, label)
    require(resolve_operand(build, build, index_command[1]) == archive,
            f"{label} indexes another archive")
    return objects, [
        _normalize_build_argv(command, build) for command in commands
    ]


def _executable_recipe_semantics(
    content: bytes,
    build: Path,
    executable: Path,
    archive: Path,
    compiler: Path,
    *,
    lexical_build_operands: bool = False,
) -> tuple[list[Path], list[str]]:
    """Validate the one canonical GNU/CMake production benchmark link argv."""
    label = "candidate executable link recipe"
    commands = _recipe_commands(content, label)
    require(len(commands) == 1,
            f"{label} must contain exactly one command")
    tokens = commands[0]
    _require_exact_driver(tokens[0], compiler, label)
    target = executable.name
    object_token = (
        f"CMakeFiles/{target}.dir/bench/leopard2/benchmark.cpp.o")
    canonical_prefix = [
        str(compiler), "-Wall", "-Wextra", "-fopenmp",
        "-O3", "-DNDEBUG", "-O3", object_token,
        "-o", target, archive.name,
    ]
    require(len(tokens) == len(canonical_prefix) + 2 and
            tokens[:len(canonical_prefix)] == canonical_prefix,
            f"{label} differs from the canonical CMake compiler, flags, "
            "benchmark object, output, or archive operands")

    def system_runtime(token: str, lexical_name: str) -> Path:
        operand = Path(token)
        require(operand.is_absolute() and operand.name == lexical_name,
                f"{label} contains a non-canonical system runtime operand")
        try:
            resolved = operand.resolve(strict=True)
            metadata = resolved.stat()
            allowed_roots = tuple(
                root.resolve(strict=True)
                for root in (Path("/usr/lib"), Path("/lib"))
                if root.exists())
        except OSError as error:
            raise BuildProvenanceError(
                f"{label} system runtime operand is invalid: {error}") \
                from error
        require(stat.S_ISREG(metadata.st_mode) and
                any(resolved.is_relative_to(root)
                    for root in allowed_roots),
                f"{label} system runtime operand is outside system "
                "library roots")
        return resolved

    system_runtime(tokens[-2], "libgomp.so")
    system_runtime(tokens[-1], "libpthread.a")
    resolve_operand = (
        _lexical_build_operand
        if lexical_build_operands else _resolve_build_operand)
    benchmark_object = resolve_operand(build, build, object_token)
    require(resolve_operand(build, build, target) == executable and
            resolve_operand(build, build, archive.name) == archive,
            f"{label} writes or links another artifact")
    return [benchmark_object], _normalize_build_argv(tokens, build)


def _compile_tokens(entry: Mapping[str, Any]) -> list[str]:
    has_arguments = "arguments" in entry
    has_command = "command" in entry
    require(has_arguments != has_command,
            "compile command must have exactly one argv representation")
    if has_arguments:
        value = entry["arguments"]
        require(isinstance(value, list) and value and
                all(isinstance(item, str) and item for item in value),
                "compile command arguments are malformed")
        return list(value)
    value = entry["command"]
    require(isinstance(value, str) and value,
            "compile command string is malformed")
    # Unlike the JSON ``arguments`` representation, this retained string is a
    # shell recipe.  Reject delimiters and shell operators in the raw text
    # before shlex can erase a command-separating newline or an escape.
    require(not any(character in value for character in "\0\r\n") and
            not any(character in value for character in SHELL_RECIPE_META),
            "compile command string contains shell control, substitution, "
            "or expansion syntax")
    try:
        return shlex.split(value, posix=True)
    except ValueError as error:
        raise BuildProvenanceError(
            f"cannot parse compile command: {error}") from error


def _compile_recipe_identity(
    entry: Mapping[str, Any], build: Path, source: Path | None = None,
) -> dict[str, str]:
    """Retain raw text plus a root-independent canonical argv rendering.

    CMake changes shell quoting around an otherwise identical ``-D`` token
    when a build root contains spaces.  Comparing raw text after pathname
    replacement therefore rejects a semantically identical clean replay.
    The command was already parsed with the strict no-shell-control grammar;
    bind its exact parsed argv and retain the original string for diagnostics.
    """
    if "arguments" in entry:
        return {"representation": "arguments"}
    command = entry.get("command")
    require(isinstance(command, str) and command,
            "compile command string is malformed")
    tokens = _compile_tokens({"command": command})
    normalized = _normalize_build_argv(tokens, build, source)
    return {
        "representation": "command",
        "command": command,
        "normalized_command": shlex.join(normalized),
    }


def _require_exact_compile_object_operand(
    entry: Mapping[str, Any], normalized_recipe_operand: str,
    build: Path, label: str,
) -> None:
    retained_output = entry.get("output")
    require(isinstance(retained_output, str) and retained_output and
            _normalize_build_token(retained_output, build) ==
            normalized_recipe_operand,
            f"{label} does not exactly name its compile output")


def _compile_output(
    entry: Mapping[str, Any], tokens: Sequence[str], build: Path,
) -> Path:
    directory_value = entry.get("directory")
    require(isinstance(directory_value, str) and directory_value,
            "compile command has no directory")
    directory = Path(directory_value).resolve(strict=True)
    require(directory == build,
            "compile command did not run at the exact candidate build root")
    output_positions = [index for index, token in enumerate(tokens)
                        if token == "-o"]
    require(len(output_positions) == 1 and
            output_positions[0] + 1 < len(tokens),
            "compile command has no unique output operand")
    output_value = tokens[output_positions[0] + 1]
    retained_output = entry.get("output")
    require(isinstance(retained_output, str) and
            retained_output == output_value,
            "compile command output metadata differs from its argv")
    return _resolve_build_operand(build, directory, output_value)


_GNU_CXX_DRIVER_BASENAME = re.compile(
    r"(?:g\+\+|(?:[A-Za-z0-9_.+]+-)+"
    r"(?:gnu|gnueabi(?:hf)?|eabi|elf|musl|mingw32)-g\+\+)"
    r"(?:-[0-9]+(?:\.[0-9]+)*)?")


def _resolved_compiler_is_gnu(
    compiler_path: Path | str | None,
) -> bool:
    """Classify only an unambiguous resolved GNU C++ driver pathname."""
    if compiler_path is None:
        return False
    path = Path(compiler_path)
    return path.is_absolute() and \
        _GNU_CXX_DRIVER_BASENAME.fullmatch(path.name) is not None


def _validate_compile_flags(
    tokens: Sequence[str], source: Path, *,
    source_root: Path | None = None,
    build_root: Path | None = None,
    cache: Mapping[str, Any] | None = None,
    library_sources: set[Path] | None = None,
    benchmark_source: bool | None = None,
    lexical_build_output: bool = False,
    compiler_path: Path | str | None = None,
) -> str:
    """Accept only the direct GCC/Clang argv emitted by this CMake graph.

    A denylist cannot bind the compilation closure: ``-B``, ``-specs``,
    plugins, wrappers, response files, sysroots, and driver-specific escape
    hatches can all select unretained code without spelling an ``-m`` flag.
    This grammar instead enumerates every admissible option and definition.
    """
    require(len(tokens) >= 5 and
            tokens[-4] == "-o" and tokens[-2] == "-c" and
            tokens[-1] == str(source) and
            tokens.count("-o") == 1 and tokens.count("-c") == 1,
            f"compile command has a non-canonical output/source tail for "
            f"{source.name}")
    require(all(
                token and not token.startswith(("@", "-Wl,@"))
                for token in tokens),
            f"compile command uses an empty or response-file operand for "
            f"{source.name}")
    if build_root is not None:
        output_operand = Path(tokens[-3])
        if lexical_build_output:
            lexical_output = Path(os.path.abspath(os.fspath(
                output_operand if output_operand.is_absolute() else
                build_root / output_operand)))
            output_is_contained = lexical_output.is_relative_to(build_root)
        else:
            output_is_contained = _resolve_build_operand(
                build_root, build_root, tokens[-3]).is_relative_to(build_root)
        require(build_root.is_absolute() and output_is_contained,
                f"compile command output escapes its build for {source.name}")
    cache = {} if cache is None else cache
    library_sources = set() if library_sources is None else library_sources

    name = source.name
    if name.startswith("Leopard2BackendAVX2"):
        require(name in PRODUCTION_AVX2_LIBRARY_SOURCES,
                "compile source is not an allowlisted production AVX2 "
                f"translation unit: {name}")
    p32_source_name = "Leopard2LowP32B64AVX2.cpp"
    enhanced_backend = (
        name.startswith((
            "Leopard2BackendSSSE3", "Leopard2BackendAVX2",
            "Leopard2BackendGFNI")) or
        name == p32_source_name)
    if name.startswith("Leopard2BackendGFNI"):
        profile = "avx2-gfni-no-avx512"
        expected_target_flags = (
            "-mavx2", "-mgfni", "-mno-avx512f")
    elif name.startswith("Leopard2BackendAVX2"):
        profile = "avx2-no-avx512"
        expected_target_flags = ("-mavx2", "-mno-avx512f")
    elif name == p32_source_name:
        profile = "avx2-no-avx512"
        expected_target_flags = ("-mavx2", "-mno-avx512f")
    elif name == "Leopard2BackendSSSE3.cpp":
        profile = "ssse3-no-avx"
        expected_target_flags = ("-mssse3", "-mno-avx")
    else:
        profile = "portable-core"
        expected_target_flags = ()

    expected_options = Counter({
        "-Wall": 1,
        "-Wextra": 1,
        "-fopenmp": 1 if enhanced_backend else 2,
        "-O3": 2,
        "-std=gnu++11": 1,
    })
    expected_options.update(expected_target_flags)
    if (name == "Leopard2BackendAVX2T8K8B1024.cpp" and
            cache.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") in {
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            } and
            _resolved_compiler_is_gnu(compiler_path)):
        expected_options["-flive-range-shrinkage"] = 1
    if ((name.startswith(("Leopard2BackendAVX2",
                          "Leopard2BackendGFNI")) or
         name == p32_source_name) and
            _cmake_true(cache.get("LEO2_FLAG_FALIGN_FUNCTIONS_64"))):
        expected_options["-falign-functions=64"] = 1
    if source_root is not None:
        require(source_root.is_absolute(),
                "compile source root is not absolute")
        expected_options[f"-I{source_root}"] = 1

    definitions = [
        token for token in tokens[1:-4] if token.startswith("-D")]
    observed_options = Counter(
        token for token in tokens[1:-4] if not token.startswith("-D"))
    require(observed_options == expected_options,
            f"{profile} source has a non-canonical or indirect compile "
            f"option: {name}: observed={sorted(observed_options.elements())!r} "
            f"expected={sorted(expected_options.elements())!r}")
    require(len(definitions) == len(set(definitions)),
            f"{profile} source has duplicate compile definitions: {name}")

    expected_definitions: set[str] = {"-DNDEBUG"}
    if source_root is not None:
        current_configuration = (
            cache.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") in {
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            })
        if current_configuration:
            expected_definitions.update((
                "-DLEO2_EXPERIMENT_CAUCHY_LOG_REUSE=1",
                "-DLEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT=1",
            ))
            if name in {
                    "leopard2.cpp", "Leopard2Backend.cpp",
                    "Leopard2BackendAVX2.cpp"}:
                expected_definitions.add(
                    "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=1")
        available_names = {path.name for path in library_sources}
        t8_k8_b1024_source_name = "Leopard2BackendAVX2T8K8B1024.cpp"
        t8_k8_b1024_available = (
            name == t8_k8_b1024_source_name or
            t8_k8_b1024_source_name in available_names)
        t16_source_name = "Leopard2BackendAVX2T16B64.cpp"
        t16_available = (
            name == t16_source_name or t16_source_name in available_names)
        if "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED" in cache:
            t16_generated = _cmake_true(
                cache.get("LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED"))
        else:
            # Reduced historical replay caches predate retention of this
            # default-on selector.  Archive/source membership is the only
            # authenticated signal available in that schema.
            t16_generated = t16_available
        t32_source_name = "Leopard2BackendAVX2T32B256.cpp"
        t32_available = (
            name == t32_source_name or t32_source_name in available_names)
        t32_generated = _cmake_true(
            cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED"))
        if "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK" in cache:
            t32_two_block = _cmake_true(
                cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK"))
        else:
            # Current replay records predate retention of this default-on
            # selector.  With the generated selector forced off, presence of
            # the isolated object in the audited archive proves that the
            # two-block selector supplied it.
            t32_two_block = t32_available and not t32_generated
        p32_available = (
            name == p32_source_name or p32_source_name in available_names)
        t16_q2_source_name = "Leopard2BackendAVX2T16Q2.cpp"
        t16_q2_enabled = _cmake_true(
            cache.get("LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED"))
        if "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL" in cache:
            low_p32_terminal = _cmake_true(
                cache.get("LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL"))
        else:
            # The selector is default-on but is not retained by the current
            # replay schema.  Archive membership is the authenticated signal
            # that it was active for the production core compilation.
            low_p32_terminal = p32_available
        if source in library_sources and not enhanced_backend:
            expected_definitions.update((
                "-DLEO2_DISABLE_AVX2_CODEGEN=1",
                "-DLEO2_DISABLE_SSSE3_CODEGEN=1",
            ))
            if current_configuration:
                expected_definitions.update((
                    "-DLEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING=1",
                    "-DLEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING=1",
                    "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1",
                ))
            for backend_name, definition in (
                ("Leopard2BackendSSSE3.cpp",
                 "-DLEO2_HAVE_SSSE3_BACKEND=1"),
                ("Leopard2BackendAVX2.cpp",
                 "-DLEO2_HAVE_AVX2_BACKEND=1"),
                ("Leopard2BackendGFNI.cpp",
                 "-DLEO2_HAVE_GFNI_BACKEND=1"),
            ):
                if backend_name in available_names:
                    expected_definitions.add(definition)
            if t8_k8_b1024_available:
                expected_definitions.add(
                    "-DLEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1")
            if name == "leopard2.cpp" and t16_generated:
                expected_definitions.add(
                    "-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=1")
            if name == "leopard2.cpp" and t32_generated:
                expected_definitions.add(
                    "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1")
            if name == "leopard2.cpp" and t32_two_block:
                expected_definitions.add(
                    "-DLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=1")
            if low_p32_terminal:
                expected_definitions.add(
                    "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1")
            if cache.get("LEO2_BACKEND_VARIANT") == "avx2":
                expected_definitions.add("-DLEO2_BACKEND_FORCE_AVX2=1")
            if (name == "leopard2.cpp" and
                    cache.get("LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE")
                    in {"1", "2"}):
                expected_definitions.add(
                    "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=" +
                    str(cache["LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"]))
            configuration_schema = cache.get(
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA")
            if (name == "leopard2.cpp" and configuration_schema in {
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA}):
                dual_direct = cache.get(
                    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT")
                require(dual_direct in {"ON", "OFF"},
                        "current GF8 small dual-direct selector is not "
                        "exactly ON or OFF")
                expected_definitions.add(
                    "-DLEO2_ENABLE_GF8_SMALL_DUAL_DIRECT=" +
                    ("1" if dual_direct == "ON" else "0"))
            if (name == "leopard2.cpp" and configuration_schema in {
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA}):
                for selector in (
                        "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                        "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
                    value = cache.get(selector)
                    require(value in {"ON", "OFF"},
                            f"current {selector} is not exactly ON or OFF")
                    expected_definitions.add(
                        f"-D{selector}=" + ("1" if value == "ON" else "0"))
        elif name.startswith("Leopard2BackendGFNI"):
            expected_definitions.update((
                "-DLEO2_HAVE_AVX2_BACKEND=1",
                "-DLEO2_HAVE_GFNI_BACKEND=1",
            ))
        elif name == "Leopard2BackendAVX2T2K4.cpp":
            expected_definitions.add("-DLEO2_HAVE_AVX2_BACKEND=1")
        elif name == t8_k8_b1024_source_name:
            expected_definitions.update((
                "-DLEO2_HAVE_AVX2_BACKEND=1",
                "-DLEO2_HAVE_AVX2_T8_K8_B1024_DIRECT=1",
            ))
        elif name == t16_source_name:
            require(t16_generated,
                    "T16/B64 source contradicts its disabled selector")
            expected_definitions.update((
                "-DLEO2_HAVE_AVX2_BACKEND=1",
                "-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=1",
            ))
        elif name == t16_q2_source_name:
            require(t16_q2_enabled,
                    "T16/Q2 source contradicts its disabled selector")
            expected_definitions.add(
                "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1")
        elif name == t32_source_name:
            require(t32_generated or t32_two_block,
                    "T32/B256 source contradicts its disabled selectors")
            expected_definitions.add("-DLEO2_HAVE_AVX2_BACKEND=1")
            if t32_generated:
                expected_definitions.update((
                    "-DLEO2_EXPERIMENT_HIGH_T32_B256_GENERATED=1",
                    "-DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=" +
                    ("1" if _cmake_true(cache.get(
                        "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED"))
                     else "0"),
                ))
            if t32_two_block:
                expected_definitions.update((
                    "-DLEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK=1",
                    "-DLEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK=" +
                    ("1" if _cmake_true(cache.get(
                        "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK"))
                     else "0"),
                ))
        elif name == p32_source_name:
            require(low_p32_terminal,
                    "low P32/B64 source contradicts its disabled selector")
            expected_definitions.update((
                "-DLEO2_HAVE_AVX2_BACKEND=1",
                "-DLEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1",
            ))
            if not _cmake_true(cache.get("LEOPARD_ENABLE_GF16")):
                expected_definitions.add("-DNO_LEO_HAS_FF16=1")
        elif name.startswith("Leopard2BackendAVX2"):
            expected_definitions.add(
                "-DLEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING=1")
        elif (
                benchmark_source is True or
                (benchmark_source is None and source == (
                    source_root / "bench/leopard2/benchmark.cpp").resolve(
                        strict=True))):
            expected_header = (
                build_root /
                "generated/leopard2-benchmark-attestation/"
                "leopard2_benchmark_source_attestation.h"
                if build_root is not None else None)
            require(expected_header is not None,
                    "benchmark compile grammar lacks a build root")
            expected_definitions.update((
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION=1",
                '-DLEO2_BENCHMARK_BUILD_TYPE="Release"',
                "-DLEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER=" +
                f'"{expected_header}"',
            ))
            hashes = [
                token for token in definitions
                if token.startswith(
                    "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=")]
            configuration_digest = cache.get(
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256")
            expected_hash = (
                "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=" +
                f'"{configuration_digest}"')
            require(
                isinstance(configuration_digest, str) and
                re.fullmatch(r"[0-9a-f]{64}", configuration_digest)
                    is not None and
                hashes == [expected_hash],
                "benchmark compile configuration digest is not bound to "
                "the retained CMake cache")
            expected_definitions.add(expected_hash)

        is_benchmark_source = (
            benchmark_source is True or
            (benchmark_source is None and source == (
                source_root / "bench/leopard2/benchmark.cpp").resolve(
                    strict=True)))
        if (t16_q2_enabled and
                ((source in library_sources and not enhanced_backend) or
                 name in {
                    "Leopard2BackendAVX2.cpp",
                    "Leopard2BackendAVX2Xor.cpp",
                    t16_q2_source_name,
                 } or is_benchmark_source)):
            expected_definitions.add(
                "-DLEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED=1")

    require(set(definitions) == expected_definitions,
            f"{profile} source has non-canonical compile definitions: "
            f"{name}: observed={sorted(definitions)!r} "
            f"expected={sorted(expected_definitions)!r}")
    return profile


def _compiler_subtool_identities(
    compiler: Path, language: str,
    inherited_descriptors: Sequence[int],
) -> list[dict[str, Any]]:
    """Bind GCC/Clang external frontend, assembler, and linker helpers."""
    require(language in ("c", "c++"), "compiler subtool language is invalid")
    names = (
        ("cc1", "as", "collect2", "ld")
        if language == "c" else
        ("cc1plus", "as", "collect2", "ld"))
    discovered: dict[Path, str] = {}
    with _RetainedFileSnapshot(
            compiler, f"candidate {language} compiler subtool driver"
            ) as driver:
        for name in names:
            raw = _run(
                (str(compiler), f"-print-prog-name={name}"),
                f"candidate {language} compiler {name} lookup",
                maximum_bytes=MAX_METADATA_BYTES,
                inherited_descriptors=inherited_descriptors,
                executable_descriptor=driver.executable_descriptor)
            try:
                value = raw.decode("utf-8", errors="strict").strip()
            except UnicodeDecodeError as error:
                raise BuildProvenanceError(
                    f"candidate {language} compiler {name} path is not "
                    "strict UTF-8") from error
            _require_safe_unicode(
                value, f"candidate {language} compiler {name} path")
            require(value and "\n" not in value and "\r" not in value,
                    f"candidate {language} compiler {name} lookup is "
                    "malformed")
            candidate = Path(value)
            if not candidate.is_absolute():
                located = shutil.which(value, path=GIT_ENVIRONMENT["PATH"])
                if located is None:
                    # Integrated frontends/assemblers legitimately report an
                    # unresolved lexical name.  They have no pathname to race.
                    continue
                candidate = Path(located)
            try:
                resolved = candidate.resolve(strict=True)
            except OSError as error:
                raise BuildProvenanceError(
                    f"candidate {language} compiler {name} subtool is "
                    f"invalid: {error}") from error
            require(resolved.is_file(),
                    f"candidate {language} compiler {name} subtool is not "
                    "a file")
            discovered.setdefault(resolved, name)
    records = []
    for path, name in sorted(
            discovered.items(), key=lambda item: str(item[0])):
        records.append({
            "language": language,
            "role": name,
            "identity": file_identity(
                path, f"candidate {language} compiler {name} subtool"),
        })
    require(records,
            f"candidate {language} compiler exposed no bindable subtool")
    return records


def _tracked_files(
    source: Path, *, inherited_descriptors: Sequence[int] = (),
) -> set[Path]:
    git = Path("/usr/bin/git").resolve(strict=True)
    raw = _run((str(git), "-C", str(source), "ls-files", "-z"),
               "candidate tracked-file inventory", maximum_bytes=32 << 20,
               inherited_descriptors=inherited_descriptors)
    result = set()
    for item in raw.split(b"\0"):
        if not item:
            continue
        try:
            relative = item.decode("utf-8", errors="strict")
        except UnicodeDecodeError as error:
            raise BuildProvenanceError(
                "candidate tracked path is not strict UTF-8") from error
        path = (source / relative).resolve(strict=True)
        require(path.is_relative_to(source),
                "candidate tracked path escapes the source root")
        result.add(path)
    require(result, "candidate tracked-file inventory is empty")
    return result


def _expected_library_sources(
    source: Path, cache: Mapping[str, str], tracked: set[Path],
    *,
    expected_configuration_schema: str =
        BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
) -> set[Path]:
    require(expected_configuration_schema in {
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            },
            "candidate source-closure configuration schema is unsupported")
    names = set(CORE_LIBRARY_SOURCES)
    names.add("LeopardFF8.cpp")
    if _cmake_true(cache.get("LEOPARD_ENABLE_GF16")):
        names.add("LeopardFF16.cpp")
    names.add("Leopard2BackendSSSE3.cpp")
    if expected_configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        names.update(UNCONDITIONAL_AVX2_LIBRARY_SOURCES)
    elif (expected_configuration_schema ==
          BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6):
        names.update(UNCONDITIONAL_AVX2_LIBRARY_SOURCES_V6)
    else:
        names.update(UNCONDITIONAL_AVX2_LIBRARY_SOURCES_V5)
    if _cmake_true(cache.get("LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED")):
        names.add("Leopard2BackendAVX2T16B64.cpp")
    if _cmake_true(cache.get("LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED")):
        names.add("Leopard2BackendAVX2T16Q2.cpp")
    if (_cmake_true(
            cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED")) or
            _cmake_true(
                cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK"))):
        names.add("Leopard2BackendAVX2T32B256.cpp")
    if _cmake_true(cache.get("LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL")):
        names.add("Leopard2LowP32B64AVX2.cpp")
    # The GFNI member joins the archive whenever its compiler probe passes.
    # It contains no AVX-512 encoding, so a pure-AVX2 candidate retains it.
    if _cmake_true(cache.get("LEO2_FLAG_MGFNI")):
        names.add("Leopard2BackendGFNI.cpp")
    expected = {(source / name).resolve(strict=True) for name in names}
    require(expected.issubset(tracked),
            "candidate build expects a library source not tracked at HEAD")
    return expected


def _require_compile_driver(
    tokens: Sequence[str],
    source: Path,
    c_compiler: Path,
    cxx_compiler: Path,
) -> None:
    """Bind each configured compile entry to its exact language driver.

    CMake emits entries for the C ABI tests when LEO2_BUILD_TESTS is enabled,
    even though those objects do not enter the production archive.  They are
    valid metadata only when a .c source uses the cached C compiler; every
    supported C++ suffix must use the cached C++ compiler.  Looking at argv[0]
    also rejects unbound compiler launchers.
    """
    suffix = source.suffix
    if suffix == ".c":
        expected = c_compiler
    elif suffix in (".C", ".cc", ".cpp", ".cxx"):
        expected = cxx_compiler
    else:
        raise BuildProvenanceError(
            f"candidate compile command uses unsupported source language: "
            f"{source.name}")
    require(bool(tokens),
            "candidate compile command omits its language compiler")
    require(not any(token.startswith("@") for token in tokens),
            "candidate compile command uses an unbound response file")
    _require_exact_driver(
        tokens[0], expected, "candidate compile command")


def _validate_candidate_required_cache(
    cache: Mapping[str, Any],
    *,
    expected_configuration_schema: str =
        BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
) -> dict[str, str]:
    """Require every production-evidence cache value that cannot vary."""
    require(isinstance(cache, Mapping),
            "candidate CMake cache is malformed")
    required_exact_v5 = {
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_GENERATOR": "Unix Makefiles",
        "ENABLE_OPENMP": "ON",
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
        "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_ENABLE_CUDA": "OFF",
        "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
        "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED": "OFF",
        "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE": "ON",
        "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
        "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "ON",
        "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
        "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING": "ON",
        "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING": "ON",
        "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
        "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "OFF",
        "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT": "ON",
        "LEOPARD_ENABLE_GF8": "ON",
        "LEOPARD_ENABLE_GF16": "ON",
    }
    require(expected_configuration_schema in {
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            },
            "candidate CMake cache configuration schema is unsupported")
    required_exact = dict(required_exact_v5)
    if expected_configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        required_exact.update({
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                expected_configuration_schema,
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK": "OFF",
            "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED": "ON",
            "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK": "ON",
            "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL": "ON",
        })
    if expected_configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        required_exact["LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT"] = "ON"
    else:
        require("LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT" not in cache,
                "historical candidate CMake cache unexpectedly retains "
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT")
    if expected_configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        required_exact.update({
            "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS": "ON",
            "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK":
                "ON" if expected_configuration_schema ==
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA else "OFF",
        })
    else:
        for selector in (
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"):
            require(selector not in cache,
                    "historical candidate CMake cache unexpectedly retains " +
                    selector)
    if expected_configuration_schema == BENCHMARK_BUILD_CONFIGURATION_SCHEMA:
        required_exact["LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED"] = "ON"
    else:
        require("LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED" not in cache,
                "historical candidate CMake cache unexpectedly retains "
                "LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED")
    for name, expected in required_exact.items():
        require(cache.get(name) == expected,
                f"candidate CMake cache {name} is {cache.get(name)!r}, "
                f"expected {expected!r}")
    configuration_digest = cache.get(
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256")
    require(isinstance(configuration_digest, str) and
            re.fullmatch(r"[0-9a-f]{64}", configuration_digest) is not None,
            "candidate CMake cache effective configuration digest is "
            "malformed")
    return {
        **required_exact,
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
            configuration_digest,
    }


def _reproducible_replay_contract(
    candidate: Mapping[str, Any],
) -> dict[str, str]:
    """Select one explicit, non-downgradable replay evidence contract.

    The historical v1 closure remains readable only with the exact metadata
    matrix it originally emitted.  Current candidates cannot infer an older
    transport from a missing selector or a forged proof schema.
    """
    require(isinstance(candidate, Mapping),
            "reproducible-build closure is malformed")
    cache = candidate.get("validated_cache")
    require(isinstance(cache, Mapping),
            "reproducible-build closure lacks retained CMake semantics")
    configuration_digest = cache.get(
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256")
    require(isinstance(configuration_digest, str) and
            re.fullmatch(r"[0-9a-f]{64}", configuration_digest) is not None,
            "reproducible-build closure configuration digest is malformed")
    closure_schema = candidate.get("schema")
    if closure_schema == PRODUCTION_BUILD_CLOSURE_SCHEMA:
        configuration_schema = cache.get(
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA")
        if configuration_schema in {
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
            require(
                cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT") ==
                    "ON" and
                cache.get(
                    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT") ==
                    "ON" and
                cache.get("LEO2_EXPERIMENT_CAUCHY_LOG_REUSE") == "ON" and
                cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED") ==
                    "ON" and
                cache.get(
                    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED") ==
                    "OFF" and
                cache.get("LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED") ==
                    "ON" and
                cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK") ==
                    "ON" and
                cache.get(
                    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK") ==
                    "OFF" and
                cache.get("LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL") ==
                    "ON" and
                ((configuration_schema in {
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA} and
                  cache.get("LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT") == "ON" and
                  cache.get("LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS") ==
                    "ON" and
                  ((configuration_schema ==
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA and
                    cache.get(
                        "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK") ==
                        "ON" and
                    cache.get(
                        "LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED") == "ON") or
                   (configuration_schema ==
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9 and
                    cache.get(
                        "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK") ==
                        "OFF" and
                    "LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED" not in cache))) or
                 (configuration_schema ==
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8 and
                  cache.get("LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT") == "ON" and
                  "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS" not in cache and
                  "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK" not in
                    cache) or
                 (configuration_schema in {
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7} and
                  "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT" not in cache and
                  "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS" not in cache and
                  "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK" not in
                    cache)),
                "v6/v7/v8/v9/current reproducible-build closure configuration "
                "contract differs")
            return {
                "closure_schema": PRODUCTION_BUILD_CLOSURE_SCHEMA,
                "configuration_schema": configuration_schema,
                "proof_schema": REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                "recipe_schema": CANONICAL_REPLAY_RECIPE_SCHEMA,
                "invocation_schema": REPLAY_INVOCATION_SCHEMA,
            }
        if configuration_schema == BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5:
            require(
                cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT") ==
                    "ON" and
                cache.get(
                    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT") ==
                    "ON" and
                cache.get("LEO2_EXPERIMENT_CAUCHY_LOG_REUSE") == "ON" and
                cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED") ==
                    "OFF" and
                cache.get(
                    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED") ==
                    "OFF" and
                all(name not in cache for name in (
                    "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED",
                    "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK",
                    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK",
                    "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL",
                    "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT",
                    "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS",
                    "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK")),
                "v5 reproducible-build closure configuration contract "
                "differs")
            return {
                "closure_schema": PRODUCTION_BUILD_CLOSURE_SCHEMA,
                "configuration_schema":
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
                "proof_schema": REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                "recipe_schema": CANONICAL_REPLAY_RECIPE_SCHEMA,
                "invocation_schema": REPLAY_INVOCATION_SCHEMA,
            }
        if configuration_schema == BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4:
            require(
                cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT") ==
                    "ON" and
                cache.get(
                    "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT") ==
                    "ON" and
                cache.get("LEO2_EXPERIMENT_CAUCHY_LOG_REUSE") == "ON" and
                "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED" not in cache and
                "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED" not in
                    cache and
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT" not in cache and
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS" not in cache and
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK" not in cache,
                "v4 reproducible-build closure configuration contract "
                "differs")
            return {
                "closure_schema": PRODUCTION_BUILD_CLOSURE_SCHEMA,
                "configuration_schema":
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4,
                "proof_schema": REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                "recipe_schema": CANONICAL_REPLAY_RECIPE_SCHEMA,
                "invocation_schema": REPLAY_INVOCATION_SCHEMA,
            }
        if configuration_schema == BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3:
            require(
                cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT") ==
                    "OFF" and
                "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT" not in
                    cache and
                "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE" not in cache and
                "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT" not in cache and
                "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS" not in cache and
                "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK" not in cache,
                "v3 reproducible-build closure configuration contract "
                "differs")
            return {
                "closure_schema": PRODUCTION_BUILD_CLOSURE_SCHEMA,
                "configuration_schema":
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3,
                "proof_schema": REPRODUCIBLE_BUILD_PROOF_SCHEMA,
                "recipe_schema": CANONICAL_REPLAY_RECIPE_SCHEMA,
                "invocation_schema": REPLAY_INVOCATION_SCHEMA,
            }
        raise BuildProvenanceError(
            "current reproducible-build closure configuration schema is "
            "unsupported")
    if closure_schema == PRODUCTION_BUILD_CLOSURE_SCHEMA_V1:
        require(
            cache.get("LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") ==
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V2 and
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT" not in cache and
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT" not in cache and
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE" not in cache and
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT" not in cache and
            "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS" not in cache and
            "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK" not in cache,
            "historical reproducible-build closure configuration contract "
            "differs")
        return {
            "closure_schema": PRODUCTION_BUILD_CLOSURE_SCHEMA_V1,
            "configuration_schema": BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V2,
            "proof_schema": REPRODUCIBLE_BUILD_PROOF_SCHEMA_V2,
            "recipe_schema": LEGACY_REPLAY_RECIPE_SCHEMA,
            "invocation_schema": REPLAY_INVOCATION_SCHEMA_V1,
        }
    raise BuildProvenanceError(
        "reproducible-build closure schema is unsupported")


def _require_reproducible_replay_artifact_contract(
    candidate: Mapping[str, Any], proof_schema: Any, recipe_schema: Any,
) -> dict[str, str]:
    contract = _reproducible_replay_contract(candidate)
    require(
        proof_schema == contract["proof_schema"] and
        recipe_schema == contract["recipe_schema"],
        "reproducible-build proof/recipe transport is incompatible with "
        "its closure contract")
    return contract


def candidate_build_provenance(
    build_root: Path | str,
    source_root: Path | str,
    executable: Path | str,
    executable_target: str,
    *,
    inherited_descriptors: Sequence[int] = (),
    tracked_source_manifest: Mapping[str, Any] | None = None,
    logical_source_root: Path | str | None = None,
    _expected_configuration_schema: str =
        BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
    sealed_artifacts: Mapping[
        Path | str, _ReplayArtifactSink] | None = None,
) -> dict[str, Any]:
    """Capture and validate one pure-AVX2 production benchmark build."""
    build = Path(build_root).resolve(strict=True)
    source = Path(source_root).resolve(strict=True)
    logical_source = (
        source if logical_source_root is None else
        Path(logical_source_root).resolve(strict=True))
    require(build.is_dir(), "candidate build root is not a directory")
    require(source.is_dir(), "candidate source root is not a directory")
    require(executable_target in {"bench_leopard2", "bench_leopard2_allk"},
            "unsupported candidate benchmark target")
    artifact_overrides: dict[Path, _ReplayArtifactSink] = {}
    if sealed_artifacts is not None:
        require(isinstance(sealed_artifacts, Mapping) and sealed_artifacts,
                "candidate sealed-artifact map is malformed")
        for raw_path, artifact in sealed_artifacts.items():
            path = Path(os.path.abspath(os.fspath(raw_path)))
            require(
                path.is_relative_to(build) and
                isinstance(artifact, _ReplayArtifactSink) and
                artifact.logical_path == path and
                path not in artifact_overrides,
                "candidate sealed-artifact identity is malformed")
            artifact.verify()
            artifact_overrides[path] = artifact

    def artifact_identity(path: Path, label: str) -> dict[str, Any]:
        lexical = Path(os.path.abspath(os.fspath(path)))
        artifact = artifact_overrides.get(lexical)
        if artifact is None:
            return file_identity(path, label)
        artifact.verify()
        return dict(artifact.identity)
    source_manifest = (
        _capture_tracked_source_tree(
            source, inherited_descriptors=inherited_descriptors)
        if tracked_source_manifest is None else
        _verify_tracked_source_manifest(source, tracked_source_manifest)
    )
    manifest_files = source_manifest.get("files")
    require(isinstance(manifest_files, list) and manifest_files,
            "candidate tracked source manifest is empty")
    manifest_by_path = {
        record["path"]: record for record in manifest_files
        if isinstance(record, dict) and isinstance(record.get("path"), str)
    }
    require(len(manifest_by_path) == len(manifest_files),
            "candidate tracked source manifest is malformed or duplicated")

    expected_executable = (build / executable_target).resolve(strict=True)
    requested_executable = Path(executable).resolve(strict=True)
    require(requested_executable == expected_executable,
            "candidate executable is not the declared CMake build target")
    archive = (build / "libleopard.a").resolve(strict=True)
    cache_path = build / "CMakeCache.txt"
    commands_path = build / "compile_commands.json"
    archive_recipe_path = build / "CMakeFiles/leopard.dir/link.txt"
    executable_recipe_path = \
        build / f"CMakeFiles/{executable_target}.dir/link.txt"

    cache_identity, cache_bytes = file_snapshot(
        cache_path, "candidate CMake cache", maximum_bytes=MAX_METADATA_BYTES)
    cache = parse_cmake_cache(cache_bytes)
    required_exact = _validate_candidate_required_cache(
        cache,
        expected_configuration_schema=_expected_configuration_schema)
    benchmark_git_value = cache.get("LEO2_BENCHMARK_GIT_EXECUTABLE")
    locator_git_value = cache.get("LEO2_LOCATOR_GIT_EXECUTABLE")
    require(isinstance(benchmark_git_value, str) and benchmark_git_value and
            isinstance(locator_git_value, str) and locator_git_value and
            Path(benchmark_git_value).is_absolute() and
            Path(locator_git_value).is_absolute() and
            benchmark_git_value == locator_git_value and
            Path(benchmark_git_value).resolve(strict=True) ==
            Path(locator_git_value).resolve(strict=True),
            "candidate CMake cache does not bind one exact Git executable")
    benchmark_git_identity = file_identity(
        benchmark_git_value, "candidate benchmark Git executable")
    if executable_target == "bench_leopard2_allk":
        require(cache.get("LEO2_BUILD_ALLK_DIAGNOSTIC") == "ON",
                "all-K candidate target was not explicitly enabled")
    variant = cache.get("LEO2_BACKEND_VARIANT")
    require(variant in {"auto", "avx2"},
            "candidate backend variant is not production auto or AVX2")
    small_direct_mode = cache.get("LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE")
    require(small_direct_mode in {"0", "1", "2"},
            "candidate small-direct selector is not exactly 0, 1, or 2")
    dual_direct = cache.get("LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT")
    if _expected_configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        require(dual_direct in {"ON", "OFF"},
                "candidate GF8 small dual-direct selector is not exactly "
                "ON or OFF")
    else:
        require("LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT" not in cache,
                "historical candidate unexpectedly retains the GF8 small "
                "dual-direct selector")
    require(_cmake_true(cache.get("LEO2_FLAG_MAVX2")) and
            _cmake_true(cache.get("LEO2_FLAG_MNO_AVX512F")),
            "candidate compiler lacks the required pure-AVX2 flags")
    for name in ("LEO2_FLAG_MAVX512F", "LEO2_FLAG_MAVX512BW",
                 "LEO2_FLAG_MAVX512VL"):
        require(_cmake_false(cache.get(name)),
                f"candidate pure-AVX2 build did not disable {name}")
    home = Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve(strict=True)
    require(home == source,
            "candidate CMake cache points at another source root")
    global_flags = " ".join((cache.get("CMAKE_CXX_FLAGS", ""),
                             cache.get("CMAKE_CXX_FLAGS_RELEASE", "")))
    require(cache.get("CMAKE_CXX_FLAGS") == "" and
            cache.get("CMAKE_CXX_FLAGS_RELEASE") == "-O3 -DNDEBUG",
            "candidate CMake flags differ from the canonical Release profile")
    for name in ("CMAKE_C_COMPILER_LAUNCHER", "CMAKE_C_COMPILER_ARG1",
                 "CMAKE_C_COMPILER_TARGET",
                 "CMAKE_CXX_COMPILER_LAUNCHER", "CMAKE_CXX_COMPILER_ARG1",
                 "CMAKE_CXX_COMPILER_TARGET", "CMAKE_SYSROOT",
                 "CMAKE_TOOLCHAIN_FILE"):
        require(not cache.get(name),
                f"candidate CMake cache uses unsupported {name}")
    require(_cmake_false(cache.get("CMAKE_POSITION_INDEPENDENT_CODE")),
            "candidate CMake cache changes position-independent code policy")
    require(not any(flag in global_flags for flag in (
                "-march=native", "-mtune=native", "-mavx", "-mssse3",
                "-mavx512")),
            "candidate CMake cache raises the global ISA floor")

    compiler_invocation = cache.get("CMAKE_CXX_COMPILER", "")
    require(bool(compiler_invocation),
            "candidate CMake cache has no C++ compiler")
    compiler = Path(compiler_invocation)
    require(compiler.is_absolute(),
            "candidate C++ compiler invocation is not absolute")
    with _RetainedFileSnapshot(
            compiler, "candidate C++ compiler") as compiler_snapshot:
        compiler_identity = compiler_snapshot.identity
        compiler_version = _run(
            (str(compiler), "--version"), "candidate compiler version",
            inherited_descriptors=inherited_descriptors,
            executable_descriptor=compiler_snapshot.executable_descriptor)
    c_compiler_invocation = cache.get("CMAKE_C_COMPILER", "")
    require(bool(c_compiler_invocation),
            "candidate CMake cache has no C compiler")
    c_compiler = Path(c_compiler_invocation)
    require(c_compiler.is_absolute(),
            "candidate C compiler invocation is not absolute")
    c_compiler_identity = file_identity(
        c_compiler, "candidate C compiler")
    make_invocation = cache.get("CMAKE_MAKE_PROGRAM", "")
    linker_invocation = cache.get("CMAKE_LINKER", "")
    require(bool(make_invocation) and bool(linker_invocation),
            "candidate CMake cache lacks its make program or linker")
    make_program = Path(make_invocation)
    linker = Path(linker_invocation)
    require(make_program.is_absolute() and linker.is_absolute(),
            "candidate make program or linker invocation is not absolute")
    make_program_identity = file_identity(
        make_program, "candidate make program")
    linker_identity = file_identity(linker, "candidate linker")
    compiler_subtools = (
        _compiler_subtool_identities(
            c_compiler, "c", inherited_descriptors) +
        _compiler_subtool_identities(
            compiler, "c++", inherited_descriptors)
    )

    commands_identity, commands_bytes = file_snapshot(
        commands_path, "candidate compile commands",
        maximum_bytes=MAX_METADATA_BYTES)
    entries = _strict_json_loads(
        commands_bytes, "candidate compile_commands.json")
    require(isinstance(entries, list) and entries,
            "candidate compile command closure is empty")
    by_output: dict[
        Path, tuple[
            Mapping[str, Any], list[str], list[str], Path, dict[str, str],
        ]
    ] = {}
    for entry in entries:
        require(isinstance(entry, dict) and
                isinstance(entry.get("file"), str),
                "candidate compile command entry is malformed")
        tokens = _compile_tokens(entry)
        source_operand = Path(entry["file"]).resolve(strict=True)
        _require_compile_driver(
            tokens, source_operand, c_compiler, compiler)
        normalized_arguments = _canonical_compile_argv(
            entry, tokens, source_operand, build, source)
        recipe_identity = _compile_recipe_identity(entry, build, source)
        output = _compile_output(entry, tokens, build)
        require(output not in by_output,
                f"duplicate candidate compile output: {output}")
        by_output[output] = (
            entry, tokens, normalized_arguments, source_operand,
            recipe_identity)

    archive_recipe_identity, archive_recipe_bytes = file_snapshot(
        archive_recipe_path, "candidate archive link recipe",
        maximum_bytes=MAX_METADATA_BYTES)
    executable_recipe_identity, executable_recipe_bytes = file_snapshot(
        executable_recipe_path, "candidate executable link recipe",
        maximum_bytes=MAX_METADATA_BYTES)

    ar_invocation = cache.get("CMAKE_AR", "")
    require(bool(ar_invocation), "candidate CMake cache has no archiver")
    archiver = Path(ar_invocation)
    require(archiver.is_absolute(),
            "candidate archiver invocation is not absolute")
    ranlib_invocation = cache.get("CMAKE_RANLIB", "")
    require(bool(ranlib_invocation), "candidate CMake cache has no ranlib")
    ranlib = Path(ranlib_invocation)
    require(ranlib.is_absolute(),
            "candidate ranlib invocation is not absolute")
    with _RetainedFileSnapshot(
            archiver, "candidate archiver") as archiver_snapshot, \
            _RetainedFileSnapshot(
                ranlib, "candidate ranlib") as ranlib_snapshot, \
            _RetainedFileSnapshot(
                compiler,
                "candidate executable-link compiler") as link_compiler:
        archiver_identity = archiver_snapshot.identity
        ranlib_identity = ranlib_snapshot.identity
        require(link_compiler.identity == compiler_identity,
                "candidate C++ compiler changed between version and link "
                "recipe validation")
        archive_objects, archive_link_commands = _archive_recipe_semantics(
            archive_recipe_bytes, build, archive, archiver, ranlib)
        executable_objects, executable_link_command = \
            _executable_recipe_semantics(
                executable_recipe_bytes, build, expected_executable, archive,
                compiler)
    require(len(executable_objects) == 1,
            "candidate executable does not have one benchmark object")
    archive_object_operands = dict(zip(
        archive_objects, archive_link_commands[0][3:]))
    executable_object_operands = {
        executable_objects[0]: executable_link_command[7],
    }

    tracked = {
        (source / relative).resolve(strict=True)
        for relative in manifest_by_path
    }
    expected_library_sources = _expected_library_sources(
        source, cache, tracked,
        expected_configuration_schema=_expected_configuration_schema)
    closure_records: list[dict[str, Any]] = []
    archive_sources: set[Path] = set()
    for role, objects in (("archive", archive_objects),
                          ("benchmark", executable_objects)):
        for obj in objects:
            require(obj in by_output,
                    f"{role} object has no exact compile command: {obj}")
            (entry, tokens, normalized_arguments, source_operand,
             recipe_identity) = \
                by_output[obj]
            normalized_recipe_operand = (
                archive_object_operands[obj] if role == "archive"
                else executable_object_operands[obj])
            _require_exact_compile_object_operand(
                entry, normalized_recipe_operand, build,
                f"candidate {role} recipe object")
            require(source_operand.is_relative_to(source) and
                    source_operand in tracked,
                    f"{role} object source is not tracked at candidate HEAD")
            require(source_operand.suffix != ".c",
                    f"{role} production closure contains a C object")
            relative_source = source_operand.relative_to(source).as_posix()
            source_record = manifest_by_path.get(relative_source)
            require(source_record is not None,
                    f"{role} object source is absent from the stable tracked "
                    "source manifest")
            source_identity = {
                "path": str(logical_source / relative_source),
                "sha256": source_record["sha256"],
                "size": source_record["size"],
                "mode": source_record["mode"],
            }
            object_identity = artifact_identity(
                obj, f"candidate {role} object")
            flag_profile = _validate_compile_flags(
                tokens, source_operand,
                source_root=source, build_root=build, cache=cache,
                library_sources=expected_library_sources,
                compiler_path=compiler_identity["path"])
            if role == "archive":
                archive_sources.add(source_operand)
            else:
                require(source_operand ==
                        (source / "bench/leopard2/benchmark.cpp").resolve(
                            strict=True),
                        "candidate executable was not compiled from the exact "
                        "Leopard2 benchmark source")
            closure_records.append({
                "role": role,
                "source": source_identity,
                "object": object_identity,
                "compile_entry": {
                    "directory": entry["directory"],
                    "file": entry["file"],
                    "output": entry["output"],
                    "arguments": tokens,
                    "normalized_arguments": normalized_arguments,
                    **recipe_identity,
                },
                "flag_profile": flag_profile,
            })
    require(archive_sources == expected_library_sources,
            "candidate archive source closure differs from the production "
            "pure-AVX2 source set")

    archive_records = [record for record in closure_records
                       if record["role"] == "archive"]
    benchmark_record = next(record for record in closure_records
                            if record["role"] == "benchmark")
    object_identity_by_name = {
        Path(record["object"]["path"]).name: record["object"]
        for record in archive_records
    }
    archive_override = artifact_overrides.get(
        Path(os.path.abspath(os.fspath(archive))))
    executable_override = artifact_overrides.get(
        Path(os.path.abspath(os.fspath(expected_executable))))
    archive_context = (
        nullcontext(archive_override) if archive_override is not None else
        _RetainedFileSnapshot(archive, "candidate production archive"))
    executable_context = (
        nullcontext(executable_override)
        if executable_override is not None else
        _RetainedFileSnapshot(
            expected_executable, "candidate benchmark executable"))
    with archive_context as archive_snapshot, \
            _RetainedFileSnapshot(
                archiver,
                "candidate archive inventory tool") as inventory_archiver, \
            executable_context as executable_snapshot:
        require(archive_snapshot is not None and
                executable_snapshot is not None,
                "candidate sealed artifact context is empty")
        archive_snapshot.verify()
        executable_snapshot.verify()
        require(inventory_archiver.identity == archiver_identity,
                "candidate archiver changed before archive inspection")
        archive_identity = archive_snapshot.identity
        executable_identity = executable_snapshot.identity
        require(stat.S_ISREG(executable_identity["mode"]) and
                executable_identity["mode"] & 0o111,
                "candidate benchmark is not executable")
        require(all(
                    archive_identity["mtime_ns"] >=
                    record["object"]["mtime_ns"]
                    for record in archive_records),
                "candidate archive predates one of its exact objects")
        require(
            executable_identity["mtime_ns"] >=
            archive_identity["mtime_ns"] and
            executable_identity["mtime_ns"] >=
            benchmark_record["object"]["mtime_ns"],
            "candidate benchmark predates one of its exact link inputs")

        inherited_for_archive = tuple(inherited_descriptors) + (
            archive_snapshot.descriptor,)
        members = _run(
            (str(archiver), "t", archive_snapshot.proc_path),
            "candidate archive inventory",
            inherited_descriptors=inherited_for_archive,
            executable_descriptor=
            inventory_archiver.executable_descriptor).decode(
                "utf-8", errors="strict").splitlines()
        expected_members = [path.name for path in archive_objects]
        require(
            members == expected_members and
            len(members) == len(set(members)),
            "candidate archive members differ from its exact object recipe")
        require(set(object_identity_by_name) == set(members),
                "candidate archive object basenames are ambiguous")
        archive_member_identities = []
        for member in members:
            content = _run(
                (str(archiver), "p", archive_snapshot.proc_path, member),
                f"candidate archive member {member}",
                maximum_bytes=MAX_FILE_BYTES,
                inherited_descriptors=inherited_for_archive,
                executable_descriptor=
                inventory_archiver.executable_descriptor)
            identity = {
                "member": member, "size": len(content),
                "sha256": hashlib.sha256(content).hexdigest(),
            }
            expected_object = object_identity_by_name[member]
            require(
                identity["size"] == expected_object["size"] and
                identity["sha256"] == expected_object["sha256"],
                f"candidate archive member bytes differ from object: "
                f"{member}")
            archive_member_identities.append(identity)

    return {
        "schema": PRODUCTION_BUILD_CLOSURE_SCHEMA,
        "build_root": str(build),
        "physical_source_root": str(source),
        "source_root": str(logical_source),
        "tracked_source_manifest": source_manifest,
        "executable_target": executable_target,
        "validated_cache": {
            **required_exact,
            "LEO2_BUILD_ALLK_DIAGNOSTIC": cache.get(
                "LEO2_BUILD_ALLK_DIAGNOSTIC"),
            "LEO2_BUILD_TESTS": cache.get("LEO2_BUILD_TESTS"),
            "ENABLE_OPENMP": cache.get("ENABLE_OPENMP"),
            "LEO2_BACKEND_VARIANT": variant,
            "LEO2_BENCHMARK_GIT_EXECUTABLE": benchmark_git_value,
            "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN":
                cache.get("LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN"),
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR":
                cache.get("LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR"),
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED":
                cache.get(
                    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED"),
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE":
                cache.get("LEO2_EXPERIMENT_CAUCHY_LOG_REUSE"),
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT":
                cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"),
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE":
                cache.get("LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE"),
            "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING":
                cache.get("LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING"),
            "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING":
                cache.get("LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING"),
            "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING":
                cache.get("LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING"),
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED":
                cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED"),
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT":
                cache.get("LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT"),
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": small_direct_mode,
            "LEOPARD_ENABLE_GF16": cache.get("LEOPARD_ENABLE_GF16"),
            "LEO2_FLAG_FALIGN_FUNCTIONS_64":
                cache.get("LEO2_FLAG_FALIGN_FUNCTIONS_64"),
            "LEO2_FLAG_MAVX2": cache.get("LEO2_FLAG_MAVX2"),
            "LEO2_FLAG_MNO_AVX512F": cache.get("LEO2_FLAG_MNO_AVX512F"),
            "LEO2_FLAG_MAVX512F": cache.get("LEO2_FLAG_MAVX512F"),
            "LEO2_FLAG_MAVX512BW": cache.get("LEO2_FLAG_MAVX512BW"),
            "LEO2_FLAG_MAVX512VL": cache.get("LEO2_FLAG_MAVX512VL"),
            "LEO2_LOCATOR_GIT_EXECUTABLE": locator_git_value,
            "CMAKE_CXX_FLAGS": cache.get("CMAKE_CXX_FLAGS"),
            "CMAKE_CXX_FLAGS_RELEASE": cache.get("CMAKE_CXX_FLAGS_RELEASE"),
            "CMAKE_AR": cache.get("CMAKE_AR"),
            "CMAKE_C_COMPILER": cache.get("CMAKE_C_COMPILER"),
            "CMAKE_CXX_COMPILER": cache.get("CMAKE_CXX_COMPILER"),
            "CMAKE_LINKER": cache.get("CMAKE_LINKER"),
            "CMAKE_MAKE_PROGRAM": cache.get("CMAKE_MAKE_PROGRAM"),
            "CMAKE_RANLIB": cache.get("CMAKE_RANLIB"),
        },
        "cmake_cache": cache_identity,
        "compile_commands": commands_identity,
        "archive_link_recipe": archive_recipe_identity,
        "executable_link_recipe": executable_recipe_identity,
        "archive_link_commands": archive_link_commands,
        "executable_link_command": executable_link_command,
        "compiler": compiler_identity,
        "compiler_version_sha256": hashlib.sha256(compiler_version).hexdigest(),
        "c_compiler": c_compiler_identity,
        "make_program": make_program_identity,
        "linker": linker_identity,
        "compiler_subtools": compiler_subtools,
        "benchmark_git": benchmark_git_identity,
        "archiver": archiver_identity,
        "ranlib": ranlib_identity,
        "archive_members": members,
        "archive_member_identities": archive_member_identities,
        "archive": archive_identity,
        "executable": executable_identity,
        "source_object_compile_closure": sorted(
            closure_records, key=lambda record: (
                record["role"], record["source"]["path"],
                record["object"]["path"])),
    }


def compare_reproducible_builds(
    candidate: Mapping[str, Any], rebuilt: Mapping[str, Any],
) -> dict[str, Any]:
    """Require a clean rebuild to reproduce every linked project byte."""
    candidate_schema = candidate.get("schema")
    require(
            candidate_schema in {
                PRODUCTION_BUILD_CLOSURE_SCHEMA_V1,
                PRODUCTION_BUILD_CLOSURE_SCHEMA,
            } and
            rebuilt.get("schema") == candidate_schema,
            "candidate/rebuild provenance schema differs")
    require(candidate.get("source_root") == rebuilt.get("source_root") and
            candidate.get("executable_target") ==
            rebuilt.get("executable_target"),
            "candidate/rebuild source or target differs")
    require(candidate.get("tracked_source_manifest") ==
            rebuilt.get("tracked_source_manifest"),
            "clean rebuild tracked source tree differs")
    require(candidate.get("validated_cache") ==
            rebuilt.get("validated_cache"),
            "clean rebuild effective CMake semantics differ")
    source = Path(str(candidate["source_root"]))

    def closure_map(record: Mapping[str, Any]) -> dict[tuple[str, str], Any]:
        result = {}
        closure = record.get("source_object_compile_closure")
        require(isinstance(closure, list) and closure,
                "candidate/rebuild compile closure is empty")
        for item in closure:
            require(isinstance(item, dict) and
                    isinstance(item.get("role"), str) and
                    isinstance(item.get("source"), dict) and
                    isinstance(item.get("object"), dict),
                    "candidate/rebuild compile closure is malformed")
            source_path = Path(str(item["source"].get("path", "")))
            require(source_path.is_relative_to(source),
                    "candidate/rebuild source escapes the bound source root")
            key = (item["role"], source_path.relative_to(source).as_posix())
            require(key not in result,
                    "candidate/rebuild compile closure contains duplicates")
            result[key] = item
        return result

    candidate_closure = closure_map(candidate)
    rebuilt_closure = closure_map(rebuilt)
    require(set(candidate_closure) == set(rebuilt_closure),
            "clean rebuild source/object closure differs")
    object_comparisons = []
    for key in sorted(candidate_closure):
        original = candidate_closure[key]
        reproduction = rebuilt_closure[key]
        require(original["source"]["sha256"] ==
                reproduction["source"]["sha256"] and
                original["source"]["size"] ==
                reproduction["source"]["size"],
                f"clean rebuild source bytes differ: {key[1]}")
        require(original["flag_profile"] == reproduction["flag_profile"],
                f"clean rebuild ISA flag profile differs: {key[1]}")
        original_entry = original.get("compile_entry")
        rebuilt_entry = reproduction.get("compile_entry")
        require(isinstance(original_entry, dict) and
                isinstance(rebuilt_entry, dict) and
                isinstance(original_entry.get("arguments"), list) and
                original_entry.get("arguments") and
                all(isinstance(argument, str) and argument
                    for argument in original_entry["arguments"]) and
                isinstance(rebuilt_entry.get("arguments"), list) and
                rebuilt_entry.get("arguments") and
                all(isinstance(argument, str) and argument
                    for argument in rebuilt_entry["arguments"]) and
                isinstance(original_entry.get("normalized_arguments"), list)
                and original_entry.get("normalized_arguments") and
                all(isinstance(argument, str) and argument
                    for argument in
                    original_entry["normalized_arguments"]) and
                isinstance(rebuilt_entry.get("normalized_arguments"), list)
                and rebuilt_entry.get("normalized_arguments") and
                all(isinstance(argument, str) and argument
                    for argument in
                    rebuilt_entry["normalized_arguments"]),
                f"clean rebuild compile argv is malformed: {key[1]}")
        candidate_build_root = Path(
            str(candidate.get("build_root", "")))
        rebuilt_build_root = Path(str(rebuilt.get("build_root", "")))
        candidate_physical_source = Path(
            str(candidate.get("physical_source_root", "")))
        rebuilt_physical_source = Path(
            str(rebuilt.get("physical_source_root", "")))
        require(candidate_build_root.is_absolute() and
                rebuilt_build_root.is_absolute() and
                candidate_physical_source.is_absolute() and
                rebuilt_physical_source.is_absolute() and
                _normalize_build_argv(
                    original_entry["arguments"], candidate_build_root,
                    candidate_physical_source) ==
                original_entry["normalized_arguments"] and
                _normalize_build_argv(
                    rebuilt_entry["arguments"], rebuilt_build_root,
                    rebuilt_physical_source) ==
                rebuilt_entry["normalized_arguments"],
                f"clean rebuild compile argv normalization is invalid: "
                f"{key[1]}")
        require(original_entry["normalized_arguments"] ==
                rebuilt_entry["normalized_arguments"],
                f"clean rebuild compile argv differs: {key[1]}")
        original_representation = original_entry.get("representation")
        rebuilt_representation = rebuilt_entry.get("representation")
        require(original_representation in {"arguments", "command"} and
                rebuilt_representation in {"arguments", "command"} and
                original_representation == rebuilt_representation,
                f"clean rebuild compile recipe representation differs: "
                f"{key[1]}")
        if original_representation == "command":
            original_command = original_entry.get("command")
            rebuilt_command = rebuilt_entry.get("command")
            original_normalized_command = original_entry.get(
                "normalized_command")
            rebuilt_normalized_command = rebuilt_entry.get(
                "normalized_command")
            require(isinstance(original_command, str) and original_command and
                    isinstance(rebuilt_command, str) and rebuilt_command and
                    isinstance(original_normalized_command, str) and
                    original_normalized_command and
                    isinstance(rebuilt_normalized_command, str) and
                    rebuilt_normalized_command and
                    _compile_tokens({"command": original_command}) ==
                    original_entry["arguments"] and
                    _compile_tokens({"command": rebuilt_command}) ==
                    rebuilt_entry["arguments"] and
                    shlex.join(original_entry["normalized_arguments"]) ==
                    original_normalized_command and
                    shlex.join(rebuilt_entry["normalized_arguments"]) ==
                    rebuilt_normalized_command,
                    f"clean rebuild compile recipe normalization is invalid: "
                    f"{key[1]}")
            require(original_normalized_command ==
                    rebuilt_normalized_command,
                    f"clean rebuild raw compile recipe differs: {key[1]}")
        else:
            require("command" not in original_entry and
                    "normalized_command" not in original_entry and
                    "command" not in rebuilt_entry and
                    "normalized_command" not in rebuilt_entry,
                    f"clean rebuild arguments-form compile recipe contains "
                    f"string-form metadata: {key[1]}")
        require(original["object"]["sha256"] ==
                reproduction["object"]["sha256"] and
                original["object"]["size"] ==
                reproduction["object"]["size"],
                f"clean rebuild object bytes differ: {key[1]}")
        object_comparisons.append({
            "role": key[0], "source": key[1],
            "sha256": original["object"]["sha256"],
            "size": original["object"]["size"],
            "flag_profile": original["flag_profile"],
        })

    for recipe in ("archive_link_recipe", "executable_link_recipe"):
        original = candidate.get(recipe)
        reproduction = rebuilt.get(recipe)
        require(isinstance(original, dict) and
                isinstance(reproduction, dict) and
                original.get("sha256") == reproduction.get("sha256") and
                original.get("size") == reproduction.get("size"),
                f"clean rebuild {recipe.replace('_', ' ')} identity differs")
    candidate_archive_commands = candidate.get("archive_link_commands")
    rebuilt_archive_commands = rebuilt.get("archive_link_commands")
    require(isinstance(candidate_archive_commands, list) and
            len(candidate_archive_commands) == 2 and
            all(isinstance(command, list) and command and
                all(isinstance(argument, str) and argument
                    for argument in command)
                for command in candidate_archive_commands) and
            isinstance(rebuilt_archive_commands, list) and
            len(rebuilt_archive_commands) == 2 and
            all(isinstance(command, list) and command and
                all(isinstance(argument, str) and argument
                    for argument in command)
                for command in rebuilt_archive_commands),
            "clean rebuild archive link recipe semantics are malformed")
    require(candidate_archive_commands == rebuilt_archive_commands,
            "clean rebuild archive link recipe semantics differ")
    candidate_executable_command = candidate.get(
        "executable_link_command")
    rebuilt_executable_command = rebuilt.get("executable_link_command")
    require(isinstance(candidate_executable_command, list) and
            candidate_executable_command and
            all(isinstance(argument, str) and argument
                for argument in candidate_executable_command) and
            isinstance(rebuilt_executable_command, list) and
            rebuilt_executable_command and
            all(isinstance(argument, str) and argument
                for argument in rebuilt_executable_command),
            "clean rebuild executable link recipe semantics are malformed")
    require(candidate_executable_command == rebuilt_executable_command,
            "clean rebuild executable link recipe semantics differ")

    require(candidate.get("archive_member_identities") ==
            rebuilt.get("archive_member_identities"),
            "clean rebuild archive member bytes differ")
    for artifact in ("archive", "executable"):
        original = candidate.get(artifact)
        reproduction = rebuilt.get(artifact)
        require(isinstance(original, dict) and
                isinstance(reproduction, dict) and
                original.get("sha256") == reproduction.get("sha256") and
                original.get("size") == reproduction.get("size"),
                f"clean rebuild {artifact} bytes differ")
    require(candidate.get("compiler", {}).get("sha256") ==
            rebuilt.get("compiler", {}).get("sha256") and
            candidate.get("compiler_version_sha256") ==
            rebuilt.get("compiler_version_sha256"),
            "clean rebuild compiler identity differs")
    for tool in (
            "c_compiler", "archiver", "ranlib", "make_program", "linker",
            "benchmark_git"):
        original = candidate.get(tool)
        reproduction = rebuilt.get(tool)
        require(isinstance(original, dict) and
                isinstance(reproduction, dict) and
                original.get("sha256") == reproduction.get("sha256") and
                original.get("size") == reproduction.get("size"),
                f"clean rebuild {tool.replace('_', ' ')} identity differs")
    require(candidate.get("compiler_subtools") ==
            rebuilt.get("compiler_subtools"),
            "clean rebuild compiler subtool identities differ")
    return {
        "schema": (
            REPRODUCIBLE_BUILD_PROOF_SCHEMA
            if candidate_schema == PRODUCTION_BUILD_CLOSURE_SCHEMA else
            REPRODUCIBLE_BUILD_PROOF_SCHEMA_V2),
        "method": "runner-owned-empty-directory-configure-build-byte-compare",
        "source_root": str(source),
        "executable_target": candidate["executable_target"],
        "compiler_sha256": candidate["compiler"]["sha256"],
        "objects": object_comparisons,
        "archive_members": candidate["archive_member_identities"],
        "archive_sha256": candidate["archive"]["sha256"],
        "executable_sha256": candidate["executable"]["sha256"],
    }


def _reproducible_configure_argv(
    source: Path, build: Path, cache: Mapping[str, Any], *,
    cmake_path: Path | str | None = None,
) -> list[str]:
    """Recreate the candidate's exact retained CMake cache semantics.

    An untyped ``-Dname=value`` initially creates an ``UNINITIALIZED`` cache
    entry.  CMake does not necessarily retag an already-defined compiler or
    compiler-probe variable later: in particular, CMAKE_CXX_COMPILER remains
    STRING and preseeded check_cxx_compiler_flag results remain UNINITIALIZED.
    The production validator deliberately requires the canonical FILEPATH and
    INTERNAL types, so the replay must preserve every declared type explicitly.
    """
    if cmake_path is None:
        cmake = Path("/usr/bin/cmake").resolve(strict=True)
    else:
        cmake = Path(cmake_path)
        require(cmake.is_absolute() and "\0" not in str(cmake),
                "retained clean-rebuild CMake path is invalid")
        _require_safe_unicode(
            str(cmake), "retained clean-rebuild CMake path")
    c_compiler = cache.get("CMAKE_C_COMPILER")
    require(isinstance(c_compiler, str) and c_compiler,
            "candidate build lacks a C compiler needed for clean rebuild")
    compiler = cache.get("CMAKE_CXX_COMPILER")
    require(isinstance(compiler, str) and compiler,
            "candidate build lacks a C++ compiler needed for clean rebuild")
    # CMake deliberately retags an explicitly -D-seeded compiler as STRING.
    # The caller supplies these drivers through a bounded CC/CXX environment
    # while executing this direct CMake argv from a sealed descriptor.
    configure = [
        str(cmake),
        "-S", str(source), "-B", str(build),
        "-G", "Unix Makefiles",
    ]
    cmake_values = {
        "CMAKE_AR": cache.get("CMAKE_AR"),
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_CXX_FLAGS": cache.get("CMAKE_CXX_FLAGS"),
        "CMAKE_CXX_FLAGS_RELEASE": cache.get("CMAKE_CXX_FLAGS_RELEASE"),
        "CMAKE_LINKER": cache.get("CMAKE_LINKER"),
        "CMAKE_MAKE_PROGRAM": cache.get("CMAKE_MAKE_PROGRAM"),
        "ENABLE_OPENMP": cache.get("ENABLE_OPENMP"),
        "LEO2_BACKEND_VARIANT": cache.get("LEO2_BACKEND_VARIANT"),
        "LEO2_BENCHMARK_GIT_EXECUTABLE":
            cache.get("LEO2_BENCHMARK_GIT_EXECUTABLE"),
        "LEO2_BUILD_ALLK_DIAGNOSTIC": cache.get(
            "LEO2_BUILD_ALLK_DIAGNOSTIC"),
        "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_BUILD_TESTS": cache.get("LEO2_BUILD_TESTS"),
        "LEO2_ENABLE_CUDA": "OFF",
        "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN":
            cache.get("LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN"),
        "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE":
            cache.get("LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE"),
        "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE":
            cache.get("LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"),
        "LEOPARD_ENABLE_GF8": "ON",
        "LEOPARD_ENABLE_GF16": cache.get("LEOPARD_ENABLE_GF16"),
        "LEO2_FLAG_FALIGN_FUNCTIONS_64":
            cache.get("LEO2_FLAG_FALIGN_FUNCTIONS_64"),
        "LEO2_FLAG_MAVX2": cache.get("LEO2_FLAG_MAVX2"),
        "LEO2_FLAG_MNO_AVX512F": cache.get("LEO2_FLAG_MNO_AVX512F"),
        "LEO2_FLAG_MAVX512F": "FALSE",
        "LEO2_FLAG_MAVX512BW": "FALSE",
        "LEO2_FLAG_MAVX512VL": "FALSE",
        "LEO2_LOCATOR_GIT_EXECUTABLE":
            cache.get("LEO2_LOCATOR_GIT_EXECUTABLE"),
        "CMAKE_RANLIB": cache.get("CMAKE_RANLIB"),
    }
    configuration_schema = cache.get(
        "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA")
    if configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V3,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        cmake_values["LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"] = \
            cache.get("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT")
    if configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V4,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        cmake_values.update({
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR":
                cache.get("LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR"),
            "LEO2_EXPERIMENT_CAUCHY_LOG_REUSE":
                cache.get("LEO2_EXPERIMENT_CAUCHY_LOG_REUSE"),
            "LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING":
                cache.get("LEO2_EXPERIMENT_HIGH_T8_PARTIAL_BINDING"),
            "LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING":
                cache.get("LEO2_EXPERIMENT_HIGH_T8_RAGGED_BINDING"),
            "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING":
                cache.get("LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING"),
            "LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT":
                cache.get("LEO2_EXPERIMENT_ONE_SHOT_EQUAL_ROUNDED_DIRECT"),
        })
    if configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        cmake_values.update({
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED":
                cache.get(
                    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED"),
            "LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED":
                cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_GENERATED"),
        })
    if configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        cmake_values.update({
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK":
                cache.get(
                    "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_TWO_BLOCK"),
            "LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED":
                cache.get("LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED"),
            "LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK":
                cache.get("LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK"),
            "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL":
                cache.get("LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL"),
        })
    if configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        cmake_values["LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT"] = cache.get(
            "LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT")
    if configuration_schema in {
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
            BENCHMARK_BUILD_CONFIGURATION_SCHEMA}:
        cmake_values.update({
            "LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS":
                cache.get("LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS"),
            "LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK":
                cache.get("LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK"),
        })
    if configuration_schema == BENCHMARK_BUILD_CONFIGURATION_SCHEMA:
        cmake_values["LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED"] = cache.get(
            "LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED")
    require(all(value is not None for value in cmake_values.values()),
            "candidate build lacks a CMake value needed for clean rebuild")
    for name, value in cmake_values.items():
        allowed_types = CMAKE_CACHE_REQUIRED_ENTRY_TYPES.get(name)
        require(allowed_types is not None and len(allowed_types) == 1,
                f"clean rebuild has no unique cache type for {name}")
        entry_type = next(iter(allowed_types))
        configure.append(f"-D{name}:{entry_type}={value}")
    return configure


def _replay_git_identity_bytes(
    source: Path, manifest: Mapping[str, Any], *,
    interpreter_descriptor: int | None = None,
) -> tuple[bytes, dict[str, Any]]:
    """Generate the exact retained Git emulator and its semantic identity."""
    git = manifest.get("git")
    require(isinstance(git, dict) and
            isinstance(git.get("commit"), str) and
            isinstance(git.get("tree"), str) and
            type(git.get("dirty")) is bool,
            "clean rebuild tracked Git identity is malformed")
    commit = git["commit"]
    tree = git["tree"]
    require(re.fullmatch(r"[0-9a-f]{40,64}", commit) is not None and
            re.fullmatch(r"[0-9a-f]{40,64}", tree) is not None,
            "clean rebuild tracked Git identity is non-canonical")
    # The generated benchmark header consumes only commit, tree, and the
    # status-empty predicate.  This direct argv emulator prevents CMake from
    # consulting the mutable original .git directory during replay.
    if interpreter_descriptor is None:
        interpreter = "/bin/sh"
    else:
        require(type(interpreter_descriptor) is int and
                interpreter_descriptor >= 0,
                "clean rebuild Git interpreter descriptor is invalid")
        interpreter = f"/proc/self/fd/{interpreter_descriptor}"
    script = (
        f"#!{interpreter}\n"
        "set -eu\n"
        "if [ \"${1-}\" = \"-C\" ]; then shift 2; fi\n"
        "if [ \"${1-}\" = \"rev-parse\" ]; then\n"
        "  shift\n"
        "  if [ \"${1-}\" = \"--show-toplevel\" ]; then\n"
        f"    printf '%s\\n' {shlex.quote(str(source))}; exit 0\n"
        "  fi\n"
        "  if [ \"${1-}\" = \"--verify\" ]; then shift; fi\n"
        "  case \"${1-}\" in\n"
        f"    HEAD) printf '%s\\n' {shlex.quote(commit)} ;;\n"
        f"    'HEAD^{{tree}}') printf '%s\\n' {shlex.quote(tree)} ;;\n"
        "    *) exit 2 ;;\n"
        "  esac\n"
        "  exit 0\n"
        "fi\n"
        "if [ \"${1-}\" = \"status\" ]; then\n" +
        ("  printf '%s\\n' ' M .leo2-retained-dirty'\n"
         if git["dirty"] else ": clean status\n") +
        "  exit 0\n"
        "fi\n"
        "exit 2\n"
    )
    encoded = script.encode("utf-8", errors="strict")
    return encoded, {
        "source_root": str(source),
        "commit": commit,
        "tree": tree,
        "dirty": git["dirty"],
        "interpreter_descriptor": interpreter_descriptor,
    }


def _write_replay_git_identity(
    path: Path, source: Path, manifest: Mapping[str, Any],
    *, interpreter_descriptor: int | None = None,
) -> tuple[bytes, dict[str, Any]]:
    if interpreter_descriptor is not None:
        os.fstat(interpreter_descriptor)
    encoded, specification = _replay_git_identity_bytes(
        source, manifest,
        interpreter_descriptor=interpreter_descriptor)
    require(not path.exists(), "clean rebuild Git identity path already exists")
    _write_private_executable(path, encoded)
    return encoded, specification


def _retarget_replay_git_transport(
    build: Path, identity_git: Path, transport_git: Path,
) -> None:
    """Change only the generated refresh command, not its attested identity.

    The benchmark configuration digest intentionally includes the configured
    Git pathname.  Reconfiguring with the private emulator would therefore
    change benchmark.cpp.  Configure with the candidate's retained pathname,
    then retarget the one runner-owned Make recipe to the private transport.
    """
    build = build.resolve(strict=True)
    require(identity_git.is_absolute() and transport_git.is_absolute() and
            identity_git.is_file() and transport_git.is_file(),
            "clean replay Git identity/transport path is invalid")
    for value in (str(identity_git), str(transport_git)):
        _require_safe_unicode(value, "clean replay Git transport path")
        require(re.fullmatch(r"[/A-Za-z0-9_.+-]+", value) is not None,
                "clean replay Git transport path requires Make escaping")
    recipe = (
        build /
        "CMakeFiles/leopard2_benchmark_source_attestation_refresh.dir/"
        "build.make")
    require(recipe.resolve(strict=True).is_relative_to(build) and
            not recipe.is_symlink(),
            "clean replay benchmark refresh recipe is unsafe")
    identity, content = file_snapshot(
        recipe, "clean replay benchmark refresh recipe",
        maximum_bytes=MAX_METADATA_BYTES)
    needle = (
        "-DLEO2_GIT_EXECUTABLE=" + str(identity_git)).encode("utf-8")
    replacement = (
        "-DLEO2_GIT_EXECUTABLE=" + str(transport_git)).encode("utf-8")
    require(content.count(needle) == 1 and replacement not in content,
            "clean replay benchmark refresh recipe does not contain one "
            "exact Git transport")
    rewritten = content.replace(needle, replacement)
    temporary = recipe.with_name(recipe.name + ".leo2-git-transport.tmp")
    descriptor = -1
    try:
        descriptor = os.open(
            temporary,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL |
            getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0),
            stat.S_IMODE(identity["mode"]))
        view = memoryview(rewritten)
        while view:
            written = os.write(descriptor, view)
            require(written > 0,
                    "clean replay Git transport rewrite stalled")
            view = view[written:]
        os.fsync(descriptor)
        consumed = [False]
        try:
            _close_raw_descriptor_with_precedence(
                descriptor, None,
                "clean replay Git transport temporary descriptor",
                "cannot close clean replay Git transport temporary "
                "descriptor",
                ownership_consumed=consumed)
        finally:
            if consumed[0]:
                descriptor = -1
        os.replace(temporary, recipe)
        directory_descriptor = os.open(
            recipe.parent,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_DIRECTORY", 0) |
            getattr(os, "O_NOFOLLOW", 0))
        try:
            os.fsync(directory_descriptor)
        finally:
            _close_raw_descriptor_with_precedence(
                directory_descriptor, sys.exception(),
                "clean replay Git transport directory descriptor",
                "cannot close clean replay Git transport directory "
                "descriptor")
    except BaseException as primary:
        cleanup_failure: BaseException | None = None
        if descriptor >= 0:
            consumed = [False]
            try:
                _close_raw_descriptor_with_precedence(
                    descriptor, None,
                    "interrupted clean replay Git transport descriptor",
                    "cannot close interrupted clean replay Git transport "
                    "descriptor",
                    ownership_consumed=consumed)
            except BaseException as error:
                cleanup_failure = error
            finally:
                if consumed[0]:
                    descriptor = -1
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        except BaseException as error:
            cleanup_failure = _owner_exception_precedence(
                cleanup_failure, error,
                "clean replay Git transport cleanup")
        if cleanup_failure is not None:
            _raise_owner_exit_failure(
                primary, cleanup_failure,
                "clean replay Git transport owner",
                "clean replay Git transport cleanup failed: "
                f"{cleanup_failure}; primary failure: "
                f"{type(primary).__name__}: {primary}")
        raise
    rewritten_identity, observed = file_snapshot(
        recipe, "retargeted clean replay benchmark refresh recipe",
        maximum_bytes=MAX_METADATA_BYTES)
    require(observed == rewritten and
            stat.S_IMODE(rewritten_identity["mode"]) ==
            stat.S_IMODE(identity["mode"]),
            "clean replay Git transport rewrite was not retained exactly")


def _write_private_executable(path: Path, content: bytes) -> None:
    """Create one private file without following its parent or final name.

    Replay paths below the configured build tree are attacker-controlled
    outputs until they are independently reduced to candidate semantics.
    Opening ``path`` directly with ``O_NOFOLLOW`` protects only its final
    component: a generated parent symlink could still redirect the write.
    Retain and verify the literal parent directory first, then create the
    basename relative to that descriptor.
    """
    require(content, "private replay executable content is empty")
    requested = Path(os.path.abspath(os.fspath(path)))
    parent = requested.parent
    name = requested.name
    require(requested.is_absolute() and name not in ("", ".", "..") and
            "/" not in name and "\0" not in name,
            "private replay executable destination is invalid")
    directory = _OwnedDescriptor()
    output = _OwnedDescriptor()
    created = False
    primary: BaseException | None = None
    cleanup_failure: BaseException | None = None
    cleanup_details: list[str] = []

    def retain_cleanup_failure(error: BaseException, action: str) -> None:
        nonlocal cleanup_failure
        cleanup_details.append(
            f"{action}: {type(error).__name__}: {error}")
        observed = _ordinary_owner_cleanup_error(
            error, f"{action}: {error}")
        cleanup_failure = _owner_exception_precedence(
            cleanup_failure, observed,
            "private replay executable cleanup")

    try:
        try:
            resolved_parent = parent.resolve(strict=True)
        except (OSError, RuntimeError) as error:
            raise BuildProvenanceError(
                f"private replay executable parent is invalid: {error}") \
                from error
        require(resolved_parent == parent,
                "private replay executable parent follows a symlink")
        directory_descriptor = directory.open(
            parent,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_DIRECTORY", 0) |
            getattr(os, "O_NOFOLLOW", 0))
        opened_parent = os.fstat(directory_descriptor)
        lexical_parent = parent.lstat()
        require(stat.S_ISDIR(opened_parent.st_mode) and
                stat.S_ISDIR(lexical_parent.st_mode) and
                (opened_parent.st_dev, opened_parent.st_ino) ==
                (lexical_parent.st_dev, lexical_parent.st_ino) and
                parent.resolve(strict=True) == parent,
                "private replay executable parent changed while opened")
        descriptor = output.open(
            name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL |
            getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
            0o700, dir_fd=directory_descriptor)
        created = True
        os.fchmod(descriptor, 0o700)
        view = memoryview(content)
        while view:
            written = os.write(descriptor, view)
            require(written > 0,
                    "private replay executable write stalled")
            view = view[written:]
        os.fsync(descriptor)
        os.fsync(directory_descriptor)
    except BaseException as error:
        primary = error
        if created and directory.descriptor >= 0:
            try:
                os.unlink(name, dir_fd=directory.descriptor)
                os.fsync(directory.descriptor)
            except FileNotFoundError:
                pass
            except BaseException as cleanup_error:
                retain_cleanup_failure(
                    cleanup_error,
                    "cannot remove interrupted private replay executable")

    for owner, description in (
            (output, "cannot close private replay executable"),
            (directory,
             "cannot close private replay executable parent directory")):
        try:
            owner.close()
        except BaseException as cleanup_error:
            retain_cleanup_failure(cleanup_error, description)

    if cleanup_failure is not None:
        details = "; ".join(cleanup_details)
        if primary is not None:
            _raise_owner_exit_failure(
                primary, cleanup_failure,
                "private replay executable owner",
                f"private replay executable cleanup failed: {details}; "
                f"primary failure: {type(primary).__name__}: {primary}")
        if isinstance(cleanup_failure, _OWNER_TERMINAL_EXCEPTIONS):
            raise cleanup_failure.with_traceback(
                cleanup_failure.__traceback__)
        raise BuildProvenanceError(
            "private replay executable cleanup failed: " + details) from (
                cleanup_failure)
    if primary is not None:
        raise primary.with_traceback(primary.__traceback__)


def _replay_exec_wrapper_bytes(
    interpreter_descriptor: int, lexical_argv0: str,
    executable_descriptor: int, extra_arguments: Sequence[str] = (),
) -> bytes:
    """Return one exact Bash wrapper for a sealed executable descriptor."""
    require(type(interpreter_descriptor) is int and
            interpreter_descriptor >= 0 and
            type(executable_descriptor) is int and
            executable_descriptor >= 0,
            "private replay wrapper descriptors are invalid")
    require(isinstance(lexical_argv0, str) and lexical_argv0 and
            "\0" not in lexical_argv0,
            "private replay wrapper argv[0] is invalid")
    _require_safe_unicode(
        lexical_argv0, "private replay wrapper argv[0]")
    require(isinstance(extra_arguments, Sequence) and
            not isinstance(extra_arguments, (str, bytes)) and
            all(isinstance(argument, str) and argument and
                "\0" not in argument for argument in extra_arguments),
            "private replay wrapper extra arguments are invalid")
    for argument in extra_arguments:
        _require_safe_unicode(
            argument, "private replay wrapper extra argument")
    command = [
        "exec", "-a", lexical_argv0,
        f"/proc/self/fd/{executable_descriptor}",
        *extra_arguments,
    ]
    return (
        f"#!/proc/self/fd/{interpreter_descriptor}\n"
        "set -eu\n" +
        " ".join(shlex.quote(argument) for argument in command) +
        ' "$@"\n'
    ).encode("utf-8", errors="strict")


def _write_replay_exec_wrapper(
    path: Path, interpreter_descriptor: int, lexical_argv0: str,
    executable_descriptor: int, extra_arguments: Sequence[str] = (),
) -> bytes:
    """Write and return a Bash wrapper around one sealed executable."""
    os.fstat(interpreter_descriptor)
    os.fstat(executable_descriptor)
    script = _replay_exec_wrapper_bytes(
        interpreter_descriptor, lexical_argv0, executable_descriptor,
        extra_arguments)
    _write_private_executable(path, script)
    return script


def _private_replay_inventory_bytes(
    interpreter_descriptor: int, expected_files: Sequence[str],
) -> bytes:
    """Render a sealed Bash-only exact inventory checker for private tmpfs."""
    require(
        type(interpreter_descriptor) is int and interpreter_descriptor >= 0 and
        isinstance(expected_files, Sequence) and
        not isinstance(expected_files, (str, bytes)) and expected_files and
        len(expected_files) == len(set(expected_files)) and
        all(
            isinstance(relative, str) and relative and
            not Path(relative).is_absolute() and
            re.fullmatch(r"[A-Za-z0-9_.+/-]+", relative) is not None and
            all(component not in ("", ".", "..")
                for component in Path(relative).parts)
            for relative in expected_files),
        "private replay expected output inventory is malformed")
    assignments = "\n".join(
        f"expected[{shlex.quote(relative)}]=1"
        for relative in sorted(expected_files))
    return (
        f"#!/proc/self/fd/{interpreter_descriptor}\n"
        "set -euo pipefail\n"
        "shopt -s globstar nullglob dotglob\n"
        "declare -A expected=()\n"
        f"{assignments}\n"
        "observed=0\n"
        "for path in **; do\n"
        "  if [[ -d \"$path\" && ! -L \"$path\" ]]; then continue; fi\n"
        "  [[ -f \"$path\" && ! -L \"$path\" && "
        "${expected[$path]+present} == present ]] || exit 91\n"
        "  ((observed += 1))\n"
        "done\n"
        f"(( observed == {len(expected_files)} )) || exit 92\n"
        "for path in \"${!expected[@]}\"; do\n"
        "  [[ -f \"$path\" && ! -L \"$path\" ]] || exit 93\n"
        "done\n"
    ).encode("utf-8", errors="strict")


def _retain_exact_generated_executable(
    path: Path, expected_content: bytes, label: str,
) -> _RetainedFileSnapshot:
    """Retain a generated executable only when its exact intended bytes remain."""
    require(isinstance(expected_content, bytes) and expected_content,
            f"{label} expected executable bytes are invalid")
    snapshot = _RetainedFileSnapshot(path, label)
    try:
        require(snapshot.content == expected_content,
                f"{label} changed between generation and retention")
        snapshot.executable_descriptor
        return snapshot
    except BaseException as primary:
        try:
            snapshot._close_without_verification()
        except BaseException as cleanup_error:
            _raise_owner_exit_failure(
                primary, cleanup_error,
                f"{label} exact generated executable retention",
                f"{label} exact generated executable cleanup failed: "
                f"{cleanup_error}; primary failure: "
                f"{type(primary).__name__}: {primary}")
        raise


def _replace_private_file(path: Path, content: bytes, mode: int) -> None:
    """Atomically replace one runner-owned generated recipe."""
    temporary = path.with_name(path.name + ".leo2-transport.tmp")
    owner = _OwnedDescriptor()
    try:
        descriptor = owner.open(
            temporary,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL |
            getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0),
            mode)
        # open(2) applies the process umask even when restoring a previously
        # captured mode.  Reinstate the exact retained mode before publishing.
        os.fchmod(descriptor, mode)
        view = memoryview(content)
        while view:
            written = os.write(descriptor, view)
            require(written > 0,
                    "private replay recipe write stalled")
            view = view[written:]
        os.fsync(descriptor)
        owner.close()
        os.replace(temporary, path)
    except BaseException as primary:
        cleanup_error: BaseException | None = None
        try:
            owner.close()
        except BaseException as observed:
            cleanup_error = observed
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        except BaseException as observed:
            cleanup_error = _owner_exception_precedence(
                cleanup_error, observed,
                "private replay recipe replacement cleanup")
        if cleanup_error is not None:
            _raise_owner_exit_failure(
                primary, cleanup_error,
                "private replay recipe replacement owner",
                "private replay recipe replacement cleanup failed: "
                f"{cleanup_error}; primary failure: "
                f"{type(primary).__name__}: {primary}")
        raise


_REPLAY_RECIPE_REQUIRED_PATHS = frozenset((
    "CMakeCache.txt",
    "Makefile",
    "CMakeFiles/Makefile2",
    "CMakeFiles/Makefile.cmake",
    "CMakeFiles/cmake.check_cache",
))
_REPLAY_DERIVED_MUTABLE_RECIPE_NAMES = frozenset((
    "compiler_depend.make",
    "depend.make",
))
REPLAY_EMPTY_DEPENDENCY_INCLUDE = (
    b"# Leopard2 immutable empty dependency include\n")
_REPLAY_DERIVED_DEPENDENCY_TOKEN = re.compile(
    rb"(?:compiler_)?depend\.make")
_REPLAY_SPECIAL_MAKE_VARIABLE = re.compile(
    rb"(?:^|[ \t:])"
    rb"(?:(?:override|export|unexport)[ \t]+)*"
    rb"(?:SHELL|MAKE|CMAKE_COMMAND|MAKESILENT|MAKEFILES|MAKEFLAGS|"
    rb"GNUMAKEFLAGS|MAKEOVERRIDES|"
    rb"\.SHELLFLAGS|\.RECIPEPREFIX)"
    rb"(?:[ \t]*(?:[:+?!]?=|$))")
_REPLAY_SPECIAL_MAKE_REFERENCE = re.compile(
    rb"\$(?:\((?P<paren>SHELL|MAKE|MAKESILENT|MAKEFILES|MAKEFLAGS|GNUMAKEFLAGS|"
    rb"MAKEOVERRIDES|\.SHELLFLAGS|\.RECIPEPREFIX)\)|"
    rb"\{(?P<brace>SHELL|MAKE|MAKESILENT|MAKEFILES|MAKEFLAGS|GNUMAKEFLAGS|"
    rb"MAKEOVERRIDES|\.SHELLFLAGS|\.RECIPEPREFIX)\})")


def _is_replay_recipe_relative_path(value: Any) -> bool:
    if not isinstance(value, str) or not value:
        return False
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != value:
        return False
    if path.name in _REPLAY_DERIVED_MUTABLE_RECIPE_NAMES:
        # These files are generated by the sealed cmake_depends invocation.
        # Recording either one as an input (or rewrite) would let a forged
        # proof authenticate mutable derived Make code as retained source.
        return False
    if value in _REPLAY_RECIPE_REQUIRED_PATHS:
        return True
    return (
        len(path.parts) >= 2 and path.parts[0] == "CMakeFiles" and
        (value.endswith(".make") or value.endswith(".cmake") or
         path.name == "link.txt")
    )


def _replay_recipe_candidates(build: Path) -> list[Path]:
    """Enumerate the complete code-bearing Unix-Makefiles input closure."""
    lexical = Path(os.path.abspath(os.fspath(build)))
    candidates = {
        lexical / "CMakeCache.txt",
        lexical / "Makefile",
        lexical / "CMakeFiles/Makefile2",
        lexical / "CMakeFiles/Makefile.cmake",
        lexical / "CMakeFiles/cmake.check_cache",
    }
    for pattern in (
            "CMakeFiles/**/*.make",
            "CMakeFiles/**/*.cmake",
            "CMakeFiles/**/link.txt"):
        for path in lexical.glob(pattern):
            candidates.add(path)
            require(len(candidates) <= MAX_REPLAY_RECIPE_FILES,
                    "private replay recipe inventory exceeds its file bound")
    result = []
    for path in sorted(candidates):
        try:
            status = path.lstat()
        except FileNotFoundError:
            continue
        relative = path.relative_to(lexical).as_posix()
        if path.name in _REPLAY_DERIVED_MUTABLE_RECIPE_NAMES:
            # cmake_depends creates these from retained *.cmake inputs before
            # the recursive build Make parses them.  They are outputs of the
            # sealed CMake tool, not pre-existing replay inputs.
            continue
        require(_is_replay_recipe_relative_path(relative) and
                stat.S_ISREG(status.st_mode) and
                path.resolve(strict=True).is_relative_to(lexical),
                "private replay recipe path is unsafe")
        parent = lexical
        for component in Path(relative).parts[:-1]:
            parent = parent / component
            require(not parent.is_symlink(),
                    "private replay recipe parent is a symlink")
        result.append(path)
    return result


def _replay_recipe_inventory_sha256(
    records: Sequence[Mapping[str, Any]],
) -> str:
    encoded = json.dumps(
        list(records), sort_keys=True, separators=(",", ":"),
        ensure_ascii=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _validate_replay_make_include_closure(
    path: Path, build: Path, content: bytes, candidates: set[str],
    cmake_logical_path: str,
    source_logical_path: str | None = None,
) -> None:
    """Accept only CMake's direct, same-target generated Make includes."""
    require(
        b"\0" not in content and b"\r" not in content and
        (not content or content.endswith(b"\n")),
        "private replay Make input has non-canonical line framing")
    relative = path.relative_to(build).as_posix()
    parent = Path(relative).parent.as_posix()
    allowed_names = {
        b"depend.make", b"compiler_depend.make",
        b"progress.make", b"flags.make",
    }
    observed: set[bytes] = set()
    dangerous_functions = (
        b"$(eval", b"${eval", b"$(shell", b"${shell",
        b"$(file", b"${file", b"$(call", b"${call",
        b"$(guile", b"${guile", b"$(value", b"${value",
    )
    allowed_reference = re.compile(
        rb"(?:ARGS|CMAKE_BINARY_DIR|CMAKE_COMMAND|CMAKE_SOURCE_DIR|COLOR|"
        rb"CXX_DEFINES|CXX_FLAGS|CXX_INCLUDES|C_DEFINES|C_FLAGS|C_INCLUDES|"
        rb"MAKE|MAKESILENT|VERBOSE|CMAKE_PROGRESS_[0-9]+)")
    value_variable = (
        rb"(?:ARGS|CMAKE_BINARY_DIR|CMAKE_COMMAND|CMAKE_SOURCE_DIR|COLOR|"
        rb"CXX_DEFINES|CXX_FLAGS|CXX_INCLUDES|C_DEFINES|C_FLAGS|C_INCLUDES|"
        rb"VERBOSE|CMAKE_PROGRESS_[0-9]+)")
    value_assignment_probe = re.compile(
        rb"(?:^|[ \t:])(?P<name>" + value_variable + rb")"
        rb"[ \t]*(?P<operator>:::=|::=|:=|\+=|\?=|!=|=)")
    canonical_value_assignment = re.compile(
        rb"(?P<name>" + value_variable + rb")[ \t]+=[ \t]*"
        rb"(?P<value>.*)")
    make_reference = re.compile(
        rb"\$\((?P<paren>[^()]*)\)|\$\{(?P<brace>[^{}]*)\}")
    make_alias = re.compile(
        rb"(?<![A-Za-z0-9_.+-])"
        rb"(?:[A-Za-z0-9_.+/-]*/)?g?make"
        rb"(?![A-Za-z0-9_.+-])")
    make_literal = re.compile(
        r"(?:^|[^A-Za-z0-9_.+-])"
        r"(?:[A-Za-z0-9_.+/-]*/)?g?make"
        r"(?:$|[^A-Za-z0-9_.+-])")
    safe_make_data = re.compile(
        r"-DBUILD_PROGRAM=(?:[A-Za-z0-9_.+/-]*/)?g?make")
    _require_safe_unicode(
        cmake_logical_path, "private replay CMake identity path")
    canonical_cmake_assignment = (
        b"CMAKE_COMMAND = " +
        cmake_logical_path.encode("utf-8", errors="strict"))
    canonical_binary_assignment = (
        b"CMAKE_BINARY_DIR = " +
        str(build).encode("utf-8", errors="strict"))
    canonical_source_assignment = (
        None if source_logical_path is None else
        b"CMAKE_SOURCE_DIR = " +
        source_logical_path.encode("utf-8", errors="strict"))
    if source_logical_path is not None:
        _require_safe_unicode(
            source_logical_path, "private replay source identity path")

    def validate_value_assignment(line: bytes) -> None:
        """Reject code-bearing values hidden behind allowed references."""
        probe = value_assignment_probe.search(line)
        if probe is None:
            return
        assignment = canonical_value_assignment.fullmatch(line)
        require(
            not line.startswith(b"\t") and assignment is not None and
            probe.group("operator") == b"=",
            "private replay Make input contains a non-canonical "
            "value-bearing variable assignment")
        name = assignment.group("name")
        value = assignment.group("value")
        # These variables are intentionally empty in the isolated replay
        # environment.  CMake emits references to them but no file assignment;
        # allowing an assignment would turn a harmless display/test argument
        # into arbitrary shell syntax after Make expansion.
        require(
            name not in (b"ARGS", b"COLOR", b"VERBOSE"),
            "private replay Make input assigns an external-only variable")
        if name == b"CMAKE_COMMAND":
            require(
                line == canonical_cmake_assignment,
                "private replay Make input changes an execution-control "
                "variable")
            return
        if name == b"CMAKE_BINARY_DIR":
            require(
                line == canonical_binary_assignment,
                "private replay Make input changes its canonical build "
                "directory")
            return
        require(
            not any(ord(character) in value
                    for character in SHELL_RECIPE_META),
            "private replay Make variable value contains shell code or "
            "indirection")
        try:
            decoded = value.decode("utf-8", errors="strict")
            _require_safe_unicode(
                decoded, "private replay Make variable value")
            tokens = shlex.split(decoded, posix=True)
        except (UnicodeDecodeError, ValueError) as error:
            raise BuildProvenanceError(
                "private replay Make variable value has non-canonical shell "
                f"quoting: {error}") from error
        require(
            all(not token.startswith("@") for token in tokens),
            "private replay Make variable value contains response-file "
            "indirection")
        if name == b"CMAKE_SOURCE_DIR":
            require(
                len(tokens) == 1 and tokens[0] == decoded and
                Path(decoded).is_absolute() and
                os.path.abspath(decoded) == decoded and
                (canonical_source_assignment is None or
                 line == canonical_source_assignment),
                "private replay Make source directory is non-canonical")
        elif name.startswith(b"CMAKE_PROGRESS_"):
            require(
                re.fullmatch(rb"[0-9]+", value) is not None,
                "private replay Make progress value is non-canonical")

    def shell_tokens(line: bytes) -> list[str]:
        command = line[1:]
        while command[:1] in (b"@", b"+", b"-"):
            command = command[1:]
        try:
            decoded = command.decode("utf-8", errors="strict")
            return shlex.split(decoded, posix=True)
        except (UnicodeDecodeError, ValueError) as error:
            raise BuildProvenanceError(
                "private replay Make recipe has non-canonical shell "
                f"quoting: {error}") from error

    def is_make_executable(token: str) -> bool:
        return (
            token and not token.startswith("-") and "=" not in token and
            token.rsplit("/", 1)[-1] in ("make", "gmake")
        )

    def has_make_control_option(tokens: Sequence[str]) -> bool:
        return any(
            token in ("-f", "-C", "--file", "--makefile", "--directory") or
            token.startswith(("-f", "-C", "--file=", "--makefile=",
                              "--directory="))
            for token in tokens
        )

    def invokes_noncanonical_make(line: bytes) -> bool:
        if not line.startswith(b"\t"):
            return (
                b"=" in line and
                make_alias.search(line.split(b"=", 1)[1]) is not None
            )
        tokens = shell_tokens(line)
        for token in tokens:
            if safe_make_data.fullmatch(token) is not None:
                continue
            if make_literal.search(token) is not None:
                return True
        while (tokens and
               re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", tokens[0])):
            tokens.pop(0)
        if not tokens:
            return False
        if is_make_executable(tokens[0]):
            return True
        # These command-position wrappers all execute a later argv token.
        # Searching only their tail avoids mistaking CMake data such as
        # -DBUILD_PROGRAM=/usr/bin/gmake for an invocation.
        wrapper = tokens[0].rsplit("/", 1)[-1]
        if wrapper in (
                "command", "env", "exec", "nohup", "nice", "stdbuf",
                "time"):
            return any(is_make_executable(token) for token in tokens[1:])
        variable = re.fullmatch(
            r"\$\([A-Za-z0-9_.+-]+\)", tokens[0]) is not None
        if variable and tokens[0] not in ("$(CMAKE_COMMAND)", "$(MAKE)"):
            return True
        return variable and has_make_control_option(tokens[1:])

    for line in content.splitlines():
        stripped = line.lstrip()
        if stripped.startswith(b"#"):
            continue
        require(
            re.match(
                rb"^(?:(?:override|export|unexport)[ \t]+)*"
                rb"(?:define|endef|undefine|-?load)(?:[ \t]|$)",
                stripped) is None,
            "private replay Make input contains a dynamic definition or "
            "load directive")
        require(not any(
                    token in line
                    for token in (
                        b"$($", b"$(${", b"${$(", b"${${",
                        b"$()", b"${}")),
                "private replay Make input contains a computed variable "
                "reference")
        if not line.startswith(b"\t") and b"=" in line:
            lhs = line.split(b"=", 1)[0]
            require(
                b"$" not in lhs or
                line == b"$(VERBOSE)MAKESILENT = -s",
                "private replay Make input contains a computed assignment")
        require(not (
                    stripped.startswith(
                        (b"export ", b"unexport ", b"override ")) and
                    b"$" in stripped),
                "private replay Make input contains a computed directive")
        require(not any(token in line for token in dangerous_functions),
                "private replay Make input contains a dynamic code or shell "
                "function")
        control = _REPLAY_SPECIAL_MAKE_VARIABLE.search(line)
        require(
            control is None or
            line in (b"SHELL = /bin/sh", canonical_cmake_assignment),
            "private replay Make input changes an execution-control "
            "variable")
        validate_value_assignment(line)
        special_references = list(
            _REPLAY_SPECIAL_MAKE_REFERENCE.finditer(line))
        require(
            not special_references or (
                line.startswith(b"\t") and
                all((match.group("paren") or match.group("brace")) in
                    {b"MAKE", b"MAKESILENT"}
                    for match in special_references) and
                b"${MAKE}" not in line and
                b"${MAKESILENT}" not in line and
                b"$(MAKE)" in line),
            "private replay Make input aliases an execution-control "
            "variable")
        recursive = re.fullmatch(
            rb"\t\$\(MAKE\) \$\(MAKESILENT\) -f "
            rb"(?P<recipe>CMakeFiles/[A-Za-z0-9_.+/-]+) "
            rb"(?P<target>[A-Za-z0-9_.+/][A-Za-z0-9_.+/-]*)",
            line)
        # Command-line MAKE is rebound to the immutable wrapper.  A literal
        # make/gmake, an alias assigned to one, or another command using Make's
        # recipe/directory switches would bypass that binding.  The canonical
        # CMake form is deliberately the sole exception.
        noncanonical_make = invokes_noncanonical_make(line)
        require(
            not noncanonical_make or recursive is not None,
            "private replay Make input contains a non-canonical recursive "
            "Make command")
        # GNU Make also accepts single-character references such as $M and
        # automatic variables such as $@.  Merely iterating over the
        # parenthesized/braced references therefore leaves unparsed dollar
        # syntax that can synthesize an executable or an unretained include.
        # Consume every dollar byte and accept only the exact generated forms
        # below; unmatched, escaped and one-character references all fail
        # closed.
        reference_offset = 0
        while True:
            reference_offset = line.find(b"$", reference_offset)
            if reference_offset < 0:
                break
            match = make_reference.match(line, reference_offset)
            require(
                match is not None,
                "private replay Make input contains a non-canonical variable "
                "reference")
            reference = match.group("paren") or match.group("brace")
            require(
                allowed_reference.fullmatch(reference) is not None,
                "private replay Make input contains a non-canonical variable "
                "reference")
            reference_offset = match.end()
        directive = re.match(
            rb"^[ \t]*(?:-?include|sinclude)(?:[ \t]|$)", line)
        if directive is not None:
            require(path.name == "build.make",
                    "private replay Make input contains an unexpected include")
            match = re.fullmatch(
                rb"include[ \t]+" +
                re.escape(parent.encode("utf-8")) +
                rb"/(?P<name>[A-Za-z0-9_.+-]+)",
                line)
            require(match is not None and
                    match.group("name") in allowed_names and
                    match.group("name") not in observed,
                    "private replay Make input contains a non-canonical "
                    "include")
            name = match.group("name")
            observed.add(name)
            decoded_name = name.decode("ascii")
            if decoded_name not in _REPLAY_DERIVED_MUTABLE_RECIPE_NAMES:
                included = parent + "/" + decoded_name
                require(included in candidates,
                        "private replay Make include is outside the retained "
                        "recipe closure")
        if b"$(MAKE)" in line:
            require(line.startswith(b"\t"),
                    "private replay Make input aliases recursive Make")
            require(recursive is not None,
                    "private replay Make input contains a non-canonical "
                    "recursive Make command")
            recipe = recursive.group("recipe").decode("ascii")
            require(
                recipe == "CMakeFiles/Makefile2" or
                (recipe.endswith("/build.make") and
                 recipe in candidates),
                "private replay recursive Make selects an unretained recipe")


def _replay_candidate_compile_closure_sha256(
    expected_closure: Sequence[Mapping[str, Any]],
) -> str:
    """Digest the candidate-normalized production compile semantics."""
    require(isinstance(expected_closure, (list, tuple)) and expected_closure,
            "private replay expected compile closure is empty")
    records = []
    outputs = set()
    for item in expected_closure:
        require(isinstance(item, Mapping) and
                item.get("role") in ("archive", "benchmark") and
                isinstance(item.get("compile_entry"), Mapping),
                "private replay expected compile closure is malformed")
        arguments = item["compile_entry"].get("normalized_arguments")
        require(isinstance(arguments, list) and len(arguments) >= 5 and
                all(isinstance(token, str) and token for token in arguments)
                and arguments[-4] == "-o" and arguments[-2] == "-c",
                "private replay expected normalized compile argv is "
                "malformed")
        output = arguments[-3]
        require(isinstance(output, str) and output not in outputs,
                "private replay expected compile closure has duplicate "
                "outputs")
        outputs.add(output)
        records.append({
            "role": item["role"],
            "output": output,
            "normalized_arguments": arguments,
        })
    encoded = json.dumps(
        sorted(records, key=lambda record: record["output"]),
        sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode(
            "ascii")
    return hashlib.sha256(encoded).hexdigest()


def _replay_candidate_target_variables_sha256(
    expected_closure: Sequence[Mapping[str, Any]],
) -> str:
    """Digest the canonical shared-variable argv for each CMake target."""
    by_target: dict[str, list[list[str]]] = {}
    for item in expected_closure:
        require(isinstance(item, Mapping) and
                isinstance(item.get("compile_entry"), Mapping),
                "private replay expected compile closure is malformed")
        arguments = item["compile_entry"].get("normalized_arguments")
        require(isinstance(arguments, list) and len(arguments) >= 5 and
                all(isinstance(token, str) and token for token in arguments)
                and arguments[-4] == "-o" and arguments[-2] == "-c",
                "private replay expected normalized compile argv is "
                "malformed")
        output = Path(arguments[-3])
        require(not output.is_absolute() and len(output.parts) >= 3 and
                output.parts[0] == "CMakeFiles" and
                output.parts[1].endswith(".dir"),
                "private replay expected compile output is non-canonical")
        target = "/".join(output.parts[:2])
        by_target.setdefault(target, []).append(arguments[1:-4])
    require(by_target,
            "private replay expected target-variable closure is empty")
    records = [{
        "target": target,
        # Source-specific CMake options are inserted around these shared
        # variables.  The target-wide member is therefore the shortest exact
        # production argv member; lexical order makes ties deterministic.
        "normalized_tokens": min(values, key=lambda item: (len(item), item)),
    } for target, values in sorted(by_target.items())]
    encoded = json.dumps(
        records, sort_keys=True, separators=(",", ":"),
        ensure_ascii=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _require_replay_candidate_compile_closure_binding(
    semantic: Mapping[str, Any],
    expected_closure: Sequence[Mapping[str, Any]],
    *,
    label: str,
) -> None:
    """Recompute the proof's semantic closure binding from the candidate."""
    require(
        semantic.get("candidate_compile_closure_sha256") ==
            _replay_candidate_compile_closure_sha256(expected_closure) and
        semantic.get("candidate_target_variables_sha256") ==
            _replay_candidate_target_variables_sha256(expected_closure),
        f"{label} replay candidate compile closure digest differs")


def _validate_replay_compile_variable_closure(
    build: Path,
    source: Path,
    compile_commands_content: bytes,
    recipe_contents: Mapping[Path, bytes],
    expected_closure: Sequence[Mapping[str, Any]],
    cmake_logical_path: str,
) -> dict[str, Any]:
    """Bind generated Make compilation inputs to the audited object closure.

    ``compile_commands.json`` describes the complete compiler argv, whereas
    CMake's Unix-Makefiles generator stores the target-wide portion in
    ``flags.make`` and may add source-specific options in ``build.make``.
    Require every linked production object's generated argv to equal its
    already-audited candidate argv, then require each executed target's three
    shared C++ variables to equal an exact retained argv member for that
    target.  Thus a post-configure edit such as ``-B``, ``-fplugin``, or
    ``-include`` cannot execute merely because it contains no shell operator.
    """
    build = build.resolve(strict=True)
    source = source.resolve(strict=True)
    require(build.is_dir() and source.is_dir(),
            "private replay semantic roots are not directories")
    _require_safe_unicode(
        cmake_logical_path, "private replay CMake identity path")
    require(isinstance(expected_closure, (list, tuple)) and expected_closure,
            "private replay expected compile closure is empty")

    expected_by_output: dict[str, list[str]] = {}
    expected_by_target: dict[str, list[list[str]]] = {}
    for item in expected_closure:
        require(isinstance(item, Mapping) and
                item.get("role") in ("archive", "benchmark"),
                "private replay expected compile closure is malformed")
        entry = item.get("compile_entry")
        require(isinstance(entry, Mapping),
                "private replay expected compile entry is malformed")
        arguments = entry.get("normalized_arguments")
        require(isinstance(arguments, list) and len(arguments) >= 5 and
                all(isinstance(token, str) and token for token in arguments),
                "private replay expected normalized compile argv is "
                "malformed")
        output_positions = [
            index for index, token in enumerate(arguments) if token == "-o"
        ]
        require(output_positions == [len(arguments) - 4] and
                arguments[-2] == "-c",
                "private replay expected compile argv has a non-canonical "
                "tail")
        output = arguments[-3]
        output_path = Path(output)
        require(not output_path.is_absolute() and
                ".." not in output_path.parts and
                len(output_path.parts) >= 3 and
                output_path.parts[0] == "CMakeFiles" and
                output_path.parts[1].endswith(".dir") and
                output_path.as_posix() == output,
                "private replay expected compile output is non-canonical")
        require(output not in expected_by_output,
                "private replay expected compile closure has duplicate "
                "outputs")
        expected_by_output[output] = list(arguments)
        target = "/".join(output_path.parts[:2])
        expected_by_target.setdefault(target, []).append(arguments[1:-4])

    commands = _strict_json_loads(
        compile_commands_content,
        "private replay compile_commands.json")
    require(isinstance(commands, list) and commands,
            "private replay compile command closure is empty")
    observed_by_output: dict[str, list[str]] = {}
    observed_tokens_by_output: dict[str, list[str]] = {}
    for raw_entry in commands:
        require(isinstance(raw_entry, Mapping),
                "private replay compile command entry is malformed")
        tokens = _compile_tokens(raw_entry)
        output = raw_entry.get("output")
        if output not in expected_by_output:
            continue
        require(isinstance(output, str),
                "private replay compile output metadata is malformed")
        require(output not in observed_by_output,
                "private replay compile command closure has a duplicate "
                "production output")
        file_value = raw_entry.get("file")
        require(isinstance(file_value, str) and file_value,
                "private replay compile command source is malformed")
        try:
            source_operand = Path(file_value).resolve(strict=True)
        except (OSError, ValueError) as error:
            raise BuildProvenanceError(
                "private replay compile command source is invalid") from error
        require(source_operand.is_relative_to(source),
                "private replay compile command source escapes its retained "
                "source root")
        normalized = _canonical_compile_argv(
            raw_entry, tokens, source_operand, build, source)
        require(normalized == expected_by_output[output],
                "private replay generated compile argv differs from the "
                "audited candidate closure")
        observed_by_output[output] = normalized
        observed_tokens_by_output[output] = tokens
    require(set(observed_by_output) == set(expected_by_output),
            "private replay compile commands omit an audited production "
            "object")

    flags_assignment = re.compile(
        rb"(?P<name>CXX_DEFINES|CXX_INCLUDES|CXX_FLAGS)"
        rb"[ \t]+=[ \t]*(?P<value>.*)")
    observed_targets: set[str] = set()
    shared_tokens_by_target: dict[str, dict[bytes, list[str]]] = {}
    for target, expected_middles in expected_by_target.items():
        flags_path = build / target / "flags.make"
        content = recipe_contents.get(flags_path)
        require(isinstance(content, bytes),
                "private replay recipe closure omits production flags.make")
        assignments: dict[bytes, list[str]] = {}
        for line in content.splitlines():
            match = flags_assignment.fullmatch(line)
            if match is None:
                continue
            name = match.group("name")
            require(name not in assignments,
                    "private replay production flags.make duplicates a "
                    "C++ variable")
            value = match.group("value")
            require(not any(ord(character) in value
                            for character in SHELL_RECIPE_META),
                    "private replay production flags.make contains shell "
                    "code or indirection")
            try:
                decoded = value.decode("utf-8", errors="strict")
                _require_safe_unicode(
                    decoded, "private replay production flags.make value")
                tokens = shlex.split(decoded, posix=True)
            except (UnicodeDecodeError, ValueError) as error:
                raise BuildProvenanceError(
                    "private replay production flags.make has "
                    f"non-canonical quoting: {error}") from error
            require(all(token and not token.startswith("@")
                        for token in tokens),
                    "private replay production flags.make contains an empty "
                    "operand or response-file indirection")
            assignments[name] = tokens
        require(set(assignments) == {
                    b"CXX_DEFINES", b"CXX_INCLUDES", b"CXX_FLAGS"},
                "private replay production flags.make omits a canonical "
                "C++ variable")
        shared_tokens = (
            assignments[b"CXX_DEFINES"] +
            assignments[b"CXX_INCLUDES"] +
            assignments[b"CXX_FLAGS"])
        normalized_shared = _normalize_build_argv(
            shared_tokens, build, source)
        canonical_shared = min(
            expected_middles, key=lambda item: (len(item), item))
        require(normalized_shared == canonical_shared,
                "private replay production flags.make variables differ from "
                "the audited compile closure")
        shared_tokens_by_target[target] = assignments
        observed_targets.add(target)
    require(observed_targets == set(expected_by_target),
            "private replay did not bind every production flags target")

    for output, compile_tokens in observed_tokens_by_output.items():
        output_path = Path(output)
        target = "/".join(output_path.parts[:2])
        build_make_path = build / target / "build.make"
        build_make = recipe_contents.get(build_make_path)
        require(isinstance(build_make, bytes),
                "private replay recipe closure omits production build.make")
        output_bytes = output.encode("utf-8", errors="strict")
        recipe_lines = []
        lines = build_make.splitlines()
        rule_indices = [
            index for index, line in enumerate(lines)
            if line.startswith(output_bytes + b":")
        ]
        expected_rules = [
            output_bytes + b": " +
                (target + "/flags.make").encode("utf-8"),
            output_bytes + b": " +
                compile_tokens[-1].encode("utf-8", errors="strict"),
            output_bytes + b": " +
                (target + "/compiler_depend.ts").encode("utf-8"),
        ]
        require([lines[index] for index in rule_indices] == expected_rules,
                "private replay production object dependency rules differ "
                "from canonical CMake output")
        for index in rule_indices:
            cursor = index + 1
            while cursor < len(lines) and lines[cursor].startswith(b"\t"):
                recipe_lines.append(lines[cursor])
                cursor += 1
        require(len(recipe_lines) == 2,
                "private replay production object rule has a non-canonical "
                "recipe count")
        echo_pattern = re.compile(
            rb"\t@\$\(CMAKE_COMMAND\) -E cmake_echo_color "
            rb"\"--switch=\$\(COLOR\)\" --green --progress-dir=" +
            re.escape(str(build).encode("utf-8")) +
            rb"/CMakeFiles --progress-num=\$\(CMAKE_PROGRESS_[0-9]+\) "
            rb"\"Building CXX object " + re.escape(output_bytes) + rb"\"")
        require(echo_pattern.fullmatch(recipe_lines[0]) is not None,
                "private replay production object status recipe differs "
                "from canonical CMake output")
        try:
            command = recipe_lines[1][1:].decode(
                "utf-8", errors="strict")
            template_tokens = shlex.split(command, posix=True)
        except (UnicodeDecodeError, ValueError) as error:
            raise BuildProvenanceError(
                "private replay production compile recipe has "
                f"non-canonical quoting: {error}") from error
        references = (
            "$(CXX_DEFINES)", "$(CXX_INCLUDES)", "$(CXX_FLAGS)")
        require([token for token in template_tokens
                 if token in references] == list(references) and
                all(token in references or
                    ("$" not in token and
                     not any(character in token
                             for character in SHELL_RECIPE_META))
                    for token in template_tokens),
                "private replay production compile recipe contains a "
                "non-canonical variable or shell operand")
        shared = shared_tokens_by_target[target]
        expansion = {
            "$(CXX_DEFINES)": shared[b"CXX_DEFINES"],
            "$(CXX_INCLUDES)": shared[b"CXX_INCLUDES"],
            "$(CXX_FLAGS)": shared[b"CXX_FLAGS"],
        }
        expanded = []
        for token in template_tokens:
            if token in expansion:
                expanded.extend(expansion[token])
            else:
                expanded.append(token)
        expected_recipe = (
            compile_tokens[:-4] +
            ["-MD", "-MT", output, "-MF", output + ".d"] +
            compile_tokens[-4:])
        require(expanded == expected_recipe,
                "private replay production compile recipe differs from the "
                "audited compile command")

    control_expected = {
        b"CMAKE_COMMAND": cmake_logical_path.encode(
            "utf-8", errors="strict"),
        b"CMAKE_BINARY_DIR": str(build).encode(
            "utf-8", errors="strict"),
        b"CMAKE_SOURCE_DIR": str(source).encode(
            "utf-8", errors="strict"),
    }
    control_assignment = re.compile(
        rb"(?P<name>CMAKE_COMMAND|CMAKE_BINARY_DIR|CMAKE_SOURCE_DIR)"
        rb"[ \t]+=[ \t]*(?P<value>.*)")
    control_counts = Counter()
    for content in recipe_contents.values():
        for line in content.splitlines():
            match = control_assignment.fullmatch(line)
            if match is None:
                continue
            name = match.group("name")
            require(match.group("value") == control_expected[name],
                    "private replay Make control variable differs from its "
                    "retained configure identity")
            control_counts[name] += 1
    require(all(control_counts[name] > 0 for name in control_expected),
            "private replay recipe closure omits a required CMake control "
            "binding")

    cache_content = recipe_contents.get(build / "CMakeCache.txt")
    require(isinstance(cache_content, bytes),
            "private replay recipe closure omits CMakeCache.txt")
    replay_cache = parse_cmake_cache(cache_content)
    require(replay_cache.get("CMAKE_HOME_DIRECTORY") == str(source),
            "private replay CMake cache points at another source root")
    return {
        "schema": "leopard2-replay-compile-variable-closure/v2",
        "compile_commands_sha256":
            hashlib.sha256(compile_commands_content).hexdigest(),
        "candidate_compile_closure_sha256":
            _replay_candidate_compile_closure_sha256(expected_closure),
        "candidate_target_variables_sha256":
            _replay_candidate_target_variables_sha256(expected_closure),
        "object_count": len(expected_by_output),
        "target_count": len(expected_by_target),
        "control_assignment_counts": {
            name.decode("ascii"): control_counts[name]
            for name in sorted(control_counts)
        },
    }


_CANONICAL_REPLAY_PLAN_RELATIVE = (
    "CMakeFiles/leopard2-canonical-replay.make")
_CANONICAL_REPLAY_ATTESTATION_SOURCE = (
    "cmake/GenerateBenchmarkSourceAttestation.cmake")
_CANONICAL_REPLAY_ATTESTATION_TARGET = (
    "leopard2-canonical-source-attestation")
_CANONICAL_REPLAY_FORCE_TARGET = "leopard2-canonical-force"
_CANONICAL_REPLAY_FREEZE_TARGET = "leopard2-canonical-freeze"


def _canonical_replay_source_record(
    candidate: Mapping[str, Any], relative: str,
) -> Mapping[str, Any]:
    """Return one exact candidate-bound tracked-source record."""
    relative_path = Path(relative)
    require(
        isinstance(relative, str) and relative and
        not relative_path.is_absolute() and
        relative_path.as_posix() == relative and
        all(component not in ("", ".", "..")
            for component in relative_path.parts),
        f"canonical replay tracked source path is unsafe: {relative!r}")
    manifest = candidate.get("tracked_source_manifest")
    require(isinstance(manifest, Mapping) and
            isinstance(manifest.get("files"), list),
            "canonical replay candidate source manifest is malformed")
    matches = [
        record for record in manifest["files"]
        if isinstance(record, Mapping) and record.get("path") == relative
    ]
    require(len(matches) == 1,
            f"canonical replay source closure does not contain {relative!r}")
    record = matches[0]
    require(set(record) == {"path", "sha256", "size", "mode"} and
            isinstance(record["sha256"], str) and
            re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is not None and
            type(record["size"]) is int and
            0 <= record["size"] <= MAX_TRACKED_SOURCE_FILE_BYTES and
            type(record["mode"]) is int and
            0 <= record["mode"] <= 0o7777,
            f"canonical replay tracked source record is malformed: "
            f"{relative!r}")
    return record


def _restore_canonical_replay_token(
    token: str, build: Path, source: Path, *,
    execution_build_root: str | None = None,
) -> str:
    """Materialize one candidate-normalized argv token in the private replay."""
    require(isinstance(token, str) and token and
            "\0" not in token and "\r" not in token and "\n" not in token,
            "canonical replay argv token is malformed")
    materialized_build = (
        str(build) if execution_build_root is None else execution_build_root)
    require(materialized_build and "\0" not in materialized_build,
            "canonical replay execution build root is malformed")
    restored = token.replace("${BUILD_ROOT}", materialized_build).replace(
        "${SOURCE_ROOT}", str(source))
    require("${BUILD_ROOT}" not in restored and
            "${SOURCE_ROOT}" not in restored and "$" not in restored,
            "canonical replay argv retains an unbound placeholder")
    _require_safe_unicode(restored, "canonical replay argv token")
    return restored


def _canonical_replay_build_operand(
    value: str, build: Path, source: Path,
) -> str:
    """Return one canonical build-relative target from normalized argv data."""
    restored = _restore_canonical_replay_token(value, build, source)
    operand = Path(restored)
    require(build.is_absolute(),
            "canonical replay build root is not absolute")
    lexical = Path(os.path.abspath(os.fspath(
        operand if operand.is_absolute() else build / operand)))
    require(lexical.is_relative_to(build),
            "canonical replay build operand escapes its private build")
    relative = lexical.relative_to(build).as_posix()
    require(relative and ".." not in Path(relative).parts and
            re.fullmatch(r"[A-Za-z0-9_.+/-]+", relative) is not None,
            "canonical replay build operand is not Make-safe")
    return relative


def _canonical_replay_output_paths(
    candidate: Mapping[str, Any], source: Path, build: Path,
) -> list[Path]:
    """Derive every file path written by the sealed canonical replay."""
    closure = candidate.get("source_object_compile_closure")
    archive_commands = candidate.get("archive_link_commands")
    target = candidate.get("executable_target")
    require(isinstance(closure, list) and closure and
            isinstance(archive_commands, list) and archive_commands and
            isinstance(archive_commands[0], list) and
            len(archive_commands[0]) >= 3 and
            isinstance(target, str) and target,
            "canonical replay output topology input is malformed")
    relative_outputs = []
    for item in closure:
        require(isinstance(item, Mapping) and
                isinstance(item.get("compile_entry"), Mapping),
                "canonical replay output compile record is malformed")
        arguments = item["compile_entry"].get("normalized_arguments")
        require(isinstance(arguments, list) and len(arguments) >= 4 and
                all(isinstance(argument, str) for argument in arguments),
                "canonical replay output compile argv is malformed")
        relative_outputs.append(_canonical_replay_build_operand(
            arguments[-3], build, source))
    relative_outputs.extend((
        _canonical_replay_build_operand(
            str(archive_commands[0][2]), build, source),
        _canonical_replay_build_operand(target, build, source),
        _CANONICAL_REPLAY_PLAN_RELATIVE,
        "generated/leopard2-benchmark-attestation/"
        "leopard2_benchmark_source_attestation.h",
        "generated/leopard2-benchmark-attestation/"
        "leopard2_benchmark_source_attestation.h.lock",
    ))
    require(len(relative_outputs) == len(set(relative_outputs)),
            "canonical replay output paths are duplicated")
    return [build / relative for relative in sorted(relative_outputs)]


def _canonical_replay_attestation_header_bytes(
    candidate: Mapping[str, Any],
) -> bytes:
    """Render the exact unchanged header expected from the retained Git state."""
    manifest = candidate.get("tracked_source_manifest")
    git = manifest.get("git") if isinstance(manifest, Mapping) else None
    require(
        isinstance(git, Mapping) and
        isinstance(git.get("commit"), str) and
        re.fullmatch(r"[0-9a-f]+", git["commit"]) is not None and
        isinstance(git.get("tree"), str) and
        re.fullmatch(r"[0-9a-f]+", git["tree"]) is not None and
        type(git.get("dirty")) is bool,
        "canonical replay source identity is malformed")
    dirty = 1 if git["dirty"] else 0
    return (
        "#ifndef LEOPARD2_BENCHMARK_SOURCE_ATTESTATION_GENERATED_H\n"
        "#define LEOPARD2_BENCHMARK_SOURCE_ATTESTATION_GENERATED_H\n"
        "\n"
        "#undef LEO2_BENCHMARK_SOURCE_COMMIT\n"
        "#undef LEO2_BENCHMARK_SOURCE_TREE\n"
        "#undef LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY\n"
        f"#define LEO2_BENCHMARK_SOURCE_COMMIT \"{git['commit']}\"\n"
        f"#define LEO2_BENCHMARK_SOURCE_TREE \"{git['tree']}\"\n"
        f"#define LEO2_BENCHMARK_SOURCE_TRACKED_DIRTY {dirty}\n"
        "\n"
        "#endif\n"
    ).encode("ascii")


def _validate_sealed_replay_artifact_bindings(
    candidate: Mapping[str, Any],
    artifacts: Mapping[str, Mapping[str, Any]],
    source: Path,
    build: Path,
    *,
    label: str,
) -> None:
    """Bind every private-output memfd to the independently known result."""
    require(source.is_absolute() and build.is_absolute() and
            isinstance(artifacts, Mapping) and artifacts,
            f"{label} sealed replay artifact binding input is malformed")
    expected: dict[str, tuple[str, int, int]] = {}

    def add(
        relative: str, sha256: Any, size: Any, permissions: int,
    ) -> None:
        require(
            relative not in expected and
            isinstance(sha256, str) and
            re.fullmatch(r"[0-9a-f]{64}", sha256) is not None and
            type(size) is int and 0 <= size <= MAX_FILE_BYTES,
            f"{label} expected sealed replay artifact identity is malformed")
        expected[relative] = (sha256, size, permissions)

    closure = candidate.get("source_object_compile_closure")
    require(isinstance(closure, list) and closure,
            f"{label} sealed replay object closure is malformed")
    for item in closure:
        require(isinstance(item, Mapping) and
                isinstance(item.get("compile_entry"), Mapping) and
                isinstance(item.get("object"), Mapping),
                f"{label} sealed replay object identity is malformed")
        arguments = item["compile_entry"].get("normalized_arguments")
        require(isinstance(arguments, list) and len(arguments) >= 4 and
                all(isinstance(argument, str) for argument in arguments),
                f"{label} sealed replay compile output is malformed")
        relative = _canonical_replay_build_operand(
            arguments[-3], build, source)
        identity = item["object"]
        add(relative, identity.get("sha256"), identity.get("size"), 0o600)

    archive = candidate.get("archive")
    executable = candidate.get("executable")
    target = candidate.get("executable_target")
    require(isinstance(archive, Mapping) and
            isinstance(executable, Mapping) and
            isinstance(target, str) and target,
            f"{label} sealed replay linked-artifact identity is malformed")
    add("libleopard.a", archive.get("sha256"), archive.get("size"), 0o600)
    add(target, executable.get("sha256"), executable.get("size"), 0o700)

    header_relative = (
        "generated/leopard2-benchmark-attestation/"
        "leopard2_benchmark_source_attestation.h")
    header = _canonical_replay_attestation_header_bytes(candidate)
    add(
        header_relative, hashlib.sha256(header).hexdigest(),
        len(header), 0o600)
    add(
        header_relative + ".lock", hashlib.sha256(b"").hexdigest(),
        0, 0o600)

    require(set(artifacts) == set(expected),
            f"{label} sealed replay artifact path closure differs")
    for relative, (sha256, size, permissions) in expected.items():
        artifact = artifacts[relative]
        require(
            isinstance(artifact, Mapping) and
            artifact.get("path") == relative and
            artifact.get("sha256") == sha256 and
            artifact.get("size") == size and
            type(artifact.get("mode")) is int and
            stat.S_ISREG(artifact["mode"]) and
            stat.S_IMODE(artifact["mode"]) == permissions,
            f"{label} sealed replay artifact differs from its independent "
            f"identity: {relative}")


def _retain_canonical_replay_output_topology(
    candidate: Mapping[str, Any], source: Path, build: Path, stack: ExitStack,
) -> dict[str, Any]:
    """Reserve and retain every exact replay output inode before execution."""
    outputs = _canonical_replay_output_paths(candidate, source, build)
    tree = _OpenDirectoryTree(
        build, stack, "canonical replay output topology")
    parents: dict[Path, int] = {}
    for output in outputs:
        relative = output.relative_to(build)
        parent_parts = relative.parts[:-1]
        parent_descriptor = tree.directory(parent_parts)
        parent = output.parent
        parents[parent] = parent_descriptor

    guard = _InotifyMutationGuard("canonical replay output topology")
    entry_guard = _InotifyMutationGuard(
        "canonical replay exact output entries")
    try:
        for parent in sorted(parents, key=str):
            guard.add_directory_path(parent)
            entry_guard.add_exact_directory_entries(
                parent, [output.name for output in outputs
                         if output.parent == parent])

        plan = build / _CANONICAL_REPLAY_PLAN_RELATIVE
        header = (
            build / "generated/leopard2-benchmark-attestation/"
            "leopard2_benchmark_source_attestation.h")
        archive = build / "libleopard.a"
        executable = build / str(candidate.get("executable_target", ""))
        initial_content = {
            archive: b"!<arch>\n",
            header: _canonical_replay_attestation_header_bytes(candidate),
        }
        retained_outputs = []
        for output in outputs:
            relative = output.relative_to(build)
            parent_descriptor = parents[output.parent]
            name = relative.parts[-1]
            descriptor = -1
            retained_file = None
            try:
                if output == plan:
                    descriptor = os.open(
                        name,
                        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                        getattr(os, "O_NONBLOCK", 0) |
                        getattr(os, "O_NOFOLLOW", 0),
                        dir_fd=parent_descriptor)
                    retained_file = os.fdopen(
                        descriptor, "rb", buffering=0)
                else:
                    mode = 0o700 if output == executable else 0o600
                    descriptor = os.open(
                        name,
                        os.O_RDWR | os.O_CREAT | os.O_EXCL |
                        getattr(os, "O_CLOEXEC", 0) |
                        getattr(os, "O_NOFOLLOW", 0),
                        mode, dir_fd=parent_descriptor)
                    retained_file = os.fdopen(
                        descriptor, "r+b", buffering=0)
                    content = initial_content.get(output, b"")
                    view = memoryview(content)
                    while view:
                        written = os.write(descriptor, view)
                        require(written > 0,
                                f"canonical replay output reservation write "
                                f"stalled: {output}")
                        view = view[written:]
                    os.fsync(descriptor)
                    entry_guard.accept_exact_regular_creation(name)
                metadata = os.fstat(descriptor)
                lexical = os.stat(
                    name, dir_fd=parent_descriptor,
                    follow_symlinks=False)
                require(
                    stat.S_ISREG(metadata.st_mode) and
                    metadata.st_nlink == 1 and
                    (metadata.st_dev, metadata.st_ino,
                     stat.S_IFMT(metadata.st_mode)) ==
                    (lexical.st_dev, lexical.st_ino,
                     stat.S_IFMT(lexical.st_mode)),
                    f"canonical replay output is a symlink, hard link, or "
                    f"non-regular file: {output}")
                stack.callback(retained_file.close)
                retained_outputs.append({
                    "path": output,
                    "parent_descriptor": parent_descriptor,
                    "name": name,
                    "file": retained_file,
                    "identity": (
                        metadata.st_dev, metadata.st_ino,
                        stat.S_IFMT(metadata.st_mode),
                        metadata.st_uid, metadata.st_gid, metadata.st_nlink,
                    ),
                })
            except BaseException as primary:
                cleanup_failure: BaseException | None = None
                try:
                    if retained_file is not None:
                        retained_file.close()
                    elif descriptor >= 0:
                        consumed = [False]
                        try:
                            _close_raw_descriptor_with_precedence(
                                descriptor, None,
                                "canonical replay output raw descriptor",
                                "cannot close canonical replay output raw "
                                "descriptor",
                                ownership_consumed=consumed)
                        finally:
                            if consumed[0]:
                                descriptor = -1
                except BaseException as error:
                    cleanup_failure = error
                if cleanup_failure is not None:
                    _raise_owner_exit_failure(
                        primary, cleanup_failure,
                        "canonical replay output reservation",
                        "canonical replay output reservation cleanup failed: "
                        f"{cleanup_failure}; primary failure: "
                        f"{type(primary).__name__}: {primary}")
                if not isinstance(primary, OSError):
                    raise
                raise BuildProvenanceError(
                    f"canonical replay output is preexisting, unstable or "
                    f"unsafe: {output}: {primary}") from primary
    except BaseException:
        guard._close_without_verification()
        entry_guard._close_without_verification()
        raise
    stack.callback(guard.close)
    stack.callback(entry_guard.close)
    retained = {
        "guard": guard,
        "entry_guard": entry_guard,
        "parents": [
            (parent, parents[parent]) for parent in sorted(parents, key=str)
        ],
        "outputs": retained_outputs,
    }
    _verify_canonical_replay_output_topology(retained, require_outputs=False)
    return retained


def _verify_canonical_replay_output_topology(
    retained: Mapping[str, Any], *, require_outputs: bool,
) -> None:
    """Verify retained output parents and optional completed regular outputs."""
    guard = retained.get("guard")
    parents = retained.get("parents")
    entry_guard = retained.get("entry_guard")
    outputs = retained.get("outputs")
    require(isinstance(guard, _InotifyMutationGuard) and
            isinstance(entry_guard, _InotifyMutationGuard) and
            isinstance(parents, list) and isinstance(outputs, list),
            "canonical replay retained output topology is malformed")
    guard.verify()
    entry_guard.verify()
    for parent, descriptor in parents:
        try:
            lexical = os.stat(parent, follow_symlinks=False)
            retained_status = os.fstat(descriptor)
        except OSError as error:
            raise BuildProvenanceError(
                f"canonical replay output parent changed: {parent}: {error}") \
                from error
        require(
            stat.S_ISDIR(lexical.st_mode) and
            (lexical.st_dev, lexical.st_ino, stat.S_IFMT(lexical.st_mode)) ==
            (retained_status.st_dev, retained_status.st_ino,
             stat.S_IFMT(retained_status.st_mode)),
            f"canonical replay output parent changed: {parent}")
    for record in outputs:
        require(isinstance(record, Mapping) and
                set(record) == {
                    "path", "parent_descriptor", "name", "file", "identity",
                },
                "canonical replay retained output record is malformed")
        output = record["path"]
        parent_descriptor = record["parent_descriptor"]
        name = record["name"]
        retained_file = record["file"]
        require(isinstance(output, Path) and
                type(parent_descriptor) is int and
                isinstance(name, str) and name and
                hasattr(retained_file, "fileno") and
                isinstance(record["identity"], tuple),
                "canonical replay retained output identity is malformed")
        try:
            metadata = os.stat(
                name, dir_fd=parent_descriptor, follow_symlinks=False)
            retained_status = os.fstat(retained_file.fileno())
        except OSError as error:
            raise BuildProvenanceError(
                f"canonical replay output is missing or unsafe: "
                f"{output}: {error}") from error
        identity = (
            retained_status.st_dev, retained_status.st_ino,
            stat.S_IFMT(retained_status.st_mode),
            retained_status.st_uid, retained_status.st_gid,
            retained_status.st_nlink,
        )
        lexical_identity = (
            metadata.st_dev, metadata.st_ino, stat.S_IFMT(metadata.st_mode),
            metadata.st_uid, metadata.st_gid, metadata.st_nlink,
        )
        require(
                identity == record["identity"] and
                lexical_identity == record["identity"] and
                stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
                f"canonical replay output is a symlink, hard link, or "
                f"non-regular/replaced file: {output}")
    guard.verify()
    entry_guard.verify()


def _canonical_replay_makefile_bytes(
    candidate: Mapping[str, Any],
    source: Path,
    build: Path,
    transports: Mapping[str, str],
    artifact_transports: Mapping[str, str] | None = None,
    execution_build_root: str | None = None,
) -> tuple[bytes, dict[str, int]]:
    """Render the sole Make input used by a clean compiler replay.

    CMake is still run in the private empty workspace so its effective
    configuration is reproduced.  None of its generated Make/CMake programs
    are executed, however.  This plan is derived only from the already-audited
    candidate compile/archive/link semantics and the one tracked attestation
    script.  Consequently its exact bytes can be regenerated by the offline
    verifier rather than authenticated by a self-rooted inventory digest.
    """
    require(source.is_absolute() and build.is_absolute() and
            source != build,
            "canonical replay roots are malformed")
    if execution_build_root is not None:
        require(
            re.fullmatch(
                r"/proc/self/fd/[0-9]+", execution_build_root) is not None,
            "canonical replay execution root is malformed")
    required_transports = {"cxx", "archiver", "ranlib", "cmake", "git"}
    if artifact_transports is not None:
        required_transports.update(("copy", "inventory"))
    require(isinstance(transports, Mapping) and
            set(transports) == required_transports,
            "canonical replay transport roles differ")
    for label, value in transports.items():
        require(isinstance(value, str) and
                re.fullmatch(r"/proc/self/fd/[0-9]+", value) is not None,
                f"canonical replay {label} transport is malformed")
    if artifact_transports is not None:
        require(
            isinstance(artifact_transports, Mapping) and
            artifact_transports and
            all(
                isinstance(relative, str) and relative and
                not Path(relative).is_absolute() and
                all(component not in ("", ".", "..")
                    for component in Path(relative).parts) and
                isinstance(destination, str) and
                re.fullmatch(r"/proc/self/fd/[0-9]+", destination)
                    is not None
                for relative, destination in artifact_transports.items()),
            "canonical replay artifact transports are malformed")

    target = candidate.get("executable_target")
    require(isinstance(target, str) and
            re.fullmatch(r"[A-Za-z0-9_.+][A-Za-z0-9_.+-]*", target)
            is not None,
            "canonical replay executable target is malformed")
    closure = candidate.get("source_object_compile_closure")
    require(isinstance(closure, list) and closure,
            "canonical replay compile closure is empty")
    cache = candidate.get("validated_cache")
    require(isinstance(cache, Mapping) and
            isinstance(cache.get("CMAKE_CXX_COMPILER"), str) and
            cache["CMAKE_CXX_COMPILER"],
            "canonical replay candidate cache is malformed")
    compiler_identity = candidate.get("compiler")
    compiler_path_value = (
        compiler_identity.get("path")
        if isinstance(compiler_identity, Mapping) else None)
    if not (isinstance(compiler_path_value, str) and
            Path(compiler_path_value).is_absolute()):
        # Current production records retain the resolved compiler artifact,
        # while historical and adversarial replay fixtures predate that
        # redundant identity.  Their already-validated cache still binds the
        # exact absolute driver used by every restored compile command.  An
        # ambiguous spelling such as /usr/bin/compiler remains fail-closed in
        # the GNU-only option classifier.
        compiler_path_value = cache["CMAKE_CXX_COMPILER"]
    require(Path(compiler_path_value).is_absolute(),
            "canonical replay compiler path is malformed")
    compiler_path = Path(compiler_path_value)

    prepared_records: list[
        tuple[Mapping[str, Any], str, list[str], str, Path]] = []
    outputs: set[str] = set()
    archive_outputs: set[str] = set()
    benchmark_outputs: list[str] = []
    for item in closure:
        require(isinstance(item, Mapping) and
                item.get("role") in ("archive", "benchmark") and
                isinstance(item.get("compile_entry"), Mapping),
                "canonical replay compile record is malformed")
        arguments = item["compile_entry"].get("normalized_arguments")
        require(isinstance(arguments, list) and len(arguments) >= 5 and
                all(isinstance(argument, str) and argument
                    for argument in arguments) and
                arguments[-4] == "-o" and arguments[-2] == "-c",
                "canonical replay compile argv is malformed")
        restored = [
            _restore_canonical_replay_token(argument, build, source)
            for argument in arguments
        ]
        _require_shell_literal_tokens(
            restored, "canonical replay raw compile command")
        require(restored[0] == cache["CMAKE_CXX_COMPILER"],
                "canonical replay compile command uses another driver")
        output = _canonical_replay_build_operand(
            arguments[-3], build, source)
        require(output not in outputs,
                "canonical replay compile outputs are duplicated")
        outputs.add(output)
        role = str(item["role"])
        if role == "archive":
            archive_outputs.add(output)
        else:
            benchmark_outputs.append(output)
        source_operand = Path(restored[-1])
        lexical_source = Path(os.path.abspath(os.fspath(source_operand)))
        require(source_operand.is_absolute() and
                lexical_source == source_operand and
                source_operand.is_relative_to(source),
                "canonical replay compile source escapes its private source")
        relative_source = source_operand.relative_to(source).as_posix()
        require(
            relative_source and
            all(component not in ("", ".", "..")
                for component in Path(relative_source).parts),
            "canonical replay compile source path is unsafe")
        require(
            (role == "benchmark") ==
            (relative_source == "bench/leopard2/benchmark.cpp"),
            "canonical replay compile source role differs")
        _canonical_replay_source_record(candidate, relative_source)
        prepared_records.append(
            (item, output, restored, role, source_operand))
    require(len(benchmark_outputs) == 1,
            "canonical replay must compile exactly one benchmark object")

    library_sources = {
        source_operand for _item, _output, _restored, role, source_operand
        in prepared_records if role == "archive"
    }
    require(library_sources,
            "canonical replay archive compile source closure is empty")
    compile_records: list[tuple[str, list[str], str]] = []
    for item, output, restored, role, source_operand in prepared_records:
        profile = _validate_compile_flags(
            restored, source_operand, source_root=source, build_root=build,
            cache=cache, library_sources=library_sources,
            benchmark_source=(role == "benchmark"),
            lexical_build_output=True, compiler_path=compiler_path)
        require(item.get("flag_profile") == profile,
                "canonical replay compile flag profile differs")
        arguments = item["compile_entry"]["normalized_arguments"]
        restored = [
            _restore_canonical_replay_token(
                argument, build, source,
                execution_build_root=execution_build_root)
            for argument in arguments
        ]
        restored[0] = transports["cxx"]
        # The canonical plan does not consume generated dependency files.
        # Omitting -MD/-MT/-MF avoids an otherwise unused writable pathname
        # and does not change compiler object bytes.
        recipe = restored
        _require_shell_literal_tokens(
            recipe, f"canonical replay compile recipe for {output}")
        compile_records.append((output, recipe, role))

    archive_commands = candidate.get("archive_link_commands")
    require(isinstance(archive_commands, list) and
            len(archive_commands) == 2 and
            all(isinstance(command, list) and command and
                all(isinstance(argument, str) and argument
                    for argument in command)
                for command in archive_commands),
            "canonical replay archive recipe is malformed")
    archive_command = [
        _restore_canonical_replay_token(
            argument, build, source,
            execution_build_root=execution_build_root)
        for argument in archive_commands[0]
    ]
    index_command = [
        _restore_canonical_replay_token(
            argument, build, source,
            execution_build_root=execution_build_root)
        for argument in archive_commands[1]
    ]
    _require_shell_literal_tokens(
        archive_command, "canonical replay raw archive command")
    _require_shell_literal_tokens(
        index_command, "canonical replay raw archive index command")
    require(len(archive_command) >= 4 and
            archive_command[1] == "qc" and
            len(index_command) == 2,
            "canonical replay archive recipe shape differs")
    archive_target = _canonical_replay_build_operand(
        archive_command[2], build, source)
    archive_dependencies = [
        _canonical_replay_build_operand(value, build, source)
        for value in archive_command[3:]
    ]
    require(archive_target == "libleopard.a" and
            len(archive_dependencies) == len(set(archive_dependencies)) and
            set(archive_dependencies) == archive_outputs and
            _canonical_replay_build_operand(
                index_command[1], build, source) == archive_target,
            "canonical replay archive object/target closure differs")
    require(isinstance(cache.get("CMAKE_AR"), str) and
            isinstance(cache.get("CMAKE_RANLIB"), str),
            "canonical replay archive tool cache is malformed")
    archive_recipe = (
        shlex.join(archive_command) + "\n" +
        shlex.join(index_command) + "\n").encode("utf-8", errors="strict")
    exact_archive_objects, _exact_archive_commands = \
        _archive_recipe_semantics(
            archive_recipe, build, build / archive_target,
            Path(cache["CMAKE_AR"]), Path(cache["CMAKE_RANLIB"]),
            lexical_build_operands=True)
    require(
        exact_archive_objects ==
        [build / dependency for dependency in archive_dependencies],
        "canonical replay archive recipe semantic closure differs")
    archive_command[0] = transports["archiver"]
    index_command[0] = transports["ranlib"]
    _require_shell_literal_tokens(
        archive_command, "canonical replay archive command")
    _require_shell_literal_tokens(
        index_command, "canonical replay archive index command")

    raw_link = candidate.get("executable_link_command")
    require(isinstance(raw_link, list) and raw_link and
            all(isinstance(argument, str) and argument
                for argument in raw_link),
            "canonical replay executable link command is malformed")
    link_command = [
        _restore_canonical_replay_token(
            argument, build, source,
            execution_build_root=execution_build_root)
        for argument in raw_link
    ]
    _require_shell_literal_tokens(
        link_command, "canonical replay raw executable link command")
    output_positions = [
        index for index, argument in enumerate(link_command)
        if argument == "-o"
    ]
    require(output_positions == [len(link_command) - 5] and
            output_positions[0] + 1 < len(link_command) and
            _canonical_replay_build_operand(
                link_command[output_positions[0] + 1],
                build, source) == target and
            archive_target in link_command,
            "canonical replay executable link closure differs")
    benchmark_output = benchmark_outputs[0]
    require(benchmark_output in {
                _canonical_replay_build_operand(argument, build, source)
                for argument in link_command
                if (argument.endswith((".o", ".obj")) and
                    (not Path(argument).is_absolute() or
                     Path(argument).resolve(strict=False).is_relative_to(
                         build)))
            },
            "canonical replay executable omits its benchmark object")
    exact_link_objects, _exact_link_command = \
        _executable_recipe_semantics(
            (shlex.join(link_command) + "\n").encode(
                "utf-8", errors="strict"),
            build, build / target, build / archive_target,
            Path(cache["CMAKE_CXX_COMPILER"]),
            lexical_build_operands=True)
    require(exact_link_objects == [build / benchmark_output],
            "canonical replay executable semantic closure differs")
    link_command[0] = transports["cxx"]
    _require_shell_literal_tokens(
        link_command, "canonical replay executable link command")

    script_record = _canonical_replay_source_record(
        candidate, _CANONICAL_REPLAY_ATTESTATION_SOURCE)
    del script_record
    script = source / _CANONICAL_REPLAY_ATTESTATION_SOURCE
    materialized_build = (
        str(build) if execution_build_root is None else execution_build_root)
    header = (
        materialized_build +
        "/generated/leopard2-benchmark-attestation/"
        "leopard2_benchmark_source_attestation.h")
    attestation_command = [
        transports["cmake"], "-E", "env",
        "GIT_CONFIG_GLOBAL=/dev/null",
        "GIT_CONFIG_NOSYSTEM=1",
        "GIT_CONFIG_SYSTEM=/dev/null",
        "GIT_NO_REPLACE_OBJECTS=1",
        "GIT_OPTIONAL_LOCKS=0",
        transports["cmake"],
        f"-DLEO2_SOURCE_DIR={source}",
        f"-DLEO2_OUTPUT_FILE={header}",
        f"-DLEO2_GIT_EXECUTABLE={transports['git']}",
        "-P", str(script),
    ]
    _require_shell_literal_tokens(
        attestation_command, "canonical replay attestation command")
    def recipe(tokens: Sequence[str]) -> str:
        rendered = shlex.join(list(tokens))
        require("$" not in rendered and "\n" not in rendered and
                "\r" not in rendered,
                "canonical replay recipe rendering is unsafe")
        return "\t" + rendered

    selected = f"CMakeFiles/{target}.dir/replay"
    selected_dependency = (
        _CANONICAL_REPLAY_FREEZE_TARGET
        if artifact_transports is not None else target)
    lines = [
        "# Leopard2 canonical candidate-bound compiler replay.",
        "# Generated by tools/leopard2_build_provenance.py; do not edit.",
        ".DELETE_ON_ERROR:",
        ".SUFFIXES:",
        "",
        f".PHONY: {selected} {_CANONICAL_REPLAY_ATTESTATION_TARGET} "
        f"{_CANONICAL_REPLAY_FORCE_TARGET}" +
        (f" {_CANONICAL_REPLAY_FREEZE_TARGET}"
         if artifact_transports is not None else ""),
        f"{selected}: {selected_dependency}",
        "",
        f"{_CANONICAL_REPLAY_FORCE_TARGET}:",
        "",
        f"{_CANONICAL_REPLAY_ATTESTATION_TARGET}:",
        recipe(attestation_command),
        "",
    ]
    for output, command, role in sorted(
            compile_records, key=lambda record: record[0]):
        # The source path is already carried as one shell-quoted compiler argv
        # token and the complete tracked tree is retained independently.
        # Do not parse that pathname as Make syntax: TMPDIR/source paths may
        # legitimately contain spaces or other Make-special separators.
        prerequisites = [_CANONICAL_REPLAY_FORCE_TARGET]
        if role == "benchmark":
            prerequisites.append(_CANONICAL_REPLAY_ATTESTATION_TARGET)
        lines.append(f"{output}: {' '.join(prerequisites)}")
        lines.append(recipe(command))
        lines.append("")
    lines.append(
        f"{archive_target}: {_CANONICAL_REPLAY_FORCE_TARGET} "
        f"{' '.join(archive_dependencies)}")
    lines.append(recipe(archive_command))
    lines.append(recipe(index_command))
    lines.append("")
    lines.append(
        f"{target}: {_CANONICAL_REPLAY_FORCE_TARGET} "
        f"{benchmark_output} {archive_target}")
    lines.append(recipe(link_command))
    lines.append("")
    freeze_count = 0
    if artifact_transports is not None:
        header_relative = (
            "generated/leopard2-benchmark-attestation/"
            "leopard2_benchmark_source_attestation.h")
        lock_relative = header_relative + ".lock"
        expected_freeze = {
            *(output for output, _command, _role in compile_records),
            archive_target, target, header_relative, lock_relative,
        }
        require(set(artifact_transports) == expected_freeze,
                "canonical replay artifact transport closure differs")
        archive_compile_outputs = sorted(
            output for output, _command, role in compile_records
            if role == "archive")
        benchmark_compile_outputs = sorted(
            output for output, _command, role in compile_records
            if role == "benchmark")
        freeze_order = [
            header_relative, lock_relative,
            *archive_compile_outputs, *benchmark_compile_outputs,
            archive_target, target,
        ]
        require(len(freeze_order) == len(set(freeze_order)) and
                set(freeze_order) == expected_freeze,
                "canonical replay artifact freeze order differs")
        lines.append(f"{_CANONICAL_REPLAY_FREEZE_TARGET}: {target}")
        lines.append(recipe([transports["inventory"]]))
        for relative in freeze_order:
            copy_command = [transports["copy"], "--", relative]
            _require_shell_literal_tokens(
                copy_command, "canonical replay artifact freeze command")
            destination = artifact_transports[relative]
            destination_descriptor = int(destination.rsplit("/", 1)[1])
            rendered = shlex.join(copy_command)
            require("$" not in rendered and "\n" not in rendered and
                    "\r" not in rendered,
                    "canonical replay artifact freeze rendering is unsafe")
            # The shell duplicates an already-open memfd onto cat's stdout.
            # Landlock therefore never needs to authorize reopening an
            # anonymous inode through procfs.
            lines.append(
                "\t" + rendered + f" 1>&{destination_descriptor}")
        lines.append("")
        freeze_count = len(freeze_order)
    encoded = ("\n".join(lines)).encode("utf-8", errors="strict")
    require(encoded.endswith(b"\n") and
            len(encoded) <= MAX_METADATA_BYTES,
            "canonical replay plan exceeds its byte/framing contract")
    return encoded, {
        "object_count": len(compile_records),
        "archive_object_count": len(archive_dependencies),
        "command_count": len(compile_records) + 4 + freeze_count +
            (1 if artifact_transports is not None else 0),
        "target_count": len(compile_records) + 5 +
            (1 if artifact_transports is not None else 0),
        "recursive_make_count": 0,
        "freeze_output_count": freeze_count,
    }


def _canonical_replay_inventory_sha256(
    records: Sequence[Mapping[str, Any]],
) -> str:
    encoded = json.dumps(
        list(records), sort_keys=True, separators=(",", ":"),
        ensure_ascii=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _canonical_replay_plan_manifest(
    candidate: Mapping[str, Any],
    plan_path: Path,
    plan_content: bytes,
    plan_snapshot: _RetainedFileSnapshot,
    counts: Mapping[str, int],
    artifact_transports: Mapping[str, str] | None = None,
    execution_root_descriptor: int | None = None,
) -> dict[str, Any]:
    """Publish the exact independently-recomputable replay input inventory."""
    require(plan_path.name == Path(_CANONICAL_REPLAY_PLAN_RELATIVE).name and
            plan_snapshot.content == plan_content,
            "canonical replay plan snapshot differs")
    descriptor = plan_snapshot.executable_descriptor
    execution = plan_snapshot.executable_identity
    plan_record = {
        "role": "generated-plan",
        "root": "build",
        "path": _CANONICAL_REPLAY_PLAN_RELATIVE,
        "sha256": hashlib.sha256(plan_content).hexdigest(),
        "size": len(plan_content),
        "mode": stat.S_IMODE(plan_snapshot.identity["mode"]),
    }
    source = _canonical_replay_source_record(
        candidate, _CANONICAL_REPLAY_ATTESTATION_SOURCE)
    source_record = {
        "role": "tracked-cmake",
        "root": "source",
        "path": source["path"],
        "sha256": source["sha256"],
        "size": source["size"],
        "mode": source["mode"],
    }
    retained = sorted(
        (plan_record, source_record),
        key=lambda record: (record["root"], record["path"], record["role"]))
    count_keys = {
        "object_count", "archive_object_count", "command_count",
        "target_count", "recursive_make_count", "freeze_output_count",
    }
    require(isinstance(counts, Mapping) and set(counts) == count_keys and
            all(type(value) is int and value >= 0
                for value in counts.values()) and
            counts["recursive_make_count"] == 0,
            "canonical replay plan counts are malformed")
    output_sinks = [] if artifact_transports is None else [
        {
            "path": relative,
            "descriptor": int(destination.rsplit("/", 1)[1]),
        }
        for relative, destination in sorted(artifact_transports.items())
    ]
    require(counts["freeze_output_count"] == len(output_sinks),
            "canonical replay output-sink count differs")
    require(
        (execution_root_descriptor is None and not output_sinks) or
        (type(execution_root_descriptor) is int and
         execution_root_descriptor >= 0 and bool(output_sinks)),
        "canonical replay execution-root descriptor differs")
    return {
        "schema": CANONICAL_REPLAY_RECIPE_SCHEMA,
        "selected_target":
            f"CMakeFiles/{candidate['executable_target']}.dir/replay",
        "file_count": len(retained),
        "total_bytes": sum(record["size"] for record in retained),
        "counts": {key: counts[key] for key in sorted(counts)},
        "retained_inputs": retained,
        "inventory_sha256":
            _canonical_replay_inventory_sha256(retained),
        "output_sinks": output_sinks,
        "execution_root_descriptor": execution_root_descriptor,
        "plan_execution": {
            "sealed_descriptor": descriptor,
            "sealed_sha256": execution["sha256"],
            "sealed_size": execution["size"],
            "sealed_seals": execution["seals"],
        },
    }


class _TemporaryReplayRecipeTransport:
    """Route generated Make recipes through immutable sealed tool FDs."""

    def __init__(
        self, build: Path,
        replacements: Mapping[str, tuple[str, str]],
        required_labels: set[str],
        *,
        expected_compile_closure: Sequence[Mapping[str, Any]] | None = None,
        expected_source: Path | None = None,
    ) -> None:
        self.build = build.resolve(strict=True)
        semantic_enabled = expected_compile_closure is not None
        require(
            semantic_enabled == (expected_source is not None),
            "private replay semantic closure and source must be supplied "
            "together")
        semantic_source = (
            expected_source.resolve(strict=True)
            if expected_source is not None else None)
        self.originals: list[tuple[Path, bytes, int]] = []
        self.manifest: dict[str, Any] = {}
        self._guards: list[_InotifyMutationGuard] = []
        self._expected: dict[
            Path, tuple[dict[str, Any], bytes]] = {}
        counts = Counter()
        rewrites = []
        try:
            candidates = _replay_recipe_candidates(self.build)
            relative_candidates = {
                path.relative_to(self.build).as_posix()
                for path in candidates
            }
            require(
                _REPLAY_RECIPE_REQUIRED_PATHS.issubset(
                    relative_candidates),
                "private replay recipe inventory omits a required CMake "
                "input")
            capture_guard = _InotifyMutationGuard(
                "private replay recipe capture")
            try:
                for path in candidates:
                    capture_guard.add_file_path(path)
                commands_path = self.build / "compile_commands.json"
                if semantic_enabled:
                    capture_guard.add_file_path(commands_path)
                captured = {}
                captured_total = 0
                for path in candidates:
                    snapshot = file_snapshot(
                        path, "private replay generated recipe",
                        maximum_bytes=MAX_METADATA_BYTES)
                    captured[path] = snapshot
                    captured_total += len(snapshot[1])
                    require(
                        captured_total <= MAX_REPLAY_RECIPE_TOTAL_BYTES,
                        "private replay recipes exceed their aggregate byte "
                        "bound")
                commands_snapshot = None
                if semantic_enabled:
                    commands_snapshot = file_snapshot(
                        commands_path, "private replay compile commands",
                        maximum_bytes=MAX_METADATA_BYTES)
                capture_guard.verify()
                require(_replay_recipe_candidates(self.build) == candidates,
                        "private replay recipe inventory changed while "
                        "captured")
            finally:
                capture_guard.close()

            cmake_replacement = replacements.get("cmake")
            require(
                isinstance(cmake_replacement, tuple) and
                len(cmake_replacement) == 2 and
                isinstance(cmake_replacement[0], str),
                "private replay CMake transport is malformed")
            for path, (_identity, original) in captured.items():
                if path.name in ("Makefile", "Makefile2") or \
                        path.suffix == ".make":
                    _validate_replay_make_include_closure(
                        path, self.build, original, relative_candidates,
                        cmake_replacement[0],
                        None if semantic_source is None
                        else str(semantic_source))
            semantic_manifest = None
            if semantic_enabled:
                require(commands_snapshot is not None and
                        semantic_source is not None and
                        expected_compile_closure is not None,
                        "private replay semantic capture is incomplete")
                semantic_manifest = \
                    _validate_replay_compile_variable_closure(
                        self.build, semantic_source, commands_snapshot[1],
                        {path: content
                         for path, (_identity, content) in captured.items()},
                        expected_compile_closure, cmake_replacement[0])

            prepared: dict[
                Path, tuple[dict[str, Any], bytes, bytes]] = {}
            prepared_total = 0
            for path in candidates:
                identity, original = captured[path]
                if "make" in replacements:
                    # Command-line MAKE is inherited by recursive GNU Make.
                    # Count only actual recipe expansions, not CMake's
                    # explanatory comments mentioning $(MAKE).
                    counts["make"] += sum(
                        line.count(b"$(MAKE)")
                        for line in original.splitlines()
                        if line.startswith(b"\t"))
                rewritten = original
                path_counts = Counter()
                for label, (old, new) in replacements.items():
                    _require_safe_unicode(
                        old, f"private replay {label} identity path")
                    _require_safe_unicode(
                        new, f"private replay {label} transport path")
                    if label == "dependency-includes":
                        require(
                            old ==
                            "CMakeFiles/**/{compiler_,}depend.make" and
                            re.fullmatch(r"/proc/self/fd/[0-9]+", new)
                            is not None,
                            "private replay dependency-include transport "
                            "is malformed")
                        relative_parent = path.relative_to(
                            self.build).parent.as_posix().encode("utf-8")
                        pattern = re.compile(
                            rb"(?m)^include[ \t]+" +
                            re.escape(relative_parent) +
                            rb"/(?P<name>(?:compiler_)?depend\.make)$")
                        matches = list(pattern.finditer(rewritten))
                        tokens = list(
                            _REPLAY_DERIVED_DEPENDENCY_TOKEN.finditer(
                                rewritten))
                        names = [
                            match.group("name") for match in matches
                        ]
                        require(
                            len(matches) == len(tokens) and
                            len(names) == len(set(names)) and
                            (not names or (
                                path.name == "build.make" and
                                set(names).issubset({
                                    b"depend.make",
                                    b"compiler_depend.make",
                                }))),
                            "private replay generated recipe contains an "
                            "indirect, duplicate, or non-canonical derived "
                            "dependency include")
                        transport = new.encode("ascii")
                        rewritten, count = pattern.subn(
                            lambda unused_match:
                                b"include " + transport,
                            rewritten)
                        require(
                            _REPLAY_DERIVED_DEPENDENCY_TOKEN.search(
                                rewritten) is None,
                            "private replay generated recipe retains a "
                            "derived dependency include")
                    else:
                        pattern = re.compile(
                            rb"(?<![A-Za-z0-9_./+:-])" +
                            re.escape(old.encode("utf-8")) +
                            rb"(?![A-Za-z0-9_./+:-])")
                        rewritten, count = pattern.subn(
                            new.encode("utf-8"), rewritten)
                    counts[label] += count
                    path_counts[label] += count
                require(
                    len(rewritten) <= MAX_METADATA_BYTES,
                    "transported private replay recipe exceeds its per-file "
                    "byte bound")
                prepared_total += len(rewritten)
                require(
                    prepared_total <= MAX_REPLAY_RECIPE_TOTAL_BYTES,
                    "transported private replay recipes exceed their "
                    "aggregate byte bound")
                prepared[path] = (identity, original, rewritten)
                if rewritten != original:
                    mode = stat.S_IMODE(identity["mode"])
                    rewrites.append({
                        "path": path.relative_to(self.build).as_posix(),
                        "mode": mode,
                        "original_sha256":
                            hashlib.sha256(original).hexdigest(),
                        "transport_sha256":
                            hashlib.sha256(rewritten).hexdigest(),
                        "size": len(rewritten),
                        "replacement_counts": {
                            label: path_counts[label]
                            for label in sorted(path_counts)
                            if path_counts[label] != 0
                        },
                    })
            require(all(counts[label] > 0 for label in required_labels),
                    "private replay generated recipes omit a required "
                    "immutable tool transport: " +
                    repr({label: counts[label]
                          for label in sorted(required_labels)}))

            rewritten_paths = {
                path for path, (_identity, original, rewritten)
                in prepared.items() if rewritten != original
            }
            unrewritten_guard = _InotifyMutationGuard(
                "unrewritten clean replay recipes")
            self._guards.append(unrewritten_guard)
            for path in candidates:
                if path not in rewritten_paths:
                    unrewritten_guard.add_file_path(path)
            if semantic_enabled:
                unrewritten_guard.add_file_path(commands_path)
            unrewritten_guard.verify()

            for path in candidates:
                identity, original, rewritten = prepared[path]
                if path not in rewritten_paths:
                    continue
                mode = stat.S_IMODE(identity["mode"])
                self.originals.append((path, original, mode))
                _replace_private_file(path, rewritten, mode)

            rewritten_guard = _InotifyMutationGuard(
                "rewritten clean replay recipes")
            self._guards.append(rewritten_guard)
            for path in candidates:
                if path in rewritten_paths:
                    rewritten_guard.add_file_path(path)
            rewritten_guard.verify()

            retained_recipes = []
            for path in candidates:
                expected = prepared[path][2]
                identity, observed = file_snapshot(
                    path, "transported clean replay recipe",
                    maximum_bytes=MAX_METADATA_BYTES)
                require(observed == expected and
                        identity["path"] == str(path),
                        "transported clean replay recipe differs from its "
                        "captured bytes")
                self._expected[path] = (identity, expected)
                retained_recipes.append({
                    "path": path.relative_to(self.build).as_posix(),
                    "sha256": identity["sha256"],
                    "size": identity["size"],
                    "mode": identity["mode"],
                    "rewritten": path in rewritten_paths,
                })
            semantic_inputs = []
            if semantic_enabled:
                require(commands_snapshot is not None,
                        "private replay semantic snapshot is missing")
                commands_identity, commands_content = file_snapshot(
                    commands_path, "retained private replay compile commands",
                    maximum_bytes=MAX_METADATA_BYTES)
                require(commands_identity == commands_snapshot[0] and
                        commands_content == commands_snapshot[1],
                        "private replay compile commands changed before "
                        "execution")
                self._expected[commands_path] = (
                    commands_identity, commands_content)
                semantic_inputs.append({
                    "path": "compile_commands.json",
                    "sha256": commands_identity["sha256"],
                    "size": commands_identity["size"],
                    "mode": commands_identity["mode"],
                })
            require(_replay_recipe_candidates(self.build) == candidates,
                    "private replay recipe inventory changed before "
                    "execution")
            for guard in self._guards:
                guard.verify()
            retained_total = sum(
                record["size"] for record in retained_recipes)
            require(
                retained_total <= MAX_REPLAY_RECIPE_TOTAL_BYTES,
                "transported private replay recipes exceed their aggregate "
                "byte bound")
            self.manifest = {
                "schema": "leopard2-replay-recipe-transport/v2",
                "file_count": len(retained_recipes),
                "total_bytes": retained_total,
                "required_labels": sorted(required_labels),
                "replacement_counts": {
                    label: counts[label] for label in sorted(counts)
                },
                "rewrites": sorted(
                    rewrites, key=lambda record: record["path"]),
                "retained_recipes": retained_recipes,
                "inventory_sha256":
                    _replay_recipe_inventory_sha256(retained_recipes),
                "semantic_compile_closure": semantic_manifest,
                "semantic_inputs": semantic_inputs,
            }
        except BaseException as primary:
            try:
                self.close()
            except BaseException as cleanup_error:
                _raise_owner_exit_failure(
                    primary, cleanup_error,
                    "private replay recipe transport constructor",
                    "private replay recipe transport constructor cleanup "
                    f"failed: {cleanup_error}; primary failure: "
                    f"{type(primary).__name__}: {primary}")
            raise

    def _cleanup_owned_state(self) -> list[str]:
        """Release every retained resource without dropping retry ownership."""
        failures: list[str] = []
        failure: BaseException | None = None
        for guard in tuple(self._guards):
            try:
                guard.close()
            except BaseException as error:
                failures.append(
                    "cannot close transported recipe guard: " + str(error))
                failure = _owner_exception_precedence(
                    failure, error,
                    "transported recipe guard cleanup")
                try:
                    guard._close_without_verification()
                except BaseException as close_error:
                    failures.append(
                        "cannot force-close transported recipe guard: " +
                        str(close_error))
                    failure = _owner_exception_precedence(
                        failure, close_error,
                        "transported recipe forced guard cleanup")
            if guard.descriptor < 0 and guard in self._guards:
                self._guards.remove(guard)

        for record in reversed(tuple(self.originals)):
            path, content, mode = record
            try:
                _replace_private_file(path, content, mode)
                identity, observed = file_snapshot(
                    path, "restored private replay generated recipe",
                    maximum_bytes=MAX_METADATA_BYTES)
                require(observed == content and
                        stat.S_IMODE(identity["mode"]) == mode,
                        "private replay generated recipe was not restored")
            except BaseException as error:
                failures.append(f"{path}: {error}")
                failure = _owner_exception_precedence(
                    failure, error,
                    "transported recipe restoration cleanup")
            else:
                if record in self.originals:
                    self.originals.remove(record)
        if not self.originals and not self._guards:
            self._expected = {}
        if isinstance(failure, _OWNER_TERMINAL_EXCEPTIONS):
            raise failure.with_traceback(failure.__traceback__)
        return failures

    def close(self) -> None:
        failures: list[str] = []
        verification_error: BaseException | None = None
        cleanup_error: BaseException | None = None
        try:
            for guard in self._guards:
                try:
                    guard.verify()
                except BaseException as error:
                    failures.append(
                        "transported recipe changed while retained: " +
                        str(error))
                    verification_error = _owner_exception_precedence(
                        verification_error, error,
                        "transported recipe guard verification")
            for path, (expected_identity, expected) in self._expected.items():
                try:
                    identity, observed = file_snapshot(
                        path, "retained transported clean replay recipe",
                        maximum_bytes=MAX_METADATA_BYTES)
                    require(identity == expected_identity and
                            observed == expected,
                            "transported recipe changed while retained")
                except BaseException as error:
                    failures.append(f"{path}: {error}")
                    verification_error = _owner_exception_precedence(
                        verification_error, error,
                        "transported recipe file verification")
        except BaseException as error:
            verification_error = _owner_exception_precedence(
                verification_error, error,
                "transported recipe verification")
        finally:
            # Cleanup is a separate retryable state machine.  In particular,
            # an asynchronous BaseException anywhere in verification must not
            # strand rewritten recipes or mutation-guard descriptors.
            try:
                observed_failures = self._cleanup_owned_state()
                failures.extend(observed_failures)
                if observed_failures:
                    cleanup_error = _owner_exception_precedence(
                        cleanup_error,
                        BuildProvenanceError("; ".join(observed_failures)),
                        "transported recipe cleanup diagnostics")
            except BaseException as error:
                cleanup_error = _owner_exception_precedence(
                    cleanup_error, error,
                    "transported recipe cleanup")
            # _cleanup_owned_state records per-resource BaseExceptions so it
            # can continue releasing unrelated resources.  Such a caught
            # interruption therefore does not reach the exception handler
            # above.  Retry the bounded remaining ownership set once whether
            # the first pass returned or raised, while retaining its original
            # diagnostic.
            if self.originals or self._guards:
                try:
                    observed_failures = self._cleanup_owned_state()
                    failures.extend(observed_failures)
                    if observed_failures:
                        cleanup_error = _owner_exception_precedence(
                            cleanup_error,
                            BuildProvenanceError(
                                "; ".join(observed_failures)),
                            "transported recipe cleanup retry diagnostics")
                except BaseException as retry_error:
                    failures.append(
                        "interrupted replay transport cleanup retry failed: "
                        f"{type(retry_error).__name__}: {retry_error}")
                    cleanup_error = _owner_exception_precedence(
                        cleanup_error, retry_error,
                        "transported recipe cleanup retry")

        if self.originals or self._guards:
            failures.append(
                "replay transport still owns rewritten recipes or guards")
        authoritative = _owner_exception_precedence(
            verification_error, cleanup_error,
            "private replay recipe transport close")
        if isinstance(authoritative, _OWNER_TERMINAL_EXCEPTIONS):
            raise authoritative.with_traceback(authoritative.__traceback__)
        if failures:
            contexts = []
            for error in (verification_error, cleanup_error):
                if error is not None:
                    contexts.append(
                        f"{type(error).__name__}: {error}")
            if contexts:
                failures.append("interrupted phase: " + "; ".join(contexts))
        elif verification_error is not None:
            raise verification_error
        elif cleanup_error is not None:
            raise cleanup_error
        require(not failures,
                "cannot verify/restore private replay generated recipes: " +
                "; ".join(failures))

    def __enter__(self) -> "_TemporaryReplayRecipeTransport":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, tb
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc is None:
                raise
            _raise_owner_exit_failure(
                exc, cleanup_error, "private replay recipe transport owner",
                "private replay recipe transport cleanup failed: "
                f"{cleanup_error}; primary failure: "
                f"{type(exc).__name__}: {exc}")


def _validate_canonical_replay_plan_manifest(
    raw: Any,
    candidate: Mapping[str, Any],
    source: Path,
    build: Path,
    transports: Mapping[str, str],
    *,
    label: str,
) -> int:
    """Bind a canonical replay plan to candidate semantics offline."""
    top_keys = {
        "schema", "selected_target", "file_count", "total_bytes", "counts",
        "retained_inputs", "inventory_sha256", "output_sinks",
        "execution_root_descriptor", "plan_execution",
    }
    require(isinstance(raw, Mapping) and set(raw) == top_keys and
            raw.get("schema") == CANONICAL_REPLAY_RECIPE_SCHEMA,
            f"{label} canonical replay plan manifest is malformed")
    execution = raw["plan_execution"]
    execution_keys = {
        "sealed_descriptor", "sealed_sha256", "sealed_size", "sealed_seals",
    }
    require(isinstance(execution, Mapping) and
            set(execution) == execution_keys and
            type(execution["sealed_descriptor"]) is int and
            execution["sealed_descriptor"] >= 0 and
            isinstance(execution["sealed_sha256"], str) and
            re.fullmatch(
                r"[0-9a-f]{64}", execution["sealed_sha256"]) is not None and
            type(execution["sealed_size"]) is int and
            execution["sealed_size"] > 0 and
            type(execution["sealed_seals"]) is int,
            f"{label} canonical replay sealed plan identity is malformed")
    required_seals = (
        getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
        getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
        getattr(fcntl, "F_SEAL_GROW", 0x0004) |
        getattr(fcntl, "F_SEAL_WRITE", 0x0008))
    require(execution["sealed_seals"] & required_seals == required_seals,
            f"{label} canonical replay plan is not fully sealed")

    raw_sinks = raw["output_sinks"]
    require(
        isinstance(raw_sinks, list) and
        all(
            isinstance(item, Mapping) and
            set(item) == {"path", "descriptor"} and
            isinstance(item["path"], str) and item["path"] and
            type(item["descriptor"]) is int and item["descriptor"] >= 0
            for item in raw_sinks) and
        raw_sinks == sorted(raw_sinks, key=lambda item: item["path"]) and
        len({item["path"] for item in raw_sinks}) == len(raw_sinks) and
        len({item["descriptor"] for item in raw_sinks}) == len(raw_sinks),
        f"{label} canonical replay output sinks are malformed")
    artifact_transports = ({
        item["path"]: f"/proc/self/fd/{item['descriptor']}"
        for item in raw_sinks
    } if raw_sinks else None)
    execution_root_descriptor = raw["execution_root_descriptor"]
    require(
        (not raw_sinks and execution_root_descriptor is None) or
        (raw_sinks and type(execution_root_descriptor) is int and
         execution_root_descriptor >= 0),
        f"{label} canonical replay execution root is malformed")
    execution_root = (
        f"/proc/self/fd/{execution_root_descriptor}"
        if execution_root_descriptor is not None else None)
    content, counts = _canonical_replay_makefile_bytes(
        candidate, source, build, transports, artifact_transports,
        execution_build_root=execution_root)
    content_digest = hashlib.sha256(content).hexdigest()
    source_record = _canonical_replay_source_record(
        candidate, _CANONICAL_REPLAY_ATTESTATION_SOURCE)
    retained = sorted((
        {
            "role": "generated-plan",
            "root": "build",
            "path": _CANONICAL_REPLAY_PLAN_RELATIVE,
            "sha256": content_digest,
            "size": len(content),
            "mode": 0o700,
        },
        {
            "role": "tracked-cmake",
            "root": "source",
            "path": source_record["path"],
            "sha256": source_record["sha256"],
            "size": source_record["size"],
            "mode": source_record["mode"],
        },
    ), key=lambda record: (
        record["root"], record["path"], record["role"]))
    expected = {
        "schema": CANONICAL_REPLAY_RECIPE_SCHEMA,
        "selected_target":
            f"CMakeFiles/{candidate.get('executable_target')}.dir/replay",
        "file_count": len(retained),
        "total_bytes": sum(record["size"] for record in retained),
        "counts": {key: counts[key] for key in sorted(counts)},
        "retained_inputs": retained,
        "inventory_sha256":
            _canonical_replay_inventory_sha256(retained),
        "output_sinks": raw_sinks,
        "execution_root_descriptor": execution_root_descriptor,
        "plan_execution": {
            "sealed_descriptor": execution["sealed_descriptor"],
            "sealed_sha256": content_digest,
            "sealed_size": len(content),
            "sealed_seals": execution["sealed_seals"],
        },
    }
    require(raw == expected,
            f"{label} canonical replay retained input inventory differs")
    return int(execution["sealed_descriptor"])


def _validate_legacy_replay_recipe_transport(
    raw_recipe: Any,
    candidate: Mapping[str, Any],
    replacement_by_label: Mapping[str, Mapping[str, Any]],
    *,
    label: str,
) -> None:
    """Preserve exact validation of historical generated-recipe v2 proofs."""

    def exact(value: Any, keys: set[str], what: str) -> Mapping[str, Any]:
        require(isinstance(value, Mapping) and set(value) == keys,
                f"{label} {what} is malformed")
        return value

    def digest(value: Any, what: str) -> str:
        require(isinstance(value, str) and
                re.fullmatch(r"[0-9a-f]{64}", value) is not None,
                f"{label} {what} digest is malformed")
        return value

    recipe_keys = {
        "schema", "file_count", "total_bytes", "required_labels",
        "replacement_counts", "rewrites", "retained_recipes",
        "inventory_sha256", "semantic_compile_closure",
        "semantic_inputs",
    }
    recipe = exact(
        raw_recipe, recipe_keys, "replay recipe transport")
    require(recipe["schema"] == "leopard2-replay-recipe-transport/v2" and
            recipe["required_labels"] ==
                sorted({"cxx-compiler", "archiver", "ranlib",
                        "cmake", "make", "git",
                        "dependency-includes"}) and
            isinstance(recipe["replacement_counts"], Mapping) and
            bool(recipe["replacement_counts"]) and
            set(recipe["replacement_counts"]) ==
                set(replacement_by_label) and
            all(type(count) is int and count >= 0
                for count in recipe["replacement_counts"].values()),
            f"{label} replay recipe replacement counts are malformed")
    semantic_keys = {
        "schema", "compile_commands_sha256", "object_count",
        "target_count", "control_assignment_counts",
        "candidate_compile_closure_sha256",
        "candidate_target_variables_sha256",
    }
    semantic = exact(
        recipe["semantic_compile_closure"], semantic_keys,
        "replay semantic compile closure")
    semantic_inputs = recipe["semantic_inputs"]
    require(isinstance(semantic_inputs, list) and
            len(semantic_inputs) == 1,
            f"{label} replay semantic input set is malformed")
    semantic_input = exact(
        semantic_inputs[0], {"path", "sha256", "size", "mode"},
        "replay semantic input")
    candidate_compile_closure = candidate.get(
        "source_object_compile_closure")
    require(isinstance(candidate_compile_closure, list) and
            candidate_compile_closure,
            f"{label} closure lacks its production compile objects")
    candidate_targets = set()
    for item in candidate_compile_closure:
        require(isinstance(item, Mapping) and
                isinstance(item.get("compile_entry"), Mapping) and
                isinstance(
                    item["compile_entry"].get("normalized_arguments"), list),
                f"{label} production compile closure is malformed")
        arguments = item["compile_entry"]["normalized_arguments"]
        require(len(arguments) >= 5 and arguments[-4] == "-o" and
                isinstance(arguments[-3], str),
                f"{label} production compile argv is malformed")
        output_parts = Path(arguments[-3]).parts
        require(len(output_parts) >= 3 and
                output_parts[0] == "CMakeFiles" and
                output_parts[1].endswith(".dir"),
                f"{label} production compile output is malformed")
        candidate_targets.add("/".join(output_parts[:2]))
    control_counts = semantic["control_assignment_counts"]
    _require_replay_candidate_compile_closure_binding(
        semantic, candidate_compile_closure, label=label)
    require(
        semantic["schema"] ==
            "leopard2-replay-compile-variable-closure/v2" and
        semantic_input["path"] == "compile_commands.json" and
        semantic_input["sha256"] ==
            semantic["compile_commands_sha256"] and
        digest(semantic_input["sha256"],
               "replay compile commands") ==
            semantic_input["sha256"] and
        type(semantic_input["size"]) is int and
        0 < semantic_input["size"] <= MAX_METADATA_BYTES and
        type(semantic_input["mode"]) is int and
        stat.S_ISREG(semantic_input["mode"]) and
        semantic["object_count"] == len(candidate_compile_closure) and
        semantic["target_count"] == len(candidate_targets) and
        isinstance(control_counts, Mapping) and
        set(control_counts) == {
            "CMAKE_COMMAND", "CMAKE_BINARY_DIR", "CMAKE_SOURCE_DIR"} and
        all(type(count) is int and count > 0
            for count in control_counts.values()),
        f"{label} replay semantic compile closure differs")
    require(all(recipe["replacement_counts"][item] > 0
                for item in recipe["required_labels"]),
            f"{label} replay recipe omits a required immutable transport")
    rewrite_keys = {
        "path", "mode", "original_sha256", "transport_sha256", "size",
        "replacement_counts",
    }
    rewrites = recipe["rewrites"]
    retained = recipe["retained_recipes"]
    require(isinstance(rewrites, list) and rewrites and
            isinstance(retained, list) and len(retained) >= len(rewrites) and
            type(recipe["file_count"]) is int and
            recipe["file_count"] == len(retained) and
            recipe["file_count"] <= MAX_REPLAY_RECIPE_FILES and
            type(recipe["total_bytes"]) is int and
            0 <= recipe["total_bytes"] <=
                MAX_REPLAY_RECIPE_TOTAL_BYTES,
            f"{label} replay recipe retained set is malformed")
    rewrite_by_path: dict[str, Mapping[str, Any]] = {}
    sums = Counter()
    for raw in rewrites:
        rewrite = exact(raw, rewrite_keys, "replay recipe rewrite")
        path = rewrite["path"]
        require(_is_replay_recipe_relative_path(path) and
                path not in rewrite_by_path and
                type(rewrite["mode"]) is int and
                0 <= rewrite["mode"] <= 0o7777 and
                type(rewrite["size"]) is int and
                0 <= rewrite["size"] <= MAX_METADATA_BYTES and
                isinstance(rewrite["replacement_counts"], Mapping) and
                bool(rewrite["replacement_counts"]) and
                set(rewrite["replacement_counts"]).issubset(
                    replacement_by_label) and
                all(type(count) is int and count > 0
                    for count in rewrite["replacement_counts"].values()),
                f"{label} replay recipe rewrite is malformed")
        digest(rewrite["original_sha256"], "recipe original")
        digest(rewrite["transport_sha256"], "recipe transport")
        require(rewrite["original_sha256"] !=
                rewrite["transport_sha256"],
                f"{label} replay recipe rewrite does not change bytes")
        sums.update(rewrite["replacement_counts"])
        rewrite_by_path[path] = rewrite
    require(rewrites == sorted(rewrites, key=lambda item: item["path"]),
            f"{label} replay recipe rewrites are not canonical")
    for replacement_label, count in recipe["replacement_counts"].items():
        if replacement_label == "make":
            require(count >= sums[replacement_label],
                    f"{label} recursive Make transport count differs")
        else:
            require(count == sums[replacement_label],
                    f"{label} recipe {replacement_label} count differs")
    retained_paths = set()
    retained_keys = {"path", "sha256", "size", "mode", "rewritten"}
    for raw in retained:
        record = exact(raw, retained_keys, "retained replay recipe")
        path = record["path"]
        require(_is_replay_recipe_relative_path(path),
                f"{label} retained replay recipe path is malformed")
        rewrite = rewrite_by_path.get(path)
        require(path not in retained_paths and
                type(record["rewritten"]) is bool and
                record["rewritten"] is (rewrite is not None) and
                type(record["size"]) is int and
                0 <= record["size"] <= MAX_METADATA_BYTES and
                type(record["mode"]) is int and record["mode"] >= 0 and
                stat.S_ISREG(record["mode"]),
                f"{label} retained replay recipe identity differs")
        digest(record["sha256"], "retained recipe")
        if rewrite is not None:
            require(
                record["sha256"] == rewrite["transport_sha256"] and
                record["size"] == rewrite["size"] and
                stat.S_IMODE(record["mode"]) == rewrite["mode"],
                f"{label} rewritten replay recipe identity differs")
        retained_paths.add(path)
    require(sum(record["size"] for record in retained) ==
                recipe["total_bytes"] and
            set(rewrite_by_path).issubset(retained_paths) and
            _REPLAY_RECIPE_REQUIRED_PATHS.issubset(retained_paths) and
            retained == sorted(retained, key=lambda item: item["path"]) and
            digest(recipe["inventory_sha256"],
                   "retained recipe inventory") ==
                _replay_recipe_inventory_sha256(retained),
            f"{label} retained replay recipe path set differs")


def _capture_replayed_candidate_provenance(
    build_root: Path,
    source_root: Path,
    target: str,
    *,
    replay_contract: Mapping[str, str],
    inherited_descriptors: Sequence[int],
    tracked_source_manifest: Mapping[str, Any],
    logical_source_root: Path,
    sealed_artifacts: Mapping[Path | str, _ReplayArtifactSink],
) -> dict[str, Any]:
    """Capture a rebuilt candidate under its authenticated schema."""
    configuration_schema = replay_contract.get("configuration_schema")
    require(configuration_schema in {
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            },
            "replayed candidate capture configuration schema is unsupported")
    return candidate_build_provenance(
        build_root, source_root, build_root / target, target,
        inherited_descriptors=inherited_descriptors,
        tracked_source_manifest=tracked_source_manifest,
        logical_source_root=logical_source_root,
        _expected_configuration_schema=configuration_schema,
        sealed_artifacts=sealed_artifacts)


def verify_reproducible_candidate_build(
    candidate: Mapping[str, Any], *, jobs: int | None = None,
    inherited_descriptors: Sequence[int] = (),
) -> dict[str, Any]:
    """Reconfigure/build in an empty directory and compare all linked bytes."""
    replay_contract = _reproducible_replay_contract(candidate)
    require(
        replay_contract["closure_schema"] ==
            PRODUCTION_BUILD_CLOSURE_SCHEMA,
        "clean replay generation requires a current production-build closure")
    replay_configuration_schema = replay_contract["configuration_schema"]
    require(replay_configuration_schema in {
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V5,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V6,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V7,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V8,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA_V9,
                BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            },
            "clean replay generation supports only v5 through v9 and current "
            "configuration contracts")
    source = Path(str(candidate.get("source_root", ""))).resolve(strict=True)
    target = str(candidate.get("executable_target", ""))
    require(re.fullmatch(r"[A-Za-z0-9_.+][A-Za-z0-9_.+-]*", target)
            is not None,
            "candidate executable target is unsafe")
    cache = candidate.get("validated_cache")
    require(isinstance(cache, dict),
            "candidate build has no retained CMake semantics")
    cmake = Path("/usr/bin/cmake").resolve(strict=True)
    shell = Path("/bin/sh").resolve(strict=True)
    bash = Path("/bin/bash").resolve(strict=True)
    copy_tool = Path("/usr/bin/cat").resolve(strict=True)
    retained_manifest = candidate.get("tracked_source_manifest")
    require(isinstance(retained_manifest, dict),
            "candidate build has no retained tracked source tree")
    if jobs is None:
        jobs = min(128, len(os.sched_getaffinity(0)))
    require(type(jobs) is int and 1 <= jobs <= 128,
            "clean rebuild job count must be in 1..128")

    workspace = Path(tempfile.mkdtemp(prefix="leopard2-proof-workspace-"))
    temporary = workspace / "build"
    private_source = workspace / "source"
    private_tools = workspace / "tools"
    private_tmp = workspace / "tmp"
    private_stage = workspace / "stage"
    try:
        private_tmp.mkdir(mode=0o700)
        private_stage.mkdir(mode=0o700)
        with _RetainedPrivateSourceTree(
                source, private_source,
                inherited_descriptors=inherited_descriptors) as source_tree, \
                ExitStack() as retained_tools:
            require(source_tree.manifest == retained_manifest,
                    "candidate tracked source tree changed before clean "
                    "rebuild")
            private_tools.mkdir(mode=0o700)
            tool_records = (
                ("CMAKE_C_COMPILER", "c_compiler", "candidate C compiler"),
                ("CMAKE_CXX_COMPILER", "compiler",
                 "candidate C++ compiler"),
                ("CMAKE_AR", "archiver", "candidate archiver"),
                ("CMAKE_RANLIB", "ranlib", "candidate ranlib"),
                ("CMAKE_MAKE_PROGRAM", "make_program",
                 "candidate make program"),
                ("CMAKE_LINKER", "linker", "candidate linker"),
                ("LEO2_BENCHMARK_GIT_EXECUTABLE", "benchmark_git",
                 "candidate benchmark Git executable"),
            )
            tool_snapshots: dict[str, _RetainedFileSnapshot] = {}
            for cache_name, record_name, label in tool_records:
                path = cache.get(cache_name)
                record = candidate.get(record_name)
                require(isinstance(path, str) and path and
                        isinstance(record, dict),
                        f"candidate build lacks retained {label}")
                snapshot = retained_tools.enter_context(
                    _RetainedFileSnapshot(path, label))
                require(snapshot.identity == record,
                        f"{label} changed before clean rebuild")
                snapshot.executable_descriptor
                tool_snapshots[record_name] = snapshot
            subtools = candidate.get("compiler_subtools")
            require(isinstance(subtools, list) and subtools,
                    "candidate build lacks retained compiler subtools")
            subtool_snapshots: list[
                tuple[Mapping[str, Any], _RetainedFileSnapshot]
            ] = []
            for index, record in enumerate(subtools):
                require(isinstance(record, dict) and
                        record.get("language") in ("c", "c++") and
                        record.get("role") in (
                            "cc1", "cc1plus", "as", "collect2", "ld") and
                        isinstance(record.get("identity"), dict) and
                        isinstance(record["identity"].get("path"), str),
                        "candidate compiler subtool record is malformed")
                snapshot = retained_tools.enter_context(
                    _RetainedFileSnapshot(
                        record["identity"]["path"],
                        f"candidate compiler subtool {index}"))
                require(snapshot.identity == record["identity"],
                        "candidate compiler subtool changed before clean "
                        "rebuild")
                snapshot.executable_descriptor
                subtool_snapshots.append((record, snapshot))
            cmake_snapshot = retained_tools.enter_context(
                _RetainedFileSnapshot(cmake, "clean replay CMake"))
            shell_snapshot = retained_tools.enter_context(
                _RetainedFileSnapshot(shell, "clean replay command shell"))
            bash_snapshot = retained_tools.enter_context(
                _RetainedFileSnapshot(bash, "clean replay Bash"))
            copy_snapshot = retained_tools.enter_context(
                _RetainedFileSnapshot(
                    copy_tool, "clean replay artifact stream tool"))
            for snapshot in (
                    cmake_snapshot, shell_snapshot, bash_snapshot,
                    copy_snapshot):
                snapshot.executable_descriptor

            configure = _reproducible_configure_argv(
                private_source, temporary, cache)
            configure_environment_overrides = {
                "CC": str(cache["CMAKE_C_COMPILER"]),
                "CXX": str(cache["CMAKE_CXX_COMPILER"]),
                "TMPDIR": str(private_tmp),
            }
            configure_inherited_descriptors = (
                *inherited_descriptors,
                cmake_snapshot.executable_descriptor,
            )
            configure_invocation = _replay_invocation_record(
                configure, executable_label="cmake",
                executable=cmake_snapshot,
                inherited_descriptors=configure_inherited_descriptors,
                environment_overrides=configure_environment_overrides,
                maximum_bytes=16 << 20, timeout=300,
                write_sandbox_root=workspace)
            _run(
                configure, "runner-owned candidate configure",
                maximum_bytes=16 << 20, timeout=300,
                inherited_descriptors=configure_inherited_descriptors,
                executable_descriptor=cmake_snapshot.executable_descriptor,
                environment_overrides=configure_environment_overrides,
                write_sandbox_root=workspace)
            for snapshot in (
                    *tool_snapshots.values(),
                    *(item[1] for item in subtool_snapshots),
                    cmake_snapshot, shell_snapshot, bash_snapshot,
                    copy_snapshot):
                snapshot.verify()

            replay_git = private_tools / "git"
            replay_git_bytes, replay_git_specification = \
                _write_replay_git_identity(
                    replay_git, private_source, source_tree.manifest,
                    interpreter_descriptor=
                        bash_snapshot.executable_descriptor)
            replay_git_snapshot = retained_tools.enter_context(
                _retain_exact_generated_executable(
                    replay_git, replay_git_bytes,
                    "clean replay Git identity executable"))
            empty_dependency_include = (
                private_tools / "empty-dependency-include")
            empty_dependency_bytes = REPLAY_EMPTY_DEPENDENCY_INCLUDE
            _write_private_executable(
                empty_dependency_include, empty_dependency_bytes)
            empty_dependency_snapshot = retained_tools.enter_context(
                _retain_exact_generated_executable(
                    empty_dependency_include, empty_dependency_bytes,
                    "clean replay empty dependency include"))

            prefix_guards: dict[str, _RetainedDirectoryTree] = {}
            prefix_records = []
            for language, directory_name in (
                    ("c", "c-prefix"), ("c++", "cxx-prefix")):
                prefix = private_tools / directory_name
                prefix.mkdir(mode=0o700)
                roles: dict[str, _RetainedFileSnapshot] = {}
                for record, snapshot in subtool_snapshots:
                    if record["language"] != language:
                        continue
                    role = str(record["role"])
                    previous = roles.get(role)
                    require(previous is None or
                            previous.identity == snapshot.identity,
                            f"candidate {language} compiler exposed "
                            f"conflicting {role} subtools")
                    roles.setdefault(role, snapshot)
                roles.setdefault("ld", tool_snapshots["linker"])
                mappings = []
                expected_targets = {}
                for role, snapshot in sorted(roles.items()):
                    target_path = (
                        f"/proc/self/fd/{snapshot.executable_descriptor}")
                    os.symlink(
                        target_path, prefix / role)
                    require(os.readlink(prefix / role) == target_path,
                            f"clean replay {language} compiler prefix "
                            f"mapping changed for {role}")
                    expected_targets[role] = target_path
                    mappings.append({
                        "role": role,
                        "target_descriptor":
                            snapshot.executable_descriptor,
                        "source_path": snapshot.identity["path"],
                        "source_sha256": snapshot.identity["sha256"],
                        "sealed_sha256":
                            snapshot.executable_identity["sha256"],
                        "sealed_size": snapshot.executable_identity["size"],
                        "sealed_seals": snapshot.executable_identity["seals"],
                    })
                guard = retained_tools.enter_context(
                    _retain_exact_symlink_directory(
                        prefix, expected_targets,
                        f"clean replay {language} compiler prefix"))
                prefix_guards[language] = guard
                prefix_records.append({
                    "language": language,
                    "path": str(prefix),
                    "descriptor": guard.descriptor,
                    "mode": os.fstat(guard.descriptor).st_mode,
                    "mappings": mappings,
                })

            wrapper_records = []
            def retain_wrapper(
                name: str, lexical: str,
                executable: _RetainedFileSnapshot,
                extra_arguments: Sequence[str] = (),
                underlying_label: str = "",
            ) -> _RetainedFileSnapshot:
                path = private_tools / name
                expected_script = _write_replay_exec_wrapper(
                    path, bash_snapshot.executable_descriptor, lexical,
                    executable.executable_descriptor, extra_arguments)
                wrapper = retained_tools.enter_context(
                    _retain_exact_generated_executable(
                        path, expected_script,
                        f"clean replay {name} wrapper"))
                wrapper_records.append({
                    "label": name,
                    "underlying_label": underlying_label,
                    "logical_argv0": lexical,
                    "injected_arguments": list(extra_arguments),
                    "interpreter_descriptor":
                        bash_snapshot.executable_descriptor,
                    "interpreter_sha256":
                        bash_snapshot.executable_identity["sha256"],
                    "underlying_descriptor":
                        executable.executable_descriptor,
                    "underlying_source_path": executable.identity["path"],
                    "underlying_source_sha256":
                        executable.identity["sha256"],
                    "underlying_sealed_sha256":
                        executable.executable_identity["sha256"],
                    "script_sha256":
                        hashlib.sha256(expected_script).hexdigest(),
                    "script_size": len(expected_script),
                    "sealed_descriptor": wrapper.executable_descriptor,
                    "sealed_sha256":
                        wrapper.executable_identity["sha256"],
                    "sealed_size": wrapper.executable_identity["size"],
                    "sealed_seals": wrapper.executable_identity["seals"],
                })
                return wrapper

            c_wrapper = retain_wrapper(
                "cc", str(cache["CMAKE_C_COMPILER"]),
                tool_snapshots["c_compiler"],
                (f"-B/proc/self/fd/"
                 f"{prefix_guards['c'].descriptor}/",),
                "c_compiler")
            cxx_wrapper = retain_wrapper(
                "cxx", str(cache["CMAKE_CXX_COMPILER"]),
                tool_snapshots["compiler"],
                (f"-B/proc/self/fd/"
                 f"{prefix_guards['c++'].descriptor}/",),
                "compiler")
            archiver_wrapper = retain_wrapper(
                "ar", str(cache["CMAKE_AR"]),
                tool_snapshots["archiver"],
                underlying_label="archiver")
            ranlib_wrapper = retain_wrapper(
                "ranlib", str(cache["CMAKE_RANLIB"]),
                tool_snapshots["ranlib"],
                underlying_label="ranlib")
            linker_wrapper = retain_wrapper(
                "ld", str(cache["CMAKE_LINKER"]),
                tool_snapshots["linker"],
                underlying_label="linker")
            cmake_wrapper = retain_wrapper(
                "cmake", str(cmake), cmake_snapshot,
                underlying_label="cmake")
            make_wrapper = retain_wrapper(
                "make", str(cache["CMAKE_MAKE_PROGRAM"]),
                tool_snapshots["make_program"],
                underlying_label="make_program")
            copy_wrapper = retain_wrapper(
                "copy", str(copy_tool), copy_snapshot,
                underlying_label="copy")

            replacements = {
                "c-compiler": (
                    str(cache["CMAKE_C_COMPILER"]),
                    f"/proc/self/fd/{c_wrapper.executable_descriptor}"),
                "cxx-compiler": (
                    str(cache["CMAKE_CXX_COMPILER"]),
                    f"/proc/self/fd/{cxx_wrapper.executable_descriptor}"),
                "archiver": (
                    str(cache["CMAKE_AR"]),
                    f"/proc/self/fd/{archiver_wrapper.executable_descriptor}"),
                "ranlib": (
                    str(cache["CMAKE_RANLIB"]),
                    f"/proc/self/fd/{ranlib_wrapper.executable_descriptor}"),
                "linker": (
                    str(cache["CMAKE_LINKER"]),
                    f"/proc/self/fd/{linker_wrapper.executable_descriptor}"),
                "cmake": (
                    str(cmake),
                    f"/proc/self/fd/{cmake_wrapper.executable_descriptor}"),
                # CMake's Unix Makefiles normally spell recursive invocations
                # as $(MAKE), but retain this exact lexical replacement for
                # generator versions that emit CMAKE_MAKE_PROGRAM directly.
                "make": (
                    str(cache["CMAKE_MAKE_PROGRAM"]),
                    f"/proc/self/fd/{make_wrapper.executable_descriptor}"),
                "git": (
                    str(cache["LEO2_BENCHMARK_GIT_EXECUTABLE"]),
                    f"/proc/self/fd/"
                    f"{replay_git_snapshot.executable_descriptor}"),
                "dependency-includes": (
                    "CMakeFiles/**/{compiler_,}depend.make",
                    f"/proc/self/fd/"
                    f"{empty_dependency_snapshot.executable_descriptor}"),
            }
            stage_outputs = _canonical_replay_output_paths(
                candidate, private_source, private_stage)
            stage_owner = _OwnedDescriptor()
            stage_descriptor = stage_owner.open(
                private_stage,
                getattr(os, "O_PATH", os.O_RDONLY) |
                getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_NOFOLLOW", 0))
            retained_tools.callback(stage_owner.close)
            execution_stage_root = f"/proc/self/fd/{stage_descriptor}"
            stage_artifacts = [
                path for path in stage_outputs
                if path.relative_to(private_stage).as_posix() !=
                    _CANONICAL_REPLAY_PLAN_RELATIVE
            ]
            require(stage_artifacts,
                    "canonical replay private artifact closure is empty")
            artifact_sinks: dict[str, _ReplayArtifactSink] = {}
            for stage_path in stage_artifacts:
                relative = stage_path.relative_to(private_stage).as_posix()
                final_path = temporary / relative
                artifact_sinks[relative] = retained_tools.enter_context(
                    _ReplayArtifactSink(
                        final_path, f"canonical replay artifact {relative}"))
            inventory_path = private_tools / "verify-output-inventory"
            inventory_bytes = _private_replay_inventory_bytes(
                bash_snapshot.executable_descriptor,
                sorted(artifact_sinks))
            _write_private_executable(inventory_path, inventory_bytes)
            inventory_snapshot = retained_tools.enter_context(
                _retain_exact_generated_executable(
                    inventory_path, inventory_bytes,
                    "clean replay private-output inventory verifier"))
            artifact_transports = {
                relative: sink.proc_path
                for relative, sink in artifact_sinks.items()
            }
            transport_descriptors = {
                *inherited_descriptors,
                *(snapshot.executable_descriptor
                  for snapshot in tool_snapshots.values()),
                *(snapshot.executable_descriptor
                  for _record, snapshot in subtool_snapshots),
                cmake_snapshot.executable_descriptor,
                shell_snapshot.executable_descriptor,
                bash_snapshot.executable_descriptor,
                copy_snapshot.executable_descriptor,
                replay_git_snapshot.executable_descriptor,
                empty_dependency_snapshot.executable_descriptor,
                inventory_snapshot.executable_descriptor,
                *(guard.descriptor for guard in prefix_guards.values()),
                c_wrapper.executable_descriptor,
                cxx_wrapper.executable_descriptor,
                archiver_wrapper.executable_descriptor,
                ranlib_wrapper.executable_descriptor,
                linker_wrapper.executable_descriptor,
                cmake_wrapper.executable_descriptor,
                make_wrapper.executable_descriptor,
                copy_wrapper.executable_descriptor,
                *(sink.descriptor for sink in artifact_sinks.values()),
                stage_descriptor,
            }
            canonical_transports = {
                "cxx":
                    f"/proc/self/fd/{cxx_wrapper.executable_descriptor}",
                "archiver":
                    f"/proc/self/fd/{archiver_wrapper.executable_descriptor}",
                "ranlib":
                    f"/proc/self/fd/{ranlib_wrapper.executable_descriptor}",
                "cmake":
                    f"/proc/self/fd/{cmake_wrapper.executable_descriptor}",
                "git":
                    f"/proc/self/fd/"
                    f"{replay_git_snapshot.executable_descriptor}",
                "copy":
                    f"/proc/self/fd/{copy_wrapper.executable_descriptor}",
                "inventory":
                    f"/proc/self/fd/{inventory_snapshot.executable_descriptor}",
            }
            plan_content, plan_counts = _canonical_replay_makefile_bytes(
                candidate, private_source, private_stage,
                canonical_transports, artifact_transports,
                execution_build_root=execution_stage_root)
            plan_path = temporary / _CANONICAL_REPLAY_PLAN_RELATIVE
            _write_private_executable(plan_path, plan_content)
            plan_snapshot = retained_tools.enter_context(
                _RetainedFileSnapshot(
                    plan_path, "canonical clean replay build plan"))
            plan_descriptor = plan_snapshot.executable_descriptor
            output_topology = _retain_canonical_replay_output_topology(
                candidate, private_source, temporary, retained_tools)
            recipe_transport_manifest = _canonical_replay_plan_manifest(
                candidate, plan_path, plan_content, plan_snapshot, plan_counts,
                artifact_transports, stage_descriptor)
            transport_descriptors.add(plan_descriptor)
            build_command = (
                str(cache["CMAKE_MAKE_PROGRAM"]),
                "-C", execution_stage_root,
                "-f", f"/proc/self/fd/{plan_descriptor}",
                "-j", str(jobs),
                f"SHELL=/proc/self/fd/"
                f"{bash_snapshot.executable_descriptor}",
                "--", f"CMakeFiles/{target}.dir/replay")
            stage_directories = sorted({
                parent.as_posix()
                for relative in artifact_sinks
                for parent in (Path(relative).parent,)
                if parent.as_posix() != "."
            })
            build_invocation = _replay_invocation_record(
                build_command, executable_label="make_program",
                executable=tool_snapshots["make_program"],
                inherited_descriptors=tuple(transport_descriptors),
                environment_overrides={"TMPDIR": str(private_tmp)},
                maximum_bytes=32 << 20, timeout=1800,
                write_sandbox_root=workspace,
                write_sandbox_descriptors=tuple(
                    sink.descriptor for sink in artifact_sinks.values()),
                private_tmpfs_root=private_stage,
                private_tmpfs_directories=stage_directories,
                private_tmpfs_descriptor=stage_descriptor)
            _run(
                build_command,
                "runner-owned candidate build",
                maximum_bytes=32 << 20, timeout=1800,
                inherited_descriptors=tuple(transport_descriptors),
                executable_descriptor=tool_snapshots[
                    "make_program"].executable_descriptor,
                environment_overrides={"TMPDIR": str(private_tmp)},
                write_sandbox_root=workspace,
                write_sandbox_descriptors=tuple(
                    sink.descriptor for sink in artifact_sinks.values()),
                private_tmpfs_root=private_stage,
                private_tmpfs_directories=stage_directories,
                private_tmpfs_descriptor=stage_descriptor)
            header_relative = (
                "generated/leopard2-benchmark-attestation/"
                "leopard2_benchmark_source_attestation.h")
            for relative, sink in artifact_sinks.items():
                sink.seal(
                    executable=(relative == target),
                    allow_empty=(relative == header_relative + ".lock"))
            plan_snapshot.verify()
            _verify_canonical_replay_output_topology(
                output_topology, require_outputs=True)
            final_sinks = {
                sink.logical_path: sink for sink in artifact_sinks.values()
            }
            for record in output_topology["outputs"]:
                output = record["path"]
                if output == plan_path:
                    continue
                sink = final_sinks.get(output)
                require(sink is not None,
                        f"canonical replay output has no sealed artifact: "
                        f"{output}")
                descriptor = record["file"].fileno()
                os.ftruncate(descriptor, 0)
                os.lseek(descriptor, 0, os.SEEK_SET)
                view = memoryview(sink.content)
                while view:
                    written = os.write(descriptor, view)
                    require(written > 0,
                            f"canonical replay output materialization "
                            f"stalled: {output}")
                    view = view[written:]
                os.fsync(descriptor)
            _verify_canonical_replay_output_topology(
                output_topology, require_outputs=True)
            rebuilt = _capture_replayed_candidate_provenance(
                temporary, private_source, target,
                replay_contract=replay_contract,
                inherited_descriptors=inherited_descriptors,
                tracked_source_manifest=source_tree.manifest,
                logical_source_root=source,
                sealed_artifacts=final_sinks)
            result = compare_reproducible_builds(candidate, rebuilt)
            _verify_canonical_replay_output_topology(
                output_topology, require_outputs=True)
            result["schema"] = REPRODUCIBLE_BUILD_PROOF_SCHEMA
            immutable_tools = []
            for label, snapshot in (
                    *tuple(tool_snapshots.items()),
                    ("cmake", cmake_snapshot),
                    ("command_shell", shell_snapshot),
                    ("wrapper_shell", bash_snapshot),
                    ("copy", copy_snapshot),
                    ("replay_git_identity", replay_git_snapshot),
                    ("empty_dependency_include",
                     empty_dependency_snapshot),
                    ("private_output_inventory", inventory_snapshot)):
                immutable_tools.append({
                    "label": label,
                    "source_path": snapshot.identity["path"],
                    "sealed_descriptor": snapshot.executable_descriptor,
                    **snapshot.executable_identity,
                })
            for record, snapshot in subtool_snapshots:
                immutable_tools.append({
                    "label": (
                        f"{record['language']}-compiler-"
                        f"{record['role']}"),
                    "source_path": snapshot.identity["path"],
                    "sealed_descriptor": snapshot.executable_descriptor,
                    **snapshot.executable_identity,
                })
            result["immutable_replay"] = {
                "schema": "leopard2-immutable-replay-tools/v2",
                "tools": sorted(
                    immutable_tools,
                    key=lambda item: (item["label"], item["source_path"])),
                "wrappers": sorted(
                    wrapper_records, key=lambda item: item["label"]),
                "compiler_prefixes": sorted(
                    prefix_records, key=lambda item: item["language"]),
                "recipe_transport": recipe_transport_manifest,
                "artifacts": [
                    {
                        "path": relative,
                        "sealed_descriptor": sink.descriptor,
                        "sha256": sink.identity["sha256"],
                        "size": sink.identity["size"],
                        "mode": sink.identity["mode"],
                        "seals": fcntl.fcntl(
                            sink.descriptor,
                            getattr(fcntl, "F_GET_SEALS", 1034)),
                    }
                    for relative, sink in sorted(artifact_sinks.items())
                ],
                "git_identity": {
                    **replay_git_specification,
                    "script_sha256":
                        hashlib.sha256(replay_git_bytes).hexdigest(),
                    "script_size": len(replay_git_bytes),
                    "sealed_descriptor":
                        replay_git_snapshot.executable_descriptor,
                    "sealed_sha256":
                        replay_git_snapshot.executable_identity["sha256"],
                    "sealed_size":
                        replay_git_snapshot.executable_identity["size"],
                    "sealed_seals":
                        replay_git_snapshot.executable_identity["seals"],
                },
                "invocations": {
                    "configure": configure_invocation,
                    "build": build_invocation,
                },
                "replacements": [
                    {
                        "label": label,
                        "logical_path": logical,
                        "transport_path": transport,
                    }
                    for label, (logical, transport) in sorted(
                        replacements.items())
                ],
            }
            source_tree.verify()
            validate_reproducible_build_proof(
                result, candidate, label="generated candidate")
            # ExitStack and source_tree now revalidate every held tool and
            # source pathname, including the inotify event history that makes
            # an A -> B -> A replacement observable.
            return result
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def validate_reproducible_build_proof(
    proof: Any, candidate: Mapping[str, Any], *,
    label: str = "candidate",
) -> None:
    """Strictly validate one exact current or historical replay contract."""

    def exact(value: Any, keys: set[str], what: str) -> Mapping[str, Any]:
        require(isinstance(value, Mapping) and set(value) == keys,
                f"{label} {what} is malformed")
        return value

    def digest(value: Any, what: str) -> str:
        require(isinstance(value, str) and
                re.fullmatch(r"[0-9a-f]{64}", value) is not None,
                f"{label} {what} digest is malformed")
        return value

    contract = _reproducible_replay_contract(candidate)
    expected_core = compare_reproducible_builds(candidate, candidate)
    require(expected_core.get("schema") == contract["proof_schema"],
            f"{label} clean reproducible-build core schema differs")
    proof_map = exact(
        proof, {*expected_core, "immutable_replay"},
        "clean reproducible-build proof")
    for key, expected in expected_core.items():
        require(proof_map[key] == expected,
                f"{label} clean reproducible-build proof field {key!r} "
                "is not bound to its closure")

    replay_keys = {
        "schema", "tools", "wrappers", "compiler_prefixes",
        "recipe_transport", "git_identity", "invocations",
        "replacements",
    }
    if contract["closure_schema"] == PRODUCTION_BUILD_CLOSURE_SCHEMA:
        replay_keys.add("artifacts")
    replay = exact(
        proof_map["immutable_replay"],
        replay_keys,
        "immutable replay proof")
    require(replay["schema"] == "leopard2-immutable-replay-tools/v2",
            f"{label} immutable replay schema is unsupported")

    tool_keys = {
        "label", "source_path", "sealed_descriptor", "sha256", "size",
        "seals", "source_sha256",
    }
    tools = replay["tools"]
    require(isinstance(tools, list) and tools,
            f"{label} immutable replay tool list is empty")
    tool_by_label: dict[str, Mapping[str, Any]] = {}
    required_seals = (
        getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
        getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
        getattr(fcntl, "F_SEAL_GROW", 0x0004) |
        getattr(fcntl, "F_SEAL_WRITE", 0x0008))
    for index, raw in enumerate(tools):
        tool = exact(raw, tool_keys, f"immutable tool {index}")
        tool_label = tool["label"]
        require(isinstance(tool_label, str) and tool_label and
                tool_label not in tool_by_label and
                isinstance(tool["source_path"], str) and
                Path(tool["source_path"]).is_absolute() and
                type(tool["sealed_descriptor"]) is int and
                tool["sealed_descriptor"] >= 0 and
                type(tool["size"]) is int and tool["size"] > 0 and
                type(tool["seals"]) is int and
                tool["seals"] & required_seals == required_seals,
                f"{label} immutable tool {index} identity is malformed")
        digest(tool["sha256"], f"immutable tool {tool_label} sealed")
        digest(
            tool["source_sha256"],
            f"immutable tool {tool_label} source")
        require(tool["sha256"] == tool["source_sha256"],
                f"{label} immutable tool {tool_label} differs from source")
        tool_by_label[tool_label] = tool
    require(tools == sorted(
                tools, key=lambda item: (item["label"], item["source_path"])),
            f"{label} immutable replay tools are not canonical")

    base_tool_labels = {
        "c_compiler": "c_compiler",
        "compiler": "compiler",
        "archiver": "archiver",
        "ranlib": "ranlib",
        "make_program": "make_program",
        "linker": "linker",
        "benchmark_git": "benchmark_git",
    }
    expected_tool_labels = set(base_tool_labels)
    for tool_label, closure_key in base_tool_labels.items():
        identity = candidate.get(closure_key)
        require(isinstance(identity, Mapping),
                f"{label} closure lacks {closure_key}")
        tool = tool_by_label.get(tool_label)
        require(tool is not None and
                tool["source_path"] == identity.get("path") and
                tool["source_sha256"] == identity.get("sha256") and
                tool["size"] == identity.get("size"),
                f"{label} immutable {tool_label} is not bound to closure")
    subtool_roles: dict[str, dict[str, str]] = {"c": {}, "c++": {}}
    subtools = candidate.get("compiler_subtools")
    require(isinstance(subtools, list) and subtools,
            f"{label} closure lacks compiler subtools")
    for record in subtools:
        require(isinstance(record, Mapping) and
                record.get("language") in subtool_roles and
                isinstance(record.get("role"), str) and
                isinstance(record.get("identity"), Mapping),
                f"{label} compiler subtool record is malformed")
        language = str(record["language"])
        role = str(record["role"])
        tool_label = f"{language}-compiler-{role}"
        require(role not in subtool_roles[language],
                f"{label} compiler subtool role is duplicated")
        subtool_roles[language][role] = tool_label
        expected_tool_labels.add(tool_label)
        identity = record["identity"]
        tool = tool_by_label.get(tool_label)
        require(tool is not None and
                tool["source_path"] == identity.get("path") and
                tool["source_sha256"] == identity.get("sha256") and
                tool["size"] == identity.get("size"),
                f"{label} immutable {tool_label} is not bound to closure")
    expected_tool_labels.update({
        "cmake", "command_shell", "wrapper_shell", "replay_git_identity",
        "empty_dependency_include",
    })
    if contract["closure_schema"] == PRODUCTION_BUILD_CLOSURE_SCHEMA:
        expected_tool_labels.update(("copy", "private_output_inventory"))
    require(set(tool_by_label) == expected_tool_labels,
            f"{label} immutable replay tool roles differ")
    empty_dependency_tool = tool_by_label["empty_dependency_include"]
    require(
        empty_dependency_tool["sha256"] ==
            hashlib.sha256(REPLAY_EMPTY_DEPENDENCY_INCLUDE).hexdigest() and
        empty_dependency_tool["source_sha256"] ==
            hashlib.sha256(REPLAY_EMPTY_DEPENDENCY_INCLUDE).hexdigest() and
        empty_dependency_tool["size"] ==
            len(REPLAY_EMPTY_DEPENDENCY_INCLUDE),
        f"{label} immutable dependency include is not inert")
    prefix_keys = {
        "language", "path", "descriptor", "mode",
        "mappings",
    }
    mapping_keys = {
        "role", "target_descriptor", "source_path", "source_sha256",
        "sealed_sha256", "sealed_size", "sealed_seals",
    }
    prefixes = replay["compiler_prefixes"]
    require(isinstance(prefixes, list) and len(prefixes) == 2,
            f"{label} compiler-prefix proof is malformed")
    prefix_by_language: dict[str, Mapping[str, Any]] = {}
    for raw in prefixes:
        prefix = exact(raw, prefix_keys, "compiler-prefix record")
        language = prefix["language"]
        require(isinstance(language, str) and
                language in ("c", "c++") and
                language not in prefix_by_language and
                isinstance(prefix["path"], str) and
                Path(prefix["path"]).is_absolute() and
                type(prefix["descriptor"]) is int and
                prefix["descriptor"] >= 0 and
                type(prefix["mode"]) is int and prefix["mode"] >= 0 and
                stat.S_ISDIR(prefix["mode"]),
                f"{label} compiler-prefix identity is malformed")
        role_tools = dict(subtool_roles[language])
        role_tools.setdefault("ld", "linker")
        mappings = prefix["mappings"]
        require(isinstance(mappings, list) and
                len(mappings) == len(role_tools),
                f"{label} {language} compiler-prefix mappings differ")
        observed_roles = set()
        for raw_mapping in mappings:
            mapping = exact(
                raw_mapping, mapping_keys,
                f"{language} compiler-prefix mapping")
            role = mapping["role"]
            require(isinstance(role, str),
                    f"{label} {language} compiler-prefix role differs")
            tool_label = role_tools.get(role)
            require(tool_label is not None and role not in observed_roles,
                    f"{label} {language} compiler-prefix role differs")
            observed_roles.add(role)
            tool = tool_by_label[tool_label]
            require(
                mapping["target_descriptor"] ==
                    tool["sealed_descriptor"] and
                mapping["source_path"] == tool["source_path"] and
                mapping["source_sha256"] == tool["source_sha256"] and
                mapping["sealed_sha256"] == tool["sha256"] and
                mapping["sealed_size"] == tool["size"] and
                mapping["sealed_seals"] == tool["seals"],
                f"{label} {language} compiler-prefix role {role} is "
                "not bound to its immutable tool")
        require(observed_roles == set(role_tools),
                f"{label} {language} compiler-prefix roles differ")
        require(mappings == sorted(
                    mappings, key=lambda item: item["role"]),
                f"{label} {language} compiler-prefix mappings are not "
                "canonical")
        prefix_by_language[language] = prefix
    require(prefixes == sorted(prefixes, key=lambda item: item["language"]),
            f"{label} compiler-prefix records are not canonical")

    wrapper_keys = {
        "label", "underlying_label", "logical_argv0",
        "injected_arguments", "interpreter_descriptor",
        "interpreter_sha256", "underlying_descriptor",
        "underlying_source_path", "underlying_source_sha256",
        "underlying_sealed_sha256", "script_sha256", "script_size",
        "sealed_descriptor", "sealed_sha256", "sealed_size",
        "sealed_seals",
    }
    wrapper_roles = {
        "cc": "c_compiler",
        "cxx": "compiler",
        "ar": "archiver",
        "ranlib": "ranlib",
        "ld": "linker",
        "cmake": "cmake",
        "make": "make_program",
    }
    if contract["closure_schema"] == PRODUCTION_BUILD_CLOSURE_SCHEMA:
        wrapper_roles["copy"] = "copy"
    wrappers = replay["wrappers"]
    require(isinstance(wrappers, list) and len(wrappers) == len(wrapper_roles),
            f"{label} immutable replay wrappers are malformed")
    wrapper_by_label: dict[str, Mapping[str, Any]] = {}
    cache = candidate.get("validated_cache")
    require(isinstance(cache, Mapping),
            f"{label} closure lacks retained CMake cache semantics")
    logical_wrapper_paths = {
        "cc": cache.get("CMAKE_C_COMPILER"),
        "cxx": cache.get("CMAKE_CXX_COMPILER"),
        "ar": cache.get("CMAKE_AR"),
        "ranlib": cache.get("CMAKE_RANLIB"),
        "ld": cache.get("CMAKE_LINKER"),
        "cmake": tool_by_label["cmake"]["source_path"],
        "make": cache.get("CMAKE_MAKE_PROGRAM"),
    }
    if "copy" in wrapper_roles:
        logical_wrapper_paths["copy"] = tool_by_label["copy"]["source_path"]
    for raw in wrappers:
        wrapper = exact(raw, wrapper_keys, "immutable replay wrapper")
        wrapper_label = wrapper["label"]
        require(isinstance(wrapper_label, str),
                f"{label} immutable wrapper role is malformed")
        underlying_label = wrapper_roles.get(wrapper_label)
        require(underlying_label is not None and
                wrapper_label not in wrapper_by_label and
                wrapper["underlying_label"] == underlying_label and
                wrapper["logical_argv0"] ==
                    logical_wrapper_paths[wrapper_label] and
                isinstance(wrapper["injected_arguments"], list) and
                all(isinstance(argument, str) and argument
                    for argument in wrapper["injected_arguments"]),
                f"{label} immutable wrapper role is malformed")
        expected_arguments: list[str] = []
        if wrapper_label == "cc":
            expected_arguments = [
                f"-B/proc/self/fd/"
                f"{prefix_by_language['c']['descriptor']}/"]
        elif wrapper_label == "cxx":
            expected_arguments = [
                f"-B/proc/self/fd/"
                f"{prefix_by_language['c++']['descriptor']}/"]
        require(wrapper["injected_arguments"] == expected_arguments,
                f"{label} immutable {wrapper_label} wrapper arguments differ")
        interpreter = tool_by_label["wrapper_shell"]
        underlying = tool_by_label[underlying_label]
        require(
            wrapper["interpreter_descriptor"] ==
                interpreter["sealed_descriptor"] and
            wrapper["interpreter_sha256"] == interpreter["sha256"] and
            wrapper["underlying_descriptor"] ==
                underlying["sealed_descriptor"] and
            wrapper["underlying_source_path"] ==
                underlying["source_path"] and
            wrapper["underlying_source_sha256"] ==
                underlying["source_sha256"] and
            wrapper["underlying_sealed_sha256"] ==
                underlying["sha256"],
            f"{label} immutable {wrapper_label} wrapper binding differs")
        script = _replay_exec_wrapper_bytes(
            wrapper["interpreter_descriptor"],
            wrapper["logical_argv0"],
            wrapper["underlying_descriptor"],
            wrapper["injected_arguments"])
        script_sha256 = hashlib.sha256(script).hexdigest()
        require(
            wrapper["script_sha256"] == script_sha256 and
            wrapper["script_size"] == len(script) and
            wrapper["sealed_sha256"] == script_sha256 and
            wrapper["sealed_size"] == len(script) and
            type(wrapper["sealed_descriptor"]) is int and
            wrapper["sealed_descriptor"] >= 0 and
            type(wrapper["sealed_seals"]) is int and
            wrapper["sealed_seals"] & required_seals == required_seals,
            f"{label} immutable {wrapper_label} wrapper bytes differ")
        wrapper_by_label[wrapper_label] = wrapper
    require(set(wrapper_by_label) == set(wrapper_roles) and
            wrappers == sorted(wrappers, key=lambda item: item["label"]),
            f"{label} immutable replay wrapper roles differ")

    replacement_keys = {"label", "logical_path", "transport_path"}
    replacements = replay["replacements"]
    replacement_wrappers = {
        "c-compiler": "cc", "cxx-compiler": "cxx", "archiver": "ar",
        "ranlib": "ranlib", "linker": "ld", "cmake": "cmake",
        "make": "make",
    }
    logical_replacements = {
        "c-compiler": cache.get("CMAKE_C_COMPILER"),
        "cxx-compiler": cache.get("CMAKE_CXX_COMPILER"),
        "archiver": cache.get("CMAKE_AR"),
        "ranlib": cache.get("CMAKE_RANLIB"),
        "linker": cache.get("CMAKE_LINKER"),
        "cmake": tool_by_label["cmake"]["source_path"],
        "make": cache.get("CMAKE_MAKE_PROGRAM"),
        "git": cache.get("LEO2_BENCHMARK_GIT_EXECUTABLE"),
        "dependency-includes":
            "CMakeFiles/**/{compiler_,}depend.make",
    }
    replacement_by_label: dict[str, Mapping[str, Any]] = {}
    require(isinstance(replacements, list),
            f"{label} immutable replacement list is malformed")
    for raw in replacements:
        replacement = exact(
            raw, replacement_keys, "immutable replay replacement")
        replacement_label = replacement["label"]
        require(isinstance(replacement_label, str) and
                replacement_label in logical_replacements and
                replacement_label not in replacement_by_label and
                replacement["logical_path"] ==
                    logical_replacements[replacement_label],
                f"{label} immutable replay replacement identity differs")
        if replacement_label == "git":
            descriptor = tool_by_label[
                "replay_git_identity"]["sealed_descriptor"]
        elif replacement_label == "dependency-includes":
            descriptor = tool_by_label[
                "empty_dependency_include"]["sealed_descriptor"]
        else:
            descriptor = wrapper_by_label[
                replacement_wrappers[replacement_label]][
                    "sealed_descriptor"]
        require(replacement["transport_path"] ==
                f"/proc/self/fd/{descriptor}",
                f"{label} immutable {replacement_label} transport differs")
        replacement_by_label[replacement_label] = replacement
    require(set(replacement_by_label) == set(logical_replacements) and
            replacements == sorted(
                replacements, key=lambda item: item["label"]),
            f"{label} immutable replay replacement roles differ")

    git_keys = {
        "source_root", "commit", "tree", "dirty",
        "interpreter_descriptor", "script_sha256", "script_size",
        "sealed_descriptor", "sealed_sha256", "sealed_size",
        "sealed_seals",
    }
    git_identity = exact(
        replay["git_identity"], git_keys, "replay Git identity")
    manifest = candidate.get("tracked_source_manifest")
    require(isinstance(manifest, Mapping) and
            isinstance(manifest.get("git"), Mapping) and
            isinstance(git_identity["source_root"], str) and
            Path(git_identity["source_root"]).is_absolute() and
            git_identity["commit"] == manifest["git"].get("commit") and
            git_identity["tree"] == manifest["git"].get("tree") and
            git_identity["dirty"] == manifest["git"].get("dirty") and
            git_identity["interpreter_descriptor"] ==
                tool_by_label["wrapper_shell"]["sealed_descriptor"],
            f"{label} replay Git semantic identity differs")
    git_script, git_specification = _replay_git_identity_bytes(
        Path(git_identity["source_root"]), manifest,
        interpreter_descriptor=git_identity["interpreter_descriptor"])
    require(all(git_identity[key] == value
                for key, value in git_specification.items()),
            f"{label} replay Git specification differs")
    git_digest = hashlib.sha256(git_script).hexdigest()
    replay_git_tool = tool_by_label["replay_git_identity"]
    require(
        git_identity["script_sha256"] == git_digest and
        git_identity["script_size"] == len(git_script) and
        git_identity["sealed_descriptor"] ==
            replay_git_tool["sealed_descriptor"] and
        git_identity["sealed_sha256"] == git_digest and
        git_identity["sealed_size"] == len(git_script) and
        git_identity["sealed_seals"] == replay_git_tool["seals"] and
        replay_git_tool["sha256"] == git_digest,
        f"{label} replay Git executable bytes differ")

    recipe_transport = replay["recipe_transport"]
    require(isinstance(recipe_transport, Mapping) and
            isinstance(recipe_transport.get("schema"), str),
            f"{label} replay recipe transport is malformed")
    _require_reproducible_replay_artifact_contract(
        candidate, proof_map["schema"], recipe_transport["schema"])
    canonical_recipe = (
        recipe_transport
        if recipe_transport["schema"] == CANONICAL_REPLAY_RECIPE_SCHEMA
        else None)
    if contract["recipe_schema"] == LEGACY_REPLAY_RECIPE_SCHEMA:
        _validate_legacy_replay_recipe_transport(
            recipe_transport, candidate, replacement_by_label, label=label)
    else:
        require(canonical_recipe is not None,
                f"{label} current replay requires a canonical build plan")

    artifact_by_path: dict[str, Mapping[str, Any]] = {}
    if contract["closure_schema"] == PRODUCTION_BUILD_CLOSURE_SCHEMA:
        artifact_keys = {
            "path", "sealed_descriptor", "sha256", "size", "mode", "seals",
        }
        artifacts = replay["artifacts"]
        require(isinstance(artifacts, list) and artifacts,
                f"{label} sealed replay artifact list is empty")
        for raw_artifact in artifacts:
            artifact = exact(
                raw_artifact, artifact_keys, "sealed replay artifact")
            path = artifact["path"]
            require(
                isinstance(path, str) and path and
                not Path(path).is_absolute() and
                all(component not in ("", ".", "..")
                    for component in Path(path).parts) and
                path not in artifact_by_path and
                type(artifact["sealed_descriptor"]) is int and
                artifact["sealed_descriptor"] >= 0 and
                type(artifact["size"]) is int and artifact["size"] >= 0 and
                type(artifact["mode"]) is int and
                stat.S_ISREG(artifact["mode"]) and
                type(artifact["seals"]) is int and
                artifact["seals"] & required_seals == required_seals,
                f"{label} sealed replay artifact identity is malformed")
            digest(artifact["sha256"], f"sealed replay artifact {path}")
            artifact_by_path[path] = artifact
        require(
            artifacts == sorted(artifacts, key=lambda item: item["path"]) and
            canonical_recipe is not None and
            canonical_recipe["output_sinks"] == [
                {
                    "path": path,
                    "descriptor":
                        artifact_by_path[path]["sealed_descriptor"],
                }
                for path in sorted(artifact_by_path)
            ] and
            canonical_recipe["counts"]["freeze_output_count"] ==
                len(artifact_by_path),
            f"{label} sealed artifacts are not bound to canonical sinks")
        inventory_script = _private_replay_inventory_bytes(
            tool_by_label["wrapper_shell"]["sealed_descriptor"],
            sorted(artifact_by_path))
        inventory_tool = tool_by_label["private_output_inventory"]
        inventory_digest = hashlib.sha256(inventory_script).hexdigest()
        require(
            inventory_tool["sha256"] == inventory_digest and
            inventory_tool["source_sha256"] == inventory_digest and
            inventory_tool["size"] == len(inventory_script),
            f"{label} private replay inventory verifier bytes differ")

    invocation_keys = {
        "schema", "argv", "environment", "pass_fds",
        "maximum_output_bytes", "timeout_seconds", "executable",
    }
    if contract["invocation_schema"] == REPLAY_INVOCATION_SCHEMA:
        invocation_keys.update({
            "write_sandbox_root", "write_sandbox_descriptors",
            "private_tmpfs",
        })
    executable_keys = {
        "label", "logical_argv0", "sealed_descriptor", "sealed_sha256",
        "sealed_size", "sealed_seals",
    }

    def invocation(
        raw: Any, name: str, tool_label: str, environment: Mapping[str, str],
        maximum_bytes: int, timeout: int,
        *, write_sandbox_root: Path | None = None,
        write_sandbox_descriptors: Sequence[int] = (),
        private_tmpfs: Mapping[str, Any] | None = None,
    ) -> Mapping[str, Any]:
        value = exact(raw, invocation_keys, f"{name} invocation")
        executable = exact(
            value["executable"], executable_keys,
            f"{name} invocation executable")
        tool = tool_by_label[tool_label]
        pass_fds = value["pass_fds"]
        require(
            value["schema"] == contract["invocation_schema"] and
            isinstance(value["argv"], list) and value["argv"] and
            all(isinstance(argument, str) and "\0" not in argument
                for argument in value["argv"]) and
            value["environment"] == {
                key: environment[key] for key in sorted(environment)
            } and
            isinstance(pass_fds, list) and
            all(type(descriptor) is int and descriptor >= 0
                for descriptor in pass_fds) and
            pass_fds == sorted(set(pass_fds)) and
            value["maximum_output_bytes"] == maximum_bytes and
            value["timeout_seconds"] == timeout and
            executable["label"] == tool_label and
            executable["logical_argv0"] == value["argv"][0] and
            executable["sealed_descriptor"] ==
                tool["sealed_descriptor"] and
            executable["sealed_sha256"] == tool["sha256"] and
            executable["sealed_size"] == tool["size"] and
            executable["sealed_seals"] == tool["seals"] and
            executable["sealed_descriptor"] in pass_fds,
            f"{label} {name} replay invocation differs")
        if contract["invocation_schema"] == REPLAY_INVOCATION_SCHEMA:
            require(
                write_sandbox_root is not None and
                value["write_sandbox_root"] ==
                    str(write_sandbox_root) and
                value["write_sandbox_descriptors"] ==
                    sorted(set(write_sandbox_descriptors)) and
                value["private_tmpfs"] == private_tmpfs,
                f"{label} {name} replay isolation differs")
        return value

    invocations = exact(
        replay["invocations"], {"configure", "build"},
        "replay invocation set")
    replay_source_root = Path(str(git_identity["source_root"]))
    replay_workspace = replay_source_root.parent
    configure_environment = dict(GIT_ENVIRONMENT)
    configure_environment.update({
        "CC": str(cache.get("CMAKE_C_COMPILER")),
        "CXX": str(cache.get("CMAKE_CXX_COMPILER")),
    })
    if contract["invocation_schema"] == REPLAY_INVOCATION_SCHEMA:
        configure_environment["TMPDIR"] = str(replay_workspace / "tmp")
    configure_record = invocation(
        invocations["configure"], "configure", "cmake",
        configure_environment, 16 << 20, 300,
        write_sandbox_root=(
            replay_workspace
            if contract["invocation_schema"] == REPLAY_INVOCATION_SCHEMA
            else None))
    configure_argv = configure_record["argv"]
    require(len(configure_argv) >= 7 and
            configure_argv[1] == "-S" and
            configure_argv[2] == git_identity["source_root"] and
            configure_argv[3] == "-B",
            f"{label} replay configure roots are malformed")
    replay_build_root = Path(configure_argv[4])
    require(replay_source_root.is_absolute() and
            replay_build_root.is_absolute() and
            replay_source_root.name == "source" and
            replay_build_root.name == "build" and
            replay_source_root.parent == replay_build_root.parent,
            f"{label} replay workspace roots are malformed")
    replay_tools_root = replay_workspace / "tools"
    require(
        Path(prefix_by_language["c"]["path"]) ==
            replay_tools_root / "c-prefix" and
        Path(prefix_by_language["c++"]["path"]) ==
            replay_tools_root / "cxx-prefix" and
        stat.S_IMODE(prefix_by_language["c"]["mode"]) == 0o700 and
        stat.S_IMODE(prefix_by_language["c++"]["mode"]) == 0o700 and
        Path(tool_by_label["replay_git_identity"]["source_path"]) ==
            replay_tools_root / "git" and
        Path(tool_by_label["empty_dependency_include"]["source_path"]) ==
            replay_tools_root / "empty-dependency-include" and
        (
            contract["closure_schema"] != PRODUCTION_BUILD_CLOSURE_SCHEMA or
            (
                Path(tool_by_label[
                    "private_output_inventory"]["source_path"]) ==
                    replay_tools_root / "verify-output-inventory" and
                Path(tool_by_label["copy"]["source_path"]).name == "cat"
            )
        ),
        f"{label} replay private-tool topology differs")
    expected_configure = _reproducible_configure_argv(
        Path(git_identity["source_root"]), replay_build_root, cache,
        cmake_path=tool_by_label["cmake"]["source_path"])
    require(configure_argv == expected_configure,
            f"{label} replay configure argv differs")

    canonical_plan_descriptor: int | None = None
    replay_stage_root = replay_workspace / "stage"
    if canonical_recipe is not None:
        canonical_transports = {
            "cxx":
                f"/proc/self/fd/"
                f"{wrapper_by_label['cxx']['sealed_descriptor']}",
            "archiver":
                f"/proc/self/fd/"
                f"{wrapper_by_label['ar']['sealed_descriptor']}",
            "ranlib":
                f"/proc/self/fd/"
                f"{wrapper_by_label['ranlib']['sealed_descriptor']}",
            "cmake":
                f"/proc/self/fd/"
                f"{wrapper_by_label['cmake']['sealed_descriptor']}",
            "git":
                f"/proc/self/fd/"
                f"{tool_by_label['replay_git_identity']['sealed_descriptor']}",
            "copy":
                f"/proc/self/fd/"
                f"{wrapper_by_label['copy']['sealed_descriptor']}",
            "inventory":
                f"/proc/self/fd/"
                f"{tool_by_label[
                    'private_output_inventory']['sealed_descriptor']}",
        }
        canonical_plan_descriptor = \
            _validate_canonical_replay_plan_manifest(
                canonical_recipe, candidate, replay_source_root,
                replay_stage_root, canonical_transports, label=label)
        _validate_sealed_replay_artifact_bindings(
            candidate, artifact_by_path, replay_source_root,
            replay_stage_root, label=label)

    build_environment = dict(GIT_ENVIRONMENT)
    artifact_descriptors = [
        artifact["sealed_descriptor"]
        for artifact in artifact_by_path.values()
    ]
    stage_directories = sorted({
        Path(path).parent.as_posix()
        for path in artifact_by_path
        if Path(path).parent.as_posix() != "."
    })
    expected_private_tmpfs = None
    if contract["invocation_schema"] == REPLAY_INVOCATION_SCHEMA:
        build_environment["TMPDIR"] = str(replay_workspace / "tmp")
        expected_private_tmpfs = {
            "root": str(replay_stage_root),
            "descriptor":
                canonical_recipe["execution_root_descriptor"]
                if canonical_recipe is not None else None,
            "bytes": PRIVATE_REPLAY_TMPFS_BYTES,
            "directories": stage_directories,
        }
    build_record = invocation(
        invocations["build"], "build", "make_program",
        build_environment, 32 << 20, 1800,
        write_sandbox_root=(
            replay_workspace
            if contract["invocation_schema"] == REPLAY_INVOCATION_SCHEMA
            else None),
        write_sandbox_descriptors=artifact_descriptors,
        private_tmpfs=expected_private_tmpfs)
    build_argv = build_record["argv"]
    if canonical_plan_descriptor is None:
        require(len(build_argv) == 11 and
                build_argv[0] == cache.get("CMAKE_MAKE_PROGRAM") and
                build_argv[1:6] == [
                    "-C", str(replay_build_root),
                    "-f", "CMakeFiles/Makefile2", "-j"] and
                re.fullmatch(
                    r"[1-9][0-9]{0,2}", build_argv[6]) is not None and
                1 <= int(build_argv[6]) <= 128 and
                build_argv[7] ==
                    f"SHELL=/proc/self/fd/"
                    f"{tool_by_label['command_shell']['sealed_descriptor']}"
                and
                build_argv[8] ==
                    f"MAKE=/proc/self/fd/"
                    f"{wrapper_by_label['make']['sealed_descriptor']}" and
                build_argv[9:] == [
                    "--",
                    f"CMakeFiles/"
                    f"{candidate.get('executable_target')}.dir/all"],
                f"{label} replay build argv differs")
    else:
        require(len(build_argv) == 10 and
                build_argv[0] == cache.get("CMAKE_MAKE_PROGRAM") and
                build_argv[1:6] == [
                    "-C",
                    f"/proc/self/fd/"
                    f"{canonical_recipe['execution_root_descriptor']}",
                    "-f", f"/proc/self/fd/{canonical_plan_descriptor}", "-j"]
                and
                re.fullmatch(
                    r"[1-9][0-9]{0,2}", build_argv[6]) is not None and
                1 <= int(build_argv[6]) <= 128 and
                build_argv[7] ==
                    f"SHELL=/proc/self/fd/"
                    f"{tool_by_label['wrapper_shell']['sealed_descriptor']}"
                and
                build_argv[8:] == [
                    "--",
                    f"CMakeFiles/"
                    f"{candidate.get('executable_target')}.dir/replay"],
                f"{label} canonical replay build argv differs")
    required_build_fds = {
        *(tool["sealed_descriptor"] for tool in tool_by_label.values()),
        *(wrapper["sealed_descriptor"]
          for wrapper in wrapper_by_label.values()),
        *(prefix["descriptor"] for prefix in prefix_by_language.values()),
    }
    if canonical_plan_descriptor is not None:
        required_build_fds.add(canonical_plan_descriptor)
        required_build_fds.add(
            canonical_recipe["execution_root_descriptor"])
        required_build_fds.update(artifact_descriptors)
    require(required_build_fds.issubset(set(build_record["pass_fds"])),
            f"{label} replay build descriptor transport is incomplete")

    all_internal_descriptors = [
        *(tool["sealed_descriptor"] for tool in tool_by_label.values()),
        *(wrapper["sealed_descriptor"]
          for wrapper in wrapper_by_label.values()),
        *(prefix["descriptor"] for prefix in prefix_by_language.values()),
    ]
    if canonical_plan_descriptor is not None:
        all_internal_descriptors.append(canonical_plan_descriptor)
        all_internal_descriptors.append(
            canonical_recipe["execution_root_descriptor"])
        all_internal_descriptors.extend(artifact_descriptors)
    require(len(all_internal_descriptors) ==
            len(set(all_internal_descriptors)),
            f"{label} immutable replay descriptor identities collide")
